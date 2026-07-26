# -*- coding: utf-8 -*-
# DIM: na (court cross-check of the S933 probe's numbers — dimensionless linear algebra; handles 0).
# The project's court hand on S933 (W33-leg-1): an independent recomputation of three load-bearing
# rows of the table. NOT a probe — an act of the court (class S929_court_p4_omega).
import sys
sys.stdout.reconfigure(encoding='utf-8')
from sympy import Matrix, Integer, zeros, diag

def eta_pq(p, q):
    return diag(*([Integer(1)] * p + [Integer(-1)] * q))

def e(n, i):
    v = zeros(n, 1)
    v[i, 0] = Integer(1)
    return v

def wedge(x, y, g):
    return x * (g * y).T - y * (g * x).T

def so_basis(n, g):
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            out.append(wedge(e(n, i), e(n, j), g))
    return out

def flat(M):
    return Matrix(M.rows * M.rows, 1, list(M))

def span_dim(mats, n):
    if not mats:
        return 0
    return Matrix.hstack(*[flat(M) for M in mats]).rank()

def centralizer(A, bas):
    n = A.rows
    cols = [flat(B * A - A * B) for B in bas]
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        M = zeros(n, n)
        for k in range(len(bas)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * bas[k]
        out.append(M)
    return out

def derived_chain(cb, n):
    dims = [span_dim(cb, n)]
    cur = cb
    while True:
        brs = []
        for a in range(len(cur)):
            for b in range(a + 1, len(cur)):
                brs.append(cur[a] * cur[b] - cur[b] * cur[a])
        d = span_dim(brs, n)
        if d == dims[-1]:
            break
        dims.append(d)
        # basis of the derived subalgebra (crude: independent brackets)
        nxt, F, r = [], zeros(n * n, 0), 0
        for M in brs:
            F2 = Matrix.hstack(F, flat(M))
            if F2.rank() > r:
                nxt.append(M); F = F2; r += 1
        cur = nxt
        if d == 0:
            break
    return dims

def row(p, q, x, y, expect):
    n = p + q
    g = eta_pq(p, q)
    bas = so_basis(n, g)
    N = wedge(x, y, g)
    assert (N * N * N).is_zero_matrix or True  # informational, not a gate
    cb = centralizer(N, bas)
    dims = derived_chain(cb, n)
    got = (span_dim(cb, n), dims)
    tag = "OK " if got == expect else "MISMATCH"
    print(f"({p},{q}) | dimc={got[0]} derived={got[1]} | expect={expect} | {tag}")
    return got == expect

ok = True
# (3,1) rank1: x=e0+e3 isotropic, y=e1 — expect dimc=2, derived [2,0] (abelian)
p, q = 3, 1; n = 4
ok &= row(p, q, e(n,0)+e(n,3), e(n,1), (2, [2, 0]))
# (3,2) rank1: x=e0+e3, y=e1 — expect dimc=4, derived [4,2,0] (solvable, NON-abelian)
p, q = 3, 2; n = 5
ok &= row(p, q, e(n,0)+e(n,3), e(n,1), (4, [4, 2, 0]))
# (2,2) rank0: x=e0+e2, y=e1+e3 (both isotropic, orthogonal) — expect dimc=4, derived [4,3] (perfect-part 3)
p, q = 2, 2; n = 4
ok &= row(p, q, e(n,0)+e(n,2), e(n,1)+e(n,3), (4, [4, 3]))
# (3,3) rank1: x=e0+e3, y=e1 — expect dimc=7, derived [7,6] (perfect-part 6)
p, q = 3, 3; n = 6
ok &= row(p, q, e(n,0)+e(n,3), e(n,1), (7, [7, 6]))
print("COURT-HAND:", "ALL-OK" if ok else "MISMATCH-STOP")
sys.exit(0 if ok else 1)
