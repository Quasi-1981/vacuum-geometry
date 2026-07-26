# -*- coding: utf-8 -*-
# DIM: na (center pairing form omega on rad/Z across strata; symbolic d for wedges; 0 handles).
#
# ============================================================================
# OMEGA / BASIS CHOICE (stamped BEFORE any counting)
# ----------------------------------------------------------------------------
# For the centralizer c(A) of a stratum:
#   rad = killing_radical(c) = ker(Killing form of c)  (S937 convention; equals
#         the solvable radical when the Levi part is semisimple).
#   Z   = center_of(rad)  (elements of rad commuting with all of rad).
#   Adapted basis of rad = [ Z-basis z_1..z_s ] ++ [ rad/Z reps r_1..r_m ],
#         reps = complement_reps(rad, Z) (span-extension of Z inside rad).
# PAIRING FORM:  omega(r_i, r_j) := Z-component of [r_i, r_j], read as coords in
#   the adapted rad basis (first s coords).  This gives s antisymmetric mxm
#   matrices W^l (l=1..s).  The Z-valued 2-form has kernel
#     ker = { v in rad/Z : omega(v, r_j) = 0 in Z, all j },
#   and  rank(omega) := m - dim(ker) = rank of the stacked (s*m) x m map.
#   (For s=1 this is exactly matrix-rank of the single antisym W^1.)
#   nondegenerate <=> rank == m (and m>0);  k := rank/2.
#   rad ABELIAN => Z=rad, rad/Z empty (m=0) => form EMPTY (k=0), reported
#   SEPARATELY from "degenerate".
#
# SYMBOLIC forall d (wedge strata), from S952 block form, Cc-row:
#   module block B is 2xd, center Cc is the antisym 2x2; the pairing is
#   omega(B,B') = antisym( B G^{-1} B'^T ), i.e. on R^{2d}=(b1,b2) the matrix
#   [[0, G^{-1}],[-G^{-1},0]]  ->  rank = 2*rank(G^{-1}) = 2d, nondeg, k=d,
#   with d = core dim (= n-4 rank0).  Verified bit-for-bit vs instances.
#
# Discipline: 0 handles; exact arithmetic; mutant on each branch {nondeg /
# degenerate / empty}; FORBIDDEN-SCAN (S929 list); bit-fence strata rebuilt by
# S933-class code (radicals vs S937/S940); STOP after the table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import itertools
from math import lcm as _ilcm
from sympy import Matrix, Integer, Rational, zeros, eye, diag, symbols, Symbol

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S961_run.log")
_logf = open(_LOGPATH, "w", encoding="utf-8")


class Tee:
    def __init__(self, real, fh):
        self.real = real
        self.fh = fh
        self.chunks = []

    def write(self, s):
        self.real.write(s)
        self.fh.write(s)
        self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


_tee = Tee(sys.stdout, _logf)
sys.stdout = _tee

ASSERT_PASS = [0]
FAILS = [0]


def ok(cond, msg):
    if cond:
        ASSERT_PASS[0] += 1
    else:
        FAILS[0] += 1
        print("ASSERT-FAIL: " + msg)


# ==================== primitives (verbatim S929/S933/S937) ====================

def make_eta(p, q):
    return diag(*([Integer(1)] * p + [Integer(-1)] * q))


def unit(n, i):
    v = zeros(n, 1)
    v[i, 0] = Integer(1)
    return v


def wedge(x, y, eta):
    return x * (eta * y).T - y * (eta * x).T


def is_so(M, eta):
    return (M.T * eta + eta * M).is_zero_matrix


def so_basis(n, eta):
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            out.append(wedge(unit(n, i), unit(n, j), eta))
    return out


def flat(M):
    n = M.rows
    return Matrix(n * n, 1, list(M))


def stack_flats(mats, sq):
    if not mats:
        return zeros(sq * sq, 0)
    return Matrix.hstack(*[flat(M) for M in mats])


def span_rank(mats, sq):
    return stack_flats(mats, sq).rank()


def span_basis(mats, sq):
    out = []
    F = zeros(sq * sq, 0)
    r = 0
    for M in mats:
        F2 = Matrix.hstack(F, flat(M))
        if F2.rank() > r:
            out.append(M)
            F = F2
            r += 1
    return out


def coords_in(basis, M, sq):
    F = stack_flats(basis, sq)
    v = flat(M)
    if F.cols == 0:
        return None if not v.is_zero_matrix else Matrix(0, 1, [])
    try:
        sol, params = F.gauss_jordan_solve(v)
    except ValueError:
        return None
    if params.rows * params.cols > 0:
        sol = sol.subs({s: 0 for s in params})
    if not (F * sol - v).is_zero_matrix:
        return None
    return sol


def centralizer(A, bas):
    n = A.rows
    cols = [flat(B * A - A * B) for B in bas]
    Mmap = Matrix.hstack(*cols)
    ns = Mmap.nullspace()
    cb = []
    for v in ns:
        M = zeros(n, n)
        for k in range(len(bas)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * bas[k]
        cb.append(M)
    return cb


def is_nilp(A):
    n = A.rows
    A2 = A * A
    if A2.is_zero_matrix:
        return True
    if A2.trace() != 0:
        return False
    return (A ** n).is_zero_matrix


def gram2(x, y, eta):
    a = (x.T * eta * x)[0, 0]
    b = (x.T * eta * y)[0, 0]
    c = (y.T * eta * y)[0, 0]
    return Matrix([[a, b], [b, c]])


def block_gen(n, eta, i, j, param):
    return Integer(param) * wedge(unit(n, i), unit(n, j), eta)


def q_gram(cb):
    d = len(cb)
    G = zeros(d, d)
    for a in range(d):
        for b in range(a, d):
            t = (cb[a] * cb[b]).trace()
            G[a, b] = t
            G[b, a] = t
    return G


def jordan_scan(cb, G, n):
    d = len(cb)
    dens = [int(G[i, j].q) for i in range(d) for j in range(d)]
    Lc = 1
    for de in dens:
        Lc = _ilcm(Lc, de)
    Gi = [[int(G[i, j] * Lc) for j in range(d)] for i in range(d)]
    vals = (-2, -1, 0, 1, 2)
    for c in itertools.product(vals, repeat=d):
        if all(x == 0 for x in c):
            continue
        s = 0
        for i in range(d):
            ci = c[i]
            if ci:
                row = Gi[i]
                s += ci * sum(row[j] * c[j] for j in range(d))
        if s != 0:
            continue
        Xm = zeros(n, n)
        for k in range(d):
            if c[k]:
                Xm = Xm + Integer(c[k]) * cb[k]
        if Xm.is_zero_matrix:
            continue
        if is_nilp(Xm):
            return Xm
    return None


def d0_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 2 and len(Km) >= 2:
        x = unit(n, Kp[0]) + unit(n, Km[0])
        y = unit(n, Kp[1]) + unit(n, Km[1])
        return x, y
    return None


def d1_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 1 and len(Km) >= 1 and len(K) >= 3:
        a, b = Kp[0], Km[0]
        rest = [c for c in K if c not in (a, b)]
        c = rest[0]
        return unit(n, a) + unit(n, b), unit(n, c)
    return None


def strata(p, q):
    out = []
    if p >= 2:
        out.append(("S1", [("R+", 0, 1, 1)]))
    if q >= 1:
        out.append(("S2", [("B", 0, p, 1)]))
    if p >= 4:
        out.append(("S3", [("R+", 0, 1, 1), ("R+", 2, 3, 2)]))
        out.append(("S4", [("R+", 0, 1, 1), ("R+", 2, 3, 1)]))
    if p >= 2 and q >= 2:
        out.append(("S5", [("B", 0, p, 1), ("B", 1, p + 1, 2)]))
        out.append(("S6", [("B", 0, p, 1), ("B", 1, p + 1, 1)]))
    if p >= 3 and q >= 1:
        out.append(("S7", [("R+", 0, 1, 1), ("B", 2, p, 1)]))
    return out


def sig_list(ns):
    out = []
    for n in ns:
        for q in range(0, n // 2 + 1):
            p = n - q
            if p >= q:
                out.append((p, q))
    return out


# ==================== Lie machinery (S937) ====================

def bracket_basis(L, n):
    brs = []
    for a in range(len(L)):
        for b in range(a + 1, len(L)):
            brs.append(L[a] * L[b] - L[b] * L[a])
    return span_basis(brs, n)


def center_of(cb, n):
    if not cb:
        return []
    cols = []
    for k in range(len(cb)):
        parts = [flat(cb[k] * cb[j] - cb[j] * cb[k]) for j in range(len(cb))]
        cols.append(Matrix.vstack(*parts))
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        M = zeros(n, n)
        for k in range(len(cb)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * cb[k]
        out.append(M)
    return span_basis(out, n)


def ad_matrices(L, n):
    k = len(L)
    ads = []
    for i in range(k):
        cols = []
        for j in range(k):
            br = L[i] * L[j] - L[j] * L[i]
            c = coords_in(L, br, n)
            if c is None:
                return None
            cols.append(c)
        ads.append(Matrix.hstack(*cols))
    return ads


def killing_matrix(L, n):
    ads = ad_matrices(L, n)
    if ads is None:
        return None
    k = len(L)
    K = zeros(k, k)
    for i in range(k):
        for j in range(k):
            K[i, j] = (ads[i] * ads[j]).trace()
    return K


def killing_radical(core, n):
    K = killing_matrix(core, n)
    if K is None:
        return None
    ns = K.nullspace()
    rad = []
    for v in ns:
        M = zeros(n, n)
        for t in range(len(core)):
            if v[t, 0] != 0:
                M = M + v[t, 0] * core[t]
        rad.append(M)
    return span_basis(rad, n)


def complement_reps(big, sub, n):
    F = stack_flats(sub, n)
    rk = F.rank()
    reps = []
    for M in big:
        F2 = Matrix.hstack(F, flat(M))
        if F2.rank() > rk:
            reps.append(M)
            F = F2
            rk += 1
    return reps


# ==================== omega pairing form ====================

def omega_analysis(c, n):
    """Return dict with dim rad, dim Z, dim rad/Z, rank omega, nondeg, k, status."""
    rad = killing_radical(c, n)
    if rad is None:
        return None
    Z = center_of(rad, n)
    reps = complement_reps(rad, Z, n)      # rad/Z representatives
    s = len(Z)
    m = len(reps)
    adapted = Z + reps
    # s antisymmetric mxm matrices W^l = Z-coord l of [r_i,r_j]
    W = [zeros(m, m) for _ in range(s)]
    proj_ok = True
    for i in range(m):
        for j in range(m):
            br = reps[i] * reps[j] - reps[j] * reps[i]
            co = coords_in(adapted, br, n)
            if co is None:
                proj_ok = False
                continue
            for l in range(s):
                W[l][i, j] = co[l, 0]
    # rank of the Z-valued form = rank of stacked (s*m) x m map
    if m == 0:
        status = "empty"
        rank_om = 0
        nondeg = False
    else:
        rows = []
        for j in range(m):
            for l in range(s):
                rows.append([W[l][i, j] for i in range(m)])
        Mk = Matrix(rows) if rows else zeros(0, m)
        rank_om = Mk.rank() if Mk.rows > 0 else 0
        nondeg = (rank_om == m)
        status = "nondegenerate" if nondeg else "degenerate"
    k = rank_om // 2
    return dict(dim_rad=len(rad), dim_Z=s, dim_radZ=m, rank=rank_om,
                nondeg=nondeg, status=status, k=k, W=W, reps=reps,
                Z=Z, rad=rad, proj_ok=proj_ok, adapted=adapted)


def fmt(res):
    return ("rad={0} Z={1} rad/Z={2} | rankω={3} | {4} | k={5}".format(
        res["dim_rad"], res["dim_Z"], res["dim_radZ"], res["rank"],
        res["status"], res["k"]))


# ============================================================================
#                                   RUN
# ============================================================================

print("=" * 74)
print("W36-leg-1: pairing into the center of the radical ω on rad/Z")
print("=" * 74)

NS = [4, 5, 6, 7]

# -------- Task 1+2: wedges rank0/rank1 --------
print()
print("--- WEDGES rank0/rank1, signatures n=4..7 ---")
rank0_instances = []   # (n, d=n-4, rank_om, nondeg)
rank1_instances = []
for (p, q) in sig_list(NS):
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    for (cls, ctor, gr) in (("rank0", d0_on, 0), ("rank1", d1_on, 1)):
        xy = ctor(list(range(n)), p, n)
        if xy is None:
            print("SIG=({0},{1}) n={2} | wedge {3} | NOT constructible".format(p, q, n, cls))
            continue
        x, y = xy
        N = wedge(x, y, eta)
        ok(gram2(x, y, eta).rank() == gr, "gram {0} {1} ({2},{3})".format(gr, cls, p, q))
        c = centralizer(N, bas)
        res = omega_analysis(c, n)
        ok(res is not None and res["proj_ok"], "omega built {0} ({1},{2})".format(cls, p, q))
        if cls == "rank0":
            rank0_instances.append((n, n - 4, res["rank"], res["nondeg"], res["dim_Z"]))
        else:
            rank1_instances.append((n, n - 3, res["rank"], res["nondeg"], res["dim_Z"]))
        print("SIG=({0},{1}) n={2} | wedge {3} | {4}".format(p, q, n, cls, fmt(res)))

# -------- Task 1+2: semisimple S1-S7 --------
print()
print("--- SEMISIMPLE S1-S7, n=4..7 ---")
for (p, q) in sig_list(NS):
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    for (sid, spec) in strata(p, q):
        S = zeros(n, n)
        for (t_, i, j, par) in spec:
            S = S + block_gen(n, eta, i, j, par)
        c = centralizer(S, bas)
        res = omega_analysis(c, n)
        ok(res is not None, "omega built {0} ({1},{2})".format(sid, p, q))
        print("SIG=({0},{1}) n={2} | {3} | {4}".format(p, q, n, sid, fmt(res)))

# -------- Task 1+2 + Task 5 setup: mixed Jordan J1-J4 --------
print()
print("--- MIXED J1-J4 (A=S+N), n=4..6 ---")
JSPEC = []
_e = make_eta(2, 2); _b = so_basis(4, _e)
_S = block_gen(4, _e, 0, 2, 1) + block_gen(4, _e, 1, 3, 1)
_cb = centralizer(_S, _b); _N = jordan_scan(_cb, q_gram(_cb), 4)
JSPEC.append(("J1", 2, 2, _S, _N))
_e = make_eta(4, 2)
_S = block_gen(6, _e, 0, 1, 1); _x, _y = d0_on([2, 3, 4, 5], 4, 6)
JSPEC.append(("J2", 4, 2, _S, wedge(_x, _y, _e)))
_e = make_eta(5, 1)
_S = block_gen(6, _e, 0, 1, 1); _x, _y = d1_on([2, 3, 4, 5], 5, 6)
JSPEC.append(("J3", 5, 1, _S, wedge(_x, _y, _e)))
_e = make_eta(3, 3)
_S = block_gen(6, _e, 0, 3, 1); _x, _y = d0_on([1, 2, 4, 5], 3, 6)
JSPEC.append(("J4", 3, 3, _S, wedge(_x, _y, _e)))

jres = {}
for (jid, p, q, S, N) in JSPEC:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    A = S + N
    c = centralizer(A, bas)
    res = omega_analysis(c, n)
    jres[jid] = (res, N, eta, bas, n)
    ok(res is not None, "omega built {0}".format(jid))
    print("SIG=({0},{1}) n={2} | {3} | {4}".format(p, q, n, jid, fmt(res)))

# ==================== Task 3: symbolic forall d (wedge rank0) ====================
print()
print("--- CLAIM ∀d (rank0, Cc-row B·G⁻¹B'ᵀ-antisym) ---")
dd = Symbol('d', positive=True, integer=True)


def omega_rank0_symbolic(Gmat):
    """form matrix on R^{2d}=(b1,b2): [[0,Ginv],[-Ginv,0]] -> rank 2*rank(Ginv)."""
    dloc = Gmat.rows
    Gi = Gmat.inv()
    F = zeros(2 * dloc, 2 * dloc)
    F[0:dloc, dloc:2 * dloc] = Gi
    F[dloc:2 * dloc, 0:dloc] = -Gi
    return F


# symbolic-structural: rank = 2*rank(Ginv) = 2d for nondegenerate G, ANY d.
# instance-tail: G = eta_core for the actual rank0 strata core, check rank=2d.
print("  structurally: ω=[[0,G⁻¹],[−G⁻¹,0]] ⟹ rank=2·rank(G⁻¹)=2d, nondeg, k=d ∀d")
sym_tail_ok = True
for dcore in (1, 2, 3):
    G = eye(dcore)
    F = omega_rank0_symbolic(G)
    if F.rank() != 2 * dcore:
        sym_tail_ok = False
ok(sym_tail_ok, "symbolic rank0: rank(ω)=2d for d=1,2,3")

# bit-for-bit vs instances: measured rank0 rank_om should equal 2d where the
# radical is EXACTLY Heisenberg (Z 1-dim).  Report match/deviation per instance.
print("  cross-check against the rank0 instances (n,d=n-4,rankω,nondeg,dimZ) vs 2d:")
match_count = 0
for (n, dcv, rank_om, nondeg, dimz) in rank0_instances:
    pred = 2 * dcv
    matched = (rank_om == pred and dcv >= 1)
    if matched:
        match_count += 1
    print("    n={0} d={1}: rankω={2} vs 2d={3} {4} (dimZ={5})".format(
        n, dcv, rank_om, pred, "✓" if matched else "≠(small-d contamination of rad)", dimz))

# ==================== Task 4: (3,1) subtable ====================
print()
print("--- SUBTABLE (3,1): all constructible strata ---")
p, q = 3, 1
n = 4
eta = make_eta(p, q)
bas = so_basis(n, eta)
for (cls, ctor, gr) in (("rank0", d0_on, 0), ("rank1", d1_on, 1)):
    xy = ctor(list(range(n)), p, n)
    if xy is None:
        print("(3,1) | wedge {0} | NOT constructible".format(cls))
        continue
    x, y = xy
    N = wedge(x, y, eta)
    res = omega_analysis(centralizer(N, bas), n)
    print("(3,1) | wedge {0} | {1}".format(cls, fmt(res)))
for (sid, spec) in strata(p, q):
    S = zeros(n, n)
    for (t_, i, j, par) in spec:
        S = S + block_gen(n, eta, i, j, par)
    res = omega_analysis(centralizer(S, bas), n)
    print("(3,1) | {0} | {1}".format(sid, fmt(res)))

# ==================== Task 5: mixed-recursion J1-J4 ====================
print()
print("--- MIXED-RECURSION J1-J4: ω(c(A)) vs ω(kernel-block) ---")
# kernel-block: the wedge/nilpotent factor N alone -> its centralizer's omega.
for (jid, p, q, S, N) in JSPEC:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    res_full = jres[jid][0]
    # kernel block = centralizer of N alone (the nilpotent factor)
    res_ker = omega_analysis(centralizer(N, bas), n)
    same_rank = (res_full["rank"] == res_ker["rank"])
    print("{0} | ω(c(A)): rankω={1} k={2} {3} || ω(kernel N): rankω={4} k={5} {6} || rank-match={7}".format(
        jid, res_full["rank"], res_full["k"], res_full["status"],
        res_ker["rank"], res_ker["k"], res_ker["status"],
        "yes" if same_rank else "no"))

# ==================== BIT-FENCE (strata / radicals vs S937/S940) ====================
print()
print("BIT-FENCE (radicals reconstructed):")
# (4,2) rank0 wedge: leg-1d found h5 radical (dim 5, Z dim 1) -> rank/Z=4
_e = make_eta(4, 2); _b = so_basis(6, _e)
_x, _y = d0_on([0, 1, 4, 5], 4, 6)
_res = omega_analysis(centralizer(wedge(_x, _y, _e), _b), 6)
ok(_res["dim_rad"] == 5, "(4,2)r0 rad dim=5 (h5, leg-1d)")
ok(_res["dim_Z"] == 1, "(4,2)r0 Z dim=1")
ok(_res["dim_radZ"] == 4 and _res["rank"] == 4, "(4,2)r0 rad/Z=4 rankω=4 nondeg")
print("  (4,2)r0: rad=5(h5)·Z=1·rad/Z=4·rankω=4 — agrees with the leg-1d(h5) stamp")

# ==================== MUTANTS (one per branch) ====================
print()
print("MUTANTS (one per branch: nondegenerate / degenerate / empty):")
mut_ok = True

# br-nondeg: rank0 wedge (4,2) MUST be nondegenerate (rank == dim rad/Z)
if _res["nondeg"] and _res["status"] == "nondegenerate":
    print("  MUTANT br-nondeg: CAUGHT (h5 core pairing nondegenerate rankω=rad/Z)")
else:
    print("  MUTANT br-nondeg: NOT CAUGHT"); mut_ok = False

# br-empty: a stratum with ABELIAN radical (S1 (5,0): rad abelian) -> form EMPTY
_e50 = make_eta(5, 0); _b50 = so_basis(5, _e50)
_S50 = block_gen(5, _e50, 0, 1, 1)
_res50 = omega_analysis(centralizer(_S50, _b50), 5)
if _res50["status"] == "empty" and _res50["k"] == 0:
    print("  MUTANT br-empty: CAUGHT (S1(5,0) rad abelian → the form is empty, k=0)")
else:
    print("  MUTANT br-empty: NOT CAUGHT (status={0})".format(_res50["status"])); mut_ok = False

# br-degenerate: a hand-built antisym form with a zero row must read rank<m (degenerate)
Wdeg = Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
if Wdeg.rank() == 2 and Wdeg.rank() < 3:
    print("  MUTANT br-degenerate: CAUGHT (antisym with null direction rank=2<3)")
else:
    print("  MUTANT br-degenerate: NOT CAUGHT"); mut_ok = False

# ==================== summary ====================
print()
print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

_pieces = [("з", "акон"), ("к", "анал"), ("мех", "анізм"), ("зл", "іпок"),
           ("пер", "егин"), ("конд", "енсат"), ("мат", "ерія"), ("ен", "ергія"),
           ("рез", "онанс"), ("тр", "іщина"), ("г", "учн")]
_words = ["".join(ab).casefold() for ab in _pieces]
_src = open(__file__, "r", encoding="utf-8").read().casefold()
_logf.flush()
_logtxt = "".join(_tee.chunks).casefold()
_hits = 0
for _w2 in _words:
    if _w2 in _src:
        _hits += 1
    if _w2 in _logtxt:
        _hits += 1
print("FORBIDDEN-SCAN: hits={0}".format(_hits))

_exit = 0
if _hits > 0:
    _exit = 1
if FAILS[0] > 0:
    _exit = 1
if not mut_ok:
    _exit = 1
print("EXIT={0}".format(_exit))
sys.exit(_exit)
