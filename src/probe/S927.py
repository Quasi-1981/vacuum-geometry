# -*- coding: utf-8 -*-
# DIM: na (a measurement table of centralizers in so(p,q); handles 0 — the probe postulates nothing).
import sys
sys.stdout.reconfigure(encoding='utf-8')

import itertools
from math import lcm as _ilcm

from sympy import Matrix, Integer, zeros, eye, diag, symbols, Poly, gcd

X_ = symbols('x')


class Tee:
    def __init__(self, real):
        self.real = real
        self.chunks = []

    def write(self, s):
        self.real.write(s)
        self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()


_tee = Tee(sys.stdout)
sys.stdout = _tee

ASSERT_PASS = [0]
FAILS = [0]


def ok(cond, msg):
    if cond:
        ASSERT_PASS[0] += 1
    else:
        FAILS[0] += 1
        print("ASSERT-FAIL: " + msg)


# ---------- exact linear-algebra primitives (S925 conventions) ----------

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


def commutes_all(cb, A):
    return all((M * A - A * M).is_zero_matrix for M in cb)


def is_nilp(A):
    n = A.rows
    A2 = A * A
    if A2.is_zero_matrix:
        return True
    if A2.trace() != 0:
        return False
    A4 = A2 * A2
    if A4.is_zero_matrix:
        return True
    if A4.trace() != 0:
        return False
    return (A ** n).is_zero_matrix


def derived_and_closed(cb, n):
    if not cb:
        return 0, True
    fc = Matrix.hstack(*[flat(M) for M in cb])
    rc = fc.rank()
    brs = []
    for a in range(len(cb)):
        for b in range(a + 1, len(cb)):
            brs.append(flat(cb[a] * cb[b] - cb[b] * cb[a]))
    if not brs:
        return 0, True
    fb = Matrix.hstack(*brs)
    dimcc = fb.rank()
    closed = (Matrix.hstack(fc, fb).rank() == rc)
    return dimcc, closed


def minpoly_poly(A):
    n = A.rows
    pows = [eye(n)]
    for _k in range(n):
        pows.append(pows[-1] * A)
    flats = [flat(P) for P in pows]
    for d in range(1, n + 1):
        Mstack = Matrix.hstack(*flats[:d])
        v = flats[d]
        if Matrix.hstack(Mstack, v).rank() == Mstack.rank():
            sol, params = Mstack.gauss_jordan_solve(v)
            if params.rows * params.cols > 0:
                sol = sol.subs({s: 0 for s in params})
            expr = X_ ** d
            for k in range(d):
                expr = expr - sol[k, 0] * X_ ** k
            return Poly(expr, X_)
    return Poly(X_ ** n, X_)  # unreachable by Cayley-Hamilton


def squarefree(mp):
    me = mp.as_expr()
    g = gcd(me, me.diff(X_))
    return Poly(g, X_).degree() == 0


def gram(x, y, eta):
    a = (x.T * eta * x)[0, 0]
    b = (x.T * eta * y)[0, 0]
    c = (y.T * eta * y)[0, 0]
    return Matrix([[a, b], [b, c]])


def eig_str(A):
    ev = A.eigenvals()
    items = sorted(ev.items(), key=lambda t: str(t[0]).replace(" ", ""))
    return "[" + ",".join("(" + str(e).replace(" ", "") + "," + str(m) + ")" for e, m in items) + "]"


def block(n, eta, i, j, param):
    return Integer(param) * wedge(unit(n, i), unit(n, j), eta)


# ---------- quadratic form Q(X)=tr(X*X) on a centralizer basis ----------

def q_gram(cb):
    d = len(cb)
    G = zeros(d, d)
    for a in range(d):
        for b in range(a, d):
            t = (cb[a] * cb[b]).trace()
            G[a, b] = t
            G[b, a] = t
    return G


def cong_signature(Gin):
    # exact rational congruence diagonalization (pivoted completion of squares)
    G = Gin.copy()
    pos = neg = zer = 0
    while G.rows > 0:
        m = G.rows
        k = None
        for i in range(m):
            if G[i, i] != 0:
                k = i
                break
        if k is None:
            pair = None
            for i in range(m):
                for j in range(i + 1, m):
                    if G[i, j] != 0:
                        pair = (i, j)
                        break
                if pair is not None:
                    break
            if pair is None:
                zer += m
                break
            i, j = pair
            # congruence b_i <- b_i + b_j creates a nonzero diagonal entry
            G[i, :] = G[i, :] + G[j, :]
            G[:, i] = G[:, i] + G[:, j]
            continue
        dv = G[k, k]
        if dv > 0:
            pos += 1
        else:
            neg += 1
        v = G[:, k]
        G = G - (v * v.T) / dv
        G.row_del(k)
        G.col_del(k)
    return (pos, neg, zer)


# ---------- deterministic nilpotent scan in span of c(S) basis ----------

def jordan_scan(cb, G, n):
    d = len(cb)
    dens = []
    for i in range(d):
        for j in range(d):
            dens.append(int(G[i, j].q))
    L = 1
    for de in dens:
        L = _ilcm(L, de)
    Gi = [[int(G[i, j] * L) for j in range(d)] for i in range(d)]
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
        if s != 0:  # tr(X*X) != 0 => not nilpotent
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


# ---------- collision builders ----------

def collisions(p, q):
    out = []
    if p >= 4:
        out.append(("RR+", [("R+", 0, 1), ("R+", 2, 3)]))
    if q >= 4:
        out.append(("RR-", [("R-", p, p + 1), ("R-", p + 2, p + 3)]))
    if p >= 2 and q >= 2:
        out.append(("BB", [("B", 0, p), ("B", 1, p + 1)]))
    return out


def d1_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 1 and len(Km) >= 1 and len(K) >= 3:
        a, b = Kp[0], Km[0]
        rest = [c for c in K if c not in (a, b)]
        c = rest[0]
        return unit(n, a) + unit(n, b), unit(n, c)
    return None


# ---------- measurement ----------

ROWCOUNT = {"coll-pure": 0, "coll-jordan": 0, "coll-wedge": 0,
            "NONE": 0, "NONE-CERT": 0, "NONE-SCAN": 0}


def prow(p, q, typ, Sspec, Nclass, Kstr, m, qsig=None):
    ROWCOUNT[typ] += 1
    if m is None:
        tail = "dimc=- | dimcc=- | abelian=- | eig=- | dimker=- | A2zero=- | minpolydeg=-"
    else:
        tail = ("dimc={dimc} | dimcc={dimcc} | abelian={abelian} | eig={eig} | "
                "dimker={dimker} | A2zero={A2zero} | minpolydeg={mpd}").format(**m)
    if qsig is not None:
        tail = tail + " | Qsig=({0},{1},{2})".format(*qsig)
    print("SIG=({0},{1}) | TYPE={2} | S={3} | N={4} | K={5} | ".format(p, q, typ, Sspec, Nclass, Kstr) + tail)


def measure(A, p, q, eta, bas, kind, Sm=None, Nm=None):
    n = p + q
    tag = kind + " (" + str(p) + "," + str(q) + ")"
    ok(is_so(A, eta), "A in so " + tag)
    cb = centralizer(A, bas)
    ok(all(is_so(M, eta) for M in cb), "c-basis in so " + tag)
    ok(commutes_all(cb, A), "c-basis commutes with A " + tag)
    dimcc, closed = derived_and_closed(cb, n)
    ok(closed, "c closed under bracket " + tag)
    mp = minpoly_poly(A)
    if kind == "coll-pure":
        ok(squarefree(mp), "collision A diagonalizable over C " + tag)
    else:
        ok(squarefree(minpoly_poly(Sm)), "S diagonalizable " + tag)
        ok(is_nilp(Nm), "N nilpotent " + tag)
        ok((Sm * Nm - Nm * Sm).is_zero_matrix, "[S,N]=0 " + tag)
        ok(is_so(Nm, eta), "N in so " + tag)
    m = dict(dimc=len(cb), dimcc=dimcc, abelian=1 if dimcc == 0 else 0,
             eig=eig_str(A), dimker=n - A.rank(),
             A2zero=1 if (A * A).is_zero_matrix else 0, mpd=mp.degree())
    return m, cb


# ---------- main enumeration ----------

L1 = [(2, 2), (3, 2), (3, 3), (4, 2)]
L2 = [(5, 2), (6, 1), (4, 3)]

for (p, q) in L1 + L2:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    ok(len(bas) == n * (n - 1) // 2, "basis size so({0},{1})".format(p, q))
    ok(all(is_so(B, eta) for B in bas), "basis membership so({0},{1})".format(p, q))
    ok(Matrix.hstack(*[flat(B) for B in bas]).rank() == n * (n - 1) // 2,
       "basis independent so({0},{1})".format(p, q))

    for (cname, blist) in collisions(p, q):
        S = zeros(n, n)
        used = []
        for (t, i, j) in blist:
            S = S + block(n, eta, i, j, 1)
            used.extend([i, j])
        sstr = ",".join("{0}({1},{2})x1".format(t, i, j) for (t, i, j) in blist)

        # 1. coll-pure: A = S alone, plus Qsig of tr-form on c(S)
        m, cb = measure(S, p, q, eta, bas, "coll-pure")
        G = q_gram(cb)
        sig = cong_signature(G)
        ok(sig[0] + sig[1] + sig[2] == len(cb), "Qsig counts sum to dimc ({0},{1}) {2}".format(p, q, cname))
        ok(sig[0] + sig[1] == G.rank(), "Qsig rank cross-check ({0},{1}) {2}".format(p, q, cname))
        prow(p, q, "coll-pure", sstr, "-", "-", m, qsig=sig)

        # 2. coll-jordan
        definite = (sig[2] == 0 and len(cb) > 0 and (sig[0] == 0 or sig[1] == 0))
        if definite:
            prow(p, q, "NONE-CERT", sstr, "JORD", "-", None)
        else:
            Xn = jordan_scan(cb, G, n)
            if Xn is None:
                prow(p, q, "NONE-SCAN", sstr, "JORD", "-", None)
            else:
                tag = "coll-jordan " + cname + " ({0},{1})".format(p, q)
                ok(is_so(Xn, eta), "X in so " + tag)
                ok((Xn * S - S * Xn).is_zero_matrix, "[X,S]=0 " + tag)
                ok(is_nilp(Xn), "X nilpotent " + tag)
                ok(not Xn.is_zero_matrix, "X nonzero " + tag)
                A = S + Xn
                mj, _cbj = measure(A, p, q, eta, bas, "coll-jordan", Sm=S, Nm=Xn)
                prow(p, q, "coll-jordan", sstr, "JORD", "-", mj)

        # 3. coll-wedge on complement axes K
        K = [a for a in range(n) if a not in used]
        p0 = len([a for a in K if a < p])
        q0 = len(K) - p0
        Kstr = "({0},{1})".format(p0, q0)
        pr = d1_on(K, p, n)
        if pr is None:
            prow(p, q, "NONE", sstr, "D1", Kstr, None)
        else:
            x, y = pr
            N = wedge(x, y, eta)
            tag = "coll-wedge " + cname + "+D1 ({0},{1})".format(p, q)
            ok(Matrix.hstack(x, y).rank() == 2, "x,y linearly independent " + tag)
            ok(gram(x, y, eta).rank() == 1, "Gram rank 1 " + tag)
            ok(is_so(N, eta), "N in so " + tag)
            ok((S * N - N * S).is_zero_matrix, "[S,N]=0 " + tag)
            ok(is_nilp(N), "N nilpotent " + tag)
            A = S + N
            mw, _cbw = measure(A, p, q, eta, bas, "coll-wedge", Sm=S, Nm=N)
            prow(p, q, "coll-wedge", sstr, "D1", Kstr, mw)

# ---------- mutants (real code paths; each must be CAUGHT) ----------

print("--- mutants ---")
mut_ok = True

p, q = 2, 2
n = p + q
eta = make_eta(p, q)
bas = so_basis(n, eta)
Sm = block(n, eta, 0, 2, 1) + block(n, eta, 1, 3, 1)  # BB collision in (2,2)
cbm = centralizer(Sm, bas)

# m1: inject a non-commuting so-generator into a computed c-basis; commute-check must reject
bad = None
for B in bas:
    if not (B * Sm - Sm * B).is_zero_matrix:
        bad = B
        break
if bad is None:
    print("MUTANT m1: NOT CAUGHT")
    mut_ok = False
else:
    tampered = cbm + [bad]
    if not commutes_all(tampered, Sm):
        print("MUTANT m1: CAUGHT")
    else:
        print("MUTANT m1: NOT CAUGHT")
        mut_ok = False

# m2: membership check with one-entry sign-flipped eta must fail on a legitimate generator
Gm2 = wedge(unit(n, 0), unit(n, 2), eta)
ok(is_so(Gm2, eta), "m2 baseline: generator is legitimate under true eta")
eta_bad = eta.copy()
eta_bad[0, 0] = -eta_bad[0, 0]
if not is_so(Gm2, eta_bad):
    print("MUTANT m2: CAUGHT")
else:
    print("MUTANT m2: NOT CAUGHT")
    mut_ok = False

# m3: semisimple element (single B block) fed to the nilpotency check must fail it
Bsemi = block(n, eta, 0, 2, 1)
ok(squarefree(minpoly_poly(Bsemi)), "m3 baseline: single block is diagonalizable over C")
if not is_nilp(Bsemi):
    print("MUTANT m3: CAUGHT")
else:
    print("MUTANT m3: NOT CAUGHT")
    mut_ok = False

# m4: corrupt the Gram matrix of Q (zero one off-diagonal pair) on a case with known
# mixed signs; the congruence signature must differ from the uncorrupted one
Gm4 = q_gram(cbm)
sig0 = cong_signature(Gm4)
ok(sig0[0] > 0 and sig0[1] > 0, "m4 baseline: Q has mixed signs on BB (2,2)")
m4_caught = False
d4 = Gm4.rows
for i in range(d4):
    if m4_caught:
        break
    for j in range(i + 1, d4):
        if Gm4[i, j] != 0:
            Gb = Gm4.copy()
            Gb[i, j] = Integer(0)
            Gb[j, i] = Integer(0)
            if cong_signature(Gb) != sig0:
                m4_caught = True
                break
if not m4_caught:
    # Gram may be diagonal in the computed basis: apply a signature-preserving
    # basis mix b_i <- b_i + b_j first, then zero the created off-diagonal pair
    for i in range(d4):
        if m4_caught:
            break
        for j in range(d4):
            if i == j:
                continue
            cb2 = list(cbm)
            cb2[i] = cbm[i] + cbm[j]
            G2 = q_gram(cb2)
            if cong_signature(G2) != sig0:
                continue  # keep only genuinely congruent rebasings
            for a in range(d4):
                if m4_caught:
                    break
                for b in range(a + 1, d4):
                    if G2[a, b] != 0:
                        Gb = G2.copy()
                        Gb[a, b] = Integer(0)
                        Gb[b, a] = Integer(0)
                        if cong_signature(Gb) != sig0:
                            m4_caught = True
                            break
            if m4_caught:
                break
if m4_caught:
    print("MUTANT m4: CAUGHT")
else:
    print("MUTANT m4: NOT CAUGHT")
    mut_ok = False

# ---------- summary ----------

total_rows = sum(ROWCOUNT.values())
print("SUMMARY: rows={0} | coll-pure={1} | coll-jordan={2} | coll-wedge={3} | NONE={4} | NONE-CERT={5} | NONE-SCAN={6}".format(
    total_rows, ROWCOUNT["coll-pure"], ROWCOUNT["coll-jordan"], ROWCOUNT["coll-wedge"],
    ROWCOUNT["NONE"], ROWCOUNT["NONE-CERT"], ROWCOUNT["NONE-SCAN"]))
print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

# ---------- self-scan of source and produced log text for banned substrings ----------

_pieces = [("з", "акон"), ("к", "анал"), ("мех", "анізм"), ("зл", "іпок"),
           ("пер", "егин"), ("конд", "енсат"), ("мат", "ерія"), ("ен", "ергія"),
           ("рез", "онанс"), ("тр", "іщина"), ("г", "учн")]
_words = ["".join(ab).casefold() for ab in _pieces]
_src = open(__file__, "r", encoding="utf-8").read().casefold()
_log = "".join(_tee.chunks).casefold()
_hits = 0
for _w in _words:
    if _w in _src:
        _hits += 1
    if _w in _log:
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
