# -*- coding: utf-8 -*-
# DIM: na (a measurement table of centralizers in so(p,q); handles 0 — the probe postulates nothing).
import sys
sys.stdout.reconfigure(encoding='utf-8')

from sympy import Matrix, Integer, zeros, eye, diag, symbols, Poly, gcd

X = symbols('x')


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


# ---------- exact linear-algebra primitives ----------

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


def nilp_check(A):
    return (A ** A.rows).is_zero_matrix


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
            expr = X ** d
            for k in range(d):
                expr = expr - sol[k, 0] * X ** k
            return Poly(expr, X)
    return Poly(X ** n, X)  # unreachable by Cayley-Hamilton


def squarefree(mp):
    me = mp.as_expr()
    g = gcd(me, me.diff(X))
    return Poly(g, X).degree() == 0


def gram(x, y, eta):
    a = (x.T * eta * x)[0, 0]
    b = (x.T * eta * y)[0, 0]
    c = (y.T * eta * y)[0, 0]
    return Matrix([[a, b], [b, c]])


def eig_str(A):
    ev = A.eigenvals()
    items = sorted(ev.items(), key=lambda t: str(t[0]).replace(" ", ""))
    return "[" + ",".join("(" + str(e).replace(" ", "") + "," + str(m) + ")" for e, m in items) + "]"


# ---------- strata builders ----------

def block(n, eta, i, j, param):
    return Integer(param) * wedge(unit(n, i), unit(n, j), eta)


def single_block_specs(p, q):
    out = []
    if p >= 2:
        out.append([("R+", 0, 1, 1)])
    if q >= 2:
        out.append([("R-", p, p + 1, 1)])
    if p >= 1 and q >= 1:
        out.append([("B", 0, p, 1)])
    return out


def two_block_specs(p, q):
    out = []
    for pars in ((1, 1), (1, 2)):
        if p >= 4:
            out.append([("R+", 0, 1, pars[0]), ("R+", 2, 3, pars[1])])
        if p >= 2 and q >= 2:
            out.append([("R+", 0, 1, pars[0]), ("R-", p, p + 1, pars[1])])
        if q >= 4:
            out.append([("R-", p, p + 1, pars[0]), ("R-", p + 2, p + 3, pars[1])])
        if p >= 2 and q >= 2:
            out.append([("B", 0, p, pars[0]), ("B", 1, p + 1, pars[1])])
    if p >= 3 and q >= 1:
        out.append([("R+", 0, 1, 1), ("B", 2, p, 1)])
    if p >= 1 and q >= 3:
        out.append([("R-", p, p + 1, 1), ("B", 0, p + 2, 1)])
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


def d0_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 2 and len(Km) >= 2:
        x = unit(n, Kp[0]) + unit(n, Km[0])
        y = unit(n, Kp[1]) + unit(n, Km[1])
        return x, y
    return None


# ---------- measurement ----------

ROWCOUNT = {"pure-block": 0, "pure-wedge": 0, "mixed": 0, "NONE": 0}


def prow(p, q, typ, Sspec, Nclass, Kstr, m):
    ROWCOUNT[typ] += 1
    if m is None:
        tail = "dimc=- | dimcc=- | abelian=- | eig=- | dimker=- | A2zero=- | minpolydeg=-"
    else:
        tail = ("dimc={dimc} | dimcc={dimcc} | abelian={abelian} | eig={eig} | "
                "dimker={dimker} | A2zero={A2zero} | minpolydeg={mpd}").format(**m)
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
    if kind == "pure-block":
        ok(squarefree(mp), "block A diagonalizable over C " + tag)
    elif kind == "pure-wedge":
        ok(nilp_check(A), "wedge A nilpotent " + tag)
    elif kind == "mixed":
        ok(squarefree(minpoly_poly(Sm)), "S diagonalizable " + tag)
        ok(nilp_check(Nm), "N nilpotent " + tag)
        ok((Sm * Nm - Nm * Sm).is_zero_matrix, "[S,N]=0 " + tag)
    return dict(dimc=len(cb), dimcc=dimcc, abelian=1 if dimcc == 0 else 0,
                eig=eig_str(A), dimker=n - A.rank(),
                A2zero=1 if (A * A).is_zero_matrix else 0, mpd=mp.degree())


def wedge_class_asserts(x, y, eta, cls, tag):
    ok(Matrix.hstack(x, y).rank() == 2, "x,y linearly independent " + tag)
    g = gram(x, y, eta)
    if cls == "D1":
        ok(g.rank() == 1, "Gram rank 1 " + tag)
    else:
        ok(g.is_zero_matrix, "Gram rank 0 " + tag)


# ---------- main enumeration ----------

def sig_list():
    out = []
    for n in (4, 5, 6):
        for q in range(0, n // 2 + 1):
            p = n - q
            if p >= q:
                out.append((p, q))
    return out


for (p, q) in sig_list():
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    ok(len(bas) == n * (n - 1) // 2, "basis size so({0},{1})".format(p, q))
    ok(all(is_so(B, eta) for B in bas), "basis membership so({0},{1})".format(p, q))
    ok(Matrix.hstack(*[flat(B) for B in bas]).rank() == n * (n - 1) // 2,
       "basis independent so({0},{1})".format(p, q))

    # 1. pure-block strata
    for spec in single_block_specs(p, q) + two_block_specs(p, q):
        A = zeros(n, n)
        for (t, i, j, par) in spec:
            A = A + block(n, eta, i, j, par)
        sstr = ",".join("{0}({1},{2})x{3}".format(t, i, j, par) for (t, i, j, par) in spec)
        m = measure(A, p, q, eta, bas, "pure-block")
        prow(p, q, "pure-block", sstr, "-", "-", m)

    # 2. pure-wedge strata
    full = list(range(n))
    for cls, builder in (("D1", d1_on), ("D0", d0_on)):
        pr = builder(full, p, n)
        if pr is None:
            prow(p, q, "NONE", "-", cls, "-", None)
            continue
        x, y = pr
        W = wedge(x, y, eta)
        tag = "pure-wedge " + cls + " (" + str(p) + "," + str(q) + ")"
        ok(is_so(W, eta), "W(x,y) in so " + tag)
        wedge_class_asserts(x, y, eta, cls, tag)
        m = measure(W, p, q, eta, bas, "pure-wedge")
        prow(p, q, "pure-wedge", "-", cls, "-", m)

    # 3. mixed strata
    stypes = []
    if p >= 2:
        stypes.append(("R+", 0, 1))
    if q >= 2:
        stypes.append(("R-", p, p + 1))
    if p >= 1 and q >= 1:
        stypes.append(("B", 0, p))
    for (t, i, j) in stypes:
        Sm = block(n, eta, i, j, 1)
        K = [a for a in range(n) if a not in (i, j)]
        p0 = len([a for a in K if a < p])
        q0 = len(K) - p0
        Kstr = "({0},{1})".format(p0, q0)
        sstr = "{0}({1},{2})x1".format(t, i, j)
        for cls, builder in (("D1", d1_on), ("D0", d0_on)):
            pr = builder(K, p, n)
            if pr is None:
                prow(p, q, "NONE", sstr, cls, Kstr, None)
                continue
            x, y = pr
            Nm = wedge(x, y, eta)
            tag = "mixed " + t + "+" + cls + " (" + str(p) + "," + str(q) + ")"
            ok(is_so(Nm, eta), "N in so " + tag)
            wedge_class_asserts(x, y, eta, cls, tag)
            A = Sm + Nm
            m = measure(A, p, q, eta, bas, "mixed", Sm, Nm)
            prow(p, q, "mixed", sstr, cls, Kstr, m)

# ---------- mutants (real code paths; each must be CAUGHT) ----------

print("--- mutants ---")
mut_ok = True

# m1: inject a non-commuting matrix into a computed c-basis; commute-check must reject
p, q = 2, 2
n = p + q
eta = make_eta(p, q)
bas = so_basis(n, eta)
A = block(n, eta, 0, 2, 1)
cb = centralizer(A, bas)
bad = None
for B in bas:
    if not (B * A - A * B).is_zero_matrix:
        bad = B
        break
if bad is None:
    print("MUTANT m1: NOT CAUGHT")
    mut_ok = False
else:
    tampered = cb + [bad]
    if not commutes_all(tampered, A):
        print("MUTANT m1: CAUGHT")
    else:
        print("MUTANT m1: NOT CAUGHT")
        mut_ok = False

# m2: membership check with one-entry sign-flipped eta on a legitimate generator must fail
G = wedge(unit(n, 0), unit(n, 2), eta)
ok(is_so(G, eta), "m2 baseline: generator is legitimate under true eta")
eta_bad = eta.copy()
eta_bad[0, 0] = -eta_bad[0, 0]
if not is_so(G, eta_bad):
    print("MUTANT m2: CAUGHT")
else:
    print("MUTANT m2: NOT CAUGHT")
    mut_ok = False

# m3: rank-2-Gram wedge (class D2) fed to the nilpotency check must fail it
x = unit(n, 0)
y = unit(n, 1)
W2 = wedge(x, y, eta)
g2 = gram(x, y, eta)
ok(g2.rank() == 2, "m3 baseline: Gram rank is 2 (class D2)")
if not nilp_check(W2):
    print("MUTANT m3: CAUGHT")
else:
    print("MUTANT m3: NOT CAUGHT")
    mut_ok = False

# ---------- summary ----------

total_rows = sum(ROWCOUNT.values())
print("SUMMARY: rows={0} | pure-block={1} | pure-wedge={2} | mixed={3} | NONE={4}".format(
    total_rows, ROWCOUNT["pure-block"], ROWCOUNT["pure-wedge"], ROWCOUNT["mixed"], ROWCOUNT["NONE"]))
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
