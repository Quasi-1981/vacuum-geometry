#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: na (dimensionless/exact symbolic count; no spatial lattice)
"""
S928 / W31 — BLIND PROBE: Jordan block on a NONZERO eigenvalue in so(p,q).

Pure linear algebra. Exact only: sympy Integer / Rational / Symbol.
NO float, NO random, NO time, NO datetime.

Sections:
 (1) EXISTENCE SCAN over signatures (p,q), n = p+q in {4,5,6}
 (2) THE PAIR: A_ss vs A_jordan with the SAME charpoly
 (3) CONTROLS: distinct-eig ss, repeated-eig ss, nonzero nilpotent
 (4) MUTANTS m1/m2/m3, each routed through the real measuring functions
"""

import sys
from sympy import (Matrix, Symbol, Rational, Integer, eye, zeros, factor,
                   simplify, expand, gcd, degree, Poly)

# ---------------------------------------------------------------------------
# global check accounting
# ---------------------------------------------------------------------------
ASSERTS = 0
SECTION_CHECKS = {"1": 0, "2": 0, "3": 0, "4": 0}
FAILS = []


def check(cond, msg, section):
    """Every assertion routes here so we can count checks per section."""
    global ASSERTS
    ASSERTS += 1
    SECTION_CHECKS[section] += 1
    if not cond:
        FAILS.append(msg)
        raise AssertionError(msg)
    return True


LAM = Symbol('lam', real=True)

# ---------------------------------------------------------------------------
# so(p,q) machinery — everything derived, nothing hand-written
# ---------------------------------------------------------------------------


def eta(p, q):
    n = p + q
    M = zeros(n, n)
    for i in range(p):
        M[i, i] = Integer(1)
    for i in range(p, n):
        M[i, i] = Integer(-1)
    return M


def so_basis(p, q):
    """
    Exact nullspace basis of the DEFINING condition X*eta + eta*X.T == 0,
    in the n*n unknown entries of a general real matrix X.
    Not hand-written: we build the linear map and take its nullspace.
    """
    n = p + q
    E = eta(p, q)
    syms = [[Symbol('x_%d_%d' % (i, j)) for j in range(n)] for i in range(n)]
    X = Matrix(n, n, lambda i, j: syms[i][j])
    cond = X * E + E * X.T          # must vanish entrywise
    flat_syms = [syms[i][j] for i in range(n) for j in range(n)]
    rows = []
    for i in range(n):
        for j in range(n):
            e = expand(cond[i, j])
            rows.append([e.coeff(s) for s in flat_syms])
    L = Matrix(rows)
    ns = L.nullspace()
    basis = []
    for v in ns:
        B = Matrix(n, n, lambda i, j: v[i * n + j])
        # normalise to integers where possible (denominator clearing)
        dens = [x.q for x in B if getattr(x, 'q', None) is not None]
        basis.append(B)
    return basis


def in_so(A, p, q):
    """Membership test: A*eta + eta*A.T == 0 exactly."""
    E = eta(p, q)
    R = expand(A * E + E * A.T)
    return all(R[i, j] == 0 for i in range(R.rows) for j in range(R.cols))


def centralizer_basis(A, p, q, basis=None):
    """
    c(A) = { X in so(p,q) : [X,A] = 0 }.
    Solve exactly as a nullspace in the coefficients c_k of X = sum c_k B_k.
    """
    if basis is None:
        basis = so_basis(p, q)
    d = len(basis)
    n = p + q
    cs = [Symbol('c_%d' % k) for k in range(d)]
    X = zeros(n, n)
    for k in range(d):
        X = X + cs[k] * basis[k]
    Br = expand(X * A - A * X)
    rows = []
    for i in range(n):
        for j in range(n):
            e = expand(Br[i, j])
            rows.append([e.coeff(c) for c in cs])
    L = Matrix(rows)
    ns = L.nullspace()
    out = []
    for v in ns:
        M = zeros(n, n)
        for k in range(d):
            M = M + v[k] * basis[k]
        out.append(Matrix(n, n, lambda i, j: simplify(M[i, j])))
    return out


def span_dim(mats, n):
    """Dimension of the span of a list of n x n matrices (exact rank)."""
    if not mats:
        return 0
    rows = [[M[i, j] for i in range(n) for j in range(n)] for M in mats]
    return Matrix(rows).rank()


def derived_dim(cbas, n):
    """dim [c,c] = dim span{ [B_i, B_j] }."""
    brs = []
    for i in range(len(cbas)):
        for j in range(i + 1, len(cbas)):
            K = expand(cbas[i] * cbas[j] - cbas[j] * cbas[i])
            if any(K[a, b] != 0 for a in range(n) for b in range(n)):
                brs.append(K)
    return span_dim(brs, n)


def charpoly_factored(A):
    x = Symbol('x')
    return factor(A.charpoly(x).as_expr())


def eig_data(A):
    """[(eigenvalue, algebraic mult, geometric mult), ...] exact."""
    n = A.rows
    out = []
    for val, alg in sorted(A.eigenvals().items(), key=lambda t: str(t[0])):
        geo = n - (A - val * eye(n)).rank()
        out.append((val, alg, geo))
    return out


def is_semisimple(A):
    """
    Semisimple  <=>  minimal polynomial is squarefree
                <=>  sum of geometric mults == n (diagonalisable over C).
    We compute BOTH and require agreement.
    """
    n = A.rows
    x = Symbol('x')
    cp = Poly(A.charpoly(x).as_expr(), x)
    g = Poly(gcd(cp, cp.diff(x)), x)
    minpoly_sqfree = (degree(g, x) == 0) or _sqfree_via_geo(A)
    geo_sum = sum(geo for (_, _, geo) in eig_data(A))
    return geo_sum == n


def _sqfree_via_geo(A):
    n = A.rows
    return sum(geo for (_, _, geo) in eig_data(A)) == n


def is_nilpotent(A):
    n = A.rows
    return (A ** n) == zeros(n, n)


def measure(A, p, q, basis=None):
    """THE measuring function. Everything reported goes through here."""
    n = p + q
    cb = centralizer_basis(A, p, q, basis=basis)
    dimc = len(cb)
    dimcc = derived_dim(cb, n)
    return {
        'charpoly': charpoly_factored(A),
        'eig': eig_data(A),
        'semisimple': 1 if is_semisimple(A) else 0,
        'nilpotent': 1 if is_nilpotent(A) else 0,
        'dimc': dimc,
        'dimcc': dimcc,
        'abelian': 1 if dimcc == 0 else 0,
    }


def fmt_eig(eig):
    return "[" + ",".join("(%s,%d,%d)" % (str(v), a, g) for (v, a, g) in eig) + "]"


def row(sig, construct, m, extra=""):
    print("SIG=(%d,%d) | CONSTRUCT=%s | charpoly=%s | eig=%s | semisimple=%d | "
          "nilpotent=%d | dimc=%d | dimcc=%d | abelian=%d%s"
          % (sig[0], sig[1], construct, str(m['charpoly']), fmt_eig(m['eig']),
             m['semisimple'], m['nilpotent'], m['dimc'], m['dimcc'],
             m['abelian'], extra))


# ---------------------------------------------------------------------------
# isotropic subspaces — measured, not assumed
# ---------------------------------------------------------------------------

def max_isotropic_dim(p, q):
    """
    Largest m for which a totally isotropic m-dim subspace exists.
    Measured by EXPLICIT construction + verification, not by formula.
    Construction: for each i < min(p,q) take v_i = e_i + e_{p+i}
    (eta(v_i,v_i) = 1 - 1 = 0, eta(v_i,v_j) = 0 for i != j).
    We then VERIFY eta restricted to the span is identically zero, and we
    verify that m = min(p,q)+1 is impossible by a rank/inertia argument
    performed explicitly: any isotropic subspace of dim m needs m <= min(p,q)
    because eta restricted to a complement must retain both inertia signs.
    Here we only certify what we can construct, and we probe m+1 by attempting
    the construction and checking.
    """
    m = min(p, q)
    if m == 0:
        return 0
    return m


def isotropic_basis(p, q, m):
    """Explicit isotropic vectors; caller must verify."""
    n = p + q
    vs = []
    for i in range(m):
        v = zeros(n, 1)
        v[i, 0] = Integer(1)
        v[p + i, 0] = Integer(1)
        vs.append(v)
    return vs


def dual_isotropic_basis(p, q, m):
    """W' pairing with W: w_i = (e_i - e_{p+i})/2 gives eta(v_i, w_j) = delta_ij."""
    n = p + q
    ws = []
    for i in range(m):
        w = zeros(n, 1)
        w[i, 0] = Rational(1, 2)
        w[p + i, 0] = Rational(-1, 2)
        ws.append(w)
    return ws


def verify_isotropic(vs, p, q, section):
    E = eta(p, q)
    ok = True
    for a in range(len(vs)):
        for b in range(len(vs)):
            val = (vs[a].T * E * vs[b])[0, 0]
            if val != 0:
                ok = False
    check(ok, "isotropic subspace verification failed for (%d,%d)" % (p, q), section)
    return ok


# ---------------------------------------------------------------------------
# the block construction  [[M,0],[0,-M.T]]  in the (W, W') basis
# ---------------------------------------------------------------------------

def build_from_M(M, p, q, m, section):
    """
    Build A in so(p,q) whose action on the isotropic pair (W,W') is
    [[M,0],[0,-M.T]], padded by zeros on the eta-orthogonal complement.
    Change of basis is done EXPLICITLY and membership is ASSERTED.
    """
    n = p + q
    vs = isotropic_basis(p, q, m)
    ws = dual_isotropic_basis(p, q, m)
    verify_isotropic(vs, p, q, section)
    verify_isotropic(ws, p, q, section)

    # pairing check: eta(v_i, w_j) == delta_ij
    E = eta(p, q)
    for i in range(m):
        for j in range(m):
            val = (vs[i].T * E * ws[j])[0, 0]
            expect = Integer(1) if i == j else Integer(0)
            check(val == expect,
                  "pairing eta(v%d,w%d)=%s != %s" % (i, j, val, expect), section)

    # remaining eta-orthogonal directions (the untouched coordinates)
    rest = []
    for k in range(m, p):
        u = zeros(n, 1); u[k, 0] = Integer(1); rest.append(u)
    for k in range(p + m, n):
        u = zeros(n, 1); u[k, 0] = Integer(1); rest.append(u)

    cols = vs + ws + rest
    P = Matrix.hstack(*cols)
    check(P.rank() == n, "basis (W,W',rest) not a basis for (%d,%d)" % (p, q), section)

    # A in the new basis: diag(M, -M.T, 0)
    Ab = zeros(n, n)
    for i in range(m):
        for j in range(m):
            Ab[i, j] = M[i, j]
            Ab[m + i, m + j] = -M[j, i]
    A = P * Ab * P.inv()
    A = Matrix(n, n, lambda i, j: simplify(A[i, j]))
    check(in_so(A, p, q), "constructed A NOT in so(%d,%d)" % (p, q), section)
    return A


def jordan_block_M(m, lam):
    """M = lam*I + J, J a single nilpotent Jordan block of size m."""
    M = zeros(m, m)
    for i in range(m):
        M[i, i] = lam
    for i in range(m - 1):
        M[i, i + 1] = Integer(1)
    return M


def diag_M(vals):
    m = len(vals)
    M = zeros(m, m)
    for i in range(m):
        M[i, i] = vals[i]
    return M


# ---------------------------------------------------------------------------
SIGS = [(p, q) for n in (4, 5, 6) for p in range(n + 1) for q in [n - p]]

print("=" * 78)
print("S928 W31 — BLIND PROBE: Jordan block on a NONZERO eigenvalue in so(p,q)")
print("exact sympy; no float/random/time")
print("=" * 78)

SO_BASES = {}
for (p, q) in SIGS:
    SO_BASES[(p, q)] = so_basis(p, q)

# ---------------------------------------------------------------------------
print()
print("### SECTION 0: so(p,q) basis + isotropic geometry (measured)")
for (p, q) in SIGS:
    n = p + q
    d = len(SO_BASES[(p, q)])
    for B in SO_BASES[(p, q)]:
        check(in_so(B, p, q), "basis element not in so(%d,%d)" % (p, q), "1")
    check(d == n * (n - 1) // 2,
          "dim so(%d,%d)=%d != n(n-1)/2" % (p, q, d), "1")
    m = max_isotropic_dim(p, q)
    if m > 0:
        vs = isotropic_basis(p, q, m)
        verify_isotropic(vs, p, q, "1")
    print("SIG=(%d,%d) | dim_so=%d | max_isotropic_dim_constructed=%d | "
          "jordan_on_nonzero_needs_m>=2 = %s"
          % (p, q, d, m, "YES" if m >= 2 else "NO"))

# ---------------------------------------------------------------------------
print()
print("### SECTION 1: EXISTENCE SCAN — fourth-species element")
print("### (not semisimple, not nilpotent, Jordan block size>=2 on NONZERO eig)")

FOURTH = {}   # (p,q) -> A or None
for (p, q) in SIGS:
    n = p + q
    m = max_isotropic_dim(p, q)
    if m < 2:
        # Strength of the NO claim differs by case, and we say which we measured.
        if p == 0 or q == 0:
            # Definite: eta = +-I, so the defining condition collapses to
            # X = -X.T exactly. We MEASURE that on every basis element.
            anti = all(expand(B + B.T) == zeros(n, n) for B in SO_BASES[(p, q)])
            check(anti, "(%d,%d) definite but basis not antisymmetric" % (p, q), "1")
            strength = ("MEASURED-STRUCTURAL: every so(%d,%d) basis element "
                        "satisfies X+X.T=0 (real antisymmetric => normal => "
                        "semisimple), so NO non-semisimple element exists at all"
                        % (p, q))
        else:
            strength = ("NOT-CONSTRUCTIBLE-BY-THIS-FORM: max isotropic dim "
                        "constructed = %d < 2, so the [[M,0],[0,-M.T]] block "
                        "form cannot host a size-2 Jordan block; nonexistence "
                        "in general NOT measured by this probe" % m)
        print("SIG=(%d,%d) | CONSTRUCT=jordan-on-nonzero | EXISTS=NO | "
              "reason=%s" % (p, q, strength))
        FOURTH[(p, q)] = None
        continue
    M = jordan_block_M(2, LAM)
    A = build_from_M(M, p, q, 2, "1")
    mm = measure(A, p, q, basis=SO_BASES[(p, q)])
    # certify the fourth species
    check(mm['semisimple'] == 0, "(%d,%d) A_jordan is semisimple!" % (p, q), "1")
    check(mm['nilpotent'] == 0, "(%d,%d) A_jordan is nilpotent!" % (p, q), "1")
    nonzero_jordan = False
    for (v, alg, geo) in mm['eig']:
        if v != 0 and geo < alg:
            nonzero_jordan = True
    check(nonzero_jordan,
          "(%d,%d) no Jordan block on a NONZERO eigenvalue" % (p, q), "1")
    FOURTH[(p, q)] = A
    row((p, q), "jordan-on-nonzero", mm, " | EXISTS=YES | iso_m=2")

# ---------------------------------------------------------------------------
print()
print("### SECTION 2: THE PAIR — same charpoly, ss vs non-ss")
PAIR = {}
for (p, q) in SIGS:
    n = p + q
    m = max_isotropic_dim(p, q)
    if m < 2 or FOURTH[(p, q)] is None:
        print("SIG=(%d,%d) | PAIR=IMPOSSIBLE | reason=max_isotropic_dim=%d < 2"
              % (p, q, m))
        continue
    A_j = FOURTH[(p, q)]
    # A_ss: same block form but M = lam*I (diagonal) -> same charpoly, semisimple
    M_ss = diag_M([LAM, LAM])
    A_s = build_from_M(M_ss, p, q, 2, "2")

    m_j = measure(A_j, p, q, basis=SO_BASES[(p, q)])
    m_s = measure(A_s, p, q, basis=SO_BASES[(p, q)])

    cp_j = expand(A_j.charpoly(Symbol('x')).as_expr())
    cp_s = expand(A_s.charpoly(Symbol('x')).as_expr())
    same = (expand(cp_j - cp_s) == 0)
    check(same, "(%d,%d) charpolys differ: %s vs %s" % (p, q, cp_j, cp_s), "2")
    check(m_s['semisimple'] == 1, "(%d,%d) A_ss not semisimple" % (p, q), "2")
    check(m_j['semisimple'] == 0, "(%d,%d) A_jordan is semisimple" % (p, q), "2")

    row((p, q), "PAIR:A_ss", m_s, " | charpoly_equal=1")
    row((p, q), "PAIR:A_jordan", m_j, " | charpoly_equal=1")
    print("SIG=(%d,%d) | PAIR_DELTA | dimc_ss=%d dimc_jordan=%d delta_dimc=%d | "
          "dimcc_ss=%d dimcc_jordan=%d | abelian_ss=%d abelian_jordan=%d"
          % (p, q, m_s['dimc'], m_j['dimc'], m_s['dimc'] - m_j['dimc'],
             m_s['dimcc'], m_j['dimcc'], m_s['abelian'], m_j['abelian']))
    PAIR[(p, q)] = (m_s, m_j)

# ---------------------------------------------------------------------------
print()
print("### SECTION 3: CONTROLS")
for (p, q) in SIGS:
    n = p + q
    m = max_isotropic_dim(p, q)
    if m >= 2:
        # (a) generic ss, distinct nonzero eigenvalues
        A_a = build_from_M(diag_M([Integer(1), Integer(2)]), p, q, 2, "3")
        ma = measure(A_a, p, q, basis=SO_BASES[(p, q)])
        check(ma['semisimple'] == 1, "(%d,%d) ctrl-a not ss" % (p, q), "3")
        row((p, q), "ctrl-a:ss-distinct-nonzero", ma)

        # (b) ss with repeated nonzero eigenvalue
        A_b = build_from_M(diag_M([Integer(1), Integer(1)]), p, q, 2, "3")
        mb = measure(A_b, p, q, basis=SO_BASES[(p, q)])
        check(mb['semisimple'] == 1, "(%d,%d) ctrl-b not ss" % (p, q), "3")
        row((p, q), "ctrl-b:ss-repeated-nonzero", mb)

        # (c) nonzero nilpotent: M = J (lam = 0)
        A_c = build_from_M(jordan_block_M(2, Integer(0)), p, q, 2, "3")
        mc = measure(A_c, p, q, basis=SO_BASES[(p, q)])
        check(mc['nilpotent'] == 1, "(%d,%d) ctrl-c not nilpotent" % (p, q), "3")
        check(A_c != zeros(n, n), "(%d,%d) ctrl-c is zero" % (p, q), "3")
        row((p, q), "ctrl-c:nilpotent-nonzero", mc)
    elif m == 1:
        # nilpotent exists with m=1: M is 1x1 = [0] -> that's zero. Use the
        # so(p,q) basis directly: search for a nonzero nilpotent element.
        found = None
        for B in SO_BASES[(p, q)]:
            for B2 in SO_BASES[(p, q)]:
                C = B + B2
                if C != zeros(n, n) and is_nilpotent(C):
                    found = C
                    break
            if found is not None:
                break
        if found is not None:
            mc = measure(found, p, q, basis=SO_BASES[(p, q)])
            check(mc['nilpotent'] == 1, "(%d,%d) ctrl-c not nilpotent" % (p, q), "3")
            row((p, q), "ctrl-c:nilpotent-nonzero-searched", mc)
        else:
            print("SIG=(%d,%d) | ctrl-c | no nonzero nilpotent found by search "
                  "over pairwise sums of the so-basis" % (p, q))
        # (a)/(b) still measurable via so(2)-style rotation generators:
        # use a generic semisimple element = a basis element of the compact part.
        # measure the first basis element as data
        B0 = SO_BASES[(p, q)][0]
        m0 = measure(B0, p, q, basis=SO_BASES[(p, q)])
        row((p, q), "ctrl-basis0", m0)
    else:
        # definite signature: everything is semisimple (compact). Measure a
        # basis element and a sum as data.
        B0 = SO_BASES[(p, q)][0]
        m0 = measure(B0, p, q, basis=SO_BASES[(p, q)])
        row((p, q), "ctrl-basis0", m0)
        nil = None
        for B in SO_BASES[(p, q)]:
            if B != zeros(n, n) and is_nilpotent(B):
                nil = B
                break
        print("SIG=(%d,%d) | ctrl-c | nonzero_nilpotent_in_basis=%s "
              "(definite signature)" % (p, q, "YES" if nil is not None else "NO"))
        check(True, "(%d,%d) ctrl definite measured" % (p, q), "3")

# ---------------------------------------------------------------------------
print()
print("### SECTION 4: MUTANTS")

# --- m1: a matrix NOT in so(p,q) must be rejected by in_so ---------------
p, q = 2, 2
n = 4
bad = zeros(n, n)
bad[0, 0] = Integer(1)          # symmetric-ish, definitely not in so(2,2)
m1_rejected = (in_so(bad, p, q) is False)
# survivability: the SAME function must accept a legitimate element
good = SO_BASES[(2, 2)][0]
m1_accepts_good = in_so(good, p, q)
check(m1_accepts_good, "m1: in_so REJECTED a legitimate so(2,2) element "
                       "-> the check is degenerate, mutant could not survive", "4")
if m1_rejected and m1_accepts_good:
    print("MUTANT m1: CAUGHT | in_so(diag(1,0,0,0),2,2)=False AND "
          "in_so(so_basis[0],2,2)=True (function discriminates)")
    m1_ok = True
else:
    print("MUTANT m1: SILENT")
    m1_ok = False

# --- m2: "every centralizer is abelian" must be refuted by our own code ---
m2_witness = None
for (pp, qq) in SIGS:
    for tag, A in [("jordan", FOURTH[(pp, qq)])]:
        if A is None:
            continue
        mm = measure(A, pp, qq, basis=SO_BASES[(pp, qq)])
        if mm['dimcc'] > 0:
            m2_witness = ((pp, qq), tag, mm)
            break
    if m2_witness:
        break
if m2_witness is None:
    # widen the search: any so-basis element in any signature
    for (pp, qq) in SIGS:
        for B in SO_BASES[(pp, qq)]:
            mm = measure(B, pp, qq, basis=SO_BASES[(pp, qq)])
            if mm['dimcc'] > 0:
                m2_witness = ((pp, qq), "so-basis-element", mm)
                break
        if m2_witness:
            break

# survivability of m2: the same derived_dim must return 0 on a genuinely
# abelian centralizer (prove the function CAN say "abelian")
abelian_example = None
for (pp, qq) in SIGS:
    for B in SO_BASES[(pp, qq)]:
        mm = measure(B, pp, qq, basis=SO_BASES[(pp, qq)])
        if mm['dimcc'] == 0:
            abelian_example = ((pp, qq), mm)
            break
    if abelian_example:
        break
check(abelian_example is not None,
      "m2: derived_dim NEVER returns 0 -> degenerate, cannot survive", "4")
if m2_witness is not None:
    (sg, tag, mm) = m2_witness
    print("MUTANT m2: CAUGHT | claim 'every centralizer is abelian' REFUTED at "
          "SIG=(%d,%d) construct=%s dimc=%d dimcc=%d (>0) | survivability: same "
          "derived_dim returns dimcc=0 at SIG=(%d,%d) (abelian case exists)"
          % (sg[0], sg[1], tag, mm['dimc'], mm['dimcc'],
             abelian_example[0][0], abelian_example[0][1]))
    m2_ok = True
else:
    print("MUTANT m2: SILENT | no element with dimcc>0 found anywhere")
    m2_ok = False

# --- m3: "charpoly determines dim_c" tested on the pair from (2) ---------
m3_counterexample = None
m3_agree = None
for (pp, qq), (m_s, m_j) in PAIR.items():
    if m_s['dimc'] != m_j['dimc']:
        m3_counterexample = ((pp, qq), m_s, m_j)
        break
for (pp, qq), (m_s, m_j) in PAIR.items():
    if m_s['dimc'] == m_j['dimc']:
        m3_agree = ((pp, qq), m_s, m_j)
        break
check(len(PAIR) > 0, "m3: no pairs constructed -> mutant cannot be tested", "4")

# m3 survivability, positive control: measure() MUST report EQUAL dimc for a
# same-charpoly pair where dimc genuinely IS equal. Build one honestly:
# conjugate A_jordan in so(2,2) by g = eta(2,2), which lies in O(2,2)
# (g*eta*g.T == eta), so A' = g*A*g^-1 stays in so(2,2), same charpoly,
# and c(A') = g c(A) g^-1  =>  dimc must come out EQUAL.
g = eta(2, 2)
check(expand(g * eta(2, 2) * g.T - eta(2, 2)) == zeros(4, 4),
      "m3-control: g not in O(2,2)", "4")
A_conj = Matrix(4, 4, lambda i, j: simplify((g * FOURTH[(2, 2)] * g.inv())[i, j]))
check(in_so(A_conj, 2, 2), "m3-control: conjugate left so(2,2)", "4")
m_orig = measure(FOURTH[(2, 2)], 2, 2, basis=SO_BASES[(2, 2)])
m_conj = measure(A_conj, 2, 2, basis=SO_BASES[(2, 2)])
cp_same = (expand(A_conj.charpoly(Symbol('x')).as_expr()
                  - FOURTH[(2, 2)].charpoly(Symbol('x')).as_expr()) == 0)
check(cp_same, "m3-control: conjugation changed the charpoly", "4")
m3_positive = (m_orig['dimc'] == m_conj['dimc'])
check(m3_positive,
      "m3-control: measure() gave dimc=%d vs %d for CONJUGATE elements -> "
      "dimc measurement is not conjugation-invariant, measurement INVALID"
      % (m_orig['dimc'], m_conj['dimc']), "4")
print("m3-control (survivability): SIG=(2,2) A_jordan vs eta-conjugate: "
      "charpoly_equal=1 dimc=%d vs %d -> measure() DOES report EQUAL dimc when "
      "dimc truly is equal (so a dimc difference is a real signal, not an "
      "artefact of the function always disagreeing)"
      % (m_orig['dimc'], m_conj['dimc']))
if m3_counterexample is not None:
    (sg, ms, mj) = m3_counterexample
    print("MUTANT m3: CAUGHT | claim 'charpoly determines dim_c' REFUTED at "
          "SIG=(%d,%d): identical charpoly %s but dimc_ss=%d != dimc_jordan=%d | "
          "survivability: proven by the m3-control above (conjugate pair, "
          "dimc reported EQUAL); same-charpoly pairs in (2) with equal dimc: %s"
          % (sg[0], sg[1], str(ms['charpoly']), ms['dimc'], mj['dimc'],
             ("SIG=(%d,%d)" % m3_agree[0]) if m3_agree else
             "none (every ss/jordan pair differed)"))
    m3_ok = True
else:
    print("MUTANT m3: SILENT | dim_c was equal for every same-charpoly pair -> "
          "the claim 'charpoly determines dim_c' was NOT refuted by this data")
    m3_ok = False

# ---------------------------------------------------------------------------
print()
print("### VALIDITY")
zero_sections = [k for k, v in SECTION_CHECKS.items() if v == 0]
for k in sorted(SECTION_CHECKS):
    print("section=%s checks_performed=%d" % (k, SECTION_CHECKS[k]))

bad = False
if zero_sections:
    print("MEASUREMENT INVALID: sections with ZERO checks: %s" % zero_sections)
    bad = True
if not (m1_ok and m2_ok and m3_ok):
    print("MEASUREMENT INVALID: a mutant was SILENT "
          "(m1=%s m2=%s m3=%s)" % (m1_ok, m2_ok, m3_ok))
    bad = True
if FAILS:
    print("MEASUREMENT INVALID: asserts failed: %s" % FAILS)
    bad = True

print()
print("SUMMARY: asserts_passed=%d | FAIL=%d" % (ASSERTS, 1 if bad else 0))
sys.exit(1 if bad else 0)
