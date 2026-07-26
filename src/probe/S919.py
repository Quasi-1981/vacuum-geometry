# -*- coding: utf-8 -*-
# DIM: na (dimensionless/exact symbolic count; no spatial lattice)
"""
S919_w29_nilpotent_strata.py — nilpotent strata of so(p,q).

Signature (p,q), eta=diag(+1 x p, -1 x q), so(p,q)={X: X*eta + eta*X^T = 0}.
"Wedge" construction: A = x*(eta*y)^T - y*(eta*x)^T.

Sections:
  N0 — non-existence of a nonzero nilpotent in the definite signatures so(n,0)/so(0,n), n=3,4,5;
  N1 — construction and taxonomy (soft/deep) in 10 signatures + machine proofs of non-existence;
  N2 — orbits: exp(tA)=1 <=> t=0 (the N-orbit does not close);
  N3 — centralizers: exact nullspace, recheck, closure, abelian-ness, dim[c,c],
       honest analysis of the pattern with counterexamples;
  N4 — isometric consistency (conjugation by a rational isometry);
  M  — mutants m1/m2/m3 (each MUST be caught);
  SUMMARY — verdicts + the full N3 table.

EXACT arithmetic: sympy Integer/Rational. NO float/random/time/datetime.
Run: python S919_w29_nilpotent_strata.py > S919_w29_run.log 2>&1
Exit 0 only if all checks OK, all mutants caught, no section empty.
"""

import sys
import sympy as sp
from sympy import Matrix, Rational, Integer, symbols, eye, zeros

sys.stdout.reconfigure(encoding="utf-8")
sys.stderr.reconfigure(encoding="utf-8")

# ---------------------------------------------------------------- bookkeeping
CHECKS = 0
FAILS = []
SECTION_CHECKS = {}
VERDICTS = []
CUR = "?"
MUTANTS_CAUGHT = 0
MUTANTS_TOTAL = 3


def begin(name, title):
    global CUR
    CUR = name
    SECTION_CHECKS[name] = 0
    print()
    print("=" * 78)
    print(f"SECTION {name} — {title}")
    print("=" * 78)


def check(label, cond):
    global CHECKS
    CHECKS += 1
    SECTION_CHECKS[CUR] = SECTION_CHECKS.get(CUR, 0) + 1
    if cond:
        print(f"  [OK]   {label}")
    else:
        print(f"  [FAIL] {label}")
        FAILS.append(f"{CUR}: {label}")
    return bool(cond)


def end_section(verdict):
    n = SECTION_CHECKS.get(CUR, 0)
    if n == 0:
        FAILS.append(f"{CUR}: ZERO checks — the count did not happen")
        print(f"  [FAIL] {CUR}: ZERO checks — the count did not happen")
    print(f"  -- checks in section {CUR}: {n}")
    print(f"VERDICT {CUR}: {verdict}")
    VERDICTS.append((CUR, verdict))


# ---------------------------------------------------------------- linear algebra (exact)
def eta_mat(p, q):
    return sp.diag(*([Integer(1)] * p + [Integer(-1)] * q))


def ip(et, u, v):
    """<u,v> = u^T eta v"""
    return sp.expand((u.T * et * v)[0, 0])


def wedge(et, x, y):
    """A = x (eta y)^T - y (eta x)^T"""
    return (x * (et * y).T - y * (et * x).T).expand()


def is_zero_mat(M):
    return M.expand() == zeros(M.rows, M.cols)


def in_so(et, X):
    return is_zero_mat(X * et + et * X.T)


def e_(n, i):
    """i — 1-based unit vector"""
    v = zeros(n, 1)
    v[i - 1, 0] = Integer(1)
    return v


def vecm(M):
    return Matrix(M.rows * M.cols, 1, list(M))


def so_basis_nullspace(et):
    """Basis of so(p,q) = exact nullspace of the defining system X*eta+eta*X^T=0."""
    n = et.rows
    xs = list(symbols(f"u0:{n * n}"))
    X = Matrix(n, n, xs)
    E = (X * et + et * X.T).expand()
    eqs = [E[i, j] for i in range(n) for j in range(n)]
    M, _ = sp.linear_eq_to_matrix(eqs, xs)
    return [Matrix(n, n, list(v)) for v in M.nullspace()]


def centralizer(et, A, basis):
    """c(A) = {X in so: XA-AX=0} — exact nullspace; returns (basis of c, system matrix)."""
    m = len(basis)
    cs = list(symbols(f"c0:{m}"))
    X = zeros(A.rows, A.cols)
    for c, B in zip(cs, basis):
        X = X + c * B
    Cm = (X * A - A * X).expand()
    eqs = [Cm[i, j] for i in range(A.rows) for j in range(A.cols)]
    M, _ = sp.linear_eq_to_matrix(eqs, cs)
    cent = [coeffs_to_mat(v, basis, A.rows) for v in M.nullspace()]
    return cent, M


def coeffs_to_mat(v, basis, n):
    Z = zeros(n, n)
    for coeff, B in zip(list(v), basis):
        Z = Z + coeff * B
    return Z.expand()


def bra(X, Y):
    return (X * Y - Y * X).expand()


def in_span(mats, Z):
    if not mats:
        return is_zero_mat(Z)
    B = Matrix.hstack(*[vecm(m) for m in mats])
    return Matrix.hstack(B, vecm(Z)).rank() == B.rank()


def span_dim(mats):
    nz = [m for m in mats if not is_zero_mat(m)]
    if not nz:
        return 0
    return Matrix.hstack(*[vecm(m) for m in nz]).rank()


def charpoly_is_lambda_n(A):
    lam = symbols("lam_cp")
    return sp.expand(A.charpoly(lam).as_expr() - lam ** A.rows) == 0


def nilp_exp(A, t):
    """exp(tA) for nilpotent A — a finite polynomial, exact Rational coefficients."""
    n = A.rows
    S = eye(n)
    Ak = eye(n)
    fact = Integer(1)
    k = 0
    while True:
        Ak = (Ak * A).expand()
        k += 1
        if is_zero_mat(Ak):
            break
        if k > n:
            raise RuntimeError("A is not nilpotent — exp is not a polynomial")
        fact = fact * k
        S = S + (t ** k) * Ak / fact
    return S.expand()


def prove_no_isotropic_partner(et, x):
    """
    Machine proof: all y with {<y,y>=0, <x,y>=0} are collinear with x (for a fixed null x).
    Method: linear equation -> pivot; residual quadratic form q2 in the free variables;
    check +-PSD (exact, rational); q2(v)=0 with PSD => v in the kernel => lift the kernel
    into y-space, check collinearity with x. Returns (ok, description).
    """
    n = et.rows
    ys = list(symbols(f"y0:{n}", real=True))
    y = Matrix(n, 1, ys)
    lin = ip(et, x, y)
    pivot = None
    for v in ys:
        if lin.coeff(v) != 0:
            pivot = v
            break
    if pivot is None:
        return False, "the linear equation <x,y>=0 is trivial"
    sols = sp.solve(sp.Eq(lin, Integer(0)), pivot)
    if len(sols) != 1:
        return False, "the pivot did not solve uniquely"
    sub = {pivot: sols[0]}
    q2 = sp.expand(ip(et, y, y).subs(sub))
    free = [v for v in ys if v is not pivot]
    H = sp.hessian(q2, free) / Integer(2)
    if H.is_positive_semidefinite:
        Mp, sgn = H, "+"
    elif (-H).is_positive_semidefinite:
        Mp, sgn = -H, "-"
    else:
        return False, "the residual form is not sign-definite — the proof path does not work"
    ker = Mp.nullspace()
    if len(ker) < 1:
        return False, "the kernel is empty (not even y=x was found — a contradiction)"
    for kv in ker:
        smap = {fv: kv[i] for i, fv in enumerate(free)}
        yfull = Matrix(n, 1, [sp.expand(sp.sympify(s).subs(sub).subs(smap)) for s in ys])
        if is_zero_mat(yfull):
            continue
        if Matrix.hstack(x, yfull).rank() > 1:
            return False, f"a non-collinear solution y={list(yfull)}"
    return True, (f"form {sgn}PSD, dim ker={len(ker)}, "
                  f"all solutions (linear lift of the kernel) lie in span(x)")


# ================================================================ SECTION N0
begin("N0", "NON-EXISTENCE of a nonzero nilpotent in the definite signatures")

for n in (3, 4, 5):
    # general antisymmetric X with symbolic super-diagonal entries
    syms = {}
    X = zeros(n, n)
    varlist = []
    for i in range(n):
        for j in range(i + 1, n):
            a = symbols(f"a_{i}{j}", real=True)
            X[i, j] = a
            X[j, i] = -a
            varlist.append(a)
    sum_sq = sum(X[i, j] ** 2 for i in range(n) for j in range(i + 1, n))
    ident = sp.expand(sp.trace(X * X) + 2 * sum_sq)
    check(f"n={n}: identity tr(X^2) + 2*Sum_(i<j) X[i,j]^2 == 0 (symbolic)", ident == 0)
    # membership: an antisymmetric X satisfies the so-definition for both eta=+I and eta=-I
    check(f"n={n}: X*(+I)+(+I)*X^T=0 for a general antisymmetric X",
          in_so(eye(n), X))
    check(f"n={n}: X*(-I)+(-I)*X^T=0 for a general antisymmetric X",
          in_so(-eye(n), X))
    # defining system: dim so(n,0) = n(n-1)/2, all basis elements antisymmetric
    bpos = so_basis_nullspace(eye(n))
    check(f"n={n}: dim so(n,0) (nullspace of the defining system) = {len(bpos)} = n(n-1)/2 = {n*(n-1)//2}",
          len(bpos) == n * (n - 1) // 2)
    check(f"n={n}: every basis element of so(n,0) is antisymmetric",
          all(is_zero_mat(B + B.T) for B in bpos))
    bneg = so_basis_nullspace(-eye(n))
    check(f"n={n}: dim so(0,n) = {len(bneg)} = n(n-1)/2 and all antisymmetric",
          len(bneg) == n * (n - 1) // 2 and all(is_zero_mat(B + B.T) for B in bneg))
    # -tr(X^2) = 2*Sum X_ij^2 — a positive DEFINITE form in the variables X_ij (exact check)
    Hn = sp.hessian(-sp.trace(X * X), varlist) / Integer(2)
    check(f"n={n}: the form -tr(X^2) is positive DEFINITE (Hessian/2 PD, exact) "
          f"=> tr(X^2)=0 <=> X=0", bool(Hn.is_positive_definite))

print("""
  Derivation chain (machine-checked links marked):
    X nilpotent (X^k=0) => X^2 nilpotent => tr(X^2)=0                 [trace of a nilpotent = 0, standard]
    tr(X^2) = -2*Sum_(i<j) X[i,j]^2                                    [MACHINE: identity]
    -tr(X^2) — a positive definite form                                [MACHINE: PD check]
    => tr(X^2)=0 only at X=0 (over R)                                 => nilpotent != 0 does not exist.
  For so(0,n) (eta=-I) the algebra is THE SAME (antisymmetric) — the same conclusion.""")

end_section("nilpotent X != 0 does NOT exist in the definite signatures so(n,0), so(0,n) "
            "(n=3,4,5 machine-checked; the argument goes exclusively through tr(X^2) as a "
            "+-definite form).")

# ================================================================ SECTION N1
begin("N1", "\"WEDGE\" CONSTRUCTION AND TAXONOMY (soft / deep)")

SIGS = [(2, 1), (1, 2), (3, 1), (1, 3), (2, 2), (3, 2), (2, 3), (4, 1), (1, 4), (3, 3)]
NILP = []          # {'p','q','n','typ','A','x','y','et','rankA','rankA2','dimker'}
TAX = {}           # sig -> {'soft': str, 'deep': str}

for (p, q) in SIGS:
    n = p + q
    et = eta_mat(p, q)
    print(f"\n--- signature ({p},{q}), n={n}, eta=diag({','.join(str(et[i,i]) for i in range(n))}) ---")

    # the wedge lies in so(p,q) for FULLY symbolic x,y
    gx = Matrix(n, 1, symbols(f"gx0:{n}"))
    gy = Matrix(n, 1, symbols(f"gy0:{n}"))
    check(f"({p},{q}): wedge A(x,y) in so(p,q) for symbolic x,y (identity)",
          in_so(et, wedge(et, gx, gy)))

    # ---- SOFT type: x null, y orthogonal, y*y != 0
    x = e_(n, 1) + e_(n, p + 1)                       # first (+) axis + first (-) axis
    y = e_(n, 2) if p >= 2 else e_(n, p + 2)          # another axis, orthogonal to x
    xx, xy, yy = ip(et, x, x), ip(et, x, y), ip(et, y, y)
    check(f"({p},{q}) soft: x={list(x)} null (<x,x>={xx})", xx == 0)
    check(f"({p},{q}) soft: y={list(y)} orthogonal to x (<x,y>={xy})", xy == 0)
    check(f"({p},{q}) soft: <y,y>={yy} != 0", yy != 0)
    check(f"({p},{q}) soft: x,y linearly independent", Matrix.hstack(x, y).rank() == 2)
    A = wedge(et, x, y)
    A2 = (A * A).expand()
    A3 = (A2 * A).expand()
    check(f"({p},{q}) soft: A in so(p,q)", in_so(et, A))
    check(f"({p},{q}) soft: A^2 != 0", not is_zero_mat(A2))
    check(f"({p},{q}) soft: A^3 == 0", is_zero_mat(A3))
    rA, rA2 = A.rank(), A2.rank()
    dk = n - rA
    check(f"({p},{q}) soft: rank A={rA} (expect 2), rank A^2={rA2} (expect 1), dim ker A={dk}",
          rA == 2 and rA2 == 1)
    check(f"({p},{q}) soft: dim ker A via nullspace = {len(A.nullspace())} = n - rank",
          len(A.nullspace()) == dk)
    check(f"({p},{q}) soft: charpoly(A) = lam^{n}", charpoly_is_lambda_n(A))
    check(f"({p},{q}) soft: tr(A)=0 and tr(A^2)=0 (consistent with the N0 chain)",
          sp.trace(A) == 0 and sp.trace(A2) == 0)
    NILP.append(dict(p=p, q=q, n=n, typ="soft", A=A, x=x, y=y, et=et,
                     rankA=rA, rankA2=rA2, dimker=dk))
    TAX.setdefault((p, q), {})["soft"] = "Yes"

    # ---- DEEP type: x,y both null, <x,y>=0, independent (a totally isotropic plane)
    if p >= 2 and q >= 2:
        x2 = e_(n, 1) + e_(n, p + 1)
        y2 = e_(n, 2) + e_(n, p + 2)
        c1 = ip(et, x2, x2) == 0
        c2 = ip(et, y2, y2) == 0
        c3 = ip(et, x2, y2) == 0
        c4 = Matrix.hstack(x2, y2).rank() == 2
        check(f"({p},{q}) deep: x={list(x2)}, y={list(y2)} — totally isotropic plane "
              f"(<x,x>=<y,y>=<x,y>=0, independent)", c1 and c2 and c3 and c4)
        Ad = wedge(et, x2, y2)
        Ad2 = (Ad * Ad).expand()
        check(f"({p},{q}) deep: A in so(p,q)", in_so(et, Ad))
        check(f"({p},{q}) deep: A^2 == 0", is_zero_mat(Ad2))
        check(f"({p},{q}) deep: A != 0", not is_zero_mat(Ad))
        rAd = Ad.rank()
        dkd = n - rAd
        check(f"({p},{q}) deep: rank A={rAd} (expect 2), rank A^2=0, dim ker A={dkd}", rAd == 2)
        check(f"({p},{q}) deep: charpoly(A) = lam^{n}", charpoly_is_lambda_n(Ad))
        NILP.append(dict(p=p, q=q, n=n, typ="deep", A=Ad, x=x2, y=y2, et=et,
                         rankA=rAd, rankA2=0, dimker=dkd))
        TAX[(p, q)]["deep"] = "Yes"
    else:
        # min(p,q)=1: proof of non-existence of a totally isotropic plane
        ok1, msg1 = prove_no_isotropic_partner(et, x)
        check(f"({p},{q}) deep: NON-EXISTENCE for x={list(x)}: {msg1}", ok1)
        x_alt = e_(n, 1) - e_(n, p + 1)
        ok2, msg2 = prove_no_isotropic_partner(et, x_alt)
        check(f"({p},{q}) deep: NON-EXISTENCE for alt. x={list(x_alt)}: {msg2}", ok2)
        print(f"         (coverage: proved for 2 representatives of a null x; an arbitrary null x "
              f"reduces to a representative by transitivity of O(p,q) on null vectors — Witt, 📖)")
        TAX[(p, q)]["deep"] = "NONE (proved)"

print("\n  TAXONOMY:")
print(f"  {'signature':<10} {'soft':<8} {'deep'}")
for (p, q) in SIGS:
    print(f"  ({p},{q}){'':<5} {TAX[(p,q)]['soft']:<8} {TAX[(p,q)]['deep']}")

end_section("the soft type (A^2!=0, A^3=0) exists in ALL 10 signatures; the deep type (A^2=0, a totally "
            "isotropic plane) exists EXACTLY where min(p,q)>=2: (2,2),(3,2),(2,3),(3,3); "
            "for min(p,q)=1 non-existence is proved by machine (PSD reduction: every isotropic "
            "partner is collinear with x).")

# ================================================================ SECTION N2
begin("N2", "ORBITS: exp(tA) = 1 <=> t = 0")

t = symbols("t", real=True)
for item in NILP:
    p, q, n, typ, A, et = item["p"], item["q"], item["n"], item["typ"], item["A"], item["et"]
    E = nilp_exp(A, t)
    check(f"({p},{q}) {typ}: exp(tA) — an isometry: exp^T*eta*exp - eta == 0 (symbolic in t)",
          is_zero_mat((E.T * et * E - et).expand()))
    check(f"({p},{q}) {typ}: exp(0*A) = 1", (E.subs(t, Integer(0))).expand() == eye(n))
    D = (E - eye(n)).expand()
    lin_entry = None
    for i in range(n):
        for j in range(n):
            if D[i, j] == 0:
                continue
            poly = sp.Poly(D[i, j], t)
            if poly.degree() == 1:
                lin_entry = (i, j, poly.LC())
                break
        if lin_entry:
            break
    ok = lin_entry is not None and lin_entry[2] != 0
    if ok:
        i, j, c = lin_entry
        check(f"({p},{q}) {typ}: entry [{i},{j}] of exp(tA)-1 = ({c})*t (purely linear, coeff != 0) "
              f"=> exp(tA)=1 requires t=0", True)
    else:
        check(f"({p},{q}) {typ}: a purely linear entry of exp(tA)-1 was found", False)

end_section("for EVERY constructed nilpotent (both types, all signatures) exp(tA)=1 <=> t=0 "
            "(a concrete matrix entry is linear in t with a nonzero coefficient); "
            "the one-parameter N-orbit {exp(tA)} ~ R — never closes.")

# ================================================================ SECTION N3
begin("N3", "CENTRALIZERS (the main section)")

SO_CACHE = {}
ROWS = []

for item in NILP:
    p, q, n, typ, A, et = item["p"], item["q"], item["n"], item["typ"], item["A"], item["et"]
    key = (p, q)
    if key not in SO_CACHE:
        b = so_basis_nullspace(et)
        check(f"so({p},{q}): dim (nullspace of the defining system) = {len(b)} = n(n-1)/2 = {n*(n-1)//2}",
              len(b) == n * (n - 1) // 2)
        check(f"so({p},{q}): all basis elements satisfy X*eta+eta*X^T=0",
              all(in_so(et, B) for B in b))
        SO_CACHE[key] = b
    basis = SO_CACHE[key]

    cent, _Meq = centralizer(et, A, basis)
    dimc = len(cent)
    check(f"({p},{q}) {typ}: c(A) — exact nullspace, dim c = {dimc} >= 1", dimc >= 1)
    check(f"({p},{q}) {typ}: EACH of the {dimc} elements of c(A) rechecked by commutation [Z,A]=0",
          all(is_zero_mat(bra(Z, A)) for Z in cent))
    check(f"({p},{q}) {typ}: every element of c(A) lies in so(p,q)",
          all(in_so(et, Z) for Z in cent))
    check(f"({p},{q}) {typ}: A in c(A)", in_span(cent, A))
    pairs = [(i, j) for i in range(dimc) for j in range(i + 1, dimc)]
    brackets = [bra(cent[i], cent[j]) for (i, j) in pairs]
    check(f"({p},{q}) {typ}: closure under the bracket — all {len(brackets)} pairwise [Zi,Zj] in span(c(A))",
          all(in_span(cent, B) for B in brackets))
    abel = all(is_zero_mat(B) for B in brackets)
    dimcc = span_dim(brackets)
    if abel:
        check(f"({p},{q}) {typ}: ABELIAN (all {len(brackets)} brackets = 0), dim[c,c]=0", dimcc == 0)
    else:
        check(f"({p},{q}) {typ}: NOT abelian, dim[c,c] = {dimcc}", dimcc >= 1)
    ROWS.append(dict(p=p, q=q, n=n, typ=typ, rankA=item["rankA"], rankA2=item["rankA2"],
                     dimker=item["dimker"], dimc=dimc, abel=abel, dimcc=(0 if abel else dimcc)))


def print_table(rows):
    print(f"  {'(p,q)':<7} {'type':<10} {'rkA':<4} {'rkA2':<5} {'dim ker A':<10} "
          f"{'dim c':<6} {'abelian':<9} {'dim[c,c]':<9} mark")
    for r in rows:
        mark = "<== (3,1)/(1,3)" if (r["p"], r["q"]) in ((3, 1), (1, 3)) else ""
        print(f"  ({r['p']},{r['q']}){'':<2} {r['typ']:<10} {r['rankA']:<4} {r['rankA2']:<5} "
              f"{r['dimker']:<10} {r['dimc']:<6} {('yes' if r['abel'] else 'NO'):<9} "
              f"{r['dimcc']:<9} {mark}")


print("\n  TABLE N3:")
print_table(ROWS)
check("table N3 is nonempty and covers every constructed nilpotent",
      len(ROWS) == len(NILP) and len(ROWS) >= 14)

# ---- honest analysis of the pattern
print("\n  PATTERN ANALYSIS (candidates against the data; counterexamples listed):")


def nonab(r):
    return not r["abel"]


BICONDS = [
    ("K1: non-abelian <=> type=deep", lambda r: r["typ"] == "deep"),
    ("K2: non-abelian <=> min(p,q)>=2", lambda r: min(r["p"], r["q"]) >= 2),
    ("K3: non-abelian <=> dim ker A >= 3", lambda r: r["dimker"] >= 3),
    ("K4: non-abelian <=> algebraic multiplicity of zero >= 5 (multiplicity = n always, charpoly=lam^n)",
     lambda r: r["n"] >= 5),
    ("K5: non-abelian <=> dim c >= 4", lambda r: r["dimc"] >= 4),
    ("K6: non-abelian <=> (deep OR n>=5)",
     lambda r: r["typ"] == "deep" or r["n"] >= 5),
]
IMPLS = [
    ("I1: deep => non-abelian", lambda r: (r["typ"] != "deep") or nonab(r)),
    ("I2: soft => abelian", lambda r: (r["typ"] != "soft") or r["abel"]),
]

survivors = []
for name, pred in BICONDS:
    mism = [r for r in ROWS if pred(r) != nonab(r)]
    if mism:
        exs = "; ".join(f"({r['p']},{r['q']}) {r['typ']}: predicate={pred(r)}, non-abelian={nonab(r)}"
                        for r in mism)
        print(f"    {name}: REJECTED. Counterexamples: {exs}")
    else:
        print(f"    {name}: CONSISTENT with all {len(ROWS)} rows")
        survivors.append(name)
for name, pred in IMPLS:
    mism = [r for r in ROWS if not pred(r)]
    if mism:
        exs = "; ".join(f"({r['p']},{r['q']}) {r['typ']} (abelian={r['abel']})" for r in mism)
        print(f"    {name}: REJECTED. Counterexamples: {exs}")
    else:
        print(f"    {name}: CONSISTENT with all {len(ROWS)} rows")
        survivors.append(name)

same_n_diff = [(a, b) for a in ROWS for b in ROWS
               if a["n"] == b["n"] and a["dimker"] == b["dimker"] and a["abel"] != b["abel"]]
if same_n_diff:
    a, b = same_n_diff[0]
    print(f"    Note: the algebraic multiplicity of zero = n in ALL rows (charpoly=lam^n) and cannot "
          f"distinguish rows with the same n; the pair ({a['p']},{a['q']}) {a['typ']} (abelian={a['abel']}) "
          f"vs ({b['p']},{b['q']}) {b['typ']} (abelian={b['abel']}) has the same n={a['n']} and "
          f"dim ker A={a['dimker']} — the multiplicity and dim ker do NOT determine abelian-ness.")
check("pattern analysis performed (candidates evaluated, survivors listed)",
      len(BICONDS) + len(IMPLS) == 8)

verdict_n3 = ("pattern from the data: " +
              ("; ".join(survivors) if survivors else "NO candidate survived") +
              ". Neither the type itself, nor dim ker A, nor the algebraic multiplicity (=n), nor the "
              "signature alone determine abelian-ness — counterexamples listed above.")
end_section(verdict_n3)

# ================================================================ SECTION N4
begin("N4", "ISOMETRIC CONSISTENCY (conjugation by a rational isometry)")


def rot_block(n, i, j):
    """rotation 3/5-4/5 on axes i,j (0-based), same-sign pair"""
    S = eye(n)
    S[i, i] = Rational(3, 5); S[i, j] = Rational(4, 5)
    S[j, i] = Rational(-4, 5); S[j, j] = Rational(3, 5)
    return S


def boost_block(n, i, j):
    """boost 5/4-3/4 on axes i,j (0-based), mixed-sign pair"""
    S = eye(n)
    S[i, i] = Rational(5, 4); S[i, j] = Rational(3, 4)
    S[j, i] = Rational(3, 4); S[j, j] = Rational(5, 4)
    return S


def find_nilp(p, q, typ):
    return next(it for it in NILP if it["p"] == p and it["q"] == q and it["typ"] == typ)


def row_of(p, q, typ):
    return next(r for r in ROWS if r["p"] == p and r["q"] == q and r["typ"] == typ)


N4_CASES = [
    # ((p,q), type, S constructor, description)
    ((3, 1), "soft", lambda: rot_block(4, 0, 1) * boost_block(4, 2, 3),
     "R(axes 1,2; both +) * B(axes 3,4; signs +,-)"),
    ((2, 2), "deep", lambda: rot_block(4, 0, 1) * boost_block(4, 1, 2),
     "R(axes 1,2; both +) * B(axes 2,3; signs +,-)"),
]

for (p, q), typ, mkS, descr in N4_CASES:
    item = find_nilp(p, q, typ)
    et, A, n = item["et"], item["A"], item["n"]
    S = mkS()
    check(f"({p},{q}) {typ}: S = {descr} — an isometry: S^T*eta*S == eta (exact)",
          is_zero_mat((S.T * et * S - et).expand()))
    Ac = (S * A * S.inv()).expand()
    check(f"({p},{q}) {typ}: A' = S*A*S^(-1) in so(p,q)", in_so(et, Ac))
    check(f"({p},{q}) {typ}: A' != A (the conjugation is nontrivial)", not is_zero_mat((Ac - A).expand()))
    cent2, _ = centralizer(et, A, SO_CACHE[(p, q)])
    centc, _ = centralizer(et, Ac, SO_CACHE[(p, q)])
    r = row_of(p, q, typ)
    check(f"({p},{q}) {typ}: dim c(A') = {len(centc)} = dim c(A) = {r['dimc']} (invariant)",
          len(centc) == r["dimc"])
    pairs_c = [(i, j) for i in range(len(centc)) for j in range(i + 1, len(centc))]
    brs = [bra(centc[i], centc[j]) for (i, j) in pairs_c]
    abel_c = all(is_zero_mat(B) for B in brs)
    check(f"({p},{q}) {typ}: abelian-ness of c(A') = {abel_c} = abelian-ness of c(A) = {r['abel']} (invariant)",
          abel_c == r["abel"])
    dimcc_c = span_dim(brs) if not abel_c else 0
    check(f"({p},{q}) {typ}: dim[c(A'),c(A')] = {dimcc_c} = dim[c,c] = {r['dimcc']} (invariant)",
          dimcc_c == r["dimcc"])

end_section("dim c, abelian-ness and dim[c,c] are invariant under a rational isometric conjugation "
            "(checked on the (3,1)-soft and (2,2)-deep cases; S^T*eta*S=eta exact).")

# ================================================================ SECTION M
begin("M", "MUTANTS (each MUST be caught)")

# ---- m1: broken centralizer solver (drops equations until the solution changes)
item22 = find_nilp(2, 2, "deep")
A22, et22 = item22["A"], item22["et"]
basis22 = SO_CACHE[(2, 2)]
cent_h, Meq = centralizer(et22, A22, basis22)
r0 = Meq.rank()
Mb = Meq.copy()
removed = 0
while Mb.rows > 0 and Mb.rank() >= r0:
    Mb.row_del(0)
    removed += 1
broken = [coeffs_to_mat(v, basis22, 4) for v in Mb.nullspace()]
print(f"  m1: (2,2)-deep; the true rank of the system = {r0}; dropped {removed} equations; "
      f"broken dim 'c' = {len(broken)} vs true {len(cent_h)}")
check("m1: the broken solver GAVE a different (larger) space", len(broken) > len(cent_h))
bad1 = [Z for Z in broken if not is_zero_mat(bra(Z, A22))]
caught1 = len(bad1) >= 1
print(f"  m1: recheck by commutation [Z,A]=0 on each of the {len(broken)} elements: "
      f"{len(bad1)} FALSE ones found -> {'CAUGHT' if caught1 else 'NOT-CAUGHT'}")
if caught1:
    MUTANTS_CAUGHT += 1
check("m1 CAUGHT: false elements of the broken centralizer caught by commutation", caught1)

# ---- m2: broken abelian-ness test (always 'yes') on the so(3) reference
def honest_abelian(mats):
    return all(is_zero_mat(bra(mats[i], mats[j]))
               for i in range(len(mats)) for j in range(i + 1, len(mats)))

def broken_abelian(mats):
    return True  # MUTANT: always 'yes'

so3 = SO_CACHE.get((3, 0)) or so_basis_nullspace(eye(3))
h_ab = honest_abelian(so3)
b_ab = broken_abelian(so3)
wit = None
for i in range(len(so3)):
    for j in range(i + 1, len(so3)):
        B = bra(so3[i], so3[j])
        if not is_zero_mat(B):
            wit = (i, j, B)
            break
    if wit:
        break
print(f"  m2: so(3) reference: honest test = {h_ab} (witness: [B{wit[0]},B{wit[1]}] != 0), "
      f"broken test = {b_ab}")
check("m2: the so(3) reference is honestly NOT abelian (witness bracket != 0)", h_ab is False and wit is not None)
caught2 = (b_ab != h_ab)
print(f"  m2: discrepancy of the broken test vs the honest one -> {'CAUGHT' if caught2 else 'NOT-CAUGHT'}")
if caught2:
    MUTANTS_CAUGHT += 1
check("m2 CAUGHT: the broken abelian-ness test disagrees with the honest one on so(3)", caught2)

# ---- m3: the false claim "the soft type has A^2=0"
item31 = find_nilp(3, 1, "soft")
A31 = item31["A"]
A31sq = (A31 * A31).expand()
wit3 = None
for i in range(4):
    for j in range(4):
        if A31sq[i, j] != 0:
            wit3 = (i, j, A31sq[i, j])
            break
    if wit3:
        break
caught3 = wit3 is not None
if caught3:
    print(f"  m3: the claim \"soft => A^2=0\" REFUTED exactly: (3,1)-soft, "
          f"A^2[{wit3[0]},{wit3[1]}] = {wit3[2]} != 0 -> CAUGHT")
    MUTANTS_CAUGHT += 1
else:
    print("  m3: no nonzero entry of A^2 was found -> NOT-CAUGHT")
check("m3 CAUGHT: the false claim refuted by a concrete nonzero entry of A^2", caught3)

check(f"all mutants caught: {MUTANTS_CAUGHT}/{MUTANTS_TOTAL}", MUTANTS_CAUGHT == MUTANTS_TOTAL)
end_section(f"mutants caught {MUTANTS_CAUGHT}/{MUTANTS_TOTAL} "
            "(m1 — commutation recheck; m2 — discrepancy with the honest bracket on so(3); "
            "m3 — a concrete nonzero entry of A^2).")

# ================================================================ SECTION SUMMARY
begin("SUMMARY", "VERDICTS OF ALL SECTIONS + THE FULL N3 TABLE")

for name, v in VERDICTS:
    print(f"  VERDICT {name}: {v}")

print("\n  FULL TABLE N3 (repeated):")
print_table(ROWS)

check("SUMMARY: all previous sections had > 0 checks",
      all(SECTION_CHECKS.get(s, 0) > 0 for s in ("N0", "N1", "N2", "N3", "N4", "M")))
check(f"SUMMARY: mutants {MUTANTS_CAUGHT}/{MUTANTS_TOTAL} caught", MUTANTS_CAUGHT == MUTANTS_TOTAL)

per_sec = ", ".join(f"{k}={v}" for k, v in SECTION_CHECKS.items())
print(f"\n  COUNTERS: total checks = {CHECKS}; by section: {per_sec}; FAIL = {len(FAILS)}")
if FAILS:
    print("  FAIL LIST:")
    for f in FAILS:
        print(f"    - {f}")

end_section(f"checks={CHECKS}, FAIL={len(FAILS)}, mutants={MUTANTS_CAUGHT}/{MUTANTS_TOTAL}; "
            + ("ALL OK" if not FAILS else "THERE ARE FAILURES"))

# ---------------------------------------------------------------- exit
if CHECKS == 0:
    print("EXIT=1 (zero checks — the count did not happen)")
    sys.exit(1)
if FAILS or MUTANTS_CAUGHT != MUTANTS_TOTAL:
    print("EXIT=1")
    sys.exit(1)
print("EXIT=0")
sys.exit(0)
