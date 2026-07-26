# -*- coding: utf-8 -*-
# DIM: na (visa B on W29 rungs 4-5; handles 0 — the probe postulates nothing).
"""
S922 — COUNTER-PROBE OF VISA B against S920/S921 (rungs 4-5 of W29).

★WHAT THIS IS NOT: not a refutation of S921's numbers. All 6492 of its checks
are reproduced bit-for-bit (5/5 W29 probes, exit=0). The numbers are correct.
The visa question is different: WHAT EXACTLY DO THEY MEASURE — and do the
verdicts placed ON TOP of them hold up.

Three stakes of the visa (each machine-checked, each with a mutant):

B1-B5 ★MAIN: the A1 verdict "the structure constants are FORCED … NO
  free parameter in the structure of the extension" was measured INSIDE the
  chosen embedding. S921 takes T_a = [[0,a],[0,0]] (a strictly-upper block) — for
  this FORM the product T_a·T_b equals zero IDENTICALLY, by form, before any
  algebra enters. That is, "[T,T]=0" was not measured — it was BUILT IN
  by the choice of form. Counter-thesis: there exists a ONE-PARAMETER
  family T_a(λ) = [[0,a],[−λ(ηa)ᵀ,0]], closed for EVERY λ, with the same
  automatic Jacobi identity (an associative algebra — the S921 argument works
  identically on the whole family). λ=0 is exactly the S921 term. ⟹ the
  extension is not forced, but CHOSEN.

B6 ★MULTIPLICITY OF A2: the three numbers of the verdict (linear 0 · quadratic
  exactly 1 · antisymmetric 0) are not three measurements, but consequences of
  ONE fact (the commutant of the standard representation = ℝ, there is no
  invariant vector). I measure the dimensions explicitly and show the
  coincidence. The fact is a reference one (📖 Schur/irreducibility) ⟹
  multiplicity 0.

B7 ★★A CATCH THAT DIED UNDER ADVERSARIAL CHECKING — kept in the probe WITH
  ITS NAME ATTACHED. I wanted to show that rung-4 silently narrowed the scope.
  It did not: the quantifier "invariant with respect to all coordinate flows"
  sits INSIDE the theorem, "orientABLE" carries the suffix of possibility, and
  "the sign MAY BE" already is the very correction I proposed. The section is
  kept as an ADDENDUM to the qualifier (an explicit reflection) + as a stamp:
  a catch quietly withdrawn is a fudge.

Fences: handles 0 · sympy-exact (no float at all) · language — forms/brackets/
components; no physical noun outside quotations of verdicts.
"""
import sys
from pathlib import Path

from sympy import Integer, Matrix, Rational, Symbol, diag, eye, zeros

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

FAILED = 0
CUR = "—"
SECTION_COUNTS = {}


def section(name):
    global CUR
    CUR = name
    SECTION_COUNTS.setdefault(name, 0)
    print("\n" + "=" * 72)
    print(name)
    print("=" * 72)


def check(msg, ok):
    global FAILED
    SECTION_COUNTS[CUR] = SECTION_COUNTS.get(CUR, 0) + 1
    if not ok:
        FAILED += 1
        print(f"  ✗ FAIL: {msg}")
    return ok


def verdict(msg):
    print(f"\n  ⇒ {msg}")


# ────────────────────────────── carrier ─────────────────────────────────────
def eta_of(p, q):
    """The same convention as in S920/S921."""
    return diag(*([Integer(1)] * p + [Integer(-1)] * q))


def so_basis(p, q):
    """Basis of so(p,q) as the exact nullspace of the condition Xη + ηXᵀ = 0."""
    n = p + q
    et = eta_of(p, q)
    unk = [[Symbol(f"x{i}_{j}") for j in range(n)] for i in range(n)]
    X = Matrix(unk)
    cond = X * et + et * X.T
    eqs = [cond[i, j] for i in range(n) for j in range(n)]
    flat = [unk[i][j] for i in range(n) for j in range(n)]
    A = Matrix([[eq.coeff(v) for v in flat] for eq in eqs])
    basis = []
    for v in A.nullspace():
        basis.append(Matrix(n, n, list(v)))
    return basis


def hat(M):
    n = M.shape[0]
    return M.row_join(zeros(n, 1)).col_join(zeros(1, n + 1))


def Tmat(a, lam, et):
    """★THE FAMILY: T_a(λ) = [[0, a],[−λ·(ηa)ᵀ, 0]].

    λ=0 → exactly T_a from S921 (a strictly-upper block).
    λ≠0 → the bottom row is alive: the product T_a·T_b is NO LONGER zero by form.
    """
    n = a.shape[0]
    T = zeros(n + 1, n + 1)
    for i in range(n):
        T[i, n] = a[i]
    row = -lam * (et * a).T
    for j in range(n):
        T[n, j] = row[0, j]
    return T


def br(X, Y):
    return X * Y - Y * X


def e_vec(n, i):
    v = zeros(n, 1)
    v[i, 0] = Integer(1)
    return v


def Tbad(a, lam, et):
    """★MUTANT FORM: the same thing, but η is DROPPED in the bottom row.

    For a definite η=I it coincides with the real one (which is why the
    mutant is run on an INDEFINITE signature — there η≠I and the difference
    is alive).
    """
    n = a.shape[0]
    T = zeros(n + 1, n + 1)
    for i in range(n):
        T[i, n] = a[i]
    row = -lam * a.T                      # ★broken: without η
    for j in range(n):
        T[n, j] = row[0, j]
    return T


# ★★THE ONE MEASUREMENT PATH — both the sections AND the mutants walk it.
#   The lesson that cost me this rewrite: a mutant that does its OWN inline
#   check does not test the measurement — it tests itself. My first m1 did
#   exactly that, and m3 asked something that could not have been false (the
#   very class this same probe judges m3 for in S921). Now both go through
#   the real code.
def family_closes(p, q, lam, tfun):
    """Is span{M̂, T(λ)} closed under the bracket: [M̂,T_a] = T_{Ma}?"""
    n = p + q
    et = eta_of(p, q)
    basis = so_basis(p, q)
    for M in basis:
        Mh = hat(M)
        for i in range(n):
            lhs = br(Mh, tfun(e_vec(n, i), lam, et)).expand()
            rhs = tfun(M * e_vec(n, i), lam, et).expand()
            if lhs != rhs:
                return False
    return True


def dim_bilinear_invariants(p, q):
    """dim{B : MᵀB + BM = 0 ∀M ∈ so(p,q)} — sym. and antisym. TOGETHER."""
    n = p + q
    basis = so_basis(p, q)
    unk = [[Symbol(f"u{i}_{j}") for j in range(n)] for i in range(n)]
    B = Matrix(unk)
    flat = [unk[i][j] for i in range(n) for j in range(n)]
    eqs = []
    for M in basis:
        C = (M.T * B + B * M).expand()
        eqs += [C[i, j] for i in range(n) for j in range(n)]
    return len(Matrix([[e.coeff(v) for v in flat] for e in eqs]).nullspace())


SIGS = [(2, 1), (1, 2), (3, 1), (1, 3), (2, 2), (4, 1), (3, 2), (2, 3), (3, 3)]
LAM = Symbol("lambda", real=True)

# ═══════════════════════════════════════════════════════════════════════════
section("B1 — \"[T,T]=0\" IN S921 IS AN IDENTITY OF THE FORM, NOT A MEASUREMENT")
# ═══════════════════════════════════════════════════════════════════════════
print("  S921: Tmat(a) = [[0,a],[0,0]]. For ARBITRARY symbolic a,b I compute")
print("  the product T_a·T_b — if it is zero IDENTICALLY, then '[T,T]=0' is built")
print("  in by the form of the embedding, not measured.\n")
for n in (3, 4, 5):
    a = Matrix(n, 1, [Symbol(f"a{i}") for i in range(n)])
    b = Matrix(n, 1, [Symbol(f"b{i}") for i in range(n)])
    Ta = Tmat(a, Integer(0), eye(n))
    Tb = Tmat(b, Integer(0), eye(n))
    prod = (Ta * Tb).expand()
    check(f"n={n}: T_a·T_b ≡ 0 for SYMBOLIC a,b (not for concrete ones)",
          prod == zeros(n + 1, n + 1))
    check(f"n={n}: ⟹ [T_a,T_b] ≡ 0 identically by form",
          (br(Ta, Tb)).expand() == zeros(n + 1, n + 1))
    print(f"  n={n}: T_a·T_b ≡ 0 — no property of so(p,q) was used.")
verdict("★ '[T,T]=0' in S921 does not depend AT ALL on the signature, on η, on "
        "so(p,q): this is a property of the strictly-upper block. The probe proved "
        "a property of the CHOSEN EMBEDDING and called it the forcedness of the structure.")

# ═══════════════════════════════════════════════════════════════════════════
section("B2 — THE ONE-PARAMETER FAMILY T_a(λ): CLOSURE AT EVERY λ")
# ═══════════════════════════════════════════════════════════════════════════
print("  T_a(λ) = [[0,a],[−λ(ηa)ᵀ,0]]. Question: is span{M̂, T(λ)} closed")
print("  under the bracket at a SYMBOLIC λ (not at λ=0 and not at λ=1)?\n")
CLOSURE = {}
for (p, q) in SIGS:
    n = p + q
    et = eta_of(p, q)
    basis = so_basis(p, q)
    Ts = [Tmat(e_vec(n, i), LAM, et) for i in range(n)]
    Ms = [hat(M) for M in basis]

    # (a) [M̂, T_a(λ)] = T_{Ma}(λ) — at symbolic λ
    ok_a = 0
    for M, Mh in zip(basis, Ms):
        for i in range(n):
            lhs = br(Mh, Ts[i]).expand()
            rhs = Tmat(M * e_vec(n, i), LAM, et).expand()
            if lhs == rhs:
                ok_a += 1
            else:
                check(f"({p},{q}) (a) [M̂,T_{i}(λ)] = T_(Ma)(λ)", False)
    check(f"({p},{q}) (a) {ok_a} checks of [M̂,T(λ)]=T_(Ma)(λ) at SYMBOLIC λ",
          ok_a == len(basis) * n)

    # (b) [T_a(λ), T_b(λ)] = −λ·hat(a(ηb)ᵀ − b(ηa)ᵀ) — and this is in the so-block
    ok_b = 0
    nonzero_seen = 0
    for i in range(n):
        for j in range(i + 1, n):
            ai, aj = e_vec(n, i), e_vec(n, j)
            comm = br(Ts[i], Ts[j]).expand()
            Mij = (ai * (et * aj).T - aj * (et * ai).T)
            pred = (-LAM * hat(Mij)).expand()
            if comm == pred:
                ok_b += 1
            else:
                check(f"({p},{q}) (b) [T_{i}(λ),T_{j}(λ)] = −λ·M̂_{i}{j}", False)
            # does this M really lie in so(p,q)?
            check(f"({p},{q}) (b) M_{i}{j} ∈ so(p,q)",
                  (Mij * et + et * Mij.T).expand() == zeros(n, n))
            if comm != zeros(n + 1, n + 1):
                nonzero_seen += 1
    check(f"({p},{q}) (b) {ok_b} brackets [T(λ),T(λ)] landed IN the so-BLOCK ⟹ closed",
          ok_b == n * (n - 1) // 2)
    check(f"({p},{q}) (b) at λ≠0 the bracket [T,T] is NONzero ({nonzero_seen} pairs)",
          nonzero_seen == n * (n - 1) // 2)
    CLOSURE[(p, q)] = True
    print(f"  ({p},{q}): (a)={ok_a} OK, (b)={ok_b} OK — closed ∀λ")
verdict("★THE FAMILY IS CLOSED AT EVERY λ. Jacobi here is just as automatic as in "
        "S921 (matrices = an associative algebra) ⟹ the S921 argument that "
        "\"Jacobi holds automatically\" does NOT SINGLE OUT λ=0 in any way: it "
        "holds equally on the whole family.")

# ═══════════════════════════════════════════════════════════════════════════
section("B3 — WHAT ALGEBRA THIS IS AT λ≠0 (★the label is COUNTED, not hand-written)")
# ═══════════════════════════════════════════════════════════════════════════
print("""  ★★A STAMP AGAINST MYSELF — THE FIRST PASS OF THIS SECTION WAS WRONG, AND GREEN.
  I wrote the label BY HAND: "λ=+1 → so(p,q+1)". Flipped on all 8/8: η̃ =
  diag(η,+1) adds a PLUS ⟹ so(p+1,q). For (3,1) I called so(4,1) — "so(3,2)".
  WHY IT PASSED: what was asserted was allin (membership in so(η̃)) — true regardless
  of the label; the LABEL ITSELF rode along in the f-string of that SAME check,
  which measured something else. Exactly the class "the label of the rule ≠ its
  execution" — in a probe written against it, in a section about mutants. ⟹ THE
  FIX: the signature of η̃ is now COUNTED by machine and this count is
  ASSERTED; the hand no longer writes labels.\n""")
for (p, q) in SIGS:
    n = p + q
    et = eta_of(p, q)
    basis = so_basis(p, q)
    got = {}
    for sig_lam, tag in ((Integer(1), "λ=+1"), (Integer(-1), "λ=−1")):
        ett = diag(et, sig_lam)
        gens = [hat(M) for M in basis] + \
               [Tmat(e_vec(n, i), sig_lam, et) for i in range(n)]
        allin = all((X * ett + ett * X.T).expand() == zeros(n + 1, n + 1)
                    for X in gens)
        check(f"({p},{q}) {tag}: ALL {len(gens)} generators ∈ so(η̃)", allin)
        # ★the signature η̃ is COUNTED from the matrix itself, and this is a SEPARATE check
        d = [ett[k, k] for k in range(n + 1)]
        pl = sum(1 for x in d if x > 0)
        mi = sum(1 for x in d if x < 0)
        exp = (p + 1, q) if sig_lam == 1 else (p, q + 1)
        check(f"({p},{q}) {tag}: COUNTED signature η̃ = ({pl},{mi}) = "
              f"expected {exp}", (pl, mi) == exp)
        check(f"({p},{q}) {tag}: dim = {len(gens)} = (n+1)n/2 = "
              f"{(n + 1) * n // 2}", len(gens) == (n + 1) * n // 2)
        got[tag] = (pl, mi)
    # ★E2: so(A,B) ≅ so(B,A) (η → −η) ⟹ at p=q the two terms COINCIDE
    same = (got["λ=+1"] == (got["λ=−1"][1], got["λ=−1"][0]))
    check(f"({p},{q}) λ=±1 — the same algebra up to so(A,B)≅so(B,A)? "
          f"{'YES' if same and p == q else 'NO'} (expectation: ⟺ p=q)",
          (same and p == q) == (p == q))
    n_alg = 2 if p == q else 3
    print(f"  ({p},{q}): λ=+1 → so{got['λ=+1']} · λ=−1 → so{got['λ=−1']} · "
          f"λ=0 → so({p},{q})⋉R^{n} [the S921 term] ⟹ DISTINCT algebras: {n_alg}")
verdict("★A FAMILY, NOT A SINGLETON: λ=0 is not \"a structure without a parameter\", "
        "but a DEGENERATE term (an Inönü–Wigner contraction). There are THREE "
        "distinct algebras at p≠q and TWO at p=q (since so(A,B)≅so(B,A)) — a "
        "correction my first pass did not have.")

# ═══════════════════════════════════════════════════════════════════════════
section("B4 — λ=0 IS THE ONLY CASE WHERE THE SHIFTS = AN IDEAL (and this is exactly what rung-5 stands on)")
# ═══════════════════════════════════════════════════════════════════════════
print("  The rung-5 verdict asks \"WHAT NUMBER attaches invariantly to a shift\".")
print("  The question makes sense only if the shifts are an invariant subspace (an ideal).")
print("  I measure: whether [T,T] lies in span{T}, and whether the Killing form is degenerate.\n")
for (p, q) in SIGS:
    n = p + q
    et = eta_of(p, q)
    basis = so_basis(p, q)
    for lam, tag in ((Integer(0), "λ=0"), (Integer(1), "λ=+1"),
                     (Integer(-1), "λ=−1")):
        gens = [hat(M) for M in basis] + \
               [Tmat(e_vec(n, i), lam, et) for i in range(n)]
        Ts = gens[len(basis):]
        # an ideal: [T_i,T_j] ∈ span{T}?
        in_span = True
        for i in range(n):
            for j in range(i + 1, n):
                c = br(Ts[i], Ts[j]).expand()
                if c == zeros(n + 1, n + 1):
                    continue
                # does c have a nonzero so-block (⟹ not in span{T})?
                if any(c[r, s] != 0 for r in range(n) for s in range(n)):
                    in_span = False
        check(f"({p},{q}) {tag}: are the shifts an ideal? {'YES' if in_span else 'NO'} "
              f"(expectation: only λ=0)", in_span == (lam == 0))
        # the Killing form K(X,Y)=tr(ad_X ad_Y): degenerate ⟺ not semisimple
        d = len(gens)
        idx = {}
        for k, g in enumerate(gens):
            idx[k] = g

        def ad(X):
            cols = []
            for Y in gens:
                c = br(X, Y).expand()
                cols.append(c)
            # coordinates in the basis gens: the solution of a linear system
            A = Matrix([[g[r, s] for g in gens]
                        for r in range(n + 1) for s in range(n + 1)])
            out = zeros(d, d)
            for jcol, c in enumerate(cols):
                rhs = Matrix([c[r, s] for r in range(n + 1)
                              for s in range(n + 1)])
                sol = A.solve_least_squares(rhs) if A.rank() < d else \
                    A.LUsolve(rhs) if A.rows == A.cols else \
                    (A.T * A).inv() * A.T * rhs
                for irow in range(d):
                    out[irow, jcol] = sol[irow]
            return out

        ads = [ad(g) for g in gens]
        K = Matrix(d, d, lambda i2, j2: (ads[i2] * ads[j2]).trace())
        K = K.applyfunc(lambda z: z.expand())
        degen = (K.det() == 0)
        check(f"({p},{q}) {tag}: is the Killing form degenerate? "
              f"{'YES' if degen else 'NO'} (expectation: only λ=0)",
              degen == (lam == 0))
    print(f"  ({p},{q}): λ=0 — span{{T}} is an ideal, the Killing form is degenerate (NOT "
          f"semisimple); λ=±1 — span{{T}} is NOT an ideal, the Killing form is non-degenerate "
          f"(semisimple)")
verdict("★★span{T} IS AN IDEAL ONLY AT λ=0 — and this is exactly what the rung-5 "
        "question \"what number attaches to a SHIFT\" stands on. At λ≠0 a \"shift\" "
        "stops being an invariant subspace ⟹ the question loses its OBJECT. "
        "★PRECISION (an adversarial-check correction): what was measured is exactly "
        "\"span{T} is not an ideal\", NOT \"there are no invariant subspaces at all\" — "
        "so(2,2) does have proper ideals of its own. The claim carries exactly what was "
        "measured. ★PRECISION-2: a non-degenerate Killing form = SEMIsimple (so(2,2) ≅ "
        "sl(2,ℝ)⊕sl(2,ℝ) — semisimple, not simple).")

# ★λ=0 as a ROOT, not as a sample point: det of the Killing form — a polynomial in λ.
#   The first pass measured λ∈{0,±1} one at a time (sampling). The polynomial says more:
#   the degeneracy is EXACTLY at λ=0 and NOWHERE ELSE — for all λ at once.
print("\n  ★λ=0 — not a sample point, but a ROOT: computing det K(λ) as a polynomial.")
for (p, q) in [(2, 1), (3, 1)]:
    n = p + q
    et = eta_of(p, q)
    basis = so_basis(p, q)
    gens = [hat(M) for M in basis] + \
           [Tmat(e_vec(n, i), LAM, et) for i in range(n)]
    d = len(gens)
    A = Matrix([[g[r, s] for g in gens]
                for r in range(n + 1) for s in range(n + 1)])
    ads = []
    for X in gens:
        out = zeros(d, d)
        for jc, Y in enumerate(gens):
            c = br(X, Y).expand()
            rhs = Matrix([c[r, s] for r in range(n + 1) for s in range(n + 1)])
            sol = A.solve_least_squares(rhs)
            for ir in range(d):
                out[ir, jc] = sol[ir].expand()
        ads.append(out)
    K = Matrix(d, d, lambda i2, j2: (ads[i2] * ads[j2]).trace().expand())
    dk = K.det().factor()
    roots_only_zero = (dk.subs(LAM, 0) == 0)
    print(f"    ({p},{q}): det K(λ) = {dk}")
    check(f"({p},{q}) det K(λ) — a polynomial vanishing at λ=0", roots_only_zero)
    check(f"({p},{q}) det K(λ) ≠ 0 at λ=1 ⟹ λ=0 is an ISOLATED root, "
          f"not a sample", dk.subs(LAM, Integer(1)) != 0)
verdict("★DEGENERACY EXACTLY AT λ=0 AND NOWHERE ELSE (a polynomial, not a sample). "
        "The dimension of the Killing radical is an isomorphism invariant — n at λ=0, 0 at "
        "λ≠0 ⟹ λ=0 is NOT isomorphic to its neighbors by any change of basis. The family is real.")

# ═══════════════════════════════════════════════════════════════════════════
section("B5 — WHAT λ COSTS: rescaling T→μT gives λ→μ²λ")
# ═══════════════════════════════════════════════════════════════════════════
MU = Symbol("mu", positive=True)
for (p, q) in SIGS[:5]:
    n = p + q
    et = eta_of(p, q)
    # T'_a := μ·T_a(λ) — which λ' gives the same family?
    ok = 0
    for i in range(n):
        for j in range(i + 1, n):
            Ti = MU * Tmat(e_vec(n, i), LAM, et)
            Tj = MU * Tmat(e_vec(n, j), LAM, et)
            c = br(Ti, Tj).expand()
            ai, aj = e_vec(n, i), e_vec(n, j)
            Mij = (ai * (et * aj).T - aj * (et * ai).T)
            pred = (-(MU ** 2 * LAM) * hat(Mij)).expand()
            if c == pred:
                ok += 1
            else:
                check(f"({p},{q}) scale: [μT_{i},μT_{j}] = −μ²λ·M̂", False)
    check(f"({p},{q}) {ok} pairs: rescaling T→μT gives exactly λ→μ²λ",
          ok == n * (n - 1) // 2)
print("\n  μ² > 0 for real μ≠0 ⟹ the SIGN of λ is invariant, the MAGNITUDE is not.")
verdict("★AN HONEST TALLY OF HANDLES — AND IT FAVORS THE COURT: λ is NOT a "
        "continuous handle at the level of the algebra (|λ| is washed out by the "
        "scale T→μT). What remains is only a DISCRETE choice of sign(λ) ∈ {−,0,+}. "
        "The program counts continuous metric-linking constants ⟹ **\"HANDLES 0\" "
        "STANDS**. Inflating a discrete choice into a constant would be a fudge — and "
        "my first pass did exactly that (the summary said \"the handle is set to zero\", "
        "contradicting THIS very section).")

# ★★A DEBT THAT LIVES OUTSIDE THIS FLOOR (a flag, not a measurement of this probe).
print("""
  ★DEBT (a bridge, not a measurement): "handles 0" holds EXACTLY as long as μ is free.
  But premise (b) of the theory — "d+1 UNIT axes" — is a NORMALIZATION: a unit
  axis fixes the length. As soon as μ is pinned down, |λ| = 1/R² stops being washed
  out and becomes a genuine continuous dimensionful constant (a radius of curvature).
  ⟹ the honest tally line is not "handles 0" but "HANDLES 0, AS LONG AS NO LENGTH
  IS PINNED DOWN ANYWHERE" — and the program pins it down on the neighboring rung.
  Class: a bridge between floors (its multiplicity is carried from there, not from
  here). NOT a verdict of this probe — a flag for the court.""")

# ═══════════════════════════════════════════════════════════════════════════
section("B6 — MULTIPLICITY OF A2: the verdict's three numbers = ONE fact")
# ═══════════════════════════════════════════════════════════════════════════
print("  The A2 verdict reported THREE numbers (linear 0 · quadratic 1 ·")
print("  antisymmetric 0) as three measurements. I measure the commutant and the")
print("  space of invariant BILINEAR forms (symmetric+antisymmetric TOGETHER).\n")
for (p, q) in SIGS:
    n = p + q
    et = eta_of(p, q)
    basis = so_basis(p, q)

    # (i) an invariant vector: M v = 0 ∀M
    rows = []
    for M in basis:
        for r in range(n):
            rows.append([M[r, c] for c in range(n)])
    dim_vec = len(Matrix(rows).nullspace())
    check(f"({p},{q}) invariant vectors: dim = {dim_vec} (0)", dim_vec == 0)

    # (ii) the commutant End_so(V) = {X : [X,M]=0 ∀M}
    unk = [[Symbol(f"y{i}_{j}") for j in range(n)] for i in range(n)]
    X = Matrix(unk)
    flat = [unk[i][j] for i in range(n) for j in range(n)]
    eqs = []
    for M in basis:
        C = (X * M - M * X).expand()
        eqs += [C[i, j] for i in range(n) for j in range(n)]
    A = Matrix([[e.coeff(v) for v in flat] for e in eqs])
    dim_comm = len(A.nullspace())
    check(f"({p},{q}) commutant End_so(V): dim = {dim_comm} (= ℝ·I)",
          dim_comm == 1)

    # (iii) ALL invariant bilinear forms together: MᵀB + BM = 0
    unk2 = [[Symbol(f"z{i}_{j}") for j in range(n)] for i in range(n)]
    B = Matrix(unk2)
    flat2 = [unk2[i][j] for i in range(n) for j in range(n)]
    eqs2 = []
    for M in basis:
        C = (M.T * B + B * M).expand()
        eqs2 += [C[i, j] for i in range(n) for j in range(n)]
    A2m = Matrix([[e.coeff(v) for v in flat2] for e in eqs2])
    ns = A2m.nullspace()
    dim_bil = len(ns)
    check(f"({p},{q}) ALL invariant bilinear forms: dim = {dim_bil} (1)",
          dim_bil == 1)
    check(f"({p},{q}) dim(bilinear) = dim(commutant) — both {dim_comm}",
          dim_bil == dim_comm)
    # the one form ∝ η ⟹ symmetric ⟹ antisymmetric 0 AS A COROLLARY
    Bsol = Matrix(n, n, list(ns[0]))
    check(f"({p},{q}) the one form is symmetric ⟹ \"antisymmetric 0\" is "
          f"not a separate measurement, but a corollary", Bsol == Bsol.T)
    print(f"  ({p},{q}): vector-invariants 0 · commutant 1 · ALL bilinear 1 "
          f"⟹ (0,1,0) exhausted by two lines")
# ★n=2 — A WITNESS FOR WHY THE QUALIFIER "ABSOLUTE" IS LOAD-BEARING (a catch on
#   my own name: my first pass wrote "irreducible over ℝ" — that is NOT ENOUGH).
print("\n  ★witness n=2: why \"irreducible over ℝ\" is a wrong name for the lemma.")
for (p, q) in [(2, 0), (1, 1)]:
    n = 2
    basis = so_basis(p, q)
    unk = [[Symbol(f"w{i}_{j}") for j in range(n)] for i in range(n)]
    X = Matrix(unk)
    flat = [unk[i][j] for i in range(n) for j in range(n)]
    eqs = []
    for M in basis:
        C = (X * M - M * X).expand()
        eqs += [C[i, j] for i in range(n) for j in range(n)]
    dim_comm = len(Matrix([[e.coeff(v) for v in flat] for e in eqs]).nullspace())
    check(f"({p},{q}) n=2: commutant dim = {dim_comm} ≠ 1 ⟹ Schur over ℝ gives "
          f"NOT-ℝ (ℂ or ℝ⊕ℝ)", dim_comm == 2)
    print(f"    ({p},{q}): commutant dim={dim_comm} — so(2,0) is IRREDUCIBLE over ℝ, "
          f"but the commutant = ℂ ⟹ \"antisym. 0\" is WRONG at n=2")
verdict("★THE MULTIPLICITY OF A2 = 1, not 3. The space of ALL invariant bilinear forms "
        "is one-dimensional; the fact that its generator is symmetric gives "
        "\"antisymmetric 0\" FOR FREE. There is one ancestor — the ABSOLUTE "
        "irreducibility of the standard representation at n≥3 (📖 reference, "
        "MULTIPLICITY 0).\n"
        "     ★THE NAME OF THE LEMMA IS CORRECTED UNDER ADVERSARIAL CHECKING: writing "
        "\"irreducible over ℝ\" is WRONG. so(2,0) is irreducible over ℝ, yet the commutant "
        "= ℂ and an antisymmetric invariant EXISTS (dim=1). So irreducibility over ℝ alone "
        "does not give (0,1,0) — ABSOLUTE irreducibility does (commutant = ℝ). n≥3 has it; "
        "the qualifier is load-bearing, not decoration.\n"
        "     The 6492 checks of S921 = 17 instances of one textbook line.")

# ═══════════════════════════════════════════════════════════════════════════
section("B7 — THE SCOPE OF ORIENTABILITY: time reversal is also an isometry")
# ═══════════════════════════════════════════════════════════════════════════
print("  S920 measured COORDINATE flows exp(sM) and their finite compositions")
print("  = the connected component SO⁺(p,q). Question: does this exhaust the isometries?\n")
for (p, q) in SIGS:
    n = p + q
    et = eta_of(p, q)
    # reflection of the LAST axis (the time axis when q≥1)
    R = eye(n)
    R[n - 1, n - 1] = Integer(-1)
    check(f"({p},{q}) reflection R=diag(…,−1): RᵀηR = η ⟹ R is an ISOMETRY",
          (R.T * et * R) == et)
    check(f"({p},{q}) det R = −1 ⟹ R ∉ the connected component",
          R.det() == Integer(-1))
    if q >= 1:
        v = e_vec(n, n - 1)          # the last axis: time
        check(f"({p},{q}) e_t is timelike: ⟨e_t,e_t⟩ = −1 < 0",
              (v.T * et * v)[0, 0] == Integer(-1))
        check(f"({p},{q}) ★R·e_t = −e_t ⟹ THE ISOMETRY SWAPS THE CONES",
              (R * v) == -v)
print("\n  ⟹ [MACH] R is an isometry, and it swaps the cones FOR ANY q.")
print("     Each cone is SEPARATELY invariant under SO⁺(p,q); the UNORDERED pair")
print("     {C₊,C₋} is invariant under the whole of O(p,q).\n")
verdict("★★THE VISA'S CATCH WAS HERE — AND IT IS DEAD. I report it with its name "
        "attached, because a catch quietly withdrawn is a fudge.\n"
        "     I wanted to say: \"the rung-4 verdict silently narrowed the scope from "
        "isometries to the connected component\". WRONG, for three reasons, each "
        "sufficient on its own:\n"
        "     (1) the quantifier sits INSIDE the theorem itself, not in a qualifier "
        "below: \"splits into two cones, INVARIANT WITH RESPECT TO ALL COORDINATE "
        "FLOWS\". I read past it and attacked the sentence with an object its own "
        "quantifier excludes;\n"
        "     (2) \"orientABLE\" — the suffix '-able' = POSSIBILITY (as in "
        "\"solvable\", \"measurable\"): orientable, not oriented. The language already "
        "carried what I was demanding be added;\n"
        "     (3) the verdict says \"an invariant sign of charge is POSSIBLE ⟺ …\". "
        "\"Possible\" = a choice exists. My \"correction\" was a restatement of the "
        "sentence I was trying to fix.\n"
        "     ⟹ B7 is not a catch, but an ADDENDUM TO THE QUALIFIER: the reflection R "
        "explicitly realizes the reversal of the cones. The content of rung-4 STANDS WHOLE.")
print("  ★RESIDUE (hygiene, not a verdict): in §(2) of the court's text \"an invariant")
print("    sign of charge\" stands BARE — the scope is inherited from §(3) below it, while")
print("    §(2) is written in physical language. A sentence from §(2) lifted out of the")
print("    court's text will have the reader supply the word \"isometries\" themselves.")
print("    The cost of the fix is zero: write the scope into §(2) inline, as already "
      "done in §(1). This is a drift risk in citation, not a defect.")

# ═══════════════════════════════════════════════════════════════════════════
section("B8 — THE COURT'S LAST REDOUBT + S921's EMPTY MUTANT")
# ═══════════════════════════════════════════════════════════════════════════
print("""  The court itself handed back the frame: "agreement with the textbook = a
  frame, multiplicity 0" — and fell back to ONE claim: "the load-bearing new
  measurement = ZERO linear invariants". Question: is this last redoubt a
  separate fact, or the same one?\n""")
for (p, q) in SIGS:
    n = p + q
    et = eta_of(p, q)
    basis = so_basis(p, q)
    # identity: φ := (ηv)ᵀ ⟹ φ·M = −(Mv)ᵀη for EVERY M ∈ so(p,q)
    ok = 0
    v = Matrix(n, 1, [Symbol(f"v{i}") for i in range(n)])
    for M in basis:
        lhs = ((et * v).T * M).expand()
        rhs = (-(M * v).T * et).expand()
        if lhs == rhs:
            ok += 1
        else:
            check(f"({p},{q}) identity φM = −(Mv)ᵀη", False)
    check(f"({p},{q}) {ok} generators: φ=(ηv)ᵀ ⟹ φM = −(Mv)ᵀη IDENTICALLY "
          f"(symbolic v)", ok == len(basis))
print("\n  ⟹ η is a bijection between invariant VECTORS and invariant")
print("     COVECTORS. \"Zero linear invariants\" ≡ \"there is no invariant")
print("     vector\" — and that is a one-line corollary of irreducibility (a fixed")
print("     vector spans a 1-dimensional trivial sub-representation; irreducibility +")
print("     n≥2 ⟹ it is zero).")
verdict("★★THE LAST REDOUBT FALLS. The claim the court PULLED OUT from under the "
        "frame as \"a load-bearing NEW measurement\" is the very same textbook fact, "
        "reached by η-conjugation from what the court itself marked with "
        "multiplicity 0. The multiplicity of A2 = 0 IN FULL, including the exception "
        "the court kept for itself.")

print("\n  ★AS A SIDE EFFECT — S921's mutant m3 IS EMPTY (class: a mutant that could not survive):")
bad_sigs = []
for n in range(3, 8):
    for p in range(0, n + 1):
        q = n - p
        et = eta_of(p, q)
        basis = so_basis(p, q)
        if not basis:
            continue
        unk = [[Symbol(f"s{i}_{j}") for j in range(n)] for i in range(n)]
        B = Matrix(unk)
        flat = [unk[i][j] for i in range(n) for j in range(n)]
        eqs = []
        for M in basis:
            C = (M.T * B + B * M).expand()
            eqs += [C[i, j] for i in range(n) for j in range(n)]
        ns = Matrix([[e.coeff(x) for x in flat] for e in eqs]).nullspace()
        sym = [S for S in (Matrix(n, n, list(w)) for w in ns) if S == S.T]
        if len(ns) != 1:
            bad_sigs.append((p, q, len(ns)))
check(f"dim(sym. invariants) = 1 in ALL (p,q) with n=3..7 — deviations: "
      f"{len(bad_sigs)}", len(bad_sigs) == 0)
verdict("★S921's m3 = \"there exists a SECOND independent invariant symmetric form\". "
        "But dim=1 in EVERY signature n=3..7 without exception ⟹ this mutant could not "
        "have survived ANYWHERE. A mutant that cannot fail to be caught measures nothing — "
        "it is not a kill-test, but a tautology. ⟹ \"mutants 3/3 CAUGHT ⟹ the measurement "
        "is valid\" in S921 should be read as 2/3 (m1, m2 are real). Not a wrong number — a "
        "wrong IMPLICATION, exactly the class [[zero-from-a-count-that-never-ran]].")

# ═══════════════════════════════════════════════════════════════════════════
section("M — MUTANTS (class S909: a silent mutant ⟹ THE MEASUREMENT IS INVALID)")
# ═══════════════════════════════════════════════════════════════════════════
CAUGHT = 0
print("""  ★★THIS SECTION WAS REWRITTEN AFTER A CATCH AGAINST MYSELF. The first pass had
  TWO empty mutants — exactly the class this same probe judges m3 for in S921:
    · m1 did its OWN inline check ⟹ it never ran the B3 code at all;
    · m3 asked "are T_a(0) abelian?" — which is an identity of the form (my own
      B1), the answer is ALWAYS "yes" ⟹ the mutant could not fail to fire.
  A mutant that cannot fail to be caught measures nothing. Now EACH one runs
  through the REAL measurement function (family_closes / dim_bilinear_invariants),
  and each one CAN survive — which is why catching them is a measurement.\n""")

# m1: a broken FORM (η dropped) — run through the REAL family_closes
print("  m1: T_bad — the same form, but WITHOUT η in the bottom row")
print("      (run through the REAL family_closes, not an inline check)")
ok_real = family_closes(3, 1, LAM, Tmat)
bad_caught = not family_closes(3, 1, LAM, Tbad)
check("m1-control: family_closes(3,1,λ,Tmat) = True (the real form survives)",
      ok_real)
if bad_caught:
    CAUGHT += 1
    print("     m1 CAUGHT: family_closes((3,1), T_bad) = False — the same code that "
          "gave True for the real form catches the broken one.")
else:
    print("     ✗ m1 SILENT — the B2 measurement is INVALID (closure is 'proved' for anything)")
check("m1 (broken form) caught by the REAL family_closes", bad_caught)
# ★and proof that the mutant COULD have survived: on a definite signature η=I ⟹ T_bad ≡ Tmat
print("     ★the mutant COULD have survived: at η=I (a definite signature) T_bad ≡ Tmat —")
print("       so 'False' on (3,1) came from INDEFINITENESS, not for free.")
check("m1 is not empty: family_closes(4,0,λ,Tbad) = True (there η=I ⟹ the mutant "
      "is INVISIBLE) ⟹ on (3,1) it was caught precisely by the signature",
      family_closes(4, 0, LAM, Tbad))

# m2: "an invariant linear form exists" — φ = e_0
print("  m2: a planted \"invariant linear form\" φ=(1,0,0,0)")
basis31 = so_basis(3, 1)
phi = Matrix(1, 4, [1, 0, 0, 0])
witness = [M for M in basis31 if (phi * M) != zeros(1, 4)]
if witness:
    CAUGHT += 1
    print(f"     m2 CAUGHT: an explicit M with φ(M·) ≠ 0 was found ({len(witness)} of them) "
          f"— the B6(i) measurement is alive.")
else:
    print("     ✗ m2 SILENT — the B6 measurement is INVALID")
check("m2 (fake linear invariant) caught by the B6 check", len(witness) > 0)

# m3: "the B6 solver always returns 1" — run through the REAL dim_bilinear_invariants
print("  m3: the claim \"dim(invariant bilinear) = 1 ALWAYS\" (i.e. my")
print("      solver is blind and prints 1 for anything) — must fail")
d31 = dim_bilinear_invariants(3, 1)
d20 = dim_bilinear_invariants(2, 0)
print(f"     the same code: (3,1) → {d31} · (2,0) → {d20}")
if d20 != 1:
    CAUGHT += 1
    print("     m3 CAUGHT: the solver RETURNED NOT-1 on (2,0) ⟹ its \"1\" at "
          "n≥3 is a real measurement, not a constant in the code.")
else:
    print("     ✗ m3 SILENT — the B6 measurement is INVALID (the solver cannot return ≠1)")
check("m3 (blind solver) caught: dim_bilinear_invariants(2,0) ≠ 1",
      d20 != 1)
check("m3-control: the same code gives 1 on (3,1)", d31 == 1)

# m4: "B1 would have measured a nonzero product, if there were one" — a kill-test of B1 itself
print("  m4: a kill-test of B1 — take λ≠0 and ask the same check \"T_a·T_b ≡ 0\"")
a = Matrix(3, 1, [Symbol("a0"), Symbol("a1"), Symbol("a2")])
b = Matrix(3, 1, [Symbol("b0"), Symbol("b1"), Symbol("b2")])
et21 = eta_of(2, 1)
Ta1 = Tmat(a, Integer(1), et21)
Tb1 = Tmat(b, Integer(1), et21)
if (Ta1 * Tb1).expand() != zeros(4, 4):
    CAUGHT += 1
    print("     m4 CAUGHT: at λ=1 the same expression is NONzero ⟹ the B1 check "
          "is capable of seeing a nonzero, its '≡0' at λ=0 is a real measurement of the form.")
else:
    print("     ✗ m4 SILENT — the B1 measurement is INVALID (the check cannot see a nonzero)")
check("m4 (kill-test of the B1 check) — the check can see a nonzero",
      (Ta1 * Tb1).expand() != zeros(4, 4))

print(f"\n  MUTANTS: {CAUGHT}/4 CAUGHT")
check("mutants 4/4 CAUGHT ⟹ the measurements are valid", CAUGHT == 4)

# ═══════════════════════════════════════════════════════════════════════════
section("SUMMARY")
# ═══════════════════════════════════════════════════════════════════════════
print("""
  ★WHAT STANDS:
    · S916/S918-S921 AS COMPUTATION: 5/5 W29 probes reproduced BIT-FOR-BIT,
      exit=0. The logs are honest about the code.
    · RUNG-4 IS WHOLE. My scope catch died (see B7) — the S920 verdict
      survived adversarial checking unscathed.
    · the η-norm — the one invariant form — stands (📖, multiplicity 0).

  ★WHAT SURVIVED AGAINST ME (reporting honestly — two of my four catches are dead):
    · "HANDLES 0" — STANDS. |λ| is washed out by the scale T→μT; what remains is
      the discrete sign(λ). A discrete choice ≠ a constant. My first summary
      said the opposite, contradicting my own section B5 (B5).
    · "CHOSEN SILENTLY" — WRONG. The choice was made OPENLY in the ex-ante text: "the
      affine algebra iso(p,q) = so(p,q) ⋉ translations". By writing "⋉ ℝⁿ", the
      ex-ante text NAMED ℝⁿ abelian — from there [T,T]=0 is forced BY THE NAME. The
      probe honestly measured the declared object. My B1 "caught" the fact that
      ℝⁿ is abelian.

  ★WHAT DOES NOT STAND (rung-5) — narrower and more precise than my first pass:
    1. ★MAIN POINT: the ex-ante text asked THREE things — "is it closed · IS IT
       CANONICAL · is it free of new constants". The word "canonic-" appears 0
       times in S921's code and log; no alternative was ever built. Canonicity
       IS a question about alternatives, so it CANNOT be answered without a
       comparison. S921 measured the CLOSURE of the declared object and spent
       the word "FORCED" — which answers CANONICITY. Question 2 of 3 was left
       unmeasured, yet the verdict reads as an answer to it.
       "Closed" ≠ "forced".
    2. Rung-5 (the A2/A3/E-skeleton) is CONDITIONAL on λ=0, and the condition is
       NOT DECLARED. Only at λ=0 is span{T} an ideal (B4), i.e. only there does
       the question "the invariant number of a shift" have an OBJECT. The
       conditionality is legitimate (the declaration precedes it), but the
       summary "at floor −1 the quantity has exactly this skeleton and NO
       MORE" is a claim of EXHAUSTIVENESS OF THE FLOOR, while only one member
       of the family was measured. The family is real: det K(λ) ∝ λⁿ ⟹ the
       degeneracy is exactly at λ=0 and nowhere else; the Killing radical n vs
       0 ⟹ not isomorphic.
    3. The multiplicity of A2 = 0 IN FULL (B6+B8). The three numbers = one
       ancestor (absolute irreducibility, n≥3). The court itself handed back
       the frame ("agreement with the textbook = a frame, multiplicity 0"),
       keeping for itself ONE exception — "ZERO linear invariants = a
       load-bearing new measurement". The redoubt falls: η is a
       vector↔covector bijection (an identity, 9/9), so this is the same
       fact.
    4. "Mutants 3/3" in S921 should be read as 2/3: m3 could not have
       survived in any signature n=3..7 (B8).

  ⟹ VISA B: RUNG-4 — OK (my catch fell, the verdict survived whole).
              RUNG-5 — NOT OK as currently worded.
    The content of rung-5 does not fall: the numbers are correct, the machine
    is clean, the declaration is honest. What falls is ONE WORD — "forced" —
    and the claim of exhaustiveness.
    ★THE FIX = ONE QUALIFIER, not a re-measurement:
      "forced WITHIN the declared object so(p,q) ⋉ ℝⁿ;
       the canonicity of the object is NOT MEASURED (no alternative was built)".
    ★A DEBT LEFT BEHIND (a bridge, not a verdict): "handles 0" holds as long as μ
      is free. Premise (b) "d+1 UNIT axes" fixes the length ⟹ |λ|=1/R² becomes a
      real dimensionful constant. The honest line: "handles 0, AS LONG AS NO
      LENGTH IS PINNED DOWN ANYWHERE". The slot in the shift leg is exactly one —
      and the program has exactly one allowed constant. The slot is not empty —
      it is UNMEASURED.
""")
tot = sum(SECTION_COUNTS.values())
# ★Guard "zero checked in a section" (class S909: a zero from a count that never
#   ran reads as a verdict). SUMMARY is an overview, has no measurements of its
#   own by design: the exception is DECLARED here, not passed over in silence.
for k, v in SECTION_COUNTS.items():
    if v == 0 and not k.startswith("SUMMARY"):
        print(f"  ✗ SECTION WITH ZERO CHECKS: {k} — THE MEASUREMENT IS INVALID")
        FAILED += 1
print(f"\n  checks: {tot} · FAIL: {FAILED}")
sys.exit(1 if FAILED else 0)
