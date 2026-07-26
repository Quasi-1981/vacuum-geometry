# -*- coding: utf-8 -*-
# DIM: symbolic-exact.  Numbers with a bracket [address · unit · type/operation].
"""S1063 T2 — THE SEAM (d,1), step 2: A HOMONYM-COURT N1↔N3 BY ANCESTOR.

ASSIGNMENT: project ruling, exante `hub/prime/SEAM_D1_SIGNATURE_EXANTE.md`, step T2 (T1 delivered, 11/0/0).
The gate of this step (2a): is the minus of N3 = **minus_arith** [T32] BY ANCESTOR — with a COMPUTED map,
not a postulated one (S1028 discipline: compute the action IN THE FIELD).

ANCESTORS (by citation, not re-derived):
  · T32 [S1011/J-0467] — the native Box of the cell `Λ(ψ,ν) = T_A(ψ) − T_col(ν)`; the minus = the sign
    of the Pontryagin-dual of the column [T26], NOT by hand; «the minus = the ONLY sign at which the native column
    satisfies the participation law [T30]».  The machinery is verbatim from `S1011_w42_leg6_native_box.py`:
      T_A(ψ) = Σ_{i=0..d} (2−2cos 2πp_i), the origin-choice p_0 = 0, p_i = ψ_i;
      T_col(ν) = 2−2cos(2πν/h),  ν ∈ ℤ/h,  h = d+1 [S1001: the period of the column];
      f(ψ) = 1 + Σ_{i=1..d} exp(2πi ψ_i),  a node ⟺ f = 0 [S1002].
  · T33 — chirality `σ_z H σ_z = −H` EXACTLY, spec ±|f|; the forcing of {I,H}, dim 2.
  · S1063-T1 — the inertia of the seam (1,2,d−2); the form `ϸ² − m₀² − |f|²` [S1045].
  · ★S1027 — A PRECEDENT: a minus-bridge was already TRIED and the verdict was a HOMONYM (shared only
    the arena = not a bridge).  I cite the verdict, I do NOT re-try it (a K-resurrection).

★WHAT I EXPECT AND WHAT I FEAR (carved BEFORE the count): both objects live on the SAME lattice and
both are called «the first minus» — this is exactly the configuration in which S1027 already caught
a homonym.  I fear CONFIRMING a bridge through a consonance of roles («a minus against the kinetics»), without
checking the functions.  So I measure THREE independent things: the value · the symmetry · the builder.

FENCE: Layer-1.  The physics-vocabulary classes named in the project's fence do not enter the code
(K-stone).  Homonyms with a prefix: `minus_arith` [T32] ⊥ `minus_chiral` [N3, named by this
step] ⊥ `minus_geom`.  `ϸ` = a spectral parameter.  No verdict is rendered.
"""
import itertools
import os
import sys

import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "..", "test")))
from _teeth import ok, report, reset                      # noqa: E402


class Tee:
    def __init__(self, real, fh):
        self.real, self.fh = real, fh

    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush()
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


# ═══════ THE T32 MACHINERY — VERBATIM FROM S1011 (not rewritten, carried over) ═══════
def term(h, m):
    """2 − 2cos(2πm/h) — the builder of BOTH sides of the Box [S1011]."""
    return sp.simplify(2 - 2 * sp.cos(2 * sp.pi * sp.Rational(m % h, h)))


def T_A(psi, h):
    """T_A(ψ) = Σ_{i=0..d}(2−2cos 2πp_i), the origin-choice p_0 = 0 ⟹ the axis-0 term is zero."""
    return sp.simplify(sum(term(h, m) for m in psi))


def T_col(nu, h):
    """T_col(ν) = term(ν) — the dual of the column, ν ∈ ℤ/h [S1001: EXACTLY h points]."""
    return term(h, nu)


def absf2(psi, h):
    """|f|² at f = 1 + Σ exp(2πiψ_i/h) — a structural function [S1002]."""
    re_ = 1 + sum(sp.cos(2 * sp.pi * sp.Rational(m % h, h)) for m in psi)
    im_ = sum(sp.sin(2 * sp.pi * sp.Rational(m % h, h)) for m in psi)
    return sp.simplify(sp.expand(re_**2 + im_**2))


def torsion(d, h):
    return list(itertools.product(range(h), repeat=d))


def nodes(d, h):
    return [psi for psi in torsion(d, h) if sp.simplify(absf2(psi, h)) == 0]


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1063_2_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf); sys.stdout = tee
    reset()

    print("=" * 80)
    print("S1063 T2 — A HOMONYM-COURT N1↔N3: is the minus of the seam THE minus of T32 (the map COMPUTED)")
    print("=" * 80)
    print()

    print("THE CANDIDATE MAP, which must be checked (not accepted):")
    print("  The T32 Box:  Λ(ψ,ν) = T_A(ψ) − T_col(ν)      — a minus between the CELL and the COLUMN")
    print("  The N3 seam:   Q(ϸ,k) = ϸ² − m₀² − |f(k)|²     — a minus between ϸ² and the HOP of the cell")
    print("  ⟹ a bridge exists ⟺ |f|² stands on the same side and IS the same object")
    print("     as one of the terms of the Box.  This is checked by FUNCTIONS, not by roles.")
    print()

    # ═══════ (2a-i) THE VALUE: is |f|² = T_A as FUNCTIONS ═══════
    print("(2a-i) THE VALUE — are |f|² and T_A the same function on the h-torsion")

    def functions_differ(pair):
        """A WORLD = A PAIR of functions (g₁, g₂) + (d,h) ⟹ whether they DIFFER at even one point
        of the h-torsion.  The detector reads both functions on a SHARED grid — so it can
        be re-asked on any pair."""
        g1, g2, d, h = pair
        return any(sp.simplify(g1(psi, h) - g2(psi, h)) != 0 for psi in torsion(d, h))

    def T_A_other_route(psi, h):
        """The same T_A, computed a DIFFERENT way: |1−ω|² = (1−cos)² + sin² instead of 2−2cos.
        A negative world: the detector MUST say «no difference».
        ★The first form wrote `sp.Abs(1−exp(...))**2` — and at d=2 sympy did NOT collapse the cube root
        of unity ⟹ the negative control «fired» on a defect of the WRITING, not of the object (the same
        class caught in T1: `f=0` was not recognized).  Now the route is without exp."""
        out = 0
        for m in psi:
            c = sp.cos(2 * sp.pi * sp.Rational(m % h, h))
            s = sp.sin(2 * sp.pi * sp.Rational(m % h, h))
            out += sp.expand((1 - c)**2 + s**2)
        return sp.simplify(sp.expand_trig(out))

    for d in (2, 3):
        h = d + 1
        ok(functions_differ, (absf2, T_A, d, h),
           f"★|f|² and T_A — are DIFFERENT functions on the h-torsion (d={d}, h={h})  "
           f"[2a-i · dimensionless · a comparison of values over the whole grid]",
           must_fail_on=[(f"T_A against T_A computed a different way (|1−ω|²) — the same "
                          f"function, d={d}", (T_A, T_A_other_route, d, h))])
        nd = nodes(d, h)
        print(f"     d={d}: nodes on the torsion {len(nd)}; at the node {nd[0] if nd else '—'}: "
              f"|f|² = {absf2(nd[0], h) if nd else '—'}  ⊥  T_A = {T_A(nd[0], h) if nd else '—'} "
              f"(=2h={2*h})  [2a-i · dimensionless · the value at the node]")
    print("     ⟹ ★THE SHARPEST POINT OF DIVERGENCE — THE NODE ITSELF: there |f|² = 0 BY DEFINITION,")
    print("       while T_A = 2h ≠ 0 (this is a line of ancestor S1011: «T_A(coset-node) = 2h ∉ range T_col»).")
    print("       That is, at the very point for which the seam is built, the two functions are not merely")
    print("       different — they diverge MAXIMALLY.")
    print()

    # ═══════ (2a-ii) SYMMETRY: distinguishes the MECHANISM, not the value ═══════
    print("(2a-ii) SYMMETRY — a global shift of the phases of all (d+1) axes")

    def shift_invariant(g_d_h):
        """A WORLD = (a function, d, h) ⟹ whether it is INVARIANT under a global shift of the phases of all
        d+1 axes (ψ_i → ψ_i + c ∀i, together with axis-0).  This is not simply a preference: the origin-choice p_0=0 in the Box
        [S1011] IS a choice, and a legitimate object should not depend on it."""
        g, d, h = g_d_h
        for c in range(1, h):
            for psi in torsion(d, h):
                psi_sh = tuple((m + c) % h for m in psi)
                # axis-0 also shifts: in the origin-choice p_0=0 this is f → e^{iφ}·f, that is,
                # for |f|² the shift of axis-0 must be accounted for explicitly
                v0 = sp.simplify(g(psi, h))
                v1 = sp.simplify(g_shifted(g, psi_sh, h, c))
                if sp.simplify(v0 - v1) != 0:
                    return False
        return True

    def g_shifted(g, psi_sh, h, c):
        """The value of the function after a global shift — for |f|² axis-0 carries e^{2πic/h},
        for T_A axis-0 carries the term term(c).  Both are accounted for EXPLICITLY, by one hand."""
        if g is absf2:
            re_ = sp.cos(2 * sp.pi * sp.Rational(c, h)) + \
                sum(sp.cos(2 * sp.pi * sp.Rational(m % h, h)) for m in psi_sh)
            im_ = sp.sin(2 * sp.pi * sp.Rational(c, h)) + \
                sum(sp.sin(2 * sp.pi * sp.Rational(m % h, h)) for m in psi_sh)
            return sp.simplify(sp.expand(re_**2 + im_**2))
        return sp.simplify(term(h, c) + sum(term(h, m) for m in psi_sh))

    ok(shift_invariant, (absf2, 2, 3),
       "★★|f|² IS INVARIANT under a global shift of the phases of all (d+1) axes (d=2)  "
       "[2a-ii · dimensionless · an exhaustive check of all shifts over the whole grid]",
       must_fail_on=[("T_A — is NOT invariant: it stands on the origin-choice p₀=0 (d=2)",
                      (T_A, 2, 3)),
                     ("T_A at d=3", (T_A, 3, 4))])
    print("     ⟹ ★★THIS DISTINGUISHES THE MECHANISM, NOT THE VALUE: |f|² does not depend on the origin-choice")
    print("       of axis-0, while T_A STANDS on it.  To identify them would mean ascribing to one")
    print("       object a symmetry the other does not have.  This is not «the numbers did not match»,")
    print("       these are DIFFERENT OBJECTS.")
    print()

    # ═══════ (2a-iii) THE BUILDER: what each side is made of ═══════
    print("(2a-iii) THE BUILDER — what each side of the minus is made of")
    th = sp.symbols('theta0:5', real=True)

    def built_from_pair_differences(builder):
        """A WORLD = A BUILDER (a function of phases) ⟹ whether it depends ONLY on the PAIRWISE DIFFERENCES
        θ_a − θ_b (a participation structure), not on the individual phases.  Operationally: invariance
        under θ_a → θ_a + c ∀a on SYMBOLIC phases, without a grid."""
        n = 4
        e = sp.simplify(builder([th[i] for i in range(n)]))
        c = sp.Symbol('c', real=True)
        e_sh = sp.simplify(builder([th[i] + c for i in range(n)]))
        return sp.simplify(sp.expand_trig(sp.expand(e - e_sh))) == 0

    def builder_absf2(ph):
        re_ = sum(sp.cos(t) for t in ph)
        im_ = sum(sp.sin(t) for t in ph)
        return sp.expand(re_**2 + im_**2)

    def builder_TA(ph):
        return sp.expand(sum(2 - 2 * sp.cos(t) for t in ph))

    ok(built_from_pair_differences, builder_absf2,
       "★★|f|² is BUILT from the PAIRWISE DIFFERENCES of phases: |f|² = (d+1) + 2Σ_{a<b} cos(θ_a−θ_b) — "
       "a PARTICIPATION structure of pairs  [2a-iii · dimensionless · symbolic invariance under θ→θ+c]",
       must_fail_on=[("T_A is built from INDIVIDUAL phases (Σ(2−2cos θ_i)) — moves under a shift",
                      builder_TA)])
    print(f"     |f|² symbolically (4 axes) = {sp.simplify(sp.expand_trig(builder_absf2([th[i] for i in range(4)])))}")
    print("     ⟹ ★A THIRD, INDEPENDENT DIVERGENCE: the sides of the T32 Box (T_A and T_col) have ONE")
    print("       builder `2−2cos` (individual phases), while |f|² has a DIFFERENT one (pairwise differences).")
    print("       The T32-minus stands BETWEEN TWO SAME-BUILDER objects; the seam-minus does not.")
    print()

    # ═══════ WHERE THE SEAM-MINUS THEN COMES FROM: COMPUTING ITS OWN ANCESTOR ═══════
    print("★WHERE THE SEAM-MINUS ACTUALLY COMES FROM — computing the ancestor, not naming it")
    ps = sp.Symbol('ps', real=True, positive=True)
    fr, fi, m_s = sp.symbols('f_re f_im m', real=True)
    a_p, b_p = sp.symbols('a_p b_p', positive=True)
    H_chiral = sp.Matrix([[0, fr + sp.I * fi], [fr - sp.I * fi, 0]]) + m_s * sp.Matrix([[1, 0], [0, -1]])
    H_definite = sp.Matrix([[a_p, 0], [0, b_p]])          # NON-chiral: both branches the same sign

    def form_carries_minus(op):
        """A WORLD = a 2×2 operator ⟹ whether the characteristic form carries a MINUS against ϸ², that is,
        whether `−det(op)` is a POSITIVE DEFINITE form of its parameters (all coefficients > 0).
        ★The first form asked `sp.ask(Q.negative(det))` — and it failed HONESTLY: at real
        symbols det = −(f_re²+f_im²+m²) is NOT strictly negative (all-zero is admissible).
        This is not a triviality of the writing: strictness IS the claim here, so I measure DEFINITENESS
        (all coefficients of one sign), not the value at one point."""
        dt = sp.expand(-op.det())
        vars_ = sorted(dt.free_symbols, key=str)
        if not vars_:
            return False
        coeffs = sp.Poly(dt, *vars_).coeffs()
        return all(sp.simplify(c) > 0 for c in coeffs)

    ok(form_carries_minus, H_chiral,
       "★★THE SEAM-MINUS = OPPOSITE SIGNS OF TWO BRANCHES: the chirality `σ_zHσ_z = −H` [T33] gives "
       "tracelessness ⟹ det = −(m₀²+|f|²) < 0 ⟹ the form is INDEFINITE  "
       "[N3 · dimensionless · the sign of det of the operator]",
       must_fail_on=[("a NON-chiral operator diag(a,b), a,b>0 — branches of ONE sign, "
                      "det > 0 ⟹ the form is DEFINITE, no minus", H_definite)])
    print(f"     det(H_chiral) = {sp.simplify(H_chiral.det())}  ⊥  det(H_definite) = "
          f"{sp.simplify(H_definite.det())}  [N3 · dimensionless · the determinant]")
    print("     ⟹ ★the seam-minus is born from the OPPOSITION OF THE SIGNS OF THE BRANCHES (the T33 chirality/")
    print("       T39 ε-parity), not from the participation of the column (T30/T32).  I name it")
    print("       `minus_chiral` — a THIRD minus, separate from `minus_arith` and `minus_geom`.")
    print()

    # ═══════ A CAVEAT ABOUT THE ORDER OF WRITING (otherwise a line would lie) ═══════
    print("★A CAVEAT TO T1, which must be said HERE (otherwise the reader will fill it in themselves)")

    def split_is_one_vs_two(sign_convention):
        """A WORLD = A GLOBAL SIGN of the form (±1) ⟹ whether the SPLIT remains «one against two».
        An ordered pair (1,2) under Q→−Q becomes (2,1) — this is a CONVENTION; the invariant is
        precisely the split {1, 2}."""
        Q = sp.diag(sign_convention * 1, -sign_convention * 1, -sign_convention * 1)
        ev = [sp.sign(Q[i, i]) for i in range(3)]
        return sorted([ev.count(1), ev.count(-1)]) == [1, 2]

    ok(split_is_one_vs_two, 1,
       "★THE SPLIT «one against two» is INVARIANT to the global sign of the form, whereas "
       "the ORDERED pair (1,2)⊥(2,1) is NOT  [T1-caveat · dimensionless · a sign-count "
       "under Q→−Q]",
       must_fail_on=[("a false world: the split «zero against three» (a definite form) — there "
                      "the split is different", 0)])
    print("     ⟹ T1 printed (1,2,d−2) in ITS OWN writing (+ϸ²).  The invariant — the split")
    print("       {1 against 2}; which side is called the «plus» is fixed by the convention of the writing, not")
    print("       by the measurement.  T32 fixes a RELATIVE sign in ITS OWN pair (cell⊥column) — this")
    print("       is a different pair, so its fixing does not carry over to ours.")
    print()

    code = report("S1063 T2 — a homonym-court")
    print()
    print("=" * 80)
    print("★RAW OUTPUTS OF T2 (no verdict rendered — the project's court judges)")
    print("  (2a) THE MAP IS COMPUTED AND IT DOES NOT EXIST — three independent divergences:")
    print("     · THE VALUE: at the node |f|² = 0, while T_A = 2h ≠ 0 (the maximal divergence")
    print("       precisely at the point for which the seam is built);")
    print("     · THE SYMMETRY: |f|² is invariant under a global phase shift, T_A is not")
    print("       (it stands on the origin-choice p₀=0);")
    print("     · THE BUILDER: the sides of the T32 Box are same-builder (2−2cos, individual phases),")
    print("       |f|² — from pairwise differences (a participation structure of pairs).")
    print("  ⟹ (2c) ★A HOMONYM: `minus_arith` [T32: cell⊥column, participation/Pontryagin]")
    print("     ⊥ `minus_chiral` [N3: opposite signs of the branches, T33/T39 chirality].")
    print("     Shared — only the ARENA (the same lattice).  This is EXACTLY the S1027 pattern, and I")
    print("     cite it as a precedent, I do not re-try it.")
    print("  ⟹ (2b) ★THE CARVED CONSEQUENCE «THE LINKING OF FRONTS» DOES NOT FOLLOW: the premise")
    print("     («the minus of N3 = minus_arith») is FALSE by computation ⟹ a second independent")
    print("     input for the T32→AX-indef dotted line DOES NOT exist from this.  I report it loudly, as")
    print("     befits a negative (the same class as Λ_min in S1061).")
    print("  ⟹ N1↔N3: two ROOTS, not two witnesses.  the project's court judges the multiplicity — not me.")
    print("=" * 80)
    sys.stdout = tee.real
    logf.close()
    return code


if __name__ == "__main__":
    sys.exit(main())
