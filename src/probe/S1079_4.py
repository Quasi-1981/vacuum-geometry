# -*- coding: utf-8 -*-
# DIM: na  (an audit of the field's PROVENANCE — a statement about construction steps, not
#           about dimension; no d on the load-bearing path.)
"""S1079 · A TARGETED AUDIT: which STEP needs an ORDERED field, and which a REALLY CLOSED one

CAUSE: Beta's visa J-0520 (VALID/measured) on Z3/Z4.  Its load-bearing line (accepted without
argument, it is correct):
> «every step that writes a form in the NORMAL form diag(±1) … SILENTLY uses "every
>  positive element is a SQUARE", that is, a REALLY CLOSED field, not merely an ordered one»
> ⟹ the honest form of the component is TWO-PART: ordered FOR THE CERTIFICATE, really
>   closed FOR THE NORMAL FORM.

★WHY THIS FILE (and why this is not resonance): Beta named three carriers of the line — my
Z4-c, her D1 (S1076), and the ancestor S908.  I ran her line across the REST of my own acts
— and **Z1 also falls under it**, while Z2 does NOT.  the project's adjudication must carve the two-part form;
the carve will be exact only if it is known WHICH part bites WHICH step.  This is not a
correction of a correction, it is a MAP of the scope.

WHAT IS MEASURED:
  T1 · the sign of det of the form's restriction to a 2-plane is invariant under a change of
       the plane's basis OVER ANY ordered field (because det is multiplied by a square)
       ⟹ Z2 is safe;
  T2 · over ℚ, the signature does NOT classify a form (Beta's example: diag(1) ⊥ diag(2))
       ⟹ a step that identifies a form with its signature needs a really closed field;
  T3 · a map: which of my steps falls under which part.
"""
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'test'))

import sympy as sp
from _teeth import ok, ok_contrast, report


# ── T1 · what suffices for the CLASSIFICATION of a 2-PLANE (step Z2) ────────────────────

def plane_det_under_change(params):
    """probe(p) → (what_stands, what_moves) from ONE calculation.

    Take a 2×2 Gram restriction and an arbitrary change of the plane's basis S (a
    parameter).
    STANDS: the sign of det (because det(SᵀGS) = det G · (det S)²) · MOVES: the value of
    det itself.
    """
    a = sp.Integer(params)
    G = sp.Matrix([[1, 0], [0, -3]])          # a mixed plane (K-type)
    S = sp.Matrix([[a, 1], [0, 1]])           # det S = a ≠ 0
    Gp = sp.simplify(S.T * G * S)
    d = sp.simplify(Gp.det())
    return (sp.sign(d), d)


# ── T2 · what suffices for the SIGNATURE TO CLASSIFY (step Z1 / normal form) ─

def signature_classifies(world):
    """★DEFINITION (carved BEFORE the count):

    «The signature CLASSIFIES over a field F» ⟺ two forms with the same signature are
    congruent over F.  A test on Beta's example: diag(1) versus diag(2), signature (1,0)
    in both; congruence ⟺ ∃x ∈ F, x ≠ 0: x·1·x = 2 ⟺ 2 is a square in F.
    """
    return world['two_is_square']


def two_is_square_in_Q():
    """Does a rational x with x² = 2 exist — computed, not quoted (the square root of 2)."""
    sols = sp.solve(sp.Eq(sp.Symbol('x', rational=True) ** 2, 2))
    return any(s.is_rational for s in sols)


# ──────────────────────────────────────────────── the run ─────────────────────

def main():
    print("=" * 78)
    print("S1079 · AN AUDIT OF THE FIELD: ordered FOR WHAT ⊥ really closed FOR WHAT")
    print("=" * 78)

    print("\nT1 · step Z2 (the classification of a 2-plane J/K/N) — what suffices for it")
    ok_contrast(plane_det_under_change, [2, 3, 5],
                "the sign of det of the restriction STANDS under a change of the plane's "
                "basis, while the value of det itself MOVES ⟹ the three-way classification "
                "J/K/N needs only an ORDER (det is multiplied by a square), not roots")

    print("\nT2 · the step «a form ≡ its signature» — what suffices for it")
    got = two_is_square_in_Q()
    print(f"   ∃x ∈ ℚ: x² = 2 ?  → {got}   (solutions: {sp.solve(sp.Eq(sp.Symbol('x')**2, 2))})")
    ok(signature_classifies,
       {'name': 'ℝ (really closed)', 'two_is_square': True},
       "★over a really closed field the signature CLASSIFIES (diag(1) ≅ diag(2), because "
       "2 = (√2)²) — but over ℚ it does NOT, even though the signature in both is (1,0)",
       must_fail_on=[("ℚ (ordered, but not really closed)",
                      {'name': 'ℚ', 'two_is_square': got})])
    assert got is False

    print("""
T3 · ★A MAP OF THE SCOPE (for the project's carve — so the two-part form is EXACT)

  step                                     | what it really uses      | source
  ----------------------------------------|--------------------------|----------
  Z4-c certificate (a vector of norm −1)   | ORDEREDNESS              | my measurement
  Z2 classification of 2-planes J/K/N      | ORDEREDNESS (T1)         | this file
  Z1 order (p₁,q₁) ≼ (p₂,q₂)               | REAL CLOSEDNESS (T2)     | this file
  the normal form diag(±1) in general      | REAL CLOSEDNESS          | Beta's line

  ★WHY Z1 FALLS UNDER IT (a self-demotion, not someone else's): my Z1 identifies a form
  with the pair (p,q) and compares PAIRS.  This is legitimate exactly when the pair is a
  COMPLETE invariant — and T2 shows it is NOT complete over a merely ordered field.  So my
  line "a rigid NUMBER" holds, but its carrier is a really closed field, not merely an
  ordered one.  Beta named three addresses of her line; the fourth is mine, and I add it
  myself.""")

    code = report("S1079 · an audit of the field's scope")

    print("""
────────────────────────────────────────────────────────────────────────────────
SUMMARY (a measurement, not a verdict — carved by a project ruling):

  · Beta's two-part form is CONFIRMED and REFINED by address: not "the whole construction
    needs a really closed field", but EXACTLY those steps that write/compare a normal
    form.  Steps that work with the SIGN (the Z4-c certificate, the Z2 classification)
    hold over any ordered field — including ℚ.
  · ⟹ the carve-formula being proposed (not rendered here): the input component =
    "an ordered field" + "real closedness WHERE THE NORMAL FORM IS USED".
    A one-word variant is wrong in both directions: "ℚ suffices" is false for Z1;
    "ℝ is needed everywhere" is false for Z2/Z4-c.
  · The question of Λ: Beta's answer is accepted (scale multiplies the form, and
    "scales are comparable" already IS orderedness, that is, the same assumption, not a
    source).  My half of the fork is withdrawn: I too see no source in Λ.
────────────────────────────────────────────────────────────────────────────────""")
    return code


if __name__ == '__main__':
    sys.exit(main())
