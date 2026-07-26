# -*- coding: utf-8 -*-
# DIM: symbolic-exact.  Numbers with a bracket [address · unit · type/operation].
"""S1063 T3 — 2↔4-COMPONENT STITCHING: is the 4-component-ness NATIVE, or an effective splice. KILL-FIRST.

ASSIGNMENT: project ruling, exante `SEAM_D1_SIGNATURE_EXANTE.md` step T3 + the T2 verdict («GO on T3; it does NOT
depend on the fate of the bridge»).  The Dirac-remainder (2) from S1060.

★THE NULL HYPOTHESIS, CARVED IN THE EXANTE (I strike precisely at it, not build up the 4-component case):
   «4-components are NOT needed — the native Dirac of the seam is EXACTLY 2-component [the T33 forcing], and 4-components =
    an EFFECTIVE splice of a continuum-reading, the same class as i∂ₜ [S1060-T4-B]».
If it withstands — a rhyme with time (the continuum and the 4-components fail by one mechanism).
If it fails — the assignment demands a COUNT: where do EXACTLY two chiralities come from, and the divergence at d≠2
must be NAMED, not papered over.

ANCESTORS (by citation):
  · T33 — the invariant algebra {I,H} dim 2 ⟹ 2 components AT A POINT k; `σ_zHσ_z = −H` ⟹
    exactly two branches ±|f| (two chiralities at the node);
  · S1002 — COSET-NODES: ψ_i = (i+1)·j mod h, j = 1..d — exactly d of them;
  · T32/T26.7 — a bijection «the nodal set ↔ d nontrivial characters of the dual group ℤ/h»;
  · S1005/S1011 — the criterion of TIME-ACTIVITY of a character: sin(2πν/h) ≠ 0, that is, ν ∉ {0, h/2};
  · S1001 — h = d+1 (the period of the column), EXACTLY h points of the dual, not a choice;
  · S1045 — the node-set has dimension d−2 (at d=3 — a LINE, not points).

★WHAT I FEAR (carved BEFORE the count): there are THREE counts (all solutions on the torsion · coset-nodes ·
time-active), and they give DIFFERENT numbers.  The temptation is to take the one that gives 4, and call it native.
So I print ALL THREE and say WHICH ANCESTOR selects each; if 4 comes out only after two
selections — I say so, and do not present it as «forced».

FENCE: Layer-1.  The physics-vocabulary classes named in the project's fence do not enter the code
(K-stone).  Zero-tests via (Σcos, Σsin) — a requirement of the project's court after a recurrence of the writing-defect.
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


def f_is_zero(psi, h):
    """A zero-test via (Σcos, Σsin) — WITHOUT exp (an order of the project's court: an sympy-simplification is not
    the definition of zero; this class was caught twice this session)."""
    re_ = 1 + sum(sp.cos(2 * sp.pi * sp.Rational(m % h, h)) for m in psi)
    im_ = sum(sp.sin(2 * sp.pi * sp.Rational(m % h, h)) for m in psi)
    return sp.simplify(re_) == 0 and sp.simplify(im_) == 0


def all_nodes(d, h):
    """COUNT-1: all solutions f=0 on the h-torsion (the whole nodal variety)."""
    return [psi for psi in itertools.product(range(h), repeat=d) if f_is_zero(psi, h)]


def coset_nodes(d, h):
    """COUNT-2: COSET-nodes [S1002]: ψ_i = (i+1)·j mod h, i=0..d−1, j=1..d — d of them."""
    return [tuple(((i + 1) * j) % h for i in range(d)) for j in range(1, d + 1)]


def is_time_active(nu, h):
    """The criterion of TIME-ACTIVITY of a character [S1005/S1011]: sin(2πν/h) ≠ 0 ⟺ ν ∉ {0, h/2}."""
    return sp.simplify(sp.sin(2 * sp.pi * sp.Rational(nu % h, h))) != 0


def active_characters(d, h):
    """COUNT-3: time-active nontrivial characters of ℤ/h."""
    return [j for j in range(1, d + 1) if is_time_active(j, h)]


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1063_3_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf); sys.stdout = tee
    reset()

    print("=" * 80)
    print("S1063 T3 — 2↔4-COMPONENT: KILL-FIRST on the null hypothesis «4-component = an effective splice»")
    print("=" * 80)
    print()

    # ═══════ STEP 1: 2 COMPONENTS AT A POINT — AN ANCESTOR, NOT A DERIVATION ═══════
    print("STEP 1 — how many components AT A POINT k (ancestor T33, checking the carry-over)")
    fr, fi = sp.symbols('f_re f_im', real=True)
    H = sp.Matrix([[0, fr + sp.I * fi], [fr - sp.I * fi, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])

    def two_branches_opposite(op):
        """A WORLD = an operator ⟹ whether it has EXACTLY two branches of opposite sign
        (σ_z-chirality: σ_z·op·σ_z = −op ⟹ a spectrum ±)."""
        if sp.simplify(sz * op * sz + op) != sp.zeros(2, 2):
            return False
        ev = list(sp.Matrix(op).eigenvals().keys())
        return len(ev) == 2 and sp.simplify(sum(ev)) == 0

    ok(two_branches_opposite, H,
       "★the T33-carry-over: `σ_zHσ_z = −H` ⟹ EXACTLY two branches ±|f| at every point k — TWO "
       "chiralities at the node, not four  [step1 · dimensionless · the spectrum of the 2×2 chiral H]",
       must_fail_on=[("H + a diagonal (non-chiral) — the chirality breaks",
                      H + sp.Rational(1, 3) * sz),
                     ("the identity (branches of ONE sign)", sp.eye(2))])
    print("     ⟹ 2 components at a point — this is an ANCESTOR [the T33 forcing {I,H} dim 2], not my derivation.")
    print("       So 4 could appear ONLY as 2 × (however many INDEPENDENT nodes).")
    print()

    # ═══════ STEP 2: THREE NODE-COUNTS — ALL PRINTED, NONE CHOSEN ═══════
    print("STEP 2 — ★THREE COUNTS OF NODES, and they give DIFFERENT numbers (all three printed)")
    print("  count-1: all solutions f=0 on the torsion · count-2: coset-nodes [S1002] ·")
    print("  count-3: TIME-ACTIVE characters [S1005/S1011: ν ∉ {0, h/2}]")
    table = {}
    for d in (2, 3, 4, 5):
        h = d + 1
        n_all = len(all_nodes(d, h)) if d <= 4 else None
        n_cos = len(coset_nodes(d, h))
        act = active_characters(d, h)
        table[d] = (n_all, n_cos, len(act), 2 * len(act))
        print(f"     d={d} (h={h}): all solutions = {n_all if n_all is not None else '—'} · "
              f"coset-nodes = {n_cos} · TIME-ACTIVE = {len(act)} {act} ⟹ components "
              f"2×{len(act)} = **{2*len(act)}**  [step2 · dimensionless · node counts]")
    print()
    print("  ★I NAME THE DIVERGENCE OF THE COUNTS, I DO NOT PAPER OVER IT (an order of the assignment):")
    print("     at d=3 the three counts give 9 / 3 / 2 — and only the THIRD gives 4 components.")
    print("     Each selection has a NAMED ancestor: coset-nodes = S1002, time-activity")
    print("     = S1005/S1011 (the same criterion by which T32 counted time-active zeros).")
    print("     ⟹ «4» is NATIVE RELATIVE TO TWO NAMED SELECTIONS, and I say it exactly this way.")
    print()

    # ═══════ STEP 3: WHY AT d=3 EXACTLY ONE CHARACTER DROPS OUT ═══════
    print("STEP 3 — ★THE MECHANISM: why at d=3, of the three characters, exactly two are time-active")

    def self_conjugate_exists(d):
        """A WORLD = d ⟹ whether among the nontrivial characters of ℤ/h (h=d+1) there is a SELF-CONJUGATE one
        (ν = h/2, that is, ν ≡ −ν).  It is precisely this one, and only this one, that has sin = 0 ⟹ time-INACTIVE."""
        h = d + 1
        return any((2 * j) % h == 0 for j in range(1, d + 1))

    ok(self_conjugate_exists, 3,
       "★at d=3 (h=4) among the characters there IS a SELF-CONJUGATE ν=h/2=2 (ν ≡ −ν) — and it alone "
       "has sin(2πν/h)=0 ⟹ time-INACTIVE ⟹ drops out  [step3 · dimensionless · the existence of ν=h/2]",
       must_fail_on=[("d=2 (h=3) — h is ODD, no self-conjugate exists, BOTH are active", 2),
                     ("d=4 (h=5) — h is odd, all four are active", 4)])
    for d in (2, 3, 4, 5):
        h = d + 1
        sc = [j for j in range(1, d + 1) if (2 * j) % h == 0]
        print(f"     d={d} (h={h}): self-conjugate characters = {sc if sc else 'none'} "
              f"⟹ time-active {len(active_characters(d, h))} of {d}  "
              f"[step3 · dimensionless · the conjugation ν ↔ h−ν]")
    print("     ⟹ ★the mechanism is NAMED: active characters go in CONJUGATE PAIRS (ν, h−ν);")
    print("       a self-conjugate ν=h/2 exists ⟺ h is EVEN ⟺ d is ODD, and it drops out.")
    print("       So at d=3 exactly ONE pair remains = 2 ⟹ 2×2 = 4.")
    print()

    # ═══════ STEP 4: ★KILL-FIRST — THE VERDICT ON THE NULL HYPOTHESIS ═══════
    print("STEP 4 — ★KILL-FIRST: did the null hypothesis «4-component = an effective splice» withstand")

    def four_components_native(d):
        """A WORLD = d ⟹ whether the NATIVE count (2 × time-active characters) gives EXACTLY 4.
        This is precisely the strike against the null hypothesis: if 4 comes out WITHOUT any continuum-reading
        and without a manual choice — the «splice» hypothesis fails."""
        h = d + 1
        return 2 * len(active_characters(d, h)) == 4

    ok(four_components_native, 3,
       "★★THE NULL HYPOTHESIS FAILS AT d=3: the native count gives EXACTLY 4 components (2 branches × 2 "
       "time-active characters) — WITHOUT a continuum-reading and without a manual choice  "
       "[step4 · dimensionless · 2×|active| at d=3]",
       must_fail_on=[("d=4 (h=5) — there 2×4 = 8, not 4", 4),
                     ("d=5 (h=6) — there 2×4 = 8", 5)])

    def four_is_universal(dlist):
        """A WORLD = a list of d ⟹ whether «4 components» is a law ∀d (rather than a fact at d≤3)."""
        return all(four_components_native(d) for d in dlist)

    ok(lambda dl: not four_is_universal(dl), [2, 3, 4, 5],
       "★★BUT «4» IS NOT A ∀d-LAW: at d≥4 the native count gives 8 ⟹ 4-component-ness is a fact "
       "precisely at d ≤ 3 (h ≤ 4)  [step4 · dimensionless · 2×|active| over d=2..5]",
       must_fail_on=[("a list of only d=2,3 — there 4 comes out always, and the law would seem ∀d",
                      [2, 3])])
    print("     ⟹ ★★THE KILL-FIRST VERDICT: the null hypothesis FAILS — 4-components are NOT an effective splice")
    print("       of a continuum-reading: it comes out of the COUNT OF CHARACTERS, that is, from the dual of")
    print("       the COLUMN, and this is a native object (S1001/S1005), not an import.")
    print("     ⟹ ★BUT there is NO rhyme with time that the exante hoped for: the continuum fails ∀d,")
    print("       while 4-components hold only at d ≤ 3.  The two mechanisms are DIFFERENT — and I say this directly,")
    print("       because the rhyme would have been prettier than the truth.")
    print()

    # ═══════ STEP 5: THE BOUNDARY d+1 ≤ 4 — A NAMED OBSERVATION, NOT A CLAIM ═══════
    print("STEP 5 — the boundary at which this breaks (a named observation, NOT a claim)")

    def boundary_is_h_le_4(count_rule):
        """A WORLD = A COUNTING RULE for nodes ⟹ whether «gives exactly 4 components» ⟺ h = d+1 ≤ 4
        over the whole range d=2..5.
        ★The first form was EMPTY, and this is a RECURRENCE of my own class from T1: as a negative world I
        gave d=5, that is, THE SAME LAW at a different parameter — and the law holds there.
        A negative world must break the LAW, that is, be a DIFFERENT COUNTING RULE."""
        return all((2 * len(count_rule(d)) == 4) == ((d + 1) <= 4) for d in (2, 3, 4, 5))

    def rule_active(d):
        return active_characters(d, d + 1)

    def rule_all_characters(d):
        """A negative world: a count WITHOUT the time-activity criterion (all nontrivial characters).
        There «4 components» ⟺ 2d=4 ⟺ d=2, that is, the boundary would stand at h ≤ 3, and at d=3
        the equivalence with h≤4 BREAKS."""
        return list(range(1, d + 1))

    ok(boundary_is_h_le_4, rule_active,
       "★the boundary of 4-component-ness is EXACTLY h = d+1 ≤ 4 — the SAME boundary at which "
       "the node-mark closed [S1059B: «combinatorially forced at d+1 ≤ 4»]  [step5 · dimensionless · "
       "a coincidence of the boundary over d=2..5]",
       must_fail_on=[("a count WITHOUT time-activity (all characters) — there the boundary would stand at "
                      "h ≤ 3, and at d=3 the equivalence breaks", rule_all_characters)])
    print("     ⟹ ★AND THIS IMMEDIATELY SHOWS THE COST OF THE SELECTION: the boundary `h ≤ 4` holds PRECISELY BECAUSE")
    print("       the count goes over TIME-ACTIVE characters.  Without this ancestor the boundary would be")
    print("       `h ≤ 3` (that is, 4-components only at d=2).  The selection is not cosmetic — it")
    print("       is what carries the whole line about d=3.")
    print("     ⟹ ★I NAME THE COINCIDENCE AND NO MORE: the same boundary `d+1 ≤ 4` already stood in S1059B")
    print("       (the forcing of the decomposition).  Whether this is ONE ancestor, or two different mechanisms with")
    print("       the same boundary — I did NOT measure.  I do not declare a multiplicity: this is a named-question")
    print("       for the project's court, and presenting it as a witness would be exactly the inflation that")
    print("       our own registry forbids.")
    print()

    code = report("S1063 T3 — 2↔4-component")
    print()
    print("=" * 80)
    print("★RAW OUTPUTS OF T3 (no verdict rendered)")
    print("  A TABLE OF COUNTS (all three, none hidden):")
    print("     d | all solutions | coset-nodes | TIME-ACTIVE | components 2×act")
    for d in (2, 3, 4, 5):
        n_all, n_cos, n_act, comp = table[d]
        print(f"     {d} | {str(n_all) if n_all is not None else '—':>15} | {n_cos:>12} | "
              f"{n_act:>12} | {comp:>15}")
    print("  ★THE NULL HYPOTHESIS («4-component = an effective splice, the i∂ₜ class») — FAILS at d=2,3:")
    print("     4 = 2 branches [T33] × 2 time-active characters [S1005], zero continuum-reading.")
    print("  ★BUT: «4» is NOT a ∀d-law — at d≥4 the native count gives 8 ⟹ the rhyme «the continuum and the 4-component")
    print("     fail by one mechanism» is NOT confirmed: the continuum fails ∀d, the 4-component case")
    print("     holds only at d+1 ≤ 4.  The mechanisms are different.")
    print("  ★THE MECHANISM OF THE COUNT: characters go in conjugate pairs (ν, h−ν); a self-conjugate")
    print("     ν=h/2 exists ⟺ h is even ⟺ d is odd and it is time-INACTIVE ⟹ it drops out.")
    print("     At d=3: {1,2,3} → the pair {1,3} + the dropped {2} ⟹ exactly 2 ⟹ 4 components.")
    print("  ★HONESTLY ABOUT THE SELECTIONS: «4» is native RELATIVE TO TWO named ancestors (coset-nodes")
    print("     S1002 · time-activity S1005).  Without the second selection, at d=3 it would come out 6.")
    print("=" * 80)
    sys.stdout = tee.real
    logf.close()
    return code


if __name__ == "__main__":
    sys.exit(main())
