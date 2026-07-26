# -*- coding: utf-8 -*-
# DIM: symbolic-exact.  Numbers with a bracket [address · unit · type/operation].
"""S1061 — T38-PRODUCTIVITY: does a split need to exist so the TACT DOES NOT LOSE DISTINCTION.

ASSIGNMENT: Omega, exante `hub/prime/T38_PRODUCTIVITY_EXANTE.md` (carved BEFORE the count).
The gate — an open line of S1060 («BOUNDARIES» item 3): T38 is measured for EXISTENCE, not for PRODUCTIVITY.
The author's word: «the conserved quantity of the break becomes productive only under a STABLE one-vector
splitting in time».

ANCESTORS (addresses, not paraphrases):
  · S1045 Component 1 — the chain `G = Σₙ G₀·(H·G₀)ⁿ`, `H=[[0,f],[f̄,0]]` (the T33-hop) ⟹
    ONE TACT = ONE factor `(H·G₀)`;
  · S1045 Component 3 — `G₀ = diag(1/(ϸ−m₀), 1/(ϸ+m₀))` = the SEGMENT WEIGHT (the address of the weight,
    which S1060-T4 did not build: there the shift was pure topology of the column);
  · S1060-T2 — the exponent `w=1/ϸ` IS the tact (the coefficient of `w^{n+1}` is exactly `Hⁿ`);
  · S1060-T3 — `kappa_conv = 1` (a bijection tact↔character, 0 handles);
  · S1060-T4-A — the dual of the tact is COMPACT (the phase ceiling does not grow with N);
  · S1060-D — at the node without a split `H≡0`, the commutant 2→4;
  · T38(3) `kappa_stab = Λ` (the SPLIT-equation 2κ=J(m₀)) · T39 (a split, an alternation period of 2).

★THE HOMONYM κ — DECLARED BEFORE THE COUNT (the exante carries a bare «κ=0» in P3, and it is two-valued):
  `kappa_conv` = the tact↔character converter [S1060-T3] = 1;
  `kappa_stab` = the stabilization stiffness = Λ [T38(3)] — THIS is the κ in P3 («κ=0 ⟹ a runaway m»).
  A bare `κ` is not used in this probe.

★WHAT I EXPECT AND WHAT I FEAR (carved BEFORE the count, in the probe's body): I fear that all three readings
will turn out fixed by the exante (a K-tautology) — and then «productivity» is empty. The second
fear is the opposite: that I will find a «ceiling» that comes from IMPORTING the continuum
form `U = e^{−iϸ}`, which T4-B precisely refused to derive. Both worlds are reported as they are.

FENCE: Layer-1. There is no bare «E» (ϸ = a spectral parameter, an E-base with a hyphen).
{Schrödinger, Matsubara, ħ, discretization-talk} — behind the fence. `m₀` = a representative of the regime `kappa_stab<a`,
not a canon-number (J-0486). Bond=1 [codex/lexicon.yaml].
"""
import os
import sys

import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "..", "test")))
from _teeth import ok, ok_contrast, report, reset          # noqa: E402


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


# ═══════════ FIELD 0(a): THE ONE-TACT OPERATOR — FROM AN ANCESTOR, NOT INVENTED ═══════════
def hop(f):
    """H = [[0,f],[f̄,0]] — the bipartite-hop of T33/S1045 (off-diagonal, chirality-odd)."""
    return sp.Matrix([[0, f], [sp.conjugate(f), 0]])


def seg_weight(ps, m0):
    """G₀ = diag(1/(ϸ−m₀), 1/(ϸ+m₀)) — the SEGMENT WEIGHT, S1045 Component 3.
    This is exactly the address the weight comes from: S1060-T4 did not build it."""
    return sp.diag(1 / (ps - m0), 1 / (ps + m0))


def U_tact(world):
    """A WORLD = (f, m₀, ϸ) ⟹ the ONE-TACT operator `U = H·G₀` [S1045: one factor of the chain].
    The tact-index n = the power of `w` [S1060-T2], kappa_conv=1 [S1060-T3]."""
    f, m0, ps = world
    return sp.simplify(hop(f) * seg_weight(ps, m0))


def shell_ps(f, m0):
    """ϸ on the shell: det(ϸ−H−m₀σ_z)=0 ⟹ ϸ² = m₀²+|f|²  [S1045: a two-band dispersion carrying a splitter]."""
    return sp.sqrt(m0**2 + sp.Abs(f)**2)


def world_on_shell(f, m0):
    return (f, m0, shell_ps(f, m0))


def eig_of(world):
    U = U_tact(world)
    return [sp.simplify(v) for v in sp.Matrix(U).eigenvals().keys()]


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1061_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf); sys.stdout = tee
    reset()

    print("=" * 80)
    print("S1061 — T38-PRODUCTIVITY: is a split needed so the tact does not lose distinction")
    print("=" * 80)
    print()

    f_gen = sp.Integer(2) + 3 * sp.I        # a generic f off the node [S1045: f(k)≠0]
    absf = sp.sqrt(13)

    # ═══════════════ FIELD 0 — MANDATORY (an anti-tautology gate) ═══════════════
    print("FIELD 0 — THE ONE-TACT OPERATOR AND ITS UNITARITY (an empty field = a refusal of the assignment)")
    print("  (0a) THE ADDRESS OF THE WEIGHT: `G₀ = diag(1/(ϸ−m₀), 1/(ϸ+m₀))` — S1045 Component 3")
    print("       («the splitter = the weight of the segment between zigzag-hops»); ONE TACT = one factor")
    print("       `(H·G₀)` — S1045 Component 1, the series `G = Σₙ G₀(H·G₀)ⁿ`.")
    print(f"       U(f={f_gen}, m₀=5, ϸ on the shell) = {U_tact(world_on_shell(f_gen, sp.Integer(5)))}")

    def not_unitary(world):
        """A WORLD ⟹ whether U†U ≠ I.  The detector reads the operator of the world ITSELF, not a ready-made boolean."""
        U = U_tact(world)
        return sp.simplify(U.H * U - sp.eye(2)) != sp.zeros(2, 2)

    ok(not_unitary, world_on_shell(f_gen, sp.Integer(5)),
       "★★(0b) BINARY: U is NOT unitary by construction — U†U = diag(|f|²/(ϸ−m₀)², |f|²/(ϸ+m₀)²), "
       "and with a split these two numbers are DIFFERENT  [0b · dimensionless · U†U−I on the shell m₀=5]",
       must_fail_on=[("the splitter-free shell (m₀=0, ϸ=|f|) — THERE U really is unitary",
                      world_on_shell(f_gen, sp.Integer(0)))])

    def is_unitary(world):
        U = U_tact(world)
        return sp.simplify(U.H * U - sp.eye(2)) == sp.zeros(2, 2)

    ok(is_unitary, world_on_shell(f_gen, sp.Integer(0)),
       "★unitarity holds EXACTLY on the splitter-free shell (m₀=0 and ϸ²=|f|²) — that is, "
       "the predicate is non-empty on BOTH sides  [0b · dimensionless · U†U at m₀=0]",
       must_fail_on=[("the same shell with a split m₀=5", world_on_shell(f_gen, sp.Integer(5))),
                     ("off the shell, splitter-free (ϸ=2|f|)", (f_gen, sp.Integer(0), 2 * absf))])
    print("  ⟹ (0c) SO P2 IS NOT AN EXERCISE: the question of growth has a subject, because unitarity")
    print("       is ABSENT by construction.  There will be no tautological PASS here — measured honestly.")
    print()

    # ═══════════════════════ P1 — NONTRIVIALITY ═══════════════════════
    print("P1 — NONTRIVIALITY: does the tact leave a RELATIVE trace (Δθ between the blades)")
    print("  ★HONESTY OF ANCESTRY (the exante, carved BEFORE the count): P1 stands on the VERY SAME")
    print("    ancestor as section D (the degeneracy at H≡0), read on a different property")
    print("    ⟹ the WITNESS MULTIPLICITY is 1, NOT 2.  This is a naming of the property, not a new witness.")

    def leaves_trace(world):
        """A WORLD ⟹ whether the spectrum of U carries a RELATIVE splitting of the blades (Δθ = arg λ₊ − arg λ₋ ≠ 0).
        A global phase marks nothing — only the difference marks anything."""
        ev = eig_of(world)
        if len(ev) < 2 or all(sp.simplify(v) == 0 for v in ev):
            return False
        args = sorted(sp.nsimplify(sp.arg(v)) for v in ev)
        return sp.simplify(args[-1] - args[0]) != 0

    ok(leaves_trace, (f_gen, sp.Integer(0), 2 * absf),
       "★the tact leaves a trace OFF the node and WITHOUT a split: Δθ = π ≠ 0  [P1 · a tact · the difference "
       "of the arguments of the eigenvalues of U]",
       must_fail_on=[("★AT THE NODE without a split (f=0, m₀=0) — H≡0 ⟹ U≡0, no trace",
                      (sp.Integer(0), sp.Integer(0), sp.Integer(7))),
                     ("★AT THE NODE WITH A SPLIT (f=0, m₀=5) — U≡0 TOO: a split does NOT restore the trace",
                      (sp.Integer(0), sp.Integer(5), sp.Integer(7)))])
    print("     ⟹ ★A NEW LINE AGAINST THE EXANTE'S EXPECTATION: at the node a split does NOT restore the trace")
    print("       of the tact (U = H·G₀ ≡ 0 at f=0 regardless of m₀), whereas in section D it")
    print("       restored the FORCING OF THE SLOT.  These are again DIFFERENT properties — and line D is not")
    print("       shaken by this: there the commutant was computed for (H+m₀σ_z), here — for the STEP per tact.")
    print()

    # ═════════ P1-addendum: THE ALPHABET OF PHASES (here new content, and it is structural) ═════════
    print("P1-addendum — WHICH PHASES PER TACT ARE POSSIBLE AT ALL (an alphabet, not values)")
    grid = [(f_gen, sp.Integer(0), 2 * absf),          # above the zone, split-free
            (f_gen, sp.Integer(5), sp.Rational(1, 2)),  # inside the split
            (f_gen, sp.Integer(5), shell_ps(f_gen, sp.Integer(5))),   # the shell with a split
            (f_gen, sp.Integer(0), absf)]               # the splitter-free shell

    def phases_quantized(weight_rule):
        """A WORLD = A WEIGHT RULE (a function (ϸ,m₀)→G₀) ⟹ whether all phases per tact fall
        into the finite alphabet {0, ±π/2, π}.  The detector rebuilds U with the world's RULE —
        so it can be re-asked on a different weight."""
        alphabet = {sp.Integer(0), sp.pi, -sp.pi / 2, sp.pi / 2, -sp.pi}
        for f, m0, ps in grid:
            U = sp.simplify(hop(f) * weight_rule(ps, m0))
            for v in sp.Matrix(U).eigenvals().keys():
                a = sp.nsimplify(sp.simplify(sp.arg(sp.simplify(v))))
                if a not in alphabet:
                    return False
        return True

    def weight_with_phase(ps, m0):
        """A negative world: the same weight, but with a phase on the segment (a hypothetical weight
        that is NOT real) — if the segment weight carried a phase, the alphabet would break."""
        return sp.exp(sp.I * sp.Rational(3, 10)) * seg_weight(ps, m0)

    ok(phases_quantized, seg_weight,
       "★★THE PHASE PER TACT IS DISCRETIZED: the segment weight is REAL (S1045) ⟹ the alphabet of phases is exactly "
       "{0, ±π/2, π} — zone/split/node, and NOTHING in between  [P1-addendum · a tact · the set of arg λ over 4 worlds]",
       must_fail_on=[("a hypothetical segment weight with a phase e^{i·3/10} — the alphabet breaks",
                      weight_with_phase)])
    for f, m0, ps in grid:
        ev = eig_of((f, m0, ps))
        print(f"     [f={f}, m₀={m0}, ϸ={ps}] λ = {ev}  ⟹ arg = "
              f"{[sp.nsimplify(sp.arg(v)) for v in ev]}, |λ| = {[sp.N(sp.Abs(v), 6) for v in ev]}")
    print("     ⟹ ★THE TACT IS A BIT-REGISTER, NOT A DIAL: a continuous phase per tact in this")
    print("       construction does NOT exist.  ★And this is the SAME fact as T4-B (there is no i∂ₜ as")
    print("       an identity): a continuous phase would be given by `U = e^{−iϸ}` — exactly the form that")
    print("       T4-B refused to derive.  ⟹ ANCESTOR MULTIPLICITY 1, declared (a K-ancestor).")
    print()

    # ═══════════════════════ P2 — BOUNDEDNESS ═══════════════════════
    print("P2 — BOUNDEDNESS: is growth over tacts possible at all")

    def squares_to_identity(world):
        """A WORLD ⟹ whether U² = I exactly (then growth is impossible BY FORM, ∀n)."""
        U = U_tact(world)
        return sp.simplify(U * U - sp.eye(2)) == sp.zeros(2, 2)

    ok(squares_to_identity, world_on_shell(f_gen, sp.Integer(5)),
       "★★ON THE SHELL U² = I EXACTLY, WITH A SPLIT TOO ⟹ ρ(U)=1, divergence is IMPOSSIBLE ∀n "
       "regardless of m₀; this is the ALTERNATION PERIOD 2 [T39] in tact form  "
       "[P2 · a tact · U²−I on the shell m₀=5]",
       must_fail_on=[("off the shell (ϸ=2|f|, m₀=0) — there U² ≠ I", (f_gen, sp.Integer(0), 2 * absf)),
                     ("off the shell inside the split (ϸ=1/2, m₀=5)", (f_gen, sp.Integer(5), sp.Rational(1, 2)))])
    print("     ⟹ ★THE READING «PRODUCTIVITY VIA BOUNDEDNESS» IS EMPTY, AND I REPORT IT AS A NEGATIVE:")
    print("       a split is NOT NEEDED for boundedness (U²=I even at m₀=0).  T38 does not work here")
    print("       either FOR or AGAINST — there is no subject.")

    def non_normal(world):
        """A WORLD ⟹ whether U is NON-normal (U†U ≠ UU†): ρ=1, but ‖U‖>1 — a distortion without growth."""
        U = U_tact(world)
        return sp.simplify(U.H * U - U * U.H) != sp.zeros(2, 2)

    ok(non_normal, world_on_shell(f_gen, sp.Integer(5)),
       "★A SPLIT BREAKS THE NORMALITY of U on the shell (U†U ≠ UU†): ‖U‖ = √((ϸ+m₀)/(ϸ−m₀)) > 1 while "
       "ρ(U)=1  [P2 · dimensionless · the commutator U†U−UU† on the shell]",
       must_fail_on=[("the splitter-free shell — there U is NORMAL (unitary), ‖U‖=1",
                      world_on_shell(f_gen, sp.Integer(0)))])
    chi2 = sp.simplify((shell_ps(f_gen, sp.Integer(5)) + 5) / (shell_ps(f_gen, sp.Integer(5)) - 5))
    print(f"     χ² = (ϸ+m₀)/(ϸ−m₀) at m₀=5, |f|=√13: {sp.nsimplify(chi2)} ≈ {sp.N(chi2, 6)}")
    print("     ⟹ ★THE CARVED RISK OF AN INVERSION WAS REALIZED IN A WEAKER FORM, and I say this directly:")
    print("       a split genuinely makes the step per tact NON-unitary — but this is a DISTORTION")
    print("       (a bounded χ), NOT growth: U²=I kills accumulation over tacts.")
    print()

    # ═══════════════════ P3 — INJECTIVITY (a step of the assignment) ═══════════════════
    print("P3 — INJECTIVITY: does `m₀ ↦ (a readout per tact)` remain distinguishing")
    print("  (the homonym is declared: here `kappa_stab=0` [T38(3)], NOT `kappa_conv` [S1060-T3])")

    def polarity(m0v):
        """p = m₀/ϸ on the shell — the RELATIVE weight of the blades per tact (the single channel through which m₀ enters)."""
        return sp.simplify(m0v / shell_ps(f_gen, m0v))

    def probe_mass(m0v):
        """ONE body, ONE object: the same U(m₀) on the shell yields BOTH quantities —
        (what stands still) the set of phases per tact  ⊥  (what moves) the polarity of the blades."""
        U = U_tact(world_on_shell(f_gen, m0v))
        ph = tuple(sorted(str(sp.nsimplify(sp.arg(sp.simplify(v))))
                          for v in sp.Matrix(U).eigenvals().keys()))
        return ph, sp.nsimplify(polarity(m0v))

    ok_contrast(probe_mass, [sp.Integer(0), sp.Rational(1, 10), sp.Integer(1), sp.Integer(10),
                            sp.Integer(1000)],
                "★★THE PHASE PER TACT IS BLIND TO m₀ (stands at {0,π} ∀m₀), whereas the POLARITY OF THE BLADES "
                "moves — the same machinery moves the control  [P3 · a tact · (arg λ) ⊥ (m₀/ϸ)]")
    print("     ⟹ ★★THIS IS THE ANSWER TO P3 BY TYPE, NOT BY MAGNITUDE: m₀ DOES NOT LIVE IN THE PHASE.")
    print("       There is no winding of the phase from m₀ — hence «winding ⟹ indistinguishability» also")
    print("       cannot occur.  The exante's positive world (κ_stab=0 ⟹ the phase winds)")
    print("       is NOT REALIZED, and the reason is structural, not numerical.")

    _RANGE = [sp.Integer(1), sp.Integer(10), sp.Integer(10**3), sp.Integer(10**6),
              sp.Integer(10**6) + 1]

    def injective_on_range(mapping):
        """A WORLD = A MAPPING m₀ ↦ a readout ⟹ whether it is injective on the range that
        admits an unstabilized regime (m₀ up to 10⁶+1).  ★The comparison is EXACT (the difference
        symbolically ≠ 0) — see the tooth below: any APPROXIMATE comparison introduces a THRESHOLD
        and collapses pairs, that is, it counterfeits «indistinguishability»."""
        vals = [mapping(m) for m in _RANGE]
        for i in range(len(vals)):
            for j in range(i + 1, len(vals)):
                if sp.simplify(vals[i] - vals[j]) == 0:
                    return False
        return True

    def winding_readout(m0v):
        """Negative world 1: a readout that WINDS (as it would if m₀ lived in the phase) —
        m₀ ↦ m₀ mod 1.  Here injectivity GENUINELY dies."""
        return sp.Mod(m0v, 1)

    def polarity_at_finite_precision(m0v):
        """★Negative world 2 — THIS TOOTH GREW OUT OF MY OWN FAILURE (the first run, ✗):
        THE SAME polarity, read with a THRESHOLD (8 significant digits) — and it collapses
        m₀=10⁶ with m₀=10⁶+1 into the same number.  That is, «indistinguishability» appears
        EXACTLY when a resolution threshold is introduced, and NOT BEFORE."""
        return sp.N(polarity(m0v), 8)

    ok(injective_on_range, polarity,
       "★★m₀ ↦ THE POLARITY OF THE BLADES IS INJECTIVE EXACTLY over the entire range that admits "
       "an unstabilized regime (m₀ = 1 … 10⁶+1) ⟹ THERE IS NO CEILING  "
       "[P3 · dimensionless · a pairwise difference of m₀/ϸ, exact]",
       must_fail_on=[("a wound readout (m₀ mod 1) — there injectivity dies",
                      winding_readout),
                     ("★THE SAME polarity with an 8-significant threshold — collapses 10⁶ and 10⁶+1",
                      polarity_at_finite_precision)])
    print(f"     [threshold control] p(10⁶) and p(10⁶+1) exactly: different "
          f"(the difference = {sp.simplify(polarity(sp.Integer(10**6)) - polarity(sp.Integer(10**6) + 1)) != 0}); "
          f"with an 8-significant threshold: {sp.N(polarity(sp.Integer(10**6)), 8)} = "
          f"{sp.N(polarity(sp.Integer(10**6) + 1), 8)}  [P3 · dimensionless · the same "
          "mapping under two resolutions]")

    def unbounded_in_mass(mapping):
        """A WORLD = A MAPPING ⟹ whether its image is UNBOUNDED (grows without a ceiling from m₀).
        The compact dual of the tact [T4-A] gives a ceiling ONLY because it lives in the phase."""
        vals = [sp.N(mapping(m)) for m in [sp.Integer(1), sp.Integer(10**3), sp.Integer(10**6)]]
        return vals[-1] > 10**5 * vals[0]

    def blade_ratio(m0v):
        """χ² = (ϸ+m₀)/(ϸ−m₀) — the same channel, in the dimensionless form of a ratio of weights."""
        ps = shell_ps(f_gen, m0v)
        return sp.simplify((ps + m0v) / (ps - m0v))

    ok(unbounded_in_mass, blade_ratio,
       "★★THE CHANNEL THROUGH WHICH THE TACT READS m₀ IS NONCOMPACT: χ² grows without a ceiling (1 → 10¹²) ⟹ "
       "the compactness of the dual [T4-A] does NOT GIVE a ceiling on m₀  [P3 · dimensionless · χ²(m₀) at "
       "m₀=1,10³,10⁶]",
       must_fail_on=[("the phase per tact itself (compact: |arg λ| ≤ π ∀m₀)",
                      lambda m: sp.Abs(sp.arg(sp.simplify(
                          eig_of(world_on_shell(f_gen, m))[0]))) + sp.Rational(1, 1000))])
    print("     ⟹ ★★A CONSEQUENCE FOR THE CARVED CONDITIONAL CONCLUSION (I report it JUST AS LOUDLY")
    print("       as a positive would be reported): the exante carved «if P3 gives a CEILING on the phase per")
    print("       tact ⟹ a ceiling on m₀ ⟹ via 2·kappa_stab=J(m₀), kappa_stab=Λ — a LOWER bound")
    print("       on Λ ⟹ a TWO-SIDED WINDOW».  THE PREMISE DOES NOT HOLD: m₀ is not a phase quantity,")
    print("       and no change of variable saves this (the polarity p=m₀/ϸ ∈ [0,1) and χ²")
    print("       are related by a BIJECTION ⟹ the injectivity is the same).  ⟹ Λ_min does NOT FOLLOW from this.")
    print("     ⟹ ★AND WHAT REALLY HAPPENS AT kappa_stab=0 (m₀→∞): the contrast CONVERGES")
    print("       to zero as 1−p ≈ |f|²/(2m₀²), but this is CONVERGENCE, not COLLISION.  Calling")
    print("       it «indistinguishability» is possible only against a resolution THRESHOLD, and there is")
    print("       no threshold in the model, and introducing one = a NEW CONSTANT ⟹ a K-handle ⟹ STOP.")
    print("     ★AND THIS IS NOT RHETORIC — IT IS SHOWN CONSTRUCTIVELY (the tooth above): the same")
    print("       mapping with an 8-significant threshold ALREADY collapses m₀=10⁶ and m₀=10⁶+1.")
    print("       A threshold — and only a threshold — produces indistinguishability.")
    for mv in [sp.Integer(0), sp.Integer(1), sp.Integer(10), sp.Integer(1000)]:
        print(f"       m₀={mv}: p = {sp.N(polarity(mv), 8)}  ·  χ² = {sp.N(blade_ratio(mv), 8)}"
              f"  [P3 · dimensionless · the polarity and the ratio of the blade weights]")
    print()

    # ═══════ ★WHAT T38 DOES GIVE: THE SINGLE CHANNEL THROUGH WHICH THE TACT READS THE SPLITTER ═══════
    print("★A SUMMARY MEASUREMENT: in HOW MANY channels the tact reads m₀")

    def mass_enters_only_polarity(world_pair):
        """A WORLD = A PAIR (m₀=0, m₀>0) on the shell ⟹ whether the DIFFERENCE between them sits EXACTLY in
        the polarity of the blades: the phases are the same, |λ| are the same (=1), and U†U differ."""
        w0, w1 = world_pair
        U0, U1 = U_tact(w0), U_tact(w1)
        same_phase = ([sp.nsimplify(sp.arg(v)) for v in sorted(sp.Matrix(U0).eigenvals(), key=str)]
                      == [sp.nsimplify(sp.arg(v)) for v in sorted(sp.Matrix(U1).eigenvals(), key=str)])
        same_modulus = sp.simplify(sp.Abs(eig_of(w0)[0]) - sp.Abs(eig_of(w1)[0])) == 0
        diff_polarity = sp.simplify(U0.H * U0 - U1.H * U1) != sp.zeros(2, 2)
        return same_phase and same_modulus and diff_polarity

    ok(mass_enters_only_polarity,
       (world_on_shell(f_gen, sp.Integer(0)), world_on_shell(f_gen, sp.Integer(5))),
       "★★★m₀ ENTERS THE TACT THROUGH EXACTLY ONE CHANNEL — THE POLARITY OF THE BLADES: the phases are the same, "
       "the moduli are the same (ρ=1), the difference is EXACTLY in U†U  [summary · a tact · a comparison of "
       "three invariants of U at m₀=0 and m₀=5]",
       must_fail_on=[("a pair off the shell (ϸ=2|f|) — there already the MODULI differ, not only the polarity",
                      ((f_gen, sp.Integer(0), 2 * absf), (f_gen, sp.Integer(5), 2 * absf)))])
    print("     ⟹ ★THE OPERATIONAL CONTENT OF THE AUTHOR'S WORD («a stable ONE-VECTOR splitting»):")
    print("       a split = the SOLE source of UNEQUAL weight of the two blades per tact (p = m₀/ϸ > 0).")
    print("       Without a split the blades weigh EQUALLY (p=0, U unitary) — a splitting exists, but")
    print("       it is SYMMETRIC, not one-vector.  ⟹ T38 is needed NOT for the existence of the tact")
    print("       (T4-C) and NOT for its boundedness (P2 above — empty), but for its")
    print("       POLARIZATION.  This is the THIRD property, and it is separate from existence and forcing.")
    print()

    print("★A POSITIVE CONTROL — THE HARNESS CAUGHT ME HERE TOO, AND THE CATCH TURNED OUT SUBSTANTIVE")
    print("  The first run: 9 ✓ / 1 ✗.  It was exactly the injectivity assert that failed — and NOT because")
    print("  the object is not injective, but because my detector compared values via")
    print("  `nsimplify`, which COLLAPSED p(10⁶) to exactly 1.  That is, the tool carried its own")
    print("  RESOLUTION THRESHOLD and would have reported «indistinguishability» as a measurement of the object —")
    print("  ★and this is EXACTLY the positive world the exante was waiting for.  A false POSITIVE")
    print("  nearly slipped through disguised as a confirmation of the wager.  Rewritten to an EXACT")
    print("  pairwise comparison, and the artifact itself was turned into a TOOTH (negative world 2).")
    print()
    code = report("S1061 T38-productivity")
    print()
    print("=" * 80)
    print("★THE BINARY OUTPUTS OF S1061")
    print("  FIELD 0. U = H·G₀ [S1045]; NOT unitary by construction; unitary EXACTLY on the")
    print("     splitter-free shell ⟹ P2 was a genuine question, not an exercise.")
    print("  P1. A trace of the tact EXISTS ⟺ there is a hop (f≠0), WITHOUT a split too.  ★At the node a split does NOT")
    print("     restore the trace (unlike the forcing of the slot in D) — ancestor multiplicity 1.")
    print("  P1-addendum. ★The alphabet of phases per tact is DISCRETIZED {0,±π/2,π}: a tact = a bit-register, not")
    print("     a dial.  The same fact as T4-B (no i∂ₜ) ⟹ multiplicity 1.")
    print("  P2. ★EMPTY, reported as a negative: U²=I on the shell ∀m₀ ⟹ divergence is impossible even without")
    print("     a split.  The carved risk of an inversion was realized MORE WEAKLY: a split breaks")
    print("     normality (‖U‖=χ>1), but this is a distortion without accumulation.")
    print("  P3. ★★A NEGATIVE BY TYPE: m₀ does not live in the phase ⟹ there is no winding ⟹ no ceiling ⟹")
    print("     ★THE EXANTE'S CONDITIONAL CONSEQUENCE (Λ_min, a two-sided window) DOES NOT FOLLOW.  The channel by which")
    print("     the tact reads m₀ is NONCOMPACT (χ² without a ceiling); injectivity holds ∀m₀.")
    print("     «Indistinguishability» would require a resolution THRESHOLD = a new constant = a K-handle.")
    print("  ★A POSITIVE. m₀ enters the tact through EXACTLY one channel — the polarity of the blades p=m₀/ϸ.")
    print("     ⟹ T38 carries a THIRD property: the POLARIZATION of the tact (not existence, not forcing).")
    print("     The author's wager gains operational content precisely here — and only here.")
    print("=" * 80)
    sys.stdout = tee.real
    logf.close()
    return code


if __name__ == "__main__":
    sys.exit(main())
