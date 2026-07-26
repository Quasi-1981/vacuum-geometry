# -*- coding: utf-8 -*-
# DIM: symbolic-d  (forms (d,1)/(d+1,0)/degenerate — symbolic in d, a run over d=2..5;
#                   no d=3 instance on the load-bearing path.)
"""S1075 · SIEGE OF AX-dimer · component **D1 (EXISTENCE)** — a lane-A probe.

EXANTE: `hub/prime/AX_DIMER_SIEGE_EXANTE.md` (carved by Ω, S1074).
The D1 question verbatim: «there is ≥1 marked axis (q≥1) — does even one native law FORCE it?»
Ω's wager: **IT IS NOT FORCED** ⟹ an irreducibility theorem ⟹ AX-dimer collapses to a single
existential bit, «time exists».

★THE CYCLE GUARD (§3 of the exante) — and exactly how it is upheld here:
  the legal sources for D1 = {AX-lambda · AX-alphabet · AX-closure · the prime directive};
  `AX-indef` and its descendants are FORBIDDEN as a FORCING source.
  In §2 of this probe `AX-indef` IS USED — but **as a REINFORCEMENT OF THE NEGATIVE (a
  fortiori)**: we freely GIFT the siege the forbidden root and show that the axis is **still
  not forced**.  A cycle needs FORCING; a negative reinforced by a gift does not create a
  cycle (it makes the verdict stronger, not weaker).  This is stated EXPLICITLY — for Beta/Ω
  to judge.

★A HOMONYM this probe is obliged to resolve (otherwise the verdict is spurious):
  «the axis is FORCED» has TWO different definitions, and they give different answers —
   · **RELATIVE** (the language of S1066-B2): does an isometry of the form carry the
     candidate axis `u₀` into another axis of the **ALPHABET** (out of d+1 named ones)? — a
     question about a LABEL inside the alphabet;
   · **ABSOLUTE** (the language of this probe): does there exist a line, invariant under
     **THE ENTIRE** group of isometries of the form? — a question of whether the form can
     single out an axis at all.
  Conflating them buys a spurious forcing.  The probe measures BOTH and prints the contrast.

WHAT IS MEASURED:
  §1  A RE-MEASUREMENT CHECK: what of D1 the ancestors have already measured (do not spend
      the count twice).
  §2  THE CORE (new): in (d,1) the form forces a CONE, not an AXIS — even with AX-indef
      gifted.
  §3  the contrast of the two definitions of «forced» + the degenerate world as a tooth.
"""
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'test'))

import sympy as sp
from _teeth import ok, report


# ───────────────────────────────── worlds = quadratic forms ──────────────────

def W_definite(n):
    """A definite arena (n,0) — the world q_sig=0: «the legal ground-state neighbor» W42."""
    return {'name': f'definite ({n},0)', 'G': sp.eye(n), 'n': n}


def W_lorentz(d):
    """(d,1) — an arena with a break.  Here AX-indef IS GIFTED (a fortiori, see the header)."""
    n = d + 1
    return {'name': f'indefinite ({d},1) [AX-indef gifted]',
            'G': sp.diag(*([sp.Integer(1)] * d + [sp.Integer(-1)])), 'n': n}


def W_degenerate(d):
    """★A NEGATIVE WORLD: a degenerate form with a 1-dimensional radical.

    The radical is invariant under ALL isometries ⟹ the axis is genuinely FORCED by the
    form.  (The same third world that S1066-B2 found outside its own exante — taken here
    deliberately as a control: the detector must say «forced» HERE.)
    """
    n = d + 1
    return {'name': f'DEGENERATE (radical dim 1), n={n}',
            'G': sp.diag(*([sp.Integer(1)] * d + [sp.Integer(0)])), 'n': n}


# ───────────────────────── detector: does the form force an AXIS (absolute) ─────────

def iso_algebra(G, n):
    """A basis of the Lie algebra of isometries of the form: {X : XᵀG + GX = 0}."""
    ent = sp.symbols(f'x0:{n * n}')
    X = sp.Matrix(n, n, ent)
    eqs = list(X.T * G + G * X)
    sol = sp.solve(eqs, ent, dict=True)
    assert sol, "empty solution of the isometry system"
    X0 = X.subs(sol[0])
    free = sorted({s for s in X0.free_symbols if s in set(ent)}, key=lambda s: s.name)
    return [sp.Matrix(n, n, lambda i, j: sp.diff(X0[i, j], f)) for f in free]


def invariant_lines(world):
    """All lines span(v) invariant under THE WHOLE algebra of isometries — a direct search.

    ★WHY NOT THE COMMUTANT (recorded as a caught defect, not as a style choice): the first
    draft of this detector counted the dimension of the commutant and read «dim=1 ⟹
    irreducible ⟹ no line».  This is FALSE for a NON-reductive algebra: in the degenerate
    world the isometry algebra is not semisimple, the representation is INDECOMPOSABLE but
    REDUCIBLE — the commutant equals 1, yet an invariant line (the radical) EXISTS.  The
    defect was caught by `_teeth` (8 ☠ EMPTIES: the detector said YES even on the negative
    world).  Hence here — a definition head-on, without representation theory.

    v defines an invariant line ⟺ Xv ∥ v for EVERY generator X
    ⟺ rank[Xv | v] ≤ 1 ∀X ⟺ all 2×2 minors vanish.
    Solved chart-by-chart with v_i = 1 (the line has some nonzero coordinate).
    """
    G, n = world['G'], world['n']
    gens = iso_algebra(G, n)
    found = []
    for i in range(n):
        ts = sp.symbols(f't0:{n}')
        v = sp.Matrix([sp.Integer(1) if k == i else ts[k] for k in range(n)])
        eqs = set()
        for X in gens:
            Xv = X * v
            for a in range(n):
                for b in range(a + 1, n):
                    e = sp.expand(Xv[a] * v[b] - Xv[b] * v[a])
                    if e != 0:
                        eqs.add(e)
        unknowns = [ts[k] for k in range(n) if k != i]
        sols = sp.solve(list(eqs), unknowns, dict=True) if eqs else [{}]
        for s in sols:
            vv = sp.simplify(v.subs(s))
            if any(x.free_symbols for x in vv):          # a whole parametric pencil
                found.append(('pencil', vv))
            elif not any(sp.simplify(vv - w[1]) == sp.zeros(n, 1) for w in found):
                found.append(('line', vv))
    return found


def axis_not_forced(world):
    """★DEFINITION (carved BEFORE the count, absolute):

    «The axis is FORCED by the form» ⟺ there exists a line, invariant under THE ENTIRE group
    of isometries.  The detector returns TRUE when no such line EXISTS (that is, the axis is
    NOT forced).
    """
    return len(invariant_lines(world)) == 0


# ───────────────────────── the relative definition (the language of S1066-B2) ─────────────────

def alphabet_axes(d):
    """d+1 candidate axes of the alphabet: cell weights (Gram 1 diag / −1/d)."""
    n = d + 1
    return sp.Matrix(n, n, lambda i, j: sp.Integer(1) if i == j else sp.Rational(-1, d))


def mark_is_label_only(gram):
    """RELATIVE: does a permutation of the alphabet's axes exist that preserves the Gram?

    TRUE ⟹ the mark on the axis is a LABEL (no invariant carrier) ⟹ the axis is not singled
    out.  We check the transposition (0↔1): it preserves the Gram iff the norms and cross
    products are symmetric with respect to it.
    """
    n = gram.shape[0]
    P = sp.eye(n)
    P[0, 0] = P[1, 1] = 0
    P[0, 1] = P[1, 0] = 1
    return sp.simplify(P.T * gram * P - gram) == sp.zeros(n, n)


# ───────────────────────────────────────────────── the run ─────────────────────

def main():
    print("=" * 78)
    print("S1075 · AX-dimer · D1 (EXISTENCE): does even one native law force the mark?")
    print("=" * 78)

    print("""
§1 · A RE-MEASUREMENT CHECK (the lesson of B1/B3: first look for a re-measurement, only then
     spend the count)
────────────────────────────────────────────────────────────────────────────────
  · AX-closure — THE NEGATIVE IS ALREADY CARVED IN THE CANON, not re-measured.  CODEX §0
      verbatim: «the existence of a time-axis = a postulated break (NOT from closure: the
       minimal carrier of closure = definite (3,0))».   ⟹ multiplicity by ancestor = 1, the
      witness = the ancestor.
  · the prime directive — not a law about an OBJECT, but a meta-rule for counting handles;
      moreover it presses IN THE OPPOSITE DIRECTION (q=0 is cheaper than q≥1).  It cannot
      force existence, by its very type.  No probe is built here: building a "measurement"
      of a meta-rule would be fabricating a witness.
  · AX-lambda — the scale ruler.  A homothety acts identically on ALL d+1 directions
      (it commutes with S_{d+1}) ⟹ by construction it cannot single out an axis.  The check
      below enters §3 as part of the relative measurement (the Gram is invariant under
      scaling).
  · AX-alphabet — the ONLY source where the question is not vacuous.  Measured in §3.
────────────────────────────────────────────────────────────────────────────────""")

    print("\n§2 · THE CORE: does the form itself force an AXIS (the absolute definition)?")
    print("     [AX-indef gifted for free — a fortiori, see the header]")
    for d in range(2, 5):
        lor, deg, defi = W_lorentz(d), W_degenerate(d), W_definite(d + 1)
        print(f"\n   d={d}:  invariant lines — ({d},1): {len(invariant_lines(lor))} · "
              f"({d + 1},0): {len(invariant_lines(defi))} · "
              f"degenerate: {[str(list(v.T)) for _, v in invariant_lines(deg)]}")
        ok(axis_not_forced, lor,
           f"d={d}: in ({d},1) there is NO line invariant under all isometries ⟹ "
           f"even the break does NOT force an axis",
           must_fail_on=[("a degenerate form (the radical = an invariant line)", deg)])
        ok(axis_not_forced, defi,
           f"d={d}: in ({d + 1},0) there is also NO invariant line ⟹ the definite "
           f"ground-state neighbor has no axis either",
           must_fail_on=[("a degenerate form (the radical = an invariant line)", deg)])

    print("\n   ★WHAT THE BREAK ACTUALLY IS, IF NOT AN AXIS: an explicit witness of transitivity (a boost)")
    for d in range(2, 6):
        n = d + 1
        G = W_lorentz(d)['G']
        w = sp.symbols('omega', real=True)
        R = sp.eye(n)
        R[0, 0] = sp.cosh(w); R[0, n - 1] = sp.sinh(w)
        R[n - 1, 0] = sp.sinh(w); R[n - 1, n - 1] = sp.cosh(w)
        assert sp.simplify(R.T * G * R - G) == sp.zeros(n, n), f"the boost is not an isometry, d={d}"
        u = sp.Matrix([0] * (n - 1) + [1])              # a timelike ray
        Ru = sp.simplify(R * u)
        para = sp.simplify(Ru[0] * u[n - 1] - Ru[n - 1] * u[0])   # a ∥-test
        norm_u = sp.simplify((u.T * G * u)[0, 0])
        norm_Ru = sp.simplify((Ru.T * G * Ru)[0, 0])
        assert norm_u == norm_Ru == -1
        print(f"     d={d}: RᵀGR=G ✓ · the norm is preserved ({norm_Ru}) ✓ · "
              f"R·u ∦ u for ω≠0 (the determinant test = {sp.simplify(para)}) ⟹ "
              f"the isometry MOVES the timelike ray")

    print("\n   ★A CONE, NOT A LINE: the timelike set has FULL dimension")
    t = sp.symbols('t', real=True)
    for d in range(2, 6):
        n, G = d + 1, W_lorentz(d)['G']
        dims = []
        for i in range(n - 1):
            v = sp.Matrix([0] * n); v[n - 1] = 1; v[i] = t
            dims.append(sp.simplify((v.T * G * v)[0, 0]))          # t² − 1 < 0 for |t|<1
        assert all(sp.simplify(x - (t ** 2 - 1)) == 0 for x in dims)
        print(f"     d={d}: perturbing a timelike vector in EVERY one of the {n - 1} "
              f"transverse directions leaves the norm {dims[0]} < 0 for |t|<1 ⟹ "
              f"the set of timelike vectors has dimension {n}, not 1")

    print("\n§3 · THE CONTRAST OF THE TWO DEFINITIONS OF «forced» (resolving the homonym)")
    print("     ★A DETECTOR FINDING (not suppressed — examined): in the alphabet's Gram an invariant")
    print("     line EXISTS.  But it is NOT an axis: the Gram is degenerate (Σ u_j = 0), and its radical —")
    print("     is the DIAGONAL (1,…,1), that is, democracy itself.  A mark must BREAK S_{d+1};")
    print("     S_{d+1} leaves the diagonal fixed ⟹ an invariant carrier exists, but it carries no mark.")
    for d in range(2, 5):
        n = d + 1
        gram = alphabet_axes(d)
        rel = mark_is_label_only(gram)
        lines = invariant_lines({'name': 'alphabet', 'G': gram, 'n': n})
        diag = sp.Matrix([sp.Integer(1)] * n)
        is_diag = [sp.simplify(v - diag) == sp.zeros(n, 1) for _, v in lines]
        # is even ONE coordinate axis of the alphabet among the invariant lines?
        axes_hit = [k for k in range(n)
                    for _, v in lines
                    if sp.simplify(v - sp.Matrix([sp.Integer(1) if i == k else sp.Integer(0)
                                                  for i in range(n)])) == sp.zeros(n, 1)]
        print(f"   d={d}: RELATIVE (a permutation of the axes preserves the Gram) = {rel} ⟹ the mark = a LABEL"
              f"   ·   invariant lines: {len(lines)}, among them the diagonal: {any(is_diag)},"
              f" coordinate axes of the alphabet: {len(axes_hit)}")
        assert rel and any(is_diag) and not axes_hit

    ok(lambda w: mark_is_label_only(w['G']), {'name': 'alphabet d=3', 'G': alphabet_axes(3)},
       "democracy S_{d+1}: a transposition of axes is an isometry of the Gram ⟹ the mark has "
       "no invariant carrier IN THE ALPHABET (this is a RE-MEASUREMENT of S1066-B2, row 1)",
       must_fail_on=[("a Gram with an unequal axis norm (norm 2)",
                      {'name': 'unequal norms', 'G': sp.diag(2, 1, 1, 1)})])

    code = report("S1075 · D1-existence of AX-dimer")

    print("""
────────────────────────────────────────────────────────────────────────────────
WHAT THIS MEANS FOR D1 (a measurement, not a verdict):

  (1) No LEGAL source forces the mark: closure — the negative is already in the canon (a
      re-measurement, not my own line); the prime directive presses the other way; Λ
      commutes with democracy; the alphabet S_{d+1} makes the mark a LABEL (a
      re-measurement of S1066-B2).
  (2) ★NEW: even the GIFTED AX-indef does not deliver an axis.  In (d,1) the group of
      isometries acts irreducibly (commutant = 1) ⟹ no invariant line exists; the timelike
      set has FULL dimension (a cone, not a line); a boost explicitly MOVES the timelike
      ray.
      ⟹ the break gives a CLASS (a cone); choosing a RAY within it is a separate bit.
  (2b) ★A SIDE FINDING (brought by the detector, not sought): the only invariant line of
      the alphabet is the DIAGONAL (1,…,1) = the radical of the degenerate Gram = the
      relation Σ u_j = 0.  It is S_{d+1}-FIXED ⟹ by construction it carries no mark (a mark
      must BREAK democracy).  That is, the alphabet has exactly one invariant carrier, and
      it IS democracy itself.  This reinforces (1): the alphabet does not simply "stay
      silent" — it has a single singled-out direction, and that direction is symmetric.

  (3) ⟹ D1 does not reduce even to the forbidden root.  So AX-dimer and AX-indef are NOT
      duplicates (important: if D1 ≡ AX-indef, the demotion of S1066 would be a cycle);
      the letter `q` in each denotes DIFFERENT quantities (marked axes ⊥ minus signs of
      the signature).
  (4) The form of the verdict handed to Ω (not rendered here): the existence of the mark =
      a SINGLE existential bit, not derivable from what is accepted — an irreducibility
      theorem, a rhyme to T39 (there the undeliverable was the SIGN, here it is the
      EXISTENCE).
────────────────────────────────────────────────────────────────────────────────""")
    return code


if __name__ == '__main__':
    sys.exit(main())
