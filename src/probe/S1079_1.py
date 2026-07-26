# -*- coding: utf-8 -*-
# DIM: symbolic-d  (an order on signatures — combinatorial, a scan p,q ≤ N without instances;
#                   no d=3 numbers on the load-bearing path.)
"""S1079 · SIEGE OF AX-closure · component **Z1 (DEFINITION)** — a lane-A probe.

EXANTE: `hub/prime/AX_CLOSURE_SIEGE_EXANTE.md` (carved by Ω, S1078).
The Z1 question verbatim: «does **(a)** give the SAME verdict {a (3,0)-minimum · (3,1) after
the break}?»  Ω's wager: INVARIANCE; the named candidate for the wager's failure — «(a),
without coercion, may NOT have a notion of a forced successor ⟹ then (a) is not an
alternative, but a WEAKER order — and this too is a verdict about rigidity (the alternative
is eliminated BY TYPE, not by taste)».

★DEFINITIONS pulled up VERBATIM (step 0, `W28_ORDER_PRECEDENCE_EXANTE.md` §2):
  (a) EMBEDDING:  (p₁,q₁) ≼ (p₂,q₂) ⟺ p₁ ≤ p₂ ∧ q₁ ≤ q₂   [a closed form, Sylvester/Witt]
  (b) COMPLETION: C1 ≺ C2 ⟺ Lie(V₁) ⊋ V₁ and C2 = (the minimal carrier, Lie(V₁))

★WHAT IS MEASURED (and what is deliberately NOT measured):
  §1 · does (a) reproduce the NUMBER of the verdict — the minimal elements under the stated
       premises;
  §2 · does (a) have a notion of a FORCED SUCCESSOR — a count of the minimal strict upper
       bounds (in (b) the successor is the VALUE of the closure function, that is, one);
  §3 · A RE-MEASUREMENT CHECK: P1 («in (a) there are no maximal elements») was carved by Ω
       back on 2026-07-16 — NOT re-measured, quoted and its multiplicity handed to the
       ancestor.

★THE GUARD OF §4: `AX-indef` is not involved — the premise `q≥1` enters here as a SCAN
PARAMETER (I scan q₀ = 0 and 1 on equal footing), not as a forcing source.
★THE FENCE OF §10: frozen words are not used.
"""
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'test'))

from _teeth import ok, report

N = 8          # the ceiling of the ladder scan (p,q ≤ N) — the verdict does not depend on it, see the §1 control


# ─────────────────────────────── (a) as an object ────────────────────────────────

def leq_a(x, y):
    """(a) EMBEDDING, closed form: (p₁,q₁) ≼ (p₂,q₂) ⟺ p₁≤p₂ ∧ q₁≤q₂."""
    return x[0] <= y[0] and x[1] <= y[1]


def ladder(n=N):
    return [(p, q) for p in range(n + 1) for q in range(n + 1)]


def minimal_elements(cond, n=N):
    """The minimal elements of the set {(p,q) : cond} in the order (a)."""
    S = [x for x in ladder(n) if cond(x)]
    return sorted(x for x in S if not any(y != x and leq_a(y, x) for y in S))


def minimal_strict_upper_bounds(x, cond=lambda z: True, n=N):
    """★The minimal STRICT upper bounds of an element — the "successors" in (a).

    How many there are is exactly the question about forcing: one ⟹ the successor is
    coerced; more than one ⟹ the notion of a «forced successor» is not defined in this
    order.
    """
    up = [y for y in ladder(n) if cond(y) and leq_a(x, y) and y != x]
    return sorted(y for y in up if not any(z != y and leq_a(z, y) and z != x
                                           and leq_a(x, z) for z in up))


# ─────────────────────────── detectors (functions of a WORLD) ────────────────────

def verdict_is_31(world):
    """Does the world's order give the minimum {(3,1)} under the world's stated premises?"""
    p0, q0 = world['p_min'], world['q_min']
    return minimal_elements(lambda x: x[0] >= p0 and x[1] >= q0) == [(3, 1)]


def has_forced_successor(world):
    """Is the notion of a FORCED successor (exactly one) defined in the order?"""
    return all(len(minimal_strict_upper_bounds(x, world['cond'])) == 1
               for x in world['probe_points'])


# ───────────────────────────────────────────────── the run ─────────────────────

def main():
    print("=" * 78)
    print("S1079 · AX-closure · Z1 (DEFINITION): does (a) give the same verdict as (b)?")
    print("=" * 78)

    # ── §1 · the NUMBER of the verdict under (a) ────────────────────────────────────────────
    print("\n§1 · the NUMBER of the verdict in the order (a) — minimal elements under the premises")
    cases = [("{p≥3, q=0} — a definite carrier", 3, 0, [(3, 0)]),
             ("{p≥3, q≥1} — after the break", 3, 1, [(3, 1)])]
    for name, p0, q0, want in cases:
        got = minimal_elements(lambda x: x[0] >= p0 and x[1] >= q0)
        print(f"   {name}: minimal = {got}  {'✓' if got == want else '✗'}")
        assert got == want

    print("   a control on the scan ceiling (the verdict must not depend on it):")
    for n in (5, 8, 12):
        got = minimal_elements(lambda x: x[0] >= 3 and x[1] >= 1, n)
        print(f"     N={n}: {got}")
        assert got == [(3, 1)]

    ok(verdict_is_31, {'name': 'the canonical premises {p≥3, q≥1}', 'p_min': 3, 'q_min': 1},
       "(a) reproduces the NUMBER of the verdict exactly: the minimum under {p≥3, q≥1} = (3,1)",
       must_fail_on=[("a weakened premise {p≥2, q≥1}",
                      {'name': 'p≥2', 'p_min': 2, 'q_min': 1}),
                     ("a weakened premise {p≥3, q≥0}",
                      {'name': 'q≥0', 'p_min': 3, 'q_min': 0})])

    print("""   ★READING §1: the verdict-NUMBER in (a) is the same — and the tooth shows WHY: it
   is a function of the PREMISES {p≥3 · q≥1}, not of the order's machinery.  Weaken a
   premise — the number moves; replace the order's definition — it does not.  This is a
   RE-MEASUREMENT of the S908-F.4 line («the entire content of T5 lives in two premises»),
   a different road to the same place; the multiplicity is handed to the ancestor, no line
   of my own is minted here.""")

    # ── §2 · FORCING: does (a) have a notion of a coerced successor ────────────────
    print("\n§2 · ★DOES (a) HAVE A FORCED SUCCESSOR (a discriminating measurement)")
    pts = [(3, 1), (3, 0), (2, 2), (4, 1)]
    for x in pts:
        s = minimal_strict_upper_bounds(x)
        print(f"   (a): the minimal strict upper bounds of {x} = {s}  ⟹ {len(s)} of them.")
    print("   (b): the successor = the VALUE of the closure function Lie(V) ⟹ exactly 1 by construction")

    chain = {'name': 'a chain (a total order q=0)', 'cond': lambda z: z[1] == 0,
             'probe_points': [(3, 0), (4, 0), (5, 0)]}
    plane = {'name': '(a) on the full ladder', 'cond': lambda z: True,
             'probe_points': pts}
    ok(has_forced_successor, chain,
       "control: in a TOTAL order the successor really is one — the detector can "
       "say TRUE, so its «FALSE» on (a) is a measurement, not blindness",
       must_fail_on=[("(a) on the full ladder", plane)])

    print("""   ★READING §2: in (a) every element has TWO minimal strict upper
   neighbors ((p+1,q) and (p,q+1)) ⟹ the notion of a «forced successor» is NOT
   DEFINED in (a).  In (b) the successor is the value of the closure function, that
   is, one.
   ⟹ (a) is not an alternative definition of THE SAME object: it is STRUCTURALLY
   weaker.  This is exactly the branch that Ω named BEFORE the probes («eliminated
   by TYPE»).""")

    # ── §3 · a re-measurement check ──────────────────────────────────────────────────
    print("""
§3 · A RE-MEASUREMENT CHECK (not spending the count on what is already carved)
   P1 (Omega, W28 §4, carved on 2026-07-16 BEFORE the numbers, verbatim):
     «in pure (a) at level-1 there are NO maximal elements (the ℕ²-lattice is infinite
      upward) ⟹ (a) itself cannot be a selector — its role = a skeleton of the order
      for Z2.  Carved so that this is not later "discovered" as a result.»
   ⟹ NOT re-measured.  The multiplicity goes to the ancestor; my §2 measures something
   DIFFERENT (not "is there a maximum", but "is a successor defined") — and precisely
   for that reason it is not a duplicate of P1.""")

    code = report("S1079 · Z1-definition of AX-closure")

    print("""
────────────────────────────────────────────────────────────────────────────────
WHAT THIS MEANS FOR Z1 (a measurement, not a verdict):

  (1) ★THE VERDICT-NUMBER IS INVARIANT: (a) hands back {(3,0)} and {(3,1)} under the same
      premises — exactly, and independent of the scan ceiling.  Ω's wager on the number
      holds.

  (2) ★BUT THE INVARIANCE IS NOT WHAT IT SEEMS: the number agrees not because the two
      definitions "agreed", but because it NEVER USED their difference.  The verdict is a
      function of the PREMISES {p≥3 · q≥1} + minimality; the machinery of (b) (closure,
      forcing) does not enter the NUMBER.  The tooth of §1 shows this by hand: the number
      moves when a premise is weakened and does NOT move when the definition is replaced.

  (3) ★WHAT EXACTLY IS LOST BY "TWISTING IT" — and this answers the author's fear
      («someone will twist it»): on the transition (b)→(a) the NUMBER DOES NOT MOVE, but
      the FORCING DISAPPEARS.  In (a) every element has TWO minimal strict successors
      ⟹ the notion of a "coerced successor" is not defined; in (b) it is the value of a
      function.  ⟹ (a) is not an alternative, but a STRUCTURALLY WEAKER order:
      the alternative is eliminated BY TYPE, not by taste — exactly the branch Ω carved
      BEFORE the probes.

  (4) ⟹ the form of the verdict handed to Ω (not rendered here): on Z1 — RIGIDITY, but with
      a precise qualifier: a rigid NUMBER (it does not move under a legal twist) ⊥ a
      NON-rigid PROVENANCE (the strength of the claim moves: without (b) the verdict
      becomes a consequence of two premises, not a partly-forced one).  To the reader of a
      preprint who "twists the order", we owe both lines.

  A BOUNDARY: measured at level-1 (classes (p,q)), because that is exactly where (a) has a
  closed form.  Level-2 ((carrier, V)) for (a) was not measured — and no conclusion is drawn
  there.
────────────────────────────────────────────────────────────────────────────────""")
    return code


if __name__ == '__main__':
    sys.exit(main())
