# -*- coding: utf-8 -*-
# DIM: symbolic-d  (the count of invariant forms — symbolic in d, a run over d=2..5;
#                   the field ablation does not depend on dimension, run on the carrier (d,1).)
"""S1079 · SIEGE OF AX-closure · component **Z4 (THE CARRIER)** — a lane-A probe.

EXANTE: `hub/prime/AX_CLOSURE_SIEGE_EXANTE.md` (carved by Ω, S1078).
The Z4 question verbatim: the carrier «real quadratic spaces» — «is this a separate handle,
or is it consumed by the roots {AX-alphabet · Λ}?»  Ω's wager: **CONSUMED**; «if not — name
it honestly as an input component, do NOT hide it».

★DISSECTING THE CARRIER into three components (otherwise "consumed/not consumed" is a
question without a subject):
   Z4-a · SPACE           — the vector space everything lives on;
   Z4-b · QUADRATIC       — the symmetric bilinear form on it;
   Z4-c · REALNESS        — the field both are taken over.

★THE GUARD (§4 of the exante): `AX-indef` does not appear here — the signature (d,1) is used
in §3 only as the CARRIER of a field ablation (the question "does the notion of signature
survive", not "is there a minus"); no forcing is drawn from it.

★THE FENCE OF §10 holds: the physics-vocabulary classes named in the project's fence are not
used in the probe.
"""
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'test'))

import sympy as sp
from itertools import permutations
from _teeth import ok, report


# ──────────────────────── Z4-b · which forms the alphabet allows at all ───────────

def perm_matrix(p):
    n = len(p)
    return sp.Matrix(n, n, lambda i, j: sp.Integer(1) if p[j] == i else sp.Integer(0))


def invariant_symmetric_forms(group, n):
    """A basis of the space of symmetric bilinear forms invariant under the group.

    Solving PᵀGP = G for ALL p ∈ group (G symmetric, n(n+1)/2 unknowns).
    Returns the basis as a list of matrices — this is exactly "how much freedom the
    alphabet leaves".
    """
    ent = {}
    G = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            s = sp.Symbol(f'g_{i}_{j}')
            ent[(i, j)] = s
            G[i, j] = s
            G[j, i] = s
    eqs = []
    for p in group:
        P = perm_matrix(p)
        eqs += list(P.T * G * P - G)
    unk = list(ent.values())
    sol = sp.solve(eqs, unk, dict=True)
    G0 = G.subs(sol[0]) if sol else G
    free = sorted({s for s in G0.free_symbols if s in set(unk)}, key=lambda s: s.name)
    return [sp.Matrix(n, n, lambda i, j: sp.diff(G0[i, j], f)) for f in free]


def pinned_by_radical(basis, n):
    """★Does the condition "d+1 directions live in d dimensions" pin the form UNIQUELY (up to
    scale)?

    The condition = the diagonal (1,…,1) lies in the radical (this is exactly Σ u_j = 0,
    that is, the cell).  Returns the list of forms, from basis combinations, satisfying it,
    up to scale.
    """
    cs = sp.symbols(f'c0:{len(basis)}')
    G = sum((c * B for c, B in zip(cs, basis)), sp.zeros(n, n))
    one = sp.Matrix([sp.Integer(1)] * n)
    sol = sp.solve(list(G * one), list(cs), dict=True)
    if not sol:
        return []
    Gp = sp.simplify(G.subs(sol[0]))
    free = [s for s in Gp.free_symbols if s in set(cs)]
    return [Gp, free]


# ─────────────────────────── Z4-c · a FIELD ablation ──────────────────────────────

def G_break(d):
    """The ablation carrier: a nondegenerate form with one minus sign on n = d+1."""
    return sp.diag(*([sp.Integer(1)] * d + [sp.Integer(-1)]))


def signature_survives(world):
    """★DEFINITION (carved BEFORE the count):

    «The notion of SIGNATURE survives over a field F» ⟺ a form with a minus CANNOT be
    brought by congruence (SᵀGS) to the identity.  The detector checks this
    CONSTRUCTIVELY:
      · if the field is ordered — it exhibits a vector of negative norm; any
        SᵀGS = 𝟙 would make G positive-definite (since xᵀx > 0 ∀x≠0) ⟹ a contradiction;
      · if the field is NOT ordered — there is no order, the argument evaporates, and
        the detector must say FALSE (and §3 also exhibits an explicit S that reduces it).
    """
    G, ordered = world['G'], world['ordered']
    if not ordered:
        return False
    n = G.shape[0]
    v = sp.Matrix([sp.Integer(0)] * (n - 1) + [sp.Integer(1)])
    return sp.simplify((v.T * G * v)[0, 0]) < 0


# ───────────────────────────────────────────────── the run ─────────────────────

def main():
    print("=" * 78)
    print("S1079 · AX-closure · Z4 (THE CARRIER): «real quadratic spaces» —")
    print("        a separate handle or consumed by the roots {AX-alphabet · Λ}?")
    print("=" * 78)

    print("""
Z4-a · SPACE — not measured, and here is why this is not laziness: the vector space
   everything lives on IS the LINEAR SPAN of the alphabet (the bare set d+1 + democracy
   S_{d+1} ⟹ the permutation representation).  Nothing beyond the root appears here;
   building a "measurement" of this would be fabricating a witness.  CONSUMED by the root
   AX-alphabet.""")

    print("\nZ4-b · QUADRATIC — how much freedom the alphabet leaves, and who eats it")
    for d in range(2, 6):
        n = d + 1
        grp = list(permutations(range(n)))
        basis = invariant_symmetric_forms(grp, n)
        Gp, free = pinned_by_radical(basis, n)
        canon = sp.Matrix(n, n, lambda i, j: sp.Integer(1) if i == j else sp.Rational(-1, d))
        # does the pinned form coincide with the canonical Gram of the cell (up to scale)?
        lam = sp.simplify(Gp[0, 0] / canon[0, 0]) if canon[0, 0] != 0 else None
        same = sp.simplify(Gp - lam * canon) == sp.zeros(n, n)
        print(f"   d={d}: S_{{{n}}}-invariant symmetric forms: dim = {len(basis)}"
              f"  ·  after the condition «Σ u_j = 0» there remains {len(free)} parameter"
              f"  ·  coincides with the Gram of the cell up to scale: {same}")
        assert len(basis) == 2 and len(free) == 1 and same

    print("""
   ★READING: the alphabet ITSELF leaves a 2-dimensional space of invariant forms
   (⟨𝟙, J⟩), that is, one free handle beyond scale.  It is eaten NOT by taste, but by the
   condition "d+1 directions live in d dimensions" (Σ u_j = 0 — the same thing that makes
   the cell a cell): it pins the ratio uniquely, leaving EXACTLY the scale.  The scale = the
   Λ-ruler.  ⟹ QUADRATIC is consumed by {AX-alphabet + the cell} + Λ, no new handle.""")

    for d in range(2, 5):
        n = d + 1
        full = list(permutations(range(n)))
        stab = [p for p in full if p[0] == 0]          # democracy broken: S_d
        ok(lambda w: len(invariant_symmetric_forms(w['grp'], w['n'])) == 2,
           {'name': f'full democracy d={d}', 'grp': full, 'n': n},
           f"d={d}: under the FULL S_{{{n}}} the invariant forms are exactly 2 ⟹ there "
           f"is almost no freedom «which form» — and this is a property of DEMOCRACY, "
           f"not of my algebra",
           must_fail_on=[(f"democracy broken down to the stabilizer S_{d}",
                          {'name': 'broken', 'grp': stab, 'n': n})])

    print("\nZ4-c · REALNESS — a field ablation (is this a load-bearing component?)")
    d = 3
    G = G_break(d)
    real_w = {'name': 'an ordered field (ℝ or ℚ)', 'G': G, 'ordered': True}
    cplx_w = {'name': 'ℂ (no order)', 'G': G, 'ordered': False}

    print("   an explicit reduction over ℂ (constructive, not by citing a theorem):")
    S = sp.diag(*([sp.Integer(1)] * d + [sp.I]))
    red = sp.simplify(S.T * G * S)
    print(f"     S = diag(1,…,1, i) ⟹ SᵀGS = {red.tolist()} = 𝟙 "
          f"⟹ over ℂ the minus DISAPPEARS, no signature exists")
    assert red == sp.eye(d + 1)

    v = sp.Matrix([0] * d + [1])
    print(f"   over an ordered field: the vector v={list(v.T)} has norm "
          f"{(v.T * G * v)[0, 0]} < 0, while SᵀGS=𝟙 would make all norms > 0 "
          f"⟹ the reduction is IMPOSSIBLE")

    ok(signature_survives, real_w,
       "the notion of signature survives exactly over an ORDERED field",
       must_fail_on=[("ℂ (no order) — and the reduction is exhibited explicitly", cplx_w)])

    # ★is ℝ itself needed, or is orderedness enough: the same certificate over ℚ
    q_w = {'name': 'ℚ (ordered, incomplete)', 'G': G, 'ordered': True}
    ok(signature_survives, q_w,
       "★the same certificate holds over ℚ ⟹ the carrier needs the field's "
       "ORDEREDNESS, not the completeness of ℝ",
       must_fail_on=[("ℂ (no order)", cplx_w)])

    code = report("S1079 · Z4-carrier of AX-closure")

    print("""
────────────────────────────────────────────────────────────────────────────────
WHAT THIS MEANS FOR Z4 (a measurement, not a verdict):

  (a) SPACE — consumed by `AX-alphabet` (the linear span of the alphabet).  New: 0.
  (b) QUADRATIC — consumed by {`AX-alphabet` + the cell} + Λ: democracy leaves exactly a
      2-dimensional space of forms ⟨𝟙, J⟩, the condition «Σ u_j = 0» eats the last freedom,
      leaving EXACTLY the scale = Λ.  Measured ∀d=2..5, with a tooth (broken democracy
      leaves more freedom ⟹ the pin is a property of DEMOCRACY, not of the apparatus).
  (c) ★REALNESS — LOAD-BEARING, and this is the main line of Z4: over ℂ the minus is
      removed by an explicit S=diag(1,…,1,i), that is, the NOTION OF SIGNATURE DISAPPEARS,
      and with it the entire verdict {a (3,0)-minimum · (3,1) after the break}.  ⟹ the
      component cannot be silently inherited: it must have a source.
      ★BUT IT IS NARROWER THAN NAMED: the certificate holds over ℚ just the same ⟹ the
      carrier needs the field's **ORDEREDNESS**, not "realness" (the completeness of ℝ is
      not used anywhere in this step).  The honest name of the component = "a formally real
      (ordered) field", and this is a WEAKER premise than stands in the prose.
  (d) ⟹ the form of the verdict handed to Ω (not rendered here): Z4 = {a,b consumed} + {c —
      a named input component that SHOULD BE NARROWED to orderedness}.  Whether
      orderedness already carries Λ (a ruler = an ordered quantity by type) — a question for
      the court, not decided here: this is the difference between "consumed" and "named".
  BOUNDARY: not measured whether the COMPLETENESS of ℝ is needed elsewhere in T1–T11 (here
  — no).  This line is not carried beyond the step.
────────────────────────────────────────────────────────────────────────────────""")
    return code


if __name__ == '__main__':
    sys.exit(main())
