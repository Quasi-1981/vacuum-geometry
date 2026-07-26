# -*- coding: utf-8 -*-
# DIM: symbolic-d  (the orbit count is exact ∀d as a statement; the run d=2..4 is exhaustive
#                   over (d+1)! permutations and 2^(d+1) subsets — no instances.)
"""S1075 · SIEGE OF AX-dimer · component **D3 (CHOICE)** — a lane-A probe.

EXANTE: `hub/prime/AX_DIMER_SIEGE_EXANTE.md` (carved by Ω, S1074).
The D3 question verbatim: «WHICH bond is marked (out of d+1 equal-rights ones) — is there an
irreducibility theorem for the choice?»  Ω's wager: **AN IRREDUCIBILITY THEOREM**, the
candidate mechanism — the democracy `S_{d+1}` (verify with a FIELD, do not just believe it).

★THE CYCLE GUARD: `AX-indef` and its descendants do not appear here at all — the whole probe
lives in the alphabet (`AX-alphabet`) and its Gram (`AX-cell`/T19).  There is no signature in
the code.

★WHAT IS MEASURED (and why this is not a retelling of an orbit triviality):
  §1  the input really is `S_{d+1}`-invariant (I check ALL (d+1)! permutations, not just
      believe it);
  §2  the action on the axes is TRANSITIVE (a permutation for every pair is CONSTRUCTED);
  §3  ★a table: how many invariant subsets of axes exist at each q = 0..d+1;
  §4  ★THE PRICE of derivability: the largest subgroup under which one axis is invariant is
      a stabilizer of index exactly d+1, that is, exactly the choice itself (a structural
      circularity).
"""
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'test'))

import sympy as sp
from itertools import permutations
from _teeth import ok, report


# ────────────────────────────────────────────── worlds ────────────────────────

def gram_democratic(d):
    """The alphabet's Gram: 1 on the diagonal, −1/d off it (the cell A_d, T19)."""
    n = d + 1
    return sp.Matrix(n, n, lambda i, j: sp.Integer(1) if i == j else sp.Rational(-1, d))


def gram_broken(d):
    """★A NEGATIVE WORLD: democracy is already BROKEN by a norm (axis 0 has norm 2).

    The mechanism is shared with S1066-B2, row 3 («a marked axis has norm 2 → it READS»)
    — stated here as a SHARED ANCESTOR, not as my own line.  Here the world is needed only
    as a tooth.
    """
    G = gram_democratic(d)
    G[0, 0] = sp.Integer(2)
    return G


def perm_matrix(p):
    n = len(p)
    return sp.Matrix(n, n, lambda i, j: sp.Integer(1) if p[j] == i else sp.Integer(0))


def symmetry_group(G):
    """MEASURING the symmetry group of the Gram: all permutations with PᵀGP = G (an exhaustive search, not belief)."""
    n = G.shape[0]
    return [p for p in permutations(range(n))
            if sp.simplify(perm_matrix(p).T * G * perm_matrix(p) - G) == sp.zeros(n, n)]


# ──────────────────────────────────────────── detectors ──────────────────────

def invariant_subsets(group, n):
    """All subsets of axes invariant under THE ENTIRE group (an exhaustive search over 2ⁿ)."""
    out = []
    for mask in range(1 << n):
        S = frozenset(k for k in range(n) if mask >> k & 1)
        if all(frozenset(p[k] for k in S) == S for p in group):
            out.append(S)
    return out


def only_trivial_invariants(world):
    """★DEFINITION (carved BEFORE the count):

    «The choice is DERIVABLE» ⟺ there exists a PROPER NONEMPTY subset of axes, invariant
    under the symmetry group of the input (because any construction from an invariant input
    is equivariant, and its output is an invariant subset).
    The detector says TRUE when only ∅ and the WHOLE set are invariant ⟹ the choice is NOT
    derivable.
    """
    G = world['G']
    n = G.shape[0]
    inv = invariant_subsets(symmetry_group(G), n)
    return len(inv) == 2 and frozenset() in inv and frozenset(range(n)) in inv


# ───────────────────────────────────────────────── the run ────────────────────

def main():
    print("=" * 78)
    print("S1075 · AX-dimer · D3 (CHOICE): is there an irreducibility theorem for the choice?")
    print("=" * 78)

    for d in range(2, 5):
        n = d + 1
        G = gram_democratic(d)
        grp = symmetry_group(G)

        print(f"\n───── d={d} ({n} axes) ─────")
        print(f"§1 · the symmetry group of the Gram MEASURED exhaustively: |Aut| = {len(grp)} "
              f"= {n}! ⟹ the input really is S_{{{n}}}-invariant "
              f"({'✓' if len(grp) == sp.factorial(n) else '✗'})")
        assert len(grp) == sp.factorial(n)

        print("§2 · transitivity — a permutation is constructed for EVERY pair of axes:")
        pairs = 0
        for i in range(n):
            for j in range(n):
                wit = [p for p in grp if p[i] == j]
                assert wit, f"no permutation {i}→{j}"
                pairs += 1
        print(f"     {pairs} pairs (i→j) — a witness found for every one ⟹ the action is TRANSITIVE")

        print("§3 · ★A TABLE of invariant subsets of axes by size q:")
        inv = invariant_subsets(grp, n)
        by_q = {q: len([S for S in inv if len(S) == q]) for q in range(n + 1)}
        for q in range(n + 1):
            tag = ("  ← the canonical q=1: NO invariants" if q == 1 else
                   ("  (trivial)" if q in (0, n) else ""))
            print(f"     q={q}: invariant subsets = {by_q[q]}{tag}")
        assert by_q[0] == 1 and by_q[n] == 1
        assert all(by_q[q] == 0 for q in range(1, n))

        ok(only_trivial_invariants, {'name': f'democracy d={d}', 'G': G},
           f"d={d}: only ∅ and the WHOLE set are invariant ⟹ no proper nonempty "
           f"mark (any 0<q<d+1) is the output of an equivariant construction",
           must_fail_on=[("democracy broken by a norm (axis 0 has norm 2)",
                          {'name': 'broken', 'G': gram_broken(d)})])

        print("§4 · ★THE PRICE OF DERIVABILITY: under which subgroup does one axis become invariant?")
        stab = [p for p in grp if p[0] == 0]
        keep0 = [p for p in grp if frozenset({p[0]}) == frozenset({0})]
        assert stab == keep0, "invariance of {0} ≠ the stabilizer — check"
        idx = len(grp) // len(stab)
        print(f"     the largest subgroup with an invariant axis {{0}} = the stabilizer, "
              f"|H| = {len(stab)} = {d}!  ·  INDEX = {idx} = d+1")
        assert idx == n
        print(f"     ⟹ to DERIVE the choice, one must already hand over exactly {n} = d+1 options,")
        print(f"       that is, exactly the choice itself.  A structural circularity, the index exact.")

        brk = symmetry_group(gram_broken(d))
        print(f"     control (the broken world): |Aut| = {len(brk)} = {d}! ⟹ there {{0}} "
              f"IS invariant, and the choice really IS DERIVABLE — the mechanism shared "
              f"with S1066-B2 row 3 (stated as an ancestor, not as my own line)")

    code = report("S1075 · D3-choice of AX-dimer")

    print("""
────────────────────────────────────────────────────────────────────────────────
WHAT THIS MEANS FOR D3 (a measurement, not a verdict):

  (1) THE IRREDUCIBILITY THEOREM OF THE CHOICE HOLDS, and the mechanism is exactly the one Ω
      wagered on: the transitivity of S_{d+1} on d+1 axes ⟹ the invariant subsets are only
      {∅, all} ⟹ no equivariant construction can hand back a single axis.
      ★This is STRONGER than «4 failures = a theorem» (a rhyme to A-nonderiv): this is not a
      list of failed attempts, but a structural impossibility for ANY construction from this
      input.

  (2) ★STRONGER THAN THE EXANTE'S FORMULATION: the theorem forbids not only «which axis
      exactly», but ANY proper nonempty mark — the table of §3 gives zero invariant subsets
      at EVERY 0 < q < d+1, not only at q=1.  ⟹ an equivariant construction can only hand
      back exactly q=0 or q=d+1; the canonical q=1 is not its output, BY TYPE.

  (3) ★A CONSEQUENCE FOR THE WITNESS COUNT (multiplicity by ancestor, not by line):
      D1 (existence) and D3 (choice) fall from ONE theorem — equivariance.
      ⟹ they are NOT two independent witnesses of irreducibility, but ONE.  What differs
      between them is only the INFORMATION (existence = 1 bit ⊥ choice = log₂(d+1) bits),
      not the mechanism.  The only thing that remains independent is the cone-measurement
      D1-§2 (it lives in the arena and knows nothing of S_{d+1}) — THAT is the second
      witness, not D3.

  (4) THE PRICE, named precisely (§4): for the choice to become derivable, democracy must be
      broken down to the stabilizer — a subgroup of INDEX exactly d+1.  That is, one must
      pay exactly the information one wants to obtain.  This is not "hard", it is circular.

  (5) ⟹ the choice = SPONTANEITY, the same class as ±m₀ and the 1-bit of T36-(v) (as Ω
      wagered).  But with the correction of (2): what is spontaneous is not only WHICH, but
      also WHETHER.
────────────────────────────────────────────────────────────────────────────────""")
    return code


if __name__ == '__main__':
    sys.exit(main())
