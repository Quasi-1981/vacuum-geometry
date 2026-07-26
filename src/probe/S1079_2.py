# -*- coding: utf-8 -*-
# DIM: symbolic-d  (the classification of 2-planes — symbolic ∀(p,q); a run over
#                   (2,1)·(3,1)·(2,2)·(3,0)·(4,0).  No d=3 instances on the load-bearing path.)
"""S1079 · SIEGE OF AX-closure · component **Z2 (THE CLASS OF MOVES)** — a lane-A probe.

EXANTE: `hub/prime/AX_CLOSURE_SIEGE_EXANTE.md` (carved by Ω, S1078).
The Z2 question verbatim: «does the verdict hold under a variation of the class of moves, or
is the class a hidden switch?»  Ω's wager: the INHERITANCE of direction; «this is the ONE
place where I expect possible new content».

★THE CARVED CLASS OF MOVES being twisted (W28 §1, verbatim):
  «generators of isometries of the COORDINATE 2-planes of the carrier: J-type (definite-sign
   pairs (+,+)/(−,−); a compact plane, closes) and K-type (mixed pairs (+,−);
   noncompact, does not close)»

★THE LEGALITY OF THE TWIST (the fishing guard of §2 of the exante — check me against it):
  twisted with ONE tool — the classification of a 2-plane by the RESTRICTION OF THE FORM to
  it (Gram/Sylvester) + bracket-closure.  New constants: 0, new entities: 0, no targeted
  tuning toward a desired verdict (the verdict (3,1) does not even appear in this probe).

★THE CYCLE GUARD OF §4: `AX-indef` forces nothing — indefinite carriers enter as a SCAN
PARAMETER alongside definite ones ((3,0)/(4,0) are measured by the same hands).
★THE FENCE OF §10: frozen words are not used.
"""
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'test'))

import sympy as sp
from itertools import combinations
from _teeth import ok, report


# ─────────────────────────────── the machinery of the carrier ──────────────────────────────

def eta(p, q):
    return sp.diag(*([sp.Integer(1)] * p + [sp.Integer(-1)] * q))


def wedge(x, y, e):
    """A generator of isometries of the plane span{x,y}: X = x(ηy)ᵀ − y(ηx)ᵀ ∈ so(η)."""
    return x * (e * y).T - y * (e * x).T


def in_so(X, e):
    return sp.simplify(X.T * e + e * X) == sp.zeros(X.rows, X.rows)


def unit(n, i):
    v = sp.zeros(n, 1)
    v[i] = sp.Integer(1)
    return v


def plane_type(x, y, e):
    """★THE CLASSIFICATION OF A 2-PLANE by the RESTRICTION of the form (a 2×2 Gram + Sylvester).

    Three — and exactly three — types for a nondegenerate carrier:
      'J' definite-sign (compact)  ·  'K' mixed nondegenerate  ·  'N' DEGENERATE
    """
    G = sp.Matrix(2, 2, lambda i, j: ((([x, y][i]).T * e * [x, y][j])[0, 0]))
    G = sp.simplify(G)
    det = sp.simplify(G.det())
    if det == 0:
        return 'N'
    if det > 0:
        return 'J'
    return 'K'


def killing_signature(basis):
    """The signature of the Killing form of the subalgebra (n₊, n₋, n₀) — a class invariant (the S900 rig)."""
    m = len(basis)
    K = sp.zeros(m, m)
    for i in range(m):
        for j in range(m):
            K[i, j] = sp.trace(basis[i] * basis[j])
    K = sp.simplify(K)
    evs = K.eigenvals()
    pos = sum(mult for val, mult in evs.items() if sp.simplify(val) > 0)
    neg = sum(mult for val, mult in evs.items() if sp.simplify(val) < 0)
    zer = sum(mult for val, mult in evs.items() if sp.simplify(val) == 0)
    return (int(pos), int(neg), int(zer))


def bracket_closure(gens, e):
    """Lie(V): closure under the bracket (the same mechanism as in (b))."""
    def flat(M):
        return sp.Matrix(M.rows * M.rows, 1, list(M))
    basis = []
    M = sp.zeros(gens[0].rows ** 2, 0)
    for X in gens + []:
        cand = M.row_join(flat(X))
        if cand.rank() > M.rank():
            M = cand
            basis.append(X)
    changed = True
    while changed:
        changed = False
        for A, B in combinations(basis, 2):
            C = sp.simplify(A * B - B * A)
            cand = M.row_join(flat(C))
            if cand.rank() > M.rank():
                M = cand
                basis.append(C)
                changed = True
    assert all(in_so(X, e) for X in basis)
    return basis


# ──────────────────────────── detectors (functions of a WORLD) ───────────────────

def carved_class_is_complete(world):
    """★Does the CARVED class of moves exhaust all types of 2-planes of the carrier?

    TRUE ⟺ among all 2-planes of the carrier there is NOT A SINGLE degenerate one ('N'),
    that is, {J,K} is a complete list of types.  The detector looks at the WORLD (the
    carrier), not at my own list.
    """
    p, q = world['p'], world['q']
    e = eta(p, q)
    n = p + q
    if q == 0:
        return True
    # searching for a degenerate 2-plane explicitly: an isotropic vector + one orthogonal to it
    u = unit(n, 0) + unit(n, n - 1)                       # η(u,u) = 1 − 1 = 0
    for k in range(1, n - 1):
        w = unit(n, k)
        if plane_type(u, w, e) == 'N':
            return False
    return True


def coordinate_planes_are_nondegenerate(world):
    """A control: COORDINATE 2-planes in the diagonal basis — never 'N'."""
    p, q = world['p'], world['q']
    e = eta(p, q)
    n = p + q
    return all(plane_type(unit(n, i), unit(n, j), e) != 'N'
               for i, j in combinations(range(n), 2))


# ───────────────────────────────────────────────── the run ─────────────────────

def main():
    print("=" * 78)
    print("S1079 · AX-closure · Z2 (THE CLASS OF MOVES): is the class a hidden switch?")
    print("=" * 78)

    WORLDS = [(2, 1), (3, 1), (2, 2), (3, 0), (4, 0)]

    # ── §1 · how many plane types the carrier has, and how many the carve sees ─────
    print("\n§1 · TYPES OF 2-PLANES of the carrier (Gram restriction + Sylvester) vs the carved class")
    for (p, q) in WORLDS:
        e, n = eta(p, q), p + q
        coord = {}
        for i, j in combinations(range(n), 2):
            t = plane_type(unit(n, i), unit(n, j), e)
            coord[t] = coord.get(t, 0) + 1
        deg = "YES" if not carved_class_is_complete({'p': p, 'q': q}) else "none"
        print(f"   ({p},{q}): coordinate planes → {coord}"
              f"   ·   degenerate (N) planes in the carrier: {deg}")

    for (p, q) in [(3, 0), (4, 0)]:
        ok(carved_class_is_complete, {'name': f'({p},{q}) definite', 'p': p, 'q': q},
           f"({p},{q}): the carved class {{J,K}} is COMPLETE — no degenerate 2-plane "
           f"exists in a definite carrier",
           must_fail_on=[("(3,1) — an indefinite carrier", {'name': '(3,1)', 'p': 3, 'q': 1}),
                         ("(2,1)", {'name': '(2,1)', 'p': 2, 'q': 1})])

    for (p, q) in WORLDS:
        assert coordinate_planes_are_nondegenerate({'p': p, 'q': q})
    print("""   ✓ a control: COORDINATE 2-planes in the diagonal basis are never
     degenerate (checked on all 5 carriers) ⟹ the carve «coordinate planes»
     EXCLUDES the N-type NOT by content, but BY THE CHOICE OF BASIS.""")

    # ── §2 · the N-move exists as an object of the same machinery ──────────────────────
    print("\n§2 · THE N-MOVE (a generator of isometries of a DEGENERATE plane) — built explicitly")
    for (p, q) in [(2, 1), (3, 1)]:
        e, n = eta(p, q), p + q
        u = unit(n, 0) + unit(n, n - 1)
        w = unit(n, 1)
        assert sp.simplify((u.T * e * u)[0, 0]) == 0
        N = wedge(u, w, e)
        assert in_so(N, e)
        pw = [sp.simplify(N ** k) for k in (2, 3)]
        nilp = next((k for k, M in zip((2, 3), pw) if M == sp.zeros(n, n)), None)
        print(f"   ({p},{q}): N ∈ so(η) ✓ · the plane type = "
              f"{plane_type(u, w, e)} · nilpotent: N^{nilp} = 0 · "
              f"Killing span{{N}} = {killing_signature([N])} · "
              f"bracket-closed (1-dim abelian) ⟹ a TERMINAL in (b)")

    # ── §3 · does the N-move add NEW classes (i.e. is the class of moves a switch) ─────
    print("\n§3 · ★DOES THE N-MOVE GIVE A CLASS THE COORDINATE LIST DOES NOT SEE")
    p, q = 2, 1
    e, n = eta(p, q), p + q
    u = unit(n, 0) + unit(n, n - 1)
    N = wedge(u, unit(n, 1), e)
    D = wedge(unit(n, 0), unit(n, n - 1), e)        # a K-move (a boost) in the same plane
    pair = bracket_closure([N, D], e)
    kill = killing_signature(pair)
    print(f"   (2,1): Lie{{N, K}} → dim = {len(pair)} · Killing = {kill}")
    print(f"           [D,N] = {sp.simplify(D * N - N * D).tolist()}")

    coord_classes_21 = {(1, (0, 0, 1)), (3, (2, 1, 0))}     # S901, block A for (2,1)
    print(f"   the coordinate classes of (2,1) from S901-A: {sorted(coord_classes_21)}")
    print(f"   our class: {(len(pair), kill)} — "
          f"{'OUTSIDE the coordinate list' if (len(pair), kill) not in coord_classes_21 else 'inside'}")

    ok(lambda w: (w['dim'], w['kill']) not in w['coord'],
       {'name': 'Lie{N,K} in (2,1)', 'dim': len(pair), 'kill': kill,
        'coord': coord_classes_21},
       "★the closure {N-move, K-move} gives a class OUTSIDE the coordinate list S901-A "
       "⟹ the carved class of moves does NOT exhaust the terminals of its own carrier",
       must_fail_on=[("a coordinate class {1·(0,0,1)}",
                      {'name': 'coordinate', 'dim': 1, 'kill': (0, 0, 1),
                       'coord': coord_classes_21}),
                     ("a coordinate class {3·(2,1,0)} (the full so(2,1))",
                      {'name': 'full', 'dim': 3, 'kill': (2, 1, 0),
                       'coord': coord_classes_21})])

    code = report("S1079 · Z2-class of moves of AX-closure")

    print("""
────────────────────────────────────────────────────────────────────────────────
WHAT THIS MEANS FOR Z2 (a measurement, not a verdict):

  (1) ★THE CLASS OF MOVES IS A SWITCH — but not the one that was feared.  It does not
      change the VERDICT-NUMBER (that lives on the premises, Z1-§1); it changes the LIST OF
      TERMINALS.  The carve «COORDINATE 2-planes» excludes a third plane type — DEGENERATE
      (N) — and does so NOT by content, but BY THE CHOICE OF BASIS: in the diagonal basis a
      coordinate plane is never degenerate (measured on 5 carriers).

  (2) ★THERE ARE EXACTLY THREE TYPES, AND THE THIRD EXISTS ONLY AT THE BREAK: restricting
      the form to a 2-plane gives {definite-sign J · mixed nondegenerate K · DEGENERATE N},
      and N exists ⟺ the carrier is indefinite.  The tooth is exactly this: on (3,0)/(4,0)
      the detector says «the class is complete», on (2,1)/(3,1) — «incomplete».  ⟹ the
      carved class {J,K} is COMPLETE for a definite carrier and INCOMPLETE exactly where our
      arena lives.

  (3) ★THIS CLOSES A NAMED LIMIT OF THE ANCESTOR WITH A MECHANISM (not a new claim):
      S901 honestly carved that block-B (random subspaces) finds classes OUTSIDE the
      coordinate list (a Borel-type), and that a full classification of subalgebras was not
      done.  Measured: Lie{N-move, K-move} in (2,1) gives exactly a class outside the
      coordinate list.  ⟹ the source of those "extra" classes is NOT a scan accident, but a
      SYSTEMATICALLY omitted type of move.

  (4) ⟹ the form of the verdict handed to Ω (not rendered here): on Z2 the wager
      "inheritance of direction" is NOT refuted by this probe (the J/K direction was not
      touched), but the CARVED phrase itself, «coordinate 2-planes», is measured to be
      BASIS-DEPENDENT and incomplete on an indefinite carrier.  This is either (i) a
      narrowing of the carve that must be stated honestly in the preprint, or (ii) an
      assignment to fill in the N-row of table O1 (the direction for {N,N}/{N,J}/{N,K}) —
      but that is already a NEW assignment, not my choice.

  BOUNDARIES: the direction (who forces whom) for the N-type was NOT measured — that is a
  separate table.  It is not claimed that the verdict (3,1) moves: it does not even appear
  in this probe.
────────────────────────────────────────────────────────────────────────────────""")
    return code


if __name__ == '__main__':
    sys.exit(main())
