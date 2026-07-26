# -*- coding: utf-8 -*-
# DIM: symbolic-exact (rational matrices; exact ranks).  Numbers with a bracket
#      [address · unit · type/operation].
"""S1066 B1 — THE AX-indef SIEGE, step 1: KILL-FIRST «EUCLID IS SELF-SUFFICIENT».

ASSIGNMENT: project ruling, exante `hub/prime/AX_INDEF_SIEGE_EXANTE.md`, task B1 (order B1 → STOP).
THE NULL HYPOTHESIS (I strike precisely at it): «the full Euclidean (q=0) version of the closure/socket
machinery T1–T4 exists WITHOUT an internal contradiction».  Binary: does something break STRUCTURALLY (to be named
a theorem) — or does everything hold.

ANCESTORS (by citation, not re-derived; addresses from `codex/graph.json`):
  · AX-closure — the CLOSURE construction of the order ≺ on real quadratic spaces
    (a method-input, NOT a numeric handle) [§3.arch-0 · preprint-1 · S899+];
  · T1 — the minimal carrier of the closure = SIGN-DEFINITE (3,0); there is no 2-dim non-abelian so(3);
  · T2 — the direction of the order: the J-sector is bracket-closed; a full K forces so(p,q);
  · T3 — the antisymmetry of ≺ (idempotency of the closure, no cycles);
  · T4 — the enumeration of sockets: the coordinate terminal is finite (2→14 classes);
  · T5 — minimality → (3,1): `p≥3 ∧ q≥1 ⟹ n≥4`; deps = {T4, AX-cell, **AX-indef**};
  · AX-indef — the POSTULATE q≥1 (a yellow root), consumed by — T5.

★The project's WAGER, CARVED BEFORE THE COUNT (I strike at it, I do not confirm it): «Euclid is ALIVE as
existence — time in the program is an OUTPUT, not a precondition; a system without a minus is not contradictory, it
is simply CLOCK-LESS».

★WHAT I FEAR (carved BEFORE the count): that I will mistake a DEGENERACY for a CONTRADICTION.  «Empty» and
«breaks» — are different outcomes, and the first does NOT demote AX-indef.  So every failure must be
presented as a THEOREM («exactly what is incompatible»), not as an empty enumeration.

A K-WINNING-DRIFT (the exante): even if Euclid falls — EXACTLY AX-indef is demoted;
the exclusivity of (3,1) is NOT proved, the p-axis does not move.

FENCE: Layer-1.  The physics-vocabulary classes named in the project's fence do not enter the code
(K-stone: the signature is measured, the GROUP is not).  Homonyms only with prefixes.  No verdict is rendered.
"""
import itertools
import json
import os
import sys

import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "..", "test")))
from _teeth import ok, report, reset                      # noqa: E402

_GRAPH = os.path.abspath(os.path.join(_HERE, "..", "..", "codex", "graph.json"))


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


# ═══════════════ THE MACHINERY: so(p,q) EXACTLY, RATIONALLY ═══════════════
def metric(p, q):
    return sp.diag(*([sp.Integer(1)] * p + [sp.Integer(-1)] * q))


def so_basis(p, q):
    """A basis of so(p,q) = {X : XᵀG + GX = 0}, split into SECTORS:
    J — rotations INSIDE a same-sign block (compact),
    K — «mixed» generators BETWEEN blocks of different sign (exactly p·q of them).
    These are precisely the two sectors between which the order ≺ [T2] lives."""
    n = p + q
    J, K = [], []
    for i, j in itertools.combinations(range(n), 2):
        E = sp.zeros(n, n)
        same = (i < p) == (j < p)
        if same:
            E[i, j], E[j, i] = 1, -1          # antisymmetric ⟹ a rotation
            J.append(E)
        else:
            E[i, j], E[j, i] = 1, 1           # a symmetric block-bridge
            K.append(E)
    return J, K


def in_algebra(X, p, q):
    G = metric(p, q)
    return sp.simplify(X.T * G + G * X) == sp.zeros(p + q, p + q)


def vec(X):
    return sp.Matrix([X[i, j] for i in range(X.rows) for j in range(X.cols)])


def span_rank(mats):
    if not mats:
        return 0
    return sp.Matrix.hstack(*[vec(M) for M in mats]).rank()


def bracket(A, B):
    return sp.expand(A * B - B * A)


def closure(gens, cap=6):
    """The CLOSURE of a subspace under the bracket: iterate span(V ∪ [V,V]) to a fixed point.
    Returns (the dimension, whether it CONVERGED within cap iterations) — finiteness here is
    precisely the subject of T3/T4 (no cycles, the terminal is finite)."""
    cur = list(gens)
    r = span_rank(cur)
    for _ in range(cap):
        new = cur + [bracket(A, B) for A, B in itertools.combinations(cur, 2)]
        r2 = span_rank(new)
        if r2 == r:
            return r, True
        cur, r = new, r2
    return r, False


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1066_1_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf); sys.stdout = tee
    reset()

    print("=" * 80)
    print("S1066 B1 — KILL-FIRST: is EUCLID (q=0) SELF-SUFFICIENT for T1–T4")
    print("=" * 80)
    print()

    # ═══════ PART 1: AN AUDIT OF THE GRAPH — WHERE EXACTLY q≥1 ENTERS ═══════
    print("PART 1 — AN AUDIT OF THE ACCEPTED GRAPH: at WHICH node is AX-indef first needed")
    graph = json.load(open(_GRAPH, encoding='utf-8'))
    nodes = {n['id']: n for n in graph['nodes']}

    def ancestors(nid, seen=None):
        seen = seen if seen is not None else set()
        for d in nodes.get(nid, {}).get('deps', []):
            if d not in seen:
                seen.add(d)
                ancestors(d, seen)
        return seen

    def indef_free(nid):
        """A WORLD = a graph node ⟹ whether its FULL closure of ancestors does NOT contain AX-indef.
        The detector reads the ACCEPTED graph, not my paraphrase."""
        return 'AX-indef' not in ancestors(nid)

    for nid in ('T1', 'T2', 'T3', 'T4', 'T5', 'T7'):
        anc = ancestors(nid)
        print(f"     {nid}: ancestors = {sorted(anc)}  ⟹ AX-indef {'YES' if 'AX-indef' in anc else 'NO'}"
              f"  [pt.1 · dimensionless · the closure of deps in graph.json]")
    ok(indef_free, 'T4',
       "★★T1–T4 DO NOT DEPEND on AX-indef: the entire closure of their ancestors = {AX-closure} ⟹ "
       "Euclid is self-sufficient ON THIS STRETCH by the construction of the graph  [pt.1 · dimensionless · "
       "the closure of deps]",
       must_fail_on=[("T5 — there AX-indef IS present (it is the first consumer)", 'T5'),
                     ("T7 — inherits it through T5", 'T7')])
    print("     ⟹ ★q≥1 enters the arc EXACTLY ONCE and EXACTLY at T5.  This is not my paraphrase —")
    print("       this is the closure of deps of the accepted graph, read by the machine.")
    print()

    # ═══════ PART 2: DOES THE MACHINERY BREAK AT q=0 ═══════
    print("PART 2 — does anything break STRUCTURALLY at a sign-definite form (q=0)")

    def naive_basis(p, q):
        """A DIFFERENT LAW of construction: always an antisymmetric generator, WITHOUT regard to the metric
        (that is, the «Euclidean» rule, applied even where a minus exists)."""
        n = p + q
        out = []
        for i, j in itertools.combinations(range(n), 2):
            E = sp.zeros(n, n)
            E[i, j], E[j, i] = 1, -1
            out.append(E)
        return out

    def basis_is_algebra(world):
        """A WORLD = (p, q, A CONSTRUCTION RULE) ⟹ whether the basis lies in so(p,q) and whether dim = n(n−1)/2.
        This is a check of MY OWN machinery BEFORE the verdict.  ★A negative world — a DIFFERENT RULE (naive,
        metric-blind), not different (p,q): my own pattern `negative-world-is-a-different-law`."""
        p, q, rule = world
        n = p + q
        B = (lambda: (lambda JK: JK[0] + JK[1])(so_basis(p, q))) () if rule == 'metric' \
            else naive_basis(p, q)
        return (all(in_algebra(X, p, q) for X in B)
                and span_rank(B) == n * (n - 1) // 2)

    ok(basis_is_algebra, (3, 0, 'metric'),
       "★the machinery holds in the SIGN-DEFINITE sector: the basis of so(3,0) lies in the algebra, "
       "dim = 3 = n(n−1)/2  [pt.2 · dimensionless · a check of XᵀG+GX=0 and the rank]",
       must_fail_on=[("★A DIFFERENT LAW: a metric-blind (always-antisymmetric) construction at "
                      "q=1 — the generators cease to lie in the algebra", (3, 1, 'naive')),
                     ("the same naive rule at (2,2)", (2, 2, 'naive'))])
    for pq in ((3, 0), (4, 0), (3, 1), (2, 2)):
        p, q = pq
        J, K = so_basis(p, q)
        print(f"     so({p},{q}): dim J = {len(J)} · dim K = {len(K)} (= p·q = {p*q}) · "
              f"together {len(J)+len(K)} = n(n−1)/2  [pt.2 · dimensionless · a count of generators]")

    def closure_stays_and_terminates(world):
        """A WORLD = (p, q, a CLOSURE OPERATION) ⟹ whether the iteration converges to a fixed point
        AND whether everything generated STAYS in the algebra.  The finiteness of the terminal is the subject of T3/T4;
        if it broke at q=0, this WOULD BE the sought CONTRADICTION.
        ★A negative world — a DIFFERENT OPERATION (an ordinary product instead of the bracket): it EXITS the
        algebra, that is, it breaks precisely what the bracket holds."""
        p, q, op = world
        J, K = so_basis(p, q)
        gens = (J + K)[:2]
        mul = (lambda A, B: sp.expand(A * B - B * A)) if op == 'bracket' \
            else (lambda A, B: sp.expand(A * B))
        cur, r = list(gens), span_rank(list(gens))
        for _ in range(6):
            new = cur + [mul(A, B) for A, B in itertools.combinations(cur, 2)]
            r2 = span_rank(new)
            if r2 == r:
                return all(in_algebra(X, p, q) for X in cur)
            cur, r = new, r2
        return False

    ok(closure_stays_and_terminates, (3, 0, 'bracket'),
       "★★THE BRACKET-CLOSURE CONVERGES EVEN AT q=0, without exiting the algebra ⟹ the T3/T4-finiteness "
       "does NOT rest on the minus  [pt.2 · dimensionless · a fixed point of the iteration span(V∪[V,V])]",
       must_fail_on=[("★A DIFFERENT OPERATION: an ordinary product instead of the bracket — exits "
                      "the algebra (symmetric parts)", (3, 0, 'product')),
                     ("the same product at (3,1)", (3, 1, 'product'))])

    def minimal_carrier_is_three(subspace_kind):
        """A WORLD = A CLASS OF SUBSPACE (how it is chosen) ⟹ whether the minimal bracket-closed
        carrier has dim 3.  ★A negative world — a DIFFERENT LAW of choice (an abelian plane), not different
        numbers: there the closure = the subspace itself, dim 2 (my own pattern
        `negative-world-is-a-different-law`)."""
        J, _ = so_basis(4, 0)
        if subspace_kind == 'generic':
            gens = [J[0], J[1]]                       # two axes that do NOT commute
        else:
            gens = [J[0], J[5]]                       # a Cartan pair: [J01, J23] = 0
        r, done = closure(gens)
        return done and r == 3

    ok(minimal_carrier_is_three, 'generic',
       "★T1 in the sign-definite sector: a generic 2-dimensional plane closes EXACTLY to "
       "dim 3 (so(3)) — «there is no 2-dim non-abelian one»  [pt.2 · dimensionless · the rank of the closure]",
       must_fail_on=[("★A DIFFERENT LAW of choosing the plane: an ABELIAN (Cartan) pair — "
                      "the closure = the plane itself, dim 2", 'cartan')])
    print()

    # ═══════ PART 3: WHAT EXACTLY DEGENERATES (not breaks) ═══════
    print("PART 3 — ★what EXACTLY fades at q=0: a degeneracy ⊥ a contradiction")

    def second_sector_exists(pq):
        """A WORLD = (p,q) ⟹ whether a SECOND sector (K, «mixed» generators) exists.
        The order ≺ [T2] — a relation BETWEEN the two sectors; at q=0 the second sector
        does not exist, and the relation becomes EMPTY (not contradictory — empty)."""
        p, q = pq
        _, K = so_basis(p, q)
        return len(K) > 0

    ok(second_sector_exists, (3, 1),
       "★THE SECOND SECTOR (K) exists ⟺ q≥1: |K| = p·q  [pt.3 · dimensionless · a count of "
       "mixed generators]",
       must_fail_on=[("q=0, (3,0) — K is EMPTY", (3, 0)),
                     ("q=0, (4,0) — K is empty there too", (4, 0))])

    def killing_is_definite(pq):
        """A WORLD = (p,q) ⟹ whether the Killing form `K(X,Y) = tr(ad_X ad_Y)` ON THE ALGEBRA ITSELF
        is sign-definite.  `ad_X` is built HONESTLY: the coordinates `[X, B_k]` in the basis B —
        the solution of an exact linear system, not a textbook formula."""
        p, q = pq
        J, K = so_basis(p, q)
        B = J + K
        m = len(B)
        Mb = sp.Matrix.hstack(*[vec(X) for X in B])

        def coords(Z):
            sol = Mb.solve_least_squares(vec(Z)) if Mb.rows != Mb.cols else Mb.solve(vec(Z))
            return sp.Matrix([sp.simplify(c) for c in sol])

        ads = []
        for X in B:
            cols = [coords(bracket(X, Y)) for Y in B]
            ads.append(sp.Matrix.hstack(*cols))
        Kf = sp.zeros(m, m)
        for i in range(m):
            for j in range(i, m):
                Kf[i, j] = Kf[j, i] = sp.simplify((ads[i] * ads[j]).trace())
        ev = [sp.simplify(v) for v in Kf.eigenvals()]
        signs = {int(sp.sign(sp.N(v, 25))) for v in ev if sp.simplify(v) != 0}
        return len(signs) == 1

    ok(killing_is_definite, (3, 0),
       "★★THE INDEFINITENESS OF THE ARENA IS READ EVEN ON THE ALGEBRA: the Killing form is SIGN-DEFINITE "
       "at q=0 and INDEFINITE at q≥1  [pt.3 · dimensionless · the signs of the eigenvalues of "
       "tr(ad·ad)]",
       must_fail_on=[("(3,1) — there the Killing form is INdefinite", (3, 1)),
                     ("(2,2) — likewise indefinite", (2, 2))])
    print()

    code = report("S1066 B1 — is Euclid self-sufficient?")
    print()
    print("=" * 80)
    print("★RAW OUTPUTS OF B1 (no verdict rendered — a B-visa, then the project's court)")
    print("  PART 1 (a graph audit): AX-indef is ABSENT from the closure of ancestors of T1–T4 and")
    print("     PRESENT in T5/T7 ⟹ q≥1 enters the arc EXACTLY ONCE, at T5.")
    print("  PART 2 (the machinery at q=0): nothing breaks —")
    print("     · the basis of so(n,0) lies in the algebra, dim = n(n−1)/2;")
    print("     · the closure CONVERGES (a fixed point) ⟹ the T3/T4-finiteness does not")
    print("       rest on the minus;")
    print("     · T1 holds: a generic plane closes exactly to dim 3.")
    print("  PART 3 (what fades): the SECOND SECTOR K is empty at q=0 (|K| = p·q) ⟹")
    print("     the order ≺ [T2] — a relation BETWEEN the sectors — becomes EMPTY.")
    print("     The Killing form at q=0 is SIGN-DEFINITE; at q≥1 it is indefinite.")
    print("  ⟹ ★THE NULL HYPOTHESIS B1 WITHSTOOD: there is no contradiction, there is a DEGENERACY.")
    print("     «Empty» ≠ «breaks» — and I do NOT present the former as the latter.")
    print("     ⟹ there is NO unconditional forcing of q≥1 from this stretch; the siege must go")
    print("       into B2 (an implication between roots), just as the project's wager carved.")
    print("=" * 80)
    sys.stdout = tee.real
    logf.close()
    return code


if __name__ == "__main__":
    sys.exit(main())
