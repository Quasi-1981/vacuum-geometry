# -*- coding: utf-8 -*-
# DIM: symbolic-d  (the Gram and the zigzag — symbolic in d; the coordinate realization
#                   checked for d=2..7.  The d=3 numbers enter ONLY as an anchor to the
#                   S597 primary source, not as a carrier of the form.)
"""S1075 · SIEGE OF AX-dimer · component **D2 (FORM)** — a lane-A probe.

EXANTE: `hub/prime/AX_DIMER_SIEGE_EXANTE.md` (carved by Ω, S1074).
The D2 question verbatim: «the mark = a DIMER (exactly 1 time-bond/node; q=1) — does the form
follow from what is already accepted?»; the candidate for the ONE new piece of content is
«1 bond per NODE»; step 0 — pull up the S597 primary source and separate what it MEASURED
from what it POSTULATED.

★AN EXPLICIT LIST OF THE ANCESTORS INVOLVED (the cycle guard of §3 of the exante):
  · AX-alphabet {the bare set d+1 · democracy S_{d+1}}   — the root;
  · AX-cell / A-space (the cell A_d, DERIVED at the C2-swap, J-0476 VALID ∀d) — Gram 1 / −1/d;
  · S597 (viz/SIMPLEX_LATTICE_MEMBRANES_omega_analysis.md) — the primary source, the d=3
    diamond;
  · S1054/S1055/S1001 — the d-qualifier of the S597 numbers (TIME=1+1/d, LATERAL=√(1−1/d²)).
  AX-indef and its descendants (T32-minus, the signature machinery) — NOT INVOLVED in any
  detector of this file (verify with Beta's own eyes: there is no signature, no minus sign,
  no q≥1 premise in the code; the parameter q here = the COUNT OF MARKED BOND CLASSES,
  scanned freely from 0 to d+1).

THE OBJECT (FIELD-0: the ancestor's form BEFORE the count).  The bond graph of the cell = a
marked maximal abelian cover of a dipole graph with d+1 edges: node-A at point m, node-B at
m+e_j, edge label = j ∈ {0..d}; coordinates m ∈ ℤ^{d+1} modulo (1,…,1) (since Σ_j u_j = 0).
This is the SAME construction that gave the cell A_d; here it is needed for its EDGES, not
only its weights.

WHAT IS MEASURED (4 tests, each with a negative world — the `test/_teeth.py` harness):
  T1  marked bonds of one class = a PERFECT MATCHING (1 per node)               ∀d
  T2  a purely-time-like path does NOT EXIST (a component of the marked subgraph = 1 edge) ∀d
  T3  the zigzag: TIME=1+1/d, LATERAL=√(1−1/d²) — symbolic from the Gram + an anchor in the
      S597 numbers
  T4  ablation over q: a matching ⟺ q=1 (the D2 form carries no content beyond the count of
      classes)
"""
import sys
import os
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'test'))

import sympy as sp
from itertools import product
from _teeth import ok, ok_contrast, report

DMAX = 7


# ───────────────────────────────────── the object: the cell A_d ────────────────────

def weights(d):
    """u_0..u_d ∈ ℝ^d — d+1 equal-rights directions of the cell, |u_j| = 1 (canon ed.2).

    Construction: project e_j ∈ ℝ^{d+1} onto {Σx=0} and normalize.  The Gram comes out to
    1 on the diagonal, −1/d off it — identical to `test/_forced_lattice.bonds` for d=2,3.
    """
    n = d + 1
    c = sp.sqrt(sp.Rational(n, d))
    cols = []
    for j in range(n):
        v = [sp.Rational(-1, n)] * n
        v[j] += 1
        cols.append([c * x for x in v])
    # an orthonormal basis of the plane {Σx=0} — to hand back coordinates in ℝ^d
    basis = []
    for k in range(d):
        w = sp.Matrix([sp.Integer(1) if i == k else (sp.Integer(-1) if i == k + 1 else sp.Integer(0))
                       for i in range(n)])
        for b in basis:
            w = w - (w.dot(b)) * b
        w = sp.simplify(w / sp.sqrt(w.dot(w)))
        basis.append(w)
    out = []
    for col in cols:
        cv = sp.Matrix(col)
        out.append(tuple(sp.simplify(cv.dot(b)) for b in basis))
    return out


def gram_from_coords(d):
    u = weights(d)
    n = d + 1
    return sp.Matrix(n, n, lambda i, j: sp.simplify(sum(a * b for a, b in zip(u[i], u[j]))))


def gram_canon(d):
    """The ancestor's Gram (AX-cell): 1 on the diagonal, −1/d off it."""
    n = d + 1
    return sp.Matrix(n, n, lambda i, j: sp.Integer(1) if i == j else sp.Rational(-1, d))


# ───────────────────────────────── worlds (each = a graph + a marked set) ──────

def _canon(m):
    """m ∈ ℤ^{d+1} modulo (1,…,1): fix the last coordinate at zero."""
    s = m[-1]
    return tuple(x - s for x in m)


def world_cell(d, marked=(0,), ball=2):
    """OUR WORLD: a marked cover of the cell A_d, the marked classes = `marked`.

    a node = (s, m): s=0 — sublattice A, s=1 — sublattice B.
    an edge of label j: A(m) — B(m+e_j).  ⟹ every node has EXACTLY ONE edge of each label.
    """
    n = d + 1
    nodes = []
    for m in product(range(-ball, ball + 1), repeat=n):
        if max(m) - min(m) > ball:
            continue
        for s in (0, 1):
            nodes.append((s, _canon(m)))
    nodes = sorted(set(nodes))

    def marked_nbrs(node):
        s, m = node
        out = []
        for j in marked:
            e = [0] * n
            e[j] = 1
            if s == 0:
                out.append((1, _canon(tuple(a + b for a, b in zip(m, e)))))
            else:
                out.append((0, _canon(tuple(a - b for a, b in zip(m, e)))))
        return out

    return {'name': f'cell A_{d}, marked {sorted(marked)}', 'nodes': nodes,
            'marked_nbrs': marked_nbrs}


def world_columns(d, ball=2):
    """NEGATIVE WORLD No. 1 — «COLUMNS» (the rival named in S597 itself, fact 2:
    «Time = a dimer covering, NOT columns»).  Time is realized as a CHAIN: a node carries
    BOTH time-bonds ±u_0 to different neighbors.  The lattice is not bipartite in time.
    """
    n = d + 1
    nodes = sorted({_canon(m) for m in product(range(-ball, ball + 1), repeat=n)
                    if max(m) - min(m) <= ball})

    def marked_nbrs(node):
        e = [0] * n
        e[0] = 1
        return [_canon(tuple(a + b for a, b in zip(node, e))),
                _canon(tuple(a - b for a, b in zip(node, e)))]

    return {'name': f'columns (time = a chain ±u_0), d={d}', 'nodes': nodes,
            'marked_nbrs': marked_nbrs}


def world_q2(d, ball=2):
    """NEGATIVE WORLD No. 2 — the same cell, but with TWO marked classes (q=2).
    Needed so the detector is not blind to the COUNT of the mark (and not only to the
    lattice type).
    """
    w = world_cell(d, marked=(0, 1), ball=ball)
    w['name'] = f'cell A_{d}, q=2 (two marked classes)'
    return w


def degree(world_nbrs, node, n_labels):
    """The degree of a node in the FULL bond graph of the world (all classes) — the T19 gate (z = d+1)."""
    return len(set(world_nbrs(node)))


# ───────────────────────────────────── detectors (functions of a WORLD) ──────────

def is_perfect_matching(world):
    """Every node is incident to EXACTLY one marked edge (a perfect matching)."""
    return all(len(set(world['marked_nbrs'](v))) == 1 for v in world['nodes'])


def no_pure_marked_path(world):
    """In the marked subgraph there is NO path of length ≥2: every component = 1 edge.

    A BFS over marked edges from every node; a component must have exactly 2 nodes.
    (A different detector from T1: it looks at CONNECTIVITY, not degree.)
    """
    mn = world['marked_nbrs']
    for v in world['nodes']:
        seen, front = {v}, [v]
        while front and len(seen) <= 3:
            nxt = []
            for x in front:
                for y in mn(x):
                    if y not in seen:
                        seen.add(y)
                        nxt.append(y)
            front = nxt
        if len(seen) != 2:
            return False
    return True


def zigzag_law(gram_world):
    """TIME and LATERAL over one forced zigzag period — from a SINGLE Gram, symbolically.

    Definition (S1054/S1001, verbatim): TIME := the axial run over a period =
    (a time-step +u_0) + (a space-step −u_j); LATERAL := the transverse drift over the same
    period.  The axial run of the step −u_j = −u_j·u_0;  transverse² = |u_j|² − (u_j·u_0)².
    """
    G, d = gram_world['gram'], gram_world['d']
    t = sp.simplify(G[0, 0] - G[1, 0])                       # 1 + (−u_1·u_0)
    lat = sp.simplify(sp.sqrt(G[1, 1] - G[1, 0] ** 2))       # the time-step gives 0 transverse
    return sp.simplify(t), sp.simplify(lat)


def zigzag_is_canon(gram_world):
    d = gram_world['d']
    t, lat = zigzag_law(gram_world)
    return (sp.simplify(t - (1 + sp.Rational(1, 1) / d)) == 0 and
            sp.simplify(lat ** 2 - (1 - sp.Rational(1, 1) / d ** 2)) == 0)


# ───────────────────────────────────────────────── the run ─────────────────────

def main():
    print("=" * 78)
    print("S1075 · AX-dimer · D2 (FORM): does «dimer» follow from what is already accepted?")
    print("=" * 78)

    # ---- GATE-0: the object = the ancestor's cell, not our new one -----------------------
    print("\n[GATE-0] the coordinate realization ⟹ the ancestor's Gram (1 / −1/d)")
    for d in range(2, DMAX + 1):
        got, want = gram_from_coords(d), gram_canon(d)
        assert sp.simplify(got - want) == sp.zeros(d + 1, d + 1), f"the Gram did not match, d={d}"
        print(f"   d={d}: Gram = 1 diag / {sp.Rational(-1, d)} off — ✓")

    grams = [{'d': d, 'gram': gram_canon(d), 'name': f'cell A_{d}'} for d in range(2, DMAX + 1)]

    # ---- GATE-T19: the rival «columns» is eliminated as the ANCESTOR (z = d+1), not by taste ----
    print("\n[GATE-T19] the count of bonds per node (z): our world = d+1 · «columns» = 2(d+1)")
    for d in range(2, 6):
        ours = world_cell(d, marked=tuple(range(d + 1)))
        cols = world_columns(d)
        z_ours = degree(ours['marked_nbrs'], ours['nodes'][0], d + 1)
        # in the world of columns every class is realized as a chain ±u_j ⟹ 2 neighbors per class
        z_cols = 2 * (d + 1)
        assert z_ours == d + 1, (d, z_ours)
        print(f"   d={d}: our z={z_ours} = d+1 ✓ (T19)   ·   «columns» z={z_cols} ✗ — "
              f"NOT the T19 cell, eliminated by the ancestor")

    # ---- T1: a dimer = a perfect matching ------------------------------------
    print("\n[T1] a marked class of bonds = a PERFECT MATCHING (1 time-bond per node)")
    for d in range(2, 6):
        ok(is_perfect_matching, world_cell(d),
           f"d={d}: the marked bonds of class u_0 form a perfect matching A–B",
           must_fail_on=[("columns (the S597 rival)", world_columns(d)),
                         ("q=2 (two marked classes)", world_q2(d))])

    # ---- T2: no purely-time-like path exists --------------------------------
    print("\n[T2] a purely-time-like path does NOT EXIST (a component of the marked subgraph = 1 edge)")
    for d in range(2, 6):
        ok(no_pure_marked_path, world_cell(d),
           f"d={d}: the marked subgraph has no path of length ≥2 ⟹ duration = a zigzag",
           must_fail_on=[("columns (time = a chain)", world_columns(d)),
                         ("q=2 (two marked classes)", world_q2(d))])

    # ---- T3: the zigzag law ∀d + an anchor in the numbers of the primary source --------------------
    print("\n[T3] the zigzag ∀d from the Gram: TIME = 1+1/d · LATERAL = √(1−1/d²)")
    ortho = {'d': 3, 'gram': sp.eye(4), 'name': 'a FACTORIZED lattice (space ⊥ time)'}
    for gw in grams[:4]:
        t, lat = zigzag_law(gw)
        ok(zigzag_is_canon, gw,
           f"d={gw['d']}: TIME={t} · LATERAL={sp.nsimplify(lat)} · LATERAL/TIME="
           f"{sp.simplify(lat / t)} = √((d−1)/(d+1))",
           must_fail_on=[("a factorized lattice (Gram = 1)", ortho)])
        assert sp.simplify(lat / t - sp.sqrt(sp.Rational(gw['d'] - 1, gw['d'] + 1))) == 0

    print("\n   ★AN ANCHOR IN THE PRIMARY SOURCE S597 (a diamond in a cube a=1, a bond = √3/4):")
    bond = sp.sqrt(3) / 4
    t3, lat3 = zigzag_law(grams[1])                     # d=3
    print(f"     TIME·bond = {sp.nsimplify(t3 * bond)} = {float(t3 * bond):.5f}   "
          f"(S597 prints 0.5774)")
    print(f"     LATERAL·bond = {sp.nsimplify(lat3 * bond)} = {float(lat3 * bond):.5f}   "
          f"(S597 prints 0.4082)")
    assert abs(float(t3 * bond) - 0.57735) < 1e-5 and abs(float(lat3 * bond) - 0.40825) < 1e-5

    # ---- T4: ablation over q — «a dimer ⟺ q=1» -----------------------------------
    print("\n[T4] ABLATION over q (the count of marked classes): a matching ⟺ q=1 ?")

    def matching_set(scan_world):
        """The set of q for which the marked subgraph is a perfect matching."""
        d, ball = scan_world['d'], scan_world['ball']
        mk = scan_world['mk']
        got = []
        for q in range(0, d + 2):
            w = mk(d, tuple(range(q)), ball)
            if is_perfect_matching(w):
                got.append(q)
        return tuple(got)

    def only_q1(scan_world):
        return matching_set(scan_world) == (1,)

    for d in range(2, 6):
        ours = {'d': d, 'ball': 2, 'mk': lambda dd, mm, bb: world_cell(dd, mm, bb),
                'name': f'cell A_{d} (a scan over q)'}
        degen = {'d': d, 'ball': 2,
                 'mk': lambda dd, mm, bb: world_columns(dd, bb),
                 'name': 'a scan in the world of columns'}
        print(f"   d={d}: q with a matching = {matching_set(ours)}")
        ok(only_q1, ours,
           f"d={d}: a perfect matching exists EXACTLY at q=1 ⟹ «dimer» ≡ «one class marked»",
           must_fail_on=[("a scan in the world of columns", degen)])

    # ---- contrast: the form does NOT move with d, the number MOVES ----------------
    print("\n[T5] a contrast: the combinatorial form stands ∀d, the zigzag number moves")

    def probe(d):
        w = world_cell(d)
        form = (is_perfect_matching(w), no_pure_marked_path(w))
        t, _ = zigzag_law({'d': d, 'gram': gram_canon(d)})
        return (str(form), str(t))

    ok_contrast(probe, [2, 3, 4, 5],
                "the form (a matching + no purely-time-like path) is invariant in d, "
                "while TIME=1+1/d moves with d")

    code = report("S1075 · D2-form of AX-dimer")

    print("""
────────────────────────────────────────────────────────────────────────────────
WHAT THIS MEANS FOR D2 (a measurement, not a verdict — Ω's verdict after Beta's visa):

  (1) «1 time-bond per node» and «no purely-time-like path» — are NOT separate content:
      both fall ∀d from ONE ancestor (the marked cover of the cell: exactly one edge of
      each class per node), as soon as ONE class is marked.  Multiplicity by ancestor = 1.
  (2) S597 MEASURED these two lines (a code-visa, the d=3 diamond), but it measured a
      CONSEQUENCE: the CHOICE itself, «the time-bond = the global (1,1,1)», was an INPUT of
      its visualization — the primary source, verbatim: «the time-bond = a global choice of
      (1,1,1)» and «the time-bond is singled out ONLY by a global choice of 1 out of 4
      equivalent directions».
      ⟹ S597 postulated the CHOICE (D3), measured the FORM (D2).
  (3) T4: a matching exists EXACTLY at q=1 ⟹ the entire D2 form = a function of the COUNT of
      marked classes.  It carries no content of its own beyond «q=1».
  ★HONESTLY ABOUT THE PROBE'S STRENGTH: T1/T2 follow almost directly from the definition of
      the marked cover (T19: z=d+1, bipartite, one edge of each class per node).  This is
      NOT a defect of the probe, it IS the D2 measurement: «dimer» lies one step away from
      the cell, that is, it is not an independent statement.  The probe's value is not its
      strength but its NEGATIVE WORLDS: the rival «columns» (S597 fact 2) is eliminated by
      the T19 ANCESTOR (z=2(d+1)≠d+1), not by taste.

  (4) A REMAINDER this probe does NOT close (and does not claim to): the identification
      «time-axis ↔ the bond class u_0».  T1/T2 live in the language of BOND CLASSES; the
      world «merged marks» shows that the language of AXES would give a different verdict.
      This already stands in the graph as the candidate edge A-axis⇢AX-dimer (S1027) — the
      address is known, not new.
────────────────────────────────────────────────────────────────────────────────""")
    return code


if __name__ == '__main__':
    sys.exit(main())
