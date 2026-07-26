# -*- coding: utf-8 -*-
# DIM: symbolic-exact.  Numbers with a bracket [address · unit · type/operation].
"""S1066 B2 — THE AX-indef SIEGE, THE HEART: does {AX-dimer} ⟹ q≥1 (is the mark readable without a minus?).

ASSIGNMENT: exante `AX_INDEF_SIEGE_EXANTE.md`, task B2 + the author's word «go on B2».
B1 is delivered (6/0/0): Euclid is self-sufficient for T1–T4, there is no unconditional forcing.

★(2a) THE DEFINITION «THE MARK IS READABLE» — CARVED BEFORE THE COUNT, FROM AN ANCESTOR:
  T27 (the scale of clock-ness: q=0 — «no mark/time») and T34 (the uniqueness of the clock) demand
  a DISTINGUISHING CAPACITY, not a label.  Operationally, at the level of the arena:
      **the mark `u₀` IS READABLE ⟺ no ISOMETRY of the form carries `u₀` into another
      candidate-axis** (that is, a group-invariant exists that singles out `u₀`).
  If an isometry exists — «which axis is marked» is not a property of the structure, it is a LABEL.
  The detector = a function of the WORLD (the form + the mark + the candidates), not a ready-made boolean.

★WHAT I EXPECT AND WHAT I FEAR (BEFORE the count): I expect that in a definite arena the mark is NOT readable
(a democracy), while the minus reads it ⟹ an implication.  I fear PROVING TOO MUCH: readability
can be bought by something OTHER than the sign (a different norm · a degenerate form).  So I measure ALL worlds, and
name each PREMISE that does the work as a separate line.

★THE FORK, named by me before this step (CODEX: «indefinite ⟺ q≥1 = a Witt tautology»):
the question stands NOT against the postulate, but against the DEFINITION.  So I formulate the verdict as
an IMPLICATION WITH EXPLICIT PREMISES, not as «q≥1 was derived from nothing».

(2c) A SURFACE GUARD: the column/dimer — an object of the LATTICE level; the arena — a quadratic
space.  I stay ON THE ARENA (where the consumer T5 lives) and I do NOT build a bridge; the precedents
S1027 / S1063-T2 by citation — I do not make a third corpse.

FENCE: Layer-1.  The physics-vocabulary classes named in the project's fence do not enter the code.
Isometries are BUILT explicitly (constructively), the Witt theorem is NOT imported as a formula.
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


def norm(G, v):
    return sp.simplify((v.T * G * v)[0, 0])


def is_isometry(G, R):
    return sp.simplify(R.T * G * R - G) == sp.zeros(G.rows, G.cols)


def build_isometry(G, u, v):
    """CONSTRUCTIVELY: a reflection that carries `u` to `v` (the classical construction via the
    difference vector).  Returns R or None.  ★I do NOT import the Witt theorem as a FORMULA —
    I build the matrix and CHECK `RᵀGR = G` and `Ru = v` exactly."""
    w = sp.simplify(u - v)
    nw = norm(G, w)
    if sp.simplify(nw) == 0:
        return None
    R = sp.eye(G.rows) - (2 / nw) * (w * (w.T * G))
    R = sp.simplify(R)
    if is_isometry(G, R) and sp.simplify(R * u - v) == sp.zeros(G.rows, 1):
        return R
    return None


def mark_is_label(world):
    """A WORLD = (G, the index of the mark, candidate-axes) ⟹ whether the mark is a LABEL, that is, whether an
    isometry EXISTS that carries the marked axis into some other candidate-axis.
    A LABEL = NOT readable.  Readable = the negation of this."""
    G, i0, cand = world
    n = G.rows
    e = [sp.Matrix([1 if r == c else 0 for r in range(n)]) for c in range(n)]
    u = e[i0]
    for j in cand:
        if j == i0:
            continue
        if build_isometry(G, u, e[j]) is not None:
            return True
    return False


def mark_is_readable(world):
    """The negation of the previous — precisely the definition of (2a)."""
    return not mark_is_label(world)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1066_2_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf); sys.stdout = tee
    reset()

    print("=" * 80)
    print("S1066 B2 — THE HEART OF THE SIEGE: is the mark READABLE without a minus")
    print("=" * 80)
    print()
    print("(2a) THE DEFINITION, carved BEFORE the count [T27/T34 — a distinguishing capacity, not a label]:")
    print("     the mark u₀ IS READABLE ⟺ NO isometry of the form carries u₀ into another candidate-axis.")
    print()

    n = 4
    cand = list(range(n))
    W_def = (sp.eye(n), 0, cand)                                    # definite, equal norms
    W_indef = (sp.diag(1, 1, 1, -1), 3, cand)                       # q=1: the mark = the minus-axis
    W_scaled = (sp.diag(2, 1, 1, 1), 0, cand)                       # definite, BUT a different norm
    W_degen = (sp.diag(1, 1, 1, 0), 3, cand)                        # DEGENERATE, q=0
    names = {"definite, equal norms (a democracy)": W_def,
             "indefinite q=1 (the mark = the minus-axis)": W_indef,
             "definite, the marked axis has NORM 2": W_scaled,
             "DEGENERATE (a radical), q=0": W_degen}
    for nm, w in names.items():
        G, i0, _ = w
        e0 = sp.Matrix([1 if r == i0 else 0 for r in range(n)])
        print(f"     [{nm}] the norm of the mark = {norm(G, e0)} · the norms of the others = "
              f"{[norm(G, sp.Matrix([1 if r == j else 0 for r in range(n)])) for j in cand if j != i0]}"
              f" ⟹ {'READABLE' if mark_is_readable(w) else 'A LABEL'}"
              f"  [2b · dimensionless · the existence of an isometry u₀↦another axis]")
    print()

    # ═══════ (2b) THE CORE: IS IT READABLE IN EUCLID ═══════
    print("(2b) THE CORE — computed IN THE FIELD (isometries are BUILT, the Witt theorem is not imported)")

    ok(mark_is_label, W_def,
       "★★IN A DEFINITE ARENA WITH EQUAL NORMS THE MARK IS A LABEL: an isometry that carries "
       "the marked axis into any other, EXISTS and is built explicitly (RᵀGR=G, Ru₀=u_j checked) "
       "[2b · dimensionless · the existence of an isometry]",
       must_fail_on=[("★indefinite (q=1): the norm of the mark −1 ≠ +1 ⟹ there is no isometry", W_indef),
                     ("definite, but the norm of the mark = 2 ⟹ there is no isometry", W_scaled)])
    R = build_isometry(sp.eye(n), sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 1, 0, 0]))
    print(f"     an explicit isometry (Euclid, e₀↦e₁):\n{sp.pretty(R)}")
    print("     ⟹ the democracy of the axes in a definite arena is NOT broken by the mark: the mark has no")
    print("       invariant carrier, that is, it IS A LABEL, not a property of the structure.")

    ok(mark_is_readable, W_indef,
       "★★★THE MINUS MAKES THE MARK READABLE: an isometry preserves the form ⟹ norms of −1 and +1 are not "
       "carried into one another ⟹ the marked axis is singled out INVARIANTLY  [2b · dimensionless · "
       "the non-existence of an isometry between different norms]",
       must_fail_on=[("definite with equal norms — there the mark IS a label", W_def)])
    print()

    # ═══════ ★WHICH PREMISES DO THE WORK (so as not to prove too much) ═══════
    print("★THE PREMISES THAT DO THE WORK — each with its own alternative world")

    ok(mark_is_readable, W_scaled,
       "★PREMISE-1 «EQUAL NORMS» (a democracy, the AX-alphabet) DOES THE WORK: in a definite "
       "arena the mark CAN be made readable WITHOUT a minus — with a different norm  [premise-1 · "
       "dimensionless · readability at q=0, norm 2]",
       must_fail_on=[("the same definite arena with EQUAL norms — the mark IS a label", W_def)])
    print("     ⟹ so the implication does NOT follow from definiteness alone: without a democracy")
    print("       (of equal norms) readability is bought by scale.  The democracy is given by the ROOT")
    print("       AX-alphabet {a bare set of d+1 · S_{d+1}} — and it is precisely this that works here,")
    print("       while the T19-cell (unit axes, a Gram −1/d) forbids rescaling.")

    ok(mark_is_readable, W_degen,
       "★★PREMISE-2 «NON-DEGENERACY» DOES THE WORK: in a DEGENERATE form (a radical) the mark "
       "is readable at q=0 — a zero norm is invariant  [premise-2 · dimensionless · "
       "readability at q=0 with a radical]",
       must_fail_on=[("a non-degenerate definite form with equal norms — the mark IS a label", W_def)])
    print("     ⟹ ★A THIRD WORLD that the exante did not name: a degenerate form reads the mark")
    print("       WITHOUT a minus.  So the implication needs NON-DEGENERACY as an EXPLICIT premise;")
    print("       it cannot be taken silently — that would be a hidden axiom.")
    print()

    # ═══════ A CHECK OF THE MACHINERY ITSELF ═══════
    print("A CHECK OF THE MACHINERY — are my «isometries» genuinely isometries")

    def built_map_is_isometry(world):
        """A WORLD = (G, u, v) ⟹ whether the BUILT matrix genuinely preserves the form AND carries
        u to v.  Without this line the whole of (2b) would rest on my word."""
        G, u, v = world
        R = build_isometry(G, u, v)
        return R is not None and is_isometry(G, R) and sp.simplify(R * u - v) == sp.zeros(G.rows, 1)

    e0 = sp.Matrix([1, 0, 0, 0]); e1 = sp.Matrix([0, 1, 0, 0]); e3 = sp.Matrix([0, 0, 0, 1])
    ok(built_map_is_isometry, (sp.eye(n), e0, e1),
       "★the built matrix IS an isometry and genuinely carries e₀↦e₁ (checked exactly) "
       "[machinery · dimensionless · RᵀGR−G and Ru−v]",
       must_fail_on=[("★A DIFFERENT LAW: an attempt to carry DIFFERENT-NORM vectors (an indefinite "
                      "arena, e₃↦e₀) — the construction does not give an isometry", (sp.diag(1, 1, 1, -1), e3, e0)),
                     ("the same attempt in a degenerate form", (sp.diag(1, 1, 1, 0), e3, e0))])
    print()

    code = report("S1066 B2 — the implication {a mark} ⟹ q≥1")
    print()
    print("=" * 80)
    print("★RAW OUTPUTS OF B2 (no verdict rendered — a B-visa, then the Ω court)")
    print("  (2a) the definition carved BEFORE the count: readable ⟺ there is no isometry u₀↦another axis.")
    print("  (2b) THE CORE:")
    print("     · a definite arena + EQUAL norms ⟹ the mark = A LABEL (an isometry is built explicitly);")
    print("     · indefinite (q=1) ⟹ the mark IS READABLE (norms −1 and +1 are not carried into one another).")
    print("  ⟹ ★AN IMPLICATION EXISTS, BUT WITH TWO EXPLICIT PREMISES:")
    print("       {EQUAL NORMS (a democracy, the AX-alphabet) ∧ NON-DEGENERACY ∧ the mark is READABLE}")
    print("       ⟹ q≥1.")
    print("     Both premises are MEASURED as doing work — for each, a world is exhibited where it")
    print("     fails and readability is bought WITHOUT a minus (norm 2 · a degenerate form).")
    print("  ★A THIRD WORLD (the exante did not name it): a DEGENERATE form reads the mark at q=0 ⟹")
    print("     non-degeneracy must stand as an EXPLICIT premise, otherwise it is a hidden axiom.")
    print("  (2c) THE SURFACE GUARD IS OBSERVED: computed ON THE ARENA (where the consumer T5 lives);")
    print("     the column/dimer of the lattice level was NOT identified with it — there is no bridge.")
    print("=" * 80)
    sys.stdout = tee.real
    logf.close()
    return code


if __name__ == "__main__":
    sys.exit(main())
