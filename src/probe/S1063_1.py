# -*- coding: utf-8 -*-
# DIM: symbolic-exact.  Numbers with a bracket [address · unit · type/operation].
"""S1063 T1 — THE SEAM (d,1), step 1: BUILD N3 and read its signature.

ASSIGNMENT: project ruling, exante `hub/prime/SEAM_D1_SIGNATURE_EXANTE.md` (carved BEFORE the count), T1.
The order of the assignment: T1 → STOP.  No verdict is rendered.

ANCESTORS (addresses, not paraphrases):
  · S1045 Component 1 — `det(ϸ·I − H − m₀σ_z) = ϸ² − m₀² − |f|²`, `H=[[0,f],[f̄,0]]`;
  · S1045 Component 2 — the node-set `f(k)=0` has dimension **d−2** (2 real conditions);
  · S1047 — at d=3 a generic node gives the Hessian |f|² rank-2 (+1 flat along the line);
    a special point (∇Re f=0) — rank-1 (2 flat);
  · T19 — the cell: d+1 unit axes, a pairwise Gram −1/d, **the sum of the axes = 0**;
    the node at the symmetric point: the phases = the full (d+1)-th roots of unity (a Gauss-sum = 0);
  · T33 — the transverse cone, `v² = ½·trM` is TRACE-based (and it is precisely this that shifted under
    relabeling — the errata: S1054/S1055/S1057/S1058);
  · T39 — `m₀` = a split (a chirality-even `σ_z` insertion);
  · S1061 — `U†U` carries `m₀` in the NORM (a polarity), not in the phase.

★WHAT I EXPECT AND WHAT I FEAR (carved BEFORE the count, in the probe's body): the exante's wager is named —
«|f|² collects d transverses ⟹ a (1,d)-candidate».  But ancestor S1045 says the node-set
has dimension d−2, and S1047 — that the Hessian there is rank-2.  These two lines are INCOMPATIBLE with
the wager at d≥3, and I fear I will recognize this only after writing «(1,d)» by hand.  So I measure
the INERTIA (a sign-count), rather than describe it.

FENCE: Layer-1.  The physics-vocabulary classes named in the project's fence do not enter the code
(K-stone `KAMIN.md`: the signature is measured, the GROUP is not).  Homonyms with a prefix: `minus_arith`
(T32, a participant) ⊥ `minus_geom`; «tact» is not used here at all.  `ϸ` = a spectral
parameter (there is no bare «E»).  `m₀` = a representative of the regime (J-0486), not a canon-number.
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


# ═══════ (1a) THE FORM — EXACTLY FROM ANCESTORS, WITHOUT NEW TERMS ═══════
def simplex_bonds(d):
    """The T19-cell: d+1 unit axes, a pairwise Gram −1/d, the SUM = 0.
    A direct construction (no expansions): δ_a = e_a − 𝟙/(d+1), normalized.
    Then δ_a·δ_b = (δ_ab − 1/(d+1))/(d/(d+1)) ⟹ −1/d at a≠b, and Σδ_a = 0 exactly.
    The axes lie in the d-dimensional plane ⟨𝟙⟩^⊥ ⊂ ℝ^{d+1} — writing in (d+1) coordinates
    changes nothing: the Gram and the sum are invariant to the choice of basis of the plane."""
    n = d + 1
    rows = []
    for a in range(n):
        v = sp.Matrix([sp.Rational(1, 1) if i == a else sp.Integer(0) for i in range(n)]) \
            - sp.Matrix([sp.Rational(1, n)] * n)
        rows.append((v / sp.sqrt(sp.Rational(d, n))).T)
    return sp.simplify(sp.Matrix.vstack(*rows))


def gram_of(rows):
    """The Gram of the normalized rows — an invariant that does not depend on the scale of the writing."""
    n = rows.rows
    out = sp.zeros(n, n)
    for a in range(n):
        na = sp.sqrt(sum(rows[a, i]**2 for i in range(rows.cols)))
        for b in range(n):
            nb = sp.sqrt(sum(rows[b, i]**2 for i in range(rows.cols)))
            out[a, b] = sp.simplify(sum(rows[a, i] * rows[b, i] for i in range(rows.cols)) / (na * nb))
    return out


def node_phases(d):
    """The EXACT phases of the node θ_a = k·δ_a  (Σθ_a = 0 because Σδ_a = 0 [T19]):
    · d EVEN (d+1 odd) ⟹ the full (d+1)-th roots of unity [T19: a Gauss-sum = 0],
      representatives chosen so that Σθ = 0 exactly;
    · d ODD (d+1 even) ⟹ PAIRINGS (α_j, α_j+π): each pair cancels itself in f
      ∀α_j ⟹ the free α = this IS the node-SET of S1045 (at d=3 — a line).
      The last α is shifted by −S/2, so that Σθ = 0 EXACTLY (f does not depend on this)."""
    if d % 2 == 0:
        raw = [2 * sp.pi * a / (d + 1) for a in range(d + 1)]
        return [sp.nsimplify(th if th <= sp.pi else th - 2 * sp.pi) for th in raw]
    # α on a SHARED grid of π/12: a generic choice (not 0 and not π/2 — there the rank drops,
    # see the special point below), and a common denominator keeps the algebraic degree
    # low.  The grid is a writing convenience; the node remains a node ∀α by construction.
    alphas = [sp.pi * sp.Rational(2 * j + 1, 12) for j in range((d + 1) // 2)]
    ph = []
    for al in alphas:
        ph += [al, al + sp.pi]
    S = sum(ph)
    alphas[-1] = sp.nsimplify(alphas[-1] - S / 2)        # the pair remains a pair: f = 0 ∀α
    ph = []
    for al in alphas:
        ph += [sp.nsimplify(al), sp.nsimplify(al + sp.pi)]
    return ph


def f_at(phases):
    """f = Σ_a e^{iθ_a}  [S1045/b1_resolvent_chain_probe: f(k)=Σ_δ e^{ik·δ}].
    Computed via (Σcos, Σsin) — otherwise sympy leaves exp/(-1)^(p/q) unresolved
    and «f = 0» is not recognized (this is a defect of the WRITING, not of the object)."""
    re_ = sp.simplify(sp.expand_trig(sum(sp.cos(t) for t in phases)))
    im_ = sp.simplify(sp.expand_trig(sum(sp.sin(t) for t in phases)))
    return sp.nsimplify(sp.radsimp(re_)) + sp.I * sp.nsimplify(sp.radsimp(im_))


def hess_absf2(phases):
    """The Hessian of |f|² in phase-coordinates AT THE NODE: at f=0 the term 2Re(f̄·∂²f) vanishes exactly,
    leaving M_ab = 2·Re(∂_a f · conj(∂_b f)) = 2·cos(θ_a − θ_b)
    [S1047 · dimensionless · the second derivative of |f|² at the node]."""
    n = len(phases)
    return sp.Matrix(n, n, lambda a, b: sp.simplify(2 * sp.cos(phases[a] - phases[b])))


def hyperplane_basis(n):
    """A basis of the plane Σθ_a = 0 — precisely where k lives (because Σδ_a = 0 [T19]).
    The transition (k)→(θ) is LINEAR and injective ⟹ the inertia will not shift (Sylvester)."""
    return sp.Matrix.hstack(*[sp.Matrix([1 if i == 0 else (-1 if i == j else 0)
                                         for i in range(n)]) for j in range(1, n)])


def _sign_changes(coeffs):
    s = [c for c in coeffs if sp.N(sp.Abs(c), 40) > sp.Float("1e-30")]
    return sum(1 for a, b in zip(s, s[1:]) if sp.N(a, 40) * sp.N(b, 40) < 0)


def inertia(Msym):
    """THE SIGN-COUNT (inertia) of a symmetric matrix: (n₊, n₋, n₀).

    ★WHY NOT VIA EIGENVALUES: sympy computes the radical eigenvalues of a 6×6 for
    hours, and this is a cost WITHOUT a gain — the inertia does not need them.  I take two EXACT
    facts: (1) n₀ = n − rank (a symbolic rank) · (2) for a symmetric matrix ALL
    roots of the characteristic polynomial are real ⟹ Descartes' rule of signs gives the
    EXACT (not an estimated) count of positive and negative roots.
    The sign of the coefficient itself (an algebraic number) is read at 40 significant digits —
    zero coefficients are filtered out BEFORE that, SYMBOLICALLY."""
    A = sp.Matrix(Msym).applyfunc(lambda e: sp.radsimp(sp.simplify(e)))
    n, lam = A.rows, sp.Symbol('lam')
    zer = n - A.rank()
    coeffs = [sp.simplify(sp.expand(c)) for c in sp.Poly(A.charpoly(lam).as_expr(), lam).all_coeffs()]
    while coeffs and sp.simplify(coeffs[-1]) == 0:       # a root of 0 — already counted in zer
        coeffs.pop()
    deg = len(coeffs) - 1
    pos = _sign_changes(coeffs)
    neg = _sign_changes([c * (-1)**(deg - i) for i, c in enumerate(coeffs)])
    return (pos, neg, zer)


def seam_form(d, phases):
    """★(1a) THE QUADRATIC PART OF THE SEAM on (ϸ, k₁…k_d), assembled EXACTLY from ancestors:
        Q = ϸ² − |f(k)|²      [S1045 · det = ϸ² − m₀² − |f|²]
    where |f|² near the node = k_⊥ᵀ·M·k_⊥ [S1047].  `m₀²` — A TERM WITHOUT A COORDINATE (a shift
    of the quadric), it does NOT enter the quadratic part; this is named, not forgotten.
    No new terms were introduced (K2)."""
    M = hess_absf2(phases)
    Bh = hyperplane_basis(len(phases))
    Mk = sp.simplify(Bh.T * M * Bh)                      # restricted to the plane Σθ=0
    Q = sp.zeros(d + 1, d + 1)
    Q[0, 0] = 1                                          # +ϸ²  [S1045 · dimensionless · char. poly]
    for i in range(d):
        for j in range(d):
            Q[i + 1, j + 1] = -Mk[i, j]                  # −|f|²  [S1047 · dimensionless · the Hessian]
    return sp.simplify(Q), Mk


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1063_1_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf); sys.stdout = tee
    reset()

    print("=" * 80)
    print("S1063 T1 — THE SEAM (d,1): build N3 and COMPUTE its signature")
    print("=" * 80)
    print()

    # ═══════ THE ANCHOR: THE CELL — T19, AND IT MATCHES THE ANCESTOR'S CELL ═══════
    print("THE ANCHOR — the cell from T19 (d+1 axes, Gram −1/d, Σδ=0), checked against ancestor S1045")

    def gram_is_T19(rows_d):
        """A WORLD = (a set of axes, d) ⟹ whether the pairwise Gram of the normalized axes = −1/d [T19]."""
        rows, d = rows_d
        G = gram_of(rows)
        return all(sp.simplify(G[a, b] + sp.Rational(1, d)) == 0
                   for a in range(G.rows) for b in range(G.cols) if a != b)

    anc_diamond = sp.Matrix([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])   # S1045 bonds(3)
    anc_honey = sp.Matrix([[1, 0], [-sp.Rational(1, 2), sp.sqrt(3) / 2],
                           [-sp.Rational(1, 2), -sp.sqrt(3) / 2]])                # S1045 bonds(2)
    ok(gram_is_T19, (anc_diamond, 3),
       "★the ANCESTOR's cell (S1045 bonds(3), a tetrahedron) has a Gram −1/3 ⟹ my T19-construction "
       "stands on the SAME cell, not a new one  [anchor · dimensionless · a pairwise Gram]",
       must_fail_on=[("orthogonal axes (a Gram of 0 — a different cell)",
                      (sp.eye(3).col_join(-sp.eye(3)[0, :]), 3)),
                     ("a square cell d=2 (a Gram of 0)", (sp.eye(2), 2))])
    for d in (2, 3, 4):
        B = simplex_bonds(d)
        G = gram_of(B)
        print(f"     d={d}: axes {B.rows}, the off-diagonal of the Gram = {sp.simplify(G[0, 1])} "
              f"(T19: −1/{d}), Σδ = {sp.simplify(sp.Matrix([sum(B[a, i] for a in range(B.rows)) for i in range(B.cols)]).T)}"
              f"  [T19 · dimensionless · the Gram and the sum of the axes]")
    print(f"     checked against the ancestor: Gram(S1045 bonds(2)) off-diag = {sp.simplify(gram_of(anc_honey)[0, 1])}")
    print()

    # ═══════ (1a) THE NODE AND THE FORM ═══════
    print("(1a) THE NODE EXACTLY, AND THE FORM FROM ANCESTORS (no new term — K2)")

    def is_node(phases):
        """A WORLD = a set of phases ⟹ whether f = Σ e^{iθ} = 0 EXACTLY (and whether Σθ = 0, that is, the phases are legitimate)."""
        return sp.simplify(f_at(phases)) == 0 and sp.simplify(sum(phases)) == 0

    for d in (2, 3, 4, 5):
        ph = node_phases(d)
        ok(is_node, ph,
           f"★d={d}: the node is EXACT — f = Σ e^{{iθ}} = 0 at Σθ = 0  [1a · dimensionless · "
           f"the value of f at the node]",
           must_fail_on=[(f"a shifted phase (θ₀ → θ₀+π/6) — f ≠ 0",
                          [ph[0] + sp.pi / 6] + [t - sp.pi / (6 * d) for t in ph[1:]])])
        print(f"     d={d}: θ = {ph},  f = {f_at(ph)}  [1a · dimensionless · the phases of the node]")
    print()

    # ═══════ ★THE MAIN MEASUREMENT: THE INERTIA OF THE SEAM ═══════
    print("★THE MAIN MEASUREMENT — THE INERTIA (n₊, n₋, n₀) of the quadratic part of the seam on (ϸ, k₁…k_d)")
    print("  the exante's wager: a (1,d)-candidate.  Measured by a sign-count, not by description.")
    table = {}
    for d in (2, 3, 4, 5):
        ph = node_phases(d)
        Q, Mk = seam_form(d, ph)
        table[d] = (inertia(Q), inertia(Mk), sp.Matrix(Mk).rank())
        print(f"     d={d}: the inertia Q(ϸ,k) = {table[d][0]}  ·  the rank of the Hessian |f|² = {table[d][2]}"
              f"  ·  the inertia of the Hessian itself = {table[d][1]}"
              f"  [T1 · dimensionless · a sign-count of eigenvalues]")

    def inertia_is_one_two(d):
        """A WORLD = d ⟹ whether the inertia of the seam is (1, 2, d−2): one plus from ϸ², TWO minuses from
        transverses, the rest — EXACT zeros along the node-set."""
        ph = node_phases(d)
        Q, _ = seam_form(d, ph)
        return inertia(Q) == (1, 2, d - 2)

    def inertia_is_one_d(d):
        """The exante's wager detector: whether the inertia is (1, d, 0)."""
        ph = node_phases(d)
        Q, _ = seam_form(d, ph)
        return inertia(Q) == (1, d, 0)

    def inertia_law_holds(world):
        """A WORLD = (d, A HESSIAN RULE) ⟹ whether the inertia of the seam = (1, 2, d−2).
        ★The first form of this assert was EMPTY: I gave d=2 and d=5 as negative worlds,
        and the law (1,2,d−2) holds there too — that is, the «negative worlds» were the same
        law at other d.  A negative world must break the LAW, not change d."""
        d, hess_rule = world
        ph = node_phases(d)
        M = hess_rule(ph)
        Bh = hyperplane_basis(len(ph))
        Mk = sp.simplify(Bh.T * M * Bh)
        Q = sp.zeros(d + 1, d + 1)
        Q[0, 0] = 1
        Q[1:, 1:] = -Mk
        return inertia(Q) == (1, 2, d - 2)

    def isotropic_hess(phases):
        """Negative world 1 — the world the exante's wager IMAGINED: if |f|² bent
        along ALL d directions (an isotropic cone), the Hessian would be full rank ⟹ (1,d,0)."""
        n = len(phases)
        return 2 * sp.eye(n)

    ok(inertia_law_holds, (3, hess_absf2),
       "★★THE INERTIA OF THE SEAM = (1, 2, d−2) at d=3: ONE plus (ϸ²) + TWO minuses (transverses) + "
       "ONE exact zero along the node-LINE  [T1 · dimensionless · the inertia of Q at d=3]",
       must_fail_on=[("★the exante's wager's world: an isotropic full-rank Hessian ⟹ (1,d,0)",
                      (3, isotropic_hess)),
                     ("★the special point (the rank drops 2→1) ⟹ (1,1,d−1) — the law does NOT hold there",
                      (3, lambda ph: hess_absf2([sp.Integer(0), sp.pi, -sp.pi, sp.Integer(0)])))])

    ok(inertia_is_one_d, 2,
       "★the exante's wager «(1,d)» holds EXACTLY at d=2 (because there d−2=0 ⟹ (1,2,0)=(1,d,0)) "
       "— that is, it is not wrong, but SPECIAL  [T1 · dimensionless · the inertia of Q at d=2]",
       must_fail_on=[("d=3 — (1,2,1) ≠ (1,3,0), the exante's wager FAILS", 3),
                     ("d=4 — (1,2,2) ≠ (1,4,0)", 4),
                     ("d=5 — (1,2,3) ≠ (1,5,0)", 5)])
    print("     ⟹ ★★THE EXANTE'S WAGER FAILS AT d≥3, AND I SAY THIS DIRECTLY: «|f|² collects d")
    print("       transverses» is true only at d=2.  At d≥3 the node-set has dimension d−2")
    print("       [S1045] and along it |f|≡0 EXACTLY ⟹ (d−2) directions give EXACT zeros,")
    print("       not minuses.  The live (non-degenerate) block of the seam = (1,2) ∀d.")
    print()

    # ═══════ THE SPECIAL POINT: A DROP IN RANK (ancestor S1047) ═══════
    print("THE SPECIAL POINT — ancestor S1047 (a drop in rank): measured, not cited")
    spec3 = [sp.Integer(0), sp.pi, -sp.pi, sp.Integer(0)]        # α=0 on the line d=3

    def rank_drops(phases):
        """A WORLD = node phases ⟹ whether the rank of the Hessian |f|² DROPS below 2 (all directions
        collinear: θ_a ∈ {φ, φ+π})."""
        M = hess_absf2(phases)
        Bh = hyperplane_basis(len(phases))
        return sp.Matrix(sp.simplify(Bh.T * M * Bh)).rank() < 2

    ok(rank_drops, spec3,
       "★the special point EXISTS at d=3 too (α=0 on the node-line): the rank of the Hessian drops 2→1 ⟹ there "
       "the inertia of the seam is (1,1,2), not (1,2,1)  [S1047 · dimensionless · the rank at the special point]",
       must_fail_on=[("a generic point of the same line (α=π/5) — rank 2", node_phases(3))])
    Qs, Mks = seam_form(3, spec3)
    print(f"     the special point d=3: θ = {spec3}, f = {f_at(spec3)}, rank = {sp.Matrix(Mks).rank()}, "
          f"the inertia of Q = {inertia(Qs)}  [S1047 · dimensionless · the inertia at the special point]")
    print("     ⟹ the signature of the seam is NOT constant over the node-set: (1,2,d−2) generically ⊥ (1,1,d−1)")
    print("       at the special points.  This is ancestor S1047, named, not a new witness (multiplicity 1).")
    print()

    # ═══════ (1b) BINARY: FORCED BY STRUCTURE OR BY THE WRITING ═══════
    print("(1b) BINARY — is the signature forced by STRUCTURE, or is it a property of the WRITING?")
    ps, m_s = sp.symbols('ps m', real=True, positive=True)
    fr, fi = sp.symbols('f_re f_im', real=True)
    f_sym = fr + sp.I * fi
    H_hop = sp.Matrix([[0, f_sym], [sp.conjugate(f_sym), 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])
    sx = sp.Matrix([[0, 1], [1, 0]])

    def three_independent_squares(mass_channel):
        """A WORLD = A SPLITTER CHANNEL (a matrix X in the insertion m·X) ⟹ whether the characteristic polynomial
        is a PURE sum of squares over (ϸ, m, f_re, f_im) — that is, NO cross
        term.  ★The first form was EMPTY: I was checking «is there a term m²», and there is one even in
        absorbing channels (|f+m|² also has m²).  It is precisely the CROSS term that distinguishes:
        it is exactly what means the splitter can be removed by a shift of ϸ (channel I) or by redefining
        f (channel σ_x), that is, there is NO separate third square."""
        Hm = H_hop + m_s * mass_channel
        c = sp.expand((ps * sp.eye(2) - Hm).det())
        vars_ = [ps, m_s, fr, fi]
        for i in range(len(vars_)):
            for j in range(i + 1, len(vars_)):
                if sp.simplify(sp.diff(c, vars_[i], vars_[j])) != 0:
                    return False
        return sp.simplify(sp.diff(c, m_s, 2)) != 0

    ok(three_independent_squares, sz,
       "★★THE THREE SQUARES ARE INDEPENDENT PRECISELY BECAUSE THE SPLITTER ANTICOMMUTES WITH THE HOP: {H, σ_z} = 0 ⟹ "
       "(H+m₀σ_z)² = (m₀²+|f|²)·I ⟹ det = ϸ² − m₀² − |f|²  [1b · dimensionless · "
       "the characteristic form by splitter-channel]",
       must_fail_on=[("the splitter ∝ I (commutes) — absorbed by a SHIFT of ϸ, no separate square",
                      sp.eye(2)),
                     ("the splitter ∝ σ_x (the same channel as the hop) — absorbed INTO f, "
                      "no separate square", sx)])
    for name, X in (("σ_z [T39]", sz), ("I", sp.eye(2)), ("σ_x", sx)):
        Hm = H_hop + m_s * X
        print(f"     channel {name}: det(ϸ−H−mX) = {sp.factor(sp.expand((ps * sp.eye(2) - Hm).det()))}"
              f"  [1b · dimensionless · the characteristic polynomial]")
    print("     ⟹ ★BINARY: THE SIGNS ARE FORCED BY STRUCTURE, not by the writing.  One plus against")
    print("       the minuses — is not a choice of convention, but the ANTICOMMUTATION of the channels {hop, splitter};")
    print("       in the two legitimate alternative channels the third square VANISHES.")
    print()

    # ═══════ ★SYLVESTER: WHY THIS IS A MEASUREMENT, NOT A WRITING (the T33 erratum lesson) ═══════
    print("★SYLVESTER — does the inertia withstand what broke the T33 number (relabeling)")

    def probe_recoord(P_gen):
        """ONE body, ONE object: the same form d=3 in RELABELED coordinates
        k → P·k yields (what stands still) the inertia  ⊥  (what moves) the trace-number ½·trM,
        that is, exactly the quantity that already shifted once under relabeling
        [the T33 erratum: S1054/S1057]."""
        ph = node_phases(3)
        Q, Mk = seam_form(3, ph)
        P = sp.Matrix(3, 3, lambda i, j: P_gen**(i * 3 + j) if i != j else 1)
        Qp = sp.zeros(4, 4)
        Qp[0, 0] = Q[0, 0]
        Qp[1:, 1:] = sp.simplify(P.T * Q[1:, 1:] * P)
        trace_num = sp.nsimplify(sp.simplify(sp.Rational(1, 2) * (P.T * Mk * P).trace()))
        return inertia(Qp), trace_num

    ok_contrast(probe_recoord, [sp.Rational(1, 3), sp.Rational(1, 2), sp.Integer(2)],
                "★★THE INERTIA HOLDS under a relabeling of coordinates, whereas the trace-number "
                "½·trM MOVES — the same machinery moves the control  [T1 · dimensionless · "
                "(the inertia) ⊥ (½·trM) under k→P·k]")
    print("     ⟹ ★THIS IS PRECISELY THE REASON T1 MEASURES THE INERTIA: a quantity that shifts under")
    print("       a free relabeling, does not measure the object (the class S1057/S1059, the T33")
    print("       erratum).  The inertia is a Sylvester invariant, and it is the only legitimate ruler here.")
    print()

    # ═══════ (1c) THE FENCE OF SQUARES: PROVENANCE ═══════
    print("(1c) THE FENCE OF SQUARES — where does ϸ² come from: from `U†U` (S1061) or from the dispersion?")

    def square_is_spectral_free(op):
        """A WORLD = an operator ⟹ whether its square is scalar AND whether that scalar is a POLYNOMIAL in
        (m, f) WITHOUT the spectral parameter ϸ.
        ★The first form was EMPTY, and the mistake was SUBSTANTIVE: I wrote that the square
        of the one-tact `U` is not scalar — but it IS scalar always (an off-diagonal
        2×2: [[0,a],[b,0]]² = ab·I).  That is, scalarity distinguishes nothing.
        What distinguishes is the PROVENANCE: in `H+m₀σ_z` the scalar = (m₀²+|f|²) — a polynomial, no ϸ;
        in `U` the scalar = |f|²/(ϸ²−m₀²) — rational, with ϸ and poles.  The squares
        of the signature are taken only from the former."""
        sq = sp.simplify(sp.expand(op * op))
        if sp.simplify(sq - sq[0, 0] * sp.eye(2)) != sp.zeros(2, 2):
            return False
        c = sp.simplify(sq[0, 0])
        return ps not in c.free_symbols and sp.denom(sp.together(c)) == 1

    U_tact = sp.simplify(H_hop * sp.diag(1 / (ps - m_s), 1 / (ps + m_s)))   # S1061 U = H·G₀
    ok(square_is_spectral_free, H_hop + m_s * sz,
       "★★THE PROVENANCE OF ALL THREE SQUARES IS ONE: `(H+m₀σ_z)² = (m₀²+|f|²)·I` — the scalar is "
       "POLYNOMIAL and WITHOUT ϸ ⟹ the squares come from the operator itself, and ϸ² — from the "
       "characteristic polynomial [S1045]  [1c · dimensionless · the scalar of the square]",
       must_fail_on=[("the one-tact U = H·G₀ [S1061] — its square is ALSO scalar, but "
                      "the scalar = |f|²/(ϸ²−m₀²): rational, with ϸ and poles ⟹ a different "
                      "object, squares are not taken from it", U_tact)])
    UU = sp.simplify(U_tact.H * U_tact)
    print(f"     U†U [S1061] = diag({sp.simplify(UU[0, 0])}, {sp.simplify(UU[1, 1])})")
    print(f"     the scalar of the square: (H+m₀σ_z)² = {sp.simplify((H_hop + m_s * sz)**2)[0, 0]}"
          f"  ⊥  U² = {sp.simplify(U_tact * U_tact)[0, 0]}"
          "  [1c · dimensionless · the scalar of the square of two objects]")
    print("     ⟹ in `U†U` the splitter sits in the NORM (a polarity of the blades, S1061) — this is not a square")
    print("       of the signature.  The square `ϸ²` comes from the CHARACTERISTIC POLYNOMIAL")
    print("       (a spectral parameter), `m₀²` and `|f|²` — from the Clifford square.")
    print("       ★An exponent-smuggling does not sneak in here: there is no `exp` in the form at all.")
    print()

    code = report("S1063 T1 — the seam (d,1)")
    print()
    print("=" * 80)
    print("★RAW OUTPUTS OF T1 (no verdict rendered — the assignment: T1 → STOP)")
    print("  a table of the seam inertias Q(ϸ, k₁…k_d) at a generic point of the node-set:")
    for d in (2, 3, 4, 5):
        n_p, n_m, n_0 = table[d][0]
        print(f"     d={d}: (n₊, n₋, n₀) = ({n_p}, {n_m}, {n_0})   the rank of the Hessian |f|² = {table[d][2]}"
              f"   [T1 · dimensionless · the inertia]")
    print("  the special point d=3 (α=0): (1, 1, 2), rank 1  [S1047 · dimensionless · the inertia]")
    print("  (1a) the form is assembled from ancestors, no new terms (K2 did not fire).")
    print("  (1b) ★FORCED BY STRUCTURE: three squares exist ⟺ {hop, splitter} = 0;")
    print("       in the channels I and σ_x the third square vanishes (both worlds were computed).")
    print("  (1c) provenance: all squares from the Clifford (H+m₀σ_z)² ∝ I; NOT from U†U,")
    print("       NOT from an exponent.")
    print("  ★THE EXANTE'S WAGER «(1,d)»: holds EXACTLY at d=2, fails at d≥3.")
    print("     The live block of the seam = (1,2) ∀d; (d−2) directions — EXACT zeros along the node.")
    print("=" * 80)
    sys.stdout = tee.real
    logf.close()
    return code


if __name__ == "__main__":
    sys.exit(main())
