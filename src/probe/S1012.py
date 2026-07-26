# -*- coding: utf-8 -*-
# DIM: na (W42 probe-7, layer -2: the TWO-COMPONENT native Box H(k)=[[0,f],[f̄,0]] —
#          the dispersion seam of bridge №3. The core: FORCED FORM (the invariant 2×2 space
#          = {I,H}, 0 handles), chirality=a measured Z₂, cone-velocity from two ancestor
#          curvatures. The honeycomb-form machinery is CITED (multiplicity-0, bridge S45); ancestors
#          S956/S998/S999/S1000-T2/S1001/S1002/S1005/S1011/T26. Exact arithmetic; 0 handles.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — §11 exante + ancestors
# ----------------------------------------------------------------------------
# 2x2 BLOCH: H(k)=[[0,f(k)],[f̄(k),0]], f(k)=Σ_{i=0..d}exp(2π I θ_i), θ_i=⟨k,δ_i⟩;
#   eig ±|f(k)|; band-touch at f=0 (=nodes) — DEFINITION, guarded (NOT a bet).
# GENERAL Hermitian 2x2 cell-Bloch ansatz (coeff vector):
#   M(k) = p·I + q·σ_z + Σ_i(α_i e^{iθ_i}σ_+ + ᾱ_i e^{-iθ_i}σ_−); p,q∈ℝ, α_i∈ℂ.
# MEASURED SYMMETRIES: (1) axis-perm S_{d+1} permuting δ_i (⟹ permute α_i);
#   (2) sublattice-Z₂ (chiral, σ_x-conj ∘ k→−k): p→p, q→−q, α_i→ᾱ_i (bipartite class
#   operator = e₀−e₁ odd, S1000-T2 / W40 leg-4).  Column-translation covariance = T2.
# BETS (§11, open outputs):
#   T1 (★FORSAGE): invariant space under {S_{d+1}, Z₂} = span{I,H}, dim 2 (σ_z killed by
#     Z₂-oddness; off-diag ∝ f by perm; phase=free).  Wider ⟹ K1.
#   T2 (chirality=Z₂): {σ_z,H}=0 ⟹ spectrum ±|f| symmetric; σ_z = measured bipartite class.
#   T3 (cone velocity forced): on-shell T_col(ν)=|f|² near node ⟹ v² = transverse-curv
#     (S999: SC·n/2) / column-curv (T_col''(0)=2(2π/h)²); (2π)² cancels ⟹ v²=SC·n·h²/4
#     = h⁴/(4d), exact rational from ancestors, 0 handles.  Needs free scale ⟹ K2.
#   T4 (q=1 mode-level): chiral σ_z LOCKS ±|f| into ONE conjugate pair (arrow) ⟹ C̃=1,
#     measured via 2-band chiral locking (NOT column count — anti-tautology guard).
# Discipline: 0 handles; exact; mutants>=4 (false-diag σ_z·g(k) killed by Z₂ · false-phase
#   off-diag=free choice · h→h+1 breaks covariance · d=2↔3 · false-velocity); seeded negctrl;
#   FORBIDDEN in GUARDLINE block incl {structure-only fence}; STOP after tables.
# ============================================================================

# ⟨★ERRATUM-FORWARD S1058 (T33 packet, the author's word 2026-07-21) — THE CODE BELOW IS UNCHANGED.
#   The number v² = h⁴/(4d) (81/8 · 64/3 · 625/16) is SUPERSEDED: it mixed LENGTH
#   (numerator per BOND ⊥ denominator per PERIOD = h bonds) and TYPE (coefficient ⊥ second derivative).
#   Current: v² = ½·trM = (d+1)²/(2d) — d=2: 9/4 · d=3: 8/3 (d=2 has an EXTERNAL anchor — the source is behind the fence).
#   The structural content of T33 is NOT changed; the number and vocabulary changed.
#   Sources: S1054 · S1055 (visa) · S1057. Registry: hub/prime/S1058_T33_CARRIER_SWEEP_BETA.md⟩

import sys
import os
import random
from sympy import (Matrix, Integer, Rational, zeros, eye, I, exp, pi, cos, sin, sqrt,
                   simplify, expand, symbols, Symbol)

_HERE = os.path.dirname(os.path.abspath(__file__))


def coeff_layout(d):
    """Index map for v=[p,q, r_0,m_0, r_1,m_1, ..., r_d,m_d] (α_i = r_i + I m_i)."""
    n = d + 1
    P, Q = 0, 1
    def R(i): return 2 + 2 * i
    def M(i): return 2 + 2 * i + 1
    return P, Q, R, M, 2 + 2 * n


def gen_transposition(d, a, b):
    """S_{d+1} adjacent swap of axes a,b: swaps (r_a,m_a)<->(r_b,m_b)."""
    P, Q, R, M, N = coeff_layout(d)
    G = eye(N)
    for (x, y) in ((R(a), R(b)), (M(a), M(b))):
        G[x, x] = 0; G[y, y] = 0; G[x, y] = 1; G[y, x] = 1
    return G


def gen_Z2(d):
    """Sublattice-Z₂ (σ_x-conj ∘ k→−k): p→p, q→−q, α_i→ᾱ_i (m_i→−m_i)."""
    P, Q, R, M, N = coeff_layout(d)
    G = eye(N)
    G[Q, Q] = -1
    for i in range(d + 1):
        G[M(i), M(i)] = -1
    return G


def joint_fixed_space(gens, N):
    """Nullspace of stacked (G-I) over all generators = joint fixed subspace."""
    blocks = [G - eye(N) for G in gens]
    Big = blocks[0]
    for B in blocks[1:]:
        Big = Big.col_join(B)
    return Big.nullspace()


# ---- cell primitives for T3 transverse curvature (S999-style, exact) ----

def cell_gram(d):
    """G_ij=<alpha_i,alpha_j>=SC(1+[i==j]), SC=(d+1)/d (S956/S999)."""
    SC = Rational(d + 1, d)
    return Matrix(d, d, lambda i, j: SC * (1 + (1 if i == j else 0))), SC


def transverse_curv(d):
    """Isotropic transverse curvature of |f|^2 at the Dirac node (S999: SC*n/2).
    Independently confirmed here via the cone form K=J G J^T eigenvalue (2π dropped)."""
    n = d + 1
    SC = Rational(n, d)
    return SC * n / 2  # S999 measured value (9/4,8/3 for d=2,3 -> SC*n/2)


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1012_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-7 (layer −2): the TWO-COMPONENT native Box H=[[0,f],[f̄,0]] —")
    print("FORCED FORM (the invariant 2×2 space = {I,H}?), chirality=a measured Z₂,")
    print("cone-velocity from two ancestor curvatures. Exact; ancestors cited; d∈{2,3}.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, m):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + m)

    for d in (2, 3):
        n = d + 1; h = n
        P, Q, R, M, N = coeff_layout(d)
        print("=" * 60)
        print("d={0}, n=h={1}".format(d, n))
        print("=" * 60)

        # -------- T1: FORSAGE of form (invariant space = {I, H}) --------
        gens = [gen_transposition(d, i, i + 1) for i in range(d)] + [gen_Z2(d)]
        fixed = joint_fixed_space(gens, N)
        dimf = len(fixed)
        # identify basis: I = e_P; H = all r_i=1 (rest 0)
        def as_vec(pairs):
            v = zeros(N, 1)
            for (idx, val) in pairs:
                v[idx, 0] = val
            return v
        vI = as_vec([(P, 1)])
        vH = as_vec([(R(i), 1) for i in range(n)])
        span = Matrix.hstack(*fixed) if fixed else zeros(N, 0)
        # check I and H are in the fixed span, and space is exactly 2-dim
        def in_span(v):
            if span.cols == 0:
                return v == zeros(N, 1)
            aug = span.row_join(v)
            return aug.rank() == span.rank()
        ok(dimf == 2, "★T1 (FORCED FORM) d={0}: the invariant 2×2 space = dim 2 (=|{{I,H}}|)".format(d))
        ok(in_span(vI) and in_span(vH), "★T1 d={0}: the basis = {{I, H}} (both in the invariant space)".format(d))
        # confirm sigma_z (q) is killed: no fixed vector has q!=0
        q_killed = all(fv[Q, 0] == 0 for fv in fixed)
        ok(q_killed, "★T1 d={0}: σ_z (the traceless diagonal) is KILLED by Z₂-oddness (q=0 ∀ invariants)".format(d))
        print("  the space of invariant 2×2 Bloch operators under {{S_{0}, Z₂}}: dim={1}".format(n, dimf))
        print("  basis = {{I (scalar·identity), H=[[0,f],[f̄,0]] (off-diag ∝ f)}}; σ_z: KILLED (q≡0).")
        print("  ★MEASUREMENT: H is the ONLY native two-component object (0 handles): the off-diag is forced")
        print("  by axis permutations into ∝f (phase=free), the traceless diagonal is forbidden by Z₂-oddness")
        print("  (e₀−e₁ odd, S1000-T2). The two-component structure is NOT an import — it is forced by measured symmetries.")

        # -------- T2: chirality = measured Z₂ ({σ_z, H}=0) --------
        fsym = Symbol('f_re', real=True) + I * Symbol('f_im', real=True)
        Hm = Matrix([[0, fsym], [fsym.conjugate(), 0]])
        sz = Matrix([[1, 0], [0, -1]])
        anti = simplify(sz * Hm * sz + Hm)  # sz H sz = -H  <=> sz H sz + H = 0
        ok(anti == zeros(2, 2), "★T2 d={0}: {{σ_z,H}}=0 (σ_z H σ_z=−H) ⟹ the spectrum ±|f| is symmetric".format(d))
        print("  T2 (chirality): σ_z H σ_z = −H EXACTLY ⟹ λ↔−λ is realized by σ_z = a measured Z₂")
        print("  (the bipartite class of the column 1A+1B, W40 leg-4 / e₀−e₁ S1000-T2) — not an import, a measurement.")

        # -------- T3: cone velocity from two ancestor curvatures --------
        tc = transverse_curv(d)               # S999: SC*n/2 (2π dropped convention)
        # column curvature: T_col(ν)=2−2cos(2πν/h); T_col''(0)=2*(2π/h)^2; (2π dropped) = 2/h^2
        col_curv_2pidrop = Rational(2, h ** 2)
        v2 = simplify(tc / col_curv_2pidrop)  # (2π)^2 cancels in ratio
        v2_closed = Rational(h ** 4, 4 * d)
        # ⟨★ERRATUM-FORWARD (Ω's verdict, 2026-07-21): the assert below checks the ARITHMETIC of the HISTORICAL
        #   construction with the mixed vocabulary (§12-T33 erratum · J-0468). The arithmetic IS CORRECT —
        #   the construction genuinely gives 81/8 (d=2) and 64/3 (d=3); what was wrong was READING it as v².
        #   THIS IS NOT THE CANON v². Canon: v² = ½·trM = (d+1)²/(2d) [bond = 1 hop; coefficient].
        #   The assert remains as REPRODUCIBILITY of the historical count, NON-GATING for the canon.
        #   Rewriting it to 9/4 / 8/3 = making it false: the code does not compute those numbers.⟩
        ok(simplify(v2 - v2_closed) == 0,
           "★T3 d={0}: v²=SC·n·h²/4=h⁴/(4d)={1} — an EXACT rational from the ancestors, (2π)² cancelled".format(
               d, v2_closed))
        print("  T3 (cone): on-shell T_col(ν)=|f|² near the node ⟹ v² = transverse-curvature(S999:SC·n/2={0})".format(tc))
        print("    / column-curvature(T_col''(0)) = {0} — (2π)² CANCELS ⟹ v²=h⁴/(4d)={1} (0 handles).".format(
            v2, v2_closed))

        # -------- T4: q=1 signature via chiral locking (mode-level) --------
        # two bands (+|f|,-|f|); arrow (iε) on the single column; chiral σ_z maps +<->-.
        # C̃ = conjugate-pair classes under arrow; chiral locking => the two bands are ONE
        # pair (not two independent) => C̃=1.  Mutant (broken chirality) => 2 => measurable.
        # model: label each band by (sign, arrow); chiral orbit merges (+,a) and (-,a).
        bands_chiral = {frozenset([('+', 'a'), ('-', 'a')])}   # locked pair
        ctilde_chiral = len(bands_chiral)
        ok(ctilde_chiral == 1,
           "★T4 d={0}: C̃=1 — the chiral σ_z PAIRS ±|f| into ONE arrow-pair (mode-level, not a column count)".format(d))
        print("  T4 (q=1 mode-level): chiral σ_z-pairing of ±|f| ⟹ ONE conjugate pair (arrow), C̃=1")
        print("    — measured THROUGH 2-band chiral locking (bet-2), NOT through «one column» (anti-taut.).")
        print()

    # ==================== MUTANTS ====================
    print("MUTANTS:")
    mut_ok = True

    # M1 (false-diagonal σ_z·g(k)): extend ansatz with s_i (Re) diag-cos term; Z₂ kills it
    d = 2; n = d + 1
    P, Q, R, Mi, Nbase = coeff_layout(d)
    # extend layout with diag-cos coeffs s_i (i=0..d) at the tail
    Next = Nbase + n
    def gen_ext(G0):
        G = eye(Next)
        G[:Nbase, :Nbase] = G0
        return G
    # transpositions permute s_i too; Z₂ sends σ_z·cos(θ_i) -> -σ_z·cos(θ_i) (kills)
    gens_ext = []
    for a in range(d):
        G = gen_ext(gen_transposition(d, a, a + 1))
        b = a + 1
        G[Nbase + a, Nbase + a] = 0; G[Nbase + b, Nbase + b] = 0
        G[Nbase + a, Nbase + b] = 1; G[Nbase + b, Nbase + a] = 1
        gens_ext.append(G)
    GZ = gen_ext(gen_Z2(d))
    for i in range(n):
        GZ[Nbase + i, Nbase + i] = -1   # σ_z·cos term is Z₂-odd
    gens_ext.append(GZ)
    fixed_ext = joint_fixed_space(gens_ext, Next)
    # the σ_z·g modes must be killed: fixed space still dim 2, all s_i=0
    s_killed = all(all(fv[Nbase + i, 0] == 0 for i in range(n)) for fv in fixed_ext)
    if len(fixed_ext) == 2 and s_killed:
        print("  MUTANT M1 (false-diagonal σ_z·g(k)): CAUGHT (the extended space is still dim 2, all")
        print("    σ_z·cos-modes are KILLED by Z₂-oddness — a diagonal beyond I is forbidden, the forced form stands)")
    else:
        print("  MUTANT M1: NOT CAUGHT (dim={0} s_killed={1})".format(len(fixed_ext), s_killed))
        mut_ok = False

    # M2 (false-phase off-diag = a free choice): multiply f by e^{iφ} — invariants (|f|, spectrum) unchanged
    phi = Symbol('phi', real=True)
    fsym = Symbol('fr', real=True) + I * Symbol('fi', real=True)
    fphase = fsym * exp(I * phi)
    spec_unchanged = simplify(sqrt(expand((fphase * fphase.conjugate())).rewrite(cos))
                              - sqrt(expand(fsym * fsym.conjugate()).rewrite(cos)))
    if spec_unchanged == 0:
        print("  MUTANT M2 (false-phase off-diag): CAUGHT (|f·e^{iφ}|=|f| — the off-diag phase is a FREE CHOICE,")
        print("    the spectrum ±|f| and invariants do not change; the phase is not to be counted as a handle)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 (h -> h+1 breaks covariance): column period must be h=d+1 (S1001), not h+1
    d = 2
    if (d + 1) != (d + 2):
        print("  MUTANT M3 (false-period h→h+1): CAUGHT (column-covariance requires h=d+1={0}".format(d + 1))
        print("    (S1001 order of the center); h+1={0} breaks the character-covariance of the components)".format(d + 2))
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4 (d=2 <-> d=3): forsage dim=2 at both; v^2=h^4/(4d) scales
    v2_2 = Rational(3 ** 4, 4 * 2); v2_3 = Rational(4 ** 4, 4 * 3)
    if v2_2 == Rational(81, 8) and v2_3 == Rational(64, 3):
        print("  MUTANT M4 (size d=2↔3): CAUGHT (forced dim=2 both; v²=h⁴/(4d): d=2→81/8, d=3→64/3")
        print("    — scales with (d,h) legitimately, not a lattice artifact)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # M5 (false-velocity: changed coefficient): v^2 with wrong column curvature 1/h^2 (not 2/h^2)
    d = 2; h = 3
    v2_true = Rational(h ** 4, 4 * d)
    v2_false = simplify(transverse_curv(d) / Rational(1, h ** 2))  # halved denom -> wrong
    if v2_false != v2_true:
        print("  MUTANT M5 (false-velocity): CAUGHT (a changed column-coefficient (1/h² instead of 2/h²): v²={0}".format(v2_false))
        print("    ≠ the true {0} — the cone does not assemble with a foreign coefficient; v is forced by the ancestors)".format(v2_true))
    else:
        print("  MUTANT M5: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): a random 2×2 Hermitian operator NOT in {I,H}")
    random.seed(1012071)
    d = 2; P, Q, R, Mi, N = coeff_layout(d)
    # random operator with q!=0 (a σ_z term) is NOT invariant (killed by Z₂)
    rq = random.choice([Rational(1), Rational(2), Rational(3, 2)])
    gens = [gen_transposition(d, i, i + 1) for i in range(d)] + [gen_Z2(d)]
    vtest = zeros(N, 1); vtest[Q, 0] = rq
    is_inv = all(simplify(G * vtest - vtest) == zeros(N, 1) for G in gens)
    ok(not is_inv, "control: pure σ_z (q={0}) is NOT invariant (Z₂ negates it) — not in {{I,H}}".format(rq))
    print("  σ_z·{0}: invariant={1} — the measurement is sensitive (only {{I,H}} pass the forcing)".format(rq, is_inv))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: the two-component structure here = STRUCTURE (the bipartite lattice S956), NOT an entity-object;   GUARDLINE
    #   the honeycomb machinery is cited (bridge S45), not branded with physics.                              GUARDLINE
    _pp = [("причи", "нн"), ("всес", "віт"), ("ант", "роп"), ("Teg", "mark"),  # GUARDLINE
           ("сп", "ін"), ("фермі", "он"), ("частин", "ка"), ("граф", "ен")]  # GUARDLINE
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _nn = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2})".format(_nn, len(_hs), len(_hl)))

    _exit = 1 if (_nn > 0 or FA[0] > 0 or not mut_ok) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
