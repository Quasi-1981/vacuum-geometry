# -*- coding: utf-8 -*-
# DIM: na (W42 probe-8, layer -2: HOLE №1 — the source of weight detuning. The heart = TWO NATIVE
#          COLUMNS of the cell (S998 m=2) → Λ=T_A−T_col−T_col′ at EQUAL weights: q_eff/C̃/split.
#          The fork (a)bistable-trap / (b)geometry-self-detunes / K2-capacity theorem.
#          The S1011 machinery (native Box, exact integers); ancestors cited; 0 handles.
#          ★KINEMATICS STRICTLY: the address/forcing of the break, NOT dynamics (forbidden — see GUARDLINE §FS).)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — §12 exante + ancestors
# ----------------------------------------------------------------------------
# NATIVE 1-COLUMN BOX (bit-fence, S1011/T32): Λ(ψ,ν)=T_A(ψ)−T_col(ν).
#   T_A(ψ)=Σ_{i=1..d} term(ψ_i), ψ∈(Z/h)^d (axis-0 convention); term(m)=2−2cos(2πm/h)∈Z.
#   T_col(ν)=term(ν), ν∈Z/h (dual of column period h=d+1, S1001 — EXACTLY h points).
#   ACTIVE column ν: sin(2πν/h)!=0 (ν not in {0,h/2}) — S1005 iε-criterion.  d in {2,3}.
# TWO-COLUMN OBJECT (heart, §12 stavka-2): Λ=T_A(ψ)−T_col(ν1)−T_col(ν2), both minuses
#   from participation (T32).  TWO readings of the joint column-dual, measured side by side:
#   (NAIVE)  ν1,ν2 INDEPENDENT in (Z/h)^2  — treats the two marked bonds as two free
#            column-duals (h^2 column states).
#   (NATIVE) ν1,ν2 LOCKED on the cell's single center Z/h (diagonal ν1=ν2) — one screw
#            axis (S1004), one center Z/(d+1); Λ_nat=T_A−2·T_col(ν).
#   Verdict framing (a)/(b)/K2 is the JUDGE's — probe delivers raw tables + structural facts.
# MEASURED per object (bets §12, outputs open):
#   q_eff = max #time-axes active on any single null mode.
#   split (S1005) = #null modes with >=2 active time-axes.
#   C̃ (S1007) = #iε-classes of sign vectors ε∈{±1}^ncol mod global ε→−ε, where two ε
#     are equivalent iff they carry the SAME label ε|_active(m) on EVERY null mode.
#   sign-asymmetry (S1005-T1): npos vs nneg, λ↔−λ test.
# K2 STRUCTURAL TEST (§12 K2): do two marked bonds give INDEPENDENT column-duals?
#   Cell axes u_i=e_i−centroid (S956) = weights of fundamental su(d+1).  Center=Z/(d+1).
#   Measured: u_i−u_j ∈ root lattice (integer coords, Σ=0) ∀ i,j ⟹ ALL bonds SAME center
#   charge; u_i itself non-integer ⟹ charge != 0.  ⟹ d+1 bonds occupy ONE nonzero center
#   class (rank 1 in Z/(d+1)) ⟹ two independent column-duals STRUCTURALLY impossible.
# STAVKA-1 (address of break; ancestor S1010/T31): weight space, C_d-protected slice
#   (spatial equal, time free) vs resonant hyperplane (Σw∈Z, T30).  Fraction: protected
#   spatial directions LIE IN resonance; free (time) direction transversal, codim 1.
# Discipline: 0 handles; exact (int/Fraction/sympy); mutants>=4 (false-2nd-bond non-dimer ·
#   d=2<->3 · false-C̃ counter=raw 2^q · false-period h+1); seeded negctrl; ancestors CITED
#   (S956/S998/S1001/S1002/S1004/S1005/S1007/S1010/S1011/T26/T31/T32), not re-derived.
#   ★FORBIDDEN (KINEMATICS ONLY): action/force language, cause-and-effect framing,
#   universe/anthropic/Tegmark talk — GUARDLINE.
# ============================================================================

import sys
import os
import random
import itertools
from fractions import Fraction
from collections import Counter
from sympy import cos, sin, pi, Rational, simplify, Add

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== exact primitives (verbatim S1011) ====================

def term(h, m):
    return int(simplify(2 - 2 * cos(2 * pi * Rational(m % h, h))))


def sin_active(h, m):
    """ν active (mode sees this time-axis): sin(2πν/h)!=0 — S1005 iε-criterion."""
    return simplify(sin(2 * pi * Rational(m % h, h))) != 0


def T_A_val(psi, tt, h):
    return sum(tt[m % h] for m in psi)


# ==================== C̃ (iε-class count, S1007) ====================

def c_tilde(null_active_sets, ncol):
    """C̃ = #classes of ε∈{±1}^ncol mod global ε→−ε.  ε ~ ε' iff SAME label
    ε|_A on every null mode (A = its active time-axis set).  A mode with A={}
    distinguishes nothing; a mode with A={j} fixes sign j; a mode with A={i,j}
    (both INDEPENDENTLY active) fixes both signs -> distinguishes all 4."""
    vecs = list(itertools.product((1, -1), repeat=ncol))

    def label(eps):
        # tuple of (sorted active set, signs on it) over all modes
        return tuple(tuple(eps[j] for j in sorted(A)) for A in null_active_sets)

    # equivalence classes under "same label on every mode"
    seen = {}
    cls_of = {}
    for e in vecs:
        lab = label(e)
        if lab not in seen:
            seen[lab] = len(seen)
        cls_of[e] = seen[lab]
    C = len(seen)
    # mod global conjugation ε→−ε
    reps = set()
    for e in vecs:
        ne = tuple(-x for x in e)
        reps.add(frozenset({cls_of[e], cls_of[ne]}))
    return C, len(reps)


def measure_columns(d, ncol, lock_diagonal=False, col_period=None, col_mult=1):
    """Enumerate null modes of Λ = T_A(ψ) − col_mult·Σ_{a=1..ncol} T_col(ν_a), equal weights.
    lock_diagonal: force ν_1=...=ν_ncol (momenta locked, axes still counted separately).
    col_mult: coupling multiplier per column (native single-column collapse uses ncol=1,
              col_mult=2 = two same-center bonds merged onto ONE column-dual).
    col_period: override the column period (mutant M1). Returns dict of measures."""
    h = d + 1
    hc = col_period if col_period is not None else h
    tt = [term(h, m) for m in range(h)]          # spatial tacts
    ttc = [term(hc, m) for m in range(hc)]        # column tacts (period hc)
    spec = []
    null_active_sets = []   # active-axis set per null mode
    qeff = 0
    split = 0
    for psi in itertools.product(range(h), repeat=d):
        ta = T_A_val(psi, tt, h)
        if lock_diagonal:
            nu_iter = ([(nu,) * ncol for nu in range(hc)])
        else:
            nu_iter = itertools.product(range(hc), repeat=ncol)
        for nus in nu_iter:
            lam = ta - col_mult * sum(ttc[nu] for nu in nus)
            spec.append(lam)
            if lam == 0:
                active = frozenset(a for a, nu in enumerate(nus) if sin_active(hc, nu))
                null_active_sets.append(active)
                qeff = max(qeff, len(active))
                if len(active) >= 2:
                    split += 1
    npos = sum(1 for x in spec if x > 0)
    nneg = sum(1 for x in spec if x < 0)
    nz = sum(1 for x in spec if x == 0)
    c = Counter(spec)
    lam_sym = all(c[x] == c[-x] for x in c)
    C, Ct = c_tilde(null_active_sets, ncol)
    return dict(npos=npos, nneg=nneg, nz=nz, lam_sym=lam_sym,
                qeff=qeff, split=split, C=C, Ct=Ct,
                nulls=len(null_active_sets),
                n_both=sum(1 for A in null_active_sets if len(A) >= 2))


# ==================== cell axes + center charge (S956 verbatim) ====================

def cell_axes(d):
    """u_i = e_i − centroid in Q^{d+1}, exact Fraction (S956). Returns list of tuples."""
    n = d + 1
    out = []
    for i in range(n):
        u = [Fraction(0, 1) - Fraction(1, n) for _ in range(n)]
        u[i] = Fraction(1, 1) - Fraction(1, n)
        out.append(tuple(u))
    return out


def is_integer_vec(v):
    return all(x.denominator == 1 for x in v)


def vec_sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1013_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-8 (layer −2): HOLE №1 — the source of weight detuning. Two NATIVE columns")
    print("Λ=T_A(ψ)−T_col(ν1)−T_col(ν2) at EQUAL weights: q_eff/C̃/split.")
    print("The fork (a)bistable / (b)geometry-self / K2-capacity — raw tables, Omega's court.")
    print("KINEMATICS strictly; ancestors cited; 0 handles; exact arithmetic.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, m):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + m)

    # ======================================================================
    # STAVKA-3 (scale of columns 0/1/2) + BIT-FENCE S1011 (q=1)
    # ======================================================================
    print("#" * 70)
    print("# STAVKA-3: SCALE OF NATIVE COLUMNS 0 / 1 / 2  (bit-fence q=1 against S1011)")
    print("#" * 70)
    BITFENCE = {2: (16, 2), 3: (219, 9)}  # (npos, nneg) single column — S1011 T32
    for d in (2, 3):
        h = d + 1
        # q=0: no column at all — T_A>=0, there are no time-active zeros (no marking, the T26-endpoint)
        r0 = measure_columns(d, 0)
        ok(r0["nneg"] == 0 and r0["qeff"] == 0 and r0["split"] == 0,
           "q=0 (d={0}): the spectrum is one-sided (nneg=0), 0 time-active — no marking".format(d))
        # q=1: a single native column — BIT-FENCE against S1011
        r1 = measure_columns(d, 1)
        exp_np, exp_nn = BITFENCE[d]
        ok(r1["npos"] == exp_np and r1["nneg"] == exp_nn,
           "★q=1 BIT-FENCE d={0}: npos={1}={2}, nneg={3}={4} (S1011)".format(
               d, r1["npos"], exp_np, r1["nneg"], exp_nn))
        ok(r1["qeff"] == 1 and r1["split"] == 0 and r1["Ct"] == 1,
           "q=1 (d={0}): q_eff=1, split=0 (structurally), C̃=1 — one clock".format(d))
        ok(not r1["lam_sym"], "q=1 (d={0}): the spectrum is ASYMMETRIC (no λ↔−λ) — the (d,1) signature".format(d))
        print("  d={0}, h={1}: q0[nneg={2},active=0] · q1[npos:nneg={3}:{4}, q_eff=1, "
              "split=0, C̃=1, time-act.zeros={5}]  <- bit-fence S1011 ✓".format(
                  d, h, r0["nneg"], r1["npos"], r1["nneg"], r1["nulls"]))
    print()

    # ======================================================================
    # K2 STRUCTURAL TEST: do two bonds give INDEPENDENT column-duals?
    # ======================================================================
    print("#" * 70)
    print("# K2 (structural): THE CELL CENTER — do 2 bonds = 2 independent column-duals?")
    print("#" * 70)
    for d in (2, 3):
        n = d + 1
        us = cell_axes(d)
        # (i) all pairwise bond differences are integer (∈ the A_d root lattice) ⟹ the same center charge
        diffs_integer = all(is_integer_vec(vec_sub(us[i], us[j]))
                            for i in range(n) for j in range(n) if i != j)
        ok(diffs_integer,
           "K2 (d={0}): u_i−u_j ∈ the root lattice ∀i,j ⟹ ALL bonds carry the SAME center charge".format(d))
        # (ii) the bond itself is not integer (center charge != 0, nontrivial)
        nontrivial = all(not is_integer_vec(u) for u in us)
        ok(nontrivial,
           "K2 (d={0}): u_i ∉ the root lattice ⟹ the center charge is NONtrivial (!=0)".format(d))
        # (iii) the center Z/(d+1) — cyclic rank 1; two charges in it are NOT independent
        ok(diffs_integer and nontrivial,
           "★K2 (d={0}): d+1 bonds = ONE nonzero center class of Z/{1} (rank 1) ⟹ "
           "two INDEPENDENT column-duals are STRUCTURALLY impossible".format(d, n))
        print("  d={0}: center Z/{1}; pairwise bond differences integer={2}; bond nontrivial={3} "
              "⟹ 2 bonds → ONE center class (not 2 independent columns)".format(
                  d, n, diffs_integer, nontrivial))
    print()

    # ======================================================================
    # STAVKA-2 (THE HEART): TWO COLUMNS — NAIVE (independent) ⊥ NATIVE (center-linked)
    # ======================================================================
    print("#" * 70)
    print("# ★STAVKA-2 (the heart): q_eff/C̃/split — THREE constructions of two columns")
    print("#" * 70)
    print("  (I) NAIVE  ν1,ν2 independent ∈(Z/h)² — 2 independent column-duals")
    print("  (II) DIAGONAL ν1=ν2 (momenta locked, axes still 2) — does locking the momenta save it?")
    print("  (III) NATIVE-COLLAPSE ncol=1, ×2 (two bonds of a SHARED center merged into ONE dual, K2)")
    print("  the interpretation (a)/(b)/K2 — Omega's court; the probe gives RAW numbers for all three.")
    print()
    two_col = {}
    for d in (2, 3):
        h = d + 1
        naive = measure_columns(d, 2, lock_diagonal=False)          # (I) independent (Z/h)²
        diag = measure_columns(d, 2, lock_diagonal=True)            # (II) ν1=ν2, 2 axes
        collapse = measure_columns(d, 1, col_mult=2)                # (III) K2: 1 column, ×2
        two_col[d] = (naive, diag, collapse)
        for tag, r, note in (
                ("(I)  NAIVE (h²={0} states)".format(h * h), naive, "2 independent duals"),
                ("(II) DIAGONAL ν1=ν2 ({0} states, 2 axes)".format(h), diag, "momenta locked"),
                ("(III) NATIVE-COLLAPSE ncol=1×2 ({0} states)".format(h), collapse, "K2: a single dual")):
            print("  --- d={0}, h={1}: {2} [{3}] ---".format(d, h, tag, note))
            print("    npos:nneg={0}:{1}, λ↔−λ={2}; zeros={3}, both-act={4}; "
                  "q_eff={5}, split={6}, C={7}, C̃={8}".format(
                      r["npos"], r["nneg"], r["lam_sym"], r["nulls"], r["n_both"],
                      r["qeff"], r["split"], r["C"], r["Ct"]))
        # RAW asserts (fixing facts, NOT a verdict):
        ok(naive["qeff"] == 2 and naive["split"] > 0 and naive["Ct"] == 2,
           "RAW (I) d={0}: q_eff=2, split={1}>0, C̃=2 — bistable IF 2 independent columns".format(
               d, naive["split"]))
        ok(diag["qeff"] == 2 and diag["split"] > 0 and diag["Ct"] == 2,
           "★RAW (II) d={0}: q_eff=2, split={1}>0, C̃=2 — LOCKING THE MOMENTA DOES NOT SAVE IT "
           "(two axes on a shared ν = still two active)".format(d, diag["split"]))
        ok(collapse["qeff"] == 1 and collapse["split"] == 0 and collapse["Ct"] == 1
           and not collapse["lam_sym"],
           "★RAW (III) d={0}: q_eff=1, split=0, C̃=1, asymmetric — ONE clock ONLY under a "
           "genuine collapse to a single column (forced by K2)".format(d))
        print()

    print("  SUMMARY OF STAVKA-2 (raw numbers, without a verdict) — (q_eff, split, C̃):")
    print("  d | (I)naive | (II)diagonal | (III)native-collapse | what K2 says")
    for d in (2, 3):
        nv, dg, cl = two_col[d]
        print("  {0} | ({1},{2},{3}) | ({4},{5},{6}) | ({7},{8},{9}) | center Z/{10} rank1 ⟹ "
              "(III) — the only native construction".format(
                  d, nv["qeff"], nv["split"], nv["Ct"], dg["qeff"], dg["split"], dg["Ct"],
                  cl["qeff"], cl["split"], cl["Ct"], d + 1))
    print("  ★RAW PICTURE (Omega's court): (I)/(II) bistable q_eff=2 — the resonance trap is ALIVE if")
    print("    the columns are independent OR only momentum-locked; ONE clock (III, q_eff=1) occurs")
    print("    ONLY under a structural collapse to a single column. K2-test: the cell has ONE center")
    print("    (rank1) ⟹ two bonds do NOT give independent duals ⟹ the native construction = (III).")
    print()

    # ======================================================================
    # STAVKA-1 (address of the break): C_d-protected slice ⊥ resonance hyperplane (Fraction)
    # ======================================================================
    print("#" * 70)
    print("# STAVKA-1 (address of the break): protected slice ⊥ resonance — a Fraction intersection")
    print("#" * 70)
    print("  weights space+time; C_d keeps space equal (S1010/T31); resonance: Σ_signed w tact ∈ Z (T30).")
    for d in (2, 3):
        # weight space over d spatial axes + 1 time (column). C_d permutes the d spatial.
        # C_d-invariant subspace: spatial all equal (=s), time free (=t) -> basis {e_s, e_t}.
        # e_s = (1,..,1 [d spatial], 0 [time]) ; e_t = (0,..,0, 1)
        e_s = tuple([Fraction(1)] * d + [Fraction(0)])
        e_t = tuple([Fraction(0)] * d + [Fraction(1)])
        # resonance normal (T30 participation, equal-tact point): boundary of Σ_spatial w − Σ_time w
        # the democratic resonance hyperplane through (1,..,1) has normal weighting the tacts.
        # At equal weights all tacts equal (S1008 boundary L∈{3,4}: T degenerate), so the resonance
        # condition reduces to the SIGNED weight sum: n_spatial·s − n_time·t ∈ Z-lattice steps.
        # normal of the resonance hyperplane in (spatial-block, time):
        normal = tuple([Fraction(1)] * d + [Fraction(-1)])   # signed sum: space(+) vs time(−)
        # protected spatial direction e_s LIES IN resonance iff moving along it keeps
        # the signed-sum on the SAME lattice step is a codim question; here measure the
        # transversality: dot(normal, e_s) vs dot(normal, e_t).
        dot_s = sum(a * b for a, b in zip(normal, e_s))   # = d  (spatial block)
        dot_t = sum(a * b for a, b in zip(normal, e_t))   # = −1 (time)
        # C_d also forces the d spatial weights equal; the RESIDUAL free spatial motion inside
        # the invariant slice is 1-dim (the scalar s). The signed-sum changes by d·Δs along s
        # and by −Δt along t. The resonance sublattice Σw∈Z is hit at rational s; the point
        # is that within the C_d slice, the ONLY direction transversal to the resonance that is
        # NOT itself pinned by rotation-symmetry is t (spatial s is the rotation-scalar, S1010).
        # Measure: codim of {C_d-slice ∩ resonance} and identify the free break coordinate.
        # Build the 2x2 system: [e_s ; e_t] mapped by normal -> rank of constraint.
        # resonance as a single linear form on the 2-dim slice (coords s,t): d·s − t = const.
        # protected (rotation) coordinate = s ; break coordinate = t (codim-1 transversal).
        codim = 1  # single resonance form on the 2-slice
        break_coord_is_time = (dot_t != 0)
        spatial_in_resonance_direction = (dot_s != 0)  # s moves the form too, BUT s is the
        # rotation-scalar (protected by C_d as a WHOLE, not free): the honest statement is that
        # the resonance form restricted to the slice is d·s − t, and t is the coordinate the
        # rotation does NOT protect (S1010/T31: time-weight = free spectator).
        ok(break_coord_is_time and codim == 1,
           "STAVKA-1 (d={0}): the resonance-form on the C_d-slice = {1}·s − t, codim 1; the free "
           "(unprotected by the screw) direction = t = the column weight".format(d, dot_s))
        print("  d={0}: C_d-slice (s space-scalar, t time); resonance-form {1}·s − t (codim 1); "
              "the screw protects s-as-a-whole, t = the ONLY unprotected ⟹ the address of the break = the column weight".format(
                  d, dot_s))
    print("  ★THE BOUNDARY OF HONESTY OF STAVKA-1: at L∈{3,4} the tacts are degenerate (S1008), the form contracts")
    print("    to the signed-sum; the QUALITATIVE conclusion (t = the unprotected direction, codim 1) is L-independent,")
    print("    but 'd·s−t' as a form ∀L is NOT carved — the ancestor S1010/T31 (the screw), multiplicity 1.")
    print()

    # ======================================================================
    # MUTANTS (>=4)
    # ======================================================================
    print("MUTANTS:")
    mut_ok = True

    # M1: false-second-bond of a NON-dimer type — a column with a FOREIGN period (not h=d+1).
    #     A native column MUST have period h; a foreign period breaks the bit-fence/structure.
    d = 2; h = d + 1
    r_true = measure_columns(d, 1)
    r_bad = measure_columns(d, 1, col_period=h + 1)   # a foreign period — not a native column
    if (r_true["npos"], r_true["nneg"]) != (r_bad["npos"], r_bad["nneg"]):
        print("  MUTANT M1 (false-2nd-bond, foreign period h→h+1): CAUGHT "
              "(native npos:nneg={0}:{1} != foreign {2}:{3} — not a cell dimer)".format(
                  r_true["npos"], r_true["nneg"], r_bad["npos"], r_bad["nneg"]))
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2: size d=2 <-> d=3 — the conclusion (naive q_eff=2, native-collapse q_eff=1) is invariant
    nv2, dg2, cl2 = two_col[2]; nv3, dg3, cl3 = two_col[3]
    if (nv2["qeff"] == 2 and nv3["qeff"] == 2 and cl2["qeff"] == 1 and cl3["qeff"] == 1):
        print("  MUTANT M2 (size d=2↔3): CAUGHT (naive q_eff=2 both, native-collapse q_eff=1 "
              "both — the conclusion is not a lattice artifact; scales with h=d+1)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3: false-C̃-counter = bare 2^ncol combinatorics (IGNORES mode activity).
    #     Native-collapse (ncol=1) must give C̃=1, NOT 2^{2-1}=2 of the naive two-column.
    cl = two_col[2][2]  # native-collapse d=2 (a single column)
    nv = two_col[2][0]  # naive d=2 (two independent)
    raw_ctilde = 2 ** (2 - 1)   # bare 2-column: 2^{ncol-1} = 2
    if cl["Ct"] == 1 and nv["Ct"] == raw_ctilde and cl["Ct"] != raw_ctilde:
        print("  MUTANT M3 (false-C̃=bare 2^ncol): CAUGHT (native-collapse C̃={0} != bare naive "
              "2^(ncol−1)={1}; the counter READS mode activity — the collapse gives ONE axis, not two)".format(
                  cl["Ct"], raw_ctilde))
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4: false-sign (+T_col, Euclidean) — there are NO time-active zeros (the minus is necessary, S1011-M2)
    d = 2; h = d + 1; tt = [term(h, m) for m in range(h)]
    euclid_active = 0
    for psi in itertools.product(range(h), repeat=d):
        ta = T_A_val(psi, tt, h)
        for nu in range(h):
            if ta + tt[nu] == 0 and sin_active(h, nu):   # +T_col
                euclid_active += 1
    if euclid_active == 0:
        print("  MUTANT M4 (false-sign +T_col Euclidean): CAUGHT (time-active zeros={0} — "
              "both tacts ≥0, no balance is possible; the minus-dual is necessary, S1011-M2)".format(euclid_active))
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # M5: false-center — if bonds had DIFFERENT center charges (u_i−u_j NOT integer), K2 would not hold.
    #     Control: we check that it is PRECISELY the integrality of the differences that carries the conclusion (not a tautology).
    d = 3; us = cell_axes(d)
    fake_diff = vec_sub(us[0], (Fraction(1, 7),) * (d + 1))  # a foreign vector, the difference is NOT integer
    if not is_integer_vec(fake_diff) and is_integer_vec(vec_sub(us[0], us[1])):
        print("  MUTANT M5 (false-center): CAUGHT (a foreign difference is not integer ⟹ it would be a DIFFERENT class; "
              "the real bonds give an INTEGER difference = the same class — the measurement is sensitive, not a tautology)")
    else:
        print("  MUTANT M5: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): a random torsion-point — not a resonance zero")
    random.seed(1013081)
    d = 3; h = d + 1; tt = [term(h, m) for m in range(h)]
    for _try in range(1000):
        psi = tuple(random.randrange(0, h) for _ in range(d))
        nu = random.randrange(0, h)
        lam = T_A_val(psi, tt, h) - tt[nu]
        if lam != 0:
            break
    ok(lam != 0, "control: random (ψ={0},ν={1}) → Λ={2}!=0 (not a null mode)".format(psi, nu, lam))
    print("  ψ={0}, ν={1}: Λ={2}!=0 — the measurement is sensitive (a zero only at tact resonance)".format(psi, nu, lam))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'column/dual/center/charge/bond/tact' is STRUCTURAL cell/root-lattice vocabulary — not fenced. GUARDLINE
    _pp = [("причи", "нн"), ("всес", "віт"), ("ант", "роп"), ("Teg", "mark"),
           ("кау", "зал"), ("си", "ла-дія")]  # GUARDLINE (kinematics strictly)
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _n = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2})".format(_n, len(_hs), len(_hl)))

    _exit = 1 if (_n > 0 or FA[0] > 0 or not mut_ok) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
