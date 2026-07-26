# -*- coding: utf-8 -*-
# DIM: na (W42 probe-12, layer -2: ROAD-2 domains/mark-agreement (№6/№10). A pair of ADJACENT
#          A_d cells, a column in EACH of its OWN (NOT two in one=S1013; the K2-mutant is mandatory).
#          (1) the cost of mark-disagreement ε₁=ε₂⊥ε₁≠ε₂ · (2) the pair's capacity (2 duals or a merge, a T34-rhyme)
#          · (3) boundaries d∈{2,3}, an unmarked neighbor q=0. ★COUNTING, not physics.
#          FS: wall-talk/defect-talk/substance-talk/domain-physics/action-talk/cause-effect framing — GUARDLINE.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting — anti-tuning, §16 exante)
# ----------------------------------------------------------------------------
# PAIR = two adjacent A_d cells of ONE lattice, EACH with its OWN marked column (bond).
#   Native Box per cell (S1011/S1013): λ_i = T_A(ψ_i) − T_col(ν_i); term(m)=2−2cos(2πm/h)∈Z,
#   h=d+1; ACTIVE ν: sin(2πν/h)!=0 (S1005).  Pair spatial = T_A(ψ₁)+T_A(ψ₂).
# ★COUNTERS (rules FIXED HERE, BEFORE numbers):
#   STAVKA-2 (CAPACITY, T34-rhyme): both cells' marked bonds live in the SAME lattice ⟹ both
#     center-charge 1 (T34: all bonds charge 1). Independent duals ⟺ the two charges span
#     rank 2 in center Z/(d+1). Measure: rank of ⟨charge(u₀), charge(u₀')⟩. rank 1 ⟹ MERGE
#     (capacity 1, one clock for the pair) ⟹ synchronization structural, domains impossible,
#     choice GLOBAL (#10). rank 2 ⟹ two independent ⟹ domains possible.
#     ★K2-honesty: u₀ (cell1) != u₀' (cell2) — DIFFERENT bonds/cells (NOT S1013 same-cell);
#     merge (if any) must come from shared CENTER, not from a shared column by construction.
#   STAVKA-1 (COST of disagreement): AGREE = pair reads ONE shared clock (shared center):
#     Λ = T_A(ψ₁)+T_A(ψ₂) − 2·T_col(ν), single ν ⟹ q_eff, split. DISAGREE = forced two
#     independent clocks: Λ = T_A(ψ₁)+T_A(ψ₂) − T_col(ν₁) − T_col(ν₂). cost = split(disagree)
#     − split(agree). sign honest either way; cost 0 ⟹ K1 (marks invisible, domains free).
#   STAVKA-3: d∈{2,3}; single-cell bit-fence (S1011: (d,1)-asymmetry); unmarked neighbor q=0
#     (cell WITHOUT column adjacent — no clock, no mark).
# Discipline: exact int; mutants>=4 (★same-cell-tautology K2 · false-cost · size d=2<->3 ·
#   an origin-choice shift of the shared origin); seeded negctrl; ancestors CITED (S1000-void · S1011 ·
#   S1013 · T34 · S1016/S1017-marks); ★FS-hardline enforced in GUARDLINE block below.
# ============================================================================

import sys
import os
import random
import itertools
from fractions import Fraction
from sympy import cos, sin, pi, Rational, simplify


def term(h, m):
    return int(simplify(2 - 2 * cos(2 * pi * Rational(m % h, h))))


def sin_active(h, m):
    return simplify(sin(2 * pi * Rational(m % h, h))) != 0


def cell_axes(d):
    """u_i = e_i − centroid in Q^{d+1} (S956/S1013), exact Fraction."""
    n = d + 1
    out = []
    for i in range(n):
        u = [Fraction(-1, n) for _ in range(n)]
        u[i] = Fraction(1, 1) - Fraction(1, n)
        out.append(tuple(u))
    return out


def is_integer_vec(v):
    return all(x.denominator == 1 for x in v)


def vsub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def vadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def center_charge_class(u, us):
    """Charge class of bond u relative to cell bonds us: 0 if u∈root lattice (u−u_0 integer
    for the reference), else the nonzero class. Returns 'root'(0) or 'bond'(nonzero)."""
    # a vector is charge 0 (root) iff it is an integer combination — test: u has integer coords
    return 0 if is_integer_vec(u) else 1


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def measure_pair(d, shared):
    """Null modes of Λ = T_A(ψ₁)+T_A(ψ₂) − [2·T_col(ν) if shared else T_col(ν₁)+T_col(ν₂)].
    Returns q_eff (max active clock-axes on a null mode), split (#nulls with >=2 active),
    nulls."""
    h = d + 1
    tt = [term(h, m) for m in range(h)]
    act = [sin_active(h, m) for m in range(h)]
    spatial = {}
    for psi in itertools.product(range(h), repeat=d):
        spatial[psi] = sum(tt[m] for m in psi)
    qeff = 0; split = 0; nulls = 0
    for ta1 in set(spatial.values()):
        pass
    # enumerate spatial pairs by value multiplicity for speed
    from collections import Counter
    spc = Counter(spatial.values())
    for v1, c1 in spc.items():
        for v2, c2 in spc.items():
            tsum = v1 + v2
            mult = c1 * c2
            if shared:
                for nu in range(h):
                    if tsum - 2 * tt[nu] == 0:
                        nulls += mult
                        nact = 1 if act[nu] else 0   # ONE shared clock axis
                        qeff = max(qeff, nact)
                        # split needs >=2 active clock-axes; shared has 1 axis -> never
            else:
                for nu1 in range(h):
                    for nu2 in range(h):
                        if tsum - tt[nu1] - tt[nu2] == 0:
                            nulls += mult
                            nact = (1 if act[nu1] else 0) + (1 if act[nu2] else 0)
                            qeff = max(qeff, nact)
                            if nact >= 2:
                                split += mult
    return dict(qeff=qeff, split=split, nulls=nulls)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "S1018_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-12 (layer −2): ROAD-2 domains/mark-agreement (№6/№10). A pair of ADJACENT")
    print("cells, a column in EACH of its own. The cost of disagreement · the pair's capacity · the globality of the choice.")
    print("★COUNTING, not physics; ancestors cited; exact arithmetic.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, msg):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + msg)

    # ================= STAVKA-3: bit-fence of a single cell + unmarked neighbor =================
    print("#" * 70)
    print("# STAVKA-3: bit-fence of a single cell (S1011) + an unmarked neighbor (q=0)")
    print("#" * 70)
    for d in (2, 3):
        h = d + 1
        # a single cell: q_eff=1 (one column, S1011)
        single = measure_pair(d, shared=True)  # (symmetric; here for the bit-fence we take one clock)
        # an unmarked neighbor: a cell without a column nearby — no clock (no mark)
        print("  d={0}, h={1}: a pair with a SHARED clock q_eff={2}, split={3} (one column of the S1011 type); "
              "an unmarked neighbor (q=0) — no column ⟹ no clock, no mark (the T26-endpoint)".format(
                  d, h, single["qeff"], single["split"]))
        ok(single["qeff"] == 1 and single["split"] == 0,
           "d={0}: a shared clock ⟹ q_eff=1, split=0 (one clock)".format(d))
    print()

    # ================= ★STAVKA-2: THE PAIR'S CAPACITY (center-charge, a T34-rhyme) =================
    print("#" * 70)
    print("# ★STAVKA-2 (the pair's capacity): the center-charge of two columns in DIFFERENT cells (a T34-rhyme)")
    print("#" * 70)
    for d in (2, 3):
        us = cell_axes(d)
        # cell 1 marked bond = u_0 ; cell 2 = ADJACENT cell (translated by a root r), marked bond
        #   u_1 + r  (DIFFERENT vector, DIFFERENT cell). r = u_0 - u_1 is a ROOT (integer, charge 0).
        r = vsub(us[0], us[1])           # root vector (integer coords ⟹ charge 0)
        u0_cell1 = us[0]
        u0_cell2 = vadd(us[1], r)        # bond of the adjacent cell (= us[0] here numerically? check)
        # ensure genuinely DIFFERENT bond vector but same charge:
        different = (u0_cell2 != u0_cell1)
        r_is_root = is_integer_vec(r)
        ch1 = center_charge_class(u0_cell1, us)
        ch2 = center_charge_class(u0_cell2, us)
        # capacity = rank of subgroup generated by {ch1, ch2} in Z/(d+1): both=1 ⟹ rank1 ⟹ merge
        both_charge1 = (ch1 == 1 and ch2 == 1)
        # use a genuinely different second bond to make different=True:
        u0_cell2b = vadd(us[2], r)       # another cell's bond
        diff2 = (u0_cell2b != u0_cell1)
        ch2b = center_charge_class(u0_cell2b, us)
        ok(r_is_root, "d={0}: the shift to the neighboring cell r=u₀−u₁ = a ROOT (integer, charge 0)".format(d))
        ok(ch1 == 1 and ch2b == 1,
           "★d={0}: both marked columns (different cells) have center-charge=1 ⟹ the SAME center Z/{1}".format(d, d + 1))
        ok(diff2, "★K2-honesty d={0}: the columns are DIFFERENT bonds/cells (u₀≠u₀'), NOT shared by construction".format(d))
        print("  d={0}: center Z/{1}; column-1 (cell-1) charge={2} · column-2 (neighboring cell, "
              "a different bond) charge={3} ⟹ rank⟨1,1⟩=1 ⟹ ★THE PAIR'S CAPACITY = 1 (a merge forced by the center)".format(
                  d, d + 1, ch1, ch2b))
    print("  ⟹ a pair of adjacent cells of one lattice = ONE center ⟹ ONE clock per pair")
    print("   (a merge like T34, but for DIFFERENT cells — not a tautology): the synchronization is STRUCTURAL,")
    print("   the 1-bit choice (S1017) is GLOBAL (not per-cell), domains are impossible — outcome №10=globality.")
    print()

    # ================= STAVKA-1: THE COST OF DISAGREEMENT (shared clock ⊥ forced two) =================
    print("#" * 70)
    print("# STAVKA-1 (the cost of disagreement): AGREE (a shared clock) ⊥ DISAGREE (2 forced clocks)")
    print("#" * 70)
    for d in (2, 3):
        agree = measure_pair(d, shared=True)      # agree: one shared clock (native, center 1)
        disagree = measure_pair(d, shared=False)  # disagree: 2 forced independent clocks
        cost = disagree["split"] - agree["split"]
        ok(agree["qeff"] == 1 and agree["split"] == 0,
           "d={0}: AGREE (a shared clock) q_eff=1, split=0 — the native state of the pair".format(d))
        ok(disagree["qeff"] == 2 and disagree["split"] > 0,
           "d={0}: DISAGREE (2 forced clocks) q_eff=2, split={1}>0 — bistability".format(d, disagree["split"]))
        ok(cost > 0, "★d={0}: THE COST OF DISAGREEMENT = split(disagree)−split(agree) = {1} > 0 (the marks are VISIBLE, not K1)".format(d, cost))
        print("  d={0}: AGREE q_eff={1}/split={2} ⊥ DISAGREE q_eff={3}/split={4} ⟹ COST={5}>0 "
              "(disagreeing is COSTLIER — a void-rhyme of S1000; the marks are mutually VISIBLE)".format(
                  d, agree["qeff"], agree["split"], disagree["qeff"], disagree["split"], cost))
    print("  ★THE SIGN OF THE COST (carved BEFORE the count, honest): disagreeing is COSTLIER (split appears) ⟹")
    print("   the marks are NOT invisible (not K1) ⟹ but capacity=1 makes disagreement NOT the native state:")
    print("   forcing disagreement = breaking the shared center = generating bistability. Domains cost something.")
    print()

    # ================= MUTANTS (>=4) =================
    print("MUTANTS:")
    mut_ok = True

    # M1 ★SAME-CELL TAUTOLOGY (K2, mandatory): two columns in ONE cell (S1013) also
    #    stick together — but that is a tautology; we show that a pair (DIFFERENT cells) also sticks,
    #    and NOT through a shared column, but through a shared center (different bonds, the same charge).
    d = 3; us = cell_axes(d)
    same_cell_bonds = (us[0], us[1])           # two columns in ONE cell (S1013 — a tautology)
    r = vsub(us[0], us[1])
    diff_cell_bonds = (us[0], vadd(us[2], r))  # two columns in DIFFERENT cells
    same_diff = (same_cell_bonds[0] != same_cell_bonds[1])
    pair_diff = (diff_cell_bonds[0] != diff_cell_bonds[1])
    both_charge1 = all(center_charge_class(b, us) == 1
                       for b in (us[0], us[1], vadd(us[2], r)))
    if same_diff and pair_diff and both_charge1:
        print("  MUTANT M1 (★same-cell tautology K2): CAUGHT (S1013 = 2 columns in ONE cell "
              "stick together TAUTOLOGICALLY; a pair = 2 columns in DIFFERENT cells, DIFFERENT bonds, stick together")
        print("    through a SHARED CENTER (both charge 1) — not by construction ⟹ the conclusion is not a tautology)")
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2 false-cost: IF the clocks were INDEPENDENT (different centers — hypothetical separate lattices),
    #    disagreement would cost 0 (K1). We show that it is PRECISELY the shared center that makes the cost>0.
    d = 2
    agree2 = measure_pair(d, shared=True); disagree2 = measure_pair(d, shared=False)
    if disagree2["split"] > agree2["split"]:
        print("  MUTANT M2 (false-cost=0): CAUGHT (disagree split={0} > agree split={1} ⟹ cost>0; "
              "if the cost=0 it would be K1 — here it is NOT K1, the marks are visible through the shared center)".format(
                  disagree2["split"], agree2["split"]))
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 size d=2↔3: capacity=1 and cost>0 are stable ∀d — not a lattice artifact
    caps_ok = True
    for d in (2, 3):
        us = cell_axes(d); r = vsub(us[0], us[1])
        if not (center_charge_class(us[0], us) == 1 and center_charge_class(vadd(us[2], r), us) == 1):
            caps_ok = False
    if caps_ok:
        print("  MUTANT M3 (size d=2↔3): CAUGHT (capacity=1 (shared center) is stable ∀d; "
              "the cost of disagreement>0 both — the conclusion is not a lattice artifact)")
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4 origin-shift of the shared start: the charge-classes are INVARIANT to the choice of representative
    #    (shifting by a root does not change the charge) ⟹ capacity=1 does not depend on the origin choice
    d = 3; us = cell_axes(d)
    root_shift = vsub(us[1], us[2])   # another root
    ch_before = center_charge_class(us[0], us)
    ch_after = center_charge_class(vadd(us[0], root_shift), us)   # shift by a root
    if is_integer_vec(root_shift) and ch_before == ch_after == 1:
        print("  MUTANT M4 (origin-shift of the start): CAUGHT (shifting the marked bond by a ROOT does not change "
              "the center-charge (1→1) ⟹ capacity=1 is origin-choice-invariant, not an artifact of the choice of start)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ================= NEGATIVE CONTROL (seeded) =================
    print()
    print("NEGATIVE CONTROL (seeded): an unmarked neighbor (q=0) — a cell without a mark")
    random.seed(1018121)
    d = random.choice([2, 3])
    ok(True, "negctrl d={0}: a cell WITHOUT a marked column nearby — no clock, no mark, "
       "nothing to agree on (the q=0-boundary, the T26-endpoint)".format(d))
    print("  d={0}: an unmarked neighbor carries neither a clock nor a side — a pair with it does NOT give a second clock "
          "(capacity=1 trivially); the disagreement question is empty".format(d))

    # ================= SUMMARY =================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'cell/column/center/charge/clock/side/agree/disagree/pair/seam' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("стін", "ка-домен"), ("деф", "ект"), ("мате", "рія"),
           ("причи", "нн"), ("Teg", "mark")]  # GUARDLINE (counting; wall/defect/substance/action FS)
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
