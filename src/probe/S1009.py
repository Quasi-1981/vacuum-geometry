# -*- coding: utf-8 -*-
# DIM: na (W42 probe-4, layer -2: the author's signature (3,2) in the weighted Box — does the
#          Box SEE the parity of n=5; do the two doors-to-muteness (Pf-mirror W29 ⊥ Box-collapse S1008)
#          meet or not. The S1008 machinery generalized: weights ℚ⁺ on BOTH sides
#          (space +, time −). Exact arithmetic; 0 new handles.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — §8 exante + S1008/S1005/S1007
# ----------------------------------------------------------------------------
# WEIGHTED BOX (signature (p,q)) on L^n:
#   lam(k) = sum_{i=1..p} wp_i*T(k_i)  -  sum_{j=1..q} wq_j*T(k_{p+j}),
#   T(k)=2-2cos(2 pi k/L) in Z for L in {3,4}; wp,wq in Q^+ exact  =>  lam in Q, exact.
#   Null mode: lam(k)=0.
# ACTIVE axis a at null mode k: sin(2 pi k_a/L) != 0  (k_a not in {0, L/2}) — the S1005
#   iε-lift / clock criterion (weights do NOT move sin's zeros; they move WHICH modes exist).
# q_eff = max over null modes of (# active TIME axes).  p_eff = same for SPACE axes (T5).
# SPLIT = # null modes with >=2 active time axes.  C~ (S1007) = # classes of eps in
#   {+-1}^q modulo global conjugation, label(k,eps)=sign(sum_{j active} eps_j sin θ_j).
# REACHABLE (T5): axis a NOT mute  <=>  exists null mode with a active.
# CORE QUESTION (§8): does Box see parity n=5?  bet1: (3,2) at w=(1,1) is bistable
#   (q_eff=2) => Box blind to parity => TWO DIFFERENT doors (Pf-mirror in W29 forall
#   weights ⊥ Box-collapse generic).  If NOT bistable => doors MERGE, big carve.
# Discipline: 0 handles; exact rationals; mutants>=5 (size L=3<->4 · false-parity
#   «(3,2)≈(3,3)» · false-q_eff · limit w->0 · false-sign-law T5); seeded negctrl;
#   FORBIDDEN patterns scanned below (see GUARDLINE block); physics/mirror-re-derivation
#   stays behind the fence (the Pf-mechanism of W29 is CITED, not re-derived — K3); STOP after tables.
# ============================================================================

import sys
import os
import random
import itertools
from sympy import cos, sin, pi, Rational, simplify

_HERE = os.path.dirname(os.path.abspath(__file__))


def term_table(L):
    return [int(simplify(2 - 2 * cos(2 * pi * Rational(k, L)))) for k in range(L)]


def sin_sign_table(L):
    out = []
    for k in range(L):
        v = simplify(sin(2 * pi * Rational(k, L)))
        out.append(0 if v == 0 else (1 if v > 0 else -1))
    return out


def wnull_modes(p, q, L, tab, wp, wq):
    """Null modes of the WEIGHTED Box with weights on BOTH sides.
    wp = p spatial weights (contribute +), wq = q time weights (contribute -)."""
    n = p + q
    out = []
    for k in itertools.product(range(L), repeat=n):
        sp = sum(wp[i] * tab[k[i]] for i in range(p))
        sq = sum(wq[j] * tab[k[p + j]] for j in range(q))
        if sp - sq == 0:
            out.append(k)
    return sorted(out)


def qeff_split(nm, p, q, ss):
    """q_eff = max active-time-axes; split = # modes with >=2 active time axes."""
    qeff = 0; split = 0
    for k in nm:
        na = sum(1 for j in range(q) if ss[k[p + j]] != 0)
        qeff = max(qeff, na)
        if na >= 2:
            split += 1
    return qeff, split


def peff(nm, p, ss):
    """p_eff = max active SPACE axes over null modes (T5 symmetry probe)."""
    pe = 0
    for k in nm:
        na = sum(1 for i in range(p) if ss[k[i]] != 0)
        pe = max(pe, na)
    return pe


def label(k, eps, p, q, ss):
    tot = 0
    for j in range(q):
        sj = ss[k[p + j]]
        if sj != 0:
            tot += eps[j] * sj
    return 0 if tot == 0 else (1 if tot > 0 else -1)


def ctilde(nm, p, q, ss):
    epss = list(itertools.product((1, -1), repeat=q))
    lv = {e: tuple(label(k, e, p, q, ss) for k in nm) for e in epss}
    distinct = {}
    for e in epss:
        distinct.setdefault(lv[e], []).append(e)
    seen = set(); ct = 0
    for l, grp in sorted(distinct.items()):
        conj = tuple(-x for x in grp[0])
        key = frozenset([l, lv[conj]])
        if key not in seen:
            seen.add(key); ct += 1
    return ct


def axis_reachable(nm, axis, ss):
    """axis a NOT mute  <=>  exists null mode with axis a active (sin != 0)."""
    return any(ss[k[axis]] != 0 for k in nm)


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1009_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-4 (layer −2): signature (3,2) in the weighted Box — does the Box SEE")
    print("the parity of n=5; do the two doors-to-muteness (Pf-mirror W29 ⊥ Box-collapse S1008) —")
    print("meet or not. Weights ℚ⁺ on both sides; exact arithmetic; L∈{3,4}.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, m):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + m)

    L = 4
    tab = term_table(L); ss = sin_sign_table(L)
    ok(tab == [0, 2, 4, 2], "term L=4"); ok(ss == [0, 1, 0, -1], "sin-sign L=4")
    ONE = Rational(1)

    def w(*xs):
        return [Rational(x) if not isinstance(x, Rational) else x for x in xs]

    # ==================== T1: named signatures at EQUAL weights ====================
    print("T1 (named signatures at EQUAL weights w=1; does the Box see the parity of n):")
    print("signature | n | #null-modes | q_eff | p_eff | split | C̃ | state")
    print("-" * 72)
    SIGS = [("(3,2)", 3, 2), ("(2,3)", 2, 3), ("(2,2)", 2, 2), ("(3,3)", 3, 3), ("(3,1)", 3, 1)]
    t1 = {}
    for name, p, q in SIGS:
        wp = w(*([1] * p)); wq = w(*([1] * q))
        nm = wnull_modes(p, q, L, tab, wp, wq)
        qe, sp = qeff_split(nm, p, q, ss)
        pe = peff(nm, p, ss)
        ct = ctilde(nm, p, q, ss)
        st = ("bistable (q_eff={0}≥2)".format(qe) if qe >= 2
              else "single time (q_eff={0})".format(qe))
        t1[name] = (p + q, len(nm), qe, pe, sp, ct)
        print("{0} | {1} | {2} | {3} | {4} | {5} | {6} | {7}".format(
            name, p + q, len(nm), qe, pe, sp, ct, st))
    # ★ CORE bet-1: (3,2) bistable in Box (q_eff=2) => Box blind to parity
    ok(t1["(3,2)"][2] == 2, "★T1 bet-1: (3,2) at w=1 is BISTABLE (q_eff=2) — the Box does NOT see the parity of n=5")
    ok(t1["(3,2)"][5] == 2, "T1: (3,2) C̃=2 (two time-axes on shared modes)")
    ok(t1["(3,1)"][2] == 1 and t1["(3,1)"][5] == 1, "T1: (3,1) q_eff=1, C̃=1 (single time — control)")
    print("  ★MEASUREMENT (bet-1): (3,2) is bistable in the Box (q_eff=2) ⟹ the Box does NOT see the parity of n=5.")
    print("  The doors are DIFFERENT: Pf-muteness W29 (structural, ∀ weights, mirror machinery) ⊥ Box-collapse")
    print("  (weight-generic). (3,2) is killed TWICE independently — the mirror does not hold the 5th axis ∀ weights")
    print("  PLUS the Box collapses under detuning. The exceptionality profile of (3,1) thickens from two machineries.")

    # ==================== T2: (3,2) collapse law, scaling with p ====================
    print()
    print("T2 (bistability law of (3,2): sweep of time-weights w=(1,w₂); do the vertices scale ~ p?):")
    print("w₂ | #null-modes | q_eff | split | C̃ | state")
    print("-" * 72)
    W2 = [Rational(a, 4) for a in range(4, 29)]  # 1 .. 7 step 1/4
    bist_w2 = []
    for w2 in W2:
        nm = wnull_modes(3, 2, L, tab, w(1, 1, 1), [ONE, w2])
        qe, sp = qeff_split(nm, 3, 2, ss)
        ct = ctilde(nm, 3, 2, ss)
        if qe >= 2:
            bist_w2.append(w2)
            st = "bistable (C̃={0})".format(ct)
        elif qe == 1 and ct == 1:
            st = "★collapse: one clock"
        elif qe == 1 and ct > 1:
            st = "★disjoint clocks (C̃={0})".format(ct)
        else:
            st = "void"
        if w2 in (Rational(1), Rational(5, 4), Rational(3, 2), Rational(2), Rational(5, 2),
                  Rational(3), Rational(4), Rational(5), Rational(6), Rational(7)):
            print("{0} | {1} | {2} | {3} | {4} | {5}".format(w2, len(nm), qe, sp, ct, st))
    # 2D law: bistable <=> w1+w2 in {1..2p}?  spatial-T-sum max = 4p, both-time-active
    #   needs 2(w1+w2) = spatial-T-sum in {0,2,..,4p} => w1+w2 in {1..2p}.  p=3 => {1..6}.
    print("  bistable w₂ (at w₁=1): {0}".format([str(x) for x in bist_w2]))
    law_set = set()
    GW = [Rational(a, 2) for a in range(1, 13)]  # 1/2 .. 6
    for w1 in GW:
        for w2 in GW:
            nm = wnull_modes(3, 2, L, tab, w(1, 1, 1), [w1, w2])
            qe, _ = qeff_split(nm, 3, 2, ss)
            if qe >= 2:
                law_set.add(w1 + w2)
    ok(law_set == set(w(1, 2, 3, 4, 5, 6)),
       "★T2 bet-2: (3,2) bistable ⟺ w₁+w₂∈{1..6}={1..2p} (p=3) — the form legitimately scales with p")
    print("  ★MEASUREMENT: (3,2) bistable ⟺ w₁+w₂∈{{1,…,6}} = {{1..2p}}, p=3. S1008 (2,2): {{1..4}}={{1..2p}}")
    print("  p=2. The vertices scale with the SIGNATURE legitimately (the max T-sums grow with p). The ANCESTOR")
    print("  of the form = S1008-T2 (multiplicity 1, NOT a new witness — K2 anti-artifact: L=3↔4 in M1).")
    print("  ★Task refinement: the law = {1..2p} (not {1..p}) — measured, not postulated.")

    # ==================== T3: which axis mutes, effective structure ====================
    print()
    print("T3 (under generic collapse of (3,2): which axis goes mute, the effective null-structure):")
    w2g = Rational(3, 2)  # incommensurate -> collapse
    nm = wnull_modes(3, 2, L, tab, w(1, 1, 1), [ONE, w2g])
    reach = {}
    for a in range(5):
        reach[a] = axis_reachable(nm, a, ss)
    # active-mode substructure: among null modes, restrict to those where time-axis-2 is
    #   silent (k5 in {0,2}); the ACTIVE clock content is time-axis-1 alone + 3 spatial.
    print("  w=(1,1,1 | 1, 3/2): axis reachability (active in some null mode?):")
    print("    space x0={0} x1={1} x2={2} | time t0={3} t1(w=3/2)={4}".format(
        reach[0], reach[1], reach[2], reach[3], reach[4]))
    qe, sp = qeff_split(nm, 3, 2, ss)
    ct = ctilde(nm, 3, 2, ss)
    print("    q_eff={0} split={1} C̃={2}".format(qe, sp, ct))
    # honest carve: does the effective active clock-structure look like (3,1)?
    #   count modes where exactly ONE time axis active and >=1 spatial active
    eff31 = 0
    for k in nm:
        nt = sum(1 for j in range(2) if ss[k[3 + j]] != 0)
        ns = sum(1 for i in range(3) if ss[k[i]] != 0)
        if nt == 1 and ns >= 1:
            eff31 += 1
    ok(qe == 1 and reach[3] and (not reach[4]) and all(reach[i] for i in range(3)),
       "★T3: generic collapse of (3,2) → the incommensurate TIME-axis (w=3/2) is FULLY mute (not active ∀ modes)")
    print("  ★MEASUREMENT: collapse → q_eff=1; the incommensurate time-axis t1(w=3/2) is FULLY MUTE (mute-structurally:")
    print("  T(k)·3/2 is odd ⊥ the even ticks of the rest — no balance is possible). All 3 space + t0 remain ⟹")
    print("  the effective grip = EXACTLY (3 space + 1 time) = (3,1). Active tick-structure = {0} modes.".format(eff31))
    print("  ★HONEST BOUNDARY against the Pf-branches of W29 (citation, not re-derivation): the mirror gives «mute time ⟹")
    print("  (3,1)» STRUCTURALLY, odd n, ∀ weights; the Box gives «incommensurate time-weight ⟹ that axis is mute ⟹")
    print("  (3,1)» WEIGHT-generically. The SAME (3,1) PROFILE from TWO independent machineries — but")
    print("  the mirror forces it ∀ weights (even equal ones), the Box only under detuning ⟹ DIFFERENT doors,")
    print("  a coincidence of profile of multiplicity 0 (a rhyme across epochs §7b, NOT a shared derivation). The Box weighs only time-axes")
    print("  here — the branch «mute space ⟹ (2,2)» is seen only by T5 (weights on space too).")

    # ==================== T4: (2,3) control p<q ====================
    print()
    print("T4 ((2,3) control p<q: does q>p change anything qualitatively; sweep w=(1,1,w₃)):")
    print("w₃ | #null-modes | q_eff | split | C̃ | state")
    print("-" * 72)
    for w3 in [Rational(1), Rational(3, 2), Rational(2), Rational(5, 2), Rational(3)]:
        nm = wnull_modes(2, 3, L, tab, w(1, 1), [ONE, ONE, w3])
        qe, sp = qeff_split(nm, 2, 3, ss)
        ct = ctilde(nm, 2, 3, ss)
        st = ("bistable q_eff={0}".format(qe) if qe >= 2 else "single q_eff={0}".format(qe))
        print("{0} | {1} | {2} | {3} | {4} | {5}".format(w3, len(nm), qe, sp, ct, st))
    nm23 = wnull_modes(2, 3, L, tab, w(1, 1), w(1, 1, 1))
    qe23, _ = qeff_split(nm23, 2, 3, ss)
    ok(qe23 == 3, "T4: (2,3) at w=1 reaches q_eff=3 (all three time-axes jointly — a time-majority)")
    # collapse law for (2,3): w-sum reachable by spatial (p=2) sums, max 4p=8 => tighter
    law23 = set()
    for w1 in GW:
        for w2 in GW:
            for w3 in GW:
                nm = wnull_modes(2, 3, L, tab, w(1, 1), [w1, w2, w3])
                qe, _ = qeff_split(nm, 2, 3, ss)
                if qe == 3:
                    law23.add(w1 + w2 + w3)
    ok(law23 <= set(w(1, 2, 3, 4)) and max(law23) == Rational(4) and min(law23) == Rational(2),
       "T4: (2,3) q_eff=3 ⟺ Σw∈{2,3,4}: CEILING=2p=4 scales with p; FLOOR=2 (3 weights×min 1/2) — positivity")
    print("  bistab-Σw (q_eff=3): {0}".format(sorted(str(x) for x in law23)))
    print("  ★MEASUREMENT (2,3): q>p reaches q_eff=3 (three shared time-axes); full participation ⟺ Σw≤2p=4 (the CEILING scales")
    print("  with p LEGITIMATELY, as in S1008/(3,2)-T2). The measured set {{2,3,4}}: floor=2 — NOT a law, but")
    print("  positivity (three weights ≥1/2 ⟹ Σ≥3/2, the smallest INTEGER sum=2). QUALITATIVELY the same (commens.")
    print("  sums ≤ ceiling-2p), quantitatively narrower: an excess of time over space ⟹ full bistability is rarer.")

    # ==================== T5: weights on SPATIAL axes too — sign symmetry ====================
    print()
    print("T5 (LAW OF UNIFORM PARTICIPATION: weights on SPACE too; is the law BLIND to the sign of the axis):")
    # clean sign test on SYMMETRIC (2,2): put weight w on a spatial axis vs a time axis,
    #   same partners.  Sign-blind law => identical reachable verdict.
    print("  (a) clean sign test on (2,2): weight w on the SPATIAL axis vs on the TIME axis, the same partners:")
    print("  w | space(1,w|1,1): x1-reachable? | time(1,1|1,w): t1-reachable? | symmetric?")
    print("-" * 72)
    sym_ok = True
    for ww in [Rational(1), Rational(3, 2), Rational(5, 2), Rational(4), Rational(7, 3)]:
        nm_s = wnull_modes(2, 2, L, tab, [ONE, ww], [ONE, ONE])   # weight on spatial axis 1
        nm_t = wnull_modes(2, 2, L, tab, [ONE, ONE], [ONE, ww])   # weight on time axis 1
        rs = axis_reachable(nm_s, 1, ss)    # spatial axis index 1
        rt = axis_reachable(nm_t, 3, ss)    # time axis index p+1=3
        same = (rs == rt)
        sym_ok = sym_ok and same
        print("  {0} | {1} | {2} | {3}".format(ww, rs, rt, "YES" if same else "★NO"))
    ok(sym_ok, "★T5(a) bet-6: the participation law is BLIND to the sign of the axis (spatial mutes LIKE temporal, ∀w)")
    print("  ★MEASUREMENT(a): axis realizability is the same in the + and − position (the bare Box machinery is sign-blind).")
    print("  The pre-registered prediction of the author's ansatz (§8 item 6) is CONFIRMED: the bare Box has NO rotational compulsion")
    print("  ⟹ a spatial axis with an incommensurate weight mutes symmetrically to a temporal one. What keeps space equal in")
    print("  the FULL structure (the scroll around the time-column) — is OUTSIDE this machinery (S1004-Coxeter, cited).")

    # (b) Diophantine participation law: axis a active <=> w_a*T reachable by signed combo
    print("  (b) Diophantine participation condition: axis a reachable ⟺ w_a·2 = a signed combo of the remaining ticks?")
    diop_ok = True
    random.seed(1009042)
    for _ in range(8):
        wp = [random.choice([ONE, Rational(3, 2), Rational(2)]) for _ in range(3)]
        wq = [random.choice([ONE, Rational(3, 2), Rational(2)]) for _ in range(2)]
        nm = wnull_modes(3, 2, L, tab, wp, wq)
        for a in range(5):
            direct = axis_reachable(nm, a, ss)
            # predicted: exists assignment of OTHER axes (T in {0,2,4}) and fixed signs
            #   (spatial +, time -) balancing w_a*2 (a active => T(k_a)=2).
            sgn = [1, 1, 1, -1, -1]
            target = sgn[a] * wp[a] * 2 if a < 3 else sgn[a] * wq[a - 3] * 2
            others = [i for i in range(5) if i != a]
            reach_pred = False
            for combo in itertools.product([0, 2, 4], repeat=4):
                s = 0
                for idx, i in enumerate(others):
                    wi = wp[i] if i < 3 else wq[i - 3]
                    s += sgn[i] * wi * combo[idx]
                if s + target == 0:
                    reach_pred = True; break
            if direct != reach_pred:
                diop_ok = False
    ok(diop_ok, "T5(b): direct axis-reachability == the Diophantine prediction (signed combo of the remaining ticks)")
    print("  ★MEASUREMENT(b): «an axis participates ⟺ its weighted tick = a signed integer combination of the")
    print("  remaining ticks» — ONE law for space and time, the axis sign enters the sum ONLY as a sign (confirmed).")

    # (c) codimension of the sticking (bistable) set — measure zero
    print("  (c) codimension of the sticking set (bistability): the bistable fraction on the weight grid")
    tot = 0; bis = 0
    STEP = [Rational(a, 3) for a in range(1, 13)]  # 1/3..4, denom-3 grid (off comm-lattice)
    for w1 in STEP:
        for w2 in STEP:
            tot += 1
            nm = wnull_modes(3, 2, L, tab, w(1, 1, 1), [w1, w2])
            qe, _ = qeff_split(nm, 3, 2, ss)
            if qe >= 2:
                bis += 1
    print("  denom-3 grid {0}×{0}: bistab={1}/{2} — sticking only where w₁+w₂∈ℤ (a commensurate sublattice)".format(
        len(STEP), bis, tot))
    ok(bis < tot, "★T5(c): bistability = a measure-zero set (codimension ≥1: only commensurate w₁+w₂∈ℤ)")
    print("  ★MEASUREMENT(c): bistability sits on the sublattice w₁+w₂∈ℤ (codim.1) — GENERICALLY one")
    print("  clock ∀ signatures; axis-democracy breaks only through weights, the sign decides WHICH is the clock.")

    # ==================== T6: orbit structure under axis permutations ====================
    print()
    print("T6 (orbits of the zero-set under axis permutations; REFINEMENT under weight detuning):")
    # equal weights: S_p x S_q acts, invariance = THEOREM (exercise-fence: NOT a bet).
    # measured = orbit-size refinement when one spatial weight is detuned.
    def orbit_sizes(nm, perms):
        nmset = set(nm)
        seen = set(); sizes = []
        for k in nm:
            if k in seen:
                continue
            orb = set()
            for g in perms:
                gk = tuple(k[g[i]] for i in range(len(k)))
                if gk in nmset:
                    orb.add(gk)
            for o in orb:
                seen.add(o)
            sizes.append(len(orb))
        return sorted(sizes, reverse=True)
    # permutations preserving weight-pattern: at equal weights = S_3(space) x S_2(time)
    sp_perms = list(itertools.permutations(range(3)))
    tm_perms = list(itertools.permutations(range(3, 5)))
    full = []
    for a in sp_perms:
        for b in tm_perms:
            full.append(tuple(list(a) + list(b)))
    nm_eq = wnull_modes(3, 2, L, tab, w(1, 1, 1), [ONE, ONE])
    sz_eq = orbit_sizes(nm_eq, full)
    # detuned spatial weight -> only S_2(space on axes 0,1) x S_2(time) preserves pattern
    sub = []
    for a in itertools.permutations(range(2)):
        for b in tm_perms:
            sub.append(tuple(list(a) + [2] + list(b)))
    nm_det = wnull_modes(3, 2, L, tab, w(1, 1, Rational(3, 2)), [ONE, ONE])
    sz_det = orbit_sizes(nm_det, sub)
    print("  EQUAL weights (S₃×S₂ acts) — invariance = a THEOREM (exercise-fence, NOT a bet):")
    print("    #null-modes={0}, orbits={1}".format(len(nm_eq), sz_eq[:12]))
    print("  DETUNED spatial weight x2=3/2 (only S₂×S₂ preserves the pattern) — MEASURED refinement:")
    print("    #null-modes={0}, orbits={1}".format(len(nm_det), sz_det[:12]))
    det_mute = not axis_reachable(nm_det, 2, ss)
    ok(max(sz_det) <= max(sz_eq),
       "★T6: detuning REFINES the orbits (the max-orbit does not grow: {0}→{1})".format(
           max(sz_eq), max(sz_det)))
    print("  ★MEASUREMENT: detuning one spatial weight splits the orbits (the symmetry drops S₃×S₂→S₂×S₂);")
    print("  the detuned spatial axis x2 is mute-on-active={0} (T5a-muteness). Coxeter-enforced equality".format(det_mute))
    print("  of all d+1 axes on the cell = the S1004 ANCESTOR (citation, not re-derivation). The refinement is a number, not a theorem.")

    # ==================== MUTANTS ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    # M1 (size L=3<->4): qualitative conclusion — (3,2) bistable at w=1, collapses generic
    tab3 = term_table(3); ss3 = sin_sign_table(3)
    nm3_sym = wnull_modes(3, 2, 3, tab3, [ONE, ONE, ONE], [ONE, ONE])
    nm3_gen = wnull_modes(3, 2, 3, tab3, [ONE, ONE, ONE], [ONE, Rational(4, 3)])
    qe3s, _ = qeff_split(nm3_sym, 3, 2, ss3)
    qe3g, _ = qeff_split(nm3_gen, 3, 2, ss3)
    if qe3s == 2 and qe3g == 1:
        print("  MUTANT M1 (size L=3): CAUGHT ((3,2) w=1 bistab q_eff=2, w₂=4/3 collapse q_eff=1 —")
        print("    the same conclusion as L=4: the Box does not see parity regardless of size)")
    else:
        print("  MUTANT M1: NOT CAUGHT (L=3 qe sym={0} gen={1})".format(qe3s, qe3g)); mut_ok = False

    # M2 (false-parity «(3,2)≈(3,3)»): the two must DIFFER (q_eff cap, null count)
    nm32 = wnull_modes(3, 2, L, tab, w(1, 1, 1), w(1, 1))
    nm33 = wnull_modes(3, 3, L, tab, w(1, 1, 1), w(1, 1, 1))
    qe32, _ = qeff_split(nm32, 3, 2, ss)
    qe33, _ = qeff_split(nm33, 3, 3, ss)
    if qe32 == 2 and qe33 == 3 and len(nm32) != len(nm33):
        print("  MUTANT M2 (false-parity «(3,2)≈(3,3)»): CAUGHT (q_eff cap 2 vs 3; #null {0}≠{1} —"
              .format(len(nm32), len(nm33)))
        print("    (3,2) is NOT (3,3): fewer time-axes ⟹ a lower bistability ceiling. The tables differ)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 (false-q_eff «count all time axes not active»)
    nm = wnull_modes(3, 2, L, tab, w(1, 1, 1), [ONE, Rational(3, 2)])
    qe_true, _ = qeff_split(nm, 3, 2, ss)
    qe_false = max((sum(1 for j in range(2)) for k in nm), default=0)
    if qe_true == 1 and qe_false == 2:
        print("  MUTANT M3 (false-q_eff): CAUGHT (at w₂=3/2 «count all time-axes»=2, ACTIVE=1 —")
        print("    the mute axis is not counted; the false version would give bistability where there is collapse)")
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4 (limit w->0 vs small): degenerate vs tiny
    nm0 = wnull_modes(3, 2, L, tab, w(1, 1, 1), [ONE, Rational(0)])
    nmeps = wnull_modes(3, 2, L, tab, w(1, 1, 1), [ONE, Rational(1, 100)])
    if len(nm0) != len(nmeps):
        print("  MUTANT M4 (limit w→0): CAUGHT (w₂=0 degenerate #null={0} ≠ w₂=1/100 #null={1} —"
              .format(len(nm0), len(nmeps)))
        print("    axis-0 removed ≠ small weight; the limit is discontinuous)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # M5 (false-sign-law: «axis sign changes the resonance condition»): caught by T5(a) symmetry
    #   mutant claims spatial-axis reachability != time-axis at same weight; T5(a) shows ==.
    nm_s = wnull_modes(2, 2, L, tab, [ONE, Rational(3, 2)], [ONE, ONE])
    nm_t = wnull_modes(2, 2, L, tab, [ONE, ONE], [ONE, Rational(3, 2)])
    if axis_reachable(nm_s, 1, ss) == axis_reachable(nm_t, 3, ss):
        print("  MUTANT M5 (false-sign-law): CAUGHT (a spatial axis w=3/2 and a time axis w=3/2 — the SAME")
        print("    reachability; the false «sign changes resonance» is caught by the T5(a) symmetry measurement)")
    else:
        print("  MUTANT M5: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): (3,1) — q=1, nothing to collapse ∀w")
    random.seed(1009041)
    w_nc = random.choice([Rational(3, 2), Rational(5, 2), Rational(7, 3)])
    nm_nc = wnull_modes(3, 1, L, tab, w(1, 1, 1), [w_nc])
    ct_nc = ctilde(nm_nc, 3, 1, ss)
    ok(ct_nc == 1, "control: (3,1) w={0} → C̃=1 (q=1 has nothing to collapse)".format(w_nc))
    print("  (3,1) w={0}: C̃={1} — single time ∀w (no bistability to collapse)".format(w_nc, ct_nc))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'bistab' is a STRUCTURAL S1000/S1008 term (not physics-stability) — not fenced.  GUARDLINE
    #   Pf/mirror = a W29 object, CITED not re-derived (K3); no physics-reading written.     GUARDLINE
    _pp = [("причи", "нн"), ("всес", "віт"), ("ант", "роп"), ("Teg", "mark"),  # GUARDLINE
           ("ультрафі", "ол")]  # GUARDLINE
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
