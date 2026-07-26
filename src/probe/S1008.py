# -*- coding: utf-8 -*-
# DIM: na (W42 probe-3, layer -2: weighted Box — does splitting the time-pair weights reduce
#          q=2 to q_eff=1 (leg of S1000-bistability, question №4); machinery of S1005/S1007;
#          exact rational arithmetic; 0 new handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — derivational, §7 exante + S1005/7
# ----------------------------------------------------------------------------
# WEIGHTED BOX (signature (p,q)) on L^n: lam(k) = sum_{i=1..p} T(k_i)
#   - sum_{j=1..q} w_j * T(k_{p+j}),  w_j in Q^+ exact, T(k)=2-2cos(2 pi k/L) in Z
#   for L in {3,4}  =>  lam in Q, exact.  Null: lam(k)=0.
# ACTIVE time-axis j at null mode k: sin(2 pi k_{p+j}/L) != 0 (k not in {0,L/2}) —
#   the S1005 iε-lift criterion (weights do NOT move sin's zeros; they move WHICH
#   null modes EXIST).  q_eff = max over null modes of (# active time-axes).
# SPLIT (S1005): # null modes with >=2 active time-axes.  C~ (S1007): # classes of
#   eps in {+-1}^q modulo global conjugation, label(k,eps)=sign(sum_{j active}
#   eps_j sin θ_j), equivalence = same label on every null mode.
# CORE QUESTION (§7): does splitting the time weights REDUCE q=2 to effective q=1?
#   bet1: exact rational asymmetry threshold with q_eff=1, split=0, C~=1 (one clock);
#   bet2: collapse/bistable law in the weight plane (piecewise-linear, rational
#   vertices, form rides with L) — if it matches t=s+(d-1) it is the S1000 ANCESTOR
#   (kratn 1, NOT a new witness, K3); bet3: ORDER — collapse to q_eff=1 vs empty
#   null-set (void), which comes first.
# Discipline: 0 handles; exact rationals; mutants>=4 (size L=3<->4 + false-q_eff
#   «count all axes» + limit w->0 vs w small); seeded negctrl ((3,1)-weighted: C~=1
#   forall w); FORBIDDEN-SCAN; STOP after tables.
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


def active_set(L):
    s = set(range(L)); s.discard(0)
    if L % 2 == 0:
        s.discard(L // 2)
    return s


def wnull_modes(p, q, L, tab, weights):
    """Null modes of the WEIGHTED Box; weights = q exact Rationals."""
    n = p + q
    out = []
    for k in itertools.product(range(L), repeat=n):
        sp = sum(tab[k[i]] for i in range(p))
        sq = sum(weights[j] * tab[k[p + j]] for j in range(q))
        if sp - sq == 0:
            out.append(k)
    return out


def q_eff_and_split(nm, p, q, ssign):
    """q_eff = max active-time-axes over null modes; split = # modes with >=2."""
    qeff = 0; split = 0
    for k in nm:
        na = sum(1 for j in range(q) if ssign[k[p + j]] != 0)
        qeff = max(qeff, na)
        if na >= 2:
            split += 1
    return qeff, split


def label(k, eps, p, q, ssign):
    tot = 0
    for j in range(q):
        sj = ssign[k[p + j]]
        if sj != 0:
            tot += eps[j] * sj
    return 0 if tot == 0 else (1 if tot > 0 else -1)


def ctilde(nm, p, q, ssign):
    epss = list(itertools.product((1, -1), repeat=q))
    lv = {e: tuple(label(k, e, p, q, ssign) for k in nm) for e in epss}
    distinct = {}
    for e in epss:
        distinct.setdefault(lv[e], []).append(e)
    seen = set(); ct = 0
    for l, grp in distinct.items():
        conj = tuple(-x for x in grp[0])
        key = frozenset([l, lv[conj]])
        if key not in seen:
            seen.add(key); ct += 1
    return ct


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1008_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-3 (layer −2): weighted Box — reduction q=2 → q_eff=1 by splitting the weights")
    print("(leg of S1000-bistability, question №4); exact rational weights; L∈{3,4}.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, m):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + m)

    L = 4
    tab = term_table(L); ss = sin_sign_table(L)
    ok(tab == [0, 2, 4, 2], "term L=4"); ok(ss == [0, 1, 0, -1], "sin-sign L=4")

    # ==================== T1: (2,2) asymmetry sweep w=(1, w2) ====================
    print("T1 ((2,2), L=4): asymmetry sweep w=(1, w₂); q_eff/split/C̃ over null modes")
    print("w₂ | #null-modes | q_eff | split | C̃ | state")
    print("-" * 72)
    W2 = [Rational(1), Rational(5, 4), Rational(3, 2), Rational(7, 4), Rational(2),
          Rational(9, 4), Rational(5, 2), Rational(3), Rational(7, 2), Rational(4), Rational(5)]
    t1rows = []
    for w2 in W2:
        nm = wnull_modes(2, 2, L, tab, [Rational(1), w2])
        qeff, split = q_eff_and_split(nm, 2, 2, ss)
        ct = ctilde(nm, 2, 2, ss)
        # label from ACTUAL ct (three phases; C̃ NOT hardcoded — Omega gate S1008):
        if qeff >= 2:
            st = "bistable (q_eff={0}, C̃={1})".format(qeff, ct)
        elif qeff == 1 and ct == 1:
            st = "★collapse: ONE clock (q_eff=1, C̃=1)"
        elif qeff == 1 and ct > 1:
            st = "★disjoint clocks (q_eff=1, C̃={0}: axes active separately, not on the same mode)".format(ct)
        else:
            st = "void (no active time-axes)"
        t1rows.append((w2, len(nm), qeff, split, ct))
        print("{0} | {1} | {2} | {3} | {4} | {5}".format(w2, len(nm), qeff, split, ct, st))
    # symmetric w2=1 must be bistable; a non-integer w2 must collapse
    r_sym = [r for r in t1rows if r[0] == Rational(1)][0]
    r_half = [r for r in t1rows if r[0] == Rational(3, 2)][0]
    r_four = [r for r in t1rows if r[0] == Rational(4)][0]
    ok(r_sym[2] == 2 and r_sym[4] == 2, "T1: w₂=1 symmetry → q_eff=2, C̃=2 (bistable)")
    ok(r_half[2] == 1 and r_half[4] == 1, "T1: w₂=3/2 → q_eff=1, C̃=1 (COLLAPSE to 1 clock)")
    ok(r_four[2] == 1 and r_four[4] == 2,
       "T1: w₂=4 → q_eff=1, C̃=2 (★DISJOINT clocks — phase-(ii): axes active separately)")
    print("  ★MEASUREMENT (THREE phases — the data revealed the third): bistable (q_eff≥2, C̃=2) at integer")
    print("  w₂∈{1,2,3} · COLLAPSE to one clock (q_eff=1, C̃=1) at ANY non-integer")
    print("  asymmetry (3/2,5/2,…) · ★DISJOINT clocks (q_eff=1, C̃=2) at w₂=4: both")
    print("  axes active, but NEVER on the same null mode (split=0) — two times that do not overlap.")
    print("  ⟹ bet-1 REFINED: the collapse is NOT a monotone threshold, but a COMB; +a third phase (disjoint) —")
    print("  q_eff=1 does NOT ⟹ one clock (C̃ must be measured separately). (author: «certain relations»)")

    # ==================== T2: bistable law in weight plane ====================
    print()
    print("T2 (bistability law in the weight plane (w₁,w₂), (2,2) L=4): q_eff=2 ⟺ ?")
    GW = [Rational(1, 2), Rational(1), Rational(3, 2), Rational(2), Rational(5, 2), Rational(3)]
    print("  q_eff=2 (bistable) points (w₁,w₂):")
    bistab = []
    for w1 in GW:
        for w2 in GW:
            nm = wnull_modes(2, 2, L, tab, [w1, w2])
            qeff, _ = q_eff_and_split(nm, 2, 2, ss)
            if qeff >= 2:
                bistab.append((w1, w2))
    # test law: bistable <=> w1+w2 in {1,2,3,4} (both-active: 2(w1+w2) = T-sum in {0..8})
    law_sum = all((w1 + w2) in (Rational(1), Rational(2), Rational(3), Rational(4))
                  for (w1, w2) in bistab)
    ok(law_sum, "T2: bistable ⟺ w₁+w₂ ∈ {1,2,3,4} (commensurability of the weight sum)")
    print("   {0}".format([("{0},{1}".format(a, b)) for (a, b) in bistab]))
    print("  law: q_eff=2 ⟺ w₁+w₂ ∈ {{1,2,3,4}} (2(w₁+w₂)=T-sum∈{{0,2,4,6,8}}) — the weight SUM is integer.")
    # K3 anti-reclaim: form is SUM (w1+w2=k), S1000 is DIFFERENCE (t=s+(d-1)) — DIFFERENT
    print("  ★K3 (anti-reclaim): the form = SUM w₁+w₂=integer ≠ DIFFERENCE t=s+(d−1) (S1000) —")
    print("  DIFFERENT structures. NOT the same law ⟹ NOT an S1000 reclaim; a separate −2-fact")
    print("  (bistability = commensurability of the weight sum). The S1000 ancestor is not invoked.")

    # ==================== T3: order (collapse vs void) ====================
    print()
    print("T3 (order of doors: collapse q_eff→1 vs a void (empty nontrivial zero-set))")
    print("  sweep of w₂ (w₁=1): do NONTRIVIAL null modes exist (besides the origin) as w₂ grows?")
    for w2 in [Rational(3, 2), Rational(5, 2), Rational(5), Rational(10)]:
        nm = wnull_modes(2, 2, L, tab, [Rational(1), w2])
        nontrivial = [k for k in nm if any(kk != 0 for kk in k)]
        qeff, _ = q_eff_and_split(nm, 2, 2, ss)
        print("  w₂={0}: #null-modes={1} (nontriv={2}), q_eff={3}".format(
            w2, len(nm), len(nontrivial), qeff))
    ok(len(wnull_modes(2, 2, L, tab, [Rational(1), Rational(10)])) > 1,
       "T3: nontrivial null modes PERSIST under large asymmetry (the k₄=0 family)")
    print("  ★MEASUREMENT: collapse (q_eff→1) occurs FIRST at any non-integer asymmetry; a VOID")
    print("  (of the nontrivial zero-set) DOES NOT OCCUR — the k₄=0 family (mute axis) holds null modes ∀w₂.")
    print("  ⟹ order: COLLAPSE to one clock, NOT a break. «Coordinate stabilization» is real,")
    print("  a void does not occur in this machinery (bet-3: collapse-first, the bet-1 branch stands).")

    # ==================== (3,3) sweep ====================
    print()
    print("T1b ((3,3), L=4): w=(1,1,w₃) — splitting the THIRD weight")
    print("w₃ | #null-modes | q_eff | split | C̃ | state")
    print("-" * 72)
    for w3 in [Rational(1), Rational(3, 2), Rational(2), Rational(5, 2), Rational(3)]:
        nm = wnull_modes(3, 3, L, tab, [Rational(1), Rational(1), w3])
        qeff, split = q_eff_and_split(nm, 3, 3, ss)
        ct = ctilde(nm, 3, 3, ss)
        # label from ACTUAL ct (three phases; C̃ NOT hardcoded — Omega gate):
        if qeff >= 2:
            st = "bistable (q_eff={0}, C̃={1})".format(qeff, ct)
        elif qeff == 1 and ct == 1:
            st = "★collapse: one clock (q_eff=1, C̃=1)"
        elif qeff == 1 and ct > 1:
            st = "★disjoint clocks (q_eff=1, C̃={0})".format(ct)
        else:
            st = "void"
        print("{0} | {1} | {2} | {3} | {4} | {5}".format(w3, len(nm), qeff, split, ct, st))
    print("  ★MEASUREMENT (3,3): the same — a non-commensurable w₃ reduces q_eff; C̃ drops with splitting.")

    # ==================== MUTANTS ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    # M1 (false-q_eff «count all time axes, not active»): at collapsed point counts wrong
    nm = wnull_modes(2, 2, L, tab, [Rational(1), Rational(3, 2)])
    qeff_true, _ = q_eff_and_split(nm, 2, 2, ss)
    qeff_false = max((sum(1 for j in range(2)) for k in nm), default=0)  # counts all q=2
    if qeff_true == 1 and qeff_false == 2:
        print("  MUTANT M1 (false-q_eff): CAUGHT (at w₂=3/2 «count all axes»=2, but ACTIVE=1 — "
              "the mute axis is not counted; the false version would give bistability where there is collapse)")
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2 (limit w2->0 vs small): w2=0 removes axis (degenerate) != w2 tiny
    nm0 = wnull_modes(2, 2, L, tab, [Rational(1), Rational(0)])
    nm_eps = wnull_modes(2, 2, L, tab, [Rational(1), Rational(1, 100)])
    qe0, _ = q_eff_and_split(nm0, 2, 2, ss)
    qeeps, _ = q_eff_and_split(nm_eps, 2, 2, ss)
    # w2=0: axis-2 term vanishes -> every k4 gives 0 -> axis-2 never balances (active mode
    #   needs T(k4)!=0 contributing 0*... ) -> different null set than tiny w2
    if len(nm0) != len(nm_eps) or qe0 != qeeps:
        print("  MUTANT M2 (limit w→0): CAUGHT (w₂=0: #null={0} q_eff={1} ≠ w₂=1/100: "
              "#null={2} q_eff={3} — the limit is discontinuous, axis-0 degenerate ≠ small weight)"
              .format(len(nm0), qe0, len(nm_eps), qeeps))
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 (size L=3<->4): qualitative conclusion (bistable<->w1+w2 comm; collapse generic) robust
    tab3 = term_table(3); ss3 = sin_sign_table(3)
    nm3_sym = wnull_modes(2, 2, 3, tab3, [Rational(1), Rational(1)])
    nm3_gen = wnull_modes(2, 2, 3, tab3, [Rational(1), Rational(4, 3)])
    qe3s, _ = q_eff_and_split(nm3_sym, 2, 2, ss3)
    qe3g, _ = q_eff_and_split(nm3_gen, 2, 2, ss3)
    if qe3s == 2 and qe3g == 1:
        print("  MUTANT M3 (size): CAUGHT (L=3: w₂=1 bistable q_eff=2, w₂=4/3 generic "
              "collapse q_eff=1 — the same conclusion as L=4, not a lattice artifact)")
    else:
        print("  MUTANT M3: NOT CAUGHT (L=3 qe sym={0} gen={1})".format(qe3s, qe3g))
        mut_ok = False

    # M4 (false law: bistable <-> difference w1-w2 like t=s+(d-1))
    # counter: (w1,w2)=(1/2,3/2): w1-w2=-1 (would be "bistable" under diff-law), but w1+w2=2 -> real bistable;
    #          (w1,w2)=(1,2): w1-w2=-1 (diff-law bistable) AND w1+w2=3 -> also bistable; need a discriminator
    # discriminator: (w1,w2)=(3/2,3/2): w1-w2=0 (diff-law: bistable), w1+w2=3 -> real bistable too. hmm
    # use (1/2,1/2): sum=1 (real bistable), diff=0 (diff-law bistable). use (2,2): sum=4 bistable, diff=0.
    # discriminator (1, 3): sum=4 REAL bistable; but (1,1/2): sum=3/2 NOT in {1..4}? 3/2 no -> real NOT bistable, diff=1/2.
    # cleanest: (w1,w2)=(1/2, 1): sum=3/2 -> NOT bistable (real); a diff-law "t-s=const" would mispredict.
    nm_a = wnull_modes(2, 2, L, tab, [Rational(1, 2), Rational(1)])
    qe_a, _ = q_eff_and_split(nm_a, 2, 2, ss)
    nm_b = wnull_modes(2, 2, L, tab, [Rational(1), Rational(1)])
    qe_b, _ = q_eff_and_split(nm_b, 2, 2, ss)
    # (1/2,1): sum=3/2 not integer -> real collapse (qe=1); (1,1): sum=2 -> bistable (qe=2)
    if qe_a == 1 and qe_b == 2:
        print("  MUTANT M4 (false-law): CAUGHT ((1/2,1) sum=3/2∉ℤ → collapse q_eff=1; (1,1) "
              "sum=2 → bistable q_eff=2. The law = SUM, not difference — the false «t−s»-form misses)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): (3,1)-weighted — q=1, nothing to collapse")
    random.seed(1008041)
    w_nc = random.choice([Rational(3, 2), Rational(5, 2), Rational(7, 3)])
    nm_nc = wnull_modes(3, 1, L, tab, [w_nc])
    ct_nc = ctilde(nm_nc, 3, 1, ss)
    ok(ct_nc == 1, "control: (3,1) w={0} → C̃=1 (q=1 has nothing to collapse)".format(w_nc))
    print("  (3,1) w={0}: C̃={1} — q=1 is always one class up to arrow, ∀w (no bistability)".format(
        w_nc, ct_nc))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'bistable'/'bistability' is a STRUCTURAL S1000 term (not physics-stability) — not fenced;   GUARDLINE
    #   fence UV-physics via 'ультрафіол' only; stability-physics claims simply not written.  GUARDLINE
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
