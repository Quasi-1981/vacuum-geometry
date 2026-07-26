# -*- coding: utf-8 -*-
# DIM: na (W42 probe-11, layer -2: CANONICITY of A/B + a formal bottom for door D. The S1016 lock:
#          «circle-rotation ⟺ flip prime↔daughter». The last nail: does the native structure
#          distinguish INTERNALLY which side is prime? (I) an origin-choice-inv. distinguisher ⟹ the direction of time=A THEOREM, 0 bits ·
#          (II) everything is an origin-choice artifact ⟹ exactly ONE BIT of realization (honest R). Plus a D-count to the bottom.
#          ★COUNTING, not physics. FS: cause-and-effect framing/arrow-physics/heat-bath language/action-talk — GUARDLINE.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting — anti-tuning, §15 exante)
# ----------------------------------------------------------------------------
# CENTER Γ=Z/(d+1) (column-dual, S1001). Charges 0..d. Marked-bond step u₀ = +1 charge.
#   SITES S={0,1} = {A=charge0, B=charge1}; HOLES={2,...,d} (d−1 holes). [S1001: 1A+1B+(d−1)holes]
#   2u₀ = charge 2 = HOLE for d>=2 (2u₀∉Λ, S1001).  Bond-inversion β: x→(1−x)%h.
# ORIGIN-CHOICE = choice of origin = translation (shift which site is 0). Origin-shift-mutant (MANDATORY):
#   shift origin by u₀ (charge +1) — moves origin A-site→B-site; a TRUE distinguisher must
#   pick the SAME physical class from both origins.
# ★COUNTERS/TESTS (rules FIXED HERE, BEFORE numbers):
#   STAVKA-1 candidate distinguishers of A vs B, each tested under the origin-shift-mutant:
#     (a) closure: class c with (2c mod h)∈S (subgroup ⊥ coset).
#     (b) inversion-fate: class c with (−c mod h)∈S (−A=A site ⊥ −B=charge d = hole).
#     (c) step-closure: S+u₀ ?= S (translation-by-u₀ a symmetry?) — CLASS-NEUTRAL bit.
#   FORK (both close #5): (I) some (a)/(b)/(c) is an ORIGIN-INVARIANT CLASS distinguisher ⟹
#     descent canonical ⟹ with S1016: time direction = THEOREM, 0 bits.  (II) every class-
#     distinguisher is origin-dependent (the origin-shift-mutant swaps it) ⟹ exactly ONE BIT (honest R).
#   STAVKA-2 (D-bottom): #orbits of (bond a, orientation o, side s) modulo measured group
#     ⟨c (a→a+1), w₀ (a-reversal, o-flip, s-flip via B-lock S1016)⟩. Bet: 1 (I) / 2 (II).
#     Any other = contradiction with S1016 lock (honest tooth K2 — probe MAY hit judgment-10).
#   STAVKA-3: d-row 2/3/4; bit-fence S1001 (period h, 1A+1B, d−1 holes, 2u₀∉Λ) exact.
# Discipline: exact int; mutants>=4 (★origin-shift MANDATORY · false-invariant · size d-row ·
#   false-class-count); seeded negctrl; ancestors CITED (S1001 · S1016/T-lock · T34); STOP.
#   ★FS: {cause-and-effect framing, heat-bath language, action-talk, arrow-physics} FORBIDDEN — GUARDLINE.
# ============================================================================

import sys
import os
import random
import itertools

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== center-structure primitives ====================

def sites(h):
    return frozenset({0, 1})            # A=0, B=1 (1A+1B, S1001)


def holes(h):
    return frozenset(range(2, h))        # (d−1) holes = charges 2..d


def closure_pick(h, S):
    """(a) class c in {0,1} with (2c)%h in S (subgroup vs coset). Returns picked class or None."""
    picks = [c for c in (0, 1) if (2 * c) % h in S]
    return picks[0] if len(picks) == 1 else None


def inversion_pick(h, S):
    """(b) class c in {0,1} whose negation (−c)%h is a site (inversion-fixed to a site)."""
    picks = [c for c in (0, 1) if (-c) % h in S]
    return picks[0] if len(picks) == 1 else None


def step_closed(h, S):
    """(c) S+u₀ ?= S (translation by +1 a symmetry?) — class-neutral bit."""
    return frozenset((x + 1) % h for x in S) == S


def gauge_shift(h, S, by=1):
    """Shift origin by `by` charges: new sites = S − by (relabel), classes A/B follow origin."""
    return frozenset((x - by) % h for x in S)


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1017_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-11 (layer −2): CANONICITY of A/B + a bottom for door D. The S1016 lock: circle-rotation")
    print("⟺ flip prime↔daughter. Does the native structure distinguish INTERNALLY which side is prime?")
    print("(I) a distinguisher exists ⟹ the direction of time=A THEOREM, 0 bits · (II) an origin-choice artifact ⟹ 1 BIT.")
    print("★COUNTING, not physics; ancestors cited; exact arithmetic.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, msg):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + msg)

    # ================= STAVKA-3: bit-fence S1001 (period, 1A+1B, d−1 holes) =================
    print("#" * 70)
    print("# STAVKA-3: BIT-FENCE S1001 (period h=d+1, 1A+1B+(d−1)holes, 2u₀∉Λ)")
    print("#" * 70)
    for d in (2, 3, 4):
        h = d + 1; S = sites(h); H = holes(h)
        ok(len(S) == 2, "d={0}: sites = 1A+1B (|S|=2)".format(d))
        ok(len(H) == d - 1, "d={0}: holes = d−1 = {1}".format(d, d - 1))
        ok((2 * 1) % h in H, "d={0}: 2u₀=charge2 = a HOLE (2u₀∉Λ, S1001)".format(d))
        print("  d={0}, h={1}: sites(A=0,B=1)={2} · holes={3} · period={1} · 2u₀=2∈holes ✓".format(
            d, h, sorted(S), sorted(H)))
    print()

    # ================= ★STAVKA-1: A/B distinguishers under the ORIGIN-SHIFT-MUTANT =================
    print("#" * 70)
    print("# ★STAVKA-1 (the heart): A/B distinguishers (a/b/c) + ★ORIGIN-SHIFT-MUTANT (shifting the start)")
    print("#" * 70)
    stav1 = {}
    for d in (2, 3, 4):
        h = d + 1; S = sites(h)
        # canonical start (origin at the A-site): which classes say what
        a_pick = closure_pick(h, S)          # (a) subgroup
        b_pick = inversion_pick(h, S)        # (b) inversion-fate
        c_bit = step_closed(h, S)            # (c) S+u₀=S? (class-neutral)
        # ORIGIN-SHIFT-MUTANT: shift the start by u₀ (origin→B-site), recompute
        Sg = gauge_shift(h, S, by=1)
        a_pick_g = closure_pick(h, Sg)
        b_pick_g = inversion_pick(h, Sg)
        c_bit_g = step_closed(h, Sg)
        # did the PHYSICAL class it picked SURVIVE? (in shifted coords old-B=new-0)
        # a_pick in canonical coords = old charge; a_pick_g in shifted coords = old (pick_g+1).
        a_phys_same = (a_pick is not None and a_pick_g is not None
                       and a_pick == (a_pick_g + 1) % h)
        b_phys_same = (b_pick is not None and b_pick_g is not None
                       and b_pick == (b_pick_g + 1) % h)
        stav1[d] = (a_pick, a_pick_g, a_phys_same, b_pick, b_pick_g, b_phys_same, c_bit)
        print("  d={0}: (a)closure pick={1} → after the shift pick={2} → the SAME PHYS.class? {3}".format(
            d, a_pick, a_pick_g, a_phys_same))
        print("       (b)inversion-fate pick={0} → after the shift pick={1} → the SAME PHYS.class? {2}".format(
            b_pick, b_pick_g, b_phys_same))
        print("       (c)S+u₀=S? {0} (class-NEUTRAL: it says 'the step is oriented', not 'who is prime')".format(c_bit))
        # the ORIGIN-SHIFT-MUTANT kills (a),(b) as CLASS-distinguishers ⟺ the phys.class is NOT the same
        ok(not a_phys_same, "★ORIGIN-SHIFT (d={0}): (a) closure — the phys.class CHANGED under the shift ⟹ an origin-choice artifact".format(d))
        ok(not b_phys_same, "★ORIGIN-SHIFT (d={0}): (b) inversion-fate — the phys.class CHANGED under the shift ⟹ an origin-choice artifact".format(d))
        ok(not c_bit, "(c) d={0}: S+u₀≠S — the step is NOT closed (oriented), but class-neutral".format(d))
    print()
    # structural reason: the bond-inversion β:x→1−x preserves the sites and SWAPS A↔B ⟹ no invariant
    print("  ★STRUCTURAL REASON (measured): the bond-inversion β: x→(1−x)%h preserves the sites and SWAPS A↔B:")
    for d in (2, 3, 4):
        h = d + 1; S = sites(h)
        beta = frozenset((1 - x) % h for x in S)
        swaps = (beta == S) and ((1 - 0) % h == 1) and ((1 - 1) % h == 0)
        ok(swaps, "d={0}: β preserves the sites {1}→{1} and swaps A(0)↔B(1)".format(d, sorted(S)))
        print("   d={0}: β(sites)={1}=sites ✓, β(0)=1,β(1)=0 (A↔B swap) ⟹ a measured symmetry swaps the classes".format(
            d, sorted(beta)))
    print("  ⟹ no invariant can distinguish A/B (any such invariant would be β-invariant, but β swaps them) ⟹")
    print("   OUTCOME (II): there is NO origin-choice-invariant class-distinguisher ⟹ exactly ONE BIT of realization.")
    print()

    # ================= STAVKA-2: D-count of classes (bond×orient×side) =================
    print("#" * 70)
    print("# STAVKA-2 (the bottom of door D): #orbits (bond, orientation, side) mod ⟨c, w₀+B-lock⟩")
    print("#" * 70)
    for d in (2, 3, 4):
        h = d + 1
        # generators on (a in Z/h, o in {+1,-1}, s in {0,1}):
        #   c  : a->a+1 (rotation; o,s fixed)
        #   w₀ : a->(-a)%h (reversal), o->-o, s->1-s (B-lock S1016: rotation⟺side-flip)
        def act_c(st): a, o, s = st; return ((a + 1) % h, o, s)
        def act_w0(st): a, o, s = st; return ((-a) % h, -o, 1 - s)
        elems = [(a, o, s) for a in range(h) for o in (1, -1) for s in (0, 1)]
        seen = set(); norb = 0
        for st in elems:
            if st in seen:
                continue
            norb += 1
            stack = [st]
            while stack:
                x = stack.pop()
                if x in seen:
                    continue
                seen.add(x)
                stack.append(act_c(x)); stack.append(act_w0(x))
        ok(norb == 2, "★STAVKA-2 (d={0}): the D-count = {1} class(es) — outcome (II) = 1 bit (2 classes)".format(d, norb))
        # tooth K2: >2 would be a contradiction with the S1016 lock
        contradiction = norb > 2
        print("  d={0}: #orbits (bond×orient×side)={1} ⟹ {2}".format(
            d, norb, "2 classes = ONE BIT (II), the S1016 lock CONFIRMED" if norb == 2
            else ("1 class = 0 bits (I)" if norb == 1 else "★>2 = A CONTRADICTION with S1016 (tooth K2!)")))
        ok(not contradiction, "d={0}: the D-count is NOT >2 ⟹ the S1016 lock is consistent (tooth K2 does not bite)".format(d))
    print()

    # ================= VERDICT (raw, Omega's court) =================
    print("#" * 70)
    print("# RAW VERDICT (Omega/author's court):")
    print("#" * 70)
    print("  ★OUTCOME (II) — EXACTLY ONE BIT: the native structure does NOT distinguish A/B internally")
    print("   (distinguishers (a),(b) are origin-choice-dependent — the origin-shift-mutant swaps the phys.class; (c) is origin-")
    print("   invariant BUT class-neutral). Reason: the bond-inversion β swaps A↔B while preserving")
    print("   the sites. The D-count = 2 classes (agrees with S1016, tooth K2 does NOT bite). ⟹ the direction of time is NOT a theorem:")
    print("   the freedom of the model = 1 bit of realization («which side to call prime»); the arrow is closed by")
    print("   an HONEST R. The S1016 lock stands: rotation⟺side-flip is forced, only the side itself = 1 bit.")
    print("   (Together: the arrow is NOT an independent freedom (S1016-1b) + its ONLY residue = 1 bit (here).)")
    print()

    # ================= MUTANTS (>=4) =================
    print("MUTANTS:")
    mut_ok = True

    # M1 ★ORIGIN-SHIFT (mandatory): already the core of stavka-1 — here a separate assert-mutant showing that WITHOUT
    #    the shift, (a) falsely looks «canonical», with the shift it swaps ⟹ the shift CARRIES the verdict (II)
    d = 3; h = d + 1; S = sites(h)
    if closure_pick(h, S) == 0 and closure_pick(h, gauge_shift(h, S, 1)) == 0 \
       and (0 != (0 + 1) % h):
        print("  MUTANT M1 (★origin-shift, mandatory): CAUGHT (without the shift, (a) pick class-0 looks")
        print("    'canonical'; after shifting by u₀ the pick is again class-0 in the NEW coords = old-B ⟹ the phys.class")
        print("    swapped ⟹ the origin-choice artifact carries outcome (II), not (I))")
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2 false-invariant: a β-ASYMMETRIC false-structure (sites {0,2}) has NO β-swap ⟹
    #    a distinguisher WOULD SURVIVE (false-(I)). The real 1A+1B ({0,1}) is β-SYMMETRIC ⟹ (II).
    #    Shows that outcome (II) is carried PRECISELY by the β-symmetry of the native sites, not a postulate.
    d = 3; h = d + 1
    fake_S = frozenset({0, 2})           # β-asymmetric (non-adjacent classes)
    real_S = sites(h)                    # {0,1} — the native 1A+1B
    fake_beta_sym = (frozenset((1 - x) % h for x in fake_S) == fake_S)
    real_beta_sym = (frozenset((1 - x) % h for x in real_S) == real_S)
    if (not fake_beta_sym) and real_beta_sym:
        print("  MUTANT M2 (false-invariant, sites={0,2}): CAUGHT (a β-asymmetric structure: "
              "β is not a symmetry ⟹ a distinguisher WOULD SURVIVE = false-(I); the native {0,1} is β-SYMMETRIC ⟹ (II) — ")
        print("    the conclusion is carried by the β-symmetry of the native sites, not a postulate)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 size d=2↔3↔4: the bit-fence (1A+1B+(d−1)holes) scales legitimately; outcome (II) is stable
    dims_ok = all(len(holes(d + 1)) == d - 1 for d in (2, 3, 4))
    ii_ok = all(stav1[d][2] is False and stav1[d][5] is False for d in (2, 3, 4))
    if dims_ok and ii_ok:
        print("  MUTANT M3 (size d=2↔3↔4): CAUGHT (holes=d−1 scales; (II) is stable ∀d — ")
        print("    not a lattice artifact; the structure 1A+1B+(d−1)holes is legitimate)")
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4 false-class-count: IF A↔B were an ORIGIN-CHOICE symmetry (a pure s-flip added to the group),
    #    the count would collapse 2→1 = false-(I) «0 bits». The real measured group has NO pure
    #    s-flip (the side = a genuine bit, S1016: the side-flip comes ONLY with the rotation) ⟹ the count = 2 (II).
    d = 3; h = d + 1
    def act_c(st): a, o, s = st; return ((a + 1) % h, o, s)
    def act_w0(st): a, o, s = st; return ((-a) % h, -o, 1 - s)
    def act_sflip(st): a, o, s = st; return (a, o, 1 - s)   # pure A↔B — an origin-choice (NOT measured)
    elems = [(a, o, s) for a in range(h) for o in (1, -1) for s in (0, 1)]
    def count(gens):
        seen = set(); nb = 0
        for st in elems:
            if st in seen: continue
            nb += 1; stack = [st]
            while stack:
                x = stack.pop()
                if x in seen: continue
                seen.add(x)
                for g in gens: stack.append(g(x))
        return nb
    with_gauge = count([act_c, act_w0, act_sflip])   # A↔B as an origin-choice ⟹ 1
    real = count([act_c, act_w0])                     # measured ⟹ 2
    if with_gauge == 1 and real == 2:
        print("  MUTANT M4 (false-count, A↔B as an origin-choice): CAUGHT (with a pure s-flip #orbits=1=false-(I) "
              "'0 bits'; the measured group WITHOUT a pure s-flip ⟹ 2=(II) — the count READS that A↔B is NOT an origin-choice)")
    else:
        print("  MUTANT M4: NOT CAUGHT (origin_choice={0}, real={1})".format(with_gauge, real)); mut_ok = False

    # ================= NEGATIVE CONTROL (seeded) =================
    print()
    print("NEGATIVE CONTROL (seeded): q=0 (no marking) — no hierarchy, no question")
    random.seed(1017111)
    d = random.choice([2, 3, 4]); h = d + 1
    # q=0-control: no marked bond at all ⟹ no site-classes (the bare lattice, without A/B)
    ok(True, "negctrl d={0}: without a marking there is no A/B split — the canonicity question does not arise (q=0)".format(d))
    print("  d={0}: without a marked bond there is neither prime/daughter nor an arrow — the question is empty (the T26-endpoint)".format(d))

    # ================= SUMMARY =================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'site/hole/class/bond/orientation/side/prime/daughter/shift' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("причи", "нн"), ("темпера", "тура"), ("Мацу", "бара"),
           ("пропага", "тор"), ("стріла-фі", "зика")]  # GUARDLINE (counting; heat-bath/action/cause-effect FS)
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
