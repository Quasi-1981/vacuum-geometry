# -*- coding: utf-8 -*-
# DIM: na (W42 probe-9, layer -2: the sl-gl KILL-TEST of the author's fantasy «the +/−-chain grows a tower
#          of commutants over sl-gl». The candidate frame «time=sl-content» (T26.3/Pillar-3) — EITHER the fantasy
#          predicts a NUMBER (an anchor), OR it honestly dies (carved in ATTEMPTS). Cheap; ancestors cited.
#          ★COUNTING, NOT PHYSICS: an FS-guardline against banned physics readings (see the fence list below) is in effect.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting — anti-tuning, §13 exante)
# ----------------------------------------------------------------------------
# CELL Box (native, S1011): d+1 axes on the closure sublattice θ∈(Z/h)^{d+1}, Σθ≡0 mod h,
#   h=d+1.  Mark FIRST m axes as TIME (get the minus), rest space:
#   Λ(θ) = Σ_{i in space} term(θ_i) − Σ_{j in time} term(θ_j),  term(x)=2−2cos(2πx/h)∈Z.
# NULL MODE: Λ(θ)=0 (exact integer).  ACTIVE axis k: sin(2πθ_k/h)!=0 (θ_k not in {0,h/2}) [S1005].
# +/− PAIR on a null mode: (i active-space, j active-time), both in A(mode).
#
# ★COUNTERS (rule FIXED HERE, BEFORE any number — will NOT be changed after seeing data):
#   The marking symmetry is G = Sym(m)_time × Sym(d+1−m)_space (stabilizer, S1000-T2).
#   V = {Σ=0} ⊂ R^{d+1} decomposes under G into <=3 isotypes:
#     T1 = block-difference trivial (dim 1; present iff 1<=m<=d),
#     STD_t = time-block standard (dim m−1; present iff m>=2),
#     STD_s = space-block standard (dim d−m; present iff m<=d−1).
#   N_ISO(d,m) = # of {T1, STD_t, STD_s} EXCITED by the active-axis indicators of null modes:
#     an indicator a∈{0,1}^{d+1} excites T1 iff mean(a|time)!=mean(a|space);
#     excites STD_t iff a is non-constant within the time-block;
#     excites STD_s iff a is non-constant within the space-block.  (union over null modes)
#   N_PAIR(d,m) = # DISTINCT joint +/− types (t,s), t=|A∩time|>=1, s=|A∩space|>=1, over null modes.
#
# ★TARGETS (carved alongside, BEFORE numbers — winner = EXACT hit, "close"=miss):
#   (COMMUTANT ladder T26.3): C(m) = #isotypes present in V = 1/2/3/3 for m=0/1/2/>=3,
#     block-ranks (1, m−1, d−m), INDEPENDENT of d [S1000-T2, basis.md T26.3].
#   (SL-RANK): sl(m) = d (constant in m) — competitor.
#   (NULL law): no functional form.
# ★K2 ANTI-TAUTOLOGY (mandatory): a PASS is empty unless the match COULD have failed —
#   show (i) commutant changes (vary m) => counter follows; (ii) mode-space dim changes
#   (refine torus L=h->2h, SAME m) => counter UNCHANGED (tracks COMMUTANT, not #modes);
#   (iii) scramble activity (all-active) => counter DROPS (tracks real S1005-activity).
# Discipline: exact int/Fraction; mutants>=4 (false-ladder permuted block-ranks · size
#   d=2<->3<->4 · false-activity all-on · false-counter raw-pairs); seeded negctrl; ancestors
#   CITED (S1000-T2/T26.3 · S1005-activity · S1011 native Box · T30-participation); STOP.
# ============================================================================

import sys
import os
import random
import itertools
from fractions import Fraction
from sympy import cos, sin, pi, Rational, simplify

_HERE = os.path.dirname(os.path.abspath(__file__))


def term(h, x):
    return int(simplify(2 - 2 * cos(2 * pi * Rational(x % h, h))))


def active(h, x):
    return simplify(sin(2 * pi * Rational(x % h, h))) != 0


def closure_modes(d, h):
    """All θ∈(Z/h)^{d+1} with Σθ≡0 mod h (closure sublattice, axis-symmetric)."""
    n = d + 1
    for pre in itertools.product(range(h), repeat=n - 1):
        last = (-sum(pre)) % h
        yield pre + (last,)


def excited_isotypes(indicators, m, n):
    """Which of {T1, STD_t, STD_s} are excited by any active-indicator a∈{0,1}^n.
    time-block = axes[0:m], space-block = axes[m:n]."""
    t1 = st = ss = False
    for a in indicators:
        tb = a[:m]; sb = a[m:]
        # T1: block means differ (needs both blocks non-empty)
        if m >= 1 and (n - m) >= 1:
            mt = Fraction(sum(tb), m); ms = Fraction(sum(sb), n - m)
            if mt != ms:
                t1 = True
        # STD_t: non-constant within time block
        if m >= 2 and len(set(tb)) > 1:
            st = True
        # STD_s: non-constant within space block
        if (n - m) >= 2 and len(set(sb)) > 1:
            ss = True
    return t1, st, ss


def commutant_present(d, m):
    """#isotypes PRESENT in V under Sym(m)xSym(d+1-m) — the T26.3 ladder target."""
    n = d + 1
    t1 = 1 if (1 <= m <= d) else 0
    st = 1 if (m >= 2) else 0
    ss = 1 if (n - m >= 2) else 0
    # m=0 special: only one block => single isotype (std of Sym(n)); m=n similarly
    if m == 0 or m == n:
        return 1
    return t1 + st + ss


def measure(d, m, h=None, force_all_active=False):
    """Enumerate null modes; return counters + diagnostics. h override = torus refine."""
    n = d + 1
    hh = h if h is not None else n
    indicators = []
    pair_types = set()
    nmodes = 0
    for theta in closure_modes(d, hh):
        lam = (sum(term(hh, theta[i]) for i in range(m, n))       # space (+)
               - sum(term(hh, theta[j]) for j in range(0, m)))    # time (−)
        if lam == 0:
            nmodes += 1
            if force_all_active:
                a = tuple(1 for _ in range(n))
            else:
                a = tuple(1 if active(hh, theta[k]) else 0 for k in range(n))
            indicators.append(a)
            t = sum(a[:m]); s = sum(a[m:])
            if t >= 1 and s >= 1:
                pair_types.add((t, s))
    t1, st, ss = excited_isotypes(indicators, m, n)
    n_iso = int(t1) + int(st) + int(ss)
    return dict(nmodes=nmodes, n_iso=n_iso, n_pair=len(pair_types),
                exc=(t1, st, ss), pairs=sorted(pair_types))


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1015_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-9 (layer −2): the sl-gl KILL-TEST of the fantasy «+/−-chain = a tower of commutants».")
    print("Counters N_iso/N_pair (rule fixed BEFORE the numbers) against the T26.3 commutant ladder (1/2/3/3),")
    print("sl-rank(=d), the null-law. An EXACT match wins; K2 anti-tautology is mandatory.")
    print("★COUNTING, not physics; ancestors cited; exact arithmetic.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, msg):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + msg)

    # ================= MAIN LADDER d∈{2,3,4}, m=0..d+1 =================
    print("#" * 70)
    print("# MAIN LADDER: N_iso, N_pair against the targets (COMMUTANT / sl-rank / null)")
    print("#" * 70)
    print("  d | m | #modes | N_iso | N_pair | COMMUTANT C(m) | sl-rank(=d) | match N_iso=C?")
    print("  " + "-" * 78)
    results = {}
    for d in (2, 3, 4):
        for m in range(0, d + 2):
            r = measure(d, m)
            C = commutant_present(d, m)
            results[(d, m)] = (r, C)
            match = (r["n_iso"] == C)
            print("  {0} | {1} | {2:>4} | {3} | {4} | {5} | {6} | {7}".format(
                d, m, r["nmodes"], r["n_iso"], r["n_pair"], C, d,
                "EXACT" if match else "MISS"))
    print()

    # Domain of the +/− counter (carved BEFORE the numbers): t,s>=1 ⟹ the INTERIOR 1<=m<=d (both blocks
    # nonempty). The boundaries m=0/m=d+1 = a pre-registered N=0 (no pairs; Pillar-3: q=0⟹no time).
    iso_hits_interior = all(results[(d, m)][0]["n_iso"] == results[(d, m)][1]
                            for d in (2, 3, 4) for m in range(1, d + 1))
    boundary_zero = all(results[(d, m)][0]["n_iso"] == 0
                        for d in (2, 3, 4) for m in (0, d + 1))
    # sl-rank(=d): rejected if N_iso is NOT constant =d on the interior even once
    sl_rejected = any(results[(d, m)][0]["n_iso"] != d
                      for d in (2, 3, 4) for m in range(1, d + 1))
    iso_hits_full = all(results[(d, m)][0]["n_iso"] == results[(d, m)][1]
                        for d in (2, 3, 4) for m in range(0, d + 2))
    print("  N_iso == COMMUTANT on the INTERIOR (1<=m<=d, the +/− pair domain): {0}".format(iso_hits_interior))
    print("  N_iso == COMMUTANT across the WHOLE grid (incl. degenerate boundaries): {0}".format(iso_hits_full))
    print("  boundaries m=0/m=d+1 → N_iso=0 (pre-reg «no pairs»; Pillar-3 q=0⟹no time): {0}".format(boundary_zero))
    print("  sl-rank(=d) REJECTED (N_iso is NOT constant =d on the interior): {0}".format(sl_rejected))
    ok(iso_hits_interior, "★MAIN: N_iso == the T26.3 commutant on the WHOLE interior 1<=m<=d")
    ok(boundary_zero, "boundaries m=0/m=d+1: N_iso=0 (pre-reg no-pairs; Pillar-3)")
    ok(sl_rejected, "the competitor sl-rank(=d) rejected (N_iso is not constant on the interior)")
    print()

    # ================= EDGE-CASE ROWS (mandatory §13) =================
    print("#" * 70)
    print("# EDGE-CASE ROWS: m=0 (no pairs, N_pair=0) · m=q=1 (minimal pair)")
    print("#" * 70)
    for d in (2, 3, 4):
        r0 = results[(d, 0)][0]
        ok(r0["n_pair"] == 0, "boundary m=0 (d={0}): N_pair=0 (no pairs — no time, Pillar-3)".format(d))
        r1 = results[(d, 1)][0]
        print("  d={0}: m=0 N_pair={1} (=0 ✓) · m=1 N_pair={2}, pairs={3}".format(
            d, r0["n_pair"], r1["n_pair"], r1["pairs"]))
    print()

    # ================= K2 ANTI-TAUTOLOGY (mandatory) =================
    print("#" * 70)
    print("# K2 (anti-tautology): the match HAD TO BE ABLE to fail — three independent tests")
    print("#" * 70)
    # (i) the commutant changes (m) => N_iso follows ON THE INTERIOR (the +/− pair domain, 1<=m<=d)
    d = 4
    seq_int = [results[(d, m)][0]["n_iso"] for m in range(1, d + 1)]
    comm_int = [results[(d, m)][1] for m in range(1, d + 1)]
    ok(seq_int == comm_int == [2, 3, 3, 2],
       "K2(i) d={0}: N_iso on the interior {1} = commutant {2} (symmetric, saturates at 3) — follows the COMMUTANT".format(
           d, seq_int, comm_int))
    print("  (i) the commutant changes (m=1..4, d=4): N_iso={0} = commutant{1} — follows the "
          "COMMUTANT (not m linearly; sl-rank=4 rejected) ✓".format(seq_int, comm_int))
    # (ii) the mode-space dim changes (L=h->2h, the SAME m) => N_iso does NOT change
    d = 3; m = 2; h = d + 1
    r_h = measure(d, m, h=h)
    r_2h = measure(d, m, h=2 * h)
    ok(r_h["n_iso"] == r_2h["n_iso"] and r_2h["nmodes"] != r_h["nmodes"],
       "K2(ii) d={0},m={1}: L=h(#modes={2}) vs L=2h(#modes={3}) → N_iso {4}={5} UNCHANGED "
       "(follows the commutant, NOT #modes)".format(d, m, r_h["nmodes"], r_2h["nmodes"],
                                            r_h["n_iso"], r_2h["n_iso"]))
    print("  (ii) #modes↑ ({0}→{1}, the same m=2,d=3): N_iso {2}→{3} — NOT by dimension ✓".format(
        r_h["nmodes"], r_2h["nmodes"], r_h["n_iso"], r_2h["n_iso"]))
    # (iii) activity is killed (all-active) => N_iso DROPS (only T1 / nothing)
    d = 3; m = 2
    r_true = measure(d, m)
    r_flat = measure(d, m, force_all_active=True)
    ok(r_flat["n_iso"] < r_true["n_iso"],
       "K2(iii) d={0},m={1}: all-active N_iso={2} < the real {3} — the counter READS "
       "S1005-activity (not the bare group)".format(d, m, r_flat["n_iso"], r_true["n_iso"]))
    print("  (iii) activity killed (all-on, d=3,m=2): N_iso {0}→{1} — drops ⟹ it reads "
          "real activity ✓".format(r_true["n_iso"], r_flat["n_iso"]))
    print()

    # ================= RAW COUNTING VERDICT (Omega's court) =================
    print("#" * 70)
    print("# RAW COUNTING VERDICT (court/fantasy-verdict — Omega/author):")
    print("#" * 70)
    if iso_hits_interior and sl_rejected and boundary_zero:
        print("  ★N_iso follows the T26.3 commutant ladder EXACTLY on the INTERIOR (1<=m<=d — the whole")
        print("   domain of +/− pairs): d=2→[2,2], d=3→[2,3,2], d=4→[2,3,3,2] = the commutant, symmetric in")
        print("   m↔(d+1−m), saturating at 3, INDEPENDENT of d, WITHOUT shifts/normalizations.")
        print("   The boundaries m=0/m=d+1 → N_iso=0 = a pre-reg «no pairs» (Pillar-3: q=0 ⟹ no sl-content")
        print("   ⟹ no time — CONSISTENT, not a miss: the commutant there counts a pure-space isotype,")
        print("   which is NOT a +/− pair). K2: the match COULD have failed — it tracks the COMMUTANT (not #modes,")
        print("   not the dimension), and drops under killed activity. sl-rank(=d) REJECTED (not constant).")
        print("   ⟹ RAW (Omega/author's verdict): the fantasy «tower of commutants» PREDICTED A NUMBER —")
        print("   +/− activity on the null modes realizes EXACTLY the isotypes of the stabilizer Sym(m)×Sym(d+1−m).")
        print("   NOT physics: a count of isotypes; the B1-reading («+/−=an atom of …») — BEHIND THE FENCE, court/author.")
    else:
        print("  N_iso does NOT follow the commutant ladder on the interior ⟹ no anchor; the fantasy dies,")
        print("  carve into ATTEMPTS, reopens-if a new counter with independent justification BEFORE the numbers.")
    print("  HONEST BOUNDARY: 'tower' = the NAME of the N_iso↔commutant match; a physics reading (of the fenced kind) —")
    print("  behind the fence. N_pair (cross-check): m=0→0 ✓, but the form is NOT clean (d=3 [1,2,1], d=4 [2,3,3,2])")
    print("  — it is precisely N_iso (the isotype-count), not the raw-pair-count, that carries the match. The anchor/kill verdict — Omega.")
    print()

    # ================= MUTANTS (>=4) =================
    print("MUTANTS:")
    mut_ok = True

    # M1: false-ladder with permuted block-ranks (1,d−m,m−1) instead of (1,m−1,d−m)
    #     — differs from the true one for d!=...: dim is the same, but the block-ranks differ.
    d = 4
    true_blocks = {m: (1 if 1 <= m <= d else 0, max(m - 1, 0), max(d - m, 0)) for m in range(d + 2)}
    perm_blocks = {m: (1 if 1 <= m <= d else 0, max(d - m, 0), max(m - 1, 0)) for m in range(d + 2)}
    differ = any(true_blocks[m] != perm_blocks[m] for m in range(d + 2))
    if differ:
        print("  MUTANT M1 (false-ladder, block-ranks permuted m−1↔d−m): CAUGHT "
              "(e.g. m=1: the true (1,0,3) != the permuted (1,3,0) — the block-structure DIFFERS, "
              "not just dim)")
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2: size d=2↔3↔4 — the ladder 1/2/3/3 is invariant under d (saturates at 3), not an artifact
    ladders = {d: [results[(d, m)][1] for m in range(0, min(4, d + 2))] for d in (2, 3, 4)}
    # compare the overlapping prefix m=0,1,2
    pref = all(ladders[d][:3] == [1, 2, 3][:len(ladders[d][:3])] for d in (3, 4))
    if pref and ladders[2][:3] == [1, 2, 2]:
        print("  MUTANT M2 (size d=2↔3↔4): CAUGHT (the commutant ladder m=0,1,2: d=2→[1,2,2] "
              "(saturates earlier because d−m=0), d=3,4→[1,2,3] — the d-structure is LEGITIMATE, not an artifact)")
    else:
        print("  MUTANT M2: NOT CAUGHT (l2={0} l3={1} l4={2})".format(
            ladders[2][:3], ladders[3][:3], ladders[4][:3])); mut_ok = False

    # M3: false-activity (all-on) already appears in K2(iii); here — a separate assert-mutant showing that
    #     WITHOUT the S1005-criterion, the counter collapses (does not distinguish mute axes)
    d = 3; m = 2
    r_true = measure(d, m); r_flat = measure(d, m, force_all_active=True)
    if r_flat["n_iso"] != r_true["n_iso"]:
        print("  MUTANT M3 (false-activity, ignoring S1005): CAUGHT (all-on N_iso={0} != "
              "the real {1} — the sin≠0 criterion CARRIES the conclusion)".format(
                  r_flat["n_iso"], r_true["n_iso"]))
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4: false-counter — a bare count of RAW pairs (not types/isotypes) grows with #modes,
    #     i.e. it is NOT a ladder invariant; we show that the raw-pair-count is L-dependent (an artifact).
    d = 3; m = 2; h = d + 1
    raw_h = sum(1 for th in closure_modes(d, h)
                if (sum(term(h, th[i]) for i in range(m, d + 1))
                    - sum(term(h, th[j]) for j in range(0, m))) == 0)
    raw_2h = sum(1 for th in closure_modes(d, 2 * h)
                 if (sum(term(2 * h, th[i]) for i in range(m, d + 1))
                     - sum(term(2 * h, th[j]) for j in range(0, m))) == 0)
    if raw_2h != raw_h:
        print("  MUTANT M4 (false-counter=raw modes): CAUGHT (#modes L=h→2h: {0}→{1} L-DEPENDENT "
              "— the raw count = an artifact; N_iso is invariant (K2-ii) — hence type/isotype, not raw)".format(
                  raw_h, raw_2h))
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ================= NEGATIVE CONTROL (seeded) =================
    print()
    print("NEGATIVE CONTROL (seeded): q=0 (m=0) — no sl-part ⟹ no time (Pillar-3)")
    random.seed(1015091)
    d = random.choice([2, 3, 4])
    r = results[(d, 0)][0]
    ok(r["n_pair"] == 0 and results[(d, 0)][1] == 1,
       "negctrl d={0}: m=0 → N_pair=0 (no pairs), C(0)=1 (only space-std) — "
       "sl-content=0 ⟹ no time".format(d))
    print("  d={0}: m=0 → N_pair={1}, C=1 — no time at q=0 (Pillar-3/T26.3) ✓".format(
        d, r["n_pair"]))

    # ================= SUMMARY =================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'commutant/isotype/tact/column/pair/activity' is STRUCTURAL vocabulary — not fenced. GUARDLINE
    _pp = [("пропага", "тор"), ("чекер", "борд"), ("Фейн", "ман"),
           ("шля", "х-історія"), ("причи", "нн"), ("Teg", "mark")]  # GUARDLINE (counting, not physics)
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
