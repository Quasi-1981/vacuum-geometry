# -*- coding: utf-8 -*-
# DIM: na (W42 probe-1, layer -2 ultrahyperbolic: a discrete ultrahyperbolic Box on small
#          exact lattices; a kill-test of seed ALTER_UHP; exact integer arithmetic; 0 handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — derivational, seed+exante known
# ----------------------------------------------------------------------------
# DISCRETE BOX (signature (p,q)) on a periodic lattice L^n, n=p+q, diagonal in the
#   Fourier basis: eigenvalue at momentum k=(k_1..k_n), k_i in {0..L-1}:
#     lam(k) = sum_{i=1..p} T(k_i) - sum_{j=1..q} T(k_{p+j}),  T(k)=2-2cos(2 pi k/L).
#   For L in {3,4} T(k) is an EXACT INTEGER (L=3: {0,3}; L=4: {0,2,4}) => the whole
#   spectrum is exact integers; no matrix diagonalization needed.
# T1 (spectrum near zero): counts n0/npos/nneg, near-zero density on BOTH signs,
#   and the sign-symmetry lam<->-lam (swap a + and a - axis) — present iff p=q.
#   Claim tested: (2,2)/(3,3) dense-both-sides (symmetric) vs (3,1) one-sided/exclusion-zone.
# T2 (ln det branch around zero): a null mode k (lam=0) can be encircled through a
#   time-axis j iff that axis is "active" (sin(2 pi k_j/L) != 0, i.e. k_j not in
#   {0, L/2}).  BRANCH is unique iff every null mode has <=1 active time-axis (one
#   iε-direction); it SPLITS iff some null mode has >=2 active time-axes (>=2
#   inequivalent iε-directions to go around the zero).  q=1 => 0 split-modes
#   (single time) ; q>=2 => count > 0.  Exact combinatorial measure of the seed's
#   "q=1 one traversal ⊥ q>=2 splitting".
# T3 (measure leg, separate): sqrt(-g) phase; det g = (-1)^q, sqrt(-(-1)^q) is REAL
#   for q ODD, IMAGINARY for q EVEN.  => (3,1) real, (2,2) imag, (3,3) real.  Ш1:
#   (3,3) has REAL measure yet the T1/T2 pathology is ALIVE => measure != selector.
# T4: honest table across (3,1)/(2,2)/(3,3).  Size-mutant: L=3 vs L=4, conclusion
#   invariant (box+1 does not flip the verdict).
# Discipline: 0 handles; exact integers; mutants>=4 + size-mutant; seeded negctrl;
#   FORBIDDEN-SCAN {physics-interpretation}; log bit-reproducible; STOP after tables.
# ============================================================================

import sys
import os
import random
import itertools
from sympy import cos, pi, Rational, simplify, Integer

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== exact term table ====================

def term_table(L):
    """T(k) = 2 - 2cos(2 pi k/L) for k=0..L-1, exact integer for L in {3,4,6}."""
    tab = []
    for k in range(L):
        v = simplify(2 - 2 * cos(2 * pi * Rational(k, L)))
        tab.append(int(v))
    return tab


def active_set(L):
    """k with sin(2 pi k/L) != 0  <=>  k not in {0, L/2}."""
    s = set(range(L))
    s.discard(0)
    if L % 2 == 0:
        s.discard(L // 2)
    return s


# ==================== spectrum + measures ====================

def spectrum(p, q, L, tab):
    """List of exact integer eigenvalues lam(k) over all k in {0..L-1}^(p+q)."""
    n = p + q
    out = []
    for k in itertools.product(range(L), repeat=n):
        sp = sum(tab[k[i]] for i in range(p))
        sq = sum(tab[k[p + j]] for j in range(q))
        out.append(sp - sq)
    return out


def null_modes(p, q, L, tab):
    """All k with lam(k)=0."""
    n = p + q
    out = []
    for k in itertools.product(range(L), repeat=n):
        sp = sum(tab[k[i]] for i in range(p))
        sq = sum(tab[k[p + j]] for j in range(q))
        if sp - sq == 0:
            out.append(k)
    return out


def branch_split_count(p, q, L, tab):
    """# null modes with >=2 active time-axes (>=2 iε-directions to encircle)."""
    act = active_set(L)
    cnt = 0
    maxact = 0
    for k in null_modes(p, q, L, tab):
        n_active_time = sum(1 for j in range(q) if k[p + j] in act)
        maxact = max(maxact, n_active_time)
        if n_active_time >= 2:
            cnt += 1
    return cnt, maxact


def near_zero(spec, radius=2):
    npos = sum(1 for x in spec if 0 < x <= radius)
    nneg = sum(1 for x in spec if -radius <= x < 0)
    nz = sum(1 for x in spec if x == 0)
    return nz, npos, nneg


def is_sign_symmetric(spec):
    from collections import Counter
    c = Counter(spec)
    return all(c[x] == c[-x] for x in c)


def measure_real(q):
    """sqrt(-g) real  <=>  q odd."""
    return (q % 2 == 1)


# ==================== master ====================

class Tee:
    def __init__(self, real, fh):
        self.real = real
        self.fh = fh
        self.chunks = []

    def write(self, s):
        self.real.write(s)
        self.fh.write(s)
        self.fh.flush()
        self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1005_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf)
    sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-1 (layer −2, ultrahyperbolic): a discrete ultrahyperbolic Box")
    print("signatures (3,1)/(2,2)/(3,3) [Ш1-control]; exact integer eigenvalues; a kill-test of the seed.")
    print("=" * 74)
    print()

    ASSERT_PASS = [0]
    FAILS = [0]

    def ok(cond, msg):
        if cond:
            ASSERT_PASS[0] += 1
        else:
            FAILS[0] += 1
            print("ASSERT-FAIL: " + msg)

    SIGS = [(3, 1), (2, 2), (3, 3)]
    L = 4
    tab4 = term_table(L)
    ok(tab4 == [0, 2, 4, 2], "term-table L=4 = [0,2,4,2] (exact integers)")
    tab3 = term_table(3)
    ok(tab3 == [0, 3, 3], "term-table L=3 = [0,3,3] (exact integers)")

    # precompute spectra
    SPEC = {(p, q): spectrum(p, q, L, tab4) for (p, q) in SIGS}

    # ==================== T1: spectrum near zero ====================
    print("T1 (spectrum near zero, L=4): density on both sides + the λ↔−λ symmetry")
    print("signature | n0 (λ=0) | npos | nneg | near-zero ±2 (0/+/−) | λ↔−λ symmetric | one-sided?")
    print("-" * 100)
    t1 = {}
    for (p, q) in SIGS:
        spec = SPEC[(p, q)]
        n0 = spec.count(0)
        npos = sum(1 for x in spec if x > 0)
        nneg = sum(1 for x in spec if x < 0)
        nz, np2, nn2 = near_zero(spec, 2)
        sym = is_sign_symmetric(spec)
        onesided = (nneg == 0 or npos == 0)
        # symmetric <=> p==q
        ok(sym == (p == q), "T1 ({0},{1}): λ↔−λ symmetry = (p==q)".format(p, q))
        t1[(p, q)] = dict(n0=n0, npos=npos, nneg=nneg, sym=sym)
        print("({0},{1}) | {2} | {3} | {4} | 0:{5} +:{6} −:{7} | {8} | {9}".format(
            p, q, n0, npos, nneg, nz, np2, nn2, "yes" if sym else "no",
            "yes (one-sided)" if onesided else "no (two-sided)"))
    print("  ★MEASUREMENT: (2,2)/(3,3) p=q ⟹ the spectrum is EXACTLY SYMMETRIC λ↔−λ (dense on BOTH sides,")
    print("  npos=nneg). (3,1) p≠q ⟹ STRONG ASYMMETRY (leaning to the + side: npos≫nneg), BUT")
    print("  HONESTLY: nneg≠0 (not a hard exclusion-zone/one-sidedness) — the distinction = symmetry⊥asymmetry,")
    print("  not exclusion-zone⊥density. Refinement of stake-1: the qualitative difference holds, «exclusion-zone» does not.")

    # ==================== T2: ln det branch ====================
    print()
    print("T2 (the ln det branch when encircling zero, L=4): null modes with ≥2 active time-axes")
    print("signature | # null modes | max active time-axes | # split modes (≥2) | branch")
    print("-" * 92)
    t2 = {}
    for (p, q) in SIGS:
        nm = null_modes(p, q, L, tab4)
        split, maxact = branch_split_count(p, q, L, tab4)
        branch = "SINGLE (q=1: 1 iε-traversal)" if q == 1 else \
                 ("SPLIT ({0} split-modes, ≥2 iε-directions)".format(split) if split > 0
                  else "single (no ≥2-active)")
        ok((q == 1 and split == 0) or (q >= 2 and split > 0),
           "T2 ({0},{1}): q=1→split=0 / q≥2→split>0".format(p, q))
        t2[(p, q)] = dict(nm=len(nm), split=split, maxact=maxact)
        print("({0},{1}) | {2} | {3} | {4} | {5}".format(
            p, q, len(nm), maxact, split, branch))
    print("  ★MEASUREMENT: q=1 (3,1) → 0 split-modes (a single iε-traversal, the ln det branch is unambiguous);")
    print("  q≥2 (2,2)/(3,3) → >0 split-modes (≥2 independent traversal time-directions ⟹ the branch")
    print("  SPLITS). The seed's stake-2 (q=1 single ⊥ q≥2 splitting) — measured exactly.")

    # ==================== T3: measure leg (separate) ====================
    print()
    print("T3 (the measure-leg, separately): the √−g phase (real ⟺ q odd) vs the T1/T2 pathology")
    print("signature | q | √−g | pathology (T1 symmetric-both / T2 split) | measure=selector?")
    print("-" * 96)
    for (p, q) in SIGS:
        real_m = measure_real(q)
        patho = (t1[(p, q)]["sym"] and t1[(p, q)]["nneg"] > 0) or (t2[(p, q)]["split"] > 0)
        # measure selects q=1?  measure real for q odd (1 AND 3) -> does NOT isolate q=1
        print("({0},{1}) | {2} | {3} | {4} | {5}".format(
            p, q, q, "real" if real_m else "imaginary",
            "alive" if patho else "none",
            "—" if not real_m else ("NO (real measure + pathology)" if patho else "possible")))
    # Ш1 assertion: (3,3) real measure BUT pathology alive
    real33 = measure_real(3)
    patho33 = t2[(3, 3)]["split"] > 0
    ok(real33 and patho33,
       "Ш1: (3,3) the measure is REAL (q=3 odd) BUT the ultrahyperbolic pathology is ALIVE (split>0) ⟹ measure≠selector")
    print("  ★MEASUREMENT (Ш1 CONFIRMED): (3,3) — the measure √−g is REAL (q=3 odd, same as (3,1)), BUT")
    print("  the ultrahyperbolic pathology is ALIVE (split-modes exist) ⟹ the measure-leg does NOT select q=1 alone (it is real for q=3 too);")
    print("  the SELECTOR = the iε-leg (the branch), not the measure. The measure⊥iε separation — done.")

    # ==================== T4: honest table ====================
    print()
    print("T4 (an honest table over the three signatures):")
    print("signature | min(p,q) | spectrum-symmetric | ln det branch | √−g | verdict")
    print("-" * 96)
    for (p, q) in SIGS:
        mpq = min(p, q)
        sym = t1[(p, q)]["sym"]
        split = t2[(p, q)]["split"]
        real_m = measure_real(q)
        if q == 1:
            verd = "HEALTHY: a single branch + a real measure (canonical time)"
        elif mpq >= 2 and not real_m:
            verd = "ultrahyperbolic+imaginary measure (a double break)"
        else:
            verd = "ultrahyperbolic+real measure (Ш1: pathology with no measure-signal)"
        print("({0},{1}) | {2} | {3} | {4} | {5} | {6}".format(
            p, q, mpq, "yes" if sym else "no",
            "split({0})".format(split) if split > 0 else "single",
            "real" if real_m else "imaginary", verd))
    print("  ★SUMMARY: (3,1) — the ONLY one with a healthy branch (q=1); (2,2)/(3,3) ultrahyperbolic (split>0);")
    print("  the measure separates (2,2) imaginary vs (3,3) real, but the branch catches BOTH ⟹ the iε-leg =")
    print("  the true selector of q=1; the measure is only a reinforcement (imaginary) at q=2 (Ш1).")

    # ==================== MUTANTS ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    # M1 (size-mutant): conclusion (symmetry + branch-split sign) invariant L=3 vs L=4
    stable = True
    for (p, q) in SIGS:
        sym4 = is_sign_symmetric(SPEC[(p, q)])
        sym3 = is_sign_symmetric(spectrum(p, q, 3, tab3))
        sp4, _ = branch_split_count(p, q, 4, tab4)
        sp3, _ = branch_split_count(p, q, 3, tab3)
        concl4 = (sym4, sp4 > 0)
        concl3 = (sym3, sp3 > 0)
        if concl4 != concl3:
            stable = False
    if stable:
        print("  MUTANT M1 (size): CAUGHT (the hypothesis «the conclusion changes at box+1» is REFUTED —"
              " symmetry+split-sign are IDENTICAL for L=3 and L=4 across all signatures; not a lattice-artifact)")
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2: broken signature -> Euclidean (4,0), q=0: one-sided, no split, no both-sides
    spec40 = spectrum(4, 0, 4, tab4)
    onesided40 = all(x >= 0 for x in spec40)
    split40, _ = branch_split_count(4, 0, 4, tab4)
    if onesided40 and split40 == 0:
        print("  MUTANT M2: CAUGHT (broken signature (4,0) Euclidean q=0: the spectrum is one-sided λ≥0, "
              "split=0 — qualitatively DIFFERENT from ultrahyperbolic (no two-sidedness/splitting))")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3: wrong measure claim (q=3 imaginary) -> actually real
    if measure_real(3) and not measure_real(2):
        print("  MUTANT M3: CAUGHT (false claim «√−g is imaginary at q=3»: it is actually REAL "
              "(q=3 is odd); imaginary only for even q)")
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4: false branch (claim q=1 splits) -> (3,1) split=0
    sp31, _ = branch_split_count(3, 1, 4, tab4)
    if sp31 == 0:
        print("  MUTANT M4: CAUGHT (false-branch «q=1 splits»: (3,1) split=0 — "
              "a single iε-traversal, the branch is unambiguous)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): Euclidean (n,0) — a non-ultrahyperbolic reference")
    random.seed(1005031)
    p_nc = random.choice([3, 4])
    spec_nc = spectrum(p_nc, 0, 4, tab4)
    nneg_nc = sum(1 for x in spec_nc if x < 0)
    sp_nc, _ = branch_split_count(p_nc, 0, 4, tab4)
    ok(nneg_nc == 0 and sp_nc == 0,
       "control: Euclidean ({0},0) λ≥0, split=0 (non-ultrahyperbolic)".format(p_nc))
    print("  ({0},0): nneg={1} (no negatives), split-modes={2} — purely one-sided, no "
          "ultrahyperbolic signal; the measurement is sensitive (signature q governs two-sidedness and splitting)"
          .format(p_nc, nneg_nc, sp_nc))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

    # ==================== FORBIDDEN-SCAN ====================
    _pp = [("причи", "нн"), ("всес", "віт"), ("ант", "роп"), ("Teg", "mark")]  # GUARDLINE
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hits_src = scan_forbidden(__file__, _PATTERNS)
    _logf.flush()
    _logtxt = "".join(_tee.chunks)
    _hits_log = scan_forbidden(_logtxt, _PATTERNS)
    _nhits = len(_hits_src) + len(_hits_log)
    print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2})".format(
        _nhits, len(_hits_src), len(_hits_log)))

    _exit = 0
    if _nhits > 0:
        _exit = 1
    if FAILS[0] > 0:
        _exit = 1
    if not mut_ok:
        _exit = 1
    print("EXIT={0}".format(_exit))
    print("PROC_EXIT={0}".format(_exit))
    print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
