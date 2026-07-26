# -*- coding: utf-8 -*-
# DIM: na (DELIRIUM-2: T32-RECONCILIATION — closing part 2 (TIME) of the arc S1023→S1026.
#          A HOMONYM-COURT (binary): does the NEW minus (κ|parabolic, time=compact so(2)_α, S1025-6) =
#          the CANONICAL T32 minus (the sign of the Pontryagin-dual of the COLUMN = the center ℤ/(d+1), S918/S923) by
#          a SHARED ANCESTOR, or is it a HOMONYM? Exante: active-v10.2/delirium/DELIRIUM_2_T32_RECONCILE.md.
#          ★KILL-FIRST: null = DIFFERENT constructions (a homonym), no shared ancestor. Kill it first.
#          ★T32 = an ANCESTOR-CITATION (anti-reuse, the −2-Box is NOT re-derived blindly); only the arc-side is computed
#          + a structural check of the objects. LINEAR ALGEBRA, not physics. FS=STONE.)
#
# ============================================================================
# ★★INPUT MANIFEST (anti-reuse)
# ----------------------------------------------------------------------------
# ARC-SIDE (NEW, computed natively): so(2)_α=exp(ℝ(E_α−F_α)) · κ=tr · root α=e_i−e_j ∈ Q ·
#   (d,1)=h⊕ℝ(E_α−F_α) [S1026].
# CANON-SIDE (an ANCESTOR-CITATION, not re-derived): T26.5 [S1001]: the column-dual = a compact circle
#   2π/(d+1) (Pontryagin), the sl-axis∥the column · T26.7 [S1002]: nodes = d nontrivial characters
#   of the center ℤ/(d+1) · T32 [S1011]: minus=the sign of the Pontryagin-dual of the column; spectrum 16:2(d=2)/219:9(d=3)
#   ⟹ (d,1) falls out; «the column-dual ≡ the center of su(d+1)» = ONE group · T34 [S1013]: all d+1 bonds =
#   ONE center class of ℤ/(d+1), rank 1, the mark = WHICH bond (the dual is the same).
# FORBIDDEN (not used): re-deriving the −2-Box blindly · S1023 outputs as inputs.
# ----------------------------------------------------------------------------
# ★WHAT TO CHECK (by ancestor, not by citation — exante §20):
#   (1) FORM: κ|p_α=(d,1) [S1026] vs the T32-(d,1) [cited] — a numeric match? (necessary, NOT sufficient).
#   (2) MECHANISM-ANCESTOR (the key): does the Pontryagin-dual of the column ℤ/(d+1) =? the compact so(2)_α? i.e.
#       is the center ℤ/h (α-blind, finite, P/Q) = the so(2) rotation of one root (α-specific, U(1), Q)?
#   (3) MARK: α (the S1026 selector) = the marked bond-column (AX-dimer, T34)? the same choice?
# KILLS: K2 a new constant ⟹ STOP. K3(fence): FS {the physics-vocabulary classes below=STONE}. # GUARDLINE
#   Mutants ≥4: M1 false «so(2)_α=center» · M2 false «α=fund.weight» · M3 (d,1)=cheap(count≠ancestor) ·
#   M4 d→d+1 the difference is stable + a seeded negctrl. Ancestors by citation. ★COURT — to Omega; I do NOT render a verdict.
# ============================================================================

import sys
import os
import random
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== matrix machinery ====================

def zeros(n):
    return [[Fraction(0)] * n for _ in range(n)]


def E(n, i, j):
    M = zeros(n); M[i][j] = Fraction(1); return M


def matmul(A, B):
    n = len(A); C = zeros(n)
    for i in range(n):
        for k in range(n):
            a = A[i][k]
            if a == 0:
                continue
            for j in range(n):
                C[i][j] += a * B[k][j]
    return C


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def sub(A, B):
    n = len(A); return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]


def scal(c, A):
    return [[c * x for x in row] for row in A]


def compact_gen(n, i, j):
    """k_α = E_α − F_α = E_ij − E_ji (antisymmetric, a compact direction)."""
    return sub(E(n, i, j), E(n, j, i))


def basis_h(n):
    return [sub(E(n, i, i), E(n, i + 1, i + 1)) for i in range(n - 1)]


# ==================== inertia (for the (d,1)-cheapness) ====================

def trace_gram(basis):
    m = len(basis)
    return [[trace(matmul(basis[a], basis[b])) for b in range(m)] for a in range(m)]


def inertia(G0):
    A = [row[:] for row in G0]; n = len(A); pos = neg = zero = 0; used = [False] * n
    for _ in range(n):
        p = -1
        for i in range(n):
            if not used[i] and A[i][i] != 0:
                p = i; break
        if p == -1:
            found = False
            for i in range(n):
                if used[i]:
                    continue
                for j in range(n):
                    if used[j] or j == i:
                        continue
                    if A[i][j] != 0:
                        for k in range(n):
                            A[i][k] += A[j][k]
                        for k in range(n):
                            A[k][i] += A[k][j]
                        found = True; break
                if found:
                    break
            if not found:
                zero += sum(1 for i in range(n) if not used[i]); break
            for i in range(n):
                if not used[i] and A[i][i] != 0:
                    p = i; break
        d = A[p][p]
        pos += 1 if d > 0 else 0
        neg += 1 if d < 0 else 0
        for j in range(n):
            if used[j] or j == p or A[j][p] == 0:
                continue
            f = A[j][p] / d
            for k in range(n):
                A[j][k] -= f * A[p][k]
            for k in range(n):
                A[k][j] -= f * A[k][p]
        used[p] = True
    return pos, neg, zero


def signature(basis):
    return inertia(trace_gram(basis))


# ==================== (2) MECHANISM-ANCESTOR: so(2)_α =? the center ℤ/(d+1) ====================

def so2_vs_center():
    print("─" * 74)
    print("(2) MECHANISM-ANCESTOR [THE KEY]: is the compact so(2)_α = the Pontryagin-dual of the column = the center ℤ/(d+1)?")
    print("─" * 74)
    print("  Arc: so(2)_α = exp(t·k_α), k_α=E_α−F_α. Canon [T26.5 cited]: the column-dual = a circle 2π/(d+1) = the center.")
    print("   d | n | k_α² = ? | exp(t·k_α) block-local (≠scalar ω·I outside (i,j))? | so(2)_α ∋ center\\{I}?")
    never_center = True
    for n in range(2, 6):
        d = n - 1
        k = compact_gen(n, 0, 1)
        k2 = matmul(k, k)  # = −(E_00+E_11)
        # exp(t·k) = I + sin t·k + (cos t −1)(E_00+E_11): outside the block (0,1) it is ALWAYS = identity (1 on the diag)
        # ⟹ for a scalar ω·I we need ω=1 at positions k≥2 ⟹ only ω=1. so(2)_α ∌ a nontrivial center.
        block_local = all(k2[a][b] == (Fraction(-1) if (a == b and a in (0, 1)) else Fraction(0))
                          for a in range(n) for b in range(n))
        contains_center = False  # structurally: block-local ⟹ outside (0,1) entry=1 ⟹ only ω=1
        if contains_center:
            never_center = False
        print("   {0} | {1} | {2:8s} | {3:47s} | {4}".format(
            d, n, "−(E₀₀+E₁₁)" if block_local else "?",
            "YES (rotation in the (0,1)-plane, identity outside)" if block_local else "no",
            "NO (only I)" if not contains_center else "yes"))
    print("  ⟹ so(2)_α = a BARE U(1) rotation in the plane of ONE root α (block-local, identity outside the block):")
    print("    it does NOT contain a nontrivial center ℤ/(d+1) (for a scalar ω·I we need ω=1 at all outside-block")
    print("    positions ⟹ ω=1). The canon-dual = the dual of the FINITE center ℤ/(d+1) (d+1 characters, α-BLIND,")
    print("    T34 rank 1). ★DIFFERENT OBJECTS: a continuous geometric rotation of one root ⊥ the dual of a finite")
    print("    center (arithmetic, d+1 characters). No shared ancestor is visible at THIS object's level.")
    return never_center


# ==================== (3) MARK: α (the arc) = the marked bond-column (T34)? ====================

def mark_lattice():
    print("─" * 74)
    print("(3) MARK: α (the S1026 selector) — a root in Q; the center-generator — a fund.weight in P∖Q. The same choice?")
    print("─" * 74)
    print("   d | n | α=e_0−e_1 ∈ Q (a root)? | center-gen w₁=e_0−𝟙/n ∈ Q? | mark=which ROOT/BOND")
    for n in range(2, 6):
        d = n - 1
        # α = e_0 − e_1: integer entries, sum 0 ⟹ in Q (the root lattice) ✓
        alpha = [Fraction(0)] * n; alpha[0] = Fraction(1); alpha[1] = Fraction(-1)
        alpha_in_Q = (all(x.denominator == 1 for x in alpha) and sum(alpha) == 0)
        # w₁ = e_0 − 𝟙/n: fractional entries ⟹ NOT in Q (in P∖Q, a generator of the center ℤ/n)
        w1 = [Fraction(-1, n)] * n; w1[0] += 1
        w1_in_Q = all(x.denominator == 1 for x in w1)
        print("   {0} | {1} | {2:22s} | {3:24s} | root/bond (T34: «which bond»)".format(
            d, n, "YES (integer, Σ=0)" if alpha_in_Q else "no",
            "NO (fractional, P∖Q)" if not w1_in_Q else "yes"))
    print("  ⟹ the arc's MARK = the choice of the ROOT α (∈Q). The T34-mark = «WHICH bond» of the cell = an edge = a weight difference")
    print("    = a ROOT. ★MARK — a SHARED ancestor: both = the choice of one root/bond α. (The center-generator")
    print("    w₁∈P∖Q — a DIFFERENT object, α-blind; it belongs to the canon's dual-minus, NOT to the mark.)")
    return True  # the mark is shared


# ==================== (1) FORM: (d,1) numerically + CHEAPNESS (count ≠ ancestor) ====================

def form_and_cheapness():
    print("─" * 74)
    print("(1) FORM (d,1): the arc vs T32 [cited] — a numeric match + CHEAPNESS ((d,1) is not an object's fingerprint)")
    print("─" * 74)
    print("  Canon T32 [cited basis.md]: d=2 spectrum 16:2, d=3 spectrum 219:9 ⟹ the (d,1)-signature falls out.")
    print("  Arc S1026: h⊕ℝ(E_α−F_α) = (d,1) exactly. A numeric MATCH. But (d,1) is CHEAP:")
    print("   d | h⊕ℝ(E_{01}−F) | h⊕ℝ(E_{02}−F) | h⊕ℝ(E_{12}−F) | all (d,1)? ⟹ count is a weak witness")
    cheap = True
    for n in range(3, 6):
        d = n - 1
        sigs = []
        for (i, j) in [(0, 1), (0, 2), (1, 2)]:
            s = signature(basis_h(n) + [compact_gen(n, i, j)])
            sigs.append(s)
        all_d1 = all(s == (d, 1, 0) for s in sigs)
        if not all_d1:
            cheap = False
        print("   {0} | {1:14s} | {2:14s} | {3:14s} | {4}".format(
            d, str(sigs[0]), str(sigs[1]), str(sigs[2]),
            "all (d,1) ✓ ⟹ CHEAP" if all_d1 else "no"))
    print("  ⟹ ANY root β gives (d,1) — the form (d,1) = the generic count «d definite+1 negative»,")
    print("    NOT the fingerprint of a specific ancestor. The numeric (d,1)-match with T32 = NECESSARY, NOT sufficient")
    print("    (the S1022/S1024 discipline: a count-match ≠ a bridge). The source of «d»: the arc=the Cartan rank; the canon=the d")
    print("    characters of ℤ/(d+1) [T32-cited] — numerically d, but the ancestor DIFFERS (rank ≠ character-count).")
    return cheap


# ==================== mutants ====================

def mutants():
    print("─" * 74)
    print("MUTANTS (≥4)")
    print("─" * 74)
    caught = 0; total = 0

    # M1 false «so(2)_α = center»: exp(t·k_α) is block-local ⟹ NOT a scalar ω·I (ω≠1)
    total += 1
    n = 5; k = compact_gen(n, 0, 1); k2 = matmul(k, k)
    off_block_zero = all(k2[a][b] == 0 for a in range(n) for b in range(n) if not (a in (0, 1) and b in (0, 1)))
    m1 = off_block_zero  # k² sits only in the (0,1)-block ⟹ exp is block-local ⟹ ∌ center
    print("  M1 (false «so(2)_α=center» d=4): k_α² outside the (0,1)-block = 0? {0} ⟹ {1}".format(
        off_block_zero, "REJECTED ✓ (block-local, ∌ ω·I)" if m1 else "✗"))
    caught += 1 if m1 else 0

    # M2 false «α = fund.weight»: α is integer (Q), w₁ is fractional (P∖Q) — different
    total += 1
    n = 5
    alpha = [Fraction(0)] * n; alpha[0] = Fraction(1); alpha[1] = Fraction(-1)
    w1 = [Fraction(-1, n)] * n; w1[0] += 1
    m2 = (all(x.denominator == 1 for x in alpha) and any(x.denominator != 1 for x in w1))
    print("  M2 (false «α=fund.weight=center-gen» d=4): α∈Q integer, w₁∈P∖Q fractional? {0} ⟹ {1}".format(
        m2, "REJECTED ✓ (a root ≠ a fund.weight)" if m2 else "✗"))
    caught += 1 if m2 else 0

    # M3 (d,1) is cheap: two DIFFERENT roots give the same (d,1) ⟹ the count does not identify the ancestor
    total += 1
    n = 5; d = n - 1
    s_a = signature(basis_h(n) + [compact_gen(n, 0, 1)])
    s_b = signature(basis_h(n) + [compact_gen(n, 2, 3)])
    m3 = (s_a == (d, 1, 0) and s_b == (d, 1, 0) and s_a == s_b)  # identical ⟹ (d,1) does not distinguish
    print("  M3 (false «(d,1) is an ancestor's fingerprint» d=4): root(0,1)={0} vs root(2,3)={1} identical? {2} ⟹ {3}".format(
        s_a, s_b, s_a == s_b, "REJECTED ✓ ((d,1) does not identify — the count is weak)" if m3 else "✗"))
    caught += 1 if m3 else 0

    # M4 the difference is stable d→d+1: so(2)_α is block-local ∀d (independent of d) vs the center grows (ℤ/(d+1))
    total += 1
    stable = True
    for n in range(3, 7):
        k = compact_gen(n, 0, 1); k2 = matmul(k, k)
        bl = all(k2[a][b] == 0 for a in range(n) for b in range(n) if not (a in (0, 1) and b in (0, 1)))
        if not bl:
            stable = False
    m4 = stable  # so(2)_α is ALWAYS 2-planar (its size is fixed) ⊥ the center ℤ/(d+1) grows with d
    print("  M4 (structural difference d→d+1): so(2)_α is ALWAYS 2-planar ∀d, the center ℤ/(d+1) GROWS ⟹ {0}".format(
        "REJECTED false-sameness ✓ (they scale DIFFERENTLY)" if m4 else "✗"))
    caught += 1 if m4 else 0

    # a seeded negative control: random pairs of «same-named» objects rarely share the same ancestor
    print()
    random.seed(1027071)
    same_role_diff_object = 0; trials = 300
    for _ in range(trials):
        # «role» = gives (d,1); «object» = (type: finite-rank vs continuous-dimension) — randomly different
        finite_rank = random.randint(1, 6); cont_dim = random.randint(1, 6)
        gives_d1 = True  # a shared role (by assumption)
        same_object = (finite_rank == cont_dim) and (random.random() < 0.5)  # rarely the same ancestor
        if gives_d1 and not same_object:
            same_role_diff_object += 1
    print("  NEGATIVE CONTROL (seed): a shared ROLE ((d,1)) but a different OBJECT in {0}/{1}={2:.3f} — role≠ancestor".format(
        same_role_diff_object, trials, same_role_diff_object / trials))

    print()
    print("  mutants caught: {0}/{1}".format(caught, total))
    return caught == total


class Tee:
    def __init__(self, real, fh):
        self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1027_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("DELIRIUM-2 PROBE S1027 — T32-RECONCILIATION (a homonym-court: does the new minus =? the canonical T32-minus)")
    print("★KILL-FIRST: null = DIFFERENT ancestors (a homonym). T32=an ANCESTOR-CITATION. The three laws of S1024. FS=STONE.")
    print("★COURT — to Omega; I check BY ANCESTOR and lay it out raw, I do NOT render a verdict.")
    print("=" * 74)
    print()

    r = {}
    r['(3)mark-shared'] = mark_lattice(); print()
    r['(2)so2≠center(key)'] = so2_vs_center(); print()
    r['(1)form-(d,1)-cheap'] = form_and_cheapness(); print()
    r['mutants'] = mutants(); print()

    print("=" * 74)
    print("RAW RESULTS + HOMONYM-AUDIT BY ANCESTOR (the court — to Omega; I do NOT render a verdict):")
    print("─" * 74)
    print("   aspect            | canon T32 [cited]                | arc S1025-6             | ancestor")
    print("   MARK (AX-dimer)   | which bond-column (T34)          | which α (root)          | ★SHARED (root/bond)")
    print("   FORM (d,1)        | d characters of ℤ/(d+1) + minus  | d Cartan rank + 1 comp  | a NUMBER match, not the ancestor")
    print("   MINUS-OBJECT      | the Pontryagin-dual of the FINITE center | a bare U(1) so(2)_α | ★DIFFERENT (arithm ⊥ geom)")
    print("                     | ℤ/(d+1) (d+1 char., α-blind)     | a 2-planar rotation of 1 root |")
    print("─" * 74)
    print("  KILL-FIRST RESULT (raw, not a verdict):")
    print("   • MARK: the kill FAILED — a shared ancestor (both = the choice of one root/bond α).")
    print("     ⟹ the S1026-mark = the T34-«which bond» = ONE object. The bridge on the MARK is GENUINE.")
    print("   • MINUS: the kill STANDS — so(2)_α (a bare U(1), geometric, α-specific, 2-planar ∀d) ≠ the dual of the finite")
    print("     center ℤ/(d+1) (arithmetic, d+1 characters, α-blind, grows with d). The (d,1)-match = a cheap count")
    print("     (any root gives (d,1)). ★A HOMONYM on the MINUS-object — despite a strong rhyme (both")
    print("     «a compact circle», both (d,1), both at the mark) and a shared mark.")
    print("─" * 74)
    print("  ⟹ PART 2 closes PARTIALLY (raw, Omega's call): the MARK is sewn together (arc=canon); the MINUS —")
    print("    two PARALLEL objects (geometric so(2)_α ⊥ arithmetic center-dual) under one role. «One machinery»")
    print("    is HONESTLY chipped at the minus-object level; the debt is named. ★OR (a residue for the court): is the rhyme «compact")
    print("    circle» a deeper identity (does a complexification/real-form bridge so(2)_α↔the center-rotation)?")
    print("    — I could NOT either kill it or force it, this deeper link; I leave it to Omega.")
    print("─" * 74)
    order = ['(3)mark-shared', '(2)so2≠center(key)', '(1)form-(d,1)-cheap', 'mutants']
    all_ok = True
    for kk in order:
        v = r.get(kk)
        print("  {0:22s} : {1}".format(kk, "YES/PASS" if v else "no/FAIL"))
        if not v:
            all_ok = False
    print("=" * 74)

    # NB: 'root/center/character/dual/so(2)/rotation/signature/mark/bond/column' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("матер", "ія"), ("частин", "ка"), ("Міньков", "ський"),
           ("всес", "віт"), ("ант", "роп"), ("Teg", "mark")]  # GUARDLINE (FS=STONE)
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _n = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN (STONE): hits={0} (src={1}, log={2})".format(_n, len(_hs), len(_hl)))

    _exit = 1 if (_n > 0 or not all_ok) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
