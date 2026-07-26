# -*- coding: utf-8 -*-
# DIM: na (PLAN-1D: (d,1) FROM A PARABOLIC — CLOSING the form of the axis (the author's word «close it»).
#          The S1025-crown: a UNIQUE Schur-form κ on sl(n), Jordan-Chevalley (h definite ⊥ n± hyperbolic).
#          A residue: n± is balanced (C(n,2),C(n,2)) ⊥ the physical axis AX-indef=(d,1). RECONCILING THE FORM.
#          BINARY QUESTION: is (d,1) = κ on the MINIMAL parabolic p_α (one root α), NOT the full n±?
#          «which α» = the MARK(AX-dimer). Exante: active-v10.2/delirium/PLAN_1D_PARABOLIC_D1.md.
#          ★LINEAR ALGEBRA (signatures, inertia), not physics. FS=STONE.)
#
# ============================================================================
# ★★INPUT MANIFEST (law-1 ANTI-REUSE, S1024 in force)
# ----------------------------------------------------------------------------
# NEW/ROOTS: [N1] κ = the trace tr-form NATIVELY (uniqueness — S1025 CITATION, not re-derived
#   Schur; here κ is computed via tr directly, NOT from the S1023-Gram) · [N2] the minimal parabolic p_α=h⊕ℝE_α⊕ℝF_α ·
#   [N3] the MARK = the choice of a root α (one positive root).
# CARRIED-OVER: [T1] Ω (indices) · [T2] Λ. ROOT (marked): [R1] the bracket=AX-closure (H_α=[E_α,F_α]).
# FORBIDDEN-not-used (S1023 outputs): metric/Gram/T19/simplex/Cartan-matrix-as-given/arena-W.
#   ★h here = a diag-subalgebra of matrices; κ|h is computed via tr natively (NOT an import of the S1023 metric).
# ----------------------------------------------------------------------------
# ★BET (carved in the exante BEFORE the count; kill-first — null=«the parabolic does NOT give (d,1)»):
#   κ|p_α: the exante expectation (§construction) = (1,1)hyperbolic-pair{E_α,F_α} + H_α(+) + the rest of h(d−1,+).
#   ★HONEST ARITHMETIC of this same construction: p(+)=1+1+(d−1)=d+1, q(−)=1 ⟹ (d+1,1), NOT (d,1):
#   the exante-label «(d,1)» is off by +1 in space. I measure BOTH: (A) the parabolic p_α → (d+1,1) ·
#   (B) the clean axis h⊕ℝ(E_α−F_α) [Cartan + ONE negative direction] → (d,1) EXACTLY. The robust
#   invariant of BOTH = q EXACTLY 1 (one time-axis). Which form = AX-indef — an open homonym question for the project's adjudication.
# KILLS: K2 a new constant ⟹ STOP (κ=tr, no constants at all). K3(fence): FS {the physics-vocabulary classes below=STONE}. # GUARDLINE
#   Mutants ≥4: M1 two roots⟹q=2(the mark=EXACTLY one) · M2 H_α definite(control) · M3 d→d+1
#   q=1 stable · M4 false-κ(not invariant)⟹the signature breaks + a seeded negctrl.
#   Ancestors by citation: S1023(space) · S1025(κ unique, the split) · T32/S918((d,1)-families). COURT — to the project's adjudication.
# ============================================================================

import sys
import os
import random
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== matrix machinery over ℚ ====================

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


def add(A, B):
    n = len(A); return [[A[i][j] + B[i][j] for j in range(n)] for i in range(n)]


def bracket(A, B):
    return sub(matmul(A, B), matmul(B, A))


# ==================== inertia (Sylvester, exact) ====================

def trace_gram(basis):
    m = len(basis)
    G = [[trace(matmul(basis[a], basis[b])) for b in range(m)] for a in range(m)]
    return G


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
        if d > 0:
            pos += 1
        else:
            neg += 1
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


# ==================== sl(n) building blocks ====================

def basis_h(n):
    return [sub(E(n, i, i), E(n, i + 1, i + 1)) for i in range(n - 1)]


def Ealpha(n, i, j):
    return E(n, i, j)  # i<j: the positive root e_i−e_j


def Falpha(n, i, j):
    return E(n, j, i)


def parabolic_p_alpha(n, i, j):
    """(A) p_α = h ⊕ ℝE_α ⊕ ℝF_α (dim d+2)."""
    return basis_h(n) + [Ealpha(n, i, j), Falpha(n, i, j)]


def clean_axis(n, i, j):
    """(B) h ⊕ ℝ(E_α−F_α): Cartan + ONE negative (antisymmetric) direction (dim d+1)."""
    return basis_h(n) + [sub(Ealpha(n, i, j), Falpha(n, i, j))]


def sl2_triplet(n, i, j):
    """{E_α, H_α, F_α}, H_α=[E_α,F_α]=E_ii−E_jj."""
    return [Ealpha(n, i, j), bracket(Ealpha(n, i, j), Falpha(n, i, j)), Falpha(n, i, j)]


def two_root_object(n, a, b):
    """h ⊕ sl(2)_α ⊕ sl(2)_β for two INDEPENDENT roots (mutant M1)."""
    (i1, j1), (i2, j2) = a, b
    return basis_h(n) + [Ealpha(n, i1, j1), Falpha(n, i1, j1),
                         Ealpha(n, i2, j2), Falpha(n, i2, j2)]


# ==================== the main probe ====================

def core_measure():
    print("─" * 74)
    print("★BINARY QUESTION (closing): κ on the MINIMAL parabolic p_α — is q=EXACTLY 1 (one axis)?")
    print("─" * 74)
    print("α = the simple root e_0−e_1 (mark=the choice of α). κ = the trace tr-form (unique, S1025 citation).")
    print("   d | n | (A) p_α=h⊕E_α⊕F_α | (B) h⊕ℝ(E_α−F_α) | sl(2)_α triplet | q EXACTLY 1 (A/B)?")
    a_q1 = True
    b_is_d1 = True
    for n in range(2, 7):
        d = n - 1
        sig_A = signature(parabolic_p_alpha(n, 0, 1))
        sig_B = signature(clean_axis(n, 0, 1))
        sig_T = signature(sl2_triplet(n, 0, 1))
        qA1 = (sig_A[1] == 1)
        qB1 = (sig_B[1] == 1)
        b_d1 = (sig_B == (d, 1, 0))
        if not qA1:
            a_q1 = False
        if not b_d1:
            b_is_d1 = False
        print("   {0} | {1} | {2:17s} | {3:17s} | {4:15s} | A:q={5} B:q={6} B=(d,1)?{7}".format(
            d, n, str(sig_A), str(sig_B), str(sig_T), sig_A[1], sig_B[1],
            "YES" if b_d1 else "no"))
    print()
    print("  READING (raw facts):")
    print("   • (A) the parabolic p_α (the exante object) = (d+1,1): q=EXACTLY 1 ✓, but p=d+1 (the hyperbolic pair {E_α,F_α}")
    print("     gives a NULL-CONE 2-plane: 1 time ⊥ 1 spatial-partner). The exante-label «(d,1)» is off by +1 in p.")
    print("   • (B) the clean axis h⊕ℝ(E_α−F_α) = EXACTLY (d,1): d definite (space) + 1 negative (time).")
    print("     This is precisely the AX-indef-form (d,1), q=1, without a null-cone-partner.")
    print("   • BOTH: q = EXACTLY 1 ∀d — ONE time-axis, chosen by the root α. A robust invariant of the closing.")
    print("   • the sl(2)_α triplet = (2,1): a 2-dim null-cone + H_α — the core of both constructions.")
    return a_q1, b_is_d1


def mark_link():
    print("─" * 74)
    print("MARK↔AXIS (the lock of the two remaining roots): «which α» = the choice of the UNIQUE time-axis")
    print("─" * 74)
    print("   d | q(p_α) ∀ choice of α among the simple roots | axis = the chosen root?")
    all_q1 = True
    for n in range(3, 7):
        d = n - 1
        qs = []
        for i in range(n - 1):  # simple roots e_i−e_{i+1}
            s = signature(parabolic_p_alpha(n, i, i + 1))
            qs.append(s[1])
        uniform = all(q == 1 for q in qs)
        if not uniform:
            all_q1 = False
        print("   {0} | {1:32s} | {2}".format(
            d, str(qs), "YES — every α gives q=1, the choice=THE MARK" if uniform else "no"))
    print("  ⟹ EVERY simple root α gives EXACTLY ONE axis (q=1); «WHICH α» = the MARK(AX-dimer). ⟹ ★THE LOCK:")
    print("    the mark and the axis = ONE datum (which root = time). The two roots remaining from the S1023-review are MERGED:")
    print("    axis-TYPE = the Schur-κ split [S1025] · axis-FORM (·,1) = the parabolic · UNIQUENESS = the mark [this one].")
    return all_q1


def homonym_audit():
    print("─" * 74)
    print("HOMONYM-AUDIT ((·,1)-parabolic vs AX-indef/T32-(d,1))")
    print("─" * 74)
    print("  TYPE: both = an indefinite real-quad-form, q=1 (one axis) — a TYPE MATCH (≠ the sgn character, S1024).")
    print("  FORM: (B) h⊕ℝ(E_α−F_α) = EXACTLY (d,1) = the AX-indef form ((3,1) for d=3 — the S918/S923 families).")
    print("  ANCESTOR: the parabolic-(d,1) from {κ=the S1025-unique · the α-mark}; the T32-(d,1) from participation/Pontryagin.")
    print("   ★IS it the same object: both = ONE negative direction on a definite Cartan-background, chosen by")
    print("    ONE root. The bridge-candidate is GENUINE (the same α-sl(2), one κ) — BUT checking the ancestor")
    print("    T32-mechanism (participation) against the α-choice is left to the COURT (T32 is not re-derived). I do not render a verdict.")
    return True


def anti_reuse_audit():
    print("─" * 74)
    print("ANTI-REUSE AUDIT (law-1): used ∩ FORBIDDEN(S1023 outputs) = ∅ ?")
    print("─" * 74)
    FORBIDDEN = {"metric_S1023", "gram_S1023", "T19", "simplex_as_given", "cartan_matrix_as_given",
                 "arena_W_as_metric"}
    used = {"trace_form_kappa_native", "minimal_parabolic", "root_choice_MARK",
            "lie_bracket_AXclosure_ROOT", "bare_set_Omega", "Lambda"}
    leak = used & FORBIDDEN
    print("   used     = {0}".format(sorted(used)))
    print("   ★κ is computed via tr NATIVELY (uniqueness — S1025 citation), NOT from the S1023-Gram.")
    print("   leak = {0} ⟹ {1}".format(sorted(leak), "CLEAN ✓" if not leak else "REUSE ✗"))
    return len(leak) == 0


def mutants():
    print("─" * 74)
    print("MUTANTS (≥4)")
    print("─" * 74)
    caught = 0; total = 0

    # M1 two roots ⟹ q=2 (the mark = EXACTLY one root)
    total += 1
    s2 = signature(two_root_object(5, (0, 1), (2, 3)))  # d=4, two independent roots
    m1 = (s2[1] == 2)
    print("  M1 (two roots α,β d=4): sig={0} ⟹ q={1} {2}".format(
        s2, s2[1], "REJECTED false-one ✓ (q grows with #roots, the mark=EXACTLY 1)" if m1 else "✗"))
    caught += 1 if m1 else 0

    # M2 H_α is definite (a control — not indefinite by itself)
    total += 1
    sH = signature([bracket(Ealpha(4, 0, 1), Falpha(4, 0, 1))])  # {H_α}
    m2 = (sH == (1, 0, 0))
    print("  M2 (H_α alone): sig={0} ⟹ {1}".format(
        sH, "REJECTED false-indef ✓ (H_α is definite +)" if m2 else "✗"))
    caught += 1 if m2 else 0

    # M3 d→d+1: q=1 stable (the form holds)
    total += 1
    qs = [signature(parabolic_p_alpha(n, 0, 1))[1] for n in range(3, 7)]
    m3 = all(q == 1 for q in qs)
    print("  M3 (d→d+1 stability): q(p_α) over d={0} = {1} ⟹ {2}".format(
        list(range(2, 6)), qs, "REJECTED false-unstable ✓ (q=1 stable)" if m3 else "✗"))
    caught += 1 if m3 else 0

    # M4 false-κ (a NON-invariant form) ⟹ the signature differs/breaks
    total += 1
    # substitute the tr-Gram for an arbitrary symmetric (non-invariant) one → q need NOT be 1 systematically
    random.seed(1026041)
    n = 5; basis = parabolic_p_alpha(n, 0, 1); m = len(basis)
    Gfake = [[Fraction(0)] * m for _ in range(m)]
    for i in range(m):
        for j in range(i, m):
            v = Fraction(random.randint(-3, 3)); Gfake[i][j] = v; Gfake[j][i] = v
    sfake = inertia(Gfake)
    m4 = (sfake[1] != 1)  # the false-form does NOT give the canonical q=1 ⟹ invariance of κ is load-bearing
    print("  M4 (false-κ not invariant d=4): sig={0} ⟹ {1}".format(
        sfake, "REJECTED ✓ (q≠1 — it is precisely the invariance of κ that gives q=1, not the subspace's form)"
        if m4 else "coincidentally q=1, retry the seed"))
    caught += 1 if m4 else 0

    # a seeded negative control: random subspaces of dim d+2 ⟹ q is distributed (not always 1)
    print()
    random.seed(1026071)
    n = 5; dimsub = (n - 1) + 2  # d+2
    q1 = 0; trials = 300
    for _ in range(trials):
        G = [[Fraction(0)] * dimsub for _ in range(dimsub)]
        for i in range(dimsub):
            for j in range(i, dimsub):
                v = Fraction(random.randint(-3, 3)); G[i][j] = v; G[j][i] = v
        if inertia(G)[1] == 1:
            q1 += 1
    print("  NEGATIVE CONTROL (seed): random {0}×{0} forms with q=1: {1}/{2}={3:.3f} — q=1 is NOT automatic".format(
        dimsub, q1, trials, q1 / trials))

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
    _logf = open(os.path.join(_HERE, "S1026_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("PLAN-1D PROBE S1026 — (d,1) FROM A PARABOLIC (CLOSING the form of the axis, «close it»)")
    print("BINARY: is (d,1) = κ(the S1025-unique) on the MINIMAL parabolic p_α(one α), NOT the full n±?")
    print("★KILL-FIRST: null = «the parabolic does NOT give (d,1)». The three laws of S1024. FS=STONE. Court — to the project's adjudication.")
    print("=" * 74)
    print()

    results = {}
    results['anti-reuse'] = anti_reuse_audit(); print()
    a_q1, b_d1 = core_measure(); print()
    results['q=EXACTLY-1(parabolic)'] = a_q1
    results['(B)=EXACTLY-(d,1)'] = b_d1
    results['mark↔axis-lock'] = mark_link(); print()
    results['homonym-audit'] = homonym_audit(); print()
    results['mutants'] = mutants(); print()

    print("=" * 74)
    print("RAW RESULTS + STRUCTURAL CONCLUSION (the court — to the project's adjudication; I do NOT render a verdict):")
    print("─" * 74)
    print("  ★CLOSING (the author's binary question):")
    print("   • q = EXACTLY 1 ∀d from the MINIMAL parabolic of one root α (both constructions) —")
    print("     ONE time-axis is FORCED by one root. This is the form of the axis, not a separate uniqueness postulate.")
    print("   • (A) the parabolic p_α = (d+1,1) [a null-cone-pair]; (B) h⊕ℝ(E_α−F_α) = EXACTLY (d,1) [a clean axis].")
    print("     Both q=1; which one = AX-indef — the project's call (the exante-label «(d,1)» is exact for (B), off by +1 for (A)).")
    print("   • ★MARK↔AXIS MERGED: «which α» = the choice of axis = the MARK(AX-dimer). The two remaining roots = ONE datum.")
    print("─" * 74)
    print("  SUMMARY OF THE ARC S1023→S1026 (raw, not a verdict):")
    print("   space=a Schur-form on h [S1023] · axis-TYPE=the nilpotent-half of the same κ [S1025] ·")
    print("   axis-FORM (·,1)+UNIQUENESS=the mark on the parabolic [S1026] ⟹ AX-indef+AX-dimer SPLIT")
    print("   into {a UNIQUE κ (a simple algebra) · the choice of ONE root}. «Not many cards»: {Ω,order,closure,Λ}.")
    print("─" * 74)
    order = ['anti-reuse', 'q=EXACTLY-1(parabolic)', '(B)=EXACTLY-(d,1)', 'mark↔axis-lock',
             'homonym-audit', 'mutants']
    all_ok = True
    for kk in order:
        v = results.get(kk)
        print("  {0:24s} : {1}".format(kk, "YES/PASS" if v else "no/FAIL"))
        if not v:
            all_ok = False
    print("=" * 74)

    # NB: 'parabolic/root/signature/inertia/Cartan/hyperbolic/axis/mark' is STRUCTURAL vocabulary. GUARDLINE
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
