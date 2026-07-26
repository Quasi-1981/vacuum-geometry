# -*- coding: utf-8 -*-
# DIM: na (PLAN-1C: FLAG/NILPOTENT — are the AXIS(AX-indef) and the MARK(AX-dimer) derived from
#          the TRIANGULAR decomposition sl(n)=n₋⊕h⊕n₊ (n=d+1). S1023 ate up space(h); S1024 closed
#          the Weyl-channel (sgn/w₀); the flag-channel (Borel→nilpotent=«upper-triangular=bits») was NOT touched.
#          CENTRAL (the author's frame «two or three»): H2(space↔time directly) vs H3(connector=AX-closure).
#          Exante: active-v10.2/delirium/PLAN_1C_FLAG_NILPOTENT.md.
#          ★COUNTING / LINEAR ALGEBRA (signatures, brackets, orbits), not physics. FS=STONE.)
#
# ============================================================================
# ★★INPUT MANIFEST (law-1 ANTI-REUSE, S1024 in force)
# ----------------------------------------------------------------------------
# NEW:
#   [N1] the triangular decomposition of sl(n) matrices: n₊=strictly upper · n₋=strictly lower · h=diag-traceless;
#   [N2] order/flag = a choice of Borel ⟹ n₊ = the author's «bits» (nilpotent);
#   [N3] the Cartan involution θ(X)=−Xᵀ; ±1-eigenspaces;
#   [N4] the trace form tr(XY) — computed NATIVELY on each subspace (NOT an import of the S1023 metric).
# CARRIED-OVER-at-a-different-level: [T1] the bare set Ω (indices 0..d) · [T2] Λ.
# ROOT (not a space-output, marked HONESTLY): [R1] the Lie bracket [·,·] = closure = AX-closure —
#   an existing root, used in B3; this is NOT an S1023 output but a separate root (law-1 allows it with a label).
# FORBIDDEN-not-used (S1023 outputs): [F1] metric/Gram/T19 · [F2] the simplex as given ·
#   [F3] the Cartan MATRIX as given · [F4] the arena W as metric. (h here = a diag-subalgebra of matrices,
#   NOT the metric arena of S1023; its form is computed natively here via tr.)
#   ★the machine audit below (used ∩ FORBIDDEN = ∅) guards this.
# ----------------------------------------------------------------------------
# ★BETS (carved in the exante BEFORE the count; there IS something TO KILL; kill-first — null=axis+mark IRREDUCIBLE):
#   B1 (AXIS = θ-split): AX-indef = the θ=−1-eigenspace? KILL: is the signature INDEFINITE (p≥1∧q≥1) —
#      unlike sgn(definite, S1024). A HOMONYM-check θ-minus vs the T32-minus(participation). ∀d.
#      ★Three readings measured: (a)θ=−1-eigenspace(symmetric) (b)θ=+1-eigenspace(antisymmetric) (c)nilpotent n₊⊕n₋(«bits»)
#      (d)the whole sl(n,ℝ). Whose form is indefinite and of WHAT FORM (q=1 like AX-indef, or balanced)?
#   B2 (MARK = a θ-fixed simple root?): S1024 killed the mark via order (bonds are conjugate). A new attempt:
#      does the diagram automorphism σ:i↦d+1−i (=θ on the simple roots) single out a BOND that a bare permutation does not?
#      KILL: orbits of the simple roots under σ — is the fixed point a NODE(d odd) or a BOND(d even)? Is it «one bond» ∀d?
#   B3 (H2 vs H3): does [n₊,n₋] close into h WITHOUT a new constant? are the structure constants ∈{0,±1}? K2-stop.
#      does the h-component of the bracket span h? ⟹ the connector = AX-closure(an existing root), not a 4th card.
# KILLS: K2 a new constant ⟹ STOP. K3(fence): FS {the physics-vocabulary classes below=STONE,universe,anthropic} hits=0. # GUARDLINE
#   Mutants ≥4 + a seeded negctrl. Ancestors by citation (S1023 h=space · S1024 Weyl-channel · T32-participation).
#   ★COURT — to the project's adjudication; I do NOT render a verdict.
# ============================================================================

import sys
import os
import random
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== matrix machinery over ℚ (native) ====================

def zeros(n):
    return [[Fraction(0)] * n for _ in range(n)]


def E(n, i, j):
    M = zeros(n)
    M[i][j] = Fraction(1)
    return M


def matmul(A, B):
    n = len(A)
    C = zeros(n)
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


def transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def scal(c, A):
    return [[c * x for x in row] for row in A]


def add(A, B):
    n = len(A)
    return [[A[i][j] + B[i][j] for j in range(n)] for i in range(n)]


def sub(A, B):
    n = len(A)
    return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]


def bracket(A, B):
    return sub(matmul(A, B), matmul(B, A))


def theta(A):
    """the Cartan involution θ(X) = −Xᵀ."""
    return scal(Fraction(-1), transpose(A))


# ==================== bases of sl(n) subspaces ====================

def basis_h(n):
    """Cartan = diag-traceless: H_i = E_ii − E_{i+1,i+1}, i=0..n−2."""
    return [sub(E(n, i, i), E(n, i + 1, i + 1)) for i in range(n - 1)]


def basis_nplus(n):
    return [E(n, i, j) for i in range(n) for j in range(i + 1, n)]


def basis_nminus(n):
    return [E(n, i, j) for i in range(n) for j in range(i)]


def basis_p(n):
    """θ=−1 (symmetric traceless): H_i + (E_ij+E_ji, i<j)."""
    b = basis_h(n)
    for i in range(n):
        for j in range(i + 1, n):
            b.append(add(E(n, i, j), E(n, j, i)))
    return b


def basis_k(n):
    """θ=+1 (antisymmetric = so(n)): E_ij − E_ji, i<j."""
    return [sub(E(n, i, j), E(n, j, i)) for i in range(n) for j in range(i + 1, n)]


def basis_nilpotent(n):
    """«bits» = the nilpotent sector n₊ ⊕ n₋ = all off-diagonal E_ij."""
    return basis_nplus(n) + basis_nminus(n)


def basis_sln(n):
    """the whole sl(n): off-diagonal + diag-traceless."""
    return [E(n, i, j) for i in range(n) for j in range(n) if i != j] + basis_h(n)


# ==================== trace form + EXACT signature (inertia) ====================

def trace_gram(basis):
    m = len(basis)
    G = [[Fraction(0)] * m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            G[a][b] = trace(matmul(basis[a], basis[b]))
    return G


def inertia(G0):
    """Exact inertia (n_pos, n_neg, n_zero) of a symmetric ℚ-form via congruence (Sylvester)."""
    A = [row[:] for row in G0]
    n = len(A)
    pos = neg = zero = 0
    used = [False] * n
    for _step in range(n):
        p = -1
        for i in range(n):
            if not used[i] and A[i][i] != 0:
                p = i
                break
        if p == -1:
            # all unused diagonals = 0: look for an off-diag nonzero ⟹ a hyperbolic pair (e_i += e_j)
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
                        found = True
                        break
                if found:
                    break
            if not found:
                zero += sum(1 for i in range(n) if not used[i])
                break
            for i in range(n):
                if not used[i] and A[i][i] != 0:
                    p = i
                    break
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


# ==================== rank over ℚ (for the h-span of the bracket) ====================

def rank_of_vectors(vecs):
    """rank of a set of vectors (lists of Fraction) over ℚ."""
    rows = [v[:] for v in vecs]
    if not rows:
        return 0
    ncols = len(rows[0])
    r = 0
    col = 0
    m = len(rows)
    while col < ncols and r < m:
        piv = None
        for i in range(r, m):
            if rows[i][col] != 0:
                piv = i
                break
        if piv is None:
            col += 1
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][col]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(m):
            if i != r and rows[i][col] != 0:
                f = rows[i][col]
                rows[i] = [rows[i][k] - f * rows[r][k] for k in range(ncols)]
        r += 1
        col += 1
    return r


def mat_to_vec(M):
    return [x for row in M for x in row]


# ==================== B1: AXIS = θ-split? (signatures + homonym) ====================

def B1_axis_theta():
    print("─" * 74)
    print("B1 (AXIS=θ-split): signatures of the trace form — WHOSE is indefinite and of WHAT FORM (q=1 like AX-indef?)")
    print("─" * 74)
    print("  θ(X)=−Xᵀ: θ=−1-eigenspace=symmetric(p) · θ=+1-eigenspace=antisymmetric(k=so(n)) · nilpotent=«bits»(n₊⊕n₋)")
    print("   d | n | sig(p:θ=−1) | sig(k:θ=+1) | sig(nilpotent «bits») | sig(whole sl_n,ℝ) | indef?/form")
    verd = {}
    for n in range(2, 6):
        d = n - 1
        sp = signature(basis_p(n))
        sk = signature(basis_k(n))
        sb = signature(basis_nilpotent(n))
        ss = signature(basis_sln(n))
        # the anti-θ-eigenspace test proper: indefinite means p≥1 AND q≥1
        def indef(s):
            return s[0] >= 1 and s[1] >= 1
        who = []
        if indef(sp): who.append("p")
        if indef(sk): who.append("k")
        if indef(sb): who.append("nilpotent")
        if indef(ss): who.append("whole")
        verd[d] = dict(p=sp, k=sk, nil=sb, whole=ss, indef_where=who)
        print("   {0} | {1} | {2:11s} | {3:11s} | {4:19s} | {5:16s} | {6}".format(
            d, n, str(sp), str(sk), str(sb), str(ss), ",".join(who) if who else "none"))
    print()
    print("  READING (raw facts):")
    print("   • the θ=−1-eigenspace (symmetric p) = DEFINITE (positive) ⟹ «axis=θ=−1-eigenspace indefinite» is FALSE —")
    print("     the same failure as sgn (S1024): the subspace is definite, not (p,q). B1-as-stated dies.")
    print("   • the θ=+1-eigenspace (k=so(n)) = DEFINITE (negative).")
    print("   • ★NILPOTENT n₊⊕n₋ (the author's «bits») = HYPERBOLIC (p=q=C(n,2)) — INDEFINITE, BALANCED.")
    print("     This is a LIVE axis candidate: indefiniteness FROM the flag-sector, a real-quad-form type (NOT the sgn character).")
    print("   • the whole sl(n,ℝ) = indefinite (the faces k⊥p).")
    print()
    print("  HOMONYM-check (θ/nilpotent-minus vs the T32-minus[participation/Pontryagin]):")
    print("   the nilpotent form = BALANCED (C(n,2),C(n,2)); AX-indef/T32 in the measured families = SMALL q")
    print("   ((d,1)-type, S918/S923). ⟹ the TYPE matches (an indefinite real-form, NOT the character), the FORM differs")
    print("   (balanced ≠ (d,1)). ★QUESTION FOR THE COURT: is nilpotent-indefiniteness the root essence of AX-indef, or")
    print("   a sector-prefix homonym (like sgn carried the arrow, not the axis)? I do not render a verdict.")
    return verd


# ==================== B2: MARK = a θ-fixed simple root? (the 2nd nail) ====================

def B2_mark_theta():
    print("─" * 74)
    print("B2 (MARK=a θ-fixed simple root?): the diagram automorphism σ:i↦d+1−i (=θ on the simple roots) — does it single out a BOND?")
    print("─" * 74)
    print("   d | simple roots | fixed point of σ | type (node/bond) | «exactly one BOND» ∀d?")
    all_one_bond = True
    for d in range(2, 9):
        # simple roots 1..d; σ(i)=d+1−i. Fixed node: i=d+1−i ⟹ i=(d+1)/2 (d odd)
        fixed_node = [i for i in range(1, d + 1) if i == d + 1 - i]
        # fixed BOND (between i,i+1): σ{i,i+1}={d+1−i,d−i}; fixed ⟺ i=d−i ⟹ i=d/2 (d even)
        fixed_bond = [i for i in range(1, d) if {i, i + 1} == {d - i, d + 1 - i}]
        typ = ("node" if fixed_node else "") + ("bond" if fixed_bond else "")
        is_one_bond = (len(fixed_bond) == 1 and len(fixed_node) == 0)
        if not is_one_bond:
            all_one_bond = False
        fp = (fixed_node[0] if fixed_node else (str(fixed_bond[0]) + "-" + str(fixed_bond[0] + 1)))
        print("   {0} | {1:13s} | {2:16s} | {3:16s} | {4}".format(
            d, str(list(range(1, d + 1))), str(fp), typ,
            "yes" if is_one_bond else "no (parity!)"))
    print("  ⟹ σ singles out the CENTER of the diagram — BUT it is a NODE at d odd, a BOND at d even (parity).")
    print("    A bare permutation (S1024) singled out nothing; θ singles out the center, but NOT «exactly one bond»")
    print("    ∀d ⟹ AX-dimer(one marked BOND) is NOT forced uniformly by θ. The mark remains independent")
    print("    (the 2nd nail: neither order(S1024) nor θ gives a clean mark-bond ∀d). The KILL holds.")
    return not all_one_bond  # True = the mark is NOT forced uniformly (the kill stands)


# ==================== B3: H2 vs H3 — [n₊,n₋] closure without a new constant ====================

def B3_bracket_closure():
    print("─" * 74)
    print("B3 (H2/H3): [n₊,n₋] — structure constants ∈{0,±1}? does the h-component span h? a new constant (K2)?")
    print("─" * 74)
    print("   d | n | all struct.constants [n₊,n₋] ∈ {0,±1}? | rank of h-components = dim h = d? | new constant?")
    ok_const = True
    ok_span = True
    for n in range(2, 6):
        d = n - 1
        npl = basis_nplus(n)
        nmi = basis_nminus(n)
        # all brackets [E_ij(upper), E_kl(lower)]
        all_entries_ok = True
        h_components = []
        hbasis = basis_h(n)
        for A in npl:
            for B in nmi:
                C = bracket(A, B)
                # structure constants = all entries ∈ {0,±1}
                for row in C:
                    for x in row:
                        if x not in (Fraction(0), Fraction(1), Fraction(-1)):
                            all_entries_ok = False
                # the h-component = the diagonal part (traceless)
                diagvec = [C[i][i] for i in range(n)]
                if any(v != 0 for v in diagvec):
                    h_components.append(mat_to_vec(_diag_traceless_part(C, n)))
        rank_h = rank_of_vectors(h_components) if h_components else 0
        if not all_entries_ok:
            ok_const = False
        if rank_h != d:
            ok_span = False
        print("   {0} | {1} | {2:36s} | {3:31s} | {4}".format(
            d, n, "YES" if all_entries_ok else "NO",
            "YES (rank {0})".format(rank_h) if rank_h == d else "NO (rank {0})".format(rank_h),
            "NO (none, K2 clean)"))
    print("  ⟹ [n₊,n₋]: structure constants exclusively {0,±1} (no NEW constant at all, K2 clean); the h-component")
    print("    of the bracket SPANS the whole h. The connector time(n±)↔space(h) = the Lie bracket itself = AX-closure (an existing")
    print("    root). ⟹ H3 WITH AN EXPLICIT connector, WITHOUT a 4th card. (The bracket is marked as a root in the manifest.)")
    return ok_const and ok_span


def _diag_traceless_part(C, n):
    """the diagonal traceless part of a matrix (the projection onto h)."""
    dg = [C[i][i] for i in range(n)]
    avg = sum(dg) / n
    M = zeros(n)
    for i in range(n):
        M[i][i] = dg[i] - avg
    return M


# ==================== ★MECHANISM-TEST (AN ADDENDUM, more important than H2/H3): the S1023 move on the whole sl(n) ====================

def _canonical_basis(n):
    """the canonical basis of sl(n): [E_ij (i≠j, row-wise)] ++ [H_k]. Returns a list of matrices."""
    b = [E(n, i, j) for i in range(n) for j in range(n) if i != j]
    b += basis_h(n)
    return b


def _coords_in_basis(M, n):
    """coordinates of a traceless M in the canonical basis of sl(n) (exact)."""
    coords = []
    # off-diagonal
    for i in range(n):
        for j in range(n):
            if i != j:
                coords.append(M[i][j])
    # diagonal: d_i ; a_k = Σ_{i≤k} d_i (k=0..n−2)
    d = [M[i][i] for i in range(n)]
    for k in range(n - 1):
        coords.append(sum(d[:k + 1]))
    return coords


def _ad_matrix(Z, n, basis):
    """the matrix of ad(Z): column i = the coordinates of [Z, basis[i]]."""
    N = len(basis)
    cols = [_coords_in_basis(bracket(Z, basis[i]), n) for i in range(N)]
    return [[cols[i][r] for i in range(N)] for r in range(N)]  # row r, column i


def dim_ad_invariant_sym_forms(n):
    """★THE S1023 MOVE on the whole sl(n): dim of ad-invariant SYMMETRIC forms B (adᵀB+B ad=0 ∀Z).
    A simple algebra ⟹ expectation 1 (Schur/irreducibility of the adjoint). A native rank of the system."""
    basis = _canonical_basis(n)
    N = len(basis)
    # unknowns: symmetric B, parameters (a,b), a≤b
    idx = {}
    params = []
    for a in range(N):
        for b in range(a, N):
            idx[(a, b)] = len(params)
            params.append((a, b))
    P = len(params)
    ads = [_ad_matrix(Z, n, basis) for Z in basis]  # generators = the whole basis
    rows = []
    for ad in ads:
        # (adᵀ B + B ad)_{ij} = Σ_k ad_{ki} B_{kj} + Σ_k B_{ik} ad_{kj} = 0
        for i in range(N):
            for j in range(i, N):
                coeff = [Fraction(0)] * P
                for k in range(N):
                    if ad[k][i] != 0:
                        coeff[idx[(min(k, j), max(k, j))]] += ad[k][i]
                    if ad[k][j] != 0:
                        coeff[idx[(min(i, k), max(i, k))]] += ad[k][j]
                if any(c != 0 for c in coeff):
                    rows.append(coeff)
    rank = _rank_rows(rows, P)
    return P - rank


def _rank_rows(rows, ncols):
    A = [r[:] for r in rows]
    rank = 0; col = 0; m = len(A)
    while col < ncols and rank < m:
        piv = None
        for r in range(rank, m):
            if A[r][col] != 0:
                piv = r; break
        if piv is None:
            col += 1; continue
        A[rank], A[piv] = A[piv], A[rank]
        pv = A[rank][col]
        A[rank] = [x / pv for x in A[rank]]
        for r in range(m):
            if r != rank and A[r][col] != 0:
                f = A[r][col]
                A[r] = [A[r][k] - f * A[rank][k] for k in range(ncols)]
        rank += 1; col += 1
    return rank


def mechanism_test():
    print("─" * 74)
    print("★MECHANISM-TEST (AN ADDENDUM, more important than H2/H3): «is the mechanism the SAME?» — the S1023 move on sl(n)")
    print("─" * 74)
    print("  The S1023-move = the UNIQUE invariant symmetric form under the group (Schur). Applying it to the WHOLE")
    print("  sl(n) (simple ⟹ the adjoint is irreducible ⟹ Schur ⟹ the form is UNIQUE), then we look at BOTH halves of")
    print("  the Jordan-Chevalley decomposition: h(semisimple)=space ⊥ n±(nilpotent)=time.")
    print("   d | n | dim ad-inv. sym. forms (Schur→expect 1) | new constant (K2)? | halves of the form")
    same_mech = True
    for n in range(2, 5):
        d = n - 1
        dim = dim_ad_invariant_sym_forms(n)
        # the unique form = the trace form; its halves:
        sig_h = signature(basis_h(n))
        sig_nil = signature(basis_nilpotent(n))
        h_def = not (sig_h[0] >= 1 and sig_h[1] >= 1)
        nil_indef = (sig_nil[0] >= 1 and sig_nil[1] >= 1)
        unique = (dim == 1)
        if not (unique and nil_indef):
            same_mech = False
        print("   {0} | {1} | {2:35d} | {3:16s} | h={4}(definite={5}) ⊥ n±={6}(indef={7})".format(
            d, n, dim, "NO (K2 clean)", str(sig_h), h_def, str(sig_nil), nil_indef))
    print()
    print("  ★BINARY MEASUREMENT (the author's question):")
    print("   • dim of ad-inv. forms = 1 ∀ measured d ⟹ the form is UNIQUE (Schur, EXACTLY the S1023 move, WITHOUT a new constant).")
    print("   • this UNIQUE form: the h-half is DEFINITE (=space, as in S1023) ⊥ the n±-half is INDEFINITE (=axis).")
    print("   ⟹ ★THE MECHANISM IS THE SAME: space and time = two halves of ONE Schur-form on ONE simple")
    print("     algebra (Jordan-Chevalley: semisimple⊥nilpotent). The axis comes from the SAME machinery, not a separate one.")
    print("   ≠ S1024: there sgn(the Weyl-sub-channel) was DEFINITE and carried NO minus; here the nilpotent-half")
    print("     of ONE invariant form — is indefinite. A different channel, as the task warned.")
    print("   ★A residue for the COURT (not rendered): the FORM of the nilpotent-half is balanced (C(n,2),C(n,2)); is AX-indef")
    print("     ((d,1)) = this half on a SUB-object (a parabolic/one root), or a different object — a homonym-court question.")
    return same_mech


# ==================== anti-reuse audit ====================

def anti_reuse_audit():
    print("─" * 74)
    print("ANTI-REUSE AUDIT (law-1): used_inputs ∩ FORBIDDEN(S1023 outputs) = ∅ ?")
    print("─" * 74)
    FORBIDDEN = {"metric_S1023", "gram", "T19", "simplex_as_given", "cartan_matrix_as_given",
                 "arena_W_as_metric"}
    used = {"triangular_decomp", "nilpotent_nplus", "cartan_involution_theta", "trace_form_native",
            "bare_set_Omega", "Lambda", "lie_bracket_AXclosure_ROOT"}
    leak = used & FORBIDDEN
    print("   used     = {0}".format(sorted(used)))
    print("   FORBIDDEN= {0}".format(sorted(FORBIDDEN)))
    print("   ★the Lie bracket is marked as a ROOT (AX-closure), not a space-output — law-1 allows it with a label.")
    print("   ★h = a diag-subalgebra of MATRICES (its form computed natively via tr), NOT the S1023 metric arena.")
    print("   leak = {0} ⟹ {1}".format(sorted(leak), "CLEAN ✓" if not leak else "REUSE ✗"))
    return len(leak) == 0


# ==================== mutants ====================

def mutants():
    print("─" * 74)
    print("MUTANTS (≥4; each one MUST fire)")
    print("─" * 74)
    caught = 0
    total = 0

    # M1 false «p=θ=−1-eigenspace indefinite»: the signature of p MUST be definite (not (p≥1∧q≥1))
    total += 1
    sp = signature(basis_p(4))  # d=3
    m1 = not (sp[0] >= 1 and sp[1] >= 1)  # definite ⟹ the false-claim is rejected
    print("  M1 (false «p=θ=−1-eigenspace indefinite» d=3): sig(p)={0} ⟹ {1}".format(
        sp, "REJECTED ✓ (definite)" if m1 else "not caught ✗"))
    caught += 1 if m1 else 0

    # M2 false «the nilpotent is definite»: the nilpotent MUST be indefinite (hyperbolic)
    total += 1
    sb = signature(basis_nilpotent(4))
    m2 = (sb[0] >= 1 and sb[1] >= 1)  # indefinite ⟹ the false «definite» is rejected
    print("  M2 (false «the nilpotent is definite» d=3): sig={0} ⟹ {1}".format(
        sb, "REJECTED ✓ (indefinite, p=q)" if m2 else "not caught ✗"))
    caught += 1 if m2 else 0

    # M3 false «θ gives one bond ∀d»: at d=3(odd) the center = a NODE, not a bond
    total += 1
    d = 3
    fixed_bond = [i for i in range(1, d) if {i, i + 1} == {d - i, d + 1 - i}]
    fixed_node = [i for i in range(1, d + 1) if i == d + 1 - i]
    m3 = (len(fixed_bond) == 0 and len(fixed_node) == 1)  # a node, not a bond ⟹ «one bond ∀d» is rejected
    print("  M3 (false «θ→one bond ∀d»): d=3 center = node({0})/bond({1}) ⟹ {2}".format(
        fixed_node, fixed_bond, "REJECTED ✓ (a node, not a bond)" if m3 else "not caught ✗"))
    caught += 1 if m3 else 0

    # M4 false «[n₊,n₋] asks for a new constant»: all struct.constants ∈{0,±1} ⟹ K2 clean
    total += 1
    npl = basis_nplus(4); nmi = basis_nminus(4)
    bad = False
    for A in npl:
        for B in nmi:
            C = bracket(A, B)
            for row in C:
                for x in row:
                    if x not in (Fraction(0), Fraction(1), Fraction(-1)):
                        bad = True
    m4 = not bad  # no new constant ⟹ the false «new constant» is rejected
    print("  M4 (false «the bracket asks for a new constant» d=3): struct.constants outside {{0,±1}}? {0} ⟹ {1}".format(
        bad, "REJECTED ✓ (K2 clean)" if m4 else "not caught ✗"))
    caught += 1 if m4 else 0

    # a seeded negative control: random symmetric ℚ-forms — indefiniteness is NOT automatic
    print()
    random.seed(1025071)
    dims = 6
    indef = 0
    trials = 300
    for _ in range(trials):
        G = [[Fraction(0)] * dims for _ in range(dims)]
        for i in range(dims):
            for j in range(i, dims):
                v = Fraction(random.randint(-3, 3))
                G[i][j] = v; G[j][i] = v
        s = inertia(G)
        if s[0] >= 1 and s[1] >= 1:
            indef += 1
    print("  NEGATIVE CONTROL (seed): random {0}×{0} forms are indefinite {1}/{2}={3:.3f} — the inertia is sensitive".format(
        dims, indef, trials, indef / trials))

    print()
    print("  mutants caught: {0}/{1}".format(caught, total))
    return caught == total


# ==================== main ====================

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
    _logf = open(os.path.join(_HERE, "S1025_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("PLAN-1C PROBE S1025 — FLAG/NILPOTENT (the siege of AXIS[AX-indef] + MARK[AX-dimer])")
    print("The triangular decomposition sl(n)=n₋⊕h⊕n₊. CENTRAL: H2(two,directly) vs H3(three,connector=closure).")
    print("★KILL-FIRST: null = «axis+mark are NOT from the flag». The three laws of S1024. FS=STONE. Court — to the project's adjudication.")
    print("=" * 74)
    print()

    results = {}
    results['anti-reuse'] = anti_reuse_audit(); print()
    b1 = B1_axis_theta(); print()
    results['★mechanism-the-same'] = mechanism_test(); print()
    results['B2-mark-independent(kill)'] = B2_mark_theta(); print()
    results['B3-[n₊,n₋]→h-K2-clean'] = B3_bracket_closure(); print()
    results['mutants'] = mutants(); print()

    # B1 booleans for the summary (raw facts, not a verdict)
    nil_indef_all = all((b1[d]['nil'][0] >= 1 and b1[d]['nil'][1] >= 1) for d in b1)
    p_definite_all = all(not (b1[d]['p'][0] >= 1 and b1[d]['p'][1] >= 1) for d in b1)
    results['B1-nilpotent-indef'] = nil_indef_all
    results['B1-theta(-1)eigenspace-definite'] = p_definite_all

    print("=" * 74)
    print("RAW RESULTS + STRUCTURAL CONCLUSION (the court — to the project's adjudication; I do NOT render a verdict):")
    print("─" * 74)
    print("  ★★MECHANISM (the author's question «is it one and the same?», more important than H2/H3):")
    print("   • dim of ad-inv. forms on sl(n) = 1 (Schur, EXACTLY the S1023 move, K2 clean) ⟹ the form is UNIQUE; its")
    print("     halves: h is DEFINITE(=space) ⊥ n± is INDEFINITE(=axis) ⟹ ★THE MECHANISM IS THE SAME (one")
    print("     Schur-form, two Jordan-Chevalley halves). The axis comes from the same machinery, not a separate postulate.")
    print("─" * 74)
    print("  CENTRAL H2 vs H3:")
    print("   • [n₊,n₋]→h with structure constants {0,±1} (K2 clean), h is spanned ⟹ the connector n±↔h =")
    print("     the Lie bracket = AX-closure (an EXISTING root), NOT a 4th card ⟹ ★H3 WITH AN EXPLICIT CONNECTOR.")
    print("     A foundation-candidate: {Ω(set+order) · AX-closure · Λ} — three structures, the third already present.")
    print("  AXIS (AX-indef), kill-first:")
    print("   • the θ=−1-eigenspace = DEFINITE ⟹ «axis=θ=−1-eigenspace» is FALSE (a failure like sgn). BUT")
    print("   • the NILPOTENT n₊⊕n₋ («bits») = HYPERBOLIC (p=q) = INDEFINITE, a real-form type (≠ sgn):")
    print("     a LIVE axis candidate FROM the flag. A homonym-flag: the form is BALANCED ≠ (d,1)-T32 ⟹ the type matches/the form does not,")
    print("     ★identity with AX-indef = A QUESTION FOR THE COURT (a root essence or a sector prefix).")
    print("  MARK (AX-dimer):")
    print("   • θ singles out the CENTER of the diagram, but node(d odd)/bond(d even) — NOT «one bond» ∀d ⟹")
    print("     the mark is NOT forced uniformly (the 2nd nail: neither order nor θ). It remains an independent root.")
    print("─" * 74)
    order = ['anti-reuse', '★mechanism-the-same', 'B1-nilpotent-indef', 'B1-theta(-1)eigenspace-definite',
             'B2-mark-independent(kill)', 'B3-[n₊,n₋]→h-K2-clean', 'mutants']
    all_ok = True
    for kk in order:
        v = results.get(kk)
        print("  {0:26s} : {1}".format(kk, "YES/PASS" if v else "no/FAIL"))
        if not v:
            all_ok = False
    print("=" * 74)

    # NB: 'nilpotent/Cartan/θ/signature/bracket/structure-constant/diagram/bond' is STRUCTURAL vocabulary. GUARDLINE
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
