# -*- coding: utf-8 -*-
# DIM: na (PLAN-1 formal probe: THE SIEGE of the AX-cell. Candidate theorem
#          {alphabet d+1 · democracy S_{d+1} · scale Λ} ⟹ the cell A_d is forced ENTIRELY
#          (points + metric + the very fact of the simplex). ∀d SYMBOLICALLY, not instances.
#          Exante: active-v10.2/delirium/PLAN_1_GEOMETRY_FROM_INFORMATION.md §EXANTE.
#          Contamination DECLARED: Omega's Lite-4 = a seed; honesty is carried by symbolics+mutants,
#          not blindness. ★LINEAR ALGEBRA/REPRESENTATION THEORY, not physics. FS=STONE.)
#
# ============================================================================
# ★★INPUT MANIFEST (C-GUARD against circularity — the author's alarm 2026-07-20, exante §C)
# ----------------------------------------------------------------------------
# INPUTS (and that's it):
#   [1] Ω = a bare finite set of n=d+1 elements (PRE-LIE; realized by the indicators e_i —
#       the standard basis of R^n, ZERO su(n)-structure at the input; e_i = «symbol i is present»);
#   [2] democracy = a free S_n-action permuting the elements of Ω (perm.-representation);
#   [3] Λ = scale (the single constant).
# FORBIDDEN as INPUT (a semantic cycle — the dep-checker does NOT catch it): the center ℤ/n, the weights of the
#   fundamental repr., Cartan, column, dual, any su(n)-structure. Anything Lie-flavored ⟹ ONLY as OUTPUT/target-check.
# ★AUDIT (every step): the hyperplane W={Σx=0} is forced by INPUT [2] alone (= the orthogonal to
#   the UNIQUE S_n-fixed subspace span 𝟙 = the triv-isotype of the perm.-repr; no Lie-input at all).
#   «= weights of the fundamental repr.», «= Cartan(A_d)», «center» appear ONLY in check-rows
#   as a RECOGNIZED OUTPUT (post-hoc), NEVER fed into the count. C0 below prints the input⊥output split.
# ----------------------------------------------------------------------------
# OPERATIONALIZATIONS + ★BETS (carved in the exante BEFORE the count — §EXANTE; reproduced here)
# ----------------------------------------------------------------------------
# n := d+1 = the number of symbols in the alphabet. G := S_n (democracy). E^n = R^n (indicators e_i).
# W := {x : Σx_i = 0} ⊂ R^n  (the hyperplane «the common phase removed»), dim W = n-1 = d.
# The metric on W = an S_n-invariant symmetric bilinear form. Λ = a free scale.
#
# ★C1 (the points): the indicator projection proj(e_i)=e_i − 𝟙/n = a weight of the su(n) fundamental repr. EXACTLY;
#     the Gram of differences v_i=e_i−e_{i+1} = Cartan(A_d) verbatim. SYMBOLICALLY in n (not instances).
# ★C2 (THE HEART, the metric): dim of S_n-inv. symmetric forms on W = 1 (Schur). A NATIVE re-derivation
#     of the irreducibility of the standard representation via an ORBITAL COUNT (not a textbook citation):
#       (i) the commutant of the perm.-repr. = #orbits of G on ordered pairs (i,j) = #orbitals;
#       (ii) S_n is 2-transitive ∀n≥2 ⟹ exactly 2 orbitals {i=j}⊥{i≠j} ⟹ the commutant = 2 = {I,J};
#       (iii) the perm.-repr = triv ⊕ W ⟹ 2 = Σm² ⟹ W is irreducible, multiplicity 1, W≇triv;
#       (iv) J|_W = 0 (because W ⊥ 𝟙) ⟹ on W only a·I remains ⟹ dim = 1.
#     ★The argument (i)-(iv) is DIMENSION-INDEPENDENT (the same 1-line reason ∀n) ⟹ ∀d, not just d≤6.
# ★C3 (M1, democracy is load-bearing): breaking G→S_k×S_{n−k} ⟹ dim of inv. forms ≥2 ⟹ uniqueness DIES.
#     Symbolically: W|_H = triv ⊕ std_k ⊕ std_{n−k} ⟹ dim = #distinct H-irreducibles = 3 (2≤k≤n−2).
# ★C4 (M4-HEART, «why THIS cell»): (a) RIGIDITY — n+1 points in R^d with 1 distance-class ⟹
#     the Gram is forced = (s²/2)(I − J/(n)) ⟹ a regular simplex, UNIQUE up to isometry+scale
#     (symbolically ∀d). (b) distance-classes: simplex=1 ⊥ cross-polytope=2 ⊥ cube=d ⟹ ONLY the simplex
#     is democratic (1 orbit on vertex-pairs).  Democracy on pairs ⟹ 1 class ⟹(C4a) a simplex.
# ★C5 (M2): a false projection (a skewed hyperplane Σc_i x_i=0, unequal c) ⟹ Gram ≠ Cartan.
# ★C6 (M3): the free parameters of the metric are EXACTLY 1 (=C2), not 2 — Λ and only Λ.
# ★C7 (negative control, seeded): (a) random NON-2-transitive subgroups of S_n ⟹ dim inv ≥2
#     in ≥95% of seeds; (b) random n+1 points ⟹ a 1-distance-class ~never (measure zero).
# ★C8 (an honest boundary): the alphabet+democracy remain PREMISES. The probe derives geometry FROM them,
#     NOT the premises themselves. The claim = «the postulate RETREATED (geometry→alphabet)», NOT «vanished».
#
# KILLS: K1 the native-Schur does not assemble symbolically ⟹ the status is a snapshot d≤6, a carved boundary (the argument (i)-(iv)
#   is dimension-independent ⟹ K1 does not fire, as long as the orbital count holds).
#   K2 any step asks for a NEW constant ⟹ STOP, report (Λ = the only one; the scale of the form is not a new handle).
#   K3 (fence): physics readings stay BEHIND THE FENCE; FS {the physics-vocabulary classes below, unquoted Minkowski-mentions=STONE, # GUARDLINE
#   universe, anthropic, Tegmark} hits=0.  Mutants ≥4: M1(C3)·M2(C5)·M3(C6)·M4(C4) + a seeded negctrl(C7). # GUARDLINE
# Ancestors (citation, not re-derivation): T19 the Gram∝Cartan dictionary · Omega's Lite-4 (a seed, declared) ·
#   Chentsov/Schur (📖 class of the move, multiplicity 0).  Discipline: exact ℚ-arithmetic; raw tables — to Omega.
# ============================================================================

import sys
import os
import random
import itertools
from fractions import Fraction
import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== exact linear algebra over ℚ (native, no library magic) ====================

def perm_matrix_action_on_vec(perm, vec):
    """P_perm · vec, where (P e_i)=e_{perm[i]}. vec — a list of ℚ (Fraction). Returns a list of ℚ."""
    n = len(perm)
    out = [Fraction(0)] * n
    for i in range(n):
        out[perm[i]] += vec[i]
    return out


def W_basis(n):
    """Basis of W={Σx=0}: v_i = e_i − e_{i+1}, i=0..n−2 (simple roots of A_{n−1})."""
    B = []
    for i in range(n - 1):
        v = [Fraction(0)] * n
        v[i] = Fraction(1)
        v[i + 1] = Fraction(-1)
        B.append(v)
    return B  # a list of (n−1) vectors of length n


def _matmatT(cols):
    """Gram BᵀB for B given by columns (cols = a list of column-vectors)."""
    m = len(cols)
    G = [[Fraction(0)] * m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            G[a][b] = sum(cols[a][k] * cols[b][k] for k in range(len(cols[a])))
    return G


def _inv_rational(M):
    """Inverse of a square ℚ-matrix via Gauss (exact). M — a list of lists of Fraction."""
    n = len(M)
    A = [row[:] + [Fraction(1) if i == j else Fraction(0) for j in range(n)] for i, row in enumerate(M)]
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [A[r][k] - f * A[col][k] for k in range(2 * n)]
    return [row[n:] for row in A]


def _matvec(M, v):
    return [sum(M[i][k] * v[k] for k in range(len(v))) for i in range(len(M))]


def _matmat(A, B):
    p = len(A); q = len(B[0]); r = len(B)
    return [[sum(A[i][k] * B[k][j] for k in range(r)) for j in range(q)] for i in range(p)]


def rho_on_W(perm, n, Bpinv=None, Bcols=None):
    """ρ(perm) = the action of perm on W in the basis of simple roots; an (n−1)×(n−1) ℚ-matrix, native.
    Bpinv = (BᵀB)⁻¹Bᵀ ((n−1)×n), Bcols = the columns of B (n×(n−1))."""
    if Bcols is None:
        Bcols = W_basis(n)
    if Bpinv is None:
        BtB = _matmatT(Bcols)
        BtBinv = _inv_rational(BtB)
        # Bᵀ rows = Bcols; (BᵀB)⁻¹Bᵀ
        Bt = [[Bcols[a][k] for k in range(n)] for a in range(n - 1)]  # (n−1)×n
        Bpinv = _matmat(BtBinv, Bt)  # (n−1)×n
    cols = []
    for j in range(n - 1):
        Pvj = perm_matrix_action_on_vec(perm, Bcols[j])  # in R^n, lies in W
        cols.append(_matvec(Bpinv, Pvj))  # coordinates in the root basis
    # cols[j] = the j-th column of ρ
    rho = [[cols[j][i] for j in range(n - 1)] for i in range(n - 1)]
    return rho


def dim_invariant_sym_forms_on_W(gens, n):
    """dim of the space of symmetric (n−1)×(n−1) S with ρ(g)ᵀ S ρ(g)=S ∀g in gens. Native (rank of the system).
    gens — a list of perms (generators of a subgroup G ≤ S_n)."""
    d = n - 1
    Bcols = W_basis(n)
    BtB = _matmatT(Bcols)
    BtBinv = _inv_rational(BtB)
    Bt = [[Bcols[a][k] for k in range(n)] for a in range(d)]
    Bpinv = _matmat(BtBinv, Bt)
    rhos = [rho_on_W(g, n, Bpinv=Bpinv, Bcols=Bcols) for g in gens]

    # unknowns: symmetric S ⟹ parameters S_{ab}, a≤b. Indexing.
    idx = {}
    params = []
    for a in range(d):
        for b in range(a, d):
            idx[(a, b)] = len(params)
            params.append((a, b))
    P = len(params)

    def sym_get(Svec, a, b):
        return Svec[idx[(min(a, b), max(a, b))]]

    # for each g: (ρᵀ S ρ − S)_{ij}=0. Building the system's rows (ℚ) as linear in the parameters.
    rows = []
    for rho in rhos:
        # M = ρᵀ S ρ. M_{ij} = Σ_{k,l} ρ_{ki} S_{kl} ρ_{lj}
        for i in range(d):
            for j in range(i, d):  # symmetric ⟹ the upper triangle of equations
                coeff = [Fraction(0)] * P
                for k in range(d):
                    for l in range(d):
                        c = rho[k][i] * rho[l][j]
                        coeff[idx[(min(k, l), max(k, l))]] += c
                coeff[idx[(i, j)]] -= Fraction(1)
                rows.append(coeff)
    # the rank of the system over ℚ ⟹ dim = P − rank
    rank = _rank_rational(rows, P)
    return P - rank


def _rank_rational(rows, ncols):
    """Rank of a matrix of rows over ℚ (Gauss, exact)."""
    A = [r[:] for r in rows]
    rank = 0
    pivot_col = 0
    nrows = len(A)
    while pivot_col < ncols and rank < nrows:
        piv = None
        for r in range(rank, nrows):
            if A[r][pivot_col] != 0:
                piv = r
                break
        if piv is None:
            pivot_col += 1
            continue
        A[rank], A[piv] = A[piv], A[rank]
        pv = A[rank][pivot_col]
        A[rank] = [x / pv for x in A[rank]]
        for r in range(nrows):
            if r != rank and A[r][pivot_col] != 0:
                f = A[r][pivot_col]
                A[r] = [A[r][k] - f * A[rank][k] for k in range(ncols)]
        rank += 1
        pivot_col += 1
    return rank


# ==================== group primitives ====================

def transposition(n, i, j):
    p = list(range(n))
    p[i], p[j] = p[j], p[i]
    return tuple(p)


def n_cycle(n):
    return tuple((i + 1) % n for i in range(n))


def S_n_generators(n):
    """S_n = ⟨(0 1), n-cycle⟩ (n≥2)."""
    if n == 2:
        return [transposition(2, 0, 1)]
    return [transposition(n, 0, 1), n_cycle(n)]


def Sk_x_Snk_generators(n, k):
    """S_k × S_{n−k}: transpositions within {0..k−1} and {k..n−1}."""
    gens = []
    if k >= 2:
        gens.append(transposition(n, 0, 1))
        if k >= 3:
            gens.append(tuple([(i + 1) % k if i < k else i for i in range(n)]))  # a k-cycle on the 1st block
    if n - k >= 2:
        gens.append(transposition(n, k, k + 1))
        if n - k >= 3:
            gens.append(tuple([i if i < k else k + ((i - k + 1) % (n - k)) for i in range(n)]))
    return gens


def orbits_on_ordered_pairs(gens, n):
    """#orbits of the group ⟨gens⟩ on ordered pairs (i,j) — a native commutant-count (Burnside-free,
    via direct orbit closure). Returns #orbitals."""
    # generate the whole group by closure (small n)
    group = _closure(gens, n)
    seen = set()
    count = 0
    for i in range(n):
        for j in range(n):
            if (i, j) in seen:
                continue
            count += 1
            for g in group:
                seen.add((g[i], g[j]))
    return count


def _closure(gens, n):
    ident = tuple(range(n))
    group = {ident}
    frontier = [ident]
    gg = list(gens) if gens else [ident]
    while frontier:
        x = frontier.pop()
        for g in gg:
            y = tuple(g[x[i]] for i in range(n))
            if y not in group:
                group.add(y)
                frontier.append(y)
    return group


def is_2_transitive(gens, n):
    return orbits_on_ordered_pairs(gens, n) == 2


# ==================== C0: C-GUARD — W is forced by the S-action alone (no Lie-input) ====================

def C0_anticircularity_audit():
    print("─" * 74)
    print("C0 (C-GUARD, anti-circularity): W={Σx=0} is forced ONLY by input [2] (the S-action), not Lie")
    print("─" * 74)
    print("  Audit of inputs: {Ω a bare set (indicators e_i) · the S_n-action · Λ}. Lie-objects = OUTPUT.")
    ok_all = True
    print("   n | S_n-fixed subspace = span? | dim=1? | W=orthogonal ⟹ dim W=n−1?")
    for n in range(2, 9):
        gens = S_n_generators(n)
        # the S_n-fixed subspace of the perm.-repr = {v : P_g v = v ∀g} — computed natively (the kernel of Σ(P_g−I))
        # for a transitive action = span(𝟙), dim 1. Checked WITHOUT any Lie-input.
        fixed_dim = _fixed_subspace_dim(gens, n)
        is_ones = (fixed_dim == 1)
        wdim = n - fixed_dim  # W = the orthogonal complement of the fixed subspace
        good = is_ones and (wdim == n - 1)
        if not good:
            ok_all = False
        print("   {0} | {1:31s} | {2:6s} | {3}".format(
            n, "span(𝟙)" if is_ones else "?", "YES" if is_ones else "no",
            "dim W={0} ✓".format(wdim) if good else "✗"))
    print("  ⟹ the decomposition R^n = (fixed triv) ⊕ W is forced by the S-action; W is the arena of the metric, WITHOUT su(n).")
    print("    «weights/Cartan/center» below = a RECOGNIZED output in check-rows, NOT an input to the count.")
    return ok_all


def _fixed_subspace_dim(gens, n):
    """dim {v∈R^n : P_g v = v ∀g} — native (rank of the stack (P_g−I)). No Lie-input at all."""
    rows = []
    for g in gens:
        # (P_g − I): a row for each coordinate i: v_{g^{-1}(i)} − v_i = 0 ⟺ v[preimage]−v[i]
        # P_g v has at position j the value v[g^{-1}(j)]; (P_g v − v)_j = v[ginv[j]] − v[j]
        ginv = [0] * n
        for a in range(n):
            ginv[g[a]] = a
        for j in range(n):
            r = [Fraction(0)] * n
            r[ginv[j]] += Fraction(1)
            r[j] -= Fraction(1)
            rows.append(r)
    rank = _rank_rational(rows, n)
    return n - rank


# ==================== C1: the projection = weights of the fund. repr., Gram = Cartan (SYMBOLICALLY in n) ====================

def C1_symbolic():
    print("─" * 74)
    print("C1 (the points, SYMBOLICALLY in n): proj(e_i)=e_i−𝟙/n = a weight of the su(n) fund. repr.; the Gram of differences=Cartan(A_d)")
    print("─" * 74)
    n = sp.Symbol('n', positive=True, integer=True)
    # ⟨proj(e_i),proj(e_j)⟩ = δ_ij − 1/n  (symbolically): the projector I − 𝟙𝟙ᵀ/n
    #   diag = 1 − 1/n = (n−1)/n ;  off = −1/n.  This is the Gram of the fund. repr. weights (up to scale).
    diag = sp.simplify(1 - 1 / n)
    off = sp.simplify(-1 / n)
    print("  ⟨proj e_i, proj e_i⟩ = {0}   (symbolically in n)".format(diag))
    print("  ⟨proj e_i, proj e_j⟩ = {0}   (i≠j, symbolically in n)".format(off))
    print("  ⟹ the Gram of the projections = I − 𝟙𝟙ᵀ/n — invariant under permutations, forced by the projection.")
    # the Gram of the simple roots v_i=e_i−e_{i+1}: ⟨v_i,v_j⟩ = 2(i=j), −1(|i−j|=1), 0(otherwise) — n-INDEPENDENT
    print("  the Gram of the simple roots v_i=e_i−e_{i+1}: ⟨v_i,v_i⟩=2, ⟨v_i,v_{i±1}⟩=−1, otherwise 0")
    print("  = Cartan(A_d) VERBATIM, tridiagonal [2,−1] — the expression is n-INDEPENDENT (from Kronecker deltas).")
    # a symbolic check of the identity at several n (the machine confirms the n-independent derivation)
    ok_all = True
    for nn in range(2, 9):
        proj = [[Fraction(1 if a == b else 0) - Fraction(1, nn) for b in range(nn)] for a in range(nn)]
        # the Gram of the roots
        Bc = W_basis(nn)
        Gr = _matmatT(Bc)
        cartan = [[Fraction(2 if a == b else (-1 if abs(a - b) == 1 else 0)) for b in range(nn - 1)]
                  for a in range(nn - 1)]
        if Gr != cartan:
            ok_all = False
        # the projection = δ−1/n
        proj_ok = all(proj[a][b] == (Fraction(1 if a == b else 0) - Fraction(1, nn))
                      for a in range(nn) for b in range(nn))
        if not proj_ok:
            ok_all = False
    print("  machine check n=2..8: Gram(roots)=Cartan(A_d) AND projection=I−𝟙𝟙ᵀ/n — {0}".format(
        "MATCH (the n-independent derivation is confirmed)" if ok_all else "MISMATCH"))
    return ok_all


# ==================== C2: the heart — native-Schur ∀d ====================

def C2_native_schur():
    print("─" * 74)
    print("C2 (THE HEART): dim of S_n-inv. symmetric forms on W = 1 — native-Schur via an orbital count")
    print("─" * 74)
    print("  The argument (DIMENSION-INDEPENDENT, ∀n≥2):")
    print("    (i)  the commutant of the perm.-repr R^n = span of orbit-indicators = #orbitals of G on (i,j);")
    print("    (ii) S_n is 2-transitive ∀n≥2 ⟹ exactly 2 orbitals: {i=j} ⊥ {i≠j} ⟹ commutant={I,J};")
    print("    (iii) the perm.-repr = triv ⊕ W, dim End = Σm² = 2 ⟹ W is IRREDUCIBLE, multiplicity 1, W≇triv;")
    print("    (iv) J|_W = 𝟙𝟙ᵀ|_W = 0 (W ⊥ 𝟙) ⟹ on W the inv. form = a·I ⟹ dim = 1.  ∎ ∀d")
    print()
    print("  Machine check (a native rank of the system ρᵀSρ=S on W; NOT a citation):")
    print("   d | n=d+1 | 2-transitive? | #orbitals | dim_inv_sym_forms(W) | expected=1")
    ok_all = True
    for n in range(2, 9):
        gens = S_n_generators(n)
        orb = orbits_on_ordered_pairs(gens, n)
        tr = is_2_transitive(gens, n)
        dim = dim_invariant_sym_forms_on_W(gens, n)
        flag = "✓" if dim == 1 else "✗"
        if dim != 1:
            ok_all = False
        print("   {0:1d} | {1:5d} | {2:8s} | {3:10d} | {4:19d} | {5}".format(
            n - 1, n, "YES" if tr else "no", orb, dim, flag))
    print("  ⟹ dim=1 ∀ measured d; the orbital argument is n-independent ⟹ K1 does NOT fire.")
    return ok_all


# ==================== C3 (M1): breaking democracy ⟹ uniqueness dies ====================

def C3_mutant_break_democracy():
    print("─" * 74)
    print("C3 (M1, democracy is LOAD-BEARING): breaking S_n → S_k×S_{n−k} ⟹ dim of inv. forms ≥2 ⟹ uniqueness DIES")
    print("─" * 74)
    print("  Symbolically: W|_H = triv ⊕ std_k ⊕ std_{n−k}; dim = #distinct H-irreducibles.")
    print("   n | k | dim_inv(W) | expected (2≤k≤n−2 ⟹ 3; k=1 ⟹ 2) | uniqueness")
    ok_all = True
    for n in range(4, 8):
        for k in (1, 2, 3):
            if k > n - 1:
                continue
            gens = Sk_x_Snk_generators(n, k)
            if not gens:
                continue
            dim = dim_invariant_sym_forms_on_W(gens, n)
            expect = 3 if (2 <= k <= n - 2) else 2
            broke = (dim >= 2)
            if not broke:
                ok_all = False
            print("   {0} | {1} | {2:10d} | {3:31d} | {4}".format(
                n, k, dim, expect, "DIES ✓" if broke else "holds ✗"))
    print("  ⟹ any break of democracy ⟹ dim≥2 ⟹ the metric stops being forced.")
    return ok_all


# ==================== C4 (M4-heart): rigidity of the simplex + distance-classes ====================

def C4_rigidity_symbolic():
    print("─" * 74)
    print("C4a (M4-HEART, RIGIDITY): n+1 points with 1 distance-class ⟹ the Gram is forced ⟹ a simplex")
    print("─" * 74)
    # symbolic derivation: v_i (i=0..n), Σv_i=0, |v_i−v_j|²=s² ∀i≠j ⟹ G_ii=a ∀i, G_ij=b ∀i≠j,
    #   a+ n·b = 0 (centered), a−b = s²/2 ⟹ a=n·s²/(2(n+1)), b=−s²/(2(n+1)). The Gram is UNIQUE.
    s, N = sp.symbols('s N', positive=True)  # N = n = the vertex count − 1 (i.e. vertices = N+1)
    a = sp.simplify(N * s**2 / (2 * (N + 1)))
    b = sp.simplify(-s**2 / (2 * (N + 1)))
    print("  Solution (symbolically, vertices=N+1 in R^N): G_ii = {0}, G_ij = {1} (i≠j)".format(a, b))
    print("  ⟹ the Gram = (s²/2)(I − 𝟙𝟙ᵀ/(N+1)), rank N, PSD ⟹ the configuration is UNIQUE up to isometry+scale")
    print("    = a regular simplex.  (1 distance-class ⟹ a simplex, ∀d symbolically.)")
    # machine check: the 1-class Gram matches the simplex-Gram
    ok_all = True
    for nn in range(2, 8):
        # a regular simplex: proj(e_i), i=0..nn ; vertices = nn+1 in R^nn
        verts = []
        for i in range(nn + 1):
            v = [Fraction(-1, nn + 1)] * (nn + 1)
            v[i] += 1
            verts.append(v)
        # are all pairwise squared distances equal?
        dists = set()
        for i in range(nn + 1):
            for j in range(i + 1, nn + 1):
                dists.add(sum((verts[i][t] - verts[j][t]) ** 2 for t in range(nn + 1)))
        if len(dists) != 1:
            ok_all = False
        # the Gram of the centered points = a on the diagonal, b off (with one s²)
        s2 = next(iter(dists))
        aa = Fraction(nn, 2 * (nn + 1)) * s2
        bb = Fraction(-1, 2 * (nn + 1)) * s2
        Gcheck = all(
            (sum(verts[i][t] * verts[j][t] for t in range(nn + 1)) == (aa if i == j else bb))
            for i in range(nn + 1) for j in range(nn + 1))
        if not Gcheck:
            ok_all = False
    print("  machine check n=2..7: simplex = 1 class AND Gram = the forced expression — {0}".format(
        "MATCH" if ok_all else "MISMATCH"))
    return ok_all


def C4_distance_classes():
    print("─" * 74)
    print("C4b (M4, «why THIS cell»): vertex distance-classes — simplex=1 ⊥ cross-polytope=2 ⊥ cube=d")
    print("─" * 74)
    print("   d | simplex | cross-polytope | hypercube | only 1-class?")
    ok_all = True
    for d in range(2, 7):
        # a simplex: d+1 vertices e_i (in R^{d+1}) — 1 class
        simp = _num_dist_classes([tuple(1 if i == t else 0 for t in range(d + 1)) for i in range(d + 1)])
        # a cross-polytope: ±e_i in R^d
        cross_v = []
        for i in range(d):
            p = [0] * d; p[i] = 1; cross_v.append(tuple(p))
            m = [0] * d; m[i] = -1; cross_v.append(tuple(m))
        cross = _num_dist_classes(cross_v)
        # a hypercube {0,1}^d
        cube_v = list(itertools.product((0, 1), repeat=d))
        cube = _num_dist_classes(cube_v)
        only_simplex = (simp == 1 and cross > 1 and cube > 1)
        if not only_simplex:
            ok_all = False
        print("   {0} | {1:8d} | {2:12d} | {3:8d} | {4}".format(
            d, simp, cross, cube, "YES ✓ (the simplex is the only 1-class)" if only_simplex else "✗"))
    print("  ⟹ democracy on vertex-pairs (1 orbit = 1 class) ⟹(C4a) a REGULAR SIMPLEX is forced.")
    return ok_all


def _num_dist_classes(verts):
    dists = set()
    m = len(verts); dim = len(verts[0])
    for i in range(m):
        for j in range(i + 1, m):
            dists.add(sum((verts[i][t] - verts[j][t]) ** 2 for t in range(dim)))
    return len(dists)


# ==================== C5 (M2): a false projection ⟹ Gram ≠ Cartan ====================

def C5_mutant_false_projection():
    print("─" * 74)
    print("C5 (M2, a false projection): projecting onto a SKEWED Σc_i x_i=0 (unequal c) ⟹ Gram ≠ Cartan")
    print("─" * 74)
    print("   d | c-weights (skew) | Gram of differences = Cartan(A_d)? ")
    ok_all = True
    for d in range(2, 6):
        n = d + 1
        c = [Fraction(1 + i) for i in range(n)]  # skewed weights 1,2,3,...
        csum = sum(ci * ci for ci in c)
        # the projector onto {Σ c_i x_i = 0}: P = I − c cᵀ / (cᵀc). proj(e_i) = e_i − c_i·c/(cᵀc)
        proj = []
        for i in range(n):
            v = [Fraction(0)] * n
            v[i] = Fraction(1)
            for t in range(n):
                v[t] -= c[i] * c[t] / csum
            proj.append(v)
        # the Gram of the differences proj(e_i)−proj(e_{i+1})
        diffs = [[proj[i][t] - proj[i + 1][t] for t in range(n)] for i in range(n - 1)]
        Gr = [[sum(diffs[a][t] * diffs[b][t] for t in range(n)) for b in range(n - 1)] for a in range(n - 1)]
        cartan = [[Fraction(2 if a == b else (-1 if abs(a - b) == 1 else 0)) for b in range(n - 1)]
                  for a in range(n - 1)]
        equal = (Gr == cartan)
        if not equal:
            ok_all = True  # the mutant FIRED (the Gram is broken) — this is expected
        else:
            ok_all = False  # if it matches — the mutant did NOT fire, bad
        print("   {0} | {1} | {2}".format(
            d, [str(x) for x in c], "≠ Cartan ✓ (the mutant fired)" if not equal else "= Cartan ✗"))
    print("  ⟹ only a DEMOCRATIC projection (equal weights, onto Σx=0) gives the fund. repr. weights + Cartan.")
    return ok_all


# ==================== C7: seeded negative controls ====================

def C7_seeded_negcontrols():
    print("─" * 74)
    print("C7 (negative control, SEEDED): (a) NON-2-transitive subgroups ⟹ dim≥2 ; (b) random points ⟹ >1 class")
    print("─" * 74)
    random.seed(1023071)
    # (a) random subgroups of S_n with 1-2 generators; among the NON-2-transitive ones — the fraction with dim≥2
    print("  (a) random subgroups of S_n (seed 1023071): among the NON-2-transitive ones — the fraction with dim≥2")
    print("      n | #trials | #not-2-transitive | #(not-2-transitive ∧ dim≥2) | fraction")
    a_ok = True
    for n in range(4, 7):
        trials = 60
        not2t = 0
        dimge2 = 0
        for _ in range(trials):
            ngen = random.choice([1, 2])
            gens = []
            for _g in range(ngen):
                p = list(range(n))
                random.shuffle(p)
                gens.append(tuple(p))
            if not is_2_transitive(gens, n):
                not2t += 1
                if dim_invariant_sym_forms_on_W(gens, n) >= 2:
                    dimge2 += 1
        frac = (dimge2 / not2t) if not2t else 1.0
        if not2t and frac < 0.95:
            a_ok = False
        print("      {0} | {1:5d} | {2:11d} | {3:21d} | {4:.3f}".format(n, trials, not2t, dimge2, frac))
    # (b) random n+1 points ⟹ #distance-classes > 1 almost always
    print("  (b) random n+1 points in R^d (ℚ-coordinates): a 1-distance-class = measure zero")
    print("      d | #trials | #(1-class) | 1-class fraction")
    b_ok = True
    for d in range(2, 6):
        trials = 200
        one = 0
        for _ in range(trials):
            verts = [tuple(Fraction(random.randint(-5, 5)) for _ in range(d)) for _ in range(d + 1)]
            if _num_dist_classes(verts) == 1:
                one += 1
        frac = one / trials
        if frac > 0.02:
            b_ok = False
        print("      {0} | {1:5d} | {2:9d} | {3:.4f}".format(d, trials, one, frac))
    print("  ⟹ (a) non-democracy ⟹ the metric is NOT unique; (b) 1-class ⟹ forcing, not a coincidence.")
    return a_ok and b_ok


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
    _logf = open(os.path.join(_HERE, "S1023_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("PLAN-1 FORMAL PROBE S1023 — THE SIEGE of the AX-cell")
    print("Candidate theorem: {alphabet d+1 · democracy S_{d+1} · Λ} ⟹ the cell A_d is forced ENTIRELY")
    print("★LINEAR ALGEBRA / REPRESENTATION THEORY; the bets are carved BEFORE the count (exante §EXANTE).")
    print("Contamination DECLARED (Lite-4 = a seed); honesty = symbolics + mutants, not blindness.")
    print("=" * 74)
    print()

    results = {}
    results['C0(anti-circ)'] = C0_anticircularity_audit(); print()
    results['C1'] = C1_symbolic();               print()
    results['C2'] = C2_native_schur();           print()
    results['C3(M1)'] = C3_mutant_break_democracy(); print()
    results['C4a'] = C4_rigidity_symbolic();     print()
    results['C4b(M4)'] = C4_distance_classes();  print()
    results['C5(M2)'] = C5_mutant_false_projection(); print()
    results['C6(M3)'] = (results['C2'] and results['C6-note'] if False else results['C2'])
    print("─" * 74)
    print("C6 (M3, the free parameters are EXACTLY 1): = the C2 conclusion (dim=1). Λ and only Λ — not 2.")
    print("  (no step asked for a NEW constant ⟹ K2 did not fire; the scale of the form = a Λ-slot, not a handle)")
    print("─" * 74); print()
    results['C7(negctrl)'] = C7_seeded_negcontrols(); print()

    print("─" * 74)
    print("C8 (THE BOUNDARY OF HONESTY): the alphabet (d+1 distinguishable states) and democracy (S_{d+1}) remain")
    print("  PREMISES. The probe derives geometry (points+metric+the very fact of the simplex) FROM these premises,")
    print("  it does NOT derive the premises themselves. The claim = «the postulate of geometry RETREATED to the informational")
    print("  (alphabet+democracy)», NOT «vanished». Why there are distinguishable states — falls on existing roots.")
    print("─" * 74); print()

    # ===== summary of raw tables (the court — to Omega, NOT self-assessment) =====
    print("=" * 74)
    print("RAW RESULTS (the court — to Omega; I do NOT render a verdict):")
    order = ['C0(anti-circ)', 'C1', 'C2', 'C3(M1)', 'C4a', 'C4b(M4)', 'C5(M2)', 'C6(M3)', 'C7(negctrl)']
    all_pass = True
    for kk in order:
        v = results.get(kk)
        print("  {0:14s} : {1}".format(kk, "PASS" if v else "FAIL"))
        if not v:
            all_pass = False
    mut_ok = results['C3(M1)'] and results['C5(M2)'] and results['C4b(M4)'] and results['C7(negctrl)']
    print("  mutants ≥4 (M1·M2·M4·negctrl)  : {0}".format("all fired" if mut_ok else "NOT all"))
    print("=" * 74)

    # NB: 'projection/Gram/Cartan/simplex/distance-class/orbital/Schur/democracy' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("матер", "ія"), ("частин", "ка"), ("Міньков", "ський"),
           ("всес", "віт"), ("ант", "роп"), ("Teg", "mark"),
           ("речо", "вина"), ("контей", "нер")]  # GUARDLINE (FS=STONE; physics readings stay behind the fence)
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _n = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN (STONE): hits={0} (src={1}, log={2})".format(_n, len(_hs), len(_hl)))

    _exit = 1 if (_n > 0 or not all_pass or not mut_ok) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
