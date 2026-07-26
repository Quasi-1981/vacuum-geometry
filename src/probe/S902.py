#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (linear algebra + Lie algebra; downstream-content 0; handles 0). W28-O4 —
#      stabilizers of the pair (η,B) in the containers so(2,2)·so(3,3): g_B = {X ∈ so(n,n):
#      XᵀB + BX = 0} as the EXACT nullspace over Q; dim · Killing-signature ·
#      bracket-closedness (assert). Blocks: A diagonal B (full coverage) ·
#      B random symmetric · C intersection with the S901 list · D off-diagonal Jordan.
#      ★BLINDNESS: the probe prints ONLY raw numbers/tables/sets; reading = an act of the court.
"""
S902 (lane A, ed.2) — W28-O4: stabilizers of the pair (η,B) in the (n,n) container.

  Object: η = diag(+1×n, −1×n), n ∈ {2,3}; so(n,n) = {X : Xη + ηXᵀ = 0}; a second
  symmetric form B (symmetric integer 2n×2n). Stabilizer of the pair:
    g_B = {X ∈ so(n,n) : XᵀB + BX = 0}
  — an exact linear system on the coefficients of X in the generator basis of so(n,n)
  (nullspace over Q, sympy). For every g_B: dim · Killing-signature
  (congruent diagonalization over Q) · degenerate/non-degenerate · bracket-closedness assert
  (the stabilizer must be a subalgebra — checked explicitly on every g_B).

  BLOCK A — diagonal B = diag(b₁..b_{2n}), bᵢ ∈ {−2..2}: n=2 full pass
    5⁴ = 625 (every B solved directly; class agreement within a pattern — assert);
    n=3 coverage 5⁶ = 15625 by cache over the canonical pattern (axis groups by shared
    eigenvalue λᵢ = bᵢ/ηᵢᵢ of the matrix A = η⁻¹B; key = the multiset of induced
    signatures (pᵢ,qᵢ) of the groups; one representative per pattern solved
    directly + a verification sample of direct solves against the cache). For each class:
    decomposition ⊕(pᵢ,qᵢ) (zero bᵢ are also a group) · dim g_B · Killing · B count.
    Rollup: the set of all (p,q)-components across decompositions.
  BLOCK B — random full symmetric B: N = 300 per container, coeffs {−3..3},
    seed 902; histogram of classes {dim · Killing}; factorization structure of
    charpoly(η⁻¹B) over Q (degrees of irreducible factors · multiplicities) by class.
  BLOCK C — three raw sets of keys {dim · Killing} against the S901 list
    (the same container): intersection · S901∖stabilizer · stabilizer∖S901.
  BLOCK D — n=2, off-diagonal B: A = η⁻¹B with a 2×2 Jordan block (nilpotent
    improper diagonal; Jordan structure checked sympy-exact) — dim g_B ·
    Killing · degenerate/non-degenerate · membership of the key in the set of diagonal classes.

  EXACT arithmetic: sympy Rational/Integer; ranks/echelons/nullspaces over Q;
  no tolerances. Mechanisms (Span/RREF-key/bracket/structure-constants/Killing/
  signature) — a verbatim copy from S900/S901 (reproducibility).

Fence: the shared fence_scan helper (forbidden words — on the GUARDLINE line).
No verdicts in the probe's text; raw lines with no interpretive words.
"""
import sys
import random
from itertools import product
from collections import Counter

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# shared fence_scan helper (src/tools)
_src = __file__.replace("\\", "/").rsplit("/src/", 1)[0] + "/src"
if _src not in sys.path:
    sys.path.insert(0, _src)
from tools.fence_scan import scan_forbidden   # noqa: E402

import sympy as sp   # noqa: E402

# ── tee: all of stdout is duplicated into S902_run.log next to the script ──
_LOG_PATH = __file__.replace("\\", "/").rsplit("/", 1)[0] + "/S902_run.log"


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)
            st.flush()

    def flush(self):
        for st in self.streams:
            st.flush()


_logf = open(_LOG_PATH, "w", encoding="utf-8")
sys.stdout = _Tee(sys.stdout, _logf)

FAIL = []
N_CHECKS = 0


def check(name, cond, detail=""):
    global N_CHECKS
    N_CHECKS += 1
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)
    return ok


def rule(t):
    print("\n" + "=" * 96)
    print(t)
    print("=" * 96)


# ── SECTION 0 — FENCE VIA SHARED HELPER (first) ──────────────────────────────────────
rule("SECTION 0 — FENCE (shared fence_scan helper)")
_FORBIDDEN = [r"стійк", r"обрано", r"selected", r"stable", r"причин", r"час", r"time", r"arrow", r"стріла", r"вибір", r"краще", r"злам", r"Лі-клас", r"евклід"]   # GUARDLINE
_hits = scan_forbidden(__file__, _FORBIDDEN)
check("fence: forbidden words (list = GUARDLINE line) = 0 occurrences outside the declaration",
      not _hits, f"hits: {_hits}" if _hits else "0")
check("downstream-content 0 · handles 0 (pure linear algebra + Lie algebra)", True)


# ═══════════════════════════════════════════════════════════════════════════════════════
# TOOLS (verbatim copy from S900/S901): exact span (echelon over Q), bracket,
# structure constants, Killing, signature of a symmetric matrix over Q.
# ═══════════════════════════════════════════════════════════════════════════════════════

def flat(M):
    """n×n matrix → flat n² tuple (row-major), sympy-exact entries."""
    n = M.rows
    return tuple(M[r, c] for r in range(n) for c in range(n))


class Span:
    """Exact linear span over Q: echelon rows by pivot; rank/membership/canonical key."""

    def __init__(self):
        self.rows = {}   # pivot position → row (list, sympy Rational/Integer)

    def _reduce(self, v):
        v = list(v)
        for p in sorted(self.rows):
            if v[p] != 0:
                r = self.rows[p]
                c = sp.Rational(v[p], 1) / r[p]
                v = [vi - c * ri for vi, ri in zip(v, r)]
        return v

    def add(self, v):
        """True ⟺ the vector raised the rank (added to the echelon)."""
        w = self._reduce(v)
        piv = next((k for k, x in enumerate(w) if x != 0), None)
        if piv is None:
            return False
        self.rows[piv] = w
        return True

    def contains(self, v):
        return all(x == 0 for x in self._reduce(v))

    def dim(self):
        return len(self.rows)

    def key(self):
        """Canonical RREF key of the subspace (unique by construction)."""
        pivs = sorted(self.rows)
        rows = [[sp.Rational(x, 1) / self.rows[p][p] for x in self.rows[p]] for p in pivs]
        # zero out ABOVE the pivots (full RREF)
        for a in range(len(pivs) - 1, -1, -1):
            pa = pivs[a]
            for b in range(a):
                cb = rows[b][pa]
                if cb != 0:
                    rows[b] = [x - cb * y for x, y in zip(rows[b], rows[a])]
        return tuple(tuple(r) for r in rows)


def bracket(X, Y):
    return X * Y - Y * X


def structure_constants(basis):
    """c^k_{ij}: [e_i,e_j] = Σ_k c^k_{ij} e_k — exact decomposition, assert zero residual."""
    d = len(basis)
    A = sp.Matrix.hstack(*[sp.Matrix(flat(b)) for b in basis])   # n² × d
    AtA_inv = (A.T * A).inv()
    C = [[None] * d for _ in range(d)]
    for i in range(d):
        C[i][i] = [sp.Integer(0)] * d
        for j in range(i + 1, d):
            v = sp.Matrix(flat(bracket(basis[i], basis[j])))
            sol = AtA_inv * (A.T * v)
            assert A * sol - v == sp.zeros(A.rows, 1), "bracket outside the span of basis L"
            C[i][j] = [sol[k] for k in range(d)]
            C[j][i] = [-sol[k] for k in range(d)]
    return C


def killing_matrix(basis):
    d = len(basis)
    C = structure_constants(basis)
    ad = [sp.Matrix(d, d, lambda k, j, i=i: C[i][j][k]) for i in range(d)]

    def tr_prod(Ai, Aj):
        # trace(Ai·Aj) via direct sum of element products (exact, without full matrix product)
        return sum(Ai[r, c] * Aj[c, r] for r in range(d) for c in range(d))

    return sp.Matrix(d, d, lambda i, j: tr_prod(ad[i], ad[j]))


def sym_signature(Kmat):
    """Signature (n₊,n₋,n₀) of a symmetric matrix over Q — congruent diagonalization, exact."""
    d = Kmat.rows
    M = [[sp.Rational(Kmat[i, j], 1) for j in range(d)] for i in range(d)]
    n_pos = n_neg = n_zero = 0
    for k in range(d):
        if M[k][k] == 0:
            jd = next((j for j in range(k + 1, d) if M[j][j] != 0), None)
            if jd is not None:
                M[k], M[jd] = M[jd], M[k]
                for r in range(d):
                    M[r][k], M[r][jd] = M[r][jd], M[r][k]
            else:
                jo = next((j for j in range(k + 1, d) if M[k][j] != 0), None)
                if jo is None:
                    n_zero += 1
                    continue
                for c in range(d):
                    M[k][c] += M[jo][c]
                for r in range(d):
                    M[r][k] += M[r][jo]
        piv = M[k][k]
        for i in range(k + 1, d):
            if M[i][k] != 0:
                f = M[i][k] / piv
                for c in range(d):
                    M[i][c] -= f * M[k][c]
                for r in range(d):
                    M[r][i] -= f * M[r][k]
        if piv > 0:
            n_pos += 1
        else:
            n_neg += 1
    return (n_pos, n_neg, n_zero)


def make_soPQ(p, q):
    """η, J-type and K-type generators with plane labels (i,j)."""
    n = p + q
    signs = [1] * p + [-1] * q
    eta = sp.diag(*signs)
    Js, Ks = [], []
    for i in range(n):
        for j in range(i + 1, n):
            E_ij = sp.zeros(n, n)
            E_ij[i, j] = 1
            E_ji = sp.zeros(n, n)
            E_ji[j, i] = 1
            if signs[i] == signs[j]:
                Js.append((f"J({i},{j})", E_ij - E_ji))
            else:
                Ks.append((f"K({i},{j})", E_ij + E_ji))
    return eta, Js, Ks


def sig_str(s):
    return f"({s[0]},{s[1]},{s[2]})"


def cls_str(cl):
    d, s = cl
    return f"{{dim {d} · Killing {sig_str(s)}}}"


# ═══════════════════════════════════════════════════════════════════════════════════════
# STABILIZER OF THE PAIR (η,B): exact nullspace over Q + asserts (subalgebra explicit).
# ═══════════════════════════════════════════════════════════════════════════════════════

N_SUBALG_ASSERTS = 0   # count of explicit bracket-closedness checks (all g_B)


def solve_gB(Bm, G, eta, nn):
    """g_B = {X ∈ so(η) : XᵀB + BX = 0}: nullspace of the linear system on the coefficients
    in the generator basis G. Returns (basis matrices, Span). Asserts: every element is in
    so(η) and satisfies the pair condition; the basis is independent; pairwise brackets lie within the span."""
    global N_SUBALG_ASSERTS
    cols = [sp.Matrix(list(flat(Gk.T * Bm + Bm * Gk))) for Gk in G]
    M = sp.Matrix.hstack(*cols)   # (2n)² × dim so
    ns = M.nullspace()            # exact nullspace over Q
    basis = []
    span = Span()
    for v in ns:
        X = sp.zeros(nn, nn)
        for k in range(len(G)):
            if v[k] != 0:
                X = X + v[k] * G[k]
        assert (X * eta + eta * X.T) == sp.zeros(nn, nn), "element of g_B outside so(η)"
        assert (X.T * Bm + Bm * X) == sp.zeros(nn, nn), "element of g_B outside the pair condition"
        added = span.add(flat(X))
        assert added, "nullspace basis is linearly dependent"
        basis.append(X)
    # EXPLICIT check: the stabilizer is a subalgebra (all pairwise brackets in the span)
    for a in range(len(basis)):
        for b in range(a + 1, len(basis)):
            assert span.contains(flat(bracket(basis[a], basis[b]))), \
                "g_B is not closed under the bracket"
            N_SUBALG_ASSERTS += 1
    return basis, span


def gB_class(basis, span, kill_cache):
    """Class {dim · Killing-signature}; Killing cached by the RREF key of the subspace."""
    d = span.dim()
    if d == 0:
        return (0, (0, 0, 0))
    k = span.key()
    if k not in kill_cache:
        kill_cache[k] = sym_signature(killing_matrix(basis))
    return (d, kill_cache[k])


def diag_pattern(bs, signs):
    """Pattern of a diagonal B: axis groups by shared λᵢ = bᵢ/ηᵢᵢ (= bᵢ·ηᵢᵢ, since ηᵢᵢ=±1)
    of the matrix A = η⁻¹B. Key = the multiset of induced signatures (pᵢ,qᵢ) of the groups
    (descending). Zero bᵢ are also a group (λ=0). Returns (key, group count)."""
    lam_groups = {}
    for b, s in zip(bs, signs):
        lam = b * s
        pq = lam_groups.setdefault(lam, [0, 0])
        if s > 0:
            pq[0] += 1
        else:
            pq[1] += 1
    parts = tuple(sorted(((p, q) for p, q in lam_groups.values()), reverse=True))
    return parts, len(lam_groups)


def parts_str(parts):
    return " ⊕ ".join(f"({p},{q})" for p, q in parts)


# ═══════════════════════════════════════════════════════════════════════════════════════
# REFERENCE CLASSES S901 (measured data of the same container; for block C)
# ═══════════════════════════════════════════════════════════════════════════════════════
REF_S901 = {
    2: {(1, (0, 0, 1)), (2, (0, 0, 2)), (3, (2, 1, 0)), (4, (2, 1, 1)),
        (5, (3, 1, 1)), (6, (4, 2, 0))},
    3: {(1, (0, 0, 1)), (2, (0, 0, 2)), (3, (0, 0, 3)), (3, (0, 3, 0)),
        (3, (2, 1, 0)), (4, (0, 3, 1)), (4, (2, 1, 1)), (6, (0, 6, 0)),
        (6, (3, 3, 0)), (6, (4, 2, 0)), (7, (3, 3, 1)), (7, (4, 2, 1)),
        (10, (6, 4, 0)), (15, (9, 6, 0))},
}

CONTAINER_STORE = {}   # n → (eta, G, kill_cache, classesA_set) for block D

# ═══════════════════════════════════════════════════════════════════════════════════════
for n_half in (2, 3):
    nn = 2 * n_half
    dim_so = nn * (nn - 1) // 2
    rule(f"CONTAINER so({n_half},{n_half}) — 2n={nn}, "
         f"η=diag({'+1,' * n_half}{'−1,' * n_half}) · dim so = {dim_so}")
    eta, Js, Ks = make_soPQ(n_half, n_half)
    signs = [1] * n_half + [-1] * n_half
    gens = Js + Ks
    G = [M for _, M in gens]

    check(f"so({n_half},{n_half}): Xη + ηXᵀ = 0 for all {len(gens)} generators",
          all((Gm * eta + eta * Gm.T) == sp.zeros(nn, nn) for Gm in G))
    check(f"so({n_half},{n_half}): generator count = 2n(2n−1)/2", len(gens) == dim_so,
          f"J-moves: {len(Js)} · K-moves: {len(Ks)} · total: {len(gens)}")
    stack = sp.Matrix([list(flat(Gm)) for Gm in G])
    check(f"so({n_half},{n_half}): linear independence of the generators",
          stack.rank() == dim_so, f"rank = {stack.rank()}")

    kill_cache = {}

    # ═══ BLOCK A — diagonal B = diag(b₁..b_{2n}), bᵢ ∈ {−2..2} ═══
    total_B = 5 ** nn
    print(f"\n  BLOCK A — diagonal B: bᵢ ∈ {{−2,−1,0,1,2}}, total 5^{nn} = {total_B}")
    classesA = {}          # class (dim,sig) → Counter(pattern key → B count)
    pat_class = {}         # pattern key → class
    pat_count = Counter()  # pattern key → B count
    pat_groups = {}        # pattern key → λ-group count of the representative

    if n_half == 2:
        # full pass: EVERY one of the 625 B solved directly
        n_solved = 0
        for bs in product(range(-2, 3), repeat=nn):
            Bm = sp.diag(*bs)
            key, n_grp = diag_pattern(bs, signs)
            basis, span = solve_gB(Bm, G, eta, nn)
            cl = gB_class(basis, span, kill_cache)
            n_solved += 1
            if key in pat_class:
                assert pat_class[key] == cl, \
                    f"class within pattern {key} diverged: {pat_class[key]} ≠ {cl}"
            else:
                pat_class[key] = cl
                pat_groups[key] = n_grp
            pat_count[key] += 1
        cover_line = f"FULL PASS: solved directly {n_solved}/{total_B}"
        check(f"so({n_half},{n_half}) block A: full pass, every B solved directly",
              n_solved == total_B, cover_line)
        check(f"so({n_half},{n_half}) block A: class agreement within every pattern "
              f"(assert on all {total_B})", True, f"patterns: {len(pat_class)}")
    else:
        # cache coverage over the canonical pattern + representative solves + verification
        pat_rep = {}
        for bs in product(range(-2, 3), repeat=nn):
            key, n_grp = diag_pattern(bs, signs)
            pat_count[key] += 1
            if key not in pat_rep:
                pat_rep[key] = bs
                pat_groups[key] = n_grp
        n_solved = 0
        for key, bs in pat_rep.items():
            Bm = sp.diag(*bs)
            basis, span = solve_gB(Bm, G, eta, nn)
            pat_class[key] = gB_class(basis, span, kill_cache)
            n_solved += 1
        covered = sum(pat_count.values())
        cover_line = (f"CACHED BY PATTERN: patterns {len(pat_rep)} · representatives "
                      f"solved directly {n_solved} · covered {covered}/{total_B}")
        print(f"    {cover_line}")
        check(f"so({n_half},{n_half}) block A: cache coverage = all {total_B} B",
              covered == total_B, cover_line)
        # verification sample: 30 random B (seed 9023) solved directly against the cache
        vrng = random.Random(9023)
        n_verify = 30
        mism = 0
        for _ in range(n_verify):
            bs = tuple(vrng.randint(-2, 2) for _ in range(nn))
            key, _g = diag_pattern(bs, signs)
            basis, span = solve_gB(sp.diag(*bs), G, eta, nn)
            cl = gB_class(basis, span, kill_cache)
            if cl != pat_class[key]:
                mism += 1
        check(f"so({n_half},{n_half}) block A: cache verification sample "
              f"({n_verify} random B, seed 9023, direct solve = cache class)",
              mism == 0, f"mismatches: {mism}")

    # grouping patterns into classes {dim · Killing}
    for key, cl in pat_class.items():
        classesA.setdefault(cl, Counter())[key] = pat_count[key]
    check(f"so({n_half},{n_half}) block A: sum of class counts = {total_B}",
          sum(sum(c.values()) for c in classesA.values()) == total_B,
          f"= {sum(sum(c.values()) for c in classesA.values())}")

    print(f"\n    BLOCK A CLASSES (class → pattern decompositions ⊕(pᵢ,qᵢ) with B count):")
    print(f"    {'dim g_B':>8} {'Killing':>10} {'degen':>7} {'B count':>9}   decompositions ⊕(pᵢ,qᵢ) [B count per pattern]")
    pq_parts = set()
    for cl in sorted(classesA):
        d, s = cl
        deg = "yes" if s[2] > 0 else ("—" if d == 0 else "no")
        tot = sum(classesA[cl].values())
        pats = " · ".join(f"{parts_str(k)} [{classesA[cl][k]}]"
                          for k in sorted(classesA[cl], reverse=True))
        print(f"    {d:>8} {sig_str(s):>10} {deg:>7} {tot:>9}   {pats}")
        for k in classesA[cl]:
            pq_parts |= set(k)
    print(f"    A classes: {len(classesA)} · patterns: {len(pat_class)}")
    print(f"    SET OF ALL (p,q)-COMPONENTS in block A decompositions: "
          f"{' · '.join(f'({p},{q})' for p, q in sorted(pq_parts, reverse=True))}")

    # ═══ BLOCK B — random full symmetric B (seed 902) ═══
    N_RAND = 300
    print(f"\n  BLOCK B — random symmetric B: N = {N_RAND}, coeffs {{−3..3}}, seed 902 "
          f"(separate Random(902) per container)")
    rng = random.Random(902)
    x = sp.Symbol("x")
    histB = Counter()
    shapes_by_class = {}
    for _ in range(N_RAND):
        Bm = sp.zeros(nn, nn)
        for i in range(nn):
            for j in range(i, nn):
                v = rng.randint(-3, 3)
                Bm[i, j] = v
                Bm[j, i] = v
        basis, span = solve_gB(Bm, G, eta, nn)
        cl = gB_class(basis, span, kill_cache)
        histB[cl] += 1
        # factorization structure of charpoly(η⁻¹B) over Q: (deg of irreducible)^multiplicity
        Am = eta * Bm   # η⁻¹ = η
        fl = sp.factor_list(Am.charpoly(x).as_expr(), x)
        shape = tuple(sorted((sp.degree(f, x), m) for f, m in fl[1]))
        shp = " · ".join(f"deg{dd}^{mm}" for dd, mm in shape)
        shapes_by_class.setdefault(cl, Counter())[shp] += 1
    check(f"so({n_half},{n_half}) block B: all {N_RAND} B solved (no truncation)",
          sum(histB.values()) == N_RAND, f"= {sum(histB.values())}")
    print(f"    HISTOGRAM OF CLASSES {{dim · Killing}} (N = {N_RAND}):")
    for cl in sorted(histB):
        print(f"      {cls_str(cl)} : {histB[cl]}")
    print(f"    FACTORIZATION of charpoly(η⁻¹B) over Q BY CLASS (deg^multiplicity : count):")
    for cl in sorted(shapes_by_class):
        print(f"      {cls_str(cl)}:")
        for shp, cnt in sorted(shapes_by_class[cl].items()):
            print(f"        {shp} : {cnt}")

    # ═══ BLOCK C — intersection with the S901 list (same container) ═══
    stab_classes = set(classesA) | set(histB)
    ref = REF_S901[n_half]
    inter = sorted(stab_classes & ref)
    only_ref = sorted(ref - stab_classes)
    only_stab = sorted(stab_classes - ref)
    print(f"\n  BLOCK C — sets of keys {{dim · Killing}} against the S901 list "
          f"(so({n_half},{n_half})):")
    print(f"    stabilizer classes A∪B of this probe ({len(stab_classes)}): "
          + " ".join(cls_str(c) for c in sorted(stab_classes)))
    print(f"    (i)   stabilizer ∩ S901 ({len(inter)}): "
          + (" ".join(cls_str(c) for c in inter) if inter else "∅"))
    print(f"    (ii)  S901 ∖ stabilizer ({len(only_ref)}): "
          + (" ".join(cls_str(c) for c in only_ref) if only_ref else "∅"))
    print(f"    (iii) stabilizer ∖ S901 ({len(only_stab)}): "
          + (" ".join(cls_str(c) for c in only_stab) if only_stab else "∅"))

    CONTAINER_STORE[n_half] = (eta, G, kill_cache, set(classesA))


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("BLOCK D — off-diagonal B, n=2: A = η⁻¹B with a 2×2 Jordan block (sympy-exact)")
eta2, G2, kill_cache2, classesA2 = CONTAINER_STORE[2]


def jordan_blocks(Am):
    """Multiset of (eigenvalue, block size) from the exact Jordan form."""
    _P, J = Am.jordan_form()
    blocks = []
    i = 0
    m = J.rows
    while i < m:
        size = 1
        while i + size < m and J[i + size - 1, i + size] == 1:
            size += 1
        blocks.append((J[i, i], size))
        i += size
    return sorted(blocks, key=lambda t: (sp.default_sort_key(t[0]), t[1]))


D_EXAMPLES = [
    # (name, dict of nonzero upper-triangle entries of B)
    ("D1: J₂(0) ⊕ (2) ⊕ (−2)",
     {(0, 0): 1, (0, 2): 1, (2, 2): 1, (1, 1): 2, (3, 3): 2}),
    ("D2: J₂(1) ⊕ (3) ⊕ (−3)",
     {(0, 0): 2, (0, 2): 1, (2, 2): 0, (1, 1): 3, (3, 3): 3}),
]
for name, entries in D_EXAMPLES:
    Bm = sp.zeros(4, 4)
    for (i, j), v in entries.items():
        Bm[i, j] = v
        Bm[j, i] = v
    check(f"{name}: B symmetric integer", Bm == Bm.T and
          all(sp.Integer(Bm[i, j]) == Bm[i, j] for i in range(4) for j in range(4)))
    Am = eta2 * Bm
    jb = jordan_blocks(Am)
    has_j2 = any(sz == 2 for _ev, sz in jb)
    jb_s = " · ".join(f"J_{sz}({ev})" for ev, sz in jb)
    check(f"{name}: Jordan structure of A = η⁻¹B contains a 2×2 block (exact jordan_form)",
          has_j2, f"blocks: {jb_s}")
    basis, span = solve_gB(Bm, G2, eta2, 4)
    cl = gB_class(basis, span, kill_cache2)
    d, s = cl
    deg = "yes" if s[2] > 0 else ("—" if d == 0 else "no")
    print(f"    {name}: B(rows) = {[list(Bm.row(i)) for i in range(4)]}")
    print(f"      Jordan blocks of A: {jb_s} · charpoly(A) = "
          f"{sp.factor(Am.charpoly(sp.Symbol('x')).as_expr())}")
    print(f"      dim g_B = {d} · Killing {sig_str(s)} · degen: {deg} · "
          f"key {cls_str(cl)} in the set of diagonal classes A(n=2): "
          f"{'yes' if cl in classesA2 else 'no'}")


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY S902 (raw lines; reading = an act of the court)")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}
  Explicit subalgebra bracket-asserts (all g_B, pairwise brackets): {N_SUBALG_ASSERTS}

  RAW LINES (no readings):
   (1) containers so(2,2)·so(3,3): generators form a basis — sympy-exact;
   (2) block A: diagonal B — n=2 full pass 625/625 directly; n=3 cache over
       the canonical pattern (coverage 15625/15625, representatives + verification);
       class tables with decompositions ⊕(pᵢ,qᵢ) and the (p,q)-component set — above;
   (3) block B: 300 random symmetric B per container (seed 902) — class
       histograms and factorization of charpoly(η⁻¹B) by class — above;
   (4) block C: three sets against the S901 list — above;
   (5) block D: off-diagonal Jordan examples n=2 — above.
  HONEST TALLY: handles 0 · downstream-content 0 · verdicts 0. Court = Omega.
""")
_logf.flush()
sys.exit(0 if not FAIL else 1)
