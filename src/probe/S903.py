#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (Lie algebra + linear algebra; handles 0). W28-O5 — centralizers
#      c(A) = {X ∈ so(η) : [X,A]=0} of the pair (η,Ω), A = η⁻¹Ω, in so(2,2)·so(3,1)·so(3,3):
#      single coordinate generators · random integer A · A² = ±𝟙 · nilpotent A;
#      sets of classes {dim·Killing} against the S901/S902 reference lists.
#      ★BLINDNESS: the probe prints ONLY raw identities/numbers/histograms/sets; reading = an act of the court.
"""
S903 (lane A, ed.2) — W28-O5: centralizers of the pair (η,Ω) in so(2,2)·so(3,1)·so(3,3).

  Object: η = diag(+1×p, −1×q); so(η) = {X : Xη + ηXᵀ = 0}; a second form Ω
  antisymmetric on the carrier ⟺ A = η⁻¹Ω ∈ so(η) (both directions by assert:
  on the generator basis — ηX is antisymmetric; reverse — generic symbolic Ω).
  Centralizer:
    c(A) = {X ∈ so(η) : [X,A] = 0}
  — an exact nullspace of the linear system on the coefficients of X in the generator basis.
  For every c(A): dim · Killing-signature in its own basis (degenerate/non-degenerate) ·
  subalgebra assert (all pairwise brackets in the span).

  BLOCK A — c(a single coordinate generator): every J and every K individually;
    grouping of identical classes {dim · Killing}, count.
  BLOCK B — generic A: N = 300 per signature, A = random integer combination of
    generators (coeffs {−3..3}, seed 903, separate Random per signature);
    histogram of classes {dim c · Killing c}; factorization structure of
    charpoly(A) over Q by class (degrees of irreducibles · multiplicities).
  BLOCK C — A² = −𝟙 (minimal polynomial x²+1) and A² = +𝟙 (x²−1):
    systematic search within stamped bounds (nn=4: full enumeration of coeffs
    {−2..2}^6 = 15625; nn=6: coeffs {−1,0,1}^15 with the exact necessary condition
    tr(A²) = tr(±𝟙) = ±6 — the tr-Gram of generators is diagonal ±2, assert +
    verification sample of the formula on random coeffs) + symbolic identities
    (ΩᵀηΩ = −ηA² · AᵀηA = −ηA² · xᵀΩx = 0, all expanded to zero) and exact
    signatures of η and −η (Sylvester congruence). For those found: c(A) in full;
    for those not found — a raw line "not found within bounds …" with the bound stamp.
  BLOCK D — nilpotent A: rank-2 constructions A = v·(ηw)ᵀ − w·(ηv)ᵀ with
    η(v,v)=0 · η(v,w)=0 (+ a Witt-block example on so(3,3)): nilpotency
    degree (A^k = 0 and A^{k−1} ≠ 0, assert) · c(A).
  BLOCK E — sets of keys {dim · Killing} A∪B∪C∪D against the S901-terminal
    reference lists (all three signatures) and the S902-stabilizer lists
    (containers so(2,2)·so(3,3)): intersection · differences both ways.

  EXACT arithmetic: sympy Rational/Integer over Q; no tolerances. Mechanisms
  (Span/RREF-key/bracket/structure-constants/Killing/signature) — a verbatim
  copy from S900/S901/S902 (reproducibility). No silent truncation:
  every search bound is stamped in the printed lines.

Fence: the shared fence_scan helper (forbidden words — on the GUARDLINE line).
No verdicts in the probe's text; raw lines with no interpretive words.
"""
import sys
import random
from math import comb
from itertools import combinations, product
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

# ── tee: all of stdout is duplicated into S903_run.log next to the script ──
_LOG_PATH = __file__.replace("\\", "/").rsplit("/", 1)[0] + "/S903_run.log"


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
_FORBIDDEN = [r"стійк", r"обрано", r"selected", r"stable", r"причин", r"час", r"time", r"arrow", r"стріла", r"вибір", r"краще", r"злам", r"Лі-клас", r"евклід", r"поле", r"бівектор", r"сила", r"енергія", r"матерія"]   # GUARDLINE
_hits = scan_forbidden(__file__, _FORBIDDEN)
check("fence: forbidden words (list = GUARDLINE line) = 0 occurrences outside the declaration",
      not _hits, f"hits: {_hits}" if _hits else "0")
check("handles 0 (pure algebra: pair (η,Ω) · centralizer · rank · Killing · nilpotent)", True)


# ═══════════════════════════════════════════════════════════════════════════════════════
# TOOLS (verbatim copy from S900/S901/S902): exact span (echelon over Q), bracket,
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


def deg_flag(cl):
    d, s = cl
    return "yes" if s[2] > 0 else ("—" if d == 0 else "no")


# ═══════════════════════════════════════════════════════════════════════════════════════
# CENTRALIZER c(A) = {X ∈ so(η): [X,A]=0}: exact nullspace over Q + asserts.
# ═══════════════════════════════════════════════════════════════════════════════════════

N_SUBALG_ASSERTS = 0   # count of explicit bracket-closedness checks (all c(A))


def solve_cA(Am, G, eta, nn):
    """c(A): nullspace of the linear system [X,A]=0 on the coefficients of X in the generator basis G.
    Returns (basis matrices, Span). Asserts: every element is in so(η) and commutes with A;
    the basis is independent; all pairwise brackets are in the span (subalgebra — explicit)."""
    global N_SUBALG_ASSERTS
    cols = [sp.Matrix(list(flat(Gk * Am - Am * Gk))) for Gk in G]
    M = sp.Matrix.hstack(*cols)   # nn² × dim so
    ns = M.nullspace()            # exact nullspace over Q
    basis = []
    span = Span()
    for v in ns:
        X = sp.zeros(nn, nn)
        for k in range(len(G)):
            if v[k] != 0:
                X = X + v[k] * G[k]
        assert (X * eta + eta * X.T) == sp.zeros(nn, nn), "element of c(A) outside so(η)"
        assert (X * Am - Am * X) == sp.zeros(nn, nn), "element of c(A) does not commute with A"
        added = span.add(flat(X))
        assert added, "nullspace basis is linearly dependent"
        basis.append(X)
    # EXPLICIT check: c(A) is a subalgebra (all pairwise brackets in the span)
    for a in range(len(basis)):
        for b in range(a + 1, len(basis)):
            assert span.contains(flat(bracket(basis[a], basis[b]))), \
                "c(A) is not closed under the bracket"
            N_SUBALG_ASSERTS += 1
    return basis, span


def cent_class(basis, span, kill_cache):
    """Class {dim · Killing-signature}; Killing cached by the RREF key of the subspace."""
    d = span.dim()
    if d == 0:
        return (0, (0, 0, 0))
    k = span.key()
    if k not in kill_cache:
        kill_cache[k] = sym_signature(killing_matrix(basis))
    return (d, kill_cache[k])


def charpoly_shape(Am):
    """Factorization structure of charpoly(A) over Q: sorted (deg of irreducible, multiplicity)."""
    xsym = sp.Symbol("x")
    fl = sp.factor_list(Am.charpoly(xsym).as_expr(), xsym)
    parts = sorted((int(sp.degree(f, xsym)), int(m)) for f, m in fl[1])
    return " · ".join(f"deg{dd}^{mm}" for dd, mm in parts)


def nilp_degree(Am, cap):
    """Smallest k ≤ cap with A^k = 0 (None if none); exact arithmetic."""
    P = sp.eye(Am.rows)
    for k in range(1, cap + 1):
        P = P * Am
        if P == sp.zeros(Am.rows, Am.rows):
            return k
    return None


def rank2_A(v, w, eta):
    """A = v·(ηw)ᵀ − w·(ηv)ᵀ (rank ≤ 2, ηA antisymmetric by construction)."""
    return v * (eta * w).T - w * (eta * v).T


def int_mat(M):
    """sympy matrix with integer entries → tuple of tuples of python int."""
    return tuple(tuple(int(M[i, j]) for j in range(M.cols)) for i in range(M.rows))


def gen_nonzeros(Mint):
    """Nonzero entries of an integer matrix: tuple (i,j,value)."""
    return tuple((i, j, v) for i, row in enumerate(Mint) for j, v in enumerate(row) if v)


def sq_is(A, s, nn):
    """A² == s·𝟙 (python-int matrix A as a list of rows; early exit)."""
    for r in range(nn):
        Ar = A[r]
        for c in range(nn):
            t = s if r == c else 0
            acc = 0
            for k in range(nn):
                a = Ar[k]
                if a:
                    acc += a * A[k][c]
            if acc != t:
                return False
    return True


# ═══════════════════════════════════════════════════════════════════════════════════════
# REFERENCE LISTS (measured data; for block E)
# ═══════════════════════════════════════════════════════════════════════════════════════
REF_S901 = {
    (2, 2): {(1, (0, 0, 1)), (2, (0, 0, 2)), (3, (2, 1, 0)), (4, (2, 1, 1)),
             (5, (3, 1, 1)), (6, (4, 2, 0))},
    (3, 1): {(1, (0, 0, 1)), (2, (0, 0, 2)), (3, (0, 3, 0)), (3, (2, 1, 0)),
             (4, (1, 1, 2)), (6, (3, 3, 0))},
    (3, 3): {(1, (0, 0, 1)), (2, (0, 0, 2)), (3, (0, 0, 3)), (3, (0, 3, 0)),
             (3, (2, 1, 0)), (4, (0, 3, 1)), (4, (2, 1, 1)), (6, (0, 6, 0)),
             (6, (3, 3, 0)), (6, (4, 2, 0)), (7, (3, 3, 1)), (7, (4, 2, 1)),
             (10, (6, 4, 0)), (15, (9, 6, 0))},
}
REF_S902 = {
    (2, 2): {(0, (0, 0, 0)), (1, (0, 0, 1)), (2, (0, 0, 2)), (3, (2, 1, 0)),
             (6, (4, 2, 0))},
    (3, 3): REF_S901[(3, 3)] | {(0, (0, 0, 0))},
}

SIGS = [(2, 2), (3, 1), (3, 3)]
CLASS_SETS = {}          # signature → set of classes A∪B∪C∪D
SEARCH_STAMPS = []       # search-bound stamps for block C (for the summary)

# ═══════════════════════════════════════════════════════════════════════════════════════
for (p, q) in SIGS:
    nn = p + q
    dim_so = nn * (nn - 1) // 2
    rule(f"SECTION so({p},{q}) — n={nn}, η=diag({'+1,' * p}{'−1,' * q}) · dim so = {dim_so}")
    eta, Js, Ks = make_soPQ(p, q)
    gens = Js + Ks
    G = [M for _, M in gens]
    gnames = [nm for nm, _ in gens]

    check(f"so({p},{q}): Xη + ηXᵀ = 0 for all {len(gens)} generators (sympy-exact)",
          all((Gm * eta + eta * Gm.T) == sp.zeros(nn, nn) for Gm in G))
    check(f"so({p},{q}): generator count = n(n−1)/2", len(gens) == dim_so,
          f"J-moves: {len(Js)} · K-moves: {len(Ks)} · total: {len(gens)}")
    stack = sp.Matrix([list(flat(Gm)) for Gm in G])
    check(f"so({p},{q}): linear independence of the generators", stack.rank() == dim_so,
          f"rank = {stack.rank()}")

    # ── EQUIVALENCE (η,Ω) ⟺ A = η⁻¹Ω ∈ so(η) ──
    check(f"so({p},{q}): ηX antisymmetric for EVERY generator X (direction so(η) ⟹ Ω)",
          all((eta * Gm).T == -(eta * Gm) for Gm in G))
    _wsyms = [sp.Symbol(f"w{i}_{j}") for i in range(nn) for j in range(i + 1, nn)]
    OmegaS = sp.zeros(nn, nn)
    _t = 0
    for i in range(nn):
        for j in range(i + 1, nn):
            OmegaS[i, j] = _wsyms[_t]
            OmegaS[j, i] = -_wsyms[_t]
            _t += 1
    AS = eta * OmegaS   # η⁻¹ = η (diagonal ±1)
    check(f"so({p},{q}): generic antisymmetric Ω ⟹ A = η⁻¹Ω ∈ so(η) "
          f"(symbolic, {len(_wsyms)} symbols)",
          sp.expand(AS * eta + eta * AS.T) == sp.zeros(nn, nn))

    kill_cache = {}
    cls_set = set()

    # ═══ BLOCK A — centralizers of single coordinate generators ═══
    print(f"\n  BLOCK A — c(a single generator): {len(Js)} J-moves + {len(Ks)} K-moves "
          f"(each individually)")
    groupsA = {}
    for nm, X in gens:
        basis, span = solve_cA(X, G, eta, nn)
        assert span.contains(flat(X)), "generator not in its own centralizer"
        cl = cent_class(basis, span, kill_cache)
        groupsA.setdefault(cl, []).append(nm)
    check(f"so({p},{q}) block A: all {len(gens)} single generators solved "
          f"(X ∈ c(X) — assert on each)",
          sum(len(v) for v in groupsA.values()) == len(gens))
    print(f"    CLASSES (grouping identical):")
    print(f"    {'dim c':>6} {'Killing':>10} {'degen':>6} {'count':>6}   generators")
    for cl in sorted(groupsA):
        names = groupsA[cl]
        print(f"    {cl[0]:>6} {sig_str(cl[1]):>10} {deg_flag(cl):>6} {len(names):>6}   "
              + " ".join(names))
    cls_set |= set(groupsA)

    # ═══ BLOCK B — generic A (random integer combinations, seed 903) ═══
    N_RAND = 300
    print(f"\n  BLOCK B — generic A: N = {N_RAND}, A = Σ cᵢ·Gᵢ, cᵢ ∈ {{−3..3}}, seed 903 "
          f"(separate Random(903) per signature)")
    rng = random.Random(903)
    histB = Counter()
    shapes_by_class = {}
    for _ in range(N_RAND):
        cs = [rng.randint(-3, 3) for _ in range(dim_so)]
        Am = sp.zeros(nn, nn)
        for c, Gk in zip(cs, G):
            if c:
                Am = Am + c * Gk
        basis, span = solve_cA(Am, G, eta, nn)
        if any(cs):
            assert span.contains(flat(Am)), "A not in its own centralizer"
        cl = cent_class(basis, span, kill_cache)
        histB[cl] += 1
        shapes_by_class.setdefault(cl, Counter())[charpoly_shape(Am)] += 1
    check(f"so({p},{q}) block B: all {N_RAND} A solved (no truncation)",
          sum(histB.values()) == N_RAND, f"= {sum(histB.values())}")
    print(f"    HISTOGRAM OF CLASSES {{dim c · Killing c}} (N = {N_RAND}):")
    for cl in sorted(histB):
        print(f"      {cls_str(cl)} · degen: {deg_flag(cl)} : {histB[cl]}")
    print(f"    FACTORIZATION of charpoly(A) over Q BY CLASS (deg^multiplicity : count):")
    for cl in sorted(shapes_by_class):
        print(f"      {cls_str(cl)}:")
        for shp, cnt in sorted(shapes_by_class[cl].items()):
            print(f"        {shp} : {cnt}")
    cls_set |= set(histB)

    # ═══ BLOCK C — special A: A² = −𝟙 (x²+1) and A² = +𝟙 (x²−1), exact ═══
    print(f"\n  BLOCK C — A² = ±𝟙 on so({p},{q}):")
    # symbolic identities (generic antisymmetric Ω, A = ηΩ; expand to zero):
    check(f"so({p},{q}) block C: identity ΩᵀηΩ + ηA² = 0 (symbolic, generic Ω)",
          sp.expand(OmegaS.T * eta * OmegaS + eta * AS * AS) == sp.zeros(nn, nn))
    check(f"so({p},{q}) block C: identity AᵀηA + ηA² = 0 (symbolic, generic Ω)",
          sp.expand(AS.T * eta * AS + eta * AS * AS) == sp.zeros(nn, nn))
    _xs = sp.Matrix([sp.Symbol(f"x{i}") for i in range(nn)])
    check(f"so({p},{q}) block C: identity xᵀΩx = 0 (symbolic; xᵀηAx = xᵀΩx)",
          sp.expand((_xs.T * OmegaS * _xs)[0, 0]) == 0)
    sig_eta = sym_signature(sp.Matrix(eta))
    sig_meta = sym_signature(sp.Matrix(-eta))
    print(f"    sym_signature(η) = {sig_str(sig_eta)} · sym_signature(−η) = {sig_str(sig_meta)} "
          f"· match: {'yes' if sig_eta == sig_meta else 'no'}")
    print(f"    LINE (from the identities above): A² = +𝟙 ⟹ ΩᵀηΩ = −η (Ω = ηA invertible, since A² = 𝟙); "
          f"Sylvester: congruence by an invertible matrix preserves the signature")
    print(f"    LINE (from the identities above): A² = −𝟙 ⟹ AᵀηA = η and xᵀηAx = 0 ∀x ⟹ "
          f"the Gram of the pair (x, Ax) = [[t,0],[0,t]], t = xᵀηx; p mod 2 = {p % 2} · q mod 2 = {q % 2}")

    # ── systematic search within stamped bounds ──
    G_int = [int_mat(Gm) for Gm in G]
    G_nz = [gen_nonzeros(Gi) for Gi in G_int]
    found = {}   # s (−1/+1) → (n_cand, n_found, first coeff vector or None)

    if nn == 4:
        stamp = (f"so({p},{q}): SEARCH BOUNDARY (stamp): coeffs ∈ {{−2..2}}^{dim_so} — "
                 f"full enumeration 5^{dim_so} = {5 ** dim_so}")
        print(f"    {stamp}")
        SEARCH_STAMPS.append(stamp)
        n_enum = 0
        res = {-1: [0, None], 1: [0, None]}
        for cs in product(range(-2, 3), repeat=dim_so):
            A = [[0] * nn for _ in range(nn)]
            for c, nz in zip(cs, G_nz):
                if c:
                    for (i, j, v) in nz:
                        A[i][j] += c * v
            n_enum += 1
            for s in (-1, 1):
                if sq_is(A, s, nn):
                    res[s][0] += 1
                    if res[s][1] is None:
                        res[s][1] = list(cs)
        check(f"so({p},{q}) block C: enumerated {n_enum} = 5^{dim_so} coeff vectors "
              f"(full enumeration, no truncation)", n_enum == 5 ** dim_so)
        for s in (-1, 1):
            found[s] = (n_enum, res[s][0], res[s][1])
    else:
        # nn=6: exact necessary condition tr(A²) = tr(s·𝟙) = 6s; tr-Gram of the generators
        T = [[sum(G_int[a][i][j] * G_int[b][j][i] for i in range(nn) for j in range(nn))
              for b in range(dim_so)] for a in range(dim_so)]
        nJ_tot = len(Js)
        check(f"so({p},{q}) block C: tr-Gram of the generators is diagonal; diagonal −2 (J) / +2 (K)",
              all(T[a][b] == 0 for a in range(dim_so) for b in range(dim_so) if a != b)
              and all(T[a][a] == -2 for a in range(nJ_tot))
              and all(T[a][a] == 2 for a in range(nJ_tot, dim_so)))
        vrng = random.Random(9031)
        n_ver = 20
        ver_ok = True
        for _ in range(n_ver):
            cs = [vrng.randint(-1, 1) for _ in range(dim_so)]
            A = [[0] * nn for _ in range(nn)]
            for c, nz in zip(cs, G_nz):
                if c:
                    for (i, j, v) in nz:
                        A[i][j] += c * v
            tr_direct = sum(sum(A[r][k] * A[k][r] for k in range(nn)) for r in range(nn))
            tr_formula = sum(c * c * T[k][k] for k, c in enumerate(cs))
            ver_ok &= (tr_direct == tr_formula)
        check(f"so({p},{q}) block C: formula tr(A²) = Σ cᵢ²·tr(Gᵢ²) — verification sample "
              f"of {n_ver} random coeff vectors (seed 9031)", ver_ok)
        for s in (-1, 1):
            # condition: −2·(count of nonzero J) + 2·(count of nonzero K) = 6s (coeffs ∈ {−1,0,1})
            need = (f"(count of nonzero J) − (count of nonzero K) = 3" if s == -1
                    else f"(count of nonzero K) − (count of nonzero J) = 3")
            n_cand = 0
            n_found = 0
            first = None
            expect = 0
            for nJ in range(0, len(Js) + 1):
                nK = nJ - 3 if s == -1 else nJ + 3
                if nK < 0 or nK > len(Ks):
                    continue
                expect += comb(len(Js), nJ) * comb(len(Ks), nK) * 2 ** (nJ + nK)
                for Jsup in combinations(range(len(Js)), nJ):
                    entsJ = [G_nz[t] for t in Jsup]
                    for Ksup in combinations(range(len(Ks)), nK):
                        ents = entsJ + [G_nz[nJ_tot + t] for t in Ksup]
                        m = nJ + nK
                        for sgn in product((1, -1), repeat=m):
                            A = [[0] * nn for _ in range(nn)]
                            for sv, ent in zip(sgn, ents):
                                for (i, j, v) in ent:
                                    A[i][j] += sv * v
                            n_cand += 1
                            if sq_is(A, s, nn):
                                n_found += 1
                                if first is None:
                                    coeffs = [0] * dim_so
                                    for sv, t in zip(sgn[:nJ], Jsup):
                                        coeffs[t] = sv
                                    for sv, t in zip(sgn[nJ:], Ksup):
                                        coeffs[nJ_tot + t] = sv
                                    first = coeffs
            stamp = (f"so({p},{q}) target A²={'−' if s == -1 else '+'}𝟙: SEARCH BOUNDARY (stamp): "
                     f"coeffs ∈ {{−1,0,1}}^{dim_so} — full domain 3^{dim_so} = {3 ** dim_so}; "
                     f"exact necessary condition tr(A²) = {6 * s} ⟹ {need}; candidates after "
                     f"the condition: {n_cand} (all enumerated); the rest of the domain is cut by the exact condition")
            print(f"    {stamp}")
            SEARCH_STAMPS.append(stamp)
            check(f"so({p},{q}) block C target A²={'−' if s == -1 else '+'}𝟙: candidate count "
                  f"= combinatorial formula", n_cand == expect, f"{n_cand} = {expect}")
            found[s] = (n_cand, n_found, first)

    # ── search report + analysis of what was found ──
    for s, label in ((-1, "A² = −𝟙 (minimal polynomial x²+1)"),
                     (1, "A² = +𝟙 (x²−1)")):
        n_cand, n_found, first = found[s]
        print(f"\n    TARGET {label}: enumerated {n_cand} · found {n_found}")
        if first is None:
            box = ("{−2..2}^" + str(dim_so) + f" (full enumeration {5 ** dim_so})") if nn == 4 \
                else ("{−1,0,1}^" + str(dim_so) + f" (domain 3^{dim_so}, tr-condition, "
                      f"{n_cand} candidates)")
            print(f"    not found within coeffs ∈ {box}")
        else:
            Asp = sp.zeros(nn, nn)
            for c, Gk in zip(first, G):
                if c:
                    Asp = Asp + c * Gk
            assert (Asp * eta + eta * Asp.T) == sp.zeros(nn, nn), "found A outside so(η)"
            assert Asp * Asp == s * sp.eye(nn), "found A: A² ≠ s·𝟙"
            if s == 1:
                assert Asp != sp.eye(nn) and Asp != -sp.eye(nn), \
                    "minimal polynomial is not x²−1"
            print(f"      first found (lexicographic enumeration order): "
                  f"coeffs = {first}")
            print(f"      A(rows) = {[list(Asp.row(i)) for i in range(nn)]}")
            print(f"      charpoly(A) = {sp.factor(Asp.charpoly(sp.Symbol('x')).as_expr())}")
            basis, span = solve_cA(Asp, G, eta, nn)
            cl = cent_class(basis, span, kill_cache)
            print(f"      dim c(A) = {cl[0]} · Killing {sig_str(cl[1])} · degen: {deg_flag(cl)}")
            cls_set.add(cl)
    check(f"so({p},{q}) block C: both targets processed (found, or bound stamp + raw line)",
          set(found) == {-1, 1})

    # ═══ BLOCK D — nilpotent A ═══
    print(f"\n  BLOCK D — nilpotent A on so({p},{q}):")
    e = [sp.eye(nn).col(i) for i in range(nn)]
    if (p, q) == (2, 2):
        D_EX = [
            ("D1: rank-2, v=e0+e2 · w=e1", rank2_A(e[0] + e[2], e[1], eta)),
            ("D2: rank-2, v=e0+e2 · w=e1+e3", rank2_A(e[0] + e[2], e[1] + e[3], eta)),
            ("D3: rank-2, v=e1+e3 · w=e0", rank2_A(e[1] + e[3], e[0], eta)),
        ]
    elif (p, q) == (3, 1):
        D_EX = [
            ("D1: rank-2, v=e2+e3 · w=e0", rank2_A(e[2] + e[3], e[0], eta)),
            ("D2: rank-2, v=e2+e3 · w=e0+e1", rank2_A(e[2] + e[3], e[0] + e[1], eta)),
            ("D3: rank-2, v=e0+e3 · w=e1", rank2_A(e[0] + e[3], e[1], eta)),
        ]
    else:
        aW = sp.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
        bW = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])
        Xw = sp.Matrix(sp.BlockMatrix([[aW, bW], [sp.zeros(3, 3), -aW.T]]))
        SW = sp.Matrix([[1, 0, 0, 1, 0, 0],
                        [0, 1, 0, 0, 1, 0],
                        [0, 0, 1, 0, 0, 1],
                        [1, 0, 0, -1, 0, 0],
                        [0, 1, 0, 0, -1, 0],
                        [0, 0, 1, 0, 0, -1]])
        etaW = SW.T * eta * SW
        assert Xw.T * etaW + etaW * Xw == sp.zeros(6, 6), "Witt block outside so(η_w)"
        A_witt = SW * Xw * SW.inv()
        D_EX = [
            ("D1: rank-2, v=e0+e3 · w=e2", rank2_A(e[0] + e[3], e[2], eta)),
            ("D2: rank-2, v=e0+e3 · w=e1+e4", rank2_A(e[0] + e[3], e[1] + e[4], eta)),
            ("D3: Witt block [[a,b],[0,−aᵀ]], a=E12+E23 · b=E23−E32 (conjugation S)", A_witt),
        ]
    for nm, Am in D_EX:
        check(f"so({p},{q}) {nm.split(':')[0]}: A ∈ so(η) (ηA antisymmetric) and A ≠ 0",
              (eta * Am).T == -(eta * Am) and Am != sp.zeros(nn, nn))
        deg_n = nilp_degree(Am, nn)
        assert deg_n is not None, "A is not nilpotent within k ≤ nn"
        Pk = sp.eye(nn)
        for _ in range(deg_n - 1):
            Pk = Pk * Am
        check(f"so({p},{q}) {nm.split(':')[0]}: A^{deg_n} = 0 and A^{deg_n - 1} ≠ 0 "
              f"(nilpotency degree {deg_n})",
              Pk != sp.zeros(nn, nn))
        basis, span = solve_cA(Am, G, eta, nn)
        assert span.contains(flat(Am)), "A not in its own centralizer"
        cl = cent_class(basis, span, kill_cache)
        print(f"    {nm}")
        print(f"      A(rows) = {[list(Am.row(i)) for i in range(nn)]}")
        print(f"      nilpotency degree = {deg_n} · dim c(A) = {cl[0]} · "
              f"Killing {sig_str(cl[1])} · degen: {deg_flag(cl)}")
        cls_set.add(cl)

    CLASS_SETS[(p, q)] = cls_set

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("BLOCK E — sets of keys {dim · Killing} against the S901/S902 reference lists (raw)")
for (p, q) in SIGS:
    cent = CLASS_SETS[(p, q)]
    ref1 = REF_S901[(p, q)]
    print(f"\n  so({p},{q}):")
    print(f"    (i)   centralizer classes A∪B∪C∪D ({len(cent)}): "
          + " ".join(cls_str(c) for c in sorted(cent)))
    inter = sorted(cent & ref1)
    only_c = sorted(cent - ref1)
    only_r = sorted(ref1 - cent)
    print(f"    (ii)  ∩ S901 ({len(inter)}): "
          + (" ".join(cls_str(c) for c in inter) if inter else "∅"))
    print(f"    (iii) centralizer ∖ S901 ({len(only_c)}): "
          + (" ".join(cls_str(c) for c in only_c) if only_c else "∅"))
    print(f"    (iv)  S901 ∖ centralizer ({len(only_r)}): "
          + (" ".join(cls_str(c) for c in only_r) if only_r else "∅"))
    if (p, q) in REF_S902:
        ref2 = REF_S902[(p, q)]
        inter2 = sorted(cent & ref2)
        only_c2 = sorted(cent - ref2)
        only_r2 = sorted(ref2 - cent)
        print(f"    (ii')  ∩ S902 ({len(inter2)}): "
              + (" ".join(cls_str(c) for c in inter2) if inter2 else "∅"))
        print(f"    (iii') centralizer ∖ S902 ({len(only_c2)}): "
              + (" ".join(cls_str(c) for c in only_c2) if only_c2 else "∅"))
        print(f"    (iv')  S902 ∖ centralizer ({len(only_r2)}): "
              + (" ".join(cls_str(c) for c in only_r2) if only_r2 else "∅"))
check("block E: three signatures against S901; containers so(2,2)·so(3,3) additionally against S902",
      True)

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY S903 (raw lines; reading = an act of the court)")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}
  Explicit subalgebra bracket-asserts (all c(A), pairwise brackets): {N_SUBALG_ASSERTS}

  SEARCH BOUND STAMPS (block C):""")
for st in SEARCH_STAMPS:
    print(f"   · {st}")
print(f"""
  RAW LINES (no readings):
   (1) so(2,2)·so(3,1)·so(3,3): generators form a basis + equivalence
       (η,Ω-antisymmetric) ⟺ A = η⁻¹Ω ∈ so(η) — both directions, sympy-exact;
   (2) block A: centralizers of all single J/K — classes grouped, count above;
   (3) block B: 300 random integer A per signature (seed 903) — class
       histograms and factorization of charpoly(A) by class — above;
   (4) block C: A² = ±𝟙 — searches within stamped bounds + symbolic identities +
       signatures of η/−η; A found with its c(A), those not found — raw bound lines — above;
   (5) block D: nilpotent A (rank-2 + Witt block) with degrees — above;
   (6) block E: sets of keys against S901/S902 — above.
  HONEST TALLY: handles 0 · verdicts 0. Court = Omega.
""")
_logf.flush()
sys.exit(0 if not FAIL else 1)
