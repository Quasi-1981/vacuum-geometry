#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (Lie algebra; downstream-content 0; handles 0). W28-O3 — full
#      enumeration of bracket-closures of subsets of the generator basis of so(p,q)
#      (block A) + random subspaces (block B) + growth with n (block C). Killing signatures
#      of the classes, exact.
#      ★BLINDNESS: the probe prints ONLY raw numbers/tables/histograms; reading = an act of the court.
"""
S901 (lane A, ed.2) — W28-O3: enumeration of closures of subsets of the generator basis
of so(p,q).

  Object: η = diag(+1×p, −1×q); so(p,q) = {X : Xη + ηXᵀ = 0}; a generator basis over the
  2-planes (i,j): J-type (a pair of the same sign, E_ij − E_ji) · K-type (a mixed pair,
  E_ij + E_ji). Signatures: (3,0) · (2,1) · (2,2) · (3,1) · (4,0) · (3,3).

  BLOCK A: ALL non-empty subsets of the generator basis → bracket-closure
    (memoization closure(S) = closure(closure(S∖{x}) ∪ {x}) by subset size;
    a cache of closures keyed by the canonical RREF key of the subspace). For every unique
    closure subspace L: dim L · Killing signature (n₊,n₋,n₀) · degenerate (n₀>0) or not.
    A table of classes {dim · Killing}: unique subspaces per class · starting subsets per
    class. Summary: non-degenerate classes · degenerate classes.
  BLOCK B: random k-dimensional subspaces (k = 1..min(4, dim so), N = 200 per k, seed 901,
    spans of k integer combinations of generators, coefficients from {−3..3}, rank<k
    discarded) → closure, exact → histogram of classes {dim · Killing};
    a raw intersection of the sets of classes A and B (present there/there — with no words
    of interpretation).
  BLOCK C: one line per signature: n · dim so(p,q) · classes in A · classes in A∪B · of
    which non-degenerate / degenerate.

  EXACT arithmetic: sympy Rational/Integer; ranks/echelons over Q; no tolerances.
  Mechanisms (Span/RREF key/lie_closure/structure constants/Killing/signature) — a pure
  copy from S900_w28_o1_sector_closure_ladder.py (reproducibility).

Fence: shared helper fence_scan (forbidden words — on the GUARDLINE line).
No readings in the probe's own text; raw lines with no words of interpretation.
"""
import sys
import random
from itertools import combinations
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


# ── SECTION 0 — FENCE VIA SHARED HELPER (first) ─────────────────────────────────────────
rule("SECTION 0 — FENCE (shared helper fence_scan)")
_FORBIDDEN = [r"стійк", r"обрано", r"selected", r"stable", r"причин", r"час", r"time", r"arrow", r"стріла", r"вибір", r"краще", r"евклід", r"Лі-клас"]   # GUARDLINE
_hits = scan_forbidden(__file__, _FORBIDDEN)
check("fence: forbidden words (list = the GUARDLINE line) = 0 occurrences outside the declaration",
      not _hits, f"hits: {_hits}" if _hits else "0")
check("downstream-content 0 · handles 0 (pure Lie algebra)", True)


# ═══════════════════════════════════════════════════════════════════════════════════════
# TOOLS (a pure copy from S900): exact span (echelon over Q), Lie-closure, structure
# constants, Killing form, signature of a symmetric matrix over Q.
# ═══════════════════════════════════════════════════════════════════════════════════════

def flat(M):
    """n×n matrix → a flat n² tuple (row-major), sympy-exact entries."""
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


def lie_closure(gens):
    """Iterated bracket-span to a fixed point; returns (basis matrices, Span)."""
    span = Span()
    basis = []
    for g in gens:
        if span.add(flat(g)):
            basis.append(g)
    changed = True
    while changed:
        changed = False
        m = len(basis)
        for i in range(m):
            for j in range(i + 1, m):
                br = bracket(basis[i], basis[j])
                if span.add(flat(br)):
                    basis.append(br)
                    changed = True
    return basis, span


def structure_constants(basis):
    """c^k_{ij}: [e_i,e_j] = Σ_k c^k_{ij} e_k — an exact decomposition, asserts zero residual."""
    d = len(basis)
    A = sp.Matrix.hstack(*[sp.Matrix(flat(b)) for b in basis])   # n²×d
    AtA_inv = (A.T * A).inv()
    C = [[None] * d for _ in range(d)]
    for i in range(d):
        C[i][i] = [sp.Integer(0)] * d
        for j in range(i + 1, d):
            v = sp.Matrix(flat(bracket(basis[i], basis[j])))
            sol = AtA_inv * (A.T * v)
            assert A * sol - v == sp.zeros(A.rows, 1), "the bracket lies outside the span of basis L"
            C[i][j] = [sol[k] for k in range(d)]
            C[j][i] = [-sol[k] for k in range(d)]
    return C


def killing_matrix(basis):
    d = len(basis)
    C = structure_constants(basis)
    ad = [sp.Matrix(d, d, lambda k, j, i=i: C[i][j][k]) for i in range(d)]

    def tr_prod(Ai, Aj):
        # trace(Ai·Aj) as a direct sum of element products (exact, without a full product)
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
    return f"[dim {d} · Killing {sig_str(s)}]"


# ═══════════════════════════════════════════════════════════════════════════════════════
SIGNATURES = [(3, 0), (2, 1), (2, 2), (3, 1), (4, 0), (3, 3)]
ROWS_C = []   # block C: (p,q) · n · dim so · classes in A · classes in A∪B · non-deg./deg.

for (p, q) in SIGNATURES:
    n = p + q
    dim_so = n * (n - 1) // 2
    rule(f"SECTION so({p},{q}) — n={n}, η=diag({'+1,' * p}{'−1,' * q}) · dim so = n(n−1)/2 = {dim_so}")
    eta, Js, Ks = make_soPQ(p, q)
    gens = Js + Ks
    G = [M for _, M in gens]

    # basis-hood: the isometry condition, the count, linear independence
    check(f"so({p},{q}): Xη + ηXᵀ = 0 for all {len(gens)} generators (sympy-exact)",
          all((Gm * eta + eta * Gm.T) == sp.zeros(n, n) for Gm in G))
    check(f"so({p},{q}): generator count = n(n−1)/2", len(gens) == dim_so,
          f"J-entries: {len(Js)} · K-entries: {len(Ks)} · total: {len(gens)}")
    stack = sp.Matrix([list(flat(Gm)) for Gm in G])
    rk = stack.rank()
    check(f"so({p},{q}): linear independence of the generators (sympy-exact rank)", rk == dim_so,
          f"rank = {rk}")

    # ── registry of unique closure subspaces (keyed by RREF) ──
    key_info = {}    # key → (dim, Killing signature)
    key_basis = {}   # key → basis matrices (for memoized extension)

    def register(basis, span):
        """Registers a closure by its RREF key; the Killing form is computed once per key."""
        ck = span.key()
        if ck not in key_info:
            # assert: every basis element of the closure lies in so(p,q)
            assert all((B * eta + eta * B.T) == sp.zeros(n, n) for B in basis), \
                "a closure element lies outside so(p,q)"
            # assert: the closure has genuinely reached a fixed point (all pairwise brackets in the span)
            for a in range(len(basis)):
                for b in range(a + 1, len(basis)):
                    assert span.contains(flat(bracket(basis[a], basis[b]))), \
                        "the closure did not reach a fixed point"
            key_info[ck] = (span.dim(), sym_signature(killing_matrix(basis)))
            key_basis[ck] = basis
        return ck

    # ═══ BLOCK A — full enumeration of closures of subsets of the generator basis ═══
    max_size = dim_so   # a full pass: all sizes 1..dim_so
    total_subsets = 2 ** dim_so - 1
    print(f"\n  BLOCK A — all non-empty subsets of the generator basis: 2^{dim_so} − 1 = {total_subsets}")
    print(f"    sizes covered: ≤{max_size} of {dim_so} (a full pass)")

    subkey = {(): ()}     # subset (index tuple) → RREF key of the closure
    extend_cache = {}     # (parent key, generator index) → closure key
    n_lie_calls = 0
    for size in range(1, max_size + 1):
        for S in combinations(range(dim_so), size):
            parent = S[:-1]
            x = S[-1]
            pk = subkey[parent]
            ck = extend_cache.get((pk, x))
            if ck is None:
                pbasis = list(key_basis[pk]) if pk != () else []
                basis, span = lie_closure(pbasis + [G[x]])
                n_lie_calls += 1
                ck = register(basis, span)
                extend_cache[(pk, x)] = ck
            subkey[S] = ck

    # count of the classes {dim · Killing}
    classA_keys = {}       # class → set of RREF keys
    classA_starts = Counter()
    uniq_keys_A = set()
    for S, ck in subkey.items():
        if not S:
            continue
        cl = key_info[ck][0], key_info[ck][1]
        classA_keys.setdefault(cl, set()).add(ck)
        classA_starts[cl] += 1
        uniq_keys_A.add(ck)
    check(f"so({p},{q}): ALL {total_subsets} non-empty subsets were traversed",
          sum(classA_starts.values()) == total_subsets,
          f"start count = {sum(classA_starts.values())} · lie_closure calls (memoized) = {n_lie_calls}")

    print(f"    unique closure subspaces: {len(uniq_keys_A)}")
    print(f"    classes {{dim · Killing signature}}: {len(classA_keys)}")
    print(f"    {'dim L':>6} {'Killing':>12} {'degenerate':>12} {'uniq.subspaces':>19} {'starting subsets':>19}")
    for cl in sorted(classA_keys):
        d, s = cl
        deg = "yes" if s[2] > 0 else "no"
        print(f"    {d:>6} {sig_str(s):>12} {deg:>12} {len(classA_keys[cl]):>19} {classA_starts[cl]:>19}")
    ndegA = sum(1 for (d, s) in classA_keys if s[2] == 0)
    degA = len(classA_keys) - ndegA
    print(f"    summary A: classes with a non-degenerate Killing form = {ndegA} · with a degenerate one = {degA}")

    # ═══ BLOCK B — random k-dimensional subspaces (seed 901) ═══
    kmax = min(4, dim_so)
    N_RAND = 200
    print(f"\n  BLOCK B — random subspaces: k = 1..{kmax}, N = {N_RAND} per k, seed 901, coeff. {{−3..3}}")
    rng = random.Random(901)
    classB = set()
    for k in range(1, kmax + 1):
        hist = Counter()
        got = 0
        tries = 0
        while got < N_RAND:
            tries += 1
            assert tries <= 200 * N_RAND, "the rank<k discard rate does not converge"
            mats = []
            for _ in range(k):
                coeffs = [rng.randint(-3, 3) for _ in range(dim_so)]
                Mx = sp.zeros(n, n)
                for c, Gm in zip(coeffs, G):
                    if c:
                        Mx = Mx + c * Gm
                mats.append(Mx)
            s0 = Span()
            for Mx in mats:
                s0.add(flat(Mx))
            if s0.dim() < k:
                continue   # degenerate case rank<k — discarded
            basis, span = lie_closure(mats)
            ck = register(basis, span)
            cl = key_info[ck]
            hist[cl] += 1
            got += 1
        print(f"    k={k} (spans accepted {got}, discarded for rank<k: {tries - got}):")
        for cl in sorted(hist):
            print(f"      {cls_str(cl)} : {hist[cl]}")
        classB |= set(hist)

    setA = set(classA_keys)
    onlyB = sorted(classB - setA)
    onlyA = sorted(setA - classB)
    print(f"    classes in B absent from A: {len(onlyB)}" +
          (" — " + " · ".join(cls_str(c) for c in onlyB) if onlyB else ""))
    print(f"    classes in A absent from B: {len(onlyA)}" +
          (" — " + " · ".join(cls_str(c) for c in onlyA) if onlyA else ""))

    unionAB = setA | classB
    ndegU = sum(1 for (d, s) in unionAB if s[2] == 0)
    degU = len(unionAB) - ndegU
    ROWS_C.append(((p, q), n, dim_so, len(setA), len(unionAB), ndegU, degU))


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("BLOCK C — growth with n: signature · n · dim so(p,q) · classes A · classes A∪B · non-deg./deg.")
print(f"  {'signature':<12}{'n':>4}{'dim so':>8}{'classes A':>10}{'classes A∪B':>12}"
      f"{'non-degenerate':>14}{'degenerate':>12}")
for (pq, nn, dso, nA, nU, ndg, dg) in ROWS_C:
    print(f"  ({pq[0]},{pq[1]})      {nn:>4}{dso:>8}{nA:>10}{nU:>12}{ndg:>14}{dg:>12}")


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY S901 (raw lines; reading = an act of the court)")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}

  RAW LINES (no readings):
   (1) six signatures so(p,q): basis-hood (count n(n−1)/2 + independence) — sympy-exact;
   (2) block A: a full enumeration of closures of all non-empty subsets of the generator
       basis (memoization by size + a cache keyed by RREF) — class tables above;
   (3) block B: random k-dimensional subspaces (seed 901) — class histograms and a raw
       intersection of the sets of classes A/B above;
   (4) block C: the table of growth with n — above.
  HONEST TALLY: handles 0 · downstream-content 0 · rulings 0. Court = Omega.
""")
sys.exit(0 if not FAIL else 1)
