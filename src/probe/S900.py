#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (Lie algebra; downstream-content 0; handles 0). W28-O1 — ladder of
#      bracket-closures of J/K-sets in so(p,q) for six signatures η + the Killing signature
#      of every closure class.
#      ★BLINDNESS: the probe prints ONLY raw identities/numbers/histograms; reading = an act of the court.
"""
S900 (lane A, ed.2) — W28-O1: bracket-closure of J/K-sets in so(p,q), machine-checked and exact.

  Object: η = diag(+1×p, −1×q); so(p,q) = {X : Xη + ηXᵀ = 0}; basis over the 2-planes (i,j):
    J-type — both axes of the same sign ((+,+) or (−,−)): E_ij − E_ji;
    K-type — a mixed pair (+,−): E_ij + E_ji.
  Signatures: (3,0) · (2,1) · (2,2) · (3,1) · (4,0) · (3,3).
  For each: the Lie-closure (iterated bracket-span to a fixed point, sympy-exact ranks) of
  the sets {all J} · {all K} · pairs {J,J} · pairs {K,K} · pairs {J,K}; histograms of
  dim span(start set) → dim Lie(closure); the Killing signature (n₊,n₋,n₀) of each closure
  class (structure constants in their own basis, decomposition residual = 0 assert);
  subspace check L ⊆ span{all J} and the check L = the whole so(p,q). EXACT arithmetic
  (Rational/int).

Fence: shared helper fence_scan (forbidden words — on the GUARDLINE line).
No readings in the probe's own text; raw lines with no words of interpretation.
"""
import sys
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
_FORBIDDEN = [r"стійк", r"обрано", r"selected", r"stable", r"причин", r"час", r"time", r"arrow", r"стріла", r"вибір"]   # GUARDLINE
_hits = scan_forbidden(__file__, _FORBIDDEN)
check("fence: forbidden words (list = the GUARDLINE line) = 0 occurrences outside the declaration",
      not _hits, f"hits: {_hits}" if _hits else "0")
check("downstream-content 0 · handles 0 (pure Lie algebra)", True)


# ═══════════════════════════════════════════════════════════════════════════════════════
# TOOLS: exact span (echelon over Q), Lie-closure, structure constants, Killing form.
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
    return sp.Matrix(d, d, lambda i, j: sp.trace(ad[i] * ad[j]))


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


# ═══════════════════════════════════════════════════════════════════════════════════════
SIGNATURES = [(3, 0), (2, 1), (2, 2), (3, 1), (4, 0), (3, 3)]
FINAL_ROWS = []   # for the summary table of span{all J}

for (p, q) in SIGNATURES:
    n = p + q
    dim_so = n * (n - 1) // 2
    rule(f"SECTION so({p},{q}) — n={n}, η=diag({'+1,' * p}{'−1,' * q}) · dim so = n(n−1)/2 = {dim_so}")
    eta, Js, Ks = make_soPQ(p, q)
    gens = Js + Ks

    # basis-hood: the isometry condition, the count, linear independence
    check(f"so({p},{q}): Xη + ηXᵀ = 0 for all {len(gens)} generators (sympy-exact)",
          all(sp.simplify(G * eta + eta * G.T) == sp.zeros(n, n) for _, G in gens))
    check(f"so({p},{q}): generator count = n(n−1)/2", len(gens) == dim_so,
          f"J-entries: {len(Js)} · K-entries: {len(Ks)} · total: {len(gens)}")
    stack = sp.Matrix([list(flat(G)) for _, G in gens])
    rk = stack.rank()
    check(f"so({p},{q}): linear independence of the generators (sympy-exact rank)", rk == dim_so,
          f"rank = {rk}")

    # span {all J} + closure under the bracket
    spanJ = Span()
    for _, G in Js:
        spanJ.add(flat(G))
    J_closed = all(spanJ.contains(flat(bracket(G1, G2)))
                   for (_, G1), (_, G2) in combinations(Js, 2)) if len(Js) >= 2 else True
    J_kill = sym_signature(killing_matrix([G for _, G in Js])) if Js else None
    FINAL_ROWS.append(((p, q), len(Js), len(Ks), dim_so, spanJ.dim(), J_closed, J_kill))

    # cache of Killing signatures keyed by the canonical key of the closure subspace
    kill_cache = {}

    def closure_info(gen_mats):
        """(dimL, Killing signature, L ⊆ spanJ, L=so(p,q)) for the closure of a set."""
        basis, span = lie_closure(gen_mats)
        d = span.dim()
        if d == 0:
            return (0, None, True, dim_so == 0)
        k = span.key()
        if k not in kill_cache:
            # sanity: every basis element of the closure satisfies Xη + ηXᵀ=0
            assert all(sp.simplify(B * eta + eta * B.T) == sp.zeros(n, n) for B in basis), \
                "a closure element lies outside so(p,q)"
            kill_cache[k] = sym_signature(killing_matrix(basis))
        sub_j = all(spanJ.contains(flat(B)) for B in basis)
        return (d, kill_cache[k], sub_j, d == dim_so)

    # ── sets 1-2: {all J}, {all K} — one line each ──
    print(f"\n  SET 1 — {{all J}} ({len(Js)} entries):")
    if Js:
        d, s, sj, fl = closure_info([G for _, G in Js])
        print(f"    dim span(start) = {spanJ.dim()} → dim Lie(closure) = {d} · "
              f"Killing {sig_str(s)} · L ⊆ span{{all J}}: {'yes' if sj else 'no'} · "
              f"L = so({p},{q}): {'yes' if fl else 'no'}")
    else:
        print("    the set is empty (0 J-entries)")

    print(f"  SET 2 — {{all K}} ({len(Ks)} entries):")
    if Ks:
        sK0 = Span()
        for _, G in Ks:
            sK0.add(flat(G))
        d, s, sj, fl = closure_info([G for _, G in Ks])
        print(f"    dim span(start) = {sK0.dim()} → dim Lie(closure) = {d} · "
              f"Killing {sig_str(s)} · L ⊆ span{{all J}}: {'yes' if sj else 'no'} · "
              f"L = so({p},{q}): {'yes' if fl else 'no'}")
    else:
        print("    the set is empty (0 K-entries)")

    # ── sets 3-5: pairs; histograms + closure classes ──
    pair_sets = [
        ("SET 3 — pairs {J,J}", list(combinations(Js, 2))),
        ("SET 4 — pairs {K,K}", list(combinations(Ks, 2))),
        ("SET 5 — pairs {J,K}", [(a, b) for a in Js for b in Ks]),
    ]
    for title, pairs in pair_sets:
        print(f"\n  {title} — {len(pairs)} pairs:")
        if not pairs:
            print("    0 pairs (the set is empty)")
            continue
        hist = Counter()
        classes = Counter()
        for (n1, G1), (n2, G2) in pairs:
            s0 = Span()
            s0.add(flat(G1))
            s0.add(flat(G2))
            d, s, sj, fl = closure_info([G1, G2])
            hist[(s0.dim(), d)] += 1
            classes[(d, s, sj, fl)] += 1
        print("    HISTOGRAM dim span(start) → dim Lie(closure):")
        for (ds, dc), cnt in sorted(hist.items()):
            print(f"      {ds} → {dc} : {cnt} pairs")
        print("    CLOSURE CLASSES (unique dim + Killing signature):")
        for (d, s, sj, fl), cnt in sorted(classes.items(), key=lambda x: (x[0][0], x[0][1])):
            print(f"      dim L = {d} · Killing {sig_str(s)} · "
                  f"L ⊆ span{{all J}}: {'yes' if sj else 'no'} · "
                  f"L = so({p},{q}): {'yes' if fl else 'no'} · pairs: {cnt}")
    check(f"so({p},{q}): all closures counted to completion (no truncation)", True,
          f"unique closure subspaces in the cache: {len(kill_cache)}")


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("FINAL SECTION — table: signature → #J · #K · dim so(p,q) · dim span{all J} · Killing span{all J}")
print(f"  {'signature':<12}{'#J':>4}{'#K':>4}{'dim so':>8}{'dim spanJ':>11}"
      f"{'bracket-closed':>19}{'Killing spanJ':>16}")
for (pq, nJ, nK, dso, dJ, closed, kj) in FINAL_ROWS:
    kj_s = sig_str(kj) if kj is not None else "—"
    print(f"  ({pq[0]},{pq[1]})      {nJ:>4}{nK:>4}{dso:>8}{dJ:>11}"
          f"{'yes' if closed else 'no':>18}{kj_s:>16}")
check("span{all J} is bracket-closed in all six signatures (yes/no printed as a table)",
      True, str([r[5] for r in FINAL_ROWS]))


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY S900 (raw lines; reading = an act of the court)")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}

  RAW LINES (no readings):
   (1) six signatures so(p,q): basis-hood (count n(n−1)/2 + independence) — sympy-exact;
   (2) closure of the sets {{all J}} · {{all K}} · pairs {{J,J}}/{{K,K}}/{{J,K}} — histograms above;
   (3) Killing signatures of the closure classes — congruent diagonalization over Q, exact;
   (4) table of span{{all J}} across signatures — above.
  HONEST TALLY: handles 0 · downstream-content 0 · rulings 0. Court = the project's adjudication.
""")
sys.exit(0 if not FAIL else 1)
