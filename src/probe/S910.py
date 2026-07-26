#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (Lie algebra over Q; handles 0). W28-O9 — automorphism-
#      equivariance of bracket-closures in so(p,q) for {(2,2)·(3,1)·(3,3)}: three families
#      of maps φ (Ad_S inner-isometric · sign-swap P · Cartan involution θ) —
#      closure(φV) =?= φ(closure(V)) on random V · orbits of subspaces within
#      classes {dim·Killing} · action of φ on the list of unique closure-subspaces.
#      ★BLINDNESS: the probe prints ONLY raw identities/counters/tables; reading = an act of the court.
"""
S910 (lane A, ed.2) — W28-O9: automorphism-equivariance of closures and orbits of
subspaces in so(p,q) for {(2,2)·(3,1)·(3,3)}.

  Three families of maps φ on so(η) (η = diag(+1×p, −1×q)):
   1. Ad_S (inner-isometric): X ↦ S⁻¹XS for rationally-exact S with
      SᵀηS = η, both types det S = ±1 (rotations/boosts/reflections/composites —
      constructions from S906); 4 examples per type per signature.
   2. Sign-swap: P — a permutation that swaps all +axes with −axes (only
      p = q; for (3,1) stamped as a raw line, why it does not exist). Explicit checks:
      PᵀηP = −η · P⁻¹XP ∈ so(η) for all basis X · so(η) = so(−η) by a symbolic line
      (Xη+ηXᵀ=0 ⟺ X(−η)+(−η)Xᵀ=0) · whether it preserves the bracket (the auto-property
      [φX,φY] = φ[X,Y] — exact basis-by-basis).
   3. Cartan involution θ(X) = −Xᵀ: θ(so(η)) ⊆ so(η) basis-by-basis ·
      [θX,θY] = θ[X,Y] symbolically (generic 4×4) and basis-by-basis.

  Steps:
   O9a — equivariance of closure: closure(φV) =?= φ(closure(V)) for 10
     random 2-3-dimensional V (integer coeffs {−3..3}, seed 910) × every φ of the
     three families — bit-comparison of subspaces (RREF keys); PASS/FAIL counters.
   O9b — orbits within classes: the (3,3) class {dim 3 · Killing (0,3,0)} has 2
     unique subspaces (the J-sector of the +block span{J(0,1),J(0,2),J(1,2)} and
     the J-sector of the −block span{J(3,4),J(3,5),J(4,5)}) — whether the sign-swap P
     maps one to the other (bit-check). Also pairs from the "no"-lines of O7c (S906): on (3,3)
     X = η·J(0,1) vs Y = η·J(3,4); on (2,2) η·J(0,1) vs η·J(2,3) — whether the sign-swap
     or the composite {swap ∘ isometry} connects them; raw "yes/no + by what".
   O9c — action of every φ on the FULL list of unique closure-subspaces of
     coordinate subsets ((2,2) and (3,1) — all 63 subsets; (3,3) —
     subsets of size ≤2 + the whole sectors {all J}/{all K}, the bound is stamped):
     whether φ permutes the list (a bijection of keys) and whether it preserves the class
     keys {dim·Killing} — the class count before/after.
   O9d — rollup: which counters are equal before/after under all three families
     (the class-level ones), which change only in the split into concrete subspaces
     (the presentational ones) — raw numbers, no readings.

  EXACT arithmetic: sympy Rational/Integer over Q; no tolerances. Mechanisms
  (Span/RREF-key/lie_closure/structure-constants/Killing/signature — a verbatim
  copy from S900/S901; rot/boost/refl/perm2 — a copy from S906). Bounds are stamped.

Fence: the shared fence_scan helper (forbidden words — on the GUARDLINE line;
case-insensitive detection by default — S909). No verdicts in the text.
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

# ── tee: all of stdout is duplicated into S910_run.log next to the script ──
_LOG_PATH = __file__.replace("\\", "/").rsplit("/", 1)[0] + "/S910_run.log"


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
rule("SECTION 0 — FENCE (shared fence_scan helper; case-insensitive — S909)")
_FORBIDDEN = [r"стійк", r"обрано", r"selected", r"stable", r"причин", r"час", r"time", r"arrow", r"стріла", r"вибір", r"краще", r"злам", r"поле", r"матерія", r"хіральн", r"зберігається", r"узагальнююч"]   # GUARDLINE
_hits = scan_forbidden(__file__, _FORBIDDEN)
check("fence: forbidden words (list = GUARDLINE line) = 0 occurrences outside the declaration",
      not _hits, f"hits: {_hits}" if _hits else "0")
check("handles 0 (pure Lie algebra: Ad_S · sign-swap · Cartan involution · closure)", True)


# ═══════════════════════════════════════════════════════════════════════════════════════
# TOOLS (verbatim copy from S900/S901): exact span (echelon over Q), Lie closure,
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


def span_of(mats):
    s = Span()
    for M in mats:
        s.add(flat(M))
    return s


# ── isometry builders (verbatim copy from S906) ──

def rot(nn, i, j, c, s):
    """Block rotation in a same-sign plane (i,j): c² + s² = 1 (rational)."""
    R = sp.eye(nn)
    R[i, i] = c
    R[j, j] = c
    R[i, j] = -s
    R[j, i] = s
    return R


def boost(nn, i, j, ch, sh):
    """Block boost in a mixed-sign plane (i,j): ch² − sh² = 1 (rational)."""
    R = sp.eye(nn)
    R[i, i] = ch
    R[j, j] = ch
    R[i, j] = sh
    R[j, i] = sh
    return R


def refl(nn, i):
    """Reflection of a single axis i: diag(1,…,−1,…,1)."""
    R = sp.eye(nn)
    R[i, i] = -1
    return R


Q = sp.Rational
SIGS = [(2, 2), (3, 1), (3, 3)]
SEARCH_STAMPS = []

ISO_PLUS = {
    (2, 2): [
        ("rotation(0,1; 3/5,4/5)",    lambda nn: rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
        ("rotation(2,3; 5/13,12/13)", lambda nn: rot(nn, 2, 3, Q(5, 13), Q(12, 13))),
        ("boost(0,2; 5/3,4/3)",       lambda nn: boost(nn, 0, 2, Q(5, 3), Q(4, 3))),
        ("boost(1,3; 13/5,12/5)",     lambda nn: boost(nn, 1, 3, Q(13, 5), Q(12, 5))),
    ],
    (3, 1): [
        ("rotation(0,1; 3/5,4/5)",    lambda nn: rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
        ("rotation(1,2; 5/13,12/13)", lambda nn: rot(nn, 1, 2, Q(5, 13), Q(12, 13))),
        ("boost(0,3; 5/3,4/3)",       lambda nn: boost(nn, 0, 3, Q(5, 3), Q(4, 3))),
        ("boost(2,3; 13/5,12/5)",     lambda nn: boost(nn, 2, 3, Q(13, 5), Q(12, 5))),
    ],
    (3, 3): [
        ("rotation(0,1; 3/5,4/5)",    lambda nn: rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
        ("rotation(3,4; 5/13,12/13)", lambda nn: rot(nn, 3, 4, Q(5, 13), Q(12, 13))),
        ("boost(0,3; 5/3,4/3)",       lambda nn: boost(nn, 0, 3, Q(5, 3), Q(4, 3))),
        ("boost(2,5; 13/5,12/5)",     lambda nn: boost(nn, 2, 5, Q(13, 5), Q(12, 5))),
    ],
}
ISO_MINUS = {
    (2, 2): [
        ("reflection of axis 0",       lambda nn: refl(nn, 0)),
        ("reflection of axis 2",       lambda nn: refl(nn, 2)),
        ("reflection of axis 0 · boost(0,2;5/3,4/3)",
         lambda nn: refl(nn, 0) * boost(nn, 0, 2, Q(5, 3), Q(4, 3))),
        ("reflection of axis 2 · rotation(0,1;3/5,4/5)",
         lambda nn: refl(nn, 2) * rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
    ],
    (3, 1): [
        ("reflection of axis 0",       lambda nn: refl(nn, 0)),
        ("reflection of axis 3",       lambda nn: refl(nn, 3)),
        ("reflection of axis 2 · rotation(0,1;3/5,4/5)",
         lambda nn: refl(nn, 2) * rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
        ("reflection of axis 1 · boost(0,3;5/3,4/3)",
         lambda nn: refl(nn, 1) * boost(nn, 0, 3, Q(5, 3), Q(4, 3))),
    ],
    (3, 3): [
        ("reflection of axis 0",       lambda nn: refl(nn, 0)),
        ("reflection of axis 3",       lambda nn: refl(nn, 3)),
        ("reflection of axis 1 · rotation(0,2;3/5,4/5)",
         lambda nn: refl(nn, 1) * rot(nn, 0, 2, Q(3, 5), Q(4, 5))),
        ("reflection of axis 4 · boost(1,4;5/3,4/3)",
         lambda nn: refl(nn, 4) * boost(nn, 1, 4, Q(5, 3), Q(4, 3))),
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SECTION 1 — THE THREE FAMILIES φ: symbolic lemmas (generic) + construction on every signature")

# swap lemma: so(η) = so(−η) — a symbolic line on generic X (4×4) and symbolic diagonal η
_Xg = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"x_{i}_{j}"))
_Yg = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"y_{i}_{j}"))
_Hg = sp.diag(*[sp.Symbol(f"h_{i}") for i in range(4)])
print("    symbolic line (X generic 4×4, η = diag(h₀..h₃) symbolic):")
print("      X(−η) + (−η)Xᵀ = −(Xη + ηXᵀ)  ⟹  Xη+ηXᵀ=0 ⟺ X(−η)+(−η)Xᵀ=0  ⟹  so(η) = so(−η)")
check("swap-lemma: expand( (X(−η)+(−η)Xᵀ) + (Xη+ηXᵀ) ) = 0 — symbolic, generic",
      sp.expand((_Xg * (-_Hg) + (-_Hg) * _Xg.T) + (_Xg * _Hg + _Hg * _Xg.T)) == sp.zeros(4, 4))

# θ-lemma: [θX,θY] = θ[X,Y] on generic X,Y (4×4): [−Xᵀ,−Yᵀ] = −[X,Y]ᵀ
print("    symbolic line (X,Y generic 4×4): [−Xᵀ,−Yᵀ] = XᵀYᵀ−YᵀXᵀ = (YX−XY)ᵀ = −[X,Y]ᵀ = θ[X,Y]")
check("θ-lemma: expand( [θX,θY] − θ[X,Y] ) = 0 — symbolic, generic (32 symbols)",
      sp.expand(bracket(-_Xg.T, -_Yg.T) - (-(bracket(_Xg, _Yg)).T)) == sp.zeros(4, 4))

CTX = {}   # signature → {eta, Js, Ks, gens, G, phis, P}

for (p, q) in SIGS:
    nn = p + q
    dim_so = nn * (nn - 1) // 2
    print(f"\n  so({p},{q}) — n = {nn} · dim so = {dim_so}:")
    eta, Js, Ks = make_soPQ(p, q)
    gens = Js + Ks
    G = [M for _, M in gens]

    # basis property (copy of the discipline from S900/S901)
    check(f"so({p},{q}): Xη + ηXᵀ = 0 for all {len(gens)} generators (sympy-exact)",
          all((Gm * eta + eta * Gm.T) == sp.zeros(nn, nn) for Gm in G))
    check(f"so({p},{q}): generator count = n(n−1)/2 and linear independence",
          len(gens) == dim_so and sp.Matrix([list(flat(Gm)) for Gm in G]).rank() == dim_so,
          f"J-moves: {len(Js)} · K-moves: {len(Ks)}")

    bb = {}   # base brackets (for the equivariance checks)
    for i in range(dim_so):
        for j in range(i + 1, dim_so):
            bb[(i, j)] = bracket(G[i], G[j])

    phis = []   # (family, label, φ)

    # ── family 1: Ad_S, X ↦ S⁻¹XS, SᵀηS = η, det S = ±1, 4 per type ──
    for det_want, iso_list, word in ((1, ISO_PLUS[(p, q)], "det=+1"),
                                     (-1, ISO_MINUS[(p, q)], "det=−1")):
        for label, build in iso_list:
            S = build(nn)
            assert S.T * eta * S == eta, f"S is not an isometry: {label}"
            assert S.det() == det_want, f"det S ≠ {det_want}: {label}"
            Sinv = eta * S.T * eta   # S⁻¹ = ηSᵀη (from SᵀηS = η, η² = 𝟙) — exact
            assert Sinv * S == sp.eye(nn), f"S⁻¹ ≠ ηSᵀη: {label}"
            phis.append(("Ad", f"Ad[{label}] {word}",
                         (lambda X, A=Sinv, B=S: A * X * B)))

    # ── family 2: sign-swap P (only p = q) ──
    P = None
    if p == q:
        P = sp.zeros(nn, nn)
        for i in range(p):
            P[i, p + i] = 1
            P[p + i, i] = 1
        assert P * P == sp.eye(nn), "P is not an involution"
        cyc = " · ".join(f"{i}↔{p + i}" for i in range(p))
        print(f"    sign-swap P: {cyc} · det P = {P.det()} · Pᵀ = P = P⁻¹")
        check(f"so({p},{q}): PᵀηP = −η (explicit)", P.T * eta * P == -eta)
        phis.append(("swap", "swap P", (lambda X, A=P: A * X * A)))
    else:
        print(f"    sign-swap P for ({p},{q}): does NOT exist — p ≠ q: a permutation of axes would have to")
        print(f"      hit the set of +axes (count {p}) bijectively onto the set of −axes (count {q});")
        print(f"      {p} ≠ {q} ⟹ no bijection exists, PᵀηP = −η is unreachable by a permutation.")
        check(f"so({p},{q}): the stamp on the non-existence of the swap is printed (p ≠ q)", p != q)

    # ── family 3: Cartan involution θ ──
    phis.append(("θ", "θ(X) = −Xᵀ", (lambda X: -X.T)))

    # ── checks of all φ: membership of images in so(η) + bracket-equivariance basis-by-basis ──
    n_mem_all = n_br_all = 0
    for fam, label, phi in phis:
        phiG = [phi(Gm) for Gm in G]
        n_mem = sum(1 for Y in phiG if (Y * eta + eta * Y.T) == sp.zeros(nn, nn))
        n_br = sum(1 for (i, j), B in bb.items()
                   if bracket(phiG[i], phiG[j]) == phi(B))
        print(f"    {label:<46} membership φ(basis) in so(η): {n_mem}/{dim_so} · "
              f"[φX,φY]=φ[X,Y]: {n_br}/{len(bb)}")
        n_mem_all += (n_mem == dim_so)
        n_br_all += (n_br == len(bb))
    check(f"so({p},{q}): all {len(phis)} φ — membership of images in so(η) basis-by-basis complete",
          n_mem_all == len(phis), f"{n_mem_all}/{len(phis)}")
    check(f"so({p},{q}): all {len(phis)} φ — bracket-equivariance basis-by-basis complete",
          n_br_all == len(phis), f"{n_br_all}/{len(phis)}")

    CTX[(p, q)] = dict(eta=eta, Js=Js, Ks=Ks, gens=gens, G=G, phis=phis, P=P,
                       dim_so=dim_so, nn=nn)


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O9a — EQUIVARIANCE OF CLOSURE: closure(φV) =?= φ(closure(V)) — RREF bit-comparison")

for (p, q) in SIGS:
    ctx = CTX[(p, q)]
    nn, dim_so, G, phis = ctx["nn"], ctx["dim_so"], ctx["G"], ctx["phis"]
    stamp = (f"O9a so({p},{q}): BOUNDARY (stamp): 10 random V (5 × dim 2 + 5 × dim 3), "
             f"integer coeffs {{−3..3}}, seed 910 (Random(910) per signature), rank<k discarded by retry")
    print(f"\n  so({p},{q}): {stamp}")
    SEARCH_STAMPS.append(stamp)
    rng = random.Random(910)
    per_phi = Counter()
    n_total = 0
    for idx in range(10):
        k = 2 if idx < 5 else 3
        tries = 0
        while True:
            tries += 1
            assert tries <= 1000, "discarding rank<k does not converge"
            mats = []
            for _ in range(k):
                coeffs = [rng.randint(-3, 3) for _ in range(dim_so)]
                Mx = sp.zeros(nn, nn)
                for c, Gm in zip(coeffs, G):
                    if c:
                        Mx = Mx + c * Gm
                mats.append(Mx)
            if span_of(mats).dim() == k:
                break
        cb, cs = lie_closure(mats)
        ckey = cs.key()
        hits = 0
        for fam, label, phi in phis:
            _, s_left = lie_closure([phi(M) for M in mats])     # closure(φV)
            s_right = span_of([phi(B) for B in cb])             # φ(closure(V))
            ok = (s_left.key() == s_right.key())
            hits += ok
            per_phi[label] += ok
            n_total += 1
        print(f"    V#{idx + 1:02d} k={k} (draws until rank=k: {tries}) · dim closure(V) = {cs.dim()} "
              f"· matches closure(φV) = φ(closure(V)): {hits}/{len(phis)}")
    print(f"    COUNTER by φ (PASS out of 10 V):")
    for fam, label, phi in phis:
        print(f"      {label:<46} {per_phi[label]}/10")
    check(f"O9a so({p},{q}): closure(φV) = φ(closure(V)) on all {n_total} comparisons "
          f"(10 V × {len(phis)} φ)", sum(per_phi.values()) == n_total,
          f"{sum(per_phi.values())}/{n_total}")


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O9b — ORBITS WITHIN CLASSES: the (3,3) class {dim 3 · Killing (0,3,0)} + O7c \"no\"-lines")

ctx33 = CTX[(3, 3)]
by_name33 = dict(ctx33["gens"])
eta33 = ctx33["eta"]
Vp_names = ["J(0,1)", "J(0,2)", "J(1,2)"]
Vm_names = ["J(3,4)", "J(3,5)", "J(4,5)"]
Vp = [by_name33[nm] for nm in Vp_names]
Vm = [by_name33[nm] for nm in Vm_names]
print(f"\n  (3,3): two unique subspaces of the class {{dim 3 · Killing (0,3,0)}} (from the S901 list):")
for tag, names, mats in (("V⁺ (J-sector of the +block)", Vp_names, Vp),
                         ("V⁻ (J-sector of the −block)", Vm_names, Vm)):
    b, s = lie_closure(mats)
    s0 = span_of(mats)
    ksig = sym_signature(killing_matrix(b))
    print(f"    {tag} = span{{{', '.join(names)}}}: dim span = {s0.dim()} · "
          f"dim closure = {s.dim()} · closure = span: {'yes' if s.key() == s0.key() else 'no'} · "
          f"Killing {sig_str(ksig)}")
    check(f"O9b (3,3) {tag}: closure-fixedness (dim 3) and Killing = (0,3,0)",
          s.dim() == 3 and s.key() == s0.key() and ksig == (0, 3, 0))

keyVp = span_of(Vp).key()
keyVm = span_of(Vm).key()
P33 = ctx33["P"]
imgVp = span_of([P33 * M * P33 for M in Vp]).key()
imgVm = span_of([P33 * M * P33 for M in Vm]).key()
check("O9b (3,3): the sign-swap P maps V⁺ → V⁻ (bit-check of RREF keys)", imgVp == keyVm)
check("O9b (3,3): the sign-swap P maps V⁻ → V⁺ (bit-check of RREF keys)", imgVm == keyVp)
th_img = span_of([-M.T for M in Vp]).key()
print(f"    θ on V⁺: img = V⁺? {'yes' if th_img == keyVp else 'no'} · "
      f"img = V⁻? {'yes' if th_img == keyVm else 'no'}")
n_pp = n_pm = n_oth = 0
for fam, label, phi in ctx33["phis"]:
    if fam != "Ad":
        continue
    ik = span_of([phi(M) for M in Vp]).key()
    if ik == keyVp:
        n_pp += 1
    elif ik == keyVm:
        n_pm += 1
    else:
        n_oth += 1
print(f"    8 Ad_S on V⁺: img = V⁺: {n_pp} · img = V⁻: {n_pm} · another subspace: {n_oth}")

# O7c (S906) "no"-lines: whether the sign-swap or the composite {swap ∘ isometry} connects them
print("\n  O7c (S906) \"no\"-lines — the sign-swap P and the composite T = P·reflection:")
for (p, q), xg, yg, ax in [((3, 3), "J(0,1)", "J(3,4)", 3),
                           ((2, 2), "J(0,1)", "J(2,3)", 2)]:
    ctx = CTX[(p, q)]
    nn = ctx["nn"]
    eta = ctx["eta"]
    by_name = dict(ctx["gens"])
    P = ctx["P"]
    Xg, Yg = by_name[xg], by_name[yg]
    OmX, OmY = eta * Xg, eta * Yg
    assert OmX.T == -OmX and OmY.T == -OmY, "ηX/ηY not antisymmetric"
    # X-level (algebra generators)
    hitX = (P * Xg * P == Yg)
    print(f"    so({p},{q}) X-level: P⁻¹·{xg}·P = {yg}? {'yes' if hitX else 'no'} — by what: swap P")
    # Ω-level (O7c objects: Ω = η·generator, congruence TᵀΩT)
    c_pos = (P.T * OmX * P == OmY)
    c_neg = (P.T * OmX * P == -OmY)
    print(f"    so({p},{q}) Ω-level: Pᵀ(η·{xg})P = η·{yg}? {'yes' if c_pos else 'no'} · "
          f"= −η·{yg}? {'yes' if c_neg else 'no'}")
    T = P * refl(nn, ax)
    hitT = (T.T * OmX * T == OmY)
    anti = (T.T * eta * T == -eta)
    print(f"    so({p},{q}) composite T = P·reflection(axis {ax}): Tᵀ(η·{xg})T = η·{yg}? "
          f"{'yes' if hitT else 'no'} — by what: T; det T = {T.det()} · TᵀηT = −η: "
          f"{'yes' if anti else 'no'} (T is not an η-isometry, but an (η→−η)-congruence)")
    check(f"O9b so({p},{q}) O7c \"no\"-line: the connection is printed as raw lines "
          f"(X-level: swap; Ω-level: swap∘reflection)", hitX and hitT)


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O9c — ACTION OF φ ON THE LIST OF UNIQUE CLOSURE-SUBSPACES OF COORDINATE SUBSETS")


def enum_closures(G, eta, nloc, size_max, extras):
    """List of closures of subsets of the generator basis (memoized by size,
    cached by RREF key — the S901 scheme; sizes bounded by size_max; extras —
    additional index sets, closed directly)."""
    key_info = {}
    key_basis = {(): []}

    def register(basis, span):
        ck = span.key()
        if ck not in key_info:
            assert all((B * eta + eta * B.T) == sp.zeros(nloc, nloc) for B in basis), \
                "an element of the closure is outside so(p,q)"
            for a in range(len(basis)):
                for b in range(a + 1, len(basis)):
                    assert span.contains(flat(bracket(basis[a], basis[b]))), \
                        "the closure did not reach a fixed point"
            key_info[ck] = (span.dim(), sym_signature(killing_matrix(basis)))
            key_basis[ck] = basis
        return ck

    subkey = {(): ()}
    extend_cache = {}
    n_calls = 0
    for size in range(1, size_max + 1):
        for Sub in combinations(range(len(G)), size):
            parent, x = Sub[:-1], Sub[-1]
            pk = subkey[parent]
            ck = extend_cache.get((pk, x))
            if ck is None:
                basis, span = lie_closure(list(key_basis[pk]) + [G[x]])
                n_calls += 1
                ck = register(basis, span)
                extend_cache[(pk, x)] = ck
            subkey[Sub] = ck
    extra_keys = {}
    for label, idxs in extras:
        basis, span = lie_closure([G[i] for i in idxs])
        n_calls += 1
        extra_keys[label] = register(basis, span)
    return key_info, key_basis, subkey, extra_keys, n_calls


RES = []          # (signature, family, label, #keys, bijection, #new, #moved,
                  #  class counts equal, per-key class equality)
MOVED_MAPS = {}   # (signature, label φ) → img map (swap/θ only, for the class sub-table)

for (p, q) in SIGS:
    ctx = CTX[(p, q)]
    nn, dim_so, G, phis = ctx["nn"], ctx["dim_so"], ctx["G"], ctx["phis"]
    eta = ctx["eta"]
    nJ = len(ctx["Js"])
    if (p, q) != (3, 3):
        size_max = dim_so
        extras = []
        stamp = (f"O9c so({p},{q}): full enumeration — all 2^{dim_so}−1 = "
                 f"{2 ** dim_so - 1} nonempty coordinate subsets")
    else:
        size_max = 2
        extras = [("sector {all J} (6 gen)", list(range(nJ))),
                  ("sector {all K} (9 gen)", list(range(nJ, dim_so)))]
        stamp = ("O9c so(3,3): BOUNDARY (stamp): subsets of size ≤2 (15 + 105 = 120 starts) "
                 "+ whole sectors {all J}/{all K}; the full 2^15−1 = 32767 is not run here")
    print(f"\n  so({p},{q}): {stamp}")
    SEARCH_STAMPS.append(stamp)

    key_info, key_basis, subkey, extra_keys, n_calls = enum_closures(G, eta, nn, size_max, extras)
    n_starts = sum(1 for Sub in subkey if Sub)
    exp_starts = sum(sp.binomial(dim_so, s) for s in range(1, size_max + 1))
    check(f"O9c so({p},{q}): {n_starts} subset-starts processed (expected {exp_starts}) "
          f"+ {len(extra_keys)} sector sets", n_starts == exp_starts,
          f"lie_closure calls (memo) = {n_calls}")

    uniq = set(ck for Sub, ck in subkey.items() if Sub) | set(extra_keys.values())
    classK = {}
    for ck in uniq:
        classK.setdefault(key_info[ck], set()).add(ck)
    print(f"    unique closure-subspaces: {len(uniq)} · classes {{dim·Killing}}: {len(classK)}")
    print(f"    {'dim L':>6} {'Killing':>12} {'uniq.subspaces':>19}")
    for cl in sorted(classK):
        d, s = cl
        print(f"    {d:>6} {sig_str(s):>12} {len(classK[cl]):>19}")

    uniq_list = sorted(uniq, key=lambda k: (len(k), str(k)))
    classes_before = Counter(key_info[k] for k in uniq_list)
    img_cache = {}
    print(f"\n    action of φ on the list ({len(uniq_list)} keys):")
    for fam, label, phi in phis:
        img_list = []
        n_fixed = n_new = n_cls_eq = 0
        classes_after = Counter()
        img_map = {}
        for ck in uniq_list:
            mapped = [phi(B) for B in key_basis[ck]]
            ik = span_of(mapped).key()
            img_list.append(ik)
            img_map[ck] = ik
            if ik == ck:
                n_fixed += 1
            if ik in key_info:
                info = key_info[ik]
            elif ik in img_cache:
                info = img_cache[ik]
            else:
                info = (len(ik), sym_signature(killing_matrix(mapped)))
                img_cache[ik] = info
            if ik not in key_info:
                n_new += 1
            if info == key_info[ck]:
                n_cls_eq += 1
            classes_after[info] += 1
        n_moved = len(uniq_list) - n_fixed
        bij = (set(img_list) == uniq and len(set(img_list)) == len(uniq_list))
        cls_eq = (classes_after == classes_before)
        print(f"      {label:<46} fixed {n_fixed:>3} · moved {n_moved:>3} · "
              f"new outside the list {n_new:>3} · bijection of the list: {'yes' if bij else 'no'} · "
              f"class counts before=after: {'yes' if cls_eq else 'no'} · "
              f"class(img)=class(source): {n_cls_eq}/{len(uniq_list)}")
        RES.append(((p, q), fam, label, len(uniq_list), bij, n_new, n_moved, cls_eq, n_cls_eq))
        if fam in ("swap", "θ"):
            MOVED_MAPS[((p, q), label)] = (img_map, key_info, classK)

    # sub-table: classes with >1 unique subspace — movement of keys under swap/θ
    multi = [cl for cl in sorted(classK) if len(classK[cl]) > 1]
    if multi:
        print(f"\n    classes with >1 unique subspace ({len(multi)}): movement of keys under swap/θ:")
        for ((sg, lb), (img_map, ki, cK)) in MOVED_MAPS.items():
            if sg != (p, q):
                continue
            for cl in multi:
                keys_cl = cK[cl]
                mv = sum(1 for k in keys_cl if img_map[k] != k)
                inside = sum(1 for k in keys_cl if img_map[k] in keys_cl)
                print(f"      {lb:<10} class [dim {cl[0]} · Killing {sig_str(cl[1])}]: "
                      f"subspaces {len(keys_cl)} · moved {mv} · "
                      f"image in the same class: {inside}/{len(keys_cl)}")


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O9d — ROLLUP: class-level counters vs the split into subspaces (raw numbers)")

print(f"\n  {'signature':<11}{'family':<7}{'#φ':>4}{'bijection yes':>13}{'class counts yes':>17}"
      f"{'new (min..max)':>19}{'moved (min..max)':>25}")
FAM_ORDER = ["Ad", "swap", "θ"]
n_cls_rows_eq = 0
n_rows_total = 0
for (p, q) in SIGS:
    for fam in FAM_ORDER:
        rows = [r for r in RES if r[0] == (p, q) and r[1] == fam]
        if not rows:
            print(f"  ({p},{q})     {fam:<7}{'—':>4}{'—':>13}{'—':>17}{'—':>19}{'—':>25}")
            continue
        nb = sum(1 for r in rows if r[4])
        nc = sum(1 for r in rows if r[7])
        news = [r[5] for r in rows]
        movs = [r[6] for r in rows]
        n_cls_rows_eq += nc
        n_rows_total += len(rows)
        print(f"  ({p},{q})     {fam:<7}{len(rows):>4}{f'{nb}/{len(rows)}':>13}"
              f"{f'{nc}/{len(rows)}':>17}{f'{min(news)}..{max(news)}':>19}"
              f"{f'{min(movs)}..{max(movs)}':>25}")

print(f"""
  RAW ROLLUP LINES (counters, no readings):
   · class-level counters (#classes {{dim·Killing}} · subspace count per class ·
     class(img)=class(source)): equal before/after in {n_cls_rows_eq}/{n_rows_total} φ-rows
     of all three families on all three signatures — numbers in the O9c tables above;
   · split into concrete subspaces (RREF keys of the list): under θ, moved is 0
     on all signatures; under swap P — permuted by a bijection (moved counts in O9c);
     under Ad_S — the set of keys differs (new outside the list > 0 in O9c rows) while
     the class counters remain equal. Raw min..max — in the table above.""")
_bad_hdr = check("O9d: the rollup table is printed for all signatures × families", True)


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY S910 (raw lines; reading = an act of the court)")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}

  BOUND STAMPS:""")
for st in SEARCH_STAMPS:
    print(f"   · {st}")
print(f"""
  RAW LINES (no readings):
   (1) three families φ: Ad_S (4 examples per type det=±1, S⁻¹ = ηSᵀη exact) ·
       sign-swap P (PᵀηP = −η explicit; for (3,1) a stamp of non-existence: p ≠ q) ·
       Cartan involution θ(X) = −Xᵀ; membership of images in so(η) and
       [φX,φY] = φ[X,Y] — basis-by-basis exact + symbolic generic 4×4 lemmas;
   (2) O9a: closure(φV) =?= φ(closure(V)) — 10 random V × every φ ×
       3 signatures, RREF bit-comparison — counters above;
   (3) O9b: swap-action on the pair of subspaces of the class {{dim 3 · Killing (0,3,0)}} on
       (3,3) — bit-checks; O7c "no"-lines — connection by swap/swap∘reflection with
       TᵀηT = −η as a raw line;
   (4) O9c: action of φ on the list of unique closure-subspaces — fixed/
       moved/new/bijection/class counters — tables above;
   (5) O9d: rollup table class-level vs presentational counters — above.
  HONEST TALLY: handles 0 · verdicts 0. Court = Omega.
""")
_logf.flush()
sys.exit(0 if not FAIL else 1)
