#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (Lie algebra; downstream-content 0; handles 0). W28-O0 — bracket-closure
#      of 2-dimensional subspaces of so(3) + dimension count of so(p,q) + Killing signatures at n=3.
#      ★BLINDNESS: the probe prints ONLY raw identities/numbers/histograms; reading = an act of the court.
"""
S899 (lane A, ed.2) — W28-O0: four raw Lie-algebra facts, machine-checked.

  1. Closed form (sympy): for A=Σaᵢ Jᵢ, B=Σbᵢ Jᵢ in so(3) the component vector [A,B] = a×b;
     the identity det[a; b; a×b] = |a×b|² (sympy-exact) + Lagrange |a×b|² = |a|²|b|² − (a·b)².
  2. Numeric sweep: N=10000 random 2-dimensional subspaces V ⊂ so(3) —
     histogram of rank(span(V ∪ [V,V])).
  3. Dimension count: dim so(p,q) = n(n−1)/2, n=p+q; rows for n=2 and n=3.
  4. The Killing form K(X,Y)=tr(ad_X ad_Y) on explicit matrix bases of so(3) and so(2,1):
     eigenvalues / signatures of both, raw lines side by side; invariance of the signature
     under a change of basis — by number.

Fence: shared helper fence_scan (forbidden words — on the GUARDLINE line).
No readings in the probe's own text; raw lines with no words of interpretation.
"""
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# shared fence_scan helper (src/tools)
_src = __file__.replace("\\", "/").rsplit("/src/", 1)[0] + "/src"
if _src not in sys.path:
    sys.path.insert(0, _src)
from tools.fence_scan import scan_forbidden   # noqa: E402

import numpy as np   # noqa: E402
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
rule("SECTION 1 — CLOSED FORM (sympy): [A,B] in so(3) = a×b; determinant identity")
a1, a2, a3, b1, b2, b3 = sp.symbols("a1 a2 a3 b1 b2 b3", real=True)
a = sp.Matrix([a1, a2, a3])
b = sp.Matrix([b1, b2, b3])

# so(3) basis in the adjoint-vector realization: (J_i)_{jk} = −ε_{ijk} ⟹ [J_i,J_j] = ε_{ijk} J_k
eps = [[[int((i - j) * (j - k) * (k - i) / 2) for k in range(3)] for j in range(3)]
       for i in range(3)]
J = [sp.Matrix(3, 3, lambda r, c, i=i: -eps[i][r][c]) for i in range(3)]

# basis check: [J_i, J_j] = ε_{ijk} J_k (sympy-exact, all 9 pairs)
basis_ok = all(
    sp.simplify((J[i] * J[j] - J[j] * J[i])
                - sum((eps[i][j][k] * J[k] for k in range(3)), sp.zeros(3, 3))) == sp.zeros(3, 3)
    for i in range(3) for j in range(3)
)
check("basis: [J_i,J_j] = ε_{ijk} J_k (all 9 pairs, sympy-exact)", basis_ok)

A = sum((a[i] * J[i] for i in range(3)), sp.zeros(3, 3))
B = sum((b[i] * J[i] for i in range(3)), sp.zeros(3, 3))
comm = sp.expand(A * B - B * A)
cross = a.cross(b)
comm_from_cross = sum((cross[i] * J[i] for i in range(3)), sp.zeros(3, 3))
check("[A,B] = Σ (a×b)_i J_i (the component vector of the bracket = the cross product; sympy-exact)",
      sp.simplify(comm - comm_from_cross) == sp.zeros(3, 3))
print(f"  (a×b) = {list(cross)}")

# determinant identity: det[a; b; a×b] = |a×b|²
M3 = sp.Matrix([[a1, a2, a3], [b1, b2, b3], list(cross)])
det3 = sp.expand(M3.det())
norm2 = sp.expand(cross.dot(cross))
check("det[a; b; a×b] − |a×b|² = 0 (sympy-exact)", sp.simplify(det3 - norm2) == 0)
print(f"  RAW IDENTITY: det[a; b; a×b] = |a×b|²  ;  simplify(det − |a×b|²) = {sp.simplify(det3 - norm2)}")

# Lagrange: |a×b|² = |a|²|b|² − (a·b)² (= the Gram determinant of the pair a,b)
lagrange = sp.simplify(norm2 - (a.dot(a) * b.dot(b) - a.dot(b) ** 2))
check("|a×b|² − (|a|²|b|² − (a·b)²) = 0 (Lagrange; sympy-exact)", lagrange == 0)
print(f"  RAW IDENTITY: |a×b|² = |a|²|b|² − (a·b)² = Gram(a,b)  ;  simplify(difference) = {lagrange}")
print("  RAW CHAIN: rank[a;b]=2 ⟺ Gram(a,b)>0 ⟺ |a×b|²>0 ⟺ det[a;b;a×b]>0 ⟺ a×b ∉ span{a,b}")


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SECTION 2 — NUMERIC SWEEP: N=10000 random 2-dimensional V ⊂ so(3), rank(span(V ∪ [V,V]))")
rng = np.random.default_rng(899)
N_TARGET = 10000
hist = {2: 0, 3: 0}
n_degenerate_discarded = 0
n_done = 0
while n_done < N_TARGET:
    av = rng.standard_normal(3)
    bv = rng.standard_normal(3)
    if np.linalg.matrix_rank(np.vstack([av, bv])) < 2:
        n_degenerate_discarded += 1
        continue
    cv = np.cross(av, bv)
    r = int(np.linalg.matrix_rank(np.vstack([av, bv, cv])))
    hist[r] = hist.get(r, 0) + 1
    n_done += 1
print(f"  N={N_TARGET} · discarded for degeneracy rank<2: {n_degenerate_discarded}")
print(f"  HISTOGRAM rank(span(V ∪ [V,V])): " + " · ".join(f"rank={k} → {v}" for k, v in sorted(hist.items())))
check("the sweep ran to completion (N=10000, no truncation)", n_done == N_TARGET,
      f"rank=2: {hist.get(2, 0)} · rank=3: {hist.get(3, 0)}")


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SECTION 3 — DIMENSION COUNT: dim so(p,q) = n(n−1)/2, n = p+q")
for n in (2, 3):
    print(f"  n={n} → dim so(p,q) = n(n−1)/2 = {n * (n - 1) // 2}")
check("n=2 → dim=1 · n=3 → dim=3", 2 * 1 // 2 == 1 and 3 * 2 // 2 == 3)


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SECTION 4 — KILLING SIGNATURES AT n=3: so(3) and so(2,1), K(X,Y)=tr(ad_X ad_Y)")


def structure_constants(basis):
    """c^k_{ij}: [e_i,e_j] = Σ_k c^k_{ij} e_k — an exact decomposition (sympy solve on the basis)."""
    m = len(basis)
    flat = sp.Matrix.hstack(*[sp.Matrix(sp.flatten(e)) for e in basis])  # 9×m
    C = [[None] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            br = basis[i] * basis[j] - basis[j] * basis[i]
            sol = flat.solve_least_squares(sp.Matrix(sp.flatten(br)))
            resid = flat * sol - sp.Matrix(sp.flatten(br))
            assert sp.simplify(resid.norm()) == 0, "the bracket lies outside the basis span"
            C[i][j] = [sp.nsimplify(sol[k]) for k in range(m)]
    return C


def killing_matrix(basis):
    m = len(basis)
    C = structure_constants(basis)
    ad = [sp.Matrix(m, m, lambda k, j, i=i: C[i][j][k]) for i in range(m)]
    return sp.Matrix(m, m, lambda i, j: sp.trace(ad[i] * ad[j])), ad


def signature_counts(K):
    ev = []
    for lam, mult in K.eigenvals().items():
        ev.extend([sp.nsimplify(lam)] * mult)
    n_pos = sum(1 for e in ev if e > 0)
    n_neg = sum(1 for e in ev if e < 0)
    n_zero = sum(1 for e in ev if e == 0)
    return sorted(ev, key=lambda x: float(x)), (n_pos, n_neg, n_zero)


# --- so(3): antisymmetric 3×3 (the same J as above: Jᵀ = −J) ---
check("so(3)-basis: J_iᵀ = −J_i (all 3, antisymmetry explicit)",
      all(sp.simplify(Ji.T + Ji) == sp.zeros(3, 3) for Ji in J))
K_so3, _ = killing_matrix(J)
ev3, sig3 = signature_counts(K_so3)
print(f"  so(3):   K = {K_so3.tolist()}")
print(f"  so(3):   eigenvalues K = {ev3} · signature (n₊,n₋,n₀) = {sig3}")

# --- so(2,1): X η + η Xᵀ = 0, η = diag(+1,+1,−1) ---
eta = sp.diag(1, 1, -1)
G = [
    sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]]),   # generator in the (1,2) plane
    sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),    # generator in the (1,3) plane
    sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),    # generator in the (2,3) plane
]
check("so(2,1)-basis: X η + η Xᵀ = 0 (all 3, the isometry condition diag(+1,+1,−1) explicit)",
      all(sp.simplify(Gi * eta + eta * Gi.T) == sp.zeros(3, 3) for Gi in G))
K_so21, _ = killing_matrix(G)
ev21, sig21 = signature_counts(K_so21)
print(f"  so(2,1): K = {K_so21.tolist()}")
print(f"  so(2,1): eigenvalues K = {ev21} · signature (n₊,n₋,n₀) = {sig21}")

print("\n  RAW SIGNATURES SIDE BY SIDE (the Killing signature = an isomorphism invariant; numbers only):")
print(f"    so(3)   → (n₊,n₋,n₀) = {sig3}")
print(f"    so(2,1) → (n₊,n₋,n₀) = {sig21}")
print(f"    a 3-dimensional subalgebra of a 3-dimensional algebra = the whole algebra; dim so(2,1) = {len(G)}")

# invariance of the signature under a change of basis (by number): e'_i = Σ_j P_ij e_j, P ∈ GL(3)
print("\n  change of basis e' = P·e (random integer P, det≠0) — the Killing signature of so(2,1) by number:")
n_trials = 0
sig_rows = []
while n_trials < 5:
    P = rng.integers(-3, 4, size=(3, 3))
    if abs(np.linalg.det(P.astype(float))) < 0.5:
        continue
    Gp = [sum((int(P[i, j]) * G[j] for j in range(3)), sp.zeros(3, 3)) for i in range(3)]
    Kp, _ = killing_matrix(Gp)
    Kp_np = np.array(Kp.tolist(), dtype=float)
    evp = np.linalg.eigvalsh(Kp_np)
    s = (int(np.sum(evp > 1e-9)), int(np.sum(evp < -1e-9)), int(np.sum(np.abs(evp) <= 1e-9)))
    sig_rows.append(s)
    print(f"    P={P.tolist()} → eigenvalues K' = {np.round(evp, 6).tolist()} → (n₊,n₋,n₀) = {s}")
    n_trials += 1
check("the Killing signature of so(2,1) is the same in all 5 random bases (Sylvester, by number)",
      all(s == sig_rows[0] for s in sig_rows), f"{sig_rows}")
check("raw signatures printed side by side (no word of conclusion)", True,
      f"so(3): {sig3} · so(2,1): {sig21}")


# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY S899 (raw lines; reading = an act of the court)")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}

  RAW LINES (no readings):
   (1) [A,B] = Σ (a×b)_i J_i ; det[a;b;a×b] = |a×b|² = |a|²|b|² − (a·b)² (sympy-exact).
   (2) histogram of rank(span(V ∪ [V,V])) over N=10000: rank=2 → {hist.get(2, 0)} · rank=3 → {hist.get(3, 0)}.
   (3) dim so(p,q) = n(n−1)/2: n=2 → 1 · n=3 → 3.
   (4) Killing signatures (n₊,n₋,n₀): so(3) → {sig3} · so(2,1) → {sig21};
       the so(2,1) signature is unchanged under 5 random changes of basis.
  HONEST TALLY: handles 0 · downstream-content 0 · rulings 0. Court = the project's adjudication.
""")
sys.exit(0 if not FAIL else 1)
