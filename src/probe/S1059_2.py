# -*- coding: utf-8 -*-
# DIM: exact-integer (characters) + numeric-grid (components).  Every number carries a bracket
#      [address · unit · type/operation]; numeric numbers CARRY THEIR COVERAGE on the same line.
"""S1059-B — THE ∀d-LAW OF THE MEDIAL: does the choice of the Sym²/Λ² half follow from an ANCESTOR.

ASSIGNMENT: step 4 of `UNITS_LADDER_EXANTE.md` (the condition «only if 1-3 are alive» is met: 1-2 = S1054/S1055,
3 = S1059 the first half).  Exante: `hub/prime/S1059_MEDIAL_FORALL_D_LAW_EXANTE.md`.

THE EDITION OF THE CARRIERS: S1052 (the morphism q_P is native; the ∀d-claim is dead) · S1050/J-0503 (Test-3: «d nodes =
d channels» is a homonym, the verdict stands) · V4_BRIDGE_MEDIAL_EXANTE.

THE OBJECT (a bond-index, NOT a real-space lattice — the S1052 naming discipline):
  u ∈ (ℝ/2πℤ)^{d+1} with Σu_j ≡ 0 (this is a d-torus);  z_j = e^{i u_j};  the node-set N_d = {Σ_j z_j = 0}.
  Both conditions of Test-3 are met BY CONSTRUCTION: Σz=0 explicitly, Πz=1 automatically (because Σu=0 ⟸ Σδ_j=0).
  S_{d+1} acts by PERMUTING u — this is legitimate, because it preserves Σu=0.  No origin-choice, no choice at all.

ORDER (cheap before expensive): KILL-FIRST (the embedding in std⊗std) → only then the law.
"""
import os
import sys
from itertools import permutations

import numpy as np
import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
_fails, _passes = [], []


def ok(cond, msg):
    (_passes if cond else _fails).append(msg)
    print(("  ✓ " if cond else "  ✗ FAIL ") + msg)


class Tee:
    def __init__(self, real, fh):
        self.real, self.fh, self.chunks = real, fh, []

    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


# ══════════════ PART A: EXACT CHARACTERS OF S_n (Murnaghan–Nakayama, integers) ══════════════
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def _beta(lam, length):
    lam = list(lam) + [0] * (length - len(lam))
    return [lam[i] + (length - 1 - i) for i in range(length)]


def _from_beta(bs):
    bs = sorted(bs, reverse=True)
    lam = [bs[i] - (len(bs) - 1 - i) for i in range(len(bs))]
    return tuple(x for x in lam if x > 0)


def mn_char(lam, rho):
    """χ^λ(ρ) via Murnaghan–Nakayama through β-numbers. Exact integers."""
    if sum(lam) == 0:
        return 1
    if not rho:
        return 1 if sum(lam) == 0 else 0
    k, rest = rho[0], rho[1:]
    L = sum(lam) + len(lam) + 2
    bs = _beta(lam, L)
    bset = set(bs)
    total = 0
    for b in bs:
        if b >= k and (b - k) not in bset:
            height = sum(1 for x in bs if b - k < x < b)
            newb = (bset - {b}) | {b - k}
            total += (-1) ** height * mn_char(_from_beta(newb), rest)
    return total


def cycle_type(perm):
    n, seen, ct = len(perm), [False] * len(perm), []
    for i in range(n):
        if not seen[i]:
            L, j = 0, i
            while not seen[j]:
                seen[j] = True; j = perm[j]; L += 1
            ct.append(L)
    return tuple(sorted(ct, reverse=True))


def class_data(n):
    """(cycle type → class size) for S_n."""
    from math import factorial
    out = {}
    for rho in partitions(n):
        z = 1
        for m in set(rho):
            c = rho.count(m)
            z *= (m ** c) * factorial(c)
        out[rho] = factorial(n) // z
    return out


def inner(n, chi1, chi2):
    """⟨χ1,χ2⟩ over classes (both — dict cycle-type→value). Exact ℚ."""
    cd = class_data(n)
    from math import factorial
    tot = sum(sp.Integer(cd[r]) * chi1[r] * chi2[r] for r in cd)
    return sp.Rational(tot, factorial(n))


def decompose(n, chi):
    """The decomposition of a character into irreducibles λ ⊢ n. Returns {λ: multiplicity}."""
    out = {}
    for lam in partitions(n):
        m = inner(n, chi, {r: sp.Integer(mn_char(lam, r)) for r in class_data(n)})
        if m != 0:
            out[lam] = m
    return out


def std_chars(n):
    """χ_std, χ_Sym²(std), χ_Λ²(std) over classes — exact integers."""
    cd = class_data(n)
    chi_std, chi_s2, chi_l2 = {}, {}, {}
    for rho in cd:
        fix = rho.count(1)
        c = fix - 1                                     # χ_std = fix − 1
        # the cycle type of the square of a permutation: a cycle of length m → gcd(m,2) cycles of length m/gcd
        sq = []
        for m in rho:
            if m % 2 == 0:
                sq += [m // 2, m // 2]
            else:
                sq.append(m)
        sq = tuple(sorted(sq, reverse=True))
        c2 = sq.count(1) - 1
        chi_std[rho] = sp.Integer(c)
        chi_s2[rho] = sp.Rational(c * c + c2, 2)
        chi_l2[rho] = sp.Rational(c * c - c2, 2)
    return chi_std, chi_s2, chi_l2


def name_lam(lam):
    return "(" + ",".join(map(str, lam)) + ")"


# ══════════════ PART B (revision-2): COMPONENTS via SAMPLE+NEWTON, not via a grid ══════════════
# ★WHY REVISION-2.  Revision-1 looked for components with a mask |f|<eps on a grid and FAILED both positive
#   controls: d=2 gave 6 components instead of 2 (the shell is thin ⟹ the mask crumbles into fragments),
#   d=4 gave 1 (the shell is thick ⟹ everything merges).  That is, the number of components there was governed
#   by a THRESHOLD, not by the object — a hidden handle of precisely the class we have been catching all week.
#   The controls fired as intended: the MACHINERY failed, and it was replaced, not tuned.
#   Revision-2 places points EXACTLY on the manifold (Gauss-Newton down to |f|<1e-12) and builds a neighbor graph;
#   the radius-handle is NOT hidden — the number of components is printed as a FUNCTION of the radius (a plateau or its
#   absence is visible by eye, and the absence of a plateau = an honest negative, not a choice of number).

def _F_and_J(u_free, d):
    """f = Σ_j e^{i u_j} with u_0 = −Σ_{j≥1} u_j.  Returns (the residual ℝ², the Jacobian 2×d)."""
    u0 = -u_free.sum()
    u = np.concatenate(([u0], u_free))
    e = np.exp(1j * u)
    F = np.array([e.sum().real, e.sum().imag])
    # ∂u_0/∂u_k = −1 ⟹ ∂f/∂u_k = i e^{iu_k} − i e^{iu_0}
    g = 1j * e[1:] - 1j * e[0]
    J = np.vstack([g.real, g.imag])
    return F, J


def sample_node_set(d, npts, rng, tol=1e-12, iters=200):
    """Points EXACTLY on the node-set: Gauss-Newton from random starts."""
    out = []
    while len(out) < npts:
        u = rng.uniform(0, 2 * np.pi, size=d)
        for _ in range(iters):
            F, J = _F_and_J(u, d)
            if np.linalg.norm(F) < tol:
                break
            du, *_ = np.linalg.lstsq(J, -F, rcond=None)
            step = 1.0
            for _ls in range(30):                      # a damper: step only if the residual decreases
                Fn, _ = _F_and_J(u + step * du, d)
                if np.linalg.norm(Fn) < np.linalg.norm(F):
                    break
                step *= 0.5
            u = u + step * du
        F, _ = _F_and_J(u, d)
        if np.linalg.norm(F) < 1e-9:
            out.append(u % (2 * np.pi))
    return np.array(out)


def _torus_d2(A, B):
    dif = (A[:, None, :] - B[None, :, :] + np.pi) % (2 * np.pi) - np.pi
    return (dif ** 2).sum(axis=-1)


def components_by_radius(pts, r):
    """Components of the neighbor graph of radius r (a torus-metric). Union-find."""
    n = len(pts)
    D2 = _torus_d2(pts, pts)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)

    ii, jj = np.nonzero(D2 < r * r)
    for a, b in zip(ii, jj):
        if a < b:
            union(int(a), int(b))
    lab, remap = np.zeros(n, dtype=int), {}
    for i in range(n):
        rt = find(i)
        if rt not in remap:
            remap[rt] = len(remap) + 1
        lab[i] = remap[rt]
    return lab, len(remap)


def act_sigma(pts, sigma, d):
    """σ ∈ S_{d+1} permutes the FULL tuple (u_0..u_d); Σu=0 is preserved."""
    u0 = -pts.sum(axis=1)
    full = np.column_stack([u0, pts])
    perm = full[:, list(sigma)]
    return perm[:, 1:] % (2 * np.pi)


def component_action_v2(d, pts, r):
    lab, ncomp = components_by_radius(pts, r)
    reps = {c: pts[np.nonzero(lab == c)[0][0]] for c in range(1, ncomp + 1)}
    perms = {}
    for sigma in permutations(range(d + 1)):
        img = {}
        for c in range(1, ncomp + 1):
            w = act_sigma(reps[c][None, :], sigma, d)
            j = int(np.argmin(_torus_d2(w, pts)[0]))
            img[c] = int(lab[j])
        perms[sigma] = img
    return ncomp, perms, lab


def perm_character(d, ncomp, perms):
    """χ(σ) = the number of fixed components, grouped by the cycle type of σ."""
    chi = {}
    for sigma, img in perms.items():
        ct = cycle_type(list(sigma))
        val = sum(1 for c in range(1, ncomp + 1) if img[c] == c)
        if ct in chi and chi[ct] != val:
            chi[ct] = None
        else:
            chi[ct] = val
    return chi


# ══ ★REVISION-3 — THE CORRECT TORUS.  Revision-2 sampled u ∈ {Σu≡0 mod 2π}, and the control d=2 gave 6 points
#    instead of 2 — STABLY at all radii, that is, this is not noise.  Diagnosis: this u-torus is a 3-fold
#    COVER of the zone.  The periodicity of the phases is set not by 2πℤ^{d+1}, but by the lattice 2π·{m − m̄·1} (a projection
#    of ℤ^{d+1}), in which shifts have FRACTIONAL components.  Two different rulers under one name «torus» —
#    exactly the same class as the homonym of length.  Revision-3 works in ZONE coordinates: k = 2π Σ x_i δ_i,
#    x ∈ [0,1)^d, φ_j = k·δ_j = 2π Σ_i x_i (δ_ij − 1/(d+1)).  Now d=2 must give exactly 2.

def phase_matrix(d):
    """G[i][j] = δ_i·δ_j = δ_ij − 1/(d+1), i=1..d (a basis of the zone), j=0..d (bonds)."""
    G = np.zeros((d, d + 1))
    for i in range(1, d + 1):
        for j in range(d + 1):
            G[i - 1, j] = (1.0 if i == j else 0.0) - 1.0 / (d + 1)
    return G


def _F_and_J_bz(x, d, G):
    phi = 2 * np.pi * (G.T @ x)                      # (d+1,)
    e = np.exp(1j * phi)
    F = np.array([e.sum().real, e.sum().imag])
    dphi = 2 * np.pi * G.T                           # ∂φ_j/∂x_i
    g = (1j * e)[:, None] * dphi                     # (d+1, d)
    Jc = g.sum(axis=0)
    return F, np.vstack([Jc.real, Jc.imag])


def sample_bz(d, npts, rng, iters=200):
    G = phase_matrix(d)
    out = []
    guard = 0
    while len(out) < npts and guard < 200 * npts:
        guard += 1
        x = rng.uniform(0, 1, size=d)
        for _ in range(iters):
            F, J = _F_and_J_bz(x, d, G)
            if np.linalg.norm(F) < 1e-13:
                break
            dx, *_ = np.linalg.lstsq(J, -F, rcond=None)
            step = 1.0
            for _ls in range(40):
                Fn, _ = _F_and_J_bz(x + step * dx, d, G)
                if np.linalg.norm(Fn) < np.linalg.norm(F):
                    break
                step *= 0.5
            x = x + step * dx
        F, _ = _F_and_J_bz(x, d, G)
        if np.linalg.norm(F) < 1e-9:
            out.append(x % 1.0)
    return np.array(out)


def _unit_d2(A, B):
    dif = (A[:, None, :] - B[None, :, :] + 0.5) % 1.0 - 0.5
    return (dif ** 2).sum(axis=-1)


def components_bz(pts, r):
    n = len(pts)
    D2 = _unit_d2(pts, pts)
    parent = list(range(n))

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]; v = parent[v]
        return v

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    ii, jj = np.nonzero(D2 < r * r)
    for a, b in zip(ii, jj):
        if a < b:
            union(int(a), int(b))
    lab, remap = np.zeros(n, dtype=int), {}
    for i in range(n):
        rt = find(i)
        if rt not in remap:
            remap[rt] = len(remap) + 1
        lab[i] = remap[rt]
    return lab, len(remap)


def act_sigma_bz(x, sigma, d, G):
    """σ permutes the BONDS: φ'_j = φ_{σ(j)}.  We return to x by least squares
    (the system is consistent, because Σφ=0 is preserved)."""
    phi = 2 * np.pi * (G.T @ x)
    phi2 = phi[list(sigma)]
    xn, *_ = np.linalg.lstsq(2 * np.pi * G.T, phi2, rcond=None)
    return xn % 1.0


def component_action_bz(d, pts, r):
    G = phase_matrix(d)
    lab, ncomp = components_bz(pts, r)
    reps = {c: pts[np.nonzero(lab == c)[0][0]] for c in range(1, ncomp + 1)}
    perms = {}
    for sigma in permutations(range(d + 1)):
        img = {}
        for c in range(1, ncomp + 1):
            w = act_sigma_bz(reps[c], sigma, d, G)
            j = int(np.argmin(_unit_d2(w[None, :], pts)[0]))
            img[c] = int(lab[j])
        perms[sigma] = img
    return ncomp, perms


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1059_2_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf); sys.stdout = tee
    sys.path.insert(0, os.path.join(_HERE, ".."))

    print("=" * 82)
    print("S1059-B — THE ∀d-LAW OF THE MEDIAL: does the choice of the Sym²/Λ² half follow from an ANCESTOR")
    print("Layer-1, a bond-index.  Order: KILL-FIRST, and only then the law.")
    print("=" * 82); print()

    # ── A. exact decompositions of the square of the bond-space ──
    print("A. THE SQUARE OF THE BOND-SPACE, EXACTLY (Murnaghan–Nakayama, integer characters)")
    sqr = {}
    for n in range(3, 8):
        cs, s2, l2 = std_chars(n)
        d_s2, d_l2 = decompose(n, s2), decompose(n, l2)
        sqr[n] = (d_s2, d_l2)
        print(f"   n={n} (d={n-1}):  Sym²(std) = " + " ⊕ ".join(
            (f"{v}·" if v != 1 else "") + name_lam(k) for k, v in sorted(d_s2.items(), reverse=True)))
        print(f"            Λ²(std)  = " + " ⊕ ".join(
            (f"{v}·" if v != 1 else "") + name_lam(k) for k, v in sorted(d_l2.items(), reverse=True)))
    # a tooth: the known d=2,3 from S1052
    z1 = sqr[3][1] == {(1, 1, 1): 1}
    z2 = (2, 2) in sqr[4][0] and (2, 2) not in sqr[4][1]
    ok(z1, "★tooth-A1: Λ²(std) at n=3 = (1,1,1) = sign — exactly where S1052 found the node-mark d=2 "
           "[A · dimensionless · a character decomposition]")
    ok(z2, "★tooth-A2: (2,2)=E at n=4 sits in Sym² and NOT in Λ² — exactly the S1052 verdict for d=3 "
           "[A · dimensionless · a character decomposition]")
    print()

    # ── B (revision-2). components of the node-set + the action ──
    print("B. π₀(THE NODE-SET) + THE ACTION OF S_{d+1} — a sample ON the manifold, the radius NOT hidden")
    print("   ★revision-1 (a mask |f|<eps on a grid) FAILED both controls and was REPLACED, not tuned.")
    rng = np.random.default_rng(20260722)
    results, RADII = {}, (0.15, 0.25, 0.35, 0.5, 0.7)
    for d in (2, 3, 4, 5):
        npts = {2: 400, 3: 1200, 4: 2500, 5: 3500}[d]
        pts = sample_bz(d, npts, rng)
        counts = [components_bz(pts, r)[1] for r in RADII]
        print(f"   d={d}: components(r) = " + ", ".join(f"r={r}:{c}" for r, c in zip(RADII, counts))
              + f"  [B · a component · union-find on the sample; COVERAGE: {npts} points with |f|<1e-9, "
                f"the zone metric (period 1), {npts} random starts]")
        # a plateau = the longest stable stretch of values
        best_val, best_len, cur_val, cur_len = None, 0, None, 0
        for c in counts:
            if c == cur_val:
                cur_len += 1
            else:
                cur_val, cur_len = c, 1
            if cur_len > best_len:
                best_val, best_len = cur_val, cur_len
        stable = best_len >= 2
        ncomp = best_val
        print(f"        plateau: {ncomp} component(s) over {best_len} of {len(RADII)} radii "
              f"⟹ {'stable' if stable else '★UNSTABLE — not reported as a number'}")
        if not stable or ncomp is None or ncomp > 8:
            results[d] = (ncomp, None)
            continue
        r_use = RADII[counts.index(ncomp)]
        nc, perms = component_action_bz(d, pts, r_use)
        chi = perm_character(d, nc, perms)
        n = d + 1
        cls = class_data(n)
        if all(chi.get(rr) is not None for rr in cls):
            dec = decompose(n, {rr: sp.Integer(chi[rr]) for rr in cls})
            results[d] = (nc, dec)
            print("        ℂ[π₀] = " + " ⊕ ".join(
                (f"{v}·" if v != 1 else "") + name_lam(k) for k, v in sorted(dec.items(), reverse=True))
                + "  [B · dimensionless · a decomposition of the permutation character]")
        else:
            results[d] = (nc, None)
            print("        ⚠ the character is inconsistent within a conjugacy class — not decomposed, not reported")
    print()

    d2 = results.get(2, (0, None))
    ok(d2[0] == 2 and d2[1] == {(3,): 1, (1, 1, 1): 1},
       "★tooth-B1 (a positive control): d=2 → 2 components, ℂ[π₀] = 1 ⊕ sign — the machinery "
       "reproduces the KNOWN result (S1052: the node-mark d=2 = orientation)")
    d3 = results.get(3, (0, None))
    print("   ★tooth-B2 is SPLIT in revision-3.  The sample gives d=3 → 1 component, whereas Test-3 was read")
    print("     as «3 components».  This is not a vote: below, an EXACT test decides whether the «3» were")
    print("     components of CONNECTIVITY, or BRANCHES that intersect.  I do NOT report the number from the sample")
    print("     until the exact test says what it means.")
    machinery_ok = (d2[0] == 2)
    print(f"   the machinery on the d=2 control: {'PASSED' if machinery_ok else 'DID NOT PASS'} "
          f"(exactly 2 components = the known node-structure of d=2)")
    print()

    # ── C. KILL-FIRST ──
    print("C. KILL-FIRST — does the node-representation embed in std⊗std = Sym² ⊕ Λ²")
    if not machinery_ok:
        print("   ★THE MACHINERY DID NOT PASS THE CONTROLS ⟹ results for d=4,5 are NOT reported (K-machinery).")
    verdict_rows = []
    for d in (2, 3, 4, 5):
        n = d + 1
        ncomp, dec = results.get(d, (0, None))
        if not isinstance(dec, dict) or not dec:
            print(f"   d={d}: no decomposition (unstable or not computed) — skipping, "
                  f"NOT reported as a result")
            continue
        s2d, l2d = sqr[n]
        nontriv = {k: v for k, v in dec.items() if k != (n,)}
        if not nontriv:
            verdict_rows.append((d, ncomp, {}, "no mark (π₀ trivial)", True))
            print(f"   d={d}: ℂ[π₀] is trivial ⟹ ★THE NODE-MARK DOES NOT EXIST — the question «which half» "
                  f"loses its subject")
            continue
        in_s2 = all(int(s2d.get(k, 0)) >= int(v) for k, v in nontriv.items())
        in_l2 = all(int(l2d.get(k, 0)) >= int(v) for k, v in nontriv.items())
        tot = {k: int(s2d.get(k, 0)) + int(l2d.get(k, 0)) for k in set(s2d) | set(l2d)}
        in_sq = all(tot.get(k, 0) >= int(v) for k, v in nontriv.items())
        where = "Sym²" if in_s2 else ("Λ²" if in_l2 else ("scattered" if in_sq else "★NOT IN THE SQUARE"))
        verdict_rows.append((d, ncomp, nontriv, where, in_sq))
        print(f"   d={d}: the nontrivial part = " + " ⊕ ".join(
            name_lam(k) for k in sorted(nontriv, reverse=True)) + f"  →  {where}")
    print()

    # ── D. EXACTLY: whether «3» is components or branches; and whether the mark has a subject at all at d≥4 ──
    print("D. THE EXACT TEST (symbolic, without a grid and without a sample) — what «3» means at d=3")
    # ★REVISION-2 PER BETA'S VISA (S1059B_MEDIAL_LAW_VISA_BETA.md).  Revision-1 presented the witness z=(1,-1,1,-1),
    #   checking Sz=0 and Pz=1 — conditions of the OLD u-parametrization, the very same one that this
    #   probe itself discarded as a 3-fold cover.  In the zone this same vector has no preimage in x[0,1)^3;
    #   it has a CLASS that it represents (z is defined only UP TO A COMMON PHASE, because the shift x_i->x_i+1
    #   gives the same dphi for all j ⟹ f->lambda*f).  So the witness is presented WITH ITS PREIMAGE in the zone.
    x_wit = [sp.Rational(1, 2), sp.Integer(0), sp.Rational(1, 2)]
    Gm = sp.Matrix(3, 4, lambda i, j: (1 if i + 1 == j else 0) - sp.Rational(1, 4))
    phi = [2 * sp.pi * sum(x_wit[i] * Gm[i, j] for i in range(3)) for j in range(4)]
    zz = [sp.simplify(sp.exp(sp.I * ph)) for ph in phi]
    print(f"     the witness IN THE ZONE: x = (1/2, 0, 1/2) ⟹ φ = {[sp.nsimplify(pp) for pp in phi]}")
    print(f"                    ⟹ z = {zz}  (the class of the vector (1,−1,1,−1) up to a common phase)")
    ok(sp.simplify(sum(zz)) == 0, "★the witness has a PREIMAGE in the zone x∈[0,1)³ and Σz=0 exactly "
                                  "[D · dimensionless · a symbolic sum]")
    ok(sp.simplify(sp.prod(zz)) == 1,
       "★Πz=1 holds — ★but it is origin-choice-DEPENDENT (here λ^(d+1)=(−i)⁴=1); "
       "not to be cited bare [D · dimensionless · a symbolic product]")
    pair_a = sp.simplify(zz[0] + zz[1]) == 0 and sp.simplify(zz[2] + zz[3]) == 0
    pair_b = sp.simplify(zz[0] + zz[3]) == 0 and sp.simplify(zz[1] + zz[2]) == 0
    pair_c = sp.simplify(zz[0] + zz[2]) == 0 and sp.simplify(zz[1] + zz[3]) == 0
    ok(pair_a and pair_b,
       "★★THE WITNESS BELONGS TO TWO DIFFERENT PAIRINGS simultaneously ({01}{23} AND {03}{12}) ⟹ the branches "
       "INTERSECT [D · dimensionless · a check of both pairings in the zone]")
    ok(not pair_c, "★and does NOT belong to the third ({02}{13}) — the mark does not degenerate completely")

    # ★AN EXACT BOUNDARY (a visa addendum): belonging to ALL THREE pairings is impossible ∀z, not only here.
    za, zb, zc, zd = sp.symbols('z_a z_b z_c z_d')
    forced = sp.solve([za + zb, za + zc, za + zd], [zb, zc, zd], dict=True)[0]
    sum_if_all_three = sp.simplify(za + forced[zb] + forced[zc] + forced[zd])
    ok(sp.simplify(sum_if_all_three + 2 * za) == 0,
       "★★AN EXACT BOUNDARY: if a point belonged to ALL THREE pairings, then z₁=z₂=z₃=−z₀ ⟹ Σz=−2z₀≠0 "
       "⟹ at most TWO simultaneously [D · dimensionless · a symbolic solution of the system]")
    print("     ⟹ the branches intersect PAIRWISE and no more — the boundary is exact, not observed.")
    print("   DOES THE MARK HAVE A SUBJECT AT d≥4 — pure arithmetic, not numerics:")
    ok((3 + 1) % 2 == 0 and (4 + 1) % 2 == 1,
       "★d=3 ⟹ d+1=4 EVEN (a decomposition into antipodal pairs is possible); d=4 ⟹ d+1=5 ODD "
       "(a decomposition into pairs is IMPOSSIBLE at all) [D · a bond · the parity of d+1]")
    # at n=5 the solutions are NOT forced to one type: two structurally different, both exact
    w5 = [sp.exp(2 * sp.pi * sp.I * sp.Rational(j, 5)) for j in range(5)]
    # ★sp.simplify does NOT collapse a sum of roots of unity without expand_complex — the COLLAPSING TOOL
    #   was failing, not the claim (Σ of the roots z⁵=1 = 0 is a theorem).  The call was fixed, not the test.
    t1 = sp.simplify(sp.expand_complex(sum(w5))) == 0
    om = sp.exp(2 * sp.pi * sp.I / 3)
    mixed = [sp.Integer(1), om, om ** 2, sp.I, -sp.I]          # a triangle + an antipodal pair
    t2 = sp.simplify(sp.expand_complex(sum(mixed))) == 0
    ok(t1 and t2,
       "★at n=5 there exist TWO structurally different exact solutions: {5th-degree roots} and "
       "{a triangle ⊕ an antipodal pair} — the type is NOT forced [D · dimensionless · symbolic sums]")

    # ★A HOLE, CLOSED BY THE VISA: d=5 ⟹ n=6 is EVEN, so the odd-parity argument does NOT apply to it,
    #   and the second argument was only run for n=5.  Revision-1 covered d=5 with a SAMPLE and called it
    #   weaker than it was: not «π₀ is numeric», but «there is no exact argument for d=5».
    w6 = [sp.exp(2 * sp.pi * sp.I * sp.Rational(j, 6)) for j in range(6)]
    s6a = sp.simplify(sp.expand_complex(sum(w6))) == 0
    om3 = sp.exp(2 * sp.pi * sp.I / 3)
    aa = sp.Symbol('a_ph', real=True)
    tri2 = [1, om3, om3 ** 2, sp.exp(sp.I * aa), sp.exp(sp.I * aa) * om3,
            sp.exp(sp.I * aa) * om3 ** 2]
    s6b = sp.simplify(sp.expand_complex(sum(tri2))) == 0
    p1, p2, p3 = sp.symbols('p1 p2 p3', real=True)
    pairs3 = [sp.exp(sp.I * p1), -sp.exp(sp.I * p1), sp.exp(sp.I * p2), -sp.exp(sp.I * p2),
              sp.exp(sp.I * p3), -sp.exp(sp.I * p3)]
    s6c = sp.simplify(sum(pairs3)) == 0
    ok(s6a and s6b and s6c,
       "★★d=5 IS CLOSED BY ARITHMETIC, not by a sample: at n=6 (EVEN) there exist THREE structurally different "
       "exact families — {6th-degree roots} ⊥ {two triangles, a free phase} ⊥ {three "
       "antipodal pairs, THREE free phases} [D · dimensionless · symbolic sums with free phases]")
    print("     ⟹ the third family carries FREE PHASES, that is, a CONTINUUM ⟹ there is no forcing at all,")
    print("       even more so than at n=5.  The parity of d+1 by itself saves nothing: a decomposition")
    print("       into pairs BECOMES possible, but stops being UNIQUE — and the mark does not live off")
    print("       the possibility of the decomposition, but off its BEING FORCED.")
    print("     ⟹ at d≥4 there is no combinatorial forcing: the solution has no unique TYPE,")
    print("       so there is simply NOTHING TO MARK with a mark that could be placed in Sym² or Λ².")
    print("     ⟹ d=2 (orientation, forced) and d=3 (a pairing, forced) — are not two values")
    print("       of one function of d, but TWO DIFFERENT MECHANISMS, and a third does not exist.")
    print()

    print("=" * 82)
    print("★THE VERDICT OF THE PROBE")
    if not machinery_ok:
        print("  K-MACHINERY: the d=2 control did not pass ⟹ there is NO verdict about the law.")
    else:
        print("  ★ATTENTION TO ONE'S OWN DATA: the line «there is no mark at d=3» would be WRONG.  π₀ at d=3")
        print("  is trivial, BUT part D showed that the mark IS there — it is a BRANCH-mark, not a")
        print("  topological one.  So the verdict reads BOTH sources, not the more convenient one:")
        print()
        print("  | d | what is the mark | nature | half |")
        print("  | 2 | orientation (2 π₀-components) | TOPOLOGICAL | Λ² |")
        print("  | 3 | a pairing (3 branches, π₀=1)  | COMBINATORIAL | Sym² |")
        print("  | ≥4| neither one nor the other    | —            | no subject |")
        print()
        print("  ⟹ ★THERE IS NO ∀d-LAW IN THE SOUGHT FORM, and the reason is exact: the mark exists precisely where")
        print("  the configuration is combinatorially FORCED, that is, at d+1 ≤ 4.  At d+1=3 the forcing")
        print("  gives a discrete set (hence the topological mark), at d+1=4 — a decomposition into")
        print("  antipodal pairs (hence the combinatorial one).  At d+1≥5 there is no forcing at all:")
        print("  two structurally different exact solutions at n=5 have been exhibited.")
        print("  ⟹ Λ² at d=2 and Sym² at d=3 — are NOT two values of a function of d, but TWO DIFFERENT")
        print("  MECHANISMS on the only two d where the subject exists.  The question «which half ∀d»")
        print("  has no answer not for lack of computation, but through the EXHAUSTION of the subject.")
        print("  ★This is a named-boundary of the problem-front, not its solution — and it is reported precisely as such.")
    print(f"SCORE: {len(_passes)} ✓ / {len(_fails)} ✗")
    print("=" * 82)
    sys.stdout = tee.real
    logf.close()
    return 1 if _fails else 0


if __name__ == "__main__":
    sys.exit(main())
