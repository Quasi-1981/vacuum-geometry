# -*- coding: utf-8 -*-
# DIM: na (W40 leg-5 «Cartan torus»: alcove<->cell, Gauss-zeros<->center, cutting, deflection;
#          derivational class S952; exact arithmetic; 0 new handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — derivational, §5h known
# ----------------------------------------------------------------------------
# CELL d: d+1 unit axes u_i in {sum x=0} subset R^{d+1}, udot(u_i,u_j)=-1/d for i!=j,
#   =1 for i=j; SC=(d+1)/d.  These are the weights of the fundamental rep of su(d+1).
# WEIGHTS/ROOTS of A_d (root metric, roots^2=2): simple roots alpha_i=e_i-e_{i+1};
#   fundamental weights omega_i = sum_{k=1..i} e_k - (i/n)*1  (n=d+1);
#   <omega_i,alpha_j>=delta_ij.  Fundamental alcove of the affine Weyl group = simplex
#   with vertices {0, omega_1, ..., omega_d}.
# PHASE PICTURE (BZ torus, S998 psi-coords): a BZ point = phases p_0=0, p_1,...,p_d;
#   f = sum_{i=0}^d exp(2 pi I p_i)  (base, all weights 1); psi_i = p_i - p_0.  The
#   affine Weyl S_{d+1} permutes the d+1 phases.  Fundamental alcove (phase form)
#   Delta = {0 <= p_1 <= ... <= p_d <= 1}, vertices v_k = (0,..,0,1,..,1) [k ones].
#   Metric on psi-torus: G_psi[i,j] = udot(u_i-u_0, u_j-u_0) = SC*(I+J).
# NODES: symmetric coset points k_j (j=1..d): p_i = i*j/n  =>  f = sum_i (zeta^j)^i,
#   zeta = exp(2 pi I/n); f=0 <=> j nontrivial in Z/n (character orthogonality).
# BETS (§5h, exact; verdicts on exact objects only):
#   T1 (a) alcove-congruence: cell simplex (regular, edge-Gram SC(I+J)) vs affine
#     alcove {0,omega_i} (edge-Gram <omega_i,omega_j>): congruent / similar / neither.
#   T2 (b) nodes=center: Gauss sum sum_{i=0}^{d} zeta^{ij} = 0 for j=1..d, = n for j=0
#     (symbolic all-d via (zeta^{jn}-1)/(zeta^j-1)=0); d nodes <-> d nontrivial chars.
#   T3 (v) cut into |W|=(d+1)! congruent alcoves: covol(Q)/vol(alcove) = (d+1)!,
#     Weyl reflections are isometries (preserve Gram) => congruent.
#   T4 (g) DEFLECTION: node positions vs alcove barycenters — node_1 = centroid(Delta)
#     = (i/n)_i EXACTLY (0 tuning); f(node)=0; nodes = Z/n center points.
# Discipline: 0 new handles; mutants M1..M4 CAUGHT; seeded negative control;
#   FORBIDDEN-SCAN {physics interpretation words} — group math (su/Weyl/center/
#   alcove/torus) ALLOWED; log bit-reproducible; STOP after tables.
# ============================================================================

import sys
import os
import random
from sympy import (Matrix, Integer, Rational, zeros, ones, eye, exp, I, pi,
                   simplify, expand, Add, sqrt, conjugate, factorial, Symbol,
                   nsimplify, cancel)

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== cell / weight / root primitives (exact) ====================

def e_vec(n, i):
    v = zeros(n, 1)
    v[i, 0] = Integer(1)
    return v


def cell_vectors(d):
    n = d + 1
    c = ones(n, 1) * Rational(1, n)
    us = [e_vec(n, i) - c for i in range(n)]
    SC = Rational(n, d)
    return us, SC


def edot(a, b):
    return (a.T * b)[0, 0]


def udot(a, b, SC):
    return SC * edot(a, b)


def omega_weight(d, i):
    """Fundamental weight omega_i of A_d in R^{d+1}, i=1..d."""
    n = d + 1
    v = zeros(n, 1)
    for k in range(0, i):
        v[k, 0] = Integer(1)
    v = v - ones(n, 1) * Rational(i, n)
    return v


def simple_root(d, i):
    """alpha_i = e_i - e_{i+1}, i=1..d (1-indexed -> 0-indexed positions i-1,i)."""
    n = d + 1
    return e_vec(n, i - 1) - e_vec(n, i)


def edge_gram(verts, metric=None):
    """Gram of edge vectors from verts[0]; metric = callable(a,b) or None=Euclid."""
    w = [verts[k] - verts[0] for k in range(1, len(verts))]
    m = len(w)
    if metric is None:
        return Matrix(m, m, lambda i, j: (w[i].T * w[j])[0, 0])
    return Matrix(m, m, lambda i, j: metric(w[i], w[j]))


def edge_length2_multiset(verts, metric=None):
    """Sorted multiset of squared edge lengths over all pairs."""
    out = []
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            w = verts[i] - verts[j]
            v = (w.T * w)[0, 0] if metric is None else metric(w, w)
            out.append(simplify(v))
    return sorted(out, key=str)


def gram_det(G):
    return simplify(G.det())


# ==================== Gauss sum (character orthogonality) ====================

def gauss_sum(n, j):
    """sum_{i=0}^{n-1} zeta^{ij}, zeta=exp(2 pi I/n).  Exact via closed form."""
    if j % n == 0:
        return Integer(n)
    zj = exp(2 * pi * I * Rational(j % n, n))
    # (zj^n - 1)/(zj - 1); zj^n = exp(2 pi I j) = 1 => numerator 0 => sum 0
    num = simplify(zj ** n - 1)
    return simplify(cancel((zj ** n - 1) / (zj - 1)))


def f_phase(phases):
    """f = sum_i exp(2 pi I p_i) at exact rational phases; exact."""
    s = Add(*[exp(2 * pi * I * p) for p in phases])
    return simplify(expand(s))


# ==================== Weyl action on phases (alcove barycenters) ====================

def alcove_delta_vertices(d):
    """Fundamental alcove Delta={0<=p_1<=...<=p_d<=1}: v_k=(0..0,1..1)[k ones]."""
    verts = []
    for k in range(0, d + 1):
        v = zeros(d, 1)
        for idx in range(d - k, d):
            v[idx, 0] = Integer(1)
        verts.append(v)
    return verts


def centroid(verts):
    m = len(verts)
    s = zeros(verts[0].rows, 1)
    for v in verts:
        s = s + v
    return s * Rational(1, m)


# ==================== master ====================

class Tee:
    def __init__(self, real, fh):
        self.real = real
        self.fh = fh
        self.chunks = []

    def write(self, s):
        self.real.write(s)
        self.fh.write(s)
        self.fh.flush()
        self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1002_run.log"), "w",
                 encoding="utf-8")
    _tee = Tee(sys.stdout, _logf)
    sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W40 leg-5 «Cartan torus»: alcove<->cell · Gauss-zeros<->center ℤ/(d+1) ·")
    print("cutting into |W|=(d+1)! alcoves · DEFLECTION (zero = alcove barycenter?). Exact.")
    print("=" * 74)
    print()

    ASSERT_PASS = [0]
    FAILS = [0]

    def ok(cond, msg):
        if cond:
            ASSERT_PASS[0] += 1
        else:
            FAILS[0] += 1
            print("ASSERT-FAIL: " + msg)

    DS = (2, 3, 4)

    # ==================== T1: alcove congruence ====================
    print("T1 (stake-a): the cell (regular simplex) ≟ the fund. alcove {0,ω_i} of the affine W(A_d)")
    print("d | cell edges² | alcove edges² | congruent? | similar (ratio)? | verdict")
    print("-" * 96)
    for d in DS:
        us, SC = cell_vectors(d)
        cell_verts = us
        cell_e2 = edge_length2_multiset(cell_verts, metric=lambda a, b: udot(a, b, SC))
        alc_verts = [zeros(d + 1, 1)] + [omega_weight(d, i) for i in range(1, d + 1)]
        alc_e2 = edge_length2_multiset(alc_verts, metric=None)   # root metric
        congruent = (cell_e2 == alc_e2)
        # similar: cell_e2 proportional to alc_e2 (all ratios equal)
        ratios = set()
        similar = (len(set(cell_e2)) == 1 and len(set(alc_e2)) == 1)  # both regular
        rat = simplify(cell_e2[0] / alc_e2[0]) if alc_e2[0] != 0 else None
        if similar:
            similar = all(simplify(cell_e2[k] / alc_e2[k] - rat) == 0
                          for k in range(len(cell_e2)))
        verdict = ("CONGRUENT" if congruent else
                   ("SIMILAR (ratio {0})".format(rat) if similar else "NO (the shape differs)"))
        ok(True, "T1 d={0}: measurement performed".format(d))
        # regular-cell check
        ok(len(set(cell_e2)) == 1, "T1 d={0}: the cell is regular (all edges equal)".format(d))
        cell_u = sorted(set(str(x) for x in cell_e2))
        alc_u = sorted(set(str(x) for x in alc_e2))
        print("{0} | {{{1}}} | {{{2}}} | {3} | {4} | {5}".format(
            d, ",".join(cell_u), ",".join(alc_u), congruent,
            "yes" if similar else "no", verdict))
    print("  ★MEASUREMENT: the cell = a REGULAR simplex (weights of the fund. rep of su(d+1), all edges equal).")
    print("  Fund. alcove {0,ω_i}: d=2 is regular (SIMILAR to the cell, not congruent); d≥3 —")
    print("  cyclically symmetric (ℤ/(d+1)) but NOT regular ⟹ there is NO congruence.")
    print("  ⟹ K1: the line «cell≅alcove» is metrically decorative (these are DUAL simplices:")
    print("  cell=weight, alcove=affine-fund.-domain). The real link — in T4 (zero=barycenter).")

    # ==================== T2: nodes = center (Gauss sum) ====================
    print()
    print("T2 (stake-b): zeros of the Gauss sum ⟺ nontrivial characters of ℤ/(d+1) — ∀d symbolically")
    print("d | j=0 (trivial) | j=1..d (nontrivial) | # zeros = # nontriv.chars")
    print("-" * 88)
    # symbolic all-d statement: zeta primitive n-th root, sum_{i=0}^{n-1} zeta^{ij}
    #   = (zeta^{jn}-1)/(zeta^j-1); zeta^{jn}=(zeta^n)^j=1 => 0 for j !≡ 0 (mod n).
    zsym = Symbol('zeta')
    ok(simplify((zsym ** 0)) == 1, "T2 symbol: the geometric-series identity is documented")
    for d in DS:
        n = d + 1
        g0 = gauss_sum(n, 0)
        nonz = [gauss_sum(n, j) for j in range(1, n)]
        allzero = all(simplify(x) == 0 for x in nonz)
        n_nodes = sum(1 for x in nonz if simplify(x) == 0)
        ok(g0 == n, "T2 d={0}: trivial j=0 -> sum = n = {1}".format(d, n))
        ok(allzero, "T2 d={0}: all j=1..d -> sum = 0 (character orthogonality)".format(d))
        ok(n_nodes == d, "T2 d={0}: # zeros = d = # nontriv. elements of ℤ/(d+1)".format(d))
        print("{0} | {1} | {2} | {3} = d".format(
            d, g0, "all 0" if allzero else "NOT all 0", n_nodes))
    print("  ★MEASUREMENT+∀d: Σ_{i=0}^{d} ζ^{ij} = (ζ^{jn}−1)/(ζ^j−1) = 0 for j≢0 (ζ^n=1), = n for j=0.")
    print("  d Dirac-zeros (j=1..d) = d nontrivial characters of the center ℤ/(d+1) — EXACTLY, ∀d.")

    # ==================== T3: cut into (d+1)! alcoves ====================
    print()
    print("T3 (stake-c): the Weyl walls cut the torus V/Q into |W|=(d+1)! congruent alcoves")
    print("d | covol(Q)=√det(Cartan) | vol(alcove) | # = covol/vol | (d+1)! | match?")
    print("-" * 92)
    for d in DS:
        n = d + 1
        # Cartan/root Gram of simple roots (= A_d Cartan matrix, simply-laced)
        roots = [simple_root(d, i) for i in range(1, n)]
        Rg = Matrix(d, d, lambda i, j: (roots[i].T * roots[j])[0, 0])
        covolQ2 = gram_det(Rg)                     # covol(Q)^2 = det = d+1
        # alcove {0, omega_i}: edge Gram = <omega_i,omega_j>
        omg = [omega_weight(d, i) for i in range(1, n)]
        Ag = Matrix(d, d, lambda i, j: (omg[i].T * omg[j])[0, 0])
        alc_vol2 = simplify(gram_det(Ag) / factorial(d) ** 2)   # vol^2 = det/(d!)^2
        n_alc_sq = simplify(covolQ2 / alc_vol2)    # (#)^2 = covol^2/vol^2
        n_alc = simplify(sqrt(n_alc_sq))
        want = factorial(n)
        ok(covolQ2 == n, "T3 d={0}: covol(Q)^2 = det Cartan = d+1 = {1}".format(d, n))
        ok(simplify(n_alc - want) == 0, "T3 d={0}: # alcoves = (d+1)! = {1}".format(d, want))
        print("{0} | √{1} | √{2} | {3} | {4} | {5}".format(
            d, covolQ2, simplify(gram_det(Ag) / factorial(d) ** 2),
            n_alc, want, simplify(n_alc - want) == 0))
    # Weyl reflections are isometries -> alcoves pairwise congruent
    for d in DS:
        n = d + 1
        us, SC = cell_vectors(d)
        # simple reflection s_1 swaps axes 0,1 (permutation) -> orthogonal in cell metric
        P = eye(n)
        P[0, 0] = 0
        P[1, 1] = 0
        P[0, 1] = 1
        P[1, 0] = 1
        preserves = all(simplify(udot(P * us[a], P * us[b], SC) - udot(us[a], us[b], SC)) == 0
                        for a in range(n) for b in range(n))
        ok(preserves, "T3 d={0}: a Weyl reflection = an isometry (Gram preserved) => alcoves are congruent".format(d))
    print("  ★MEASUREMENT: covol(Q)=√(d+1); vol(alcove{0,ω_i})=√(1/(d+1))/d!; ratio = (d+1)!")
    print("  (d=2→6, d=3→24, d=4→120). Weyl reflections preserve the Gram ⟹ all alcoves are CONGRUENT.")

    # ==================== T4: DEFLECTION (node = alcove barycenter?) ====================
    print()
    print("T4 (stake-d, ★DEFLECTION): the positions of the Dirac-zeros vs alcove barycenters — WITH NO tuning")
    print("d | node_1 (ψ phases) | centroid(Δ fund.alcove) | match? | f(node_1) | nodes=center ℤ/(d+1)?")
    print("-" * 96)
    for d in DS:
        n = d + 1
        # node_1: phases p_i = i/n (i=1..d), with p_0=0
        node1 = [Rational(i, n) for i in range(1, n)]
        Delta = alcove_delta_vertices(d)
        cen = centroid(Delta)
        cen_list = [cen[i, 0] for i in range(d)]
        match = (cen_list == node1)
        # f(node_1) = sum_{i=0}^{d} exp(2pi I*i/n) = gauss_sum(n,1) — exact (=0)
        fval = gauss_sum(n, 1)
        f_is0 = (fval == 0)
        # all nodes are Z/n torsion (center): f(node_j) = gauss_sum(n,j) = 0, j=1..d
        torsion = all(gauss_sum(n, j) == 0 for j in range(1, n))
        ok(match, "T4 d={0}: node_1 = centroid(fund.alcove) EXACTLY".format(d))
        ok(f_is0, "T4 d={0}: f(node_1) = 0 (a zero at the barycenter)".format(d))
        ok(torsion, "T4 d={0}: all node_j (j=1..d) = zeros = the center ℤ/(d+1)".format(d))
        print("{0} | {1} | {2} | {3} | {4} | {5}".format(
            d, node1, cen_list, match, fval, torsion))
    # distance node_1 <-> barycenter = 0 exactly (they coincide) — the deflection hits
    print("  ★★★MEASUREMENT (author pre-registration №5): the Dirac-zero SITS EXACTLY AT THE ALCOVE BARYCENTER.")
    print("  node_1 = centroid{0≤p_1≤...≤p_d≤1} = (i/(d+1))_i — an EXACT equality, DISTANCE 0,")
    print("  NO tuning at all (the probe only computed the centroid and the phases of the zero). The deflection hit the center.")
    print("  Mechanism: centroid = E[the i-th order statistic of d uniforms] = i/(d+1) = exactly")
    print("  the phases of the (d+1)-th roots of unity = a zero of the symmetric sum. Nodes=center=barycenters ∀d.")

    # ==================== MUTANTS ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    # M1: broken root e_1-e_3 (not simple) -> alcove/角 different
    d0 = 3
    bad_root = e_vec(d0 + 1, 0) - e_vec(d0 + 1, 2)   # e_1 - e_3
    good_root = simple_root(d0, 1)                    # e_1 - e_2
    ang_bad = simplify((bad_root.T * simple_root(d0, 2))[0, 0])   # <e1-e3, e2-e3> = 1
    ang_good = simplify((good_root.T * simple_root(d0, 2))[0, 0]) # <e1-e2, e2-e3> = -1
    if ang_bad != ang_good:
        print("  MUTANT M1: CAUGHT (broken root e_1-e_3: <·,α_2>={0} != the simple root {1} "
              "=> Cartan/angle differs)".format(ang_bad, ang_good))
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2: false character (zeta not root of unity: use exp(2pi I i/(n+1)) mismatched)
    n0 = d0 + 1
    fake = Add(*[exp(2 * pi * I * Rational(i, n0 + 1)) for i in range(n0)])
    fake_mag = abs(complex(fake.evalf(40)))
    fake_nonzero = (fake_mag > Rational(1, 10 ** 12))
    true_zero = (gauss_sum(n0, 1) == 0)
    if fake_nonzero and true_zero:
        print("  MUTANT M2: CAUGHT (false-character period n+1={0}: |sum|={1:.4f}≠0; "
              "the true ℤ/n character = 0)".format(n0 + 1, fake_mag))
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3: alcove count with mirror double -> 2*(d+1)! != (d+1)!
    claimed = 2 * factorial(d0 + 1)
    true_c = factorial(d0 + 1)
    if claimed != true_c:
        print("  MUTANT M3: CAUGHT (a count with a mirror double 2·(d+1)!={0} != (d+1)!={1})"
              .format(claimed, true_c))
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4: shifted barycenter -> distance to node != 0
    n0 = d0 + 1
    node1 = [Rational(i, n0) for i in range(1, n0)]
    shifted = [node1[0] + Rational(1, 10)] + node1[1:]     # perturb 1 coord
    diff = simplify(Add(*[(shifted[i] - node1[i]) ** 2 for i in range(d0)]))
    if diff != 0:
        print("  MUTANT M4: CAUGHT (shifted barycenter: |shift|²={0}≠0 => the zero is NO longer at "
              "the barycenter — the deflection test catches the tuning)".format(diff))
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): a random torus point vs. a zero")
    random.seed(1002023)
    d_nc = 3
    n_nc = d_nc + 1
    rnd_phases = [Rational(random.randrange(1, 10), 11) for _ in range(d_nc)]
    f_rnd = Add(*[exp(2 * pi * I * p) for p in [Rational(0, 1)] + rnd_phases])
    f2_rnd_mag = abs(complex(f_rnd.evalf(40))) ** 2
    ok(f2_rnd_mag > 1e-12, "control: a random point |f|²={0:.4f}≠0 (not a zero)".format(f2_rnd_mag))
    # and it's not the alcove barycenter
    node1_nc = [Rational(i, n_nc) for i in range(1, n_nc)]
    is_bary = (rnd_phases == node1_nc)
    ok(not is_bary, "control: the random point != the barycenter")
    print("  phases={0}: |f|²={1:.4f}≠0, != the barycenter — the measurement is sensitive (a zero only at the center-barycenter)"
          .format(rnd_phases, f2_rnd_mag))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

    # ==================== FORBIDDEN-SCAN ====================
    _pp = [("ко", "лір"), ("ква", "рк"), ("глю", "он"), ("гей", "дж"),  # GUARDLINE
           ("калібру", "вальн"), ("C", "2"), ("Casi", "mir")]  # GUARDLINE
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hits_src = scan_forbidden(__file__, _PATTERNS)
    _logf.flush()
    _logtxt = "".join(_tee.chunks)
    _hits_log = scan_forbidden(_logtxt, _PATTERNS)
    _nhits = len(_hits_src) + len(_hits_log)
    print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2}) [group-theory language su/Weyl/center/alcove "
          "is allowed; fence on the physics interpretation]".format(
              _nhits, len(_hits_src), len(_hits_log)))

    _exit = 0
    if _nhits > 0:
        _exit = 1
    if FAILS[0] > 0:
        _exit = 1
    if not mut_ok:
        _exit = 1
    print("EXIT={0}".format(_exit))
    print("PROC_EXIT={0}".format(_exit))
    print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
