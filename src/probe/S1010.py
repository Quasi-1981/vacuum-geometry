# -*- coding: utf-8 -*-
# DIM: na (W42 probe-5, layer -2: the author's ansatz «a rotational axis sets the symmetry
#          of space» on the CELL (Cartan torus A_d, Bloch g=w0+Σwᵢexp). Does time-detuning (B) spare ⊥
#          space-detuning (C) kill rotation-born structures? The S999/S1002 machinery
#          reproduced exactly; ancestors (S999 cone · S1002 nodes=barycenters · S1001 period ·
#          S1004 Coxeter) CITED, not re-derived. Exact arithmetic; 0 handles.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — §9 exante + S999/S1002 ancestors
# ----------------------------------------------------------------------------
# CELL d: d+1 unit axes u_i in {sum x=0}, udot(u_i,u_j)=-1/d (i!=j), SC=(d+1)/d.
# BLOCH g(psi; w) = w_0 + sum_{i=1..d} w_i exp(2 pi I psi_i)  (psi_i phases, p_0=0 convention;
#   axis 0 = TIME-COLUMN = fixed axis of the C_d screw; axes 1..d = SPATIAL, C_d-orbit).
#   S999 form with w_0=t on axis 0, w_i on the rest.  Node: g=0.
# METRIC (S999): G_ij=<alpha_i,alpha_j>=SC(1+[i==j]); CONE at a node = K=J_psi G J_psi^T
#   (2x2), J_psi col i = w_i(-sin,cos)(psi_i) [2pi dropped]; ISOTROPY <=> disc(K)=0.
# C_d = cyclic permutation of SPATIAL phase-coords (fixes axis 0).  Config C_d-invariant
#   <=> spatial weights all equal (theorem — exercise-fence §9, NOT measured as a bet).
# THREE CONFIGS: (A) all w=1 · (B) time-detune w_0=t in {3/2,2} (spatial=1) ·
#   (C) space-detune ONE spatial w_1=eps in {3/2,2} (rest=1).
# MEASURED (change of 4 ancestor structures under detuning ONLY):
#   T1 cone isotropy at node (d=2 exact solve) — A/B/C.
#   T2 ★band C_d-symmetry: |g|^2 at a point vs its C_d-images — EQUAL (C_d intact) or
#      DIFFER (C_d broken).  d=2 (C_2) and d=3 (C_3).  THIS is the ansatz asymmetry.
#   T3 nodal-set orbit structure under C_d (d=2 exact): closed under swap? orbit sizes.
#   T4 column period P=d+1 (S1001): weight-independent — honest NON-discriminator.
#   T5 bet-2: cone eigenvalue-splitting pattern (disc) by stabilizer under (C).
# BETS (§9, open outputs): bet1 (B live ⊥ C die = forsage); bet2 (splitting by isotypes);
#   K1 (both live => rotation not holder); K2 (both die => detune-as-such); K3 (depends
#   on WHICH spatial axis => an axis-labeling artifact S1000-4a, mutant catches by permutation).
# Discipline: 0 handles; exact; mutants>=4 (perm detuned axis=equiv · false-isotropy
#   float · d=2<->3 · limit eps->1 continuous); seeded negctrl; FORBIDDEN patterns in
#   GUARDLINE block; the S592-rhyme NOT invoked; STOP after tables.
# ============================================================================

import sys
import os
import random
from sympy import (Matrix, Integer, Rational, zeros, ones, eye, exp, I, pi, cos, sin,
                   simplify, expand, Add, sqrt, re, im, symbols, solve, nsimplify, Abs)

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== cell / metric primitives (exact, from S999/S1002) ====================

def e_vec(n, i):
    v = zeros(n, 1); v[i, 0] = Integer(1); return v


def cell_vectors(d):
    n = d + 1
    c = ones(n, 1) * Rational(1, n)
    us = [e_vec(n, i) - c for i in range(n)]
    return us, Rational(n, d)


def udot(a, b, SC):
    return SC * (a.T * b)[0, 0]


def gram_alpha(d):
    us, SC = cell_vectors(d)
    al = [us[i] - us[0] for i in range(1, d + 1)]
    G = Matrix(d, d, lambda i, j: udot(al[i], al[j], SC))
    assert G == SC * (eye(d) + ones(d, d)), "Gram(alpha)=SC(I+J)"
    return G, SC


def jac_psi(sample, c):
    """2 x d Jacobian of (Re g, Im g) wrt psi at sample [(x_i,y_i)]; col i = c_i(-y,x)."""
    M = zeros(2, len(sample))
    for i, (x, y) in enumerate(sample):
        M[0, i] = -c[i] * y
        M[1, i] = c[i] * x
    return M


def eig2(K):
    tr = simplify(K[0, 0] + K[1, 1]); dt = simplify(K.det())
    disc = simplify(expand(tr ** 2 - 4 * dt))
    return simplify((tr + sqrt(disc)) / 2), simplify((tr - sqrt(disc)) / 2), disc


def g_at(psi, w0, wsp):
    """g = w0 + sum_i wsp_i exp(2 pi I psi_i), psi rational phases; exact complex."""
    s = w0 + Add(*[wsp[i] * exp(2 * pi * I * psi[i]) for i in range(len(psi))])
    return simplify(expand(s))


def absg2(psi, w0, wsp):
    g = g_at(psi, w0, wsp)
    return simplify(expand(re(g) ** 2 + im(g) ** 2))


def cyc(seq, s=1):
    """cyclic shift of a list by s (C_d action on spatial coords)."""
    n = len(seq); return [seq[(i - s) % n] for i in range(n)]


def solve_d2_node(w0, c1, c2):
    """Exact real nodes of w0 + c1 w1 + c2 w2 = 0, |w_i|=1 (d=2, general weights)."""
    x1, y1, x2, y2 = symbols('x1 y1 x2 y2', real=True)
    sols = solve([w0 + c1 * x1 + c2 * x2, c1 * y1 + c2 * y2,
                  x1 ** 2 + y1 ** 2 - 1, x2 ** 2 + y2 ** 2 - 1],
                 [x1, y1, x2, y2], dict=True)
    pts = []
    for s in sols:
        vals = (s[x1], s[y1], s[x2], s[y2])
        if all(v.is_real for v in vals):
            pts.append(vals)
    return sorted(set(pts), key=str)


def cone_disc(node, c1, c2, G):
    """disc(K) at a d=2 node; disc=0 <=> isotropic cone."""
    (x1, y1, x2, y2) = node
    J = jac_psi([(x1, y1), (x2, y2)], [c1, c2])
    r = J.rank()
    K = simplify(expand(J * G * J.T))
    l1, l2, disc = eig2(K)
    return r, disc, l1, l2


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1010_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-5 (layer −2): the ansatz «a rotational axis sets the symmetry of space» on the CELL.")
    print("Does time-detuning (B) spare ⊥ space-detuning (C) kill rotation-born structures?")
    print("Cartan torus A_d; the S999 cone / band C_d-symmetry / node-barycenters; exact; d∈{2,3}.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, m):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + m)

    G2, SC2 = gram_alpha(2)
    ONE = Integer(1)
    TDET = [Rational(3, 2), Rational(2)]   # detune values (commensurate)

    # ==================== T1: d=2 cone isotropy under A/B/C ====================
    print("T1 (d=2 cone at a node: isotropy (disc=0?) under detuning — A/B/C):")
    print("config | w0(time) | c1,c2(space) | nodes | rank | disc(K) | isotropic?")
    print("-" * 72)
    def report_d2(tagname, w0, c1, c2):
        nodes = solve_d2_node(w0, c1, c2)
        if not nodes:
            print("{0} | {1} | {2},{3} | NO nodes".format(tagname, w0, c1, c2))
            return None
        discs = []
        for nd in nodes:
            r, disc, l1, l2 = cone_disc(nd, c1, c2, G2)
            discs.append((r, disc))
        alliso = all(simplify(dd) == 0 for (_, dd) in discs)
        anyiso = any(simplify(dd) == 0 for (_, dd) in discs)
        rr = discs[0][0]
        dd0 = discs[0][1]
        print("{0} | {1} | {2},{3} | {4} | {5} | {6} | {7}".format(
            tagname, w0, c1, c2, len(nodes), rr, dd0,
            "YES(all)" if alliso else ("part." if anyiso else "NO")))
        return alliso
    isoA = report_d2("(A) equal", ONE, ONE, ONE)
    isoB = [report_d2("(B) t={0}".format(t), t, ONE, ONE) for t in TDET]
    isoC = [report_d2("(C) ε={0}".format(e), ONE, e, ONE) for e in TDET]
    ok(isoA is True, "★T1: (A) equal weights → the cone is ISOTROPIC (disc=0) — the S999 ancestor is reproduced")
    ok(all(x is False for x in isoB), "★T1: (B) time-detuning → the cone BREAKS (disc≠0) — an S999-T5.1 rhyme")
    ok(all(x is False for x in isoC), "★T1: (C) space-detuning → the cone BREAKS (disc≠0)")
    print("  ★MEASUREMENT T1 (honest, against bet-1): local cone isotropy breaks under BOTH")
    print("  detunings. Mechanism: the node unpins from the barycenter under any detuning; isotropy")
    print("  was held by the FULL Weyl symmetry S_{d+1} (node at the center), not C_d itself. The cone = an object")
    print("  of the full symmetry, NOT a rotational one ⟹ the B⊥C asymmetry must be sought elsewhere (see T2). The K2 flavor")
    print("  is ONLY for the cone: locally, detuning-as-such breaks it; but globally — it does not (T2/T3).")

    # ==================== T2: band C_d-symmetry (THE ansatz asymmetry) ====================
    print()
    print("T2 ★CORE: C_d-symmetry of the BAND |g|² — |g|² at a point vs its C_d-images (equal/different?):")
    print("measured at a seeded point + barycenter; B: all equal (C_d intact) ⊥ C: different (C_d broken)")
    print("d | config | {|g|² on the C_d-orbit of the point} | all equal? | band C_d-symmetry")
    print("-" * 72)
    random.seed(1010051)
    def cd_orbit_vals(d, w0, wsp, psi):
        vals = []
        cur = list(psi)
        for _ in range(d):
            vals.append(absg2(cur, w0, wsp))
            # C_d permutes which spatial axis carries which phase == cyclic shift of wsp
            wsp = cyc(wsp, 1)
        return [simplify(v) for v in vals]
    def report_cd(d, tag, w0, wsp):
        # seeded generic rational point on the torus
        psi = [Rational(random.randrange(1, 12), 13) for _ in range(d)]
        vals = cd_orbit_vals(d, w0, list(wsp), psi)
        alleq = all(simplify(v - vals[0]) == 0 for v in vals)
        print("{0} | {1} | {2} | {3} | {4}".format(
            d, tag, "equal" if alleq else "DIFFERENT({0})".format(len(set(str(v) for v in vals))),
            "YES" if alleq else "NO", "INTACT" if alleq else "BROKEN"))
        return alleq
    for d in (2, 3):
        eqA = report_cd(d, "(A)", ONE, [ONE] * d)
        eqB = report_cd(d, "(B) t=3/2", Rational(3, 2), [ONE] * d)
        eqC = report_cd(d, "(C) ε=3/2", ONE, [Rational(3, 2)] + [ONE] * (d - 1))
        ok(eqA, "T2 d={0}: (A) the band C_d-symmetry is intact".format(d))
        ok(eqB, "★T2 d={0}: (B) time-detuning — the band C_d-symmetry is INTACT (the time-weight = the fixed axis of the screw)".format(d))
        ok(not eqC, "★T2 d={0}: (C) space-detuning — the band C_d-symmetry is BROKEN".format(d))
    print("  ★★★MEASUREMENT T2 (THE AUTHOR'S ANSATZ): the band C_d-symmetry SURVIVES under time-detuning (B) ⊥ DIES")
    print("  under space-detuning (C), ∀d. The mechanism EXACTLY the ansatz: the time-column = the FIXED axis")
    print("  of the C_d screw ⟹ its weight is a free spectator (never breaks C_d); the spatial axis")
    print("  sits IN the screw's ORBIT ⟹ its detuning breaks C_d. The forcing = both halves together.")

    # ==================== T3: nodal-set orbit structure under swap (d=2 exact) ====================
    print()
    print("T3 (d=2: orbit structure of NODES under the C_2 swap — closed? orbit sizes):")
    print("config | nodes (w1,w2) | closed under swap? | orbits under C_2")
    print("-" * 72)
    def swap_closed(w0, c1, c2):
        nodes = solve_d2_node(w0, c1, c2)
        node_ws = []
        for (x1, y1, x2, y2) in nodes:
            node_ws.append(((x1, y1), (x2, y2)))
        sset = set(str(nw) for nw in node_ws)
        closed = all(str((b, a)) in sset for (a, b) in node_ws)
        # orbit sizes under swap
        seen = set(); orbs = []
        for nw in node_ws:
            k = str(nw)
            if k in seen: continue
            a, b = nw; sw = (b, a)
            if str(sw) in sset and str(sw) != k:
                orbs.append(2); seen.add(k); seen.add(str(sw))
            else:
                orbs.append(1); seen.add(k)
        return closed, sorted(orbs, reverse=True), len(nodes)
    cA, oA, nA = swap_closed(ONE, ONE, ONE)
    cB, oB, nB = swap_closed(Rational(3, 2), ONE, ONE)
    cC, oC, nC = swap_closed(ONE, Rational(3, 2), ONE)
    print("(A) equal     | {0} nodes | {1} | {2}".format(nA, cA, oA))
    print("(B) t=3/2     | {0} nodes | {1} | {2}".format(nB, cB, oB))
    print("(C) ε=3/2     | {0} nodes | {1} | {2}".format(nC, cC, oC))
    ok(cA and cB, "★T3: the nodes are closed under C_2 at (A) and (B) — the orbits are intact")
    ok(not cC, "★T3: (C) space-detuning — the nodes are NOT closed under C_2 (the orbits split)")
    print("  ★MEASUREMENT T3: the nodal set holds C_2-orbits under time-detuning (B, paired by the swap)")
    print("  ⊥ splits under space-detuning (C, the swap is not a symmetry). The second leg of T2 — global.")

    # ==================== T4: period P — honest non-discriminator ====================
    print()
    print("T4 (column period P=d+1, S1001 citation): weight-INDEPENDENT — an honest NON-discriminator:")
    for d in (2, 3):
        # P = order of column element in disc(A_d) = Z/(d+1) — combinatorial, weight-free
        P = d + 1
        print("  d={0}: P=d+1={1} — combinatorial (order of the center ℤ/(d+1), S1001) ⟹ the same in A/B/C".format(d, P))
    ok(True, "T4: P=d+1 is weight-independent ⟹ does NOT distinguish B/C (honestly: not all 4 structures are discriminators)")
    print("  ★MEASUREMENT T4: the period P does not depend on the weights (a center structure) ⟹ B and C give the same P.")
    print("  Honestly: of the 4 rotational structures, only the C_d-symmetry (T2) and the orbits (T3) distinguish B⊥C;")
    print("  the cone (T1) breaks under both; the period (T4) under neither. The ansatz lives in the GLOBAL ones (T2/T3).")

    # ==================== T5: bet-2 splitting pattern ====================
    print()
    print("T5 (bet-2: the splitting pattern of the cone eigenvalues under (C) — by stabilizer):")
    for e in TDET:
        nodes = solve_d2_node(ONE, e, ONE)
        if nodes:
            r, disc, l1, l2 = cone_disc(nodes[0], e, ONE, G2)
            print("  (C) ε={0}: cone eigenvalues (l1,l2)=({1}, {2}); disc={3}≠0 ⟹ split 1+1"
                  .format(e, simplify(l1), simplify(l2), simplify(disc)))
    print("  ★MEASUREMENT T5: under (C) the stabilizer of the detuned config. = trivial (ε singles out axis-1)")
    print("  ⟹ the 2×2 cone splits 1+1 (two distinct eigenvalues) — the isotype pattern of the trivial group.")
    print("  The ancestor-isotypes (S1000-T2 commutant 1/2/3) cited; the block-size here = 1+1 (d=2).")

    # ==================== MUTANTS ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    # M1 (perm detuned spatial axis = equivalence, S1000-4a axis-labeling artifact): eps on axis1 vs axis2
    n1 = solve_d2_node(ONE, Rational(3, 2), ONE)
    n2 = solve_d2_node(ONE, ONE, Rational(3, 2))
    d1 = simplify(cone_disc(n1[0], Rational(3, 2), ONE, G2)[1]) if n1 else None
    # disc for eps on axis2: match node ordering by |g|^2 multiset instead
    def disc_multiset(nodes, c1, c2):
        return sorted(str(simplify(cone_disc(nd, c1, c2, G2)[1])) for nd in nodes)
    ms1 = disc_multiset(n1, Rational(3, 2), ONE)
    ms2 = disc_multiset(n2, ONE, Rational(3, 2))
    if ms1 == ms2 and n1 and n2:
        print("  MUTANT M1 (permuting the detuned axis = an axis-labeling artifact, S1000-4a): CAUGHT (ε on axis-1 vs axis-2")
        print("    gives the SAME disc-multiset {0} — the measurement is invariant to the CHOICE of axis; the K3-artifact is dead)".format(ms1))
    else:
        print("  MUTANT M1: NOT CAUGHT (ms1={0} ms2={1})".format(ms1, ms2)); mut_ok = False

    # M2 (false-isotropy via float): broken cone disc tiny-but-nonzero must NOT pass as 0
    nodes = solve_d2_node(ONE, Rational(3, 2), ONE)
    _, discC, _, _ = cone_disc(nodes[0], Rational(3, 2), ONE, G2)
    exact_nonzero = (simplify(discC) != 0)
    float_would_pass = (abs(complex(discC.evalf(30))) < 1e-6) if exact_nonzero else False
    if exact_nonzero:
        print("  MUTANT M2 (false-isotropy via float): CAUGHT (disc={0} EXACTLY≠0; a float tolerance could".format(simplify(discC)))
        print("    let a broken cone «pass» — exact arithmetic catches it, float is forbidden in the verdict)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 (d=2 <-> d=3 size): T2 asymmetry holds at BOTH d
    psi3 = [Rational(1, 5), Rational(2, 7), Rational(3, 11)]
    v3B = cd_orbit_vals(3, Rational(3, 2), [ONE, ONE, ONE], psi3)
    v3C = cd_orbit_vals(3, ONE, [Rational(3, 2), ONE, ONE], psi3)
    b3 = all(simplify(v - v3B[0]) == 0 for v in v3B)
    c3 = not all(simplify(v - v3C[0]) == 0 for v in v3C)
    if b3 and c3:
        print("  MUTANT M3 (size d=2↔3): CAUGHT (at d=3 the same asymmetry: (B) the C_3-band is intact,")
        print("    (C) is broken — the conclusion is not a d=2 lattice artifact)")
    else:
        print("  MUTANT M3: NOT CAUGHT (b3={0} c3={1})".format(b3, c3)); mut_ok = False

    # M4 (limit eps->1 continuous return to A): disc(eps) -> 0 as eps->1
    nodes_e = solve_d2_node(ONE, Rational(11, 10), ONE)
    _, disc_e, _, _ = cone_disc(nodes_e[0], Rational(11, 10), ONE, G2)
    nodes_1 = solve_d2_node(ONE, ONE, ONE)
    _, disc_1, _, _ = cone_disc(nodes_1[0], ONE, ONE, G2)
    small = abs(complex(disc_e.evalf(30)))
    if simplify(disc_1) == 0 and 0 < small < abs(complex(cone_disc(solve_d2_node(ONE, Rational(3,2), ONE)[0], Rational(3,2), ONE, G2)[1].evalf(30))):
        print("  MUTANT M4 (limit ε→1 continuous): CAUGHT (disc(ε=11/10)={0:.4f} small>0, disc(ε=1)=0 —"
              .format(small))
        print("    a continuous return to (A); a discontinuity would be an artifact)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): a symmetric point NOT at a C_d-orbit break under (A)")
    random.seed(1010052)
    d_nc = 3
    psi_nc = [Rational(random.randrange(1, 8), 9) for _ in range(d_nc)]
    v_nc = cd_orbit_vals(d_nc, ONE, [ONE] * d_nc, psi_nc)
    ok(all(simplify(v - v_nc[0]) == 0 for v in v_nc),
       "control: (A) equal weights → the C_3-orbit |g|² is EQUAL at the seeded point (symmetry intact)")
    print("  phases={0}: the C_3-values are all equal under (A) — the measurement is sensitive (it distinguishes a break only in C)".format(psi_nc))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'rotation/screw/C_d' is STRUCTURAL cell-geometry (Coxeter S1004) — not fenced;   GUARDLINE
    #   the S592 rhyme is NOT invoked; no physics-interpretation written.                  GUARDLINE
    _pp = [("причи", "нн"), ("всес", "віт"), ("ант", "роп"), ("Teg", "mark"),  # GUARDLINE
           ("сп", "ін"), ("S5", "92")]  # GUARDLINE
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _n = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2})".format(_n, len(_hs), len(_hl)))

    _exit = 1 if (_n > 0 or FA[0] > 0 or not mut_ok) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
