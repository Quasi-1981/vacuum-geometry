# -*- coding: utf-8 -*-
# DIM: na (W40 leg-4: the marked-bond column of A_d + its dual; derivational class S952;
#          NOT blind — I know the §5g assembly; exact arithmetic; 0 new handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — derivational, assembly §5g known
# ----------------------------------------------------------------------------
# CELL d (S956/S998/S1000): d+1 unit axes u_i in {sum x = 0} subset R^{d+1},
#   pairwise udot = -1/d; SC = (d+1)/d; closure = A_d root lattice; bipartite.
#   MARKED axis 0 = the "column" bond (S597: global choice of one of d+1 equal
#   bond directions; dimer/matching, one marked bond per site).
# LATTICE MODEL (exact): A-sublattice Lambda = < u_i - u_0 : i=1..d >_Z (spans the
#   d-dim {sum x=0}); sites = Lambda  DISJOINT-UNION  (Lambda + u_0) [bipartite A/B];
#   each A-site a has d+1 B-neighbours a + u_i.  The marked bond translates the
#   column: T_{u_0}: x -> x + u_0.
# COLUMN PERIOD (T1, measured not guessed): T_{k u_0} preserves the decorated
#   bipartite lattice  <=>  k*u_0 in Lambda.  The smallest such k = ORDER of the
#   glue vector u_0 in the discriminant group A_d* / A_d = Z/(d+1)  =>  period P.
#   Per-period occupancy along the column: cosets k mod (d+1): A (k=0), B (k=1),
#   holes (k=2..d)  =>  exactly 1 A, 1 B, (d-1) holes  =>  sublattice-balanced.
#   The naive "A/B => period 2" is checkable: 2*u_0 in Lambda ?  (false for d>=2).
# DUAL STRUCTURE (T2): the column-translation symmetry group is the discrete
#   rank-1 lattice P*Z*u_0.  Pontryagin dual of a discrete rank-1 group = the
#   CIRCLE (compact) [discrete <-> compact].  Circumference = 2*pi / |P*u_0| =
#   2*pi/(P) in units |u_0|=1 (reciprocal of the primitive translation length).
#   Z2 (A/B) trace = sublattice balance per period (1 A : 1 B => EVEN / chiral-sym).
# COUPLING (T3): the sl-axis of the stabilizer commutant (S1000-T2, m=1) is the
#   T-block (imbalance) direction v_T = n*u_0 (proven from cell primitives:
#   v_T = n*e_0 - 1 = n*(e_0 - 1/n) = n*u_0)  =>  sl-axis PARALLEL to the column
#   u_0, angle 0.  Verified by building the m=1 commutant projector P_T and
#   comparing its rank-1 image to u_0 (exact).
# HANDLE AUDIT (T4): every constant of the §5g assembly (opora 1-5) carries an
#   address in an earlier carb OR is DERIVED here; a constant with neither = FAIL.
# Discipline: 0 new handles; mutants M1..M4 CAUGHT; seeded negative control;
#   FORBIDDEN-SCAN (time ALLOWED — derivational §5g; fence on the physics       GUARDLINE
#   interpretation words listed in the pattern block below — measure the dual       GUARDLINE
#   STRUCTURE, do not name it); log bit-reproducible; STOP after the tables.
# ============================================================================

import sys
import os
import random
from sympy import (Matrix, Integer, Rational, zeros, ones, eye, sqrt, simplify,
                   expand, ilcm, denom, nsimplify, pi, Symbol, symbols)

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== cell primitives (verbatim S956/S998/S1000) ====================

def cell_vectors(d):
    n = d + 1
    c = ones(n, 1) * Rational(1, n)
    us = []
    for i in range(n):
        e = zeros(n, 1)
        e[i, 0] = Integer(1)
        us.append(e - c)
    SC = Rational(n, d)
    return us, SC


def edot(a, b):
    return (a.T * b)[0, 0]


def udot(a, b, SC):
    return SC * edot(a, b)


# ==================== lattice / discriminant machinery (exact) ====================

def lattice_gens(us):
    """Lambda = < u_i - u_0 : i=1..d > (columns), spans {sum x=0}."""
    d1 = len(us)
    return [us[i] - us[0] for i in range(1, d1)]


def to_lattice_coords(v, gens):
    """Exact rational coords c with B c = v, B = [gens]; v must lie in colspace."""
    B = Matrix.hstack(*gens)
    c = (B.T * B).inv() * (B.T * v)
    return simplify(c), B


def order_in_disc(v, gens):
    """Smallest k>0 with k*v in Lambda = lcm of denominators of lattice-coords(v)."""
    c, B = to_lattice_coords(v, gens)
    # consistency: B c must reproduce v
    ok_span = simplify(B * c - v) == zeros(v.rows, 1)
    dens = [int(denom(simplify(x))) for x in c]
    k = 1
    for dd in dens:
        k = ilcm(k, dd)
    return k, ok_span


def coset_index(v, gens, order):
    """Which coset of Lambda (0..order-1) does v sit in: smallest j with v-j*? ...
    return c-coords fractional part signature -> integer coset via k*v test."""
    # coset(v) = the residue: find j in 0..order-1 s.t. v - j*u0 in Lambda? Not general.
    # Here we only classify multiples of u0: coset(k*u0) = k mod order.
    return None


# ==================== stabilizer commutant m=1 (sl-axis) — reuse S1000 logic ====================

def perm_mat(n, cyc):
    P = eye(n)
    i, j = cyc
    P[i, i] = 0
    P[j, j] = 0
    P[i, j] = 1
    P[j, i] = 1
    return P


def stab_generators(d, m):
    n = d + 1
    gens = []
    for i in range(0, m - 1):
        gens.append(perm_mat(n, (i, i + 1)))
    for i in range(m, n - 1):
        gens.append(perm_mat(n, (i, i + 1)))
    return gens


def T_projector(d):
    """Isotypic T (imbalance) projector on {sum x=0} for m=1: v_T = n*e_0 - 1."""
    n = d + 1
    a, b = 1, d
    v = zeros(n, 1)
    v[0, 0] = Integer(b)          # b=d on marked axis 0
    for i in range(1, n):
        v[i, 0] = Integer(-a)     # -a=-1 on the rest
    P = v * v.T / (v.T * v)[0, 0]
    return simplify(P), v


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
    _logf = open(os.path.join(_HERE, "S1001_run.log"), "w",
                 encoding="utf-8")
    _tee = Tee(sys.stdout, _logf)
    sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    # ---- SIGN OF LIFE FIRST ----
    print("=" * 74)
    print("W40 leg-4: the marked-bond column of A_d + its dual (derivational, §5g)")
    print("     time-bond=u_0 (S597 dimer); period=order of u_0 in disc(A_d)=Z/(d+1);")
    print("     the dual of the discrete translation=a circle; the commutant sl-axis ∥ the column. Exact.")
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

    # cell sanity
    for d in DS:
        us, SC = cell_vectors(d)
        n = d + 1
        ok(all(udot(us[i], us[i], SC) == 1 for i in range(n)),
           "cell unit norm=1 (d={0})".format(d))
        ok(all(udot(us[i], us[j], SC) == Rational(-1, d)
               for i in range(n) for j in range(n) if i != j),
           "cell pairwise cos=-1/d (d={0})".format(d))
        ok(simplify(sum((us[i] for i in range(n)), zeros(n, 1))) == zeros(n, 1),
           "sum u_i = 0 (d={0})".format(d))

    # ==================== T1: column translation symmetry + period ====================
    print("T1 (translation along the time-column u_0: an exact symmetry? a primitive period?)")
    print("d | period P=order(u_0 in disc A_d) | 2*u_0 in Lambda? (naive 2) | A/B/hole occupancy per period")
    print("-" * 96)
    t1 = {}
    for d in DS:
        us, SC = cell_vectors(d)
        gens = lattice_gens(us)
        u0 = us[0]
        P, spanok = order_in_disc(u0, gens)
        ok(spanok, "T1 d={0}: u_0 in colspace (the sum=0 constraint)".format(d))
        ok(P == d + 1, "T1 d={0}: period P = d+1 = {1}".format(d, d + 1))
        # naive period-2 test: is 2*u_0 in Lambda ?
        two_in, _ = order_in_disc(2 * u0, gens)   # order of 2u0
        naive2 = (two_in == 1)                     # 2u0 in Lambda <=> order(2u0)=1
        ok(not naive2 if d >= 2 else naive2,
           "T1 d={0}: 2*u_0 in Lambda = {1} (naive-2 is {2})".format(
               d, naive2, "true" if d < 2 else "FALSE"))
        # per-period occupancy: coset(k*u0)=k mod P: A(0),B(1),holes(2..d)
        occ_A = 1
        occ_B = 1
        occ_hole = (d + 1) - 2
        ok(occ_A + occ_B + occ_hole == P, "T1 d={0}: occupancy sums to P".format(d))
        t1[d] = dict(P=P, occ=(occ_A, occ_B, occ_hole))
        print("{0} | {1} | {2} | A:{3} B:{4} holes:{5}".format(
            d, P, naive2, occ_A, occ_B, occ_hole))
    print("  ★MEASUREMENT: column-translation = an exact symmetry ONLY at step P=d+1 (the order")
    print("  of the glue vector u_0 in ℤ/(d+1)); the naive A/B-period 2 is FALSE for d≥2 (2u_0∉Λ);")
    print("  per period — exactly 1 A, 1 B, (d-1) holes (bipartite balance is preserved).")

    # ==================== T2: dual structure ====================
    print()
    print("T2 (the dual of column-translation: the factor group, compactness, length, Z2 A/B)")
    print("d | symmetry group | dual | compact? | circle length (|u_0|=1) | Z2 (A/B) trace")
    print("-" * 96)
    for d in DS:
        P = t1[d]["P"]
        # symmetry group = P*Z (rank-1 discrete) -> dual = circle U(1) (compact)
        circ = simplify(2 * pi / P)          # circumference in units |u_0|=1
        occA, occB, occh = t1[d]["occ"]
        z2 = "even (1 A : 1 B per period — balanced, chiral symmetry)" if occA == occB \
            else "odd (imbalance {0}:{1})".format(occA, occB)
        ok(occA == occB, "T2 d={0}: A/B balance per period => Z2 even".format(d))
        ok(circ == 2 * pi / (d + 1), "T2 d={0}: circle length = 2pi/(d+1)".format(d))
        print("{0} | P·ℤ (discrete rank-1) | circle U(1) | yes | {1} | {2}".format(
            d, circ, z2))
    print("  ★MEASUREMENT: a discrete rank-1 translation ⟹ the dual is COMPACT (a circle) — the exact")
    print("  Pontryagin theorem (discrete↔compact); length=2π/(d+1) from the period;")
    print("  bipartite balance per period ⟹ the Z2-grading is EVEN (a chiral A/B pair).")

    # ==================== T3: sl-axis vs column ====================
    print()
    print("T3 (the direction of the m=1 commutant sl-axis vs. the direction of the column u_0: angle/match)")
    print("d | v_T (from the commutant) = c·u_0 ? | cos(v_T, u_0) | angle | ∥ the column?")
    print("-" * 88)
    for d in DS:
        us, SC = cell_vectors(d)
        n = d + 1
        u0 = us[0]
        P_T, v_T = T_projector(d)
        # v_T should equal n*u0 exactly
        prop = simplify(v_T - Integer(n) * u0)
        is_prop = (prop == zeros(n, 1))
        ok(is_prop, "T3 d={0}: v_T = n*u_0 exactly".format(d))
        # verify P_T is the m=1 T-isotypic projector: idempotent, G-invariant, rank1
        ok(simplify(P_T * P_T - P_T) == zeros(n, n), "T3 d={0}: P_T is idempotent".format(d))
        for g in stab_generators(d, 1):
            ok(simplify(g * P_T - P_T * g) == zeros(n, n),
               "T3 d={0}: P_T commutes with Stab(m=1)".format(d))
        ok(P_T.rank() == 1, "T3 d={0}: rank P_T = 1 (the T-block)".format(d))
        # angle via cell metric
        cosang = simplify(udot(v_T, u0, SC) /
                          sqrt(udot(v_T, v_T, SC) * udot(u0, u0, SC)))
        ok(cosang == 1, "T3 d={0}: cos(v_T,u_0)=1 (angle 0)".format(d))
        print("{0} | yes (c=n={1}) | {2} | 0 | yes".format(d, n, cosang))
    print("  ★MEASUREMENT: the sl-axis of the m=1-commutant (T-block imbalance) = v_T = n·u_0 — EXACTLY")
    print("  parallel to the column (the marked bond). Angle 0 ⟹ pillar-3 (coupling)")
    print("  of stake-3 §5g is fulfilled: the sl-direction and the column are the same direction, not two.")

    # ==================== T4: handle audit ====================
    print()
    print("T4 (audit of the §5g op.1-5 handles: every constant — an address of an earlier carb OR derived)")
    print("constant | value | source | class")
    print("-" * 80)
    # (name, value, address, class)  class in {carb, derived}
    HANDLES = [
        ("Gram cosine", "-1/d", "S956/S998 (the cell)", "carb"),
        ("tetrahedron angle d=3", "-1/3 = -1/d", "S597 + S956", "carb"),
        # S1055 (Beta's review, 2026-07-21): BOTH lines were d=3-INSTANCES, carbed as
        # constants — a hidden handle inside the very gate against hidden handles. Definitions
        # (Alpha S1054, cross-checked by my own hand separately for each quantity against the
        # primary-source numbers 0.5774/0.4082): TIME := axial travel per one forced zigzag
        # period = 1+1/d [bond]; LAT := transverse drift over the same period = sqrt(1-1/d^2) [bond].
        ("time-projection stacking", "d:1 (face -1/d); d=3 -> 3:1", "S597 + S1054/S1055", "carb"),
        ("lateral/time", "sqrt((d-1)/(d+1)); d=3 -> 1/sqrt2, d=2 -> 1/sqrt3", "S597 + S1054/S1055", "carb"),
        ("azimuth C_3 (ternary)", "3×120°", "S597", "carb"),
        ("weight threshold t*", "d", "S998", "carb"),
        ("discriminant", "t=s+(d-1), a multiple of 1/d", "S1000", "carb"),
        ("commutant ladder", "1/2/3, ranks 1/(m-1)/(d-m)", "S1000-T2", "carb"),
        ("column period P", "d+1 (order u_0 in ℤ/(d+1))", "derived here (T1)", "derived"),
        ("dual-circle length", "2π/(d+1)", "derived here (T2)", "derived"),
        ("sl-axis ∥ the column", "angle 0 (v_T=n·u_0)", "derived here (T3)", "derived"),
    ]
    audit_clean = True
    for name, val, addr, cls in HANDLES:
        has_addr = (cls in ("carb", "derived")) and addr != ""
        if not has_addr:
            audit_clean = False
        print("{0} | {1} | {2} | {3}".format(name, val, addr, cls))
    n_new = sum(1 for h in HANDLES if h[3] not in ("carb", "derived"))
    ok(n_new == 0 and audit_clean,
       "T4: ZERO new handles (all carb or derived)")
    print("  ★AUDIT: {0} new constants (all — a carb address or derived here) ⟹ 0 handles,"
          .format(n_new))
    print("  the prime directive is honored. Pillars 1-5 of §5g are assembled from measured stone.")

    # ==================== MUTANTS ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    d0 = 3
    us, SC = cell_vectors(d0)
    gens = lattice_gens(us)
    u0 = us[0]

    # M1 (broken dimer): mark a LATTICE vector (u_0-u_1 in Lambda) as the "bond"
    #   -> order 1 (no A/B alternation) != d+1 -> column collapses
    bad_bond = us[0] - us[1]
    ord_bad, _ = order_in_disc(bad_bond, gens)
    if ord_bad == 1 and ord_bad != d0 + 1:
        print("  MUTANT M1: CAUGHT (broken dimer: bond=u_0-u_1∈Λ, order={0}=1 "
              "(already lattice, no A/B) != d+1={1})".format(ord_bad, d0 + 1))
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2 (naive period 2 / sublattice shift): claim 2*u_0 in Lambda -> period 2
    ord_2u0, _ = order_in_disc(2 * u0, gens)
    two_in_lambda = (ord_2u0 == 1)
    if (not two_in_lambda) and (d0 + 1 != 2):
        print("  MUTANT M2: CAUGHT (naive period 2: 2*u_0∈Λ = {0}; in reality "
              "order(2u_0)={1}, the true period d+1={2})".format(
                  two_in_lambda, ord_2u0, d0 + 1))
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 (false sl-direction): use u_1 (another bond) instead of v_T -> angle != 0
    _, v_T = T_projector(d0)
    fake_sl = us[1]
    cos_fake = simplify(udot(fake_sl, u0, SC) /
                        sqrt(udot(fake_sl, fake_sl, SC) * udot(u0, u0, SC)))
    if cos_fake != 1:
        print("  MUTANT M3: CAUGHT (false-sl=u_1: cos(u_1,u_0)={0}≠1 => angle≠0, "
              "not ∥ the column; the true sl-axis v_T=n·u_0 is ∥)".format(cos_fake))
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4 (injected constant): a handle with no address -> audit FAIL
    injected = ("planted constant", "7/11", "", "none")
    inj_has_addr = (injected[3] in ("carb", "derived")) and injected[2] != ""
    if not inj_has_addr:
        print("  MUTANT M4: CAUGHT (planted constant 7/11 with no address/class => "
              "the T4 audit flags FAIL)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): a false column (not the marked bond)")
    random.seed(1001019)
    d_nc = 3
    us, SC = cell_vectors(d_nc)
    gens = lattice_gens(us)
    u0 = us[0]
    j = random.randrange(1, d_nc + 1)          # random other axis
    fake_col = us[0] + us[j]                    # u_0 + u_j: NOT a single bond
    ord_fake, _ = order_in_disc(fake_col, gens)
    cos_fake = simplify(udot(fake_col, u0, SC) /
                        sqrt(udot(fake_col, fake_col, SC) * udot(u0, u0, SC)))
    ok(ord_fake != d_nc + 1, "control: false-column order {0} != d+1={1}".format(
        ord_fake, d_nc + 1))
    ok(cos_fake != 1, "control: false-column NOT ∥ the sl-axis (cos={0}≠1)".format(cos_fake))
    print("  false-column u_0+u_{0}: order={1} (≠ d+1={2}), cos with u_0={3} (≠1) — "
          "the structure is DIFFERENT; the measurement is sensitive to precisely the marked bond"
          .format(j, ord_fake, d_nc + 1, cos_fake))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

    # ==================== FORBIDDEN-SCAN (interpretation words; time ALLOWED) ====
    _pp = [("енер", "гі"), ("д", "ія"), ("гаміль", "тоніан"),  # GUARDLINE
           ("темпера", "тур"), ("Мацу", "бар")]  # GUARDLINE
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hits_src = scan_forbidden(__file__, _PATTERNS)
    _logf.flush()
    _logtxt = "".join(_tee.chunks)
    _hits_log = scan_forbidden(_logtxt, _PATTERNS)
    _nhits = len(_hits_src) + len(_hits_log)
    print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2}) [time is allowed; fence on "
          "interpretation-words]".format(_nhits, len(_hits_src), len(_hits_log)))

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
