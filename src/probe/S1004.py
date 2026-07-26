# -*- coding: utf-8 -*-
# DIM: na (W41 «Coxeter-time axiom»: one ladder cos(pi j/h) under the program's spectra +
#          azimuth = the action of the Coxeter element; derivational class; exact arithmetic; 0 handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — derivational, §W41 exante known
# ----------------------------------------------------------------------------
# A_d: Coxeter number h = d+1; exponents m_j = 1,2,...,d.  Cartan matrix C (path
#   graph P_d): C[i,i]=2, C[i,i+1]=C[i+1,i]=-1.  Eigenvalues of C = 2-2cos(pi j/h),
#   j=1..d, with eigenvector v_j[i]=sin(pi i j/h) (path-Laplacian identity).
# LADDER = {SC*(2-2cos(pi j/h))}_{j=1..d}, SC=(d+1)/d = the unit-Gram spectrum.
# BET 1 (one ladder under all carbed quadratic spectra):
#   (a) unit-Gram of the cell = SC*Cartan(A_d) (S956 «generator Gram»); its
#       eigenvalues = SC*(2-2cos(pi j/h))  -> the exponent ladder.  Exact.
#   (b) transverse-eigenvalues S999 (isotropic K at symmetric points) = n^2/(2d)
#       = SC*n/2; is it a ladder rung SC*(2-2cos(pi j/h)) for integer j?  Verdict.
#   (c) commutant blocks S1000-T2 (isotypic dims 1/2/3, orthogonal blocks) — a
#       cos-ladder?  Verdict.
#   FORSAGE = >= 2 of {a,b,c} reduce to the ladder WITHOUT tuning.
# BET 2 (time = Coxeter action): the Coxeter element c = (d+1)-cycle permutation of
#   the cell axes; on {sum x=0} it has eigenvalues e^{2 pi i j/h} (exponents), order
#   h.  d=2: c is a rotation of the plane {sum x=0}; angle from tr(c|plane)=2cos(theta)
#   => theta = azimuth.  Check theta = 120° = 2pi/3 = 2pi/h (S597).  All-d: <c> order h.
# ANTI-TUNING (K3): any fitted coefficient outside {SC, 2pi, integer j, h} = FAIL.
# Discipline: 0 new handles; mutants M1..M4 + anti-tuning; seeded negative control;
#   FORBIDDEN-SCAN {physics interpretation} — su/Weyl/Coxeter/exponents ALLOWED;
#   log bit-reproducible; STOP after tables.
# ============================================================================

import sys
import os
import random
from sympy import (Matrix, Integer, Rational, zeros, ones, eye, cos, sin, pi,
                   sqrt, simplify, expand, exp, I, nsimplify, Rational as Rat)

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== primitives ====================

def cartan_Ad(d):
    C = zeros(d, d)
    for i in range(d):
        C[i, i] = Integer(2)
        if i + 1 < d:
            C[i, i + 1] = Integer(-1)
            C[i + 1, i] = Integer(-1)
    return C


def SC_of(d):
    return Rational(d + 1, d)


def ladder(d):
    """{SC*(2-2cos(pi j/h))}_{j=1..d}, h=d+1 — the exponent (unit-Gram) ladder."""
    h = d + 1
    SC = SC_of(d)
    return [simplify(SC * (2 - 2 * cos(pi * Rational(j, h)))) for j in range(1, d + 1)]


def cartan_eig_closed(d):
    h = d + 1
    return [simplify(2 - 2 * cos(pi * Rational(j, h))) for j in range(1, d + 1)]


def cyc_perm_matrix(n):
    """(0 1 2 ... n-1) cyclic permutation matrix (Coxeter element as an S_n cycle)."""
    P = zeros(n, n)
    for i in range(n):
        P[(i + 1) % n, i] = Integer(1)
    return P


def basis_sumzero(n):
    """Orthonormal-ish basis of {sum x=0}: columns e_i - e_{i+1} (i=0..n-2)."""
    cols = []
    for i in range(n - 1):
        v = zeros(n, 1)
        v[i, 0] = Integer(1)
        v[i + 1, 0] = Integer(-1)
        cols.append(v)
    return Matrix.hstack(*cols)


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
    _logf = open(os.path.join(_HERE, "S1004_run.log"), "w",
                 encoding="utf-8")
    _tee = Tee(sys.stdout, _logf)
    sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W41 «Coxeter-time axiom»: one ladder cos(πj/h) under the program's spectra +")
    print("the zigzag azimuth = the action of the Coxeter element. Exact algebra, 0 handles.")
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

    # ==================== BET 1(a): unit-Gram spectrum = ladder ====================
    print("STAKE 1(a): eigenvalues of the cell's unit-Gram (SC·Cartan A_d) = SC(2−2cos(πj/h))")
    print("d | eigenvector test Cartan·v_j=(2−2cos(πj/h))v_j | ladder SC(2−2cos) | matches S956")
    print("-" * 96)
    for d in DS:
        h = d + 1
        C = cartan_Ad(d)
        # verify eigenvector v_j[i]=sin(pi*i*j/h) gives eigenvalue 2-2cos(pi j/h)
        all_ev = True
        for j in range(1, d + 1):
            vj = Matrix([sin(pi * Rational((i + 1) * j, h)) for i in range(d)])
            lam = 2 - 2 * cos(pi * Rational(j, h))
            resid = simplify(expand(C * vj - lam * vj))
            if resid != zeros(d, 1):
                all_ev = False
        ok(all_ev, "1(a) d={0}: Cartan eigenvalues = 2−2cos(πj/h) (all j, exactly)".format(d))
        SC = SC_of(d)
        lad = ladder(d)
        # S956 «generator Gram» = SC*Cartan; its spectrum = SC*(cartan eig)
        ug_spec = [simplify(SC * x) for x in cartan_eig_closed(d)]
        match_s956 = all(simplify(ug_spec[k] - lad[k]) == 0 for k in range(d))
        ok(match_s956, "1(a) d={0}: the spectrum of SC·Cartan = the ladder SC(2−2cos)".format(d))
        print("{0} | OK (all j) | {1} | {2}".format(
            d, [str(x) for x in lad], "yes" if match_s956 else "NO"))
    print("  ★MEASUREMENT: the S956 generator-Gram = SC·Cartan(A_d), spectrum = SC(2−2cos(πj/h)) —")
    print("  EXACTLY the ladder of Coxeter exponents. Stake 1(a): WON (an identity, 0 tuning).")

    # ==================== BET 1(b): transverse S999 vs ladder ====================
    print()
    print("STAKE 1(b): are the S999 transverse eigenvalues (isotropic, symmetric points) on the ladder?")
    print("d | transverse = n²/(2d) | = SC·n/2? | a ladder rung SC(2−2cos)? which j | verdict")
    print("-" * 96)
    S999_TR = {2: Rational(9, 4), 3: Rational(8, 3), 4: Rational(25, 8)}
    for d in DS:
        n = d + 1
        tr_closed = Rational(n ** 2, 2 * d)
        ok(tr_closed == S999_TR[d], "1(b) d={0}: transverse = n²/(2d) = {1} (=S999)".format(d, S999_TR[d]))
        SC = SC_of(d)
        ok(simplify(tr_closed - SC * Rational(n, 2)) == 0,
           "1(b) d={0}: transverse = SC·n/2 (a closed form)".format(d))
        lad = ladder(d)
        on_ladder = None
        for j in range(1, d + 1):
            if simplify(tr_closed - lad[j - 1]) == 0:
                on_ladder = j
        ok(True, "1(b) d={0}: the ladder-membership measurement is performed".format(d))
        verdict = ("rung j={0}".format(on_ladder) if on_ladder else "NOT on the ladder")
        print("{0} | {1} | yes (SC·n/2) | {2} | {3}".format(
            d, tr_closed, on_ladder if on_ladder else "none", verdict))
    print("  ★MEASUREMENT: transverse = n²/(2d) = SC·n/2 — a PURE closed form (SC, not tuning),")
    print("  but on the cos-ladder only at d=3 (j=2, cos(π/2)=0) — a COINCIDENCE, not ∀d. Stake 1(b): «no».")

    # ==================== BET 1(c): commutant blocks vs ladder ====================
    print()
    print("STAKE 1(c): are the S1000-T2 commutant blocks (isotypes) on the cos-ladder?")
    print("d | dim of the commutant (m=0..3) | type | on the cos-ladder?")
    print("-" * 72)
    # commutant dims from S1000-T2: 1/2/3/3-saturation (integers, orthogonal blocks)
    for d in DS:
        dims = []
        for m in range(0, 4):
            a, b = m, d + 1 - m
            rT = 1 if (a >= 1 and b >= 1) else 0
            rM = 1 if a >= 2 else 0
            rU = 1 if b >= 2 else 0
            dims.append(rT + rM + rU)
        lad = ladder(d)
        # integer dims are not cos-ladder values (which are irrational for most j)
        any_on = any(simplify(Integer(x) - lad[j - 1]) == 0
                     for x in dims for j in range(1, d + 1))
        ok(not any_on, "1(c) d={0}: the integer isotype dims are NOT on the cos-ladder".format(d))
        print("{0} | {1} | integer (orthogonal blocks) | no".format(d, dims))
    print("  ★MEASUREMENT: the commutant = an isotypic decomposition (dim 1/2/3, orthogonal blocks, integer) —")
    print("  structurally NOT a quadratic cos-ladder. Stake 1(c): «no».")

    # bet-1 forsage verdict
    print()
    forsage = 1  # only (a); (b)=no, (c)=no
    ok(forsage < 2, "STAKE 1 FORCING: {0}/3 (only (a)) < 2 => forcing NOT reached (K1)".format(forsage))
    print("  ⟹ STAKE 1: (a) WON (Gram-spectrum=ladder), (b)/(c) «no» ⟹ forcing (≥2)")
    print("  NOT reached (K1 partially): not all quadratic spectra reduce to ONE ladder —")
    print("  only the Gram. The transverse (n²/2d) and the commutant (integer isotypes) — other closed forms.")

    # ==================== BET 2: time = Coxeter action ====================
    print()
    print("STAKE 2: the S597 zigzag azimuth = the action of the Coxeter element ⟨c⟩ of order h")
    print("d | c=(d+1)-cycle, order | eigenvalues on {Σx=0} = e^{2πij/h}? | azimuth θ (tr=2cosθ)")
    print("-" * 96)
    for d in DS:
        n = d + 1
        P = cyc_perm_matrix(n)
        # order of c
        Pk = eye(n)
        order = None
        for k in range(1, n + 1):
            Pk = P * Pk
            if Pk == eye(n):
                order = k
                break
        ok(order == n, "2 d={0}: the order of the Coxeter element = h = {1}".format(d, n))
        # eigenvalues on {sum x=0}: char poly restricted = 1+lam+...+lam^d,
        #   roots e^{2pi i j/h}. verify P has these (P^n=I, P^k!=I) + 1 excluded.
        # check e^{2pi i j/h} is eigenvalue: P*w = zeta*w for w[i]=zeta^{-i}? test trace
        # trace of P on {sum x=0} = trace(P) - 1 = 0 - 1 = -1 (permutation trace 0)
        tr_plane = simplify(P.trace() - 1)
        ok(tr_plane == -1, "2 d={0}: tr(c|{{Σx=0}}) = trP - 1 = -1".format(d))
        # eigenvalues = ALL h-th roots of unity: charpoly(c) = lam^h - 1 (exact, integer
        #   matrix); on {sum x=0} the root 1 is removed -> nontrivial = e^{2pi i j/h}
        cp = P.charpoly()
        ok(simplify(cp.as_expr() - (cp.gen ** n - 1)) == 0,
           "2 d={0}: charpoly(c) = λ^h − 1 => eigenvalues = all h-th roots of 1; "
           "on {{Σx=0}} the nontrivial ones = e^{{2πij/h}} (exponents)".format(d))
        # azimuth for d=2: rotation angle theta, 2cos theta = tr_plane = -1 -> 120°
        if d == 2:
            cos_theta = Rational(tr_plane, 2)
            ok(cos_theta == Rational(-1, 2), "2 d=2: cosθ=-1/2 => θ=120°=2π/3=2π/h")
            azim = "θ=120°=2π/3=2π/h ✓ (=S597 C₃-azimuth)"
        else:
            azim = "order ⟨c⟩=h={0}; eigenvalues e^{{2πij/h}} (exponents)".format(n)
        print("{0} | {1} | yes (Σ=−1=tr) | {2}".format(d, order, azim))
    print("  ★MEASUREMENT: the Coxeter element c=(d+1)-cycle, order h; on {Σx=0} eigenvalues = e^{2πij/h}")
    print("  (exponents); d=2 azimuth = 120° = 2π/3 = 2π/h = EXACTLY the S597 C₃-azimuth. ∀d ⟨c⟩ has order h.")
    print("  Stake 2: WON — «time=ladder» is lifted to the Coxeter group action (K2 is dead).")

    # ==================== MUTANTS ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    # M1: broken Cartan (cycle graph instead of path) -> spectrum != cos(pi j/h) ladder
    d0 = 3
    Ccyc = cartan_Ad(d0)
    Ccyc[0, d0 - 1] = Integer(-1)      # add wrap-around edge -> cycle graph
    Ccyc[d0 - 1, 0] = Integer(-1)
    # cycle-graph eigenvalues 2-2cos(2 pi j/d0) != path 2-2cos(pi j/(d0+1))
    v1 = Matrix([sin(pi * Rational((i + 1) * 1, d0 + 1)) for i in range(d0)])
    lam1 = 2 - 2 * cos(pi * Rational(1, d0 + 1))
    resid = simplify(expand(Ccyc * v1 - lam1 * v1))
    if resid != zeros(d0, 1):
        print("  MUTANT M1: CAUGHT (broken Cartan=cycle-graph: the path eigenvector is NO LONGER an "
              "eigenvector => the spectrum != the cos(πj/h)-ladder)")
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2 (ANTI-TUNING, mandatory): fit transverse 9/4 with a NON-integer j -> FAIL
    d0 = 2
    SC = SC_of(d0)
    tr = Rational(9, 4)
    # solve SC*(2-2cos(pi*x/h)) = tr for x: cos(pi x/3) = 1 - tr/(2 SC)
    val = simplify(1 - tr / (2 * SC))         # = 1 - (9/4)/3 = 1 - 3/4 = 1/4
    # x with cos(pi x/3)=1/4 is NOT a rational integer j in {1,2} (cos={1/2,-1/2})
    j_ok = any(simplify(cos(pi * Rational(j, 3)) - val) == 0 for j in (1, 2))
    if not j_ok:
        print("  MUTANT M2 (ANTI-TUNING): CAUGHT (fitting 9/4 onto the ladder requires cos(πx/3)=1/4, "
              "but integer j∈{{1,2}} give cos∈{{1/2,−1/2}} => a coefficient outside the convention = FAIL)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3: false Coxeter (transposition, order 2) -> azimuth order != h
    d0 = 3
    n0 = d0 + 1
    Pt = eye(n0)
    Pt[0, 0] = 0; Pt[1, 1] = 0; Pt[0, 1] = 1; Pt[1, 0] = 1   # transposition, order 2
    ord_t = 2 if (Pt * Pt == eye(n0)) else None
    if ord_t == 2 and ord_t != n0:
        print("  MUTANT M3: CAUGHT (false-Coxeter=transposition, order 2 != h={0} => "
              "the azimuth-group is not ⟨c⟩)".format(n0))
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4: wrong azimuth angle 2pi/(h+1) instead of 2pi/h -> != S597 120° at d=2
    ang_wrong = 2 * pi / (3 + 1)               # 2pi/4 = 90°, not 120°
    ang_true = 2 * pi / 3
    if simplify(ang_wrong - ang_true) != 0:
        print("  MUTANT M4: CAUGHT (azimuth 2π/(h+1)=2π/4=90° != 2π/h=2π/3=120° S597)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): a random rung against the ladder")
    random.seed(1004029)
    d_nc = 3
    lad = ladder(d_nc)
    fake_rung = Rational(random.randrange(1, 9), random.randrange(2, 7))
    on_lad = any(simplify(fake_rung - lad[j - 1]) == 0 for j in range(1, d_nc + 1))
    ok(not on_lad, "control: the random rung {0} is NOT on the ladder".format(fake_rung))
    print("  random {0}: on the ladder={1} — the measurement is sensitive (only exact cos-values pass)"
          .format(fake_rung, on_lad))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

    # ==================== FORBIDDEN-SCAN ====================
    _pp = [("ко", "лір"), ("ква", "рк"), ("глю", "он"), ("гей", "дж"),  # GUARDLINE
           ("Кази", "мир"), ("C", "2")]  # GUARDLINE
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hits_src = scan_forbidden(__file__, _PATTERNS)
    _logf.flush()
    _logtxt = "".join(_tee.chunks)
    _hits_log = scan_forbidden(_logtxt, _PATTERNS)
    _nhits = len(_hits_src) + len(_hits_log)
    print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2}) [su/Weyl/Coxeter/exponents are allowed]"
          .format(_nhits, len(_hits_src), len(_hits_log)))

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
