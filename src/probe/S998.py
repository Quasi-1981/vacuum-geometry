# -*- coding: utf-8 -*-
# DIM: na (W40: marked cell axes — spectral trace of marking; exact sympy; 0 handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — blind probe, reuse of S956
# ----------------------------------------------------------------------------
# CELL d (as S956): d+1 unit axes, pairwise product -1/d; realization u_i =
#   e_i - c in {sum x = 0} of R^{d+1}; SC = (d+1)/d rational scale; closure =
#   A_d root lattice; bipartite; bond vectors delta_i = the d+1 axes.
# MARKING: first m of the d+1 axes get weight t (rational), the rest weight 1.
#   Bloch: H(k) = [[0, f], [f*, 0]], f_t(k) = sum_i w_i exp(2 pi I <k, delta_i>).
#   Base t=1 = S956 TASK 4 verbatim (bit-fence T1 against its run.log).
# BZ REDUCTION (exact): psi_i = <k, delta_i - delta_0> mod 1 (i = 1..d) are
#   coordinates on the BZ torus (delta_i - delta_0 = e_i - e_0 = basis of A_d;
#   the pairing L x L* -> Z gives the iso).  f = z_0 * g with |z_0| = 1, so
#   f = 0  <=>  g = 0, where with w_i = exp(2 pi I psi_i):
#     m=1:  g_t(psi) = t + sum_{i=1..d} w_i;
#     m=2:  g_t(psi) = t + t*w_1 + sum_{i=2..d} w_i.
# EXACT CLASSIFICATION of {g_t = 0} (verdicts on sympy exact objects only):
#   (E) empty:    m=1, t > d  =>  |g| >= t - d > 0 (triangle bound, rational).
#   (P) point:    m=1, t = d  =>  Re g = t + sum cos(2 pi psi_i) >= t - d = 0,
#                 equality iff every w_i = -1: unique psi = (1/2, ..., 1/2).
#   (M) variety:  explicit exact sample (g(sample) = 0 by rational cancel) +
#                 rank of the real 2 x d Jacobian J (columns c_i*(-Im w_i, Re w_i))
#                 at the sample; local dim = d - rank(J) at a regular sample.
#   Singular (collinear) configurations: all w_i in {+1, -1}  =>  m=1: t = q - p
#   with q + p = d;  m=2: t*(1 + w_1) + (p - q) = 0 — enumerated + shown exactly.
# POLYGON CRITERION: a zero exists on the torus  <=>  the moduli multiset
#   {t x m, 1 x (d+1-m)} closes a polygon  <=>  2*max <= sum.  The product-of-
#   phases tie is killed by a global phase: z_i -> z_i * e^{I a} rescales the
#   product by e^{I (d+1) a}, so any polygon solution lands on the torus.
#   Critical t* = saturation of the inequality.
# Discipline: 0 handles; mutants M1..M5 on every classification branch; seeded
# negative control; FORBIDDEN-SCAN via tools.fence_scan (source + log);
# bit-fence T1 vs S956 run.log read verbatim; STOP after the tables.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import random
from sympy import (Matrix, Integer, Rational, zeros, ones, exp, I, pi,
                   simplify, Add, sqrt, conjugate, expand, symbols, solve)

_HERE = os.path.dirname(os.path.abspath(__file__))
_LOGPATH = os.path.join(_HERE, "S998_run.log")
_logf = open(_LOGPATH, "w", encoding="utf-8")


class Tee:
    def __init__(self, real, fh):
        self.real = real
        self.fh = fh
        self.chunks = []

    def write(self, s):
        self.real.write(s)
        self.fh.write(s)
        self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


_tee = Tee(sys.stdout, _logf)
sys.stdout = _tee

sys.path.insert(0, os.path.join(_HERE, ".."))
from tools.fence_scan import scan_forbidden

ASSERT_PASS = [0]
FAILS = [0]


def ok(cond, msg):
    if cond:
        ASSERT_PASS[0] += 1
    else:
        FAILS[0] += 1
        print("ASSERT-FAIL: " + msg)


# ==================== rational primitives (verbatim from S956) ====================

def cell_vectors(d):
    """u_i = e_i - centroid in R^{d+1}, exact rational; and unit-scale SC."""
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


def cartan_Ad(d):
    M = zeros(d, d)
    for i in range(d):
        M[i, i] = Integer(2)
        if i + 1 < d:
            M[i, i + 1] = Integer(-1)
            M[i + 1, i] = Integer(-1)
    return M


def root_basis(d):
    n = d + 1
    out = []
    for p in range(d):
        v = zeros(n, 1)
        v[p, 0] = Integer(1)
        v[p + 1, 0] = Integer(-1)
        out.append(v)
    return out


def coset_search(d, us, SC, box=2):
    basis = root_basis(d)
    u0 = us[0]
    best = None
    minima = []
    ranges = [range(-box, box + 1)] * d
    import itertools
    for combo in itertools.product(*ranges):
        w = u0.copy()
        for coef, b in zip(combo, basis):
            if coef:
                w = w + Integer(coef) * b
        nrm = udot(w, w, SC)
        if best is None or nrm < best:
            best = nrm
            minima = [w]
        elif nrm == best:
            minima.append(w)
    return best, minima


def bloch_symmetric_point(d, j):
    """f(k_j) at the symmetric coset point, exact (geometric-series identity).
    Verbatim S956 TASK 4 (t=1 base)."""
    n = d + 1
    w = zeros(n, 1)
    for i in range(n):
        w[i, 0] = Rational(2 * i - d, 2 * n)
    kj = Integer(j) * w
    us, _ = cell_vectors(d)
    ok(edot(w, ones(n, 1)) == 0, "w in {sum=0} hyperplane")
    phases_ok = True
    for i in range(n):
        ph = edot(kj, us[i])
        want = Rational(j * (2 * i - d), 2 * n)
        if ph != want:
            phases_ok = False
    ok(phases_ok, "phases <k_j,delta_i> concrete-consistent (d={0},j={1})".format(d, j))
    omega = exp(2 * pi * I * Rational(j % n, n))
    num = simplify(omega ** n - 1)
    ok(num == 0, "omega^n = 1 exact (d={0},j={1})".format(d, j))
    if (j % n) == 0:
        f = Integer(n)
    else:
        f = Integer(0)
    direct = Add(*[exp(2 * pi * I * Rational((j * i) % n, n)) for i in range(n)])
    dval = complex(direct.evalf(30))
    ok(abs(dval - complex(f)) < 1e-20,
       "direct sum ~= geom-series numerically (d={0},j={1})".format(d, j))
    f2 = f * conjugate(f)
    return f, f2


# ==================== marked machinery (new, exact) ====================

def sym_f2_weighted(d, j, marked, t):
    """|f|^2 at symmetric point k_j with weight t on axes in `marked`, 1 elsewhere.
    Exact: f_j = geom_j + (t-1)*sum_{a in marked} omega^{a j}, geom_j = n*[j==0]
    (geometric-series identity established per (d,j) in bloch_symmetric_point)."""
    n = d + 1
    base = Integer(n) if (j % n) == 0 else Integer(0)
    if marked:
        dev = Add(*[(t - 1) * exp(2 * pi * I * Rational((a * j) % n, n))
                    for a in marked])
    else:
        dev = Integer(0)
    f = base + dev
    f2 = simplify(expand(f * conjugate(f)))
    return f2


def sym_node_count(d, marked, t):
    cnt = 0
    for j in range(d + 1):
        if sym_f2_weighted(d, j, marked, t) == 0:
            cnt += 1
    return cnt


def g_val(sample, c, t):
    """g = t + sum_i c_i * w_i at sample [(x_i, y_i)]; returns (Re, Im) exact."""
    gx = t + Add(*[c[i] * sample[i][0] for i in range(len(sample))])
    gy = Add(*[c[i] * sample[i][1] for i in range(len(sample))])
    return simplify(gx), simplify(gy)


def jac_rank(sample, c):
    """Rank of the real 2 x d Jacobian of (Re g, Im g) wrt psi at the sample:
    column i = c_i * (-y_i, x_i)   (2 pi factor dropped: rank-invariant)."""
    M = zeros(2, len(sample))
    for i, (x, y) in enumerate(sample):
        M[0, i] = -c[i] * y
        M[1, i] = c[i] * x
    return M.rank()


def poly_ok(d, t, m, strict=False):
    """Polygon criterion on moduli {t x m, 1 x (d+1-m)}: 2*max <= sum."""
    sides = [t] * m + [Integer(1)] * (d + 1 - m)
    mx = max(sides)
    total = Add(*sides)
    return (2 * mx < total) if strict else (2 * mx <= total)


def build_sample_m1(d, t):
    """Exact zero of g_t = t + sum_{i=1..d} w_i for 0 <= t < d:
    conj pair (x +- I y) + extras (+-1) with x = -(t+s)/2, s = d-2-2j, |x|<1."""
    for j in range(d - 1):
        s = Integer(d - 2 - 2 * j)
        x = -(t + s) / 2
        if abs(x) < 1:
            y = sqrt(1 - x ** 2)
            extras = ([(Integer(1), Integer(0))] * (d - 2 - j)
                      + [(Integer(-1), Integer(0))] * j)
            return [(x, y), (x, -y)] + extras
    return None


def solve_points_d2_m1(t):
    """d=2, m=1, t>0: full exact solve of t + w_1 + w_2 = 0 over the torus."""
    x1, y1, x2, y2 = symbols('x1 y1 x2 y2', real=True)
    sols = solve([x1 + x2 + t, y1 + y2,
                  x1 ** 2 + y1 ** 2 - 1, x2 ** 2 + y2 ** 2 - 1],
                 [x1, y1, x2, y2], dict=True)
    pts = set()
    for s in sols:
        vals = (s[x1], s[y1], s[x2], s[y2])
        if all(v.is_real for v in vals):
            pts.add(vals)
    return sorted(pts, key=str)


def classify_m1(d, t):
    """Exact verdict on {g_t = 0} (m=1) in the d-torus BZ."""
    c = [Integer(1)] * d
    if t > d:                                   # (E) empty
        ok(t - d > 0, "E-cert: t-d>0 (d={0}, t={1})".format(d, t))
        return dict(typ="empty", dim="-", npts=0, rank="-",
                    note="|g| >= t-d = {0} > 0".format(t - d))
    if t == d:                                  # (P) merge point
        sample = [(Integer(-1), Integer(0))] * d
        gx, gy = g_val(sample, c, t)
        ok(gx == 0 and gy == 0, "P: g(all psi=1/2) = 0 exact (d={0}, t={1})".format(d, t))
        r = jac_rank(sample, c)
        ok(r == 1, "P: rank(J)=1 at merge point (d={0})".format(d))
        return dict(typ="point", dim=0, npts=1, rank=r,
                    note="unique psi=(1/2,..,1/2): Re g >= t-d=0, equality <=> all w=-1")
    if t == 0 and d == 2:                       # (M) degenerate: circle
        s = symbols('s', real=True)
        gsym = exp(2 * pi * I * s) + exp(2 * pi * I * s) * exp(2 * pi * I * Rational(1, 2))
        ok(simplify(gsym) == 0, "d=2 t=0: g == 0 identically on psi_2 = psi_1 + 1/2")
        sample = [(Integer(0), Integer(1)), (Integer(0), Integer(-1))]
        r = jac_rank(sample, c)
        ok(r == 1, "d=2 t=0: rank(J)=1 on the circle")
        return dict(typ="variety", dim=1, npts="-", rank=r,
                    note="circle psi_2 - psi_1 = 1/2 (saturation 1 = t + 1)")
    # (M) generic 0 <= t < d
    sample = build_sample_m1(d, t)
    ok(sample is not None, "M: sample exists (d={0}, t={1})".format(d, t))
    gx, gy = g_val(sample, c, t)
    ok(gx == 0 and gy == 0, "M: g(sample) = 0 exact (d={0}, t={1})".format(d, t))
    r = jac_rank(sample, c)
    ok(r == 2, "M: rank(J)=2 at regular sample (d={0}, t={1})".format(d, t))
    dim = d - r
    if d == 2:
        pts = solve_points_d2_m1(t)
        npts = len(pts)
        ok(npts == 2, "d=2: exactly 2 solved points (t={0}) got {1}".format(t, npts))
        note = "; ".join("(x={0},y={1})".format(p[0], p[1]) for p in pts)
        return dict(typ="points", dim=0, npts=npts, rank=r, note="solve: " + note)
    return dict(typ="variety", dim=dim, npts="-", rank=r,
                note="dim = d - rank(J) = {0}".format(dim))


def classify_m2(d, t):
    """Exact verdict on {g_t = 0} (m=2) in the d-torus BZ, d in {3,4}."""
    c = [t] + [Integer(1)] * (d - 1)
    ok(poly_ok(d, t, 2), "m=2 polygon criterion holds (d={0}, t={1})".format(d, t))
    if d == 3:
        sample = [(Integer(-1), Integer(0)), (Integer(0), Integer(1)),
                  (Integer(0), Integer(-1))]
    else:
        s3 = sqrt(3)
        sample = [(Integer(-1), Integer(0)), (Integer(1), Integer(0)),
                  (Rational(-1, 2), s3 / 2), (Rational(-1, 2), -s3 / 2)]
    gx, gy = g_val(sample, c, t)
    ok(gx == 0 and gy == 0, "m=2: g(sample) = 0 exact (d={0}, t={1})".format(d, t))
    r = jac_rank(sample, c)
    dim = d - r
    note = "dim = d - rank(J)"
    if t == 0:
        note = "marked axes excluded: d-1 bonds; psi_1 free => factor T^1"
    return dict(typ="variety", dim=dim, npts="-", rank=r, note=note)


# ============================================================================
#                                   RUN
# ============================================================================

print("=" * 74)
print("W40: marked cell axes — spectral trace of marking (blind probe)")
print("=" * 74)

# ==================== T1: bit-fence against S956 run.log ====================
print()
print("T1 (bit-fence, base t=1, d=2..4): cross-check against S956 run.log BIT-FOR-BIT")
S956LOG = os.path.join(_HERE, "S956_run.log")
_s956txt = open(S956LOG, encoding="utf-8").read().splitlines()
s956rows = {}
for _ln in _s956txt:
    _parts = [p.strip() for p in _ln.split("|")]
    if len(_parts) == 11 and _parts[1] == "SC*Cartan(A_d)":
        s956rows[int(_parts[0])] = _parts
ok(all(dd in s956rows for dd in (2, 3, 4)), "S956 log rows d=2,3,4 parsed")

t1_nodes = {}
for d in (2, 3, 4):
    n = d + 1
    us, SC = cell_vectors(d)
    ok(all(udot(us[i], us[i], SC) == 1 for i in range(n)),
       "cell unit norm =1 (d={0})".format(d))
    ok(all(udot(us[i], us[j], SC) == Rational(-1, d)
           for i in range(n) for j in range(n) if i != j),
       "cell pairwise cos = -1/d (d={0})".format(d))
    basis = root_basis(d)
    G_eucl = Matrix(d, d, lambda i, j: edot(basis[i], basis[j]))
    G_unit = SC * G_eucl
    detE = G_eucl.det()
    detU = G_unit.det()
    ok(G_eucl == cartan_Ad(d), "Euclid Gram = Cartan(A_d) (d={0})".format(d))
    ok(detE == d + 1, "det Gram = d+1 (d={0})".format(d))
    dmin, minima = coset_search(d, us, SC, box=2)
    z = len(minima)
    ok(dmin == 1, "dist^2 = 1 (d={0})".format(d))
    ok(z == d + 1, "z = d+1 (d={0})".format(d))
    node_js = []
    for j in range(n):
        f, f2 = bloch_symmetric_point(d, j)
        if f2 == 0:
            node_js.append(j)
    nnodes = len(node_js)
    ok(nnodes == d, "zeros of f at symmetric points = d (d={0}) got {1}".format(d, nnodes))
    if d == 2:
        ok(nnodes == 2, "d=2: exactly 2 zeros")
    t1_nodes[d] = nnodes
    # bit-fence: string-exact against the S956 log row fields
    row = s956rows[d]
    mine = {2: str(SC), 3: str(detE), 4: str(detU), 5: str(z),
            6: str(udot(us[0], us[1], SC)), 9: str(nnodes), 10: str(dmin)}
    for idx, val in mine.items():
        ok(row[idx] == val,
           "bit-fence d={0} field[{1}]: '{2}' == '{3}'".format(d, idx, row[idx], val))
    print("  d={0}: SC={1} detE={2} detU={3} z={4} bondcos={5} zeros={6} dist2={7}"
          "  <- matches S956 run.log".format(d, SC, detE, detU, z,
                                            udot(us[0], us[1], SC), nnodes, dmin))

# ==================== T2: zero-set of the marked sum, m=1 ====================
print()
print("T2 (m=1): {g_t = 0} on the BZ torus, exact classification")
TVALS = [Integer(0), Rational(1, 2), Integer(1), Rational(3, 2), Integer(2),
         Rational(5, 2), Integer(3), Rational(7, 2), Integer(4)]

EXPECT_T2 = {}  # (d,t) -> (typ, dim, npts)
for tv in TVALS:
    for d in (2, 3, 4):
        if tv > d + 1:
            continue
        if tv > d:
            EXPECT_T2[(d, tv)] = ("empty", "-", 0)
        elif tv == d:
            EXPECT_T2[(d, tv)] = ("point", 0, 1)
        elif d == 2 and tv == 0:
            EXPECT_T2[(d, tv)] = ("variety", 1, "-")
        elif d == 2:
            EXPECT_T2[(d, tv)] = ("points", 0, 2)
        else:
            EXPECT_T2[(d, tv)] = ("variety", d - 2, "-")

hdr2 = "d | t | type | dim | #points | rank(J) | note"
print(hdr2)
print("-" * 100)
t2rows = []
for d in (2, 3, 4):
    for tv in [x for x in TVALS if x <= d + 1]:
        v = classify_m1(d, tv)
        exp_typ, exp_dim, exp_npts = EXPECT_T2[(d, tv)]
        ok(v["typ"] == exp_typ and v["dim"] == exp_dim and v["npts"] == exp_npts,
           "T2 verdict (d={0}, t={1}): {2}/{3}/{4}".format(
               d, tv, v["typ"], v["dim"], v["npts"]))
        t2rows.append((d, tv, v))
        print("{0} | {1} | {2} | {3} | {4} | {5} | {6}".format(
            d, tv, v["typ"], v["dim"], v["npts"], v["rank"], v["note"]))

print()
print("T2 critical t* (symbolic, polygon saturation 2*max = sum):")
for d in (2, 3, 4):
    line = "  d={0}: t* = d = {0} (the zeros merge into psi=(1/2,..,1/2) and vanish)".format(d)
    if d == 2:
        line += " ; lower t* = 2-d = 0 (points -> circle, smeared)"
    sing = [Integer(d - 2 * j) for j in range(d // 2 + 1) if d - 2 * j >= 0]
    line += " ; collinear t (singular points on the zero-set): {0}".format(sing)
    print(line)
# demonstrations of singular points (rank(J)=1 on the zero-set), exact:
_c3 = [Integer(1)] * 3
_sngl3 = [(Integer(1), Integer(0)), (Integer(-1), Integer(0)), (Integer(-1), Integer(0))]
_gx, _gy = g_val(_sngl3, _c3, Integer(1))
ok(_gx == 0 and _gy == 0, "d=3 t=1: collinear configuration (1,-1,-1) is a zero")
ok(jac_rank(_sngl3, _c3) == 1, "d=3 t=1: rank(J)=1 there (singular point)")
_c4 = [Integer(1)] * 4
_sngl4a = [(Integer(1), Integer(0))] * 2 + [(Integer(-1), Integer(0))] * 2
_gx, _gy = g_val(_sngl4a, _c4, Integer(0))
ok(_gx == 0 and _gy == 0, "d=4 t=0: (1,1,-1,-1) is a zero")
ok(jac_rank(_sngl4a, _c4) == 1, "d=4 t=0: rank(J)=1 there")
_sngl4b = [(Integer(-1), Integer(0))] * 3 + [(Integer(1), Integer(0))]
_gx, _gy = g_val(_sngl4b, _c4, Integer(2))
ok(_gx == 0 and _gy == 0, "d=4 t=2: (-1,-1,-1,1) is a zero")
ok(jac_rank(_sngl4b, _c4) == 1, "d=4 t=2: rank(J)=1 there")
print("  demonstrations: d=3 t=1 (1,-1,-1) rank=1 ; d=4 t=0 (1,1,-1,-1) rank=1 ;"
      " d=4 t=2 (-1,-1,-1,1) rank=1 — all exact")

# t=1: symmetric points lie ON the zero-set (geometric sum = 0, T1)
for d in (2, 3, 4):
    ok(t1_nodes[d] == d,
       "t=1: all j=1..d symmetric points are zeros of g (d={0})".format(d))

# ==================== T3: m=2 for d in {3,4} ====================
print()
print("T3 (m=2, weight t on two axes): {g_t = 0}, d in {3,4}, t in {0,1/2,1,2}")
EXPECT_T3 = {(3, Integer(0)): 2, (3, Rational(1, 2)): 1, (3, Integer(1)): 1,
             (3, Integer(2)): 1,
             (4, Integer(0)): 2, (4, Rational(1, 2)): 2, (4, Integer(1)): 2,
             (4, Integer(2)): 2}
print(hdr2)
print("-" * 100)
t3rows = []
for d in (3, 4):
    for tv in [Integer(0), Rational(1, 2), Integer(1), Integer(2)]:
        v = classify_m2(d, tv)
        ok(v["dim"] == EXPECT_T3[(d, tv)],
           "T3 dim (d={0}, t={1}) = {2}".format(d, tv, v["dim"]))
        t3rows.append((d, tv, v))
        print("{0} | {1} | {2} | {3} | {4} | {5} | {6}".format(
            d, tv, v["typ"], v["dim"], v["npts"], v["rank"], v["note"]))
# m=2 collinear singular points (exact):
_cm2_3 = [symbols('tt', positive=True)] * 1  # placeholder unused
for tv in [Rational(1, 2), Integer(1), Integer(2)]:
    _c = [tv, Integer(1), Integer(1)]
    _s = [(Integer(-1), Integer(0)), (Integer(1), Integer(0)), (Integer(-1), Integer(0))]
    _gx, _gy = g_val(_s, _c, tv)
    ok(_gx == 0 and _gy == 0,
       "d=3 m=2 t={0}: (-1,1,-1) is a zero (independent of t)".format(tv))
    ok(jac_rank(_s, _c) == 1, "d=3 m=2 t={0}: rank(J)=1 there".format(tv))
_c = [Rational(1, 2), Integer(1), Integer(1), Integer(1)]
_s = [(Integer(1), Integer(0)), (Integer(1), Integer(0)),
      (Integer(-1), Integer(0)), (Integer(-1), Integer(0))]
_gx, _gy = g_val(_s, _c, Rational(1, 2))
ok(_gx == 0 and _gy == 0, "d=4 m=2 t=1/2: (1,1,-1,-1) is a zero")
ok(jac_rank(_s, _c) == 1, "d=4 m=2 t=1/2: rank(J)=1 (singular)")
print("  singular points m=2: d=3 — (-1,1,-1) for ALL t (rank=1); d=4 — at t=1/2 (rank=1)")

# ==================== T4: table of two prescriptions for the pair (p,q) ====================
print()
print("T4: pairs (p,q) — C1: d=p, m=0  ·  C2: d=p+q-1, m=q, t=1")
hdr4 = ("pair (p,q) | cell | d | m | t | rank | det Gram(unit) | z | "
        "zeros(symm) | dim zero-set")
print(hdr4)
print("-" * len(hdr4))
DIM_FULL = {}  # d -> dim of the full zero-set at t=1 (from T2 classification)
for d in (2, 3, 4):
    v = classify_m1(d, Integer(1))
    DIM_FULL[d] = v["dim"]
ok(DIM_FULL == {2: 0, 3: 1, 4: 2}, "t=1 full zero-set: dim = {2:0, 3:1, 4:2}")

t4rows = []
for (p, q) in [(2, 1), (3, 1), (2, 2), (3, 2)]:
    for label, dd, m in [("C1", p, 0), ("C2", p + q - 1, q)]:
        us, SC = cell_vectors(dd)
        basis = root_basis(dd)
        G_unit = SC * Matrix(dd, dd, lambda i, j: edot(basis[i], basis[j]))
        rk = G_unit.rank()
        detU = G_unit.det()
        dmin, minima = coset_search(dd, us, SC, box=2)
        z = len(minima)
        marked = list(range(m))
        nodes = sym_node_count(dd, marked, Integer(1))
        ok(rk == dd, "T4 {0}({1},{2}): rank = d".format(label, p, q))
        ok(nodes == t1_nodes[dd],
           "T4 {0}({1},{2}): zeros(symm) at t=1 = base d={3}".format(label, p, q, dd))
        row = ((p, q), label, dd, m, 1, rk, detU, z, nodes, DIM_FULL[dd])
        t4rows.append(row)
        print("({0},{1}) | {2} | {3} | {4} | {5} | {6} | {7} | {8} | {9} | {10}".format(
            p, q, label, dd, m, 1, rk, detU, z, nodes, DIM_FULL[dd]))

# ==================== MUTANTS ====================
print()
print("MUTANTS (each classification branch):")
mut_ok = True

# M1: Gram scale making pairwise product -1/(d+1) instead of -1/d -> T1 breaks
d0 = 3
us0, SC0 = cell_vectors(d0)
SC_bad = Integer(1)
pair_bad = udot(us0[0], us0[1], SC_bad)
norm_bad = udot(us0[0], us0[0], SC_bad)
if (pair_bad == Rational(-1, d0 + 1) and pair_bad != Rational(-1, d0)
        and norm_bad != 1
        and SC_bad ** d0 * (d0 + 1) != SC0 ** d0 * (d0 + 1)):
    print("  MUTANT M1: CAUGHT (pair={0} != -1/d; norm={1} != 1; det shifted)"
          .format(pair_bad, norm_bad))
else:
    print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

# M2: t = 1 + 1/100 in the T1 cross-check -> bit-fence breaks (zeros d -> 0)
nodes_m2 = sym_node_count(3, [0], Rational(101, 100))
if nodes_m2 == 0 and nodes_m2 != int(s956rows[3][9]):
    print("  MUTANT M2: CAUGHT (t=101/100: zeros(symm)={0} != {1} from S956 log)"
          .format(nodes_m2, s956rows[3][9]))
else:
    print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

# M3: polygon criterion <= -> < : at t=d the mutant says "empty",
#     but an exact zero exists (all w=-1) -> contradiction
t_sat = Integer(3)
mut_empty = not poly_ok(3, t_sat, 1, strict=True)
gx_wit, gy_wit = g_val([(Integer(-1), Integer(0))] * 3, [Integer(1)] * 3, t_sat)
if mut_empty and gx_wit == 0 and gy_wit == 0 and poly_ok(3, t_sat, 1):
    print("  MUTANT M3: CAUGHT (strict criterion: 'empty' at t=d=3, but "
          "an exact zero psi=(1/2,1/2,1/2) exists; t* would shift)")
else:
    print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

# M4: off-by-one in the count of symmetric zeros (loop j=0..d-1 misses j=d)
d0 = 3
cnt_off = 0
for j in range(0, d0):  # off-by-one: j=d is skipped
    if sym_f2_weighted(d0, j, [], Integer(1)) == 0:
        cnt_off += 1
if cnt_off == d0 - 1 and cnt_off != t1_nodes[d0]:
    print("  MUTANT M4: CAUGHT (off-by-one count={0} != {1})".format(cnt_off, t1_nodes[d0]))
else:
    print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

# M5: substituting m (2 -> 1) in T3: the m=2 sample does NOT satisfy the m=1 equation
tv = Rational(1, 2)
s_m2 = [(Integer(-1), Integer(0)), (Integer(0), Integer(1)), (Integer(0), Integer(-1))]
gx2, gy2 = g_val(s_m2, [tv, Integer(1), Integer(1)], tv)      # m=2: zero
gx1, gy1 = g_val(s_m2, [Integer(1)] * 3, tv)                  # m=1: NOT zero
if gx2 == 0 and gy2 == 0 and (gx1 != 0 or gy1 != 0):
    print("  MUTANT M5: CAUGHT (m=2 sample: g_m2=0, but g_m1={0}+I*{1} != 0)"
          .format(gx1, gy1))
else:
    print("  MUTANT M5: NOT CAUGHT"); mut_ok = False

# ==================== NEGATIVE CONTROL (seeded) ====================
print()
print("NEGATIVE CONTROL (seeded): weight 7/5 on a NON-symmetric axis")
random.seed(940031)
d_nc = 3
axis_nc = random.randrange(1, d_nc + 1)   # axis != 0 (non-central)
nodes_nc = sym_node_count(d_nc, [axis_nc], Rational(7, 5))
f2_j1 = sym_f2_weighted(d_nc, 1, [axis_nc], Rational(7, 5))
ok(f2_j1 == Rational(4, 25), "control: |f|^2(j=1) = 4/25 exact (axis {0})".format(axis_nc))
ok(nodes_nc != t1_nodes[d_nc],
   "control: the zero-set is DIFFERENT from the base ({0} != {1})".format(nodes_nc, t1_nodes[d_nc]))
print("  d={0}, axis={1}, w=(1,..,7/5,..): zeros(symm)={2} vs base {3}; "
      "|f|^2(j=1)=4/25 != 0 — the control fired (the measurement is sensitive)"
      .format(d_nc, axis_nc, nodes_nc, t1_nodes[d_nc]))

# ==================== SUMMARY ====================
print()
print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

# ==================== FORBIDDEN-SCAN (tools.fence_scan) ====================
_pp = [("ч", "ас"), ("t", "ime"), ("тем", "пор"), ("кау", "зал"),
       ("сві", "тло"), ("Мінков", "ськ"), ("сигна", "тур")]
_PATTERNS = ["".join(ab) for ab in _pp]
_hits_src = scan_forbidden(__file__, _PATTERNS)
_logf.flush()
_logtxt = "".join(_tee.chunks)
_hits_log = scan_forbidden(_logtxt, _PATTERNS)
_nhits = len(_hits_src) + len(_hits_log)
print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2})".format(
    _nhits, len(_hits_src), len(_hits_log)))

_exit = 0
if _nhits > 0:
    _exit = 1
if FAILS[0] > 0:
    _exit = 1
if not mut_ok:
    _exit = 1
print("EXIT={0}".format(_exit))
print("STOP")
sys.exit(_exit)
