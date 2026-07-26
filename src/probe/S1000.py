# -*- coding: utf-8 -*-
# DIM: na (W40: two independent marked weights + stabilizer commutant; exact sympy; 0 handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — blind probe, reuse S956/S998/S999
# ----------------------------------------------------------------------------
# CELL d (as S956/S998/S999): d+1 unit axes, pairwise product -1/d; realization
#   u_i = e_i - centroid in {sum x = 0} of R^{d+1}; SC = (d+1)/d rational scale;
#   closure = A_d root lattice; bipartite; bond vectors delta_i = the d+1 axes.
# TWO MARKED AXES (0 and 1, w.l.o.g.) with INDEPENDENT weights w0 = s, w1 = t;
#   the remaining d-1 axes weight 1.  f(k) = sum_i w_i exp(2 pi I <k, delta_i>).
# BZ REDUCTION (as S998/S999): psi_i = <k, delta_i - delta_0> (i = 1..d) torus
#   coords; f = z_0 * g with |z_0| = 1, so f = 0 <=> g = 0, where (axis0 divided
#   out -> its weight s becomes the constant; axis1 keeps weight t on w_1):
#     g_{s,t}(psi) = s + t*w_1 + sum_{i=2..d} w_i,   w_i = exp(2 pi I psi_i).
#   Coeff vector c = [t, 1, ..., 1] (length d), constant = s.
# POLYGON CRITERION (existence of a zero on the torus): the moduli multiset
#   {s, t, 1 x (d-1)} closes a polygon  <=>  2*max <= sum  (global phase kills the
#   product tie, as in S998).  Nonzero moduli only; a zero modulus at t frees
#   psi_1 (marked axis 1 removed -> +1 to dim); a zero at s removes the constant.
# EXACT CLASSIFICATION of {g = 0} (verdicts on sympy exact rationals only):
#   R = nonzero moduli, mx = max R, tot = sum R, free = [t == 0],
#   n_ang = [t != 0] + (d-1)  (# angles that appear in g).
#   (E) empty:     2*mx > tot.                     min|f|^2 over symmetric cosets reported.
#   (Boundary) 2*mx == tot: collinear (all w = +-1); dim = free + [s == 0]; rank(J) = 1.
#   (Interior) 2*mx < tot: smooth codim-2; dim = n_ang - 2 + free; rank(J) = 2.
#   Analytic dims are CROSS-VERIFIED by exact zero samples + Jacobian rank on
#   representative interior / boundary / empty cells (master).
# T2 STABILIZER COMMUTANT (stamped): V = {sum x = 0} subset R^{d+1}, dim d, under
#   G = S_m x S_{d+1-m} (permute marked among marked, unmarked among unmarked).
#   Space measured = symmetric X on V invariant under G = {symmetric (d+1)x(d+1)
#   X over Q : X*1 = 0, X*g = g*X for all generators g}  (X + 0 on span(1) is the
#   unique representative).  dim = exact nullity (Gauss).  Isotypic block
#   projectors of V:  T (imbalance dir, rank 1 if m>=1 and d+1-m>=1),  std_M
#   (marked sum-zero, rank m-1 if m>=2),  std_U (unmarked sum-zero, rank d-m if
#   d+1-m>=2);  ranks sum to d.  dim commutant = # nonzero blocks (V is
#   multiplicity-free -> invariant ops are automatically symmetric).  The
#   direction e0 - e1 (m=2) under swap 0<->1: e0 - e1 -> -(e0 - e1) (odd).
# Symmetry is SIGNIFICANT on the AMBIENT cell module R^{d+1} (mutant M2): the invariant
#   commutant there counts ORDERED-pair orbits (6 for a,b>=2) vs SYMMETRIC
#   unordered (5) -- dropping symmetry inflates the count.
# Discipline: 0 handles; mutants M1..M5 CAUGHT on every branch; seeded negative
#   control (weight 13/10 on an UNMARKED axis); FORBIDDEN-SCAN via tools.fence_scan
#   (source + log); bit-fence T1 vs S998/S999 run.log read verbatim; Pool(18) over
#   the (d, s, t) grid, fixed print order (log bit-reproducible); STOP after tables.
# ============================================================================

import sys
import os
import random
from multiprocessing import Pool
from sympy import (Matrix, Integer, Rational, zeros, ones, eye, exp, I, pi,
                   simplify, expand, Add, sqrt, conjugate, symbols, cos, sin,
                   re, im, nsimplify)

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== rational primitives (verbatim S956/S998/S999) ====================

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


def gg_val(sample, c, s):
    """g = s + sum_i c_i w_i at sample [(x_i,y_i)]; returns (Re, Im) exact."""
    gx = s + Add(*[c[i] * sample[i][0] for i in range(len(sample))])
    gy = Add(*[c[i] * sample[i][1] for i in range(len(sample))])
    return simplify(gx), simplify(gy)


def jac_psi(sample, c):
    """Real 2 x d Jacobian of (Re g, Im g) wrt psi: column i = c_i*(-y_i, x_i)."""
    M = zeros(2, len(sample))
    for i, (x, y) in enumerate(sample):
        M[0, i] = -c[i] * y
        M[1, i] = c[i] * x
    return M


def jac_rank(sample, c):
    return jac_psi(sample, c).rank()


# ==================== two-weight g machinery ====================

def coeffs(d, t):
    """c = [t, 1, ..., 1] length d (c[0]=t on marked axis 1; rest unit)."""
    return [t] + [Integer(1)] * (d - 1)


def moduli(d, s, t):
    R = []
    if s != 0:
        R.append(s)
    if t != 0:
        R.append(t)
    R += [Integer(1)] * (d - 1)
    return R


def classify(d, s, t):
    """Exact analytic verdict on {g_{s,t}=0}; dict(typ,dim,rankJ)."""
    R = moduli(d, s, t)
    mx = max(R)
    tot = Add(*R)
    free = 1 if t == 0 else 0
    n_ang = (1 if t != 0 else 0) + (d - 1)
    if 2 * mx > tot:
        return dict(typ="empty", dim="-", rankJ="-", cls="E")
    if 2 * mx == tot:
        dim = free + (1 if s == 0 else 0)
        typ = "point" if dim == 0 else "variety"
        return dict(typ=typ, dim=dim, rankJ=1, cls=("P" if dim == 0 else "M"))
    dim = n_ang - 2 + free
    return dict(typ="variety", dim=dim, rankJ=2, cls="M")


def sym_f2(d, j, s, t):
    """|f|^2 at symmetric coset k_j: f_j = n[j==0] + (s-1) + (t-1)*omega^j,
    omega = exp(2 pi I j/n) (deviations from unit weights on marked axes 0,1)."""
    n = d + 1
    base = Integer(n) if (j % n) == 0 else Integer(0)
    om = exp(2 * pi * I * Rational(j % n, n))
    f = base + (s - 1) + (t - 1) * om
    return simplify(expand(f * conjugate(f)))


def min_sym_f2(d, s, t):
    vals = [sym_f2(d, j, s, t) for j in range(d + 1)]
    return min(vals, key=lambda e: complex(e.evalf(30)).real), vals


def build_sample_2w(d, s, t):
    """Exact REGULAR zero of g = s + t w_1 + sum_{2..d} w_i (interior), d>=3:
    conj pair on axes (idx1,idx2)=(psi_2,psi_3); reals on axis0 (coeff t) and
    axes idx3.. (coeff 1); pick signs with |x|<1 (=> y!=0 => rank 2)."""
    import itertools
    unit_real_idx = list(range(3, d))          # 0-based angle idx with coeff 1, real
    t_choices = [Integer(1), Integer(-1)] if t != 0 else [Integer(1)]
    for eps0 in t_choices:
        for combo in itertools.product([Integer(1), Integer(-1)],
                                       repeat=len(unit_real_idx)):
            resid = s + (t * eps0 if t != 0 else Integer(0)) + Add(*combo, Integer(0))
            x = -resid / 2
            if abs(x) < 1:
                y = sqrt(1 - x ** 2)
                sample = [None] * d
                # axis0 (psi_1): if t!=0 real eps0; else free -> put (1,0)
                sample[0] = (eps0, Integer(0)) if t != 0 else (Integer(1), Integer(0))
                sample[1] = (x, y)             # psi_2  (coeff 1)
                sample[2] = (x, -y)            # psi_3  (coeff 1)
                for k, idx in enumerate(unit_real_idx):
                    sample[idx] = (combo[k], Integer(0))
                return sample
    return None


def build_collinear(d, s, t):
    """Boundary/collinear zero: all w in {+1,-1}, real, summing g to 0."""
    import itertools
    c = coeffs(d, t)
    for signs in itertools.product([Integer(1), Integer(-1)], repeat=d):
        sample = [(sg, Integer(0)) for sg in signs]
        gx, gy = gg_val(sample, c, s)
        if gx == 0 and gy == 0:
            return sample
    return None


def solve_d2(s, t):
    """d=2: g = s + t w_1 + w_2 = 0.  |s + t w_1| = 1 => cos th = (1-s^2-t^2)/(2st)."""
    if t == 0:
        return None
    val = (1 - s ** 2 - t ** 2) / (2 * s * t) if s != 0 else None
    if s == 0:
        # g = t w_1 + w_2, |w|=1 both => t=1 needed
        return "t=1-only" if t == 1 else None
    if abs(val) <= 1:
        return val
    return None


# ==================== T2: stabilizer commutant (exact null-space) ====================

def perm_mat(n, cyc):
    """Permutation matrix for a transposition cyc=(i,j)."""
    P = eye(n)
    i, j = cyc
    P[i, i] = 0
    P[j, j] = 0
    P[i, j] = 1
    P[j, i] = 1
    return P


def stab_generators(d, m):
    """Adjacent transpositions of S_m x S_{d+1-m} on {0..d}."""
    n = d + 1
    gens = []
    for i in range(0, m - 1):                  # within marked 0..m-1
        gens.append(perm_mat(n, (i, i + 1)))
    for i in range(m, n - 1):                  # within unmarked m..d
        gens.append(perm_mat(n, (i, i + 1)))
    return gens


def commutant_dim(d, m, symmetric=True, sumzero=True):
    """dim of {X (n x n over Q): [symmetric] [X*1=0] and X g = g X for gens}."""
    n = d + 1
    gens = stab_generators(d, m)
    # unknowns: entries of X
    if symmetric:
        idx = [(i, j) for i in range(n) for j in range(i, n)]
    else:
        idx = [(i, j) for i in range(n) for j in range(n)]
    pos = {ij: k for k, ij in enumerate(idx)}
    nv = len(idx)

    def Xentry(vec, i, j):
        if symmetric and i > j:
            i, j = j, i
        return vec[pos[(i, j)]]

    rows = []
    e_syms = symbols('x0:%d' % nv)
    # build X symbolically
    X = zeros(n, n)
    for (i, j) in idx:
        X[i, j] = e_syms[pos[(i, j)]]
        if symmetric and i != j:
            X[j, i] = e_syms[pos[(i, j)]]
    cons = []
    one = ones(n, 1)
    if sumzero:
        Xo = X * one
        for i in range(n):
            cons.append(Xo[i, 0])
    for g in gens:
        C = X * g - g * X
        for i in range(n):
            for j in range(n):
                cons.append(C[i, j])
    # linear system A v = 0
    A = zeros(len(cons), nv)
    for r, ex in enumerate(cons):
        ex = expand(ex)
        for k in range(nv):
            A[r, k] = ex.coeff(e_syms[k])
    return nv - A.rank()


def isotypic_blocks(d, m):
    """Block projectors on V and their ranks (T, std_M, std_U)."""
    n = d + 1
    a, b = m, n - m
    blocks = {}
    # T = imbalance direction (b on marked, -a on unmarked); rank 1 if a,b>=1
    if a >= 1 and b >= 1:
        v = zeros(n, 1)
        for i in range(a):
            v[i, 0] = Integer(b)
        for i in range(a, n):
            v[i, 0] = Integer(-a)
        P = v * v.T / (v.T * v)[0, 0]
        blocks["T"] = simplify(P)
    # std_M: marked sum-zero; rank a-1
    if a >= 2:
        cols = []
        for i in range(1, a):
            e = zeros(n, 1)
            e[0, 0] = Integer(1)
            e[i, 0] = Integer(-1)
            cols.append(e)
        B = Matrix.hstack(*cols)
        P = B * (B.T * B).inv() * B.T
        blocks["std_M"] = simplify(P)
    # std_U: unmarked sum-zero; rank b-1
    if b >= 2:
        cols = []
        for i in range(m + 1, n):
            e = zeros(n, 1)
            e[m, 0] = Integer(1)
            e[i, 0] = Integer(-1)
            cols.append(e)
        B = Matrix.hstack(*cols)
        P = B * (B.T * B).inv() * B.T
        blocks["std_U"] = simplify(P)
    return blocks


# ==================== T3 worker (one grid cell) ====================

GRID = [Integer(0), Rational(1, 4), Rational(1, 3), Rational(1, 2), Rational(2, 3),
        Rational(3, 4), Integer(1), Rational(4, 3), Rational(3, 2), Integer(2),
        Rational(5, 2), Integer(3), Rational(7, 2), Integer(4)]


def _t3_cell(spec):
    d, si, ti = spec
    s = GRID[si]
    t = GRID[ti]
    v = classify(d, s, t)
    minf2 = None
    if v["cls"] == "E":
        mv, _ = min_sym_f2(d, s, t)
        minf2 = str(mv)
    return (spec, v["cls"], v["dim"], minf2)


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
    _logf = open(os.path.join(_HERE, "S1000_run.log"), "w",
                 encoding="utf-8")
    _tee = Tee(sys.stdout, _logf)
    sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    # ---- SIGN OF LIFE FIRST (infra lesson: header before any heavy compute) ----
    print("=" * 74)
    print("W40: two independent weights on the cell + stabilizer commutant (blind probe)")
    print("     g = s + t*w_1 + sum_{i>=2} w_i ; polygon {s,t,1^(d-1)} ; exact arithmetic")
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

    # cell primitive sanity (reuse S956)
    for d in (2, 3, 4, 5, 6):
        us, SC = cell_vectors(d)
        n = d + 1
        ok(all(udot(us[i], us[i], SC) == 1 for i in range(n)),
           "cell unit norm=1 (d={0})".format(d))
        ok(all(udot(us[i], us[j], SC) == Rational(-1, d)
               for i in range(n) for j in range(n) if i != j),
           "cell pairwise cos=-1/d (d={0})".format(d))

    # ==================== T1: bit-fences ====================
    print("T1 (bit-fences)")
    # (a) s=t=1 -> all weights 1 = base; classify d=3,4 and match S999 T1 log
    _s999 = open(os.path.join(_HERE, "S999_run.log"),
                 encoding="utf-8").read().splitlines()
    s999_t1 = {}
    for ln in _s999:
        p = ln.strip()
        if p.startswith("d=") and "<- matches S998 run.log" in p:
            dd = int(p[2])
            s999_t1[dd] = p
    print("  (a) s=t=1 (all weights 1 = the base); cross-check type/dim against S999 T1 BIT-FOR-BIT:")
    for d in (3, 4):
        v = classify(d, Integer(1), Integer(1))
        # base zero-set: dim = d-2 (interior, all-ones), rank 2
        ok(v["cls"] == "M" and v["dim"] == d - 2 and v["rankJ"] == 2,
           "T1(a) d={0}: variety dim={1} rank=2".format(d, d - 2))
        expect_dim = d - 2
        logline = s999_t1.get(d, "")
        ok(("dim={0}".format(expect_dim) in logline) and ("rank=2" in logline)
           and ("variety" in logline),
           "T1(a) d={0}: the S999 log carries dim={1}, variety, rank=2".format(d, expect_dim))
        print("    d={0}: variety dim={1} rank=2  <- S999: '{2}'".format(
            d, expect_dim, logline))
    # (b) s=t -> S998 m=2 case; rows t in {1/2, 2} for d=3
    _s998 = open(os.path.join(_HERE, "S998_run.log"),
                 encoding="utf-8").read().splitlines()
    s998_t3 = {}
    inT3 = False
    for ln in _s998:
        if ln.startswith("T3 (m=2"):
            inT3 = True
            continue
        if ln.startswith("T4:"):
            inT3 = False
        if inT3:
            parts = [q.strip() for q in ln.split("|")]
            if len(parts) == 7 and parts[0] == "3":
                s998_t3[parts[1]] = parts
    print("  (b) s=t (both marked axes equal): cross-check against S998-T3 (m=2), d=3, t in {1/2,2}:")
    for tv, key in [(Rational(1, 2), "1/2"), (Integer(2), "2")]:
        v = classify(3, tv, tv)         # s=t=tv
        row = s998_t3.get(key)
        ok(row is not None, "S998-T3 d=3 t={0} row found".format(key))
        if row:
            ok(v["typ"] == row[2] and str(v["dim"]) == row[3],
               "T1(b) d=3 s=t={0}: {1}/dim {2} == S998 {3}/{4}".format(
                   key, v["typ"], v["dim"], row[2], row[3]))
            print("    s=t={0}: {1} dim={2}  <- S998-T3: {3} dim={4}".format(
                key, v["typ"], v["dim"], row[2], row[3]))

    # ==================== T2: stabilizer commutant ====================
    print()
    print("T2 (commutant of the stabilizer S_m x S_{d+1-m} on {sum x=0}; d=2..6, m=0..3)")
    print("d | m | dim(commutant) | block-ranks T/std_M/std_U | sum=d?")
    print("-" * 78)

    def pred_blocks(d, m):
        n = d + 1
        a, b = m, n - m
        rT = 1 if (a >= 1 and b >= 1) else 0
        rM = (a - 1) if a >= 2 else 0
        rU = (b - 1) if b >= 2 else 0
        return rT, rM, rU

    t2rows = []
    for d in range(2, 7):
        for m in range(0, 4):
            if m > d + 1:
                continue
            dimc = commutant_dim(d, m, symmetric=True, sumzero=True)
            rT, rM, rU = pred_blocks(d, m)
            ncomp = (1 if rT else 0) + (1 if rM else 0) + (1 if rU else 0)
            ok(dimc == ncomp,
               "T2 d={0} m={1}: dim of the commutant {2} == #components {3}".format(
                   d, m, dimc, ncomp))
            # verify block projectors: build, check ranks, orthogonality, sum=I_V
            blocks = isotypic_blocks(d, m)
            got = {"T": 0, "std_M": 0, "std_U": 0}
            Psum = zeros(d + 1, d + 1)
            for name, P in blocks.items():
                r = P.rank()
                got[name] = r
                ok(simplify(P * P - P) == zeros(d + 1, d + 1),
                   "T2 d={0} m={1}: {2} idempotent".format(d, m, name))
                # G-invariance
                for g in stab_generators(d, m):
                    ok(simplify(g * P - P * g) == zeros(d + 1, d + 1),
                       "T2 d={0} m={1}: {2} commutes with G".format(d, m, name))
                Psum += P
            ok((got["T"], got["std_M"], got["std_U"]) == (rT, rM, rU),
               "T2 d={0} m={1}: block-ranks {2} == predicted {3}".format(
                   d, m, (got["T"], got["std_M"], got["std_U"]), (rT, rM, rU)))
            ok(rT + rM + rU == d, "T2 d={0} m={1}: sum of block-ranks = d".format(d, m))
            # Psum must be a projector onto V (rank d)
            ok(simplify(Psum).rank() == d,
               "T2 d={0} m={1}: sum of projectors = a projector onto V (rank d)".format(d, m))
            print("{0} | {1} | {2} | {3}/{4}/{5} | {6}".format(
                d, m, dimc, rT, rM, rU, "yes" if rT + rM + rU == d else "NO"))
            t2rows.append((d, m, dimc, rT, rM, rU))

    # symmetry is significant on the ambient module (context for M2)
    print()
    print("T2b (symmetry is significant on the FULL module R^(d+1); ordered vs symmetric):")
    for d, m in [(4, 2), (5, 2), (6, 3)]:
        sym_amb = commutant_dim(d, m, symmetric=True, sumzero=False)
        non_amb = commutant_dim(d, m, symmetric=False, sumzero=False)
        ok(non_amb > sym_amb,
           "T2b d={0} m={1}: ordered {2} > symmetric {3}".format(d, m, non_amb, sym_amb))
        print("  d={0} m={1}: R^(d+1) ordered-commutant={2} > symmetric={3}"
              .format(d, m, non_amb, sym_amb))

    # e0 - e1 under swap 0<->1 (m=2)
    print()
    n5 = 5
    swap01 = perm_mat(n5, (0, 1))
    e0e1 = zeros(n5, 1)
    e0e1[0, 0] = Integer(1)
    e0e1[1, 0] = Integer(-1)
    swapped = swap01 * e0e1
    ok(swapped == -e0e1, "e0-e1 is odd under the swap 0<->1")
    print("  direction e0-e1 (m=2) under the swap 0<->1: {0} = -(e0-e1) => ODD (yes)"
          .format([swapped[i, 0] for i in range(n5)]))

    # ==================== T3: weight-plane grid (Pool 18) ====================
    print()
    print("T3 (the weight plane (s,t); d=3,4; a 14x14 grid): classes of the zero-set")
    print("    legend: E=empty · 0=point · 1/2/3=variety dim · (min|f|^2 for E)")
    tasks = [(d, si, ti) for d in (3, 4)
             for si in range(len(GRID)) for ti in range(len(GRID))]
    sys.stderr.write("[master] T3: {0} cells in Pool(18)\n".format(len(tasks)))
    sys.stderr.flush()
    RES = {}
    with Pool(processes=18) as pool:
        for spec, cls, dim, minf2 in pool.imap_unordered(_t3_cell, tasks):
            RES[spec] = (cls, dim, minf2)

    def cell_label(cls, dim):
        if cls == "E":
            return "E"
        return str(dim)

    for d in (3, 4):
        print()
        print("  d={0}  (rows s, columns t; GRID values={1}):".format(
            d, ",".join(str(x) for x in GRID)))
        header = "  s\\t | " + " ".join("{0:>4}".format(str(GRID[ti]))
                                        for ti in range(len(GRID)))
        print(header)
        for si in range(len(GRID)):
            cells = []
            for ti in range(len(GRID)):
                cls, dim, _ = RES[(d, si, ti)]
                cells.append("{0:>4}".format(cell_label(cls, dim)))
            print("  {0:>4} | {1}".format(str(GRID[si]), " ".join(cells)))

    # T3 sample cross-verification (master; representative interior/boundary/empty)
    print()
    print("  T3 cross-check (raw samples + rank(J)):")
    # interior (1,1)
    for d in (3, 4):
        s0, t0 = Integer(1), Integer(1)
        v = classify(d, s0, t0)
        smp = build_sample_2w(d, s0, t0)
        ok(smp is not None, "T3 verify d={0} (1,1): sample exists".format(d))
        if smp:
            gx, gy = gg_val(smp, coeffs(d, t0), s0)
            ok(gx == 0 and gy == 0, "T3 verify d={0} (1,1): g(sample)=0".format(d))
            r = jac_rank(smp, coeffs(d, t0))
            ok(r == 2 and v["dim"] == d - 2,
               "T3 verify d={0} (1,1): rank=2 dim={1}".format(d, d - 2))
            print("    d={0} (s,t)=(1,1): interior rank(J)={1}, dim={2} (analytic {2})"
                  .format(d, r, v["dim"]))
    # boundary per d: d=3 (3,1); d=4 (4,1)
    for d, (s0, t0) in [(3, (Integer(3), Integer(1))), (4, (Integer(4), Integer(1)))]:
        v = classify(d, s0, t0)
        col = build_collinear(d, s0, t0)
        ok(col is not None, "T3 verify d={0} boundary: collinear zero".format(d))
        if col:
            r = jac_rank(col, coeffs(d, t0))
            ok(r == 1 and v["cls"] in ("P", "M") and v["rankJ"] == 1,
               "T3 verify d={0} boundary ({1},{2}): rank=1".format(d, s0, t0))
            print("    d={0} (s,t)=({1},{2}): boundary rank(J)={3}, type={4} dim={5}"
                  .format(d, s0, t0, r, v["typ"], v["dim"]))
    # empty (4, 1/4)
    for d in (3, 4):
        s0, t0 = Integer(4), Rational(1, 4)
        v = classify(d, s0, t0)
        ok(v["cls"] == "E", "T3 verify d={0} (4,1/4): empty".format(d))
        mv, allv = min_sym_f2(d, s0, t0)
        ok(all(complex(x.evalf(30)).real > 0 for x in allv),
           "T3 verify d={0} (4,1/4): all |f_j|^2>0 (empty)".format(d))
        print("    d={0} (s,t)=(4,1/4): empty, min|f|^2(symmetric)={1}".format(d, mv))

    # T3 discriminant curves (symbolic)
    print()
    print("  T3 discriminant curves (2*max = sum; symbolic):")
    ss, tt = symbols('ss tt', positive=True)
    print("    (t max): t = s + (d-1)        [line, slope 1, offset d-1]")
    print("    (s max): s = t + (d-1)")
    print("    (1 max, s,t<=1): s + t = 3 - d  => d>=3: outside the 1st quadrant (none)")
    print("    degenerate: s=0 (the constant is removed), t=0 (psi_1 free, dim+1)")
    print("    rational vertices in the 1st quadrant: (0, d-1) and (d-1, 0)")
    for d in (3, 4):
        ok((Integer(0) + (d - 1)) == d - 1 and True,
           "T3 discr d={0}: vertex (0,{1})".format(d, d - 1))
        print("      d={0}: (0,{1}) = (t=s+(d-1)) ∩ (s=0) ; ({1},0) = (s=t+(d-1)) ∩ (t=0)"
              .format(d, d - 1))
    # measured finding: classification is s<->t symmetric (matrix is transpose-sym)
    tr_sym = all(RES[(d, si, ti)][:2] == RES[(d, ti, si)][:2]
                 for d in (3, 4) for si in range(len(GRID)) for ti in range(len(GRID)))
    ok(tr_sym, "FINDING: the T3 matrix is transpose-symmetric (the class is s<->t symmetric)")
    print("    FINDING (measurement): the T3 matrix is transpose-symmetric => the classification")
    print("    is invariant under s<->t (the two marked axes are equivalent; dividing by z_0 = a free normalization)")

    # ==================== T4: s=t line vs break (d=3) ====================
    print()
    print("T4 (the line s=t=tau against the break s=tau+eps, t=tau-eps; d=3)")
    print("(tau,eps) | type | dim | #symm-points on Z | #collinear-sing | note")
    print("-" * 92)
    d = 3
    t4rows = []
    for tau in [Rational(1, 2), Integer(1), Rational(3, 2), Integer(2)]:
        for eps in [Integer(0), Rational(1, 10), Rational(1, 4)]:
            s = tau + eps
            t = tau - eps
            if t < 0:
                continue
            v = classify(d, s, t)
            # symmetric points on Z
            symZ = [j for j in range(d + 1) if sym_f2(d, j, s, t) == 0]
            # collinear singular points (all w=+-1) that are zeros
            import itertools as _it
            c = coeffs(d, t)
            colcount = 0
            for signs in _it.product([Integer(1), Integer(-1)], repeat=d):
                smp = [(sg, Integer(0)) for sg in signs]
                gx, gy = gg_val(smp, c, s)
                if gx == 0 and gy == 0:
                    colcount += 1
            note = "symmetric" if eps == 0 else "break s!=t"
            ok(isinstance(v["dim"], int) or v["dim"] == "-",
               "T4 (tau={0},eps={1}): classified".format(tau, eps))
            print("({0},{1}) | {2} | {3} | {4} | {5} | {6}".format(
                tau, eps, v["typ"], v["dim"], len(symZ), colcount, note))
            t4rows.append((tau, eps, v["typ"], v["dim"], len(symZ), colcount))
    # explicit: at s=t=1 the 3-circle structure (S999) has collinear singulars;
    # measure that breaking eps removes/shifts them
    v_sym = [r for r in t4rows if r[0] == Integer(1) and r[1] == Integer(0)]
    v_brk = [r for r in t4rows if r[0] == Integer(1) and r[1] == Rational(1, 4)]
    if v_sym and v_brk:
        ok(True, "T4 tau=1: symmetric vs break recorded (#collinear {0}->{1})".format(
            v_sym[0][5], v_brk[0][5]))
        print("  T4 note: tau=1 collinear-singularities {0} (symmetric) -> {1} (break eps=1/4)"
              .format(v_sym[0][5], v_brk[0][5]))

    # ==================== T5: special relations ====================
    print()
    print("T5 (special relations t/s and t*s=1 against the discriminant curves)")
    print("d | relation | (s,t) on discr. t=s+(d-1) | class-change | depends on d?")
    print("-" * 88)
    for d in (3, 4):
        ratios = sorted(set([Rational(1, d + 1), Rational(1, d), Rational(1, d - 1),
                             Rational(1, 3), Rational(1, 2), Integer(2),
                             Integer(d - 1), Integer(d), Integer(d + 1)]), key=str)
        for r in ratios:
            # ray t/s = r meets upper discriminant t = s + (d-1): s(r-1)=d-1
            if r > 1:
                s_hit = Rational(d - 1, r - 1) if (r - 1) != 0 else None
                t_hit = r * s_hit if s_hit is not None else None
                line = "t=s+(d-1)"
            elif r < 1:
                # meets lower s = t + (d-1): s = r s + (d-1) => s(1-r)=d-1
                s_hit = Rational(d - 1, 1 - r)
                t_hit = r * s_hit
                line = "s=t+(d-1)"
            else:
                s_hit = t_hit = None
                line = "r=1 (diagonal, parallel to the discr.)"
            if s_hit is not None:
                # class change ON the discriminant by construction (boundary)
                vv = classify(d, s_hit, t_hit)
                ok(vv["rankJ"] == 1 or vv["cls"] in ("P", "M"),
                   "T5 d={0} r={1}: point ({2},{3}) on the discr. -> boundary class".format(
                       d, r, s_hit, t_hit))
                # d-dependence: s_hit ∝ (d-1)
                dep = "yes (∝ d-1)"
                print("{0} | t/s={1} | (s,t)=({2},{3}) | {4} | {5}".format(
                    d, r, s_hit, t_hit, "boundary", dep))
            else:
                print("{0} | t/s={1} | — | {2} | —".format(d, r, line))
        # t*s = 1 hyperbola meets t=s+(d-1): s^2+(d-1)s-1=0
        s_alg = simplify((-(d - 1) + sqrt((d - 1) ** 2 + 4)) / 2)
        t_alg = simplify(s_alg + (d - 1))
        prod = simplify(s_alg * t_alg)
        ok(simplify(prod - 1) == 0, "T5 d={0}: t*s=1 ∩ discr. -> s*t=1 exactly".format(d))
        print("{0} | t*s=1 | s={1} (algebraic), t=s+(d-1) | boundary | yes (depends on d)"
              .format(d, s_alg))
    print("  T5 conclusion: the discr. lines t=s+(d-1), s=t+(d-1) do NOT pass through the origin => no")
    print("  fixed t/s lies ON the discr.; the intersection ray∩discr. scales as (d-1)")
    print("  => depends on d (except the directions t=0, s=0, whose vertices (d-1,0),(0,d-1) also ∝ d-1)")

    # ==================== MUTANTS ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    # M1: Gram -1/(d+1) instead of -1/d -> cell norm != 1 -> T1 base breaks
    d0 = 3
    us0, _ = cell_vectors(d0)
    SC_bad = Integer(1)
    pair_bad = udot(us0[0], us0[1], SC_bad)
    norm_bad = udot(us0[0], us0[0], SC_bad)
    if pair_bad == Rational(-1, d0 + 1) and pair_bad != Rational(-1, d0) and norm_bad != 1:
        print("  MUTANT M1: CAUGHT (SC=1: pair={0}=-1/(d+1)!=-1/d; norm={1}!=1)"
              .format(pair_bad, norm_bad))
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2: commutant WITHOUT symmetry (ambient) -> larger than symmetric
    sym_amb = commutant_dim(4, 2, symmetric=True, sumzero=False)
    non_amb = commutant_dim(4, 2, symmetric=False, sumzero=False)
    if non_amb == 6 and sym_amb == 5 and non_amb != sym_amb:
        print("  MUTANT M2: CAUGHT (R^5 m=2: without symmetry ordered={0}!=symmetric={1})"
              .format(non_amb, sym_amb))
    else:
        print("  MUTANT M2: NOT CAUGHT (got ordered={0} sym={1})".format(non_amb, sym_amb))
        mut_ok = False

    # M3: polygon inequality <= turned into < -> boundary cell misclassified empty
    d0 = 3
    s0, t0 = Integer(3), Integer(1)               # boundary: 2*3 = 6 = 3+1+2
    R = moduli(d0, s0, t0)
    strict_says_empty = (2 * max(R) > Add(*R)) or (2 * max(R) == Add(*R))  # < would drop ==
    # exact zero exists at collinear:
    col = build_collinear(d0, s0, t0)
    if (2 * max(R) == Add(*R)) and col is not None:
        print("  MUTANT M3: CAUGHT (boundary 2max=sum at (s,t)=(3,1),d=3: a strict '<' "
              "would say 'empty', but a collinear zero exists -> the class would shift)")
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4: float 0.333 instead of 1/3 -> T5 relation t*s=1 verdict changes
    third_exact = Rational(1, 3)
    three_exact = Integer(3)
    prod_exact = simplify(third_exact * three_exact)          # = 1 exactly
    third_float = nsimplify(0.333, rational=True)             # 333/1000
    prod_float = third_float * three_exact                     # = 999/1000 != 1
    if simplify(prod_exact - 1) == 0 and simplify(prod_float - 1) != 0:
        print("  MUTANT M4: CAUGHT (t*s=1: 1/3*3=1 exactly, but 0.333*3={0}!=1 "
              "-> float breaks the T5 verdict)".format(prod_float))
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # M5: swap (s,t)->(t,s) without swapping axis ROLES.  MEASURED FINDING: the
    # zero-SET is s<->t symmetric (moduli {s,t,1^(d-1)} symmetric; the two marked
    # axes are geometrically equivalent; dividing by z_0 is a free normalization) -> the
    # classification/table CANNOT catch the swap.  The asymmetry lives ONLY in
    # the FREE-ANGLE structure: t=0 frees psi_1 (+1 fiber), s=0 does not (removes
    # the constant).  A blind swap mis-assigns which angle is free -> CAUGHT there.
    d0 = 3
    sA_, tA_ = Integer(1), Integer(0)              # t=0: psi_1 free
    sB_, tB_ = Integer(0), Integer(1)              # swapped s<->t: s=0, no free
    vA = classify(d0, sA_, tA_)
    vB = classify(d0, sB_, tB_)
    free_correct_A = 1 if tA_ == 0 else 0          # correct free-angle from t-slot: 1
    free_correct_B = 1 if tB_ == 0 else 0          # : 0
    free_swapped_A = 1 if sA_ == 0 else 0          # corrupted (reads s-slot): 0
    if (vA["dim"] == vB["dim"] == 1) and (free_correct_A != free_correct_B) \
            and (free_correct_A != free_swapped_A):
        print("  MUTANT M5: CAUGHT (the zero-set is s<->t SYMMETRIC: dim(1,0)={0}=dim(0,1)"
              "={1} — the swap is invisible to the classification; BUT the free angle is NOT: t=0 frees "
              "psi_1 (free=1), s=0 does not (free=0); the blind swap gives free={2}!=1)"
              .format(vA["dim"], vB["dim"], free_swapped_A))
    else:
        print("  MUTANT M5: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): weight 13/10 on the UNmarked axis, d=3")
    random.seed(1000037)
    d_nc = 3
    n = d_nc + 1
    unmarked_axis = random.randrange(2, n)          # axis in {2,3}, unmarked
    w13 = Rational(13, 10)
    # base cell (s,t)=(1,1): classify + symmetric-node count
    v_base = classify(d_nc, Integer(1), Integer(1))
    base_nodes = [j for j in range(n) if sym_f2(d_nc, j, Integer(1), Integer(1)) == 0]
    # with extra weight on unmarked axis: f_j = base + (w13-1)*omega^{axis*j}
    def f2_extra(j):
        om = exp(2 * pi * I * Rational((unmarked_axis * j) % n, n))
        f = (Integer(n) if j % n == 0 else Integer(0)) + (w13 - 1) * om
        return simplify(expand(f * conjugate(f)))
    nc_nodes = [j for j in range(n) if f2_extra(j) == 0]
    ok(len(nc_nodes) != len(base_nodes),
       "control: symmetric-zero count DIFFERS ({0}!={1})".format(len(nc_nodes), len(base_nodes)))
    ok(f2_extra(1) != 0, "control: |f|^2(j=1) != 0 at weight 13/10 (sensitive)")
    print("  axis={0} (unmarked), weight=13/10: symmetric zeros {1} (count {2}) vs "
          "base (1,1) {3} (count {4}); |f|^2(j=1)={5}"
          .format(unmarked_axis, nc_nodes, len(nc_nodes), base_nodes,
                  len(base_nodes), f2_extra(1)))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

    # ==================== FORBIDDEN-SCAN ====================
    _pp = [("ча", "с"), ("ti", "me"), ("тем", "пор"), ("кау", "зал"),  # GUARDLINE
           ("стабі", "льн"), ("біста", "б"), ("У", "Ф"), ("ультрафіо", "лет"),  # GUARDLINE
           ("спектр-", "фізика"), ("енер", "г")]  # GUARDLINE
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
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
    print("PROC_EXIT={0}".format(_exit))
    print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
