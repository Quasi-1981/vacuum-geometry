# -*- coding: utf-8 -*-
# DIM: na (cell closure by dimension d=2..6; exact rational arithmetic; 0 handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting)
# ----------------------------------------------------------------------------
# CELL of dimension d:  d+1 unit vectors with pairwise inner product -1/d.
#   Realization in the hyperplane {sum x = 0} of R^{d+1}:  u_i = e_i - c,
#   c = 1/(d+1) * 1  (centroid).  |u_i|^2 = d/(d+1) (Euclid),  u_i.u_j = -1/(d+1).
#   The UNIT vectors are v_i = u_i/|u_i|; their inner products are RATIONAL:
#   v_i.v_i = 1, v_i.v_j = -1/d.  We never take a root: carry the scale
#   SC = (d+1)/d so that <a,b>_unit = SC * (a.b)_Euclid on the u-frame.  All
#   u_i, differences and duals are exact rational vectors in R^{d+1}.
#
# TASK 1 (translational closure): L = Z-span of {v_i - v_j} = Z-span{e_i - e_j}
#   = the A_d root lattice (rank d, in {sum=0} of Z^{d+1}).  Reduced basis
#   a_p = e_{p-1} - e_p, p=1..d.  Euclid Gram = Cartan(A_d) (2 on diag, -1 sub),
#   det = d+1; unit Gram = SC * Cartan, det = SC^d * (d+1).
#
# TASK 2 (two-sublattice): L u (L + u_0).  Shortest distance^2 between the two
#   sublattices = min unit-norm^2 over the coset (u_0 + L); z(d) = # coset
#   vectors on that shell; bond angles via cos = v_i.v_j (exact).
#
# TASK 3 (traversal of the simplest flat figure) -- CHOICE stamped here:
#   class = COMPOSITION OF TWO EDGE-REFLECTIONS at a cell vertex.  Two adjacent
#   bond vectors b_i,b_j (unit, b_i.b_j = -1/d); the reflection through the
#   hyperplane _|_ b then _|_ b' is a rotation by 2*theta, theta = angle(b,b'),
#   so cos(2*theta) = 2*cos^2(theta) - 1 = 2/d^2 - 1 (exact, rational).  Order =
#   finite iff 2*theta is a rational multiple of pi <=> (Niven) cos(2theta) in
#   {0,+-1/2,+-1}.  Report finite order / infinite, per d.
#
# TASK 4 (nearest-neighbour Bloch, bipartite): H(k)=[[0,f],[f*,0]],
#   f(k) = sum_i exp(2*pi*i <k,delta_i>), delta_i = bond vectors = the d+1 cell
#   vectors.  Symmetric points: the coset family k_j = j*w (j=0..d), w in the
#   {sum=0} hyperplane with w_i = (2i-d)/(2(d+1)) so <k_j,delta_i> = j*i/(d+1)
#   (mod a global phase).  |f|^2 = 0 at a symmetric point <=> the d+1 phases are
#   the full (d+1)-th roots of unity summed (Gauss sum = 0).  Count nodes per d.
#
# Discipline: 0 handles; sympy exact (Integer/Rational/roots of unity); mutant
# on every classification branch; FORBIDDEN-SCAN; bit-fence d=2 vs the program's
# known carve (triangular/honeycomb, 120 deg); STOP after the table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
from sympy import (Matrix, Integer, Rational, zeros, eye, ones, exp, I, pi,
                   simplify, Add, nsimplify, sqrt, conjugate, re, im)

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S956_run.log")
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

ASSERT_PASS = [0]
FAILS = [0]


def ok(cond, msg):
    if cond:
        ASSERT_PASS[0] += 1
    else:
        FAILS[0] += 1
        print("ASSERT-FAIL: " + msg)


# ==================== rational primitives ====================

def cell_vectors(d):
    """u_i = e_i - centroid in R^{d+1}, exact rational; and unit-scale SC."""
    n = d + 1
    c = ones(n, 1) * Rational(1, n)
    us = []
    for i in range(n):
        e = zeros(n, 1)
        e[i, 0] = Integer(1)
        us.append(e - c)
    SC = Rational(n, d)  # (d+1)/d : <.,.>_unit = SC * (.)_Euclid on u-frame
    return us, SC


def edot(a, b):
    return (a.T * b)[0, 0]


def udot(a, b, SC):
    return SC * edot(a, b)


def cartan_Ad(d):
    """Cartan matrix of A_d: 2 on diagonal, -1 on first off-diagonals."""
    M = zeros(d, d)
    for i in range(d):
        M[i, i] = Integer(2)
        if i + 1 < d:
            M[i, i + 1] = Integer(-1)
            M[i + 1, i] = Integer(-1)
    return M


def root_basis(d):
    """a_p = e_{p-1} - e_p in R^{d+1}, p=1..d (basis of A_d root lattice)."""
    n = d + 1
    out = []
    for p in range(d):
        v = zeros(n, 1)
        v[p, 0] = Integer(1)
        v[p + 1, 0] = Integer(-1)
        out.append(v)
    return out


def coset_search(d, us, SC, box=2):
    """Enumerate coset (u_0 + L) with L-coords in [-box,box]^d; return
    (min unit-norm^2, list of minimizing vectors as R^{d+1})."""
    n = d + 1
    basis = root_basis(d)
    u0 = us[0]
    best = None
    minima = []
    # iterate integer combos of the d basis vectors
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


# ==================== Niven / rotation order ====================

NIVEN = {Rational(0), Rational(1, 2), Rational(-1, 2), Integer(1), Integer(-1)}


def rotation_cos2theta(d):
    return Rational(2, d * d) - 1  # 2/d^2 - 1, exact


def order_from_cos(c2):
    """Finite order of a planar rotation with cos(angle)=c2, via Niven.
    Returns integer order or None (infinite)."""
    if c2 not in NIVEN:
        return None
    # map cos(angle) -> angle as fraction of 2pi, then order
    table = {Integer(1): Integer(1), Rational(1, 2): Integer(6),
             Rational(0): Integer(4), Rational(-1, 2): Integer(3),
             Integer(-1): Integer(2)}
    return table[c2]


# ==================== Bloch f(k) at symmetric points ====================

def bloch_symmetric_point(d, j):
    """f(k_j) = sum_{i=0}^{d} exp(2 pi I * j*i/(d+1)); exact.  Also returns the
    concrete rational w-vector so <k_j,delta_i> = j*i/(d+1) is verified."""
    n = d + 1
    # concrete k_j = j*w with w_i = (2i-d)/(2(d+1)); delta_i = u_i
    w = zeros(n, 1)
    for i in range(n):
        w[i, 0] = Rational(2 * i - d, 2 * n)
    kj = Integer(j) * w
    # verify <k_j, delta_i> == j*i/(d+1) mod 1 (up to the global -j*d/(2(d+1)))
    us, _ = cell_vectors(d)
    ok(edot(w, ones(n, 1)) == 0, "w in {sum=0} hyperplane")
    phases_ok = True
    for i in range(n):
        ph = edot(kj, us[i])
        want = Rational(j * (2 * i - d), 2 * n)  # = j*i/n - j*d/(2n)
        if ph != want:
            phases_ok = False
    ok(phases_ok, "phases <k_j,delta_i> concrete-consistent (d={0},j={1})".format(d, j))
    # f up to global phase = sum_{i=0}^{n-1} omega^i, omega = exp(2 pi I j/n).
    # Geometric series (EXACT for all n): omega^n = exp(2 pi I j) = 1, so if
    # omega != 1  ->  sum = (omega^n - 1)/(omega - 1) = 0;  if omega == 1 -> n.
    omega = exp(2 * pi * I * Rational(j % n, n))
    num = simplify(omega ** n - 1)          # exact 0 (omega^n = 1)
    ok(num == 0, "omega^n = 1 exact (d={0},j={1})".format(d, j))
    if (j % n) == 0:
        f = Integer(n)
    else:
        # omega != 1 and numerator == 0  =>  sum is exactly 0
        f = Integer(0)
    # cross-check: direct sum evaluated at high precision must match f (the
    # exact verdict rests on the geom-series identity above; this is evidence)
    direct = Add(*[exp(2 * pi * I * Rational((j * i) % n, n)) for i in range(n)])
    dval = complex(direct.evalf(30))
    ok(abs(dval - complex(f)) < 1e-20,
       "direct sum ~= geom-series numerically (d={0},j={1})".format(d, j))
    f2 = f * conjugate(f)
    return f, f2


# ============================================================================
#                                   RUN
# ============================================================================

print("=" * 74)
print("W34-leg-1: cell closure d=2..6 (A_d), rational arithmetic")
print("=" * 74)

rows = []
DS = [2, 3, 4, 5, 6]

for d in DS:
    n = d + 1
    us, SC = cell_vectors(d)

    # ---- sanity: cell inner products ----
    ok(all(udot(us[i], us[i], SC) == 1 for i in range(n)),
       "cell unit norm =1 (d={0})".format(d))
    ok(all(udot(us[i], us[j], SC) == Rational(-1, d)
           for i in range(n) for j in range(n) if i != j),
       "cell pairwise cos = -1/d (d={0})".format(d))

    # ---- TASK 1: translational closure ----
    basis = root_basis(d)
    G_eucl = Matrix(d, d, lambda i, j: edot(basis[i], basis[j]))
    G_unit = SC * G_eucl
    rank = G_eucl.rank()
    detE = G_eucl.det()
    detU = G_unit.det()
    ok(rank == d, "L rank = d (d={0})".format(d))
    ok(G_eucl == cartan_Ad(d), "Euclid Gram = Cartan(A_d) (d={0})".format(d))
    ok(detE == d + 1, "det Euclid Gram = d+1 (d={0})".format(d))
    ok(detU == SC ** d * (d + 1), "det unit Gram = SC^d (d+1) (d={0})".format(d))

    # ---- TASK 2: two-sublattice ----
    dmin, minima = coset_search(d, us, SC, box=2)
    z = len(minima)
    # bond vectors = the d+1 cell vectors; angle cos between distinct bonds
    bond_cos = udot(us[0], us[1], SC)
    ok(dmin == 1, "shortest sublattice dist^2 = 1 (d={0})".format(d))
    ok(z == d + 1, "z(d) = d+1 (d={0}) got {1}".format(d, z))
    ok(bond_cos == Rational(-1, d), "bond cos = -1/d (d={0})".format(d))

    # ---- TASK 3: traversal order ----
    c2 = rotation_cos2theta(d)
    order = order_from_cos(c2)
    ordstr = str(order) if order is not None else "inf"

    # ---- TASK 4: Bloch nodes at symmetric points ----
    node_js = []
    f_gamma2 = None
    for j in range(0, n):
        f, f2 = bloch_symmetric_point(d, j)
        if j == 0:
            f_gamma2 = f2
        if f2 == 0:
            node_js.append(j)
    n_nodes = len(node_js)
    ok(f_gamma2 == (d + 1) ** 2, "Gamma (j=0) |f|^2 = (d+1)^2 (d={0})".format(d))
    ok(n_nodes == d, "node count = d (d={0}) got {1}".format(d, n_nodes))

    rows.append(dict(d=d, gram="SC*Cartan(A_d)", SC=SC, detE=detE, detU=detU,
                     z=z, bondcos=bond_cos, order=ordstr, c2=c2,
                     nodes=n_nodes, dmin=dmin))

# ==================== TABLE ====================
print()
print("COLUMNS (raw):")
hdr = ("d | generator Gram | SC=(d+1)/d | det(Euclid) | det(unit) | "
       "z(d) | cos(bond) | traversal(order) | cos2θ | null-modes(symm) | dist²")
print(hdr)
print("-" * len(hdr))
for r in rows:
    print("{d} | {gram} | {SC} | {detE} | {detU} | {z} | {bondcos} | "
          "{order} | {c2} | {nodes} | {dmin}".format(**r))

# ==================== BIT-FENCE d=2 (known carve) ====================
print()
print("BIT-FENCE d=2 (cross-check against the program's stamps):")
r2 = rows[0]
ok(r2["detE"] == 3, "d=2: A_2 root lattice det = 3 (triangular/honeycomb)")
ok(r2["bondcos"] == Rational(-1, 2), "d=2: bond cos = -1/2 => 120° (honeycomb)")
ok(r2["z"] == 3, "d=2: z=3 (honeycomb coordination)")
ok(r2["nodes"] == 2, "d=2: 2 null-modes (K,K' Dirac points)")
ok(r2["order"] == "3", "d=2: traversal order=3 (finite)")
print("  det(A_2)=3 · cos=−1/2(120°) · z=3 · nodes=2(K,K') · order=3 — all agree")

# ==================== MUTANTS (one per classification branch) ====================
print()
print("MUTANTS (one per classification branch):")
mut_ok = True

# br-rank (Task1): a rank-deficient generator set must NOT close to rank d
d0 = 3
basis0 = root_basis(d0)
bad_basis = basis0[:-1] + [basis0[0]]  # duplicate -> rank d-1
Gbad = Matrix(d0, d0, lambda i, j: edot(bad_basis[i], bad_basis[j]))
if cartan_Ad(d0).rank() == d0 and Gbad.rank() < d0:
    print("  MUTANT br-rank: CAUGHT (good basis rank=d; duplicated gen rank<d)")
else:
    print("  MUTANT br-rank: NOT CAUGHT"); mut_ok = False

# br-shell (Task2): a coset vector off the minimal shell must be excluded from z
usm, SCm = cell_vectors(d0)
dminm, minm = coset_search(d0, usm, SCm, box=2)
# a genuine second-shell coset vector: u_0 + (e0-e1)+(e1-e2)... take u_0 + 2*a_1
far = usm[0] + Integer(2) * root_basis(d0)[0]
far_norm = udot(far, far, SCm)
if far_norm > dminm and all(udot(m, m, SCm) == dminm for m in minm):
    print("  MUTANT br-shell: CAUGHT (far coset vec norm²={0}>{1}=min; not in z)".format(far_norm, dminm))
else:
    print("  MUTANT br-shell: NOT CAUGHT"); mut_ok = False

# br-order (Task3): a non-Niven cos(2θ) must be rejected as INFINITE
c2_d3 = rotation_cos2theta(3)  # = -7/9
c2_d2 = rotation_cos2theta(2)  # = -1/2
if (order_from_cos(c2_d3) is None) and (order_from_cos(c2_d2) == 6 or order_from_cos(c2_d2) is not None):
    # d=2 finite, d=3 infinite
    if order_from_cos(c2_d2) is not None and order_from_cos(c2_d3) is None:
        print("  MUTANT br-order: CAUGHT (d=2 cos=-1/2 finite; d=3 cos=-7/9 infinite)")
    else:
        print("  MUTANT br-order: NOT CAUGHT"); mut_ok = False
else:
    print("  MUTANT br-order: NOT CAUGHT"); mut_ok = False

# br-node (Task4): the Gamma point (j=0) must NOT be a node
_, f2_gamma = bloch_symmetric_point(2, 0)
_, f2_K = bloch_symmetric_point(2, 1)
if f2_gamma != 0 and f2_K == 0:
    print("  MUTANT br-node: CAUGHT (Γ |f|²={0}≠0 no node; K |f|²=0 node)".format(f2_gamma))
else:
    print("  MUTANT br-node: NOT CAUGHT"); mut_ok = False

# ==================== summary ====================
print()
print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

_pieces = [("з", "акон"), ("к", "анал"), ("мех", "анізм"), ("зл", "іпок"),
           ("пер", "егин"), ("конд", "енсат"), ("мат", "ерія"), ("ен", "ергія"),
           ("рез", "онанс"), ("тр", "іщина"), ("г", "учн")]
_words = ["".join(ab).casefold() for ab in _pieces]
_src = open(__file__, "r", encoding="utf-8").read().casefold()
_logf.flush()
_logtxt = "".join(_tee.chunks).casefold()
_hits = 0
for _w2 in _words:
    if _w2 in _src:
        _hits += 1
    if _w2 in _logtxt:
        _hits += 1
print("FORBIDDEN-SCAN: hits={0}".format(_hits))

_exit = 0
if _hits > 0:
    _exit = 1
if FAILS[0] > 0:
    _exit = 1
if not mut_ok:
    _exit = 1
print("EXIT={0}".format(_exit))
sys.exit(_exit)
