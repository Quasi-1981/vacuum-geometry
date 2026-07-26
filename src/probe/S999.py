# -*- coding: utf-8 -*-
# DIM: na (W40: structure of f near the zero-set — linear order + eigenvalues;
#          exact sympy; 0 handles).
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — blind probe, reuse S956/S998
# ----------------------------------------------------------------------------
# CELL d (as S956/S998): d+1 unit axes, pairwise product -1/d; realization
#   u_i = e_i - c in {sum x = 0} of R^{d+1}; SC = (d+1)/d rational scale;
#   closure = A_d root lattice; bipartite; bond vectors delta_i = the d+1 axes.
# BZ REDUCTION (as S998): psi_i = <k, delta_i - delta_0> (i = 1..d) coordinates
#   on the BZ torus; f = z_0 * g with |z_0| = 1, so f = 0 <=> g = 0, where
#   g_t(psi) = t + sum_{i=1..d} w_i, w_i = exp(2 pi I psi_i)  (weight t on
#   axis 0, weight 1 on the rest; base t = 1).
# METRIC CONVENTION (Cartesian): alpha_i = delta_i - delta_0; Gram matrix
#   G_ij = <alpha_i, alpha_j>_unit = SC*(1 + [i==j])  (SPD, rational; verified
#   from the S956 primitives, not assumed).  For any orthonormal frame A with
#   A*A^T = G (here: exact Cholesky), J_k = J_psi * A, hence
#     J_k J_k^T = J_psi G J_psi^T   and   ker(J_k^T J_k) = A^{-1} ker(J_psi):
#   the kernel-vs-tangent comparison reduces to psi-space (A invertible), and
#   the nonzero eigenvalues of J^T J equal the eigenvalues of J_psi G J_psi^T.
# J CONVENTION: the 2 pi factor is dropped uniformly (J := d(Re g, Im g)/dk / 2pi);
#   eigenvalues are reported in (2 pi)^2 units; rank/kernel/equality verdicts
#   are invariant under this uniform scale.  f vs g: J_f = R(phase) J_g with R
#   orthogonal, so J^T J, |det| and all verdicts are identical for f and g.
# TANGENTS are measured from EXPLICIT parametrizations (never from ker J — no
#   circularity): d=3 circles (z,-z,-1)/(z,-1,-z)/(-1,z,-z) with an exact
#   algebraic certificate (1+w1)(1+w2)(1+w3) = 0 on Z(t=1); d=4 two-parameter
#   family w1, w2 free + pair (w3, w4) = S*(1/2 +- I*sqrt(1/|S|^2 - 1/4)),
#   S = -(1+w1+w2); conj-symmetric slice through symmetric points; derivative
#   direction d psi_i  proportional to  dw_i/(I w_i)  (real, uniform scale).
# EXACT zero-verdicts: expand/simplify chain, fallback minimal_polynomial
#   (algebraic numbers); float is FORBIDDEN in verdicts (mutant M3 shows why).
# PARALLEL HARNESS: independent (d, point) subtasks run in a
#   multiprocessing.Pool(18); workers RETURN serialized results (strings/bools),
#   the master prints in FIXED task order => the log is bit-reproducible.
#   All assert accounting (ok) + mutants live in the MASTER process.
# Discipline: 0 handles; mutants M1..M5 CAUGHT; seeded negative control;
#   FORBIDDEN-SCAN via tools.fence_scan (source + log); bit-fence T1 vs S998
#   run.log read verbatim; STOP after the tables.
# ============================================================================

import sys
import os
import random
from multiprocessing import Pool
from sympy import (Matrix, Integer, Rational, zeros, ones, eye, exp, I, pi,
                   simplify, expand, Add, sqrt, conjugate, symbols, solve,
                   cos, sin, re, im, Symbol, Poly, roots, diff, radsimp,
                   minimal_polynomial)

_HERE = os.path.dirname(os.path.abspath(__file__))

# ==================== exact zero-verdict machinery ====================

_mx = Symbol('_mx')


def is_zero_exact(e):
    """Exact verdict e == 0 for sympy objects.  Chain: expand+radsimp+simplify;
    explicit .is_zero; fallback minimal_polynomial (algebraic numbers:
    minpoly(0) = _mx).  NO float anywhere."""
    e2 = simplify(radsimp(expand(e)))
    if e2 == 0:
        return True
    if e2.free_symbols:
        return False        # symbolic: simplify said not identically zero
    z = e2.is_zero
    if z is not None:
        return bool(z)
    try:
        return minimal_polynomial(e2, _mx) == _mx
    except Exception:
        return False


def eq_verdict(diff_expr):
    return "yes" if is_zero_exact(diff_expr) else "no"


def _izf(v):
    return is_zero_exact(v)


# ==================== check collector (worker-side; ok() in master) ====================

_CHK = []


def chk(cond, msg):
    _CHK.append((bool(cond), msg))


def drain():
    out = list(_CHK)
    del _CHK[:]
    return out


# ==================== rational primitives (verbatim from S956/S998) ====================

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


def g_val(sample, c, t):
    """g = t + sum_i c_i * w_i at sample [(x_i, y_i)]; returns (Re, Im) exact."""
    gx = t + Add(*[c[i] * sample[i][0] for i in range(len(sample))])
    gy = Add(*[c[i] * sample[i][1] for i in range(len(sample))])
    return simplify(gx), simplify(gy)


def jac_psi(sample, c):
    """Real 2 x d Jacobian of (Re g, Im g) wrt psi at the sample:
    column i = c_i * (-y_i, x_i)   (2 pi dropped uniformly — see convention)."""
    M = zeros(2, len(sample))
    for i, (x, y) in enumerate(sample):
        M[0, i] = -c[i] * y
        M[1, i] = c[i] * x
    return M


def jac_rank(sample, c):
    return jac_psi(sample, c).rank()


def build_sample_m1(d, t):
    """Exact zero of g_t = t + sum_{i=1..d} w_i for 0 <= t < d (verbatim S998)."""
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
    """d=2, m=1, t>0: full exact solve of t + w_1 + w_2 = 0 (verbatim S998)."""
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
    """Exact verdict on {g_t = 0} (m=1) in the d-torus BZ (verbatim S998,
    asserts -> chk)."""
    c = [Integer(1)] * d
    if t > d:
        chk(t - d > 0, "E-cert: t-d>0 (d={0}, t={1})".format(d, t))
        return dict(typ="empty", dim="-", npts=0, rank="-",
                    note="|g| >= t-d = {0} > 0".format(t - d))
    if t == d:
        sample = [(Integer(-1), Integer(0))] * d
        gx, gy = g_val(sample, c, t)
        chk(gx == 0 and gy == 0, "P: g(all psi=1/2) = 0 exact (d={0}, t={1})".format(d, t))
        r = jac_rank(sample, c)
        chk(r == 1, "P: rank(J)=1 at merge point (d={0})".format(d))
        return dict(typ="point", dim=0, npts=1, rank=r,
                    note="unique psi=(1/2,..,1/2): Re g >= t-d=0, equality <=> all w=-1")
    sample = build_sample_m1(d, t)
    chk(sample is not None, "M: sample exists (d={0}, t={1})".format(d, t))
    gx, gy = g_val(sample, c, t)
    chk(gx == 0 and gy == 0, "M: g(sample) = 0 exact (d={0}, t={1})".format(d, t))
    r = jac_rank(sample, c)
    chk(r == 2, "M: rank(J)=2 at regular sample (d={0}, t={1})".format(d, t))
    dim = d - r
    if d == 2:
        pts = solve_points_d2_m1(t)
        npts = len(pts)
        chk(npts == 2, "d=2: exactly 2 solved points (t={0}) got {1}".format(t, npts))
        note = "; ".join("(x={0},y={1})".format(p[0], p[1]) for p in pts)
        return dict(typ="points", dim=0, npts=npts, rank=r, note="solve: " + note)
    return dict(typ="variety", dim=dim, npts="-", rank=r,
                note="dim = d - rank(J) = {0}".format(dim))


def rowfmt(d, tv, v):
    return "{0} | {1} | {2} | {3} | {4} | {5} | {6}".format(
        d, tv, v["typ"], v["dim"], v["npts"], v["rank"], v["note"])


# ==================== metric machinery (exact) ====================

def gram_alpha(d):
    """G_ij = <alpha_i, alpha_j>_unit, alpha_i = u_i - u_0 (i=1..d).  Exact."""
    us, SC = cell_vectors(d)
    al = [us[i] - us[0] for i in range(1, d + 1)]
    G = Matrix(d, d, lambda i, j: udot(al[i], al[j], SC))
    chk(G == SC * (eye(d) + ones(d, d)),
        "Gram(alpha) = SC*(I + 1*1^T) from the S956 primitives (d={0})".format(d))
    return G, SC


G2, SC2 = gram_alpha(2)
G3, SC3 = gram_alpha(3)
G4, SC4 = gram_alpha(4)
A2 = G2.cholesky()
chk(simplify(A2 * A2.T - G2) == zeros(2, 2), "A2*A2^T = G2 (exact Cholesky)")
_IMPORT_CHKS = drain()


def eig2(K):
    """Exact eigenvalues of symmetric 2x2 + discriminant (equality <=> disc = 0)."""
    tr = simplify(K[0, 0] + K[1, 1])
    dt = simplify(K.det())
    disc = simplify(expand(tr ** 2 - 4 * dt))
    s = sqrt(disc)
    l1 = simplify((tr + s) / 2)
    l2 = simplify((tr - s) / 2)
    return l1, l2, disc, tr, dt


def kfmt(K):
    return "[[{0}, {1}], [{2}, {3}]]".format(K[0, 0], K[0, 1], K[1, 0], K[1, 1])


def dlog_dir(w, param):
    """d psi direction component: dw/dparam / (I w)  (real; uniform scale
    inside one family — span/kernel comparisons invariant)."""
    return diff(w, param) / (I * w)


def sym_sample(d, j):
    """Symmetric coset point k_j (S956 TASK 4): psi_i = j*i/(d+1) mod 1."""
    n = d + 1
    return [(cos(2 * pi * Rational((i * j) % n, n)),
             sin(2 * pi * Rational((i * j) % n, n))) for i in range(1, n)]


COMP_TANGENTS = {"C1": (1, 1, 0), "C2": (1, 0, 1), "C3": (0, 1, 1)}


def comps_through(ws):
    out = []
    if is_zero_exact(ws[2] + 1) and is_zero_exact(ws[0] + ws[1]):
        out.append("C1")
    if is_zero_exact(ws[1] + 1) and is_zero_exact(ws[0] + ws[2]):
        out.append("C2")
    if is_zero_exact(ws[0] + 1) and is_zero_exact(ws[1] + ws[2]):
        out.append("C3")
    return out


def t_point_data(ws, cvec, G, d, tangents, label, t=None):
    """Full Jacobian structure at the point: rank, ker-vs-tangent, transverse eigenvalues.
    Returns (table row, dict of raw fields); checks go through chk."""
    if t is None:
        t = Integer(1)
    sample = [(simplify(re(w)), simplify(im(w))) for w in ws]
    gx, gy = g_val(sample, cvec, t)
    chk(is_zero_exact(gx) and is_zero_exact(gy), "{0}: point on Z".format(label))
    J = jac_psi(sample, cvec)
    r = J.rank(iszerofunc=_izf)
    K = simplify(expand(J * G * J.T))
    l1, l2, disc, tr, dt = eig2(K)
    tgs = [Matrix(list(tg)) for tg in tangents]
    for tg in tgs:
        prod = simplify(J * tg)
        chk(all(is_zero_exact(prod[iq]) for iq in range(2)),
            "{0}: tangent in ker J".format(label))
    span_rank = Matrix.hstack(*tgs).rank(iszerofunc=_izf) if tgs else 0
    ker_dim = d - r
    if r == 2 and span_rank == ker_dim:
        kerv = "yes (ker = span of tangents, dim {0})".format(ker_dim)
        everd = eq_verdict(disc)
    elif r < 2 and span_rank == ker_dim:
        kerv = "no: singular (rank {0}, ker dim {1} = the plane of {2} tangent branches)".format(
            r, ker_dim, len(tgs))
        everd = eq_verdict(disc) + " (degenerate: rank<2)"
    else:
        kerv = "MISMATCH (span {0} != ker {1})".format(span_rank, ker_dim)
        chk(False, "{0}: ker/tangents did not match".format(label))
        everd = eq_verdict(disc)
    line = ("{0} | rank={1} | ker=tangents: {2} | K={3} | eig=({4}, {5}) | equal: {6}"
            .format(label, r, kerv, kfmt(K), l1, l2, everd))
    return line, dict(rank=r, everd=everd, l1=str(l1), l2=str(l2))


# ==================== families (module level, cheap symbolic builds) ====================

_aa, _bb = symbols('_aa _bb', real=True)
_p = symbols('_p', real=True)
_lam = symbols('_lam')


def W_re(s):
    return (1 - s ** 2) / (1 + s ** 2)


def W_im(s):
    return 2 * s / (1 + s ** 2)


w1s = W_re(_aa) + I * W_im(_aa)
w1sb = W_re(_aa) - I * W_im(_aa)
w2s = W_re(_bb) + I * W_im(_bb)
w2sb = W_re(_bb) - I * W_im(_bb)
Sx = -(1 + w1s + w2s)
Sxb = -(1 + w1sb + w2sb)
R2ab = simplify(expand(Sx * Sxb))
tau_ab = sqrt(1 / R2ab - Rational(1, 4))
w3s = Sx * (Rational(1, 2) + I * tau_ab)
w4s = Sx * (Rational(1, 2) - I * tau_ab)
FAM4 = [w1s, w2s, w3s, w4s]

w1p = cos(2 * pi * _p) + I * sin(2 * pi * _p)
w1pb = cos(2 * pi * _p) - I * sin(2 * pi * _p)


# ==================== worker tasks (pure functions of spec) ====================

def _task_t1(d):
    v = classify_m1(d, Integer(1))
    myrow = rowfmt(d, Integer(1), v)
    line = ("  d={0}: {1} | dim={2} | #points={3} | rank={4}  <- matches S998 run.log"
            .format(d, v["typ"], v["dim"], v["npts"], v["rank"]))
    return [line], {"row": myrow}


def _t2_point_machinery(idx, t):
    pts = solve_points_d2_m1(t)
    chk(len(pts) == 2, "exactly 2 zeros (d=2, t={0})".format(t))
    (px1, py1, px2, py2) = pts[idx]
    sample = [(px1, py1), (px2, py2)]
    W1 = px1 + I * py1
    W2 = px2 + I * py2
    gx, gy = g_val(sample, [Integer(1)] * 2, t)
    chk(gx == 0 and gy == 0, "point on Z (t={0})".format(t))
    # (i) exact linear order: g(psi0 + eps q) = c0 + c1 eps + O(eps^2)
    _epsq, _q1, _q2 = symbols('_epsq _q1 _q2', real=True)
    gser = t + W1 * exp(2 * pi * I * _epsq * _q1) + W2 * exp(2 * pi * I * _epsq * _q2)
    c0 = simplify(gser.subs(_epsq, 0))
    c1 = simplify(diff(gser, _epsq).subs(_epsq, 0))
    chk(c0 == 0, "c0 = 0 (the zeroth order vanishes)")
    Jp = jac_psi(sample, [Integer(1), Integer(1)])
    Jq = Jp * Matrix([_q1, _q2])
    chk(is_zero_exact(re(c1) - 2 * pi * Jq[0]) and is_zero_exact(im(c1) - 2 * pi * Jq[1]),
        "c1 == 2 pi (J q) for the symbolic direction (L = the linear form J)")
    # (ii) J in the Cartesian frame, det, eigenvalues
    Jk = simplify(Jp * A2)
    chk(simplify(Jk * Jk.T - Jp * G2 * Jp.T) == zeros(2, 2),
        "J_k J_k^T == J_psi G J_psi^T (cross-check of the frame)")
    detJk = simplify(Jk.det())
    JtJ = simplify(Jk.T * Jk)
    l1, l2, disc, tr, dt = eig2(JtJ)
    verd = eq_verdict(disc)
    chk(is_zero_exact(dt - detJk ** 2), "det(J^T J) = (det J_k)^2")
    row = " | ".join(("w=({0}, {1})".format(W1, W2),
                      "[{0}, {1}; {2}, {3}]".format(Jk[0, 0], Jk[0, 1], Jk[1, 0], Jk[1, 1]),
                      str(detJk), "{0}, {1}".format(l1, l2), verd))
    return [row], {"verd": verd}


def _task_t2(idx):
    return _t2_point_machinery(idx, Integer(1))


def _task_t51(idx):
    return _t2_point_machinery(idx, Rational(11, 10))


def _task_t3u():
    lines = []
    _u = symbols('_u', real=True)
    xu = (1 - _u ** 2) / (1 + _u ** 2)
    yu = 2 * _u / (1 + _u ** 2)
    sample_u = [(xu, yu), (-xu, -yu), (Integer(-1), Integer(0))]
    gx_u, gy_u = g_val(sample_u, [Integer(1)] * 3, Integer(1))
    chk(simplify(gx_u) == 0 and simplify(gy_u) == 0, "T3(a): g == 0 identically on C1(u)")
    Ju = jac_psi(sample_u, [Integer(1)] * 3)
    minor02 = simplify(Ju.extract([0, 1], [0, 2]).det())
    chk(simplify(minor02 - yu) == 0, "T3(a): minor (col.1,3) = y1(u) => rank=2 at u != 0")
    tg_u = Matrix([1, 1, 0])
    prod_u = simplify(Ju * tg_u)
    chk(prod_u == zeros(2, 1), "T3(a): tangent (1,1,0) in ker J identically in u")
    Ku = simplify(expand(Ju * G3 * Ju.T))
    l1u, l2u, discu, tru, dtu = eig2(Ku)
    chk(simplify(tru - Rational(16, 3)) == 0, "T3(a): tr K = 16/3 identically")
    chk(simplify(discu - (Rational(16, 3) * xu) ** 2) == 0,
        "T3(a): disc = ((16/3) x1(u))^2 identically")
    eqsol = solve(xu, _u)
    lines.append("  rank=2 (u outside the collinear intersections of the circles); "
                 "ker = span(1,1,0) = tangent: yes")
    lines.append("  K(u) = {0}".format(kfmt(Ku)))
    lines.append("  tr=16/3; disc=((16/3)x1(u))^2 => eigenvalues = 8/3 +- (8/3)|x1(u)|")
    lines.append("  equal only at x1(u)=0, i.e. u in {0} (w1 = +-I); GENERICALLY: no"
                 .format(eqsol))
    return lines, {}


def _task_t3w():
    w_half = Rational(3, 5) + I * Rational(4, 5)
    line, data = t_point_data([w_half, -w_half, Integer(-1)],
                              [Integer(1)] * 3, G3, 3, [(1, 1, 0)],
                              "C1 u=1/2 w1=(3/5,4/5)")
    return [line], data


def _task_t3j(j):
    sm = sym_sample(3, j)
    ws = [sm[iq][0] + I * sm[iq][1] for iq in range(3)]
    comps = comps_through(ws)
    chk(len(comps) >= 1, "k_{0}: point on some component".format(j))
    tangents = [COMP_TANGENTS[cc] for cc in comps]
    line, data = t_point_data(ws, [Integer(1)] * 3, G3, 3, tangents,
                              "k_{0} (branches: {1})".format(j, ",".join(comps)))
    return [line], data


def _task_t4g(a0s, b0s, label):
    a0 = Rational(a0s)
    b0 = Rational(b0s)
    sub0 = {_aa: a0, _bb: b0}
    R2v = simplify(R2ab.subs(sub0))
    chk(R2v > 0 and R2v < 4, "{0}: 0 < |S|^2 = {1} < 4 (a regular pair)".format(label, R2v))
    ws = [simplify(w.subs(sub0)) for w in FAM4]
    for w in ws:
        chk(is_zero_exact(w * conjugate(w) - 1), "{0}: |w|=1".format(label))
    Ta = Matrix([simplify(dlog_dir(w, _aa).subs(sub0)) for w in FAM4])
    Tb = Matrix([simplify(dlog_dir(w, _bb).subs(sub0)) for w in FAM4])
    for TT, nm in ((Ta, "T_a"), (Tb, "T_b")):
        chk(all(is_zero_exact(im(TT[iq])) for iq in range(4)),
            "{0}: {1} is real".format(label, nm))
    Ta_r = Matrix([simplify(re(Ta[iq])) for iq in range(4)])
    Tb_r = Matrix([simplify(re(Tb[iq])) for iq in range(4)])
    chk(Matrix.hstack(Ta_r, Tb_r).rank(iszerofunc=_izf) == 2,
        "{0}: tangents T_a, T_b independent".format(label))
    line, data = t_point_data(ws, [Integer(1)] * 4, G4, 4,
                              [tuple(Ta_r), tuple(Tb_r)], label)
    return [line], data


def _task_t4j(j):
    sm = sym_sample(4, j)
    ws = [simplify(sm[iq][0] + I * sm[iq][1]) for iq in range(4)]
    # tangent A: the conj-symmetric slice psi = (p, q, -q, -p);
    # slice identity: g|slice = 1 + 2cos(2pi p) + 2cos(2pi q) (real)
    chk(is_zero_exact(1 + 2 * cos(2 * pi * Rational(j, 5))
                      + 2 * cos(2 * pi * Rational(2 * j, 5))),
        "k_{0}: the slice constraint holds".format(j))
    qp = simplify(-sin(2 * pi * Rational(j, 5)) / sin(2 * pi * Rational(2 * j, 5)))
    TA = Matrix([1, qp, -qp, -1])
    # tangent B: w2 frozen, w1 free, (w3,w4) = S(1/2 +- I tau)
    z2 = ws[1]
    z2b = conjugate(z2)
    Sp = -(1 + w1p + z2)
    Spb = -(1 + w1pb + z2b)
    R2p = simplify(expand(Sp * Spb))
    taup = sqrt(1 / R2p - Rational(1, 4))
    epsj = None
    for etry in (Integer(1), Integer(-1)):
        w3p_try = Sp * (Rational(1, 2) + I * etry * taup)
        if is_zero_exact(simplify(w3p_try.subs(_p, Rational(j, 5))) - ws[2]):
            epsj = etry
            break
    chk(epsj is not None, "k_{0}: the eps branch found for w3".format(j))
    if epsj is None:
        epsj = Integer(1)
    w3p = Sp * (Rational(1, 2) + I * epsj * taup)
    w4p = Sp * (Rational(1, 2) - I * epsj * taup)
    chk(is_zero_exact(simplify(w4p.subs(_p, Rational(j, 5))) - ws[3]),
        "k_{0}: the w4 branch is consistent".format(j))
    d3v = simplify(dlog_dir(w3p, _p).subs(_p, Rational(j, 5)))
    d4v = simplify(dlog_dir(w4p, _p).subs(_p, Rational(j, 5)))
    chk(is_zero_exact(im(d3v)) and is_zero_exact(im(d4v)),
        "k_{0}: the T_B components are real".format(j))
    TB = Matrix([2 * pi, 0, simplify(re(d3v)), simplify(re(d4v))])
    chk(Matrix.hstack(TA, TB).rank(iszerofunc=_izf) == 2,
        "k_{0}: T_A, T_B independent".format(j))
    line, data = t_point_data(ws, [Integer(1)] * 4, G4, 4,
                              [tuple(TA), tuple(TB)], "k_{0} (d=4)".format(j))
    return [line], data


def _task_t52(p0s, lbl):
    p0 = Rational(p0s)
    t5v = Rational(11, 10)
    Sp5 = -(t5v + w1p)
    Spb5 = -(t5v + w1pb)
    R2p5 = simplify(expand(Sp5 * Spb5))
    taup5 = sqrt(1 / R2p5 - Rational(1, 4))
    w2p5 = Sp5 * (Rational(1, 2) + I * taup5)
    w3p5 = Sp5 * (Rational(1, 2) - I * taup5)
    chk(simplify(expand(t5v + w1p + w2p5 + w3p5)) == 0,
        "T5.2 {0}: g == 0 identically on the parametrization".format(lbl))
    R2v = simplify(R2p5.subs(_p, p0))
    chk(R2v > 0 and R2v < 4, "T5.2 {0}: 0 < |S|^2 = {1} < 4".format(lbl, R2v))
    ws5 = [simplify(w.subs(_p, p0)) for w in (w1p, w2p5, w3p5)]
    for w in ws5:
        chk(is_zero_exact(w * conjugate(w) - 1), "T5.2 {0}: |w|=1".format(lbl))
    d1v = simplify(dlog_dir(w1p, _p).subs(_p, p0))
    d2v = simplify(dlog_dir(w2p5, _p).subs(_p, p0))
    d3v = simplify(dlog_dir(w3p5, _p).subs(_p, p0))
    chk(all(is_zero_exact(im(vv)) for vv in (d1v, d2v, d3v)),
        "T5.2 {0}: the tangent is real".format(lbl))
    T5t = Matrix([simplify(re(d1v)), simplify(re(d2v)), simplify(re(d3v))])
    chk(any(not is_zero_exact(T5t[iq]) for iq in range(3)),
        "T5.2 {0}: the tangent is nonzero".format(lbl))
    line, data = t_point_data(ws5, [Integer(1)] * 3, G3, 3,
                              [tuple(T5t)], "T5.2 " + lbl, t=t5v)
    return [line], data


_DISPATCH = {
    "t1": _task_t1, "t2": _task_t2, "t3u": _task_t3u, "t3w": _task_t3w,
    "t3j": _task_t3j, "t4g": _task_t4g, "t4j": _task_t4j,
    "t51": _task_t51, "t52": _task_t52,
}


def run_task(spec):
    drain()   # discard import residue in the worker
    fn = _DISPATCH[spec[0]]
    try:
        lines, data = fn(*spec[1:])
        checks = drain()
    except Exception as ex:
        lines, data = ["  [TASK FAILURE {0}: {1}]".format(spec, ex)], {}
        checks = drain() + [(False, "task {0} crashed: {1}".format(spec, ex))]
    return (spec, checks, lines, data)


TASKS = [
    ("t1", 2), ("t1", 3), ("t1", 4),
    ("t2", 0), ("t2", 1),
    ("t3u",), ("t3w",),
    ("t3j", 1), ("t3j", 2), ("t3j", 3),
    ("t4g", "2", "-2", "P1 (a,b)=(2,-2)"), ("t4g", "2", "1/2", "P2 (a,b)=(2,1/2)"),
    ("t4j", 1), ("t4j", 2), ("t4j", 3), ("t4j", 4),
    ("t51", 0), ("t51", 1),
    ("t52", "1/2", "w1=-1 (p=1/2)"), ("t52", "1/4", "w1=I (p=1/4)"),
]


# ==================== master-side: log Tee ====================

class Tee:
    def __init__(self, real, fh):
        self.real = real
        self.fh = fh
        self.chunks = []

    def write(self, s):
        self.real.write(s)
        self.fh.write(s)
        self.fh.flush()          # flush on every write: progress visible from outside
        self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


def plane_eigs(J, P, G, Ginv, okf):
    JP = J * P
    M = simplify((JP.T) * JP)
    N = simplify(P.T * Ginv * P)
    okf(not is_zero_exact(N.det()), "the plane is non-degenerate (det N != 0)")
    chp = expand((M - _lam * N).det())
    rr = roots(Poly(chp, _lam))
    out = []
    for rt, ml in rr.items():
        out += [simplify(rt)] * ml
    return sorted(out, key=str)


# ============================================================================
#                                   RUN (master)
# ============================================================================

def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S999_run.log"), "w",
                 encoding="utf-8")
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

    for cnd, msg in _IMPORT_CHKS:
        ok(cnd, msg)

    print("=" * 74)
    print("W40: structure of f near the zero-set — linear order and eigenvalues")
    print("     (blind probe; convention: J without 2pi, eigenvalues in (2pi)^2-units;")
    print("      frame A=Cholesky(G), G_ij=<alpha_i,alpha_j>_unit; verdicts exact;")
    print("      Pool(18), collected in fixed order — the log is bit-reproducible)")
    print("=" * 74)

    # ---- parallel accounting: workers return, the master prints ----
    sys.stderr.write("[master] dispatching {0} tasks to Pool(18)\n".format(len(TASKS)))
    sys.stderr.flush()
    RES = {}
    with Pool(processes=18) as pool:
        for spec, checks, lines, data in pool.imap(run_task, TASKS):
            RES[spec] = (checks, lines, data)
            sys.stderr.write("[done] {0}\n".format(spec))
            sys.stderr.flush()

    def emit(spec):
        checks, lines, data = RES[spec]
        for cnd, msg in checks:
            ok(cnd, msg)
        for ln in lines:
            print(ln)
        return data

    # ==================== T1: bit-fence against S998 run.log ====================
    print()
    print("T1 (bit-fence): type/dimension of Z (t=1, d=2,3,4) — cross-check against S998 run.log BIT-FOR-BIT")
    _s998txt = open(os.path.join(_HERE, "S998_run.log"),
                    encoding="utf-8").read().splitlines()
    _int2 = False
    s998rows = {}
    for _ln in _s998txt:
        if _ln.startswith("T2 (m=1)"):
            _int2 = True
            continue
        if _ln.startswith("T2 critical"):
            _int2 = False
            continue
        if not _int2:
            continue
        _parts = [pq.strip() for pq in _ln.split("|")]
        if len(_parts) == 7 and _parts[0] in ("2", "3", "4"):
            s998rows[(_parts[0], _parts[1])] = (_ln.strip(), _parts)
    ok(all((str(dd), "1") in s998rows for dd in (2, 3, 4)),
       "S998 log: T2 t=1 rows for d=2,3,4 parsed")
    for d in (2, 3, 4):
        data = emit(("t1", d))
        logrow = s998rows[(str(d), "1")][0]
        ok(data["row"] == logrow,
           "bit-fence d={0}: row '{1}' == log '{2}'".format(d, data["row"], logrow))

    # ==================== T2 ====================
    print()
    print("T2 (d=2, t=1): expansion g(psi0 + dk) = L(dk) + O(dk^2) at each of the 2 zeros")
    print("point w=(w1,w2) | J_k (rows) | det J_k | eigenvalues J^T J | equal?")
    print("-" * 100)
    emit(("t2", 0))
    emit(("t2", 1))
    print("  (convention: J_k = J_psi*A2, A2 = Cholesky(G2); the O(dk^2) remainder — the series is analytic)")

    # ==================== T3 ====================
    print()
    print("T3 (d=3, t=1): the zero-set = 3 circles (exact certificate) + Jacobian structure")
    _a, _b, _c = symbols('_a _b _c')
    _e1 = _a + _b + _c
    _e2 = _a * _b + _a * _c + _b * _c
    _e3 = _a * _b * _c
    ok(expand((1 + _a) * (1 + _b) * (1 + _c) - (1 + _e1 + _e2 + _e3)) == 0,
       "certificate: (1+a)(1+b)(1+c) = 1 + e1 + e2 + e3 (identically)")
    ok(simplify(_e2 / _e3 - (1 / _a + 1 / _b + 1 / _c)) == 0,
       "certificate: e2/e3 = 1/a + 1/b + 1/c (identically; for |w|=1 this is the conjugate of e1)")
    _s3 = symbols('_s3')
    ok(simplify(1 + Integer(-1) + (-_s3) + _s3) == 0,
       "certificate: 1 + e1 + e2 + e3 = 0 at e1=-1, e2=-e3")
    print("  certificate: on Z(t=1,d=3) some w_i = -1, the other two are antipodes =>")
    print("  Z = C1 (z,-z,-1) ∪ C2 (z,-1,-z) ∪ C3 (-1,z,-z)  (three circles; intersections at")
    print("  the collinear points (1,-1,-1),(-1,1,-1),(-1,-1,1))")
    _z = symbols('_z')
    for lbl, trio in (("C1", (_z, -_z, Integer(-1))), ("C2", (_z, Integer(-1), -_z)),
                      ("C3", (Integer(-1), _z, -_z))):
        ok(expand(1 + trio[0] + trio[1] + trio[2]) == 0,
           "component {0}: g == 0 identically".format(lbl))
    print()
    print("T3(a) GENERIC point on C1 (symbolic parameter u; w1 = ((1-u^2)+2Iu)/(1+u^2)):")
    emit(("t3u",))
    print()
    print("T3 table of points (witness u=1/2 + all symmetric cosets k_j on Z):")
    print("point | rank | ker-vs-tangents | K | eigenvalues (transverse) | equal?")
    print("-" * 110)
    emit(("t3w",))
    # count of ALL symmetric points on Z (master, exact arithmetic)
    t3_onZ = []
    for j in range(0, 4):
        sm = sym_sample(3, j)
        gx, gy = g_val(sm, [Integer(1)] * 3, Integer(1))
        if is_zero_exact(gx) and is_zero_exact(gy):
            t3_onZ.append(j)
    ok(t3_onZ == [1, 2, 3], "T3(b): symmetric points on Z = j=1,2,3 (count of all)")
    for j in t3_onZ:
        emit(("t3j", j))
    # cross-check of the transverse eigenvalues via the pair of forms (M, N) — base for mutant M2
    sm_j1 = sym_sample(3, 1)
    J_j1 = jac_psi(sm_j1, [Integer(1)] * 3)
    K_j1 = simplify(expand(J_j1 * G3 * J_j1.T))
    G3inv = G3.inv()
    P_true = G3 * J_j1.T
    eigs_true = plane_eigs(J_j1, P_true, G3, G3inv, ok)
    l1_j1, l2_j1, _d0, _d1, _d2 = eig2(K_j1)
    ok(sorted([l1_j1, l2_j1], key=str) == eigs_true,
       "cross-check: eigenvalues on the transverse plane (pair of forms) == eig(K) at k_1")

    # ==================== T4 ====================
    print()
    print("T4 (d=4, t=1): a 2-variety — generic points + all symmetric cosets")
    ok(simplify(expand(1 + w1s + w2s + w3s + w4s)) == 0,
       "T4 family: g == 0 identically (w3+w4 = S by construction)")
    print("point | rank | ker-vs-tangents | K | eigenvalues (transverse) | equal?")
    print("-" * 110)
    emit(("t4g", "2", "-2", "P1 (a,b)=(2,-2)"))
    emit(("t4g", "2", "1/2", "P2 (a,b)=(2,1/2)"))
    t4_onZ = []
    for j in range(0, 5):
        sm = sym_sample(4, j)
        gx, gy = g_val(sm, [Integer(1)] * 4, Integer(1))
        if is_zero_exact(gx) and is_zero_exact(gy):
            t4_onZ.append(j)
    ok(t4_onZ == [1, 2, 3, 4], "T4: symmetric points on Z = j=1..4 (count of all)")
    for j in t4_onZ:
        emit(("t4j", j))

    # ==================== T5 ====================
    print()
    print("T5: shift the weight t = 11/10 on axis 0")
    print()
    print("T5.1 (d=2, t=11/10): repeat of T2")
    print("point w=(w1,w2) | J_k (rows) | det J_k | eigenvalues J^T J | equal?")
    print("-" * 100)
    v51a = emit(("t51", 0))
    v51b = emit(("t51", 1))
    print("  verdict T5.1: the eigenvalues at t=11/10 do NOT remain equal"
          if "no" in (v51a.get("verd"), v51b.get("verd"))
          else "  verdict T5.1: remain equal")
    print()
    print("T5.2 (d=3, t=11/10): repeat of T3(a,b) for the surviving component")
    t5v = Rational(11, 10)
    xbound = simplify((3 - t5v ** 2) / (2 * t5v))
    ok(xbound == Rational(179, 220), "T5.2: arc w1: x1 <= 179/220")
    print("  the parametrization is complete: w1 free on the arc x1 <= 179/220 (|t+w1| <= 2),")
    print("  the pair (w2,w3) = S(1/2 +- I tau), S = -(t+w1); the two branches join at tau=0")
    print("  => ONE circle (there are no collinear singularities at t=11/10)")
    print("point | rank | ker-vs-tangents | K | eigenvalues (transverse) | equal?")
    print("-" * 110)
    emit(("t52", "1/2", "w1=-1 (p=1/2)"))
    emit(("t52", "1/4", "w1=I (p=1/4)"))
    print("T5.2(b) symmetric cosets k_j at t=11/10 (the value of g, exact):")
    t52_sym = []
    for j in range(0, 4):
        sm = sym_sample(3, j)
        gx, gy = g_val(sm, [Integer(1)] * 3, t5v)
        g2 = simplify(gx ** 2 + gy ** 2)
        onZ = is_zero_exact(gx) and is_zero_exact(gy)
        if onZ:
            t52_sym.append(j)
        print("  j={0}: |g|^2 = {1}  -> on Z: {2}".format(j, g2, "yes" if onZ else "no"))
    ok(t52_sym == [], "T5.2(b): NO symmetric point is on Z at t=11/10 (count of all)")

    # ==================== MUTANTS (master) ====================
    print()
    print("MUTANTS:")
    mut_ok = True

    # M1: the expansion truncated at O(dk^0): the linear form L is replaced by zero ->
    # the prediction 'g remains 0 near the zero' contradicts the exact c1 != 0
    pts2 = solve_points_d2_m1(Integer(1))
    (px1, py1, px2, py2) = pts2[0]
    W1 = px1 + I * py1
    W2 = px2 + I * py2
    J_m1 = jac_psi([(px1, py1), (px2, py2)], [Integer(1)] * 2)
    det_true = simplify((J_m1 * A2).det())
    c1_dir = simplify(2 * pi * I * (W1 * 1 + W2 * 0))   # direction q=(1,0)
    J_mut = zeros(2, 2)
    if (not is_zero_exact(c1_dir)) and (not is_zero_exact(det_true)) \
            and simplify((J_mut * A2).det()) == 0:
        print("  MUTANT M1: CAUGHT (truncation at O(dk^0): L=0 => det=0, but exactly "
              "c1(q=(1,0)) = {0} != 0, det J_k = {1} != 0)".format(c1_dir, det_true))
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2: confuse tangent/transverse in T3: in the plane with the tangent, the pair
    # of forms (M,N) gives an eigenvalue 0 — contradicts the exact transverse pair at k_1
    P_mut = Matrix.hstack(P_true[:, 0], Matrix([1, 0, 1]))
    eigs_mut = plane_eigs(J_j1, P_mut, G3, G3inv, ok)
    if (Integer(0) in eigs_mut) and (Integer(0) not in eigs_true) \
            and eigs_mut != eigs_true:
        print("  MUTANT M2: CAUGHT (the plane with the tangent: eigenvalues {0} != transverse {1};"
              " the zero gives away the substitution)".format(eigs_mut, eigs_true))
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3: comparing eigenvalues with a float tolerance of 1e-9 on a deliberately degenerate pair
    l_a = Integer(2)
    l_b = Integer(2) + Rational(1, 10 ** 12)
    float_says_equal = abs(float(l_a) - float(l_b)) < 1e-9
    exact_says_equal = is_zero_exact(l_a - l_b)
    if float_says_equal and (not exact_says_equal):
        print("  MUTANT M3: CAUGHT (pair (2, 2+10^-12): float(1e-9) says 'equal', "
              "symbolically the difference = -1/10^12 != 0)")
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4: the T1 bit-fence with a shifted t = 11/10: the fields (type,dim,#points,rank) agree,
    # but the note (solve coordinates) diverges from the log -> the fence catches it
    v_m4 = classify_m1(2, Rational(11, 10))
    for cnd, msg in drain():
        ok(cnd, msg)
    logrow_fields = s998rows[("2", "1")][1]
    m4_match_coarse = (v_m4["typ"] == logrow_fields[2]
                       and str(v_m4["dim"]) == logrow_fields[3]
                       and str(v_m4["npts"]) == logrow_fields[4]
                       and str(v_m4["rank"]) == logrow_fields[5])
    m4_note_differs = (v_m4["note"] != logrow_fields[6])
    if m4_match_coarse and m4_note_differs:
        print("  MUTANT M4: CAUGHT (t=11/10: type/dim/#/rank would agree, but the note "
              "'{0}' != log '{1}')".format(v_m4["note"], logrow_fields[6]))
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # M5: a planted 'symmetric point' j=0 is NOT on Z -> the point-on-Z assert catches it
    sm0 = sym_sample(3, 0)
    gx0, gy0 = g_val(sm0, [Integer(1)] * 3, Integer(1))
    if not (is_zero_exact(gx0) and is_zero_exact(gy0)):
        print("  MUTANT M5: CAUGHT (j=0: g = {0} != 0 — the point is not on Z, "
              "the 'point on Z' assert fails the probe)".format(simplify(gx0 + I * gy0)))
    else:
        print("  MUTANT M5: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): weight 13/10 on two axes, d=3")
    random.seed(999331)
    axes_nc = sorted(random.sample(range(4), 2))
    t_nc = Rational(13, 10)
    c_full = [Integer(1)] * 4
    for aq in axes_nc:
        c_full[aq] = t_nc
    base_cnt = []
    nc_cnt = []
    for j in range(4):
        fb = Add(*[exp(2 * pi * I * Rational((iq * j) % 4, 4)) for iq in range(4)])
        fn = Add(*[c_full[iq] * exp(2 * pi * I * Rational((iq * j) % 4, 4))
                   for iq in range(4)])
        if simplify(expand(fb * conjugate(fb))) == 0:
            base_cnt.append(j)
        if simplify(expand(fn * conjugate(fn))) == 0:
            nc_cnt.append(j)
    ok(base_cnt == [1, 2, 3], "control: the base t=1 has 3 symmetric zeros")
    ok(len(nc_cnt) != len(base_cnt),
       "control: the count of symmetric zeros is DIFFERENT ({0} != {1})".format(
           len(nc_cnt), len(base_cnt)))
    print("  axes={0}, weight=13/10: symmetric zeros j={1} (count {2}) against the base j={3} (count 3)"
          .format(axes_nc, nc_cnt, len(nc_cnt), base_cnt))
    for j in nc_cnt:
        sm = sym_sample(3, j)
        c_red = c_full[1:]
        gx, gy = g_val(sm, c_red, c_full[0])
        ok(is_zero_exact(gx) and is_zero_exact(gy),
           "control: surviving zero j={0}".format(j))
        Jnc = jac_psi(sm, c_red)
        rnc = Jnc.rank(iszerofunc=_izf)
        Knc = simplify(expand(Jnc * G3 * Jnc.T))
        l1n, l2n, discn, _t1, _t2 = eig2(Knc)
        print("  surviving j={0}: rank={1}, K={2}, eigenvalues=({3}, {4}), equal: {5}"
              "  <- the structure differs from the base (there rank=2, 8/3=8/3)"
              .format(j, rnc, kfmt(Knc), l1n, l2n, eq_verdict(discn)))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

    # ==================== FORBIDDEN-SCAN (tools.fence_scan) ====================
    _pp = [("ч", "ас"), ("t", "ime"), ("тем", "пор"), ("кау", "зал"),
           ("сві", "тло"), ("швид", "к"), ("Мінков", "ськ"), ("кон", "ус"),
           ("енер", "г")]
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


if __name__ == "__main__":
    main()
