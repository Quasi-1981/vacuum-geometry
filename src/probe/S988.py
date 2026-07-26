# -*- coding: utf-8 -*-
# DIM: na (W39 tower step: centralizers of tower-nilpotents in so(4,4); 0 handles).
#
# ============================================================================
# CARVE CHOICE (stamped BEFORE any counting)
# ----------------------------------------------------------------------------
# so(4,4), eta = diag(+^4, -^4).  Three objects, one shared pipeline (the
# S933/S937/S940/S950 primitives, VERBATIM), each measured raw:
#   (A) TWO-WEDGE  N = N1+N2, two wedges over pairwise-orthogonal totally
#       isotropic supports (rank 4, N^2=0);
#   (B) DEEP       one nilpotent with N^2!=0 (Jordan-3 block inside an
#       so(2,1) on a (2,1)-signature 3-space; built explicitly);
#   (C) CONTROL    one single wedge (rank 2, N^2=0).
# For each: full centralizer c(N) dim/basis, Levi+radical split, Levi type
# (Killing exact), radical lower-central-series + number of central
# directions, Levi-invariant bilinear forms on the module (symmetric AND
# antisymmetric, exact dims), raw [module,module] bracket table into the
# centre, and the residual so(eta|G) signature via eigenvalues of squares
# of the Levi generators.  Branches CAUGHT by mutants; negative control;
# bit-fence vs the S950 single-wedge (5,2)r0 row; FORBIDDEN-SCAN; STOP.
# No targets are assumed -- every column is a raw measurement.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
from sympy import Matrix, Integer, zeros, eye, diag, Rational

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S988_run.log")
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


# ==================== primitives (VERBATIM from S950/S940/S937/S933) ====================

def make_eta(p, q):
    return diag(*([Integer(1)] * p + [Integer(-1)] * q))


def unit(n, i):
    v = zeros(n, 1)
    v[i, 0] = Integer(1)
    return v


def wedge(x, y, eta):
    return x * (eta * y).T - y * (eta * x).T


def so_basis(n, eta):
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            out.append(wedge(unit(n, i), unit(n, j), eta))
    return out


def flat(M):
    n = M.rows
    return Matrix(n * n, 1, list(M))


def stack_flats(mats, sq):
    if not mats:
        return zeros(sq * sq, 0)
    return Matrix.hstack(*[flat(M) for M in mats])


def span_basis(mats, sq):
    out = []
    F = zeros(sq * sq, 0)
    r = 0
    for M in mats:
        F2 = Matrix.hstack(F, flat(M))
        if F2.rank() > r:
            out.append(M)
            F = F2
            r += 1
    return out


def centralizer(A, bas):
    n = A.rows
    cols = [flat(B * A - A * B) for B in bas]
    ns = Matrix.hstack(*cols).nullspace()
    cb = []
    for v in ns:
        M = zeros(n, n)
        for k in range(len(bas)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * bas[k]
        cb.append(M)
    return cb


def is_nilp(A):
    n = A.rows
    A2 = A * A
    if A2.is_zero_matrix:
        return True
    if A2.trace() != 0:
        return False
    return (A ** n).is_zero_matrix


def nilp_depth(A):
    n = A.rows
    P = eye(n)
    for k in range(1, n + 1):
        P = P * A
        if P.is_zero_matrix:
            return k
    return None


def in_so(M, eta):
    return (M.T * eta + eta * M).is_zero_matrix


def bracket_basis(L, n):
    brs = []
    for a in range(len(L)):
        for b in range(a + 1, len(L)):
            brs.append(L[a] * L[b] - L[b] * L[a])
    return span_basis(brs, n)


def derived_series(L, n):
    cur = span_basis(L, n)
    dims = [len(cur)]
    while len(cur) > 0:
        nxt = bracket_basis(cur, n)
        if len(nxt) == len(cur):
            break
        dims.append(len(nxt))
        cur = nxt
    return dims


def perfect_core(L, n):
    cur = span_basis(L, n)
    while len(cur) > 0:
        nxt = bracket_basis(cur, n)
        if len(nxt) == len(cur):
            return cur
        cur = nxt
    return []


def lower_central_series(L, n):
    L0 = span_basis(L, n)
    cur = L0
    dims = [len(cur)]
    while len(cur) > 0:
        brs = []
        for a in range(len(L0)):
            for b in range(len(cur)):
                brs.append(L0[a] * cur[b] - cur[b] * L0[a])
        nxt = span_basis(brs, n)
        if len(nxt) == len(cur):
            break
        dims.append(len(nxt))
        cur = nxt
    return dims


def coords_in(basis, M, sq):
    F = stack_flats(basis, sq)
    v = flat(M)
    if F.cols == 0:
        return None if not v.is_zero_matrix else Matrix(0, 1, [])
    try:
        sol, params = F.gauss_jordan_solve(v)
    except ValueError:
        return None
    if params.rows * params.cols > 0:
        sol = sol.subs({s: 0 for s in params})
    if not (F * sol - v).is_zero_matrix:
        return None
    return sol


def cong_signature(Gin):
    G = Gin.copy()
    pos = neg = zer = 0
    while G.rows > 0:
        m = G.rows
        k = None
        for i in range(m):
            if G[i, i] != 0:
                k = i
                break
        if k is None:
            pair = None
            for i in range(m):
                for j in range(i + 1, m):
                    if G[i, j] != 0:
                        pair = (i, j)
                        break
                if pair is not None:
                    break
            if pair is None:
                zer += m
                break
            i, j = pair
            G[i, :] = G[i, :] + G[j, :]
            G[:, i] = G[:, i] + G[:, j]
            continue
        dv = G[k, k]
        if dv > 0:
            pos += 1
        else:
            neg += 1
        v = G[:, k]
        G = G - (v * v.T) / dv
        G.row_del(k)
        G.col_del(k)
    return (pos, neg, zer)


def ad_matrices(L, n):
    k = len(L)
    ads = []
    for i in range(k):
        cols = []
        for j in range(k):
            c = coords_in(L, L[i] * L[j] - L[j] * L[i], n)
            if c is None:
                return None
            cols.append(c)
        ads.append(Matrix.hstack(*cols))
    return ads


def killing_matrix(L, n):
    ads = ad_matrices(L, n)
    if ads is None:
        return None
    k = len(L)
    K = zeros(k, k)
    for i in range(k):
        for j in range(k):
            K[i, j] = (ads[i] * ads[j]).trace()
    return K


def killing_sig(L, n):
    K = killing_matrix(L, n)
    return None if K is None else cong_signature(K)


def n_invforms(L, n):
    ads = ad_matrices(L, n)
    if ads is None:
        return None
    k = len(L)
    symb = []
    for i in range(k):
        for j in range(i, k):
            E = zeros(k, k)
            E[i, j] = Integer(1)
            E[j, i] = Integer(1)
            symb.append(E)
    cols = []
    for E in symb:
        cols.append(Matrix.vstack(*[flat(a.T * E + E * a) for a in ads]))
    ns = Matrix.hstack(*cols).nullspace()
    return len(ns)


def fingerprint(L, n):
    b = span_basis(L, n)
    return (len(b), killing_sig(b, n), n_invforms(b, n))


def killing_radical(core, n):
    K = killing_matrix(core, n)
    ns = K.nullspace()
    rad = []
    for v in ns:
        M = zeros(n, n)
        for t in range(len(core)):
            if v[t, 0] != 0:
                M = M + v[t, 0] * core[t]
        rad.append(M)
    return span_basis(rad, n)


def center_of(cb, n):
    if not cb:
        return []
    cols = []
    for k in range(len(cb)):
        cols.append(Matrix.vstack(*[flat(cb[k] * cb[j] - cb[j] * cb[k]) for j in range(len(cb))]))
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        M = zeros(n, n)
        for k in range(len(cb)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * cb[k]
        out.append(M)
    return span_basis(out, n)


def complement_reps(core, rad, n):
    F = stack_flats(rad, n)
    rk = F.rank()
    reps = []
    for M in core:
        F2 = Matrix.hstack(F, flat(M))
        if F2.rank() > rk:
            reps.append(M)
            F = F2
            rk += 1
    return reps


def end_basis(d):
    out = []
    for i in range(d):
        for j in range(d):
            E = zeros(d, d)
            E[i, j] = Integer(1)
            out.append(E)
    return out


def commutant_dim(ops, dv):
    eb = end_basis(dv)
    cols = []
    for E in eb:
        cols.append(Matrix.vstack(*[flat(E * op - op * E) for op in ops]))
    return dv * dv - Matrix.hstack(*cols).rank()


def cyclic_dim(v0, ops, dv):
    F = v0.copy()
    r = F.rank()
    frontier = [v0]
    while frontier:
        nf = []
        for w in frontier:
            for op in ops:
                u = op * w
                F2 = Matrix.hstack(F, u)
                if F2.rank() > r:
                    F = F2
                    r += 1
                    nf.append(u)
        frontier = nf
    return r


def min_cyclic_dim(ops, dv):
    best = dv
    for i in range(dv):
        d = cyclic_dim(unit(dv, i), ops, dv)
        if 0 < d < best:
            best = d
    return best


def weight_multiset(ops, dv):
    best = None
    best_score = -1
    for op in ops:
        ev = op.eigenvals()
        if sum(ev.values()) == dv and all(v.is_rational for v in ev.keys()):
            nonzero = any(k != 0 for k in ev.keys())
            score = len(ev) + (100 if nonzero else 0)
            if score > best_score:
                best_score = score
                best = {int(k) if k.is_integer else str(k): int(m) for k, m in ev.items()}
    return best


def noncentral(rad, Z, n):
    F = stack_flats(Z, n)
    rk = F.rank()
    out = []
    for M in rad:
        F2 = Matrix.hstack(F, flat(M))
        if F2.rank() > rk:
            out.append(M)
            F = F2
            rk += 1
    return out


def classify_radical(rad, n):
    dr = len(rad)
    if dr == 0:
        return "0", dict(dr=0, dz=0, drr=0)
    Z = center_of(rad, n)
    dz = len(Z)
    drr = len(bracket_basis(rad, n))
    nc = noncentral(rad, Z, n)
    allmove = all(len(span_basis([x * r - r * x for r in rad], n)) > 0 for x in nc)
    info = dict(dr=dr, dz=dz, drr=drr, allmove=allmove)
    if drr == 0:
        return "R%d" % dr, info
    if dz == 1 and drr == 1 and allmove:
        return "h%d" % dr, info
    a = dz - drr
    b = dr - a
    info["a"] = a
    info["b"] = b
    return "R%d+h%d" % (a, b), info


def model_from_constraints(d, cons):
    eb = []
    for i in range(d):
        for j in range(d):
            E = zeros(d, d)
            E[i, j] = Integer(1)
            eb.append(E)
    cols = []
    for E in eb:
        parts = [flat(fn(E)) for fn in cons]
        cols.append(Matrix.vstack(*parts))
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        M = zeros(d, d)
        for k in range(d * d):
            if v[k, 0] != 0:
                M = M + v[k, 0] * eb[k]
        out.append(M)
    return out


def sp_form(m):
    J = zeros(2 * m, 2 * m)
    for i in range(m):
        J[i, m + i] = Integer(1)
        J[m + i, i] = Integer(-1)
    return J


def model_sp_w(w):
    d = w.rows
    return model_from_constraints(d, [lambda E: E.T * w + w * E])


# ---------- Levi identification menu (Killing signature exact) ----------
LEVI_SIMPLE = [("so(3)", (3, (0, 3, 0))), ("so(2,1)", (3, (2, 1, 0))),
               ("so(4)", (6, (0, 6, 0))), ("so(3,1)", (6, (3, 3, 0))),
               ("so(2,2)", (6, (4, 2, 0))), ("so(5)", (10, (0, 10, 0))),
               ("so(4,1)", (10, (4, 6, 0))), ("so(3,2)", (10, (6, 4, 0))),
               ("sp(4,R)", (10, (6, 4, 0)))]
LEVI_MENU = list(LEVI_SIMPLE)
for _i in range(len(LEVI_SIMPLE)):
    for _j in range(_i, len(LEVI_SIMPLE)):
        _ni, _fi = LEVI_SIMPLE[_i]
        _nj, _fj = LEVI_SIMPLE[_j]
        LEVI_MENU.append((_ni + "+" + _nj,
                          (_fi[0] + _fj[0],
                           (_fi[1][0] + _fj[1][0], _fi[1][1] + _fj[1][1], _fi[1][2] + _fj[1][2]))))


def quotient_ads(reps, rad, n):
    full = reps + rad
    kq = len(reps)
    ads = []
    for i in range(kq):
        cols = []
        for j in range(kq):
            co = coords_in(full, reps[i] * reps[j] - reps[j] * reps[i], n)
            cols.append(Matrix(kq, 1, [co[t, 0] for t in range(kq)]))
        ads.append(Matrix.hstack(*cols))
    return ads


def identify_levi_quotient(reps, rad, n):
    if not reps:
        return "Levi=0", (0, (0, 0, 0))
    ads = quotient_ads(reps, rad, n)
    kq = len(reps)
    K = zeros(kq, kq)
    for i in range(kq):
        for j in range(kq):
            K[i, j] = (ads[i] * ads[j]).trace()
    ks = cong_signature(K)
    fp = (kq, ks)
    names = [nm for (nm, f) in LEVI_MENU if f == fp]
    return ("Levi=" + "&".join(names)) if names else "Levi=OTHER", fp


# ---------- module data: Levi action on rad (or rad/Z), plus centre ----------
def module_data(reps, rad, n):
    Z = span_basis(center_of(rad, n), n)
    drr = len(bracket_basis(rad, n))
    if drr == 0:
        Mbasis = list(rad)
        full = list(rad)
        dz = 0
    else:
        Mbasis = noncentral(rad, Z, n)
        full = Z + Mbasis
        dz = len(Z)
    dv = len(Mbasis)
    ops = []
    for z in reps:
        cols = []
        for b in range(dv):
            c = coords_in(full, z * Mbasis[b] - Mbasis[b] * z, n)
            cols.append(Matrix(dv, 1, [c[dz + t, 0] for t in range(dv)]))
        ops.append(Matrix.hstack(*cols) if dv > 0 else zeros(0, 0))
    return dict(Z=Z, Mbasis=Mbasis, dv=dv, ops=ops, dz=dz)


def invariant_forms_split(ops, dv):
    """Levi-invariant bilinear forms B on the dv-module: op^T B + B op = 0.
    Return (n_symmetric, n_antisymmetric) exact dims."""
    if dv == 0:
        return 0, 0
    eb = end_basis(dv)

    def count(symtype):
        cols = []
        for E in eb:
            parts = [flat(op.T * E + E * op) for op in ops]
            if symtype == 'sym':
                parts.append(flat(E - E.T))
            else:
                parts.append(flat(E + E.T))
            cols.append(Matrix.vstack(*parts))
        return len(Matrix.hstack(*cols).nullspace())

    return count('sym'), count('anti')


def module_bracket_table(Mbasis, Zbasis, n):
    """Raw [m_i, m_j] expressed in centre-basis coords; None => outside centre span."""
    rows = []
    for i in range(len(Mbasis)):
        for j in range(i + 1, len(Mbasis)):
            br = Mbasis[i] * Mbasis[j] - Mbasis[j] * Mbasis[i]
            if br.is_zero_matrix:
                rows.append((i, j, "0"))
                continue
            c = coords_in(Zbasis, br, n)
            if c is None:
                rows.append((i, j, "OUTSIDE-CENTRE"))
            else:
                rows.append((i, j, [c[t, 0] for t in range(c.rows)]))
    return rows


def gen_square_eigs(basis):
    """Eigenvalues of the square of each generator on the ambient rep (raw)."""
    out = []
    for M in basis:
        ev = (M * M).eigenvals()
        out.append({str(k): int(v) for k, v in ev.items()})
    return out


# ---------- wedge constructors ----------
def d0_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 2 and len(Km) >= 2:
        return unit(n, Kp[0]) + unit(n, Km[0]), unit(n, Kp[1]) + unit(n, Km[1])
    return None


def build_h(m2plus1):
    m = (m2plus1 - 1) // 2
    d = m + 2
    gens = []
    z = zeros(d, d); z[0, d - 1] = Integer(1)
    gens.append(z)
    for i in range(1, m + 1):
        X = zeros(d, d); X[0, i] = Integer(1); gens.append(X)
    for i in range(1, m + 1):
        Y = zeros(d, d); Y[i, d - 1] = Integer(1); gens.append(Y)
    return gens, d


# ==================== unified analyzer ====================
def analyze_object(tag, N, n, eta):
    print("=== OBJECT: {0} ===".format(tag))
    bas = so_basis(n, eta)
    ok(in_so(N, eta), tag + ": N in so(eta)")
    dep = nilp_depth(N)
    ok(dep is not None, tag + ": N nilpotent")
    print("  N: rank={0} nilp_depth(min k: N^k=0)={1}".format(N.rank(), dep))
    c = centralizer(N, bas)
    # verify every centralizer basis element commutes exactly
    ok(all((B * N - N * B).is_zero_matrix for B in c), tag + ": centralizer commutes exactly")
    dc = len(c)
    dcc = len(bracket_basis(c, n))
    dser = derived_series(c, n)
    Zc = center_of(c, n)
    print("  c(N): dim={0} dim[c,c]={1} derived={2} dim_centre(c)={3}".format(dc, dcc, dser, len(Zc)))
    core = perfect_core(c, n)
    print("  perfect_core dim={0}".format(len(core)))
    if not core:
        lcs = lower_central_series(c, n)
        Zr = center_of(c, n)
        rtag, rinfo = classify_radical(span_basis(c, n), n)
        print("  class=SOLVABLE  radical=whole c  type={0} lcs={1} dim_centre={2} info={3}".format(
            rtag, lcs, len(Zr), rinfo))
        print("  -- STOP (raw, no interpretation) --")
        return dict(tag=tag, dimc=dc, cls="SOLVABLE")
    cid, fp = identify_core(core, n) if True else (None, None)
    rad = killing_radical(core, n)
    reps = complement_reps(core, rad, n)
    ltag, lfp = identify_levi_quotient(reps, rad, n)
    print("  core: dim={0} Killing={1}".format(len(core), fp[1]))
    print("  Levi: {0}  dim={1} Killing_sig={2}".format(ltag, lfp[0], lfp[1]))
    # residual so(eta|G) signature via eigenvalues of squares of Levi generators (raw)
    sq = gen_square_eigs(reps)
    print("  Levi gen^2 eigenvalues on {0}-dim rep (raw, per generator):".format(n))
    for i, e in enumerate(sq):
        print("    g{0}^2 eig={1}".format(i, e))
    # radical structure
    rtag, rinfo = classify_radical(rad, n)
    Zr = center_of(rad, n)
    lcs = lower_central_series(rad, n)
    print("  radical: type={0} dim={1} lower_central_series={2} n_central_dirs={3} info={4}".format(
        rtag, len(rad), lcs, len(Zr), rinfo))
    # module
    md = module_data(reps, rad, n)
    dv = md["dv"]
    if dv > 0:
        cdim = commutant_dim(md["ops"], dv)
        wm = weight_multiset(md["ops"], dv)
        mc = min_cyclic_dim(md["ops"], dv) if cdim > 1 else dv
        nsym, nanti = invariant_forms_split(md["ops"], dv)
        print("  module: dv={0} commutant={1} weights={2} mincyc={3}".format(dv, cdim, wm, mc))
        print("  module Levi-invariant bilinear forms: symmetric={0} antisymmetric={1}".format(nsym, nanti))
        Zb = md["Z"] if md["Z"] else Zr
        print("  [module,module] bracket table (coords in centre basis; else flag):")
        for (i, j, val) in module_bracket_table(md["Mbasis"], Zb, n):
            print("    [m{0},m{1}] -> {2}".format(i, j, val))
    else:
        print("  module: dv=0 (radical is its own centre / no translation module)")
    print("  -- STOP (raw, no interpretation) --")
    return dict(tag=tag, dimc=dc, cls="PERFECT-NON-SS" if fp[1][2] > 0 else "SEMISIMPLE-CORE",
                levi=ltag, radical=rtag)


def identify_core(core, n):
    fp = fingerprint(core, n)
    return "dim{0}".format(fp[0]), fp


# ==================== bit-fence: shared machinery vs S950 (5,2)r0 ====================
print("--- bit-fence: shared pipeline reproduces S950 single-wedge (5,2) rank0 (n=7) ---")
_eta7 = make_eta(5, 2)
_bas7 = so_basis(7, _eta7)
_xy = d0_on(list(range(7)), 5, 7)
_Nf = wedge(_xy[0], _xy[1], _eta7)
_cf = centralizer(_Nf, _bas7)
ok(len(_cf) == 13, "fence dimc (5,2)r0 == 13 (got {0})".format(len(_cf)))
_coref = perfect_core(_cf, 7)
_radf = killing_radical(_coref, 7)
_repsf = complement_reps(_coref, _radf, 7)
_ltf, _ = identify_levi_quotient(_repsf, _radf, 7)
ok("so(3)" in _ltf and "so(2,1)" in _ltf, "fence Levi (5,2)r0 = so(3)+so(2,1) (got {0})".format(_ltf))
_rtf, _ = classify_radical(_radf, 7)
ok(_rtf == "h7", "fence radical (5,2)r0 == h7 (got {0})".format(_rtf))
print("  FENCE (5,2)r0: dimc={0} {1} radical={2}".format(len(_cf), _ltf, _rtf))


# ==================== so(4,4) objects ====================
N8 = 8
ETA8 = make_eta(4, 4)

# (A) TWO-WEDGE: pairwise orthogonal totally isotropic supports
x1 = unit(8, 0) + unit(8, 4)
y1 = unit(8, 1) + unit(8, 5)
x2 = unit(8, 2) + unit(8, 6)
y2 = unit(8, 3) + unit(8, 7)
_supports = [x1, y1, x2, y2]
for _a in range(4):
    for _b in range(_a, 4):
        val = (_supports[_a].T * ETA8 * _supports[_b])[0, 0]
        ok(val == 0, "two-wedge supports totally isotropic/orthogonal pair ({0},{1})".format(_a, _b))
N1 = wedge(x1, y1, ETA8)
N2 = wedge(x2, y2, ETA8)
NTWO = N1 + N2
ok((NTWO * NTWO).is_zero_matrix, "two-wedge N^2 == 0")
ok(NTWO.rank() == 4, "two-wedge rank == 4 (got {0})".format(NTWO.rank()))

# (C) CONTROL: single wedge (rank 2)
NCTRL = wedge(x1, y1, ETA8)
ok((NCTRL * NCTRL).is_zero_matrix, "control single-wedge N^2 == 0")
ok(NCTRL.rank() == 2, "control rank == 2 (got {0})".format(NCTRL.rank()))

# (B) DEEP: Jordan-3 nilpotent inside so(2,1) on the (2,1) 3-space {e0,e1,e4}
NDEEP = wedge(unit(8, 1), unit(8, 0) + unit(8, 4), ETA8)
ok(not (NDEEP * NDEEP).is_zero_matrix, "deep N^2 != 0")
ok((NDEEP ** 3).is_zero_matrix, "deep N^3 == 0")
ok(NDEEP.rank() == 2, "deep rank == 2 (Jordan [3,1^5])")

print()
res_ctrl = analyze_object("(C) CONTROL single-wedge so(4,4) rank2", NCTRL, N8, ETA8)
print()
res_two = analyze_object("(A) TWO-WEDGE so(4,4) rank4 (N^2=0)", NTWO, N8, ETA8)
print()
res_deep = analyze_object("(B) DEEP so(4,4) (N^2!=0, Jordan-3)", NDEEP, N8, ETA8)


# ==================== mutants — one per classification branch ====================
print()
print("--- mutants (one per branch) ---")
mut_ok = True

# br-depth: N^2=0 (two-wedge/control) vs N^2!=0 (deep)
if (NTWO * NTWO).is_zero_matrix and not (NDEEP * NDEEP).is_zero_matrix:
    print("MUTANT br-depth: CAUGHT (N^2=0 for wedges; N^2!=0 for deep)")
else:
    print("MUTANT br-depth: NOT CAUGHT"); mut_ok = False

# br-rank: rank 4 (two-wedge) vs rank 2 (control/deep)
if NTWO.rank() == 4 and NCTRL.rank() == 2 and NDEEP.rank() == 2:
    print("MUTANT br-rank: CAUGHT (two-wedge rank4; single/deep rank2)")
else:
    print("MUTANT br-rank: NOT CAUGHT"); mut_ok = False

# br-constructible: definite (8,0) has no isotropic wedge; (4,4) does
_eta80 = make_eta(8, 0)
_iso_exists_80 = any(((unit(8, i) + unit(8, j)).T * _eta80 * (unit(8, i) + unit(8, j)))[0, 0] == 0
                     for i in range(8) for j in range(8) if i != j)
_iso_exists_44 = ((unit(8, 0) + unit(8, 4)).T * ETA8 * (unit(8, 0) + unit(8, 4)))[0, 0] == 0
if (not _iso_exists_80) and _iso_exists_44:
    print("MUTANT br-constructible: CAUGHT ((8,0) no isotropic e_i+e_j; (4,4) has)")
else:
    print("MUTANT br-constructible: NOT CAUGHT"); mut_ok = False

# br-radical: h vs abelian vs mixed (canonical builds)
_h3, _dh3 = build_h(3)
_h5, _dh5 = build_h(5)
_R3 = [unit(4, 0) * unit(4, 1).T, unit(4, 0) * unit(4, 2).T, unit(4, 0) * unit(4, 3).T]
_th3, _ = classify_radical(_h3, _dh3)
_th5, _ = classify_radical(_h5, _dh5)
_tR3, _ = classify_radical(_R3, 4)
if _th3 == "h3" and _th5 == "h5" and _tR3 == "R3":
    print("MUTANT br-radical: CAUGHT (h3->h3, h5->h5, R3->R3)")
else:
    print("MUTANT br-radical: NOT CAUGHT ({0},{1},{2})".format(_th3, _th5, _tR3)); mut_ok = False

# br-module irreducible vs reducible (commutant dim)
_hh = diag(3, 1, -1, -3)
_ee = zeros(4, 4); _ee[0, 1] = Integer(1); _ee[1, 2] = Integer(2); _ee[2, 3] = Integer(3)
_ff = zeros(4, 4); _ff[1, 0] = Integer(3); _ff[2, 1] = Integer(2); _ff[3, 2] = Integer(1)
_ci = commutant_dim([_hh, _ee, _ff], 4)
_h2 = diag(1, -1, 1, -1)
_e2 = zeros(4, 4); _e2[0, 1] = Integer(1); _e2[2, 3] = Integer(1)
_f2 = zeros(4, 4); _f2[1, 0] = Integer(1); _f2[3, 2] = Integer(1)
_cr = commutant_dim([_h2, _e2, _f2], 4)
if _ci == 1 and _cr == 4:
    print("MUTANT br-module: CAUGHT (irreducible commutant=1; reducible 2+2 commutant=4)")
else:
    print("MUTANT br-module: NOT CAUGHT ({0},{1})".format(_ci, _cr)); mut_ok = False

# br-invform sym vs antisym: so(3) vector rep -> (sym=1, anti=0); sp(2,R) -> (sym=0, anti=1)
_so3 = so_basis(3, make_eta(3, 0))
_sym3, _anti3 = invariant_forms_split(_so3, 3)
_sp2 = model_sp_w(sp_form(1))
_sym2, _anti2 = invariant_forms_split(_sp2, 2)
if (_sym3, _anti3) == (1, 0) and (_sym2, _anti2) == (0, 1):
    print("MUTANT br-invform: CAUGHT (so(3):sym1/anti0 ; sp(2,R):sym0/anti1)")
else:
    print("MUTANT br-invform: NOT CAUGHT (so3={0} sp2={1})".format((_sym3, _anti3), (_sym2, _anti2))); mut_ok = False


# ==================== negative control ====================
# a specific so(eta) element that does NOT commute with the two-wedge must be
# excluded from the centralizer span (centralizer is a proper subalgebra).
_g = wedge(unit(8, 0), unit(8, 1), ETA8)  # compact rotation in +,+ plane
_noncomm = not (_g * NTWO - NTWO * _g).is_zero_matrix
_ctwo = centralizer(NTWO, so_basis(8, ETA8))
_excluded = coords_in(_ctwo, _g, 8) is None
ok(_noncomm and _excluded, "neg-control: non-commuting so element excluded from centralizer")

# sign-fence: eta -> -eta gives same two-wedge core fingerprint (isomorphic realization)
_cA = perfect_core(centralizer(NTWO, so_basis(8, ETA8)), 8)
_cB = perfect_core(centralizer(NTWO, so_basis(8, -ETA8)), 8)
ok(fingerprint(_cA, 8) == fingerprint(_cB, 8), "sign-fence: two-wedge core fingerprint eta->-eta")


# ==================== summary ====================
print()
print("SUMMARY: objects=3 (control/two-wedge/deep) + bit-fence + 6 mutants + neg-control")
print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

_pieces = [("ен", "ерг"), ("ч", "ас"), ("кв", "ант"), ("ма", "са")]
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
sys.stdout = _tee.real
_logf.close()
sys.exit(_exit)
