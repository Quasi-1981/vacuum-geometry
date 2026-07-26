# -*- coding: utf-8 -*-
# DIM: na (W39 leg2: two-wedge centralizers in so(5,4) and so(5,5); 0 handles).
#
# ============================================================================
# CARVE CHOICE (stamped BEFORE any counting)
# ----------------------------------------------------------------------------
# Same shared pipeline as S988 (J-0444-verified), VERBATIM, on four objects,
# each measured raw:
#   (C1) single wedge in so(5,4)  rank2
#   (A1) two-wedge  in so(5,4)  rank4 (pairwise-orthogonal isotropic supports)
#   (C2) single wedge in so(5,5)  rank2
#   (A2) two-wedge  in so(5,5)  rank4  <- main object
# Per object, the S988 raw package PLUS two extensions:
#   [EXT-1] invariant bilinear forms PER Levi-irreducible module block (sym/anti
#           dims) and per ordered pair of distinct blocks (cross-invariants);
#   [EXT-2] the [module,module] bracket table in the centre basis, BY BLOCKS:
#           (block_i, block_i) and (block_i, block_j), raw elementwise.
# Module block operationalization (CARVED): commutant scalar -> irreducible
# (one block = whole module); else split by the rational eigenspaces of the
# first commutant element with >=2 rational eigenvalues (each eigenspace is
# Levi-invariant since it commutes); no rational split -> NAMED-STOP (no fake).
# Branches CAUGHT by mutants; negative control; bit-fence vs S988 (4,4) control
# (dimc=18 / h9); FORBIDDEN-SCAN; raw tables; STOP.  No targets assumed.
# Workbench limit (so(5,5) dim 45): if exact arithmetic stalls -> named report.
# ============================================================================
#
# ============================================================================
# PROVENANCE / EPISTEMIC CLASS  (self-report, carved per Omega's request 2026-07-18)
# ----------------------------------------------------------------------------
# CLASS = DECLARED-CONTAMINATED CORROBORATION (NOT a blind witness).
#   The "blind" token in this filename is INACCURATE for this run (kept for the
#   committed reference; see below).  Load-bearing blind witness for W39-leg2 =
#   S991 (fresh zero-context subagent, Omega route-(A), commit 86d6efe3).
# Timeline (UTC 2026-07-18):
#   08:55:48  Omega "gate resolved -- do it" -> lane-A (route decision).
#   09:11:53  Alpha (321d3b67) DECLARED contamination to Omega (fork xmsg).
#   09:13:01  Alpha reserved S990 and ran this probe  <-- count STARTED AFTER the
#             contamination declaration (~68 s later): declared-before-count.
#   09:13:27  fresh hand reserved S991 (parallel, 26 s later).
# Knowledge at count-time: T23 (saturation-law) + W39-leg1 ex-ante = YES;
#   W39_LEG2_* ex-ante / leg-2 stakes = NO (leg-2 blindness held).
# Outcome: conclusions IDENTICAL to blind S991 (T23 refuted by empty residual;
#   (5,5) Lambda^2 cocycle) -- reported AGAINST Alpha's own T23 expectation, so
#   the contamination did not bias the sign.  Accepted by Omega as corroboration
#   (visa J-0446); S991 remains the load-bearing witness.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
from sympy import Matrix, Integer, zeros, eye, diag

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S990_run.log")
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


# ==================== primitives (VERBATIM from S988/S950) ====================

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
    if dv == 0:
        return 0
    eb = end_basis(dv)
    cols = []
    for E in eb:
        cols.append(Matrix.vstack(*[flat(E * op - op * E) for op in ops]))
    return dv * dv - Matrix.hstack(*cols).rank()


def commutant_basis(ops, dv):
    if dv == 0:
        return []
    eb = end_basis(dv)
    cols = [Matrix.vstack(*[flat(E * op - op * E) for op in ops]) for E in eb]
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        X = zeros(dv, dv)
        for k in range(dv * dv):
            if v[k, 0] != 0:
                X = X + v[k, 0] * eb[k]
        out.append(X)
    return out


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


LEVI_SIMPLE = [("so(3)", (3, (0, 3, 0))), ("so(2,1)", (3, (2, 1, 0))),
               ("so(4)", (6, (0, 6, 0))), ("so(3,1)", (6, (3, 3, 0))),
               ("so(2,2)", (6, (4, 2, 0))), ("so(5)", (10, (0, 10, 0))),
               ("so(4,1)", (10, (4, 6, 0))), ("so(3,2)", (10, (6, 4, 0))),
               ("sp(4,R)", (10, (6, 4, 0))), ("so(5,1)", (15, (5, 10, 0))),
               ("so(4,2)", (15, (8, 7, 0))), ("so(3,3)", (15, (9, 6, 0)))]
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


def gen_square_eigs(basis):
    out = []
    for M in basis:
        ev = (M * M).eigenvals()
        out.append({str(k): int(v) for k, v in ev.items()})
    return out


# ---------- module data + block decomposition ----------
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
    if dv == 0:
        return 0, 0
    eb = end_basis(dv)

    def count(symtype):
        cols = []
        for E in eb:
            parts = [flat(op.T * E + E * op) for op in ops]
            parts.append(flat(E - E.T) if symtype == 'sym' else flat(E + E.T))
            cols.append(Matrix.vstack(*parts))
        return len(Matrix.hstack(*cols).nullspace())

    return count('sym'), count('anti')


def cross_invariant_forms(ops_i, ops_j, di, dj):
    """Levi-invariant bilinear pairings B: V_i x V_j -> field, op_i^T B + B op_j = 0."""
    if di == 0 or dj == 0:
        return 0
    eb = []
    for a in range(di):
        for b in range(dj):
            E = zeros(di, dj); E[a, b] = Integer(1); eb.append(E)
    cols = []
    for E in eb:
        parts = [flat(oi.T * E + E * oj) for (oi, oj) in zip(ops_i, ops_j)]
        cols.append(Matrix.vstack(*parts))
    return len(Matrix.hstack(*cols).nullspace())


def module_blocks(ops, dv):
    """Return (blocks, note). blocks = list of column-vector lists spanning
    Levi-invariant subspaces. CARVED: commutant scalar (dim<=1) -> irreducible,
    one block = whole space; else split by rational eigenspaces of the first
    commutant element with >=2 rational eigenvalues (invariant, commutes);
    no rational split -> (None, NAMED-STOP)."""
    if dv == 0:
        return [], "empty"
    C = commutant_basis(ops, dv)
    if len(C) <= 1:
        return [[unit(dv, i) for i in range(dv)]], "irreducible(commutant<=1)"
    for X in C:
        ev = X.eigenvals()
        if len(ev) >= 2 and all(k.is_rational for k in ev.keys()):
            blocks = []
            for lam in sorted(ev.keys(), key=lambda z: (float(z))):
                blocks.append((X - lam * eye(dv)).nullspace())
            return blocks, "rational-eigenspace split (evals={0})".format(
                sorted(str(k) for k in ev.keys()))
    return None, "NAMED-STOP: commutant dim {0}, no rational eigen-split".format(len(C))


def restrict_ops(ops, blockcols):
    dk = len(blockcols)
    if dk == 0:
        return []
    B = Matrix.hstack(*blockcols)
    rops = []
    for op in ops:
        cols = []
        for j in range(dk):
            w = op * blockcols[j]
            sol, params = B.gauss_jordan_solve(w)
            if params.rows * params.cols > 0:
                sol = sol.subs({s: 0 for s in params})
            cols.append(sol)
        rops.append(Matrix.hstack(*cols))
    return rops


def block_matrices(blockcols, Mbasis, n):
    out = []
    for v in blockcols:
        M = zeros(n, n)
        for k in range(len(Mbasis)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * Mbasis[k]
        out.append(M)
    return out


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


def two_wedge(p, q):
    """N=N1+N2 over pairwise-orthogonal isotropic supports e_i+e_{p+i}, i=0..3."""
    n = p + q
    eta = make_eta(p, q)
    x1 = unit(n, 0) + unit(n, p)
    y1 = unit(n, 1) + unit(n, p + 1)
    x2 = unit(n, 2) + unit(n, p + 2)
    y2 = unit(n, 3) + unit(n, p + 3)
    return wedge(x1, y1, eta) + wedge(x2, y2, eta), [x1, y1, x2, y2], eta


def single_wedge(p, q):
    n = p + q
    eta = make_eta(p, q)
    x1 = unit(n, 0) + unit(n, p)
    y1 = unit(n, 1) + unit(n, p + 1)
    return wedge(x1, y1, eta), eta


# ==================== unified analyzer ====================
def analyze_object(tag, N, n, eta):
    print("=== OBJECT: {0} ===".format(tag))
    bas = so_basis(n, eta)
    ok(in_so(N, eta), tag + ": N in so(eta)")
    dep = nilp_depth(N)
    ok(dep is not None, tag + ": N nilpotent")
    print("  N: rank={0} nilp_depth={1}".format(N.rank(), dep))
    c = centralizer(N, bas)
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
        rtag, rinfo = classify_radical(span_basis(c, n), n)
        print("  class=SOLVABLE radical=whole-c type={0} lcs={1} info={2}".format(rtag, lcs, rinfo))
        print("  -- STOP (raw) --")
        return dict(tag=tag, dimc=dc, cls="SOLVABLE")
    fp = fingerprint(core, n)
    rad = killing_radical(core, n)
    reps = complement_reps(core, rad, n)
    ltag, lfp = identify_levi_quotient(reps, rad, n)
    print("  core: dim={0} Killing={1}".format(len(core), fp[1]))
    print("  Levi: {0} dim={1} Killing_sig={2}".format(ltag, lfp[0], lfp[1]))
    sq = gen_square_eigs(reps)
    print("  Levi gen^2 eigenvalues on {0}-dim rep (raw, reduct signature):".format(n))
    for i, e in enumerate(sq):
        print("    g{0}^2 eig={1}".format(i, e))
    rtag, rinfo = classify_radical(rad, n)
    Zr = center_of(rad, n)
    lcs = lower_central_series(rad, n)
    print("  radical: type={0} dim={1} lower_central_series={2} n_central_dirs={3} info={4}".format(
        rtag, len(rad), lcs, len(Zr), rinfo))
    md = module_data(reps, rad, n)
    dv = md["dv"]
    if dv == 0:
        print("  module: dv=0 (no translation module)")
        print("  -- STOP (raw) --")
        return dict(tag=tag, dimc=dc, levi=ltag, radical=rtag)
    cdim = commutant_dim(md["ops"], dv)
    wm = weight_multiset(md["ops"], dv)
    print("  module: dv={0} commutant={1} weights(aggregate)={2}".format(dv, cdim, wm))
    # block decomposition
    blocks, bnote = module_blocks(md["ops"], dv)
    print("  module blocks: {0}".format(bnote))
    if blocks is None:
        print("  [EXT-1/EXT-2 by-block] NAMED-STOP: block split undecided (see note)")
        # still give whole-module forms as fallback raw
        ns, na = invariant_forms_split(md["ops"], dv)
        print("  whole-module invariant forms: sym={0} anti={1}".format(ns, na))
        print("  -- STOP (raw) --")
        return dict(tag=tag, dimc=dc, levi=ltag, radical=rtag, blocks="NAMED-STOP")
    # per-block analysis
    blk_ops = []
    blk_mats = []
    for bi, bc in enumerate(blocks):
        dk = len(bc)
        rops = restrict_ops(md["ops"], bc)
        blk_ops.append(rops)
        blk_mats.append(block_matrices(bc, md["Mbasis"], n))
        bcomm = commutant_dim(rops, dk) if dk > 0 else 0
        bw = weight_multiset(rops, dk)
        ns, na = invariant_forms_split(rops, dk)
        print("  [EXT-1] block{0}: dim={1} commutant={2} weights={3} invforms(sym={4},anti={5})".format(
            bi, dk, bcomm, bw, ns, na))
    # cross-block invariant forms
    for bi in range(len(blocks)):
        for bj in range(bi + 1, len(blocks)):
            cij = cross_invariant_forms(blk_ops[bi], blk_ops[bj], len(blocks[bi]), len(blocks[bj]))
            print("  [EXT-1 cross] block{0}xblock{1}: invariant pairings dim={2}".format(bi, bj, cij))
    # [EXT-2] block bracket table into centre basis
    Zb = md["Z"] if md["Z"] else Zr
    print("  [EXT-2] [module,module] by blocks (coords in centre basis; else flag):")
    for bi in range(len(blocks)):
        for bj in range(bi, len(blocks)):
            mats_i = blk_mats[bi]
            mats_j = blk_mats[bj]
            for a in range(len(mats_i)):
                brange = range(a + 1, len(mats_j)) if bi == bj else range(len(mats_j))
                for b in brange:
                    br = mats_i[a] * mats_j[b] - mats_j[b] * mats_i[a]
                    if br.is_zero_matrix:
                        val = "0"
                    else:
                        cc = coords_in(Zb, br, n)
                        val = "OUTSIDE-CENTRE" if cc is None else [cc[t, 0] for t in range(cc.rows)]
                    print("    [blk{0}.m{1}, blk{2}.m{3}] -> {4}".format(bi, a, bj, b, val))
    print("  -- STOP (raw) --")
    return dict(tag=tag, dimc=dc, levi=ltag, radical=rtag, nblocks=len(blocks))


# ==================== bit-fence: shared pipeline reproduces S988 (4,4) control ====================
print("--- bit-fence: shared pipeline reproduces S988 (4,4) single-wedge control ---")
_N44, _eta44 = single_wedge(4, 4)
_c44 = centralizer(_N44, so_basis(8, _eta44))
ok(len(_c44) == 18, "fence (4,4) control dimc == 18 (got {0})".format(len(_c44)))
_core44 = perfect_core(_c44, 8)
_rad44 = killing_radical(_core44, 8)
_reps44 = complement_reps(_core44, _rad44, 8)
_lt44, _ = identify_levi_quotient(_reps44, _rad44, 8)
_rt44, _ = classify_radical(_rad44, 8)
ok(_rt44 == "h9", "fence (4,4) control radical == h9 (got {0})".format(_rt44))
ok("so(2,1)" in _lt44 and "so(2,2)" in _lt44, "fence (4,4) Levi so(2,1)+so(2,2) (got {0})".format(_lt44))
print("  FENCE (4,4) control: dimc={0} {1} radical={2}".format(len(_c44), _lt44, _rt44))


# ==================== objects (order: cheap before heavy) ====================
print()
_Nc1, _e54 = single_wedge(5, 4)
ok((_Nc1 * _Nc1).is_zero_matrix and _Nc1.rank() == 2, "(C1) so(5,4) single-wedge rank2 N^2=0")
r_c1 = analyze_object("(C1) so(5,4) single-wedge rank2", _Nc1, 9, _e54)

print()
_Na1, _sup1, _e54b = two_wedge(5, 4)
for _a in range(4):
    for _b in range(_a, 4):
        ok((_sup1[_a].T * _e54b * _sup1[_b])[0, 0] == 0, "(A1) supports isotropic/orth ({0},{1})".format(_a, _b))
ok((_Na1 * _Na1).is_zero_matrix and _Na1.rank() == 4, "(A1) so(5,4) two-wedge rank4 N^2=0")
r_a1 = analyze_object("(A1) so(5,4) two-wedge rank4", _Na1, 9, _e54b)

print()
_Nc2, _e55 = single_wedge(5, 5)
ok((_Nc2 * _Nc2).is_zero_matrix and _Nc2.rank() == 2, "(C2) so(5,5) single-wedge rank2 N^2=0")
r_c2 = analyze_object("(C2) so(5,5) single-wedge rank2", _Nc2, 10, _e55)

print()
_Na2, _sup2, _e55b = two_wedge(5, 5)
for _a in range(4):
    for _b in range(_a, 4):
        ok((_sup2[_a].T * _e55b * _sup2[_b])[0, 0] == 0, "(A2) supports isotropic/orth ({0},{1})".format(_a, _b))
ok((_Na2 * _Na2).is_zero_matrix and _Na2.rank() == 4, "(A2) so(5,5) two-wedge rank4 N^2=0")
r_a2 = analyze_object("(A2) so(5,5) two-wedge rank4 [MAIN]", _Na2, 10, _e55b)


# ==================== mutants — one per branch ====================
print()
print("--- mutants (one per branch) ---")
mut_ok = True

if _Na1.rank() == 4 and _Nc1.rank() == 2:
    print("MUTANT br-rank: CAUGHT (two-wedge rank4; single rank2)")
else:
    print("MUTANT br-rank: NOT CAUGHT"); mut_ok = False

_eta90 = make_eta(9, 0)
_iso90 = any(((unit(9, i) + unit(9, j)).T * _eta90 * (unit(9, i) + unit(9, j)))[0, 0] == 0
             for i in range(9) for j in range(9) if i != j)
_iso54 = ((unit(9, 0) + unit(9, 5)).T * _e54 * (unit(9, 0) + unit(9, 5)))[0, 0] == 0
if (not _iso90) and _iso54:
    print("MUTANT br-constructible: CAUGHT ((9,0) no isotropic; (5,4) has)")
else:
    print("MUTANT br-constructible: NOT CAUGHT"); mut_ok = False

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

_hh = diag(3, 1, -1, -3)
_ee = zeros(4, 4); _ee[0, 1] = Integer(1); _ee[1, 2] = Integer(2); _ee[2, 3] = Integer(3)
_ff = zeros(4, 4); _ff[1, 0] = Integer(3); _ff[2, 1] = Integer(2); _ff[3, 2] = Integer(1)
_ci = commutant_dim([_hh, _ee, _ff], 4)
_h2 = diag(1, -1, 1, -1)
_e2 = zeros(4, 4); _e2[0, 1] = Integer(1); _e2[2, 3] = Integer(1)
_f2 = zeros(4, 4); _f2[1, 0] = Integer(1); _f2[3, 2] = Integer(1)
_cr = commutant_dim([_h2, _e2, _f2], 4)
if _ci == 1 and _cr == 4:
    print("MUTANT br-module: CAUGHT (irreducible commutant=1; reducible commutant=4)")
else:
    print("MUTANT br-module: NOT CAUGHT ({0},{1})".format(_ci, _cr)); mut_ok = False

_so3 = so_basis(3, make_eta(3, 0))
_sym3, _anti3 = invariant_forms_split(_so3, 3)
_sp2 = model_sp_w(sp_form(1))
_sym2, _anti2 = invariant_forms_split(_sp2, 2)
if (_sym3, _anti3) == (1, 0) and (_sym2, _anti2) == (0, 1):
    print("MUTANT br-invform: CAUGHT (so(3):sym1/anti0 ; sp(2,R):sym0/anti1)")
else:
    print("MUTANT br-invform: NOT CAUGHT (so3={0} sp2={1})".format((_sym3, _anti3), (_sym2, _anti2))); mut_ok = False

# br-block: a reducible 2+2 module splits into >=2 blocks; an irreducible splits into 1
_blk_red, _ = module_blocks([_h2, _e2, _f2], 4)
_blk_irr, _ = module_blocks([_hh, _ee, _ff], 4)
if _blk_red is not None and len(_blk_red) >= 2 and _blk_irr is not None and len(_blk_irr) == 1:
    print("MUTANT br-block: CAUGHT (reducible->{0} blocks; irreducible->1 block)".format(len(_blk_red)))
else:
    print("MUTANT br-block: NOT CAUGHT ({0},{1})".format(
        None if _blk_red is None else len(_blk_red), None if _blk_irr is None else len(_blk_irr))); mut_ok = False


# ==================== negative control ====================
_g = wedge(unit(9, 0), unit(9, 1), _e54)
_noncomm = not (_g * _Na1 - _Na1 * _g).is_zero_matrix
_ca1 = centralizer(_Na1, so_basis(9, _e54))
_excluded = coords_in(_ca1, _g, 9) is None
ok(_noncomm and _excluded, "neg-control: non-commuting so element excluded from centralizer (A1)")


# ==================== summary ====================
print()
print("SUMMARY: objects=4 (C1/A1 so(5,4), C2/A2 so(5,5)) + bit-fence(4,4) + 6 mutants + neg-control")
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
