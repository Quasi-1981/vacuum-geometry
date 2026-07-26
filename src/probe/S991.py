# -*- coding: utf-8 -*-
# DIM: na (W39 leg-2: two-wedge centralizers in so(5,4) and so(5,5); 0 handles).
#
# ============================================================================
# CARVE CHOICE (stamped BEFORE any counting)
# ----------------------------------------------------------------------------
# Blind probe per TZ_W39_LEG2_BLIND.md.  Pipeline = S988 primitives VERBATIM
# (S933/S937/S940/S950 lineage), plus the two extensions required by the TZ:
#   [EXT-1] Levi-invariant bilinear forms PER irreducible module block
#           (sym/anti dims) + cross-block invariant pairings;
#   [EXT-2] full [mod,mod] bracket table in centre-basis coords, blockwise:
#           (block_i,block_i) and (block_i,block_j) i!=j, element by element.
# Objects, in TZ order (cheap before heavy):
#   (C1) single wedge in so(5,4), rank 2 (control);
#   (A1) two-wedge N=N1+N2 in so(5,4), rank 4, N^2=0, supports pairwise
#        orthogonal totally isotropic;
#   (C2) single wedge in so(5,5), rank 2 (control);
#   (A2) two-wedge in so(5,5), rank 4, N^2=0 (main object).
# For each: c(N) dim/derived/centre, Levi(+)radical split, Levi type (Killing
# exact), radical lower-central series + number of central directions, module
# split into Levi-irreducible blocks over Q (dims, weights, commutant), EXT-1,
# EXT-2, and the reductive-centre / residual so(eta|G) data: Levi gen^2
# eigenvalues on the ambient rep + centre-of-c(N) basis with exact eigenvalues
# (W3-type raw).  Bit-fence FIRST: the pipeline must reproduce the S988 (4,4)
# single-wedge control (dim c = 18, radical h9).  One mutant per
# classification branch; negative control; sign-fence; FORBIDDEN-SCAN; STOP.
# Exact rational arithmetic only (sympy) -- no floats anywhere.
# No targets are assumed -- every column is a raw measurement.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
from sympy import Matrix, Integer, zeros, eye, diag, Rational, symbols, Poly, factor_list

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S991_run.log")
_logf = open(_LOGPATH, "w", encoding="utf-8")


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


# ==================== primitives (VERBATIM from S988/S950/S940/S937/S933) ====================

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


def flat2(M):
    return Matrix(M.rows * M.cols, 1, list(M))


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
    # multi-RHS variant of the S988 primitive (same output, one solve per row)
    k = len(L)
    if k == 0:
        return []
    F = stack_flats(L, n)
    ads = []
    for i in range(k):
        R = Matrix.hstack(*[flat(L[i] * L[j] - L[j] * L[i]) for j in range(k)])
        try:
            sol, params = F.gauss_jordan_solve(R)
        except ValueError:
            return None
        if params.rows * params.cols > 0:
            sol = sol.subs({s: 0 for s in params})
        if not (F * sol - R).is_zero_matrix:
            return None
        ads.append(sol)
    return ads


def killing_matrix(L, n):
    ads = ad_matrices(L, n)
    if ads is None:
        return None
    k = len(L)
    K = zeros(k, k)
    for i in range(k):
        for j in range(i, k):
            v = sum(ads[i][a, b] * ads[j][b, a] for a in range(k) for b in range(k))
            K[i, j] = v
            K[j, i] = v
    return K


def killing_sig(L, n):
    K = killing_matrix(L, n)
    return None if K is None else cong_signature(K)


def killing_radical_from(K, core, n):
    ns = K.nullspace()
    rad = []
    for v in ns:
        M = zeros(n, n)
        for t in range(len(core)):
            if v[t, 0] != 0:
                M = M + v[t, 0] * core[t]
        rad.append(M)
    return span_basis(rad, n)


def killing_radical(core, n):
    K = killing_matrix(core, n)
    return killing_radical_from(K, core, n)


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
    if not ops:
        return dv * dv
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


# ---------- Levi identification menu (Killing signature exact; extended for n=9,10) ----------
LEVI_SIMPLE = [("so(3)", (3, (0, 3, 0))), ("so(2,1)", (3, (2, 1, 0))),
               ("so(4)", (6, (0, 6, 0))), ("so(3,1)", (6, (3, 3, 0))),
               ("so(2,2)", (6, (4, 2, 0))),
               ("sl(3,R)", (8, (5, 3, 0))), ("su(3)", (8, (0, 8, 0))),
               ("su(2,1)", (8, (4, 4, 0))),
               ("so(5)", (10, (0, 10, 0))), ("so(4,1)", (10, (4, 6, 0))),
               ("so(3,2)", (10, (6, 4, 0))), ("sp(4,R)", (10, (6, 4, 0))),
               ("so(6)", (15, (0, 15, 0))), ("so(5,1)", (15, (5, 10, 0))),
               ("so(4,2)", (15, (8, 7, 0))), ("so(3,3)", (15, (9, 6, 0))),
               ("sl(4,R)", (15, (9, 6, 0))),
               ("so(7)", (21, (0, 21, 0))), ("so(6,1)", (21, (6, 15, 0))),
               ("so(5,2)", (21, (10, 11, 0))), ("so(4,3)", (21, (12, 9, 0))),
               ("sp(6,R)", (21, (12, 9, 0)))]
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
    clean = True
    for z in reps:
        cols = []
        for b in range(dv):
            c = coords_in(full, z * Mbasis[b] - Mbasis[b] * z, n)
            if c is None:
                clean = False
                c = zeros(dz + dv, 1)
            cols.append(Matrix(dv, 1, [c[dz + t, 0] for t in range(dv)]))
        ops.append(Matrix.hstack(*cols) if dv > 0 else zeros(0, 0))
    return dict(Z=Z, Mbasis=Mbasis, dv=dv, ops=ops, dz=dz, clean=clean)


def invariant_forms_split(ops, dv):
    """Levi-invariant bilinear forms B on the dv-module: op^T B + B op = 0.
    Return (n_symmetric, n_antisymmetric) exact dims."""
    if dv == 0:
        return 0, 0
    if not ops:
        return dv * (dv + 1) // 2, dv * (dv - 1) // 2
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


def cross_invariant_dim(ropsA, ropsB, ka, kb):
    """Invariant pairings B: V_a x V_b -> Q with ropsA[i]^T B + B ropsB[i] = 0
    simultaneously (same Levi generator index i). Exact dim."""
    if ka == 0 or kb == 0:
        return 0
    if not ropsA or not ropsB:
        return ka * kb
    units = []
    for a in range(ka):
        for b in range(kb):
            E = zeros(ka, kb)
            E[a, b] = Integer(1)
            units.append(E)
    cols = []
    for E in units:
        parts = [flat2(ropsA[i].T * E + E * ropsB[i]) for i in range(len(ropsA))]
        cols.append(Matrix.vstack(*parts))
    return len(Matrix.hstack(*cols).nullspace())


def gen_square_eigs(basis):
    """Eigenvalues of the square of each generator on the ambient rep (raw)."""
    out = []
    for M in basis:
        ev = (M * M).eigenvals()
        out.append({str(k): int(v) for k, v in ev.items()})
    return out


def eig_dict(M):
    ev = M.eigenvals()
    return {str(k): int(v) for k, v in ev.items()}


def nz_entries(M):
    out = []
    for i in range(M.rows):
        for j in range(M.cols):
            if M[i, j] != 0:
                out.append((i, j, M[i, j]))
    return out


# ==================== module block split (exact, over Q) ====================

def restrict_ops(B, ops):
    """B: m x k full-column-rank basis of a subspace; if invariant, return the
    k x k restricted operators, else None."""
    rops = []
    for op in ops:
        R = op * B
        try:
            sol, params = B.gauss_jordan_solve(R)
        except ValueError:
            return None
        if params.rows * params.cols > 0:
            sol = sol.subs({s: 0 for s in params})
        if not (B * sol - R).is_zero_matrix:
            return None
        rops.append(sol)
    return rops


def commutant_basis(ops, m):
    eb = end_basis(m)
    cols = []
    for E in eb:
        cols.append(Matrix.vstack(*[flat2(E * op - op * E) for op in ops]))
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        M = zeros(m, m)
        for k in range(m * m):
            if v[k, 0] != 0:
                M = M + v[k, 0] * eb[k]
        out.append(M)
    return out


def mat_poly(p, A, m):
    Mres = zeros(m, m)
    for c in p.all_coeffs():
        Mres = Mres * A + c * eye(m)
    return Mres


def try_split_with(cand, m):
    """If cand (commutant element) yields invariant subspaces, return list of
    kernel bases; last element None means: [W, None] -> complement needed."""
    lam = symbols('lam')
    p = cand.charpoly(lam)
    fl = factor_list(p.as_expr())
    factors = list(fl[1])
    if len(factors) >= 2:
        subs_list = []
        for (f, e) in factors:
            fp = Poly(f, lam)
            Mf = mat_poly(fp, cand, m) ** e
            ker = Mf.nullspace()
            subs_list.append(Matrix.hstack(*ker))
        return subs_list
    (f, e) = factors[0]
    fp = Poly(f, lam)
    Mf = mat_poly(fp, cand, m)
    if Mf.is_zero_matrix:
        return None
    ker = Mf.nullspace()
    if not ker:
        return None
    W = Matrix.hstack(*ker)
    if W.cols == m:
        return None
    return [W, None]


def invariant_complement(Bw, rops, m):
    """Invariant complement of the invariant subspace with basis Bw (m x k)
    via an exact invariant projection P = Bw*X, X*Bw = I, [P, op] = 0."""
    k = Bw.cols
    units = []
    for a in range(k):
        for b in range(m):
            E = zeros(k, m)
            E[a, b] = Integer(1)
            units.append(E)
    cols = []
    for E in units:
        parts = [flat2(E * Bw)]
        PE = Bw * E
        for op in rops:
            parts.append(flat2(PE * op - op * PE))
        cols.append(Matrix.vstack(*parts))
    rhs = Matrix.vstack(flat2(eye(k)), zeros(m * m * len(rops), 1))
    F = Matrix.hstack(*cols)
    try:
        sol, params = F.gauss_jordan_solve(rhs)
    except ValueError:
        return None
    if params.rows * params.cols > 0:
        sol = sol.subs({s: 0 for s in params})
    if not (F * sol - rhs).is_zero_matrix:
        return None
    X = zeros(k, m)
    idx = 0
    for a in range(k):
        for b in range(m):
            X[a, b] = sol[idx, 0]
            idx += 1
    P = Bw * X
    ns = P.nullspace()
    if not ns:
        return None
    C = Matrix.hstack(*ns)
    if C.cols != m - k:
        return None
    return C


def split_module(ops, m):
    """Split the m-dim module into invariant blocks, irreducible over Q where
    certifiable (commutant dim 1, or field commutant note). Exact throughout.
    Returns (list of basis matrices in top coords, notes)."""
    notes = []

    def rec(B, rops, mloc):
        if mloc == 0:
            return []
        if not rops:
            return [B[:, i] for i in range(mloc)]
        comm = commutant_basis(rops, mloc)
        if len(comm) == 1:
            return [B]
        cands = list(comm)
        for i0 in range(len(comm)):
            for j0 in range(i0 + 1, len(comm)):
                cands.append(comm[i0] + comm[j0])
                cands.append(comm[i0] * comm[j0])
        for cand in cands:
            if (cand - cand[0, 0] * eye(mloc)).is_zero_matrix:
                continue
            sp = try_split_with(cand, mloc)
            if sp is None:
                continue
            if sp[-1] is None:
                W = sp[0]
                C = invariant_complement(W, rops, mloc)
                if C is None:
                    continue
                parts = [W, C]
            else:
                parts = sp
            out = []
            good = True
            for Wp in parts:
                ropsW = restrict_ops(Wp, rops)
                if ropsW is None:
                    good = False
                    break
                out.extend(rec(B * Wp, ropsW, Wp.cols))
            if good:
                return out
        comm_commutes = all((comm[i0] * comm[j0] - comm[j0] * comm[i0]).is_zero_matrix
                            for i0 in range(len(comm)) for j0 in range(i0 + 1, len(comm)))
        notes.append("RAW-NOTE: block dim {0} commutant dim {1} commutative={2} "
                     "no split found (treated as one Q-block)".format(
                         mloc, len(comm), comm_commutes))
        return [B]

    blocks = rec(eye(m), ops, m)
    return blocks, notes


# ==================== unified analyzer (S988 package + EXT-1 + EXT-2) ====================

def analyze_leg2(tag, N, n, eta, expect_rank):
    print("=== OBJECT: {0} ===".format(tag))
    bas = so_basis(n, eta)
    ok(in_so(N, eta), tag + ": N in so(eta)")
    ok((N * N).is_zero_matrix, tag + ": N^2 == 0")
    ok(N.rank() == expect_rank, tag + ": rank == {0} (got {1})".format(expect_rank, N.rank()))
    print("  N: rank={0} nilp_depth(min k: N^k=0)={1}".format(N.rank(), nilp_depth(N)))
    c = centralizer(N, bas)
    ok(all((B * N - N * B).is_zero_matrix for B in c), tag + ": centralizer commutes exactly")
    dc = len(c)
    dcc = len(bracket_basis(c, n))
    dser = derived_series(c, n)
    Zc = center_of(c, n)
    print("  c(N): dim={0} dim[c,c]={1} derived={2} dim_centre(c)={3}".format(dc, dcc, dser, len(Zc)))
    core = perfect_core(c, n)
    print("  perfect_core dim={0}".format(len(core)))
    out = dict(tag=tag, c=c, dimc=dc)
    if not core:
        lcs0 = lower_central_series(c, n)
        rtag0, rinfo0 = classify_radical(span_basis(c, n), n)
        print("  class=SOLVABLE radical=whole c type={0} lcs={1} info={2}".format(rtag0, lcs0, rinfo0))
        print("  -- STOP (raw, no interpretation) --")
        return out
    K = killing_matrix(core, n)
    ks = cong_signature(K)
    print("  core: dim={0} Killing_sig={1}".format(len(core), ks))
    rad = killing_radical_from(K, core, n)
    reps = complement_reps(core, rad, n)
    ltag, lfp = identify_levi_quotient(reps, rad, n)
    print("  Levi: {0}  dim={1} Killing_sig={2}".format(ltag, lfp[0], lfp[1]))
    # residual so(eta|G): Levi gen^2 eigenvalues on ambient rep (raw)
    sq = gen_square_eigs(reps)
    print("  Levi gen^2 eigenvalues on {0}-dim rep (raw, per generator):".format(n))
    for i2, e2 in enumerate(sq):
        print("    g{0}^2 eig={1}".format(i2, e2))
    # radical structure
    rtag, rinfo = classify_radical(rad, n)
    Zr = center_of(rad, n)
    lcs = lower_central_series(rad, n)
    print("  radical: type={0} dim={1} lower_central_series={2} n_central_dirs={3} info={4}".format(
        rtag, len(rad), lcs, len(Zr), rinfo))
    # module
    md = module_data(reps, rad, n)
    ok(md["clean"], tag + ": Levi action closes on Z+module span")
    dv = md["dv"]
    if dv > 0:
        cdim = commutant_dim(md["ops"], dv)
        wm = weight_multiset(md["ops"], dv)
        mc = min_cyclic_dim(md["ops"], dv) if cdim > 1 else dv
        print("  module (whole): dv={0} commutant={1} weights={2} mincyc={3}".format(dv, cdim, wm, mc))
        blocks, bnotes = split_module(md["ops"], dv)
        blocks = sorted(blocks, key=lambda Bb: Bb.cols)
        for note in bnotes:
            print("  " + note)
        dims = [Bb.cols for Bb in blocks]
        ok(sum(dims) == dv, tag + ": sum of block dims == dv ({0} vs {1})".format(sum(dims), dv))
        if cdim == 1:
            ok(len(blocks) == 1, tag + ": commutant=1 -> single block")
        rops_list = []
        allinv = True
        for Bb in blocks:
            rops = restrict_ops(Bb, md["ops"])
            if rops is None:
                allinv = False
                rops = []
            rops_list.append(rops)
        ok(allinv, tag + ": all blocks invariant under Levi ops")
        print("  module blocks (Levi-irreducible over Q): dims={0}".format(dims))
        for bi2 in range(len(blocks)):
            ki = blocks[bi2].cols
            bc = commutant_dim(rops_list[bi2], ki)
            bw = weight_multiset(rops_list[bi2], ki) if rops_list[bi2] else None
            print("    block{0}: dim={1} commutant={2} weights={3}".format(bi2, ki, bc, bw))
            print("    block{0} basis (module coords): {1}".format(
                bi2, [list(blocks[bi2][:, s2]) for s2 in range(ki)]))
        print("  [EXT-1] per-block invariant bilinear forms (sym/anti dims) + cross-block invariants:")
        for bi2 in range(len(blocks)):
            nsym, nanti = invariant_forms_split(rops_list[bi2], blocks[bi2].cols)
            print("    block{0}: sym={1} anti={2}".format(bi2, nsym, nanti))
        for bi2 in range(len(blocks)):
            for bj2 in range(bi2 + 1, len(blocks)):
                ncr = cross_invariant_dim(rops_list[bi2], rops_list[bj2],
                                          blocks[bi2].cols, blocks[bj2].cols)
                print("    cross block{0} x block{1}: dim={2}".format(bi2, bj2, ncr))
        Zb = md["Z"] if md["Z"] else Zr
        print("  centre basis for [EXT-2] coords ({0} dirs):".format(len(Zb)))
        for zi2, Zel in enumerate(Zb):
            print("    zc{0} nz={1}".format(zi2, nz_entries(Zel)))
        bmats = []
        for Bb in blocks:
            colm = []
            for s2 in range(Bb.cols):
                Mm = zeros(n, n)
                for t2 in range(dv):
                    if Bb[t2, s2] != 0:
                        Mm = Mm + Bb[t2, s2] * md["Mbasis"][t2]
                colm.append(Mm)
            bmats.append(colm)
        print("  [EXT-2] full [mod,mod] bracket table in centre-basis coords, blockwise:")
        for bi2 in range(len(blocks)):
            for bj2 in range(bi2, len(blocks)):
                print("    -- pairs (block{0},block{1}) --".format(bi2, bj2))
                for s2 in range(len(bmats[bi2])):
                    tstart = s2 + 1 if bi2 == bj2 else 0
                    for t2 in range(tstart, len(bmats[bj2])):
                        br = bmats[bi2][s2] * bmats[bj2][t2] - bmats[bj2][t2] * bmats[bi2][s2]
                        if br.is_zero_matrix:
                            val = "0"
                        else:
                            cc = coords_in(Zb, br, n)
                            val = "OUTSIDE-CENTRE" if cc is None else [cc[u2, 0] for u2 in range(cc.rows)]
                        print("      [b{0}.{1}, b{2}.{3}] -> {4}".format(bi2, s2, bj2, t2, val))
    else:
        print("  module: dv=0 (radical is its own centre / no translation module)")
    print("  reductive-centre (centre of c(N)) basis, exact eigenvalues (W3-type raw):")
    if not Zc:
        print("    (empty)")
    for zi2, Zel in enumerate(Zc):
        print("    z{0}: nz={1}".format(zi2, nz_entries(Zel)))
        print("    z{0}: eig={1} eig(z^2)={2}".format(zi2, eig_dict(Zel), eig_dict(Zel * Zel)))
    print("  -- STOP (raw, no interpretation) --")
    out.update(dict(core=len(core), levi=ltag, radical=rtag))
    return out


# ==================== bit-fence FIRST: reproduce S988 (4,4) single-wedge control ====================
print("--- bit-fence: pipeline reproduces S988 (4,4) single-wedge control (dimc=18, radical h9) ---")
_eta44 = make_eta(4, 4)
_bas44 = so_basis(8, _eta44)
_Nf = wedge(unit(8, 0) + unit(8, 4), unit(8, 1) + unit(8, 5), _eta44)
_cf = centralizer(_Nf, _bas44)
ok(len(_cf) == 18, "fence dimc (4,4) single-wedge == 18 (got {0})".format(len(_cf)))
_coref = perfect_core(_cf, 8)
_radf = killing_radical(_coref, 8)
_repsf = complement_reps(_coref, _radf, 8)
_ltf, _ = identify_levi_quotient(_repsf, _radf, 8)
_rtf, _ = classify_radical(_radf, 8)
ok(_rtf == "h9", "fence radical (4,4) single-wedge == h9 (got {0})".format(_rtf))
print("  FENCE (4,4) single-wedge: dimc={0} {1} radical={2}".format(len(_cf), _ltf, _rtf))


# ==================== objects (TZ order) ====================
N9 = 9
ETA9 = make_eta(5, 4)
x1_9 = unit(9, 0) + unit(9, 5)
y1_9 = unit(9, 1) + unit(9, 6)
x2_9 = unit(9, 2) + unit(9, 7)
y2_9 = unit(9, 3) + unit(9, 8)
_sup9 = [x1_9, y1_9, x2_9, y2_9]
for _a in range(4):
    for _b in range(_a, 4):
        _v = (_sup9[_a].T * ETA9 * _sup9[_b])[0, 0]
        ok(_v == 0, "so(5,4) supports totally isotropic/orthogonal pair ({0},{1})".format(_a, _b))
NC1 = wedge(x1_9, y1_9, ETA9)
NA1 = NC1 + wedge(x2_9, y2_9, ETA9)

N10 = 10
ETA10 = make_eta(5, 5)
x1_10 = unit(10, 0) + unit(10, 5)
y1_10 = unit(10, 1) + unit(10, 6)
x2_10 = unit(10, 2) + unit(10, 7)
y2_10 = unit(10, 3) + unit(10, 8)
_sup10 = [x1_10, y1_10, x2_10, y2_10]
for _a in range(4):
    for _b in range(_a, 4):
        _v = (_sup10[_a].T * ETA10 * _sup10[_b])[0, 0]
        ok(_v == 0, "so(5,5) supports totally isotropic/orthogonal pair ({0},{1})".format(_a, _b))
NC2 = wedge(x1_10, y1_10, ETA10)
NA2 = NC2 + wedge(x2_10, y2_10, ETA10)

print()
res_C1 = analyze_leg2("(C1) CONTROL single-wedge so(5,4) rank2", NC1, N9, ETA9, 2)
print()
res_A1 = analyze_leg2("(A1) TWO-WEDGE so(5,4) rank4 (N^2=0)", NA1, N9, ETA9, 4)
print()
try:
    res_C2 = analyze_leg2("(C2) CONTROL single-wedge so(5,5) rank2", NC2, N10, ETA10, 2)
except MemoryError:
    res_C2 = None
    print("MACHINE-LIMIT-NOTE(so(5,5) dim45, object C2): exact arithmetic exhausted memory; no float fallback.")
print()
try:
    res_A2 = analyze_leg2("(A2) TWO-WEDGE so(5,5) rank4 (N^2=0) [MAIN]", NA2, N10, ETA10, 4)
except MemoryError:
    res_A2 = None
    print("MACHINE-LIMIT-NOTE(so(5,5) dim45, object A2): exact arithmetic exhausted memory; no float fallback.")


# ==================== mutants — one per classification branch ====================
print()
print("--- mutants (one per branch) ---")
mut_ok = True

# br-depth: N^2=0 (our wedges) vs N^2!=0 (deep Jordan-3 in so(5,4))
_ND9 = wedge(unit(9, 1), unit(9, 0) + unit(9, 5), ETA9)
if (NA1 * NA1).is_zero_matrix and (not (_ND9 * _ND9).is_zero_matrix) and (_ND9 ** 3).is_zero_matrix:
    print("MUTANT br-depth: CAUGHT (N^2=0 for wedges; deep Jordan-3 has N^2!=0, N^3=0)")
else:
    print("MUTANT br-depth: NOT CAUGHT"); mut_ok = False

# br-rank: two-wedge rank 4 vs single-wedge rank 2 (both signatures)
if NA1.rank() == 4 and NC1.rank() == 2 and NA2.rank() == 4 and NC2.rank() == 2:
    print("MUTANT br-rank: CAUGHT (two-wedge rank4; single-wedge rank2; both signatures)")
else:
    print("MUTANT br-rank: NOT CAUGHT"); mut_ok = False

# br-constructible: definite (9,0) has no isotropic e_i+e_j; (5,4) does
_eta90 = make_eta(9, 0)
_iso_exists_90 = any(((unit(9, i) + unit(9, j)).T * _eta90 * (unit(9, i) + unit(9, j)))[0, 0] == 0
                     for i in range(9) for j in range(9) if i != j)
_iso_exists_54 = ((unit(9, 0) + unit(9, 5)).T * ETA9 * (unit(9, 0) + unit(9, 5)))[0, 0] == 0
if (not _iso_exists_90) and _iso_exists_54:
    print("MUTANT br-constructible: CAUGHT ((9,0) no isotropic e_i+e_j; (5,4) has)")
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
_h2d = diag(1, -1, 1, -1)
_e2d = zeros(4, 4); _e2d[0, 1] = Integer(1); _e2d[2, 3] = Integer(1)
_f2d = zeros(4, 4); _f2d[1, 0] = Integer(1); _f2d[3, 2] = Integer(1)
_cr = commutant_dim([_h2d, _e2d, _f2d], 4)
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

# br-blocksplit: block decomposition must find [2,3] for 2+3 and [2,2] for isotypic 2+2
_h2m = Matrix([[1, 0], [0, -1]])
_e2m = Matrix([[0, 1], [0, 0]])
_f2m = Matrix([[0, 0], [1, 0]])
_h3m = diag(2, 0, -2)
_e3m = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
_f3m = Matrix([[0, 0, 0], [2, 0, 0], [0, 2, 0]])


def _bd(A, B2):
    m1 = A.rows; m2 = B2.rows
    Z2 = zeros(m1 + m2, m1 + m2)
    Z2[:m1, :m1] = A
    Z2[m1:, m1:] = B2
    return Z2


_ops23 = [_bd(_h2m, _h3m), _bd(_e2m, _e3m), _bd(_f2m, _f3m)]
_blk23, _n23 = split_module(_ops23, 5)
_d23 = sorted(Bb.cols for Bb in _blk23)
_ops22 = [_bd(_h2m, _h2m), _bd(_e2m, _e2m), _bd(_f2m, _f2m)]
_blk22, _n22 = split_module(_ops22, 4)
_d22 = sorted(Bb.cols for Bb in _blk22)
if _d23 == [2, 3] and _d22 == [2, 2]:
    print("MUTANT br-blocksplit: CAUGHT (2+3 -> blocks [2,3]; isotypic 2+2 -> blocks [2,2])")
else:
    print("MUTANT br-blocksplit: NOT CAUGHT ({0},{1})".format(_d23, _d22)); mut_ok = False

# br-crossform: cross-invariants: standard x standard (sl2) dim 1; standard x adjoint dim 0
_blk22s = sorted(_blk22, key=lambda Bb: Bb.cols)
_blk23s = sorted(_blk23, key=lambda Bb: Bb.cols)
_cAB = None
_c23x = None
if _d22 == [2, 2] and _d23 == [2, 3]:
    _rA = restrict_ops(_blk22s[0], _ops22)
    _rB = restrict_ops(_blk22s[1], _ops22)
    _cAB = cross_invariant_dim(_rA, _rB, 2, 2)
    _r2x = restrict_ops(_blk23s[0], _ops23)
    _r3x = restrict_ops(_blk23s[1], _ops23)
    _c23x = cross_invariant_dim(_r2x, _r3x, 2, 3)
if _cAB == 1 and _c23x == 0:
    print("MUTANT br-crossform: CAUGHT (std x std cross=1; std x adj cross=0)")
else:
    print("MUTANT br-crossform: NOT CAUGHT ({0},{1})".format(_cAB, _c23x)); mut_ok = False


# ==================== negative control ====================
# a specific so(eta) element that does NOT commute with the main two-wedge must
# be excluded from the centralizer span (centralizer is a proper subalgebra).
_g10 = wedge(unit(10, 0), unit(10, 1), ETA10)  # compact rotation in +,+ plane
_noncomm = not (_g10 * NA2 - NA2 * _g10).is_zero_matrix
if res_A2 is not None:
    _excluded = coords_in(res_A2["c"], _g10, 10) is None
else:
    _excluded = coords_in(centralizer(NA2, so_basis(10, ETA10)), _g10, 10) is None
ok(_noncomm and _excluded, "neg-control: non-commuting so element excluded from centralizer")

# sign-fence: eta -> -eta gives same two-wedge (5,4) centralizer invariants
_cM = centralizer(NA1, so_basis(9, -ETA9))
ok(len(res_A1["c"]) == len(_cM), "sign-fence: dim c equal under eta->-eta (got {0} vs {1})".format(
    len(res_A1["c"]), len(_cM)))
_coreP = perfect_core(res_A1["c"], 9)
_coreM = perfect_core(_cM, 9)
ok(len(_coreP) == len(_coreM), "sign-fence: core dim equal under eta->-eta")
ok(killing_sig(_coreP, 9) == killing_sig(_coreM, 9), "sign-fence: core Killing sig equal under eta->-eta")


# ==================== summary ====================
print()
print("SUMMARY: objects=4 (C1/A1/C2/A2) + bit-fence + 8 mutants + neg-control + sign-fence")
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
