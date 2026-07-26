# -*- coding: utf-8 -*-
# DIM: na (action of the 6 non-abelian centralizer cores on the graded pieces of the radical-image
#      filtration, and match against the affine / centrally-extended construct; 0 handles).
#
# ============================================================================
# CARVE / BRANCH CHOICES (stamped BEFORE any counting)
# ----------------------------------------------------------------------------
# FILTRATION.  W := image of the radical on R^n (column span of the radical
# basis).  Alternative W' := span of the wedge support {x,y}.  Measured fact
# (stamped raw): W = ker N  and  W' = im N, so in general W != W' (mismatch).
# Because W need not be isotropic, the graded pieces are taken from the ALWAYS
# well-defined flag  0 ⊆ (W ∩ W^⊥) ⊆ W ⊆ R^n :
#      G1 := W ∩ W^⊥            (isotropic bottom = radical of η|_W)
#      G2 := W / (W ∩ W^⊥)      (nondegenerate middle; η descends nondegenerate)
#      G3 := R^n / W            (top)
# When W is co-isotropic (W^⊥ ⊆ W) this is exactly the TZ's G2 = W^⊥/W-analogue
# (middle nondegenerate piece); sig(η|G2) = (p'',q'') is reported raw.
#
# CONSTRUCT MATCH (step 4).  The translation module is M := rad (if rad abelian)
# or M := rad/Z (if [rad,rad]!=0, Z = center of rad).  The affine construct is
# A := Levi ⋉ M with M abelian ([M,M]=0), Levi acting by its module action —
# this is so(p',q') ⋉ R^{d'} (S921 shape) with (p',q') the signature of the
# Levi-invariant form on M, d' = dim M.  Match by iso-fingerprint
# (dim, Killing-sig, derived-series):
#      fp(core) == fp(A)            -> "affine"        (rad abelian)
#      fp(core) != fp(A) and
#      fp(core/Z) == fp(A), Z!=0    -> "central-extended"
#      otherwise                    -> "none"
# EVERY branch (affine / extended / none / class-1 yes-no / equivariance yes-no)
# carries its own CAUGHT mutant (lesson J-0415).
#
# Fences: blindness, exact arithmetic, FORBIDDEN-SCAN, bit-fence S937/S940
# (cores re-built by the same code), STOP after the table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import itertools

from sympy import Matrix, Integer, Rational, zeros, eye, diag, symbols, gcd, linsolve

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S943_run.log")
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


# ==================== primitives (VERBATIM) ====================

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


def bracket_basis(L, n):
    brs = []
    for a in range(len(L)):
        for b in range(a + 1, len(L)):
            brs.append(L[a] * L[b] - L[b] * L[a])
    return span_basis(brs, n)


def perfect_core(L, n):
    cur = span_basis(L, n)
    while len(cur) > 0:
        nxt = bracket_basis(cur, n)
        if len(nxt) == len(cur):
            return cur
        cur = nxt
    return []


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


def killing_sig_mats(L, n):
    ads = ad_matrices(L, n)
    if ads is None:
        return None
    k = len(L)
    K = zeros(k, k)
    for i in range(k):
        for j in range(k):
            K[i, j] = (ads[i] * ads[j]).trace()
    return cong_signature(K)


def killing_radical(core, n):
    ads = ad_matrices(core, n)
    k = len(core)
    K = zeros(k, k)
    for i in range(k):
        for j in range(k):
            K[i, j] = (ads[i] * ads[j]).trace()
    ns = K.nullspace()
    rad = []
    for v in ns:
        M = zeros(n, n)
        for t in range(k):
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


def d0_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 2 and len(Km) >= 2:
        return unit(n, Kp[0]) + unit(n, Km[0]), unit(n, Kp[1]) + unit(n, Km[1])
    return None


def d1_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 1 and len(Km) >= 1 and len(K) >= 3:
        a, b = Kp[0], Km[0]
        rest = [c for c in K if c not in (a, b)]
        c = rest[0]
        return unit(n, a) + unit(n, b), unit(n, c)
    return None


# ==================== vector-subspace helpers on R^n ====================

def vspan_basis(cols, n):
    out = []
    F = zeros(n, 0)
    r = 0
    for v in cols:
        F2 = Matrix.hstack(F, v)
        if F2.rank() > r:
            out.append(v)
            F = F2
            r += 1
    return out


def vdim(cols, n):
    return Matrix.hstack(*cols).rank() if cols else 0


def image_span(mats, n):
    cols = []
    for M in mats:
        for i in range(M.cols):
            cols.append(M[:, i])
    return vspan_basis(cols, n)


def eta_perp(Wcols, eta, n):
    if not Wcols:
        return [unit(n, i) for i in range(n)]
    A = Matrix.hstack(*[eta * w for w in Wcols]).T
    return A.nullspace()


def contained(Acols, Bcols, n):
    FB = Matrix.hstack(*Bcols) if Bcols else zeros(n, 0)
    rb = FB.rank()
    for v in Acols:
        if Matrix.hstack(FB, v).rank() > rb:
            return False
    return True


def same_subspace(A, B, n):
    return contained(A, B, n) and contained(B, A, n)


def intersect_sub(A, B, n):
    # basis of span(A) ∩ span(B) via nullspace of [A | -B]
    if not A or not B:
        return []
    FA = Matrix.hstack(*A)
    FB = Matrix.hstack(*B)
    H = Matrix.hstack(FA, -FB)
    ns = H.nullspace()
    ka = FA.cols
    out = []
    for v in ns:
        w = FA * v[0:ka, 0]
        if not w.is_zero_matrix:
            out.append(w)
    return vspan_basis(out, n)


def eta_sig_on(cols, eta, n):
    if not cols:
        return (0, 0, 0)
    G = Matrix(len(cols), len(cols), lambda a, b: (cols[a].T * eta * cols[b])[0, 0])
    return cong_signature(G)


def preserves_sub(M, cols, n):
    # M * span(cols) ⊆ span(cols)
    if not cols:
        return True
    FB = Matrix.hstack(*cols)
    rb = FB.rank()
    for v in cols:
        if Matrix.hstack(FB, M * v).rank() > rb:
            return False
    return True


# ==================== structure-constant fingerprint (abstract algebras) ====================

def sc_from_mats(basis, n):
    # {(i,j): tuple of len(basis) coeffs} for i<j, [e_i,e_j] in basis coords
    tab = {}
    k = len(basis)
    for i in range(k):
        for j in range(i + 1, k):
            c = coords_in(basis, basis[i] * basis[j] - basis[j] * basis[i], n)
            if c is None:
                return None
            tab[(i, j)] = tuple(c[t, 0] for t in range(k))
    return tab


def ads_from_sc(sc, dim):
    ads = []
    for i in range(dim):
        M = zeros(dim, dim)
        for j in range(dim):
            if i == j:
                continue
            if i < j:
                co = sc[(i, j)]
                s = Integer(1)
            else:
                co = sc[(j, i)]
                s = Integer(-1)
            for k in range(dim):
                M[k, j] = s * co[k]
        ads.append(M)
    return ads


def bracket_coords(u, v, sc, dim):
    # [u,v] in coords, u,v coord vectors
    out = [Integer(0)] * dim
    for i in range(dim):
        if u[i] == 0:
            continue
        for j in range(dim):
            if v[j] == 0 or i == j:
                continue
            if i < j:
                co = sc[(i, j)]
                s = Integer(1)
            else:
                co = sc[(j, i)]
                s = Integer(-1)
            for k in range(dim):
                out[k] += s * u[i] * v[j] * co[k]
    return out


def derived_series_sc(sc, dim):
    # coords-space derived series dims
    def rankcols(vecs):
        if not vecs:
            return 0
        return Matrix.hstack(*[Matrix(dim, 1, w) for w in vecs]).rank()

    def span_reduce(vecs):
        out = []
        F = zeros(dim, 0)
        r = 0
        for w in vecs:
            col = Matrix(dim, 1, w)
            F2 = Matrix.hstack(F, col)
            if F2.rank() > r:
                out.append(w)
                F = F2
                r += 1
        return out

    cur = span_reduce([[Integer(1) if t == i else Integer(0) for t in range(dim)] for i in range(dim)])
    dims = [len(cur)]
    while len(cur) > 0:
        brs = []
        for a in range(len(cur)):
            for b in range(a + 1, len(cur)):
                brs.append(bracket_coords(cur[a], cur[b], sc, dim))
        nxt = span_reduce(brs)
        if len(nxt) == len(cur):
            break
        dims.append(len(nxt))
        cur = nxt
    return dims


def dim_derived_sc(sc, dim):
    brs = []
    for (i, j) in sc:
        brs.append(list(sc[(i, j)]))
    if not brs:
        return 0
    return Matrix.hstack(*[Matrix(dim, 1, w) for w in brs]).rank()


def killing_sig_sc(sc, dim):
    ads = ads_from_sc(sc, dim)
    K = zeros(dim, dim)
    for i in range(dim):
        for j in range(dim):
            K[i, j] = (ads[i] * ads[j]).trace()
    return cong_signature(K)


def fp_sc(sc, dim):
    return (dim, killing_sig_sc(sc, dim), tuple(derived_series_sc(sc, dim)), dim_derived_sc(sc, dim))


# ---------- explicit Levi (from S937) ----------

def levi_explicit(core, rad, n):
    reps = complement_reps(core, rad, n)
    kq = len(reps)
    dr = len(rad)
    full = reps + rad
    cstr = {}
    for i in range(kq):
        for j in range(kq):
            co = coords_in(full, reps[i] * reps[j] - reps[j] * reps[i], n)
            cstr[(i, j)] = [co[t, 0] for t in range(kq)]
    us = [[symbols('u_%d_%d' % (i, a)) for a in range(dr)] for i in range(kq)]
    t = []
    for i in range(kq):
        Ti = zeros(n, n)
        for a in range(dr):
            Ti = Ti + us[i][a] * rad[a]
        t.append(Ti)
    eqs = []
    for i in range(kq):
        for j in range(i + 1, kq):
            resid = (reps[i] * t[j] - t[j] * reps[i]) + (t[i] * reps[j] - reps[j] * t[i]) \
                    + (reps[i] * reps[j] - reps[j] * reps[i])
            for k in range(kq):
                resid = resid - cstr[(i, j)][k] * (reps[k] + t[k])
            for e in range(n):
                for f in range(n):
                    ent = resid[e, f].expand()
                    if ent != 0:
                        eqs.append(ent)
    unknowns = [us[i][a] for i in range(kq) for a in range(dr)]
    sol = linsolve(eqs, unknowns)
    if not sol:
        return None
    solset = list(sol)[0]
    freevars = set()
    for expr in solset:
        freevars |= expr.free_symbols
    zsub = {s: Integer(0) for s in freevars}
    vals = {unknowns[idx]: solset[idx] for idx in range(len(unknowns))}
    levi = [(reps[i] + t[i]).subs(vals).subs(zsub) for i in range(kq)]
    lb = span_basis(levi, n)
    if len(lb) != kq:
        return None
    for a in range(kq):
        for b in range(kq):
            if coords_in(lb, lb[a] * lb[b] - lb[b] * lb[a], n) is None:
                return None
    return lb


# ---------- build affine construct A = Levi ⋉ M (M abelian) via SC ----------

def build_affine_sc(core, rad, n):
    # A = explicit-Levi ⋉ M with M abelian; [Levi,Levi] in Levi, [Levi,M] projected to M, [M,M]:=0
    L = levi_explicit(core, rad, n)
    if L is None:
        return None, 0
    Z = span_basis(center_of(rad, n), n)
    drr = len(bracket_basis(rad, n))
    M = list(rad) if drr == 0 else complement_reps(rad, Z, n)
    Zc = [] if drr == 0 else Z
    full = L + M + Zc
    kq = len(L)
    dm = len(M)
    dim = kq + dm
    gens = L + M
    sc = {}
    for i in range(dim):
        for j in range(i + 1, dim):
            c = coords_in(full, gens[i] * gens[j] - gens[j] * gens[i], n)
            if c is None:
                return None, 0
            if i < kq and j < kq:                       # Levi-Levi -> Levi part
                sc[(i, j)] = tuple(c[t, 0] if t < kq else Integer(0) for t in range(dim))
            elif i < kq <= j:                            # Levi-module -> module part (drop Z)
                sc[(i, j)] = tuple(c[t, 0] if kq <= t < kq + dm else Integer(0) for t in range(dim))
            else:                                        # module-module := 0
                sc[(i, j)] = tuple(Integer(0) for _ in range(dim))
    return sc, dim


# ==================== the 6 cores ====================
CORES = [
    ("N", 3, 2, "rank0", d0_on),
    ("N", 5, 1, "rank1", d1_on),
    ("N", 4, 2, "rank1", d1_on),
    ("N", 3, 3, "rank1", d1_on),
    ("N", 4, 2, "rank0", d0_on),
    ("N", 3, 3, "rank0", d0_on),
]
FENCE = {   # bit-fence S937/S940
    ("N", 3, 2, "rank0"): dict(dimc=6, core=6, rad=3, drr=1),
    ("N", 5, 1, "rank1"): dict(dimc=7, core=6, rad=3, drr=0),
    ("N", 4, 2, "rank1"): dict(dimc=7, core=6, rad=3, drr=0),
    ("N", 3, 3, "rank1"): dict(dimc=7, core=6, rad=3, drr=0),
    ("N", 4, 2, "rank0"): dict(dimc=9, core=8, rad=5, drr=1),
    ("N", 3, 3, "rank0"): dict(dimc=9, core=8, rad=5, drr=1),
}


def quotient_by_ideal_sc(core, ideal, n):
    # SC of core/ideal in a complement-reps basis
    F = stack_flats(ideal, n)
    rk = F.rank()
    reps = []
    for M in core:
        F2 = Matrix.hstack(F, flat(M))
        if F2.rank() > rk:
            reps.append(M)
            F = F2
            rk += 1
    full = reps + ideal
    kq = len(reps)
    sc = {}
    for i in range(kq):
        for j in range(i + 1, kq):
            c = coords_in(full, reps[i] * reps[j] - reps[j] * reps[i], n)
            sc[(i, j)] = tuple(c[t, 0] for t in range(kq))
    return sc, kq


print("--- W33-нога-2: graded action + construct match ---")
for (kind, p, q, cls, ctor) in CORES:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    xy = ctor(list(range(n)), p, n)
    x, y = xy
    N = wedge(x, y, eta)
    c = centralizer(N, bas)
    core = perfect_core(c, n)
    rad = killing_radical(core, n)
    reps = complement_reps(core, rad, n)
    fk = (kind, p, q, cls)
    exp = FENCE[fk]
    ok(len(c) == exp["dimc"], "bit-fence dimc " + str(fk))
    ok(len(core) == exp["core"], "bit-fence core " + str(fk))
    ok(len(rad) == exp["rad"], "bit-fence rad " + str(fk))
    ok(len(bracket_basis(rad, n)) == exp["drr"], "bit-fence [rad,rad] " + str(fk))

    # --- filtration ---
    W = image_span(rad, n)
    Wp = [x, y]
    dW = vdim(W, n)
    match_WWp = same_subspace(W, Wp, n)
    kerN = N.nullspace()
    W_is_kerN = same_subspace(W, kerN, n) if kerN else (dW == 0)
    Wperp = eta_perp(W, eta, n)
    # core-invariance of W and W^⊥
    invW = all(preserves_sub(M, W, n) for M in core)
    invWp = all(preserves_sub(M, Wperp, n) for M in core)
    ok(invW, "W core-invariant " + str(fk))
    ok(invWp, "Wperp core-invariant " + str(fk))
    # graded flag 0 ⊆ W∩Wperp ⊆ W ⊆ R^n
    G1 = intersect_sub(W, Wperp, n)
    dG1 = vdim(G1, n) if G1 else 0
    dG2 = dW - dG1
    dG3 = n - dW
    # sig of η on G2 = W / (W∩Wperp): take W basis, quotient out G1 -> nondeg part
    sigW = eta_sig_on(W, eta, n)          # (pos,neg,zer); zer = dim G1
    sigG2 = (sigW[0], sigW[1])            # nondegenerate part
    ok(sigW[2] == dG1, "sig(η|W) zero-count = dim G1 " + str(fk))

    # --- translation module M and construct ---
    Z = span_basis(center_of(rad, n), n)
    drr = len(bracket_basis(rad, n))
    if drr == 0:
        Mmod = rad
    else:
        Mmod = complement_reps(rad, Z, n)   # rad/Z reps (non-central)
    # radical action class: R·G3⊆G2, R·G2⊆G1 (class-1 nilpotent on the flag)
    # flag subspaces as column-spans: F0=G1(=W∩Wperp), F1=W, F2=R^n
    Wperp_b = list(Wperp)
    class1 = True
    for R in rad:
        # R·W ⊆ W∩Wperp ? (drop one grade)   and   R·R^n ⊆ W ?
        if not contained([R * v for v in W], G1 if G1 else [zeros(n, 1)], n):
            class1 = False
        allcols = [unit(n, i) for i in range(n)]
        if not contained([R * v for v in allcols], W, n):
            class1 = False
    # Levi-equivariance: [L,R] ∈ rad and acts consistently (machine: [L,R] preserves flag same way)
    equivar = True
    for L in reps:
        for R in rad:
            LR = L * R - R * L
            if coords_in(rad, LR, n) is None:
                equivar = False
    ok(equivar, "radical is a Levi-submodule ([L,R]∈rad) " + str(fk))

    # --- construct fingerprint match ---
    A_sc, A_dim = build_affine_sc(core, rad, n)
    fpA = fp_sc(A_sc, A_dim)
    core_sc = sc_from_mats(core, n)
    fpCore = fp_sc(core_sc, len(core))
    if drr > 0:
        cz_sc, cz_dim = quotient_by_ideal_sc(core, Z, n)
        fpCoreZ = fp_sc(cz_sc, cz_dim)
    else:
        fpCoreZ = None
    if fpCore == fpA:
        match = "affine"
    elif fpCoreZ is not None and fpCoreZ == fpA and len(Z) > 0:
        match = "central-extended"
    else:
        match = "none"

    # dimension audit
    audit = (dG1 + dG2 + dG3 == n) and (len(core) == len(reps) + len(rad))
    ok(audit, "dimension audit " + str(fk))

    print(("CORE ({0},{1}) {2} | dimW={3} W=W'?{4} W=kerN?{5} | invar W/Wperp={6}/{7} | "
           "G1={8} G2={9} G3={10} sigG2={11} | rad: class-1={12} Levi-equiv={13} | "
           "M(transl)dim={14} [rad,rad]={15} | bt-match={16} | audit={17}").format(
        p, q, cls, dW, match_WWp, W_is_kerN, invW, invWp,
        dG1, dG2, dG3, sigG2, class1, equivar, len(Mmod), drr, match, "OK" if audit else "FAIL"))


# ==================== mutants — one per classification branch ====================
print("--- mutants (one per branch) ---")
mut_ok = True


def _core_of(p, q, ctor):
    n = p + q
    eta = make_eta(p, q)
    N = wedge(*ctor(list(range(n)), p, n), eta)
    core = perfect_core(centralizer(N, so_basis(n, eta)), n)
    rad = killing_radical(core, n)
    reps = complement_reps(core, rad, n)
    return core, rad, reps, n


# branch AFFINE: an abelian-radical core matches the affine construct
_core, _rad, _reps, _n = _core_of(5, 1, d1_on)
_A, _Ad = build_affine_sc(_core, _rad, _n)
_fpA = fp_sc(_A, _Ad)
_fpC = fp_sc(sc_from_mats(_core, _n), len(_core))
if _fpC == _fpA:
    print("MUTANT br-affine: CAUGHT (abelian-rad (5,1)r1 fp == affine construct fp)")
else:
    print("MUTANT br-affine: NOT CAUGHT"); mut_ok = False

# branch EXTENDED: a Heisenberg-radical core does NOT match affine but core/Z does
_core, _rad, _reps, _n = _core_of(4, 2, d0_on)
_Z = span_basis(center_of(_rad, _n), _n)
_Mmod = complement_reps(_rad, _Z, _n)
_A, _Ad = build_affine_sc(_core, _rad, _n)
_fpA = fp_sc(_A, _Ad)
_fpC = fp_sc(sc_from_mats(_core, _n), len(_core))
_czsc, _czd = quotient_by_ideal_sc(_core, _Z, _n)
_fpCZ = fp_sc(_czsc, _czd)
if _fpC != _fpA and _fpCZ == _fpA and len(_Z) > 0:
    print("MUTANT br-extended: CAUGHT (h5-rad (4,2)r0: core != affine, core/Z == affine)")
else:
    print("MUTANT br-extended: NOT CAUGHT (fpC==fpA:{0}, fpCZ==fpA:{1})".format(_fpC == _fpA, _fpCZ == _fpA))
    mut_ok = False

# branch NONE: a deliberately-wrong construct (Levi of wrong signature) does NOT match
_core, _rad, _reps, _n = _core_of(5, 1, d1_on)     # core has Levi so(3)
_wrongLevi = so_basis(3, make_eta(2, 1))           # so(2,1) instead of so(3)
_Tr = [unit(4, 0), unit(4, 1), unit(4, 2)]         # dummy module basis (unused shape)
# build a mismatched affine so(2,1)⋉R^3 in S921 (4x4) rep
_hats = [Mm.row_join(zeros(3, 1)).col_join(zeros(1, 4)) for Mm in _wrongLevi]
_Ts = []
for i in range(3):
    T = zeros(4, 4); T[i, 3] = Integer(1); _Ts.append(T)
_wrong = _hats + _Ts
_fpWrong = fp_sc(sc_from_mats(_wrong, 4), len(_wrong))
_fpC = fp_sc(sc_from_mats(_core, _n), len(_core))
if _fpC != _fpWrong:
    print("MUTANT br-none: CAUGHT (so(3)-core fp != so(2,1)⋉R^3 construct fp)")
else:
    print("MUTANT br-none: NOT CAUGHT"); mut_ok = False

# branch CLASS-1 (yes): every radical element drops the flag by one grade (R·W ⊆ G1);
#   a genuine (explicit) Levi element preserves the middle grade, so some L·W ⊄ G1.
_core, _rad, _reps, _n = _core_of(4, 2, d0_on)
_W = image_span(_rad, _n)
_G1 = intersect_sub(_W, eta_perp(_W, make_eta(4, 2), _n), _n)
_G1c = _G1 if _G1 else [zeros(_n, 1)]
_allcols = [unit(_n, i) for i in range(_n)]
_rad_c1 = all(contained([R * v for v in _allcols], _W, _n) for R in _rad)   # R·R^n ⊆ W
_L = levi_explicit(_core, _rad, _n)
_levi_not_c1 = any(not contained([L * v for v in _allcols], _W, _n) for L in _L)   # L·R^n ⊄ W
if _rad_c1 and _levi_not_c1:
    print("MUTANT br-class1: CAUGHT (all radical drop top grade R^n->W; some Levi element does NOT)")
else:
    print("MUTANT br-class1: NOT CAUGHT (rad_c1={0}, levi_not_c1={1})".format(_rad_c1, _levi_not_c1))
    mut_ok = False

# branch EQUIVARIANCE (yes): [Levi,rad] stays in rad; [Levi,Levi] does NOT land in rad
_core, _rad, _reps, _n = _core_of(3, 2, d0_on)
_eq_in = all(coords_in(_rad, L * R - R * L, _n) is not None for L in _reps for R in _rad)
_ll_out = any(coords_in(_rad, _reps[a] * _reps[b] - _reps[b] * _reps[a], _n) is None
              for a in range(len(_reps)) for b in range(len(_reps)) if a != b)
if _eq_in and _ll_out:
    print("MUTANT br-equivar: CAUGHT ([Levi,rad]⊆rad; [Levi,Levi]⊄rad)")
else:
    print("MUTANT br-equivar: NOT CAUGHT (eq_in={0}, ll_out={1})".format(_eq_in, _ll_out))
    mut_ok = False


# ==================== sign-fence {eta,-eta} ====================
_p, _q = 4, 2; _nn = 6
_N = wedge(*d0_on(list(range(_nn)), _p, _nn), make_eta(_p, _q))
_coreA = perfect_core(centralizer(_N, so_basis(_nn, make_eta(_p, _q))), _nn)
_coreB = perfect_core(centralizer(_N, so_basis(_nn, -make_eta(_p, _q))), _nn)
_radA = killing_radical(_coreA, _nn)
_radB = killing_radical(_coreB, _nn)
ok(vdim(image_span(_radA, _nn), _nn) == vdim(image_span(_radB, _nn), _nn),
   "sign-fence: dim W invariant under eta->-eta")


# ==================== summary ====================
print("SUMMARY: cores={0}".format(len(CORES)))
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
sys.stdout = _tee.real
_logf.close()
sys.exit(_exit)
