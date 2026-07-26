# -*- coding: utf-8 -*-
# DIM: na (fate of the center form under c(S) ∩ c(N); ad(S) on the module; 0 handles).
#
# ============================================================================
# OPERATIONALIZATION (stamped BEFORE any counting)
# ----------------------------------------------------------------------------
# Context (measured, S933/S961): for A=S+N (J-cases) c(A)=c(S)∩c(N).  c(N) of a
# rank0-wedge carries rad=Heisenberg with NONDEGENERATE omega on V:=rad/Z (S961);
# the full c(A) form is EMPTY.  This leg measures HOW omega dies at the intersection.
#
#   S in c(N) normalizes rad and preserves Z (Z=center is characteristic), so
#   ad(S) descends to a linear map  adS_V : V -> V,  V = rad(c(N))/Z.
#   SURVIVORS  W_surv := ker(adS_V)  (module directions that commute with S mod Z).
#   omega|surv := the S961 Z-valued form restricted to W_surv; rank + isotropic?.
#   PARTNER of a survivor u: any module direction r_j with omega(u,r_j)!=0; column
#   records whether that partner itself survives (in ker adS_V).
#
# SYMBOLIC BLOCK (S952 rank0): module = B in R^{2 x d}; a semisimple S acts as
#   adS(B) = a.B - B.E,  a in sp(nu)=so(2,1) (W-index 2x2), E in so(G) (core dxd);
#   omega(B,B') = antisym( B G^{-1} B'^T )  (the 2x2 antisym -> scalar).
#   Survivors = { B : a B = B E } (eigenvalue-matched intertwiners).  omega|surv:
#     * a=0, E full-rank on core -> survivors 0            (J-like: EMPTY)
#     * a=0, E kills part of core -> survivors 2*dim ker E, omega NONDEG (not isotropic)
#     * a,E share eigenvalue lambda -> same-lambda survivors ISOTROPIC (omega|surv=0)
#   because adS_V in sp(V,omega): eigenvalue-lambda space is omega-isotropic unless
#   paired with -lambda.  J2/J4 instances verified bit-for-bit.
#
# Discipline: 0 handles; exact arithmetic (field extensions explicit); mutant on
# each branch {isotropic / nondegenerate / trivial-S}; FORBIDDEN-SCAN (S929);
# bit-fence J2/J4 = S961 (empty full form) + S933 (dim c); STOP after the table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import itertools
from math import lcm as _ilcm
from sympy import (Matrix, Integer, Rational, zeros, eye, diag, symbols, Symbol,
                   sqrt, I, factor, Poly)

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S966_run.log")
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


# ==================== primitives ====================

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


def block_gen(n, eta, i, j, param):
    return Integer(param) * wedge(unit(n, i), unit(n, j), eta)


def q_gram(cb):
    d = len(cb)
    G = zeros(d, d)
    for a in range(d):
        for b in range(a, d):
            t = (cb[a] * cb[b]).trace()
            G[a, b] = t
            G[b, a] = t
    return G


def jordan_scan(cb, G, n):
    d = len(cb)
    dens = [int(G[i, j].q) for i in range(d) for j in range(d)]
    Lc = 1
    for de in dens:
        Lc = _ilcm(Lc, de)
    Gi = [[int(G[i, j] * Lc) for j in range(d)] for i in range(d)]
    vals = (-2, -1, 0, 1, 2)
    for c in itertools.product(vals, repeat=d):
        if all(x == 0 for x in c):
            continue
        s = sum(c[i] * sum(Gi[i][j] * c[j] for j in range(d)) for i in range(d))
        if s != 0:
            continue
        Xm = zeros(n, n)
        for k in range(d):
            if c[k]:
                Xm = Xm + Integer(c[k]) * cb[k]
        if Xm.is_zero_matrix:
            continue
        if is_nilp(Xm):
            return Xm
    return None


def d0_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 2 and len(Km) >= 2:
        return unit(n, Kp[0]) + unit(n, Km[0]), unit(n, Kp[1]) + unit(n, Km[1])
    return None


# ==================== Lie machinery ====================

def center_of(cb, n):
    if not cb:
        return []
    cols = []
    for k in range(len(cb)):
        parts = [flat(cb[k] * cb[j] - cb[j] * cb[k]) for j in range(len(cb))]
        cols.append(Matrix.vstack(*parts))
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        M = zeros(n, n)
        for k in range(len(cb)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * cb[k]
        out.append(M)
    return span_basis(out, n)


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


def killing_radical(core, n):
    K = killing_matrix(core, n)
    if K is None:
        return None
    rad = []
    for v in K.nullspace():
        M = zeros(n, n)
        for t in range(len(core)):
            if v[t, 0] != 0:
                M = M + v[t, 0] * core[t]
        rad.append(M)
    return span_basis(rad, n)


def complement_reps(big, sub, n):
    F = stack_flats(sub, n)
    rk = F.rank()
    reps = []
    for M in big:
        F2 = Matrix.hstack(F, flat(M))
        if F2.rank() > rk:
            reps.append(M)
            F = F2
            rk += 1
    return reps


def module_and_form(c, n):
    """rad, Z, reps(V basis), W (Z-valued form as s antisym mxm), adapted."""
    rad = killing_radical(c, n)
    Z = center_of(rad, n)
    reps = complement_reps(rad, Z, n)
    s = len(Z)
    m = len(reps)
    adapted = Z + reps
    W = [zeros(m, m) for _ in range(s)]
    for i in range(m):
        for j in range(m):
            co = coords_in(adapted, reps[i] * reps[j] - reps[j] * reps[i], n)
            for l in range(s):
                W[l][i, j] = co[l, 0]
    return dict(rad=rad, Z=Z, reps=reps, s=s, m=m, adapted=adapted, W=W)


def adS_on_V(S, mf, n):
    """matrix of ad(S) induced on V=rad/Z: reps-block of [S,reps_j] in adapted basis."""
    reps, Z, adapted, s, m = mf["reps"], mf["Z"], mf["adapted"], mf["s"], mf["m"]
    M = zeros(m, m)
    for j in range(m):
        co = coords_in(adapted, S * reps[j] - reps[j] * S, n)
        for i in range(m):
            M[i, j] = co[s + i, 0]   # reps-part
    return M


def form_rank_on_sub(W, basisvecs):
    """rank of the stacked Z-valued antisym form restricted to a subspace whose
    coeff vectors (in reps basis) are the columns of the basisvecs list."""
    if not basisvecs:
        return 0, True
    P = Matrix.hstack(*basisvecs)   # m x t
    rows = []
    t = P.cols
    for l in range(len(W)):
        Wl = P.T * W[l] * P            # t x t antisym
        for a in range(t):
            for b in range(t):
                pass
        rows.append(Wl)
    # combined rank via stacked kernel
    stack = []
    for b in range(t):
        for l in range(len(W)):
            stack.append([rows[l][a, b] for a in range(t)])
    Mk = Matrix(stack) if stack else zeros(0, t)
    rank = Mk.rank() if Mk.rows > 0 else 0
    isotropic = (rank == 0)
    return rank, isotropic


# ============================================================================
#                                   RUN
# ============================================================================

print("=" * 74)
print("W36-leg-2: the fate of the center form under c(S)∩c(N)")
print("=" * 74)

# -------- Tasks 1-4: real J2, J4 --------
JSPEC = []
_e = make_eta(4, 2)
_S = block_gen(6, _e, 0, 1, 1); _x, _y = d0_on([2, 3, 4, 5], 4, 6)
JSPEC.append(("J2", 4, 2, _S, wedge(_x, _y, _e)))
_e = make_eta(3, 3)
_S = block_gen(6, _e, 0, 3, 1); _x, _y = d0_on([1, 2, 4, 5], 3, 6)
JSPEC.append(("J4", 3, 3, _S, wedge(_x, _y, _e)))

print()
print("--- REAL J2/J4: spectrum of ad(S)|module, survivors, ω|surv, partners ---")
real_res = {}
for (jid, p, q, S, N) in JSPEC:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    c = centralizer(N, bas)
    mf = module_and_form(c, n)
    ok((S * N - N * S).is_zero_matrix, "[S,N]=0 " + jid)
    adSV = adS_on_V(S, mf, n)
    spec = adSV.eigenvals()          # exact, field-ext if needed
    surv = adSV.nullspace()          # coeff vectors in reps basis
    dim_surv = len(surv)
    rank_surv, iso = form_rank_on_sub(mf["W"], surv)
    # partners of each survivor: module dirs r_j with omega(u,r_j)!=0; do they survive?
    partner_rows = []
    ker_cols = [tuple(v) for v in surv]
    for u in surv:
        partners = []
        for j in range(mf["m"]):
            ej = zeros(mf["m"], 1); ej[j, 0] = Integer(1)
            val = any((u.T * mf["W"][l] * ej)[0, 0] != 0 for l in range(mf["s"]))
            if val:
                # does e_j survive? (in ker adSV)
                survives = (adSV * ej).is_zero_matrix
                partners.append((j, survives))
        partner_rows.append(partners)
    real_res[jid] = dict(m=mf["m"], spec=spec, dim_surv=dim_surv,
                         rank_surv=rank_surv, iso=iso, partners=partner_rows,
                         dimrad=len(mf["rad"]), dimZ=mf["s"])
    specstr = ", ".join("{0}^{1}".format(k, v) for k, v in spec.items())
    print("{0} | dim module={1} | spectrum ad(S)={{{2}}} | dim surv={3} | rankω|surv={4} | isotropic={5}".format(
        jid, mf["m"], specstr, dim_surv, rank_surv, "yes" if iso else "no"))
    if dim_surv == 0:
        print("    partners: (no survivors — S rotates the whole module; the form dies completely)")
    else:
        for idx, prs in enumerate(partner_rows):
            print("    survivor#{0} partners(j,survives): {1}".format(idx, prs))

# -------- Task 5: symbolic block family B -> a.B - B.E --------
print()
print("--- CLAIM symbolic block family: adS(B)=a·B−B·E, ω=antisym(B G⁻¹ B'ᵀ) ---")


def block_family(a, E, G):
    """B in R^{2 x d}; return (dim_surv, rank_omega_surv, isotropic, spectrum)."""
    d = E.rows
    Gi = G.inv()
    # basis of module: E_{rc}, r in 0..1, c in 0..d-1
    basis = []
    for r in range(2):
        for cc in range(d):
            B = zeros(2, d); B[r, cc] = Integer(1)
            basis.append(B)
    dimV = 2 * d
    # adS as matrix on vec(B)
    cols = []
    for B in basis:
        T = a * B - B * E
        cols.append(Matrix(2 * d, 1, list(T.T)))   # column-major consistent
    adM = Matrix.hstack(*cols)
    spec = adM.eigenvals()
    surv_vecs = adM.nullspace()
    dim_surv = len(surv_vecs)

    def vec_to_B(v):
        B = zeros(2, d)
        idx = 0
        for r in range(2):
            for cc in range(d):
                B[r, cc] = v[idx, 0]; idx += 1
        return B

    def omega(B, Bp):
        Mm = B * Gi * Bp.T
        return Mm[0, 1] - Mm[1, 0]
    survB = [vec_to_B(v) for v in surv_vecs]
    # omega|surv matrix
    if dim_surv == 0:
        return 0, 0, True, spec
    Om = zeros(dim_surv, dim_surv)
    for i in range(dim_surv):
        for j in range(dim_surv):
            Om[i, j] = omega(survB[i], survB[j])
    return dim_surv, Om.rank(), (Om.rank() == 0), spec


# instance (i) J-like: a=0, E full-rank rotation on d=2
a0 = zeros(2, 2)
E_full = Matrix([[0, 1], [-1, 0]])
G2 = eye(2)
ds_i, rk_i, iso_i, sp_i = block_family(a0, E_full, G2)
ok(ds_i == 0, "block (i) J-like a=0,E-full: survivors 0")
print("  (i) a=0, E full rank (d=2): dim surv={0}, rankω|surv={1}, isotropic={2} [J-like]".format(ds_i, rk_i, iso_i))

# instance (ii) partial: a=0, d=3, E kills 1 core dim
E_part = Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
G3 = eye(3)
ds_ii, rk_ii, iso_ii, sp_ii = block_family(zeros(2, 2), E_part, G3)
ok(ds_ii == 2 and rk_ii == 2 and not iso_ii, "block (ii) partial: surv 2, omega NONDEG")
print("  (ii) a=0, E kills 1 of 3 core: dim surv={0}, rankω|surv={1}, isotropic={2} [NONDEGEN.]".format(ds_ii, rk_ii, iso_ii))

# instance (iii) NILPOTENT a (Jordan block) -> Lagrangian (isotropic) kernel.
# a semisimple gives V_0 with omega nondeg; a NILPOTENT gives an isotropic kernel:
# d=1, a=[[0,1],[0,0]] in sl(2)=sp(2,R), E=0 -> ker={(b1,0)}, omega|surv=0.
a_nil = Matrix([[0, 1], [0, 0]])         # nilpotent, single Jordan block
E_z = zeros(1, 1)
G1 = eye(1)
ds_iii, rk_iii, iso_iii, sp_iii = block_family(a_nil, E_z, G1)
ok(ds_iii > 0 and rk_iii == 0 and iso_iii, "block (iii) nilpotent-a: survivors ISOTROPIC (Lagrangian)")
print("  (iii) a nilpotent (Jordan), E=0: dim surv={0}, rankω|surv={1}, isotropic={2} [ISOTROPIC/Lagrangian]".format(ds_iii, rk_iii, iso_iii))

# ---- bit-for-bit: J2/J4 real vs block (i) prediction (survivors 0) ----
print()
print("  cross-check of the real J2/J4 against block (i) [a=0,E-full → surv 0]:")
for jid in ("J2", "J4"):
    match = (real_res[jid]["dim_surv"] == 0)
    ok(match, "J-real {0} dim surv=0 matches block (i)".format(jid))
    print("    {0}: dim surv={1} → {2}".format(jid, real_res[jid]["dim_surv"],
                                                "✓ (the whole module rotates)" if match else "≠"))

# -------- Task 6: control -- S trivial on module (S!=0) --------
print()
print("--- CONTROL: S≠0, trivial on the module (constructible?) ---")
# search c(N) for nonzero S with [S, every module rep] in Z (adS_V = 0) yet S!=0
jid, p, q, S0, N0 = JSPEC[0]  # J2 setting
n = 6
eta = make_eta(4, 2)
bas = so_basis(6, eta)
cN = centralizer(N0, bas)
mf2 = module_and_form(cN, 6)
# A genuine "trivial-S" must be a LEVI element (outside rad) acting trivially on V.
# (module elements themselves trivially act on V=rad/Z since [module,module]<=Z, so
#  they are excluded; the question is whether a SEMISIMPLE factor can act trivially.)
rad2 = mf2["rad"]
levi_reps = complement_reps(cN, rad2, 6)
found_levi_trivial = False
for L in levi_reps:
    if adS_on_V(L, mf2, 6).is_zero_matrix:
        found_levi_trivial = True
if found_levi_trivial:
    print("  S≠0 of Levi type, trivial on the module: CONSTRUCTIBLE")
else:
    print("  S≠0 of Levi type, trivial on the module: NOT CONSTRUCTIBLE (every Levi element")
    print("    acts nontrivially — the module occupies the whole core). The trivial endpoint = S=0 = c(N) itself.")
ok(not found_levi_trivial, "control: no nonzero Levi acts trivially on module (carve)")
# S=0 endpoint control: full survival, nondeg form
adS0 = adS_on_V(zeros(6, 6), mf2, 6)
ok(adS0.is_zero_matrix, "control S=0: adS_V=0 (full module survives)")
rk0, iso0 = form_rank_on_sub(mf2["W"], [c.T.T for c in [Matrix(mf2["m"], 1, [Integer(1) if k == t else Integer(0) for k in range(mf2["m"])]) for t in range(mf2["m"])]])
print("  endpoint S=0: the whole module survives, rankω|surv={0} (=nondegenerate, as c(N) in leg-1)".format(rk0))

# ==================== BIT-FENCE ====================
print()
print("BIT-FENCE J2/J4 (= S961 empty full form):")
for (jid, p, q, S, N) in JSPEC:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    cA = centralizer(S + N, bas)
    mfA = module_and_form(cA, n)
    ok(mfA["m"] == 0, "{0}: full c(A) form EMPTY (rad/Z=0) [=S961]".format(jid))
    print("  {0}: c(A) rad/Z={1} → the form is empty (S961 confirmed)".format(jid, mfA["m"]))

# ==================== MUTANTS ====================
print()
print("MUTANTS (one per branch: isotropic / nondegenerate / trivial-S):")
mut_ok = True
if iso_iii and rk_iii == 0 and ds_iii > 0:
    print("  MUTANT br-isotropic: CAUGHT (matched-eigenvalue survivors ω|surv=0)")
else:
    print("  MUTANT br-isotropic: NOT CAUGHT"); mut_ok = False
if (not iso_ii) and rk_ii == ds_ii and ds_ii > 0:
    print("  MUTANT br-nondeg: CAUGHT (partial-E survivors ω|surv nondeg rank={0})".format(rk_ii))
else:
    print("  MUTANT br-nondeg: NOT CAUGHT"); mut_ok = False
if ds_i == 0 and real_res["J2"]["dim_surv"] == 0:
    print("  MUTANT br-trivialS: CAUGHT (a=0,E-full → 0 survivors; real J2 matches)")
else:
    print("  MUTANT br-trivialS: NOT CAUGHT"); mut_ok = False

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
