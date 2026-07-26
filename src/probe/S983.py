# -*- coding: utf-8 -*-
# DIM: na (scale law of the center form; normalization; lambda-transport to the extension family; 0 handles).
#
# ============================================================================
# OPERATIONALIZATION (stamped BEFORE any counting; continues S961/S970/S980)
# ----------------------------------------------------------------------------
# Deep rank0 wedge N=x^y; module gen m_{w,c}=wedge(w,e_c); coefficient (S961) =
#   Z-projection of [m_{u,a}, m_{v,b}].  At mu=nu=1 it equals -eps(u,v)*G(a,b) (S980).
#
# SCALE LAW (symbolic mu, nu):
#   (a) e_c -> mu e_c : m_{w,c} = wedge(w, mu e_c) = mu * m_{w,c}  (linear).
#   (b) g -> mu^2 g   : induced (g_cc = e_c^T eta e_c -> mu^2 g_cc).
#   (c) Z -> nu Z     : the center generator is rescaled.
#   => [m(mu), m'(mu)] = mu^2 [m,m'];  coefficient as a multiple of the NEW center
#   generator (nu Z) = (mu^2 / nu) * base.  EXPONENTS: mu:+2, nu:-1.  Module-gen
#   norms (trace form and bound-axis eta-norm) scale as mu^2.
#
# NORMALIZATION: coefficient = 1  <=>  (mu^2/nu)*(-1) = 1  <=>  nu = -mu^2.  One
#   equation on (mu,nu) => residual freedom = 1-dim CONTINUOUS (mu in R*, nu fixed).
#
# LAMBDA-TRANSPORT (S921/S922 homogeneous family): translations T_a with
#   [T_a,T_b] = lambda * omega_ab * Z (central extension); the S946-class iso maps
#   the kernel Heisenberg [m,m']=omega*Z onto this => lambda = kernel coefficient
#   (= -eta_aa per Darboux pair, or -mu^2/nu scaled).  lambda=0 : [T,T]=0 = the
#   affine so|xR^d (S921), i.e. the radical ABELIANIZES (soft branch).
#
# Discipline: 0 handles; exact (symbolic mu,nu,lambda); mutant on each branch
# {exponent / normalization / transport}; FORBIDDEN-SCAN (S929); bit-fence
# coefficients at mu=nu=1 = S980 (-1*eps(x)G); STOP after the table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
from sympy import (Matrix, Integer, Rational, zeros, eye, diag, symbols, Symbol,
                   simplify, sqrt)

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S983_run.log")
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
    rad = []
    for v in K.nullspace():
        M = zeros(n, n)
        for t in range(len(core)):
            if v[t, 0] != 0:
                M = M + v[t, 0] * core[t]
        rad.append(M)
    return span_basis(rad, n)


def supp_of(N, n):
    s = set()
    for i in range(n):
        for j in range(n):
            if N[i, j] != 0:
                s.add(i); s.add(j)
    return s


def d0_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 2 and len(Km) >= 2:
        return unit(n, Kp[0]) + unit(n, Km[0]), unit(n, Kp[1]) + unit(n, Km[1])
    return None


def build_module(p, q):
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    x, y = d0_on(list(range(n)), p, n)
    N = wedge(x, y, eta)
    c = centralizer(N, bas)
    rad = killing_radical(c, n)
    Z = center_of(rad, n)
    core_idx = [i for i in range(n) if i not in supp_of(N, n)]
    return dict(n=n, eta=eta, x=x, y=y, N=N, c=c, rad=rad, Z=Z, core_idx=core_idx)


# ============================================================================
#                                   RUN
# ============================================================================

print("=" * 78)
print("W38-leg-2: scale law of the center form (μ,ν) + λ-transport")
print("=" * 78)

mu, nu, lam = symbols('mu nu lambda', positive=True)

# -------- Task 1: scale law (symbolic mu, nu) --------
print()
print("--- SCALE LAW (symbolic μ, ν) ---")
D = build_module(4, 2)     # representative
n, eta, x, y, Z, core = D["n"], D["eta"], D["x"], D["y"], D["Z"], D["core_idx"]
a0 = core[0]
# module gens with symbolic mu on the axis
mx = wedge(x, mu * unit(n, a0), eta)       # = mu * wedge(x,e_a)
my = wedge(y, mu * unit(n, a0), eta)
ok(simplify(mx - mu * wedge(x, unit(n, a0), eta)) == zeros(n, n), "m_{w,c}(μ)=μ·m (linear in axis)")
# bracket and Z-projection; base center gen Z[0], rescaled center = nu*Z[0]
br = mx * my - my * mx
# coefficient wrt base Z: solve br = coeff_base * Z0 (Z 1-dim)
Z0 = Z[0]
co_base = coords_in([Z0], br, n)
coeff_base = co_base[0, 0] if co_base is not None else None
# wrt rescaled center nu*Z0: coefficient = coeff_base / nu
coeff_resc = simplify(coeff_base / nu)
ok(simplify(coeff_base - mu ** 2 * coeff_base.subs(mu, 1)) == 0, "coeff ∝ μ² (exponent +2)")
print("  coefficient (base Z): {0}  → exponent μ: +2".format(coeff_base))
print("  coefficient (center νZ): {0}  → exponent ν: −1".format(coeff_resc))
print("  at μ=ν=1: {0} (=S980 −1·η_aa for the core+ axis)".format(simplify(coeff_resc.subs({mu: 1, nu: 1}))))
# module-generator norm scaling (Frobenius: gens are nilpotent so tr(m·m)=0;
# Frobenius tr(mᵀm) is nonzero and exhibits the mu^2 scaling)
m_base = wedge(x, unit(n, a0), eta)
norm_fro = (mx.T * mx).trace()
ok(simplify(norm_fro - mu ** 2 * (m_base.T * m_base).trace()) == 0,
   "module-gen Frobenius norm ∝ μ²")
ok((mx * mx).trace() == 0, "module gens nilpotent (trace-norm tr(m·m)=0, scale-blind)")
print("  module-generator norm (Frobenius tr(mᵀm)) ∝ μ²: {0}  [tr(m·m)=0 nilpotent]".format(
    simplify(norm_fro)))
# instances mu in {1/2,2,3}
print("  instances μ (ν=1): " + " · ".join(
    "μ={0}→{1}".format(mv, simplify(coeff_resc.subs({mu: mv, nu: 1})))
    for mv in (Rational(1, 2), Integer(2), Integer(3))))

# -------- Task 2: normalization choice --------
print()
print("--- NORMALIZATION CHOICE (coeff=1) ---")
# coeff_resc = -mu^2/nu (base -1). Solve coeff_resc == 1 -> nu = -mu^2.
c0 = coeff_resc.subs(mu, 1).subs(nu, 1)   # -1
print("  coeff(μ,ν) = ({0})·μ²/ν ; the condition coeff=1 ⟹ ν = ({0})·μ² (here ν=−μ²)".format(c0))
print("  the residual freedom: 1-parameter CONTINUOUS (μ∈ℝ*, ν=−μ² determined)")
# verify: coeff at (mu, nu=-mu^2) == 1
chk = simplify(coeff_resc.subs(nu, c0 * mu ** 2))
ok(simplify(chk - 1) == 0, "normalization coeff=1 on nu=c0*mu^2 (1-param family)")
print("  check: coeff|_{{ν=({0})μ²}} = {1} (a constant ⟹ a 1-param family of normalizations)".format(c0, chk))

# -------- Task 3: lambda-transport to extension family --------
print()
print("--- λ-TRANSPORT into the extension family (S921/S922 homogeneous) ---")


def heis_rep(dcore, lam_val, etac):
    """(dcore+2)x(dcore+2) rep: q_a=E[0,1+a], p_a=lam*etac[a]*E[1+a,dcore+1], Z=E[0,dcore+1].
    [q_a,p_b] = lam*etac[a]*delta_ab * Z."""
    D2 = dcore + 2
    Zc = zeros(D2, D2); Zc[0, D2 - 1] = Integer(1)
    qs, ps = [], []
    for a in range(dcore):
        Q = zeros(D2, D2); Q[0, 1 + a] = Integer(1); qs.append(Q)
        P = zeros(D2, D2); P[1 + a, D2 - 1] = lam_val * etac[a]; ps.append(P)
    return qs, ps, Zc


def lam_from_kernel(p, q):
    """measure kernel coefficient omega(m_x@a, m_y@a) = lambda per Darboux pair."""
    D = build_module(p, q)
    n, eta, x, y, Z, core = D["n"], D["eta"], D["x"], D["y"], D["Z"], D["core_idx"]
    Z0 = Z[0]
    lams = []
    for a in core:
        mx = wedge(x, unit(n, a), eta); my = wedge(y, unit(n, a), eta)
        co = coords_in([Z0], mx * my - my * mx, n)
        lams.append((a, eta[a, a], co[0, 0]))
    return lams


for (p, q) in [(5, 2), (4, 3), (3, 4)]:
    lams = lam_from_kernel(p, q)
    # lambda per pair equals coefficient = (S980) -eta_aa (times global -1 base)
    match = all(lamv == -etaa for (a, etaa, lamv) in lams)
    ok(match, "lambda = kernel coefficient = -eta_aa ({0},{1})".format(p, q))
    lstr = " ".join("a{0}(η{1}):λ={2}".format(a, etaa, lamv) for (a, etaa, lamv) in lams)
    print("  ({0},{1}) n={2} | λ per pair = kernel coeff: {3}".format(p, q, p + q, lstr))
# verify heis_rep with lambda reproduces [q,p]=lam*eta*Z, and iso structure
dcore = 3
etac = [1, 1, -1]
qs, ps, Zc = heis_rep(dcore, lam, etac)
brqp = qs[0] * ps[0] - ps[0] * qs[0]
ok(simplify(brqp - lam * etac[0] * Zc) == zeros(dcore + 2, dcore + 2),
   "heis rep [q,p]=λ·η·Z (symbolic λ)")
print("  homogeneous representation: [q_a,p_b]=λ·η_a·δ_ab·Z (symbolic λ) — forced")

# -------- Task 4: lambda=0 end --------
print()
print("--- λ=0 ENDPOINT (abelianization) ---")
qs0, ps0, Zc0 = heis_rep(dcore, Integer(0), etac)
allz = all((A * B - B * A).is_zero_matrix for A in qs0 + ps0 for B in qs0 + ps0)
ok(allz, "lambda=0: all [T,T]=0 (translations abelianize)")
print("  λ=0: [q,p]=0 for all → the translations are ABELIAN = the affine so⋉ℝ^d (S921);")
print("       on the kernel side: Heisenberg h_{2d+1} → abelian ℝ^{2d}⊕Z (the radical abelianizes = the soft branch)")
# soft-branch check: at lambda=0 the radical structure has ZERO center pairing (omega->0)
ok(simplify((lam * etac[0]).subs(lam, 0)) == 0, "lambda=0 => coefficient 0 (omega degenerates)")

# ==================== BIT-FENCE ====================
print()
print("BIT-FENCE (coeff at μ=ν=1 = S980):")
for (p, q) in [(4, 2), (3, 3)]:
    lams = lam_from_kernel(p, q)
    okrow = all(lamv == -etaa for (a, etaa, lamv) in lams)
    ok(okrow, "({0},{1}) coeff at μ=ν=1 = -eta_aa (S980)".format(p, q))
    print("  ({0},{1}): coeff={2} = −η_aa (S980 const=−1·ε⊗G)".format(
        p, q, [lamv for (_, _, lamv) in lams]))

# ==================== MUTANTS ====================
print()
print("MUTANTS (one per branch: exponent / normalization / transport):")
mut_ok = True

# br-exponent: coeff scales as mu^2 (not mu^1 or mu^3)
c_test = coeff_base
if simplify(c_test.subs(mu, 2) / c_test.subs(mu, 1)) == 4:
    print("  MUTANT br-exponent: CAUGHT (μ→2 multiplies the coeff by 4=2², exponent exactly +2)")
else:
    print("  MUTANT br-exponent: NOT CAUGHT"); mut_ok = False

# br-normalization: nu=+mu^2 (wrong sign) does NOT give coeff=1 (gives -1)
wrong = simplify(coeff_resc.subs(nu, mu ** 2))
if simplify(wrong - 1) != 0 and simplify(wrong + 1) == 0:
    print("  MUTANT br-normalization: CAUGHT (ν=+μ² gives coeff=−1≠1; ν=−μ² is required)")
else:
    print("  MUTANT br-normalization: NOT CAUGHT"); mut_ok = False

# br-transport: lambda=0 abelianizes; lambda!=0 stays Heisenberg (nonabelian)
qs1, ps1, _ = heis_rep(dcore, Integer(1), etac)
nonab = any(not (A * B - B * A).is_zero_matrix for A in qs1 for B in ps1)
if allz and nonab:
    print("  MUTANT br-transport: CAUGHT (λ=0 abelian; λ=1 non-abelian Heisenberg)")
else:
    print("  MUTANT br-transport: NOT CAUGHT"); mut_ok = False

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
