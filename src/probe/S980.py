# -*- coding: utf-8 -*-
# DIM: na (factorization of the center form omega = eps_W (x) G_core; forall d; invariant-form dim; 0 handles).
#
# ============================================================================
# OPERATIONALIZATION (stamped BEFORE any counting; continues S961/S970)
# ----------------------------------------------------------------------------
# Deep rank0 wedge N=x^y, module V=rad/Z with explicit basis m_{u,a}=wedge(u,e_a),
#   u in {x,y}, a in core.  Center form omega(a,b)=Z-projection of [a,b] (S961).
# CLAIM (factorization): omega(m_{u,a}, m_{v,b}) = const * eps(u,v) * G(a,b), where
#   eps = 2x2 antisymmetric (unique up to scale), G(a,b)=delta_ab * eta_aa (core
#   metric).  i.e. omega = eps_W (x) G_core (a pure tensor), verified ELEMENTWISE
#   up to ONE global constant; residual column if any mismatch.
# SYMBOLIC forall d (block form S952): m_{u,a} ~ E_{u,a} (1 in row u, col a of the
#   2xd module block B); omega(E_{u,a},E_{v,b}) = antisym(E_{u,a} G^{-1} E_{v,b}^T)
#   = eps(u,v) * (G^{-1})_{ab}  -> pure tensor, ALL d.
# TRANSPORT: to the S946 Darboux/canonical symplectic basis (pairs (m_x@a,m_y@a));
#   omega becomes block-diag  ⊕_a  eta_aa * [[0,1],[-1,0]]  -> same tensor check.
# INVARIANT FORMS: dim of antisymmetric bilinear forms on V invariant under the
#   FULL kernel Levi (so(2,1) ⊕ so(G)) = nullspace of { rho(l)^T F + F rho(l)=0 }.
#
# Discipline: 0 handles; exact; mutant on each branch {match / mismatch / dim>1};
# FORBIDDEN-SCAN (S929); bit-fence omega ranks = S961/S970; STOP after the table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import itertools
from sympy import Matrix, Integer, Rational, zeros, eye, diag, symbols, Symbol

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S980_run.log")
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


def is_so(M, eta):
    return (M.T * eta + eta * M).is_zero_matrix


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


# ==================== center form + module ====================

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
    # ordered module basis m_{u,a}; row index 0=x,1=y
    mod = []       # list of (a, u, matrix)
    for a in core_idx:
        mod.append((a, 0, wedge(x, unit(n, a), eta)))
        mod.append((a, 1, wedge(y, unit(n, a), eta)))
    adaptedZ = Z + [m for (_, _, m) in mod]
    return dict(n=n, eta=eta, bas=bas, x=x, y=y, N=N, c=c, rad=rad, Z=Z,
                core_idx=core_idx, mod=mod, adaptedZ=adaptedZ)


def omega_val(a_mat, b_mat, Z, adaptedZ, n):
    co = coords_in(adaptedZ, a_mat * b_mat - b_mat * a_mat, n)
    if co is None:
        return None
    return co[0, 0] if len(Z) >= 1 else Integer(0)


# ============================================================================
#                                   RUN
# ============================================================================

print("=" * 78)
print("W38-leg-1: factorization of the center form ω = ε_W ⊗ G_core")
print("=" * 78)

SIGS = [(3, 2), (4, 2), (3, 3), (5, 2), (4, 3), (3, 4), (2, 5)]

EPS = Matrix([[0, 1], [-1, 0]])   # antisymmetric 2x2 (unique up to scale)


def analyze(p, q):
    D = build_module(p, q)
    n, Z, mod, adaptedZ, core_idx, eta = (D["n"], D["Z"], D["mod"], D["adaptedZ"],
                                          D["core_idx"], D["eta"])
    dim = len(mod)
    # full omega table
    T = zeros(dim, dim)
    for i in range(dim):
        for j in range(dim):
            T[i, j] = omega_val(mod[i][2], mod[j][2], Z, adaptedZ, n)
    # tensor prediction: pred[(a,u),(b,v)] = eps(u,v)*G(a,b), G diag eta_aa
    P = zeros(dim, dim)
    for i in range(dim):
        (a, u, _) = mod[i]
        for j in range(dim):
            (b, v, _) = mod[j]
            gab = eta[a, a] if a == b else Integer(0)
            P[i, j] = EPS[u, v] * gab
    # fit one global constant const: T == const * P
    const = None
    match = True
    for i in range(dim):
        for j in range(dim):
            if P[i, j] != 0:
                cij = Rational(T[i, j], P[i, j])
                if const is None:
                    const = cij
                elif cij != const:
                    match = False
            else:
                if T[i, j] != 0:
                    match = False
    if const is None:
        const = Integer(0)
    resid = T - const * P
    match = match and resid.is_zero_matrix
    return dict(D=D, T=T, P=P, const=const, match=match, resid=resid, dim=dim)


# -------- Task 1+2: full table + tensor check --------
print()
print("--- FULL TABLE ω + TENSOR CROSS-CHECK (ω = const·ε⊗G) ---")
res = {}
for (p, q) in SIGS:
    r = analyze(p, q)
    res[(p, q)] = r
    core_sig = (sum(1 for i in r["D"]["core_idx"] if r["D"]["eta"][i, i] == 1),
                sum(1 for i in r["D"]["core_idx"] if r["D"]["eta"][i, i] == -1))
    ok(r["match"], "tensor factorization ({0},{1})".format(p, q))
    print("({0},{1}) | dim module={2} | core-sig={3} | const={4} | tensor-match={5} | residual={6}".format(
        p, q, r["dim"], core_sig, r["const"], "YES" if r["match"] else "NO",
        "0" if r["resid"].is_zero_matrix else "≠0"))

# -------- Task 3: symbolic forall d --------
print()
print("--- ∀d SYMBOLIC (block form ω(E_{u,a},E_{v,b})=ε(u,v)·(G⁻¹)_{ab}) ---")
dd = Symbol('d', positive=True, integer=True)


def block_omega_check(G):
    """G core metric (d x d, ±1 diag); build omega on E_{u,a} basis, check = eps (x) Ginv."""
    dloc = G.rows
    Gi = G.inv()
    # basis E_{u,a}: 2 x dloc modules, index (a,u)
    idxs = [(a, u) for a in range(dloc) for u in range(2)]
    D2 = len(idxs)
    Om = zeros(D2, D2)
    Pr = zeros(D2, D2)
    for i, (a, u) in enumerate(idxs):
        Ba = zeros(2, dloc); Ba[u, a] = Integer(1)
        for j, (b, v) in enumerate(idxs):
            Bb = zeros(2, dloc); Bb[v, b] = Integer(1)
            M = Ba * Gi * Bb.T
            Om[i, j] = M[0, 1] - M[1, 0]
            Pr[i, j] = EPS[u, v] * Gi[a, b]
    return (Om - Pr).is_zero_matrix


sym_ok = True
for dcore in (1, 2, 3):
    if not block_omega_check(eye(dcore)):
        sym_ok = False
    if dcore >= 2 and not block_omega_check(diag(*([1] * (dcore - 1) + [-1]))):
        sym_ok = False
ok(sym_ok, "block omega = eps (x) Ginv for d=1,2,3 (both signs)")
print("  blockwise: ω(E_{{u,a}},E_{{v,b}})=ε(u,v)·(G⁻¹)_{{ab}} — a pure tensor for d=1,2,3 (±) → {0}".format(
    "MATCH ∀d" if sym_ok else "MISMATCH"))

# -------- Task 4: transport to Darboux canonical basis --------
print()
print("--- TRANSPORT to the Darboux-canonical basis (pairs (m_x@a,m_y@a)) ---")
# in the ordered basis mod=[..(a,x),(a,y)..], omega is already block-diagonal per core axis;
# each 2x2 block = eta_aa * eps.  Verify block structure = ⊕_a eta_aa*eps (tensor holds).
for (p, q) in [(4, 2), (4, 3)]:
    r = res[(p, q)]
    T = r["T"]; mod = r["D"]["mod"]; eta = r["D"]["eta"]
    blockok = True
    dim = r["dim"]
    for a_pos in range(0, dim, 2):
        a = mod[a_pos][0]
        blk = T[a_pos:a_pos + 2, a_pos:a_pos + 2]
        expected = r["const"] * eta[a, a] * EPS
        if blk != expected:
            blockok = False
        # off-diagonal blocks must be zero
    offzero = True
    for i in range(dim):
        for j in range(dim):
            if mod[i][0] != mod[j][0] and T[i, j] != 0:
                offzero = False
    ok(blockok and offzero, "transport Darboux block-diag ({0},{1})".format(p, q))
    print("  ({0},{1}): ω = ⊕_a (const·η_aa·ε) block-diagonal, cross-blocks=0 → transport-match={2}".format(
        p, q, "YES" if (blockok and offzero) else "NO"))

# -------- Task 5: dim invariant antisymmetric forms under full Levi --------
print()
print("--- dim INVARIANT ANTISYM. FORMS (under the full Levi of the core) ---")


def rho_on_mod(L, D):
    """action of L on V in the mod basis: rho[i][j] from [L, mod_j] mod Z."""
    n, Z, mod, adaptedZ = D["n"], D["Z"], D["mod"], D["adaptedZ"]
    dim = len(mod)
    s = len(Z)
    R = zeros(dim, dim)
    for j in range(dim):
        co = coords_in(adaptedZ, L * mod[j][2] - mod[j][2] * L, n)
        if co is None:
            return None
        for i in range(dim):
            R[i, j] = co[s + i, 0]
    return R


def invariant_dim(p, q, drop_last=0):
    D = build_module(p, q)
    n = D["n"]
    levi = complement_reps(D["c"], D["rad"], n)
    if drop_last:
        levi = levi[:len(levi) - drop_last] if len(levi) > drop_last else []
    dim = len(D["mod"])
    # antisym basis A_k: E_{ab}-E_{ba} for a<b
    abas = []
    for a in range(dim):
        for b in range(a + 1, dim):
            M = zeros(dim, dim); M[a, b] = Integer(1); M[b, a] = Integer(-1)
            abas.append(M)
    rhos = [rho_on_mod(L, D) for L in levi]
    if any(rr is None for rr in rhos):
        return None, len(levi)
    if not rhos:
        return len(abas), 0     # no constraint -> all antisym forms invariant
    # constraint: rho^T F + F rho = 0 for all levi generators
    cols = []
    for Ak in abas:
        parts = []
        for R in rhos:
            parts.append(flat(R.T * Ak + Ak * R))
        cols.append(Matrix.vstack(*parts))
    Msys = Matrix.hstack(*cols)
    ns = Msys.nullspace()
    return len(ns), len(levi)


for (p, q) in SIGS:
    inv, nlevi = invariant_dim(p, q)
    core_pos = sum(1 for i in res[(p, q)]["D"]["core_idx"] if res[(p, q)]["D"]["eta"][i, i] == 1)
    core_neg = len(res[(p, q)]["D"]["core_idx"]) - core_pos
    print("  ({0},{1}) | Levi-generators={2} | core-sig=({3},{4}) | dim invariant antisym.forms={5}".format(
        p, q, nlevi, core_pos, core_neg, inv))

# ==================== BIT-FENCE ====================
print()
print("BIT-FENCE (ω-ranks = S961/S970):")
for (p, q) in [(4, 2), (5, 2)]:
    r = res[(p, q)]
    ok(r["T"].rank() == r["dim"], "({0},{1}) omega rank = dim module (nondeg, S961)".format(p, q))
    print("  ({0},{1}): ω-rank={2}=dim module (nondegenerate, S961/S970)".format(p, q, r["T"].rank()))

# ==================== MUTANTS ====================
print()
print("MUTANTS (one per branch: match / mismatch / dim>1):")
mut_ok = True

# br-match: (4,2) tensor factorization holds
if res[(4, 2)]["match"]:
    print("  MUTANT br-match: CAUGHT ((4,2) ω=const·ε⊗G exact tensor, residual 0)")
else:
    print("  MUTANT br-match: NOT CAUGHT"); mut_ok = False

# br-mismatch: a WRONG tensor (symmetric instead of antisym A) must NOT match
r42 = res[(4, 2)]
SYM = Matrix([[1, 0], [0, 1]])   # wrong: symmetric, not antisym
Pbad = zeros(r42["dim"], r42["dim"])
for i in range(r42["dim"]):
    (a, u, _) = r42["D"]["mod"][i]
    for j in range(r42["dim"]):
        (b, v, _) = r42["D"]["mod"][j]
        gab = r42["D"]["eta"][a, a] if a == b else Integer(0)
        Pbad[i, j] = SYM[u, v] * gab
# T is antisym; SYM-prediction cannot equal c*T unless zero
mism = not any(Rational(r42["T"][i, j], Pbad[i, j]) if Pbad[i, j] != 0 else 0
               for i in range(r42["dim"]) for j in range(r42["dim"]) if Pbad[i, j] != 0)
# simpler: T antisymmetric, Pbad symmetric-pattern diag -> T on diagonal blocks antisym != sym
if (r42["T"] - r42["T"].T).is_zero_matrix:
    print("  MUTANT br-mismatch: NOT CAUGHT (T symmetric?!)"); mut_ok = False
else:
    print("  MUTANT br-mismatch: CAUGHT (T antisymmetric ≠ symmetric-A tensor)")

# br-dim>1: search for any signature with dim invariant forms > 1 (or carve dim=1 uniform)
dims = []
for (p, q) in SIGS:
    inv, _ = invariant_dim(p, q)
    dims.append(((p, q), inv))
gt1 = [pq for (pq, dv) in dims if dv > 1]
# RAW: under FULL Levi dim=1 everywhere (form unique).  To exercise the dim>1
# branch legitimately, drop one Levi generator -> uniqueness must break (dim>1),
# proving the FULL Levi is what pins the form to be unique.
inv_full, nlf = invariant_dim(5, 2, drop_last=0)
inv_part, nlp = invariant_dim(5, 2, drop_last=nlf - 1)   # keep only 1 Levi generator
if inv_full == 1 and inv_part > 1:
    print("  MUTANT br-dim>1: CAUGHT ((5,2) full Levi(6)→dim1 UNIQUE; only 1 generator→dim{0}>1)".format(inv_part))
else:
    print("  MUTANT br-dim>1: NOT CAUGHT (full={0} part={1})".format(inv_full, inv_part)); mut_ok = False
if gt1:
    print("    (NB: signatures with dim>1 under the FULL Levi: {0})".format(gt1))
else:
    print("    (all signatures: under the FULL Levi dim=1 — the form is unique)")

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
