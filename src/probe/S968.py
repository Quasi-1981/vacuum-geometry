# -*- coding: utf-8 -*-
# DIM: na (constructibility of the form-fate regimes for real A=S+N; deep wedges n=6,7; 0 handles).
#
# ============================================================================
# NOTE: the referenced blind TZ file (TZ_W36_LEG3_..._BLIND.md) was NOT present
# in hub/ (commit 426c84c9 added only a hub/prime/* exante, behind the blindness
# fence).  Executed from the inline task spec verbatim; fence held (no hub/prime/*).
#
# OPERATIONALIZATION (stamped BEFORE any counting; continues S966):
# For a deep rank0 wedge N in so(p,q), c(N) = (so(2,1) + so(G)) |x Heisenberg.
# ad(S) descends to V=rad/Z as a symplectic operator; SURVIVORS = ker(adS_V).
# Leg-2 (S966) showed 3 fate regimes in an ABSTRACT block family; this leg asks
# whether they are REALIZED by actual A=S+N with real S in c(N):
#   (ii) k>0  : a SEMISIMPLE S with partial core action (ker E != 0) -> survivors
#               with omega|surv NONDEGENERATE -> full c(A) carries k = dim(ker E) live.
#   (iii) Lagrangian : a NILPOTENT S in the so(2,1) Levi -> isotropic survivors.
# Candidate S drawn from c(N):
#   * core rotations block_gen(core_i,core_j) (in c(N) since core is disjoint from
#     supp(N)); PARTIAL iff core dim d is odd (single rotation always fixes a line)
#     -> so d=3 (n=7) admits regime (ii); d=2 (n=6) core-rotation is FULL -> k=0.
#   * nilpotent Levi element: integer scan of the Levi complement of rad for a
#     nilpotent matrix with nonzero adS_V (the so(2,1) nilpotent cone).
# Per candidate A=S+N: spectrum ad(S)|V, dim 0-space, omega|surv rank + isotropy,
#   and k(c(A)) := rank(omega of full c(A))/2 (S961 operationalization).
#
# Discipline: 0 handles; exact arithmetic; mutant on each branch {k>0 / Lagrangian /
# full-core-kills}; FORBIDDEN-SCAN (S929); bit-fence vs S961/S966; STOP after table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import itertools
from math import lcm as _ilcm
from sympy import Matrix, Integer, Rational, zeros, eye, diag

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S968_run.log")
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


def is_semisimple(A):
    """diagonalizable over C: minimal poly squarefree."""
    try:
        return A.is_diagonalizable()
    except Exception:
        return False


def block_gen(n, eta, i, j, param):
    return Integer(param) * wedge(unit(n, i), unit(n, j), eta)


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
    reps, adapted, s, m = mf["reps"], mf["adapted"], mf["s"], mf["m"]
    M = zeros(m, m)
    for j in range(m):
        co = coords_in(adapted, S * reps[j] - reps[j] * S, n)
        for i in range(m):
            M[i, j] = co[s + i, 0]
    return M


def form_rank_on_sub(W, basisvecs):
    if not basisvecs:
        return 0, True
    P = Matrix.hstack(*basisvecs)
    t = P.cols
    Wls = [P.T * W[l] * P for l in range(len(W))]
    stack = []
    for b in range(t):
        for l in range(len(W)):
            stack.append([Wls[l][a, b] for a in range(t)])
    Mk = Matrix(stack) if stack else zeros(0, t)
    rank = Mk.rank() if Mk.rows > 0 else 0
    return rank, (rank == 0)


def k_of_full(A, bas, n):
    """k = rank(omega of full c(A))/2 (S961)."""
    mf = module_and_form(centralizer(A, bas), n)
    if mf["m"] == 0:
        return 0, "empty"
    rows = []
    for j in range(mf["m"]):
        for l in range(mf["s"]):
            rows.append([mf["W"][l][i, j] for i in range(mf["m"])])
    Mk = Matrix(rows)
    rk = Mk.rank()
    return rk // 2, ("nondegenerate" if rk == mf["m"] else "degenerate")


def supp_of(N, n):
    s = set()
    for i in range(n):
        for j in range(n):
            if N[i, j] != 0:
                s.add(i); s.add(j)
    return s


# ============================================================================
#                                   RUN
# ============================================================================

print("=" * 74)
print("W36-leg-3: constructibility of the form-fate regimes (real A=S+N)")
print("=" * 74)

# deep rank0 wedges: n=6 (4,2) d=2 core; n=7 (5,2) d=3 core
CASES = [("n6", 4, 2), ("n7", 5, 2)]

found_k_pos = False       # target: any A with k>0 (regime ii)
found_lagrangian = False  # target: any Lagrangian (regime iii)

for (tag, p, q) in CASES:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    xy = d0_on(list(range(n)), p, n)
    x, y = xy
    N = wedge(x, y, eta)
    ok(is_nilp(N), "N nilpotent " + tag)
    c = centralizer(N, bas)
    mf = module_and_form(c, n)
    d_core = n - 4
    supp = supp_of(N, n)
    core_idx = [i for i in range(n) if i not in supp]
    print()
    print("--- {0}: sig ({1},{2}) deep rank0 wedge, core dim d={3}, core-indices={4} ---".format(
        tag, p, q, d_core, core_idx))
    print("    c(N): dim of module V={0}, dim Z={1}".format(mf["m"], mf["s"]))

    candidates = []
    # regime-(ii) candidates: core rotations (partial iff d odd)
    for (i, j) in itertools.combinations(core_idx, 2):
        S = block_gen(n, eta, i, j, 1)
        if (S * N - N * S).is_zero_matrix:
            candidates.append(("semisimple core-rot({0},{1})".format(i, j), S))

    # regime-(iii) candidate: scan Levi complement of rad for a nilpotent with adS_V != 0
    rad = mf["rad"]
    levi = complement_reps(c, rad, n)
    nil_found = None
    for combo in itertools.product((-1, 0, 1), repeat=len(levi)):
        if all(v == 0 for v in combo):
            continue
        S = zeros(n, n)
        for coef, L in zip(combo, levi):
            if coef:
                S = S + Integer(coef) * L
        if is_nilp(S) and not adS_on_V(S, mf, n).is_zero_matrix:
            nil_found = S
            break
    if nil_found is not None:
        candidates.append(("nilpotent-Levi", nil_found))
    else:
        print("    (a nilpotent Levi with nontrivial action: not found in the {−1,0,1}-scan)")

    # analyze each candidate
    for (name, S) in candidates:
        adSV = adS_on_V(S, mf, n)
        spec = adSV.eigenvals()
        ker0 = adSV.nullspace()
        dim0 = len(ker0)
        rk_surv, iso = form_rank_on_sub(mf["W"], ker0)
        A = S + N
        kfull, kstat = k_of_full(A, bas, n)
        nilp = is_nilp(S)
        ss = is_semisimple(S)
        typ = "nilpotent" if nilp else ("semisimple" if ss else "mixed")
        specstr = ", ".join("{0}^{1}".format(kk, vv) for kk, vv in spec.items())
        if kfull > 0:
            found_k_pos = True
        if dim0 > 0 and rk_surv == 0 and not nilp is False and iso and nilp:
            found_lagrangian = True
        # also flag Lagrangian purely by: nonzero 0-space that is isotropic AND from nilpotent S
        lag = (dim0 > 0 and iso and nilp)
        if lag:
            found_lagrangian = True
        print("    A=S+N [{0}, {1}] | spectrum ad(S)|V={{{2}}} | dim 0-space={3} | ω|surv: rank={4} isotropic={5} | k(c(A))={6} [{7}]".format(
            name, typ, specstr, dim0, rk_surv, "yes" if iso else "no", kfull, kstat))

# ==================== RAW TARGET ====================
print()
print("RAW TARGET:")
print("  regime (ii) k>0 is constructible: {0}".format("YES" if found_k_pos else "NO"))
print("  regime (iii) Lagrangian is constructible: {0}".format("YES" if found_lagrangian else "NO"))

# ==================== BIT-FENCE ====================
print()
print("BIT-FENCE (S961/S966):")
# n=6 core rotation is FULL (d=2 even) -> k=0, like J2 (S966)
_e = make_eta(4, 2); _b = so_basis(6, _e)
_x, _y = d0_on(list(range(6)), 4, 6)
_N = wedge(_x, _y, _e)
_supp = supp_of(_N, 6)
_core = [i for i in range(6) if i not in _supp]
_S = block_gen(6, _e, _core[0], _core[1], 1)
_k, _st = k_of_full(_S + _N, _b, 6)
ok(_k == 0, "n=6 full core-rot: k=0 (=J2 S966, core dimension even)")
print("  n=6 core={0} full-rotation: k={1} (even core → full action → the form dies, =S966 J2)".format(_core, _k))
# c(N) alone: nondeg form (S961)
_mf = module_and_form(centralizer(_N, _b), 6)
_rows = [[_mf["W"][l][i, j] for i in range(_mf["m"])] for j in range(_mf["m"]) for l in range(_mf["s"])]
ok(Matrix(_rows).rank() == _mf["m"], "n=6 c(N) form nondeg (S961)")
print("  n=6 c(N) alone: ω nondegenerate rank={0} (S961)".format(_mf["m"]))

# ==================== MUTANTS ====================
print()
print("MUTANTS (by branch: k>0 / Lagrangian / full-core-kills):")
mut_ok = True

# br-kpos: n=7 partial core rotation MUST give k>0
_e7 = make_eta(5, 2); _b7 = so_basis(7, _e7)
_x7, _y7 = d0_on(list(range(7)), 5, 7)
_N7 = wedge(_x7, _y7, _e7)
_core7 = [i for i in range(7) if i not in supp_of(_N7, 7)]
_S7 = block_gen(7, _e7, _core7[0], _core7[1], 1)   # rotates 2 of 3 core -> partial
_k7, _ = k_of_full(_S7 + _N7, _b7, 7)
if _k7 > 0:
    print("  MUTANT br-kpos: CAUGHT (n=7 partial core-rot → k={0}>0, regime ii is alive)".format(_k7))
else:
    print("  MUTANT br-kpos: NOT CAUGHT"); mut_ok = False

# br-fullkill: n=6 full core rotation MUST give k=0
if _k == 0:
    print("  MUTANT br-fullkill: CAUGHT (n=6 full core-rot → k=0, the form is dead)")
else:
    print("  MUTANT br-fullkill: NOT CAUGHT"); mut_ok = False

# br-lagrangian: a nilpotent S with nonzero isotropic 0-space (if constructed)
if found_lagrangian:
    print("  MUTANT br-lagrangian: CAUGHT (a nilpotent S gave an isotropic 0-space)")
else:
    print("  MUTANT br-lagrangian: NOT CAUGHT (a Lagrangian was not constructed in the scan)")
    mut_ok = False

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
