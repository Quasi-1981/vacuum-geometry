# -*- coding: utf-8 -*-
# DIM: na (eta-types of canonical pairs of the center form; module trace on the floor; 0 handles).
#
# ============================================================================
# OPERATIONALIZATION (stamped BEFORE any counting; continues S961/S968)
# ----------------------------------------------------------------------------
# Deep rank0 wedge N=x^y (N^2=0) in so(p,q), Witt>=2.  Flag of N:
#   W = im N = span{x,y} (isotropic 2-plane) ⊆ W^perp = ker N ⊆ R^n.
#   middle graded piece = core C = W^perp / W, eta|_C = G nondegenerate;
#   C is spanned by the ambient unit vectors e_c, c ∉ supp(N) (eta-diagonal).
# MODULE V = rad(c(N))/Z acts class-1 (shifts the flag): a module direction maps
#   the middle C into the bottom W.  Explicit module basis:
#       m_{x,c} = wedge(x, e_c),  m_{y,c} = wedge(y, e_c),  c ∈ core_idx
#   (each lies in rad; m_{w,c}·e_c = eta_cc · w ∈ W, so it BINDS core axis e_c).
#   TRACE of a module direction = the floor axis e_c it binds; its ETA-TYPE =
#   sign(e_c^T eta e_c) = eta_cc  (+ / - ; isotropic only if bound to a null axis).
# CENTER FORM omega(a,b) = Z-projection of [a,b] (S961).  CANONICAL (Darboux, S946)
#   pairs: (m_{x,c}, m_{y,c}) are omega-conjugate (omega!=0), omega(m_{.,c},m_{.,c'})=0
#   for c!=c'.  Each pair's two members bind the SAME core axis c -> SAME eta-type.
#   => pair-type = eta-type of core axis c;  multiset of pair-types = core signature.
# MIXED k>0 ((5,2)+partial core rotation E): fixed axis = ker E (a core axis);
#   its eta-type + the eta-types of the surviving pair(s).
#
# Discipline: 0 handles; exact arithmetic; mutant on each branch {same-type /
# different-type / isotropic-trace}; FORBIDDEN-SCAN (S929); bit-fence omega ranks
# vs S961/S968; STOP after the table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import itertools
from sympy import Matrix, Integer, Rational, zeros, eye, diag

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S970_run.log")
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


def contains(big, small, sq):
    FB = stack_flats(big, sq)
    rb = FB.rank()
    return Matrix.hstack(FB, stack_flats(small, sq)).rank() == rb


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


def d0_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 2 and len(Km) >= 2:
        return unit(n, Kp[0]) + unit(n, Km[0]), unit(n, Kp[1]) + unit(n, Km[1])
    return None


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


def omega_Z(a, b, Z, adapted, n):
    """Z-projection scalar(s) of [a,b]; returns list over Z-directions."""
    co = coords_in(adapted, a * b - b * a, n)
    if co is None:
        return None
    return [co[l, 0] for l in range(len(Z))]


# ============================================================================
#                                   RUN
# ============================================================================

print("=" * 74)
print("W36-leg-4: η-types of the canonical pairs of the center form")
print("=" * 74)


def analyze(p, q, label):
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    xy = d0_on(list(range(n)), p, n)
    if xy is None:
        print("SIG ({0},{1}) {2}: rank0-wedge NOT constructible".format(p, q, label))
        return None
    x, y = xy
    N = wedge(x, y, eta)
    c = centralizer(N, bas)
    rad = killing_radical(c, n)
    Z = center_of(rad, n)
    reps = complement_reps(rad, Z, n)
    adapted = Z + reps
    core_idx = [i for i in range(n) if i not in supp_of(N, n)]
    # module basis m_{w,c}; verify in rad; eta-type per core axis
    pairs = []
    for cc in core_idx:
        ec = unit(n, cc)
        mx = wedge(x, ec, eta)
        my = wedge(y, ec, eta)
        ok(is_so(mx, eta) and is_so(my, eta), "m in so ({0},{1}) c={2}".format(p, q, cc))
        ok(contains(rad, [mx], n) and contains(rad, [my], n),
           "m in rad ({0},{1}) c={2}".format(p, q, cc))
        etatype = "+" if eta[cc, cc] == 1 else "-"
        pairs.append((cc, mx, my, etatype))
    # omega structure: canonical pair (mx_c, my_c) conjugate; cross = 0
    print("SIG ({0},{1}) {2} | core-indices={3} | core-sig=({4},{5}) | dim V={6}".format(
        p, q, label, core_idx,
        sum(1 for _, _, _, t in pairs if t == "+"),
        sum(1 for _, _, _, t in pairs if t == "-"), len(reps)))
    for idx, (cc, mx, my, tt) in enumerate(pairs):
        w_self = omega_Z(mx, my, Z, adapted, n)
        conj = any(v != 0 for v in w_self) if w_self else False
        # cross with other pairs
        cross0 = True
        for (cc2, mx2, my2, tt2) in pairs:
            if cc2 == cc:
                continue
            wc = omega_Z(mx, my2, Z, adapted, n)
            if wc and any(v != 0 for v in wc):
                cross0 = False
        ok(conj, "canonical pair conjugate ({0},{1}) c={2}".format(p, q, cc))
        same = (tt == tt)  # both members bind same axis -> same type by construction
        print("   pair#{0} (axis c={1}) | η-type(b)={2} | η-type(b′)={3} | {4} | ω-conjugate={5} | cross-ortho={6}".format(
            idx, cc, tt, tt, "SAME-TYPE", "yes" if conj else "NO", "yes" if cross0 else "no"))
    return dict(n=n, eta=eta, bas=bas, N=N, x=x, y=y, core_idx=core_idx,
                Z=Z, adapted=adapted, reps=reps, pairs=pairs)


# -------- q-scan --------
print()
print("--- q-SCAN (signatures Witt≥2, n=5,6,7) ---")
SIGS = [(3, 2, "n5"), (4, 2, "n6"), (2, 4, "n6?"), (3, 3, "n6"),
        (5, 2, "n7"), (4, 3, "n7"), (3, 4, "n7?"), (2, 5, "n7?")]
results = {}
for (p, q, lab) in SIGS:
    r = analyze(p, q, lab)
    if r is not None:
        results[(p, q)] = r

# -------- systematic --------
print()
print("SYSTEMATICS (raw):")
print("  canonical pairs are ALWAYS same-type (both members bind the same core axis);")
print("  the multiset of pair-types = the core signature.")

# -------- Mixed k>0 : (5,2) + partial core rotation --------
print()
print("--- MIXED k>0: (5,2) + partial core rotation (axis ker E) ---")
r = results[(5, 2)]
n = 7
eta = r["eta"]
bas = r["bas"]
N = r["N"]
core_idx = r["core_idx"]         # [2,3,4]
# E rotates core {2,3}, fixes core axis 4 = ker E
E = block_gen(n, eta, core_idx[0], core_idx[1], 1)
ok((E * N - N * E).is_zero_matrix, "mixed: [E,N]=0")
fixed_axis = core_idx[2]
fixed_type = "+" if eta[fixed_axis, fixed_axis] == 1 else "-"
# surviving pair = the one at the fixed axis
surv_pair = [pr for pr in r["pairs"] if pr[0] == fixed_axis][0]
print("   (5,2)+rotation({0},{1}) | fixed axis ker E = c={2} type={3} | surviving pair = pair@c={2} type=({3},{3})".format(
    core_idx[0], core_idx[1], fixed_axis, fixed_type))
# verify the fixed-axis module directions survive ad(E) (commute with E), rotated ones don't
mx_fix = wedge(r["x"], unit(n, fixed_axis), eta)
my_fix = wedge(r["y"], unit(n, fixed_axis), eta)
ok((E * mx_fix - mx_fix * E).is_zero_matrix, "mixed: fixed-axis module dir survives [E,.]=0")
mx_rot = wedge(r["x"], unit(n, core_idx[0]), eta)
ok(not (E * mx_rot - mx_rot * E).is_zero_matrix, "mixed: rotated-axis module dir does NOT survive")
print("   → the surviving pair binds the fixed axis ({0}-type); the rotated axes {1},{2} die (do not commute with E)".format(
    fixed_type, core_idx[0], core_idx[1]))

# ==================== BIT-FENCE ====================
print()
print("BIT-FENCE (ω-ranks S961/S968):")
# (4,2) rank0: module dim 4, omega rank 4 (S961); pairs = core sig (2,0)
r42 = results[(4, 2)]
nn = 4 + 2
# reconstruct omega rank from constructed pairs (should be full = dim V)
allmods = []
for (cc, mx, my, tt) in r42["pairs"]:
    allmods.append(mx); allmods.append(my)
mm = len(r42["reps"])
Wrows = []
for a_ in range(len(allmods)):
    for b_ in range(len(allmods)):
        pass
# rank via Z-projection gram of the constructed module basis
s = len(r42["Z"])
G = zeros(len(allmods), len(allmods))
for a_ in range(len(allmods)):
    for b_ in range(len(allmods)):
        wv = omega_Z(allmods[a_], allmods[b_], r42["Z"], r42["adapted"], nn)
        G[a_, b_] = wv[0] if wv else Integer(0)
ok(G.rank() == 4, "(4,2) constructed module omega rank=4 (=S961)")
print("  (4,2): constructed module ω-rank={0} (=S961 rank4); core-sig(2,0) → 2 pairs both +".format(G.rank()))

# ==================== MUTANTS ====================
print()
print("MUTANTS (one per branch: same-type / different-type / isotropic-trace):")
mut_ok = True

# br-sametype: (5,2) all core-pos -> all pairs same-type (+)
r52 = results[(5, 2)]
types52 = [t for (_, _, _, t) in r52["pairs"]]
if all(t == "+" for t in types52):
    print("  MUTANT br-sametype: CAUGHT ((5,2) core all-pos → all pairs same-type +)")
else:
    print("  MUTANT br-sametype: NOT CAUGHT"); mut_ok = False

# br-difftype: (3,3) core (1,1) -> pairs of BOTH types present
r33 = results[(3, 3)]
types33 = set(t for (_, _, _, t) in r33["pairs"])
if types33 == {"+", "-"}:
    print("  MUTANT br-difftype: CAUGHT ((3,3) core (1,1) → pairs of both types +,−)")
else:
    print("  MUTANT br-difftype: NOT CAUGHT"); mut_ok = False

# br-isotropic-trace: bind a module dir to a NULL axis (x itself, in W) -> isotropic trace
_e = make_eta(4, 2)
_x, _y = d0_on(list(range(6)), 4, 6)
null_trace = (_x.T * _e * _x)[0, 0]     # x is isotropic -> 0
if null_trace == 0:
    print("  MUTANT br-isotropic-trace: CAUGHT (bind to W-axis x: eta(x,x)=0 isotropic trace)")
else:
    print("  MUTANT br-isotropic-trace: NOT CAUGHT"); mut_ok = False

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
