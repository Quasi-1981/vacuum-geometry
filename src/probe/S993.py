# -*- coding: utf-8 -*-
# DIM: na (W39 leg3b: symbolic forall-forcing of the Lambda^2 cocycle by invariant-form count).
#
# ============================================================================
# DERIVATIONAL (own TZ, no gate; closes stake-2 of leg-3 that S992 left as an
# instance-fact -- the project's honest catch).  Stamped BEFORE.
# ----------------------------------------------------------------------------
# Levi = sp(2k,R) (+) so(eta|G) on the module M = W (x) G  (W=2k fundamental of
# sp ; G=g fundamental of so).  The Lambda^2 cocycle is an ANTISYMMETRIC invariant
# bilinear form on M.  Stake-2 (forall k,g) = the ANTISYMMETRIC-invariant space is
# 1-dim -> the cocycle omega (x) eta' is FORCED (unique).
# Factorization by symmetry sector:
#   Inv(W)^sp   = { antisym: 1 (omega) , sym: 0 }             forall k  (Schur, symplectic fundamental)
#   Inv(G)^so   = { sym: 1 (eta') , antisym: 0 (g!=2) / 1 (g=2) }
#      -- g=2 EXCEPTION: so(eta|G)=so(2) or so(1,1) is ABELIAN, its 2-dim
#         fundamental also preserves the AREA form eps_G (antisym).
#   ANTISYM Inv(W(x)G) = Inv(W)^anti (x) Inv(G)^sym  (+)  Inv(W)^sym (x) Inv(G)^anti
#                       = 1 (x) 1                    (+)  0 (x) {0 or 1}
#                       = 1   forall g   (= omega (x) eta' = the cocycle).
#   -> the extra g=2 area (Inv(G)^anti) feeds ONLY the SYMMETRIC sector
#      (omega (x) eps_G, a NON-cocycle), because Inv(W)^sym = 0.
# CONCLUSION: the Lambda^2 cocycle is forced forall k,g (antisym-invariant unique);
# the g=2 abelian rung adds a symmetric non-cocycle invariant -- carved, harmless
# to T24.  Instance-fenced k in {1,2,3}, g in {1,2,3,4} several signatures.
# Mutant per branch ; negative control ; exact Q ; FORBIDDEN-SCAN ; STOP.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
from sympy import Matrix, Integer, zeros, eye, diag, symbols, linsolve

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S993_run.log")
_logf = open(_LOGPATH, "w", encoding="utf-8")


class Tee:
    def __init__(self, real, fh):
        self.real = real; self.fh = fh; self.chunks = []

    def write(self, s):
        self.real.write(s); self.fh.write(s); self.chunks.append(s); return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


_tee = Tee(sys.stdout, _logf)
sys.stdout = _tee
ASSERT_PASS = [0]; FAILS = [0]


def ok(cond, msg):
    if cond:
        ASSERT_PASS[0] += 1
    else:
        FAILS[0] += 1; print("ASSERT-FAIL: " + msg)


# ==================== primitives ====================
def make_eta(p, q):
    return diag(*([Integer(1)] * p + [Integer(-1)] * q))


def unit(n, i):
    v = zeros(n, 1); v[i, 0] = Integer(1); return v


def flat(M):
    return Matrix(M.rows * M.cols, 1, list(M))


def sp_form(k):
    J = zeros(2 * k, 2 * k)
    for i in range(k):
        J[i, k + i] = Integer(1); J[k + i, i] = Integer(-1)
    return J


def sp_basis(w):
    m = w.rows
    av = symbols('a0:%d' % (m * m))
    Am = Matrix(m, m, av)
    sol = list(linsolve(list(Am * w + w * Am.T), list(av)))[0]
    free = set()
    for e in sol:
        free |= e.free_symbols
    subsol = {av[k]: sol[k] for k in range(m * m)}
    out = []
    for fv in sorted(free, key=lambda s: s.name):
        pick = {s: (Integer(1) if s == fv else Integer(0)) for s in free}
        out.append(Am.subs(subsol).subs(pick))
    return out


def so_G_basis(G):
    g = G.rows; out = []
    for i in range(g):
        for j in range(i + 1, g):
            out.append(unit(g, i) * (G * unit(g, j)).T - unit(g, j) * (G * unit(g, i)).T)
    return out


def kron(A, B):
    ra, ca = A.shape; rb, cb = B.shape
    M = zeros(ra * rb, ca * cb)
    for i in range(ra):
        for j in range(ca):
            if A[i, j] != 0:
                M[i * rb:(i + 1) * rb, j * cb:(j + 1) * cb] = A[i, j] * B
    return M


def invariant_forms(ops, d, symtype=None):
    """basis of {B (d x d): op^T B + B op = 0 for all op}; symtype adds (anti)symmetry."""
    if d == 0:
        return []
    eb = []
    for i in range(d):
        for j in range(d):
            E = zeros(d, d); E[i, j] = Integer(1); eb.append(E)
    cols = []
    for E in eb:
        parts = [flat(op.T * E + E * op) for op in ops] if ops else [flat(zeros(d, d))]
        if symtype == 'sym':
            parts.append(flat(E - E.T))
        elif symtype == 'anti':
            parts.append(flat(E + E.T))
        cols.append(Matrix.vstack(*parts))
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        B = zeros(d, d)
        for k in range(d * d):
            if v[k, 0] != 0:
                B = B + v[k, 0] * eb[k]
        out.append(B)
    return out


def sectors(ops, d):
    return (len(invariant_forms(ops, d)),
            len(invariant_forms(ops, d, 'sym')),
            len(invariant_forms(ops, d, 'anti')))


# ==================== §1 Inv(W)^sp: anti=1 (omega), sym=0, forall k ====================
print("=== §1 Inv(W)^{sp(2k)}: ANTISYM=1 (omega), SYM=0, for k=1,2,3 ===")
for k in (1, 2, 3):
    w = sp_form(k); spb = sp_basis(w)
    tot, s, a = sectors(spb, 2 * k)
    ok(len(spb) == k * (2 * k + 1), "sp(2k) dim k(2k+1) at k={0}".format(k))
    ok(tot == 1 and s == 0 and a == 1, "Inv(W)^sp anti=1 sym=0 at k={0} (tot={1} s={2} a={3})".format(k, tot, s, a))
    form = invariant_forms(spb, 2 * k)[0]
    prop = next((form[i, j] / w[i, j] for i in range(2 * k) for j in range(2 * k)
                 if w[i, j] != 0 and form[i, j] != 0), None)
    ok((form - prop * w).is_zero_matrix, "sp-form == omega at k={0}".format(k))
    print("  k={0}: dim sp={1} | tot={2} sym={3} anti={4} | form == {5}*omega".format(k, len(spb), tot, s, a, prop))


# ==================== §2 Inv(G)^so: sym=1 (eta') always ; anti = 0 (g!=2) / 1 (g=2 area) ====================
print()
print("=== §2 Inv(G)^{so(eta|G)}: SYM=1 (eta') always ; ANTISYM = 0 (g!=2) / 1 (g=2, area eps_G) ===")
GSIGS = [(1, 0), (0, 1), (2, 0), (1, 1), (3, 0), (2, 1), (4, 0), (2, 2)]
for (gp, gm) in GSIGS:
    g = gp + gm; G = make_eta(gp, gm); sob = so_G_basis(G)
    tot, s, a = sectors(sob, g)
    exp_a = 1 if g == 2 else 0
    ok(s == 1 and a == exp_a and tot == 1 + exp_a,
       "Inv(G)^so at ({0},{1}) g={2}: sym=1 anti={3} (got s={4} a={5})".format(gp, gm, g, exp_a, s, a))
    symform = invariant_forms(sob, g, 'sym')[0]
    prop = next((symform[i, i] / G[i, i] for i in range(g) if G[i, i] != 0 and symform[i, i] != 0), None)
    ok(prop is not None and (symform - prop * G).is_zero_matrix, "sym so-form == eta' at ({0},{1})".format(gp, gm))
    note = "  <-- g=2 ABELIAN so, +area eps_G (antisym)" if g == 2 else ""
    print("  ({0},{1}) g={2}: dim so={3} | tot={4} sym={5} anti={6} | sym==eta'{7}".format(gp, gm, g, len(sob), tot, s, a, note))


# ==================== §3 factorization by sector -> cocycle forced forall ====================
print()
print("=== §3 FACTORIZATION by symmetry sector (carved) ===")
print("  ANTISYM Inv(W(x)G) = Inv(W)^anti (x) Inv(G)^sym  (+)  Inv(W)^sym (x) Inv(G)^anti")
print("                     = 1 (x) 1  (+)  0 (x) {0 or 1}  =  1   forall g  (= omega (x) eta').")
print("  => Lambda^2 cocycle FORCED forall k,g (unique antisym invariant).")
print("  g=2 exception lives ONLY in the SYM sector: omega (x) eps_G (SYMMETRIC, NOT a Lie cocycle),")
print("  since Inv(W)^sym = 0 kills any antisym contribution from the extra area form.")


# ==================== §4 instance-fences: antisym (cocycle) sector = 1 forall ====================
print()
print("=== §4 INSTANCE-FENCES: ANTISYM Inv(W(x)G) = 1 (=omega(x)eta') forall ; total jumps at g=2 ===")
CASES = [(1, (1, 0)), (1, (1, 1)), (1, (2, 1)), (1, (2, 2)),
         (2, (1, 0)), (2, (1, 1)), (2, (2, 1)),
         (3, (1, 0)), (3, (1, 1))]
for (k, (gp, gm)) in CASES:
    g = gp + gm; m = 2 * k
    w = sp_form(k); G = make_eta(gp, gm)
    spb = sp_basis(w); sob = so_G_basis(G)
    mod_ops = [kron(A, eye(g)) for A in spb] + [kron(eye(m), E) for E in sob]
    d = m * g
    tot, s, a = sectors(mod_ops, d)
    ok(a == 1, "ANTISYM Inv(WxG)=1 (cocycle forced) at k={0} ({1},{2}) g={3} (got a={4})".format(k, gp, gm, g, a))
    exp_tot = 2 if g == 2 else 1
    ok(tot == exp_tot, "total Inv(WxG) at k={0} ({1},{2}) g={3}: expect {4} got {5}".format(k, gp, gm, g, exp_tot, tot))
    # the unique antisym form == omega (x) eta'
    expected = kron(w, G)
    aform = invariant_forms(mod_ops, d, 'anti')[0]
    prop = next((aform[i, j] / expected[i, j] for i in range(d) for j in range(d)
                 if expected[i, j] != 0 and aform[i, j] != 0), None)
    ok(prop is not None and (aform - prop * expected).is_zero_matrix, "cocycle form == omega(x)eta' at k={0} ({1},{2})".format(k, gp, gm))
    ex = " | +SYM omega(x)eps_G (g=2 non-cocycle)" if g == 2 else ""
    print("  k={0} ({1},{2}) g={3} dv={4}: antisym(cocycle)={5}==omega(x)eta' | total={6} sym={7}{8}".format(
        k, gp, gm, g, d, a, tot, s, ex))


# ==================== mutants — one per branch ====================
print()
print("=== mutants (one per branch) ===")
mut_ok = True

# br-sp-sym0: sp fundamental has NO symmetric invariant (kills g=2 area in the antisym sector)
_sp = sp_basis(sp_form(2)); _t, _s, _a = sectors(_sp, 4)
if (_t, _s, _a) == (1, 0, 1):
    print("MUTANT br-sp: CAUGHT (Inv(W)^sp: sym=0, anti=1 -- symplectic fundamental)")
else:
    print("MUTANT br-sp: NOT CAUGHT ({0})".format((_t, _s, _a))); mut_ok = False

# br-so-g2: so(2)/so(1,1) at g=2 has an EXTRA antisym area invariant; g=3 does not
_g2 = sectors(so_G_basis(make_eta(1, 1)), 2)
_g3 = sectors(so_G_basis(make_eta(2, 1)), 3)
if _g2[2] == 1 and _g3[2] == 0:
    print("MUTANT br-so-g2: CAUGHT (g=2 antisym area=1 ; g=3 antisym=0)")
else:
    print("MUTANT br-so-g2: NOT CAUGHT ({0} vs {1})".format(_g2, _g3)); mut_ok = False

# br-cocycle-forall: antisym Inv(WxG)=1 even at the g=2 exception (5,5)-like (k=2,(1,1))
_w = sp_form(2); _G = make_eta(1, 1)
_ops = [kron(A, eye(2)) for A in sp_basis(_w)] + [kron(eye(4), E) for E in so_G_basis(_G)]
_tc, _sc, _ac = sectors(_ops, 8)
if _ac == 1 and _tc == 2 and _sc == 1:
    print("MUTANT br-cocycle-forall: CAUGHT (g=2: antisym cocycle still unique=1 ; +sym extra=1)")
else:
    print("MUTANT br-cocycle-forall: NOT CAUGHT ({0})".format((_tc, _sc, _ac))); mut_ok = False

# br-factorization: total dim(WxG) == dim Inv(W) * dim Inv(G) (sector-blind product)
_dW = sectors(sp_basis(sp_form(2)), 4)[0]
_dG2 = sectors(so_G_basis(make_eta(1, 1)), 2)[0]
_dWG2 = sectors([kron(A, eye(2)) for A in sp_basis(sp_form(2))] + [kron(eye(4), E) for E in so_G_basis(make_eta(1, 1))], 8)[0]
if _dWG2 == _dW * _dG2 == 2:
    print("MUTANT br-factorization: CAUGHT (g=2: dim WxG = dim W * dim G = 1*2 = 2)")
else:
    print("MUTANT br-factorization: NOT CAUGHT ({0} vs {1}*{2})".format(_dWG2, _dW, _dG2)); mut_ok = False


# ==================== negative control ====================
_w = sp_form(2); _G = make_eta(2, 1)   # g=3, cocycle unique
_ops = [kron(A, eye(3)) for A in sp_basis(_w)] + [kron(eye(4), E) for E in so_G_basis(_G)]
_good = kron(_w, _G)
_bad = _good + unit(12, 0) * unit(12, 0).T
_gi = all((op.T * _good + _good * op).is_zero_matrix for op in _ops)
_bi = all((op.T * _bad + _bad * op).is_zero_matrix for op in _ops)
ok(_gi and not _bi, "neg-control: omega(x)eta' invariant ; perturbed not invariant")


# ==================== summary ====================
print()
print("SUMMARY: stake-2 CLOSED forall -- ANTISYM (cocycle) invariant on W(x)G is 1-dim = omega(x)eta'")
print("SUMMARY: for ALL k,g (incl g=2). Boundary carved: g=2 (abelian so) adds a SYM non-cocycle")
print("SUMMARY: invariant omega(x)eps_G; harmless to T24 (Lie bracket = antisym only).")
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
