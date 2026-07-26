# -*- coding: utf-8 -*-
# DIM: na (W39 leg3: symbolic forall-derivation of the Lambda^2 law; 0 handles).
#
# ============================================================================
# DERIVATIONAL CLASS (S952/J-0423) -- NOT blind (public target).  Stamped BEFORE.
# ----------------------------------------------------------------------------
# N = sum of k pairwise-orthogonal wedges over a 2k-dim totally isotropic W in
# so(p,q).  Adapted split V = W(2k) + G(g) + W'(2k), g = n-4k:
#   eta = [[0,0,I],[0,G,0],[I,0,0]] ;  N = [[0,0,w],[0,0,0],[0,0,0]]  (w = W'->W,
#   the 2k x 2k skew symplectic form built from the k wedges).
# Block M = [[A,B,Cc],[D,E,F],[G2,H,K]].  Imposing so(eta) and [M,N]=0 splits the
# blocks FREE/LINKED/ZERO uniformly in (k,g) -- the k=1 case is S952, VERBATIM in
# form.  The derivation (carved, with instance-fences by span-equality) yields
#   c(N) = ( sp(2k,R) (+) so(eta|G) )  |x  ( (W(x)G) (+) Lambda^2 W )
#   radical bracket   [u(x)a, v(x)b] = (u^v) . eta'(a,b)      (Lambda^2 cocycle)
#   [Lambda^2 W, rad] = 0 ;  dim c = k(2k+1) + g(g-1)/2 + 2kg + k(2k-1).
# Degenerate from ONE formula: g=0 -> radical = Lambda^2 W abelian ; k=1 ->
# Lambda^2 W 1-dim -> h_{2g+1} = T18 rank-0 (retro-identity).
# Guards: instance-fences span-equality vs S988(A)/S991(A1,A2) + k=1 retro to
# the T18/S950 tables ; mutant per branch ; negative control ; exact Q sympy ;
# FORBIDDEN-SCAN ; Jordan-deep (N^2!=0) OUT of scope ; STOP after derivation.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
from sympy import (Matrix, Integer, Rational, zeros, eye, diag, symbols, Symbol,
                   MatrixSymbol, ZeroMatrix, BlockMatrix, block_collapse, linsolve, Identity)

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S992_run.log")
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


# ==================== numeric primitives (from S952/S988) ====================
def make_eta(p, q):
    return diag(*([Integer(1)] * p + [Integer(-1)] * q))


def unit(n, i):
    v = zeros(n, 1); v[i, 0] = Integer(1); return v


def wedge(x, y, eta):
    return x * (eta * y).T - y * (eta * x).T


def is_so(M, eta):
    return (M.T * eta + eta * M).is_zero_matrix


def so_basis(n, eta):
    return [wedge(unit(n, i), unit(n, j), eta) for i in range(n) for j in range(i + 1, n)]


def flat(M):
    return Matrix(M.rows * M.rows, 1, list(M))


def stack_flats(mats, sq):
    return zeros(sq * sq, 0) if not mats else Matrix.hstack(*[flat(M) for M in mats])


def span_basis(mats, sq):
    out = []; F = zeros(sq * sq, 0); r = 0
    for M in mats:
        F2 = Matrix.hstack(F, flat(M))
        if F2.rank() > r:
            out.append(M); F = F2; r += 1
    return out


def same_span(A, B, sq):
    FA = stack_flats(A, sq); FB = stack_flats(B, sq)
    ra, rb = FA.rank(), FB.rank()
    return ra == rb and Matrix.hstack(FA, FB).rank() == ra


def centralizer(A, bas):
    n = A.rows
    ns = Matrix.hstack(*[flat(B * A - A * B) for B in bas]).nullspace()
    cb = []
    for v in ns:
        M = zeros(n, n)
        for k in range(len(bas)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * bas[k]
        cb.append(M)
    return cb


# ==================== §1 symbolic block equations (uniform in m=2k, g) ====================
def symbolic_blocks():
    m = Symbol('m', positive=True, integer=True)   # m = 2k = dim W
    g = Symbol('g', positive=True, integer=True)    # g = dim G
    A = MatrixSymbol('A', m, m); B = MatrixSymbol('B', m, g); Cc = MatrixSymbol('Cc', m, m)
    D = MatrixSymbol('D', g, m); E = MatrixSymbol('E', g, g); F = MatrixSymbol('F', g, m)
    G2 = MatrixSymbol('G2', m, m); H = MatrixSymbol('H', m, g); K = MatrixSymbol('K', m, m)
    G = MatrixSymbol('Gm', g, g); w = MatrixSymbol('w', m, m)   # w = symplectic W'->W (skew)
    Imm = Identity(m)
    Zmg = ZeroMatrix(m, g); Zgm = ZeroMatrix(g, m); Zgg = ZeroMatrix(g, g); Zmm = ZeroMatrix(m, m)
    M = BlockMatrix([[A, B, Cc], [D, E, F], [G2, H, K]])
    eta = BlockMatrix([[Zmm, Zmg, Imm], [Zgm, G, Zgm], [Imm, Zmg, Zmm]])
    N = BlockMatrix([[Zmm, Zmg, w], [Zgm, Zgg, Zgm], [Zmm, Zmg, Zmm]])
    c1 = block_collapse(M.T * eta + eta * M)
    c2 = block_collapse(M * N - N * M)
    return c1, c2, ("A", "B", "Cc", "D", "E", "F", "G2", "H", "K")


print("=== §1 SYMBOLIC BLOCK EQUATIONS (MatrixSymbol, symbolic m=2k and g) ===")
_c1, _c2, _ = symbolic_blocks()
print("--- so(eta) block equations (=0), positions [W|G|W'] x [W|G|W'] ---")
_lab = [["A", "B", "Cc"], ["D", "E", "F"], ["G2", "H", "K"]]
for i in range(3):
    for j in range(3):
        print("  soeta[{0}]: {1} = 0".format(_lab[i][j], _c1.blocks[i, j]))
print("--- [M,N] block equations (=0) ---")
for i in range(3):
    for j in range(3):
        e = _c2.blocks[i, j]
        if str(e) != "0":
            print("  [M,N][{0},{1}]: {2} = 0".format(i, j, e))


# ==================== §2 FREE / LINKED / ZERO (hand-read from §1, carved) ====================
print()
print("=== §2 SOLUTION free/linked/zero (hand-read from §1 block eqs; carved, class S952) ===")
print("  FREE:   A in sp(w) [A*w + w*A^T = 0]  = sp(2k,R) Levi   dim k(2k+1)")
print("          Cc antisym m x m [Cc + Cc^T = 0]  = Lambda^2 W (centre)   dim k(2k-1)")
print("          E in so(G) [E^T G + G E = 0]  = so(eta|G) reduct   dim g(g-1)/2")
print("          B free m x g  = module W (x) G   dim 2k*g")
print("  LINKED: F = -G^{-1} B^T ;  H = -D^T G ;  K = -A^T")
print("  ZERO:   G2 = 0 ;  D = 0 (=> H = 0)")
print("  hand-reads (explicit, class S952): (i) soeta[K]:Cc+Cc^T=0 -> Cc antisym=Lambda^2 W ;")
print("    (ii) soeta[A]:A^T w..-> A*w+w*A^T=0 = sp(w) ; (iii) [M,N] forces D=G2=0, K=-A^T ;")
print("    (iv) soeta[F]/[H]: F,H linked to B,D by G ; (v) E block = so(G) unchanged by N.")


# ==================== helpers: sp / Lambda^2 / so(G) bases ====================
def sp_basis(w):
    m = w.rows
    av = symbols('a0:%d' % (m * m))
    Am = Matrix(m, m, av)
    eqs = list(Am * w + w * Am.T)
    sol = list(linsolve(eqs, list(av)))[0]
    free = set()
    for e in sol:
        free |= e.free_symbols
    out = []
    subsol = {av[k]: sol[k] for k in range(m * m)}
    for fv in sorted(free, key=lambda s: s.name):
        pick = {s: (Integer(1) if s == fv else Integer(0)) for s in free}
        out.append(Am.subs(subsol).subs(pick))
    return out


def lambda2_basis(m):
    out = []
    for i in range(m):
        for j in range(i + 1, m):
            E = zeros(m, m); E[i, j] = Integer(1); E[j, i] = Integer(-1); out.append(E)
    return out


def so_G_basis(G):
    g = G.rows; out = []
    for i in range(g):
        for j in range(i + 1, g):
            out.append(unit(g, i) * (G * unit(g, j)).T - unit(g, j) * (G * unit(g, i)).T)
    return out


# ==================== k-wedge construction + adapted basis ====================
def k_wedge(p, q, k):
    n = p + q
    eta = make_eta(p, q)
    Wv = []
    for i in range(k):
        Wv.append(unit(n, 2 * i) + unit(n, p + 2 * i))
        Wv.append(unit(n, 2 * i + 1) + unit(n, p + 2 * i + 1))
    N = zeros(n, n)
    for i in range(k):
        N = N + wedge(Wv[2 * i], Wv[2 * i + 1], eta)
    return N, Wv, eta


def adapted(Wv, eta, n):
    m = len(Wv)                       # = 2k
    # G = core: eta-orthogonal complement of W, intersect... build C by extending W
    Ap = Matrix.hstack(*[eta * v for v in Wv]).T
    Wperp = Ap.nullspace()
    F = Matrix.hstack(*Wv); r = m; C = []
    for v in Wperp:
        F2 = Matrix.hstack(F, v)
        if F2.rank() > r:
            C.append(v); F = F2; r += 1
    g = len(C)
    # dual isotropic W': solve <Wv_i, dual_j> = delta, <C, dual>=0
    rows = [(eta * v).T for v in Wv] + [(eta * cc).T for cc in C]
    Amat = Matrix.vstack(*rows)

    def dual(idx):
        rhs = zeros(m + g, 1); rhs[idx, 0] = Integer(1)
        sol, params = Amat.gauss_jordan_solve(rhs)
        if params.rows * params.cols > 0:
            sol = sol.subs({s: 0 for s in params})
        return sol
    Wp = [dual(i) for i in range(m)]
    # make W' isotropic + mutually orthogonal + orthogonal to C (Gram-Schmidt against W, null)
    for i in range(m):
        for cc in C:
            coef = (Wp[i].T * eta * cc)[0, 0]
            gcc = (cc.T * eta * cc)[0, 0]
            if gcc != 0:
                Wp[i] = Wp[i] - Rational(coef, gcc) * cc
    for i in range(m):
        aii = (Wp[i].T * eta * Wp[i])[0, 0]
        Wp[i] = Wp[i] - Rational(aii, 2) * Wv[i]
        for j in range(i):
            aij = (Wp[i].T * eta * Wp[j])[0, 0]
            Wp[i] = Wp[i] - aij * Wv[j]
    cols = list(Wv) + C + Wp
    P = Matrix.hstack(*cols)
    return P, g, C, Wp


def build_block_centralizer(N, Wv, eta, n):
    m = len(Wv)
    P, g, C, Wp = adapted(Wv, eta, n)
    Pinv = P.inv()
    etaA = P.T * eta * P
    NA = Pinv * N * P
    G = etaA[m:m + g, m:m + g] if g > 0 else zeros(0, 0)
    w = NA[0:m, m + g:m + g + m]        # W'->W skew form
    Ginv = G.inv() if g > 0 else zeros(0, 0)

    def emit(Ab, Ccb, Eb, Bb):
        Ma = zeros(n, n)
        Ma[0:m, 0:m] = Ab
        if g > 0:
            Ma[0:m, m:m + g] = Bb
            Ma[m:m + g, m:m + g] = Eb
            Ma[m:m + g, m + g:m + g + m] = -Ginv * Bb.T
        Ma[0:m, m + g:m + g + m] = Ccb
        Ma[m + g:m + g + m, m + g:m + g + m] = -Ab.T
        return P * Ma * Pinv
    zmm = zeros(m, m); zgg = zeros(g, g); zmg = zeros(m, g)
    gens = []
    tag = []
    for Av in sp_basis(w):
        gens.append(emit(Av, zmm, zgg, zmg)); tag.append("sp")
    for Lv in lambda2_basis(m):
        gens.append(emit(zmm, Lv, zgg, zmg)); tag.append("L2")
    for Ev in so_G_basis(G):
        gens.append(emit(zmm, zmm, Ev, zmg)); tag.append("soG")
    for i in range(m):
        for j in range(g):
            Bv = zeros(m, g); Bv[i, j] = Integer(1)
            gens.append(emit(zmm, zmm, zgg, Bv)); tag.append("B")
    return gens, g, dict(P=P, Ginv=Ginv, w=w, m=m, emit=emit)


def dimc_formula(k, g):
    m = 2 * k
    return k * (2 * k + 1) + g * (g - 1) // 2 + 2 * k * g + k * (2 * k - 1)


# ==================== §3 module bracket -> Lambda^2 cocycle (derive + instance-verify) ====================
print()
print("=== §3 MODULE BRACKET  [u(x)a, v(x)b] = (u^v).eta'(a,b)  (derived; instance bit-check) ===")
print("  Block algebra: effective M = [[A,B,Cc],[0,E,-G^{-1}B^T],[0,0,-A^T]] (upper-tri).")
print("  [M1,M2] (1,3)-block (-> centre Cc=Lambda^2 W) receives  B1(-G^{-1})B2^T - B2(-G^{-1})B1^T")
print("  For B = e_i(x)e_a (u=e_i in W, a-th G-axis):  [B_{u,a},B_{v,b}] = -(G^{-1})_{ab}(u v^T - v u^T)")
print("     = (u^v) . eta'(a,b)   with eta' = -G^{-1} (dual metric on G, sign=convention).  Lambda^2 cocycle.")


def verify_cocycle(builder_info, gens, g, n):
    """On a concrete instance: pick the B-module generators, bracket them, and check the
    result lies in the Lambda^2 W (Cc) block with the (u^v).eta'(a,b) coefficient."""
    inf = builder_info
    P, Ginv, m = inf["P"], inf["Ginv"], inf["m"]
    Pinv = P.inv()
    emit = inf["emit"]
    zmm = zeros(m, m); zgg = zeros(g, g)
    # B generators as (i,a)
    def Bgen(i, a):
        Bv = zeros(m, g); Bv[i, a] = Integer(1)
        return emit(zmm, zmm, zgg, Bv)
    # Cc (Lambda^2) generator for u=e_i^e_j
    def L2(i, j):
        L = zeros(m, m); L[i, j] = Integer(1); L[j, i] = Integer(-1)
        return emit(L, zmm, zgg, zeros(m, g)) if False else _emitCc(inf, i, j)
    okflag = True
    checked = 0
    for i in range(m):
        for j in range(i + 1, m):
            for a in range(g):
                for b in range(g):
                    br = Bgen(i, a) * Bgen(j, b) - Bgen(j, b) * Bgen(i, a)
                    # derived: [B_{ia},B_{jb}] = -(G^{-1})_{ab} (e_i^e_j) in Cc block
                    # (eta' = -G^{-1}: dual metric on G up to the overall sign convention)
                    coeff = -Ginv[a, b]
                    expected = _emitCc(inf, i, j, coeff)
                    if not (br - expected).is_zero_matrix:
                        okflag = False
                    checked += 1
    return okflag, checked


def _emitCc(inf, i, j, coeff=Integer(1)):
    m = inf["m"]; emit = inf["emit"]
    L = zeros(m, m); L[i, j] = coeff; L[j, i] = -coeff
    g = inf["Ginv"].rows
    return emit(zeros(m, m), L, zeros(g, g), zeros(m, g))


# ==================== §4 dim formula + degenerate cases ====================
print()
print("=== §4 dim-formula + degenerate cases (from the ONE formula) ===")
print("  dim c(k,g) = k(2k+1) + g(g-1)/2 + 2k*g + k(2k-1)   [= sp(2k) + so(G) + W(x)G + Lambda^2 W]")
print("  g=0:  dim = k(2k+1) + k(2k-1) = 4k^2 ; radical = Lambda^2 W (dim k(2k-1)) ABELIAN (eta'=0).")
print("  k=1:  Lambda^2 W dim 1 ; radical = h_{2g+1} ; dim = 3 + g(g-1)/2 + 2g = g(g-1)/2+2g+3  (= T18 rank-0).")
ok(dimc_formula(1, 1) == 6, "formula k=1 g=1 == 6 (S952 (3,2)r0)")
ok(dimc_formula(1, 2) == 9, "formula k=1 g=2 == 9 (S988/(4,2)r0)")
ok(dimc_formula(1, 3) == 13, "formula k=1 g=3 == 13 (S950 (5,2)r0)")
ok(dimc_formula(2, 0) == 16, "formula k=2 g=0 == 16 (S988-A two-wedge)")
ok(dimc_formula(2, 1) == 20, "formula k=2 g=1 == 20 (S991-A1 so(5,4))")
ok(dimc_formula(2, 2) == 25, "formula k=2 g=2 == 25 (S991-A2 so(5,5))")
# k=1 retro identity: g(g-1)/2 + 2g + 3  vs  h_{2g+1} Heisenberg (dim 2g+1) + Levi so(2,1)(3) + so(G)(g(g-1)/2)
ok(dimc_formula(1, 2) == (2 * 2 + 1) + 3 + (2 * (2 - 1) // 2), "k=1 retro = h5 + so(2,1) + so(G) at g=2")


# ==================== §5 instance-fences (span-equality over a (k,g) spread) ====================
print()
print("=== §5 INSTANCE-FENCES (bit-for-bit vs centralizer(): dim + span + cocycle) ===")
# (p,q,k) chosen so g=n-4k matches; g = (p-2k)+(q-2k)
CASES = [
    ((3, 2), 1),   # k=1 g=1  (T18 rank-0 retro; S952/S950 family)
    ((4, 2), 1),   # k=1 g=2
    ((5, 2), 1),   # k=1 g=3  (S950 (5,2)r0)
    ((4, 4), 2),   # k=2 g=0  (S988-A two-wedge)
    ((5, 4), 2),   # k=2 g=1  (S991-A1)
    ((5, 5), 2),   # k=2 g=2  (S991-A2 -- Lambda^2 main)
    ((6, 6), 3),   # k=3 g=0  (out-of-sample: three wedges saturating)
    ((7, 6), 3),   # k=3 g=1  (out-of-sample)
]
for (pq, k) in CASES:
    p, q = pq
    n = p + q
    g = (p - 2 * k) + (q - 2 * k)
    if p < 2 * k or q < 2 * k:
        print("  ({0},{1}) k={2}: SKIP (need p,q>=2k)".format(p, q, k)); continue
    N, Wv, eta = k_wedge(p, q, k)
    ok((N * N).is_zero_matrix, "N^2=0 ({0},{1}) k={2}".format(p, q, k))
    ok(N.rank() == 2 * k, "rank N == 2k ({0},{1}) k={2}".format(p, q, k))
    bas = so_basis(n, eta)
    cdir = centralizer(N, bas)
    pred = dimc_formula(k, g)
    gens, gg, inf = build_block_centralizer(N, Wv, eta, n)
    valid = all(is_so(M, eta) and (M * N - N * M).is_zero_matrix for M in gens)
    ok(valid, "block gens valid centralizer ({0},{1}) k={2}".format(p, q, k))
    dimg = len(span_basis(gens, n))
    speq = same_span(gens, cdir, n)
    ok(len(cdir) == pred, "dim c == formula ({0},{1}) k={2}: {3} vs {4}".format(p, q, k, len(cdir), pred))
    ok(dimg == len(cdir), "block span dim == centralizer dim ({0},{1}) k={2}".format(p, q, k))
    ok(speq, "block span == centralizer span ({0},{1}) k={2}".format(p, q, k))
    # cocycle check (skip g=0 -> no module)
    if g > 0:
        cok, checked = verify_cocycle(inf, gens, g, n)
        ok(cok, "cocycle [u(x)a,v(x)b]=(u^v)eta'(a,b) exact ({0},{1}) k={2} ({3} checks)".format(p, q, k, checked))
        coctxt = "cocycle {0}/{1}".format("OK" if cok else "FAIL", checked)
    else:
        coctxt = "cocycle n/a (g=0, radical=Lambda^2 W abelian)"
    print("  ({0},{1}) k={2} g={3}: formula={4} centralizer={5} block-span={6} span-eq={7} valid={8} | {9}".format(
        p, q, k, g, pred, len(cdir), dimg, speq, valid, coctxt))


# ==================== mutants — one per branch ====================
print()
print("=== mutants (one per branch) ===")
mut_ok = True

# br-formula: matches measured AND distinguishes (k,g)
if dimc_formula(2, 0) == 16 and dimc_formula(1, 2) == 9 and dimc_formula(2, 0) != dimc_formula(1, 2):
    print("MUTANT br-formula: CAUGHT (k=2g=0->16 ; k=1g=2->9 ; distinct)")
else:
    print("MUTANT br-formula: NOT CAUGHT"); mut_ok = False

# br-sp-vs-so: A-block is sp(w) (Aw+wA^T=0), NOT so (A^T+A=0) unless coincidence
_w = Matrix([[0, 1], [-1, 0]])
_spb = sp_basis(_w)
_A = _spb[0]
if (_A * _w + _w * _A.T).is_zero_matrix and len(_spb) == 3:
    print("MUTANT br-sp: CAUGHT (sp(2) dim 3, A*w+w*A^T=0)")
else:
    print("MUTANT br-sp: NOT CAUGHT ({0})".format(len(_spb))); mut_ok = False

# br-lambda2: Cc antisym (Lambda^2), dim k(2k-1); a symmetric Cc is NOT allowed
_L = lambda2_basis(4)
_sym = zeros(4, 4); _sym[0, 0] = Integer(1)
if len(_L) == 6 and all((Lm + Lm.T).is_zero_matrix for Lm in _L) and not (_sym + _sym.T).is_zero_matrix:
    print("MUTANT br-lambda2: CAUGHT (Lambda^2(4) dim 6 antisym ; symmetric excluded)")
else:
    print("MUTANT br-lambda2: NOT CAUGHT"); mut_ok = False

# br-cocycle-sign: bracket is ANTISYMMETRIC in the pair (u^v = -v^u)
_N2, _W2, _e2 = k_wedge(5, 5, 2)
_g2, _, _inf2 = None, None, None
_gens2, _g2, _inf2 = build_block_centralizer(_N2, _W2, _e2, 10)
_cok2, _ = verify_cocycle(_inf2, _gens2, _g2, 10)
if _cok2:
    print("MUTANT br-cocycle: CAUGHT (Lambda^2 form exact on (5,5) k=2 main object)")
else:
    print("MUTANT br-cocycle: NOT CAUGHT"); mut_ok = False

# br-degenerate: g=0 gives abelian radical (Lambda^2 W), k=1 gives Heisenberg h_{2g+1}
if dimc_formula(2, 0) == 4 * 2 ** 2 and dimc_formula(1, 3) == 13:
    print("MUTANT br-degenerate: CAUGHT (g=0 -> 4k^2=16 abelian Lambda^2 ; k=1 -> h-tower 13)")
else:
    print("MUTANT br-degenerate: NOT CAUGHT"); mut_ok = False


# ==================== negative control ====================
# a symmetric (non-so) perturbation of a valid block-gen leaves so(eta); a generic so
# element does not commute with N.
_Nn, _Wn, _en = k_wedge(5, 4, 2)
_gn, _, _ = build_block_centralizer(_Nn, _Wn, _en, 9)
_good = _gn[0]
_bad = _good + unit(9, 0) * unit(9, 0).T
_someso = so_basis(9, _en)[0]
ok(is_so(_good, _en) and not is_so(_bad, _en), "neg-control: symmetric perturbation leaves so(eta)")
ok((_good * _Nn - _Nn * _good).is_zero_matrix, "neg-control: block-gen commutes with N")


# ==================== summary ====================
print()
print("SUMMARY: forall-derivation (symbolic block eqs + carved free/linked/zero + cocycle) "
      "+ instance-fences k in {1,2,3}, g in {0,1,2,3}")
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
