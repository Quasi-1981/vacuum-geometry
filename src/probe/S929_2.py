# -*- coding: utf-8 -*-
# DIM: na (a measurement table of block decompositions and restriction images in so(p,q); handles 0 — the probe postulates nothing).
import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import itertools
from math import lcm as _ilcm

from sympy import Matrix, Integer, Rational, zeros, eye, diag, symbols, Poly, gcd, sqrt as ssqrt

X_ = symbols('x')

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S929_2_run.log")
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


# ---------- exact linear-algebra primitives (S925/S927 conventions) ----------

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


def span_rank(mats, sq):
    return stack_flats(mats, sq).rank()


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


def same_span(A, B, sq):
    FA = stack_flats(A, sq)
    FB = stack_flats(B, sq)
    ra = FA.rank()
    rb = FB.rank()
    if ra != rb:
        return False
    return Matrix.hstack(FA, FB).rank() == ra


def contains(big, small, sq):
    FB = stack_flats(big, sq)
    rb = FB.rank()
    return Matrix.hstack(FB, stack_flats(small, sq)).rank() == rb


def inter_dim(A, B, sq):
    ra = span_rank(A, sq)
    rb = span_rank(B, sq)
    ru = Matrix.hstack(stack_flats(A, sq), stack_flats(B, sq)).rank()
    return ra + rb - ru


def centralizer(A, bas):
    n = A.rows
    cols = [flat(B * A - A * B) for B in bas]
    Mmap = Matrix.hstack(*cols)
    ns = Mmap.nullspace()
    cb = []
    for v in ns:
        M = zeros(n, n)
        for k in range(len(bas)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * bas[k]
        cb.append(M)
    return cb


def commutes_all(cb, A):
    return all((M * A - A * M).is_zero_matrix for M in cb)


def is_nilp(A):
    n = A.rows
    A2 = A * A
    if A2.is_zero_matrix:
        return True
    if A2.trace() != 0:
        return False
    return (A ** n).is_zero_matrix


def derived_and_closed(cb, sq):
    if not cb:
        return 0, True
    fc = stack_flats(cb, sq)
    rc = fc.rank()
    brs = []
    for a in range(len(cb)):
        for b in range(a + 1, len(cb)):
            brs.append(flat(cb[a] * cb[b] - cb[b] * cb[a]))
    if not brs:
        return 0, True
    fb = Matrix.hstack(*brs)
    dimcc = fb.rank()
    closed = (Matrix.hstack(fc, fb).rank() == rc)
    return dimcc, closed


def minpoly_poly(A):
    n = A.rows
    pows = [eye(n)]
    for _k in range(n):
        pows.append(pows[-1] * A)
    flats = [flat(P) for P in pows]
    for d in range(1, n + 1):
        Mstack = Matrix.hstack(*flats[:d])
        v = flats[d]
        if Matrix.hstack(Mstack, v).rank() == Mstack.rank():
            sol, params = Mstack.gauss_jordan_solve(v)
            if params.rows * params.cols > 0:
                sol = sol.subs({s: 0 for s in params})
            expr = X_ ** d
            for k in range(d):
                expr = expr - sol[k, 0] * X_ ** k
            return Poly(expr, X_)
    return Poly(X_ ** n, X_)  # unreachable by Cayley-Hamilton


def squarefree(mp):
    me = mp.as_expr()
    g = gcd(me, me.diff(X_))
    return Poly(g, X_).degree() == 0


def gram2(x, y, eta):
    a = (x.T * eta * x)[0, 0]
    b = (x.T * eta * y)[0, 0]
    c = (y.T * eta * y)[0, 0]
    return Matrix([[a, b], [b, c]])


def block_gen(n, eta, i, j, param):
    return Integer(param) * wedge(unit(n, i), unit(n, j), eta)


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


def cong_diag(Gin):
    # returns (D, T) with T^T G T = D diagonal (exact congruence)
    A = Gin.copy()
    d = A.rows
    T = eye(d)
    for i in range(d):
        if A[i, i] == 0:
            js = None
            for j in range(i + 1, d):
                if A[j, j] != 0:
                    js = j
                    break
            if js is not None:
                A.col_swap(i, js)
                A.row_swap(i, js)
                T.col_swap(i, js)
            else:
                jn = None
                for j in range(i + 1, d):
                    if A[i, j] != 0:
                        jn = j
                        break
                if jn is None:
                    continue
                A[:, i] = A[:, i] + A[:, jn]
                A[i, :] = A[i, :] + A[jn, :]
                T[:, i] = T[:, i] + T[:, jn]
        piv = A[i, i]
        if piv == 0:
            continue
        for j in range(i + 1, d):
            if A[i, j] != 0:
                f = A[i, j] / piv
                A[:, j] = A[:, j] - f * A[:, i]
                A[j, :] = A[j, :] - f * A[i, :]
                T[:, j] = T[:, j] - f * T[:, i]
    return A, T


# ---------- coordinates / structure constants ----------

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


def sc_table(basis, sq):
    tab = {}
    for i in range(len(basis)):
        for j in range(i + 1, len(basis)):
            br = basis[i] * basis[j] - basis[j] * basis[i]
            c = coords_in(basis, br, sq)
            if c is None:
                return None
            tab[(i, j)] = tuple(c)
    return tab


def bt_verify(L, Mo, mtable, sq):
    # direct structure constants of L must equal those transported through the model table
    k = len(L)
    if k == 0:
        return True
    if mtable is None:
        return False
    A = [coords_in(Mo, L[i], sq) for i in range(k)]
    if any(a is None for a in A):
        return False
    B = [coords_in(L, Mo[g], sq) for g in range(len(Mo))]
    if any(b is None for b in B):
        return False
    Ltab = sc_table(L, sq)
    if Ltab is None:
        return False
    km = len(Mo)
    for i in range(k):
        for j in range(i + 1, k):
            pv = [Integer(0)] * km
            for a in range(km):
                for b in range(km):
                    if a == b:
                        continue
                    coef = A[i][a] * A[j][b]
                    if coef == 0:
                        continue
                    if a < b:
                        f = mtable[(a, b)]
                        sgn = 1
                    else:
                        f = mtable[(b, a)]
                        sgn = -1
                    for g in range(km):
                        pv[g] = pv[g] + sgn * coef * f[g]
            pl = [Integer(0)] * k
            for g in range(km):
                if pv[g] != 0:
                    for t in range(k):
                        pl[t] = pl[t] + pv[g] * B[g][t]
            direct = Ltab[(i, j)]
            if any(pl[t] != direct[t] for t in range(k)):
                return False
    return True


def bt_str(L, sq):
    tab = sc_table(L, sq)
    if tab is None:
        return "bt=UNCLOSED"
    items = []
    for (i, j) in sorted(tab.keys()):
        c = tab[(i, j)]
        if any(v != 0 for v in c):
            items.append("({0},{1}):({2})".format(i, j, ",".join(str(v).replace(" ", "") for v in c)))
    return "bt=[" + ",".join(items) + "]"


# ---------- invariant block decomposition of S ----------

def evalmat(pl, S):
    n = S.rows
    M = zeros(n, n)
    for c in pl.all_coeffs():
        M = M * S + c * eye(n)
    return M


def dual_poly(f):
    return Poly(f.as_expr().subs(X_, -X_), X_).monic()


def poly_str(pl):
    return str(pl.as_expr()).replace(" ", "").replace("**", "^")


def blocks_of(S, eta, n, tag):
    mp = minpoly_poly(S)
    ok(squarefree(mp), "S diagonalizable over C (squarefree minpoly) " + tag)
    _c, facs = mp.factor_list()
    fl = []
    for (f, mult) in facs:
        ok(mult == 1, "factor multiplicity 1 " + tag)
        fm = Poly(f, X_).monic()
        ok(fm.degree() <= 2, "factor degree <= 2 " + tag)
        if fm.degree() == 2:
            ok(fm.discriminant() < 0, "quadratic factor irreducible over R " + tag)
        fl.append(fm)
    fl = sorted(fl, key=lambda t: (t.degree(), str(t.as_expr())))
    seen = set()
    recs = []
    for f in fl:
        key = str(f.as_expr())
        if key in seen:
            continue
        fd = dual_poly(f)
        kd = str(fd.as_expr())
        if kd == key:
            bp = f
            pair = None
            seen.add(key)
        else:
            ok(kd in [str(x.as_expr()) for x in fl], "dual factor present " + tag)
            bp = Poly(f.as_expr() * fd.as_expr(), X_).monic()
            pair = (f, fd)
            seen.add(key)
            seen.add(kd)
        B = evalmat(bp, S)
        K = (B ** n).nullspace()
        P = Matrix.hstack(*K)
        d = P.cols
        ok(P.rank() == d, "block basis full column rank " + tag)
        g = P.T * eta * P
        sig = cong_signature(g)
        ok(sig[2] == 0, "block Gram nondegenerate " + tag)
        ok(sig[0] + sig[1] + sig[2] == d, "block sig counts sum to dim " + tag)
        ok(sig[0] + sig[1] == g.rank(), "block sig pos+neg = rank " + tag)
        ptpinv = (P.T * P).inv()
        Ts = ptpinv * P.T * (S * P)
        ok((P * Ts - S * P).is_zero_matrix, "restriction of S exact " + tag)
        halves = None
        if pair is not None:
            h1 = evalmat(pair[0], S).nullspace()
            h2 = evalmat(pair[1], S).nullspace()
            C1 = []
            for u in h1:
                c = ptpinv * P.T * u
                if (P * c - u).is_zero_matrix:
                    C1.append(c)
            C2 = []
            for u in h2:
                c = ptpinv * P.T * u
                if (P * c - u).is_zero_matrix:
                    C2.append(c)
            if C1 and C2 and len(C1) == len(C2) and len(C1) + len(C2) == d:
                halves = (Matrix.hstack(*C1), Matrix.hstack(*C2))
        recs.append(dict(poly=bp, pstr=poly_str(bp), P=P, d=d, g=g,
                         sig=(sig[0], sig[1]), ptpinv=ptpinv, T=Ts, halves=halves))
    recs = sorted(recs, key=lambda r: (r["poly"].degree(), r["pstr"]))
    tot = sum(r["d"] for r in recs)
    ok(tot == n, "block dims sum to n " + tag)
    allP = Matrix.hstack(*[r["P"] for r in recs])
    ok(allP.rank() == n, "blocks span R^n " + tag)
    for i in range(len(recs)):
        for j in range(i + 1, len(recs)):
            ok(Matrix.hstack(recs[i]["P"], recs[j]["P"]).rank() == recs[i]["d"] + recs[j]["d"],
               "blocks pairwise zero intersection " + tag)
            ok((recs[i]["P"].T * eta * recs[j]["P"]).is_zero_matrix,
               "blocks pairwise eta-orthogonal " + tag)
    return recs


def preserves(M, P):
    return Matrix.hstack(P, M * P).rank() == P.rank()


def restrict(M, rec, tag):
    R = rec["ptpinv"] * rec["P"].T * (M * rec["P"])
    ok((rec["P"] * R - M * rec["P"]).is_zero_matrix, "restriction exact " + tag)
    return R


# ---------- model builders (full menu) ----------

def end_basis(d):
    out = []
    for i in range(d):
        for j in range(d):
            E = zeros(d, d)
            E[i, j] = Integer(1)
            out.append(E)
    return out


def model_from_constraints(d, cons):
    eb = end_basis(d)
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


def model_so_g(g):
    d = g.rows
    return model_from_constraints(d, [lambda E: E.T * g + g * E])


def model_u_g(g, J):
    d = g.rows
    return model_from_constraints(d, [lambda E: E.T * g + g * E,
                                      lambda E: E * J - J * E])


def model_sp_w(w):
    d = w.rows
    return model_from_constraints(d, [lambda E: E.T * w + w * E])


def hyp_form(m):
    H = zeros(2 * m, 2 * m)
    for i in range(m):
        H[i, m + i] = Integer(1)
        H[m + i, i] = Integer(1)
    return H


def gl_H_basis(m):
    out = []
    for a in range(m):
        for b in range(m):
            H = zeros(2 * m, 2 * m)
            H[a, b] = Integer(1)
            H[m + b, m + a] = Integer(-1)
            out.append(H)
    return out


def hyper_from_halves(Cp, Cm, g):
    m = Cp.cols
    if not (Cp.T * g * Cp).is_zero_matrix:
        return None
    if not (Cm.T * g * Cm).is_zero_matrix:
        return None
    B = Cp.T * g * Cm
    if B.det() == 0:
        return None
    Q = Matrix.hstack(Cp, Cm * B.inv())
    if not (Q.T * g * Q - hyp_form(m)).is_zero_matrix:
        return None
    return Q


def hyper_from_split(g):
    d = g.rows
    if d % 2 != 0:
        return None
    m = d // 2
    D, T = cong_diag(g)
    for i in range(d):
        for j in range(d):
            if i != j and D[i, j] != 0:
                return None
    pos = [i for i in range(d) if D[i, i] > 0]
    neg = [i for i in range(d) if D[i, i] < 0]
    if len(pos) != m or len(neg) != m:
        return None
    used = set()
    pairs = []
    for i in pos:
        found = None
        for j in neg:
            if j in used:
                continue
            t = ssqrt(-D[i, i] / D[j, j])
            if t.is_rational:
                found = (j, t)
                break
        if found is None:
            return None
        used.add(found[0])
        pairs.append((i, found[0], found[1]))
    xs = []
    ys = []
    for (i, j, t) in pairs:
        u = T[:, i]
        v = t * T[:, j]
        xs.append(u + v)
        ys.append((u - v) / (2 * D[i, i]))
    Q = Matrix.hstack(*(xs + ys))
    if not (Q.T * g * Q - hyp_form(m)).is_zero_matrix:
        return None
    return Q


def try_gl_model(rec):
    d = rec["d"]
    if d % 2 != 0:
        return None
    m = d // 2
    if rec["sig"] != (m, m):
        return None
    Q = None
    if rec["halves"] is not None:
        Q = hyper_from_halves(rec["halves"][0], rec["halves"][1], rec["g"])
    if Q is None:
        Q = hyper_from_split(rec["g"])
    if Q is None:
        return None
    Qi = Q.inv()
    return [Q * H * Qi for H in gl_H_basis(m)]


def try_u_model(rec):
    d = rec["d"]
    T2 = rec["T"] * rec["T"]
    s = T2[0, 0]
    if not (T2 - s * eye(d)).is_zero_matrix:
        return None
    if not (s < 0):
        return None
    a = ssqrt(-s)
    if not a.is_rational:
        return None
    J = rec["T"] / a
    if not (J * J + eye(d)).is_zero_matrix:
        return None
    return model_u_g(rec["g"], J)


def antisym_basis(d):
    out = []
    for i in range(d):
        for j in range(i + 1, d):
            A = zeros(d, d)
            A[i, j] = Integer(1)
            A[j, i] = Integer(-1)
            out.append(A)
    return out


def inv_antisym_forms(L, d):
    ab = antisym_basis(d)
    if not ab:
        return []
    if not L:
        return ab
    cols = []
    for A in ab:
        parts = [flat(Xm.T * A + A * Xm) for Xm in L]
        cols.append(Matrix.vstack(*parts))
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        W = zeros(d, d)
        for k in range(len(ab)):
            if v[k, 0] != 0:
                W = W + v[k, 0] * ab[k]
        out.append(W)
    return out


def find_nondeg_w(Ws, d):
    for W in Ws:
        if W.rank() == d:
            return W
    kk = len(Ws)
    if kk >= 2 and kk <= 8:
        for cf in itertools.product((0, 1, 2), repeat=kk):
            if all(c == 0 for c in cf):
                continue
            W = zeros(d, d)
            for k in range(kk):
                if cf[k]:
                    W = W + Integer(cf[k]) * Ws[k]
            if W.rank() == d:
                return W
    return None


def try_sp_model(rec, L):
    d = rec["d"]
    Ws = inv_antisym_forms(L, d)
    w = find_nondeg_w(Ws, d)
    if w is None:
        return None
    return model_sp_w(w)


# ---------- identification against the full menu ----------

COUNTS = {"gl": 0, "u": 0, "so": 0, "sp": 0, "OTHER": 0, "sub-of": 0}


def identify(rec, L, tag):
    d = rec["d"]
    p0, q0 = rec["sig"]
    dimL = span_rank(L, d)
    _dcc, closed = derived_and_closed(L, d)
    ok(closed, "image bracket-closed " + tag)
    cand = []
    glm = try_gl_model(rec)
    if glm is not None:
        ok(span_rank(glm, d) == (d // 2) ** 2, "gl model dim " + tag)
        cand.append(("gl({0},R)".format(d // 2), glm))
    um = try_u_model(rec)
    if um is not None:
        ok(span_rank(um, d) == (d // 2) ** 2, "u model dim " + tag)
        cand.append(("u({0})".format(d // 2), um))
    som = model_so_g(rec["g"])
    ok(span_rank(som, d) == d * (d - 1) // 2, "so model dim " + tag)
    ok(contains(som, L, d), "image inside so(g) " + tag)
    cand.append(("so({0},{1})".format(p0, q0), som))
    spm = try_sp_model(rec, L)
    if spm is not None:
        ok(span_rank(spm, d) == (d // 2) * (d + 1), "sp model dim " + tag)
        cand.append(("sp({0},R)".format(d), spm))
    matches = []
    for (name, mo) in cand:
        if same_span(L, mo, d):
            matches.append((name, mo))
    if matches:
        for (name, mo) in matches:
            mob = span_basis(mo, d)
            tab = sc_table(mob, d)
            ok(bt_verify(L, mob, tab, d), "bracket table vs model " + name + " " + tag)
            COUNTS[name.split("(")[0]] += 1
        return "match=" + "&".join(nm for (nm, _m) in matches) + ":dim" + str(dimL)
    containers = []
    for idx, (name, mo) in enumerate(cand):
        if contains(mo, L, d):
            containers.append((span_rank(mo, d), idx, name))
    if containers:
        containers.sort()
        COUNTS["sub-of"] += 1
        return "sub-of=" + containers[0][2] + ":dim" + str(dimL) + ":" + bt_str(span_basis(L, d), d)
    COUNTS["OTHER"] += 1
    diags = ",".join(name + ":int" + str(inter_dim(L, mo, d)) for (name, mo) in cand)
    return "OTHER:dim" + str(dimL) + "(" + diags + ")"


# ---------- shared per-stratum pipeline ----------

def blocks_str(recs):
    return "[" + ",".join("({0},{1},({2},{3}))".format(r["pstr"], r["d"], r["sig"][0], r["sig"][1])
                          for r in recs) + "]"


def pipeline(A, S, p, q, eta, bas, tag):
    n = p + q
    ok(is_so(A, eta), "A in so " + tag)
    recs = blocks_of(S, eta, n, tag)
    cb = centralizer(A, bas)
    ok(all(is_so(M, eta) for M in cb), "c-basis in so " + tag)
    ok(commutes_all(cb, A), "c-basis commutes " + tag)
    dimcc, closed = derived_and_closed(cb, n)
    ok(closed, "c closed under bracket " + tag)
    viol = 0
    for M in cb:
        for r in recs:
            if not preserves(M, r["P"]):
                viol += 1
    pres = "ALL" if viol == 0 else "VIOL:" + str(viol)
    ok(viol == 0, "all c elements preserve all blocks " + tag)
    Rls = []
    for r in recs:
        Rls.append([restrict(M, r, tag) for M in cb])
    if cb:
        rho_cols = []
        for i in range(len(cb)):
            rho_cols.append(Matrix.vstack(*[flat(Rls[k][i]) for k in range(len(recs))]))
        ker_rho = len(cb) - Matrix.hstack(*rho_cols).rank()
    else:
        ker_rho = 0
    images = []
    sumdim = 0
    for k, r in enumerate(recs):
        L = span_basis(Rls[k], r["d"])
        sumdim += len(L)
        res = identify(r, L, tag + " f=" + r["pstr"])
        images.append(r["pstr"] + ":" + res)
    audit_ok = (len(cb) == sumdim + ker_rho)
    ok(audit_ok, "dimension audit " + tag)
    return dict(blocks=blocks_str(recs), dimc=len(cb), dimcc=dimcc, pres=pres,
                ker_rho=ker_rho, images="[" + "; ".join(images) + "]",
                audit="OK" if audit_ok else "FAIL({0} vs {1}+{2})".format(len(cb), sumdim, ker_rho))


# ---------- strata (semisimple S) ----------

def strata(p, q):
    out = []
    if p >= 2:
        out.append(("S1", [("R+", 0, 1, 1)]))
    if q >= 1:
        out.append(("S2", [("B", 0, p, 1)]))
    if p >= 4:
        out.append(("S3", [("R+", 0, 1, 1), ("R+", 2, 3, 2)]))
        out.append(("S4", [("R+", 0, 1, 1), ("R+", 2, 3, 1)]))
    if p >= 2 and q >= 2:
        out.append(("S5", [("B", 0, p, 1), ("B", 1, p + 1, 2)]))
        out.append(("S6", [("B", 0, p, 1), ("B", 1, p + 1, 1)]))
    if p >= 3 and q >= 1:
        out.append(("S7", [("R+", 0, 1, 1), ("B", 2, p, 1)]))
    return out


def sig_list():
    out = []
    for n in (4, 5, 6):
        for q in range(0, n // 2 + 1):
            p = n - q
            if p >= q:
                out.append((p, q))
    return out


ROWS_STRATA = [0]
ROWS_JORDAN = [0]

print("--- strata (semisimple S) ---")
for (p, q) in sig_list():
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    ok(len(bas) == n * (n - 1) // 2, "basis size so({0},{1})".format(p, q))
    ok(all(is_so(B, eta) for B in bas), "basis membership so({0},{1})".format(p, q))
    ok(span_rank(bas, n) == n * (n - 1) // 2, "basis independent so({0},{1})".format(p, q))
    for (sid, spec) in strata(p, q):
        S = zeros(n, n)
        for (t, i, j, par) in spec:
            S = S + block_gen(n, eta, i, j, par)
        sstr = ",".join("{0}({1},{2})x{3}".format(t, i, j, par) for (t, i, j, par) in spec)
        tag = sid + " (" + str(p) + "," + str(q) + ")"
        r = pipeline(S, S, p, q, eta, bas, tag)
        ROWS_STRATA[0] += 1
        print("SIG=({0},{1}) | S={2}[{3}] | BLOCKS={4} | dimc={5} | dimcc={6} | preserved={7} | ker_rho={8} | IMAGES={9} | audit={10}".format(
            p, q, sid, sstr, r["blocks"], r["dimc"], r["dimcc"], r["pres"], r["ker_rho"], r["images"], r["audit"]))

# ---------- Jordan-reduction strata (A = S + N) ----------

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
    dens = []
    for i in range(d):
        for j in range(d):
            dens.append(int(G[i, j].q))
    Lc = 1
    for de in dens:
        Lc = _ilcm(Lc, de)
    Gi = [[int(G[i, j] * Lc) for j in range(d)] for i in range(d)]
    vals = (-2, -1, 0, 1, 2)
    for c in itertools.product(vals, repeat=d):
        if all(x == 0 for x in c):
            continue
        s = 0
        for i in range(d):
            ci = c[i]
            if ci:
                row = Gi[i]
                s += ci * sum(row[j] * c[j] for j in range(d))
        if s != 0:  # tr(X*X) != 0 => not nilpotent
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


def d1_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 1 and len(Km) >= 1 and len(K) >= 3:
        a, b = Kp[0], Km[0]
        rest = [c for c in K if c not in (a, b)]
        c = rest[0]
        return unit(n, a) + unit(n, b), unit(n, c)
    return None


def d0_on(K, p, n):
    Kp = [a for a in K if a < p]
    Km = [a for a in K if a >= p]
    if len(Kp) >= 2 and len(Km) >= 2:
        x = unit(n, Kp[0]) + unit(n, Km[0])
        y = unit(n, Kp[1]) + unit(n, Km[1])
        return x, y
    return None


print("--- jordan (A = S + N) ---")

JCASES = []

# J1: (2,2) S = two equal B, N from deterministic scan in span of c(S)
_p, _q = 2, 2
_n = 4
_eta = make_eta(_p, _q)
_bas = so_basis(_n, _eta)
_S = block_gen(_n, _eta, 0, 2, 1) + block_gen(_n, _eta, 1, 3, 1)
_cb = centralizer(_S, _bas)
_Xn = jordan_scan(_cb, q_gram(_cb), _n)
ok(_Xn is not None, "J1 scan finds nilpotent")
JCASES.append(("J1", _p, _q, _S, "B(0,2)x1,B(1,3)x1", _Xn, "scan"))

# J2: (4,2) S = single R on (0,1), N = D0 wedge on complement
_p, _q = 4, 2
_n = 6
_eta = make_eta(_p, _q)
_S = block_gen(_n, _eta, 0, 1, 1)
_K = [2, 3, 4, 5]
_x, _y = d0_on(_K, _p, _n)
ok(gram2(_x, _y, _eta).is_zero_matrix, "J2 Gram rank 0")
JCASES.append(("J2", _p, _q, _S, "R+(0,1)x1", wedge(_x, _y, _eta), "D0"))

# J3: (5,1) S = single R on (0,1), N = D1 wedge on complement
_p, _q = 5, 1
_n = 6
_eta = make_eta(_p, _q)
_S = block_gen(_n, _eta, 0, 1, 1)
_K = [2, 3, 4, 5]
_x, _y = d1_on(_K, _p, _n)
ok(gram2(_x, _y, _eta).rank() == 1, "J3 Gram rank 1")
JCASES.append(("J3", _p, _q, _S, "R+(0,1)x1", wedge(_x, _y, _eta), "D1"))

# J4: (3,3) S = single B on (0,3), N = D0 wedge on complement
_p, _q = 3, 3
_n = 6
_eta = make_eta(_p, _q)
_S = block_gen(_n, _eta, 0, _p, 1)
_K = [1, 2, 4, 5]
_x, _y = d0_on(_K, _p, _n)
ok(gram2(_x, _y, _eta).is_zero_matrix, "J4 Gram rank 0")
JCASES.append(("J4", _p, _q, _S, "B(0,3)x1", wedge(_x, _y, _eta), "D0"))

for (jid, p, q, S, sstr, N, ncls) in JCASES:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    tag = jid + " (" + str(p) + "," + str(q) + ")"
    ok(is_so(N, eta), "N in so " + tag)
    ok(is_nilp(N), "N nilpotent " + tag)
    ok(not N.is_zero_matrix, "N nonzero " + tag)
    ok((S * N - N * S).is_zero_matrix, "[S,N]=0 " + tag)
    A = S + N
    r = pipeline(A, S, p, q, eta, bas, tag)
    ROWS_JORDAN[0] += 1
    print("SIG=({0},{1}) | A=S+N {2}[S={3};N={4}] | BLOCKS={5} | dimc={6} | dimcc={7} | preserved={8} | ker_rho={9} | IMAGES={10} | audit={11}".format(
        p, q, jid, sstr, ncls, r["blocks"], r["dimc"], r["dimcc"], r["pres"], r["ker_rho"], r["images"], r["audit"]))

# ---------- mutants (real code paths; each must be CAUGHT) ----------

print("--- mutants ---")
mut_ok = True

# m1: corrupt one structure constant of a computed model bracket table -> MISMATCH
_g3 = eye(3)
_mod3 = span_basis(model_so_g(_g3), 3)
_tab3 = sc_table(_mod3, 3)
ok(_tab3 is not None, "m1 baseline: model table exists")
ok(bt_verify(_mod3, _mod3, _tab3, 3), "m1 baseline: uncorrupted table verifies")
_tab_bad = dict(_tab3)
_k0 = sorted(_tab_bad.keys())[0]
_row = list(_tab_bad[_k0])
_row[0] = _row[0] + 1
_tab_bad[_k0] = tuple(_row)
if not bt_verify(_mod3, _mod3, _tab_bad, 3):
    print("MUTANT m1: CAUGHT (bracket-table MISMATCH reported)")
else:
    print("MUTANT m1: NOT CAUGHT")
    mut_ok = False

# m2: inject a matrix that does NOT preserve a block into a c-basis -> flagged
p, q = 2, 2
n = 4
eta = make_eta(p, q)
bas = so_basis(n, eta)
S_m2 = block_gen(n, eta, 0, 1, 1)
recs_m2 = blocks_of(S_m2, eta, n, "m2")
cb_m2 = centralizer(S_m2, bas)
bad_m2 = wedge(unit(n, 0), unit(n, 2), eta)
tampered = cb_m2 + [bad_m2]
viol_m2 = 0
for M in tampered:
    for r in recs_m2:
        if not preserves(M, r["P"]):
            viol_m2 += 1
ok(all(preserves(M, r["P"]) for M in cb_m2 for r in recs_m2), "m2 baseline: true c-basis preserves")
if viol_m2 > 0:
    print("MUTANT m2: CAUGHT (block-preservation flagged {0} violations)".format(viol_m2))
else:
    print("MUTANT m2: NOT CAUGHT")
    mut_ok = False

# m3: kernel-block image tested against gl of wrong dim -> subspace equality must reject
p, q = 3, 3
n = 6
eta = make_eta(p, q)
bas = so_basis(n, eta)
S_m3 = block_gen(n, eta, 0, p, 1)
recs_m3 = blocks_of(S_m3, eta, n, "m3")
cb_m3 = centralizer(S_m3, bas)
rk = None
for r in recs_m3:
    if r["pstr"] == "x":
        rk = r
ok(rk is not None and rk["d"] == 4 and rk["sig"] == (2, 2), "m3 baseline: kernel block (2,2) dim 4")
L_m3 = span_basis([restrict(M, rk, "m3") for M in cb_m3], rk["d"])
ok(span_rank(L_m3, rk["d"]) == 6, "m3 baseline: kernel image is so(2,2) dim 6")
glm_m3 = try_gl_model(rk)
ok(glm_m3 is not None, "m3 baseline: gl(2) model constructible on split block")
if not same_span(L_m3, glm_m3, rk["d"]):
    print("MUTANT m3: CAUGHT (gl(2,R) equality rejected for so(2,2)-image)")
else:
    print("MUTANT m3: NOT CAUGHT")
    mut_ok = False

# m4 (survivability control): each tester PASSES on a hand-built canonical instance
m4_ok = True

# m4-gl: hand-built gl(2) inside so(2,2) in hyperbolic coordinates
_gh = hyp_form(2)
_handL = gl_H_basis(2)
_rec_gl = dict(poly=None, pstr="hand", P=eye(4), d=4, g=_gh, sig=(2, 2),
               ptpinv=eye(4), T=zeros(4, 4), halves=None)
_glm = try_gl_model(_rec_gl)
if _glm is not None and same_span(_handL, _glm, 4):
    print("MUTANT m4: gl tester PASS on canonical gl(2,R)")
else:
    print("MUTANT m4: gl tester FAIL")
    m4_ok = False

# m4-u: hand-built u(1) inside so(2)
_J2 = Matrix([[0, -1], [1, 0]])
_rec_u = dict(poly=None, pstr="hand", P=eye(2), d=2, g=eye(2), sig=(2, 0),
              ptpinv=eye(2), T=_J2, halves=None)
_um = try_u_model(_rec_u)
if _um is not None and same_span([_J2], _um, 2):
    print("MUTANT m4: u tester PASS on canonical u(1)")
else:
    print("MUTANT m4: u tester FAIL")
    m4_ok = False

# m4-so: hand-built so(2,1)
_g21 = make_eta(2, 1)
_hand_so = so_basis(3, _g21)
_som = model_so_g(_g21)
if same_span(_hand_so, _som, 3) and bt_verify(_hand_so, span_basis(_som, 3), sc_table(span_basis(_som, 3), 3), 3):
    print("MUTANT m4: so tester PASS on canonical so(2,1)")
else:
    print("MUTANT m4: so tester FAIL")
    m4_ok = False

# m4-sp: hand-built sp(2,R) = traceless 2x2
_hand_sp = [Matrix([[1, 0], [0, -1]]), Matrix([[0, 1], [0, 0]]), Matrix([[0, 0], [1, 0]])]
_Ws = inv_antisym_forms(_hand_sp, 2)
_w = find_nondeg_w(_Ws, 2)
if _w is not None:
    _spm = model_sp_w(_w)
    if same_span(_hand_sp, _spm, 2) and bt_verify(_hand_sp, span_basis(_spm, 2), sc_table(span_basis(_spm, 2), 2), 2):
        print("MUTANT m4: sp tester PASS on canonical sp(2,R)")
    else:
        print("MUTANT m4: sp tester FAIL")
        m4_ok = False
else:
    print("MUTANT m4: sp tester FAIL (no invariant form found)")
    m4_ok = False

if m4_ok:
    print("MUTANT m4: CAUGHT (all four testers pass on their own canonical models)")
else:
    print("MUTANT m4: NOT CAUGHT")
    mut_ok = False

# ---------- summary ----------

print("SUMMARY: rows_strata={0} | rows_jordan={1}".format(ROWS_STRATA[0], ROWS_JORDAN[0]))
print("SUMMARY: matches gl={0} | u={1} | so={2} | sp={3} | OTHER={4} | sub-of={5}".format(
    COUNTS["gl"], COUNTS["u"], COUNTS["so"], COUNTS["sp"], COUNTS["OTHER"], COUNTS["sub-of"]))
print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

# ---------- self-scan of source and produced log text for banned substrings ----------

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
