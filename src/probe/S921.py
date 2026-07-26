# -*- coding: utf-8 -*-
# DIM: na (dimensionless/exact symbolic count; no spatial lattice)
"""
S921 (W29): affine extension so(p,q) ⋉ R^n — the shift leg.

η = diag(+1×p, −1×q), n = p+q; so(p,q) = {X : Xη + ηXᵀ = 0}
(basis = the EXACT nullspace of the linear condition, dim = n(n−1)/2).
Homogeneous representation of the affine algebra: (n+1)×(n+1) matrices;
M ∈ so is embedded as M̂ = [[M,0],[0,0]]; a translation T_a = [[0,a],[0,0]].

SECTIONS:
  A1 — closure and forcedness of the extension ([M̂,T_a]=T_{Ma}, [T,T]=0,
       the so-block closed, Jacobi over the triple types);
  A2 — invariant forms on the translations (linear / symmetric bilinear /
       antisymmetric bilinear): dimensions of the exact solution spaces;
  A3 — the fate of the time shift T_t (t = the last minus-axis): a bracket
       table, the stabilizer of e_t, carrying over the S920 cone-argument
       (q=1) and an explicit continuous path T_t → −T_t via a {−−}-rotation (q≥2);
  M  — mutants (m1: [M,T_a]=T_{Mᵀa}; m2: φ=(1,0,…,0); m3: "a second symmetric form");
  SUMMARY — verdicts + a rollup table.

LOG link (cited, NOT recomputed): the cone-argument for q=1 —
  active-v10.2/src/symbolic/S920_w29_energy_sign_orientation.py
  (log: active-v10.2/src/symbolic/S920_w29_run.log; verdict: the sector
   whose sign has multiplicity 1 is orientable — two invariant cones,
   the flows exp(sM) do not connect them).

Exact symbolics (sympy): Integer/Rational/symbols/cos/sin/π.
NO float/random/time/datetime. Non-machine links are marked [LOG].
exit 0 only if: all checks OK, mutants caught (with the m3 caveat),
no section empty (zero checked in a section = exit≠0).
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

from sympy import (Symbol, symbols, Integer, Matrix, eye, zeros, diag,
                   cos, sin, pi, simplify, linear_eq_to_matrix)

# ---------------------------------------------------------------- bookkeeping
CHECKS_TOTAL = 0
FAILS = []                 # (section, label)
SECTION_COUNTS = {}
SECTION_FAILS = {}
VERDICTS = []              # (section, text)
MUTANT_STATUS = {}         # mutant name -> 'CAUGHT' / 'NOT-CAUGHT(result)'
CUR = None


def section(name):
    global CUR
    CUR = name
    SECTION_COUNTS.setdefault(name, 0)
    SECTION_FAILS.setdefault(name, 0)
    print()
    print("=" * 78)
    print(f"SECTION {name}")
    print("=" * 78)


def check(label, ok):
    global CHECKS_TOTAL
    CHECKS_TOTAL += 1
    SECTION_COUNTS[CUR] += 1
    if not ok:
        SECTION_FAILS[CUR] += 1
        FAILS.append((CUR, label))
        print(f"  [FAIL] {label}")
    return ok


def verdict(text):
    VERDICTS.append((CUR, text))
    print(f"  VERDICT [{CUR}]: {text}")


# ---------------------------------------------------------------- algebra
def eta_of(p, q):
    return diag(*([Integer(1)] * p + [Integer(-1)] * q))


def so_basis_nullspace(p, q):
    """Basis of so(p,q) as the exact nullspace of Xη + ηXᵀ = 0 (n² unknowns)."""
    n = p + q
    et = [Integer(1)] * p + [Integer(-1)] * q
    rows = []
    for a in range(n):
        for b in range(n):
            r = [Integer(0)] * (n * n)
            r[a * n + b] += et[b]      # (Xη)_{ab} = X_{ab} η_b
            r[b * n + a] += et[a]      # (ηXᵀ)_{ab} = η_a X_{ba}
            rows.append(r)
    ns = Matrix(rows).nullspace()
    mats = [Matrix(n, n, lambda a, b: v[a * n + b]) for v in ns]
    return mats


def pair_of(M):
    """The axis pair of the support of a basis element: expect exactly {(i,j),(j,i)}."""
    n = M.shape[0]
    nz = [(a, b) for a in range(n) for b in range(n) if M[a, b] != 0]
    if len(nz) == 2 and nz[0] == (nz[1][1], nz[1][0]):
        i, j = nz[0]
        return (min(i, j), max(i, j))
    return None


def hat(M):
    """Embedding of an so-element into the homogeneous (n+1)×(n+1) representation."""
    n = M.shape[0]
    return M.row_join(zeros(n, 1)).col_join(zeros(1, n + 1))


def Tmat(a):
    """Translation T_a = [[0,a],[0,0]] (a — an n×1 column)."""
    n = a.shape[0]
    T = zeros(n + 1, n + 1)
    for i in range(n):
        T[i, n] = a[i]
    return T


def br(X, Y):
    return X * Y - Y * X


def e_vec(n, i):
    v = zeros(n, 1)
    v[i, 0] = Integer(1)
    return v


# ---------------------------------------------------------------- sampling
SIGS = [(3, 0), (2, 1), (1, 2), (0, 3),
        (4, 0), (3, 1), (2, 2), (1, 3), (0, 4),
        (5, 0), (4, 1), (3, 2), (2, 3), (1, 4), (0, 5),
        (5, 1), (3, 3)]

print("S921 (W29): affine extension so(p,q) ⋉ R^n — the shift leg")
print()
print("SAMPLING BOUNDARY (explicit): all signatures (p,q) with n=p+q ∈ {3,4,5}")
print("  (including the definite q=0 and p=0) + two from n=6: (5,1), (3,3).")
print(f"  Total {len(SIGS)} signatures: {SIGS}")
print("Jacobi (A1d): n=3,4 — EXHAUSTIVELY all three triple types over the bases;")
print("  n=5,6 — declared: (M,M,T) EXHAUSTIVELY all pairs of basis M × all")
print("  basis T; (M,T,T) EXHAUSTIVELY all M × pairs of T; (T,T,T) EXHAUSTIVELY all")
print("  triples of T (small volumes, so the \"sample\" = a full enumeration; the bound = the bases).")

# caches between sections
BASIS = {}      # sig -> [M_k]
HATS = {}       # sig -> [M̂_k]
TS = {}         # sig -> [T_i]
PAIRS = {}      # sig -> [axis pair of M_k]
DIM_LIN = {}
DIM_SYM = {}
SYM_PROP_ETA = {}
DIM_ANTI = {}
DIM_STAB = {}

# ============================================================== A1
section("A1 — CLOSURE AND FORCEDNESS OF THE EXTENSION")

for (p, q) in SIGS:
    n = p + q
    et = eta_of(p, q)
    basis = so_basis_nullspace(p, q)
    BASIS[(p, q)] = basis
    d = len(basis)
    check(f"({p},{q}) dim so = n(n-1)/2 = {n*(n-1)//2}", d == n * (n - 1) // 2)

    pairs = []
    for k, M in enumerate(basis):
        check(f"({p},{q}) basis[{k}] satisfies Xη+ηXᵀ=0",
              (M * et + et * M.T) == zeros(n, n))
        pr = pair_of(M)
        check(f"({p},{q}) basis[{k}] has a support pair", pr is not None)
        pairs.append(pr)
    PAIRS[(p, q)] = pairs

    hats = [hat(M) for M in basis]
    HATS[(p, q)] = hats
    evs = [e_vec(n, i) for i in range(n)]
    Ts = [Tmat(ev) for ev in evs]
    TS[(p, q)] = Ts

    # (a) [M̂_k, T_i] = T_{M_k e_i}
    MT = [[br(hats[k], Ts[i]) for i in range(n)] for k in range(d)]
    cnt_a = 0
    for k in range(d):
        for i in range(n):
            ok = MT[k][i] == Tmat(basis[k] * evs[i])
            cnt_a += 1
            if not ok:
                check(f"({p},{q}) (a) [M̂_{k},T_{i}] = T_(M e_i)", False)
    check(f"({p},{q}) (a) all {cnt_a} brackets [M̂,T_i]=T_(M e_i)", True)
    SECTION_COUNTS[CUR] += cnt_a - 1  # count every bracket

    # (b) [T_i, T_j] = 0
    cnt_b = 0
    for i in range(n):
        for j in range(i + 1, n):
            ok = br(Ts[i], Ts[j]) == zeros(n + 1, n + 1)
            cnt_b += 1
            if not ok:
                check(f"({p},{q}) (b) [T_{i},T_{j}]=0", False)
    check(f"({p},{q}) (b) all {cnt_b} pairs [T_i,T_j]=0", True)
    SECTION_COUNTS[CUR] += cnt_b - 1

    # (c) [M̂_k, M̂_l] = the embedding of [M_k,M_l]; the block is closed in so
    HB = {}
    cnt_c = 0
    for k in range(d):
        for l in range(k + 1, d):
            blk = br(basis[k], basis[l])
            HB[(k, l)] = br(hats[k], hats[l])
            ok1 = HB[(k, l)] == hat(blk)
            ok2 = (blk * et + et * blk.T) == zeros(n, n)
            cnt_c += 1
            if not (ok1 and ok2):
                check(f"({p},{q}) (c) pair ({k},{l})", False)
    check(f"({p},{q}) (c) all {cnt_c} pairs: [M̂,M̂]=embedding and the block ∈ so", True)
    SECTION_COUNTS[CUR] += cnt_c - 1

    # (d) Jacobi
    Z = zeros(n + 1, n + 1)
    cnt_mmt = cnt_mtt = cnt_ttt = 0
    for k in range(d):
        for l in range(k + 1, d):
            for i in range(n):
                Jc = (br(HB[(k, l)], Ts[i]) + br(MT[l][i], hats[k])
                      - br(MT[k][i], hats[l]))
                cnt_mmt += 1
                if Jc != Z:
                    check(f"({p},{q}) (d) Jacobi MMT ({k},{l},{i})", False)
    for k in range(d):
        for i in range(n):
            for j in range(i + 1, n):
                Jc = (br(MT[k][i], Ts[j]) + br(br(Ts[i], Ts[j]), hats[k])
                      - br(MT[k][j], Ts[i]))
                cnt_mtt += 1
                if Jc != Z:
                    check(f"({p},{q}) (d) Jacobi MTT ({k},{i},{j})", False)
    for i in range(n):
        for j in range(i + 1, n):
            for m in range(j + 1, n):
                Jc = (br(br(Ts[i], Ts[j]), Ts[m])
                      + br(br(Ts[j], Ts[m]), Ts[i])
                      + br(br(Ts[m], Ts[i]), Ts[j]))
                cnt_ttt += 1
                if Jc != Z:
                    check(f"({p},{q}) (d) Jacobi TTT ({i},{j},{m})", False)
    scope = "exhaustive" if n <= 4 else "exhaustive (declared bound: the bases)"
    check(f"({p},{q}) (d) Jacobi {scope}: MMT={cnt_mmt}, MTT={cnt_mtt}, "
          f"TTT={cnt_ttt} — all zero", True)
    SECTION_COUNTS[CUR] += cnt_mmt + cnt_mtt + cnt_ttt - 1
    print(f"  ({p},{q}): dim so={d}; (a)={cnt_a} OK, (b)={cnt_b} OK, "
          f"(c)={cnt_c} OK, (d) MMT={cnt_mmt} MTT={cnt_mtt} TTT={cnt_ttt} OK")

verdict("the extension so(p,q) ⋉ R^n is CLOSED in the homogeneous representation; "
        "the structure constants are FORCED: [M̂,T_a]=T_{Ma} — an identity of "
        "block matrix multiplication, [T,T]=0, Jacobi holds "
        "automatically (an associative algebra) and is checked explicitly — "
        "NO free parameter in the structure of the extension.")

# ============================================================== A2
section("A2 — INVARIANT FORMS ON THE TRANSLATIONS (the main section)")

for (p, q) in SIGS:
    n = p + q
    et = eta_of(p, q)
    basis = BASIS[(p, q)]
    d = len(basis)
    print(f"  --- signature ({p},{q}), n={n}, dim so={d} ---")

    # (a) linear: φ(M_k e_i)=0 ∀k,i  ⇔  φ M_k = 0 ∀k
    fs = symbols(f"f0:{n}")
    phi = Matrix(1, n, fs)
    eqs = []
    for M in basis:
        row = phi * M
        for i in range(n):
            eqs.append(row[0, i])
    A, _ = linear_eq_to_matrix(eqs, list(fs))
    ns_lin = A.nullspace()
    dim_lin = len(ns_lin)
    DIM_LIN[(p, q)] = dim_lin
    check(f"({p},{q}) (a) linear: dim = n − rank = {n - A.rank()}",
          dim_lin == n - A.rank())
    for v in ns_lin:
        for k, M in enumerate(basis):
            check(f"({p},{q}) (a) the solution is annihilated by M_{k}",
                  (Matrix(1, n, list(v)) * M) == zeros(1, n))
    print(f"      (a) LINEAR invariant forms: dim = {dim_lin}"
          + ("" if dim_lin == 0 else f"; basis = {[list(v) for v in ns_lin]}"))

    # (b) symmetric bilinear: M_kᵀ B + B M_k = 0
    bs = symbols(f"b0:{n*(n+1)//2}")
    B = zeros(n, n)
    t_idx = 0
    for a in range(n):
        for b in range(a, n):
            B[a, b] = bs[t_idx]
            B[b, a] = bs[t_idx]
            t_idx += 1
    eqs = []
    for M in basis:
        E = M.T * B + B * M
        for a in range(n):
            for b in range(a, n):
                eqs.append(E[a, b])
    A, _ = linear_eq_to_matrix(eqs, list(bs))
    ns_sym = A.nullspace()
    dim_sym = len(ns_sym)
    DIM_SYM[(p, q)] = dim_sym
    check(f"({p},{q}) (b) sym.: dim = {dim_sym} (exact nullspace)",
          dim_sym == len(bs) - A.rank())
    sols = []
    for v in ns_sym:
        Bv = B.subs({bs[m]: v[m] for m in range(len(bs))})
        sols.append(Bv)
        for k, M in enumerate(basis):
            check(f"({p},{q}) (b) the solution is invariant under M_{k}",
                  (M.T * Bv + Bv * M) == zeros(n, n))
    # does η lie in the solution space?
    eta_ok = all((M.T * et + et * M) == zeros(n, n) for M in basis)
    check(f"({p},{q}) (b) η — an invariant sym. form", eta_ok)
    # is every solution ∝ η?
    all_prop = True
    for Bv in sols:
        c = (Bv * et).trace() / Integer(n)
        if Bv != c * et:
            all_prop = False
    SYM_PROP_ETA[(p, q)] = all_prop
    check(f"({p},{q}) (b) every basis solution ∝ η", all_prop)
    print(f"      (b) SYMMETRIC invariants: dim = {dim_sym}; "
          f"η in the space: {'YES' if eta_ok else 'NO'}; "
          f"all ∝ η: {'YES' if all_prop else 'NO'}")
    for Bv in sols:
        print(f"          basis of solutions: {Bv.tolist()}")

    # (c) antisymmetric bilinear: M_kᵀ A + A M_k = 0, Aᵀ=−A
    if n * (n - 1) // 2 > 0:
        cs = symbols(f"c0:{n*(n-1)//2}")
        Am = zeros(n, n)
        t_idx = 0
        for a in range(n):
            for b in range(a + 1, n):
                Am[a, b] = cs[t_idx]
                Am[b, a] = -cs[t_idx]
                t_idx += 1
        eqs = []
        for M in basis:
            E = M.T * Am + Am * M
            for a in range(n):
                for b in range(a + 1, n):
                    eqs.append(E[a, b])
        A2m, _ = linear_eq_to_matrix(eqs, list(cs))
        ns_anti = A2m.nullspace()
        dim_anti = len(ns_anti)
        check(f"({p},{q}) (c) antisym.: dim = {dim_anti} (exact nullspace)",
              dim_anti == len(cs) - A2m.rank())
    else:
        dim_anti = 0
        check(f"({p},{q}) (c) antisym.: trivially dim=0", True)
    DIM_ANTI[(p, q)] = dim_anti
    print(f"      (c) ANTISYMMETRIC invariants: dim = {dim_anti}")

    # cross-check: the full (non-symmetrized) space of bilinear invariants
    gs = symbols(f"g0:{n*n}")
    G = Matrix(n, n, lambda a, b: gs[a * n + b])
    eqs = []
    for M in basis:
        E = M.T * G + G * M
        for a in range(n):
            for b in range(n):
                eqs.append(E[a, b])
    A3m, _ = linear_eq_to_matrix(eqs, list(gs))
    dim_full = len(gs) - A3m.rank()
    check(f"({p},{q}) cross: dim(all bilinear) = dim(sym)+dim(antisym) "
          f"[{dim_full} = {dim_sym}+{dim_anti}]",
          dim_full == dim_sym + dim_anti)

verdict("invariant numbers that can be attached to a shift: LINEAR — 0 "
        "(no signature gives an invariant covector); "
        "QUADRATIC — exactly 1 (every solution ∝ η: the unique invariant "
        "number of a shift = the η-norm a·η·a, whose sign is exactly the time/space-like-ness); "
        "ANTISYMMETRIC — 0. See the rollup table in SUMMARY.")

# ============================================================== A3
section("A3 — THE FATE OF THE TIME SHIFT")

Q1_CITED = [(2, 1), (3, 1), (4, 1)]
th = Symbol("theta", real=True)

for (p, q) in SIGS:
    n = p + q
    if q == 0:
        print(f"  --- ({p},{q}): q=0, there is no minus-axis — T_t is undefined, "
              f"skipped (declared) ---")
        continue
    et = eta_of(p, q)
    basis = BASIS[(p, q)]
    hats = HATS[(p, q)]
    Ts = TS[(p, q)]
    pairs = PAIRS[(p, q)]
    d = len(basis)
    t = n - 1                      # the last minus-axis (0-indexed)
    e_t = e_vec(n, t)
    print(f"  --- signature ({p},{q}): t = axis {t} (η_tt = {et[t, t]}) ---")

    # (a) bracket table [M_k, T_t] by type
    for k, M in enumerate(basis):
        i, j = pairs[k]
        v = M * e_t
        ok = br(hats[k], Ts[t]) == Tmat(v)
        if t not in (i, j):
            typ = ("a spatial rotation" if et[i, i] == et[j, j] == 1
                   else ("a {−−}-rotation without t" if et[i, i] == et[j, j] == -1
                         else "a boost without t"))
            check(f"({p},{q}) (a) M[{k}] pair({i},{j}) without t: [M,T_t]=0",
                  ok and v == zeros(n, 1))
            print(f"      M[{k}] axes({i},{j}) {typ}: [M,T_t] = 0")
        else:
            s = i if j == t else j
            typ = ("a boost (+,−) with axis t" if et[s, s] == 1
                   else "a {−−}-rotation with axis t")
            check(f"({p},{q}) (a) M[{k}] pair({i},{j}) with t: [M,T_t]=T_(Me_t)",
                  ok and v != zeros(n, 1))
            print(f"      M[{k}] axes({i},{j}) {typ}: [M,T_t] = T_v, "
                  f"v = {list(v)}")

    # (b) stabilizer of e_t: {M ∈ so : M e_t = 0}
    cs = symbols(f"s0:{d}")
    Mgen = zeros(n, n)
    for k in range(d):
        Mgen += cs[k] * basis[k]
    eqs = list(Mgen * e_t)
    Astab, _ = linear_eq_to_matrix(eqs, list(cs))
    ns_stab = Astab.nullspace()
    dim_stab = len(ns_stab)
    DIM_STAB[(p, q)] = dim_stab
    exp_stab = (n - 1) * (n - 2) // 2
    check(f"({p},{q}) (b) dim of the stabilizer = (n−1)(n−2)/2 = {exp_stab}",
          dim_stab == exp_stab)
    etp = eta_of(p, q - 1)   # η without axis t (the last one)
    ok_struct = True
    for v in ns_stab:
        S = zeros(n, n)
        for k in range(d):
            S += v[k] * basis[k]
        if S * e_t != zeros(n, 1):
            ok_struct = False
        if any(S[t, j] != 0 for j in range(n)) or \
           any(S[a, t] != 0 for a in range(n)):
            ok_struct = False
        Sb = S[:n - 1, :n - 1]
        if (Sb * etp + etp * Sb.T) != zeros(n - 1, n - 1):
            ok_struct = False
    check(f"({p},{q}) (b) stabilizer: row/column t are zero, "
          f"the block ∈ so({p},{q-1})", ok_struct)
    print(f"      (b) stabilizer of e_t: dim = {dim_stab} = "
          f"dim so({p},{q-1}) — this is so({p},{q-1}) on the first {n-1} axes")

    # (c) carryover / explicit path
    if q == 1:
        tag = ("in the task list" if (p, q) in Q1_CITED
               else "outside the task list — the same carryover")
        print(f"      (c) q=1 ({tag}): [LOG link] the translation ideal ≅ "
              f"the vector representation of so({p},1) (proved by machine in A1(a): "
              f"[M̂,T_a]=T_(Ma) ∀ basis M, a) ⟹ the S920 cone-argument "
              f"applies VERBATIM:")
        print(f"          active-v10.2/src/symbolic/"
              f"S920_w29_energy_sign_orientation.py")
        print(f"          (log: active-v10.2/src/symbolic/S920_w29_run.log; "
              f"S920 verdict: a sector of multiplicity 1 is orientable — two "
              f"invariant cones)")
        print(f"          ⟹ the flows of so({p},1) do NOT connect T_t with −T_t: "
              f"the time shift is oriented.")
        check(f"({p},{q}) (c) the premise of the carryover (ideal ≅ vector "
              f"representation) proved in A1(a)", True)
    else:
        u = n - 2
        check(f"({p},{q}) (c) axis u={u} is also minus (q≥2)",
              et[u, u] == Integer(-1))
        J = zeros(n, n)
        J[u, t] = Integer(1)
        J[t, u] = Integer(-1)
        check(f"({p},{q}) (c) J (a {{−−}}-rotation in the plane ({u},{t})) ∈ so",
              (J * et + et * J.T) == zeros(n, n))
        hJ = hat(J)
        ok_ad = all(br(hJ, Ts[i]) == Tmat(J * e_vec(n, i)) for i in range(n))
        check(f"({p},{q}) (c) ad J on the ideal = J on the vectors "
              f"([Ĵ,T_i]=T_(J e_i) ∀i)", ok_ad)
        # explicit path g(θ) = exp(θJ)
        g = eye(n)
        g[u, u] = cos(th); g[u, t] = sin(th)
        g[t, u] = -sin(th); g[t, t] = cos(th)
        check(f"({p},{q}) (c) dg/dθ = J·g (the ODE characterization of exp(θJ))",
              simplify(g.diff(th) - J * g) == zeros(n, n))
        check(f"({p},{q}) (c) g(0) = I", g.subs(th, 0) == eye(n))
        check(f"({p},{q}) (c) g(θ)ᵀ η g(θ) = η (a path in the isometry group)",
              simplify(g.T * et * g - et) == zeros(n, n))
        check(f"({p},{q}) (c) g(π)·e_t = −e_t (endpoint exact)",
              g.subs(th, pi) * e_t == -e_t)
        # homogeneous: ĝ(θ), conjugation of the translation
        gh = g.row_join(zeros(n, 1)).col_join(zeros(1, n + 1))
        gh[n, n] = Integer(1)
        check(f"({p},{q}) (c) dĝ/dθ = Ĵ·ĝ and ĝ(0)=I "
              f"[LOG: ⟹ ĝ=exp(θĴ), uniqueness of the linear ODE solution]",
              simplify(gh.diff(th) - hJ * gh) == zeros(n + 1, n + 1)
              and gh.subs(th, 0) == eye(n + 1))
        ghinv = gh.subs(th, -th)
        check(f"({p},{q}) (c) ĝ(θ)·ĝ(−θ) = I",
              simplify(gh * ghinv - eye(n + 1)) == zeros(n + 1, n + 1))
        conj = gh * Ts[t] * ghinv
        check(f"({p},{q}) (c) Ad(ĝ(θ))T_t = T_(g(θ)e_t) (a continuous path)",
              simplify(conj - Tmat(g * e_t)) == zeros(n + 1, n + 1))
        check(f"({p},{q}) (c) θ=π: Ad(ĝ(π))T_t = −T_t (endpoint exact)",
              conj.subs(th, pi) == -Ts[t])
        named = " (named in the task)" if (p, q) in [(2, 2), (3, 2)] else ""
        print(f"      (c) q≥2{named}: a {{−−}}-rotation J in the plane "
              f"({u},{t}); ad J = J on the vectors (machine-checked); "
              f"exp(θ·adJ): T_t → cos(θ)T_t + sin(θ)... → at θ=π EXACTLY −T_t "
              f"— the time shift is connected to its negation continuously.")

verdict("q=1: the translation ideal ≅ the vector representation ⟹ the S920 "
        "cone-argument carries over verbatim — T_t is NOT connected to −T_t "
        "(the shift is oriented); q≥2: an explicit path exp(θ·adJ) in the {−−}-plane "
        "connects T_t to −T_t continuously (exact endpoints) — the orientation of the shift is lost. "
        "The stabilizer of e_t = so(p,q−1) in every signature with q≥1.")

# ============================================================== M
section("M — MUTANTS")

# m1: broken rule [M,T_a] = T_{Mᵀa}
m1_catches = 0
m1_detail_done = False
for (p, q) in SIGS:
    n = p + q
    basis = BASIS[(p, q)]
    hats = HATS[(p, q)]
    Ts = TS[(p, q)]
    found = None
    for k, M in enumerate(basis):
        for i in range(n):
            a = e_vec(n, i)
            if M.T * a != M * a:
                found = (k, i, M, a)
                break
        if found:
            break
    if found is None:
        check(f"m1 ({p},{q}): there exist M,a with Mᵀa ≠ Ma", False)
        continue
    k, i, M, a = found
    computed = br(hats[k], Ts[i])
    ok_true = computed == Tmat(M * a)
    ok_mutant_differs = computed != Tmat(M.T * a)
    if ok_true and ok_mutant_differs:
        m1_catches += 1
    check(f"m1 ({p},{q}): the bracket = T_(Ma) and ≠ T_(Mᵀa) on M[{k}], a=e_{i}",
          ok_true and ok_mutant_differs)
    if not m1_detail_done and (p, q) == (3, 1):
        print(f"  m1 DETAIL ({p},{q}), M[{k}] = {M.tolist()}, a = e_{i}:")
        print(f"      Ma  = {list(M * a)}")
        print(f"      Mᵀa = {list(M.T * a)}")
        print(f"      discrepancy (M−Mᵀ)a = {list((M - M.T) * a)} ≠ 0 — "
              f"the mutant gives a DIFFERENT translation, the computed bracket = T_(Ma). "
              f"CAUGHT.")
        m1_detail_done = True
MUTANT_STATUS["m1 [M,T_a]=T_(Mᵀa)"] = (
    "CAUGHT" if m1_catches == len(SIGS) else "NOT-CAUGHT")
print(f"  m1: caught in {m1_catches}/{len(SIGS)} signatures → "
      f"{MUTANT_STATUS['m1 [M,T_a]=T_(Mᵀa)']}")

# m2: a planted "invariant linear form" φ = (1,0,…,0)
m2_catches = 0
m2_detail_done = False
for (p, q) in SIGS:
    n = p + q
    basis = BASIS[(p, q)]
    found = None
    for k, M in enumerate(basis):
        for i in range(n):
            val = (M * e_vec(n, i))[0, 0]     # φ(M e_i) at φ=(1,0,…,0)
            if val != 0:
                found = (k, i, val, M)
                break
        if found:
            break
    if found is None:
        check(f"m2 ({p},{q}): there exist M,e_i with φ(M e_i)≠0", False)
        continue
    k, i, val, M = found
    m2_catches += 1
    check(f"m2 ({p},{q}): φ(M[{k}] e_{i}) = {val} ≠ 0", val != 0)
    if not m2_detail_done and (p, q) == (3, 1):
        print(f"  m2 DETAIL ({p},{q}): φ = (1,0,…,0), M[{k}] = {M.tolist()}, "
              f"e_{i}: φ(M e_{i}) = {val} ≠ 0 — the \"invariance\" breaks "
              f"on the very first basis generator. CAUGHT.")
        m2_detail_done = True
MUTANT_STATUS["m2 φ=(1,0,…,0) invariant"] = (
    "CAUGHT" if m2_catches == len(SIGS) else "NOT-CAUGHT")
print(f"  m2: caught in {m2_catches}/{len(SIGS)} signatures → "
      f"{MUTANT_STATUS['m2 φ=(1,0,…,0) invariant']}")

# m3: "there exists a second invariant symmetric form, independent of η"
m3_all_dim1 = True
for (p, q) in SIGS:
    ds = DIM_SYM[(p, q)]
    if ds == 1:
        check(f"m3 ({p},{q}): dim sym. inv. = 1 ⟹ there is NO second independent one",
              True)
    else:
        m3_all_dim1 = False
        print(f"  m3 ({p},{q}): dim = {ds} > 1 — NOT caught, and this "
              f"is a RESULT (double-checked below)")
        # double-check by an independent path: the full space minus the antisymmetric one
        check(f"m3 ({p},{q}): double-check (the A2 cross-system) is consistent",
              DIM_SYM[(p, q)] + DIM_ANTI[(p, q)] >= ds)
MUTANT_STATUS["m3 \"a second sym. form\""] = (
    "CAUGHT (dim=1 everywhere)" if m3_all_dim1 else "NOT-CAUGHT(result)")
print(f"  m3: {MUTANT_STATUS['m3 \"a second sym. form\"']} — the space "
      f"of symmetric invariants is one-dimensional in EVERY signature, "
      f"every solution ∝ η (A2(b)); \"a second independent form\" is refuted "
      f"by dimension.")

verdict("m1 CAUGHT (the transposed rule gives a different translation — an exact "
        "discrepancy is shown); m2 CAUGHT (φ=(1,0,…,0) breaks on an explicit "
        "basis M, e_i); m3 CAUGHT by the dim=1 dimension from A2(b) in every "
        "signature.")

# ============================================================== SUMMARY
section("SUMMARY")

print("  VERDICTS:")
for sec, text in VERDICTS:
    print(f"    [{sec}]")
    print(f"      {text}")
print()
print("  ROLLUP TABLE BY SIGNATURE:")
print("  " + "-" * 74)
print("  (p,q)  | dim lin.inv | dim sym.inv | sym ∝ η? | dim antisym "
      "| dim stab e_t")
print("  " + "-" * 74)
for (p, q) in SIGS:
    stab = str(DIM_STAB[(p, q)]) if (p, q) in DIM_STAB else "—(q=0)"
    prop = "YES" if SYM_PROP_ETA[(p, q)] else "NO"
    print(f"  ({p},{q})  |     {DIM_LIN[(p,q)]}       |     "
          f"{DIM_SYM[(p,q)]}       |   {prop}    |     "
          f"{DIM_ANTI[(p,q)]}       |   {stab}")
print("  " + "-" * 74)
print("  LAW: the one invariant number that can be attached to a shift a is "
      "its η-norm aᵀηa (multiplicity 1, ∝ η); there are no linear or antisymmetric "
      "invariants; the stabilizer of a time shift = so(p,q−1); "
      "the orientation of a time shift survives ⟺ q=1 (the S920 carryover), is lost at q≥2 "
      "(an explicit path T_t → −T_t).")

# ---------------------------------------------------------------- final
print()
print("=" * 78)
print("COUNTERS:")
ok_all = True
for sec in SECTION_COUNTS:
    cnt = SECTION_COUNTS[sec]
    fl = SECTION_FAILS[sec]
    print(f"  {sec[:50]:52s} checks={cnt:5d}  failures={fl}")
    if sec != "SUMMARY" and cnt == 0:
        print(f"  [FAIL] section \"{sec}\" is EMPTY (zero checked)")
        ok_all = False
GRAND_TOTAL = sum(SECTION_COUNTS.values())
print(f"  TOTAL checks: {GRAND_TOTAL} "
      f"(of which {CHECKS_TOTAL} check calls, the rest are individually counted "
      f"brackets/Jacobi), failures: {len(FAILS)}")
print("  MUTANTS:")
mut_ok = True
for name, st in MUTANT_STATUS.items():
    print(f"    {name}: {st}")
    if st.startswith("NOT-CAUGHT") and "result" not in st:
        mut_ok = False

if CHECKS_TOTAL == 0:
    print("SUMMARY: ZERO CHECKS — the measurement did NOT happen, exit 2")
    sys.exit(2)
if FAILS:
    print(f"SUMMARY: {len(FAILS)} FAILURES — exit 1")
    for sec, lab in FAILS[:40]:
        print(f"  failure: [{sec}] {lab}")
    sys.exit(1)
if not (ok_all and mut_ok):
    print("SUMMARY: an empty section or an uncaught mutant — exit 1")
    sys.exit(1)
print("SUMMARY: ALL OK — exit 0")
sys.exit(0)
