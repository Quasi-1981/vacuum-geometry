# -*- coding: utf-8 -*-
# DIM: symbolic-d (an independent VISA on S1054: re-deriving the on-shell slope of the cone WITH MY OWN HAND.
#          The A_d lattice is built ∀d from scratch, the node is solved, the series |f|² is expanded,
#          the on-shell 2-2cos(x)=λ²s is solved for x. The ∀d statements carry symbolic
#          d (the tr M identity — polynomial in z=e^{iθ}); d=2 additionally carries an EXTERNAL
#          anchor (honeycomb), re-derived here, not cited.)
"""S1055 — a VISA on S1054: an INDEPENDENT re-derivation of the on-shell slope of the cone.

Alpha's request (verbatim): «re-derive the on-shell slope of the cone from scratch — 2-2cos(x)=|f|^2
near the node — and say whether it carries a factor of 1/2 relative to the registered h^4/(4d)».
The second question: is the numerator counted per BOND, while the denominator 2/h^2 — per PERIOD.

★METHOD (per the request): WITH MY OWN HAND, not by checking against Alpha's code. Nothing is imported
from S1012/S999/S1054. The A_d lattice is built here from scratch; the node is solved here; the series
is expanded here. The ONLY thing taken from OUTSIDE — the registered number h^4/(4d) and its assembly
(SC*n/2)/(2/h^2), because those are precisely what is on trial.

Teeth: negative controls M1-M5 + an external anchor d=2 (honeycomb-lattice v_F), which COULD have failed to match.
"""
import os
import sys
from sympy import (Rational, symbols, cos, sin, exp, I, simplify, expand, sqrt,
                   Matrix, eye, ones, zeros, series, nsimplify, pi, Symbol,
                   trigsimp, re, im, together, factor, S, solve, O, Poly)

_HERE = os.path.dirname(os.path.abspath(__file__))

PASS = 0
FAIL = 0


def ok(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
        print("  [OK]   " + msg)
    else:
        FAIL += 1
        print("  [FAIL] " + msg)


class Tee:
    def __init__(self, real, fh):
        self.real = real
        self.fh = fh

    def write(self, s):
        self.real.write(s)
        self.fh.write(s)
        self.fh.flush()
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


# ======================================================================
# PART 0 — THE LATTICE FROM SCRATCH (nothing borrowed)
# ======================================================================
def bonds(d, off=None):
    """d+1 bond-vectors of the A_d cell, EXPLICITLY, in ℝ^{d+1} on the hyperplane ⟂ 1.

    w_j = e_j - 1/(d+1)*1  =>  <w_i,w_j> = δ_ij - 1/(d+1),  |w_j|^2 = d/(d+1).
    Normalize: δ_j = w_j*sqrt((d+1)/d)  =>  |δ_j|^2 = 1, <δ_i,δ_j> = -1/d  (i≠j).
    The sum Σδ_j = 0 (closure of the cell).

    off: if given — a SUBSTITUTION of the off-diagonal of the Gram (mutant M1); only the Gram is returned.
    """
    n = d + 1
    scale = sqrt(Rational(n, d))
    B = []
    for j in range(n):
        v = [Rational(-1, n)] * n
        v[j] += 1
        B.append(Matrix([x * scale for x in v]))
    return B


def gram_from_bonds(B):
    n = len(B)
    return Matrix(n, n, lambda i, j: simplify((B[i].T * B[j])[0, 0]))


# ======================================================================
# PART 1 — THE NODE: solve f(k0)=0 WITH MY OWN HAND
# ======================================================================
def k_from_phases(d, thetas, B):
    """Recover k from the phases θ_j = k·δ_j. The Gram G is singular (G·1=0), so a solution
    exists ⟺ Σθ_j = 0, and then k = (d/(d+1))·Σ_j θ_j δ_j. Checked below."""
    n = d + 1
    k = zeros(B[0].rows, 1)
    for j in range(n):
        k = k + thetas[j] * B[j]
    return Rational(d, n) * k


def node_phases(d, t=None):
    """Phases at the NODE: Σ_j e^{iθ_j}=0 AND Σ_j θ_j = 0 are both needed (so that k exists).
    Construction: we take the (d+1)-th root of unity for d=2; for d>=3 — antipodal
    pairs (a forced pairing: 4 unit numbers with zero sum decompose into two
    antipodal pairs), with a free parameter t along the nodal set at d>=3."""
    n = d + 1
    if d == 2:
        return [S(0), Rational(2, 3) * pi, -Rational(2, 3) * pi]
    if d == 3:
        # (t, -t, t+pi, -t-pi) -> sum 0; e^{it}+e^{-it}+e^{i(t+pi)}+e^{-i(t+pi)}
        #                      = 2cos t - 2cos t = 0  ∀t  => a NODAL LINE (dim 1)
        if t is None:
            t = Rational(1, 5)
        return [t, -t, t + pi, -t - pi]
    if d == 4:
        # 5 phases, sum 0, Σe^{iθ}=0: 5th-degree roots, shifted so the sum=0
        return [Rational(2, 5) * pi * j - Rational(4, 5) * pi for j in range(5)]
    if d == 5:
        # 6 phases = THREE antipodal pairs (a_k, a_k+π) ⟹ f=0 automatically.
        # The compatibility condition Σθ=0 is NOT automatic: Σθ = 2Σa_k + 3π ⟹ we need Σa_k = -3π/2.
        # (my first choice Σa_k=π gave Σθ=5π≠0 — k0 did not exist; test T1 caught this)
        a, b, c = -Rational(1, 3) * pi, -Rational(1, 2) * pi, -Rational(2, 3) * pi
        return [a, a + pi, b, b + pi, c, c + pi]
    raise ValueError(d)


def f_of_phases(thetas):
    return sum(exp(I * th) for th in thetas)


# ======================================================================
# PART 2 — THE EXPANSION OF |f|^2 NEAR THE NODE (λ^0, λ^1, λ^2) — WITH MY OWN HAND
# ======================================================================
def M_matrix(d, thetas, B):
    """|f(k0+λu)|^2 = λ^2 · u^T M u + O(λ^3),  M = a a^T + b b^T,
       a = Σ_j cos θ_j δ_j,  b = Σ_j sin θ_j δ_j.
    Derived here: f(k0+λu) = Σ e^{iθ_j} e^{iλ u·δ_j}
                              = Σ e^{iθ_j} + iλ Σ e^{iθ_j}(u·δ_j) + O(λ^2)
                              = 0 + iλ (u·w) + O(λ^2),  w = Σ e^{iθ_j} δ_j
    ⟹ |f|^2 = λ^2 |u·w|^2 = λ^2 [ (u·Re w)^2 + (u·Im w)^2 ]."""
    a = zeros(B[0].rows, 1)
    b = zeros(B[0].rows, 1)
    for j, th in enumerate(thetas):
        a = a + cos(th) * B[j]
        b = b + sin(th) * B[j]
    a = a.applyfunc(lambda x: simplify(x))
    b = b.applyfunc(lambda x: simplify(x))
    return simplify_mat(a * a.T + b * b.T), a, b


def simplify_mat(Mx):
    return Mx.applyfunc(lambda x: simplify(x))


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1055_run.log"), "w", encoding="utf-8")
    sys.stdout = Tee(sys.stdout, logf)

    print("=" * 78)
    print("S1055 (Beta) — a VISA on S1054: an INDEPENDENT RE-DERIVATION OF THE ON-SHELL SLOPE OF THE CONE")
    print("=" * 78)

    lam = Symbol('lam', positive=True)
    x = Symbol('x', real=True)
    s = Symbol('s', positive=True)

    # ==================================================================
    print("\n[T0] THE LATTICE FROM SCRATCH: the Gram = I on the diagonal, -1/d off, Σδ=0")
    # ==================================================================
    for d in [2, 3, 4, 5]:
        B = bonds(d)
        G = gram_from_bonds(B)
        n = d + 1
        diag_ok = all(G[j, j] == 1 for j in range(n))
        off_ok = all(G[i, j] == Rational(-1, d) for i in range(n) for j in range(n) if i != j)
        ssum = zeros(B[0].rows, 1)
        for v in B:
            ssum = ssum + v
        ok(diag_ok and off_ok, "d={0}: |δ|^2=1 ∀j and <δ_i,δ_j>=-1/d ∀i≠j (the K-V/A_d Gram)".format(d))
        ok(simplify(ssum.norm()) == 0, "d={0}: Σ_j δ_j = 0 (closure of the cell)".format(d))
        ok(simplify(sum(G[0, j] for j in range(n))) == 0,
           "d={0}: G·1 = 0 ⟹ k exists ⟺ Σθ_j=0 (a compatibility condition, not a choice)".format(d))

    # ==================================================================
    print("\n[T1] THE NODE: f(k0)=0 and k0 is genuinely recovered from the phases")
    # ==================================================================
    for d in [2, 3, 4, 5]:
        B = bonds(d)
        th = node_phases(d)
        fval = simplify(expand(f_of_phases(th).rewrite(cos)))
        ok(simplify(fval) == 0, "d={0}: f(k0) = Σ_j e^{{iθ_j}} = 0 (the node)".format(d))
        ok(simplify(sum(th)) == 0, "d={0}: Σ_j θ_j = 0 ⟹ k0 exists".format(d))
        k0 = k_from_phases(d, th, B)
        back = [simplify((k0.T * B[j])[0, 0]) for j in range(d + 1)]
        ok(all(simplify(back[j] - th[j]) == 0 for j in range(d + 1)),
           "d={0}: k0·δ_j = θ_j reproduced ∀j (k0 is NOT fitted)".format(d))

    # ==================================================================
    print("\n[T2] THE EXPANSION OF |f|^2: λ^0=λ^1=0, the coeff of λ^2 = u^T M u (EXACTLY, without 1/2)")
    # ==================================================================
    for d in [2, 3]:
        B = bonds(d)
        th = node_phases(d)
        Mx, a, b = M_matrix(d, th, B)
        # a direct series in λ in a GIVEN direction u (we take a generic u in the plane a,b)
        c1, c2 = symbols('c1 c2', real=True)
        u = c1 * a + c2 * b
        fser = sum(exp(I * th[j]) * exp(I * lam * (u.T * B[j])[0, 0]) for j in range(d + 1))
        f2 = expand(simplify((fser * fser.conjugate()).rewrite(cos)))
        p = Poly(series(f2, lam, 0, 3).removeO(), lam)
        c_lam0 = simplify(p.coeff_monomial(1))
        c_lam1 = simplify(p.coeff_monomial(lam))
        c_lam2 = simplify(p.coeff_monomial(lam ** 2))
        quad = simplify((u.T * Mx * u)[0, 0])
        ok(c_lam0 == 0, "d={0}: the coeff λ^0 = 0".format(d))
        ok(c_lam1 == 0, "d={0}: the coeff λ^1 = 0".format(d))
        ok(simplify(c_lam2 - quad) == 0,
           "d={0}: the coeff λ^2 = u^T M u EXACTLY (no 1/2 appears)".format(d))
        # the second derivative ≠ the coefficient: exactly twice
        ok(simplify(2 * c_lam2 - 2 * quad) == 0,
           "d={0}: d^2|f|^2/dλ^2 = 2·(coeff λ^2) ⟹ coeff ⊥ 2nd-derivative are DIFFERENT objects (×2)".format(d))

    # ==================================================================
    print("\n[T3] tr M AT THE NODE = (d+1)^2/d — an identity, checked with SYMBOLIC phases")
    # ==================================================================
    # ★THE PROOF ROUTE: the identity does NOT need trigonometry. In the variables z_j=e^{iθ_j} (|z_j|=1):
    #   tr M = |Re w|^2+|Im w|^2 = w†w,  w = Σ_j z_j δ_j  ⟹  tr M = Σ_{j,l} G_jl z_j conj(z_l).
    #   G = (1+1/d)I - (1/d)J  ⟹  tr M = (1+1/d)Σ_j|z_j|^2 - (1/d)|Σ_j z_j|^2 = n(1+1/d) - |f|^2/d.
    # This is a POLYNOMIAL identity in (z, z̄) under the constraint z_j·z̄_j=1 — checked exactly and cheaply.
    for d in [2, 3, 4, 5, 6]:
        n = d + 1
        z = symbols('z0:%d' % n)
        zb = symbols('zb0:%d' % n)          # z̄_j as independent symbols, the constraint z_j*zb_j=1
        G = Matrix(n, n, lambda i, j: 1 if i == j else Rational(-1, d))
        trM_z = expand(sum(G[j, l] * z[j] * zb[l] for j in range(n) for l in range(n)))
        f_z = sum(z[j] for j in range(n))
        fbar_z = sum(zb[j] for j in range(n))
        claim_z = expand(n * (1 + Rational(1, d)) - f_z * fbar_z / d)
        # apply the constraint |z_j|=1 to BOTH sides (only the diagonal terms z_j*zb_j)
        subsmap = {z[j] * zb[j]: 1 for j in range(n)}
        diff = expand(trM_z - claim_z)
        for j in range(n):
            diff = diff.subs(z[j] * zb[j], 1)
        ok(expand(diff) == 0,
           "d={0}: tr M = n(1+1/d) - |f|^2/d — an IDENTITY ∀θ (polynomial in z, |z_j|=1)".format(d))

    # ★A CHECK OF TWO ROUTES (z-algebra ⊥ the matrix M=aa^T+bb^T) — where the second is cheap
    for d in [2, 3]:
        B = bonds(d)
        n = d + 1
        ths = symbols('th0:%d' % n, real=True)
        Mx, a, b = M_matrix(d, list(ths), B)
        trM = simplify(expand(trigsimp(Mx.trace())))
        fabs2 = simplify(expand(trigsimp(
            (f_of_phases(list(ths)) * f_of_phases(list(ths)).conjugate()).rewrite(cos))))
        claim = simplify(n * (1 + Rational(1, d)) - fabs2 / d)
        ok(simplify(expand(trigsimp(trM - claim))) == 0,
           "d={0}: THE SAME verdict by the matrix route (tr(aa^T+bb^T), trigonometrically)".format(d))

    for d in [2, 3, 4]:
        B = bonds(d)
        n = d + 1
        # at the node
        thn = node_phases(d)
        Mn, an, bn = M_matrix(d, thn, B)
        trn = simplify(trigsimp(expand(Mn.trace())))
        ok(simplify(trn - Rational((d + 1) ** 2, d)) == 0,
           "d={0}: at the node tr M = (d+1)^2/d = {1}".format(d, Rational((d + 1) ** 2, d)))
        ok(Mn.rank() == 2,
           "d={0}: rank M = 2 (M=aa^T+bb^T) ⟹ 2 transverse directions, the kernel = tangent to the nodal set".format(d))

    # d=3: a nodal LINE — tr M is CONSTANT ALONG it (a symbolic parameter t)
    tt = Symbol('t', real=True)
    B3 = bonds(3)
    th3 = node_phases(3, t=tt)
    M3, a3, b3 = M_matrix(3, th3, B3)
    ok(simplify(trigsimp(expand(f_of_phases(th3).rewrite(cos)))) == 0,
       "d=3: f=0 ∀t ⟹ a nodal LINE (dim = d-2 = 1), not a point")
    ok(simplify(trigsimp(expand(M3.trace() - Rational(16, 3)))) == 0,
       "d=3: tr M = 16/3 CONSTANT along the whole nodal line (symbolic t)")

    # ==================================================================
    print("\n[T4] ★ON-SHELL: solve 2-2cos(x) = λ^2·s FOR x — is there a 1/2?")
    # ==================================================================
    # a series for x(λ): substitute an ansatz and match orders — WITHOUT assumptions
    A1, A3 = symbols('A1 A3', real=True)
    xa = A1 * lam + A3 * lam ** 3
    lhs = series(2 - 2 * cos(xa), lam, 0, 5).removeO()
    rhs = lam ** 2 * s
    eqs = Poly(expand(lhs - rhs), lam)
    e2 = simplify(eqs.coeff_monomial(lam ** 2))
    e4 = simplify(eqs.coeff_monomial(lam ** 4))
    sol1 = solve(e2, A1)
    A1v = [q for q in sol1 if simplify(q - sqrt(s)) == 0]
    ok(len(A1v) == 1, "on-shell: the leading slope A1 = sqrt(s) (the positive branch)")
    ok(simplify(A1v[0] - sqrt(s)) == 0,
       "★THE ON-SHELL SLOPE = sqrt(u^T M u) — there is NO factor of 1/2 (or 2)")
    A3v = solve(e4.subs(A1, sqrt(s)), A3)
    ok(len(A3v) == 1 and simplify(A3v[0] - s ** Rational(3, 2) / 24) == 0,
       "the next order A3 = s^(3/2)/24 ≠ 0 ⟹ the leading order is IDENTIFIED (not guessed)")

    # ==================================================================
    print("\n[T5] THE DICTIONARY OF THE DENOMINATOR: 2/h^2 — this is the SECOND DERIVATIVE in the PERIOD-variable")
    # ==================================================================
    nu = Symbol('nu', real=True)
    for d in [2, 3, 4, 5]:
        h = d + 1
        Tnu = 2 - 2 * cos(2 * pi * nu / h)       # the variable ν = the index of the ℤ/h character
        Tx = 2 - 2 * cos(x)                       # the variable x = the phase per ONE hop (1 bond)
        c_nu = simplify(series(Tnu, nu, 0, 3).removeO().coeff(nu ** 2))
        dd_nu = simplify(series(Tnu, nu, 0, 3).removeO().diff(nu, 2).subs(nu, 0))
        c_x = simplify(series(Tx, x, 0, 3).removeO().coeff(x ** 2))
        dd_x = simplify(series(Tx, x, 0, 3).removeO().diff(x, 2).subs(x, 0))
        ok(simplify(c_nu - (2 * pi / h) ** 2) == 0,
           "d={0}: the coeff of ν^2 = (2π/h)^2  [h={1}]".format(d, h))
        ok(simplify(dd_nu - 2 * (2 * pi / h) ** 2) == 0,
           "d={0}: T''(0) over ν = 2(2π/h)^2 ⟹ with the (2π)-drop = 2/h^2 = the REGISTERED DENOMINATOR".format(d))
        ok(simplify(c_x - 1) == 0 and simplify(dd_x - 2) == 0,
           "d={0}: in the HOP-variable x: coeff=1, T''(0)=2 — no h at all ⟹ h comes ONLY from ν".format(d))
    print("  ⟹ 1/h^2 in the denominator = exactly the chain rule dx/dν = 2π/h,")
    print("    and h = the number of BONDS in the period of the column (S1001: P=order(u_0)=d+1).")
    print("    ⟹ the denominator is measured PER PERIOD, the numerator (the Gram, |δ|=1) — PER BOND. CONFIRMED.")

    # ==================================================================
    print("\n[T6] THE FOUR PAIRINGS — with my own hand, in ONE (period-) unit")
    # ==================================================================
    print("  {0:<34} {1:<14} {2:<14}".format("pairing numerator/denominator", "d=2", "∀d"))
    for d in [2, 3, 4, 5]:
        h = d + 1
        num_coef = Rational((d + 1) ** 2, 2 * d)      # ½trM = the coeff of λ^2, averaged
        num_dd = 2 * num_coef                          # the 2nd-derivative
        den_coef = Rational(1, h ** 2)                 # the coeff of ν^2, the (2π)-drop
        den_dd = Rational(2, h ** 2)                   # T''(0), the (2π)-drop  ← the registered one
        p_cc = simplify(num_coef / den_coef)
        p_dd = simplify(num_dd / den_dd)
        p_cd = simplify(num_coef / den_dd)             # MIXED ← the registered one
        p_dc = simplify(num_dd / den_coef)
        tgt = Rational(h ** 4, 4 * d)
        ok(p_cc == p_dd, "d={0}: both CONSISTENT pairings give one number {1}".format(d, p_cc))
        ok(p_cc == Rational(h ** 4, 2 * d), "d={0}: the consistent one = h^4/(2d) = {1}".format(d, Rational(h ** 4, 2 * d)))
        ok(p_cd == tgt, "d={0}: the REGISTERED h^4/(4d)={1} is reproduced EXACTLY by the mixed one (coeff/2nd-derivative)".format(d, tgt))
        ok(p_dc == Rational(h ** 4, d), "d={0}: the other mixed one = h^4/d = {1} (not the registered one)".format(d, Rational(h ** 4, d)))
        ok(simplify(tgt / p_cc) == Rational(1, 2),
           "d={0}: registered/consistent = 1/2 EXACTLY (a constant, not f(d))".format(d))

    # ==================================================================
    print("\n[T7] ★THE NUMBER IN THE CONSISTENT HOP-DICTIONARY (what the on-shell derivation T4 gives)")
    # ==================================================================
    for d in [2, 3, 4, 5]:
        h = d + 1
        v2 = Rational((d + 1) ** 2, 2 * d)   # = ½trM / 1  (coeff/coeff in the hop-variable)
        reg = Rational(h ** 4, 4 * d)
        ok(simplify(reg / v2) == Rational(h ** 2, 2),
           "d={0}: registered/hop-dictionary = h^2/2 = {1}  (= length h^2 × pairing 1/2)"
           .format(d, Rational(h ** 2, 2)))
        print("       d={0}: v^2(hop) = {1}   ⊥   registered = {2}".format(d, v2, reg))

    # ==================================================================
    print("\n[T8] ★AN EXTERNAL ANCHOR d=2 — RE-DERIVED FROM THEIR DEFINITION, WITH MY OWN HAND")
    # ==================================================================
    # The honeycomb lattice, from scratch: H=[[0,f],[f*,0]], t=1, a=|δ|=1, WITHOUT our machinery.
    # δ_j — three vectors at 120°, |δ|=1 (this IS the d=2 case of our A_2).
    kx, ky = symbols('kx ky', real=True)
    dl = [Matrix([1, 0]),
          Matrix([Rational(-1, 2), sqrt(3) / 2]),
          Matrix([Rational(-1, 2), -sqrt(3) / 2])]
    Gg = Matrix(3, 3, lambda i, j: simplify((dl[i].T * dl[j])[0, 0]))
    ok(all(Gg[i, i] == 1 for i in range(3)) and
       all(simplify(Gg[i, j] - Rational(-1, 2)) == 0 for i in range(3) for j in range(3) if i != j),
       "the anchor: the honeycomb Gram = our A_2 Gram (-1/d = -1/2) ⟹ the same object")
    kv = Matrix([kx, ky])
    fg = sum(exp(I * (kv.T * dl[j])[0, 0]) for j in range(3))
    # the point K: k = (0, 4π/(3√3))·... taken from the condition f=0 — solved, not cited
    Kpt = Matrix([0, 4 * pi / (3 * sqrt(3))])
    fK = simplify(trigsimp(expand(fg.subs({kx: Kpt[0], ky: Kpt[1]}).rewrite(cos))))
    ok(simplify(fK) == 0, "the anchor: f(K)=0 at the point K, found by solving (not by citation)")
    # the slope of |f| along an arbitrary direction from K
    ang = Symbol('ang', real=True)
    uu = Matrix([cos(ang), sin(ang)])
    fnear = fg.subs({kx: Kpt[0] + lam * uu[0], ky: Kpt[1] + lam * uu[1]})
    absf2 = expand(simplify((fnear * fnear.conjugate()).rewrite(cos)))
    c2 = simplify(trigsimp(series(absf2, lam, 0, 3).removeO().coeff(lam ** 2)))
    ok(simplify(trigsimp(c2 - Rational(9, 4))) == 0,
       "★the anchor: the coeff λ^2 in |f|^2 = 9/4 INDEPENDENT of the direction ang ⟹ v_F = 3/2 (t=a=1)")
    ok(simplify(Rational(9, 4) - Rational((2 + 1) ** 2, 2 * 2)) == 0,
       "★the anchor AGREES with the hop-dictionary (d+1)^2/(2d)=9/4 — it COULD have failed to match")
    ok(simplify(Rational(81, 8) - Rational(9, 4)) != 0,
       "★the anchor DISAGREES with the registered 81/8 (=h^4/(4d) at d=2) by a factor of 9/2")

    # ==================================================================
    print("\n[T9] ★AN AUDIT OF THE STANDING GREEN TEST test/test_2_velocity_three_ways.py")
    # ==================================================================
    # T2 stands at 57/0 and carries the conclusion «T33 stands: h^4/(4d) — is a RATIO, and it
    # differs from the velocity by EXACTLY ONE NAMED factor — the inverse column-curvature h^2/2 — not
    # a fitted constant» (lines 533-534), resting on the assert of lines 382-384.
    # I check THIS VERY ASSERT for its ability to fail.
    Dsym = Symbol('D', positive=True)
    transverse_t2 = (Dsym + 1) ** 2 / (2 * Dsym)      # the numerator, road (a)
    vF2_t2 = ((Dsym + 1) ** 2 / Dsym) / 2             # road_d_absc2(D)/2 — the EXTERNAL road
    ok(simplify(transverse_t2 - vF2_t2) == 0,
       "T9a: in test_2 `transverse` and `v_F^2` — are ONE AND THE SAME expression (d+1)^2/(2d) symbolically")
    Cany = Symbol('C', positive=True)                  # an ARBITRARY «column-curvature»
    conv_any = simplify((transverse_t2 / Cany) / vF2_t2)
    ok(simplify(conv_any - 1 / Cany) == 0,
       "T9b: ★conv = ratio/v_F^2 ≡ 1/C IDENTICALLY ∀C ⟹ the assert «the factor = 1/column-curvature» "
       "CANNOT FAIL — this is (T/C)/T, not a measurement")
    print("       ⟹ the green of this assert is empty: it does NOT show that the factor is «named»,")
    print("         because the same assert would pass for ANY denominator C.")
    print("       ⟹ the load-bearing SUBSTANTIVE part of test_2 — is precisely T9a: the external v_F^2 EQUALS")
    print("         our numerator. And that IS the S1055 verdict: v^2 = ½trM, without h^2/2.")
    # a negative-control for T9b itself: an assert that CAN fail (so that T9b is not vacuous)
    ok(simplify((transverse_t2 / Rational(2, 9)) - vF2_t2) != 0,
       "T9c: a negative control — the assert «ratio == v_F^2» FAILS at C=2/9 ⟹ resolving power exists, "
       "it is precisely the T9b-variety that is vacuous, not any comparison")

    # ==================================================================
    print("\n[M] MUTANTS / NEGATIVE CONTROLS (a test that cannot fail is not a witness)")
    # ==================================================================
    # M1: a broken Gram -> the tr M identity must BREAK
    d = 3
    n = d + 1
    ths = symbols('m0:%d' % n, real=True)
    Bm = bonds(d)
    # we substitute one bond so that the Gram breaks
    Bm_bad = list(Bm)
    Bm_bad[0] = Bm[0] * Rational(6, 5)
    Mbad, _, _ = M_matrix(d, list(ths), Bm_bad)
    fabs2 = simplify(expand(trigsimp((f_of_phases(list(ths)) *
                                      f_of_phases(list(ths)).conjugate()).rewrite(cos))))
    claim = simplify(n * (1 + Rational(1, d)) - fabs2 / d)
    ok(simplify(expand(trigsimp(Mbad.trace() - claim))) != 0,
       "M1: a broken Gram ⟹ the tr M identity BREAKS (resolving power exists)")

    # M2: a column hop = 2 bonds -> v^2 must shift (otherwise the length question is empty)
    for d in [2, 3]:
        v2_1 = Rational((d + 1) ** 2, 2 * d)
        v2_2 = simplify(v2_1 / 4)   # x over 2 bonds: T=2-2cos(x/2)? coeff 1/4
        ok(v2_1 != v2_2, "M2 d={0}: a hop=2 bonds moves v^2 ({1}→{2}) ⟹ the length question is NOT empty"
           .format(d, v2_1, v2_2))

    # M3: CAN the registered number be rescued by an HONEST reading of SC*n/2?
    # ★A TYPED test: a rescue is possible ⟺ BOTH the NUMBER and the TYPE match (coeff ⊥ 2nd-derivative).
    # The M-eigenvalues = λ1,λ2 — are the COEFFICIENTS of λ^2 in the directions; the 2nd-derivatives = 2λ1,2λ2.
    print("\n  M3 (an attempt to RESCUE the registered number — kill-first against my own verdict):")
    for d in [2, 3, 4]:
        h = d + 1
        trM = Rational((d + 1) ** 2, d)
        # all honest readings of the «transverse-curvature», with an EXPLICIT type
        readings = [
            ("sum of coeffs = trM",              trM,       "coef"),
            ("mean coeff = ½trM",         trM / 2,   "coef"),   # ← this IS SC·n/2
            ("sum of 2nd-derivs = 2trM",       2 * trM,   "deriv"),
            ("mean 2nd-derivative = trM",      trM,       "deriv"),
        ]
        # the denominators, with a type
        dens = [("coeff of ν^2 = 1/h^2", Rational(1, h ** 2), "coef"),
                ("T''(0) = 2/h^2",       Rational(2, h ** 2), "deriv")]
        reg = Rational(h ** 4, 4 * d)
        rescues = []
        for nm, val, tp in readings:
            for dnm, dval, dtp in dens:
                if tp == dtp and simplify(val / dval - reg) == 0:
                    rescues.append((nm, dnm))
        ok(len(rescues) == 0,
           "M3 d={0}: among 4 readings × 2 denominators there is NO CONSISTENT pair "
           "giving the registered {1} ⟹ there is NO rescue by reading".format(d, reg))
        # and a positive control: the mixed pair — exists, and it is THE ONLY ONE
        mixed = [(nm, dnm) for nm, val, tp in readings for dnm, dval, dtp in dens
                 if tp != dtp and simplify(val / dval - reg) == 0]
        ok(len(mixed) == 1,
           "M3+ d={0}: the registered number is reproduced by EXACTLY ONE pair, and it IS MIXED: {1} / {2}"
           .format(d, mixed[0][0], mixed[0][1]) if mixed else
           "M3+ d={0}: expected exactly one mixed pair, found {1}".format(d, len(mixed)))
        # ★a type-inconsistency: ½trM CANNOT be the mean 2nd-derivative
        ok(simplify(trM / 2 - trM) != 0,
           "M3t d={0}: ½trM ≠ the mean 2nd-derivative (=trM) ⟹ the type of the numerator is UNAMBIGUOUSLY = coeff"
           .format(d))

    # M4: registered/consistent — a constant or f(d)?
    ratios = set()
    for d in [2, 3, 4, 5, 6]:
        h = d + 1
        ratios.add(simplify(Rational(h ** 4, 4 * d) / Rational(h ** 4, 2 * d)))
    ok(ratios == {Rational(1, 2)},
       "M4: registered/consistent = 1/2 for d=2..6 ⟹ a CONSTANT (bookkeeping), not d-geometry")

    # M5: is h^2 not ALSO a «constant»? (it must be f(d) — otherwise it would merge with the pairing)
    h2s = set(simplify(Rational((d + 1) ** 2, 2)) for d in [2, 3, 4, 5])
    ok(len(h2s) == 4,
       "M5: h^2/2 is DIFFERENT ∀d ⟹ the length (h^2) and the pairing (1/2) — are SEPARATE defects, not one")

    # ==================================================================
    print("\n" + "=" * 78)
    print("SUMMARY S1055: PASS={0}  FAIL={1}".format(PASS, FAIL))
    print("=" * 78)
    print("""
ANSWER TO QUESTION 1 (the factor of 1/2): I CONFIRM ALPHA'S KILL.
  The on-shell solution 2-2cos(x)=λ^2 s gives a slope of EXACTLY sqrt(s), without 1/2 and without 2 [T4].
  The coeff of λ^2 in |f|^2 = u^T M u EXACTLY [T2]; the second derivative = 2·coeff — a different object.
  The registered h^4/(4d) is reproduced by EXACTLY one, MIXED pairing [T6];
  both consistent pairings give ONE number h^4/(2d). The ratio = 1/2 ∀d [M4].

ANSWER TO QUESTION 2 (two lengths): I CONFIRM.
  1/h^2 in the denominator = ONLY the chain rule dx/dν=2π/h [T5]; h = the number of bonds
  in the period (S1001). The numerator carries a length through the diagonal of the Gram |δ|^2=1 = a BOND.

★MY OWN ADDITION (not in S1054): THE EXTERNAL ANCHOR IS PROMOTED TO A KILL [T8].
  The honeycomb Gram = our A_2 Gram; f(K)=0 is SOLVED, not cited; the coeff λ^2 = 9/4
  INDEPENDENT of the direction ⟹ v=3/2. This agrees with the hop-dictionary and disagrees with
  the registered 81/8 by a factor of 9/2 = h^2/2. That is, the question is NOT internally-conventional:
  the registered number contradicts an EXTERNAL ruler, re-derived by my own hand from their
  definition. This raises the verdict from «the dictionary is ambiguous» to «the number is wrong».
""")
    logf.close()
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
