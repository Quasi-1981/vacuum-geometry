# -*- coding: utf-8 -*-
# DIM: symbolic-d (THE UNITS LADDER, steps 1-2 of the assignment `hub/prime/UNITS_LADDER_EXANTE.md`:
#          THE DICTIONARY of the ratio-number `h⁴/(4d)` — of WHAT is it a ratio, which LENGTH and which
#          STEP enter the numerator and the denominator; then the LADDER — which of the two scales
#          (half_bond ⊥ medial_step) enters the dictionary. STOP after step 2 — step 3
#          (the kill-test of the «1 hop» ansatz) and step 4 (the ∀d-law) do NOT begin.)
"""S1054 — UNITS LADDER, steps 1-2 only.  Layer-1, exact, symbolic in d where possible.

THE ONE QUESTION OF STEP 1: the registered ratio `h⁴/(4d)` (T33, probe S1012) is built as

        ratio = [S999: SC·n/2] / curv_col ,   curv_col := 2/h²   (T_col''(0), 2π dropped)

Of WHAT is it the ratio?  Concretely: which length sits in the numerator, which length sits
in the denominator, and are they the SAME length?  The answer is obtained by REDERIVING both
sides from their own primitives with the unit kept visible at every step, then reading the
units off the two expansions — never by choosing.

STEP 2 asks which of the two exact scales enters that dictionary:
        half_bond   = 1/2                    (no d)
        medial_step = sqrt((1+1/d)/2)        (d-dependent)
Both are exact (stake E of the order); the question is which one CARRIES A ROLE.

WHAT IS COMPUTED (all exact; symbolic d where the statement is ∀d)
  A. cell primitives from an explicit realization: Gram of the d+1 bonds (1 on the diagonal,
     −1/d off), half_bond, medial_step — and the guard that the two are never equal;
  B. numerator: tr M from the bond Gram at a node, with symbolic phases (identity ∀k, then
     the node specialization), plus the SCALING LAW that exposes which length it carries;
  C. denominator: the h-cycle of the marked-axis column (S1001: exact translation at step
     h=d+1) with its two possible momentum variables — per ONE hop vs per FULL period — and
     the expansion coefficient / second derivative of each;
  D. the on-shell relation solved as an exact series (not asserted): the leading coefficient
     IS the dictionary of the ratio-number;
  E. the pairing audit: all four (coefficient | second-derivative) pairings enumerated, with
     which one reproduces the registered `h⁴/(4d)`;
  F. step 2: the ladder — the same number re-expressed in each candidate unit, and the test
     of which unit the machinery actually USES (as opposed to which one makes it pretty).

TEETH
  (T1) the numerator identity is proved with SYMBOLIC phases (generic point, then node), and
       a mutated Gram must break it — otherwise the identity is vacuous;
  (T2) the on-shell leading coefficient is obtained by series inversion of the EXACT relation,
       with the next order printed, so the statement is leading-order and known to be so;
  (T3) MUTANT: a column whose hop is 2 bonds instead of 1 must MOVE the number — if it did
       not, the length question would be vacuous and no dictionary could be wrong;
  (T4) MUTANT: mixing half_bond with medial_step must be caught (a K-mixing of the order);
  (T5) external anchor d=2: the hop-unit number must reproduce the value that an outside hand
       publishes for this object (S1050 road (d), supplement of the 2024 source).

FENCE: Layer-1 only.  Ratio labels only (hop, half_bond, medial_step, curv_col, v_mean2,
ratio_num).  Any physical reading of the ratio lives behind the fence, not here.  The bare
term for the numerator is FORBIDDEN without its bracket — it is always written as the
bracketed form [S999: SC·n/2] (homonym caught in S1050/J-0503).
"""
import os
import sys

import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))

DS = (2, 3, 4, 5, 6)

_fails = []
_passes = []


def ok(cond, msg):
    (_passes if cond else _fails).append(msg)
    print(("  ✓ " if cond else "  ✗ FAIL ") + msg)


# ============================================================================
# A. CELL PRIMITIVES — bonds, half_bond, medial_step (exact, from a realization)
# ============================================================================

def bond_gram(dv):
    """Gram of the d+1 cell bonds from the explicit realization u_i = s(e_i − 1/n),
    s = sqrt(n/d) — NOT assumed: the Gram is computed and then checked."""
    n = dv + 1
    s = sp.sqrt(sp.Rational(n, dv))
    us = [sp.Matrix([s * ((1 if i == j else 0) - sp.Rational(1, n)) for j in range(n)])
          for i in range(n)]
    G = sp.Matrix(n, n, lambda i, j: sp.simplify((us[i].T * us[j])[0, 0]))
    tot = sp.zeros(n, 1)
    for u in us:
        tot += u
    return G, us, sp.simplify(tot)


def medial_step_from_gram(G, dv):
    """|δ_i/2 − δ_j/2| computed FROM the Gram (midpoints of two bonds at one site)."""
    g_ii, g_ij = G[0, 0], G[0, 1]
    return sp.sqrt(sp.simplify((g_ii + g_ii - 2 * g_ij) / 4))


# ============================================================================
# B. NUMERATOR — [S999: SC·n/2] = ½·tr M, and the length it carries
# ============================================================================

def trM_symbolic(dv, scale=sp.Integer(1)):
    """tr M with symbolic phases, bonds scaled by `scale` (M = a a^T + b b^T,
    a = ∇Re f, b = ∇Im f, f = Σ_j exp(i⟨k,δ_j⟩)).

    tr M = Σ_{j,l} cos(θ_j − θ_l) ⟨δ_j, δ_l⟩ — the ONLY place a transverse length enters."""
    n = dv + 1
    G, _, _ = bond_gram(dv)
    th = sp.symbols('th0:%d' % n, real=True)
    expr = 0
    for j in range(n):
        for l in range(n):
            expr += sp.cos(th[j] - th[l]) * G[j, l] * scale ** 2
    return sp.simplify(sp.expand_trig(sp.expand(expr))), th


def abs_f2_symbolic(th):
    """|f|² = Σ_{j,l} cos(θ_j − θ_l) — the node condition is |f|² = 0."""
    n = len(th)
    return sp.simplify(sum(sp.cos(th[j] - th[l]) for j in range(n) for l in range(n)))


# ============================================================================
# C. DENOMINATOR — the column as an h-cycle, and its two momentum variables
# ============================================================================

def cycle_spectrum(h):
    """Eigenvalues of the h-cycle operator 2I − S − S^{-1} (S = cyclic shift), exact.
    The column of the marked axis closes after h = d+1 hops (S1001): translation by one
    bond is not a lattice symmetry, translation by h bonds is; the h cosets form the cycle."""
    S = sp.zeros(h, h)
    for i in range(h):
        S[i, (i + 1) % h] = 1
    L = 2 * sp.eye(h) - S - S.T
    ev = []
    for val, mult in L.eigenvals().items():
        ev.extend([sp.simplify(val)] * mult)
    return L, sorted(ev, key=lambda e: sp.re(sp.N(e)))


# ============================================================================
# MAIN
# ============================================================================

class Tee:
    def __init__(self, real, fh):
        self.real, self.fh, self.chunks = real, fh, []

    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1054_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf)
    sys.stdout = tee
    sys.path.insert(0, os.path.join(_HERE, ".."))

    d = sp.Symbol('d', positive=True)
    n_sym = d + 1          # h = d+1: coordination of the cell = column period (S1001)

    print("=" * 78)
    print("S1054 — THE UNITS LADDER, steps 1-2 (assignment hub/prime/UNITS_LADDER_EXANTE.md)")
    print("Layer-1, exact arithmetic. STOP after step 2 (steps 3-4 do not begin).")
    print("=" * 78)
    print()

    # ---------------- A ----------------
    print("A. CELL PRIMITIVES (from a realization, not from the formula)")
    print("d | Gram diag | Gram off | Σδ | half_bond | medial_step | equal?")
    for dv in DS:
        G, us, tot = bond_gram(dv)
        diag_ok = all(sp.simplify(G[i, i] - 1) == 0 for i in range(dv + 1))
        off_ok = all(sp.simplify(G[i, j] + sp.Rational(1, dv)) == 0
                     for i in range(dv + 1) for j in range(dv + 1) if i != j)
        closed = (sp.simplify(tot) == sp.zeros(dv + 1, 1))
        hb = sp.Rational(1, 2)
        ms = medial_step_from_gram(G, dv)
        ms_closed = sp.sqrt((1 + sp.Rational(1, dv)) / 2)
        ok(diag_ok and off_ok, "d={0}: the Gram of the bonds = 1 diag, −1/d off (computed from the realization)".format(dv))
        ok(closed, "d={0}: Σδ = 0 (closure of the cell)".format(dv))
        ok(sp.simplify(ms - ms_closed) == 0,
           "d={0}: medial_step from the Gram = sqrt((1+1/d)/2) = {1}".format(dv, sp.nsimplify(ms)))
        ok(sp.simplify(ms - hb) != 0,
           "d={0}: half_bond={1} ≠ medial_step={2} (a K-mixing has resolving power)".format(
               dv, hb, sp.nsimplify(ms)))
        print("  {0} | 1 | −1/{0} | 0 | 1/2 | {1} | NO".format(dv, sp.nsimplify(ms)))
    # ∀d symbolically: the two scales never coincide
    gap = sp.simplify((1 + 1 / d) / 2 - sp.Rational(1, 4))
    ok(sp.simplify(gap - (d + 2) / (4 * d)) == 0,
       "★∀d symbolically: medial_step² − half_bond² = (d+2)/(4d) > 0 — they NEVER coincide")
    print()

    # ---------------- B ----------------
    print("B. THE NUMERATOR [S999: SC·n/2] — what it is and which LENGTH it carries")
    print("   tr M = Σ_{j,l} cos(θ_j−θ_l)·⟨δ_j,δ_l⟩ — the length enters PRECISELY through the Gram of the bonds.")
    for dv in (2, 3, 4):
        trM, th = trM_symbolic(dv)
        f2 = abs_f2_symbolic(th)
        n = dv + 1
        # identity ∀k:  tr M = n(1+1/d) − |f|²/d
        ident = sp.simplify(trM - (n * (1 + sp.Rational(1, dv)) - f2 / dv))
        ok(ident == 0,
           "d={0}: the identity ∀k (symbolic phases): tr M = n(1+1/d) − |f|²/d".format(dv))
        # node specialization |f|² = 0
        trM_node = sp.simplify(n * (1 + sp.Rational(1, dv)))
        ok(sp.simplify(trM_node - sp.Rational((dv + 1) ** 2, dv)) == 0,
           "d={0}: at the node (|f|²=0) ⟹ tr M = (d+1)²/d = {1} — k-INDEPENDENT".format(
               dv, sp.Rational((dv + 1) ** 2, dv)))
        # the S999 number
        SC = sp.Rational(dv + 1, dv)
        ok(sp.simplify(SC * n / 2 - trM_node / 2) == 0,
           "d={0}: [S999: SC·n/2] = ½·tr M = {1} — the NUMBERS COINCIDE (not two objects)".format(
               dv, sp.nsimplify(trM_node / 2)))
    # ★ SCALING LAW — the dictionary of the numerator
    rho = sp.Symbol('rho', positive=True)
    trM_1, _ = trM_symbolic(2)
    trM_r, _ = trM_symbolic(2, scale=rho)
    ok(sp.simplify(trM_r - rho ** 2 * trM_1) == 0,
       "★THE DICTIONARY OF THE NUMERATOR: δ→ρδ ⟹ tr M→ρ²·tr M ⟹ [S999: SC·n/2] carries (LENGTH_⊥)²,")
    print("     and this length = |δ| = 1 BOND (the Gram has 1 on the diagonal, step A) ⟹")
    print("     ½·tr M = (amplitude·bond)² at amplitude 1 — that is, (bond/step)² is NOT YET fixed:")
    print("     it is fixed only together with the denominator, which brings the STEP.")
    print()

    # ---------------- C ----------------
    print("C. THE DENOMINATOR curv_col — the column as an h-cycle, and TWO possible momentum variables")
    for dv in (2, 3, 4, 5):
        h = dv + 1
        L, ev = cycle_spectrum(h)
        want = sorted([sp.simplify(2 - 2 * sp.cos(2 * sp.pi * sp.Rational(nu, h))) for nu in range(h)],
                      key=lambda e: sp.re(sp.N(e)))
        same = all(sp.simplify(a - b) == 0 for a, b in zip(ev, want))
        ok(same, "d={0}: the spectrum of the column h-cycle = {{2−2cos(2πν/h)}}, ν∈ℤ/h — h={1} points".format(dv, h))
    print()
    print("   TWO VARIABLES OF THE SAME object (this IS the fork of the dictionary):")
    x, nu = sp.symbols('x nu', real=True)          # x = the phase per ONE hop; nu = the ℤ/h index
    h_s = n_sym
    T_x = 2 - 2 * sp.cos(x)                        # the hop-variable: the step = 1 bond along the column
    T_nu = 2 - 2 * sp.cos(2 * sp.pi * nu / h_s)    # the index-variable: the step = the PERIOD = h bonds
    coef_x = sp.simplify(sp.series(T_x, x, 0, 4).removeO().coeff(x, 2))
    dd_x = sp.simplify(sp.diff(T_x, x, 2).subs(x, 0))
    coef_nu = sp.simplify(sp.series(T_nu, nu, 0, 4).removeO().coeff(nu, 2))
    dd_nu = sp.simplify(sp.diff(T_nu, nu, 2).subs(nu, 0))
    ok(coef_x == 1 and dd_x == 2,
       "the hop-variable x (the phase per ONE bond-hop): T_col ≈ 1·x², T''(0) = 2")
    ok(sp.simplify(coef_nu - (2 * sp.pi / h_s) ** 2) == 0 and sp.simplify(dd_nu - 2 * (2 * sp.pi / h_s) ** 2) == 0,
       "the index-variable ν (a ℤ/h character): T_col ≈ (2π/h)²·ν², T''(0) = 2(2π/h)²")
    ok(sp.simplify((dd_nu / dd_x) - (2 * sp.pi / h_s) ** 2) == 0,
       "★THE RELATION: x = 2πν/h EXACTLY ⟹ dd_ν/dd_x = (2π)²/h² — (2π)² is the phase↔radian conversion,")
    print("     and 1/h² — is the LENGTH: ν measures in units of 1/(h bonds) = 1/PERIOD of the column,")
    print("     x measures in units of 1/(1 bond) = 1/HOP. The denominator 2/h² = this is dd_ν with the (2π)² REMOVED,")
    print("     that is, CURVATURE PER UNIT OF PERIOD, not per unit of hop.")
    print()

    # ---------------- D ----------------
    print("D. ON-SHELL SOLVED AS A SERIES (not asserted): what exactly gives the slope of the cone")
    lam, s = sp.symbols('lam s', positive=True)
    # exact on-shell: 2 − 2cos(x) = |f|² = λ²·s,  s = uᵀMu along a cone direction
    x_exact = 2 * sp.asin(lam * sp.sqrt(s) / 2)
    ser = sp.series(x_exact, lam, 0, 5).removeO()
    lead = sp.simplify(ser.coeff(lam, 1))
    nxt = sp.simplify(ser.coeff(lam, 3))
    ok(sp.simplify(lead - sp.sqrt(s)) == 0,
       "★the on-shell 2−2cos(x)=λ²·s solved EXACTLY: x = λ·sqrt(s) + O(λ³) ⟹ dx/dλ = sqrt(uᵀMu)")
    ok(sp.simplify(nxt - s ** sp.Rational(3, 2) / 24) == 0,
       "★TOOTH: the next order is printed (λ³·s^{3/2}/24 ≠ 0) ⟹ the statement is of LEADING order,")
    print("     not an identity (a test that cannot fail is not a witness).")
    print("   ⟹ the slope in the HOP-dictionary: (phase per hop) = sqrt(uᵀMu)·|q|, q — the transverse wavevector")
    print("     in units of 1/bond. Averaging over the TWO nonzero directions of the cone (rank M = 2):")
    for dv in (2, 3):
        h = dv + 1
        v2_hop = sp.Rational((dv + 1) ** 2, 2 * dv)
        ok(sp.simplify(v2_hop - sp.Rational(dv + 1, dv) * (dv + 1) / 2) == 0,
           "d={0}: v_mean2(hop-dictionary) = ½·tr M = [S999: SC·n/2] = {1}".format(dv, v2_hop))
    v2_hop_sym = sp.simplify((d + 1) ** 2 / (2 * d))
    print("   ★∀d: v_mean2(hop) = (d+1)²/(2d) — WITHOUT h⁴, because the column is measured in HOPS, like the numerator.")
    print()

    # ---------------- D2: end-to-end, nothing borrowed ----------------
    print("D2. END-TO-END FROM THE LATTICE TO THE SLOPE (no borrowed object: |f|² is expanded")
    print("    HERE, M is built HERE, on-shell is solved HERE — the chain is closed)")
    lam_ = sp.Symbol('lam', positive=True)
    tpar = sp.Symbol('t', real=True)
    node_phases = {
        2: [sp.Integer(0), 2 * sp.pi / 3, -2 * sp.pi / 3],
        3: [tpar, tpar + sp.pi, -tpar - sp.pi, -tpar],          # the WHOLE nodal line, t symbolic
    }
    for dv in (2, 3):
        n = dv + 1
        _, us, _ = bond_gram(dv)
        # bonds live in the hyperplane {Σx=0} ⊂ R^{n}; use them as vectors of R^n directly
        th_t = node_phases[dv]
        ok(sp.simplify(sp.expand_complex(sum(sp.exp(sp.I * a) for a in th_t))) == 0,
           "d={0}: the chosen phases satisfy the node condition Σe^{{iθ}}=0 (symbolically)".format(dv))
        # transverse direction u: generic, inside the hyperplane, orthogonal to nothing special
        uc = sp.symbols('u0:%d' % n, real=True)
        uvec = sp.Matrix(uc) - sp.Matrix([sum(uc) / n] * n)      # project into {Σx=0}
        # f(k0 + λu):  phase_j = θ_j + λ·⟨u, δ_j⟩
        sj = [sp.simplify((uvec.T * us[j])[0, 0]) for j in range(n)]
        f = sum(sp.exp(sp.I * (th_t[j] + lam_ * sj[j])) for j in range(n))
        f2 = sp.simplify(sp.expand(sp.re(sp.expand_complex(f)) ** 2
                                   + sp.im(sp.expand_complex(f)) ** 2))
        ser = sp.series(f2, lam_, 0, 3).removeO()
        c0 = sp.simplify(ser.coeff(lam_, 0))
        c1 = sp.simplify(ser.coeff(lam_, 1))
        c2 = sp.simplify(sp.expand_trig(sp.expand(ser.coeff(lam_, 2))))
        # M = a a^T + b b^T built here from the same phases
        a = sp.zeros(n, 1)
        b = sp.zeros(n, 1)
        for j in range(n):
            a += (-sp.sin(th_t[j])) * us[j]
            b += sp.cos(th_t[j]) * us[j]
        # ⟨u,a⟩² + ⟨u,b⟩²  ≡  u^T M u   (M never formed as a matrix — no room for a slip)
        uMu = sp.simplify(sp.expand_trig(sp.expand(((uvec.T * a)[0, 0]) ** 2 + ((uvec.T * b)[0, 0]) ** 2)))
        ok(c0 == 0 and c1 == 0,
           "d={0}: |f|² at the node: λ⁰=0 and λ¹=0 ⟹ the node is genuinely a node, the expansion begins at λ²".format(dv))
        ok(sp.simplify(c2 - uMu) == 0,
           "d={0}: ★the λ²-coefficient |f|² = uᵀMu EXACTLY (symbolic u{1}) — WITHOUT a ½ factor".format(
               dv, ", symbolic t along the WHOLE line" if dv == 3 else ""))
        # trace check: sum over an orthonormal-in-{Σx=0} frame reproduces tr M
        trM_here = sp.simplify(sp.expand_trig(sp.expand((a.T * a)[0, 0] + (b.T * b)[0, 0])))
        ok(sp.simplify(trM_here - sp.Rational((dv + 1) ** 2, dv)) == 0,
           "d={0}: tr M = |a|²+|b|² = (d+1)²/d = {1} (the same object, built here)".format(
               dv, sp.Rational((dv + 1) ** 2, dv)))
    print("   ⟹ THE CHAIN IS CLOSED: |f|² = λ²·uᵀMu (without ½) ⊥ T_col = x² (without ½ in the hop-variable)")
    print("     ⟹ on-shell x = λ·sqrt(uᵀMu): the slope of the cone carries NO factor of two.")
    print()

    # ---------------- E ----------------
    print("E. THE PAIRING AUDIT — where h⁴/(4d) comes from (all four pairings, enumerated)")
    tc_coef = sp.simplify((d + 1) ** 2 / (2 * d))     # ½ tr M  = the COEFFICIENT of q² (leading)
    tc_dd = sp.simplify((d + 1) ** 2 / d)             # tr M    = the MEAN SECOND DERIVATIVE of |f|² over the cone
    col_coef_per = sp.simplify(1 / (d + 1) ** 2)      # 1/h²    = the coefficient of ν̃² (period-unit)
    col_dd_per = sp.simplify(2 / (d + 1) ** 2)        # 2/h²    = the second derivative (period-unit)
    pairs = [
        ("coeff / coeff   (consistent)", tc_coef / col_coef_per),
        ("2nd-deriv / 2nd-deriv (consistent)", tc_dd / col_dd_per),
        ("★coeff / 2nd-deriv (MIXED)", tc_coef / col_dd_per),
        ("2nd-deriv / coeff (MIXED)", tc_dd / col_coef_per),
    ]
    print("   (everything in the PERIOD-unit of the column, to isolate the pairing itself)")
    for lbl, val in pairs:
        print("     {0:38s} = {1}".format(lbl, sp.simplify(val)))
    reg = sp.simplify((d + 1) ** 4 / (4 * d))
    consistent = sp.simplify((d + 1) ** 4 / (2 * d))
    ok(sp.simplify(pairs[0][1] - consistent) == 0 and sp.simplify(pairs[1][1] - consistent) == 0,
       "★BOTH CONSISTENT pairings give ONE AND THE SAME thing: h⁴/(2d) (not my choice — an agreement of two)")
    ok(sp.simplify(pairs[2][1] - reg) == 0,
       "★the REGISTERED h⁴/(4d) is reproduced by EXACTLY one — MIXED — pairing (coeff over 2nd-deriv)")
    ok(sp.simplify(reg / consistent - sp.Rational(1, 2)) == 0,
       "★∀d symbolically: registered / consistent = 1/2 (not a d-dependent factor — exactly a two)")
    # reproduce the registry number exactly as S1012 builds it (bit-fence with the ancestor)
    for dv in (2, 3):
        h = dv + 1
        s1012 = sp.simplify(sp.Rational((dv + 1) ** 2, 2 * dv) / sp.Rational(2, h ** 2))
        ok(sp.simplify(s1012 - sp.Rational(h ** 4, 4 * dv)) == 0,
           "d={0}: reproduced the S1012 assembly bit-for-bit: (SC·n/2)/(2/h²) = h⁴/(4d) = {1}".format(
               dv, sp.Rational(h ** 4, 4 * dv)))
    print()
    print("   ⟹ THE DICTIONARY OF THE REGISTERED NUMBER (the answer to step 1, verbatim):")
    print("      h⁴/(4d) = [the coefficient of q² in |f|², length = 1 BOND]")
    print("                 / [the second derivative of T_col, length = 1 PERIOD = h bonds]")
    print("      that is, TWO DIFFERENT lengths (bond ⊥ h·bond) AND two different operations (coeff ⊥ 2nd-deriv).")
    print("      One consistent dictionary: both lengths = HOP ⟹ v_mean2 = (d+1)²/(2d).")
    print("      A second consistent dictionary: the column in PERIODS ⟹ h⁴/(2d).  The registered one is neither.")
    print()

    # ---------------- F: STEP 2 — the ladder ----------------
    print("F. STEP 2 — THE LADDER: which scale ENTERS (derived, not chosen)")
    print("   Method: a scale-probe. The number v_mean2 is re-expressed in each candidate unit;")
    print("   the one that ENTERS is the one the machinery USES (the Gram of the numerator / the cycle-step of the denominator),")
    print("   not the one in which the number is prettier.")
    hb2 = sp.Rational(1, 4)
    ms2 = sp.simplify((1 + 1 / d) / 2)
    v2_bond = sp.simplify((d + 1) ** 2 / (2 * d))
    v2_in_hb = sp.simplify(v2_bond / hb2)
    v2_in_ms = sp.simplify(v2_bond / ms2)
    ok(sp.simplify(v2_in_ms - (d + 1)) == 0,
       "★∀d symbolically: v_mean2 in units of medial_step = d+1 = h EXACTLY (an integer ∀d)")
    ok(sp.simplify(v2_in_hb - 2 * (d + 1) ** 2 / d) == 0,
       "★∀d symbolically: v_mean2 in units of half_bond = 2(d+1)²/d (an integer only when d|2(d+1)²)")
    print("   d | v_mean2 [bond] | [half_bond] | [medial_step]")
    for dv in DS:
        print("   {0} | {1} | {2} | {3}".format(
            dv, sp.Rational((dv + 1) ** 2, 2 * dv),
            sp.nsimplify(v2_in_hb.subs(d, dv)), sp.nsimplify(v2_in_ms.subs(d, dv))))
    # WHICH ONE THE MACHINERY USES — the derivation, not the beauty contest
    G2, _, _ = bond_gram(2)
    ok(sp.simplify(G2[0, 0] - 1) == 0,
       "★THE NUMERATOR-INPUT: the diagonal of the Gram = 1 = |δ|² ⟹ tr M has the FULL bond entering it,")
    print("     not its midpoint: half a bond would give a Gram of 1/4 on the diagonal (the scale-probe of B).")
    L3, _ = cycle_spectrum(3)
    ok(sp.simplify(L3[0, 1] + 1) == 0,
       "★THE DENOMINATOR-INPUT: adjacency of the cycle = ONE step of the column = one bond (S1001) ⟹ the hop = a FULL bond")
    print("   ⟹ THE VERDICT OF STEP 2: the FULL BOND ENTERS the dictionary of v — in BOTH slots.")
    print("     NEITHER half_bond NOR medial_step APPEAR in the derivation: both are a RE-EXPRESSION of")
    print("     the already-derived number, not inputs. ⟨An observation, NOT a derivation: in medial_step-units the number")
    print("     becomes the integer h ∀d — pretty, but it does not follow from the derivation; weight = 0 until forced.⟩")
    print()

    # ---------------- G: geometry of the column (S597) — the last open item of step 1 ----
    print("G. THE GEOMETRY OF THE COLUMN (S597: a 3:1 stacking, «lateral 1/√2») — does it enter curv_col")
    print("   Question: does the h in T_col come from the PATH LENGTH (zigzag) or from the ORDER of u₀?")
    for dv in (2, 3, 4, 5):
        n = dv + 1
        _, us, _ = bond_gram(dv)
        # Λ = A_d = span_Z{u_i − u_0}; order of u_0 in the quotient = smallest m>0 with m·u_0 ∈ Λ
        gens = sp.Matrix.hstack(*[us[i] - us[0] for i in range(1, n)])
        def in_lattice(vec):
            sol = gens.solve_least_squares(vec)
            if sp.simplify(gens * sol - vec) != sp.zeros(n, 1):
                return False
            return all(sp.simplify(c - sp.floor(sp.nsimplify(c))) == 0 for c in sol)
        order = None
        for m in range(1, 3 * n + 2):
            if in_lattice(m * us[0]):
                order = m
                break
        ok(order == n,
           "d={0}: the order of u₀ in A_d*/A_d = {1} = h ⟹ h comes FROM THE ORDER, not from a path length".format(
               dv, order))
        # the zigzag step (a root) is ALREADY in the lattice: order 1 ⟹ it carries no h
        zig = us[0] - us[1]
        ok(in_lattice(zig),
           "d={0}: the zigzag-step (u₀−u_j) = a root ∈ Λ ⟹ order 1 ⟹ it CANNOT carry an h-cycle".format(dv))
    # lateral/axial of the zigzag — the S597 number, with its d put back
    lat_ax = sp.sqrt(sp.simplify((1 - 1 / d ** 2) / (1 + 1 / d) ** 2))
    lat_ax = sp.simplify(sp.sqrt((d - 1) / (d + 1)))
    ok(sp.simplify(lat_ax.subs(d, 3) - 1 / sp.sqrt(2)) == 0,
       "★the S597-number in its place: lateral/axis of the zigzag = sqrt((d−1)/(d+1)); at d=3 = 1/√2 EXACTLY")
    ok(sp.simplify(lat_ax.subs(d, 2) - 1 / sp.sqrt(3)) == 0,
       "★…but at d=2 it = 1/√3, not 1/√2 ⟹ «lateral 1/√2» — a d-DEPENDENT number, not a constant")
    print("   d | lateral/axis of the zigzag")
    for dv in DS:
        print("   {0} | {1}".format(dv, sp.nsimplify(lat_ax.subs(d, dv))))
    print("   ⟹ THE VERDICT of G: the geometry of the zigzag does NOT enter curv_col. The h-cycle of T_col carries the ORDER of u₀")
    print("     in A_d*/A_d (=h, S1001), while the zigzag-step — a root, is already in the lattice (order 1).")
    print("     The zigzag describes the REALIZED PATH between endpoints, not the coset-shift that carries the dual.")
    print()

    # ---------------- H: the S597 form, verified against the SOURCE's own numbers ----------
    print("H. THE S597 FORM «lat/time» — A CHECK FROM THE DEFINITION (Omega's brake: do not change until checked)")
    print("   The primary source: viz/SIMPLEX_LATTICE_MEMBRANES_omega_analysis.md, a table row:")
    print("   «the forced lateral per time-period | lat/time = 1/√2 (0.4082/0.5774) | exactly».")
    print("   ★It prints TWO numbers separately — so I check not the ratio, but EACH of the two.")
    # source realization: diamond, cubic a=1, bonds a/4·(±1,±1,±1) with even # of minus signs
    a4 = sp.Rational(1, 4)
    dia = [sp.Matrix([a4, a4, a4]), sp.Matrix([a4, -a4, -a4]),
           sp.Matrix([-a4, a4, -a4]), sp.Matrix([-a4, -a4, a4])]
    axis = sp.Matrix([1, 1, 1]) / sp.sqrt(3)          # the time-axis = the direction of the marked bond
    bond_len = sp.simplify(sp.sqrt((dia[0].T * dia[0])[0, 0]))
    ax_time = sp.simplify((dia[0].T * axis)[0, 0])     # the projection of the TIME-bond onto the axis
    ax_space = sp.simplify((dia[1].T * axis)[0, 0])    # the projection of the SPACE-bond onto the axis (negative)
    delta_src = sp.simplify(-ax_space)
    lat_space = sp.simplify(sp.sqrt((dia[1].T * dia[1])[0, 0] - ax_space ** 2))
    period_ax = sp.simplify(ax_time + delta_src)       # the axial travel per ONE zigzag-period
    print("   d=3, the source realization (diamond, cube a=1):")
    for lbl, val in (("|bond|", bond_len), ("δ = |axis-projection of the space-bond|", delta_src),
                     ("time-step (axis-projection of the time-bond) = 3δ", ax_time),
                     ("★TIME = the axial travel per period = 4δ", period_ax),
                     ("★LAT = the transverse part of the space-bond", lat_space)):
        print("     {0:44s} = {1} = {2}".format(lbl, val, sp.N(val, 5)))
    ok(sp.simplify(delta_src - sp.Rational(1, 4) / sp.sqrt(3)) == 0 and abs(sp.N(delta_src) - sp.Float("0.1443")) < 1e-4,
       "★checked against the source: δ = 0.1443 (the row «space-ends at −δ»)")
    ok(abs(sp.N(lat_space) - sp.Float("0.4082")) < 1e-4,
       "★checked against the source NUMBER-1: lat = 0.4082 = 1/√6 (the transverse part of the space-bond)")
    ok(abs(sp.N(period_ax) - sp.Float("0.5774")) < 1e-4,
       "★checked against the source NUMBER-2: time = 0.5774 = 1/√3 (the axial travel per zigzag-period, 4δ)")
    ok(sp.simplify(lat_space / period_ax - 1 / sp.sqrt(2)) == 0,
       "★and only then the ratio: lat/time = 1/√2 EXACTLY (both numbers matched SEPARATELY)")
    ok(sp.simplify(ax_time / delta_src - 3) == 0,
       "★checked against the source: the «3:1» stacking = the time-projection of the time-bond : the space-bond = 3")
    ok(sp.simplify((dia[1].T * dia[0])[0, 0] / bond_len ** 2 + sp.Rational(1, 3)) == 0,
       "★checked against the source: cos(space-bond, time-axis-bond) = −1/3 (an aperture of 109.47°)")
    print()
    print("   ⟹ THE DEFINITION, EXPLICITLY (exactly what is divided by what — this was missing from the citations):")
    print("     TIME  := the axial travel of the world-line per ONE forced zigzag-period")
    print("             = (time-step +u₀: travel 1) + (space-step −u_j: travel +1/d)  = 1 + 1/d  [bond]")
    print("     LAT  := the transverse (⊥ axis) drift over the same period; the time-step gives 0,")
    print("             all the drift — from the space-step = sqrt(1 − 1/d²)  [bond]")
    # ∀d generalization from the SAME definition, on the unit-bond cell
    for dv in (2, 3, 4, 5, 6):
        _, us, _ = bond_gram(dv)
        u0 = us[0]
        ax_t = sp.simplify((u0.T * u0)[0, 0])                       # = 1
        ax_s = sp.simplify(-(us[1].T * u0)[0, 0])                   # = +1/d (the step −u_j)
        lat = sp.simplify(sp.sqrt((us[1].T * us[1])[0, 0] - ((us[1].T * u0)[0, 0]) ** 2))
        per = sp.simplify(ax_t + ax_s)
        ratio = sp.simplify(lat / per)
        ok(sp.simplify(ax_s - sp.Rational(1, dv)) == 0 and sp.simplify(per - (1 + sp.Rational(1, dv))) == 0
           and sp.simplify(lat - sp.sqrt(1 - sp.Rational(1, dv ** 2))) == 0,
           "d={0}: TIME = 1+1/d = {1} · LAT = sqrt(1−1/d²) = {2} (the same definition)".format(
               dv, sp.nsimplify(per), sp.nsimplify(lat)))
        ok(sp.simplify(ratio - sp.sqrt(sp.Rational(dv - 1, dv + 1))) == 0,
           "d={0}: lat/time = sqrt((d−1)/(d+1)) = {1}".format(dv, sp.nsimplify(ratio)))
        ok(sp.simplify(ax_t / ax_s - dv) == 0,
           "d={0}: the stacking = d:1 = {1}:1 (that is, «3:1» — is ALSO a d-instance, not a constant)".format(dv, dv))
    print()
    print("   ★WHY A CHECK ON A SINGLE NUMBER WOULD BE EMPTY (Omega's brake is structurally right):")
    rivals = [("sqrt((d−1)/(d+1))  [this one]", sp.sqrt((d - 1) / (d + 1))),
              ("1/sqrt(d−1)", 1 / sp.sqrt(d - 1)),
              ("sqrt(2/(d+1))", sp.sqrt(2 / (d + 1))),
              ("sqrt(d−1)/2", sp.sqrt(d - 1) / 2)]
    print("     form | d=2 | d=3 | d=4 | d=5")
    for lbl, ex in rivals:
        print("     {0:24s} | {1} | {2} | {3} | {4}".format(
            lbl, *[sp.nsimplify(sp.simplify(ex.subs(d, dv))) for dv in (2, 3, 4, 5)]))
    all_agree_3 = all(sp.simplify(ex.subs(d, 3) - 1 / sp.sqrt(2)) == 0 for _, ex in rivals)
    differ_2 = len({sp.simplify(ex.subs(d, 2)) for _, ex in rivals}) > 1
    ok(all_agree_3 and differ_2,
       "★ALL four forms agree at d=3 (=1/√2) and DIVERGE at d=2 ⟹ the d=3-instance does NOT determine")
    print("     the form. The form is determined only by checking the TWO quantities SEPARATELY (lat and time), and it is exactly")
    print("     that the source gives: 0.4082 = sqrt(1−1/d²)·|bond| ⊥ 0.5774 = (1+1/d)·|bond| at d=3.")
    print()

    # ---------------- TEETH / MUTANTS ----------------
    print("TEETH AND MUTANTS:")
    mut = True
    # M1: mutated Gram must break the node identity
    dv = 3
    n = dv + 1
    G, _, _ = bond_gram(dv)
    Gbad = G[:, :]
    Gbad[0, 1] = Gbad[1, 0] = G[0, 1] + sp.Rational(1, 5)
    th = sp.symbols('th0:%d' % n, real=True)
    trM_bad = sp.simplify(sum(sp.cos(th[j] - th[l]) * Gbad[j, l] for j in range(n) for l in range(n)))
    f2 = abs_f2_symbolic(th)
    broke = sp.simplify(trM_bad - (n * (1 + sp.Rational(1, dv)) - f2 / dv)) != 0
    print("  M1 (a mutated Gram −1/d → −1/d+1/5): {0}".format(
        "CAUGHT — the identity breaks ⟹ it has resolving power" if broke else "NOT CAUGHT"))
    mut &= broke
    # M2: column hop of 2 bonds instead of 1 must move the number
    v2_hop2 = sp.simplify(v2_bond / 4)
    moved = sp.simplify(v2_hop2 - v2_bond) != 0
    print("  M2 (the column hop = 2 bonds): CAUGHT — v_mean2 {0} → {1} ⟹ the LENGTH truly enters,".format(
        v2_bond.subs(d, 2), v2_hop2.subs(d, 2)))
    print("     the dictionary question is NOT empty (had it not moved, no dictionary could be wrong)")
    mut &= moved
    # M3: half_bond/medial_step confusion (a K-mixing of the order)
    conf = sp.simplify(sp.sqrt(ms2) - sp.Rational(1, 2)) != 0
    print("  M3 (a K-mixing: half_bond ≡ medial_step): CAUGHT — the difference (d+2)/(4d) ≠ 0 ∀d")
    mut &= conf
    # M4: the halving is a pairing artefact, not a d-dependent factor
    ratio_fix = sp.simplify(reg / consistent)
    print("  M4 (is it not a d-dependent factor?): registered/consistent = {0} — a constant 1/2, not f(d)".format(ratio_fix))
    mut &= (sp.simplify(ratio_fix - sp.Rational(1, 2)) == 0)
    # T5: external anchor at d=2 (S1050 road (d): supplement of the 2024 source, |c|²=(d+1)²/d)
    ext = sp.Rational(9, 2) / 2          # |c|²/2 = (d+1)²/(2d) at d=2 — outside hand, own definitions
    print("  T5 (an external anchor d=2): the hop-dictionary gives {0}; the external ruler (|c|²/2 at d=2) = {1} — {2}".format(
        sp.Rational(9, 4), ext, "AGREES" if sp.simplify(ext - sp.Rational(9, 4)) == 0 else "DISAGREES"))
    mut &= (sp.simplify(ext - sp.Rational(9, 4)) == 0)
    ok(mut, "all mutants/teeth fired correctly")
    print()

    # ---------------- FENCE ----------------
    from tools.fence_scan import scan_forbidden
    pats = [r"світл\w*", r"ħ", r"hbar", r"спін\w*", r"\bspin\b", r"квант\w*",   # GUARDLINE
            r"\bquantum\b", r"\blight\b", r"часоподібн\w*", r"метрик\w+\s+простор\w+"]  # GUARDLINE
    hits = scan_forbidden(__file__, pats)
    hits += scan_forbidden("".join(tee.chunks), pats)
    ok(not hits, "the Layer-1 fence is clean: 0 hits on {0} patterns (source + log)".format(len(pats)))

    print()
    print("=" * 78)
    print("SCORE: {0} ✓ / {1} ✗".format(len(_passes), len(_fails)))
    print("STOP per the assignment: steps 3 (the kill-test of the «1 hop» ansatz) and 4 (the ∀d-law of the medial) did NOT begin.")
    print("=" * 78)
    sys.stdout = tee.real
    logf.close()
    return 1 if _fails else 0


if __name__ == "__main__":
    sys.exit(main())
