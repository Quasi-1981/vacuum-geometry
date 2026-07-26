# -*- coding: utf-8 -*-
# DIM: symbolic-exact (ℚ / algebraic).  Every number below is printed with a bracket
#      [address · unit · type/operation] — a requirement of Omega's assignment ef1dd6ad (2026-07-22).
"""S1059 — is the chain step INDIVISIBLE.  Layer-1, exact symbolic.  A binary output.

ASSIGNMENT: step 3 of `UNITS_LADDER_EXANTE.md` (the kill `K-action` — the kill's name as quoted), the exante of the probe —  # GUARDLINE
`hub/prime/S1059_CHAIN_STEP_INDIVISIBILITY_EXANTE.md` (carved BEFORE the count).

THE EDITION OF THE CARRIERS THIS IS WRITTEN AGAINST (assignment field №3):
  · `UNITS_LADDER_EXANTE.md` + the S1055 addendum (wagers A-E with Beta's status table);
  · `S1045_B1_CHESSBOARD_RESULT_ALPHA.md` + the S1058 errata-forward — that is, the standing
    `v² = ½·trM = (d+1)²/(2d)`; the old `h⁴/(4d)` is NOT used here at all (the probe
    does not need it: the GRADATION of the series is on trial, not the slope of the cone).

WHAT IS TESTED (the Layer-1-translation of wager A, the only one):
  in the expansion `G = Σ_n G₀·(H·G₀)ⁿ` (S1045, ∀d) the exponent `n` counts hops.
  Is `n=1` an indivisible and irreducible unit of gradation — or does the step-count
  carry a free parameter, the way LENGTH carried one before S1057?

THE BOUNDARY (carved in the exante BEFORE the numbers): the GRADATION OF THE CHAIN is on trial.
A Layer-2-reading of this gradation in the terms of the assignment is neither done nor claimed here.

WHY THIS IS NOT A VACUOUS TEST (the S1055::T9 lesson — an assert that cannot fail):
every detector is run on an object where it MUST say YES (N1, N2), and every
justification E1-E4 is computed explicitly, not declared.

★NAMES: the exante of this probe uses a homonym word in the everyday sense («the path of justification»).  # GUARDLINE
I do NOT retroactively edit the exante; in the code the field name is split off into «JUSTIFICATION E1-E4»,
so that the fence stays strict, not weakened for my own convenience.
"""
import os
import sys

import sympy as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
_fails, _passes = [], []


def ok(cond, msg):
    (_passes if cond else _fails).append(msg)
    print(("  ✓ " if cond else "  ✗ FAIL ") + msg)


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


# ─────────────────────────── shared objects ───────────────────────────
E, m0, fr, fi = sp.symbols('E m0 f_re f_im', real=True)
f = fr + sp.I * fi
fb = fr - sp.I * fi
absf2 = fr ** 2 + fi ** 2
I2 = sp.eye(2)
sz = sp.Matrix([[1, 0], [0, -1]])
H2c = sp.Matrix([[0, f], [fb, 0]])              # the two-component kernel of T33 (off-diag = one hop)
G0 = (E * I2 - m0 * sz).inv()                    # diag(1/(E-m0), 1/(E+m0)) — the segment weight


def f_lattice(d, k):
    """f(k) = Σ_{j=0..d} e^{i k_j} in the lattice basis, k_0 ≡ 0 (z = d+1 bonds, A_d/T19)."""
    return 1 + sum(sp.exp(sp.I * k[a]) for a in range(d))


def node_k(d):
    """the EXACT node ∀d: k_a = 2πa/(d+1) ⟹ Σ_{j=0}^{d} e^{2πij/(d+1)} = 0 (the sum of all roots)."""
    return [2 * sp.pi * sp.Rational(a + 1, d + 1) for a in range(d)]


def chain_circulant(N, t=1):
    """★the S1057 machinery, verbatim: a nearest-neighbor chain, closed."""
    S = sp.zeros(N, N)
    for i in range(N):
        S[i, (i + 1) % N] = 1
    return 2 * t * sp.eye(N) - t * (S + S.T)


def hop_matrix(N, reach=1):
    """The bare hop matrix on a closed chain: edges i ↔ i+reach."""
    A = sp.zeros(N, N)
    for i in range(N):
        A[i, (i + reach) % N] += 1
        A[(i + reach) % N, i] += 1
    return A


def min_power_connecting(A, u, v, nmax=12):
    """The minimal power n≥1 with (Aⁿ)[u,v] ≠ 0.  None ⟺ not connected up to nmax."""
    P = sp.eye(A.rows)
    for n in range(1, nmax + 1):
        P = P * A
        if sp.simplify(P[u, v]) != 0:
            return n
    return None


def supercell_bloch(P, phi, t=1):
    """★the S1057 machinery, verbatim (including the ACCUMULATE-fix): the same chain,
    described with P nodes per cell, a Bloch phase phi across the cell."""
    H = sp.zeros(P, P)
    for i in range(P):
        H[i, i] = 2 * t
    for i in range(P - 1):
        H[i, i + 1] += -t
        H[i + 1, i] += -t
    H[0, P - 1] += -t * sp.exp(-sp.I * phi)
    H[P - 1, 0] += -t * sp.exp(sp.I * phi)
    if P == 1:
        H[0, 0] = 2 * t - t * (sp.exp(sp.I * phi) + sp.exp(-sp.I * phi))
    return H


# ─────────────────────────── T1: SCALE (justification E4) ───────────────────────────
def T1_scale():
    print("T1. SCALE — does a REAL parameter enter the step-count (justification E4)")
    c = sp.Symbol('c', positive=True)
    G0_c = (c * E * I2 - c * m0 * sz).inv()
    step_plain = sp.simplify(H2c * G0)
    step_scaled = sp.simplify((c * H2c) * G0_c)
    inv = sp.simplify(step_scaled - step_plain) == sp.zeros(2, 2)
    ok(inv, "a COMMON scale: (H·G₀) is identical under (H,E,m₀)→c·(H,E,m₀) — ★this is an IDENTITY BY "
            "CONSTRUCTION (weight 0), NOT a measurement  [T1 · dimensionless · a 2×2 identity]")
    print("     ★BETA'S VISA (the S1059-visa): this line by itself is NOT A WITNESS — H·G₀ has weight 0,")
    print("       so invariance cannot fail to occur.  The substantive form is below, and it")
    print("       gives a MORE PRECISE support for the main thesis than a cancelling of scale.")
    # ★THE SUBSTANTIVE FORM (Beta's route): scale H INDEPENDENTLY of E.
    t = sp.Symbol('t', positive=True)
    term_n = []
    cur = G0
    for n in range(4):
        cur_t = sp.simplify(G0 * (t * H2c * G0) ** n)
        cur_1 = sp.simplify(G0 * (H2c * G0) ** n)
        ratio = sp.simplify(cur_t[0, 1] / cur_1[0, 1]) if cur_1[0, 1] != 0 else                 sp.simplify(cur_t[0, 0] / cur_1[0, 0])
        term_n.append(sp.simplify(ratio - t ** n) == 0)
    grades_split = all(term_n)
    ok(grades_split,
       "★SUBSTANTIVELY: under H→t·H (E,m₀ left UNTOUCHED) the term of index n goes EXACTLY as t^n, n=0..3 ⟹ "
       "the parameter SEPARATES the gradations, it does not mix them  [T1 · dimensionless · a ratio of series terms]")
    print("     ⟹ n is an integer NOT because a scale cancelled, but because the parameter enters as")
    print("       an EXPONENT.  A continuum handle cannot substitute for the integer, because it sits elsewhere.")
    # a negative control: under a COMMON scale the separation vanishes — the detector must see this.
    # ★the first form of this line COMPARED THE TERMS and failed: the whole resolvent carries a common weight
    #   1/c, so the term n scales as 1/c regardless of n.  The separation of gradations lives not in
    #   the terms, but in the RATIO member_n/member_0 — that is what is measured.  The FORMULATION of the test
    #   was fixed, not the claim: a common weight is a known fact about G, not news about the gradation.
    def ratio_n0(Hop, Gz, n):
        top = sp.simplify(Gz * (Hop * Gz) ** n)
        bot = sp.simplify(Gz)
        return sp.simplify((top[0, 1] if top[0, 1] != 0 else top[0, 0]) /
                           (bot[0, 0] if bot[0, 0] != 0 else bot[0, 1]))
    r_indep = sp.simplify(ratio_n0(t * H2c, G0, 2) / ratio_n0(H2c, G0, 2))
    r_common = sp.simplify(ratio_n0(c * H2c, G0_c, 2) / ratio_n0(H2c, G0, 2))
    common_blind = (sp.simplify(r_indep - t ** 2) == 0) and (sp.simplify(r_common - 1) == 0)
    ok(common_blind,
       "★N3 A NEGATIVE CONTROL for T1: the ratio member₂/member₀ under an INDEPENDENT scale = t², under "
       "a COMMON one = 1 ⟹ the detector DISTINGUISHES the two scalings, it does not always say «invariant»")
    return inv and grades_split and common_blind


# ─────────────────────── T2: FROM BELOW, does the step divide (E1) ───────────────────────
def T2_below():
    print("T2. FROM BELOW — does a HALF-STEP exist: an operator K with K² = H (justification E1)")

    # (0) A LIVENESS CONTROL: without constraints a root EXISTS ⟹ it is precisely {Hermiticity, locality}
    #     that kill it, not «no root at all».  Without this line the test would be vacuous.
    lam = sp.sqrt(sp.Symbol('s', positive=True))
    Hnum = sp.Matrix([[0, 1], [1, 0]])                       # f=1 (fr=1, fi=0)
    Kbare = sp.Matrix([[sp.Rational(1, 2) * (1 + sp.I), sp.Rational(1, 2) * (1 - sp.I)],
                       [sp.Rational(1, 2) * (1 - sp.I), sp.Rational(1, 2) * (1 + sp.I)]])
    bare_ok = sp.simplify(Kbare * Kbare - Hnum) == sp.zeros(2, 2)
    ok(bare_ok, "★A LIVENESS CONTROL: a bare (non-Hermitian) root K²=H EXISTS at f=1 — "
                "the test is not vacuous  [T2.0 · dimensionless · a matrix identity]")
    ok(sp.simplify(Kbare - Kbare.conjugate().T) != sp.zeros(2, 2),
       "   and it is NOT Hermitian (K ≠ K†) — that is, it is precisely the Hermiticity requirement that kills it, not existence")

    # (i) MECHANISM-1: Hermiticity.  K†=K ⟹ K² is positive-semidefinite ⟹ spec(H) ⊆ [0,∞).
    ev = list(H2c.eigenvals().keys())
    ev_s = sorted([sp.simplify(e) for e in ev], key=lambda e: sp.N(e.subs({fr: 1, fi: 0})))
    neg_exists = sp.simplify(ev_s[0] + sp.sqrt(absf2)) == 0
    ok(neg_exists, "★MECHANISM-1 (Hermiticity): spec(H) = ±√(f_re²+f_im²) — a NEGATIVE eigen"
                   "value EXISTS at f≠0  [T2.i · dimensionless · eigenvalues of a 2×2]")
    print("     ⟹ there is no Hermitian K with K²=H: K†=K would give K² ⪰ 0, but H has −|f| < 0.")

    # (ii) MECHANISM-2: locality.  A finite radius ⟺ the entries of K are trigonometric
    #      polynomials in k ⟹ det K is analytic.  But (det K)² = det H = −|f|² ⟹ det K = ±i|f|,
    #      and |f| is NOT differentiable at the node f=0 (it exists ∀d — S1045, dim=d−2).
    detH = sp.simplify(H2c.det())
    ok(sp.simplify(detH + absf2) == 0,
       "★det H = −(f_re²+f_im²)  [T2.ii · dimensionless · a 2×2 determinant]")
    # ★REVISION-2 PER BETA'S VISA.  Revision-1 measured |f| along a ray with the expression √((as)²+(bs)²) — this
    #   is a MODEL expression that ASSUMES linearity of the node; at the degenerate point (S1047:
    #   (0,π/2,π/2), ∇Re f=0, |f|~s²) the mechanism falls silent.  The replacement removes the ansatz ENTIRELY:
    #   K local ⟹ the entries are Laurent polynomials ⟹ so is det K.  (det K)² = det H = −f·f̄,
    #   so −f·f̄ must be a SQUARE in C[X^{±1}].  A UFD ring ⟹ FACTORIZATION decides,
    #   with no analysis at all and — ★crucially — WITHOUT HERMITICITY.
    def laurent_detH(d):
        """−f·f̄ as a polynomial (up to a unit ∏X): f = 1+ΣX_a, f̄·∏X = ∏X + Σ_a ∏_{b≠a}X_b."""
        X = sp.symbols(f'X1:{d+1}', positive=True)
        f_pol = 1 + sum(X)
        g_pol = sp.prod(X) + sum(sp.prod([X[b] for b in range(d) if b != a]) for a in range(d))
        return sp.expand(-f_pol * g_pol), f_pol, g_pol, X

    all_notsquare = []
    for dd in (2, 3, 4):
        prod_pol, f_pol, g_pol, X = laurent_detH(dd)
        fl = sp.factor_list(prod_pol)
        odd_factors = [(fac, mult) for fac, mult in fl[1] if mult % 2 == 1]
        # a square in a UFD ⟺ EVERY irreducible factor has even multiplicity (the unit −1 = i² does not interfere)
        is_square = len(odd_factors) == 0
        all_notsquare.append(not is_square and len(odd_factors) >= 2)
        print(f"     d={dd}: −f·f̄ = a unit · " + " · ".join(
            f"({sp.factor(fac)})^{mult}" for fac, mult in fl[1])
            + f"  ⟹ factors of odd multiplicity: {len(odd_factors)}")
    ok(all(all_notsquare),
       "★★T2.ii REVISION-2 (the visa route): −f·f̄ has TWO DISTINCT irreducible factors of ODD multiplicity "
       "∀d∈{2,3,4} ⟹ not a square in C[X^{±1}] ⟹ NO LOCAL root EXISTS  "
       "[T2.ii · dimensionless · factorization in a UFD]")
    # a negative control: on H² the multiplicities are even ⟹ the detector must say YES
    prod_pol, f_pol, g_pol, X = laurent_detH(2)
    sq_pol = sp.expand(prod_pol ** 2)
    fl_sq = sp.factor_list(sq_pol)
    ok(all(m % 2 == 0 for _, m in fl_sq[1]),
       "★N4 A NEGATIVE CONTROL for T2.ii: on the square (det H)² all multiplicities are EVEN ⟹ the detector says "
       "YES and is able to fail; this agrees with T2.0 (at f=1 we have −1=i², a square — it does not block)")
    print("     ★THE CLAIM IS RAISED (per the visa's pointer): this argument NOWHERE uses Hermiticity ⟹")
    print("       the correct statement — «no local root exists AT ALL», not «within the class")
    print("       {Hermitian, local}».  Revision-1 under-claimed its own result; fixed.")
    print("     ★AND the mechanism rests on the LATTICE ITSELF (the factorization of f), not the analysis of the node.")

    # the node exists ∀d — checked EXACTLY, not cited
    node_ok = []
    for d in (2, 3, 4, 5):
        k = node_k(d)
        val = sp.simplify(sp.expand(sp.re(f_lattice(d, k))) + sp.I * sp.simplify(sp.expand(sp.im(f_lattice(d, k)))))
        node_ok.append(sp.simplify(val) == 0)
    ok(all(node_ok), "★the node f=0 EXISTS ∀d∈{2,3,4,5} at the point k_a=2πa/(d+1) (the sum of all "
                     "(d+1)-th roots of 1)  [T2.ii · dimensionless · an exact symbolic substitution]")

    # (iii) NEGATIVE CONTROL N2: the «Hermitian root» detector MUST say YES on a positive-definite matrix
    Ppos = sp.Matrix([[2, 1], [1, 2]])
    evp = sorted([sp.simplify(e) for e in Ppos.eigenvals().keys()], key=sp.N)
    Kpos_exists = all(sp.N(e) > 0 for e in evp)
    ok(Kpos_exists, "★N2 A NEGATIVE CONTROL: on the positive-definite [[2,1],[1,2]] (spec={1,3}) the detector "
                    "says YES — a Hermitian root DOES exist ⟹ the detector is able to fail")

    return neg_exists and all(node_ok) and Kpos_exists and bare_ok


# ─────────────────────── T3: FROM ABOVE, is the step not 2+ (E2) ───────────────────────
def T3_above():
    print("T3. FROM ABOVE — are there OMISSIONS in the series (justification E2)")
    term = G0
    empties, shown = [], []
    for n in range(0, 7):
        nonzero = sp.simplify(term.norm()) != 0
        blk = 'AB≠0' if sp.simplify(term[0, 1]) != 0 else 'AB=0'
        dia = 'AA≠0' if sp.simplify(term[0, 0]) != 0 else 'AA=0'
        shown.append(f"n={n}: {blk}, {dia}")
        if not nonzero:
            empties.append(n)
        term = sp.simplify(term * H2c * G0)
    print("     " + " | ".join(shown))
    print("     ★BETA'S VISA: the non-vanishing of the series terms CANNOT BE A WITNESS for E2 — it")
    print("       is structural: n counts the powers of THAT operator, whatever its reach, and")
    print("       an empty term would require the product of nonzero 2×2 matrices to be zero.")
    print("       Running it on an operator of reach 2 gives the same result — the detector is blind.")
    print("       ⟹ the line below remains as an OBSERVATION (an alternation, S1045), not as a witness.")
    print("     (the AB↔AA alternation — is binary; this is the PARITY of the block, not an omitted term)")

    # NEGATIVE CONTROL N1: an object where the unit IS GENUINELY 2 base edges — the detector must catch this
    N = 8
    A1 = hop_matrix(N, reach=1)
    A2 = hop_matrix(N, reach=2)
    u1 = min_power_connecting(A1, 0, 1)
    u2_neighbour = min_power_connecting(A2, 0, 1)
    u2_second = min_power_connecting(A2, 0, 2)
    print(f"     the r-space detector: reach=1 → min power(0→1) = {u1} · "
          f"reach=2 → (0→1) = {u2_neighbour}, (0→2) = {u2_second}")
    ok(u1 == 1 and u2_neighbour is None and u2_second == 1,
       "★★THE E2 WITNESS (relocated per the visa): reachability in r-space — at reach=1 the neighbor "
       "is reached by EXACTLY one power; at reach=2 the neighbor is unreachable by ANY power, and the unit = 2 "
       "base edges  [T3 · a base edge · the minimal power with a nonzero entry]")
    print("     ⟹ it is PRECISELY THIS detector that is able to fail and it is the one that carries E2; the non-vanishing test")
    print("       of the series terms no longer carries a witness (the empty safeguard is removed, not hidden).")
    return u1 == 1 and u2_neighbour is None and u2_second == 1


# ─────────────────── T4: RELABELING — the same S1057 machinery (E3) ───────────────────
def T4_relabel():
    print("T4. RELABELING — the same supercell machinery that SHIFTED the length (justification E3)")
    N = 12
    ref = chain_circulant(N)
    ref_spec = sorted([sp.nsimplify(sp.simplify(v)) for v, m in ref.eigenvals().items()
                       for _ in range(m)], key=lambda e: sp.N(e))
    same_all = True
    for Pv in (2, 3, 4, 6):
        cells = N // Pv
        got = []
        for j in range(cells):
            phi_j = 2 * sp.pi * sp.Rational(j, cells)
            Hb = supercell_bloch(Pv, phi_j)
            for v, m in Hb.eigenvals().items():
                got.extend([sp.nsimplify(sp.simplify(v))] * m)
        got = sorted(got, key=lambda e: sp.N(e))
        same = len(got) == len(ref_spec) and all(
            abs(sp.N(x - y)) < sp.Float("1e-25") for x, y in zip(got, ref_spec))
        same_all &= bool(same)
        ok(same, f"P={Pv}: the supercell spectrum is IDENTICAL to the direct chain (N={N}) ⟹ THE SAME object")

    # ★REVISION-2 PER BETA'S VISA — and this was the most expensive finding.  Revision-1 took A = hop_matrix(N,1)
    #   OUTSIDE the loop and ran `for _ in (2,3,4,6)`: four iterations counted the same quantity,
    #   and the assert checked that four identical computations gave the same thing.  That is, what
    #   happened to the LENGTH under relabeling was never reproduced in the test at all — and it is
    #   precisely on this that the thesis «it happened to the length, it did not happen to the count» rests.  Now P enters
    #   the CONSTRUCTION of the matrix, and the SAME two physical nodes are compared in both descriptions.
    def hop_supercell(Ntot, P):
        """The same chain, BUILT block-by-block: the node (c,s) → the index c·P+s."""
        cells = Ntot // P
        A = sp.zeros(Ntot, Ntot)
        for c in range(cells):
            for sidx in range(P - 1):                     # edges INSIDE the cell
                i, j = c * P + sidx, c * P + sidx + 1
                A[i, j] += 1; A[j, i] += 1
            i, j = c * P + (P - 1), ((c + 1) % cells) * P  # the edge BETWEEN cells
            A[i, j] += 1; A[j, i] += 1
        return A

    per_bond, per_cell = [], []
    for Pv in (2, 3, 4, 6):
        Asup = hop_supercell(N, Pv)
        per_bond.append(min_power_connecting(Asup, 0, 1))          # two PHYSICAL neighbors
        per_cell.append(min_power_connecting(Asup, 0, Pv))         # the same node of the next cell
    print("     P | hops per BOND | hops per CELL")
    for Pv, hb, hc in zip((2, 3, 4, 6), per_bond, per_cell):
        print(f"     {Pv} | {hb} | {hc}")
    inv_bond = len(set(per_bond)) == 1 and per_bond[0] == 1
    moves_cell = per_cell == [2, 3, 4, 6]
    ok(inv_bond and moves_cell,
       "★the count «per BOND» = 1 and does NOT move under relabeling · the count «per CELL» moves "
       "as P — ★on a matrix BUILT from P (revision-2 per the visa)  "
       "[T4 · a hop · the minimal power in a supercell construction]")
    print("     ⟹ the same free choice of description that shifted the LENGTH as 1/P² (S1057), does NOT")
    print("       move the hop-count.  The asymmetry is measured on ONE object by the same machinery.")
    print("     ★A LINK WITH T2: to take the CELL as the unit, the series would need a power of 1/P")
    print("       (the neighbor 0→1 = «half a cell» at P=2) — but T2 proved that a LOCAL fractional")
    print("       power does not exist at all.  The cell-unit is not inconvenient, it is IMPOSSIBLE.")
    return same_all and inv_bond and moves_cell


# ─────── T5: THE MOST OBVIOUS OBJECTION — just SPLIT THE BOND IN HALF (E1-bis) ───────
def T5_subdivision():
    """The step can be «subdivided» by adding a node in the middle of each bond (a subdivision).
    The question is not «can it be drawn», but: is this the SAME object in a new description (like
    a supercell) — or a DIFFERENT object?  The criterion is the same as in T4: the spectrum.

    ★REVISION-2 PER BETA'S VISA, three corrections, and all THREE STRENGTHEN the verdict:
      (1) revision-1 measured on C6 — this is r=2 and m=n, the ONLY regular case where subdivision exactly
          DOUBLES the states.  Our object has z=d+1, that is r=3.  The word «doubles» is dropped: at r=3
          the excess = m−n, that is, subdivision costs MORE than revision-1 showed.
      (2) the form all(any(...)) saw neither the multiplicities nor the m−n extra zeros — replaced by
          the STRONG form: equality of MULTISETS with multiplicities.
      (3) 📖 THE PILLAR IS NAMED: μ²=λ+r for the subdivision of an r-regular graph together with a multiplicity
          of m−n zeros — a textbook result of spectral graph theory
          (Cvetković–Doob–Sachs).  It is CHECKED here with my own hand, but the status grows
          from the pillar, not from the check — that is why the 📖 mark sits on that line, not in a footnote.
    """
    print("T5. THE MOST OBVIOUS OBJECTION — «just add a node in the middle of a bond» (E1-bis)")
    print("     📖 the pillar: μ²=λ+r + (m−n) zeros = spectral graph theory (Cvetković–Doob–Sachs);")
    print("        checked with my own hand, but the status is carried by the PILLAR, not my check.")

    def subdivision(edges, n):
        """The subdivision graph: n original vertices + one vertex per edge."""
        m = len(edges)
        A = sp.zeros(n + m, n + m)
        for idx, (a, b) in enumerate(edges):
            e = n + idx
            A[a, e] += 1; A[e, a] += 1
            A[b, e] += 1; A[e, b] += 1
        return A

    def edges_of(A, n):
        return [(i, j) for i in range(n) for j in range(i + 1, n) if A[i, j] != 0]

    def spec(A):
        return sorted([sp.nsimplify(sp.simplify(v)) for v, mult in A.eigenvals().items()
                       for _ in range(mult)], key=lambda e: sp.N(e))

    # three objects: C6 (r=2, m=n — the revision-1 case) · K33 (r=3) · a honeycomb-torus 2×2 (r=3, our z=d+1)
    C6 = hop_matrix(6, reach=1)
    K33 = sp.zeros(6, 6)
    for a in range(3):
        for b in range(3, 6):
            K33[a, b] = 1; K33[b, a] = 1
    HC = sp.zeros(8, 8)                       # a honeycomb on a 2×2 torus: 8 vertices, 12 edges, r=3
    hc_edges = [(0, 4), (0, 5), (0, 6), (1, 5), (1, 6), (1, 7),
                (2, 4), (2, 6), (2, 7), (3, 4), (3, 5), (3, 7)]
    for a, b in hc_edges:
        HC[a, b] = 1; HC[b, a] = 1

    results = []
    for name, A, r in (("C6 (r=2, m=n — the revision-1 case)", C6, 2),
                       ("K33 (r=3)", K33, 3),
                       ("a honeycomb-torus 2×2 (r=3 = our z=d+1)", HC, 3)):
        n = A.rows
        ed = edges_of(A, n)
        m = len(ed)
        S = subdivision(ed, n)
        sp_small, sp_big = spec(A), spec(S)
        # the STRONG form: the multiset spec(S) = {±√(λ+r)} ∪ {0 × (m−n)}
        predicted = []
        for lam in sp_small:
            val = sp.simplify(lam + r)
            root = sp.sqrt(val)
            predicted += [sp.nsimplify(root), sp.nsimplify(-root)]
        predicted += [sp.Integer(0)] * (m - n)
        predicted = sorted(predicted, key=lambda e: sp.N(e))
        same = len(predicted) == len(sp_big) and all(
            abs(sp.N(x - y)) < sp.Float("1e-20") for x, y in zip(predicted, sp_big))
        results.append((name, n, m, len(sp_big), m - n, same))
        print(f"     {name}: n={n}, m={m} ⟹ the subdivision has {len(sp_big)} states "
              f"(the excess m−n = {m - n})  [T5 · a state · the count of eigenvalues]")
        ok(same, f"★the STRONG form of the law on «{name}»: the multiset spec(subdivision) = "
                 f"{{±√(λ+r)}} ∪ {{0×(m−n)}} EXACTLY, with multiplicities  "
                 f"[T5 · dimensionless · equality of multisets]")

    doubling_only_c6 = (results[0][4] == 0) and all(r[4] > 0 for r in results[1:])
    ok(doubling_only_c6,
       "★★«DOUBLING» — A PROPERTY OF C6, NOT OF OUR LATTICE (a visa finding): at m=n the excess is 0, "
       "while at r=3 the excess m−n>0 ⟹ subdivision costs EVEN MORE than revision-1 showed")
    print("     ⟹ ★THE KEY (stands unchanged): subdivision gives a root NOT of H, but of (H + r·I) — of a")
    print("       SHIFTED operator.  The shift by r is exactly what removes the obstruction of T2.i.")
    print("       Subdivision BUYS the half-step at exactly the price T2 named: by changing the object, not the description.")
    print("     ⟹ splitting a bond = a DIFFERENT LATTICE (more states, a different spectrum), NOT a relabeling")
    print("       of ours.  The step is indivisible RELATIVE TO the lattice; the lattice is set by T19, not this probe.")
    return all(r[5] for r in results) and doubling_only_c6


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    logf = open(os.path.join(_HERE, "S1059_1_run.log"), "w", encoding="utf-8")
    tee = Tee(sys.stdout, logf)
    sys.stdout = tee
    sys.path.insert(0, os.path.join(_HERE, ".."))

    print("=" * 78)
    print("S1059 — A KILL-TEST OF ANSATZ A: is the chain step INDIVISIBLE (step 3 of the assignment)")
    print("Layer-1: the GRADATION of the series G = Σ_n G₀(H·G₀)ⁿ is on trial, not its Layer-2-reading.")
    print("=" * 78)
    print()

    t1 = T1_scale(); print()
    t2 = T2_below(); print()
    t3 = T3_above(); print()
    t4 = T4_relabel(); print()
    t5 = T5_subdivision(); print()

    # ───────────── JUSTIFICATION E1-E4 (a mandatory field of assignment №1) ─────────────
    print("JUSTIFICATION E1-E4 — what WOULD FORCE ansatz A to FAIL (all four computed, not declared)")
    print(f"  E1 «a LOCAL K, K²=H exists» ⟹ the step divides, A is demoted: "
          f"{'DID NOT OCCUR' if t2 else 'OCCURRED'}")
    print("     ★THE CLAIM IS RAISED (revision-2 per the visa): the «Hermitian» condition is DROPPED from the E1 formulation —")
    print("     the argument via the factorization of det kills ANY local root, Hermitian or not.")
    print("     (two mechanisms: spec(H) ∌ [0,∞) — requires Hermiticity ⊥ −f·f̄ not a square in")
    print("      C[X^±1] — does NOT require it; the second is stronger, and it is the one that carries E1)")
    print(f"  E2 «the neighbor is reached by more than one hop» ⟹ the unit > 1, A is wrong: "
          f"{'DID NOT OCCUR' if t3 else 'OCCURRED'}")
    print("     ★THE WITNESS WAS RELOCATED (revision-2 per the visa) from the non-vanishing of the series terms (an empty")
    print("     detector — is silent even where the unit is genuinely 2) to REACHABILITY in r-space.")
    print(f"  E3 «the count moves under a free relabeling» ⟹ a convention, as the length once was: "
          f"{'DID NOT OCCUR' if t4 else 'OCCURRED'}")
    print(f"  E4 «a real parameter enters the count» ⟹ a handle was added, a STOP per K2: "
          f"{'DID NOT OCCUR' if t1 else 'OCCURRED'}")
    print()

    # ───────────── AN HONEST MULTIPLICITY COUNT (the rule: witnesses by ANCESTOR) ─────────────
    print("★THE MULTIPLICITY COUNT, HONESTLY (the rule «witnesses by ANCESTOR, not by citation»)")
    print("  The mechanisms T2.i (Hermiticity) and T2.ii (locality+node) both read ONE")
    print("  structure — the off-diagonal form H=[[0,f],[f̄,0]] (T19-bipartite + T33).  A shared")
    print("  ancestor ⟹ MULTIPLICITY 1, not 2.  The exante named this risk BEFORE the count; it came true.")
    print("  T4 rests on a different ancestor (the S1057-machinery of relabeling) — but it measures")
    print("  the INVARIANCE of the count, not its indivisibility: this is the other half of the question, not a second")
    print("  witness of the same one.  ⟹ on indivisibility the witness is SINGLE.")
    print()

    # ───────────── FENCE ─────────────
    from tools.fence_scan import scan_forbidden
    pats = [r"світл\w*", r"спін\w*", r"\bspin\b", r"квант\w*", r"\bquantum\b",          # GUARDLINE
            r"\blight\b", r"\bдія\b", r"\bдії\b", r"\bдією\b", r"\baction\b",           # GUARDLINE
            r"Фейнман\w*", r"\bFeynman\b", r"\bшлях\w*", r"пропагатор\w*",              # GUARDLINE
            r"чекерборд\w*"]                                                            # GUARDLINE
    hits = scan_forbidden(__file__, pats) + scan_forbidden("".join(tee.chunks), pats)
    ok(not hits, f"the Layer-1 fence: 0 hits on {len(pats)} patterns (source + log)")
    print("     (English `path` is deliberately NOT scanned — a technical homonym of os.path, not a Layer-2 word)")

    print()
    print("=" * 78)
    verdict = t1 and t2 and t3 and t4 and t5
    print("★A BINARY OUTPUT (as the assignment asked):")
    if verdict:
        print("  THE CHAIN STEP IS FORCED to one hop — FROM BELOW (a LOCAL fractional power")
        print("  does not exist at all, without any Hermiticity condition) and FROM ABOVE (the neighbor is reached by exactly one")
        print("  hop, and on a lattice of reach 2 — by none), and moreover")
        print("  the hop-count is INVARIANT under the same relabeling that shifted the length.")
        print("  ⟹ none of the four justifications occurred; NO HANDLE WAS ADDED:")
        print("  the unit of gradation — is an INTEGER, not a value, so the question «what value»")
        print("  does not apply to it.  Wager A stands as an OUTPUT at Layer-1.")
    else:
        print("  AT LEAST ONE JUSTIFICATION OCCURRED — wager A is demoted, named.")
    print(f"SCORE: {len(_passes)} ✓ / {len(_fails)} ✗")
    print("STOP. The registry is untouched; step 4 is not begun.")
    print("=" * 78)
    sys.stdout = tee.real
    logf.close()
    return 1 if _fails else 0


if __name__ == "__main__":
    sys.exit(main())
