#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (linear algebra over Q; handles 0). W28-O6 — Pfaffian mirror
#      of the log-pair (η,Ω): det Ω = Pf(Ω)² (even n symbolic + random integers; odd n
#      det ≡ 0 symbolic) · invariance of det A under congruence (g,Ω) → (SᵀgS,SᵀΩS) ·
#      table of det A across the signature ladder {(2,1)·(2,2)·(3,1)·(4,0)·(3,3)} ·
#      orientation flip (swap of two axes, Pf and det A by the numbers) ·
#      intersection with the A² = ±𝟙 examples of S903.
#      ★BLINDNESS: the probe prints ONLY raw identities/numbers/histograms; reading = an act of the court.
"""
S905 (lane A, ed.2) — W28-O6: Pfaffian mirror of the log-pair (η,Ω).

  Objects: signatures {(2,1)·(2,2)·(3,1)·(4,0)·(3,3)}, η = diag(±1);
  antisymmetric Ω (integer); A = η⁻¹Ω (η⁻¹ = η, diagonal ±1).

  O6a — Pfaffian root: Pf(Ω) by recursion on the first row
        (Pf = Σ_{j≥2} (−1)^j ω_{1j} Pf(Ω without rows/columns 1,j); Pf(2×2)=ω₁₂;
        Pf(0×0)=1). Verification det Ω − Pf(Ω)² = 0: symbolic (full antisymmetric
        Ω on symbols) for n = 4 and n = 6, expanded to zero; + 50 random integer
        Ω per n (entries {−9..9}, seed 905). For n = 3 and n = 5 — symbolic
        det Ω ≡ 0 (expanded to zero).
  O6b — invariance: for n = 4 (signatures (2,2)·(3,1)·(4,0)) and n = 6 ((3,3)):
        20 random integer S (entries {−3..3}, det S ≠ 0, seed per signature):
        det A′ where g′ = SᵀgS, Ω′ = SᵀΩS, A′ = g′⁻¹Ω′ — bit-match with det A (raw
        line); separately det g·det Ω multiplied by (det S)⁴ — factors by the numbers for
        3 examples.
  O6c — ladder table: for every signature det A for EACH single
        coordinate generator (J and K separately) · for the sums {all J} and {all K} ·
        for 100 random integer A (coeffs {−3..3}, seed 905): histogram of the
        sign of det A {positive · zero · negative} and min/max; det g alongside.
        For odd n ((2,1)) — det A ≡ 0 symbolic (separate line).
  O6d — orientation flip: n = 4 ((2,2)) and n = 6 ((3,3)): swap of two
        axes P — a same-sign pair and a mixed-sign pair separately: Pf(PᵀΩP) against
        Pf(Ω) (sign by the numbers, 5 examples of each type) and det A′ = det(η·PᵀΩP)
        against det A (by the numbers alongside). Raw lines with no readings.
  O6e — intersection with the square structures of S903 (address: S903_run.log, block C):
        A² = −𝟙 on (2,2), coeffs [−1,−1,0,0,0,0] in the S903 generator basis;
        A² = +𝟙 on (2,2), first example of S903, coeffs [−2,−2,−2,−1,−1,2];
        A² = +𝟙 on (3,3), first example of S903: K(0,3)+K(1,4)+K(2,5).
        det A by the numbers; the formula line det(A²) = det(±𝟙) = (±1)ⁿ (sympy) alongside.

  EXACT arithmetic: sympy Rational/Integer over Q; no tolerances. Mechanisms
  (make_soPQ/flat) — a verbatim copy from S903 (reproducibility). No silent
  truncation: every bound is stamped in the printed lines.

Fence: the shared fence_scan helper (forbidden words — on the GUARDLINE line;
"orientation" allowed as a math term). No verdicts in the probe's text.
"""
import sys
import random

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# shared fence_scan helper (src/tools)
_src = __file__.replace("\\", "/").rsplit("/src/", 1)[0] + "/src"
if _src not in sys.path:
    sys.path.insert(0, _src)
from tools.fence_scan import scan_forbidden   # noqa: E402

import sympy as sp   # noqa: E402

# ── tee: all of stdout is duplicated into S905_run.log next to the script ──
_LOG_PATH = __file__.replace("\\", "/").rsplit("/", 1)[0] + "/S905_run.log"


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)
            st.flush()

    def flush(self):
        for st in self.streams:
            st.flush()


_logf = open(_LOG_PATH, "w", encoding="utf-8")
sys.stdout = _Tee(sys.stdout, _logf)

FAIL = []
N_CHECKS = 0


def check(name, cond, detail=""):
    global N_CHECKS
    N_CHECKS += 1
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" — {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)
    return ok


def rule(t):
    print("\n" + "=" * 96)
    print(t)
    print("=" * 96)


# ── SECTION 0 — FENCE VIA SHARED HELPER (first) ──────────────────────────────────────
rule("SECTION 0 — FENCE (shared fence_scan helper)")
_FORBIDDEN = [r"стійк", r"обрано", r"selected", r"stable", r"причин", r"час", r"time", r"arrow", r"стріла", r"вибір", r"краще", r"злам", r"поле", r"матерія", r"хіральн"]   # GUARDLINE
_hits = scan_forbidden(__file__, _FORBIDDEN)
check("fence: forbidden words (list = GUARDLINE line) = 0 occurrences outside the declaration",
      not _hits, f"hits: {_hits}" if _hits else "0")
check("handles 0 (pure algebra: Pfaffian · det · congruence · signature ladder)", True)


# ═══════════════════════════════════════════════════════════════════════════════════════
# TOOLS (make_soPQ/flat — verbatim copy from S903; Pfaffian — recursion on the first row)
# ═══════════════════════════════════════════════════════════════════════════════════════

def flat(M):
    """n×n matrix → flat n² tuple (row-major), sympy-exact entries."""
    n = M.rows
    return tuple(M[r, c] for r in range(n) for c in range(n))


def make_soPQ(p, q):
    """η, J-type and K-type generators with plane labels (i,j)."""
    n = p + q
    signs = [1] * p + [-1] * q
    eta = sp.diag(*signs)
    Js, Ks = [], []
    for i in range(n):
        for j in range(i + 1, n):
            E_ij = sp.zeros(n, n)
            E_ij[i, j] = 1
            E_ji = sp.zeros(n, n)
            E_ji[j, i] = 1
            if signs[i] == signs[j]:
                Js.append((f"J({i},{j})", E_ij - E_ji))
            else:
                Ks.append((f"K({i},{j})", E_ij + E_ji))
    return eta, Js, Ks


def pfaffian(M):
    """Pf(M) by recursion on the first row (M antisymmetric, even size):
    Pf = Σ_{j=1..n−1} (−1)^{j+1}·M[0,j]·Pf(M without rows/columns 0,j);
    Pf(0×0) = 1; Pf(2×2) = M[0,1]. (0-based indices; the sign is the classic 1-based convention.)"""
    n = M.rows
    if n == 0:
        return sp.Integer(1)
    if n == 2:
        return M[0, 1]
    tot = sp.Integer(0)
    for j in range(1, n):
        a = M[0, j]
        if a == 0:
            continue
        idx = [k for k in range(n) if k != 0 and k != j]
        tot = tot + sp.Integer(-1) ** (j + 1) * a * pfaffian(M.extract(idx, idx))
    return tot


def sym_omega(nn, tag):
    """Full antisymmetric Ω on symbols w{tag}_{i}_{j} (upper triangle)."""
    O = sp.zeros(nn, nn)
    for i in range(nn):
        for j in range(i + 1, nn):
            w = sp.Symbol(f"w{tag}_{i}_{j}")
            O[i, j] = w
            O[j, i] = -w
    return O


def rand_omega(nn, rng, lo, hi):
    """Random integer antisymmetric Ω, upper-triangle entries ∈ {lo..hi}."""
    O = sp.zeros(nn, nn)
    for i in range(nn):
        for j in range(i + 1, nn):
            v = rng.randint(lo, hi)
            O[i, j] = v
            O[j, i] = -v
    return O


def rand_S(nn, rng, lo, hi):
    """Random integer S with det S ≠ 0 (retry until success; attempt count is returned)."""
    tries = 0
    while True:
        tries += 1
        S = sp.Matrix(nn, nn, lambda i, j: rng.randint(lo, hi))
        if S.det() != 0:
            return S, tries


def sign_word(v):
    if v > 0:
        return "positive"
    if v < 0:
        return "negative"
    return "zero"


SEARCH_STAMPS = []

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O6a — PFAFFIAN ROOT: det Ω = Pf(Ω)² (even n) · det Ω ≡ 0 (odd n)")

# base recursion identities (stamp: explicit reference values for small n)
_O2 = sym_omega(2, "b2")
check("O6a: Pf(2×2) = ω₀₁ and det(2×2) = ω₀₁² (symbolic)",
      pfaffian(_O2) == _O2[0, 1] and sp.expand(_O2.det() - _O2[0, 1] ** 2) == 0)
_O4t = sym_omega(4, "b4")
_pf4_ref = _O4t[0, 1] * _O4t[2, 3] - _O4t[0, 2] * _O4t[1, 3] + _O4t[0, 3] * _O4t[1, 2]
check("O6a: Pf(4×4) = ω₀₁ω₂₃ − ω₀₂ω₁₃ + ω₀₃ω₁₂ (recursion = classical form, symbolic)",
      sp.expand(pfaffian(_O4t) - _pf4_ref) == 0)

for nn in (4, 6):
    Osym = sym_omega(nn, f"a{nn}")
    n_syms = nn * (nn - 1) // 2
    pf_s = pfaffian(Osym)
    det_s = Osym.det(method="berkowitz")
    check(f"O6a n={nn}: expand(det Ω − Pf(Ω)²) = 0 — symbolic, full antisymmetric Ω "
          f"({n_syms} symbols)", sp.expand(det_s - pf_s ** 2) == 0)
    stamp = (f"O6a n={nn}: BOUNDARY (stamp): 50 random integer Ω, entries {{−9..9}}, "
             f"seed 905 (separate Random(905) per n)")
    print(f"    {stamp}")
    SEARCH_STAMPS.append(stamp)
    rng = random.Random(905)
    n_ok = 0
    for _ in range(50):
        Om = rand_omega(nn, rng, -9, 9)
        if Om.det() == pfaffian(Om) ** 2:
            n_ok += 1
    check(f"O6a n={nn}: det Ω = Pf(Ω)² on 50/50 random integer Ω (bit-exact)",
          n_ok == 50, f"{n_ok}/50")

for nn in (3, 5):
    Osym = sym_omega(nn, f"o{nn}")
    check(f"O6a n={nn} (odd): expand(det Ω) = 0 — symbolic, full antisymmetric Ω",
          sp.expand(Osym.det(method="berkowitz")) == 0)

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O6b — INVARIANCE OF det A UNDER CONGRUENCE (g,Ω) → (SᵀgS, SᵀΩS)")

O6B_SIGS = [(2, 2), (3, 1), (4, 0), (3, 3)]
for (p, q) in O6B_SIGS:
    nn = p + q
    eta, Js, Ks = make_soPQ(p, q)
    stamp = (f"O6b so({p},{q}): BOUNDARY (stamp): 1 random integer Ω (entries {{−5..5}}) + "
             f"20 random integer S (entries {{−3..3}}, det S ≠ 0 by retry), "
             f"seed 905{nn}{p} (separate Random per signature)")
    print(f"\n  {stamp}")
    SEARCH_STAMPS.append(stamp)
    rng = random.Random(int(f"905{nn}{p}"))
    n_draw = 0
    while True:
        n_draw += 1
        Om = rand_omega(nn, rng, -5, 5)
        if Om.det() != 0:
            break
        print(f"    Ω-attempt {n_draw}: det Ω = 0 (degenerate; printed and discarded) · "
              f"Ω(rows) = {[list(Om.row(i)) for i in range(nn)]}")
    print(f"    Ω-attempts until det Ω ≠ 0 (stamp): {n_draw}")
    A = eta * Om            # η⁻¹ = η
    detA = A.det()
    detg = eta.det()
    detOm = Om.det()
    print(f"    Ω(rows) = {[list(Om.row(i)) for i in range(nn)]}")
    print(f"    det g = {detg} · det Ω = {detOm} · det A = {detA} · "
          f"det g·det Ω = {detg * detOm}")
    n_bit = 0
    tot_tries = 0
    for t in range(20):
        S, tries = rand_S(nn, rng, -3, 3)
        tot_tries += tries
        gp = S.T * eta * S
        Op = S.T * Om * S
        Ap = gp.inv() * Op
        detAp = Ap.det()
        if detAp == detA and sp.srepr(sp.nsimplify(detAp)) == sp.srepr(sp.nsimplify(detA)):
            n_bit += 1
        if t < 3:
            dS = S.det()
            lhs = gp.det() * Op.det()
            rhs = (dS ** 4) * detg * detOm
            print(f"    example {t + 1}: det S = {dS} · (det S)⁴ = {dS ** 4} · "
                  f"det g′·det Ω′ = {lhs} · (det S)⁴·(det g·det Ω) = {rhs} · "
                  f"match: {'yes' if lhs == rhs else 'no'} · det A′ = {detAp} "
                  f"(det A = {detA})")
    check(f"O6b so({p},{q}): det A′ bit-matches det A on 20/20 random S "
          f"(sympy-equality + srepr-equality)", n_bit == 20, f"{n_bit}/20")
    fac_ok = True
    rng2 = random.Random(int(f"905{nn}{p}") + 1)
    for _ in range(20):
        S, _t = rand_S(nn, rng2, -3, 3)
        gp = S.T * eta * S
        Op = S.T * Om * S
        fac_ok &= (gp.det() * Op.det() == (S.det() ** 4) * detg * detOm)
    check(f"O6b so({p},{q}): det g′·det Ω′ = (det S)⁴·det g·det Ω on 20/20 independent S",
          fac_ok, f"attempts for non-degenerate S (first series): {tot_tries}")

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O6c — SIGNATURE LADDER TABLE: det A (generators · sums · 100 random)")

O6C_SIGS = [(2, 1), (2, 2), (3, 1), (4, 0), (3, 3)]
for (p, q) in O6C_SIGS:
    nn = p + q
    dim_so = nn * (nn - 1) // 2
    eta, Js, Ks = make_soPQ(p, q)
    gens = Js + Ks
    G = [M for _, M in gens]
    detg = eta.det()
    print(f"\n  so({p},{q}) — n = {nn} · det g = {detg} · J-moves: {len(Js)} · "
          f"K-moves: {len(Ks)}")
    check(f"O6c so({p},{q}): ηX antisymmetric for all {len(gens)} generators "
          f"(A = η⁻¹Ω correspondence)",
          all((eta * Gm).T == -(eta * Gm) for Gm in G))
    if nn % 2 == 1:
        Osym = sym_omega(nn, f"c{p}{q}")
        Asym = eta * Osym
        check(f"O6c so({p},{q}) n={nn} ODD: det A ≡ 0 — symbolic "
              f"(full antisymmetric Ω)",
              sp.expand(Asym.det(method="berkowitz")) == 0)
    print(f"    det A over the single generators (Ω = ηX, A = X):")
    for nm, X in gens:
        print(f"      {nm:>8}: det A = {X.det()}")
    if Js:
        AJ = sp.zeros(nn, nn)
        for _, X in Js:
            AJ = AJ + X
        print(f"    sum {{all J}} ({len(Js)} terms): det A = {AJ.det()}")
    else:
        print(f"    sum {{all J}}: J-moves 0 — the sum is empty (zero matrix), det A = 0")
    if Ks:
        AK = sp.zeros(nn, nn)
        for _, X in Ks:
            AK = AK + X
        print(f"    sum {{all K}} ({len(Ks)} terms): det A = {AK.det()}")
    else:
        print(f"    sum {{all K}}: K-moves 0 — the sum is empty (zero matrix), det A = 0")
    stamp = (f"O6c so({p},{q}): BOUNDARY (stamp): 100 random integer A = Σ cᵢ·Gᵢ, "
             f"cᵢ ∈ {{−3..3}}, seed 905 (separate Random(905) per signature)")
    print(f"    {stamp}")
    SEARCH_STAMPS.append(stamp)
    rng = random.Random(905)
    n_pos = n_zero = n_neg = 0
    d_min = None
    d_max = None
    for _ in range(100):
        cs = [rng.randint(-3, 3) for _ in range(dim_so)]
        Am = sp.zeros(nn, nn)
        for c, Gk in zip(cs, G):
            if c:
                Am = Am + c * Gk
        d = Am.det()
        if d > 0:
            n_pos += 1
        elif d < 0:
            n_neg += 1
        else:
            n_zero += 1
        d_min = d if d_min is None or d < d_min else d_min
        d_max = d if d_max is None or d > d_max else d_max
    check(f"O6c so({p},{q}): 100 random A processed (no truncation)",
          n_pos + n_zero + n_neg == 100)
    print(f"    HISTOGRAM of the sign of det A (100 random): positive {n_pos} · "
          f"zero {n_zero} · negative {n_neg} · min {d_min} · max {d_max} · "
          f"det g of the signature = {detg}")

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O6d — ORIENTATION FLIP: swap of two axes P — Pf(PᵀΩP) and det A by the numbers")

O6D = [((2, 2), (0, 1), (0, 2)), ((3, 3), (0, 1), (0, 3))]
for (p, q), same_pair, mixed_pair in O6D:
    nn = p + q
    eta, Js, Ks = make_soPQ(p, q)
    signs = [1] * p + [-1] * q
    for label, (i, j) in (("same-sign", same_pair), ("mixed-sign", mixed_pair)):
        P = sp.eye(nn)
        P[i, i] = 0
        P[j, j] = 0
        P[i, j] = 1
        P[j, i] = 1
        assert P.T * P == sp.eye(nn), "P is not an orthogonal permutation"
        assert (signs[i] == signs[j]) == (label == "same-sign"), "pair type mismatch"
        print(f"\n  so({p},{q}) — P = permutation of axes ({i},{j}), type: {label} "
              f"(η signs: {signs[i]:+d},{signs[j]:+d}) · det P = {P.det()} · "
              f"PᵀηP = η: {'yes' if P.T * eta * P == eta else 'no'}")
        stamp = (f"O6d so({p},{q}) pair ({i},{j}) [{label}]: BOUNDARY (stamp): 5 random "
                 f"integer Ω (entries {{−5..5}}), seed 9052{nn}{i}{j}")
        print(f"    {stamp}")
        SEARCH_STAMPS.append(stamp)
        rng = random.Random(int(f"9052{nn}{i}{j}"))
        n_flip = 0
        n_det_eq = 0
        for t in range(5):
            Om = rand_omega(nn, rng, -5, 5)
            pf0 = pfaffian(Om)
            Om1 = P.T * Om * P
            pf1 = pfaffian(Om1)
            dA0 = (eta * Om).det()
            dA1 = (eta * Om1).det()
            if pf1 == -pf0:
                n_flip += 1
            if dA1 == dA0:
                n_det_eq += 1
            print(f"    example {t + 1}: Pf(Ω) = {pf0} · Pf(PᵀΩP) = {pf1} · "
                  f"Pf-ratio = {sp.Rational(pf1, pf0) if pf0 != 0 else 'Pf(Ω)=0'} · "
                  f"det A = {dA0} · det A′ = {dA1}")
        print(f"    raw tally: Pf(PᵀΩP) = −Pf(Ω) in {n_flip}/5 · "
              f"det A′ = det A in {n_det_eq}/5")
        check(f"O6d so({p},{q}) pair ({i},{j}) [{label}]: 5 examples printed by the numbers",
              True)

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O6e — INTERSECTION WITH A² = ±𝟙 (S903 examples, address: S903_run.log block C)")

O6E = [
    ((2, 2), [-1, -1, 0, 0, 0, 0], -1,
     "A² = −𝟙 on (2,2): coeffs [−1,−1,0,0,0,0] (S903, first found)"),
    ((2, 2), [-2, -2, -2, -1, -1, 2], 1,
     "A² = +𝟙 on (2,2): coeffs [−2,−2,−2,−1,−1,2] (S903, first found)"),
    ((3, 3), None, 1,
     "A² = +𝟙 on (3,3): K(0,3)+K(1,4)+K(2,5) (S903, first found)"),
]
for (p, q), coeffs, s, title in O6E:
    nn = p + q
    eta, Js, Ks = make_soPQ(p, q)
    gens = Js + Ks
    if coeffs is not None:
        Am = sp.zeros(nn, nn)
        used = []
        for c, (nm, X) in zip(coeffs, gens):
            if c:
                Am = Am + c * X
                used.append(f"{c:+d}·{nm}")
    else:
        want = {"K(0,3)", "K(1,4)", "K(2,5)"}
        Am = sp.zeros(nn, nn)
        used = []
        for nm, X in gens:
            if nm in want:
                Am = Am + X
                used.append(f"+1·{nm}")
        assert len(used) == 3, "not all three K-generators found by name"
    assert (Am * eta + eta * Am.T) == sp.zeros(nn, nn), "A outside so(η)"
    assert Am * Am == s * sp.eye(nn), "A² ≠ s·𝟙 (reconstruction of the S903 example)"
    dA = Am.det()
    dA2 = (Am * Am).det()
    formula = sp.Integer(s) ** nn
    print(f"\n  {title}")
    print(f"    A = {' '.join(used)}")
    print(f"    A(rows) = {[list(Am.row(i)) for i in range(nn)]}")
    print(f"    det A = {dA}")
    print(f"    formula line: det(A²) = det({'+' if s == 1 else '−'}𝟙) = "
          f"({'+' if s == 1 else '−'}1)^{nn} = {formula} · sympy det(A·A) = {dA2} · "
          f"(det A)² = {dA ** 2}")
    check(f"O6e {title.split(':')[0]}: det(A²) = ({'+' if s == 1 else '−'}1)ⁿ = "
          f"(det A)² (sympy-exact)",
          dA2 == formula and dA ** 2 == formula)

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY S905 (raw lines; reading = an act of the court)")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}

  BOUND STAMPS:""")
for st in SEARCH_STAMPS:
    print(f"   · {st}")
print(f"""
  RAW LINES (no readings):
   (1) O6a: det Ω − Pf(Ω)² = 0 symbolic (n = 4, 6; full antisymmetric Ω) +
       50/50 random integers per n; det Ω ≡ 0 symbolic (n = 3, 5);
   (2) O6b: det A′ = det A bit-match 20/20 on every signature {{(2,2)·(3,1)·(4,0)·(3,3)}};
       det g′·det Ω′ = (det S)⁴·det g·det Ω — factors by the numbers for 3 examples — above;
   (3) O6c: det A ladder table {{(2,1)·(2,2)·(3,1)·(4,0)·(3,3)}} — single
       generators · sums of J/K · sign histograms on 100 random · det g alongside;
       (2,1): det A ≡ 0 symbolic — separate line above;
   (4) O6d: axis permutations (same-sign · mixed-sign) — Pf and det A by the numbers,
       5 examples per type — above;
   (5) O6e: det A by the numbers for three S903 examples A² = ±𝟙 + formula line
       det(A²) = (±1)ⁿ — above.
  HONEST TALLY: handles 0 · verdicts 0. Court = Omega.
""")
_logf.flush()
sys.exit(0 if not FAIL else 1)
