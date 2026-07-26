#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (linear algebra over Q; handles 0). W28-O7 — Pfaffian leg
#      of orientation: the transformation law Pf(SᵀΩS) = det(S)·Pf(Ω) (symbolic n=4 full +
#      50 random integers n=6) · explicit isometries SᵀηS = η of both types det S = ±1
#      (rational rotations/boosts · reflections) — Pf under det=+1 and under det=−1 by the numbers ·
#      Pf table (sums J/K · S903 square structures · 100 random, sign histogram) ·
#      orientation doubling of classes (single generators · pairs of the same centralizer class
#      with explicit det=−1 isometries) · tie-in with O6: det A = Pf(ηA)²/det g by bit-check.
#      ★BLINDNESS: the probe prints ONLY raw identities/numbers/histograms; reading = an act of the court.
"""
S906 (lane A, ed.2) — W28-O7: Pfaffian leg of orientation for the log-pair (η,Ω).

  Objects: signatures {(2,2)·(3,1)·(3,3)}, η = diag(±1); A ∈ so(η) (Xη + ηXᵀ = 0);
  Ω = ηA — antisymmetric (assert on every object). Pf — recursion on the first
  row (verbatim copy from S905; reference values Pf(2×2) and Pf(4×4) — stamped below).

  O7a — transformation law: Pf(SᵀΩS) = det(S)·Pf(Ω):
        n = 4 — symbolic FULL (full antisymmetric Ω on 6 symbols × full
        symbolic S on 16 symbols; expand of the difference to zero);
        n = 6 — 50 random integer pairs (Ω, S), seed 906, bit-equality.
        Then ISOMETRIES: explicit S with SᵀηS = η of both types on every signature
        (det = +1: block-rotation (3/5,4/5)-family and boost (5/3,4/3)-family,
        rationally exact; det = −1: single-axis reflection + composites) —
        5 examples per type per signature: Pf(SᵀΩS) against Pf(Ω) by the numbers
        (invariance under det=+1 · flip under det=−1 — raw lines).
  O7b — Pf table: Pf(ηA) for the sum {all J} · the sum {all K} · the square
        structures of S903 (on (2,2): A² = −𝟙 coeffs [−1,−1,0,0,0,0] and A² = +𝟙
        coeffs [−2,−2,−2,−1,−1,2] in the S903 generator basis; on (3,3): A² = +𝟙 =
        K(0,3)+K(1,4)+K(2,5)) · 100 random integer A (coeffs {−3..3}, seed 906)
        per signature: histogram of the sign of Pf {positive·zero·negative}; for
        Pf = 0 cases — rank of A as a raw line.
  O7c — orientation doubling of classes: Pf(ηX) for EACH single coordinate
        generator (J and K separately, all). Then pairs of generators of the SAME
        centralizer class (from S903: on (2,2)/(3,1) all singles — one class
        {2·(0,0,2)}; on (3,3) J-class {7·(3,3,1)} and K-class {7·(4,2,1)}):
        raw lines "Pf(X)=…, Pf(Y)=…, isometry R: RᵀηR=η, det R=−1,
        X→RᵀXR=Y? yes/no" — at least 3 explicit attempts per signature
        (axis reflection · permutation of same-sign axes).
  O7d — tie-in with O6: for ALL O7b objects, bit-check det A = Pf(ηA)²/det g
        (sympy-equality; every line PASS/FAIL).

  EXACT arithmetic: sympy Rational/Integer over Q; no tolerances. Mechanisms
  (make_soPQ/flat/pfaffian) — a verbatim copy from S903/S905 (reproducibility).
  No silent truncation: every bound is stamped in the printed lines.

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

# ── tee: all of stdout is duplicated into S906_run.log next to the script ──
_LOG_PATH = __file__.replace("\\", "/").rsplit("/", 1)[0] + "/S906_run.log"


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
check("handles 0 (pure algebra: Pfaffian · isometries SᵀηS=η · congruence · tie-in with det)",
      True)


# ═══════════════════════════════════════════════════════════════════════════════════════
# TOOLS (make_soPQ/flat/pfaffian — verbatim copy from S903/S905; isometry builders)
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


def rot(nn, i, j, c, s):
    """Block rotation in a same-sign plane (i,j): c² + s² = 1 (rational)."""
    R = sp.eye(nn)
    R[i, i] = c
    R[j, j] = c
    R[i, j] = -s
    R[j, i] = s
    return R


def boost(nn, i, j, ch, sh):
    """Block boost in a mixed-sign plane (i,j): ch² − sh² = 1 (rational)."""
    R = sp.eye(nn)
    R[i, i] = ch
    R[j, j] = ch
    R[i, j] = sh
    R[j, i] = sh
    return R


def refl(nn, i):
    """Reflection of a single axis i: diag(1,…,−1,…,1)."""
    R = sp.eye(nn)
    R[i, i] = -1
    return R


def perm2(nn, i, j):
    """Permutation of two axes (i,j)."""
    R = sp.eye(nn)
    R[i, i] = 0
    R[j, j] = 0
    R[i, j] = 1
    R[j, i] = 1
    return R


def sign_word(v):
    if v > 0:
        return "positive"
    if v < 0:
        return "negative"
    return "zero"


Q = sp.Rational
SIGS = [(2, 2), (3, 1), (3, 3)]
SEARCH_STAMPS = []
OBJS = {}   # signature → list of (label, A, Pf(ηA)) for O7d

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O7a — TRANSFORMATION LAW: Pf(SᵀΩS) = det(S)·Pf(Ω)")

# recursion reference values (stamp: explicit small n, verbatim copy of S905 reference values)
_O2 = sym_omega(2, "b2")
check("O7a reference: Pf(2×2) = ω₀₁ and det(2×2) = ω₀₁² (symbolic)",
      pfaffian(_O2) == _O2[0, 1] and sp.expand(_O2.det() - _O2[0, 1] ** 2) == 0)
_O4t = sym_omega(4, "b4")
_pf4_ref = _O4t[0, 1] * _O4t[2, 3] - _O4t[0, 2] * _O4t[1, 3] + _O4t[0, 3] * _O4t[1, 2]
check("O7a reference: Pf(4×4) = ω₀₁ω₂₃ − ω₀₂ω₁₃ + ω₀₃ω₁₂ (recursion = classical form)",
      sp.expand(pfaffian(_O4t) - _pf4_ref) == 0)

# n = 4: FULL symbolic (full Ω on 6 symbols, full S on 16 symbols)
_Om4 = sym_omega(4, "t4")
_S4 = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"s_{i}_{j}"))
_C4 = _S4.T * _Om4 * _S4
check("O7a n=4: (SᵀΩS)ᵀ = −(SᵀΩS) — symbolic (congruence preserves antisymmetry)",
      sp.expand(_C4.T + _C4) == sp.zeros(4, 4))
check("O7a n=4: expand(Pf(SᵀΩS) − det(S)·Pf(Ω)) = 0 — symbolic FULL "
      "(Ω: 6 symbols · S: 16 symbols)",
      sp.expand(pfaffian(_C4) - _S4.det(method="berkowitz") * pfaffian(_Om4)) == 0)

# n = 6: 50 random integer pairs (Ω, S)
stamp = ("O7a n=6: BOUNDARY (stamp): 50 random integer pairs (Ω entries {−9..9}, "
         "S entries {−4..4}, det S unrestricted), seed 906 (Random(906))")
print(f"    {stamp}")
SEARCH_STAMPS.append(stamp)
_rng = random.Random(906)
_n_ok = 0
for _ in range(50):
    Om = rand_omega(6, _rng, -9, 9)
    S = sp.Matrix(6, 6, lambda i, j: _rng.randint(-4, 4))
    if pfaffian(S.T * Om * S) == S.det() * pfaffian(Om):
        _n_ok += 1
check("O7a n=6: Pf(SᵀΩS) = det(S)·Pf(Ω) on 50/50 random integer pairs (bit-exact)",
      _n_ok == 50, f"{_n_ok}/50")

# ── ISOMETRIES: explicit S with SᵀηS = η, det = ±1, 5 examples per type per signature ──
print("\n  ISOMETRIES SᵀηS = η (rationally exact: rotations (c,s) with c²+s²=1 · "
      "boosts (ch,sh) with ch²−sh²=1 · reflections · composites):")

ISO_PLUS = {
    (2, 2): [
        ("rotation(0,1; 3/5,4/5)",              lambda nn: rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
        ("rotation(2,3; 5/13,12/13)",           lambda nn: rot(nn, 2, 3, Q(5, 13), Q(12, 13))),
        ("boost(0,2; 5/3,4/3)",                 lambda nn: boost(nn, 0, 2, Q(5, 3), Q(4, 3))),
        ("boost(1,3; 13/5,12/5)",               lambda nn: boost(nn, 1, 3, Q(13, 5), Q(12, 5))),
        ("rotation(0,1;8/17,15/17)·boost(1,2;17/8,15/8)",
         lambda nn: rot(nn, 0, 1, Q(8, 17), Q(15, 17)) * boost(nn, 1, 2, Q(17, 8), Q(15, 8))),
    ],
    (3, 1): [
        ("rotation(0,1; 3/5,4/5)",              lambda nn: rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
        ("rotation(1,2; 5/13,12/13)",           lambda nn: rot(nn, 1, 2, Q(5, 13), Q(12, 13))),
        ("boost(0,3; 5/3,4/3)",                 lambda nn: boost(nn, 0, 3, Q(5, 3), Q(4, 3))),
        ("boost(2,3; 13/5,12/5)",               lambda nn: boost(nn, 2, 3, Q(13, 5), Q(12, 5))),
        ("rotation(0,2;8/17,15/17)·boost(1,3;17/8,15/8)",
         lambda nn: rot(nn, 0, 2, Q(8, 17), Q(15, 17)) * boost(nn, 1, 3, Q(17, 8), Q(15, 8))),
    ],
    (3, 3): [
        ("rotation(0,1; 3/5,4/5)",              lambda nn: rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
        ("rotation(3,4; 5/13,12/13)",           lambda nn: rot(nn, 3, 4, Q(5, 13), Q(12, 13))),
        ("boost(0,3; 5/3,4/3)",                 lambda nn: boost(nn, 0, 3, Q(5, 3), Q(4, 3))),
        ("boost(2,5; 13/5,12/5)",               lambda nn: boost(nn, 2, 5, Q(13, 5), Q(12, 5))),
        ("rotation(1,2;7/25,24/25)·boost(0,4;25/7,24/7)",
         lambda nn: rot(nn, 1, 2, Q(7, 25), Q(24, 25)) * boost(nn, 0, 4, Q(25, 7), Q(24, 7))),
    ],
}
ISO_MINUS = {
    (2, 2): [
        ("reflection of axis 0",                lambda nn: refl(nn, 0)),
        ("reflection of axis 1",                lambda nn: refl(nn, 1)),
        ("reflection of axis 2",                lambda nn: refl(nn, 2)),
        ("reflection of axis 3",                lambda nn: refl(nn, 3)),
        ("reflection of axis 0 · boost(0,2;5/3,4/3)",
         lambda nn: refl(nn, 0) * boost(nn, 0, 2, Q(5, 3), Q(4, 3))),
    ],
    (3, 1): [
        ("reflection of axis 0",                lambda nn: refl(nn, 0)),
        ("reflection of axis 1",                lambda nn: refl(nn, 1)),
        ("reflection of axis 2",                lambda nn: refl(nn, 2)),
        ("reflection of axis 3",                lambda nn: refl(nn, 3)),
        ("reflection of axis 2 · rotation(0,1;3/5,4/5)",
         lambda nn: refl(nn, 2) * rot(nn, 0, 1, Q(3, 5), Q(4, 5))),
    ],
    (3, 3): [
        ("reflection of axis 0",                lambda nn: refl(nn, 0)),
        ("reflection of axis 3",                lambda nn: refl(nn, 3)),
        ("reflection of axis 5",                lambda nn: refl(nn, 5)),
        ("reflection of axis 1 · rotation(0,2;3/5,4/5)",
         lambda nn: refl(nn, 1) * rot(nn, 0, 2, Q(3, 5), Q(4, 5))),
        ("reflection of axis 4 · boost(1,4;5/3,4/3)",
         lambda nn: refl(nn, 4) * boost(nn, 1, 4, Q(5, 3), Q(4, 3))),
    ],
}

for (p, q) in SIGS:
    nn = p + q
    eta, Js, Ks = make_soPQ(p, q)
    stamp = (f"O7a-isometries so({p},{q}): BOUNDARY (stamp): 5 explicit S per type det=±1; "
             f"random Ω entries {{−5..5}}, Pf(Ω) ≠ 0 by retry, seed 906{nn}{p}")
    print(f"\n  so({p},{q}) — n = {nn}: {stamp}")
    SEARCH_STAMPS.append(stamp)
    rng = random.Random(int(f"906{nn}{p}"))
    for det_want, iso_list, word in ((1, ISO_PLUS[(p, q)], "det=+1"),
                                     (-1, ISO_MINUS[(p, q)], "det=−1")):
        n_hold = 0
        for t, (label, build) in enumerate(iso_list):
            S = build(nn)
            assert S.T * eta * S == eta, f"S is not an isometry: {label}"
            dS = S.det()
            assert dS == det_want, f"det S ≠ {det_want}: {label}"
            n_draw = 0
            while True:
                n_draw += 1
                Om = rand_omega(nn, rng, -5, 5)
                pf0 = pfaffian(Om)
                if pf0 != 0:
                    break
            pf1 = pfaffian(S.T * Om * S)
            if pf1 == det_want * pf0:
                n_hold += 1
            print(f"    [{word}] example {t + 1}: S = {label} · det S = {dS} · "
                  f"SᵀηS = η: yes · Ω-draws until Pf≠0: {n_draw} · Pf(Ω) = {pf0} · "
                  f"Pf(SᵀΩS) = {pf1} · ratio = {sp.Rational(pf1, pf0)}")
        check(f"O7a-isometries so({p},{q}) {word}: Pf(SᵀΩS) = ({det_want})·Pf(Ω) "
              f"on 5/5 explicit examples", n_hold == 5, f"{n_hold}/5")

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O7b — Pf TABLE: Pf(ηA) for sums J/K · S903 square structures · 100 random")

SQ_TABLE = {
    (2, 2): [
        ([-1, -1, 0, 0, 0, 0], -1,
         "A² = −𝟙 (2,2): coeffs [−1,−1,0,0,0,0] in the S903 generator basis"),
        ([-2, -2, -2, -1, -1, 2], 1,
         "A² = +𝟙 (2,2): coeffs [−2,−2,−2,−1,−1,2] in the S903 generator basis"),
    ],
    (3, 1): [],
    (3, 3): [
        (None, 1, "A² = +𝟙 (3,3): K(0,3)+K(1,4)+K(2,5) (S903)"),
    ],
}

for (p, q) in SIGS:
    nn = p + q
    dim_so = nn * (nn - 1) // 2
    eta, Js, Ks = make_soPQ(p, q)
    gens = Js + Ks
    G = [M for _, M in gens]
    detg = eta.det()
    print(f"\n  so({p},{q}) — n = {nn} · det g = {detg} · J-moves: {len(Js)} · "
          f"K-moves: {len(Ks)}")
    objs = []

    # sums {all J} and {all K}
    AJ = sp.zeros(nn, nn)
    for _, X in Js:
        AJ = AJ + X
    AK = sp.zeros(nn, nn)
    for _, X in Ks:
        AK = AK + X
    objs.append((f"sum {{all J}} ({len(Js)} terms)", AJ))
    objs.append((f"sum {{all K}} ({len(Ks)} terms)", AK))

    # S903 square structures (reconstruction with asserts)
    for coeffs, s, title in SQ_TABLE[(p, q)]:
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
        objs.append((title + f" [A = {' '.join(used)}]", Am))

    # Pf for the named objects (Ω = ηA antisymmetric — assert on each)
    named_pf = []
    for label, Am in objs:
        Om = eta * Am
        assert Om.T == -Om, f"Ω = ηA not antisymmetric: {label}"
        pf = pfaffian(Om)
        named_pf.append((label, Am, pf))
        print(f"    {label}: Pf(ηA) = {pf} ({sign_word(pf)})")

    # 100 random integer A
    stamp = (f"O7b so({p},{q}): BOUNDARY (stamp): 100 random integer A = Σ cᵢ·Gᵢ, "
             f"cᵢ ∈ {{−3..3}}, seed 906 (separate Random(906) per signature)")
    print(f"    {stamp}")
    SEARCH_STAMPS.append(stamp)
    rng = random.Random(906)
    n_pos = n_zero = n_neg = 0
    zero_rows = []
    rand_pf = []
    for k in range(100):
        cs = [rng.randint(-3, 3) for _ in range(dim_so)]
        Am = sp.zeros(nn, nn)
        for c, Gk in zip(cs, G):
            if c:
                Am = Am + c * Gk
        Om = eta * Am
        assert Om.T == -Om, "Ω = ηA not antisymmetric (random A)"
        pf = pfaffian(Om)
        rand_pf.append((f"rand#{k + 1:03d}", Am, pf))
        if pf > 0:
            n_pos += 1
        elif pf < 0:
            n_neg += 1
        else:
            n_zero += 1
            zero_rows.append((k + 1, Am.rank()))
    check(f"O7b so({p},{q}): 100 random A processed (no truncation)",
          n_pos + n_zero + n_neg == 100)
    print(f"    HISTOGRAM of the sign of Pf(ηA) (100 random): positive {n_pos} · "
          f"zero {n_zero} · negative {n_neg}")
    if zero_rows:
        for k, rk in zero_rows:
            print(f"      Pf = 0 at rand#{k:03d}: rank A = {rk}")
    else:
        print(f"      no Pf = 0 cases (empty)")
    OBJS[(p, q)] = named_pf + rand_pf

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O7c — ORIENTATION DOUBLING OF CLASSES: single generators · pairs of one class")

# centralizer classes of single generators — reference from S903 (block A, S903_run.log):
REF_CLASS = {
    (2, 2): "all 6 singles: one class {2·(0,0,2)}",
    (3, 1): "all 6 singles: one class {2·(0,0,2)}",
    (3, 3): "J-class {7·(3,3,1)} (all 6 J) · K-class {7·(4,2,1)} (all 9 K)",
}

# attempts: (X-label, X-generator expression, Y-label, Y-generator expression, R-label, R-builder)
# X, Y — Ω-objects (η·generator); the expressions below are on generator names.
ATTEMPTS = {
    (2, 2): [
        ("η·K(0,2)", ("K(0,2)", 1), "η·K(1,2)", ("K(1,2)", 1),
         "permutation of axes (0,1)", lambda nn: perm2(nn, 0, 1)),
        ("η·K(0,2)", ("K(0,2)", 1), "η·K(0,3)", ("K(0,3)", 1),
         "permutation of axes (2,3)", lambda nn: perm2(nn, 2, 3)),
        ("η·J(0,1)", ("J(0,1)", 1), "η·(−J(0,1))", ("J(0,1)", -1),
         "reflection of axis 0", lambda nn: refl(nn, 0)),
        ("η·J(0,1)", ("J(0,1)", 1), "η·J(2,3)", ("J(2,3)", 1),
         "reflection of axis 0", lambda nn: refl(nn, 0)),
    ],
    (3, 1): [
        ("η·J(0,1)", ("J(0,1)", 1), "η·J(0,2)", ("J(0,2)", 1),
         "permutation of axes (1,2)", lambda nn: perm2(nn, 1, 2)),
        ("η·K(0,3)", ("K(0,3)", 1), "η·K(1,3)", ("K(1,3)", 1),
         "permutation of axes (0,1)", lambda nn: perm2(nn, 0, 1)),
        ("η·J(0,1)", ("J(0,1)", 1), "η·(−J(0,1))", ("J(0,1)", -1),
         "reflection of axis 1", lambda nn: refl(nn, 1)),
        ("η·K(2,3)", ("K(2,3)", 1), "η·(−K(2,3))", ("K(2,3)", -1),
         "reflection of axis 3", lambda nn: refl(nn, 3)),
    ],
    (3, 3): [
        ("η·J(0,1)", ("J(0,1)", 1), "η·J(0,2)", ("J(0,2)", 1),
         "permutation of axes (1,2)", lambda nn: perm2(nn, 1, 2)),
        ("η·J(3,4)", ("J(3,4)", 1), "η·J(3,5)", ("J(3,5)", 1),
         "permutation of axes (4,5)", lambda nn: perm2(nn, 4, 5)),
        ("η·K(0,3)", ("K(0,3)", 1), "η·K(1,3)", ("K(1,3)", 1),
         "permutation of axes (0,1)", lambda nn: perm2(nn, 0, 1)),
        ("η·K(0,3)", ("K(0,3)", 1), "η·K(0,4)", ("K(0,4)", 1),
         "permutation of axes (3,4)", lambda nn: perm2(nn, 3, 4)),
        ("η·J(0,1)", ("J(0,1)", 1), "η·(−J(0,1))", ("J(0,1)", -1),
         "reflection of axis 0", lambda nn: refl(nn, 0)),
        ("η·J(0,1)", ("J(0,1)", 1), "η·J(3,4)", ("J(3,4)", 1),
         "reflection of axis 2", lambda nn: refl(nn, 2)),
    ],
}

for (p, q) in SIGS:
    nn = p + q
    eta, Js, Ks = make_soPQ(p, q)
    gens = Js + Ks
    by_name = {nm: X for nm, X in gens}
    print(f"\n  so({p},{q}) — n = {nn}:")
    print(f"    Pf(ηX) for EACH single generator "
          f"(S903 class reference: {REF_CLASS[(p, q)]}):")
    for nm, X in gens:
        Om = eta * X
        assert Om.T == -Om, f"ηX not antisymmetric: {nm}"
        print(f"      {nm:>8}: Pf(ηX) = {pfaffian(Om)}")
    n_att = 0
    for (xl, (xn, xc), yl, (yn, yc), rl, rb) in ATTEMPTS[(p, q)]:
        X = eta * (xc * by_name[xn])
        Y = eta * (yc * by_name[yn])
        R = rb(nn)
        assert X.T == -X and Y.T == -Y, "X or Y not antisymmetric"
        iso_ok = (R.T * eta * R == eta)
        dR = R.det()
        conj = R.T * X * R
        hit = (conj == Y)
        n_att += 1
        print(f"    attempt {n_att}: Pf(X={xl}) = {pfaffian(X)}, Pf(Y={yl}) = {pfaffian(Y)}, "
              f"isometry R = {rl}: RᵀηR = η? {'yes' if iso_ok else 'no'}, det R = {dR}, "
              f"X→RᵀXR = Y? {'yes' if hit else 'no'}")
        assert iso_ok and dR == -1, "R is not a det=−1 isometry"
    check(f"O7c so({p},{q}): {n_att} explicit attempts (minimum 3) printed as raw lines",
          n_att >= 3, f"attempts: {n_att}")

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("O7d — TIE-IN WITH O6: det A = Pf(ηA)²/det g — bit-check on ALL O7b objects")

for (p, q) in SIGS:
    nn = p + q
    eta, Js, Ks = make_soPQ(p, q)
    detg = eta.det()
    print(f"\n  so({p},{q}) — det g = {detg} · O7b objects: {len(OBJS[(p, q)])}")
    n_pass = 0
    for label, Am, pf in OBJS[(p, q)]:
        dA = Am.det()
        rhs = pf ** 2 / detg
        ok = (dA == rhs)
        if ok:
            n_pass += 1
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}: det A = {dA} · Pf(ηA) = {pf} · "
              f"Pf²/det g = {rhs}")
    check(f"O7d so({p},{q}): det A = Pf(ηA)²/det g on {n_pass}/{len(OBJS[(p, q)])} "
          f"O7b objects (sympy-equality)", n_pass == len(OBJS[(p, q)]))

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY S906 (raw lines; reading = an act of the court)")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}

  BOUND STAMPS:""")
for st in SEARCH_STAMPS:
    print(f"   · {st}")
print(f"""
  RAW LINES (no readings):
   (1) O7a: Pf(SᵀΩS) = det(S)·Pf(Ω) — symbolic FULL n=4 (Ω 6 symbols × S 16
       symbols) + 50/50 random integer pairs n=6; isometries SᵀηS = η of both types
       det = ±1, 5 explicit examples per type per signature — Pf by the numbers above;
   (2) O7b: Pf(ηA) table {{sums J/K · S903 square structures · 100 random}}
       on {{(2,2)·(3,1)·(3,3)}} — values, sign histograms, ranks of Pf=0 cases above;
   (3) O7c: Pf(ηX) for every single coordinate generator + explicit attempts
       to connect pairs of one centralizer class by det=−1 isometries — raw lines above;
   (4) O7d: det A = Pf(ηA)²/det g — per-object bit-check PASS/FAIL on all
       O7b objects — above.
  HONEST TALLY: handles 0 · verdicts 0. Court = Omega.
""")
_logf.flush()
sys.exit(0 if not FAIL else 1)
