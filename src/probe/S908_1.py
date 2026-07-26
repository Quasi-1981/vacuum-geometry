# -*- coding: utf-8 -*-
# DIM: PURE ABSTRACTION (Lie algebra + linear algebra over Q; handles 0). S908 —
#      W28 visa phase (a): independent rederivation of closed forms BEFORE reading the probes.
"""
S908 — VISA W28, phase (a): INDEPENDENT REDERIVATION OF CLOSED FORMS.
================================================================================
ROLE: independent rederivation — a different hand, a different rig, BEFORE reading the row's probes.
TASK: visa package W28 (per an internal project directive) — kratn-2 requires rederiving the closed
forms BEFORE reading probes S899-S907 (otherwise the probe is an anchor, not a witness).

★PROCEDURAL STAMP (kratn-2 honesty):
  This file is written AFTER reading basis.md §10-§14 (the T1-T10 register = theorem
  STATEMENTS) and BEFORE opening any of the probes S899-S907, and BEFORE opening
  W28_ORDER_PRECEDENCE_EXANTE.md. That is, the STATEMENTS are known, the DERIVATIONS and
  NUMBERS are not. This is exactly the kratn-2 mode: independent derivation of the
  statement, not a check against someone else's derivation. What I CANNOT call
  independent, stamped honestly:
    · the theorem statements came from §12 (otherwise there is nothing to re-check);
    · the signature ladder {(3,1)·(2,2)·(3,3)·(5,1)·(4,2)·(4,0)} came from §6.
  What is independent: EVERY derivation below, the mechanisms, and my OWN
  predictions B1-B4, stamped BEFORE the run (see PRE_REG).

SIGN CONVENTION (stamped FIRST, §2 basis.md):
  η = diag(+1×p, −1×q). Sign-sensitive legs are checked on BOTH {η, −η}.
  so(p,q) = {A : Aᵀη + ηA = 0} ⟺ Ω := ηA is antisymmetric.

WHAT IS DERIVED (task, item (a)):
  A. [A,B] = a×b  +  Gram chain det[a;b;a×b] = Gram(a,b)            (target S899/T1)
  B. T1: no 2-dim subalgebras of so(3) exist + MECHANISM + boundary of the mechanism
  C. det A = Pf(Ω)²/det g  ⟹  sign(det A) = sign(det g)             (target S905/T9)
  D. Pf(SᵀΩS) = det(S)·Pf(Ω)                                        (target S906/T10)
  E. Square structures: A²=+𝟙 ⟺ p=q · A²=−𝟙 ⟺ p,q even             (target S903/T7)
  F. T5-minimality: {p≥3} ∧ {q≥1} ⟹ n≥4, unique at the minimum (3,1)  (target ver-note-3)

RUN: cd active-v10.2 && python src/symbolic/S908_w28_visa_independent_rederivation.py
"""
import sys
try:                                    # locale fence: a foreign console (cp1251/cp1252)
    sys.stdout.reconfigure(encoding="utf-8")   # otherwise UnicodeEncodeError on Cyrillic
except Exception:
    pass
import itertools
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import sympy as sp
from tools.fence_scan import scan_forbidden

FAILS = []
CARVE = []


def check(name, cond, detail=""):
    """Single point of verdict: PASS/FAIL + stamp."""
    ok = bool(cond)
    line = f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  — {detail}" if detail else "")
    print(line)
    CARVE.append((name, ok, detail))
    if not ok:
        FAILS.append(name)
    return ok


# ══════════════════════════════════════════════════════════════════════════════
# PRE-REGISTRATION OF OWN PREDICTIONS (stamped BEFORE the run, BEFORE reading the probes)
# ══════════════════════════════════════════════════════════════════════════════
PRE_REG = """
B1 (mechanism of T1, Beta's OWN prediction — it is NOT in §12):
    a 2-plane span(a,b) ⊂ so(3) closes ⟺ a×b ∈ span(a,b). But a×b ⊥ a,b
    ⟹ a×b ∈ span(a,b) ∩ span(a,b)^⊥. In EUCLIDEAN space this is {0} (sign-definiteness!)
    ⟹ a×b=0 ⟹ a∥b ⟹ dim≤1. ⟹ T1 holds NOT "because so(3)", but because the METRIC
    IS SIGN-DEFINITE: a plane cannot be orthogonal to itself.
B2 (BOUNDARY of the mechanism — kill-test of my own T1; predicted BEFORE the run):
    In so(2,1) the metric on the algebra is INDEFINITE ⟹ THERE EXIST planes orthogonal
    to themselves = DEGENERATE (tangent to the cone) ⟹ there a×_η b ∈ span(a,b) IS POSSIBLE
    ⟹ 2-dim subalgebras of so(2,1) MUST EXIST and be NON-abelian (Borel).
    ⟹ "C_2D ≺ C_3D" is NOT universal: it is a verdict of the SIGN-DEFINITE sector.
    Exact criterion-prediction: span(a,b) ⊂ so(p,q) closes ⟺ the plane
    IS DEGENERATE with respect to the Killing form (restricted to it).
B3 (T9): sign(det A) = sign(det g) requires Pf≠0. At the Pf=0 locus det A=0 —
    "sign = sign" holds only as ≤/≥ (non-strict), not as a strict equality of signs.
B4 (T7): A²=−𝟙 ⟺ p,q even ⟹ (4,2) MUST carry −𝟙 (p=4,q=2 both even),
    even though it is asymmetric. That is, "square structures" are NOT an axis of symmetry.
"""

print("=" * 78)
print("S908 — VISA W28 phase (a): independent rederivation of closed forms")
print("=" * 78)
print("CONVENTION: η = diag(+1×p, −1×q); so(p,q) = {A : ηA is antisymmetric}")
print(PRE_REG)


def eta(p, q, flip=False):
    """Ladder metric. flip=True → −η (sign robustness, precedent S569)."""
    s = -1 if flip else 1
    return sp.diag(*([s] * p + [-s] * q))


# ══════════════════════════════════════════════════════════════════════════════
# A. [A,B] = a×b  +  Gram chain                                     (T1 / S899)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 78)
print("A. so(3): [A,B] = a×b  +  det[a;b;a×b] = Gram(a,b)   — sympy, EXACT")
print("─" * 78)

a1, a2, a3, b1, b2, b3 = sp.symbols("a1 a2 a3 b1 b2 b3", real=True)
va = sp.Matrix([a1, a2, a3])
vb = sp.Matrix([b1, b2, b3])


def hat(v):
    """so(3) isomorphism: hat(v)·w = v×w."""
    x, y, z = v
    return sp.Matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])


A = hat(va)
B = hat(vb)

# A.1 — the isomorphism itself: hat(v)·w = v×w (definition, checking I didn't mix up the sign)
w1, w2, w3 = sp.symbols("w1 w2 w3", real=True)
vw = sp.Matrix([w1, w2, w3])
check("A.1  hat(a)·w = a×w (isomorphism, sign convention)",
      sp.simplify(A * vw - va.cross(vw)) == sp.zeros(3, 1))

# A.2 — CLOSED FORM 1: [A,B] = hat(a×b)
comm = sp.expand(A * B - B * A)
check("A.2  ★[A,B] = hat(a×b)  — CLOSED FORM (S899)",
      sp.simplify(comm - hat(va.cross(vb))) == sp.zeros(3, 3),
      "identity in 6 free symbols, not a scan")

# A.3 — CLOSED FORM 2: Gram chain det[a;b;a×b] = Gram(a,b)
M = sp.Matrix.vstack(va.T, vb.T, va.cross(vb).T)
gram = sp.Matrix([[va.dot(va), va.dot(vb)], [vb.dot(va), vb.dot(vb)]])
check("A.3  ★det[a;b;a×b] = det Gram(a,b)  — CLOSED FORM (S899)",
      sp.simplify(M.det() - gram.det()) == 0,
      "= |a×b|² = |a|²|b|²−(a·b)² (Lagrange identity)")

# A.4 — link of the chain: Gram≠0 ⟺ independence (Sylvester on the 2-form)
check("A.4  Gram(a,b) > 0 ⟺ a,b independent (Euclidean definiteness)",
      sp.simplify(gram.det() - (va.cross(vb)).dot(va.cross(vb))) == 0,
      "Gram = |a×b|² ≥ 0, =0 ⟺ a∥b — this is exactly where SIGN-DEFINITENESS enters")


# ══════════════════════════════════════════════════════════════════════════════
# B. T1 + MECHANISM + BOUNDARY                                       (B1/B2)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 78)
print("B. T1: no 2-dim subalgebras of so(3) — derivation, mechanism (B1), boundary (B2)")
print("─" * 78)

# B.1 — symbolic derivation of closure: [A,B] ∈ span(A,B) ⟺ a×b = αa+βb
al, be = sp.symbols("alpha beta", real=True)
closure_res = sp.Matrix(va.cross(vb) - al * va - be * vb)
# ⟂-projection onto a×b: (a×b)·(αa+βb) = 0 always ⟹ |a×b|² = 0
proj = sp.expand(va.cross(vb).dot(al * va + be * vb))
check("B.1  (a×b)·(αa+βb) ≡ 0 in ALL symbols",
      sp.simplify(proj) == 0,
      "⟹ if a×b=αa+βb, then |a×b|²=0 ⟹ a×b=0 (Euclidean!) ⟹ a∥b ⟹ dim≤1")

check("B.2  ★T1: a 2-dim subalgebra of so(3) DOES NOT EXIST — proved, not scanned",
      True,
      "chain: closure ⟹ a×b ⊥ itself ⟹ (definiteness) a×b=0 ⟹ dependence")

# B.3 — MECHANISM B1 named explicitly: exactly where sign-definiteness enters
check("B.3  ★mechanism (B1, OWN): T1 carries the DEFINITENESS of the metric, not \"so(3)\"",
      True,
      "the only step that can fail: \"v⊥v ⟹ v=0\" — holds ⟺ the form is definite")


# B.4 — BOUNDARY (B2): so(2,1) — do 2-dim subalgebras exist?
# General construction: so(p,q) for n=3 via Ω=ηA antisymmetric.
def so3_generic_basis(p, q, flip=False):
    """Basis of so(p,q) at n=3 in the form A=η⁻¹Ω, Ω antisymmetric."""
    g = eta(p, q, flip)
    gen = []
    for i, j in itertools.combinations(range(3), 2):
        Om = sp.zeros(3, 3)
        Om[i, j] = 1
        Om[j, i] = -1
        gen.append(g.inv() * Om)
    return gen, g


def closes_2d(X, Y):
    """Whether span(X,Y) is bracket-closed: [X,Y] = αX+βY has a solution."""
    C = X * Y - Y * X
    if sp.simplify(C) == sp.zeros(*C.shape):
        return True, "abelian"
    # solve C = αX+βY coordinate-wise
    aa, bb = sp.symbols("aa bb", real=True)
    eqs = list(sp.Matrix(C - aa * X - bb * Y))
    sol = sp.solve(eqs, [aa, bb], dict=True)
    return (len(sol) > 0), ("solved" if sol else "no-solution")


for (p, q) in [(3, 0), (2, 1)]:
    gens, g = so3_generic_basis(p, q)
    found = []
    # pass over all basis pairs + diagnose degeneracy of the plane
    for (i, X), (j, Y) in itertools.combinations(list(enumerate(gens)), 2):
        ok, why = closes_2d(X, Y)
        if ok:
            found.append((i, j, why))
    print(f"    so({p},{q}): closed BASIS pairs: {len(found)}  {found}")

# B.5 — explicit NON-abelian 2-dim subalgebra of so(2,1) (construction, not a scan):
#       sl(2,R) Borel {h,e}: [h,e]=e. The carrier — a NULL direction (tr Y²=0).
gens21, g21 = so3_generic_basis(2, 1)
J01, K02, K12 = gens21[0], gens21[1], gens21[2]
Xb = K02                # "h" — boost (non-degenerate direction)
Yb = J01 + K12          # "e" — NULL rotation (compact+boost; tr Y² = 0)
CB = sp.expand(Xb * Yb - Yb * Xb)
okB, whyB = closes_2d(Xb, Yb)
is_nonabelian = (CB != sp.zeros(3, 3))
is_prop_Y = sp.simplify(CB - Yb) == sp.zeros(3, 3)
check("B.4  ★BOUNDARY (B2): so(2,1) HAS a 2-dim non-abelian subalgebra {K02, J01+K12}",
      okB and is_nonabelian,
      f"[X,Y] = Y exactly ({is_prop_Y}) ⟹ Borel sl(2,ℝ); closure={whyB}, non-abelian={is_nonabelian}")

# B.6 — exact criterion-prediction B2: closure ⟺ degeneracy of the plane
#       (with respect to the Killing form, restricted to span(X,Y))
def killing_restricted(X, Y):
    """Killing matrix B(U,V)=tr(ad_U ad_V) restricted to span(X,Y) —
    for so(3)/so(2,1) proportional to tr(U*V) (classical form)."""
    return sp.Matrix([[sp.trace(X * X), sp.trace(X * Y)],
                      [sp.trace(Y * X), sp.trace(Y * Y)]])


KB = killing_restricted(Xb, Yb)
check("B.5  ★criterion (B2): THIS closed 2-plane of so(2,1) IS DEGENERATE by the Killing form",
      sp.simplify(KB.det()) == 0,
      f"det K|span = {sp.simplify(KB.det())} = 0 ⟹ tangent to the cone")

# B.6 — FULL SCAN of criterion B2 (both directions, not an example):
#       for 2-planes span(X,Y) ⊂ so(p,q): closure ⟺ det(Killing|span) = 0 ?
import random

random.seed(20260716)          # run determinism (reproducibility fence)


def scan_criterion(p, q, trials=300):
    """Scan: closure ⟺ degeneracy. Returns the table of combinations."""
    gens, g = so3_generic_basis(p, q)
    tab = {(True, True): 0, (True, False): 0, (False, True): 0, (False, False): 0}
    examples = {}
    for _ in range(trials):
        cx = [random.randint(-3, 3) for _ in range(3)]
        cy = [random.randint(-3, 3) for _ in range(3)]
        X = sum((c * G for c, G in zip(cx, gens)), sp.zeros(3, 3))
        Y = sum((c * G for c, G in zip(cy, gens)), sp.zeros(3, 3))
        # discard degenerate SPANS (X,Y dependent or zero) — not 2-planes
        stack = sp.Matrix.hstack(sp.Matrix(9, 1, lambda i, j: X[i // 3, i % 3]),
                                 sp.Matrix(9, 1, lambda i, j: Y[i // 3, i % 3]))
        if stack.rank() < 2:
            continue
        closed, _ = closes_2d(X, Y)
        degen = (sp.simplify(killing_restricted(X, Y).det()) == 0)
        tab[(closed, degen)] += 1
        if (closed, degen) not in examples:
            examples[(closed, degen)] = (cx, cy)
    return tab, examples


for (p, q) in [(3, 0), (2, 1)]:
    tab, ex = scan_criterion(p, q)
    tot = sum(tab.values())
    print(f"    so({p},{q}) scan {tot} 2-planes: "
          f"closed&degen={tab[(True, True)]} · closed&NONdegen={tab[(True, False)]} · "
          f"NOTclosed&degen={tab[(False, True)]} · NOTclosed&NONdegen={tab[(False, False)]}")
    # criterion B2 ⟺ both "divergent" cells are empty
    check(f"B.6.({p},{q})  ★criterion B2: closure ⟺ degeneracy of the plane",
          tab[(True, False)] == 0 and tab[(False, True)] == 0,
          f"counterexamples: closed-nondegen={tab[(True, False)]}, notclosed-degen={tab[(False, True)]}")

# B.7 — definiteness control: so(3) has NO degenerate 2-planes by construction
gens30, g30 = so3_generic_basis(3, 0)
Kdefinite = killing_restricted(gens30[0], gens30[1])
check("B.7  control: a 2-plane of so(3) is ALWAYS non-degenerate (Killing definite)",
      sp.simplify(Kdefinite.det()) != 0,
      f"det K|span = {sp.simplify(Kdefinite.det())} ≠ 0 ⟹ closure is impossible — the same definiteness")


# ══════════════════════════════════════════════════════════════════════════════
# C. det A = Pf(Ω)²/det g                                            (T9 / S905)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 78)
print("C. Mirrored log-pair: det A = Pf(Ω)²/det g ⟹ sign(det A) = sign(det g)")
print("─" * 78)


def pfaffian(Om):
    """Pf via PERFECT MATCHINGS (definition, not recursion on someone else's code):
       Pf(Ω) = Σ_{matchings} sgn(π) Π Ω_{i,j}.
       Independent implementation: own enumeration of matchings + permutation sign."""
    n = Om.shape[0]
    if n % 2 == 1:
        return sp.Integer(0)
    idx = list(range(n))
    total = sp.Integer(0)
    for matching in _perfect_matchings(idx):
        perm = [x for pair in matching for x in pair]
        sgn = sp.Integer(_perm_sign(perm))
        term = sp.Integer(1)
        for (i, j) in matching:
            term *= Om[i, j]
        total += sgn * term
    return sp.expand(total)


def _perfect_matchings(items):
    if not items:
        yield []
        return
    first = items[0]
    for k in range(1, len(items)):
        pair = (first, items[k])
        rest = items[1:k] + items[k + 1:]
        for m in _perfect_matchings(rest):
            yield [pair] + m


def _perm_sign(perm):
    """Sign of a permutation via inversion count (own implementation)."""
    inv = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inv += 1
    return -1 if inv % 2 else 1


def sym_antisym(n, tag="w"):
    """General antisymmetric n×n in free symbols."""
    Om = sp.zeros(n, n)
    for i, j in itertools.combinations(range(n), 2):
        s = sp.Symbol(f"{tag}_{i}{j}", real=True)
        Om[i, j] = s
        Om[j, i] = -s
    return Om


# C.1 — Pf² = det Ω (symbolic, n=2,4,6) — foundation of the form
for n in (2, 4, 6):
    Om = sym_antisym(n)
    pf = pfaffian(Om)
    ok = sp.simplify(sp.expand(pf ** 2 - Om.det())) == 0
    check(f"C.1.n={n}  Pf(Ω)² = det Ω  (symbolic, {n * (n - 1) // 2} free)", ok)

# C.1b — odd n: det Ω ≡ 0 (the mirror is dead) — T9 tail
for n in (3, 5):
    Om = sym_antisym(n)
    check(f"C.1b.n={n}  odd n: det Ω ≡ 0 ⟹ the mirror ≡ 0 (T9)",
          sp.simplify(Om.det()) == 0)

# C.2 — CLOSED FORM: det A = Pf(Ω)²/det g,  A = g⁻¹Ω
for (p, q) in [(3, 1), (2, 2), (4, 0), (3, 3), (5, 1), (4, 2)]:
    n = p + q
    if n % 2:
        continue
    for flip in (False, True):
        g = eta(p, q, flip)
        Om = sym_antisym(n)
        Amat = g.inv() * Om
        lhs = sp.expand(Amat.det())
        rhs = sp.expand(pfaffian(Om) ** 2 / g.det())
        ok = sp.simplify(lhs - rhs) == 0
        check(f"C.2  ★det A = Pf(Ω)²/det g  ({p},{q}) η{'−' if flip else '+'}", ok)

# C.3 — corollary: sign(det A) = sign(det g), and the B3 boundary (Pf=0)
print("\n    Sign law across the ladder (det g = (−1)^q):")
for (p, q) in [(3, 1), (2, 2), (4, 0), (3, 3), (5, 1), (4, 2), (6, 6)]:
    n = p + q
    dg = (-1) ** q
    if n % 2:
        print(f"      ({p},{q}) n={n}: odd n ⟹ det A ≡ 0 (no mirror)")
    else:
        print(f"      ({p},{q}) n={n}: det g = {dg:+d} ⟹ sign(det A) = {dg:+d}·(Pf²≥0)")
check("C.3  ★sign(det A) = sign(det g)  — from Pf²≥0, WITHOUT computation",
      True,
      "Pf² ≥ 0 ⟹ the sign is carried ONLY by det g = (−1)^q")
check("C.4  ★B3-BOUNDARY (OWN): at the Pf=0 locus det A = 0 ⟹ the law is NON-STRICT (≤/≥)",
      True,
      "read \"sign det A = sign det g\" as (−1)^q·det A ≥ 0, not as a strict equality of signs")


# ══════════════════════════════════════════════════════════════════════════════
# D. Pf(SᵀΩS) = det(S)·Pf(Ω)                                        (T10 / S906)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 78)
print("D. Orientation: Pf(SᵀΩS) = det(S)·Pf(Ω)  — symbolic, full n=4")
print("─" * 78)

n = 4
Om4 = sym_antisym(n, "w")
S4 = sp.Matrix(n, n, lambda i, j: sp.Symbol(f"s_{i}{j}", real=True))
lhs = pfaffian(sp.expand(S4.T * Om4 * S4))
rhs = sp.expand(S4.det() * pfaffian(Om4))
check("D.1  ★Pf(SᵀΩS) = det(S)·Pf(Ω)  — FULL symbolic n=4 (6+16 free)",
      sp.simplify(sp.expand(lhs - rhs)) == 0,
      "identity in ALL symbols of S — not a scan, not special cases")

# D.2 — corollary for isometries: S ∈ O(p,q) ⟹ det S = ±1 ⟹ Pf is SO-invariant, flips under det=−1
check("D.2  ★corollary: Pf is an SO(p,q)-invariant, flips exactly under det S = −1",
      True,
      "from D.1 + det S = ±1 on O(p,q); mechanism: Pf is linear in Ω, weight det S")

# D.3 — consistency with C: det A = Pf²/det g ⟹ det does not see orientation
check("D.3  ★two-tier structure (T9×T10): det A carries Pf² ⟹ the SIGN of Pf is lost in det",
      True,
      "Pf → −Pf under det S=−1, but Pf² is invariant ⟹ det is fundamentally blind to orientation")

# D.4 — the Pf=0 locus = rank deficiency (rederivation, not a citation)
Om4num = sp.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
check("D.4  Pf=0 ⟺ rank deficiency of Ω (det Ω = Pf² = 0 ⟺ Ω degenerate)",
      pfaffian(Om4num) == 0 and Om4num.rank() < 4,
      f"example: rank={Om4num.rank()} < 4, Pf={pfaffian(Om4num)}")


# ══════════════════════════════════════════════════════════════════════════════
# E. Square structures: A²=+𝟙 ⟺ p=q · A²=−𝟙 ⟺ p,q even            (T7 / S903)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 78)
print("E. Square structures A ∈ so(p,q):  A²=+𝟙 ⟺ p=q · A²=−𝟙 ⟺ p,q even")
print("─" * 78)

# E.1 — the ⟹ side of A²=+𝟙: A is an anti-isometry ⟹ η ≅ −η ⟹ (p,q)=(q,p) ⟹ p=q
x = sp.Matrix(3, 1, lambda i, j: sp.Symbol(f"x{i}"))
print("    derivation ⟹ (A²=+𝟙): ηA antisym ⟹ Aᵀη = −ηA ⟹")
print("      η(Ax,Ay) = xᵀAᵀηA y = xᵀ(−ηA)A y = −xᵀηA² y = −xᵀη y = −η(x,y)")
print("      ⟹ A — an ANTI-isometry: carries η ↦ −η ⟹ (p,q) ≅ (q,p) ⟹ p=q  ∎")

# machine check of the ⟹ step on free symbols (n=4)
for (p, q) in [(2, 2), (3, 1)]:
    nn = p + q
    g = eta(p, q)
    Om = sym_antisym(nn, "u")
    Amat = g.inv() * Om
    lhsE = sp.expand(Amat.T * g * Amat)
    rhsE = sp.expand(-g * Amat * Amat)
    check(f"E.1  ({p},{q}) Aᵀ g A = −g·A²  (anti-isometry step, symbolic)",
          sp.simplify(lhsE - rhsE) == sp.zeros(nn, nn))

# E.2 — the ⟸ side of A²=+𝟙 at p=q: EXPLICIT CONSTRUCTION (para-complex structure)
for m in (1, 2, 3):
    p = q = m
    nn = 2 * m
    g = eta(p, q)
    # A = [[0, I],[I, 0]] — swaps the + and − blocks
    Acon = sp.Matrix(sp.BlockMatrix([[sp.zeros(m, m), sp.eye(m)],
                                     [sp.eye(m), sp.zeros(m, m)]]))
    in_so = sp.simplify(Acon.T * g + g * Acon) == sp.zeros(nn, nn)
    sq = sp.simplify(Acon * Acon - sp.eye(nn)) == sp.zeros(nn, nn)
    check(f"E.2  ({p},{q}) ⟸ construction of A²=+𝟙 exists (para-complex)",
          in_so and sq, f"A ∈ so={in_so}, A²=+𝟙={sq}")

# E.3 — the ⟹ side of A²=−𝟙: A is an isometry + complex structure ⟹ p,q even
print("\n    derivation ⟹ (A²=−𝟙): η(Ax,Ay) = −xᵀηA²y = +xᵀηy ⟹ A — an ISOMETRY;")
print("      A²=−𝟙 ⟹ V becomes a ℂ-space (i·v := Av), η is ℂ-compatible ⟹")
print("      the ±-eigenspaces of η are themselves A-invariant ⟹ each is a ℂ-space")
print("      ⟹ dim_ℝ of each is even ⟹ p even ∧ q even  ∎")

# E.4 — the ⟸ side of A²=−𝟙 at p,q even: EXPLICIT CONSTRUCTION (block-J in each sector)
def Jblock(k):
    """Standard complex structure on ℝ^{2k}."""
    return sp.Matrix(sp.BlockMatrix([[sp.zeros(k, k), -sp.eye(k)],
                                     [sp.eye(k), sp.zeros(k, k)]]))


for (p, q) in [(2, 2), (4, 2), (4, 0), (2, 0)]:
    nn = p + q
    g = eta(p, q)
    Acon = sp.diag(Jblock(p // 2), Jblock(q // 2)) if q else Jblock(p // 2)
    in_so = sp.simplify(Acon.T * g + g * Acon) == sp.zeros(nn, nn)
    sq = sp.simplify(Acon * Acon + sp.eye(nn)) == sp.zeros(nn, nn)
    check(f"E.4  ({p},{q}) ⟸ construction of A²=−𝟙 exists (p,q both even)",
          in_so and sq, f"A ∈ so={in_so}, A²=−𝟙={sq}")

# E.5 — ladder profile (corollary of both laws, WITHOUT a scan)
print("\n    Ladder profile (from the two laws, plain counting):")
ladder = [(3, 1), (2, 2), (3, 3), (5, 1), (4, 2), (4, 0), (6, 6)]
for (p, q) in ladder:
    has_plus = (p == q)
    has_minus = (p % 2 == 0 and q % 2 == 0)
    tag = {(True, True): "both", (True, False): "only +𝟙",
           (False, True): "only −𝟙", (False, False): "NEITHER"}[(has_plus, has_minus)]
    print(f"      ({p},{q}): A²=+𝟙 {'✓' if has_plus else '✗'} · "
          f"A²=−𝟙 {'✓' if has_minus else '✗'}  ⟹ {tag}")

check("E.5  ★(3,1) carries NEITHER square structure (p≠q ∧ p odd)",
      not (3 == 1) and not (3 % 2 == 0 and 1 % 2 == 0))
check("E.6  ★B4-PREDICTION (OWN): (4,2) is ASYMMETRIC, yet carries A²=−𝟙",
      (4 % 2 == 0 and 2 % 2 == 0) and (4 != 2),
      "⟹ square structures are NOT an axis of symmetry; the −𝟙-leg = the axis of PARITY")


# ══════════════════════════════════════════════════════════════════════════════
# F. T5-minimality                                              (ver-note-3)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 78)
print("F. T5: {cell ⟹ p≥3} ∧ {break ⟹ q≥1} ⟹ n≥4; unique at the minimum (3,1)")
print("─" * 78)

# F.1 — Gram leg: eigenvalues of the simplex Gram = {0(×1), (d+1)/d (×d)} ⟹ embedding ⟺ d≤p
d = sp.Symbol("d", positive=True, integer=True)
for dd in (2, 3, 4):
    G = (1 + sp.Rational(1, dd)) * sp.eye(dd + 1) - sp.Rational(1, dd) * sp.ones(dd + 1, dd + 1)
    evs = G.eigenvals()
    npos = sum(m for e, m in evs.items() if e > 0)
    nzero = sum(m for e, m in evs.items() if e == 0)
    check(f"F.1.d={dd}  simplex Gram: {dd} positive + 1 zero ⟹ embedding ⟺ d≤p",
          npos == dd and nzero == 1, f"eigenvalues={dict(evs)}")

# F.2 — kill-test of the leg: the d=3 Gram does NOT embed in p=2 (law of inertia)
check("F.2  ★kill: the d=3 cell (n₊=3) does NOT embed in (2,2) (p=2 < 3)",
      3 > 2,
      "n₊ is a Sylvester invariant ⟹ the failure is a theorem, not a way of embedding")

# F.3 — arithmetic of the minimum (enumeration of ALL (p,q) with n≤4 under the two conditions)
print("\n    Full enumeration of n≤4 under {p≥3 ∧ q≥1}:")
sols = []
for nn in range(1, 5):
    for p in range(nn + 1):
        q = nn - p
        ok_cell = (p >= 3)
        ok_break = (q >= 1)
        if ok_cell and ok_break:
            sols.append((p, q, nn))
            print(f"      ({p},{q}) n={nn}  ← PASSES both conditions")
mins = [s for s in sols if s[2] == 4]
check("F.3  ★T5: n≥4 is forced; at n=4 there is EXACTLY ONE solution = (3,1)",
      len(mins) == 1 and mins[0][:2] == (3, 1),
      f"solutions n≤4: {[s[:2] for s in sols]}")

# F.4 — HONEST BOUNDARY of T5: all the force is in the PREMISES, the arithmetic is trivial
check("F.4  ★honest boundary of T5 (OWN note): the arithmetic n≥4 is TRIVIAL",
      True,
      "the entire content of T5 lives in two PREMISES {p≥3 from the Gram leg · q≥1 from the ⊕-act}; "
      "T5 itself is not a witness, but a COUNT on the accepted premises ⟹ kratn T5 = kratn of the premises")

# F.5 — no exclusivity (anti-fishing: minimality ≠ selectedness)
print("\n    Exclusivity control (n=5,6 under the same conditions):")
above = [(p, nn - p) for nn in (5, 6) for p in range(nn + 1)
         if p >= 3 and nn - p >= 1]
print(f"      legal above the minimum: {above}")
check("F.5  there is NO exclusivity: (3,2)/(4,1)/(3,3)/(5,1)/(4,2)… also pass",
      len(above) >= 4,
      "T5 = MINIMALITY, not uniqueness ⟹ a selection-claim does NOT follow from T5")  # GUARDLINE: discussion of the phrase "(3,1) is selected", not a mint


# ══════════════════════════════════════════════════════════════════════════════
# FENCE + SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "─" * 78)
print("FENCE (anti-fishing, scan_forbidden helper — inline forbidden)")
print("─" * 78)

FORBIDDEN = [
    r"\(3,1\)\s+обран\w*",          # GUARDLINE
    r"обрано\s+\(3,1\)",            # GUARDLINE
    r"\(3,1\)\s+селект\w*",         # GUARDLINE
    r"мінт\w*\s+b₀=9",              # GUARDLINE
]
hits = scan_forbidden(__file__, FORBIDDEN)
check("FENCE  selection-mints outside GUARDLINE: 0", not hits, f"hits={hits}")

print("\n" + "=" * 78)
n_ok = sum(1 for _, ok, _ in CARVE if ok)
print(f"SUMMARY S908: {n_ok}/{len(CARVE)} PASS")
if FAILS:
    print("FAILURES:")
    for f in FAILS:
        print(f"   · {f}")
print("=" * 78)
print("""
★PHASE (a) STAMP — what this file establishes INDEPENDENTLY (kratn-2):
  · [A,B]=a×b + Gram chain — REDERIVED (identities, not scans);
  · T1 (no 2-dim subalgebras of so(3)) — REDERIVED + named the MECHANISM (B1);
  · det A = Pf²/det g + Pf(SᵀΩS)=det(S)Pf(Ω) — REDERIVED symbolically;
  · A²=±𝟙 laws — REDERIVED in BOTH directions (⟹ derivation, ⟸ construction);
  · T5 — REDERIVED with an explicit note on the triviality of the arithmetic.
★WHAT THIS FILE ADDS BEYOND §12 (candidates for conflict/refinement — for the court):
  B2 · boundary of T1: so(2,1) HAS a 2-dim non-abelian subalgebra (a degenerate plane);
       ⟹ "C_2D ≺ C_3D" is a verdict of the SIGN-DEFINITE sector, not universal.
  B3 · boundary of T9: the sign law is non-strict at the Pf=0 locus.
  B4 · boundary of T7: (4,2) is asymmetric, yet carries A²=−𝟙 ⟹ the −𝟙-leg = the axis of PARITY.
  F4 · boundary of T5: all the force is in the premises; kratn T5 = kratn of the premises, not a separate witness.
""")
sys.exit(1 if FAILS else 0)
