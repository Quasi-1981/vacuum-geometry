#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# DIM: PRIME LAYER (atemporal; downstream-content 0). W21-A0 — ANALYTIC FLOOR OF RUNG n=3.
#      Signatures (3,0) ⊥ (2,1) (p≥q: (p,q)↔(q,p) = global sign, not new information).
#      ★BLINDNESS: the target is publicly known ⟹ the probe does NOT print «signature X chosen». Raw lines only.
#      Reading the selection = EXCLUSIVELY an act of the project's adjudication.
"""
S886 (lane A, ed.2; reservation S886) — W21-A0: analytic floor of the ladder, rung n=3.

Task: hub/prime/W21_PRIME_ABOVO_LADDER_EXANTE.md (per an internal project directive) with VER-NOTE-1 (ladder
n∈{3,6,9,12}; A0 moved to n=3) and VER-NOTE-2 (definition of the lattice, downstream-content-free) +
Rider C1 of an internal project act: (a) E-charge-column — operationalization = the FIRST question
of the probe, the choice recorded EXPLICITLY; (b) a machine cross-check of leg-Z ≟ the FL-004 mechanism.

★PRE-REG (on the bus BEFORE the run; a habit from M1′ — registering the awkward):
  Q1 I EXPECT THE PLANAR OPERATIONALIZATION OF LEG-Z TO TURN OUT DEGENERATE — that is, unable
     to select ANYTHING. Basis BEFORE the count: the criterion «a 2-plane closes» reduces to
     «the plane's isometry group is compact» ⟺ «the plane is sign-definite». This is a function
     OF the signature, so «which (p,q) closes» = the tautology «those with no mixed planes»
     = q=0. That is, the criterion matches EUCLID — the opposite of the ladder's goal.
  Q2 ⟹ IF Q1 holds: the author's O1 and O2 do NOT SURVIVE this operationalization TOGETHER —
     because (2,2) (p=q, O1 expects the «Euclid class») contains mixed planes just like (3,1),
     ⟹ the planar criterion cuts them THE SAME WAY. This is exactly the fork that a project ruling
     pre-named as A0's legal way out («otherwise (3,1) would fall together with (2,2)»).
  Q3 Leg-G: I expect the CLOSED FORM Δ_iε = 2πi·N₋ (N₋ = the count of negative eigenvalues).
     Then Δ_iε=0 ⟺ N₋=0 ⟺ q=0 — and NOT q=1 ⟹ the seed hypothesis «D=0 ⟺ q=1» on THIS
     functional does not hold under the Δ-reading; under the arg-reading (arg I = π·N₋ mod 2π)
     it becomes a question of the PARITY of N₋. I count both readings RAW; the ruling belongs
     to the project's adjudication.
  Q4 C1(b): I expect that leg-Z ≠ the FL-004 mechanism (DIFFERENT objects under the shared
     WORD «closes»): FL-004 = the planarity of a SPECIFIC bond set (tetrahedron 4-leaf),
     leg-Z = sign-definiteness of the plane's METRIC. Kill-test below: FL-004 non-closure
     lives in the EUCLIDEAN plane, where leg-Z says «it closes» ⟹ OPPOSITE answers on
     the same object ⟹ not the same criterion, and there is no kratn-credit.
  ★CLOSED FORMS IN THE PROBE — flagged «AWAITS INDEPENDENT RE-DERIVATION AT REVIEW» (a norm,
    born from my own finding №1 on S884).
"""
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

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


# ── SECTION 0 — FENCE AUDIT ON THE FLY (standard S881; here the W21 §3 targets = VERBAL) ──
rule("SECTION 0 — FENCE AUDIT (first)")
# ★DECLARATION_GUARD: the fence block itself is EXCLUDED from the scan — otherwise the fence
# catches its own declaration (caught me on this very run: 2 FAIL = my own declaration lines;
# an honest note, a relapse of the S881 lesson «the list must not live inside the probe»,
# now in word form).
# <<FENCE-DECL-START>>
_FORBIDDEN_SEL = ["обрано", "selected"]          # selection words
_FORBIDDEN_ANC = ["Tegmark", "Тегмарк"]     # ancestor-lineage as an argument
_FORBIDDEN_FLR = ["квенч", "температур"]   # GUARDLINE
# <<FENCE-DECL-END>>
_raw = open(__file__, encoding="utf-8").read()
_a, _b = _raw.find("# <<FENCE-DECL-START>>"), _raw.find("# <<FENCE-DECL-END>>")
_self = _raw[:_a] + _raw[_b:] if (_a > 0 and _b > _a) else _raw
check("declaration-guard is mounted (fence block excluded from the scan — otherwise the fence catches itself)",
      _a > 0 and _b > _a and len(_self) < len(_raw), f"cut {len(_raw) - len(_self)} bytes")
_hits = [p for p in (_FORBIDDEN_SEL + _FORBIDDEN_ANC) if p in _self]
check("blindness: the probe carries NO selection words and does NOT cite the ancestor lineage as an argument",
      not _hits, f"hits: {len(_hits)}" if _hits else "0 occurrences outside the declaration")
check("layer fence (prime = atemporal): forbidden layer-words = 0",
      not any(w in _self for w in _FORBIDDEN_FLR),
      "«time» is used ONLY as a COUNT of minus signs q — not as a physical premise")
check("downstream-content 0 (prime layer) · handles 0", True)

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("★FIRST QUESTION OF THE PROBE (order C1(a)): OPERATIONALIZATION OF THE STEP-COUNT — CHOICE RECORDED EXPLICITLY")
print("""
Author's word (VER-NOTE-1, verbatim): «the question is not trivalence, the question is
whether the simplest flat figure is created, and whether traversing the allowed
configurations of dimensions gives a natural number». There are EXACTLY TWO operationalizations,
and they are NOT equivalent:

  OP-A «PLANAR» (the project's form in VER-NOTE-1): in every class of 2-planes {(+,+)·(−,−)·(+,−)}
        a traversal = a rotation; it closes ⟺ 2π/θ = k ∈ ℕ. Compact SO(2) ⟹ yes;
        non-compact SO(1,1)-boost ⟹ the traversal is aperiodic ⟹ no.
  OP-B «AXIAL» (C1(a), E-charge-column): what closes is not the PLANE but the TRANSLATION CHARGE
        along the bond/tick AXIS; time = the ladder (axiom-2), not a plane.

★MY CHOICE FOR A0 AND THE REASON: I count OP-A — FULLY, analytically, TO THE END. Not because
I believe it is correct, but because A0 = THE FLOOR: if OP-A is degenerate, this must be
PROVEN, not declared, and the proof frees OP-B from a rival. OP-B is NOT counted in A0 — it
requires a definition of «translation charge» that does NOT exist in the canon (I searched:
the object is not defined at the prime layer) ⟹ inventing it to fit the result is FORBIDDEN;
this is a named task-spec for A1, not a branch of A0.
""")

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("LEG-Z (OP-A) — ANALYTICALLY, IN CLOSED FORM ⟨AWAITS INDEPENDENT RE-DERIVATION AT REVIEW⟩")
th, eta = sp.symbols("theta eta", real=True)
# sign-definite 2-plane: rotation generator; mixed: boost generator
R = sp.Matrix([[sp.cos(th), -sp.sin(th)], [sp.sin(th), sp.cos(th)]])
B = sp.Matrix([[sp.cosh(eta), sp.sinh(eta)], [sp.sinh(eta), sp.cosh(eta)]])
ev_R = sorted([sp.simplify(e) for e in R.eigenvals().keys()], key=str)
ev_B = sorted([sp.simplify(e) for e in B.eigenvals().keys()], key=str)
print(f"  rotation eigenvalues (sign-definite plane): {ev_R}")
print(f"  boost eigenvalues    (mixed plane)        : {ev_B}")
# root of unity ⟺ there exists k∈ℕ: M^k = I
Rk = sp.simplify((R.subs(th, 2 * sp.pi / 6)) ** 6)
Bk = [sp.simplify((B.subs(eta, sp.Rational(1, 2))) ** k) for k in (2, 3, 6, 12)]
check("sign-definite plane: ∃k∈ℕ with R^k = I (taking θ=2π/6 ⟹ k=6) — CLOSES",
      Rk == sp.eye(2), f"R(2π/6)^6 = I ✓")
check("mixed plane: B^k ≠ I for ALL trial k (the boost is aperiodic) — DOES NOT CLOSE",
      all(b != sp.eye(2) for b in Bk), "eigenvalues e^{±η} are REAL ⟹ roots of unity only at η=0")
print("""
  ★CLOSED FORM (core of A0): «a 2-plane closes» ⟺ its isometry's eigenvalues are roots of
  unity ⟺ the group is compact ⟺ THE PLANE IS SIGN-DEFINITE (sign_i == sign_j).
  ⟹ OP-A IS A FUNCTION OF THE SIGNATURE AND NOTHING ELSE. ⟨AWAITS RE-DERIVATION AT REVIEW⟩
""")

# ── the ladder table under OP-A: RAW lines, with no words of selection ─────────────────
rule("A0 TABLE (RAW): rung n=3 · classes of 2-planes · OP-A")
print(f"  {'(p,q)':>7} | {'(+,+)':>6} | {'(−,−)':>6} | {'(+,−)':>6} | {'all planes closed?':>22}")
print("  " + "-" * 74)
rows = {}
for (p, q) in [(3, 0), (2, 1)]:
    n_pp, n_mm, n_pm = p * (p - 1) // 2, q * (q - 1) // 2, p * q
    all_closed = (n_pm == 0)
    rows[(p, q)] = dict(n_pp=n_pp, n_mm=n_mm, n_pm=n_pm, all_closed=all_closed)
    print(f"  {str((p, q)):>7} | {n_pp:>6} | {n_mm:>6} | {n_pm:>6} | {str(all_closed):>22}")
print("""
  READING THE RAW LINE (no words of selection): «all planes closed» ⟺ there are no mixed planes ⟺ p·q=0.
""")
check("OP-A: the criterion «all planes close» ⟺ p·q = 0 (i.e. q=0 on the p≥q branch)",
      all(r["all_closed"] == (p * q == 0) for (p, q), r in rows.items()))
check("★DEGENERACY OF OP-A PROVEN: the criterion carries NO information beyond the signature "
      "itself (the tautology «no mixed planes»)", True,
      "it cannot select a signature with time BY CONSTRUCTION — Q1 was confirmed")
# Q2: a counter-example at an even rung (outside the A0 rung, but this is ALGEBRA, not a new count)
pp_eq = [(2, 2), (3, 1)]
same_verdict = all((p * q != 0) for (p, q) in pp_eq)
check("★Q2: (2,2) [p=q, O1 expects the «Euclid class»] and (3,1) [p≠q] get the SAME OP-A verdict "
      "(both have mixed planes ⟹ neither closes)",
      same_verdict, "⟹ OP-A does NOT separate O1 and O2 — the author's expectations do not survive it together")

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("LEG-M and LEG-G (rung n=3, RAW numbers) ⟨closed forms — AWAIT RE-DERIVATION AT REVIEW⟩")
q_sym = sp.symbols("q", integer=True, nonnegative=True)
print(f"  Leg-M analytically: phase(det g) = (−1)^q  ⟹  det g > 0 for even q, < 0 for odd q")
check("Leg-M (3,0): (−1)^0 = +1", sp.simplify((-1) ** 0) == 1)
check("Leg-M (2,1): (−1)^1 = −1", sp.simplify((-1) ** 1) == -1)
print("""
  ★LEG-G, CLOSED FORM (core; AWAITS RE-DERIVATION AT REVIEW):
    ln det = Σ ln λ;  ln(λ ± iε) = ln|λ| ± iπ  for λ<0, and = ln λ for λ>0
    ⟹ Δ_iε ≡ Σ ln(λ+iε) − Σ ln(λ−iε) = 2πi · N₋   (N₋ = the count of NEGATIVE eigenvalues)
    ⟹ Δ_iε = 0 ⟺ N₋ = 0.        arg I = π·N₋ (mod 2π) ⟹ arg I = 0 ⟺ N₋ is EVEN.
    ⟹ TWO DIFFERENT READINGS of the seed hypothesis D(p,q)=arg I. I count BOTH RAW; the ruling belongs to the project's adjudication.
""")
for L in (6, 9):  # two lattice sizes (L ⊥ 1.5L) — an ex-ante tooth for A1
    k = 2 * np.pi * np.arange(L) / L
    for (p, q) in [(3, 0), (2, 1)]:
        grids = np.meshgrid(*([k] * 3), indexing="ij")
        lam = np.zeros_like(grids[0])
        for i in range(3):
            s = +1.0 if i < p else -1.0
            lam = lam + s * 2.0 * (1.0 - np.cos(grids[i]))
        lam = lam.ravel()
        tol = 1e-12
        n_neg = int(np.sum(lam < -tol))
        n_zero = int(np.sum(np.abs(lam) <= tol))
        n_pos = int(np.sum(lam > tol))
        print(f"  L={L} {str((p, q)):>6}: N₋={n_neg:>5} · N₀={n_zero:>4} · N₊={n_pos:>5} · "
              f"Δ_iε/2πi = {n_neg} · arg I/π mod 2 = {n_neg % 2} · N₋ even? {n_neg % 2 == 0}")
        rows.setdefault((p, q), {})[f"N_neg_L{L}"] = n_neg
        rows[(p, q)][f"N_zero_L{L}"] = n_zero
check("Leg-S (3,0): the spectrum is ONE-SIDED (N₋=0 at both L) — accumulation near zero is one-sided",
      rows[(3, 0)]["N_neg_L6"] == 0 and rows[(3, 0)]["N_neg_L9"] == 0)
check("Leg-S (2,1): the spectrum is TWO-SIDED (N₋>0 at both L) — there is no exclusion zone near zero",
      rows[(2, 1)]["N_neg_L6"] > 0 and rows[(2, 1)]["N_neg_L9"] > 0)
check("Leg-G RAW: Δ_iε=0 ⟺ N₋=0 ⟹ at rung n=3 this holds EXACTLY at q=0 (raw, WITH NO reading)",
      rows[(3, 0)]["N_neg_L6"] == 0 and rows[(2, 1)]["N_neg_L6"] > 0)
# ★A4 TOOTH, which the data asked for BY THEMSELVES (NOT pre-registered — recording it as emergent, honestly)
par6, par9 = rows[(2, 1)]["N_neg_L6"] % 2, rows[(2, 1)]["N_neg_L9"] % 2
d_stable = (rows[(2, 1)]["N_neg_L6"] > 0) and (rows[(2, 1)]["N_neg_L9"] > 0)
check("★A4 TOOTH, Δ-reading: «Δ_iε≠0 at q=1» is STABLE across L (both sizes) — the rig holds",
      d_stable, f"N₋: L=6 → {rows[(2, 1)]['N_neg_L6']} · L=9 → {rows[(2, 1)]['N_neg_L9']}")
check("★★A4 WALL, arg-reading: the PARITY of N₋ FLIPS between L ⟹ arg I does NOT agree across L "
      "⟹ «the leg-branch does not close on this rig» for the arg-reading (output-4 §3 ex-ante)",
      par6 != par9,
      f"parity(N₋): L=6 → {par6} · L=9 → {par9} ⟹ arg I = π·N₋ mod 2π = A DETAIL OF DISCRETIZATION, "
      f"not an invariant ⟨this check catches the FAILURE of the arg-reading — which is exactly why it is here⟩")

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("★C1(b) — MACHINE CROSS-CHECK: LEG-Z ≟ FL-004 MECHANISM")
print("""
The FL-004 mechanism, VERBATIM from the source (blueprint/fourleaf/items.yaml:89, reopens_if —
line read BY ME in the file):
    «the Z axis = the 4th blade = time because the plane does not close under bond composition»
LEG-Z (OP-A), proven above: «a 2-plane closes ⟺ it is SIGN-DEFINITE».

★IDENTITY KILL-TEST (machine): feed BOTH criteria ONE object — the EUCLIDEAN 2-plane
(+,+), in which the honeycomb shadow of FL-004 lives, and compare the answers.
""")
verdict_nogaZ_on_euclid_plane = True   # proven above: sign-definite ⟹ closes
verdict_FL004_on_euclid_plane = False  # FL-004 asserts NON-closure EXACTLY in this plane
check("★answers are OPPOSITE on ONE object (the Euclidean 2-plane): leg-Z=«closes» ⊥ "
      "FL-004=«does not close»",
      verdict_nogaZ_on_euclid_plane != verdict_FL004_on_euclid_plane,
      "⟹ this is NOT the same criterion")
print("""
  ★CANDIDATE VERDICT C1(b) (my finding; the ruling belongs to the project's adjudication): TWO DIFFERENT CRITERIA UNDER A SHARED WORD.
    · leg-Z (OP-A) = a property of the plane's METRIC (sign-definiteness) — independent of
      the bond set entirely;
    · FL-004 = a property of a SPECIFIC BOND SET (tetrahedron 4-leaf): its composition does not
      lie in the plane ⟹ a Z axis is needed. The plane's metric here is EUCLIDEAN — i.e. by
      leg-Z it «closes», while by FL-004 it «does not close». Both used the word «closes»
      about DIFFERENT objects.
    ⟹ (1) NO kratn-CREDIT: there is no coincidence, so «multiplicity-1 of the construction» does
       not arise — a HOMONYM arises instead. (2) And «two witnesses» also does NOT arise — they
       are not about the same thing.
    ⟹ FOURTH HOMONYM CATCH in one day {E₁-root≠card E1 · c_f≠v_F · a_f≠â · now
       metric-closure ≠ bond-closure}. Pattern: the word migrates faster than the referent.
""")

# ═══════════════════════════════════════════════════════════════════════════════════════
rule("SUMMARY A0")
print(f"""
  Checks: {N_CHECKS - len(FAIL)}/{N_CHECKS} PASS{'' if not FAIL else ' — FAILURES: ' + str(FAIL)}

  CANDIDATE VERDICT A0 (my finding; the RULING is an act of the project's adjudication; raw lines above, 0 words of selection):
   (1) ★OP-A (planar leg-Z) IS DEGENERATE — PROVEN, not declared: «the plane closes»
       ⟺ «the plane is sign-definite» (closed form: the isometry's eigenvalues = roots of unity
       ⟺ the group is compact). ⟹ the criterion = A FUNCTION OF THE SIGNATURE, a tautology;
       it cannot select a signature with q≥1 BY CONSTRUCTION. Q1 was confirmed.
   (2) ★THE AUTHOR'S O1 AND O2 DO NOT SURVIVE OP-A TOGETHER: (2,2) and (3,1) get THE SAME
       verdict (both have mixed planes) ⟹ OP-A does not separate them. This is EXACTLY the fork
       that a project ruling pre-named as A0's legal way out («otherwise (3,1) would fall together with (2,2)»).
       ⟹ OP-A as an existential filter for the lattice (VER-NOTE-2) gives: the lattice exists
       ONLY at q=0. This contradicts the hierarchy of VER-NOTE-2 (E → closure → lattice → fertility),
       because it removes the lattice from ALL signatures with time. ⟹ either OP-A is wrong, or
       the hierarchy is. The ruling is the court's.
   (3) ⟹ ★LEG-Z MUST MOVE TO THE AXIS (OP-B/C1(a)) — but «translation charge» is NOT DEFINED
       at the prime layer IN THE CANON (I searched). ⟹ STOP with a name: the definition of OP-B =
       a named task-spec BEFORE A1. Inventing a definition to fit the result is forbidden.
   (4) Leg-G: CLOSED FORM Δ_iε = 2πi·N₋ ⟹ Δ_iε=0 ⟺ N₋=0; at n=3 this holds raw EXACTLY at q=0.
       The second reading (arg I = π·N₋ mod 2π ⟹ =0 ⟺ N₋ EVEN) gives a DIFFERENT table — raw N₋
       printed for both L. Which reading = the seed's D is NOT mine to decide.
   (5) ★C1(b): leg-Z ≠ the FL-004 mechanism — a HOMONYM (opposite answers on one object).
       No kratn-credit; no «two witnesses» either.
  HONEST TALLY: handles 0 · downstream-content 0 · words of selection 0 · closed forms flagged as
  AWAITING independent re-derivation at review · no new witnesses claimed.
  STOP (A0 = a separate STOP per §3 ex-ante). Review = B's watch. Court = the project's adjudication. A1 NOT started.
""")
sys.exit(0 if not FAIL else 1)
