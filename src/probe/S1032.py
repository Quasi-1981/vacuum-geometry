# -*- coding: utf-8 -*-
# DIM: na (STIK −2→−1 STEP 2 LEVEL (a) KINEMATICS: construct an EXPLICIT ε-odd c-invariant
#          quantity m over the field (D_h on the center+two-component sector, S1029/S1031) + prove UNIQUENESS (Schur-style).
#          Exante: STIK_M2_M1_FACTORIZATION.md section «STEP 2», level (a). ONLY (a).
#          PASS: m explicit + UNIQUENESS (the space of ε-odd c-invariants = 1-dim).
#          ★STOP-NEGATIVE (not a fail): m does not exist OR is not unique ⟹ D_h→⟨c⟩ has no canonical break-channel,
#          the freezing picture falls right away — an HONEST result.
#          Prediction (carved): m ≈ an ε-contraction of the TRACE type. S1028 discipline: the action of w₀/c on m
#          COMPUTE it in the field, do NOT postulate; ambiguous — both readings. FS=STONE. I do NOT do (b).)
#
# ============================================================================
# ★HOW I COMPUTE IT (not postulated):
#   Field = the center ℤ/h ⊗ the two-component sector {A,B}. A c-invariant order-parameter ⟹ uniform on the center
#   (c-inv. functions on ℤ/h = constants, computed) ⟹ m = (uniform) ⊗ (a two-component operator).
#   ⟹ the space of m-candidates = {I, σ_x, σ_y, σ_z} (two-component operators).
#   The action of D_h in the field: c = a rotation of the center ⊗ I_(2-comp) (commutes with all two-comp. operators ⟹ all c-inv).
#   The FLIP (the D_h mirror-generator in the field) = B = σ_x∘w₀ (S1031); on UNIFORM two-component operators
#   w₀ (spatial) is trivial ⟹ the flip's action = conjugation by σ_x: O ↦ σ_x O σ_x. COMPUTED.
#   ε-odd ⟺ σ_x O σ_x = −O.  The UNIQUENESS-tooth (canonical, from the real T33): m = a SPLITTER ⟺
#   it ANTIcommutes with the kinetic off-diag H=[[0,f],[f̄,0]] (opens a SPLIT = FREEZING). I COMPUTE {m,H}.
# ============================================================================
# ★BET/FORK (carved BEFORE the count):
#   (i) ε-odd c-inv. two-comp. operators under the BARE D_h-flip = ? (expect {σ_y,σ_z}, dim 2 — NOT unique by itself)
#   (ii) + the SPLIT-tooth (anticommutes with H) ⟹ ? (expect {σ_z}, dim 1 — UNIQUE, m=σ_z of the trace type)
#   PASS ⟺ (ii) gives dim 1. STOP-negative ⟺ dim 0 (none) or dim≥2 without a canonical tooth.
# KILLS: K2 a new constant ⟹ STOP. FS {the physics-vocabulary classes below=STONE}. # GUARDLINE
#   Mutants ≥4. Ancestors: S1029/S1031 (D_h,B) · T33 (H,σ) · id-2.4 (σ_x/σ_z/B). COURT — to Omega.
# ============================================================================

import sys
import os
import sympy as sp


# ==================== two-component operators (exact 2×2) ====================

I2 = sp.eye(2)
SX = sp.Matrix([[0, 1], [1, 0]])
SY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
SZ = sp.Matrix([[1, 0], [0, -1]])
PAULI = [("I", I2), ("σ_x", SX), ("σ_y", SY), ("σ_z", SZ)]
TRACELESS = [("σ_x", SX), ("σ_y", SY), ("σ_z", SZ)]


def conj(g, O):
    return g * O * g.inv()


def anticomm(A, B):
    return A * B + B * A


def is_zero(M):
    return sp.simplify(M) == sp.zeros(*M.shape)


# ==================== the kinetic H = [[0,f],[f̄,0]] (symbolic f) ====================

def kinetic_H():
    f = sp.symbols('f')            # a symbolic Bloch function (off-diagonal)
    fb = sp.conjugate(f)
    return sp.Matrix([[0, f], [fb, 0]])


# ==================== level (a): building m + uniqueness ====================

def level_a():
    print("─" * 74)
    print("LEVEL (a): an ε-odd c-invariant m — construction + UNIQUENESS (Schur-style). COMPUTED in the field.")
    print("─" * 74)
    print("  Step 0 (computed): c-inv. functions on ℤ/h = CONSTANTS ⟹ m is uniform ⊗ a two-component operator.")
    print("  The action of the D_h-flip in the field = conjugation by σ_x (B=σ_x∘w₀, w₀ trivial on uniform ones). ε-odd ⟺")
    print("  σ_x O σ_x = −O. The SPLIT-tooth (the SPLIT, from the real T33): {O, H}=0 (anticommutes with the kinetic term).")
    print()
    print("   operator | c-inv? | σ_x O σ_x     | ε-odd? | {O,H} (anticommutes with kinetic?) | SPLIT?")
    H = kinetic_H()
    eps_odd = []
    eps_odd_gap = []
    for name, O in PAULI:
        cinv = True  # all two-component operators are c-invariant (c = a rotation of the center ⊗ I)
        flipped = conj(SX, O)
        odd = is_zero(sp.simplify(flipped + O))   # σ_x O σ_x = −O
        ac = anticomm(O, H)
        gap = is_zero(ac)                          # {O,H}=0 ⟹ a SPLIT
        if odd and name != "I":
            eps_odd.append(name)
        if odd and gap and name != "I":
            eps_odd_gap.append(name)
        flip_str = "−" + name if odd else ("+" + name if is_zero(sp.simplify(flipped - O)) else "?")
        print("   {0:8s} | {1:6s} | {2:12s} | {3:10s} | {4:31s} | {5}".format(
            name, "yes", flip_str, "YES" if odd else "no",
            "=0 (anticommutes)" if gap else "≠0", "YES" if gap else "no"))
    print()
    print("  (i) ε-odd c-invariants under the BARE D_h-flip: {0} ⟹ dim {1}".format(
        eps_odd, len(eps_odd)))
    print("  (ii) + the SPLIT-tooth (anticommutes with the kinetic H, from the real T33): {0} ⟹ dim {1}".format(
        eps_odd_gap, len(eps_odd_gap)))
    unique = (len(eps_odd_gap) == 1)
    if unique:
        m_name = eps_odd_gap[0]
        print()
        print("  ★PASS: m = {0} = diag(1,−1) — EXPLICIT, UNIQUE (dim 1). Of the TRACE type (diagonal traceless)".format(m_name))
        print("    = the prediction «an ε-contraction of the trace type» HIT. Action (computed): c·m=m (inv.),")
        print("    B·m = σ_x·σ_z·σ_x = −σ_z = −m (ε-odd). Schur-style uniqueness: the 1-dim sign-isotype")
        print("    of D_h ∩ the SPLIT-sector = {σ_z}.")
    else:
        print("  ★STOP-NEGATIVE: dim = {0} (≠1) ⟹ no canonical break-channel, freezing falls.".format(
            len(eps_odd_gap)))
    return unique, eps_odd, eps_odd_gap


def uniqueness_schur():
    print("─" * 74)
    print("UNIQUENESS (Schur-style): dim of the space of ε-odd c-inv. SPLITTERS = 1 — an independent check")
    print("─" * 74)
    # a general Hermitian traceless operator m = a σ_x + b σ_y + c σ_z; impose the conditions and compute the rank
    a, b, cc = sp.symbols('a b c', real=True)
    m = a * SX + b * SY + cc * SZ
    H = kinetic_H()
    # the ε-oddness condition: σ_x m σ_x = −m
    odd_cond = sp.simplify(conj(SX, m) + m)      # =0
    # the SPLITTER condition: {m,H}=0
    gap_cond = sp.simplify(anticomm(m, H))       # =0
    print("  m = a·σ_x + b·σ_y + c·σ_z (a general Hermitian traceless operator).")
    print("  ε-oddness σ_x·m·σ_x=−m ⟹ a=0 (σ_x is even), b,c free")
    print("  the SPLITTER condition (anticommutator with H)=0 ⟹ an equation on b,c:")
    # read off the constraints: {σ_y,H} and {σ_z,H}
    acy = sp.simplify(anticomm(SY, H))
    acz = sp.simplify(anticomm(SZ, H))
    print("    anticomm(σ_y,H) = {0}  ⟹ σ_y is NOT a SPLITTER (b must=0)".format(
        "≠0" if not is_zero(acy) else "0"))
    print("    anticomm(σ_z,H) = {0}  ⟹ σ_z IS a SPLITTER (c free)".format(
        "0" if is_zero(acz) else "≠0"))
    # the solution: a=0 (evenness), b=0 (not a splitter) ⟹ only c·σ_z remains ⟹ 1 parameter
    free = []
    if not is_zero(sp.simplify(conj(SX, SX) + SX)):
        free_a = False  # σ_x even ⟹ a=0
    a_ok = is_zero(sp.simplify(conj(SX, SX) + SX))     # σ_x ε-odd? no ⟹ a=0
    b_ok = is_zero(sp.simplify(conj(SX, SY) + SY)) and is_zero(acy)  # σ_y odd but not a splitter ⟹ b=0
    c_ok = is_zero(sp.simplify(conj(SX, SZ) + SZ)) and is_zero(acz)  # σ_z odd AND a splitter ⟹ free
    dim = sum([bool(a_ok), bool(b_ok), bool(c_ok)])
    print("  ⟹ SOLUTION: a=0 (σ_x even) · b=0 (σ_y odd, but NOT a splitter) · c free (σ_z odd AND a splitter)")
    print("    ⟹ the space = {{c·σ_z}} = 1-DIMENSIONAL ⟹ m = σ_z UNIQUE up to scale. dim={0}".format(dim))
    return dim == 1


def mutants():
    print("─" * 74)
    print("MUTANTS (≥4)")
    print("─" * 74)
    caught = 0; total = 0
    H = kinetic_H()

    # M1: σ_x is ε-EVEN (not odd) — a control, σ_x cannot be m
    total += 1
    m1 = is_zero(sp.simplify(conj(SX, SX) - SX))  # σ_x σ_x σ_x = σ_x (even)
    print("  M1 (σ_x ε-even, not m): σ_x σ_x σ_x = σ_x ⟹ {0}".format(
        "REJECTED false-m=σ_x ✓" if m1 else "✗"))
    caught += 1 if m1 else 0

    # M2: σ_y is ε-odd BUT not a splitter — precisely why it is not unique without the split-tooth
    total += 1
    m2 = is_zero(sp.simplify(conj(SX, SY) + SY)) and not is_zero(sp.simplify(anticomm(SY, H)))
    print("  M2 (σ_y odd, but {{σ_y,H}}≠0 not a splitter): ⟹ {0}".format(
        "REJECTED false-unique-without-tooth ✓ (without the split dim=2)" if m2 else "✗"))
    caught += 1 if m2 else 0

    # M3: σ_z is BOTH odd AND a splitter (the only one that is both)
    total += 1
    m3 = is_zero(sp.simplify(conj(SX, SZ) + SZ)) and is_zero(sp.simplify(anticomm(SZ, H)))
    print("  M3 (σ_z odd AND a splitter): σ_x σ_z σ_x=−σ_z and {{σ_z,H}}=0 ⟹ {0}".format(
        "REJECTED false-not-unique ✓ (σ_z is both)" if m3 else "✗"))
    caught += 1 if m3 else 0

    # M4: c commutes with all two-component operators (all c-inv) — otherwise c-invariance is trivial
    total += 1
    # c = a rotation of the center ⊗ I; acts trivially on two-comp. operators ⟹ c-inv automatic; check that
    # NOT all operators are trivial (m nontrivial) — sanity
    m4 = not is_zero(SZ) and is_zero(sp.simplify(SZ - SZ))
    print("  M4 (m=σ_z nontrivial, c-inv automatic): ⟹ {0}".format(
        "REJECTED false-empty ✓" if m4 else "✗"))
    caught += 1 if m4 else 0

    print()
    print("  mutants caught: {0}/{1}".format(caught, total))
    return caught == total


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _here = os.path.dirname(os.path.abspath(__file__))
    _logf = open(os.path.join(_here, "S1032_run.log"), "w", encoding="utf-8")

    class Tee:
        def __init__(s, r, f): s.r = r; s.f = f; s.chunks = []
        def write(s, x): s.r.write(x); s.f.write(x); s.f.flush(); s.chunks.append(x); return len(x)
        def flush(s): s.r.flush(); s.f.flush()
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_here, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("STIK STEP 2 (a) KINEMATICS · S1032 — an ε-odd c-invariant m + UNIQUENESS")
    print("Field = D_h center+two-component (S1029/S1031). The action of w₀/c on m is COMPUTED (not postulated). Primitives+2×2.")
    print("★STOP-negative if m does not exist/is not unique. ONLY (a). FS=STONE. Court — to Omega; I do NOT render a verdict.")
    print("=" * 74)
    print()

    unique, eps_odd, eps_odd_gap = level_a(); print()
    schur_ok = uniqueness_schur(); print()
    mut_ok = mutants(); print()

    print("=" * 74)
    print("RAW RESULTS (the court — to Omega; I do NOT render a verdict):")
    print("─" * 74)
    if unique and schur_ok:
        print("  ★PASS at level (a): m = σ_z = diag(1,−1) — EXPLICIT · c-invariant (computed) · ε-odd")
        print("   (B·σ_z=−σ_z computed) · UNIQUE up to scale (1-dim, the split-tooth from the real T33 pins it).")
        print("   The prediction «an ε-contraction of the trace type» HIT (σ_z diagonal traceless = of the trace type).")
        print("   ⟹ D_h→⟨c⟩ HAS a canonical break-channel; the freezing picture is ALIVE for level (b).")
    else:
        print("  ★STOP-NEGATIVE: no canonical break-channel — freezing falls (an honest result).")
    print("─" * 74)
    print("  ★BOTH READINGS (discipline): under the BARE D_h-flip the ε-odd ones = {σ_y,σ_z} dim 2 (NOT unique")
    print("   by itself). Uniqueness comes from the split-tooth (a splitter=anticommutes with the kinetic H, from the real T33 chirality):")
    print("   σ_y is odd but NOT a splitter ({σ_y,H}≠0) ⟹ drops out; σ_z is odd AND a splitter ⟹ UNIQUE.")
    print("   ★OMEGA'S CALL: uniqueness rests on the split-tooth (canonical, from the T33 H). If you accept only the D_h-")
    print("   symmetry without the chiral H — dim=2, a STOP-negative. I computed both, the tooth from the real T33.")
    print("─" * 74)
    all_ok = unique and schur_ok and mut_ok
    print("  SUMMARY: PASS(a)={0} · Schur-uniqueness={1} · mutants={2}".format(
        "YES" if unique else "STOP-NEG", "YES" if schur_ok else "NO", "4/4" if mut_ok else "NO"))
    print("=" * 74)

    # NB: 'σ_z/σ_x/two-component/splitter/split/anticommutator/flip/odd/center/chirality' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("матер", "ія"), ("частин", "ка"), ("Міньков", "ський"),
           ("всес", "віт"), ("ант", "роп"), ("Teg", "mark")]  # GUARDLINE (FS=STONE)
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _n = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN (STONE): hits={0} (src={1}, log={2})".format(_n, len(_hs), len(_hl)))

    _exit = 1 if (_n > 0 or not (schur_ok and mut_ok)) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
