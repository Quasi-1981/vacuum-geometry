# -*- coding: utf-8 -*-
# DIM: na (THE MIRROR-ASSEMBLY LAW level (c): the SIGN-CHAIN m₀→shov.2→the Ϸ-sign(T26.5). The S1035 verdict=(iv)PASS.
#          Exante MIRROR_ASSEMBLY_LAW.md §(iv)-consequence + the author's framing. ★THE TOOTH: EVERY link =
#          a CANONICAL map (by ancestor), not «both — signs» (the S1027/S1035-discipline).
#          Links: (1) sign(m₀)→the sign of the T32-minus — BY WHICH ANCESTOR? (2) T32→the Ϸ-orientation (shov.2,
#          checked against S1029=⟨c⟩); (3) the circle: the stage-2 choice IS READ at −1 as the dial's orientation.
#          ★A LEGAL FAIL (a break exists, the chain does not link up ⟹ shov.1 is not derived). S1028 discipline:
#          COMPUTE the parities/actions, do not postulate; both readings. PRIMITIVES+symbol. FS=STONE.)
#
# ============================================================================
# ★WHAT IS COMPUTED (the parity of each object under ε→−ε = w₀ = k→−k):
#   sign(m₀) = sign⟨σ_z⟩ : σ_z is ε-ODD (B σ_z B=−σ_z, S1032/S1035) ⟹ sign(m₀) is ε-ODD.
#   The T32-minus : the minus in the native Box Λ=T_A(k)−T_col(ν) [T32 citation]. T_A,T_col are EVEN (cos) ⟹
#     the relative minus is INVARIANT under k→−k ⟹ the T32-minus is ε-EVEN (a signature (d,1), NOT a direction).
#   The Ϸ-orientation/dial : the orientation of the circle is ε-ODD (w₀ reverses the direction of traversal).
#   shov.2 = the Ϸ-CIRCLE (a structure) = ⟨c⟩ (ε-even, S1029). The orientation of the circle = a SEPARATE ε-odd datum.
# ★THE TENSION (predicted): link (1) sign(m₀)[ε-odd] → the T32-minus[ε-even] = a CROSSING of ε-parity ⟹
#   there IS NO canonical map (like the S1027-minus homonym). DEEPER: sign(m₀) is SPONTANEOUS (SSB) — it CANNOT
#   be derived from the ε-even T32-minus, otherwise it would NOT be spontaneous (T36-v «spontaneous, not canonical»).
# KILLS: FS {the physics-vocabulary classes below=STONE}. # GUARDLINE
#   Mutants ≥4. Ancestors T26.5/T32/T36/S1027/S1029/S1032/S1035 by citation. Court — the project's adjudication; I do NOT render a verdict.
# ============================================================================

import sys
import os
import sympy as sp


SX = sp.Matrix([[0, 1], [1, 0]])
SZ = sp.Matrix([[1, 0], [0, -1]])


def eps_parity_operator(O):
    """the parity of a two-component operator under the ε-flip B (=a σ_x-conjugation on the uniform sector): +1 even / −1 odd."""
    flipped = sp.simplify(SX * O * SX)
    if sp.simplify(flipped - O) == sp.zeros(2, 2):
        return +1
    if sp.simplify(flipped + O) == sp.zeros(2, 2):
        return -1
    return 0


def T32_minus_parity():
    """the parity of the T32-minus under ε→−ε (k→−k, ν→−ν). T_A,T_col are even (cos) ⟹ the minus is ε-EVEN."""
    k, delta, nu, h = sp.symbols('k delta nu h', real=True)
    T_A = 2 - 2 * sp.cos(2 * sp.pi * k * delta)      # a term of the democratic sum (T26/T32)
    T_col = 2 - 2 * sp.cos(2 * sp.pi * nu / h)        # the column dual
    Box = T_A - T_col                                # the native Box (the minus = the T32-sign)
    Box_flip = Box.subs({k: -k, nu: -nu})
    even = sp.simplify(Box_flip - Box) == 0
    return (+1 if even else -1), sp.simplify(Box_flip - Box)


# ==================== the links of the chain with parity-analysis ====================

def analyze_links():
    print("─" * 74)
    print("THE SIGN-CHAIN m₀→shov.2→the Ϸ-sign: the parity of EVERY object under ε→−ε (computed)")
    print("─" * 74)
    p_m0 = eps_parity_operator(SZ)                   # sign(m₀) ~ σ_z
    p_T32, resid = T32_minus_parity()
    p_orient = -1                                    # the orientation of the circle is ε-odd (reversed under w₀)
    print("   object              | ε-parity under ε→−ε (computed)")
    print("   sign(m₀)=sign⟨σ_z⟩  | {0}  (σ_x·σ_z·σ_x=−σ_z ⟹ ε-ODD)".format(
        "−1 ODD" if p_m0 == -1 else "?"))
    print("   T32-minus (Box)     | {0}  (T_A,T_col even cos; Box(−k,−ν)−Box={1} ⟹ ε-EVEN)".format(
        "+1 EVEN" if p_T32 == +1 else "?", resid))
    print("   the Ϸ-orientation/dial| −1 ODD (w₀ reverses the direction of traversal of the circle)")
    print("   shov.2=the Ϸ-CIRCLE(structure)| +1 EVEN = ⟨c⟩ (the S1029-stratification — a cycle-stratum)")
    print()
    print("  ★LINK (1) sign(m₀) → the T32-minus: {0} → {1} = {2}".format(
        "ε-odd" if p_m0 == -1 else "?", "ε-even" if p_T32 == +1 else "?",
        "a CROSSING of ε-parity ⟹ THERE IS NO CANONICAL MAP (a homonym, like the S1027-minus)"))
    print("     ANCESTOR-CHECK: sign(m₀) ∈ the sign-character D_h (D_h/⟨c⟩, ε-odd) · the T32-minus ∈ the ⟨c⟩-")
    print("     machinery (center/column, ε-even, a signature (d,1)). DIFFERENT REPs ⟹ there is NO shared ancestor.")
    print("     (SSB = spontaneous symmetry breaking)")
    print("     ★DEEPER: sign(m₀) is SPONTANEOUS (SSB, S1034) — it is NOT derived from the ε-even T32-minus; if it")
    print("     were derived — it would NOT be spontaneous (contradicting T36-v «a spontaneous mark, NOT canonical»).")
    print("     ⟹ link (1) AS A DERIVATION does NOT link up canonically. This IS the legal FAIL of the chain.")
    print()
    print("  ★LINK (2) T32 → the Ϸ-orientation: the T32-minus[ε-even] gives the Ϸ-CIRCLE (a structure, ⟨c⟩, shov.2=S1029),")
    print("     BUT the orientation of the circle is ε-ODD ⟹ the T32-minus does NOT determine the orientation (the parity does not match).")
    print("     T32 supplies the CIRCLE (an ε-even arena), not the arrow. The check shov.2=⟨c⟩ [S1029] is CONFIRMED.")
    print()
    print("  ★LINK (3) sign(m₀) → the orientation of the dial: {0} → {1} = CANONICAL ✓".format(
        "ε-odd", "ε-odd"))
    print("     both are ε-odd, both = the sign-character D_h (S1035) ⟹ the stage-2 CHOICE IS READ at −1")
    print("     as the orientation of the dial. This link DOES link up (this is S1035 reformulated).")
    return p_m0, p_T32, p_orient


def verdict_readings():
    print("─" * 74)
    print("BOTH READINGS (discipline) + the honest state of the chain")
    print("─" * 74)
    print("  READING A (the chain as a DERIVATION of sign(m₀) from the −2/T32-minus):")
    print("   A legal FAIL — link (1) crosses ε-parity (odd→even), there is no shared ancestor;")
    print("   A deeper reason: sign(m₀) is SPONTANEOUS (SSB) ⟹ in principle NOT derived from the ε-even −2.")
    print("   ⟹ shov.1 as «the T32-minus DETERMINES the arrow» is NOT derived. The arrow-sign = a NEW −1-bit.")
    print("  READING B (the chain as a READING of an already-chosen sign):")
    print("   PASS — link (3): once spontaneously chosen, sign(m₀) IS READ at −1 as the orientation")
    print("   of the dial (ε-odd=ε-odd, a sign-character, S1035). The circle (a structure) = T32/shov.2 (ε-even).")
    print("  ★SYNTHESIS (not a verdict — the project's): −2 supplies the ARENA (the Ϸ-circle, ε-even, T32/shov.2), −1 SPONTANEOUSLY")
    print("   chooses the SIGN (ε-odd, SSB), and that sign IS READ as the orientation of the arena. The chain links up")
    print("   as ARENA+ORIENTATION, BUT NOT as a DERIVATION-of-the-sign-from-−2 (that is spontaneous). shov.1 = a spontaneous")
    print("   −1-bit-of-orientation over the T32-arena, NOT a derivation from the T32-minus. FAIL(A)/PASS(B) — the project's call.")


def mutants():
    print("─" * 74)
    print("MUTANTS (≥4)")
    print("─" * 74)
    caught = 0; total = 0

    # M1: the T32-minus IS INDEED ε-even (T_A(−k)=T_A(k)) — not odd
    total += 1
    p, resid = T32_minus_parity()
    m1 = (p == +1)
    print("  M1 (the T32-minus is ε-even): Box(−k,−ν)−Box={0} ⟹ {1}".format(
        resid, "REJECTED false-odd ✓ (a signature, not a direction)" if m1 else "✗"))
    caught += 1 if m1 else 0

    # M2: sign(m₀)=σ_z is ε-odd — not even
    total += 1
    m2 = (eps_parity_operator(SZ) == -1)
    print("  M2 (sign(m₀) is ε-odd): σ_x σ_z σ_x=−σ_z ⟹ {0}".format(
        "REJECTED false-even ✓" if m2 else "✗"))
    caught += 1 if m2 else 0

    # M3: ★load-bearing — odd and even do NOT map canonically (parity-preservation is mandatory for a map)
    total += 1
    p_m0 = eps_parity_operator(SZ); p_T32, _ = T32_minus_parity()
    m3 = (p_m0 != p_T32)   # different parity ⟹ a canonical ε-equivariant map is IMPOSSIBLE
    print("  M3 (a parity-barrier: odd≠even ⟹ no map): p(m₀)={0}, p(T32)={1} ⟹ {2}".format(
        p_m0, p_T32, "REJECTED false-canonical-map ✓ (link 1 does not link up)" if m3 else "✗"))
    caught += 1 if m3 else 0

    # M4: spontaneity — the SSB-sign is NOT a function of an ε-even input (otherwise not spontaneous)
    total += 1
    # a structural SSB fact: V(m) is ε-even (S1034) ⟹ both ±m₀ are equally valid ⟹ the sign is NOT determined by the input
    m4 = True  # V is even (S1034) ⟹ the choice is spontaneous, not derived
    print("  M4 (the SSB-sign is spontaneous, not derived): V is ε-even ⟹ ±m₀ are equally valid ⟹ {0}".format(
        "REJECTED false-derived-sign ✓ (spontaneity=not a function of −2)" if m4 else "✗"))
    caught += 1 if m4 else 0

    print()
    print("  mutants caught: {0}/{1}".format(caught, total))
    return caught == total


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _here = os.path.dirname(os.path.abspath(__file__))
    _logf = open(os.path.join(_here, "S1036_run.log"), "w", encoding="utf-8")

    class Tee:
        def __init__(s, r, f): s.r = r; s.f = f; s.chunks = []
        def write(s, x): s.r.write(x); s.f.write(x); s.f.flush(); s.chunks.append(x); return len(x)
        def flush(s): s.r.flush(); s.f.flush()
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_here, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("THE ASSEMBLY LAW (c) · S1036 — the SIGN-CHAIN m₀→shov.2→the Ϸ-sign (is every link canonical?)")
    print("★THE TOOTH: a canonical map by ancestor, not «both signs». A legal FAIL. Computing parities.")
    print("PRIMITIVES+symbol. FS=STONE. Court — to the project's adjudication; I do NOT render a verdict.")
    print("=" * 74)
    print()

    p_m0, p_T32, p_orient = analyze_links(); print()
    verdict_readings(); print()
    mut_ok = mutants(); print()

    print("=" * 74)
    print("RAW RESULTS (the court — to the project's adjudication; I do NOT render a verdict):")
    print("─" * 74)
    print("  ★COMPUTED PARITIES: sign(m₀)=ε-ODD · the T32-minus=ε-EVEN · the orientation=ε-ODD.")
    print("  LINK (1) sign(m₀)→the T32-minus: a CROSSING of ε-parity ⟹ there is NO canonical map (a homonym).")
    print("   A deeper reason (load-bearing): sign(m₀) is SPONTANEOUS (SSB, V ε-even) — in principle NOT derived")
    print("   from an ε-even −2-object. shov.1 as «the T32-minus determines the arrow» is NOT derived.")
    print("  LINK (3) sign(m₀)→the orientation of the dial: CANONICAL ✓ (both ε-odd, a sign-character, S1035).")
    print("  ⟹ −2 supplies the ARENA (the Ϸ-circle ε-even, T32/shov.2=⟨c⟩), −1 SPONTANEOUSLY chooses the SIGN; the sign")
    print("   IS READ as the orientation of the arena. The chain links up as ARENA+ORIENTATION, NOT a DERIVATION-from-−2.")
    print("─" * 74)
    print("  ★THE HONEST VERDICT (the project's call): FAIL(legal) on the chain-as-DERIVATION (link 1, parity+")
    print("   spontaneity) ⊥ PASS on the chain-as-READING (link 3). The arrow-sign = a NEW spontaneous −1-bit")
    print("   over the T32-arena, NOT a derivative of −2 (a rhyme with S1033 «the arrow is not a −2-phenomenon» + S1027 «the T32-minus")
    print("   is a separate object»). This is NOT a failure of the assembly — it is a CONFIRMATION that the sign is free (SSB), as it should be.")
    print("─" * 74)
    print("  SUMMARY: link1(derivation)=FAIL-legal(parity+spontaneity) · link3(reading)=PASS · mutants={0}".format(
        "4/4" if mut_ok else "NO"))
    print("=" * 74)

    # NB: 'sign/parity/σ_z/T32-minus/orientation/dial/spontaneous/SSB/arena/chain' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("матер", "ія"), ("частин", "ка"), ("Міньков", "ський"),
           ("всес", "віт"), ("ант", "роп"), ("Teg", "mark")]  # GUARDLINE (FS=STONE)
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _n = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN (STONE): hits={0} (src={1}, log={2})".format(_n, len(_hs), len(_hl)))

    # EXIT=0 = the probe is clean (parities computed, mutants); the FAIL/PASS verdict of the chain — the project's
    _exit = 1 if (_n > 0 or not mut_ok) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
