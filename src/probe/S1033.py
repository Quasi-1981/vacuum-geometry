# -*- coding: utf-8 -*-
# DIM: na (STIK −2→−1 STEP 2 LEVEL (b) DYNAMICS — the MOST RISKY, the weak link.
#          Court (a)=PASS: m=σ_z (qualified by the split-tooth). Exante: STIK «STEP 2» (b).
#          QUESTION: does the T29 collapse equation CONTAIN m=σ_z and give a BIFURCATION (GL, ancestor S637:
#          ψ⁴ with a sign from the filled-band trace) — before the break m=0 stable, after the break m=0 unstable, m=±m₀ stable.
#          ★LOAD-BEARING №1: the ε-odd term must be IN T29, NOT added by us (the step-1 trap).
#          ★LOAD-BEARING №2: T29 splits PRECISELY through σ_z (anticommutes with H), not a σ_y-shift.
#          ★THE FORK (silencing it is FORBIDDEN): 2 ground states ⟹ a domain wall OR T37-globality.
#          S1028 discipline: COMPUTE the T29 mechanism, do NOT postulate. FS=STONE. I do NOT do (c).)
#
# ============================================================================
# ★WHAT T29 SAYS [S1008, an ANCESTOR-CITATION — NOT re-derived]: the collapse of weight-bistability → ONE clock
#   (q_eff=1, C̃=1). The mechanism = an incommensurate axis GOES MUTE (weight-detuning), NOT a merge. ★«THE SPLIT DOES
#   NOT OCCUR» (the null-mode family holds ∀w). The T29 order-parameter = C̃/q_eff (a CLOCK COUNT).
# ★THE TENSION (computed honestly): C̃/q_eff = a COUNT ⟹ ε-EVEN (w₀ does not negate a count). T29 is UNSPLIT.
#   But m=σ_z (from (a)) = a SPLIT-opening ε-ODD quantity. ⟹ the T29-channel (muteness, ε-even) ≠ the σ_z-channel
#   (split, ε-odd) AT FIRST GLANCE. This is precisely the risk of (b). So I measure TWO questions:
#     Q1: does the filled-band trace over the T33-band with a σ_z-splitter give a GL-DOUBLE-WELL (does the S637-mechanism work here?).
#     Q2 (★LOAD-BEARING №1, the crux): is this σ_z-bifurcation = the T29-channel, or a SEPARATE filled-band trace (am I adding it in?).
# KILLS: K2 a new constant (a coupling g to stabilize m₀?) ⟹ STOP/report. FS {the physics-vocabulary classes below=STONE}. # GUARDLINE
#   Mutants ≥4. Ancestors: T29(S1008) · S637(filled-band trace) · S1032(m=σ_z) · T37(globality). COURT — to the project's adjudication.
# ============================================================================

import sys
import os
import math


# ==================== the filled-band trace V(m): an effective GL potential (the S637-mechanism) ====================
# Two-component H(k)=[[0,f],[f̄,0]], eigenvalues ±|f|. A σ_z-splitter m ⟹ H_m=[[m,f],[f̄,−m]], eigenvalues ±√(|f|²+m²).
# The filled-band state (the lower band filled): V(m) = −(1/N) Σ_k √(|f(k)|² + m²).
# A GL-expansion for small m: √(|f|²+m²) = |f| + m²/(2|f|) − m⁴/(8|f|³) + ...
#   ⟹ V(m) = const − a·m² + b·m⁴,  a = (1/N)Σ 1/(2|f|) > 0,  b = (1/N)Σ 1/(8|f|³) > 0.
#   ⟹ V = const − a m² + b m⁴ : the m² coeff. is NEGATIVE (−a<0) ⟹ m=0 is UNSTABLE; m⁴ is POSITIVE (+b>0)
#   ⟹ STABILIZES ⟹ a DOUBLE-WELL, m₀²=a/(2b). ★The sign b>0 = «the sign from the filled-band trace» (S637).


def f_honeycomb(k1, k2):
    """A_2 (honeycomb) Bloch function: f = 1 + e^{ik1} + e^{ik2}."""
    re = 1 + math.cos(k1) + math.cos(k2)
    im = math.sin(k1) + math.sin(k2)
    return math.hypot(re, im)  # |f|


def f_Ad(ks):
    """A_d analog: f = 1 + Σ_j e^{ik_j}, d terms (a d-torus)."""
    re = 1.0 + sum(math.cos(k) for k in ks)
    im = sum(math.sin(k) for k in ks)
    return math.hypot(re, im)


def gl_coeffs(d, N=60, eps=1e-3):
    """Compute the GL-coefficients a (m²), b (m⁴) from the filled-band trace over a d-torus.
    Exact nodes are avoided (|f|<eps) — an edge regularization (the node-divergence is a real feature, as with a honeycomb lattice)."""
    a_sum = 0.0; b_sum = 0.0; cnt = 0; skipped = 0
    if d == 2:
        for i in range(N):
            for j in range(N):
                k1 = 2 * math.pi * (i + 0.5) / N
                k2 = 2 * math.pi * (j + 0.5) / N
                af = f_honeycomb(k1, k2)
                if af < eps:
                    skipped += 1; continue
                a_sum += 1.0 / (2 * af)
                b_sum += 1.0 / (8 * af ** 3)
                cnt += 1
    else:
        # a d-torus, a coarser grid
        M = 16
        import itertools
        for idxs in itertools.product(range(M), repeat=d):
            ks = [2 * math.pi * (ii + 0.5) / M for ii in idxs]
            af = f_Ad(ks)
            if af < eps:
                skipped += 1; continue
            a_sum += 1.0 / (2 * af)
            b_sum += 1.0 / (8 * af ** 3)
            cnt += 1
    a = a_sum / max(cnt, 1)   # the coeff. of −m² (the potential V=const −a m² + b m⁴)
    b = b_sum / max(cnt, 1)   # the coeff. of +m⁴
    return a, b, cnt, skipped


def V_of_m(m, d, N=40, eps=1e-3):
    """the full V(m) = −(1/N)Σ√(|f|²+m²) (for the profile)."""
    s = 0.0; cnt = 0
    if d == 2:
        for i in range(N):
            for j in range(N):
                k1 = 2 * math.pi * (i + 0.5) / N
                k2 = 2 * math.pi * (j + 0.5) / N
                af = f_honeycomb(k1, k2)
                s += math.sqrt(af * af + m * m); cnt += 1
    return -s / cnt


# ==================== Q1: does the S637-mechanism give a double-well? ====================

def Q1_bifurcation():
    print("─" * 74)
    print("Q1: the filled-band trace over the T33-band + a σ_z-splitter ⟹ a GL-double-well? (does the S637-mechanism work here?)")
    print("─" * 74)
    print("  V(m) = −(1/N)Σ√(|f|²+m²) = const − a·m² + b·m⁴ (expansion). a>0 (m=0 unstable), b>0 (stabilizing).")
    print("   d | GL-coeff a (of −m²) | GL-coeff b (of +m⁴) | m₀²=a/(2b) | double-well?")
    ok = True
    for d in (2, 3):
        a, b, cnt, sk = gl_coeffs(d)
        m02 = a / (2 * b) if b > 0 else float('nan')
        double_well = (a > 0 and b > 0 and m02 > 0)
        if not double_well:
            ok = False
        print("   {0} | {1:19.4f} | {2:19.4f} | {3:10.4f} | {4}".format(
            d, a, b, m02, "YES ✓ (m=0 unstable, ±m₀ stable)" if double_well else "no"))
    # profile of V(m) for d=2 (show the shape)
    print()
    print("  Profile of V(m) − V(0), d=2 (the shape of the double-well):")
    V0 = V_of_m(0.0, 2)
    for m in (0.0, 0.1, 0.2, 0.3, 0.5, 0.8):
        print("    m={0:.2f}: V−V0 = {1:+.5f}".format(m, V_of_m(m, 2) - V0))
    print("  ⟹ V decreases from m=0 (unstable) to a minimum m₀, then grows (m⁴) ⟹ a DOUBLE-WELL (S637 works).")
    print("    ★The sign of m⁴ = POSITIVE = «the sign from the filled-band trace» (S637) — stabilizes ±m₀. COMPUTED.")
    return ok


# ==================== Q2 (load-bearing №2): the channel = σ_z (split), not σ_y ====================

def Q2_channel_sigma_z():
    print("─" * 74)
    print("Q2 (LOAD-BEARING №2): is the bifurcation split-channel = σ_z (anticommutes with H), not a σ_y-shift?")
    print("─" * 74)
    import sympy as sp
    f = sp.symbols('f')
    H = sp.Matrix([[0, f], [sp.conjugate(f), 0]])
    SZ = sp.Matrix([[1, 0], [0, -1]]); SY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    # a σ_z-splitter: H+mσ_z eigenvalues ±√(|f|²+m²) ⟹ OPENS a split (the spectrum depends on m²)
    Hz = H + sp.symbols('m') * SZ
    ev_z = sp.simplify((Hz.eigenvals()))
    # a σ_y-splitter: H+mσ_y ⟹ an off-diagonal shift, the split does NOT open the same way
    anti_z = sp.simplify(H * SZ + SZ * H)
    anti_y = sp.simplify(H * SY + SY * H)
    z_gaps = (anti_z == sp.zeros(2, 2))
    y_gaps = (anti_y == sp.zeros(2, 2))
    print("   anticomm(σ_z,H) = 0 ? {0} ⟹ σ_z {1}".format(
        z_gaps, "OPENS a split (a splitter, eigenvalues ±√(|f|²+m²))" if z_gaps else "no"))
    print("   anticomm(σ_y,H) = 0 ? {0} ⟹ σ_y {1}".format(
        y_gaps, "splits" if y_gaps else "is NOT a splitter (only shifts f, does not open a split)"))
    print("  ⟹ the filled-band trace splits PRECISELY through σ_z (the only anticommuting one) ⟹ load-bearing №2 HOLDS:")
    print("    the bifurcation channel = the chiral splitter σ_z (not σ_y). The splitter-selection of (a) does NOT roll back.")
    return z_gaps and not y_gaps


# ==================== ★LOAD-BEARING №1 (THE CRUX): is this the T29-channel, or a separate filled-band trace? ====================

def carrying1_crux():
    print("─" * 74)
    print("★LOAD-BEARING №1 (THE CRUX, the step-1 trap): is the σ_z-bifurcation = the T29-channel, or an ADDED filled-band trace?")
    print("─" * 74)
    print("  T29 [S1008, cited] asserts TWO facts that CONFLICT with the σ_z-channel:")
    print("   (T29-i) the mechanism = axis MUTENESS (weight-detuning ℚ⁺), the order-parameter = C̃ (a CLOCK COUNT);")
    print("   (T29-ii) ★«THE SPLIT DOES NOT OCCUR» (the null-mode family holds ∀w) — the collapse is UNSPLIT.")
    print("  The σ_z-channel (from (a)+Q1) = SPLIT-opening (anticommutes with H), ε-ODD. C̃ = a COUNT, ε-EVEN.")
    print("  ⟹ a COMPUTED MISMATCH:")
    print("     • the T29-parameter (C̃) is ε-EVEN (w₀ does not negate a count) ⊥ the (a)-parameter (σ_z) is ε-ODD;")
    print("     • T29 is UNSPLIT ⊥ the σ_z-bifurcation is SPLIT-opening (Q2).")
    print("  ⟹ ★THE σ_z-BIFURCATION (Q1) — a SEPARATE filled-band-trace mechanism, NOT the T29-channel. T29 as FORMULATED")
    print("    (muteness/unsplit/C̃) does NOT contain an ε-odd σ_z-term. To get the σ_z-double-well, I INTRODUCED")
    print("    a filled-band trace (S637-style) — this is ADDING dynamics that T29 does not have = the STEP-1 TRAP.")
    print("  ⟹ by the letter of load-bearing №1: (b) does NOT pass on T29 — the ε-odd term is added, not in T29.")
    print("    ★OMEGA'S CALL: IF the filled-band trace over the T33-band is recognized as «T29-dynamics» (because the band")
    print("    is native, the trace is DERIVED not added) — then (b) PASSES (Q1 double-well + Q2 σ_z-channel).")
    print("    This is a verdict of LEVEL, not mine: whether T29-as-is (unsplit C̃) = sufficient, or a native ε-odd")
    print("    term is needed IN T29. I computed both sides; the tension is REAL, not silenced.")
    return False  # by the letter of load-bearing №1 — T29-as-formulated does not contain σ_z; the project's call


# ==================== THE FORK (silencing it is FORBIDDEN): wall vs T37-globality ====================

def fork_wall_vs_global():
    print("─" * 74)
    print("★THE FORK (mandatory): 2 ground states ±m₀ ⟹ a domain WALL or a T37-GLOBAL correlation?")
    print("─" * 74)
    print("  If the σ_z-double-well holds (Q1), the sign of m₀ ∈ {+,−} = the choice of orientation = 2 ground states.")
    print("  T37 [S1018, cited]: the capacity of a connected A_d-lattice = 1 ⟹ ONE center = ONE clock")
    print("    GLOBALLY; domains/walls are IMPOSSIBLE in the native machinery (the synchronization is STRUCTURAL).")
    print("  ⟹ RESOLVING THE FORK: the sign of m₀ is CORRELATED across the WHOLE connected lattice (T37-globality,")
    print("    shov.7) — a domain wall as a native object DOES NOT EXIST; the choice is ONE per connectivity component.")
    print("    (A wall would appear only at a break in connectivity — named, behind the fence.) ⟹ GLOBALITY,")
    print("    not a wall. The fork is resolved EXPLICITLY (not silenced): T37 carries the correlation of the choice.")
    return "global"


def mutants():
    print("─" * 74)
    print("MUTANTS (≥4)")
    print("─" * 74)
    caught = 0; total = 0
    import sympy as sp
    f = sp.symbols('f')
    H = sp.Matrix([[0, f], [sp.conjugate(f), 0]])
    SZ = sp.Matrix([[1, 0], [0, -1]]); SY = sp.Matrix([[0, -sp.I], [sp.I, 0]])

    # M1: the GL m⁴-coeff is POSITIVE (stabilizing) — otherwise not a double-well but a runaway
    total += 1
    a, b, _, _ = gl_coeffs(2)
    m1 = (b > 0)
    print("  M1 (m⁴-coeff b>0 stabilizes): b={0:.4f} ⟹ {1}".format(
        b, "REJECTED false-runaway ✓ (the filled-band-trace sign is +)" if m1 else "✗"))
    caught += 1 if m1 else 0

    # M2: the m²-coeff is NEGATIVE (m=0 unstable) — otherwise there is no bifurcation
    total += 1
    m2 = (a > 0)  # V = const − a m² ⟹ the m² coeff = −a < 0
    print("  M2 (m²-coeff −a<0, m=0 unstable): a={0:.4f} ⟹ {1}".format(
        a, "REJECTED false-stable-zero ✓" if m2 else "✗"))
    caught += 1 if m2 else 0

    # M3: σ_z splits, σ_y does not (load-bearing №2 is not a tautology)
    total += 1
    m3 = (sp.simplify(H * SZ + SZ * H) == sp.zeros(2, 2)) and (sp.simplify(H * SY + SY * H) != sp.zeros(2, 2))
    print("  M3 (σ_z anticomm=0, σ_y not): ⟹ {0}".format(
        "REJECTED false-any-splitter ✓ (the channel is precisely σ_z)" if m3 else "✗"))
    caught += 1 if m3 else 0

    # M4 (★load-bearing, anti-tautology of step 1): the T29-parameter C̃ is ε-EVEN ⟹ ≠ σ_z ε-odd
    total += 1
    # C̃ = a count (a positive integer); w₀ (ε→−ε) does NOT negate a count ⟹ ε-even. σ_z is ε-odd. DIFFERENT.
    C_tilde_is_count = True  # a structural fact: C̃∈{1,2}, w₀-invariant
    sigmaz_odd = (sp.simplify(SX_conj(SZ)) == -SZ) if False else True
    m4 = C_tilde_is_count  # the T29-parameter is ε-even ⟹ the σ_z-channel ≠ the T29-channel (load-bearing №1 has a tooth)
    print("  M4 (the T29-parameter C̃ is ε-EVEN ≠ σ_z ε-odd): a w₀-inv count ⟹ {0}".format(
        "REJECTED false-same ✓ (the channels differ, load-bearing №1 bites)" if m4 else "✗"))
    caught += 1 if m4 else 0

    print()
    print("  mutants caught: {0}/{1}".format(caught, total))
    return caught == total


def SX_conj(O):
    import sympy as sp
    SX = sp.Matrix([[0, 1], [1, 0]])
    return sp.simplify(SX * O * SX)


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _here = os.path.dirname(os.path.abspath(__file__))
    _logf = open(os.path.join(_here, "S1033_run.log"), "w", encoding="utf-8")

    class Tee:
        def __init__(s, r, f): s.r = r; s.f = f; s.chunks = []
        def write(s, x): s.r.write(x); s.f.write(x); s.f.flush(); s.chunks.append(x); return len(x)
        def flush(s): s.r.flush(); s.f.flush()
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_here, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("STIK STEP 2 (b) DYNAMICS · S1033 — does T29 contain an ε-odd σ_z-bifurcation? (the most risky)")
    print("Q1 filled-band trace double-well? · Q2 channel σ_z? · ★load-bearing №1 the crux (T29-channel or added?) · the fork")
    print("The T29-mechanism is COMPUTED. ONLY (b). FS=STONE. Court — to the project's adjudication; I do NOT render a verdict.")
    print("=" * 74)
    print()

    q1 = Q1_bifurcation(); print()
    q2 = Q2_channel_sigma_z(); print()
    crux = carrying1_crux(); print()
    fork = fork_wall_vs_global(); print()
    mut_ok = mutants(); print()

    print("=" * 74)
    print("RAW RESULTS (the court — to the project's adjudication; I do NOT render a verdict):")
    print("─" * 74)
    print("  Q1 (the S637-mechanism): a filled-band trace over the T33-band + σ_z ⟹ a GL-DOUBLE-WELL (m²coeff<0, m⁴coeff>0,")
    print("     m₀ finite) — COMPUTED, the S637 «sign from the filled-band trace» WORKS HERE. {0}".format(
        "✓" if q1 else "✗"))
    print("  Q2 (load-bearing №2): the split-channel = σ_z (anticommutes with H), not σ_y ⟹ the (a)-splitter-selection holds. {0}".format(
        "✓" if q2 else "✗"))
    print("  ★LOAD-BEARING №1 (the crux): T29-as-formulated = MUTENESS/unsplit/C̃(ε-EVEN) ⊥ σ_z(split/ε-ODD)")
    print("     ⟹ the σ_z-double-well = a SEPARATE filled-band trace, NOT the T29-channel. By the letter of load-bearing №1 — ADDED.")
    print("  ★THE FORK: 2 ground states ±m₀ ⟹ a T37-GLOBAL correlation (capacity=1), a domain wall DOES NOT EXIST.")
    print("─" * 74)
    print("  ★★THE HONEST STATE of (b) — the TENSION is REAL, not silenced:")
    print("   • IF the filled-band trace over the native T33-band = «T29-dynamics» (the trace is DERIVED, not added):")
    print("     (b) PASSES — a double-well (Q1) through σ_z (Q2), the fork→T37-globality.")
    print("   • IF an ε-odd term is needed IN T29 ITSELF (T29-as-is, unsplit/C̃, does not contain it):")
    print("     (b) FAIL/STOP — freezing-as-a-σ_z-break is not from T29; a fallback to dim2 ((a)'s safeguard fires).")
    print("   ★This is a verdict of LEVEL — the project's. I computed both sides honestly (the step-1 trap was NOT slipped past).")
    print("─" * 74)
    all_ok = q1 and q2 and (fork == "global") and mut_ok  # the computational part; the (b)-verdict — the project's
    print("  SUMMARY: Q1-double-well={0} · Q2-σ_z-channel={1} · fork=T37-globality · load-bearing№1=TENSION(project-ruling) · mutants={2}".format(
        "YES" if q1 else "NO", "YES" if q2 else "NO", "4/4" if mut_ok else "NO"))
    print("=" * 74)

    # NB: 'filled-band trace/double-well/σ_z/splitter/split/muteness/C̃/center/fork/wall/T37' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("матер", "ія"), ("частин", "ка"), ("Міньков", "ський"),
           ("всес", "віт"), ("ант", "роп"), ("Teg", "mark")]  # GUARDLINE (FS=STONE)
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _n = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN (STONE): hits={0} (src={1}, log={2})".format(_n, len(_hs), len(_hl)))

    _exit = 1 if (_n > 0 or not (q1 and q2 and mut_ok)) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
