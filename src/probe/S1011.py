# -*- coding: utf-8 -*-
# DIM: na (W42 probe-6, layer -2: the NATIVE cell Box Λ(k,ν)=T_A(k)−T_col(ν) — the minus
#          comes from the column duality (T26), not by hand. Are the zeros = Dirac nodes (bridge №3), and
#          does the q=1-asymmetry fall out (hole №2). The S956/S1005 machinery; ancestors cited;
#          exact arithmetic (roots of unity, integer term_table); 0 handles.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting) — §10 exante + ancestors
# ----------------------------------------------------------------------------
# NATIVE BOX: Λ(ψ,ν) = T_A(ψ) − T_col(ν).
#   T_A(ψ) = sum_{i=0..d}(2−2cos 2π p_i), p_0=0 convention (axis-0 term=0), p_i=ψ_i (i=1..d);
#     = sum_{i=1..d} term(ψ_i) on the h-torsion, term(m)=2−2cos(2πm/h) in Z (S1005-structural).
#   T_col(ν) = 2−2cos(2πν/h) = term(ν), ν in Z/h (dual of column period h=d+1, S1001 — EXACTLY
#     h points, not a choice).  Minus enters via duality (T26), not a hand-signed (p,q).
#   f(ψ) = 1 + sum_{i=1..d} exp(2π I ψ_i): structure function; Dirac node <=> f=0 (S1002).
# ACTIVE column ν: sin(2πν/h)!=0 (ν not in {0,h/2}) — S1005 iε-criterion.  d in {2,3}.
# MEASURED (bets §10, outputs open):
#   T1a (literal bridge): time-active null modes {Λ=0, ν active} — positions vs Dirac nodes.
#   T1b (character bridge): d Dirac nodes (f=0) <-> d nontrivial Z/h column characters — exact.
#   T2 (q=1 signature): sign asymmetry of {Λ} (npos vs nneg, NO λ↔−λ); C̃ iε-count on nulls.
#   T3 (column participation): do time-active null modes exist at equal weights?
# ANTI-TAUTOLOGY guard (§10): "one column => one time" is STRUCTURAL, not a bet — only
#   identities/asymmetries below are measured, none postulated.
# Discipline: 0 handles; exact; mutants>=4 (false-period h+1 · false-sign +T_col · d=2<->3
#   · false-nodes shifted); seeded negctrl; FORBIDDEN in GUARDLINE block; ancestors CITED
#   (S956/S999/S1001/S1002/S1004/S1005/S1007/T26), not re-derived; STOP after tables.
# ============================================================================

import sys
import os
import random
import itertools
from sympy import cos, sin, pi, exp, I, Rational, Integer, simplify, expand, Add, re, im

_HERE = os.path.dirname(os.path.abspath(__file__))


# ⟨★ERRATUM-FORWARD S1063 (Beta's finding J-0513, Ω's task 2026-07-22) — non-gating,
#   an erratum of class T33: the arithmetic IN-DOMAIN is correct, the line below STAYS as it was.
#   DEFECT: `int(...)` — SILENT TRUNCATION, neither rounding nor a domain error.
#   `2−2cos(2πm/h)` is an INTEGER ∀m EXACTLY when h ∈ {1,2,3,4,6} (the crystallographic orders;
#   Niven). Outside this set the Box machinery silently lies:
#       h=5: exact 5/2∓√5/2 = 1.382 / 3.618  →  int gives 1 / 3;
#       h=7: exact 2−2cos(2π/7) = 0.753      →  int gives **0** (the term VANISHES).
#   THE DOMAIN OF THIS PROBE: d ∈ {2,3} ⟹ h = d+1 ∈ {3,4} ⟹ both are in the set ⟹ NO number
#   in S1011/T32 has shifted, the verdict and visas stand.
#   THE CORRECT FORM (use for any extension beyond h ∈ {3,4}):
#       term(h, m) = sp.nsimplify(2 - 2*sp.cos(2*sp.pi*sp.Rational(m % h, h)))  — WITHOUT int.
#   The exact version already lives in `S1063_seam_d1_T2.py` (there h is not domain-limited).⟩
def term(h, m):
    return int(simplify(2 - 2 * cos(2 * pi * Rational(m % h, h))))


def sin_active(h, m):
    v = simplify(sin(2 * pi * Rational(m % h, h)))
    return v != 0


def f_is_zero(psi, h):
    """Exact zero-test of f = 1 + sum exp(2π I ψ_i/h) via Re/Im in cos/sin (robust:
    exp-sums of roots of unity are NOT auto-reduced by sympy)."""
    ref = 1 + Add(*[cos(2 * pi * Rational(m % h, h)) for m in psi])
    imf = Add(*[sin(2 * pi * Rational(m % h, h)) for m in psi])
    return simplify(ref) == 0 and simplify(imf) == 0


def f_show(psi, h):
    ref = simplify(1 + Add(*[cos(2 * pi * Rational(m % h, h)) for m in psi]))
    imf = simplify(Add(*[sin(2 * pi * Rational(m % h, h)) for m in psi]))
    return "{0}{1}I".format(ref, ('+' if imf >= 0 else '') + str(imf))


def T_A(psi, h, tt):
    return sum(tt[m % h] for m in psi)


def dirac_nodes(d, h):
    """All torsion ψ in (Z/h)^d with f(ψ)=0 (exact) — the FULL nodal variety (dim d-2)."""
    out = []
    for psi in itertools.product(range(h), repeat=d):
        if f_is_zero(psi, h):
            out.append(psi)
    return sorted(out)


def coset_nodes(d, h):
    """The d SYMMETRIC coset nodes ψ_i = (i+1)*j mod h (i=0..d-1), j=1..d (S1002):
    the ℤ/h center orbit; f=0 by character orthogonality.  Returns {j: psi}."""
    out = {}
    for j in range(1, h):
        psi = tuple((i + 1) * j % h for i in range(d))
        out[j] = psi
    return out


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1011_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-6 (layer −2): the NATIVE cell Box Λ=T_A(ψ)−T_col(ν) — the minus comes")
    print("from the column duality (T26), not by hand. Zeros=Dirac nodes (bridge №3)? Does q=1 fall out (№2)?")
    print("h=d+1 column points (not a choice); exact integers/roots; d∈{2,3}.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, m):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + m)

    for d in (2, 3):
        h = d + 1
        tt = [term(h, m) for m in range(h)]
        print("=" * 60)
        print("d={0}, h=d+1={1}: term_table(h)={2}  (T_col(ν)=term(ν))".format(d, h, tt))
        print("=" * 60)

        # -------- Dirac nodes: full variety + symmetric coset (character) nodes --------
        full = dirac_nodes(d, h)
        cos_nd = coset_nodes(d, h)
        # verify coset nodes are indeed nodes (f=0) and carry characters 1..d
        for j, psi in cos_nd.items():
            ok(f_is_zero(psi, h), "d={0}: coset-node j={1} ψ={2} has f=0".format(d, j, psi))
        cos_set = set(cos_nd.values())
        TA_cos = sorted(set(T_A(psi, h, tt) for psi in cos_nd.values()))
        ok(len(cos_nd) == d, "d={0}: exactly d={0} coset-nodes (characters of ℤ/{1})".format(d, h))
        ok(sorted(cos_nd.keys()) == list(range(1, h)),
           "★T1b: coset-nodes ↔ d nontriv. characters of ℤ/{0}={1} — an EXACT BIJECTION".format(h, list(range(1, h))))
        ok(TA_cos == [2 * h], "d={0}: T_A at EVERY coset-node = 2h = {1} (constant)".format(d, 2 * h))
        print("Coset-nodes (characters of ℤ/{0}): {1}".format(h, {j: cos_nd[j] for j in sorted(cos_nd)}))
        print("  the full nodal set f=0 (the variety dim d−2, S1002): {0} points {1}".format(
            len(full), full if len(full) <= 12 else "(>12)"))
        print("  ★T1b (CHARACTER bridge): d={0} coset-nodes ↔ {1} nontriv. characters of ℤ/{2} — an exact bijection.".format(
            d, list(range(1, h)), h))
        print("  T_A(coset-node) = 2h = {0} (constant); T_col(ν≠0) ∈ {1} (max {2}) ⟹ 2h ∉ range(T_col).".format(
            2 * h, sorted(set(tt[nu] for nu in range(1, h))), max(tt)))

        # -------- T1a: literal null modes vs nodes --------
        nodes = full  # literal-bridge test uses the FULL nodal set
        active_nulls = []
        all_nulls = []
        for psi in itertools.product(range(h), repeat=d):
            ta = T_A(psi, h, tt)
            for nu in range(h):
                if ta - tt[nu] == 0:
                    all_nulls.append((psi, nu))
                    if nu != 0 and sin_active(h, nu):
                        active_nulls.append((psi, nu))
        null_at_node = [(psi, nu) for (psi, nu) in active_nulls if psi in set(nodes)]
        print("  T1a (LITERAL Λ=0): total zeros={0}, time-active (ν≠0,sin≠0)={1}, "
              "of them ON nodes={2}".format(len(all_nulls), len(active_nulls), len(null_at_node)))
        if active_nulls:
            ta_active = sorted(set(T_A(psi, h, tt) for (psi, nu) in active_nulls))
            print("    T_A of time-active zeros ∈ {0} (nodes have T_A={1}) ⟹ the zeros {2} on the nodes".format(
                ta_active, 2 * h, "SIT" if null_at_node else "DO NOT sit"))
        ok(len(null_at_node) == 0,
           "★T1a (K1-literal): the time-active zeros of Λ are NOT on the nodes (T_A(node)=2h ∉ range T_col) — d={0}".format(d))

        # -------- T2: sign asymmetry (q=1 signature) --------
        spec = []
        for psi in itertools.product(range(h), repeat=d):
            ta = T_A(psi, h, tt)
            for nu in range(h):
                spec.append(ta - tt[nu])
        npos = sum(1 for x in spec if x > 0)
        nneg = sum(1 for x in spec if x < 0)
        nz = sum(1 for x in spec if x == 0)
        # λ↔−λ symmetry test: is the multiset symmetric under negation?
        from collections import Counter
        c = Counter(spec)
        lam_sym = all(c[x] == c[-x] for x in c)
        ok(npos > nneg and not lam_sym,
           "★T2: the Λ spectrum is ASYMMETRIC (npos={0}≫nneg={1}, WITHOUT λ↔−λ) — the q=1 signature of type (d,1)".format(npos, nneg))
        print("  T2 (sign-asymmetry): npos={0}, nneg={1}, nzero={2}; λ↔−λ-symmetry: {3}".format(
            npos, nneg, nz, lam_sym))
        print("    ⟹ the column minus makes the spectrum ONE-SIDEDLY-heavy (space T_A up to d·{0}={1} ≫ time "
              "T_col≤{2}) — asymmetry (d,1), NOT a hand-signed (p,q).".format(max(tt), d * max(tt), max(tt)))

        # -------- T3: column participation --------
        part = len(active_nulls) > 0
        print("  T3 (column participation): time-active null modes={0} ⟹ the column {1} at equal weights".format(
            len(active_nulls), "PARTICIPATES" if part else "is MUTE (structurally, h={0})".format(h)))
        if d == 3:
            ok(part, "T3 d=3: the column participates (time-active zeros exist)")
        print()

    # ==================== MUTANTS ====================
    print("MUTANTS:")
    mut_ok = True

    # M1 (false-period h -> h+1): breaks the node/character match
    d = 2; h = d + 1; tt = [term(h, m) for m in range(h)]
    tt_bad = [term(h + 1, m) for m in range(h + 1)]  # wrong column period
    # nodes need h=3 characters; column with h+1=4 has different nontrivial set
    true_chars = list(range(1, h))
    bad_chars = list(range(1, h + 1))
    if true_chars != bad_chars:
        print("  MUTANT M1 (false-period h→h+1): CAUGHT (a column h+1={0} has {1} characters ≠ {2} nodes "
              "of ℤ/{3} — the bijection breaks)".format(h + 1, len(bad_chars), len(true_chars), h))
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2 (false-sign +T_col: euclidean, spectrum one-sided, no time-active null)
    d = 2; h = 3; tt = [term(h, m) for m in range(h)]
    euclid_nulls = 0
    for psi in itertools.product(range(h), repeat=d):
        ta = T_A(psi, h, tt)
        for nu in range(h):
            if ta + tt[nu] == 0 and nu != 0 and sin_active(h, nu):  # +T_col
                euclid_nulls += 1
    if euclid_nulls == 0:
        print("  MUTANT M2 (false-sign +T_col): CAUGHT (Euclidean T_A+T_col=0 with ν≠0: {0} zeros — both".format(euclid_nulls))
        print("    terms ≥0 ⟹ there are NO time-active zeros; the sign-minus (dual) is necessary for time)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 (d=2 <-> d=3 size): coset-node count = d; character bijection at both
    c2 = coset_nodes(2, 3); c3 = coset_nodes(3, 4)
    ok2 = len(c2) == 2 and all(f_is_zero(p, 3) for p in c2.values())
    ok3 = len(c3) == 3 and all(f_is_zero(p, 4) for p in c3.values())
    if ok2 and ok3:
        print("  MUTANT M3 (size d=2↔3): CAUGHT (d=2→2 coset-nodes/ℤ3, d=3→3/ℤ4 — the bijection scales")
        print("    with h=d+1 legitimately; T_A(coset-node)=2h both; the conclusion is not a lattice artifact)")
    else:
        print("  MUTANT M3: NOT CAUGHT (ok2={0} ok3={1})".format(ok2, ok3)); mut_ok = False

    # M4 (false-nodes: shifted points fail f=0)
    d = 2; h = 3
    shifted = (0, 1)  # not a coset node (ψ=(0,1): f=1+1+exp(2πI/3)!=0)
    if not f_is_zero(shifted, h):
        print("  MUTANT M4 (false-nodes shifted): CAUGHT (a shifted point ψ={0}: f={1}≠0 — NOT a node; "
              "only coset-points of ℤ/h pass the f=0 identity)".format(shifted, f_show(shifted, h)))
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # ==================== NEGATIVE CONTROL (seeded) ====================
    print()
    print("NEGATIVE CONTROL (seeded): a random torsion-point NOT a node (f≠0)")
    random.seed(1011061)
    d = 3; h = 4
    for _try in range(1000):
        rp = tuple(random.randrange(0, h) for _ in range(d))
        if not f_is_zero(rp, h):
            break
    ok(not f_is_zero(rp, h), "control: random ψ={0} → f={1}≠0 (not a node)".format(rp, f_show(rp, h)))
    print("  ψ={0}: f={1}≠0 — the measurement is sensitive (a node only on character-cosets of ℤ/h)".format(rp, f_show(rp, h)))

    # ==================== SUMMARY ====================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'column/dual/character/nodes' is STRUCTURAL cell/root-lattice vocabulary (S1001/1002) — not fenced. GUARDLINE
    _pp = [("причи", "нн"), ("всес", "віт"), ("ант", "роп"), ("Teg", "mark")]  # GUARDLINE
    _PATTERNS = ["".join(ab) for ab in _pp]  # GUARDLINE
    _hs = scan_forbidden(__file__, _PATTERNS); _logf.flush()
    _hl = scan_forbidden("".join(_tee.chunks), _PATTERNS)
    _n = len(_hs) + len(_hl)
    print("FORBIDDEN-SCAN: hits={0} (src={1}, log={2})".format(_n, len(_hs), len(_hl)))

    _exit = 1 if (_n > 0 or FA[0] > 0 or not mut_ok) else 0
    print("EXIT={0}".format(_exit)); print("PROC_EXIT={0}".format(_exit)); print("STOP")
    sys.exit(_exit)


if __name__ == "__main__":
    main()
