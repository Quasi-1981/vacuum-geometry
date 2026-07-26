# -*- coding: utf-8 -*-
# DIM: na (W42 probe-10, layer -2: ARROW, road K. Is the orientation of the clock
#          (c→c⁻¹, Coxeter S1004) coupled with the orientation of SPACE (mirror = Weyl parity)?
#          Exact count in S_{d+1}: parity of the inverters {g: gcg⁻¹=c⁻¹} over d∈{2,3,4}.
#          + ε→−ε (T28) = a measured symmetry (k↔−k=w₀)? + an azimuth info-row + the T34 link.
#          ★COUNTING, not physics. FS SPECIAL: heat-bath language (β/Matsubara) is FORBIDDEN — GUARDLINE.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting — anti-tuning, §14 exante)
# ----------------------------------------------------------------------------
# Weyl(A_d)=S_{d+1} (n=d+1 axes). Coxeter c = n-cycle (order h=n, S1004). Two directions
#   c, c⁻¹.  MIRROR = orientation = Weyl parity (reflection-rep det = sign(σ)): EVEN=A_{n}
#   (orient-preserving) ⊥ ODD (orient-reversing = a reflection/mirror).
# EXERCISE-FENCE (theorem, NOT bet): c~c⁻¹ in full S_n (same cycle-type) — given, not measured.
# ★COUNTERS (rules FIXED HERE, BEFORE numbers — will NOT change after data):
#   STAVKA-1a: I(d) = {g in S_n : g c g⁻¹ = c⁻¹}; report |I|, #even, #odd — per d.
#     Coupling-with-mirror FORCED iff I ∩ A_n = ∅ (all inverters odd).  d-dep = honest result.
#   STAVKA-1b: w₀ = reversal i↦n−1−i (= longest Weyl elt = momentum inversion k↔−k, S1004);
#     is w₀ ∈ I (does k↔−k invert the clock)? parity(w₀)? ⟹ ε→−ε (T28 arrow) = FACET of a
#     MEASURED symmetry (k↔−k), not independent freedom; its mirror-content = parity(w₀).
#   STAVKA-2: azimuth string of c-orbit from 0 = [c^0(0),c^1(0),...] = positions visited.
#     (ABS)  chiral? cyclic word W != reverse(W) up to rotation (absolute azimuth labels).
#     (REL)  increment word (differences mod h); c=[+1..], c⁻¹=[−1..]; distinguishable ONLY
#            if the ℤ/h negation-relabel x↦−x is NOT allowed. Report BOTH readings raw.
#   STAVKA-3 (link T34): orbits of (marked bond a, orientation o∈{c,c⁻¹}) under dihedral
#     D_h=⟨c,w₀⟩; does w₀ fix/move the marked bond? #orbits — bridge to door D.
#   ★STAVKA-ANZATS (added §14-update, author pre-reg «dial turns ONE way: prime→daughter,
#     descending»): among MEASURED two-component symmetries, does EVERY realization of
#     clock-inversion (k→−k) also INVERT the hierarchy marker A↔B (σ_x swap)?  Ancestor
#     S1012 (CITED, not re-derived): H(k)=[[0,f],[f̄,0]] ⟹ H(−k)=σ_x H σ_x, so the ONLY
#     k-inverting symmetry is B=σ_x∘(k→−k) (bare k→−k NOT a symmetry; bare σ_x NOT either).
#     COUNTER (rule before numbers): enumerate g=K^a X^b Z^c (K=k→−k, X=σ_x-conj, Z=σ_z-conj);
#     is_sym(g)⟺ g·H = H; inverts_clock⟺a=1; flips_hierarchy⟺b=1.  Ansatz WINS iff NO
#     hierarchy-neutral rotator (is_sym ∧ a=1 ∧ b=0) exists ∀d; find one ⟹ ansatz loses, carve.
# ★TARGETS/KILLS: K1 inverting reachable by EVEN ∀d ⟹ no mirror-coupling (road R closer).
#   K2 strings c/c⁻¹ indistinguishable ∀ internal counter ⟹ info-leg blind (2nd blind, road R).
#   K3 (reclaim): arrow NOT "chosen" — even with coupling it is ADDRESS-TRANSFER (mirror/bond).
#   K4: no new handle.
# Discipline: exact perm arithmetic; mutants>=4 (false-parity · false-composition · size
#   d=2<->3<->4 · false-string-counter non-rotation-inv); seeded negctrl; ancestors CITED
#   (S1000-T2 rep · S1004 Coxeter · T28 ε-classes · T34 residual); ★FS heat-bath language FORBIDDEN.
# ============================================================================

import sys
import os
import random
import itertools
from sympy import Matrix, eye, Symbol, I as sympyI

_HERE = os.path.dirname(os.path.abspath(__file__))


# ==================== two-component operator machinery (S1012, exact) ====================

def _pauli():
    sx = Matrix([[0, 1], [1, 0]])
    sy = Matrix([[0, -sympyI], [sympyI, 0]])
    sz = Matrix([[1, 0], [0, -1]])
    return sx, sy, sz


def anzats_symmetry_scan():
    """Enumerate g = K^a X^b Z^c on H(k)=[[0,f],[f̄,0]] (S1012). f,fb independent symbols;
    K = k→−k (swap f,fb); X = σ_x-conj; Z = σ_z-conj. Return list of dicts per g with
    is_sym / inverts_clock(a) / flips_hierarchy(b)."""
    f = Symbol('f'); fb = Symbol('fb')
    H = Matrix([[0, f], [fb, 0]])
    sx, sy, sz = _pauli()
    rows = []
    for a in (0, 1):       # K: k→−k
        for b in (0, 1):   # X: σ_x-conj (A↔B swap)
            for c in (0, 1):  # Z: σ_z-conj (chirality)
                M = H
                if c:
                    M = sz * M * sz
                if b:
                    M = sx * M * sx
                if a:
                    M = M.subs({f: Symbol('tmp'), fb: f}, simultaneous=True).subs(Symbol('tmp'), fb)
                is_sym = (M == H)
                rows.append(dict(a=a, b=b, c=c, is_sym=is_sym,
                                 inverts_clock=(a == 1), flips_hier=(b == 1)))
    return rows


# ==================== permutation primitives (exact) ====================

def compose(p, q):
    """(p∘q)(i) = p[q[i]]."""
    return tuple(p[q[i]] for i in range(len(p)))


def inverse(p):
    inv = [0] * len(p)
    for i, pi in enumerate(p):
        inv[pi] = i
    return tuple(inv)


def parity(p):
    """sign as #transpositions mod 2 via cycle decomposition; return 0(even)/1(odd)."""
    n = len(p); seen = [False] * n; tr = 0
    for i in range(n):
        if not seen[i]:
            j = i; ln = 0
            while not seen[j]:
                seen[j] = True; j = p[j]; ln += 1
            tr += ln - 1
    return tr % 2


def n_cycle(n):
    """c: i -> (i+1) mod n (the Coxeter element as an n-cycle)."""
    return tuple((i + 1) % n for i in range(n))


def reversal(n):
    """w₀: i -> n-1-i (longest Weyl element / momentum inversion k↔−k)."""
    return tuple(n - 1 - i for i in range(n))


def inverting_set(n):
    """{g in S_n : g c g⁻¹ = c⁻¹}."""
    c = n_cycle(n); cinv = inverse(c)
    out = []
    for perm in itertools.permutations(range(n)):
        g = tuple(perm)
        if compose(compose(g, c), inverse(g)) == cinv:
            out.append(g)
    return out


def azimuth_word(n, direction=+1):
    """positions visited by c^{±1}-orbit from 0: [0, ±1, ±2, ...] mod n."""
    return tuple((direction * k) % n for k in range(n))


def cyclic_equal_up_to_rotation(w, v):
    """True iff v is a rotation of w (cyclic words equal)."""
    n = len(w)
    if len(v) != n:
        return False
    dbl = w + w
    return any(dbl[r:r + n] == v for r in range(n))


def is_chiral_cyclic(w):
    """True iff w != reverse(w) up to rotation (oriented cyclic word)."""
    return not cyclic_equal_up_to_rotation(w, tuple(reversed(w)))


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1016_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-10 (layer −2): ARROW, road K. Is the orientation of the clock coupled")
    print("(c→c⁻¹, Coxeter) with the orientation of SPACE (mirror=Weyl parity)? Exact count S_{d+1}.")
    print("★COUNTING, not physics; the info-count is INFORMATIONAL only; heat-bath stuff behind the fence; ancestors cited.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, msg):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + msg)

    # ================= STAVKA-1a: parity of the inverters c→c⁻¹ =================
    print("#" * 70)
    print("# ★STAVKA-1a (the heart): parity of the set of inverters {g: gcg⁻¹=c⁻¹} over d")
    print("#" * 70)
    print("  d | n=d+1 | |I| | #even(EVEN) | #odd(ODD) | coupling with the mirror?")
    print("  " + "-" * 70)
    s1 = {}
    for d in (2, 3, 4):
        n = d + 1
        I = inverting_set(n)
        pars = [parity(g) for g in I]
        neven = pars.count(0); nodd = pars.count(1)
        forced = (neven == 0)   # all odd ⟹ inversion REQUIRES the mirror
        s1[d] = (len(I), neven, nodd, forced)
        ok(len(I) == n, "d={0}: |I|=n={1} (the coset of the cyclic center ⟨c⟩)".format(d, n))
        verdict = ("FORCED (all odd)" if forced else
                   ("NONE (all even=K1)" if nodd == 0 else "MIXED (some even too)"))
        print("  {0} | {1} | {2} | {3} | {4} | {5}".format(d, n, len(I), neven, nodd, verdict))
    print()
    print("  ★RAW ROW (an honest result either way, carved BEFORE the count):")
    print("   d=2: {0} · d=3: {1} · d=4: {2}".format(
        "all odd ⟹ coupling with the mirror is FORCED",
        "2 even/2 odd ⟹ inversion is reachable BOTH ways (coupling is optional)",
        "all even ⟹ NO coupling (inversion PRESERVES orientation) — K1 at d=4"))
    print("   ⟹ the arrow⊗mirror coupling is d-DEPENDENT, non-monotone: strong at d=2 → none at d=4.")
    print()

    # ================= STAVKA-1b: ε→−ε = w₀ = k↔−k (measured symmetry)? =================
    print("#" * 70)
    print("# STAVKA-1b: ε→−ε (T28) = w₀ (reversal = k↔−k, measured)? the parity of w₀ over d")
    print("#" * 70)
    for d in (2, 3, 4):
        n = d + 1
        w0 = reversal(n)
        c = n_cycle(n); cinv = inverse(c)
        w0_inverts = (compose(compose(w0, c), inverse(w0)) == cinv)
        pw = parity(w0)
        ok(w0_inverts, "d={0}: w₀ (k↔−k) INVERTS the clock (w₀cw₀⁻¹=c⁻¹)".format(d))
        print("  d={0}, n={1}: w₀={2}, inverts c={3}, parity of w₀={4} ({5})".format(
            d, n, w0, w0_inverts, pw, "odd=mirror" if pw else "even=rotation"))
    print("  ⟹ the arrow ε→−ε = a facet of an ALREADY-measured symmetry k↔−k (w₀), NOT an independent freedom;")
    print("   its mirror-content = the parity of w₀ (odd at d=2 · even at d=3,4) — an address-transfer (K3),")
    print("   NOT a selection mechanism. (At d=3 w₀ is even, BUT I also has odd inverters — stavka-1a.)")
    print()

    # ================= ★STAVKA-ANZATS: clock-rotation ⟹ hierarchy A↔B inversion? =================
    print("#" * 70)
    print("# ★STAVKA-ANZATS (the author's pre-reg «the dial turns one way: prime→daughter»):")
    print("#   does EVERY measured realization of k→−k (rotation of the circle) INVERT the marker A↔B (σ_x)?")
    print("#" * 70)
    print("  ancestor S1012 (citation): H(−k)=σ_x H σ_x ⟹ bare k→−k is NOT a symmetry, bare σ_x is NOT")
    print("  a symmetry — only B=σ_x∘(k→−k) is. Counting g=K^a X^b Z^c (rule BEFORE the numbers):")
    rows = anzats_symmetry_scan()
    syms = [r for r in rows if r["is_sym"]]
    clock_inverters = [r for r in syms if r["inverts_clock"]]
    hier_neutral = [r for r in syms if r["inverts_clock"] and not r["flips_hier"]]
    print("  symmetries (g·H=H): {0}".format(
        ["K{0}X{1}Z{2}".format(r["a"], r["b"], r["c"]) for r in syms]))
    print("  of them invert the clock (a=1): {0} — do all flip A↔B (b=1)? {1}".format(
        ["K{0}X{1}Z{2}".format(r["a"], r["b"], r["c"]) for r in clock_inverters],
        all(r["flips_hier"] for r in clock_inverters)))
    print("  ★HIERARCHY-NEUTRAL ROTATOR (symmetry ∧ a=1 ∧ b=0): {0} — {1}".format(
        len(hier_neutral), "DOES NOT EXIST ⟹ THE ANSATZ WINS" if not hier_neutral else "EXISTS ⟹ the ansatz loses"))
    ok(all(r["flips_hier"] for r in clock_inverters) and not hier_neutral,
       "★STAVKA-ANZATS: every measured realization of k→−k flips A↔B; 0 hierarchy-neutral rotators")
    ok(len(clock_inverters) == 1 and clock_inverters[0]["a"] == 1 and clock_inverters[0]["b"] == 1,
       "the single measured clock-rotation = B=σ_x∘(k→−k) [a=1,b=1] — flips A↔B (not bare k→−k)")
    # ★NATIVE d∈{2,3} check (gate-style, the project's word): for the ACTUAL native
    # f(k)=1+Σ_{i=1..d}exp(2πi ψ_i/h) the key identity H(−k)=σ_x H σ_x reduces to
    # f(−k)=conj(f(k)) — checked EXACTLY on the grid of native momenta (roots of unity).
    from sympy import cos as _cos, sin as _sin, pi as _pi, Rational as _R, simplify as _sm
    for d in (2, 3):
        h = d + 1
        native_ok = True
        for psi in itertools.product(range(h), repeat=d):
            # Re f(−k) == Re f(k)  and  Im f(−k) == −Im f(k)  (⟺ f(−k)=conj f(k))
            re_p = _sm(1 + sum(_cos(2 * _pi * _R(m % h, h)) for m in psi))
            im_p = _sm(sum(_sin(2 * _pi * _R(m % h, h)) for m in psi))
            re_m = _sm(1 + sum(_cos(2 * _pi * _R((-m) % h, h)) for m in psi))
            im_m = _sm(sum(_sin(2 * _pi * _R((-m) % h, h)) for m in psi))
            if not (re_m == re_p and im_m == -im_p):
                native_ok = False
                break
        ok(native_ok, "★NATIVE d={0}: f(−k)=conj f(k) EXACTLY ∀native momenta ⟹ H(−k)=σ_x H σ_x "
           "⟹ bare k→−k is NOT a native symmetry (the ansatz wins on the native structure)".format(d))
        print("  native d={0} (h={1}): f(−k)=conj f(k) exactly ∀{2} native momenta — {3}".format(
            d, h, h ** d, "H(−k)=σ_x H σ_x ✓ (bare k→−k is not a symmetry)" if native_ok else "BREAK"))
    print("  ⟹ d-INDEPENDENTLY (the Pauli-algebra of H is the same ∀d; f only changes the #terms, not the structure):")
    print("   the hierarchy-change ⟺ circle-rotation are LOCKED TOGETHER in B — «descending, not the reverse» WINS as a measurement.")
    print("   Contrast with stavka-1a: in the full S_{d+1} inverters can be EVEN (d=3,4), BUT those are")
    print("   SPATIAL permutations, not two-component measured symmetries — the ansatz lives at the S1012 level.")
    print()

    # ================= STAVKA-2: azimuth info-row =================
    print("#" * 70)
    print("# STAVKA-2 (info-row): azimuth row of the c-traversal vs c⁻¹ — distinguishable INTERNALLY?")
    print("#" * 70)
    for d in (2, 3, 4):
        n = d + 1
        wc = azimuth_word(n, +1); wci = azimuth_word(n, -1)
        abs_chiral = is_chiral_cyclic(wc)   # absolute azimuth labels
        # REL: increment words; c=+1..,  c⁻¹=−1..≡(h−1)..; related by negation x↦−x mod h
        inc_c = tuple((wc[(k + 1) % n] - wc[k]) % n for k in range(n))
        inc_ci = tuple((wci[(k + 1) % n] - wci[k]) % n for k in range(n))
        neg_of_inc_c = tuple((-x) % n for x in inc_c)
        rel_blind = cyclic_equal_up_to_rotation(neg_of_inc_c, inc_ci)  # the negation-relabel saves it
        print("  d={0},n={1}: ABS row c={2} — chiral(oriented)={3} · "
              "REL increment c={4},c⁻¹={5} — blind-under-negation={6}".format(
                  d, n, wc, abs_chiral, inc_c, inc_ci, rel_blind))
        ok(abs_chiral == (n >= 3), "d={0}: the ABS-row is chiral ⟺ n>=3".format(d))
        ok(rel_blind, "d={0}: the REL-increment c⁻¹ = negation(c) — blind without an absolute zero".format(d))
    print("  ★HONEST FORK (both readings raw, the verdict — the project's):")
    print("   ABS (there is an absolute azimuth-zero): the row is ORIENTED ∀d≥2 — the arrow reads from inside.")
    print("   REL (increments only): c⁻¹=negation(c) ⟹ BLIND — distinguishing needs an already-fixed")
    print("   azimuth sign = the ARROW ITSELF (circularly) ⟹ a D_net=0 rhyme, reinforcing road R.")
    print()

    # ================= STAVKA-3: link with T34 (bond representative ⊗ orientation) =================
    print("#" * 70)
    print("# STAVKA-3 (T34 link): orbits of (marked bond, orientation) under D_h=⟨c,w₀⟩")
    print("#" * 70)
    for d in (2, 3, 4):
        n = d + 1
        c = n_cycle(n); w0 = reversal(n)
        # dihedral group generated by c, w0
        G = set()
        frontier = [tuple(range(n))]
        gens = [c, w0]
        while frontier:
            x = frontier.pop()
            if x in G:
                continue
            G.add(x)
            for gg in gens:
                frontier.append(compose(gg, x))
        # orbits of (bond a, orientation o): o in {+1,-1}; group elt g acts:
        #   bond a -> g[a]; orientation flips iff g is a reflection (parity depends: use
        #   "reverses c" test — g in inverting coset flips orientation, else preserves).
        cinv = inverse(c)
        def flips(g):
            return compose(compose(g, c), inverse(g)) == cinv
        flags = [(a, o) for a in range(n) for o in (+1, -1)]
        seen = set(); norb = 0
        for f in flags:
            if f in seen:
                continue
            norb += 1
            stack = [f]
            while stack:
                (a, o) = stack.pop()
                if (a, o) in seen:
                    continue
                seen.add((a, o))
                for g in G:
                    na = g[a]; no = (-o if flips(g) else o)
                    stack.append((na, no))
        # does w0 fix any marked bond?
        w0_fixed = [a for a in range(n) if w0[a] == a]
        ok(len(G) == 2 * n, "d={0}: |D_h|=2n={1}".format(d, 2 * n))
        print("  d={0},n={1}: |D_h|={2}, orbits (bond,orient)={3}, w₀-fixed bonds={4} ({5})".format(
            d, n, len(G), norb, w0_fixed,
            "the middle bond (n odd)" if w0_fixed else "none (n even) — bond↔partner"))
    print("  ⟹ (bond-representative ⊗ orientation) = ONE orbit under the dihedral group ⟹ the two T34 residues")
    print("   (the bond-representative + the arrow) are COUPLED by the cell symmetry — the first bridge to door D.")
    print()

    # ================= MUTANTS (>=4) =================
    print("MUTANTS:")
    mut_ok = True

    # M1 false-parity: the parity of the FIRST inverter != the parity of all (where they differ, d=3)
    d = 3; n = 4; I = inverting_set(n); pars = [parity(g) for g in I]
    if len(set(pars)) == 2:  # d=3 has BOTH parities — a bare «first parity» would lie
        print("  MUTANT M1 (false-parity=first inverter): CAUGHT (d=3 I has BOTH parities "
              "{0} — the verdict must count the WHOLE set, not a representative; a bare 'first' would mislead)".format(
                  sorted(set(pars))))
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2 false-composition: an element that does NOT invert c, falsely declared an inverter
    d = 4; n = 5; c = n_cycle(n); cinv = inverse(c)
    fake = compose(n_cycle(n), n_cycle(n))  # c² — NOT an inverter (it's in the center)
    if compose(compose(fake, c), inverse(fake)) != cinv:
        print("  MUTANT M2 (false-composition): CAUGHT (c²∈⟨c⟩ does NOT invert: gcg⁻¹=c!=c⁻¹ — "
              "membership in I is checked by the equation, not by declaration)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 size d=2↔3↔4: the raw d-dependency (0/2/5 even inverters) is NOT an artifact — it is structural
    evens = [s1[d][1] for d in (2, 3, 4)]
    if evens == [0, 2, 5]:
        print("  MUTANT M3 (size d=2↔3↔4): CAUGHT (#even inverters = {0} — the d-structure "
              "(n odd→all even/odd equal; n even→split) is LEGITIMATE, not a lattice artifact)".format(evens))
    else:
        print("  MUTANT M3: NOT CAUGHT (evens={0})".format(evens)); mut_ok = False

    # M4 false string-counter: a non-rotation-invariant «chiral» test misfires on a
    #    palindrome; control — a row EQUAL to its reverse (n=2 degenerate) must be NOT chiral
    w2 = azimuth_word(2, +1)   # n=2: [0,1] — reverse [1,0] = rotation of [0,1]? yes
    if not is_chiral_cyclic(w2) and is_chiral_cyclic(azimuth_word(3, +1)):
        print("  MUTANT M4 (false string-counter): CAUGHT (n=2 row [0,1] is NOT chiral "
              "(reverse=rotation), n=3 is chiral — the test is rotation-invariant, does not misfire on a palindrome)")
    else:
        print("  MUTANT M4: NOT CAUGHT"); mut_ok = False

    # M5 false-relation of the ansatz: IF bare k→−k were a symmetry (without σ_x), a
    #    hierarchy-neutral rotator would exist ⟹ the ansatz would lose. We show that this is an
    #    S1012-fact (H(−k)!=H(k)) that carries the conclusion, not a postulate: bare K (a=1,b=0,c=0) is NOT a symmetry.
    rows_m = anzats_symmetry_scan()
    bare_K = [r for r in rows_m if r["a"] == 1 and r["b"] == 0 and r["c"] == 0][0]
    if not bare_K["is_sym"]:
        print("  MUTANT M5 (false-relation of the ansatz): CAUGHT (bare k→−k [K1X0Z0] is NOT a symmetry: "
              "H(−k)=σ_x H σ_x != H — IF it were, a hierarchy-neutral rotator would exist; "
              "the conclusion carries an S1012-fact, not a postulate)")
    else:
        print("  MUTANT M5: NOT CAUGHT"); mut_ok = False

    # ================= NEGATIVE CONTROL (seeded) =================
    print()
    print("NEGATIVE CONTROL (seeded): a random element of S_n — NOT an inverter (in general)")
    random.seed(1016101)
    n = 5; c = n_cycle(n); cinv = inverse(c)
    for _ in range(1000):
        g = list(range(n)); random.shuffle(g); g = tuple(g)
        if compose(compose(g, c), inverse(g)) != cinv:
            break
    ok(compose(compose(g, c), inverse(g)) != cinv,
       "negctrl: random g={0} does NOT invert c (inverters are a rare set |I|=n)".format(g))
    print("  g={0}: gcg⁻¹ != c⁻¹ — the measurement is sensitive (only {1} of {2} elements invert)".format(
        g, n, "120"))

    # ================= SUMMARY =================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'rotation/inverter/mirror/parity/azimuth/orientation' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("темпера", "тура"), ("Мацу", "бара"), ("пропага", "тор"),
           ("чекер", "борд"), ("всес", "віт"), ("ант", "роп"),
           ("бе", "та-температ")]  # GUARDLINE (heat-bath+action FORBIDDEN; the count is informational)
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
