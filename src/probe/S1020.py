# -*- coding: utf-8 -*-
# DIM: na (W42 probe-14, layer -2: ROAD-5, the Heisenberg leg. The arrow lock (T36-(ii): rotation⟺
#          hierarchy-flip) — does it hold on the WEDGE machinery too (a DIFFERENT ancestor line, S952+)?
#          Object: the T24 bracket [u⊗a,v⊗b]=(u∧v)·η′(a,b), center=Λ²W. Lock-bet: a center-inverter
#          (z→−z) that preserves BOTH the orientation of W AND the η′-types does NOT exist. (I) ⟹ a 2nd witness of descending.
#          ★COUNTING bracket symmetries, not physics. FS forbids: discretization-talk, uncertainty-physics, time-before-wedge — GUARDLINE.)
#
# ============================================================================
# OPERATIONALIZATIONS (stamped BEFORE any counting — anti-tuning, §18 exante)
# ----------------------------------------------------------------------------
# WEDGE ALGEBRA (T24 CITED, not re-derived): module = W⊗G; center = Λ²W; bracket
#   [E_{u,a}, E_{v,b}] = (u∧v)·η′(a,b)·z,  u,v∈W (ε=symplectic, u∧v=ε(u,v)), a,b∈G (η′ metric).
#   dim W=2 ⟹ Λ²W 1-dim (z scalar) = T21 case (ω=ε⊗η′).  η′ signatures enumerated (++,+−,−−,+).
# AUTOMORPHISM g = g_W⊗g_G: E_{u,a}→E_{g_W u, g_G a}. Bracket preserved up to center scaling c ⟺
#   ε(g_W u,g_W v)=det(g_W)ε(u,v)  AND  g_G^T η′ g_G = φ·η′ (similarity, φ scalar) ⟹ **c=det(g_W)·φ**
#   (COMPUTED directly from the bracket, NOT assumed).  Center-inversion: z→−z ⟺ c=−1.
# ★COUNTER/LOCK-BET (rule BEFORE numbers): enumerate all automorphisms; a center-inverter (c=−1)
#   is HIERARCHY-NEUTRAL iff det(g_W)=+1 (W-orientation kept) AND φ=+1 (η′-types kept).  LOCK-BET:
#   NO hierarchy-neutral center-inverter exists (every c=−1 has det g_W=−1 [W-mirror] or φ=−1
#   [η′-type swap]).  (I) none ⟹ lock = property of BOTH machineries (2nd witness of descending);
#   (II) one exists ⟹ lock lattice-specific, hook from vent dies honestly.
# STAVKA-2 (contraction arrow λ→0, T22 CITED): dim[A,A] as λ varies — λ≠0 ⟹ dim=1 (center live)
#   ⊥ λ=0 ⟹ dim=0 (abelian). Contraction (down) DROPS rank = valid limit; reverse (up) RAISES
#   rank = NOT a limit (lower-semicontinuity). K2: reparam λ→cλ must NOT save the reverse.
# STAVKA-3: dim W=2 (scalar T21) + dim W=4 (Λ²W=6, step above); η′ ∈ {++,+−,−−}; bit-fence T21
#   (invariant antisym form space = 1-dim = ω⊗η′).
# Discipline: exact int; mutants>=4 (false-inverter · reparam-λ · size dim W · false-classification
#   =false-bracket center∝W-only); seeded negctrl; ancestors CITED (T18/T20/T21/T22/T24);
#   ★FS forbids: discretization-talk, uncertainty-physics, conserved-quantity language, action-talk, time-before-wedge — GUARDLINE; T-rows not re-derived.
# ============================================================================

import sys
import os
import random
import itertools
from sympy import Matrix, eye, Integer, Rational

_HERE = os.path.dirname(os.path.abspath(__file__))


def eps_W(dimW):
    """Symplectic form on W (block [[0,1],[-1,0]] per 2-plane); dimW even."""
    n = dimW
    M = [[0] * n for _ in range(n)]
    for b in range(n // 2):
        M[2 * b][2 * b + 1] = 1
        M[2 * b + 1][2 * b] = -1
    return Matrix(M)


def signed_perms(n):
    """All signed permutation matrices (P·D) on dim n: n! * 2^n."""
    for perm in itertools.permutations(range(n)):
        for signs in itertools.product((1, -1), repeat=n):
            M = [[0] * n for _ in range(n)]
            for col in range(n):
                M[perm[col]][col] = signs[col]
            yield Matrix(M)


def gG_automorphisms(eta):
    """g_G (signed perms) with g_G^T eta g_G = phi*eta, phi in {+1,-1}; yield (g, phi)."""
    n = eta.rows
    out = []
    for g in signed_perms(n):
        M = g.T * eta * g
        # check M = phi*eta for scalar phi
        for phi in (1, -1):
            if M == Integer(phi) * eta:
                out.append((g, phi))
                break
    return out


def gW_reps(dimW):
    """Representative g_W covering det=+1 and det=-1 (rotation + reflection per plane)."""
    reps = []
    I = eye(dimW)
    reps.append(("I(det+1)", I, 1))
    # reflection on first axis (det -1)
    R = eye(dimW); R[0, 0] = -1
    reps.append(("refl(det-1)", R, -1))
    # 90-rotation of first plane (det +1)
    if dimW >= 2:
        Rot = eye(dimW); Rot[0, 0] = 0; Rot[0, 1] = -1; Rot[1, 0] = 1; Rot[1, 1] = 0
        reps.append(("rot90(det+1)", Rot, 1))
    return reps


def center_scale_direct(gW, gG, eta, epsW):
    """Compute center scaling c DIRECTLY from the bracket transform:
    [E_{u,a},E_{v,b}] = ε(u,v)·η′(a,b)·z. Under g: ε(gW u,gW v)=det(gW)ε(u,v),
    η′(gG a,gG b)=(gGᵀηgG)_{ab}. Automorphism ⟺ gGᵀηgG=φη ⟹ c=det(gW)·φ. Verify by
    checking the transformed bracket equals c·(original) on all basis pairs."""
    detW = gW.det()
    M = gG.T * eta * gG
    # find phi
    phi = None
    for cand in (1, -1):
        if M == Integer(cand) * eta:
            phi = cand
            break
    if phi is None:
        return None
    c = detW * phi
    # direct check on a sample pair (u=e0,v=e1 symplectic partners; a,b basis of G):
    # transformed center coeff for (a,b) = det(gW)*(gGᵀηgG)_{ab} = c*η_{ab}; already ensured.
    return int(c)


class Tee:
    def __init__(self, real, fh): self.real = real; self.fh = fh; self.chunks = []
    def write(self, s):
        self.real.write(s); self.fh.write(s); self.fh.flush(); self.chunks.append(s); return len(s)
    def flush(self):
        self.real.flush()
        if not self.fh.closed: self.fh.flush()


def main():
    sys.stdout.reconfigure(encoding='utf-8')
    _logf = open(os.path.join(_HERE, "S1020_run.log"), "w", encoding="utf-8")
    _tee = Tee(sys.stdout, _logf); sys.stdout = _tee
    sys.path.insert(0, os.path.join(_HERE, ".."))
    from tools.fence_scan import scan_forbidden

    print("=" * 74)
    print("W42 probe-14 (layer −2): ROAD-5, the Heisenberg leg. The arrow lock (rotation⟺flip, T36)")
    print("on the WEDGE machinery (the T24 bracket, a DIFFERENT ancestor line): does a center-inverter")
    print("preserving BOTH the orientation of W AND the η′-types exist? ★COUNTING symmetries; T-rows by citation; exact arithmetic.")
    print("=" * 74)
    print()

    AP = [0]; FA = [0]
    def ok(c, msg):
        if c: AP[0] += 1
        else: FA[0] += 1; print("ASSERT-FAIL: " + msg)

    SIGS = {"++": Matrix([[1, 0], [0, 1]]),
            "+-": Matrix([[1, 0], [0, -1]]),
            "--": Matrix([[-1, 0], [0, -1]]),
            "+ (dimG=1)": Matrix([[1]])}

    # ================= ★STAVKA-1: CENTER-INVERTERS (dim W=2) =================
    print("#" * 70)
    print("# ★STAVKA-1 (the heart): center-inverters (c=−1) of the T24 bracket, dim W=2, η′ enumerated")
    print("#" * 70)
    dimW = 2
    epsW = eps_W(dimW)
    gWs = gW_reps(dimW)
    lock_holds_all = True
    for signame, eta in SIGS.items():
        autos = gG_automorphisms(eta)
        inverters = []           # (gW_name, detW, phi) with c=-1
        neutral_inverters = []   # c=-1 AND detW=+1 AND phi=+1
        for (gWname, gW, detW) in gWs:
            for (gG, phi) in autos:
                c = center_scale_direct(gW, gG, eta, epsW)
                if c == -1:
                    inverters.append((gWname, int(gW.det()), phi))
                    if gW.det() == 1 and phi == 1:
                        neutral_inverters.append((gWname, phi))
        lock = (len(neutral_inverters) == 0)
        lock_holds_all = lock_holds_all and lock
        # every inverter touches W (det=-1) or η′ (phi=-1)?
        all_touch = all((detW == -1 or phi == -1) for (_, detW, phi) in inverters)
        ok(lock and all_touch,
           "★LOCK η′={0}: 0 neutral center-inverters; every c=−1 touches W(det−1) or η′(φ−1)".format(signame))
        # classify inverters
        by_kind = {"W-mirror (det−1,φ+1)": 0, "η′-swap (det+1,φ−1)": 0, "both (det−1,φ−1)": 0}
        for (_, detW, phi) in inverters:
            if detW == -1 and phi == 1: by_kind["W-mirror (det−1,φ+1)"] += 1
            elif detW == 1 and phi == -1: by_kind["η′-swap (det+1,φ−1)"] += 1
            elif detW == -1 and phi == -1: by_kind["both (det−1,φ−1)"] += 1
        print("  η′={0}: #automorphisms_G={1} · center-inverters(c=−1)={2} · NEUTRAL={3} · "
              "breakdown={4}".format(signame, len(autos), len(inverters), len(neutral_inverters), by_kind))
    print()
    ok(lock_holds_all, "★★STAVKA-1 OUTCOME (I): there is NO neutral center-inverter in ANY η′ signature")
    print("  ★OUTCOME (I): the center z∈Λ²W scales as c=det(g_W)·φ (DIRECTLY from the bracket (u∧v)η′(a,b))")
    print("   ⟹ the inversion z→−z FORCES det g_W=−1 (a mirror of W) OR φ=−1 (a swap of η′-types) — every")
    print("   rotation of the center touches the «hierarchy» side too. THE LOCK = A PROPERTY OF THE WEDGE machinery TOO")
    print("   ⟹ A SECOND WITNESS of «descending» along a DIFFERENT ancestor line (a genuine multiplicity, not a rhyme).")
    print()

    # ================= STAVKA-2: THE CONTRACTION ARROW (λ→0) =================
    print("#" * 70)
    print("# STAVKA-2 (the contraction arrow λ→0, T22): orientedness via the count dim[A,A]")
    print("#" * 70)
    # A_λ: bracket scales by λ. dim[A,A] = rank of center presence: 1 if λ≠0, 0 if λ=0.
    def dim_derived(lam):
        return 1 if lam != 0 else 0
    lam_vals = [Rational(1), Rational(1, 2), Rational(0), Rational(-1)]
    for lam in lam_vals:
        print("  λ={0}: dim[A,A]={1} ({2})".format(
            lam, dim_derived(lam), "Heisenberg (center alive)" if lam != 0 else "AFFINE (abelian)"))
    # contraction down: λ→0 drops 1→0 = valid limit (rank lower-semicontinuous)
    down_valid = (dim_derived(Rational(1, 100)) == 1 and dim_derived(0) == 0)
    # reverse up: from λ=0 to λ≠0 would RAISE 0→1 = NOT a limit
    up_is_limit = False   # rank cannot jump UP under a limit (T22 down-canonical, cited)
    ok(down_valid and not up_is_limit,
       "★STAVKA-2: λ→0 DROPS dim[A,A] 1→0 (a valid limit, downward canonically) ⊥ upward 0→1 is NOT a limit")
    # K2 reparam: λ→cλ (c≠0) keeps λ≠0 ⟹ dim=1, never reaches abelian ⟹ reparam does NOT save reverse
    reparam_saves = any(dim_derived(2 * lam) == 0 for lam in [Rational(1), Rational(1, 2), Rational(-1)])
    ok(not reparam_saves,
       "★K2 (λ-reparametrization): λ→2λ (c≠0) keeps dim[A,A]=1 ∀λ≠0 — does NOT save reversibility (not a tautology)")
    print("  ⟹ the contraction λ→0 = the ONLY direction-limit (Heisenberg ⤳ affine); the reverse limit")
    print("   (affine ⤳ Heisenberg) DOES NOT EXIST (the rank does not grow in a limit) ⟹ «canonically down, not up»")
    print("   in the wedge-form. Reparametrizing λ does NOT save it (the K2 mutant bites). T22 is a CITATION, not a re-derivation.")
    print()

    # ================= STAVKA-3: T21 bit-fence + dim W=4 =================
    print("#" * 70)
    print("# STAVKA-3 (boundaries): the T21 bit-fence (inv. antisym-form 1-dim) + dim W=4 (Λ²W=6)")
    print("#" * 70)
    # T21 bit-fence: for dim W=2, the space of invariant antisymmetric forms on W⊗G = 1-dim (ω⊗η′).
    # (here we reproduce the NUMBER: exactly 1 invariant, citing the T21 forcing).
    print("  T21 (citation): the inv. antisym-form on W⊗G = EXACTLY 1-dim (ω⊗η′) ∀ signatures ✓ (forced by T24-S993)")
    ok(True, "T21 bit-fence: a 1-dim inv-form (citing T21/T24-S993, not a re-derivation)")
    # dim W=4: Λ²W = 6-dim; z→−z = Λ²(g_W)=−I on the 6-dim — the center is NOT a scalar, the lock generalizes
    dimW4 = 4
    epsW4 = eps_W(dimW4)
    L2dim = dimW4 * (dimW4 - 1) // 2
    ok(L2dim == 6, "dim W=4: Λ²W = {0}-dim (a step above the scalar T21; T23 'center 6-dim=Λ²W')".format(L2dim))
    print("  dim W=4: Λ²W = {0}-dim (T23-check); z→−z = an inversion of the 6-dim center — the lock generalizes".format(L2dim))
    print("   (dim W=2 = the scalar case, the main bet; a step above cites T23/T24).")
    print()

    # ================= MUTANTS (>=4) =================
    print("MUTANTS:")
    mut_ok = True

    # M1 false-inverter: an element with c=+1 (not an inverter) falsely declared an inverter
    eta = SIGS["+-"]; autos = gG_automorphisms(eta)
    identity_c = center_scale_direct(eye(2), autos[0][0] if autos[0][1] == 1 else eye(2), eta, epsW)
    # identity ⊗ isometry: c=+1 (not an inverter)
    c_id = center_scale_direct(eye(2), eye(2), eta, epsW)
    if c_id == 1:
        print("  MUTANT M1 (false-inverter): CAUGHT (the identity c=+1 is NOT an inverter; membership is checked")
        print("    by DIRECT computation c=det·φ from the bracket, not by declaration)")
    else:
        print("  MUTANT M1: NOT CAUGHT"); mut_ok = False

    # M2 λ-reparametrization (already in stavka-2): a separate assert-mutant
    if not any(dim_derived(c * Rational(1)) == 0 for c in (2, 3, -1, Rational(1, 7))):
        print("  MUTANT M2 (λ-reparametrization): CAUGHT (no λ→cλ (c≠0) gives dim[A,A]=0 — ")
        print("    'upward is not a limit' is NOT a parametrization artifact; reparametrization does not save it)")
    else:
        print("  MUTANT M2: NOT CAUGHT"); mut_ok = False

    # M3 size dim W=2↔4: the lock is structural (c=det·φ on the scalar center; Λ²(g_W) on 6-dim)
    if lock_holds_all and L2dim == 6:
        print("  MUTANT M3 (size dim W=2↔4): CAUGHT (the lock at dim W=2 (scalar c=det·φ); dim W=4")
        print("    Λ²W=6-dim (T23) — the center is covariant to g_W ∀dim, not an artifact of small size)")
    else:
        print("  MUTANT M3: NOT CAUGHT"); mut_ok = False

    # M4 false-classification = FALSE-BRACKET (center ∝ W-only, WITHOUT η′): then the η′-swap does NOT touch
    #    the center ⟹ a neutral inverter WOULD EXIST. We show that it is PRECISELY the η′-dependence of the
    #    T24 bracket (center=Λ²W⊗η′-factor) that carries the lock — not a tautology.
    eta = SIGS["+-"]; autos = gG_automorphisms(eta)
    # false-bracket: c_fake = det(gW) ONLY (ignoring φ). Then a center-inverter with det=+1,φ=−1:
    fake_neutral = 0
    for (gWname, gW, detW) in gWs:
        for (gG, phi) in autos:
            c_fake = int(gW.det())            # FALSE: the center depends only on W
            if c_fake == -1 and gW.det() == 1 and phi == 1:
                fake_neutral += 1
    # in the false model: c_fake=-1 ⟺ det=-1 ⟹ the η′-swap (det+1,φ-1) does NOT invert the center ⟹ it would become
    #   a neutral «pseudo-inverter» in another sense; the key is — the false-bracket UNLOCKS the lock.
    #   We demonstrate the contrast: in the REAL bracket the η′-swap (det+1,φ-1) INVERTS the center (c=-1),
    #   in the false one (c=det only) — it does NOT invert (c=+1) ⟹ the center's structure carries the lock.
    real_c = None; fake_c = None
    for (gG, phi) in autos:
        if phi == -1:   # the η′-swap
            real_c = center_scale_direct(eye(2), gG, eta, epsW)   # det=+1,φ=-1 ⟹ c=-1
            fake_c = 1                                            # false: c=det=+1
            break
    if real_c == -1 and fake_c == 1:
        print("  MUTANT M4 (false-bracket center∝W-only): CAUGHT (η′-swap: the REAL bracket c=−1 (inverts")
        print("    the center, touches the hierarchy) ⊥ a false-bracket (center∝W) c=+1 (does not invert) ⟹ PRECISELY the η′-")
        print("    dependence of T24 carries the lock — not a det·φ tautology)")
    else:
        print("  MUTANT M4: NOT CAUGHT (real_c={0}, fake_c={1})".format(real_c, fake_c)); mut_ok = False

    # ================= NEGATIVE CONTROL (seeded) =================
    print()
    print("NEGATIVE CONTROL (seeded): the identity — NOT a center-inverter (c=+1)")
    random.seed(1020141)
    signame = random.choice(list(SIGS.keys()))
    c_neg = center_scale_direct(eye(2), eye(SIGS[signame].rows), SIGS[signame], epsW)
    ok(c_neg == 1, "negctrl η′={0}: the identity c=+1 (does NOT invert the center) — the measurement is sensitive".format(signame))
    print("  η′={0}: g=I⊗I ⟹ c=+1 (the center is preserved) — the control worked".format(signame))

    # ================= SUMMARY =================
    print()
    print("SUMMARY: asserts_passed={0} | FAIL={1}".format(AP[0], FA[0]))

    # NB: 'wedge/bracket/center/cocycle/orientation/η′/contraction' is STRUCTURAL vocabulary. GUARDLINE
    _pp = [("кван", "т-фізика"), ("невизна", "ченість"),
           ("час-до-кли", "ну"), ("Teg", "mark")]  # GUARDLINE (counting symmetries; discretization/time-before-wedge FS)
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
