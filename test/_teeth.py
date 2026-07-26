# -*- coding: utf-8 -*-
"""`_teeth.py` — SHARED ASSERTION HARNESS WITH TEETH (author's word "stands", 2026-07-22).

WHAT THIS FIXES
-----------
The `S1055::T9` class — an **assert that cannot fail**. The statement is TRUE, the test is
GREEN, and the evidence is empty, because the detector says YES the same way in a world
where the object is correct and in a world where it's broken. Such a test is not a witness:
it does not distinguish worlds.

Instances the harness is built on (not hypothetical — both caught live in the repo):
  · `S1055::T9` — the assert "factor = 1/column-curvature": `(T/C)/T ≡ 1/C ∀C`, an identity;
  · `S1059::T3` — "no term of the series is empty": true for ANY non-degenerate
    operator, including the one where the unit really is 2 hops ⟹ it distinguishes nothing;
  · `S1059::T4` — "the per-bond count does not move under P": `P` never reached the
    computation, four iterations counted the same quantity.

WHY THIS EXACT SIGNATURE (and not "a rule in the README")
---------------------------------------------------
Carved lesson S1056: **a lesson becomes a gate only by becoming a REQUIRED FIELD.** So here
the detector is a **function of the WORLD**, not a ready-made boolean: a ready-made boolean
cannot be re-asked on another object, and that is exactly why the local `def ok(cond, msg)`
in test_1..test_3 are structurally unable to catch this class. Forgetting the negative
control is impossible not from discipline, but because `ok()` cannot be called without one.

BOUNDARIES — READ BEFORE RELYING ON THIS (otherwise the line lies)
------------------------------------------------------------
1. The harness measures **SENSITIVITY, not correctness.** A test can have teeth and still
   be about the wrong thing.
2. It **does not check whether the negative world is honest.** One can name a world where
   the detector fails trivially and buy a fake tooth. Cutting corners becomes more
   expensive, not impossible.
3. It **does not catch "stated instead of derived"** — an absent derivation is not an
   empty test. That is a different layer (proof / Lean), and this file neither replaces
   it nor competes with it.
4. `ok_contrast` requires a PAIR, because true invariance and blindness of the machinery
   give the same output. Without a quantity that MUST move, "didn't move" is worth nothing.

RELATION TO OMEGA'S DIRECTION ("a source of values + an assert against it", 2026-07-22)
--------------------------------------------------------------------------------
Different defects, not competitors: her mechanism catches **DESYNC** (a number diverged
from its source), this one catches **EMPTINESS** (a test does not distinguish worlds).
On 2026-07-22 the repo carried both.

SELF-TEST
-------------
`python _teeth.py` runs the harness on FOUR real S1059 detectors and requires
two to be flagged empty and two to pass. A harness that checks teeth must
have its own: if it stops catching, its own self-test fails with exit 1.
"""
import sys

_pass, _fail, _void = [], [], []


def _say(s=''):
    print(s)


def ok(detector, world, msg, *, must_fail_on):
    """An assert that is required to be able to fail.

    detector    — a FUNCTION OF THE WORLD returning bool. Not a ready-made value: a
                  ready-made value cannot be re-asked on another object.
    world       — our object (the one the statement is about).
    must_fail_on— (name, world) or a list of such pairs: object(s) where the detector
                  MUST say NO. If it says YES there too, the assert self-destructs (☠),
                  because it distinguishes nothing.

    ★The negative control here runs the SAME `detector` that carries the evidence — this
    is not a convention, it's the only thing the signature allows. Alpha's class (S1059
    rev.2): it is not enough to HAVE a negative control, it must hit the same detector;
    otherwise we get two greens that know nothing about each other (exactly S1059::T3:
    the test measured non-zero-ness of a series term, the negative control measured
    reachability in r-space).

    Returns the detector's verdict on `world`; ☠ is counted separately from ✗ and also
    gives exit≠0.
    """
    negatives = must_fail_on if isinstance(must_fail_on, list) else [must_fail_on]
    if not negatives:
        raise ValueError("must_fail_on is empty: an assert without a negative world is forbidden")
    blind = [name for name, w in negatives if detector(w)]
    if blind:
        _void.append(msg)
        _say(f"  ☠ EMPTY  {msg}")
        for name in blind:
            _say(f"       the detector said YES on «{name}» too ⟹ it does NOT distinguish worlds")
        return False
    verdict = bool(detector(world))
    (_pass if verdict else _fail).append(msg)
    _say(f"  {'✓' if verdict else '✗ FAIL'} {msg}")
    _say(f"       (tooth: said NO on «{negatives[0][0]}» ⟹ able to fail)")
    return verdict


def ok_contrast(probe, params, msg):
    """An assert about INVARIANCE under a parameter — with a mandatory contrast.

    probe(p) → (what_stands, what_moves) — ONE call, ONE body, ONE object.
    params   — parameter values (≥2).

    Why a pair: true invariance and "the parameter never reached the computation" give
    the same output. A moving quantity proves the machinery is CAPABLE of seeing motion —
    and only then is "didn't move" a result, not blindness.

    ★WHY ONE FUNCTION, NOT TWO (Alpha's class, S1059 rev.2, 2026-07-22): "it is not
    enough to HAVE a negative control — it must hit the SAME detector that carries the
    evidence; otherwise we get two greens that know nothing about each other." That is
    exactly what happened in S1059::T3, where the main test measured non-zero-ness of a
    series term, and the negative control measured reachability in r-space. Two separate
    functions would let the contrast be computed by DIFFERENT machinery; returning the
    pair from one body makes that impossible by form.
    """
    if len(params) < 2:
        raise ValueError("contrast requires ≥2 parameter values")
    got = [probe(p) for p in params]
    if any((not isinstance(g, tuple)) or len(g) != 2 for g in got):
        raise ValueError("probe(p) must return a PAIR (what_stands, what_moves) "
                         "from ONE computation — separate functions are not accepted by form")
    inv = [g[0] for g in got]
    mov = [g[1] for g in got]
    inv_holds = len(set(map(str, inv))) == 1
    mov_moves = len(set(map(str, mov))) > 1
    if not mov_moves:
        _void.append(msg)
        _say(f"  ☠ EMPTY  {msg}")
        _say(f"       the control quantity ALSO didn't move ({mov}) ⟹ the parameter never reached")
        _say(f"       the computation; the test reformulated the definition instead of measuring a response")
        return False
    (_pass if inv_holds else _fail).append(msg)
    _say(f"  {'✓' if inv_holds else '✗ FAIL'} {msg}")
    _say(f"       invariant {inv} ⊥ control {mov} — the same machinery MOVES the control,")
    _say(f"       so 'the invariant stands' is a measurement, not blindness")
    return inv_holds


def report(name="harness"):
    """Summary + exit code. ☠ gives exit≠0 on par with ✗ — an empty assert IS a failure."""
    _say()
    _say(f"SCORE [{name}]: {len(_pass)} ✓ / {len(_fail)} ✗ / {len(_void)} ☠ EMPTY")
    if _void:
        _say("★EMPTY ASSERTS (not witnesses — remove or give them a negative world):")
        for m in _void:
            _say(f"   ☠ {m}")
    return 1 if (_fail or _void) else 0


def reset():
    """For the self-test: start the count from zero."""
    _pass.clear(); _fail.clear(); _void.clear()


# ───────────────────────────── SELF-TEST ─────────────────────────────
def _selftest():
    """A run on FOUR real S1059 detectors. The expectations are carved in here:
    T3 and T4-as-it-stood MUST be flagged empty, T2 and T4-meaningful must pass."""
    import sympy as sp

    _say("=" * 78)
    _say("SELF-TEST `_teeth.py` — on real detectors from probe S1059")
    _say("=" * 78)

    g1, g2 = sp.symbols('g1 g2', nonzero=True)
    G0 = sp.Matrix([[g1, 0], [0, g2]])
    f, fb, q = sp.symbols('f fbar q', nonzero=True)

    def no_gaps(H):
        term = G0
        for _ in range(7):
            if sp.simplify(term.norm()) == 0:
                return False
            term = sp.simplify(term * H * G0)
        return True

    _say("\n[1] S1059::T3 «no omissions in the series» — GREEN as it stood in the probe")
    reset()
    ok(no_gaps, sp.Matrix([[0, f], [fb, 0]]),
       "T3: no term of the series n=0..6 is empty ⟹ the unit is not p≥2",
       must_fail_on=("a reach-2 operator — there the unit REALLY IS 2 edges",
                     sp.Matrix([[q, 0], [0, q]])))
    t3_void = len(_void) == 1

    _say("\n[2] S1059::T2 (algebraic route) — has teeth")
    reset()

    def no_local_root(pair):
        fp, fbp = pair
        num, _ = sp.fraction(sp.cancel(sp.together(-fp * fbp)))
        _, facs = sp.factor_list(sp.expand(num))
        return any(m % 2 == 1 and not b.is_Symbol for b, m in facs)

    X1, X2 = sp.symbols('X1 X2', nonzero=True)
    f2, fb2 = 1 + X1 + X2, 1 + 1 / X1 + 1 / X2
    ok(no_local_root, (f2, fb2),
       "T2: no local root K²=H exists (d=2)",
       must_fail_on=[("H² — there a local root DOES exist (H itself)", (f2 * fb2, f2 * fb2)),
                     ("f=1 — no node, a root is possible", (sp.Integer(1), sp.Integer(1)))])
    t2_pass = len(_pass) == 1 and not _void

    def hop_matrix(N, reach=1):
        A = sp.zeros(N, N)
        for i in range(N):
            A[i, (i + reach) % N] += 1
            A[(i + reach) % N, i] += 1
        return A

    def min_power(A, u, v, nmax=14):
        P = sp.eye(A.rows)
        for n in range(1, nmax + 1):
            P = P * A
            if sp.simplify(P[u, v]) != 0:
                return n
        return None

    _say("\n[3] S1059::T4 as it stood in the probe — P never reached the computation")
    reset()
    A_fixed = hop_matrix(12, 1)

    def probe_as_was(P):
        """A literal reconstruction: the matrix is built OUTSIDE the parameter, so both
        quantities read an object that P never touched."""
        return min_power(A_fixed, 0, 1), min_power(A_fixed, 0, 1)

    ok_contrast(probe_as_was, [2, 3, 4, 6], "T4-as-it-stood: per-bond count = 1 ∀P")
    t4_void = len(_void) == 1

    _say("\n[4] S1059::T4 in a meaningful form — P enters the CONSTRUCTION of the description")
    reset()

    def supercell(N, P):
        perm = [(i % P) * (N // P) + (i // P) for i in range(N)]
        base = hop_matrix(N, 1)
        A = sp.zeros(N, N)
        for i in range(N):
            for j in range(N):
                A[perm[i], perm[j]] = base[i, j]
        return A, perm

    def probe_meaningful(P):
        """ONE construction, ONE object: the same supercell matrix yields both
        quantities — otherwise 'invariant' and 'control' could read different worlds."""
        A, perm = supercell(12, P)
        return min_power(A, perm[0], perm[1]), min_power(A, perm[0], perm[P])

    ok_contrast(probe_meaningful, [2, 3, 4, 6],
                "T4-meaningful: the per-BOND count is invariant under relabelling P")
    t4_pass = len(_pass) == 1 and not _void

    _say()
    _say("=" * 78)
    _say("SELF-TEST VERDICT (expectations carved into the code, not read off the output):")
    checks = [("T3 flagged empty", t3_void),
              ("T2 passed with teeth", t2_pass),
              ("T4-as-it-stood flagged empty", t4_void),
              ("T4-meaningful passed with contrast", t4_pass)]
    for name, good in checks:
        _say(f"  {'✓' if good else '✗ FAIL'} {name}")
    bad = [n for n, g in checks if not g]
    if bad:
        _say(f"★HARNESS IS BROKEN: {bad} — it no longer catches what it exists to catch.")
    else:
        _say("★the harness catches exactly the two that passed green in the probe, and lets")
        _say("  through the two that have teeth. With no knowledge of physics whatsoever.")
    _say("=" * 78)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding='utf-8')
    sys.exit(_selftest())
