# -*- coding: utf-8 -*-
# DIM: na (dimensionless/exact symbolic count; no spatial lattice)
"""
S920 (W29): sign convention / orientability of the time-sector in signature (p,q).

η = diag(+1×p, −1×q); ⟨u,v⟩ = uᵀ η v; time-sector = minus-axes (indices p..n−1);
v is timelike ⟺ ⟨v,v⟩ < 0. Coordinate flows exp(sM), M∈so(p,q)
(defining equation Mη + ηMᵀ = 0): a same-sign axis pair → rotation,
a mixed-sign pair → boost.

QUESTION: does the set of timelike directions split into two components
that the flows exp(sM) do not connect — depending on q.

Exact symbolics (sympy): Integer/Rational/symbols/cos/sin/cosh/sinh/π/log.
NO float/random/time/datetime. Logical (non-machine) links are marked [LOG].
exit 0 only if: all checks OK, mutants 3/3 CAUGHT, no section empty.
"""

import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

from sympy import (Symbol, symbols, Integer, Matrix, eye, zeros, diag,
                   cos, sin, cosh, sinh, exp, pi, log, simplify, expand)

# ---------------------------------------------------------------- bookkeeping
CHECKS_TOTAL = 0
FAILS = []                 # (section, label)
SECTION_COUNTS = {}        # section -> check count
SECTION_FAILS = {}         # section -> failure count
VERDICTS = []              # (section, text, ok)
MUTANTS_CAUGHT = 0
CUR = None


def section(name):
    global CUR
    CUR = name
    SECTION_COUNTS[name] = 0
    SECTION_FAILS[name] = 0
    print("\n" + "=" * 78)
    print("SECTION " + name)
    print("=" * 78)


def check(label, ok):
    global CHECKS_TOTAL
    CHECKS_TOTAL += 1
    SECTION_COUNTS[CUR] += 1
    if not ok:
        SECTION_FAILS[CUR] += 1
        FAILS.append((CUR, label))
    print("  [%s] %s" % ("OK " if ok else "FAIL", label))
    return ok


def verdict(text):
    ok = (SECTION_FAILS[CUR] == 0) and (SECTION_COUNTS[CUR] > 0)
    VERDICTS.append((CUR, text, ok))
    print("  SECTION COUNTER: checks=%d, failures=%d"
          % (SECTION_COUNTS[CUR], SECTION_FAILS[CUR]))
    print("  VERDICT: %s — %s" % (text, "OK" if ok else "FAILURE"))


# ---------------------------------------------------------------- geometry
def eta_mat(p, q):
    return diag(*([Integer(1)] * p + [Integer(-1)] * q))


def axis_sign(p, i):
    return Integer(1) if i < p else Integer(-1)


def generator(p, q, i, j):
    """Coordinate generator of so(p,q) on the axis pair (i,j), i<j.
    A same-sign pair -> rotation, a mixed-sign pair -> boost."""
    n = p + q
    M = zeros(n)
    if axis_sign(p, i) * axis_sign(p, j) > 0:
        M[i, j] = Integer(-1)
        M[j, i] = Integer(1)
        kind = "rot"
    else:
        M[i, j] = Integer(1)
        M[j, i] = Integer(1)
        kind = "boost"
    return M, kind


def flow(p, q, i, j, s):
    """Explicit E(s): a 2×2 block cos/sin (rot) or cosh/sinh (boost) on (i,j).
    Machine-checked below: E(0)=1, dE/ds = M·E, EᵀηE = η."""
    n = p + q
    M, kind = generator(p, q, i, j)
    E = eye(n)
    if kind == "rot":
        E[i, i] = cos(s); E[j, j] = cos(s)
        E[i, j] = -sin(s); E[j, i] = sin(s)
    else:
        E[i, i] = cosh(s); E[j, j] = cosh(s)
        E[i, j] = sinh(s); E[j, i] = sinh(s)
    return E, M, kind


def is_zero_mat(A):
    return all(simplify(x) == 0 for x in A)


def basis_vec(n, k):
    e = zeros(n, 1)
    e[k] = Integer(1)
    return e


# ================================================================ T1
section("T1 — SPLIT AT q=1")
for (p, q) in [(2, 1), (3, 1), (4, 1)]:
    n = p + q
    t = n - 1                      # the single time axis (last)
    eta = eta_mat(p, q)
    vs = symbols("v0:%d" % n, real=True)
    v = Matrix(n, 1, vs)
    v_t0 = v.subs(vs[t], Integer(0))
    Q = expand((v_t0.T * eta * v_t0)[0, 0])
    sq = expand(sum(vs[i] ** 2 for i in range(p)))
    print("  (%d,%d): ⟨v,v⟩|_(v_t=0) = %s" % (p, q, Q))
    check("(%d,%d): ⟨v,v⟩|_(v_t=0) − Σ_{i<p} v_i² == 0 (a sum of squares of real "
          "symbols ⟹ PSD exactly)" % (p, q), expand(Q - sq) == 0)
    et = basis_vec(n, t)
    check("(%d,%d): ⟨e_t,e_t⟩ = −1 < 0 (e_t is timelike)" % (p, q),
          (et.T * eta * et)[0, 0] == Integer(-1))
    check("(%d,%d): ⟨−e_t,−e_t⟩ = −1 < 0 (−e_t is timelike)" % (p, q),
          ((-et).T * eta * (-et))[0, 0] == Integer(-1))
print("  CHAIN:")
print("   [MACH] v_t=0 ⟹ ⟨v,v⟩ = Σ v_i² ≥ 0 (identity above)")
print("   ⟹ timelike (⟨v,v⟩<0) ⟹ v_t ≠ 0 (contrapositive of the machine identity)")
print("   [LOG] v ↦ v_t is continuous and does not vanish on the timelike set ⟹ sign(v_t) is defined")
print("        and locally constant ⟹ {v_t>0} and {v_t<0} — open, disjoint, cover")
print("        all timelike vectors ⟹ ≥2 connected components")
print("   [MACH] both nonempty: e_t (v_t=+1) and −e_t (v_t=−1) are timelike")
verdict("at q=1 the set of timelike vectors SPLITS into ≥2 components (v_t>0 ⊥ v_t<0)")

# ================================================================ T2
section("T2 — PRESERVATION AT q=1")
s = Symbol("s", real=True)
for (p, q) in [(2, 1), (3, 1), (4, 1)]:
    n = p + q
    t = n - 1
    eta = eta_mat(p, q)
    vs = symbols("v0:%d" % n, real=True)
    v = Matrix(n, 1, vs)
    n_rot, n_boost = 0, 0
    for i in range(n):
        for j in range(i + 1, n):
            E, M, kind = flow(p, q, i, j, s)
            check("(%d,%d) pair(%d,%d) %s: Mη + ηMᵀ = 0 (M ∈ so(p,q))"
                  % (p, q, i, j, kind), is_zero_mat(M * eta + eta * M.T))
            check("(%d,%d) pair(%d,%d) %s: E(0)=𝟙 ∧ dE/ds = M·E "
                  "(⟹ E=exp(sM) by uniqueness of the ODE solution [LOG])"
                  % (p, q, i, j, kind),
                  is_zero_mat(E.subs(s, Integer(0)) - eye(n))
                  and is_zero_mat(E.diff(s) - M * E))
            check("(%d,%d) pair(%d,%d) %s: EᵀηE = η symbolically in s (isometry)"
                  % (p, q, i, j, kind), is_zero_mat(E.T * eta * E - eta))
            if kind == "rot":
                n_rot += 1
                w = E * v
                check("(%d,%d) rot(%d,%d): (E(s)·v)_t = v_t identically (time is untouched)"
                      % (p, q, i, j), simplify(w[t] - v[t]) == 0)
            else:
                n_boost += 1
                et = basis_vec(n, t)
                w = E * et
                check("(%d,%d) boost(%d,%d): (E(s)·e_t)_t = cosh(s)"
                      % (p, q, i, j), simplify(w[t] - cosh(s)) == 0)
                check("(%d,%d) boost(%d,%d): cosh(s) − 1 − 2·sinh²(s/2) = 0 "
                      "⟹ cosh(s) = 1 + 2·sinh²(s/2) ≥ 1 > 0 (sinh²≥0 for real s)"
                      % (p, q, i, j),
                      simplify(cosh(s) - 1 - 2 * sinh(s / 2) ** 2) == 0)
    print("  (%d,%d): total pairs=%d, rotations=%d, boosts=%d"
          % (p, q, n_rot + n_boost, n_rot, n_boost))
print("  CHAIN:")
print("   [MACH] EᵀηE = η ⟹ the flow preserves ⟨v,v⟩ ⟹ timelike stays timelike")
print("   [MACH, T1] a timelike vector with v_t=0 is impossible (⟨v,v⟩≥0 there)")
print("   [LOG] s ↦ (E(s)v)_t is continuous; a sign change would require a zero (the")
print("        intermediate value theorem) — and a zero is forbidden on a timelike orbit [T1]")
print("   ⟹ every coordinate flow preserves the component (the sign of v_t)")
verdict("at q=1 the two cones (v_t>0 / v_t<0) are INVARIANT under all coordinate flows")

# ================================================================ T3
section("T3 — MERGING AT q≥2")
theta = Symbol("theta", real=True)
for (p, q) in [(2, 2), (3, 2), (2, 3), (3, 3)]:
    n = p + q
    eta = eta_mat(p, q)
    a, b = p, p + 1                # the two time (minus) axes
    E, J, kind = flow(p, q, a, b, theta)
    check("(%d,%d): pair(%d,%d) ∈ {−−} ⟹ the generator J is a ROTATION" % (p, q, a, b),
          kind == "rot")
    check("(%d,%d): Jη + ηJᵀ = 0" % (p, q), is_zero_mat(J * eta + eta * J.T))
    ea = basis_vec(n, a)
    g = E * ea                     # γ(θ) = exp(θJ)·e_a
    Qg = (g.T * eta * g)[0, 0]
    check("(%d,%d): ⟨γ(θ),γ(θ)⟩ = −1 IDENTICALLY in θ (the whole path is timelike)" % (p, q),
          simplify(Qg + 1) == 0)
    check("(%d,%d): γ(0) = e_a exactly" % (p, q),
          is_zero_mat(g.subs(theta, Integer(0)) - ea))
    check("(%d,%d): γ(π) = −e_a exactly (cos π=−1, sin π=0)" % (p, q),
          is_zero_mat(g.subs(theta, pi) + ea))
print("  CHAIN:")
print("   [MACH] γ(θ)=exp(θJ)e_a: the norm ≡ −1 along the whole path, γ(0)=e_a, γ(π)=−e_a")
print("   [LOG] θ ↦ γ(θ) is continuous ⟹ e_a and −e_a lie in ONE component of the timelike set")
print("   ⟹ no continuous \"orientation of the time-sector\" is invariant: the flow in the")
print("     {−−}-plane carries e_a to −e_a without leaving the timelike set")
verdict("at q≥2 e_a and −e_a ARE CONNECTED by a path within the timelike set — there is "
        "no invariant orientation of the time-sector")

# ================================================================ T4
section("T4 — DUALITY (plus-sector)")
# (i) p=1: split of the plus-sector (⟨v,v⟩>0)
for (p, q) in [(1, 2), (1, 3), (1, 4)]:
    n = p + q
    sx = 0                          # the single plus-axis
    eta = eta_mat(p, q)
    vs = symbols("v0:%d" % n, real=True)
    v = Matrix(n, 1, vs)
    v_s0 = v.subs(vs[sx], Integer(0))
    Q = expand((v_s0.T * eta * v_s0)[0, 0])
    msq = expand(-sum(vs[i] ** 2 for i in range(1, n)))
    print("  (%d,%d): ⟨v,v⟩|_(v_s=0) = %s" % (p, q, Q))
    check("(%d,%d): ⟨v,v⟩|_(v_s=0) − (−Σ_{i≥1} v_i²) == 0 (a minus-sum of squares "
          "⟹ ≤ 0 exactly)" % (p, q), expand(Q - msq) == 0)
    e0 = basis_vec(n, 0)
    check("(%d,%d): ⟨e_0,e_0⟩ = +1 > 0 and ⟨−e_0,−e_0⟩ = +1 > 0 (both components "
          "nonempty)" % (p, q),
          (e0.T * eta * e0)[0, 0] == Integer(1)
          and ((-e0).T * eta * (-e0))[0, 0] == Integer(1))
print("   ⟹ [MACH] ⟨v,v⟩>0 requires v_s≠0; [LOG] sign(v_s) continuous ⟹ two components")
# (ii) p>=2: merging of the plus-sector
for (p, q) in [(2, 2), (3, 2)]:
    n = p + q
    eta = eta_mat(p, q)
    E, J, kind = flow(p, q, 0, 1, theta)   # the {++}-plane (0,1)
    check("(%d,%d): pair(0,1) ∈ {++} ⟹ J is a rotation; Jη+ηJᵀ=0" % (p, q),
          kind == "rot" and is_zero_mat(J * eta + eta * J.T))
    e0 = basis_vec(n, 0)
    g = E * e0
    Qg = (g.T * eta * g)[0, 0]
    check("(%d,%d): ⟨γ(θ),γ(θ)⟩ = +1 IDENTICALLY (the whole path is among plus-vectors)" % (p, q),
          simplify(Qg - 1) == 0)
    check("(%d,%d): γ(0)=e_0, γ(π)=−e_0 exactly" % (p, q),
          is_zero_mat(g.subs(theta, Integer(0)) - e0)
          and is_zero_mat(g.subs(theta, pi) + e0))
verdict("MIRROR: the plus-sector splits at p=1 (two components), at p≥2 it is "
        "connected by a flow in the {++}-plane. A sector whose sign has multiplicity-1 "
        "is orientable, with ≥2 it is not")

# ================================================================ M
section("M — MUTANTS (each MUST be caught)")

# --- m1: "at q=1 the timelike set is connected: e_t and −e_t are joined by a path within the timelike set"
print("  m1: FALSE claim — \"at q=1 e_t and −e_t are joined by a path within the timelike set\".")
print("      Refutation: a path v(τ), v_t(0)=+1, v_t(1)=−1, continuous ⟹ [LOG: the intermediate")
print("      value theorem] ∃τ*: v_t(τ*)=0; at that point [MACH, the T1 identity]:")
m1_ok = True
for (p, q) in [(2, 1), (3, 1), (4, 1)]:
    n = p + q
    t = n - 1
    eta = eta_mat(p, q)
    vs = symbols("v0:%d" % n, real=True)
    v = Matrix(n, 1, vs)
    Q = expand((v.subs(vs[t], Integer(0)).T * eta * v.subs(vs[t], Integer(0)))[0, 0])
    ok = expand(Q - sum(vs[i] ** 2 for i in range(p))) == 0
    m1_ok = m1_ok and ok
    check("m1 (%d,%d): the barrier v_t=0 ⟹ ⟨v,v⟩ = Σ v_i² ≥ 0 — timelike-ness (<0) "
          "is IMPOSSIBLE at the crossing point" % (p, q), ok)
if m1_ok:
    MUTANTS_CAUGHT += 1
    print("  m1 CAUGHT: any path e_t → −e_t would have to pass through v_t=0, where ⟨v,v⟩≥0 — "
          "a contradiction.")

# --- m2: "a boost flips the cone: ∃s, (exp(sK)e_t)_t < 0"
print("  m2: FALSE claim — \"there exists s with (exp(sK)·e_t)_t < 0\".")
p, q = 3, 1
n = p + q
t = n - 1
E, K, kind = flow(p, q, 0, t, s)
et = basis_vec(n, t)
w = E * et
c1 = check("m2 (3,1) boost(0,3): (exp(sK)·e_t)_t = cosh(s) symbolically",
           kind == "boost" and simplify(w[t] - cosh(s)) == 0)
c2 = check("m2: cosh(s) = 1 + 2·sinh²(s/2) identically ⟹ cosh(s) ≥ 1 > 0 ∀ real s",
           simplify(cosh(s) - 1 - 2 * sinh(s / 2) ** 2) == 0)
c3 = check("m2 (a control exact point, no float): s=log(2) ⟹ cosh(log 2) = 5/4 > 0",
           simplify(w[t].subs(s, log(2)).rewrite(exp)) == Integer(5) / Integer(4))
if c1 and c2 and c3:
    MUTANTS_CAUGHT += 1
    print("  m2 CAUGHT: the t-component of the boost orbit = cosh(s) ≥ 1, the sign cannot be flipped.")

# --- m3: broken isometry check: EᵀE = 𝟙 instead of EᵀηE = η
print("  m3: BROKEN isometry check (Eᵀ·E = 𝟙 instead of Eᵀ·η·E = η) on a boost.")
eta = eta_mat(3, 1)
honest = is_zero_mat(E.T * eta * E - eta)          # honest: must pass
D = E.T * E - eye(n)                                # broken: residual
d00 = simplify(D[0, 0])
c1 = check("m3: the HONEST check EᵀηE=η PASSES on the boost (a boost is an isometry)", honest)
c2 = check("m3: the broken residual (EᵀE−𝟙)[0,0] = cosh²+sinh²−1 = 2·sinh²(s) ≠ 0 "
           "(as a symbolic expression)", simplify(d00 - 2 * sinh(s) ** 2) == 0 and d00 != 0)
c3 = check("m3: exact point s=log(2): (EᵀE−𝟙)[0,0] = 9/8 ≠ 0 — the broken check FALSELY "
           "fails a genuine isometry",
           simplify(d00.subs(s, log(2)).rewrite(exp)) == Integer(9) / Integer(8))
if c1 and c2 and c3:
    MUTANTS_CAUGHT += 1
    print("  m3 CAUGHT: the discrepancy of the broken check vs the honest one is shown exactly "
          "(honest OK, broken gives 9/8≠0 at s=log 2).")

check("MUTANTS: caught 3/3", MUTANTS_CAUGHT == 3)
verdict("all 3 mutants caught (m1 barrier-identity, m2 cosh≥1, m3 discrepancy of the "
        "broken check)")

# ================================================================ SUMMARY
section("SUMMARY")
for name, text, ok in VERDICTS:
    print("  [%s] %s: %s" % ("OK " if ok else "FAIL", name.split(" —")[0], text))
print()
print("  ROLLUP TABLE (signature | sector | minority-1? | orientable?):")
print("  " + "-" * 66)
rows = [
    ("(2,1)", "time (−, q=1)",  "YES", "YES"),
    ("(3,1)", "time (−, q=1)",  "YES", "YES"),
    ("(4,1)", "time (−, q=1)",  "YES", "YES"),
    ("(1,2)", "plus (+, p=1)", "YES", "YES"),
    ("(1,3)", "plus (+, p=1)", "YES", "YES"),
    ("(1,4)", "plus (+, p=1)", "YES", "YES"),
    ("(2,2)", "time (−, q=2)",  "NO",  "NO"),
    ("(2,2)", "plus (+, p=2)", "NO",  "NO"),
    ("(3,2)", "time (−, q=2)",  "NO",  "NO"),
    ("(3,2)", "plus (+, p=3)", "NO",  "NO"),
    ("(2,3)", "time (−, q=3)",  "NO",  "NO"),
    ("(3,3)", "time (−, q=3)",  "NO",  "NO"),
]
print("  %-8s| %-16s| %-12s| %s" % ("sig.", "sector", "minority-1?", "orientable?"))
print("  " + "-" * 66)
for r in rows:
    print("  %-8s| %-16s| %-12s| %s" % r)
print("  " + "-" * 66)
print("  LAW: a sector whose sign has multiplicity 1 is orientable (two invariant cones);")
print("       multiplicity ≥2 — a rotation in the same-sign plane connects ±e within")
print("       the sector, the orientation is lost.")

# ---------------------------------------------------------------- final
print("\n" + "=" * 78)
print("COUNTERS:")
# SUMMARY — an overview section with no checks of its own; the "zero checked" guard
# applies to the measurement sections T1/T2/T3/T4/M.
empty_sections = [k for k, c in SECTION_COUNTS.items()
                  if c == 0 and not k.startswith("SUMMARY")]
for k in SECTION_COUNTS:
    print("  %-40s checks=%3d  failures=%d"
          % (k.split(" —")[0], SECTION_COUNTS[k], SECTION_FAILS[k]))
print("  TOTAL checks: %d, failures: %d, mutants caught: %d/3"
      % (CHECKS_TOTAL, len(FAILS), MUTANTS_CAUGHT))

ok_all = (len(FAILS) == 0 and MUTANTS_CAUGHT == 3
          and CHECKS_TOTAL > 0 and not empty_sections)
if empty_sections:
    print("  EMPTY SECTIONS (zero checked = an invalid measurement): %s" % empty_sections)
if FAILS:
    print("  FAILURES:")
    for sec, lab in FAILS:
        print("    - [%s] %s" % (sec, lab))
print("SUMMARY: %s" % ("ALL OK — exit 0" if ok_all else "THERE ARE FAILURES — exit 1"))
sys.exit(0 if ok_all else 1)
