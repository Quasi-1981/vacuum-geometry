# -*- coding: utf-8 -*-
# DIM: na (dimensionless/exact symbolic count; no spatial lattice)
"""
S918 (w29): "residue buildup on mirror cracks" — residual structure of the
antisymmetric form Omega in so(p,q).

Setup:
  eta = diag(+1 x p, -1 x q);  so(p,q) = {X : X*eta + eta*X^T = 0};
  Omega^T = -Omega, block-structured (consecutive 2-blocks [[0,w],[-w,0]], an odd axis left unpaired);
  A = eta^{-1} * Omega  (eta^{-1} = eta, since diag +-1);
  c(A) = centralizer = {X in so(p,q) : X*A - A*X = 0}.

Question: when is c(A) abelian; the pattern
  "non-abelian <=> degenerate (collision of eigenvalues of A, or a rank drop of Omega)".

EXACT arithmetic: sympy Integer/Rational. NO float, NO random, NO time/datetime.
Run: python S918_w29_condensate_on_mirror_cracks.py > S918_w29_run.log 2>&1
"""
import sys
from sympy import Matrix, Integer, Rational, Symbol, Poly, factor, gcd, diff, eye, zeros

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

lam = Symbol("lam")
FAIL = []


def fail(msg):
    FAIL.append(msg)
    print("  !! FAIL: " + msg)


# ---------------------------------------------------------------- building blocks
def eta_matrix(p, q):
    n = p + q
    return Matrix(n, n, lambda i, j: Integer(0) if i != j else (Integer(1) if i < p else Integer(-1)))


def build_omega(n, profile):
    """Block Omega: 2-blocks on axes (2b, 2b+1) with value profile[b]."""
    assert len(profile) == n // 2, "the profile must have n//2 values"
    Om = zeros(n, n)
    for b, w in enumerate(profile):
        a = 2 * b
        Om[a, a + 1] = Integer(w)
        Om[a + 1, a] = -Integer(w)
    return Om


def pfaffian(M):
    """Pfaffian by recursion on the first row (exact arithmetic).
    Pf(A) = sum_{j=2..n} (-1)^j a_{1j} Pf(A_without rows/columns 1,j) (1-based).
    Odd n -> 0 by definition."""
    n = M.rows
    if n % 2 == 1:
        return Integer(0)
    if n == 0:
        return Integer(1)
    total = Integer(0)
    for j in range(1, n):
        a = M[0, j]
        if a == 0:
            continue
        idx = [k for k in range(n) if k != 0 and k != j]
        sub = M.extract(idx, idx)
        total += Integer(-1) ** (j + 1) * a * pfaffian(sub)
    return total


def vec(M):
    n = M.rows
    return Matrix(n * n, 1, lambda i, j: M[i // n, i % n])


_SO_CACHE = {}


def so_basis(p, q):
    """Basis of so(p,q) via the exact solution of the linear system X*eta + eta*X^T = 0."""
    if (p, q) in _SO_CACHE:
        return _SO_CACHE[(p, q)]
    n = p + q
    eta = eta_matrix(p, q)
    rows = []
    for i in range(n):
        for j in range(n):
            row = [Integer(0)] * (n * n)
            # (X*eta + eta*X^T)[i,j] = X[i,j]*eta[j,j] + eta[i,i]*X[j,i]
            row[i * n + j] += eta[j, j]
            row[j * n + i] += eta[i, i]
            rows.append(row)
    C = Matrix(rows)
    null = C.nullspace()
    basis = [Matrix(n, n, lambda a, b, v=v: v[a * n + b]) for v in null]
    # basis sanity
    for B in basis:
        if not (B * eta + eta * B.T).is_zero_matrix:
            fail("so(%d,%d): a basis element does not satisfy the defining equation" % (p, q))
    if len(basis) != n * (n - 1) // 2:
        fail("dim so(%d,%d) = %d != n(n-1)/2 = %d" % (p, q, len(basis), n * (n - 1) // 2))
    _SO_CACHE[(p, q)] = basis
    return basis


def centralizer(basis, A, drop_last_rows=0):
    """c(A) in span(basis): exact nullspace of the system vec(B_k*A - A*B_k)*c = 0.
    drop_last_rows > 0 -> DELIBERATELY broken solver (mutant m1)."""
    n = A.rows
    cols = [vec(B * A - A * B) for B in basis]
    M = Matrix.hstack(*cols)
    if drop_last_rows:
        M = M[: M.rows - drop_last_rows, :]
    if M.rows == 0:
        null = [eye(len(basis)).col(k) for k in range(len(basis))]
    else:
        null = M.nullspace()
    cb = []
    for c in null:
        Z = zeros(n, n)
        for k, B in enumerate(basis):
            if c[k] != 0:
                Z = Z + c[k] * B
        cb.append(Z)
    return cb


def rank_of(mats):
    cols = [vec(Z) for Z in mats]
    if not cols:
        return 0
    return Matrix.hstack(*cols).rank()


def in_span(mats, W):
    if W.is_zero_matrix:
        return True
    if not mats:
        return False
    P = Matrix.hstack(*[vec(Z) for Z in mats])
    return Matrix.hstack(P, vec(W)).rank() == P.rank()


def closed_under_bracket(cb):
    for a in range(len(cb)):
        for b in range(a + 1, len(cb)):
            W = cb[a] * cb[b] - cb[b] * cb[a]
            if not in_span(cb, W):
                return False
    return True


def bracket_stats(cb):
    """(abelian?, dim [c,c])"""
    brs = []
    for a in range(len(cb)):
        for b in range(a + 1, len(cb)):
            W = cb[a] * cb[b] - cb[b] * cb[a]
            if not W.is_zero_matrix:
                brs.append(W)
    return (len(brs) == 0), rank_of(brs)


def eig_info(A):
    """Exact charpoly. Returns (m0 = multiplicity of zero, coll = whether there is a REPEATED
    NONzero eigenvalue (non-squarefree nonzero part), factorization)."""
    cp = Poly(A.charpoly(lam).as_expr(), lam)
    q = cp
    m0 = 0
    while q.degree() > 0 and q.eval(0) == 0:
        q = q.quo(Poly(lam, lam))
        m0 += 1
    if q.degree() > 0:
        g = gcd(q.as_expr(), diff(q.as_expr(), lam))
        coll = Poly(g, lam).degree() > 0
    else:
        coll = False
    return m0, coll, factor(cp.as_expr())


# ---------------------------------------------------------------- enumeration of cases
def make_cases():
    cases = []

    def sigs(n):
        return [(p, n - p) for p in range(n, -1, -1)]

    for (p, q) in sigs(4):
        for fam, prof in [("GEN", (1, 2)), ("COLL", (1, 1)), ("CRACK", (1, 0)), ("CRACK2", (0, 0))]:
            cases.append((4, p, q, fam, prof))
    for (p, q) in sigs(3):
        for fam, prof in [("GEN", (1,)), ("CRACK", (0,))]:
            cases.append((3, p, q, fam, prof))
    for (p, q) in sigs(5):
        for fam, prof in [("GEN", (1, 2)), ("COLL", (1, 1))]:
            cases.append((5, p, q, fam, prof))
    for (p, q) in [(3, 3), (5, 1)]:
        for fam, prof in [("GEN", (1, 2, 3)), ("COLL", (1, 1, 2)), ("CRACK", (1, 2, 0))]:
            cases.append((6, p, q, fam, prof))
    return cases


def process_case(n, p, q, fam, prof):
    eta = eta_matrix(p, q)
    Om = build_omega(n, prof)
    A = eta * Om  # eta^{-1} = eta
    # A must lie in so(p,q)
    if not (A * eta + eta * A.T).is_zero_matrix:
        fail("A = eta^{-1}Omega not in so(%d,%d) [%s %s]" % (p, q, fam, prof))
    pf_rec = pfaffian(Om)
    detOm = Om.det()
    rkOm = Om.rank()
    if n % 2 == 0:
        pf_block = Integer(1)
        for w in prof:
            pf_block *= Integer(w)
    else:
        pf_block = Integer(0)
    basis = so_basis(p, q)
    cb = centralizer(basis, A)
    # sanity: every element of the centralizer really commutes with A
    for Z in cb:
        if not (Z * A - A * Z).is_zero_matrix:
            fail("an element of c(A) does not commute with A [%s (%d,%d) %s]" % (fam, p, q, prof))
    closed = closed_under_bracket(cb)
    abelian, ddim = bracket_stats(cb)
    m0, coll, cpf = eig_info(A)
    rank_drop = rkOm < 2 * (n // 2)
    return dict(n=n, p=p, q=q, fam=fam, prof=prof, pf_rec=pf_rec, pf_block=pf_block,
                det=detOm, rk=rkOm, dimso=len(basis), dimc=len(cb), closed=closed,
                abelian=abelian, ddim=ddim, m0=m0, coll=coll, cpf=cpf,
                rank_drop=rank_drop, degenerate=(coll or rank_drop),
                kernel=n - rkOm)


def table_row(r):
    return "| %d | %d | %d | %-6s | %-9s | %d | %-3s | %-3s | %2d | %-3s | %2d |" % (
        r["n"], r["p"], r["q"], r["fam"], str(tuple(r["prof"])), r["rk"],
        ("yes" if r["pf_rec"] == 0 else "no"),
        ("yes" if r["coll"] else "no"),
        r["dimc"],
        ("yes" if r["abelian"] else "NO"),
        r["ddim"])


TABLE_HEAD = ("| n | p | q | family | profile   | rkOm | Pf=0 | coll | dim c | abel | dim[c,c] |\n"
              "|---|---|---|--------|-----------|------|------|-------|-------|------|----------|")


def print_table(recs):
    print("| n | p | q | family | profile   | rkOm | Pf=0 | coll | dimc | abel | dim[c,c] |")
    print("|---|---|---|--------|-----------|------|------|-------|------|------|----------|")
    for r in recs:
        print(table_row(r))


# ---------------------------------------------------------------- S2: isometries
def rot2():
    return Matrix([[Rational(3, 5), Rational(4, 5)], [Rational(-4, 5), Rational(3, 5)]])


def boost2():
    return Matrix([[Rational(5, 4), Rational(3, 4)], [Rational(3, 4), Rational(5, 4)]])


def embed2(n, i, B):
    S = eye(n)
    S[i, i] = B[0, 0]
    S[i, i + 1] = B[0, 1]
    S[i + 1, i] = B[1, 0]
    S[i + 1, i + 1] = B[1, 1]
    return S


def isometry_for(p, q):
    """Rational S in O(p,q) for n=4: rotation on same-sign pairs,
    boost on the mixed-sign pair."""
    if (p, q) == (3, 1):
        # axes (0,1): (+,+) -> rotation; axes (2,3): (+,-) -> boost
        return embed2(4, 0, rot2()) * embed2(4, 2, boost2()), "rot(0,1)*boost(2,3)"
    if (p, q) == (2, 2):
        # rotation on the (+,+)-pair (0,1) and on the (-,-)-pair (2,3), boost on the (+,-)-pair (1,2)
        return (embed2(4, 0, rot2()) * embed2(4, 2, rot2()) * embed2(4, 1, boost2()),
                "rot(0,1)*rot(2,3)*boost(1,2)")
    raise ValueError("signature has no isometry recipe")


# ================================================================ MAIN
def main():
    print("S918 (w29): residual structure of Omega in so(p,q) — centralizer c(A), A = eta^{-1}*Omega")
    print("Exact sympy arithmetic; no float/random/time.")
    print()
    print("SAMPLING BOUNDARY (explicit): n=4 — ALL signatures x {GEN,COLL,CRACK,CRACK2};")
    print("  n=3 — all signatures x {(1),(0)}; n=5 — all signatures x {(1,2),(1,1)};")
    print("  n=6 — ONLY signatures (3,3) and (5,1), families GEN/COLL/CRACK (sampling deliberately limited).")
    print()

    cases = make_cases()
    recs = [process_case(*c) for c in cases]
    verdicts = {}

    # ---------------------------------------------------------- S0
    print("=" * 100)
    print("SECTION S0 — SANITY: Pf(Omega)^2 = det(Omega) (Pf by recursion + block formula)")
    print("=" * 100)
    s0_checked = 0
    s0_bad = 0
    for r in recs:
        ok1 = (r["pf_rec"] == r["pf_block"])
        ok2 = (r["pf_rec"] ** 2 == r["det"])
        ok3 = True
        if r["n"] % 2 == 1:
            ok3 = (r["det"] == 0 and r["pf_rec"] == 0)
        ok = ok1 and ok2 and ok3
        s0_checked += 1
        if not ok:
            s0_bad += 1
            fail("S0: (%d,%d) %s %s: Pf_rec=%s Pf_block=%s det=%s" % (
                r["p"], r["q"], r["fam"], r["prof"], r["pf_rec"], r["pf_block"], r["det"]))
        print("  (%d,%d) %-6s profile=%-9s: Pf=%s, det=%s, rank=%d  [Pf_rec==Pf_block: %s; Pf^2==det: %s%s]" % (
            r["p"], r["q"], r["fam"], str(tuple(r["prof"])), r["pf_rec"], r["det"], r["rk"],
            "ok" if ok1 else "FAIL", "ok" if ok2 else "FAIL",
            ("; odd n: det==0: " + ("ok" if ok3 else "FAIL")) if r["n"] % 2 == 1 else ""))
    print("  Omega checked: %d (zero => invalid)" % s0_checked)
    v = "OK" if (s0_checked > 0 and s0_bad == 0) else "FAIL"
    verdicts["S0"] = "%s — %d/%d Omega passed Pf^2=det (recursion agreed with the block formula everywhere)" % (
        v, s0_checked - s0_bad, s0_checked)
    print("  VERDICT S0: " + verdicts["S0"])

    # ---------------------------------------------------------- S1
    print()
    print("=" * 100)
    print("SECTION S1 — CENTRALIZERS BY STRATUM (signature x family)")
    print("=" * 100)
    s1_checked = 0
    s1_bad = 0
    for r in recs:
        s1_checked += 1
        line = ("  n=%d (%d,%d) %-6s profile=%-9s: dim so=%d (=n(n-1)/2 ok), dim c(A)=%d, "
                "closed under [,]: %s, abelian: %s, dim[c,c]=%d") % (
            r["n"], r["p"], r["q"], r["fam"], str(tuple(r["prof"])), r["dimso"], r["dimc"],
            "yes" if r["closed"] else "NO", "yes" if r["abelian"] else "NO", r["ddim"])
        print(line)
        print("      charpoly(A) = %s  [multiplicity of zero m0=%d, collision of nonzero: %s; rank Omega=%d%s]" % (
            r["cpf"], r["m0"], "yes" if r["coll"] else "no", r["rk"],
            ", RANK DROP" if r["rank_drop"] else ""))
        if not r["closed"]:
            s1_bad += 1
            fail("S1: c(A) not closed under the bracket: (%d,%d) %s %s" % (r["p"], r["q"], r["fam"], r["prof"]))
        if r["dimso"] != r["n"] * (r["n"] - 1) // 2:
            s1_bad += 1
    print()
    print("  SUMMARY TABLE S1:")
    print_table(recs)
    # pattern: non-abelian <=> degenerate (collision OR rank drop)
    viol_A = [r for r in recs if (not r["abelian"]) and (not r["degenerate"])]  # non-abelian => degenerate
    viol_B = [r for r in recs if r["degenerate"] and r["abelian"]]              # degenerate => non-abelian
    print()
    print("  Checking the pattern \"non-abelian <=> degenerate (collision of eigenvalues of A OR rank drop of Omega)\":")
    print("    direction 1 (non-abelian => degenerate): violations %d" % len(viol_A))
    for r in viol_A:
        print("      BREAKS: n=%d (%d,%d) %s %s — non-abelian WITHOUT degeneracy" % (
            r["n"], r["p"], r["q"], r["fam"], tuple(r["prof"])))
    print("    direction 2 (degenerate => non-abelian): violations %d" % len(viol_B))
    for r in viol_B:
        print("      BREAKS: n=%d (%d,%d) %-6s %-9s — degenerate (coll=%s, rank_drop=%s), but ABELIAN (dim c=%d)" % (
            r["n"], r["p"], r["q"], r["fam"], str(tuple(r["prof"])), r["coll"], r["rank_drop"], r["dimc"]))
    # pattern refined from the data
    viol_A2 = [r for r in recs if (not r["abelian"]) and not (r["coll"] or r["kernel"] >= 3)]
    viol_B2 = [r for r in recs if (r["coll"] or r["kernel"] >= 3) and r["abelian"]]
    print("    REFINEMENT from the data: \"non-abelian <=> (collision of NONzero eigenvalues OR dim ker(Omega) >= 3)\":")
    print("      violations direction 1: %d, direction 2: %d" % (len(viol_A2), len(viol_B2)))
    for r in viol_A2 + viol_B2:
        print("      BREAKS the refinement: n=%d (%d,%d) %s %s (coll=%s, ker=%d, abel=%s)" % (
            r["n"], r["p"], r["q"], r["fam"], tuple(r["prof"]), r["coll"], r["kernel"], r["abelian"]))
    nonab = sum(1 for r in recs if not r["abelian"])
    if s1_checked > 0 and s1_bad == 0:
        vtext = ("OK (%d cases, %d non-abelian). The pattern from the statement does NOT hold as an equivalence: "
                 "the direction \"non-abelian => degenerate\" holds (violations %d), but \"degenerate => non-abelian\" "
                 "BREAKS on %d cases (a rank drop with ker=2 gives an ABELIAN c). "
                 "The exact pattern in the data: non-abelian <=> (collision of nonzero eigenvalues of A "
                 "OR dim ker >= 3): violations %d+%d.") % (
            s1_checked, nonab, len(viol_A), len(viol_B), len(viol_A2), len(viol_B2))
        if len(viol_A) > 0:
            vtext = ("OK-with-violations (%d cases): both directions have counterexamples: %d and %d." %
                     (s1_checked, len(viol_A), len(viol_B)))
        if len(viol_A2) + len(viol_B2) > 0:
            vtext += " NOTE: the refined pattern also has counterexamples."
    else:
        vtext = "FAIL (closure/dimensions broken)"
    verdicts["S1"] = vtext
    print("  VERDICT S1: " + verdicts["S1"])

    # ---------------------------------------------------------- S2
    print()
    print("=" * 100)
    print("SECTION S2 — INVARIANCE UNDER AN O(p,q) ISOMETRY (rational rotations/boosts)")
    print("=" * 100)
    s2_checked = 0
    s2_bad = 0
    for (p, q) in [(2, 2), (3, 1)]:
        eta = eta_matrix(p, q)
        S, Sdesc = isometry_for(p, q)
        iso_ok = (S.T * eta * S - eta).is_zero_matrix
        print("  (%d,%d): S = %s; checking S^T*eta*S = eta: %s" % (p, q, Sdesc, "ok" if iso_ok else "FAIL"))
        if not iso_ok:
            fail("S2: S is not an isometry in (%d,%d)" % (p, q))
            s2_bad += 1
            continue
        basis = so_basis(p, q)
        for fam, prof in [("GEN", (1, 2)), ("CRACK", (1, 0))]:
            Om = build_omega(4, prof)
            A = eta * Om
            Om2 = S.T * Om * S
            if not (Om2.T + Om2).is_zero_matrix:
                fail("S2: Omega' is not antisymmetric (%d,%d) %s" % (p, q, fam))
            A2 = eta * Om2
            cb1 = centralizer(basis, A)
            cb2 = centralizer(basis, A2)
            ab1 = bracket_stats(cb1)[0]
            ab2 = bracket_stats(cb2)[0]
            same = (len(cb1) == len(cb2)) and (ab1 == ab2)
            s2_checked += 1
            if not same:
                s2_bad += 1
                fail("S2: invariance broken (%d,%d) %s" % (p, q, fam))
            print("    %-6s profile=%-7s: dim c(A)=%d, dim c(A')=%d; abelian: %s / %s  => %s" % (
                fam, str(tuple(prof)), len(cb1), len(cb2),
                "yes" if ab1 else "NO", "yes" if ab2 else "NO",
                "MATCH" if same else "MISMATCH"))
    print("  Pairs checked (Omega, Omega'): %d (zero => invalid)" % s2_checked)
    v = "OK" if (s2_checked > 0 and s2_bad == 0) else "FAIL"
    verdicts["S2"] = "%s — dim c and abelian-ness are invariant under the exact isometry in %d/%d cases" % (
        v, s2_checked - s2_bad, s2_checked)
    print("  VERDICT S2: " + verdicts["S2"])

    # ---------------------------------------------------------- M: mutants
    print()
    print("=" * 100)
    print("SECTION M — MUTANTS (each MUST be caught, otherwise the measurement is invalid)")
    print("=" * 100)
    mut_caught = 0
    mut_total = 3

    # m1: broken centralizer solver (dropping equations from the end of the system)
    print("  m1: centralizer solver with dropped equations (from the end of the system).")
    m1_cases = [(4, 4, 0, "GEN", (1, 2)), (4, 2, 2, "GEN", (1, 2)), (4, 3, 1, "CRACK", (1, 0))]
    m1_hits = 0
    for (n, p, q, fam, prof) in m1_cases:
        eta = eta_matrix(p, q)
        Om = build_omega(n, prof)
        A = eta * Om
        basis = so_basis(p, q)
        cb_true = centralizer(basis, A)
        # minimal k >= 1 for which dropping the last k equations changes the solution
        k = 1
        cb_bad = centralizer(basis, A, drop_last_rows=k)
        while len(cb_bad) == len(cb_true) and k < n * n:
            k += 1
            cb_bad = centralizer(basis, A, drop_last_rows=k)
        dim_diff = (len(cb_bad) != len(cb_true))
        closed_bad = closed_under_bracket(cb_bad)
        resid = sum(1 for Z in cb_bad if not (Z * A - A * Z).is_zero_matrix)
        caught = dim_diff or (not closed_bad) or (resid > 0)
        if caught:
            m1_hits += 1
        print("    (%d,%d) %-6s %-7s: dropped last equations k=%d; dim: %d (true) vs %d (broken) [%s]; "
              "broken closed under [,]: %s%s; elements with [X,A]!=0: %d => %s" % (
                  p, q, fam, str(tuple(prof)), k, len(cb_true), len(cb_bad),
                  "DIFFERENT" if dim_diff else "match",
                  "yes" if closed_bad else "NO",
                  " (CATCHES)" if not closed_bad else "",
                  resid, "CAUGHT" if caught else "NOT-CAUGHT"))
    m1_ok = m1_hits > 0
    if m1_ok:
        mut_caught += 1
    else:
        fail("m1 NOT caught — the measurement is invalid")
    print("    m1: caught in %d/%d cases => %s" % (m1_hits, len(m1_cases), "CAUGHT" if m1_ok else "NOT-CAUGHT"))

    # m2: broken abelian-ness test (always "abelian") on the so(3) reference
    print("  m2: an abelian-ness test that always says \"abelian\", on the so(3) reference in (3,0).")
    J12 = zeros(3, 3); J12[0, 1] = Integer(1); J12[1, 0] = Integer(-1)
    J13 = zeros(3, 3); J13[0, 2] = Integer(1); J13[2, 0] = Integer(-1)
    J23 = zeros(3, 3); J23[1, 2] = Integer(1); J23[2, 1] = Integer(-1)
    so3 = [J12, J13, J23]
    eta3 = eta_matrix(3, 0)
    for J in so3:
        if not (J * eta3 + eta3 * J.T).is_zero_matrix:
            fail("m2: the reference J is not in so(3)")
    true_ab = bracket_stats(so3)[0]      # the real test
    mutant_ab = True                      # mutant: always "abelian"
    m2_ok = (true_ab is False) and (mutant_ab is True)
    if m2_ok:
        mut_caught += 1
    else:
        fail("m2 NOT caught — the measurement is invalid")
    print("    real test: is so(3) abelian? %s (must be NO); mutant: %s => discrepancy %s" % (
        "yes" if true_ab else "NO", "yes" if mutant_ab else "NO", "CAUGHT" if m2_ok else "NOT-CAUGHT"))

    # m3: false claim Pf = w1 + w2 instead of the product, n=4, profile (1,2)
    print("  m3: the false claim Pf(Omega) = w1 + w2 for a block Omega, n=4, profile (1,2).")
    Om = build_omega(4, (1, 2))
    pf_fake = Integer(1) + Integer(2)
    d = Om.det()
    m3_ok = (pf_fake ** 2 != d) and (pfaffian(Om) ** 2 == d)
    if m3_ok:
        mut_caught += 1
    else:
        fail("m3 NOT caught — the measurement is invalid")
    print("    Pf_fake = %s, Pf_fake^2 = %s, det = %s: %s != %s => %s (real Pf=%s, Pf^2=det: ok)" % (
        pf_fake, pf_fake ** 2, d, pf_fake ** 2, d, "CAUGHT" if m3_ok else "NOT-CAUGHT", pfaffian(Om)))

    verdicts["M"] = "%s — caught %d/%d mutants" % ("OK" if mut_caught == mut_total else "FAIL", mut_caught, mut_total)
    print("  VERDICT M: " + verdicts["M"])

    # ---------------------------------------------------------- SUMMARY
    print()
    print("=" * 100)
    print("SECTION SUMMARY")
    print("=" * 100)
    print("  VERDICTS:")
    for k in ["S0", "S1", "S2", "M"]:
        print("    %s: %s" % (k, verdicts[k]))
    print()
    print("  COUNTERS: S0=%d Omega checked; S1=%d strata; S2=%d isometry pairs; M=%d/%d mutants caught."
          % (s0_checked, s1_checked, s2_checked, mut_caught, mut_total))
    print()
    print("  FULL TABLE S1 (repeated):")
    print_table(recs)
    print()

    counters_ok = (s0_checked > 0 and s1_checked > 0 and s2_checked > 0 and mut_caught == mut_total)
    if not counters_ok:
        fail("counters: zero checked, or a mutant not caught")
    if FAIL:
        print("SUMMARY: FAIL (%d issues):" % len(FAIL))
        for m in FAIL:
            print("  - " + m)
        print("EXIT: 1")
        return 1
    print("SUMMARY: OK — all sections passed, all mutants caught, counters nonzero.")
    print("EXIT: 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
