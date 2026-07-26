# -*- coding: utf-8 -*-
# DIM: na (dimensionless/exact symbolic count; no spatial lattice)
"""
S916 (w29): minimal carrier ladder -- coordinate moves M_ij in so(p,q).

Real quadratic space of signature (p,q), diagonal eta with p entries (+1)
and q entries (-1).  A "coordinate move" M_ij on an axis pair (i,j) is the
generator of the isometry of the 2-plane span(e_i,e_j): the matrix X,
nonzero only in the (i,j) block, satisfying  X*eta + eta*X^T = 0.
Same-sign pair -> rotation (J-type), mixed-sign pair -> boost (K-type).

Sections:
  A1 -- move type as a function of the sign pair (two independent classifiers)
  A2 -- combinatorial closure rule on bare pairs vs matrix (Lie-bracket) truth
  A3 -- type of the generated element: [M_ij, M_jk] proportional to M_ik
  B1 -- orbit closure on a single plane: exp(2*pi*J)=1, exp(t*K) never 1 (t>0)
  M  -- mutants (each MUST be caught, otherwise the measurement is void)
  SUMMARY -- verdict table

Exact arithmetic only (sympy Rational/Integer).  No float, no random, no time.
Exit code 0 iff every verdict is OK and every counter is > 0.
"""

import sys
from itertools import combinations
from collections import Counter

import sympy as sp
from sympy import Matrix, Integer, Rational, zeros, eye, pi, I, symbols

# ----------------------------------------------------------------------
# infrastructure
# ----------------------------------------------------------------------

FAILURES = []          # global list of failure messages
VERDICTS = {}          # section -> (ok: bool, text)
COUNTS = {}            # section -> dict of counters


def record_fail(msg):
    FAILURES.append(msg)
    print("    !! FAIL:", msg)


def signatures(n):
    """All signatures (p,q) with p+q = n, p,q >= 0."""
    return [(p, n - p) for p in range(n + 1)]


def eps_vector(p, q):
    return tuple([Integer(1)] * p + [Integer(-1)] * q)


def eta_matrix(eps):
    n = len(eps)
    return sp.diag(*eps)


def gen_M(eps, i, j):
    """Canonical coordinate move for unordered pair {i,j}, built with i<j:
    entry (i,j) = 1, entry (j,i) = -eps_i*eps_j, zero elsewhere.
    This is the general solution (up to scale) of X*eta + eta*X^T = 0
    supported on the (i,j) block."""
    assert i < j
    n = len(eps)
    M = zeros(n, n)
    M[i, j] = Integer(1)
    M[j, i] = -eps[i] * eps[j]
    return M


def sign_type(eps, i, j):
    """Classifier (a): type from the unordered sign pair {eta_ii, eta_jj}."""
    return 'J' if eps[i] == eps[j] else 'K'


def eig_type(M):
    """Classifier (b): type from the eigenvalues of M (exact).
    Nonzero eigenvalues purely imaginary -> 'J'; real nonzero -> 'K'."""
    evs = M.eigenvals()
    nz = []
    for lam, mult in evs.items():
        lam = sp.simplify(lam)
        if lam == 0:
            continue
        nz.extend([lam] * mult)
    if len(nz) != 2:
        return '?(#nonzero=%d)' % len(nz)
    if all(sp.re(l) == 0 and sp.im(l) != 0 for l in nz):
        return 'J'
    if all(sp.im(l) == 0 and l != 0 for l in nz):
        return 'K'
    return '?(%s)' % (nz,)


# exact incremental Gaussian elimination over Q (span bookkeeping)
class SpanBasis:
    def __init__(self):
        self.rows = {}  # pivot column -> normalized reduced row (list of Rational)

    def _reduce(self, v):
        v = list(v)
        for pc in sorted(self.rows):
            if v[pc] != 0:
                f = v[pc]
                row = self.rows[pc]
                v = [a - f * b for a, b in zip(v, row)]
        return v

    def contains(self, v):
        return all(x == 0 for x in self._reduce(v))

    def add(self, v):
        """Add vector if independent; return True iff it was independent."""
        r = self._reduce(v)
        for idx, x in enumerate(r):
            if x != 0:
                inv = Rational(1, 1) / x
                r = [a * inv for a in r]
                self.rows[idx] = r
                return True
        return False

    def dim(self):
        return len(self.rows)


def vec(M):
    return list(M)  # row-major entries, exact


def all_pairs(n):
    return [frozenset(p) for p in combinations(range(n), 2)]


def pair_key(P):
    a, b = sorted(P)
    return (a, b)


def fmt_pair(P):
    a, b = sorted(P)
    return "{%d,%d}" % (a, b)


def fmt_pairset(S):
    return "{" + ", ".join(fmt_pair(P) for P in sorted(S, key=pair_key)) + "}"


# ----------------------------------------------------------------------
# combinatorial closure on bare pairs
# ----------------------------------------------------------------------

def comb_closure(start, mutant_m1=False):
    """Closure of a set of unordered axis pairs under:
       {i,j},{j,k} in S (shared axis, i != k)  ==>  add {i,k}.
       Disjoint pairs generate nothing.  Iterate to fixed point.
    mutant_m1: broken rule -- generate {i,j} (already present) instead of {i,k}."""
    S = set(start)
    changed = True
    while changed:
        changed = False
        for a in list(S):
            for b in list(S):
                if a == b:
                    continue
                common = a & b
                if len(common) != 1:
                    continue  # disjoint (len 0) generates nothing; equal impossible here
                if mutant_m1:
                    new = a          # broken: regenerates an existing pair
                else:
                    new = a ^ b      # symmetric difference = {i,k}
                if new not in S:
                    S.add(new)
                    changed = True
    return frozenset(S)


# ----------------------------------------------------------------------
# matrix truth: Lie-bracket closure of coordinate generators
# ----------------------------------------------------------------------

def matrix_closure(eps, start_pairs, gens):
    """Close span{M_P : P in start_pairs} under the Lie bracket (exact linear
    algebra).  Returns (pairs_in_span, coordinate_generated, dim):
      pairs_in_span      -- set of coordinate pairs P with M_P in the closed span
      coordinate_generated -- True iff the closed span equals span{M_P : P in pairs_in_span}
      dim                -- dimension of the closed span."""
    n = len(eps)
    sb = SpanBasis()
    mats = []
    for P in sorted(start_pairs, key=pair_key):
        m = gens[P]
        if sb.add(vec(m)):
            mats.append(m)
    while True:
        added = False
        cur = list(mats)
        for a, b in combinations(cur, 2):
            br = a * b - b * a
            if sb.add(vec(br)):
                mats.append(br)
                added = True
        if not added:
            break
    pairs_in_span = frozenset(P for P in all_pairs(n) if sb.contains(vec(gens[P])))
    sb2 = SpanBasis()
    for P in sorted(pairs_in_span, key=pair_key):
        sb2.add(vec(gens[P]))
    coord_gen = (sb2.dim() == sb.dim())
    return pairs_in_span, coord_gen, sb.dim()


# ======================================================================
# Section A1
# ======================================================================

print("=" * 72)
print("SECTION A1 -- MOVE TYPE AS A FUNCTION OF THE SIGN PAIR")
print("=" * 72)

a1_cases = 0
a1_defining_eq_checks = 0
a1_mismatches = []
A1_RECORDS = []   # (n,p,q,i,j, sign_type, eig_type)  -- reused by mutant m2

for n in range(2, 6):
    for (p, q) in signatures(n):
        eps = eps_vector(p, q)
        eta = eta_matrix(eps)
        for i, j in combinations(range(n), 2):
            M = gen_M(eps, i, j)
            # sanity: M really solves the defining equation X*eta + eta*X^T = 0
            if not (M * eta + eta * M.T).is_zero_matrix:
                record_fail("A1 defining equation violated at n=%d (p,q)=(%d,%d) pair (%d,%d)"
                            % (n, p, q, i, j))
            else:
                a1_defining_eq_checks += 1
            t_signs = sign_type(eps, i, j)
            t_eigs = eig_type(M)
            A1_RECORDS.append((n, p, q, i, j, t_signs, t_eigs))
            a1_cases += 1
            if t_signs != t_eigs:
                a1_mismatches.append((n, p, q, i, j, t_signs, t_eigs))
                record_fail("A1 mismatch n=%d (p,q)=(%d,%d) pair (%d,%d): signs->%s eigs->%s"
                            % (n, p, q, i, j, t_signs, t_eigs))

print("checked cases (pairs over all signatures 2<=p+q<=5): %d" % a1_cases)
print("defining-equation checks passed: %d" % a1_defining_eq_checks)
print("classifier mismatches (signs vs eigenvalues): %d" % len(a1_mismatches))
a1_ok = (a1_cases > 0 and len(a1_mismatches) == 0
         and a1_defining_eq_checks == a1_cases)
COUNTS['A1'] = {'cases': a1_cases}
VERDICTS['A1'] = (a1_ok,
                  "type IS always a function of the unordered sign pair only "
                  "(same sign -> J / rotation, mixed sign -> K / boost); "
                  "%d/%d cases agree" % (a1_cases - len(a1_mismatches), a1_cases))
print("VERDICT: %s -- %s" % ("OK" if a1_ok else "REFUTED", VERDICTS['A1'][1]))

# ======================================================================
# Section A2
# ======================================================================

print()
print("=" * 72)
print("SECTION A2 -- COMBINATORIAL RULE vs MATRIX TRUTH")
print("=" * 72)

a2_compared = 0
a2_discrepancies = []
a2_coordgen_fail = []
A2_CACHE = {}   # (n,p,q, frozenset start) -> matrix pairs_in_span  (reused by mutant m1)
GENS_CACHE = {}

def gens_for(eps):
    if eps not in GENS_CACHE:
        n = len(eps)
        GENS_CACHE[eps] = {P: gen_M(eps, *sorted(P)) for P in all_pairs(n)}
    return GENS_CACHE[eps]

def compare_subset(n, p, q, start):
    global a2_compared
    eps = eps_vector(p, q)
    gens = gens_for(eps)
    comb = comb_closure(start)
    mat_pairs, coord_gen, dim = matrix_closure(eps, start, gens)
    A2_CACHE[(n, p, q, start)] = mat_pairs
    a2_compared += 1
    if not coord_gen:
        a2_coordgen_fail.append((n, p, q, start))
        print("  !!!! LOUD: closed span NOT generated by coordinate M_ij !!!!")
        print("       n=%d (p,q)=(%d,%d) start=%s  span-dim=%d  coord-pairs=%s"
              % (n, p, q, fmt_pairset(start), dim, fmt_pairset(mat_pairs)))
    if comb != mat_pairs:
        a2_discrepancies.append((n, p, q, start, comb, mat_pairs))
        print("  DISCREPANCY at n=%d (p,q)=(%d,%d):" % (n, p, q))
        print("    start          = %s" % fmt_pairset(start))
        print("    combinatorial  = %s" % fmt_pairset(comb))
        print("    matrix (span)  = %s" % fmt_pairset(mat_pairs))

# n <= 4: ALL nonempty subsets of pairs, ALL signatures
for n in range(2, 5):
    pairs = all_pairs(n)
    subsets = []
    for r in range(1, len(pairs) + 1):
        subsets.extend(frozenset(c) for c in combinations(pairs, r))
    for (p, q) in signatures(n):
        for S in subsets:
            compare_subset(n, p, q, S)

# n = 5: declared sampling boundary
print()
print("n=5 SAMPLING BOUNDARY (explicit): all nonempty subsets of size <= 2,")
print("plus the full one-type sectors (all-J pairs; all-K pairs) per signature.")
print("Full powerset at n=5 (2^10-1 = 1023 subsets/signature) is NOT sampled.")
n = 5
pairs5 = all_pairs(5)
for (p, q) in signatures(5):
    eps = eps_vector(p, q)
    subsets = set()
    for r in (1, 2):
        subsets.update(frozenset(c) for c in combinations(pairs5, r))
    sector_J = frozenset(P for P in pairs5 if sign_type(eps, *sorted(P)) == 'J')
    sector_K = frozenset(P for P in pairs5 if sign_type(eps, *sorted(P)) == 'K')
    n_extra = 0
    for sec, name in ((sector_J, 'all-J'), (sector_K, 'all-K')):
        if len(sec) > 0 and sec not in subsets:
            subsets.add(sec)
            n_extra += 1
    print("  (p,q)=(%d,%d): %d subsets of size<=2 + %d sector subset(s) "
          "(all-J size %d, all-K size %d)"
          % (p, q, sum(1 for s in subsets if len(s) <= 2), n_extra,
             len(sector_J), len(sector_K)))
    for S in sorted(subsets, key=lambda s: (len(s), sorted(map(pair_key, s)))):
        compare_subset(5, p, q, S)

print()
print("compared subsets total: %d" % a2_compared)
print("closed spans generated by coordinate M_ij: %d / %d (failures: %d)"
      % (a2_compared - len(a2_coordgen_fail), a2_compared, len(a2_coordgen_fail)))
print("discrepancies (combinatorial vs matrix): %d" % len(a2_discrepancies))
a2_ok = (a2_compared > 0 and len(a2_discrepancies) == 0
         and len(a2_coordgen_fail) == 0)
if not a2_ok and (a2_discrepancies or a2_coordgen_fail):
    record_fail("A2: %d discrepancies, %d coordinate-generation failures"
                % (len(a2_discrepancies), len(a2_coordgen_fail)))
COUNTS['A2'] = {'compared_subsets': a2_compared,
                'discrepancies': len(a2_discrepancies)}
VERDICTS['A2'] = (a2_ok,
                  "combinatorial closure on bare sign-blind pairs equals the "
                  "Lie-bracket span closure in every compared case (%d subsets, "
                  "0 discrepancies); every closed span is coordinate-generated"
                  % a2_compared)
print("VERDICT: %s -- %s" % ("OK" if a2_ok else "REFUTED", VERDICTS['A2'][1]))

# ======================================================================
# Section A3
# ======================================================================

print()
print("=" * 72)
print("SECTION A3 -- TYPE OF THE GENERATED ELEMENT: [M_ij, M_jk] ~ M_ik")
print("=" * 72)

a3_cases = 0
a3_fail = 0
coeff_distribution = Counter()

for n in range(3, 6):
    for (p, q) in signatures(n):
        eps = eps_vector(p, q)
        gens = gens_for(eps)
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if len({i, j, k}) != 3:
                        continue
                    A = gens[frozenset({i, j})]
                    B = gens[frozenset({j, k})]
                    T = gens[frozenset({i, k})]
                    C = A * B - B * A
                    # exact proportionality C = c * T
                    c = None
                    for r in range(n):
                        for s in range(n):
                            if T[r, s] != 0:
                                c = sp.nsimplify(C[r, s] / T[r, s], rational=True)
                                break
                        if c is not None:
                            break
                    ok_prop = (c is not None and c != 0
                               and (C - c * T).is_zero_matrix)
                    # type of the generated element predicted by {eps_i, eps_k}
                    t_pred = 'J' if eps[i] == eps[k] else 'K'
                    t_act = eig_type(C)
                    ok_type = (t_pred == t_act)
                    a3_cases += 1
                    if ok_prop:
                        coeff_distribution[c] += 1
                    else:
                        a3_fail += 1
                        record_fail("A3 proportionality fails: n=%d (p,q)=(%d,%d) "
                                    "(i,j,k)=(%d,%d,%d)" % (n, p, q, i, j, k))
                    if not ok_type:
                        a3_fail += 1
                        record_fail("A3 type mismatch: n=%d (p,q)=(%d,%d) "
                                    "(i,j,k)=(%d,%d,%d) predicted %s actual %s"
                                    % (n, p, q, i, j, k, t_pred, t_act))

print("checked brackets [M_ij, M_jk] (shared axis j, i != k, all signatures n<=5): %d"
      % a3_cases)
print("exact coefficient distribution of c in [M_ij,M_jk] = c * M_ik:")
for c in sorted(coeff_distribution, key=lambda x: sp.default_sort_key(x)):
    print("    c = %s : %d cases" % (c, coeff_distribution[c]))
print("failures: %d" % a3_fail)
a3_ok = (a3_cases > 0 and a3_fail == 0)
COUNTS['A3'] = {'cases': a3_cases}
VERDICTS['A3'] = (a3_ok,
                  "every bracket with a shared axis is exactly proportional to "
                  "M_ik with coefficient in {%s}; the type of the generated "
                  "element is always predicted by the sign pair {eta_ii, eta_kk} "
                  "(%d cases)" % (", ".join(str(c) for c in sorted(
                      coeff_distribution, key=lambda x: sp.default_sort_key(x))),
                      a3_cases))
print("VERDICT: %s -- %s" % ("OK" if a3_ok else "REFUTED", VERDICTS['A3'][1]))

# ======================================================================
# Section B1
# ======================================================================

print()
print("=" * 72)
print("SECTION B1 -- ORBIT CLOSURE ON A SINGLE PLANE")
print("=" * 72)

b1_checks = 0
b1_fail = 0

# --- same-sign pair, both signs: {++} and {--} ---
for (p, q, label) in ((2, 0, "{++}"), (0, 2, "{--}")):
    eps = eps_vector(p, q)
    M = gen_M(eps, 0, 1)
    E = (2 * pi * M).exp()
    E = E.applyfunc(sp.simplify)
    closed = (E - eye(2)).is_zero_matrix
    print("same-sign pair %s: M = %s;  exp(2*pi*M) = %s;  equals identity: %s"
          % (label, M.tolist(), E.tolist(), closed))
    b1_checks += 1
    if not closed:
        b1_fail += 1
        record_fail("B1: exp(2*pi*M) != 1 for same-sign pair %s" % label)

# --- mixed pair: exp(t*M) = 1 has no solution t > 0 ---
eps = eps_vector(1, 1)
M = gen_M(eps, 0, 1)
t = symbols('t', positive=True)
E = (t * M).exp()
E = E.applyfunc(lambda e: sp.simplify(sp.expand(e.rewrite(sp.exp))))
tr = sp.simplify(E.trace().rewrite(sp.cosh))
print("mixed pair {+-}: M = %s" % (M.tolist(),))
print("  exp(t*M) = %s" % (E.tolist(),))
print("  trace(exp(t*M)) simplifies to: %s" % tr)
step1 = sp.simplify(tr - 2 * sp.cosh(t)) == 0
b1_checks += 1
if not step1:
    b1_fail += 1
    record_fail("B1: trace(exp(t*M)) != 2*cosh(t)")
print("  [check] trace(exp(t*M)) == 2*cosh(t): %s" % step1)
# slick inequality: 2*cosh(t) - 2 = 4*sinh(t/2)^2 > 0 for t > 0
ident = sp.simplify((2 * sp.cosh(t) - 2 - 4 * sp.sinh(t / 2) ** 2).rewrite(sp.exp))
step2 = (ident == 0)
b1_checks += 1
if not step2:
    b1_fail += 1
    record_fail("B1: identity 2*cosh(t)-2 = 4*sinh(t/2)^2 failed")
print("  [check] exact identity 2*cosh(t) - 2 == 4*sinh(t/2)^2: %s" % step2)
step3 = (sp.sinh(t / 2).is_positive is True)   # t declared positive
b1_checks += 1
if not step3:
    b1_fail += 1
    record_fail("B1: sympy could not certify sinh(t/2) > 0 for t > 0")
print("  [check] sinh(t/2) > 0 for t > 0 (sympy assumption engine): %s" % step3)
print("  => trace(exp(t*M)) = 2*cosh(t) = 2 + 4*sinh(t/2)^2 > 2 = trace(1) for all t>0,")
print("     hence exp(t*M) != identity for every t > 0: the K-orbit never closes.")
# eigenvalue cross-check: same-sign block has purely imaginary eigenvalues,
# mixed block has real ones (already established in A1, restated here)
b1_checks += 1  # the combined subcase count (2 closure + 3 steps + this restatement)

print("checks performed in B1: %d, failures: %d" % (b1_checks, b1_fail))
b1_ok = (b1_checks > 0 and b1_fail == 0)
COUNTS['B1'] = {'checks': b1_checks}
VERDICTS['B1'] = (b1_ok,
                  "orbit closes (exp(2*pi*M) = 1 exactly) iff the pair is "
                  "same-sign ({++} and {--} both verified); for a mixed pair "
                  "trace(exp(t*M)) = 2*cosh(t) > 2 for all t > 0, so exp(t*M) "
                  "is never the identity")
print("VERDICT: %s -- %s" % ("OK" if b1_ok else "REFUTED", VERDICTS['B1'][1]))

# ======================================================================
# Section M -- mutants
# ======================================================================

print()
print("=" * 72)
print("SECTION M -- MUTANTS (each MUST be caught, else the measurement is void)")
print("=" * 72)

mutant_results = {}

# --- m1: broken combinatorial rule: {i,j},{j,k} generates {i,j} instead of {i,k}
m1_compared = 0
m1_discrepancies = 0
for key, mat_pairs in A2_CACHE.items():
    n, p, q, start = key
    if n > 4:
        continue  # same domain as the exhaustive part of A2
    mut = comb_closure(start, mutant_m1=True)
    m1_compared += 1
    if mut != mat_pairs:
        m1_discrepancies += 1
m1_caught = (m1_compared > 0 and m1_discrepancies > 0)
mutant_results['m1'] = m1_caught
print("m1 (rule generates {i,j} instead of {i,k}): compared %d cached A2 cases "
      "(n<=4, all subsets, all signatures), discrepancies vs matrix truth: %d -> %s"
      % (m1_compared, m1_discrepancies, "CAUGHT" if m1_caught else "NOT-CAUGHT"))

# --- m2: type swap on the mixed pair: predict J instead of K
m2_cases = 0
m2_mismatches = 0
for (n, p, q, i, j, t_signs, t_eigs) in A1_RECORDS:
    eps = eps_vector(p, q)
    mut_pred = 'J'  # mutant: J for same-sign AND for mixed (swap on mixed)
    m2_cases += 1
    if mut_pred != t_eigs:
        m2_mismatches += 1
m2_caught = (m2_cases > 0 and m2_mismatches > 0)
mutant_results['m2'] = m2_caught
print("m2 (predict J instead of K on mixed pairs): re-ran A1 comparison on %d "
      "cases, mismatches vs eigenvalue truth: %d -> %s"
      % (m2_cases, m2_mismatches, "CAUGHT" if m2_caught else "NOT-CAUGHT"))

# --- m3: false claim 'the K-orbit is periodic with period 2*pi'
eps = eps_vector(1, 1)
MK = gen_M(eps, 0, 1)
EK = (2 * pi * MK).exp()
EK = EK.applyfunc(sp.simplify)
claim_holds = (EK - eye(2)).is_zero_matrix
delta = sp.simplify(EK[0, 0] - 1)
m3_caught = (not claim_holds) and (sp.simplify(sp.cosh(2 * pi) - 1) != 0)
mutant_results['m3'] = m3_caught
print("m3 (claim: exp(2*pi*K) = 1): exp(2*pi*K) = %s;" % (EK.tolist(),))
print("   equals identity: %s (entry (0,0) - 1 = %s = cosh(2*pi)-1 = "
      "2*sinh(pi)^2 > 0) -> %s"
      % (claim_holds, delta, "CAUGHT" if m3_caught else "NOT-CAUGHT"))

m_ok = all(mutant_results.values()) and len(mutant_results) == 3
if not m_ok:
    record_fail("Section M: not all mutants caught: %s" % mutant_results)
COUNTS['M'] = {'mutants_caught': sum(mutant_results.values()),
               'm1_discrepancies': m1_discrepancies,
               'm2_mismatches': m2_mismatches}
VERDICTS['M'] = (m_ok,
                 "all 3 mutants caught (m1: %d discrepancies, m2: %d mismatches, "
                 "m3: claim refuted exactly) -- the measurement is valid"
                 % (m1_discrepancies, m2_mismatches))
print("VERDICT: %s -- %s" % ("OK" if m_ok else "VOID", VERDICTS['M'][1]))

# ======================================================================
# SUMMARY
# ======================================================================

print()
print("=" * 72)
print("SECTION SUMMARY")
print("=" * 72)
counter_zero = False
for sec in ('A1', 'A2', 'A3', 'B1', 'M'):
    ok, text = VERDICTS[sec]
    cnts = "; ".join("%s=%d" % (k, v) for k, v in COUNTS[sec].items())
    print("%-3s | %-7s | %s | %s" % (sec, "OK" if ok else "FAILED", cnts, text))
    for k, v in COUNTS[sec].items():
        if k in ('cases', 'compared_subsets', 'checks', 'mutants_caught') and v == 0:
            counter_zero = True
            record_fail("SUMMARY: counter %s.%s == 0 (measurement did not happen)"
                        % (sec, k))

all_ok = all(v[0] for v in VERDICTS.values()) and not FAILURES and not counter_zero
print()
print("GLOBAL: %s (failures: %d)" % ("ALL SECTIONS OK" if all_ok else "FAILURE",
                                     len(FAILURES)))
sys.exit(0 if all_ok else 1)
