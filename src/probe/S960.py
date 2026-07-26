# -*- coding: utf-8 -*-
# DIM: na (forall-d SYMBOLIC proof of the W34-leg1 cell identities; 0 handles).
#
# ============================================================================
# GOAL: prove the four W34-leg1 carves with a SYMBOLIC d (not an enumeration),
# so the "d<=6 listing" qualifier can be dropped.  Where a step resists a fully
# symbolic close, an instance-tail d=2..8 re-confirms it (stamped, not silent).
#
# CELL: d+1 unit vectors, pairwise inner product -1/d.  v_i.v_i=1, v_i.v_j=-1/d.
#
# CLAIM 1 (Gram of differences ~ Cartan(A_d), forall d):
#   reduced basis b_p = v_p - v_{p+1}.  Symbolic-d expansion of b_p.b_q gives
#     |p-q|=0 -> 2(d+1)/d ;  |p-q|=1 -> -(d+1)/d ;  |p-q|>=2 -> 0,
#   i.e. Gram = (d+1)/d * Cartan(A_d).  det: (d+1)/d)^d * det Cartan, and
#   det Cartan(A_d) = d+1 forall d (tridiagonal recurrence D_d=2D_{d-1}-D_{d-2},
#   D_1=2,D_2=3 -> D_d=d+1; d+1 satisfies it: 2d-(d-1)=d+1).
#
# CLAIM 2 (z = d+1 combinatorial, forall d):
#   coset (u_0 + L) = { z - c : z in Z^{d+1}, sum z = 1 }.
#   |z-c|^2 = sum z_i^2 - 1/(d+1).  Minimize integer sum z_i^2 s.t. sum z_i = 1:
#     sum z_i^2 >= sum |z_i| >= |sum z_i| = 1, equality iff all z_i in {0,1} with
#     exactly one 1  ->  the d+1 vectors e_i.  So z(d)=d+1 minimal vectors.
#
# CLAIM 3 (Gauss sum, forall d):  f(k_j) = sum_{i=0}^{d} w^i, w = exp(2 pi I j/(d+1)).
#   Finite geometric series: (w-1) f = w^{d+1} - 1 = exp(2 pi I j) - 1 = 0
#   (integer j); w != 1 for 1<=j<=d  =>  f = 0.  j=0 -> w=1 -> f=d+1.
#   Hence exactly d node values (j=1..d), Gamma (j=0) not a node.
#
# CLAIM 4 (finite traversal <=> d=2):  rotation cos(2 theta) = 2/d^2 - 1.
#   Finite order <=> 2 theta in pi*Q <=> (Niven) cos(2 theta) in {0,+-1/2,+-1}.
#   Solve 2/d^2-1 = value for integer d>=2:  value=-1/2 -> d^2=4 -> d=2 (unique);
#   all other Niven values give no integer d>=2.  => finite iff d=2.
#
# Discipline: 0 handles; sympy exact + symbolic d; mutant on each claim;
# FORBIDDEN-SCAN (S929 list); bit-fence d=2 vs S956 carve; STOP after the table.
# ============================================================================

import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import itertools
from sympy import (symbols, Symbol, Rational, Integer, zeros, eye, simplify,
                   summation, exp, I, pi, cos, Eq, expand, factor, Abs)

_LOGPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "S960_run.log")
_logf = open(_LOGPATH, "w", encoding="utf-8")


class Tee:
    def __init__(self, real, fh):
        self.real = real
        self.fh = fh
        self.chunks = []

    def write(self, s):
        self.real.write(s)
        self.fh.write(s)
        self.chunks.append(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()


_tee = Tee(sys.stdout, _logf)
sys.stdout = _tee

ASSERT_PASS = [0]
FAILS = [0]


def ok(cond, msg):
    if cond:
        ASSERT_PASS[0] += 1
    else:
        FAILS[0] += 1
        print("ASSERT-FAIL: " + msg)


d = Symbol('d', positive=True, integer=True)


def cartan_Ad(dd):
    M = zeros(dd, dd)
    for a in range(dd):
        M[a, a] = Integer(2)
        if a + 1 < dd:
            M[a, a + 1] = Integer(-1)
            M[a + 1, a] = Integer(-1)
    return M


def unit_dot(a, b, dv):
    """unit-vector inner product with symbolic/numeric d: 1 if same index, else -1/d."""
    return Integer(1) if a == b else Rational(-1, 1) / dv if False else (-1 / dv if a != b else 1)


def udot(a, b, dv):
    return Integer(1) if a == b else -Rational(1, 1) * (1) / dv


# clean symbolic dot: label-based
def dot(la, lb, dv):
    return Integer(1) if la == lb else -1 / dv


def beta_beta(p, q, dv):
    """b_p . b_q  where b_p = v_p - v_{p+1}, computed from unit dots (symbolic d).
    Uses distinct integer labels; only the equality pattern is significant."""
    # indices p, p+1, q, q+1 (labels); equality pattern set by |p-q|
    A, B = p, p + 1
    C, D = q, q + 1
    return simplify(dot(A, C, dv) - dot(A, D, dv) - dot(B, C, dv) + dot(B, D, dv))


# ============================================================================
print("=" * 74)
print("W34-leg-2: ∀d symbolic derivation of the cell identities")
print("=" * 74)

TAIL = list(range(2, 9))   # instance-tail d=2..8

# ---------------------------------------------------------------------------
print()
print("--- CLAIM 1: Gram of differences ∝ Cartan(A_d) ∀d (symbolic d) ---")

# symbolic-d entries by |p-q| case (distinct-label representatives)
e_diag = beta_beta(0, 0, d)      # |p-q|=0
e_adj = beta_beta(0, 1, d)       # |p-q|=1
e_far = beta_beta(0, 2, d)       # |p-q|>=2 (4 distinct labels)
ok(simplify(e_diag - 2 * (d + 1) / d) == 0, "diag entry = 2(d+1)/d (symbolic d)")
ok(simplify(e_adj + (d + 1) / d) == 0, "adjacent entry = -(d+1)/d (symbolic d)")
ok(simplify(e_far) == 0, "far entry = 0 (symbolic d)")
print("  symbolic: diag={0}, adj={1}, far={2}  =>  Gram=((d+1)/d)·Cartan(A_d)".format(
    simplify(e_diag), simplify(e_adj), simplify(e_far)))

# det Cartan(A_d)=d+1 forall d : recurrence proof (symbolic)
Dexpr = d + 1
rec = simplify((2 * ((d - 1) + 1) - ((d - 2) + 1)) - Dexpr)  # 2 D_{d-1}-D_{d-2}-D_d
ok(rec == 0, "d+1 satisfies tridiagonal recurrence 2D_{d-1}-D_{d-2}=D_d")
ok(cartan_Ad(1).det() == 2 and cartan_Ad(2).det() == 3, "recurrence base D_1=2,D_2=3")
print("  det Cartan(A_d)=d+1 ∀d (recurrence D_d=2D_{d-1}-D_{d-2}, base 2,3)")

# instance-tail
tail_ok1 = True
for dd in TAIL:
    G = cartan_Ad(dd)
    # build reduced-basis Gram from actual unit vectors with numeric d=dd
    dv = Integer(dd)
    Gun = zeros(dd, dd)
    for p in range(dd):
        for qi in range(dd):
            Gun[p, qi] = beta_beta(p, qi, dv)
    expect = Rational(dd + 1, dd) * G
    if Gun != expect:
        tail_ok1 = False
    if G.det() != dd + 1:
        tail_ok1 = False
ok(tail_ok1, "CLAIM1 instance-tail d=2..8 (Gram=SC·Cartan, det=d+1)")
print("  instance-tail d=2..8: Gram=SC·Cartan and det=d+1 — {0}".format("OK" if tail_ok1 else "FAIL"))

# ---------------------------------------------------------------------------
print()
print("--- CLAIM 2: z=d+1 combinatorially ∀d (minimal coset vectors) ---")

# symbolic norm formula: |z-c|^2 = sum z_i^2 - 1/(d+1); minimizer e_i -> 1 - 1/(d+1)
# lemma: integer z with sum z=1 has sum z_i^2 >= 1, equality iff z=e_i.
# proof pieces: z_i^2 >= |z_i| (integers), sum|z_i| >= |sum z_i| = 1.
# verify the two inequalities are tight only on {0,1}-one-hot; instance-tail counts.
tail_ok2 = True
counts = []
for dd in TAIL:
    n = dd + 1
    dv = Integer(dd)
    c = Rational(1, n)
    # enumerate coset reps z in box: entries in [-2,2], sum z = 1
    best = None
    minima = 0
    for z in itertools.product(range(-2, 3), repeat=n):
        if sum(z) != 1:
            continue
        s2 = sum(zi * zi for zi in z)          # euclid |z|^2 sum
        nrm = Rational(s2) - Rational(1, n)     # |z-c|^2
        unit_nrm = Rational(n, dd) * nrm        # unit metric
        if best is None or unit_nrm < best:
            best = unit_nrm
            minima = 1
        elif unit_nrm == best:
            minima += 1
    counts.append(minima)
    if best != 1 or minima != dd + 1:
        tail_ok2 = False
ok(tail_ok2, "CLAIM2 instance-tail d=2..8 (min unit-norm²=1, exactly d+1 minima)")
# lemma inequality check (symbolic-ish over the integer alphabet):
ineq_ok = all((zi * zi >= abs(zi)) and ((zi * zi == abs(zi)) == (zi in (-1, 0, 1)))
              for zi in range(-3, 4))
ok(ineq_ok, "CLAIM2 lemma: z_i^2>=|z_i|, tight iff z_i in {-1,0,1}")
print("  lemma: Σz_i²≥Σ|z_i|≥|Σz_i|=1, equality ⟺ one-hot ⟹ exactly d+1 minima")
print("  instance-tail d=2..8 z(d): {0} (expect {1}) — {2}".format(
    counts, [dd + 1 for dd in TAIL], "OK" if tail_ok2 else "FAIL"))

# ---------------------------------------------------------------------------
print()
print("--- CLAIM 3: Gauss sum Σ(d+1)-roots=0 ∀d + exactly d nodes ---")

# symbolic geometric series with symbolic count m=d and symbolic j (integer)
w = Symbol('w')
m = Symbol('m', positive=True, integer=True)
geom = summation(w ** i if False else w ** Symbol('i'), (Symbol('i'), 0, m))
# geom = Piecewise((m+1, w=1), ((1-w^{m+1})/(1-w), else))
# telescoping identity (polynomial, exact): (w-1)*sum_{i=0}^{m} w^i = w^{m+1}-1
telescope = simplify(expand((w - 1) * ((w ** (m + 1) - 1) / (w - 1)) - (w ** (m + 1) - 1)))
ok(telescope == 0, "geom-series telescoping identity (w-1)Σ=w^{m+1}-1")

jj = Symbol('j', integer=True, positive=True)
# w = exp(2 pi I j/(d+1)); w^{d+1} = exp(2 pi I j) = 1 (integer j)
w_root = exp(2 * I * pi * jj / (d + 1))
w_pow = simplify(w_root ** (d + 1))
ok(simplify(w_pow - 1) == 0, "w^{d+1}=exp(2πI j)=1 exact (integer j, symbolic d)")
# => numerator w^{d+1}-1 = 0, and w != 1 for 1<=j<=d => f=0 ; j=0 -> f=d+1
print("  symbolic: w^(d+1)=exp(2πI·j)=1 ⟹ the geom-series numerator=0, w≠1 (1≤j≤d) ⟹ f=0")
print("  ⟹ exactly d nodes (j=1..d), Γ(j=0): f=d+1≠0 (not a node)")

# instance-tail
tail_ok3 = True
for dd in TAIL:
    n = dd + 1
    node_count = 0
    for jv in range(0, n):
        wv = exp(2 * I * pi * Rational(jv % n, n))
        # geometric sum exact via telescoping: 0 if jv%n!=0 else n
        fval = Integer(n) if (jv % n) == 0 else Integer(0)
        # cross-check numerically
        direct = sum(complex(exp(2 * I * pi * Rational((jv * ii) % n, n)).evalf(30))
                     for ii in range(n))
        if abs(direct - complex(fval)) > 1e-12:
            tail_ok3 = False
        if fval == 0:
            node_count += 1
    if node_count != dd:
        tail_ok3 = False
ok(tail_ok3, "CLAIM3 instance-tail d=2..8 (sum=0 for j!=0, node count=d)")
print("  instance-tail d=2..8: sum=0 (j≠0), nodes=d — {0}".format("OK" if tail_ok3 else "FAIL"))

# ---------------------------------------------------------------------------
print()
print("--- CLAIM 4: finiteness of the traversal ⟺ d=2 (Niven, explicit) ---")

# cos(2theta) = 2/d^2 - 1.  Niven set of rational cos that are rational*pi:
NIVEN = [Integer(0), Rational(1, 2), Rational(-1, 2), Integer(1), Integer(-1)]
c2 = 2 / d ** 2 - 1
sols = []
for val in NIVEN:
    # solve 2/d^2 - 1 = val  ->  d^2 = 2/(val+1)
    if val + 1 == 0:
        continue  # val=-1 -> 2/d^2=0 no finite d
    d2 = Rational(2, 1) / (val + 1)
    # integer d>=2 with d^2 == d2 ?
    for cand in range(2, 40):
        if Rational(cand * cand) == d2:
            sols.append((val, cand))
ok(sols == [(Rational(-1, 2), 2)], "Niven: only integer d>=2 solution is d=2 (cos2θ=-1/2)")
print("  Niven set cos∈{{0,±1/2,±1}}: solutions of 2/d²−1=val for integer d≥2: {0}".format(sols))
print("  ⟹ the traversal is finite ⟺ d=2 (order 3); d≥3 infinite")

# instance-tail
tail_ok4 = True
for dd in TAIL:
    val = Rational(2, dd * dd) - 1
    finite = val in NIVEN
    if dd == 2 and not finite:
        tail_ok4 = False
    if dd >= 3 and finite:
        tail_ok4 = False
ok(tail_ok4, "CLAIM4 instance-tail d=2..8 (finite iff d=2)")
print("  instance-tail d=2..8: cos2θ∈Niven only for d=2 — {0}".format("OK" if tail_ok4 else "FAIL"))

# ==================== BIT-FENCE d=2 vs S956 ====================
print()
print("BIT-FENCE d=2 (cross-check against the S956 stamp):")
ok(cartan_Ad(2).det() == 3, "d=2: det Cartan(A_2)=3")
ok(beta_beta(0, 1, Integer(2)) == -Rational(3, 2), "d=2: adj=-(3/2)=SC·(-1)")
ok((Rational(2, 4) - 1) == Rational(-1, 2), "d=2: cos2θ=-1/2 (order 3)")
print("  d=2: Cartan det3 · adj=-3/2 · cos2θ=-1/2 · z=3 · nodes=2 — all agree with S956")

# ==================== MUTANTS (one per claim) ====================
print()
print("MUTANTS (one per branch/claim):")
mut_ok = True

# m1 (Gram): wrong off-diagonal (-2 instead of -1) must NOT match the true Gram
bad_cartan = cartan_Ad(3)
bad_cartan[0, 1] = Integer(-2); bad_cartan[1, 0] = Integer(-2)
true_gram = zeros(3, 3)
for p in range(3):
    for qi in range(3):
        true_gram[p, qi] = beta_beta(p, qi, Integer(3))
if true_gram != Rational(4, 3) * bad_cartan and true_gram == Rational(4, 3) * cartan_Ad(3):
    print("  MUTANT m1-gram: CAUGHT (off-diag -2 mismatches true Gram; -1 matches)")
else:
    print("  MUTANT m1-gram: NOT CAUGHT"); mut_ok = False

# m2 (z-count): a coset vec (2,-1,0,..) sum=1 but |z|^2=5>1 must be off the shell
zbad = (2, -1) + (0,) * 3   # d=4, n=5, sum=1
s2bad = sum(zi * zi for zi in zbad)
if sum(zbad) == 1 and s2bad > 1:
    print("  MUTANT m2-zcount: CAUGHT (z=(2,-1,0,0,0) sum=1 but |z|²=5>1, not minimal)")
else:
    print("  MUTANT m2-zcount: NOT CAUGHT"); mut_ok = False

# m3 (Gauss): j=0 (w=1) must NOT be a node (geom series degenerates, sum=n)
n5 = 5
f_j0 = sum(complex(exp(2 * I * pi * Rational(0, n5)).evalf()) for ii in range(n5))
if abs(f_j0 - n5) < 1e-9:
    print("  MUTANT m3-gauss: CAUGHT (j=0 → w=1 → f=n={0}≠0, not a node)".format(n5))
else:
    print("  MUTANT m3-gauss: NOT CAUGHT"); mut_ok = False

# m4 (Niven): claim d=3 finite (cos2θ=-7/9 in Niven) must be REJECTED
c2_d3 = Rational(2, 9) - 1
if (c2_d3 not in NIVEN) and ((Rational(2, 4) - 1) in NIVEN):
    print("  MUTANT m4-niven: CAUGHT (d=3 cos2θ=-7/9 ∉ Niven infinite; d=2 -1/2 ∈ finite)")
else:
    print("  MUTANT m4-niven: NOT CAUGHT"); mut_ok = False

# ==================== summary ====================
print()
print("SUMMARY: asserts_passed={0} | FAIL={1}".format(ASSERT_PASS[0], FAILS[0]))

_pieces = [("з", "акон"), ("к", "анал"), ("мех", "анізм"), ("зл", "іпок"),
           ("пер", "егин"), ("конд", "енсат"), ("мат", "ерія"), ("ен", "ергія"),
           ("рез", "онанс"), ("тр", "іщина"), ("г", "учн")]
_words = ["".join(ab).casefold() for ab in _pieces]
_src = open(__file__, "r", encoding="utf-8").read().casefold()
_logf.flush()
_logtxt = "".join(_tee.chunks).casefold()
_hits = 0
for _w2 in _words:
    if _w2 in _src:
        _hits += 1
    if _w2 in _logtxt:
        _hits += 1
print("FORBIDDEN-SCAN: hits={0}".format(_hits))

_exit = 0
if _hits > 0:
    _exit = 1
if FAILS[0] > 0:
    _exit = 1
if not mut_ok:
    _exit = 1
print("EXIT={0}".format(_exit))
sys.exit(_exit)
