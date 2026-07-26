# -*- coding: utf-8 -*-
# DIM: na (exact sympy linear algebra on so(p,q) blocks, n≤6; the court's P4 measurement — handles 0).
# COURT HAND W32 (per an internal project directive): P4 measurement + a liveness check of the gl-candidate
# on the kernel blocks J2/J4. Primitives are pulled by a SLICE from the source of hand S929
# (exec of the definitions segment) — the court measures ITS constructions, not a retelling.
import sys, os
sys.stdout.reconfigure(encoding='utf-8')

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "S929_2.py")
src = open(SRC, encoding="utf-8").read()

a = src.index("# ---------- exact linear-algebra primitives")
b = src.index("# ---------- shared per-stratum pipeline")
c = src.index("def q_gram(")
d = src.index('print("--- jordan')
seg = src[a:b] + "\n" + src[c:d]

ns = {}
exec("import itertools\n"
     "from math import lcm as _ilcm\n"
     "from sympy import Matrix, Integer, Rational, zeros, eye, diag, symbols, Poly, gcd, sqrt as ssqrt\n"
     "X_ = symbols('x')\n"
     "PASS=[0]; FAILS=[0]\n"
     "def ok(cond, msg):\n"
     "    if cond: PASS[0]+=1\n"
     "    else:\n"
     "        FAILS[0]+=1\n"
     "        print('ASSERT-FAIL: '+msg)\n"
     + seg, ns)

ok = ns["ok"]; PASS = ns["PASS"]; FAILS = ns["FAILS"]
make_eta=ns["make_eta"]; so_basis=ns["so_basis"]; block_gen=ns["block_gen"]
centralizer=ns["centralizer"]; blocks_of=ns["blocks_of"]; restrict=ns["restrict"]
span_basis=ns["span_basis"]; span_rank=ns["span_rank"]; same_span=ns["same_span"]
contains=ns["contains"]; stack_flats=ns["stack_flats"]; flat=ns["flat"]
try_gl_model=ns["try_gl_model"]; is_nilp=ns["is_nilp"]; is_so=ns["is_so"]
wedge=ns["wedge"]; unit=ns["unit"]; gram2=ns["gram2"]
q_gram=ns["q_gram"]; jordan_scan=ns["jordan_scan"]; d0_on=ns["d0_on"]; d1_on=ns["d1_on"]
derived_and_closed=ns["derived_and_closed"]
Matrix=ns["Matrix"]; zeros=ns["zeros"]; Integer=ns["Integer"]

def centralizer_in_span(F, T, d):
    # elements of span(F) with [X, T] = 0: a basis — exactly as in centralizer(), but within subspace F
    cols = [flat(B*T - T*B) for B in F]
    M = Matrix.hstack(*cols)
    out = []
    for v in M.nullspace():
        X = zeros(d, d)
        for k in range(len(F)):
            if v[k, 0] != 0:
                X = X + v[k, 0]*F[k]
        out.append(X)
    return out

# J-cases — the construction is VERBATIM per hand S929
JC = []
_eta = make_eta(2,2); _bas = so_basis(4,_eta)
_S = block_gen(4,_eta,0,2,1) + block_gen(4,_eta,1,3,1)
_cb = centralizer(_S,_bas)
_Xn = jordan_scan(_cb, q_gram(_cb), 4)
ok(_Xn is not None, "J1 scan")
JC.append(("J1",2,2,_S,_Xn))
_eta = make_eta(4,2); _S = block_gen(6,_eta,0,1,1)
_x,_y = d0_on([2,3,4,5],4,6); ok(gram2(_x,_y,_eta).is_zero_matrix,"J2 gram0")
JC.append(("J2",4,2,_S,wedge(_x,_y,_eta)))
_eta = make_eta(5,1); _S = block_gen(6,_eta,0,1,1)
_x,_y = d1_on([2,3,4,5],5,6); ok(gram2(_x,_y,_eta).rank()==1,"J3 gram1")
JC.append(("J3",5,1,_S,wedge(_x,_y,_eta)))
_eta = make_eta(3,3); _S = block_gen(6,_eta,0,3,1)
_x,_y = d0_on([1,2,4,5],3,6); ok(gram2(_x,_y,_eta).is_zero_matrix,"J4 gram0")
JC.append(("J4",3,3,_S,wedge(_x,_y,_eta)))

print("--- COURT P4: image of c(A) per block == centralizer of N|block inside factor(=image of c(S)|block) ---")
eq_cnt = 0; tot = 0
for (jid,p,q,S,N) in JC:
    n = p+q
    eta = make_eta(p,q); bas = so_basis(n,eta)
    ok(is_so(N,eta) and is_nilp(N) and not N.is_zero_matrix, jid+" N legal")
    ok((S*N-N*S).is_zero_matrix, jid+" [S,N]=0")
    A = S+N
    recs = blocks_of(S, eta, n, jid)
    cS = centralizer(S,bas); cA = centralizer(A,bas)
    for r in recs:
        d = r["d"]
        F  = span_basis([restrict(M,r,jid) for M in cS], d)   # factor: image of c(S)|block
        LA = span_basis([restrict(M,r,jid) for M in cA], d)   # image of c(A)|block
        Nb = restrict(N,r,jid)
        CF = centralizer_in_span(F, Nb, d)                    # centralizer of N|block IN the factor
        same = same_span(LA, CF, d)
        _dcc, closed = derived_and_closed(LA, d)
        ok(closed, jid+" image closed f="+r["pstr"])
        tot += 1
        if same: eq_cnt += 1
        print("{0} f={1} d={2} sig={3} | dim F={4} | dim cF(N)={5} | dim im c(A)={6} | span-EQUAL={7}".format(
            jid, r["pstr"], d, r["sig"], len(F), len(CF), len(LA), "YES" if same else "NO"))
print("COURT P4: {0}/{1} span-equal".format(eq_cnt, tot))

print("--- COURT: gl(2)-liveness on kernel (2,2) blocks of J2/J4 (dim-vs-bracket discriminant) ---")
for (jid,p,q,S,N) in JC:
    if jid not in ("J2","J4"):
        continue
    n = p+q
    eta = make_eta(p,q); bas = so_basis(n,eta)
    A = S+N
    recs = blocks_of(S, eta, n, jid+"g")
    cA = centralizer(A,bas)
    for r in recs:
        if r["pstr"] != "x" or r["d"] != 4:
            continue
        LA = span_basis([restrict(M,r,jid+"g") for M in cA], r["d"])
        glm = try_gl_model(r)
        live = glm is not None
        dim_eq = live and (span_rank(glm,4) == span_rank(LA,4))
        span_eq = live and same_span(LA, glm, 4)
        print("{0} kernel(2,2): gl(2)-model constructible={1} | dim gl=4 vs dim im={2} -> dim-equal={3} | span-equal={4}".format(
            jid, live, span_rank(LA,4), dim_eq, span_eq))
        ok(live, jid+" gl model must be constructible on split (2,2)")
        ok(not span_eq, jid+" image must NOT be gl(2) by brackets")

print("SUMMARY: asserts={0} FAIL={1}".format(PASS[0], FAILS[0]))
sys.exit(0 if FAILS[0]==0 else 1)
