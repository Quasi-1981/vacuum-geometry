# -*- coding: utf-8 -*-
# DIM: na (Beta's review of probe S929: re-measure P4 (equality of the image with the
#      canonical centralizer of the image of N in the factor) + count the DISCRIMINATING
#      power of the menu, block by block. Handles 0.)
import sys
sys.stdout.reconfigure(encoding='utf-8')

import os
import importlib.util

_HERE = os.path.dirname(os.path.abspath(__file__))
_LOGPATH = os.path.join(_HERE, "S930_run.log")
_logf = open(_LOGPATH, "w", encoding="utf-8")


class Tee:
    def __init__(self, real, fh):
        self.real = real
        self.fh = fh

    def write(self, s):
        self.real.write(s)
        self.fh.write(s)
        return len(s)

    def flush(self):
        self.real.flush()
        if not self.fh.closed:
            self.fh.flush()

    def reconfigure(self, **kw):
        return None


sys.stdout = Tee(sys.stdout, _logf)

PASS = [0]
FAIL = [0]


def ok(cond, msg):
    if cond:
        PASS[0] += 1
    else:
        FAIL[0] += 1
        print("ASSERT-FAIL: " + msg)


# --- import S929's hand WITHOUT running it: pull only the primitives, not the script ---
# (the hand runs everything on import; so we read the source and execute ONLY up to the marker)
_SRC = open(os.path.join(_HERE, "S929_2.py"),
            "r", encoding="utf-8").read()
_CUT = _SRC.index("# ---------- strata (semisimple S) ----------")
# + a clean slice of the Jordan-block defs (q_gram/jordan_scan/d1_on/d0_on) VERBATIM from
#   the hand, so the review measures ITS constructions, not my copy
_J0 = _SRC.index("def q_gram(cb):")
_J1 = _SRC.index('print("--- jordan (A = S + N) ---")')
_HEAD = _SRC[:_CUT] + "\n\n" + _SRC[_J0:_J1]
# strip the hand's file tee, so as not to overwrite its log
_HEAD = _HEAD.replace('_logf = open(_LOGPATH, "w", encoding="utf-8")',
                      'import io\n_logf = io.StringIO()')
_HEAD = _HEAD.replace('sys.stdout = _tee', 'pass')
_ns = {"__file__": os.path.join(_HERE, "S929_2.py"),
       "__name__": "s929_head"}
exec(compile(_HEAD, "S929_head", "exec"), _ns)

for _n in ("make_eta", "so_basis", "wedge", "unit", "block_gen", "centralizer",
           "blocks_of", "restrict", "span_basis", "span_rank", "same_span",
           "contains", "preserves", "try_gl_model", "try_u_model", "model_so_g", "try_sp_model",
           "d0_on", "d1_on", "gram2", "q_gram", "jordan_scan", "is_nilp", "flat",
           "stack_flats", "derived_and_closed", "sc_table", "bt_verify", "is_so"):
    globals()[_n] = _ns[_n]

from sympy import Matrix, Integer, zeros, eye

ok(True, "primitives imported from S929 head without running its body")


def centralizer_in_model(Nr, model, d):
    """canonical centralizer of the element Nr INSIDE the subalgebra model ⊂ gl(d)"""
    mb = span_basis(model, d)
    if not mb:
        return []
    cols = [flat(B * Nr - Nr * B) for B in mb]
    ns = Matrix.hstack(*cols).nullspace()
    out = []
    for v in ns:
        M = zeros(d, d)
        for k in range(len(mb)):
            if v[k, 0] != 0:
                M = M + v[k, 0] * mb[k]
        out.append(M)
    return out


# ---------------- P4: Jordan rows ----------------
print("--- P4: image vs canonical centralizer of N-image inside the factor ---")

JC = []

# J1: (2,2), S = two equal B's, N from the scan (as in S929)
p, q = 2, 2
n = 4
eta = make_eta(n - q, q) if False else make_eta(p, q)
S = block_gen(n, eta, 0, 2, 1) + block_gen(n, eta, 1, 3, 1)
cb = centralizer(S, so_basis(n, eta))
N = jordan_scan(cb, q_gram(cb), n)
JC.append(("J1", p, q, S, N))

# J2: (4,2), S = R+(0,1), N = a D0-wedge on the complement
p, q = 4, 2
n = 6
eta = make_eta(p, q)
S = block_gen(n, eta, 0, 1, 1)
x, y = d0_on([2, 3, 4, 5], p, n)
JC.append(("J2", p, q, S, wedge(x, y, eta)))

# J3: (5,1)
p, q = 5, 1
n = 6
eta = make_eta(p, q)
S = block_gen(n, eta, 0, 1, 1)
x, y = d1_on([2, 3, 4, 5], p, n)
JC.append(("J3", p, q, S, wedge(x, y, eta)))

# J4: (3,3)
p, q = 3, 3
n = 6
eta = make_eta(p, q)
S = block_gen(n, eta, 0, p, 1)
x, y = d0_on([1, 2, 4, 5], p, n)
JC.append(("J4", p, q, S, wedge(x, y, eta)))

P4_ROWS = [0]
P4_EQ = [0]
P4_NE = [0]

for (jid, p, q, S, N) in JC:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    A = S + N
    ok(is_so(N, eta) and is_nilp(N) and not N.is_zero_matrix, jid + " N is nonzero nilpotent in so")
    ok((S * N - N * S).is_zero_matrix, jid + " [S,N]=0")
    recs = blocks_of(S, eta, n, jid)
    cbA = centralizer(A, bas)
    for r in recs:
        d = r["d"]
        L = span_basis([restrict(M, r, jid) for M in cbA], d)
        Nr = restrict(N, r, jid)
        # the block's factor by the MENU (the same pick as in the hand)
        fac = None
        fname = None
        glm = try_gl_model(r)
        um = try_u_model(r)
        som = model_so_g(r["g"])
        # the factor = the menu item that equals the image of c(S) (not c(A)) on the block
        cbS = centralizer(S, bas)
        LS = span_basis([restrict(M, r, jid) for M in cbS], d)
        for (nm, mo) in (("gl", glm), ("u", um), ("so", som)):
            if mo is not None and same_span(LS, mo, d):
                fac, fname = mo, nm
                break
        ok(fac is not None, jid + " factor of c(S) identified on block " + r["pstr"])
        if fac is None:
            continue
        cN = centralizer_in_model(Nr, fac, d)
        eq = same_span(L, cN, d) if (L or cN) else True
        P4_ROWS[0] += 1
        if eq:
            P4_EQ[0] += 1
        else:
            P4_NE[0] += 1
        print("{0} | block f={1} d={2} sig={3} | factor={4} dim={5} | N|block {6} | dim c(A)|block={7} | dim c_factor(N)={8} | EQUAL={9}".format(
            jid, r["pstr"], d, r["sig"], fname, span_rank(fac, d),
            "ZERO" if Nr.is_zero_matrix else "nonzero",
            span_rank(L, d), span_rank(cN, d), "YES" if eq else "NO"))
        ok(eq, jid + " P4 equality on block " + r["pstr"])

print("P4-SUMMARY: rows={0} | EQUAL={1} | NOT-EQUAL={2}".format(P4_ROWS[0], P4_EQ[0], P4_NE[0]))

# ---------------- MENU DISCRIMINATING POWER ----------------
# For each block: how many menu items were even ALIVE (constructible) at all?
# menu-discriminating = >=2 live candidates, of which one is picked.
print("--- menu power per block (semisimple strata) ---")


def strata_local(p, q):
    out = []
    if p >= 2:
        out.append(("S1", [(0, 1, 1)]))
    if q >= 1:
        out.append(("S2", [(0, p, 1)]))
    if p >= 4:
        out.append(("S3", [(0, 1, 1), (2, 3, 2)]))
        out.append(("S4", [(0, 1, 1), (2, 3, 1)]))
    if p >= 2 and q >= 2:
        out.append(("S5", [(0, p, 1), (1, p + 1, 2)]))
        out.append(("S6", [(0, p, 1), (1, p + 1, 1)]))
    if p >= 3 and q >= 1:
        out.append(("S7", [(0, 1, 1), (2, p, 1)]))
    return out


TOT = [0]
DISCR = [0]
TRIV = [0]
AMBIG = [0]

for nn in (4, 5, 6):
    for q in range(0, nn // 2 + 1):
        p = nn - q
        if p < q:
            continue
        eta = make_eta(p, q)
        bas = so_basis(p + q, eta)
        for (sid, spec) in strata_local(p, q):
            S = zeros(p + q, p + q)
            for (i, j, par) in spec:
                S = S + block_gen(p + q, eta, i, j, par)
            recs = blocks_of(S, eta, p + q, sid)
            cb = centralizer(S, bas)
            for r in recs:
                d = r["d"]
                L = span_basis([restrict(M, r, sid) for M in cb], d)
                live = []
                glm = try_gl_model(r)
                if glm is not None:
                    live.append(("gl", span_rank(glm, d), glm))
                um = try_u_model(r)
                if um is not None:
                    live.append(("u", span_rank(um, d), um))
                som = model_so_g(r["g"])
                live.append(("so", span_rank(som, d), som))
                spm = try_sp_model(r, L)
                if spm is not None:
                    live.append(("sp", span_rank(spm, d), spm))
                matched = [nm for (nm, _dm, mo) in live if same_span(L, mo, d)]
                # are the live candidates DIFFERENT as SUBSPACES? (gl(1)=so(1,1), u(1)=so(2,0) coincide)
                distinct = []
                for (nm, dm, mo) in live:
                    dup = False
                    for (nm2, _d2, mo2) in distinct:
                        if same_span(mo, mo2, d):
                            dup = True
                            break
                    if not dup:
                        distinct.append((nm, dm, mo))
                TOT[0] += 1
                if len(distinct) >= 2:
                    DISCR[0] += 1
                    cls = "DISCRIMINATING"
                else:
                    TRIV[0] += 1
                    cls = "menu-of-1"
                if len(matched) >= 2:
                    AMBIG[0] += 1
                print("({0},{1}) {2} f={3} d={4} sig={5} | live={6} | distinct-as-subspaces={7} | matched={8} | {9}".format(
                    p, q, sid, r["pstr"], d, r["sig"],
                    ",".join(nm + ":" + str(dm) for (nm, dm, _m) in live),
                    ",".join(nm for (nm, _dm, _m) in distinct),
                    ",".join(matched) if matched else "-none-", cls))

print("POWER-SUMMARY: blocks={0} | >=2-distinct-live={1} | menu-of-1={2} | rows-with-ambiguous-match={3}".format(
    TOT[0], DISCR[0], TRIV[0], AMBIG[0]))

# ---------------- REVIEW ADD-ON: multiplicity m=3 (the hand's strata do not produce it) ----------------
# In the hand max m = 2 ⟹ Σm² is checked only at 1 and 4. m=3 must give 9. n<=6 — within the declared scope.
print("--- visa add-on: multiplicity m=3 (absent from the hand's strata; Sum m^2 must give 9) ---")

M3 = []
# (6,0): three identical R+'s ⟹ one block x^2+1, d=6, sig (6,0) ⟹ predicts u(3), dim 9
M3.append(("M3a", 6, 0, [(0, 1, 1), (2, 3, 1), (4, 5, 1)], "u", 9))
# (3,3): three identical B's ⟹ one block x^2-1, d=6, sig (3,3) ⟹ predicts gl(3,R), dim 9
M3.append(("M3b", 3, 3, [(0, 3, 1), (1, 4, 1), (2, 5, 1)], "gl", 9))
# (6,0): two identical + one different ⟹ u(2) + u(1) ⟹ dim c = 4 + 1 = 5
M3.append(("M3c", 6, 0, [(0, 1, 1), (2, 3, 1), (4, 5, 2)], "mixed", 5))

for (mid, p, q, spec, want, wantdim) in M3:
    n = p + q
    eta = make_eta(p, q)
    bas = so_basis(n, eta)
    S = zeros(n, n)
    for (i, j, par) in spec:
        S = S + block_gen(n, eta, i, j, par)
    recs = blocks_of(S, eta, n, mid)
    cb = centralizer(S, bas)
    viol = sum(1 for M in cb for r in recs if not preserves(M, r["P"]))
    ok(viol == 0, mid + " P1: all c elements preserve all blocks")
    outs = []
    tot = 0
    for r in recs:
        d = r["d"]
        L = span_basis([restrict(M, r, mid) for M in cb], d)
        tot += len(L)
        live = []
        glm = try_gl_model(r)
        if glm is not None:
            live.append(("gl({0},R)".format(d // 2), glm))
        um = try_u_model(r)
        if um is not None:
            live.append(("u({0})".format(d // 2), um))
        som = model_so_g(r["g"])
        live.append(("so{0}".format(r["sig"]), som))
        spm = try_sp_model(r, L)
        if spm is not None:
            live.append(("sp({0},R)".format(d), spm))
        matched = []
        for (nm, mo) in live:
            if same_span(L, mo, d):
                mob = span_basis(mo, d)
                ok(bt_verify(L, mob, sc_table(mob, d), d), mid + " bracket table vs " + nm)
                matched.append(nm)
        outs.append("f={0} d={1} sig={2} live=[{3}] dim(image)={4} matched={5}".format(
            r["pstr"], d, r["sig"], ",".join(nm + ":" + str(span_rank(mo, d)) for (nm, mo) in live),
            span_rank(L, d), ",".join(matched) if matched else "-NONE(OTHER/sub-of)-"))
    ok(len(cb) == tot, mid + " P1: dim c = sum of block images (no cross-terms)")
    ok(len(cb) == wantdim, mid + " P3: dim c = {0} (predicted Sum m^2 arithmetic)".format(wantdim))
    print("{0} ({1},{2}) | dimc={3} (want {4}) | {5}".format(mid, p, q, len(cb), wantdim, " ; ".join(outs)))
print("VISA-SUMMARY: asserts_passed={0} | FAIL={1}".format(PASS[0], FAIL[0]))
print("EXIT={0}".format(1 if FAIL[0] > 0 else 0))
sys.stdout = sys.stdout.real
_logf.close()
sys.exit(1 if FAIL[0] > 0 else 0)
