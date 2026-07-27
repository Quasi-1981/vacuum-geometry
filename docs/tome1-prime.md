# TOME I — PRIME FLOOR: pure mathematics

> **Status of this text.** This English text is the original, written directly from the project's
> canonical sources. Qualifiers are rendered conservatively and close to the letter: a qualifier is
> exactly what gets lost in transfer.
>
> **Scope (by the author's word, 2026-07-25).** This is the first tome published. Tomes II and III —
> the bridge, $(3,1)$, and the identifications — are not yet ready; they will follow. Tome I therefore
> has to stand alone, so §7 lists everything it does not prove.
>
> **What it is.** One chain, already proved, each link carrying the address of its primary act.
> Distillation moves no status: where the source says `candidate`, `candidate` stands here. This is
> representation, not derivation — no new claim and no new number.
> **Floor fence:** a word is admissible only if its referent is derived or axiomatic on this floor.

## Why any of this

<!--program-preface-->
*This preface speaks of the programme one level above this tome's floor, so it may use words the
floor itself has not earned; the floor fence applies from §0 onward.*

The original goal was rather modest: to determine whether a topology could be constructed at all if
**E** (energy) were taken as the fundamental quantity and **Λ** as the characteristic unit of vacuum
stiffness. The answer proved to be yes.

The difficulty emerged elsewhere. Modern physics offers a wide variety of theories, each describing
only its own limiting regime—a diverse collection of mathematical frameworks for projecting
interactions and processes toward the sub-Planck scale. A broader analysis suggested that these
families of topologies are structurally overdetermined. Together with the heterogeneity of the
existing approaches, this observation motivated a step one level deeper, leading to the present
volume.

**Tome I** is therefore a purely mathematical abstraction. Its purpose is to construct a more robust
mathematical framework for these topologies and to place them within a unified space of solutions
based on Lie algebras.

**A note on confidence.** The author does not regard personal intuition as sufficient
evidence—including the possibility that the entire construction may be fundamentally mistaken. This
is not a rhetorical disclaimer but a design principle. Nothing in this work is intended to rest on
the author's confidence alone. Every claim is accompanied by its current status
(measured / ∀d-symbolic / candidate), the probe that produced it, and the conditions under which it
would be considered falsified.

The reader is therefore invited to approach this work exactly as it was written: assume the author
may be wrong, and verify every step.

*How the work was organized — the roles, the ids, and what is public — is described on
[Process & roles](process.md).*
<!--/program-preface-->

---

## 0 · READ THIS FIRST — the boundaries

Six boundaries, stated before the body. An unwritten boundary works as a cancelled negative.

1. **Six status words, and they are not equal.** The graph's `status` field takes exactly six
   values: `axiom` (declared premise) · `⚓ measured` (probe plus independent visa) · `candidate`
   (a form fitted to data — neither a derivation nor a test; rendered `candidate/cast` in this tome
   where the fit is a constant) · `measured-negative` (proved absence) · `named` (identified, not
   yet derived) · `author-word` (the author's own addition, carried as data). Most of floor −2,
   T26–T39, is `candidate`. **`∀d-symbolic` is not a status** — it is a qualifier in verdict prose
   saying a `measured` node was derived for symbolic $d$ rather than enumerated; it never stands
   alone in the status column. Calling a row a theorem does not raise its grade; the status column
   does.
2. **Enumeration is not ∀d.** Symbolic derivations say so (`∀d-symb`). Enumerations state their
   range: $d\in\{2,3\}$, $n\le 6$, $n\le 7$, up to $n=10$. A $d=2$ instance is never a witness on its
   own, since the discriminants vanish identically there; the norm is symbolic $d$, or the pair
   $d\in\{2,3\}$.
3. **Classical frames count for nothing.** Root lattices $A_d$, Weyl alcoves, Pontryagin duality,
   Gauss sums, Schur's lemma, Levi–Malcev, İnönü–Wigner, Niven — all 📖 reference, multiplicity 0.
   They add no weight to anything below.
4. **Witnesses count by ancestor, not by citation.** Two descendants of one ancestor agreeing is
   multiplicity 1.
5. **One handle.** The whole tome declares exactly one dimensional input, $\Lambda$. Nothing below
   is a new constant; where a number looks new, its ancestor is named.
6. **Absent by construction.** Dynamics and the variational principle do not live on this floor
   (fence of floor −1: −1 is the factory of dials). The seam to the action formula, the two-component operator,
   the spectrum and $(3,1)$-non-splitting are Tome II. Interpretation is Tome III.

---

## 1 · INPUTS

### 1.1 The handle

$\Lambda$ is the ruler: a unit of length, and the only dimensional handle. It enters the central
coefficient exactly once. The debt "the slot $|\lambda| = 1/R^2$ is unmeasured" (S922) closes by
identification — the weight $\mu^2/\nu$ *is* the ruler, not a handle.
`AX-lambda` · **T22** · S983 · **J-0442**.

> **Boundary:** in the prime graph $\Lambda$ is `AX-lambda`; in the conveyor graph the same quantity
> is the node `AX-001`, and no edge joins them (hole **E-1**). Tome I uses the prime branch only.
> Tome II writes the seam out; unifying the graphs is a separate pipeline.

### 1.2 Four declared premises (inputs, not numbers)

| # | premise | content | address |
|:-:|:--|:--|:--|
| 1 | **`AX-alphabet`** | a bare set of $d+1$ elements, with democracy $S_{d+1}$ (the full symmetric group) | §0 · arc S1023–S1026 · **J-0476** |
| 2 | **`AX-dimer`** | the dimer marker $u_0$ on a bond | §0 · S597 |
| 3 | **`AX-closure`** | the closure **method**: construction of the order $\prec$ on real quadratic spaces. A method input, not a numerical handle | §3.arc-0 · preprint-1 · S899+ |
| 4 | **`AX-lambda`** | §1.1 | — |

Four roots, no more. Two former inputs are now theorems:

- **`AX-cell`, the cell $A_d$: input → theorem**, `deps = [A-space]`. Derived from the alphabet at
  the C2 swap by three independent routes — rank, character, Molien — with circularity killed.
  **J-0476 VALID ∀d.** Geometry has left the inputs of the programme.
- **`AX-indef`, the indefiniteness break: input → theorem**, `deps = [AX-alphabet, AX-dimer]`, at
  0 handles. It was the postulate "at least one break axis exists". What carries the demotion: it
  consumes exactly the prohibition $q=0$; the monotone laws are robust; the bound $n\le9$ drops.
  **Permanent qualifier: at 0 handles, sign is the only channel.**
  §0 · S1066 · **J-0515…J-0518** · graph operation 5→4 **S1093/J-0527**.

> **Boundary:** 5→4 counts roots. It strengthens no statement, and no status rose with it.

---

## 2 · DEMOCRACY → ARENA

Layer `alg` takes a simple algebra from the alphabet of §1.2-1.

| link | statement | status | address |
|:--|:--|:--|:--|
| **`OBJ-sln`** | $\mathfrak{sl}(n)$, $n=d+1$, simple, with a unique invariant form $\kappa$ (Schur) | ⚓ measured | §3.arc · S1025 · **J-0478** |
| **`A-space`** | the arena is a Cartan subalgebra; Gram $\propto \mathrm{Cartan}(A_d)$, **∀d** | ⚓ measured | **J-0476** |
| **`A-kappa`** | arena and axis split by one mechanism, Jordan–Chevalley | ⚓ measured | **J-0478** |
| **`A-axis`** | signature $(d,1)$; axis $=\mathfrak{so}(2)_\alpha$; the marker selects a root | ⚓ measured | **J-0479** |
| **`A-time-neg`** | **⚓✗ measured negative:** the axis does not follow from the order. The Weyl channel yields the arrow (the $\mathbb{Z}/2$ orientation bit; no temporal content is claimed) [T28] and nothing more | ⚓✗ negative | arc S1023–S1026 |

> **Boundary (carved homonym):** the minus of the arc, geometric $\mathfrak{so}(2)_\alpha$, and the
> minus of T32, the arithmetic centre dual, share the arena and nothing else. There is no bridge.
> The edge `A-axis ⇢ AX-indef` is blocked for that reason [S1027]. Do not read them as one.

---

## 3 · THE TERMINAL STEP OF THE DERIVED SERIES

Take the derived series of the democracy $S_{d+1}$. **A terminal abelian step exists iff the series
reaches $\{e\}$ — that is, iff $d+1 \le 4$.**

| $d$ | derived series of $S_{d+1}$ | terminal step |
|:-:|:--|:--|
| 2 | $S_3 \rhd A_3 \rhd \{e\}$ | $\mathbb{Z}/3$ |
| 3 | $S_4 \rhd A_4 \rhd V_4 \rhd \{e\}$ | $V_4$ |
| 4 | $S_5 \rhd A_5 \rhd A_5 \rhd \dots$ — stalls | none |

Two bounds $d+1\le4$ were measured independently: the node marker [S1059B] and the four-component
reading/clock [S1063-T3]. One object realises both — this step.

- Status: **candidate**. The common ancestor is carved as a candidate, not as a theorem.
- Multiplicity: T1+T2 is **one** ancestor, multiplicity 1. The Galois rhyme is reference,
  multiplicity 0.
- Statuses S1059B and S1063-T3 are **not** moved; raising their weight through the ancestor takes a
  separate decision.
- Address: S1101 → **S1102** · visa **J-0531** · an independent re-run, 23/0/0 · seal
  `LADDER_ANCESTOR_EXANTE.md`.

> **Boundary (essential):** the verdict rests on reading the derived series by hand. Identifying the
> Galois group is an instrument inside the declared scope, not the support. The series is computed
> from our own democracy; taking a Galois form as input would have killed the probe.

**The writing rule** (same seal; installed as a point of enforcement, below). A legitimate number of
the canon is a citizen of a solvable tower over $\mathbb{Q}$. A tower with an unsolvable ($A_5$)
factor is a writing artefact. The register passed whole; $x^5-x-1$ and $x^5-x+1$ were cut.
**Qualifiers:** (a) the rule measures algebraic numbers — minimal polynomial over $\mathbb{Q}$;
$\pi$ and the transcendentals fall outside it; (b) the verdict rests on reading the series by hand.

> **Installation (P4, S1136):** the rule is no longer a paragraph. Every published number carries
> obligatory fields — value, address, unit, type, algebraic form — and a checker with an exit code
> computes the Galois group of the minimal polynomial and rejects an unsolvable tower. An empty
> field is visible; a forgotten principle is not. The $x^5-x-1$ demonstration (group of order 120,
> rejected) lives as data, not as a hidden test.

---

## 4 · THE CELL $A_d$

### 4.1 The cell

$d+1$ unit axes, pairwise Gram $-1/d$; Gram $= SC\cdot\left(I - \dfrac{\mathbf{1}\mathbf{1}^\mathsf{T}}{d+1}\right)$ with
$SC = (d+1)/d$ (unit diagonal, off-diagonal $-1/d$). These are the weights of the fundamental representation of $\mathfrak{su}(d+1)$.
`OBJ-cell` · T26/T19 · S956.

### 4.2 T19 — the cell generates a lattice

Under translational closure the cell generates the root lattice $A_d$:

- Gram of differences $= \dfrac{d+1}{d}\,\mathrm{Cartan}(A_d)$, with $\det\mathrm{Cartan}(A_d) = d+1$ by the
  recursion $D_d = 2D_{d-1}-D_{d-2}$ (so $\det(\text{Gram of differences}) = SC^d\cdot(d+1)$, not $d+1$ itself);
- the two-sublattice construction gives a bipartite lattice with $z = d+1$, by the integer lemma
  $\sum z^2 \ge \sum|z| \ge |\sum z| = 1$ with equality iff one-hot; bond angles are $\arccos(-1/d)$;
- the nearest-neighbour operator has exactly $d$ zero modes at the symmetric points, the phases
  being the full $(d+1)$-th roots of unity; the Gauss sum telescopes,
  $(w-1)\sum w^i = w^{d+1}-1 = 0$;
- finiteness of the rotational circuit is exclusive to $d=2$ — by Niven's theorem, with
  $\cos\theta = -1/d$, finiteness requires $\cos 2\theta \in \{0,\pm1/2,\pm1\}$, which among
  integers $d\ge2$ holds only at $d=2$ (giving $\cos2\theta=-1/2$).

**Status: ⚓ measured; mechanism derived ∀d symbolically** [ver:2, leg 2].
Address: S956 + **S960** · **J-0425** / **J-0430**.

> **Boundaries (essential, verbatim from the verdict):**
> (a) the 📖 frames — root-lattice theory, $\det\mathrm{Cartan} = d+1$, the Gauss sum — carry
> multiplicity 0;
> (b) the circuit was measured under one operationalisation: two edge reflections at a vertex;
> (c) read the $d=2$ exclusivity with care. A planarity criterion does not select whether the lattice
> exists — lattices for $d\ge3$ exist translationally, and that is measured;
> (d) the zero-mode sector is independent of the circuit: the $d$ zero modes exist ∀d, including
> under an infinite circuit at $d\ge3$ *(the source words this line in the register of another floor;
> here it is restated on this floor — same content, different word)*;
> (e) "hexagonal trivalent $d=3$" is demoted: native $d=3$ has $z=4$, and trivalence does not
> survive ∀d.

### 4.3 The zero set

$f(k)=0$ imposes two real conditions, $\mathrm{Re}\,f=0$ and $\mathrm{Im}\,f=0$, so the zero set has
codimension 2 generically: points at $d=2$, a line at $d=3$, a surface at $d=4$.
$d=2,3$ are exact; $d\ge4$ is extrapolation.
S1045 · **J-0498** → correction S1047 · **J-0500**.

> **Boundary (a correction, not cosmetics):** the witness for the line at $d=3$ is the generic node,
> rank 2 and nullity 1. The special point $(0,\pi/2,\pi/2)$ is degenerate — rank 1, nullity 2 — and
> is not a witness; J-0498 cited it in error. More degenerate points exist, so codimension 2 holds
> generically, not everywhere.

### 4.4 Torus and alcoves

The Gauss zeros are the non-trivial characters of the centre $\mathbb{Z}/(d+1)$ (∀d symbolically),
and they sit exactly at the barycentres of the Weyl alcoves, with no tuning. The torus cuts into
$(d+1)!$ congruent alcoves. Cell and alcove are dual, and not congruent at any measured $d$: at
$d=2$ the alcove is similar to the cell (ratio $9/2$) but still not congruent; at $d\ge3$ the alcove
is cyclically symmetric but not even similar — a pre-registered expectation died here, an honest kill.
T26-(7) · S1002 · **J-0457**. Status: **candidate/cast**.

### 4.5 Two ladders

| ladder | content | status | address |
|:--|:--|:--|:--|
| Coxeter | eigenvalues $\mathrm{Cartan}(A_d) = 2-2\cos(\pi j/h)$, azimuth $2\pi/h$ | candidate | T26 · S1004 |
| commutant | $\dim$ of the stabiliser's invariant symmetric operators (ladder over $m$) $= 1/2/2/1$ ($d=2$) $\cdot\ 1/2/3/2$ ($d=3$) $\cdot\ 1/2/3/3$ ($d\ge4$, Lean debt above); block ranks $(1, m-1, d-m)$ for $1\le m\le d$ | candidate/cast | T26.3 · S1000-T2 · **Lean S1090/J-0526** |

> **Boundary on the Lean stamp:** machine checking hardens the verification; it does not raise the
> status, and a cast stays a cast. Scope $d\in\{2,3\}$ — $d\ge4$ is a debt, Mathlib having no
> `det_fin_four`. `ladderXX` is an exercise by ruling, and isotypy is imported from the ancestor
> S1000-T2.

---

## 5 · COLUMN AND CLOCK

### 5.1 Column, period, dual

The dimer marker $u_0$ of §1.2-2 generates the column. Translating the column is an exact symmetry
of period $P = d+1 = |\mathrm{disc}(A_d)|$ — the naive $P=2$ fails, since $2u_0$ is not in the
lattice. The dual is a compact circle $2\pi/(d+1)$ (Pontryagin), and the $\mathfrak{sl}$ axis runs
parallel to the column, at angle 0. One period holds $1A + 1B + (d-1)$ holes, with $\mathbb{Z}_2$
even. **0 new handles.**
`OBJ-column`, `OBJ-dual` · T26.5 · S1001 · **J-0456**. Status: **candidate**.

### 5.2 T34 — the cell holds one clock

All $d+1$ bonds of the cell fall in a single non-zero class of the centre $\mathbb{Z}/(d+1)$: the
pairwise differences $u_i-u_j = e_i-e_j$ lie in the root lattice exactly, computed in `Fraction`,
while the bond itself is non-integral. Hence the centre has rank 1, hence there is no second
independent column dual — **the cell admits exactly one clock**.

Three constructions of two columns at equal weights: independent gives $q_{\text{eff}}=2$;
locking $\nu_1=\nu_2$ still gives $q_{\text{eff}}=2$, so the lock does not save it — both are
counterfactual, the centre forbidding independence; the native merged construction, where the merge
is forced by the centre rather than chosen, gives $q_{\text{eff}}=1$ and $\mathrm{split}=0$.
So the democratic point of the native cell is not bistable.

T34 · S1013 · **J-0469** · **Lean S1088/J-0525** (10 theorems ∀n symbolically, not a cast).
Status: **candidate/cast** — the Lean stamp does not promote.

> **Boundaries:** scope is the cell $A_d$, whose centre is cyclic of rank 1; non-cyclic
> discriminants (D type, $\mathbb{Z}_2\times\mathbb{Z}_2$) are a named contrast and are not
> measured. Capacity is a structural count of independent duals, not a budgetary or dynamical
> quantity — that reading sits behind the fence. One freedom remains: which bond carries the marker.
> The dual is the same either way. A subtlety worth keeping: the locked and the merged constructions
> have identical spectra and differ only in mode activity, so the criterion for one clock is
> $\tilde{C}$, not the spectrum.

### 5.3 T37 — the clock is global

Two marked columns in adjacent cells of one $A_d$ lattice — different bonds, different cells, the
same-cell tautology killed by a mutant — share the same centre $\mathbb{Z}/(d+1)$. The common charge
has rank 1, so the pair has capacity 1: one connected lattice is one centre and **one clock
globally**. Separated domains are impossible in the native machinery, and synchronisation is
structural — geometry, not dynamics. Forced disagreement costs
$q_{\text{eff}}=2$ and $\mathrm{split} = 96\ (d=2)\ /\ 264\ (d=3)$.
T37 · S1018 · **J-0473**. Status: **candidate/cast**.

> **Boundary (essential):** at capacity 1, disagreement is not a native state. The cost prices a
> counterfactual rupture of the shared centre; it is not a quantity of the native state. Scope is a connected lattice —
> disconnected components are separate centres, hence separate bits, naturally.

---

## 6 · SIGN AND SIGNATURE

### 6.1 The order of closure (T1–T11)

`AX-closure` builds the order $\prec$ on real quadratic spaces.

| # | statement | address |
|:--|:--|:--|
| **T1** | $C_{2D} \prec C_{3D}$: $\mathfrak{so}(3)$ has no two-dimensional subalgebras, so the minimal carrier of closure is the definite $(3,0)$ | S899+S908 |
| T2 | the J sector is bracket-closed in every signature; a full K set forces the full $\mathfrak{so}(p,q)$ | S900 |
| T3 | $\prec$ is antisymmetric by construction — closure is idempotent, and there are no cycles | S900-ruling |
| T4 | the coordinate terminal is finite: $2\to14$ classes per rung, with $(3,3)$ passed in full, 32767 | S901 |
| **T5** | $\{$cell $\Rightarrow p\ge3\} \wedge \{$H-base $\Rightarrow q\ge1\} \Rightarrow n\ge4$, and at the minimum $(3,1)$ is the unique solution; $(2,2)$ falls by its own arithmetic | ver-note-3 + S908/F4 |
| T6 | the terminal of the order is the stabilisers of diagonal breaks, class by class 14/14 on $(3,3)$, up to the trivial class | S901×S902 |
| **T7** | $A^2=+\mathbb{1} \iff p=q$; $A^2=-\mathbb{1} \iff p,q$ even. Hence $(3,1)$ is rigid: it admits none | S903 |
| T8 | three generics: configuration→full algebra, symmetric break→zero symmetry, centraliser→Cartan | S901/S902/S903 |
| T9 | the mirror log pair $\ln\det\Omega - \ln\det g = \ln\det A$ with $\det A = \mathrm{Pf}^2/\det g$, so the sign law is non-strict: $(-1)^q \det A \ge 0$, with equality on the $\mathrm{Pf}=0$ stratum | S905 |
| T10 | $\mathrm{Pf}(S^\mathsf{T}\Omega S) = \det(S)\,\mathrm{Pf}(\Omega)$; the flip happens exactly under $\det=-1$ isometries | S906 |
| T11 | $\prec$ is equivariant under $\{\mathrm{Ad}$, sign swap, Cartan involution$\}$, 290/290 | S910+S912 |

All under ruling, visa **J-0410** (core VALID) with S912. Status: **⚓ measured**.

> **Boundaries — the qualifiers that die in transfer, hence verbatim:**
> **T1** rules on the **definite** sector. $\mathfrak{so}(2,1)$ does have two-dimensional non-abelian
> subalgebras (Borel, $[X,Y]=Y$); what carries the result is definiteness of the metric, not
> $\mathfrak{so}(3)$ as such.
> **T5** rests on a declared premise for the cell, plus an accepted definition, plus a $d\ge3$
> import. Only the implications are derived — embedding iff $d\le p$, by Gram and Sylvester — and the
> arithmetic $n\ge4$. So the multiplicity of T5 equals the multiplicity of its premises, and T5 is
> not an independent witness.
> **T6** without its qualifier is drift from act to register. It reads: class by class, up to the
> trivial class, the stabiliser minus terminal being 1.

### 6.2 The signature is an output

Given the cell ($p\ge3$) and the indefiniteness break ($q\ge1$, `AX-indef`), $(3,1)$ is forced —
by **T5** through minimality and uniqueness, and by **T7** through rigidity. In the graph `AX-sig`
descends from [T5, T7]; it is not a root. The composition — **3** definite arena
directions against **1** break axis — and the dimension $n=4$ are both derived. *(The source prints
this as a "3:1 ratio"; a bare "3:1" is forbidden by the homonym register, §8.1, so the composition
is written out by sides.)*
§0 · §3.arc-0 · a finding of the author, 2026-07-19, tightened by a project ruling.

> **Boundary:** the break axis does not come from closure. The minimal carrier of closure is the
> definite $(3,0)$ [T1]. It comes from `AX-indef` (§1.2), itself a theorem at 0 handles.
>
> There is no circularity, and the check is named [ruling, addendum to J-0549]. Closure does not
> presuppose $(3,1)$: T2 closes brackets in every signature, T4 passes $(3,3)$ in full, T7 states its
> condition in general $(p,q)$, and T5 selects afterwards rather than before. No link sees $(3,1)$
> before minimality picks it.
>
> What does remain is a different question. The premise $p\ge3$ carries a declared $d\ge3$ import
> (qualifier F4, visa **J-0410**). Declared, not hidden — so it is a named front, the $d{=}3$
> disposition of §7, and not circularity. The difference is essential: a hidden assumption devalues
> a derivation, a declared premise does not.

### 6.3 T32 — the first minus comes from participation

The native Box of the cell is $\Lambda(\psi,\nu) = T_A(\psi) - T_{\text{col}}(\nu)$, where
$T(k) = 2 - 2\cos(2\pi k/h) \in \mathbb{Z}$ (only for $h\in\{1,2,3,4,6\}$, i.e. the measured window
$d\in\{2,3\}$) is the lattice form — $T_A$ on the cell, $T_{\text{col}}$
on the column — and **`Box`** is the source's name for **their difference of two forms** (§8.3).
The minus is the sign of the Pontryagin dual of the column [T26], not a hand. The spectrum is strongly
asymmetric — mode count over the full space, $n_{\text{pos}}{:}n_{\text{neg}}{:}n_{\text{zero}} =
16{:}2{:}9$ at $d=2$, $219{:}9{:}28$ at $d=3$, with no $\lambda \leftrightarrow -\lambda$
symmetry. The derivation came out stronger than the expectation: with "+" the column does not
participate at all, since no zero has it active — both summands of the Box are $\ge 0$ and nothing
balances. So the minus is the only sign under which the cell's own column satisfies the law of
participation [T30].
T32 · S1011 · **J-0467**. Status: **candidate/cast**.

> **Boundaries:** "the minus is necessary" holds within the machinery of the native Box, by the M2
> contrast, and not as a ∀ theorem. The ancestors S956/S999/S1001/S1002/S1005/T26 are cited, not
> re-derived.

### 6.4 T39 — the sign is non-derivable

**Notation.** $m_0$ is the quantity of the §12 arc whose **sign** is the subject of this link (§8.3).
It carries **no number**: the value $m_0$ is a representative of the regime, **not canon**, and this
tome does not derive it [`A-stab`.scope · CODEX: "the graph nowhere holds the number $m_0^2$"].

$\mathrm{sign}(m_0)$ is $\varepsilon$-odd while the whole geometry of the machinery is
$\varepsilon$-even, so the sign does not follow from the construction. Otherwise $\pm m_0$ would be
inequivalent, the potential would carry an $\varepsilon$-odd term, and that contradicts the proved
evenness [S1034-i].

- **Witness multiplicity is 1, not 6.** All six parity measurements were made through one lens,
  $\varepsilon$-parity; by our own ancestor rule that is a single witness.
- The pattern of four failures — S1027 the minus homonym, S1033, S1036 leg 1 — was a theorem
  knocking, not a run of defeats.
- **Inheritance:** the choice is spontaneous exactly once per connected lattice [T37]. Local choices
  read the descending frame; they are not new bits.
- **The transmission mechanism is found:** a symmetric exchange $m_i \cdot m_j$ on the bipartite
  $A_d$ lattice, forced by symmetry and locality. Directionality from parent to child is withdrawn —
  the exchange is symmetric in $i \leftrightarrow j$, and "which came first" is an artefact of the
  order of construction. The debt closes by dissolution, not by postulation.
- **The through-going law: spontaneity ⟺ non-derivability.**

T39 · S1035·S1036·S1037·S1038·**S1043** · J-0488…**J-0497** · `MIRROR_ASSEMBLY_LAW.md`.
Status: **candidate/cast**.

> **Boundaries:** casts at $d=2,3$; the law holds within the $\varepsilon$-even machinery;
> bipartiteness ∀d is a property of the two-sublattice construction [T19], not a machine identity;
> **0 new tuned constants**, $J$ scaling in magnitude with the $\Lambda$ ruler, and the sign of $J$
> being a datum. The expectation "a third force is the relief" lost honestly: the relief exists —
> ridges at $d\ge3$, flat at $d=2$ — but it is $\varepsilon$-even and does not select the sign.
> There is no arbiter inside the geometry.

---

## 7 · WHAT THIS CHAIN DOES NOT CLOSE

Stated so that silence does not work as a claim.

| open, or behind the fence | where it stands |
|:--|:--|
| the seam between the two graph identifiers of $\Lambda$ (**E-1**) | Tome II; unifying the graphs is a separate pipeline |
| the two-component operator $\{I,H\}$, chirality, the cone | **Tome II** (T33) |
| $(3,1)$-non-splitting and centralisers | **Tome II** |
| dynamics, the variational principle | not on this floor (fence of floor −1) |
| non-cyclic discriminants, D type | a named contrast, not measured |
| the declared $d\ge3$ import in the premise $p\ge3$ [F4 · J-0410] | the $d{=}3$ disposition front (§6.2); no circularity |
| $d\ge4$ in the Lean core | a debt (§4.5) |
| any interpretive identification | **Tome III** |

---

## 8 · LEXICON (a word is admissible in exactly this sense)

> **Column 3 is not a status.** It answers one question of the floor fence: does the word have a
> referent on this floor. `⚓` is deliberately absent here — across the project `⚓` means measured,
> probe plus independent visa, and the same glyph in two roles would build exactly the homonym §8.1
> stands against. The referent marker is **`⊙`**: "the referent lives here".
> **`⊙` says nothing about the grade of proof** — grades live in §0 and in the status columns.

| word | sense **in this tome** | referent (⊙ = lives on this floor) |
|:--|:--|:--|
| **clock** | the class of the centre $\mathbb{Z}/(d+1)$ and its compact dual [T26.5/T34] | ⊙ derived here |
| **marker / dimer** | $u_0$, `AX-dimer` | ⊙ axiomatic input |
| **arena** | the Cartan subalgebra $A_d$ (`A-space`) | ⊙ derived here |
| **axis** | $\mathfrak{so}(2)_\alpha$ (`A-axis`) | ⊙ derived here |
| **sign** | $\mathrm{sign}(m_0)$, an $\varepsilon$-odd quantity [T39] | ⊙ derived here, as non-derivability |
| **arrow** | the unique free $\mathbb{Z}/2$ of the sign character of $D_h$ [T36] — one bit of realisation | ⊙ derived here |
| **democracy** | the full symmetric group $S_{d+1}$ | ⊙ axiomatic input |
| **cell** | $d+1$ axes, Gram $-1/d$ [T19] | ⊙ derived here (`AX-cell` is a theorem) |
| **ruler** | $\Lambda$, the only dimensional handle | ⊙ declared handle |

### 8.1 Separations against the homonym register

The register marks these words as forbidden bare. Tome I uses each in one sense only; no other sense
exists here.

| word | register knows | **here — exactly this object** | not used here | addr |
|:--|:-:|:--|:--|:--|
| **minus** | 3 | `minus_arith` — the first minus from participation: sign of the Pontryagin dual of the column, Box $T_A - T_{\text{col}}$ [T32, §6.3] | `minus_geom`, the arc of bits, $\mathfrak{so}(2)_\alpha$ — a homonym, §2 · `minus_chiral`, the seam minus [T33] — Tome II | S1027 · J-0481 · S1063-T2 |
| **sign** | 2 + ours | **the sign of the arrow** — the value of the free $\mathbb{Z}/2$ bit [T36, row "arrow" above], that is $\mathrm{sign}(m_0)$, $\varepsilon$-odd [T39, §6.4]. The tome coins **no new word**: both parts are already carved here | sign-MARKER, a convention · sign-CONVEXITY, another floor | ruling 2026-07-16 |
| **$\kappa$** | 2 + ours | **Schur $\kappa$**, the unique invariant form on $\mathfrak{sl}(n)$ [`A-kappa`, §2] | `kappa_conv` $=1$ · `kappa_stab` $=\Lambda$ [T38] | S1061 · §12-T38 |
| **closure** | 3 + ours | two uses, both qualified: the closure **method** for the order $\prec$ [§1.2], and **translational** closure [§4.2] | closure of bonds (FL-004) · closure of the metric (S886) · closure of a quantity in a volume, referent not derived | ruling 2026-07-16 |
| **beat** | 2 | not used at all — both objects of the register, the shift $S$ and the multiplier $U$, lie outside this tome | — | S1061 · J-0511 |
| **"3:1"** | 4 | not used as a ratio; the signature is written as $(3,1)$ with the side named | stacking projection $d{:}1$ · bond count $1{:}d$ · anisotropy parameter | S1054 · visa S1055 |

> **A correction kept visible.** The first edition of this tome coined its own word for the class of
> the centre — manufacturing a third homonym against a pair the register already carved. Withdrawn;
> the source's word, "clock", is used. No reader caught it: the check against the register did.
>
> **A second correction, same class, caught a different way.** Edition ver:2 named our object with a
> **three-letter abbreviation from another floor**: it (a) coined, **again**, a name the register does
> not hold — the register knows two objects — and (b) did so with a word whose expansion lives behind
> the fence. The abbreviation itself is **not reproduced here**: the fence grants no exemption for
> "mentioned in order to withdraw it", since that is the very loophole a word travels through.
> Withdrawn; the object is named from parts already carved in this tome. ★No machine check caught
> this one: the probe matched word stems, and an abbreviation has no stem — the class was invisible
> **by construction of the method**, not by oversight. A person found it (line B). The probe has
> carried a separate layer for abbreviations ever since (§10).

### 8.2 How numbers are written — `[address · unit · type]`

| number | form |
|:--|:--|
| $\det \mathrm{Cartan}(A_d) = d+1$ | [T19 · dimensionless · determinant] |
| $z = d+1$ | [T19 · dimensionless · neighbour count] |
| column period $= d+1$ | [T26.5 · **in bonds** · translation length] |
| dual circle $2\pi/(d+1)$ | [T26.5 · angular · period of the compact dual] |
| Gram $-1/d$; angle $\arccos(-1/d)$ | [T19 · dimensionless / angle · pairwise form] |
| $\mathrm{split} = 96\,(d{=}2)\,/\,264\,(d{=}3)$ | [T37 · dimensionless · count of the **counterfactual** rupture, not a quantity of the native state] |
| spectral asymmetry (mode count over full space, npos:nneg:nzero) $16{:}2{:}9$; $219{:}9{:}28$ | [T32 · dimensionless · mode count] |

### 8.3 Short carriers — expansions

The tome uses a few **short** names taken from the sources. Each is expanded here, because otherwise
the **reader** supplies the expansion, and the only one available to them lives on another floor.
Renaming a short name is not this tome's business: the name of a canon object is changed by a
ruling, not by an exposition.

| short | expansion in this tome | why the short name is dangerous | addr |
|:--|:--|:--|:--|
| `Box` | the **difference of two lattice forms** $\Lambda(\psi,\nu) = T_A(\psi) - T_{\text{col}}(\nu)$ [§6.3] | the same symbol on another floor denotes an operator, not a difference of forms | T32 · S1011 |
| $T_A$ | the lattice form **on the cell**: $T(k) = 2 - 2\cos(2\pi k/h) \in \mathbb{Z}$, only for $h\in\{1,2,3,4,6\}$ (measured window $d\in\{2,3\}$) | the letter $T$ denotes something else on another floor | T32 · S1011 |
| $T_{\text{col}}$ | the same form **on the column** | the same | T32 · S1011 |
| $m_0$ | the quantity of the §12 arc whose **sign** is the subject [T39]; the value $m_0$ is a representative of the regime, **not canon**, and is not derived here | the letter $m$ denotes something else on another floor | `A-stab` · T39 · §6.4 |

**Forbidden words.** The list lives in the checking probe; each has a referent on another floor.
Tome I does not use them, and any appearance fails the check.

---

## 9 · SOURCES

**The order of sources, by project ruling, 2026-07-25: the structured codex view first.**

| source | what is taken | coverage |
|:--|:--|:--|
| the structured codex view | §0 inputs and handles · §1 objects · §2 identities · §3.arc-0 / arc / −1 / −2 · §3.seam · §5 open items and fences | the § references above are its sections |
| the machine graph | 90 nodes, 4 roots, `deps` — machine topology of layers 0–7 | §1–§6 |
| the homonym register | separations of homonyms, and the form for writing numbers | **§8.1, §8.2** |
| seals of rulings | links absent from the graphs — hole **R-1**: cited directly | §3 (`LADDER_ANCESTOR_EXANTE.md`) · §6.4 (`MIRROR_ASSEMBLY_LAW.md`) |
| the register of verdicts | primary rulings T1–T39 with qualifiers verbatim — a qualifier dies in paraphrase, so it is quoted, not retold | every boundary block |

### 9.1 Published parts, with verified concept DOIs

| # | title | concept DOI | coverage |
|:-:|:--|:--|:--|
| 1 | *Order from closure on real quadratic spaces: terminal, minimality, orientation* | `10.5281/zenodo.21412389` | **T1–T11** — §6.1 |
| 2 | *Centralizers of nilpotent wedges in real orthogonal Lie algebras $\mathfrak{so}(p,q)$: a signature-resolved Levi-radical map* | `10.5281/zenodo.21429542` | T17–T19 — support for §4.2 |
| 3 | *The central charge of nilpotent-wedge strata in real $\mathfrak{so}(p,q)$* | `10.5281/zenodo.21437068` | T20–T25 — support for §1.1 |

> **DOI verification, closed 2026-07-25.** Opus 2 sat for a while in the memory index marked "concept DOI not
> yet verified", so it could not be cited as the address of the work. Two hands closed it: the
> Zenodo API (`conceptdoi` of record …602) and the carving in the citation file. Current is
> **concept `10.5281/zenodo.21429542`**; versions 1.0.0 = …543 and **1.0.1 = …602**, current. Cite
> the concept, not the version — a concept always resolves to the newest.

### 9.2 Not woven in

The interactive graph map goes out as a separate page, as is. It is deliberately not embedded in the
body of Tome I; only a pointer stands here.

---

## 10 · GATE

No status moved; no new numbers; every link carries an address.

The text stands under a machine fence: the **floor fence** (word stems, with separate UA and EN
rules) · **abbreviations** (exact token, case-sensitive — a word stem cannot see this class) ·
**expansions** of short carriers (§8.3) · carved terms and register homonyms (the register is read
as data) · the form in which numbers are written (a solvable tower over $\mathbb{Q}$). A check that
cannot fail is not a witness.

**Scope of the `_teeth` harness (honest, not universal):** the shared negative-control harness
covers the 15 probes at S1061 and later, from the 2026-07-22 carving onward. Earlier probes —
including every `measured` node of §4/§5/§6.1 — carry their own bespoke negative controls, built
before the harness existed, and are **not retrofitted**: a probe is the record of an act and is not
edited after the fact. The harness is a forward-looking policy, not a blanket guarantee over the
whole corpus; migrating an older probe under it requires its own ed.2 re-verification.
