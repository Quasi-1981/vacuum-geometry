# Symbolic output of the probes

Each block below is printed **verbatim from the run log of the probe named beside it** — no line here was written by hand, and nothing is restated in prose. The probe and its log are published under `src/probe/`, so every block can be reproduced by running the file.

## A-axis

```text
   space=a Schur-form on h [S1023] · axis-TYPE=the nilpotent-half of the same κ [S1025] ·
```

*Source: probe `S1026`, run log lines 78 — see `src/probe/S1026.py` and `src/probe/S1026_run.log`. Node: [A-axis](theorems/A-axis.md).*

## A-collapse-pre

```text
  V(m)=const −a·m²+b·m⁴, a,b = Σ over the native band |f(k)| (H=[[0,f],[f̄,0]], native momenta).
```

*Source: probe `S1034`, run log lines 18 — see `src/probe/S1034.py` and `src/probe/S1034_run.log`. Node: [A-collapse-pre](theorems/A-collapse-pre.md).*

## A-instability

```text
  V(m) = −(1/N)Σ√(|f|²+m²) = const − a·m² + b·m⁴ (expansion). a>0 (m=0 unstable), b>0 (stabilizing).
```

*Source: probe `S1033`, run log lines 10 — see `src/probe/S1033.py` and `src/probe/S1033_run.log`. Node: [A-instability](theorems/A-instability.md).*

## A-inherit

```text
  ⟹ d=2: 2π/60°=6 an integer ⟹ FLAT (deficit 0) · d≥3: 2π/θ is NOT an integer ⟹ deficit ≠0 (BUMPS exist).
```

*Source: probe `S1037`, run log lines 15 — see `src/probe/S1037.py` and `src/probe/S1037_run.log`. Node: [A-inherit](theorems/A-inherit.md).*

## A-kappa

```text
    axis-TYPE = the Schur-κ split [S1025] · axis-FORM (·,1) = the parabolic · UNIQUENESS = the mark [this one].
```

*Source: probe `S1026`, run log lines 43 — see `src/probe/S1026.py` and `src/probe/S1026_run.log`. Node: [A-kappa](theorems/A-kappa.md).*

## A-nonderiv

```text
   sign(m₀)=sign⟨σ_z⟩  | −1 ODD  (σ_x·σ_z·σ_x=−σ_z ⟹ ε-ODD)
```

*Source: probe `S1036`, run log lines 11 — see `src/probe/S1036.py` and `src/probe/S1036_run.log`. Node: [A-nonderiv](theorems/A-nonderiv.md).*

## A-space

```text
  ⟹ the Gram = (s²/2)(I − 𝟙𝟙ᵀ/(N+1)), rank N, PSD ⟹ the configuration is UNIQUE up to isometry+scale
```

*Source: probe `S1023`, run log lines 76 — see `src/probe/S1023.py` and `src/probe/S1023_run.log`. Node: [A-space](theorems/A-space.md).*

## A-ssb-bit

```text
  σ_z: the c-action = INVARIANT (c→+1) [S1032] · the w₀/B-action = NEGATION (w₀→−1) ⟹ REP = the sign-character D_h.
```

*Source: probe `S1035`, run log lines 18 — see `src/probe/S1035.py` and `src/probe/S1035_run.log`. Node: [A-ssb-bit](theorems/A-ssb-bit.md).*

## A-stab

```text
(β) THE BALANCE V_tot=κm²−Σ√ + ★A REFINEMENT-TEST L→2L→4L (mandatory, the J-0486 lesson)
```

*Source: probe `S1039`, run log lines 18 — see `src/probe/S1039.py` and `src/probe/S1039_run.log`. Node: [A-stab](theorems/A-stab.md).*

## A-time-neg

```text
    (iii) the perm.-repr = triv ⊕ W, dim End = Σm² = 2 ⟹ W is IRREDUCIBLE, multiplicity 1, W≇triv;
```

*Source: probe `S1023`, run log lines 39 — see `src/probe/S1023.py` and `src/probe/S1023_run.log`. Node: [A-time-neg](theorems/A-time-neg.md).*

## AX-cell

```text
    (iv) J|_W = 𝟙𝟙ᵀ|_W = 0 (W ⊥ 𝟙) ⟹ on W the inv. form = a·I ⟹ dim = 1.  ∎ ∀d
```

*Source: probe `S1023`, run log lines 40 — see `src/probe/S1023.py` and `src/probe/S1023_run.log`. Node: [AX-cell](theorems/AX-cell.md).*

## AX-indef

```text
     so(3,0): dim J = 3 · dim K = 0 (= p·q = 0) · together 3 = n(n−1)/2  [pt.2 · dimensionless · a count of generators]
```

*Source: probe `S1066_1`, run log lines 20 — see `src/probe/S1066_1.py` and `src/probe/S1066_1_run.log`. Node: [AX-indef](theorems/AX-indef.md).*

## AX-sig

```text
    LINE (from the identities above): A² = −𝟙 ⟹ AᵀηA = η and xᵀηAx = 0 ∀x ⟹ the Gram of the pair (x, Ax) = [[t,0],[0,t]], t = xᵀηx; p mod 2 = 1 · q mod 2 = 1
```

*Source: probe `S903`, run log lines 114 — see `src/probe/S903.py` and `src/probe/S903_run.log`. Node: [AX-sig](theorems/AX-sig.md).*

## OBJ-Dh

```text
   d=2 (even): prediction shov.4/6=dihedral · actual=dihedral ⟹ MATCH ✓
```

*Source: probe `S1028`, run log lines 19 — see `src/probe/S1028.py` and `src/probe/S1028_run.log`. Node: [OBJ-Dh](theorems/OBJ-Dh.md).*

## OBJ-H

```text
  T3 (cone): on-shell T_col(ν)=|f|² near the node ⟹ v² = transverse-curvature(S999:SC·n/2=9/4)
```

*Source: probe `S1012`, run log lines 17 — see `src/probe/S1012.py` and `src/probe/S1012_run.log`. Node: [OBJ-H](theorems/OBJ-H.md).*

## OBJ-cell

```text
  det(A_2)=3 · cos=−1/2(120°) · z=3 · nodes=2(K,K') · order=3 — all agree
```

*Source: probe `S956`, run log lines 15 — see `src/probe/S956.py` and `src/probe/S956_run.log`. Node: [OBJ-cell](theorems/OBJ-cell.md).*

## OBJ-sln

```text
  ⟹ dim=1 ∀ measured d; the orbital argument is n-independent ⟹ K1 does NOT fire.
```

*Source: probe `S1023`, run log lines 51 — see `src/probe/S1023.py` and `src/probe/S1023_run.log`. Node: [OBJ-sln](theorems/OBJ-sln.md).*

## T1

```text
   (1) [A,B] = Σ (a×b)_i J_i ; det[a;b;a×b] = |a×b|² = |a|²|b|² − (a·b)² (sympy-exact).
```

*Source: probe `S899`, run log lines 65 — see `src/probe/S899.py` and `src/probe/S899_run.log`. Node: [T1](theorems/T1.md).*

## T10

```text
  ISOMETRIES SᵀηS = η (rationally exact: rotations (c,s) with c²+s²=1 · boosts (ch,sh) with ch²−sh²=1 · reflections · composites):
```

*Source: probe `S906`, run log lines 18 — see `src/probe/S906.py` and `src/probe/S906_run.log`. Node: [T10](theorems/T10.md).*

## T11

```text
    sign-swap P: 0↔3 · 1↔4 · 2↔5 · det P = -1 · Pᵀ = P = P⁻¹
```

*Source: probe `S910`, run log lines 57 — see `src/probe/S910.py` and `src/probe/S910_run.log`. Node: [T11](theorems/T11.md).*

## T15

```text
  m3: BROKEN isometry check (Eᵀ·E = 𝟙 instead of Eᵀ·η·E = η) on a boost.
```

*Source: probe `S920`, run log lines 195 — see `src/probe/S920.py` and `src/probe/S920_run.log`. Node: [T15](theorems/T15.md).*

## T18

```text
          Q with N3·Q=0 [only ker N3 row ⟹ dim d' = module]
```

*Source: probe `S952`, run log lines 61 — see `src/probe/S952.py` and `src/probe/S952_run.log`. Node: [T18](theorems/T18.md).*

## T19

```text
  d=2: Cartan det3 · adj=-3/2 · cos2θ=-1/2 · z=3 · nodes=2 — all agree with S956
```

*Source: probe `S960`, run log lines 25 — see `src/probe/S960.py` and `src/probe/S960_run.log`. Node: [T19](theorems/T19.md).*

## T2

```text
      dim L = 3 · Killing (0,3,0) · L ⊆ span{all J}: yes · L = so(3,0): yes · pairs: 3
```

*Source: probe `S900`, run log lines 24 — see `src/probe/S900.py` and `src/probe/S900_run.log`. Node: [T2](theorems/T2.md).*

## T20

```text
  structurally: ω=[[0,G⁻¹],[−G⁻¹,0]] ⟹ rank=2·rank(G⁻¹)=2d, nondeg, k=d ∀d
```

*Source: probe `S961`, run log lines 109 — see `src/probe/S961.py` and `src/probe/S961_run.log`. Node: [T20](theorems/T20.md).*

## T21

```text
  (4,2): ω = ⊕_a (const·η_aa·ε) block-diagonal, cross-blocks=0 → transport-match=YES
```

*Source: probe `S980`, run log lines 18 — see `src/probe/S980.py` and `src/probe/S980_run.log`. Node: [T21](theorems/T21.md).*

## T22

```text
  coeff(μ,ν) = (-1)·μ²/ν ; the condition coeff=1 ⟹ ν = (-1)·μ² (here ν=−μ²)
```

*Source: probe `S983`, run log lines 13 — see `src/probe/S983.py` and `src/probe/S983_run.log`. Node: [T22](theorems/T22.md).*

## T26

```text
  (d=2→6, d=3→24, d=4→120). Weyl reflections preserve the Gram ⟹ all alcoves are CONGRUENT.
```

*Source: probe `S1002`, run log lines 34 — see `src/probe/S1002.py` and `src/probe/S1002_run.log`. Node: [T26](theorems/T26.md).*

## T28

```text
  ★MEASUREMENT: C̃(0)=0 (degenerate) · C̃(1)=1 · C̃(2)=2 · C̃(3)=4 — GROWS. Only q=1 has
```

*Source: probe `S1007`, run log lines 13 — see `src/probe/S1007.py` and `src/probe/S1007_run.log`. Node: [T28](theorems/T28.md).*

## T29

```text
  law: q_eff=2 ⟺ w₁+w₂ ∈ {{1,2,3,4}} (2(w₁+w₂)=T-sum∈{{0,2,4,6,8}}) — the weight SUM is integer.
```

*Source: probe `S1008`, run log lines 30 — see `src/probe/S1008.py` and `src/probe/S1008_run.log`. Node: [T29](theorems/T29.md).*

## T3

```text
      dim L = 3 · Killing (2,1,0) · L ⊆ span{all J}: no · L = so(2,1): yes · pairs: 1
```

*Source: probe `S900`, run log lines 52 — see `src/probe/S900.py` and `src/probe/S900_run.log`. Node: [T3](theorems/T3.md).*

## T33

```text
  T3 (cone): on-shell T_col(ν)=|f|² near the node ⟹ v² = transverse-curvature(S999:SC·n/2=8/3)
```

*Source: probe `S1012`, run log lines 32 — see `src/probe/S1012.py` and `src/probe/S1012_run.log`. Node: [T33](theorems/T33.md).*

## T34

```text
  d=2: center Z/3; pairwise bond differences integer=True; bond nontrivial=True ⟹ 2 bonds → ONE center class (not 2 independent columns)
```

*Source: probe `S1013`, run log lines 17 — see `src/probe/S1013.py` and `src/probe/S1013_run.log`. Node: [T34](theorems/T34.md).*

## T35

```text
  boundaries m=0/m=d+1 → N_iso=0 (pre-reg «no pairs»; Pillar-3 q=0⟹no time): True
```

*Source: probe `S1015`, run log lines 31 — see `src/probe/S1015.py` and `src/probe/S1015_run.log`. Node: [T35](theorems/T35.md).*

## T36

```text
  native d=2 (h=3): f(−k)=conj f(k) exactly ∀9 native momenta — H(−k)=σ_x H σ_x ✓ (bare k→−k is not a symmetry)
```

*Source: probe `S1016`, run log lines 39 — see `src/probe/S1016.py` and `src/probe/S1016_run.log`. Node: [T36](theorems/T36.md).*

## T37

```text
  d=2: center Z/3; column-1 (cell-1) charge=1 · column-2 (neighboring cell, a different bond) charge=1 ⟹ rank⟨1,1⟩=1 ⟹ ★THE PAIR'S CAPACITY = 1 (a merge forced by the center)
```

*Source: probe `S1018`, run log lines 16 — see `src/probe/S1018.py` and `src/probe/S1018_run.log`. Node: [T37](theorems/T37.md).*

## T38

```text
    B·m = σ_x·σ_z·σ_x = −σ_z = −m (ε-odd). Schur-style uniqueness: the 1-dim sign-isotype
```

*Source: probe `S1032`, run log lines 25 — see `src/probe/S1032.py` and `src/probe/S1032_run.log`. Node: [T38](theorems/T38.md).*

## T39

```text
   2. REP: both = the sign-character D_h (c→+1,w₀→−1, 1-dim unique) — THE SAME. ✓
```

*Source: probe `S1035`, run log lines 62 — see `src/probe/S1035.py` and `src/probe/S1035_run.log`. Node: [T39](theorems/T39.md).*

## T5

```text
    ⟹ a×b=0 ⟹ a∥b ⟹ dim≤1. ⟹ T1 holds NOT "because so(3)", but because the METRIC
```

*Source: probe `S908_1`, run log lines 9 — see `src/probe/S908_1.py` and `src/probe/S908_1_run.log`. Node: [T5](theorems/T5.md).*

## T6

```text
      dim g_B = 0 · Killing (0,0,0) · degen: — · key {dim 0 · Killing (0,0,0)} in the set of diagonal classes A(n=2): yes
```

*Source: probe `S902`, run log lines 105 — see `src/probe/S902.py` and `src/probe/S902_run.log`. Node: [T6](theorems/T6.md).*

## T7

```text
    LINE (from the identities above): A² = −𝟙 ⟹ AᵀηA = η and xᵀηAx = 0 ∀x ⟹ the Gram of the pair (x, Ax) = [[t,0],[0,t]], t = xᵀηx; p mod 2 = 0 · q mod 2 = 0
```

*Source: probe `S903`, run log lines 47 — see `src/probe/S903.py` and `src/probe/S903_run.log`. Node: [T7](theorems/T7.md).*

## T8

```text
      dim g_B = 0 · Killing (0,0,0) · degen: — · key {dim 0 · Killing (0,0,0)} in the set of diagonal classes A(n=2): yes
```

*Source: probe `S902`, run log lines 110 — see `src/probe/S902.py` and `src/probe/S902_run.log`. Node: [T8](theorems/T8.md).*

## T9

```text
    A = -2·J(0,1) -2·J(2,3) -2·K(0,2) -1·K(0,3) -1·K(1,2) +2·K(1,3)
```

*Source: probe `S905`, run log lines 209 — see `src/probe/S905.py` and `src/probe/S905_run.log`. Node: [T9](theorems/T9.md).*

## id-2.1

```text
  d=3: center Z/4; pairwise bond differences integer=True; bond nontrivial=True ⟹ 2 bonds → ONE center class (not 2 independent columns)
```

*Source: probe `S1013`, run log lines 18 — see `src/probe/S1013.py` and `src/probe/S1013_run.log`. Node: [id-2.1](theorems/id-2.1.md).*

## id-2.2

```text
  native d=3 (h=4): f(−k)=conj f(k) exactly ∀64 native momenta — H(−k)=σ_x H σ_x ✓ (bare k→−k is not a symmetry)
```

*Source: probe `S1016`, run log lines 40 — see `src/probe/S1016.py` and `src/probe/S1016_run.log`. Node: [id-2.2](theorems/id-2.2.md).*

## id-2.3

```text
  basis = {{I (scalar·identity), H=[[0,f],[f̄,0]] (off-diag ∝ f)}}; σ_z: KILLED (q≡0).
```

*Source: probe `S1012`, run log lines 11 — see `src/probe/S1012.py` and `src/probe/S1012_run.log`. Node: [id-2.3](theorems/id-2.3.md).*

## id-2.4

```text
  basis = {{I (scalar·identity), H=[[0,f],[f̄,0]] (off-diag ∝ f)}}; σ_z: KILLED (q≡0).
```

*Source: probe `S1012`, run log lines 26 — see `src/probe/S1012.py` and `src/probe/S1012_run.log`. Node: [id-2.4](theorems/id-2.4.md).*

## id-2.5

```text
  d=2,n=3: ABS row c=(0, 1, 2) — chiral(oriented)=True · REL increment c=(1, 1, 1),c⁻¹=(2, 2, 2) — blind-under-negation=True
```

*Source: probe `S1016`, run log lines 49 — see `src/probe/S1016.py` and `src/probe/S1016_run.log`. Node: [id-2.5](theorems/id-2.5.md).*

## id-2.6

```text
  d=2, h=3: sites(A=0,B=1)=[0, 1] · holes=[2] · period=3 · 2u₀=2∈holes ✓
```

*Source: probe `S1017`, run log lines 11 — see `src/probe/S1017.py` and `src/probe/S1017_run.log`. Node: [id-2.6](theorems/id-2.6.md).*

## seam.3

```text
  d=2: C_d-slice (s space-scalar, t time); resonance-form 2·s − t (codim 1); the screw protects s-as-a-whole, t = the ONLY unprotected ⟹ the address of the break = the column weight
```

*Source: probe `S1013`, run log lines 55 — see `src/probe/S1013.py` and `src/probe/S1013_run.log`. Node: [seam.3](theorems/seam.3.md).*

## seam.4

```text
    / column-curvature(T_col''(0)) = 81/8 — (2π)² CANCELS ⟹ v²=h⁴/(4d)=81/8 (0 handles).
```

*Source: probe `S1012`, run log lines 18 — see `src/probe/S1012.py` and `src/probe/S1012_run.log`. Node: [seam.4](theorems/seam.4.md).*

## seam.6

```text
  d=3,n=4: ABS row c=(0, 1, 2, 3) — chiral(oriented)=True · REL increment c=(1, 1, 1, 1),c⁻¹=(3, 3, 3, 3) — blind-under-negation=True
```

*Source: probe `S1016`, run log lines 50 — see `src/probe/S1016.py` and `src/probe/S1016_run.log`. Node: [seam.6](theorems/seam.6.md).*

## seam.7

```text
  d=3: center Z/4; column-1 (cell-1) charge=1 · column-2 (neighboring cell, a different bond) charge=1 ⟹ rank⟨1,1⟩=1 ⟹ ★THE PAIR'S CAPACITY = 1 (a merge forced by the center)
```

*Source: probe `S1018`, run log lines 17 — see `src/probe/S1018.py` and `src/probe/S1018_run.log`. Node: [seam.7](theorems/seam.7.md).*


### 8.3 Notation — short carriers used above

| short | expansion | why the short name is dangerous |
|:--|:--|:--|
| `T_col` | the same lattice form on the column | the same |
| `m_0` | the quantity of the arc whose sign is the subject; its value is a representative of the regime, not canon, and is not derived here | the letter m denotes something else on another floor |

---

[← all nodes](theorems/index.md)
