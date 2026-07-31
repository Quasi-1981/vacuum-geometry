# Symbolic output of the probes

Each block below is printed **verbatim from the run log of the probe named beside it** — no line here was written by hand, and nothing is restated in prose. The probe and its log are published under `src/probe/`, so every block can be reproduced by running the file.

## T18 — Universal mechanism derived for all n: c=(sp2+so(eta|G)) semidirect h_heis; center from [module,module]

```text
--- rank-0 : so(η) block equations (=0) ---
  soη[A,]: G2 + G2.T = 0
  soη[B,]: H + D.T*G = 0
  soη[Cc,]: K + A.T = 0
```

*Source: probe `S952`, run log lines 2, 3, 4, 5 — see `src/probe/S952.py` and `src/probe/S952_run.log`.*

## T19 — Cell to A_d native map (Gram proportional to Cartan, z=d+1, d nodes); symbolic for all d (ver:2)

```text
BIT-FENCE d=2 (cross-check against the program's stamps):
  det(A_2)=3 · cos=−1/2(120°) · z=3 · nodes=2(K,K') · order=3 — all agree
```

*Source: probe `S956`, run log lines 14, 15 — see `src/probe/S956.py` and `src/probe/S956_run.log`.*

## T21 — Junction of T16 and T20: the cocycle w=eps_W (x) eta'_core is forced; Heisenberg part required to be metric

```text
--- ∀d SYMBOLIC (block form ω(E_{u,a},E_{v,b})=ε(u,v)·(G⁻¹)_{ab}) ---
  blockwise: ω(E_{u,a},E_{v,b})=ε(u,v)·(G⁻¹)_{ab} — a pure tensor for d=1,2,3 (±) → MATCH ∀d
```

*Source: probe `S980`, run log lines 14, 15 — see `src/probe/S980.py` and `src/probe/S980_run.log`.*

## T22 — Lambda-slot of the charge mu^2/nu equals the Lambda ruler; lambda_ext=-eta'; contraction unity of branches

```text
--- SCALE LAW (symbolic μ, ν) ---
  coefficient (base Z): -mu**2  → exponent μ: +2
```

*Source: probe `S983`, run log lines 5, 6 — see `src/probe/S983.py` and `src/probe/S983_run.log`.*

---

[← all nodes](theorems/index.md)
