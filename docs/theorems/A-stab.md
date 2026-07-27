# Stabilization intrinsic: κ=Λ implies m₀>0 ⟺ Λ<a

**Kind:** theorem · **Status:** candidate · **Address:** §3.seam · MIRROR_ASSEMBLY_LAW.md · S1039 · J-0492

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [a>0 ⟹ m=0 unstable; monotonic trace forces runaway, no m₀² term](A-instability.md) | theorem | measured-negative |
| [Λ-scale: the sole dimensional free handle](AX-lambda.md) | input | axiom |
| [Threshold-weight law: t*=d (symbolic for all d); the Dirac case is live iff 0<t<d](T26.1.md) | theorem | candidate |
| [Two-weight discriminant law: t=s+(d-1) (symbolic for all d); Z2 symmetry exchanges s and t](T26.2.md) | theorem | candidate |

## FUNCTION

Stabilization is native: the detuning cost of the weights is κm² (T26.1/T26.2, second order in κ>0, ε-even), with κ=Λ as the single free handle; the separation equation 2κ=J(m₀) is refinement-convergent, giving the criterion m₀>0 if and only if Λ<a.

> **Boundary:** Restricted to the −2↔−1 joint; the value m₀ is a representative of the regime, not canonical; since κ=Λ, the reachability of Λ is honest (depends on AX-lambda).

## OUTPUT

Provides the stabilization mechanism and the criterion m₀>0 ⟺ Λ<a, resolving the instability recorded by the negative node upstream.

**Consumed by:** [Stage 2: the free ±m₀ choice is the single free bit (T36-bit)](A-ssb-bit.md) · [Registry wrapper (Section 12) over the A-* arc, not a separate derivation: mirror-assembly law (two-stage freezing + kappa=Lambda stabilization)](T38.md)

### 8.3 Notation — short carriers used above

| short | expansion | why the short name is dangerous |
|:--|:--|:--|
| `m_0` | the quantity of the arc whose sign is the subject; its value is a representative of the regime, not canon, and is not derived here | the letter m denotes something else on another floor |

## Measurements

**S1039**

??? note "Probe code (`S1039.py`)"
    ```python
    --8<-- "src/probe/S1039.py"
    ```

**Result:** exit=0

??? note "Full run log (`S1039_run.log`)"
    ```text
    --8<-- "src/probe/S1039_run.log"
    ```


---

[← all nodes](index.md)
