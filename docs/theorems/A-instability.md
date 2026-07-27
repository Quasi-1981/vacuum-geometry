# a>0 ⟹ m=0 unstable; monotonic trace forces runaway, no m₀² term

**Kind:** theorem · **Status:** measured-negative · **Address:** §3.seam · MIRROR_ASSEMBLY_LAW.md · S1033-Q1 + correction J-0486

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Stage 1: T29 collapse (ε-even, zero separation) as precondition](A-collapse-pre.md) | theorem | candidate |
| [Two-component structure forced (0 free parameters) {I,H}; v^2=(1/2)*trM=(d+1)^2/(2d) (d=2->9/4, d=3->8/3); sub-bridge #3(a) sealed](T33.md) | theorem | candidate |

## FUNCTION

For a>0 the m=0 state is unstable because the trace favors a nonzero lowest value; since the full trace V=−Σ√(|f|²+m²) is strictly monotonic, the system runs away without a stabilizing term, and this negative result carries no m₀² contribution (a lattice artifact of the truncated series, per J-0486).

> **Boundary:** Negative result: instability without a number; the finding 'negative carries no m₀²' is load-bearing for J-0486 and is kept separate from A-stab.

## OUTPUT

Establishes a negative result — instability without producing a numeric value — and rules out any m₀² contribution, to be kept separate from the A-stab result.

**Consumed by:** [Stabilization intrinsic: κ=Λ implies m₀>0 ⟺ Λ<a](A-stab.md)

### 8.3 Notation — short carriers used above

| short | expansion | why the short name is dangerous |
|:--|:--|:--|
| `m_0` | the quantity of the arc whose sign is the subject; its value is a representative of the regime, not canon, and is not derived here | the letter m denotes something else on another floor |

## Measurements

**S1033**

??? note "Probe code (`S1033.py`)"
    ```python
    --8<-- "src/probe/S1033.py"
    ```

**Result:** exit=0

??? note "Full run log (`S1033_run.log`)"
    ```text
    --8<-- "src/probe/S1033_run.log"
    ```


---

[← all nodes](index.md)
