# κ unique: Jordan-Chevalley split — space (h) ⊥ axis (n±)

**Kind:** theorem · **Status:** measured · **Address:** §3.arc · J-0478 · arc S1023–S1026 · heat MIGRATION_ARC_BITS S1

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Space = Cartan algebra A_d via Schur-κ, for all d](A-space.md) | theorem | measured |
| [Simple algebra sl(n), n=d+1, with unique invariant form κ](OBJ-sln.md) | object | measured |

## FUNCTION

The form κ is unique; the Jordan-Chevalley decomposition splits it into h (semisimple, definite, identified with space) orthogonal to n± (nilpotent, indefinite, identified as the axis type), both arising from a single mechanism.

> **Boundary:** Single form κ; the space/axis split follows from one and the same mechanism only.

## OUTPUT

Yields the space-versus-axis orthogonal splitting used to build the space object (A-space) from the sl(n) algebra.

**Consumed by:** [κ: minimal parabolic p_α=(d,1), q=1, variant B; time=so(2)_α](A-axis.md)

## Measurements

**S1023**

??? note "Probe code (`S1023.py`)"
    ```python
    --8<-- "src/probe/S1023.py"
    ```

**Result:** exit=0

??? note "Full run log (`S1023_run.log`)"
    ```text
    --8<-- "src/probe/S1023_run.log"
    ```

**S1026**

??? note "Probe code (`S1026.py`)"
    ```python
    --8<-- "src/probe/S1026.py"
    ```

**Result:** exit=0

??? note "Full run log (`S1026_run.log`)"
    ```text
    --8<-- "src/probe/S1026_run.log"
    ```


---

[← all nodes](index.md)
