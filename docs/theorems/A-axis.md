# κ: minimal parabolic p_α=(d,1), q=1, variant B; time=so(2)_α

**Kind:** theorem · **Status:** measured · **Address:** §3.arc · J-0479 · arc S1023–S1026 · heat MIGRATION_ARC_BITS stage 1

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [κ unique: Jordan-Chevalley split — space (h_cart) ⊥ axis (n±)](A-kappa.md) | theorem | measured |

## FUNCTION

For the minimal parabolic case p_α=(d,1) with q=1 (variant B), time is realized as the compact so(2)_α, and the label serves as a selector of the root α.

> **Boundary:** Restricted to the minimal parabolic p_α=(d,1), q=1, variant B.

## OUTPUT

Fixes the (d,1) signature on the parabolic with time identified as so(2)_α, giving the root-selector label used downstream.

**Consumed by:** nothing in this graph — a terminal node.

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
