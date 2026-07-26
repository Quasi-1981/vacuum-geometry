# (marked bond tensor orientation) equals a single D_h orbit, D_h=<c,w_0>

**Kind:** identity · **Status:** measured · **Address:** §2.5 · T34 · S1016 · J-0471

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Cell-locality theorem: rank-1 center implies a single clock; hole #1 dissolved](T34.md) | theorem | candidate |

## FUNCTION

Shows that the combination of a marked bond and an orientation forms a single orbit under D_h=<c,w_0>.

> **Boundary:** Cell-tensor-orientation level.

## OUTPUT

Grounds the cell/orientation orbit structure reused in the Coxeter-action seam seam.6.

**Consumed by:** nothing in this graph — a terminal node.

## Measurements

**S1016**

??? note "Probe code (`S1016.py`)"
    ```python
    --8<-- "src/probe/S1016.py"
    ```

**Result:** asserts_passed=20 · FAIL=0

??? note "Full run log (`S1016_run.log`)"
    ```text
    --8<-- "src/probe/S1016_run.log"
    ```


---

[← all nodes](index.md)
