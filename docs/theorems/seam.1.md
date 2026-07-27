# Seam joint of the -2 -> -1 break (T29 collapse) leads to loss of column-tick equality

**Kind:** seam · **Status:** measured · **Address:** §3.seam.1 · S1008 · J-0463 · ruling §7

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Collapse and three regimes: generic detuning yields 1 clock; criterion is Ctilde=1](T29.md) | theorem | candidate |

## FUNCTION

Records that, at the -2 to -1 joint (the T29 collapse), column-tick equality fails across the column.

> **Boundary:** Levels -2 to -1.

## OUTPUT

Marks the opening event of the -2/-1 arc that the later seams (seam.2 through seam.7) build on.

**Consumed by:** nothing in this graph — a terminal node.

## Measurements

**S1008**

??? note "Probe code (`S1008.py`)"
    ```python
    --8<-- "src/probe/S1008.py"
    ```

**Result:** asserts_passed=8 · FAIL=0

??? note "Full run log (`S1008_run.log`)"
    ```text
    --8<-- "src/probe/S1008_run.log"
    ```


---

[← all nodes](index.md)
