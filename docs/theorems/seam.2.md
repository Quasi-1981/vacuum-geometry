# The sign of the first minus (T32) leads to the sign of the circle dual (Pontryagin dual)

**Kind:** seam · **Status:** measured · **Address:** §3.seam.2 · S1011 · J-0467 · ruling §10

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Period P=d+1 plus dual column equals the circle 2π/(d+1) (Pontryagin); 1A+1B+(d-1) holes](T26.5.md) | theorem | candidate |
| [First minus from participation + character bridge: minus equals the sign of the Pontryagin dual; hole #2 closed kinematically](T32.md) | theorem | candidate |

## FUNCTION

Traces how the sign established in T32 determines the sign of the circle dual, the Pontryagin dual.

> **Boundary:** Levels -2 to -1.

## OUTPUT

Feeds the center/dual identification recorded in id-2.1.

**Consumed by:** nothing in this graph — a terminal node.

## Measurements

**S1011**

??? note "Probe code (`S1011.py`)"
    ```python
    --8<-- "src/probe/S1011.py"
    ```

**Result:** asserts_passed=17 · FAIL=0

??? note "Full run log (`S1011_run.log`)"
    ```text
    --8<-- "src/probe/S1011_run.log"
    ```


---

[← all nodes](index.md)
