# Center Z/(d+1), the column dual, and Z/P_col are identified

**Kind:** identity · **Status:** measured · **Address:** §2.1 · T26.5/T32/T34 · S1001/S1011/S1013 · J-0456/J-0467/J-0469

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Period P_col=d+1 plus dual column equals the circle 2π/(d+1) (Pontryagin); 1A+1B+(d-1) holes](T26.5.md) | theorem | candidate |
| [Cartan torus: nodes equal the center Z/(d+1), the barycenters of the alcoves (symbolic for all d)](T26.7.md) | theorem | candidate |
| [First minus from participation + character bridge: minus equals the sign of the Pontryagin dual; hole #2 closed kinematically](T32.md) | theorem | candidate |
| [Cell-locality theorem: rank-1 center implies a single clock; hole #1 dissolved](T34.md) | theorem | candidate |

## FUNCTION

Establishes the three-way identity between the center Z/(d+1), the dual of the column, and Z/P_col.

> **Boundary:** Cell level, A_d.

## OUTPUT

Establishes a reusable three-way identification (center = column dual = Z/h) underlying the center-based reasoning in the seam analyses (seam.2, seam.3).

**Consumed by:** no node lists this one among its premises. (That is a statement about the `deps` edges, which record ancestry to the forcing root rather than a chain — not a claim that nothing follows from it.)

## Measurements

**S1001**

??? note "Probe code (`S1001.py`)"
    ```python
    --8<-- "src/probe/S1001.py"
    ```

**Result:** asserts_passed=48 · FAIL=0

??? note "Full run log (`S1001_run.log`)"
    ```text
    --8<-- "src/probe/S1001_run.log"
    ```

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

**S1013**

??? note "Probe code (`S1013.py`)"
    ```python
    --8<-- "src/probe/S1013.py"
    ```

**Result:** asserts_passed=23 · FAIL=0

??? note "Full run log (`S1013_run.log`)"
    ```text
    --8<-- "src/probe/S1013_run.log"
    ```


---

[← all nodes](index.md)
