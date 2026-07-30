# Center Z/(d+1) (T32 nodes + T34 locality) leads to nodes = center = barycenters of T26.7

**Kind:** seam · **Status:** measured · **Address:** §3.seam.3 · S1011/S1013 · J-0467/J-0469 · rulings §10/§12-ruling

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Cartan torus: nodes equal the center Z/(d+1), the barycenters of the alcoves (symbolic for all d)](T26.7.md) | theorem | candidate |
| [First minus from participation + character bridge: minus equals the sign of the Pontryagin dual; hole #2 closed kinematically](T32.md) | theorem | candidate |
| [Cell-locality theorem: rank-1 center implies a single clock; hole #1 dissolved](T34.md) | theorem | candidate |

## FUNCTION

Shows that combining the T32 node structure with the T34 locality result yields the identification of nodes with the center and with the barycenters from T26.7.

> **Boundary:** Levels -2 to -1.

## OUTPUT

Feeds the center identification recorded in id-2.1.

**Consumed by:** nothing in this graph — a terminal node.

## FACTORIZATION

A part of [Dihedral hinge D_h=⟨c,w₀⟩ of the −2↔−1 joint](OBJ-Dh.md) — component ⟨c⟩ (measured).

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
