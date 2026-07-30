# Coxeter ladder: Cartan(A_d) eigenvalues 2−2cos(πj/h_cox)

**Kind:** object · **Status:** measured · **Address:** §1 · T26 · S1004

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [A_d cell as weights of the fundamental su(d+1) representation](OBJ-cell.md) | object | measured |

## FUNCTION

The Coxeter ladder consists of the eigenvalues of Cartan(A_d), given by 2−2cos(πj/h_cox).

> **Boundary:** Restricted to the A_d lattice.

## OUTPUT

Provides the eigenvalue ladder used elsewhere in the A_d structure analysis.

**Consumed by:** nothing in this graph — a terminal node.

## Measurements

**S1004**

??? note "Probe code (`S1004.py`)"
    ```python
    --8<-- "src/probe/S1004.py"
    ```

**Result:** asserts_passed=30 · FAIL=0

??? note "Full run log (`S1004_run.log`)"
    ```text
    --8<-- "src/probe/S1004_run.log"
    ```


---

[← all nodes](index.md)
