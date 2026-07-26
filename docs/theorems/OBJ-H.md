# Two-component matrix H with invariants {I,H}

**Kind:** object · **Status:** measured · **Address:** §1 · T33 · S1012

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Scalar Box: Λ(ψ,ν)=T_A(ψ)−T_col(ν)](OBJ-box.md) | object | measured |

## FUNCTION

H is the two-component object [[0,f],[f̄,0]]; together with the identity I, {I,H} are the only 2×2 invariants.

> **Boundary:** Restricted to A_d with d∈{2,3}.

## OUTPUT

Provides the two-component structure used to define the nodal set (OBJ-nodal) via f=0.

**Consumed by:** nothing in this graph — a terminal node.

### 8.3 Notation — short carriers used above

| short | expansion | why the short name is dangerous |
|:--|:--|:--|
| `Box` | the difference of two lattice forms, Lambda(psi,nu) = T_A(psi) - T_col(nu) | the same symbol on another floor denotes an operator, not a difference of forms |
| `T_A` | the lattice form on the cell, T(k) = 2 - 2cos(2 pi k / h), an integer | the letter T denotes something else on another floor |
| `T_col` | the same lattice form on the column | the same |

## Measurements

**S1012**

??? note "Probe code (`S1012.py`)"
    ```python
    --8<-- "src/probe/S1012.py"
    ```

**Result:** asserts_passed=13 · FAIL=0

??? note "Full run log (`S1012_run.log`)"
    ```text
    --8<-- "src/probe/S1012_run.log"
    ```


---

[← all nodes](index.md)
