# Nodal set f=0 corresponds to d characters of ℤ/h

**Kind:** object · **Status:** measured · **Address:** §1 · T32 · S1011

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Scalar Box: Λ(ψ,ν)=T_A(ψ)−T_col(ν)](OBJ-box.md) | object | measured |
| [Dual circle ℤ/h, Pontryagin dual of the column translation](OBJ-dual.md) | object | measured |

## FUNCTION

The nodal set defined by f=0 corresponds bijectively to d characters of the dual group ℤ/h.

> **Boundary:** Restricted to the A_d lattice.

## OUTPUT

Provides the nodal-set structure used in constructing the two-component object H (OBJ-H).

**Consumed by:** nothing in this graph — a terminal node.

### 8.3 Notation — short carriers used above

| short | expansion | why the short name is dangerous |
|:--|:--|:--|
| `Box` | the difference of two lattice forms, Lambda(psi,nu) = T_A(psi) - T_col(nu) | the same symbol on another floor denotes an operator, not a difference of forms |
| `T_A` | the lattice form on the cell, T(k) = 2 - 2cos(2 pi k / h), an integer only for h in {1,2,3,4,6} (measured window d in {2,3}) | the letter T denotes something else on another floor |
| `T_col` | the same lattice form on the column | the same |

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
