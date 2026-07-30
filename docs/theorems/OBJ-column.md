# Column/dimer: time-bond u₀ with translation symmetry P_col=d+1

**Kind:** object · **Status:** measured · **Address:** §1 · T26.5 · S1001

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Dimer axiom: the time-bond marked axis is 1 bit of 'time exists'](AX-dimer.md) | input | axiom |
| [A_d cell as weights of the fundamental su(d+1) representation](OBJ-cell.md) | object | measured |

## FUNCTION

The column (dimer) object is defined by the time-bond u₀, with translation as an exact symmetry of period P_col=d+1.

> **Boundary:** Restricted to the A_d lattice.

## OUTPUT

Supplies the dual circle object (OBJ-dual) via Pontryagin duality of the translation.

**Consumed by:** [Scalar Box: Λ(ψ,ν)=T_A(ψ)−T_col(ν)](OBJ-box.md) · [Dual circle ℤ/P_col, Pontryagin dual of the column translation](OBJ-dual.md)

### 8.3 Notation — short carriers used above

| short | expansion | why the short name is dangerous |
|:--|:--|:--|
| `Box` | the difference of two lattice forms, Lambda(psi,nu) = T_A(psi) - T_col(nu) | the same symbol on another floor denotes an operator, not a difference of forms |
| `T_A` | the lattice form on the cell, T(k) = 2 - 2cos(2 pi k / h), an integer only for h in {1,2,3,4,6} (measured window d in {2,3}) | the letter T denotes something else on another floor |
| `T_col` | the same lattice form on the column | the same |

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


---

[← all nodes](index.md)
