# Scalar Box: Λ(ψ,ν)=T_A(ψ)−T_col(ν)

**Kind:** object · **Status:** measured · **Address:** §1 · T32 · S1011

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [A_d cell as weights of the fundamental su(d+1) representation](OBJ-cell.md) | object | measured |
| [Column/dimer: time-bond u₀ with translation symmetry P_col=d+1](OBJ-column.md) | object | measured |

## FUNCTION

The scalar Box function is Λ(ψ,ν)=T_A(ψ)−T_col(ν), with T(k)=2−2cos(2πk/P_col) an integer.

> **Boundary:** Restricted to L∈{3,4}.

## OUTPUT

Feeds into the two-component H object (OBJ-H) and the nodal set (OBJ-nodal).

**Consumed by:** [Two-component matrix H with invariants {I,H}](OBJ-H.md) · [Nodal set f=0 corresponds to d characters of ℤ/P_col](OBJ-nodal.md)

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
