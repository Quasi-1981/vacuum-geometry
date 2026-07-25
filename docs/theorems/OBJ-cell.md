# A_d cell as weights of the fundamental su(d+1) representation

**Kind:** object · **Status:** measured · **Address:** §1 · T26/T19 · S956

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [A_d cell derived from alphabet (∀d, three independent routes)](AX-cell.md) | theorem | measured |
| [Cell to A_d native map (Gram proportional to Cartan, z=d+1, d nodes); symbolic for all d (ver:2)](T19.md) | theorem | measured |

## FUNCTION

The A_d cell consists of the weights of the fundamental representation of su(d+1), with Gram matrix proportional to SC(I+11ᵀ).

> **Boundary:** Holds for A_d, symbolically for all d.

## OUTPUT

Provides the base cell object used to build the commutant ladder, the Coxeter ladder, the column/dimer, and the Box object.

**Consumed by:** [Scalar Box: Λ(ψ,ν)=T_A(ψ)−T_col(ν)](OBJ-box.md) · [Column/dimer: time-bond u₀ with translation symmetry P=d+1](OBJ-column.md) · [Commutant ladder with block ranks (1,m−1,d−m)](OBJ-commutant.md) · [Coxeter ladder: Cartan(A_d) eigenvalues 2−2cos(πj/h)](OBJ-coxeter.md)

### 8.3 Notation — short carriers used above

| short | expansion | why the short name is dangerous |
|:--|:--|:--|
| `Box` | the difference of two lattice forms, Lambda(psi,nu) = T_A(psi) - T_col(nu) | the same symbol on another floor denotes an operator, not a difference of forms |
| `T_A` | the lattice form on the cell, T(k) = 2 - 2cos(2 pi k / h), an integer | the letter T denotes something else on another floor |
| `T_col` | the same lattice form on the column | the same |

---

[← all nodes](index.md)
