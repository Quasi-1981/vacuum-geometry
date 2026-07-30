# Space = Cartan algebra A_d via Schur-κ, for all d

**Kind:** theorem · **Status:** measured · **Address:** §3.arc · J-0476 · arc S1023–S1026 · heat MIGRATION_ARC_BITS S1†

† the act behind this id is recorded internally; no probe is published under it, so the id is a reference, not a link.

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Λ-scale: the sole dimensional free handle](AX-lambda.md) | input | axiom |
| [Simple algebra sl(n), n=d+1, with unique invariant form κ](OBJ-sln.md) | object | measured |

## FUNCTION

From {Ω, the S_{d+1} democracy, Λ} via the Schur form of κ, space is identified with h_cart = A_d, with the Gram matrix proportional to Cartan(A_d) for all d.

> **Boundary:** Space realized as the Cartan algebra A_d; stands as a candidate reduction for AX-cell pending the C2 Lean gate ∀d.

## OUTPUT

Identifies space with the A_d Cartan algebra, providing the candidate reduction feeding AX-cell (pending the C2 Lean gate for all d).

**Consumed by:** [κ unique: Jordan-Chevalley split — space (h_cart) ⊥ axis (n±)](A-kappa.md) · [A_d cell derived from alphabet (∀d, three independent routes)](AX-cell.md)

## Measurements

**S1023**

??? note "Probe code (`S1023.py`)"
    ```python
    --8<-- "src/probe/S1023.py"
    ```

**Result:** exit=0

??? note "Full run log (`S1023_run.log`)"
    ```text
    --8<-- "src/probe/S1023_run.log"
    ```

**S1026**

??? note "Probe code (`S1026.py`)"
    ```python
    --8<-- "src/probe/S1026.py"
    ```

**Result:** exit=0

??? note "Full run log (`S1026_run.log`)"
    ```text
    --8<-- "src/probe/S1026_run.log"
    ```


---

[← all nodes](index.md)
