# No A/B discriminator exists (beta-inversion swaps them); the residual is 1 bit (FORK II)

**Kind:** identity · **Status:** measured · **Address:** §2.6 · T36 · S1017 · J-0472 · ruling §15

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [The arrow is a single bit of realization (5 legs + wedge witness ver:2): no A/B discriminator exists](T36.md) | theorem | candidate |

## FUNCTION

Confirms that beta-inversion swaps cases A and B so no discriminator between them can exist, leaving a residual of exactly 1 bit, labeled FORK II.

> **Boundary:** Level: side of the cell, modulo the lock.

## OUTPUT

Provides the 1-bit residual (FORK II) that underlies T36's arrow-as-single-bit conclusion.

**Consumed by:** no node lists this one among its premises. (That is a statement about the `deps` edges, which record ancestry to the forcing root rather than a chain — not a claim that nothing follows from it.)

## Measurements

**S1017**

??? note "Probe code (`S1017.py`)"
    ```python
    --8<-- "src/probe/S1017.py"
    ```

**Result:** asserts_passed=28 · FAIL=0

??? note "Full run log (`S1017_run.log`)"
    ```text
    --8<-- "src/probe/S1017_run.log"
    ```


---

[← all nodes](index.md)
