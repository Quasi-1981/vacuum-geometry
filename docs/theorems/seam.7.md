# Globality (S1018): pair-locality=1 leads to a single dial on the connected lattice

**Kind:** seam · **Status:** measured · **Address:** §3.seam.7 · S1018 · J-0473 · ruling §16-ruling

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Global clock (T34 to the lattice): pair-locality=1 implies a single dial on the connected lattice](T37.md) | theorem | candidate |

## FUNCTION

Traces how the T37 result of pair-locality equal to 1 produces a single dial across the connected lattice.

> **Boundary:** Levels -2 to -1; the same center Z/(d+1) applies across adjacent cells.

## OUTPUT

Feeds the global-clock conclusion used in T39's inheritance clause.

**Consumed by:** nothing in this graph — a terminal node.

## FACTORIZATION

A part of [Dihedral hinge D_h=⟨c,w₀⟩ of the −2↔−1 joint](OBJ-Dh.md) — component ⟨c⟩ (measured).

## Measurements

**S1018**

??? note "Probe code (`S1018.py`)"
    ```python
    --8<-- "src/probe/S1018.py"
    ```

**Result:** asserts_passed=15 · FAIL=0

??? note "Full run log (`S1018_run.log`)"
    ```text
    --8<-- "src/probe/S1018_run.log"
    ```


---

[← all nodes](index.md)
