# Commutant ladder with block ranks (1,m−1,d−m)

**Kind:** object · **Status:** measured · **Address:** §1 · T26.3 · S1000-T2

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [A_d cell as weights of the fundamental su(d+1) representation](OBJ-cell.md) | object | measured |

## FUNCTION

The commutant ladder has block ranks (1,m−1,d−m), with dimensions 1,2,3,3.

> **Boundary:** Restricted to the A_d lattice.

## OUTPUT

Provides the block-rank decomposition used in the structure analysis of the A_d cell.

**Consumed by:** nothing in this graph — a terminal node.

## Measurements

**S1000**

??? note "Probe code (`S1000.py`)"
    ```python
    --8<-- "src/probe/S1000.py"
    ```

**Result:** asserts_passed=327 · FAIL=0

??? note "Full run log (`S1000_run.log`)"
    ```text
    --8<-- "src/probe/S1000_run.log"
    ```


---

[← all nodes](index.md)
