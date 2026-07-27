# Commutant ladder with block ranks (1,m−1,d−m)

**Kind:** object · **Status:** measured · **Address:** §1 · T26.3 · S1000-T2

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [A_d cell as weights of the fundamental su(d+1) representation](OBJ-cell.md) | object | measured |

## FUNCTION

The commutant ladder has block ranks (1,m−1,d−m) for 1≤m≤d, with dimensions 1/2/2/1 (d=2), 1/2/3/2 (d=3), 1/2/3/3 (d≥4) — read off S1000_run.log.

> **Boundary:** Restricted to the A_d lattice; dimension sequence is d-dependent (not a single constant across d).

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

??? note "Run log (paraphrased under the word fence — see provenance) (`S1000_run.log`)"
    ```text
    --8<-- "src/probe/S1000_run.log"
    ```


---

[← all nodes](index.md)
