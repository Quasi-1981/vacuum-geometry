# A_d cell derived from alphabet (∀d, three independent routes)

**Kind:** theorem · **Status:** measured · **Address:** §0 · §3.arc · J-0476 · arc S1023–S1026 · C2_CELL_GATE_EXANTE.md (Path FAST)

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Space = Cartan algebra A_d via Schur-κ, for all d](A-space.md) | theorem | measured |

## FUNCTION

The A_d cell is derived from the alphabet via the Schur form of κ (J-0476, valid for all d), confirmed by three independent routes (rank, character, Molien), eliminating circularity; it was demoted from root-premise status at gate C2 (previously a premise as of S943/S946).

> **Boundary:** Holds for A_d, derived for all d; Lean confirmation steps L1–L3 proceed in parallel and do not block the root swap.

## OUTPUT

Provides the derived A_d cell object used throughout the object layer (OBJ-cell) as a consequence rather than an input.

**Consumed by:** [A_d cell as weights of the fundamental su(d+1) representation](OBJ-cell.md) · [Cell to A_d native map (Gram proportional to Cartan, z=d+1, d nodes); symbolic for all d (ver:2)](T19.md) · [Minimality implies (3,1): p>=3 and q>=1 implies n>=4; the minimum is uniquely (3,1); (2,2) fails](T5.md)

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
