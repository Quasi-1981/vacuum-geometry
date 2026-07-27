# The sl-gl tower (T35) leads to the T26.3 commutant ladder

**Kind:** seam · **Status:** measured · **Address:** §3.seam.5 · S1015 · J-0470 · ruling §13-ruling

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Commutant ladder of labels: block ranks (1,m-1,d-m), dimensions 1/2/2/1 (d=2), 1/2/3/2 (d=3), 1/2/3/3 (d>=4)](T26.3.md) | theorem | candidate |
| [Tower of commutatives, anchored sl-gl ladder: N_iso is the ladder of the T26.3 commutant](T35.md) | theorem | candidate |

## FUNCTION

Traces how the sl-gl tower result of T35 produces the commutant ladder defined in T26.3.

> **Boundary:** Levels -2 to -1.

## OUTPUT

Establishes the commutant-ladder count (N_iso) used by T35.

**Consumed by:** nothing in this graph — a terminal node.

## Measurements

**S1015**

??? note "Probe code (`S1015.py`)"
    ```python
    --8<-- "src/probe/S1015.py"
    ```

**Result:** asserts_passed=10 · FAIL=0

??? note "Full run log (`S1015_run.log`)"
    ```text
    --8<-- "src/probe/S1015_run.log"
    ```


---

[← all nodes](index.md)
