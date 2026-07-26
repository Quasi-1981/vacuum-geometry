# w_0, the map k <-> -k, and the epsilon-flip are identified

**Kind:** identity · **Status:** measured · **Address:** §2.2 · T28 · S1007/S1016 · J-0462/J-0471

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Count of ie classes: C(q)=2^q, Ctilde=0/1/2/4; q=1 gives a unique class up to the arrow](T28.md) | theorem | candidate |

## FUNCTION

Establishes the identity between w_0, the momentum-reversal map k <-> -k, and the epsilon-flip.

> **Boundary:** Scalar/Weyl level (S_{d+1}).

## OUTPUT

Grounds the scalar/Weyl-level (S_{d+1}) structure used in T28 and reused in the class-count seam seam.6.

**Consumed by:** nothing in this graph — a terminal node.

## Measurements

**S1007**

??? note "Probe code (`S1007.py`)"
    ```python
    --8<-- "src/probe/S1007.py"
    ```

**Result:** asserts_passed=13 · FAIL=0

??? note "Full run log (`S1007_run.log`)"
    ```text
    --8<-- "src/probe/S1007_run.log"
    ```

**S1016**

??? note "Probe code (`S1016.py`)"
    ```python
    --8<-- "src/probe/S1016.py"
    ```

**Result:** asserts_passed=20 · FAIL=0

??? note "Full run log (`S1016_run.log`)"
    ```text
    --8<-- "src/probe/S1016_run.log"
    ```


---

[← all nodes](index.md)
