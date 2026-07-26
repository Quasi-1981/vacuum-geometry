# B = sigma_x composed with (k <-> -k) is the unique reflector; bare k<->-k or sigma_x alone are not symmetries

**Kind:** identity · **Status:** measured · **Address:** §2.3 · T33 · S1012/S1016 · J-0468/J-0471

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Two-component structure forced (0 free parameters) {I,H}; v^2=(1/2)*trM=(d+1)^2/(2d) (d=2->9/4, d=3->8/3); sub-bridge #3(a) sealed](T33.md) | theorem | candidate |

## FUNCTION

Shows that only the composite B=sigma_x o (k<->-k) acts as a genuine reflector, while the bare individual maps k<->-k and sigma_x are not symmetries on their own.

> **Boundary:** Native two-component level.

## OUTPUT

Grounds the two-component reflector structure reused in id-2.4's Z_2 classification.

**Consumed by:** nothing in this graph — a terminal node.

## Measurements

**S1012**

??? note "Probe code (`S1012.py`)"
    ```python
    --8<-- "src/probe/S1012.py"
    ```

**Result:** asserts_passed=13 · FAIL=0

??? note "Full run log (`S1012_run.log`)"
    ```text
    --8<-- "src/probe/S1012_run.log"
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
