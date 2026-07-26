# sigma_x (bipartite), sigma_z (chirality), and B (mirror) are three distinct Z_2 involutions

**Kind:** identity · **Status:** measured · **Address:** §2.4 · T33 · S1012 · J-0468

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Two-component structure forced (0 free parameters) {I,H}; v^2=(1/2)*trM=(d+1)^2/(2d) (d=2->9/4, d=3->8/3); sub-bridge #3(a) sealed](T33.md) | theorem | candidate |

## FUNCTION

Establishes that the three Z_2 involutions sigma_x, sigma_z, and B are mutually distinct on the two-component Hamiltonian H.

> **Boundary:** Two-component H level.

## OUTPUT

Provides the three-way Z_2 classification used in later two-component and arrow-realization arguments.

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


---

[← all nodes](index.md)
