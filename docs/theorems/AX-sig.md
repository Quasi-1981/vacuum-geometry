# Signature (3,1) derived, forced by closure of T5 and T7

**Kind:** conclusion · **Status:** measured · **Address:** §3.arc-0 · §0 · T5(S908·J-0410) / T7(S903)

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Minimality implies (3,1): p>=3 and q>=1 implies n>=4; the minimum is uniquely (3,1); (2,2) fails](T5.md) | theorem | measured |
| [Square structures: A^2=+1 iff p=q; A^2=-1 iff p,q both even; implies (3,1) is rigid](T7.md) | theorem | measured |

## FUNCTION

The signature (3,1) is derived, forced by the closure combining T5 and T7; it is not a free input, as found by the author on 2026-07-19.

> **Boundary:** Derived result, not a premise.

## OUTPUT

Fixes the signature (3,1) as a downstream consequence rather than an assumed input.

**Consumed by:** nothing in this graph — a terminal node.

## PROOF COMPOSITION

The proof consumes one group of premises:

- [Minimality implies (3,1): p>=3 and q>=1 implies n>=4; the minimum is uniquely (3,1); (2,2) fails](T5.md)

**Independent witness groups:** 1

## Measurements

**S903**

??? note "Probe code (`S903.py`)"
    ```python
    --8<-- "src/probe/S903.py"
    ```

??? note "Full run log (`S903_run.log`)"
    ```text
    --8<-- "src/probe/S903_run.log"
    ```

**S908_1**

??? note "Probe code (`S908_1.py`)"
    ```python
    --8<-- "src/probe/S908_1.py"
    ```

??? note "Full run log (`S908_1_run.log`)"
    ```text
    --8<-- "src/probe/S908_1_run.log"
    ```

**S908_2**

??? note "Probe code (`S908_2.py`)"
    ```python
    --8<-- "src/probe/S908_2.py"
    ```

??? note "Full run log (`S908_2_run.log`)"
    ```text
    --8<-- "src/probe/S908_2_run.log"
    ```


---

[← all nodes](index.md)
