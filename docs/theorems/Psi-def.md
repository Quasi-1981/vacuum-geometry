# Ϸ-definition of time and its dual quantity via the time-step↔circle pair

**Kind:** definition · **Status:** author-word · **Address:** §3.-1 · S998/S1000/S1000-T2/S597†/S1001/S1002 · J-0453…J-0457

† the act behind this id is recorded internally; no probe is published under it, so the id is a reference, not a link.

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [W40 re-reading of the cell's time axis: umbrella of 7 laws (see T26.1-T26.7)](T26.md) | theorem | candidate |
| [Period P_col=d+1 plus dual column equals the circle 2π/(d+1) (Pontryagin); 1A+1B+(d-1) holes](T26.5.md) | theorem | candidate |

## FUNCTION

Ϸ defines time together with its Pontryagin-dual quantity, via the pair time-step↔circle (per the author's directive to cut it out).

> **Boundary:** Restricted to A_d; when q=0, neither time nor its dual quantity exists.

## OUTPUT

Supplies the time/dual-quantity Pontryagin pair definition used downstream wherever q=0 forces both to vanish.

**Consumed by:** nothing in this graph — a terminal node.

## PROOF COMPOSITION

The proof consumes one group of premises:

- [Period P_col=d+1 plus dual column equals the circle 2π/(d+1) (Pontryagin); 1A+1B+(d-1) holes](T26.5.md)

**Independent witness groups:** 1

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

**S1001**

??? note "Probe code (`S1001.py`)"
    ```python
    --8<-- "src/probe/S1001.py"
    ```

**Result:** asserts_passed=48 · FAIL=0

??? note "Full run log (`S1001_run.log`)"
    ```text
    --8<-- "src/probe/S1001_run.log"
    ```

**S1002**

??? note "Probe code (`S1002.py`)"
    ```python
    --8<-- "src/probe/S1002.py"
    ```

**Result:** asserts_passed=36 · FAIL=0

??? note "Full run log (`S1002_run.log`)"
    ```text
    --8<-- "src/probe/S1002_run.log"
    ```

**S998**

??? note "Probe code (`S998.py`)"
    ```python
    --8<-- "src/probe/S998.py"
    ```

**Result:** asserts_passed=253 · FAIL=0

??? note "Full run log (`S998_run.log`)"
    ```text
    --8<-- "src/probe/S998_run.log"
    ```


---

[← all nodes](index.md)
