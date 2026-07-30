# Arrow/class-count (T28+T36 from W41): Coxeter action c leads to epsilon -> -epsilon, 1 bit, D_h orbit

**Kind:** seam · **Status:** measured · **Address:** §3.seam.6 · S1004/S1016/S1017 · J-0460/J-0471/J-0472 · rulings §14/§15

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Count of ie classes: C(q)=2^q, Ctilde=0/1/2/4; q=1 gives a unique class up to the arrow](T28.md) | theorem | candidate |
| [The arrow is a single bit of realization (5 legs + wedge witness ver:2): no A/B discriminator exists](T36.md) | theorem | candidate |
| [W41 Coxeter-time axiom: time as the action of <c> of order h_cox, won for all d; spectral monism refuted (3 closed forms)](W41.md) | theorem | measured |

## FUNCTION

Traces how the Coxeter action c, feeding from W41 through T28 and T36, produces the epsilon -> -epsilon flip, a single bit, and the D_h orbit structure.

> **Boundary:** Levels -2 to -1.

## OUTPUT

Feeds the epsilon-flip identity id-2.2 and the orbit identity id-2.5.

**Consumed by:** nothing in this graph — a terminal node.

## FACTORIZATION

A part of [Dihedral hinge D_h=⟨c,w₀⟩ of the −2↔−1 joint](OBJ-Dh.md) — component w0 (measured).

## Measurements

**S1004**

??? note "Probe code (`S1004.py`)"
    ```python
    --8<-- "src/probe/S1004.py"
    ```

**Result:** asserts_passed=30 · FAIL=0

??? note "Full run log (`S1004_run.log`)"
    ```text
    --8<-- "src/probe/S1004_run.log"
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

**S1017**

??? note "Probe code (`S1017.py`)"
    ```python
    --8<-- "src/probe/S1017.py"
    ```

**Result:** asserts_passed=28 · FAIL=0

??? note "Full run log (`S1017_run.log`)"
    ```text
    --8<-- "src/probe/S1017_run.log"
    ```


---

[← all nodes](index.md)
