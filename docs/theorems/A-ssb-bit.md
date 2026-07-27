# Stage 2: the free ±m₀ choice is the single free bit (T36-bit)

**Kind:** theorem · **Status:** candidate · **Address:** §3.seam · MIRROR_ASSEMBLY_LAW.md · S1035 · J-0482..0492

*A node of the prime graph. Status is carried by the graph and is not restated here.*

## INPUT

| premise | kind | status |
|:--|:--|:--|
| [Stabilization intrinsic: κ=Λ implies m₀>0 ⟺ Λ<a](A-stab.md) | theorem | candidate |
| [Dihedral hinge D_h=⟨c,w₀⟩ of the −2↔−1 joint](OBJ-Dh.md) | object | measured |
| [The arrow is a single bit of realization (5 legs + wedge witness ver:2): no A/B discriminator exists](T36.md) | theorem | candidate |

## FUNCTION

The free choice of sign in ±m₀ is the single free bit (the T36-bit), given by the sign-character of D_h (c→+1, w₀→−1, unique for all n, with the realizer σ_z carrying values (1,−1)); it is a canonical map determined by ancestor, not by role.

> **Boundary:** Restricted to the −2↔−1 joint; the spontaneous ±m₀ choice is the free ℤ/2 sign-character of D_h.

## OUTPUT

Supplies the single free bit consumed downstream by the stabilization result (A-stab).

**Consumed by:** [Crown result: spontaneity ⟺ non-derivability of sign(m₀)](A-nonderiv.md) · [Registry wrapper (Section 12) over the A-* arc, not a separate derivation: mirror-assembly law (two-stage freezing + kappa=Lambda stabilization)](T38.md)

### 8.3 Notation — short carriers used above

| short | expansion | why the short name is dangerous |
|:--|:--|:--|
| `m_0` | the quantity of the arc whose sign is the subject; its value is a representative of the regime, not canon, and is not derived here | the letter m denotes something else on another floor |

## Measurements

**S1035**

??? note "Probe code (`S1035.py`)"
    ```python
    --8<-- "src/probe/S1035.py"
    ```

**Result:** exit=0

??? note "Full run log (`S1035_run.log`)"
    ```text
    --8<-- "src/probe/S1035_run.log"
    ```


---

[← all nodes](index.md)
