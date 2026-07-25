# Graph map

The chain of [Tome I](tome1-prime.md) has a machine representation: a **graph of nodes and
dependencies**, where each node carries `id · kind · status · address · ancestors`. An edge is an
**ancestor**, not a citation: agreement between two descendants of a common ancestor counts as
**one** witness, not two.

**[Open the node index →](theorems/index.md)**

*(the graph is published as node pages rather than as a single dump: one page per node, and every
dependency is a link between pages, so the structure is walked rather than viewed)*

---

## How to read it

| node field | meaning |
|:--|:--|
| `kind` | input · theorem · object · identity · seam · conclusion · fence |
| `status` | `axiom` (declared premise) · `measured` (probe + verification) · `candidate` (a cast fitted to data) · `measured-negative` (**proved absence**) |
| `address` | where the primary act lies: section · probe `S…` · verification `J-…` · seal |
| `deps` | ancestors — witness multiplicity is counted along these |

**Read `measured-negative` carefully:** it does not say "it did not work out", it says the object was
**proved not to exist**. Such a node carries as much weight as a positive one.

---

## A boundary of the map, named explicitly

The graph covers the **prime floor**. The most recent layer of the chain (the terminal step of the
derived series and everything after it) is **not yet represented in the graph** — those links are
cited in Tome I **directly by the seals of the rulings**. So a clean state of the graph checker
applies to the covered part, not to the whole exposition. This is said here so that a "one green
tick" reading is not taken more widely than it holds.
