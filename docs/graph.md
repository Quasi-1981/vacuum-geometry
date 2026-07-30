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
| `status` | `axiom` (declared premise) · `measured` (probe + verification) · `candidate` (a cast fitted to data) · `measured-negative` (**proved absence**) · `named` (identified, not yet derived) · `author-word` (the author's own addition, carried as data) |
| `address` | where the primary act lies: section · probe `S…` · verification `J-…` · seal |
| `deps` | ancestors — witness multiplicity is counted along these |

**Read `measured-negative` carefully:** it does not say "it did not work out", it says the object was
**proved not to exist**. Such a node carries as much weight as a positive one.

**`∀d-symbolic` is not a seventh status:** it is a qualifier used in verdict prose to say a
`measured` node was derived for symbolic $d$ rather than enumerated; no node carries it in the
status field.

**A `measured` node may stand on `candidate` ancestors, and this is not a contradiction.** The
status grades *the act of measuring that node* — was its own claim probed and independently
verified — and it does not grade the truth of what the node stands on. So `measured` reads as
"this step was carried out and checked", not as "everything below it is settled". A chain is
therefore only as load-bearing as its **weakest** link: to judge how much a result carries,
follow `deps` and read the statuses along the way, rather than reading the status of the last
node alone. The graph deliberately keeps the two questions apart instead of collapsing them
into one number — a node whose own act is clean does not inherit certainty from below, and
does not lose it either.

---

## The map

<!--graph-map-start-->

<div class="gm-wrap">
<div class="gm-controls">
<input class="gm-search" type="text" placeholder="Search a node…" aria-label="Search a node">
</div>
<div class="gm-legend"><span class="gm-legend-item"><svg width="14" height="14" aria-hidden="true"><rect x="1" y="1" width="12" height="12" rx="3" fill="#3b6fb5" stroke="#20406e" stroke-width="1.2"/></svg> axiom</span> &nbsp; <span class="gm-legend-item"><svg width="14" height="14" aria-hidden="true"><rect x="1" y="1" width="12" height="12" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2"/></svg> measured</span> &nbsp; <span class="gm-legend-item"><svg width="14" height="14" aria-hidden="true"><rect x="1" y="1" width="12" height="12" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2"/></svg> candidate</span> &nbsp; <span class="gm-legend-item"><svg width="14" height="14" aria-hidden="true"><rect x="1" y="1" width="12" height="12" rx="3" fill="#b8433f" stroke="#7a2b28" stroke-width="1.2"/></svg> measured-negative</span> &nbsp; <span class="gm-legend-item"><svg width="14" height="14" aria-hidden="true"><rect x="1" y="1" width="12" height="12" rx="3" fill="#8a8a8a" stroke="#5a5a5a" stroke-width="1.2"/></svg> named</span> &nbsp; <span class="gm-legend-item"><svg width="14" height="14" aria-hidden="true"><rect x="1" y="1" width="12" height="12" rx="3" fill="#7a4fb0" stroke="#4f3373" stroke-width="1.2"/></svg> author-word</span></div>
<svg class="gm-svg" viewBox="0 0 2110 1502" width="2110" height="1502" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Prime graph map">
<g class="gm-edges">
<line class="gm-edge" x1="770.0" y1="209.0" x2="550.0" y2="305.0" />
<line class="gm-edge" x1="250.0" y1="65.0" x2="550.0" y2="329.0" />
<line class="gm-edge" x1="250.0" y1="113.0" x2="550.0" y2="329.0" />
<line class="gm-edge" x1="250.0" y1="89.0" x2="550.0" y2="353.0" />
<line class="gm-edge" x1="250.0" y1="89.0" x2="550.0" y2="617.0" />
<line class="gm-edge" x1="250.0" y1="89.0" x2="550.0" y2="1049.0" />
<line class="gm-edge" x1="250.0" y1="89.0" x2="550.0" y2="1313.0" />
<line class="gm-edge" x1="770.0" y1="1313.0" x2="550.0" y2="1337.0" />
<line class="gm-edge" x1="770.0" y1="305.0" x2="550.0" y2="1337.0" />
<line class="gm-edge" x1="770.0" y1="329.0" x2="550.0" y2="1337.0" />
<line class="gm-edge" x1="250.0" y1="161.0" x2="550.0" y2="1337.0" />
<line class="gm-edge" x1="770.0" y1="1313.0" x2="550.0" y2="1361.0" />
<line class="gm-edge" x1="770.0" y1="1337.0" x2="550.0" y2="1385.0" />
<line class="gm-edge" x1="250.0" y1="89.0" x2="550.0" y2="1409.0" />
<line class="gm-edge" x1="250.0" y1="89.0" x2="550.0" y2="1433.0" />
<line class="gm-edge" x1="250.0" y1="89.0" x2="550.0" y2="377.0" />
<line class="gm-edge" x1="770.0" y1="377.0" x2="550.0" y2="401.0" />
<line class="gm-edge" x1="770.0" y1="1337.0" x2="1590.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="1385.0" x2="1590.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="305.0" x2="290.0" y2="137.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="290.0" y2="137.0" />
<line class="gm-edge" x1="510.0" y1="137.0" x2="290.0" y2="161.0" />
<line class="gm-edge" x1="250.0" y1="113.0" x2="290.0" y2="161.0" />
<line class="gm-edge" x1="510.0" y1="161.0" x2="290.0" y2="233.0" />
<line class="gm-edge" x1="510.0" y1="137.0" x2="290.0" y2="113.0" />
<line class="gm-edge" x1="510.0" y1="161.0" x2="290.0" y2="113.0" />
<line class="gm-edge" x1="510.0" y1="113.0" x2="290.0" y2="89.0" />
<line class="gm-edge" x1="510.0" y1="113.0" x2="290.0" y2="257.0" />
<line class="gm-edge" x1="510.0" y1="233.0" x2="290.0" y2="257.0" />
<line class="gm-edge" x1="510.0" y1="137.0" x2="290.0" y2="209.0" />
<line class="gm-edge" x1="510.0" y1="137.0" x2="290.0" y2="185.0" />
<line class="gm-edge" x1="250.0" y1="89.0" x2="550.0" y2="425.0" />
<line class="gm-edge" x1="770.0" y1="425.0" x2="550.0" y2="449.0" />
<line class="gm-edge" x1="770.0" y1="425.0" x2="550.0" y2="473.0" />
<line class="gm-edge" x1="770.0" y1="425.0" x2="550.0" y2="497.0" />
<line class="gm-edge" x1="770.0" y1="425.0" x2="550.0" y2="521.0" />
<line class="gm-edge" x1="770.0" y1="449.0" x2="550.0" y2="545.0" />
<line class="gm-edge" x1="770.0" y1="473.0" x2="550.0" y2="545.0" />
<line class="gm-edge" x1="770.0" y1="521.0" x2="550.0" y2="545.0" />
<line class="gm-edge" x1="770.0" y1="473.0" x2="550.0" y2="569.0" />
<line class="gm-edge" x1="770.0" y1="521.0" x2="550.0" y2="569.0" />
<line class="gm-edge" x1="770.0" y1="305.0" x2="550.0" y2="593.0" />
<line class="gm-edge" x1="770.0" y1="569.0" x2="550.0" y2="641.0" />
<line class="gm-edge" x1="770.0" y1="521.0" x2="550.0" y2="665.0" />
<line class="gm-edge" x1="770.0" y1="641.0" x2="550.0" y2="665.0" />
<line class="gm-edge" x1="770.0" y1="665.0" x2="550.0" y2="689.0" />
<line class="gm-edge" x1="250.0" y1="137.0" x2="550.0" y2="689.0" />
<line class="gm-edge" x1="770.0" y1="569.0" x2="550.0" y2="713.0" />
<line class="gm-edge" x1="770.0" y1="737.0" x2="550.0" y2="713.0" />
<line class="gm-edge" x1="770.0" y1="641.0" x2="550.0" y2="737.0" />
<line class="gm-edge" x1="770.0" y1="737.0" x2="550.0" y2="761.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="785.0" />
<line class="gm-edge" x1="770.0" y1="809.0" x2="550.0" y2="785.0" />
<line class="gm-edge" x1="770.0" y1="833.0" x2="550.0" y2="785.0" />
<line class="gm-edge" x1="770.0" y1="857.0" x2="550.0" y2="785.0" />
<line class="gm-edge" x1="770.0" y1="881.0" x2="550.0" y2="785.0" />
<line class="gm-edge" x1="770.0" y1="905.0" x2="550.0" y2="785.0" />
<line class="gm-edge" x1="770.0" y1="929.0" x2="550.0" y2="785.0" />
<line class="gm-edge" x1="770.0" y1="953.0" x2="550.0" y2="785.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="809.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="833.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="857.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="881.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="905.0" />
<line class="gm-edge" x1="250.0" y1="113.0" x2="550.0" y2="905.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="929.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="953.0" />
<line class="gm-edge" x1="770.0" y1="785.0" x2="1330.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="905.0" x2="1330.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="785.0" x2="550.0" y2="1457.0" />
<line class="gm-edge" x1="770.0" y1="953.0" x2="550.0" y2="1457.0" />
<line class="gm-edge" x1="770.0" y1="1457.0" x2="1850.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="785.0" x2="1850.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="1433.0" x2="550.0" y2="977.0" />
<line class="gm-edge" x1="770.0" y1="977.0" x2="550.0" y2="1001.0" />
<line class="gm-edge" x1="770.0" y1="977.0" x2="550.0" y2="1025.0" />
<line class="gm-edge" x1="770.0" y1="1025.0" x2="550.0" y2="1073.0" />
<line class="gm-edge" x1="770.0" y1="1433.0" x2="550.0" y2="1073.0" />
<line class="gm-edge" x1="770.0" y1="881.0" x2="550.0" y2="1097.0" />
<line class="gm-edge" x1="770.0" y1="905.0" x2="550.0" y2="1097.0" />
<line class="gm-edge" x1="770.0" y1="953.0" x2="550.0" y2="1097.0" />
<line class="gm-edge" x1="770.0" y1="1457.0" x2="550.0" y2="1097.0" />
<line class="gm-edge" x1="770.0" y1="905.0" x2="550.0" y2="1121.0" />
<line class="gm-edge" x1="770.0" y1="1073.0" x2="550.0" y2="1121.0" />
<line class="gm-edge" x1="770.0" y1="977.0" x2="550.0" y2="1121.0" />
<line class="gm-edge" x1="770.0" y1="593.0" x2="550.0" y2="1121.0" />
<line class="gm-edge" x1="770.0" y1="1121.0" x2="550.0" y2="1145.0" />
<line class="gm-edge" x1="770.0" y1="857.0" x2="550.0" y2="1145.0" />
<line class="gm-edge" x1="770.0" y1="881.0" x2="550.0" y2="1145.0" />
<line class="gm-edge" x1="770.0" y1="953.0" x2="550.0" y2="1169.0" />
<line class="gm-edge" x1="770.0" y1="1025.0" x2="550.0" y2="1169.0" />
<line class="gm-edge" x1="770.0" y1="1097.0" x2="550.0" y2="1169.0" />
<line class="gm-edge" x1="770.0" y1="857.0" x2="550.0" y2="1193.0" />
<line class="gm-edge" x1="770.0" y1="977.0" x2="550.0" y2="1193.0" />
<line class="gm-edge" x1="770.0" y1="1001.0" x2="550.0" y2="1217.0" />
<line class="gm-edge" x1="770.0" y1="1145.0" x2="550.0" y2="1217.0" />
<line class="gm-edge" x1="770.0" y1="1169.0" x2="550.0" y2="1217.0" />
<line class="gm-edge" x1="770.0" y1="737.0" x2="550.0" y2="1217.0" />
<line class="gm-edge" x1="770.0" y1="689.0" x2="550.0" y2="1217.0" />
<line class="gm-edge" x1="770.0" y1="1169.0" x2="550.0" y2="1241.0" />
<line class="gm-edge" x1="770.0" y1="1217.0" x2="550.0" y2="1241.0" />
<line class="gm-edge" x1="770.0" y1="905.0" x2="810.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="953.0" x2="810.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="1121.0" x2="810.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="1169.0" x2="810.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="1001.0" x2="810.0" y2="89.0" />
<line class="gm-edge" x1="770.0" y1="1145.0" x2="810.0" y2="113.0" />
<line class="gm-edge" x1="770.0" y1="1145.0" x2="810.0" y2="137.0" />
<line class="gm-edge" x1="770.0" y1="1169.0" x2="810.0" y2="161.0" />
<line class="gm-edge" x1="770.0" y1="1217.0" x2="810.0" y2="185.0" />
<line class="gm-edge" x1="770.0" y1="1193.0" x2="810.0" y2="209.0" />
<line class="gm-edge" x1="770.0" y1="1025.0" x2="1070.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="1121.0" x2="1070.0" y2="89.0" />
<line class="gm-edge" x1="770.0" y1="905.0" x2="1070.0" y2="89.0" />
<line class="gm-edge" x1="770.0" y1="1121.0" x2="1070.0" y2="113.0" />
<line class="gm-edge" x1="770.0" y1="1169.0" x2="1070.0" y2="113.0" />
<line class="gm-edge" x1="770.0" y1="953.0" x2="1070.0" y2="113.0" />
<line class="gm-edge" x1="770.0" y1="1145.0" x2="1070.0" y2="137.0" />
<line class="gm-edge" x1="770.0" y1="905.0" x2="1070.0" y2="137.0" />
<line class="gm-edge" x1="770.0" y1="1193.0" x2="1070.0" y2="161.0" />
<line class="gm-edge" x1="770.0" y1="857.0" x2="1070.0" y2="161.0" />
<line class="gm-edge" x1="770.0" y1="1457.0" x2="1070.0" y2="185.0" />
<line class="gm-edge" x1="770.0" y1="1001.0" x2="1070.0" y2="185.0" />
<line class="gm-edge" x1="770.0" y1="1217.0" x2="1070.0" y2="185.0" />
<line class="gm-edge" x1="770.0" y1="1241.0" x2="1070.0" y2="209.0" />
<line class="gm-edge" x1="770.0" y1="953.0" x2="290.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="1001.0" x2="290.0" y2="65.0" />
<line class="gm-edge" x1="770.0" y1="1025.0" x2="550.0" y2="89.0" />
<line class="gm-edge" x1="770.0" y1="977.0" x2="550.0" y2="89.0" />
<line class="gm-edge" x1="770.0" y1="1145.0" x2="550.0" y2="89.0" />
<line class="gm-edge-neg" x1="770.0" y1="1145.0" x2="550.0" y2="113.0" />
<line class="gm-edge-neg" x1="770.0" y1="89.0" x2="550.0" y2="113.0" />
<line class="gm-edge" x1="770.0" y1="809.0" x2="550.0" y2="257.0" />
<line class="gm-edge" x1="770.0" y1="833.0" x2="550.0" y2="257.0" />
<line class="gm-edge" x1="770.0" y1="113.0" x2="550.0" y2="257.0" />
<line class="gm-edge" x1="250.0" y1="137.0" x2="550.0" y2="257.0" />
<line class="gm-edge" x1="770.0" y1="257.0" x2="550.0" y2="233.0" />
<line class="gm-edge" x1="770.0" y1="1217.0" x2="550.0" y2="233.0" />
<line class="gm-edge" x1="510.0" y1="65.0" x2="550.0" y2="233.0" />
<line class="gm-edge" x1="770.0" y1="233.0" x2="550.0" y2="185.0" />
<line class="gm-edge" x1="770.0" y1="1121.0" x2="550.0" y2="185.0" />
<line class="gm-edge" x1="770.0" y1="185.0" x2="550.0" y2="137.0" />
<line class="gm-edge" x1="770.0" y1="1241.0" x2="550.0" y2="137.0" />
<line class="gm-edge" x1="770.0" y1="1217.0" x2="550.0" y2="137.0" />
<line class="gm-edge" x1="770.0" y1="89.0" x2="550.0" y2="1265.0" />
<line class="gm-edge" x1="770.0" y1="233.0" x2="550.0" y2="1265.0" />
<line class="gm-edge" x1="770.0" y1="257.0" x2="550.0" y2="1265.0" />
<line class="gm-edge" x1="770.0" y1="185.0" x2="550.0" y2="1289.0" />
<line class="gm-edge" x1="770.0" y1="137.0" x2="550.0" y2="1289.0" />
<line class="gm-edge" x1="250.0" y1="65.0" x2="290.0" y2="281.0" />
<line class="gm-edge" x1="510.0" y1="281.0" x2="550.0" y2="209.0" />
<line class="gm-edge" x1="250.0" y1="137.0" x2="550.0" y2="209.0" />
<line class="gm-edge-neg" x1="510.0" y1="281.0" x2="550.0" y2="281.0" />
<line class="gm-edge" x1="510.0" y1="281.0" x2="550.0" y2="161.0" />
<line class="gm-edge" x1="770.0" y1="209.0" x2="550.0" y2="161.0" />
<line class="gm-edge" x1="770.0" y1="161.0" x2="550.0" y2="65.0" />
</g>
<g class="gm-colheads">
<text class="gm-colhead" x="30.0" y="44">input</text>
<text class="gm-colhead" x="290.0" y="44">object</text>
<text class="gm-colhead" x="550.0" y="44">theorem</text>
<text class="gm-colhead" x="810.0" y="44">identity</text>
<text class="gm-colhead" x="1070.0" y="44">seam</text>
<text class="gm-colhead" x="1330.0" y="44">definition</text>
<text class="gm-colhead" x="1590.0" y="44">conclusion</text>
<text class="gm-colhead" x="1850.0" y="44">fence</text>
</g>
<g class="gm-nodes">
<a href="../theorems/AX-alphabet/" class="gm-node" data-id="AX-alphabet" data-label="alphabet axiom: bare set of d+1 with s_{d+1} democracy (root swap)">
<rect x="30.0" y="56.0" width="18" height="18" rx="3" fill="#3b6fb5" stroke="#20406e" stroke-width="1.2" />
<title>AX-alphabet — Alphabet axiom: bare set of d+1 with S_{d+1} democracy (root swap)</title>
<text class="gm-label" x="52.0" y="69.0">Alphabet axiom: bare set of d…</text>
</a>
<a href="../theorems/AX-closure/" class="gm-node" data-id="AX-closure" data-label="closure axiom on real quadratic spaces forces (3,0)/(3,1)">
<rect x="30.0" y="80.0" width="18" height="18" rx="3" fill="#3b6fb5" stroke="#20406e" stroke-width="1.2" />
<title>AX-closure — Closure axiom on real quadratic spaces forces (3,0)/(3,1)</title>
<text class="gm-label" x="52.0" y="93.0">Closure axiom on real quadrat…</text>
</a>
<a href="../theorems/AX-dimer/" class="gm-node" data-id="AX-dimer" data-label="dimer axiom: the time-bond marked axis is 1 bit of 'time exists'">
<rect x="30.0" y="104.0" width="18" height="18" rx="3" fill="#3b6fb5" stroke="#20406e" stroke-width="1.2" />
<title>AX-dimer — Dimer axiom: the time-bond marked axis is 1 bit of 'time exists'</title>
<text class="gm-label" x="52.0" y="117.0">Dimer axiom: the time-bond ma…</text>
</a>
<a href="../theorems/AX-lambda/" class="gm-node" data-id="AX-lambda" data-label="λ-scale: the sole dimensional free handle">
<rect x="30.0" y="128.0" width="18" height="18" rx="3" fill="#3b6fb5" stroke="#20406e" stroke-width="1.2" />
<title>AX-lambda — Λ-scale: the sole dimensional free handle</title>
<text class="gm-label" x="52.0" y="141.0">Λ-scale: the sole dimensional…</text>
</a>
<a href="../theorems/IMP-d3/" class="gm-node" data-id="IMP-d3" data-label="import d≥3: spatial dimensionality is a declared boundary, not derived">
<rect x="30.0" y="152.0" width="18" height="18" rx="3" fill="#8a8a8a" stroke="#5a5a5a" stroke-width="1.2" />
<title>IMP-d3 — IMPORT d≥3: spatial dimensionality is a declared boundary, not derived</title>
<text class="gm-label" x="52.0" y="165.0">IMPORT d≥3: spatial dimension…</text>
</a>
<a href="../theorems/OBJ-Dh/" class="gm-node" data-id="OBJ-Dh" data-label="dihedral hinge d_h=⟨c,w₀⟩ of the −2↔−1 joint">
<rect x="290.0" y="56.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-Dh — Dihedral hinge D_h=⟨c,w₀⟩ of the −2↔−1 joint</title>
<text class="gm-label" x="312.0" y="69.0">Dihedral hinge D_h=⟨c,w₀⟩ of …</text>
</a>
<a href="../theorems/OBJ-H/" class="gm-node" data-id="OBJ-H" data-label="two-component matrix h with invariants {i,h}">
<rect x="290.0" y="80.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-H — Two-component matrix H with invariants {I,H}</title>
<text class="gm-label" x="312.0" y="93.0">Two-component matrix H with i…</text>
</a>
<a href="../theorems/OBJ-box/" class="gm-node" data-id="OBJ-box" data-label="scalar box: λ(ψ,ν)=t_a(ψ)−t_col(ν)">
<rect x="290.0" y="104.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-box — Scalar Box: Λ(ψ,ν)=T_A(ψ)−T_col(ν)</title>
<text class="gm-label" x="312.0" y="117.0">Scalar Box: Λ(ψ,ν)=T_A(ψ)−T_c…</text>
</a>
<a href="../theorems/OBJ-cell/" class="gm-node" data-id="OBJ-cell" data-label="a_d cell as weights of the fundamental su(d+1) representation">
<rect x="290.0" y="128.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-cell — A_d cell as weights of the fundamental su(d+1) representation</title>
<text class="gm-label" x="312.0" y="141.0">A_d cell as weights of the fu…</text>
</a>
<a href="../theorems/OBJ-column/" class="gm-node" data-id="OBJ-column" data-label="column/dimer: time-bond u₀ with translation symmetry p_col=d+1">
<rect x="290.0" y="152.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-column — Column/dimer: time-bond u₀ with translation symmetry P_col=d+1</title>
<text class="gm-label" x="312.0" y="165.0">Column/dimer: time-bond u₀ wi…</text>
</a>
<a href="../theorems/OBJ-commutant/" class="gm-node" data-id="OBJ-commutant" data-label="commutant ladder with block ranks (1,m−1,d−m)">
<rect x="290.0" y="176.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-commutant — Commutant ladder with block ranks (1,m−1,d−m)</title>
<text class="gm-label" x="312.0" y="189.0">Commutant ladder with block r…</text>
</a>
<a href="../theorems/OBJ-coxeter/" class="gm-node" data-id="OBJ-coxeter" data-label="coxeter ladder: cartan(a_d) eigenvalues 2−2cos(πj/h_cox)">
<rect x="290.0" y="200.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-coxeter — Coxeter ladder: Cartan(A_d) eigenvalues 2−2cos(πj/h_cox)</title>
<text class="gm-label" x="312.0" y="213.0">Coxeter ladder: Cartan(A_d) e…</text>
</a>
<a href="../theorems/OBJ-dual/" class="gm-node" data-id="OBJ-dual" data-label="dual circle ℤ/p_col, pontryagin dual of the column translation">
<rect x="290.0" y="224.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-dual — Dual circle ℤ/P_col, Pontryagin dual of the column translation</title>
<text class="gm-label" x="312.0" y="237.0">Dual circle ℤ/P_col, Pontryag…</text>
</a>
<a href="../theorems/OBJ-nodal/" class="gm-node" data-id="OBJ-nodal" data-label="nodal set f=0 corresponds to d characters of ℤ/p_col">
<rect x="290.0" y="248.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-nodal — Nodal set f=0 corresponds to d characters of ℤ/P_col</title>
<text class="gm-label" x="312.0" y="261.0">Nodal set f=0 corresponds to …</text>
</a>
<a href="../theorems/OBJ-sln/" class="gm-node" data-id="OBJ-sln" data-label="simple algebra sl(n), n=d+1, with unique invariant form κ">
<rect x="290.0" y="272.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>OBJ-sln — Simple algebra sl(n), n=d+1, with unique invariant form κ</title>
<text class="gm-label" x="312.0" y="285.0">Simple algebra sl(n), n=d+1, …</text>
</a>
<a href="../theorems/A-axis/" class="gm-node" data-id="A-axis" data-label="κ: minimal parabolic p_α=(d,1), q=1, variant b; time=so(2)_α">
<rect x="550.0" y="56.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>A-axis — κ: minimal parabolic p_α=(d,1), q=1, variant B; time=so(2)_α</title>
<text class="gm-label" x="572.0" y="69.0">κ: minimal parabolic p_α=(d,1…</text>
</a>
<a href="../theorems/A-collapse-pre/" class="gm-node" data-id="A-collapse-pre" data-label="stage 1: t29 collapse (ε-even, zero separation) as precondition">
<rect x="550.0" y="80.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>A-collapse-pre — Stage 1: T29 collapse (ε-even, zero separation) as precondition</title>
<text class="gm-label" x="572.0" y="93.0">Stage 1: T29 collapse (ε-even…</text>
</a>
<a href="../theorems/A-instability/" class="gm-node" data-id="A-instability" data-label="a&gt;0 ⟹ m=0 unstable; monotonic trace forces runaway, no m₀² term">
<rect x="550.0" y="104.0" width="18" height="18" rx="3" fill="#b8433f" stroke="#7a2b28" stroke-width="1.2" />
<title>A-instability — a&gt;0 ⟹ m=0 unstable; monotonic trace forces runaway, no m₀² term</title>
<text class="gm-label" x="572.0" y="117.0">a&gt;0 ⟹ m=0 unstable; monotonic…</text>
</a>
<a href="../theorems/A-inherit/" class="gm-node" data-id="A-inherit" data-label="inheritance: arrow spontaneous exactly once per connected lattice">
<rect x="550.0" y="128.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>A-inherit — Inheritance: arrow spontaneous exactly once per connected lattice</title>
<text class="gm-label" x="572.0" y="141.0">Inheritance: arrow spontaneou…</text>
</a>
<a href="../theorems/A-kappa/" class="gm-node" data-id="A-kappa" data-label="κ unique: jordan-chevalley split — space (h_cart) ⊥ axis (n±)">
<rect x="550.0" y="152.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>A-kappa — κ unique: Jordan-Chevalley split — space (h_cart) ⊥ axis (n±)</title>
<text class="gm-label" x="572.0" y="165.0">κ unique: Jordan-Chevalley sp…</text>
</a>
<a href="../theorems/A-nonderiv/" class="gm-node" data-id="A-nonderiv" data-label="crown result: spontaneity ⟺ non-derivability of sign(m₀)">
<rect x="550.0" y="176.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>A-nonderiv — Crown result: spontaneity ⟺ non-derivability of sign(m₀)</title>
<text class="gm-label" x="572.0" y="189.0">Crown result: spontaneity ⟺ n…</text>
</a>
<a href="../theorems/A-space/" class="gm-node" data-id="A-space" data-label="space = cartan algebra a_d via schur-κ, for all d">
<rect x="550.0" y="200.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>A-space — Space = Cartan algebra A_d via Schur-κ, for all d</title>
<text class="gm-label" x="572.0" y="213.0">Space = Cartan algebra A_d vi…</text>
</a>
<a href="../theorems/A-ssb-bit/" class="gm-node" data-id="A-ssb-bit" data-label="stage 2: the free ±m₀ choice is the single free bit (t36-bit)">
<rect x="550.0" y="224.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>A-ssb-bit — Stage 2: the free ±m₀ choice is the single free bit (T36-bit)</title>
<text class="gm-label" x="572.0" y="237.0">Stage 2: the free ±m₀ choice …</text>
</a>
<a href="../theorems/A-stab/" class="gm-node" data-id="A-stab" data-label="stabilization intrinsic: κ=λ implies m₀&gt;0 ⟺ λ&lt;a">
<rect x="550.0" y="248.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>A-stab — Stabilization intrinsic: κ=Λ implies m₀&gt;0 ⟺ Λ&lt;a</title>
<text class="gm-label" x="572.0" y="261.0">Stabilization intrinsic: κ=Λ …</text>
</a>
<a href="../theorems/A-time-neg/" class="gm-node" data-id="A-time-neg" data-label="time not derivable from order: weyl channel gives only an arrow">
<rect x="550.0" y="272.0" width="18" height="18" rx="3" fill="#b8433f" stroke="#7a2b28" stroke-width="1.2" />
<title>A-time-neg — Time not derivable from order: Weyl channel gives only an arrow</title>
<text class="gm-label" x="572.0" y="285.0">Time not derivable from order…</text>
</a>
<a href="../theorems/AX-cell/" class="gm-node" data-id="AX-cell" data-label="a_d cell derived from alphabet (∀d, three independent routes)">
<rect x="550.0" y="296.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>AX-cell — A_d cell derived from alphabet (∀d, three independent routes)</title>
<text class="gm-label" x="572.0" y="309.0">A_d cell derived from alphabe…</text>
</a>
<a href="../theorems/AX-indef/" class="gm-node" data-id="AX-indef" data-label="indefiniteness (q≥1) demoted from input to theorem">
<rect x="550.0" y="320.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>AX-indef — Indefiniteness (q≥1) demoted from input to theorem</title>
<text class="gm-label" x="572.0" y="333.0">Indefiniteness (q≥1) demoted …</text>
</a>
<a href="../theorems/T1/" class="gm-node" data-id="T1" data-label="closure preserves definiteness: the carrier is definite exactly when the start's carrier is (r-b)">
<rect x="550.0" y="344.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T1 — Closure preserves definiteness: the carrier is definite exactly when the start's carrier is (R-B)</title>
<text class="gm-label" x="572.0" y="357.0">Closure preserves definitenes…</text>
</a>
<a href="../theorems/T10/" class="gm-node" data-id="T10" data-label="orientation identity pf(sᵀωs)=det(s)·pf(ω)">
<rect x="550.0" y="368.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T10 — Orientation identity Pf(SᵀΩS)=det(S)·Pf(Ω)</title>
<text class="gm-label" x="572.0" y="381.0">Orientation identity Pf(SᵀΩS)…</text>
</a>
<a href="../theorems/T11/" class="gm-node" data-id="T11" data-label="aut-invariance of the preorder (290/290 on a carved grid); no swap on (3,1) when p!=q">
<rect x="550.0" y="392.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T11 — Aut-invariance of the preorder (290/290 on a carved grid); no swap on (3,1) when p!=q</title>
<text class="gm-label" x="572.0" y="405.0">Aut-invariance of the preorde…</text>
</a>
<a href="../theorems/T12/" class="gm-node" data-id="T12" data-label="the carrier is three-layered (combinatorial closure equals lie closure)">
<rect x="550.0" y="416.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T12 — The carrier is three-layered (combinatorial closure equals Lie closure)</title>
<text class="gm-label" x="572.0" y="429.0">The carrier is three-layered …</text>
</a>
<a href="../theorems/T13/" class="gm-node" data-id="T13" data-label="the 'channel law' form is refuted; leg-1 measured across families, dim c=sum m^2+k(k-1)/2">
<rect x="550.0" y="440.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T13 — The 'channel law' form is refuted; leg-1 measured across families, dim c=Sum m^2+k(k-1)/2</title>
<text class="gm-label" x="572.0" y="453.0">The 'channel law' form is ref…</text>
</a>
<a href="../theorems/T14/" class="gm-node" data-id="T14" data-label="taxonomy of nilpotents; deep iff witt index&gt;=2; (3,1) immune beyond the imported witt boundary">
<rect x="550.0" y="464.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T14 — Taxonomy of nilpotents; deep iff Witt index&gt;=2; (3,1) immune beyond the imported Witt boundary</title>
<text class="gm-label" x="572.0" y="477.0">Taxonomy of nilpotents; deep …</text>
</a>
<a href="../theorems/T15/" class="gm-node" data-id="T15" data-label="sector law (orientability): multiplicity 1 gives two cones; multiplicity&gt;=2 non-orientable">
<rect x="550.0" y="488.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T15 — Sector law (orientability): multiplicity 1 gives two cones; multiplicity&gt;=2 non-orientable</title>
<text class="gm-label" x="572.0" y="501.0">Sector law (orientability): m…</text>
</a>
<a href="../theorems/T16/" class="gm-node" data-id="T16" data-label="shift: brackets forced within so(p,q) semidirect r^n; construct-class lineage measured (s943)">
<rect x="550.0" y="512.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T16 — Shift: brackets forced within so(p,q) semidirect R^n; construct-class lineage measured (S943)</title>
<text class="gm-label" x="572.0" y="525.0">Shift: brackets forced within…</text>
</a>
<a href="../theorems/T17/" class="gm-node" data-id="T17" data-label="the map of nonabelian kernels is closed: 4 classes; radical dichotomy abelian vs heisenberg">
<rect x="550.0" y="536.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T17 — The map of nonabelian kernels is closed: 4 classes; radical dichotomy abelian vs Heisenberg</title>
<text class="gm-label" x="572.0" y="549.0">The map of nonabelian kernels…</text>
</a>
<a href="../theorems/T18/" class="gm-node" data-id="T18" data-label="universal mechanism derived for all n: c=(sp2+so(eta|g)) semidirect h_heis; center from [module,module]">
<rect x="550.0" y="560.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T18 — Universal mechanism derived for all n: c=(sp2+so(eta|G)) semidirect h_heis; center from [module,module]</title>
<text class="gm-label" x="572.0" y="573.0">Universal mechanism derived f…</text>
</a>
<a href="../theorems/T19/" class="gm-node" data-id="T19" data-label="cell to a_d native map (gram proportional to cartan, z=d+1, d nodes); symbolic for all d (ver:2)">
<rect x="550.0" y="584.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T19 — Cell to A_d native map (Gram proportional to Cartan, z=d+1, d nodes); symbolic for all d (ver:2)</title>
<text class="gm-label" x="572.0" y="597.0">Cell to A_d native map (Gram …</text>
</a>
<a href="../theorems/T2/" class="gm-node" data-id="T2" data-label="direction of order: the j-sector is bracket-closed; the full k forces so(p,q)">
<rect x="550.0" y="608.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T2 — Direction of order: the J-sector is bracket-closed; the full K forces so(p,q)</title>
<text class="gm-label" x="572.0" y="621.0">Direction of order: the J-sec…</text>
</a>
<a href="../theorems/T20/" class="gm-node" data-id="T20" data-label="central-charge laws: k=n-4; density holds iff q&gt;=2, zero on (3,1); type=(p-2,q-2)">
<rect x="550.0" y="632.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T20 — Central-charge laws: k=n-4; density holds iff q&gt;=2, zero on (3,1); type=(p-2,q-2)</title>
<text class="gm-label" x="572.0" y="645.0">Central-charge laws: k=n-4; d…</text>
</a>
<a href="../theorems/T21/" class="gm-node" data-id="T21" data-label="junction of t16 and t20: the cocycle w=eps_w (x) eta'_core is forced; heisenberg part required to be metric">
<rect x="550.0" y="656.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T21 — Junction of T16 and T20: the cocycle w=eps_W (x) eta'_core is forced; Heisenberg part required to be metric</title>
<text class="gm-label" x="572.0" y="669.0">Junction of T16 and T20: the …</text>
</a>
<a href="../theorems/T22/" class="gm-node" data-id="T22" data-label="lambda-slot of the charge mu^2/nu equals the lambda ruler; lambda_ext=-eta'; contraction unity of branches">
<rect x="550.0" y="680.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T22 — Lambda-slot of the charge mu^2/nu equals the Lambda ruler; lambda_ext=-eta'; contraction unity of branches</title>
<text class="gm-label" x="572.0" y="693.0">Lambda-slot of the charge mu^…</text>
</a>
<a href="../theorems/T23/" class="gm-node" data-id="T23" data-label="child-metric law (out-of-sample corrected): the layer is live iff g=(p-2,q-2)!=0; shell orthogonal to orbit">
<rect x="550.0" y="704.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T23 — Child-metric law (out-of-sample corrected): the layer is live iff G=(p-2,q-2)!=0; shell orthogonal to orbit</title>
<text class="gm-label" x="572.0" y="717.0">Child-metric law (out-of-samp…</text>
</a>
<a href="../theorems/T24/" class="gm-node" data-id="T24" data-label="lambda^2-valued cocycle proved for all cases: center=lambda^2 w; t21 is the case dim w=2">
<rect x="550.0" y="728.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T24 — Lambda^2-valued cocycle proved for all cases: center=Lambda^2 W; T21 is the case dim W=2</title>
<text class="gm-label" x="572.0" y="741.0">Lambda^2-valued cocycle prove…</text>
</a>
<a href="../theorems/T25/" class="gm-node" data-id="T25" data-label="depth discriminant n^2=0 vs n^2!=0 switches antisymmetric to symmetric; the t24 boundary is sharp">
<rect x="550.0" y="752.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T25 — Depth discriminant N^2=0 vs N^2!=0 switches antisymmetric to symmetric; the T24 boundary is sharp</title>
<text class="gm-label" x="572.0" y="765.0">Depth discriminant N^2=0 vs N…</text>
</a>
<a href="../theorems/T26/" class="gm-node" data-id="T26" data-label="w40 re-reading of the cell's time axis: umbrella of 7 laws (see t26.1-t26.7)">
<rect x="550.0" y="776.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T26 — W40 re-reading of the cell's time axis: umbrella of 7 laws (see T26.1-T26.7)</title>
<text class="gm-label" x="572.0" y="789.0">W40 re-reading of the cell's …</text>
</a>
<a href="../theorems/T26.1/" class="gm-node" data-id="T26.1" data-label="threshold-weight law: t*=d (symbolic for all d); the dirac case is live iff 0&lt;t&lt;d">
<rect x="550.0" y="800.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T26.1 — Threshold-weight law: t*=d (symbolic for all d); the Dirac case is live iff 0&lt;t&lt;d</title>
<text class="gm-label" x="572.0" y="813.0">Threshold-weight law: t*=d (s…</text>
</a>
<a href="../theorems/T26.2/" class="gm-node" data-id="T26.2" data-label="two-weight discriminant law: t=s+(d-1) (symbolic for all d); z2 symmetry exchanges s and t">
<rect x="550.0" y="824.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T26.2 — Two-weight discriminant law: t=s+(d-1) (symbolic for all d); Z2 symmetry exchanges s and t</title>
<text class="gm-label" x="572.0" y="837.0">Two-weight discriminant law: …</text>
</a>
<a href="../theorems/T26.3/" class="gm-node" data-id="T26.3" data-label="commutant ladder of labels: block ranks (1,m-1,d-m), dimensions 1/2/2/1 (d=2), 1/2/3/2 (d=3), 1/2/3/3 (d&gt;=4)">
<rect x="550.0" y="848.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T26.3 — Commutant ladder of labels: block ranks (1,m-1,d-m), dimensions 1/2/2/1 (d=2), 1/2/3/2 (d=3), 1/2/3/3 (d&gt;=4)</title>
<text class="gm-label" x="572.0" y="861.0">Commutant ladder of labels: b…</text>
</a>
<a href="../theorems/T26.4/" class="gm-node" data-id="T26.4" data-label="nodal threads plus forced transversality (schur argument, isotropic cone)">
<rect x="550.0" y="872.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T26.4 — Nodal threads plus forced transversality (Schur argument, isotropic cone)</title>
<text class="gm-label" x="572.0" y="885.0">Nodal threads plus forced tra…</text>
</a>
<a href="../theorems/T26.5/" class="gm-node" data-id="T26.5" data-label="period p_col=d+1 plus dual column equals the circle 2π/(d+1) (pontryagin); 1a+1b+(d-1) holes">
<rect x="550.0" y="896.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T26.5 — Period P_col=d+1 plus dual column equals the circle 2π/(d+1) (Pontryagin); 1A+1B+(d-1) holes</title>
<text class="gm-label" x="572.0" y="909.0">Period P_col=d+1 plus dual co…</text>
</a>
<a href="../theorems/T26.6/" class="gm-node" data-id="T26.6" data-label="reading-b: time equals a labeled coordinate (invariant at q=1; branches at q&gt;=2)">
<rect x="550.0" y="920.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T26.6 — Reading-B: time equals a labeled coordinate (invariant at q=1; branches at q&gt;=2)</title>
<text class="gm-label" x="572.0" y="933.0">Reading-B: time equals a labe…</text>
</a>
<a href="../theorems/T26.7/" class="gm-node" data-id="T26.7" data-label="cartan torus: nodes equal the center z/(d+1), the barycenters of the alcoves (symbolic for all d)">
<rect x="550.0" y="944.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T26.7 — Cartan torus: nodes equal the center Z/(d+1), the barycenters of the alcoves (symbolic for all d)</title>
<text class="gm-label" x="572.0" y="957.0">Cartan torus: nodes equal the…</text>
</a>
<a href="../theorems/T27/" class="gm-node" data-id="T27" data-label="clock-ness scale: none at q=0; canonical time at q=1; split at min(p,q)&gt;=2">
<rect x="550.0" y="968.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T27 — Clock-ness scale: none at q=0; canonical time at q=1; split at min(p,q)&gt;=2</title>
<text class="gm-label" x="572.0" y="981.0">Clock-ness scale: none at q=0…</text>
</a>
<a href="../theorems/T28/" class="gm-node" data-id="T28" data-label="count of ie classes: c(q)=2^q, ctilde=0/1/2/4; q=1 gives a unique class up to the arrow">
<rect x="550.0" y="992.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T28 — Count of ie classes: C(q)=2^q, Ctilde=0/1/2/4; q=1 gives a unique class up to the arrow</title>
<text class="gm-label" x="572.0" y="1005.0">Count of ie classes: C(q)=2^q…</text>
</a>
<a href="../theorems/T29/" class="gm-node" data-id="T29" data-label="collapse and three regimes: generic detuning yields 1 clock; criterion is ctilde=1">
<rect x="550.0" y="1016.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T29 — Collapse and three regimes: generic detuning yields 1 clock; criterion is Ctilde=1</title>
<text class="gm-label" x="572.0" y="1029.0">Collapse and three regimes: g…</text>
</a>
<a href="../theorems/T3/" class="gm-node" data-id="T3" data-label="antisymmetry of the preorder (idempotency of closure, no cycles)">
<rect x="550.0" y="1040.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T3 — Antisymmetry of the preorder (idempotency of closure, no cycles)</title>
<text class="gm-label" x="572.0" y="1053.0">Antisymmetry of the preorder …</text>
</a>
<a href="../theorems/T30/" class="gm-node" data-id="T30" data-label="participation law: an axis participates iff there is diophantine resonance; blind to the axis sign">
<rect x="550.0" y="1064.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T30 — Participation law: an axis participates iff there is Diophantine resonance; blind to the axis sign</title>
<text class="gm-label" x="572.0" y="1077.0">Participation law: an axis pa…</text>
</a>
<a href="../theorems/T31/" class="gm-node" data-id="T31" data-label="rotation-holder: the time-column is fixed, orthogonal to the space that moves in orbit">
<rect x="550.0" y="1088.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T31 — Rotation-holder: the time-column is fixed, orthogonal to the space that moves in orbit</title>
<text class="gm-label" x="572.0" y="1101.0">Rotation-holder: the time-col…</text>
</a>
<a href="../theorems/T32/" class="gm-node" data-id="T32" data-label="first minus from participation + character bridge: minus equals the sign of the pontryagin dual; hole #2 closed kinematically">
<rect x="550.0" y="1112.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T32 — First minus from participation + character bridge: minus equals the sign of the Pontryagin dual; hole #2 closed kinematically</title>
<text class="gm-label" x="572.0" y="1125.0">First minus from participatio…</text>
</a>
<a href="../theorems/T33/" class="gm-node" data-id="T33" data-label="two-component structure forced (0 free parameters) {i,h}; v^2=(1/2)*trm=(d+1)^2/(2d) (d=2-&gt;9/4, d=3-&gt;8/3); sub-bridge #3(a) sealed">
<rect x="550.0" y="1136.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T33 — Two-component structure forced (0 free parameters) {I,H}; v^2=(1/2)*trM=(d+1)^2/(2d) (d=2-&gt;9/4, d=3-&gt;8/3); sub-bridge #3(a) sealed</title>
<text class="gm-label" x="572.0" y="1149.0">Two-component structure force…</text>
</a>
<a href="../theorems/T34/" class="gm-node" data-id="T34" data-label="cell-locality theorem: rank-1 center implies a single clock; hole #1 dissolved">
<rect x="550.0" y="1160.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T34 — Cell-locality theorem: rank-1 center implies a single clock; hole #1 dissolved</title>
<text class="gm-label" x="572.0" y="1173.0">Cell-locality theorem: rank-1…</text>
</a>
<a href="../theorems/T35/" class="gm-node" data-id="T35" data-label="tower of commutatives, anchored sl-gl ladder: n_iso is the ladder of the t26.3 commutant">
<rect x="550.0" y="1184.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T35 — Tower of commutatives, anchored sl-gl ladder: N_iso is the ladder of the T26.3 commutant</title>
<text class="gm-label" x="572.0" y="1197.0">Tower of commutatives, anchor…</text>
</a>
<a href="../theorems/T36/" class="gm-node" data-id="T36" data-label="the arrow is a single bit of realization (5 legs + wedge witness ver:2): no a/b discriminator exists">
<rect x="550.0" y="1208.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T36 — The arrow is a single bit of realization (5 legs + wedge witness ver:2): no A/B discriminator exists</title>
<text class="gm-label" x="572.0" y="1221.0">The arrow is a single bit of …</text>
</a>
<a href="../theorems/T37/" class="gm-node" data-id="T37" data-label="global clock (t34 to the lattice): pair-locality=1 implies a single dial on the connected lattice">
<rect x="550.0" y="1232.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T37 — Global clock (T34 to the lattice): pair-locality=1 implies a single dial on the connected lattice</title>
<text class="gm-label" x="572.0" y="1245.0">Global clock (T34 to the latt…</text>
</a>
<a href="../theorems/T38/" class="gm-node" data-id="T38" data-label="registry wrapper (section 12) over the a-* arc, not a separate derivation: mirror-assembly law (two-stage freezing + kappa=lambda stabilization)">
<rect x="550.0" y="1256.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T38 — Registry wrapper (Section 12) over the A-* arc, not a separate derivation: mirror-assembly law (two-stage freezing + kappa=Lambda stabilization)</title>
<text class="gm-label" x="572.0" y="1269.0">Registry wrapper (Section 12)…</text>
</a>
<a href="../theorems/T39/" class="gm-node" data-id="T39" data-label="registry wrapper (section 12) over the a-* arc, not a separate derivation: the sign is non-derivable (spontaneity equals non-derivability) + inheritance">
<rect x="550.0" y="1280.0" width="18" height="18" rx="3" fill="#c98a2c" stroke="#8a5c16" stroke-width="1.2" />
<title>T39 — Registry wrapper (Section 12) over the A-* arc, not a separate derivation: the sign is non-derivable (spontaneity equals non-derivability) + inheritance</title>
<text class="gm-label" x="572.0" y="1293.0">Registry wrapper (Section 12)…</text>
</a>
<a href="../theorems/T4/" class="gm-node" data-id="T4" data-label="enumeration of niches: the coordinate terminal is finite (2 -&gt; 14 classes)">
<rect x="550.0" y="1304.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T4 — Enumeration of niches: the coordinate terminal is finite (2 -&gt; 14 classes)</title>
<text class="gm-label" x="572.0" y="1317.0">Enumeration of niches: the co…</text>
</a>
<a href="../theorems/T5/" class="gm-node" data-id="T5" data-label="minimality implies (3,1): p&gt;=3 and q&gt;=1 implies n&gt;=4; the minimum is uniquely (3,1); (2,2) fails">
<rect x="550.0" y="1328.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T5 — Minimality implies (3,1): p&gt;=3 and q&gt;=1 implies n&gt;=4; the minimum is uniquely (3,1); (2,2) fails</title>
<text class="gm-label" x="572.0" y="1341.0">Minimality implies (3,1): p&gt;=…</text>
</a>
<a href="../theorems/T6/" class="gm-node" data-id="T6" data-label="coincidence of records: terminal equals stabilizers, 14/14 down to the trivial (3,3) case">
<rect x="550.0" y="1352.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T6 — Coincidence of records: terminal equals stabilizers, 14/14 down to the trivial (3,3) case</title>
<text class="gm-label" x="572.0" y="1365.0">Coincidence of records: termi…</text>
</a>
<a href="../theorems/T7/" class="gm-node" data-id="T7" data-label="square structures: a^2=+1 iff p=q; a^2=-1 iff p,q both even; implies (3,1) is rigid">
<rect x="550.0" y="1376.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T7 — Square structures: A^2=+1 iff p=q; A^2=-1 iff p,q both even; implies (3,1) is rigid</title>
<text class="gm-label" x="572.0" y="1389.0">Square structures: A^2=+1 iff…</text>
</a>
<a href="../theorems/T8/" class="gm-node" data-id="T8" data-label="three generic cases (full algebra, zero symmetry, cartan centralizer)">
<rect x="550.0" y="1400.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T8 — Three generic cases (full algebra, zero symmetry, Cartan centralizer)</title>
<text class="gm-label" x="572.0" y="1413.0">Three generic cases (full alg…</text>
</a>
<a href="../theorems/T9/" class="gm-node" data-id="T9" data-label="mirror log-pair ln det a; sign law (-1)^q*det a&gt;=0, equivalent to pf=0">
<rect x="550.0" y="1424.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>T9 — Mirror log-pair ln det A; sign law (-1)^q*det A&gt;=0, equivalent to Pf=0</title>
<text class="gm-label" x="572.0" y="1437.0">Mirror log-pair ln det A; sig…</text>
</a>
<a href="../theorems/W41/" class="gm-node" data-id="W41" data-label="w41 coxeter-time axiom: time as the action of &lt;c&gt; of order h_cox, won for all d; spectral monism refuted (3 closed forms)">
<rect x="550.0" y="1448.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>W41 — W41 Coxeter-time axiom: time as the action of &lt;c&gt; of order h_cox, won for all d; spectral monism refuted (3 closed forms)</title>
<text class="gm-label" x="572.0" y="1461.0">W41 Coxeter-time axiom: time …</text>
</a>
<a href="../theorems/id-2.1/" class="gm-node" data-id="id-2.1" data-label="center z/(d+1), the column dual, and z/p_col are identified">
<rect x="810.0" y="56.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>id-2.1 — Center Z/(d+1), the column dual, and Z/P_col are identified</title>
<text class="gm-label" x="832.0" y="69.0">Center Z/(d+1), the column du…</text>
</a>
<a href="../theorems/id-2.2/" class="gm-node" data-id="id-2.2" data-label="w_0, the map k &lt;-&gt; -k, and the epsilon-flip are identified">
<rect x="810.0" y="80.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>id-2.2 — w_0, the map k &lt;-&gt; -k, and the epsilon-flip are identified</title>
<text class="gm-label" x="832.0" y="93.0">w_0, the map k &lt;-&gt; -k, and th…</text>
</a>
<a href="../theorems/id-2.3/" class="gm-node" data-id="id-2.3" data-label="b = sigma_x composed with (k &lt;-&gt; -k) is the unique reflector; bare k&lt;-&gt;-k or sigma_x alone are not symmetries">
<rect x="810.0" y="104.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>id-2.3 — B = sigma_x composed with (k &lt;-&gt; -k) is the unique reflector; bare k&lt;-&gt;-k or sigma_x alone are not symmetries</title>
<text class="gm-label" x="832.0" y="117.0">B = sigma_x composed with (k …</text>
</a>
<a href="../theorems/id-2.4/" class="gm-node" data-id="id-2.4" data-label="sigma_x (bipartite), sigma_z (chirality), and b (mirror) are three distinct z_2 involutions">
<rect x="810.0" y="128.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>id-2.4 — sigma_x (bipartite), sigma_z (chirality), and B (mirror) are three distinct Z_2 involutions</title>
<text class="gm-label" x="832.0" y="141.0">sigma_x (bipartite), sigma_z …</text>
</a>
<a href="../theorems/id-2.5/" class="gm-node" data-id="id-2.5" data-label="(marked bond tensor orientation) equals a single d_h orbit, d_h=&lt;c,w_0&gt;">
<rect x="810.0" y="152.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>id-2.5 — (marked bond tensor orientation) equals a single D_h orbit, D_h=&lt;c,w_0&gt;</title>
<text class="gm-label" x="832.0" y="165.0">(marked bond tensor orientati…</text>
</a>
<a href="../theorems/id-2.6/" class="gm-node" data-id="id-2.6" data-label="no a/b discriminator exists (beta-inversion swaps them); the residual is 1 bit (fork ii)">
<rect x="810.0" y="176.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>id-2.6 — No A/B discriminator exists (beta-inversion swaps them); the residual is 1 bit (FORK II)</title>
<text class="gm-label" x="832.0" y="189.0">No A/B discriminator exists (…</text>
</a>
<a href="../theorems/id-2.7/" class="gm-node" data-id="id-2.7" data-label="b1 '+/- as the elementary unit of the two-point structural link', a seed, not a measurement">
<rect x="810.0" y="200.0" width="18" height="18" rx="3" fill="#8a8a8a" stroke="#5a5a5a" stroke-width="1.2" />
<title>id-2.7 — B1 '+/- as the elementary unit of the two-point structural link', a seed, not a measurement</title>
<text class="gm-label" x="832.0" y="213.0">B1 '+/- as the elementary uni…</text>
</a>
<a href="../theorems/seam.1/" class="gm-node" data-id="seam.1" data-label="seam joint of the -2 -&gt; -1 break (t29 collapse) leads to loss of column-tick equality">
<rect x="1070.0" y="56.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>seam.1 — Seam joint of the -2 -&gt; -1 break (T29 collapse) leads to loss of column-tick equality</title>
<text class="gm-label" x="1092.0" y="69.0">Seam joint of the -2 -&gt; -1 br…</text>
</a>
<a href="../theorems/seam.2/" class="gm-node" data-id="seam.2" data-label="the sign of the first minus (t32) leads to the sign of the circle dual (pontryagin dual)">
<rect x="1070.0" y="80.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>seam.2 — The sign of the first minus (T32) leads to the sign of the circle dual (Pontryagin dual)</title>
<text class="gm-label" x="1092.0" y="93.0">The sign of the first minus (…</text>
</a>
<a href="../theorems/seam.3/" class="gm-node" data-id="seam.3" data-label="center z/(d+1) (t32 nodes + t34 locality) leads to nodes = center = barycenters of t26.7">
<rect x="1070.0" y="104.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>seam.3 — Center Z/(d+1) (T32 nodes + T34 locality) leads to nodes = center = barycenters of T26.7</title>
<text class="gm-label" x="1092.0" y="117.0">Center Z/(d+1) (T32 nodes + T…</text>
</a>
<a href="../theorems/seam.4/" class="gm-node" data-id="seam.4" data-label="two-component structure (t33) leads to classes 1a+1b of the t26.5 column">
<rect x="1070.0" y="128.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>seam.4 — Two-component structure (T33) leads to classes 1A+1B of the T26.5 column</title>
<text class="gm-label" x="1092.0" y="141.0">Two-component structure (T33)…</text>
</a>
<a href="../theorems/seam.5/" class="gm-node" data-id="seam.5" data-label="the sl-gl tower (t35) leads to the t26.3 commutant ladder">
<rect x="1070.0" y="152.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>seam.5 — The sl-gl tower (T35) leads to the T26.3 commutant ladder</title>
<text class="gm-label" x="1092.0" y="165.0">The sl-gl tower (T35) leads t…</text>
</a>
<a href="../theorems/seam.6/" class="gm-node" data-id="seam.6" data-label="arrow/class-count (t28+t36 from w41): coxeter action c leads to epsilon -&gt; -epsilon, 1 bit, d_h orbit">
<rect x="1070.0" y="176.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>seam.6 — Arrow/class-count (T28+T36 from W41): Coxeter action c leads to epsilon -&gt; -epsilon, 1 bit, D_h orbit</title>
<text class="gm-label" x="1092.0" y="189.0">Arrow/class-count (T28+T36 fr…</text>
</a>
<a href="../theorems/seam.7/" class="gm-node" data-id="seam.7" data-label="globality (s1018): pair-locality=1 leads to a single dial on the connected lattice">
<rect x="1070.0" y="200.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>seam.7 — Globality (S1018): pair-locality=1 leads to a single dial on the connected lattice</title>
<text class="gm-label" x="1092.0" y="213.0">Globality (S1018): pair-local…</text>
</a>
<a href="../theorems/Psi-def/" class="gm-node" data-id="Psi-def" data-label="ϸ-definition of time and its dual quantity via tact↔circle pair">
<rect x="1330.0" y="56.0" width="18" height="18" rx="3" fill="#7a4fb0" stroke="#4f3373" stroke-width="1.2" />
<title>Psi-def — Ϸ-definition of time and its dual quantity via tact↔circle pair</title>
<text class="gm-label" x="1352.0" y="69.0">Ϸ-definition of time and its …</text>
</a>
<a href="../theorems/AX-sig/" class="gm-node" data-id="AX-sig" data-label="signature (3,1) derived, forced by closure of t5 and t7">
<rect x="1590.0" y="56.0" width="18" height="18" rx="3" fill="#2f8f5b" stroke="#1b5c3a" stroke-width="1.2" />
<title>AX-sig — Signature (3,1) derived, forced by closure of T5 and T7</title>
<text class="gm-label" x="1612.0" y="69.0">Signature (3,1) derived, forc…</text>
</a>
<a href="../theorems/boundary-1/" class="gm-node" data-id="boundary-1" data-label="boundary of level -1: 'the dial exists, there is no action'; level -1 is a dial-factory, action lives at level 0/field">
<rect x="1850.0" y="56.0" width="18" height="18" rx="3" fill="#7a4fb0" stroke="#4f3373" stroke-width="1.2" />
<title>boundary-1 — Boundary of level -1: 'the dial exists, there is no action'; level -1 is a dial-factory, action lives at level 0/field</title>
<text class="gm-label" x="1872.0" y="69.0">Boundary of level -1: 'the di…</text>
</a>
</g>
</svg>
</div>
<style>
.gm-wrap{max-width:100%;overflow-x:auto;border:1px solid var(--md-default-fg-color--lightest,#ccc);
  border-radius:6px;padding:8px;}
.gm-controls{display:flex;gap:8px;align-items:center;margin-bottom:8px;flex-wrap:wrap;}
.gm-search{padding:4px 8px;border:1px solid #888;border-radius:4px;font-size:0.9em;min-width:220px;}
.gm-legend{font-size:0.8em;margin-bottom:6px;}
.gm-legend-item{display:inline-flex;align-items:center;gap:4px;margin-right:6px;}
.gm-svg{background:transparent;}
.gm-colhead{font-size:11px;font-weight:600;fill:currentColor;opacity:0.6;
  text-transform:uppercase;letter-spacing:0.04em;}
.gm-edge{stroke:currentColor;stroke-opacity:0.18;stroke-width:1;fill:none;}
.gm-edge-neg{stroke:#b8433f;stroke-opacity:0.35;stroke-width:1;fill:none;stroke-dasharray:2,2;}
.gm-label{font-size:10px;fill:currentColor;}
.gm-node rect{cursor:pointer;}
.gm-node:hover rect{stroke-width:2.4;}
.gm-node.gm-dim{opacity:0.18;}
.gm-node.gm-hit rect{stroke-width:2.6;}
a.gm-node, a.gm-node:visited{text-decoration:none;color:inherit;}
</style>
<script>
(function(){
  var box = document.currentScript.previousElementSibling;
  while (box && !box.classList.contains('gm-wrap')) { box = box.previousElementSibling; }
  if (!box) return;
  var input = box.querySelector('.gm-search');
  var nodes = box.querySelectorAll('a.gm-node');
  if (!input) return;
  input.addEventListener('input', function(){
    var q = input.value.trim().toLowerCase();
    var first = null;
    nodes.forEach(function(a){
      var hay = (a.getAttribute('data-id') + ' ' + a.getAttribute('data-label')).toLowerCase();
      var hit = q.length > 0 && hay.indexOf(q) !== -1;
      a.classList.toggle('gm-hit', hit);
      a.classList.toggle('gm-dim', q.length > 0 && !hit);
      if (hit && !first) { first = a; }
    });
    if (first) { first.scrollIntoView({block: 'center', inline: 'center'}); }
  });
})();
</script>

### Notation

- `Box` — the difference of two lattice forms, Lambda(psi,nu) = T_A(psi) - T_col(nu)
- `T_A` — the lattice form on the cell, T(k) = 2 - 2cos(2 pi k / h), an integer only for h in {1,2,3,4,6} (measured window d in {2,3})
- `T_col` — the same lattice form on the column
- `m_0` — the quantity of the arc whose sign is the subject; its value is a representative of the regime, not canon, and is not derived here

<!--graph-map-close-->

---

## A boundary of the map, named explicitly

The graph covers the **prime floor**. The most recent layer of the chain (the terminal step of the
derived series and everything after it) is **not yet represented in the graph** — those links are
cited in Tome I **directly by the seals of the rulings**. So a clean state of the graph checker
applies to the covered part, not to the whole exposition. This is said here so that a "one green
tick" reading is not taken more widely than it holds.
