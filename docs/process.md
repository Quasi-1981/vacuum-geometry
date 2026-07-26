# Process & roles

*This page describes how the work was organized. It speaks of the programme's process, not of its
mathematics; nothing here is a claim of the chain.*

## The plain statement first

The computations and texts behind this site were produced by AI instances working under human
authorship and direction. That is visible in the artifacts anyway, so it is stated here plainly
instead of being left for the reader to infer. What this page explains is the part that is *not*
visible at a glance: the organization — because in our experience it is the organization, not the
model, that carries the reliability.

## The pyramid of roles

A pyramidal system of three core instances, plus an engineering department:

| Role | Function |
|:--|:--|
| **Omega** | project architect: plans, adjudicates, sole writer of the shared project memory |
| **Alpha** | hypothesis computation and blind probes |
| **Beta** | independent verification and journaling |
| **Sigma** | engineering executor: production, tooling, publication |

These are AI instances cast into **roles**; the roles are permanent, the instances are not. All of
it runs under a single human author, whose word gates every step that leaves the workshop.

## Why the roles are separated

Earlier iterations showed that combining different roles within a single instance produces a purely
psychological imbalance inside that instance: the work destabilizes, and the method of blind probes
stops being possible — an instance that computed a hypothesis cannot blindly grade it. Role
separation is therefore not an aesthetic choice; it is what makes a probe's PASS worth anything.

The same separation applies to verification: the instance that produced a result never writes its
own verdict. A separate instance re-derives or re-runs the claim and records the outcome
independently.

## The rig

Architecturally the solution is fairly trivial: a shared memory bank for the project, individual
memory banks for each type of instance, and dependency graphs within the mathematical models. The
one component that is perhaps non-trivial is the communication bus between the instances. A minimal
set of working rules, hierarchy, and methodology is mandatory: without it the system — like any
system — drifts toward simpler solutions.

## What the reader sees in the artifacts

- **`S###`** — a work id reserved on the communication bus before a probe is run; probe files carry
  it in their names. These ids resolve to the published probe sources.
- **Attribution lines** (e.g. *probe: Sigma · verification: Beta · adjudication: Omega*) — the
  collective work is left visible, not erased.
- **Verification records** are internal to the project's review process; public pages state the
  fact of independent verification rather than the record number, since the internal journal is not
  published.

## The obvious question

Everyone is trying AI in research work, and mostly it does not hold up. In our experience it fails
exactly where a single instance holds all the roles at once: generator, grader, and archivist in
one context, with nothing blind and nothing independent. The pyramid above is the answer we
arrived at, and the site is built so that the reader does not have to take any of this on trust:
every claim carries its status, its probe, and the probe can be re-run.

## Feedback & errata

Found a broken link in the chain, an error in a probe, or a claim that does not survive your own
re-derivation? Open an issue in the
[GitHub repository](https://github.com/Quasi-1981/vacuum-geometry/issues). Refutation attempts are
welcome — a failed one is worth more to this project than praise.
