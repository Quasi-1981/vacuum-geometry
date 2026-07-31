import VGLean.B1Resolvent

set_option linter.style.header false

/-!
# Axiom-hygiene gate for the B1 preprint-4 §2 resolvent identity

Every B1 §2 theorem must depend ONLY on the mathlib triple `[propext, Classical.choice, Quot.sound]`
— NO `sorryAx`, NO `Lean.ofReduceBool` (no `native_decide`).  Same discipline as A1–A3 / C2.
Run `lake build VGLean.B1AxiomCheck` and read the emitted axiom sets by hand.
-/

/-- info: 'VGLean.B1.Resolvent.O_eq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.B1.Resolvent.O_eq
/-- info: 'VGLean.B1.Resolvent.H_sq' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.B1.Resolvent.H_sq
/-- info: 'VGLean.B1.Resolvent.resolvent_id' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.B1.Resolvent.resolvent_id
/-- info: 'VGLean.B1.Resolvent.Adj_mul_O' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.B1.Resolvent.Adj_mul_O
