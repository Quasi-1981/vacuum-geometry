import VGLean.So51r1
import VGLean.So42r1
import VGLean.So33r1
import VGLean.So32r0
import VGLean.So42r0
import VGLean.So33r0

set_option linter.style.header false

-- Axiom-hygiene gate for A2 (all 6 FULL-centralizer isomorphisms). Every A2 theorem must depend
-- ONLY on the mathlib triple [propext, Classical.choice, Quot.sound] — NO sorryAx, NO ofReduceBool
-- (no native_decide leaking Lean.ofReduceBool into the trust base). A1 discipline (J-0428).
/-- info: 'VGLean.A2.So51r1.brC_2_3' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So51r1.brC_2_3
/-- info: 'VGLean.A2.So51r1.B_indep' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So51r1.B_indep
/-- info: 'VGLean.A2.So51r1.psi0_anchor' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So51r1.psi0_anchor
/-- info: 'VGLean.A2.So42r1.brC_2_4' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So42r1.brC_2_4
/-- info: 'VGLean.A2.So33r1.brB_1_4' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So33r1.brB_1_4
/-- info: 'VGLean.A2.So32r0.C_indep' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So32r0.C_indep
/-- info: 'VGLean.A2.So42r0.brC_3_5' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So42r0.brC_3_5
/-- info: 'VGLean.A2.So42r0.B8_comm' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So42r0.B8_comm
/-- info: 'VGLean.A2.So33r0.brC_4_6' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So33r0.brC_4_6
/-- info: 'VGLean.A2.So33r0.B_indep' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So33r0.B_indep
/-- info: 'VGLean.A2.So33r0.N_nilpotent' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.A2.So33r0.N_nilpotent
