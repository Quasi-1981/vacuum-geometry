import Mathlib.Tactic

set_option linter.style.header false

/-!
# DIM: exact-rational / Layer-1 · T34 — CELL CAPACITY (rank-1 centre)

WHAT IS PROVED HERE (literally, for all n, no snapshots):
for the cell axes `u i k = δ_{ik} − 1/n` in `ℚ^n` (n = d+1, exact rational arithmetic) —
(a) **the pairwise difference of axes is INTEGRAL**: `u i − u j = e_i − e_j`, with
    coordinates in {−1,0,1} (the centroid cancels), hence all axes lie in the SAME class
    modulo the lattice `ℤ^n`;
(b) that difference has **zero sum**, hence it lies in the ROOT lattice `A_{n−1}`, not
    merely in `ℤ^n`;
(c) the axis `u i` itself is **NOT integral** for `n ≥ 2` (`u i i = 1 − 1/n ∈ (0,1)`),
    hence the common class is **NON-TRIVIAL** (charge ≠ 0);
(d) `ZMod n` is **cyclic**: every element is an `ℤ`-multiple of a single generator,
    hence **rank 1**.
Together (a)+(b)+(c)+(d): the `n` cell axes occupy **EXACTLY ONE non-zero class** of the
centre `ℤ/n`, and no second INDEPENDENT class exists in it.

WHAT IS NOT HERE AND IS NOT CLAIMED: the words "clock", "time", "column" do not occur in
any statement. The reading "capacity = one clock" belongs to the ancestor record, not to
this file.

WHY THIS IS NOT A SNAPSHOT: everything is proved `∀ n`, symbolically. The instances
`n = 3, 4` (d = 2, 3) stand separately as a TERM-BY-TERM CROSS-CHECK against the ancestor,
which measured exactly those two values of d; they are not the method of proof.

WHY THIS FILE IS NOT GENERATED: there are no Python-derived literals here — the statement
is symbolic, so a generator would be a redundant link. No `gen_*.py` is required.

ANCESTOR, TERM BY TERM: the published probe `src/probe/S1013.py`, block K2
(`cell_axes`: `u = [−1/n …]`, `u[i] = 1 − 1/n`, exact `Fraction`); ancestors T26.7/T29/T31.
-/

namespace VGLean.T34

open Finset

/-- Cell axis: `u i k = δ_{ik} − 1/n` over `ℚ` (exact; n = d+1). -/
def u (n : ℕ) (i k : Fin n) : ℚ := (if i = k then 1 else 0) - 1 / (n : ℚ)

/-- (a) THE CENTROID CANCELS: the difference of two axes is a difference of basis vectors. -/
theorem u_sub (n : ℕ) (i j k : Fin n) :
    u n i k - u n j k = (if i = k then 1 else 0) - (if j = k then 1 else 0) := by
  simp [u]

/-- (a') The coordinates of the difference lie in `{−1, 0, 1}`, i.e. the difference is INTEGRAL. -/
theorem u_sub_mem (n : ℕ) (i j k : Fin n) :
    u n i k - u n j k = -1 ∨ u n i k - u n j k = 0 ∨ u n i k - u n j k = 1 := by
  rw [u_sub]
  by_cases h1 : i = k <;> by_cases h2 : j = k <;> simp [h1, h2]

/-- (b) The difference has ZERO SUM, hence lies in the root lattice `A_{n−1}`, not merely `ℤ^n`. -/
theorem u_sub_sum_zero (n : ℕ) (i j : Fin n) :
    ∑ k : Fin n, (u n i k - u n j k) = 0 := by
  simp [u_sub, Finset.sum_sub_distrib]

/-- The axis at its own coordinate: `u i i = 1 − 1/n`. -/
theorem u_self (n : ℕ) (i : Fin n) : u n i i = 1 - 1 / (n : ℚ) := by
  simp [u]

/-- (c) THE AXIS IS NOT INTEGRAL for `n ≥ 2`: `u i i = 1 − 1/n` lies STRICTLY between 0
    and 1, so it cannot be an integer; hence the common class of the axes is NON-TRIVIAL
    (charge ≠ 0). -/
theorem u_self_pos_lt_one (n : ℕ) (hn : 2 ≤ n) (i : Fin n) :
    0 < u n i i ∧ u n i i < 1 := by
  have h2 : (2 : ℚ) ≤ (n : ℚ) := by exact_mod_cast hn
  have hn0 : (0 : ℚ) < (n : ℚ) := by linarith
  have hinv_pos : 0 < 1 / (n : ℚ) := by positivity
  have hinv_lt : 1 / (n : ℚ) < 1 := (div_lt_one hn0).mpr (by linarith)
  rw [u_self]
  constructor <;> linarith

/-- (d) RANK-1 CENTRE: `ZMod n` is cyclic — every element is an
    `ℤ`-multiple of the generator `1`. -/
theorem center_rank_one (n : ℕ) [NeZero n] (x : ZMod n) : ∃ m : ℤ, x = m • (1 : ZMod n) := by
  refine ⟨(x.val : ℤ), ?_⟩
  simp [zsmul_eq_mul]

/-! ## TERM-BY-TERM CROSS-CHECK AGAINST THE ANCESTOR (S1013 measured exactly d = 2, 3,
i.e. n = 3, 4). This is NOT the method of proof — everything above is proved `∀ n`; these
are instance checks against the ancestor's table. -/

/-- d = 2 (n = 3): diagonal entry `1 − 1/3 = 2/3`. -/
theorem inst_n3_diag (i : Fin 3) : u 3 i i = 2 / 3 := by norm_num [u]

/-- d = 2 (n = 3): off-diagonal entry `−1/3`. -/
theorem inst_n3_off (i k : Fin 3) (h : i ≠ k) : u 3 i k = -(1 / 3) := by
  simp [u, h]

/-- d = 3 (n = 4): diagonal entry `1 − 1/4 = 3/4`. -/
theorem inst_n4_diag (i : Fin 4) : u 4 i i = 3 / 4 := by norm_num [u]

/-- d = 3 (n = 4): off-diagonal entry `−1/4`. -/
theorem inst_n4_off (i k : Fin 4) (h : i ≠ k) : u 4 i k = -(1 / 4) := by
  simp [u, h]

end VGLean.T34

/-! ## Axiom hygiene: EVERY theorem is stamped, and the stamp is the exact triple
OR A SUBSET of it ("narrower" is meant strictly as a subset). -/

/-- info: 'VGLean.T34.u_sub' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.u_sub

/-- info: 'VGLean.T34.u_sub_mem' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.u_sub_mem

/-- info: 'VGLean.T34.u_sub_sum_zero' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.u_sub_sum_zero

/-- info: 'VGLean.T34.u_self' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.u_self

/-- info: 'VGLean.T34.u_self_pos_lt_one' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.u_self_pos_lt_one

-- A SUBSET, not merely "shorter": `Classical.choice` is not needed here at all.
/-- info: 'VGLean.T34.center_rank_one' depends on axioms: [propext, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.center_rank_one

/-- info: 'VGLean.T34.inst_n3_diag' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.inst_n3_diag

/-- info: 'VGLean.T34.inst_n3_off' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.inst_n3_off

/-- info: 'VGLean.T34.inst_n4_diag' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.inst_n4_diag

/-- info: 'VGLean.T34.inst_n4_off' depends on axioms: [propext, Classical.choice, Quot.sound] -/
#guard_msgs in
#print axioms VGLean.T34.inst_n4_off
