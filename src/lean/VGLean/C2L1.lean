import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Matrix.Mul
import Mathlib.Logic.Equiv.Fintype
import Mathlib.Tactic

set_option linter.style.header false
set_option linter.style.longLine false

/-!
# C2-cell L1 — "geometry from bits", milestone L1, STAGE 1 (orbit characterization)

The heart of C2: the space of `S_n`-invariant symmetric bilinear forms on the standard
representation is **1-dimensional** (∀ n ≥ 2, char 0).

STAGE 1 (this file): the matrix core.  An `S_n`-invariant form on `ℚ^(n+2)` (coordinate
permutation) corresponds to a matrix `M` with `M (σ i) (σ j) = M i j` ∀σ.  Such an `M` is
constant on the **two** orbits of the diagonal `S_n`-action on ordered pairs — the diagonal
`{(i,i)}` and the off-diagonal `{(i,j) : i ≠ j}` (2-transitivity) — hence
`M = (a - b)·I + b·O`, i.e. `span{I, O}`, `dim = 2`.  (`O` = all-ones = the `J`-form.)
The 1-dimensionality on `W = ker(augmentation)` (where `O` restricts to `0`) is STAGE 2.

`ι = Fin (n + 2)` encodes `n ≥ 2` (cardinality ≥ 2, so `0 ≠ 1` always available).
Char-0 field `ℚ`.  Layer-2 (form = geometry of an alphabet) NOT taken here.  2-transitivity
witness: `Equiv.Perm.exists_extending_pair` (no `native_decide`, kernel triple only).
-/

namespace VGLean.C2.L1

open Matrix

variable {n : ℕ}

/-- Index type with cardinality `≥ 2` (`n ≥ 2` sites). -/
abbrev ι (n : ℕ) := Fin (n + 2)

/-- All-ones matrix (the `J`/`O` form: `O x y = (∑ xᵢ)(∑ yⱼ)`). -/
def O : Matrix (ι n) (ι n) ℚ := Matrix.of fun _ _ => 1

/-- `S_n`-invariance of a form-matrix: `M (σ i) (σ j) = M i j` for every permutation `σ`. -/
def Inv (M : Matrix (ι n) (ι n) ℚ) : Prop :=
  ∀ σ : Equiv.Perm (ι n), M.submatrix σ σ = M

/-- Entrywise form of invariance. -/
theorem Inv.apply {M : Matrix (ι n) (ι n) ℚ} (h : Inv M) (σ : Equiv.Perm (ι n)) (i j : ι n) :
    M (σ i) (σ j) = M i j := by
  have := congrFun (congrFun (h σ) i) j
  simpa [Matrix.submatrix_apply] using this

/-- The identity form (`I`, the dot product) is invariant. -/
theorem one_inv : Inv (1 : Matrix (ι n) (ι n) ℚ) := by
  intro σ
  ext i j
  simp [Matrix.submatrix_apply, Matrix.one_apply, σ.injective.eq_iff]

/-- The all-ones form (`O`) is invariant. -/
theorem O_inv : Inv (O : Matrix (ι n) (ι n) ℚ) := by
  intro σ
  ext i j
  simp [O, Matrix.submatrix_apply]

/-- **2-transitivity witness.** For `i ≠ j` there is a permutation sending `0 ↦ i`, `1 ↦ j`. -/
theorem exists_perm_pair {i j : ι n} (hij : i ≠ j) :
    ∃ σ : Equiv.Perm (ι n), σ 0 = i ∧ σ 1 = j := by
  have hf : Function.Injective (![ (0 : ι n), 1 ]) := by
    intro a b hab
    fin_cases a <;> fin_cases b <;> simp_all
  have hg : Function.Injective (![ i, j ]) := by
    intro a b hab
    fin_cases a <;> fin_cases b <;> simp_all
  obtain ⟨σ, hσ⟩ := Equiv.Perm.exists_extending_pair (![ (0 : ι n), 1 ]) (![ i, j ]) hf hg
  exact ⟨σ, by simpa using hσ 0, by simpa using hσ 1⟩

/-- **Orbit characterization (pointwise).** An invariant `M` takes only two values: `M 0 0`
    on the diagonal, `M 0 1` off the diagonal. -/
theorem inv_two_valued {M : Matrix (ι n) (ι n) ℚ} (h : Inv M) (i j : ι n) :
    M i j = if i = j then M 0 0 else M 0 1 := by
  by_cases hij : i = j
  · -- diagonal: move 0 ↦ i by a swap
    subst hij
    rw [if_pos rfl]
    have := h.apply (Equiv.swap 0 i) 0 0
    simpa using this
  · -- off-diagonal: 2-transitivity 0,1 ↦ i,j
    rw [if_neg hij]
    obtain ⟨σ, h0, h1⟩ := exists_perm_pair hij
    have := h.apply σ 0 1
    rw [h0, h1] at this
    exact this

/-- **STAGE 1 result: `Inv M ⟹ M = (a − b)·I + b·O`**, `a = M 0 0`, `b = M 0 1`.
    The invariant forms are exactly `span{I, O}` (dimension 2 on the full space). -/
theorem inv_decomp {M : Matrix (ι n) (ι n) ℚ} (h : Inv M) :
    M = (M 0 0 - M 0 1) • (1 : Matrix (ι n) (ι n) ℚ) + (M 0 1) • O := by
  ext i j
  rw [inv_two_valued h i j]
  simp only [Matrix.add_apply, Matrix.smul_apply, O, Matrix.of_apply, Matrix.one_apply,
    smul_eq_mul, mul_one]
  by_cases hij : i = j
  · subst hij; simp
  · rw [if_neg hij]; simp [hij]

/-! ## STAGE 2 — the bilinear form, restriction to `W = ker(augmentation)`, and `dim = 1`. -/

/-- The bilinear form of a form-matrix: `bF M x y = ∑ᵢ ∑ⱼ xᵢ · Mᵢⱼ · yⱼ`. -/
def bF (M : Matrix (ι n) (ι n) ℚ) (x y : ι n → ℚ) : ℚ := ∑ i, ∑ j, x i * M i j * y j

/-- **Form expansion.** An invariant form splits into a diagonal (`I`) part and an all-ones
    (`O`) part:  `bF M x y = (a − b)·⟨x,y⟩ + b·(∑x)(∑y)`, `a = M 0 0`, `b = M 0 1`. -/
theorem bF_inv {M : Matrix (ι n) (ι n) ℚ} (h : Inv M) (x y : ι n → ℚ) :
    bF M x y = (M 0 0 - M 0 1) * (∑ i, x i * y i)
             + M 0 1 * ((∑ i, x i) * (∑ i, y i)) := by
  have e1 : bF M x y = ∑ i, ∑ j,
      ((M 0 0 - M 0 1) * (if i = j then x i * y i else 0) + M 0 1 * (x i * y j)) := by
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    rw [inv_two_valued h i j]
    by_cases hij : i = j <;> simp [hij] <;> ring
  rw [e1]
  simp only [Finset.sum_add_distrib, ← Finset.mul_sum]
  congr 1
  · congr 1
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.sum_ite_eq Finset.univ i (fun _ => x i * y i)]
    simp
  · rw [Finset.sum_mul_sum]
    refine congrArg _ (Finset.sum_congr rfl fun i _ => ?_)
    rw [Finset.mul_sum]

/-- **Restriction to `W = ker(augmentation)`.** On the subspace `∑ xᵢ = 0` the all-ones part
    vanishes, so every invariant form is a scalar multiple of the dot product:
    `bF M x y = (a − b)·⟨x,y⟩`. -/
theorem bF_on_W {M : Matrix (ι n) (ι n) ℚ} (h : Inv M) {x y : ι n → ℚ}
    (hx : ∑ i, x i = 0) (hy : ∑ i, y i = 0) :
    bF M x y = (M 0 0 - M 0 1) * (∑ i, x i * y i) := by
  rw [bF_inv h, hx, hy]; ring

/-- Witness vector in `W`: `e₀ − e₁` has zero coordinate sum and self-dot `2 ≠ 0` (needs `n ≥ 2`). -/
def wWit : ι n → ℚ := fun i => (if i = 0 then 1 else 0) - (if i = 1 then 1 else 0)

theorem wWit_sum : ∑ i, (wWit : ι n → ℚ) i = 0 := by
  simp [wWit, Finset.sum_sub_distrib, Finset.sum_ite_eq' Finset.univ (0 : ι n),
    Finset.sum_ite_eq' Finset.univ (1 : ι n)]

theorem wWit_dot : ∑ i, (wWit : ι n → ℚ) i * wWit i = 2 := by
  have h01 : (0 : ι n) ≠ 1 := Fin.zero_ne_one
  have pt : ∀ i : ι n, (wWit : ι n → ℚ) i * wWit i
      = (if i = 0 then (1:ℚ) else 0) + (if i = 1 then 1 else 0) := by
    intro i
    by_cases h0 : i = 0
    · subst h0; simp [wWit, h01]
    · by_cases h1 : i = 1
      · subst h1; simp [wWit, Ne.symm h01]
      · simp [wWit, h0, h1]
  rw [Finset.sum_congr rfl fun i _ => pt i]
  simp only [Finset.sum_add_distrib, Finset.sum_ite_eq' Finset.univ (0 : ι n),
    Finset.sum_ite_eq' Finset.univ (1 : ι n), Finset.mem_univ, if_true]
  norm_num

/-- **L1 (dimension one, ∀ n ≥ 2, char 0).** For every `S_n`-invariant form `M`, its restriction
    to `W = ker(augmentation)` is a *unique* scalar multiple of the dot product.  The space of
    invariant symmetric forms on the standard representation is therefore exactly 1-dimensional. -/
theorem L1_dim_one {M : Matrix (ι n) (ι n) ℚ} (h : Inv M) :
    ∃! c : ℚ, ∀ x y : ι n → ℚ, (∑ i, x i = 0) → (∑ i, y i = 0) →
      bF M x y = c * (∑ i, x i * y i) := by
  refine ⟨M 0 0 - M 0 1, fun x y hx hy => bF_on_W h hx hy, ?_⟩
  intro c hc
  have := hc wWit wWit wWit_sum wWit_sum
  rw [bF_on_W h wWit_sum wWit_sum, wWit_dot] at this
  have h2 : (2 : ℚ) ≠ 0 := by norm_num
  field_simp at this
  linarith [this]

end VGLean.C2.L1
