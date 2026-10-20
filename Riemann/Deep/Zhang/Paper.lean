import Riemann.Deep.DH.ZeroFinal

/-!
# Zhang, Preprints.org 202108.0146 v60 (24 July 2026): chain laid down, failing step pinned and refuted

Paper: pair `ρᵢ, ρ̄ᵢ` in the Hadamard product, `ξ(s) = ξ(0) ∏ (βᵢ²/|ρᵢ|² + (s-αᵢ)²/|ρᵢ|²)^{mᵢ}` (eq. 25),
so up to a constant `ξ(s) = c ∏ gᵢ(s)^{mᵢ}`, `gᵢ(s) = 1 + (s-αᵢ)²/βᵢ²` (eq. 29). Chain:
* (22)–(25), Lemma 7: convergence, `∑ 1/βᵢ² < ∞`. Correct.
* Lemma 3 + (16) ⇒ (17): every `gᵢ(s)` equals some `gₗ(1-s)`. Correct (the zero set is symmetric).
* **(19) (fails, `Step19`)**: "due to the uniqueness of the multiplicity mᵢ, the only solution to
  Equation (17) is" `l = i`, since otherwise (19c) `αᵢ + αₗ = 1, βᵢ² = βₗ²` "means the i-th and l-th
  quadruplets of zeros coincide ..., violating the uniqueness of mᵢ".
  But `gᵢ` carries only the conjugate PAIR `αᵢ ± iβᵢ`, not a quadruplet. An off-line quadruplet
  `{α ± iβ, 1-α ± iβ}` is two distinct factors `gᵢ, gₗ` with `αₗ = 1-αᵢ`; (19c) is exactly that case,
  not a contradiction. `l = i` for all `i` is equivalent to RH itself: the step is circular.
* (20)/(30): `αᵢ = 1/2`. Lemma 6 (⇒ of eq. 13) is false.

Lean (all sorry-free, standard axioms):
* `chain`: (17) + `Step19` ⇒ `Lemma6` (the paper's deduction is valid; the input is false).
* `hyp0`, `symm0`: the family `α = (3/10, 7/10, 1/2, 1/2, …)`, `β = (1, 1, 2, 3, …)`, `m ≡ 1` meets every
  hypothesis of Lemma 6 (strip, `β ≠ 0`, `|β|` monotone, `∑ 1/β² < ∞`, distinct pairs, convergence) and
  `F(s) = F(1-s)`. `not_lemma6`, `not_step19`: refuted.
* `dh_zero_in_strip`: the paper's DH remark ("DH zeros are not all in the strip") does not save it:
  the counterexample is in the strip, and DH has a certified off-line zero inside the strip too.
-/

namespace Riemann.Deep.Zhang

open Complex

/-- Paired conjugate factor `gᵢ(s) = 1 + (s - αᵢ)²/βᵢ²` (eq. 14, 29), zeros `αᵢ ± iβᵢ`. -/
noncomputable def g (α β : ℕ → ℝ) (i : ℕ) (s : ℂ) : ℂ := 1 + (s - α i) ^ 2 / (β i : ℂ) ^ 2

/-- `f(s) = ∏ gᵢ(s)^{mᵢ}` (eq. 11, 28). -/
noncomputable def F (α β : ℕ → ℝ) (m : ℕ → ℕ) (s : ℂ) : ℂ := ∏' i, g α β i s ^ m i

/-- Hypotheses of Lemma 6. `distinct`: each factor appears once, with its multiplicity `mᵢ`. -/
structure Hyp (α β : ℕ → ℝ) (m : ℕ → ℕ) : Prop where
  strip : ∀ i, 0 < α i ∧ α i < 1
  ne : ∀ i, β i ≠ 0
  mono : Monotone fun i => |β i|
  summ : Summable fun i => 1 / β i ^ 2
  mult : ∀ i, 1 ≤ m i
  distinct : Function.Injective fun i => (α i, β i ^ 2)
  conv : ∀ s, Multipliable fun i => g α β i s ^ m i

/-- The functional equation `f(s) = f(1-s)`. -/
def Symm (α β : ℕ → ℝ) (m : ℕ → ℕ) : Prop := ∀ s, F α β m s = F α β m (1 - s)

/-- Lemma 6, direction ⇒ of eq. (13). -/
def Lemma6 : Prop := ∀ α β m, Hyp α β m → Symm α β m → ∀ i, α i = 1 / 2

/-- Eq. (17), via (16) and Lemma 3 (true for symmetric products). -/
def Div17 (α β : ℕ → ℝ) : Prop := ∀ i, ∃ l, ∀ s, g α β i s = g α β l (1 - s)

/-- Eq. (19): "the only solution to Equation (17) is" `l = i`. -/
def Step19 : Prop :=
  ∀ α β m, Hyp α β m → Symm α β m → ∀ i l, (∀ s, g α β i s = g α β l (1 - s)) → l = i

/-- Comparing like terms in (19a): `gᵢ(s) = gᵢ(1-s)` forces `αᵢ = 1/2`. -/
lemma alpha_half {α β : ℕ → ℝ} {i : ℕ} (hβ : β i ≠ 0) (h : ∀ s, g α β i s = g α β i (1 - s)) :
    α i = 1 / 2 := by
  have h1 := h (α i)
  simp only [g, sub_self] at h1
  have hb : (β i : ℂ) ^ 2 ≠ 0 := pow_ne_zero _ (ofReal_ne_zero.mpr hβ)
  have h2 : (1 - (α i : ℂ) - α i) ^ 2 / (β i : ℂ) ^ 2 = 0 := by linear_combination -h1
  rw [div_eq_zero_iff, or_iff_left hb] at h2
  have h3 : ((1 - α i - α i : ℝ) : ℂ) = 0 := by push_cast; exact pow_eq_zero_iff two_ne_zero |>.mp h2
  have h4 : 1 - α i - α i = 0 := by exact_mod_cast h3
  linarith

/-- The paper's chain is a valid deduction: (17) + step (19) ⇒ Lemma 6. -/
theorem chain (h17 : ∀ α β m, Hyp α β m → Symm α β m → Div17 α β) (h19 : Step19) : Lemma6 := by
  intro α β m hH hS i
  obtain ⟨l, hl⟩ := h17 α β m hH hS i
  have := h19 α β m hH hS i l hl
  rw [this] at hl
  exact alpha_half (hH.ne i) hl

/-! ## Counterexample: one off-line quadruplet `3/10 ± i, 7/10 ± i`, the rest on the line -/

/-- `α = (3/10, 7/10, 1/2, 1/2, …)`. -/
noncomputable def α0 (i : ℕ) : ℝ := if i = 0 then 3 / 10 else if i = 1 then 7 / 10 else 1 / 2

/-- `β = (1, 1, 2, 3, …)`. -/
noncomputable def β0 (i : ℕ) : ℝ := if i ≤ 1 then 1 else (i : ℝ)

/-- `m ≡ 1`. -/
def m0 (_ : ℕ) : ℕ := 1

lemma β0_add_two (n : ℕ) : β0 (n + 2) = ((n + 2 : ℕ) : ℝ) := by
  simp [β0]

lemma β0_pos (i : ℕ) : 0 < β0 i := by
  unfold β0; split_ifs with h
  · norm_num
  · exact Nat.cast_pos.mpr (by omega)

lemma α0_mem (i : ℕ) : 0 < α0 i ∧ α0 i < 1 := by
  unfold α0; split_ifs <;> norm_num

lemma summ0 : Summable fun i => 1 / β0 i ^ 2 := by
  apply (summable_nat_add_iff 2).mp
  have := (summable_nat_add_iff 2).mpr
    (Real.summable_one_div_nat_pow.mpr one_lt_two : Summable fun n : ℕ => 1 / (n : ℝ) ^ 2)
  refine this.congr fun n => ?_
  rw [β0_add_two]

lemma conv0 (s : ℂ) : Multipliable fun i => g α0 β0 i s ^ m0 i := by
  simp only [m0, pow_one, g]
  apply Complex.multipliable_one_add_of_summable
  refine Summable.of_norm_bounded (summ0.mul_left ((‖s‖ + 1) ^ 2)) fun i => ?_
  have hα := α0_mem i
  have hs : ‖s - (α0 i : ℂ)‖ ≤ ‖s‖ + 1 := by
    calc ‖s - (α0 i : ℂ)‖ ≤ ‖s‖ + ‖(α0 i : ℂ)‖ := norm_sub_le _ _
      _ ≤ ‖s‖ + 1 := by
        rw [Complex.norm_real, Real.norm_eq_abs, abs_of_pos hα.1]; linarith [hα.2]
  have hb := β0_pos i
  calc ‖(s - (α0 i : ℂ)) ^ 2 / (β0 i : ℂ) ^ 2‖ = ‖s - (α0 i : ℂ)‖ ^ 2 / β0 i ^ 2 := by
        rw [norm_div, norm_pow, norm_pow, Complex.norm_real, Real.norm_eq_abs, sq_abs]
    _ ≤ (‖s‖ + 1) ^ 2 / β0 i ^ 2 := by gcongr
    _ = (‖s‖ + 1) ^ 2 * (1 / β0 i ^ 2) := by ring

lemma distinct0 : Function.Injective fun i => (α0 i, β0 i ^ 2) := by
  intro a b h
  simp only [Prod.mk.injEq] at h
  obtain ⟨h1, h2⟩ := h
  have off : ∀ i, 1 < i → α0 i = 1 / 2 := fun i hi => by
    simp [α0, show i ≠ 0 by omega, show i ≠ 1 by omega]
  by_cases ha : a ≤ 1 <;> by_cases hb : b ≤ 1
  · interval_cases a <;> interval_cases b <;> first | rfl | norm_num [α0] at h1
  · rw [off b (by omega)] at h1
    interval_cases a <;> norm_num [α0] at h1
  · rw [off a (by omega)] at h1
    interval_cases b <;> norm_num [α0] at h1
  · have ea : β0 a = a := by simp [β0, ha]
    have eb : β0 b = b := by simp [β0, hb]
    rw [ea, eb] at h2
    have h3 := (pow_left_inj₀ (by positivity) (by positivity) two_ne_zero).mp h2
    exact_mod_cast h3

lemma mono0 : Monotone fun i => |β0 i| := by
  intro a b hab
  simp only [abs_of_pos (β0_pos _)]
  unfold β0
  split_ifs with ha hb hb
  · exact le_rfl
  · exact Nat.one_le_cast.mpr (by omega)
  · omega
  · exact Nat.cast_le.mpr hab

theorem hyp0 : Hyp α0 β0 m0 where
  strip := α0_mem
  ne i := (β0_pos i).ne'
  mono := mono0
  summ := summ0
  mult _ := le_rfl
  distinct := distinct0
  conv := conv0

/-- `s ↦ 1 - s` swaps the factors `g₀, g₁` and fixes the on-line ones. -/
lemma g0_reflect (i : ℕ) (s : ℂ) : g α0 β0 i (1 - s) = g α0 β0 (Equiv.swap 0 1 i) s := by
  rcases i with _ | _ | n
  · simp only [Equiv.swap_apply_left, g, α0, β0]; norm_num; ring
  · simp only [Equiv.swap_apply_right, g, α0, β0]; norm_num; ring
  · rw [Equiv.swap_apply_of_ne_of_ne (by omega) (by omega)]
    simp only [g, α0, β0]; norm_num; ring

theorem symm0 : Symm α0 β0 m0 := by
  intro s
  unfold F
  have : (fun i => g α0 β0 i (1 - s) ^ m0 i) =
      fun i => (fun j => g α0 β0 j s ^ m0 j) (Equiv.swap 0 1 i) := by
    funext i; rw [g0_reflect]; rfl
  rw [this, Equiv.tprod_eq (Equiv.swap 0 1) (fun j => g α0 β0 j s ^ m0 j)]

/-- **Lemma 6 is false**: all its hypotheses and `f(s) = f(1-s)` hold, yet `α₀ = 3/10`. -/
theorem not_lemma6 : ¬ Lemma6 := by
  intro h
  have := h α0 β0 m0 hyp0 symm0 0
  norm_num [α0] at this

/-- **Step (19) is false**: (17) holds with `i = 0, l = 1`, `l ≠ i`. -/
theorem not_step19 : ¬ Step19 := by
  intro h
  have := h α0 β0 m0 hyp0 symm0 0 1 fun s => by rw [g0_reflect]; rfl
  exact one_ne_zero this

/-- (17) itself holds for the counterexample (so the failure is (19), not (17)). -/
theorem div17_0 : Div17 α0 β0 := fun i => ⟨Equiv.swap 0 1 i, fun s => by
  rw [g0_reflect, Equiv.swap_apply_self]⟩

/-! ## Davenport–Heilbronn -/

/-- The paper's Remark says Lemma 6 does not apply to DH because DH's zeros are not all in the
strip. Irrelevant: `not_lemma6` lives in the strip, and DH has a certified off-line zero
`0.808517 + 85.699348 i` in the strip as well. The argument uses only the Hadamard product and the
functional equation, so it fails the DH test. -/
theorem dh_zero_in_strip :
    ∃ s : ℂ, Fayed.DH.DH s = 0 ∧ 0 < s.re ∧ s.re < 1 ∧ s.re ≠ 1 / 2 := by
  obtain ⟨s, hs, hre, _⟩ := Fayed.DH.Zero.offLineZero_holds
  rw [abs_lt] at hre
  norm_num at hre
  obtain ⟨h1, h2⟩ := hre
  refine ⟨s, hs, by linarith, by linarith, fun h => ?_⟩
  rw [h] at h1
  norm_num at h1

end Riemann.Deep.Zhang

#print axioms Riemann.Deep.Zhang.chain
#print axioms Riemann.Deep.Zhang.not_lemma6
#print axioms Riemann.Deep.Zhang.not_step19
#print axioms Riemann.Deep.Zhang.div17_0
#print axioms Riemann.Deep.Zhang.dh_zero_in_strip
