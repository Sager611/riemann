import Mathlib

/-!
# The Davenport–Heilbronn function and its functional equation

Counter-model "Dirichlet series + Riemann-type functional equation, no Euler product, yet zeros
off the critical line" (Davenport–Heilbronn 1936; Titchmarsh, *Theory of the Riemann
zeta-function*, §10.25).

* `χ`: the primitive Dirichlet character mod 5 with `χ 2 = I` (order 4, odd).
* `κ = (√(10 - 2√5) - 2)/(√5 - 1)`, `DH s = (1 - iκ)/2 · L(χ, s) + (1 + iκ)/2 · L(χ̄, s)`.
* `completedDH s = 5^(s/2) · ((1 - iκ)/2 · Λ(χ, s) + (1 + iκ)/2 · Λ(χ̄, s))` with Mathlib's
  `completedLFunction` (`Λ(χ, s) = Γℝ(s + 1) L(χ, s)`, `Γℝ(s) = π^(-s/2) Γ(s/2)`).

Main theorem: `completedDH_one_sub : completedDH (1 - s) = completedDH s` (all `s : ℂ`).
Key input: the root number `ε(χ) = τ(χ)/(i√5) = (1 + iκ)/(1 - iκ)`, computed from the Gauss sum
via `cos(2π/5) = (√5 - 1)/4`, `sin²(2π/5) = (10 + 2√5)/16`, which gives `κ = 2 sin(2π/5) - (√5 + 1)/2`.
-/

open Complex DirichletCharacter

namespace Fayed.DH

/-! ### The character -/

/-- Values of `χ` on `ZMod 5 = {0,1,2,3,4}`: `0, 1, i, -i, -1` (`2` generates `(ZMod 5)ˣ`). -/
noncomputable def chiFun : ZMod 5 → ℂ := ![0, 1, I, -I, -1]

lemma chiFun_mul (a b : ZMod 5) : chiFun (a * b) = chiFun a * chiFun b := by
  fin_cases a <;> fin_cases b <;> simp [chiFun] <;> rfl

lemma chiFun_nonunit (a : ZMod 5) (h : ¬ IsUnit a) : chiFun a = 0 := by
  have : Fact (Nat.Prime 5) := ⟨by norm_num⟩
  rw [isUnit_iff_ne_zero, not_not] at h
  subst h; rfl

/-- The Dirichlet character mod 5 with `χ 2 = I`. -/
noncomputable def χ : DirichletCharacter ℂ 5 where
  toFun := chiFun
  map_one' := rfl
  map_mul' := chiFun_mul
  map_nonunit' := chiFun_nonunit

lemma χ_apply (a : ZMod 5) : χ a = chiFun a := rfl

lemma χ_zero : χ 0 = 0 := rfl
lemma χ_one' : χ 1 = 1 := rfl
lemma χ_two : χ 2 = I := rfl
lemma χ_three : χ 3 = -I := rfl
lemma χ_four : χ 4 = -1 := rfl

lemma χ_neg_one : χ (-1) = -1 := rfl

lemma χ_odd : χ.Odd := χ_neg_one

lemma χ_not_even : ¬ χ.Even := by
  intro h
  have h' : χ (-1) = 1 := h
  rw [χ_neg_one] at h'
  norm_num at h'

lemma χinv_apply (a : ZMod 5) : χ⁻¹ a = (chiFun a)⁻¹ := MulChar.inv_apply_eq_inv' χ a

lemma χinv_neg_one : χ⁻¹ (-1) = -1 := by
  rw [χinv_apply]; norm_num [chiFun]

lemma χinv_not_even : ¬ (χ⁻¹).Even := by
  intro h
  have h' : χ⁻¹ (-1) = 1 := h
  rw [χinv_neg_one] at h'
  norm_num at h'

lemma isUnit_two : IsUnit (2 : ZMod 5) := ⟨⟨2, 3, by decide, by decide⟩, rfl⟩

lemma χ_ne_one : χ ≠ 1 := by
  intro h
  have := congrArg (fun ψ : DirichletCharacter ℂ 5 ↦ ψ 2) h
  simp only [χ_two, MulChar.one_apply isUnit_two] at this
  simp [Complex.ext_iff] at this

/-- **Step 1.** `χ` is primitive: its conductor is `5`. -/
theorem χ_isPrimitive : χ.IsPrimitive := by
  rw [isPrimitive_def]
  rcases (Nat.prime_five).eq_one_or_self_of_dvd _ (conductor_dvd_level χ) with h | h
  · exact absurd (eq_one_iff_conductor_eq_one.mpr h) χ_ne_one
  · exact h

theorem χinv_isPrimitive : (χ⁻¹).IsPrimitive := by
  rw [isPrimitive_def, conductor_inv]; exact χ_isPrimitive

/-! ### Trigonometric constants -/

/-- `S = sin(2π/5)`. -/
noncomputable def S : ℝ := Real.sin (2 * Real.pi / 5)

lemma cos_two_pi_div_five : Real.cos (2 * Real.pi / 5) = (√5 - 1) / 4 := by
  have hr : (√5 : ℝ) ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  rw [show 2 * Real.pi / 5 = 2 * (Real.pi / 5) by ring, Real.cos_two_mul, Real.cos_pi_div_five]
  linear_combination (1 / 8 : ℝ) * hr

lemma S_pos : 0 < S := Real.sin_pos_of_pos_of_lt_pi (by positivity) (by linarith [Real.pi_pos])

lemma S_sq : 16 * S ^ 2 = 10 + 2 * √5 := by
  have hr : (√5 : ℝ) ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  have h := Real.sin_sq_add_cos_sq (2 * Real.pi / 5)
  rw [cos_two_pi_div_five] at h
  unfold S
  linear_combination 16 * h - hr

/-- `κ = (√(10 - 2√5) - 2)/(√5 - 1)`. -/
noncomputable def κ : ℝ := (√(10 - 2 * √5) - 2) / (√5 - 1)

lemma sqrt5_gt_one : 1 < √5 := by
  rw [show (1 : ℝ) = √1 by simp]
  exact Real.sqrt_lt_sqrt (by norm_num) (by norm_num)

lemma κ_eq : κ = 2 * S - (√5 + 1) / 2 := by
  have hr : (√5 : ℝ) ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  have hS := S_sq
  have h1 : 0 < √5 - 1 := by linarith [sqrt5_gt_one]
  have hT : √(10 - 2 * √5) = 2 * S * (√5 - 1) := by
    rw [show 10 - 2 * √5 = (2 * S * (√5 - 1)) ^ 2 by
      linear_combination (-(2 * √5 + 6) / 4) * hr + (-(√5 - 1) ^ 2 / 4) * hS]
    exact Real.sqrt_sq (by have := S_pos; positivity)
  rw [κ, hT]
  field_simp
  linear_combination hr

/-! ### Gauss sums and root numbers -/

/-- `ζ₅ = exp(2πi/5)`. -/
noncomputable def z : ℂ := exp (2 * (Real.pi : ℂ) * I / 5)

lemma z_eq : z = ((√5 : ℝ) - 1) / 4 + (S : ℂ) * I := by
  rw [z, show 2 * (Real.pi : ℂ) * I / 5 = ((2 * Real.pi / 5 : ℝ) : ℂ) * I by push_cast; ring,
    exp_mul_I, ← ofReal_cos, ← ofReal_sin, cos_two_pi_div_five, S]
  push_cast; ring

lemma stdAddChar_five (k : ℕ) : (ZMod.stdAddChar (k : ZMod 5) : ℂ) = z ^ k := by
  have := ZMod.stdAddChar_coe (N := 5) (k : ℤ)
  push_cast at this
  rw [this, z, ← exp_nat_mul]
  congr 1; ring

lemma univ_zmod5 : (Finset.univ : Finset (ZMod 5)) = {0, 1, 2, 3, 4} := by decide

lemma gaussSum_χ : gaussSum χ ZMod.stdAddChar = z + I * z ^ 2 - I * z ^ 3 - z ^ 4 := by
  rw [gaussSum, univ_zmod5]
  simp only [Finset.sum_insert (show (0 : ZMod 5) ∉ ({1, 2, 3, 4} : Finset (ZMod 5)) by decide),
    Finset.sum_insert (show (1 : ZMod 5) ∉ ({2, 3, 4} : Finset (ZMod 5)) by decide),
    Finset.sum_insert (show (2 : ZMod 5) ∉ ({3, 4} : Finset (ZMod 5)) by decide),
    Finset.sum_insert (show (3 : ZMod 5) ∉ ({4} : Finset (ZMod 5)) by decide),
    Finset.sum_singleton]
  have h1 := stdAddChar_five 1
  have h2 := stdAddChar_five 2
  have h3 := stdAddChar_five 3
  have h4 := stdAddChar_five 4
  push_cast at h1 h2 h3 h4
  rw [h1, h2, h3, h4]
  rw [χ_zero, χ_one', χ_two, χ_three, χ_four]
  ring


lemma χinv_zero : χ⁻¹ 0 = 0 := by rw [χinv_apply]; simp [chiFun]
lemma χinv_one : χ⁻¹ 1 = 1 := by rw [χinv_apply, show chiFun 1 = 1 from rfl, inv_one]
lemma χinv_two : χ⁻¹ 2 = -I := by rw [χinv_apply, show chiFun 2 = I from rfl, inv_I]
lemma χinv_three : χ⁻¹ 3 = I := by rw [χinv_apply, show chiFun 3 = -I from rfl, inv_neg, inv_I, neg_neg]
lemma χinv_four : χ⁻¹ 4 = -1 := by rw [χinv_apply, show chiFun 4 = -1 from rfl, inv_neg, inv_one]

lemma gaussSum_χinv : gaussSum χ⁻¹ ZMod.stdAddChar = z - I * z ^ 2 + I * z ^ 3 - z ^ 4 := by
  rw [gaussSum, univ_zmod5]
  simp only [Finset.sum_insert (show (0 : ZMod 5) ∉ ({1, 2, 3, 4} : Finset (ZMod 5)) by decide),
    Finset.sum_insert (show (1 : ZMod 5) ∉ ({2, 3, 4} : Finset (ZMod 5)) by decide),
    Finset.sum_insert (show (2 : ZMod 5) ∉ ({3, 4} : Finset (ZMod 5)) by decide),
    Finset.sum_insert (show (3 : ZMod 5) ∉ ({4} : Finset (ZMod 5)) by decide),
    Finset.sum_singleton]
  have h1 := stdAddChar_five 1
  have h2 := stdAddChar_five 2
  have h3 := stdAddChar_five 3
  have h4 := stdAddChar_five 4
  push_cast at h1 h2 h3 h4
  rw [h1, h2, h3, h4]
  rw [χinv_zero, χinv_one, χinv_two, χinv_three, χinv_four]
  ring

lemma five_cpow_half : ((5 : ℕ) : ℂ) ^ (1 / 2 : ℂ) = ((√5 : ℝ) : ℂ) := by
  rw [Real.sqrt_eq_rpow, ofReal_cpow (by norm_num)]
  push_cast; rfl

/-- Pure algebra behind `ε(χ) = (1 + iκ)/(1 - iκ)`. -/
lemma keyA (S r : ℂ) (hr : r ^ 2 = 5) (hS : 16 * S ^ 2 = 10 + 2 * r) :
    (1 - I * (2 * S - (r + 1) / 2)) / 2 *
      (((r - 1) / 4 + S * I) + I * ((r - 1) / 4 + S * I) ^ 2 - I * ((r - 1) / 4 + S * I) ^ 3
        - ((r - 1) / 4 + S * I) ^ 4) = (1 + I * (2 * S - (r + 1) / 2)) / 2 * (I * r) := by
  have hI : I ^ 2 = -1 := I_sq
  linear_combination (-(-256*S^5*I^3 + 256*S^5*I + 64*S^4*I^3*r - 192*S^4*I^3 - 256*S^4*I^2*r + 384*S^4*I^2 - 64*S^4*I*r + 192*S^4*I + 256*S^4*r - 384*S^4 + 64*S^3*I^3*r + 64*S^3*I^3 + 64*S^3*I^2*r^2 - 192*S^3*I^2*r + 512*S^3*I^2 - 96*S^3*I*r^2 + 256*S^3*I*r - 288*S^3*I - 64*S^3*r^2 + 192*S^3*r - 512*S^3 + 48*S^2*I^2*r^2 - 64*S^2*I^2*r - 112*S^2*I^2 + 24*S^2*I*r^3 - 72*S^2*I*r^2 + 296*S^2*I*r - 376*S^2*I - 16*S^2*r^3 + 48*S^2*r^2 - 80*S^2*r + 432*S^2 + 12*S*I*r^3 - 44*S*I*r^2 - 12*S*I*r + 44*S*I + 4*S*r^4 - 12*S*r^3 + 52*S*r^2 + 44*S*r + 40*S + r^4 - 6*r^3 - 60*r^2 - 58*r - 5)/256) * hI + ((32*S^3*I - 8*S^2*I*r + 24*S^2*I + 32*S^2*r - 48*S^2 - 12*S*I*r^2 + 36*S*I*r - 16*S*I - 8*S*r^2 + 24*S*r - 64*S + 3*I*r^3 - 10*I*r^2 + 35*I*r - 32*I - 2*r^3 + 10*r^2 + 4*r + 24)/512) * hS + (-(-4*S*I*r^2 + 48*S*I*r + 132*S*I - 16*S*r^2 + 80*S*r - 224*S + I*r^3 - 15*I*r^2 - 5*I*r - 133*I + 6*r^2 + 16*r + 66)/1024) * hr

/-- Pure algebra behind `ε(χ̄) = (1 - iκ)/(1 + iκ)`. -/
lemma keyB (S r : ℂ) (hr : r ^ 2 = 5) (hS : 16 * S ^ 2 = 10 + 2 * r) :
    (1 + I * (2 * S - (r + 1) / 2)) / 2 *
      (((r - 1) / 4 + S * I) - I * ((r - 1) / 4 + S * I) ^ 2 + I * ((r - 1) / 4 + S * I) ^ 3
        - ((r - 1) / 4 + S * I) ^ 4) = (1 - I * (2 * S - (r + 1) / 2)) / 2 * (I * r) := by
  have hI : I ^ 2 = -1 := I_sq
  linear_combination ((-256*S^5*I^3 + 256*S^5*I + 64*S^4*I^3*r + 320*S^4*I^3 - 256*S^4*I^2*r + 128*S^4*I^2 - 64*S^4*I*r - 320*S^4*I + 256*S^4*r - 128*S^4 - 64*S^3*I^3*r - 64*S^3*I^3 + 64*S^3*I^2*r^2 + 192*S^3*I^2*r - 384*S^3*I^2 - 96*S^3*I*r^2 + 128*S^3*I*r + 96*S^3*I - 64*S^3*r^2 - 192*S^3*r + 384*S^3 - 48*S^2*I^2*r^2 + 64*S^2*I^2*r + 112*S^2*I^2 + 24*S^2*I*r^3 + 24*S^2*I*r^2 - 152*S^2*I*r - 24*S^2*I - 16*S^2*r^3 + 48*S^2*r^2 - 16*S^2*r + 112*S^2 - 12*S*I*r^3 + 44*S*I*r^2 + 12*S*I*r - 44*S*I + 4*S*r^4 - 4*S*r^3 - 4*S*r^2 + 132*S*r - r^4 + 6*r^3 - 68*r^2 - 70*r + 5)/256) * hI + (-(32*S^3*I - 8*S^2*I*r - 40*S^2*I + 32*S^2*r - 16*S^2 - 12*S*I*r^2 + 20*S*I*r + 32*S*I - 8*S*r^2 - 24*S*r + 48*S + 3*I*r^3 + 2*I*r^2 - 29*I*r - 28*I - 2*r^3 + 10*r^2 + 16*r + 4)/512) * hS + ((-4*S*I*r^2 + 80*S*I*r + 36*S*I - 16*S*r^2 + 48*S*r + 192*S + I*r^3 - 15*I*r^2 - 53*I*r - 117*I + 10*r^2 - 16*r + 46)/1024) * hr

/-! ### The Davenport–Heilbronn function -/

/-- Weight of `L(χ, s)`. -/
noncomputable def a : ℂ := (1 - I * κ) / 2
/-- Weight of `L(χ̄, s)`. -/
noncomputable def b : ℂ := (1 + I * κ) / 2

lemma hr_c : ((√5 : ℝ) : ℂ) ^ 2 = 5 := by
  rw [← ofReal_pow, Real.sq_sqrt (by norm_num)]; push_cast; rfl

lemma hS_c : 16 * (S : ℂ) ^ 2 = 10 + 2 * ((√5 : ℝ) : ℂ) := by
  have := congrArg (fun x : ℝ ↦ (x : ℂ)) S_sq
  push_cast at this; exact this

lemma κ_c : (κ : ℂ) = 2 * (S : ℂ) - (((√5 : ℝ) : ℂ) + 1) / 2 := by
  rw [κ_eq]; push_cast; ring

lemma sqrt5_ne : ((√5 : ℝ) : ℂ) ≠ 0 := by
  rw [ofReal_ne_zero]; positivity

/-- **Step 3a.** `(1 - iκ)/2 · ε(χ) = (1 + iκ)/2`. -/
lemma a_mul_rootNumber : a * rootNumber χ = b := by
  rw [rootNumber]
  simp only [χ_not_even, ↓reduceIte]
  rw [pow_one, gaussSum_χ, five_cpow_half, div_div,
    mul_div_assoc', div_eq_iff (mul_ne_zero I_ne_zero sqrt5_ne), a, b, κ_c, z_eq]
  exact keyA _ _ hr_c hS_c

/-- **Step 3b.** `(1 + iκ)/2 · ε(χ̄) = (1 - iκ)/2`. -/
lemma b_mul_rootNumber_inv : b * rootNumber χ⁻¹ = a := by
  rw [rootNumber]
  simp only [χinv_not_even, ↓reduceIte]
  rw [pow_one, gaussSum_χinv, five_cpow_half, div_div,
    mul_div_assoc', div_eq_iff (mul_ne_zero I_ne_zero sqrt5_ne), a, b, κ_c, z_eq]
  exact keyB _ _ hr_c hS_c

set_option linter.dupNamespace false in
/-- **Step 2.** The Davenport–Heilbronn function. -/
noncomputable def DH (s : ℂ) : ℂ := a * LFunction χ s + b * LFunction χ⁻¹ s

/-- Completed DH function `5^(s/2) (a Λ(χ,s) + b Λ(χ̄,s))` (entire). -/
noncomputable def completedDH (s : ℂ) : ℂ :=
  (5 : ℂ) ^ (s / 2) * (a * completedLFunction χ s + b * completedLFunction χ⁻¹ s)

/-- **Step 3. Functional equation** of the Davenport–Heilbronn function, for all `s`. -/
theorem completedDH_one_sub (s : ℂ) : completedDH (1 - s) = completedDH s := by
  have h1 := χ_isPrimitive.completedLFunction_one_sub s
  have h2 := χinv_isPrimitive.completedLFunction_one_sub s
  rw [inv_inv] at h2
  simp only [Nat.cast_ofNat] at h1 h2
  have hpow : (5 : ℂ) ^ ((1 - s) / 2) * (5 : ℂ) ^ (s - 1 / 2) = (5 : ℂ) ^ (s / 2) := by
    rw [← cpow_add _ _ (by norm_num)]; congr 1; ring
  unfold completedDH
  rw [h1, h2]
  calc _ = ((5 : ℂ) ^ ((1 - s) / 2) * (5 : ℂ) ^ (s - 1 / 2)) *
        ((a * rootNumber χ) * completedLFunction χ⁻¹ s +
          (b * rootNumber χ⁻¹) * completedLFunction χ s) := by ring
    _ = _ := by rw [hpow, a_mul_rootNumber, b_mul_rootNumber_inv]; ring

lemma χinv_odd : (χ⁻¹).Odd := χinv_neg_one

lemma completed_eq_of_odd {ψ : DirichletCharacter ℂ 5} (hψ : ψ.Odd) {s : ℂ}
    (h : Gammaℝ (s + 1) ≠ 0) : completedLFunction ψ s = Gammaℝ (s + 1) * LFunction ψ s := by
  rw [LFunction_eq_completed_div_gammaFactor ψ s (Or.inr (by norm_num)), hψ.gammaFactor_def]
  field_simp

/-- Away from the poles of the Gamma factor, `completedDH = 5^(s/2) Γℝ(s+1) DH`. -/
lemma completedDH_eq {s : ℂ} (h : Gammaℝ (s + 1) ≠ 0) :
    completedDH s = (5 : ℂ) ^ (s / 2) * Gammaℝ (s + 1) * DH s := by
  rw [completedDH, DH, completed_eq_of_odd χ_odd h, completed_eq_of_odd χinv_odd h]
  ring

/-- Functional equation in Gamma-factor form (where neither Gamma factor has a pole):
`5^((1-s)/2) Γℝ(2-s) DH(1-s) = 5^(s/2) Γℝ(s+1) DH(s)`, `Γℝ(s) = π^(-s/2) Γ(s/2)`. -/
theorem DH_functional_equation {s : ℂ} (h1 : Gammaℝ (s + 1) ≠ 0) (h2 : Gammaℝ (2 - s) ≠ 0) :
    (5 : ℂ) ^ ((1 - s) / 2) * Gammaℝ (2 - s) * DH (1 - s) =
      (5 : ℂ) ^ (s / 2) * Gammaℝ (s + 1) * DH s := by
  have h := completedDH_eq (s := 1 - s) (by rwa [show 1 - s + 1 = 2 - s by ring])
  rw [show 1 - s + 1 = 2 - s by ring] at h
  rw [← h, ← completedDH_eq h1, completedDH_one_sub]

/-- The classical normalisation `Λ_DH(s) = (5/π)^(s/2) Γ((s+1)/2) DH(s)`. -/
noncomputable def classicalDH (s : ℂ) : ℂ :=
  ((5 / Real.pi : ℝ) : ℂ) ^ (s / 2) * Gamma ((s + 1) / 2) * DH s

lemma classicalDH_eq (s : ℂ) :
    classicalDH s = ((√Real.pi : ℝ) : ℂ) * ((5 : ℂ) ^ (s / 2) * Gammaℝ (s + 1) * DH s) := by
  have hπ : (Real.pi : ℂ) ≠ 0 := ofReal_ne_zero.mpr Real.pi_ne_zero
  have e1 : ((5 / Real.pi : ℝ) : ℂ) ^ (s / 2) =
      (5 : ℂ) ^ (s / 2) * ((Real.pi⁻¹ : ℝ) : ℂ) ^ (s / 2) := by
    rw [div_eq_mul_inv, ofReal_mul, mul_cpow_ofReal_nonneg (by norm_num) (by positivity)]
    norm_num
  have e2 : ((Real.pi⁻¹ : ℝ) : ℂ) ^ (s / 2) = (Real.pi : ℂ) ^ (-(s / 2)) := by
    have harg : (Real.pi : ℂ).arg ≠ Real.pi := by
      rw [arg_ofReal_of_nonneg Real.pi_pos.le]; exact Real.pi_ne_zero.symm
    rw [ofReal_inv, inv_cpow _ _ harg, cpow_neg]
  have e3 : ((√Real.pi : ℝ) : ℂ) = (Real.pi : ℂ) ^ (1 / 2 : ℂ) := by
    rw [Real.sqrt_eq_rpow, ofReal_cpow Real.pi_pos.le]; push_cast; rfl
  have e4 : (Real.pi : ℂ) ^ (1 / 2 : ℂ) * (Real.pi : ℂ) ^ (-(s + 1) / 2) =
      (Real.pi : ℂ) ^ (-(s / 2)) := by
    rw [← cpow_add _ _ hπ]; congr 1; ring
  rw [classicalDH, e1, e2, e3, Gammaℝ_def, ← e4]
  ring

/-- Functional equation in the classical form `Λ_DH(1 - s) = Λ_DH(s)`,
`Λ_DH(s) = (5/π)^(s/2) Γ((s+1)/2) DH(s)`, away from the poles of `Γ((s+1)/2)`, `Γ((2-s)/2)`. -/
theorem classicalDH_one_sub {s : ℂ} (h1 : Gammaℝ (s + 1) ≠ 0) (h2 : Gammaℝ (2 - s) ≠ 0) :
    classicalDH (1 - s) = classicalDH s := by
  rw [classicalDH_eq, classicalDH_eq, show 1 - s + 1 = 2 - s by ring,
    DH_functional_equation h1 h2]

/-! ### Step 4: real Dirichlet coefficients -/

/-- Dirichlet coefficients of DH. -/
noncomputable def coeff (n : ℕ) : ℂ := a * χ n + b * χ⁻¹ n

lemma zmod5_cases (j : ZMod 5) : j = 0 ∨ j = 1 ∨ j = 2 ∨ j = 3 ∨ j = 4 := by
  revert j; decide

/-- The coefficients are `1, κ, -κ, -1, 0` for `n ≡ 1, 2, 3, 4, 0 (mod 5)`. -/
lemma coeff_values (n : ℕ) :
    ((n : ZMod 5) = 0 ∧ coeff n = 0) ∨ ((n : ZMod 5) = 1 ∧ coeff n = 1) ∨
    ((n : ZMod 5) = 2 ∧ coeff n = κ) ∨ ((n : ZMod 5) = 3 ∧ coeff n = -κ) ∨
    ((n : ZMod 5) = 4 ∧ coeff n = -1) := by
  unfold coeff
  rcases zmod5_cases (n : ZMod 5) with h | h | h | h | h <;> rw [h]
  · left; simp [χ_zero, χinv_zero]
  · right; left; simp only [χ_one', χinv_one, a, b]; ring_nf; simp
  · right; right; left; simp only [χ_two, χinv_two, a, b]; ring_nf; simp [I_sq]
  · right; right; right; left; simp only [χ_three, χinv_three, a, b]; ring_nf; simp [I_sq]
  · right; right; right; right; simp only [χ_four, χinv_four, a, b]; ring_nf; simp

/-- **Step 4.** The Dirichlet coefficients of DH are real. -/
theorem coeff_real (n : ℕ) : coeff n = ((coeff n).re : ℂ) := by
  rcases coeff_values n with ⟨-, h⟩ | ⟨-, h⟩ | ⟨-, h⟩ | ⟨-, h⟩ | ⟨-, h⟩ <;> rw [h] <;> simp

theorem coeff_im (n : ℕ) : (coeff n).im = 0 := by
  rw [coeff_real n]; simp

/-- For `Re s > 1`, DH is the Dirichlet series with the (real) coefficients `coeff`. -/
theorem DH_eq_LSeries {s : ℂ} (hs : 1 < s.re) : DH s = LSeries coeff s := by
  rw [DH, LFunction_eq_LSeries _ hs, LFunction_eq_LSeries _ hs, ← LSeries_smul, ← LSeries_smul,
    ← LSeries_add ((LSeriesSummable_of_one_lt_re χ hs).smul a)
      ((LSeriesSummable_of_one_lt_re χ⁻¹ hs).smul b)]
  congr 1

/-! ### Off-line zeros (numerical, cited only) -/

/-- **Cited numerical fact, not proved here.** DH has a zero near `0.808517 + 85.699348 i`
(R. Spira, *Some zeros of the Titchmarsh counterexample*, Math. Comp. 63 (1994); see also
Balanzario–Sánchez-Ortiz, Math. Comp. 76 (2007)), i.e. off the critical line. -/
def OffLineZeroClaim : Prop :=
  ∃ s : ℂ, DH s = 0 ∧ |s.re - 0.808517| < 1e-5 ∧ |s.im - 85.699348| < 1e-5

end Fayed.DH

#print axioms Fayed.DH.completedDH_one_sub
#print axioms Fayed.DH.DH_functional_equation
#print axioms Fayed.DH.classicalDH_one_sub
#print axioms Fayed.DH.coeff_real
