module

public import Riemann.Deep.LeanLi.IbpS

/-!
# LeanLi: ingredients of Theorem 4 (smoothed mean `7/4`)
-/

@[expose] public section

open Real Filter Set Finset MeasureTheory Topology

namespace LeanLi

/-- A height below the first ordinate. -/
lemma exists_T0 : ∃ T₀ : ℝ, 0 < T₀ ∧ T₀ ≤ 1 ∧ zset T₀ = ∅ := by
  by_cases h : (zset 1).Nonempty
  · obtain ⟨ρ₀, hρ₀, hmin⟩ := (zset 1).exists_min_image Complex.im h
    have hpos : 0 < ρ₀.im := (mem_zset.1 hρ₀).2.1
    have hle : ρ₀.im ≤ 1 := (mem_zset.1 hρ₀).2.2
    refine ⟨ρ₀.im / 2, by positivity, by linarith, ?_⟩
    ext ρ
    simp only [Finset.notMem_empty, iff_false]
    intro hρ
    have h1 := mem_zset.1 hρ
    have hρ1 : ρ ∈ zset 1 := mem_zset.2 ⟨h1.1, h1.2.1, by linarith [h1.2.2]⟩
    have := hmin ρ hρ1
    linarith [h1.2.2]
  · exact ⟨1, one_pos, le_rfl, Finset.not_nonempty_iff_eq_empty.1 h⟩

lemma ZK_eq_empty {T : ℝ} (h : zset T = ∅) : ZK T = ∅ := by
  have h1 := card_ZK T
  rw [zeroCount_eq, h, Finset.sum_empty] at h1
  exact Finset.card_eq_zero.1 (by exact_mod_cast h1)

lemma zeroCount_eq_zero {T : ℝ} (h : zset T = ∅) : (zeroCount T : ℝ) = 0 := by
  rw [zeroCount_eq, h, Finset.sum_empty]

lemma Ej_re (φ : ℝ → ℝ) (B N j : ℕ) (x : ℝ) :
    (Ej φ B N j x).re = ∑ n ∈ range (B * N + 1), φ (n / N) * (n : ℝ) ^ j * Real.cos (n * x) := by
  unfold Ej
  rw [Complex.re_sum]
  refine Finset.sum_congr rfl fun n _ => ?_
  have e : Complex.I * (((n : ℝ) : ℂ) * (x : ℂ)) = (((n : ℝ) * x : ℝ) : ℂ) * Complex.I := by push_cast; ring
  rw [Complex.re_ofReal_mul, e, Complex.exp_ofReal_mul_I_re]

lemma Ej_im (φ : ℝ → ℝ) (B N j : ℕ) (x : ℝ) :
    (Ej φ B N j x).im = ∑ n ∈ range (B * N + 1), φ (n / N) * (n : ℝ) ^ j * Real.sin (n * x) := by
  unfold Ej
  rw [Complex.im_sum]
  refine Finset.sum_congr rfl fun n _ => ?_
  have e : Complex.I * (((n : ℝ) : ℂ) * (x : ℂ)) = (((n : ℝ) * x : ℝ) : ℂ) * Complex.I := by push_cast; ring
  rw [Complex.im_ofReal_mul, e, Complex.exp_ofReal_mul_I_im]

/-- `F₀(t) = Σ φ(n/N) cos nθ(t)`, `F₁ = Σ φ n sin nθ`, `F₂ = Σ φ n² cos nθ`. -/
noncomputable def F0 (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : ℝ :=
  ∑ n ∈ range (B * N + 1), φ (n / N) * Real.cos (n * liTheta t)
noncomputable def F1 (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : ℝ :=
  ∑ n ∈ range (B * N + 1), φ (n / N) * n * Real.sin (n * liTheta t)
noncomputable def F2 (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : ℝ :=
  ∑ n ∈ range (B * N + 1), φ (n / N) * (n : ℝ) ^ 2 * Real.cos (n * liTheta t)

lemma F0_eq (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : F0 φ B N t = (Ej φ B N 0 (liTheta t)).re := by
  rw [Ej_re]; unfold F0; simp

lemma F1_eq (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : F1 φ B N t = (Ej φ B N 1 (liTheta t)).im := by
  rw [Ej_im]; unfold F1; simp

lemma F2_eq (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : F2 φ B N t = (Ej φ B N 2 (liTheta t)).re := by
  rw [Ej_re]; unfold F2; simp

/-- The kernel `G = 2F₁/(t² + ¼)` and its derivative. -/
noncomputable def Gf (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : ℝ := 2 / (t ^ 2 + 1 / 4) * F1 φ B N t
noncomputable def Gd (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : ℝ :=
  -(4 * t) / (t ^ 2 + 1 / 4) ^ 2 * F1 φ B N t - 2 / (t ^ 2 + 1 / 4) ^ 2 * F2 φ B N t

lemma hasDerivAt_F1 (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) :
    HasDerivAt (F1 φ B N) (-(1 / (t ^ 2 + 1 / 4)) * F2 φ B N t) t := by
  unfold F1 F2
  rw [Finset.mul_sum]
  apply HasDerivAt.fun_sum
  intro n _
  have h := (((hasDerivAt_liTheta t).const_mul (n : ℝ)).sin).const_mul (φ (n / N) * n)
  convert h using 1
  ring

lemma hasDerivAt_Gf (φ : ℝ → ℝ) (B N : ℕ) (t : ℝ) : HasDerivAt (Gf φ B N) (Gd φ B N t) t := by
  have hw : t ^ 2 + 1 / 4 ≠ 0 := by positivity
  have h1 : HasDerivAt (fun t : ℝ => 2 / (t ^ 2 + 1 / 4)) (-(4 * t) / (t ^ 2 + 1 / 4) ^ 2) t := by
    have := ((hasDerivAt_pow 2 t).add_const (1 / 4 : ℝ)).inv hw
    have h2 := this.const_mul (2 : ℝ)
    convert h2 using 1
    · funext s; simp [div_eq_mul_inv]
    · field_simp; push_cast; ring
  have h := h1.mul (hasDerivAt_F1 φ B N t)
  unfold Gf Gd
  convert h using 1
  field_simp
  ring

lemma continuous_F1 (φ : ℝ → ℝ) (B N : ℕ) : Continuous (F1 φ B N) := by
  unfold F1; have := continuous_liTheta; fun_prop

lemma continuous_F2 (φ : ℝ → ℝ) (B N : ℕ) : Continuous (F2 φ B N) := by
  unfold F2; have := continuous_liTheta; fun_prop

lemma continuous_F0 (φ : ℝ → ℝ) (B N : ℕ) : Continuous (F0 φ B N) := by
  unfold F0; have := continuous_liTheta; fun_prop

lemma continuous_Gd (φ : ℝ → ℝ) (B N : ℕ) : Continuous (Gd φ B N) := by
  unfold Gd
  have h1 := continuous_F1 φ B N
  have h2 := continuous_F2 φ B N
  have hw : ∀ t : ℝ, (t ^ 2 + 1 / 4) ^ 2 ≠ 0 := fun t => by positivity
  exact ((continuous_const.mul continuous_id).neg.div (by fun_prop) hw |>.mul h1).sub
    ((continuous_const.div (by fun_prop) hw).mul h2)

section
variable {φ : ℝ → ℝ} {B : ℕ} (h : TestFn φ B)
include h

/-- `|F₀(t)| ≤ C (2t + 1)`. -/
lemma F0_bound : ∃ C : ℝ, 0 ≤ C ∧ ∀ N : ℕ, 1 ≤ N → ∀ t : ℝ, 0 ≤ t → |F0 φ B N t| ≤ C * (2 * t + 1) := by
  obtain ⟨C, hC0, hC⟩ := Ej_theta h 0
  refine ⟨C, hC0, fun N hN t ht => ?_⟩
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  rw [F0_eq]
  refine (Complex.abs_re_le_norm _).trans ((hC N hN t ht).trans ?_)
  have hu : 1 ≤ 2 * t + 1 := by linarith
  have hv : 0 < 2 * t + 1 + N := by linarith
  rw [div_le_iff₀ (by positivity), pow_one]
  have h1 : (2 * t + 1) ^ 3 ≤ (2 * t + 1 + N) ^ 3 := pow_le_pow_left₀ (by linarith) (by linarith) 3
  have h2 : (N : ℝ) ≤ 2 * t + 1 + N := by linarith
  have h3 : (N : ℝ) * (2 * t + 1) ^ 3 ≤ (2 * t + 1 + N) * (2 * t + 1 + N) ^ 3 :=
    mul_le_mul h2 h1 (by positivity) (by positivity)
  have e1 : C * (N : ℝ) * (2 * t + 1) ^ 4 = C * (2 * t + 1) * ((N : ℝ) * (2 * t + 1) ^ 3) := by ring
  have e2 : C * (2 * t + 1) * (2 * t + 1 + N) ^ 4 = C * (2 * t + 1) * ((2 * t + 1 + N) * (2 * t + 1 + N) ^ 3) := by ring
  rw [e1, e2]
  exact mul_le_mul_of_nonneg_left h3 (by positivity)

/-- `|G(t)| ≤ 16 C`. -/
lemma Gf_bound : ∃ C : ℝ, 0 ≤ C ∧ ∀ N : ℕ, 1 ≤ N → ∀ t : ℝ, 0 ≤ t → |Gf φ B N t| ≤ C := by
  obtain ⟨C, hC0, hC⟩ := Ej_theta h 1
  refine ⟨16 * C, by positivity, fun N hN t ht => ?_⟩
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have hw : 0 < t ^ 2 + 1 / 4 := by positivity
  have hF := (Complex.abs_im_le_norm _).trans (hC N hN t ht)
  rw [← F1_eq] at hF
  unfold Gf
  rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 2 / (t ^ 2 + 1 / 4))]
  have hv : 0 < 2 * t + 1 + N := by linarith
  have hu2 : (2 * t + 1) ^ 2 ≤ 8 * (t ^ 2 + 1 / 4) := by nlinarith [sq_nonneg (2 * t - 1)]
  have hNv : (N : ℝ) ≤ 2 * t + 1 + N := by linarith
  have huv : 2 * t + 1 ≤ 2 * t + 1 + N := by linarith
  calc 2 / (t ^ 2 + 1 / 4) * |F1 φ B N t| ≤ 2 / (t ^ 2 + 1 / 4) * (C * (N : ℝ) ^ (1 + 1) * (2 * t + 1) ^ 4 / (2 * t + 1 + N) ^ 4) :=
        mul_le_mul_of_nonneg_left hF (by positivity)
    _ ≤ 16 * C := by
        rw [div_mul_div_comm, div_le_iff₀ (by positivity)]
        have h1 : (N : ℝ) ^ 2 * (2 * t + 1) ^ 2 ≤ (2 * t + 1 + N) ^ 4 := by
          have := mul_le_mul (pow_le_pow_left₀ (Nat.cast_nonneg _) hNv 2) (pow_le_pow_left₀ (by linarith) huv 2)
            (by positivity) (by positivity)
          calc (N : ℝ) ^ 2 * (2 * t + 1) ^ 2 ≤ (2 * t + 1 + N) ^ 2 * (2 * t + 1 + N) ^ 2 := this
            _ = (2 * t + 1 + N) ^ 4 := by ring
        have h2 : 2 * (C * (N : ℝ) ^ (1 + 1) * (2 * t + 1) ^ 4) ≤ 16 * C * (t ^ 2 + 1 / 4) * ((N : ℝ) ^ 2 * (2 * t + 1) ^ 2) := by
          have : (2 * t + 1) ^ 4 = (2 * t + 1) ^ 2 * (2 * t + 1) ^ 2 := by ring
          rw [this]
          have := mul_le_mul_of_nonneg_left hu2 (show 0 ≤ 2 * C * (N : ℝ) ^ 2 * (2 * t + 1) ^ 2 by positivity)
          nlinarith
        calc 2 * (C * (N : ℝ) ^ (1 + 1) * (2 * t + 1) ^ 4) ≤ 16 * C * (t ^ 2 + 1 / 4) * ((N : ℝ) ^ 2 * (2 * t + 1) ^ 2) := h2
          _ ≤ 16 * C * ((t ^ 2 + 1 / 4) * (2 * t + 1 + N) ^ 4) := by
              rw [mul_assoc (16 * C)]; exact mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left h1 hw.le) (by positivity)

/-- `|G'(t)| ≤ C N²/(t + N)³`. -/
lemma Gd_bound : ∃ C : ℝ, 0 ≤ C ∧ ∀ N : ℕ, 1 ≤ N → ∀ t : ℝ, 0 ≤ t →
    |Gd φ B N t| ≤ C * (N : ℝ) ^ 2 / (t + N) ^ 3 := by
  obtain ⟨C₁, hC₁0, hC₁⟩ := Ej_theta h 1
  obtain ⟨C₂, hC₂0, hC₂⟩ := Ej_theta h 2
  refine ⟨64 * (4 * C₁ + 2 * C₂), by positivity, fun N hN t ht => ?_⟩
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have hw : 0 < t ^ 2 + 1 / 4 := by positivity
  have hF1 := (Complex.abs_im_le_norm _).trans (hC₁ N hN t ht)
  rw [← F1_eq] at hF1
  have hF2 := (Complex.abs_re_le_norm _).trans (hC₂ N hN t ht)
  rw [← F2_eq] at hF2
  set u := 2 * t + 1 with hu
  set v := 2 * t + 1 + (N : ℝ) with hv
  have hv0 : 0 < v := by linarith
  have hu0 : 0 < u := by linarith
  have hu4 : u ^ 4 ≤ 64 * (t ^ 2 + 1 / 4) ^ 2 := by
    have hu2 : u ^ 2 ≤ 8 * (t ^ 2 + 1 / 4) := by rw [hu]; nlinarith [sq_nonneg (2 * t - 1)]
    have := pow_le_pow_left₀ (by positivity) hu2 2
    calc u ^ 4 = (u ^ 2) ^ 2 := by ring
      _ ≤ (8 * (t ^ 2 + 1 / 4)) ^ 2 := this
      _ = 64 * (t ^ 2 + 1 / 4) ^ 2 := by ring
  have htN : 0 < t + N := by linarith
  have htv : t + N ≤ v := by rw [hv]; linarith
  -- |G'| ≤ (4t|F₁| + 2|F₂|)/(t²+¼)²
  have hG : |Gd φ B N t| ≤ (4 * t * |F1 φ B N t| + 2 * |F2 φ B N t|) / (t ^ 2 + 1 / 4) ^ 2 := by
    unfold Gd
    have e : -(4 * t) / (t ^ 2 + 1 / 4) ^ 2 * F1 φ B N t - 2 / (t ^ 2 + 1 / 4) ^ 2 * F2 φ B N t =
        -((4 * t * F1 φ B N t + 2 * F2 φ B N t) / (t ^ 2 + 1 / 4) ^ 2) := by ring
    rw [e, abs_neg, abs_div, abs_of_pos (by positivity : (0 : ℝ) < (t ^ 2 + 1 / 4) ^ 2)]
    apply div_le_div_of_nonneg_right _ (by positivity)
    calc |4 * t * F1 φ B N t + 2 * F2 φ B N t| ≤ |4 * t * F1 φ B N t| + |2 * F2 φ B N t| := abs_add_le _ _
      _ = 4 * t * |F1 φ B N t| + 2 * |F2 φ B N t| := by
          rw [abs_mul, abs_mul, abs_mul, abs_of_nonneg ht, abs_of_pos (by norm_num : (0 : ℝ) < 4),
            abs_of_pos (by norm_num : (0 : ℝ) < 2)]
  refine hG.trans ?_
  have hnum : 4 * t * |F1 φ B N t| + 2 * |F2 φ B N t| ≤
      (4 * t * C₁ * (N : ℝ) ^ 2 + 2 * C₂ * (N : ℝ) ^ 3) * u ^ 4 / v ^ 4 := by
    have e1 : 4 * t * |F1 φ B N t| ≤ 4 * t * (C₁ * (N : ℝ) ^ (1 + 1) * u ^ 4 / v ^ 4) :=
      mul_le_mul_of_nonneg_left hF1 (by positivity)
    have e2 : 2 * |F2 φ B N t| ≤ 2 * (C₂ * (N : ℝ) ^ (2 + 1) * u ^ 4 / v ^ 4) :=
      mul_le_mul_of_nonneg_left hF2 (by positivity)
    have e3 : 4 * t * (C₁ * (N : ℝ) ^ (1 + 1) * u ^ 4 / v ^ 4) + 2 * (C₂ * (N : ℝ) ^ (2 + 1) * u ^ 4 / v ^ 4) =
        (4 * t * C₁ * (N : ℝ) ^ 2 + 2 * C₂ * (N : ℝ) ^ 3) * u ^ 4 / v ^ 4 := by ring
    linarith
  have hA : 4 * t * C₁ * (N : ℝ) ^ 2 + 2 * C₂ * (N : ℝ) ^ 3 ≤ (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (t + N) := by
    nlinarith [sq_nonneg (N : ℝ), mul_nonneg hC₁0 (sq_nonneg (N : ℝ)), mul_nonneg hC₂0 (sq_nonneg (N : ℝ))]
  rw [div_le_div_iff₀ (by positivity) (by positivity)]
  have hv4 : (t + N) ^ 4 ≤ v ^ 4 := pow_le_pow_left₀ htN.le htv 4
  -- (num) * (t+N)^3 ≤ C N² (t²+¼)²
  have h1 : (4 * t * |F1 φ B N t| + 2 * |F2 φ B N t|) ≤ (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (t + N) * u ^ 4 / v ^ 4 := by
    refine hnum.trans ?_
    apply div_le_div_of_nonneg_right _ (by positivity)
    exact mul_le_mul_of_nonneg_right hA (by positivity)
  have h2 : (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (t + N) * u ^ 4 / v ^ 4 * (t + N) ^ 3 ≤
      64 * (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (t ^ 2 + 1 / 4) ^ 2 := by
    rw [div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
    have hK : 0 ≤ (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 := by positivity
    calc (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (t + N) * u ^ 4 * (t + N) ^ 3 =
          (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * u ^ 4 * (t + N) ^ 4 := by ring
      _ ≤ (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (64 * (t ^ 2 + 1 / 4) ^ 2) * v ^ 4 := by
          apply mul_le_mul (mul_le_mul_of_nonneg_left hu4 hK) hv4 (by positivity) (by positivity)
      _ = 64 * (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (t ^ 2 + 1 / 4) ^ 2 * v ^ 4 := by ring
  calc (4 * t * |F1 φ B N t| + 2 * |F2 φ B N t|) * (t + N) ^ 3 ≤
        (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (t + N) * u ^ 4 / v ^ 4 * (t + N) ^ 3 :=
        mul_le_mul_of_nonneg_right h1 (by positivity)
    _ ≤ 64 * (4 * C₁ + 2 * C₂) * (N : ℝ) ^ 2 * (t ^ 2 + 1 / 4) ^ 2 := h2

end

/-- `∫_a^b N²/(t+N)³ ≤ 1/2`. -/
lemma integral_Nsq_div_cube_le {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) {N : ℝ} (hN : 1 ≤ N) :
    ∫ t in a..b, N ^ 2 / (t + N) ^ 3 ≤ 1 / 2 := by
  have hd : ∀ t ∈ Set.uIcc a b, HasDerivAt (fun t => -(N ^ 2 / 2) * ((t + N) ^ 2)⁻¹) (N ^ 2 / (t + N) ^ 3) t := by
    intro t ht
    rw [Set.uIcc_of_le hab] at ht
    have htN : t + N ≠ 0 := by linarith [ht.1]
    have h0 : HasDerivAt (fun s : ℝ => (s + N) ^ 2) (2 * (t + N)) t := by
      have := ((hasDerivAt_id' t).add_const N).pow 2
      convert this using 1; push_cast; ring
    have h3 := (h0.inv (pow_ne_zero 2 htN)).const_mul (-(N ^ 2 / 2))
    convert h3 using 1
    field_simp
  have hc : ContinuousOn (fun t => N ^ 2 / (t + N) ^ 3) (Set.uIcc a b) := by
    intro t ht
    rw [Set.uIcc_of_le hab] at ht
    have hp : 0 < t + N := by linarith [ht.1]
    have hne : (t + N) ^ 3 ≠ 0 := pow_ne_zero 3 hp.ne'
    exact (continuousAt_const.div ((continuousAt_id.add continuousAt_const).pow 3) hne).continuousWithinAt
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd (hc.intervalIntegrable)]
  have h1 : 0 ≤ N ^ 2 / 2 * ((b + N) ^ 2)⁻¹ := by
    have : 0 < b + N := by linarith
    positivity
  have h2 : N ^ 2 / 2 * ((a + N) ^ 2)⁻¹ ≤ 1 / 2 := by
    have : 0 < a + N := by linarith
    rw [← div_eq_mul_inv, div_le_iff₀ (by positivity)]
    nlinarith
  linarith

end LeanLi
