import Riemann.Deep.GL2Ceil.Lean.Defs
import Mathlib.Analysis.SpecialFunctions.Trigonometric.DerivHyp
import Mathlib.Analysis.Real.Sqrt

/-!
# GL2Ceil (gl2lean): the inner minimum over P, monotonicity in λ and in ξ

DERIVATION §1.4: for a, b ≥ 0 with b = aλ², every C¹ P with P(0) = 0, P(1) = 1 has
a C_P + b B_P ≥ a C_λ + b B_λ, where P_λ = sinh(λx)/sinh λ (P_0 = x), C_λ = ∫ P_λ'², B_λ = ∫ P_λ²
(proof: P = P_λ + η, η(0) = η(1) = 0, the cross term vanishes by P_λ'' = λ² P_λ).
§1.5: B_λ is non-increasing and C_λ non-decreasing on λ ≥ 0 (exchange argument).
§1.7: C_λ + λ² B_λ = P_λ'(1) = λ coth λ, so the optimal value √(W₀W₁) coth(ξ√(W₁/W₀)) decreases in ξ.
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

/-- Inner minimum, general form: P* with P*' = Pd, Pd' = t P*, P*(0) = 0, P*(1) = 1. -/
theorem quad_ge {P Ps Pd : ℝ → ℝ} {t a b : ℝ} (hP : PAdm P) (hs0 : Ps 0 = 0) (hs1 : Ps 1 = 1)
    (hsd : ∀ x, HasDerivAt Ps (Pd x) x) (hdd : ∀ x, HasDerivAt Pd (t * Ps x) x)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : b = a * t) :
    a * (∫ x in (0:ℝ)..1, Pd x ^ 2) + b * (∫ x in (0:ℝ)..1, Ps x ^ 2) ≤ a * CP P + b * BP P := by
  obtain ⟨hPd, hPc⟩ := contDiff_one_iff_deriv.mp hP.smooth
  have hP0 : Continuous P := hPd.continuous
  have hsc : Continuous Ps := continuous_iff_continuousAt.2 fun x => (hsd x).continuousAt
  have hdc : Continuous Pd := continuous_iff_continuousAt.2 fun x => (hdd x).continuousAt
  have hibp : ∫ x in (0:ℝ)..1, Pd x * (deriv P x - Pd x)
      = Pd 1 * (P 1 - Ps 1) - Pd 0 * (P 0 - Ps 0) - ∫ x in (0:ℝ)..1, t * Ps x * (P x - Ps x) :=
    intervalIntegral.integral_mul_deriv_eq_deriv_mul (fun x _ => hdd x)
      (fun x _ => (hPd x).hasDerivAt.fun_sub (hsd x))
      ((continuous_const.fun_mul hsc).intervalIntegrable 0 1)
      ((hPc.fun_sub hdc).intervalIntegrable 0 1)
  have hcross : ∫ x in (0:ℝ)..1, Pd x * (deriv P x - Pd x)
      = -(t * ∫ x in (0:ℝ)..1, Ps x * (P x - Ps x)) := by
    rw [hibp, hP.zero, hP.one, hs0, hs1]
    have : ∫ x in (0:ℝ)..1, t * Ps x * (P x - Ps x) = t * ∫ x in (0:ℝ)..1, Ps x * (P x - Ps x) := by
      rw [← intervalIntegral.integral_const_mul]
      congr 1
      funext x
      ring
    rw [this]
    ring
  have e1 : CP P = 1 * (∫ x in (0:ℝ)..1, Pd x ^ 2) + 2 * (∫ x in (0:ℝ)..1, Pd x * (deriv P x - Pd x))
      + 1 * ∫ x in (0:ℝ)..1, (deriv P x - Pd x) ^ 2 := by
    rw [← integral_lin3 1 2 1 (hdc.fun_pow 2) (hdc.fun_mul (hPc.fun_sub hdc))
      ((hPc.fun_sub hdc).fun_pow 2)]
    unfold CP
    apply intervalIntegral.integral_congr
    intro x _
    simp only
    ring
  have e2 : BP P = 1 * (∫ x in (0:ℝ)..1, Ps x ^ 2) + 2 * (∫ x in (0:ℝ)..1, Ps x * (P x - Ps x))
      + 1 * ∫ x in (0:ℝ)..1, (P x - Ps x) ^ 2 := by
    rw [← integral_lin3 1 2 1 (hsc.fun_pow 2) (hsc.fun_mul (hP0.fun_sub hsc))
      ((hP0.fun_sub hsc).fun_pow 2)]
    unfold BP
    apply intervalIntegral.integral_congr
    intro x _
    simp only
    ring
  have h3 : a * (∫ x in (0:ℝ)..1, Pd x * (deriv P x - Pd x))
      + b * (∫ x in (0:ℝ)..1, Ps x * (P x - Ps x)) = 0 := by
    rw [hcross, hab]
    ring
  have hZ : 0 ≤ ∫ x in (0:ℝ)..1, (deriv P x - Pd x) ^ 2 :=
    intervalIntegral.integral_nonneg zero_le_one (fun _ _ => sq_nonneg _)
  have hZ' : 0 ≤ ∫ x in (0:ℝ)..1, (P x - Ps x) ^ 2 :=
    intervalIntegral.integral_nonneg zero_le_one (fun _ _ => sq_nonneg _)
  have key : a * CP P + b * BP P - (a * (∫ x in (0:ℝ)..1, Pd x ^ 2) + b * (∫ x in (0:ℝ)..1, Ps x ^ 2))
      = a * (∫ x in (0:ℝ)..1, (deriv P x - Pd x) ^ 2) + b * ∫ x in (0:ℝ)..1, (P x - Ps x) ^ 2 := by
    rw [e1, e2]
    linear_combination (2:ℝ) * h3
  nlinarith [mul_nonneg ha hZ, mul_nonneg hb hZ']

/-- The minimiser P_λ = sinh(λx)/sinh λ, P_0 = x (DERIVATION §1.4). -/
def Pl (lam x : ℝ) : ℝ := if lam = 0 then x else Real.sinh (lam * x) / Real.sinh lam
/-- Its derivative. -/
def Pld (lam x : ℝ) : ℝ := if lam = 0 then 1 else lam * Real.cosh (lam * x) / Real.sinh lam
/-- C_λ = ∫₀¹ P_λ'². -/
def Cl (lam : ℝ) : ℝ := ∫ x in (0:ℝ)..1, Pld lam x ^ 2
/-- B_λ = ∫₀¹ P_λ². -/
def Bl (lam : ℝ) : ℝ := ∫ x in (0:ℝ)..1, Pl lam x ^ 2

lemma Pl_hasDerivAt (lam x : ℝ) : HasDerivAt (Pl lam) (Pld lam x) x := by
  by_cases h : lam = 0
  · have e : Pl lam = fun x => x := by funext y; simp [Pl, h]
    have e' : Pld lam x = 1 := by simp [Pld, h]
    rw [e, e']
    exact hasDerivAt_id' (x := x)
  · have e : Pl lam = fun y => Real.sinh (lam * y) / Real.sinh lam := by funext y; simp [Pl, h]
    have e' : Pld lam x = lam * Real.cosh (lam * x) / Real.sinh lam := by simp [Pld, h]
    rw [e, e']
    have := (((hasDerivAt_id' (x := x)).const_mul lam).sinh).div_const (Real.sinh lam)
    convert this using 1
    ring

lemma Pld_hasDerivAt (lam x : ℝ) : HasDerivAt (Pld lam) (lam ^ 2 * Pl lam x) x := by
  by_cases h : lam = 0
  · have e : Pld lam = fun _ => 1 := by funext y; simp [Pld, h]
    rw [e, h]
    simpa using hasDerivAt_const (x := x) (c := (1:ℝ))
  · have e : Pld lam = fun y => lam * Real.cosh (lam * y) / Real.sinh lam := by
      funext y; simp [Pld, h]
    have e' : Pl lam x = Real.sinh (lam * x) / Real.sinh lam := by simp [Pl, h]
    rw [e, e']
    have := ((((hasDerivAt_id' (x := x)).const_mul lam).cosh).const_mul lam).div_const (Real.sinh lam)
    convert this using 1
    ring

lemma Pl_zero (lam : ℝ) : Pl lam 0 = 0 := by simp [Pl]

lemma Pl_one (lam : ℝ) : Pl lam 1 = 1 := by
  by_cases h : lam = 0
  · simp [Pl, h]
  · rw [show Pl lam 1 = Real.sinh lam / Real.sinh lam by simp [Pl, h]]
    exact div_self (Real.sinh_ne_zero.2 h)

lemma Pl_cont (lam : ℝ) : Continuous (Pl lam) :=
  continuous_iff_continuousAt.2 fun x => (Pl_hasDerivAt lam x).continuousAt

lemma Pld_cont (lam : ℝ) : Continuous (Pld lam) :=
  continuous_iff_continuousAt.2 fun x => (Pld_hasDerivAt lam x).continuousAt

lemma deriv_Pl (lam : ℝ) : deriv (Pl lam) = Pld lam :=
  funext fun x => (Pl_hasDerivAt lam x).deriv

lemma Pl_adm (lam : ℝ) : PAdm (Pl lam) where
  smooth := contDiff_one_iff_deriv.mpr
    ⟨fun x => (Pl_hasDerivAt lam x).differentiableAt, by rw [deriv_Pl]; exact Pld_cont lam⟩
  zero := Pl_zero lam
  one := Pl_one lam

lemma CP_Pl (lam : ℝ) : CP (Pl lam) = Cl lam := by
  unfold CP Cl
  rw [deriv_Pl]

lemma BP_Pl (lam : ℝ) : BP (Pl lam) = Bl lam := rfl

/-- Inner minimum over P (DERIVATION §1.4). -/
theorem inner_min {P : ℝ → ℝ} (hP : PAdm P) {a b lam : ℝ} (ha : 0 ≤ a) (hb : 0 ≤ b)
    (hab : b = a * lam ^ 2) : a * Cl lam + b * Bl lam ≤ a * CP P + b * BP P :=
  quad_ge hP (Pl_zero lam) (Pl_one lam) (Pl_hasDerivAt lam) (Pld_hasDerivAt lam) ha hb hab

lemma Cl_zero : Cl 0 = 1 := by simp [Cl, Pld]

/-- C_λ + λ² B_λ = P_λ'(1) (integration by parts). -/
lemma Cl_add (lam : ℝ) : Cl lam + lam ^ 2 * Bl lam = Pld lam 1 := by
  have h : ∫ x in (0:ℝ)..1, Pld lam x * Pld lam x
      = Pld lam 1 * Pl lam 1 - Pld lam 0 * Pl lam 0 - ∫ x in (0:ℝ)..1, lam ^ 2 * Pl lam x * Pl lam x :=
    intervalIntegral.integral_mul_deriv_eq_deriv_mul (fun x _ => Pld_hasDerivAt lam x)
      (fun x _ => Pl_hasDerivAt lam x)
      ((continuous_const.fun_mul (Pl_cont lam)).intervalIntegrable 0 1)
      ((Pld_cont lam).intervalIntegrable 0 1)
  rw [Pl_one, Pl_zero] at h
  have e1 : Cl lam = ∫ x in (0:ℝ)..1, Pld lam x * Pld lam x := by
    unfold Cl; congr 1; funext x; ring
  have e2 : (∫ x in (0:ℝ)..1, lam ^ 2 * Pl lam x * Pl lam x) = lam ^ 2 * Bl lam := by
    unfold Bl; rw [← intervalIntegral.integral_const_mul]; congr 1; funext x; ring
  rw [e1, h, e2]
  ring

/-- §1.5: B_λ is non-increasing on λ ≥ 0. -/
theorem Bl_anti {l1 l2 : ℝ} (h1 : 0 ≤ l1) (h12 : l1 ≤ l2) : Bl l2 ≤ Bl l1 := by
  have m1 := inner_min (Pl_adm l2) (a := 1) (b := l1 ^ 2) (lam := l1) zero_le_one (sq_nonneg _)
    (by ring)
  have m2 := inner_min (Pl_adm l1) (a := 1) (b := l2 ^ 2) (lam := l2) zero_le_one (sq_nonneg _)
    (by ring)
  rw [CP_Pl, BP_Pl] at m1 m2
  by_contra hc0
  have hc := not_le.mp hc0
  have hl : l1 ^ 2 < l2 ^ 2 := by
    rcases eq_or_lt_of_le h12 with h | h
    · rw [h] at hc
      exact absurd hc (lt_irrefl _)
    · nlinarith
  nlinarith [mul_pos (sub_pos.2 hl) (sub_pos.2 hc)]

/-- §1.5: C_λ is non-decreasing on λ ≥ 0. -/
theorem Cl_mono {l1 l2 : ℝ} (h1 : 0 ≤ l1) (h12 : l1 ≤ l2) : Cl l1 ≤ Cl l2 := by
  have m1 := inner_min (Pl_adm l2) (a := 1) (b := l1 ^ 2) (lam := l1) zero_le_one (sq_nonneg _)
    (by ring)
  rw [CP_Pl, BP_Pl] at m1
  have hB := Bl_anti h1 h12
  nlinarith [mul_nonneg (sq_nonneg l1) (sub_nonneg.2 hB)]

lemma Cl_ge_one {lam : ℝ} (h : 0 ≤ lam) : 1 ≤ Cl lam := Cl_zero ▸ Cl_mono le_rfl h

lemma Bl_pos (lam : ℝ) : 0 < Bl lam := by
  unfold Bl
  apply intervalIntegral.integral_pos (by norm_num : (0:ℝ) < 1) ((Pl_cont lam).fun_pow 2).continuousOn
  · intro x _
    positivity
  · exact ⟨1, ⟨zero_le_one, le_refl 1⟩, by simp [Pl_one]⟩

lemma coth_anti {x y : ℝ} (hx : 0 < x) (hxy : x ≤ y) :
    Real.cosh y / Real.sinh y ≤ Real.cosh x / Real.sinh x := by
  have hsx := Real.sinh_pos_iff.2 hx
  have hsy := Real.sinh_pos_iff.2 (hx.trans_le hxy)
  rw [div_le_div_iff₀ hsy hsx]
  have h := Real.sinh_nonneg_iff.2 (sub_nonneg.2 hxy)
  rw [Real.sinh_sub] at h
  nlinarith [h]

lemma val_eq {W0 W1 ξ : ℝ} (h0 : 0 < W0) (h1 : 0 ≤ W1) (hξ : 0 < ξ) :
    Cl (ξ * Real.sqrt (W1 / W0)) / ξ * W0 + ξ * Bl (ξ * Real.sqrt (W1 / W0)) * W1
      = W0 / ξ * Pld (ξ * Real.sqrt (W1 / W0)) 1 := by
  rw [← Cl_add]
  have hk : Real.sqrt (W1 / W0) ^ 2 = W1 / W0 := Real.sq_sqrt (div_nonneg h1 h0.le)
  have : (ξ * Real.sqrt (W1 / W0)) ^ 2 = ξ ^ 2 * (W1 / W0) := by rw [mul_pow, hk]
  rw [this]
  have hξ' := hξ.ne'
  have h0' := h0.ne'
  field_simp

/-- §1.7: the optimal value is non-increasing in ξ. -/
theorem xi_mono {W0 W1 ξ ξ0 : ℝ} (h0 : 0 < W0) (h1 : 0 ≤ W1) (hξ : 0 < ξ) (hξξ : ξ ≤ ξ0) :
    Cl (ξ0 * Real.sqrt (W1 / W0)) / ξ0 * W0 + ξ0 * Bl (ξ0 * Real.sqrt (W1 / W0)) * W1
      ≤ Cl (ξ * Real.sqrt (W1 / W0)) / ξ * W0 + ξ * Bl (ξ * Real.sqrt (W1 / W0)) * W1 := by
  have hξ0 : 0 < ξ0 := hξ.trans_le hξξ
  rw [val_eq h0 h1 hξ, val_eq h0 h1 hξ0]
  have hk0 : 0 ≤ Real.sqrt (W1 / W0) := Real.sqrt_nonneg _
  rcases eq_or_lt_of_le hk0 with hk' | hk'
  · rw [← hk', mul_zero, mul_zero]
    have : Pld 0 1 = 1 := by simp [Pld]
    rw [this, mul_one, mul_one]
    gcongr
  · have hne : ∀ z, 0 < z → Pld z 1 = z * Real.cosh z / Real.sinh z := fun z hz => by
      simp [Pld, hz.ne']
    have ha := mul_pos hξ hk'
    have hb := mul_pos hξ0 hk'
    rw [hne _ hb, hne _ ha]
    have hc := coth_anti ha (mul_le_mul_of_nonneg_right hξξ hk0)
    have hsa := (Real.sinh_pos_iff.2 ha).ne'
    have hsb := (Real.sinh_pos_iff.2 hb).ne'
    have e1 : W0 / ξ0 * (ξ0 * Real.sqrt (W1 / W0) * Real.cosh (ξ0 * Real.sqrt (W1 / W0))
        / Real.sinh (ξ0 * Real.sqrt (W1 / W0)))
        = W0 * Real.sqrt (W1 / W0) * (Real.cosh (ξ0 * Real.sqrt (W1 / W0))
          / Real.sinh (ξ0 * Real.sqrt (W1 / W0))) := by
      field_simp
    have e2 : W0 / ξ * (ξ * Real.sqrt (W1 / W0) * Real.cosh (ξ * Real.sqrt (W1 / W0))
        / Real.sinh (ξ * Real.sqrt (W1 / W0)))
        = W0 * Real.sqrt (W1 / W0) * (Real.cosh (ξ * Real.sqrt (W1 / W0))
          / Real.sinh (ξ * Real.sqrt (W1 / W0))) := by
      field_simp
    rw [e1, e2]
    exact mul_le_mul_of_nonneg_left hc (mul_nonneg h0.le hk0)

/-- J(R, λ, Q) at ξ: Bernard's c after the sinh step (DERIVATION §1.4). -/
def J (ξ R lam : ℝ) (Q : ℝ → ℝ) : ℝ := Jag R (Cl lam / ξ) (ξ * Bl lam) Q

/-- §1.4 + §1.7: c(P,Q,R,ξ) ≥ J(ξ₀, R, λ₀, Q) with λ₀ = ξ₀ √(W₁/W₀), for 0 < ξ ≤ ξ₀. -/
theorem c_ge_J {P Q : ℝ → ℝ} (hP : PAdm P) (hQ : ContDiff ℝ 1 Q) (hQ0 : Q 0 = 1) (R : ℝ)
    {ξ ξ0 : ℝ} (hξ : 0 < ξ) (hξξ : ξ ≤ ξ0) :
    J ξ0 R (ξ0 * Real.sqrt (W1 R Q / W0 R Q)) Q ≤ c P Q R ξ := by
  obtain ⟨hQd, hQc⟩ := contDiff_one_iff_deriv.mp hQ
  have h0 := W0_pos (R := R) hQd.continuous hQ0
  have h1 := W1_nonneg R Q
  have hm := xi_mono h0 h1 hξ hξξ
  have hk2 : Real.sqrt (W1 R Q / W0 R Q) ^ 2 = W1 R Q / W0 R Q := Real.sq_sqrt (div_nonneg h1 h0.le)
  have hξ' := hξ.ne'
  have h0' := h0.ne'
  have hi := inner_min hP (a := W0 R Q / ξ) (b := ξ * W1 R Q)
    (lam := ξ * Real.sqrt (W1 R Q / W0 R Q)) (div_nonneg h0.le hξ.le) (mul_nonneg hξ.le h1)
    (by rw [mul_pow, hk2]; field_simp)
  rw [c_expand hP hQ hQ0 R hξ']
  unfold J Jag
  have e1 : Cl (ξ * Real.sqrt (W1 R Q / W0 R Q)) / ξ * W0 R Q
      + ξ * Bl (ξ * Real.sqrt (W1 R Q / W0 R Q)) * W1 R Q
      = W0 R Q / ξ * Cl (ξ * Real.sqrt (W1 R Q / W0 R Q))
        + ξ * W1 R Q * Bl (ξ * Real.sqrt (W1 R Q / W0 R Q)) := by ring
  have e2 : CP P / ξ * W0 R Q + ξ * BP P * W1 R Q = W0 R Q / ξ * CP P + ξ * W1 R Q * BP P := by
    ring
  linarith

end GL2Ceil


