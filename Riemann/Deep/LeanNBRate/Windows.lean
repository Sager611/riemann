import Riemann.Deep.LeanNBRate.Gevrey
import Mathlib.Analysis.Calculus.IteratedDeriv.Lemmas

/-!
# Gevrey-2 windows

`GevB C A f`: `f` smooth with `|f^{(n)}| ≤ C Aⁿ n^{2n}` for all `n` (`0⁰ = 1`). Stable under affine
changes of variable and products. Windows built from the smooth step `Gevrey.rhoTwo`:
* `cutoff a b η` = 1 on `[a − η/2, b + η/2]`, 0 off `(a − η, b + η)` (frequency cutoff);
* `twin ℓ` = 1 on `[ℓ/4, 3ℓ/4]`, 0 off `(0, ℓ)` (time window).
-/

open Filter Topology Set

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- Gevrey-2 bounds with constants `C`, `A`. -/
structure GevB (C A : ℝ) (f : ℝ → ℝ) : Prop where
  smooth : ContDiff ℝ (⊤ : ℕ∞) f
  bound : ∀ (n : ℕ) (x : ℝ), |iteratedDeriv n f x| ≤ C * A ^ n * (n : ℝ) ^ (2 * n)

/-- `B₀ = 2e⁸`, `A₀ = 36/e`: the constants of `rhoTwo`. -/
def B₀ : ℝ := 2 * Real.exp 8
def A₀ : ℝ := 36 / Real.exp 1

theorem B₀_pos : 0 < B₀ := by unfold B₀; positivity
theorem A₀_pos : 0 < A₀ := by unfold A₀; positivity
theorem one_le_B₀ : 1 ≤ B₀ := by
  unfold B₀
  have : 1 ≤ Real.exp 8 := Real.one_le_exp (by norm_num)
  linarith

theorem rhoTwo_gevB : GevB B₀ A₀ Gevrey.rhoTwo := by
  have hP := Gevrey.gevreyProfile_rhoTwo
  refine ⟨hP.smooth, fun n x => ?_⟩
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · simp only [iteratedDeriv_zero, pow_zero, mul_one, Nat.cast_zero, mul_zero]
    have h0 := hP.taper.monotone (show (0:ℝ) ≤ max 0 (min x 1) from le_max_left _ _)
    have h1 : 0 ≤ Gevrey.rhoTwo x ∧ Gevrey.rhoTwo x ≤ 1 := by
      constructor
      · rcases le_total x 0 with hx | hx
        · rw [hP.taper.eq_zero x hx]
        · rw [← hP.taper.eq_zero 0 le_rfl]; exact hP.taper.monotone hx
      · rcases le_total x 1 with hx | hx
        · rw [← hP.taper.eq_one 1 le_rfl]; exact hP.taper.monotone hx
        · rw [hP.taper.eq_one x hx]
    rw [abs_of_nonneg h1.1]
    exact h1.2.trans one_le_B₀
  · have := hP.bound n hn x
    have e : ((n : ℝ) ^ ((2 : ℝ) * n)) = (n : ℝ) ^ (2 * n) := by
      rw [show ((2 : ℝ) * n) = ((2 * n : ℕ) : ℝ) by push_cast; ring, Real.rpow_natCast]
    rw [e] at this
    simpa [B₀, A₀] using this

namespace GevB

theorem nonneg_C {C A : ℝ} {f : ℝ → ℝ} (h : GevB C A f) : 0 ≤ C := by
  have := h.bound 0 0
  simp only [pow_zero, mul_one, Nat.cast_zero, mul_zero] at this
  exact (abs_nonneg _).trans this

/-- Affine change of variable `x ↦ αx + β`. -/
theorem comp_affine {C A : ℝ} {f : ℝ → ℝ} (h : GevB C A f) (α β : ℝ) :
    GevB C (|α| * A) (fun x => f (α * x + β)) := by
  have hs : ContDiff ℝ (⊤ : ℕ∞) (fun y => f (y + β)) := h.smooth.comp (contDiff_id.add contDiff_const)
  refine ⟨hs.comp (contDiff_const.mul contDiff_id), fun n x => ?_⟩
  have e1 : iteratedDeriv n (fun x => f (α * x + β)) x =
      α ^ n * iteratedDeriv n (fun y => f (y + β)) (α * x) := by
    have := iteratedDeriv_comp_const_mul (n := n) (hs.of_le (WithTop.coe_le_coe.mpr le_top)) α
    exact congrFun this x
  rw [e1, iteratedDeriv_comp_add_const, abs_mul, abs_pow]
  calc |α| ^ n * |iteratedDeriv n f (α * x + β)| ≤ |α| ^ n * (C * A ^ n * (n : ℝ) ^ (2 * n)) := by
        gcongr; exact h.bound n _
    _ = C * (|α| * A) ^ n * (n : ℝ) ^ (2 * n) := by ring

/-- Products (Leibniz). -/
theorem mul {C₁ A₁ C₂ A₂ : ℝ} {f g : ℝ → ℝ} (hf : GevB C₁ A₁ f) (hg : GevB C₂ A₂ g)
    (hA₁ : 0 ≤ A₁) (hA₂ : 0 ≤ A₂) :
    GevB (C₁ * C₂) (2 * max A₁ A₂) (fun x => f x * g x) := by
  refine ⟨hf.smooth.mul hg.smooth, fun n x => ?_⟩
  have hC₁ := hf.nonneg_C
  have hC₂ := hg.nonneg_C
  have hL := iteratedDeriv_mul (n := n) (x := x)
    (hf.smooth.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt
    (hg.smooth.of_le (WithTop.coe_le_coe.mpr le_top)).contDiffAt
  have hfg : iteratedDeriv n (fun x => f x * g x) x = iteratedDeriv n (f * g) x := rfl
  rw [hfg, hL]
  set A := max A₁ A₂
  have hA : 0 ≤ A := le_max_of_le_left hA₁
  calc |∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * iteratedDeriv i f x * iteratedDeriv (n - i) g x|
      ≤ ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * (C₁ * A ^ i * (n : ℝ) ^ (2 * i)) *
          (C₂ * A ^ (n - i) * (n : ℝ) ^ (2 * (n - i))) := by
        refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun i hi => ?_)
        rw [Finset.mem_range] at hi
        rw [abs_mul, abs_mul, Nat.abs_cast]
        gcongr
        · refine (hf.bound i x).trans ?_
          gcongr
          · exact le_max_left _ _
          · exact_mod_cast (Nat.le_of_lt_succ hi)
        · refine (hg.bound (n - i) x).trans ?_
          gcongr
          · exact le_max_right _ _
          · exact_mod_cast (Nat.sub_le n i)
    _ = ∑ i ∈ Finset.range (n + 1), (n.choose i : ℝ) * (C₁ * C₂ * A ^ n * (n : ℝ) ^ (2 * n)) := by
        refine Finset.sum_congr rfl fun i hi => ?_
        rw [Finset.mem_range] at hi
        have h1 : A ^ i * A ^ (n - i) = A ^ n := by rw [← pow_add]; congr 1; omega
        have h2 : (n : ℝ) ^ (2 * i) * (n : ℝ) ^ (2 * (n - i)) = (n : ℝ) ^ (2 * n) := by
          rw [← pow_add]; congr 1; omega
        calc (n.choose i : ℝ) * (C₁ * A ^ i * (n : ℝ) ^ (2 * i)) *
              (C₂ * A ^ (n - i) * (n : ℝ) ^ (2 * (n - i)))
            = (n.choose i : ℝ) * (C₁ * C₂ * (A ^ i * A ^ (n - i)) *
                ((n : ℝ) ^ (2 * i) * (n : ℝ) ^ (2 * (n - i)))) := by ring
          _ = _ := by rw [h1, h2]
    _ = C₁ * C₂ * (2 * A) ^ n * (n : ℝ) ^ (2 * n) := by
        rw [← Finset.sum_mul, ← Nat.cast_sum, Nat.sum_range_choose]; push_cast; ring

end GevB

/-- The smooth step and its basic values. -/
theorem rhoTwo_eq_zero {x : ℝ} (hx : x ≤ 0) : Gevrey.rhoTwo x = 0 :=
  Gevrey.gevreyProfile_rhoTwo.taper.eq_zero x hx

theorem rhoTwo_eq_one {x : ℝ} (hx : 1 ≤ x) : Gevrey.rhoTwo x = 1 :=
  Gevrey.gevreyProfile_rhoTwo.taper.eq_one x hx

theorem rhoTwo_nonneg (x : ℝ) : 0 ≤ Gevrey.rhoTwo x := by
  rcases le_total x 0 with hx | hx
  · rw [rhoTwo_eq_zero hx]
  · rw [← rhoTwo_eq_zero (le_refl 0)]; exact Gevrey.gevreyProfile_rhoTwo.taper.monotone hx

theorem rhoTwo_le_one (x : ℝ) : Gevrey.rhoTwo x ≤ 1 := by
  rcases le_total x 1 with hx | hx
  · rw [← rhoTwo_eq_one (le_refl 1)]; exact Gevrey.gevreyProfile_rhoTwo.taper.monotone hx
  · rw [rhoTwo_eq_one hx]

/-- Frequency cutoff: `1` on `[a − η/2, b + η/2]`, `0` off `(a − η, b + η)`. -/
def cutoff (a b η : ℝ) (θ : ℝ) : ℝ :=
  Gevrey.rhoTwo ((2 / η) * θ + (-(2 / η) * (a - η))) *
    Gevrey.rhoTwo ((-(2 / η)) * θ + (2 / η) * (b + η))

theorem cutoff_gevB (a b η : ℝ) (hη : 0 < η) :
    GevB (B₀ * B₀) (2 * (2 / η * A₀)) (cutoff a b η) := by
  have h1 := rhoTwo_gevB.comp_affine (2 / η) (-(2 / η) * (a - η))
  have h2 := rhoTwo_gevB.comp_affine (-(2 / η)) ((2 / η) * (b + η))
  have hp : 0 ≤ 2 / η * A₀ := by have := A₀_pos; positivity
  have := h1.mul h2 (by rw [abs_of_pos (by positivity)]; exact hp)
    (by rw [abs_neg, abs_of_pos (by positivity)]; exact hp)
  rw [abs_neg, abs_of_pos (by positivity : (0:ℝ) < 2 / η), max_self] at this
  exact this

theorem one_le_two_div_mul {η x : ℝ} (hη : 0 < η) (hx : η / 2 ≤ x) : 1 ≤ 2 / η * x := by
  have : 2 / η * (η / 2) = 1 := by field_simp
  calc (1 : ℝ) = 2 / η * (η / 2) := this.symm
    _ ≤ 2 / η * x := by gcongr

theorem cutoff_eq_one {a b η θ : ℝ} (hη : 0 < η) (h1 : a - η / 2 ≤ θ) (h2 : θ ≤ b + η / 2) :
    cutoff a b η θ = 1 := by
  unfold cutoff
  rw [rhoTwo_eq_one, rhoTwo_eq_one, one_mul]
  · rw [show -(2 / η) * θ + 2 / η * (b + η) = 2 / η * (b + η - θ) by ring]
    exact one_le_two_div_mul hη (by linarith)
  · rw [show 2 / η * θ + -(2 / η) * (a - η) = 2 / η * (θ - (a - η)) by ring]
    exact one_le_two_div_mul hη (by linarith)

theorem cutoff_eq_zero {a b η θ : ℝ} (hη : 0 < η) (h : θ ≤ a - η ∨ b + η ≤ θ) :
    cutoff a b η θ = 0 := by
  unfold cutoff
  rcases h with h | h
  · rw [rhoTwo_eq_zero, zero_mul]
    rw [show 2 / η * θ + -(2 / η) * (a - η) = 2 / η * (θ - (a - η)) by ring]
    exact mul_nonpos_of_nonneg_of_nonpos (by positivity) (by linarith)
  · rw [rhoTwo_eq_zero (x := -(2 / η) * θ + 2 / η * (b + η)), mul_zero]
    rw [show -(2 / η) * θ + 2 / η * (b + η) = 2 / η * (b + η - θ) by ring]
    exact mul_nonpos_of_nonneg_of_nonpos (by positivity) (by linarith)

theorem cutoff_nonneg (a b η θ : ℝ) : 0 ≤ cutoff a b η θ :=
  mul_nonneg (rhoTwo_nonneg _) (rhoTwo_nonneg _)

theorem cutoff_le_one (a b η θ : ℝ) : cutoff a b η θ ≤ 1 :=
  by unfold cutoff; nlinarith [rhoTwo_le_one ((2 / η) * θ + (-(2 / η) * (a - η))),
    rhoTwo_nonneg ((2 / η) * θ + (-(2 / η) * (a - η))),
    rhoTwo_le_one ((-(2 / η)) * θ + (2 / η) * (b + η)), rhoTwo_nonneg ((-(2 / η)) * θ + (2 / η) * (b + η))]

theorem cutoff_support {a b η : ℝ} (hη : 0 < η) :
    Function.support (cutoff a b η) ⊆ Ioo (a - η) (b + η) := by
  intro θ hθ
  by_contra h
  rw [mem_Ioo, not_and_or, not_lt, not_lt] at h
  exact hθ (cutoff_eq_zero hη h)

/-- `cutoff = 1` on a neighbourhood of every point of `[a, b]`. -/
theorem cutoff_eventuallyEq_one {a b η θ : ℝ} (hη : 0 < η) (hθ : θ ∈ Icc a b) :
    cutoff a b η =ᶠ[𝓝 θ] fun _ => 1 := by
  have : Ioo (a - η / 2) (b + η / 2) ∈ 𝓝 θ :=
    Ioo_mem_nhds (by linarith [hθ.1]) (by linarith [hθ.2])
  filter_upwards [this] with x hx
  exact cutoff_eq_one hη hx.1.le hx.2.le

/-- Time window on `[0, ℓ]`: `1` on `[ℓ/4, 3ℓ/4]`, `0` off `(0, ℓ)`. -/
def twin (ℓ : ℝ) (u : ℝ) : ℝ := cutoff (ℓ / 2) (ℓ / 2) (ℓ / 2) u

theorem twin_gevB {ℓ : ℝ} (hℓ : 0 < ℓ) :
    GevB (B₀ * B₀) (2 * (2 / (ℓ / 2) * A₀)) (twin ℓ) :=
  cutoff_gevB _ _ _ (by positivity)

theorem twin_eq_one {ℓ u : ℝ} (hℓ : 0 < ℓ) (hu : u ∈ Icc (ℓ / 4) (3 * ℓ / 4)) : twin ℓ u = 1 :=
  cutoff_eq_one (by positivity) (by linarith [hu.1]) (by linarith [hu.2])

theorem twin_eq_zero {ℓ u : ℝ} (hℓ : 0 < ℓ) (hu : u ≤ 0 ∨ ℓ ≤ u) : twin ℓ u = 0 :=
  cutoff_eq_zero (by positivity) (by rcases hu with h | h <;> [left; right] <;> linarith)

theorem twin_nonneg (ℓ u : ℝ) : 0 ≤ twin ℓ u := cutoff_nonneg _ _ _ _

theorem twin_le_one (ℓ u : ℝ) : twin ℓ u ≤ 1 := cutoff_le_one _ _ _ _

end Riemann.Deep.LeanNBRate
