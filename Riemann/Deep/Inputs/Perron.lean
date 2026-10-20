import Riemann.Deep.Inputs.PerronKernel
import Riemann.Deep.Inputs.PerronShift
import Riemann.Deep.Inputs.All

/-!
# perron33 is a theorem; hence Fayed.Inputs holds unconditionally

Assembly of the truncated Perron formula (c = 2, T = u⁸) and the contour shift:

* On Re z = 2, ζ(s+z) u^z / z = Σ_n (n+1)^{-s} e^{L_n z}/z with L_n = log u - log(n+1);
  sum and integral commute (absolutely summable), and each term is a Perron kernel
  (`perron_kernel`): error ≤ 2(u/(n+1))² / (π |L_n| T).
* On D_γ, |u - n| ≥ γ for every integer n, so 1/|L_n| ≤ (u + n + 1)/γ: total O(u³/T) = O(u⁻⁵).
* Contour shift (`shift_identity`) to Re z = 1-σ-η picks up u^{1-s}/(1-s); the horizontal sides
  are O((1 + log u) u^{E-8}), E = max(1-σ+7η, 2) < 4 (`horiz_bound`).
-/

open Complex Set MeasureTheory intervalIntegral Filter Topology
open scoped Real Interval

namespace Fayed.InputsProof

namespace PerronAux

noncomputable section

/-- log(n+1). -/
def lg (n : ℕ) : ℝ := Real.log ((n : ℝ) + 1)

lemma succ_pos' (n : ℕ) : (0 : ℝ) < (n : ℝ) + 1 := by positivity

lemma exp_neg_mul_lg (r : ℝ) (n : ℕ) : Real.exp (-(r * lg n)) = ((n : ℝ) + 1) ^ (-r) := by
  rw [Real.rpow_def_of_pos (succ_pos' n), lg]; ring_nf

lemma summable_rpow_succ {r : ℝ} (hr : 1 < r) :
    Summable (fun n : ℕ => ((n : ℝ) + 1) ^ (-r)) := by
  have h := (summable_nat_add_iff 1).mpr (Real.summable_nat_rpow_inv.mpr hr)
  refine h.congr (fun n => ?_)
  push_cast
  rw [Real.rpow_neg (by positivity)]

lemma natsucc_cpow (n : ℕ) (w : ℂ) : ((n : ℂ) + 1) ^ w = cexp ((lg n : ℂ) * w) := by
  have h : ((n : ℂ) + 1) = (((n : ℝ) + 1 : ℝ) : ℂ) := by push_cast; ring
  rw [h, cpow_def_of_ne_zero (by exact_mod_cast (succ_pos' n).ne'),
    ← ofReal_log (succ_pos' n).le, lg]

lemma real_cpow {u : ℝ} (hu : 0 < u) (w : ℂ) : (u : ℂ) ^ w = cexp ((Real.log u : ℂ) * w) := by
  rw [cpow_def_of_ne_zero (by exact_mod_cast hu.ne'), ← ofReal_log hu.le]

/-- Dirichlet coefficient (n+1)^{-s}. -/
def coef (s : ℂ) (n : ℕ) : ℂ := cexp (-((lg n : ℂ) * s))

/-- Perron frequency L_n = log u - log(n+1). -/
def Lf (u : ℝ) (n : ℕ) : ℝ := Real.log u - lg n

lemma norm_coef (s : ℂ) (n : ℕ) : ‖coef s n‖ = Real.exp (-(s.re * lg n)) := by
  rw [coef, Complex.norm_exp]; congr 1; simp; ring

lemma Fz_eq_tsum (s : ℂ) {u : ℝ} (hu : 0 < u) (hs : 0 < s.re) (y : ℝ) :
    Fz s u ((2 : ℝ) + y * I) = ∑' n : ℕ, coef s n * kern (Lf u n) ((2 : ℝ) + y * I) := by
  have hre : 1 < (s + ((2 : ℝ) + y * I)).re := by simp; linarith
  rw [Fz, zeta_eq_tsum_one_div_nat_add_one_cpow hre, ← tsum_mul_right, ← tsum_div_const]
  congr 1; funext n
  rw [natsucc_cpow, real_cpow hu, coef, kern, one_div, ← Complex.exp_neg]
  rw [← mul_div_assoc, ← Complex.exp_add, ← Complex.exp_add, Lf]
  congr 2; push_cast; ring

lemma kern_cont (L : ℝ) : Continuous (fun y : ℝ => kern L ((2 : ℝ) + y * I)) := by
  unfold kern
  refine Continuous.div (by fun_prop) (by fun_prop) (fun y => ?_)
  intro h; have := congrArg Complex.re h; simp at this

lemma norm_kern_two (L y : ℝ) : ‖kern L ((2 : ℝ) + y * I)‖ ≤ Real.exp (2 * L) / 2 := by
  have h := norm_kern_le_re (L := L) (x := 2) (y := y) (by norm_num)
  rwa [abs_of_pos (by norm_num : (0 : ℝ) < 2), mul_comm L 2] at h

lemma exp_combo (σ : ℝ) {u : ℝ} (hu : 0 < u) (n : ℕ) :
    Real.exp (-(σ * lg n)) * Real.exp (2 * Lf u n) = u ^ 2 * Real.exp (-((σ + 2) * lg n)) := by
  have hu2 : u ^ 2 = Real.exp (2 * Real.log u) := by
    rw [show (2 : ℝ) * Real.log u = ((2 : ℕ) : ℝ) * Real.log u by norm_num, Real.exp_nat_mul,
      Real.exp_log hu]
  rw [hu2, ← Real.exp_add, ← Real.exp_add, Lf]; congr 1; ring

/-- Truncated Perron on Re z = 2: the line integral as a convergent sum of kernels. -/
lemma hasSum_P (s : ℂ) {u T : ℝ} (hu : 0 < u) (hs : 0 < s.re) (hT : 0 < T) :
    HasSum (fun n => coef s n *
        (((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, kern (Lf u n) ((2 : ℝ) + y * I)))
      (((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, Fz s u ((2 : ℝ) + y * I)) := by
  set c : ℂ := ((1 / (2 * π) : ℝ) : ℂ)
  set F : ℕ → ℝ → ℂ := fun n y => coef s n * kern (Lf u n) ((2 : ℝ) + y * I) with hF
  have hint : ∀ n, Integrable (F n) (volume.restrict (Ioc (-T) T)) := fun n =>
    ((continuous_const.mul (kern_cont _)).integrableOn_Icc).mono_set Ioc_subset_Icc_self
  have hbd : ∀ n, ∫ y in Ioc (-T) T, ‖F n y‖ ≤ T * u ^ 2 * ((n : ℝ) + 1) ^ (-(s.re + 2)) := by
    intro n
    set B := Real.exp (-(s.re * lg n)) * (Real.exp (2 * Lf u n) / 2)
    have h1 : ∀ y, ‖F n y‖ ≤ B := by
      intro y
      simp only [F]
      rw [norm_mul, norm_coef]
      exact mul_le_mul_of_nonneg_left (norm_kern_two _ _) (Real.exp_pos _).le
    have hBint : IntegrableOn (fun _ : ℝ => B) (Ioc (-T) T) :=
      (continuous_const.integrableOn_Icc).mono_set Ioc_subset_Icc_self
    calc ∫ y in Ioc (-T) T, ‖F n y‖ ≤ ∫ y in Ioc (-T) T, B :=
          setIntegral_mono_on (hint n).norm hBint measurableSet_Ioc (fun y _ => h1 y)
      _ = (T - -T) * B := by
          rw [setIntegral_const, Real.volume_real_Ioc_of_le (by linarith), smul_eq_mul]
      _ = T * u ^ 2 * ((n : ℝ) + 1) ^ (-(s.re + 2)) := by
          rw [← exp_neg_mul_lg]
          have := exp_combo s.re hu n
          simp only [B]
          rw [show (T - -T) * (Real.exp (-(s.re * lg n)) * (Real.exp (2 * Lf u n) / 2)) =
            T * (Real.exp (-(s.re * lg n)) * Real.exp (2 * Lf u n)) by ring, this]
          ring
  have hsum : Summable (fun n => ∫ y in Ioc (-T) T, ‖F n y‖) :=
    Summable.of_nonneg_of_le (fun n => integral_nonneg (fun y => norm_nonneg _)) hbd
      ((summable_rpow_succ (by linarith)).mul_left (T * u ^ 2))
  have hswap := integral_tsum_of_summable_integral_norm hint hsum
  have hsumF : Summable (fun n => ∫ y in Ioc (-T) T, F n y) :=
    Summable.of_norm_bounded hsum (fun n => norm_integral_le_integral_norm _)
  have hf : ∀ n, coef s n * (c * ∫ y in (-T)..T, kern (Lf u n) ((2 : ℝ) + y * I)) =
      c * ∫ y in Ioc (-T) T, F n y := by
    intro n
    rw [intervalIntegral.integral_of_le (by linarith)]
    simp only [F]
    rw [MeasureTheory.integral_const_mul]
    ring
  have hFz : ∫ y in (-T)..T, Fz s u ((2 : ℝ) + y * I) = ∫ y in Ioc (-T) T, ∑' n, F n y := by
    rw [intervalIntegral.integral_of_le (by linarith)]
    congr 1; funext y; exact Fz_eq_tsum s hu hs y
  rw [hFz, ← hswap]
  simp_rw [hf]
  exact hsumF.hasSum.mul_left c

/-- D_γ stays γ away from the integers. -/
lemma Dγ_far {γ : ℝ} (hγ : 0 < γ) {M₀ : ℕ} {u : ℝ} (hu : u ∈ Dγ γ M₀) (n : ℕ) : γ ≤ |u - n| := by
  simp only [Dγ, unionIoo, mem_iUnion, mem_Ioo, aγ, bγ] at hu
  obtain ⟨M, h1, h2⟩ := hu
  rcases le_or_gt n (M₀ + M) with h | h
  · have : (n : ℝ) ≤ ((M₀ + M : ℕ) : ℝ) := by exact_mod_cast h
    rw [abs_of_pos (by linarith)]; linarith
  · have : ((M₀ + M : ℕ) : ℝ) + 1 ≤ n := by exact_mod_cast h
    rw [abs_of_neg (by linarith)]; linarith

lemma log_sub_ge {x y γ : ℝ} (hx : 0 < x) (hy : 0 < y) (hγ : 0 < γ) (h : γ ≤ |x - y|) :
    γ / (x + y) ≤ |Real.log x - Real.log y| := by
  rcases le_total y x with hxy | hxy
  · have hd : γ ≤ x - y := by rwa [abs_of_nonneg (by linarith)] at h
    have h1 := Real.log_le_sub_one_of_pos (div_pos hy hx)
    rw [Real.log_div hy.ne' hx.ne'] at h1
    have h2 : x * (Real.log y - Real.log x) ≤ y - x := by
      have e : x * (y / x - 1) = y - x := by field_simp
      calc x * (Real.log y - Real.log x) ≤ x * (y / x - 1) := mul_le_mul_of_nonneg_left h1 hx.le
        _ = y - x := e
    have h4 : 0 ≤ Real.log x - Real.log y := by linarith [Real.log_le_log hy hxy]
    rw [abs_of_nonneg h4, div_le_iff₀ (by linarith)]
    nlinarith
  · have hd : γ ≤ y - x := by rwa [abs_of_nonpos (by linarith), neg_sub] at h
    have h1 := Real.log_le_sub_one_of_pos (div_pos hx hy)
    rw [Real.log_div hx.ne' hy.ne'] at h1
    have h2 : y * (Real.log x - Real.log y) ≤ x - y := by
      have e : y * (x / y - 1) = x - y := by field_simp
      calc y * (Real.log x - Real.log y) ≤ y * (x / y - 1) := mul_le_mul_of_nonneg_left h1 hy.le
        _ = x - y := e
    have h4 : Real.log x - Real.log y ≤ 0 := by linarith [Real.log_le_log hx hxy]
    rw [abs_of_nonpos h4, div_le_iff₀ (by linarith)]
    nlinarith

lemma S_eq (s : ℂ) (u : ℝ) : S u s = ∑ n ∈ Finset.range ⌊u⌋₊, coef s n := by
  rw [S]
  symm
  refine Finset.sum_nbij' (fun n => n + 1) (fun k => k - 1) ?_ ?_ ?_ ?_ ?_
  · intro a ha; simp only [Finset.mem_range, Finset.mem_Icc] at ha ⊢; omega
  · intro a ha; simp only [Finset.mem_range, Finset.mem_Icc] at ha ⊢; omega
  · intro a _; simp
  · intro a ha; simp only [Finset.mem_Icc] at ha; omega
  · intro a _
    rw [coef]
    have : ((a + 1 : ℕ) : ℂ) = (a : ℂ) + 1 := by push_cast; ring
    rw [this, natsucc_cpow]; congr 1; ring

lemma hasSum_S (s : ℂ) {γ u : ℝ} {M₀ : ℕ} (hγ : 0 < γ) (hu : u ∈ Dγ γ M₀) (hu1 : 1 ≤ u) :
    HasSum (fun n => coef s n * (if 0 < Lf u n then 1 else 0)) (S u s) := by
  rw [S_eq]
  have hpos : ∀ n ∈ Finset.range ⌊u⌋₊, 0 < Lf u n := by
    intro n hn
    have hn' : n + 1 ≤ ⌊u⌋₊ := Finset.mem_range.mp hn
    have h1 : ((n + 1 : ℕ) : ℝ) ≤ u :=
      (Nat.cast_le.mpr hn').trans (Nat.floor_le (by linarith))
    have h2 := Dγ_far hγ hu (n + 1)
    push_cast at h1 h2
    have h3 : (n : ℝ) + 1 < u := by
      rcases h1.lt_or_eq with h | h
      · exact h
      · rw [h, sub_self, abs_zero] at h2; linarith
    simp only [Lf, lg]; linarith [Real.log_lt_log (succ_pos' n) h3]
  have hneg : ∀ n ∉ Finset.range ⌊u⌋₊, Lf u n < 0 := by
    intro n hn
    have hn' : ⌊u⌋₊ < n + 1 := by simp only [Finset.mem_range] at hn; omega
    have h1 : u < ((n + 1 : ℕ) : ℝ) := (Nat.floor_lt (by linarith)).mp hn'
    push_cast at h1
    simp only [Lf, lg]; linarith [Real.log_lt_log (by linarith) h1]
  have h : HasSum (fun n => coef s n * (if 0 < Lf u n then 1 else 0))
      (∑ n ∈ Finset.range ⌊u⌋₊, coef s n * (if 0 < Lf u n then 1 else 0)) :=
    hasSum_sum_of_ne_finset_zero (fun n hn => by simp [not_lt.mpr (hneg n hn).le])
  convert h using 1
  refine Finset.sum_congr rfl (fun n hn => ?_)
  simp [hpos n hn]

lemma term_bound (s : ℂ) {γ u T : ℝ} {M₀ : ℕ} (hγ : 0 < γ) (hu : u ∈ Dγ γ M₀) (hu1 : 1 ≤ u)
    (hT : 0 < T) (n : ℕ) :
    ‖coef s n * (if 0 < Lf u n then 1 else 0) - coef s n *
        (((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, kern (Lf u n) ((2 : ℝ) + y * I))‖ ≤
      4 * u ^ 3 / (π * γ * T) * ((n : ℝ) + 1) ^ (-(s.re + 1)) := by
  have hu0 : 0 < u := by linarith
  have hπ : 0 < π := Real.pi_pos
  set L := Lf u n with hLdef
  have hfar := Dγ_far hγ hu (n + 1)
  push_cast at hfar
  have hLabs : γ / (u + ((n : ℝ) + 1)) ≤ |L| := log_sub_ge hu0 (succ_pos' n) hγ hfar
  have hLpos : 0 < |L| := lt_of_lt_of_le (by positivity) hLabs
  have hL : L ≠ 0 := abs_pos.mp hLpos
  have hk := perron_kernel (L := L) (c := 2) (T := T) hL two_pos hT
  have e1 := exp_combo s.re hu0 n
  set X := Real.exp (-((s.re + 2) * lg n))
  set W := Real.exp (-((s.re + 1) * lg n))
  have e2 : X * ((n : ℝ) + 1) = W := by
    have : (n : ℝ) + 1 = Real.exp (lg n) := (Real.exp_log (succ_pos' n)).symm
    rw [this, ← Real.exp_add]; congr 1; ring
  have hX : 0 < X := Real.exp_pos _
  rw [← mul_sub, norm_mul, norm_coef, norm_sub_rev, ← exp_neg_mul_lg]
  have hγL : γ ≤ |L| * (u + ((n : ℝ) + 1)) := by rwa [div_le_iff₀ (by positivity)] at hLabs
  have hun : u + ((n : ℝ) + 1) ≤ 2 * u * ((n : ℝ) + 1) := by
    have : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    nlinarith
  calc Real.exp (-(s.re * lg n)) * ‖((1 / (2 * π) : ℝ) : ℂ) *
          (∫ y in (-T)..T, kern L ((2 : ℝ) + y * I)) - (if 0 < L then 1 else 0)‖
      ≤ Real.exp (-(s.re * lg n)) * (2 * Real.exp (2 * L) / (π * |L| * T)) :=
        mul_le_mul_of_nonneg_left hk (Real.exp_pos _).le
    _ = 2 * u ^ 2 * X / (π * |L| * T) := by
        rw [show 2 * u ^ 2 * X = 2 * (u ^ 2 * X) by ring, ← e1, ← hLdef]; ring
    _ ≤ 4 * u ^ 3 / (π * γ * T) * W := by
        rw [← e2, div_le_iff₀ (by positivity)]
        have e3 : 4 * u ^ 3 / (π * γ * T) * (X * ((n : ℝ) + 1)) * (π * |L| * T) =
            4 * u ^ 3 * X * ((n : ℝ) + 1) * |L| / γ := by
          field_simp
        rw [e3, le_div_iff₀ hγ]
        calc 2 * u ^ 2 * X * γ ≤ 2 * u ^ 2 * X * (|L| * (u + ((n : ℝ) + 1))) := by gcongr
          _ ≤ 2 * u ^ 2 * X * (|L| * (2 * u * ((n : ℝ) + 1))) := by gcongr
          _ = 4 * u ^ 3 * X * ((n : ℝ) + 1) * |L| := by ring

end

end PerronAux

open PerronAux

/-- **perron33** (truncated Perron, c = 2, T = u⁸, plus contour shift): a theorem. -/
theorem perron33_holds : Perron33 := by
  intro s hσ0 hσ1 η hη0 hη1 hη2 γ hγ0 _hγ1 ε hε
  set Z := ∑' n : ℕ, ((n : ℝ) + 1) ^ (-(s.re + 1)) with hZ
  set E := max (1 - s.re + 7 * η) 2 with hE
  have hZ0 : 0 ≤ Z := tsum_nonneg (fun n => by positivity)
  have hE4 : E ≤ 4 := max_le (seven_eta hη1 hη2).le (by norm_num)
  refine ⟨4 * Z / (π * γ) + 52 / π * (2 + 8 / ε), |s.im| + 2, ?_⟩
  intro M₀ u hu hU
  have hπ : 0 < π := Real.pi_pos
  have hu1 : 1 ≤ u := by linarith [abs_nonneg s.im]
  have hu0 : 0 < u := by linarith
  have huT : u ≤ u ^ 8 := le_self_pow₀ hu1 (by norm_num)
  have hT1 : |s.im| + 1 ≤ u ^ 8 := by linarith
  have hT0 : 0 < u ^ 8 := by positivity
  have htT : |s.im| < u ^ 8 := by linarith
  set c : ℂ := ((1 / (2 * π) : ℝ) : ℂ) with hc
  have hJ : J η u s = c * ∫ y in (-(u ^ 8))..(u ^ 8),
      Fz s u (((1 - s.re - η : ℝ) : ℂ) + y * I) := by
    rw [J]
    congr 1
    apply intervalIntegral.integral_congr
    intro y _
    have hw : s + (((1 - s.re - η : ℝ) : ℂ) + y * I) = ((1 - η : ℝ) : ℂ) + ((s.im + y : ℝ) : ℂ) * I := by
      apply Complex.ext <;> simp <;> ring
    simp only [Φ, Fz, hw]
    push_cast
    ring
  have hshift := shift_identity s hσ0 hη0 hη2 hu0 htT
  set P := c * ∫ y in (-(u ^ 8))..(u ^ 8), Fz s u ((2 : ℝ) + y * I)
  have hSP : ‖S u s - P‖ ≤ 4 * u ^ 3 / (π * γ * u ^ 8) * Z := by
    have hd := (hasSum_S s hγ0 hu hu1).sub (hasSum_P s hu0 hσ0 hT0)
    have hg : HasSum (fun n : ℕ => 4 * u ^ 3 / (π * γ * u ^ 8) * ((n : ℝ) + 1) ^ (-(s.re + 1)))
        (4 * u ^ 3 / (π * γ * u ^ 8) * Z) :=
      (summable_rpow_succ (by linarith)).hasSum.mul_left _
    exact hd.norm_le_of_bounded hg (fun n => term_bound s hγ0 hu hu1 hT0 n)
  have hH1 := horiz_bound s hη0 hη1 hη2 hu1 hT1 (y := -(u ^ 8)) (by rw [abs_neg, abs_of_pos hT0])
  have hH2 := horiz_bound s hη0 hη1 hη2 hu1 hT1 (y := u ^ 8) (by rw [abs_of_pos hT0])
  set Hb := ∫ x in (1 - s.re - η)..2, Fz s u (x + ((-(u ^ 8) : ℝ) : ℂ) * I)
  set Ht := ∫ x in (1 - s.re - η)..2, Fz s u (x + ((u ^ 8 : ℝ) : ℂ) * I)
  set K := 52 * (1 + Real.log (2 * u ^ 8)) * u ^ E / u ^ 8
  have hmain : ‖S u s - (J η u s + (u : ℂ) ^ (1 - s) / (1 - s))‖ ≤
      4 * u ^ 3 / (π * γ * u ^ 8) * Z + 1 / (2 * π) * (K + K) := by
    have e : S u s - (J η u s + (u : ℂ) ^ (1 - s) / (1 - s)) =
        (S u s - P) + (P - c * (∫ y in (-(u ^ 8))..(u ^ 8),
          Fz s u (((1 - s.re - η : ℝ) : ℂ) + y * I)) - (u : ℂ) ^ (1 - s) / (1 - s)) := by
      rw [hJ]; ring
    rw [e, hshift]
    refine (norm_add_le _ _).trans (add_le_add hSP ?_)
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
    refine mul_le_mul_of_nonneg_left ?_ (by positivity)
    refine (norm_sub_le _ _).trans ?_
    simp only [norm_mul, Complex.norm_I, one_mul]
    exact add_le_add hH1 hH2
  refine hmain.trans ?_
  -- numerics
  set R := u ^ (-4 + ε) with hR
  have h1 : 4 * u ^ 3 / (π * γ * u ^ 8) * Z ≤ 4 * Z / (π * γ) * R := by
    have e : 4 * u ^ 3 / (π * γ * u ^ 8) * Z = 4 * Z / (π * γ) * (u ^ 3 / u ^ 8) := by
      field_simp
    rw [e]
    refine mul_le_mul_of_nonneg_left ?_ (by positivity)
    have e5 : u ^ 3 / u ^ 8 = u ^ (-5 : ℝ) := by
      rw [Real.rpow_neg hu0.le, show (5 : ℝ) = ((5 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
      field_simp
    rw [e5, hR]
    exact Real.rpow_le_rpow_of_exponent_le hu1 (by linarith)
  have hlog : 1 + Real.log (2 * u ^ 8) ≤ (2 + 8 / ε) * u ^ ε := by
    rw [Real.log_mul (by norm_num) (by positivity), Real.log_pow]
    have hl2 : Real.log 2 ≤ 1 := by linarith [Real.log_two_lt_d9]
    have hlu : Real.log u ≤ u ^ ε / ε := Real.log_le_rpow_div hu0.le hε
    have hue : 1 ≤ u ^ ε := Real.one_le_rpow hu1 hε.le
    have : 8 / ε * u ^ ε = 8 * (u ^ ε / ε) := by ring
    push_cast
    nlinarith
  have h2 : 1 / (2 * π) * (K + K) ≤ 52 / π * (2 + 8 / ε) * R := by
    have eK : 1 / (2 * π) * (K + K) =
        52 / π * ((1 + Real.log (2 * u ^ 8)) * (u ^ E / u ^ 8)) := by
      simp only [K]; field_simp; ring
    rw [eK, mul_assoc]
    refine mul_le_mul_of_nonneg_left ?_ (by positivity)
    have hEu : u ^ ε * (u ^ E / u ^ 8) ≤ R := by
      have e8 : (u ^ 8 : ℝ) = u ^ (8 : ℝ) := by
        rw [show (8 : ℝ) = ((8 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
      rw [e8, ← mul_div_assoc, ← Real.rpow_add hu0, ← Real.rpow_sub hu0, hR]
      exact Real.rpow_le_rpow_of_exponent_le hu1 (by linarith)
    have hq : 0 ≤ u ^ E / u ^ 8 := by positivity
    calc (1 + Real.log (2 * u ^ 8)) * (u ^ E / u ^ 8) ≤ (2 + 8 / ε) * u ^ ε * (u ^ E / u ^ 8) :=
          mul_le_mul_of_nonneg_right hlog hq
      _ = (2 + 8 / ε) * (u ^ ε * (u ^ E / u ^ 8)) := by ring
      _ ≤ (2 + 8 / ε) * R := mul_le_mul_of_nonneg_left hEu (by positivity)
  rw [add_mul]
  linarith

/-- **Fayed.Inputs holds**: all six literature inputs are theorems. -/
theorem inputs_holds : Inputs := inputs_of_perron33 perron33_holds

/-- Fayed's reduction, now unconditional on the literature inputs: Theorem 1C ⇒ RH.
(Theorem 1C itself is false; see Riemann/Fayed/Theorem1.lean.) -/
theorem RH_of_theorem1C : Theorem1C → RiemannHypothesis := fayed_RH inputs_holds

end Fayed.InputsProof

#print axioms Fayed.InputsProof.perron33_holds
#print axioms Fayed.InputsProof.inputs_holds
#print axioms Fayed.InputsProof.RH_of_theorem1C
