module

public import Mathlib.Analysis.Calculus.MeanValue
public import Mathlib.Analysis.Calculus.IteratedDeriv.Defs
public import Mathlib.Analysis.Calculus.ContDiff.Deriv
public import Mathlib.Analysis.Calculus.Deriv.Support
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
public import Mathlib.Topology.Algebra.InfiniteSum.Real

/-!
# LeanLi: finite differences and decay of smooth exponential sums

`|Σ_{n∈ℤ} φ(n/N) e^{inx}| ≤ C N (1 + N|x|)^{-K}` for `|x| ≤ π`, `φ` smooth with compact support,
by `K`-fold summation by parts.
-/

@[expose] public section

open Complex Filter Topology Set
open scoped ContDiff

namespace LeanLi

/-- Iterated backward difference with step `h`: `D⁰g = g`, `D^{K+1}g(x) = D^K g(x) - D^K g(x - h)`. -/
noncomputable def fdiff (h : ℝ) (g : ℝ → ℝ) : ℕ → ℝ → ℝ
  | 0 => g
  | K + 1 => fun x => fdiff h g K x - fdiff h g K (x - h)

lemma fdiff_zero (h : ℝ) (g : ℝ → ℝ) : fdiff h g 0 = g := rfl

lemma fdiff_succ (h : ℝ) (g : ℝ → ℝ) (K : ℕ) (x : ℝ) :
    fdiff h g (K + 1) x = fdiff h g K x - fdiff h g K (x - h) := rfl

lemma hasDerivAt_fdiff (h : ℝ) {g g' : ℝ → ℝ} (hg : ∀ x, HasDerivAt g (g' x) x) (K : ℕ) :
    ∀ x, HasDerivAt (fdiff h g K) (fdiff h g' K x) x := by
  induction K with
  | zero => exact hg
  | succ K ih =>
    intro x
    have h1 := ih x
    have h2 := (ih (x - h)).comp_sub_const x h
    exact h1.sub h2

lemma abs_fdiff_le {h : ℝ} (hh : 0 ≤ h) (K : ℕ) :
    ∀ (g : ℝ → ℝ), ContDiff ℝ K g → ∀ M : ℝ, (∀ x, |iteratedDeriv K g x| ≤ M) →
      ∀ x, |fdiff h g K x| ≤ M * h ^ K := by
  induction K with
  | zero =>
    intro g _ M hM x
    simpa [fdiff_zero] using hM x
  | succ K ih =>
    intro g hg M hM x
    have hg2 : ContDiff ℝ ((K : WithTop ℕ∞) + 1) g := by exact_mod_cast hg
    have hg' := (contDiff_succ_iff_deriv.1 hg2)
    have hdiff : Differentiable ℝ g := hg'.1
    have hcd : ContDiff ℝ K (deriv g) := by exact_mod_cast hg'.2.2
    have hM' : ∀ y, |iteratedDeriv K (deriv g) y| ≤ M := by
      intro y; rw [← iteratedDeriv_succ']; exact hM y
    have hIH := ih (deriv g) hcd M hM'
    have hd : ∀ y, HasDerivAt (fdiff h g K) (fdiff h (deriv g) K y) y :=
      hasDerivAt_fdiff h (fun y => (hdiff y).hasDerivAt) K
    have hmvt := norm_image_sub_le_of_norm_deriv_le_segment' (f := fdiff h g K)
      (f' := fdiff h (deriv g) K) (C := M * h ^ K) (a := x - h) (b := x)
      (fun y _ => (hd y).hasDerivWithinAt) (fun y _ => by rw [Real.norm_eq_abs]; exact hIH y)
      x ⟨by linarith, le_rfl⟩
    rw [fdiff_succ, ← Real.norm_eq_abs]
    calc ‖fdiff h g K x - fdiff h g K (x - h)‖ ≤ M * h ^ K * (x - (x - h)) := hmvt
      _ = M * h ^ (K + 1) := by ring

/-- `D^K g` vanishes outside `[-R, R + K h]` if `g` vanishes outside `[-R, R]`. -/
lemma fdiff_eq_zero {h : ℝ} (hh : 0 ≤ h) {g : ℝ → ℝ} {R : ℝ} (hg : ∀ y, R < |y| → g y = 0) (K : ℕ) :
    ∀ y, (y < -R ∨ R + K * h < y) → fdiff h g K y = 0 := by
  induction K with
  | zero =>
    intro y hy
    simp only [Nat.cast_zero, zero_mul, add_zero] at hy
    refine hg y ?_
    rcases hy with hy | hy
    · exact lt_of_lt_of_le (by linarith) (neg_le_abs y)
    · exact lt_of_lt_of_le hy (le_abs_self y)
  | succ K ih =>
    intro y hy
    rw [fdiff_succ]
    have hK : (((K + 1 : ℕ) : ℝ)) * h = K * h + h := by push_cast; ring
    rcases hy with hy | hy
    · rw [ih y (Or.inl hy), ih (y - h) (Or.inl (by linarith)), sub_zero]
    · rw [hK] at hy
      rw [ih y (Or.inr (by nlinarith)), ih (y - h) (Or.inr (by linarith)), sub_zero]

/-- `Σ_{n ∈ ℤ} a(n) e^{inx}`. -/
noncomputable def esum (a : ℤ → ℝ) (x : ℝ) : ℂ := ∑' n : ℤ, (a n : ℂ) * Complex.exp (I * ((n : ℝ) * x))

/-- Backward difference of a sequence on `ℤ`. -/
def bdiff (a : ℤ → ℝ) : ℤ → ℝ := fun n => a n - a (n - 1)

/-- Finitely supported sequences. -/
def FinSupp (a : ℤ → ℝ) : Prop := ∃ S : Finset ℤ, ∀ n ∉ S, a n = 0

lemma finSupp_bdiff {a : ℤ → ℝ} (ha : FinSupp a) : FinSupp (bdiff a) := by
  obtain ⟨S, hS⟩ := ha
  refine ⟨S ∪ S.image (· + 1), fun n hn => ?_⟩
  simp only [Finset.mem_union, Finset.mem_image, not_or, not_exists, not_and] at hn
  have h1 : a n = 0 := hS n hn.1
  have h2 : a (n - 1) = 0 := hS (n - 1) (fun h => hn.2 (n - 1) h (by ring))
  simp [LeanLi.bdiff, h1, h2]

lemma finSupp_iter {a : ℤ → ℝ} (ha : FinSupp a) (K : ℕ) : FinSupp (bdiff^[K] a) := by
  induction K with
  | zero => simpa using ha
  | succ K ih => rw [Function.iterate_succ_apply']; exact finSupp_bdiff ih

lemma summable_of_finSupp {a : ℤ → ℝ} (ha : FinSupp a) (g : ℤ → ℂ) :
    Summable (fun n => (a n : ℂ) * g n) := by
  obtain ⟨S, hS⟩ := ha
  exact summable_of_ne_finset_zero fun n hn => by rw [hS n hn]; simp

/-- Summation by parts: `Σ (a(n) - a(n-1)) e^{inx} = (1 - e^{ix}) Σ a(n) e^{inx}`. -/
lemma esum_bdiff {a : ℤ → ℝ} (ha : FinSupp a) (x : ℝ) :
    esum (bdiff a) x = (1 - Complex.exp (I * x)) * esum a x := by
  unfold esum bdiff
  have hshift : ∑' n : ℤ, ((a (n - 1) : ℝ) : ℂ) * Complex.exp (I * ((n : ℝ) * x)) =
      Complex.exp (I * x) * ∑' n : ℤ, (a n : ℂ) * Complex.exp (I * ((n : ℝ) * x)) := by
    rw [← tsum_mul_left, ← (Equiv.addRight (1 : ℤ)).tsum_eq]
    refine tsum_congr fun m => ?_
    simp only [Equiv.coe_addRight, add_sub_cancel_right]
    have e : Complex.exp (I * ((((m + 1 : ℤ) : ℝ) : ℂ) * x)) =
        Complex.exp (I * x) * Complex.exp (I * (((m : ℝ) : ℂ) * x)) := by
      rw [← Complex.exp_add]; congr 1; push_cast; ring
    rw [e]; ring
  have hs1 := summable_of_finSupp ha (fun n => Complex.exp (I * ((n : ℝ) * x)))
  have hs2 : Summable (fun n : ℤ => ((a (n - 1) : ℝ) : ℂ) * Complex.exp (I * ((n : ℝ) * x))) := by
    obtain ⟨S, hS⟩ := ha
    exact summable_of_ne_finset_zero (s := S.image (· + 1)) fun n hn => by
      rw [hS (n - 1) (fun h => hn (Finset.mem_image.2 ⟨n - 1, h, by ring⟩))]; simp
  simp only [Complex.ofReal_sub, sub_mul]
  rw [hs1.tsum_sub hs2, hshift]
  ring

lemma esum_bdiff_iter {a : ℤ → ℝ} (ha : FinSupp a) (x : ℝ) (K : ℕ) :
    esum (bdiff^[K] a) x = (1 - Complex.exp (I * x)) ^ K * esum a x := by
  induction K with
  | zero => simp
  | succ K ih =>
    rw [Function.iterate_succ_apply', esum_bdiff (finSupp_iter ha K), ih, pow_succ]
    ring

lemma norm_esum_le {a : ℤ → ℝ} (S : Finset ℤ) (hS : ∀ n ∉ S, a n = 0) (B : ℝ)
    (hB : ∀ n, |a n| ≤ B) (x : ℝ) : ‖esum a x‖ ≤ S.card * B := by
  unfold esum
  rw [tsum_eq_sum (s := S) (fun n hn => by rw [hS n hn]; simp)]
  calc ‖∑ n ∈ S, (a n : ℂ) * Complex.exp (I * ((n : ℝ) * x))‖
      ≤ ∑ n ∈ S, ‖(a n : ℂ) * Complex.exp (I * ((n : ℝ) * x))‖ := norm_sum_le _ _
    _ ≤ ∑ n ∈ S, B := by
        refine Finset.sum_le_sum fun n _ => ?_
        have h1 : ‖Complex.exp (I * (((n : ℝ) : ℂ) * (x : ℂ)))‖ = 1 := by
          rw [Complex.norm_exp]; simp
        rw [norm_mul, h1, mul_one, Complex.norm_real, Real.norm_eq_abs]
        exact hB n
    _ = S.card * B := by simp

/-- Differences of `φ(n/N)` are scaled differences of `φ`. -/
lemma bdiff_iter_eq (φ : ℝ → ℝ) (N : ℝ) (K : ℕ) (n : ℤ) :
    (bdiff^[K] (fun m : ℤ => φ (m / N))) n = fdiff (1 / N) φ K (n / N) := by
  induction K generalizing n with
  | zero => rfl
  | succ K ih =>
    rw [Function.iterate_succ_apply', fdiff_succ]
    change (bdiff^[K] (fun m : ℤ => φ (m / N))) n - (bdiff^[K] (fun m : ℤ => φ (m / N))) (n - 1) = _
    rw [ih n, ih (n - 1)]
    congr 2
    push_cast; ring

/-- Integers in `[a, b]`. -/
lemma card_Icc_le (a b : ℝ) (hab : a ≤ b) :
    ((Finset.Icc ⌈a⌉ ⌊b⌋ : Finset ℤ).card : ℝ) ≤ b - a + 1 := by
  rw [Int.card_Icc]
  have h1 : (⌊b⌋ : ℝ) ≤ b := Int.floor_le b
  have h2 : a ≤ (⌈a⌉ : ℝ) := Int.le_ceil a
  rcases le_or_gt 0 (⌊b⌋ + 1 - ⌈a⌉) with h | h
  · have : (((⌊b⌋ + 1 - ⌈a⌉).toNat : ℕ) : ℝ) = ((⌊b⌋ + 1 - ⌈a⌉ : ℤ) : ℝ) := by
      rw [show (((⌊b⌋ + 1 - ⌈a⌉).toNat : ℕ) : ℝ) = (((⌊b⌋ + 1 - ⌈a⌉).toNat : ℤ) : ℝ) by norm_cast,
        Int.toNat_of_nonneg h]
    rw [this]; push_cast; linarith
  · have : (⌊b⌋ + 1 - ⌈a⌉).toNat = 0 := Int.toNat_eq_zero.2 h.le
    rw [this]; push_cast; linarith

lemma exists_R_of_compact {φ : ℝ → ℝ} (hc : HasCompactSupport φ) :
    ∃ R : ℝ, 0 ≤ R ∧ ∀ y, R < |y| → φ y = 0 := by
  obtain ⟨R, hR⟩ := (Metric.isBounded_iff_subset_closedBall 0).1 hc.isCompact.isBounded
  refine ⟨max R 0, le_max_right _ _, fun y hy => ?_⟩
  apply image_eq_zero_of_notMem_tsupport
  intro hy'
  have := hR hy'
  rw [Metric.mem_closedBall, dist_zero_right, Real.norm_eq_abs] at this
  linarith [le_max_left R 0]

lemma exists_bound_iteratedDeriv {φ : ℝ → ℝ} (hφ : ContDiff ℝ ∞ φ) (hc : HasCompactSupport φ) (K : ℕ) :
    ∃ M : ℝ, ∀ y, |iteratedDeriv K φ y| ≤ M := by
  have hcs : HasCompactSupport (iteratedDeriv K φ) := by
    have h := hc.iteratedFDeriv (𝕜 := ℝ) K
    have : iteratedDeriv K φ = (fun L => L (fun _ => (1 : ℝ))) ∘ iteratedFDeriv ℝ K φ := by
      funext y; simp [iteratedDeriv_eq_iteratedFDeriv]
    rw [this]
    exact h.comp_left (by simp)
  have hcont : Continuous (iteratedDeriv K φ) := hφ.continuous_iteratedDeriv K (by exact_mod_cast le_top)
  obtain ⟨M, hM⟩ := hcs.exists_bound_of_continuous hcont
  exact ⟨M, fun y => by rw [← Real.norm_eq_abs]; exact hM y⟩

/-- **Decay of smooth exponential sums**: `|Σ_n φ(n/N) e^{inx}| (1 + N|x|)^K ≤ C N` for
`N ≥ 1`, `|x| ≤ π`. -/
theorem esum_decay {φ : ℝ → ℝ} (hφ : ContDiff ℝ ∞ φ) (hc : HasCompactSupport φ) (K : ℕ) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ N : ℝ, 1 ≤ N → ∀ x : ℝ, |x| ≤ Real.pi →
      ‖esum (fun n => φ (n / N)) x‖ * (1 + N * |x|) ^ K ≤ C * N := by
  obtain ⟨R, hR0, hR⟩ := exists_R_of_compact hc
  obtain ⟨M₀, hM₀⟩ := exists_bound_iteratedDeriv hφ hc 0
  obtain ⟨MK, hMK⟩ := exists_bound_iteratedDeriv hφ hc K
  simp only [iteratedDeriv_zero] at hM₀
  have hM₀0 : 0 ≤ M₀ := (abs_nonneg _).trans (hM₀ 0)
  have hMK0 : 0 ≤ MK := (abs_nonneg _).trans (hMK 0)
  set C₀ := (2 * R + 1) * M₀
  set C₁ := (2 * R + K + 1) * MK * (Real.pi / 2) ^ K
  refine ⟨2 ^ K * (C₀ + C₁), by positivity, fun N hN x hx => ?_⟩
  have hN0 : 0 < N := by linarith
  set a : ℤ → ℝ := fun m => φ (m / N) with ha
  -- support of a
  set S₀ : Finset ℤ := Finset.Icc ⌈-R * N⌉ ⌊R * N⌋
  have hS₀ : ∀ n ∉ S₀, a n = 0 := by
    intro n hn
    simp only [S₀, Finset.mem_Icc, not_and_or, not_le] at hn
    apply hR
    rcases hn with h | h
    · have h' : (n : ℝ) < -R * N := Int.lt_ceil.1 h
      have : (n : ℝ) / N < -R := by rw [div_lt_iff₀ hN0]; linarith
      calc R < -((n : ℝ) / N) := by linarith
        _ ≤ |(n : ℝ) / N| := neg_le_abs _
    · have h' : R * N < (n : ℝ) := Int.floor_lt.1 h
      have : R < (n : ℝ) / N := by rw [lt_div_iff₀ hN0]; linarith
      exact lt_of_lt_of_le this (le_abs_self _)
  have hfin : FinSupp a := ⟨S₀, hS₀⟩
  -- support of the K-th difference
  set SK : Finset ℤ := Finset.Icc ⌈-R * N⌉ ⌊R * N + K⌋
  have hSK : ∀ n ∉ SK, (bdiff^[K] a) n = 0 := by
    intro n hn
    rw [ha, bdiff_iter_eq]
    apply fdiff_eq_zero (by positivity) hR K
    simp only [SK, Finset.mem_Icc, not_and_or, not_le] at hn
    rcases hn with h | h
    · have h' : (n : ℝ) < -R * N := Int.lt_ceil.1 h
      left; rw [div_lt_iff₀ hN0]; linarith
    · have h' : R * N + K < (n : ℝ) := Int.floor_lt.1 h
      right; rw [lt_div_iff₀ hN0]
      have : (R + K * (1 / N)) * N = R * N + K := by field_simp
      linarith
  have hcard0 : (S₀.card : ℝ) ≤ (2 * R + 1) * N := by
    have := card_Icc_le (-R * N) (R * N) (by nlinarith)
    nlinarith
  have hcardK : (SK.card : ℝ) ≤ (2 * R + K + 1) * N := by
    have := card_Icc_le (-R * N) (R * N + K) (by nlinarith)
    have hK0 : (0 : ℝ) ≤ K := Nat.cast_nonneg K
    nlinarith
  have hA0 : ‖esum a x‖ ≤ C₀ * N := by
    have := norm_esum_le S₀ hS₀ M₀ (fun n => hM₀ _) x
    calc ‖esum a x‖ ≤ S₀.card * M₀ := this
      _ ≤ (2 * R + 1) * N * M₀ := mul_le_mul_of_nonneg_right hcard0 hM₀0
      _ = C₀ * N := by ring
  have hAK : ‖esum (bdiff^[K] a) x‖ ≤ (2 * R + K + 1) * MK * N * (1 / N) ^ K := by
    have hb : ∀ n, |(bdiff^[K] a) n| ≤ MK * (1 / N) ^ K := by
      intro n; rw [ha, bdiff_iter_eq]
      exact abs_fdiff_le (by positivity) K φ (hφ.of_le (by exact_mod_cast le_top)) MK hMK _
    have := norm_esum_le SK hSK _ hb x
    calc ‖esum (bdiff^[K] a) x‖ ≤ SK.card * (MK * (1 / N) ^ K) := this
      _ ≤ (2 * R + K + 1) * N * (MK * (1 / N) ^ K) :=
          mul_le_mul_of_nonneg_right hcardK (by positivity)
      _ = (2 * R + K + 1) * MK * N * (1 / N) ^ K := by ring
  -- |1 - e^{ix}| ≥ 2|x|/π
  have hsin : 2 * |x| / Real.pi ≤ ‖1 - Complex.exp (I * x)‖ := by
    rw [norm_sub_rev, Complex.norm_exp_I_mul_ofReal_sub_one, Real.norm_eq_abs, abs_mul, abs_two]
    have hx2 : |x / 2| ≤ Real.pi / 2 := by rw [abs_div, abs_two]; linarith
    have hJ := Real.mul_le_sin (abs_nonneg (x / 2)) hx2
    have hs : Real.sin |x / 2| ≤ |Real.sin (x / 2)| := by
      rcases le_total 0 x with h | h
      · rw [abs_of_nonneg (by linarith : (0 : ℝ) ≤ x / 2)]; exact le_abs_self _
      · rw [abs_of_nonpos (by linarith : x / 2 ≤ 0), Real.sin_neg]; exact neg_le_abs _
    have e : |x / 2| = |x| / 2 := by rw [abs_div, abs_two]
    rw [e] at hJ hs
    have e2 : 2 / Real.pi * (|x| / 2) = |x| / Real.pi := by ring
    rw [e2] at hJ
    have : 2 * |x| / Real.pi = 2 * (|x| / Real.pi) := by ring
    rw [this]
    linarith
  have hkey : ‖esum a x‖ * (N * |x|) ^ K ≤ C₁ * N := by
    have h1 := esum_bdiff_iter hfin x K
    have h2 : ‖esum (bdiff^[K] a) x‖ = ‖1 - Complex.exp (I * x)‖ ^ K * ‖esum a x‖ := by
      rw [h1, norm_mul, norm_pow]
    have h3 : (2 * |x| / Real.pi) ^ K * ‖esum a x‖ ≤ (2 * R + K + 1) * MK * N * (1 / N) ^ K := by
      calc (2 * |x| / Real.pi) ^ K * ‖esum a x‖ ≤ ‖1 - Complex.exp (I * x)‖ ^ K * ‖esum a x‖ :=
            mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (by positivity) hsin K) (norm_nonneg _)
        _ ≤ _ := h2 ▸ hAK
    have hNK : N ^ K * (1 / N) ^ K = 1 := by rw [← mul_pow, mul_one_div_cancel hN0.ne', one_pow]
    have e : ‖esum a x‖ * (N * |x|) ^ K = (Real.pi / 2) ^ K * N ^ K * ((2 * |x| / Real.pi) ^ K * ‖esum a x‖) := by
      simp only [mul_pow, div_pow]
      field_simp
    rw [e]
    calc (Real.pi / 2) ^ K * N ^ K * ((2 * |x| / Real.pi) ^ K * ‖esum a x‖)
        ≤ (Real.pi / 2) ^ K * N ^ K * ((2 * R + K + 1) * MK * N * (1 / N) ^ K) :=
          mul_le_mul_of_nonneg_left h3 (by positivity)
      _ = C₁ * N * (N ^ K * (1 / N) ^ K) := by ring
      _ = C₁ * N := by rw [hNK, mul_one]
  -- combine
  have hbin : (1 + N * |x|) ^ K ≤ 2 ^ K * (1 + (N * |x|) ^ K) := by
    have hy : 0 ≤ N * |x| := by positivity
    rcases le_total (N * |x|) 1 with h | h
    · calc (1 + N * |x|) ^ K ≤ 2 ^ K := pow_le_pow_left₀ (by positivity) (by linarith) K
        _ ≤ 2 ^ K * (1 + (N * |x|) ^ K) := by
            have : 0 ≤ (N * |x|) ^ K := by positivity
            nlinarith [pow_pos (by norm_num : (0:ℝ) < 2) K]
    · calc (1 + N * |x|) ^ K ≤ (2 * (N * |x|)) ^ K := pow_le_pow_left₀ (by positivity) (by linarith) K
        _ = 2 ^ K * (N * |x|) ^ K := by rw [mul_pow]
        _ ≤ 2 ^ K * (1 + (N * |x|) ^ K) := by
            nlinarith [pow_pos (by norm_num : (0:ℝ) < 2) K]
  calc ‖esum a x‖ * (1 + N * |x|) ^ K ≤ ‖esum a x‖ * (2 ^ K * (1 + (N * |x|) ^ K)) :=
        mul_le_mul_of_nonneg_left hbin (norm_nonneg _)
    _ = 2 ^ K * (‖esum a x‖ + ‖esum a x‖ * (N * |x|) ^ K) := by ring
    _ ≤ 2 ^ K * (C₀ * N + C₁ * N) := by gcongr
    _ = 2 ^ K * (C₀ + C₁) * N := by ring

end LeanLi
