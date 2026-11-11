import Riemann.Deep.LiRate.Lean.Core
import Riemann.Deep.LiRate.Lean.Cor

/-!
# LiRate Theorem 1 from named inputs (proofs.md Theorem 1, step 7: constants, `u = 3/5`)

Inputs (hypotheses of `thm1`):
* (H_Li) Li's formula with real phases (RH): `λ_n = Σ_k 2(1 - cos nθ_k)`;
* (H_N) Riemann–von Mangoldt `RvM`: `#Z(T) = (T/2π) log(T/2πe) + O(log T)`;
* (H_w) for every `η > 0` a weight family with `0 ≤ φ ≤ 1`, `Φ ≥ 0`, `Φ(0) ≥ (1-η)N` (`GoodWeights`);
* (H_τ) for the same weights, the trend/tail estimate of steps 4, 5, 7 (`TrendTail`):
  `τ' ≤ (1+η) T³L²/(6π)`, `L = log(T/2π)`, uniformly in `√N/log N ≤ T ≤ √N/2`, `N ≥ N₀(η)`.
The proof chooses `T = √(6N/(5 log N))`, i.e. `T² · (log N)/2 = (3/5) N` (the optimum `u = 3/5`).
-/

open Finset Real Filter

namespace LiRate

/-- (H_N) Riemann–von Mangoldt for `Z T` = the zeros with `0 < γ ≤ T` (with multiplicity). -/
def RvM {κ : Type*} (Z : ℝ → Finset κ) : Prop :=
  ∃ C T₀ : ℝ, ∀ T ≥ T₀,
    |((Z T).card : ℝ) - T / (2 * π) * Real.log (T / (2 * π * exp 1))| ≤ C * Real.log T

/-- (H_w) weights `w N n = φ(n/N)` of quality `η`. -/
def GoodWeights (w : ℕ → ℕ → ℝ) (η : ℝ) : Prop :=
  (∀ N n, 0 ≤ w N n) ∧ (∀ N n, w N n ≤ 1) ∧ (∀ N x, 0 ≤ Phi N (w N) x) ∧
    ∀ᶠ N : ℕ in atTop, (1 - η) * N ≤ Phi N (w N) 0

/-- (H_τ) trend and tails (proofs.md steps 4, 5 and 7: `τ' = τ_T + O(N(T) log N) + tiny`,
`τ_T ~ T³L²/(6π)`). -/
def TrendTail {κ : Type*} (w : ℕ → ℕ → ℝ) (θ : κ → ℝ) (lam t : ℕ → ℝ) (Z : ℝ → Finset κ)
    (η : ℝ) : Prop :=
  ∀ᶠ N : ℕ in atTop, ∀ T : ℝ, √N / Real.log N ≤ T → T ≤ √N / 2 →
    tauPrime N (w N) θ lam t (Z T) ≤ (1 + η) * (T ^ 3 * Real.log (T / (2 * π)) ^ 2 / (6 * π))

/-! ## Pure inequalities -/

/-- `(M - τ)²/M ≥ M₀ (1 - v)²` for `M ≥ M₀ > 0`, `0 ≤ τ ≤ v M₀`, `v < 1`. -/
lemma opt (M M₀ τ v : ℝ) (hM₀ : 0 < M₀) (hM : M₀ ≤ M) (hτ0 : 0 ≤ τ)
    (hv : v < 1) (hτ : τ ≤ v * M₀) : M₀ * (1 - v) ^ 2 ≤ (M - τ) ^ 2 / M := by
  have hMpos : 0 < M := lt_of_lt_of_le hM₀ hM
  rw [le_div_iff₀ hMpos]
  have h1 : 0 ≤ M - v * M₀ := by nlinarith
  have h2 : M - v * M₀ ≤ M - τ := by linarith
  have h3 : (M - v * M₀) ^ 2 ≤ (M - τ) ^ 2 := pow_le_pow_left₀ h1 h2 2
  have hv2 : v ^ 2 * M₀ ≤ M₀ := by nlinarith
  have h4 : M₀ * (1 - v) ^ 2 * M ≤ (M - v * M₀) ^ 2 := by
    nlinarith [mul_nonneg (sub_nonneg.2 hM) (show 0 ≤ M - v ^ 2 * M₀ by linarith)]
  linarith

/-- `(1-η)³ (1 - (1+η)/(5(1-η)²))² ≥ 16/25 - (16/5) η` for `0 < η ≤ 1/10`. -/
lemma const_bound (η : ℝ) (hη0 : 0 < η) (hη1 : η ≤ 1 / 10) :
    16 / 25 - 16 / 5 * η ≤ (1 - η) ^ 3 * (1 - (1 + η) / (5 * (1 - η) ^ 2)) ^ 2 := by
  have h1η : 0 < 1 - η := by linarith
  have hv : (1 + η) / (5 * (1 - η) ^ 2) ≤ 1 / 5 + 3 / 4 * η := by
    rw [div_le_iff₀ (mul_pos (by norm_num) (pow_pos h1η 2))]
    nlinarith [mul_nonneg hη0.le (by linarith : (0 : ℝ) ≤ 3 / 4 - 13 / 2 * η), pow_nonneg hη0.le 3]
  have hs0 : 0 ≤ 4 / 5 - 3 / 4 * η := by linarith
  have hsq : (4 / 5 - 3 / 4 * η) ^ 2 ≤ (1 - (1 + η) / (5 * (1 - η) ^ 2)) ^ 2 :=
    pow_le_pow_left₀ hs0 (by linarith) 2
  have hc : 1 - 3 * η ≤ (1 - η) ^ 3 := by nlinarith [sq_nonneg η, pow_nonneg hη0.le 3]
  calc 16 / 25 - 16 / 5 * η ≤ (1 - 3 * η) * (4 / 5 - 3 / 4 * η) ^ 2 := by
        nlinarith [mul_nonneg (sq_nonneg η) (by linarith : (0 : ℝ) ≤ 333 / 80 - 27 / 16 * η)]
    _ ≤ (1 - η) ^ 3 * (1 - (1 + η) / (5 * (1 - η) ^ 2)) ^ 2 :=
        mul_le_mul hc hsq (sq_nonneg _) (pow_nonneg h1η.le 3)

/-- Step 7 at one `N`: from the core bound and the lower bounds on `#Z(T)`, `Φ(0)`, `L`. -/
lemma step7 {Nr T L L₀ Zc Φ0 η S : ℝ} (hη0 : 0 < η) (hη1 : η ≤ 1 / 10)
    (hNr : 0 < Nr) (hT : 0 < T) (hL₀ : 0 < L₀) (hTL : T ^ 2 * L₀ = 3 / 5 * Nr)
    (hL : (1 - η) * L₀ ≤ L) (hLu : L ≤ L₀)
    (hZ : (1 - η) * (T / (2 * π) * L) ≤ Zc) (hΦ : (1 - η) * Nr ≤ Φ0)
    (hcore : (1 + η) * (T ^ 3 * L ^ 2 / (6 * π)) < Zc * Φ0 →
      (Zc * Φ0 - (1 + η) * (T ^ 3 * L ^ 2 / (6 * π))) ^ 2 / (Zc * Φ0) ≤ S) :
    (16 / 25 - 16 / 5 * η) * (T * L₀ * Nr / (2 * π)) ≤ S := by
  have hπ : 0 < π := Real.pi_pos
  have h2π : 0 < 2 * π := by positivity
  have h1η : 0 < 1 - η := by linarith
  have h1η' : 0 < 1 + η := by linarith
  have hT2π : 0 < T / (2 * π) := div_pos hT h2π
  have hLpos : 0 < L := lt_of_lt_of_le (mul_pos h1η hL₀) hL
  set v : ℝ := (1 + η) / (5 * (1 - η) ^ 2) with hv
  set τ : ℝ := (1 + η) * (T ^ 3 * L ^ 2 / (6 * π)) with hτdef
  set M₀ : ℝ := (1 - η) ^ 2 * (T / (2 * π) * L * Nr) with hM₀def
  have hM₀ : 0 < M₀ := mul_pos (pow_pos h1η 2) (mul_pos (mul_pos hT2π hLpos) hNr)
  have hZ0 : 0 ≤ (1 - η) * (T / (2 * π) * L) := (mul_pos h1η (mul_pos hT2π hLpos)).le
  have hM : M₀ ≤ Zc * Φ0 := by
    have := mul_le_mul hZ hΦ (mul_pos h1η hNr).le (hZ0.trans hZ)
    calc M₀ = (1 - η) * (T / (2 * π) * L) * ((1 - η) * Nr) := by rw [hM₀def]; ring
      _ ≤ Zc * Φ0 := this
  have hτ0 : 0 ≤ τ :=
    (mul_pos h1η' (div_pos (mul_pos (pow_pos hT 3) (pow_pos hLpos 2)) (by positivity))).le
  have h5 : 0 < 5 * (1 - η) ^ 2 := mul_pos (by norm_num) (pow_pos h1η 2)
  have hv1 : v < 1 := by
    rw [hv, div_lt_one h5]
    have h81 : (9 / 10 : ℝ) ^ 2 ≤ (1 - η) ^ 2 := pow_le_pow_left₀ (by norm_num) (by linarith) 2
    norm_num at h81
    linarith
  have hkey : T ^ 2 * L ≤ 3 / 5 * Nr := by
    rw [← hTL]; exact mul_le_mul_of_nonneg_left hLu (sq_nonneg T)
  have h1ηne : 1 - η ≠ 0 := h1η.ne'
  have hπne : π ≠ 0 := hπ.ne'
  have hτv : τ ≤ v * M₀ := by
    have e1 : v * M₀ = (1 + η) * T * L / (2 * π) * (Nr / 5) := by
      rw [hv, hM₀def]; field_simp
    have e2 : τ = (1 + η) * T * L / (2 * π) * (T ^ 2 * L / 3) := by
      rw [hτdef]; field_simp; ring
    rw [e1, e2]
    apply mul_le_mul_of_nonneg_left _ (div_pos (mul_pos (mul_pos h1η' hT) hLpos) h2π).le
    linarith
  have hτM : τ < Zc * Φ0 := by
    have := mul_lt_mul_of_pos_right hv1 hM₀
    linarith
  have hopt := opt (Zc * Φ0) M₀ τ v hM₀ hM hτ0 hv1 hτv
  have hS := hcore hτM
  have hM₀' : (1 - η) ^ 3 * (T * L₀ * Nr / (2 * π)) ≤ M₀ := by
    have e3 : (1 - η) ^ 3 * (T * L₀ * Nr / (2 * π)) =
        (1 - η) ^ 2 * (T * Nr / (2 * π)) * ((1 - η) * L₀) := by ring
    have e4 : M₀ = (1 - η) ^ 2 * (T * Nr / (2 * π)) * L := by rw [hM₀def]; ring
    rw [e3, e4]
    exact mul_le_mul_of_nonneg_left hL (mul_pos (pow_pos h1η 2) (div_pos (mul_pos hT hNr) h2π)).le
  have hc := const_bound η hη0 hη1
  have hX : 0 ≤ T * L₀ * Nr / (2 * π) := (div_pos (mul_pos (mul_pos hT hL₀) hNr) h2π).le
  calc (16 / 25 - 16 / 5 * η) * (T * L₀ * Nr / (2 * π))
      ≤ (1 - η) ^ 3 * (1 - v) ^ 2 * (T * L₀ * Nr / (2 * π)) :=
        mul_le_mul_of_nonneg_right hc hX
    _ = (1 - v) ^ 2 * ((1 - η) ^ 3 * (T * L₀ * Nr / (2 * π))) := by ring
    _ ≤ (1 - v) ^ 2 * M₀ := mul_le_mul_of_nonneg_left hM₀' (sq_nonneg _)
    _ = M₀ * (1 - v) ^ 2 := by ring
    _ ≤ (Zc * Φ0 - τ) ^ 2 / (Zc * Φ0) := hopt
    _ ≤ S := hS

/-! ## The choice `T = √(6N/(5 log N))` and its asymptotics -/

/-- `T(N) = √(6/5) √N / √(log N)`. -/
noncomputable def Tch (N : ℕ) : ℝ := √(6 / 5) * √(N : ℝ) / √(Real.log N)

/-- `k = √(6/5)/(4π)`, so that `T · (log N / 2) · N / (2π) = k N^{3/2} (log N)^{1/2}`. -/
noncomputable def kk : ℝ := √(6 / 5) / (4 * π)

lemma c₁_eq : c₁ = 16 / 25 * kk := by
  have h : √(2 : ℝ) * √(6 / 5) = 2 * √(3 / 5) := by
    rw [← Real.sqrt_mul (by norm_num), show (2 : ℝ) * (6 / 5) = 2 ^ 2 * (3 / 5) by norm_num,
      Real.sqrt_mul (by norm_num), Real.sqrt_sq (by norm_num)]
  have h3 : √(3 / 5 : ℝ) = √2 * √(6 / 5) / 2 := by rw [h]; ring
  have h2 : (√2 : ℝ) ≠ 0 := by positivity
  have hπ : π ≠ 0 := Real.pi_ne_zero
  unfold c₁ kk
  rw [h3]
  field_simp
  ring

lemma kk_le : kk ≤ 1 / 10 := by
  have h1 : √(6 / 5 : ℝ) ≤ 1.1 := by rw [Real.sqrt_le_left (by norm_num)]; norm_num
  unfold kk
  rw [div_le_iff₀ (by positivity)]
  nlinarith [Real.pi_gt_three]

/-- Basic asymptotics in `N`: `log N → ∞`, `log log N = o(log N)`, `log N = o(√N)`. -/
lemma eventually_basic (A c' η : ℝ) (hc' : 0 < c') (hη : 0 < η) :
    ∀ᶠ N : ℕ in atTop, A ≤ Real.log N ∧ Real.log (Real.log N) ≤ η / 4 * Real.log N ∧
      Real.log N ≤ c' * √(N : ℝ) := by
  have h1 : ∀ᶠ x : ℝ in atTop, A ≤ Real.log x := Real.tendsto_log_atTop.eventually_ge_atTop A
  have h2 : ∀ᶠ x : ℝ in atTop, Real.log (Real.log x) ≤ η / 4 * Real.log x := by
    have hlo := (Real.isLittleO_log_id_atTop.comp_tendsto Real.tendsto_log_atTop).bound
      (show 0 < η / 4 by linarith)
    filter_upwards [hlo, Real.tendsto_log_atTop.eventually_ge_atTop 0] with x hx hx0
    simp only [Function.comp, id, Real.norm_eq_abs] at hx
    rw [abs_of_nonneg hx0] at hx
    exact (le_abs_self _).trans hx
  have h3 : ∀ᶠ x : ℝ in atTop, Real.log x ≤ c' * √x := by
    have hlo := (isLittleO_log_rpow_atTop (by norm_num : (0 : ℝ) < 1 / 2)).bound hc'
    filter_upwards [hlo, eventually_gt_atTop 0] with x hx hx0
    rw [Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos (Real.rpow_pos_of_pos hx0 _),
      ← Real.sqrt_eq_rpow] at hx
    exact (le_abs_self _).trans hx
  exact tendsto_natCast_atTop_atTop.eventually (h1.and (h2.and h3))

set_option maxHeartbeats 1000000 in
/-- The facts about `T = √(6/5) √x/√ℓ`, `L = log(T/2π)` used in step 7, at one `x`. -/
lemma facts_at {x ℓ η C T₀ T L : ℝ} (hη0 : 0 < η) (hη1 : η ≤ 1 / 10) (hx : 1 ≤ x)
    (hℓ : ℓ = Real.log x) (h5 : 5 ≤ ℓ) (h4 : 4 / η ≤ ℓ) (hpi : 4 * Real.log (2 * π) / η ≤ ℓ)
    (hll : Real.log ℓ ≤ η / 4 * ℓ) (hs1 : ℓ ≤ 1 / (|T₀| + 1) * √x)
    (hs2 : ℓ ≤ η / (5 * π * (|C| + 1)) * √x)
    (hT : T = √(6 / 5) * √x / √ℓ) (hL : L = Real.log (T / (2 * π))) :
    0 < T ∧ √x / ℓ ≤ T ∧ T ≤ √x / 2 ∧ T₀ ≤ T ∧ T ^ 2 * (ℓ / 2) = 3 / 5 * x ∧
      (1 - η) * (ℓ / 2) ≤ L ∧ L ≤ ℓ / 2 ∧ C * Real.log T ≤ T / (2 * π) * (η * L - 1) := by
  have hℓ0 : 0 < ℓ := by linarith
  have hℓne : ℓ ≠ 0 := hℓ0.ne'
  have hx0 : 0 < x := by linarith
  have hsx : 0 < √x := Real.sqrt_pos.2 hx0
  have hsℓ : 0 < √ℓ := Real.sqrt_pos.2 hℓ0
  have hs65 : 1 ≤ √(6 / 5 : ℝ) := Real.one_le_sqrt.2 (by norm_num)
  have hs65' : √(6 / 5 : ℝ) ≤ 1.1 := by rw [Real.sqrt_le_left (by norm_num)]; norm_num
  have hsq65 : √(6 / 5 : ℝ) ^ 2 = 6 / 5 := Real.sq_sqrt (by norm_num)
  have h2π : 0 < 2 * π := by positivity
  have hTpos : 0 < T := by rw [hT]; exact div_pos (mul_pos (by positivity) hsx) hsℓ
  have hT1 : √x / ℓ ≤ T := by
    rw [hT, div_le_div_iff₀ hℓ0 hsℓ]
    have h : √ℓ ≤ ℓ := by
      rw [Real.sqrt_le_left hℓ0.le, sq]; exact le_mul_of_one_le_right hℓ0.le (by linarith)
    linarith [mul_le_mul_of_nonneg_left h hsx.le,
      mul_le_mul_of_nonneg_right hs65 (mul_nonneg hsx.le hℓ0.le)]
  have hT2 : T ≤ √x / 2 := by
    rw [hT, div_le_div_iff₀ hsℓ (by norm_num)]
    have h : (2.2 : ℝ) ≤ √ℓ := by
      rw [Real.le_sqrt (by norm_num) hℓ0.le]; norm_num; linarith
    linarith [mul_le_mul_of_nonneg_left (show √(6 / 5) * 2 ≤ √ℓ by linarith) hsx.le]
  have hb1 : (|T₀| + 1) * ℓ ≤ √x := by
    have h := mul_le_mul_of_nonneg_left hs1 (by positivity : (0 : ℝ) ≤ |T₀| + 1)
    have e : (|T₀| + 1) * (1 / (|T₀| + 1) * √x) = √x := by field_simp
    linarith
  have hbig : |T₀| + 1 ≤ √x / ℓ := by rw [le_div_iff₀ hℓ0]; linarith
  have hT0 : T₀ ≤ T := by linarith [le_abs_self T₀]
  have hT1' : 1 ≤ T := by linarith [abs_nonneg T₀]
  have hTL : T ^ 2 * (ℓ / 2) = 3 / 5 * x := by
    rw [hT, div_pow, mul_pow, hsq65, Real.sq_sqrt hx0.le, Real.sq_sqrt hℓ0.le]
    field_simp
    ring
  have hLu : L ≤ ℓ / 2 := by
    rw [hL, hℓ, ← Real.log_sqrt hx0.le]
    apply Real.log_le_log (div_pos hTpos h2π)
    rw [div_le_iff₀ h2π]
    linarith [mul_le_mul_of_nonneg_left
      (show (1 : ℝ) / 2 ≤ 2 * π by linarith [Real.pi_gt_three]) hsx.le]
  have hLl : (1 - η) * (ℓ / 2) ≤ L := by
    have hlow : Real.log (√x / (2 * π * ℓ)) ≤ L := by
      rw [hL]
      apply Real.log_le_log (div_pos hsx (mul_pos h2π hℓ0))
      rw [div_le_div_iff₀ (mul_pos h2π hℓ0) h2π]
      have h := hT1
      rw [div_le_iff₀ hℓ0] at h
      linarith [mul_le_mul_of_nonneg_right h h2π.le]
    have heq : Real.log (√x / (2 * π * ℓ)) = ℓ / 2 - Real.log (2 * π) - Real.log ℓ := by
      rw [Real.log_div hsx.ne' (mul_pos h2π hℓ0).ne', Real.log_mul h2π.ne' hℓne,
        Real.log_sqrt hx0.le, ← hℓ]
      ring
    rw [heq] at hlow
    have h : Real.log (2 * π) ≤ η / 4 * ℓ := by
      rw [div_le_iff₀ hη0] at hpi; linarith
    linarith
  have hRv : C * Real.log T ≤ T / (2 * π) * (η * L - 1) := by
    have hlogT0 : 0 ≤ Real.log T := Real.log_nonneg hT1'
    have hlogT : Real.log T ≤ ℓ / 2 := by
      rw [hℓ, ← Real.log_sqrt hx0.le]; exact Real.log_le_log hTpos (by linarith)
    have hηL : η * ℓ / 5 ≤ η * L - 1 := by
      have h : 1 ≤ η * ℓ / 4 := by rw [div_le_iff₀ hη0] at h4; linarith
      linarith [mul_le_mul_of_nonneg_left hLl hη0.le,
        mul_nonneg (mul_nonneg hη0.le hℓ0.le) (by linarith : (0 : ℝ) ≤ 1 / 10 - η)]
    have hA : C * Real.log T ≤ |C| * (ℓ / 2) :=
      calc C * Real.log T ≤ |C| * Real.log T := mul_le_mul_of_nonneg_right (le_abs_self C) hlogT0
        _ ≤ |C| * (ℓ / 2) := mul_le_mul_of_nonneg_left hlogT (abs_nonneg C)
    have hb2 : 5 * π * (|C| + 1) * ℓ ≤ η * √x := by
      have h := mul_le_mul_of_nonneg_left hs2 (by positivity : (0 : ℝ) ≤ 5 * π * (|C| + 1))
      have e : 5 * π * (|C| + 1) * (η / (5 * π * (|C| + 1)) * √x) = η * √x := by field_simp
      linarith
    have hB : |C| * (ℓ / 2) ≤ η * √x / (10 * π) := by
      rw [le_div_iff₀ (by positivity)]
      linarith [mul_pos Real.pi_pos hℓ0]
    have hT2π' : √x / (2 * π * ℓ) ≤ T / (2 * π) := by
      rw [div_le_div_iff₀ (mul_pos h2π hℓ0) h2π]
      have h := hT1
      rw [div_le_iff₀ hℓ0] at h
      linarith [mul_le_mul_of_nonneg_right h h2π.le]
    have hD : η * √x / (10 * π) ≤ T / (2 * π) * (η * L - 1) :=
      calc η * √x / (10 * π) = √x / (2 * π * ℓ) * (η * ℓ / 5) := by field_simp; ring
        _ ≤ T / (2 * π) * (η * L - 1) :=
          mul_le_mul hT2π' hηL (div_nonneg (mul_nonneg hη0.le hℓ0.le) (by norm_num))
            (div_pos hTpos h2π).le
    linarith
  exact ⟨hTpos, hT1, hT2, hT0, hTL, hLl, hLu, hRv⟩

/-- **Theorem 1** (RH) from (H_Li), (H_N), (H_w), (H_τ):
`∀ ε > 0`, eventually `Σ_{n≤N} R(n)² ≥ (c₁ - ε) N^{3/2} (log N)^{1/2}`, `R(n) = λ_n - t(n)`. -/
theorem thm1 {κ : Type*} (θ : κ → ℝ) (lam t : ℕ → ℝ)
    (hlam : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)) (ht0 : t 0 = 0)
    (Z : ℝ → Finset κ) (hN : RvM Z)
    (hW : ∀ η > 0, ∃ w : ℕ → ℕ → ℝ, GoodWeights w η ∧ TrendTail w θ lam t Z η) :
    Thm1Bound (fun n => lam n - t n) := by
  intro ε hε
  obtain ⟨C, T₀, hC⟩ := hN
  obtain ⟨η, hη0, hη1, hηε⟩ : ∃ η : ℝ, 0 < η ∧ η ≤ 1 / 10 ∧ η ≤ ε :=
    ⟨min (1 / 10) ε, lt_min (by norm_num) hε, min_le_left _ _, min_le_right _ _⟩
  obtain ⟨w, ⟨hw0, hw1, hΦ, hΦ0⟩, hTT⟩ := hW η hη0
  have hc' : 0 < min (1 / (|T₀| + 1)) (η / (5 * π * (|C| + 1))) :=
    lt_min (by positivity) (div_pos hη0 (by positivity))
  filter_upwards [hΦ0, hTT, eventually_basic (max (max 5 (4 / η)) (4 * Real.log (2 * π) / η))
    _ η hc' hη0, eventually_ge_atTop 1] with N h1 h2 h3 h4
  obtain ⟨hA, hll, hsq⟩ := h3
  have hx : (1 : ℝ) ≤ N := by exact_mod_cast h4
  obtain ⟨hTpos, hT1, hT2, hT0, hTL, hLl, hLu, hRv⟩ :=
    facts_at (x := N) (ℓ := Real.log N) (T := Tch N) (L := Real.log (Tch N / (2 * π)))
      (C := C) (T₀ := T₀) hη0 hη1 hx rfl
      (le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hA)
      (le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hA)
      (le_trans (le_max_right _ _) hA) hll
      (hsq.trans (mul_le_mul_of_nonneg_right (min_le_left _ _) (Real.sqrt_nonneg _)))
      (hsq.trans (mul_le_mul_of_nonneg_right (min_le_right _ _) (Real.sqrt_nonneg _))) rfl rfl
  have hℓ5 : 5 ≤ Real.log N := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hA
  have hℓ0 : 0 < Real.log N := by linarith
  have hZ : (1 - η) * (Tch N / (2 * π) * Real.log (Tch N / (2 * π))) ≤ ((Z (Tch N)).card : ℝ) := by
    have h := hC (Tch N) hT0
    have hlog : Real.log (Tch N / (2 * π * exp 1)) = Real.log (Tch N / (2 * π)) - 1 := by
      rw [show Tch N / (2 * π * exp 1) = Tch N / (2 * π) / exp 1 by ring,
        Real.log_div (div_pos hTpos (by positivity)).ne' (Real.exp_pos 1).ne', Real.log_exp]
    rw [hlog] at h
    have := (abs_le.1 h).1
    linarith [hRv]
  have hτ0 : 0 ≤ (1 + η) * (Tch N ^ 3 * Real.log (Tch N / (2 * π)) ^ 2 / (6 * π)) :=
    mul_nonneg (by linarith) (div_nonneg (mul_nonneg (pow_nonneg hTpos.le 3) (sq_nonneg _))
      (by positivity))
  have hS := step7 (Nr := N) (S := ∑ n ∈ Icc 1 N, (lam n - t n) ^ 2) hη0 hη1 (by linarith)
    hTpos (by linarith) hTL hLl hLu hZ h1
    (fun hτM => core_bound N (w N) (hw0 N) (hw1 N) (hΦ N) θ lam t hlam ht0 (Z (Tch N)) _
      (h2 (Tch N) hT1 hT2) hτ0 hτM)
  have hsl : √(Real.log N) * √(Real.log N) = Real.log N := Real.mul_self_sqrt hℓ0.le
  have hsl0 : √(Real.log N) ≠ 0 := (Real.sqrt_pos.2 hℓ0).ne'
  have e : Tch N * (Real.log N / 2) * N / (2 * π) = kk * (N * √N * √(Real.log N)) := by
    unfold Tch kk
    calc √(6 / 5) * √(N : ℝ) / √(Real.log N) * (Real.log N / 2) * N / (2 * π) =
        √(6 / 5) * √(N : ℝ) / √(Real.log N) * ((√(Real.log N) * √(Real.log N)) / 2) * N /
          (2 * π) := by rw [hsl]
      _ = √(6 / 5) / (4 * π) * (N * √N * √(Real.log N)) := by field_simp; ring
  rw [e] at hS
  have hk := kk_le
  have hX : 0 ≤ (N : ℝ) * √N * √(Real.log N) := by positivity
  have hmain : 16 / 25 * kk - ε ≤ (16 / 25 - 16 / 5 * η) * kk := by
    linarith [mul_le_mul_of_nonneg_left hk hη0.le]
  rw [c₁_eq]
  calc (16 / 25 * kk - ε) * (N * √N * √(Real.log N))
      ≤ (16 / 25 - 16 / 5 * η) * kk * (N * √N * √(Real.log N)) :=
        mul_le_mul_of_nonneg_right hmain hX
    _ = (16 / 25 - 16 / 5 * η) * (kk * (N * √N * √(Real.log N))) := by ring
    _ ≤ _ := hS

end LiRate
