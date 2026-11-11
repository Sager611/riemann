import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Algebra.QuadraticDiscriminant
import Mathlib.Topology.Algebra.InfiniteSum.Real

/-!
# LiRate, Theorem 1 core: Fourier positivity + Cauchy–Schwarz (proofs.md, Theorem 1, steps 1–6)

Everything is written for even sequences on `ℤ`, stored on `ℕ`:
`Σ_{|n|≤N} w(|n|) F(|n|) G(|n|) = w 0 F 0 G 0 + 2 Σ_{n=1}^N w n F n G n`.
`w n` plays `φ(n/N)`, `Phi N w x = Σ_{|n|≤N} φ(n/N) e^{inx}` (proofs.md `Φ`).
The zeros are an arbitrary family `θ : κ → ℝ` of real phases (RH: `1 - 1/ρ = e^{iθ_γ}`), and Li's
formula is the hypothesis `λ_n = Σ_k 2(1 - cos nθ_k)` (`HasSum`).
-/

open Finset Real

namespace LiRate

/-- `Φ(x) = Σ_{|n|≤N} w(|n|) e^{inx}`. -/
noncomputable def Phi (N : ℕ) (w : ℕ → ℝ) (x : ℝ) : ℝ :=
  w 0 + 2 * ∑ n ∈ Icc 1 N, w n * cos (n * x)

/-- `⟨F, G⟩ = Σ_{|n|≤N} w(|n|) F(|n|) G(|n|)` (proofs.md step 2). -/
noncomputable def ip (N : ℕ) (w : ℕ → ℝ) (F G : ℕ → ℝ) : ℝ :=
  w 0 * F 0 * G 0 + 2 * ∑ n ∈ Icc 1 N, w n * F n * G n

/-- `B(θ, θ') = 2Φ(θ) - Φ(θ - θ') - Φ(θ + θ')` (proofs.md step 3). -/
noncomputable def B (N : ℕ) (w : ℕ → ℝ) (a b : ℝ) : ℝ :=
  2 * Phi N w a - Phi N w (a - b) - Phi N w (a + b)

/-- Dual vector `W(n) = Σ_{γ ≤ T} cos nθ_γ`, `J` = the zeros with `γ ≤ T` (proofs.md step 5). -/
noncomputable def Wfun {κ : Type*} (θ : κ → ℝ) (J : Finset κ) : ℕ → ℝ :=
  fun n => ∑ j ∈ J, cos (n * θ j)

/-- Trend transform `τ(θ) = Σ_{|n|≤N} φ(n/N) e^{-inθ} t(|n|)` (proofs.md step 4). -/
noncomputable def trend (N : ℕ) (w : ℕ → ℝ) (t : ℕ → ℝ) (a : ℝ) : ℝ :=
  ip N w (fun n => cos (n * a)) t

/-- Contribution of the zeros outside `J` to `Σ_n φ(n/N) e^{-inθ} λ_{|n|}`:
`Σ_{k ∉ J} B(θ, θ_k)` (see `tail_eq_tsum_compl`; proofs.md step 5, `γ' > T`). -/
noncomputable def tail {κ : Type*} (N : ℕ) (w : ℕ → ℝ) (θ : κ → ℝ) (lam : ℕ → ℝ)
    (J : Finset κ) (a : ℝ) : ℝ :=
  ip N w (fun n => cos (n * a)) lam - ∑ k ∈ J, B N w a (θ k)

/-- The exact residual `τ'` with `⟨R̃, W⟩ = τ' - 2‖W‖²` (`ip_R_W`). proofs.md's `τ'` (step 5,
line 53) is an upper bound for it, obtained by dropping nonpositive terms (`tail_le`). -/
noncomputable def tauPrime {κ : Type*} (N : ℕ) (w : ℕ → ℝ) (θ : κ → ℝ) (lam t : ℕ → ℝ)
    (J : Finset κ) : ℝ :=
  2 * J.card * ∑ j ∈ J, Phi N w (θ j) + ∑ j ∈ J, tail N w θ lam J (θ j)
    - ∑ j ∈ J, trend N w t (θ j)

/-! ## Bilinear algebra of `ip` -/

lemma ip_comm (N : ℕ) (w F G : ℕ → ℝ) : ip N w F G = ip N w G F := by
  unfold ip
  have h : ∀ n ∈ Icc 1 N, w n * F n * G n = w n * G n * F n := fun n _ => by ring
  rw [sum_congr rfl h]; ring

lemma ip_sub_right (N : ℕ) (w F G H : ℕ → ℝ) :
    ip N w F (fun n => G n - H n) = ip N w F G - ip N w F H := by
  unfold ip
  have h : ∀ n ∈ Icc 1 N, w n * F n * (G n - H n) = w n * F n * G n - w n * F n * H n :=
    fun n _ => by ring
  rw [sum_congr rfl h, sum_sub_distrib]; ring

lemma ip_sum_right {ι : Type*} (N : ℕ) (w F : ℕ → ℝ) (J : Finset ι) (G : ι → ℕ → ℝ) :
    ip N w F (fun n => ∑ j ∈ J, G j n) = ∑ j ∈ J, ip N w F (G j) := by
  simp only [ip, mul_sum, sum_add_distrib]
  rw [sum_comm (s := Icc 1 N)]

lemma ip_self_nonneg (N : ℕ) (w : ℕ → ℝ) (hw : ∀ n, 0 ≤ w n) (F : ℕ → ℝ) : 0 ≤ ip N w F F := by
  unfold ip
  have h1 : 0 ≤ w 0 * F 0 * F 0 := by rw [mul_assoc]; exact mul_nonneg (hw 0) (mul_self_nonneg _)
  have h2 : 0 ≤ ∑ n ∈ Icc 1 N, w n * F n * F n :=
    sum_nonneg fun n _ => by rw [mul_assoc]; exact mul_nonneg (hw n) (mul_self_nonneg _)
  linarith

lemma ip_quad (N : ℕ) (w F G : ℕ → ℝ) (s : ℝ) :
    ip N w (fun n => F n + s * G n) (fun n => F n + s * G n) =
      ip N w G G * (s * s) + 2 * ip N w F G * s + ip N w F F := by
  unfold ip
  have h : ∀ n ∈ Icc 1 N, w n * (F n + s * G n) * (F n + s * G n) =
      w n * F n * F n + (2 * s) * (w n * F n * G n) + (s * s) * (w n * G n * G n) :=
    fun n _ => by ring
  rw [sum_congr rfl h, sum_add_distrib, sum_add_distrib, ← mul_sum, ← mul_sum]
  ring

/-- Weighted Cauchy–Schwarz (proofs.md step 6): `⟨F,G⟩² ≤ ‖F‖²‖G‖²` for `w ≥ 0`. -/
lemma ip_sq_le (N : ℕ) (w : ℕ → ℝ) (hw : ∀ n, 0 ≤ w n) (F G : ℕ → ℝ) :
    ip N w F G ^ 2 ≤ ip N w F F * ip N w G G := by
  have h : ∀ s : ℝ, 0 ≤ ip N w G G * (s * s) + 2 * ip N w F G * s + ip N w F F := by
    intro s; rw [← ip_quad]; exact ip_self_nonneg N w hw _
  have := discrim_le_zero h
  unfold discrim at this
  nlinarith [this]

/-! ## Trigonometric identities -/

lemma ip_cos_cos (N : ℕ) (w : ℕ → ℝ) (a b : ℝ) :
    ip N w (fun n => cos (n * a)) (fun n => cos (n * b)) =
      (Phi N w (a - b) + Phi N w (a + b)) / 2 := by
  have h : ∀ n ∈ Icc 1 N, w n * cos (n * a) * cos (n * b) =
      1 / 2 * (w n * cos (n * (a - b))) + 1 / 2 * (w n * cos (n * (a + b))) := by
    intro n _
    rw [mul_sub, mul_add, cos_sub, cos_add]; ring
  simp only [ip, Phi, Nat.cast_zero, zero_mul, cos_zero]
  rw [sum_congr rfl h, sum_add_distrib, ← mul_sum, ← mul_sum]
  ring

/-- Li pairing (proofs.md step 3): `Σ_n φ(n/N) e^{-inθ} · 2(1 - cos nθ') = B(θ, θ')`. -/
lemma li_pair (N : ℕ) (w : ℕ → ℝ) (a b : ℝ) :
    ip N w (fun n => cos (n * a)) (fun n => 2 * (1 - cos (n * b))) = B N w a b := by
  have h : ∀ n ∈ Icc 1 N, w n * cos (n * a) * (2 * (1 - cos (n * b))) =
      2 * (w n * cos (n * a)) - (w n * cos (n * (a - b)) + w n * cos (n * (a + b))) := by
    intro n _
    rw [show (n : ℝ) * (a - b) = n * a - n * b by ring,
      show (n : ℝ) * (a + b) = n * a + n * b by ring, cos_sub, cos_add]; ring
  simp only [ip, B, Phi, Nat.cast_zero, zero_mul, cos_zero]
  rw [sum_congr rfl h, sum_sub_distrib, sum_add_distrib, ← mul_sum]
  ring

/-! ## The dual vector -/

lemma ip_W_W {κ : Type*} (N : ℕ) (w : ℕ → ℝ) (θ : κ → ℝ) (J : Finset κ) :
    ip N w (Wfun θ J) (Wfun θ J) =
      ∑ j ∈ J, ∑ k ∈ J, (Phi N w (θ j - θ k) + Phi N w (θ j + θ k)) / 2 := by
  unfold Wfun
  rw [ip_sum_right N w _ J (fun j n => cos (n * θ j))]
  refine sum_congr rfl fun j _ => ?_
  rw [ip_comm, ip_sum_right N w _ J (fun k n => cos (n * θ k))]
  exact sum_congr rfl fun k _ => ip_cos_cos N w (θ j) (θ k)

/-- Diagonal bound (proofs.md step 6): `Φ ≥ 0 ⇒ ‖W‖² ≥ D = ½ Φ(0) #J`. -/
lemma diag_le {κ : Type*} (N : ℕ) (w : ℕ → ℝ) (hΦ : ∀ x, 0 ≤ Phi N w x) (θ : κ → ℝ)
    (J : Finset κ) : J.card * Phi N w 0 / 2 ≤ ip N w (Wfun θ J) (Wfun θ J) := by
  rw [ip_W_W]
  have h : ∀ j ∈ J, Phi N w 0 / 2 ≤
      ∑ k ∈ J, (Phi N w (θ j - θ k) + Phi N w (θ j + θ k)) / 2 := by
    intro j hj
    have hnn : ∀ k ∈ J, 0 ≤ (Phi N w (θ j - θ k) + Phi N w (θ j + θ k)) / 2 :=
      fun k _ => div_nonneg (add_nonneg (hΦ _) (hΦ _)) (by norm_num)
    have h1 := single_le_sum hnn hj
    simp only [sub_self] at h1
    have h2 := hΦ (θ j + θ j)
    linarith
  calc (J.card : ℝ) * Phi N w 0 / 2 = ∑ j ∈ J, Phi N w 0 / 2 := by
        rw [sum_const, nsmul_eq_mul]; ring
    _ ≤ _ := sum_le_sum h

/-! ## Li's formula -/

section Li

variable {κ : Type*} (N : ℕ) (w : ℕ → ℝ) (θ : κ → ℝ) (lam : ℕ → ℝ)

lemma lam_zero (hlam : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)) :
    lam 0 = 0 := by
  have h := hlam 0
  simp only [Nat.cast_zero, zero_mul, cos_zero, sub_self, mul_zero] at h
  exact h.unique hasSum_zero

/-- proofs.md step 3: `Σ_n φ(n/N) e^{-inθ} λ_{|n|} = Σ_{γ'} B(θ, θ_{γ'})` (finite `n`-sum,
interchange free). -/
lemma hasSum_B (hlam : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)) (a : ℝ) :
    HasSum (fun k => B N w a (θ k)) (ip N w (fun n => cos (n * a)) lam) := by
  have key : (fun k => B N w a (θ k)) = fun k =>
      w 0 * cos ((0 : ℕ) * a) * (2 * (1 - cos ((0 : ℕ) * θ k))) +
        2 * ∑ n ∈ Icc 1 N, w n * cos (n * a) * (2 * (1 - cos (n * θ k))) := by
    funext k; rw [← li_pair]; rfl
  rw [key]
  exact ((hlam 0).mul_left _).add ((hasSum_sum fun n _ => (hlam n).mul_left _).mul_left 2)

/-- `tail` is the sum over the zeros outside `J`. -/
lemma tail_eq_tsum_compl (hlam : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n))
    (J : Finset κ) (a : ℝ) :
    tail N w θ lam J a = ∑' k : ↑((J : Set κ)ᶜ), B N w a (θ k) := by
  have h := hasSum_B N w θ lam hlam a
  unfold tail
  rw [← h.tsum_eq, ← h.summable.sum_add_tsum_compl (s := J)]
  ring

end Li

/-- The exact identity `⟨R̃, W⟩ = τ' - 2‖W‖²` (proofs.md step 5 with every term kept; `τ'` here is
the exact residual `tauPrime`). -/
lemma ip_R_W {κ : Type*} (N : ℕ) (w : ℕ → ℝ) (θ : κ → ℝ) (lam t : ℕ → ℝ) (J : Finset κ) :
    ip N w (fun n => lam n - t n) (Wfun θ J) =
      tauPrime N w θ lam t J - 2 * ip N w (Wfun θ J) (Wfun θ J) := by
  have hj : ∀ j ∈ J, ip N w (fun n => lam n - t n) (fun n => cos (n * θ j)) =
      2 * J.card * Phi N w (θ j) + tail N w θ lam J (θ j) - trend N w t (θ j)
        - 2 * ∑ k ∈ J, (Phi N w (θ j - θ k) + Phi N w (θ j + θ k)) / 2 := by
    intro j _
    have hB : ∑ k ∈ J, B N w (θ j) (θ k) = ∑ k ∈ J, (2 * Phi N w (θ j) -
        2 * ((Phi N w (θ j - θ k) + Phi N w (θ j + θ k)) / 2)) :=
      sum_congr rfl fun k _ => by unfold B; ring
    rw [sum_sub_distrib, sum_const, nsmul_eq_mul, ← mul_sum] at hB
    rw [ip_comm, ip_sub_right]
    unfold tail trend
    rw [hB]
    ring
  rw [ip_W_W]
  unfold Wfun
  rw [ip_sum_right N w _ J (fun j n => cos (n * θ j)), sum_congr rfl hj]
  unfold tauPrime
  rw [sum_sub_distrib, sum_sub_distrib, sum_add_distrib, ← mul_sum, ← mul_sum]

/-- **Theorem 1 core** (proofs.md Theorem 1, steps 2–6). Hypotheses: `0 ≤ φ ≤ 1`, `Φ ≥ 0`
(Fourier positivity of `φ = ψ * ψ`), Li's formula with real phases (RH), `t(0) = 0`; `J` any
finite set of zeros, `τ` any upper bound of the residual `τ'` (proofs.md's `τ'` is one). Conclusion:
`Σ_{n≤N} R(n)² ≥ Φ(0)#J (1 - τ/(Φ(0)#J))² = (Φ(0)#J - τ)²/(Φ(0)#J)`. -/
theorem core_bound {κ : Type*} (N : ℕ) (w : ℕ → ℝ) (hw0 : ∀ n, 0 ≤ w n) (hw1 : ∀ n, w n ≤ 1)
    (hΦ : ∀ x, 0 ≤ Phi N w x) (θ : κ → ℝ) (lam t : ℕ → ℝ)
    (hlam : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)) (ht0 : t 0 = 0)
    (J : Finset κ) (τ : ℝ) (hτ : tauPrime N w θ lam t J ≤ τ) (hτ0 : 0 ≤ τ)
    (hτM : τ < J.card * Phi N w 0) :
    (J.card * Phi N w 0 - τ) ^ 2 / (J.card * Phi N w 0) ≤
      ∑ n ∈ Icc 1 N, (lam n - t n) ^ 2 := by
  set M : ℝ := J.card * Phi N w 0 with hM
  set x : ℝ := ip N w (Wfun θ J) (Wfun θ J) with hx
  set A : ℝ := ip N w (fun n => lam n - t n) (Wfun θ J) with hA
  have hMpos : 0 < M := by linarith
  have hxD : M / 2 ≤ x := diag_le N w hΦ θ J
  have hAeq : A = tauPrime N w θ lam t J - 2 * x := ip_R_W N w θ lam t J
  have hCS : A ^ 2 ≤ ip N w (fun n => lam n - t n) (fun n => lam n - t n) * x :=
    ip_sq_le N w hw0 _ _
  have hQsum : ip N w (fun n => lam n - t n) (fun n => lam n - t n) ≤
      2 * ∑ n ∈ Icc 1 N, (lam n - t n) ^ 2 := by
    have h0 : lam 0 - t 0 = 0 := by rw [lam_zero θ lam hlam, ht0, sub_self]
    simp only [ip, h0, mul_zero, zero_add]
    have : ∑ n ∈ Icc 1 N, w n * (lam n - t n) * (lam n - t n) ≤
        ∑ n ∈ Icc 1 N, (lam n - t n) ^ 2 :=
      sum_le_sum fun n _ => by
        nlinarith [mul_nonneg (sub_nonneg.2 (hw1 n)) (sq_nonneg (lam n - t n))]
    linarith
  set Q := ip N w (fun n => lam n - t n) (fun n => lam n - t n) with hQ
  have hxpos : 0 < x := by linarith
  have h1 : (2 * x - τ) ^ 2 ≤ A ^ 2 := by
    have ha : 0 ≤ 2 * x - τ := by linarith
    have hb : 2 * x - τ ≤ -A := by linarith
    nlinarith
  have h2 : (2 * x - τ) ^ 2 ≤ Q * x := h1.trans hCS
  have h3 : (M - τ) ^ 2 * x ≤ (2 * x - τ) ^ 2 * (M / 2) := by
    have ha : 0 ≤ x - M / 2 := by linarith
    have hb : 0 ≤ 2 * x * M - τ ^ 2 := by
      nlinarith [mul_nonneg ha hMpos.le, mul_nonneg (sub_nonneg.2 hτM.le) (add_nonneg hMpos.le hτ0)]
    nlinarith [mul_nonneg ha hb]
  have h4 : (M - τ) ^ 2 * x ≤ (Q * (M / 2)) * x := by
    nlinarith [mul_le_mul_of_nonneg_right h2 (by linarith : (0 : ℝ) ≤ M / 2)]
  have h5 : (M - τ) ^ 2 ≤ Q * (M / 2) := le_of_mul_le_mul_right h4 hxpos
  rw [div_le_iff₀ hMpos]
  nlinarith [mul_le_mul_of_nonneg_right hQsum hMpos.le]

end LiRate
