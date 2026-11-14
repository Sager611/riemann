/-
E1Lean.Arch — Theorem E1 (E1.md §1) from the named analytic inputs H1–H4 (sources: LEAN.md, Hypotheses).
-/
import Riemann.Deep.SantiAudit.Lean.Kappa
import Riemann.Deep.SantiAudit.Lean.Detector

namespace E1Lean

/-- The objects the analytic inputs speak about; abstract, their meaning is in these docstrings. -/
structure Inputs where
  /-- N(T, H): zeros β + iγ of ζ with T < γ ≤ T + H, counted with multiplicity. -/
  N : ℝ → ℝ → ℝ
  /-- S(T, H): simple zeros with β = 1/2 and T < γ ≤ T + H. -/
  S : ℝ → ℝ → ℝ
  /-- O(T, H): distinct sign changes of Hardy's Z(t) in (T, T + H]. -/
  O : ℝ → ℝ → ℝ
  /-- `Moment θ ν R P Q`: SL Thm 1.2's conclusion at H = T^θ, y = T^ν, σ₀ = 1/2 − R/L:
    ∫ w |Vψ(σ₀ + it)|² dt = c(P, Q, R, ν) ŵ(0) + O(H/L) for every admissible weight w. -/
  Moment : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop
  /-- `cdet P Q R ν` = c(P, Q, R, ν) = 1 + ν⁻¹ ∫₀¹∫₀¹ (w_R(v)P′(u) + ν w_R′(v)P(u))² du dv, w_R(v) = e^{Rv}Q(v). -/
  cdet : (ℝ → ℝ) → (ℝ → ℝ) → ℝ → ℝ → ℝ

variable (F : Fns) (I : Inputs)

/-- SL Thm 1.2's mollifier data: R > 0; P, Q real polynomials; P(0) = 0, P(1) = 1, Q(0) = 1. -/
def Adm (R : ℝ) (P Q : ℝ → ℝ) : Prop := 0 < R ∧ IsPoly P ∧ IsPoly Q ∧ P 0 = 0 ∧ P 1 = 1 ∧ Q 0 = 1

/-- SL Thm 1.2's conclusion at (θ, ν) for all admissible (R, P, Q). -/
def MomentRange (θ ν : ℝ) : Prop := ∀ R P Q, Adm R P Q → I.Moment θ ν R P Q

/-- (H1) Theorem E1-M (E1.md §1–§2, the new input): the moment holds for every 0 < ν < (3θ − 1)/4. -/
def H1 : Prop := ∀ θ ν : ℝ, 1 / 2 < θ → θ < 1 → 0 < ν → ν < (3 * θ - 1) / 4 → MomentRange I θ ν

/-- (H2) SL Thm 1.3 (Zenodo 23134787 v0.03, PDF text l. 125–129), at every (θ, ν) where Thm 1.2's conclusion
  holds: Q(x) + Q(1 − x) a nonzero constant ⟹ liminf O/N ≥ κ(P, Q, R, ν) = 1 − R⁻¹ log c(P, Q, R, ν). -/
def H2 : Prop := ∀ (θ ν R : ℝ) (P Q : ℝ → ℝ), 1 / 2 < θ → θ < 1 → 0 < ν → ν < 1 / 2 → Adm R P Q →
  (∃ C : ℝ, C ≠ 0 ∧ ∀ x, Q x + Q (1 - x) = C) → MomentRange I θ ν →
  LiminfGe (ratio F I.O I.N θ) (kappa F (I.cdet P Q R ν) R)

/-- (H3) SL Thm 5.1 (PDF text l. 844–845; parity product [SL ref 4, Prop 2.5] + Wang, arXiv 2609.07918):
  liminf O/N ≥ k ≥ 0 ⟹ liminf S/N ≥ h(θ; k). -/
def H3 : Prop := ∀ θ k : ℝ, 1 / 2 < θ → θ < 1 → 0 ≤ k → k < 1 →
  LiminfGe (ratio F I.O I.N θ) k → LiminfGe (ratio F I.S I.N θ) (h F θ k)

/-- (H4) Detector certificate (e1_constants.py/.log, Arb 4000 bits): c(P₁, Q₁, 177/20, 3/25) ∈ [c_lo, c_hi]. -/
def H4 : Prop := (3302286923204427565638 / 1000000000000000000 : ℝ) ≤ I.cdet P₁ Q₁ (177 / 20) (3 / 25) ∧
  I.cdet P₁ Q₁ (177 / 20) (3 / 25) ≤ 3302286923204427565640 / 1000000000000000000

/-- κ₁ = 1 − R⁻¹ log c(P₁, Q₁, R, ν) at R = 177/20, ν = 3/25 (E1.md §4). -/
noncomputable def κ₁ : ℝ := kappa F (I.cdet P₁ Q₁ (177 / 20) (3 / 25)) (177 / 20)

/-- **Theorem E1** (E1.md §1): for every θ ∈ (1/2, 1), liminf S(T, T^θ)/N(T, T^θ) ≥ h(θ; κ₁) ≥ h(1/2; κ₁) > 0.0055. -/
theorem theoremE1 (h1 : H1 I) (h2 : H2 F I) (h3 : H3 F I) (h4 : H4 I) (θ : ℝ) (hθ1 : 1 / 2 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₁ F I)) ∧ h F (1 / 2) (κ₁ F I) ≤ h F θ (κ₁ F I) ∧
      (0.0055 : ℝ) < h F (1 / 2) (κ₁ F I) := by
  obtain ⟨hk, hk1⟩ := kappa_bounds F h4.1 h4.2
  have hν : (3 / 25 : ℝ) < (3 * θ - 1) / 4 := (nu_range_iff θ).2 (by linarith)
  have hM : MomentRange I θ (3 / 25) := h1 θ (3 / 25) hθ1 hθ2 (by norm_num) hν
  obtain ⟨hP, hQ, hP0, hP1, hQ0, -, hsym⟩ := detector_admissible
  have hadm : Adm (177 / 20) P₁ Q₁ := ⟨by norm_num, hP, hQ, hP0, hP1, hQ0⟩
  have hO : LiminfGe (ratio F I.O I.N θ) (κ₁ F I) :=
    h2 θ (3 / 25) (177 / 20) P₁ Q₁ hθ1 hθ2 (by norm_num) (by norm_num) hadm ⟨1, one_ne_zero, hsym⟩ hM
  have hS := h3 θ (κ₁ F I) hθ1 hθ2 (le_trans (by norm_num) hk) hk1 hO
  exact ⟨hS, h_mono_theta F (by norm_num) hθ1.le hθ2.le hk1, h_half_gt F hk hk1⟩

/-- E1.md §4: from H4 alone, κ₁ − k*(1/2) > 0 and κ₁ > k*(θ) for every θ ∈ [1/2, 1]. -/
theorem kappa1_gt_kstar (h4 : H4 I) {θ : ℝ} (hθ1 : 1 / 2 ≤ θ) (hθ2 : θ ≤ 1) :
    0 < κ₁ F I - kstar F (1 / 2) ∧ kstar F θ < κ₁ F I := by
  obtain ⟨hk, hk1⟩ := kappa_bounds F h4.1 h4.2
  have hpos : 0 < h F θ (κ₁ F I) :=
    lt_of_lt_of_le (lt_trans (by norm_num) (h_half_gt F hk hk1)) (h_mono_theta F (by norm_num) hθ1 hθ2 hk1)
  exact ⟨sub_pos.2 (kstar_half_lt F hk hk1), (h_pos_iff F (by linarith) hθ2 hk1).1 hpos⟩

/-- A model of the inputs: N = S = O = 1, Moment ≡ True, c ≡ c_mid (used only for non-vacuity). -/
noncomputable def trivialInputs : Inputs :=
  ⟨fun _ _ => 1, fun _ _ => 1, fun _ _ => 1, fun _ _ _ _ _ => True,
    fun _ _ _ _ => 3302286923204427565639 / 1000000000000000000⟩

/-- Non-vacuity: for every `F`, H1–H4 are jointly satisfiable (by `trivialInputs`). -/
theorem inputs_consistent : H1 trivialInputs ∧ H2 F trivialInputs ∧ H3 F trivialInputs ∧ H4 trivialInputs := by
  have hc0 : (0 : ℝ) < 3302286923204427565639 / 1000000000000000000 := by norm_num
  have hlog : 0 < F.log (3302286923204427565639 / 1000000000000000000) := by
    by_contra hneg
    rw [not_lt] at hneg
    have h1 := (F.log_le_iff_le_exp _ 0 hc0).1 hneg
    rw [F.exp_zero] at h1
    norm_num at h1
  refine ⟨fun _ _ _ _ _ _ _ _ _ _ => trivial, ?_, ?_, ?_⟩
  · intro θ ν R P Q _ _ _ _ hadm _ _ ε hε
    refine ⟨0, fun T _ => ?_⟩
    have hk : kappa F (3302286923204427565639 / 1000000000000000000) R ≤ 1 := by
      unfold kappa
      have := mul_pos (inv_pos.2 hadm.1) hlog
      linarith
    simp only [ratio, trivialInputs]
    norm_num
    linarith
  · intro θ k _ _ _ hk _ ε hε
    refine ⟨0, fun T _ => ?_⟩
    have hh : h F θ k ≤ 1 := by
      unfold h hc
      have := F.sqrt_nonneg ((1 - k) * (9 - k - 8 * cW F θ))
      linarith
    simp only [ratio, trivialInputs]
    norm_num
    linarith
  · simp only [H4, trivialInputs]
    norm_num

end E1Lean
