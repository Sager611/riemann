/-
E1PLean.Arch — Theorem E1′ (Riemann/Deep/B3Cert/REPORT.md, Theorem; FULLPROOF.md §5 l. 105–118, §8) from the named
analytic inputs H1′–H3′ and the Arb certificate H4′ (sources: Riemann/Deep/B3Cert/LEAN.md, Hypotheses).
Reuses E1Lean's definitions unchanged (Fns, cot, cW, hc, h, kstar, kappa, LiminfGe, ratio, Inputs, Adm, MomentRange, P₁, Q₁).
-/
import Riemann.Deep.B3Cert.Lean.Num

namespace E1PLean
open E1Lean

variable (F : Fns) (I : Inputs)

/-- (H1′) Theorem E1-M′ (REPORT 'Moment form'; FULLPROOF §§1–5, Lemma B, Lemma R, §4 branches): SL Thm 1.2's
  conclusion for every θ ∈ (1/3, 1), 0 < ν < (3θ − 1)/4 and every admissible (R, P, Q). -/
def H1P : Prop := ∀ θ ν : ℝ, 1 / 3 < θ → θ < 1 → 0 < ν → ν < (3 * θ - 1) / 4 → MomentRange I θ ν

/-- (H2′) SL Thm 1.3 (Zenodo 23134787 v0.03, pdftotext l. 125–129) at every (θ, ν) where Thm 1.2's conclusion holds,
  θ ∈ (0, 1): θ-free per FULLPROOF §5 l. 111–112, 116 (SL TeX l. 823–928: Lemma 4.2, Thm 1.3). -/
def H2P : Prop := ∀ (θ ν R : ℝ) (P Q : ℝ → ℝ), 0 < θ → θ < 1 → 0 < ν → ν < 1 / 2 → Adm R P Q →
  (∃ C : ℝ, C ≠ 0 ∧ ∀ x, Q x + Q (1 - x) = C) → MomentRange I θ ν →
  LiminfGe (ratio F I.O I.N θ) (kappa F (I.cdet P Q R ν) R)

/-- (H3′) SL Thm 5.1 (pdftotext l. 844–845; TeX l. 932–963), θ ∈ (0, 1): θ-free per FULLPROOF §5 l. 112–114, 116
  (Wang 2609.07918 Thm local-pair, its TeX l. 231–232, 'Let 0 < λ < θ < 1 be fixed'). -/
def H3P : Prop := ∀ θ k : ℝ, 0 < θ → θ < 1 → 0 ≤ k → k < 1 →
  LiminfGe (ratio F I.O I.N θ) k → LiminfGe (ratio F I.S I.N θ) (h F θ k)

/-- (H4′) Detector certificate (code/b3c_arb.py, logs/b3c_arb.log, Arb 4000 bits): c(P₁, Q₁, 44/5, 121/1000) ∈
  [3138.997965019943478277, 3138.997965019943478278] ⊇ [3138.99796501994347827760614809 ± 1.36·10⁻²⁷]. -/
def H4P : Prop := (3138997965019943478277 / 1000000000000000000 : ℝ) ≤ I.cdet P₁ Q₁ (44 / 5) (121 / 1000) ∧
  I.cdet P₁ Q₁ (44 / 5) (121 / 1000) ≤ 3138997965019943478278 / 1000000000000000000

/-- κ₂ = 1 − R⁻¹ log c(P₁, Q₁, R, ν) at R = 44/5, ν = 121/1000 (REPORT, Theorem). -/
noncomputable def κ₂ : ℝ := kappa F (I.cdet P₁ Q₁ (44 / 5) (121 / 1000)) (44 / 5)

/-- The moment edge: ν = 121/1000 < (3θ − 1)/4 ⟺ θ > 371/750. -/
theorem nu_edge_iff (θ : ℝ) : (121 / 1000 : ℝ) < (3 * θ - 1) / 4 ↔ 371 / 750 < θ := by
  constructor <;> intro hθ <;> linarith

/-- Equality at the edge: (3·371/750 − 1)/4 = 121/1000. -/
theorem nu_edge_eq : (3 * (371 / 750 : ℝ) - 1) / 4 = 121 / 1000 := by norm_num

/-- FULLPROOF §2 bookkeeping at ν = 121/1000: the relative error T^{1/2 + 2ν − 3θ/2 + ε} is power-small iff θ > 371/750. -/
theorem exponent_edge (θ : ℝ) : 1 / 2 + 2 * (121 / 1000 : ℝ) - 3 * θ / 2 < 0 ↔ 371 / 750 < θ :=
  (exponent_iff θ (121 / 1000)).trans (nu_edge_iff θ)

/-- **Theorem E1′** (REPORT): for every θ ∈ (371/750, 1),
  liminf S(T, T^θ)/N(T, T^θ) ≥ h(θ; κ₂) ≥ h(371/750; κ₂) > 3.35·10⁻⁵. -/
theorem theoremE1P (h1 : H1P I) (h2 : H2P F I) (h3 : H3P F I) (h4 : H4P I) (θ : ℝ) (hθ1 : 371 / 750 < θ)
    (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₂ F I)) ∧ h F (371 / 750) (κ₂ F I) ≤ h F θ (κ₂ F I) ∧
      (0.0000335 : ℝ) < h F (371 / 750) (κ₂ F I) := by
  obtain ⟨hk, hk1⟩ := kappa2_bounds F h4.1 h4.2
  have hν : (121 / 1000 : ℝ) < (3 * θ - 1) / 4 := (nu_edge_iff θ).2 hθ1
  have hM : MomentRange I θ (121 / 1000) := h1 θ (121 / 1000) (by linarith) hθ2 (by norm_num) hν
  obtain ⟨hP, hQ, hP0, hP1, hQ0, -, hsym⟩ := detector_admissible
  have hadm : Adm (44 / 5) P₁ Q₁ := ⟨by norm_num, hP, hQ, hP0, hP1, hQ0⟩
  have hO : LiminfGe (ratio F I.O I.N θ) (κ₂ F I) :=
    h2 θ (121 / 1000) (44 / 5) P₁ Q₁ (by linarith) hθ2 (by norm_num) (by norm_num) hadm ⟨1, one_ne_zero, hsym⟩ hM
  have hS := h3 θ (κ₂ F I) (by linarith) hθ2 (le_trans (by norm_num) hk) hk1 hO
  exact ⟨hS, h_mono_theta F (by norm_num) hθ1.le hθ2.le hk1, lt_trans (by norm_num) (h_edge_gt F hk hk1)⟩

/-- From H4′ alone: κ₂ − k*(371/750) > 0 and κ₂ > k*(θ) for every θ ∈ [371/750, 1]. -/
theorem kappa2_gt_kstar (h4 : H4P I) {θ : ℝ} (hθ1 : 371 / 750 ≤ θ) (hθ2 : θ ≤ 1) :
    0 < κ₂ F I - kstar F (371 / 750) ∧ kstar F θ < κ₂ F I := by
  obtain ⟨hk, hk1⟩ := kappa2_bounds F h4.1 h4.2
  have h0 : 0 < h F (371 / 750) (κ₂ F I) := lt_trans (by norm_num) (h_edge_gt F hk hk1)
  have hpos : 0 < h F θ (κ₂ F I) := lt_of_lt_of_le h0 (h_mono_theta F (by norm_num) hθ1 hθ2 hk1)
  exact ⟨sub_pos.2 ((h_pos_iff F (by norm_num) (by norm_num) hk1).1 h0), (h_pos_iff F (by linarith) hθ2 hk1).1 hpos⟩

/-- **E1′ constants** (FULLPROOF §8), from H4′: c(371/750), κ₂, k*(371/750), h(371/750; κ₂), κ₂ − k*(371/750) enclosed
  (Arb: −0.18577184013932, 0.08503876039177, 0.08499141434976, 3.354573606·10⁻⁵, 4.7346042·10⁻⁵). -/
theorem e1p_enclosures (h4 : H4P I) :
    (-1857718405 / 10000000000 : ℝ) ≤ cW F (371 / 750) ∧ cW F (371 / 750) ≤ (-18577184 / 100000000 : ℝ) ∧
    (850387602 / 10000000000 : ℝ) ≤ κ₂ F I ∧ κ₂ F I ≤ 850387611 / 10000000000 ∧
    (849914142 / 10000000000 : ℝ) ≤ kstar F (371 / 750) ∧ kstar F (371 / 750) ≤ 849914146 / 10000000000 ∧
    (0.000033545 : ℝ) < h F (371 / 750) (κ₂ F I) ∧ h F (371 / 750) (κ₂ F I) ≤ 0.000033547 ∧
    (0.000047345 : ℝ) < κ₂ F I - kstar F (371 / 750) ∧ κ₂ F I - kstar F (371 / 750) ≤ 0.000047347 := by
  obtain ⟨hk, hk1⟩ := kappa2_bounds F h4.1 h4.2
  have hku := kappa2_upper F h4.1
  obtain ⟨ks1, ks2⟩ := kstar_edge_bounds F
  have g1 : (0.000047345 : ℝ) < 850387602 / 10000000000 - 849914146 / 10000000000 := by norm_num
  have g2 : (850387611 / 10000000000 : ℝ) - 849914142 / 10000000000 ≤ 0.000047347 := by norm_num
  unfold κ₂
  exact ⟨cW_edge_ge F, cW_edge_le F, hk, hku, ks1, ks2, h_edge_gt F hk hk1, h_edge_le F hku,
    by linarith, by linarith⟩

/-- A model of the inputs: N = S = O = 1, Moment ≡ True, c ≡ c_mid (used only for non-vacuity). -/
noncomputable def trivialInputsP : Inputs :=
  ⟨fun _ _ => 1, fun _ _ => 1, fun _ _ => 1, fun _ _ _ _ _ => True,
    fun _ _ _ _ => 6277995930039886956555 / 2000000000000000000⟩

/-- Non-vacuity: for every F, H1′–H4′ are jointly satisfiable (by trivialInputsP). -/
theorem inputs_consistentP : H1P trivialInputsP ∧ H2P F trivialInputsP ∧ H3P F trivialInputsP ∧ H4P trivialInputsP := by
  have hc0 : (0 : ℝ) < 6277995930039886956555 / 2000000000000000000 := by norm_num
  have hlog : 0 < F.log (6277995930039886956555 / 2000000000000000000) := by
    by_contra hneg
    rw [not_lt] at hneg
    have h1 := (F.log_le_iff_le_exp _ 0 hc0).1 hneg
    rw [F.exp_zero] at h1
    norm_num at h1
  refine ⟨fun _ _ _ _ _ _ _ _ _ _ => trivial, ?_, ?_, ?_⟩
  · intro θ ν R P Q _ _ _ _ hadm _ _ ε hε
    refine ⟨0, fun T _ => ?_⟩
    have hk : kappa F (6277995930039886956555 / 2000000000000000000) R ≤ 1 := by
      unfold kappa
      have := mul_pos (inv_pos.2 hadm.1) hlog
      linarith
    simp only [ratio, trivialInputsP]
    norm_num
    linarith
  · intro θ k _ _ _ hk _ ε hε
    refine ⟨0, fun T _ => ?_⟩
    have hh : h F θ k ≤ 1 := by
      unfold h hc
      have := F.sqrt_nonneg ((1 - k) * (9 - k - 8 * cW F θ))
      linarith
    simp only [ratio, trivialInputsP]
    norm_num
    linarith
  · simp only [H4P, trivialInputsP]
    norm_num

end E1PLean
