import Riemann.Deep.GL2Ceil.Lean.InnerMin
import Riemann.Deep.GL2Ceil.Lean.Duality
import Riemann.Deep.GL2Ceil.Lean.Crude
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.Topology.Algebra.Polynomial

/-!
# GL2Ceil (gl2lean): the ceiling theorem from the certificate interface (DERIVATION §3-§5)

ξ₀ = 25/256 (ν < 25/128, AT21 Thm thm:nu l. 227-231 with θ = 7/64), κ_up = 0.07005633.
Named hypotheses (checked outside Lean):
* Cover cells: the cells tile [1/2, 28] × [0, 3] (covcheck.py).
* CertCells cells: each dual cell carries an admissible dual (σ, μ) with D(R, C_λ/ξ₀, ξ₀B_λ) ≥ e^{(1-κ_up)R}
  on the whole cell (ceil.py check_lam, Arb; replay.py); each crude cell [R_a, R_b] × [λ_a, λ_b] satisfies
  the number inequality e^{(1-κ_up)R_b} ≤ crude(R_a, (C_{λ_a}/ξ₀) G(R_a)) (§4.3 (i), Arb).
* TailSmallR (R ≤ 1/2, all λ) and the first half of TailLargeLam (λ ≥ 3, R ∈ [1/2, 19/32]):
  M(R,λ) ≥ e^{(1-κ_up)R} (§4.1-4.2, closed form at R = 0 and slope; not formalized).
* Second half of TailLargeLam: e^{(1-κ_up)R} ≤ crude(R, (C_3/ξ₀) G(R)) on [19/32, 28] (§4.3 (ii));
  TailLargeR: e^{(1-κ_up)R} ≤ crude(R, G(R)/ξ₀) for R ≥ 28 (§4.3 (iii)): one-variable number inequalities.
Proved in Lean: expansion (E), the sinh inner minimum, monotonicity of C_λ, B_λ, the ξ-monotonicity,
Lemma D, the crude bound J ≥ crude(R, αG(R)) with its monotonicity, and the assembly: κ(P,Q,R,ξ) ≤ 0.07005633 for every C¹ (in particular polynomial) P, Q of the
class, R > 0 and 0 < ξ ≤ 25/256.
-/

open Real MeasureTheory Set

noncomputable section

namespace GL2Ceil

/-- ξ₀ = ν₀/2 = 25/256. -/
def xi0 : ℝ := 25 / 256
/-- κ_up = 0.07005633. -/
def kup : ℝ := 7005633 / 100000000
/-- e^{(1-κ_up)R}. -/
def target (R : ℝ) : ℝ := Real.exp ((1 - kup) * R)

/-- M(R,λ) ≥ t at ξ: J(ξ,R,λ,Q) ≥ t for every Q in the real class. -/
def MGe (ξ R lam t : ℝ) : Prop := ∀ Q, QAdm Q → t ≤ J ξ R lam Q

/-- A cell of logs/ceil_25_256_0.07005633.json: [Rc - dR, Rc + dR] × [l0, l1]·2^{-20}. -/
structure Cell where
  Rc : ℚ
  dR : ℚ
  l0 : ℕ
  l1 : ℕ
  dual : Bool

def Cell.mem (c : Cell) (R lam : ℝ) : Prop :=
  ((c.Rc : ℝ) - c.dR ≤ R ∧ R ≤ (c.Rc : ℝ) + c.dR) ∧
    ((c.l0 : ℝ) / 2 ^ 20 ≤ lam ∧ lam ≤ (c.l1 : ℝ) / 2 ^ 20)

/-- §3: an admissible dual whose value beats the target on the whole cell. -/
def DualCellCert (c : Cell) : Prop :=
  ∃ σ σ' μ : ℝ → ℝ, DualAdm σ σ' μ ∧
    ∀ R lam, c.mem R lam → target R ≤ Dual R (Cl lam / xi0) (xi0 * Bl lam) σ σ' μ

/-- §4.3 (i): the crude bound closes the cell: e^{(1-κ_up)R_b} ≤ crude(R_a, (C_{λ_a}/ξ₀) G(R_a)). -/
def CrudeCellCert (c : Cell) : Prop :=
  target ((c.Rc : ℝ) + c.dR) ≤ crude ((c.Rc : ℝ) - c.dR) (Cl ((c.l0 : ℝ) / 2 ^ 20) / xi0 * G ((c.Rc : ℝ) - c.dR))

def CertCells (cells : List Cell) : Prop :=
  ∀ c ∈ cells, (c.dual = true → DualCellCert c) ∧ (c.dual = false → CrudeCellCert c)

def Cover (cells : List Cell) : Prop :=
  ∀ R lam : ℝ, 1 / 2 ≤ R → R ≤ 28 → 0 ≤ lam → lam ≤ 3 → ∃ c ∈ cells, c.mem R lam

def TailSmallR : Prop := ∀ R lam : ℝ, 0 < R → R ≤ 1 / 2 → 0 ≤ lam → MGe xi0 R lam (target R)
def TailLargeLam : Prop :=
  (∀ R lam : ℝ, 1 / 2 ≤ R → R ≤ 19 / 32 → 3 ≤ lam → MGe xi0 R lam (target R)) ∧
    ∀ R : ℝ, 19 / 32 ≤ R → R ≤ 28 → target R ≤ crude R (Cl 3 / xi0 * G R)
def TailLargeR : Prop := ∀ R : ℝ, 28 ≤ R → target R ≤ crude R (1 / xi0 * G R)

lemma target_mono {R R' : ℝ} (h : R ≤ R') : target R ≤ target R' :=
  Real.exp_le_exp.2 (mul_le_mul_of_nonneg_left h (by norm_num [kup]))

/-- A crude inequality at (R, K) with K ≤ (C_λ/ξ₀) G(R) gives M(R,λ) ≥ target. -/
lemma mge_of_crude {R lam K : ℝ} (hlam : 0 ≤ lam) (hK0 : 0 ≤ K) (hK : K ≤ Cl lam / xi0 * G R)
    (h : target R ≤ crude R K) : MGe xi0 R lam (target R) := by
  intro Q hQ
  have hα : 0 ≤ Cl lam / xi0 := div_nonneg (zero_le_one.trans (Cl_ge_one hlam)) (by norm_num [xi0])
  have hγ : 0 ≤ xi0 * Bl lam := mul_nonneg (by norm_num [xi0]) (Bl_pos lam).le
  exact h.trans ((crude_mono le_rfl hK0 hK).trans (crude_le_Jag (R := R) hα hγ hQ))

/-- Cover + cells + tails: M(R,λ) ≥ e^{(1-κ_up)R} for all R > 0, λ ≥ 0. -/
theorem M_ge_all (cells : List Cell) (hcov : Cover cells) (hcert : CertCells cells)
    (hS : TailSmallR) (hL : TailLargeLam) (hR : TailLargeR) {R lam : ℝ} (hR0 : 0 < R)
    (hlam : 0 ≤ lam) : MGe xi0 R lam (target R) := by
  rcases le_or_gt R (1 / 2) with h1 | h1
  · exact hS R lam hR0 h1 hlam
  rcases le_or_gt R 28 with h2 | h2
  · rcases le_or_gt lam 3 with h3 | h3
    · obtain ⟨c, hc, hmem⟩ := hcov R lam h1.le h2 hlam h3
      obtain ⟨hd, hcr⟩ := hcert c hc
      cases hb : c.dual
      · have hct := hcr hb
        obtain ⟨⟨hRa, hRb⟩, ⟨hla, _⟩⟩ := hmem
        have hla0 : (0:ℝ) ≤ (c.l0 : ℝ) / 2 ^ 20 := by positivity
        have hC : Cl ((c.l0 : ℝ) / 2 ^ 20) / xi0 ≤ Cl lam / xi0 := by
          have := Cl_mono hla0 hla
          unfold xi0
          linarith
        have hC0 : 0 ≤ Cl ((c.l0 : ℝ) / 2 ^ 20) / xi0 :=
          div_nonneg (zero_le_one.trans (Cl_ge_one hla0)) (by norm_num [xi0])
        refine mge_of_crude hlam (mul_nonneg hC0 (G_nonneg _))
          (mul_le_mul hC (G_mono hRa) (G_nonneg _) (hC0.trans hC)) ?_
        exact ((target_mono hRb).trans hct).trans
          (crude_mono hRa (mul_nonneg hC0 (G_nonneg _)) le_rfl)
      · obtain ⟨σ, σ', μ, hadm, hD⟩ := hd hb
        intro Q hQ
        have hα : 0 < Cl lam / xi0 :=
          div_pos (lt_of_lt_of_le one_pos (Cl_ge_one hlam)) (by norm_num [xi0])
        have hγ : 0 < xi0 * Bl lam := mul_pos (by norm_num [xi0]) (Bl_pos lam)
        exact (hD R lam hmem).trans (weak_duality hα hγ hadm hQ)
    · rcases le_or_gt R (19 / 32) with h4 | h4
      · exact hL.1 R lam h1.le h4 h3.le
      · have hC3 : Cl 3 / xi0 ≤ Cl lam / xi0 := by
          have := Cl_mono (by norm_num : (0:ℝ) ≤ 3) h3.le
          unfold xi0
          linarith
        have hC0 : 0 ≤ Cl 3 / xi0 :=
          div_nonneg (zero_le_one.trans (Cl_ge_one (by norm_num))) (by norm_num [xi0])
        exact mge_of_crude hlam (mul_nonneg hC0 (G_nonneg _))
          (mul_le_mul_of_nonneg_right hC3 (G_nonneg _)) (hL.2 R h4.le h2)
  · have hC1 : 1 / xi0 ≤ Cl lam / xi0 := by
      have := Cl_ge_one hlam
      unfold xi0
      linarith
    exact mge_of_crude hlam (mul_nonneg (by norm_num [xi0]) (G_nonneg _))
      (mul_le_mul_of_nonneg_right hC1 (G_nonneg _)) (hR R h2.le)

/-- Main theorem (C¹ form): κ(P,Q,R,ξ) ≤ 0.07005633. -/
theorem ceiling (cells : List Cell) (hcov : Cover cells) (hcert : CertCells cells)
    (hS : TailSmallR) (hL : TailLargeLam) (hR : TailLargeR)
    {P Q : ℝ → ℝ} (hP : PAdm P) (hQ : QAdm Q) {R ξ : ℝ} (hR0 : 0 < R) (hξ : 0 < ξ)
    (hξ0 : ξ ≤ 25 / 256) : kappa P Q R ξ ≤ 0.07005633 := by
  have hk0 : 0 ≤ xi0 * Real.sqrt (W1 R Q / W0 R Q) :=
    mul_nonneg (by norm_num [xi0]) (Real.sqrt_nonneg _)
  have hM := M_ge_all cells hcov hcert hS hL hR hR0 hk0 Q hQ
  have hc := c_ge_J hP hQ.smooth hQ.zero R hξ (show ξ ≤ xi0 by unfold xi0; exact hξ0)
  have hle : target R ≤ c P Q R ξ := hM.trans hc
  have hlog : (1 - kup) * R ≤ Real.log (c P Q R ξ) := by
    have := Real.log_le_log (Real.exp_pos ((1 - kup) * R)) hle
    rwa [Real.log_exp] at this
  unfold kup at hlog
  have h2 : 1 - 7005633 / 100000000 ≤ 1 / R * Real.log (c P Q R ξ) := by
    rw [div_mul_eq_mul_div, one_mul, le_div_iff₀ hR0]
    exact hlog
  unfold kappa
  norm_num at h2 ⊢
  linarith

/-- Polynomials are C¹. -/
lemma poly_C1 (p : Polynomial ℝ) : ContDiff ℝ 1 (fun x => p.eval x) := by
  rw [contDiff_one_iff_deriv]
  refine ⟨p.differentiable, ?_⟩
  have : deriv (fun x => p.eval x) = fun x => p.derivative.eval x :=
    funext fun x => (p.hasDerivAt x).deriv
  rw [this]
  exact p.derivative.continuous

/-- Main theorem, literal polynomial form (B l. 161, l. 770-795; AT21 l. 192-218): real P, Q ∈ ℝ[x]
with P(0) = 0, P(1) = 1, Q(0) = 1, Q(y) + Q(1-y) constant. -/
theorem ceiling_poly (cells : List Cell) (hcov : Cover cells) (hcert : CertCells cells)
    (hS : TailSmallR) (hL : TailLargeLam) (hR : TailLargeR)
    (P Q : Polynomial ℝ) (hP0 : P.eval 0 = 0) (hP1 : P.eval 1 = 1) (hQ0 : Q.eval 0 = 1)
    (hQs : ∃ β : ℝ, ∀ y : ℝ, Q.eval y + Q.eval (1 - y) = β) {R ξ : ℝ} (hR0 : 0 < R) (hξ : 0 < ξ)
    (hξ0 : ξ ≤ 25 / 256) :
    kappa (fun x => P.eval x) (fun x => Q.eval x) R ξ ≤ 0.07005633 := by
  obtain ⟨β, hβ⟩ := hQs
  exact ceiling cells hcov hcert hS hL hR ⟨poly_C1 P, hP0, hP1⟩
    ⟨poly_C1 Q, hQ0, β, fun y _ => hβ y⟩ hR0 hξ hξ0

end GL2Ceil





