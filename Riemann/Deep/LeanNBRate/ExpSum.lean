import Riemann.Deep.LeanNBRate.LitInputs
import Mathlib.Analysis.SpecialFunctions.Exponential
import Mathlib.Analysis.Complex.Exponential
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds

/-!
# Exponential sums with jets and the Turán–Nazarov lemma for them

`r(u) = Σ_j Σ_{k < m_j} a_{jk} (iu)^k e^{iδ_j u}` with real nodes `δ_j`. The Turán–Nazarov input
(`Nazarov1993_Turan`, pure exponentials) extends to jets by the forward-difference limit
`e^{iδu}((e^{iεu} − 1)/ε)^k → (iu)^k e^{iδu}` (ours; [F] §0 "continuity argument").
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- A finite exponential sum with jets and real nodes. -/
structure ExpSum where
  J : ℕ
  node : Fin J → ℝ
  mult : Fin J → ℕ
  coef : Fin J → ℕ → ℂ

namespace ExpSum

variable (r : ExpSum)

/-- `r(u) = Σ_j Σ_{k < m_j} a_{jk} (iu)^k e^{iδ_j u}`. -/
def eval (u : ℝ) : ℂ :=
  ∑ j, ∑ k ∈ Finset.range (r.mult j), r.coef j k * (I * u) ^ k * cexp (I * r.node j * u)

/-- Total multiplicity `M = Σ_j m_j`. -/
def M : ℕ := ∑ j, r.mult j

/-- Forward-difference approximation by pure exponentials. -/
def approx (ε : ℝ) (u : ℝ) : ℂ :=
  ∑ j, ∑ k ∈ Finset.range (r.mult j),
    r.coef j k * ((cexp (I * ε * u) - 1) / ε) ^ k * cexp (I * r.node j * u)

/-- Coefficients of `approx` in the pure exponential basis `e^{i(δ_j + lε)u}`, `l < m_j`. -/
def acoef (ε : ℝ) (j : Fin r.J) (l : ℕ) : ℂ :=
  ∑ k ∈ Finset.range (r.mult j),
    if l ≤ k then r.coef j k * (ε : ℂ)⁻¹ ^ k * (Nat.choose k l : ℂ) * (-1) ^ (k - l) else 0

theorem approx_eq (ε : ℝ) (u : ℝ) :
    r.approx ε u = ∑ j, ∑ l ∈ Finset.range (r.mult j),
      r.acoef ε j l * cexp (I * (r.node j + l * ε) * u) := by
  unfold approx acoef
  refine Finset.sum_congr rfl fun j _ => ?_
  have hpow : ∀ k : ℕ, ((cexp (I * ε * u) - 1) / ε) ^ k =
      ∑ l ∈ Finset.range (k + 1), (ε : ℂ)⁻¹ ^ k * (Nat.choose k l : ℂ) * (-1) ^ (k - l) *
        cexp (I * (l * ε) * u) := by
    intro k
    rw [div_eq_mul_inv, mul_pow, sub_eq_add_neg, add_pow, Finset.sum_mul]
    refine Finset.sum_congr rfl fun l _ => ?_
    rw [← Complex.exp_nat_mul]
    ring_nf
  simp_rw [hpow, Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm' (t' := Finset.range (r.mult j))
    (s' := fun l => (Finset.range (r.mult j)).filter (fun k => l ≤ k))]
  · refine Finset.sum_congr rfl fun l _ => ?_
    rw [Finset.sum_filter]
    refine Finset.sum_congr rfl fun k _ => ?_
    split_ifs with h
    · have : cexp (I * (l * ε) * u) * cexp (I * r.node j * u) =
          cexp (I * (r.node j + l * ε) * u) := by rw [← Complex.exp_add]; ring_nf
      rw [← this]; ring
    · simp
  · intro k l
    simp only [Finset.mem_range, Finset.mem_filter]
    constructor <;> intro h <;> omega

end ExpSum

/-- Turán–Nazarov for pure exponential sums with real frequencies over any finite index type. -/
theorem turan_real (hT : Nazarov1993_Turan) : ∃ c : ℝ, 0 < c ∧
    ∀ {ι : Type} [Fintype ι] (coef : ι → ℂ) (θ : ι → ℝ) (a b : ℝ) (Ω : Set ℝ) (K : ℝ),
      1 ≤ Fintype.card ι → a < b → MeasurableSet Ω → Ω ⊆ Icc a b → 0 < volume Ω →
      (∀ t ∈ Ω, ‖∑ i, coef i * cexp (I * θ i * t)‖ ≤ K) →
      ∀ t ∈ Icc a b, ‖∑ i, coef i * cexp (I * θ i * t)‖ ≤
        (c * (b - a) / (volume Ω).toReal) ^ (Fintype.card ι - 1) * K := by
  obtain ⟨c, hc, hT⟩ := hT
  refine ⟨c, hc, fun {ι} _ coef θ a b Ω K hcard hab hΩ hΩs hΩpos hK t ht => ?_⟩
  set n := Fintype.card ι with hn
  have hn' : n - 1 + 1 = n := Nat.sub_add_cancel hcard
  let e : Fin (n - 1 + 1) ≃ ι := (finCongr hn').trans (Fintype.equivFin ι).symm
  have hsum : ∀ u : ℝ, ∑ i, coef i * cexp (I * θ i * u) =
      ∑ k : Fin (n - 1 + 1), coef (e k) * cexp ((I * θ (e k)) * u) := by
    intro u
    rw [← e.sum_comp]
  have hre : (⨆ k : Fin (n - 1 + 1), |(I * (θ (e k) : ℂ)).re|) = 0 := by
    simp
  have := hT (n - 1) (fun k => coef (e k)) (fun k => I * θ (e k)) a b Ω K hab hΩ hΩs hΩpos
    (fun u hu => by rw [← hsum u]; exact hK u hu) t ht
  rw [hre, mul_zero, Real.exp_zero, one_mul, ← hsum t] at this
  exact this

namespace ExpSum

variable (r : ExpSum)

/-- `approx` written as a pure exponential sum over the index type `Σ j, Fin (m_j)`. -/
theorem approx_eq_sigma (ε : ℝ) (u : ℝ) :
    r.approx ε u = ∑ x : (Σ j : Fin r.J, Fin (r.mult j)),
      r.acoef ε x.1 x.2 * cexp (I * ((r.node x.1 + (x.2 : ℕ) * ε : ℝ) : ℂ) * u) := by
  rw [approx_eq, Fintype.sum_sigma]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.sum_range (fun l => r.acoef ε j l * cexp (I * (r.node j + l * ε) * u))]
  refine Finset.sum_congr rfl fun l _ => ?_
  push_cast; ring_nf

theorem card_sigma : Fintype.card (Σ j : Fin r.J, Fin (r.mult j)) = r.M := by
  simp [M, Fintype.card_sigma]

/-- `|X^k − Y^k| ≤ k R^{k−1} |X − Y|` for `|X|, |Y| ≤ R`. -/
theorem norm_pow_sub_pow_le {X Y : ℂ} {R : ℝ} (hX : ‖X‖ ≤ R) (hY : ‖Y‖ ≤ R) (k : ℕ) :
    ‖X ^ k - Y ^ k‖ ≤ k * R ^ (k - 1) * ‖X - Y‖ := by
  have hR : 0 ≤ R := (norm_nonneg _).trans hX
  induction k with
  | zero => simp
  | succ k ih =>
    have h : X ^ (k + 1) - Y ^ (k + 1) = X * (X ^ k - Y ^ k) + (X - Y) * Y ^ k := by ring
    rw [h]
    calc ‖X * (X ^ k - Y ^ k) + (X - Y) * Y ^ k‖
        ≤ ‖X‖ * ‖X ^ k - Y ^ k‖ + ‖X - Y‖ * ‖Y‖ ^ k := by
          refine (norm_add_le _ _).trans ?_
          rw [norm_mul, norm_mul, norm_pow]
      _ ≤ R * (k * R ^ (k - 1) * ‖X - Y‖) + ‖X - Y‖ * R ^ k := by
          gcongr
      _ = (k + 1 : ℕ) * R ^ (k + 1 - 1) * ‖X - Y‖ := by
          rcases k with _ | k
          · simp
          · simp only [Nat.add_sub_cancel]; push_cast; ring

/-- Uniform error of the forward-difference approximation on `|u| ≤ R`, `0 < ε ≤ 1/(R+1)`. -/
theorem norm_approx_sub_eval_le {R ε u : ℝ} (hR : 0 ≤ R) (hε : 0 < ε) (hεR : ε * (R + 1) ≤ 1)
    (hu : |u| ≤ R) :
    ‖r.approx ε u - r.eval u‖ ≤
      ε * ∑ j, ∑ k ∈ Finset.range (r.mult j), ‖r.coef j k‖ * (k * (2 * R + 1) ^ (k + 1)) := by
  unfold approx eval
  rw [← Finset.sum_sub_distrib, Finset.mul_sum]
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun j _ => ?_)
  rw [← Finset.sum_sub_distrib, Finset.mul_sum]
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun k _ => ?_)
  set X : ℂ := (cexp (I * ε * u) - 1) / ε
  set Y : ℂ := I * u
  have hεu : ‖(I * ε * u : ℂ)‖ ≤ 1 := by
    rw [norm_mul, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Complex.norm_real,
      Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos hε]
    nlinarith [abs_nonneg u]
  have hY : ‖Y‖ ≤ 2 * R + 1 := by
    simp only [Y, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Real.norm_eq_abs]
    linarith
  have hX : ‖X‖ ≤ 2 * R + 1 := by
    have h1 : ‖cexp (I * ε * u) - 1‖ ≤ ‖(ε * u : ℝ)‖ := by
      have := Real.norm_exp_I_mul_ofReal_sub_one_le (x := ε * u)
      simpa [mul_assoc] using this
    simp only [X, norm_div, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hε]
    rw [div_le_iff₀ hε]
    calc ‖cexp (I * ε * u) - 1‖ ≤ |ε * u| := by simpa using h1
      _ = ε * |u| := by rw [abs_mul, abs_of_pos hε]
      _ ≤ (2 * R + 1) * ε := by nlinarith [abs_nonneg u]
  have hXY : ‖X - Y‖ ≤ ε * u ^ 2 := by
    have h2 := Complex.norm_exp_sub_one_sub_id_le hεu
    have hε' : (ε : ℂ) ≠ 0 := by exact_mod_cast hε.ne'
    have : X - Y = (cexp (I * ε * u) - 1 - I * ε * u) / ε := by
      simp only [X, Y]; field_simp
    rw [this, norm_div, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hε, div_le_iff₀ hε]
    calc ‖cexp (I * ε * u) - 1 - I * ε * u‖ ≤ ‖(I * ε * u : ℂ)‖ ^ 2 := h2
      _ = (ε * |u|) ^ 2 := by
          rw [norm_mul, norm_mul, Complex.norm_I, one_mul, Complex.norm_real, Complex.norm_real,
            Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos hε]
      _ = ε * u ^ 2 * ε := by rw [mul_pow, sq_abs]; ring
  have hexp : ‖cexp (I * r.node j * u)‖ = 1 := by
    rw [show (I * r.node j * u : ℂ) = ((r.node j * u : ℝ) : ℂ) * I by push_cast; ring]
    exact Complex.norm_exp_ofReal_mul_I _
  have hdiff : r.coef j k * X ^ k * cexp (I * r.node j * u) -
      r.coef j k * Y ^ k * cexp (I * r.node j * u) =
      r.coef j k * (X ^ k - Y ^ k) * cexp (I * r.node j * u) := by ring
  rw [hdiff, norm_mul, norm_mul, hexp, mul_one]
  have hpow := norm_pow_sub_pow_le hX hY k
  have hu2 : u ^ 2 ≤ (2 * R + 1) ^ 2 := by
    have : |u| ≤ 2 * R + 1 := by linarith
    calc u ^ 2 = |u| ^ 2 := (sq_abs u).symm
      _ ≤ (2 * R + 1) ^ 2 := by gcongr
  calc ‖r.coef j k‖ * ‖X ^ k - Y ^ k‖
      ≤ ‖r.coef j k‖ * (k * (2 * R + 1) ^ (k - 1) * (ε * u ^ 2)) := by
        gcongr
        exact hpow.trans (by gcongr)
    _ ≤ ‖r.coef j k‖ * (k * (2 * R + 1) ^ (k - 1) * (ε * (2 * R + 1) ^ 2)) := by gcongr
    _ ≤ ε * (‖r.coef j k‖ * (k * (2 * R + 1) ^ (k + 1))) := by
        rcases k with _ | k
        · simp
        · simp only [Nat.add_sub_cancel]
          apply le_of_eq; ring

end ExpSum

/-- **Turán–Nazarov for exponential sums with jets** (real nodes; total multiplicity `M ≥ 1`). -/
theorem turan_expSum (hT : Nazarov1993_Turan) : ∃ c : ℝ, 0 < c ∧
    ∀ (r : ExpSum) (a b : ℝ) (Ω : Set ℝ) (K : ℝ),
      1 ≤ r.M → a < b → MeasurableSet Ω → Ω ⊆ Icc a b → 0 < volume Ω →
      (∀ t ∈ Ω, ‖r.eval t‖ ≤ K) →
      ∀ t ∈ Icc a b, ‖r.eval t‖ ≤ (c * (b - a) / (volume Ω).toReal) ^ (r.M - 1) * K := by
  obtain ⟨c, hc, hTr⟩ := turan_real hT
  refine ⟨c, hc, fun r a b Ω K hM hab hΩ hΩs hΩpos hK t ht => ?_⟩
  set F := (c * (b - a) / (volume Ω).toReal) ^ (r.M - 1) with hF
  have hF0 : 0 ≤ F := by
    apply pow_nonneg; apply div_nonneg (by nlinarith) ENNReal.toReal_nonneg
  set R : ℝ := |a| + |b| with hRdef
  have hR : 0 ≤ R := by positivity
  have hmem : ∀ u ∈ Icc a b, |u| ≤ R := by
    intro u hu; rw [abs_le]; constructor <;> cases abs_cases a <;> cases abs_cases b <;>
      linarith [hu.1, hu.2]
  set C₀ := ∑ j, ∑ k ∈ Finset.range (r.mult j), ‖r.coef j k‖ * (k * (2 * R + 1) ^ (k + 1))
  have hC₀ : 0 ≤ C₀ := by positivity
  -- for every small ε: ‖eval t‖ ≤ F (K + C₀ ε) + C₀ ε
  have key : ∀ ε : ℝ, 0 < ε → ε * (R + 1) ≤ 1 → ‖r.eval t‖ ≤ F * (K + ε * C₀) + ε * C₀ := by
    intro ε hε hεR
    have herr : ∀ u ∈ Icc a b, ‖r.approx ε u - r.eval u‖ ≤ ε * C₀ :=
      fun u hu => r.norm_approx_sub_eval_le hR hε hεR (hmem u hu)
    have hcard : 1 ≤ Fintype.card (Σ j : Fin r.J, Fin (r.mult j)) := by rw [r.card_sigma]; exact hM
    have hΩ' : ∀ u ∈ Ω, ‖∑ x : (Σ j : Fin r.J, Fin (r.mult j)), r.acoef ε x.1 x.2 *
        cexp (I * ((r.node x.1 + (x.2 : ℕ) * ε : ℝ) : ℂ) * u)‖ ≤ K + ε * C₀ := by
      intro u hu
      rw [← r.approx_eq_sigma]
      have h0 := norm_add_le (r.eval u) (r.approx ε u - r.eval u)
      rw [add_sub_cancel] at h0
      exact h0.trans (add_le_add (hK u hu) (herr u (hΩs hu)))
    have h := hTr (ι := Σ j : Fin r.J, Fin (r.mult j))
      (fun x => r.acoef ε x.1 x.2) (fun x => r.node x.1 + (x.2 : ℕ) * ε) a b Ω
      (K + ε * C₀) hcard hab hΩ hΩs hΩpos hΩ' t ht
    rw [← r.approx_eq_sigma, r.card_sigma] at h
    have h0 := norm_sub_le (r.approx ε t) (r.approx ε t - r.eval t)
    rw [sub_sub_cancel] at h0
    exact h0.trans (add_le_add h (herr t ht))
  -- let ε → 0
  have hlim : Tendsto (fun ε : ℝ => F * (K + ε * C₀) + ε * C₀) (𝓝[>] 0) (𝓝 (F * K)) := by
    have : Tendsto (fun ε : ℝ => F * (K + ε * C₀) + ε * C₀) (𝓝 0) (𝓝 (F * (K + 0 * C₀) + 0 * C₀)) :=
      (by fun_prop : Continuous fun ε : ℝ => F * (K + ε * C₀) + ε * C₀).tendsto 0
    simpa using this.mono_left nhdsWithin_le_nhds
  refine ge_of_tendsto hlim ?_
  filter_upwards [Ioo_mem_nhdsGT (show (0:ℝ) < 1 / (R + 1) by positivity)] with ε hε
  have h2 : ε * (R + 1) ≤ 1 := by
    have := hε.2
    rw [lt_div_iff₀ (by positivity)] at this
    linarith
  exact key ε hε.1 h2

end Riemann.Deep.LeanNBRate
