module

public import Riemann.Deep.LeanLi.LiCoeff
public import Mathlib.Analysis.SpecialFunctions.Trigonometric.Arctan

/-!
# LeanLi: Li's formula under RH: `λ_n = Σ_{γ>0} 2 m_γ (1 - cos nθ_γ)`, `θ(t) = π - 2 arctan 2t`
-/

@[expose] public section

open Complex Filter Topology Set ComplexConjugate MeasureTheory

namespace LeanLi

/-- `θ(t) = 2 arctan(1/(2t))` for `t > 0`, written as `π - 2 arctan(2t)` (smooth, `θ(0) = π`). -/
noncomputable def liTheta (t : ℝ) : ℝ := Real.pi - 2 * Real.arctan (2 * t)

lemma re_eq_half (hRH : RiemannHypothesis) {ρ : ℂ} (hρ : ρ ∈ NontrivialZeros) : ρ.re = 1 / 2 :=
  hRH ρ hρ.1 hρ.2.1 hρ.2.2

/-- Under RH, `1 - 1/ρ = e^{iθ(γ)}` for `ρ = ½ + iγ`. -/
lemma one_sub_inv_eq_exp {γ : ℝ} :
    (1 : ℂ) - 1 / ((1 / 2 : ℂ) + γ * I) = Complex.exp (liTheta γ * I) := by
  have hpos : (0 : ℝ) < 1 + (2 * γ) ^ 2 := by positivity
  have hs := Real.sqrt_pos.2 hpos
  have hss : √(1 + (2 * γ) ^ 2) * √(1 + (2 * γ) ^ 2) = 1 + (2 * γ) ^ 2 := Real.mul_self_sqrt hpos.le
  have hd : (0 : ℝ) < 4 * γ ^ 2 + 1 := by positivity
  have hc : Real.cos (liTheta γ) = (4 * γ ^ 2 - 1) / (4 * γ ^ 2 + 1) := by
    unfold liTheta
    rw [Real.cos_pi_sub, Real.cos_two_mul, Real.cos_arctan, div_pow, one_pow, Real.sq_sqrt hpos.le]
    field_simp
    ring
  have hsn : Real.sin (liTheta γ) = 4 * γ / (4 * γ ^ 2 + 1) := by
    unfold liTheta
    rw [Real.sin_pi_sub, Real.sin_two_mul, Real.sin_arctan, Real.cos_arctan]
    rw [show (2 : ℝ) * (2 * γ / √(1 + (2 * γ) ^ 2)) * (1 / √(1 + (2 * γ) ^ 2)) =
        4 * γ / (√(1 + (2 * γ) ^ 2) * √(1 + (2 * γ) ^ 2)) by field_simp; ring, hss]
    congr 1; ring
  rw [Complex.exp_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin, hc, hsn]
  set z : ℂ := (1 / 2 : ℂ) + γ * I with hz
  have hz0 : z ≠ 0 := by
    intro h
    have := congrArg Complex.re h
    simp [hz] at this
  have h12re : ((1 : ℂ) / 2).re = 1 / 2 := by norm_num
  have h12im : ((1 : ℂ) / 2).im = 0 := by norm_num
  have key : (((4 * γ ^ 2 - 1) / (4 * γ ^ 2 + 1) : ℝ) + ((4 * γ / (4 * γ ^ 2 + 1) : ℝ) : ℂ) * I) * z =
      z - 1 := by
    apply Complex.ext
    · simp only [hz, Complex.mul_re, Complex.add_re, Complex.add_im, Complex.ofReal_re,
        Complex.ofReal_im, Complex.I_re, Complex.I_im, Complex.mul_im, Complex.sub_re, Complex.one_re,
        h12re, h12im]
      field_simp
      ring
    · simp only [hz, Complex.mul_re, Complex.add_re, Complex.add_im, Complex.ofReal_re,
        Complex.ofReal_im, Complex.I_re, Complex.I_im, Complex.mul_im, Complex.sub_im, Complex.one_im,
        h12re, h12im]
      field_simp
      ring
  have h1 : (1 : ℂ) - 1 / z = (z - 1) / z := by field_simp
  rw [h1, ← key, mul_div_assoc, div_self hz0, mul_one]

/-- Under RH the paired term is `2 m (1 - cos nθ(γ))`. -/
lemma upperTerm_eq (hRH : RiemannHypothesis) (n : ℕ) {ρ : ℂ} (hρ : ρ ∈ NontrivialZeros) :
    upperTerm n ρ = 2 * mult ρ * (1 - Real.cos (n * liTheta ρ.im)) := by
  have hρeq : ρ = (1 / 2 : ℂ) + (ρ.im : ℂ) * I := by
    apply Complex.ext <;> simp [re_eq_half hRH hρ]
  have h1 : (1 : ℂ) - 1 / ρ = Complex.exp (liTheta ρ.im * I) := by
    have := one_sub_inv_eq_exp (γ := ρ.im)
    rwa [← hρeq] at this
  unfold upperTerm
  congr 1
  rw [h1, ← Complex.exp_nat_mul]
  have : (n : ℂ) * ((liTheta ρ.im : ℂ) * I) = ((n * liTheta ρ.im : ℝ) : ℂ) * I := by
    push_cast; ring
  rw [this, Complex.sub_re, Complex.one_re, Complex.exp_ofReal_mul_I_re]

/-- The zeros with `γ > 0` listed with multiplicity. -/
abbrev ZeroIdx : Type := Σ ρ : UpperZeros, Fin (mult (ρ : ℂ))

/-- Ordinate of an indexed zero. -/
noncomputable def ordK (k : ZeroIdx) : ℝ := (k.1 : ℂ).im

/-- Phase `θ_γ` of an indexed zero. -/
noncomputable def thK (k : ZeroIdx) : ℝ := liTheta (ordK k)

/-- The upper zeros with `γ ≤ T` as a finset of the subtype. -/
noncomputable def upperUpTo (T : ℝ) : Finset UpperZeros :=
  ((finite_zerosUpTo T).preimage (f := fun ρ : UpperZeros => (ρ : ℂ)) Subtype.val_injective.injOn).toFinset

lemma mem_upperUpTo {T : ℝ} {ρ : UpperZeros} : ρ ∈ upperUpTo T ↔ (ρ : ℂ).im ≤ T := by
  simp [upperUpTo, zerosUpTo, ρ.2.1, ρ.2.2]

/-- `Z(T)`: the indexed zeros with `0 < γ ≤ T`. -/
noncomputable def ZK (T : ℝ) : Finset ZeroIdx := (upperUpTo T).sigma fun _ => Finset.univ

lemma mem_ZK {T : ℝ} {k : ZeroIdx} : k ∈ ZK T ↔ ordK k ≤ T := by
  simp [ZK, ordK, mem_upperUpTo]

lemma upperUpTo_map (T : ℝ) : (upperUpTo T).map (Function.Embedding.subtype _) = zset T := by
  ext z
  simp only [Finset.mem_map, Function.Embedding.coe_subtype, mem_zset, mem_upperUpTo]
  constructor
  · rintro ⟨ρ, hρ, rfl⟩; exact ⟨ρ.2.1, ρ.2.2, hρ⟩
  · rintro ⟨h1, h2, h3⟩; exact ⟨⟨z, h1, h2⟩, h3, rfl⟩

lemma card_ZK (T : ℝ) : ((ZK T).card : ℝ) = zeroCount T := by
  rw [zeroCount_eq, ← upperUpTo_map, Finset.sum_map]
  simp [ZK, Finset.card_sigma]

/-- **Li's formula under RH**, in the form `λ_n = Σ_k 2(1 - cos nθ_k)` over the zeros with `γ > 0`
counted with multiplicity. -/
theorem hasSum_lam (L : LitCore) (hRH : RiemannHypothesis) (n : ℕ) :
    HasSum (fun k : ZeroIdx => 2 * (1 - Real.cos (n * thK k))) (lam n) := by
  have hnn : ∀ k : ZeroIdx, 0 ≤ 2 * (1 - Real.cos (n * thK k)) := fun k => by
    have := Real.cos_le_one (n * thK k); linarith
  have hfib : ∀ ρ : UpperZeros, HasSum (fun i : Fin (mult (ρ : ℂ)) =>
      2 * (1 - Real.cos (n * thK ⟨ρ, i⟩))) (upperTerm n ρ) := by
    intro ρ
    rw [upperTerm_eq hRH n ρ.2.1]
    have := hasSum_fintype (fun i : Fin (mult (ρ : ℂ)) => 2 * (1 - Real.cos (n * thK ⟨ρ, i⟩)))
    convert this using 1
    simp only [thK, ordK, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
    ring
  have hsum : Summable (fun k : ZeroIdx => 2 * (1 - Real.cos (n * thK k))) := by
    rw [summable_sigma_of_nonneg hnn]
    refine ⟨fun ρ => (hfib ρ).summable, ?_⟩
    refine (summable_upperTerm L n).congr fun ρ => ?_
    exact ((hfib ρ).tsum_eq).symm
  have h1 := hsum.hasSum
  have h2 := h1.sigma hfib
  rw [lam_eq L n]
  have h3 := (summable_upperTerm L n).hasSum
  rw [h2.unique h3] at h1
  exact h1

end LeanLi
