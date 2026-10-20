import Mathlib.Analysis.InnerProductSpace.l2Space
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.PosDef

/-!
# Yakaboylu, arXiv 2408.15135 ("On the Existence of the Hilbert–Pólya Hamiltonian",
# v7–v11; retitled "Nontrivial Riemann Zeros as Spectrum", v12–v17), laid down and pinned

Setting (all versions): non-normal operator `R = -D - i μ(T)` on `L²[0,∞)`; eigenvectors
`Ψ_ρ` of `R` (eigenvalue `i(1/2-ρ)`), `Φ_ρ` of `R†` (eigenvalue `-i(1/2-ρ̄)`), biorthogonal
`⟨Φ_ρ'|Ψ_ρ⟩ = δ`. Claimed RH route: an intertwiner `W ≥ 0`, `W R = R† W`, then
`h = W^{1/2} R W^{-1/2}` self-adjoint with eigenvalues `i(1/2-ρ)` ⇒ `Re ρ = 1/2`.

* v11 (eqs. (19), (23), (24)): `W := Σ_ρ |Φ_ρ⟩⟨Φ_ρ|`. In any biorthogonal eigen-system
  `W R Ψ_ρ = λ_ρ Φ_ρ` but `R† W Ψ_ρ = conj λ_ρ Φ_ρ`, so `W R = R† W` iff every eigenvalue is
  real, i.e. iff RH (`v11_intertwine_iff`). The paper's only justification, that the eigenvalue
  SETS of `R` and `R†` coincide, does not give it (`v11_lemma_false`). Circular.
* v16 (Lemma "Positive Semidefiniteness of W"): `W = V K`, `V = ω(T)^{-2} > 0`, `V K = K† V`,
  rewritten as `W = (K^{1/2})† V K^{1/2} ≥ 0` using a "conjugation-symmetric branch" of `√`.
  No such branch exists (`no_conj_symm_sqrt`), and `V`-symmetric `K` need not give `W ≥ 0`
  (`v16_lemma_false`); with `V = 1` no square root `M` of `K` has `Mᴴ V M = K`.
* v11 "complete basis / resolution of identity on the closed span" and v17 Cor. 4.2
  ("`P = Σ |Ψ_ρ⟩⟨Φ_ρ|`, strong limit, identity on the closed span"): false for general
  biorthogonal systems in ℓ² (`resolution_false`).
v17 (latest) no longer claims RH: it states `W > 0 ⟺ RH`.
-/

namespace Riemann.Deep.Yakaboylu

open Matrix
open scoped ComplexOrder

/-! ## v11: the intertwining identity is equivalent to reality of the spectrum -/

section V11

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- Finite model of v11. Columns of `P` are `Ψ_k`, columns of `(P⁻¹)ᴴ` are `Φ_k`
(biorthogonal: `((P⁻¹)ᴴ)ᴴ * P = 1`), `R = P D P⁻¹` has `R Ψ_k = d_k Ψ_k`, `Rᴴ Φ_k = conj d_k Φ_k`,
and `W = Σ_k Φ_k Φ_kᴴ = (P⁻¹)ᴴ P⁻¹` (v11 eq. (23)). Then `W R = Rᴴ W ↔ ∀ k, d_k ∈ ℝ`. -/
theorem v11_intertwine_iff (P : Matrix n n ℂ) (hP : IsUnit P.det) (d : n → ℂ) :
    ((P⁻¹)ᴴ * P⁻¹) * (P * diagonal d * P⁻¹) = (P * diagonal d * P⁻¹)ᴴ * ((P⁻¹)ᴴ * P⁻¹) ↔
      ∀ i, star (d i) = d i := by
  have h1 : P⁻¹ * P = 1 := nonsing_inv_mul P hP
  have h2 : Pᴴ * (P⁻¹)ᴴ = 1 := by rw [← conjTranspose_mul, h1, conjTranspose_one]
  have hL : ((P⁻¹)ᴴ * P⁻¹) * (P * diagonal d * P⁻¹) = (P⁻¹)ᴴ * diagonal d * P⁻¹ := by
    simp only [Matrix.mul_assoc]; rw [← Matrix.mul_assoc P⁻¹ P, h1, Matrix.one_mul]
  have hR : (P * diagonal d * P⁻¹)ᴴ * ((P⁻¹)ᴴ * P⁻¹) = (P⁻¹)ᴴ * diagonal (star d) * P⁻¹ := by
    simp only [conjTranspose_mul, diagonal_conjTranspose, Matrix.mul_assoc]
    rw [← Matrix.mul_assoc Pᴴ, h2, Matrix.one_mul]
  rw [hL, hR]
  have key : ∀ D : Matrix n n ℂ, Pᴴ * ((P⁻¹)ᴴ * D * P⁻¹) * P = D := by
    intro D
    rw [show Pᴴ * ((P⁻¹)ᴴ * D * P⁻¹) * P = (Pᴴ * (P⁻¹)ᴴ) * D * (P⁻¹ * P) by
      simp only [Matrix.mul_assoc], h2, h1, Matrix.one_mul, Matrix.mul_one]
  constructor
  · intro h i
    have hd : diagonal d = diagonal (star d) := by rw [← key (diagonal d), h, key]
    have := congrFun (congrFun hd i) i
    simpa [diagonal_apply_eq] using this.symm
  · intro h
    have : star d = d := funext h
    rw [this]

/-- The general lemma v11 relies on: biorthogonal eigen-systems whose eigenvalue set is closed
under conjugation (`spec R = spec R†`, v11 eq. (21) "the sets intersect") admit the intertwiner
`W = Σ|Φ⟩⟨Φ|`. -/
def V11Lemma : Prop :=
  ∀ (P : Matrix (Fin 2) (Fin 2) ℂ) (d : Fin 2 → ℂ), IsUnit P.det →
    (∀ i, ∃ j, star (d i) = d j) →
    ((P⁻¹)ᴴ * P⁻¹) * (P * diagonal d * P⁻¹) = (P * diagonal d * P⁻¹)ᴴ * ((P⁻¹)ᴴ * P⁻¹)

/-- Counterexample: `d = (i, -i)` (a reflection pair `ρ₀`, `1-ρ̄₀` off the line), `P = 1`. -/
theorem v11_lemma_false : ¬ V11Lemma := by
  intro h
  have hd : ∀ i, ∃ j, star (![Complex.I, -Complex.I] i) = ![Complex.I, -Complex.I] j := by
    intro i; fin_cases i
    · exact ⟨1, by simp [Complex.conj_I]⟩
    · exact ⟨0, by simp [Complex.conj_I]⟩
  have := (v11_intertwine_iff (1 : Matrix (Fin 2) (Fin 2) ℂ) (by simp) _).1
    (h 1 _ (by simp) hd) 0
  have h2 := congrArg Complex.im this
  simp [Complex.conj_I] at h2
  norm_num at h2

end V11

/-! ## v16: the positivity lemma -/

/-- No branch of `√` commutes with conjugation (needed for `(K^{1/2})† = (K†)^{1/2}`):
it fails already at `-1`. -/
theorem no_conj_symm_sqrt :
    ¬ ∃ f : ℂ → ℂ, (∀ z, f z ^ 2 = z) ∧ ∀ z, f (star z) = star (f z) := by
  rintro ⟨f, hsq, hc⟩
  have hw : star (f (-1)) = f (-1) := by simpa using (hc (-1)).symm
  have him : (f (-1)).im = 0 := Complex.conj_eq_iff_im.1 hw
  have hre := congrArg Complex.re (hsq (-1))
  simp only [pow_two, Complex.mul_re, him, mul_zero, sub_zero] at hre
  norm_num at hre
  nlinarith [mul_self_nonneg (f (-1)).re]

/-- v16's general claim: `V > 0`, `V K = Kᴴ V` ⇒ `W = V K ⪰ 0`. -/
def V16Lemma : Prop :=
  ∀ (V K : Matrix (Fin 2) (Fin 2) ℂ), V.PosDef → V * K = Kᴴ * V →
    ∀ x : Fin 2 → ℂ, 0 ≤ (star x ⬝ᵥ (V * K) *ᵥ x).re

theorem v16_lemma_false : ¬ V16Lemma := by
  intro h
  have := h 1 (diagonal ![1, -1]) PosDef.one (by simp [diagonal_conjTranspose])
    (Pi.single 1 1)
  simp [dotProduct, mulVec] at this
  norm_num at this

/-- The rewriting `W = Mᴴ V M` (`M = K^{1/2}`, `V = 1`) is impossible for every matrix `M`. -/
theorem v16_sqrt_rewrite_false (M : Matrix (Fin 2) (Fin 2) ℂ) :
    Mᴴ * M ≠ diagonal ![1, -1] := by
  intro h
  have := congrArg Complex.re (congrFun (congrFun h 1) 1)
  simp [Matrix.mul_apply, Fin.sum_univ_two] at this
  nlinarith [sq_nonneg (M 0 1).re, sq_nonneg (M 0 1).im, sq_nonneg (M 1 1).re,
    sq_nonneg (M 1 1).im]

/-! ## Resolution of identity from biorthogonality alone (v11 after (19); v17 Cor. 4.2) -/

section L2

open scoped InnerProductSpace
open Filter Topology Finset

/-- `ℓ²(ℕ)`, a separable Hilbert space, unitarily equivalent to `L²[0,∞)`. -/
noncomputable abbrev E := lp (fun _ : ℕ => ℂ) 2

/-- Standard basis. -/
noncomputable def e (k : ℕ) : E := lp.single 2 k (1 : ℂ)

theorem inner_e (a b : ℕ) : ⟪e a, e b⟫_ℂ = if a = b then 1 else 0 := by
  rw [e, lp.inner_single_left]
  simp only [e, lp.single_apply]
  by_cases h : a = b
  · subst h; simp
  · simp [h]

/-- `Ψ_n = e_n - e_{n+1}`. -/
noncomputable def Ψ (n : ℕ) : E := e n - e (n + 1)

/-- `Φ_n = e_0 + ⋯ + e_n`. -/
noncomputable def Φ (n : ℕ) : E := ∑ k ∈ range (n + 1), e k

theorem biorth (m n : ℕ) : ⟪Φ m, Ψ n⟫_ℂ = if m = n then 1 else 0 := by
  simp only [Φ, Ψ, sum_inner, inner_sub_right, inner_e, sum_ite_eq',
    mem_range]
  by_cases h1 : m = n
  · subst h1; simp
  · by_cases h2 : n < m
    · simp [h1, show n + 1 < m + 1 by omega, show n < m + 1 by omega]
    · simp [h1, show ¬ n + 1 < m + 1 by omega, show ¬ n < m + 1 by omega]

/-- Partial sums of the claimed resolution of identity `Σ |Ψ_n⟩⟨Φ_n|`. -/
noncomputable def S (N : ℕ) (x : E) : E := ∑ n ∈ range N, ⟪Φ n, x⟫_ℂ • Ψ n

theorem S_e0 (N : ℕ) : S N (e 0) = e 0 - e N := by
  have : ∀ n, ⟪Φ n, e 0⟫_ℂ = 1 := by
    intro n; simp [Φ, sum_inner, inner_e]
  simp only [S, this, one_smul, Ψ]
  exact sum_range_sub' (fun k => e k) N

theorem not_cauchy : ¬ CauchySeq (fun N => S N (e 0)) := by
  rw [Metric.cauchySeq_iff]
  push Not
  refine ⟨1, one_pos, fun N => ⟨N, le_rfl, N + 1, by omega, ?_⟩⟩
  rw [S_e0, S_e0, dist_eq_norm, show e 0 - e N - (e 0 - e (N + 1)) = e (N + 1) - e N by abel]
  calc (1 : ℝ) = ‖(e (N + 1) - e N) (N + 1)‖ := by
        simp [e, lp.single_apply]
    _ ≤ ‖e (N + 1) - e N‖ := lp.norm_apply_le_norm (by norm_num) _ _

theorem norm_sq_sum (N : ℕ) : ‖∑ j ∈ range N, e (j + 1)‖ ^ 2 = N := by
  have h : ⟪∑ j ∈ range N, e (j + 1), ∑ j ∈ range N, e (j + 1)⟫_ℂ = (N : ℂ) := by
    simp only [sum_inner, inner_sum, inner_e, add_left_inj]
    rw [sum_congr rfl fun j hj => (sum_ite_eq' (range N) j fun _ => (1 : ℂ)).trans (ite_eq_left hj)]
    simp
  have h' := @inner_self_eq_norm_sq_to_K ℂ _ _ _ _ (∑ j ∈ range N, e (j + 1))
  rw [h] at h'
  apply Complex.ofReal_injective
  push_cast
  exact h'.symm

theorem e0_mem_closure : e 0 ∈ closure (Submodule.span ℂ (Set.range Ψ) : Set E) := by
  rw [Metric.mem_closure_iff]
  intro ε hε
  obtain ⟨N, hN⟩ := exists_nat_gt (1 / ε ^ 2)
  have hN0 : (0 : ℝ) < N := lt_trans (by positivity) hN
  refine ⟨((N : ℂ)⁻¹) • ∑ j ∈ range N, (e 0 - e (j + 1)), ?_, ?_⟩
  · refine Submodule.smul_mem _ _ (Submodule.sum_mem _ fun j _ => ?_)
    rw [← sum_range_sub' (fun k => e k) (j + 1)]
    exact Submodule.sum_mem _ fun k _ => Submodule.subset_span ⟨k, rfl⟩
  · have hsum : e 0 - (N : ℂ)⁻¹ • ∑ j ∈ range N, (e 0 - e (j + 1)) =
        (N : ℂ)⁻¹ • ∑ j ∈ range N, e (j + 1) := by
      rw [sum_sub_distrib, sum_const, card_range, smul_sub, ← Nat.cast_smul_eq_nsmul ℂ,
        smul_smul, inv_mul_cancel₀ (by exact_mod_cast hN0.ne'), one_smul]
      abel
    rw [dist_eq_norm, hsum, norm_smul]
    have h2 : (‖(N : ℂ)⁻¹‖ * ‖∑ j ∈ range N, e (j + 1)‖) ^ 2 = 1 / N := by
      rw [mul_pow, norm_sq_sum, norm_inv, Complex.norm_natCast]; field_simp
    have h3 : (1 : ℝ) / N < ε ^ 2 := by
      rw [div_lt_iff₀ hN0]; rw [div_lt_iff₀ (by positivity)] at hN; linarith
    nlinarith [norm_nonneg ((N : ℂ)⁻¹), norm_nonneg (∑ j ∈ range N, e (j + 1)),
      mul_nonneg (norm_nonneg ((N : ℂ)⁻¹)) (norm_nonneg (∑ j ∈ range N, e (j + 1)))]

/-- v11/v17 claim: for a biorthogonal pair, `Σ_n |Ψ_n⟩⟨Φ_n|` converges strongly to the identity
on the closed span of the `Ψ_n`. -/
def ResolutionOfIdentity : Prop :=
  ∀ (A B : ℕ → E), (∀ m n, ⟪B m, A n⟫_ℂ = if m = n then 1 else 0) →
    ∀ x ∈ closure (Submodule.span ℂ (Set.range A) : Set E),
      Tendsto (fun N => ∑ n ∈ range N, ⟪B n, x⟫_ℂ • A n) atTop (𝓝 x)

theorem resolution_false : ¬ ResolutionOfIdentity := fun h =>
  not_cauchy (h Ψ Φ biorth (e 0) e0_mem_closure).cauchySeq

end L2

end Riemann.Deep.Yakaboylu

#print axioms Riemann.Deep.Yakaboylu.v11_intertwine_iff
#print axioms Riemann.Deep.Yakaboylu.v11_lemma_false
#print axioms Riemann.Deep.Yakaboylu.no_conj_symm_sqrt
#print axioms Riemann.Deep.Yakaboylu.v16_lemma_false
#print axioms Riemann.Deep.Yakaboylu.v16_sqrt_rewrite_false
#print axioms Riemann.Deep.Yakaboylu.resolution_false
