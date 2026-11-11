import Riemann.Deep.KimEff.Lean.Theorem
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.Analysis.Complex.TaylorSeries

/-!
# Explicit Kim–Lee wedge, part 5: Riemann's ξ

ξ(s) = s(s-1)π^{-s/2}Γ(s/2)ζ(s)/2 from Mathlib's completed zeta; the Taylor coefficients γ(m) of
8ξ(1/2+z) = Σ γ(m) z^{2m}/m!; Θ^{(n)}(w) = Σ_k γ(n+k) w^k/k!; proof that Θ(z^2) = 8ξ(1/2+z) and
that (I2), in the form of Bebel's half-plane, puts every zero of Θ at z^2 with |Re z| ≤ a.
Main theorem: xi_jensen_hyperbolic.
-/

namespace KimEff

noncomputable section

/-- Riemann's ξ(s) = s(s-1)π^{-s/2}Γ(s/2)ζ(s)/2 = (s(s-1)Λ₀(s) + 1)/2 (entire). -/
def riemannXi (s : ℂ) : ℂ := (s * (s - 1) * completedRiemannZeta₀ s + 1) / 2

lemma riemannXi_eq_completed {s : ℂ} (hs0 : s ≠ 0) (hs1 : s ≠ 1) :
    riemannXi s = s * (s - 1) * completedRiemannZeta s / 2 := by
  rw [completedRiemannZeta_eq, riemannXi]
  have h1 : s - 1 ≠ 0 := sub_ne_zero.mpr hs1
  have h2 : 1 - s ≠ 0 := sub_ne_zero.mpr (Ne.symm hs1)
  field_simp
  ring

lemma riemannXi_one_sub (s : ℂ) : riemannXi (1 - s) = riemannXi s := by
  unfold riemannXi
  rw [completedRiemannZeta₀_one_sub]
  ring

lemma differentiable_riemannXi : Differentiable ℂ riemannXi := by
  unfold riemannXi
  have h := differentiable_completedZeta₀
  fun_prop

/-- (I2): Bebel's Lean-checked half-plane (QuasiRiemannTracker): ζ(s) ≠ 0 for Re s > 1/2 + a. -/
def HalfPlane (a : ℝ) : Prop := ∀ s : ℂ, 1 / 2 + a < s.re → riemannZeta s ≠ 0

lemma riemannXi_ne_zero_of_re_gt {a : ℝ} (ha : 0 ≤ a) (hH : HalfPlane a) {s : ℂ}
    (hs : 1 / 2 + a < s.re) : riemannXi s ≠ 0 := by
  have hre : 0 < s.re := by linarith
  have hs0 : s ≠ 0 := by
    intro h
    rw [h, Complex.zero_re] at hre
    exact lt_irrefl _ hre
  by_cases hs1 : s = 1
  · rw [hs1]
    unfold riemannXi
    norm_num
  rw [riemannXi_eq_completed hs0 hs1]
  have hz : riemannZeta s ≠ 0 := hH s hs
  rw [riemannZeta_def_of_ne_zero hs0] at hz
  have hΛ : completedRiemannZeta s ≠ 0 := by
    intro h
    rw [h, zero_div] at hz
    exact hz rfl
  have h1 : s - 1 ≠ 0 := sub_ne_zero.mpr hs1
  exact div_ne_zero (mul_ne_zero (mul_ne_zero hs0 h1) hΛ) two_ne_zero

/-- (I2) + functional equation: every zero of ξ has |Re s - 1/2| ≤ a. -/
lemma abs_re_sub_half_le {a : ℝ} (ha : 0 ≤ a) (hH : HalfPlane a) {s : ℂ}
    (hs : riemannXi s = 0) : |s.re - 1 / 2| ≤ a := by
  rw [abs_le]
  constructor
  · by_contra h
    have h := not_le.mp h
    have h1 : 1 / 2 + a < (1 - s).re := by
      rw [Complex.sub_re, Complex.one_re]
      linarith
    exact riemannXi_ne_zero_of_re_gt ha hH h1 (by rw [riemannXi_one_sub]; exact hs)
  · by_contra h
    have h := not_le.mp h
    exact riemannXi_ne_zero_of_re_gt ha hH (by linarith) hs

/-- f(z) = 8 ξ(1/2 + z). -/
def xiF (z : ℂ) : ℂ := 8 * riemannXi (1 / 2 + z)

lemma differentiable_xiF : Differentiable ℂ xiF := by
  unfold xiF
  have h := differentiable_riemannXi
  fun_prop

lemma xiF_neg (z : ℂ) : xiF (-z) = xiF z := by
  unfold xiF
  have h := riemannXi_one_sub (1 / 2 + z)
  rw [show (1 : ℂ) - (1 / 2 + z) = 1 / 2 + -z by ring] at h
  rw [h]

lemma iteratedDeriv_xiF_odd (k : ℕ) : iteratedDeriv (2 * k + 1) xiF 0 = 0 := by
  have h := iteratedDeriv_comp_neg (2 * k + 1) xiF 0
  have hf : (fun x => xiF (-x)) = xiF := funext xiF_neg
  have hp : ((-1 : ℂ)) ^ (2 * k + 1) = -1 := by
    rw [pow_succ, pow_mul, neg_one_sq, one_pow, one_mul]
  rw [hf, neg_zero, hp, smul_eq_mul] at h
  linear_combination (1 / 2 : ℂ) * h

/-- γ_C(m) = m!/(2m)! f^{(2m)}(0), so that 8ξ(1/2+z) = Σ γ(m) z^{2m}/m!. -/
def xiGammaC (m : ℕ) : ℂ :=
  (m.factorial : ℂ) / ((2 * m).factorial : ℂ) * iteratedDeriv (2 * m) xiF 0

/-- γ(m) := Re γ_C(m) (γ_C(m) is real since ξ is real on the real axis; this is not used). -/
def xiGamma (m : ℕ) : ℝ := (xiGammaC m).re

/-- Θ^{(n)}(w) = Σ_k γ(n+k) w^k/k!, the n-th derivative of Θ(w) = 8ξ(1/2 + √w). -/
def ThetaDeriv (n : ℕ) (w : ℂ) : ℂ := ∑' k : ℕ, xiGammaC (n + k) * w ^ k / (k.factorial : ℂ)

/-- Θ(z^2) = 8ξ(1/2 + z) (Taylor series of the entire even function f). -/
lemma ThetaDeriv_zero_sq (z : ℂ) : ThetaDeriv 0 (z ^ 2) = xiF z := by
  have h := Complex.hasSum_taylorSeries_of_entire differentiable_xiF 0 z
  have hodd : ∀ x ∉ Set.range (fun k : ℕ => 2 * k),
      (fun n : ℕ => ((n.factorial : ℂ))⁻¹ • (z - 0) ^ n • iteratedDeriv n xiF 0) x = 0 := by
    intro x hx
    rcases Nat.even_or_odd x with ⟨k, hk⟩ | ⟨k, hk⟩
    · exact absurd (show x ∈ Set.range (fun k : ℕ => 2 * k) from ⟨k, by show 2 * k = x; omega⟩) hx
    · simp only
      rw [hk, iteratedDeriv_xiF_odd, smul_zero, smul_zero]
  have hinj : Function.Injective (fun k : ℕ => 2 * k) := fun a b hab => by
    simp only at hab
    omega
  have h2 := (hinj.hasSum_iff hodd).mpr h
  unfold ThetaDeriv
  rw [← h2.tsum_eq]
  congr 1
  funext k
  simp only [Function.comp_apply, xiGammaC, zero_add, sub_zero, smul_eq_mul]
  have hk : (k.factorial : ℂ) ≠ 0 := by exact_mod_cast k.factorial_ne_zero
  have h2k : ((2 * k).factorial : ℂ) ≠ 0 := by exact_mod_cast (2 * k).factorial_ne_zero
  rw [← pow_mul]
  field_simp

/-- (2.1) for Riemann's Θ: from (I2), every zero of Θ is z^2 with |Re z| ≤ a. -/
theorem zerosSq_xi {a : ℝ} (ha : 0 ≤ a) (hH : HalfPlane a) : ZerosSq a ThetaDeriv := by
  intro w hw
  obtain ⟨z, hz⟩ := IsAlgClosed.exists_pow_nat_eq w (by norm_num : 0 < 2)
  refine ⟨z, hz, ?_⟩
  rw [← hz, ThetaDeriv_zero_sq] at hw
  have hxi : riemannXi (1 / 2 + z) = 0 := by
    unfold xiF at hw
    exact (mul_eq_zero.mp hw).resolve_left (by norm_num)
  have h := abs_re_sub_half_le ha hH hxi
  have hre : (1 / 2 + z).re - 1 / 2 = z.re := by
    rw [Complex.add_re]
    norm_num
  rwa [hre] at h

/-- Main theorem for Riemann's ξ (a = 0.37338196661320752, n_* = 998,616,377,771): the named
inputs imply (ii) J^{d,n} hyperbolic for d ≤ D_n, (iii) for d ≤ 6.45 n^2/log^2 n, and the
power form n ≥ max(n_*, 1.1853e-10 d). -/
theorem xi_jensen_hyperbolic (hH : HalfPlane a0) (hO : ObreschkoffKimLee)
    (hHad : ∀ n, HadamardApprox xiGamma ThetaDeriv n)
    (hGL : ∀ n, GaussLucasGenus0 ThetaDeriv n) (hpos : ∀ m, 0 < xiGamma m)
    (hi : ZeroFreeRegion a0 ThetaDeriv) :
    (∀ n : ℕ, nStar ≤ n → ∀ u : ℝ, 0 < u → Zfun u = n → ∀ d : ℕ,
        (d : ℝ) ≤ Dfun a0 (Rfun n u) → Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n : ℕ, nStar ≤ n → ∀ d : ℕ, (d : ℝ) ≤ 6.45 * (n : ℝ) ^ 2 / Real.log n ^ 2 →
        Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 1.1853e-10 * (d : ℝ) ≤ n → Hyperbolic (jensen xiGamma d n)) := by
  have H : Inputs xiGamma ThetaDeriv a0 :=
    ⟨hO, hHad, hGL, zerosSq_xi a0_pos.le hH, hpos, hi⟩
  exact ⟨fun n hn u hu hZ d hd => jensen_hyperbolic_ii a0_pos H hn hu hZ hd,
    fun n hn d hd => jensen_hyperbolic_iii H hn hd,
    fun n d hn hd => jensen_hyperbolic_power_one H hn hd⟩

end

end KimEff

