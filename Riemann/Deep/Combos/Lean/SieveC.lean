import Riemann.Deep.Combos.Lean.SieveB

/-! # Lemma 2 (NOTE l.76–89), core form: weight W and ζ-type constant Z₀ abstract, L = log(3|D|)
abstract, |D|^{1/2} = √|D|; Lemma2.lean specialises W(m) = m^{1+2ε}, Z₀ = 1/ε, L = log(3|D|). -/

open Finset ComplexConjugate

namespace RingClass

variable {G ι : Type*} [AddCommGroup G]

/-- Lemma 2, core form, for every finite set 𝒞 of characters (𝒞 = Ĝ in the NOTE). -/
theorem sieve_core [Fintype G] [DecidableEq G] [DecidableEq ι] (d : RCData G ι) (hF12 : F12 d)
    (L : ℝ) (hF3 : (Fintype.card G : ℝ) ≤ 3 * Real.sqrt |(d.D : ℝ)| * L ^ 2)
    (W : ℕ → ℝ) (hW : ∀ m, 0 < m → 0 < W m) (Z₀ : ℝ)
    (hZ : ∀ K : Finset ℕ, (∀ k ∈ K, 0 < k) → ∑ k ∈ K, (W k)⁻¹ ≤ Z₀)
    (𝒞 : Finset (AddChar G ℂ)) (N : ℝ) (hN : 1 ≤ N) (S : Finset ι)
    (hS : ∀ A ∈ S, (d.nrm A : ℝ) ≤ N) (c : ι → ℂ) :
    ∑ χ ∈ 𝒞, Complex.normSq (∑ A ∈ S, c A * χ (d.cls A)) ≤
      Z₀ * (30 * L ^ 2 * (N + Real.sqrt |(d.D : ℝ)|)) *
        ∑ A ∈ S, W (d.cont A) * Complex.normSq (c A) := by
  set K := 30 * L ^ 2 * (N + Real.sqrt |(d.D : ℝ)|) with hK
  have hK0 : 0 ≤ K := by
    have : 0 ≤ N + Real.sqrt |(d.D : ℝ)| := by have := Real.sqrt_nonneg |(d.D : ℝ)|; linarith
    positivity
  calc ∑ χ ∈ 𝒞, Complex.normSq (∑ A ∈ S, c A * χ (d.cls A))
      = ∑ χ ∈ 𝒞, Complex.normSq (∑ g, χ g * vfun d S c g) :=
        Finset.sum_congr rfl (fun χ _ => by rw [sum_decomp])
    _ ≤ (Fintype.card G : ℝ) * ∑ g, Complex.normSq (vfun d S c g) := bessel 𝒞 _
    _ ≤ (Fintype.card G : ℝ) * ∑ g, (((Pcl d S g).card : ℝ) *
          ∑ B ∈ Pcl d S g, Complex.normSq (Bco d S c B)) := by
        apply mul_le_mul_of_nonneg_left _ (Nat.cast_nonneg _)
        exact Finset.sum_le_sum (fun g _ => normSq_sum_le_card _ _)
    _ = ∑ g, ((Fintype.card G : ℝ) * ((Pcl d S g).card : ℝ)) *
          ∑ B ∈ Pcl d S g, Complex.normSq (Bco d S c B) := by
        rw [Finset.mul_sum]; exact Finset.sum_congr rfl (fun g _ => by ring)
    _ ≤ ∑ g, K * ∑ B ∈ Pcl d S g, Complex.normSq (Bco d S c B) := by
        apply Finset.sum_le_sum
        intro g _
        exact mul_le_mul_of_nonneg_right (class_count d hF12 L hF3 N hN S hS g)
          (Finset.sum_nonneg (fun B _ => Complex.normSq_nonneg _))
    _ = K * ∑ B ∈ S.image d.prim, Complex.normSq (Bco d S c B) := by
        rw [← Finset.mul_sum]
        congr 1
        exact Finset.sum_fiberwise (S.image d.prim) d.cls _
    _ ≤ K * ∑ B ∈ S.image d.prim,
          Z₀ * ∑ A ∈ S with d.prim A = B, W (d.cont A) * Complex.normSq (c A) := by
        apply mul_le_mul_of_nonneg_left _ hK0
        exact Finset.sum_le_sum (fun B _ => fiber_bound d S c W hW Z₀ hZ B)
    _ = Z₀ * K * ∑ A ∈ S, W (d.cont A) * Complex.normSq (c A) := by
        rw [← Finset.mul_sum, Finset.sum_fiberwise_of_maps_to
          (fun A hA => Finset.mem_image_of_mem _ hA)]
        ring

end RingClass

