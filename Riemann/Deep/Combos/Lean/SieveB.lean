module

public import Riemann.Deep.Combos.Lean.SieveA

@[expose] public section

/-! # Lemma 2, steps: decomposition 𝔞 = m𝔟 (l.80–81), π_C ≤ P_Q(N) and the constant (l.85–88),
|B_𝔟|² ≤ ζ-bound · Σ_m weight·|c_{m𝔟}|² (l.88–89) -/

open Finset ComplexConjugate

namespace RingClass

variable {G ι : Type*} [AddCommGroup G]

/-- B_𝔟 = Σ_m c_{m𝔟} (NOTE l.81) -/
def Bco [DecidableEq ι] (d : RCData G ι) (S : Finset ι) (c : ι → ℂ) (B : ι) : ℂ :=
  ∑ A ∈ S with d.prim A = B, c A

/-- the primitive 𝔟 (of support elements) in the class C; π_C counts them -/
def Pcl [DecidableEq ι] [DecidableEq G] (d : RCData G ι) (S : Finset ι) (g : G) : Finset ι :=
  (S.image d.prim).filter (fun B => d.cls B = g)

/-- v(C) = Σ_{𝔟 ∈ C} B_𝔟 -/
def vfun [DecidableEq ι] [DecidableEq G] (d : RCData G ι) (S : Finset ι) (c : ι → ℂ) (g : G) : ℂ :=
  ∑ B ∈ Pcl d S g, Bco d S c B

/-- Σ_𝔞 c_𝔞 χ(𝔞) = Σ_𝔟 χ(𝔟) B_𝔟 = Σ_C χ(C) v(C) (NOTE l.80–81) -/
lemma sum_decomp [Fintype G] [DecidableEq G] [DecidableEq ι] (d : RCData G ι) (S : Finset ι)
    (c : ι → ℂ) (χ : AddChar G ℂ) :
    ∑ A ∈ S, c A * χ (d.cls A) = ∑ g, χ g * vfun d S c g := by
  rw [← Finset.sum_fiberwise_of_maps_to (g := d.prim) (t := S.image d.prim)
    (fun A hA => Finset.mem_image_of_mem _ hA)]
  have h1 : ∀ B ∈ S.image d.prim, ∑ A ∈ S with d.prim A = B, c A * χ (d.cls A)
      = χ (d.cls B) * Bco d S c B := by
    intro B _
    rw [Bco, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun A hA => ?_)
    rw [Finset.mem_filter] at hA
    rw [← hA.2, d.cls_prim]; ring
  rw [Finset.sum_congr rfl h1, ← Finset.sum_fiberwise (S.image d.prim) d.cls]
  refine Finset.sum_congr rfl (fun g _ => ?_)
  rw [vfun, Pcl, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun B hB => ?_)
  rw [Finset.mem_filter] at hB
  rw [hB.2]

/-- h π_C ≤ 30 L² (N + |D|^{1/2}) by (F1), (F2), Lemma 1 and (F3) h ≤ 3|D|^{1/2}L²
(NOTE l.85–88; L plays log(3|D|)) -/
lemma class_count [Fintype G] [DecidableEq G] [DecidableEq ι] (d : RCData G ι) (hF12 : F12 d)
    (L : ℝ) (hF3 : (Fintype.card G : ℝ) ≤ 3 * Real.sqrt |(d.D : ℝ)| * L ^ 2)
    (N : ℝ) (hN : 1 ≤ N) (S : Finset ι) (hS : ∀ A ∈ S, (d.nrm A : ℝ) ≤ N) (g : G) :
    (Fintype.card G : ℝ) * ((Pcl d S g).card : ℝ) ≤ 30 * L ^ 2 * (N + Real.sqrt |(d.D : ℝ)|) := by
  obtain ⟨a, b, c', ideal, ha, hb, hac, hD, hrep⟩ := hF12 g
  have hb' : |(-b)| ≤ a := by rwa [abs_neg]
  have hsub : (↑(Pcl d S g) : Set ι) ⊆ ideal '' PQset a (-b) c' N := by
    intro B hB
    rw [Finset.mem_coe, Pcl, Finset.mem_filter, Finset.mem_image] at hB
    obtain ⟨⟨A, hA, rfl⟩, hcl⟩ := hB
    obtain ⟨x, y, hxy, hQ, hid⟩ := hrep (d.prim A) (d.cont_prim A) hcl
    refine ⟨(x, y), ⟨hxy, ?_⟩, hid⟩
    have h2 : d.nrm (d.prim A) ≤ d.nrm A := by
      rw [d.nrm_eq A]; exact Nat.le_mul_of_pos_left _ (pow_pos (d.cont_pos A) 2)
    have h3 : (d.nrm (d.prim A) : ℝ) ≤ d.nrm A := by exact_mod_cast h2
    have h4 : ((Qf a (-b) c' x y : ℤ) : ℝ) = (d.nrm (d.prim A) : ℝ) := by exact_mod_cast hQ
    show ((Qf a (-b) c' x y : ℤ) : ℝ) ≤ N
    rw [h4]; exact h3.trans (hS A hA)
  have hfin : (PQset a (-b) c' N).Finite := PQset_finite a (-b) c' ha hb' hac N (by linarith)
  have hPg : ((Pcl d S g).card : ℝ) ≤ PQ a (-b) c' N := by
    have h1 := Set.ncard_le_ncard hsub (hfin.image ideal)
    rw [Set.ncard_coe_finset] at h1
    have h2 := Set.ncard_image_le (f := ideal) hfin
    exact_mod_cast h1.trans h2
  have hL1 := lemma1 a (-b) c' ha hb' hac N hN
  have hDD : ((-b : ℤ) : ℝ) ^ 2 - 4 * a * c' = (d.D : ℝ) := by
    rw [← hD]; push_cast; ring
  rw [hDD] at hL1
  have hDpos : (0 : ℝ) < |(d.D : ℝ)| := abs_pos.mpr (by exact_mod_cast d.D_neg.ne)
  set u := Real.sqrt (Real.sqrt |(d.D : ℝ)|) with hu
  have hu0 : 0 < u := Real.sqrt_pos.mpr (Real.sqrt_pos.mpr hDpos)
  have hu2 : u ^ 2 = Real.sqrt |(d.D : ℝ)| := Real.sq_sqrt (Real.sqrt_nonneg _)
  have hx2 : Real.sqrt N ^ 2 = N := Real.sq_sqrt (by linarith)
  rw [← hu2] at hL1 hF3 ⊢
  exact sieve_arith L _ _ N (Real.sqrt N) u hu0 hx2 hF3 (hPg.trans hL1) (Nat.cast_nonneg _)

/-- |B_𝔟|² ≤ Z₀ Σ_m W(m)|c_{m𝔟}|² when Σ_{m ∈ K} W(m)⁻¹ ≤ Z₀ for all finite K ⊂ ℕ_{>0}
(NOTE l.88–89 with W(m) = m^{1+2ε}, Z₀ = ζ(1+2ε)) -/
lemma fiber_bound [DecidableEq ι] (d : RCData G ι) (S : Finset ι) (c : ι → ℂ) (W : ℕ → ℝ)
    (hW : ∀ m, 0 < m → 0 < W m) (Z₀ : ℝ)
    (hZ : ∀ K : Finset ℕ, (∀ k ∈ K, 0 < k) → ∑ k ∈ K, (W k)⁻¹ ≤ Z₀) (B : ι) :
    Complex.normSq (Bco d S c B) ≤
      Z₀ * ∑ A ∈ S with d.prim A = B, W (d.cont A) * Complex.normSq (c A) := by
  have h1 := normSq_sum_le (S.filter (fun A => d.prim A = B)) (fun A => W (d.cont A))
    (fun A _ => hW _ (d.cont_pos A)) c
  have hinj : Set.InjOn d.cont (S.filter (fun A => d.prim A = B)) := by
    intro A hA A' hA' h
    rw [Finset.mem_coe, Finset.mem_filter] at hA hA'
    exact d.decomp_inj A A' (hA.2.trans hA'.2.symm) h
  have e : ∑ k ∈ (S.filter (fun A => d.prim A = B)).image d.cont, (W k)⁻¹
      = ∑ A ∈ S with d.prim A = B, (W (d.cont A))⁻¹ := Finset.sum_image hinj
  have hz : ∑ A ∈ S with d.prim A = B, (W (d.cont A))⁻¹ ≤ Z₀ := by
    rw [← e]
    apply hZ
    intro k hk
    rw [Finset.mem_image] at hk
    obtain ⟨A, _, rfl⟩ := hk
    exact d.cont_pos A
  have h0 : 0 ≤ ∑ A ∈ S with d.prim A = B, W (d.cont A) * Complex.normSq (c A) :=
    Finset.sum_nonneg (fun A _ => mul_nonneg (hW _ (d.cont_pos A)).le (Complex.normSq_nonneg _))
  exact h1.trans (mul_le_mul_of_nonneg_right hz h0)

end RingClass

