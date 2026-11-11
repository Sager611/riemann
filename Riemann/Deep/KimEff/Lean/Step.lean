import Riemann.Deep.KimEff.Lean.Core

/-!
# Explicit Kim–Lee wedge, part 2: the named hypotheses and the Obreschkoff step (DERIVATION §3)
-/

open Polynomial Filter Topology

namespace KimEff

noncomputable section

/-! ## Named hypotheses (sources: LEAN.md) -/

/-- (I3) Obreschkoff's sector theorem as Kim–Lee (arXiv 2105.05386) Thm 3 and its corollary:
if the real polynomial P has all zeros in S(δ) then J(P; d) is hyperbolic for d ≤ δ^{-2}
(here also P(0) ≠ 0, which only weakens the statement). -/
def ObreschkoffKimLee : Prop :=
  ∀ δ : ℝ, 0 < δ → ∀ P : ℝ[X], P.eval 0 ≠ 0 →
    (∀ z : ℂ, (P.map (algebraMap ℝ ℂ)).eval z = 0 → z ∈ sector δ) →
    ∀ d : ℕ, (d : ℝ) * δ ^ 2 ≤ 1 → Hyperbolic (jensenOf P d)

/-- (2.2) + §3 (Hadamard, genus 0, and Weierstrass): Θ^{(n)} = F n is the locally uniform
limit of real partial products P_k whose zeros are zeros of F n; in particular
P_k^{(j)}(0) → γ(n+j) for every j. -/
def HadamardApprox (γ : ℕ → ℝ) (F : ℕ → ℂ → ℂ) (n : ℕ) : Prop :=
  ∃ P : ℕ → ℝ[X], (∀ k z, ((P k).map (algebraMap ℝ ℂ)).eval z = 0 → F n z = 0) ∧
    ∀ j, Tendsto (fun k => (derivative^[j] (P k)).eval 0) atTop (𝓝 (γ (n + j)))

/-- (2.3) Gauss–Lucas for genus 0: Z(Θ^{(n)}) lies in the closed convex hull of Z(Θ). -/
def GaussLucasGenus0 (F : ℕ → ℂ → ℂ) (n : ℕ) : Prop :=
  ∀ w, F n w = 0 → w ∈ closure (convexHull ℝ {v | F 0 v = 0})

/-- (2.1) input: every zero of Θ = F 0 is z^2 with |Re z| ≤ a
(proved for Riemann's Θ from (I2) in Xi.lean). -/
def ZerosSq (a : ℝ) (F : ℕ → ℂ → ℂ) : Prop :=
  ∀ w, F 0 w = 0 → ∃ z : ℂ, z ^ 2 = w ∧ |z.re| ≤ a

/-! ## §3: zero-free region ⇒ hyperbolicity -/

/-- DERIVATION §3: if Θ^{(n)} has no zero in W_R then J^{d,n} is hyperbolic for every
d ≤ D(R). -/
theorem hyperbolic_of_zeroFree_W {γ : ℕ → ℝ} {F : ℕ → ℂ → ℂ} {a R : ℝ} {n : ℕ}
    (ha : 0 < a) (hR : 0 < R) (hO : ObreschkoffKimLee) (hH : HadamardApprox γ F n)
    (hG : GaussLucasGenus0 F n) (hZ : ZerosSq a F) (hpos : ∀ m, 0 < γ m)
    (hW : ∀ w ∈ wRegion a R, F n w ≠ 0) {d : ℕ} (hd : (d : ℝ) ≤ Dfun a R) :
    Hyperbolic (jensen γ d n) := by
  have hδ : 0 < sectorDelta a R := sectorDelta_pos ha hR
  have hpar : closure (convexHull ℝ {v | F 0 v = 0}) ⊆ parabola a := by
    refine closure_minimal (convexHull_min ?_ (convex_parabola a)) (isClosed_parabola a)
    intro v hv
    obtain ⟨z, rfl, hz⟩ := hZ v hv
    exact sq_mem_parabola hz
  have hsec : ∀ w, F n w = 0 → w ∈ sector (sectorDelta a R) := by
    intro w hw
    have hwP := hpar (hG w hw)
    have hre : w.re < -R := by
      by_contra h
      exact hW w ⟨hwP, not_lt.mp h⟩ hw
    exact mem_sector_of_parabola ha hR hwP hre
  obtain ⟨P, hPz, hPlim⟩ := hH
  have hdδ : (d : ℝ) * sectorDelta a R ^ 2 ≤ 1 := by
    calc (d : ℝ) * sectorDelta a R ^ 2 ≤ Dfun a R * sectorDelta a R ^ 2 :=
          mul_le_mul_of_nonneg_right hd (sq_nonneg _)
      _ = 1 := Dfun_mul_sq ha hR
  have hev0 : ∀ᶠ k in atTop, (P k).eval 0 ≠ 0 := by
    have h := (hPlim 0).eventually (eventually_ne_nhds (hpos (n + 0)).ne')
    simpa using h
  have hhyp : ∀ᶠ k in atTop, Hyperbolic (jensenOf (P k) d) := by
    filter_upwards [hev0] with k hk
    exact hO _ hδ (P k) hk (fun z hz => hsec z (hPz k z hz)) d hdδ
  refine hyperbolic_of_tendsto (p := fun k => jensenOf (P k) d) (d := d)
    (fun k => natDegree_jensenPoly_le _ _) (natDegree_jensenPoly_le _ _) ?_ ?_ hhyp
  · rw [jensen, coeff_jensenPoly]
    simp [(hpos (n + d)).ne']
  · intro j
    simp only [jensenOf, jensen, coeff_jensenPoly]
    split_ifs
    · exact (hPlim j).const_mul _
    · exact tendsto_const_nhds

end

end KimEff

