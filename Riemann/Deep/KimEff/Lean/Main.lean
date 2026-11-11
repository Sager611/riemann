import Riemann.Deep.KimEff.Lean.Power
import Riemann.Deep.KimEff.Lean.Half
import Riemann.Deep.KimEff.Lean.NStar
import Riemann.Deep.KimEff.Lean.Real
import Riemann.Deep.KimEff.Lean.TwoThirds

/-!
# Explicit Kim–Lee wedge: main theorem (DERIVATION §0)

kimEff_main: for Riemann's ξ, a = 0.37338196661320752, n_* = 998,616,377,771, the named inputs
(Bebel's half-plane, Obreschkoff–Kim–Lee, Hadamard partial products, Gauss–Lucas for genus 0,
γ > 0, and the zero-free region (i) of DERIVATION §4–§9) imply (ii) and (iii) with all four
power forms. xi_jensen_hyperbolic_half: (ii), (iii) with a = 1/2 and no zero-location hypothesis.
Also imported: Zfun_12_bounds, twelve_le_u_iff (n_* = ceil Z(12), u(n) ≥ 12 iff n ≥ n_*) and
xiGammaC_eq, xiF_eq_tsum (γ is the real Taylor coefficient of 8ξ(1/2+z)), and
jensen_hyperbolic_twoThirds (the previous two-thirds wedge, for d > 1.04710e8, lies inside).
-/

namespace KimEff

noncomputable section

/-- DERIVATION §0, (ii) and (iii), for Riemann's ξ and a = 0.37338196661320752. -/
theorem kimEff_main (hH : HalfPlane a0) (hO : ObreschkoffKimLee)
    (hHad : ∀ n, HadamardApprox xiGamma ThetaDeriv n)
    (hGL : ∀ n, GaussLucasGenus0 ThetaDeriv n) (hpos : ∀ m, 0 < xiGamma m)
    (hi : ZeroFreeRegion a0 ThetaDeriv) :
    (∀ n : ℕ, nStar ≤ n → ∀ u : ℝ, 0 < u → Zfun u = n → ∀ d : ℕ,
        (d : ℝ) ≤ Dfun a0 (Rfun n u) → Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n : ℕ, nStar ≤ n → ∀ d : ℕ, (d : ℝ) ≤ 6.45 * (n : ℝ) ^ 2 / Real.log n ^ 2 →
        Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 7.6890 * (d : ℝ) ^ (0.51 : ℝ) ≤ n →
        Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 0.87160 * (d : ℝ) ^ (0.55 : ℝ) ≤ n →
        Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 0.069829 * (d : ℝ) ^ (0.6 : ℝ) ≤ n →
        Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n d : ℕ, nStar ≤ n → 1.1853e-10 * (d : ℝ) ≤ n → Hyperbolic (jensen xiGamma d n)) := by
  have H : Inputs xiGamma ThetaDeriv a0 :=
    ⟨hO, hHad, hGL, zerosSq_xi a0_pos.le hH, hpos, hi⟩
  exact ⟨fun n hn u hu hZ d hd => jensen_hyperbolic_ii a0_pos H hn hu hZ hd,
    fun n hn d hd => jensen_hyperbolic_iii H hn hd,
    fun n d hn hd => jensen_hyperbolic_power_051 H hn hd,
    fun n d hn hd => jensen_hyperbolic_power_055 H hn hd,
    fun n d hn hd => jensen_hyperbolic_power_06 H hn hd,
    fun n d hn hd => jensen_hyperbolic_power_one H hn hd⟩

end

end KimEff

