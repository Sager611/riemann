import Riemann.Deep.ExplorePointProc.Lean.Cert
import Riemann.Deep.ExplorePointProc.Lean.Liminf

/-!
# Final forms with Lamzouri's D(λ) = liminf (TeX l.183)
-/

namespace PPLean

/-- **Theorem 1** (DERIVATION §4) with Lamzouri's D: 0 < D(λ). -/
theorem theorem1_D {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z)
    (h12 : Steps12 Z C) (hFC : FourierCc2 C) (hT : LemmaT' C)
    {θ lam U β c : ℝ} {f A : ℝ → ℝ}
    (hθ1 : 1 / 2 < θ) (hθ2 : θ < 4 / 7) (hlam0 : 0 < lam) (hlam1 : lam ≤ 1)
    (hU1 : lam < U) (hU3 : U ≤ 3 * lam)
    (hfS : C.Sch f) (hfe : ∀ x, f (-x) = f x) (hfn : ∀ x, 0 ≤ f x)
    (hfneg : ∀ u, U ≤ u → C.hat f u ≤ 0)
    (hA : ∃ K : ℝ, ∀ y, lam ≤ y → |A y| ≤ K) (hβ : 0 ≤ β) (hH : HH lam β c (C.hat f) A)
    (hJ : C.hat f 0 + 2 * c < C.J θ f) : 0 < Dlam Z lam :=
  (dpos_iff (nToInf hsc hRvM) lam).1
    (theorem1 hZ hsc hRvM h12 hFC hT hθ1 hθ2 hlam0 hlam1 hU1 hU3 hfS hfe hfn hfneg hA hβ hH hJ)

/-- **Corollary** (DERIVATION §5): D(λ) > 0 for every λ ≥ 0.768, Lamzouri's D, exact certificate form. -/
theorem D_pos_ge_768 {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z)
    (h12 : Steps12 Z C) (hFC : FourierCc2 C) (hT : LemmaT' C) {f : ℝ → ℝ}
    (hS : C.Sch f) (he : ∀ x, f (-x) = f x) (hrep : HermRep f (C.hat f)) (hP : CertP) (hQ : CertQ)
    (hJ : CertJ C f) (hR1 : CertR1 C f Acert) (hR2 : CertR2 C f Acert) :
    ∀ lam : ℝ, 96 / 125 ≤ lam → 0 < Dlam Z lam := fun lam hlam =>
  (dpos_iff (nToInf hsc hRvM) lam).1
    (dpos_ge_768_exact hZ hsc hRvM h12 hFC hT hS he hrep hP hQ hJ hR1 hR2 lam hlam)

end PPLean
