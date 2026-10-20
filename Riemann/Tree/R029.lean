import Riemann.Deep.Backbone
import Riemann.Deep.Lagarias.Suff

/-!
# Round 29 — program blueprints: the open core of each program

| program | open core | chain in Lean |
|---------|-----------|---------------|
| Weil positivity | ∀ L, WeilPos L | criterion cited (Riemann/Deep/Weil/Defs.lean) |
| Connes trace formula | positivity of the trace form ≡ Weil | cited |
| Deninger | arithmetic H¹ + self-adjoint A with det = ξ | conjectural objects |
| de Branges | Hermite–Biehler E with A_E = Ξ | cited |
| Nyman–Beurling / Báez-Duarte | d_N → 0 | sufficiency machine-checked (open_cores_imply_RH) |
| de Bruijn–Newman | Λ ≤ 0 | cited (Rodgers–Tao Λ ≥ 0) |
| Li / Keiper | λ_n ≥ 0 ∀ n | cited |
| automorphic GRH | GRH for all cuspidal π | GL₁ case = RH |
| Mertens / Redheffer | M(x) = O(x^{1/2+ε}) | sufficiency machine-checked |
| Lagarias | Re ξ'/ξ > 0 on Re s > 1/2 | sufficiency machine-checked |
-/

namespace Fayed.R29

/-- Any one of the machine-checked open cores implies RH. -/
theorem open_cores_imply_RH
    (h : Fayed.NB.NBApprox ∨ Fayed.Backbone.BDApprox ∨ Fayed.R9.LagariasPositivity ∨
      (∀ ε > 0, Fayed.R2.MertensBound (1 / 2 + ε))) : RiemannHypothesis := by
  rcases h with h | h | h | h
  · exact Fayed.NB.RH_of_NBApprox h
  · exact Fayed.Backbone.RH_of_BDApprox h
  · exact Fayed.Lagarias.RH_of_LagariasPositivity h
  · exact Fayed.Backbone.RH_of_mertens h

end Fayed.R29

#print axioms Fayed.R29.open_cores_imply_RH
