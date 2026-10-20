/-!
# Round 54 — adversarial self-audit of our Lean statements

| file | statement | finding |
|------|-----------|---------|
| Riemann/Fayed/Theorem1.lean | Theorem1, Theorem1C, theorem1_false, theorem1C_false | none |
| Riemann/Deep/Backbone.lean | MertensBound, BDApprox, RH_of_* | none |
| Riemann/Deep/NymanBeurling/Suff.lean | rhoNB, NBApprox | none |
| Riemann/Tree/R009.lean + Riemann/Deep/Lagarias/Suff.lean | LagariasPositivity, RH_of_LagariasPositivity | SERIOUS: at a zero of ξ, Lean's deriv ξ / ξ = deriv ξ / 0 = 0, so "0 < Re(ξ'/ξ)" already excludes zeros by junk value; the proof uses exactly that. Fix: quantify only over ξ s ≠ 0 and prove via the pole m/(s-ρ) of ξ'/ξ. |
| Riemann/Deep/Riesz/Suff.lean | RH_of_rieszBound | none |
| Riemann/Fayed/Paper.lean | fayed_RH : Inputs → Theorem1C → RH | SERIOUS (fidelity): Inputs is never inhabited, so the reduction could be vacuous if a field is misstated; strip, reflect, convexity are proved, em17, perron33, leibniz44 are not. Fix: prove them (queued). |

Lesson: every RH criterion written with a quotient or log must exclude the zeros explicitly, or the junk
value 0 silently turns the criterion into "no zeros".
-/
