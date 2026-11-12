/-!
# Round 446 — node N12 (Selberg class and degree-two analogues), round 15, state of the art around today's ceiling of the one-piece Levinson–Conrey method for GL(2): no other method beats our 7.0056% (Rezvyakova's proportion is unspecified); no mollified second moment for a fixed GL(2) L-function beyond Andersen–Thorner's ν < 25/128, and Kim–Sarnak's 7/64 is still the best exponent, so the length cannot move now; the only multi-piece GL(2) mollifier, Kühn–Robles–Zeindler 1609.03738, gives 2.97607% unconditionally and 6.93872% under Ramanujan at shorter lengths and without optimisation; Bernard's admissible class is exactly the one the ceiling covers; network flows, Markov chain Monte Carlo and Szegő/Fisher–Hartwig transfer no tool

* Node: N12 (15th round; R439 before). Round type: [check] the escape routes from a fresh barrier of ours
  (discoveries/gl2-levinson-ceiling.md: the one-piece Levinson–Conrey functional for GL(2) is capped at 0.07005633 for
  ν < 25/128, 0.08971026 under Ramanujan; our record 0.0700562390).
* [known] Other methods: no proportion above 0.0700562 for every GL(2) newform; Andersen–Thorner 2004.03581 (Thm
  main_theorem2: > 0.0696, > 0.0896 under Ramanujan) was the best; Rezvyakova 1212.2948 gives an unspecified c > 0: HOLDS.
* [check] Length: no 2024-2026 mollified second moment for a fixed GL(2) L-function beyond ν < 1/4 − θ/2 = 25/128; the closest,
  2601.15268, is a half-integral-weight family result with mollifier scale X^ε (TeX l.153-156, 301-307): HOLDS.
* [check] θ_KS: 7/64 (Kim–Sarnak) is still the best exponent over Q (2402.11761, introduction: λ_f ≥ 1/4 − (7/64)²);
  2011.09054 is withdrawn: HOLDS.
* [known] Multi-piece mollifiers for GL(2): Kühn–Robles–Zeindler 1609.03738 (On mean values of mollifiers and L-functions
  associated to primitive cusp forms), Thm theoremofkappa (TeX l.293-299): ψ = ψ₁ + ψ₂ (l.202-206), κ_f ≥ 2.97607%
  unconditionally and 6.93872% under Ramanujan, at θ = 7/64, ν₁ = 5/27, ν₂ = 25/149, R = 3.21 (l.1765): KNOWN, outside our
  one-piece ceiling, below our record and unoptimised (Bernard's one-piece 2.97% at 5/27 became 6.64% with optimisation).
* [falsify] Admissible class (Bernard 1404.6432 l.161, l.770, l.785-795): real P with P(0) = 0, P(1) = 1, R > 0, Q(0) = 1 in
  Bernard's form, zero bound with c(P, Q, R, ν/2); complex Q appears only in the mean value theorem: exactly the class the
  ceiling covers: HOLDS.
* Exploration (tool transfers, tested): Hoffman's circulation theorem (no finite network for Weil's form without re-encoding
  Q(f)), coupling from the past (no canonical kernel on zero ordinates; real parts unconstrained), strong Szegő/Fisher–Hartwig
  (leads only to the Toeplitz positivity of Li coefficients, already false at 2 × 2, R395): no transfer.
-/

namespace Fayed.R446
end Fayed.R446
