/-!
# Round 187 — the matrix-norm conjecture as stated in 2610.30538 fails trivially; other hunts clean

* 2610.30538 (Thompson-type inequalities), TeX l. 730–733: "if (A, B) is a monotone pair in M_n⁺, it is natural to
  conjecture [Bou99] that ‖AZB‖ ≤ √2 ‖ZAB‖ for every symmetric norm", with monotone meaning A = f(C), B = g(C) for
  non-decreasing f, g (l. 494–497) and no condition on Z. Two scouts confirm A = diag(1, 2), B = I (g ≡ 1), Z = E₂₁:
  ‖AZB‖ = 2 > √2 = √2 ‖ZAB‖ in operator and trace norms, and the ratio is unbounded (A = diag(10⁻¹², 1)). So the
  statement as written is false; most likely a hypothesis of Bourin's original (e.g. Z ≥ 0) was dropped. To check.
* No counterexample: 2610.30493 (1,176 Gaussian water-filling cases), 2610.30800 (finite X), 2610.30164 (k^n, M_n(Q)),
  2610.25661 (|n! − m^k| > n for n ≤ 500), 2610.31557 (γ(F₂(K_{a,b,c})) for a + b + c ≤ 9). Gale R⁴/R⁵ QA timed out.
-/

namespace Fayed.R187
end Fayed.R187
