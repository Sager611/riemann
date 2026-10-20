import Mathlib

/-!
# Round 9 — gap map: RH-equivalents vs best unconditional results

| id | RH-equivalent | best unconditional | gap / status |
|----|---------------|--------------------|--------------|
| R9.1 | M(x) = O(x^{1/2+ε}) | M(x) ≪ x exp(-c (log x)^{3/5} (log log x)^{-1/5}) (Korobov–Vinogradov) | exponent 1/2: DEAD |
| R9.2 | |ψ(x) - x| ≤ √x log²x / 8π (Schoenfeld) | same shape as R9.1, explicit (Fiori–Kadiri–Swidinsky 2023) | constants only, not exponent |
| R9.3 | Re ξ'/ξ > 0 on Re s > 1/2 (Lagarias 1999) | Re s ≥ 1 | ≡ quasi-RH strip: open |
| R9.4 | ζ' ≠ 0 on 0 < Re s < 1/2 (Speiser, Levinson–Montgomery) | counts only (Soundararajan 1998, Zhang 2001) | needs exact 0: DEAD |
| R9.5 | Weil positivity on dense classes | windows ≤ 1.6 (Zhu 2026); RH to 3·10¹² | see R3.1 deep dive |
| R9.6 | Λ_DBN ≤ 0 | Λ ≤ 0.2 (Platt–Trudgian 2021), Λ ≥ 0 | 0.1 reachable by computation, 0 not |
| R9.7 | Báez-Duarte d_N → 0 | d_N ≤ d_1 ≈ 0.561 | no o(1): open |
| R9.8 | density hypothesis + simplicity ⇒ RH? | Huxley exponent; Maynard–Pratt 2024 | no known implication: DEAD |
-/

open Real ArithmeticFunction

namespace Fayed.R9

/-- Chebyshev ψ(x) = Σ_{n ≤ x} Λ(n). -/
noncomputable def psi (x : ℝ) : ℝ := ∑ n ∈ Finset.Icc 1 ⌊x⌋₊, vonMangoldt n

/-- R9.2: Schoenfeld's form of RH (≡ RH, cited). -/
def SchoenfeldBound : Prop := ∀ x : ℝ, 73 ≤ x → |psi x - x| ≤ (8 * π)⁻¹ * √x * (log x) ^ 2

/-- R9.3: Riemann's ξ, entire: ξ(s) = s(s-1)/2 · Λ₀(s) + 1/2. -/
noncomputable def xi (s : ℂ) : ℂ := s * (s - 1) / 2 * completedRiemannZeta₀ s + 1 / 2

/-- R9.3: Lagarias' positivity (≡ RH, cited): Re ξ'(s)/ξ(s) > 0 at every non-zero of ξ in
Re s > 1/2. The hypothesis `xi s ≠ 0` keeps Lean's junk value ξ'/0 = 0 out of the statement. -/
def LagariasPositivity : Prop :=
  ∀ s : ℂ, 1 / 2 < s.re → xi s ≠ 0 → 0 < (deriv xi s / xi s).re

/-- R9.4: Speiser's criterion (≡ RH, cited). -/
def Speiser : Prop := ∀ s : ℂ, 0 < s.re → s.re < 1 / 2 → deriv riemannZeta s ≠ 0

end Fayed.R9
