/-!
# Round 437 — node N10 (explicit ψ and π − li bounds), round 13, falsification of the technical lemmas of Chirre–Helfgott 2511.14736 under our MertensCH record: Lemmas arborio, barbar and artanor (constants), zetinvbound, paulin, pommedupe, benvenuto and Prop. shiftbounded hold at random and corner parameters; arborio's 1/(16πu²) and zetinvbound's first bound are sharp in the limit (0.99992 at λ = −0.001, 0.99958 at s = −200); two harmless misprints; both papers still v1, no pre-emption; finite elements, S-matrix dispersion relations and tidal analysis transfer no tool

* Node: N10 explicit ψ and π − li bounds (13th round; R428 before). Round type: falsify the unrefereed technical lemmas of our record's input
  (chaudit's open item 7; TeX /tmp copy of the 2511.14736v1 e-print). The four [falsify] scouts timed out; the orchestrator ran each check.
* [falsify] Lemma arborio (l. 5513–5521): the proof's representation of φ̂_λ − I_λ (Lemma 5 of zbMATH06384942, eq. monthmayending)
  agrees with the direct Fourier integral to 1.2·10⁻¹¹ (8 λ, ~700 points); on |u| ∈ [1/2, 200] (step 0.004), 24 λ from ±10⁻⁶ to ±40:
  max 16πu²|f| = 0.999918 (λ = −0.001, u = −199.752; sharp as λ → 0 at quarter-integers), max (22/3)u²|f′| = 0.91666 (limit 11/12);
  for 0 < |u| < 1/2 the sign and monotonicity claims hold (max u·sgn(λ)·f = −1·10⁻¹³, min sgn(λ)·f′ = 2·10⁻⁹), jump −sgn(λ), f(±1/2) = 0:
  HOLDS (arb2's outlier at λ = ±10⁻⁴ was a 500-node quadrature artefact). Misprint: the proof's chain h(u) ≤ 2 sin πu/(πu) ≤ 4/(3π)
  on [1, 3/2] is false as written, but max h = 4/(3π) (u = 3/2) and max κ = 1 hold there, and κ ≤ 1.055 < 12/11 on [3/2, 50].
* [falsify] Lemma barbar and artanor's constants (l. 5602–5643): Ei(6.4) − e^{6.4}·5.4/(6.4·4.4) = 2.524241 ≤ 2.525; c₁(1/2) = 3.396086 < 17/5,
  c₀(1/2) = 0.754603 < 16/21, c₂(1/2) = 0.212816 < 2/9, c₁(1/6) = 2.835678 ≤ 2.85, c₀(1/6) = 0.637721 ≤ 2/3: HOLDS.
* [falsify] Lemma zetinvbound (l. 6138–6150) by the exact functional-equation reduction (equal to direct ζ evaluation to 10⁻²⁹ in log):
  20,000 random s and 182 corners: claim (i) max ratio 0.99958 (s = −200 + 10⁻⁶i: Stirling's lower bound for |Γ(1 − s)|, sharp on the
  negative real axis), claim (ii) max 0.384 (s = i): HOLDS.
* [falsify] Prop. shiftbounded (l. 2697–2705) end to end with F = 1/ζ, T = 20, 30, x = 50, 10³, 10⁴, λ = ±0.5, 2, −2, −2π/T (residues at
  the zeros, 40 trivial zeros and s* = 1 + λT/2π): max |LHS − ΣRes|/bound = 0.864 (T = 30, x = 10³, λ = −2π/T): HOLDS (the scout's
  ratio 1.75 used only the ξ = +1 half of the bound).
* [falsify] Lemmas paulin, pommedupe, benvenuto (l. 3061–3121, 2160): |w_{δ,σ}(t)| ≤ 1/(δ(σ − t)) + 2 on 3,000 random cases; pommedupe's
  trivial-zero sum stays below its right side on 424 cases, max 0.9999996 (σ → −2, δ → 0, x = 2: sharp), its "=" is a misprint for
  O*(…) (the proof bounds it; Cor. jolene repeats it); |Φ_λ(±1 + it)| = |tanh(πt/2 + λ/4) − tanh(λ/4)|/4 ≤ π|t|/8: HOLDS.
* [check] 2511.14736 and 2512.15709 still v1; OpenAlex adds only Zuniga-Alterman 2610.17720 (no error reported); no newer explicit
  |M(x)| ≤ c·x^θ, θ < 1, beating mertensch: HOLDS. Open (unchanged): Platt's private sums above 1.9·10⁶, Props. gould/qmediano.
* Exploration (tool transfers, tested): finite-element Rayleigh–Ritz gives only upper bounds on the Weil form's minimum; twice-subtracted
  dispersion positivity of −log(ξ(1/2 + iz)/ξ(1/2)) (c₂ = Σ1/γ² = 0.02310499, c₄ = 1.85863·10⁻⁵) is an RH consequence blind to zeros off
  the line; the tidal Rayleigh resolution criterion is weaker than known gap bounds (R408): no transfer.
-/

namespace Fayed.R437
end Fayed.R437
