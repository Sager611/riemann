/-!
# Round 442 — node N8 (critical-line zeros: gaps), round 15, falsification of the small-gap method: Lamzouri's criterion (arXiv 2611.14486 Thm Main2, TeX l. 392–400: 1/2 < θ < 4/7, f ∈ C_c² even, f ≥ 0, f^(0) = 1, f^ ≤ 0 on |u| > λ, J_θ(f) = f(0) + 2∫_0^θ f(α)α dα > 1 ⇒ D(λ) > 0) is a linear program whose duals, pseudo pair-correlation lattice measures with a hard core λ and structure factor ≥ |α| on [−θ, θ], certify (Arb) that no admissible f exists for λ ≤ 0.775 (every θ ≤ 4/7) nor, in the RH case θ = 1, for λ ≤ 0.55 (explore34; the 1/(1+θ)-discrete sine DPP alone gives 1/(1+θ)); our record 0.7851 is within 0.0101 of the floor; no prior floor for this LP (Chirre–Gonçalves–de Laat, Lagarias–Rodgers' 1/2, Carneiro–Chandee–Chirre–Milinovich); graph limits, Remez inequalities and Azuma–Hoeffding transfer no tool

* Node: N8 (15th round; R441 before). Round type: falsify the method behind our own record (small gaps: D(λ) > 0 for λ ≥ 0.7851,
  smallgap, 11 Nov 2026) by its LP dual; the listing of 12 Nov was 2 h late and goes to R443.
* [falsify] explore34 (Orca, Riemann/Deep/ExploreDPP): Lemma 1 (weak duality): a measure σ = Σ_{k≠0} h(1 + η_|k|)δ_{kh} ≥ 0 with no
  mass in (−λ, λ) whose structure factor P(α) = (1 − h) + 2hΣ η_k cos(2παkh) is ≥ |α|1[|α| ≤ θ] on [0, 1/(2h)] gives J_θ(f) ≤ f^(0)
  for every admissible f; Theorem A: the 1/(1+θ)-discrete sine (η_k = −sinc²(ka)) has P = 1 − Σ_m Δ(α − m/a) = |α| on [−θ, θ];
  Theorem B: Arb certificates (rational η, cells with |P''| ≤ M2) for θ = 4/7, λ = 0.775 (min G ≥ 0.000204, 90,324 cells; rerun by the
  orchestrator: PASS, identical, 80 MB) and θ = 1, λ = 0.55 (min ≥ 0.0322). Reviewer: core CORRECT, SIGNIFICANT as a method-specific
  barrier; text repairs requested (θ ≤ 1/(2h) at θ = 1, a false intermediate inequality, interval upper end 0.7851, not 0.7849).
* [check] Lamzouri 2611.14486: Thm Main2 as stated above (l. 392–400), D(λ) = liminf N(T)⁻¹ #{γ_{n+1} − γ_n ≤ 2πλ/log T} (l. 183);
  his own best λ > 0.7864, 'further optimization may improve' (l. 195–199): HOLDS.
* [check] novelty: Chirre–Gonçalves–de Laat 1810.08843 Thm 4 (l. 317–325): 0.6039 under RH and N*/N → 1, no dual floor; Lagarias–
  Rodgers 1905.12123: band-limited correlations are compatible with gaps ≥ 1/2 (the AH obstruction = Theorem A at θ = 1);
  Carneiro–Chandee–Chirre–Milinovich 2108.09258: no small-gap limit above 1/2; 2604.05733 (μ < 0.50895 under RH by a
  resonance-correlation method) quotes the Montgomery–Odlyzko method's limit 0.508 for μ (Inoue–Kobayashi–Toma 2510.14309), a
  different problem; a scout's KNOWN on it was a false positive. No prior floor for the D(λ) LP.
* [falsify] Theorem A's identity numerically (K = 200,000, 20,001 points, θ = 4/7 and 1): max error 1.0·10⁻⁶ (the truncation tail): HOLDS.
* Exploration (tool transfers, tested): Frieze–Kannan weak regularity for graphons (no bounded graphon preserves Weil's form for all
  test functions), the Remez inequality (Taylor truncations of Ξ carry uncontrolled remainders), Azuma–Hoeffding for Li coefficients
  (the zero tail has no bounded differences or filtration): no transfer.
-/

namespace Fayed.R442
end Fayed.R442
