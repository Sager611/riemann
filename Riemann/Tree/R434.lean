/-!
# Round 434 — node N1 (Weil/Connes positivity), round 9, falsification of new claims: Veltri's archimedean ceiling (Zenodo 24124207: the archimedean part is positive on odd test functions only below support a ≈ 0.3716) reproduced (a* = 0.37162); Mori's rigorous e^{−4πμ} upper bounds agree with every certified Chekmachev ε_c (margin ≥ 1.1·10^15, slope −12.32 vs −4π); Kastoris's prolate bound (Zenodo 22864295) has non-explicit constants and the form Cλ^A δ_4(λ) (R419's author and form corrected); no new certified positivity; quality control, epidemiology and sequence design transfer no tool

* Node: N1 Weil/Connes positivity (ninth round; R420 before). Round type: falsify new claims.
* [falsify] Veltri, Zenodo 24124207 (11 Nov 2026), §7: A(g) = (1/π)∫|ĝ|²θ′ dr + 2ĝ(i/2)ĝ(−i/2) on real odd g supported in
  [−a, a] is ≥ 0 for a ≤ 0.349 (Thm 7.2, elementary), a ≤ 0.360 (his earlier Temple certificate) and has a negative
  eigenvalue for a > 0.3715. Orchestrator (the scout timed out): odd-Legendre Galerkin (degree ≤ 5, 9, 13, 17) gives
  thresholds a* = 0.37183, 0.37166, 0.37163, 0.37162 (upper bounds, converging) and λ_min = 0.0659, 0.0334, 0.00033 at
  a = 0.349, 0.360, 0.3715 (sine basis converges more slowly to the same value): HOLDS (barrier added to the N1 line).
* [falsify] Mori, Zenodo 24116425 v4 (λ_min ≤ 1.362·10⁹ μ²³ e^{−4πμ}, μ = λ² ≥ 5) against Chekmachev, Zenodo 23848597
  (certified ε_c on L = ln c, c = 2, …, 31, 45, 47; same log window, μ = c): every c ≥ 5 is consistent, smallest ratio
  bound/ε_c = 1.124·10^15 at c = 5 (bound 8.3745·10⁻³, ε_5 = 7.4475·10⁻¹⁸); least-squares slope of log ε_c is −12.32
  (all c) and −12.36 (c ≥ 9) against −4π = −12.566: HOLDS, the certified data follow Mori's rate.
* [falsify] Kastoris, Zenodo 22864295 (single version, 21 Sep 2026; author Kastoris alone): Thm 6.2 (eqs. 87–88) gives
  |QW_λ(k_λ, k_λ)|/‖k_λ‖² ≤ Cλ^A δ_4(λ) with C, A, c₀, λ₀ existential ('A is not optimized'), δ_4(λ) = 1 − χ_4(λ) ≍
  λ⁵e^{−4πλ²}; no (log λ)³ form: NOT TESTABLE; R419's 'Connes–Kastoris, ε_λ ≤ C δ_4(λ)(log λ)³' corrected on the N13 line
  (orchestrator read the PDF lines).
* [check] state of the art since R420: new Zenodo 24120157 v5 (Wiktorowicz: A401172 predictions; its a = 2.3029 LDL
  certificate is self-declared not fully rigorous), 24126058 (near-zero spectral clusters, no positivity), 24136767
  (disclaims a proof); 2611.12152, 2611.13129, 2610.32704, 2611.02944 and Zenodo 24076153 seen before; HAL nothing newer:
  HOLDS, Chekmachev's 1.925 (2.131 even) stands.
* Veltri's Li-truncation claim ('the prime side converges to λ_n at rate X^{−1/2} for all n iff RH'): the scout timed
  out; an asymptotic equivalence that finite X can only illustrate: not run.
* Exploration (tool transfers, tested): Wald's SPRT on 'zero on the line' indicators (assumes independence; no bound),
  the epidemic final-size relation 1 − z = e^{−R₀z} (no positive renewal kernel for ζ), the Welch bound for phase vectors
  of three ordinates (tautological; max squared correlation 0.911 ≥ 1/4): no transfer.
-/

namespace Fayed.R434
end Fayed.R434
