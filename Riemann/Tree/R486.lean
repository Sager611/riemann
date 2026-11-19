/-!
# Round 486 — node N1 (Weil positivity), round 12, [check] of new versions and notes: Mori v9 keeps both explicit λ_min bounds and makes Theorem 5 two-sided, which holds for every τ > 0.10755 (routine); Ursachi and Niiranen add nothing beyond the node line; Zeng's 11a1 zero table is wrong from the second zero (REFUTATIONS); chromatic roots, continued fractions and Ford circles transfer no tool

* Node: N1 (12th round; R468 before). Round type: [check] new versions and new Zenodo notes next to the node line (lessons 181, 201).
* [check] R. Mori v9 (Zenodo 24285436, 17 Nov 2026) against v6 (24147570, R445), pdftotext + difflib on the statements: Theorem 1's c₃
  now depends on (C₀, c_ht) (v6 l.156 → v9 l.155); Theorem 5 gains a lower bound, −log λ_min = 4πμ − τ(log μ)² − C′_τ log μ log(2 + log μ)
  ≤ …, so its limit 4π is now two-sided (v9 l.1578–1589); everything else, including both explicit bounds on our line
  (1.362·10⁹μ²³e^{−4πμ}, μ ≥ 5; 1.24·10²³μ⁵(log μ)⁵e^{−4πμ}, μ ≥ 50; v9 l.63–64), unchanged: HOLDS (scout moridiff).
* [check] Theorem 5's hypothesis: 'Fix τ > Λ'; the proof uses it only for real simple zeros of H_τ (l.1582), spacing comes from fixed-τ
  asymptotics (T_gap(τ), l.1621), and Λ ≤ 0.2 is only cited (l.1538, 2184): with Λ ≤ 0.107550869585 Theorem 5 holds for every
  τ > 0.10755 (scout morithm5; COROLLARIES, routine).
* [falsify] Mori's explicit upper bound against certified two-sided data (window [λ⁻¹, λ], μ = λ² = e^{2L}): at μ = 25, e^{3.2}, e^{3.6}
  his bounds 7.07·10⁻⁹⁶, 1.63·10⁻⁹³, 2.28·10⁻¹⁵⁵ lie above Chekmachev's ε₂₅ = 1.17·10⁻¹²³ and Zhu's U(1.6) = 8.61·10⁻¹²¹,
  U(1.8) = 7.86·10⁻¹⁸⁶; his unconditional O(log μ) is log(1.362·10⁹) + 23 log μ: HOLDS (scout morifals).
* [check] Ursachi (Zenodo 24313137, 17 Nov): Thm 1.1 proves local positivity only on [−1.1, 1.1] (below the node's 1.925); Thms 5.2/6.1
  are RH-equivalent all-scale criteria; RH stated as unproved (l.151–152): HOLDS (scout ursachi).
* [check] Niiranen (Zenodo 24259150, 15 Nov): code and data only; no note found on arXiv or Zenodo; the README's detection law
  L_det ≈ log(QT) − 0.2 (offsets δ = 0.003–0.3, T ≈ 15–101) is numerical: NOT TESTABLE (scout niiranen).
* [falsify] Zeng (Zenodo 24348441, 18 Nov), Numerical Certificate 12.6, the 14 'zeros' of L(E, s) for 11a1 below 27: wrong from the
  second entry on (PARI lfunzeros 6.3626138947, 8.6035396193, 10.0355090972, … against 6.3624972, 8.6177760, 10.1817428, …; PARI's
  Hardy function at the listed points +2.9·10⁻⁴, +0.027, −0.277, +2.43, …; orchestrator) and an independent mpmath Λ(1 + it) from the
  η-product (Λ(1) = A·L(E, 1) to 2·10⁻²⁰) is nonzero at all 14 (scout zengz): REFUTATIONS (the Li layer built on it is unsupported).
* Exploration: chromatic-root density (Thomassen, Sokal: ζ is zero-free for Re s > 1, no analogue), Stieltjes/S-fraction positivity
  (automatic for positive atomic measures; off-line ordinates break the hypothesis), Ford circles / Hurwitz approximation of zero
  ratios (no spacing information): no transfer.
* Scouts: 9, all answered.
-/

namespace Fayed.R486
end Fayed.R486
