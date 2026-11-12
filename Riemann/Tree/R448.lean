/-!
# Round 448 — node N8 (small gaps between consecutive zeros), round 18, escape routes from our small-gap floor (Lamzouri's criterion cannot reach λ ≤ 0.775): the only method below it, Durkan–Pearce-Crump's two-line moments (D(λ) > 0 for λ > 1/2 under RH), uses RH essentially at Rodgers' zero moments, the contour shift through Re s > 1/2 and the boundary values of log ξ; no unconditional weighted moments of logarithmic differences exist; with Conrey's slope 8/7 (still the best) a Lamzouri-type removal would cap its length at θ < 4/7, i.e. λ > 7/11, exactly the floor 1/(1 + θ) of our Theorem A (heuristic); no unconditional result below 0.7851 in print; hypergraph containers, Bourgain–Gamburd expansion and Gaussian entire functions transfer no tool

* Node: N8 (18th round; R447 before). Round type: escape routes from a fresh barrier of ours (discoveries/small-gap-floor.md:
  Lamzouri 2611.14486 Thm Main2 cannot prove D(λ) > 0 for λ ≤ 0.775 with θ < 4/7; record λ ≥ 0.7851, ours).
* [check] Where θ < 4/7 enters Lamzouri's method (TeX l.377-389, read by the orchestrator): removing the exceptional multiset E
  needs Σ_{ρ∈E, β≥1/2} T^{θ(2β−1)} ≪ δ^{1/q}N(T), bounded by Hölder and Lemma ZeroDensity (Σ x^{β−1/2} ≪ N(T) for x ≤ T^{θ'},
  θ' < 8/7, from Conrey's N(σ, T) ≪ T^{1−(8/7−ε)(σ−1/2)} log T); the square forces 2θ < 8/7. BGST's Lemma 5 itself holds for
  supp f ⊂ [−1, 1], so the support cap is the density slope's, not the pair correlation's.
* [check] Slope: Conrey's 8/7 (Baluyot's thesis Thm 1.2.2) is still the largest unconditional slope near σ = 1/2 (2015-2026
  arXiv; 2609.39706 is degree 2): HOLDS.
* [check] Durkan–Pearce-Crump 2610.28527 Thm thm:positive-proportion (RH): RH enters at l.607 and l.629-641 (Rodgers' zero
  moments, Hadamard/Hilbert bound for the absolute moments of V_φ), l.341 (contour shift of the logarithmic difference through
  the zero-free Re s > 1/2), l.164/194 (boundary values of log ξ; Im V_φ counts zeros) and l.774-785 (integer zero counts);
  explicit formula, prime sums, Gaussian averaging, operator algebra and the Hölder counting are unconditional. Constraints:
  1/φ − 1 < θ < 1, sign of p₀(0) = sin(πφθ)sin(πφ(1+θ))/(⋯) < 0 (l.750-753), i.e. φ > 1/(1 + θ); at θ = 4/7 − 10⁻⁹,
  φ > 0.636363636769 (7/11).
* [escape] An off-line zero would carry the weight T^{θ(2β−1)} against |R_T|², as in Lamzouri's removal, so slope a needs
  2θ < a: Conrey's 8/7 would allow θ < 4/7 and λ > 7/11 = 0.6364, exactly explore34's Theorem A floor 1/(1 + θ) for the
  LP at the same support; but unconditional weighted moments of complex log-differences are unavailable (closest: Fujii's
  unweighted moments of S(t + h) − S(t), 2611.14486 Lem Fujii), so this route is open, not a theorem (candidate dpcuncond).
* [known] No 2025-2026 unconditional positive proportion of gaps below 0.7851 and no unconditional version of DPC's or
  Bondarenko–Heap's (2610.26038) method (2611.14486, 2610.22927, 2610.28527, 2610.26038, 2604.05733 checked; Zenodo none): HOLDS.
* Exploration (tool transfers, tested): Saxton–Thomason containers (no canonical forbidden-configuration hypergraph on zero
  ordinates), Balog–Szemerédi–Gowers / flattening (no group structure or Cayley operator), Offord-type hole bounds for
  Gaussian entire functions (need random coefficients; ζ's holes are deterministic): no transfer.
-/

namespace Fayed.R448
end Fayed.R448
