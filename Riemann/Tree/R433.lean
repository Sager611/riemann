/-!
# Round 433 — node N9 (moments, Lindelöf, subconvexity), round 12, falsification of the load-bearing inputs of our 0.568 t^{1/6} log t: Lemmas 2.3 and 2.5–2.7 of 2602.05614 hold at random and corner parameters (also at m = 1, outside the stated hypotheses of 2.5 and 2.6); 2602.05614 still v1, no erratum; Rafiq's Zenodo 22952309 confirms Hiary–Patel–Yang's interval comparisons; Zenodo 22265217's claimed RH proof breaks at eq. (10) (REFUTATIONS); filter design, geostatistics and demography transfer no tool

* Node: N9 moments, Lindelöf, subconvexity (twelfth round; R412 before). Round type: falsify the unrefereed inputs of our record.
* [falsify] 2602.05614 v1 Lemma 2.6(b) (TeX l. 355–398; any t > 0): 23,580 random and corner cases (t ∈ [1, 10^9], K ≤ 200,
  r ≤ 60, K0 = K included, Δ ∈ (0.01, 0.49)): min RHS − LHS = 0.8756 (t = 189,639.8, K = 15, r = 33, m = 2, L = 12);
  with m = 1 (outside its m ≥ 2; Lemma 2.7 allows m = 1): min 0.5821 (orchestrator rerun): HOLDS.
* [falsify] Lemma 2.7 (l. 400–419; t ≥ 5·10^6): 5,000 admissible draws to t = 10^11, half with K0 = K:
  min (RHS − LHS)/RHS = 0.579 (t = 3.13·10^10, K = 3275, r = 21, m = 1, L = 2065): HOLDS.
* [falsify] Lemma 2.5 (l. 329–346), exact rational test (t cancels, λ^{1/3} = (r+1)/r): K ≤ 160, r ≤ 40, every m < K
  incl. m = 1, K0 = K and (K+1)/2: max LHS/RHS = 0.98689 (K = 160, r = 2, m = 1, K0 = K), tending to 1 as K grows
  (sharp at K0 = K), no failure (orchestrator): HOLDS.
* [falsify] Lemma 2.3 (refined Kusmin–Landau = Hiary–Patel–Yang 2207.02366 Lemma 2.3, l. 272–285): 41,723 log,
  quadratic and near-linear phases: max |Σ|/bound = 0.96634 (near-linear, U ≈ V ≈ 0.0338, 44 terms; mpmath agrees);
  the scout timed out, the orchestrator reran its script with a 2·10^5-term cap: HOLDS.
* [check] 2602.05614: v1 only (5 Feb 2026); OpenAlex lists 3 citing works (2603.16023, 2607.10053, 2609.22624), none
  reports an error: HOLDS.
* [check] Zenodo 22952309 (Rafiq, 25 Sep 2026) re-verifies Hiary–Patel–Yang's interval comparisons (breakpoints 10^8,
  8.5·10^10, 10^12): no failure, smallest margin D2(8.5·10^10) = 0.00747, independent mpmath recomputation agrees: HOLDS.
* R432's open item (orchestrator, PDF lines and own re-derivation): Zenodo 22265217 (Aksman, claims RH) breaks at eq. (10)
  of Theorem 2: the zero mode of the completed Eisenstein series for PSL(2, Z) is ξ(2s)y^s + ξ(2s−1)y^{1−s}, not
  ξ(s)(y^s + y^{1−s}); at s = 1/2 + iγ₁ both coefficients have modulus 2.509612451·10⁻⁷ and the true zero mode at
  y = √163/2 is −4.48·10⁻⁸ ≠ 0; its Lemma 2 (C(s) ≠ 0 for all s) is false too: C vanishes at s = 1/2 + 0.451286199i,
  1/2 + 1.838024182i, … (mpmath): REFUTATIONS.
* Exploration (tool transfers, tested): Chebyshev alternation / Parks–McClellan (cos(γt) on the first three ordinates is
  not a Haar system: residuals (−1, −1, −1)), Schoenberg's conditionally negative definite variograms (CND for any real
  ordinates, no constraint; 3-point tests −0.0171, −0.0421) and Perron–Frobenius / Leslie renewal roots (controls only
  the principal root, as R432's Krein–Rutman): no transfer.
-/

namespace Fayed.R433
end Fayed.R433
