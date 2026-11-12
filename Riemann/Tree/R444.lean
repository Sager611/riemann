/-!
# Round 444 — node N8 (critical-line zeros, large gaps), round 16, check of the published inputs of our new record Λ > 3.30 (gapceil: Bui–Milinovich's criterion with r = 2 and a degree-60 amplifier): BM's Theorem smoothed and Lemma main hold for every fixed r ∈ ℕ with the weight x^{r²−1}(x1x2x3x4)^{r−1} and any polynomial P; BBLR's twisted fourth moment needs only a_h ≪ h^ε and ϑ < 1/4; Saker's Λ ≥ 3.5555 (arXiv 1002.2695) was withdrawn by its author; no other unconditional bound above 3.2648 and no ceiling of the Wirtinger method in print; BHSJ's Λ > 2.64 is unoptimised; reaction–diffusion, generic chaining and Littlewood–Offord transfer no tool

* Node: N8 (16th round; R442 before). Round type: [check] the published inputs of a fresh record of ours (gapceil, 12 Nov 2026,
  Riemann/Deep/GapCeil: certified Λ > 3.30, h1 = 1.00000835431 ± 2.8·10⁻¹² at r = 2, ϑ = 2499/10000, v = 1.353909, degree 60;
  float ceilings of BM's criterion over all amplifiers 3.27210 (r = 1) and 3.30279 (r = 2)).
* [check] Bui–Milinovich 1410.3635v3, Thm 2.3 (smoothed, TeX l.300-313), Lemma 2.4 (main, l.370-384), definitions l.243-254: r ∈ ℕ
  arbitrary, 0 < ϑ < 1/4, P any polynomial (no P(0) = 0 condition), v bounded real; the integrand carries x^{r²−1}(x1x2x3x4)^{r−1}
  and the common factor A_{r+2}(log y)^{r²+4r}/[2(r²−1)!((r−1)!)⁴] cancels in h1: HOLDS for r = 2, degree 60.
* [check] Bettin–Bui–Li–Radziwiłł 1609.02539, Thm 1 (label mt) = BM's Thm BB (l.500-517): coefficients a_h ≪ h^ε on h ≤ T^ϑ,
  ϑ < 1/4, no smoothness condition; d_2(h)P(log(y/h)/log y) qualifies: HOLDS.
* [known] Unconditional bounds above 3.2648 by others: arXiv 1002.2695 (S. H. Saker, 2010) claims Λ ≥ 3.5555, but its v2
  (22 Jun 2010) is withdrawn by the author ('the inequality applied in this paper is not sharp'); BM's history (l.100: Selberg 1,
  Bredberg 2.76, BM 3.18) does not cite it; arXiv searches and the 22 OpenAlex citers of BM (W2609434580) give nothing above
  3.2648: NEW.
* [known] Ceiling of the Wirtinger method: none in print (BM l.98-113 conjecture Λ = λ = ∞; Hall 2.63, Bredberg 2.76); gapceil's
  float ceilings are the first: NEW.
* [check] Bui–Hall–Subirà Jorge 2511.14415, Thm 1.2: Λ > 2.64 unconditionally (amplified fourth moment plus Hall's Wirtinger
  inequality; r = 1, P = 1, ϑ = 0.1249; l.126: no attempt to optimise): HOLDS, below 3.30.
* Exploration (tool transfers, tested): Turing linear stability (needs a self-adjoint operator, RH-equivalent), Dudley's entropy
  bound (needs a centred subgaussian process; ζ on the line is deterministic), Erdős–Littlewood–Offord (no independent signs;
  anti-concentration recorded in R346): no transfer.
-/

namespace Fayed.R444
end Fayed.R444
