/-!
# Round 466 — node N7 (zero-free regions), round 16, the 13 Nov listing and the tracker: QuasiRiemannTracker PR 141 (Bebel, merged 13 Nov 2026) proves Re s > 0.8733458561067185768 for every Dirichlet L-function by continuum marked anchors (admitted, standard axioms, three kernels, our static scan clean; 9.44·10⁻⁷ below PR 137); open PR 143's exact continuum contact is the degree-9 root that PR 141's rational matches to 60 digits; arXiv 2611.16665's quartic Q(i) route to Re s > 599/600 has matching inputs and its own constraint system caps it at 1 − 1/597.7; no other listing hit; Nevanlinna–Pick, Hoffman's ratio bound and Alexandrov–Fenchel transfer no tool

* Node: N7 (16th round; R463 before). Round type: listing round (arXiv Friday 13 Nov 2026, two days of submissions after the 11 Nov holiday;
  the Mac slept 22:31 UTC 12 Nov to 17:12 UTC 14 Nov) plus [check] of new QuasiRiemannTracker (QRT) PRs.
* Listing (orchestrator scan of math.NT, math.CV, math.SP /new, header 'Friday, 13 November 2026', 142 entries, max id 2611.17255, 33 keyword
  hits): the only RH-node hit is 2611.16665 (below); 2611.12438 and 2611.12467 are seen replacements; nothing beats a node record.
* [check] QRT PR 141 (Joseph Bebel, GitHub josephbebel, opened 13 Nov 04:04 UTC, merged 08:28 UTC), 'Dirichlet nonvanishing at
  0.8733458561067185768 via continuum marked anchors', submission 2026-11-13-continuum-anchor: L(s, χ) ≠ 0 for every modulus, every χ
  and Re s > θ₁₄₁ (theta_exact = 9908119238306612133369374975000622988619121777878229799765263366946606045502455625911569003364955967645745531862430068266267717
  / 11345012023616779600791582323373494582287250283647959965785871011074347975117040113347744503437234376916434871293982689424818700
  = 0.873345856106718576801184796065…), except the pole. Admission record (branch admission-records, records/141/03a1e5ad…/continuum-anchor,
  run 43382115499): status PASS, judge exit 0, target QRH.canonical_bound, allowed axioms propext, Classical.choice, Quot.sound, kernels
  Lean default, nanoda, con-ron, verified 2026-11-13T08:25:53Z; builds_on marked-anchor-induction (Sil, PR 137), arvi-repeated-marks-on-sil-
  face-20261025, arvi-selected-leaf, hailey-collet-1252270103, openai-baseline, sil-density-split-20261026. Our static scan (all 632 .lean
  files of the submission at head 03a1e5ad, 117,201 lines, comments stripped): no sorry, admit, axiom, native_decide, implemented_by,
  extern or unsafe; 1252 decide +kernel and 231 set_option heartbeat/recursion lines only. Gain over PR 137 (2183367/2500000):
  9.438932814·10⁻⁷; below the old B_2 contact by 3.611·10⁻⁵. RECORD (external, Lean-checked) → TREE N7 line; routine half-plane updates
  → COROLLARIES (Jensen d ≤ 16,143,979,969,895,705,813,461,097; Λ_χ ≤ 0.27877425654; p = 1.145022; A(0.87) ≤ 1.692204).
* [check] QRT PR 143 (Bebel, open) 'Exact continuum contact' (needs PR 142, open: integer-polynomial roots with rational isolating
  intervals as exact bounds): θ = the unique root in (13/15, 7/8) of 3888·19⁵(15t − 13)⁴(1 − t)⁵ − 7(7 − 6t)⁴(9t − 7)⁵ = 0. Orchestrator:
  one sign change on a 4000-point grid of (13/15, 7/8), root 0.873345856106718576801184796065 (mpmath, 60 digits) = θ₁₄₁ to 60 digits.
  Scout (files pre-fetched): AlgebraicContact.lean l.16–18, 48–51 and PolynomialContact.lean l.21–32 state this polynomial and interval;
  the two sides are the source ideal/gain curve and the physical slot-capacity floor (ContactCapacity.lean l.9–17); the files call it the
  closest supported contact, not a ceiling of the method; no sorry, axiom or native_decide (only #print axioms, decide +kernel): HOLDS.
* [check] arXiv 2611.16665 (Xin Hang Ji, 12 Nov 2026; section 8: idea by the author, initial proof by 'GPT6.1Sol', revised and checked by
  the author), 'Proof of the Quasi Riemann Hypothesis using Q[i]': ζ_{Q(i)}(s) = ζ(s)L(s, χ_{−4}) has no zeros in Re s > 599/600, hence ζ.
  Route: quartic Gauss sums, the angular-weight-1 quartic Kubota series (no residues), Poisson summation twice, de Faveri's zero density
  for quartic characters, a lattice average over fourth powers. (a) Input DDHL (arXiv 2306.11875v6) Prop 4.3 (TeX l.1008–1013): the ℓ ≠ 0
  completed cusp series are entire with the stated functional equation: HOLDS (scout). (b) Input de Faveri (arXiv 2610.04045v1) Corollary
  'Improved zero density estimate' (TeX l.246–262): with n = 4, d = 2 its 1 − δ_4(σ) equals the paper's E(σ) in all three ranges (1/2–3/4,
  3/4–7/8, 7/8–1; orchestrator algebra) and the height exponent 1 + 2(1 − σ)/(3 − 2σ) ≤ 2: HOLDS (de Faveri is an unrefereed v1).
  (c) [falsify] Lemma angular (l.761–806): a scout FLAW ('the integrand is not radial: P(r) and the Gauss sum are not invariant under
  r → ir') is refuted by the proof: on r = a₀t the Gauss sum is Cψ(t), the periodic factor P(a₀t)ψ(t) goes into the finite Fourier sum
  mod lcm(bλ⁷, f), and the continuous part is α(a₀)³α(t)³h(Na₀·Nt), whose zero-frequency integral vanishes: HOLDS. (d) [rate] With
  H = X^h, Y = X^y and the w-contour Re w = c, its own constraints are: top group E(1 − δ) < c/2 − 1/4; main and top errors δ ≤ (1 − c)y,
  δ ≤ h/8; dual bound δ ≤ 9h/4 − 1 − 4y; Lemma density max_σ σ + h(max(E − 1/4, 0) + Λ(σ, 1)) ≤ 1 − δ on [5/8, 1 − δ]. The paper's
  (65/144, 1/288, 13/25) satisfies all with δ = 1/600 (main term with equality, dual margin 6.9·10⁻⁵, top margin 6.6·10⁻⁵); the optimum over
  (h, y, c) is δ = 0.0016732 = 1/597.7 at h = 0.451384, c = 0.51994, y = 0.003485 (orchestrator bisection, 400,001-point σ grid): the
  quartic route is capped at Re s > 0.998327 by its own propositions (binding: density margin, the dual bound's H > X^{4/9}, the top
  group); far above 7/8 and θ₁₄₁. Not a record: TREE N7 line (second route); method cap → NOGO.
* Also checked: HAL (newest 8 Nov, all seen); Riemann.fail arena (newest record 12 Nov 11:37 UTC, Knausgård 0.673639287955…, unchanged).
* Exploration (N1): Nevanlinna–Pick/Carathéodory–Fejér: the Schur normalization F(z) = 2ξ((1 + z)/(1 − z)) of ξ saturates at F(0) = 1
  (2×2 CF determinant −|a₁|² ≈ −0.0021), so interpolation positivity carries no Weil positivity; Hoffman's ratio bound controls independent
  sets, not positive semidefiniteness (KG(4,2) bound 3 against a 2×2 form with eigenvalue −1), no canonical graph encodes Weil truncations;
  Alexandrov–Fenchel: Weil's signed form is no mixed volume (2-D test 0 ≥ 1 fails; Brunn–Minkowski was R352): no transfer.
* Scouts: 8 (listing, contact, rate, de Faveri timed out or partial; the listing, rate and de Faveri checks were redone by the orchestrator).
-/

namespace Fayed.R466
end Fayed.R466
