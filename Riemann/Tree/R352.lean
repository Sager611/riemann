/-!
# Round 352 — node N10 (explicit prime-error bounds), round 6: Mertens78 power saving stands as new; explicit Liouville and Moebius-mean power savings (routine)

* Node: N10 explicit psi, pi - li and Mertens bounds (sixth round; R336 before).
* Processed this stretch: mertens78: |M(x)| <= 1.84e6 x^0.9382 for all x >= 1, the first explicit power saving for the Mertens
  function (from 7/8; Harnack and LP-certified Poisson-kernel bounds for 1/zeta; Arb checker rerun on a copy; reviewer: correct,
  significant). DISCOVERIES (N10).
* [check] novelty recheck: no explicit C x^theta (theta < 1) Mertens bound since October (2610.27591 is RH-conditional and
  non-explicit; scout timed out after the arXiv and Zenodo passes). The best previous explicit bound (2208.06141, TeX l.135-161)
  is of Korobov-Vinogradov shape, 25.85 x log x exp(-A (log x)^{3/5} (log log x)^{-1/5}) for x >= e^{e^{10.01}}; its l.168-173 say
  power-saving forms were unavailable: settled.
* Adversarial: Zenodo's 'Kernel-checked explicit Chebyshev bounds' deposit (19 Oct, 700 x^{61/64}) not found by the API; at
  x = 2^256, 10^100, 10^1000 Fiori-Jaskari's 0.04934 x^{3/4} log^2 x + 0.02636425 x^{7/8} is smaller by 10^10.4 to 10^82.5,
  so it is not a record (NOT TESTABLE as stated).
* Computation (routine corollaries, COROLLARIES): |L(x)| <= 3.28e6 x^0.9382 (L(x) = sum_{d <= sqrt x} M(x/d^2), zeta(1.8764) =
  1.7781647014) and |sum_{n<=x} mu(n)/n| <= 3.17e7 x^{-0.0618} for all x >= 1; sieve to 10^7: all three ratios peak at x = 1.
* Exploration: information geometry (2104.10548: Fisher metric of the zeta family only for real s > 1), Brunn-Minkowski /
  Prekopa-Leindler (log-concavity, NOGO R38, R56, R214), zeros of derivatives (2610.12499 Lemma 1: differentiation contracts the
  zero variance by (m-1)/m; a weak quantitative-margin idea for N6, overlapping the R34 real-rootedness toolbox): no new route.
-/

namespace Fayed.R352
end Fayed.R352
