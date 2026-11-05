/-!
# Round 323 — node N10 (explicit psi and pi - li bounds), round 4: Skewes and Mertens records updated; a prime in (x, x + x^{7/8}] for every x >= 1 (routine); differenced xi real-rooted for h > 0.37338 (routine, explore5)

* Node: N10 explicit psi and pi - li bounds (fourth round; R286, R295, R307 before).
* [check] Skewes: TREE's crossing 'below 1.3982*10^316' (Bays-Hudson 2000) was stale: the first x with pi(x) > li(x) is below
  e^{727.9513468} = 1.397182*10^316 unconditionally (e^{727.9513386} = 1.397171*10^316 under RH), Saouter-Trudgian-Demichel,
  Math. Comp. 2015; a crossing near e^{727.95133}. The scout's 'e^{727.9513} = 1.3971167*10^316' truncated the exponent
  (wrong in the 5th digit); corrected by the orchestrator. TREE N10 updated.
* [check] Mertens: the first counterexample to |M(x)| < sqrt(x) is below about exp(1.96*10^19) (Kim-Nguyen, arXiv 2502.21021,
  lattice reduction; known here since R118, now in TREE N10). Hurst's limsup/liminf unchanged.
* [check] pi(x) < li(x) for x <= 10^19 (Buthe) and Johnston's pi - li range to 1.101*10^26: HOLDS (nothing newer found).
* Computation (routine corollary, COROLLARIES): for every real x >= 1 there is a prime in (x, x + x^{7/8}]. Inputs: the 7/8
  half-plane, Fiori-Jaskari 2609.23222 Thm 2 (|psi(x) - x| < 0.04934 x^{3/4} log^2 x + 0.02636425 x^{7/8}), used from
  x = e^{41.56} (r = 0.99992 < 1), psi - theta < 1.4262 sqrt(x) (Rosser-Schoenfeld), the maximal prime gap 1476 below 4*10^18
  (Oliveira e Silva-Herzog-Pardi) for 4186 <= x < e^{41.56}, and a sieve to 10^7. Reviewer: correct given the inputs, routine;
  the scout's 'best explicit interval from 10^43' was wrong: Trudgian's x/(111 log^2 x) is shorter below about 10^49.
  [ShortPrimes](Riemann/Deep/ShortPrimes/REPORT.md). The scout timed out; the orchestrator finished the computation.
* Harvest: 2609.23222 states no open computational question and no short-interval result (TeX checked).
* Exploration: explore5 (Orca; Fourier quasicrystals / Lee-Yang) NO IDEA (NOGO); its routine corollary: Lagarias' differenced
  xi, A_h = (xi(s+h) + xi(s-h))/2, has only simple interlacing critical-line zeros for h > 0.37338196661320752 (was 1/2),
  since Xi(z + ih) is strict Hermite-Biehler once xi has no zeros on Re s >= 1/2 + h (reviewer: correct after writing Xi as
  its even product; routine; DH-blind). Scouts: AlphaEvolve's problem list (2511.02864): the sign-uncertainty constant
  0.3279 (Problem 11) constrains explicit-formula test functions for discriminant bounds, no mechanism toward an RH node
  (weak idea); Nevanlinna / a-point theory (2510.07710, 2008.11071): no idea (uniqueness theorems do not constrain zeros;
  DH test fails); spin glasses / log-correlated fields: scout timed out (1-RSB conjecture for log|zeta| noted).
-/

namespace Fayed.R323
end Fayed.R323
