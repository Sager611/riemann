/-!
# Round 353 — node N4 (Robin / Lagarias), round 5: the 10^(10^46.79) record stands; adversarial pass on its inputs (Chirre-Helfgott 2512.15709, C_T, the k-box MAD_k) finds no flaw

* Node: N4 Robin/Lagarias (fifth round; R334 before). Record under test: Robin for 5041 <= n <= 10^(10^46.79), Nicolas for
  p_k <= 1.78e48 (robinch, Riemann/Deep/RobinCH).
* [check] state of the art: no newer verified range on arXiv (2610.39120, 2606.15096, 2511.02106 SEEN or below ours) or Zenodo
  (22014299 certifies only X <= 1000): HOLDS.
* Adversarial, hypotheses: 2512.15709v1 Corollaries cor:psi / cor:psiexpl (TeX l.389-406) assume RH for zeros with 0 < Im s <= T,
  T >= 10^7, and x > max(T, 10^9) (cor:psiexpl: every x >= 1); no upper limit on x. Our use (T = 3e12 + 1 + pi/3 below
  Platt-Trudgian's 3,000,175,332,800; 10^19 <= x <= 10^48) satisfies all of them: HOLDS.
* Adversarial, constant: TeX C_T = log^2(T/2 pi)/(2 pi) - log(T/2 pi)/(6 pi) = 113.66888055282766 at our T (independent 50-digit
  mpmath), the formula coded in check_all2.py; < 113.67: HOLDS.
* Adversarial, k-box ramp (Lemma 11_k): MAD_k is computed as an exact Fraction, then 192-bit arb (certify3.py l.108-116), so the
  k = 120 cancellation is harmless; exact values k = 2, 7, 60, 120: 1/3, 0.613873, 1.785616, 2.524186; Monte Carlo (10^6) agrees: HOLDS.
* Harvest: 2511.02106 asks whether the threshold kappa >= 2 of its effective Ramanujan theorem can be lowered to 3/2; asymptotic,
  not settled by finite verification (OPEN-HARD).
* Exploration: billiards (2507.08788's prime-labelled BKL oscillator trace is the primon gas, NOGO R7.9, R48, R330); mock modular
  forms (2506.17178 Thm 2: Hecke polynomials of j real-rooted in [0, 1728] by Rankin-Swinnerton-Dyer sign counting; the argument
  needs a remainder uniformly below the cosine main term on a compact arc, which xi lacks, so for N6 it gives only GORZ-type
  fixed-d results); large deviations (2606.27516, 2604.25579: conditional Gaussian tails only): no new route.
-/

namespace Fayed.R353
end Fayed.R353
