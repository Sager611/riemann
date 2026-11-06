/-!
# Round 405 — node N10 (explicit psi and pi - li bounds), round: the inputs of today's prime-between-powers discovery authenticated in the TeX (Kadiri-Lumley-Ng Lemma 4.14, Theorem 1.1 and Table 1; Hasanalizade-Shen-Wong; Platt-Trudgian; Visser states k = 8 only, a wording fix); an elementary HAL proof of primes between all consecutive cubes refuted (circular); a newer explicit Ingham-type density (2507.15184) is a candidate refinement; cellular automata, robot motion planning and nuclear structure give no route

* Node: N10 explicit psi and pi - li bounds. Round type: authenticate the inputs of a fresh discovery
  (discoveries/primes-between-powers.md, ingham3).
* Input authenticated: Kadiri-Lumley-Ng 2101.12263: Lemma 4.14 holds for delta > 0 and every sigma > 1/2 + d/log H_0,
  H_0 = 3.0610046e10; Theorem 1.1 is stated with delta >= 1; Table 1 row 0.60: k = 0.5, mu = 1.251, alpha = 0.288, delta = 0.3140,
  d = 0.341, H = H_0 - 1, as quoted: HOLDS.
* Input authenticated: Hasanalizade-Shen-Wong 2107.06506 Cor. 1.2 (T >= e; 0.1038 log T + 0.2573 log log T + 9.3675) implies
  N(t) <= t log t for all t >= 14 (minimum margin 27.5 at t = 14; N(t) = 0 below 14.1347); Platt-Trudgian 2004.09765: RH to 3e12: HOLDS.
* Before authenticated: 2601.15564 (l. 122-156): cubes computed to n^3 <= 1.649e40, the all-n problem open, the explicit route
  'futile', Cully-Hugill's thesis bound exp(exp(32.537)); Visser 2610.20058's theorem (l. 221-276) states k = 8 only, from the
  7/8 half-plane; all k >= 8 is our ShortPrimes: wording fixed on the page (FLAW in our quote, not in the result).
* Adversarial: Diouf, HAL hal-02680213v1 (2020), 'Primes between two consecutive cubes for all positive integers': Corollary
  2.2.2 (p. 117) replaces c_{n+1} by c_n without proof (Lemma 2.1 gives only c_{n+1} <= c_n); with c_n = floor(g_n^2/(8 p_n)) both
  vanish whenever g_n^2 < 8 p_n, and the key bound c_n < p_n/g_n is essentially g_n^3 < 8 p_n^2, the claim itself: circular
  (REFUTATIONS). Our check: all 17,982 consecutive prime pairs below 2e5 have c_{n+1} = c_n, as expected: FLAW.
* Harvest: explicit zero-density papers after KLN: 2507.15184 (an explicit Ingham-type bound B_1 T^{3(1 - sigma)/(2 - sigma)}
  (log T)^{(7 - 5 sigma)/(2 - sigma)} for sigma in [0.5, 0.875], T >= 3e12) gives about 0.77 times KLN's count at sigma = 0.85,
  T = 1e13 (not enough for k = 6 for every n, which needs about 1/6): CANDIDATE (refined thresholds, after auditing its
  constants); 2311.05136 (sigma >= 0.98) and 2405.12545 (constants unsupported, REFUTATIONS) do not apply.
* Exploration: cellular automata (Moore-Myhill's Garden of Eden theorem concerns amenable groups; de Bruijn graph zeta functions
  are rational with finite-matrix spectra), robot motion planning (Farber's topological complexity is cup-length data) and
  nuclear structure (shell-model spectra are finite self-adjoint matrices with random-matrix statistics, no Euler product):
  no route.
-/

namespace Fayed.R405
end Fayed.R405
