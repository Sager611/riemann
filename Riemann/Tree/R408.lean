/-!
# Round 408 — node N11 (GRH families), round 11: four more routine corollaries from harvest4 checked in the TeX: Kim's Goldbach numbers with missing digits, Mangerel's Littlewood-type characterisation, Bahnson et al.'s polynomial saving and Akbary-Maarefparvar's 161 imaginary quadratic fields of Polya index 2 (exact recount), all unconditional; Fuselier et al. still need GRH; collective motion, paired-comparison ranking and tidal harmonic analysis give no route

* Node: N11 GRH families and Siegel zeros (eleventh round; R407 before). Round type: check harvest4's routine corollaries
  (Riemann/Deep/Harvest4/HARVEST.md), four more of sixteen (nine done, seven left); recorded in COROLLARIES.
* Kim 2412.19975 Thm 1.2 (l. 93): the Assumption (l. 79) is a zero-free half-plane Re s > 1 - c2 for every Dirichlet L-function;
  Bebel's half-plane gives c2 = 1 - B: almost all even 2n in [X, X + H], X^{2/3+eps} << H << X^{1-eps}, with a missing digit are
  Goldbach numbers, unconditionally: HOLDS.
* Mangerel 2502.16014 Thm thm:iff(b): A-Littlewood type (l. 189-194: no zero in Re s > 1 - A log log N/log N, |Im s| <= N) holds
  once A log log N/log N <= 1 - B: the characterisation for c in (0, 1) is unconditional for large primes N: HOLDS.
* Bahnson-Daskalakis-Dohadwala-Shah 2412.07055 Rem. Polsaving (l. 1281-1291): a zeta-only strip of width delta gives the saving
  t^{-min(delta/2, 1/31)}; delta = 1/8 already gives t^{-1/31}: HOLDS.
* Akbary-Maarefparvar 2508.11125 Thm (l. 379-381): GRH only through Ihara's bound (l. 1279-1333). Index 2 means h = 2^t, so the
  class group is (Z/2)^{t-2} x Z/4 of exponent 4 and lies in EKN's exponent-4 list, complete without GRH (Exp48, largest |D| 435,435);
  exact class numbers by counting reduced forms for every d <= 435,435 (orchestrator, 2.4 s) give 161 fields, largest d = 40,755,
  equal to the paper's table, and the classical 65 with one class per genus as a sanity check: HOLDS (polya_index2.py).
* Fuselier-Iezzi-Kozek-Morrison-Namoijam 2306.03051: the scout's FLAW stands: besides Littlewood's bound (l. 808-846) GRH enters
  through Wesolowski's Thm 6.4 (l. 997-999); not routine, and Page-Wesolowski 2309.10432 is already unconditional.
* Exploration: collective motion (Cucker-Smale flocking, Vicsek, graph consensus), paired-comparison ranking (Bradley-Terry, Elo,
  Rank Centrality, HodgeRank) and tidal harmonic analysis (Doodson sums, Rayleigh criterion): dissipative or finite-graph
  Laplacian and Perron positivity, or almost-periodic least squares, with no arithmetic kernel: no route.
-/

namespace Fayed.R408
end Fayed.R408
