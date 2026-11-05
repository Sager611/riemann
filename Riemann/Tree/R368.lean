/-!
# Round 368 — node N3 (Li/Keiper), round 7: records stand; a sharper explicit zero count (2606.31441) leaves the Li ranges unchanged; 2509.06390 asks about Li coefficients averaged over characters; generalized Li criteria give no larger range

* Node: N3 Li/Keiper (seventh round; R356 before).
* [check] state of the art: Fabbian's lambda_n > 0 for n <= 6.93e26 = 77.03 T^2 (Zenodo 22839423, latest version 19 Sep 2026)
  and our 106.25 T^2 stand; the Oct-Nov listings, Zenodo and HAL give no new Li range or certified numerics beyond n <= 10^5: HOLDS.
* Computation: arXiv 2606.31441 (Counting zeros of Artin L-functions, June 2026, TeX Cor. improvedzeta) gives
  |N(T) - (T/2pi) log(T/(2pi e))| <= 0.097 log T + 4.954 (and 0.0946 log T + 5.746) for T >= 1, sharper than
  Hasanalizade-Shen-Wong; our Li certificate (Riemann/Deep/Li200k/verify.py) rerun with it on a copy gives the identical
  ranges 78.714, 105.774, 106.245 T^2 (the range is insensitive to the counting error, as in R333): HOLDS.
* Adversarial: the newest Zenodo item on Li coefficients (23370435) itself disclaims an RH proof (E_n <= B_n for every n is
  missing, p. 7): NOT TESTABLE.
* Wild: generalized Li criteria (Sekatskii 1304.7895, 1410.4384 Prop. 22, tau-Li coefficients) use the verified height only to
  bound computational error; no family member has a range beyond c T^2: NO IDEA.
* Harvest: 2509.06390 (TeX l. 191, Motivating Question) asks what can be said about Li coefficients averaged over Dirichlet
  characters mod a fixed q; for q <= 400,000 our DirLi ranges make every Re lambda_n(chi) > 0 for n <= 26.12 t0(q)^2, so the
  averages are positive there, and its lambda_2(Q(zeta_m)) > 0 for large primes m holds for every prime m <= 400,000 (routine,
  by additivity); 2508.08407 states no open question.
* Exploration: sum-product and additive energy (2403.13157 proves the converse implication, zero density implies large-value
  bounds; the scout's Balog-Szemeredi-Gowers route is not in it), sampling and frame theory (no zeta sampling or Riesz-basis
  handle found), Donaldson-Thomas wall-crossing (1408.2671 is an exposition; no zeta construction): no new route.
-/

namespace Fayed.R368
end Fayed.R368
