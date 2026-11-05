/-!
# Round 374 — node N6 (Jensen-polynomial hyperbolicity), round 8: records stand; Holland's absolute constant K (2608.08682) is not explicit although its localization constants are (candidate); Chasse's sector lemma is sharp; concentration of measure, communication complexity and Bernoulli convolutions give no route

* Node: N6 Jensen-polynomial hyperbolicity (eighth round; R364 before).
* [check] Zenodo and HAL since 20 Oct: Zenodo 23987495 (Chung, a revision claiming every J^{d,k} hyperbolic) keeps the
  false Lemma 4.0(2) (fails at (d,k) = (9,0): z_8 = 10.6458329 > y*_turn = 5.4267842), already in REFUTATIONS; nothing
  beats d <= 1.6141e25 for all n: KNOWN / HOLDS.
* Harvest: O'Sullivan 2007.13582 (xi-revised.tex l. 188): n/log^2 n >= d^{3/4}/2 implies hyperbolicity only for all
  sufficiently large d; its expansions are explicit (l. 213-219, 279-285, 636-643) but the error terms carry unspecified
  implied constants: no computable threshold (partial).
* Harvest: Holland 2608.08682 (l. 81-96): an absolute K with n^3 log^2(n+2) >= K d^5 implies J^{d,n} hyperbolic; l. 2039-2044
  says the localization constants are explicit but K is never optimized, and the proof (l. 1814-1820) leaves the thresholds
  unnamed. An explicit K plus Chasse's d <= 1.6141e25 for all n would give the first explicit hyperbolic region unbounded in d
  (candidate hollandk, N6).
* Sharpness: Chasse's sector lemma (1910.01227 Thm 3.6) is sharp in its own setting (the zero pair -5 +- i gives J^26
  hyperbolic and J^27 not); its proof gives d <= (T^2 + a^2)^2/(4 a^2 T^2), within 1 of our (T^2 - 1/4)/(4 a^2): SHARP.
* Wild: an effective GORZ saddle-point expansion with a Hermite root-separation bound and interval certification of the
  finite cases is the same route as an explicit K in Holland's wedge (candidate above).
* Exploration: concentration of measure (subgaussian tails of Dirichlet polynomials are classical moment bounds; transfer to
  zeta needs unproved independence), communication complexity (nonnegative rank of Weil Gram matrices has no link to RH),
  Bernoulli convolutions (Fourier decay of the Selberg / Bohr-Jessen law is known): no route.
-/

namespace Fayed.R374
end Fayed.R374
