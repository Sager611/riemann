/-!
# Round 393 — node N6 (Jensen-polynomial hyperbolicity), round 9: records stand (arXiv, Zenodo, HAL nothing new; Savva v1.1 and the RH height 3e12 unchanged); refinements of Holland's constant cannot beat the two-thirds wedge; Polya-Wiman gives no unconditional n0 for all d; real Waring rank, Hagedorn thermal duality and traffic flow give no route

* Node: N6 Jensen-polynomial hyperbolicity (ninth round; R374 before). arXiv's 6 Nov listing had not appeared by 04:22 UTC.
* [check] arXiv state of the art (axq.py: Jensen polynomials; Turan or Laguerre inequalities; Laguerre-Polya and zeta) since 20 Oct:
  nothing beyond d <= 1.6141e25 for all n or the wedge n^2 log(n+2) >= 24 d^3; Holland 2608.08682 and 2611.02944 still v1: HOLDS.
* [check] Zenodo and HAL since 28 Oct: 23863376 (higher Turan inequalities for fixed order and large n only), 23861931 (numerics
  for n <= 1000), 23849651 (Laguerre inequalities only); Savva 22293642 still v1.1 (4 Sep); Chung's newest record 23987495 (5 Nov,
  19:04) is the revision already in REFUTATIONS: HOLDS.
* [check] input height: no rigorous RH verification beyond Platt-Trudgian's T = 3,000,175,332,800 (2021) in arXiv or Zenodo
  2022-2026 (Gourdon-Demichel 2004 is not rigorous): HOLDS.
* Harvest: Holland's localization A >= 64(B + d) (2608.08682, l. 2041) could be relaxed, but no constant K helps: Holland's wedge
  has exponent 5/3 against 3/2 for the two-thirds wedge, which at d = 1.6141e25 starts at n = 3.4e37; Holland's would beat it
  only with K < 2.7e-10 (1.4e-17 at d = 1e40): DEAD (superseded by savvak).
* Wild (orchestrator's check; the scout timed out): the Polya-Wiman theorem (Craven-Csordas-Smith 1987, Kim 1990) gives F^(n) in
  the Laguerre-Polya class for large n, i.e. J^{d,n} hyperbolic for every d once n >= n0, only for functions with finitely many
  nonreal zeros; for F(x) = sum gamma(k) x^k / k! that means all but finitely many zeta zeros on the line (unknown). The
  hypothesis cannot simply be dropped: sum x^{2k} / (k!)^4 has genus 0 and every derivative violates a Turan inequality. Ki-Kim's
  Fourier-Polya theorem (Duke 2000) counts nonreal zeros of genus-0 functions through real-axis derivative data, which for F is
  an RH restatement; R072's large-k result is for t-derivatives of Xi, not derivatives in x = z^2: NOT TESTABLE (no unconditional
  n0, consistent with the d-dependent wedges of GORZ, O'Sullivan, Kim-Lee, Holland and Savva).
* Exploration: real Waring rank (Comon-Ottaviani conjecture, proved by Causa-Re and Blekherman: a real binary form of degree d has
  real rank d iff all roots are real) restates hyperbolicity of each J^{d,n}; a certificate would be the nonexistence of a shorter
  real decomposition, and a degree-d apolar form with d real roots only shows rank <= d, true for every real binary form; no
  uniformity in d. Hagedorn growth and thermal duality (Atick-Witten) give radius and symmetry only; total positivity of string
  state counts would concern another sequence. Traffic flow (TASEP and ASEP stationary weights, Lee-Yang zeros of their partition
  functions; closest 2609.23942, a fermion Lee-Yang theorem): no Weil, Robin or Jensen input. No route.
-/

namespace Fayed.R393
end Fayed.R393
