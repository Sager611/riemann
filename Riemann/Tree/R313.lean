/-!
# Round 313 — node N6 (Jensen polynomials), third round: all large-n thresholds are ineffective, so the explicit region stays d <= 1.6141*10^25 for all n; Dirichlet analogue blocked by Chasse's positivity hypothesis; Laguerre inequalities to height T found (2611.02944)

* Node: N6 Jensen-polynomial hyperbolicity (third round; open questions weaker than RH).
* Effectivity: GORZ (1902.07321, Theorem XiTheorem, 'sufficiently large n'), Kim-Lee (2105.05386, unknown n_1) and
  Holland (2608.08682, K only 'sufficiently large') are all ineffective, so the fully explicit region is all n with
  d <= 1.6141*10^25 (R310 update): KNOWN. GORZ's Theorem 'effective' (d <= 8, all n) is a complete claim with
  compressed numerics; no gap published.
* Dirichlet analogue (Chasse's lemma with Platt's GRH heights 1305.3087 and the Lean-checked half-plane): the degree
  bounds would be 1.99*10^15 (q = 3), 1.79*10^12 (q = 100), 112,076 (q = 400,000), but Chasse's lemma needs positive
  Taylor coefficients, which fail in general for signed chi (real primitive chi give real coefficients only): HARD.
* Laguerre: 2611.02944 proves the generalized Laguerre inequalities of every order on |x| <= T - 1/2
  (T = 3,000,175,332,800), a height-tied result (added to TREE.md N6); Jensen hyperbolicity up to degree D implies only
  the finite Newton/Turan inequalities, not all Laguerre orders.
* Wild: the forward heat flow preserves global real-rootedness, but finite-height reality of Xi does not give reality
  of H_{0.1}; HB's walls reach only |Re z| <= X, |Im z| < 1/5: no all-n range for H_t, 0 < t < Lambda.
* Exploration: Riemann's Phi is not PF_infinity in the claimed Zenodo sense (23788542 already in REFUTATIONS, R214);
  a real-stable trivariate F(y, t, X) = e^t sum gamma(n+j) y^n (Xt)^j/j! would give all J^{d,n} at once (equivalent to
  PF_infinity of gamma: known, N6); uniform saddle expansions with explicit relative errors (2207.00717 Prop.
  pr:partial) could make n_0(d) explicit (N6, idea for an Orca).
-/

namespace Fayed.R313
end Fayed.R313
