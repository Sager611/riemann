/-!
# Round 378 — node N13 (Hilbert-Polya and spectral approaches), round 8: no new operator claim; CCM's second missing step is a proximity of k_lambda to xi_lambda with no stated tolerance, and by Hurwitz plus the Lean-checked half-plane a rate O(lambda^-0.3734) suffices (candidate ccmstep2, low); the ground state stays simple and even numerically at c = 25, 30, 40; crystal bases, renewal equations and point vortices give no route

* Node: N13 Hilbert-Polya and spectral approaches (eighth round; R359 before).
* [check] state of the art: the Oct-Nov listings (math.NT, math.SP, math-ph, math.FA), Zenodo ('"Hilbert-Polya"', '"spectral realization"')
  and HAL give no new operator claim; the scout's 'new' item, Suzuki 2610.32704 (assuming RH, the positive eigenvalues of a finite-interval
  self-adjoint extension converge to the zeta ordinates), is still v1 and was verified in R260: HOLDS.
* Harvest: 2511.22755 sec. 8 (TeX l. 1388-1392) names exactly two missing steps. The second: k_lambda = E(h_lambda) (h_lambda the combination
  of the prolate functions h_{0,lambda}, h_{4,lambda} with vanishing integral; its Fourier transform tends to Xi uniformly on closed substrips
  of |Im z| < 1/2, Lemma hermfact1, l. 1353) must approximate a multiple of xi_lambda well enough that the zeros of xi_lambda-hat tend to the
  zeta zeros; no norm or tolerance is stated. Our remark (routine): |f-hat(z)| <= ||f||_{L^2(d*u)} (2 ln lambda)^{1/2} lambda^{|Im z|} on
  [1/lambda, lambda]; all zeros of Xi lie in |Im z| < 0.37338196661320752 (Lean-checked half-plane), so by Hurwitz, step 1 along a sequence
  lambda_n -> infinity plus d(lambda_n) = inf_c ||k_lambda - c xi_lambda|| = O(lambda_n^{-0.3734}) implies RH (without the half-plane one
  needs o(lambda^{-y}) for every y < 1/2). A certified d(lambda) at our 52 weilsimple windows would test the rate numerically only
  (candidate ccmstep2, low; CCM's Fig. fpro1 already shows the trend).
* Harvest (nearby): Zhu 2608.24827 open problems (i) (TeX l. 2416-2419, a conditional proof of its law via multi-interval Landau-Widom
  asymptotics) and (ii) (l. 2420-2422, the O(1/ln N) correction, an old candidate); 2601.12133 Conj. 4.1 (lim E(kappa) ln kappa exists):
  asymptotic, not checkable at fixed lambda.
* Computation (orientation, not rigorous; Riemann/Deep/WeilSimple/eigs_quick.py): lowest even eigenvalue e0 = 1.88e-98 (c = 25, N = 70),
  1.2467e-114 (c = 30, N = 80, rerun here), 1.42e-128 (c = 40, N = 80, not converged in N); odd0/e0 = 8131, 10981, 10300 and
  even1/e0 = 3.7e7, 6.4e7, 5.5e7: the ground state stays simple and even beyond the certified c <= 20: HOLDS.
* Wild: a de Branges / Hermite-Biehler bridge E_lambda with A_lambda = xi_lambda-hat closes under locally uniform limits exactly as
  Hurwitz does, so it needs the same convergence as CCM's second step: no new route.
* Exploration: crystal and canonical bases (Kazhdan-Lusztig and Soergel positivity are finite coefficient positivity, no bridge to Weil,
  Li or Jensen positivity), renewal and delay equations (Euler-Lotka, Pontryagin, Krein-Rutman give one real dominant root, not all
  roots on a line; 2107.11348 concerns rightmost roots only), point-vortex equilibria and Adler-Moser / Calogero-Moser zero dynamics
  (1709.05632: finite real-rootedness, no bridge to the heat flow or Jensen polynomials): no route.
-/

namespace Fayed.R378
end Fayed.R378
