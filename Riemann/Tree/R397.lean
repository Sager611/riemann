/-!
# Round 397 — node N13 (Hilbert-Polya and spectral approaches), round 10: records stand; CCM's first missing step authenticated in the TeX (smallest eigenvalue of QW_lambda simple with an inversion-invariant eigenvector, fixed lambda with N to infinity), so our 'every lambda in [sqrt 2, sqrt 20]' is exactly their statement on that range; Zhu 2608.24827 v2 covers lambda = e^0.8 only; OEIS A401172 still lacks our 8 certified terms; no new operator claim; Suzuki's Fredholm equation is a numerical harvest candidate; mathematical finance, reaction networks and crystallography give no route

* Node: N13 Hilbert-Polya and spectral approaches (tenth round; R388 before), used to authenticate the inputs of our N13 discovery.
* CCM 2511.22755: the first missing step (TeX l. 1389-1391) asks that the smallest eigenvalue of QW_lambda on L^2([1/lambda,
  lambda], d*u) be simple with an eigenvector invariant under u -> 1/u; Theorem 5.10 (l. 1085-1088) then gives self-adjoint
  D_{lambda,N} for every large N; lambda is fixed while N -> infinity, then lambda -> infinity (l. 1219-1224). Our certified 'every
  lambda in [sqrt 2, sqrt 20]' (WeilSimple REPORT l. 8-24) is exactly that statement on the range: HOLDS.
* Zhu 2608.24827 (v2, 2 Sep 2026), Theorem thm:parity: lambda = e^0.8 only (even and odd gaps 2.0857e-12 and 8.2065e-15), no
  all-lambda claim: HOLDS (our 'before').
* OEIS A401172: its b-file terms n = 2..9 (2, 7, 12, 17, 22, 27, 32, 37) agree with ours; our certified n = 10, 11, 12, 13, 15, 16,
  18, 19 (42, 48, 53, 58, 69, 74, 85, 90) are still absent from OEIS: HOLDS (new).
* [check] Zenodo, HAL and arXiv since 4 Nov: Zenodo 23441946 and 23442026 are conditional symmetry or determinant criteria; 23555113
  and 23788542 are already refuted or known; no new certified CCM window and no unconditional operator: HOLDS.
* Harvest (candidate, low): solve Suzuki's Fredholm equation for v_+-(a, x) numerically and test his EQ_112 at theta = pi
  (2606.09096, TeX l. 3693-3699), and probe whether J_a(z, w) develops off-diagonal singularities as a grows (l. 3372-3378).
* Exploration: mathematical finance (Breeden-Litzenberger and Dupire positivity are one-dimensional convexity statements, no
  backward-heat or Weil-kernel transfer), chemical reaction networks (already R370: positive steady states, not zero dynamics) and
  crystallography (Patterson and Karle-Hauptman positivity is phase-blind autocorrelation, no arithmetic constraint on Weil's signed
  prime kernel): no route.
-/

namespace Fayed.R397
end Fayed.R397
