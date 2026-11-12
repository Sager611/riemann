/-!
# Round 461 — node N5 (Nyman–Beurling / Báez-Duarte), round 13, the unseen Zenodo N5 records: none breaks a node record or our certified d_N²; a digamma closed form of the Báez-Duarte Gram (Zenodo 22831678) matches our Vasyunin entries to 10⁻¹⁶; a sum-zero 'scaling transition' note (22980911) works in L²(0,1) and is consistent with our certified L²(0,∞) values; the rest are conditional frameworks, off-node or cranks; Kesten's theorem, self-concordant barriers and percolation transfer no tool

* Node: N5 (13th round; R429 before). Round type: [check]/[falsify] of the unseen Zenodo records next to the node (Zenodo without a date
  filter on "Nyman-Beurling", "Baez-Duarte", "Beurling criterion": 31 records, 15 unseen), plus three exploration scouts. All five node scouts
  timed out on PDF work; the decisive checks were rerun by the orchestrator.
* [falsify] Zenodo 22980911 (26 Sep 2026), 'Scaling Law Transition in the Beurling-Nyman-Baez-Duarte Criterion under Sum-Zero Constraints':
  its Table (d^2 unconstrained / constrained with sum c_n = 0) works in L^2(0,1) with rho_n(x) = {1/(nx)} (TeX l. 40-45); the L^2(0,1)
  distance is at most Baez-Duarte's L^2(0,infinity) distance, and its values 1.006e-2 (N = 100), 7.442e-3 (500), 6.662e-3 (1000), 6.285e-3 (2000)
  respect that bound against our certified d_N^2 (Riemann/Deep/NBRate/results.txt: 0.0101887625 at N = 100; 0.0079537 at 320, 0.0070171 at 640,
  0.0063611 at 1280; d_N^2 is nonincreasing): HOLDS (consistent; the constrained ratios are numerical).
* [check] Zenodo 22831678 (16 Sep 2026), 'Proof of the closed-form Nyman-Beurling Gram formula': G_jk = (1/L) sum_{i<L} c_i [psi((i+2)/L) -
  psi((i+1)/L)], L = lcm(j,k), c_i = r_j(i) r_k(i)/(jk), r_n(i) = (i+1) mod n, the Gram of rho(1/(nx)) - rho(1/x)/n (G_11 = 0): equals our
  Vasyunin-based entries A[2,2], A[2,3], A[3,5], A[4,6], A[5,7], A[6,10] to within 5e-16 (mpmath, 30 digits; /tmp check): HOLDS (an equivalent
  closed form; its novelty against Vasyunin 1995 not judged).
* [check] Kline, Zenodo 23462535 (19 Oct 2026), 'The secondary eigenvalues of a totient Gram matrix and their limiting operator': Q_N = R_N^T R_N,
  R_N(j,k) = #{1 <= a <= j/k : gcd(a,j) = 1}; lambda_{r+1}(Q_N)/N^2 tends to the r-th eigenvalue of an explicit integral operator; a
  reproducible release, not tied to an RH criterion: off the node.
* [check] Zenodo 22680292 (9 Sep 2026), 'A Unified Framework for the Computational Verification of Equivalent Formulations of the RH':
  empirical scaling indices of M(x), L(x), theta(x) to x = 10^8 only: below every node record.
* Zenodo 21337373 (a 'partial obstruction map' for a variance route; its Theorem 4.2b 'does not refute H3') and 20924818 (a residual-exclusion
  framework 'conditional on A-E'): conditional frameworks, no checkable bound; 20405591 (a Q_N 'route' operator), 21693904 ('saturation law'),
  23459912 ('six routes'), 19449236, 20367106, 23696594, 23701874: not read beyond the description (cranks or expository).
* Exploration: Kesten's amenability criterion (N13): needs a group or Markov operator; a diagonal operator built from ordinates has spectral
  radius 1 from unbounded ordinates alone, blind to real parts: no transfer. Self-concordant barriers (N13): -log|xi(1/2+it)| is not convex
  between zeros, so no Dikin-ellipsoid certificate: no transfer. Percolation (Harris-FKG, N7): needs a product measure; the zero set is
  deterministic: no transfer.
-/

namespace Fayed.R461
end Fayed.R461
