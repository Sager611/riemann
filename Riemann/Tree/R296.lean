/-!
# Round 296 — node N3 (Li coefficients), second round: lambda_n > 0 for n <= 77 T^2 from the RH height is already in Fabbian's preprint; the 7/8 version is routine

* Node: N3 Li / Keiper-Li (second round; supports the li200k Orca agent).
* State of the art: Zenodo 22839423 (G. Fabbian, 19 Sep 2026, Thm 1 and Remark B.7): RH to T = 3*10^12 gives
  lambda_n > 0 for n <= 6.93*10^26 = 77.03 T^2 (continuation 7.0858*10^26 = 78.72 T^2); Oesterle (unpublished) and
  Brown (2005) noted T^2-scale ranges. li200k's a = 1/2 bound (c* ~ 78) is therefore KNOWN; its 7/8 half-plane version
  (a = 3/8, c* ~ 105) is routine (agent told to certify it as a corollary). TREE.md N3 updated.
* Adversarial: a scout's counterexample to lambda_n >= M_n - Q_n (quartet rho = 0.01 + 3*10^12 i, n ~ 1.8*10^25:
  C = -2.0092 < -1.0630) counted the quartet's two zeros above T once; with both, M - Q = -2.1260 <= C, so the bound
  holds there (rechecked with mpmath, 60 digits).
* Computation: a smooth-density estimate giving n_max/T^2 ~ 2.6 (a = 1/2) and 3.9 (a = 3/8) contradicts Fabbian's 77.03
  and our balance check (Q_n ~ e^{c/2} (T/c) log(T/2pi)/(2pi) against M_n ~ (c T^2/2) log(c T^2) meet near c ~ 78):
  rejected. The lambda_n sanity scout timed out (zetazero too slow for 2 minutes).
* Wild: Dirichlet L-functions: Omar-Ouni-Mazhouda (1507.03431, Prop. 1) claim positivity for n <= T^2 from GRH to
  height T but dismiss the tail by conditional convergence; with Platt's GRH verification (1305.3087) only small ranges.
* Exploration: murmurations - the murmuration density is the n = 2g Frobenius-trace correction, a lower-order
  one-level-density term (2608.01337, 2504.12295): family averages only (N11/N12, weak). Homological stability - stable
  cohomology plus exponential Betti bounds give point-count errors O(q^{-1/2}) uniformly in height (2503.03861); lifting
  to twisted local systems could give uniform family moment bounds (N12/N9, idea). Lehmer/Hankel - under RH the zero
  power sums p_k = sum gamma^{-2k} are Stieltjes moments, so their Hankel determinants are positive (first four
  1.35e-2, 6.37e-10, 2.69e-23, 3.34e-43); a Hankel-positivity criterion (N3/N6), probably classical, unchecked.
-/

namespace Fayed.R296
end Fayed.R296
