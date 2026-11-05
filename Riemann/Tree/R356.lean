/-!
# Round 356 — node N3 (Li / Keiper), round 6: Li records stand; Platt heights behind the Dirichlet ranges rechecked; explicit zero density lifts the range only to about 114 T^2

* Node: N3 Li/Keiper (sixth round; R346 before).
* [check] state of the art: no newer certified Li range (Zenodo 23070590 computes n <= 120 only; 23370435 has no Li computation);
  Fabbian's Zenodo 22839423 unchanged (updated 19 Sep 2026): HOLDS.
* Adversarial (our Dirichlet ranges, DirLi): Platt (arXiv 1305.3087, final theorem, TeX) verified GRH for q <= 400,000 to height
  t0 = max(10^8/q, 3.75e7/q + 200) (odd q); DirLi uses t0(3) = 3.33e7, so q = 3 reaches 6.05e16 = 54.42 t0^2 and the uniform
  constant 26.12 t0(q)^2 is the block minimum: HOLDS.
* Computation (feasibility): Kadiri-Lumley-Ng's explicit zero density (arXiv 2101.12263, Thm 1 / Table 1, valid for sigma > 1/2 +
  d/log H0, so at B) at T = 3e12 bounds the zeros in [T, 2T] with beta >= 0.6, 0.75, B by 3.54e19, 2.41e14, 1.23e10; balancing
  e^{c delta} N(1/2 + delta, 2T) against (c T^2/2) log(c T^2) gives c <= 212.5, 130.6, 113.6: the half-plane end binds, so zero density
  could lift lambda_n > 0 from about 106 T^2 to at most about 114 T^2 (estimate, not certified): NO (routine, small).
* Harvest: no cheap open Li question 2023-2026; Sekatskii's generalized coefficients (1405.7354) use Re lambda_{n,a} <= 0 for a > 1/2,
  and the apparent T^2/(2a - 1) gain near a = 1/2 is a normalisation artifact.
* Wild (degree 2): no published Li range and no rigorous GRH height for L(Delta, s) or a fixed elliptic curve found (1308.6704 is a
  generic criterion): NO.
* Exploration: control theory (KYP / positive-real lemma: no arithmetic application; 2608.11520's Hankel positivity is an RH
  reformulation, not a certificate), Hopf algebras of renormalization / Rota-Baxter (QFT factorization only), Gromov-Witten / Gamma
  conjecture (zeta values in periods, no zero or positivity mechanism): no new route.
-/

namespace Fayed.R356
end Fayed.R356
