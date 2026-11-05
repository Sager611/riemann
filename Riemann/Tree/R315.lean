/-!
# Round 315 — node N3 (Li coefficients), third round: Fabbian's 77.03 T^2 rechecked exactly; Dirichlet-L height-tied ranges estimated (not certified); normalization note on Johansson's table

* Node: N3 Li / Keiper-Li (third round).
* Adversarial: Fabbian (Zenodo 22839423) Theorem 1's number rechecked exactly (rationals and 100-digit mpmath):
  floor(3,000,175,332,800^2 (201/200)^871) = 693,318,936,098,417,738,514,687,162 = 77.0264 T^2: HOLDS.
* Computation (Dirichlet L-functions, routine but uncertified): Omar-Ouni-Mazhouda (1507.03431 Prop. 1) need 'T >= T_1
  sufficiently large' and dismiss the tail by conditional convergence; with Platt's GRH heights (1305.3087 Thm 1) a
  li200k-type argument would give roughly n ~ 6.1*10^16 (q = 3, t_0 = 3.33*10^7), 4.1*10^11 (q = 1000, t_0 = 10^5),
  4.1*10^6 (q = 400,000, t_0 = 387.5); explicit M_n - Q_n constants, T_1 and parity handling are missing (next Orca
  candidate). L(s, Delta): LMFDB lists only 20 zeros and no verified GRH height, so no range.
* Normalization: Johansson (1309.2877 Sec. 5) tabulates Keiper's normalization, lambda_100000 = 4.62580782406902...
  (= Li's lambda_n / n); the main term (1/2)(log n - 1 - log 2 pi + gamma) = 4.62613 matches to 7*10^-5. No flaw.
* Exploration: Pontryagin-space index bounds the number of off-line pairs by kappa but not conversely (1706.09103;
  N3/N1, weak); no PSLQ/experimental closed form for Li coefficients or Xi's Taylor data beyond heuristics
  (1907.01903) (none). Two scouts (centres, umbral) returned nothing.
-/

namespace Fayed.R315
end Fayed.R315
