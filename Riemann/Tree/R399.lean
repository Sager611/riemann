/-!
# Round 399 — node N11 (GRH families and Siegel zeros), round 8: the quoted inputs of today's two N11 discoveries (class numbers up to 10,000; numeric Chebotarev least primes) authenticated in the sources: Cureton-Thorner's Lean theorem covers every number field and every s != 1 with Re s > 31/32 (no sorry), Kadiri-Wong's Delta_L^310 is their statement and still the best fully numeric bound, Frei-Sofos-Widmer give no numeric constant, Mosunov-Jacobson's tables are unconditional to 2^40, Kimura's classno_upto is unconditional; seismic inversion, photonic crystals and polar codes give no route

* Node: N11 GRH families and Siegel zeros (eighth round; R370 before). Round type: authenticate the quoted inputs of today's
  two N11 discoveries (imqclass2: h <= 10,000; chebexp: numeric polylogarithmic Chebotarev least primes) in the sources.
* Cureton-Thorner (github.com/AxiomMath/DedekindZetaZeroFree, Zenodo 23373237): Main.lean l. 9-14, dedekindZeta_ne_zero: every
  number field K, every complex s with Re s > 31/32 and s != 1, no height restriction; no sorry or axiom: HOLDS.
* Kadiri-Wong (2105.14181) Thm 1.1: Galois L/K, L != Q, a degree-one unramified prime with Frobenius C and norm <= d_L^310; the
  introduction records Ahn-Kwon's 12577; the arXiv API shows nothing better that is fully numeric (2021-2026): HOLDS.
* Frei-Sofos-Widmer (2610.18491) Cor. 1.2 (main.tex l. 221-231): x >= c(delta) (log Delta_L)^{1/delta} (log log Delta_L)^{2/delta}
  under QRH(delta), delta in (0, 1/2]; no numeric c(delta): HOLDS.
* Mosunov-Jacobson (arXiv 1502.07953, Math. Comp. 85): unconditional tabulation (Eichler-Selberg verification removes ERH) of all
  negative fundamental discriminants |Delta| < 2^40 (TeX l. 80-88, 157-160, 879-884): HOLDS (we use it below 2.98e10).
* Kimura's classno_upto (cln-imq-2k v1.4, postproc.gp l. 25-50): prime forms up to sqrt(d/3) with qfbprimeform and qfbred, no
  GRH routine; classno_upto(28154696707, 10000) = 9819 = qfbclassno: HOLDS.
* Exploration: seismic inversion (Gel'fand-Levitan would make xi'/xi a Weyl-Titchmarsh m-function; its Herglotz property is
  Lagarias's RH equivalent Re xi'/xi > 0 on Re s > 1/2 (1999), and testing it on prolate approximants is CCM's programme, so no new
  mechanism), photonic crystals (band gaps of local self-adjoint operators, no arithmetic trace formula linking a Hill or
  Fibonacci discriminant to xi) and polar codes (Arikan's positive recursion Z^- = 2Z - Z^2, Z^+ = Z^2 discards the Moebius and
  Euler-product phases): no route.
-/

namespace Fayed.R399
end Fayed.R399
