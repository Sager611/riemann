/-!
# Round 295 — node N10 (explicit prime bounds), second round: Buthe's partial-RH range is height-tied and its tail uses only Re rho < 1; the 7/8 half-plane should enlarge it about 1.8-fold

* Node: N10 explicit psi and pi - li bounds (second round; height-tied ranges).
* State of the art: |theta(x) - x| < sqrt(x) log^2 x/(8 pi) for 599 <= x <= 2.169*10^25 (Buthe, arXiv 1410.7015,
  Theorem t:partial-prime, with T = 3*10^12; Platt-Trudgian 2109.02249 Prop. 2). No 2025-2026 paper (math.NT listings
  Sep-Nov 2026 scanned) uses the 7/8 half-plane or T to extend this range; Fiori-Jaskari 2609.23222 Thm 2 (7/8) is a
  different bound (it beats the classical one from x ~ e^191 on).
* Mechanism (TeX, 1410.7015): the Logan kernel is not band-limited; zeros with |gamma| > c/eps enter only through the
  remainder (e:rem1) with x^{Re rho} <= x + 1. With c = (1/2) log x + 5 and eps = log^{3/2} x/(8 sqrt x),
  c/eps = (4 + 40/log x) sqrt(x/log x) <= 4.92 sqrt(x/log x) for x >= 10^19: that is the 4.92.
* Computation: with Re rho <= 7/8 above c/eps one may take c = (3/8) log x + 5, so T >= (3 + 40/log x) sqrt(x/log x)
  and the range grows about (4.92/3.69)^2 = 1.78-fold, to about 3.9*10^25 (estimate; the other error terms must be
  rechecked; one scout kept c and found no gain). Orca agent launched to certify it (routine, COROLLARIES).
* Adversarial: the Skewes-region recheck (Lehman sum near e^727.95) and the mixed-volume (N6) scout timed out.
* Exploration: Turan power sums (2606.24823, 2607.17505) keep classical constants, so R_0 and explicit N(sigma,T) do
  not move (N7): none. Proof mining: no application to Siegel zeros, Littlewood oscillation or Linnik's constant found
  (N10/N11): none.
-/

namespace Fayed.R295
end Fayed.R295
