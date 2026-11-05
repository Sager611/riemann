/-!
# Round 288 — node N2 (de Bruijn-Newman): 0.13963 still the record; its dependence on the half-plane parameter quantified

* Node: N2 de Bruijn-Newman constant.
* State of the art: Lambda <= 0.13963 (Hakimoglu-Brown, 2610.12803, from the 7/8 half-plane plus numerics) is still the best;
  2611.06245 (< 0.156821) and 2610.32694 (<= 0.1732) are weaker; no erratum found.
* Dependence on theta (scout recomputation from the paper's Theorem 1 / Remark 4: only q_0 = (2 theta - 1)^2 changes): the
  stepped barrier of its Table 2 gives Lambda <= 0.13962647 at theta = 7/8 (the stated 0.13963), 0.13958053 at
  0.874886647 (Liu 2610.26835), 0.13940839 at 0.87446205769353 (Chang-Wang), 0.13924469 at 0.874058504 (Liu Hexagon v3).
  Claiming these needs a re-run of Table 2's q_i, the gate inequalities, the polynomial-positivity certificates and the final
  numerical verification; all would be conditional on the unrefereed half-plane claims (Hexagon audit running: Orca hexaudit).
* Adversarial: the paper's density proposition constants re-checked (Z < 1.34149 < 1.345, J < 0.38734 < 0.39, S > 21/4).
* Pure de Bruijn strip bound Delta^2/2 with Delta = 2 theta - 1: 0.28125 at 7/8 (finite-height numerics are what bring it to
  0.13963). A moving-strip splice using RH up to T = 3*10^12 gives no meaningful gain (heat leakage ~10^{-3.2*10^25}).
* Exploration ideas: Hamilton-Jacobi zero dynamics (2512.17808) with certified contraction on finite zero boxes (N2);
  equilibrium potentials of polynomial iterates as obstructions to hyperbolicity (N6).
-/

namespace Fayed.R288
end Fayed.R288
