/-!
# Round 141 — counterexample claims (Fried, Sard, Tóth–Varolin, Plotnick, bipanchromatic), Kozma–Nitzan, E₂ zeros, a Lean dBN upload

Scouts (8/8 ran):
* 2610.30063 (exposition of the proof of the Kozma–Nitzan conjecture attributed to OpenAI): candidate false step.
  The operator at l. 562, 𝒢_r F(C) = E_(off C)[E_(off C_(X_(r−1))) F(C_(x_r))], does not represent the transition of
  l. 543–550 ("off C" closes every edge at C, the transition leaves them open). Example: path x_r − u, X_(r−1) = {y}
  isolated, edge probability 1/2, F = 1_(|C| > 1), C = {x_r, u}: transition gives 1/2, the operator gives 0. The
  invariance/convergence (l. 599–604) used at l. 621–650 is then unjustified. → verify + check OpenAI's original.
* 2610.27660 (Tóth–Varolin counterexample): X_(6,7) is supported (McLean's computation); the claim for the whole
  family X_(d−1,d), odd d ≥ 3, rests on "proceeding similarly" (l. 1090–1093): an overclaim.
* 2610.27133 (bipanchromatic colourings): the 6-vertex, 3-edge counterexample and its minimality verified by brute
  force; the recolouring step of the characterization (bounds.tex l. 105–123) is unjustified.
* 2610.22281 (Fried conjecture counterexamples near hyperbolic metrics), 2610.24080 (equiregular Sard
  counterexample; determinants −8 and 32e(e−1)⁴(e+1) recomputed), 2610.21051 (Plotnick: the augmentation-49 unit in
  Z[SL₂(F₅)] verified on all 120 elements), 2610.30025 (zeros of E₂: N(Y) = 1/(2πY) + O(…); roots 0.5235217, 0.1309190
  by q-expansion): none found.
* Zenodo 10.5281/zenodo.23382724 (Lean de Bruijn–Newman bounds): conditional; Liu's zero-free region, the 7/8
  half-plane and a Polymath step are hypotheses, not proofs (4 sorry only in challenge files).
-/

namespace Fayed.R141
end Fayed.R141
