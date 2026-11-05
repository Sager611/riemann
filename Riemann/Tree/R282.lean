/-!
# Round 282 — node N3 (Li/Keiper): the n <= 10^5 record re-checked; no extension found

* Node: N3 Li / Keiper-Li coefficients.
* State of the art: Johansson's Arb computation (arXiv 1309.2877, section on Keiper-Li coefficients) certifies lambda_n > 0
  for 1 <= n <= 100,000; no 2020-2026 extension found; RH-to-height statements do not give a larger unconditional range.
  TREE.md N3 updated ([check] cleared).
* Adversarial: Zenodo 23884485 claims only an RH-equivalence plus finite checks (no proof for all n); nothing to refute.
* Wild (heuristic): with zeros verified on the line up to T, the tail contribution to lambda_n is about
  n (log(T/2 pi) + 1)/(2 pi T), negligible even at n = T^2 for T = 3*10^12; a rigorous version needs an explicit tail bound
  (not done).
* Exploration ideas: Lorentzian-polynomial Hessians as Weil-form certificates (they only force signature (1, *), so they
  cannot give positive semidefiniteness alone; N1); Lee-Yang zero sequences of Ising/hardcore models as finite analogues for
  Jensen tests (N6/N13); tropical Newton polygons of scaled Jensen coefficients (N6/N2).
* Two computations (lambda_n for n <= 300 from Stieltjes constants; generalized Li coefficients of L(s, chi_{-3})) timed out.
-/

namespace Fayed.R282
end Fayed.R282
