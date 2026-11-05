/-!
# Round 316 — node N4 (Robin/Lagarias), third round: Morrill-Platt's range stands; Nicolas' height-tied range needs a tail bound; an explicit Jensen threshold traced to GORTTW Thm 1

* Node: N4 Robin / Lagarias (third round).
* State of the art: nothing beyond Morrill-Platt 1809.10813 (5041 <= n <= 10^{10^{13.11485}}). 2510.23889 (least colossally
  abundant exception) states 5040 < n <= 10^{10^{13.099}} (Thm 7.1), a smaller range; 2602.15905 (MacArevey) proves that a
  least counterexample to Lagarias' inequality must be superabundant (B_n increasing for n >= 55), no new range. Zenodo 2026
  Robin deposits 23765626, 22014299, 21725356, 20931450, 20733714 seen by title only.
* Adversarial: 2610.39120 Q1 (CA2 doubling) is already recorded (R206): KNOWN; it holds for the tested steps (log N <= 6*10^8).
* Computation (height-tied Nicolas): no published finite verified range for Nicolas' primorial criterion was found
  (abstract search; 2509.11182 only discusses the sign of the Mertens remainder). Buthe's |theta(x) - x| < sqrt(x) log^2 x/(8 pi)
  up to B = 2.169*10^25 does not by itself give Nicolas at p_k <= B, because the Mertens-product remainder at x depends on
  theta(t) - t for all t > x. With RH to T and Re rho <= 0.87338 above T a smoothed explicit formula might reach x ~ T^{8/3}
  (heuristic, R294): Orca candidate (enlarge Buthe's B with the half-plane; feeds Nicolas, Axler's t-free family and N10).
* Open questions: 2510.23889 and 2602.15905 contain only RH-equivalent ones (a least colossally abundant counterexample).
* Wild: the sigma_{-a} analogue of Robin's criterion (zeros with Re s > a) is still unchecked (two scouts timed out).
* Exploration: (N1, computation) off-diagonal blocks of the smooth archimedean Weil kernel have numerical rank 4-14 at
  tolerance 1e-12 (n = 300, L = 1.6), so hierarchical-matrix compression could cut the memory of Weil certificates (idea; prime
  terms and rigorous multiprecision errors not addressed). (N1) KYP/IQC certificates need a rational symbol; Weil's symbol
  (digamma term, infinitely many delays) is not rational: no idea. (N6) GORTTW 1910.01227 Thm 1: J^{d,n} is hyperbolic for
  n >= c e^{d/2} with an inexplicit c (checked in the TeX); tracking its constants would give the first explicit large-n
  threshold (Orca candidate). Six of twelve scouts timed out.
-/

namespace Fayed.R316
end Fayed.R316
