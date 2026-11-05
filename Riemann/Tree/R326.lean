/-!
# Round 326 — node N5 (Nyman-Beurling / Baez-Duarte): 2611.03924's optimal d_N reaches N = 5120, not 10^6 ([check] corrected); no newer criterion; weighted variants give no new unconditional statement

* Node: N5 Nyman-Beurling / Baez-Duarte (fifth round; R312 before).
* [check] 2611.03924 (Table tab:dn, sec. 4): the optimal distance d_N is computed by Gram-matrix quadrature only up to N = 5120,
  with d_N^2 log N = 0.0470, 0.0475, 0.0454, 0.0455, 0.0466, 0.0451 for N = 40 ... 5120 (ours: 0.0525, 0.0456, 0.0467 at
  N = 10, 50, 200; consistent); their N = 10^6 run is the Moebius-coefficient (Baez-Duarte natural) error, not d_N. TREE N5
  said 'to 10^6': corrected.
* State of the art: no 2026 paper improves a d_N bound or the criterion (the only NEW hit, 2610.32672, is about prime pairs).
* Harvest: 1211.5191 (TeX l.61) leaves open the RH-conditional convergence rate of d_N; its asymptotic needs RH plus a
  zero-sum condition (l.90-99); Gonek's stronger conjectural rate (l.100-103).
* Wild: weighted Bergman / weighted Hilbert versions (1911.04029, 2310.03972) are RH-equivalent reformulations, not
  half-plane criteria; no weaker-norm statement turns the 7/8 half-plane into a new density result (no idea).
* Adversarial scout (Burnol tail check) timed out; the recorded bound C - 6e-12 (R301) stands.
* Exploration: prime-state entanglement (Latorre-Sierra 1302.6245) encodes pi(x) by QFT estimates, nothing on pi - li or RH
  (no idea); persistent homology of primes (2512.22079, 2605.06352): no zero statistic (no idea); machine learning of mu
  (2604.02480 is about network expressivity): no idea.
-/

namespace Fayed.R326
end Fayed.R326
