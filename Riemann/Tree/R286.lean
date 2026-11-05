/-!
# Round 286 — node N10 (explicit psi, pi - li): explicit 7/8 consequence already published (Fiori-Jaskari); records re-checked

* Node: N10 explicit psi and pi - li bounds.
* State of the art: RH verified to height 3*10^12 (Platt-Trudgian, 2109.02249; no later extension found); classical explicit
  bound |psi(x) - x| < 9.22022 x (log x)^{3/2} e^{-0.8476836 sqrt(log x)} for x > 2 (Fiori-Kadiri-Swidinsky 2204.02588) and
  the matching pi - Li bound (2206.12557); first sign change of pi - li near 1.3982*10^316 (Bays-Hudson; rigorous work by
  Saouter-Demichel-Trudgian). The explicit consequence of the 7/8 half-plane is already published: Fiori-Jaskari
  (arXiv 2609.23222, Thm 2) |psi(x) - x| < 0.04934 x^{3/4} log^2 x + 0.02636425 x^{7/8} for x > e^3, which beats the
  classical bound from about x = e^{191}. TREE.md N10 updated.
* Computation: the explicit formula with 100 zeros (T = 236.52) reproduces psi(x) for x = 10^3..10^7 within errors
  -2.7..380, far below the x log^2(xT)/T truncation bound.
* Exploration: Maynard-basis prime-tuple constants H_2 <= 31,140 etc. (2610.41113), primes in [x - x^0.52, x] (2610.18923),
  Mobius cancellation in intervals H >= X^{19/35+eps} (2610.09567); idea: feed short-interval Mobius cancellation into an
  explicit Perron zero detector (N10/N11). Proposed computational record: certify Li coefficients lambda_n > 0 up to
  n = 2*10^5 (doubling Johansson's 10^5; N3) with streamed ball arithmetic: given to Orca "li200k".
-/

namespace Fayed.R286
end Fayed.R286
