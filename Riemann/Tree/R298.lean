/-!
# Round 298 — node N1 (Weil positivity), second round: partial RH + the 7/8 strip + explicit S(t) bounds certify only small windows (heuristic 0.72-1.26), below our numeric 1.55

* Node: N1 Weil / Connes positivity (second round; height-tied idea).
* Idea tested: with RH up to T = 3*10^12 and |beta - 1/2| <= 3/8 above T, an off-line pair gamma0 -+ i delta contributes
  2 Re F(gamma0 - i delta)^2 = 2(C^2 - S^2) with ||S|| <= tanh(delta A) ||C|| in L^2 (F of exponential type A); if the
  zeros above T were a sampling set for the Paley-Wiener space with good frame bounds, positivity on small windows would
  follow without numerics.
* Computation (zero gaps from explicit S(t) bounds): the largest zero-free gap not excluded is h = 3.11 at t = 3*10^12
  (Trudgian 0.112 log t + 0.278 log log t + 2.510), 2.48 at 10^20, 1.86 at 10^50, 1.64 at 10^100; the newer
  Bellotti-Wong bound (2412.15470 Cor. 1.5, constant 7.20844) is worse at these heights (h = 5.11 at 3*10^12). Up to
  13 zeros may sit in an interval of length 0.1 near 3*10^12.
* Sufficient condition (frame scout): gap-only sampling bounds fail when clusters are allowed (Bessel bound unbounded);
  a toy frame model m = (1 - A g/pi)^2 against tanh(3A/8) gives A = 0.722 (g = 3.2) and 1.258 (g = 1.4): below the
  numeric record half-width 1.55, and not rigorous. The explicit negativity construction scout timed out.
* State of the art: no height-tied Weil positivity result (Suzuki 2606.09096 Thms 3-4; Connes-Consani 2006.13771:
  multiplicative support [1/2, 2]; Connes-Consani-Moscovici 2511.22755 Thm finmain: prime cutoff lambda^2 and numerical
  zeros, not positivity from zeros to height T). The general novelty scout timed out.
* Exploration: explicit finite PSD pencils for strictly hyperbolic cones (2609.24542 Thm 1) suggest finite certificates
  for single Jensen polynomials (N6/N13, weak); Connes-Kreimer renormalization: nothing acting on zeros (none);
  quantitative unique continuation applied to the heat-deformed Xi_t in (x, t) could control zero creation (2007.07034,
  2012.07169; N2/N7 idea; static Xi gives only Jensen-scale bounds).
-/

namespace Fayed.R298
end Fayed.R298
