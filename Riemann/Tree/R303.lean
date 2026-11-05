/-!
# Round 303 — node N12 (Selberg class, function fields), second round: no half-plane beyond Hecke; the GL(2) gap is a fourth-moment input; Siegel-zero exclusion makes Binyamini's Andre-Oort for Y(1)^n effective

* Node: N12 function-field and Selberg-class analogues (second round).
* State of the art: no fixed zero-free half-plane Re s > theta < 1 is known beyond finite-order Hecke L-functions
  (Liu 2610.26835; Chang-Wang 2611.00924 conditional 0.87446205769353); new: a GL(2) converse theorem over number and
  function fields (2611.02217 Thm 1.1); no Selberg-class extension.
* Wild (GL(2) transfer): Euler product, prime coefficients and the approximate functional equation hold for a newform
  L(s, f), but there is no pole at s = 1, Lambda_f(n) is not nonnegative, the degree is 2, and the amplified
  fourth-moment/product-moment estimate behind Liu's Prop. v:prop:refined-numerator (OpenAI Lemma 18.1) has no GL(2)
  analogue: that estimate (plus cubic-theta positivity) is the exact missing input.
* Computation: the Davenport-Heilbronn off-line zero 0.808517182457 + 85.699348485378 i confirmed (mpmath, two seeds,
  |f| = 3.8*10^-21): our R2 filter stands.
* Adversarial: no 2026 GRH-for-Selberg/Dirichlet/Dedekind claim beyond the known-flawed 2607.04338.
* Exploration: amoebas/Ronkin functions give nothing beyond the classical half-plane for zeta (none). Unlikely
  intersections: Binyamini's effective Andre-Oort for Y(1)^n (1809.05302 Thm 1) is ineffective only through the
  Siegel-Tatuzawa field K_*; the 2026 Siegel-zero exclusion makes it fully effective (COROLLARIES; routine, node N11).
  Discrepancy: Tao's theorem gives only unbounded partial sums of lambda (1509.05363), far from RH-scale control (N10,
  weak). Function fields: no Weil-free uniform-in-q zero-free result (none).
-/

namespace Fayed.R303
end Fayed.R303
