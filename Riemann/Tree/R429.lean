/-!
# Round 429 — node N5 (Nyman–Beurling/Báez-Duarte), round 13, falsification and rate questions: no second-order term for d_N² is known (Bettin–Conrey–Farmer 1211.5191 give only C/log N under RH and a moment condition; our line had it as BDLS, corrected); 2611.03924's optimal d_N² log N (double precision, non-rigorous) falls below C = 2 + γ − log 4π from N = 640 to 5120 (0.04538, 0.04512): a negative second-order term or a numerical artefact, an open rate question; no published finite-N claim refuted; nothing new since 1 Nov; naval architecture, cartography and genomics transfer no tool

* Node: N5 Nyman–Beurling/Báez-Duarte (thirteenth round; R417 before). Round type (new methodology, AGENTS.md 11 Nov): falsify
  stated claims and ask rate questions instead of extending ranges.
* Rate: arXiv 1211.5191 is Bettin–Conrey–Farmer (not Báez-Duarte–Balazard–Landreau–Saias, whose papers are math.NT/0306251 and
  Adv. Math. 149 (2000)): d_N² ~ C/log N, C = 2 + γ − log 4π = 0.0461914179, under RH and their condition; no second-order
  term, sign or rate for the optimal d_N is stated there or found later: OPEN. TREE N5 attribution corrected.
* [check] 2611.03924, Table tab:dn: d_N(1/2)² (natural log, optimal Gram error) = 0.012740, 0.007023, 0.005283 at N = 40, 640,
  5120 (Möbius coefficients: 0.1442, 0.03659, 0.03851), so d_N² log N = 0.04700, 0.04538, 0.04512: below C from N = 640 on;
  exact Gram data to N = 1280, quadrature and Cholesky to 5120, IEEE double (non-rigorous): HOLDS as quoted. Below C at finite N
  contradicts no theorem (Burnol's bound is a liminf); under BCF the sequence must return to C, so either the second-order term
  is negative or the double-precision Cholesky of an ill-conditioned Gram matrix underestimates d_N². Candidate: certify d_N² at
  N ≤ 640 in Arb and derive the second-order term (nbrate).
* Falsify: no published or posted finite-N inequality for d_N found that these numerics could refute (2312.00211's Nyman–
  Beurling–Milman equivalences, Zenodo 24054105's reciprocal-dilate certificates: NOT TESTABLE). State of the art: nothing since
  1 Nov (arXiv, Zenodo, HAL; 2610.03568, 2606.22562, 2510.18132, hal-05852361, hal-05718412 add no bound).
* Computation: the float Gram recompute (our scout) ran out of time; the 0.0547 at N = 10 stands from R289/explore13.
* Exploration (tool transfers, tested): Kelvin's wake-angle theorem needs the gravity-wave dispersion relation (no transfer);
  Tissot's indicatrix gives K = 1 at every simple zero of ζ (checked at two zeros to 2e-20; automatic for holomorphic maps, no
  bound); Karlin–Altschul Gumbel laws need independent scores, zero spacings are correlated (no transfer).
-/

namespace Fayed.R429
end Fayed.R429
