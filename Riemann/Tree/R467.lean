/-!
# Round 467 — node N2 (de Bruijn–Newman), round 13, the unseen Zenodo N2 records (no date filter, 80 records, 70 unseen): Koide's kernel results hold numerically — (log F)'' < −200 on t ≥ 0 for F = S'² − SS'', S(t) = Φ(√t) (our value −205.0616 at 0) and the ninth derivative Laguerre shift as the first failure (zero ±0.0248971); his Open Problem 4.14 'resolution' was claimed first by Kielhorn (Aug 2026, certificates), uncited; Koide's correlation hierarchies and Spark's dynamics paper restate known criteria and bounds; no N2 record moves; Helly, Brouwer and cluster expansions transfer no tool

* Node: N2 (13th round; R449 before). Round type: [falsify]/[check] of the unseen Zenodo records next to the node (lessons 216, 227; Zenodo
  without a date filter on "de Bruijn-Newman", "de Bruijn–Newman", "de Bruijn Newman", "Newman constant", "Newman conjecture", "Laguerre
  inequalities", "Turan inequalities", "Csordas", "Polya-Jensen", "heat flow" "Riemann": 80 records, 70 unseen, most off topic); texts
  pre-extracted (lesson 228); a cluster of new records (10–14 Nov) on Csordas' open problems for the Jacobi-theta kernel Φ of Ξ (Csordas,
  CMFT 15 (2015) 373–391).
* [falsify] A. Koide, Zenodo 24215600 (14 Nov 2026), 'Strict Log-Concavity of the Laguerre Quantity for the Riemann Xi Kernel': for
  S(t) = Φ(√t) and f = S'² − SS'', f > 0 and (log f)'' < −200 on [0, ∞). Scout (definitions from the TeX, mpmath 50/80 digits, θ series
  n ≤ 12–16): (log f)'' = −205.0616 at 0, −205.616 at 10⁻⁴, … −2962.95 at 2: HOLDS. Orchestrator (independent, Taylor coefficients of Φ at
  0, S(t) = Σ a_k t^k): f(0) = 38.0492425862565, (log f)''(0) = −205.0616020377808663 (identical).
* [falsify] Koide, Zenodo 24215798 (14 Nov 2026), 'The First Failure of the Derivative Laguerre Inequalities' (Open Problem 4.13): J_1, …, J_8
  > 0 everywhere and J_9 has exactly two real zeros ±t₉, negative between them; Zenodo 24100362 (10 Nov): J_9 < 0 for |t| ≤ 1/50. Scout
  (TeX definitions, mpmath 60/90 digits, mpmath.diff and the derivative polynomials agree): J_9(0, 0.01, 0.02) = (−7.8376, −6.4698,
  −2.6101)·10¹⁹, positive zero t₉ = 0.024897075143065087 inside the paper's bracket (and > 1/50); J_1..J_8 > 0 at t = 0, 0.05, 0.2: HOLDS.
  Siblings 24215540 (center failures of high order) and 24215545 (power-saving failure counts, via the irrationality exponent of ζ(2))
  read by abstract only.
* [check] priority for Open Problem 4.14 ((log F)'' < 0, F = (S')² − SS''): R. Kielhorn, Zenodo 21894122 (11 Aug 2026), claims it with
  exact rational interval/Bernstein certificates (his proof-status section and README); same F as Koide's f (scout: Kielhorn (2), Koide
  l.84–111); Koide cites neither Kielhorn nor C. Junquero Vila (Zenodo 22728597, 12 Sep 2026: the third level (log f₂)'' < 0, computer-
  assisted). KNOWN: Koide's qualitative statement is Kielhorn's; his uniform −200 is a strengthening. Side results on the kernel, not N2
  records.
* [check] Koide, Zenodo 24215625 (14 Nov 2026), 'Positive-Definite Correlation Hierarchies, Zero Strips, and Derivative de Bruijn–Newman
  Constants': its positive-definiteness ⟺ Laguerre–Pólya theorem is the known Laguerre/Jensen–Turán criterion (Csordas–Varga 1990); for
  'derivative constants' only 0 ≤ Λ^{(m+1)} ≤ Λ^{(m)} ≤ Λ and lim Λ^{(m)} = 0, no numeric bound: HOLDS (framework, no record).
* [check] L. Spark, Zenodo 24201014 v1 (13 Nov 2026), 'The de Bruijn–Newman Constant and the Dynamics of Riemann Zeros': restates
  Rodgers–Tao (Λ ≥ 0), Polymath (Λ ≤ 0.2), CSV and the zero ODE as theorems (TeX l.661–662); no new bound: KNOWN.
* Not read beyond titles: 20 further unseen records (Missiaen's conditional heat-flow hierarchy 24141074, Spark's Jensen trace matrices
  24222576, the 'shifted Laguerre converse' series 22242139–22763987, Turán–Jacobi recurrences 22181599/22181680, Polya–Jensen and
  Wright-circle Jensen polynomials 15283690, claimed proofs 15399563, 24048053, frameworks). Record Λ ≤ 0.107559925318 stands.
* Exploration: Helly's theorem (interval Helly number 2; ordinates give no constraint on Re ρ), Brouwer's fixed-point theorem (a
  self-map built from ξ has only boundary fixed points; existence carries no zero information), Sokal/Fernández–Procacci cluster
  expansions (need polymer weights; the test reproduces only Re s > 1): no transfer.
* Scouts: 8, all answered; the decisive value (log f)''(0) recomputed by the orchestrator.
-/

namespace Fayed.R467
end Fayed.R467
