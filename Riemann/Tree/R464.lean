/-!
# Round 464 — node N8 (critical-line proportions and simple zeros), round 27, falsification and novelty of Theorem E1 (santiaudit): a positive proportion of simple critical zeros in every window (T, T + T^θ], θ > 1/2, from an absolute Bessel/Hankel bound for Santibáñez-Leal's dual sum (mollifier length ν < (3θ − 1)/4); complex shifts, the Hankel identity and three block referees hold; no prior result below θ = 0.527; DISCOVERY; Aubry–Mather theory, submodularity and differential privacy transfer no tool

* Node: N8 (27th round; R459 before). Round type: the full gate of lesson 198 for a fresh theorem of ours (Riemann/Deep/SantiAudit/E1.md,
  12 Nov 2026): orchestrator rerun, falsification scouts aimed at the new estimate, three block referees, adversarial novelty.
* Theorem E1: for every fixed θ in (1/2, 1), liminf S(T,T^θ)/N(T,T^θ) >= h(θ; κ1) >= h(1/2; κ1) = 0.0055353042386, κ1 = 0.0844779079697
  (detector P, Q of Santibáñez-Leal's EXP-010 set at ν = 3/25, R = 177/20); E1-M: his short-window moment for every ν < (3θ − 1)/4.
* Rerun (orchestrator): e1_constants.py on a copy in /tmp (Arb, 4000 bits): output identical.
* [falsify] complex shifts α = 0.01 + 0.1i, β = 0.02 − 0.1i in the per-term size test (T up to 25600, pq up to 256): S/pred 0.073–0.665,
  r 0.923–8.362, decay at 4N0 <= 1.58e-3, the same scaling as with real shifts: HOLDS.
* [falsify] Step 2 (Hankel form of B = sin(πd/2)J_d + cos(πd/2)Y_d, d in {0, ±0.2, 0.1+0.15i, 0.05−0.2i}, X in {1, ..., 100}): identity error
  <= 3.0e-8, |a_±| <= 1, |X a_±'| <= 0.0977, |X² a_±''| <= 0.159: HOLDS.
* Referees (exact blocks pasted): architecture CORRECT, ROUTINE (η with θ − η > 1/2 and ν < (3(θ − η) − 1)/4; fixed J; residues need only
  ν < 1/2; Wang with λ ↑ θ); Steps 1–5 CORRECT, SIGNIFICANT (no poles crossed; uniform symbol bounds; each integration by parts costs
  T/(H X0); K > 3/2 + 2ε); constants and conclusion CORRECT, ROUTINE (frozen (P, Q) admissible at any ν; h > 0 iff κ > k*(θ); uniform lower
  bound, liminf for fixed θ).
* [novelty] GXL evidence API (one claim), arXiv abstracts ('short intervals' + 'simple zeros' / Levinson / 'critical line' + proportion),
  Zenodo without a date filter ('"simple zeros" "short intervals"', '"Levinson" "short intervals"'), HAL: the closest results are Wang
  2609.07918 Thm 1.1 (θ > 0.5501940, TeX l. 115–145), Steuding 2002 (θ >= 0.552) and Santibáñez-Leal 23134787 (θ >= 0.527); none reaches
  θ <= 0.527 (the scouts' two 'KNOWN' verdicts name these earlier results, not E1). DISCOVERY (N8): discoveries/simple-zeros-short-windows.md.
* Exploration: Aubry–Mather theory (N13): no twist map or action functional on zeros: no transfer. Submodularity (N3): Li-zero contributions
  are modular; greedy adds nothing: no transfer. Differential privacy (N1): Laplace sensitivity of single Weil terms gives only noise
  calibration: no transfer.
-/

namespace Fayed.R464
end Fayed.R464
