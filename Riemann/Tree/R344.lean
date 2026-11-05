/-!
# Round 344 — node N8 (critical-line proportion), round 5: records stand; Johns's Lean headline is hypothesis-free; 2610.08965 is an energy-majorant scheme (theta = 17/20) that no 2026 input lengthens; N4 record moved to the Lean-checked half-plane (10^(10^38.71))

* Node: N8 critical-line proportion (fifth round; R331 before).
* [check] Zenodo and arXiv since R331: 67.353% (2610.08965 thm:line, 1669159/2478195) stands; Johns's Zenodo 23070760 is still
  v1.0.1 (30 Sep 2026): HOLDS.
* Adversarial (faithfulness of Johns's Lean headline): Families/Headline.lean thmMain (l.149-167) takes no hypotheses; its
  statement (Families/Main.lean l.189-193) counts actual zeros, Nchi chi T = finsum over zerosI chi T of multiplicities (l.38),
  bridged to the ported zero count (Ported/Zero/Bridge.lean l.110): HOLDS (the ported definitions were not audited).
* Wild: 2610.08965 is not a mollifier method; it fixes theta = 17/20 for a Fourier-supported energy majorant (TeX l.337-347,
  442-467), and its own ceilings cap the scheme at 0.67383 (N8Opt), so the Lean-checked half-plane, the Siegel-zero exclusion,
  new exponent pairs or Kloosterman bounds cannot move it: KNOWN.
* Related records: arXiv abstracts give no lower bound for the proportion of zeros of xi' on the line (Conrey 1983 is in the
  tree only qualitatively, R72): not tracked (NOT TESTABLE this round).
* Harvest: Alpoge-Furman 2608.13637 leaves only conditional HL*(4) and higher-moment questions (TeX l.589); Cobeli-Zaharescu
  2610.23467 is in NOGO: nothing usable.
* Exploration: phase retrieval (zero flips keep |F| on the line but leave the Xi and H_t classes; no invariant), geometry of
  numbers (Siegel and Rogers averaging give identities, never pointwise positivity of Weil, Li or Nyman-Beurling forms),
  algorithmic randomness (no theorem ties the Kolmogorov complexity of mu or lambda to M(x)): no idea.
* Processed (N4): the Robin/Nicolas record restated with the Lean-checked half-plane Re s > 0.87338196661320752 instead of 7/8;
  RobinHeight check_all.py Level 2 rerun on a copy: Nicolas on [1e7, 1.31e40] PASS (372 boxes, min 4.2e-05), Robin for
  largest prime up to 1.19e39 PASS (1240 boxes, min 5.81e-06), so Robin for 5041 <= n <= 10^(10^38.7133) (was 10^(10^38.6095)).
-/

namespace Fayed.R344
end Fayed.R344
