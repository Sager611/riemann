/-!
# Round 463 — node N7 (zero-free regions), round 15, [check] of QuasiRiemannTracker PR 137 (Sil, merged 12 Nov 2026): every Dirichlet L-function, ζ included, has no zero with Re s > 2183367/2500000 = 0.8733468, the first Lean-checked bound below the contact B_2 = 0.8733819666; statement, Comparator, axioms, kernels, builds_on chain, author and a static scan check out; the gain comes from improved marked moments proved with an exact second Poisson summation and a Gauss-root phase identity, which tests our moment barrier R310 (mechanism and reach: b2below); routine half-plane updates; Lieb–Robinson bounds, smoothed analysis and Kolmogorov widths transfer no tool

* Node: N7 (15th round; R456 before). Round type: [check] of a new Lean record next to the node (lesson 136 formulation traps) and
  [falsify] of our own barrier lines against it, plus three exploration scouts. Source: github.com/GTimG/QuasiRiemannTracker PR 137,
  head 9d81a4c0b0e8338ef458f3ba5139657ae5404fde, merged 2026-11-12 20:35:22 UTC; texts in .scratch/r463 (deleted after the round).
* [check] Statement: src/Candidate/Main.lean, theorem Candidate.bound {q} [NeZero q] (χ : DirichletCharacter ℂ q) {s}
  (hs : (2183367 / 2500000 : ℝ) < s.re) (hpole : ¬ (χ = 1 ∧ s = 1)) : DirichletCharacter.LFunction χ s ≠ 0. It is the pinned
  Comparator template verifier/upstream-challenge.lean (sha256 8fb13ad977bca104…, equal to pins.json challenge_sha256; OpenAI's
  statement with 7/8 replaced by the bound and renamed QRH.canonical_bound). Mathlib's LFunction is the Hurwitz-zeta continuation and
  LFunction_modOne_eq gives riemannZeta at q = 1 (HPAudit, definitions): the bound covers ζ. HOLDS.
* [check] Admission record (branch admission-records, records/137/9d81a4c0.../43254868968-1.json): targets [QRH.canonical_bound], allowed
  axioms [Classical.choice, Quot.sound, propext], kernels [Lean default, nanoda, con-ron], judge exit code 0, verified 19:27:29 UTC; check
  runs verify (3 h 22 min) and qrh/pr-admission success. HOLDS (third-party CI; no local replay).
* [check] builds_on: exact-contact-point (PR 129), closest-fraction-contact (103), positive-reserve-selected-gap (99), selected-leaf-refined-return
  (98), arvi-selected-leaf (97), joint-detector-fine-chain (93), arvi-density (94), arvi-repeated-marks-on-sil-face (76), sil-density-split (83),
  density-split-fine-chain (86), hailey-collet (69): each has an admission record; openai-baseline is the pinned openai/math fd4aeeb2 (HPAudit A).
  Author in metadata.json: Sil (github.com/TheSil), AI-assisted. HOLDS.
* [check] Static scan (scout, comments and strings stripped): 521 Lean files, 91,255 lines; 0 hits for axiom, sorry, admit, native_decide,
  implemented_by, extern, unsafe, opaque, skipKernelTC, ofReduceBool, def riemannZeta, def ...LFunction, macro_rules, elab_rules. HOLDS.
* [check] PR conversation: no comments or reviews. Prior art (arXiv queries, lst.py, Zenodo "zero-free half-plane"): no other bound below B_2.
* [falsify our barrier] The driver (MarkedAnchor/Driver/Final.lean): 176 anchor stages -> Contract.ImprovedMoments -> ImprovedFiberMoments ->
  Return.dirichlet_of_improved. The improved inverse and ordinary moments (widths r in [1467/2000, 311/400], slot capacity (1 - w)/2 +
  552429/10^7 - (7227847/10^8) w - nu) are proved with a second Poisson summation over the Gauss rows, evaluated EXACTLY with a phase identity
  and the primitive Gauss root (Inverse/Poisson/PairPoisson, Phase, GaussRoot, Grouping: 'the exact physical-row identity'). So the first
  step below B_2 is phase-sensitive, as beyondb2 required; R310's 'stronger moments ... do not help; only cross-row cancellation' is tested
  (a scout read it as false; whether false or only too narrow is the b2below agent's task, with the mechanism's own ceiling).
* Routine updates (exact Fractions; COROLLARIES): Jensen d <= 16,143,898,339,760,043,974,392,331 (a = 933367/2500000), Dirichlet d >= 154,765
  for q <= 400,000, Lambda_chi <= 0.27877566614, p = 1.145020, h > 0.3733468, sigma* = 0.856069; certificate reruns pending.
* Exploration: Lieb-Robinson bounds (N13): a zero-count 'light cone' passes 100 zeros but ζ has no local tensor-product structure: no transfer.
  Smoothed analysis (N13): Wilkinson's root sensitivity under xi + eps g matches (7.232e-6 i) but bounds only perturbation stability: no
  transfer. Kolmogorov widths (N13): zero modes are not extremal; Chebyshev alternation was tried in R433: no transfer.
-/

namespace Fayed.R463
end Fayed.R463
