/-!
# Round 491 — node N5 (Nyman–Beurling / Báez-Duarte), round 17, [check] of Munawer's Zenodo records of 19 Nov 21:34 UTC: RiemannBlueprint v3.0.15's 'riemannHypothesis_iff_collapse_zero' is proved by its own hypothesis (with baezDuarteDistSq D N := D N for an arbitrary D) and WaveletRayleighRHCriterion concludes only a finite Rayleigh inequality, so the claimed machine-checked RH ⟺ Báez-Duarte equivalence is not checked (REFUTATIONS); Paper L's 'unconditional Mellin floor' is the classical d_N² ≥ (2β − 1)/|ρ|²; no RH content

* Node: N5 (17th round; R489 before). Round type: [check] of the Zenodo sweep from 2026-11-19 at 22:35 UTC (o44_zsweep.py: 10 records,
  one new since R489: Munawer's 24380524 and 24380525 at 21:34 UTC; Mohammed 24315880 updated). Texts and the Lean tree pre-extracted in
  .scratch/r491 (104 .lean files of RiemannBlueprint_Verified_Lean4_v3.0.zip). Load ≈ 3: all eight scouts answered.
* [check] RiemannBlueprint v3.0.15 (Zenodo 24380524): CMIAuditBridge.lean l.52–55, theorem riemannHypothesis_iff_collapse_zero
  (D : ℕ → ℝ) (h_bd : RiemannHypothesisStrip ↔ Tendsto (baezDuarteDistSq D) atTop (nhds 0)) : the same iff := h_bd, with
  baezDuarteDistSq D N := D N (l.47, 'abstract'): the equivalence is assumed, not proved (orchestrator read l.44–56). FLAW → REFUTATIONS.
* [check] Modules/WaveletRayleighRHCriterion.lean l.122–124: its main theorem is SCot vasyuninGCot c N ≤ tangentRayleighQuotient
  (c − cBaseAnchor) N for c ∈ A_N, a finite inequality; no RH conclusion despite the file name (orchestrator read l.118–125).
  PurePlaqueRHCriterion (l.216–239) and LogarithmicPlaqueDecayCriterion (l.65–84) are sufficiency directions conditional on
  TitchmarshMellinProperty; OpenStripRouche's theorem2_rouche_exhaustion_implies_rh_unconditional still assumes a boundary-exhaustion
  hypothesis h_exhaust (∀ R, ∃ N, c with boundary bounds; l.546–560). MatrixFactorizationAndHalfObservation l.852, l.2059: sufficiency
  only, no converse (as in R489). Static scan: no sorry or admit outside comments, no axiom declarations.
* [check] Paper L (Zenodo 24380525, v3.0.26): §7 item 2 (text l.758–759) calls (2 Re s₀ − 1)/|s₀|² ≤ αᵀM_Nα for an off-line zero s₀ an
  'unconditional Mellin right half-plane floor': this is the classical Mellin–Cauchy–Schwarz bound d_N² ≥ (2β − 1)/|ρ|² of the
  Báez-Duarte criterion (KNOWN); §7 item 4 (l.783–785) is the sufficiency direction of R489. No RH claim.
* Mohammed 24315880 (interval package for effective Selberg–Delange expansions and crossing windows): off-tree; not refereed here.
* Exploration (tool transfers, tested): Nehari's Schwarzian univalence criterion → N13 (centred at s = 2 it certifies only the zero-free
  Re s > 3/2; at s = 1/2, ξ′(1/2) = 0 breaks it), Fomin–Zelevinsky Laurent positivity → N6 (coefficient signs do not give
  real-rootedness: 1 + x + x², R456), Mittag-Leffler / Cartan B → N1 (principal parts impose no constraint on real parts): no transfer.
-/

namespace Fayed.R491
end Fayed.R491
