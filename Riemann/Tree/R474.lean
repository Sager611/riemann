/-!
# Round 474 — node N9 (moments, Lindelöf, subconvexity), round 15, the unseen Zenodo N9 records (no date filter, 114 records): Casañá-Giner's earlier S(T) paper is superseded and its zero-distance corollary gives lim sup |γ_k − τ_k| ≤ 0.52666 with his later constant (COROLLARIES); Santibáñez-Leal's sharp ratio holds and his proportion rests on Wang's withdrawn 2609.24167 (our E1 uses 2609.07918); Sarnowski's corridor is weaker than Ki–Kim–Lee; miruka's explicit threshold carries an unspecified O(1/log T); Nagy's Hankel positivity is automatic; N9 records stand; modular theory, property (T) and network flows transfer no tool

* Node: N9 (15th round; R462 before). Round type: [check]/[falsify] of the unseen Zenodo N9 records, no date filter (lessons 216, 227):
  Zenodo on "Lindelof hypothesis", "subconvexity", "subconvex", "fourth moment" zeta, "moments of the Riemann zeta", "exponent pair",
  "critical line" "explicit bound", "Weyl bound", "S(T)" zeta, "argument of the Riemann zeta", "mean value" "Riemann zeta", "sixth
  moment": 114 records, about 12 unseen on topic; 8 texts pre-extracted (lesson 228).
* [check] V. Casañá-Giner, Zenodo 23815260 (30 Oct 2026, a journal submission): |S(T)| ≤ 0.09783 log T + 0.7214 log log T + 10.6
  (Thm 1.3) is superseded by his 23866074 (0.08382 log T + 4.444, N9 line); Cor. 1.4(a) lim sup |τ_k − γ_k| ≤ 2πΛ(e_B) = 0.6083129
  (τ_k the Riemann–von Mangoldt points, N0(τ_k) = k − 1/2; first 100 τ_k checked, max ratio 0.263) is the S(T) bound read at the
  zeros: with 23866074's constant it gives 0.52666 (routine, COROLLARIES); Thm 1.7 (mean square of τ_k − γ_k) is Selberg's
  central limit scale: HOLDS.
* [check] F. Santibáñez-Leal, Zenodo 24227401 v0.03 (14 Nov 2026), N8: Thm 2.1 R(α, β) ≤ √2 with equality only at (0, 1), (1, 0)
  holds (2001² grid + 121 local maximisations on [0, 20]²); Thm 1.1 (0.6725008 simple critical) and Cor. 3.1 rest on Wang's arXiv
  2609.24167, withdrawn by its author on 6 Oct 2026; our Theorem E1 cites only Wang's 2609.07918 (v2, 6 Oct 2026, standing): HOLDS,
  E1 unaffected.
* [check] M. Sarnowski, Zenodo 22146266 (28 Aug 2026), Paper 134DN, N2: Thm 7.1 (high real zeros of H_t simple on the corridor
  84/29 + δ ≤ t log(x/4π) ≤ 24) is weaker than Ki–Kim–Lee 2009 (for every t > 0 all but finitely many zeros of H_t are real and
  simple); Thm 4.2 is conditional on an admissible θ (Bourgain's 13/84 supplies one): KNOWN.
* [check] 'miruka', Zenodo 17202483 (25 Sep 2025), short-window certificate (windows of height 0.25/log T, strip width 2/log T,
  T ≥ 2.4·10¹²): constants check (−ζ′(2)/ζ(2) = 0.5699609930945, C_thin ≤ 56/3), but (7.1) carries o(1) = O(1/log T) with no
  explicit constant, so the explicit threshold log T0 ≈ 6.03 (and 'explicit constants') is not proved; qualitatively fine for T
  beyond an ineffective T0: GAP (rigour), no record.
* [check] T. Nagy, Zenodo 21294742 (10 Jul 2026, draft): Hankel determinants of the zeta moments are positive because the moments
  m_2k(T) are moments of a positive measure (the pushforward of |ζ(1/2 + it)|² dt), whatever the growth (exact check for ν_k =
  2^(k²): H₀..H₇ > 0); the 'final bridge' to RH is open by the author's own account: KNOWN (trivial).
* Not read beyond titles: Chakhkiev 23284687 (critical mean values on the binary moment cone, off the node), Xu 24126280 (sixth-
  moment decoupling for Frostman measures on the parabola), Fox 20584852 (Lindelöf-type numerics at level 14), conditional
  'small solutions under Lindelöf' 10834877, frameworks. N9 records stand.
* Exploration: Tomita–Takesaki modular positivity (no map from Δ ≥ 0 to the Li coefficients; only λ₁ > 0 rechecked), Kazhdan
  property (T) via Zuk's link criterion (no canonical group action or link graph for ζ), max-flow min-cut (no canonical flow
  network for Λ): no transfer.
* Scouts: 8, all answered.
-/

namespace Fayed.R474
end Fayed.R474
