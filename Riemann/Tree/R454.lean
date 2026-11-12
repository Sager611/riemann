/-!
# Round 454 — node N8 (small gaps), round 22, falsification of dpcbeyond's Theorem U′ (D(λ) > 0 unconditionally for every λ > 7/11: Durkan–Pearce-Crump's two-line method at every mollifier length θ < κ/2 = 4/7, the fourth moment replaced by Lemma Z, a pointwise bound of the rough Möbius polynomial through the zeros, |R_T(t − iy)| ≤ C(1 + W_{−y})(1 + On)^{k−1}, plus Conrey's uniform density slope 8/7): Lemma Z holds on actual zeros (max |R_T|/(1 + On)^k = 1.47 at T = 2000, 1.53 at T = 6000), the coupled cutoff avoids the lossy |I_1|^j (767× at depth a = 20) while |I_j|/e^{θa} decays, DPC's coefficients and Newton corrections are covered, Conrey's density gives ∫W^q ≪ T without an extra log, no prior art; three referees: Lemma Z CORRECT, Lemma M/6e′ and the assembly repairable GAPs (an index |m| + 1; the range for κ > 2 and an unproved sharpness remark), all SIGNIFICANT; repairs sent to dpcbeyond, CLOSED re-reviews pending; Lieb–Thirring, the volume conjecture and Cannon's growth transfer no tool

* Node: N8 (22nd round; R453 before). Round type: falsification of a fresh theorem (lesson 198; first orchestrator's advice:
  planted configurations aimed at the new pointwise Lemma Z, three referees, then a line-by-line audit) before any DISCOVERY.
* Orchestrator reruns: n1_mechanism.py and n2_tails.py identical; n3_identity.py (zeros 408–889, t = 1000, 1003.7): the
  multivariate explicit formula behind Lemma Z holds to |LHS − RHS| ≤ 4.1·10⁻¹² (j = 1), 7.8·10⁻¹³ (j = 2); peaks ≤ 146 MB.
* [falsify] Lemma Z on actual zeros (DPC's rough coefficients, θ = 0.55, η = 0.2): max |R_T|/(1 + On)^k = 1.466 (T = 2000),
  1.530 (T = 6000), 99.9% quantiles 1.407, 1.428: HOLDS. Off-line mechanism (own code): |I_j(a)|/e^{θa} = 0.189 → 0.006 (j = 1),
  5.4·10⁻⁷ → 1.9·10⁻⁷ (j = 2) for a = 0 → 20, while the lossy |I_1|^j/e^{θa} reaches 2.14 (j = 2) and 767.6 (j = 3): HOLDS.
* [check] DPC 2610.28527 TeX l.678–683 match PROOF.md l.16; Newton corrections are p_h/h, bounded by the packet sum (l.25):
  HOLDS. Density: Lamzouri l.223–226, N(σ, T) ≪ T^{1−(8/7−ε)(σ−1/2)} log T uniformly in σ ≥ 1/2 (including σ − 1/2 ≪ 1/log T);
  Σ X^{q(β−1/2)} ≪ T log T for qθ < κ, so ∫W^q ≪ T with no extra log: HOLDS. Prior art: none (2010.10675 is unconditional
  but only near μ = 1).
* Referees (orchestrator, independent, PROOF.md readable): A (Lemma Z) CORRECT, SIGNIFICANT: coupled cutoff kept, partition
  inversion and prime-power packets correct, one weighted zero sum, uniform constants. B (Lemma M, 6e′) GAP repairable: M4's
  radius (m + 1)/L + h must read (|m| + 1)/L + h; no exponent obstruction; SIGNIFICANT. C (7′, choices, assembly) GAP
  repairable: φ > max(2/(2 + κ), 1/2) and θ < min(κ/2, 1) (or κ ≤ 2); 'Lemma Z is sharp there' unproved; 7′ valid;
  SIGNIFICANT. Repairs sent to dpcbeyond; CLOSED re-reviews of B and C, then DISCOVERY, audit and Lean agents: next orchestrator.
* Exploration: Lieb–Thirring/CLR (N13: −Ξ″/Ξ has nonintegrable multiple-zero singularities), Kashaev/Garoufalidis–Zagier
  (N13: one knot sequence, no positivity), Cannon's rational growth (N3: no finite-state encoding of ζ): no transfer.
-/

namespace Fayed.R454
end Fayed.R454
