/-!
# Round 452 — node N8 (critical-line proportion), round 20, mapping the Riemann.fail arena's TEN method for a ceiling: 53 accepted records since 1 Oct 2026 from κ₀ = 2 − 1/c_MT = 0.6725007 (Montgomery–Taylor) to 0.673639288, increments saturating (largest 1.49·10⁻⁴, last 3.1·10⁻⁷); TEN is a joint LP over a 16-harmonic cosine window, ten-point weights and a floor c with κ = (H_W⁻ − B)/(1 − c) on a root-covered certificate; zeta23's 2/3 is sup_{λ≤1} (2 − 1/λ − λ/3) from a rank–trace inequality; no ceiling is stated anywhere; Chirre–Gonçalves–de Laat's semidefinite pair-correlation bounds under RH (Σm ≤ 1.3208N, simple ≥ 0.6792, distinct ≥ 0.8477) lie above the arena's class, so their import into zeta23's unconditional frame is the candidate; Schäfer's Zenodo 22837912 (0.673425, 0.836713) is superseded; Kato smoothing, Stein–Tomas restriction and Bochner–Riesz transfer no tool

* Node: N8 (20th round; R451 aside updated the arena leader). Round type: map a record-producing method for its ceiling
  (lessons 172, 183: LP-duality method floors pay), before spending an Orca slot.
* [check] Arena history (GitHub commits 'record(...)'): 53 accepted critical-line records, first κ = 0.6728525 (1 Oct 2026),
  last κ = 0.673639287955 (knausgard-ten-win10-lp, 12 Nov); largest Δκ = 1.49·10⁻⁴, last 3.09·10⁻⁷: saturating, not monotone;
  README l.8: κ₀ = 2 − 1/c_MT = 0.672500703679 (the Montgomery–Taylor value); no ceiling or limit stated (rate limits only).
* [check] TEN (submission.json l.19–20 of beiyi-ten-d250-t022 and knausgard-ten-win10-lp; Solution.lean l.100–113): a joint
  LP in 16 cosine-harmonic window coefficients, ten-point weights a_ij, b_r, a floor c and auxiliaries, subject to window
  positivity with root cover, reflection, tangent/leaf and H_W-floor constraints; κ = (H_W⁻ − B)/(1 − c). No ceiling stated.
* [check] zeta23 (our port Riemann/Deep/WeilEF/Z23): the core is rank_trace_ineq (RankTrace.lean l.160–170:
  c·tr P − c²r/4 + 2c·tr Q − c²b ≤ ‖P + Q‖²_F) with a thresholded Cauchy–Schwarz (Weyl.lean l.82–94); the profile
  H(λ) = 2 − 1/λ − λ/3 has sup_{λ≤1} H = H(1) = 2/3 (Defs.lean l.57–73); certificates tighten the Block/Tail trace and
  Frobenius inputs (Assembly/Inputs.lean l.45–69); beyond 2/3 needs a new profile, which TEN's window supplies.
* [known] Chirre–Gonçalves–de Laat 1810.08843 (TeX l.255, 278, 305): under RH Σ_γ m_γ ≤ 1.3208 N, simple ≥ 0.6792, distinct
  ≥ 0.8477 (before 4/3, 1.3275, 67.27%); under GRH 1.3155, 0.6845, 0.8486; no method limit stated (only Montgomery–Taylor's
  'optimal function in this class', l.247). So the arena's κ sits between Montgomery–Taylor's 2 − 1.3275 and CGL's
  2 − 1.3208 = 0.6792: importing CGL's semidefinite class into zeta23's unconditional frame, with a dual bound for the class,
  is the next Orca candidate (record or barrier).
* [check] Zenodo: 22837912 (D. Schäfer, 18 Sep 2026): unconditional simple on-line > 0.67342504425, distinct > 0.83671252212
  (rank–trace with a pair/pressure certificate, no Lean): below 2610.08965's 0.6735382 and 0.83900, superseded; 23027645
  (G. Fabbian, 29 Sep): proportions conditional on a hypothesis (W), θ > 0.2554312, ceilings 0.8876200, 0.8362503 as θ → 1.
* Exploration: Kato smoothing (N13) needs absolutely continuous spectrum, eigenvectors are annihilated; Stein–Tomas
  restriction (N7) is weaker than Guth–Maynard and blind to Euler-product phases; Bochner–Riesz / Fefferman's ball
  multiplier (N1) gives no ζ-dependent form: no transfer.
-/

namespace Fayed.R452
end Fayed.R452
