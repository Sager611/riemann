RECORD: for every 6 < A < 12, M_A(T) = ∫_0^T |ζ(1/2+it)|^A dt ≪ T^{μ*(A)+ε} with μ*(A) below the best published bound (e.g. M_8 ≪ T^{988/671+ε}, saving 37/1342 ≈ 0.0276 against Heath-Brown's 3/2 instead of 0.0115); the sixth (37/30) and twelfth (2 − 520/445121) exponents of de Faveri–Pandey are already optimal for their method.

## Claim
With de Faveri–Pandey's large-value theorems (2610.01035 v1, Thm 8.1 for all k ≥ 2) combined optimally with the 4th and
12th moments, Heath-Brown's auxiliary bound, Ivić's exponent-pair bound over all known exponent pairs and Bourgain's 13/84:
- μ*(A) = 7A/60 + 8/15 on [4, 128/17]; 173A/1342 + 296/671 on [128/17, 1340/157]; 235A/1793 + 758/1793 on
  [1340/157, 2564/279]; 443A/3361 + 1398/3361 on [2564/279, 4172/433]; 715A/5409 + 2230/5409 on [4172/433, 6164/619];
  (58907A + 182838)/445121 on [6164/619, 3587/288].
- Table (δ = 1 + (A−4)/8 − μ): A = 7: 27/20, δ = 1/40; 8: 988/671, 37/1342; 9: 2873/1793, 0.02266; 10: 771908/445121,
  0.01585; 11: 830815/445121, 0.00851. A = 5, 6, 12 unchanged: 67/60, 37/30, 2 − 520/445121.
- Optimal: no other choice within their chain (W per N, N, k, G, the Υ-capped E_2 sup, the σ-ranges) and no other
  published input gives less (exact vertex certificate over 13823 bounds plus tails).
- Range above 12: the improvement over classical bounds holds exactly for 12 ≤ A < 3587/288 (ν = 131/288, paper: "a small
  ν"); there it equals Hölder between their M_12 and the classical M_{3587/288} ≤ T^{593/288}, so it is explicit, not new.
- Blocking inputs. Sixth: Iwaniec's short-interval fourth moment (their X-term) at the phase-truncation limit
  N^{1/2}W² ≤ T^{3/2}; a saving T^{−η} there gives 37/30 − η/15, never below 16/13. Twelfth: their Thm 5.1 itself (middle
  and diagonal terms) and Heath-Brown's exponent pair (p_10,q_10); dropping the X-term entirely gives only 2 − 0.0011720.

## Before
Classical: μ(A) ≤ 1 + (A−4)/8 (Hardy–Littlewood/Ingham + Heath-Brown). DFP printed μ(6) ≤ 37/30, μ(12) ≤ 2 − 520/445121,
δ_k > 0 "by interpolation"; with Hölder: δ_7 = 0.01408, δ_8 = 0.01150, δ_9 = 0.00892, δ_10 = 0.00633, δ_11 = 0.00375.

## Method
Exact duality μ(A) = sup_s (A s + E*(s)) over the pointwise minimum E* of all large-value bounds; exact re-derivation and
LP re-optimisation of DFP's final exponent system (TeX §9) in W, N, k, G and an unused cap D ≤ Υ in Thm 5.1's E_2 sup;
exact rationals throughout (`DERIVATION.md`).

## Certificate
`uv run --no-project --with numpy python3 Riemann/Deep/DFPMoments/dfp_moments.py --all` → `ALL PASS` (982 checks,
≈ 85 s, one thread, peak 111 MB; log `logs/all.log`). Part 1 reproduces every printed exponent first (811 checks).

## Checks
- Reproduction: 37/30, 24/17 and 128/17 on [7/60, 5/38], E_k and σ_±(k) for k = 2..12, (1070 − 6164σ)/619 on
  [451/3688, 21/158], pair (64/135, 122/243), 593/288 and 3587/288, σ_2 = 58907/445121, 2 − 520/445121.
- Every curve vertex: all four terms and all constraints with slack; DFP's τ is the LP optimum (53 checks); Υ cap and
  out-of-range use give 0 gain (exact grids).
- Vertex certificate: at each maximising vertex all 13823 bounds ≥ E*; tails k > 61 and m > 120 by explicit inequalities.
- Sensitivity: a saving T^{−η} in Iwaniec's short fourth moment gives sixth = 37/30 − η/15 (η ≤ 1/26), floor 16/13 where
  the phase-truncation limit meets Thm 5.1's |η| ≤ cW/N²; the twelfth barely moves (1.9988280 even without that term).
- Guth–Maynard and Huxley via the approximate functional equation: worse than the fourth moment (T^{7/10}, T^{3/4} vs
  T^{1/2} at V = T^{1/8}).

## Novelty
arXiv listings Oct–Nov 2026, arXiv API, OpenAlex citations (only 2506.16766, log savings), HAL, Zenodo, ANTEDB (no DFP,
no commit since 26 Sep): not found. See `NOVELTY.md`.

## Caveat
Routine optimisation the authors explicitly skipped (TeX l.74); likely COROLLARIES-grade. Rests on DFP v1 (unrefereed,
partly AI-assisted per their note) as printed; we re-derived only their exponent bookkeeping. The sixth and twelfth
moment records are unchanged. "Optimal" means over the listed inputs and the chain's free parameters, not over methods.

RH relevance: N9

## Files
`Riemann/Deep/DFPMoments/`: `REPORT.md`, `DERIVATION.md`, `NOVELTY.md`, `dfp_moments.py`, `logs/all.log` (certificate),
`logs/novelty_*.log`.
