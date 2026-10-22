# Working notes: DBN6, Λ below 0.13963 under the 7/8 theorem (2026-10-21)

## Sources
- arXiv e-prints read from the TeX: Hakimoglu-Brown 2610.12803v1 (sha256 141d758a…), Polymath 1904.12438v2, and Planat 2610.18214 (sha256 4599b663…, 48 MB with its ancillary files).
- Raymond's deposit 23272725 was read from the copy already in .scratch/raymond (downloaded this morning; not touched).

## Hakimoglu-Brown's chain as a function of the wall X (evidence/hbmodel.py)
- Rate: k_i = 4{λC/p − μG − 2λ/(B−q) − 2μ(B+q)/(B−q)²}, with p = 61/20 and B = p².
- C and G come from Prop. "two-density": u = σ_p + (L/4)t⁻, u₈ = σ_p + (L/8)t⁻, the δ_j sums up to N₁ = 65536 plus the integral tail, ρ_j taken from N_* = ⌊√(X/4π)⌋, and η_j kept at the paper's values.
- (λ, μ) per step: 2-variable LP over the constraints Φ(D, r) ≥ 0, with D compactified and r ∈ [0, q_i].
- Gates are taken at equality and the final step is cut where Q = b.
- At the paper's X the model gives 0.13943 against the paper's 0.13963. The model is slightly optimistic: finer LP, no rounding.
- Results (evidence/hb_model.out):

| Case | Model bound |
|---|---|
| (a) dt = 1/400 | 0.13921 |
| (b) X = 6·10¹², b = 1/25 | 0.12408 |
| (c) Polymath Table 1 wall at 5·10¹²+194858, b = 0.16733² | 0.12350 |
| Raymond's line wall, b → 0 | 0.12179 |
| X = 6·10¹⁰ with b → 0 | 0.13759 |

- The cap comes from p = 61/20. Lowering p leaves the approximation band 3 ≤ y ≤ 6. There the comparison "F = ζ(s_B)(1+ε)" fails: at σ ≈ 1.2–1.5 the (log n)²-weighted zeta tail and the heat deformation make ε₂ of order 1.
- Polymath Table 1, row 2: X = 5·10¹²+194858, t₀ = 0.186, y₀ = 0.16733, winding number 0.
  - It is published (Res. Math. Sci. 2019) and covered by Platt–Trudgian's RH height. This is the input Platt–Trudgian use for Λ ≤ 0.2.
  - Its canopy condition (ii) is only "expected" in Polymath's text. The barrier (iii) is stated as computed.

## Raymond's chain (Zenodo 23272725)
- Rate: c = (2λ/h)(ℓ_{i,r} − M/h), on cells of 1/4000 with 32 substeps each. The seed is 9/16 + 10⁻⁶ and the chain ends at q_e = 10⁻⁸.
- Planat's ancillary SOURCE_RETUNED.json and SOURCE_OLD_ENHANCED.json are exactly the two providers, and they include the 57 kernels.
- Taking the best provider and kernel at each step reproduces 0.116397726408291 (evidence/rmodel_lib.py; inputs at DBN6_SRC).
- Kernel catalogue lever: a dense ρ grid (step 0.005) on all 10 shifts adds only 2.6·10⁻⁷ (rcat.py).
- Derivative lever: 2K_y ≥ λ'(y/h)J − μ'y∂_hJ, with the displacement held fixed.
  - Own-pair terms: M = max s_r(w) and g = min h²∂J_own over w ≤ ρ².
  - Rate: c' = (2λ'/h)(ℓ − M/h) − 2μ'(G − g/h²).
  - The far field forces λ' − μ' ≤ 2. The near field limits μ' at large y/h, so the lever adds little early (y/h ≈ 0.4; this agrees with Planat's 0.1% pilot) and more late (y/h ≤ 0.3).
  - Model bounds (rmodel_runs.py, rmodel_g.py):

| Assumed G | Model bound |
|---|---|
| ψ₂(σ_eff)/4 | 0.11233 |
| 1.2·ψ₂(σ_eff)/4 + 0.1 | 0.11325 |
| 1.5·ψ₂(σ_eff)/4 + 0.2 | 0.11434 |
| 2·ψ₂(σ_eff)/4 + 0.3 | 0.11546 |
| G_ideal | 0.11474 |
| G_ideal + 0.1 | 0.11505 |

  - Here G_ideal = ¼ Σ_{n≤690988} (D log_* b_t)(n) log n n^{−σ_B} (gideal.py). It is phase-free and covers the first branch only. It is much larger than ψ₂(σ_eff)/4 at small heights and late times, because large n see an effective σ = (1+h)/2 + (t/4)(L − log n) ≈ 1.2.

## Certifying G (not attempted)
- G must hold uniformly in x ≥ X₀ − 2, on every used cell (t ≤ 0.117, h = 0.69–1.88), for every frequency class and translated pair.
- Ingredients:
  1. Second-derivative error bounds on Raymond's Cauchy disks.
  2. Reflected-branch second derivatives: weights (L − log n)²/4 ≈ 45–180 on a mass of about 2·10⁻³.
  3. Second-order matched-convolution identities, (log n)² c_t(n) = Σ d̃(d)[log d + log(n/d)]c_t(n/d), with unmatched masses divided by |Q| ≥ 1 − Ā.
  4. Exact Bernstein certificates for the new (r, ρ, λ', μ') kernels and the own-pair bounds.
  5. An exact rational chain and an independent auditor.
- Estimated gain: 0.0013–0.0017 if G comes close to the ideal, and less once real certification losses are paid.

## Hygiene
- No commits. Only Riemann/Deep/DBN6 and .scratch/dbn6 were edited, and .scratch/dbn6 has been deleted. No orphaned processes.
