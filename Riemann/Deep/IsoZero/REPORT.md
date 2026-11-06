RECORD: at least 0.6724 of the zeta zeros are c-isolated for c = 0.004 (unconditional)

## Claim
For every c ≤ 0.00476, liminf_{T→∞} B_0(T,c)/N(T) ≥ 0.6724; at c = 0.001 the certified bound is 0.672487. B_0(T,c) counts zeros whose two neighbouring ordinates are more than 2πc/log T away (Lamzouri's Thm 1.1); such zeros are simple and on the critical line. Zeros counted with multiplicity.
More explicit constants (all certified): ≥ 0.672 for c ≤ 0.01106 (so for c = 1/100: 0.672090); ≥ 0.67 for c ≤ 0.02485; ≥ 0.665 for c ≤ 0.04282; ≥ 0.66 for c ≤ 0.05495; Lamzouri's 0.6552, for which he gives no c, for every c ≤ 0.06427; at c = 1/6: ≥ 0.5247 (paper: 0.50875).
So, up to 1.4e-5, every simple zero on the line that Lamzouri's Hilbert-space method produces (proportion 2 − C_MT = 0.672501, his 2609.02882) is 1/1000-isolated.
Key new lemma (DERIVATION §3.7): for every fixed ρ > 0, Σ_{0<γ≤T, β>1/2} ((β − ½)(log T)/2π − ρ)^+ ≤ (B(ρ) + o(1)) N(T), B(ρ) = (1/4π) log(1 + (1 − e^{−4πρ})/(e^{16πρ/7} − 1)); B(1.2) ≤ 1.441e-5. It comes from Littlewood's lemma at σ = ½ + 2πρ/log T with Conrey's published mollified mean value (Crelle 1989, Theorem 2, moved right of the line by the functional equation). Off-line zeros below normalised height 1.22 cannot spoil isolation (they are charged like double zeros); higher ones are exponentially rare. Lamzouri instead used only Σ(β − ½) ≤ 0.0806·T, an announced, never-printed constant of Conrey, against a linear deficiency bound.

## Before
- Lamzouri, arXiv 2610.26841v2 (26 Oct 2026, latest version): ≥ 0.6552 for some non-explicit c; explicit only 0.50875 at c = 1/6.
- Lee, Zenodo 23719902 / 23807772 v1.1 (28/30 Oct 2026, unrefereed): explicit 0.655208 at c = 0.00577, 0.655200 at c = 0.0059665 (Conrey's 0.0806, forty-mode kernel). Lee, Zenodo 23922845 v3.1 (3 Nov), conditional on his local-block hypotheses: 0.5329 at c = 1/6, 0.65155 at c = 0.05. Lee, Zenodo 23688612: > 1/2 for c = 0.17661804.
- Angón–Ross–Tang–Zhao, arXiv 2610.22927: ≥ 0.3450 for every c ≤ 1/4.
- The first explicit c for the 65.5% statement is Lee's, not ours: with the paper's own inputs (column A) we get 0.655 at c = 0.00801 and 0.6552 at c = 0.00096, weaker than Lee.
- Our N8 record for simple zeros on the line (67.353%, 2610.08965) is unchanged; isolated zeros are a subset.

## Method
1. Pro:IsolatedSMALL (TeX l.511–573) with the O(τ²N) made explicit (DERIVATION §2): cluster loss λ = max(2(1 − K(2c)), β_T) ≈ 12.2c² per clustered point, charged only on the clustered fraction; the witness error θ̄ = π²c²/2 folded into the deficiency H_τ(y) = 1 − (1 − θ̄)²·4K(iy)²/(1 + K(2iy)); one inequality for all c, reproducing both parts of Thm 1.1 (c → 0, and c = τ_0 = 0.169642).
2. Hinge lemma (above, DERIVATION §3.2, §3.3, §3.7); certified ρ_max(c) with H_τ − 3λ < 0 on (0, ρ_max] and α_ρ = sup (H_τ − 3λ)^+/(y − ρ); off-line term α_ρ B(ρ) ≤ 9.5e-6 at c = 0.001.
3. Result: p(c) = (2 − λ − C_MT − α_ρ B(ρ))/(1 − λ); kernel cos(wx) optimised per c (only helps for c ≥ 0.05).
4. By-products: Lamzouri's Lem:Conrey (0.0806) is (1/4π)log(11/4) = 0.0805006 rounded up, and follows from Conrey's Theorem 2 by the same argument with ρ → 0 (DERIVATION §3.4); the intermediate first-moment routes (columns B, C in logs/cert.log) are superseded.

## Certificate
cd Riemann/Deep/IsoZero && uv run --no-project --with python-flint --with mpmath python3 cert.py
(Arb, 160 bits, one thread, about 6 minutes, under 100 MB; output logs/cert.log. Every bound is rounded in the safe direction, including α_0 ∈ [0.214612853119, 0.214612853122] (the paper's Mathematica step) and C_MT = 1.3274992963… (closed form and Arb quadrature).)
Independent float check of the headline (mpmath, kernel by quadrature): uv run --no-project --with mpmath python3 sanity_mp.py → ρ_max 1.22365, α_ρ 0.6701238, B 1.40382e-5, p 0.6724873 (logs/sanity.log).

## Table (liminf B_0(T,c)/N(T) ≥)
A: the paper's inputs (Σ(β − ½) ≤ 0.0806T), Montgomery–Taylor kernel f_0, our explicit error terms. H: new (hinge lemma), f_0. H+w: new, cos(wx) kernel optimised (w² in brackets).

| c | A | H | H+w |
|---|---|---|---|
| 0.001 | 0.655200 | 0.672487 | 0.672487 (2.0020) |
| 0.002 | 0.655190 | 0.672475 | |
| 0.003 | 0.655175 | 0.672455 | |
| 0.005 | 0.655124 | 0.672391 | |
| 0.0075 | 0.655025 | 0.672266 | |
| 0.01 | 0.654887 | 0.672090 | 0.672090 (2.0020) |
| 0.015 | 0.654491 | 0.671588 | |
| 0.02 | 0.653934 | 0.670882 | |
| 0.03 | 0.652329 | 0.668852 | |
| 0.05 | 0.647037 | 0.662211 | 0.662237 (2.1257) |
| 0.075 | 0.635983 | 0.648565 | |
| 0.1 | 0.618700 | 0.627814 | 0.628354 (2.5315) |
| 0.125 | 0.592766 | 0.597867 | |
| 0.15 | 0.553786 | 0.555171 | |
| 0.16 | 0.532872 | 0.533258 | |
| 1/6 | 0.516640 | 0.516706 | 0.524702 (3.7500) |

Largest c (to 1e-5) with the target, column H: 0.6724 → 0.00476; 0.672 → 0.01106; 0.67 → 0.02485; 0.665 → 0.04282; 0.66 → 0.05495; 0.6552 → 0.06427. Column A: 0.655 → 0.00801; 0.6552 → 0.00096. (B_0(T,c) decreases in c, so each value holds for all smaller c.)

## Checks
- Reproduces the paper: c → 0 with the paper's inputs gives 2 − C_MT − 0.0806α_0 = 0.655202908 (l.816 claims > 0.655202; margin 9e-7); λ = 1/3 (c = τ_0 = 0.169641840) gives 5/2 − 3/2·C_MT = 0.508751056 (Pro:IsolatedZeros, (1.6)); H_0(2) = 0.4092557 > 0.4092.
- τ–c relation: Z_T elements are γ log T/2π, so separation τ in Z_T is 2πτ/log T between ordinates; the proposition is applied with τ = c; B_0(T,c) ≥ |I_1| − 2.
- Hidden c-dependence: none beyond λ, θ̄, ρ_max and the side conditions 2K(c)² ≥ 1, 4c ≤ 1, λ < 1 (asserted in cert.py); the pair-correlation constant and B(ρ) do not depend on c.
- liminf versus finite T: all statements are liminf; the o(1) terms are not explicit; no finite-T claim.
- Multiplicity: as in the paper; a multiple or off-line witness is charged by its multiplicity; sums over distinct off-line pairs are bounded by sums over zeros with multiplicity.
- Limits ε → 0 (smoothing), θ → 4/7, P → sinh(bx)/sinh b: each admissible choice gives a valid liminf bound, and the bound is continuous (DERIVATION §4).
- Independent verdicts (verdict-only reviewer runs on inline summaries, logs/review.log): first-moment lemma SOUND; Hilbert-space budgets SOUND; hinge lemma: steps SOUND, one prose rounding (B(1.2) quoted as ≤ 1.44e-5, true 1.4406e-5) fixed; the certificate was unaffected.
- AUDIT.md: no flaw in 2610.26841v2; no input needs RH or Montgomery's conjecture; the one unproved input (Lem:Conrey) now has a proof.

## Novelty
- arXiv listings math.NT and math.CA, 2026-10 and 2026-11 (lst.py, 29 title hits): only 2610.26841 and 2610.22927 give isolated-zero proportions; nothing on off-line moments.
- Zenodo '"isolated zeros"': Lee's four records (above); the three at the 65.5% level use Conrey's 0.0806 with the linear deficiency bound, the fourth is the c = 1/6 route. HAL: nothing relevant. Lamzouri: v2 is the latest.
- GXL evidence API, two claims on the constant in Σ(β − ½) ≪ T and its improvements: 0 supporting, 0 refuting of 256 passages each; arXiv web search: nothing.
- Our files: Riemann/Deep/ZetaZeros, N8Opt, N8Fam treat simple and on-line proportions, not isolation; R122 and R318 record 2610.26841 only.
- New here: the hinge lemma in this explicit form and its use against the deficiency, giving 0.6724 (from 0.6552); the explicit c-table; the unconditional 0.5247 at c = 1/6. Risk: Littlewood's lemma off the line is classical (Selberg; Conrey's density N(σ,T) ≪ T^{1−8(σ−½)/7} log T); an explicit integrated bound of this kind may exist in older literature outside the indexed sources. The isolated-zero application is new as far as these searches reach.

## Caveat
- Same footing as the paper: Lem:ConstructionKernel is Lamzouri's 2609.02882 Lemma 3.2 (2026 preprint) on top of Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh (Acta Arith. 2024). Conrey's Theorem 2 (case Q = 1, fixed R > 0, θ < 4/7) is used as quoted by Bui–Conrey–Young (arXiv 1002.4127, l.338–345); Conrey's 1989 paper was not read here. Our own derivation of the shifted main term (DERIVATION §3.1) agrees with that formula.
- liminf only. The reviewer verdicts check summaries, not the full text; DERIVATION §2–3 should go to a specialist reviewer.
- The bound cannot pass 2 − C_MT = 0.672501 with this method (the pair-correlation energy is at least C_MT for every admissible kernel, l.222–227 and CCLM); the simple-and-on-line record 0.67353 uses a different method.

## Side result: linear combinations of two Dirichlet L-functions (Lamzouri's Thm 1.4)
Claim. Let χ_1 ≠ χ_2 be primitive characters mod q ≥ 3 of the same parity, b_1, b_2 nonzero reals, e^{2iθ_j} = i^a τ̄(χ_j)/√q. At least 0.3449 of the nontrivial zeros of F(s) = b_1e^{iθ_1}L(s,χ_1) + b_2e^{iθ_2}L(s,χ_2) lie on the critical line and have odd multiplicity (certified 2p − 1 ≥ 0.344982508).
Before. Unconditionally: Selberg (1998) a positive proportion; Dousselin (thesis) 7.38·10⁻⁹; Lamzouri 2610.26841v2, Thm 1.4: 31% (= 2(0.6552 − ½)). Under GRH alone, Bombieri–Hejhal (Duke 1995, Theorem C): 34.5%; under GRH plus a well-spacing hypothesis: 100%. No other 2026 paper found (lst.py 'linear combination|Epstein|Davenport.Heilbronn|Bombieri.Hejhal|Selberg class', math.NT/CA/CV, 2026-10/11; Zenodo).
So the unconditional bound now matches Bombieri–Hejhal's GRH-conditional 34.5% to within 2e-5 (2(2 − C_MT) − 1 = 0.345001).
Method. Lamzouri's Pro:BombieriHejhal (audited, AUDIT.md, TeX l.824–964) gives the proportion δ_1 + δ_2 − 1, where δ_j is the proportion of c-isolated zeros of L(s,χ_j) in dyadic blocks (a fixed small c > 0, as the paper requires; a smaller c only increases δ_j). δ_j = p(c) of the headline, because every ingredient transfers with the same constants: pair correlation (Lem:ConstructionKernelDirichlet, from Garunkštis–Paliulionytė and La), and the hinge lemma for L(s,χ) from Wu's mollified mean value (Math. Z. 293 (2019), arXiv 1802.09704, Proposition 3.1, uniform in complex shifts a, b ≪ 1, θ < 4/7 for mollifiers μ(n)P(log(y/n)/log y), P(0) = 0, P(1) = 1), applied directly at the shift +2πρ/log T (DERIVATION §8).
Certificate. cd Riemann/Deep/IsoZero && uv run --no-project --with python-flint python3 cert8.py (Arb, about 5 s, logs/cert8.log): c = 10⁻⁴: p ≥ 0.672491254, 2p − 1 ≥ 0.344982508; c = 0.001: 0.344974576; the paper's own inputs reproduce 0.310405815.
Caveat. Rests on Bombieri–Hejhal's Theorem B and Lemmas 7, 9 as cited by Lamzouri (not re-derived; they need only Selberg's zero density for fixed χ), on La's preprint (as the paper does), and on Wu's Proposition 3.1 as stated (two harmless typos noted in DERIVATION §8.2: (h,k)^{α+β} for (h,k)^{1+α+β}, and the sign of R in his (1.6)). The evaluation of Wu's main term at the shift is the ζ computation with Euler factors at p | q removed; it was checked structurally (and at α = β = 0 against his value 1 + 1/θ), not symbolically term by term. Fixed q and χ only (no uniformity in q). Reviewer (verdict-only, pasted chain, logs/review.log): first pass GAP — the per-character transfer and the order of limits were not stated; added as DERIVATION §8.4; second pass SOUND (logical closure from the stated premises; it did not re-check the cited papers or rerun the certificate).

RH relevance: N8

## Files
- Riemann/Deep/IsoZero/REPORT.md (this), AUDIT.md (chain and verdict), DERIVATION.md (every step and constant), cert.py (certificate), cert8.py (side result), sanity_mp.py (float cross-check), logs/cert.log, logs/cert8.log, logs/sanity.log, logs/novelty.log (searches), logs/review.log (reviewer verdicts).
