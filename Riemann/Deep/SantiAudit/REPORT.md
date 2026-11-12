AUDIT PASSED — Santibáñez-Leal's Thm 1.1 (Zenodo 23134787 v0.03: positive proportion of simple critical zeros in (T, T+T^θ] for every θ ∈ [0.527,1), h(0.527;κ) > 0.0005947542001) holds as written, unconditionally, given its exact attributed sources (Wang arXiv 2609.07918, Estermann FE, Young's main-term evaluation); no error found; its new off-diagonal estimate is correct but lossy.

# SantiAudit — report (orca agent santiaudit, 12 Nov 2026)

## Short answer
- Chain: Thm 1.2 (short-window mollified moment, ν < νA = min{1/2, 2θ−1, (5θ−2)/6}) → Thm 1.3 (Levinson count, O/N ≥ κ = 1 − R⁻¹log c) → Thm 5.1 (parity product (Q−S)(N−O) ≥ 2(N−S)² + Wang's pair sum S_K/N → C(f)) → S/N ≥ h(θ;κ) = (3+κ−√((1−κ)(9−κ−8c(θ))))/4.
- θ = 0.527 is reached without c(θ) > 0: Wang's pair-sum asymptotic is proved for every 0 < λ < θ < 1 (also in v2); positivity needs only κ > 1 − 2/(2−c(θ)) = 0.0349241 at θ = 0.527, and the detector gives κ = 0.0357950. No input is used outside its proved range.
- Recomputed in Arb (exact rationals, 4000 bits): κ = 0.03579499548791817334727…, h(0.527;κ) = 0.000594754200102962… ✓.
- Lead E1 (ours, not yet gated): bounding the paper's own exact dual sum absolutely (Hankel form) gives ν < (3θ−1)/4 and, with κ1 = 0.0844779 at ν = 3/25, positivity for **every fixed θ > 1/2** (h ≥ 0.0055353).

## What was checked (details: FULLPROOF.md)
| Step | Status | How |
|---|---|---|
| (10) Gaussian twisted moment = divisor sum with kernel K_β | ✓ | re-derived; numeric (mpmath, T=40, H=6): rel. 1.4e−11, 5.0e−10, 2.0e−12 for (h,k) = (1,1),(2,1),(1,2) |
| (12) composite (Estermann) FE, (14) dual Mellin form, (15) multiplier | ✓ | re-derived from classical Estermann FE (σ_{a,b}(n) = n^{−b}σ_{b−a}(n)) |
| Lemma 3.1, (17) → (18): relative errors TM/H², TM³/H^{5/2} | ✓ | re-derived; reviewer CORRECT; exponents −51/12500, −6711/40000 at θ−η = 0.52699 |
| Residues = Young's diagonal main term (ν < 1/2); compact windows (§3.7) | ✓ | Fourier-deficit identity for w_D checked |
| Sec. 2 (Young localized, ν < θ−1/2) | ✓ | not sufficient alone at θ = 0.527 (needs ν < 0.027) |
| Lemma 4.1, (20), Lemma 4.2, Thm 1.3 (short-rectangle Littlewood, horizontal O(L), concavity) | ✓ | sign of Δφ, multiplicity bookkeeping, O(L²) count error vs N ≍ HL |
| Parity product [4, Prop 2.5] for conjugation-invariant multisets (off-line zeros allowed) | ✓ | Lamzouri-type bound Q ≥ 2tN − (2t−1)S − t²d re-derived (Bessel + Gram–Schmidt ranges); 120k adversarial multisets: no violation |
| Wang Thm 2.2 + (eq:unweighted): S_K(I) = (C(f)+o(1))HL/2π, 0<λ<θ<1 | ✓ statement/range | v1 and v2 (7 Oct) TeX; errors O(H + xL³), x ≤ T^λ; not refereed line by line (attributed preprint) |
| Detector (degree 201, Q(y)+Q(1−y) ≡ 1) and all 8 Table-1 κ values | ✓ | arb_check.py; Table 1 reproduced |

## Falsification: where each key step stops
| Step | Fails at |
|---|---|
| Sec. 2 diagonal range | ν ≥ θ−1/2; alone it proves θ ≥ 0.5339856 (ν = 0.0339 detector) |
| (18) with ν = 0.0499 | θ ≤ 0.52495 |
| h(θ; 0.0357950) > 0 | θ ≤ 0.5264286 (binding bottom constraint; paper rounds to 0.527) |
| best possible with νA and κ ≈ 0.71733ν | θ ≈ 0.52573 |
| hidden RH / pair-correlation conjecture / simplicity | none used |
Minor: Thm 5.1 writes η ∈ C_c (Wang needs C_c^∞; density argument); cites "[14, Thm 2.2]" for what is Wang's unweighted lemma.

## What survives
Everything: Thm 1.1 exactly as written. Range by moment route (honest tradeoffs):
| Moment input | Least θ proved | Status |
|---|---|---|
| Sec. 2 only (classical) | 0.5339856 | proved (paper) |
| Sec. 3 νA (paper) | 0.5264286 (stated 0.527) | proved (paper, audited) |
| absolute Hankel bound, ν < (3θ−1)/4 | every θ > 1/2 | Lead E1: sketch + reviewer CORRECT on the key bound; needs full write-up, second check, novelty gate |

## Lead E1 (for the orchestrator)
Only §3.5 changes: per (g,p,q) the dual is ≪ T^ε g⁻¹(T/H)^{1/2} (X ≍ √(nT/pq), integration by parts against the width-H bump), total relative error T^{1/2+2ν−3θ/2+ε}. Steuding (2002) had this line only with linear Q (θ ≥ 0.552); combined with the degree-201 detector and the parity transfer it reaches every θ > 1/2, the natural barrier of the method. Constants: nu12.log, lead_h.log. Candidate N8 threshold result: run the AGENTS.md gate (own proof write-up, reviewer, arXiv/Zenodo novelty) before any DISCOVERIES entry.

## Files
FULLPROOF.md (referee notes); arb_check.py/.log, arb_all.log (detector constants); onsets.py/.log; parity_fals.py/.log; check10.py/.log; nu12.py/.log, lead_h.py/.log; params_exp010.json (two frozen detector sets, CC BY 4.0, Zenodo 23135354). Peaks ≤ 25 MB.

RH relevance: N8
