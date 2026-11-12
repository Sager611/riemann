AUDIT PASSED — every step of Proposition N, Theorem N^R and Lemma M has a complete proof (FULLPROOF.md) and no step is false, so lincomb3's 2.061% (N = 3) and 34.707% (N = 2) stand on their quoted inputs, of which only Selberg's mean square of S(t+h,χ) − S(t,χ) could not be read at its source.

# LC3Audit — referee report on lincomb3 (DERIVATION §3–§5, §8, §9)

Referee: lc3audit (Orca), 12 Nov 2026, 13:41–14:25 UTC. Proofs: [FULLPROOF.md](FULLPROOF.md). Code/logs: lemmaM_hs.py, logs/lemmaM_hs.log.

## Short answer
- Proposition N, the capacity step, Theorem N and Theorem N^R: proved line by line (FULLPROOF §1–§2). Correct.
- Lemma M: proved (FULLPROOF §3, M1–M9). The merged multiset satisfies every hypothesis of 2610.08965's argument; the HS change is an exact
  identity bounded by 2κ²ΣΨ₁ = O_ε(κ²N); long chains cost ≤ 4C√κN. Two verdict-only reviewers, text pasted inline:
  Lemma M chain (M4–M8) CORRECT; Theorem N^R capacity-free (P2–P8, R2, R0) CORRECT.
- Numerics found nothing: the termwise centroid bound holds with ratio ≤ 7.3·10⁻⁴ (predicted 1/1152 by the leading term), never > 1.

## Findings
| # | kind | finding | effect |
|---|---|---|---|
| G1 | slip | DERIVATION 3.3 uses BH Lemma 9 "with δ/N" but keeps the √δ windows; Lemma 9's window is tied to its δ | none: use σ = √(δ/(N−1)) throughout (FULLPROOF P2, P4) |
| G2 | structure | with equal δ_j, 1 − Σ min(1 − δ, 1/N) = max(Nδ − (N−1), 0): the capacity step never helps | the record needs only the pairwise CLT (Lamzouri's Lem:SelbergCLT), not BH Thm B for N functions, nor 3.4/4.6 |
| G3 | precision | pair counts (P) hold uniformly in R ≥ 1 (FULLPROOF M6) | as DERIVATION intends |
| G4 | wording | the abstract Gram bound is a Lean theorem (GramEnergyStatement, any finite injective real family), resting on the native_decide axiom; Lemma M itself is outside Lean | accurate as stated in DERIVATION |
| G5 | inputs | Selberg 1946 / Fujii 1974 (mean square of S(t+h,χ) − S(t,χ)) unread: AMS browser check, no open copy | only unread load-bearing input; classical |
| G6 | numerics | at t ≤ 2400 the q = 7 dominance shares for unequal b are 0.11–0.65, far from 1/3 | asymptotic claim unaffected; not load-bearing (G2) |

## Inputs, as read
| input | source read | used for | load-bearing |
|---|---|---|---|
| bivariate Selberg CLT (Lem:SelbergCLT) | Lamzouri l.873–881 (cites BH Thm B + Selberg density, l.878–880) | separation set E, all pairs | yes |
| BH Lemma 9 (stability of zero counts, Eq:BH6.10) | Lamzouri l.898–908 | blocks, B_ℓ, k_ℓ | yes |
| N_F ~ N_T; N(2T,χ) − N(T,χ) = N_T + O(T) | Lamzouri l.856–861 | normalisation | yes |
| dyadic Garunkštis–Paliulionytė pair correlation | Lamzouri l.968–975 | energy tr A² for L(s, χ) | yes (Lemma M) |
| eq:key2, lem:threshold, certificate (x, b) = (8722, 4389)·10⁻⁶ | 2610.08965 l.517–535, 1072–1104, 1149–1235; Lean Spec.lean | Lemma M | yes |
| Selberg/Fujii mean square of S(t+h,χ) − S(t,χ) | not read (G5) | pair counts (P) | yes (Lemma M) |
| BH Thm B for N functions | not read (Duke 1995) | capacity 3.4 only | no (G2) |

## Numerical falsification (logs/lemmaM_hs.log; 91 s, peak 189 MB, one thread)
800 unfolded ζ zeros near t = 18 000 with planted tight pairs, spread clusters (diameter ≤ √κ), double points and long chains, κ = 1/16 … 0.001,
two windows (smooth bump; 2610.08965's second window, edge-smoothed). Exact HS norms on the Fourier side (x-space cross-check agrees to ≤ 1.6·10⁻⁸ relative).
- Termwise |tr(D_C D_C')| ≤ 2ψ(Δ)m_C m_C' κ'⁴: max ratio 6.7–7.3·10⁻⁴, all runs. Sums: E2 ≤ B2 ≤ B2'. ‖A'−A‖_HS/(κ√N) ≤ 0.90 vs proved constant ≥ 17.
- Centroid merging is a refinement only: endpoint merging (E2 up to 2.4·10⁴ times larger) obeys its first-order bound, O(√κ√N) under (P), which also suffices (F3).
- Unplanted zeros at κ = 0.2, 0.3 (7 and 15 natural clusters): same bounds hold.
- Simple points of Z' = Y* ⊔ L exactly; |L| ≤ 4√κ·#pairs; triangle inequality for √tr A²; real zeros alone: no cluster at κ ≤ 1/16.
- q = 7 logs (not rerun): the per-block inequality holds in all 304 strongly dominated blocks, with isolation measured from on-line zeros (the Y* notion).

## What survives if an input falls
| input fails | N = 3 bound left | why |
|---|---|---|
| none | 3s₀ − 2 = 0.0206146 | this audit |
| 2610.08965 certificate, or (P) | 3p − 2 = 0.0174739 | Theorem N with IsoZero's c-isolated p = 0.672491294 (no Lemma M) |
| dyadic GaPa (energy for L(s, χ)) | Dousselin's 4.91·10⁻⁹ only | both s₀ and p use it |
| BH Thm B for N functions | unchanged | G2 |
| pairwise CLT or BH Lemma 9 | Dousselin's 4.91·10⁻⁹ only | the dominance mechanism needs both |

By-product (Lemma M for ζ and fixed L(s, χ)): ≥ s₀ − η simple zeros on the line isolated from the other on-line zeros at a non-explicit κ(η).
It does not supersede discoveries/isolated-zeros.md (0.6724, isolation from all zeros, explicit c ≤ 0.00476): different notion, neither implies the other.

RH relevance: N8, N12
