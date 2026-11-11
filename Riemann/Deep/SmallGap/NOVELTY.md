# SmallGap: novelty check (11 Nov 2026)

Claim checked: an unconditional positive proportion of gaps γ_{n+1} − γ_n ≤ 2πλ/log T for λ < 0.7892 (here 0.7851).

| Source | Query | Result |
|---|---|---|
| arXiv listing math.NT 2026-11 (280 entries, ids parsed per entry) | titles: gap, spacing, pair correlation, zeros, critical line, density | only 2611.14486 (Lamzouri, the baseline); others off-topic (2611.06041 Poissonian pair correlation, 2611.12438 RH-proof claim, 2611.13129 Weil form) |
| arXiv listing math.NT 2026-10 (854 entries) | same | 2610.26038 (Bondarenko-Heap: μ ≤ 1/2 under RH, a liminf), 2610.28527 (Durkan-Pearce-Crump: unconditional liminf < 1/2 + ε for each Dirichlet L, positive proportion only under RH), 2610.22927 (moderate gaps ≥ 2πc/log T, large side), 2610.31744 (gap biases, numerics), 2610.08965 / 2610.23467 (proportions of simple / critical zeros), 2610.15989 (density in short strips, σ ≥ 7/10): none gives an unconditional small-gap proportion |
| arXiv API, sorted by date | abs:"small gaps" AND abs:zeta | 7 hits; newest relevant 2604.05733 (Inoue, RH, μ < 0.50895) |
| | abs:"small gaps" AND abs:zeros | 40 hits; only 2604.05733 relevant |
| | au:Lamzouri | newest 2611.14486v1 (no v2) |
| | abs:"Selberg minorant"; abs:semidefinite AND abs:zeta | 1702.04579; 1810.08843 (CGdL), 1908.04876 |
| | all:2611.14486 | only the paper itself |
| | mollifier AND "zero density"; "longer mollifier"; "zero-density" AND "critical line" | no density slope above 8/7 near σ = 1/2 |
| OpenAlex | works citing CGdL (W2897734907, 19) | newest: critical-zero proportions (Aug-Sep 2026), BGMM "Small gaps and small spacings between zeta zeros" (2023, RH), "Fourier optimization and Montgomery's pair correlation conjecture" (2024): none unconditional for D(λ) |
| Semantic Scholar | citations of 2611.14486, 1810.08843 | HTTP 429 (rate limit); 2611.14486 is one day old |
| Zenodo | "small gaps" (25 newest); "gaps between zeros" (8); Lamzouri gaps | 24095120 (conditional r-gaps), 22067122 (numerics), 24116425 (Weil forms); nothing on D(λ) |
| HAL | "small gaps" AND zeta | 0 |
| repo | rg small gap, 0.7892, 0.7864 in *.md; NOVELTY files | TREE R421 (records 0.7892, suggests this optimisation), Harvest6 candidate 2; GapLambda is the large-gap Λ record (other quantity). No prior small-gap result of ours |

Verdict: new unconditional record 0.7851 (was 0.7892 published; 0.7864 claimed in Lamzouri's remark l.199-201 without details). Not pre-empted as of 11 Nov 2026.
Significance: small (−0.0041 in λ) and anticipated by the remark ("by … optimizing their polynomial, one might obtain slightly better results"). New content beyond the number: the certificate, the black-box truncation (Lemma T), and the method's limit (float λ* ≈ 0.7849; no further gain without a density slope > 8/7 near σ = 1/2).
