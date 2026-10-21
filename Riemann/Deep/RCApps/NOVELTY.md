NEW: nothing found states or implies Theorems 2–3 (every primitive form of discriminant D represents a prime ≤ 2|D|^A, any A > 1.6536, for all but O(X^{1−c}) fundamental D ≤ X; also in short intervals and sectors). The almost-all record is still Ditchen 2013 (|D|^{20/3+ε}, density-zero exceptions). Exception: Corollary 4's least-prime clause is already known in a stronger form (Baier–Zhao 2007: a prime m² + n < 4n for density-one squarefree n).

# Adversarial novelty check of RCApps (2026-10-21; raw data in `.scratch/rcnov/`)

| Source (statement read) | What it proves | Compared with RCApps |
|---|---|---|
| Ditchen, arXiv 1312.1502, §1 Corollary (same as his ETH thesis, 2013) | fundamental q ≢ 0 (mod 8): max_C p(q;C) ≤ \|q\|^{20/3+ε} outside a density-zero set (BV with log-power savings); ≤ \|q\|^{3+ε} for 100% of classes; 2+ε under Lindelöf for weight-one Rankin–Selberg | RCApps: 1.6536 < 2 unconditionally, every class, all fundamental D, power-saving exceptions |
| Zaman thesis 2017, (1.36); Gaudet thesis 2023, p. 5; Gaudet arXiv 2408.01610, §1 | all quote Ditchen's 20/3 (density zero) as the almost-all bound | confirms the record as of 2024 |
| Every D: Zaman thesis Cor 1.4.1 (\|D\|^{455}); Kowalski–Michel 2002 Cor 7 (n^{L₂}, L₂ unspecified); Thorner–Zhang 2412.01802 (constants unspecified) | Linnik for each D | weaker exponents |
| 2026, unrefereed: Hu, Hexagon 2610.00216, Cor. "Effective Chebotarev" (q_𝔭 ≪ 𝒟_L⁸); Frei–Sofos–Widmer 2610.18491 (x ≥ (log Δ_L)^{1/δ}); Hexagon 2610.00178 §54 | for ring class fields log Δ_L ≈ h log\|D\|, so \|D\|^{4+ε} for every D; §54 leaves A unspecified | weaker |
| Baier–Zhao, Math. Ann. 338 (2007), Thm 1 and Cor. with y = x²; extended in IJNT 5 (2009) | all squarefree k ≤ y except O(y(log y)^{−C}) have Σ_{n≤√y} Λ(n²+k) ~ 𝔖(k)√y, so there is a prime n² + k ≤ 2y: the least prime a² + kb² is < 4k for density-one k | beats Cor 4's 2(4n)^{1.7}. Cor 4 keeps only its O(N^{1−10⁻⁴}) exceptions, all-x short intervals and sectors. The REPORT row with baseline "n^{L₂} (KM)" should cite this instead (Ditchen §10 states n^{2+ε} under an unnecessary assumption) |
| Sardari, JEMS 23 (2021); Diao 2506.18065; Browning–Sofos–Teräväinen 2212.10373 | a positive proportion of classes has a prime ≤ h(D)log\|D\| (all D); random forms or polynomials in boxes, log-power exceptions | never every class |
| Fixed-group families: Pierce–Turnage-Butterbaugh–Wood 1709.09637, An 2004.14410, Thorner–Zaman 1909.01338, Lemke Oliver–Thorner–Zaman 2012.14422, Lemke Oliver–Smith 2405.08383 | Chebotarev for almost all G-fields, G fixed | ring class fields have degree 2h → ∞, so not covered |
| Tools with no forms or class-group applications (TeX searched with rg): Humphries–Thorner 2103.05634 and 2004.09488, Pascadi–Thorner 2508.14888, Thorner–Zaman 1906.07717, Brumley–Thorner–Zaman 1804.06402, Lemke Oliver–Thorner 1505.03122, Blomer–Thorner 2410.17158, Thorner 2608.12257 | density, large-sieve and zero-free-region tools | no least-prime corollary for forms |
| Sectors: Khale et al. 2008.09677, Panidapu–Thorner 2105.11093, Huang–Liu–Rudnick 1903.04005, Järviniemi–Teräväinen 2303.05822, Duke 1989 | a fixed field (often Q(i)) | not uniform in D |

Searches (2026-10-21):
- GXL evidence API, 6 claims (`gxl/`). Hits: only Ditchen, Zaman 2018 (upper bounds only), Gaudet (quotes Ditchen), Sardari and the Gaussian-sector papers.
- arXiv API: 65 new queries, id lists and the 300 newest math.NT papers (10–20 Oct 2026) (`ax/`). 22 TeX sources searched with rg (`src/`).
- Hexagon API: all 323 works, 2609.00001–2610.00419 (none new since 20 Oct); TeX of the 102 math.NT works. Zenodo API: 32 queries, 15 Aug–31 Oct 2026. openai/math: 749 preprints (none new) and 7 code searches. zbMATH API. Semantic Scholar and OpenAlex citations of Ditchen (4) and of Baier–Zhao (38). Nothing relevant in any of these.
- Residual risk: Ditchen's ETH thesis PDF was not opened (HTTP 429), but Gaudet 2023 cites it for 20/3. Journals outside arXiv, Zenodo and Hexagon were searched only through OpenAlex and zbMATH metadata.

Verdict: Theorems 2–3 and the short-interval statement are new. The exponent beats Ditchen's 20/3, and even his Lindelöf-conditional 2. The power-saving exceptional set is also new: before, the best bound with power-saving exceptions was the every-D bound |D|^{455}. The REPORT/NOTE comparison for Corollary 4 must cite Baier–Zhao.
