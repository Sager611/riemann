NEW: no source found states or implies (i) for almost all D, every primitive form of discriminant D represents a prime ≤ |D|^A with A < 1.1307 (let alone every A > 1), (ii) any average-over-D zero-density bound for class group L-functions (Σ_D N_D(σ) ≪ X^{2−σ+ε}), or (iii) a large sieve over (D, χ ∈ Cl(D)^) at the family size X^{3/2}. Caveats: GRH gives A = 1 for every D (Sardari 1802.06193); for a² + nb² alone Baier–Zhao 2007 is stronger; the pair-count device is Soundararajan's and Heath-Brown–Pierce's.

# Second adversarial novelty check of RCApps5 (2026-10-21; scratch `.scratch/rcnov5/`, deleted)

## Closest work, per claim
| Claim | Closest prior work (statements read here or in the [R3] audit, ../RCApps3/NOVELTY.md) | Verdict |
|---|---|---|
| (i) every form, almost all D, A > 1 | Unconditional, almost all D: Ditchen 1312.1502, \|D\|^{20/3+ε} (abstract: 7+ε) with density-zero exceptions, 3+ε for most forms, 2+ε under Lindelöf; Sardari 1803.03218 (JEMS 2021): only a positive proportion of classes. Every D: Kowalski–Michel 2002 Cor 7 and Michel's Park City notes (2007) Thm 3.11 (exponent unspecified), Zaman thesis 455, Thorner–Zaman 694, Gaudet 2408.01610 (sieve, L unspecified). Under GRH: Sardari 1802.06193 (every form represents a prime < \|D\|(log\|D\|)^{2+ε}), the A = 1 reached here without GRH for almost all D | new |
| corollary a² + nb² | Baier–Zhao, Math. Ann. 338 (2007): a prime m² + n < 4n for all but O(N(log N)^{−C}) squarefree n ≤ N | least-prime clause known in stronger form for this one form; power-saving exceptions, short intervals and other classes are new |
| (ii) density on average over D | nothing for class group L-functions. Blomer–Thorner 2410.17158 (TeX read): a density-hypothesis-type saving q^{−B(σ−½)} near ½, but for level-q GL(m) families (relative trace formula), no class-group content. General-family theorems (Kowalski–Michel 2002, Michel 2007, Humphries–Thorner 2103.05634, Pascadi–Thorner 2508.14888) give \|F\|^{c(1−σ)} with c too large near ½. Fixed-group families (Pierce–Turnage-Butterbaugh–Wood, Lemke Oliver–Thorner–Zaman, Lemke Oliver–Smith, An, Cho–Lemke Oliver–Zaman 2512.24963: TeX, exponent depends on a fixed G) exclude ring class fields (group of order 2h → ∞) | new |
| (iii) family large sieve | Ditchen: N log³N + N^{1/2}Q^{5/2+ε} (χ² ≠ 1); Pascadi–Thorner: N + Q^n\|S\| for any S; Michel 2007 (3.18): N + \|F\|²Q^d; Young 2208.03358: Eisenstein series of varying level; Humphries–Khan: one real quadratic field | none at X^{3/2} for N ≤ X^{1−ε}: new |
| method (Lemma P) | Soundararajan 2000, Heath-Brown–Pierce 2017 (pair counts summed over d); Fouvry–Iwaniec 2003 (GRH, low-lying zeros averaged over D via 4p = m² + Dn²) | device known; the positivity over all orders without conductor Möbius terms not found |

## Searches (2026-10-21)
- arXiv API: "Rate exceeded" on the first call, stopped. arXiv HTML search instead: 125 abstract, title and author queries; the four brief queries ("class group characters" + "large sieve"; "binary quadratic forms" + "least prime" + "almost all"; "ring class" + "zero density"; "dihedral forms" + "large sieve") return no paper at all.
- arXiv monthly listings: math.NT 2025-01 to 2026-10 (22 months, 8,333 titles, 100 keyword hits checked, abstracts read where relevant); math.GM and math.CV 2026-09/10 (627 titles). Nothing on primes in ideal classes for almost all D, on class-group densities over D, or on family sieves over D.
- Hexagon API: all 342 works (2609.00001–2610.00441, none new since the ChebST audit), titles and abstracts; relevant only 2610.00178 §54, 2610.00200, 2610.00216 (per-D or half-plane results, already in the NOTE).
- openai/math: 749 preprints, unchanged since 8 Oct; CONTENTS.md (373 entries) scanned: only the quasi-Riemann hypothesis (Dirichlet, Q(√−3)) and uniform exclusion of Siegel zeros.
- Zenodo: 15 queries, only unrelated items.
- Semantic Scholar citations: Ditchen 4, Zaman 1710.08914 5, Gaudet 0, Heath-Brown–Pierce 22, Humphries–Thorner 11, Thorner–Zaman 1604.01750 29, Kowalski–Michel 2002 67, Duke–Kowalski 2000 142, Sardari (two papers) 11. OpenAlex citations: Fouvry–Iwaniec 62, Blomer 2004 22, Baier–Zhao 25, Heath-Brown 1995 187, Humphries–Thorner 4; 10 topic searches and 4 title lookups. No citing work improves Ditchen or treats the family over D.
- Full texts searched with rg: Blomer–Thorner 2410.17158, Cho–Lemke Oliver–Zaman 2512.24963, Pierce 2206.08351 (TeX); Kowalski–Michel 2002 and Michel's Park City notes (PDF).
- GXL evidence API: HTTP 429 (daily quota used up).

## Residual risk
- Not opened: Duke–Kowalski 2000, Iwaniec–Li 2007, H. Kim, "Minimal generators of the ideal class group" (JNT 2021, no abstract online; by its title about generators, not every class), Y. Martin (JNT 2025, no abstract online). Journals outside arXiv were searched only through OpenAlex and Semantic Scholar metadata.
