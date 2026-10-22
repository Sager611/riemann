NEW: no source found states or implies, for any totally real F ≠ Q, that almost all CM K/F have a degree-one prime in every ideal class below d_K^A for some A < 1.1307 (let alone every A > 1), nor Theorem 1 (sieve over (K, χ ∈ Ĉl_K) at X^{3/2}), Theorem 2 (Σ_K N_K(σ) ≪ X^{2−σ+ε}) or a relative Lemma P. Not new: the sharpness clause. It is the key lemma of Miller–Peckner (J. Number Theory 132 (2012), arXiv 1003.5336): every degree-one principal prime of a CM field K has norm ≥ C(F)·d_K, with C(F) = 1/(4ⁿd_F²) in the case their proof writes out.

# Third adversarial novelty check of CM3 (2026-10-21/22; scratch `.scratch/cmnov3/`, deleted)
Scope: only what [CMQuartic NOVELTY](../CMQuartic/NOVELTY.md), [RCApps5 NOVELTY2](../RCApps5/NOVELTY2.md) and [NOTE](NOTE.md) §7 did not cover. Their tables stand. [CMQ] = CMQuartic (NOTE and NOVELTY), [R5] = RCApps5.

## Closest work found in this check (statements read)
| Claim | Closest work | Verdict |
|---|---|---|
| Thm 3: almost all K/F, every class, every A > 1 | Miller–Peckner 2012: the same family (all CM K with K⁺ = K₀ fixed, h(K₀) = 1), but under GRH and only the 1-level density of {L(s,χ) : χ ∈ Ĉl_K} for each K, supp φ̂ ⊂ (−1,1). No primes in classes, no zero density. Kim, JNT 222 (2021): off a density-zero subset of "some family", Cl_K is *generated* by degree-one primes of norm O((log d_K)^α), α > 4 (generation, not every class). Blomer–Brumley–Radziwiłł 2603.05609: many small split primes for all but O((log log X)^{1+o(1)}) discriminants ≤ X (F = Q; split primes, no classes). Xie 2606.30567, improving Deshouillers–Gun–Ramaré–Sivaraman: every narrow ray class mod 𝔮 of a fixed field is a product of three primes of norm ≤ (N𝔮)^{103/64+κ} | new |
| sharpness (exponent 1) | Miller–Peckner, key lemma, display (2.10): p ≥ Δ/(4^N D_{K⁺}²) for degree-one principal primes. Their proof is written for O_K = O_{K⁺}[√β] only ("the other case is similar"); [CMQ] Remark 2.3 covers every K | known |
| Thm 1: sieve over (K, χ ∈ Ĉl_K), K/F CM | Fixed field only: Duke 1989 (large sieve for Größencharaktere; hybrid zero density in the torsion-free aspect) and Schumer 1986. Varying field: only F = Q (Ditchen) and the general GL_n sieves already compared in [CMQ] and [R5]. "class group characters" + "large sieve": 0 hits in OpenAlex and Zenodo | new |
| Thm 2: density on average over K/F | Besides the general family densities compared in [CMQ] (Pascadi–Thorner, Humphries–Thorner), the only averaged zero statistics found are 1-level densities: Miller–Peckner (GRH) and quadratic Hecke families over a fixed imaginary quadratic field (Gao–Zhao 1708.01701, 1711.10557) | new |
| Lemma P_F | no relative version of the Soundararajan / Heath-Brown–Pierce pair count found (GXL claim: 0 supports; forward citations of HBP17 read) | not found |

## 2026 items that bear on CM3 without anticipating it
- Frei–Kai, Hexagon 2610.00466 (22 Oct): for a fixed number field F, (1 − β)·log N𝔮 ≥ ε_F for every real zero of every real Hecke character of F. Nothing about classes or averages. Applied to the genus characters, it would remove the Siegel-zero half of the inherited rough-𝔡 restriction ([CMQ] §7.1). It does not supply the missing 𝒬-aspect log-free density over F.
- Frei–Sofos–Widmer 2610.18491 and Hu 2610.00216: per-field bounds, about d_K^{4+ε} on their own and d_K^{8/5+ε} with [CMQ] Theorem A (already in [CMQ] NOVELTY).

## Searches (2026-10-21/22)
- arXiv math.NT listings 2025-01 … 2026-10: 8,333 entries, 800 title hits. All 8,333 abstracts were scanned with CM/relative × least-prime/zero-density/large-sieve/Chebotarev/almost-all patterns: 8,260 via the Semantic Scholar batch API (the arXiv API answered 429 to the first id_list call) and the other 73 by one arXiv id_list call on 22 Oct. About 40 abstracts read. 2026-10 re-fetched on 22 Oct: still 625 entries. arXiv API search, 8 abs: queries: 2 known items.
- Hexagon API: 381 works, 2609.00001–2610.00483 (39 newer than the brief's range); titles and abstracts scanned; TeX of the 111 math.NT works searched with rg. Relevant: only Frei–Kai (Kominers 2610.00451: weighted Chebotarev for one fixed field).
- openai/math: 749 preprints, last manuscript commit 8 Oct. CONTENTS.md (9,142 lines): no ideal-class prime, class-group density or CM-family result. Zenodo API: 28 queries, nothing.
- Forward citations (Semantic Scholar + OpenAlex) of Thorner–Zaman 1604.01750, Zaman 1508.00287 and 1502.05679, Lemke Oliver–Thorner 1505.03122, Lemke Oliver–Thorner–Zaman 2012.14422, Pierce–Turnage-Butterbaugh–Wood 1709.09637, Kowalski–Michel 2002, Ditchen 1312.1502 and Heath-Brown–Pierce 2017: 321 distinct citing works, 39 kept by a CM/relative × family/least-prime filter, abstracts read. Also on Semantic Scholar: Habegger–Pazuki (11 citing works), Kim 2021 (2), Sardari 1802.06193 (5), Cho–Kim 2014 (5), Fiorilli–Jouve 2001.05428 (15), Miller–Peckner (29). Nothing relevant.
- OpenAlex: 24 phrase queries. Crossref: 16. zbMATH: 14. Semantic Scholar relevance search: 8 (6 answered). GXL evidence API: 9 claims. The only supporting passages come from per-field or GRH-conditional work (Sardari, Zaman, Thorner–Zaman, Deshouillers–Gun–Ramaré–Sivaraman, 1807.01561) and from Miller–Peckner's lemma.

## Residual risk
- Not opened: Kim 2021 full text (abstract only; family unnamed; a generation statement), Liu–Masri–Young 2015 (zbMATH match on "least prime"; F = Q; no abstract online) and Y. Martin, JNT 2025 (binary forms; no abstract online).
- Forward citations of Yang 1801.05562 not checked: Semantic Scholar lists none, and the OpenAlex daily quota was used up on 22 Oct (429).
- Journals outside arXiv were searched only through metadata.
