NEW RECORD: K₄(10,3) ≤ 768 (was 832: Bertolo, di Nasso, Di Pasquale, Jurcovich, Santisi 2004, in Kéri's tables; lower bound 313).
Independent check PASS: own C program, two methods (BFS in the Hamming graph H(10,4); brute-force minimum distance) agree word by word. 768 lines, 768 distinct words over {0,1,2,3} (read as Z₄; Hamming distance does not care whether Z₄ or F₄). All 4¹⁰ = 1,048,576 words lie within distance 3; covering radius exactly 3; counts at distance 0..3: 768 / 21504 / 253184 / 773120. Positions 1–8 project onto a Z₄-linear code of 256 words with minimum distance 4 (the octacode).

## Table values (checked 2026-10-23)
- Kéri tables (old.sztaki.hu/~keri/codes; www.sztaki.hu/~keri is 404): 4-5_tables.pdf (2009-10-15) gives 313–832, upper key x. The update log (to 2011-11-21) names K₄(10,3) only on 2005-01-19 (the 832).
- arXiv 2610.03760 (Maharaj): one version (v1, 27 Sep 2026), TeX re-fetched. Its two tables (density for R=2; the 48-word K₄(6,2) code) have no K₄(10,3) entry, and the text derives none. No repository ("available from the author"); 0 citing works (OpenAlex; Semantic Scholar rate-limited).
- Newer tables: Gijswijt–Polak 2504.01932 (rev. 19 Jun 2026) has no K₄(10,3) row. Florath's post-Kéri table (github florath/covering-codes-lean, 2026-09-21): 313–832. Mapika/coldcase bounds.json: 832. thiagopatzdorf/Matematica ledger (2026-10-15, Zenodo 23085769): no entry. Marosi 2608.19872 v3: q ≥ 5 only.

## Searches (all negative for anything below 832)
- GXL evidence API, 6 claims (≤ 768; best = 832; octacode radius 2; new K₄ bounds by search; q-ary amalgamated direct sum; octacode-based sum): no supporting passage.
- arXiv listings 2026-09 and 2026-10 of math.CO, cs.IT, cs.DM (3,941 titles parsed): only 2610.03760 is about q-ary covering codes. The arXiv API was rate-limited (429), so arxiv.org/search was used instead: "covering codes" (108 hits), "covering radius" (254), football pool, octacode, author Maharaj: nothing newer.
- OpenAlex (5 queries, 2009 on), Crossref, GitHub (gh code and repo search), Zenodo API (8 queries): nothing on K₄(10,3). Older sources (Cohen–Honkala–Litsyn–Lobstein 1997, Östergård 1991/94/99, Kéri–Östergård 2003/05) predate the 832 and are folded into Kéri's table.

## Is it routine? (adversarial)
- 768 = 256·12/4 is an amalgamated direct sum (ADS) of the octacode (K₄(8,2) ≤ 256, from 2610.03760) with a 12-word code B of length 3 and covering radius 1. The textbook q-ary ADS (Lobstein–van Wee) needs both codes normal. B is normal (norm 7 ≤ 7); the octacode is not (norm 14 > 11 at all 8 coordinates). So 768 does not follow from the textbook theorem: B had to be found by search.
- The other obvious route, 16·48: the ADS of the published 48-word K₄(6,2) code with the [5,3]₄ Hamming code. All 720 ways to join them (6×5 coordinate pairs × 24 symbol matchings) have covering radius 4, not 3.
- Direct sums of published codes give at least 1024. F₄-linear codes cannot go below 1024: the sphere bound (V₄(10,3) = 3676) forces codimension ≤ 5.
- Verdict: a new table record, but a modest one (one cell, −7.7%) that rests on a September 2026 preprint.

## Side findings
- thiagopatzdorf/Matematica claims K₄(10,4) ≤ 192 (was 208), K₄(6,3) ≥ 12 and K₄(7,4) = 10. The repo ruturajr-raval/quaternary-covering-code-8-4-upper-bound claims K₄(8,4) ≤ 27. Any K₄(10,4) attempt of ours must beat 192, not 208.
