STILL NEW: no source found (classical papers, arXiv, Hexagon, Zenodo, OpenAI, ATLAS, GitHub) states or implies Σ_χ N(σ;t) ≪ |D|^{2σ(1−σ)+ε}T³ over the ring class characters of D. Every bound found is weaker in |D| at σ = 3/4 and at σ = 0.9.

# Second, adversarial novelty check of Theorem 1 (2026-10-21; raw data in `.scratch/novelty2/`)

Exponent of |D| at fixed T (smaller is better). Log-free bounds win only within O(1/log|D|) of σ = 1, where Theorem 1 loses |D|^ε. Pascadi–Thorner is weaker for every σ because 3/(3−2σ) > 2σ (4σ² − 6σ + 3 > 0); Kowalski–Michel is weaker because 4σ² − 2σ ≤ 2 < c₀.

| Source (statement read) | Bound for Σ_χ N(σ,T,χ) | σ = 3/4 | σ = 0.9 |
|---|---|---|---|
| Theorem 1 (ours) | \|D\|^{2σ(1−σ)+ε} | 0.375 | 0.18 |
| Pascadi–Thorner, arXiv 2508.14888, Thm 1.3 (n = 1 over K; constant depends only on [K:Q]) | (\|d_K\|^{−2}Q³\|S\|⁴T⁸)^{(1−σ)/(3−2σ)+o(1)}, Q ≤ \|D\|, \|S\| = h(D) | 0.50 (fundamental D) to 0.83 (fixed d_K) | 0.25 to 0.42 |
| Kowalski–Michel, Pacific J. Math. 207 (2002), Thm 2 (dihedral GL₂: A = 1, d > 1/2) | T^B\|D\|^{c₀(1−σ)/(2σ−1)}, c₀ > 5nA/2 + d > 5.5, σ ≥ 3/4 | > 2.75 | > 0.69 |
| Weiss 1983; Thorner–Zaman arXiv 1510.08086; Zaman thesis 2017, (5.2) | (e^{O(1)}d_K²QT⁴)^{81(1−σ)}, d_K²Q ≥ \|D\| (log-free) | ≥ 20 | ≥ 8.1 |
| Duke, Acta Arith. 52 (1989), Thm 2.1 (scan, read via OCR) | Σ_{\|m\|≤T} N(σ,T,χλ^m) ≪ T^{3n(1−σ)/(2−σ)}: fixed field and modulus, T-aspect only | no \|D\| aspect | – |

Searched in this round (the first check did ~105 arXiv queries):
- Hexagon (hexagonmath.org API): all 323 works, 2609.00001–2610.00419, titles and abstracts; full TeX of the 102 math.NT works searched with rg (ring class, class group character, binary quadratic form, zero density, large sieve). Only hit: 2610.00178 §54 (least prime of a form via Thorner–Zaman 2019, exponent unspecified, no density bound). Hu 2610.00216 and Cureton–Thorner 2610.00200 have no family density.
- arXiv API: id_list 2610.00001–2610.00400 (400 records, no number-theory hit) and the 400 newest math.NT papers (7–20 Oct). Closest: 2610.27493 (moments of angular Hecke L-functions), Frei–Sofos–Widmer 2610.18491v2 (ℓ-torsion), 2610.17675 (sextic large sieve). None has a class group density bound.
- Zenodo API (~70 queries with AND phrases; anonymous page size is capped at 25; a Mozilla user agent gets HTTP 403): nothing. Closest: Jaswal 23403019 (class numbers ≤ 3000 assuming 7/8), Baier 23386232 (anatomy of the quasi-RH method), Daugherty–Ryan–Ward 23462320.
- openai/math: 749 preprint names; GitHub code search inside the repo: "ring class" only in BSD/Selmer papers; "class group character" and "zero density" with "class group": 0 hits. atlas-lean (14,757 paths): RingClassField, QuadraticForms/IdealClass and Dirichlet large sieves, but no class group density. anthropics/formal-math (zeta23): nothing.
- GitHub-wide code search (22 queries, logged in): nothing states the bound. AI-agent repos (perfect-cuboid-game, math-street/math-st) only cite Thorner–Zaman and Sardari, and note that no class group large sieve is available.
- Classical, beyond the table: Ditchen arXiv 1312.1502 averages over D and proves no density bound for a single D. Gaudet arXiv 2408.01610 (2026 version too) avoids density bounds and has only a single-L-function lemma. Blomer 2004 has neither a large sieve nor a density bound. Blomer–Harcos (9 joint arXiv papers) prove subconvexity, not density. Heath-Brown (59 arXiv titles) has nothing on class group density; his Dedekind zeta density is for a fixed field. Conrey–Iwaniec (math/0111012) studies zero spacing, not density.
- GXL evidence API (3 claims): density claim 0 supports and 0 refutes out of 265 passages. The large sieve claim's single "support" is the power-moduli additive large sieve (irrelevant). Gaudet's text says it avoids class group density bounds.

Residual risk: Zenodo indexes metadata only; arXiv full text was searched for only 21 + 102 sources; the keyless web search returned noise. The verdict stays "new but routine" (classical method).
