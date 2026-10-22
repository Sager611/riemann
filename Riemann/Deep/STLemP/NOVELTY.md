(a) NEW: no source states Σ_{|D|~X}Σ_{χ∈Ĝ(D)}|L(½+δ+it, f⊗θ_χ)|² ≪ X^{3/2+O(δ)} (or X^{9/4} for Sym²f) unconditionally; only a GRH, per-D form is known (Blomer–Brumley, JEMS 26 (2024), arXiv 2009.07093, Lemma 13, m = 1). (b) NEW: Lemma L as a sieve tool and Theorem 2₁ (average-over-D density for f⊗θ_χ) found nowhere. (c) NEW: no external statement of the thresholds 1, (23−3√5)/8, 3, or of any explicit power of |D|, for these Frobenius conditions in every class.

# STLemP novelty audit (adversarial, 2026-10-22)
Brief `.scratch/orca/stlnov-brief.md`; scratch `.scratch/stlnov/` (deleted after the run). Earlier checks: ../ChebST/NOVELTY.md, NOVELTY2.md, ../RSFam/NOVELTY.md, ../RCApps5/NOVELTY2.md (family sieve, m = 0).

## Searched
- arXiv math.NT monthly listings 2024-01 … 2026-10: 34 pages, 11,712 entries. Title keyword scan (class group, ideal/ring class, Heegner, CM point, toric, Waldspurger, theta, dihedral; Sato–Tate, Hecke eigenvalue, Fourier coefficient × primes/forms; large sieve, zero density, Lindelöf, moments × Rankin/quadratic/class/family): 393 titles checked, 14 abstracts read on arxiv.org/abs. API id_list: 897 abstracts (2024-01 … 03), then "Rate exceeded" (batch 7 of 30): stopped, no retry. Instead 87 arXiv full-search phrase queries (all years, abstract level; 0–14 hits each, 64 papers).
- zbMATH Open API (with reviews): 36 reference-text searches = forward citations of Blomer 2004, Michel–Venkatesh 2007, Templier 2008/2011, Liu–Masri–Young 2013/2015, Feigon–Whitehouse 2009, Sugiyama–Tsuzuki, Duke–Friedlander–Iwaniec 1995, Duke–Kowalski 2000, Michel 2004, Harcos–Michel 2006, Fouvry–Iwaniec 2003, Blomer–Harcos 2008/2010, Blomer–Brumley, Blomer–Brumley–Khayutin, Khayutin 2019, Ditchen, Nordentoft, Masri 2007, Michel–Ramakrishnan, Rodriguez-Villegas–Yang, Huang–Lester, Vatsal 2002, BTZ, Humphries–Thorner, Pascadi–Thorner, Thorner 2021, Lemke Oliver–Thorner, Kowalski–Michel 2002, Thorner–Zaman, Heath-Brown–Pierce, Wong 2019, Young; 59 review-text searches. 1,909 documents, keyword-filtered, about 30 reviews read.
- TeX read (rg): Templier 0902.2496, Nordentoft 2105.09130, Blomer–Brumley 2009.07093, Blomer–Brumley–Khayutin 2212.06280 (Van Order 1903.06722: source not served, abstract only).
- Semantic Scholar: citations of 17 arXiv seeds (146 works, 16 keyword hits, none relevant). OpenAlex: refused (IP daily budget exhausted, HTTP 429), not used. Zenodo: 15 phrase queries, nothing. openai/math: 749 preprints, CONTENTS.md (372 entries): only torus packets in degree ≥ 4 fields.

## Closest results (statements read)
| Work | Proves | Relation to (a)–(c) |
|---|---|---|
| Blomer–Brumley, JEMS 26 (2024), §8, Lemma 13 | under GRH, each quadratic E: Σ_{χ∈Cl_E^∨}L(½, π×χ)² ≪ exp(O(𝒞^{1+ε}))L(1,η_E)√D, π on PGL₂ | (a) for m = 1 under GRH, per D (stronger); nothing unconditional, no Sym² |
| Templier, Duke 157 (2011) (0902.2496, §1.2) | survey of known class-group moments for one D: 2nd moment of L(½,χ) (DFI, Templier 2008), 1st moment of L(½,f×χ) (MV), 1st moment of L′ (Ricotta–Vidick over D, Templier); mixed moments via Σλ_f(n²+d) | no 2nd moment of L(½, f×θ_χ), per D or over D |
| Ricotta–Vidick, CJM 2008; Ricotta–Templier, JTNB 2009 | average over D of Heegner heights = Σ_DΣ_χL′(½, f×θ_χ) | first moment only |
| MV 2007; Blomer 2004; Nordentoft, ANT 18 (2024) (wide moments; n = 2 is the 1st moment of L); L. Yang 2610.22728 (1st/2nd moments of L(½,χ), m = 0); Van Order 1903.06722 | one D; first moments, degree-2 moments, non-vanishing | not (a) |
| LMY 2013/2015, Feigon–Whitehouse, Michel–Ramakrishnan, Nelson 2013, Sugiyama–Tsuzuki, Holowinsky–Templier | averages over forms f, χ fixed | other family |
| Ditchen 1312.1502; Pascadi–Thorner 2508.14888; Humphries–Thorner; Duke–Kowalski 2000; Ganguly 2009 | large sieves over all D (N log³N + N^{1/2}Q^{5/2+ε}) or general GL_n families | at the AFE length N ≍ X: ≥ X³, not X^{3/2} |
| Kowalski–Michel 2002; BTZ, JEMS 2022; An 2004.14410; HT; PT; Cossaboom 2609.27715 | zero density for families or single forms | no average-over-D density for f⊗θ_χ ((b)) |
| Ditchen; Zaman 1710.08914; Gaudet (2408.01610); Sardari; H. Kim, JNT 2021; Hua 2025, Vaishya et al. | primes in classes without Frobenius condition; λ_f on one fixed form | not (c) |

## Notes for the NOTE / reviewer
- (a) baseline: unconditionally, at s = ½, per-D first moment (≪ D^{1/2+ε}, Waldspurger, L ≥ 0) × subconvexity (Michel 2004, Harcos–Michel 2006) already gives X^{2−δ₀} for the family sum, not X²; the gain to X^{3/2} stands.
- Derivable but not found stated (our sketch, not checked in detail): Böcherer's formula for Yoshida lifts (Dickson–Pitale–Saha–Schmidt, J. Math. Soc. Japan 72 (2020), Prop. 3.14, as used by Blomer–Brumley Cor. 4; f, g distinct, weights 2 and 2k−2, N squarefree with an odd number of prime factors, all p | N inert in Q(√D)) plus the classical Rankin–Selberg series Σ_T|a(F,T)|²det(T)^{−s} (pole at s = k) give Σ_{D~X}Σ_χL(½,f×θ_χ)L(½,g×θ_χ) ≪ X^{3/2+ε} for these D. It does not imply (a) (it needs f ≠ g).
- (b) The algebra of Lemma L (inflation along G_{p₀²D} ↠ G_D, λ_{χ∘π}(n) = λ_χ(n) for p₀ ∤ n) is standard (ring class characters); the positivity over all orders and the sieve beyond the conductor were not found.
- (c) As in the ChebST/RSFam audits: the thresholds improve this project's own 9/4, 3, 10/3; A = 1 is the GRH value (GRH: every D, routine, not found in print); 𝒜(3) = 3 rests on the preprint Pascadi–Thorner.

## Residual risk
Bulk abstracts only for 2024-01 … 03 (rate limit); later months via titles and 87 full-text queries. OpenAlex unavailable; Semantic Scholar citation lists sparse; zbMATH withholds some reviews (LMY 2015, Ricotta–Vidick, Dittmer–Proulx–Seybert 2015, Khayutin 2020; the last two are first-moment or m = 0 non-vanishing results by their citations/abstracts). Journals outside arXiv and zbMATH were not searched.
