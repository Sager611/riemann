# NOVELTY (checked 2026-11-11)

Question: has anyone optimised P (or r) in Bui-Hall-Subirà Jorge's L⁴-Wirtinger large-gap argument, or beaten
the unconditional Λ > 3.18? Answer: no follow-up found; nothing beats 3.18 unconditionally.

| Source | Query | Result |
|---|---|---|
| arXiv API | id_list 2511.14415 | v1 only (2025-11-18), no journal ref |
| OpenAlex | works citing W7106101550 (= 2511.14415) | 0 |
| Semantic Scholar | citations of arXiv:2511.14415 | 0 |
| arXiv listings 2026-10, 2026-11 (math.NT, math.CA, math.CV), `.scratch/lst.py` | titles `gaps? between\|large gaps\|consecutive zeros` (also `Wirtinger\|amplif` on 6 Nov) | 2610.14127, 2610.15291, 2610.22927, 2610.26038, 2610.28527, 2610.31744, 2610.36293: none unconditional large gaps for ζ |
| arXiv API | abs "amplified" AND "fourth moment" | 2606.27323 (Durkan-Page, two-piece amplified moments): cites BHSJ (main.tex l.355) for the amplifier, uses it for moment lower bounds only (no "gap" hit); 2610.15291 unrelated |
| arXiv API | ("gaps between consecutive zeros" OR Wirtinger) AND zeta, 2024-2026 | 2610.26038 (RH, small gaps), 2604.05733 (small gaps) |
| arXiv API | au:Milinovich AND ti:gaps | 1410.3635 (3.18), 0910.2052 (RH, 2.69), 2208.02359, 2610.31744: no newer unconditional value |
| arXiv 2610.28527 (Durkan, Pearce-Crump, Oct 2026) | TeX l.45, l.73-81 | RH: λ > 3.63846723 for every primitive Dirichlet L incl. ζ (Q = 1); l.81: 3.18 is the unconditional Wirtinger value; no BHSJ citation |
| arXiv 2611.14486 (Lamzouri, Nov 2026) | abstract | small gaps (λ0 = 0.7892), no large-gap theorem |
| Zenodo | "gaps between zeros" (8 hits), "large gaps" AND zeta (4 hits) | none relevant (24095120 is conditional r-gaps) |
| HAL | "gaps between zeros" | 0 hits |
| repo | rg `2511.14415\|Subir\|large gaps between\|Wirtinger` | Harvest5 candidate 4; Riemann/Tree/R421.lean l.13 (mislabels 3.18 as RH-only, see REPORT) |

Verdict: no one has optimised P in this method, as far as these searches show, but the resulting 2.903 is below
the unconditional record 3.18: not a record, not a discovery entry.

Follow-up (optimising P in Bui-Milinovich's own L² argument): BM (QJM 2017, doi 10.1093/qmath/hax047, OpenAlex
W2609434580) has 22 citing works: r-gaps, small gaps, moments, explicit results (Simonič-Trudgian-Turnage-Butterbaugh
2021), random matrices; none improves 3.18. BHSJ 2511.14415 l.116 (Nov 2025) and 2610.28527 l.81 (Oct 2026) still
quote 3.18 as the best Wirtinger value. The float value Λ > 3.25 (REPORT, Follow-up) is not in the literature we found.
Record (follow-up 2): the certified Λ > 3.2648 uses Bui-Milinovich's own published criterion with a degree-14
amplifier. No source above certifies or claims any unconditional value above 3.18 (checked 2026-11-11).
