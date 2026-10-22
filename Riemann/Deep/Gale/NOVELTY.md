NEW RECORD: β₄/D ≥ 0.8906210556 (our 170 points) is above every value I could find. The best one already posted is 0.8892811624: the full-precision 101-point set in github.com/vsvor/Gale-conjecture, the code repository of arXiv 2610.28828.

# Gale's conjecture in R³: adversarial novelty check (2026-10-22, about 21:15–22:00 UTC)

## Values above 0.8892809651
- **Repo vsvor/Gale-conjecture** ([gale_repo] in the paper; created 2026-10-19, content pushed 2026-10-22 12:56 UTC; no forks, issues or PRs). best_known/3d_101p_088928/points.csv has 17-digit coordinates and the SHA-256 quoted in the paper's TeX. Our verify.py gives exactly √(3954104928943379114534178573401953/5000000000000001234799393185523051) = 0.8892811624 (README: 0.889281162393).
- The paper's 7-digit appendix table gives 0.8892809651 (rechecked). So the true previous record is 1.97·10⁻⁷ higher than REPORT/NOTE say; our gain is +1.3399·10⁻³ (still +1.340·10⁻³). Exact rational check: ours > repo set.
- Other 3-D sets in the repo, all rechecked: 89 points 0.8888923465, 81 points 0.8887119259, 65 points 0.8881115021, fixtures 0.8878496. The 0.892–0.898 values in its docs are R⁴/R⁵ sets or cover upper bounds; paper_v5.pdf (AI-written algorithms paper) gives no 3-D values.
- **Nothing anywhere is above 0.8906210556.**

## The paper itself
- arXiv 2610.28828 has only v1 (20 Oct 18:47 UTC, public 22 Oct 00:00 UTC). Source: main_en.tex, main1.bib, four PNG views (no numbers), 00README.json; no ancillary files.
- The TeX says 0.88928 everywhere (abstract, theorem, Table 1, appendix). The comment "SELECTED_3d.zip / 3d_088891" is a rough directory label for the repo CSV above. Prior bounds as given: 0.88807 ≤ β₃,₄ ≤ 0.9654 (upper bound from vsvor/Borsuk-problem, which has only covers and upper bounds).

## Searched; nothing better found
- **arXiv API** (calls ≥ 4 s apart; result counts): Gale ∧ (math.MG ∨ math.CO) 134; "Gale's/Gale conjecture" 2; Borsuk 668 (newest 200 read); "Borsuk number" 11; partition ∧ diameter 273 (newest 200); "smaller diameter" 91. Only 2610.28828 is relevant.
- **Listings** math.MG and math.CO, 2026-09/10 plus new/recent: 2818 ids, all titles scanned, 2000 abstracts scanned. The API returned "Rate exceeded" at batch 11, so I stopped; ids 2610.08605–2610.30123 were checked by title only.
- **Citing works**: none (public < 24 h). GitHub code hits for 2610.28828 are RSS dumps, a paper list and one prize-tracking note quoting 0.88928. Kalai's blog post of 2026-10-21 links the paper in a comment, no new values. Semantic Scholar gave 429 and OpenAlex's daily budget was used up.
- **GitHub code/repo search** ("Gale('s) conjecture", 2610.28828, 0.88928, coordinate strings, the SHA-256, beta_{3,4}, Gale, Borsuk, optimization constants): nothing on β₄ in R³. Checked teorth/optimizationproblems (only C28), Tao-constant certificate repos, AlphaEvolve problems, acprk/Kodaxadev/igorbat Borsuk repos.
- **Other sources**: Zenodo (7 queries) has only Borsuk-number preprints (CP², 1.37^√d). openai/math has 749 preprints; only the 9-D Borsuk one is related, and it says nothing about Gale or R³. MathOverflow/MSE and web search show only the paper.

## Caveats
Tonight's arXiv mailing (23 Oct 00:00 UTC) is not visible yet. GitHub code search had not indexed the new vsvor repo, so it could miss copies of its files.
