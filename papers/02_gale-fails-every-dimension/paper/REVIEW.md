# Review log: gale.pdf

Reader reviewers: fresh `reviewer` subagents, role "average STEM reader (BSc engineering/physics/CS), no discrete geometry", full text, captions and page images inline. Math reviewer: fresh `reviewer`, numbers and proofs.

## Re-verification (before round 1)

| Check | Result |
|---|---|
| `verify.py` (points7, points8) | β₈/D = √(8900856812318928/10000000213244525) = 0.94344351301…; β₉/D = √(9023994457279939/10000000157507010) = 0.94994706774…; both > σₙ; G(≥s) 1765 / 2195 edges, search nodes 1327; 23097 and 169807; 801496; 1.0 s and 36.8 s |
| `verify.py` X₃, X₄ (Gale/, Gale45/) | 0.8906210556…, 0.9164673596… |
| DRAT, drat-trim built from source | cert7.drat and cert8.drat: `s VERIFIED` |
| `family.py` | 996 rows n = 5…1000 pass; gap 7.9·10⁻⁴ (n = 5) … 1.43·10⁻³ (n = 7) … 2.9·10⁻⁷ (n = 1000); gap·n² = 0.290 at n = 1000 |
| `family_opt.py 3`, `4` | best 0.875161 < σ₃, 0.910647 < σ₄ |
| Integer σ test (2) | true for n = 3, 4, 7, 8 |
| X₇ vs V₇(θ₇) | every point within 0.1757 of a family point; X₈ norms 0.6166…0.6325 |
| Citations | Crossref DOIs (Gale, Borsuk, Kahn–Kalai, Hinrichs, Jenrich–Brouwer, DRAT-trim, Lean 4, mathlib); Voronov TeX (Gale's conjecture, σₙ, β₂,₃ = √3/2, Problem) |

## Rounds
