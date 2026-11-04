DISCOVERY: first counterexample to Gale's conjecture in R⁷: β₈/D = 0.9434435130 > σ₇ = 0.9419651

# Gale's conjecture in R⁷, R⁸ and every dimension n ≥ 5 (2026-11-04)
Gale's conjecture says no bounded set in Rⁿ is harder to cut into n+1 pieces of small diameter than the ball; σₙ is the ball's regular-simplex partition value. Before this: counterexamples only for n = 3–6 (Voronov, arXiv 2610.28828, computer search), n = 9 (OpenAI's Borsuk counterexample) and n ≥ 64 (Borsuk counterexamples). Details: [NOTE.md](NOTE.md).

## Results
| space | set | β_{n+1}/D, exact | σₙ |
|---|---|---|---|
| R⁷ | [points7.txt](points7.txt), 84 points | √(8900856812318928/10000000213244525) = 0.9434435130 | 0.9419651451 |
| R⁷ | [points7_sym69.txt](points7_sym69.txt), 69 points (the family member) | √(4450002248247449/5000000078753505) = 0.9433983441 | 0.9419651451 |
| R⁸ | [points8.txt](points8.txt), 85 points (the family member) | √(9023994457279939/10000000157507010) = 0.9499470677 | 0.9486832981 |
| Rⁿ, every n = 5…1000 | family Vₙ, [family_params.txt](family_params.txt) | excess over σₙ from 7.9·10⁻⁴ (n = 5), 1.43·10⁻³ (n = 7) down to 2.9·10⁻⁷ (n = 1000) | |

## What it is
- One explicit family: two regular simplices in orthogonal subspaces R^⌊n/2⌋ ⊕ R^⌈n/2⌉; their vertices (poles), points near each pole pushed away from the other simplex's vertices, and points between pole pairs; 8 parameters, (⌊n/2⌋+1) + (⌈n/2⌉+1) + 3(⌊n/2⌋+1)(⌈n/2⌉+1) points.
- A short hand proof (NOTE §3) shows that the graph of its long distances has no proper (n+1)-colouring; the distances are 36 rational numbers per n. No SAT is needed, so this answers the paper's open problem (a counterexample for every n, described geometrically, proved without enumerating colourings) for every n ≥ 5. With Borsuk counterexamples for n ≥ 64 and our DRAT-certified sets for n = 3, 4 (../Gale, ../Gale45), Gale's conjecture fails in every dimension n ≥ 3.

## Checks
- [verify.py](verify.py) (standard library, exact integers, about 1 min): for points7.txt and points8.txt finds dc², proves G(≥dc²) has no proper (n+1)-colouring by two independent exhaustive searches, checks an explicit colouring of G(>dc²) (so β = dc exactly), compares exactly with σₙ.
- DRAT proofs [cert7.drat](cert7.drat) (188 KB) and [cert8.drat](cert8.drat) (1.2 MB), "s VERIFIED" by drat-trim: python3 verify.py points7.txt --cnf g7.cnf --sym; drat-trim g7.cnf cert7.drat.
- [family.py](family.py) (standard library, exact): all 996 dimensions 5…1000 in 1.6 s; --lemma re-checks the colouring lemma by exhaustive search for small cases (CaDiCaL also confirmed it up to n = 11).

## Method and caveats
- Found by making the R⁶ simplex-product structure exactly symmetric in R³ ⊕ R⁴ and optimising its 8 parameters with SAT and LP steps (about a minute), then reading off which distances are long and proving that pattern can't be coloured. A free search in R⁷ (twins, splits) added only 4.5·10⁻⁵.
- The excess is small (≈ 0.29/n² for large n); no asymptotic proof (not needed: exact check up to n = 1000, Borsuk beyond 64). The family fails for n = 3, 4. Nothing in Lean.


## Correction after the novelty check (orchestrator, 4 Nov 2026)
Borsuk counterexamples (hence Gale counterexamples) exist for every n ≥ 9 (OpenAI's Borsuk preprint, Cor. 7.1; finite ones from dimension 63, Ji 2608.12561), not only n ≥ 64; so the family is new for n ≥ 9 only as an explicit elementary alternative, and the new dimensions are 7 and 8. The two-simplex template appears as a search seed in Voronov's repository (22 Oct 2026). See NOVELTY.md.
