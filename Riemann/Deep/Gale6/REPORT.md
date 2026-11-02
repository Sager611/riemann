RECORD: β₇/D ≥ 0.9373384468 in R⁶ (was 0.93696724796…), certified in exact arithmetic

# Gale's conjecture in R⁶ (2026-11-02)
Previous record: the 64-point set best_known/6d_64p_093696 of arXiv 2610.28828 (Voronov, Table 2) and github.com/vsvor/Gale-conjecture; exact value √(8779076237580130897884591797000439/10000000000000002601437547173996163) = 0.936967247964. Ball value σ₆ = √(7/8) = 0.9354143. Details: [NOTE.md](NOTE.md).

## Result
| | |
|---|---|
| set | [points6.txt](points6.txt), 99 points, coordinates in 10⁻⁸ Z⁶ |
| β₇/D, exact | √(8786033541902487/9999999889285618) |
| value | 0.937338446836… |
| over the old record | +3.712·10⁻⁴ |
| over σ₆ | +1.924·10⁻³ |

## Check
- [verify.py](verify.py): standard library only, exact integers, about 2 s. Finds dc², proves G(≥dc²) has no proper 7-colouring by two independent exhaustive searches, checks a 7-colouring of G(>dc²) (so β₇ = dc exactly), compares exactly with σ₆ and the old record.
- DRAT proof [cert6.drat](cert6.drat) (CaDiCaL), "s VERIFIED" by drat-trim: `python3 verify.py points6.txt --cnf g6.cnf --sym; drat-trim g6.cnf cert6.drat`.

## Method
- The old set is a perturbed product of two regular tetrahedra in orthogonal 3-spaces; it was not locally optimal (polishing alone gives 0.93721).
- Interior refinement (our R³–R⁵ method) stalled at 0.93722. What worked: SAT-screened vertex splitting (a core vertex becomes two twins that share its core edges minus one each; kept if still not 7-colourable), up to 6 splits per round, each round followed by core-SLP polish and pruning. About 3 CPU-hours.
- Side result: a fully symmetric 56-point counterexample ([points56.txt](points56.txt), 0.9368472 > σ₆, three parameters), below the record but with fewer points.
- Not converged (last rounds still gained ~10⁻⁵ each). Nothing in Lean.

