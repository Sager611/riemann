RECORD: stronger counterexamples to Gale's conjecture in R⁴ and R⁵: β₅/D ≥ 0.9164673596 (126 points) and β₆/D ≥ 0.9283700398 (149 points), certified in exact arithmetic

# Gale's conjecture in R⁴ and R⁵ (2026-10-23)
Previous records: arXiv 2610.28828 (Voronov, 20 Oct 2026, Table 2) and its repository github.com/vsvor/Gale-conjecture (`best_known/`): R⁴, 61 points, β₅/D = 0.913148536990; R⁵, 87 points, β₆/D = 0.925763634680. Ball values (regular-simplex partition): σ₄ = √(5/6) = 0.9128709, σ₅ = cos(π/8) = 0.9238795. Details: [NOTE.md](NOTE.md).

## What the repository numbers are
- β_{d+1}(V)/D(V): the least possible largest part diameter over partitions of V into d+1 parts, divided by diam V; bracketed by an UNSAT witness graph and an explicit partition, both checked in integers.
- examples/d4_m61 (0.913149, same set as best_known/4d_61p) beats σ₄: a counterexample in R⁴. examples/d5_m32 (0.892305) is below σ₅, so it is not one; the R⁵ counterexample is best_known/5d_87p (0.925764).
- Posted only in the paper's TeX ("most likely far from optimal") and the repository; nothing else on arXiv (math.MG, math.CO, Oct 2026), Zenodo or the GXL evidence index (see NOVELTY.md).

## Result
| | R⁴ | R⁵ |
|---|---|---|
| set | [points4.txt](points4.txt), 126 points | [points5.txt](points5.txt), 149 points |
| β/D, exact | √(8399124363373578/10000000179265409) | √(2154677368122705/2500000047209683) |
| value | 0.916467359637… | 0.928370039894… |
| over the old record | +3.319·10⁻³ | +2.606·10⁻³ |
| over the ball value σ_d | +3.596·10⁻³ | +4.491·10⁻³ |

## Check
- [verify.py](verify.py): Python standard library only, exact integers, about 15 s for both files. It finds the critical squared distance dc², proves G(≥dc²) has no proper (d+1)-colouring by two independent exhaustive searches, checks an explicit (d+1)-colouring of G(>dc²) (so β = dc exactly), and compares exactly with σ_d and with the old records.
- DRAT proofs [cert4.drat](cert4.drat), [cert5.drat](cert5.drat) (CaDiCaL), "s VERIFIED" by drat-trim: `python3 verify.py points4.txt --cnf g4.cnf --sym; drat-trim g4.cnf cert4.drat`.

## Method and caveats
- Same method as our R³ record (../Gale/): core-SLP polish, then interior refinement (add convex combinations, re-optimise the critical core, prune to it, polish), plus a variant that pushes new points outward. About 15 chain-hours, single-threaded.
- The published R⁵ set was not locally optimal (polishing alone gives 0.92630). The searches had not converged. R⁶ (0.93696) was not attempted. Nothing in Lean.
