# Gale's conjecture fails in R⁷ and R⁸, hence in every dimension n ≥ 3; an explicit family with a hand proof for n ≥ 5

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** (1) Explicit finite sets in R⁷ (84 points) and R⁸ (85 points) cannot be split into n + 1 parts of diameter at most σₙ·D: β₈/D = 0.9434435130 > σ₇ = 0.9419651451 and β₉/D = 0.9499470677 > σ₈ = 0.9486832981 (exact square roots of rationals). (2) An explicit 8-parameter family Vₙ (two regular simplices in orthogonal subspaces R^⌊n/2⌋ ⊕ R^⌈n/2⌉, their vertices, points near the poles and points between pole pairs) is a counterexample for every n = 5…1000, with a short hand proof that its long-distance graph has no proper (n + 1)-colouring and exact rational checks of β_{n+1}/D > σₙ. With known counterexamples in the other dimensions, Gale's conjecture (1953) fails in every dimension n ≥ 3.

**Before.** Counterexamples for n = 3–6 (Voronov, arXiv 2610.28828, computer search; our stronger ones: [R³](gale-r3.md), [R⁴–R⁶](gale-r4-r5.md)) and every n ≥ 9 (OpenAI's Borsuk counterexamples, openai/math, Cor. 7.1, public 6 Oct 2026; finite Borsuk counterexamples from dimension 63, Ji 2608.12561). R⁷ and R⁸ were the only open dimensions; Voronov asked for a geometric family proved without enumerating colourings.

**Method.** The R⁶ simplex-product structure made exactly symmetric in R³ ⊕ R⁴, its 8 parameters optimised with SAT and LP steps; then the pattern of long distances was read off and proved uncolourable by hand (NOTE §3). The two-simplex template itself appears as a search seed in Voronov's repository (22 Oct 2026, "not a certificate"); the parameters, the lemma and the sets in R⁷, R⁸ are new.

**Certificate.** [verify.py](../Riemann/Deep/Gale7/verify.py) (standard library, exact integers): finds the critical distance, proves the long-distance graph has no proper (n+1)-colouring by two independent exhaustive searches, checks an explicit colouring just above it, compares exactly with σₙ; DRAT proofs cert7.drat / cert8.drat ("s VERIFIED", drat-trim); [family.py](../Riemann/Deep/Gale7/family.py) checks all 996 dimensions 5…1000 in seconds and re-checks the lemma graphs. All re-run by the orchestrator on 4 Nov 2026: pass.

**Review.** The colouring lemma is sound and the combination with the known counterexamples is valid; significant.

**Novelty.** No Gale counterexample in R⁷ or R⁸ in any source checked (arXiv incl. all versions of 2610.28828 and the Oct–Nov 2026 lists, Zenodo, openai/math, GitHub, our own notes), 4 Nov 2026 ([novelty report](../Riemann/Deep/Gale7/NOVELTY.md)). Not new: existence for every n ≥ 9.

**Files.** [report](../Riemann/Deep/Gale7/REPORT.md), [note with the proof](../Riemann/Deep/Gale7/NOTE.md), [R⁷ points](../Riemann/Deep/Gale7/points7.txt), [R⁸ points](../Riemann/Deep/Gale7/points8.txt), [family parameters](../Riemann/Deep/Gale7/family_params.txt)
