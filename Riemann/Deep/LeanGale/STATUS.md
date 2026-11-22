WORKING
UTC: 2026-11-22T12:56:40Z
sorry count: 0
files: Cert.lean (SWAR refutation certificates), Splx.lean (unit simplices, every m ≥ 1), Row.lean (family row check), Pts.lean (explicit point sets), P7.lean + Data/P7.lean (84 points, certificate), P8.lean + Data/P8.lean (85 points, H(5,5) labels), Data/Rows00–19.lean (family rows n = 5…1000), tools/{gen,rows,fpcap,probe,closure,hom}.py/sh
main theorem: (planned) LeanGale.gale_fails : ∀ n ≥ 3, ¬ GaleConjecture n (EuclideanSpace ℝ (Fin n), Metric.diam)
literature inputs (planned): Hinrichs 2002 Thm 1 (Discrete Math. 243, 253–256: Borsuk fails in all dimensions > 322) for n ≥ 1001; Voronov arXiv 2610.28828 (unrefereed) Thm 1 (n = 3) and Table 2 row n = 4, unless our own R³/R⁴ sets get kernel certificates
NOT-FORMALIZED: none so far
last reviewer verdict: none
next step: build Rows01–19, Family.lean (∀ n ∈ [5, 1000]), Statement.lean (EuclideanSpace, Metric.diam, β, GaleConjecture), literature inputs; reviewer rounds

## Progress
- Memory cap lifted (orchestrator, 22 Nov): compiles via .scratch/lean4/leanrun.sh; peaks so far ≤ 4.0 GB.
- P7 (points7.txt, 84 points): no H(4,5)-homomorphism into G(≥ dc²) (tools/hom.py), so the lower bound is a kernel-checked refutation certificate: 18887 tokens (1204 splits, 16433 unit propagations, 45 colour-symmetry cuts), state = one natural number with 8 bits per vertex (Cert.lean; run_ok by decide +kernel in ~30 s). Upper bound: explicit 8-colouring with same-colour squared distances ≤ dc².
- P8 (points8.txt, 85 points): the labels map H(5,5) into G(≥ dc²) (min over edges = dc² exactly), lower bound by the Gale7 colouring lemma; explicit 9-colouring for the upper bound.
- Family: every row of family_params.txt checked in integers (Row.rowOK: 100 descriptors, exact σ test); Rows00 (n = 5…54) kernel-checked in 17 s.
- Literature (scouts, verbatim): Voronov 2610.28828 l. 88–107 (Gale's conjecture "the ball is extremal", σₙ bound), Thm 1 (n = 3: 101 points, β₄ > 0.88928), Table 2 (n = 4: 61 points, β₅ = 0.91314…, no stated certificate); Hinrichs 2002 Thm 1 ("Hence Borsuk's conjecture fails in all dimensions exceeding 322"); Pikhurko arXiv math/0202112 Thm 1 (all n ≥ 321); Gale 1953 PDF not retrievable (AMS serves HTML only).
