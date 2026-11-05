NO FLAW FOUND: reduction re-derived by hand; all provided checks reproduce; 5,000+ critical-stratum leaves re-excluded by my own Arb/LP/Krawczyk prover with masses unrestricted; full coverage of the n=6 critical tree verified. Not checkable: Lei–Zhang entropy theorem, n=7 N=6 tree (54 GB logs withheld), smaller-strata leaves (never logged).

# Audit: Wikström, computer-assisted Krzyż n = 6, 7 (Zenodo 23964169 v2.0; arXiv:2610.15206)
Claim: |a_n| ≤ 2/e for n = 6, 7, equality only for rotations of f_n = exp(−(1−z^n)/(1+z^n)).

## Reduction (checked step by step)
1. Extremal f = e^{−g}: Toeplitz data on the cone boundary ⇒ f = exp(−Σ_{j≤N} λ_j(1+α_j z)/(1−α_j z)), N ≤ n (HSZ 1977; Agler–McCarthy 2021, Lemmas 3.2–3.3). Quoted correctly.
2. N = n: AM Thm 5.1/Props 6.4, 6.6 + Lei–Zhang entropy inequality with equality ⇒ E(γ) = −log η, η ≤ 2/e, Vandermonde ⇒ λ_j = 1/n ⇒ rotation of f_n. Correct, given Lei–Zhang.
3. N ≤ n−2: no first-order conditions; proves sup|a_n| < 2/e on D_N (heaviest atom at 1, ordered angles, reflection, T_lo ≤ t ≤ T_hi). The KKT test A is gated to N = n−1 in code (eval.c:349) and closed 0 leaves in these strata.
4. N = n−1: h ≥ 0, h(θ_j)=0 ⇒ Fejér–Riesz ⇒ Σ_j(λ_j−1/k)e^{ikθ_j} = T_k(q)/k (k=1..n) after rotating X to q ≥ 1; on that locus |a_n| = e^{−t}C(θ,q) (I matched it to the direct |a_n| to 1e−15 at 50 KKT points). Root θ_1∈[0,π], θ_j∈[0,2π], q∈[1, 1.31 | 1.3], λ_j∈[0,T_hi]: the normalizations lose nothing and the bounds hold.
5. Constants (kz_constants.py, Arb): T_lo literal ≤ T_lo (gap 7e−15); threshold literal < 2/e (gap 9e−17); mass cutoffs n=6: 0.7354864, n=7: 0.7353315 < 2/e; T_6(1.31) = 50.3 > 46.94, T_7(1.3) = 99.7 > 63.05.
6. Kernel source read line by line (kernel/fixed.cpp, fixed.hpp, native.cpp): Taylor remainders, Neumann certificate, projected residual, mean-value C bound, Gauss–Seidel, LP rows and witness, directed rounding. No error found.

## Provided checks (venv, python-flint 0.9)
shasum all OK; verify_release PASS; smoke PASS; replay leaves 9,331/9,331; splits 3,326 certified + 618 inconclusive, then both covers 618/618; summary identical to the deposit.

## Own verification (scripts in this folder; outputs in results/)
- kz_check.py: separate leaf prover (tight trig ranges, LP relaxation with Arb-verified Farkas/dual certificates, value cut t ≤ ln C_max + 1 − ln 2, Krawczyk on (θ,q,λ)). Masses run over all of [0,T], so it uses none of the logged contractions. Soundness test (kz_soundness.py): boxes around real KKT points are never excluded once the threshold is lowered below their |a_n|.
- Leaves re-excluded: n=7 N=6: 411/411 (150 random, 30 deepest down to depth 65, 126 chart_value, 137 the author's replay had to subdivide). n=6 N=5: 325/325 exported, plus 3,740/3,740 sampled from the full log (3,000 uniform over 18.16M leaves, so the bad fraction is < 0.1% at 95% confidence). n=5 N=4: 156/156.
- Hardest leaves (kz_locate.py, kz_hard.py): all 19 KKT points of n=6, N=5 (best |a_6| = 0.700865261, margin 0.0349) lie in production leaves closed by chart_value; the runner-up's is at depth 51 of 52. No contraction on any of these paths dropped the true (θ,q,λ); a perturbed control point was flagged. My prover closes all 19.
- Coverage (kz_cover.py, own streaming audit of the complete n=6 N=5 log: 36,317,049 nodes, 18,158,525 leaves; and of the n=5 log): the root contains the normalized domain; every node is used exactly once; each parent lies inside its input box; parameter contractions follow from ordering, conjugation and q ≥ 1; bisections are exact; nothing is left pending; counts match the summary. n=7 N=6 (kz_chain.py): the 29 batch frontiers rebuilt from METADATA match the committed hashes; the chain starts at the canonical root; nodes = 2·leaves + |F_out| − |F_in| in every batch; the final frontier is empty.
- Smaller strata: the author's job-log audits re-run and pass (n=6 N≤5, n=7 N≤5). I built the C kernel from the deposited source and ran kz_vtest.py on 1,800 boxes, near-maximal ones included: 0 bound violations. My own branch-and-bound (kz_atoms.py) proves N=1 from scratch for n=6 and n=7. Floating-point multistart maxima: n=6 N=4 0.6630, n=7 N=5 0.6742 (gaps ≥ 0.06).

## Could not be checked
- Lei–Zhang entropy theorem (arXiv:2605.03949, v2): essential and not audited. The cited AM 2021 results were taken as given.
- n=7 N=6 node trees: the 54.46 GB logs are withheld ("available by arrangement"). Only the batch chain and samples could be checked.
- Smaller-strata leaves are never logged: coverage inside a job rests on the DFS code plus Kraft counts, and the arithmetic on the Taylor-model kernel (enclosures.tex not re-proved; n=7 N=5 alone is 6.7e10 nodes). My own branch-and-bound was too slow beyond N=1.
- The 18M n=6 leaves were only sampled. The 618 (n=7) and 188 (n=6) logged mass contractions that the author's reference could not certify are not directly certified by anyone; my leaf checks do not rely on them.

## Known elsewhere
n=6 was first claimed by the same author in arXiv:2610.15206v1 (10 Oct 2026). n=7 was first claimed by Xie, arXiv:2611.00397 (30 Oct 2026); this deposit gives a second, independent proof. No refutation found on arXiv or Zenodo (Stupin 2504.10223 covers only n ≤ 3; Krushkal's claims are not accepted).
