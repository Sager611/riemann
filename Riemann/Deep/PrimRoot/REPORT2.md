NO DISCOVERY: Thm 3.1 (one-sided marked Type II estimate) of OpenAI "The Poisson–Dirichlet law for prime predecessors" (24 Sep 2026) re-derived step by step; no false step, missing range or circular use found.

Source: openai/math preprints/The-Poisson-Dirichlet-Law-for-Prime-Predecessors-September-24-2026/paper.pdf (blob c46ae43, unchanged). Text `.scratch/pd/pd.txt`; checks `.scratch/pd/checks.py`. No Lean: no ComparatorChallenges file for it, CONTENTS entry 011 has no Lean link.

## The proof, step by step (with each saving)
1. Expand the marks: W(mn−1) = sum over b0 (one prime per band) dividing mn−1; put mn−1 = b0·h. F ignores band primes, so F(mn−1) = F(h). Prime-square terms cost X·e^{−L^0.1}.
2. Cauchy–Schwarz over h ≪ X/Y: need Q_Y ≪ XY·L^{−D1}. Pairs a,b sharing a prime: Y²e^{−L^0.1}. Coprime pairs ⇔ bmn − ars = b − a (checked).
3. Split 1_{t=b−a} into major arcs (denominators ≤ L^{A0}, radius 2L^{A0}/Y) and minor arcs.
4. Minor arcs: pads D (L^0.01 primes per band), positions P = (Dbm, s), operator A = GSTSG. Hölder turns the trace bound (1/UV)Σ⟨u_P,(AA*)^R u_P⟩ ≤ L^{−E0·N}, N ≈ L^0.5, into pairing ≪ UV·L^{−E0+O(1)}; the box loss H^{1/N} = O(1) since log H ≪ L^0.21.
5. Trace → model (Lemma 3.4): roots equidistribute mod S ≤ e^{L^0.98} (Kloosterman 4th moment), so each band prime gets an independent uniform line; error x^{−c}, far below L^{−AN}. Never-active primes: Bonferroni at degree L^0.76.
6. Memory identity (4.5)–(4.14): each prime's line average equals its lifespan sum (one 1/(p+1) per lifespan, V^{−1} per active run, q_j per pending hit, ghost ends −η'). Brute-force check: exact (error 4e−16).
7. Edge norms ≤ L^{−G}: clean edge = minor-arc large sieve on a lattice box A'B' ≪ Y (good test (i)) ≪ 1/d + Y^{−0.1}log Y + d log d/Y ≪ L^{0.2−A0}. Dirty raw: phase separation (good test (ii)) leaves one omission tuple, row B·M^{−K}, column C_K → L^{1−0.005K}; reversal table checked both ways (anchor line, ban, ρ-trick for repeated hits). Dirty comparison: one line congruence → e^{−cL^0.1}. Ghosts ≤ C_K since √q < ρ².
8. Birth distinctness: partition Möbius sum (verified); ranks ≥ N·log L/L^0.1 die by atoms e^{−L^0.1}; lower ranks touch o(N) edges, cost L^{o(N)}.
9. Goodness removal (Lemma 3.2: e^{−cL^0.1}), pad stripping (M²e^{−L^0.1}), O(L) pad dyads.
10. Major arcs (Prop 5.1): Mellin separation; Lemma 5.2 (primes minus rough model are unbiased against characters mod ≤ L^{A0} at all heights ≤ 2XL^B); Lemma 5.3 (one band prime is large on only e^{O(L^0.91)} unit intervals).

## Inputs and weight
- No level-of-distribution input: the minor-arc bound holds for arbitrary bounded coefficients. Prime input = Siegel–Walfisz (ineffective) below height L^{B2} + zero-free strip 1−β ≥ c(log log x)²/log x at heights 1..x³ (Lemma A.8, from a re-proved Vinogradov mean-value bound). No Siegel-zero-free hypothesis.
- Parameter order sound: E0 → A0 → K (0.005K > G+2); every L-exponent in the Hölder and pad steps is independent of K.
- Counterexample attempts failed: congruence patterns mod a minor-arc prime, characters of conductor > L^{A0}, a Siegel zero at modulus e^{L^0.2}. The mark b0, invisible to F, averages each away.
- Use in §7.2: F = clipped f_x − C_x is band-invariant (slot primes > x^ε); both scales ≥ x^γ, so no range is missing.

## Numerical checks (`.scratch/pd/checks.py`, all pass)
- Lemma 3.3: T_m ≤ 2m²τ3(m) for m ≤ 36 (max ratio 0.45); root counts mod 7 within 11% of UV·vol/(ζ(2)|SL2|) at U = V = 300.
- (4.32): ‖T(θ)‖² ≤ (A'/d+1)(B'+d log 2d); determinant-matrix norm equals ‖T‖² exactly.
- Toy clean edge (Y = 600): fiber norm 1.4–4.2 on arcs d ≤ 3, 0.22–0.26 on minor arcs.

## Residual risk (not line-checked)
- App. A Lemmas A.5–A.6 (log-phase sums, approximate formula (A.13)): structure only.
- §§6–8 (Thm 6.1 two-sided estimate, mark removal, PD law) and the block sieve of [S]: outside this brief.
