NO REFUTATION: arXiv 2610.26686 checked line by line (Lemmas 3, 4, 6, 7, Thm 2, Prop 8, Thm 1; l. 321–683) and by exact computation for k = 4; no false step.

# Sunflower conjecture (all k) audit: report (2026-10-22)
Details: [NOTE.md](NOTE.md). Scripts: qc_check.py, tester_n2.py, dimcount.py, transfer_check.py (each < 90 s CPU, < 90 MB).

## How it beats the k ≥ 4 barrier
- The barrier (the naive count; Karingula–Lovett 2609.27044) is about testers on the Boolean cube. The paper never uses one. It bounds sunflower-free sets of vectors in [D]^n with D = 16k², then averages over random block systems (Prop 8).
- Key step: Lemma 6 (l. 409–532) picks the rainbow values c = (−1)^k ∏_{1≤i<h≤k−1}(a_h−a_i)/(a_0−a_i). Then every monomial of q_c has degree ≤ k−2 in one of x_1…x_{k−1}.
- Lemma 3 (l. 321–354) needs slices only in x_1…x_{k−1}, with one common space. Hence dim P = D^n·Pr[Bin(n,(k−1)/D) ≥ n/(k−1)], a power saving in D once D > (k−1)².

## Brief items
- (1) Diagonality, the lemma, the alphabet:
  - The tester only ever sees V^k for one block system. Repeated vectors force all vectors equal (Lemma 4); mixed blocks cannot occur.
  - Lemma 3 is a correct, self-contained slice-rank lemma with a common space. A reference variable is harmless.
  - The alphabet is [D] ⊂ Q, not {0,1}^b. Lagrange interpolation needs only D ≥ k−1.
- (2) k = 4, exact over Q:
  - n = 1: the q_c values, the decomposition identity and the degree property hold for D = 4…9 (also k = 3, 5). The naive choice c ≡ 1 fails (min degree D−1).
  - n = 2, D = 4, 5: T_c has the tester property on all 65,536 / 390,625 tuples, and all 52,900 / 254,016 monomials have the required form.
  - Contraction test: Σ T_c(a,·)h⊗h⊗h = 0 for every a and every h ⊥ P on [D]². With the naive tester it is never 0.
  - Dimension count = rank (D = 4–6, n ≤ 3). Sunflower-free V: |V| = rank(P|_V) ≤ dim P (9 ≤ 15, 27 ≤ 63); n = 1 is tight (3 = 3).
- (3) Transfer and constants:
  - Set sunflower ⟺ vector sunflower on all k-subsets of G: 5 configurations, 0 mismatches.
  - Pr[S ∈ G] = D^t/C(N,s), by full permutation enumeration. The ρ^t bound and the sum identity hold.
  - Constants checked for k = 3…400 (exact rationals); Chernoff step as an exact integer inequality for n ≤ 150. k = 4: 1 + ρ^{1/255} = 1.99819 ≤ 2 − 1/1024.

## Flaws, risk, consequence
- No fatal flaw. Cosmetic: l. 648 credits the binomial theorem with an inequality. The binomial theorem gives the equality; the inequality is ρ ≤ ρ^{1/(D−1)}.
- Residual risk: general n rests on the tensor-product step (checked by hand; computed for n ≤ 2 decomposition, n ≤ 3 rank). The claimed Lean proof (l. 744) is not public. Only v1 was audited.
- Unless a gap turns up elsewhere, this settles the 1970s Erdős–Szemerédi conjecture for all k ≥ 3. It is not our discovery. R138 can be closed as verified; nothing goes into DISCOVERIES or REFUTATIONS.
