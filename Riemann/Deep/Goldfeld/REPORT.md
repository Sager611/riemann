NO DISCOVERY: I found no false or unjustified step in OpenAI's "Goldfeld's analytic density conjecture and the 2-converse for elliptic curves" (Oct 7 2026). Every step I checked holds. The Heegner/odd and irreducible parts were not checked, and nobody has machine-checked the paper.

# Audit: OpenAI, Goldfeld density + 2-converse (openai/math preprints, Oct 7 2026)
Source: `preprints/Goldfelds-analytic-density-conjecture-and-the-2-converse-for-elliptic-curves-October-7-2026/build` (TeX, 460 KB) in `.scratch/goldfeld/src-*`.
The Oct 7 version is identical to the Sep 23 version except for the date and three citation URLs (`diff -r`).
**Not machine-checked.** `lean/ComparatorChallenges` has 833 files and none matches goldfeld, converse, selmer, elliptic, heegner or twist. `lean/formalization.yaml` and `lean/docs` (242 files) have no entry. No published critique found (web search; arXiv 2610.17090 only lists the paper).

## How the argument fits together
- Density (Thm 1.2) = pointwise 2-converse (Thm 1.1) + Smith, arXiv:2503.17619, Thm 1.1. I checked Smith's statement verbatim: it counts d ∈ Z≠0, |d| ≤ H, and gives 2^∞-Selmer corank 0 and 1 each with density 1/2. That is exactly what §14 cites. Lemma 14.1 (all integers → squarefree, Möbius inversion) is correct. The only proportions claimed are 1/2, 1/2 and 0, and they follow exactly.
- Converse: put the curve at vertex 0 of a binary cube of twists h_x = h0·∏ q*^{λ_q(x)}. Arrange nonvanishing at all 2^b−1 other vertices with a normalized 2-adic valuation bounded by a constant B that does not depend on b. One determinant coordinate u ∈ O[G] (Kato class at rank 0, Heegner class at rank 1) then forces u_0(0) ≠ 0 (Lemmas 11.4–11.5, a Chevalley–Warning argument).

## Brief's targets
- **GRH, mollifiers, zero-free half-planes: none used.** The analytic input is only Prop 3.8. In (i), Hoffstein–Luo 1997 is quoted correctly (I checked the PDF): infinitely many odd fundamental d with L ≠ 0, ω(d) ≤ 4 and χ_d(p) = 1 on any finite set S. In (ii), Lemma 10.1 is a first moment of L′(1/2) in a fixed progression. I re-derived it: approximate functional equation, Poisson summation, zero frequency through sym² continuation, off-diagonal ≪ J^B X^{7/8} from convexity at Re s = 3/4. It is correct. The question whether the 2026 half-planes supply GRH-type input does not arise.
- **Uniformity over the twist family** (the steps most likely to fail). Checked:
  - Lemma 11.2: the Schur pole bound −Mh holds for any number of blocks.
  - Lemmas 11.4–11.5: proofs re-done. I also brute-forced Lemma 11.4 (`checks.py`): no failures above the threshold b > q(2^{M−1}−1); failures do occur below it.
  - Lemma 12.4 (residual concentration, reducible E[2]): the S-unit, Brauer and class-group counts over Q(ζ_{2^{n+2}})^+ give H¹ of size 2^n+O(1) with H² bounded, so C⁺ ⊗ O[G] ≅ O[G]²[−1] for every S.
  - Prop 12.3: the Poitou–Tate length bookkeeping is right. Euler factors at unused primes cancel exactly, active primes give t_q, and the result is v(u_x(0)) = D(h_x) + O(1).
  - Kato's Thm 12.5(3) in Kato's own text has no p ≠ 2 or image hypothesis (only 12.5(4) needs p ≠ 2), so its use at p = 2 over Λ[1/2] is legitimate.
  - The correction units U_q (eq. 12.9) and η (Lemma 13.21) involve no division by 2 per prime. In Lemma 13.22 the residual determinant orders are bounded: with odd Ψ-exponents they are at most 4.
  - Prop 13.3 (moving-conductor character): the Kummer compatibility argument holds.
- **Combinatorics.**
  - Thm 6.6 (address lemma): the degree-budget allocation in the two-symbol case is correct, including the ordinary-ring product of the top-degree parts.
  - Lemma 6.5: every F₂ matrix claim (cyclic and path programs) verified for b = 2..7.
  - Prop 5.5(iii) (contraction rule), the twin calculation in Prop 5.6, and the shell factor in Lemma 5.1 all re-derived and correct.
  - Lemmas 7.1–7.3 and Prop 7.4 (seeds, reserves, packets, even base) are consistent with the above.

## Why numerics cannot refute it
Every intermediate quantitative claim is a consequence of 2-primary BSD for twists: under BSD, D(h) is bounded for every rank-0 twist and the index bounds in Thms 3.5–3.7 and Prop 3.2 follow. Small-case computation therefore only tests BSD. A gap would have to be a logical one, not a false estimate. PARI/GP 2.19 was installed natively but turned out not to be needed.

## Not verified (highest risk first)
1. §13 Heegner part (109 KB):
   - ultraproduct "outer limit" behind the uniform ring bound (Lemma 13.20, Prop 13.16, Lemma 13.18);
   - height-one rank reduction (Lemma 13.14) and the local switch (Lemmas 13.11–13.13);
   - limits of complexes built from moving-conductor characters.
2. Odd coefficient system: Prop 4.2 (supersingular reduction / ternary theta series at precision 2^M) and the minima stabilized over ultrafilter stages (§4.3, Lemma 5.1 odd source).
3. Irreducible case: §§8–9 (Selmer matrices, saturation, interface repair) and Prop 5.8.
4. Reducible odd case: Lemma 7.9 "Nonempty K-split witness" and Prop 7.11 "Families with a fixed companion".

Verdict: the rank-0 converse for curves with rational two-torsion (cyclotomic route) looks sound in everything I checked. The rank-1 and irreducible routes depend on items 1–4. The paper is unrefereed and not machine-checked.
