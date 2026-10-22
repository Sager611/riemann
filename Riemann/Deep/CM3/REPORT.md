DISCOVERY: least primes in every ideal class of almost all CM fields, for every A > 1 (was (3+√69)/10 ≈ 1.1307; 1 is sharp): for a fixed totally real F and every A > 1, all but O(X^{1−c}) CM fields K/F with prime (or δ₀-rough, odd) relative discriminant of norm ≤ X have degree-one primes in every ideal class in every interval (x, x + xD^{−κ}], x ≥ D^A (D = N𝔡_{K/F}); so every class has a prime of norm ≤ 2d_K^A. The exponent cannot go below 1: for every K, the principal class has no degree-one prime of norm below d_K/(4ⁿd_F²).

# CM3 report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Checks: [checks/](checks/). Brief: `.scratch/orca/cm3-brief.md`. Paper proofs, not refereed, no Lean. Unconditional.

## Results (complete proofs in the NOTE)
- Lemma P_F (task 1): in every order O_F + 𝔣O_K of a CM field K/F, the pair count over all elements (no coprimality to 𝔣) splits into a sum over squarefree 𝔢 | 𝔣 of class-group sums of the orders O_F + (𝔣/𝔢)O_K, each with its own unit weight ω (roots of unity times the relative unit index). So it is ≥ 0, and its dual is CM2's lattice count without the Möbius terms.
- Thm 1: CM2's family large sieve over pairs (K, χ ∈ Ĉl_K) is at the family size X^{3/2} for every length N ≤ X^{1−4δ} (CM2: + N³/X).
- Thm 2: Σ_K N_K(σ) ≪ X^{2−σ+ε} for ½ < σ ≤ 9/10, the density hypothesis of each K's class group family on average (CM2: (13−8σ)/(7−2σ)).
- Thm 3: CMQuartic's Theorem B for every A > 1 (CM2: A > 1.1307). A > 1 is used only for the zeros near Re s = ½, the GRH limit of the method.

## Where F ≠ Q changes something (task 3; NOTE §5)
- Removing 𝔢 is a lattice identity, not division by an element: 𝔢 may be non-principal and may share primes with 𝔣/𝔢; shifting the lattice 𝔢𝔠𝔯 by some λ ∈ F^× to an ideal prime to 𝔣/𝔢 finishes the class bookkeeping.
- Each block has its own unit weight; extra roots of unity (Q(ζ₅): ω = 10 for O_K, 2 for proper orders) and Hasse index 2 (Q(ζ₁₂): 24) only change positive weights.
- Totally positive generators, narrow classes, classes of F, the reference primes S₀, primes above 2, the dual count, the family 𝒦(X) and its counting: checked one by one, nothing breaks.
- Thm 3 still needs rough relative discriminants (genus characters of small conductor, CMQuartic §7.1); that restriction is inherited, not caused here.
- The transfer does not fail anywhere.

## Checks (pari/gp, hard alarms, all pass, < 20 MB)
- Lemma P_F exact in 141/141 (field, conductor) tests over Q(√5), Q(√3) (U⁺ ≠ U², Q(ζ₁₂)) and Q(√10) (class number 2, non-principal 𝔢): 372 blocks, 226 nonzero with 𝔢 ≠ 1, weights ω ∈ {2, 4, 6, 8, 10, 24}; all 12073 pairs with different 𝔣-parts vanish; ring class number formula 372/372.
- CM2's Poisson main term for norms with square factors (the pairs Lemma P_F keeps), F = Q(√5): brute/main 0.989–1.016; the smallest case 0.94, then 1.01 at larger X′.

## Novelty
- arXiv math.NT listings 2026-09 (812) and 2026-10 (625): 118 keyword hits, closest abstracts read (L. Yang 2610.22728, moments of toric L-functions over number fields; Gao 2610.27493; Diao 2610.28001); 8 arXiv searches; GXL API refused (429, daily quota). Nothing states or implies an almost-all exponent below 1.1307 for CM fields over any F ≠ Q.

## Not done
- All CM fields (not only rough relative discriminants), F varying, sectors, non-CM K/F: open as before (NOTE §8).

## Review (2026-10-22)
Reviewer: correct for the degree-one theorem (Lemma P_F with non-principal 𝔢, unit weights incl. Hasse index 2, Theorem 1 restricted to detector sequences, Theorems 2–3). The sharpness sentence needed "degree-one" (in Q(√−163), (2) is a principal prime of norm 4 < 163/4); fixed here and in NOTE §4.
