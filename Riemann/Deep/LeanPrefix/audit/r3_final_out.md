VERDICT: FULLY-PROVED-AND-FAITHFUL

Remarks

- **Theorem A:** Quantifiers, range restrictions, parity, diagonal equality, gcd, divisibility and sign exception match the Claim. `G = ⊤` means generation of the full symmetric group on `Fin n`.

- **Line domain:** From `2 ≤ k`, `k < m` and `m + k = n - 1`, one obtains `n = m + k + 1`, hence `m < n` and `m + k < n`. Natural-number subtraction introduces no exceptional cases.

- **Orbits:** The subgroup action, 1-based congruences, orbit-quotient cardinality and transitivity equivalence state the claimed decomposition. Integer negation correctly represents the negative residue class. Natural division by two gives the required floor.

- **Even branch:** Transitivity, non-primitivity and a block satisfying `2 * B.ncard = n` give exactly two translated blocks. Since the domain forces `n ≥ 6`, this block is neither a singleton nor the whole set.

- **Wreath branch:** Omitting the explicit gcd hypothesis is harmless: with `α = k + 1`, the line gives `m + 1 = n + 1 - α`; therefore `α ∣ n` implies `gcd(m + 1, α) = 1`. Transitivity follows from the preceding conjunct, so non-primitivity genuinely means imprimitivity.

- **Wreath object and action:** `Wreath α q` is the permutational wreath product indexed by the vertices of the q-gon, not the regular wreath product. Its inverse-index action and target-index permutation formula agree. Transport by `e.symm` correctly places its permutation image inside `Perm (Fin n)`. Both permutation equivalence and abstract isomorphism are explicitly asserted. Here `α ≥ 3` and odd `q ≥ 3`, excluding degenerate, nonfaithful polygon actions.

- **Primitive branch:** Its hypotheses match “otherwise” within the odd, gcd-one case. Mathlib primitivity, containment of `alternatingGroup`, equality with it **iff** the stated congruence condition, and generation of `S_n` otherwise are all present. On the line, `α ∣ m ↔ α ∣ n`.

- **Literature inputs:** LZZW matches the quoted first sentence with its standing range assumptions. Jones supplies only weakened conclusions of the quoted classification. The full-support cycle condition supplies the required degree restriction. Neither input contains an additional result of this formalization.

- **Nonvacuity:** The even branch has `(10,6,3)`; the wreath branch `(9,6,2)`; the primitive alternating branch `(13,8,4)`; and the primitive symmetric branch `(7,4,2)`.

The verdict is conditional on the two explicitly permitted literature inputs and the supplied compilation/axiom evidence.