VERDICT: FULLY-PROVED-AND-FAITHFUL

### Remarks

1. **Scope:** Verdict concerns the proposed amended Claim, relative to the two permitted literature hypotheses and the supplied compilation evidence—not the original, unqualified Claim.

2. **Theorem A:** All bounds, quantifiers, divisibility conditions, parity conditions and the biconditional match exactly. `G = ⊤` correctly means generation of `S_n`.

3. **Line hypotheses:** From `2 ≤ k`, `k < m` and `m + k = n - 1`, one obtains `n = m + k + 1`, hence `m < n` and `m + k < n`. Their omission from the line conjunct introduces no exceptional cases.

4. **Orbits:** The stated orbit-membership equivalence gives exactly the classes of 1-based positions modulo sign modulo `g`. Here `0 < g ≤ k+1 < n`, so every residue class occurs among the positions. These classes have the canonical representatives `0 ≤ c ≤ ⌊g/2⌋`. The quotient cardinality and transitivity biconditional therefore express NOTE B.1 faithfully.

5. **Wreath branch:** With `α = k+1`, the line gives `g = gcd(n+1, α)`. Thus `α ∣ n` already implies `g = 1`; the missing explicit coprimality hypothesis is redundant, not a mismatch. Moreover, `α ≥ 3`, and oddness plus `n ≥ 2α` implies `q = n/α ≥ 3`. Both polygon actions are therefore faithful; the small-degree dihedral degeneracies are excluded.

6. **Wreath construction:** `Wreath α q` is the permutational—not regular—wreath product: its base is indexed by the `q` vertices, with inverse-pullback action. `wreathPerm` is its standard imprimitive action. Conjugation through `e.symm` has the correct direction. The permutation-form equality, abstract isomorphism, interval blocks and order formula all match the intended structure.

7. **Primitive branch:** Within the odd, coprime regime, the `¬ α ∣ n` branch states primitivity, containment of Mathlib’s alternating group, the exact alternating-group equality criterion, and equality with `S_n` otherwise. No additional hypothesis weakens it.

8. **Amendment:** This is a minimal scope repair: the orbit assertion remains unrestricted, while the dichotomy acquires the intended B.3/B.4 hypotheses. Neither oddness nor coprimality can simply be dropped from that common regime. “Otherwise” must remain scoped inside it. B.3’s explicit block-size classification is not separately restated in the amended Claim or theorem; the supplied permutation-action identification determines it. This is not a missing amended-Claim clause.

9. **Inputs:** The off-line hypothesis matches the quoted first sentence exactly. Jones’s hypothesis only weakens the quoted conclusions; the semilinear formulation correctly expresses the projective-semilinear upper bound. Neither input includes the new orbit, wreath, primitivity or alternating-containment conclusions. The extra B.2 conjunct causes no faithfulness gap.