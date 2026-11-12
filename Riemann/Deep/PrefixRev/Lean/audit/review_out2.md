VERDICT: FAITHFUL

The supplied signatures resolve the round-1 hypothesis blocker.

- **Lemma 2:** All three identities carry the size, ordering, and modular-inverse hypotheses. Since `n+1 = (k+1)+(m+1)`, the inverse condition also forces `gcd(k+1,m+1)=1`.
- **Coordinates:** `phi` uses the required 1-based input `x+1`. The formulas give cyclic increment, full reversal, and `n+2−j` on the prefix image, respectively. `↑x < k` correctly identifies that prefix.
- **Lemma 3:** `block_classes` classifies every nonempty block: singleton, universal, or a complete residue class of `psi = phi−1`, with `d ∣ n` and `1<d<n`. It does not assert the unjustified converse that every such class is a block.
- **(★ₙ correction: ★_d):** `star_key` excludes mixed prefix membership in every nonzero `psi` class—exactly every `phi` class except `1 mod d`. `star_of_block` explicitly assumes odd `n`; then `d ∣ n` supplies the oddness of `d` required by `star_key`. This is an odd-`n` derivation, not an unrestricted one.
- **B.1:** `G_le_cls` asserts preservation of each `cls` value, hence invariance of every corresponding class.
- **`fibreFix`:** Yes—setwise preservation of each fibre, not pointwise fixation of its elements.
- **`fibreStab`:** Yes—preservation of the fibre equivalence relation. A permutation may exchange nonempty fibres; it need not preserve their labels.
- **B.4:** The revised docstring accurately limits the claim to the alternating branch.

No tools used; file not written under the read-only restriction.