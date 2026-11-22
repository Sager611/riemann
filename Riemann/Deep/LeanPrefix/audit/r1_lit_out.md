VERDICT: FULLY-PROVED-AND-FAITHFUL

Remarks

- **Scope:** the two literature inputs. Both follow from the supplied source statements; neither imports a step claimed as ours. The main theorem remains conditional on these named literature hypotheses.

- **LZZW:** hypotheses exactly preserve `2 ≤ k < m < n` and `m + k < n - 1`. They ensure positive n, so natural subtraction agrees with integer subtraction. The identification `j = x.val + 1` converts `rev n i` exactly into the source’s reversal. All three reversal lengths are in range; the `i > n` identity convention cannot intervene. Since `G` is already a subgroup of `S_n`, `G n m k ≠ ⊤` expresses properness exactly.

- **Jones, antecedents:** `IsPreprimitive` has the supplied meaning of transitivity plus trivial blocks. `Perm.IsCycle c ∧ c.support = Finset.univ` expresses a single cycle moving every point. The missing explicit `n ≥ 2` causes no strengthening: for n = 0, 1, 2, the alternating group is trivial, making the noncontainment hypothesis impossible. For n = 3, a full cycle generates A₃; for n = 4, a primitive subgroup containing a 4-cycle is S₄. Neither gives a counterinstance.

- **Jones (a):** one common bijection followed by elementwise affine formulas expresses containment in AGL₁(p). With `n.Prime`, nonzero elements of `ZMod n` are precisely the possible invertible multipliers. Dropping the cyclic lower bound only weakens the conclusion.

- **Jones (b):** a finite field of cardinality q supplies q ≥ 2 and prime-power cardinality; consequently the natural-number degree formula agrees with the source. `IsCollineationVia` expresses a semilinear bijection inducing the required projective permutation. Its line-inclusion condition suffices: an injective semilinear map sends a nonzero line onto the target line. Elementwise choices of lifts are appropriate; no coherent choice of lifts is required.

- **Hyperplanes:** no omitted case. The action on hyperplanes becomes the projective action on the dual space through annihilators. The contragredient action is semilinear, and choosing coordinates identifies the dual with `Fin d → F`. The existential bijection therefore accommodates both point and hyperplane actions.

- **Jones (c):** retaining only `n = 11 ∨ n = 23` weakens the three listed exceptional alternatives.

- **Citations:** the supplied theorem numbers, source extracts and bibliographic identifications agree with the inputs. LZZW is explicitly marked an **unrefereed preprint**; Jones is marked published.

- **No hidden bespoke input:** neither Prop assumes prefix-group primitivity, the fixed-point bound, the constructed n-cycle, exclusion of Jones’s alternatives, B.1–B.3, or the line necessity directions. Given the supplied compilation evidence, those conclusions are derived rather than inserted into these literature hypotheses.