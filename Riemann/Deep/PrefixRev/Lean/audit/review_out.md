VERDICT: ISSUES

1. **Blocking evidence gap — omitted section hypotheses.** As displayed, `phi_P`, `phi_rn`, and `phi_rk` lack necessary assumptions. For example, θ = 0 and (n,m,k) = (7,4,2) make `phi_P` assert 0 = 1. These may correctly inherit section hypotheses; the excerpt does not establish that. Supply elaborated signatures showing the line, ordering, and modular-inverse hypotheses. Likewise check `block_classes` and `star_of_block`. This is not evidence that the actual compiled declarations are wrong.

2. **Faithful — `rev` and `G`.** Position x+1 maps to i−x, hence zero-based x maps to i−1−x. The cutoff x<i is correct. Identity for i>n does not affect the stated parameter regimes. Subgroup closure gives precisely the generated subgroup.

3. **Faithful — product and cycle coordinates.** `P = r_k r_n r_m` applies r_m first. `P_val` and `P_firstReturn` translate Lemma 1 correctly. The displayed coordinate identities match Lemma 2 when their necessary ambient hypotheses hold.

4. **Faithful — published statements.** `Conjecture4Literal` matches Conjecture 4’s first sentence with 2≤k<m<n. `OffLineLZZW` matches Theorem 3.1’s first sentence; natural subtraction causes no discrepancy. Neither encodes the respective “Moreover” sentence. All three displayed counterexamples satisfy the conjecture’s standing assumptions.

5. **Faithful — `JonesK0`, with the intended `Coll`.**
   - (a): a nonidentity affine map fixes at most one point.
   - (b): prime-power q≥2 makes the multiplied degree equation equivalent to Jones’s quotient formula.
   - (c): the bounds apply to **all** nonidentity elements, not merely involutions. Degree-11 M₁₁ is sharply 4-transitive; L₂(11) has maximum three fixed points; degree-23 M₂₃ has maximum seven.
   Thus these are valid weakenings. The absence of an explicit degree lower bound is harmless: failure to contain A_n excludes degrees 0 and 1.

6. **Faithful — `Lemma7`, conditional on `Coll`’s intended meaning.** `2 * nfix h ≤ n` expresses the real bound nfix(h)≤n/2 exactly. The exceptional pair is correct. Valid finite-field parameters and the natural projective action must come from `Coll`; its abstract type alone supplies neither.

7. **Faithful restricted form — `StarD` / `Prop5Core`.** Subtracting one makes the exempt class j≡1 become ψ≡0. Membership x.val<k correctly pulls back E. On the line, (k+1)|n iff (k+1)|m. `Prop5Core` restricts Proposition 5 to the odd, coprime setting needed here; it is not its unrestricted formulation.

8. **Faithful — primitivity.** `IsBlockOf` is the usual block condition. Transitivity plus “subsingleton or universal” correctly includes empty and singleton trivial blocks and matches the stated preprimitive convention.

9. **Faithful conclusions, conditional proofs.** `theoremA`, `theoremA_core`, and `theoremA_jones` reproduce Theorem A’s equivalence. Their named assumptions remain obligations; standard kernel axioms do not discharge them. No supplied evidence proves `Prop5Core` or `Prop5Primitive` false, but neither is established here.

10. **Minor — B.4 coverage/documentation.** `line_B4` gives both group-equality branches. `line_B4_jones` gives only the alternating branch under `h4`; its “exactly when” comment overstates that individual declaration. Neither declaration itself includes primitivity.

11. **Minor — B.1/B.3 coverage and missing definitions.** B.1’s claimed invariance/nonemptiness is appropriately weaker than exact orbit identification. B.3’s displayed results do not assert the wreath-product isomorphism. To verify the subgroup inclusions literally, supply `fibreFix` and `fibreStab`: the latter must stabilize the partition, not fix each interval individually.

12. **No identified intrinsic vacuity.** The named hypotheses are conditional inputs, not axioms. `JonesLemma7` does not fail at the excluded projective case: PGL₂(4)=A₅ already contains the required alternating group.

No tools used; requested output file not written.