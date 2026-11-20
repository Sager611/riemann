# PROBES — alpties (20 Nov 2026). Exact arithmetic; code/ and logs/ in this folder. memcap.sh 290 MB wrapper on every run (its sampler reported peak 0 MB: the runs were short, < 30 s, pure-Python/sympy on tiny inputs).

## A. InducedSimpleLowestAlcove is Jantzen's closed-bottom-alcove criterion, and it is sharp where the brief predicts (code/p_alcove.py, logs/p_alcove.log)
- Archive statement (Development/AbelianVarietyTheory/PublishedTheorems/GlobalTheta/WeylDualWeylLatticesAgreeOfLowestAlcove.lean:19-25): for p prime, k of characteristic p, a dominant weight λ of GL_N with λ_i − λ_j + (j − i) ≤ p for all i < j, the polynomial model ∇(λ) = nablaSq k N λ ⊂ k[x_ij] is ≠ ⊥ and has no GL_N-submodules except ⊥ and itself (IsGLSubmodule: stable under GL_N(A) for every k-algebra A, i.e. as an algebraic-group module).
- Method: λ = (n, 0, …, 0), ∇(λ) = S^n(V). Weight spaces are 1-dimensional (monomials), so a G-submodule is a span of monomials closed under the divided-power root operators E_ij^(k): x^a ↦ C(a_j, k)·x^{a + k(e_i − e_j)}. S^n(V) is simple iff the monomial graph with edges where C(a_j, k) ≢ 0 mod p is strongly connected. Exact integers.
- Cases: p ∈ {2, 3, 5, 7}, N ∈ {2, 3, 4}, 0 ≤ n ≤ 2p (dimension ≤ 400): 111 cases.
- Results:
  - condition ⇒ simple: 0 counterexamples.
  - GL_2, n = p (λ_1 − λ_2 = p, one step outside the closed alcove): never simple (p = 2, 3, 5, 7); the span of x^p, y^p is a submodule (Frobenius twist).
  - outside but simple: 16 cases (e.g. GL_2 n = 2p − 1; GL_3 n = p − 1, where λ_1 − λ_3 + 2 = p + 1). The condition is sufficient, not necessary.
  - This is exactly the direction of Jantzen, RAGS II.5.6 (H^0(λ) = L(λ) for λ in the closed bottom alcove).
- Verdict: PASS; the formal condition is the right one (⟨λ + ρ, α∨⟩ ≤ p for every positive root α = e_i − e_j).

## B. Toy analogue of MomAt: jets of a CM kernel at split vs inert primes (code/p_jets.py, logs/p_jets.log). NOT a faithful instance
- MomAt (BlockTheory/GlobalSignature/JetSumTable.lean:99-120) needs A.thInd, an ambient built from a Base, D7 data and a selection; no finite instance fits under the cap.
- Toy: ℘ of Ω·Z[i] in the integral normalisation of y² = 4x³ − 4x; jets D^a(℘ − z^{−2})(0) = a!·c_{a/2}, a ≤ 600.
- Split p = 5, 13, 17 (the type of the level primes T_u, probe C): min_a v_p(a!·c)/log_p a = −1.11, −1.02, −1.01. So the jets are integral up to a p^{O(log a)} denominator and show no growing divisibility. A moment bound "jets lie in the integral closure of a unit-period table" is the right shape there (Katz-type p-adic periods at ordinary CM points).
- Inert p = 3, 7, 11: v_p(a!·c) ≥ Σ_{odd i} a/p^i − 0.76…0.87·log_p a for every a. The growing divisibility lives only at inert primes (agrees with AlpCore Test A).
- Use: the m-step chain's p-adic input (MomAt, KerBar, KerP) sits at split level primes, where only integrality (no gain) is expected. This is consistent with §C and TIES §2.

## C. The Split axioms force the level primes to split in E/E⁺ (code/p_split.py, logs/p_split.log)
- Built.Split (Definitions/GlobalTwisting/KerIsTerm.lean:23-37) requires T ∈ P, T ∉ P², and nonconj: Q ≠ conjPrime₂ P for P, Q ∈ ST. With Q = P this gives P ≠ P̄.
- Probe: E = Q(ζ4), Q(ζ3), Q(ζ12), primes T < 200 not dividing n. Two independent tests of P = P̄:
  - Galois: −1 ∈ ⟨T mod n⟩;
  - factor Φ_n mod T and test each irreducible factor against its reciprocal.
- Results:
  - 0 mismatches between the two tests, over 134 (field, T) pairs.
  - ST can be nonempty only for T ≡ 1 (4), T ≡ 1 (3), and T ≡ 1, 5, 7 (12) respectively.
  - The primes T that are inert in E/E⁺ (T ≡ 3 (4), 2 (3), 11 (12)) are excluded.
- Verdict: PASS. KerBar (seenP w = false, i.e. w ↔ P̄) and KerP (w ↔ P) live at the two conjugate primes over a split level prime T_u. OverT w := resChar w = T (Definitions/GlobalTwisting/LevelChainSharpColumnValuationOwnChainIff.lean:79). They never live at an inert prime. Alprows' gloss "conjugate (inert)" is a misreading.

## D. Zero-kernel screen of the five kernel conditions (by hand from the definitions, KerIsTerm.lean:86-113)
- k = 0 satisfies KerEq.equiv (0 = ρ'∘0∘ρ⁻¹), KerDil (0 = 0), KerUK (0 = T^{⌊(t+1)/2⌋}·0), KerBar (take F' = 0 ∈ valAt) and KerP (0 ∈ L).
- Only KerEq.kmap pins k: Kmap i c c' (ev_i F) = ev_{i+1}(k … F), where Kmap = (rhoDet·rhoDet')•(divDe ∘ unf ∘ sharp) (Formalism/LocalRingCardinality/DataKernels.lean:63-65).
- So laws_true_fam is a theorem about the actual chain kernel kAt = KerT.ker (a theta-series kernel over Hom(E^{N_i}, E^{N_{i+1}}) with test function Ξ; Definitions/GlobalTwisting/BelowFiniteNonempty.lean:20-23), exactly when Kmap ≠ 0. That Kmap ≠ 0 at built data was not checked (open).

## E. hDivG / BkLaw: no small theta-lift instance fits
- DivRows3G (Formalism/LocalToGlobal/Divisibility/TupUsed.lean:112-126) relates a RowSideFamily: the indefinite realisation A3 with θ_{i,t}, lPrj, eMt, jMt, colData. This family is built from a D7 whose m lifts have ranks N_i = k + i(t+1) on host fields of degree 2^e. Even m = 2 needs a full Base (CM host, frames, spaces sp.X i), so no instance fits under the cap.
- BkLaw D := RationalMapValues D → BuRational D (NonzeroPlacesOverPrimesLifts.lean:107) is checked logically instead. If every value Bu u (jN φ)(jEE x) lies in embLL(LL), then choosing Bk u φ x as the preimage gives BuRational. So BkLaw is an identification of Bk with the LL-descent of the pairing Bu. It holds automatically when Bk is defined by descent (BkDef, WeightDegreeSizes.lean:1267-1269): PLUMBING.
