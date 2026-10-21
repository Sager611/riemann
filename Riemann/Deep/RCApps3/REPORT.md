DISCOVERY: threshold A > (3+√69)/10 = 1.1307 (was 3/2): for every such A, all but O(X^{1−c}) negative fundamental discriminants D with |D| ≤ X have primes in every ring class in every interval (x, x + x|D|^{−κ}] with x ≥ |D|^A; so every positive definite primitive binary quadratic form of discriminant D represents a prime ≤ 2|D|^A, and for all but O(N^{1−c}) squarefree n ≤ N the least prime a² + nb² is ≤ 2(4n)^A (for density-one n, Baier–Zhao 2007 already give < 4n). GRH would give A = 1.

# RCApps3 report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Brief: `.scratch/orca/rcapps3-brief.md`. Paper proofs, not refereed, no Lean. Unconditional.
Review: correct in substance; its five corrections are applied in the NOTE (header lists them). Novelty: [NOVELTY.md](./NOVELTY.md), NEW.

## Theorems (complete proofs in the NOTE)
- **Thm 1 (family large sieve, new tool).** For the coefficient sequences b of a new zero detector (squarefree n, mollified, with the pole of the trivial character removed): Σ_{X/2<|D|≤X} Σ_{χ∈Cl(D)^} |Σ_n b(n)λ_χ(n)|² ≪ X^{O(δ)}(X^{3/2} + N³/X)·N^{1−2σ} for N ≤ X^{1−4δ}. X^{3/2} is the size of the family; adding the per-D large sieves gives X·N·N^{1−2σ}, i.e. Thm 1 saves N/√X for √X < N < X^{5/6}.
- **Thm 2 (density on average over D).** Σ_D #{zeros of L(s,χ), χ ≠ 1, σ ≤ β ≤ 9/10, |γ − t| ≤ 1} ≪ X^{(13−8σ)/(7−2σ)+ε}. At σ = 5/6: 19/16, where RCApps2 and RCApps both gave 5/4 (density hypothesis on average: 7/6).
- **Thm 3 (least primes, D < 0).** RCApps Thm 2 and its corollary for a² + nb² hold for every A > A₀ = (3+√69)/10. A₀ is where Thm 2 and the old family bound X^{10(1−σ)/(3−2σ)} meet the Markov bar 1 + A(1−σ), at σ = 5 − √69/2 = 0.8467.

## How it works
- Positivity: the family over fundamental D is at most max h(Δ) ≈ √X times the same sum over all orders of all discriminants Δ ~ X.
- Dual identity (Lemma 2): over all orders, Σ_Δ Σ_classes |Σ_{𝔞∈C} b(N𝔞)|² is a plain lattice-point count {(T, y): y² | T² − 4nn′, (T² − 4nn′)/y² a discriminant in range}, twisted by Hecke relations: no class numbers, no squarefree sieve.
- Poisson in T, then Mellin in y with exact local densities. The secondary terms are bounded by Heath-Brown's quadratic large sieve.
- The main term is √X·|Σ_n b(n)ψ(n)(n/q)|². The trivial character's pole is killed by the detector factor 1 − 2^{2−2s} (the rational ideal (2) is principal in every order), and the other twists by Jutila's mean value.
- Only loss: pairs n, n′ sharing a factor w > X/N, where the T-intervals are shorter than the modulus; they are bounded trivially (the N³/X term).
- The class I coefficients of the detector are independent of D (squarefree Dirichlet series 𝒢 = L·E·H·M̃), so Thm 1 applies to them; class II uses the per-D mean value. Then Markov over D, as in RCApps.

## Checks
- Dual identities exact on 28 test pairs over all discriminants in [20, 1500]. Local densities verified exhaustively. The Poisson main term matches brute force.
- Brute-force family sums over all orders (X up to 80000, σ = 5/6):
  - Pole-killed sequences: 𝒬/(X‖b‖²) ≈ 0.12, flat for N from 500 to 1.5X, the same as random signs.
  - Unmollified sequences with the pole grow like N/√X; the mollified one grows too (its pole coefficient is a truncated Möbius sum, not estimated).
- `thresh.py` gives the crossover σ = 0.846688 and A₀ = 1.130662; a grid sup agrees.

## Novelty
- Adversarial check ([NOVELTY.md](./NOVELTY.md)): NEW. Closest: Ditchen's large sieve for varying fields (N + N^{1/2}Q^{5/2+ε}), Pascadi–Thorner 2025 (here N + X^{7/2}), Young 2024 (sharp, but Eisenstein series); all are weaker than the sum of the per-D sieves for N ≤ X. Known pieces: the pair count behind Lemma 2 (Soundararajan 2000; summed over d by Heath-Brown–Pierce 2017); for a² + nb² alone, Baier–Zhao 2007 (prime < 4n for density-one n).
- arXiv API: 20 queries answered; 22 more refused with "Rate exceeded".
- Also checked: web search, the GXL evidence API (3 claims), Sardari 1803.03218, Gaudet's thesis, and Ditchen 1312.1502 (|D|^{20/3+ε} outside a density-zero set).
- No family large sieve at the family size and no average density below RCApps were found. The best almost-all exponent was 3/2 (RCApps2).

## Not done / open
- Sectors (angular characters) and D > 0 (indefinite forms) stay at 3/2. Their dual counts need angles, or unit-reduced hyperbolic lattice points.
- A₀ → 1 (GRH strength) needs Thm 1 without N³/X up to N = X^{1−δ}. Numerically the family sum stays at the family size there (N = 1.5X). A proof would need cancellation in Salié-type sums over T² ≡ 4nn′ (mod y²) for intervals shorter than the modulus.
