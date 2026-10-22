DISCOVERY: Fix a newform f without CM (e.g. f = f_E, E/Q without CM). For all but O(X^{1−c}) fundamental discriminants D with |D| ≤ X, negative or positive, every (narrow) ideal class of Q(√D) contains, in every interval (x, x + x|D|^{−κ}] with x ≥ |D|^A, degree-one primes p with λ_f(p) > −α (any fixed α > 0) as soon as A > 1 (was 9/4; 1 is the GRH value of the method), primes with λ_f(p) > 0 and primes with λ_f(p) < 0 as soon as A > (23 − 3√5)/8 = 2.0365 (was 3), and primes with cos θ_p in any interval J admitting a degree-3 Sato–Tate minorant as soon as A > 3 (was 10/3; uses Pascadi–Thorner). In particular every primitive binary quadratic form of discriminant D represents primes p, p′ ≤ 2|D|^A with a_E(p) > 0 > a_E(p′) for every A > 2.0365. Thresholds for minorant degree ≥ 4 (11/3, (M+3)/3) and the discrepancy 10/(3A−3) are unchanged.

# STLemP report (2026-10-22)
Proofs: [NOTE.md](NOTE.md); checks: [checks/](checks/). Brief: `.scratch/orca/stlemp-brief.md`. Paper proofs, not refereed, no Lean, no 2026 half-planes.

## Results (complete proofs in the NOTE)
- **Theorem 1′.** The RCApps5 family sieve (Lemma P) holds for every sequence whose quadratic twists are controlled (class 𝔅): the detector sequences, the Sym^m f-twisted detector sequences and the approximate-functional-equation pieces of L(s, Sym^m f ⊗ θ_χ), m = 1, 2. Not for arbitrary coefficients (the pole term is exactly what the class excludes).
- **Lemma L (ring-class lift, new).** Class group characters of D lift to the ring class group of the order of conductor p₀; Lemma P makes all other orders of discriminant ≍ p₀²|D| nonnegative. Result: the family sieve beyond the conductor, ≪ max{X, N}^{3/2+o(1)}N^{1−2σ} (per-D sum X·N^{2−2σ}).
- **Theorem MS.** Σ_{|D|~X}Σ_{χ∈Ĝ(D)}|L(½+δ+it, f⊗θ_χ)|² ≪ X^{3/2+O(δ)}: Lindelöf on average for the double family (per-D sum X²). For Sym²f: X^{9/4} (was X^{5/2}).
- **Theorem 2₁.** Σ_{|D|~X}Σ_χ N(σ, f⊗θ_χ) ≪ X^{2−σ+O(δ)}: the density hypothesis of each D, on average (was X^{11(1−σ)/(3−2σ)} and the per-D Theorem C). **Proposition F₂** (Sym²f): X^{15(5−4σ)/(16(2−σ))}, ¾ ≤ σ ≤ 13/14.
- **Theorem D′, Corollary D′1, Sato–Tate table** 𝒜(1) = 1, 𝒜(2) = (23−3√5)/8, 𝒜(3) = 3, 𝒜(M ≥ 4) unchanged; D > 0 (narrow classes) included.

## What set the old thresholds (task 1)
- 9/4: RSFam's per-D Theorem C against its Proposition F at σ = 5/6; Proposition F carried the per-D second moment I₂ = X², although the AFE length |D| is the conductor, where Lemma P works. 3: Proposition F (HT threshold X^{7/2} × I₂ = X^{5/2}) meeting the trivial count at σ = 5/6. 10/3: ChebST's route H′ was evaluated at a fixed Y; with Y free it gives 3 (bookkeeping, not Lemma P).
- None of these contained RCApps3's N³X^{−1} term: neither RSFam nor ChebST used RCApps3 Theorem 1.
- The m = 2 value is where Proposition F₂ meets the bar: σ* = 1 − (1+3√5)/44 = 0.8248 (class II with I₂ = X^{9/4} against the lifted class I at the top length Y = X^{15/(8(2−σ))}, mollifier length Y^{1/4}).

## Limits (NOGO with present tools)
- 11/3 (degrees 4–8): class I at N = X^{11/4}, where Pascadi–Thorner meets Humphries–Thorner; X^{7/4} beyond the conductor; the lift costs N^{3/2} there.
- Slope ⅓ per degree (M ≥ 8): the second moment; the lift gives X^{3(m+1)/4}, no gain for m ≥ 3. Below A = 1: the explicit-formula barrier.
- Numerically (37a, X = 20000) the twisted family sums stay at the family size up to N = 4X; a proof would give 𝒜(2) = 3/2, 𝒜(M) = (M+3)/3.

## Checks
51 exact sympy checks of every exponent (0 failures); Lemma L on 43 (D, p₀); brute-force family sums over all orders for λ_E, λ_{Sym²E} (5 lengths × 5 sequence types); parameter optimiser reproducing the old values 9/4, 3, 10/3, 11/3, 13/3.

## Novelty
arXiv math.NT listings 2026-09 (812) and 2026-10 (625) scanned, 8 abstracts read, 6 arXiv searches: nothing on Hecke-eigenvalue or Sato–Tate conditions on primes in ideal classes of varying fields, on family mean values of L(½, f⊗θ_χ) over D and χ, or on class-group family sieves beyond the conductor.

## Significance (for the reviewer)
Records for this project's own theorems (RSFam, ChebST, audited NEW there): one-sided at the GRH value 1, sign from 3 to 2.0365. New tools of independent interest: the ring-class lift and the double-family Lindelöf on average for f ⊗ θ_χ. The 𝒜(3) gain is a re-optimisation of ChebST (could go to COROLLARIES.md).
