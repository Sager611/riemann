DISCOVERY: A > 1 (was (3+√69)/10 ≈ 1.1307; 1 is what GRH gives by the same method): for every A > 1, all but O(X^{1−c}) fundamental discriminants D with |D| ≤ X, negative or positive, have prime ideals in every (narrow) ideal class in every interval (x, x + x|D|^{−κ}] with x ≥ |D|^A, also in prescribed sectors for D < 0; so every primitive binary quadratic form of discriminant D, definite or indefinite, represents a prime ≤ 2|D|^A.

# RCApps5 report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Checks: [checks/](checks/). Brief: `.scratch/orca/rcapps5-brief.md`. Paper proofs, not refereed, no Lean. Unconditional (no 2026 half-plane, no Kerr).

## What fixed 1.1307 (task 1)
- Markov over D needs Σ_D N_D(σ) < X^{1+A(1−σ)}. RCApps3's bound X^{(13−8σ)/(7−2σ)} needs A > 6/(7−2σ); RCApps's X^{10(1−σ)/(3−2σ)} needs A > 10/(3−2σ) − 1/(1−σ). They cross at σ = 5 − √69/2 = 0.8467.
- The first bound came from the N³/X term of the RCApps3 family sieve at N ≈ X^{0.94}. Its only source: RCApps3 summed every non-maximal order over ideals prime to its conductor, which forces a Möbius inversion over common factors w of n, n′; for w > X/N Poisson summation fails. No zero-free input was involved.

## The fix (new, complete proofs)
- Lemma P: in every quadratic order, count pairs over all elements, with no coprimality to the conductor f. Dividing by e (γ ↦ γ/e) splits this count into a sum over squarefree e | f of class-group sums of the order of discriminant Δ/e². So it is ≥ 0, dominates what RCApps3 used, and its dual is RCApps3's lattice count without the Möbius terms.
- Thm 1: the large sieve over pairs (D, χ ∈ Cl(D)^) is at the family size X^{3/2} for every length N ≤ X^{1−ε} (RCApps3: + N³/X).
- Thm 2: Σ_D N_D(σ) ≪ X^{2−σ+ε} for ½ < σ ≤ 9/10: the density hypothesis of each D, on average over D (RCApps3: (13−8σ)/(7−2σ)).
- Thms 3, 3⁺, 3^∠: the least-prime theorems of RCApps, RCApps2 and RCApps4 hold for every A > 1: every class for D < 0; narrow classes, indefinite forms and a² − nb² for D > 0; arcs and a² + nb² for sectors. For 1 < A ≤ 1.2: δ = (A−1)/8000, κ = δ/4, c = δ/2.

## Checks
- Lemma P exact on 1830 (Δ, b) tests, 1449 of them non-maximal orders: D < 0 (blocks also recomputed from reduced forms), sectors ℓ = 2, 4, 6, D > 0 with a smooth unit weight. Deviation ≤ 6·10^{−14}, no negative block.
- RCApps3's Poisson main term matches brute force for M with square factors, the pairs Lemma P keeps.
- Brute-force family sums over all orders, X = 20000 and 40000, σ = 5/6:
  - pole-killed sums stay at the family size from N = X/16 to N = 1.5X;
  - Lemma P's sum exceeds RCApps3's by about 1%;
  - sequences that keep the trivial character's pole grow like N/√X.

## Other routes (task 2)
- DFI amplification, Heath-Brown refinements, fourth moments, per-D log-free, Huxley or hybrid inputs: not needed, and none can go below 1.
- Liu's half-plane (0.8749): no effect, since the binding point is now σ → ½.
- Kerr (no Siegel zeros): only removes the Siegel factor ν_D, giving ≫ H/(h log x) primes in every class uniformly. COROLLARY (unrefereed input).
- Below 1: NOGO for this method. Zeros near ½ enter with |x^{ρ−1}|, and h·x^{−1/2} ≤ 1 needs x ≥ |D|. GRH gives exactly A = 1 this way.

## Novelty
- Scanned the arXiv math.NT listings for 2026-09 (812) and 2026-10 (590), and ran 16 arXiv searches. The GXL evidence API refused (429, daily quota).
- Nothing states or implies an almost-all exponent below 1.1307. The best in the literature is Ditchen's 20/3 (definite forms, density-zero exceptions). For a² + nb² alone, Baier–Zhao 2007 give a prime < 4n for density-one n.

## Corrections and pointers
- RCApps3 §6 and RCApps4 §4 called the N³/X term intrinsic (needing Salié-sum cancellation). It was an artifact of the conductor inclusion–exclusion.
- CM2 (CM fields over a totally real F, threshold 1.1307) has the same Möbius terms. Lemma P looks transferable (A → 1 there too); not checked.

## Review (2026-10-22)
- Reviewer: correct as proved. Theorem 1 holds for the detector sequences only; for arbitrary coefficients it is false (b = 1 on primes: the trivial characters give ≫ Xπ(N)²).
- 8 independent scout checks passed (Lemma P tests, Poisson check, §4/§5 arithmetic, D > 0, sectors).
- Wording fix for NOTE.md l. 74–76: read "a principal ideal of norm M prime to f_e is (γ)"; principal ideals form the identity-class part, r_{Δ_e}(M) = (w_e/2)·[1](T_M).
