DISCOVERY: two thresholds drop from 3/2 to A₀ = (3+√69)/10 = 1.130662…: (1) indefinite forms / narrow classes, D > 0 (Theorem 3⁺); (2) sectors, D < 0, incl. primes a² + nb² with arg(a+ib√n) in a short arc (Theorem 3^∠). Complete paper proofs in [NOTE.md](NOTE.md). Goals 3 (a_E signs, A > 3 and 9/4) and 4 (sharper family sieve): no improvement.

# RCApps4 report (2026-10-21)

## New results
- **Theorem 3⁺ (D > 0).** For every A > A₀, all but O(X^{1−c}) fundamental D ∈ (0, X] have, in every narrow ideal class, prime ideals of
  norm in (x, x + xD^{−κ}] for every x ≥ D^A; so every primitive indefinite binary quadratic form of discriminant D represents a prime
  ≤ 2D^A, and a² − nb² represents a prime ≤ 2(4n)^A for almost all squarefree n ≡ 2, 3 (4). Was A > 3/2 ([R2] Theorem 4).
- **Theorem 3^∠ (sectors, D < 0).** Same A₀ for primes in every ring class with the angle of the generator in any arc of length ≥ |D|^{−κ};
  for almost all n, a prime a² + nb² in (x, x + x(4n)^{−κ}], x ≥ (4n)^A, with arg(a+ib√n) in any arc of length ≥ (4n)^{−κ}. Was 3/2.

## How
- **D > 0.** The obstruction named in RCApps3 §6 (units make the dual count infinite) is removed by averaging over the Hecke
  Grössencharacters λ_k (∞-type e(k·u/R⁺)) with weights ŵ(k/R⁺) ≥ 0. Poisson summation in k unfolds the unit group: the dual becomes
  h⁺R⁺·Σ_{γ≫0, Nγ=M} w(u_γ), and h⁺R⁺ = √Δ L(1,χ)Π(…) ≪ √X log²X has no regulator in it (Lemma 2⁺). The lattice count
  T² − Δy² = 4M, T ≤ 2√M cosh 1, has the shape of the D < 0 count, so RCApps3's Poisson (Lemma 6⁺, kernel K⁺) and Mellin steps go
  through: Theorem 1⁺ (family sieve X^{O(δ)}(X^{3/2} + N³/X)N^{1−2σ}), Theorem 2⁺ (density X^{(13−8σ)/(7−2σ)}), Theorem 3⁺.
- **Sectors.** For ∞-type ℓ the dual count gets the weight T_ℓ(T/(2√M)) (Chebyshev); Markov's inequality keeps Poisson and Mellin
  losses at (1+ℓ²)^{O(1)} ≤ X^{O(δ)} for |ℓ| ≤ X^δ (Lemma 2^∠, Theorems 1^∠, 2^∠); [R] Prop 1′ covers β > 9/10.

## Checks (all exact or to rounding; scripts in .scratch/rcapps4/)
- Lemma 2⁺ orbit form: 534 real quadratic orders × 16 pairs, 8544/8544 exact; wrong variants fail.
- Unfolded identity (Grössencharacters computed from ideal generators): 60 tests, max rel. deviation 8.9·10⁻¹⁶.
- Lemma 6⁺ main term vs brute force: e.g. 44.4781 vs 44.4781, 19.3945 vs 19.3948 (where X′ ≫ c√M).
- Lemma 2^∠: 399 imaginary quadratic orders × ℓ ∈ {2, 4, 10}, max rel. deviation 6.7·10⁻¹⁵.

## No improvement
- **a_E signs (A > 3; one-sided 9/4).** At σ = 5/6 the bottleneck is class II (second moments of degree-4/6 L-functions, |D|^{(m+1)/2})
  or polynomials longer than X; the family sieve only covers lengths ≤ X. Best detector bound: X^{9/4−σ} (m = 1, Y = X) or X × Theorem C;
  X^{5/3} for m = 2 against the bar X^{3/2}.
- **Sharper family sieve.** The N³/X term comes from pairs sharing w > X/N: both Poisson main term and error are N³X^{−3/2} there; a
  per-Δ large sieve gives the same; needs cancellation in sums of Salié-type sums (half-integral weight Kuznetsov). Open.

## Novelty
Evidence API (3 claims) and earlier adversarial audits ([R2], [R3] NOVELTY): nothing on indefinite forms on average over D (any
exponent), nothing below 3/2 for sectors; the closest are Ditchen 2013 (definite, 20/3, density zero) and fixed-field sector theorems.
The unfolding device is classical in spirit; its use to make the D > 0 family dual regulator-free appears new.

Status: paper proofs, not refereed, no Lean. Suggested DISCOVERIES.md edit: in entry 3, "1.1307" now also for D > 0 and for sectors.
