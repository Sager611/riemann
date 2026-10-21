NO DISCOVERY: the ring class zero-density theorem extends to three wider thin families (new statements, same method, routine); for general number fields its key lattice lemma fails in degree 3; no unconditional least-prime bound follows.

# Density2 report (2026-10-20)
Proofs: [NOTE.md](NOTE.md). Base: [Combos NOTE](../Combos/NOTE.md), Theorem 1 (Σ over ring class characters of D of N(σ;t) ≪ |D|^{2σ(1−σ)+ε}). Unconditional unless marked (L) or (Hu).

## New theorems (complete proofs; new as stated; routine)
- **A. Anticyclotomic characters** of an imaginary quadratic field (ring class characters times angular characters; their L-functions are those of CM newforms of weight ≤ M+1): Σ N(σ;t) ≪ (|D|(M+3)²)^{2σ(1−σ)+ε}T^{1+ε}. M = 0 is the old theorem, now with T^{1+ε} per window instead of T³. For the Gaussian angular characters λ_k, |k| ≤ K: K^{4σ(1−σ)+ε}, better than Ricci's K^{(10/3)(1−σ)} (1976) exactly for σ < 5/6.
- **B. Narrow ring class characters of a real quadratic order** (discriminant D > 0): Σ N(σ;t) ≪ D^{2σ(1−σ)+ε}T^{1+ε}. Beats the trivial bound when the class number h⁺ exceeds D^{2σ(1−σ)}, e.g. for D = n² + 4.
- **C. All finite-order Hecke characters mod 𝔪** of a quadratic field (q = least positive integer in 𝔪): ≪ φ(q)(q²|d_K|)^{2σ(1−σ)+ε}T^{1+ε}. For bounded q this is |d_K|^{2σ(1−σ)+ε}; the previous note had (|d_K|N𝔪)^{3(1−σ)} (imaginary K).
- Tools: every sector of angle φ ≤ π and radius R holds at most 1 + R²φ/covolume primitive lattice vectors (Lemma 1), plus a hyperbolic analogue (Lemma 3). They give weighted large sieves with constants ≍ L_D(N + |D|^{1/2}(M+6)) and L_D(N + D^{1/2}), both essentially sharp (Lemmas 4–5). A coset reduction handles ray classes.
- Checks: qa1–qa3.py (lattice counts and exact large sieve norms, all within the bounds), cubic.gp.

## Item (b), general number fields: no unconditional result (NOGO)
- The primitive-ideal lemma fails in degree 3. For K = Q(d^{1/3}), d = m³ + 1, the m + 1 principal ideals (θ − x) are not divisible by any integer, have norm ≤ d, and lie in one class, while h_K ≫ d/log²d. So any such large sieve loses a factor ≫ d^{1/3}/log²d (Prop 7; gp: the ratio reaches 1.5 at d ≈ 10⁵).
- What would suffice is an L² bound for pairs of ideals in the same class (Prop 8, conditional). It holds for quadratic fields and is open in degree ≥ 3.

## Least primes
- Unconditional: none. All bounds lose a factor 𝔔^ε (𝔔 = the analytic conductor), which is fatal within O(1/log 𝔔) of σ = 1, where Linnik-type arguments need log-free density.
- Conditional on (L): (i) every primitive indefinite binary quadratic form of discriminant D > 0 represents a prime ≪ min(D^{1.7498}, (h⁺ log D)^{7.99}); before, (L) gave up to D^{4.0}, and the unconditional bound is D^{O(1)}. (ii) For bounded q, every ray class mod 𝔪 contains a prime ideal of norm ≪ |d_K|^{1.7498} (the previous note had |d_K|^{3+ε}, imaginary K).
- Under (Hu): a ring class contains primes whose generators lie in a given angular sector; weaker than Ricci's unconditional sector bound for Q(i).

## Side result
- The "optimal" hybrid large sieve N + MT for Hecke characters of Q(i), which Stucky (arXiv 2008.11325) calls unavailable, is false. Coefficients on the ideals (a + i), a in a short interval, exceed it by a factor ≫ N^{1/4} (Remark 9; numerically 67 at N = 10¹², 661 at N = 10¹⁶). So the thin-family large sieve has no t-aspect version.

## Novelty
Pascadi–Thorner 2508.14888 Thm 1.3 (any GL₁ family) is weaker in |D| and in M for every σ. In the angular aspect: Ricci 1976, Coleman 1990/1992, Kubilius 1955, Chen–Park–Swaminathan 1506.09170. Nothing found for real quadratic class group families. Searches: arXiv API (30 queries), OpenAlex (10), GXL (4 claims); 6 papers read (NOTE §7).

## Leads
- A log-free version of Theorem A (M = 0) near σ = 1 would have exponent ≈ 2 where Zaman's thesis has 162, the stated bottleneck of the unconditional record |D|^455 for the least prime represented by a form. This needs explicit constants throughout (Zaman, Ch. 5–7).
- Prove the L² pair bound for CM fields over a fixed totally real field (rank-2 lattices over its ring of integers).

Significance: Theorems A–C are routine extensions. Suggest COROLLARIES.md, not DISCOVERIES.md.
