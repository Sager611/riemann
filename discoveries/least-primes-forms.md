# Least primes represented by binary quadratic forms, for almost all discriminants

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim (unconditional).** Fix A > 1. For all but O(X^{1−c}) non-square discriminants Δ ≡ 0, 1 (mod 4) with |Δ| ≤ X (fundamental or not, negative or positive), every ring class of the order of discriminant Δ (narrow ring class if Δ > 0) contains degree-one prime ideals in every interval (x, x + x|Δ|^{−κ}] with x ≥ |Δ|^A. For fundamental D < 0 this also holds in prescribed sectors. The constants c, κ > 0 depend on A; intervals have length H ≤ x.

**Consequences.**
- Every primitive binary quadratic form of discriminant Δ, positive definite or indefinite, represents a prime ≤ 2|Δ|^A.
- For almost all n (not only squarefree n), x² + ny² represents a prime in every such interval.

**Sharpness.** For D < 0, A > 1 cannot be improved: the principal class has no degree-one prime of norm below |D|/4. GRH gives A = 1 for every D (Sardari 2018).

**Before.**
- |D|^{20/3+ε} outside a density-zero exceptional set, definite forms only (Ditchen 2013).
- Indefinite forms: only one class on average (ny² − x² represents a prime < n for almost all n, via the prime-plus-a-square theorem), not every class.
- The clause for a² + nb² is known in stronger form (Baier–Zhao 2007).
- Our own earlier steps: A > 1.6536, then 3/2, then (3+√69)/10 ≈ 1.1307 (the last also for D > 0 and sectors).

**Inputs.** Our [ring class zero-density theorem](ring-class-zero-density.md) and its real-quadratic analogue; a long-mollifier per-D density (zeros with β < 3/2 − 1/A); the Humphries–Thorner zero detector over D (narrow class characters give tempered Maass forms). Genus characters and Siegel zeros are handled separately. No 2026 half-plane is used.

**Steps.**
1. A > 1.1307 ([RCApps3](../Riemann/Deep/RCApps3/NOTE.md)): a new large sieve over pairs (D, χ ∈ Cl(D)^) at family size X^{3/2} (all-orders positivity, a dual lattice-point identity with Hecke relations, Poisson/Mellin, Heath-Brown's quadratic large sieve) and a family density X^{(13−8σ)/(7−2σ)+ε}. Reviewed correct in substance (corrections applied); novelty checked ([novelty 3](../Riemann/Deep/RCApps3/NOVELTY.md)).
2. Every A > 1 ([RCApps5](../Riemann/Deep/RCApps5/NOTE.md)): count the pairs over all elements of each order (no conductor Möbius terms; the pair-count device is Soundararajan's and Heath-Brown–Pierce's). This removes the N³/X term and gives the family density Σ_D N_D(σ) ≪ X^{2−σ+ε}, the density hypothesis on average over D. Reviewed correct for detector sequences (the pair large sieve is false for arbitrary coefficients). Adversarial novelty: nothing states or implies any A < 1.1307 for almost all D ([novelty 5](../Riemann/Deep/RCApps5/NOVELTY2.md)).
3. All orders ([AllDisc](../Riemann/Deep/AllDisc/NOTE.md)): no new analytic input. Induced ring class characters come from unique primitive pairs; genus characters and Siegel zeros are handled per order. Reviewed: one range restriction (Proposition 1°, σ ≥ 7/8, all that is used) and a wording fix, then correct. It contains the fundamental case, which passed adversarial novelty.
4. D > 0 and sectors, at 1.1307 ([RCApps4](../Riemann/Deep/RCApps4/NOTE.md)): reviewed correct in substance (units unfolded via nonnegative Fourier weights, finite truncated counts, regulator via h⁺R⁺; Theorem 3⁺ must keep H ≤ x).

**Status.** Paper proofs; not refereed; not in Lean.

**Review.** Both versions of the main note correct; significant.

**Novelty.** Nothing found states or implies it ([novelty](../Riemann/Deep/RCApps/NOVELTY.md)). For the indefinite case and the 3/2 exponent, an adversarial re-check found nothing that states or implies them ([novelty 2](../Riemann/Deep/RCApps2/NOVELTY.md)).

**Files.** [main note](../Riemann/Deep/RCApps2/NOTE.md), [first version](../Riemann/Deep/RCApps/NOTE.md), and the step notes above.
