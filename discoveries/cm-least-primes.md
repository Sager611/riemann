# Least primes in every ideal class of almost all CM fields

[Back to DISCOVERIES](../DISCOVERIES.md)

**Setting.** F ≠ Q is a fixed totally real field of degree n; K/F runs over CM extensions.

**Claims (unconditional).**
1. Large sieve: a content-weighted large sieve over the class group characters of K with constant ≈ N + √d_K (via a relative lattice lemma; the regulator factor of F is unnecessary).
2. Zero density: Σ_χ N(σ; t, χ) ≪ (d_K Tⁿ)^{2(1−σ)/(3−2σ)+ε} per unit window.
3. Least primes: for every A > 1, almost all CM fields K/F with prime (or rough) relative discriminant have degree-one prime ideals in every ideal class, in every interval (x, x + x d_K^{−κ}] with x ≥ d_K^A.

**Sharpness.** A > 1 cannot be improved: the principal class has no degree-one prime of norm below d_K/(4ⁿ d_F²) (the key lemma of Miller–Peckner 2012).

**Before.** No almost-all result was known for F ≠ Q. For every field: d_K^{694} (Thorner–Zaman 2017); under GRH d_K^{1+ε}.

**Steps.** First A > 3/2, then (3+√69)/10 ≈ 1.1307, then every A > 1.
- 1.1307 ([CM2](../Riemann/Deep/CM2/NOTE.md)): a family large sieve over pairs (K, χ ∈ Ĉl_K) at family size X^{3/2}. Reviewed correct in substance.
- Every A > 1 ([CM3](../Riemann/Deep/CM3/NOTE.md)): count the pairs over all elements of each order O_F + 𝔣O_K (a relative Lemma P). This removes the N³/X term and gives Σ_K N_K(σ) ≪ X^{2−σ+ε}, the density hypothesis on average over K. Reviewed correct for detector sequences; 5 independent checks OK.

**Method.** Lemma 1 is classical and the method is routine (the same mechanism as our ring class theorems); the statements are new.

**Status.** Paper proofs; not refereed; not in Lean.

**Review.** Main inequalities correct.

**Novelty.** Nothing states or implies it ([novelty](../Riemann/Deep/CMQuartic/NOVELTY.md)). For every A > 1, a third adversarial check found nothing either ([novelty 3](../Riemann/Deep/CM3/NOVELTY2.md)).

**Files.** [main note](../Riemann/Deep/CMQuartic/NOTE.md), [CM2](../Riemann/Deep/CM2/NOTE.md), [CM3](../Riemann/Deep/CM3/NOTE.md)
