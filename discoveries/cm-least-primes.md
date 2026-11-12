# Least primes in every ideal class of almost all CM fields

[Back to DISCOVERIES](../DISCOVERIES.md)

**Setting.** F ≠ Q is a fixed totally real field of degree n; K/F runs over CM extensions.

**Claims (unconditional).**
1. Large sieve: a content-weighted large sieve over the class group characters of K with constant ≈ N + √d_K (via a relative lattice lemma; the regulator factor of F is unnecessary).
2. Zero density: Σ_χ N(σ; t, χ) ≪ (d_K Tⁿ)^{2(1−σ)/(3−2σ)+ε} per unit window, for 0.6 ≤ σ ≤ 1 (for σ ≤ 1/2 it follows from the trivial count, CMQuartic Theorem A(b), as 2(1−σ)/(3−2σ) ≥ 1/2 there; 1/2 < σ < 0.6 is not proved).
3. Least primes: for every A > 1, almost all CM fields K/F with prime (or odd δ₀-rough) relative discriminant have degree-one prime ideals in every ideal class, in every interval (x, x + x d_K^{−κ}] with x ≥ d_K^A (some κ = κ(A) > 0, half the κ of CMQuartic NOTE l.219–220; X₀ effective but astronomical).

**Sharpness.** A > 1 cannot be improved: the principal class has no degree-one prime of norm below d_K/(4ⁿ d_F²) (the key lemma of Miller–Peckner 2012).

**Before.** No almost-all result was known for F ≠ Q. For every field: d_K^{694} (Thorner–Zaman 2017); under GRH d_K^{1+ε}.

**Steps.** First A > 3/2, then (3+√69)/10 ≈ 1.1307, then every A > 1.
- 1.1307 ([CM2](../Riemann/Deep/CM2/NOTE.md)): a family large sieve over pairs (K, χ ∈ Ĉl_K) at family size X^{3/2}. Reviewed correct in substance.
- Every A > 1 ([CM3](../Riemann/Deep/CM3/NOTE.md)): count the pairs over all elements of each order O_F + 𝔣O_K (a relative Lemma P). This removes the N³/X term and gives Σ_K N_K(σ) ≪ X^{2−σ+ε}, the density hypothesis on average over K. Reviewed correct for detector sequences; 5 independent checks OK.

**Method.** Lemma 1 is classical and the method is routine (the same mechanism as our ring class theorems); the statements are new.

**Status.** Paper proofs, refereed by cmaudit (12 Nov 2026): AUDIT PASSED — CMQuartic Lemmas 1–4, Theorem A (0.6 ≤ σ ≤ 1), Proposition C, Theorem B, CM2 Theorems 1–3, CM3 Lemma P_F and Theorems 1–3 correct as stated; defects cosmetic (σ-range and wording, fixed above); no counterexample in 22 test fields with d_K ≤ 10⁸ or in the family sieve over all 1429 CM fields over Q(√5) with 500 < N𝔡 ≤ 16000; orchestrator rerun of its genus check identical ([CMAudit](../Riemann/Deep/CMAudit/REPORT.md)). Not in Lean.

**Review.** Main inequalities correct.

**Novelty.** Nothing states or implies it ([novelty](../Riemann/Deep/CMQuartic/NOVELTY.md)). For every A > 1, a third adversarial check found nothing either ([novelty 3](../Riemann/Deep/CM3/NOVELTY2.md)).

**Files.** [main note](../Riemann/Deep/CMQuartic/NOTE.md), [CM2](../Riemann/Deep/CM2/NOTE.md), [CM3](../Riemann/Deep/CM3/NOTE.md)
