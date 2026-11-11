# The limsup of Σ(−a)^{Ω(n)} tends to 1: Johnston–Leong–Tudzi's conjecture

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** For real a > 2 let W_a(x) = Σ_{n≤x} (−a)^{Ω(n)}, α = log₂ a and s_a = limsup_{x→∞} |W_a(x)|/x^α. Then lim_{a→∞} s_a = 1. More precisely, for every a > 2,
(a/(1+a))·(1 − Σ_{m odd > 1} a^{−j(m)}) ≤ s_a ≤ (a/(1+a))·Π_{p odd}(1 − a^{1−log₂ p})^{−1}, with j(m) = ⌈log₂ m⌉ − Ω(m) ≥ 1,
and both sides tend to 1.

**Before.** D. R. Johnston, N. Leong, S. Tudzi, arXiv 2408.04143v3 (Sep 2024; Experimental Mathematics 2025), Theorem limsupthm2 (TeX l. 203–212): s_a exists and 1/2 ≤ s_a < ((a−1)/(a+1))·((3^α − 1)/(3^α − a))·ζ(α); the final conjecture (l. 1142–1147) is lim_{a→∞} s_a = 1, from local maxima 0.9981 (a = 1000) and 0.9991 (a = 2000) on 10⁴ ≤ x ≤ 10⁶.

**Method.** Write n = 2^k m with m odd and g(m) = (−a)^{Ω(m)}. Summing the geometric series in k gives the exact identity
W_a(x) = (1/(1+a)) Σ_{m odd ≤ x} g(m) + (a/(1+a)) x^α Σ_{m odd ≤ x} g(m) m^{−α} h(log₂(x/m)), h(w) = (−1)^{⌊w⌋} a^{−{w}}, |h| ≤ 1.
Σ_{m odd} |g(m)| m^{−α} = Π_{p odd}(1 − a p^{−α})^{−1} < ∞ for a > 2, so the first sum is o(x^α) by dominated convergence (Σ_{m≤x}|g(m)| = x^α Σ_{m≤x}|g(m)|m^{−α}(m/x)^α), and the upper bound follows from |h| ≤ 1. At x = 2^k the m-term is (−1)^k (−1/a)^{j(m)} (log₂ m ∉ ℤ for odd m > 1, and Ω(m) ≤ log₃ m < log₂ m gives j(m) ≥ 1), hence |W_a(2^k)|/2^{kα} → (a/(1+a))|1 + Σ_{m odd>1}(−1/a)^{j(m)}| ≥ the lower bound. Σ_{m odd>1} a^{−j(m)} ≤ Σ a^{Ω(m)−log₂ m} → 0 as a → ∞ by dominated convergence (majorant 3^{Ω(m)} m^{−log₂ 3}, the a = 3 product), and the upper bound tends to 1 the same way.

**Certificate.** Exhaustive segmented sieve of W_a(n) for n ≤ 2^30 and 17 values 2 ≤ a ≤ 32 (sunWa.c, long double; a = 2 exact to 2^32 in sunW.c): octave suprema stable to 4–5 digits over the last 6 octaves. The expansion's limit value (a/(1+a))|1 + Σ_{m odd>1}(−1/a)^{j(m)}| (m ≤ 2·10⁶, orchestrator's independent script) equals the sieve's supremum on [2^29, 2^30) to 5 digits: 0.94391 (a = 32), 0.89918 (16), 0.87476 (12), 0.83710 (8), 0.81232 (6), 0.79457 vs 0.79456 (4). Lower bounds: 0.93440 (a = 32), 0.86177 (a = 16).

**Review.** Agent's verdict-only reviewer: correct, elementary, resolves the conjecture. Orchestrator's reviewer (proof pasted inline): CORRECT, SIGNIFICANT ("elementary does not itself mean routine: the cited upper bound alone cannot supply this lower bound"); it replaced the pointwise justification of the o(x^α) step by the dominated-convergence argument written above.

**Novelty.** arXiv: title search "summatory function" (only 2408.04143 itself on (−a)^{Ω(n)}), monthly listings 2026-10/11 (lst.py); OpenAlex: one citing work (2025, unrelated); Zenodo quoted queries ("summatory function" AND "Omega(n)", Johnston AND Leong AND Tudzi): 0; HAL: only Azzouza's Σ′2^{Ω(n)}/n (different sum); our files: nothing. The authors' own method (after Tenenbaum's Exercise 58) gives only s_a ≥ 1/2.

**Caveat.** Off the RH tree: the main term comes from the Euler factor at 2, not from zeros. Lean 4 + Mathlib: OmegaLean.tendsto_s : Tendsto s atTop (𝓝 1) with the literal limsup, plus the explicit bounds L(a) ≤ s_a ≤ U(a) for a ≥ 64, sorry-free, standard axioms only ([LEAN.md](../Riemann/Deep/ConjFals/LEAN.md); omegalean, 11 Nov 2026; a simpler finite-sum route than the paper proof above).

**RH relevance.** None directly (side result); context N10 (Pólya/Mertens-type size conjectures for Σ λ-like sums).

**Files.** [ConjFals report §2](../Riemann/Deep/ConjFals/REPORT.md), sunW.c, sunWa.c, sunWa.log, sunW_2e32.log in Riemann/Deep/ConjFals/.
