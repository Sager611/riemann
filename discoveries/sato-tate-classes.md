# Sato–Tate in every ideal class, for almost all discriminants

[Back to DISCOVERIES](../DISCOVERIES.md)

**Setting.** f ∈ S_k(Γ₀(N)) is a newform without CM, k ≥ 2, with trivial nebentypus (e.g. f = f_E for an elliptic curve E/Q without CM). Write λ_f(p) = a_f(p) p^{−(k−1)/2} = 2 cos θ_p. J ⊂ [−1, 1] is an interval with nonempty interior.

**Claim (unconditional).** For every A > A(J), all but O(X^{1−c}) fundamental discriminants D with |D| ≤ X, negative or positive, have degree-one primes 𝔭 with cos θ_p ∈ J in every (narrow) ideal class of Q(√D), in every interval (x, x + x|D|^{−κ}] with x ≥ |D|^A.

**Threshold.** A(J) = 𝒜(M*(J)), where M*(J) is the least degree of a polynomial minorant of 1_J with positive Sato–Tate mean:

| M | 1 | 2 | 3 | 4 to 8 | ≥ 8 |
|---|---|---|---|---|---|
| 𝒜(M) | 1 | (23 − 3√5)/8 ≈ 2.0365 | 3 | 11/3 | (M + 3)/3 |
| first proved | 9/4 | 3 | 10/3 | | |

Hence A(J) ≤ max{11/3, 25.2/L + 1.7}, where L is the length of the θ-interval.

**Special cases.**
- Signs: for A > 2.0365 every class has primes with a_f(p) > 0 and primes with a_f(p) < 0. So every primitive form of discriminant D represents primes p, p′ ≤ 2|D|^A with a_f(p) > 0 > a_f(p′).
- One-sided: a_f(p) > −α p^{(k−1)/2} (α > 0) needs only A > 1, the GRH value of the method.
- Discrepancy: for A > 11/3, the smooth-weight Sato–Tate discrepancy of cos θ_p over all degree-one primes of each class in that interval is at most 10/(3A − 3).

**Before.** No bound with any power of |D| was known, per D or on average.

**Inputs.** The Pascadi–Thorner 2025 preprint (for M = 3 and M ≥ 4); with published inputs only, 11/3 becomes 5. Input (J5) is pinned to Thorner, IMRN 2025.

**New tools behind the lower thresholds ([STLemP](../Riemann/Deep/STLemP/REPORT.md)).**
- A lift of class group characters to ring class groups (Lemma P positivity over all orders: the family sieve beyond the conductor).
- A mean square for the double family: Σ_{|D|~X} Σ_{χ ∈ Ĝ(D)} |L(½ + δ + it, f ⊗ θ_χ)|² ≪ X^{3/2+O(δ)}, i.e. Lindelöf on average (X^{9/4} for Sym² f), plus the density hypothesis on average for f ⊗ θ_χ.

**Status.** Paper proofs; not refereed; not in Lean.

**Review.**
- Sign proof: flawed as first written; all five repairs checked correct.
- General note: no internal gap ([fix](../Riemann/Deep/ChebST/FIX.md)).
- STLemP: one gap (an uncovered AFE block) repaired, then correct.
- Signs and the general case form one entry, not two (reviewer, 2026-10-22).

**Novelty.**
- Nothing found states it ([signs](../Riemann/Deep/RSFam/NOVELTY.md), [general](../Riemann/Deep/ChebST/NOVELTY.md), [adversarial](../Riemann/Deep/ChebST/NOVELTY2.md)).
- STLemP: nothing states the mean square (only a per-D GRH form, Blomer–Brumley 2024), the lift or the thresholds ([novelty](../Riemann/Deep/STLemP/NOVELTY.md); arXiv API coverage partial).

**Files.** [signs](../Riemann/Deep/RSFam/NOTE.md), [general](../Riemann/Deep/ChebST/REPORT.md), [lower thresholds](../Riemann/Deep/STLemP/REPORT.md)
