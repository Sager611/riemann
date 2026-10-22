DISCOVERY: for every A > 1, all but O(X^{1−c}) non-square discriminants Δ ≡ 0, 1 (mod 4) with |Δ| ≤ X (fundamental or not, negative or positive) have degree-one prime ideals in every (narrow, if Δ > 0) ring class of the order of discriminant Δ, with norm in every interval (x, x + x|Δ|^{−κ}] with x ≥ |Δ|^A; so every primitive binary quadratic form of discriminant Δ, definite or indefinite, represents a prime ≤ 2|Δ|^A (c, κ > 0 depend on A; for 1 < A ≤ 1.2, κ = (A−1)/32000, c = (A−1)/16000).

# AllDisc report (2026-10-22)
Proofs: [NOTE.md](NOTE.md). Checks: [checks/](checks/). Brief: `.scratch/orca/alldisc-brief.md`. Paper proofs, not refereed, no Lean. Unconditional; same inputs as [RCApps5](../RCApps5/NOTE.md).

## Result (task 1)
- Theorems 3^o (Δ < 0) and 3⁺^o (Δ > 0) extend RCApps5 from fundamental D (61% of all discriminants) to every quadratic order. Each class gets ≥ ν_Δ H/(16h log 2x) primes.
- Corollary: for all but O(N^{1−c}) integers n ≤ N (any n, not only squarefree n ≡ 1, 2 mod 4), x² + ny² has a prime in every (x, x + x(4n)^{−κ}] with x ≥ (4n)^A. The same holds for x² − ny² (non-square n).

## What had to be redone (tasks 2–5)
- Family sieve (Thm 1^o): no new work. Lemma P already sums over all orders; the left side just keeps the non-fundamental Δ.
- Detector and density (Lemma D, Thm 2^o): λ_χ(p) = 0 = (Δ/p) for p | f, L(s,χ) = L(s,χ*)·(zero-free factor). Result: Σ_{Δ} N_Δ(σ) ≪ X^{2−σ+ε} over all orders.
- Induced characters (Lemma M): ring class fields satisfy L_{f₁} ∩ L_{f₂} = L_{(f₁,f₂)}, so every χ comes from a unique primitive pair (Δ′, χ′). A primitive pair lies under ≤ 1 + (X/|Δ′|)^{1/2} discriminants Δ = m²Δ′, and each GL₂ form comes from ≤ 6 primitive pairs. Proposition 1^o: Σ_Δ N^c_Δ ≪ X^{10(1−σ)/(3−2σ)+ε} + X^{1/2+ε}; the X^{1/2} (small conductors, large multiplicity) is harmless.
- Genus and Siegel zeros (Lemma G): real characters of G_Δ correspond to disjoint pairs {d, d′} with L(χ*) = L(χ_d)L(χ_{d′}) and dd′ | Δ. An inert-prime argument shows the field is biquadratic for both signs. E₂ counts the Δ with d | Δ; ν_Δ ≥ c|Δ|^{−δ/2}/log|Δ| as before.
- Δ > 0: units and regulators of orders enter only through Lemma P⁺ (per order), h⁺R⁺ = √Δ L_Δ and h⁺ ≤ 5Δ^{1/2}log²(3Δ).

## Where it could fail (task 6)
- Nowhere. Each use of fundamentality in RCApps–RCApps5 has a replacement (NOTE §8).
- Not covered: sectors over non-maximal orders (not written out). A < 1 remains NOGO: the principal class has no prime below |Δ|/4.

## Checks (all under hard alarms, < 40 MB)
- Lemma G: 1462 discriminants with |Δ| ≤ 1500 (550 non-fundamental), 0 failures. Wide-class mutation: 465 of 712 fail.
- Family sum over all orders: 0.22·X^{3/2}‖b‖², stable for N from X/16 to 1.5X (X = 2·10⁴, 4·10⁴); about one third of it comes from non-fundamental orders.
- Least prime in every class, 19857 discriminants with |Δ| ≤ 20000 (7699 non-fundamental): P(Δ) ≤ 4.17|Δ| throughout; non-fundamental at most 1.34|Δ|.

## Novelty
- arXiv math.NT 2026-10 (625 titles) and 2026-09 (812) scanned; abstracts of 2610.28001, 2610.20214 and others read. Nothing on primes in classes for non-fundamental discriminants or on averages over orders.
- Baier–Zhao average over squarefree k only. Since the statement contains RCApps5's (already new per [NOVELTY2](../RCApps5/NOVELTY2.md)), it is new.
- No new analytic input: the content is that the bookkeeping over orders closes.
