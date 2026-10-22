NO REFUTATION: Deng's proof (arXiv 2609.22336v2) holds up. Every lemma was re-derived by hand and tested exactly (PARI/GP, p-adic factorizations): no false step, no gap. Granting the cited classical inputs, Stieltjes' conjecture (P_{2j} and P_{2j+1}/x irreducible over Q) is proved.

# Audit of Deng, "A proof of the irreducibility conjecture for Legendre polynomials" (2026-10-21)
Read from the v2 TeX e-print (27 Sep 2026). v1 (16 Sep) has the same skeleton, a weaker large-prime bound and leaned on Groth's computation up to n < 10⁵. Details: [NOTE.md](NOTE.md).

**Chain (every step holds).** Suppose M_n = 2^{v₂(n!)}P_n/x^ε = A(x²)B(x²) in Z[x], k = deg A ≤ e = deg B, k+e = m = ⌊n/2⌋, R = Res_t(A,B) ≠ 0, 𝓡 = odd part of |R|.
1. *2-adic.* Wahab's polygon: one even Q₂-irreducible factor per binary digit ≥ 2 of n, so a rational factor is a union of digit groups. Hence A(1) = 2^{2k−s₂(k)} < a = lc A, and v₂(t_i−u_j) is the min of two unequal values, so v₂(R) ≤ e(2k−s₂(k)).
2. *Odd p.* With δ_p(α,γ) = v_p(α−γ) − min(v_pα,0) − min(v_pγ,0) ≥ 0: 2v_p(R) = Σ_{α∈F,γ∈G} δ_p, and a full row sum equals v_p(P_n'(α)) (reciprocal version if v_pα < 0). The identity (1−α²)P_n'(α)U_n(α) = 1 (U_n has denominators ≤ n), plus orthogonality (denominators ≤ 2n−1) when v_pα < 0, caps each row at ⌊log_p(2n−1)⌋. So v_p(R) ≤ k⌊log_p(2n−1)⌋.
3. *Odd p > m.* P_n ≡ x^{λp}P_r (mod p), P_r squarefree, so unit roots have v_p(P_n') = 0. v_p(lc M_n), v_p(M_n(0)) ≤ 1 put all roots of positive (negative) valuation in one factor. Every cross δ_p = 0, so v_p(R) = 0. This step is essential: without it the upper bound is ≈ e^{4km}.
4. *Upper.* 𝓡 | D_n^k, D_n = ∏_{odd p ≤ m} p^{⌊log_p(2n−1)⌋} divides lcm(1..m)·lcm(1..⌊√(2n−1)⌋) < 3^{m+√(2n−1)} (Hanson).
5. *Lower.* |R| = a^{e−k}∏|H'(t_i)|/V_A. (1−x²)P_n'² + n(n+1)P_n² increases on [0,1), so |H'(t_i)| > 2^{v₂(n!)}γ_n with γ_n > √m. Chebyshev + Hadamard give V_A ≤ k^k/2^{2(k−1)²} for ANY k points of [0,1]: no spreading or interlacing is assumed. With step 1: log𝓡/k > m log 4 − log(4(m+1)²/√m).
6. *Close.* 4 > 3 gives a contradiction for n ≥ 256. n ≤ 121: Cullinan–Hajir Thm 1.9(b). 122–255: their Cor 3.4(b), or the paper's 15-row table.

**Exact checks (PARI/GP, peak 81 MB, all pass).**
- (3.4) and (3.6) as polynomial identities, n ≤ 40.
- Lemma 3.3 per root, over every Q_p-irreducible factor, all odd p ≤ 2n−1, n ≤ 70: 6023 factors, 0 violations, 776 equalities (the cap is sharp). Lemma 3.6 for n ≤ 90, all p ∈ (m, 2n): 0 violations. Lemma 3.7 for n ∈ [64, 4000].
- Lemma 2.1 and Wahab's polygon (vertices and hull), n ≤ 300. Q₂-factors: even, one per digit, v₂(f(1)) = deg f − 1, and v₂(R) ≤ e(2k−s₂(k)) for all 387 splits, n ≤ 110.
- Lemma 5.4: min |H'(t)|/(κ_nγ_n) = 1.16 over n ≤ 300 (at 400 digits; 80 digits gives a spurious 0.04 from cancellation). Lemma 5.5 against the Fekete maxima, k ≤ 150: holds.
- (4.4) and D_n < 3^{m+q} for n ∈ [64, 5000]; (6.1) for m ∈ [128, 2·10⁵]; the set 𝒮 and all 15 table rows reproduced.
- Hypothetical factor degrees: for every admissible k (digit subsets, k ≤ m/2) and 64 ≤ n ≤ 1200, the sharper bound (5.10) beats the exact log D_n (least margin 5.0, at n = 86, k = 8). The method alone covers n ≥ 64.
- Citations say what is used (Cullinan–Hajir, Indag. Math. 25 (2014)). Independently, PARI's polisirreducible confirms P_n and P_n/x irreducible for n ≤ 255.
- The scouts' "duplicate lemma header" is in neither v1 nor v2.

**Premises taken as given:** Wahab's polygon for n > 300 (Duke 1952, re-derived by Cullinan–Hajir), Hanson's lcm(1..N) < 3^N, Cullinan–Hajir's Magma Galois computation (re-confirmed here for n ≤ 255).
**Status:** not our result, so no DISCOVERIES or REFUTATIONS entry. Unrefereed (math.GM). A proof of an 1890 conjecture needs expert review; a Lean port is feasible but heavy (needs p-adic Newton polygons and Wahab).
