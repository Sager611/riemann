NO REFUTATION: Hakimoglu-Brown's Λ ≤ 0.13963 (arXiv 2610.12803v1) is sound under "no zeta zero with Re s > 7/8" plus published inputs. Every step was checked, the new pair inequality has a proof, and every finite certificate was re-verified exactly. Our 0.148409 and 0.147374 are not records.

# Audit of Hakimoglu-Brown, "A conditional bound Λ ≤ 0.13963" (2026-10-22)
Read from the TeX (e-print), with Polymath's TeX/PDF (arXiv 1904.12438v2) and Polymath's barrier code (`TloopSinglematv2.c`, at the commit the paper cites). Details: NOTE.md. Scripts and outputs: evidence/.

**Chain, step by step (all hold):**
1. *Hypothesis → start.* Re ρ ≤ 7/8 and ξ(s) = ξ(1−s) give |Im z| ≤ 3/4 for the zeros of H₀, so M(0) ≤ 9/16. This is the only place the hypothesis is used.
2. *Wall.* Polymath §8.4, claim (a), says exactly what is used: H_t ≠ 0 for X ≤ x ≤ X+1, 1/5 ≤ y ≤ 1, 0 ≤ t ≤ 1/5, with X = 6·10¹⁰+83951.5. Platt's RH check to height 3.06·10¹⁰ covers |x| ≤ X+1 at t = 0.
3. *Lemma 2 (finite region).* A first-contact argument shows that for t ≤ 1/5, zeros with |x| ≤ X have |y| < 1/5. The zero's own conjugate gives −1/h. Pairs with v ≤ h give K_h ≥ 0. Higher pairs lie at |x| > X+1, so d² ≥ 1 ≥ v² and K_h > 0. Multiple zeros split vertically backward in time (Polymath Prop 3.1(ii)).
4. *Lemma 4 (approximation on 3 ≤ y ≤ 6).* Inputs: Polymath (39)–(41), Prop 6.1 (any real σ) and Prop 6.2 (any real σ and K; for σ < 0 the remainder needs K+σ ≥ 2, which K = max(100, ⌊−U⌋+3) meets). The remainder rederived for σ ∈ [−3,4] is correct. All numerical side claims re-checked: d(x) = 3.1e-11; Gaussian checks 5.7e-10, 0.2548, 0.0297; ε_σ < 1/25; |γ| = Q^{−y/2}(1+O(x⁻²)). Claimed error 1e-6; true error about 1e-10.
5. *Prop 11 (S and ∂_pS at height 61/20).* F = ζ(s_B)(1+ε), Leibniz bounds for ε, Euler-product bounds for ζ'/ζ, (ζ'/ζ)' and 1/ζ. The signs check (∂_pS = −Re(H'/H)'). Independent check at t = 0, where S = ½ Re ξ'/ξ(w) holds exactly: at four points x ≥ X, S = 5.33–5.56 ≥ C₀ = 5.2902 and ∂_pS = 0.02–0.16 ≤ G₀ = 0.234.
6. *Lemma 10 (new conjugate-pair inequality): PROVED.* K_p = 2pw/Δ and ∂_pK_p = ((2w+4B)Δ−8Bw²)/Δ². K_h/h ≥ 2/(D+A+3r) is equivalent to r(A−r) ≥ 0. Then Φ ≥ 0 closes the case 𝒩 > 0. At (λ,μ) = (1,0) it reduces to the old p² ≥ 5q comparison.
7. *Theorem 12 (24 steps).* q' = −2−4hΣ'K_h ≤ −2−k_i·q. The own pair gives K_p(0,h) = 2p/(B−q) and ∂_pK_p(0,h) = −2(B+q)/(B−q)². The bracket decreases in q. The strip theorem keeps every contact height at q ≤ q_i, which is the certificate's range. Each gate is exactly Q_i(t_{i+1}) ≤ q_{i+1}, so each step's hypothesis holds at the previous step's output.

**Certificates (python-flint/Arb at 128 bits plus exact fractions; peak memory 33 MB):**
- All 24×9 entries of the final tables (Z, ψ, ψ₂, δ₀–δ₂, ρ₀–ρ₂) are true upper bounds (relative slack ≥ 2e-14).
- From the printed entries: max ε₀ = 0.01909 < 1/50. min(C^cert−C) = 2.0511e-6 (row 2). min(G−G^cert) = 1.1046e-6 (row 6). Smallest gate (E₂₀) = 2.787e-8 (row 5). 8.6377 ≤ k ≤ 10.3152. All exactly as printed.
- Bernstein: all 4800 coefficients exceed 13/10000. The least one is the printed fraction exactly (row 21, piece ℓ = 4).
- t₂₄ = 11963/100000, so Λ ≤ 13963/100000.
- Earlier tables (0.17172, 0.16408): every entry is an upper bound. Margins 8.56e-7 and 1.487e-6; gates 3.48e-7 and 1.25e-7.
- Scalars all pass: 269/200; 39/100 (0.387405); 5.25228; L/4 = 5.5716469 > 5.5716; √(X/4π) = 69098.88; J(2) ≤ 0.94287; tails; thresholds 0.90664, 0.93835, 0.88328, 0.91043; 0.1797275414 and 0.1938798367.

**The paper's footnote on Polymath (91)** is right: the time term has no 1/(2n). One extra observation: the code puts ⌈D_z⌉ points on each edge, so the mesh term D_z/(2(num−1)) is slightly above ½. With the step (min|f|−½)/D_t, the forward bound therefore does not quite close at the end of each step. The next node's backward bound (at least 0.46 − 0.04) covers that gap, so this is harmless. Only the very end of the last step, near t = 0.2, could stay uncovered, and the paper uses t ≤ 0.1197.

**Remaining premises:** the 7/8 statement (a theorem in our Lean port), plus published results: Polymath's barrier, its Theorem 1.5(i) and Props 3.1/6.1/6.2, and Platt's RH height. Polymath's barrier computation was not replayed.

**Consequence:** the removal of our de Bruijn–Newman entry from DISCOVERIES stands.
