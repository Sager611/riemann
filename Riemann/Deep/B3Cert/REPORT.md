THEOREM: a positive proportion of simple critical zeros in every (T, T + T^θ], θ > 371/750 = 0.494666… (< 1/2)

# B3Cert — certificate below √T (node N8)
b3cert (Orca), 20 Nov 2026, 04:23–05:25 UTC. Derivations FULLPROOF.md (§ numbers below); scripts code/, logs logs/; every run single-threaded under
.scratch/memcap.sh 290 (peaks: float 139 MB, Λ-route 272/233 MB, Arb 18.6 MB; one float run killed at 309 MB, fixed, nothing from it used).

## Theorem
For every fixed θ ∈ (371/750, 1):  liminf_{T→∞} S(T, T^θ)/N(T, T^θ) ≥ h(θ; κ) ≥ h(371/750; κ) = 3.3546·10⁻⁵ > 0,
S = simple zeros on β = 1/2 with T < γ ≤ T + T^θ, N = all zeros there (multiplicity), h(θ; k) = (3 + k − √((1 − k)(9 − k − 8c(θ))))/4,
c(θ) = 2 − θ/2 − cot(θ/√2)/√2 (Wang), κ = κ(P, Q, R = 44/5, ν = 121/1000) = 0.08503876039176576024 (Arb), (P, Q) = SL's frozen EXP-010 detector (deg Q = 201).
Moment form (E1-M′): SL Thm 1.2 (Gaussian-weighted shifted mollified second moment on windows T^θ, any fixed real Q with Q(0) = 1, error O(H/L))
holds for every θ ∈ (1/3, 1) and 0 < ν < (3θ − 1)/4: the condition θ > 1/2 is removed.
Before: E1 θ > 1/2; Wang 2609.07918 θ > 0.5501 (TeX l. 133); SL θ > 1/2; CeilAudit B3: θ > 0.49467 in floats only (reviewer GAP: uniform dual estimate unwritten).
Inputs (audited earlier unless NEW): SL (10) exact Gaussian reduction, (12) Estermann FE, (14) dual Mellin form, residues l. 444–447 (E1Audit §1, end-to-end
numerical Voronoi identity); E1 Steps 1, 2, 5 (E1Audit §2); Hankel bounds X^k∂^k a_± = O(1); Titchmarsh Lemma 4.4; SL l. 730–758 compact windows; SL Lemma 2.2
(Young's diagonal, ν < 1/2); SL Lemma 4.2, Thm 1.3 (θ-free); SL Thm 5.1 with Wang's Thm local-pair (fixed 0 < λ < θ < 1, TeX l. 231–232) and the parity product.
NEW: Lemma B (exact-profile dual bound, every θ), Lemma R (exact Mellin value of the exact profile), the §4 branch checks, the §8 Arb certificates.

## The five certificate steps — all hold
1. Step 4″ = Lemma B (§2). Keep SL's exact profile A_{−1}(ξ) = ξ^{−β}π^{−1/2}∫e^{−y²+(1−2β)y/H}e^{i(2Ty/H − λ(y)ξ)}dy, λ = 2π(e^{2y/H} − 1) (SL's exact kernel
   piece itself, K_β = Σ_{ε=±1}e(εx)A_ε exactly; no Taylor expansion, no series in j). Coordinates Λ = H(e^{2y/H} − 1), ξ = Tu/2π make the joint phase K[H log(1 + Λ/H) − Λu + σc_r√u], K = T/H, c_r = X0/K = 2√(2π)√(n/N0).
   Inner u-integral: complete the square in s = √u; m := e^{−iψ0}J has |m|, Var m ≪ X0^{−1/2} (second-derivative test). Outer Λ-integral on the stationary
   range I_σ: second-derivative test with the amplitude f·m kept inside, |ψ″| ≥ K/(8√2c_r) when c_r ≤ H/46; off I_σ two integrations by parts.
   (i) |𝒯_σ(n)| ≤ C1·H·X0^{−1/2} if c_r ≤ H/46; (ii) |𝒯(n)| ≤ C_N·T·X0^{−1/2}c_r^{−N} (E1's IBP with the pointwise |∂^kA| ≪ H^{−k}); 1 ≤ H ≤ T.
   Corollary B′: |𝒯(n)| ≪ H(pq/(nT))^{1/4}(1 + n/N0)^{−K′}T^ε and the dual total ≪ y²T^{1/2+ε}H^{−1/2}, uniformly in τ = T/H², n, p, q and complex
   shifts ≪ 1/L (Cauchy for shift derivatives). No θ condition: the profile's y-chirp (curvature 4τ) is never resolved; the inner gain X0^{−1/2} and the
   outer gain (c_r/K)^{1/2} multiply to K^{−1}. Relative error T^{1/2+2ν−3θ/2+ε}: power-small iff ν < (3θ − 1)/4.
2. Residues = Lemma R (§3): M(s) = ∫χ(x/T)A_{−1}(x)x^{s−1}dx = (H/(2√π))(T/2π)^{s−β−1}(1 + O(H/T)), via Parseval in v = e^{2y/H} − 1 (k = Tv/2π) and
   ∫ω̂(k)e^{2πik}dk = ω(1) = 1. Put in the Estermann residue formula it gives SL's displayed residue expression (built from the same leading values,
   "whole-v integral √π") times 1 + O(H/T); the pole-cancelling combination, maximum modulus, Cauchy and the 1/(gpq) sum change the main term by
   O(H²T^{−1}L^{O(1)}) = o(H/L): √πH·c unchanged. No Taylor profile and no sum over j is used (M − M_0 against SL's closed-form Gaussian is a direct comparison).
3. Compact windows (§5): SL l. 730–758 is profile-free; Gaussians of width T^{θ−η}, η = (θ − 371/750)/2, so θ − η > 371/750 ⟺ ν = 121/1000 < (3(θ − η) − 1)/4.
4. Branches (§4): off-band = SL's own (1 − χ)K_β term, unchanged for any θ < 1 (incl. the large-x tail); SL (11)'s Taylor remainder is gone; ε = +1 retained
   term ≪ (H/T)^A pointwise (y-phase derivative ≥ 2T/H), summed absolutely; K-Bessel branch on Re u = σ large using only sup|A| ≤ C: ≪ T(y²/T)^σ;
   E1 Steps 1–2 need only A ∈ C_c^∞. The chain (Emean) → Lemma 4.2 → Thm 1.3 → Thm 5.1 is θ-free on (0, 1).
5. Arb (§8; code/b3c_arb.py, logs/b3c_arb.log, 4000 bits, exact rationals, rerun identical): detector exactly admissible (P(0) = 0, P(1) = 1, Q(0) = 1,
   Q(y) + Q(1 − y) ≡ 1); check κ(3/25, 177/20) = E1's κ1 to 30 digits; κ(121/1000, 44/5) = 0.0850387603917657602445070238 ± 5·10⁻³²;
   (3·371/750 − 1)/4 = 121/1000 exactly; k*(371/750) = 0.0849914143497594430960, κ − k* = 4.7346·10⁻⁵ > 0; h(371/750; κ) = 3.354573606042305·10⁻⁵ > 0;
   h(0.4947) = 7.0693·10⁻⁵, h(1/2) = 0.0059310; detector threshold θ_det = 0.494636567599178 ± 7·10⁻¹⁷ < 371/750 (the moment edge binds).

## Float falsification (§6)
- F1 (code/b3c_fals.py): exact profile by FFT (vs brute y-quadrature 2·10⁻⁵), both Hankel branches, r = n/N0 ∈ [0.003, 8] (29 points; 𝒯 depends on n, p, q only
  via n/pq, so all (n, p, q) are covered), 4 shift sets incl. (0.2 − 0.1i, 0.03 + 0.2i), T ∈ {10⁶, 10⁷, 10⁸}, θ ∈ {0.40, 0.45, 0.47, 0.49, 0.495} (τ = 1.15–39.8):
  max |𝒯|/(HX0^{−1/2}) = 0.705 over all 60 rows, flat in τ and T; max |𝒯_±| = 0.3535 = 1/(2√2), the 2D stationary-phase constant of Lemma B's coordinates;
  tail r ≥ 4: ≤ 5.4·10⁻⁴; Lemma R: |M/M_pred − 1| ≤ 1.1·10⁻⁷ (bound H/T ≥ 1.6·10⁻⁵). No case approaches Lemma B's constant or breaks its scaling.
- F2 (code/b3c_lambda.py): Lemma B's own coordinates with a genuine smooth cutoff reproduce F1 to ≤ 8·10⁻⁵ (T = 10⁶, θ = 0.45, 0.40, r = 0.01, 0.3, 2);
  sup|m|X0^{1/2} = 5.7 ≤ M0 = 54.4, Var m·X0^{1/2} ≤ 15 ≤ M1 ≈ 2.8·10⁴, |J| on N_σ ≥ 10⁴× below its bound; explicit C1 ≈ 8.5·10⁴ (observed ≤ 0.36 per branch).
Smallest failing case of the key lemma: none found in the tested range; beyond its hypothesis (c_r > H/46, i.e. n/N0 > H²/53180) part (ii) bounds the tail.

## Review
reviewer, statement and proof pasted inline (FULLPROOF §7): pass 1 GAP (bookkeeping: Var f majorant at H = 1, (ii) needs H ≤ T, tail exponent) — "stationary-range
cancellation survives … not a failure of the iterated cancellation argument"; pass 2: "No further gap found in Lemma B(i)–(ii)", Cor B′ used N0 ≤ T — repaired with N0 ≤ T²;
pass 3: summation valid, N ≥ 2(A + 5)/ε — fixed; pass 4 (Lemma R + §4): CORRECT. All repairs are in FULLPROOF.
Orchestrator referee (20 Nov): Lemma B, Cor B′, Lemma R CORRECT; removable GAP — the claimed identity Σ_j A_{−1,j} = A_{−1} (false as a convergent series) and
the residue-sum identification built on it; both deleted, nothing else used them (FULLPROOF §7 Repair). Theorem unchanged.
Novelty: CeilAudit's scan this morning (lst.py, math.NT/CA/CV 2026-10, -11; Zenodo "simple zeros" "short intervals": only SL 23134787/23135354) rechecked:
nothing on simple zeros in windows below √T; 2610.26841 (isolated zeros, 13/20) works over long ranges with a simple-zero input; Wang's simple-zero range is θ > 0.5501.
Limits: the frozen detector; the one-piece optimum (CeilAudit b3_ceilnu, θ_c* ≈ 0.4940) and per-θ re-optimised ν are not certified; liminf only (constants depend on θ via η, D).

RH relevance: N8 (critical-line proportion): ends our stated barrier "1/2 is the method's barrier on both sides" (E1.md §5, E1Audit l. 81): the moment side now reaches
θ > 1/3 + 4ν/3; the binding limit is the detector.

## Files
Riemann/Deep/B3Cert/REPORT.md (this), FULLPROOF.md (§1 setup, §2 Lemma B + Cor B′, §3 Lemma R, §4 branches, §5 windows + chain, §6 floats, §7 reviews, §8 Arb);
code/b3c_fals.py, code/b3c_lambda.py, code/b3c_arb.py; logs/b3c_fals_part1.log, b3c_fals_part2.log, b3c_lambda.log (+ partial), b3c_arb.log, *.err (memcap/time peaks).
