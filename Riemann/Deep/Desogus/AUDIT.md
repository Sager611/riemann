# Independent audit: Desogus 2609.20367v3, Thm 8.6 ineq. (201)
**VERDICT: CONFIRMED.** (201) is false at k = 7 by about 10 orders of magnitude. The bare diagonal alone already breaks it: Q(e_7) = 1.2263 < g ≥ 1.878.

## (1) Weil form: no mismatch found
- The code uses the standard Guinand–Weil form, c(0)(−γ − log π − log(1−e^{−4a})) + ∫₀^{2a} 2(c(0)e^{−2t} − c(t)e^{−t/2})/(1−e^{−2t}) − 2(∫v e^{y/2})² − 2Σ_{n<e^{2a}} Λ(n)n^{−1/2}c(log n). Its spectral side Σ_all γ |v̂(γ)|² is the paper's (3).
- Checked against zeta zeros: for v = x(1−x²)^6 and x(1−x²)^10 on window log 8, the formula and 2Σ_{40 zeros}|v̂|² agree to 12 digits (3.46149513652e-7; 2.70254218392e-6).
- Checked against the paper's (63) by hand: h(12) + V(13) − K_γ(14) reproduces the scalar −log a − log 2π − γ = c₀(a) (15). Pole −2⟨sinh(y/2), v⟩² matches (11). Primes n < e^{2a} with shift r_n = log n/a match.
- Coordinates match: H_a = L²_odd(−a, a) is physical (§1.1), a_k = ½ log k (§6.13), so A_8 is on window log 8 with primes n ≤ 7. Collar K_k = odd L²(a_k < |y| < a_{k+1}) (64).
- One real ambiguity: the paper's e_k ("normalized true ground line in the new collar", q_Y = t^{−1/2}, Lemmas 6.4/6.8) is not pinned down. The code takes sgn(y)e^{−|y|/2}; the other reading, sinh(y/2) on the collar, gives ratios within 6% of it. No normalisation constant can close a 10¹⁰ gap.

## (2) Re-run of `cert/rigorous.py 80 7` (Arb)
- W(v) ∈ [9.7789091e-29 ± 2.6e-37], ⟨e_7, v⟩ ∈ [9.4217747e-10 ± 4.6e-18]. **Ratio ∈ [1.1016027e-10 ± 2.3e-18]**, identical to rigorous.log. Cost: 84 s, 39 MB.

## (3) New mpmath script `.scratch/desogus-audit/audit_mp.py` (exact-rational correlations, odd monomial basis, 50 digits)
- Q(e_7) = **1.226339147** (paper: B_7 = 2.513); Q(e_12) = 1.731997 (paper: B_12 = 3.025). min Q(v)/⟨e,v⟩² at k = 7, window log 8:

| degree | e = sgn·e^{−\|y\|/2} | e = sinh | λ_min of Galerkin matrix |
|---|---|---|---|
| 9 | 0.0528 | 0.0544 | 4.1e-10 |
| 13 | 0.00928 | 0.00968 | 3.5e-14 |
| 19 | **1.577e-4** | 1.671e-4 | 5.4e-21 (positive) |

- Direct quadrature without the matrix gives the same deg-19 value (1.576585889e-4). The ratio keeps falling with degree (deg 39: 5.6e-7; deg 79: 1.1e-10), as expected when 1/⟨e, A_8^{−1}e⟩ is driven by A_8's near-null modes.

## (4) Reduction to "Q(v) ≥ g⟨e_k,v⟩² for all odd v": valid under the paper's definitions
- H_{a_{k+1}} = H_{a_k} ⊕ K_k. Lemma 6.27 gives H_k* A_{k+1} H_k = T_k, the Schur complement, which is the infimum over the old part (Lemma 6.30).
- So Q(x ⊕ f) ≥ ⟨H_k f, A H_k f⟩ ≥ g|⟨e_k, f⟩|² + R⊥ ≥ g|⟨e_k, v⟩|², since R⊥ ≥ 0 and e_k ∈ K_k. The paper claims g_k^{M+} ≥ 1.87814156 ((210); table line 3990). Any g > 1.2e-10 is already refuted.
- A_7 ≻ 0 holds numerically; if it failed, the base certificate 6.64 would fail instead. A polynomial v with jumps at ±a lies in the form domain.
- Caveats: only run (2) is rigorous; (3) is floating point. Both readings of e_k were tested. Refuted: the inequality (201). Not refuted: the conclusion T_k ≻ 0.
