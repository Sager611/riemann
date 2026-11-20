AUDIT PASSED: the claim (2p − 1 ≥ 0.344982508 of the zeros of b_1e^{iθ_1}L(s,χ_1) + b_2e^{iθ_2}L(s,χ_2) on the line, odd multiplicity) holds; one citation gap (DERIVATION §8.2 uses Lamzouri's unproved "dyadic variant" of Garunkštis–Paliulionytė's pair-correlation lemma) is closed by a full-range route that needs only their published statement (FULLPROOF §G).

Auditor: lcaudit (Orca), 2026-11-20, 08:32–09:08 UTC. Target: discoveries/linear-combinations-line.md; Riemann/Deep/IsoZero/ REPORT.md "Side result", DERIVATION.md §8, cert8.py. Chain with source lines: [FULLPROOF.md](FULLPROOF.md).
Sources (latest versions, TeX): Wu arXiv 1802.09704v2 (Math. Z. 293 (2019)); Lamzouri arXiv 2610.26841v2; Garunkštis–Paliulionytė (GaPa), Glasnik Mat. 60(80) (2025) 207–227, accepted manuscript (journal site, pGM7479.pdf).

## Short answer
| item | result |
|---|---|
| (1) per-character transfer | every constant of p(c) is χ-free: C_MT (pair correlation, same constant), λ_C, θ̄, ρ_max, α_ρ (kernel and c only), B(ρ) (Wu's main term at the shift has the ζ leading term) ✓ |
| Wu Prop. 3.1 | statement authenticated (l.388–413: complex a, b ≪ 1 uniformly, w ∈ [T,2T], (B1) θ < 4/7); our mollifier is Wu's own Conrey part (l.265–266, P_2 = 0) ✓ |
| Wu's two typos | both authenticated in TeX (exponent (h,k)^{α+β} at l.400/202/422 vs 1/[h,k] at l.405 and l.1353; +R/L_χ at l.171 vs a = b = −R at l.203/426 and left-shift mollifier l.262) ✓; three further lower-order slips found, harmless (FULLPROOF D3(iii)) |
| main term at the shift, term by term | sympy, general quartic P: both terms, q/φ(q) cancellation, c(−R,−R) = c_1, c(R,R) = e^{−2R}c_1, FE symmetry, R → 0 limit 1 + 1/θ, inf_P = c* ✓; exact sums to y = 10⁶ converge to it at rate 1/log y, M_5/M_1 → 1 ✓ |
| (2) Pro:BombieriHejhal use, order of limits, p ∣ q | c-isolated ⇒ not in H_j for log log T > c⁻⁴ (boundary zeros O(1)); ε fixed, T → ∞, then ε → 0; Euler factors at p ∣ q cancel at leading order, B(ρ) unchanged ✓ |
| pair correlation per character | GAP in the citation chain: GaPa's Lemma 3.1 is for 0 < γ, γ′ ≤ T; Lamzouri states a dyadic form without proof (l.967), and §8.2 uses it. Repaired (FULLPROOF G2): full-range isolation + block-wise hinge + BH summed over dyadic blocks; same 2p − 1 |
| (3) falsification | no failing step found. Same parity is necessary: for χ_odd + χ_even mod 5 the phase gap is π/4, F(½+it) has no zero on [10,500] (min \|F\|/(\|L_1\| + \|L_2\|) = 0.383); mod 8 has no admissible pair. DH function mod 5: on-line fraction 0.90–0.96 in dyadic blocks to t = 2000 (b = (1,−1): 0.56–0.70), all ≥ 0.345 |
| (4) certificate | cert8.py rerun from a copy: identical log; independent recomputation (Gauss–Legendre/tanh-sinh quadrature for K, grids + golden section, 40 digits): p(10⁻⁴) = 0.672491254144, 2p − 1 = 0.344982508288 (cert ≥ 0.344982508) ✓ |
| reviewer (one, on the most delicate step: block hinge at a positive shift, full-range repair, BH summation) | GAP conditional on the isolation lemma, which was not pasted by design: DERIVATION §2 (audited) + the weight removal, now written out for L(s,χ) (FULLPROOF G4); steps 1–3 correct; its endpoint double count O(log T) adopted (logs/review.log) |

## Details

### 1. Per-character transfer (FULLPROOF §C–E)
- The Hilbert-space inequality (DERIVATION §2) uses only a finite conjugation-invariant multiset; Z_T^{(j)} is invariant under ρ ↦ 1 − ρ̄ (FE of L(s,χ_j) plus conjugation; same χ, same ordinate, same multiplicity). Nothing in §2 uses ζ.
- Pair correlation: Lem:GaPa/Lem:ConstructionKernelDirichlet (2610.26841 l.968–997) with r = f*f: r(0) + 2∫_0^1 αr(α)dα = ∫f² + ∬|u − v|f f = C(f), the ζ constant; GaPa's lemma has exactly BGST Lemma 5's form (weight 4/(4 − u²), same class of r), so La Lemma 3.2's proof (sole input BGST Lemma 5; re-derived in LamzAudit2 S3) transfers.
- Hinge B(ρ): Wu Prop. 3.1 at a = b = R_T = 2πρL_χ/log T (bounded, real, positive: allowed, l.390, l.398); |L(σ⁺+it)|²|A_χ(t)|² is the integrand; main term → c(R,R) = e^{−2R}c_1 (E1, symbolic), inf_P = c*(R); Littlewood on [σ⁺,2] × [T,2T] (right edge o(T), horizontal O(log T), Jensen) gives Σ_block(y − ρ)^+ ≤ (B(ρ) + o(1))(T/2π)log T with the ζ function B. Unsmoothing as Wu l.425–428 (E4).
- Hence p(c) = (2 − λ − C_MT − α_ρB(ρ))/(1 − λ) is the ζ number for every χ (cert8 calls cert.hinge_p) ✓.

### 2. Reduction, limits, p ∣ q (FULLPROOF §B, §G2(d))
- c-isolated (inside the block, threshold 2πc/log T) excludes H_j-membership once (log log T)^{1/4} > 1/c; only the first/last zero of a block can be lost: |H_j| ≤ (1 − δ_j + o(1))(T/2π)log T.
- Limits: for fixed (q, χ_j, kernel ε, θ < 4/7, P, ρ, c, BH's ε) let T → ∞, then P → sinh(bx)/sinh b, θ → 4/7, kernel ε → 0, BH's ε → 0; each limit is of lower bounds valid for every admissible choice (DERIVATION §4, §8.4(iv)) ✓.
- p ∣ q: L(1 ± z, χ_0) carries φ(q)/q, the Euler product of the h,k-sum over p ∤ q carries q/φ(q); they cancel in each of Wu's two terms (E1; numerically M_5/M_1 − 1 = −1.1% at y = 10³, −0.55% at 10⁶) ✓.

### 3. Falsification (FULLPROOF §F)
- Smallest failing case: opposite parity. Mod 5, χ_quartic (odd) + χ_quadratic (even): the gamma phases differ by π/4 (mod π), so F(½+it) = 0 forces L_1 = L_2 = 0; numerically no zero on [10, 500]. Mod 8 (χ_8 even, χ_−8 odd) has no admissible pair. The parity hypothesis is in the statement (Lamzouri l.249, ours) and is needed.
- Same parity, mod 5 (the smallest admissible case): Davenport–Heilbronn F (b = (1,1); its known off-line zero 0.808517 + 85.699348i reproduced, |F| = 1.3e−6) has 90–96% of the main-term count as sign changes in each dyadic block up to 2000 ((1000,2000]: 1011 of 1124.5); b = (1,−1): 56–70%, rising. Calibration (L_a alone): sign changes = main term ± 0.4 in every block. Consistent with ≥ 34.5% (a liminf claim; numerics only check the setting).

### 4. Certificate (FULLPROOF §E5–E6)
- cert8.py copied to .scratch/lcaudit/run and run under memcap: output identical to IsoZero/logs/cert8.log.
- Independent method: C(f_0) by quadrature = closed form to 15 digits; c = 10⁻⁴: λ = 1.22382e−7 (β_T ≈ −2/3, so λ_C = 2(1 − K(2c))), ρ_max = 1.2236206, α_ρ(1.2036) = 0.6702740537 (y* = 1.38857), B = 1.403822e−5, p = 0.672491254144, 2p − 1 = 0.344982508288; c = 10⁻³: 0.344974576351. Both match cert8 to all printed digits.

## New findings
1. Citation gap, repaired (FULLPROOF §G): 2610.26841 l.967 calls its Lem:GaPa "a dyadic variant of their Lemma 3.1"; GaPa prove only 0 < γ, γ′ ≤ T (normalisation log(qT)/2π). DERIVATION §8.2 inherits this. Repair: full-range isolation proportion (GaPa as stated), full-range hinge assembled from block hinges with ρ_k = ρ log T_k/log T, and Pro:BombieriHejhal summed over dyadic blocks (endpoint double count O(log T)). The same repair fixes Lamzouri's own 31%.
2. Wu, lower-order slips (no effect): (1.10) (2^{1+α+β} − 1)/(1+α+β) should be (2^{1−α−β} − 1)/(1−α−β) (only this reproduces (1.12)'s −1 + 2 log 2); (x3.4)/(1.12) "+ c_q" should be "+ 2c_q" (Laurent expansion of L(1 ± z, χ_0); classical mean square of L(½+it,χ)); (1.12) "2πHK" for 2πhk.
3. Cosmetic: DERIVATION §8.2 cites "BCHB's (1.6), l.152–155"; in Wu v2 it is l.147–149. cert8's ρ = 1.2036 is valid but not optimal (ρ = 1.210: +7e−8).


## Not re-derived here (residual risk, all inherited by Lamzouri's 31% as well)
- Bombieri–Hejhal (Duke 80 (1995)) Theorem B and Lemmas 7, 9, used inside Pro:BombieriHejhal (its reduction was checked as cited in IsoZero/AUDIT.md; not re-audited here).
- Proofs of Wu's Prop. 3.1 (Math. Z., refereed) and GaPa's Lemma 3.1 (Glasnik, refereed): statements authenticated in TeX / accepted manuscript only.
- (Done, not residual: La Lemma 3.2's weight removal redone for GaPa's normalisation log(qT)/2π, FULLPROOF G4; it uses no ζ-specific input.)
- Numerics stop at t = 2000 (DH blocks), a sanity check of the setting only.

## Files
FULLPROOF.md (chain, source lines, repair §G); mainterm_sym.py, mainterm_num.py (Wu's main term: symbolic, exact sums); indep_p.py (independent p); zeros_F.py, zeros_F2.py (zeros of F, parity test); logs/: cert8_rerun.log, mainterm_sym.log, mainterm_num.log, indep_p.log, zeros_F.log, zeros_F2.log, review.log, peaks.txt (memory/time; one indep_p attempt was killed by memcap at a sampled 1.9 GB and rerun chunked at 123 MB).


## Dependent: LinComb3 (2026-11-20, follow-up; LinComb3 DERIVATION §5, §8, §9.3; LC3Audit REPORT l.49)
- s₀ (Lemma M, §9.3): APPLIES. Lemma M can run on the full-range multiset {−i(ρ − ½)log(qT)/2π : 0 < γ ≤ T}. GaPa's published Lemma 3.1 already sums over all pairs, on and off the line, at the complex arguments i(ρ − ρ′)log(qT)/2π. The G4 weight removal 1/W(u) = 1 − u²/4 = 1 + π²z²/log²(qT) holds for complex u, so tr A² = Σ K_ε(z_ρ − z_ρ′)² = (C(f_ε) + o(1))N includes the off-line pairs. Input (P) follows by summing Selberg's S(t,χ) bound over dyadic pieces. The bound spn ≤ N + o(N) and thm:line-general hold for any multiset, so δ^R_j ≥ s₀ − η over the full range.
- p (stage 1, δ_j = p(c)): APPLIES. This is G2(a)–(c) verbatim.
- Theorem N / N^R summation (5.3, 8.2; needed for both): APPLIES, with one extra inequality.
  - Per dyadic block, let H_j(T_k) be the block's zeros minus the full-range isolated set. For T_k ≥ T^{1/2}, full-range κ-isolation implies isolation at κ/3 at the block's scale.
  - The capacity term of 4.6 then sums as Σ_k min(|H_j(T_k)|, cap_j(T_k)) ≤ min(Σ_k |H_j(T_k)|, Σ_k cap_j(T_k)) ≤ (min(1 − δ_j, 1/N) + o(1))N_T.
  - So 2s₀ − 1 = 860123/2478195 and 3s₀ − 2 = 51087/2478195 stand without the dyadic GaPa.

RH relevance: N8, N12
