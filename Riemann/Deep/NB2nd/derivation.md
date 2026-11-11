# nb2nd: second-order term of the optimal Báez-Duarte distance — derivation and tests

Notation as nbrate (Riemann/Deep/NBRate/derivation.md). V1 only: d_N² = min_c ‖χ_(0,1] − Σ_{k≤N} c_k ρ_{1/k}‖²_{L²(0,∞)} = 1 − bᵀG_N⁻¹b (Vasyunin Gram).
⟨F, H⟩ := (1/2π)∫ F(½+it) conj H(½+it) dt/(¼+t²); F_A := 1 − ζA, A = −Σ c_k k^{−s}; I(A) = ‖F_A‖². L = log N, C = 2 + γ − log 4π = 0.0461914179,
J = 3.5626634565 (nbrate). Ramp V_N = Σ_{n≤N} μ(n)(1 − log n/L) n^{−s}, R_N := 1 − ζV_N, I_N := ‖R_N‖², R(N) := L²I_N − CL − J (= D(N) + remainder).
Tags: [R] rigorous; [R-RH] rigorous under RH; [H] heuristic / proof sketch at the level of nbrate's theorem (needs RH, simple zeros, BCF's
Σ_{|γ|≤T}|ζ′(ρ)|⁻² ≪ T^{3/2−δ}, and (H′) = uniformity of the O(N^{−ε}) leftover line integrals); [N] numerical (double precision, cond G_640 = 1.6e6).

## 1. Exact reduction [R]
r := G_N c^ramp − b_N has r_n = ⟨R_N, ζn^{−s}⟩, and d_N² = I_N − rᵀG_N⁻¹r. Hence, with S_N := (Lr)ᵀG_N⁻¹(Lr),
  C₂(N) := L²d_N² − CL = J + R(N) − S_N.                                                         (1)
The question "second-order term of the optimum" is the asymptotics of S_N (nbrate: S_N ≈ 3.55, J + D₀ ≈ 3.55).

## 2. The head functional ℓ [R-RH]
ℓ_n := ⟨ζ′/ζ, ζn^{−s}⟩ = (1/2πi)∫_(½)(ζ′/ζ)(s) ζ(1−s) n^{s−1} ds/(s(1−s)). The integrand is regular on the line (ζ(1−s) vanishes at every ρ) and,
under RH, in 0 < Re s < ½. Shift to Re s → −∞ (the left line integral → 0: (ζ′/ζ)(s) ≪ log|s| there, ζ(1−s)n^{s−1} bounded). Residues:
s = 0 (double pole: ζ(1−s) = −1/s + γ + … times 1/s; (ζ′/ζ)(s) = log 2π + α₁s + … regular), s = −2k (simple):
  ℓ_n = −[log 2π·(1 + log n − γ) + α₁]/n − Σ_{k≥1} ζ(2k+1)/(2k(2k+1) n^{2k+1}),   α₁ = π²/12 − γ² − 2γ₁ = 0.6349208006.     (2)
ℓ_1..ℓ_6 = −1.754554, −1.369816, −1.151332, −0.993129, −0.875599, −0.785098.
- [N] Check: L r_n → ℓ_n (BCF Lemma 1: L R_N = ζ′/ζ − E − ζF): at N = 40…640, |L r_n − ℓ_n| ≤ 0.003 for n ≤ 6 (t1).
- [H, formal] Σ_n μ(n) log n · ℓ_n = J term by term (Σμ(n)log n/n = −1, Σμ(n)log²n/n = −2γ, Σμ(n)log n·n^{−2k−1} = ζ′/ζ²(2k+1) reproduce
  nbrate's closed form of J exactly): formally ζ′/ζ = ζ·Σμ(n) log n n^{−s} ("B = ζ′/ζ² would cancel J"); the series is conditionally
  convergent (partial sums 2.30, 2.17, 4.27, 3.39 at K = 10, 100, 300, 640) and ‖ζ·ζ′/ζ²‖ = ‖ζ′/ζ‖ = ∞, so this cancellation is not available.

## 3. Fixed (N-independent) small-k corrections [H]
Claim. For every fixed Dirichlet polynomial B = Σ_{n≤K} β_n n^{−s} (A = V_N + B/L, i.e. c ↦ c^ramp − β/L):
  L²I(V_N + B/L) − CL = J(B) + D(N) + O(N^{−ε}),   J(B) := J − 2ℓ(B) + ‖ζB‖²  (ℓ(B) = Σβ_nℓ_n, ‖ζB‖² = βᵀG_Kβ).      (3)
Proof sketch. Exactly, L²I(V_N + B/L) = L²I_N − 2⟨LR_N, ζB⟩ + ‖ζB‖² [R], and ⟨LR_N, ζB⟩ = ℓ(B) − ⟨E, ζB⟩ + O(N^{−5/2}). In
⟨E, ζB⟩ = (1/2πi)∫_(½) E(s)ζ(1−s)B(1−s) ds/(s(1−s)) move to Re s = ½ + ε: E has simple poles (residue 1) at the zeros, ζ(1−s) simple zeros
there, so no residues; on the new line |N^{ρ−s}| = N^{−ε}: O(N^{−ε}) under (H′) (as BCF's leftover integrals). Combine with nbrate's theorem. ∎
So B is invisible to the zero residues (the C L term) and to the pair term D(N); it only re-weights the head J. Best B of length K:
  min J(B) = J − Σ_K,   Σ_K := ℓ_Kᵀ G_K⁻¹ ℓ_K   (sup over B of 2ℓ(B) − ‖ζB‖²).
[N] Σ_K (t1, t5): 2.4419 (K=1), 2.9954 (2), 3.3108 (3), 3.4753 (5), 3.5073 (10), 3.5161 (20), 3.5189 (40), 3.5201 (80), 3.5215 (160),
3.5224 (320), 3.5230 (640). **J − Σ_K ≥ 0.0397 for every K ≤ 640.** Fits on K ∈ [40,640]: Σ_∞ − A/log K (Σ_∞ = 3.5287, A = 0.036, rms 7.7e−5);
Σ_∞ − A K^{−0.43} (Σ_∞ = 3.5249, rms 7.2e−5); a + b log log K (divergent, rms 1.1e−4, worse). If Σ_K converges: **J − Σ_∞ ∈ [0.034, 0.038]**.
Test of (3) [N] (t7, mean ± sd over N):
| K | J − Σ_K | L²I(V_N+B*_K/L) − CL, N∈[250,450] | N∈[450,640] | (construction) − (J − Σ_K) |
|---|---|---|---|---|
| 5 | 0.0874 | 0.0745 ± 0.0012 | 0.0744 ± 0.0015 | −0.0129, −0.0130 |
| 20 | 0.0465 | 0.0333 ± 0.0006 | 0.0337 ± 0.0005 | −0.0132, −0.0128 |
| 60 | 0.0428 | 0.0293 ± 0.0006 | 0.0298 ± 0.0005 | −0.0135, −0.0130 |
The last column is D₀ + (mean remainder): **D₀ = −0.0131 ± 0.0003**, inside nbrate's −0.0142 ± 0.003 (40 zeros + tail law), now from Gram data.
Side result: the fixed head correction also removes the ramp's scatter (sd 0.0005 vs 0.013 for R(N)): that scatter is the head–E leftover
2Re⟨(ζ′/ζ)_smooth, E⟩ = O(N^{−1/2}) (t6: R(N) ≈ 0.29 N^{−1/2} cos(γ₁ log N + φ) + …, R² 0.95), which B cancels with (ζ′/ζ)_smooth.
Shape of the head optimizer β* = G_K⁻¹ℓ_K (t2, t5): β*_1 = −0.152 (so P_opt(1) = 1 − 0.152/L: 0.979 at N = 1280; nbrate measured 0.981);
β*/μ = 0.63 (n = 2…7), 0.53 (8–15), 0.49 (16–31), 0.43 (32–63), 0.35 (64–127), 0.28, 0.21, 0.18 (n ≈ K = 640): not a boundary layer but a slowly
decaying profile with end value ≈ 0.19 (K = 160, 320, 640: 0.205, 0.190, 0.186); Σβ*_n/n = −1.018 ≈ −1, i.e. A(1) = V_N(1) + B(1)/L ≈ (1 − 1)/L:
the optimum kills the 1/x tail on x > 1, which is why V1 ≈ V2 (L²(V2 − V1) = 2.6e−4 at N = 640).
Consequence: **no N-independent correction cancels J.** The best fixed construction has second-order mean J − Σ_∞ + D₀ ∈ [0.021, 0.025] (> 0),
the K = 60 one 0.0295; the optimum has ≈ −0.002 (N ≤ 640).

## 4. Endpoint mechanism: a jump P(1) = j/L [H]
Profile P(x) = 1 − (1 − j/L)x, A_j = Σ_{n≤N} μ(n)P(log n/L) n^{−s} (a tilt with P(1) = j/L, i.e. A_j = V_N + j(M_N − V_N)/L, M_N sharp Möbius sum).
Perron with p(y) = (1 − j/L)y + j (y = log(N/n)), p̃(z) = (1 − j/L)/z² + j/z: L(1 − ζA_j) = (1 − j/L)ζ′/ζ − E_p, E_p = ζΣ_ρ N^{ρ−s}p̃(ρ−s)/ζ′(ρ).
Residue at s = ρ of E_p(s)(ζ′/ζ)(1−s)g(s) (g = 1/(s(1−s)), u = s − ρ, a = ζ″/2ζ′(ρ)): E_p = (1−j/L)/u + (1−j/L)(a − L) − j + O(u),
(ζ′/ζ)(1−s) = −1/u + ā + O(u) ⇒ Res = (1 − j/L)[g(ρ)(L + ā − a) − g′(ρ)] + g(ρ)j; real part g(ρ)L (ā − a, g′(ρ) imaginary under RH). The cross term carries
the factor (1 − j/L): total (1 − j/L)CL. So
  L²I(A_j) − CL = (1 − j/L)²J + D_p(N) − jC + O(N^{−ε}),  D_p = D(N) + j²X(N) + (oscillating j-linear terms),                (4)
X(N) := ‖1 − ζM_N‖² (sharp-cut norm; the cross term between the 1/z² and 1/z profiles is pointwise imaginary on the line, mean 0).
Independent check of the linear gain −jC (local Fourier, near each zero): F ≈ ((1−j/L)/L)(1 − e^{−iLτ})/(iτ) + (j/L)e^{−iLτ};
∫|F|² = 2π(1 − j/L)/L + (j/L)²·∞_local (the j² part is the non-local X), using ∫(e^{iLτ} − 1)/(iτ)dτ = π: the jump sits at the band edge u = L and
overlaps the ramp's spectrum by half a delta. [In Burnol's MA(1) toy Q = 1 − z the jump sits outside the band (degree N + 1): ramp exactly 1/N,
optimum 1/(N + 2), i.e. second-order 0 vs −2, both ends gaining 1 by global tilts a^opt − a^ramp = (2k − N)/(N(N+2)); ζ has no band extension.]
[N] (t3, N ∈ [400,640]): measured L²·2I₁ = −1.146 ± 0.218 vs −C − 2J/L = −1.15 at N = 640; L²I₂ = 0.161 ± 0.040 vs X + J/L² = 0.076 + 0.085.
X(N): mean 0.064–0.079, sd 0.034–0.055, range [0.033, 0.30] (N ∈ [100,640]). Jump-only asymptotic gain mean[C²/(4X(N))] = 0.009–0.010.
(The finite-N tilt also scales the head, −2jJ/L ≈ −0.55 j at L = 6.5: the separation of orders is poor because J ≈ 77 C.)

## 5. Optimum versus fixed constructions; where it stops
Write S_N = Σ_N + κ(N), Σ_N = ℓᵀG_N⁻¹ℓ (head problem with K = N), κ(N) = −2eᵀG_N⁻¹ℓ + eᵀG_N⁻¹e, e := ℓ − Lr (e_n = ⟨E, ζn^{−s}⟩ + O(N^{−5/2})).
By (1): C₂(N) = (J − Σ_N) + R(N) − κ(N)  [R]. Under (H), for N → ∞: Σ_N → Σ_∞ (if finite), mean R(N) → D₀, so
  C₂^opt (mean) = (J − Σ_∞) + D₀ − κ̄,   κ̄ := lim mean_N κ(N).                                                            (5)
- [N] κ(N) = +0.0275 ± 0.0168 (N ∈ [100,250]), +0.0311 ± 0.0109 ([250,450]), +0.0273 ± 0.0153 ([450,640]) (t8, every N; fit κ_∞ − a/L:
  κ_∞ = 0.032, no significant trend); corrected for the correlated ramp scatter, κ(N) − (R(N) − D₀) = 0.0255, 0.0281, 0.0284. For N-independent B, κ → 0 (§3,
  confirmed by t7); so κ̄ is entirely the gain of **N-dependent** corrections, located at k near N: with the head support 1..60 fixed, the extra gain
  S_N − S_N(1..60) is 0.026, 0.030, 0.0315 (same N ranges); the tilt direction alone recovers 43–60 % of it, tilt + support [N e^{−1/2}, N] 75–88 % (t4).
  The optimal profile excess w(y) = L(P_opt − ramp) at y = log(N/n) is 1.48 (y < 0.02), 1.3 (0.02–0.1), 1.2 (0.1–0.35), 1.05–1.25 (0.35–1.5),
  1.0 (2–3) at N = 640 (t2): a tilt j ≈ 1.2–1.4 plus the head profile (≈ 0.2), i.e. P(1) ≈ 0.2 (nbrate: 0.16–0.20); non-squarefree L·c_k near k = N:
  rms 0.59 (y < 0.5).
- Not controlled: κ̄. It is the mean over N of a nonlinear functional (the inverse Gram G_N⁻¹ on the N-dependent zero-sum vector e_N, coupled to the
  arithmetic vector ℓ) of almost periodic functions of log N (phases N^{iγ}); the endpoint minimisation adapts to these phases (jump cost X(N) varies
  by ×9), so κ̄ is not a diagonal zero sum like D₀, and no second-order asymptotics of G_N⁻¹ is known (Burnol's orthogonal vectors give first order
  only). **C₂^opt = 0 ⟺ κ̄ = J − Σ_∞ + D₀ ∈ [0.021, 0.027]** (J − Σ_∞ ≤ J − Σ_640 = 0.0397: Σ_K non-decreasing, Σ_640 in double precision, error < 1e−8): an identity
  between an arithmetic constant (Σ_∞, from ℓ and Vasyunin's Gram) and an endpoint zero-sum mean, for which no mechanism was found. The jump alone
  gives κ̄ ≥ 0.009 (asymptotic, heuristic), hence C₂^opt ≤ 0.018. If κ(N) has no trend beyond N = 640 (κ̄ ∈ [0.025, 0.033]):
  C₂^opt ∈ [−0.012, +0.002], central −0.005.
- Also not controlled: convergence of Σ_K (only K ≤ 640 computed; a log log K growth fits slightly worse but is not excluded). Heuristic for the
  rate: ζB must reproduce the regular part a_ρ = ζ″/2ζ′(ρ) of ζ′/ζ at each zero (Re a_ρ = −½ log(γ/2π) to 1e−4, t8), an NB-type local problem
  of cost |a_ρ|²/(|ρ|² log K) ⇒ A = Σ_ρ|a_ρ|²/|ρ|² = 0.041 (100 zeros) … ≈ 0.075 (with tail) vs fitted 0.036: right order, factor ≈ 2 open.
  The poles themselves cannot be exploited: approximating 1/(s−ρ) by ζB with the one-sided family B = 1/(ζ′(ρ)(s−ρ+η)²) has local f.p. cost
  (1/η)∫(2 − v²)/(1+v²)²dv = +π/(2η) > 0; a two-sided (principal-value) approximant would gain −¼‖Σ_ρx^{−ρ}/ρ‖² = −∞ but is, heuristically,
  not reachable by ζ·(Dirichlet polynomials) (x-space: the steps σ_m of Σc_n⌊y/n⌋ would have to track (ψ(m) − m)/2; the Franel-type energy
  Σc_n{y/n} beyond the tracked range K appears to grow like log² K against a gain (C/4) log K; numerically Σ_K grows ≤ 6e−4 per doubling).

## 6. Data and oscillating part [N]
C₂(N) = L²d_N² − CL (optimum): +0.0055 ± 0.0025 (N ∈ [50,100]), +0.0050 ± 0.0018 ([100,200]), −0.0015 ± 0.0018 ([200,400]), −0.0015 ± 0.0023
([400,640]) (t2, double; agrees with nbrate's certified values: N = 640: −0.0055); nbrate certified [800,1280]: −0.0062 ± 0.0009; 2611.03924 (V2,
quadrature, V2 − V1 ≈ 3e−4 in these units): +0.0028 (2560), −0.0091 (5120). Mean over [200, 5120] ≈ −0.003 ± 0.004.
C₂ + C₃/L fit on [100,640]: C₂ = −0.021, C₃ = 0.12 (R² 0.35, unstable; nbrate: −0.024…−0.044, 0.13…0.25). Adding cos/sin of the pair frequencies
γ₂−γ₁ = 6.887, γ₃−γ₁ = 10.876, γ₃−γ₂ = 3.989: R² 0.75, amplitudes 0.0026, 0.0019, 0.0008; single-zero terms N^{−1/2}cos(γ log N) add little (R² 0.41
alone). So the optimum's oscillation (sd 0.003) has zero-pair frequencies (D(N)-type), amplitude ≈ 10× the ramp's off-diagonal D amplitudes
(≤ 1.8e−4, nbrate), consistent with an endpoint optimisation that adapts to N^{i(γ_b−γ_a)}; log N ∈ [4.6, 6.5] resolves no more than 3 frequencies.
Finite-N couplings are of the size of C₂ itself (κ(N) sd 0.011–0.017; head–E coupling S_N(1..K) − Σ_K = 0.019 at N = 640, K = 60), so N ≤ 5120
cannot fix the sign of C₂^opt.
