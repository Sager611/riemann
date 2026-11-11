# Second-order term of the Báez-Duarte distance for the Bettin–Conrey–Farmer ramp

Notation: L = log N; s = 1/2 + it on the critical line; ρ = 1/2 + iγ zeros (RH); C = Σ_ρ 1/|ρ|² = 2 + γ − log 4π = 0.0461914179.
Distance (V1, Báez-Duarte, L²(0,∞)): for c ∈ R^N,
I(c) = ‖χ_(0,1] − Σ_{k≤N} c_k {1/(kx)}‖²_{L²(0,∞)} = (1/2π) ∫ |1 + ζ(s) Σ c_k k^{−s}|² dt/|s|²
(Mellin: χ ↦ 1/s, {1/(kx)} ↦ −k^{−s} ζ(s)/s for 0 < Re s < 1; Plancherel on Re s = 1/2).
BCF ramp: c_k = −μ(k)(1 − log k/L), i.e. 1 − ζ V_N with V_N(s) = Σ_{n≤N} μ(n)(1 − log n/L) n^{−s}; I_N := I(c).
d_N² = min_c I(c) ≤ I_N.

## Hypotheses (H)
RH; all zeros simple; BCF's condition Σ_{|γ|≤T} |ζ'(ρ)|^{−2} ≪ T^{3/2−δ} (arXiv 1211.5191 eq. (cond));
(H′) uniformity of the O(N^{−ε}) remainders: the contour shift in (B) across infinitely many zeros (horizontal segments → 0 along heights T_n → ∞) and the bounds
O(N^{−ε}) for the leftover line integrals in (B), (B′), with N-independent constants, as used (not re-proved) in BCF's proof via their Lemma 2 and Lindelöf.
The interchange Σ_{ρ1,ρ2} ↔ ∫ in (D) is NOT an extra hypothesis: it is proved below from (H) by domination.

## Theorem (conditional on H)
L² I_N = C L + J + D(N) + O(N^{−ε}) for every fixed 0 < ε < 1/2, where
- J = (1/2πi) ∫_(c) (ζ'/ζ)(s) (ζ'/ζ)(1−s) ds/(s(1−s)), any 0 < c < 1/2,
  = π²/12 − γ² − 2γ₁ + (1+γ) log 2π − Σ_{n≥1} (ζ'/ζ)(2n+1)/(2n(2n+1))
  = 0.63492080058374799 + 2.89872849930409450 + 0.02901415664852516 = **3.5626634565363677** (γ₁ = Stieltjes constant; mpmath 30 digits, s4_rate.py);
- D(N) = Σ_{ρ1,ρ2} K(ρ1,ρ2) N^{ρ1+ρ2−1}/(ζ'(ρ1) ζ'(ρ2)),  K(ρ1,ρ2) = (1/2πi) ∫_(c) ζ(s)ζ(1−s) ds / ((ρ1−s)²(ρ2−1+s)² s(1−s)), 0 < c < 1/2.
  N enters only through N^{ρ1+ρ2−1} = N^{i(γ1+γ2)}: D is real, bounded and almost periodic in log N, with mean
  D₀ = Σ_ρ K(ρ, 1−ρ)/|ζ'(ρ)|² (pairs ρ2 = ρ̄1) and oscillating part Ω(L) = Σ_{γ1+γ2≠0} … e^{iL(γ1+γ2)} (dominant frequencies γ_b − γ_a = 6.887, 3.989, 10.876, …).
Hence d_BCF²(N) := I_N = C/log N + (J + D(N))/log² N + O(N^{−ε}/log² N): **no constant C₂ exists unless Ω ≡ 0; the mean second-order coefficient is J + D₀**.

## Proof (BCF's decomposition, evaluated instead of bounded)
BCF Lemma 1 (unconditional): L(1 − ζV_N)(s) = ζ'/ζ(s) − E(s) − ζ(s)F_s(1/N), E(s) = ζ(s) Σ_ρ N^{ρ−s}/(ζ'(ρ)(ρ−s)²), F_s(1/N) ≪ N^{−5/2} on the strip.
Then L² I_N = (1/2πi) ∫_(1/2−ε) H(s)H(1−s) ds/(s(1−s)), H = ζ'/ζ − E − ζF (the full integrand is regular in 0 < Re s < 1, so the shift from 1/2 is free). Expand:
(A) ζ'/ζ · ζ'/ζ term = J (N-free; c-independent since 0 < Re s < 1/2 has no poles under RH).
    Shift left to Re s → −∞ (on far-left lines (ζ'/ζ)(1−s) ≪ 2^{Re s}, (ζ'/ζ)(s) ≪ log|s|, so the line integral → 0; avoid s = −2n).
    Poles: s = 0 double [(ζ'/ζ)(1−s) = 1/s + γ + O(s), 1/(s(1−s)) = (1 + s + …)/s, (ζ'/ζ)(s) = log 2π + α₁ s + …] → residue α₁ + (1+γ) log 2π;
    s = −2n simple [(ζ'/ζ) residue 1] → (ζ'/ζ)(2n+1)/(−2n(2n+1)). α₁ = (ζ'/ζ)'(0) = π²/12 − γ² − 2γ₁ from the functional equation
    ζ'/ζ(s) = log 2π + (π/2)cot(πs/2) − ψ(1−s) − ζ'/ζ(1−s) and ζ'/ζ(1+u) = −1/u + γ − (γ² + 2γ₁)u + O(u²).
(B) −(1/2πi)∫_(1/2−ε) (ζ'/ζ)(1−s) E(s) ds/(s(1−s)). Move to Re s = 1/2 + ε. At s = ρ, u = s − ρ: (ζ'/ζ)(1−s) = −1/u + a* + O(u),
    ζ(s)N^{ρ−s}/(ζ'(ρ)(ρ−s)²) = u^{−1}(1 + (a − L)u + …), g = 1/(s(1−s)) = g(ρ) + g'(ρ)u + …, with a = ζ''(ρ)/(2ζ'(ρ)), a* = ζ''(1−ρ)/(2ζ'(1−ρ)).
    Residue = g(ρ)(L + a* − a) − g'(ρ) = |ρ|^{−2}[L − i Im(ζ''/ζ'(ρ)) − 2iγ/|ρ|²] (RH: 1−ρ = ρ̄, a* = ā, g(ρ) = |ρ|^{−2}, −g'(ρ) = (1−2ρ)|ρ|^{−4}).
    The bracket minus L is purely imaginary and odd under ρ ↦ ρ̄, so Σ_ρ Res = L Σ 1/|ρ|² = CL exactly.
    [BCF (1211.5191, proof of Thm 1) print the residue as |ρ|^{−2}(L − ½ζ''/ζ'(ρ) + χ'/χ(ρ) + (1−2ρ)/|ρ|²). The functional equation
    ζ'/ζ(s) = χ'/χ(s) − ζ'/ζ(1−s) gives a + a* = χ'/χ(ρ), so a* − a = χ'/χ(ρ) − ζ''/ζ'(ρ): the coefficient of ζ''/ζ' is −1, not −½.
    On RH Re(ζ''/ζ'(ρ)) = 2 Re a = χ'/χ(ρ), so the correct bracket has real part exactly L; with −½ it would acquire ½χ'/χ(ρ) ≈ −½ log(|γ|/2π),
    i.e. a spurious constant ½ Σ_ρ χ'/χ(ρ)/|ρ|² ≈ −0.05. Harmless for BCF's O(1/log² N) bound; it matters for the constant.]
    The leftover line integral on Re s = 1/2 + ε has |N^{ρ−s}| = N^{−ε}: O(N^{−ε}) (BCF Lemma 2 + Lindelöf).
(B') −(1/2πi)∫_(1/2−ε) (ζ'/ζ)(s) E(1−s) ds/(s(1−s)) = (s ↦ 1−s) the same leftover integral on Re s = 1/2 + ε: O(N^{−ε}).
(D) (1/2πi)∫_(1/2−ε) E(s)E(1−s) ds/(s(1−s)): N^{ρ1−s} N^{ρ2−1+s} = N^{ρ1+ρ2−1} is s-free ⇒ D(N) as stated, once Σ_{ρ1,ρ2} and ∫ may be swapped.
    Domination (Fubini–Tonelli): put η = ε = 1/2 − c, s = c + it, S(t) = Σ_ρ |ζ'(ρ)|^{−1}/(η² + (t−γ)²) (finite: Σ 1/(|ζ'(ρ)||ρ|²) < ∞ by BCF).
    |N^{ρ1−s}| = N^{η}, |ρ1 − s|² = η² + (t−γ1)²; |N^{ρ2−(1−s)}| = N^{−η}, |ρ2 − (1−s)|² = η² + (t+γ2)²; S(−t) = S(t) (ρ ↦ ρ̄). So the double series is dominated,
    uniformly in N, by |ζ(s)ζ(1−s)| S(t)²/|s(1−s)|. RH ⇒ Lindelöf ⇒ |ζ(s)ζ(1−s)| = |χ(s)||ζ(1−s)|² ≪ (1+|t|)^{η+κ} (any κ > 0); |s(1−s)| ≍ (1+|t|)².
    Cauchy–Schwarz: S(t)² ≤ (Σ_ρ w_ρ)(Σ_ρ w_ρ/|ζ'(ρ)|²), w_ρ = 1/(η² + (t−γ)²), and local zero counting gives Σ_ρ w_ρ ≪ η^{−2} log(2+|t|).
    Hence ∫_R |ζζ| S²/|s(1−s)| dt ≪_η Σ_ρ |ζ'(ρ)|^{−2} ∫ (1+|t|)^{η+κ−2} log(2+|t|) w_ρ(t) dt ≪_η Σ_ρ (1+|γ|)^{−2+η+2κ}/|ζ'(ρ)|²
    (the Lorentzian has mass π/η near t = γ; the part |t−γ| > (1+|γ|)/2 contributes ≪ (1+|γ|)^{−2}), and by partial summation with BCF's condition
    the last sum converges iff η + 2κ < 1/2 + δ, true for η < 1/2 and κ < δ/2. So Σ_{ρ1,ρ2}|K(ρ1,ρ2)|/|ζ'(ρ1)ζ'(ρ2)| < ∞, the swap holds, and D is an absolutely
    convergent generalized trigonometric series in log N, hence Bohr almost periodic with mean D₀ (the pairs ρ1 + ρ2 = 1).
(F) terms with F_s(1/N): O(N^{−5/2+ε}). ∎ (proof sketch; the only unproved input is (H′), the uniformity of the O(N^{−ε}) remainders.)

## Consequences
1. The ramp approaches C from ABOVE: d_BCF² log N − C = (J + D(N))/log N + o(1/log N) > 0 for large N (J ≫ |D| numerically).
   It is far from asymptotic at accessible N: at N = 1280, I_N = 0.07560 vs C/L = 0.00646 (J/L² = 0.0697 dominates); C L > J only for L > J/C = 77, N > 10^33.
2. L²(0,1) variant (dropping the tail x > 1): subtract (Σ_k c_k/k)² = L^{−2}(Σ_{k≤N} μ(k) log(N/k)/k)² = L^{−2}(1 + O(N^{−1/2+ε})) ⇒ coefficient J − 1 + D(N) = 2.5627 + D(N).
   This is exactly what Veltri (Zenodo 22647730 v3, §'A numerical remark') tabulated as "λ·(1/2π)∫|1−ζV_N|²dμ": his c_N = 0.60431, 0.52511, 0.47279, 0.42631 (N = 100, 200, 400, 800)
   equal our L·(I_N − (Σc_k/k)²) to all 5 digits (logs/veltri_check.log), while L·I_N = 0.82191, 0.71358, 0.63980, 0.57582; direct quadrature on the critical line confirms
   I_100 ≥ ∫_0^400 part = 0.17767 > 0.13123 (s5_direct.py). So his K ≈ 2.5 is K − 1; for the stated integral K = J + D₀ ≈ 3.55. His open question (iv) ("identify the second-order term")
   is answered by the Theorem: J + D(N), with the residue sums he lists contributing nothing (their real parts cancel).
3. General smooth profile c_k = −μ(k)P(log k/L), P(0) = 1, P(1) = 0: Perron with F̂ gives near each zero the local profile h(v) = −∫_0^1 P'(x)e^{−ivx}dx (v = L(t−γ)),
   so the first-order constant is C∫_0^1 P'² ≥ C (Cauchy–Schwarz; equality only for the ramp), and the non-local term is −P'(0) ζ'/ζ(s)/L, which contributes P'(0)² J at second order.
   The ramp is the unique first-order-optimal fixed profile, and its slope P'(0) = −1 at k = 1 costs J/L². The optimal coefficients (s4_rate.py) are flat at small k
   (−c_k/μ(k) = 0.981, 0.969, 0.964, 0.903 for k = 1, 2, 3, 5 at N = 1280; the ramp has 1, 0.903, 0.846, 0.775) and keep weight P(1) ≈ 0.16–0.2 at k ≈ N: they remove the J-term by an
   N-dependent boundary layer, which is why d_N² is 12× smaller than I_N at N = 1280 and why the optimal second-order coefficient is ≈ 0 instead of ≈ 3.55.
   Where the derivation stops: the optimal second-order term needs the joint optimisation of an N-dependent boundary layer (k = O(1)) and the endpoint weight P(1) = o(1)
   against the zero sums; not done. Numerically, Möbius × polynomial(deg 10) + free c_1..c_12 gives L²d² − CL = +0.0077 at N = 1280 (optimum −0.0049).

## Numerical check (rigorous Gram values of I_N from s2_rig.py; pair constants from s6_pairs.py / s6b_D0.py)
- R(N) := log² N·I_N − C log N − J: mean −0.0072 (N ∈ [100,200]), −0.0147 ([200,400]), −0.0131 ([400,800]), −0.0149 ([800,1280]); sd 0.015 → 0.008 (decaying: the O(N^{−ε}) remainder).
- D₀ partial sums over the first m zeros (each term = 2 Re K(ρ,ρ̄)/|ζ′(ρ)|² for the pair ρ, ρ̄): −0.00123 (m=1), −0.00402 (6), −0.00509 (10), −0.00749 (30), −0.00803 (40, γ ≤ 123);
  term·γ²/log(γ/2π) fluctuates in [−0.41, −0.14] with mean −0.28 ⇒ tail ≈ −0.28 × Σ_{γ>123} log(γ/2π)/γ² = −0.28 × 0.0217 = −0.0061 ⇒ **D₀ = −0.0142 ± 0.003**.
- Off-diagonal (oscillating) amplitudes |K/(ζ′ζ′)| ≤ 1.8e−4 for the first 6 zeros (frequencies ±6.887, ±3.989, ±2.51, ±5.414): Ω is ≲ 1e−3, below the remainder at N ≤ 1280.
- J + D₀ = 3.5485 ± 0.003 vs data mean 3.5478 (N ∈ [800,1280]).
- Reviewer (algebra only): residues of Z at 0 and −2n, α₁, the zero residue with a + a* = χ′/χ(ρ) (BCF's −½ should be −1), and the contour orientations: all CORRECT;
  contour-limit and convergence justifications not assessed.
