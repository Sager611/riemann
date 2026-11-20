BARRIER HOLDS: with the zeros' almost-periodic terms regressed out (random-amplitude zero lines, consistent out of sample, stable from K = 50 to 100) every full-window estimate of C₂ is negative, but C₂ = 0 is excluded at ≥ 3σ only if the 1/L expansion stops at C₃/L (C₂ = −0.019 ± 0.0035, N ≤ 5500, model-conditional); with C₄/L² admitted (unconstrained by N ≤ 5500: C₄ = 0.43 ± 1.7) C₂ = −0.007 ± 0.024, and a design projection puts 3σ at N ≈ 10⁵–10⁷.

# nb2fals: falsification of the nb2nd barrier (node N5, rate) — can per-N data fix the sign of C₂?

## Answer (Pyramid)
1. **No, not model-robustly at 3σ.** Once the zero terms are regressed out with random amplitudes, the oscillation that nb2nd blamed is no longer the main uncertainty; the slow 1/L expansion is. The constant left after the regression is degenerate with C₃/L + C₄/L² on u = log N ∈ [4.6, 8.6].
2. **The literal test (free amplitudes, K = 0…50 zeros) is numerically infeasible beyond K ≈ 5.** The least-squares designs are numerically singular: cond 8e12 (K = 5, N ≤ 640), 1e16 (K = 10, N ≤ 4000). The Landau–Pollak–Slepian effective-dimension heuristic explains this (it is not a theorem about exact rank): lines with frequencies ≤ Ω have only ≈ ΩT/π well-conditioned dimensions on a window of length T. K zeros give ≈ K² pair columns with Ω ≈ γ_K ≈ 2πK/log K, so free amplitudes need T ≳ (K log K)/2, i.e. N ≳ 10⁷ already for K = 10. Identifiable fits: K = 0 + C₃/L: −0.0193 ± 0.0062 (block bootstrap, N ≤ 4000); K = 5: −0.027 ± 0.003, but the unmodelled lines leak ±0.07 into it.
3. **Random amplitudes give stable, model-conditional errors.** Treating the line coefficients as Gaussian with the theory's height/gap shape (Gaussian-process GLS, REML) gives:
   - shape B, the close-pair-enhanced shape (amplitude ∝ 1/(γ_aγ_b(γ_b − γ_a))), is consistent out of sample: window-mean z = +0.54, +0.76, −0.61 (N ≤ 4000 runs) and +0.73, −0.47, +0.48 (N ≤ 5500 runs; the last predicts the extension N ≤ 4000 → (4000,5500]). It is also consistent over 8 windows. These overlapping holdouts test the fitted curve, not the coverage of the extrapolated constant;
   - the gap-free shape A is refuted: window spread 0.04 against claimed s.e. 0.001, out-of-sample z = +13;
   - stability in K: s.e. 0.0069 → 0.0072 from K = 50 to 100 (a sensitivity check; the tail of lines beyond K = 100 is not controlled);
   - no extra 1/L-modulated lines: their fitted scale goes to 1e−5.
   With C₃/L, N ≤ 5500, a REML-weighted sensitivity mixture over the gap exponent β (not a posterior) gives **C₂ = −0.0192 ± 0.0035** (5.5σ). Shape B alone gives −0.0194 ± 0.0061 (3.2σ). Both s.e. are conditional on the covariance model and the trend basis.
4. **The trend order decides the sign.**
   - The data want C₃: C₃ = +0.108 ± 0.039, and the no-trend model is worse by 2ΔlogL = 7.6.
   - C₄ is unconstrained: adding C₄/L² improves the fit by only 2ΔlogL = 0.07 (C₄ = 0.43 ± 1.7). This detects no C₄, but it does not show C₄ ≈ 0. With C₄ admitted, C₂ = −0.007 ± 0.024 (fitted shape) or −0.008 ± 0.044 (shape B).
   - So a 3σ sign holds only under a truncation that N ≤ 5500 cannot check. With C₃ refitted, c moves by ≈ +0.026 per unit of C₄: C₄ ≈ C₃ ≈ 0.1 moves it by 0.003, the fitted C₄ = 0.43 by 0.011, and the ±1.7 uncertainty of C₄ gives the ±0.044.
5. **One-sided head route.** Σ_N is non-decreasing, so C₂ ≤ (J − Σ_5500) + Q̄, with J − Σ_5500 = 0.03839 and Q = C₂(N) − (J − Σ_N) = R − κ. The inequality is rigorous only in the head part. Q̄ is a fitted limit, so the bound is as model-dependent as c. It depends on the trend allowed in Q:
   - Q trend-free (the data prefer a trend only at 2ΔlogL = 3.7): C₂ ≤ −0.0039 ± 0.0005;
   - Q + C₃/L: C₂ ≤ −0.015 ± 0.003 (fitted shape) or ± 0.006 (shape B);
   - Q + C₃/L + C₄/L²: ± 0.024–0.044, sign open.
   Same conclusion: the sign rests on the 1/L expansion of κ(N), not on its oscillation.
6. **N needed.** Data-free projection of the posterior s.e. with C₄ free:
   - fitted shape: 0.025 (N = 5500), 0.013 (2.2·10⁴), 0.0041 (1.6·10⁵), 0.0016 (3.3·10⁶);
   - shape B: 0.045, 0.024, 0.0085, 0.0036.
   This is a design calculation with the covariance and the trend basis held fixed, not a detection forecast. 3σ needs s.e. ≤ |C₂|/3:
   - |C₂| ≈ 0.019 (the C₃-only value): N ≈ 8·10⁴ (fitted shape) to 5·10⁵ (shape B);
   - |C₂| ≈ 0.007 (the C₄-model value): N ≈ 1·10⁶ to 1·10⁷;
   - |C₂| ≈ 0.004 (the trend-free Q value): N ≈ 6·10⁶ to 5·10⁷.
   A dense double Gram at N = 10⁵ takes 80 GB, so this is out of reach of the present method under the 300 MB cap.
7. **No identity replaces the data.**
   - The problem is one-sided: S_N = max_x [2x·Lr − xᵀG_N x], so constructions bound C₂^opt from above, while lower bounds need dual (Burnol-type) certificates, known only at first order.
   - Proving C₂ < 0 therefore needs an N-dependent construction with mean end gain κ_c > (J − Σ_K) + D₀ ≥ 0.0253 (the head part is rigorous for K ≤ 5500 and D₀ = −0.0131 is from nb2nd; the optimum has κ ≈ 0.030 at N ≤ 5500, extrapolating to 0.040–0.044).
   - Candidates: the jump alone gives 0.009 (nb2nd §4); tilt plus end support [Ne^{−1/2}, N] recovers 75–88 % of κ at N ≤ 640 (nb2nd t4). Its mean is again an inverse-Gram functional of the end block, so a signed decomposition would only move the barrier.

## Data (f0_perN.py dense N ≤ 4000; f0p_packed.py packed N ≤ 5500)
- One script merges NB2nd t0 (Gram), t2 (d_N², R), t5 (ℓ, Σ) and t8 (κ), with the same formulas; per N one matrix-vector product with the leading block of the in-place Cholesky factor (U_N leading block of U): S_N = L²‖U_N c_N − U_N^{−T}b_N‖², I_N = 1 − 2c_N·b + ‖U_N c_N‖², C₂(N) = L²d_N² − CL, κ(N) = S_N − Σ_N, identity C₂ = (J − Σ_N) + R − κ to 6e−13.
- Checks: N ≤ 640 identical to NB2nd t8 (κ means 0.0275/0.0311/0.0273); d_N² at 640, 695, 1280 = NBRate certified to 5e−15; C₂(N) on [800,1280]: mean −0.00617, sd 0.00088 (NBRate certified −0.0062, 0.0009); packed = dense to 1e−15 (N = 640…4000).
- Size: N = 4000 dense 65 s, peak 204 MB; N = 5500 packed 158 s, peak 227 MB (N = 6000 would need ≈ 280 MB: stopped at 5500). λ_min(G_5500) = 5.5e−8, cond 1.1e8 (≤ for all N ≤ 5500 by interlacing): double precision suffices.
- Per-N values: data/perN_640.csv (required table), data/perN_5500.csv.gz (N ≤ 5500): C₂, κ, J − Σ_N, R, d_N².
- Shape of the data: smooth in u = log N (step |ΔC₂| rms 3e−4 at N ≈ 300, 2e−5 at N ≈ 3500): all uncertainty is low-frequency. Periodogram of C₂ (Hann, [100,4000], resolution 2π/T = 1.70): peaks at 2.6 (amp. 0.0038), 11.1, 16.4; κ: 14.0 (0.014, γ₁), 20.8, 25.1.
- Head: Σ_N non-decreasing ⇒ J − Σ_∞ ≤ J − Σ_5500 = 0.03839 (rigorous up to double rounding).

## Fits
Model (theory: NB2nd §6, NBRate D(N)): y(u) = c + Σ_{p≤P} C_{p+2}u^{−p} + Σ_{a<b≤K}[p_ab cos + q_ab sin]((γ_b − γ_a)u) + Σ_{k≤K} e^{−u/2}[a_k cos + b_k sin](γ_k u), u = log N, N ≥ 100.
A. Free amplitudes (OLS, weights 1/N; f1_fits.py), errors = moving-block bootstrap of the residual in u (block b), leakage = transfer of the unmodelled lines of zeros ≤ 100:
| target, window | K = 0 | K = 5 | K = 10 | K = 20 | K = 50 |
|---|---|---|---|---|---|
| C₂ + C₃/L, [100,640] | −0.0225 ± 0.0073 | cond 8e12: c ≈ −4e7 | cond 9e16 | cond 5e17 | p > n |
| C₂ + C₃/L, [100,4000] | −0.0193 ± 0.0062 (leak 0.008) | −0.0270 ± 0.0032 (leak 0.069) | cond 1e16 | cond 9e16 | 1225 lines, spacing 0.1 ≪ 1.7 |
| C₂ + C₃/L, [100,5500] | −0.0206 ± 0.0059 | −0.0178 ± 0.0005 (fails out of sample) | — | — | — |
Bootstrap block length matters:
- 20-N blocks give 0.0043 ([100,640]), 0.0013 ([100,4000]) and 0.0010 ([100,5500]).
- b = 0.6 in u gives 0.0073, 0.0062 and 0.0059.
- The residual is correlated over Δu ≈ 1, so 20-N blocks understate the error by 2–6×.

Out of sample (f7), the window-mean error of the prediction:
- K = 0, fit ≤ 400 → (400,640]: +0.0016 (sd of y 0.0022);
- K = 0, fit ≤ 1280 → (1280,5500]: +0.0050 (sd 0.0052);
- K = 0, fit ≤ 4000 → (4000,5500]: −0.0021 (sd 0.0013);
- K = 5, fit ≤ 4000 → (4000,5500]: +0.0156 (sd 0.0013). The K = 5 fit is overfitted; its bootstrap s.e. is meaningless.
Why OLS fails (Landau–Pollak–Slepian effective-dimension heuristic): on a window of length T, lines with frequencies ≤ Ω have only ≈ ΩT/π well-conditioned dimensions; K zeros bring ≈ K² pair columns with Ω ≈ γ_K ≈ 2πK/log K, so a free-amplitude fit needs T ≳ (K log K)/2: T ≥ 11.5 (N ≥ 10⁷) for K = 10, T ≥ 98 for K = 50. At N ≤ 5500 (T = 4.0) only K ≤ 5 is identifiable.
B. Random amplitudes (Gaussian-process / GLS; f2_gp.py, f3_robust.py): line coefficients N(0, (τ_p s_ab)²), s_ab = (γ₁γ₂/γ_aγ_b)^α (g₁₂/(γ_b − γ_a))^β, singles N(0, (τ_s/(γ_k²|ζ′(ρ_k)|))²), white nugget; c, C₃, C₄ flat; hyperparameters by REML; data interpolated to a uniform u-grid (step 0.005).

| model (C₂ unless Q) | window | K | shape | c ± s.e. | check |
|---|---|---|---|---|---|
| + C₃/L | [100,640] | 50 | B | −0.029 ± 0.021 | out-of-sample fit ≤ 400 → (400,640]: z = +0.54 |
| + C₃/L | [100,4000] | 50 | A | −0.0185 ± 0.0009 | **refuted**: z = +13.1 (fit ≤ 1280), window spread 0.04 |
| + C₃/L | [100,4000] | 10 / 20 / 50 / 100 | B | −0.0196 ± 0.0043 / −0.0207 ± 0.0059 / −0.0203 ± 0.0069 / −0.0204 ± 0.0072 | stable in K (sensitivity only) |
| + C₃/L | [100,5500] | 50 | B | −0.0194 ± 0.0061 | z = +0.73 (fit ≤ 1280), −0.47 (≤ 2000), +0.48 (≤ 4000 → (4000,5500]) |
| + C₃/L | [100,5500] | 50 | β-profile, REML-weighted (sensitivity mixture, not a posterior) | **−0.0192 ± 0.0035** | β ∈ [0.5, 1.5]: c ∈ [−0.0215, −0.0189], s.e. 0.0014–0.036; REML peaks at β = 0.8 |
| no trend | [100,5500] | 50 | B | −0.0026 ± 0.0005 | ML: worse than + C₃/L by 2ΔlogL = 7.6 |
| + C₃/L + C₄/L² | [100,5500] | 50 | B / fitted | −0.008 ± 0.044 / −0.007 ± 0.024 | ML: 2ΔlogL = 0.07 over + C₃/L (C₄ = 0.43 ± 1.7) |
| Q, no trend | [100,5500] | 50 | B | −0.04229 ± 0.0005 ⇒ C₂ ≤ −0.0039 ± 0.0005 | windows drift −0.0399 → −0.0440 ([100,640] → [640,5500]) |
| Q + C₃/L | [100,5500] | 50 | B / fitted | −0.0537 ± 0.0060 / −0.0534 ± 0.0028 ⇒ C₂ ≤ −0.015 | ML 2ΔlogL = 3.7 over no trend |
| Q + C₃/L + C₄/L² | [100,5500] | 50 | B / fitted | −0.044 ± 0.044 / −0.042 ± 0.024 | 2ΔlogL = 0.05 |
Window stability (shape B, K = 50, C₂ + C₃/L): [100,640] −0.029 ± 0.021, [100,1280] −0.028 ± 0.014, [100,2000] −0.016 ± 0.009, [100,4000] −0.020 ± 0.007, [100,5500] −0.019 ± 0.006, [300,5500] −0.012 ± 0.013, [640,5500] −0.006 ± 0.022, [1000,5500] −0.012 ± 0.029: mutually consistent. Shape A's windows give −0.034 ± 0.0013 ([100,1280]) and +0.0098 ± 0.006 ([640,4000]), which is inconsistent.
Implied κ̄. In the C₃-only model, κ̄ = (J − Σ_∞) + D₀ − C₂ = 0.040–0.044, outside nb2nd's C₂ = 0 window [0.021, 0.027]. This needs κ(N) to keep rising beyond its window mean 0.030. The Q fit gives κ(N) ≈ κ̄ − 0.073/L, consistent with that.
Head extrapolation over K ∈ [40, 5500]: Σ_∞ − A/log K gives J − Σ_∞ = 0.0344 (A = 0.0345); adding a B/log² K term gives 0.0347; Σ_∞ − A K^{−0.35} gives 0.0372; a log log K divergence fits worse (rms 1.1e−4 vs 6e−5). So J − Σ_∞ ∈ [0.034, 0.0384], where 0.0384 = J − Σ_5500 is rigorous.

## Error analysis
- **Correlation.** The data are a smooth function of u, so every error is low-frequency. 20-N block bootstraps understate the error by 2–6×. The honest bootstrap needs blocks of about 0.6 in u, and then gives ≈ 0.006 for the K = 0 + C₃/L fit, the same as the shape-B Gaussian-process s.e.
- **Unmodelled lines.** Their leakage into ĉ is Σ amp·|Σ_N a_N e^{iωu_N}| (f1). For the OLS K = 5 fit it is 0.07: every line resolved in the model makes the estimator ring at the unresolved ones.
- **Out-of-sample checks.** These overlapping holdouts test the fitted curve; they cannot validate the coverage of a constant extrapolated to L = ∞. Shape B passes:
  - 4 out-of-sample splits, |z| ≤ 0.76 for the window mean;
  - 8 windows, all consistent;
  - extension from N ≤ 4000 to N ≤ 5500, which was predicted at z = +0.48.
  Shape A fails (z = 13). The fitted gap exponent β ≈ 0.8 gives about half of shape B's s.e. (β = 1); that difference is the shape uncertainty and is carried by the β mixture.
- **Model truncation, the dominant error.** The constant is ĉ ≈ window mean − C₃⟨1/L⟩ − C₄⟨1/L²⟩ with ⟨1/L⟩ ≈ 0.16, ⟨1/L²⟩ ≈ 0.025.
  - Moving from no trend to C₃/L shifts c by −0.017, and the data support that shift at 2ΔlogL = 7.6.
  - Admitting C₄ is unsupported (2ΔlogL = 0.07) and unrefuted: C₄ = 0.43 ± 1.7, s.e.(c) ×7.
  - With nothing known about κ(N)'s 1/L expansion, the honest error at N ≤ 5500 is ±0.024 (fitted shape) to ±0.044 (shape B).
- **Numerics.** Double precision, cond ≤ 1.1e8; quadratic-form errors ≲ 1e−9, against effects ≥ 1e−4.

## Review
A statistics reviewer was given the regression model and the error bars inline (verdict only). It returned **GAP** on 13 items. All the gaps concern wording and the conditionality of the claims, none the verdict, and all are fixed above:
- Confirmed OK:
  - no model-robust 3σ (1);
  - REML for shapes, ML for trend orders (4);
  - the C₃-only result is conditional (6);
  - the Q inequality (8);
  - the projection, as a design calculation (10).
- Fixed:
  - the Gaussian-process s.e. are conditional on the covariance, the trend basis and the hyperparameters (2);
  - the holdouts test the fitted curve, not the coverage of the extrapolated constant, and stability in K is a sensitivity check (3);
  - the β mixture is not a posterior (5);
  - "C₄ untestable" became "C₄ unconstrained" (7);
  - Q̄ in the head bound is fitted (9);
  - a threshold for |C₂| ≈ 0.007 was added (11);
  - the Slepian count is a heuristic (12);
  - "every model gives C₂ < 0" became "every full-window estimate is negative", since the refuted shape A gave +0.0098 ± 0.006 on [640,4000] (13).

## Proposed NOGO line (replaces nb2nd's)
"nb2nd + nb2fals (Orca, N5 rate): C₂ = (J − Σ_∞) + D₀ − κ̄ (finite-N identity). Per-N data, N ≤ 5500, with the zero terms regressed out: free amplitudes are numerically singular beyond K ≈ 5 zeros, while random amplitudes are stable from K = 50 to 100. They give C₂ < 0, at ≥ 3σ only if the 1/L expansion stops at C₃/L (−0.019 ± 0.0035); with C₄/L², −0.007 ± 0.024. The sign needs N ≈ 10⁵–10⁷ or the 1/L expansion of κ(N) [NB2Fals](Riemann/Deep/NB2Fals/REPORT.md)"

## Files
Riemann/Deep/NB2Fals/. Run with uv run --no-project --with numpy --with scipy [--with mpmath], single-threaded, under .scratch/memcap.sh 290. Intermediates go to .scratch/nb2fals/ (perN{N}.npz, zeros100.npz).
| script | does | peak |
|---|---|---|
| code/f0_perN.py | dense Gram for N ≤ 4000 → per-N C₂, κ, J − Σ_N, R, d², cond | 204 MB (N = 4000; maxrss 219 MB) |
| code/f0p_packed.py | the same with packed storage, N ≤ 5500 | 227 MB (maxrss 243 MB) |
| code/e0_look.py | first look: step noise, periodograms, zeros ≤ 100 | 25 MB |
| code/f1_fits.py | free-amplitude OLS fits, K = 0…50, block bootstrap, leakage | 173 MB |
| code/f2_gp.py | random-amplitude GLS: main fits, windows, out-of-sample, projection | 115 MB |
| code/f3_robust.py | fitted shape, K = 100, 1/L-modulated lines, ML trend tests, C₄/L², β profile | 135 MB |
| code/f5_head.py | Σ_K extrapolations; the one-sided Q route | 26 MB |
| code/f6_proj.py | data-free s.e. projections, C₃ and C₃ + C₄ | not sampled (< 3 s) |
| code/f7_ols_oos.py | out-of-sample tests of the identifiable OLS fits | not sampled (< 3 s) |
| code/f4_export.py | writes the CSV files in data/ | not sampled (< 3 s) |
- logs/: the outputs of each script (the _5500 suffix marks the N ≤ 5500 runs).
- logs/f2_checks_A.log holds the shape-A half of a run whose shape-B half the watchdog killed at 355 MB: freed blocks of the 0.01-step projection were not returned by macOS malloc. The rerun used step 0.02 (f2_checks_B*.log, peak 111 MB).
- data/perN_640.csv, data/perN_5500.csv.gz: the per-N table.

## Caveats
- The Gaussian-process s.e. are conditional on:
  - lines from zeros ≤ 50 only (≤ 100 in the sensitivity run);
  - the shape family (γ_aγ_b)^−α(γ_b − γ_a)^−β;
  - a white nugget;
  - the trend basis;
  - REML point estimates of the hyperparameters.
  Only β is profiled, and the remaining hyperparameter uncertainty is not propagated.
- At N < 200 the data are interpolated linearly in u (the integer spacing 1/N exceeds the grid step 0.005); the step noise there (3e−4) is absorbed by the nugget.
- D₀ = −0.0131 and the C₂ = 0 window [0.021, 0.027] are nb2nd's. J − Σ_∞ ≥ 0.034 is an extrapolation; only J − Σ_∞ ≤ 0.0384 is rigorous.
- The projections assume the N ≤ 4000 hyperparameters and the trend basis hold up to the target N.
- No new theorem; no Lean.

RH relevance: N5
