RESULT: 2611.03924's dip is real, not an artefact: certified d_640² log 640 = 0.0453405057 < C = 0.0461914179 (ball width < 1e-23), but d_N² log N − C changes sign 15 times for N ≤ 1280 (first below C at N = 41), and for Bettin–Conrey–Farmer's ramp the second-order term is (J + D(N))/log² N with J = 3.5626634565 explicit (RH + BCF's hypothesis), D(N) almost periodic with mean ≈ −0.014.

# nbrate: rigorous d_N and the second-order term of the Nyman–Beurling/Báez-Duarte distance (node N5)

## Definitions (fixed by the TeX)
- V1 (Báez-Duarte, BBLS, Burnol, BCF): d_N² = min_c ‖χ_(0,1] − Σ_{k≤N} c_k {1/(kx)}‖² in **L²(0,∞)** = 1 − bᵀG⁻¹b, G_jk = ⟨ρ_{1/j}, ρ_{1/k}⟩ (Vasyunin: for coprime h,k
  g(h,k) = [(log 2π − γ)(h+k) + (k−h) log(h/k) − π(V(h,k)+V(k,h))]/(2hk), V(h,k) = Σ_{m<k} {mh/k} cot(πm/k), G_{dh,dk} = g(h,k)/d), b_k = (log k + 1 − γ)/k, ‖χ‖ = 1.
  (The brief's "L²(0,1), unconstrained" is not this problem: on x > 1 the span is (Σc_k/k)/x.)
- V2 = 2611.03924's d_n(1/2)² (their §subsec:dn): same with the constraint Σ c_k/k = 0 (span supported on (0,1]); V2 ≥ V1.
- Checks: our Vasyunin-based constrained Gram equals 2611.03924's elementary series (1/jk)Σ_m (m mod j)(m mod k)/(m(m+1)) (digamma over a period) to 15 digits at
  (2,2), (2,3), (3,5), (4,6), (5,7), (6,10) (logs/s1_double.log); G₁₁ = log 2π − γ, G₁₂, G₂₃ match explore13's direct integration.

## (A) Certified values — method
1. Cotangent sums exactly: cot(πm/k) as Arb balls (220 bits) → integers C_m = round(2^140 cot) ; V(a,k) = (1/k)Σ_{m≤(k−1)/2}(2(ma mod k) − k)C_m·2^−140 summed EXACTLY
   (26-bit limbs, integer-valued float64 products < 2^47); error ≤ (k/2)·2^−140 per sum.
2. g(h,k) in Arb, Gram entries as fixed-point integers Gt = round(2^120 G), |G − Gt·2^−120| ≤ ε_G = 7.5e−37 (all 1280² entries); b likewise.
3. Certificate (any c): d_N² = Q_N(c) − rᵀG_N⁻¹r, r = G_N c − b_N, Q_N(c) = ‖χ − Σc_kρ_{1/k}‖² ⇒ d_N² ∈ [Q − ‖r‖²/λ_min, Q]. c = double Cholesky solution rounded to
   2^−60; G_N c, bᵀc, cᵀGc computed EXACTLY in integers (20-bit limbs) plus explicit ε_G, ε_b propagation; Q evaluated as an exact rational.
4. λ_min(G_N) ≥ λ_min(G_1280) (interlacing) ≥ μ − ‖E‖_F where G − μI = L̂L̂ᵀ + E, L̂ = 50-bit dyadic rounding of a double Cholesky of G − μI, L̂L̂ᵀ exact (17-bit limbs),
   floating error of the subtraction bounded by 12u(|G| + μ + ‖ŷ_i‖‖ŷ_j‖) per entry: **λ_min(G_1280) ≥ 4.768·10⁻⁷** (μ = 2⁻²¹, ‖E‖_F ≤ 5.5·10⁻¹⁴; estimate 8.998·10⁻⁷).
5. V2: c = (c₁, c′) with c′ dyadic and c₁ = −Σc′_k/k an exact rational (≈1850-bit denominator); residual r′_k = r_k − r₁/k; λ_min(TᵀG T) ≥ λ_min(G) (TᵀT ⪰ I).
Resulting ball widths: ≤ 6.4·10⁻²⁴ for every N ≤ 1280 (V1) and ≤ 8.7·10⁻²⁴ (V2); the residual term ‖r‖²/λ dominates. Condition numbers (double eigvalsh, λmax/λmin):
5.7e3 (40), 2.4e4 (80), 9.9e4 (160), 3.9e5 (320), 1.6e6 (640), 6.2e6 (1280) — mild: IEEE double loses ≤ 7 digits, so double Cholesky is adequate here.

## (A) Certified table (V1 unconstrained; C = 2 + γ − log 4π = 0.04619141793)
| N | d_N² (certified; display rounding ≤ 4·10⁻¹⁶, true width < 10⁻²³) | d_N² log N | d_N² log N − C | cond(G_N) |
|---|---|---|---|---|
| 40 | 0.012715197195337 | 0.0469048297 | +7.134e−4 | 5.7e3 |
| 41 | 0.012326491948694 | 0.0457753162 | −4.161e−4 | |
| 48 | 0.011873207147379 | 0.0459635715 | −2.279e−4 | |
| 49 | 0.011872614669304 | 0.0462060828 | +1.467e−5 | |
| 80 | 0.010983481552538 | 0.0481299087 | +1.939e−3 | 2.4e4 |
| 160 | 0.0093455910085642 | 0.0474304988 | +1.239e−3 | 9.9e4 |
| 205 | 0.0086793734308363 | 0.0462003914 | +8.97e−6 | |
| 206 | 0.0086522007981348 | 0.0460978544 | −9.36e−5 | |
| 320 | 0.0079537374311735 | 0.0458797079 | −3.117e−4 | 3.9e5 |
| 411 | 0.0076744466976446 | 0.0461893728 | −2.045e−6 | |
| 524 | 0.0073767301615984 | 0.0461893345 | −2.083e−6 | |
| 583 | 0.0072572423228183 | 0.0462154779 | +2.406e−5 | |
| 584 | 0.0072490880608171 | 0.0461759732 | −1.545e−5 | |
| **640** | **0.0070170593529348** | **0.0453405057** | **−8.509e−4** | 1.6e6 |
| 695 | 0.0067526049812340 | 0.0441884517 | −2.003e−3 (minimum) | |
| 1280 | 0.0063611026375655 | 0.0455112426 | −6.802e−4 | 6.2e6 |
Certified sign of d_N² log N − C for EVERY 2 ≤ N ≤ 1280 (no undecided N; tightest margin 2.0e−6 at N = 411, 524): + [2,40], − [41,48], + [49,205], − [206,212],
+ [213,219], − [220,282], + [283,286], − [287,371], + [372,380], − [381,382], + [383,399], − [400,411], + [412,417], − [418,524], + [525,583], − [584,1280]
(981 of 1279 values below C). **First drop below C: N = 41.**

V2 (2611.03924's constrained d_n(1/2)²), certified vs their Table tab:dn:
| n | V2 certified | d² log n − C | 2611.03924 |
|---|---|---|---|
| 40 | 0.01273583539133112 | +7.895e−4 | 1.274e−2 |
| 80 | 0.01099879595104947 | +2.006e−3 | – |
| 160 | 0.009356604331758195 | +1.295e−3 | 9.357e−3 |
| 320 | 0.007961688783869451 | −2.658e−4 | – |
| 640 | 0.007023264380134277 | −8.108e−4 | 7.023e−3 |
| 1280 | 0.006366197363433226 | −6.437e−4 | 6.366e−3 |
Their n ≤ 1280 values (exact Gram data, double, Cholesky) are correctly rounded to all 4 printed digits: **no artefact**. Their n = 2560, 5120 values (quadrature, ≈5e−4
relative accuracy by their own check) are not certified here; they give d² log n − C = +3.6e−4 (2560, above C again, margin ≫ their error) and −1.07e−3 (5120).
So "convergence from below" is false on the accessible range: d_N² log N oscillates around C.

## (B) Second-order term
### BCF ramp (c_k = −μ(k)(1 − log k/log N), i.e. 1 − ζV_N) — derived (details: derivation.md)
Under RH, simple zeros, BCF's Σ_{|γ|≤T}|ζ′(ρ)|⁻² ≪ T^{3/2−δ}, and (H′) uniformity of the O(N^{−ε}) remainders as in BCF's proof (the pair-series
interchange is proved by domination, derivation.md):
  log² N · I_N = C log N + J + D(N) + O(N^{−ε}),  I_N = (1/2π)∫|1 − ζV_N|² dt/(¼ + t²),
  J = (1/2πi)∫_(c) (ζ′/ζ)(s)(ζ′/ζ)(1−s) ds/(s(1−s)) = π²/12 − γ² − 2γ₁ + (1+γ) log 2π − Σ_{n≥1}(ζ′/ζ)(2n+1)/(2n(2n+1)) = **3.5626634565363677**,
  D(N) = Σ_{ρ1,ρ2} K(ρ1,ρ2) N^{ρ1+ρ2−1}/(ζ′(ρ1)ζ′(ρ2)), K(ρ1,ρ2) = (1/2πi)∫_(c) ζ(s)ζ(1−s) ds/((ρ1−s)²(ρ2−1+s)² s(1−s)), 0 < c < ½.
Method: BCF's own decomposition, each piece evaluated: the ζ′/ζ·ζ′/ζ piece (BCF: "trivially O(1/log² N)") is exactly J (residues at 0 and −2n); the zero residues give exactly
C log N — their non-log parts are purely imaginary and cancel in conjugate pairs (BCF's displayed residue has −½ζ″/ζ′ where the functional equation gives −1; we checked
Re ζ″/ζ′(ρ) = χ′/χ(ρ) to 12 digits for ρ₁, ρ₂, ρ₃); the pair piece is D(N). So **C₂ does not exist as a constant unless the oscillating part vanishes; the mean coefficient is J + D₀**,
D₀ = Σ_ρ K(ρ, 1−ρ)/|ζ′(ρ)|². Numerically (s6b, 40 zeros, terms ≈ −0.28 log(γ/2π)/γ², density-tail extrapolation): D₀ = −0.0142 ± 0.003; the off-diagonal amplitudes
are ≤ 1.8e−4 (largest at frequency γ₂ − γ₁ = 6.887), so Ω is tiny; the visible scatter is the O(N^{−ε}) remainder.
Check against certified Gram values of I_N: mean of log² N·I_N − C log N = 3.5478 (N ∈ [800,1280], sd 0.008), 3.5496 ([400,800]), 3.5480 ([200,400]) vs **J + D₀ = 3.5485 ± 0.003**.
Sign: positive and large; the ramp approaches C from above and is far from asymptotic (I_1280 = 0.0756 = 12 d²_1280; C log N exceeds J only for N > e^77).
### Optimal coefficients — not derived; where it stops
- A fixed smooth profile c_k = −μ(k)P(log k/log N), P(0)=1, P(1)=0 has first-order constant C∫₀¹P′² (≥ C, equality only for the ramp) and a non-local term −P′(0)ζ′/ζ/log N,
  i.e. P′(0)²J at second order (heuristic local Fourier analysis). The optimum removes J with an N-dependent boundary layer: its profile −c_k/μ(k) is flat at small k
  (0.981, 0.969, 0.964, 0.903 at k = 1, 2, 3, 5 for N = 1280; ramp 1, 0.903, 0.846, 0.775), keeps weight ≈ 0.16–0.20 near k = N, and non-squarefree c_k up to 0.26.
  The optimal second-order term needs this boundary-layer/endpoint optimisation against the zero sums: not done.
- Numerical fit (certified V1 data): log N·(d_N² log N − C) has mean −0.0040 (sd 0.0044) on [100,1280], −0.0062 (sd 0.0009) on [800,1280], range [−0.0131, +0.0087];
  2611.03924's V2 gives +0.0028 (2560), −0.0091 (5120). Least squares d² = C/L + C₂/L² gives C₂ = −0.0012 … −0.0052 depending on the window; adding C₃/L³ gives
  C₂ = −0.024 … −0.044, C₃ = +0.13 … +0.25 (unstable). Honest summary: **C₂^opt ∈ [−0.01, +0.005], sign not determined**; the second-order coefficient of the optimum
  is ≈ 0 ± 0.01 and oscillates (≈ 700× smaller than the ramp's 3.55). Möbius × degree-10 polynomial + free c₁…c₁₂ already reaches +0.0077 at N = 1280.
- Consistency with (A): the optimal's second-order coefficient oscillates in sign, matching the certified sign changes; the ramp's is +3.55.

### Side finding: Veltri's numerics (Zenodo 22647730 v3, §"A numerical remark")
His table "c_N = λ·(1/2π)∫|1−ζV_N|²dμ" = 0.60431, 0.52511, 0.47279, 0.42631 (N = 100, 200, 400, 800) equals our λ·(I_N − (Σc_k/k)²) to all 5 digits
(results.txt, veltri_check), i.e. he computed the L²(0,1) norm, dropping the tail (Σc_k/k)² = (1 + O(N^{−1/2+ε}))/λ². The stated integral is λI_N = 0.82191, 0.71358, 0.63980, 0.57582
(direct critical-line quadrature: I_100 ≥ ∫₀^400 part = 0.17767 > his 0.13123; s5_direct.py). His K ≈ 2.5 is J − 1 + D; for the stated quantity K = J + D₀ ≈ 3.55.
This answers his open question (iv) ("identify the second-order term").

## (C) Novelty
- Read in TeX: BCF 1211.5191 (only O(1/log² N); residue typo above); 2611.03924 (double precision, no bounds, dip "compatible with the liminf");
  Veltri Zenodo 22647730 v3 (K ≈ 2.5 conjecture, open question (iv)); Zenodo 22650745 (averaged decoupling; no constant).
- Scans: lst.py '(Nyman|Beurling|B[aá]ez|Baez|Vasyunin|autocorr|fractional part)' over math.NT/CA/CV/SP/PR/FA/CO, Oct+Nov 2026 + Zenodo '"Nyman-Beurling" OR "Baez-Duarte"';
  axq.py 'abs:"Nyman-Beurling" OR abs:"Baez-Duarte"', 'abs:"Beurling" AND abs:"Riemann hypothesis" AND abs:"distance"'; HAL API ("Nyman-Beurling", "Báez-Duarte", "Baez-Duarte",
  Vasyunin): no second-order term, no rigorous d_N beyond small N. Zenodo 24054105 (8 Nov 2026): rational certificates for M ≤ 24 only; Zenodo 24088197 (9 Nov 2026): Arb Gram to
  N = 128 (its d²₁₂₈ = 0.009670 is the constrained V2 value, ours 0.0096703; V1 certified 0.0096585). BBLS, Adv. Math. 149 (2000) §8 (floating point, n ≤ 10⁴, fit 1/(21.99 log n) = 0.04548/log n < C/log n, as quoted
  by 2611.03924) was not available (math/0306251 is a different BBLS paper without numerics).
- New here (as far as found): (i) rigorous d_N² for all N ≤ 1280 and the certified sign pattern of d_N² log N − C; (ii) the explicit second-order constant J for the BCF ramp and the
  structure J + D(N) (conditional; proof sketch); (iii) Veltri's normalisation slip and BCF's residue typo. (i) is a record/falsification datum; (ii) answers a stated open question.

## Caveats
- The Theorem is a proof sketch under RH + simple zeros + BCF's hypothesis + (H′). Orchestrator reviewer: CORRECT; ROUTINE (explicit conditional refinement and erratum
  of BCF's decomposition → COROLLARIES; certified d_N table → RECORDS); H′ = remainder uniformity only (the pair-series interchange is proved by domination with
  S(t) = Σ_ρ|ζ′(ρ)|⁻¹/(η² + (t−γ)²), checked here). Earlier algebra review: residues, α₁, the zero residue and contour orientations CORRECT. D₀'s value uses a 40-zero sum plus an extrapolated tail (empirical law −0.28 log(γ/2π)/γ²).
- The rigorous part relies on Vasyunin's closed formula (literature theorem; checked here against an independent elementary series to 15 digits) and on standard IEEE
  error bounds (12u per entry) for the floating part of the λ_min certificate; all other arithmetic is exact integer arithmetic.
- 2611.03924's n > 1280 values are not certified (memory cap: a 2560 Gram does not fit in 300 MB).
- MEMORY CAP VIOLATED: s2_rig.py 1280 peaked at ~567 MB and s3_v2.py at ~462 MB (top), despite trimming; both finished (54 s, 33 s). Do not rerun at N = 1280 without
  further trimming (stream the Gram build; avoid eigvalsh).
- Under RH + BCF's hypothesis the excursions above C are ≤ (J + D(N))/log N + o(1/log N) (d_N ≤ ramp) and Burnol makes those below C o(1); no finite computation
  sees this; the data cannot discriminate RH (logarithmic rate).

RH relevance: N5 (rate of the Nyman–Beurling criterion)

## Files
Riemann/Deep/NBRate/: REPORT.md; derivation.md (second-order derivation + numerical check); results.txt (all run logs concatenated); code/ (scripts, run order
s1 → s2 → s3 → s4 → s5 → s6 → s6b; they read/write .scratch/nbrate/, run with uv run --no-project --with python-flint --with numpy --with scipy --with mpmath):
s1_double.py (double Gram, all-N d², V2, eigenvalues, ramp), s2_rig.py (rigorous V1 for all N ≤ 1280, λ_min certificate, ramp bounds), s3_v2.py (rigorous V2 at 6 n),
s4_rate.py (J to 20 digits, optimal profile, profile families), s5_direct.py (critical-line quadrature of I_N), s6_pairs.py (pair term, 6 zeros), s6b_D0.py (D₀, 40 zeros).
Working copies, rig_V1.npz (certified lo/hi for every N ≤ 1280, signs, ramp bounds) and STATE.md in .scratch/nbrate/ (gitignored).
