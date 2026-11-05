FLAW: Lemma 4.0, Part 2 (p.19; Rev.3.1 p.20): the turnover-lobe identification r*(z_m)=m for the two outermost critical roots ("z~_{m-1}<z_m<z~_m by Chebyshev–Markov–Stieltjes interlacing"), which Theorem 4.1 uses to close the edge roots at every inductive step, is false; witness (d,k)=(9,0): z_8=10.6458329, y*_turn(z_8)=5.4267842 < y_7=5.7444601 (r*=7, not 8); z_7=9.7810145, y*_turn(z_7)=4.1471169 < y_6=4.7785316 (r*=6, not 7).

RH relevance: N6 (Jensen-polynomial hyperbolicity criterion; adversarial audit of a claimed proof)

# Audit of Zenodo 23977400 (J. F. Chung, "Hyperbolicity of Jensen Polynomials for the Riemann Ξ-Function …", Rev.3.0, 5 Nov 2026)

## 1. Short answer

- **First false step:** Lemma 4.0 (pp.18–19). Part 2 claims the turnover point y*_turn(z_m) of the outermost root lies in Hermite lobe m (m = d−1, d−2). It lies one lobe lower at every pair of the first inductive anti-diagonal N = d+k = 9 ((9,0),(8,1),(7,2),(6,3),(5,4),(4,5)) and at every (d,0), 6 ≤ d ≤ 30. Part 1's "L_j decreases strictly for j > r*" (peak lobe = r*) fails at the same witness: L_7 = 1.300755e6 < L_8 = 1.350514e6 with r* = 7.
- **Consequence:** Theorem 4.1 Part 2 (p.20) starts "At z_{d−1}, Lemma 4.0(2) establishes r*(z_{d−1}) = d−1, so the turnover point … lies inside the (d−1)-st Hermite interval" and builds the lobe-dominance chain on it. The premise is false, so the edge-root closure is unproved. Its conclusion (Ω_1 < … < Ω_{d−1}) is still true for k = 0, d ≤ 30 (check C4), so this part may be repairable.
- **First load-bearing failures (interior roots, the part that carries RH):** the uniform bounds that reduce all (d,k,m) to one anchor computation are false at actual roots:
  - Theorem 5.3 Part 2 (p.28): s_m ≤ μ_{d+k}K^(3) with K^(3) < 1.8165; false at (23,0), m = 20: s/μ = 1.83601.
  - Theorem 5.3 Part 3 (pp.28–30): disk η|u_s|² ≤ 0.7945; false at (22,0), m = 19 (0.79976), rising to 0.936 at (90,0), m = 87. The printed chain itself only gives 1−(1−λ)(1−0.45λ) with λ ≤ (1−4.38167(8d−2)^{−2/3})² → 1; the stated "≤ 0.788" is evaluated at d = 13 only.
  - Lemma A.4 Part 2a (pp.55–56): self-map T(Q) ⊆ Q false at (22,0), m = 19 and (30,0), m = 27 (82% and 87% of sampled points of Q leave Q).
  - Lemma 6.2b (p.34) / Lemma 6.4 (p.40): "locking τ_s ≥ 0.420 rad, ΔΘ_0(u_s) ≥ 0.4388 rad … c_den ≥ 0.306"; false at (4,5), m = 1 (τ_s = 0.2933, ΔΘ_0 = 0.3166), where the saddle is the unique fixed point inside the paper's own Q. Zone-I floors (τ_s ≥ 0.486692, ΔΘ_0 ≥ 0.4583) fail at (8,1), m = 1. With the true ΔΘ_0 = 0.2918 at (4,8), m = 1, the paper's own envelope gives M_total ≈ 1.04 ≥ 1.
- **Lean suite does not close it:** `Main.lean` `inductive_step_concrete_hermite_d_plus_k_ge_9` (l.4420) takes as hypothesis `h_bulk`, which contains `0 < (-1)^(d+(m+1)-1) * concreteSoninMarginA …`, i.e. the interior Obreschkoff sign itself, plus τ_s ∈ [0.486692, 0.486693] for every root and the remainder bound. Its docstring lists this as "Residual Obligations (b)". `bulkSaddleRealPart` is defined as 0.470776 + clamp(σ_offset); the "Olver–Volterra contraction map" is an affine map with that fixed point by construction, unrelated to the integrand.
- **Not flaws:** the frequency floor ω ≥ 12.2757 (Thm 5.1–5.2) holds (Samuelson + induction hypothesis; check C1 on 822 pairs). The true critical values have the right signs (the polynomials are hyperbolic for these small d). So this refutes the proof, not RH.

## 2. Logical chain (statement → what it feeds)

1. Lemma 1.1 moments (γ_n = (2n−1)!!^{-1}∫_0^∞u^{2n}Φ, with Φ(u) = 2Φ_dBN(u/2), Φ_dBN the de Bruijn–Newman kernel, so ∫_0^∞Φ(u)cos(zu)du = ξ(1/2+iz)/2) → Thm 1.2 (Pólya) → Lemma 2.1 (Obreschkoff window at the d−1 roots x_m = −z_m²/2 of J^{d−1,k+1}).
2. Lemma 3.1 (κ_k increasing), Thm 3.2: J^{d,k}(x_m) ∝ A_{m,k}, z_m f'(z_m) = B_{m,k}; step ⇔ A·B > 0 for all m.
3. Edge roots m ∈ {d−2, d−1}: **Lemma 4.0** (lobes L_j, turnover r*, M_r > 0) → Thm 4.1 (Ω-lobe dominance, tail annihilation).
4. Interior roots 1 ≤ m ≤ d−3: Thm 5.1–5.2 (floor ω ≥ 12.2757) + Thm 5.3 (ceiling, disk |z| ≤ 0.7945) → Lemma A.4 (cube K, Q, Banach saddle) → Thm 6.1 (split contour, exact) → Thm 6.2 (anchor constants at (9,0), ω = 12.272) → Lemma 6.2b (phase lock, c_den ≥ 0.306 via ΔΘ_0 ≥ 0.4388) → Thm 6.3 / Lemma 6.4 (M_total < 1) → A·B > 0.
5. Thm 7.1: induction along anti-diagonals from d ≤ 3 [CNV, Dimitrov–Lucas] and the 15 base pairs d+k ≤ 8.

Every constant in step 4 is computed once at the anchor (d,k,m) = (9,0,1), ω = 12.272, and extended to all (d,k,m) by asserted monotonicity; the failures above are exactly those extensions.

## 3. Checks (all in `Riemann/Deep/ChungAudit/`, each prints PASS/FAIL; PASS = paper's claim holds there)

Main driver (8 s): `uv run --no-project --with mpmath --with numpy python3 audit_checks.py`

```
C1 Thm 5.1-5.2 frequency floor omega_1>=Omega_LS>=12.2757 on [(9,0),(8,1),(7,2),(6,3),(5,4),(4,5),(15,0),(30,0),(10,10)] : PASS
C2 Lemma 4.0(2) r*(z_{d-1})=d-1 at (9,0): z_8=10.6458328564 y*=5.4267842, need (y_7,y_8)=(5.7444601,6.8891224): FAIL
C2 Lemma 4.0(2) r*(z_{d-2})=d-2 at (9,0): z_7=9.78101453303 y*=4.1471169, need (y_6,y_7)=(4.7785316,5.7444601): FAIL
C3 Lemma 4.0(1) peak of L_j at r* at (9,0), z_8: r*=7, L_7=1300754.966, L_8=1350513.744, L_9=606895.724: FAIL
C4 Thm 4.1(2) conclusion Omega increasing at (9,0) m=8 (8 lobes): PASS
C5 Thm 5.3(2) ceiling s_{d-3}/mu_{d+k} < 1.8165 at (22,0), m=19: s=77.542788 mu=42.781995 ratio=1.81251: PASS
C5 Thm 5.3(2) ceiling s_{d-3}/mu_{d+k} < 1.8165 at (23,0), m=20: s=79.902282 mu=43.519524 ratio=1.83601: FAIL
C6 Thm 5.3(3) disk eta|u_s|^2<=0.7945 at (21,0,18): omega=55.8413 u_s=1.29366+0.27837i |F|=6e-15 alpha=0.79259: PASS
C6 Thm 5.3(3) disk eta|u_s|^2<=0.7945 at (22,0,19): omega=58.0785 u_s=1.31112+0.27487i |F|=1e-14 alpha=0.79976: FAIL
C6 Thm 5.3(3) disk eta|u_s|^2<=0.7945 at (30,0,27): omega=75.3386 u_s=1.42724+0.25234i |F|=1e-14 alpha=0.84198: FAIL
C7 Lemma 6.4 floors at (9,0,1): omega=12.8297 sigma=0.5012 tau=0.4885 dTheta0=0.4720: PASS
C7 Lemma 6.4 floors at (8,1,1): omega=13.7020 sigma=0.5375 tau=0.4484 dTheta0=0.4416: FAIL
C7 Lemma 6.4 floors at (4,5,1): omega=16.5665 sigma=0.7049 tau=0.2933 dTheta0=0.3166: FAIL
C7 Lemma 6.4 floors at (4,8,1): omega=19.6493 sigma=0.7812 tau=0.2598 dTheta0=0.2918: FAIL
C8 App. D anchor is a root of printed F (A.2): F=-0.1021-0.0532i: FAIL
C9 Thm 6.3 consequence B/A >= (1-0.853)Lambda_pos at (9,0), m=1..6: PASS
```

Supporting scripts (scratch paths are regenerable; data/gammas_dps70.json is kept):

| Script | What it checks | Result |
|---|---|---|
| gammas.py | γ_0..γ_100 by quadrature (dps 70); μ_2..μ_7 vs Lemma 1.1 | PASS (match to 1e-8) |
| crosscheck.py | γ_n from Taylor coefficients of ξ(1/2+w) vs quadrature; J^{8,1} roots by sign-change root finding vs polyroots; He_17 zeros two ways | PASS (1e-49, 2e-29, 7e-15) |
| roots_floor.py | roots of J^{d−1,k+1}, 4 ≤ d ≤ 30, k ≤ 30, d+k ≥ 9 (822 pairs); Thm 5.1–5.2 floor | PASS (min ω_1 = 12.830 at (9,0)) |
| turnover_k0.py | Lemma 4.0(2) at k = 0, d = 9..30 | FAIL for every d (r*(z_{d−1}) = d−2) |
| turnover_k.py | Lemma 4.0(2) for k ≥ 1: (8,1),(7,2),(6,3),(5,4),(4,5),(9,1),(8,2),(10,2),(6,6) | FAIL for all 9 pairs |
| lobes_k0.py / lobes_k0_small.py | Lemma 4.0(1) lobe unimodality with peak at r*; M_r > 0 | peak-at-r* FAIL at (3,0),(4,0),(5,0) m=1; (6,0) m=5; (7,0) m=3,6; (8,0) m=3,7; (9,0) m=8; (10,0) m=6,7,9; (11,0) m=5,8,10. M_r > 0 PASS everywhere tested (d ≤ 14, k = 0) |
| omega_k0.py (mpmath) / omega_np.py (numpy) | Thm 4.1(2) conclusion Ω increasing at m = d−1, k = 0 | PASS (d = 4..30; both methods agree at (9,0)) |
| saddle_scan.py | printed saddle equation (A.2/A.4) solved at all 11,683 interior roots of the grid by continuation from the anchor (residual ≤ 1e−13); 13 uniform claims | FAIL: τ_s ≥ 0.420 (9205 roots), ΔΘ_0 ≥ 0.4388 (4973), Φ_a ≤ 1.2604 (112), α ≤ 0.7945 (183); PASS: σ_s ≥ 0.4708, A_core ≥ 6.10, A_eff ≥ 6.40, κ_Z ≥ 0.72, β ≤ 1.10, Re S ≥ 0.420, |Z| ≥ 13.83; |V(u_s)| ≥ 0.22 fails only at (4,29),(4,30), Regime II, where the paper uses Olver's form (not counted) |
| disk_trend.py | disk α at m = d−3, k = 0, d up to 90 | FAIL for d ≥ 22 (α = 0.818, 0.842, 0.874, 0.909, 0.936 at d = 25, 30, 40, 60, 90) |
| q_selfmap.py | Lemma A.4(2a): saddle ∈ Q and T(Q) ⊆ Q | PASS at (4,5,1),(8,1,1); FAIL at (22,0,19),(30,0,27): saddle outside Q, T sends 82% / 87% of sampled Q-points out of Q |
| ba_ratio.py | exact B/A = 2x_m J''_{d,k}(x_m)/(d J_{d,k}(x_m)) vs (1−0.853)Λ_pos(u_s) | PASS (no contradiction; (B/A)/Λ_pos ∈ [0.961, 1.220] on the grid) |

## 4. Witness details

**Lemma 4.0 at (d,k) = (9,0).** J^{8,1}(X) = Σ_j C(8,j)γ_{1+j}X^j has roots x_m = −z_m²/2 with z = 4.3372171, 5.6384398, 6.5994351, 7.4479891, 8.2348860, 9.0001786, 9.7810145, 10.6458329 (z_8 checked: f_0(z_8)/∫|integrand| = 2e−30). Positive zeros of He_17: 0.7518, 1.5099, 2.2810, 3.0738, 3.9001, 4.7785, 5.7445, 6.8891. The turnover y* solves 1/y² + 1/2 = (2/z²)κ_0(√2y/z) (paper's own formula, p.18), confirmed as the argmax of Λ~(y) = yΦ(√2y/z)e^{y²/4}. y*(z_8) = 5.4268 ∈ (y_6, y_7) → lobe 7; y*(z_7) = 4.1471 ∈ (y_5, y_6) → lobe 6. The paper requires lobes 8 and 7. Full list r*(z_m), m = 1..8: [2,2,2,3,4,5,6,7]. In u-space (as used on p.20): u*_turn(z_8) = √2·y*/z_8 = 0.72091, nodes ξ_6 = 0.63480, ξ_7 = 0.76312, ξ_8 = 0.91517; so u* ∈ (ξ_6, ξ_7), not (ξ_7, ξ_8) as Theorem 4.1(2) states, and the peak ξ_7 of arch r = d−2 = 7 lies above u*, contradicting "every pre-peak arch 1 ≤ r ≤ d−2 lies below u*".

**Theorem 5.3(2) at (23,0), m = 20.** s_20 = z_20²/2 = 79.902282, μ_23 = γ_22/γ_23 = 43.519524, s/μ = 1.83601 > 1.8165. The ratio grows: 1.879 (25,0), 1.971 (30,0). Root data only, no saddle model. The paper's own Table A.3c ceiling ω_ceil(25) = 62.444 is below the actual ω_22(25,0) = 64.674.

**Lemma 6.4 at (4,5), m = 1.** ω = 16.5665 (actual root of J^{3,6}); the printed saddle equation F(u) = W − 13/2 − k/u − γu − iωS(u) − E(u) = 0 has root u_s = 0.7049 + 0.2933i with W_s = 21.429 + 14.245i ∈ Q, and T maps all 426 sampled points of Q into Q. So u_s is the paper's own unique saddle. ΔΘ_0 = arg κ_0(u_s) = 0.3166 < 0.4388 and τ_s = 0.2933 < 0.420.

**Anchor (App. D).** The certified anchor u = 0.4707764 + 0.4866925i is not a root of the printed equation: F = −0.1021 − 0.0532i. The actual root at ω = 12.272 is 0.47413 + 0.48433i (|Z| = 14.08, Re Z = 6.508, A_eff = 7.035, vs printed 13.830, 6.1086, 6.910). The paper's numbers come from an equation it does not print.

## 5. Page-image confirmations (pdftoppm, 100–110 dpi)

- p.18: Lemma 4.0(1): L_j(z) = (2/z²)∫Λ~_k|ψ_{2d−1}|dy, Λ~_k = yΦ_k(√2y/z)e^{y²/4}; "(1/y)∂_y ln Λ~ = 1/y² + 1/2 − (2/z²)κ_k(√2y/z)"; "proving 0 < L_1 < ⋯ < L_{r*(z)}"; "ensuring that L_j(z) decreases strictly for j > r*(z)".
- p.19: "each ordered root z_m of f_k(z) = 0 (m ∈ {d−2, d−1}) satisfies z~_{m−1} < z_m < z~_m, so y_{d−2} < y*_turn(z_{d−1}) < y_{d−1} (r*(z_{d−1}) = d−1) and y_{d−3} < y*_turn(z_{d−2}) < y_{d−2} (r*(z_{d−2}) = d−2)."
- p.20: Theorem 4.1(2): "At z_{d−1}, Lemma 4.0(2) establishes r*(z_{d−1}) = d−1, so the turnover point u*_turn,k(z_{d−1}) ∈ (ξ_{d−2}, ξ_{d−1}) lies inside the (d−1)-st Hermite interval. On (0, ξ_{d−1}], every pre-peak arch 1 ≤ r ≤ d−2 lies below u*_turn,k(z_{d−1}) …" → "0 < Ω_1 < Ω_2 < ⋯ < Ω_{d−1}". (Same page also uses κ_k'' > 0, which Lemma 3.1 does not prove; numerically true for k = 0.)
- p.28: "combined with r*(z_{d−3}) = d−3 … K^(3)_{d,k} := min(K^(3,exact), y²_{d−3}m_{d,k}/((2d−1)(u*)²)) < 1.8165 bounds s_m ≤ μ_{d+k}K^(3)_{d,k} uniformly across all d ≥ 4, k ≥ 0, 1 ≤ m ≤ d−3." (r*(z_{d−3}) = d−3 is not even claimed by Lemma 4.0; actual r* = d−4.)
- p.40: "Consequently, ω_{m,k} dominates k/|u_s| for every k ≥ 0, locking τ_s ≥ 0.420 rad, ΔΘ_0(u_s) ≥ 0.4388 rad, … c_den(d,k) ≥ 0.306 across all interior critical roots in K"; "every saddle in Zone I … τ_s ≥ τ_min = 0.486692 …, ΔΘ_0(u_s) ≥ 0.4583 rad".

## 6. Self-review

1. Definitions: γ_n, μ_n, J^{d,k}, z_m, ω_{m,k} = z_m√(4d+2k−1)/2, s_m = z_m²/2 and the saddle equation are the paper's. γ_1..γ_3 and μ_2..μ_7 reproduce Lemma 1.1 digit for digit; κ_0(0) = 18.726905 matches §1.1. L_j in u-form equals the paper's y-form (substitution y = zu/√2). Hermite zeros are probabilists' He_{2d−1}, lobes j = 1..d with y_0 = 0, y_d = ∞, as on p.18.
2. Is the failing statement used? Lemma 4.0(2) is cited verbatim by Thm 4.1 Part 2 (p.20) and, as r*(z_{d−3}) = d−3, by Thm 5.3 Part 2 (p.28). The K^(3), disk and ΔΘ_0 bounds feed Lemma A.4 → Lemma 6.2b → Lemma 6.4 → M_total < 1 (pp.37–41).
3. Second computation: γ by ξ-Taylor coefficients (agree to 1e−49); roots by sign changes (2e−29); He zeros by Golub–Welsch and numpy (7e−15); y* as root of the log-derivative and as argmax; L_7, L_8 by 40-point Gauss–Legendre and by mpmath tanh-sinh (agree to 10 digits); Ω by mpmath and numpy.
4. Caveat on the saddle-model findings (C6, C7, A.4): they solve the paper's printed equation, which its own anchor does not satisfy (C8). The (4,5,1) and (8,1,1) witnesses lie inside the paper's Q with T(Q) ⊆ Q on samples, so they are the paper's saddles by its own Lemma A.4. The Lemma 4.0 and Thm 5.3(2) witnesses use root data only and do not depend on the saddle model.
5. Part 1 wording: "L_j decreases strictly for j > r*" could be read as "decreasing among j > r*" (true here: L_8 > L_9). The Part 2 claim r*(z_{d−1}) = d−1 has no such ambiguity; it is the headline.

## 7. Other observations (not needed for the verdict)

- §1.1 / Lean `auditedStiffnessDeriv_pos_global`: the "audited lower majorant" κ'_{0,1} − 16a_2²e^{−3a_1} is not a lower bound of κ'_0 near u = 0.05 (15.73 vs κ'_0(0.05) = 2.377). κ'_0 > 0 itself is fine (Taylor on (0,0.2], n = 1 mode beyond).
- Table A.3a lists ω = 16.030 as "m = 6, d = 9, k = 0"; the actual ω_6(9,0) = 26.623.
- B/A at actual roots tracks Λ_pos(u_s) closely (ratio 0.961–1.220 on the grid), so the remainders R/A are small in reality; the proof fails in its claimed uniform bounds, not in the shape of the split-contour step.

## 8. Novelty and versions

- `rg -n '23977400' REFUTATIONS.md NOGO.md`: no entry. R345/R346 (Riemann/Tree) recorded "flaw not yet located".
- Zenodo versions (concept 23405056): 9 versions 15 Oct–5 Nov 2026; newest is 23987495 (Rev.3.1, 5 Nov). Rev.3.1 keeps the same statements: Lemma 4.0(2) CMS claim (p.20), K^(3) < 1.8165 (p.29), disk ≤ 0.788 (p.31), ΔΘ_0 ∈ [0.4388, 1.5708] (p.35), "locking τ_s ≥ 0.420" (p.43). The flaw persists. No public comments on the record were found (none sought beyond the API).

## 9. Suggested REFUTATIONS.md line (≤ 40 words)

Zenodo 23977400 (Chung, RH via Jensen hyperbolicity): Lemma 4.0(2) false — at (d,k)=(9,0) outermost turnover y*=5.4268<y_7=5.7445 (r*=7≠8), so Theorem 4.1's edge closure is unproved; interior bounds also fail (Thm 5.3 ceiling 1.836>1.8165 at (23,0); ΔΘ₀=0.317<0.4388 at (4,5)).

## 10. Reproduce

```
cd Riemann/Deep/ChungAudit
uv run --no-project --with mpmath --with numpy python3 audit_checks.py          # 8 s, all headline checks
uv run --no-project --with mpmath python3 turnover_k.py                         # 4 min, Lemma 4.0(2) for k>=1
uv run --no-project --with mpmath --with numpy python3 crosscheck.py            # 1 min, independent recomputation
# full grid (regenerates scratch files):
uv run --no-project --with mpmath python3 gammas.py 100 70 data/gammas_dps70.json            # 4 min
uv run --no-project --with mpmath python3 roots_floor.py data/gammas_dps70.json 30 30 /tmp/roots30.json   # 11 min
python3 saddle_scan.py /tmp/roots30.json /tmp/saddles30.json                                  # 3 min
uv run --no-project --with mpmath python3 turnover_k0.py /tmp/roots30.json                   # 1 min
```
