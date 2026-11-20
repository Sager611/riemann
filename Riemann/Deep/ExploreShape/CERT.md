# CERT: Λ > 3.367 (= 3367/1000), unconditionally (explore53, 20 Nov 2026; built up as Λ > 3.34 (order 4) → 3.358 (order 5) → 3.367 (order 6))

Λ = limsup (t_{n+1} − t_n)/(2π/log t_n) over consecutive ordinates of zeros on the critical line. Before: Λ > 3.324
(GapBar, Wirtinger step). Here BM's Wirtinger step is replaced by a certified 4th-order per-gap inequality; every
mean value is BM's (Lemma main), evaluated in Arb.

## Statement and proof outline
Let F(t) = e^{ivtL}ζ(½+it)ζ(½+it+iκπ/L)M(½+it+iκπ/(2L)), L = log(T/2π), M(s) = Σ_{h≤y} d_2(h)P[h]h^{−s}, y = T^ϑ,
r = 2, ϑ = 2499/10000, κ = 167/50, v = 21294/15625, P of degree 64 (coefficients A_n below).
1. **Per-gap inequality (V).** For every g ∈ (0, π] and every real f ∈ H⁴(0, g) with f(0) = f(g) = 0:
   Σ_{j≤4} a_j‖f^{(j)}‖² + φ12[f'f'']_0^g + φ32[f'''f'']_0^g ≥ 0, with
   a = (−5899179108889603/6134092019850000, 1, 14789027537/948575124080, −4037122085/176350902122,
   1982999731/771179571628) (`logs/cert4_k3.34.json`; a0 already contains the +2·10⁻⁵ safety shift),
   φ12 = 3691572375143311/104390384155500000, φ32 = 4037122085/176350902122 = −a3.
   Proof: Lemma 2 (REPORT.md), sufficiency direction: p > 0 on [1, ∞), S_odd(g) ≥ 0 and S_even(g) + (a3 + φ32)/2 ≥ 0
   on (0, π] imply (V). Arb check `cert4.py` (128 bits): p's real roots −5.441, 0.96554 < 1; n3 = 1.98·10⁻⁸ > 0;
   X_N = 1.2856683897…; g < π/√X_N termwise positive; [π/√X_N, π]: 204 adaptive cells, k ≤ 61, min S_odd ≥ 4.0·10⁻⁶,
   min S_even ≥ 0.01605. CERTIFIED (`logs/cert4_k3.34_final.log`).
2. **Mean values, j ≤ 4.** ∫_T^{2T}|F^{(j)}|² = c_j K T L^{(r+2)²+2j}(1 + o(1)) with one K > 0 and
   c_j = BM's 9-fold main-term integral (midpoint phase) × (v − S)^{2j}: BM's Theorem smoothed extended from j ≤ 2 to
   j ≤ 4 by BM's own argument (§ Step 1 below).
3. **Reduction.** If Λ ≤ κ, every gap of F's zeros in [T, 2T] is ≤ (1+o(1))κπ/L; (V) rescaled by λ = (1+ε)κ/L,
   applied to Re F and Im F on each gap, telescoped and summed, gives Σ_j a_j((1+ε)κ)^{2j}c_j ≥ −o(1), hence
   Σ_j a_jκ^{2j}c_j ≥ 0 (REPORT.md "Statement"; reviewer: CORRECT).
4. **Arb.** H := −Σ_j a_jκ^{2j}·Q[(v − S)^{2j}] = 9.3142330897·10⁻³⁵ ± 1.9·10⁻⁴⁶ > 0, H/Q[1] = 0.00136298995295 ±
   5.3·10⁻¹⁵ (Q[·] = GapBar's normalisation of c_j, a common positive factor; Q[1] = c0 > 0). So Σ a_jκ^{2j}c_j < 0,
   contradiction: **Λ > 167/50 = 3.34.**

## Step 1: differentiation in the shifts (sources quoted from the TeX)
- BBLR 1609.02539 Thm 1.1 (l.103-104): "α, β, γ, δ ∈ ℂ with α, β, γ, δ ≪ (log T)^{−1}", coefficients ≪ a^ε, error
  O_ε(T^{1/2+2ϑ+ε} + T^{3/4+ϑ+ε}).
- BM 1410.3635v3 l.366-368: I_{α,β}(M) with M(½+α3+it)M(½+β3−it), "where the shifts α_j, β_j ≪ L^{−1}";
  Lemma main (l.370-384) "uniformly for α_j, β_j ≪ L^{−1}".
- BM l.417: "Since I(α,β) and c(α,β) are holomorphic with respect to each variable α_j and β_j in a small disc centered
  at 0, the derivatives … can be expressed using Cauchy's integral formula as integrals of radii ≍ L^{−1} around the
  points … Since the asymptotic formula in Lemma main holds uniformly on these contours, each derivative adds a factor
  that is O(L) to the error term". BM l.418-433 print the operator identity
  Q(L⁻¹Σ∂_α)^j X1^{α1}X2^{α2}X3^{α3} = Q(L⁻¹Σ log X_i)^j X1^{α1}X2^{α2}X3^{α3} and the resulting integrand
  Q(−S)^{2j} = (v − S)^{2j} for general j; Theorem smoothed is printed for j ≤ 2.
- Our j ≤ 4: F^{(j)} = (iL)^j e^{ivtL}Q(L⁻¹Σ_{i≤3}∂_{α_i})^j[ζζM]; E = I − Main is holomorphic and ≪ T L^{(r+2)²−1} on
  the polydisc of radius ρ/L; Cauchy gives O_j(T L^{(r+2)²+2j−1}) after the operator; the midpoint shifts
  α3 = −β3 = iκπ/(2L) keep |h^{−α3}| = O(1) for h ≤ y. Unsmoothing by sandwiching (leading term only, which is all (3)
  needs).
- Reviewer (20 Nov): "no gap in extending the smoothed asymptotic to j = 3, 4"; complex uniformity "Yes"; mollifier
  shift "Admissible". Remark: the exact multiplier is (v − (log T/L)S)^{2j} = (v − S)^{2j} + O_j(1/L), absorbed.

## Step 2: Arb main terms (`certJ.py`, read-only import of GapBar's `gcert.py`)
- Same reduction as GapBar `certr.py` + `gtw.py` (midpoint twist): Q[(·)] = Σ_{M,N} Γ_MN⟨I^Mφ, X^{r²−1}(·)I^Nφ⟩,
  φ = I^{2r}f, f(U) = Σ_n A_n P_n(2U/ϑ − 1). New: the τ-weights (τ1 + τ2)^m, m ≤ 8, i.e. g_m = Σ_{k+l=m} C(m,k)
  conv2(u_k, w_l) over the four twisted terms (u_k from ∫_0^{1−y} t^k e^{zt}dt, any k: `h_series_K`), and the
  S-kernels H_n = Σ_j C(n,j) sh^j(g_{n−j}), ⟨S^m⟩ = Σ_i C(m,i) Q(H_{m−i}, X^i). Hermitian imaginary parts contain 0.
- **Lemma T (truncation M, N < D).** |H_n[M][N]| ≤ 4((r−1)!)⁴ C(M+2r−1, 2r−1) C(N+2r−1, 2r−1) β^{M+N}(2+M+N+4r)^n,
  β = ⌈3c/2 + K + 1⌉ (c = κπ, K = 8): |∂^q h_k(0)| ≤ (k + |z| + 1)^q, twisted phases |α| = c/2, so
  |u_ab| ≤ (3c/2 + k + 1)^{a+b}; Chu–Vandermonde for the fibre weights; Σ_k C(m,k) = 2^m; each sh multiplies by
  ≤ (M + N + 4r). With |⟨I^Mφ, X^q I^Nφ⟩| ≤ ϑ^q(ϑ^M/M!)(ϑ^N/N!)‖φ‖², the tail is summed in closed factorised form.
- Run: D = 72, 382 bits, β = 25; tail bound for n = 8: 1.9·10⁻¹⁹ (relative radii of ⟨S^m⟩ ≤ 2·10⁻¹⁵).
- Amplifier (`amp.py`): float top eigenvector of (C0, Σ_{j≥1}a_jC_{2j}) for the certified form, converted to GapBar's
  convention (A_0 = 1, 65 decimals of 16 digits, `logs/amp_d64_k3.34.json`); float H/c0 = 0.00136304, Arb
  0.00136298995295 ± 5.3·10⁻¹⁵.
- Checks: d = 8, κ = 3.25, D = 60, Wirtinger form: H/c0 = 0.0153629584984787 ± 1.3·10⁻¹⁷ = 1 − 1/h1 with GapBar's float
  h1 = 1.0156027 (optmid); the same run's standardised central moments (var 0.0932201, skew −0.00128, kurt 2.18571)
  equal `moments.py`'s float values; and `crosscheck.py` matched BM's printed integrand to 8 digits (REPORT.md).

## Commands (from Riemann/Deep/ExploreShape, OMP_NUM_THREADS=1), memory peaks
| step | command | time | peak |
|---|---|---|---|
| per-gap (V) | `bash ../../../.scratch/memcap.sh 290 600 uv run --no-project --with python-flint python3 cert4.py 200 61` | 1 s | 36 MiB RSS |
| amplifier | `… --with numpy --with scipy --with cvxpy python3 amp.py 64 3.34 72 logs/amp_d64_k3.34.json` | 18 s | 125 MB |
| main terms | `… --with python-flint python3 certJ.py logs/amp_d64_k3.34.json` | 154 s | 201 MiB RSS, 168 MB memcap |

Logs: `logs/cert4_k3.34_final.log`, `logs/amp_d64_k3.34.log`, `logs/certJ_d64_k3.34.log`, `logs/certJ_d64_k3.34.json`.

## Review (reviewer subagent, 20 Nov, two parallel referees, evidence pasted inline)
- Step 1 (shifts): no gap for j = 3, 4 (normalisation remark above; unsmoothing gives the leading term, enough here).
- Reduction (V) ⇒ criterion: CORRECT.
- Lemma 2: the representation (sine coefficients + f''(0), f''(g); admissibility ω⁴b_k − ω d_k ∈ ℓ²) and the termwise
  minimisation are valid, the minimum is attained; the "iff" needs p > 0 (counterexample with p(1) = 0). We use only
  the sufficiency direction, with p > 0 certified: no effect on the certificate.

## UPDATE (20 Nov, 00:55 UTC): Λ > 3.358 (= 1679/500), certified, order J = 5
- Same pipeline with a 5th-order per-gap inequality. **Lemma 3** (REPORT.md style, `exact5.py`): for f ∈ H⁵(0, g), f(0) = f(g) = 0,
  with free boundary values A = (f''(0), f''(g)), B = (f''''(0), f''''(g)), each sine mode contributes
  p b² + ω q_d b d + ω q_e b e + γ_dd d² + γ_de d e + γ_ee e² (p = Σ_{j≤5}a_jx^j, q_d = φ12 − (2a3+φ32)x − 2a4x² − 2a5x³,
  q_e = φ14 − φ34x + 2a5x², γ_dd = a3+φ32 + a4x + a5x², γ_de = φ34 − 2a5x, γ_ee = a5, x = ω²); minimising over b leaves
  R = [[N11, N12], [N12, N22]]/(4p), N11 = 4pγ_dd − x q_d², N12 = 2pγ_de − x q_d q_e, N22 = 4pγ_ee − x q_e² (degree ≤ 4).
  Sufficient: p > 0 on [1, ∞), S_odd(g) = Σ_{k odd}R(x_k) ⪰ 0, S_even(g) + C/2 ⪰ 0, C = [[a3+φ32, φ34/2], [φ34/2, a5]], g ∈ (0, π].
  Reviewer (20 Nov): **CORRECT** (representation complete: admissibility e_k + ω⁵b_k − ω²d_k ∈ ℓ²; grouped sums converge
  absolutely; per-mode lower bound legitimate; signs and C correct). Its Arb caveat (X* must be certified) is met: X* is the
  largest upper end of the Arb enclosures of all possibly-real roots of N11, N22, det (`fmpq_poly.complex_roots`), leading
  coefficients > 0, and every omitted mode has x_k ≥ (K+1)² > X*.
- Form: SDP with cutting planes and a 10⁻⁷ margin in the H⁵-normalised LMIs (`scanJ.py`, MARGIN=1e-7), rationalised to
  denominators ≤ 10¹², no repair needed. `cert5.py`: p(1) > 0, no real root ≥ 1; C ≻ 0 (det C = 1.9·10⁻⁹); X* = 6.4929;
  3748 adaptive cells on [π/√X*, π], k ≤ 61: min det S_odd ≥ 4.0·10⁻¹³, min det(S_even + C/2) ≥ 4.8·10⁻⁷: CERTIFIED
  (`logs/cert5_m80_3.358.log`, form in `logs/cert5_m80_3.358.json`). Cross-check (J = 5, d = 64 form): `validate.py` on a
  degree-31 polynomial basis gives +8.7·10⁻⁸; a0 − 10⁻³ and a0 − 10⁻⁴ fail in both checks at g = π (`logs/ex5_cross.log`).
- Mean values for j ≤ 5: Step 1 verbatim (the reviewer's verdict covers any fixed j).
- Arb main terms (`certJ.py`, K = 10, D = 72, β = 25): r = 2, ϑ = 2499/10000, degree 80, κ = 1679/500, v = 1.361505,
  amplifier `logs/amp_d80_k3.358_J5.json`: H = 3.534·10⁻³⁶ ± 4.3·10⁻⁴⁰ > 0, H/c0 = 0.000759 ± 2.4·10⁻⁷ (float 0.00075895);
  imaginary check passed (`logs/certJ_d80_k3.358_J5.json`). **Hence Λ > 1679/500 = 3.358.**
- Also certified: Λ > 671/200 = 3.355 with degree 64 (H/c0 = 0.0021774 ± 8.5·10⁻⁸, `logs/certJ_d64_k3.355_J5.json`,
  form `logs/cert5_m3.355.json`).

| step (J = 5, d = 80) | command (OMP_NUM_THREADS=1, under `bash ../../../.scratch/memcap.sh 290 …`) | time | peak |
|---|---|---|---|
| form | `MARGIN=1e-7 NIT=2 ROUNDS=10 JMIN=5 uv run … --with cvxpy python3 scanJ.py 2 0.2499 80 5 3.358` | 113 s | 173 MB |
| per-gap (V) | `SRC=logs/scanJ5m_r2_d80_k3.358.json uv run … --with python-flint --with numpy python3 cert5.py m80_3.358 0 0 0 61` | 10 s | 21 MB |
| amplifier | `uv run … python3 amp.py 80 3.358 72 logs/amp_d80_k3.358_J5.json logs/cert5_m80_3.358.json` | 40 s | 132 MB |
| main terms | `uv run … --with python-flint python3 certJ.py logs/amp_d80_k3.358_J5.json` | 306 s | 200 MiB RSS, 142 MB memcap |

## UPDATE (20 Nov, 01:25 UTC): Λ > 3.365 (= 673/200), certified, order J = 6
- **Lemma 4** (`exact6.py`): f ∈ H⁶(0, g), same free boundary data A, B as Lemma 3 (f⁽⁶⁾'s sine coefficients are
  −(ωe + ω⁶b − ω³d)), six telescoping terms φ12, φ14, φ32, φ34, φ52 [f⁽⁵⁾f''], φ54 [f⁽⁵⁾f''''] (f⁽⁵⁾(x) = δ_B/g + Σ(e + ω⁵b − ω²d)cos ωx):
  q_d = φ12 − (2a3+φ32)x + (φ52 − 2a4)x² − 2a5x³ − 2a6x⁴, q_e = φ14 − φ34x + (2a5+φ54)x² + 2a6x³,
  γ_dd = a3+φ32 + (a4−φ52)x + a5x² + a6x³, γ_de = φ34+φ52 − (2a5+φ54)x − 2a6x², γ_ee = a5+φ54 + a6x,
  C = [[a3+φ32, (φ34+φ52)/2], [(φ34+φ52)/2, a5+φ54]]; otherwise as Lemma 3. Setting a6 = φ52 = φ54 = 0 reproduces Lemma 3's
  numerators exactly (`ex6_test.py`). Reviewer: **CORRECT** given deg p = 6 (a6 > 0, here 1.5·10⁻⁴); increments, C,
  representation (admissibility ωe_k + ω⁶b_k − ω³d_k ∈ ℓ²) and convergence checked.
- Small gaps: as g → 0 the conditions become C ⪰ 0 and L∞ := lim x R(x) ⪰ 0, with (computer algebra)
  4a6·L∞ = [[4a2a6 + 4a6φ12 − φ52², 2a6φ14 + 2a6φ32 − φ52φ54], [·, 4a4a6 + 4a6φ34 − φ54²]], i.e. the Schur LMI
  [[a2+φ12, (φ14+φ32)/2, φ52/2], [·, a4+φ34, φ54/2], [·, ·, a6]] ⪰ 0. Without these LMIs the SDP's form failed at tiny gaps
  (lead(N11) = −7·10⁻⁸): tiny gaps do occur in {t_n} ∪ {t_n − κπ/L}, so they are now imposed in `sdpcut.py` (margin 10⁻⁸).
- `cert5.py` (JJ=6): p > 0 on [1, ∞); det C = 3.1·10⁻¹¹ > 0; X* = 6.8548; 3279 cells; min det S_odd ≥ 3.1·10⁻¹¹,
  min det(S_even + C/2) ≥ 1.65·10⁻⁸: CERTIFIED (`logs/cert6_m80_3.365.log`, form `logs/cert5_J6m80_3.365.json`).
- **Lemma T' (sharper truncation majorant, needed for n = 11, 12):** |h_k^{(n)}(0)| ≤ T_k(n+1)^k c^n, T_k = max(k!, max(1, k/c)^k),
  hence |H_n[M][N]| ≤ 4((r−1)!)⁴V_n C(M+2r−1,2r−1)C(N+2r−1,2r−1)(3c/2)^{M+N}(M+N+4r)^n, U_m = Σ_k C(m,k)T_kT_{m−k},
  V_n = Σ_j C(n,j)U_{n−j}. Reviewer: **CORRECT** for z ∈ {0, ±ic} (our case). The old base 3c/2 + K + 1 left s_12 with
  relative radius 1.6·10⁻⁵ at D = 72 (H not certified); with T' the tails are 2.6·10⁻⁵² (n = 0) to 3.1·10⁻¹⁶ (n = 12).
- Arb main terms (K = 12, D = 72): d = 80, κ = 673/200, v = 1.361827 (`logs/amp_d80_k3.365_J6.json`):
  **H = 9.332778·10⁻³⁶ ± 3.8·10⁻⁴³ > 0, H/c0 = 0.0013218642 ± 6.2·10⁻¹¹** (float 0.00132234), imaginary check passed
  (`logs/certJ_d80_k3.365_J6.json`). Mean values for j ≤ 6: Step 1 (any fixed j). **Hence Λ > 673/200 = 3.365.**

| step (J = 6, d = 80) | command (OMP_NUM_THREADS=1, under `bash ../../../.scratch/memcap.sh 290 …`) | time | peak |
|---|---|---|---|
| form | `MARGIN=1e-7 NIT=2 ROUNDS=12 JMIN=6 uv run … --with cvxpy python3 scanJ.py 2 0.2499 80 6 3.365` | 216 s | 180 MB |
| per-gap (V) | `JJ=6 SRC=logs/scanJ6m_r2_d80_k3.365b.json uv run … --with python-flint --with numpy python3 cert5.py J6m80_3.365 0 0 0 61` | 10 s | 20 MB |
| amplifier | `uv run … python3 amp.py 80 3.365 72 logs/amp_d80_k3.365_J6.json logs/cert5_J6m80_3.365.json` | 70 s | 141 MB |
| main terms | `uv run … --with python-flint python3 certJ5.py logs/amp_d80_k3.365_J6.json` | 442 s | 248 MiB RSS, 194 MB memcap |

## Certified values (all unconditional; same r = 2, ϑ = 2499/10000, midpoint mollifier)
| J | degree | κ | H/c0 (Arb) | per-gap form | log |
|---|---|---|---|---|---|
| 4 | 64 | 167/50 = 3.34 | 0.00136298995295 ± 5.3·10⁻¹⁵ | `logs/cert4_k3.34.json` | `logs/certJ_d64_k3.34.json` |
| 5 | 64 | 671/200 = 3.355 | 0.0021774 ± 8.5·10⁻⁸ | `logs/cert5_m3.355.json` | `logs/certJ_d64_k3.355_J5.json` |
| 5 | 80 | 1679/500 = 3.358 | 0.000759 ± 2.4·10⁻⁷ | `logs/cert5_m80_3.358.json` | `logs/certJ_d80_k3.358_J5.json` |
| 6 | 80 | 673/200 = 3.365 | 0.0013218642 ± 6.2·10⁻¹¹ | `logs/cert5_J6m80_3.365.json` | `logs/certJ_d80_k3.365_J6.json` |

## UPDATE (20 Nov, 01:36 UTC): Λ > 3.367 (= 3367/1000), certified, order J = 6, degree 80
- Form re-optimised at κ = 3.367 (`logs/scanJ6m_r2_d80_k3.367.json`, float h6 = 1.0000156); `cert5.py` (JJ=6): det C = 3.1·10⁻¹¹,
  X* = 6.8094, 3253 cells, min det S_odd ≥ 4.7·10⁻¹¹: CERTIFIED (`logs/cert6_m80_3.367.log`, form `logs/cert5_J6m80_3.367.json`).
- Amplifier `logs/amp_d80_k3.367_J6.json` (v = 1.361752, float H/c0 = 1.460·10⁻⁵); `certJ.py` (K = 12, D = 72, Lemma T'):
  **H = 1.14360·10⁻³⁷ ± 5.6·10⁻⁴³ > 0, H/c0 = 1.43562·10⁻⁵ ± 5.6·10⁻¹¹** (`logs/certJ_d80_k3.367_J6.json`; 337 s, 245 MiB RSS,
  189 MB memcap). **Hence Λ > 3367/1000 = 3.367.** The float crossing of order 6 at degree 80 is ≈ 3.3670: this is the ceiling
  of the order-6 class with this amplifier family; order 7 needs a third boundary pair (f⁽⁶⁾(0), f⁽⁶⁾(g)) and better-conditioned
  SDP validation (the float J = 7 run did not converge), not attempted.
| 6 | 80 | 3367/1000 = 3.367 | 1.43562·10⁻⁵ ± 5.6·10⁻¹¹ | `logs/cert5_J6m80_3.367.json` | `logs/certJ_d80_k3.367_J6.json` |

## Reproduce the record (from Riemann/Deep/ExploreShape; every command under `bash ../../../.scratch/memcap.sh 290 1750 …`, OMP_NUM_THREADS=1)
1. `MARGIN=1e-7 NIT=2 ROUNDS=12 JMIN=6 uv run --no-project --with numpy --with scipy --with cvxpy python3 scanJ.py 2 0.2499 80 6 3.367`
   (writes logs/scanJ_r2_d80.json; rename to logs/scanJ6m_r2_d80_k3.367.json)
2. `JJ=6 SRC=logs/scanJ6m_r2_d80_k3.367.json uv run --no-project --with python-flint --with numpy python3 cert5.py J6m80_3.367 0 0 0 61`
3. `uv run --no-project --with numpy --with scipy --with cvxpy python3 amp.py 80 3.367 72 logs/amp_d80_k3.367_J6.json logs/cert5_J6m80_3.367.json`
4. `uv run --no-project --with python-flint python3 certJ5.py logs/amp_d80_k3.367_J6.json` → `certified: true` (337 s, 245 MiB RSS, 190 MB memcap).
Steps 2 and 4 are the certificate (Arb); steps 1 and 3 only produce the rational inputs (any inputs that pass 2 and 4 suffice).

## Code versions and reproducibility (20 Nov, 01:54 UTC)
- `certJ.py` = the committed (HEAD) version: J = 4 and J = 5 certificates (old majorant, base ⌈3c/2 + K + 1⌉, valid: the
  same proof as Lemma T' with T_k = 1 and |h_k^{(n)}(0)| ≤ (k + |z| + 1)^n). Rerun on HEAD: `logs/certJ_d64_k3.34_headcheck.json`
  equals `logs/certJ_d64_k3.34.json` in H, H/c0, s, c, tails, flags (178 s, 209 MiB RSS, 168 MB memcap).
- `certJ5.py` = `certJ.py` + Lemma T' (sharper tail majorant; nothing else changed): the J = 6 certificates (3.365, 3.367).
  Rerun: `logs/certJ5_d80_k3.367_J6_check.json` equals `logs/certJ_d80_k3.367_J6.json` in all fields (and an earlier rerun,
  `logs/certJ_d80_k3.367_J6_rerun.json`, too). Per-gap rerun `logs/cert6_m80_3.367_rerun.log`: same numbers.
- `cert4.py` (J = 4) and `cert5.py` (J = 5 and, with JJ=6, J = 6; exact numerators from `exact5.py`, `exact6.py`) are the
  per-gap certificates. `scanJ.py`, `sdpcut.py`, `amp.py` only produce rational inputs (float).
- `logs/` is git-ignored; the J ≥ 5 record needs (force-add): `logs/scanJ6m_r2_d80_k3.367.json`, `logs/cert5_J6m80_3.367.json`,
  `logs/cert6_m80_3.367.log`, `logs/amp_d80_k3.367_J6.json`, `logs/certJ_d80_k3.367_J6.json`, `logs/certJ_d80_k3.367_J6.log`,
  `logs/certJ5_d80_k3.367_J6_check.json`.

## Reviews (reviewer subagents, evidence inline, 20 Nov)
| item | verdict |
|---|---|
| Step 1: BM Lemma main differentiated in the shifts, j ≤ 4 (any fixed j) | no gap (normalisation (log T/L) absorbed; leading term only) |
| Reduction (V_J) ⇒ criterion J | CORRECT |
| Lemma 2 (J = 4 per-gap) | sufficiency CORRECT; iff needs p > 0 (we certify p > 0) |
| Lemma 3 (J = 5 per-gap) | CORRECT (sufficiency); Arb caveat on X* met |
| Lemma 4 (J = 6 per-gap) | CORRECT given deg p = 6 (a6 > 0) |
| Lemma T' (truncation majorant) | CORRECT for z ∈ {0, ±ic} (our case); HS bound needs M ≥ 1, M = 0 by the operator norm |

**Best certified record: Λ > 3367/1000 = 3.367** (order 6, degree 80); the committed J = 4 certificate (Λ > 3.34) is unchanged.
