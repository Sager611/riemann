# CERT: Λ > 3.34 (= 167/50), unconditionally (explore53, 20 Nov 2026)

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
