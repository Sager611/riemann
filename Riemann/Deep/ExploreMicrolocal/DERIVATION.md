# ExploreMicrolocal — DERIVATION (explore47)

## 1. Choice (measured 21:07–21:55Z)
- N1 (Gårding/FP/anti-Wick on Weil's form): Q_c is Toeplitz on [−a,a], a = ½ln c, with multiplier band-limited to exactly
  [−2a, 2a] (primes n < c ⇔ log n < 2a): the symbol oscillates at the window's own Heisenberg scale, h_eff ≍ 1 for every c;
  the Weyl symbol is discontinuous in x (cut-offs |x| ≤ a − ½log n). Gårding/FP errors are O(1), margins are ε_c ≈ e^{−4πc}c⁵.
  Plane-wave anti-Wick positivity ⇔ Krein extension ⇔ Q_c ≥ 0 (exact, circular); Landau excess 2c − 7/4 and Krein are
  already explore19's dictionary. Not chosen: no small parameter, nothing beyond explore19/explore39.
- N13 (two-term Weyl law for xp): 1-degree-of-freedom Bohr–Sommerfeld ⇒ N_H − N̄_H = O(1); N_ζ − N̄ = S(T) unbounded
  (Tsang Ω±((log T/log log T)^{1/3})): prolatefals' barrier in general form. Not chosen: known (Berry–Keating 1999).
- CHOSEN: tool (b) Bérard's Ehrenfest-time cutoff of the wave trace → node N9/N7 quantity S(T); consequences for N1 (Mori's ϖ)
  and N8 (RH gap constant).

## 2. Dictionary (compact hyperbolic surface X ↔ ζ)
| X (Bérard 1977, Hejhal 1976, Randol 1978) | ζ |
|---|---|
| eigenvalues r_j, density (area/2π) r | ordinates γ, density (1/2π)log(t/2π) |
| Selberg trace formula (exact, all times) | Guinand–Weil explicit formula (exact) |
| closed geodesics ℓ, #{ℓ ≤ L} ~ e^L/L (entropy 1) | prime powers log n, #{log n ≤ L} ~ e^L/L |
| amplitude ℓ₀/(2 sinh(ℓ/2)) ~ ℓ₀e^{−ℓ/2} | Λ(n)/√n = log p·e^{−(log n)/2} |
| self-adjointness: r_j real | RH: γ real |
| Ehrenfest time T_E: orbit sum e^{T/2} ≍ density | Δ_E: e^{Δ/2} ≍ log t ⇒ Δ_E = 2 log log t − 2 log log log t |
| Tauberian at T: error ≍ density/T | Beurling–Selberg R^± at type Δ: excess (2π/Δ)·density |
| Bérard: R(λ) = O(λ/log λ) | S(t) = O(log t/log log t) (Littlewood 1924) |
| JPT 2008 / Hejhal–Randol Ω-bound | Montgomery Ω±(√(log t/log log t)) (RH), Tsang (uncond.) |

## 3. Transferred claim (Bérard–Ehrenfest bound for ζ), exact form
Units: R^±(x) = ½{H^±(δ(t+x)) + H^±(δ(t−x))}, δ = Δ/2π, H^± Beurling's sgn majorant/minorant (CCM 1503.00955 Lemma 2);
R̂^± supported in |ξ| < δ, so the explicit formula sees n < e^Δ. With η_n = log n/Δ ∈ (0,1) (Vaaler's J):
  R̂^±(log n/2π) = (2π/Δ)[(1/π + (1−η)cot πη) sin(t log n) ± (1−η) cos(t log n)].
Guinand–Weil (CCM Lemma 4) + R^− ≤ χ_[−t,t] ≤ R^+ + RH (γ real) + N(t) = ϑ(t)/π + 1 + S(t) give, for all Δ, t:
  −τ_Δ(t) + O^−(t,Δ) − p^− ≤ S(t) ≤ τ_Δ(t) + O^+(t,Δ) + p^+,
  τ_Δ(t) = (1/4π)∫|R^± − χ|(Re ψ(¼+iu/2) − log π)du = log(t/2π)/(2Δ)·(1 + O(log t/(Δt))),
  O^±(t,Δ) = −(1/2π)Σ_{n<e^Δ} Λ(n)n^{−1/2} R̂^±(log n/2π)   (the "orbit sum"),
  p^± = ½(R^±(i/2) + R^±(−i/2)) − 1 = O(e^{Δ/2}/(δt)²) (computed exactly below).
Bérard step = triangle inequality on the orbit sum: |O^±| ≤ TB(Δ) := (1/Δ)Σ_{n<e^Δ}Λ(n)n^{−1/2}√((1/π+(1−η)cot πη)² + (1−η)²)
  ~ (const) e^{Δ/2}/Δ, so |S(t)| ≤ min_Δ[log(t/2π)/(2Δ) + TB(Δ)] → (¼ + o(1)) log t/log log t at Δ = Δ_E.
This is exactly CCM 1503.00955's "new and simple proof" of the 1/4 (GG 2007 had 1/2 from S(t+h) − S(t), two endpoints).
Bérard's own kernel (Hörmander's Tauberian theorem, smooth ρ ≥ 0, ρ̂ ⊂ [−1,1]) has excess > Beurling–Selberg's (Beurling's
L¹-extremality of H^±), so the microlocal version gives ≥ 1/4; no constant below 1/4 without |O^±| ≪ e^{θΔ/2}, θ < 1, for ALL t.
Consequences of the same limit: multiplicity/cluster ≤ (½+o(1)) log t/log log t (Fejér majorant, same Δ_E) ⇒ Mori's
RH-alone rate λ_min ≥ exp(−exp(μ^ϖ)) needs ϖ > 2·(¼) = ½ (MoriRate (6)); GG's gap ≤ (π+o(1))/log log t.
With an orbit hypothesis |O| ≪ e^{θΔ/2}: Δ = (2/θ)log log t, S-constant θ/4, ϖ > θ/2, gap θπ/log log t.

## 4. Tests (scripts and logs in this directory; every process < 290 MB, peaks 18–228 MB)
Smallest case (bt_validate.py, logs/bt_validate.txt): Odlyzko zeros1 (first 10⁵ zeros), t = 100, 1000, Δ = 3–7, both signs.
  Σ_γ R^±(γ) summed over zeros vs the explicit-formula side (pole + ∫R^±·ReΓ_R'/Γ_R/π − primes): |diff| ≤ 4·10⁻³ (zeros beyond
  74920 and quadrature); pole − 2 ≤ 10⁻⁴; τ_exact within 0.1–0.4 % of log(t/2π)/(2Δ) (and below it for R⁺);
  bracket e.g. t = 100, Δ = 5: −0.0178 ≤ S = −0.0024 ≤ 0.0112; t = 1000, Δ = 7: −0.371 ≤ 0.384 ≤ 0.584.
Scale (bt_window.py, logs/bt_*.txt|json, summary logs/bt_summary.txt): six windows, 2000–2500 sampled zeros each,
  S(γ_k^±) from Odlyzko's indices and mpmath ϑ (dps 50; offsets kept exact), orbit sums for Δ = 2…14 (93,371 prime powers,
  phases T0·log n mod 2π in mpmath), random-phase model: iid θ_p per prime, p^k ↦ kθ_p, same sample count.
  Bound check B⁻(γ_k^−) ≤ S(γ_k^−), S(γ_k^+) ≤ B⁺(γ_k^+): 14,000 zeros × 2 × 13 Δ = 364,000 checks, all hold, min margin 0.0037.
| window | t | log(t/2π) | Δ_tri | B_tri | ¼log t/loglog t | max|S| | Δ_E: TB / sup|O_ζ| / sup|O_rand| | Δ=14: TB / sup|O_ζ| / rand | max B⁺(14) |
|---|---|---|---|---|---|---|---|---|---|
| zeros1 #1001–3000 | 2.5·10³ | 5.99 | 3.00 | 1.88 | 0.95 | 1.31 | 0.88 / 0.69 / 0.72 | 24.7 / 1.21 / 1.34 | 1.42 |
| zeros1 #20001–22000 | 1.9·10⁴ | 8.01 | 3.25 | 2.20 | 1.08 | 1.46 | 0.88 / 0.69 / 0.72 | 24.7 / 1.34 / 1.34 | 1.56 |
| zeros1 #90001–10⁵ | 7.2·10⁴ | 9.34 | 3.50 | 2.39 | 1.16 | 1.63 | 0.88 / 0.70 / 0.72 | 24.7 / 1.32 / 1.34 | 1.65 |
| zeros3 (#10¹²+) | 2.7·10¹¹ | 24.48 | 5.00 | 4.13 | 2.01 | 1.94 | 1.69 / 0.98 / 0.97 | 24.7 / 1.38 / 1.34 | 2.12 |
| zeros4 (#10²¹+) | 1.4·10²⁰ | 44.58 | 6.00 | 5.94 | 3.02 | 2.00 | 2.23 / 1.09 / 1.08 | 24.7 / 1.27 / 1.34 | 2.86 |
| zeros5 (#10²²+) | 1.4·10²¹ | 46.83 | 6.00 | 6.13 | 3.13 | 2.03 | 2.23 / 1.08 / 1.08 | 24.7 / 1.32 / 1.34 | 2.99 |
Reading: (i) the transferred inequality holds; (ii) optimized over Δ it gives B_tri = 1.4–3.0 × max|S| and ≈ 2× its own leading term
(bt_range.py: B_tri/(log t/loglog t) = 0.48–0.52 for 10³ ≤ t ≤ 10¹², 0.42 at 10¹⁰⁰; ¼ is reached only as t → ∞);
(iii) at Δ_E the triangle inequality loses only ×1.3–2.1 (TB vs observed sup|O|), at Δ = 14 ×18–20: the Tauberian term τ(Δ_E) =
1.0–3.9 dominates, and only Δ ≫ Δ_E (where TB explodes like e^{Δ/2}/Δ) could lower it; (iv) with exact orbit sums the same
inequality at Δ = 14 is within 0.1–1.0 of max S (Selberg's formula); (v) ζ's sup orbit sums are within −11…+5 % and rms within −3…+14 %
of the random-phase model's at every Δ (e.g. 1.08 vs 1.08 at 10²¹, Δ = 6; +14 % rms only at 2.5·10³, Δ = 14, ζ sampled at its zeros):
the observed cancellation is Bohr–Jessen-generic, not arithmetic.
Window S-values: zeros1 windows from exact indices; zeros3–5 indices from Odlyzko's headers (10¹²+1, 10²¹+1, 10²²+1).

## 5. Range statement (asked by the brief)
RH holds for |γ| ≤ H0 = 3·10¹² (Platt–Trudgian). For t ≤ 10¹², R^± decays like (δ(|x|−t))^{−2} and grows ≤ e^{Δ/2} off the axis,
so any zero above H0 changes Σ_ρR^±(γ_ρ) by ≪ e^{Δ/2} log H0/(δ²H0) ≈ 10⁻¹⁰ (Δ ≤ 8): unconditionally |S(t)| ≤ B_tri(t) + o(1)
on 10³ ≤ t ≤ 10¹²; B_tri(10¹²) = 4.26 vs Trudgian 2014's 6.53 and HSW 2022's 13.09 there (bt_range.py). Not certified (needs
explicit constants for τ, p^± and the off-line term); a range record at most (RECORDS.md), and Platt–Trudgian's own zero data
(not public) give S directly in that range.

## 6. N1 measurement behind "not chosen" (n1_symbol.py, logs/n1_symbol.txt)
Q_c(f) = (1/2π)∫|f̂|²σ_c, σ_c(ξ) = Reψ(¼+iξ/2) − log π + FT[2cosh(s/2)1_{|s|≤2a}](ξ) − Σ_{n<c}2Λ(n)n^{−1/2}cos(ξ log n), so the
plain symbol (anti-Wick-type) bound is λ_min ≥ inf σ_c. Scan |ξ| ≤ 400 (the negative set continues to ξ ~ 2π e^{Σ2Λ/√n}):
| c | 1.2 | 1.5 | 2 | 3 | 5 | 10 | 20 | 47 |
|---|---|---|---|---|---|---|---|---|
| inf σ_c | −4.64 | −3.74 | −2.54 | −1.73 | −1.62 | −2.57 | −4.26 | −5.05 |
| negative cells 2a·meas{σ_c<0}/2π | 0.18 | 0.18 | 1.13 | 2.16 | 6.79 | 27.5 | 66.6 | 116.9 |
Q_c ≥ 0 at all these c (Yoshida c ≤ 2; certificates to 47), so the symbol is negative on up to 117 quantum cells with no negative
eigenvalue: positivity is a coherent (Krein-extension) effect, not a symbol inequality; at c ≤ 1.5 the negative set is a fifth of a
cell (an uncertainty-principle regime, which is where Fefferman–Phong-type bounds could at best work, below Yoshida's c = 2).
