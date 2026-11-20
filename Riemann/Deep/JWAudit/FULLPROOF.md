# JWAudit — the chain, re-derived (line numbers: HD = HollandK/DERIVATION.md, HC = HollandK/certify.py, SD = SavvaK/DERIVATION.md, SC = SavvaK/certify6.py, G = GORTTW 1910.01227 TeX `Jensen_Revision.tex`)

## 1. Shared analytic input: explicit complex saddle on D(x0) = {|z − x0| ≤ x0/10} (HD l. 38–98, HC l. 49–90)

Setup (HD l. 21–27). φ(u) = Φ(u/2) (Φ the de Bruijn–Newman kernel), Mt(z) = ∫_0^∞ φ u^{2z} du, ξ(1/2+w) = 4∫_0^∞ φ(v) cosh(wv) dv,
so γ(m) = 4Γ(m+1)Mt(m)/Γ(2m+1) (my normalization; HD's C·4^{cz} differs by a geometric factor C r^{n+j}, which changes neither
hyperbolicity of J^{d,n} nor any Δ^{≥2} log γ, nor h^{(k≥2)}). ✓ φ > 0 on u ≥ 0 (2πk²e^{2u} > 3), so Mt(x) > 0. ✓
Γ-part: log Γ(w+1) − log Γ(2w+1) = −2w log 2 − log Γ(w+1/2) + const (duplication), hence g = const + c w + h − log Γ(w+1/2) (HD l. 152). ✓

A1 (HD l. 47–51; HC l. 50–53). k = 1 term: φ_1 u^{2z} = 2π² e^{f} a, f = 2z log u + 9u/2 − πe^{2u}, a = 1 − (3/2π)e^{−2u}. ✓
f'(u) = 0 ⟺ u = T_z(u) := ½ log((z/u + 9/4)/π) ⟺ Z(u) = z. ✓ With X := x0/a0 = πe^{2a0} − 9/4, ζ := (z/x0)(a0/u):
T_z(u) − a0 = ½[log(z/x0) − log(u/a0) + log(1 + 9/(4Xζ)) − log(1 + 9/(4X))] (all four arguments small, principal logs add). ✓
|log(1+w)| ≤ −log(1−|w|) gives −log 0.9 and −log(1 − 0.06/a0); |ζ| ≥ 0.9/(1 + 0.06/a0) gives |9/(4Xζ)| ≤ w1; |log(1+w)| ≤ 2|w|. So
|T_z(u) − a0| ≤ β = ½(−log 0.9 − log(1 − 0.06/a0) + 4w1) < 0.06 (invariance). T_z' = −1/(2u(1 + 9u/(4z))), |9u/(4z)| ≤ w1:
|T_z'| ≤ 1/(2(a0 − 0.06)(1 − w1)). β, w1 decrease in a0, so a0 = a* is the worst case (HC uses a_s). ✓ Banach: unique fixed point,
holomorphic (uniform limit of iterates), real for real z (the real root lies in the same disk). ✓ Hence |Im u0| ≤ β, Re u0 ≥ a* − β =: alo.

A2 (HD l. 52–73; HC l. 54–80). Contour [0, al/2] ∪ [al/2, al/2 + iβ'] ∪ {u0 + s: s ≥ −al/2}, u0 = al + iβ': integrand analytic in Re u > 0,
|e^{−πe^{2u}}| = e^{−πe^{2Re u}cos(2 Im u)} decays doubly exponentially for |Im u| ≤ β < π/4. ✓
Ray: F(s) := Re(f(u0+s) − f(u0)), F(0) = 0, F'(0) = Re f'(u0) = 0, F'' = Re(−2z/u²) − 4πe^{2(al+s)}cos 2β' ≤ −4π cos(2β) e^{2(al+s)},
since |arg z| ≤ asin 0.1, |arg u| ≤ atan(2β/al) on the ray (Re u ≥ al/2), so Re(z/u²) > 0 iff asin 0.1 + 2 atan(2β/al) < π/2 (HC l. 56). ✓
Core |s| ≤ s1 = 0.3: f(u0+s) − f(u0) = ½f''(u0)s² + R3, |R3| ≤ C3|s|³, C3 = sup|f'''|/6; f''' = 4z/u³ − 8πe^{2u},
|z| = |u0||πe^{2u0} − 9/4| ≤ (al+β)(πe^{2al} + 9/4), |u| ≥ al − s1, e^{2Re u} ≤ e^{2(al+s1)} ⇒ C3 ≤ c3 e^{2al}, c3 as HC l. 59. ✓ (c3 decreasing in al ✓.)
If C3 s1 ≤ κ0/4 (κ0 = Re(−f''(u0)) ≥ 4π cos 2β e^{2al}; HC l. 60 checks c3 s1 ≤ π cos 2β ✓): |e^{R3} − 1| ≤ C3|s|³e^{κ0 s²/4} and
∫|s|³e^{−κ0 s²/4} ds = 16/κ0² ⇒ core error ≤ 16C3/κ0². Relative to |√(2π/(−f''))| with |f''(u0)| ≤ 4πe^{2al}φ2,
φ2 = 1 + 1/(2alo) + 9/(8π alo e^{2alo}) (from −f''(u0) = 4πe^{2u0}(1 + 1/(2u0) − 9/(8πu0e^{2u0})), HC l. 57 ✓):
rel_core ≤ c3 √(2φ2) e^{−al}/(π² cos² 2β) — HC l. 61 ✓. Gaussian tails (Mills) 2e^{−κ0 s1²/2}/(κ0 s1) — HC l. 63 ✓. Ray tails (tangent line
of concave F at ±s1, κ_min = 4π cos 2β e^{2(al − s1)}) 2e^{−κ_min s1²/2}/(κ_min s1) — HC l. 65 ✓. Amplitude: |a − 1| ≤ (3/2π)e^{−al} on the
ray (Re u ≥ al/2), ∫_ray e^F ≤ √(2π/κ_min) + tails ⇒ rel_am = (3/2π)e^{−al}(√(φ2 e^{2s1}/cos 2β) + rel_rt) — HC l. 66 ✓.
Pieces near 0 (my derivation, HD l. 63 states only the result): for u on [0, al/2] ∪ vertical, |u/u0|² ≤ 1/4 + β²/al² < 1, Re z ≥ 0.9x0,
|Im z| ≤ 0.1x0, |arg u − arg u0| ≤ 3β/al, 9 Re(u − u0)/2 ≤ 0, −πe^{2Re u}cos(2 Im u) ≤ 0 and πe^{2al}cos 2β' ≤ |z/u0 + 9/4| ≤ 1.1x0/al + 9/4
(saddle equation) ⇒ |e^{f(u)−f(u0)}| ≤ e^{b x0 + 9/4}, b exactly HC l. 68; with |a| ≤ 1.48, length ≤ al ≤ log x0, √(|f''|/2π) ≤ √(2φ2)e^{al} ≤ √x0:
relative ≤ 20 √x0 log x0 e^{b x0} ≤ e^{5.2} x0 e^{b x0} ✓ (decreasing in x0, so HC l. 70 at x0 = N* is the sup). ✓
θ-terms k ≥ 2 (real axis): |φ_k| ≤ (2π²k⁴ + 3πk²)e^{−π(k²−4)} e^{9u/2 − 4πe^{2u}} (u ≥ 0) ⇒ Σ ≤ C_θ 2^{−9/2} e^{9v/2 − πe^{2v}}|_{v = u + log 2};
(v − log 2)^{2x} ≤ v^{2x}e^{−2x log 2/v}; gain ≥ x log 2/u0(x) ≥ 0.9 x0 log 2/(a0(1 + β/a0)) on v ≤ 2u0(x); deficit
Ψ(x) − Re Ψ(z) ≤ (y²/2) sup|Ψ''| ≤ q2max x0/(360(a0 − β)) (the vertical segment x → z lies in D(x0)) ⇒ exponent −g x0/a0, g as HC l. 76. ✓
The prefactor (e^{2β+1}√φ2 + 1) (HD l. 64, l. 172–174) I did not re-derive in full; any prefactor below e^{10^9} is irrelevant
(exponent −g x0/a0 ≤ −0.617·8.3·10^10 already at L* = 24).
Numbers (rerun, HC): δA ≤ 5.3215·10⁻²⁰ (L* = 90), 7.9042·10⁻¹⁸ (80), 1.1637·10⁻⁵ (24, SavvaK -self); dominated by rel_core + rel_am ≈ 1.11e^{−al} + 0.66e^{−al}.
N* = Z(45) = 1.72530·10^41 (L* = 90) ✓ (recomputed independently in saddle_num.py: a0 = 45.0 at x0 = N*).
By hand at alo = 44.94665: e^{−alo} = 3.0193·10⁻²⁰, rel_core = 7.63354·√2.02226·e^{−alo}/(π²·0.994313²) = 3.3590·10⁻²⁰, rel_am =
(3/2π)e^{−alo}√(1.01113 e^{0.6}/0.994313) = 1.9624·10⁻²⁰; sum 5.3214·10⁻²⁰ = certified δA ✓.

A3 (HD l. 74–78; HC l. 82). |δ| < 1 ⇒ Mt ≠ 0 on the convex D(x0): holomorphic log, real on the axis. −f''(u0) = (4z/u0)(1 + 1/(2u0) + 9u0/(4z)) ✓.
|G − G(x0)| ≤ ½(|log(z/x0)| + |log(u0/a0)| + 2|ε(z) − ε(x0)|) + 2·2δA, ε = 1/(2u0) + 9u0/(4z), |ε(z) − ε(x0)| ≤ β/(2 al a0) + 2w1 ⇒ etaG as HC l. 82. ✓
Values: etaG ≤ 0.0532866 (90), 0.0533654 (80), 0.0552241 (24). Cauchy: |G^{(k)}(w)| ≤ k! etaG/(x0/10 − |w − x0|)^k. ✓
Uniformity: every bound is a sup over D(x0) and decreasing in a0, evaluated at a0 = a* (HC uses a_s, alo). ✓

Part B (HD l. 80–86; SD l. 35–42). Ψ' = 2 log u0 (envelope), Ψ'' = 2/(u0 Z'(u0)), d/dz = Z'^{-1} d/du, Z' = (9/4)(2 + t − ts)/(ts) ✓;
Q2 = 4(1−s)/(2 + t(1−s)) ✓ (by hand). Q3..Q6 re-derived independently (qk_check.py, Sec. 2 below).

## 2. Falsification tests of the saddle input (this audit's own code)

saddle_num.py (mpmath, 70–200 digits; logs/saddle_num_12pts.log, 433 s, peak 69 MB). Mt by direct quadrature, the k = 1 term,
on the horizontal ray through the complex saddle; validated against the real-axis integral with k ≤ 3 theta terms at x0 = 2000
and 10^6 (real z: agreement to all 50 digits) and at complex z = 2000(1 + 0.1e^{iθ}), θ = π/2, 3π/4 (agreement 2·10⁻⁴⁸).
Edge |z − x0| = x0/10, 12 points, plus the centre:

| x0 | regime | true max \|δ\| | claimed δA | true max \|G(z) − G(x0)\| (at z = 0.9x0) | claimed ηG | max \|u0(z) − a0\| | β |
|---|---|---|---|---|---|---|---|
| 9.986·10^11 = N*(L* = 24) | SavvaK -self | 1.90·10⁻¹¹ | 1.164·10⁻⁵ | 0.050484 | 0.055224 | 0.050569 | 0.055187 |
| 6.963·10^36 = N*(80) | SavvaK main | 9.06·10⁻³⁶ | 7.904·10⁻¹⁸ | 0.052021 | 0.053365 | 0.052030 | 0.053431 |
| 1.725·10^41 = N*(90) | HollandK | 4.11·10⁻⁴⁰ | 5.321·10⁻²⁰ | 0.052095 | 0.053287 | 0.052101 | 0.053348 |
| 10^45, 10^50, 10^60 | HollandK | ≤ 7.8·10⁻⁴⁴ | 5.321·10⁻²⁰ | 0.052146, 0.052201, 0.052283 | 0.053287 | ≤ 0.052287 | 0.053348 |

True δ(x0) equals the first-order Laplace term (5f'''²/24κ³ + f''''/8κ² + amplitude) to 4 digits: δ ≍ e^{−2a0}, the proof's δA ≍ e^{−a0}
(odd terms bounded, not cancelled) — valid, about 20 orders loose at L* = 90. ηG is tight: margin 0.0010–0.0013 (2%) at L* ≥ 80,
0.0047 at L* = 24; the maximum sits at z = 0.9x0 where ½|log(z/x0)| = 0.05268 dominates. Ray tails beyond 40 standard deviations:
exponent ≤ −796. No counterexample.

qk_check.py (logs/qk_check.log): Q2..Q6 re-derived directly in u (Ψ' = 2 log u0, d/dz = Z'^{-1}d/du) agree with the (t,s) construction of
SC l. 106–116 to 10⁻⁵⁰ at five u; Q_k(t,0) = 2 − t + t²/2, −2 + t², 4 + t − 5t²/2, −12 − 5t + 15t²/2, 48 + 26t − 28t² + O(t³)
(SD l. 41 ✓); Ψ^{(k)}(10^8) by numerical differentiation of Ψ = f(u0(x), x) matches Q_k/(x^{k−1}2u0) to ≤ 6·10⁻⁹⁸, k = 1..6.

## 3. Spot check of one certified SavvaK cell by a different method (spotcheck.py, logs/spotcheck_1e40.log, 158 s, 41 MB)

n0 = 10^40 lies in cell 917 (ε ∈ [0.0114625, 0.011475], L = 87.18380). Method: h(n0 + j), j = 0..6, by direct quadrature at 330
digits (error estimate 8·10⁻³⁴⁸), the exact matching equations Δ²h(n0+k) = Σ ±log(1 + 1/U) (k = 0..3; no Q_k, no R_k series, no 6kL/n
ball) solved by Newton (residual 2·10⁻²¹⁰), and n0^5 L E_F^{(6)} from Δ⁶h and mpmath ψ^{(5)}. Certified box for the cell from a
single-cell rerun of a patched copy (logs/cert_one_cell917.log; its witness row equals cells_main.csv row 917):

| quantity | exact at n0 = 10^40 | certificate, cell 917 |
|---|---|---|
| a, t | 0.3246122, 1.9598015 | [0.324233, 0.324993], [1.959145, 1.960466] |
| w, δ | 5.1962473, 0.3410670 | [5.192006, 5.200518], [0.340898, 0.341235] |
| B/n | 2.019403 | ≤ 2.02015 |
| X0..X3 | 1.9773201, −1.9994885, 4.0216397, −12.1107553 | 1.9773 ± 3.5·10⁻⁵, −1.999489 ± 9·10⁻⁷, 4.0216 ± 4.9·10⁻⁵, −12.111 ± 3.2·10⁻⁴ |
| n^5 L E_F^{(6)} | −1.93878 (at z ≈ 3) | ball −1.9 ± 0.098; C6 ≤ 1.99797 |
| n^5 L h^{(6)} | 48.58157 | Q6 leading 48 + 26t − 28t² = 48.58171 |

Every exact value lies in its certified enclosure.

Second spot check, the weakest regime (SavvaK -self, L* = 24, δA = 1.2·10⁻⁵): n0 = 10^12, ε = 0.0416644 (top cell 399, last of its
8 sub-cells, ε ∈ [0.0416536, 0.0416667]); logs/spotcheck_1e12.log (160 digits, 18 s) and logs/cert_one_self_cell399.log:

| quantity | exact at n0 = 10^12 | certificate, sub-cell |
|---|---|---|
| a, t | 0.3033466, 1.8697914 | [0.302672, 0.304027], [1.868638, 1.870967] |
| w, δ | 4.8833700, 0.3623262 | [4.875773, 4.891046], [0.361976, 0.362670] |
| X0..X3 | 1.9200043, −1.9937287, 4.0667255, −12.3648729 | 1.9200 ± 2.4·10⁻⁵, −1.99373 ± 2·10⁻⁶, 4.06672 ± 8.7·10⁻⁶, −12.365 ± 2.3·10⁻⁴ |
| n^5 L E_F^{(6)} | −1.77762 | ball −2 ± 0.33; C6 ≤ 1.95717 |

Both spot checks: every exact value lies in its certified enclosure (cert_one_patch.py makes the single-cell copy; it never edits SavvaK/).


## 4. (A) SavvaK assembly (SD l. 44–139, SC l. 137–386; Savva TeX = Zenodo 22293642 v1.1 / GitHub aimathpapers/Jensen-2-Thirds v1.1.0)

Six matches (SD l. 44–62). P(y) := J^{d,n}(−y/S)/γ(n) = Σ p_j c(j) y^j, p_F = 3F2(−d, A, C; B, D; Dy/(AC)) = Σ(−1)^j C(d,j) m_j y^j ✓
((−d)_j/j! = (−1)^j C(d,j)); m_1 = 1/B, so S = Bγ(n+1)/γ(n) gives c(0) = c(1) = 1 ✓; c(2..5) = 1 ⟺ Δ² log γ(n+k) = Δ² log m_k (k = 0..3)
⟺ Δ²h(n+k) = f(A+k) + f(C+k) − f(B+k) + f(n+½+k) − f(D+k), f = log(1 + 1/U) (Δ² log Γ(U+k) = f(U+k); −log Γ(w+½) from the
duplication) ✓. Δ^j f = (−1)^j j! Φ_{j+1} (ΔΦ_q = −qΦ_{q+1}) ✓; n^qΦ_q(U) = (n/U)^q(1 − θq²/U), θ ∈ [0,1] ✓; the scaled equations
a^k ε^{k−1}/k + w R_k(t, wε) + δ R_k(1 + ν/2, δε) = Y_{k−1} ✓. Error: |ρ_k| ≤ (k L/n) Σ_U (n/U)^{k+1} ≤ 5kL/n when all five n/U ≤ 1,
which the box check SC l. 273 (δ > 0, t > 1, 0 < aε < 1) guarantees; L ≤ log n (u0 ≥ 1) and log x/x decreasing give ≤ 6k log N*/N* ✓.
Data: Δ^{j+2}h(n) is a B-spline average of h^{(j+2)} on [n, n+j+2] ✓; t_w, s_w, (n/w)^{j+1}u0(n)/u0(w) enclosures SC l. 255–258 ✓; G part
(j+2)! ηG (10/(1 − 50ν))^{j+2} L/n by Cauchy on D(n) ✓.
Brouwer/Krawczyk (SD l. 63–75, SC l. 172–227, 262–275). T(y) = y − P G_exact(y); componentwise MVT for G_approx over the box with every
datum (cell ε, ν ∈ [0, ν*], Y-intervals) as an interval and ρ as a ball in G(y_m) (SC l. 180) gives T(Y) ⊂ K(Y); K(Y) ⊂ int Y ⇒ T(Y) ⊂ Y,
T continuous (Φ_k smooth for U > 0) ⇒ fixed point (Brouwer), P nonsingular (Rump) ⇒ G_exact = 0 in K(Y) ∩ Y, for every n in the cell ✓.
Uniqueness not needed: any real (A,B,C,D) with c(0..5) = 1 serves ✓. Per-cell witnesses: cells_main.csv (1000 rows, 0 splits) and
cells_self.csv (2210 rows, 1810 splits) reproduced byte-identically; cell 917 recomputed alone and matched by the exact solution (§3).
Positivity, simplicity (SD l. 76–82; Savva TeX l. 796–816). Certified A − B − d > 0, C − D − d > 0 (SC l. 280–281), B, D > 0: both factors
of p_F = Q^{(1)}_{A,B} ⊠_d Q^{(D)}_{C,D} (Savva l. 784; sign convention re-derived: the MSS ⊠ multiplies the positive normalized
elementary symmetric functions, giving Σ(−1)^j C(d,j) m_j y^j ✓) are rescaled Jacobi P_d^{(α,β)}, α = V − 1 > −1, β = U − V − d > 0:
d simple roots in (0, U). Cited: Martínez-Finkelshtein–Morales–Perales [MMP2024] Prop. 2.7(iii) (nonnegative real roots preserved
by ⊠_d) and Prop. 2.17 (lmesh(p ⊠ q) ≥ lmesh(p) > 1 for simple positive roots, d ≥ 2) — hypotheses hold. NB: "MMP" is
Martínez-Finkelshtein–Morales–Perales (Savva TeX l. 162–165), not "Mirabelli–Marcus–Peters" as the discovery page says.
Localization (SD l. 83–95, SC l. 228–237, 283–289). Monic Jacobi recurrence (standard b_k, a_k²); Y = U(1 − x)/2 transport; diagonal
identity re-derived symbolically (qk_check.py ✓); with S + d + 1 = U exactly, k ≤ d, (2k+S)(2k+S+2) ≥ S², (2k+S)² − 1 ≥ S² − 1,
k + V − 1 ≤ V + d, k + U − V − d ≤ U − V, k + S ≤ S + d I re-derive both bounds of SD l. 89–90 / SC l. 233–235 exactly ✓; Gershgorin C_J = sup|diag − V|/√(Vd) + 2 sup off/√(Vd) ✓. Cited: MSS [MSS2022, Thm 1.6] maxroot(p ⊠_d q) ≤
maxroot(p)maxroot(q), applied to the reversals too (⊠_d commutes with reversal; needs positive lower endpoints, certified
SC l. 286–287) ⇒ every root, and by Gauss–Lucas every critical point, has |y − B| ≤ C_loc√(Bd), C_loc = C_J1(1+η) + C_J2√(B/D) ✓.
Critical radius (SD l. 96–108, SC l. 289–301). Recurrence P3T_{m+3} + P2T_{m+2} + P1T_{m+1} + P0T_m = 0 verified exactly with the
factor (m − d + 1) (qk_check.py; Savva App l. 688's (m − d) fails, as AUDIT.md says) ✓. Decomposition of P1 ✓; every term of P2/n,
|P1|/(n√(Bd)), P0/(n²d) re-derived from m ≤ d ≤ ω²n ✓. Max argument: T_0 = 1, T_1 = 0 (critical point), T_{d+1} = 0; maximizer k ∈ [2, d],
recurrence at m = k − 2 ⇒ 1 ≤ budget; certified budget ≤ 0.97277 < 1 ⇒ |T_k| ≤ R^k ✓.
Residual (SD l. 109–120, SC l. 303–325). E_F^{(6)} = h^{(6)}(n+z) − ψ^{(5)}(n+½+z) − ψ^{(5)}(A+z) − ψ^{(5)}(C+z) + ψ^{(5)}(B+z) + ψ^{(5)}(D+z) ✓;
pairings x^{−5} − y^{−5} = (y − x)Σ x^{−(i+1)}y^{−(5−i)} give term_B, term_D, term_A exactly ✓; |ψ^{(5)}(w) − 24w^{−5}| ≤ 120(Re w)^{−6}
(ψ^{(5)}(w) − 24w^{−5} = ∫t⁵e^{−wt}(1/(1−e^{−t}) − 1/t)dt, bracket in [½, 1); 240 used) ✓; the five (Re w)^{−6} terms in SC l. 320 ✓; G^{(6)} by Cauchy on radius n/20 inside D(n) ✓;
u0(n)/u0(n+z), t, s balls ✓.
HG, stability, assembly (SD l. 121–139, SC l. 345–366). E_F(0..5) = 0 (real, exp = 1) ✓; |Π(z − j)| ≤ (d + λR)^6 on Ω ✓ ⇒ sup_Ω|E_F|
≤ Cfin d³/(n²L) ✓. P/p = Σ_k (Δ^k c(0)/k!) T_k (Newton series) ✓; Δ^k c(0) = 0, k = 1..5 ✓; |c[0..k]| ≤ e/(λR)^k (Hermite–Genocchi on
[0,k] + Cauchy on radius λR, k ≤ d; T_k = 0 for k > d) ✓ ⇒ |P/p − 1| ≤ e/(λ⁵(λ−1)); λ⁵(λ−1) = 1.373317 (λ = 1.33) ✓;
|e^{E} − 1| ≤ e^{|E|} − 1 ⇒ θE = log(2.373317) − 0.001 = 0.863290 ✓. Signs at 0, the d − 1 critical points and +∞ (c(d) > 0) ⇒ d simple
positive zeros of P ⇒ d distinct negative zeros of J^{d,n} ✓. log(n+2) ≤ (1 + lamL)L (n < (πL/2)e^L) ✓. K_res = 1.06043·18.752/0.86329
≤ 23.035 from rounded displays, 23.0332 certified ✓ < 24. K_N(9e24) = N*²log(N*+2)/(9e24)³ = 5.641 ✓, K_N(d0) = 0.978 ✓ (by hand).
Regime: d > D_all and n²log(n+2) ≥ 24d³ ⇒ n > N* (else n²log(n+2) ≤ 5.65 D_all³ < 24d³) and (d/n)³ ≤ log(N*+2)/N* ✓.
-self: K = N*²log(N*+2) = 2.756·10^25 forces n ≥ N* and d/n ≤ 1/N*; K_res = 830 ≤ K ✓; d ≤ 5: P = p_F ✓.


## 5. (B) HollandK: where our constants enter Holland's chain (HD l. 114–175, HC l. 49–90, 121–270; Holland TeX = arXiv 2608.08682 jensen.tex)

- C3 s1 < κ0/4: c3 s1 = 7.63354·0.3 = 2.2901 ≤ π cos 2β = 3.1237 (HC l. 60) ✓, margin 27%; re-derived in §1, reviewer CORRECT (§8).
- δA ≤ 5.33·10⁻²⁰, ηG ≤ 0.053287: re-derived term by term (§1); true values 4.1·10⁻⁴⁰ and 0.05210 at x0 = N* (§2) ✓.
- a + b + θ = 0.97084 < 1 suffices: Holland TeX l. 1797–1807 shows |T_{m+2}|/r^{m+2} ≤ a_m|T_{m+1}|/r^{m+1} + b_m|T_m|/r^m + ϑM and at the
  maximizer M ≤ (a + b + ϑ)M; Holland's "each ≤ 1/4" is only a sufficient choice ✓. a_m, b_m (TeX l. 1748–1755): |B + m − y + (d−1−2m)y/A|
  ≤ C_loc√(Bd) + d(1 + y/A) and y/B ≤ 1 + C_loc√(d/B) give HD l. 126–127 ✓. ϑ: Holland's absolute C* (TeX l. 1719–1726) replaced by an
  explicit C_V; re-derived |β_m| ≤ A(d+1) + 5d² + 4d + 1, |γ_m| ≤ (d²/4)(A + d), A − d + 3 + 3m ≤ A + 2d + 3 ⇒ C_V ≤ 1 + r/A + (2d+3)/A
  + 2√(d/B)/K6 + 10(d/A)√(d/B)/K6 + d(1 + d/A)/(4K6²B); HC l. 228 has exactly these terms with 3d/A for (2d+3)/A, valid for
  d ≥ 3 (d ≤ 4 needs no stability: P_F = p_F) ✓; U-backward (TeX l. 1690)
  from U_{d+1} = 0 by induction: |U_m| ≤ ε_p C_V M r^{m+2}/(D(1 − r/D)) (check: (ε_p C_V M + Xr)/D = X for X := ε_p C_V M/(D − r)) ✓;
  ϑ = (y/D) ε_p C_V/((1 − y/A)(1 − r/D)) ✓. Certified a ≤ 0.81815, b ≤ 0.15075, ϑ ≤ 0.0019431 (rerun) ✓.
- Stability (HD l. 30–35, Holland Prop. 2.2): five matches, |P/p − 1| ≤ e/(λ⁴(λ−1)), λ⁴(λ−1) = 1.53664 (λ = 1.4), θE = log 2.53664 − 0.001
  = 0.929840 ✓. K_res = (Cfin/θE)²(1 + lamL)² = (70.600/0.929840)²·1.05502² = 6416.7 from rounded displays, 6416.59 certified ✓;
  K = max(6416.59, K_N(9e24) = 784.13) ≤ 7000 ✓. Self-contained: N*(80)³log²(N*+2) = 3.375·10^110·7196.7 = 2.429·10^114 ✓.
- Outward rounding (the reviewer's 6 Nov demand): fup/flo (HC l. 21–28) take the exact rational endpoint of each Arb ball and round
  with decimal ROUND_CEILING/ROUND_FLOOR plus an assert ✓. Every constant quoted in HD §2–9 equals the rerun's DOC line (beta 0.0533474↑,
  alo 44.9466↓, c3 7.63354↑, φ2 1.01113↑, cos 2β 0.994313↓, b −1.2224↑, g 0.61746↓, ηG 0.053287↑, tQ 1.9116↓, C_loc 2.0874↑, a, b, ϑ,
  sum, ce5 12.769↑, Cfin 70.600↑, 1 + lamL 1.05502↑, K_res 6416.59↑) ✓. Cosmetic: HD l. 85 says |Q5| ≤ 14.968 on |t| ≤ 1/(a* − 0.07);
  HC l. 113–116 computes it on the slightly smaller |t| ≤ 1/(alo − 0.01); the number is only reported (ce5 evaluates Q5 per cell
  on its own balls), so nothing depends on it.
- Holland's identities (factorization l. 1676–1680, V_m expansion, 3F2 ODE, Jacobi diagonal difference l. 1421): check_identities.py and
  check_jacobi.py rerun from copies: all True / max error 7·10⁻⁴⁹ ✓. Not re-audited (outside the brief): Laguerre defect, Jacobi model,
  four-coefficient matching (Brouwer replacement of the 1/30 contraction) — they had the orchestrator reviewer's CORRECT.

## 6. Reruns (logs/reruns.txt; copies in .scratch/jwaudit/runs, deleted after)

HollandK certify.py 90 2.6 1000 1.4 and 80 2.6 1000 1.4; SavvaK certify6.py 80 2.3 1000 1.33 and 24 4 400 1.33 -self: all four logs and both
witness CSVs byte-identical to the committed files; 7 s, 6 s, 37 s, 134 s; peak 85 MB RSS / 65 MB (memcap). SavvaK chk_alg, chk_rec,
chk_jacobi6: output identical to chk.log. Final lines: K_res ≤ 23.0332, K_N(9e24) ≤ 5.6413, "K <= 24: PASS"; K_SELF ≤ 2.756·10^25 PASS;
K ≤ 6416.59 ≤ 7000; K_SELF(B) ≤ 2.430·10^114.

## 7. All-n input d ≤ 9.0·10^24 (G l. 310–326, 748–777)

Source: Griffin–Ono–Rolen–Thorner–Tripp–Wagner, Adv. Math. 397 (2022) 108186, Theorem 1.2 (G l. 312–315): RH_m(T) and d ≤ ⌊T⌋² ⇒
J^{d,n} hyperbolic for all n ≥ m; m = 0 is Chasse's Theorem 1.8. Proof (G l. 772–776): zeros of ψ(√z), ψ(w) = ξ(½ + w), lie in
S(π, 2 arctan(1/(2T))) (on-line zeros up to height T map to the negative axis; |Re w| < ½ for every zero, so higher zeros have
|arg z − π| ≤ 2 arctan(1/(2T))); Lemma (G l. 762–771, Chasse Thm 3.6, positive Taylor coefficients γ(j) > 0) gives hyperbolicity for
d ≤ |sin δ|⁻² = T² + ½ + 1/(16T²) (sin(2 arctan x) = 2x/(1 + x²)) ✓. Input: Platt–Trudgian (Bull. LMS 53 (2021)) RH_0(3,000,175,332,800).
⌊T⌋² = 9.0010520·10^24 ≥ 9.0·10^24 ✓. Remarks: G's sentence "γ(j) > 0 ⇒ zeros in Re z < 0" is not needed for m = 0 (the S(0,δ) case
cannot occur since 0 < Re ρ < 1); the class C(θ,δ) should admit σ = 0 (G l. 752 writes σ ∈ S ⊂ C^×) — a definitional slip, harmless.
The alternative D_all = d0 = 1.6141·10^25 rests on a Lean-checked, unrefereed half-plane (COROLLARIES) and was not audited.

## 8. Reviewer on the most delicate step (§1 A2–A3, pasted inline with all formulas)

reviewer (gpt-6-astra, xhigh), verdict-only: **CORRECT** — horizontal ray legitimate (Re(z/(u0+s)²) > 0 throughout); absorption valid and
∫|s|³e^{−κ0s²/4} = 16/κ0²; concavity gives both tangent-line tails including the finite left tail; amplitude covers the whole ray;
evaluate at α ≥ a* − β = 44.94665 (HC does: alo), core + amplitude ≈ 5.322·10⁻²⁰; S6 covers both remaining pieces; 4δ means 4 sup|δ|
(HC uses 4δA). No missing term or uniformity loss. (logs/reviewer.txt)


## 9. Additional checks and findings (none load-bearing)

- MMP statements checked in the TeX of arXiv 2309.10970v3 (l. 711–719, 905–908): (iii) p, q ∈ P_n(R≥0) ⇒ p ⊠_n q ∈ P(R≥0) (Szegő–Walsh);
  p, q ∈ P_n(R>0) ⇒ lmesh(p ⊠ q) ≥ lmesh(p). Hypotheses hold (both factors: d simple positive roots, certified A − B − d > 0,
  C − D − d > 0, B, D > 0). MSS = Marcus–Spielman–Srivastava, PTRF 182 (2022), Thm 1.6 (Savva references.bib).
- Witness files (csv_check.py, logs/csv_check.log): cells_main.csv 1000 cells, cells_self.csv 2210 bisected cells, each covers
  [0, 1/L*] with no gap; maxima Cfin 18.7514 (cell 999) and 622.436, radius budgets 0.97277 and 0.545894 < 1, as logged.
- SavvaK chk_alg/chk_rec/chk_jacobi6 and HollandK check_identities/check_jacobi rerun from copies: outputs identical / all True.
- Misnaming: discoveries/jensen-two-thirds-wedge.md (Method, Caveat) and the brief say "Mirabelli–Marcus–Peters"; the cited source is
  Martínez-Finkelshtein–Morales–Perales, IMRN 2024 (16) 11642–11687.
- Rounding: discoveries/jensen-wedge-constant.md says "Kres < 6416.68 even from rounded displays"; the rounded displays give
  (70.600/0.929840)²·1.05502² = 6416.74 (certified 6416.59; K = 7000 unaffected).
- HD l. 85 quotes the |Q5| polydisc as |t| ≤ 1/(a* − 0.07); HC l. 113 uses 1/(alo − 0.01) (a smaller disc); reported value only.
- GORTTW l. 752 requires σ ∈ S(θ,δ) ⊂ C^× (σ = 0 must be allowed) and l. 775's "γ(j) > 0 ⇒ Re z < 0" is superfluous for m = 0;
  both harmless for Thm 1.2 with m = 0.
- Not checked (outside the brief or cited): Holland's Laguerre defect, Jacobi model and four-coefficient matching blocks (HD §4–5);
  the |a| ≤ 1.48 and e^{2β+1}√φ2 + 1 prefactors are irrelevant against exponents ≤ −10^10; Savva's Lean development.

