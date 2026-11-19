# kimaudit FULLPROOF: the chain of jensen-kim-lee-wedge, step by step, with source lines

Sources: **D** = Riemann/Deep/KimEff/DERIVATION.md (line numbers), **KL** = Kim–Lee arXiv 2105.05386 TeX (kim_lee.tex, line
numbers), **C** = Riemann/Deep/KimEff/certify.py. ✓ = re-derived by hand by kimaudit (19 Nov 2026); [x.py] = own script here.

## 0. Statement (D l.16–27) vs page
Page items 1, 2, 3 and the a = 1/2 line are D (i), (ii), (iii) verbatim: same Z, u(n), n_*, R_n, D_n, a, constants ✓.

## 1. Inputs
- (I1) D l.36–40, Riemann's integral. Titchmarsh (2.16.1–2): Ξ(t) = 2∫₀^∞ Φ(u) cos(ut) du, Φ = 2 Σ_k (2π²k⁴e^{9u/2} − 3πk²e^{5u/2}) e^{−πk²e^{2u}};
  with z = it this is ξ(1/2 + z) = 4∫₀^∞ φ(t) cosh(zt) dt, φ_k = πk²e^{5t/2}(2πk²e^{2t} − 3)e^{−πk²e^{2t}} ✓. Then 8ξ(1/2+z) = 32 Σ_m Mt(m) z^{2m}/(2m)!,
  so γ(m) = 32 m! Mt(m)/(2m)! ✓; 8ξ(1/2+z) = (4z² − 1)Λ(1/2 + z) is the GORZ normalisation ✓. Numerical: [jensen_small.py] vs mpmath ξ.
- (I2) D l.13, l.41: a = 0.37338196661320752 = 0.87338196661320752 − 1/2. TREE.md l.11 (N7) lists the Lean-checked half-planes
  (QuasiRiemannTracker): 0.87338196661320752 (Bebel, 31 Oct), then B_2 + 2.2e-254, B_2 = (777 − √1081)/852, 0.8733468 (PR 137),
  θ₁₄₁ = 0.87334585… (PR 141); each is ≤ 0.87338196661320752, and zeros are symmetric about Re s = 1/2 (ρ ↦ 1 − ρ̄), so
  |Re ρ − 1/2| ≤ a ✓. Refereed as an input by HPAudit (TREE l.11); not re-audited here. a = 1/2: classical ✓.
- (I3) D l.42–45 = KL l.168–189 verbatim: S(δ) = {|Im z| ≤ δ|z|}; P, Q real polynomials, Z(P) ⊂ S(δ), Q hyperbolic, deg Q ≤ δ^{−2}
  ⇒ P(D)Q = Σ_k P^{(k)}(0)Q^{(k)}/k! hyperbolic (KL Thm 3, proof "see [6]" = Obreschkoff 1961; the 1961 paper itself not seen).
  Corollary (KL l.183–189): J(P;d) hyperbolic for d ≤ δ^{−2}, as J(P;d)(z) = z^d (P(D)z^d)(1/z) ✓ (P(D)z^d = Σ_k binom(d,k)P^{(k)}(0)z^{d−k}).
  KL's J(f;d) = Σ_k binom(d,k) f^{(k)}(0) z^k (KL l.111) ✓, so J(Θ^{(n)};d) = J^{d,n} ✓; Θ(w) = 8Ξ₀(−w) with KL's Ξ₀ ✓.
  Sharpness example D l.45 ✓ (P(D)z^d = z^{d−2}(z² − 2d cos t z + d(d−1)), discriminant 4d(1 − d sin²t)). Random test [obresch.py].
- (I4) D l.46–48: DLMF 10.21.15 (product, ν > −1), 10.16.9 (₀F₁(;ν+1;−z²/4) = Γ(ν+1)(z/2)^{−ν}J_ν(z)), 10.21.3 (ν ≤ j'_{ν,1} < j_{ν,1}, ν ≥ 0) ✓.
  Hence G(x) = ₀F₁(;b;x) = Π(1 + x/c_k), c_k = j²_{ν,k}/4, b = ν + 1 ✓; S1 = 1/b, S2 = 1/(b²(b+1)) ✓ (the Rayleigh sums
  Σ j^{−2} = 1/(4(ν+1)), Σ j^{−4} = 1/(16(ν+1)²(ν+2))); c_1 > ν²/4 ✓.
- (I5) D l.49–50: Hadamard (order 1/2 ⇒ genus 0), Hurwitz ✓ standard.

## 2. Parabola and Gauss–Lucas (D l.52–63) ✓
z = b0 + it (|b0| ≤ a), w = z²: (Im w)² = 4b0²t² ≤ 4a²t² = 4a²(b0² − Re w) ≤ 4a²(a² − Re w) ✓; P closed, convex (hypograph of a
concave function of Im w) ✓. Θ^{(n)} order 1/2, genus 0, Θ^{(n)}(0) = γ(n) > 0 ✓. Genus-0 Gauss–Lucas: f'/f = Σ 1/(w − w_j)
absolutely convergent; strict separation of w0 from the closed hull gives Re(e f'/f(w0)) > 0 ✓; induction on n ✓.

## 3. Obreschkoff step (D l.65–74) ✓ (same structure as KL's proof of Thm 1, KL l.192–211)
No zero in W_R ⇒ every zero has s = −Re w > R and |Im w|² ≤ 4a²(a² + s), |w| ≥ s ⇒ |Im w|/|w| ≤ 2a√(a² + s)/s, decreasing in s
(derivative −(s + 2a²)/(2s²√(a²+s))) ⇒ ≤ δ = 2a√(a² + R)/R ✓. Partial products P_k (finite conjugation-closed zero sets,
multiplicity kept) are real, Z(P_k) ⊂ S(δ), P_k → Θ^{(n)} locally uniformly; P_k^{(j)}(0) → γ(n+j) (Weierstrass); limit degree d,
leading coefficient γ(n+d) > 0 ⇒ Hurwitz ✓. δ^{−2} = R²/(4a²(R + a²)) = D(R) ✓. (R > 2a² is not needed: if δ ≥ 1 then D(R) ≤ 1.)

## 4. Mixture identity (D l.76–84) ✓
d^n/dw^n cosh(t√w) = t^{2n} F_n(t²w), F_n(v) = Σ_k (n+k)! v^k/(k!(2n+2k)!); (2m)! = 4^m m! (1/2)_m and (1/2)_{n+k} = (1/2)_n (n+1/2)_k
give (n+k)!/(2n+2k)! = (n!/(2n)!)/(4^k (n+1/2)_k) ✓, so F_n(v) = (n!/(2n)!) G(v/4) ✓; Fubini bound G(y) ≤ e^{y/b} (y ≥ 0) ✓.

## 5. Facts (D l.86–112) ✓
(5.1) d/dq log(φ_k/φ_1) = −6π(k²−1)/((2πk²q−3)(2πq−3)) − π(k²−1) < 0 ✓ ⇒ sup at t = 0 ⇒ φ_1 ≤ φ ≤ (1+ε0)φ_1 ✓; tail k ≥ 12:
(2πk²−3)/(2π−3) ≤ 1.92k² ≤ 2k², (π−3)k² − π − 4 log k ≥ 7.3 > 0 at k = 12 and increasing ✓ ⇒ ≤ 4e^{−432} ✓. (5.2) ✓.
(5.3) ✓ (coefficient comparison). (5.4) |1 + x/c| ≥ Re = 1 − X/c ✓. (5.5) |1+x/c|² = (1−X/c)²(1 + Y²/(c−X)²), log(1+y) ≤ y,
c_k − X ≥ c_k(1 − τ+) ✓. (5.6) −y − y²/(2(1−y)) ≤ log(1−y) ≤ −y − y²/2 on [0,1) ✓. (5.7) geometric series with
Σ_k c_k^{−m−1} ≤ S2 c_1^{1−m} ✓. (5.8) c/(c−X)² = Σ(m+1)X^m/c^{m+1}, Σ_{j≥0}(j+2)τ^j = (2−τ)/(1−τ)² ✓. (5.9) (b)_k ≥ b^k ✓.

## 6. Tilted log-density (D l.114–127) ✓
h'(u) = 0 ⟺ Z(u) = n; both h' identities ✓ (substitute 2πe^{2u} = 2n/u + 9/2). ψ' = h' + 2βe^{−2t}/(1 − βe^{−2t}) − (ts/2)g(X) ✓
(d/dX log G(−X) = −g(X), X' = ts/2). ψ'' ✓ (d/dt of 2q/(1−q) with q = βe^{−2t} is −4q/(1−q)²); −ψ'' ≥ 2n/t² + 4πe^{2t} − a²/(2b) ✓
(for s < 0: (s/2)g ≥ −a²S1/2). Domain: X(t) < c_1; for s > 0, T_s = 2√(c_1/s) ≥ ν/√R = (u+1)/√θ > (u+1)ρ = U (θU < 1) ✓.

## 7. Parameters (D l.129–136) ✓
U²R = θU ν² ⇒ X ≤ θU c_1 on t ≤ U, X ≤ θ c_1 on t ≤ u + 1 ✓; |v| ≤ vmax, |w| ≤ Wmax (s² ≤ max(R², a⁴) = R²) ✓; Σ from (5.8)
with τ ≤ θ ((2−τ)/(1−τ)² increasing) ✓; E+ from (5.5) ✓.

## 8. Zero-free region (D l.138–189) ✓
- 8.1 (l.142–150). Q1: −h'(u_R) ≥ 2πe^{2u}(e^{2/n} − 1) ≥ 4πe^{2u}/n = 4/u + 9/n ✓; tilt ≤ u_R a²/(2b) ≤ (u+1)a²/(2n) ✓;
  β-term ≤ 2βe^{−2u}/(1−β) ✓ ⇒ −u ψ'(u_R) ≥ M1 > 0 ✓. Q2: h'(u_L) ≥ (2n/u)(1 − e^{−2Δ}) ✓; (u_L R/2)S1 ≤ (2n/u)θ/4 ⟸ ν² ≤ nb,
  u u_L ≤ (u+1)² ✓; (u_L R/2)X_L S2/(1−θ) ≤ (2n/u)θ²/(16(1−θ)) ⟸ ν⁴ ≤ nb²(b+1), u u_L³ ≤ (u+1)⁴ ✓; slack 0.012089223 ✓ [mono.py].
  Strict concavity ⇒ unique max ũ ∈ [u_L, u_R], ψ'(ũ) = 0 ✓.
- 8.2 (l.152–154) ✓: I ⊂ [u_L − 2ε_η, u + 1]; −ψ'' ≥ κ0 there; ψ'(u1) ≥ 2ηκ0 = 2K/η, ψ'(u2) ≤ −2K/η ✓.
- 8.3 (l.156–163) ✓: d/dt log G(x(t)) = (tw/2)Σ 1/(c_k + x); Im(w/(c+x)) = (sY + v(c − X))/|c+x|² = vc/|c+x|² as sY = vX ✓;
  |arg drift| ≤ η(u+1)vmax Σ = Φ ✓; |∫_I f| ≥ cos Φ ∫_I |f| (Φ ≤ π/2) ✓; ∫_I e^ψ ≥ 2η(e^{ψ(u1)} + e^{ψ(u2)}) ✓. Majorant:
  vmax ≤ 2a(√R + a), √R ≤ √θ n/(u+1), Σ ≤ (1 + Ĉ)/n (ν² ≤ b(b+1)), κ_low ≤ κ0 (4πe^{2u} ≥ 4n/u, b ≥ n) ✓.
- 8.4 (l.165–168) ✓: tangent bounds of the concave ψ on (0, u1] and [u2, U] ⊂ (0, T_s); E+ majorant: U⁴vmax²S2/(32(1−θU)²) with
  S2 ≤ n^{−3}, R(u+1)⁴/n³ ≤ θ(u+1)²/n ✓.
- 8.5 (l.170–186) ✓ + reviewer (below): (5.9) global bound, Hc concave on [U, ∞) (4πe^{2t} ≥ 4n/u), Hc'(U) ≤ −Hcp_low ≤ −1;
  lower bound of e^{ψ(u2)} via (5.6) with X(u2) ≤ X2 = θν²/4 (u2 ≤ u + 1) ✓; h(U) − h(u), h(u) − h(u2) from the h' identities
  (h concave ⇒ min over u2 ∈ [u_L, u + e] at an endpoint) ✓; U²Wmax/(4bn) ≤ θU/4 + aρ²(u+1)(√θ + a(u+1)/n)/(2n) ✓;
  −LG/n ≤ θ/4 + θ²/(32(1−θ)) ✓; −log η ≤ log(4πe^{2u}/K)/2 ≤ log(5n/K)/2 (4n/u + 9 ≤ 5n) ✓; dh1 + dh2 ≥ max ✓.
- 8.6 (l.188–189) ✓: |J| ≥ η(e^{ψ(u1)} + e^{ψ(u2)})[2cos Φ − (1+ε0)e^{E+}/(2K) − e^{logratio}], using e^{ψ(u2)} ≤ e^{ψ(u1)} + e^{ψ(u2)};
  the four pieces [0,u1], I, [u2,U], [U,∞) cover (0,∞) (u2 ≤ u+1 < U) ✓.


## 9. Uniformity in n ≥ n_* (D l.191–202) ✓
Every majorant is a function of u alone with n = Z(u) (C l.58–86) ✓; u(n_*) = 12.000000000000044 ≥ 12 [indep.py] ✓. Directions ✓:
(u+1)^k/n ↓ (k ≤ 4), n/u ↑; F' = ([2(ρ−1)u − 1]e^{2y} + 1 + 2ρ)/u² > 0 for u ≥ 1/(2(ρ−1)) = 11.4947 ✓; e = 1/n + 2η_up uses the
upper bound η_up (C l.73), and dh2 is increasing in e ✓; d(logratio)/dn ≤ 1/(2n) − c6 < 0 ✓. [mono.py]: 5766 grid points
u ∈ [12, 300] ∪ {400, 600, 10³, 3·10³, 10⁴}, both a: no monotonicity violation, every minimum at u = 12 (not a proof; the proof is §9).

## 10. (i) ⇒ (ii)–(iii) (D l.204–221) ✓
- D(R) ≥ R/(4a²) − 1/4 ⟺ R² ≥ R² − a⁴ ✓. u(n) + 1 ≤ log(n)/2 ⟸ (L/2 − 1)(π/e² − 9/(4n)) ≥ 1; its value at n_* is 5.4485
  (C l.122 prints value − 1 = H2 ≥ 4.4484549; D l.207 "(≥ 4.448)" means this margin: wording only) ✓. [indep.py]: max over
  n ∈ [n_*, 10^297] of u(n) + 1 − log(n)/2 = −0.8148 ✓.
- 3.6/(4a²) = 6.4555883891 ✓; slack (κ − 6.45)n² − κn − log²n/4 > 0 at n_*, increasing ✓. [indep.py]: min over n ∈ [n_*, 10^297]
  of D_n/(6.45 n²/log²n) = 1.0155 (a = 1/2, 3.59: 1.0175) ✓.
- Power forms ✓: recomputed sup f = 7.6889273, 0.8715955, 0.06982881, 1.1852008e-10 (a = 1/2: 10.36676, 1.2030145, 0.09924618,
  2.1293998e-10), each ≤ the printed constant; direct sampling of 2404 points (c, 1 ≤ d ≤ 10^60), n = ⌈max(n_*, C d^c)⌉: D_n ≥ d in
  every case (min D_n/d = 1.097; a = 1/2: 1.099) ✓.
- Comparisons (page "Before"): D_n = D_all at n = 46,974,172,465,412.6 (integer threshold 46,974,172,465,413 ≤ 4.6974173e13) ✓;
  D_n(10^13) = 8.1088e23; D_n − D_all = −3.09e17 at 4.6974172e13, +3.55e17 at 4.6974173e13 ✓; two-thirds wedge at d = D_all + 1:
  n ≥ 3.4172529e37 ✓; n = 10^30: D_n = 1.4616e57 vs 1.4225e20 ✓; d0 = 104,709,744.4 ≤ 1.04710e8 ✓; containment log-derivative
  4 − n/((n+2)log(n+2)) − 6/log n ≥ 4 − 1/log 10 − 6/log 8 = 0.68 ✓.

## 11. Code audit of certify.py
- Directed rounding: up()/lo() (C l.22–33) step a decimal until the Arb comparison with the ball endpoint holds ✓; checks are
  whole-ball strict (C l.17–20) ✓; 256-bit balls ✓.
- Formulas (C l.62–85, l.104) equal the D formulas term by term ✓; independently re-implemented from D's text in mono.py: same
  u = 12 values to 10 digits (η_up 6.856146951e-6, Φ_up 1.250725315e-4, E_up 7.86781902e-9, c6 0.7858282718, FINAL 1.874727976).
- Bounds certified: Laplace/peak (M1, Q2, M3a/b), phase with the largest |Im w| (vmax = 2a√(R + a²) enters Φ_up), side tails (E_up),
  crude tail (hcc, Hcp_low, c6, logratio), FINAL ✓. Not in code, analytic, checked here: the k ≥ 12 tail, 4n/u + 9 ≤ 5n, §9 directions.
- Thresholds (lesson 277): the proof of (i) has no scan: one evaluation at u_* = 12 ≤ u(n_*) plus monotone continuation ✓; H2, H3 at
  N0 = n_* exactly ✓; u_upper (C l.128–135) certifies Z(hi) > n, so R, D are lower bounds ✓; crossover bisection keeps a certified
  D > D_all point (C l.145–150), D_n increasing ✓. C l.157 prints a non-certified side of the two-thirds bisection (comparison only).
- Rerun on a copy in /tmp: logs byte-identical to logs/certify.log and logs/certify_half.log, "ALL CHECKS PASS" (logs/*_rerun.log).


## 12. Falsification (own mpmath code, not part of the proof; logs/)
- T1 [indep.py]: Z(12) = 998,616,377,770.9088…, ⌈Z(12)⌉ = 998,616,377,771 ✓; u(n_*) = 12 + 4.4·10⁻¹⁴ ≥ 12 ✓.
- T2 [indep.py]: D_n(10¹³) = 8.1088·10²³ (< D_all); D_n(4.6974173·10¹³) = D_all + 3.55·10¹⁷, D_n(4.6974172·10¹³) = D_all − 3.09·10¹⁷;
  exact crossover n = 46,974,172,465,412.6 ✓ (stated threshold rounded up correctly).
- T3 [indep.py]: the four (c, C) pairs for both a, by closed form and by a 40,000-point grid sup, plus 2404-point direct sampling ✓ (§10).
- T5 [mono.py]: §9 directions on 5766 grid points ✓.
- T7 [jensen_small.py]: γ(0), γ(1), γ(2) from the moments of φ equal 8 m! ξ^{(2m)}(1/2)/(2m)! (mpmath ξ; ratio − 1 ≤ 10⁻⁴¹) ✓;
  Taylor Θ(−50 + 10i) = 8ξ(1/2 + √w) to 12 digits; Θ(−γ₁²) ≈ 0 (2·10⁻¹⁰⁸) ✓. J^{d,n} hyperbolic, certified by d sign alternations
  at points between the polyroots approximations (110 digits), for d = ⌊D_n⌋, ⌊D_n⌋ + 1, 2⌊D_n⌋ at n = 3 (D_n = 4.31), 5 (12.45),
  8 (30.35), 10 (45.59): 12/12 hyperbolic, inside and outside (D_n is used below its range n ≥ n_*; outside holds too, as these d
  lie far inside the all-n range D_all ≈ 1.6·10²⁵).
- T8 [obresch.py]: KL Thm 3 (= (I3)) on 800 random instances (δ ∈ [0.3, 0.75]; P: 1–3 conjugate pairs mostly on ∂S(δ), both
  half-planes, plus real zeros; Q: random hyperbolic, deg Q = ⌊δ^{−2}⌋): 0 non-hyperbolic P(D)Q ✓; control deg Q = ⌊δ^{−2}⌋ + 1:
  173/800 non-hyperbolic (the test detects failures); sharp example P = z² − 2z cos 1 + 1: d = 1 hyperbolic, d = 2 not ✓.
  (A first run had my own bug, P^{(k)}(0)/k! divided by k! again; fixed before the logged run.)

- T6 [theta_n.py] zero-free parabola at moderate n, from the Taylor series Θ^{(n)}(w)/γ(n) = Σ_m γ(n+m) w^m/(γ(n) m!) (γ from
  Gauss–Legendre moments of φ, panels of width 0.1, degree 7, degree-6 cross-check; Horner), not from the mixture identity; zeros
  counted by the argument principle on ∂W (upper half by conjugate symmetry, adaptive steps ≤ 0.3 rad). (i) is claimed only for
  n ≥ n_*; these runs test the same formula R_n far below that range.
  - n = 60 (60 digits, M = 200; moments agree to 10⁻⁶⁰): W_{R_n} (R_n = 577.15, D_n = 1034.7) has 0 zeros; Θ^{(60)}(−s) > 0 on
    [0, R_n]; first zero on the negative axis at s1 = 4090.65 = 7.09 R_n; W_{0.999 s1}: 0 zeros; control W_{1.001 s1}: 1 zero ✓.
  - n = 1000 (300 digits, M = 950; moments degree 7 vs 6 agree to 6.4·10⁻¹⁴⁸): W_{R_n} (R_n = 76,042.8, vmax = 205.9) has 0 zeros
    (max arg step 0.001 rad; cancellation 10⁹⁷ < 10¹⁴⁷, truncation 10⁻⁶¹²: reliable); Θ^{(1000)}(−s) > 0 on [0, R_n] (min 7.4·10⁻⁵²).
    W_{2R_n}: 0 zeros (cancellation 10¹⁹⁸: relies on degree-7 moments beyond the 147 checked digits; indicative). W_{3R_n} and the
    sign scan near 3R_n: not reliable (truncation tail 10⁻¹¹⁵ exceeds |Θ^{(n)}/γ(n)| ≈ 10⁻¹¹⁸ there); not used.
  - Reading: no zero violates item 1's region at these n; the true zero-free reach is several R_n (θ = 9/10 and the tilted peak
    make R_n a method limit, as the page says).

- T9 [mixture.py]: D §4 at n = 60 (30 digits): Taylor series = ∫ p_n(t) ₀F₁(; n + 1/2; t²w/4) dt (Gauss–Legendre, mpmath hyp0f1)
  at the corner w = −R_n + i vmax (rel. diff 3·10⁻²⁸) and at w = −R_n/2 + 0.7 i v(R_n/2) (5·10⁻²⁹) ✓.

## 13. Reviewer: see REPORT.md (run 1 on §8.5: (A)–(H) correct, GAP only from my prompt; run 2 closure: CORRECT).
