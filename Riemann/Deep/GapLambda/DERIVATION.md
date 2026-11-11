# DERIVATION: from Bui-Hall-Subirà Jorge's main terms to a certified bound for Λ

Source: arXiv 2511.14415v1, file `BuiHallSubiraJorge_amplification.tex`. "l." means a line of that file.
Λ = limsup (t_{n+1} − t_n)/(2π/log t_n), t_n the ordinates of the zeros on the critical line (l.102, l.110).

## 1. Reduction (their §3, re-checked)
- Wirtinger-type inequality (Thm Wi, l.507-516; λ0 at l.514): for f(a) = f(b) = 0 and v ≥ 0,
  3λ0(v)∫|f|⁴ ≤ ((b−a)/π)⁴∫|f'|⁴ + 6v((b−a)/π)²∫|ff'|², with λ0(v) = (1 + 4v + √(1+8v))/8.
- f(t) = e^{iu(log T)t} ζ(½+it) A(½+it) (l.554). Its real zeros include every t_n, so Λ ≤ κ gives gaps
  b − a ≤ (1+o(1))2πκ/log T between consecutive zeros of f in [T, 2T]. Then ((b−a)/π)⁴ ≤ (1+o(1))16κ⁴/L⁴ and
  6v((b−a)/π)² ≤ (1+o(1))24vκ²/L², L = log T. Summing over the intervals (end pieces are O(1) long, the
  integrals are ≫ T; same summation as l.450-462) gives keyineq (l.575-578): if 3λ0∫|f|⁴ > 16κ⁴L⁻⁴∫|f'|⁴ + 24vκ²L⁻²∫|ff'|², then Λ > κ.
- Main terms (Thm 1.1, l.75-100, needs ϑ < 1/8, l.76; r ∈ ℕ): ∫|f|⁴, L⁻⁴∫|f'|⁴, L⁻²∫|ff'|² have main terms
  K·c(0,0), K·c1(u), K·c2(u) with the same K = a_{2r+2}T(log y)^{4r²+8r}L⁴/(2((r−1)!)⁸((r²−1)!)⁴)
  (l.589, l.610, l.644). K > 0 (a_{2r+2} is an Euler product of positive factors, l.82) cancels: Λ > κ if 3λ0(v)c > 16κ⁴c1(u) + 24vκ²c2(u) (keyineq2, l.666-668).
- Solving the quadratic in X = κ²: 16c1X² + 24vc2X − 3λ0c < 0 ⟺ X < (√(144v²c2² + 48λ0c·c1) − 12vc2)/(16c1)
  (l.671, keyineq3 l.674-676). Constants 3, 16, 24, 144, 48, 12, 16: rg-checked in `cert.py`.

## 2. c1, c2 from c(α, β) (re-derived, matches l.615-663)
f'/(ie^{iuLt}) = L·Q(L⁻¹(∂_α + ∂_β)) ζ(½+α+it)A(½+β+it)|₀ with Q(x) = u + x (l.585-600). Q(L⁻¹(∂_α+∂_β)) acting
on T^{aα}T^{bβ} multiplies by u + a + b. With y = T^ϑ, E12 = v1 − v2 + ϑ(x1+x2−z1−z2−v1Σx+v2Σz),
E34 = same with x3+x4−z3−z4, the T-exponents of α_j, β_j in the integrand of c (l.86-99) are
- α1: −ϑ(x1+x2) − v1(1−ϑΣx) + v3E12; β1: −ϑ(t1+t2+x3+z3)
- α2: −ϑ(z1+z2) − v2(1−ϑΣz) − v3E12; β2: −ϑ(t3+t4+x4+z4)
- α3: −ϑ(x3+x4) − v1(1−ϑΣx) + v4E34; β3: −ϑ(t1+t3+x1+z1)
- α4: −ϑ(z3+z4) − v2(1−ϑΣz) − v4E34; β4: −ϑ(t2+t4+x2+z2)
so c1 inserts F1F2F3F4 and c2 inserts F1F3 (pairs (α1,β1), (α3,β3): one +it and one −it factor), with
F1 = u − ϑ(t1+t2+x1+x2+x3+z3) − v1(1−ϑΣx) + v3E12, F2 = u − ϑ(t3+t4+x4+z1+z2+z4) − v2(1−ϑΣz) − v3E12,
F3 = u − ϑ(t1+t3+x1+x3+x4+z1) − v1(1−ϑΣx) + v4E34, F4 = u − ϑ(t2+t4+x2+z2+z3+z4) − v2(1−ϑΣz) − v4E34.
All four agree with the printed c1 (l.626-640); F1, F3 with c2 (l.653-663). Code: `core.poly_parts`.

## 3. Exact evaluation of the 16-fold integrals (`core.py`)
- Slacks w1 = 1−t1−t2−x3−z3, w2 = 1−t3−t4−x4−z4, w3 = 1−t1−t3−x1−z1, w4 = 1−t2−t4−x2−z2 (the P arguments).
  Each constraint holds one (x, z) pair: (x3,z3) in the triangle of size R1 = 1−t1−t2, (x4,z4) in R2 = 1−t3−t4,
  (x1,z1) in R3 = 1−t1−t3, (x2,z2) in R4 = 1−t2−t4. Since 1 − w1 = t1+t2+x3+z3 etc., the integrand without
  P and weights is a polynomial in X12 = x1+x2, Z12, X34, Z34, w1..w4, u, v1..v4.
- v1..v4 ∈ [0,1] are integrated exactly (rational polynomial arithmetic, `fmpq_mpoly`).
- Each monomial is expanded binomially in (x1,x2), (z1,z2), (x3,x4), (z3,z4); triangle integrals by Dirichlet:
  ∫_{x,z≥0,x+z≤R} x^{p+r−1}z^{q+r−1}(R−x−z)^m = R^{p+q+m+2r}(p+r−1)!(q+r−1)!m!/(p+q+m+2r)! (`core.dirichlet`),
  weights (x_jz_j)^{r−1}, P(w) = Σ P_c w^c contributes m = b + c.
- t-integral G(n;e) = ∫(t1t2t3t4)^e R1^{n1}R2^{n2}R3^{n3}R4^{n4}, e = r²−1, over t ≥ 0, R_k ≥ 0. Split t3 ≤ t2
  (then R3 = R1+δ, R2 = R4+δ, δ = t2−t3 ∈ [0,t2], and R1, R4 ≥ 0 imply the other two constraints): expand
  R3^{n3}, R2^{n2} binomially, integrate t1 ∈ [0,1−t2], t4 ∈ [0,1−t2], δ, then t2:
  H(n) = Σ_{j≤n3,k≤n2} C(n3,j)C(n2,k) B(n1+j+1,e+1) B(n4+k+1,e+1) B(Q+1,e+1) B(P+1,2e+Q+2),
  Q = n3−j+n2−k, P = n1+j+n4+k+2e+2. The swap t2 ↔ t3 maps R1↔R3, R2↔R4, so G(n) = H(n1,n2,n3,n4) + H(n3,n4,n1,n2)
  (`core.Hcase`, `core.G`; 6 cases checked against sympy, `logs/test1.txt`). `core.G_exact_table` regroups the
  same terms with j' = n3−j, k' = n2−k, N = Σn_i: H = Σ u_{n1,n3}[j'] γ_N(j'+k') u_{n4,n2}[k'],
  u_{a,c}[j'] = C(c,j')B(a+c−j'+1,e+1), γ_N(s) = B(s+1,e+1)B(N−s+2e+3,2e+s+2) (Hankel form; 50 random entries equal
  `core.G` exactly, `logs/test4.txt`).

## 4. Checks of the evaluation
- ϑ = 0, r = 1, P = 1, u = 1/2: c/c1 = 560, c2/c1 = 28/3 exactly, i.e. Λ² ≥ √(49v²+105λ0(v)) − 7v (l.680);
  v = 22/49 gives λ0 = 121/196 and Λ² ≥ 11/2 exactly (l.682). `cert.py` asserts this to 10⁻⁶⁰.
- Independent Monte Carlo of the raw 16-dim integrand as printed (l.86-99, l.626-640, l.653-663), ϑ = 0.1249,
  P = 1 − 3x + 1.5x², u = 0.6: 7.5·10⁶ and 3·10⁷ samples, |z| ≤ 2.07 and ≤ 1.21 for c, c1, c2 (`logs/mc_seed*.txt`).
- ϑ = 0.1249, r = 1, P = 1, u = 0.6, v = 0.4 (their choice, l.684): Λ ≥ 2.639532 < 2.64; tuned u = 0.5942,
  v = 0.4038: Λ ≥ 2.641023, certified Λ > 2.64. Their Theorem 1.2 holds; the printed (u, v) is slightly off.

## 5. Certificate (`cert.py`)
Inputs are exact rationals: ϑ = 1249/10000 (< 1/8), r = 1, u = 843/1250 (0.6744), v = 40741/100000, κ = 2907/1000,
P = Σ_m p_m ψ_m(x) with ψ_m the shifted Legendre polynomials on [0,1] (exact integer coefficients
(−1)^{m+k}C(m,k)C(m+k,k)) and p = (1, 3.2279, 4.6423, 9.7348, 5.6966, 21.5327, 1.8217, 36.1817, −4.141,
42.4179, −5.2206, 28.2819) as exact decimals (degree 11); the exact monomial P is printed in `logs/cert.txt`. The polynomial integrand and
all Dirichlet coefficients and G values are exact rationals; they are converted to Arb balls (320 bits) and the
matrix products and the final inequality run in ball arithmetic (outward rounding). Arb's `H > 0` is true only if
the whole ball is positive. Result: H = 3λ0c − 16κ⁴c1 − 24vκ²c2 ∈ [1.0192444·10⁻¹² ± 4.4·10⁻²⁰] > 0, so Λ > 2.907;
κ_max = 2.9071547… (ball [2.90715470829085 ± 1.5·10⁻¹⁵]).

## 6. Untrusted search (`search/`)
P enters c, c1, c2 as quartic forms (|A|⁴ contains P four times), so the quotient is a ratio of quartic forms,
not a generalised eigenproblem. Search: multi-start BFGS over (Legendre coefficients of P on the unit sphere,
u, v). The monomial-basis float evaluation is ill-conditioned (degree 5: float 3.22 vs exact 2.786,
`logs/exact_d5_monomial_float.txt`); `fast2.FastL` uses Legendre bases in w and in each R_k, with the basis
changes done exactly (rationals, and a 384-bit ball transform of the G table), and agrees with the exact values to
8 digits. A guard rejects points where c or c1 is below 10⁻⁹ of its absolute-value bound.
