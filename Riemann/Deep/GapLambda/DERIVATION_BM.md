# DERIVATION_BM: certified Λ via Bui-Milinovich's criterion with an optimised amplifier

Source: Bui-Milinovich, arXiv 1410.3635v3 (`main.tex`, Q. J. Math. 2018). "l." means a line of that file.
Λ = limsup (t_{n+1} − t_n)/(2π/log t_n) over the ordinates t_n of zeros on the critical line (l.53, l.61-67).

## 1. Taken as printed (published, refereed)
- **Wirtinger inequality** (Thm W, l.139-153): f(a) = f(b) = 0 ⇒ ∫_a^b|f|² ≤ ((b−a)/π)²∫_a^b|f'|².
- **Twisted fourth moment** of Bettin-Bui-Li-Radziwiłł (Thm 5.1, l.500-517; error O(T^{1/2+2ϑ+ε} + T^{3/4+ϑ+ε}),
  l.511): valid for Dirichlet polynomials of length T^ϑ with ϑ < 1/4 (l.246, l.259) and coefficients ≪ h^ε.
- **Lemma main** (l.370-384) and **Theorem smoothed** (l.300-313): for M(s) = Σ_{h≤y} d_r(h)P[h]h^{-s},
  y = T^ϑ, ϑ < 1/4, P[h] = P(log(y/h)/log y) (l.243-255), the mean values of |F^{(j)}|², j = 0, 1, equal
  c_j(v,κ)·A_{r+2}(log y)^{r²+4r}𝓛^{4+2j}T/(2(r²−1)!((r−1)!)⁴) + O(T𝓛^{(r+2)²+2j−1}) (unsmoothed form l.323-326),
  with c_j the 9-fold integral at l.306-312.
- **Summation step** (l.180-202): the end pieces are O(1) long and the integrals are ≫ T (l.194).

## 2. Re-derived here
- **Reduction (l.164-202, l.327-330).** F = e^{ivt𝓛}ζ(½+it)ζ(½+it+iκπ/𝓛)M(½+it) (l.169-171). Its real zeros
  contain t_n and t_n − κπ/𝓛. If Λ ≤ κ, consecutive t_n are at most (1+o(1))2πκ/𝓛 apart, so consecutive zeros
  of F are at most (1+o(1))κπ/𝓛 apart (l.180-183). Wirtinger on each gap and summing give
  ∫_T^{2T}|F|² ≤ (1+o(1))(κ/𝓛)²∫_T^{2T}|F'|². The prefactors of c0 and c1 differ only by 𝓛², hence
  h1 = c0/(κ²c1) (l.329), and h1 > 1 contradicts Λ ≤ κ, so Λ > κ.
- **c_j at the shifts.** c_0 = c(α,β) at α1 = α3 = β1 = β3 = 0, α2 = iκπ/𝓛, β2 = −iκπ/𝓛 (l.391). In c(α,β)
  (l.376-381), y = T^ϑ and these shifts give exactly the phase Φ of l.307. The operators Q(𝓛⁻¹Σ∂_α) and
  Q(𝓛⁻¹Σ∂_β) with Q(x) = v + x (l.402) each produce v − ϑ(x+x1+x2+x3+x4) − (1−ϑ(x1+x3))t1 − (1−ϑ(x2+x4))t2:
  the (α2−α1)- and (β2−β1)-terms cancel in the sum of derivatives. This gives the square (v − S)² in c1.
- **Any polynomial P is admissible.**
  - P is "a certain polynomial" with arbitrary coefficients b_j (l.249-250), and the proof uses only its Mellin
    form (l.252-255).
  - The coefficients d_r(h)P[h] are ≪ h^ε because P is bounded on [0,1], so BBLR applies.
  - |M(½+it)|² ≪ T^{ϑ+ε}, so |F^{(k)}|² ≪ T^{2/3+1/4+ε} = o(T^{1−ε}) for every P (Weyl's bound for ζ), as l.194 needs.
  - BM themselves use a non-integer-coefficient P (l.333). Here r = 1, so x^{r²−1}(x1x2x3x4)^{r−1} = 1.

## 3. Evaluation of c_j (all re-derived; implemented in `cert_bm.py`)
With a = 1−ϑ(x1+x3), b = 1−ϑ(x2+x4), s = at1 − bt2, d12 = x1−x2, d34 = x3−x4,
S = ϑ(x+Σx_i) + at1 + bt2, and the arguments of P
w1 = 1−x−x1−x2, w2 = 1−x−x3−x4 (BM l.311; also in c(α,β), l.381). w1 and w2 involve only x-variables, never
t1..t4, and lie in [0,1] on the polytope:
1. **t3, t4 in closed form.** Φ is linear in t3 and in t4, so ∫∫e^{iΦ}dt3dt4 = e^{iϑκπ(x2−x4)}E(ω3)E(ω4) with
   E(ω) = (e^{iω}−1)/(iω), ω3 = −κπ(ϑd34 + s), ω4 = κπ(ϑd12 + s). Since ω3ω4 = −(κπ)²(ϑd34+s)(ϑd12+s), the
   integrand factor (ϑd12+s)(ϑd34+s) cancels the denominator exactly:
   (ϑd12+s)(ϑd34+s)E(ω3)E(ω4) = (e^{iω3}−1)(e^{iω4}−1)/(κπ)². No division remains.
2. **Four terms.** e^{iϑκπ(x2−x4)}(e^{iω3}−1)(e^{iω4}−1) = A + B + C + D, with ω3+ω4 = ϑκπ(d12−d34):
   - A = e^{iϑκπ(x1−x3)} and B = e^{iϑκπ(x2−x4)} do not depend on t;
   - C = −e^{−iκπ(t1−t2)}e^{iℓ^C·x}, ℓ^C = κπϑ(t1, 1−t2, −(1−t1), −t2);
   - D = −e^{iκπ(t1−t2)}e^{iℓ^D·x}, ℓ^D = κπϑ(1−t1, t2, −t1, −(1−t2)).
   So c_j = (κπ)⁻²∫ ab(v−S)^{2j}P(w1)P(w2)(A+B+C+D) over the 5 x-variables and t1, t2 ∈ [0,1].
3. **t1, t2 for A and B** are exact: ab(v−S)^{2j} is a polynomial in t (exact rational t-moments).
   **For C and D** we use the n_t × n_t tensor Gauss-Legendre rule in (t1, t2) (Arb's `legendre_p_root`, balls;
   cert_bm.py l.165, sums at l.181-182). The remainder is ≤ 4Mρ^{−(2n_t−1)}/(ρ−1) per variable and term. The
   bounds are **summed**: over t1 and t2 by the telescoping I − Q = (I1−Q1)⊗I2 + Q1⊗(I2−Q2) (the Q1 weights are
   positive and sum to 1), and over the terms C and D. The total is 2·2·4Mρ^{−(2n_t−1)}/(ρ−1) (cert_bm.py
   l.183-192, `eb = [2*2*base, 2*2*base*Smax**2]` at l.192), added as a ball of that radius to Re and Im of c_j
   (l.198-199). This is the Chebyshev bound on the Bernstein ellipse E_ρ of
   [0,1], with semi-axes (ρ±1/ρ)/4. The analytic t-integrand is bounded by
   M ≤ e^{κπ b_ρ}·sup|ab(v−S)^{2j}|·‖P‖²_∞/20, with b_ρ = (ρ−1/ρ)/4, because:
   - Im Φ_C = −κπa·Im t1 (resp. κπb·Im t2) with 0 < a, b ≤ 1;
   - |v−S| ≤ |v| + 2ϑ + t_max + 1 with t_max = ½ + (ρ+1/ρ)/4;
   - P(w1)P(w2) does not depend on t, so it stays real for complex t1 or t2. Since w1, w2 ∈ [0,1],
     |P(w1)P(w2)| ≤ ‖P‖²_∞ with ‖P‖_∞ the sup over [0,1];
   - the x-polytope has volume 1/20, so ∫|P(w1)P(w2)| ≤ ‖P‖²_∞/20.
   ρ is chosen by a grid search.
4. **x-integrals for fixed real ℓ (|ℓ_i| ≤ κπϑ).** The x-polytope is x ∈ [0,1] with (x1,x2) and (x3,x4) each in
   the triangle of size R = 1 − x.
   - Expand e^{i(ℓ1x1+ℓ2x2)} and e^{i(ℓ3x3+ℓ4x4)} to total order N. For real z with |z| ≤ κπϑ the remainder is
     ≤ ε = (κπϑ)^{N+1}/(N+1)!; the product remainder ≤ 2ε+ε² is added to every result times Σ|q|·‖P‖²_∞/20.
   - Each monomial uses the Dirichlet integral ∫_{u1+u2≤R} u1^α u2^β (R−u1−u2)^c = α!β!c!/(α+β+c+2)! R^{α+β+c+2}.
   - The last step is exact: ∫_0^1 x^e(1−x)^m dx = e!m!/(e+m+1)!.
   - ‖P‖_∞ is bounded by the smaller of two rigorous bounds, both exact rationals (`sup_abs_bound`, used at
     cert_bm.py l.74): max |Bernstein coefficient|, and the max over 400 subintervals [m−h, m+h] of
     Σ_k |P^{(k)}(m)|h^k/k!.
5. **Arb balls throughout** (192 bits) from exact rationals. Checks: Im c0 and Im c1 contain 0 (the c_j are real),
   c0, c1 > 0, and H = c0 − κ²c1 > 0 for the whole ball.

## 4. Numerical consistency (not part of the proof)
- The ball values equal the earlier float quadrature (`search/bm/bmfloat.py`, t3 and t4 via E-functions):
  BM's P at κ = 3.18, v = 1.26, ϑ = 1/4 gives c0 = 1.29667198837·10⁻⁵ (float 1.296672·10⁻⁵) and h1 = 1.0021778
  (float 1.002178). The degree-8 P gives h1 = 1.004200 (float 1.004200).
- A Monte Carlo of the raw printed 9-fold integrand matches the float quadrature (`logs/bm/bmmc_raw_integrand.txt`).

## 5. Certified instance (cert_bm.py, logs/cert_bm.txt)
- **Inputs:** ϑ = 2499/10000 < 1/4, r = 1, κ = 4081/1250 (3.2648), v = 131063/100000, N = 30 (Taylor order),
  n_t = 20 (Gauss points per t-variable), Arb precision 192 bits.
- **P** = Σ_{m=0}^{14} p_mψ_m, where ψ_m(x) = Σ_{k=0}^{m} (−1)^{m+k} C(m,k) C(m+k,k) x^k is the shifted Legendre
  polynomial on [0,1] (cert_bm.py l.27-28) and p = (1, 2.295083, 5.711135, 2.4756, 13.800476, −2.126714, 25.66932,
  −11.008222, 38.228923, −19.964301, 45.41707, −22.792755, 40.630752, −14.895337, 21.791921) as exact decimals.
  The exact monomial coefficients are printed at logs/cert_bm.txt line 17 ("record P: monomial coefficients
  (exact) = …", from cert_bm.py l.247), beginning 258266243/1000000, −14106897609/500000, …
- **Result:** h1 ∈ [1.00005692466 ± 2.5·10⁻¹²] and H = c0 − κ²c1 ∈ [2.4580659·10⁻⁹ ± 4.5·10⁻¹⁷] > 0, hence Λ > 3.2648.
- **Same script:** BM's own P at ϑ = 1/4 (h1 = 1.00217781433, their l.337) and at ϑ = 2499/10000
  (h1 = 1.00201577199, i.e. Λ > 3.18 certified).
