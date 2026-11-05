# Diagonal Fasi–Higham positivity: proof

**Theorem.** Let m ≥ 1 be odd and let q(z) = q_{m,m}(z) = Σ_{j=0}^{m} (−1)^j C(m,j) (2m−j)!/(2m)! z^j be the
denominator of the diagonal [m/m] Padé approximant to e^z. Then every Taylor coefficient c_n of 1/q(z) is
positive. This is the diagonal case k = m of Fasi–Higham, SIMAX 40 (2019), Conjecture 3.2.

The proof is computer-assisted. Every finite inequality is checked in Arb ball arithmetic (python-flint), and
`verify.py` re-runs every check (§11).

- m = 1: q = 1 − z/2, so c_n = 2^{−n}.
- 3 ≤ m ≤ 201: certified computation for each m (§9).
- m ≥ 203: analytic argument (§1–§8). It uses explicit constants and bounds that hold uniformly for
  ν = m + ½ ≥ ν0 = 203.5.

Notation. ν = m + ½ and z = 2νζ. D = ζ d/dζ. a = 2l + 1 for l = 0, 1, 2, … and α = 1/a.
y_J = 0.85, and log, √ are principal branches.

- s(ζ) = √(1+ζ²)
- Ψ_e(ζ) = s − log(1+s)
- η(ζ) = Ψ_e(ζ) + log ζ
- P(ζ) = ¼ log(1+ζ²)
- A_a(ζ) = ζ + aΨ_e(ζ)
- U1(p) = (3p − 5p³)/24, with p = 1/s

All of these are analytic on Ω = {Re ζ > 0} ∪ {|ζ| < 1}. The branch cuts of s are ±i[1, ∞), and Re s > 0, so 1 + s ≠ 0.

## 1. Exact identities

**(R) Padé remainder.** Let p(z) = q(−z), so p has positive coefficients. Let Φ(t) = t^m(1−t)^m/(2m)!.

1. Integrate ∫_0^1 e^{tz}Φ(t)dt by parts 2m+1 times. Φ^{(k)}(0) = Φ^{(k)}(1) = 0 for k < m.
2. By the symmetry Φ(1−t) = Φ(t), Φ^{(k)}(1) = (−1)^kΦ^{(k)}(0), and Φ^{(m+i)}(0) = (−1)^i C(m,i)(m+i)!/(2m)!.
3. This gives z^{2m+1}∫_0^1 e^{tz}Φ = (−1)^m(e^z q(z) − p(z)).

Hence, for odd m,

  e^z q(z) = p(z) − z^{2m+1}E(z),  E(z) = Σ_j E_j z^j,  E_j = m!(m+j)!/(j!(2m)!(2m+j+1)!) > 0.

`identities.py` also checks this exactly as power series for m ≤ 11.

**(I1) Bessel form of p.** Put x = z/2. Then p(z) = θ_m(x)/θ_m(0) = √(2/π) x^ν e^x K_ν(x)/θ_m(0), with
θ_m(0) = (2m)!/(m! 2^m). This is the reverse Bessel polynomial (DLMF 10.49.12).

p has no zeros in Re z ≥ 0, because K_ν has none in |ph x| ≤ π/2 (DLMF 10.42). Hence q has none in Re z ≤ 0.

**(I2) Bessel form of X.** Let X := z^{2m+1}E/p. Then X = πI_ν(x)/K_ν(x).

- Proof: use K_ν(xe^{iπ}) = e^{−iπν}K_ν(x) − iπI_ν(x) (DLMF 10.34.2) and e^{iπν} = (−1)^m i = −i for m odd.
- This gives θ_m(−x) = √(2/π)x^ν e^{−x}(K_ν − πI_ν)(x).
- So e^z q/p = e^{2x}θ_m(−x)/θ_m(x) = 1 − πI_ν/K_ν, and X follows from (R).
- Checked numerically in `identities.py` (relative error 5e−54).

**(LG) Liouville–Green form.** Define h_K and h_I on Re ζ > 0 by

  K_ν(νζ) = √(π/(2ν)) e^{−νη}(1+ζ²)^{−1/4}(1 + h_K),   I_ν(νζ) = κ (2πν)^{−1/2} e^{νη}(1+ζ²)^{−1/4}(1 + h_I),

where κ = √(2πν)ν^ν e^{−ν}/Γ(ν+1). DLMF 5.6.1 gives e^{−1/(12ν)} < κ < 1.

With z = 2νζ, (I1) and (I2) become

  e^z/p(z) = Λ e^{ν(ζ+Ψ_e)}(1+ζ²)^{1/4}/(1+h_K),   X = κ e^{2νη}(1+h_I)/(1+h_K),   Λ = θ_m(0)√ν ν^{−ν}.

Consequences:

- **Disc.** Both 1 + h_K and 1 + h_I extend analytically to |ζ| < 1, and h_I is even.
- **Left half-plane.** For Re ζ < 0, 1/q(2νζ) = Λ e^{ν(ζ+Ψ_e(ζ))}(1+ζ²)^{1/4}/(1 + h_K(−ζ)).
- **Reflection.** For Re w > 0, 1 + h_K(−w) = (1 + h_K(w)) − κ e^{2νη(w)}(1 + h_I(w)). This is (R) at z = 2νw.

## 2. Liouville–Green error bounds (Lemma 1)

**Lemma 1.** For Re ζ > 0,

  |h_K(ζ)| ≤ exp(2V_K/ν) − 1,  where V_K = ∫|U1′(w)||dw| along a path from ζ to +∞;
  |h_I(ζ)| ≤ exp(2V_I/ν) − 1,  where V_I = ∫|U1′(w)||dw| along a path from 0 to ζ.

In both cases Re η must be nondecreasing along the path (in the direction from ζ to ∞, and from 0 to ζ).
Here |U1′(w)| = |3 − 15p²||w||p|³/24.

*Proof.* Write y(ζ) = K_ν(νζ) or I_ν(νζ). Then w = ζ^{1/2}y solves w″ = (ν²f + g)w with f = (1+ζ²)/ζ² and g = −1/(4ζ²).

1. The Liouville map ξ = η(ζ) (η′ = s/ζ = f^{1/2}) and W = f^{1/4}w = (1+ζ²)^{1/4}y give d²W/dξ² = (ν² + ψ)W.
   Here ψ = ζ²(4 − ζ²)/(4(1+ζ²)³).
2. Write W = c e^{∓νξ}(1 + h). With the upper sign (K), h solves the Volterra equation

     h(ξ) = ∫_ξ^{∞} (1 − e^{−2ν(v−ξ)})/(2ν) ψ(v)(1 + h(v)) dv.

   With the lower sign (I) the equation runs from −∞ (ζ = 0) and has kernel (1 − e^{−2ν(ξ−v)})/(2ν).
3. On a progressive path, Re(v − ξ) ≥ 0 (respectively Re(ξ − v) ≥ 0), so the kernel is bounded by 1/ν.
4. Picard iteration then gives |h| ≤ exp(∫|ψ||dξ|/ν) − 1. By sympy (`identities.py`), ψ·dξ/dζ = 2U1′(ζ).
5. The constants c are fixed by the limits:
   - K_ν(νζ) ~ √(π/(2νζ)) e^{−νζ} as ζ → +∞, with η = ζ + o(1);
   - I_ν(νζ) ~ (νζ/2)^ν/Γ(ν+1) as ζ → 0, with e^{νη} ~ e^ν(ζ/2)^ν.

This is Olver, *Asymptotics and Special Functions*, Ch. 6 Thm 11.1 and Ch. 10 §7, in the crude form with error-control function F = 2U1.

Existence of the solution, its identification with K_ν and I_ν, and independence of the path are the standard parts of Olver's theorem. ∎

Check: `identities.py` (8) compares the exact h_K and h_I (mpmath Bessel functions) with the `lg_box` bounds at 21 points for ν = 5.5, 21.5 and 203.5. The largest ratio |h|/bound is 0.495, which matches the factor 2 lost against Olver's sharper constant.

**Admissible paths** (`fh2lib.py`).

- *Rays from 0* (for I) and *rays to ∞* (for K): d/dρ Re η(ρe^{iθ}) = Re s/ρ > 0.
- *Rightward horizontal segments* in Re ζ > 0 (for K): Re(s/ζ) > 0, because arg s lies between 0 and arg ζ.
- *The imaginary segment* [0, iy], y < 1 (for I): Re η(iy) is increasing (§3).
- `lg_box(r0,r1,θ0,θ1,ν)` bounds V_I and V_K over a polar box with Riemann sums of ball suprema of |U1′|. It uses:
  - V_I: the ray from 0.
  - V_K for θ ≤ 1.2: the ray to ∞; the tail ρ ≥ 8 uses |p| ≤ (ρ²−1)^{−1/2}, contributing (3 + 15/63)/(24√63).
  - V_K for θ > 1.2: a horizontal segment to the ray of angle 1.2, then that ray.
  - The result is ε = exp(2V/ν0) − 1 and λ = −log(1 − ε).
- Since ε decreases in ν, every ε or λ computed at ν0 is valid for all ν ≥ ν0.
- On a box, |log(1+h)| ≤ λ and 1 + ε ≤ e^λ. Hence W_l := (1+h_I)^l(1+h_K)^{−l−1} satisfies

    |log W_l| ≤ λ_W := lλ_I + (l+1)λ_K ≤ a[(λ_I+λ_K)/2 + α|λ_K−λ_I|/2].

- `disc0_lg` (`goodsmall.py`) bounds h_I and h_K on the disc |ζ| ≤ 0.3 by the maximum principle on |ζ| = 0.3.
  - The right half-circle is bounded directly.
  - On the left half-circle, h_I is even. For h_K, use the reflection identity: on |w| = 0.3 with Re w ≥ 0, Re η(w) ≤ Ψ_e(0.3) + log 0.3 < 0.
- `assembly.py` caches λ_I and λ_K on a polar grid of 2448 cells (radii 0.004·1.04^k up to 1.75, 16 angle cells on [0, π/2]), rounding values up.
  - Any region is bounded by the maximum over all cells that meet it, padded by 1e−12.

## 3. Contour decomposition (Lemma 2)

Fix m ≥ 203 odd and n ≥ 0. Let N = N_l = n − l(2m+1) (an integer), t = N/ν and τ = n/ν = t + 2l.

**(F0) Monotonicity on the imaginary axis.** On [0, 1):

- Ψ_e(iy) = √(1−y²) − log(1+√(1−y²)) is decreasing, since its y-derivative is −y/(1+√(1−y²)).
- Re η(iy) = Ψ_e(iy) + log y is increasing, since its y-derivative is 1/y − y/(1+√(1−y²)) > 0.

**(X<1) X stays below 1 on the imaginary segment.** For |y| ≤ y_J, Lemma 1 and κ < 1 give

  |X(iy)| ≤ X_max := exp(2ν Re η(iy_J) + λ_I + λ_K) < 1,

with Re η(0.85i) = −0.0589. The check runs at ν0 (`neg_pieces` and `large_tau.py`: "X(i yJ) < 1"). Since Re η(iy_J) < 0, it holds for all ν ≥ ν0.

**Lemma 2.** (2ν)^n c_n/Λ = I_L + Σ_{l≥0} κ^l J_l, where:

- I_L = (1/2πi)∫ ζ^{−n−1} e^{ν(ζ+Ψ_e)}(1+ζ²)^{1/4}/(1+h_K(−ζ)) dζ, taken over the left half-circle |ζ| = y_J from iy_J to −iy_J.
- J_l = (1/2πi)∫_{P_l} e^{ν(A_a(ζ) − t log ζ)}(1+ζ²)^{1/4} W_l(ζ) dζ/ζ.
- P_l is any path from −iy_J to iy_J in {Re ζ ≥ 0, |Im ζ| ≤ y_J} that passes to the right of 0. If N_l < 0, P_l may be the segment [−iy_J, iy_J].

*Proof.*

1. c_n = (1/2πi)∮ z^{−n−1}/q over a small circle.
2. Deform that circle to the closed curve made of: the left half-circle of radius 2νy_J; the imaginary segments; and a tiny right half-circle around 0. The enclosed region contains no zero of q (by I1, the zeros have Re z > 0).
3. On the imaginary segments and the tiny half-circle, |X| < 1 (by (X<1), and X = O(z^{2m+1})). So 1/q = e^z/(p(1−X)) = (e^z/p)Σ_l X^l converges uniformly there, and the series can be integrated term by term.
4. The l-th term is e^z z^{l(2m+1)}E^l p^{−l−1} z^{−n−1}. It is analytic in Re z ≥ 0 except at 0, where there is no singularity if N_l < 0. So it can be deformed to P_l.
5. Inserting (LG) gives the stated integrands, because e^z X^l/p ∝ κ^l e^{ν(ζ + aΨ_e + 2l log ζ)} and ζ^{2lν} ζ^{−n−1} = ζ^{−N−1}. ∎

**Facts used repeatedly** (elementary; derivative identities are in `identities.py`).

- (F1) d/dφ Re Ψ_e(re^{iφ}) = −Im s = −r² sin2φ/(2 Re s) ≤ 0 on [0, π/2]. Hence Re A_a(re^{iφ}) ≤ A_a(r), and Re Ψ_e(re^{iφ}) − Ψ_e(r) = −r²∫_0^φ sin2ψ/(2Re s)dψ.
- (F3) For real ρ, νA_a(ρ) − N log ρ is convex in log ρ (D² = ν(ρ + aρ²/s) > 0). Its minimum is at the saddle r* = r*(a,t), where DA_a(r*) = r* + a(s(r*) − 1) = t. The minimum value is νΦ(a,t), with Φ = A_a(r*) − t log r*.
  - Closed forms: r* = t(t+2)/(2(t+1)) for a = 1, and r* = (−(t+a) + a√((t+a)² − a² + 1))/(a² − 1) for a ≥ 3.
- (F4) Φ(a,t) is concave in t, since dΦ/dt = −log r* and r* increases with t. So on a t-interval, Φ is minimal at an endpoint.
- (F5) D^kΨ_e = s−1, ζ²/s, ζ²(2+ζ²)/s³, ζ²(1/s+3/s⁵). D^kP = v/2, uv, −2uv(1−2u), 4uv(1−6v+6v²), with u = 1/(1+ζ²) and v = ζ²u. On (0, ∞), 0 ≤ DP < 1/2 and 0 ≤ D²P ≤ 1/4.
- (F6) Cauchy for D^k. If |F| ≤ M on a disc of radius δ around ζ, then |D^kF(ζ)| ≤ M c_k(|ζ|/δ). Here D^k = Σ_j S(k,j)ζ^j∂^j, and
  - c1 = ρ
  - c2 = ρ + 2ρ²
  - c3 = ρ + 6ρ² + 6ρ³
  - c4 = ρ + 14ρ² + 36ρ³ + 24ρ⁴

## 4. Good terms: a rigorous saddle-point lower bound (Lemma 3, `good.py`, `goodsmall.py`)

Write f = νA_a − N log ζ + P + log W_l. Then κ^l J_l = κ^l (1/2πi)∫_{P_l} e^f dζ/ζ. Note that f is real on (0, ∞), and W_l > 0 there.

**Lemma 3 (Good bound).** Take a box of saddle radii r* ∈ [r_lo, r_hi] and α ∈ [α_lo, α_hi], all ν ≥ ν0. Let `good_bound` return e < 0.9 for that box. Then for every (n, ν, l) with r*(a, t) in the box,

  κ^l J_l ≥ exp(νΦ(a,t) + G),  G = log(1−e) − λ_{W0} − l/(12ν) − ½ log(2πν a b_h),

where λ_{W0} bounds |log W_l| on [0.8r_lo, 1.25r_hi] ⊂ (0, ∞) and b_h is stored in the table. In particular J_l > 0.

*Proof.* There are six steps.

1. **Exact critical point.**
   - Let F(ρ) = Df(ρ)/(νa) for real ρ.
   - F(r*) = (DP + D log W_l)/(νa), so |F(r*)| ≤ (½/a + λ_a c1)/ν. Here λ_a bounds λ_W/a on the Cauchy region: the polar box [r_lo−δ−Δ_g, r_hi+δ+Δ_g] × [−θ, θ], with θ = φ1 + asin(δ/(r_lo−Δ_g)) < π/2.
   - The radius is δ = 0.3r_lo for r_lo < 2/3 and δ = 0.2 otherwise; Δ_g ≤ 0.2r_lo.
   - On [r_lo − Δ_g, r_hi + Δ_g], dF/dρ = D²f/(νaρ) ≥ (B(r_lo−Δ_g) − λ_a c2/ν0)/(r_hi + Δ_g), where B(ρ) = αρ + ρ²/s(ρ) is increasing and D²P ≥ 0 was used.
   - Hence F has a unique zero r_ex with |r_ex − r*| ≤ Δ, provided Δ < Δ_g (checked; Δ_g is enlarged once if needed).
   - For r* < 0.15 the Cauchy region is the disc |ζ| ≤ 0.3, and δ is the distance to its boundary.
2. **Contour P_l through r_ex.**
   - If r_ex ≤ y_J: the arc |ζ| = r_ex for |φ| ≤ π/2, then the segments [±ir_ex, ±iy_J].
   - If r_ex > y_J: the arc for |φ| ≤ φ_c = asin(y_J/r_ex), then the horizontal segments from x_c ± iy_J to ±iy_J, where x_c = √(r_ex² − y_J²).
   - Boxes straddling y_J add the bounds of both types. The check φ1 < asin(y_J/r_max) guarantees that the central arc |φ| ≤ φ1 lies on the contour.
3. **Value at r_ex.** By (F3), e^{f(r_ex)} ≥ e^{νΦ}(1+r_ex²)^{1/4}W_l(r_ex) ≥ e^{νΦ − λ_{W0}}, because r_ex ∈ [0.8r*, 1.2r*].
4. **Central arc (Gaussian lemma).**
   - Put g(φ) = f(r_ex e^{iφ}) − f(r_ex). Then g(0) = g′(0) = 0, g″ = −D²f, and g(−φ) is the complex conjugate of g(φ).
   - Let β = D²f(r_ex) > 0. Let M3 and M4 bound |D³f| and |D⁴f| on the arc |φ| ≤ φ1. Put β₋ = β − M4φ1²/12 > 0.
   - Taylor's theorem gives −βφ²/2 − M4φ⁴/24 ≤ Re g ≤ −β₋φ²/2 (the cubic term is purely imaginary) and |Im g| ≤ M3|φ|³/6.
   - Using Re e^g ≥ e^{Re g}(1 − (Im g)²/2) and Gaussian moments:

       (1/2π)∫_{|φ|≤φ1} Re e^g dφ ≥ (2πβ)^{−1/2}[1 − M4/(8β²) − (5/24)M3²β^{1/2}β₋^{−7/2} − 2e^{−βφ1²/2}/(φ1(2πβ)^{1/2})].

   - In the code β = νa·b with b ∈ [b_l, b_h] = [B(r_min) − λ_a c2/ν0, B(r_max) + (¼/a_lo + λ_a c2)/ν0], and r_min, r_max = r_lo − Δ_g, r_hi + Δ_g. Also M3 = νa·m3 and M4 = νa·m4, where m3 and m4 are sup bounds over the polar box [r_min, r_max] × [0, φ1] (angle chunks ≤ 0.02, F5, F6).
   - The quantities m3, m4, b_l, b_h are worst case over ν ≥ ν0, a ≥ a_lo. Each error term is then largest at νa = ν0a_lo.
   - Replacing β by β_hi = νa·b_h multiplies the relative error by at most √(b_h/b_l).
5. **Arc tail and segments**, relative to e^{f(r_ex)}.
   - *Arc tail* (φ1 ≤ |φ| ≤ φ_c). By (F1), |e^{g}| ≤ exp(νa[α_lo r_min(cos φ − 1) − r_min²J(φ)] + λ_W(ζ) + λ_{W0}), where J(φ) is a ball lower bound of ∫_0^φ sin2ψ/(2Re s)dψ, valid for every radius in [r_min, r_max]. Also Re P(ζ) ≤ P(r_ex), and the term −N log|ζ| is constant on the arc.
   - *Imaginary segments.* Use (F0) Ψ_e(iy) ≤ Ψ_e(ir_ex), (F1) at φ = π/2, and ∫_{r_ex}^{y_J}(r_ex/y)^N dy/y ≤ 1/N with N ≥ ν0 t_min.
   - *Horizontal segments.* The exponent per νa at ζ = x + iy_J is

       α(x − ρ) + ReΨ_e(ζ) − Ψ_e(ρ) − (t/a)log(|ζ|/ρ), with ρ = r_ex and |ζ| ≤ ρ.

     It is at most the value at t' = DA_a(r_max) ≥ t, because it is increasing in t. For t' fixed it is increasing in ρ up to r_max, since its ρ-derivative is (t' − DA_a(ρ))/(aρ) ≥ 0. So `good.py` evaluates it at ρ = r_max.
   - *Monotonicity.* For every piece the code checks (λ terms) < ν0|deficit|/2 and ν0a_lo|deficit| − a_lo(λ terms) > 1. So each piece times √(2πνa b_h) decreases in ν and in a, and its value at (ν0, a_lo) bounds all ν ≥ ν0 and all a in the box.
6. **Combination.** J_l is real (by conjugate symmetry), so J_l ≥ e^{f(r_ex)}(2πβ_hi)^{−1/2}(1 − e), with e = E√(b_h/b_l) + (tail + seg)√(2πν0a_lo b_h). Together with κ^l ≥ e^{−l/(12ν)}, this proves the lemma. ∎

**Lemma 3′ (l = 0, r* ≤ 0.15, n ≥ 8; `goodsmall.py`).** This is the same lemma, parametrised by n ≥ n0 := max(8, ν0 t(r_lo)), where t(r) = r + s(r) − 1.

- Since n ≤ ν t(r_hi), the effective ν is at least n0/t(r_hi), and the LG bounds on |ζ| ≤ 0.3 are taken there.
- *Critical point.* Work in u = log ρ.
  - F′(u) = D²f ≥ νρ − λc2 ≥ 0.99 n/(1 + r_hi/2) − λc2 for ρ ≥ 0.99r*, using D²A_1(ρ) ≥ ρ and n = νt(r*) ≤ νr*(1 + r*/2).
  - |F(u*)| ≤ r_max0²/2 + λc1, with r_max0 = 1.01r_hi.
  - Hence |log(r_ex/r*)| ≤ drel < 0.009, so r_ex/r* ∈ [1 − drel, e^{drel}] and r_ex lies inside the Cauchy radius 1.01r_hi.
- *Gaussian constants.*
  - b_l = n(1−drel)/(1+r_hi/2) − λc2
  - b_h = n e^{drel}(1+r_max) + r_max² + λc2
  - M3 = n e^{drel}(1 + r_max sup|(2+ζ²)/s³|) + sup|D³P| + λc3, and M4 analogously, using νr_ex ≤ n e^{drel}.
- *Tail.* The deficit is at most n κ_lo(cos φ − 1), with κ_lo = (1−drel)/(1+r_hi/2).
- *Segments.* The contribution is at most e^{−nκ_lo}/(πn) times the W-ratio e^{2λ}.
- *Monotonicity in n.* Each term is a quotient of affine functions of n. The signs needed are checked: the n-slope of β₋ is positive, and n0κ_lo(1 − cos φ1) > 1. So the value at n0 bounds all n ≥ n0.
- The table stores β/ν ≤ t(r_hi)e^{drel}(1+r_max) + (r_max² + λc2)/ν_eff.

**Tables** (`mktab.py` → `tab203.json`, ν0 = 203.5). Boxes of r* have relative width 3%, split at y_J.

| family | α range | r* covered | max e |
|---|---|---|---|
| l = 0, r* ≤ 0.15 (Lemma 3′) | α = 1 | [0, 0.15] | 0.089 |
| l = 0 | α = 1 | [0.15, 1.6] | 0.029 |
| l = 1 | α = 1/3 | [0.03, 1.6] | 0.320 |
| l = 2 | α = 1/5 | [0.05, 1.6] | 0.118 |
| l = 3 | α = 1/7 | [0.07, 1.6] | 0.071 |
| l = 4 | α = 1/9 | [0.08, 1.6] | 0.058 |
| l ≥ 5 | α ∈ [0, 1/11] | [0.116, 1.6] | 0.490 |

## 5. Bad terms (Lemma 4, `assembly.py`)

All bounds below are upper bounds for |κ^l J_l| or |I_L|. They have the form Σ len·exp(νR + c), where R and c are ball upper bounds valid for all ν ≥ ν0, and κ ≤ 1 is used.

- **Crude(l, contour radius r_c), for N ≥ 0 and t ∈ [t_lo, t_hi].**
  - Arc: R = A_a(r_c) − t log r_c (worst t endpoint), c = ¼log(1+r_c²) + λ_W, len = φ_c/π. This uses (F1) and |1+ζ²| ≤ 1+r_c².
  - If r_c ≤ y_J: segments with R = aΨ_e(ir_c) − t_hi log r_c, c = λ_W + log log(y_J/r_c), len = 1/π. This uses (F0), y^{−N} ≤ r_c^{−N} and ∫dy/y.
  - If r_c > y_J: 30 horizontal pieces, each bounded by ball suprema of ReA_a − t log|ζ| and ¼log|1+ζ²| − log|ζ| + λ_W.
  - r_c is the saddle at mid-box, clamped to [0.01, 1.6]. Any r_c is a valid contour.
- **Neg (all l ≥ LN with N_l ≤ −1, where LN = ⌊τ_a/2⌋ + 1).** Use the segment contour and sum the geometric series.

    Σ |κ^lJ_l| ≤ (1/π)∫_0^{y_J} e^{ν(2LN+1)Ψ_e(iy) + λ_K + LN(λ_I+λ_K)} y^{−N_{LN}−1} dy/(1−X_max).

  - It uses (F0), |X(iy)| ≤ e^{2νReη(iy)+λ_I+λ_K}, |e^z/p| ≤ e^{νΨ_e(iy)+λ_K} and −N_{LN} = ν(2LN−τ) ≥ 1.
  - It is integrated over 40 pieces [y_0, y_1] of [0, y_J], with |N| = −N_{LN}.
  - Case ν0(2LN−τ_b) ≥ 1: on the first piece, ∫_0^{y_1} y^{|N|−1}dy = y_1^{|N|}/|N| ≤ y_1^{ν(2LN−τ_b)}. On the other pieces, y^{|N|−1} ≤ y^{ν(2LN−τ_b)}/y_0, and the 1/y_0 goes into the prefactor.
  - Otherwise: y^{|N|−1} ≤ 1, since |N| ≥ 1.
  - τ-boxes never contain an even integer in their interior. A right end τ_b = 2LN is excluded from the box (it is the left end of the next box), so N_{LN} ≤ −1 holds on the box. l with N_l ≥ 0 are covered by the l-loop.
- **I_L.** On the left half-circle, the integrand modulus is y_J^{−n} e^{νReΨ}|1+ζ²|^{1/4}/|1+h_K(−ζ)|, with Ψ = ζ + Ψ_e. This gives R = ReΨ − τ_b log y_J and c = ¼log|1+ζ²| + λ_K(−ζ), over 24 angle pieces per half; together the halves have len (f1−f0)/π.

## 6. Proposition 1: n ≥ 8 and 0 ≤ τ < 60 (`assembly.py`)

`assembly.py` partitions [0, 60) into 3000 boxes [τ_a, τ_b) of width ≤ 0.02, cut at even integers. For each box it does the following:

1. For every l ≤ ⌊τ_b/2⌋ it computes the t-range [max(0, τ_a − 2l), τ_b − 2l] and the saddle range (widened by 1e−12).
   - The l is "good" if l = 0, or t_lo > 0; and r_g(l) ≤ r* ≤ 1.6; and the tables cover the r*-range with e < 0.9.
   - The thresholds are r_g = 0, 0.03, 0.05, 0.07, 0.08 for l ≤ 4, and r_g = 0.116 for l ≥ 5.
   - Otherwise l is "bad", and its Crude bound is added.
2. It takes the best good l*, with lower bound exp(ν0Φ* + G*) where Φ* = min(Φ(t_lo), Φ(t_hi)) by (F4).
3. It adds the Neg and I_L bounds.
4. It requires ratio := Σ_bad len·exp(ν0(R − Φ*) + c)/e^{G*} < ½ and ν0·min(Φ* − R) > 1.

**Validity for all ν ≥ ν0.** At ν the good lower bound has an extra factor (ν0/ν)^{1/2} relative to G*. Each bad term then changes by exp(−(ν−ν0)(Φ*−R))·(ν/ν0)^{1/2}, which is ≤ 1 because ν0(Φ*−R) > ½.

**Conclusion.** Every other good term is positive by Lemma 3. Hence

  (2ν)^n c_n/Λ ≥ κ^{l*}J_{l*}(1 − ratio) > 0.

The l = 0 rows with r* ≤ 0.15 need n ≥ 8; this is where n ≥ 8 enters.

Result: 3000/3000 boxes pass, with worst ratio 5.5e−18 (attained for τ < 12) and minimum rate gap Φ* − R = 0.183, so ν0·gap ≥ 37. For larger τ the ratios are much smaller: below 1e−124 for τ ≥ 12.

## 7. Proposition 2: τ ≥ 60 (`large_tau.py`)

**Choice of l̂.** Let l̂ be the largest l with r_l ≥ 0.5. Here r_l := r*(2l+1, τ − 2l) decreases in l.

- Since r_{25}(τ) ≥ r_{25}(60) > 0.5, l̂ ≥ 25.
- Let ρ(l) be the saddle as a function of a continuous l. Then dρ/dl = −2s²/(s + aρ), and |dρ/dl| ≤ 0.1055 < 0.135 on [0.5, 0.635] for a ≥ 51. Hence r_{l̂} ∈ [0.5, 0.635].

**Main term.** Since a = (τ + 1 − r)/s,

  Φ = r(1 − log r) + (1−r)g(r) + τg(r),  where g = (Ψ_e − (s−1)log r)/s.

So Φ̂ ≥ I0 + g_min τ with I0 = 1.0464 and g_min = 0.40154. On r ∈ [0.5, 0.635] the l ≥ 5 rows have e ≤ 0.0055, b_h ≤ 0.4153 and λ_{W0}/a ≤ 0.0014. The l ≥ 5 table rows give

  G ≥ G_0 − ½log(ν/ν0),  with G_0 := log(1−e) − a·λ − τ/(24ν0) − ½log(2πν0 a b_h),

where a ≤ τ + 1 and l ≤ τ/2. The code uses G_0. The factor (ν0/ν)^{1/2} is absorbed below by the ν-monotonicity step.

**Other l.** Every other l falls in one of these classes:

- *Far* (r_l > 1.6; checked to include every l ≤ 4): crude bound at r_c = 1.6. The rate is linear in a.
  - On the arc it is increasing in a, since η(1.6) > 0. The worst case is a = (τ − 0.6)/s(1.6).
  - On the horizontal pieces, both endpoints in a are bounded.
  - There are at most τ/(2s(1.6)) + 1 far l, and λ_W ≤ a(λ_I + λ_K) ≤ (τ+1)(λ_I+λ_K).
- *Start* (r_l < 0.116, N ≥ 0): crude bound at r_c = 0.116, using t < DA_a(0.116) and a ≤ τ + 1. The a-coefficients are positive (checked). There are at most τ/2 + 1 such l.
- *Good, r_l ∈ [0.116, 1.6], l ≥ 5*: J_l > 0, because the table rows cover this range with e < 0.9 (checked).
- *Neg* (N < 0): rate ≤ (τ+3)Ψ_e(0).
- *I_L*: rate ≤ max ReΨ − τ log y_J.

**Final inequality.**

- Every class has the form rate ≤ C0 + C1τ with gap Φ̂ − rate ≥ D0 + D1(τ − 60), where D1 > 0 and ν0D0 > 1.
- The log of (count × bad term / main term) has τ-derivative at most −ν0D1 + λ_c + 1/60 + λ + 1/(24ν0) + 1/120. This is < 0 for every class (checked).
- The exponent, including the −½log(ν/ν0) from G, has ν-derivative −gap + 1/(2ν) < 0, because ν0·gap > 1.
- So the worst case is τ = 60 and ν = ν0, and there the total is 8.9e−435 < ½.

At τ = 60 the class gaps are 11.4 (far arc), 13.4 (far horizontal), 4.97 and 5.50 (start), 5.81 (Neg) and 15.3 (I_L). The τ-derivatives of the log ratios are ≤ −15.

## 8. Proposition 3: n ≤ 12 (`smalln.py`)

c_n(m) is a rational function of m (from the recursion c_n = −Σ_j q_j c_{n−j}). For n ≤ 12, after substituting m = 203 + x, its numerator and denominator have coefficients all of one sign. Hence c_n > 0 for all real m ≥ 203.

Examples: c_1 = 1/2, c_2 = m/(4(2m−1)), c_3 = (m+1)/(24(2m−1)), c_4 = m(2m²+m−7)/(96(2m−3)(2m−1)²).

**Coverage.** Propositions 1–3 cover every n for every odd m ≥ 203:

- n ≤ 12 by Proposition 3;
- n ≥ 8 with n/ν < 60 by Proposition 1;
- n/ν ≥ 60 by Proposition 2.

## 9. Proposition 4: odd 3 ≤ m ≤ 201 (`../FasiHigham/arb_cert.py diag 201`)

For each m, do the following:

1. Compute certified isolating balls for all roots of b_0q.
2. Check that the smallest-modulus root r is strictly smaller in modulus than all the other roots. So r is real (a conjugate would have the same modulus); it is positive by (I1).
3. Check R_1 = −b_0/(rQ′(r)) > 0.
4. Partial fractions give c_n r^n = R_1 + Σ_{j≥2}R_j(r/z_j)^n, so c_n > 0 for n ≥ n_1. Here n_1 is the least n with Σ|R_j/R_1|(r/|z_j|)^n < 1; this sum is decreasing in n.
5. For n < n_1, compute the Arb power series of 1/q and certify every coefficient > 0.

All 100 values of m return POS (step G). Independently, Goodchild (Zenodo 23919117) certifies odd m ≤ 361 with the same method.

## 10. Rigor bookkeeping

- All analytic quantities are Arb balls.
- Floats appear only in the following places, and only with outward rounding:
  - table storage (values rounded up by ×(1+1e−9));
  - grid cell indices (padding ±1e−12);
  - saddle ranges (widening ±1e−12).
- The τ endpoints are floats converted exactly.
- Every "valid for all ν ≥ ν0" claim rests on one of two facts:
  - ε and λ decrease in ν;
  - an explicitly checked monotonicity inequality (stated where used).
- Cross-checks against exact values at m = 203 (`verify.py --sanity`):
  - S1: the assembled lower bound is ≤ the exact c_n for n ≤ 4000 (every third n), with min log(true/LB) = 0.044.
  - S2: the Good bound is ≤ the quadrature value of κ^lJ_l in 12 cases, with log-margins 0.056–0.31.
  - S3 and S4: the I_L, Neg and Crude bounds are ≥ the numeric |integrals|.

## 11. Verification

`python3 verify.py [--sanity]` (python-flint ≥ 0.9, mpmath, sympy) runs the following steps one after another, single-threaded:

- A: identities
- B: n ≤ 12
- C: tables
- D: LG grid
- E: five sweeps of τ in [0, 60]
- F: τ ≥ 60
- G: m ≤ 201

With --sanity it also runs the cross-checks S1–S4. Each step runs under a wall-clock limit of at most 1800 s. Run times and outputs are in REPORT.md.

## 12. Novelty

- Scardino, arXiv:2610.39051 (2026), proves k = 0 and eventual positivity under the "real zero is the unique minimal-modulus zero" condition. It states that "the diagonal case k = m … remains open".
- Goodchild (Zenodo 23919117, 2026-11-03) certifies odd m ≤ 361 one m at a time.
- No proof for all m was found (arXiv listings for Nov 2026, Zenodo, and the GitHub history checked 2026-11-05).
- New ingredients of this proof:
  - the exact expansion 1/q = (e^z/p)Σ(πI_ν/K_ν)^l, with each term bounded through its own saddle;
  - uniform Liouville–Green bounds, which need no localisation of zeros.
