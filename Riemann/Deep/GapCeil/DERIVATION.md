# DERIVATION: Bui–Milinovich's large-gap criterion as an operator pencil, its ceiling, and the Arb certificate

Sources: BM = Bui–Milinovich, arXiv 1410.3635v3 (Thm smoothed l.300-313, Lemma main l.370-384, reduction l.164-202, l.327-330);
GapLambda = `Riemann/Deep/GapLambda/` (DERIVATION_BM.md: the reduction re-derived, the 9-fold c_j and its 4-term form).
Λ > κ follows if h1 = c0/(κ²c1) > 1 for some admissible (r, P, v, ϑ < 1/4) (BM l.327-330).

## 1. The criterion in X-variables (r = 1)
GapLambda §3 writes, after integrating t3, t4 exactly,
c_j = (κπ)⁻² ∫ ab (v−S)^{2j} P(w1)P(w2)(A+B+C+D), x, x1..x4 ∈ polytope, t1, t2 ∈ [0,1],
a = 1−ϑ(x1+x3), b = 1−ϑ(x2+x4), s = at1−bt2, S = ϑ(x+Σx_i)+at1+bt2, w1 = 1−x−x1−x2, w2 = 1−x−x3−x4.
1. **C and D simplify.** With c = κπ: C = −e^{icϑ(x2−x3)}e^{−ics}, D = −e^{icϑ(x1−x4)}e^{ics}, and
   A+B+C+D = (e^{icX1} − e^{icX2}e^{−icσ})(e^{−icX3} − e^{−icX4}e^{icσ}) with X_i = ϑx_i, σ = s.
2. **Substitute** X = ϑx, X_i = ϑx_i, τ1 = at1, τ2 = bt2. ab dt1dt2 = dτ1dτ2, s = τ1−τ2, S = X+ΣX_i+τ1+τ2, and
   c_j = (κπ)⁻²ϑ⁻⁵ Q_j[f], Q_j[f] = ∫ f(X+X1+X2) f(X+X3+X4) ρ_j(X, X1..X4) dX dX1..dX4,
   ρ_j = ∫_0^{1−X1−X3}∫_0^{1−X2−X4} (v−S)^{2j} [e^{ic(X1−X3)} + e^{ic(X2−X4)} − e^{ic(X1−X4)}e^{ic(τ1−τ2)} − e^{ic(X2−X3)}e^{−ic(τ1−τ2)}] dτ1dτ2,
   f(U) = P(1−U/ϑ) on [0, ϑ], f = 0 beyond.
3. **ϑ enters only through supp f.** The prefactor is the same for j = 0, 1, so h1 = Q_0/(κ²Q_1). For ϑ' < ϑ the amplifier
   P(1−U/ϑ')·1[U ≤ ϑ'] is an L² function on [0, ϑ]. Hence the ceiling is non-decreasing in ϑ, and ϑ = 1/4 with
   f ∈ L²[0, 1/4] covers every admissible ϑ < 1/4. (The τ-box is never empty: X1+X3 ≤ 1/2.)
4. **General r** (BM's d_r(h)P[h]): extra weight X^{r²−1}(X1X2X3X4)^{r−1}; normalising constants cancel in h1.

## 2. φ-representation: multiplication + compact
Parametrise X1 = r1θ1, X2 = r1(1−θ1), X3 = r2θ2, X4 = r2(1−θ2), so U1 = X+r1 and U2 = X+r2:
Q_j = ∫_0^{1/4} dX ∫∫ k_j(X,r1,r2) f(X+r1)f(X+r2) with k_j = r1r2·ρ̄_j, ρ̄_j the θ-average (fiber length r1 r2).
ρ_j depends on X only through V = v−X−r1−r2: ρ̄_1 = (v−X)²a − 2(v−X)b + c with a = ρ̄_0 and b, c built from three
θ-averaged functions g0, g1, g2 of (r1, r2) (`gcore.py` `rho_parts`).
- Expand ρ̄ in a Taylor series in (r1, r2) (entire). With F_m(X) = ∫_X^{1/4}(u−X)^m f(u)du, φ := F_1 = I²f, I g(X) = ∫_X^{1/4} g, and F_m = m! I^{m−1}φ:
  **Q_j = Σ_{i,j≥0} Γ_ij ⟨I^iφ, X^p I^jφ⟩**, Γ_ij = (i+1)!(j+1)!·[r1^i r2^j]ρ̄. So Q_0 = ⟨φ, M0φ⟩ and Q_1(v) = ⟨φ,(v²M0 − 2vMS + MSS)φ⟩,
  with M0 = G_0^(a), MS = G_1^(a)+G_0^(b), MSS = G_2^(a)+2G_1^(b)+G_0^(c), G_p^(h) = Σ Γ^(h)_ij (I^i)^T X^p I^j.
- **Structure.** The (0,0) term is multiplication by ρ_j(X,0) = ∫∫(v−X−τ1−τ2)^{2j}·4sin²(c(τ1−τ2)/2)dτ; every other term
  contains a Volterra factor I, so M_j = multiplication + compact. The pencil is therefore well conditioned in φ, while in P
  both forms are compact (hence the ill-conditioned high-degree scans of GapLambda).
- **Closure.** {I²P} (polynomials with a double zero at 1/4) is dense in L²[0,1/4] and M_j are bounded, so
  sup over all polynomials P = sup over φ ∈ L². An L² optimiser φ* with φ*(1/4), φ*'(1/4) ≠ 0 corresponds to
  f = φ*'' + c1δ_{1/4} + c2δ'_{1/4}: point masses at h = y. Polynomials converge to it only algebraically; this is
  why κ_max(deg P) creeps up by ≈ 0.001 per degree.
- For general r: φ = I^{2r}f, Γ_ij = (i+2r−1)!(j+2r−1)![r1^i r2^j]ρ̄^(r) (θ-weights (θ(1−θ))^{r−1}), weights X^{r²−1+p} (`gcore_r.py`).

## 3. The pencil and its ceiling
κ*_r = sup{κ : sup_v λ_max(M0, M1(v); κ) > κ²}, λ_max = top generalised eigenvalue (Galerkin in shifted Legendre
polynomials of φ on [0, 1/4]; Taylor coefficients by a 2-D FFT on |r| = 0.4, D = 36 terms).
- **Checks.** BM's own P (κ = 3.18, v = 1.26, ϑ = 1/4) gives h1 = 1.002177814 (Arb certificate of GapLambda: 1.00217781433).
  Degree-2/8/14 polynomial subspaces give κ_max = 3.183185/3.255011/3.265064 (GapLambda: 3.183/3.2548/3.2648 at ϑ = 0.2499).
  The r = 1 and r = 2 polynomial values reproduce GapLambda's scans at ϑ = 0.2499 (`logs/float/check_r2.txt`).
- The Arb-exact matrices of §4 give the same pencil values (κ = 3.28: 0.9931487 both; κ = 5: 0.3709250; κ = 6: 0.2334784).

## 4. Certificate (r = 1): κ²Q_1(v) − Q_0 ≥ 0 for every φ ∈ L²[0,1/4] and every v ∈ ℝ
Implemented in `gcert.py` (operators) and `runcert.py` (covering of v); Arb balls throughout, κ exact rational.
1. **Γ arrays exactly.** ρ_j is a sum of products u(X1,X3)·w(X2,X4), u = e^{i(α1X1+α3X3)} m_k(1−X1−X3, z), m_k(A,z) = ∫_0^A τ^k e^{zτ}dτ,
   z ∈ {0, ±ic}. The θ-average turns X1^aX3^bX2^cX4^d into r1^{a+c}r2^{b+d}·a!c!b!d!/((a+c+1)!(b+d+1)!), so Γ = û∗ŵ (2-D
   convolution of derivative arrays û_ab = ∂^a_{X1}∂^b_{X3}u(0)), computed by Kronecker substitution in acb_poly.
   Bound used for the tail: |û_ab| ≤ γ^{a+b}, γ ≥ 2c+2, hence |Γ^(h)_ij| ≤ 4(i+j+4)²(i+1)(j+1)γ^{i+j}.
2. **Legendre algebra exactly.** In P̃_n = P_n(8X−1): I is tridiagonal (∫_X^{1/4}P̃_n = (P̃_{n−1}−P̃_{n+1})/(8(2n+1))), X is tridiagonal,
   Gram = diag(1/(4(2n+1))). A_kl = ⟨P̃_k, BP̃_l⟩ for k < N_w = N+2D+4, l < N (Horner in I^T on a window of size N_w+D+6).
   With 0 ≤ i, j ≤ D−1 and p ≤ 2, (I^T)^i X^p I^j raises the degree by at most 2D, so B maps degree < N into degree
   ≤ N+2D−1 < N_w: the window holds P_W B p exactly for the truncated series, and no Horner step is cut. Series tail (i or j ≥ D): operator norm ≤ ε = Σ|coef|Σ|Γ-bound|ϑ^{i+j+p}/(i!j!).
3. **High modes W = span{P̃_n : n ≥ N}.** For w ∈ W, (Iw)(0) = ⟨w,1⟩ = 0 and (Iw)(1/4) = 0, so integration by parts gives
   (0,1)+(1,0): pΓ_01⟨Iw, X^{p−1}Iw⟩; (0,j)+(j,0), j ≥ 2: 2Γ_0j⟨Iw, pX^{p−1}I^jw − X^pI^{j−1}w⟩; (i,j ≥ 1): Γ_ij⟨I^iw, X^pI^jw⟩.
   With |I^kw| ≤ t_k|w|, t_k = Π_{l<k} 2β_{N−l}, β_m = (1/8)/√((2m−1)(2m+1)) (I maps span{n ≥ m} to span{n ≥ m−1} with norm ≤ 2β_m),
   ⟨w,Bw⟩ ≥ (δ_D − C)|w|² =: δ_W|w|², δ_D = min_{[0,1/4]} d(X) of the multiplication part (256 subintervals).
4. **Schur step.** For φ = p + w: ⟨φ,Bφ⟩ ≥ ⟨p,Bp⟩ − |P_WBp|²/δ_W, and |P_WBp|² ≤ (1+η)pᵀEᵀΛEp + (1+1/η)ε²|p|²,
   E = A[N:N_w, :N], Λ = diag((2k+1)·4) (exact by step 2), η = 10⁻⁶. Test S = A_VV − (…)/δ_W ≻ 0 by Arb Cholesky.
   Run parameters: N = 24, D = max(50, ⌊2.5γ/4⌋+34), γ = ⌊2κπ+3⌋+1. Every operator's δ_D, C, δ_W, ε and Cholesky result
   is in `logs/cert_<κ>.json`.
5. **Covering v ∈ ℝ** (all operators scalar): (i) M0 ⪰ 0; (ii) MS − v_lo M0 ⪰ 0 and v_hi M0 − MS ⪰ 0 (v_lo = 0.8, v_hi = 1.6),
   so κ²c1(v) − c0 ≥ κ²c1(v_lo) − c0 for v ≤ v_lo (and likewise above v_hi); (iii) cells [v_m−h, v_m+h] tiling [v_lo, v_hi]:
   κ²c1(v) − c0 ≥ min_± (κ²c1(v_m) − c0 ± 2κ²h⟨φ,(MS − v_mM0)φ⟩) because κ²τ²c0 ≥ 0; both signs certified.
   A block relaxation treating v as an operator ([[κ²M0, −κ²MS],[−κ²MS, κ²MSS − M0]] ⪰ 0) needs no cells, but it holds only
   from κ = 3.3015 (float), so it was dropped.
6. **Precision.** Radii grow like |I|^i|Γ_i| in the Horner sum (true values decay factorially, absolute values do not),
   so prec = 128 + D·log2(γ/4) + 32 bits: 289-297 bits for κ ≤ 3.8 and 308-411 bits for κ = 4-8.

## 5. One amplifier, general r: the record certificate (`certr.py`)
- For a polynomial P, φ = I^{2r}f is a polynomial on [0, ϑ] (ϑ = 2499/10000 used as the domain), so Q_j = φᵀAφ exactly
  in the Legendre window; only the series tail (i or j ≥ D) needs a bound: |Γ^(r)_ij| ≤ 4(i+j+8)²((i+2r)(j+2r))^{2r}γ^{i+j}
  (from |û^(r)_ab| ≤ (a+r−1)!(b+r−1)!/(a!b!)·γ^{a+b} and the shifts (i+2r−1)G_{i−1,j} + (j+2r−1)G_{i,j−1}).
  D = 80 makes it 5·10⁻³⁸.
- r-weights: û^(r)_ab = (a+r−1)!(b+r−1)![X1^aX3^b]u, because the θ-average of θ^{a+r−1}(1−θ)^{c+r−1} is (a+r−1)!(c+r−1)!/(a+c+2r−1)!.
  The operator carries X^{r²−1+p}.
- **Validation of r = 2:** the best degree-0..6 amplifiers at κ = 3.25, ϑ = 0.2499 reproduce GapLambda's direct 7-D
  quadrature of BM's printed r = 2 integrand to all 6 printed digits (`logs/float/check_r2.txt`). Arb and float agree to 7-8 digits
  on the record amplifiers.
- **Results:**

| r | deg P | ϑ | κ | v | h1 (Arb) |
|---|---|---|---|---|---|
| 2 | 24 | 2499/10000 | 3.29 | 1.349501 | 1.00041599450 ± 1.5·10⁻¹² |
| 2 | 36 | 2499/10000 | 3.296 | 1.35222 | 1.00029463305 ± 1·10⁻¹² |
| 2 | 60 | 2499/10000 | 3.30 | 1.353909 | 1.00000835431 ± 2.8·10⁻¹² |

  The exact P (shifted-Legendre coefficients as exact decimals) is in `logs/record_P_r2_d60.txt` (and `_d36`),
  `logs/certr_*.txt` and `logs/r2opt_*.json`.

## 6. What is not covered
- **κ between grid points** (r = 1 and r = 2). Each run certifies one exact κ. The operators depend analytically on κ, but every
  Lipschitz bound tried (ball κ, Cauchy estimates) loses a factor 10⁵-10¹⁰ through the cancelling Taylor series, so the
  certified set is the grid, not an interval. The float pencil is smooth and decreasing there (REPORT.md table).
- **Large κ.** Certified up to κ = 8 (r = 1) and 5 (r = 2), at grid points; above that the evidence is the float pencil
  (r = 1: sup_v λ/κ² = 0.0763 at κ = 10 from the Arb-exact matrices, `logs/evid/`). There is no analytic tail bound.
- **r ≥ 3 barrier.** Float only (3.226093, 3.141892, 3.071102 for r = 3, 4, 5); §7 applies verbatim but was run for r = 2 only,
  and no uniform-in-r bound is known.
- (Correction of an earlier draft: for r ≥ 2 the Volterra parts are not "Hardy operators, bounded but not compact". Their
  kernels are homogeneous of degree 0, which gives operators of homogeneity −1 like I itself, so they are compact. The only
  obstacle was that the unweighted high-mode bound used min d = 0; §7 removes it.)
- **Inputs.** BM's main terms (Thm smoothed, Lemma main) and BBLR's twisted fourth moment, as printed.

## 7. Barrier certificate for r ≥ 2: the weighted Jacobi basis (`gcert_r.py`, `runcert_r.py`)
For r ≥ 2 every term of B carries X^{b+p}, b = r²−1, so the natural space is L²(X^b dX) on [0,1/4]. The multiplication part
is X^b d(X) with d ≥ δ > 0.
1. **Basis.** q_n ∝ R_n = P_n^{(0,b)}(2X/ϑ−1), orthonormal in L²(X^b): ‖R_n‖² = ϑ^{b+1}/(2n+b+1). Built exactly from the
   Jacobi recurrence (DLMF 18.9.2, α = 0) in the Legendre coefficient space. `verify_q.py` checks Q W_b Qᵀ ∋ I in Arb.
2. **I shifts W_m = span{q_n : n ≥ m} into W_{m−1}.** For deg P < m−1, ⟨Iw, P⟩_b = ∫w(Y)[∫_0^Y X^bP] dY = ⟨w, Y·Q(Y)⟩_b
   with ∫_0^Y X^bP = Y^{b+1}Q(Y), deg YQ < m.
3. **I is tridiagonal.** On [0,1], I R_n = A_nR_{n−1} + B_nR_n + C_nR_{n+1}. Matching the t^{n+1} and t^n coefficients
   (R_n = K_n t^n − K_n s_n t^{n−1} + …, s_n = n(n+b)/(2n+b), K_n/K_{n+1} = (n+1)(n+b+1)/((2n+b+1)(2n+b+2))) and using
   R_n(1) = 1, (IR_n)(1) = 0:
   C_n = −(n+b+1)/((2n+b+1)(2n+b+2)), B_n = b/((2n+b)(2n+b+2)), A_n = −B_n − C_n > 0 (b = 0 gives Legendre).
   Orthonormal entries: |a_n| ≤ 1/(2n+b−1), |b_n| ≤ b/(2n+b)², |c_n| ≤ 1/(2n+b+1), all decreasing. So ‖Iw‖_b ≤ τ_m‖w‖_b on W_m
   with τ_m = ϑ(1/(2m+b−1) + b/(2m+b)² + 1/(2m+b+1)), and t_k = Π_{l<k} τ_{N−l}; ‖I^k‖_b ≤ ϑ^k/k!
   (Schur test for the kernel (X/Y)^{b/2}(Y−X)^{k−1}/(k−1)!).
4. **High modes.** For w ∈ W_N: ⟨w,Bw⟩ ≥ (δ − C)‖w‖_b² with
   - (0,1)+(1,0), p = 0: 2|Γ_01|t_1;
   - (0,1)+(1,0), p ≥ 1, by parts: (b+p)ϑ^{p−1}|Γ_01|t_1² (X^{b+p}(Iw)² → 0 at 0);
   - (0,j), j ≥ 2: 2|Γ_0j|ϑ^p t_j;
   - (i,j ≥ 1): |Γ_ij|ϑ^p t_i t_j;
   - tail: the r-bound of §5.

   The p = 0 term is linear in t_1 ≈ ϑ/N, hence N = 72.
5. **Schur step.** Split B = B_T + R, with B_T the truncated series (i, j < D) and R the tail.
   - B_T p = X^b g_p with g_p a polynomial of degree < N_w (I^T(X^b h) is divisible by X^{b+1}). So ⟨w, B_T p⟩ = ⟨w, g_p⟩_b,
     and ‖P_W g_p‖_b² = Σ_{k≥N}⟨q_k, B_T p⟩² is exact on the window k < N_w.
   - The tail obeys |⟨φ, Rφ'⟩| ≤ ε‖φ‖_b‖φ'‖_b, i.e. ‖X^{−b}R‖_{L²_b→L²_b} ≤ ε. It enters every block: entrywise radius ε on
     A_VV, the term (1+1/η)ε² in the V–W coupling, and the tail sum inside C (W–W).
   - Test A_VV − ((1+η)EᵀE + (1+1/η)ε²I)/(δ − C) ≻ 0 by Arb Cholesky. gcert.py (r = 1) does the same.
6. **v-covering** as in §4 (v_lo = 0.8, v_hi = 1.7). **Results:**

| κ | float sup_v λ/κ² | certified (r = 2, N = 72, cells shrink 0.7; 0.4 at 3.305) |
|---|---|---|
| 3.305 | 0.9981910 | yes: 39 operators (18 v-cells), D = 50, 289 bits, max C/δ = 0.37, 179 MB, 206 s |
| 3.31 | 0.9941218 | yes: 21 operators (9 v-cells), D = 50, 289 bits, max C/δ = 0.37, 166 MB, 168 s |
| 3.35 | 0.9625653 | yes: 15 operators (6 v-cells), D = 50, 292 bits, max C/δ = 0.36, 159 MB, 145 s |
| 3.50 | 0.8585740 | yes: 13 operators (5 v-cells), D = 50, 292 bits, max C/δ = 0.34, 187 MB, 151 s |
| 4.00 | 0.6284830 | yes: 9 operators (3 v-cells), D = 76, 377 bits, N = 48, max C/δ = 0.63, 214 MB, 184 s |
| 5.00 | 0.3956651 | yes: 9 operators (3 v-cells), D = 76, 397 bits, N = 48, max C/δ = 0.61, 284 MB, 147 s |

   - κ = 4 and 5 ran with N = 48 and D = 76. Their first runs, with D = 52 and 55, had failed even though the high-mode
     bound itself passed (C/δ ≤ 0.36).
   - Cause: the r = 2 tail bound carries ((i+4)(j+4))⁴, so the r = 1 rule D = 52-55 left ε ≈ 0.01-0.07, and
     (1+1/η)ε²/δ_W broke the Cholesky at index 0. The window entries E were ≈ 10⁻⁹⁶ (`diag_r2.py`). κ ≥ 4 needs D ≈ 76.
