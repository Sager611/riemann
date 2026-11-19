# FULLPROOF: the chain behind Λ > 3.324 (r = 2, midpoint mollifier, real P of degree 80)

Line numbers: BM = arXiv 1410.3635v3 TeX; BBLR = arXiv 1609.02539v1 TeX; GapCeil/GapBar = our DERIVATION.md files.
Λ := limsup_{n→∞} (t_{n+1} − t_n)/(2π/log t_n), t_n the ordinates of the zeros of ζ on the critical line, with multiplicity (BM l.53, l.63-66).

## Step 1. The criterion (BM §2.2, l.164-202), for F_mid
1. Suppose Λ ≤ κ (l.166). Put L = log(T/2π) (l.59-61), y = T^ϑ, ϑ = 2499/10000, r = 2, M(s) = Σ_{h≤y} d_r(h)P[h]h^{−s},
   P[h] = P(log(y/h)/log y) (l.243-250), P real, and
   F(t) = e^{ivtL} ζ(½+it) ζ(½+it+iκπ/L) M(½+it+iκπ/(2L)), v ∈ ℝ (BM's F, l.171, with M shifted).
2. F vanishes at every t_n and every t_n − κπ/L. Under Λ ≤ κ, consecutive t_n in [T, 2T] are ≤ (1+o(1))2κπ/L apart, so this union
   has gaps ≤ (1+o(1))κπ/L (l.182); other zeros of F only shorten gaps.
3. Wirtinger (i) for complex C¹ functions vanishing at both ends (l.139-155): ∫|F|² ≤ ((b−a)/π)²∫|F'|² ≤ (1+o(1))(κ/L)²∫|F'|² per gap
   (l.184-189); summing (l.191-192).
4. End pieces (l.194-196): t̃_1 − T, 2T − t̃_N ≪ 1; |F^{(k)}|² ≪ T^{2/3+ϑ+ε} ≪ T^{1−ε} since |M(½+it+iκπ/(2L))| ≤ Σ|d_r(h)P[h]|h^{−½}
   ≪ y^{½+ε}; the integrals are ≫ T because c_0 > 0 (Step 4).
5. Hence if h1 := limsup (L²/κ²)∫_T^{2T}|F|²/∫_T^{2T}|F'|² > 1 then Λ > κ (l.198-202).

## Step 2. The mean values (BM Lemma 2.4 "main", l.364-384; Thm 2.3 l.300-313; l.321-330)
1. I_{α,β}(M) = ∫ ζ(½+α1+it)ζ(½+α2+it)ζ(½+β1−it)ζ(½+β2−it) M(½+α3+it) M(½+β3−it) w(t/T) dt, "shifts α_j, β_j ≪ L⁻¹" (l.364-368).
2. Lemma main (ϑ < 1/4): I = c(α,β)·A_{r+2}(log y)^{r²+4r}L⁴ŵ(0)T/(2(r²−1)!((r−1)!)⁴) + O(TL^{(r+2)²−1}) uniformly (l.372-384), with
   c(α,β) = ∫_{[0,1]⁹, x+x1+x2≤1, x+x3+x4≤1} y^{−(α3+β3)x−α3(x1+x2)−β3(x3+x4)−β1x1−β2x2−α1x3−α2x4−(α2−α1)(x3−x4)t3−(β2−β1)(x1−x2)t4}
   (Ty^{−x1−x3})^{…}(Ty^{−x2−x4})^{…}(1−ϑ(x1+x3))(1−ϑ(x2+x4))(ϑ(x1−x2)+D)(ϑ(x3−x4)+D) x^{r²−1}(x1x2x3x4)^{r−1}
   P(1−x−x1−x2)P(1−x−x3−x4), D = (1−ϑ(x1+x3))t1 − (1−ϑ(x2+x4))t2 (l.377-383).
3. Its proof uses BBLR Thm 1.1 (BBLR l.103-136: two independent Dirichlet polynomials of length T^ϑ, coefficients ≪ a^ε, complex
   shifts ≪ 1/log T, error T^{1/2+2ϑ+ε} + T^{3/4+ϑ+ε}); BM quote it at l.500-517. α3, β3 enter only as h^{−α3}, k^{−β3} (l.531, 538),
   become y^{−(α3+β3)x−α3(x1+x2)−β3(x3+x4)} via Lemmas 4.1-4.2 (l.440-492, l.580-605), and survive the swaps giving I_2..I_6 (l.610-666).
4. |F_mid|² ↔ α1 = β1 = 0, α2 = −β2 = iκπ/L (BM's own point, l.391), α3 = −β3 = iκπ/(2L) (real P: conj M(½+α3+it) = M(½+ᾱ3−it)).
   All shifts ≪ 1/L. The extra factor is exp(−i(κπ/2)(log y/L)(x1+x2−x3−x4)) = e^{−iκπϑ(x1+x2−x3−x4)/2} + O(1/L).
5. F_mid' (l.393-397): F'/(ie^{ivtL}) = L·Q(L⁻¹(∂_{α1}+∂_{α2}+∂_{α3}))[ζζM] at α = 0, Q(x) = v + x; Cauchy circles of radius ≍ 1/L
   (l.417); on c(α,β) both operator sums give the factor −S, S = ϑ(x+Σx_i) + (1−ϑ(x1+x3))t1 + (1−ϑ(x2+x4))t2 (l.420-433).
   So ∫|F_mid^{(j)}|²w = c_j^mid·(same constant)·L^{4+2j}T + O(TL^{(r+2)²+2j−1}), c_j^mid = the l.377 integral at the point of 2.4
   times (v−S)^{2j}. Unsmoothing (l.321-325, standard; |F|² ≥ 0) gives the same on [T, 2T], and h1 = c_0^mid/(κ²c_1^mid) (l.327-330).

## Step 3. Reduction to the twisted kernel (own derivation; agrees with GapCeil §1-2 and GapBar §3)
1. t3, t4 integrate exactly: ∫∫e^{−iat3+ibt4}·(a/κπ)(b/κπ) dt3dt4 = (1−e^{−ia})(1−e^{ib})/(κπ)², a = κπ(ϑ(x3−x4)+D), b = κπ(ϑ(x1−x2)+D).
2. X_i = ϑx_i, X = ϑx, τ1 = (1−X1−X3)t1, τ2 = (1−X2−X4)t2 (Jacobian cancels BM's (1−ϑ(x1+x3))(1−ϑ(x2+x4))), σ = τ1 − τ2, c = κπ:
   e^{ic(X2−X4)}(1−e^{−ia})(1−e^{ib}) = e^{ic(X1−X3)} + e^{ic(X2−X4)} − e^{ic(X1−X4)}e^{icσ} − e^{ic(X2−X3)}e^{−icσ} =: U_L conj U_R.
3. Midpoint: e^{ic(X3+X4−X1−X2)/2}·U_L conj U_R = 2cos(c(X1−X2−X3+X4)/2) − 2cos(c(X1−X2+X3−X4)/2 + cσ)
   = 4 sin(c(X1−X2+σ)/2)·sin(c(X3−X4+σ)/2) (a = c(X1−X2)/2, b = c(X3−X4)/2, θ = cσ/2: 2cos(a−b) − 2cos(a+b+2θ) = 4 sin(a+θ)sin(b+θ)).
   Checked at 2000 random points to 1.9·10⁻⁵⁰ (`logs/idmaj.txt`).
4. c_j^mid ∝ Q_j[f] = ∫ f(X+X1+X2) f(X+X3+X4) X^{r²−1}(X1X2X3X4)^{r−1} ρ̃_j, ρ̃_j = ∫∫_{τ-box}(v−S)^{2j}·(3.) dτ, f(U) = P(1−U/ϑ) on [0, ϑ];
   the constant is the same for j = 0, 1.
5. Validation: `indep7.py` evaluates Steps 2.2-3.1 by Gauss quadrature from BM's TeX only: BM's own data give 1.0021778143 (BM: > 1.0002);
   the degree-8 midpoint P gives 1.0156026615 = certmid's Arb 1.01560266154 ± 4.5·10⁻¹³.

## Step 4. The certificate (GapCeil §2, §4.1-4.2, §5; GapBar certmid.py → certr.py + gtw.py)
1. Taylor expansion of the fiber-averaged ρ̃_j in (r1, r2) = (X1+X2, X3+X4): Q_0 = Σ_{i,j} Γ^a_ij⟨I^iφ, X^{q0}I^jφ⟩,
   Q_1 = v²Q_0 − 2v(G_1^a + G_0^b) + (G_2^a + 2G_1^b + G_0^c), φ = I^{2r}f, I g(X) = ∫_X^ϑ g, q0 = r² − 1 = 3.
2. Γ arrays exact in Arb (prec 366) for i, j < D = 80; φ and every I^iφ exact in the shifted-Legendre window on [0, ϑ].
3. Tail i or j ≥ D: |Γ_ij| ≤ 4(i+j+8)²((i+4)(j+4))⁴·24^{i+j} (24 ≥ 2c+2), ‖I^i‖ ≤ ϑ^i/i!: adds radius ≤ 5.6·10⁻⁴⁷ to H.
   Computed arrays sit ≥ 10⁷ below the majorant; the certificate tolerates a 10¹⁰-fold error in it.
4. Result (κ = 831/250, v = 683629/500000, ϑ = 2499/10000, deg P = 80, coefficients in `GapBar/logs/certmid_r2_d80_k3.324.txt`):
   Q1 > 0 and H = Q0 − κ²Q1 = 6.5058025·10⁻³⁷ ± 5.64·10⁻⁴⁷ > 0, so h1 = 1.00056020750 ± 4.5·10⁻¹⁵ > 1. Rerun identical
   (`logs/rerun_certmid_r2_d80_k3.324.txt`).

## Conclusion
Steps 1-4: Λ > 3.324, unconditionally. Inputs: BM Lemma main and BBLR Thm 1.1 as printed (no extension needed for real P).
