RECORD: θ > 24699/50000 = 0.49398 (ν < ν₁(θ) = (696θ − 229)/934 = (3θ − 1)/4 + 9(1 − θ)/1868; exponent pair (13/502, 445/502); (1/30, 13/15) gives θ > 0.49402)

# B3Nu — the dual range of Theorem E1′ beyond ν < (3θ − 1)/4 (node N8)
b3nu (Orca), 20 Nov 2026, 06:26– UTC. Code code/, logs logs/ (every run single-threaded under .scratch/memcap.sh 290; peaks in the logs).
Notation as B3Cert FULLPROOF (= B3F) §1: H = T^θ, K = T/H, τ = T/H², y = T^ν, (p, q) = 1, N0 = pqT/H², X0 = 2√(2π)√(nT/(pq)), r = n/N0, c_r = 2√(2πr).
Dual total (E1 §2, B3F §4 last item): 𝔇 = Σ_{g,p,q} b_{gq}b_{gp}(4π/(gpq))(q/p)^{−z/2} S(p,q), S(p,q) = Σ_n n^{−z/2}σ_{−α,β}(n)e(±p̄n/q)𝒯(n), 𝒯 = 𝒯_+ + 𝒯_−.

## 1. Where Corollary B′ discards cancellation (B3F l. 55–63, l. 101–102)
Three absolute-value steps, each a separate loss:
- (P) pointwise: |𝒯(n)| ≤ 2C1·HX0^{−1/2} (Lemma B(i)). Sharp: F1 (B3F §6) max|𝒯|/(HX0^{−1/2}) = 0.705, |𝒯_±| up to 1/(2√2); here max|𝒯| = 0.52–0.66, max|𝒯_σ| = 0.3530 (logs/b3nu_Tcheck.log, b3nu_lemS.log). No loss to recover.
- (N) the n-sum per (p, q): Σ_n d(n)|𝒯(n)| ≍ H(pq/T)^{1/4}N0^{3/4}; the twist e(±p̄n/q) and the phase of 𝒯 are dropped. 𝒯_σ(n) = HX0^{−1/2}[W_σ(n)e^{iΦ_σ(n/pq)} + O(T^{Cε}/K)]
  with the explicit phase Φ_σ(x) = T[log(1 + v) + v/(1 + v)], v² = (2πx/T)(1 + v), sign v = σ (the 2D stationary value of B3F's joint phase; Lemma S, §3),
  Φ_σ(x) = σ·2√(2πTx) − πx + O(x^{3/2}T^{−1/2}): a √n-phase of size F ≍ T/H at n ≍ N0, plus a linear phase.
- (G) the (g, p, q)-sum: the Möbius-signed weights b_{gq}b_{gp}/(gpq) are summed absolutely (y² terms of equal size T^{1/2}H^{−1/2}/g each).
At the critical length the n-sum sits exactly at the edge of van der Corput's second-derivative scale: Σ_{n≤N0}e(√(2Tn/(πpq))) has N0·λ2^{1/2} = T^{1/2}H^{−1/2}
(λ2 = |f''(N0)|), which equals N0 iff pq = H^{3/2}T^{−1/2} ⟺ ν = (3θ − 1)/4 at pq = y².

## 2. Measurement (task 1): code/b3nu_T.py, b3nu_Tcheck.py, b3nu_dual.py; logs/b3nu_dual.log, b3nu_dual_040.log
Exact 𝒯 (B3Cert's FFT profile, zero shifts so B = Y_0, Hankel series to X^{−4}); 𝒯_σ e^{−iΦ_σ} x^{1/4} is smooth: a cubic spline on 1500 log-nodes reproduces the exact
quadrature to 3·10⁻¹² of HX0^{−1/2} (12 random x per case; T = 10⁶, θ = 0.40, 0.45, 0.47). Mollifier: Levinson P(x) = x, b_h = μ(h)log(y/h)/log y, all (g, p, q),
all n ≤ 4N0 (tail r ≥ 4 ≤ 5.4·10⁻⁴, B3F §6). Ratios: Rn = Σ|w||S|/Σ|w|Σ_n d|𝒯| (n-sum), Rpq = |𝔇|/Σ|w||S| ((g,p,q)-sum), Rtot = |𝔇|/Σ|w|Σ_n d|𝒯|.
| T, θ (y_crit) | y | Σ|·|/(√πH) | Rn (random-sign model) | Rpq | Rtot | |𝔇|/(√πH) |
|---|---|---|---|---|---|---|
| 10⁶, 0.45 (3.4) | 4 | 0.60 | 0.455 (0.359) | 0.050 | 0.023 | 0.014 |
| 10⁶, 0.45 | 12 | 2.79 | 0.243 (0.226) | 0.059 | 0.014 | 0.040 |
| 10⁶, 0.45 | 30 | 12.5 | 0.126 (0.127) | 0.003 | 4·10⁻⁴ | 0.005 |
| 10⁶, 0.40 (2.0) | 4 | 2.69 | 0.229 (0.197) | 0.019 | 0.004 | 0.012 |
| 10⁶, 0.40 | 12 | 11.5 | 0.135 (0.127) | 0.046 | 0.006 | 0.071 |
| 10⁶, 0.40 | 30 | 47.4 | 0.066 (0.073) | 0.001 | 5·10⁻⁵ | 0.002 |
| 10⁶, 0.47 (4.1) | 4 | 0.31 | 0.527 (0.433) | 0.078 | 0.041 | 0.013 |
| 10⁶, 0.47 | 30 | 7.13 | 0.150 (0.153) | 0.004 | 5·10⁻⁴ | 0.004 |
| 10⁷, 0.45 (4.1) | 4 | 0.37 | 0.349 (0.290) | 0.164 | 0.057 | 0.021 |
| 10⁷, 0.45 | 30 | 8.95 | 0.116 (0.110) | 0.021 | 0.002 | 0.022 |
(The twist sign e(−p̄n/q) gives the same picture: Rtot 4·10⁻⁴–0.05. Frozen SL/EXP-010 mollifier P = S(r_Px)/S(r_P), r_P = 0.7707, logs/b3nu_dual_P.log,
T = 10⁶ (θ = 0.40, 0.45), 10⁷ (θ = 0.45), y = 4, 12, 30: Rn 0.07–0.46 (random model 0.08–0.36), Rpq 3·10⁻⁴–0.19, Rtot 3·10⁻⁵–0.066, |𝔇| ≤ 0.063√πH.) Findings: the absolute sum is NOT attained. (N): the n-sum cancels like random signs
(Rn within 30% (y = 4) and 10% (y ≥ 12) of the random-sign value (Σd²|𝒯|²)^{1/2}/Σd|𝒯| ∝ N0^{−1/2}); (G): the Möbius-weighted (g,p,q)-sum cancels a further 6–1000×; |𝔇| ≤ 0.08√πH for 4 ≤ y ≤ 30
(ν ≤ 0.25, twice the critical ν) while Σ|·| grows like y². No barrier in the truth; the loss is in the tools. Only y = 2 (one term (1,1,1), n ≤ 4τ) has Rtot = Rn = 0.28–0.81.

## 3. Lemma S — phase structure of the exact-profile dual transform (new; numerical falsification test passed, §3.1)
Setting: B3F §1–2 (exact profile A_{−1}, shift disk, coprime p, q, pq ≤ y²); θ ∈ (1/3, 1), 0 < ν < θ/2, ε > 0 small; C denotes absolute constants.
By B3F Step a, 𝒯_σ(n) = κ_σπ^{−1/2}(T/2π)X0^{−1/2}𝒥_σ, 𝒥_σ = ∫∫f(Λ)α(u)e^{iKP(Λ,u)}du dΛ, P(Λ, u) := ℓ(Λ) − Λu + σc√u, c := c_r (KP = 2Ty/H − λ(y)ξ + σX exactly).
**Lemma S.** For 1 ≤ n ≤ N0T^ε and σ = ±1: 𝒯_σ(n) = HX0^{−1/2}[e^{iΦ_σ(n/pq)}W_σ(n) + E_σ(n)], |E_σ(n)| ≪_{ε} T^{Cε}K^{−1}, where
  Φ_σ(x) := T[log(1 + v) + v/(1 + v)], v = v_σ(x) the root of v² = (2πx/T)(1 + v) with sign σ;
  W_σ(n) := κ_σπ^{−1/2}|det P″(x*)|^{−1/2}f(Λ*)α(u*), x* = (Λ*, u*) the unique critical point of P: u* = 1/(1 + Λ*/H), Λ* = σc/(2√u*) (so Λ* = Hv);
  |W_σ| ≤ 2^{−3/2}e^{C0}(1 + o(1))e^{−y*²}, y* = (H/2)log(1 + v), and (n∂_n)^jW_σ(n) ≪_j T^{Cjε} (review repair: shift factor |(T/2π)^{−d/2}| ≤ e^{C0}, not 1 + O(1/L)).
Proof. (1) Critical point and value. ∇P = (ℓ′(Λ) − u, σc/(2√u) − Λ), ℓ′ = (1 + Λ/H)^{−1}; ∇P = 0 ⟺ u = (1 + Λ/H)^{−1}, Λ = σc/(2√u). With v := Λ/H: v = σc√(1+v)/(2H),
  v² = c²(1 + v)/(4H²) = (2πn/(pqT))(1 + v) (c² = 8πn/N0, N0 = pqT/H²). Value: σc√u* = 2Λ*u*, so KP(x*) = KH log(1 + v) + KΛ*u* = T log(1 + v) + Tv/(1 + v) = Φ_σ(n/pq).
(2) Nondegeneracy, uniformly. P_ΛΛ = −H^{−1}(1 + Λ/H)^{−2}, P_Λu = −1, P_uu = −σc/(4u^{3/2}); det P″ = σc/(4Hu^{3/2}(1 + Λ/H)²) − 1. Cut f smoothly to |y| ≤ L (rest ≪ Te^{−L²/2}
  ≪ T^{−A}; the cutoff amplitude at the saddle differs from f(Λ*) by e^{−Ω(L²)}); on the support u ∈ [1/2, 2], |Λ| ≤ 3L ≤ H/2, c ≤ 2√(2π)T^{ε/2}: |det P″ + 1| ≤ 2^{5/2}c/H = o(1), signature 0. The segment-averaged Hessian M = ∫₀¹P″(x* + t(x − x*))dt
  has the same form (|M_ΛΛ| ≤ 4/H, |M_uu| ≤ 2^{−1/2}c, M_Λu = −1), so det M = −1 + o(1), ‖M^{−1}‖ ≤ 2(1 + c) and |∇P(x)| = |M(x − x*)| ≥ |x − x*|/(2(1 + c)): x* is the only critical point.
(3) Stationary phase (Hörmander, ALPDO I, Thm 7.7.5, n = 2, k = 2; its error is Cω^{−k}Σ_{|γ|≤2k}sup|∂^γa|, so k = 2 is needed): ∫∫a e^{iKP} = (2π/K)|det P″(x*)|^{−1/2}e^{iKP(x*)}
  [a(x*) + K^{−1}L_1a] + O(K^{−2}Σ_{|γ|≤4}sup|∂^γa|), C bounded while P stays in a bounded C⁷-set and sup|x − x*|/|∇P| is bounded. Uniformity: partition |Λ| ≤ 3L into unit pieces (non-stationary
  pieces: Thm 7.7.1, O(K^{−N})); rescale Λ − Λ* = R·a, R := 1 + c + L, so KP = (KR)·P/R with P/R in a bounded C⁷-set on a bounded (a, u)-support and Hessian det in (a, u) equal to det P″ ≈ −1; amplitude derivatives ≪ (1 + c)^C L^C
  (∂_Λ^i f ≪ (1 + |Λ|)^i f-envelope, ‖α^{(j)}‖ ≤ A_j, B3F Step a). All constants ≪ T^{Cε}; |L_1a| ≪ T^{Cε}.
  a = fα has a(x*) = f(Λ*)α(u*) (χ̃(u*) = 1 as u* = 1 + O(T^ε/H)). Multiply by κ_σπ^{−1/2}(T/2π)X0^{−1/2}, (T/2π)(2π/K) = H: the formula with |E_σ| ≪ T^{Cε}K^{−1}.
(4) W_σ: |κ_σ|π^{−1/2} = 2^{−1/2}, |f(Λ*)| = ½e^{−y*² −Re(1+2β)y*/H}, |α(u*)| = u*^{−1/4}|(Tu*/2π)^{−d/2}||a_σ| ≤ e^{C0}(1 + o(1)) (−β − z/2 = −d/2), a_σ = 1 + O(1/X0); x* is a smooth function of c (implicit
  functions, det P″ ≈ −1) with ∂_c^j x* ≪ (1 + c)^{Cj}; n∂_n = (c/2)∂_c + (X0/2)∂_{X0} and X^k∂_X^k a_σ = O(1) (E1A §2) give (n∂_n)^jW_σ ≪ T^{Cjε}. ∎
Remark. Lemma B's |𝒯_σ| ≤ C1HX0^{−1/2} is the size of this leading term; Lemma S adds its phase Φ_σ, which Corollary B′ discarded.
### 3.1 Test (code/b3nu_lemS.py, logs/b3nu_lemS.log; exact FFT profile, zero shifts, r = x/τ ∈ [10⁻³, 3], 60 points, both σ)
| T, θ | K | max|𝒯_σ − lead_σ|/(HX0^{−1/2}) | K·err |
|---|---|---|---|
| 10⁵, 0.45 | 562 | 9.2·10⁻⁵ | 0.05 |
| 10⁶, 0.45 | 1995 | 2.6·10⁻⁵ | 0.05 |
| 10⁷, 0.45 | 7079 | 7.3·10⁻⁶ | 0.05 |
| 10⁶, 0.40 | 3981 | 1.3·10⁻⁵ | 0.05 |
| 10⁶, 0.49 | 1148 | 4.5·10⁻⁵ | 0.05 |
| 10⁷, 0.49 | 3715 | 1.4·10⁻⁵ | 0.05 |
The error is exactly ≍ 0.05/K (Lemma S's E_σ), max|𝒯_σ|/(HX0^{−1/2}) = 0.3530 → 2^{−3/2} = 0.35355 (W_σ at y* → 0). Also 𝒯_σe^{−iΦ_σ}x^{1/4} is a smooth function of log x:
a cubic spline on 1500 nodes reproduces 𝒯 to 9·10⁻¹³ (logs/b3nu_Tcheck.log). Complex shifts (code/b3nu_lemS2.py, B3Cert's b3c_fals.py imported read-only:
complex β in the profile, a_± for complex d; shifts (0.05 + 0.1i, −0.03 + 0.07i), (0.2 − 0.1i, 0.03 + 0.2i); T = 10⁶, 10⁷, θ = 0.45, 0.49; 24 r-points): K·err = 0.013–0.049, flat in T
(logs/b3nu_lemS2.log; the lead term carries e^{−y*² + (1−2β)y*/H}ξ*^{−β−z/2}a_σ(X*, d)).

## 4. Lemma C — exponent-pair dual bound (task 2; new)
Exponent pairs (k, l): Graham–Kolesnik (GK), *Van der Corput's method of exponential sums* (1991), Ch. 3: Σ_{v∈J}e(f(v)) ≪ L^kV^l + L^{−1} for J ⊂ (V, 2V] and f in GK's
monomial class (f^{(r)}(v) = (1 + O(ε0))·c(−y)_r…: |f^{(r)}| ≍ LV^{1−r}, r ≤ r0, with f′ ≍ L). Process A: (k′, l′) ↦ (k′/(2k′+2), (k′+l′+1)/(2k′+2)); B: (k, l) ↦ (l − ½, k + ½).
Linear-phase immunity: for (k, l) = A(k′, l′) and L ≥ 1 (as here) the bound holds for F = f + λv, any real λ, with f in the class for 2 ≤ r ≤ r0 only. Proof: Weyl–van der Corput inequality (GK Ch. 2),
|Σ_J e(F)|² ≤ (1 + |J|/Q)Σ_{|h|<Q}(1 − |h|/Q)Σ_{v,v+h∈J}e(F(v+h) − F(v)); F(v+h) − F(v) = [f(v+h) − f(v)] + λh, and the constant λh leaves |Σ_v| unchanged; the
A-process proof (GK Ch. 3) then uses only g_h = f(·+h) − f, whose r-th derivatives are h·f^{(r+1)}(1 + O(h/V)), r ≥ 1. Pairs used (classical): (1/30, 13/15) = A³B(0, 1) = A(1/14, 11/14) and (13/502, 445/502) = A³BA²BA²B(0, 1) = A(13/238, 97/119) (§5).
**Lemma C.** Let (k, l) = A(k′, l′) be an exponent pair, θ ∈ (1/3, 1), 0 < ν < θ/2. Uniformly in g, coprime p, q (gp, gq ≤ y), shifts in the disk and the twist sign:
  S(p, q) := Σ_n n^{−z/2}σ_{−α,β}(n)e(±p̄n/q)𝒯(n) ≪_ε T^{Cε}[(pq)^{(1+l−k)/2}T^{(k+l)/2}H^{1/2−l} + H(pq/T)^{1/4}N0^{3/4}K^{−1} + pqH^{1/2}T^{−1/2}] + T^{−A},
and the oscillatory dual total 𝔇 ≪ T^{Cε}[y^{1+l−k}T^{(k+l)/2}H^{1/2−l} + y²T^{−1/2}H^{1/2}] + T^{−A}. Relative to the main term √πH·c this bound is power-small iff ν < ν1(θ),
  ν1(θ) := [θ(1 + 2l) − (k + l)]/(2(1 + l − k)),   ν1(θ) − (3θ − 1)/4 = (1 − θ)(1 − l − 3k)/(4(1 + l − k));   (1/30, 13/15): ν1(θ) = (82θ − 27)/110 = (3θ − 1)/4 + (1 − θ)/220.
Proof. (a) n > N0T^ε: B3F Cor B′ second range (Lemma B(ii)): ≪ T^{−A}. (b) n ≤ N0T^ε: insert Lemma S for σ = ±1. E_σ-part, summed absolutely as in Cor B′ (Σ_{n≤N0T^ε}d(n)n^ε(pq/(nT))^{1/4}):
  ≪ T^{Cε}H(pq/T)^{1/4}N0^{3/4}K^{−1}. (c) Main part, dyadic n ∈ (N, 2N], N ≤ N0T^ε. σ_{−α,β}(n) = Σ_{uv=n}u^{α}v^{−β} (E1A convention); split u ≤ v / u > v (symmetric: swap roles).
  For fixed u ≤ √(2N), v runs over an interval J_u ⊂ (N/u, 2N/u], V_u := N/u ≥ √(N/2). Summand: g_u(v)e(F_u(v)), g_u(v) := u^{α−z/2}v^{−β−z/2}HX0(uv)^{−1/2}W_σ(uv),
  F_u(v) := f_u(v) + λ_uv, f_u(v) := Φ_σ(uv/pq)/(2π), λ_u := ±(p̄u mod q)/q. sup|g_u| + Var_{J_u}g_u ≪ T^{Cε}H(pq/(NT))^{1/4} (Lemma S: (n∂_n)W ≪ T^{Cε}; n∂_nX0^{−1/2} = −X0^{−1/2}/4;
  |v∂_v v^{−β−z/2}| ≪ 1/L). Partial summation: |Σ_{J_u}g_ue(F_u)| ≤ (sup|g_u| + Var g_u)·max_{J′⊂J_u}|Σ_{J′}e(F_u)|.
  Phase class. Φ_σ(x) = Tφ_σ(x/T), φ_σ(t) = log(1 + v) + v/(1 + v) with v analytic in s = √t (v = σ√(2π)s + πs² + …): φ_σ(t) = 2σ√(2πt) − πt + t^{3/2}ψ_σ(√t), ψ_σ analytic near 0.
  On J_u, t = uv/(pqT) ≤ 2T^ε/H², so for 2 ≤ r ≤ r0: f_u^{(r)}(v) = σ√(2Tu/(πpq))·(½)_r v^{1/2−r}(1 + O_r(T^{ε}/H)), (½)_r := ½(½ − 1)…(½ − r + 1) ≠ 0; the −πt term is linear in v.
  So f_u is in GK's class for r ≥ 2 with L_u := √(Tu/(pqV_u)) ≍ uF_N/N, F_N := √(TN/pq), and L_u ≥ √(T/(pqN)) ≥ H/(pqT^{ε/2}) ≥ T^{θ−2ν−ε} ≥ 1. Linear-phase immunity:
  |Σ_{J′}e(F_u)| ≪ L_u^kV_u^l + L_u^{−1} = F_N^k(N/u)^{l−k} + N/(uF_N). Sum over u ≤ √(2N) (k − l > −1; for (0, 1) one extra L): Σ_u(N/u)^{l−k} ≪ N^{(1+l−k)/2}; Σ_u N/(uF_N) ≪ (N/F_N)L.
  Hence the dyadic block is ≪ T^{Cε}H(pq/(NT))^{1/4}[F_N^kN^{(1+l−k)/2} + (N/F_N)L] = T^{Cε}[H(pq/T)^{1/4}(T/pq)^{k/2}N^{1/4+l/2} + H(pq)^{3/4}N^{1/4}T^{−3/4}].
  Both increase with N; at N = N0T^ε (N0 = pqT/H²): H(pq/T)^{1/4}(T/pq)^{k/2}N0^{1/4+l/2} = (pq)^{(1+l−k)/2}T^{(k+l)/2}H^{1/2−l} and H(pq)^{3/4}N0^{1/4}T^{−3/4} = pqH^{1/2}T^{−1/2} (review repair: was T^{−1}). O(L) blocks.
(d) Total. Weights |b_{gq}b_{gp}(4π/(gpq))(q/p)^{−z/2}| ≪ 1/(gpq); Σ_{g≤y}g^{−1}(Σ_{p≤y/g}p^{(l−k−1)/2})² ≪ Σ_g g^{−1}(y/g)^{1+l−k} ≪ y^{1+l−k} (0 < 1 + l − k < 2); the E-term gives
  Σ g^{−1}T^{1/2}H^{−1/2}K^{−1} ≪ y²T^{−1/2}H^{1/2}L, the last term the same y²H^{1/2}T^{−1/2}L. Divide by H: T^{(1+l−k)ν + (k+l)/2 − (1/2+l)θ} (< 1 ⟺ ν < ν1(θ)) and T^{2ν − (1+θ)/2} (< 1 since
  ν < θ/2 < (1 + θ)/4). ∎
Checks (algebraic consistency, outside the hypotheses θ ∈ (1/3, 1)). θ = 1: ν1(1) = 1/2 for every pair = Balasubramanian–Conrey–Heath-Brown's range (no 4/7: the (g, p, q)-sum is still absolute, Conrey's Kloosterman step not used).
ν → 0, p = q = 1: θ > (k + l)/(1 + 2l) = 27/82 for (1/30, 13/15) — Karatsuba's short-interval exponent for critical zeros (Izv. 1984, H = T^{27/82+ε}), the classical B(1/30, 13/15) = (11/30, 16/30)
exponent of E(T); the trivial pair (0, 1) returns Cor B′ exactly (gain factor 1 − l − 3k = 0). Gain needs l + 3k < 1: (1/2, 1/2), (1/6, 2/3) give none, (1/14, 11/14) zero, A³B(0, 1) gives (1 − θ)/220.
Where E1′ changes: only B3F §4 last bullet (E1 Step 5 with Cor B′) is replaced by Lemma C(d); every other input keeps its range (ν < 1/2, θ < 1; X0 ≥ 1 needs ν < 1/2; Lemma B(i) needs c_r ≤ H/46).

## 5. Theorem E1″ and the certificate (task 3)
**Theorem E1″ (moment).** SL Thm 1.2 (Gaussian-weighted shifted mollified second moment on windows T^θ, any fixed real Q, Q(0) = 1, error O(H/L)) holds for every θ ∈ (1/3, 1) and
0 < ν < ν1(θ), for every exponent pair (k, l) = A(k′, l′) (ν1 < θ/2 on (1/3, 1)). Proof: B3F §§1–5 verbatim with §4's last bullet replaced by Lemma C(d).
Pairs used: (13/502, 445/502) = A³BA²BA²B(0, 1) (classical; best A-type pair of the (0, 1) tree to depth 12 for the gain (1 − l − 3k)/(1 + l − k) = 0.019272, logs/b3nu_pairs.log):
ν1(θ) = (696θ − 229)/934 = (3θ − 1)/4 + 9(1 − θ)/1868; and the simplest (1/30, 13/15) = A³B(0, 1): ν1(θ) = (82θ − 27)/110 = (3θ − 1)/4 + (1 − θ)/220.
**Theorem (simple zeros).** liminf_{T→∞} S(T, T^θ)/N(T, T^θ) ≥ h(θ_c; κ) > 0 for every fixed θ ∈ (θ_c, 1), θ_c := 24699/50000 = 0.49398 (before: 371/750 = 0.494666…, B3Cert);
with (1/30, 13/15) alone: θ > 24701/50000 = 0.49402.
Proof: B3F §5 with ν0 := ν1(θ_c): for θ > θ_c take η = (θ − θ_c)/2, θ′ = θ − η > θ_c, so ν0 < ν1(θ′) (ν1 increasing) and the Gaussian moments of width T^{θ′} hold; chain
(Emean) → Lemma 4.2 → Thm 1.3 → Thm 5.1 unchanged (θ-free on (0, 1)); detector = SL/EXP-010 frozen (P, Q), R = 44/5 as B3Cert; k* decreasing, so h(θ; κ) ≥ h(θ_c; κ) > 0.
Arb (code/b3nu_arb.py; construction executed verbatim from code/b3c_arb_copy.py = B3Cert's b3c_arb.py; detector exactly admissible: P(0) = 0, P(1) = 1, Q(0) = 1,
Q(y) + Q(1 − y) ≡ 1, deg Q = 201; check E1's κ1 = 0.0844779079697484032032447429603 ± 4·10⁻³²):
| pair | θ_c | ν0 = ν1(θ_c) | κ(ν0, 44/5) | k*(θ_c) | κ − k* | h(θ_c; κ) | logs |
|---|---|---|---|---|---|---|---|
| (13/502, 445/502) | 24699/50000 | 717563/5837500 = 0.1229230 | 0.08609149289632616 | 0.08607188418093183 | 1.96087·10⁻⁵ | 1.39047·10⁻⁵ | b3nu_arb_b.log (4000 bits), b3nu_arb_b3000.log (3000): identical |
| (1/30, 13/15) | 24701/50000 | 337741/2750000 = 0.1228149 | 0.08603324712433759 | 0.08600892477917683 | 2.43223·10⁻⁵ | 1.72464·10⁻⁵ | b3nu_arb.log (4000), b3nu_arb_3000.log (3000): identical |
Radii ≤ 2·10⁻²⁵; ν0 − (3θ_c − 1)/4 = (1 − θ_c)(1 − l − 3k)/(4(1 + l − k)) and ν0 < θ_c/2 checked exactly. Control: Cor B′'s range ν = (3θ_c − 1)/4 gives κ − k*(θ_c) = −1.32·10⁻³
and −1.24·10⁻³ < 0 (θ_c is out of B3Cert's reach).
Float optimum (logs/b3nu_thr.log, bisection of κ(ν1(θ), 44/5) = k*(θ)): (0, 1): 0.494642881; (1/30, 13/15): 0.494007692; (13/502, 445/502): 0.493970073 (ν1 = 0.122916);
A²(Bourgain 13/84 + ε, 55/84 + ε) = (13/414, 359/414)+: 0.493908737 (ν1 = 0.123095; best from Bourgain/Huxley seeds; not certified: Bourgain's pair carries his own hypotheses on f,
unchecked here). κ′(ν) = 0.552, dk*/dθ = −1.573 near θ = 0.4947. The detector (frozen R = 44/5) was not re-optimised.

### 5.1 Range audit (orchestrator referee GAP: "§§1–5 verbatim" asserted, not shown) — no range fails, no number changes
Standing facts: ν1(θ) < θ/2 on (1/3, 1) for both pairs ((696θ − 229)/934 < θ/2 ⟺ 229θ < 229; (82θ − 27)/110 < θ/2 ⟺ 27θ < 27). So every ν ≤ ν1(θ) has ν < 1/2,
pq ≤ y² = T^{2ν} < H ≤ T and X0 ≥ 2√(2π)T^{1/2−ν} → ∞. Certified use: ν0 = 717563/5837500 = 0.12292 at θ′ = θ − η ∈ (θ_c, 1), θ_c = 0.49398 (§5). Sources: B3F = B3Cert/FULLPROOF.md,
E1 = SantiAudit/E1.md, E1A = E1Audit/FULLPROOF.md, B3A = B3Audit/FULLPROOF.md, SL = main.tex of Zenodo 23135354, Wang = arXiv 2609.07918 TeX (line numbers re-read today).
| # | import (source lines) | needs (ν, θ) | holds for ν ≤ ν1(θ) < θ/2, θ ∈ (1/3, 1) |
|---|---|---|---|
| 1 | Gaussian reduction (10), exact kernel K_β = Σ_ε e(εx)A_ε (SL l. 357–390; B3F l. 7–13) | none (identity, every x > 0) | ✓ |
| 2 | Estermann FE (12), dual Mellin form (14), residue formula (SL l. 434–475; E1A l. 8–22, numerical identity l. 87) | coprime p, q; shift disk \|d\| ≤ 1/4 (L > 8R) | ν-, θ-free ✓ |
| 3 | E1 Steps 1–2, x-space + Hankel forms (E1 l. 23–34; E1A l. 24–37; B3F l. 100) | A ∈ C_c^∞(0, ∞); X^k∂^k a_± = O(1) for all X > 0 | ν-, θ-free ✓ |
| 4 | X0 ≥ 1, hypothesis of Lemmas B, S (B3F l. 5, 17) | ν < 1/2 | ν < θ/2 < 1/2 ✓ |
| 5 | Lemma B(i) (B3F l. 18, proof l. 21–51) | c_r ≤ H/46, 1 ≤ H ≤ T | used for n ≤ N0T^ε: c_r ≤ 2√(2π)T^{ε/2} ≤ H/46 for T ≥ T0(ε, θ); ν-free ✓ |
| 6 | Lemma B(ii) + Cor B′ tail n > N0T^ε (B3F l. 19, 52–63) | N0 ≤ T² (pq ≤ T, H ≥ 1); N ≥ 2(A + 5)/ε | pq ≤ T^{2ν} < T ✓ |
| 7 | Lemma S (§3 here) | n ≤ N0T^ε; c ≤ 2√(2π)T^{ε/2}; 3L ≤ H/2; K = T^{1−θ} → ∞; X0 ≥ 1 | θ ∈ (1/3, 1), ν < 1/2 ✓ |
| 8 | Lemma C main term (§4(c)) | L_u ≥ T^{θ−2ν−ε} ≥ 1: ν < θ/2, ε < θ − 2ν | ν1(θ) < θ/2; fix ε after (θ, ν) ✓ |
| 9 | Lemma C side terms (§4(d)) | T^{2ν−(1+θ)/2} < 1: ν < (1 + θ)/4 | ν < θ/2 < (1 + θ)/4 ✓ |
| 10 | Main term: SL Lemma diagonal (SL l. 285–309; Remark l. 336–345) | ν < 1/2 (l. 304–305, 338); θ < 1 (supp w ⊂ [T/2, 2T], l. 339) | ✓ |
| 11 | Residues = Young/Conrey off-diagonal main term (SL l. 444–449, 709–729; B3F Lemma R l. 67–87) | ν < 1/2 (SL l. 728); \|w − 1\| ≤ 1/4; change O(H²T^{−1}L^{O(1)}) = o(H/L) ⟸ θ < 1 | ✓ |
| 12 | Off-band (1 − χ)K_β (SL l. 389–391; B3F l. 89–91; B3A 3.1) | O(T^{1+ν+ε}(H/T)^A) ≪ HT^{−δ}: A(1 − θ) > 1 + ν − θ + δ; tail e^{−H²/400} uses θ ≥ 1/3 | fixed A for θ < 1, ν < 1/2; θ > 1/3 ✓ |
| 13 | Representation error SL (11) (SL l. 395–398; B3F l. 92–93; B3A 3.2) | Taylor part absent (exact profile); (H/T)^A part as row 12 | ✓ |
| 14 | ε = +1 term (SL l. 541–543; B3F l. 94–96; B3A 3.3) | \|A_{+1}\| ≪ (H/T)^A, absolute sum O(T^{1+ν+ε}(H/T)^A) | as row 12 ✓ |
| 15 | K-Bessel branch (SL l. 539–542; B3F l. 97–99; B3A 3.4) | T(y²/T)^σ ≪ T^{−A}: ν < 1/2, σ > (1 + A)/(1 − 2ν) | ✓ |
| 16 | SL §2 hypothesis 0 < ν < θ − 1/2 (SL l. 199; Remark l. 337) | only in Lemma twisted, replaced by the §3 signed route (B3A 4.4) | not used ✓ |
| 17 | SL l. 536–538 (N0 < 1 at p = 1 or q = 1) | θ > 1/2 | not used: Lemma C covers p = 1 (p̄ = 1), q = 1 (no twist) and N0 ≥ 1 ✓ |
| 18 | Shift derivatives (SL l. 723–726; E1A l. 63) | Cauchy circles of radius ≍ 1/L: L^{O(1)}; all bounds uniform in the disk | Lemmas S, C uniform in the disk ✓ |
| 19 | Compact windows (SL l. 731–758; B3F l. 106–110; B3A 4.1) | moment at width s0 = T^{θ′}, centres T + O(H): ν < ν1(θ′); D > (1 + 2ν − θ + δ)/(2η) | η = (θ − θ_c)/2: ν0 = ν1(θ_c) < ν1(θ′) (ν1 increasing), θ′ > θ_c > 1/3; S, C uniform in T′ ≍ T ✓ |
| 20 | SL l. 731–732 "θ − η > 1/2" | SL's own Gaussian range (Taylor profiles, j ≥ 1 residues, N0 < 1) | not used (B3A 4.1) ✓ |
| 21 | (Emean) (SL l. 812–815) | Thm 1.2 for Q = 1, 1 + x^j, same P, R, ν | from E1″ ✓ |
| 22 | Lemma count, Thm 1.3, signs (SL l. 821–928) | θ > 0; constants depend on P, Q, R, ν (l. 921–923) | fixed ν0 ✓ |
| 23 | Thm D transfer (SL l. 930–955): Wang thm:local-pair (Wang l. 231–241), CAOSShort Prop 2.5 | 0 < λ < θ < 1, supp g ⊂ [−λ, λ], error O_g(H + T^λL²); parity inequality deterministic | ν-free; λ ↑ θ ✓ |
| 24 | c(θ), c′(θ) = ½cot²(θ/√2) > 0, monotonicity of h (Wang l. 125, 129; SL l. 957–961) | 0 < θ < 1 | ✓ |
| 25 | Detector κ = 1 − R^{−1}log c(P, Q, R, ν) (SL Remark l. 343; B3F l. 115–117) | Thm 1.2 at that ν; (P, Q) admissible | at ν0, Arb §5 ✓ |
Only B3F l. 101–102 (E1 Step 5 with Cor B′) and its window carry-over (B3F l. 108) used ν < (3θ − 1)/4; both are replaced by Lemma C(d) at θ and θ′.

## 6. Review (reviewer subagent; statement and proof of §§3–4 pasted inline with the audited context; .scratch/b3nu/review1.txt, review2.txt)
- Pass 1: GAP — (1) |W_σ| bound: the shift factor |(Tu/2π)^{−d/2}| is O(e^{C0}), not 1 + O(1/L); (2) Lemma C(c): H(pq)^{3/4}N0^{1/4}T^{−3/4} = pqH^{1/2}T^{−1/2}, not T^{−1}.
  Both repaired (marked "review repair"); neither changes ν1. Reviewer also confirmed: A-process immunity valid; partial summation with W_σ legitimate; Hörmander's O(ω^{−2}) suffices.
- Pass 2 (repaired text, plus §5's θ_c/η logic): **VERDICT: CORRECT** — "no gaps … Pair: A³BA²BA²B(0,1) = (13/502,445/502) = A(13/238,97/119), a classical A-type exponent pair applicable
  here … ν1 = (696θ − 229)/934 … ν0 = ν1(θ_c) < ν1(θ′)". Minor points applied: immunity stated for L ≥ 1; rescaling R = 1 + c + L; theorem stated on (θ_c, 1); "this bound gives a power
  saving iff". Not re-derived by the reviewer: the Arb logs (supplied as evidence; rerun identical at 3000/4000 bits here).

## 7. Novelty (adversarial)
arXiv listings math.NT/CA/CV 2026-10, 2026-11 (.scratch/lst.py, regex exponent pair|short interval|Levinson|mollif|simple zero|zeros of the Riemann|twisted; logs/novelty.log): no short-window
mollified/twisted moment below √T besides our sources; 2608.14852 (short twisted 2nd moment, reciprocity) needs δ > 1/2 and single twists (its Remark 1 expects δ > 1/3 via the quadratic
Taylor term); 2610.26841 long ranges. Zenodo ("short intervals" "simple zeros"; "exponent pair" zeta mollifier): only SL 23134787/23135354. arXiv API (abs: Levinson/mollified + short intervals,
exponent pair + mollifier): nothing. Known ingredients (not new): exponent pairs on the √n-phase dual are classical (Karatsuba 1984 critical zeros in (T, T + T^{27/82+ε}] by Selberg's
method; E(T) ≪ T^{27/82} via (11/30, 16/30)); Lemma S is standard stationary phase. New here: the twisted exact-profile dual (twist removed by the A-process uniformly in p̄/q), the joint range
ν1(θ) for Levinson/Conrey mollifiers of length T^ν in windows T^θ, and the certified simple-zero threshold below 371/750. Class: RECORD (constant/range improvement of B3Cert's node).

## 8. What remains (barrier analysis; logs/b3nu_stakes.log)
Lemma C keeps two absolute sums: (U) the divisor split's short variable u ≤ √N (worst at u ≈ v ≈ √N: 1D pairs are weakest there; 2D methods — Kolesnik/Huxley for E(T),
131/416 < 27/82 — would need uniformity in the additive twist mod q); (G) the (g, p, q)-sum (numerics §2: a further 6–1000× cancellation; needs Kloosterman-fraction bilinear bounds
à la Deshouillers–Iwaniec/Conrey (structured b) or Duke–Friedlander–Iwaniec/BCR (general b), transferred to the oscillatory weight e^{iΦ(n/pq)}).
Stakes with the frozen detector (Arb floats, R = 44/5): gain factor g in ν = (3θ − 1)/4 + g(1 − θ): g = 0.0048 (ours) → θ = 0.49397; g = 0.01 → 0.49326; g = 0.05 → 0.48846; g = 0.1 → 0.48400;
the full Lemma-C range ν = θ/2 → 0.47912; fixed ν = 0.25 → 0.47918 (the frozen (P, Q, R) saturates: κ(0.2) = 0.1078, κ(0.25) = 0.1095). Below ≈ 0.479 the detector, not the dual, binds.

## Suggested NOGO/TREE line
TREE N8 (short windows): moment range extended past Cor B′ by exponent pairs on the √n-phase dual (Lemma S phase + Lemma C, A-process kills the twist): ν < (696θ − 229)/934, simple
critical zeros in (T, T + T^θ] for θ > 0.49398 (B3Nu, record). Next: (U) 2D/Huxley-type twisted divisor sums, (G) Kloosterman bilinear bounds on the (p, q)-sum; detector floor ≈ 0.479.
NOGO: none (the absolute dual sum is not attained: Rn ≈ random-sign, Rpq 10⁻³–0.16; the limit (3θ − 1)/4 was a tool loss, not a barrier).

RH relevance: N8

## Files
Riemann/Deep/B3Nu/REPORT.md (this). code/: b3nu_T.py (exact 𝒯 by B3Cert's FFT profile, Hankel branches, demodulated spline), b3nu_Tcheck.py, b3nu_dual.py (task 1 ratios), b3nu_lemS.py,
b3nu_lemS2.py (Lemma S tests, zero and complex shifts), b3nu_pairs.py (exponent-pair search), b3nu_thr.py (float thresholds), b3nu_arb.py (Arb certificate), b3nu_stakes.py, b3c_arb_copy.py (read-only copy of B3Cert/code/b3c_arb.py).
logs/: b3nu_dual.log, b3nu_dual_040.log, b3nu_dual_P.log, b3nu_Tcheck.log, b3nu_lemS.log, b3nu_lemS2.log, b3nu_pairs.log, b3nu_thr.log, b3nu_arb.log, b3nu_arb_3000.log, b3nu_arb_b.log, b3nu_arb_b3000.log,
b3nu_stakes.log, novelty.log, footprint.log (max RSS ≤ 170 MB, memcap peaks ≤ 117 MB; one early run killed at 589 MB before a chunking fix, unused).
