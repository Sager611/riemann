# DERIVATION — selbergn4: Selberg's method for F = Σ_{j≤N} c_j ε_j q^{s/2−1/4} L(s, χ_j), explicit, re-derived and sharpened

Source: J. Dousselin, arXiv 2311.10285 (= Forum Math. 37 (2025) 1871–1912), TeX \`src/d.tex\` (line numbers below refer to it).
Setting (his §2, l.648–700): q ≥ 3, χ_1, …, χ_N distinct even primitive characters mod q (odd: "similarly", his remark l.640), c_j ∈ ℝ∖{0},
ξ = T^ϰ, H = A/log T, mollifier η_j = Σ_{n≤ξ} τ_{−1/2}(n)χ_j(n)M(n)n^{−s}, M(x) = [log(ξ/x) − log⁺(ξ^θ/x)]/((1−θ)log ξ) for x ≤ ξ (l.672–684),
I_j, J_j, M_j as in l.687–693; κ_F = liminf (N_0(2T,F) − N_0(T,F))/(N(2T,F) − N(T,F)).
"≲" means ≤ (1 + o(1)) as T → ∞ with q, χ_j, c_j, N, θ, ϰ, A fixed. Every limit below is T → ∞.

**What is reused unchanged (Dousselin's proofs, not re-derived here):** the dominance decomposition (S_j^*, l.809–885, including (I*), (J*));
Bauer's mean value theorem and its corollaries (l.1191–1362); the reduction ∫|I|² = 4T(φ(q)/q)∫_0^H (H−v) Im[K(v)/v] dv + O(T/log³T) (eqIntfacile, l.1586)
with K(v) = 𝒯^{iv/2} Σ_{d≤ξ², (d,q)=1} φ_{iv}(d) A'_d(v)² (eqExpK l.1563, eqa' l.1566); Lemma 12 (l.1107), its consequence eqintero (l.1641–1645),
the r = 2, 3 combination in the proof of coro-molli2 (l.1146–1150); eqprodeul (l.1461); the bounds for Ã_d (bd, l.1474–1484) and for Σ_{k|d}Λ(k)kφ(d/k) (l.1490–1497);
the rectangle/Gabriel argument for M (l.1363–1405). **What is new:** §1–§9 below.

## 1. Counting (replaces "pairs of intervals", l.736–741 and l.897–902)
Let G = {t ∈ (T,2T) : 𝓕(1/2+i·) changes sign in (t, t+H)} and γ_1 < … < γ_m the sign-change points of 𝓕 in (T, 2T+H) (zeros of odd order of F on the line).
If t ∈ G then t ∈ (γ_i − H, γ_i) for some i, so meas(G) ≤ mH, i.e. **m ≥ meas(G)/H = (T − meas E)log T / A**, E = (T,2T)∖G
(to count zeros below 2T restrict to t < 2T − H; this loses H in measure, i.e. O(1) zeros).
Dousselin's pairing gives m ≥ T/(2H) − meas(E)/H, i.e. half the main term. Zeros in (2T, 2T+H) are O(A) = o(T log T).

## 2. Combination per dominance set (replaces the H/(2+ε)–split, l.887–896)
For t ∈ S_j^* (Dousselin's set; Σ_j meas S_j^* = T + o(T), l.858–863) with no sign change of 𝓕 in (t, t+H), the integrand of I_j^* has constant sign, so |I_j^*| = J_j^*.
By (I*), (J*) and (min-J): |c_j||I_j| + O(H/loglogT) ≥ |c_j|J_j − O(H/loglogT) ≥ |c_j|(H − |M_j|) − O(H/loglogT), hence |I_j| + |M_j| ≥ H(1 − o(1)) on
S'_j := {t ∈ S_j^* : no sign change}. Cauchy–Schwarz on S'_j: H·meas(S'_j) ≲ ∫_{S'_j}(|I_j| + |M_j|) ≤ √meas(S'_j)·(‖I_j‖₂ + ‖M_j‖₂), so with (estiI), (estiM)
  meas(S'_j) ≲ (√C1 + √C2)² T/A²   (this is Dousselin's N = 1 argument, l.722–733, run inside S_j^*).
E ⊂ ((T,2T)∖∪S_j^*) ∪ ∪_j S'_j, so meas E ≲ N(√C1+√C2)²T/A². With §1 and (eq:N):
  **κ_F ≥ 2π (1/A − N(√C1(A) + √C2)²/A³)**   for every A > 1/ϰ, θ ∈ (0,1), ϰ ∈ (0, 1/8).    (Dousselin, Thm mainth: 2π(1/(2A) − 4N(C1+C2)/A³).)

## 3. The M-estimate (replaces the factor 6 in C2, l.1377–1405)
M(t,H) = a(t) + b(t) − c(t) (three sides of the rectangle, l.1369–1372). Minkowski instead of |a+b+c|² ≤ 3(|a|²+|b|²+|c|²): ‖a‖₂² ≪ Tξ^{−θ/2+o(1)} = o(T/log²T)
(first-esti), ‖b‖₂, ‖c‖₂ ≲ ‖Lη²−1‖₂/(θ log ξ) (his Gabriel/Cauchy–Schwarz step, l.1389–1401, with ‖Lη²(1/2+i·)−1‖₂² in place of (C3+1+2√C3)T).
First moment: ∫_T^{2T}(Lη²(1/2+it) − 1)dt = o(T). Proof: the same rectangle [1/2, 3/2] × [T, 2T] (Lη² − 1 is entire): on σ = 3/2, Lη² − 1 = Σ_{n≥ξ^θ}c(n)n^{−s}, |c(n)| ≤ τ_3(n)
(l.1380–1384), so that side is ≤ Σ_n 2τ_3(n)n^{−3/2}/log n = O(1); the horizontal sides are ≤ sup_{1/2≤σ≤3/2}(|L(σ+iT)||η(σ+iT)|² + 1) ≪ T^{1/6+ε}·4ξ = O(T^{7/24+ε})
(|η|² ≤ (Σ_{n≤ξ}n^{−1/2})² ≤ 4ξ; L(1/2+it) ≪ t^{1/6+o(1)} as used at l.1507, Phragmén–Lindelöf). Hence ∫Lη² = T + o(T) and
‖Lη²−1‖₂² = ‖Lη²‖₂² − 2Re∫Lη² + T ≤ (C3 − 1 + o(1))T. So
  **C2 = 4(C3 − 1)/(θϰ)²**   (Dousselin: 6(C3 + 1 + 2√C3)/(θϰ)²).

## 4. Exact main terms (replace Lemma A's contour bound, l.907–959)
4.1 Δ_2 in closed form. With φ(x) := −2x^{−1/2} + ∫_0^x (e^{−it}−1)t^{−3/2}dt and Δ_2(X) := ∫_0^X φ (l.976–978): φ'(x) = e^{−ix}x^{−3/2}; ψ(X) := Xφ(X) − Δ_2(X) = ∫_0^X e^{−ix}x^{−1/2}dx
  = √π e^{−iπ/4} erf(e^{iπ/4}√X); φ(x) = −2x^{−1/2}e^{−ix} − 2iψ(x); **Δ_2(X) = −2√X e^{−iX} − (1+2iX)ψ(X)**; **E(Y) := Δ_2(Y)/√Y = Σ_{j≥0} (−iY)^j/(j!(j²−1/4))**
  (his series l.1029–1040, r = 2), E(0) = −4; Δ_2(X) = c1X − ψ∞ + O(X^{−1/2}) with c1 = 2√π e^{−3iπ/4}, ψ∞ = √π e^{−iπ/4}; |R0(Y)| := |Δ_2(Y) − c1Y| ≤ √π + 4/√Y.
  Certified: series = closed form (|diff| < 10^{−30} at Y = 0.01, 0.5, 1, 2.5), = the defining double integral (quadrature, 10^{−8}), cert.log "selfcheck".
4.2 For d ≤ ξ², (d,q)=1, δ | d^∞, u := log δ/log ξ ∈ [0,1), ℓ1 := 1−u, ℓ2 := (θ−u)⁺, Π_{qd}(v) := ∏_{p|qd}(1 − p^{−1−iv})^{−1/2}, 0 < v ≤ H, X := v log ξ.
  eqintero (Lemma 12, r = 2, applied to ξ and ξ^θ through lien-molli, valid for 0 < v ≤ 1/√log ξ ⊃ (0,H]):
  S_d(δ,v) = Π_{qd}(v)·s(δ,v) + O(τ_{1/2}(δ)∏_{p|d}(1+p^{−3/4})(log T)^{−2/3}),  s(δ,v) := −τ_{−1/2}(δ)·Δ(ξ,δ,θ,v)/((1−θ)2√π √v log ξ),
  Δ(ξ,δ,θ,v) = Δ_2(Xℓ1) − Δ_2(Xℓ2) = √X·[√ℓ1 E(Xℓ1) − √ℓ2 E(Xℓ2)]. Hence, with **S_X(u) := |√ℓ1 E(Xℓ1) − √ℓ2 E(Xℓ2)|**,
  |s(δ,v)| = τ_{1/2}(δ)·S_X(u)/((1−θ)2√π √(log ξ)).
4.3 v → 0 (S_d(δ,v) is a finite sum, continuous in v; the error is uniform): B_d(δ) = S_d(δ,0) = Π_{qd}(0) b(δ) + O(same),
  b(δ) = τ_{−1/2}(δ)(2/√π)[√log(ξ/δ) − √log⁺(ξ^θ/δ)]/((1−θ)log ξ) (= the r = 2 main term at v → 0: −v^{−1/2}Δ_2(vL)/(2√π) → (2/√π)√L).
  With **B̂(u) := (√ℓ1 − √ℓ2)/√(1−θ) ∈ [0,1]** (max at u = θ): |b(δ)| = τ_{1/2}(δ)·2B̂(u)/(√π √(1−θ) √(log ξ)) ≤ τ_{1/2}(δ) **C5'** /√log T, **C5' = 2/√(πϰ(1−θ))**
  (Dousselin's C5 = (e^ϱ+e^{ϱθ})Γ(1/4)/(2(1−θ)√(πϰϱ)Γ(3/4)) ≈ 2.66·C5' at θ → 0).
  Numerical check of these normalisations (q = 1, d = δ = 1, ξ = 10⁶, θ = 0.1): sums/main terms = 1.042 (B), 1.040, 1.037, 1.027 (S at X = 0.25, 0.5, 1); .scratch mtcheck.py.
4.4 Imaginary part (replaces coro-molli2's constant C4). Im(δ^{−iv}S_d(δ,v)) = −Σ_n M(nδ)τ_{−1/2}(nδ)sin(v log nδ)/n. Write δ^{−iv}S_d(δ,v)/Π_{qd}(v) = ∫_0^v ∂_u(·)du + real.
  ∂_u of it = −iδ^{−iu}Σ_n M τ log(nδ)n^{−1−iu}/Π_{qd}(u) − (·)Π'/Π; the first sum has the asymptotic of coro-molli2's proof (l.1146–1150; error O(∏(1+p^{−3/4})(log T)^{1/3}) uniformly
  in 0 < u ≤ v), whose main term is i∂_u of the r = 2 main term (both are the Hankel integrals ∫ y^s s^{−r}√(s+iu)ds; L·M_2 − M_3 = i∂_u M_2 by integrating by parts in s);
  |Π'/Π| ≤ ½Σ_{p|qd} log p/(p−1) ≤ ½(log log d + O_q(1)) ≪ log log T (worst case: d a primorial). Explicitly: the log(nδ)-weighted sum has error O(∏_{p|d}(1+p^{−3/4})·log T·√(loglog T))
  (l.1148, from Lemma 12 with r = 2 times log ξ and r = 3), against a main term of size (log ξ)^{3/2}; with F_r := (r−1)!Π·H_r, H_r := (1/2πi)∫y^s s^{−r}√(s+iu)ds (y = ξ/δ, L = log y),
  integration by parts in s gives iH_2' = LH_2 − 2H_3, so the main term of ∂_u(exact r = 2 sum) = −i(L F_2 − F_3) equals ∂_u F_2 up to the Π'(u) term. So, uniformly, Im(δ^{−iv}S_d(δ,v)/Π_{qd}(v))/v = Im(δ^{−iv}s(δ,v))/v + o(τ√log T·∏(1+p^{−3/4})²), and exactly
  Im(δ^{−iv}s(δ,v))/v = −τ_{−1/2}(δ)√(log ξ)·F4_X(u)/((1−θ)2√π), **F4_X(u) := Im(e^{−iXu}[√ℓ1E(Xℓ1) − √ℓ2E(Xℓ2)])/X
  = cos(Xu)(ℓ1^{3/2}IE(Xℓ1) − ℓ2^{3/2}IE(Xℓ2)) − u·sinc(Xu)·Re[√ℓ1E(Xℓ1) − √ℓ2E(Xℓ2)]**, IE(Y) := Im E(Y)/Y = Σ_k (−1)^{k+1}Y^{2k}/((2k+1)!((2k+1)²−1/4)).
  F4_0(u) = 4u(√ℓ1−√ℓ2) − (4/3)(ℓ1^{3/2}−ℓ2^{3/2}), sup ≈ 4/3 (Dousselin's C4·2√π(1−θ) ≈ 5.2).

## 5. Two lemmas on the d-sums
**Lemma Φ.** For f, g ≥ 0 on [0,1], d | δδ' ⇒ u_δ + u_{δ'} ≥ t_d := log d/log ξ ∈ [0,2]; with Φ_{f,g}(t) := sup{f(u)g(u') : u, u' ∈ [0,1], u+u' ≥ t} (non-increasing),
  Σ_{d|δδ'|d^∞} τ_{1/2}(δ)τ_{1/2}(δ') f(u_δ)g(u_{δ'})/(δδ') ≤ Φ_{f,g}(t_d)·Σ τ_{1/2}(δ)τ_{1/2}(δ')/(δδ') = Φ_{f,g}(t_d)/φ(d)  (last identity: Dousselin l.1442–1446).
  Computed with u-boxes (40 on [0,θ], 600 on [θ,1]): Φ(t_m) ≤ max{f_i g_j : u_i^+ + u_j^+ ≥ t_m}, step function on 200 cells of [0,2].
**Lemma W.** For F ≥ 0 of bounded variation on [0,2]: Σ_{d≤ξ², (d,q)=1} (d²/φ(d)³) F(t_d) ≤ (φ(q)/q) P1 log ξ ∫_0^2 F(t)dt + O(sup F + Var F).
  Proof: d²/φ(d)³ = Σ_{kn=d} μ²(k)ψ(k)/n with ψ ≥ 0, Σ_k μ²ψ(k) = P1 (Dousselin l.1453–1460); for each k, partial summation against Σ_{n≤z}χ0(n)/n = (φ(q)/q)log z + O(1)
  gives (φ(q)/q)∫_{log k}^{log ξ²} F(s/log ξ)ds + O(sup F + Var F) ≤ (φ(q)/q)log ξ∫_0^2 F + O(·). (F ≡ 1 recovers eqprodeul.)
  The error is additive, O(1) for every F used here (each is a product of ≤ 200-step upper step functions bounded by 4.4² and of |sin y|/X ≤ 1/(2ϰ)+2, so sup F + Var F = O(1));
  every use divides the sum by log ξ (or log T), so the error enters the constants as O(1/log T) = o(1), whatever the size of ∫F. The sums run over (d,q) = 1 because A_d = A'_d = 0 otherwise (l.1428, l.1567).
  P1 = ∏_p(1 + (3p²−3p+1)/((p−1)³p)) ∈ [11.54290, 11.54297] (primes ≤ 10⁶ in Arb, tail ≤ exp(3.00001/(10⁶−1))).

## 6. Mollified second moment C3 (lemme-mollifer, l.1406–1501), with C5' and Lemma Φ, W
|A_d| ≤ Σ|B_d(δ1)||B_d(δ2)|/(δ1δ2) ≲ (dq/(φ(d)φ(q)))·4Φ_BB(t_d)/(π(1−θ)log ξ·φ(d)), Φ_BB := Φ_{B̂,B̂}.
Σ1 = Σφ(d)A_d² ≲ (q/φ(q))·16P1/(π²(1−θ)²log ξ)·∫_0^2Φ_BB².  Σ2'' ≤ Σ_d A_d²·2φ(d)log d (l.1490–1497) ≲ (q/φ(q))·32P1/(π²(1−θ)²)∫_0^2 tΦ_BB².
Σ2' = 2Σφ(d)A_dÃ_d; Ã_d = Σ_kΛ(k)A_{[d,k]}, |A_{[d,k]}| ≤ (same bound)·Φ_BB(t_d)/φ([d,k]) (Φ non-increasing, [d,k] ≥ d), and Dousselin's k-sum ≤ log d + 2ϰ log T + o(log T)
(l.1476–1484, before his step log d ≤ 2ϰ log T): Σ2' ≲ (q/φ(q))·32P1/(π²(1−θ)²)∫_0^2 (t+2)Φ_BB².
∫|Lη²|² = (φ(q)/q)[Σ1 T log T(1+o(1)) − Σ2 T] (eqlemme4) ≤ (φ/q)T[Σ1 log T + |Σ2'| + 2Σ2''], so
  **C3 = (16P1/(π²(1−θ)²))·[ (1/ϰ)∫_0^2 Φ_BB² + ∫_0^2 (6t+4)Φ_BB² ]**   (code: C5'⁴P1(ϰ i0 + ϰ² i1)).
All bounds in this section are absolute: |Ã_d| ≤ Σ_kΛ(k)|A_{[d,k]}| (triangle inequality), and −Σ2 is bounded by |Σ2'| + 2Σ2'' (Σ2'' ≥ 0), never by a signed upper bound.
Remark (Dousselin's Σ2'): his displayed steps give |Σ2'| ≤ 2·Σφ(d)·[C5²dq/(log Tφ²φ)]·[4ϰC5²dq/(φ²φ)] ≲ 16ϰ²C5⁴(q/φ)P1, but l.1487 prints 8ϰ² (a factor 2 is lost between l.1485 and
l.1487; the printed bound is unsupported, not necessarily false); with 16ϰ², C3 = (1/(8ϰ)+2)(…)⁴P1 instead of (1/(8ϰ)+3/2)(…)⁴P1. Effect on his table (cert.log):
N = 1: 5.4529e-8 → 5.4277e-8 (printed 5.45e-8); N ≥ 2 unchanged at 3 digits.

## 7. ∫|I|², small v: 0 < v ≤ 1/log ξ (replaces ineImaPart–petitv, l.1590–1640)
Exact identity (φ_{iv}(d) = d^{1+iv}∏_{p|d}(1−p^{−1−iv}), Π_{qd}(v)² = ∏_{p|qd}(1−p^{−1−iv})^{−1}, (d,q)=1): with Ã'_d := A'_d/(Π_{qd}(v)Π_{qd}(0)),
  𝒯^{iv/2}φ_{iv}(d)A'_d² = e^{ivΛ_d}·(q/φ(q))²(d²/φ(d))·ε_q(v)·Ã'_d²,  Λ_d := log(√𝒯 d),  ε_q(v) := (φ(q)/q)∏_{p|q}(1−p^{−1−iv})^{−1} = 1 + O_q(v)
  (one q/φ(q) from ∏_{p|q}(1−p^{−1−iv})^{−1}, one from Π_{qd}(0)²; Lemma W then returns φ(q)/q, so one q/φ(q) survives, cancelled by 4T(φ(q)/q)).
So the local factors at p | d cancel (Dousselin bounds them separately: |φ_{iv}(d)| ≤ σ(d) and (d/φ(d))² from |A'_d|², giving P2 ≈ 79.9; here only d/φ(d) remains, giving P1).
|Im(ε e^{iα}w²)| ≤ |sin α||w|² + 2|w||Im w| + |ε−1||w|², and Ã'_d = Σ s̃(δ,v)b̃(δ')/(δ^{1+iv}δ') with s̃ = S_d/Π_{qd}(v), b̃ = B_d/Π_{qd}(0) real; Im Ã'_d = Σ b̃(δ')Im(δ^{−iv}s̃)/(δδ').
By §4 and Lemma Φ: |Ã'_d| ≲ Φ6(t_d)/(π(1−θ)^{3/2} log ξ·φ(d)), |Im Ã'_d|/v ≲ Φ4(t_d)/(π(1−θ)^{3/2}φ(d)), Φ6 := Φ_{S_X,B̂}, Φ4 := Φ_{|F4_X|,B̂}.
Lemma W with F(t) = Φ6²|sin(vΛ_d)|/v, vΛ_d = X/(2ϰ) + tX + O(v) (1/v = log ξ/X), and F = 2Φ6Φ4:
  **Im(K(v)/v) ≲ (q/φ(q))·P1/(π²(1−θ)³)·W(X),  W(X) := (1/X)∫_0^2 Φ6(t)²|sin(X/(2ϰ) + tX)|dt + 2∫_0^2 Φ6Φ4 dt**
(|ε−1| term: O(v)·Σ(d²/φ³)Φ6²/log²ξ = O(1/log T)). Then 4T(φ(q)/q)∫_0^{1/log ξ}(H−v)Im(K/v)dv ≲ (T/log²T)·(4P1/(π²ϰ(1−θ)³))∫_0^1 (A − X/ϰ)W(X)dX.
Implementation: 128 X-boxes; on each, S_X, F4_X on 640 u-boxes (Arb balls), |sin y|/X ≤ min(1/(2ϰ) + t, max_{y-range}|sin y|/X_lo); upper Riemann sums.
Error terms of §4 (factors ∏(1+p^{−3/4})², (log T)^{−1/6} relative) sum to o(1) as in Dousselin's eq-Im_Kv/v (l.1617–1630).

## 8. ∫|I|², large v: 1/log ξ < v ≤ H — integration by parts (replaces gainSel, l.1679–1700; kills the A log A term)
Dousselin (intKV, lastone): this part ≤ (1/(4(1−θ)²π log²ξ))Σ_d d(q/φ(q))Σ_{δ1,δ1',δ2,δ2'} τ(δ1)τ(δ2)/(δ1δ1'δ2δ2')·|B_d(δ1')B_d(δ2')|·|J|,
J = ∫_a^H (H−v)e^{ivΛ}D1(v)D2(v)dv, a = 1/log ξ, Λ = log(√𝒯 d/(δ1δ2)) ≥ (1/2)log T + (t_d − 2)log ξ + O(1) (δ_i ≤ ξ), so log ξ/Λ ≤ ϰ r(t_d)(1+o(1)),
r(t) := 1/(1/2 + (t−2)ϰ) (r(0) = 1/(1/2−2ϰ); r is decreasing, so r(t_m) bounds it on each t-cell); D_i = Δ(ξ,δ_i,θ,v)/v = log ξ·d_i(x), x = v log ξ ∈ [1, Aϰ],
d(x;u) = [Δ_2(xℓ1) − Δ_2(xℓ2)]/x, d'(x;u) = [ψ(xℓ1) − ψ(xℓ2)]/x² (since (Xφ − Δ_2)' = Xφ' ).
He bounds |D_i| by |c1|L_i + (32/3)/v without using e^{ivΛ}; the cross term gives (128/3)√π(1−θ)Aϰ log(Aϰ), the A log A in C1.
Integrate by parts with P := (H−v)D1D2 (P(H) = 0): J = −P(a)e^{iaΛ}/(iΛ) + (1/(iΛ))∫_a^H D1D2 e^{ivΛ} − (1/(iΛ))∫_a^H (H−v)(D1D2)'e^{ivΛ}, and once more on the middle integral:
  |J| ≤ (log ξ/Λ)[(Aϰ−1)S1(u1)S1(u2) + Aϰ(Dm(u1)Dp(u2) + Dm(u2)Dp(u1))] + (log ξ/Λ)²[S1S1 + Dm Dm + Dm1Dp2 + Dm2Dp1],
S1(u) := |d(1;u)|, Dm(u) := sup_{x≥1}|d(x;u)|, Dp(u) := ∫_1^∞|d'(x;u)|dx. Bounds used: u ≥ θ: Dm ≤ √(1−u)·max(sup_{[0,1]}|E|, sup_{Y≥1}|Δ_2(Y)|/Y) = 4.348√(1−u),
Dp ≤ sup_ℓ ℓ∫_ℓ^∞|ψ(Y)|Y^{−2}dY (ℓ ∈ the box of 1−u); u < θ: Dm ≤ D_B(θ) (x-boxes on [1,10], tail |c1|(1−θ) + 2R0max/10), Dp ≤ p_B(θ) (direct, x ≤ 200, tail 2Ψmax/200).
S1, Dm, Dp depend on δ_i only through u_i ∈ [0,1] and are bounded by the certified box values (uniform in d, δ_i); the O(1) in Λ is (1/2)log(q/2π), uniform.
|B_d(δ')| ≤ √(dq/(φ(d)φ(q)))·τ(δ')·2B̂(u')/(√π√(1−θ)√log ξ) (§4.3); Lemma Φ (f = S1, Dm, Dp; g = B̂), Lemma W with log ξ/Λ ≤ ϰr(t_d) inside the t-integral:
  this part ≲ (T/log²T)·(8P1/(π²ϰ(1−θ)³))·[A·ϰ(aS1 + aX) − aS1 + ϰ·a2],  aS1 = ½∫_0^2 r Φ_{S1,B̂}², aX = ½∫ r·2Φ_{Dm,B̂}Φ_{Dp,B̂}, a2 = ½∫ r²(Φ_{S1}² + Φ_{Dm}² + 2Φ_{Dm}Φ_{Dp})
(−aS1 is legitimate: it comes from the factor (Aϰ − 1) > 0 multiplying an upper bound).

## 9. Assembly
  **C1(A) ≤ pref·[(La + Sa)A + (Lb + Sb)],  pref = 8P1/(π²ϰ(1−θ)³) (P1 → P1eff by §10),  La = ϰ(aS1+aX),  Lb = −aS1 + ϰ·a2,
  Sa = ∫_0^1 W(X)dX/2,  Sb = −∫_0^1 X W(X)dX/(2ϰ)**   (code: w = W/2; pref·(A − X/ϰ)w = (4P1/(π²ϰ(1−θ)³))(A − X/ϰ)W).
Then §2 with C2 from §3, §6. Rigour: every constant is an Arb ball or an upward-rounded float (products ×(1+10^{−12}), sums ×(1+10^{−10})); κ_F is the lower end of the ball;
the free parameters (θ, ϰ = 0.1249999 < 1/8, A) only need to be admissible, they are chosen by a float search. Large N: Kb ≤ 0 (asserted), so C1 ≤ Ka·A and with A = αN,
N·κ_F ≥ 2π[1/α − (√(Kaα) + √(C2/N))²/α³], increasing in N: κ_F ≥ c(N1)/N for all N ≥ N1; c(∞) = π/(2Ka) at α = 2Ka.

## 10. Signed local factors at the small primes (removes the Euler product P1)
All main-term δ-sums of §6–§8 have the form Σ_d := Σ_{d|δδ'|d^∞} τ_{−1/2}(δ)τ_{−1/2}(δ') f(u_δ) g(u_{δ'}) δ^{−1−iv}δ'^{−1}, f, g continuous on [0,1] (0 at u = 1):
the main terms of S_d(δ,v)/Π_{qd}(v), B_d(δ)/Π_{qd}(0) and Im(δ^{−iv}S_d(δ,v)/Π_{qd}(v))/v are τ_{−1/2}(δ) times functions of u_δ only (§4.2–4.4; the last has no phase left,
so its local sums below are taken at v = 0). §5–§8 bounded Σ_d by sup|f|sup|g|·Σ τ_{1/2}τ_{1/2}/(δδ') = …/φ(d) (an upper bound for the absolute sum, which is, e.g., 7/2 − 2√2 at d = 2).
Fix a finite set P_0 of primes (all p ≤ 10⁶) and split d = d_0d_1, δ = δ_0δ_1, δ' = δ_0'δ_1' (index 0: supported on P_0; index 1: coprime to P_0). τ_{−1/2} is multiplicative and
d | δδ' | d^∞ splits prime by prime, so Σ_d = Σ_{δ_0,δ_0'} τ(δ_0)τ(δ_0')δ_0^{−1−iv}δ_0'^{−1} · [Σ_{δ_1,δ_1'} τ(δ_1)τ(δ_1') f(u_{δ_1} + ε_0) g(u_{δ_1'} + ε_0')/(δ_1^{1+iv}δ_1')], ε_0 = log δ_0/log ξ.
For fixed d_0 ≤ D_0 and K: if δ_0, δ_0' ≤ K then ε_0, ε_0' ≤ log K/log ξ → 0 and, by uniform continuity, the bracket is I_1 + o(1/φ(d_1)) (I_1 := bracket at ε = 0); the terms with
max(δ_0, δ_0') > K are ≤ η_K(d_0) sup|f||g|/φ(d_1), η_K(d_0) → 0 (K → ∞). With σ_p(k,v) := Σ_{a+b≥k} τ_{−1/2}(p^a)τ_{−1/2}(p^b)p^{−a(1+iv)−b} and (Σ_a τ_{−1/2}(p^a)x^a)² = 1 − x:
  σ_p(1,0) = −1/p, σ_p(k,0) = 0 (k ≥ 2), σ_p(k,v) → σ_p(k,0) uniformly for v ≤ H → 0.  Hence |Σ_d| ≤ (μ²(d_0)/d_0)|I_1| + o(1/(φ(d_0)φ(d_1))),
and Lemma Φ applies to I_1 with t_{d_1} = t_d + o(1). (Reviewer, logs/review.log: Claims 1–4 SOUND up to two exact-constant statements, which are upper bounds; for the Im part
the per-δ imaginary parts are taken first, so no Im σ/v term arises.)
d-sums: the weights d(d/φ(d))|local|² (§7, §8) and φ(d)(d/φ(d))²|local|² (§6) become (μ²(d_0)/φ(d_0))·d_1²/φ(d_1)³. Lemma W over d_1 coprime to qP (P := ∏_{p∈P_0}p):
Σ_{d_1}(d_1²/φ(d_1)³)F ≤ (φ(qP)/(qP))∏_{p∉P_0, p∤q}(1 + f1(p))·log ξ∫F + O_P(1), and Σ_{d_0}μ²(d_0)/φ(d_0) ≤ ∏_{p∈P_0, p∤q} p/(p−1); the product is
≤ (φ(q)/q)·∏_{p∉P_0}(1 + f1(p))·log ξ∫F: the Euler factor of every p ∈ P_0 is exactly 1 (it was 1 + f1(p), 4.5 at p = 2, 1.79 at p = 3). The part d_0 > D_0 (absolute bounds)
is ≤ (Σ_{d_0>D_0, P_0-smooth} d_0²/φ(d_0)³)·O(log ξ) → 0 (D_0 → ∞; the series converges because P_0 is finite).
Σ2' (§6): in Ã_d = Σ_kΛ(k)A_{[d,k]}, for k = p^m with p ∈ P_0 the signed main term of A_{[d,k]} vanishes unless m = 1, and those terms add ≤ (log d_0 + Σ_{p∈P_0} log p/(p−1))·(|A_d|-bound)
= O_{P_0}(1)·(|A_d|-bound) = o(log T)·(|A_d|-bound); primes p ∉ P_0 as in §6; so Ã_d keeps the factor (log d + 2ϰ log T + o(log T)) with signed weights.
Order of limits: T → ∞, then K, D_0 → ∞, with P_0 fixed. **Conclusion: P1 is replaced by P1eff := ∏_{p>10⁶}(1 + f1(p)) ≤ exp(3.00001/(10⁶ − 1)) = 1.000003001** in pref and C3.
Plausibility (not a proof): Σ1 ≲ (q/φ(q))·16P1eff∫Φ_BB²/(π²(1−θ)² log ξ) ≈ 1.13/log ξ (θ → 0), against the exact 0.889/log U (sharp cutoff) and 0.678/log U (power 1/4) computed
for ζ in arXiv 2609.15329 (Thm in §1 there) for mollifiers of the same kind (coefficients of ζ^{−1/2}); with P1 the bound was 15 times larger.

## 11. Constants (logs/cert.log; every number below is rg-checked against the log)
Reproduction: N=4, θ=0.0013, A=425715589.6389222: 3.682008e-09 (printed 3.68e-09); N=10: 1.468218e-09; P1 = [11.5429 +/- 6.25e-5], P2 = [79.886 +/- 5.21e-4];
with the Σ2' bound doubled: N=1 5.427728e-08.
Special functions: R0max 2.338139, DA 4.348005, Emax 4.334941, Psimax 2.382235, Ipsi[0,1] 2.119367.
Ablation (N = 4): §1–2 on Dousselin's C1, C2: 4.5418e-08; + §4–§8 with sup weights, C2 = 4(√C3+1)²/(θϰ)² and constant log ξ/Λ ≤ ϰ/(1/2−2ϰ): 4.1504e-05 (Ka=6442.16); + Lemma Φ, first moment (§3), r(t) (§8):
9.7320e-05 (θ = 0.1, Ka=2398.069, C3=282.0849); + §10 (P1eff = 1.000003000): 6.8158e-04.
Final constants, ϰ = 0.1249999: θ = 0.150: C5'=3.461711 Ka=251.9215 Kb=-471.7217 La=5.87599 Sa=17.98234 C3=29.80734 C2=3.277641e+05 pref=10.55906,
aS1r 21.27826, aXr 25.72968, D_B 4.0603, p_B 1.93901, i0 0.87381, i1 6.29385. θ = 0.100: Ka=207.7522, C2=6.000112e+05. θ = 0.030: Ka=158.7670, limit pi/(2 Ka) = 9.8937e-03.
Results κ_F ≥ (θ, A): N=4 6.8158e-04 (0.150, 5265.25); N=5 5.7578e-04 (0.150, 6196.71); N=6 5.0081e-04; N=7 4.4455e-04; N=8 4.0060e-04; N=9 3.6520e-04;
N=10 3.3601e-04 (0.150, 10428.5); N=100 5.1386e-05; N=1000 6.8655e-06; N=1 1.9004e-03; N=2 1.1472e-03; N=3 8.4751e-04.
Large N (c(N1), §9): all N ≥ 11: 3.4254e-03/N (θ = 0.150); N ≥ 100: 5.1183e-03/N (θ = 0.100); N ≥ 1000: 6.8429e-03/N (θ = 0.050); N ≥ 10⁶: 9.6960e-03/N (θ = 0.030).
