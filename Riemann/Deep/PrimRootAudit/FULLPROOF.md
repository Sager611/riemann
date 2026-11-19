# PrimRootAudit — FULLPROOF: chain of discoveries/prime-primitive-root.md with source lines
D: = Riemann/Deep/PrimRootNum/DERIVATION.md, C: = cert.py, T: = tail.py, S: = sieve.py, EPD: = ExplProg/DERIVATION.md. Checked = re-derived by hand (ppraudit).

## 1. Detection and characters (D:27–47) — checked
- L1.1 (D:29): m = r^k ∈ [x e⁻⁷, x] generating G ⇒ r generates G (⟨r^k⟩ ⊆ ⟨r⟩), r ≤ x. w ≥ 0 (triangle × t^{−B}). OK.
- (1.1) (D:36): G cyclic, G/G^d cyclic of order d ⇒ 1_{G^d} = d⁻¹Σ_{χ^d=1}χ on (m,p)=1; χ ≠ χ₀ kills p | m ⇒ A_d = A_1/d + r_d, |r_d| ≤ (1−1/d)R. OK.
- (1.2) (D:41): χ ≠ χ₀ mod prime p is primitive, conductor p ≥ 3; EPD Lemma EF (D:40 needs x t₀ = e^Y > 1), Lemma P (EPD:62, pairs β ↔ 1−β), Lemma Z (EPD:71: primitive χ of conductor 3 ≤ q* ≤ q, q ≥ 9 — here q = e^{L₉}, L₉ = max(log p, log 9) ✓), K1 (EPD:117, Y > 0). Trivial zeros ≤ W(0)/(1−e^{−2Y}) incl. even-χ zero at 0 (EPD:44–46). OK.
- (1.3) (D:42): ζ explicit formula: xW(1) − Σ_ρ x^ρW(ρ) − Σ_m x^{−2m}W(−2m) (last ≥ 0, subtracted); Ω₁ ≤ 0.0231856(x^B k₂ + x^{1/2}k₂^L) (EPD:96–99; C:47 z1 = SUM_INV_G2·k2/2, SUM_INV_G2 = 0.046371 = Σ over all zeros ✓). Removed m = p^k: ≤ ‖w‖_∞ log p (1 + 7/log p). OK.
- Freezing k's at Y_m (D:97–99): K1's k₁, k₂ = n₂ + 2n₁/(eY) + (2/(eY))²n₀ decrease in Y, so the envelope G_Y ≤ G_{Y_m} POINTWISE; Lemma Z is applied to G_{Y_m} (no monotonicity of Z_G in k needed — ρ₀ < 0 near 0 would break that). OK.

## 2. Sieve lemmas (D:49–78) — checked
- S1 (D:54): unique decomposition d = e f (e shortest non-D prefix, f | n_e) ⇒ σ_n = [n=1] − Σ_{e bdry, e|n} μ(e)[n_e = 1]; μ(e) = +1 (r even) ⇒ σ_n ≤ [n=1]. Rosser rule (D:60): conditions only at even positions ⇒ boundary even. Restriction to d | P_ω (P⁺(d) ≤ p_ω) is again a DFS set (prefixes keep P⁺). OK.
- S2 (D:64–73): Σλ_d g(d) = Σ_m σ_m Π_{p|m}g Π_{p∤m}(1−g) (inclusion probability); h = g/(1−g) = 1/(l_j −1) ≤ 1/(p_j −1) as l_j ≥ p_j; σ_m ≤ 0 (m ≥ 2) ⇒ Σλ_d/d' ≥ V_Υ T/V_ω. With X ≥ 0: S ≥ X(V_Υ/V_ω)T − R Σ_{d≠1}|λ_d|. OK.
- (2.2) (D:76): V_Υ ≥ V_ω, T ≥ 0, A_1 ≥ 0 ⇒ S ≥ A_1T − R(C−1). ν(m) = Π{l: m ∈ G^l}, d'|ν(m) ⇔ m ∈ G^{d'} (cyclic). OK.

## 4. ω ≤ 10⁴ (D:91–113; C:36–53, 135–165) — checked
- x = (κL)^{1/b} ⇒ x/(x^B L) = κ; (4.1)/(4.2) = (1.3)/(1.2) divided by x^B L; (4.1) subtracted terms non-increasing in L with frozen k's; Ro bounded by sup-majorant R̂_m (C:30–34: else-branch valid for any Lm). (4.4) ⇒ for all L ≥ L_m: T·Xo(L) − (C−1)Ro(L) ≥ T·Xo(L_m) − (C−1)R̂_m > 0 (T > 0, C ≥ 1). OK.
- L_m = log(P_ω + 1) as an Arb ball of the exact integer (C:154); p − 1 ≥ Υ ≥ P_ω. OK. Admissibility Y_m > 1 (C:41).
- Thm B (D:110): c³L_ω² ≥ κ ⇒ κL ≤ c³L³ (L ≥ L_ω). C:157–158 certifies c_w³ ≥ κ/Lm². OK.
- Thm B′ (C:375–386): κ RE-certified at Lm = max(L_ω, log 10⁷) and c_w³Lm² ≥ κ certified; primes p > 10⁷ have log p > log 10⁷. OK.

## 5. Tail ω(p − 1) > 10⁴ (D:114–193; C:167–265; T:1–264) — checked
- T1 (D:123): Buchstab at w₀: 1_{(n,P⁺(z))=1} = 1_{(n,P(w₀))=1} − Σ_{w₀≤q≤z, q|n} 1_{(n,P(q))=1}, each 1_{(n,P(q))=1} ≤ (Σ_{d|(n,P(q))}ρ_d)² ⇒ σ_n ≤ [n=1] on n | P_ω; λ_1 = 1. OK.
- T2 (D:127): Selberg's minimum Σρρ/[d₁,d₂] = 1/G (h(d) = 1/φ(d)); Σ|ρ_d| = G⁻¹Σ_{k<ξ}μ²(k)σ(k)/φ(k) ≤ c_σ ξ/G (σ(k)/φ(k) = Σ_{e|k}Π_{p|e}2/(p−1)). OK.
- (5.1) (D:134): telescoping Σ_{q≥w₀} V(q)/q = V(w₀) − V⁺(z) (no primes in [50,53)); Φ_q ≥ 1 as G ≤ 1/V(q). (5.2): |λ_e| ≤ 1_{e|P(w₀)} + Σ_q Σ|ρρ| ⇒ C − 1 ≤ 2¹⁵ + c_σ²DΣ_q 1/(qG_q²), D = z^s. OK.
- T3 (D:137): (a) 1/φ(d) = Σ_{rad m = d}1/m ⇒ G ≥ Σ_{m<min(ξ,q)}1/m ≥ log min(ξ,q); (b) m = l·k, Σ_{k<y}1/k ≤ 1 + log y; (c) Rankin, Π(1 + p^α/(p−1)) = V⁻¹Π(1 + (p^α−1)/p). OK.
- (5.3) (D:146): Stieltjes by parts, f' = −1/x; ∫_{x₀}^{ξ}(log(log x/log x₀) + η₀)dx/x = log x₀(a₀log a₀ − a₀ + 1) + η₀ log(ξ/x₀). OK. (5.4) (D:148): division by log q with η₀ ≤ η_B, 1/log q ≤ 1/log Q₂ (q > Q₂). OK.
- T:G_lo (T:105–143): (a) min of lower endpoints; (b) exact sum over primes ≤ 283 in [q, ξ) (undecided l included = larger subtraction), Mertens on [max(q,288), ξ) only if ξ may exceed 293 (no primes in [284, 292]); (c) r = e^{−α lx}Π_α(q), Π_α and V over p < q (dict stores the product BEFORE p). OK.
- Part B ψ (D:153–157; T:193–214): three candidates each valid for prime q > Q₂ with a = a_q; each non-increasing in a (h_Q concave, h_Q'(a_cap) > 0 — recomputed below), so Ψ(t) = ψ((s/t − 1)/2) non-decreasing. R_B(c) (D:158–161; T:80–86): primes ≤ 283 by monotonicity of p^{c/log q} in q; Abel on [288, q) with |E| ≤ 1/(2log²x) (x ≥ 286): main ≤ Σ c^k/(k·k!), E-part ≤ (e^c − 1)(1/(2log²Q₂) + 1/(2log²288)). ein_plus tail ≤ 2c⁸⁰/80! (c ≤ 6). OK.
- (5.5) (D:163–167): re-derived: ∫F dS₀ = u⁻¹∫Ψt⁻²; |F(z)E(z)| ≤ Ψ(1)/(2u³), |F(Q₂)E(Q₂)| ≤ K₃/(2u³), ∫|E||dF| ≤ (Ψ(1) + 4J)/(2u³) (dΨ ≥ 0, ∫dΨ/t³ ≤ Ψ(1) + 3∫Ψ/t⁴). Total (2Ψ(1) + 4J + K₃)/(2u³). OK.
- Integrals (T:219–253): upper sums use Ψ(t₁) (right end, Ψ non-decreasing) with exact ∫t⁻², ∫t⁻⁴ per cell; K from Ψ(t₁)/t₀³; (0, 1/50] by one Rankin candidate with r(1/50) < 1/2, closed forms ∫₀^{t_m}e^{−k/t}t⁻²dt = e^{−kW}/k, ∫e^{−k/t}t⁻⁴ = e^{−kW}(W²/k + 2W/k² + 2/k³), sup w³e^{−kw} = W³e^{−kW} (W ≥ 3/k checked). OK.
- Cells (C:192–208): u from u₁ = log p_{10001} (ball) in steps 0.1/0.25/0.5 while u < 60 — contiguous, first cell starts AT u₁, last ends ≥ 60; Part A at left end (G ↑ in ξ ↑ in u), V⁺ lower bound e^{−γ}(1 − 1/(2u_l²))/u_r, Part B integrals over [t_i, 1] with t_i ≤ log Q₂/u_r, error terms at u_l. u ≥ 60 (C:210–223): r_q(u) = r_q(60)e^{−κ(u−60)}, κ = αs/2 = 0.138, ue^{−κ(u−60)}/(1 − 1/(2u²)) decreasing (60κ ≥ 1, certified); B_rel decreasing in u. No gap. OK.
- Cost (T:158–175, 260–262): K_A at u₁ (smallest ξ); K_B: G ≥ log q·min(1, a_q) ≥ 0.65 log q; Σ_{q>Q} 1/(q log²q) ≤ (1/(2log²Q))(1 + 3/(2log²Q)) by Abel (mine) ≤ stated (1 + 2/log²Q). OK.
- (5.6) (D:179–192; C:234–257): θ(z) ≤ log(p−1) < L, R5 ⇒ u ≤ ℓ(L); T ≥ τV⁺ ≥ τ'/ℓ(L); x/(x^B L) = c_T^β L^{s+δ}; division by L^s gives (5.6) exactly as coded. Extension: L^δ/ℓ ↑ when δℓ ≥ 1 (δℓ(L₁) = 1.165); sub/(L^sℓ) ↓; Ro ≤ R̂₁ (sup-majorant at L₁ ≥ log 9); L^{−s} ↓. OK.

## 6–7. Small primes, Theorem D, B′, Corollary E (D:195–228; C:268–352, 362–392) — checked
- Small (C:320–352): float scan is information only; proof = max g* (p ≤ 10⁷) = 211 < (0.814 log 5)^{3/b} = 601.44 (Arb) and monotonicity in p. p = 3 excluded from the 0.814 form (B with 1.4312 and D cover it).
- D (C:268–281): need_w = (κ_w/L_w^{β−1})^{1/β}, certified c^βL_w^{β−1} ≥ κ_w for all 10⁴ w and c ≥ c_T. OK. Exponents 23.70 ≥ 23.6933, 26.86 ≥ 26.8524 (bases > 1 for p ≥ 5 resp. p ≥ 3). OK.
- E (C:283–318): per w either table κ already gives (1/b)log(κL) + slack < fac·L at L_w with fac·b·L_w > 1 (difference increasing), or w ≤ 120 scanned on L_w + k/4 with κ RE-certified at each L_m; Λ = max upper endpoint; tail via L = 10⁵ < L₁ and L/2 − β/b·log(c_T L) increasing for L > 2β/b. slack 10⁻³: x e^{0.001} < √p and 0.001x > 2 (x > e⁸). OK.


## 0. Inputs as used (checked against sources)
- H (D:9–11): B = 0.8733819666132076 > B₂ + 10⁻¹⁷ certified (cert_output [core]); Bebel's half-plane audited in HPAudit (every q, χ, so ζ).
- EP lemmas (D:12–18 → EPD:40–46 EF, EPD:62–69 P, EPD:71–94 Z, EPD:96–99 Zζ, EPD:117–125 K1): hypotheses met — Y = log x − 7 > 1 > 0
  (C:41), conductor p ≤ e^{L₉} with L₉ ≥ log 9 (C:30–34), primitive χ (p prime). Kernel constants and Lemma Z coefficients recomputed by
  quadrature (indep_kappa.py, zcmp.py): equal to EP's closed forms to 10–14 digits.
- RS62 (D:19–22; T:9–13): Thm 5 Cor. |Σ1/p − log log x − M| ≤ 1/(2log²x), x ≥ 286; Thm 7 lower Π bound x ≥ 285, upper x > 1; Thm 4
  θ(x) > x(1 − 1/log x), x ≥ 41 — ranges from the quotation in arXiv 2608.01498 (TeX l.686, 806, 826, 1197–1200). Every use: x ≥ 288.
- Reviewer subagent (S1, S2, (2.2), T1, application to primitive roots): VERDICT CORRECT.

