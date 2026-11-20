# B3Cert — FULLPROOF (b3cert, Orca, 20 Nov 2026)
Sources: SL = Santibáñez-Leal, Zenodo 23134787 v0.03, TeX = main.tex of the source zip 23135354 ("l." = its lines); E1 = Riemann/Deep/SantiAudit/E1.md;
E1A = Riemann/Deep/E1Audit/FULLPROOF.md; CA = Riemann/Deep/CeilAudit/REPORT.md (B3). Notation as E1 §2: T large, L = log T, H = T^θ, K := T/H, τ := T/H²,
shifts α, β with |α|, |β| ≤ C0/L (so |d|, |z| ≤ 1/8 for L ≥ 16C0; d = α+β, z = α−β), coprime p, q, pq ≤ y² = T^{2ν}, n ≥ 1,
X0 := 2√(2π)√(nT/(pq)), N0 := pqT/H², r := n/N0, c_r := X0/K = 2√(2π)·√r. X0 ≥ 2√(2π)T^{1/2−ν} → ∞.

## 1. The exact profile (SL l. 357–390, no Taylor step)
SL (10) is exact with kernel K_β(x) = (2x^{−β}/√π)∫e^{−y²+(1−2β+2iT)y/H}cos(2πxe^{2y/H})dy. Writing cos = ½Σ_ε e^{iε2πxe^{2y/H}}:
  K_β(x) = Σ_{ε=±1} e(εx)·A_ε(x),  A_ε(x) := x^{−β}π^{−1/2}∫_ℝ e^{−y²+(1−2β)y/H}·e^{2iTy/H}·e^{iελ(y)x}dy,  λ(y) := 2π(e^{2y/H} − 1)   (exact, every x > 0).
Everything below uses A_ε itself; SL's Taylor profiles A_{ε,j} are used nowhere (their integrated series need not converge to A_ε, §7 Repair). Cutoff: χ ∈ C_c^∞(0,∞) fixed, χ̃(u) := χ(u/2π),
supp χ̃ ⊂ [1/2, 2], χ̃ = 1 on [2/3, 3/2]. Retained kernel χ(x/T)Σ_ε e(εx)A_ε(x); the rest (1 − χ(x/T))K_β(x) is SL's off-band part (§4 below).
Dual profile (ε = −1; it replaces E1's truncated Taylor profiles, no expansion): A(ξ) := χ(ξ/T)A_{−1}(ξ) ∈ C_c^∞(0,∞) (differentiate under ∫dy). Hence E1 Steps 1–2 hold verbatim
(they need only A ∈ C_c^∞): 𝒯(n) = −π∫_0^∞A(ξ)ξ^{−z/2}B(X(ξ))dξ, X(ξ) = 4π√(nξ/(pq)), and with Step 2's Hankel form 𝒯 = 𝒯_+ + 𝒯_−,
  𝒯_σ(n) = κ_σ∫_0^∞A(ξ)ξ^{−z/2}X^{−1/2}a_σ(X,d)e^{iσX}dξ,  κ_σ = iσ√(π/2)e^{−iσπ/4} (|κ_σ| = √(π/2)),  X^k∂_X^k a_σ = O_k(1) for |d| ≤ 1/4, all X > 0 (E1A §2).

## 2. Lemma B (uniform exact-profile dual bound; Step 4″)
**Lemma B.** Constants C1, C_N depend only on χ, C0, N. Assume X0 ≥ 1.
 (i) (bulk) If c_r ≤ H/46, then |𝒯_σ(n)| ≤ C1·H·X0^{−1/2} (σ = ±1), so |𝒯(n)| ≤ 2C1·H·X0^{−1/2}.
 (ii) (tail) For every N ≥ 0 and all n: |𝒯(n)| ≤ C_N·T·X0^{−1/2}·c_r^{−N}.
No condition on τ = T/H² or on θ; only 1 ≤ H ≤ T, X0 ≥ 1 and the shift disk ((ii) uses T^{−1} ≤ H^{−1}).
**Proof of (i).** Step a (coordinates). Insert A into 𝒯_σ; the double integral converges absolutely (∫e^{−y²+(1+2|β|)|y|/H}dy < ∞, ξ in supp χ(·/T)),
so Fubini applies. Put ξ = Tu/(2π) and Λ := H(e^{2y/H} − 1) ∈ (−H, ∞), y = (H/2)log(1 + Λ/H), dy = dΛ/(2(1 + Λ/H)). Then
2Ty/H = Kℓ(Λ), ℓ(Λ) := H log(1 + Λ/H); λ(y)ξ = KΛu; σX(ξ) = σX0√u = Kσc_r√u. With ρ := Λ/c_r (so KΛu = X0ρu):
  𝒯_σ = κ_σπ^{−1/2}(T/2π)X0^{−1/2}·𝒥_σ,  𝒥_σ := ∫_{−H}^∞ f(Λ)e^{iKℓ(Λ)}J(Λ)dΛ,  J(Λ) := ∫α(u)e^{iX0(σ√u − ρu)}du,
  α(u) := χ̃(u)(Tu/2π)^{−β−z/2}u^{−1/4}a_σ(X0√u, d),   f(Λ) := ½e^{−y²−(1+2β)y/H}   (y = y(Λ)).
Bounds: supp α ⊂ [1/2, 2]; ‖α^{(k)}‖_∞ ≤ A_k (k ≤ 3) uniformly (|(Tu/2π)^{−β−z/2}| ≤ e^{2C0}, ∂_u = (X/2u)∂_X on a_σ). |f| ≤ ½e^{(1+2|β|)²/(4H²)} ≤ 1,
Var_ℝ f = ∫|∂_y f|dy ≤ 2 (c := (1+2β)/H, |c| ≤ 5/4: complete the square, −y² − cy = −(y + c/2)² + c²/4, so ∫|∂_y f|dy ≤ ½e^{(Re c)²/4}∫(2|t| + |Im c|)e^{−t²}dt ≤ ½e^{0.4}(2 + 0.45) < 2; reviewer repair). It suffices to show |𝒥_σ| ≤ C·K^{−1}: then |𝒯_σ| ≤ 2^{−1/2}(T/2π)X0^{−1/2}CK^{−1} = (C/(2π√2))HX0^{−1/2}.
Step b (non-stationary Λ). Let ρ_a := 1/(4√2), ρ_b := √2, I_σ := {Λ : σρ ∈ [ρ_a, ρ_b]}, N_σ := (−H, ∞) \ I_σ. For u ∈ [1/2, 2], ∂_u(σ√u − ρu) = σ/(2√u) − ρ
and σ/(2√u) ∈ σ[1/(2√2), 1/√2]; for σρ ∉ (ρ_a, ρ_b) one checks |σ/(2√u) − ρ| ≥ (1 + |ρ|)/8 (cases σρ ≤ 0, 0 < σρ ≤ ρ_a, σρ ≥ ρ_b). The u-phase has
|∂²_u| = 1/(4u^{3/2}) ≤ 1/√2, |∂³_u| = 3/(8u^{5/2}) ≤ 2.13, independent of ρ. Two integrations by parts (α compactly supported, no boundary terms):
  |J(Λ)| ≤ X0^{−2}∫|((α/φ')'/φ')'|du ≤ B2·X0^{−2}(1 + |ρ|)^{−2}  (φ = σ√u − ρu; B2 explicit from A_0..A_2),
  ∫_{N_σ}|f||J|dΛ ≤ B2X0^{−2}c_r∫_ℝ(1 + |ρ|)^{−2}dρ = 2B2·X0^{−2}c_r = 2B2·K^{−1}X0^{−1}   (X0 = Kc_r).
Step c (stationary range: complete the square). On I_σ, |ρ| ≥ ρ_a > 0. Substitute u = s², β̃(s) := 2sα(s²) (supp ⊂ [2^{−1/2}, 2^{1/2}], ‖β̃^{(k)}‖_1 ≤ B̃_k):
  σs − ρs² = −ρ(s − s0)² + 1/(4ρ), s0 := σ/(2ρ)  ⟹  J(Λ) = e^{iψ0(Λ)}·m(Λ),  ψ0 := X0/(4ρ) = Kc_r²/(4Λ),  m(Λ) := ∫β̃(s)e^{−iX0ρ(s−s0)²}ds.
Second-derivative test (Titchmarsh Lemma 4.4: |∫_a^b e^{iF}| ≤ 8/√r0 if |F''| ≥ r0 > 0 on [a, b]) plus one partial integration in s
(|∫_a^b g e^{iF}| ≤ (|g(b)| + ∫_a^b|g'|)·sup_x|∫_a^x e^{iF}|): the s-phase has |F''| = 2X0|ρ| ≥ 2X0ρ_a, so
  (c1) |m| ≤ 8(2X0ρ_a)^{−1/2}B̃_1 =: M0·X0^{−1/2}  (B̃_k := ‖β̃^{(k)}‖_1).
∂_ρ[ρ(s − s0)²] = (s − s0)² + (σ/ρ)(s − s0) and (s − s0)e^{−iX0ρ(s−s0)²} = (i/(2X0ρ))∂_s e^{−iX0ρ(s−s0)²}. One/two partial integrations in s:
  |∫β̃(s−s0)e^{..}| ≤ 8(2X0|ρ|)^{−3/2}B̃_2,  |∫β̃(s−s0)²e^{..}| ≤ 8(2X0|ρ|)^{−3/2}B̃_1 + 8(2X0|ρ|)^{−5/2}B̃_3, hence
  (c2) |∂_ρ m| ≤ X0[|∫β̃(s−s0)²e| + |∫β̃(s−s0)e|/|ρ|] ≤ 8(2ρ_a)^{−3/2}[B̃_1 + B̃_2/ρ_a + B̃_3/(2X0ρ_a)]·X0^{−1/2}, so Var_{I_σ}m = ∫|∂_ρm|dρ ≤ M1·X0^{−1/2}.
Step d (outer Λ-integral). On I_σ, ∫f e^{iKℓ}J dΛ = ∫(fm)e^{iψ}dΛ with ψ := Kℓ + ψ0,
  ψ''(Λ) = K[−H^{−1}(1 + Λ/H)^{−2} + c_r²/(2Λ³)],  |c_r²/(2Λ³)| = 1/(2|ρ|³c_r) ≥ 1/(2ρ_b³c_r) = 1/(4√2c_r), sign σ.
On I_σ, |Λ| ≤ ρ_bc_r ≤ √2H/46 < H/2, so H^{−1}(1 + Λ/H)^{−2} ≤ 4/H. σ = −1: both terms negative, |ψ''| ≥ K/(4√2c_r). σ = +1: ψ'' ≥ K[1/(4√2c_r) − 4/H]
≥ K/(8√2c_r) because 4/H ≤ 1/(8√2c_r) ⟺ c_r ≤ H/(32√2) = H/45.25 (hypothesis c_r ≤ H/46). So |ψ''| ≥ r0 := K/(8√2c_r),
one sign on the interval I_σ. Second-derivative test with amplitude fm (C¹ on I_σ):
  |∫_{I_σ}(fm)e^{iψ}| ≤ 8r0^{−1/2}(sup|fm| + Var(fm)) ≤ 8(8√2)^{1/2}(c_r/K)^{1/2}·X0^{−1/2}[M0(1 + 2) + M1]   (|f| ≤ 1, Var f ≤ 2)
                      = 8(8√2)^{1/2}[3M0 + M1]·K^{−1}   since (c_r/K)^{1/2}X0^{−1/2} = (c_r/K)^{1/2}(Kc_r)^{−1/2} = K^{−1}.
With Step b: |𝒥_σ| ≤ {8(8√2)^{1/2}(3M0 + M1) + 2B2}·K^{−1} (X0 ≥ 1), so C1 = (2π√2)^{−1}{8(8√2)^{1/2}(3M0 + M1) + 2B2}. ∎(i)
Remark. Nothing used τ: the y-chirp of the profile (curvature 4τ) never has to be resolved; the only large parameters are X0 (inner) and K/c_r (outer),
and the two gains X0^{−1/2} (inner) and (c_r/K)^{1/2} (outer) multiply to K^{−1}: E1's per-term size HX0^{−1/2} for every θ. The Gaussian factor f is only
used through sup|f| ≤ 1 and Var f ≤ 2 (so (i) is uniform in r; the e^{−πr/2} decay is not needed, (ii) supplies the tail).
**Proof of (ii)** (E1 Step 4 with the exact profile). ∂_ξ^k A_{−1}(ξ) = Σ_{i≤k}binom(k,i)∂^{k−i}(ξ^{−β})·π^{−1/2}∫(−iλ(y))^i e^{..}dy and, for every y,
|λ(y)| ≤ (4π|y|/H)e^{2|y|/H}, so |∂_ξ^kA(ξ)| ≤ C_k H^{−k} POINTWISE on supp χ(·/T) (|ξ^{−β}| ≤ e^{C0}, T^{−1} ≤ H^{−1}); no E_K(ξ) decay is claimed.
E1 Step 4 verbatim (𝒟G = (G/X')', each 𝒟 costs ≪ (T/X0)H^{−1}) with ∫_{supp}dξ ≤ 2T gives |𝒯| ≤ C_N T X0^{−1/2}(T/(HX0))^N = C_N T X0^{−1/2}c_r^{−N}. ∎(ii)
**Corollary B′ (Step 4″ in the brief's form).** Fix ν < 1/2, θ ∈ (0, 1), H = T^θ (or H ≥ T^θ, H ≤ T^{1−δ}). For all K′ ≥ 0, ε > 0, uniformly in n, p, q, shifts:
  |𝒯(n)| ≪_{K′,ε} H(pq/(nT))^{1/4}(1 + n/N0)^{−K′}T^ε,   and  Σ_{n≥1} d(n)n^{ε}|𝒯(n)| ≪_{ε,A} H(pq/T)^{1/4}N0^{3/4}T^{3ε} + T^{−A}  (any fixed A).
Proof. X0^{−1/2} = (2√(2π))^{−1/2}(pq/(nT))^{1/4}. If n ≤ N0T^{ε}: c_r ≤ 2√(2π)T^{ε/2} ≤ H/46 for T ≥ T0(ε, θ) (ε < θ), so (i): |𝒯| ≤ 2C1HX0^{−1/2} ≤
2C1·2^{K′}HX0^{−1/2}(1 + n/N0)^{−K′}T^{εK′}. If n > N0T^{ε}: (ii) with N ≥ 2K′ + 2(1 + 10 + K′)/ε: TX0^{−1/2}c_r^{−N} ≤ HX0^{−1/2}(N0/n)^{K′}·K·(8π)^{−N/2}T^{−ε(N/2−K′)}
≤ HX0^{−1/2}(1 + n/N0)^{−K′}2^{K′}T^{−10}. Rename εK′ → ε. Sum: first range ≪ H(pq/T)^{1/4}Σ_{n≤N0T^ε}d(n)n^{−1/4+ε} ≪ H(pq/T)^{1/4}(N0T^{ε})^{3/4+ε}L
(empty if N0T^ε < 1); second range (reviewer repairs; only N0 ≤ T² is used, from pq ≤ T, H ≥ 1): TX0^{−1/2}c_r^{−N} ≤ T(N0/n)^{N/2}, and with Y := N0T^ε,
  Σ_{n>Y}d(n)n^{ε}(N0/n)^{N/2} ≪ N0^{1+2ε}T^{ε(1+2ε−N/2)} ≤ T^{2+6ε−εN/2} (Y ≥ 1),  ≪ N0^{N/2} ≤ T^{−εN/2} (Y < 1)   (N/2 > 1 + 2ε);
  so the tail is ≪ T^{3+6ε−εN/2} ≤ T^{−A−2+6ε} ≪ T^{−A} for N ≥ 2(A + 5)/ε (ε ≤ 1/5; N chosen after A, ε). ∎
Review of Lemma B: three reviewer passes (§7): CORRECT for (i)–(ii) and the cancellation mechanism; the bookkeeping points they raised (Var f, H ≤ T, tail sum) are repaired above.


## 3. Step (2): residues = exact Mellin value of the exact profile (main term unchanged)
**Lemma R.** For |w − 1| ≤ 1/4, w := s − β, let M(s) := ∫_0^∞χ(x/T)A_{−1}(x)x^{s−1}dx. Then M(s) = (H/(2√π))(T/2π)^{w−1}(1 + E), |E| ≤ C_R·H/T, uniformly
(C_R depends on χ only); M is holomorphic in (s, β). Write M_lead(s) := (H/(2√π))(T/2π)^{s−β−1}, so M = M_lead(1 + E).
Direct comparison (two explicit integrals, no series): SL's closed-form leading Gaussian A_{−1,0}(x) = x^{−β}exp(−v² − i(1−2β)v/H + (1−2β)²/(4H²)), x = T/2π + Hv/2π
(SL l. 709–711), has M_0(s) := ∫χ(x/T)A_{−1,0}(x)x^{s−1}dx = M_lead(s)(1 + O((H/T)L^{O(1)})) (SL l. 716–717), hence M − M_0 = M_lead·O((H/T)L^{O(1)}).
Proof. Substitute v = e^{2y/H} − 1: A_{−1}(x) = x^{−β}P(x), P = π^{−1/2}Ĝ, Ĝ(x) = ∫G(v)e^{−2πixv}dv, G(v) := g0(v)(1 + v)^{iT}·1_{v>−1},
g0(v) := (H/2)(1 + v)^{−1}e^{−y²+(1−2β)y/H} = (H/2)e^{−y²−(1+2β)y/H}, y = (H/2)log(1 + v); |g0| ≤ H, ∫|G|dv = ∫|e^{−y²+(1−2β)y/H}|dy < ∞.
With W(x) := χ(x/T)x^{w−1} ∈ C_c^∞, Fubini: M = π^{−1/2}∫G(v)Ŵ(v)dv, Ŵ(v) = (T/2π)^w ω̂(Tv/2π), ω(t) := χ̃(t)t^{w−1}, |ω̂(k)| ≤ C_N(1 + |k|)^{−N}.
v = 2πk/T: M = π^{−1/2}(T/2π)^{w−1}∫_{k>−T/2π}G(2πk/T)ω̂(k)dk. For |v| ≤ 1/2 (|k| ≤ T/4π): |log(1 + v) − v| ≤ v², |g0'(v)| ≤ C_gH² (|ye^{−y²}| ≤ 1), so
  |G(2πk/T) − g0(0)e^{2πik}| ≤ C_gH²·2π|k|/T + (H/2)·T(2πk/T)² ≤ C(H²/T)(|k| + k²);
for |k| > T/4π use |G| ≤ H and ∫_{|k|>T/4π}|ω̂| ≤ C_NT^{−N}. Fourier inversion: ∫_ℝω̂(k)e^{2πik}dk = ω(1) = χ̃(1) = 1. Hence
  M(s) = π^{−1/2}(T/2π)^{w−1}[g0(0) + O(H²/T)],  g0(0) = H/2.  ∎
Consequence (SL l. 444–447, 709–728), from Lemma R alone. With (a, b) = (α, −β) the Estermann residues are q^{−1+a−b}ζ(1−a+b)F̃(1−a) + q^{−1−a+b}ζ(1+a−b)F̃(1−b),
F̃(s) = (q/p)^sM(s) (F(x) = χ(xp/(qT))A_{−1}(xp/q)); s ∈ {1 − α, 1 + β} has |w − 1| ≤ |d| ≤ 1/8. Let R_lead be this residue formula (with SL's outer factor,
the same for every profile) with M replaced by M_lead: it is SL's displayed expression, which SL obtains from exactly these two leading values (whole-v
integral √π, x-power at its centred value, l. 709–718). The exact residue sum (a contour around both poles, SL l. 448–449) and R_lead (poles cancel at
α + β = 0) are holomorphic; their difference is R_lead's formula with M_lead·E, so O((H/(gpq))·(H/T)L) on the circle |α + β| = C′/L (|ζ(1 ± (α+β))| ≪ L there),
hence inside by maximum modulus (in d with β fixed; enlarge the disks slightly before Cauchy); Cauchy for shift derivatives costs L^{O(1)}; the outer weights 1/(gpq)
sum to L^{O(1)}. Total change O(H^2T^{−1}L^{O(1)}) = o(H/L) for θ < 1: SL's residue sum (H√π/g)[ζ(1+α+β)q^{−1−α}p^{−1−β} + (T/2π)^{−α−β}ζ(1−α−β)q^{−1+β}p^{−1+α}]
and Young's diagonal (SL Lemma 2.2, needs only ν < 1/2) hold verbatim: main term √πH·c(α, β) unchanged. No Taylor profile enters this step: E1's
j ≥ 1 residue bounds ≪ H(T/H²)^jL⁴ are not needed, and no sum over j is formed. Float check: |M/M_pred − 1| ≤ 1.1·10⁻⁷ (§6).

## 4. Step (4): SL §3's profile-dependent branches for the exact profile
- Off-band (1 − χ(x/T))K_β(x): this is SL's own off-band term, unchanged (it is a piece of the exact kernel, not of a profile): for u = 2πx/T ∉ [2/3, 3/2] and
  |y| ≤ H/20 (e^{2y/H} ∈ [0.904, 1.106]) both y-phases 2Ty/H ± 2πxe^{2y/H} have |derivative| ≥ (T + x)/(5H); IBP gives O_A(x^{−Re β}(H/(T + x))^A), the
  u = e^{2y/H} form gives the summable large-x tail (SL l. 389–391). Error O(T^{1+ν+ε}(H/T)^A) ≪ HT^{−δ} for A(1 − θ) > 1 + ν − θ + δ: any θ < 1.
- Representation error SL (11): only the (H/T)^A part remains; the Taylor part T^{1+ν+ε}(T/H²)^J is gone (no expansion: the retained kernel is exactly
  χ(x/T)K_β(x) = χ(x/T)Σ_ε e(εx)A_ε(x), §1).
- ε = +1 retained term: on supp χ(x/T) the y-phase 2Ty/H + λ(y)x has derivative (2/H)(T + 2πxe^{2y/H}) ≥ 2T/H; IBP (amplitude e^{−y²+...}, smooth y-cutoff
  |y| ≤ H/20, tail e^{−H²/400}) gives |A_{+1}(x)| ≪_A (H/T)^A; summed absolutely over the band (≍ Tq/p terms, weight 1/(gq), σ ≪ n^ε): O(T^{1+ν+ε}(H/T)^A). No
  transformation, no residues needed (SL transforms it and bounds it the same way, l. 541–542).
- K-Bessel (cos(πd/2)) branch of (14): move u to Re u = σ large (no poles: Γ(u ± d/2) poles are at Re u ≤ |Re d|/2; Ã entire, decaying vertically since
  A ∈ C_c^∞). |Ã(1 − u − z/2)| ≤ ∫|A|ξ^{−σ−Re z/2}dξ ≪ T·(T/4π)^{−σ}·sup|A| ≪ T^{1−σ}(4π)^σ (sup|A| ≤ C, NO L¹ gain claimed). So the branch is
  ≪_σ T(pq/(nT))^σ per n; Σ_n (σ > 1), weights 1/(gpq), g, p, q ≤ y: ≪ T(y²/T)^σ = T^{1−σ(1−2ν)} ≪ T^{−A}. Same conclusion as SL (O(H(M²/T)^σ)).
- E1 Steps 1–2 (Mellin–Parseval, contour move, Hankel form): need only A ∈ C_c^∞(0,∞) and the shift disk: unchanged. E1 Step 3 is replaced by Lemma B.
- E1 Step 5 with Corollary B′: per (g, p, q) Σ_n|D| ≪ (1/(gpq))H(pq/T)^{1/4}N0^{3/4}T^{3ε} = g^{−1}T^{1/2+3ε}H^{−1/2}; Σ_{g,p,q} ≪ y²T^{1/2+3ε}H^{−1/2}.
  Relative to the main term ≍ H: T^{1/2+2ν−3θ/2+3ε}, power-small iff ν < (3θ − 1)/4 — now for every θ ∈ (1/3, 1) (θ > 1/2 no longer used).


## 5. Step (3): compact windows, and the chain to simple zeros
SL l. 730–758 (profile-free): w_D = φ_{s0}*Σ_{j<D}(−s0²/4)^jw^{(2j)}/j!, s0 = T^{θ−η}; ‖w − w_D‖_∞ ≪_D T^{−2Dη}L^{2D+O(1)}; absolute coefficient norm O(H/s0);
Gaussian moments of width s0 (centres T + O(H)) with error O(s0/L) or O(s0T^{−δ}) give O(H/L), O(HT^{−δ}); mass preserved exactly. The only input is the
Gaussian moment at width s0 = T^{θ′}, θ′ := θ − η, which §§2–4 give whenever ν < (3θ′ − 1)/4 (θ′ > 1/3): E1-M at θ′ (all bounds of §§2–4 are uniform
in the centre T′ ≍ T). Choice: ν = 121/1000; (3θ′ − 1)/4 > 121/1000 ⟺ θ′ > 371/750 = 0.494666…; for a given θ > 371/750 take η := (θ − 371/750)/2,
D > (1 + ν − θ + δ)/(2η) fixed. SL's requirement "θ − η > 1/2" (l. 731) was only for the Taylor step and the j ≥ 1 residues: both are gone.
Chain (E1 §3, E1A §3 rows unchanged, re-read for θ < 1/2): Thm 1.2 for Q = 1 and 1 + x^j ((Emean), l. 812–815: same P, R, Q(0) = 1) → Lemma 4.2 and Thm 1.3
(l. 823–928: N(T, H) = (HL/2π)(1 + o(1)) and Backlund–Jensen O(L) hold for H = T^θ, every θ > 0; N(T+H) − N(T) = (H/2π)log(T/2π) + O(H²/T + log T)) → Thm 5.1
(l. 932–963: Wang 2609.07918 Thm "local-pair", TeX l. 231–232, "Let 0 < λ < θ < 1 be fixed"; c(θ) defined for 0 < θ < 1, l. 125; c′(θ) = ½cot²(θ/√2) > 0,
l. 129) → liminf S(T, T^θ)/N(T, T^θ) ≥ h(θ; κ), h(θ; k) = (3 + k − √((1 − k)(9 − k − 8c(θ))))/4, positive iff k > k*(θ) = 1 − 2/(2 − c(θ)).
Detector: SL/EXP-010 frozen (P, Q) (exactly admissible: P(0) = 0, P(1) = 1, Q(0) = 1, Q(y) + Q(1 − y) ≡ 1), reused at ν = 121/1000, R = 44/5 (as nu12.py):
κ = 1 − R^{−1}log c(P, Q, R, ν). Every input of this chain other than Thm 1.2 is θ-free on (0, 1).


## 6. Float falsification (own scripts; single-threaded; memcap peaks in logs/*.err)
F1 (code/b3c_fals.py, logs/b3c_fals_part1.log + part2.log; peak 139 MB; part 1 killed by memcap at 309 MB on an unchunked Bessel array at T = 1e7,
θ = 0.40 — fixed, rows redone in part 2, nothing from the killed row used). Exact profile by FFT in u = H(e^{2y/H} − 1) (vs brute y-quadrature 2.0–2.1·10⁻⁵),
complex β kept in the amplitude, both Hankel branches separately, r = n/N0 on 29 log-points in [0.003, 8] (𝒯 depends on (n, pq) only via n/pq, so this
covers every n, p, q: N0 ≤ y²τ, r ≥ 1/N0), shifts (α, β) ∈ {(0.01, 0.02), (0.05+0.1i, −0.03+0.07i), (0.1+0.2i, −0.1−0.15i), (0.2−0.1i, 0.03+0.2i)}:
| T | θ: max_{r,shifts}|𝒯|/(HX0^{−1/2}) (τ) | max|𝒯_±| | r ≥ 4 | max|M/M_pred − 1| (H/T) |
|---|---|---|---|---|
| 1e6 | 0.400: 0.695 (15.8); 0.450: 0.702 (3.98); 0.470: 0.654 (2.29); 0.490: 0.672 (1.32); 0.495: 0.686 (1.15) | 0.353 | ≤ 5.3e−4 | ≤ 1.1e−7 (≥ 2.5e−4) |
| 1e7 | 0.400: 0.692 (25.1); 0.450: 0.704 (5.01); 0.470: 0.691 (2.63); 0.490: 0.705 (1.38); 0.495: 0.620 (1.17) | 0.353 | ≤ 5.1e−4 | ≤ 1.1e−8 (≥ 6.3e−5) |
| 1e8 | 0.400: 0.695 (39.8); 0.450: 0.672 (6.31); 0.470: 0.700 (3.02); 0.490: 0.700 (1.45); 0.495: 0.671 (1.20) | 0.353 | ≤ 5.4e−4 | ≤ 1.1e−9 (≥ 1.6e−5) |
Worst case over all 60 rows: 0.705 (T = 1e7, θ = 0.49, r = 0.003): flat in τ ∈ [1.15, 39.8] and T. max|𝒯_±|/(HX0^{−1/2}) = 0.3535 in every row with Re d = 0
= 1/(2√2) = 0.35355, the 2D stationary-phase value of §2's coordinates (det ≈ −1, f(0)α(1) = ½·1: (T/2π)·(2π/K)·½·2^{−1/2} = HX0^{−1/2}/(2√2)), reached as
r → 0; with Re d ≠ 0 it is times |(T/2π)^{−d/2}| (0.294 = 0.3536·(1.6·10⁵)^{−0.015} at T = 1e6, d = 0.03). Tail: ≤ 5.4·10⁻⁴ for r ≥ 4 (≈ 2^{−3/2}e^{−2π} = 6.6·10⁻⁴).
Lemma R: |M/M_pred − 1| ≤ 1.1·10⁻⁷, far below its bound H/T. No row exceeds Lemma B's explicit constant (F2) — nor its scaling.
F2 (code/b3c_lambda.py, logs/b3c_lambda.log): the proof's OWN coordinates with a genuine smooth cutoff χ̃ (smooth steps on [1/2, 2/3] and [3/2, 2]): inner
m(ρ) on 801 points of I_σ by s-quadrature, outer ∫_{I_σ}f·m·e^{iψ}dΛ, compared with F1's ξ-route; explicit constants M0, M1, B2, C1 of §2 evaluated for χ̃:
| T | θ | r | X0 | σ | |𝒯_σ|/(HX0^{−1/2}): ξ-route / Λ-route (rel. diff) | sup|m|X0^{1/2} ≤ M0 | Var m·X0^{1/2} ≤ M1 | max|J|X0²(1+|ρ|)² on N_σ ≤ B2 | C1 |
|---|---|---|---|---|---|---|---|---|---|
| 1e6 | 0.45 | 0.01 | 1000 | ±1 | 0.31426 / 0.31427 (1.6e−6) | 5.68 ≤ 54.4 | 14.8 ≤ 2.86e4 | 5.5e−2 ≤ 374 | 8.7e4 |
| 1e6 | 0.45 | 0.3 | 5479 | ±1 | 0.19930 / 0.19930 (≤ 3.7e−5) | 5.64 ≤ 54.4 | 11.8 ≤ 2.77e4 | 6.9e−7 ≤ 374 | 8.45e4 |
| 1e6 | 0.45 | 2 | 1.4e4 | ±1 | 0.01380 / 0.01380 (≤ 8.0e−5) | 5.63 ≤ 54.4 | 11.4 ≤ 2.76e4 | 1.6e−5 ≤ 374 | 8.41e4 |
| 1e6 | 0.40 | 0.01 | 1996 | ±1 | 0.31427 / 0.31427 (1.6e−6) | 5.66 ≤ 54.4 | 13.0 ≤ 2.81e4 | 6.2e−4 ≤ 374 | 8.54e4 |
| 1e6 | 0.40 | 0.3 | 1.1e4 | ±1 | 0.19932 / 0.19932 (≤ 6.4e−5) | 5.63 ≤ 54.4 | 11.5 ≤ 2.76e4 | 5.6e−6 ≤ 374 | 8.42e4 |
(shifts 0.05+0.1i, −0.03+0.07i; σ = +1 and −1 agree to 4 digits; θ = 0.40, r = 2 skipped: its grids would pass the 290 MB cap; peaks 272/233 MB.)
The Λ-route (Fubini + both substitutions + the I_σ/N_σ split, with χ̃) reproduces the ξ-route to ≤ 8·10⁻⁵: Step a is numerically exact. Each intermediate
bound of the proof holds with room (|m| 10×, Var m 2000×, J on N_σ ≥ 10⁴×); the explicit constant C1 ≈ 8.5·10⁴ is loose (observed ≤ 0.36 per branch).


## 7. Reviews (reviewer subagent; statement and proof of §§1–2 pasted inline each time)
- Pass 1: "GAP" for bookkeeping only — (a) displayed majorant of Var f false at H = 1 (3.868 > 3), (b) (ii) needs H ≤ T, (c) Cor B′ tail gave T^{−9}; it confirmed
  Fubini, both substitutions, the three phases, α and the Jacobian; Step b; Step c (|m|, |∂_ρm| = O(X0^{−1/2})); Step d's signs and constants under c_r ≤ H/46;
  "Crucially, stationary-range cancellation survives: m remains inside the outer oscillatory integral. Absolute integration occurs only in the nonstationary region";
  "the bulk argument proves (i) … uniformly in τ, n, p, q and the stipulated shifts"; "not a failure of the iterated cancellation argument". All three repaired.
- Pass 2: repairs (a), (b) correct; "No further gap found in Lemma B(i)–(ii)"; Cor B′ used N0 ≤ T (needs pq ≤ H²) — repaired with N0 ≤ T² only.
- Pass 3 (Cor B′ only): the summation with N0 ≤ T² "is valid"; the N-choice arithmetic needed N ≥ 2(A + 5)/ε (was 2(A + 4)/ε) and an unneeded aside
  (N0 < T needs θ > ν) — both fixed. Net: Lemma B (i)–(ii) CORRECT (passes 1–2), Corollary B′ correct after the pass-3 arithmetic fix.
- Pass 4 (Lemma R and §4, pasted inline): VERDICT CORRECT — Fubini (‖G‖₁‖W‖₁), uniform holomorphic O(H/T), pole cancellation at d = 0 including for M − M_0
  (maximum modulus in d with β fixed; use slightly enlarged shift disks before Cauchy), off-band with the large-x tail, ε = +1 pointwise negligible, K-branch with
  sup|A| only. §5 (profile-free) not sent. Lemma R also checked numerically to 10⁻⁷–10⁻⁹ (F1).
- Repair (orchestrator referee, 20 Nov): an independent referee found Lemma B (i)–(ii), Corollary B′ (n-sum, dual-total exponent) and Lemma R's direct
  Fourier-inversion proof CORRECT, and one removable GAP: §1 asserted Σ_j A_{−1,j} = A_{−1}, false as a convergent series (β = 0, ξ = T/2π:
  |A_{−1,j}(ξ)| ≳ (T/2)^je^{j²/H²}/j! for large j, so the terms do not tend to 0; exp(iφ)'s Taylor series converges only pointwise in y, and
  e^{|φ|} = exp(2πξ|e^{2y/H} − 1 − 2y/H|) is not integrable against e^{−y²}, so Σ_j and ∫dy do not commute), and §3 identified the sum of E1's
  j ≥ 1 residues with M − M_0. Both statements are deleted: §1 uses A_{−1} only; §3 derives the residues from Lemma R alone (M = M_lead(1 + E) against
  SL's displayed R_lead), with M − M_0 kept only as a direct comparison of two explicit integrals. Checked that no other step uses the series: Lemma B uses
  A_{−1} and ∂_ξ^kA_{−1} (differentiation under ∫dy, dominated by |λ(y)|^ke^{−y²+O(|y|/H)}); §4 uses the exact identity K_β = Σ_ε e(εx)A_ε (two terms,
  cos = ½Σ_ε e^{±i·}, absolutely convergent integrals), pointwise bounds for A_{+1} and sup|A_{−1}|; §5 and the chain are profile-free; F1/F2 compute A_{−1}
  directly (FFT and quadrature of the y-integral), never a series. The theorem is unchanged.


## 8. Step (5): Arb certificates (code/b3c_arb.py, logs/b3c_arb.log; python-flint, 4000 bits, exact rationals for every integral; rerun identical; peak 18.6 MB)
- Detector (SantiAudit/params_exp010.json, ν = 999/10000 set, nu12.py's construction): P(0) = 0, P(1) = 1, Q(0) = 1, Q(y) + Q(1 − y) ≡ 1 exactly, deg Q = 201.
- Check: κ(3/25, R = 177/20) = 0.0844779079697484032032447429603 ± 4·10⁻³² = E1's κ1 to every digit.
- κ := κ(121/1000, R = 44/5) = 0.0850387603917657602445070238217 ± 4.8·10⁻³² (c = 3138.99796501994347827760614809 ± 1.4·10⁻²⁷).
- Moment edge exact: (3·371/750 − 1)/4 = 121/1000.
- θ0 = 371/750: c(θ0) = −0.18577184013931683646750 ± 4·10⁻²⁶, k*(θ0) = 0.084991414349759443096037 ± 4·10⁻²⁷, κ − k*(θ0) = 4.7346042006317148·10⁻⁵ ± 4·10⁻²⁵ > 0,
  h(θ0; κ) = 3.354573606042305490255·10⁻⁵ ± 1.5·10⁻³⁰ > 0.
- θ = 4947/10000: h = 7.069262334884751296·10⁻⁵ (> 0); θ = 1/2: h = 0.005931006257946888362 (E1's frozen κ1 gave 0.0055353).
- Detector threshold: k*(θ) = κ at θ_det ∈ 0.494636567599178 ± 7·10⁻¹⁷ (bisection with Arb comparisons, k* strictly decreasing since c′ = ½cot²(θ/√2) > 0)
  < 371/750: the moment edge binds.
Conclusion. For fixed θ > 371/750: η := (θ − 371/750)/2 gives ν = 121/1000 < (3(θ − η) − 1)/4 (§5), so SL Thm 1.2 holds on (T, T + T^θ] with (P, Q, R = 44/5);
SL Thm 5.1 then gives liminf S/N ≥ h(θ; κ) ≥ h(371/750; κ) = 3.3546·10⁻⁵ > 0 (h increasing in θ since c is).

