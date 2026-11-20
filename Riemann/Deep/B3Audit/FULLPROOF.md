# B3Audit — FULLPROOF (b3audit, Orca, 20 Nov 2026, 05:40– UTC): line-by-line referee of B3Cert Theorem E1′
Object: Riemann/Deep/B3Cert/FULLPROOF.md (= B3F, § and line numbers of that file), REPORT.md. SL = Santibáñez-Leal v0.03, main.tex of Zenodo 23135354
(l. = its lines, the numbering B3Cert and E1Audit use); Wang = arXiv 2609.07918 v2 TeX. Notation as B3F §1: H = T^θ, K = T/H, τ = T/H², X0 = 2√(2π)√(nT/pq),
N0 = pqT/H², r = n/N0, c_r = X0/K = 2√(2π)√r, d = α + β, z = α − β, |α|, |β| ≤ C0/L.

## 1. Lemma B (B3F §2, l. 16–65) — re-derived, CORRECT
1.1 Kernel split (B3F §1 l. 7–15 vs SL l. 357–360): cos(2πxe^{2y/H}) = ½Σ_ε e(εx)e^{iελ(y)x}, λ = 2π(e^{2y/H} − 1), gives K_β = Σ_ε e(εx)A_ε exactly; both
    y-integrals converge absolutely (|e^{−y²+(1−2β)y/H}| integrable). ✓
1.2 Step a. ξ = Tu/2π, Λ = H(e^{2y/H} − 1): λ(y)ξ = T(e^{2y/H} − 1)u = KΛu; 2Ty/H = T log(1 + Λ/H) = Kℓ(Λ); X(ξ) = 4π√(nTu/(2πpq)) = X0√u; c_r = X0H/T
    = 2√(2π)√(nH²/(pqT)) ✓. Jacobians: dξ = (T/2π)du, dy = dΛ/(2(1 + Λ/H)) = ½e^{−2y/H}dΛ, so the y-amplitude e^{−y²+(1−2β)y/H}dy = f dΛ with
    f = ½e^{−y²−(1+2β)y/H} ✓; X^{−1/2} = X0^{−1/2}u^{−1/4}; ξ^{−β−z/2} and a_σ go into α(u); χ(ξ/T) = χ̃(u) ✓. Prefactor κ_σπ^{−1/2}(T/2π)X0^{−1/2},
    |κ_σ|π^{−1/2} = 2^{−1/2} ✓. Phase K[ℓ − Λu + σc_r√u] real; all shift dependence sits in the amplitudes f, α ✓. Fubini: absolute integrability ✓.
    |f| ≤ ½e^{(1+2|β|)²/(4H²)} ≤ ½e^{0.39} < 1 for H ≥ 1, |β| ≤ 1/8 ✓; Var f ≤ 2 (complex c, completed square) ✓ (bounded variation is parametrization-free).
1.3 Step b. φ(u) = σ√u − ρu, ρ = Λ/c_r; σ/(2√u) ∈ σ[0.3536, 0.7071] on [1/2, 2]. Cases (own check): σρ ≤ 0: |φ′| ≥ 0.3536 + |ρ|; 0 < σρ ≤ 0.1768:
    |φ′| ≥ 0.1768 ≥ (1 + 0.1768)/8 = 0.147; σρ ≥ √2: |φ′| ≥ |ρ| − 0.7071 ≥ (1 + |ρ|)/8 ⟺ |ρ| ≥ 0.951 ✓. φ″, φ‴ bounded independently of ρ, so
    ((α/φ′)′/φ′)′ = α″/φ′² − 3α′φ″/φ′³ − αφ‴/φ′³ + 3αφ″²/φ′⁴ ≪ (1 + |ρ|)^{−2} ✓. ∫_{N_σ}|f||J|dΛ ≤ B2X0^{−2}c_r∫(1 + |ρ|)^{−2}dρ = 2B2X0^{−2}c_r
    = 2B2K^{−1}X0^{−1} ✓ (X0 = Kc_r).
1.4 Step c. u = s²: σs − ρs² = −ρ(s − s0)² + 1/(4ρ), s0 = σ/(2ρ) (σ² = 1) ✓ exact, so J = e^{iψ0}m with ψ0 = X0/(4ρ) = Kc_r²/(4Λ) ✓.
    (c1) partial integration against E(x) = ∫_a^x e^{iF} (β̃ vanishes at the ends) + Titchmarsh Lemma 4.4 (|F″| = 2X0|ρ|) ✓.
    (c2) ∂_ρ[ρ(s − s0)²] = (s − s0)² + 2ρ(s − s0)·σ/(2ρ²) ✓; (s − s0)e^{−iX0ρ(s−s0)²} = (i/(2X0ρ))∂_s e^{…} ✓; the three bounds
    8(2X0|ρ|)^{−3/2}B̃2, 8(2X0|ρ|)^{−3/2}B̃1 + 8(2X0|ρ|)^{−5/2}B̃3 re-derived ✓; so |∂_ρm| ≪ X0^{−1/2} on |ρ| ≥ ρ_a, Var m over I_σ (length ≤ 1.24 in ρ)
    ≪ X0^{−1/2} ✓ (total variation is invariant under Λ = c_rρ).
1.5 Step d. ψ = Kℓ + Kc_r²/(4Λ): ψ″ = K[−H^{−1}(1 + Λ/H)^{−2} + c_r²/(2Λ³)] ✓; |c_r²/(2Λ³)| = 1/(2|ρ|³c_r) ≥ 1/(4√2c_r) on |ρ| ≤ √2 ✓; sign(Λ) = σ on I_σ.
    σ = −1: both terms < 0 ✓. σ = +1: 4/H ≤ 1/(8√2c_r) ⟺ c_r ≤ H/(32√2) = H/45.25 ✓ (and |Λ| ≤ √2c_r < H/2 keeps I_σ inside (−H, ∞)). Real phase, one
    sign of ψ″ on the interval: Titchmarsh 4.4 with r0 = K/(8√2c_r) ✓; amplitude fm handled by ∫_a^b g e^{iψ} = g(b)E(b) − ∫g′E, so the bound
    (sup|fm| + Var fm)·8r0^{−1/2} ✓ (boundary term kept: the endpoints of I_σ are not zeros of fm). Var(fm) ≤ sup|f|Var m + sup|m|Var f ≤ (M1 + 2M0)X0^{−1/2} ✓.
    (c_r/K)^{1/2}X0^{−1/2} = (c_r/K)^{1/2}(Kc_r)^{−1/2} = K^{−1} ✓. Hence |𝒥_σ| ≤ C·K^{−1} and |𝒯_σ| ≤ (C/(2π√2))HX0^{−1/2} ✓.
    Mechanism check: m stays inside the outer oscillatory integral (no absolute value is taken on I_σ) ✓; no step uses τ, θ or e^{−y²} beyond sup|f|, Var f ✓.
1.6 (ii). ∂_ξ^k under ∫dy brings (−iλ(y))^k, |λ| ≤ (4π|y|/H)e^{2|y|/H} (|e^t − 1| ≤ |t|e^{|t|}) ✓, so |∂_ξ^kA| ≤ C_kH^{−k} pointwise (T^{−1} ≤ H^{−1}) ✓;
    N integrations by parts against e^{iσX}, X′ = X/(2ξ) ≍ X0/T, each costing (T/X0)·H^{−1} = 1/c_r ✓; ξ-length ≤ T, |A| ≤ C ✓.
1.7 Corollary B′ (B3F l. 55–64): first range n ≤ N0T^ε has c_r ≤ 2√(2π)T^{ε/2} ≤ H/46 for T ≥ T0(ε, θ) ✓; second range (ii) with T = HK,
    (N0/n)^{N/2} ≤ (N0/n)^{K′}T^{−ε(N/2−K′)} ✓; sums: Σ_{n≤Y}d(n)n^{−1/4+ε} ≪ Y^{3/4+ε}L ✓; tail with N0 ≤ y²T ≤ T² (y² ≤ T) ✓, exponent 3 + 6ε − εN/2 ✓.
    Dual total: H(pq/T)^{1/4}N0^{3/4} = pqT^{1/2}H^{−1/2} ✓; with weights 1/(gpq), Σ_{g≤y}g^{−1}(y/g)² ≪ y² ✓; relative to H: T^{1/2+2ν−3θ/2+ε} ✓,
    power-small ⟺ ν < (3θ − 1)/4 ✓.
1.8 Sharpness (own 2D stationary phase, B3F coordinates): stationary point ℓ′(Λ) = u, Λ = σc_r/(2√u) ⟹ Λ0 ≈ σc_r/2, y0 ≈ σc_r/4; Hessian det = ℓ″F_uu − 1 ≈ −1;
    |𝒯_σ| ≈ 2^{−1/2}(T/2π)X0^{−1/2}(2π/K)|f(Λ0)| = HX0^{−1/2}e^{−πr/2}|(T/2π)^{−d/2}|/(2√2) (c_r²/16 = πr/2), uniformly in τ. So HX0^{−1/2} is the true size
    (Lemma B cannot be improved termwise) and the bound has no hidden τ-dependence; tested numerically in §6.
1.9 Uniformity: 𝒯_σ is holomorphic in (α, β) (A_{−1} in β, ξ^{−z/2}, a_σ(X, d) in d), Lemma B holds on the enlarged disk |α|, |β| ≤ 2C0/L, so Cauchy on
    circles of radius C0/L bounds the shift derivatives of fixed order (deg Q = 201) with L^{O(1)} loss ✓; the centre enters only through T in K_β's phase,
    χ(x/T′) and X0, and every constant depends on χ, C0, N only, so T′ ≍ T is uniform ✓; θ enters only via H in the hypothesis c_r ≤ H/46 and T0(ε, θ) ✓.

## 2. Lemma R and the residues (B3F §3, l. 66–87) — re-derived, CORRECT
2.1 v = e^{2y/H} − 1: λ(y)x = 2πvx, e^{2iTy/H} = (1 + v)^{iT}, dy = (H/2)dv/(1 + v), so A_{−1}(x) = x^{−β}π^{−1/2}Ĝ(x), G = g0(v)(1 + v)^{iT}1_{v>−1},
    g0 = (H/2)e^{−y²−(1+2β)y/H} ✓, |g0| ≤ H ✓, G ∈ L¹ ✓. Fubini (‖G‖₁‖W‖₁ < ∞): M = π^{−1/2}∫GŴ ✓. x = Tt/2π: Ŵ(v) = (T/2π)^wω̂(Tv/2π), ω = χ̃(t)t^{w−1} ✓;
    v = 2πk/T ✓. For |v| ≤ 1/2: T log(1 + v) = 2πk + T(log(1 + v) − v), |log(1 + v) − v| ≤ v², |∂_vg0| ≤ C_gH² ✓, so |G − g0(0)e^{2πik}| ≤ C(H²|k| + Hk²)/T ✓;
    |k| > T/4π: |G| ≤ H, ∫|ω̂| ≤ C_NT^{−N} ✓; ∫ω̂(k)e^{2πik}dk = ω(1) = χ̃(1) = 1 ✓ (ω̂ Schwartz uniformly on |w − 1| ≤ 1/4). M = M_lead(1 + O(H/T)) ✓.
    Size check: the profile has x-width ≍ K = T/H for θ < 1/2, so x^{w−1} varies by O(1/H) across it; the relative error is nevertheless O(H/T): in the
    Fourier picture Ŵ lives at |v| ≲ 1/T where G varies by O(H/T) relative (g0 on scale 1/H; the chirp e^{−iTv²/2} is 1 + O(1/T)). Numerics §6 (F4).
2.2 Residues: SL l. 441–447 (both residues, a = α, b = −β): q^{−1+α+β}ζ(1−α−β)F̃(1−α) + q^{−1−α−β}ζ(1+α+β)F̃(1+β), F̃(s) = (q/p)^sM(s) (x = tp/q) ✓.
    Own computation with M_lead and the outer factor 2π√(p/q)·1/(g√(pq)) = 2π/(gq) (SL (10) l. 363–367, b_{gq}b̄_{gp}/(gq) l. 413):
    term ζ(1+α+β): (2π/(gq))q^{−1−α−β}(q/p)^{1+β}H/(2√π) = (H√π/g)q^{−1−α}p^{−1−β}; term ζ(1−α−β): (H√π/g)(T/2π)^{−α−β}q^{−1+β}p^{−1+α}.
    This is SL's displayed expression l. 713–719 exactly ✓ (R_lead). Both R_exact (contour around both poles, l. 449) and R_lead (pole cancels at α + β = 0:
    q^{−1−α}p^{−1+α} − q^{−1−α}p^{−1+α} = 0) are holomorphic in d; their difference is the residue formula with M_lead·E, ≪ (H/(gpq))(H/T)L on |d| = C′/L,
    so inside by maximum modulus ✓; Σ_{g,p,q}1/(gpq) ≪ L³ ✓; total O(H²T^{−1}L^{O(1)}) = o(H/L) for θ < 1 ✓. Then SL l. 724–727 (Lemma diagonal + Young
    Lemma 7, θ-free: w enters only via ∫w g_{α,β}(s,t)dt, l. 293–300; needs ν < 1/2 only, l. 337–338) gives c(α, β)√πH ✓. The ε = +1 residues use F̃ of
    A_{+1} ≪ (H/T)^A (§3.3): negligible ✓. The polar term −2πζ(α+β)G(1−α)(h/k)^{1/2−α} of SL (10) has |G(1−α)| ≤ e^{−(T²−O(1))/H²}: negligible ✓.
2.3 The deleted series identity (B3F §7 Repair): rg of B3F for A_{ε,j}, A_{−1,j}, A_{−1,0}, Σ_j hits only l. 10 (disclaimer), l. 69–70 (the "direct comparison"
    M − M_0: two explicit integrals, no series, used in no bound) and l. 159–161 (repair note) ✓. Nothing in Lemma B, Lemma R, §4 or §5 uses a Taylor profile ✓.

## 3. Branches of SL §3 (B3F §4, l. 88–104) — CORRECT
3.1 Off-band (1 − χ(x/T))K_β: SL l. 389–391 unchanged; own check of |φ′_±| ≥ (T + x)/(5H) on |y| ≤ H/20 for u ∉ [2/3, 3/2] (u ≤ 2/3: (2T/H)·0.263 ≥ 0.221T/H;
    u ≥ 3/2: (2T/H)(0.904u − 1) ≥ (T/5H)(1 + u/2π)) ✓; the minus-phase stationary point y = −(H/2)log u has |y| ≥ 0.2H > H/20 ✓; tail |y| > H/20 costs e^{−H²/400}
    ≤ e^{−T^{2/3}/400} for θ ≥ 1/3 ✓. Error O(T^{1+ν+ε}(H/T)^A) ✓ for any θ < 1.
3.2 Representation error SL (11) l. 395–398: its (T/H²)^J part came only from the Taylor expansion (l. 383–391); with the exact two-term identity it is absent ✓.
3.3 ε = +1: φ′ = (2/H)(T + 2πxe^{2y/H}) ≥ 2T/H for all y (no cutoff needed), φ^{(k)} ≤ (2/H)^{k−1}φ′, φ″/φ′² = a/(T + a)² ≤ 1/(4T) (a = 2πxe^{2y/H}), so every
    integration by parts gains O((H/T)(1 + |y|)) ✓: |A_{+1}| ≪_A (H/T)^A, summed absolutely (≍ Tq/p terms, weight 2π/(gq)) ✓.
3.4 K-Bessel branch (cos(πd/2) term of SL (12), l. 434–437): Ã entire; Γ(u ± d/2) poles at Re u ≤ |Re d|/2; on Re u = σ, |ΓΓ| ≪ e^{−π|t|}(1+|t|)^{2σ−1} integrable and
    |Ã(1 − u − z/2)| ≤ sup|A|·∫_{T/4π}^{T/π}ξ^{−σ−Re z/2}dξ ≪ T^{1−σ} ✓ → T(y²/T)^σ·L ≪ T^{−A} ✓ (also directly: X ≥ X0/√2 ≫ T^{1/2−ν}, K_d(X) ≪ e^{−X}).
3.5 E1 Steps 1–2 (Mellin–Parseval, contour, Hankel form): exact identities for any A ∈ C_c^∞(0, ∞) (Voronoi for smooth compactly supported weights; the
    numerical end-to-end identity of E1Audit) ✓; Steps 3–4 replaced by Lemma B; Step 5 unchanged with Cor B′ (§1.7) ✓.


## 4. Compact windows and the chain to simple zeros (B3F §5, l. 105–118) — CORRECT
4.1 SL l. 731–732: "choose fixed η > 0 so that θ − η > 1/2 and ν < ν_*(θ − η). Choose the estimate that admits the fixed ν at this narrower exponent".
    The step itself (l. 733–758) is profile-free: w_D = φ_{s0}*Σ_{j<D}(−s0²/4)^jw^{(2j)}/j!, deficit ≤ u^D/D!, ‖w − w_D‖_∞ ≪ T^{−2Dη}L^{2D+O(1)} (w^{(j)} ≪ (H/L)^{−j}),
    replacement error by a coarse mollified-zeta bound on |t| ≤ 2T plus Gaussian tails "without using the desired moment", coefficient norm O(H/s0),
    centres T + O(H), mass preserved. Its only θ-input is the Gaussian moment at width s0 = T^{θ−η}. So "θ − η > 1/2" is exactly SL's own Gaussian range
    (SL Thm 1.2 needs θ > 1/2 through the Taylor profiles l. 383–398, the j ≥ 1 residues l. 722 and N0 < 1 l. 537), not a property of the window step ✓.
    With E1-M′ at θ′ = θ − η (§§1–3 bounds are uniform in the centre: T enters only through K_β's phase 2iT′y/H, the band χ(x/T′) and the shift disk; all
    constants depend on χ, C0, N only) the condition becomes ν < (3(θ − η) − 1)/4; with ν = 121/1000, η = (θ − 371/750)/2: 3(θ − η) − 1 > 3·371/750 − 1 =
    121/250 = 4ν ✓ (exact).
4.2 SL Thm 1.3 / lem:count (l. 823–925): O(T,H) ≥ N(T,H) − 2N_Ṽ − O(L) uses argument variation and Riemann–von Mangoldt "Δϑ/π = N(T,H) + O(L)" (l. 864) ✓;
    the proof of Thm signs (l. 890–920) uses Littlewood's lemma, |a + b|² ≤ (1 + L^{−1/2})|a|² + (1 + L^{1/2})|b|², Thm 1.2 and (Emean) (l. 812–815: Thm 1.2 for
    Q = 1 and 1 + x^j, same window and ν) ✓. N(T, T^θ) = (H/2π)log(T/2π) + O(H²/T + L) ~ HL/2π for every θ > 0 ✓; horizontal sides O(L) (Backlund) ✓.
    No θ > 1/2 anywhere in l. 760–928 (rg of "theta" in that range: only the η = θ − 1/2 − ν remark l. 925 about SL's own constants) ✓.
4.3 SL Thm D (transfer, l. 932–954): parity product (q_T − s_T)(1 − o_T) ≥ 2(1 − s_T)² [CAOSShort Prop 2.5] is a deterministic inequality on zero multisets
    (SantiAudit: random-multiset falsification, no θ) ✓; Wang thm:local-pair (TeX l. 232–240): "Let 0 < λ < θ < 1 be fixed", error O_g(H + T^λL²) ✓ θ-free
    (WangAudit: no RH, no zero-density, no simplicity input; only |δ| < 1/2); c(θ) defined for 0 < θ < 1 (TeX l. 125), c′ = ½cot²(θ/√2) > 0 ✓.
    Wang's own θ > 0.5501 (TeX l. 133) comes from his moment input, not from the pair theorem ✓. Hence liminf S/N ≥ h(θ; κ) given O/N ≥ κ − o(1) ✓.
4.4 Diagonal: SL Lemma diagonal (l. 285–300) and Remark rem:uses (l. 336–345): "ν < θ − 1/2 is used only in Lemma twisted" (the §2 off-diagonal IBP, l. 231–262,
    which the §3 signed route replaces), "ν < 1/2 only in Lemma diagonal"; θ < 1 only places supp w in [T/2, 2T] ✓. So below 1/2 the main term is
    residues (§2.2) + Lemma diagonal; §2's section-wide hypothesis "0 < ν < θ − 1/2" (l. 199) is not needed for θ < 1/2 ✓.


## 5. Arb certificates (B3F §8, l. 169–) — rerun identical; independent recomputation agrees to 30 digits
5.1 Rerun of a copy (code/b3c_arb.py, arg SantiAudit/params_exp010.json; uv --with python-flint; peak footprint 15.3 MB): logs/rerun_b3c_arb.log is byte-identical
    to B3Cert/logs/b3c_arb.log (diff empty): admissibility P(0) = 0, P(1) = 1, Q(0) = 1, Q(y) + Q(1 − y) ≡ 1, deg Q = 201; κ(3/25, 177/20) = 0.08447790796974840320…;
    κ(121/1000, 44/5) = 0.0850387603917657602445070238217 ± 4.8·10⁻³²; k*(371/750) = 0.084991414349759443096…; κ − k* = 4.7346042006317148·10⁻⁵ > 0;
    h(371/750; κ) = 3.354573606042305490·10⁻⁵ > 0; θ_det = 0.494636567599178 ± 7·10⁻¹⁷ < 371/750.
5.2 Independent recomputation (code/b3a_kappa.py, logs/b3a_kappa_deg6.log; mpmath 50 digits, peak 17 MB): different method — floating point instead of exact
    rationals; Q, Q′ from cos(nφ), n sin(nφ)/sin φ (u = 1 − 2x = cos φ) instead of polynomial coefficients; Gauss–Legendre (24 × 96 nodes) instead of exact
    antiderivatives; c = Q(0)² + ν^{−1}[CP·∫w_R² + 2ν∫PP′·∫w_Rw_R′ + ν²BP·∫w_R′²] (SL l. 174) with the cross term integrated (Q(1) = 0 not used; it returns
    ∫w_Rw_R′ = −0.5, ∫PP′ = 0.5, Q(1−) = 4·10⁻³⁵): κ(3/25, 177/20) = 0.0844779079697484032032447429603, κ(121/1000, 44/5) = 0.0850387603917657602445070238217,
    c = 3138.99796501994347827760614809, k*(371/750) = 0.08499141434975944309603663, κ − k* = 4.734604200631714847·10⁻⁵, h(371/750; κ) = 3.354573606042305490255193·10⁻⁵,
    θ_det (findroot) = 0.49463656759917793405: all equal to the Arb values in every printed digit ✓. Edge: (3·371/750 − 1)/4 = 121/1000 exactly ✓.


## 6. Falsification (own code; single-threaded; memcap 290 MB; peaks in logs/*.err)
F3 (code/b3a_lemB.py, logs/b3a_lemB.log, peak 149 MB; check run logs/b3a_lemB_check.log). Independent route: A_{−1}(ξ) by direct trapezoid in y on [−6.5, 6.5]
   with step 2π/(max|φ′| + 30) (Poisson aliasing ≤ e^{−225}), validated against mpmath.quad (400 subintervals, 25 digits) at 10 points: rel. error ≤ 3.9·10⁻¹³;
   𝒯_± by direct trapezoid in δ = 2πξ/T − 1 over |y0| ≤ 6.8 (|A| at the ends ≤ 1.4·10⁻¹³ max|A|), Hankel functions of complex order by their asymptotic
   series (vs mpmath hankel1/2: ≤ 1.3·10⁻¹² at X = 60…2·10⁴). 24 blocks = (T, θ) ∈ {10⁴, 10⁵, 10⁶} × {0.40, 0.45} ∪ {(10⁵, 0.35), (10⁶, 0.33)} (τ = 2.5 … 110)
   × shifts {(0, 0), (0.05+0.1i, −0.03+0.07i), (0.2−0.1i, 0.03+0.2i)} × r ∈ {0.003, 0.03, 0.1, 0.3, 1, 2, 4, 8} (192 rows; 𝒯 depends on n/pq only).
   max |𝒯|/(HX0^{−1/2}) = 0.668 (T = 10⁴, θ = 0.45, r = 0.003); every branch equals the 2D stationary-phase value e^{−πr/2}|(T/2π)^{−d/2}|/(2√2) of §1.8
   (r ≤ 1: to ≤ 10⁻³ relative; r = 4: ≤ 3.6%), e.g. r = 0.003: 0.35189 at every (T, θ), zero shifts. Flat in T (10⁴ → 10⁶) and τ (2.5 → 110): the scaling
   |𝒯| ∝ HX0^{−1/2} holds with constant 1/(2√2) per branch; the explicit C1 ≈ 8.5·10⁴ of B3F F2 is never approached. Rows outside Lemma B(i)'s hypothesis
   (c_r > H/46) obey the same bound (the hypothesis is a proof artefact; (ii) covers them).
F5 (logs/b3a_lemB_tau251.log, peak 99 MB): T = 10⁶, θ = 0.30, τ = 251, K = 1.6·10⁴ (42437 × 34454 grid, 122–250 s per profile), the three shift sets,
   r = 0.003 … 8 (24 rows): zero shifts |𝒯_±|/(HX0^{−1/2}) = 0.35189, 0.33728, 0.30216, 0.22070, 0.07350, 0.01529, 0.00066 = prediction to 5 digits; all rows:
   branches within 0.52% of the prediction, max |𝒯| ratio 0.580. τ-uniformity confirmed to τ = 251 (F3 + F5: 27 blocks, 216 rows, no failing case).
F4 (code/b3a_lemR.py, logs/b3a_lemR.log, peak 132 MB): Lemma R off w = 1 (w ∈ {1.2+0.1i, 0.8−0.15i, 1+0.25i, 1.25}, β ∈ {0, −0.03+0.07i}):
   |M/M_lead − 1| = 5.9·10⁻⁶, 5.5·10⁻⁷, 5.3·10⁻⁸ at T = 10⁴, 10⁵, 10⁶ (θ = 0.40), 6.5·10⁻⁸ (10⁶, 0.45), 5.0·10⁻⁸ (10⁶, 0.33): ≈ 0.035/T, below the bound
   H/T by 10³–10⁴ (the leading correction is g0′(0)/(g0(0)T) = O(1/T) since e^{−y²} is flat at 0) ✓. At w = 1 exactly M = H/(2√π) to 10⁻¹³ (Fourier inversion).
F6 (code/b3a_mpT.py, logs/b3a_mpT.log, peak 74 MB): one dual term entirely in mpmath (25 digits; mp.expj trapezoid, mp.hankel1/2), T = 10⁴, θ = 0.40, shifts 0,
   r = 0.003, 0.3: profile = float profile to 1.9·10⁻¹⁴; |𝒯_±|, |𝒯| ratios 0.35188668, 0.56378897 (r = 0.003), 0.22069143, 0.36443289 (r = 0.3) = float route ✓.
Rerun of B3Cert's own F1 (copy code/b3c_fals.py, jobs 10⁶:0.45, 10⁶:0.40, logs/rerun_b3c_fals.log, peak 79 MB): all 8 rows identical to B3Cert's logs ✓.


## 7. Review
- reviewer subagent (20 Nov, 06:03 UTC), Lemma B(i) Steps b–d pasted inline with exact formulas (statement, N_σ/I_σ split, completed square, (c1), (c2),
  ψ″, Titchmarsh 4.4 with non-vanishing amplitude, gain bookkeeping): "VERDICT: CORRECT … nonstationary derivative bound is valid; compact support eliminates
  the integration-by-parts boundary terms"; it re-derived ∂_ρm = −(1/2ρ)∫b̃E − (σ/2ρ²)∫b̃′E + (i/(4X0ρ²))∫b̃″E (own check: equals B3F (c2) after the two
  partial integrations ✓); "the formula and signs for ψ″ are correct; 46 > 32√2 ensures the positive-sign lower bound"; "the amplitude estimate allows complex,
  nonvanishing amplitudes: its retained boundary term handles the endpoints"; "|J_σ| ≲ K^{−1}, uniformly without an additional restriction on τ or θ";
  "the bulk condition c_r ≤ H/46 remains essential to this proof" (it is a hypothesis of (i), used in Cor B′ only for n ≤ N0T^ε, where it holds for T ≥ T0).
- Own adversarial points, all resolved: (a) SL §2's section hypothesis "0 < ν < θ − 1/2" (l. 199) is impossible below 1/2 — used only in Lemma twisted
  (rem:uses l. 337), which the §3 route replaces (§4.4); (b) SL l. 731 "θ − η > 1/2" — SL's Gaussian range only (§4.1); (c) the profile's x-width K = T/H ≫ H
  below 1/2 does not spoil Lemma R's O(H/T) (Fourier argument §2.1, F4); (d) the profile has curvature 4τ → ∞ in y, but Lemma B never resolves it (the
  stationary phase is 2D in (Λ, u) with Hessian ≈ −1, §1.8; F3/F5 flat to τ = 251); (e) no zero-density or pair-correlation input with θ > 1/2 in Wang's
  local pair theorem or in lem:count (§4.2–4.3).

