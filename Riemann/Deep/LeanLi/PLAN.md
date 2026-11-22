# LeanLi PLAN (T3: DISCOVERIES #3, Li-coefficient remainder laws)

Main: `LeanLi.li_remainder_laws : LiteratureInputs → LiRemainderLaws` (Statement.lean: Part0–Part5).
Objects (Defs.lean): λ_n = lim_{T→∞} Σ_{|ρ|≤T} m_ρ[1 − (1 − 1/ρ)ⁿ] (Lagarias' convention), R, c₁, Keiper y, Voros S̄_n, N(T), Titchmarsh N(T) (averaged), θ_RS (∫ of ½Re ψ(1/4+iu/2) − ½log π), S(T), Montgomery F(α,T).
Inputs (Lit.lean): ζ ≠ 0 on (0,1); θ_RS Stirling; DLMF 5.11.2 (ψ); S = O(log T) (Titchmarsh 9.3+9.4); S₁ = O(log T) (9.9); Selberg ∫S² ≪ T log log T; Montgomery 1973 (RH, α ∈ [0,1/2]); Voros (TS).
NOT-FORMALIZED: "a = 2 first fails at n = 5083" (needs certified λ_n, n ≤ 5083).

## Files / lemma DAG
1. ZerosBasic: strip 0<Re<1, ζ(s̄)=conj ζ(s), mult conj-invariant, mult ≥ 1, 1−ρ zero, Im ρ ≠ 0, finiteness. ✔
2. Abel: Σ m_i f(x_i) = f(b)Σm_i − ∫_a^b f'·cnt. ✔
3. ZerosCount: zset, N = finset sum, N ≤ 2·N_avg, non-ordinates exist, N_avg = N off ordinates. ✔
4. CountBounds: N_avg, N = (T/2π)log(T/2πe) + O(log T); jump O(log T); N ≤ T log T. ✔
5. LiCoeff: Σ_{γ>0} m/|ρ|² < ∞; |Re(1−wⁿ)| ≤ (n + n²Kⁿ)/|ρ|²; liPartial → Σ_U 2m Re(1−wⁿ) (Part0); RH ⇒ w = e^{iθ(γ)}, θ(γ)=π−2arctan 2γ, λ_n = Σ_κ 2(1−cos nθ_k) over κ = Σ ρ, Fin m_ρ.
6. Smooth: finite differences ≤ h^K sup|f^{(K)}|; summation by parts; |Σ_n φ(n/N) n^j e^{inx}| ≤ C N^{j+1}(1+N‖x‖)^{−K}.
7. ThetaRS: A = θ_RS/π + 1, A(0)=1, A' = (2π)^{-1}log(t/2π) + r, |r| ≤ C/t² (t ≥ 1), ∫_0^∞ r = −1/8.
8. LemmaA: I(n) := ∫_0^π 2(1−cos nv)a₀(v)dv = Σ_{k odd<n}(n−k)/k − (n/2)log 4π (cos·log sin integrals); I − t = O(1) (and O(1/n)); λ⁽⁰⁾ = I + 2∫r − 2J; Abel over zeros: Σ_γ f(γ) = ∫fA' + f(0) − ∫f'S; pointwise decomposition at T.
9. Thm1: weights w = autocorrelation of v(m) = ψ(m/N) (Φ = |V|²/Σv² ≥ 0); decay of Φ, Φ''; model trend I (|I − t| ≤ C): τ_I(a) = ∫_0^π B(a,v)a₀(v)dv ≥ log(2πa)/a² − C log N; TauT by Abel + RvM; TinyTerms; LiRate core (copied) ⇒ Thm1Bound(λ − I) ⇒ Thm1Bound(R).
10. Thm2: R(n) = −2Re Z_T(n) − E(n) + O(log N) (T² log T = N); (a) Σφ|Z_T|² ≪ N·N(T)·log T (local count); Montgomery: ≪ N·N(T) (majorant |∫_{−1/4}^{1/4} e^{2πiαx}dα|²); (b) Schur + Selberg: ≪ N² log log T/T.
11. Thm4: Σφλ − Σφλ⁽⁰⁾ = G_N(0) − ∫G_N'S = G_N(0) + ∫G_N''S₁; bounds ⇒ O(log N); Σφ(λ⁽⁰⁾ − t + 1/4) = O(1) ⇒ 7/4.
12. Lemma0: ¬RH ⇒ ρ₀ with Re < 1/2; dominant |w| = M > 1 (finitely many with |w| ≥ 1+ε); mean square of Σ m e^{inα} ⇒ |λ_n| ≥ Mⁿ/4 i.o.
13. Cor (LiRate Cor/Cor3 copied): Cor1–4, Ω, RH ⟺ weighted ℓ², Part5.
14. Main: assemble; #print axioms.
