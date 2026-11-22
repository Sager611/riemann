# LeanLi PLAN (T3: DISCOVERIES #3, Li-coefficient remainder laws) — all steps ✔

Main: `LeanLi.li_remainder_laws : LiteratureInputs → LiRemainderLaws` (Main.lean; Statement.lean: Part0–Part5).
Objects (Defs.lean): λ_n = lim_{T→∞} Σ_{|ρ|≤T} m_ρ[1 − (1 − 1/ρ)ⁿ] (Lagarias' convention), R, c₁, Keiper y, Voros S̄_n, N(T), Titchmarsh N(T) (averaged), θ_RS (∫ of ½Re ψ(1/4+iu/2) − ½log π), S(T), Montgomery F(α,T).
Inputs (Lit.lean): LitCore = ζ ≠ 0 on (0,1); θ_RS Stirling; DLMF 5.11.2 (ψ); S = O(log T); S₁ = O(log T); Selberg ∫S²; Voros (TS). LiteratureInputs = LitCore + Montgomery 1973 (RH, α ∈ [0,1/2]).
NOT-FORMALIZED: "a = 2 first fails at n = 5083" (needs certified λ_n, n ≤ 5083).

## Files / lemma DAG
1. ZerosBasic, Abel, ZerosCount, CountBounds, ZeroSums: zeros, counts, N(T) = (T/2π)log(T/2πe) + O(log T). ✔
2. LiCoeff (Part0), LiRH (Li's formula under RH: λ_n = Σ_k 2(1−cos nθ_k)). ✔
3. Smooth, Weights, PhiFacts: Fourier decay of smooth sums; positive-definite weights. ✔
4. LogCos, IClosed, ITrend, IFine, ChangeVar: I(n) closed form, I − t = O(1/n), I(n) = π⁻¹∫(1 − cos nθ)log(t/2π). ✔
5. A0, TrendI, TrendLow, TrendFinal, TrendPt2, LiTheta, TauT, Tails, Tiny, LR/*, Thm1: Theorem 1 (c₁). ✔
6. ThetaRS, AbelS, LemmaA2, Decomp: Lemma A, R = −2W + 2na + Rest, Rest = O(log N). ✔
7. Bessel, SelbergT, MeanW, Thm2: Theorem 2 (N^{3/2}(log N)^{3/2}). ✔
8. PairCorr (Fejér kernel + Montgomery ⇒ Σ(1+log T|γ−γ'|)^-2 ≪ T log T), Thm2M (log log N bound). ✔
9. SmoothSum, IbpS, Thm4a, Thm4b, Thm4, Thm4c: Theorem 4 (7/4), real scales N ≥ 2. ✔
10. Lemma0: ¬RH ⇒ |λ_n| ≥ e^{cn} i.o. (max |1−1/ρ| > 1 on a finite set; cosine sums with positive coefficients). ✔
11. Corollaries, Part13, Main: Claims 1–5 assembled; #print axioms standard. ✔

