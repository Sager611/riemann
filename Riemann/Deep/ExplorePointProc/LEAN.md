FORMALIZED: step (iii′), the finite hard-core pair-sum inequality (pairSum_le; weighted form pairSum_weighted_le), from no hypotheses; Theorem 1 (theorem1, theorem1_D: 0 < D(λ) with Lamzouri's liminf D) and Lamzouri's Thm Main2 as its case A = β = c = 0 (thmMain2), from the named analytic hypotheses Steps12 (Lamzouri's (Eq:ReductionGoodSet2), TeX l.425–446), FourierCc2, LemmaT′, ScaleToInf, RvM and the zero model's Zeros.Ordered; the corollary D(λ) > 0 for every λ ≥ 96/125 = 0.768 (D_pos_ge_768) from these plus the named certificate facts CertP, CertQ (exact Sturm), CertJ, CertR1, CertR2 (Arb) and HermRep (Hermite–Gaussian form of f, f^). Proved, not assumed: the rational margin, region R3, the certificate's piecewise-linear A, f > 0, f^ < 0 on [U, ∞), f^(0) ≤ 1.0011, Z_0's hard-core property, count and reindexing, #E_λ ≤ 2·#(small gaps), DPos ⟺ 0 < liminf, monotonicity in λ. Sorry-free; axioms propext, Classical.choice, Quot.sound only; no native_decide, no decide.

# ExplorePointProc in Lean 4 + Mathlib (pplean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile). Riemann/Deep/ExplorePointProc/Lean/, 1,304 lines, namespace PPLean; not imported by Riemann.lean. No Mathlib zeta, Fourier, LSeries, log, π or exp (290 MB cap): the test-function calculus and log T/2π are abstract, every analytic input is a named hypothesis.

| file | content | lines | peak |
|---|---|---|---|
| FinSum | HH (= (H)), pairSum, (iii′) pairSum_le, pairSum_weighted_le | 152 | 254 MiB |
| Model | zero model Zeros, Zeros.Ordered; smallGaps, DPos, Elam (E_λ), Z0, xs; card_Elam_le, DPos.mono, card_Z0, Z0_sep, pairSum_Z0, xs_gap, xs_bounds; TestCalc; Steps12, FourierCc2, LemmaT′, ScaleToInf, RvM; core_arith | 260 | 264 MiB |
| NoGood | weight_le, weight_unfold, last_term, eta_small, eps_small, exists_delta; no_good_height (the contradiction at one height T) | 186 | 274 MiB |
| Theorem1 | theorem1_Cc2, theorem1, thmMain2, HH_div, margin_div | 114 | 258 MiB |
| PL | plin, plin_tail, plin_mem, nodesA (78 nodes), Acert, Acert_tail, Acert_range | 187 | 268 MiB |
| Herm | horner, pcoef, qcoef (22 exact coefficients each), CertP, CertQ, HermRep, pos_of_rep, neg_of_rep, fh0_of_rep | 120 | 245 MiB |
| Cert | CertPos, CertNeg, CertJ, CertF0, CertA, Phi, CertR1, CertR2, certHH (R3 proved), cert_margin, dpos_768, dpos_ge_768, dpos_ge_768_exact | 143 | 255 MiB |
| Liminf | Dlam (Mathlib liminf), dpos_iff, nToInf | 85 | 250 MiB |
| Final | theorem1_D, D_pos_ge_768 | 32 | 247 MiB |
| Check | #print axioms (23 theorems) | 25 | 245 MiB |

## Main statements

    -- (H) and the pair sum (FinSum)
    def HH (lam β c : ℝ) (F A : ℝ → ℝ) : Prop :=
      ∀ x y : ℝ, lam ≤ x → lam ≤ y → F x + F (x + y) + A x - A y + β * (1 - x) ≤ c
    def pairSum (F : ℝ → ℝ) (x : ℕ → ℝ) (n : ℕ) : ℝ := ∑ i ∈ range n, ∑ j ∈ range n, F (x i - x j)

    -- (1) the new step (iii′), 0-indexed: x_{n−1} − x_1 ↦ x (n-2) − x 0, A(g_{n−1}) − A(g_1), F(g_{n−1})
    theorem pairSum_le {lam U β c : ℝ} {F A : ℝ → ℝ} (hlam : 0 ≤ lam) (hU : U ≤ 3 * lam)
        (heven : ∀ u, F (-u) = F u) (hneg : ∀ u, U ≤ u → F u ≤ 0) (hH : HH lam β c F A)
        {n : ℕ} (hn : 2 ≤ n) (x : ℕ → ℝ) (hgap : ∀ i, i + 1 < n → lam ≤ x (i + 1) - x i) :
        pairSum F x n ≤ (n : ℝ) * F 0 + 2 * ((n : ℝ) - 2) * c
          + 2 * β * ((x (n - 2) - x 0) - ((n : ℝ) - 2))
          + 2 * (A (x (n - 1) - x (n - 2)) - A (x 1 - x 0)) + 2 * F (x (n - 1) - x (n - 2))
    theorem pairSum_weighted_le ... (hweven : ∀ u, w (-u) = w u) (hwnn : ∀ u, 0 ≤ w u) (hw0 : w 0 = 1)
        (hη : ∀ u, F u * w u ≤ F u + η) ... :
        pairSum (fun u => F u * w u) x n ≤ (n : ℝ) * F 0 + 2 * ((n : ℝ) - 2) * (c + 2 * η) + …
          + 2 * (F (x (n - 1) - x (n - 2)) * w (x (n - 1) - x (n - 2)))   -- same other terms

    -- zero model (Model): ordinates with multiplicity, N(T), sc T = log T/(2π)
    structure Zeros where  γ : ℕ → ℝ;  N : ℝ → ℕ;  sc : ℝ → ℝ
    structure Zeros.Ordered (Z : Zeros) : Prop where
      mono : ∀ k : ℕ, 1 ≤ k → Z.γ k ≤ Z.γ (k + 1)
      pos : 0 < Z.γ 1
      count : ∀ (T : ℝ) (k : ℕ), 1 ≤ k → (Z.γ k ≤ T ↔ k ≤ Z.N T)
    def smallGaps (Z) (lam T) : Finset ℕ := (Icc 1 (Z.N T)).filter (fun k => (Z.γ (k + 1) - Z.γ k) * Z.sc T ≤ lam)
    def DPos (Z) (lam) : Prop := ∃ δ : ℝ, 0 < δ ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → δ * (Z.N T : ℝ) ≤ ((smallGaps Z lam T).card : ℝ)
    def Dlam (Z) (lam) : ℝ := liminf (fun T : ℝ => ((smallGaps Z lam T).card : ℝ) / (Z.N T : ℝ)) atTop
    def Elam (Z) (lam T) : Finset ℕ := (Icc 1 (Z.N T)).filter
        (fun k => 2 ≤ k ∧ min (Z.γ (k + 1) - Z.γ k) (Z.γ k - Z.γ (k - 1)) * Z.sc T ≤ lam)
    def Z0 (Z) (lam T) : Finset ℕ := Icc 1 (Z.N T) \ Elam Z lam T
    lemma card_Elam_le : (Elam Z lam T).card ≤ 2 * (smallGaps Z lam T).card
    lemma card_Z0 : (Z0 Z lam T).card + (Elam Z lam T).card = Z.N T
    lemma Z0_sep (hZ : Z.Ordered) (hsc : 0 ≤ Z.sc T) (ha : a ∈ Z0 Z lam T) (hb : b ∈ Z0 Z lam T) (hab : a < b) :
        lam < (Z.γ b - Z.γ a) * Z.sc T
    lemma pairSum_Z0 (Φ : ℝ → ℝ) : ∑ a ∈ Z0 Z lam T, ∑ b ∈ Z0 Z lam T, Φ ((Z.γ a - Z.γ b) * Z.sc T)
        = pairSum Φ (xs Z lam T) (Z0 Z lam T).card          -- xs = increasing unfolded ordinates of Z_0
    theorem DPos.mono (h : lam ≤ lam') (hD : DPos Z lam) : DPos Z lam'
    theorem dpos_iff (hN : ∀ K : ℝ, ∃ T₀, ∀ T, T₀ ≤ T → K ≤ (Z.N T : ℝ)) (lam : ℝ) : DPos Z lam ↔ 0 < Dlam Z lam

    -- (3) Theorem 1 (Theorem1, Final)
    theorem theorem1 (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z) (h12 : Steps12 Z C)
        (hFC : FourierCc2 C) (hT : LemmaT' C) {θ lam U β c : ℝ} {f A : ℝ → ℝ}
        (hθ1 : 1 / 2 < θ) (hθ2 : θ < 4 / 7) (hlam0 : 0 < lam) (hlam1 : lam ≤ 1) (hU1 : lam < U) (hU3 : U ≤ 3 * lam)
        (hfS : C.Sch f) (hfe : ∀ x, f (-x) = f x) (hfn : ∀ x, 0 ≤ f x) (hfneg : ∀ u, U ≤ u → C.hat f u ≤ 0)
        (hA : ∃ K : ℝ, ∀ y, lam ≤ y → |A y| ≤ K) (hβ : 0 ≤ β) (hH : HH lam β c (C.hat f) A)
        (hJ : C.hat f 0 + 2 * c < C.J θ f) : DPos Z lam
    theorem theorem1_D … (same hypotheses) : 0 < Dlam Z lam
    theorem theorem1_Cc2 … (F ∈ C_c² directly, no LemmaT′, only U ≤ 3λ) : DPos Z lam
    theorem thmMain2 … (hfc : C.Cc2 f) … (hfneg : ∀ u, lam ≤ u → C.hat f u ≤ 0) (hJ : C.hat f 0 < C.J θ f) : DPos Z lam
    lemma HH_div (ht : 0 < t) (h : HH lam β c F A) : HH lam (β / t) (c / t) (fun u => F u / t) (fun u => A u / t)
    lemma margin_div (ht : 0 < t) (h : F0 + 2 * c < J) : F0 / t + 2 * (c / t) < J / t

    -- (4) corollary (Cert, Final); Acert = the certificate's A, pcoef/qcoef = its p, q
    theorem D_pos_ge_768 (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z) (h12 : Steps12 Z C)
        (hFC : FourierCc2 C) (hT : LemmaT' C) {f : ℝ → ℝ} (hS : C.Sch f) (he : ∀ x, f (-x) = f x)
        (hrep : HermRep f (C.hat f)) (hP : CertP) (hQ : CertQ) (hJ : CertJ C f)
        (hR1 : CertR1 C f Acert) (hR2 : CertR2 C f Acert) : ∀ lam : ℝ, 96 / 125 ≤ lam → 0 < Dlam Z lam
    theorem dpos_ge_768 … (hpos : CertPos f) (hneg : CertNeg C f) (hJ : CertJ C f) (hF0 : CertF0 C f)
        (hA : CertA A) (hR1 : CertR1 C f A) (hR2 : CertR2 C f A) : ∀ lam : ℝ, 96 / 125 ≤ lam → DPos Z lam
        -- the brief's literal form: f > 0, f^ < 0, f^(0), A as named hypotheses (all four derived in D_pos_ge_768)
    lemma certHH (hneg : CertNeg C f) (hA : CertA A) (hR1 : CertR1 C f A) (hR2 : CertR2 C f A) :
        HH (96 / 125) (416188506333 / 640725853462) (-347861653223 / 1000000000000) (C.hat f) A
    lemma cert_margin (hJ : CertJ C f) (hF0 : CertF0 C f) :
        C.hat f 0 + 2 * (-347861653223 / 1000000000000) < C.J (399993 / 700000) f      -- norm_num

## Hypotheses and sources

Lamzouri 2611.14486v1, SmallGapsPaper.tex (curl of https://arxiv.org/e-print/2611.14486, gzip/tar, 574 lines; (Lamzouri TeX, re-fetch https://arxiv.org/e-print/2611.14486; not kept)/). DERIVATION = ../DERIVATION.md. cert_aug.py = ../code/cert_aug.py.

| hypothesis | Lean statement | source |
|---|---|---|
| Zeros.Ordered | γ_n nondecreasing (n ≥ 1), γ_1 > 0, γ_k ≤ T ⟺ k ≤ N(T) | l.176–177: 0 < γ_1 ≤ γ_2 ≤ ⋯ with multiplicity, N(T) = Σ_{0<γ_n≤T} 1 |
| ScaleToInf | ∀ K, eventually K ≤ sc T | sc T = log T/(2π) |
| RvM | ∀ ε > 0, eventually T·sc T ≤ (1 + ε) N(T) | l.177 (Riemann–von Mangoldt); DERIVATION: x_{N_0−1} − x_1 ≤ T log T/2π = N(1 + O(1/log T)) |
| Steps12 | ∀ θ ∈ (1/2, 4/7), λ ∈ (0, 1], F ∈ C_c² even ≥ 0, ε > 0 ∃ δ₀ > 0 ∀ δ ∈ (0, δ₀) ∃ T₀ ∀ T ≥ T₀: #E_λ(T) < δN(T) ⇒ (J_θ(F) − ε)N(T) ≤ Σ_{a,b∈Z_0(T)} F^((γ_a − γ_b) log T/2π)·4/(4 + (γ_a − γ_b)²) | (Eq:ReductionGoodSet2) l.437–442, from Prop SmallSetZeros l.331–340 for g = ψ_εF (l.413–429) and the positivity identity l.430–436; on-line summand l.443–446 (w(u) = 4/(4 − u²), l.292); error O(δ^{1/16−7θ/64}) with T large → ε–δ–T₀ |
| FourierCc2 | F ∈ C_c² even ⇒ F^ even and ∃ M, |F^(u)|u² ≤ M | two integrations by parts, F^(u) = −(F'')^(u)/(4π²u²) |
| LemmaT′ | Schwartz f, even, ≥ 0, f^ ≤ 0 on [U, ∞), (H), J_θ(f) > f^(0) + 2c, 0 < λ ≤ U, θ ≥ 0 ⇒ ∃ F ∈ C_c², even, ≥ 0, F^ ≤ 0 on [U, ∞), (H) with the same A, β, c, J_θ(F) > F^(0) + 2c | DERIVATION §4 Lemma T′ (F = f h_R + c_r G; Lamzouri l.418–440 for its admissibility) |
| CertP | ∀ w ≥ 0, p(w) > 0 (p = pcoef) | cert_aug.py l.43–46: p(0) > 0, exact Sturm count 0 on (0, ∞) |
| CertQ | ∀ w ≥ w0 = 2·(314159265/10⁸)·s²U², q(w) < 0 (q = qcoef) | l.47–52: q(w0) < 0, Sturm count 0 on (w0, ∞) |
| HermRep | f(α) = p(κα²)G(α), f^(u) = (11/10)·q(κ′u²)G′(u); κ ≥ 0, G, G′ > 0, κ′ ≥ 2·3.14159265·s², G′(0) = 1 | DERIVATION §5 (κ = 2π/s², G = e^{−πα²/s²}, κ′ = 2πs², G′ = e^{−πs²u²}: FT of H_{2n}(√(2π)α/s)e^{−πα²/s²} is s(−1)^n H_{2n}(√(2π)su)e^{−πs²u²}) |
| CertJ | J_{θ0}(f) ≥ 153419640989/(5·10¹¹) = 0.306839281978, θ0 = 399993/700000 | l.53–66 (Arb, 1200 bits): J = 0.306839281979 ± 3.75·10⁻¹³ |
| CertR1 | Φ(x, y) ≤ c0 = −347861653223/10¹² on [λ, U]² | l.98–105, 126–127 (corners max −0.3478888283 + 3M₂δ²/8 = 2.718·10⁻⁵; M₂ l.77–96); c_cert = −0.3478616532233634 (binary64, l.130) ≤ c0 |
| CertR2 | f^(x) + A(x) − A(U) + β(1 − x) ≤ c0 on [λ, U] | l.128 (R2 = −0.5640124229) |
| Sch f, f even | Schwartz, even | Hermite polynomials of even degree × Gaussian |

Exact constants (Lean/audit/cert_exact.out, a copy of cert_aug.py with exact prints; same PASS and margin): β = 416188506333/640725853462, f^(0) = s·q(0) = 15275573730468744137181893367080725968011933/15258789062500000000000000000000000000000000 < 1.0011, min_k A_k = A_0 = −378347742845/758342134484, last node 0, R3 = −0.3481091233 ≤ c0. Margin in Lean: 0.306839281978 − 1.0011 + 2·0.347861653223 = 0.001462588424 > 0 (norm_num). pcoef, qcoef, nodesA equal the checker's python-flint p, q and the JSON's A exactly (Lean/audit/xcoef.py: all True).

## Build

Each module separately, one at a time, lean directly (no lake process), oleans to Lean/audit/(olean build output, deleted):

    lake env printenv LEAN_PATH > Lean/audit/leanpath.txt     # once
    bash Lean/audit/verify.sh    # for m in FinSum Model NoGood Theorem1 PL Herm Cert Liminf Final Check:
      # LEAN_PATH=Lean/audit/(olean build output, deleted):<lake path> bash Lean/audit/memcap.sh 290 1800 /usr/bin/time -l lean -j1 -o … $m.lean

All 10: rc = 0, 0 errors, 0 warnings; peak memory footprint (/usr/bin/time -l) 245–274 MiB (table above), 8–17 s each (Lean/audit/verify.log, final run 12:58 UTC). Kill-free; the memcap readings were 235–245 MB (0 = finished between samples). An earlier single-file Theorem1 reached 287–292 MiB and was split (NoGood). By lake instead: `lake build Riemann.Deep.ExplorePointProc.Lean.Final` (not run: lake's own process adds to the cap). No sorry, admit, axiom, native_decide, decide, skipKernelTC (rg). Check.lean:

    'PPLean.pairSum_le' depends on axioms: [propext, Classical.choice, Quot.sound]
    (same line for pairSum_weighted_le, card_Elam_le, DPos.mono, pairSum_Z0, theorem1_Cc2, theorem1, thmMain2, HH_div,
     certHH, cert_margin, dpos_768, dpos_ge_768, Acert_tail, Acert_range, pos_of_rep, neg_of_rep, fh0_of_rep,
     dpos_ge_768_exact, dpos_iff, nToInf, theorem1_D, D_pos_ge_768: 23 theorems, all [propext, Classical.choice, Quot.sound])

## Proof route

- (iii′) pairSum_le: induction on n, adding the point x_{n}: its row Σ_{i<n} F(x_n − x_i) keeps only the two nearest left neighbours (gap_sum: x_i + dλ ≤ x_{i+d}, d ≥ 3 ⇒ ≥ 3λ ≥ U, F ≤ 0); (H) at (g_{n−2}, g_{n−1}) absorbs 2F(g_{n−2}) + 2F(g_{n−2} + g_{n−1}); A telescopes, β collects x_{n−1} − x_0 − (n − 1). Exact: equality for n = 2.
- Weighted: G = F·w is even, ≤ 0 on [U, ∞) (w ≥ 0), G(0) = F(0), and satisfies (H) with c + 2η.
- Theorem 1 (no_good_height, theorem1_Cc2): μ = J_θ(F) − F^(0) − 2c > 0; Steps12 at ε = μ/6 gives δ₀; δ ≤ min(δ₀/2, 1/2, μ/(6(|F^(0)+2c| + 2β + 1))); for T ≥ T₀: if #E_λ < δN, the pair sum over Z_0 equals pairSum G xs n (pairSum_Z0, w(u) = 4s²/(4s² + u²) = 4/(4 + (γ−γ′)²), s = log T/2π); (1 − δ)N ≤ n ≤ N (card_Z0), gaps > λ (Z0_sep), x_{n−2} − x_0 ≤ T s ≤ (1 + ε₂)N (RvM), η = |M|/(4s²) ≤ μ/24 (s → ∞), boundary terms ≤ 4|K_A| + 2|M|/λ² + 4|c| + 4β ≤ μN/6: core_arith gives μN ≤ 0, false. So #E_λ ≥ δN eventually, #small gaps ≥ δN/2 (card_Elam_le, l.407–411): DPos.
- Corollary: certHH covers [λ, ∞)² by R1, R2 and R3 (x ≥ U: A(x) = 0, f^ < 0 twice, β(1 − x) ≤ β(1 − U), −A(y) ≤ −A_0; norm_num), margin by norm_num, theorem1 at (θ0, λ, U = 3λ, β, c0), DPos.mono; dpos_iff (N(T) → ∞ from ScaleToInf + RvM) turns DPos into 0 < liminf.

## Fidelity (against DERIVATION §4–§5 and Lamzouri's TeX)

| Lean | source | match |
|---|---|---|
| HH | (H), §4 | exact |
| pairSum_le | (iii′), §4; brief's form | exact boundary form, 0-indexed; gaps ≥ λ (weaker than > λ), λ ≥ 0, n ≥ 2, no condition on β, A, c (more general) |
| pairSum_weighted_le | §4 "w = 1 − O(u²/log²T), sup |f^|u² < ∞, replacing w by 1 costs o(N)"; "index distance ≥ 3: w > 0, drop" | exact (the o(N) is 4nη) |
| Zeros, Ordered, N | l.176–177 | exact; ordinates only (on-line, simple: l.443–446 is inside Steps12) |
| smallGaps, DPos, Dlam | l.183 D(λ) = liminf_{T→∞} N(T)⁻¹ #{0 < γ_n ≤ T : γ_{n+1} − γ_n ≤ 2πλ/log T} | exact (unfolded: (γ_{n+1} − γ_n)·log T/2π ≤ λ); dpos_iff proves DPos ⟺ 0 < Dlam |
| Elam, Z0 | (Eq:DefMultiSetE) l.404–406; Z_0 l.411 | exact (n ≥ 2, min of both gaps, multiplicity = indices) |
| card_Elam_le | l.407–410 (#E ≥ δN ⇒ #small gaps ≥ δN/2 + O(1)) | proved, without the O(1) |
| Z0_sep | l.443–444 "|γ − γ′| > 2πλ/log T" | proved |
| Steps12 | l.425–446 | implied by Lamzouri: his O(δ^{1/16−7θ/64}) (T large, constant depending on f, θ) gives the ε–δ–T₀ form; for Thm Main2's class (no support condition, l.392, 413–440) |
| FourierCc2, LemmaT′ | standard; DERIVATION §4 Lemma T′ | true for the intended objects (reviewer); LemmaT′ for every θ ≥ 0 (cut-off radius R ≥ θ) |
| ScaleToInf, RvM | log T → ∞; RvM l.177 | true |
| theorem1, theorem1_D | Theorem 1, §4 | exact hypotheses (1/2 < θ < 4/7, 0 < λ ≤ 1, λ < U ≤ 3λ, f even ≥ 0 Schwartz, f^ ≤ 0 on [U, ∞), A bounded on [λ, ∞), β ≥ 0, (H), J > f^(0) + 2c); conclusion D(λ) > 0 |
| theorem1_Cc2 | the C_c² core | needs only U ≤ 3λ |
| thmMain2 | Thm Main2 l.392–400 | A = β = c = 0, U = λ; f^ ≤ 0 on [λ, ∞) (Lamzouri: |u| > λ; equal by continuity of f^); without f^(0) = 1 (J > f^(0) is the homogeneous form) |
| HH_div, margin_div | §4 Normalisation | divide f, A, β, c by f^(0) > 0; not needed by theorem1 (it never normalises) |
| CertP, CertQ, HermRep, CertJ, CertR1, CertR2 | §5 items 1–4, cert_aug.py l.43–130 | one named hypothesis per checked fact; directions: J rounded down, f^(0) and c rounded up |
| Acert, nodesA | §5 "A piecewise linear on 78 rational nodes λ + k·dg, A(x) := A(U) for x ≥ U" | exact (plin, dg = 192/9625) |
| certHH | §5 item 4, R1–R3 | R1, R2 assumed, R3 proved |
| D_pos_ge_768 | Corollary §5 | exact; monotonicity proved |

Reviewer (verdict-only, statement fidelity; DERIVATION §4–§5, Lamzouri TeX l.176–446 and all Lean statements pasted inline; Lean/audit/review_out.md): FAITHFUL. Confirmed: (1) pairSum_le is the endpoint-expanded (iii′) (its extra −4c + 4β + 2(A(g_last) − A(g_first)) + 2F(g_last) is DERIVATION's O(1)); (2) Steps12 follows from (Eq:ReductionGoodSet2); (3) smallGaps/Dlam match l.183, DPos ⟺ positivity as N(T) → ∞; (4) Elam/Z0 match l.404–411; (5) the ordinate model is faithful; (6) FourierCc2 and LemmaT′ are true; theorem1_D matches §4, D_pos_ge_768 the corollary; certificate predicates are valid consequences of the checker facts. One [minor]: this is a conditional formalization (named hypotheses), not a Lean construction of ζ or Fourier analysis, nor a Lean check of the certificate.

## Not done

- The named hypotheses themselves: Lamzouri's steps (i)–(ii) (Steps12), Lemma T′, the Fourier facts, Riemann–von Mangoldt, and the checker's Sturm/Arb facts are not proved in Lean (cap: Mathlib's zeta, Fourier and Real.log imports exceed 290 MB; a Lean Sturm/Arb check of degree-21 polynomials and 386² lattice corners was out of the time box).
- .lake: left as it is (real directory, another agent's Mathlib cache; not modified: no file under .lake newer than this run's start, no lake build). My oleans are in Lean/audit/(olean build output, deleted). The orchestrator restores the placeholder.

## FORMALIZATIONS.md entry

- Hard-core augmentation of Lamzouri's small-gap criterion (our discovery [small-gap-hardcore](discoveries/small-gap-hardcore.md)) in Lean 4 + Mathlib: step (iii′) proved outright (PPLean.pairSum_le); Theorem 1 (0 < D(λ) for Lamzouri's liminf D, PPLean.theorem1_D) and D(λ) > 0 for every λ ≥ 0.768 (PPLean.D_pos_ge_768) from named inputs (Lamzouri's (Eq:ReductionGoodSet2), Lemma T′, Fourier facts, Riemann–von Mangoldt; the checker's Sturm and Arb facts), with the margin, R3, the piecewise-linear A, monotonicity in λ and DPos ⟺ liminf > 0 proved; sorry-free, standard axioms. [LEAN.md](Riemann/Deep/ExplorePointProc/LEAN.md)
