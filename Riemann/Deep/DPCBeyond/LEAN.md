FORMALIZED: Theorem U′ (DPCBLean.theoremUB, theoremUB_fine: for κ > 0 and every φ ∈ (max{2/(2+κ), 1/2}, 1) there is c > 0 with 𝒩(T) ≥ c·T·log T/(2π) for all large T) and D(λ) > 0 for every λ > 7/11 (DPCBLean.D_pos_seven_elevenths(_fine), Lamzouri's liminf D; theoremUB_conrey: 𝒩(T) ≫ T log T for every φ ∈ (7/11, 1)), from named analytic inputs: Theorem U's DPC construction with Prop 1 and Lemma 6 (MainTerm), the slit identity (Step3), Step 1's positivity (Step1Pos) and the covering bound (Covering), reused verbatim from dpclean; the uniform density slope κ (DensityK; Conrey κ = 8/7 − δ for every δ ∈ (0, 8/7)); 6e′'s class inequalities (Step6ClassesB: 6a's decomposition, Lemma Z with the y-integral, Hölder over S_a, M4 and 6b–6d's moments); Step 7′'s Hölder (HolderB); Lemma M3 (MomentM, or in the fine version M3Inputs: Lemma Z + Hölder, M1/Fujii, Minkowski, the depth bins of Off, M2); 5c (Step5cB); Riemann–von Mangoldt's upper bound (RvMUp); pplean's zero model (Zeros.Ordered, Reflect); four exp/log facts (ExpLog). Lemma Z itself is a theorem from its named inputs (i)–(iii) in both cases (lemmaZ_point1, lemmaZ_point2: |R_T(t − iy)| ≤ K(1 + W_{−y})(1 + On)^{k−1}, |R*(t − iy)| ≤ KX^y(1 + W)(1 + On)^{k−1}) and is wired into M3's first inequality (M3InputsW ⇒ M3Inputs: theoremUB_wired, D_pos_seven_elevenths_wired) and into 6e′ (a) pointwise (hairpin_pointwise); Step6ClassesB (a)–(b) stays named (see Not done). Proved, not assumed: Steps 0′/8′ and the order of choices (θ ∈ (1/φ − 1, min{κ/2, 1}) exists iff φ > max{2/(2+κ), 1/2}, then 0 < φθ < 1 < φ(1+θ) < 2 and q, s, u, w, v; κ = 8/7 − δ gives exactly φ > 7/11), Lemma Z's mechanism (∏e^{(Re z_i)u_i} ≤ (X/𝔞)^{max (Re z_i)⁺} on Σu_i ≤ log(X/𝔞)) and summation (Σ_{ρ_1..ρ_m} X^{max} ∏κ_A ≤ m·W·On^{m−1}, the packet exponents, the assembly over terms from named per-term bounds), Lemma M's bookkeeping (M2 ≪ Te^{−κa(1−1/w)}, M3's bin series and choice of w, M4's exponent κ(1−1/w)² with (|m|+2), the kernel sum ≤ 4 for A ≥ 3, M3's assembly ∫|R_T|^q ≪ T), 6e′'s class bookkeeping and Φ(ε) → 0, Step 7′ (meas E_T ≥ c₂T, c₂ independent of ε, then Step 7's contradiction), Step 2 (dpclean's #off-line ≤ 2𝒩(T)) and Theorem U′ ⇒ D(λ) > 0. Sorry-free; axioms propext, Classical.choice, Quot.sound only; no native_decide, no decide, no skipKernelTC.

# DPCBeyond (Theorem U′) in Lean 4 + Mathlib (dpcblean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile). Riemann/Deep/DPCBeyond/Lean/, 17 files, 2445 lines, namespace DPCBLean; not imported by Riemann.lean. Imports dpclean's Theorem U modules unchanged (Riemann/Deep/DPCUncond/Lean: Tail's ExpLog, Phi, Hairpin's hairpin_sum, Step7's step7_bad, Gaps' gapSet/offSet/Reflect/card_off_le, TheoremU's Data/Step3/MainTerm/Step1Pos/Covering, DLam's RvMUp/ScaleHalf/scaleToInf_of/scaleHalf_of/gapSet_half_sub, Fine's deepSet/classSet) and through them pplean's zero model and liminf D (Riemann/Deep/ExplorePointProc/Lean). No Mathlib zeta, Fourier, LSeries, complex analysis, measure theory, Real.log, Real.rpow or π (290 MB cap): integrals and norms are abstract real functions of T, exp/log abstract (ExpLog), x^e = rp ex lg x e, every analytic input a named hypothesis.

| file | content | lines | footprint (time -l) | memcap reading | time |
|---|---|---|---|---|---|
| ParamsB | (a) one_lt_of, step0B_nonempty, step0B_signs, step8B_empty, step0B_iff, exists_q_s, exists_u_w, order_of_choices, conrey_threshold, conrey_sharp | 154 | 274 MiB | 262 MB | 23 s |
| Rpow | rp (x^e) and lgex, lg_lt, lg_le, lg_nonneg, lg_nonpos, rp_of_pos, rp_of_nonpos, rp_pos, rp_nonneg, rp_one, rp_zero, rp_add, rp_mul, rp_rp, rp_neg, rp_le_rp, rp_lt_rp, rp_le_rp_exp, rp_le_max | 120 | 254 MiB | 239 MB | 31 s |
| Step7B | (e) step7B_meas, step7B | 181 | 277 MiB | 254 MB | 26 s |
| HairpinB | (d) DensityK, Step6ClassesB, class_count, class_bound | 179 | 282 MiB | 242 MB | 40 s |
| SumB | (d) step6B_of_classes, PhiB_small | 120 | 272 MiB | 245 MB | 31 s |
| DLamB | (f) UFrom, dpos_of_from, nToInf_of_from, dlam_pos_of_from | 106 | 273 MiB | 248 MB | 31 s |
| TheoremUB | TheoremUB, DataB, DataB.base, HolderB, MomentM, Step5cB, InputsB, DPCInputsB, cover_of, moments_of, split_c, theoremUB | 165 | 265 MiB | 245 MB | 31 s |
| FinalB | (f) D_pos_of_kappa, theoremUB_conrey, D_pos_seven_elevenths | 53 | 255 MiB | 246 MB | 41 s |
| LemmaZ | (b) ex_sum, sum_mul_le, prod_ex_le, prod_ex_le_sup, rp_sup_le_sum, mul_prod_eq, tuple_one, tuple_sum, tuple_bound, packet_case1, packet_case2, packet_nofree, pos_part_shift, pos_part_split, term_bound, mono_bound, lemmaZ_sum | 262 | 272 MiB | 256 MB | 33 s |
| LemmaM | (c) rp_ex, M2_bound, geom_mul, geom_le, M3_series, M3_exists_w | 115 | 262 MiB | 242 MB | 8 s |
| LemmaM4 | (c) M4_count, M4_moment, key_ineq, inv_sq_sum, kernel_sum | 173 | 278 MiB | 245 MB | 15 s |
| LemmaM3 | (c) off_bound, momentM_of, exists_r_w, M3Inputs, momentM_of_inputs | 169 | 274 MiB | 242 MB | 43 s |
| FineB | InputsFineB, DPCInputsFineB, inputsB_of_fine, theoremUB_fine, D_pos_seven_elevenths_fine | 57 | 252 MiB | 206 MB | 8 s |
| LemmaZ2 | (b) Wsum, Onsum, Offsum, W_split, tuple_case2, term_bound2, term_nofree1, term_nofree2, ZPoint1, ZPoint2, final_shape, lemmaZ_point1, lemmaZ_point2 | 263 | 276 MiB | 252 MB | 14 s |
| LemmaZW | (b) wiring: m3_first_of_lemmaZ, hairpin_pointwise | 134 | 277 MiB | 210 MB | 16 s |
| CheckB | #print axioms (72 theorems) | 79 | 250 MiB | 242 MB | 17 s |
| WiredB | M3InputsW, m3Inputs_of_wired, InputsWB, DPCInputsWB, inputsFineB_of_wired, theoremUB_wired, D_pos_seven_elevenths_wired (+ 4 #print axioms) | 115 | 262 MiB | 240 MB | 47 s |

Imported modules in the same final run (footprint MiB): FinSum 255, Model 265, Liminf 250, Params 268, Gaps 255, Tail 278, Phi 267, Step7 273, TheoremU 265, DLam 279, Hairpin 271, Fine 269. Final run: all 29 compiles rc = 0, peak footprint 250–282 MiB, no watchdog kill.

## Main statements (verbatim; proofs omitted)

    -- pplean: Z.γ n (n ≥ 1, with multiplicity), Z.N T, Z.sc T = log T/(2π); dpclean: gapSet Z φ T = 𝒩(T), Data, Step3,
    -- MainTerm, Step1Pos, Covering, Reflect, offSet, classSet, deepSet, ExpLog, Phi (Riemann/Deep/DPCUncond/LEAN.md)
    def UFrom (Z : Zeros) (φ₀ : ℝ) : Prop :=
      ∀ φ : ℝ, φ₀ < φ → φ < 1 → ∃ c : ℝ, 0 < c ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
        c * T * Z.sc T ≤ ((gapSet Z φ T).card : ℝ)
    def TheoremUB (Z : Zeros) (κ : ℝ) : Prop := UFrom Z (max (2 / (2 + κ)) (1 / 2))

    /-- The analytic quantities at fixed Step-0 data (φ, θ, P, f, η; PROOF Step 0), as functions of T: Theorem U's
    Data with ∫|R_T|⁴, ∫|P(V)|⁴ replaced by the moment families Rq q T = ∫_T^{2T}|R_T|^q, Ps s T = ∫_T^{2T}|P(V_φ)|^s. -/
    theorem theoremUB {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        {κ : ℝ} (hκ : 0 < κ) (hden : DensityK Z β ex lg κ) (hin : DPCInputsB Z β ex lg κ) :
        TheoremUB Z κ
    theorem theoremUB_fine {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        {κ : ℝ} (hκ : 0 < κ) (hden : DensityK Z β ex lg κ) (hin : DPCInputsFineB Z β ex lg κ) :
        TheoremUB Z κ :=
    theorem theoremUB_conrey {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        (hden : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DensityK Z β ex lg (8 / 7 - δ))
        (hin : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DPCInputsB Z β ex lg (8 / 7 - δ)) : UFrom Z (7 / 11)
    theorem D_pos_of_kappa {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        (hR : RvMUp Z) {κ : ℝ} (hκ : 0 < κ) (hden : DensityK Z β ex lg κ)
        (hin : DPCInputsB Z β ex lg κ) :
        ∀ lam : ℝ, max (2 / (2 + κ)) (1 / 2) < lam → 0 < Dlam Z lam
    theorem D_pos_seven_elevenths {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        (hR : RvMUp Z) (hden : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DensityK Z β ex lg (8 / 7 - δ))
        (hin : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DPCInputsB Z β ex lg (8 / 7 - δ)) :
        ∀ lam : ℝ, 7 / 11 < lam → 0 < Dlam Z lam
    theorem D_pos_seven_elevenths_fine {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        (hR : RvMUp Z) (hden : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DensityK Z β ex lg (8 / 7 - δ))
        (hin : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DPCInputsFineB Z β ex lg (8 / 7 - δ)) :
        ∀ lam : ℝ, 7 / 11 < lam → 0 < Dlam Z lam :=
    structure DataB where
      /-- Re ∫_T^{2T} G_T|R_T|² P(V_φ) -/
      I : ℝ → ℝ
      /-- Re M_T -/
      M : ℝ → ℝ
      /-- Re 𝓗_T, the hairpin total -/
      Hp : ℝ → ℝ
      /-- Σ_n |c_n|² -/
      S : ℝ → ℝ
      /-- Rq q T = ∫_T^{2T} |R_T|^q (Lemma M3) -/
      Rq : ℝ → ℝ → ℝ
      /-- Ps s T = ∫_T^{2T} |P(V_φ)|^s (5c) -/
      Ps : ℝ → ℝ → ℝ
      /-- ∫_{E_T} G_T |R_T|² |P(V_φ)| -/
      Bad : ℝ → ℝ
      /-- meas E_T -/
      mE : ℝ → ℝ

    /-- Theorem U's Data of the same construction (R4 = Rq 4, P4 = Ps 4), for reusing Step3, MainTerm, Step1Pos,
    Covering verbatim. -/
    def DataB.base (D : DataB) : Data := ⟨D.I, D.M, D.Hp, D.S, D.Rq 4, D.Ps 4, D.Bad, D.mE⟩

    /-- Step 7′, Hölder with exponents (q/2, s, 1/α), α = 1 − 2/q − 1/s, G_T bounded:
    ∫_{E_T} G_T|R_T|²|P(V)| ≤ K_G ‖R_T‖_q² ‖P(V)‖_s meas(E_T)^α. -/
    def HolderB (ex lg : ℝ → ℝ) (D : DataB) : Prop :=
      ∀ q s : ℝ, 2 < q → 1 < s → 2 / q + 1 / s < 1 → ∃ KG T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
        D.Bad T ≤ KG * rp ex lg (D.Rq q T) (2 / q) * rp ex lg (D.Ps s T) (1 / s) *
          rp ex lg (D.mE T) (1 - 2 / q - 1 / s)
    def MomentM (κ θ : ℝ) (D : DataB) : Prop :=
      ∀ q : ℝ, 2 < q → q < κ / θ → ∃ K₁ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → D.Rq q T ≤ K₁ * T
    def Step5cB (D : DataB) : Prop :=
      ∀ s : ℝ, 1 < s → ∃ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → D.Ps s T ≤ K₂ * T
    structure InputsB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ φ θ : ℝ) (D : DataB) : Prop where
      step3 : Step3 D.base
      main : MainTerm D.base
      step6 : Step6ClassesB Z β ex lg κ θ D.Hp
      step1 : Step1Pos D.base
      holder : HolderB ex lg D
      momM : MomentM κ θ D
      step5c : Step5cB D
      cover : Covering Z φ D.base
    def DPCInputsB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ : ℝ) : Prop :=
      ∀ φ θ : ℝ, 0 < φ → φ < 1 → 0 < θ → θ < min (κ / 2) 1 → 0 < φ * θ → φ * θ < 1 →
        1 < φ * (1 + θ) → φ * (1 + θ) < 2 → ∃ D : DataB, InputsB Z β ex lg κ φ θ D
    def DensityK (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ : ℝ) : Prop :=
      ∃ Cd T₀ : ℝ, 0 < Cd ∧ ∀ T : ℝ, T₀ ≤ T → ∀ a : ℕ,
        ((deepSet Z β lg T a).card : ℝ) ≤ Cd * T * lg T * ex (-(κ * a))
    def Step6ClassesB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ θ : ℝ) (Hp : ℝ → ℝ) : Prop :=
      ∀ u w : ℝ, 1 < u → 1 < w → 0 < 1 - 1 / u - 1 / w → 2 * u * θ < κ → ∃ (C : ℕ) (K : ℝ), 0 ≤ K ∧
        ∃ (Ha Y₁ Y₂ MY : ℝ → ℕ → ℝ) (MO : ℝ → ℝ) (Mt : ℝ → ℕ), ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ,
          T₀ ≤ T → |Hp T| ≤ ∑ a ∈ range (Mt T + 1), Ha T a + η * T ∧ MO T ≤ K * (T * lg T) ∧
          ∀ a : ℕ, Ha T a ≤ K * ((a : ℝ) + 2) ^ C / lg T * (ex (2 * θ) ^ (a + 1) * Y₁ T a + Y₂ T a) ∧
            MY T a ≤ K * ((a : ℝ) + 2) ^ C * (T * lg T) ∧
            Y₁ T a ≤ rp ex lg ((classSet Z β lg T a).card : ℝ) (1 - 1 / w) * rp ex lg (MY T a) (1 / w) ∧
            Y₂ T a ≤ rp ex lg ((classSet Z β lg T a).card : ℝ) (1 - 1 / u - 1 / w) *
              rp ex lg (MO T) (1 / u) * rp ex lg (MY T a) (1 / w)
    def M3Inputs (ex lg : ℝ → ℝ) (κ θ : ℝ) (D : DataB) : Prop :=
      ∀ q r w : ℝ, 2 < q → q < κ / θ → 1 < r → 1 < w → q * r * θ < κ * (1 - 1 / w) →
        ∃ (KZ K₀ K₁ T₀ : ℝ) (IW IO NOn NOff : ℝ → ℝ) (NF : ℝ → ℕ → ℝ), 0 ≤ KZ ∧ 0 ≤ K₀ ∧ 0 ≤ K₁ ∧
          (∀ T : ℝ, T₀ ≤ T → D.Rq q T ≤ KZ * rp ex lg (IW T) (1 / r) * rp ex lg (IO T) (1 - 1 / r)) ∧
          (∀ T : ℝ, T₀ ≤ T → IO T ≤ K₀ * T) ∧
          (∀ T : ℝ, T₀ ≤ T →
            rp ex lg (IW T) (1 / (q * r)) ≤ rp ex lg T (1 / (q * r)) + NOn T + NOff T) ∧
          (∀ T : ℝ, T₀ ≤ T → NOn T ≤ K₀ * rp ex lg T (1 / (q * r))) ∧
          (∀ T : ℝ, T₀ ≤ T → ∃ M : ℕ,
            NOff T ≤ ∑ a ∈ range M, ex θ ^ (a + 1) * NF T a + rp ex lg T (1 / (q * r))) ∧
          (∀ T : ℝ, T₀ ≤ T → ∀ a : ℕ,
            NF T a ≤ K₁ * rp ex lg T (1 / (q * r)) * ex (-(κ * a * (1 - 1 / w) / (q * r))))
    structure InputsFineB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ φ θ : ℝ) (D : DataB) : Prop where
      step3 : Step3 D.base
      main : MainTerm D.base
      step6 : Step6ClassesB Z β ex lg κ θ D.Hp
      step1 : Step1Pos D.base
      holder : HolderB ex lg D
      m3 : M3Inputs ex lg κ θ D
      step5c : Step5cB D
      cover : Covering Z φ D.base
    theorem step0B_iff {κ φ : ℝ} (hκ : 0 < κ) (h0 : 0 < φ) :
        (∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < min (κ / 2) 1) ↔ max (2 / (2 + κ)) (1 / 2) < φ
    theorem step0B_signs {κ φ θ : ℝ} (h1 : max (2 / (2 + κ)) (1 / 2) < φ) (h2 : φ < 1)
        (h3 : 1 / φ - 1 < θ) (h4 : θ < min (κ / 2) 1) :
        0 < θ ∧ θ < κ / 2 ∧ θ < 1 ∧ 0 < φ * θ ∧ φ * θ < 1 ∧ 1 < φ * (1 + θ) ∧ φ * (1 + θ) < 2
    theorem order_of_choices {κ φ : ℝ} (hκ : 0 < κ) (h1 : max (2 / (2 + κ)) (1 / 2) < φ) (h2 : φ < 1) :
        ∃ θ q s u w : ℝ, 1 / φ - 1 < θ ∧ θ < min (κ / 2) 1 ∧ 0 < φ * θ ∧ φ * θ < 1 ∧
          1 < φ * (1 + θ) ∧ φ * (1 + θ) < 2 ∧ 2 < q ∧ q < κ / θ ∧ 1 < s ∧ 2 / q + 1 / s < 1 ∧
          1 < u ∧ u < κ / (2 * θ) ∧ u / (u - 1) < w ∧ 2 * θ < κ * (1 - 1 / w) ∧ 0 < 1 - 1 / u - 1 / w
    theorem conrey_threshold {φ : ℝ} (h1 : 7 / 11 < φ) :
        ∃ δ : ℝ, 0 < δ ∧ δ < 8 / 7 ∧ max (2 / (2 + (8 / 7 - δ))) (1 / 2) < φ
    theorem conrey_sharp {φ δ : ℝ} (h1 : φ ≤ 7 / 11) (hδ0 : 0 ≤ δ) (hδ1 : δ < 8 / 7) :
        ¬ max (2 / (2 + (8 / 7 - δ))) (1 / 2) < φ
    theorem prod_ex_le_sup (h : ExpLog ex lg) {ι : Type*} (s : Finset ι) (hs : s.Nonempty)
        (d u : ι → ℝ) {Y : ℝ} (hY : 0 < Y) (hu : ∀ i ∈ s, 0 ≤ u i) (hl : ∑ i ∈ s, u i ≤ lg Y) :
        ∏ i ∈ s, ex (d i * u i) ≤ rp ex lg Y (s.sup' hs (fun i => max (d i) 0))
    lemma rp_sup_le_sum (h : ExpLog ex lg) {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (c : ι → ℝ)
        (X : ℝ) : rp ex lg X (s.sup' hs c) ≤ ∑ i ∈ s, rp ex lg X (c i)
    theorem tuple_bound (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] (S : Finset ι) (c k : ι → ℝ)
        (hk : ∀ x, 0 ≤ k x) (X : ℝ) (n : ℕ) :
        ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
            rp ex lg X (univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i) ≤
          (n + 1) * (∑ x ∈ S, rp ex lg X (c x) * k x) * (∑ x ∈ S, k x) ^ n
    lemma packet_case2 (h : ExpLog ex lg) {a X y D : ℝ} (ha : 1 ≤ a) (hX : 0 < X) (hD : 0 ≤ D) :
        rp ex lg a (-(1 / 2) + y) * rp ex lg (X / a) (y + D) ≤
          rp ex lg X (y + D) * rp ex lg a (-(1 / 2))
    theorem term_bound (h : ExpLog ex lg) {ι α : Type*} [DecidableEq ι] (S : Finset ι) (c k : ι → ℝ)
        (hc : ∀ x, 0 ≤ c x) (hk : ∀ x, 0 ≤ k x) {X y CI Kpk : ℝ} (hX : 0 < X) (hy : 0 ≤ y)
        (hCI : 0 ≤ CI) (n : ℕ) (P : Finset α) (v C F : α → ℝ) (hv : ∀ p ∈ P, 1 ≤ v p)
        (hC : ∀ p ∈ P, 0 ≤ C p)
        (hF : ∀ p ∈ P, |F p| ≤ CI * ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
          rp ex lg (X / v p) (univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i))
        (hP : ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) ≤ Kpk) :
        ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) - y) * |F p| ≤
          Kpk * CI * ((n + 1) * (∑ x ∈ S, rp ex lg X (c x) * k x) * (∑ x ∈ S, k x) ^ n)
    theorem lemmaZ_sum {γ : Type*} (𝒯 : Finset γ) (Bt : γ → ℝ) (nt : γ → ℕ) {Rabs E W On Kpk CI : ℝ}
        {k : ℕ} (hk : 1 ≤ k) (hW : 0 ≤ W) (hOn : 0 ≤ On) (hK : 0 ≤ Kpk * CI)
        (hR : Rabs ≤ ∑ τ ∈ 𝒯, Bt τ + E)
        (hB : ∀ τ ∈ 𝒯, nt τ + 1 ≤ k ∧
          (Bt τ ≤ Kpk * CI * (((nt τ : ℝ) + 1) * W * On ^ nt τ) ∨ Bt τ ≤ Kpk * CI)) :
        Rabs ≤ 𝒯.card * (Kpk * CI * (k * (1 + W) * (1 + On) ^ (k - 1))) + E
    theorem M2_bound (h : ExpLog ex lg) {I μ B N L C T κ a w : ℝ} (hL : 0 < L) (hC : 0 < C)
        (hT : 0 < T) (hw : 1 < w) (hI : I ≤ rp ex lg μ (1 - 1 / w) * rp ex lg N (1 / w))
        (hμ : μ ≤ B / L) (hB : B ≤ C * T * L * ex (-(κ * a))) (hN : N ≤ C * T) :
        I ≤ C * T * ex (-(κ * a * (1 - 1 / w)))
    theorem M3_series (h : ExpLog ex lg) {θ κ w p : ℝ} (hp : 0 < p) (hpθ : p * θ < κ * (1 - 1 / w))
        (M : ℕ) : ∑ a ∈ range M, ex (θ * (a + 1) - κ * a * (1 - 1 / w) / p) ≤
          ex θ / (1 - ex (θ - κ * (1 - 1 / w) / p))
    theorem M3_exists_w {θ κ p : ℝ} (hκ : 0 < κ) (hpθ0 : 0 ≤ p * θ) (hpθ : p * θ < κ) :
        ∃ w : ℝ, 1 < w ∧ p * θ < κ * (1 - 1 / w) ∧ p * θ < κ * (1 - 1 / w) ^ 2
    theorem M4_count (h : ExpLog ex lg) {Nc B Sw C T L κ a w mm : ℝ} (hC : 0 < C) (hT : 0 < T)
        (hL : 0 < L) (hw : 1 < w) (hm : 0 ≤ mm) (hN : Nc ≤ rp ex lg B (1 - 1 / w) * rp ex lg Sw (1 / w))
        (hB : B ≤ C * (T * L) * ex (-(κ * a))) (hS : Sw ≤ C * rp ex lg (mm + 2) w * (T * L)) :
        Nc ≤ C * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w)))
    theorem M4_moment (h : ExpLog ex lg) {I Nc Sp C C' T L κ a w mm : ℝ} (hC : 0 ≤ C) (hC' : 0 ≤ C')
        (hT : 0 < T) (hL : 0 < L) (hw : 1 < w) (hm : 0 ≤ mm)
        (hI : I ≤ rp ex lg Nc (1 - 1 / w) * rp ex lg Sp (1 / w))
        (hN : Nc ≤ C * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w)))) (hS : Sp ≤ C' * (T * L)) :
        I ≤ (C + 1) * (C' + 1) * (mm + 2) * (T * L) * ex (-(κ * a * (1 - 1 / w) ^ 2))
    theorem kernel_sum (A : ℕ) (hA : 3 ≤ A) (M : ℕ) :
        ∑ n ∈ range M, ((n : ℝ) + 2) / ((n : ℝ) + 1) ^ A ≤ 4
    theorem momentM_of_inputs (h : ExpLog ex lg) {κ θ : ℝ} (hθ : 0 < θ) {D : DataB}
        (hin : M3Inputs ex lg κ θ D) : MomentM κ θ D
    lemma class_bound (h : ExpLog ex lg) {Ha Y₁ Y₂ MY MO n σ T L K B G w u : ℝ} (hT : 0 < T)
        (hL : 1 ≤ L) (hK : 0 ≤ K) (hB : 1 ≤ B) (hG : 0 ≤ G) (hσ : 0 < σ) (hw : 1 < w) (hu : 1 < u)
        (hv : 0 < 1 - 1 / u - 1 / w) (hn : n ≤ σ * (T * L)) (hMO : MO ≤ K * (T * L))
        (hMY : MY ≤ K * B * (T * L)) (hY₁ : Y₁ ≤ rp ex lg n (1 - 1 / w) * rp ex lg MY (1 / w))
        (hY₂ : Y₂ ≤ rp ex lg n (1 - 1 / u - 1 / w) * rp ex lg MO (1 / u) * rp ex lg MY (1 / w))
        (hHa : Ha ≤ K * B / L * (G * Y₁ + Y₂)) :
        Ha ≤ K * (K + 1) ^ 2 * T * B ^ 2 * G * rp ex lg σ (1 - 1 / w) +
          K * (K + 1) ^ 2 * T * B ^ 2 * rp ex lg σ (1 - 1 / u - 1 / w)
    theorem step6B_of_classes (hEL : ExpLog ex lg) {Z : Zeros} {β : ℕ → ℝ} {κ θ : ℝ} {Hp : ℝ → ℝ}
        (hθ : 0 ≤ θ) (hden : DensityK Z β ex lg κ) (h6c : Step6ClassesB Z β ex lg κ θ Hp)
        {u w : ℝ} (hu : 1 < u) (hw : 1 < w) (hv : 0 < 1 - 1 / u - 1 / w) (huθ : 2 * u * θ < κ) :
        ∃ (C : ℕ) (K₁ K₂ : ℝ), 0 ≤ K₁ ∧ 0 ≤ K₂ ∧ ∀ c ε : ℝ, 0 < ε → ε < 1 → ∀ η : ℝ, 0 < η →
          ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → ((offSet Z β T).card : ℝ) ≤ 2 * ε * T * lg T →
            ∃ m : ℕ, |Hp T| ≤ (K₁ * Phi ex lg C (2 * θ) κ (1 - 1 / w) c ε m +
              K₂ * Phi ex lg C 0 κ (1 - 1 / u - 1 / w) c ε m) * T + η * T
    theorem PhiB_small (hEL : ExpLog ex lg) (C : ℕ) {κ θ w u c K₁ K₂ : ℝ} (hK₁ : 0 ≤ K₁) (hK₂ : 0 ≤ K₂)
        (hκ : 0 < κ) (hw : 1 < w) (hv : 0 < 1 - 1 / u - 1 / w) (hθw : 2 * θ < κ * (1 - 1 / w))
        (hc : 0 < c) (hcθ : c * (2 * θ) < 1 - 1 / w) :
        ∀ τ : ℝ, 0 < τ → ∃ ε₁ : ℝ, 0 < ε₁ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ ε₁ → ∀ m : ℕ,
          K₁ * Phi ex lg C (2 * θ) κ (1 - 1 / w) c ε m + K₂ * Phi ex lg C 0 κ (1 - 1 / u - 1 / w) c ε m
            ≤ τ
    lemma step7B_meas (h : ExpLog ex lg) {Bad mE T b KB α : ℝ} (hT0 : 0 < T) (hb : 0 < b) (hKB : 0 < KB)
        (hα : 0 < α) (hBad : b * T ≤ Bad) (hHol : Bad ≤ KB * rp ex lg T (1 - α) * rp ex lg mE α) :
        rp ex lg (b / KB) (1 / α) * T ≤ mE
    theorem step7B (hEL : ExpLog ex lg) {I M Hp S Rq Ps Bad mE L h Nc Off : ℝ → ℝ} {Φ : ℝ → ℕ → ℝ}
        {H₀ q s : ℝ} (hH₀ : 0 < H₀) (hq : 2 < q) (hs : 0 < s) (hqs : 2 / q + 1 / s < 1)
        (hL : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 1 ≤ L T)
        (hNc : ∀ T : ℝ, 0 ≤ Nc T)
        (h3 : ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → |I T - (M T + Hp T)| ≤ η * T)
        (hM : ∃ δ₀ S₀ : ℝ, 0 < δ₀ ∧ 0 < S₀ ∧ ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
          S₀ ≤ S T ∧ M T ≤ -3 * δ₀ * (T / 2) * S T * (1 - η))
        (h6 : ∀ ε : ℝ, 0 < ε → ε < 1 → ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
          Off T ≤ 2 * ε * T * L T → ∃ m : ℕ, |Hp T| ≤ Φ ε m * T + η * T)
        (hΦ : ∀ τ : ℝ, 0 < τ → ∃ ε₁ : ℝ, 0 < ε₁ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ ε₁ → ∀ m : ℕ, Φ ε m ≤ τ)
        (h2 : ∀ T : ℝ, Off T ≤ 2 * Nc T)
        (h1 : ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → -Bad T - η * T ≤ I T)
        (hHol : ∃ KG T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → Bad T ≤ KG * rp ex lg (Rq T) (2 / q) *
          rp ex lg (Ps T) (1 / s) * rp ex lg (mE T) (1 - 2 / q - 1 / s))
        (h5 : ∃ K₁ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → Rq T ≤ K₁ * T ∧ Ps T ≤ K₂ * T)
        (hcov : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
          0 ≤ mE T ∧ 0 ≤ h T ∧ h T * L T ≤ H₀ ∧ mE T ≤ h T * Nc T + h T) :
        ∃ ε₀ : ℝ, 0 < ε₀ ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → ε₀ * T * L T < Nc T
    theorem dlam_pos_of_from {Z : Zeros} (hZ : Z.Ordered) {φ₀ : ℝ} (hφ₀ : 0 ≤ φ₀) (hφ₀1 : φ₀ < 1)
        (hU : UFrom Z φ₀) (hs1 : ScaleToInf Z) (hs2 : ScaleHalf Z) (hR : RvMUp Z) {lam : ℝ}
        (hlam : φ₀ < lam) : 0 < Dlam Z lam :=
    noncomputable def rp (ex lg : ℝ → ℝ) (x e : ℝ) : ℝ := if 0 < x then ex (e * lg x) else 0

    -- follow-up (LemmaZ2, LemmaZW, WiredB): Lemma Z from its named inputs, wired into M3 and 6e′ (a)
    def ZPoint1 (ex lg : ℝ → ℝ) {ι : Type*} [DecidableEq ι] (X y CI Kpk : ℝ) (k NT : ℕ) (S : Finset ι)
        (d κ : ι → ℝ) (Rabs : ℝ) : Prop :=
      ∃ (𝒯 : Finset ℕ) (nt : ℕ → ℕ) (P : ℕ → Finset ℕ) (v C F : ℕ → ℕ → ℝ) (E : ℝ), 𝒯.card ≤ NT ∧ E ≤ 1 ∧
        Rabs ≤ ∑ τ ∈ 𝒯, ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) - y) * |F τ p| + E ∧
        ∀ τ ∈ 𝒯, nt τ + 1 ≤ k ∧ (∀ p ∈ P τ, 1 ≤ v τ p ∧ 0 ≤ C τ p) ∧
          ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2)) ≤ Kpk ∧
          ((∀ p ∈ P τ, |F τ p| ≤ CI * ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (nt τ + 1) => S),
              rp ex lg (X / v τ p) (univ.sup' univ_nonempty (fun i => max (d (ρ i) - y) 0)) *
                ∏ i, κ (ρ i)) ∨
            ∀ p ∈ P τ, |F τ p| ≤ CI)
    theorem lemmaZ_point1 (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] {X y CI Kpk : ℝ} {k NT : ℕ}
        {S : Finset ι} {d κ : ι → ℝ} {Rabs : ℝ} (hX : 0 < X) (hy : 0 ≤ y) (hCI : 0 ≤ CI) (hKpk : 0 ≤ Kpk)
        (hk : 1 ≤ k) (hκ : ∀ x, 0 ≤ κ x) (hz : ZPoint1 ex lg X y CI Kpk k NT S d κ Rabs) :
        Rabs ≤ (NT * (Kpk * CI * k) + 1) * ((1 + Wsum ex lg X y S d κ) * (1 + Onsum S κ) ^ (k - 1))
    theorem lemmaZ_point2 (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] {X y CI Kpk : ℝ} {k NT : ℕ}
        {S : Finset ι} {d κ : ι → ℝ} {Rabs : ℝ} (hX : 1 ≤ X) (hy : 0 ≤ y) (hCI : 0 ≤ CI) (hKpk : 0 ≤ Kpk)
        (hk : 1 ≤ k) (hκ : ∀ x, 0 ≤ κ x) (hz : ZPoint2 ex lg X y CI Kpk k NT S d κ Rabs) :
        Rabs ≤ (NT * (Kpk * CI * k) + 1) *
          (rp ex lg X y * ((1 + Wsum ex lg X 0 S d κ) * (1 + Onsum S κ) ^ (k - 1)))
    theorem m3_first_of_lemmaZ (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] {T X q CI Kpk Rq IA B : ℝ}
        {k NT : ℕ} (hX : 0 < X) (hq : 0 ≤ q) (hCI : 0 ≤ CI) (hKpk : 0 ≤ Kpk) (hk : 1 ≤ k)
        (S : ℝ → Finset ι) (d : ι → ℝ) (κ : ℝ → ι → ℝ) (hκ : ∀ t x, 0 ≤ κ t x) (Rabs : ℝ → ℝ)
        (hz : ∀ t, T ≤ t → t ≤ 2 * T → ZPoint1 ex lg X 0 CI Kpk k NT (S t) d (κ t) (Rabs t))
        (hcmp : ∀ c : ℝ, 0 ≤ c → (∀ t, T ≤ t → t ≤ 2 * T → rp ex lg (Rabs t) q ≤
          c * (rp ex lg (1 + Wsum ex lg X 0 (S t) d (κ t)) q *
            rp ex lg ((1 + Onsum (S t) (κ t)) ^ (k - 1)) q)) → Rq ≤ c * IA)
        (hHol : IA ≤ B) : Rq ≤ rp ex lg (NT * (Kpk * CI * k) + 1) q * B
    theorem hairpin_pointwise (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] {X y CI Kpk : ℝ} {k NT : ℕ}
        {S : Finset ι} {d κ : ι → ℝ} {R1 R2 : ℝ} (hX : 1 ≤ X) (hy : 0 ≤ y) (hCI : 0 ≤ CI)
        (hKpk : 0 ≤ Kpk) (hk : 1 ≤ k) (hκ : ∀ x, 0 ≤ κ x) (hR1 : 0 ≤ R1)
        (h1 : ZPoint1 ex lg X y CI Kpk k NT S d κ R1) (h2 : ZPoint2 ex lg X y CI Kpk k NT S d κ R2) :
        R1 * R2 ≤ (NT * (Kpk * CI * k) + 1) ^ 2 * ((rp ex lg X y + Offsum ex lg X S d κ) *
          (1 + Offsum ex lg X S d κ) * ((1 + Onsum S κ) ^ k) ^ 2)
    def M3InputsW (ex lg : ℝ → ℝ) (κ θ : ℝ) (D : DataB) : Prop :=
      ∀ q r w : ℝ, 2 < q → q < κ / θ → 1 < r → 1 < w → q * r * θ < κ * (1 - 1 / w) →
        ∃ (CI Kpk K₀ K₁ T₀ : ℝ) (k NT : ℕ) (X : ℝ → ℝ) (S : ℝ → ℝ → Finset ℕ) (d : ℕ → ℝ)
          (κk : ℝ → ℝ → ℕ → ℝ) (Rabs : ℝ → ℝ → ℝ) (IA IW IO NOn NOff : ℝ → ℝ) (NF : ℝ → ℕ → ℝ),
          0 ≤ CI ∧ 0 ≤ Kpk ∧ 1 ≤ k ∧ 0 ≤ K₀ ∧ 0 ≤ K₁ ∧
          ∀ T : ℝ, T₀ ≤ T → 0 < X T ∧ (∀ t x, 0 ≤ κk T t x) ∧
            (∀ t, T ≤ t → t ≤ 2 * T →
              ZPoint1 ex lg (X T) 0 CI Kpk k NT (S T t) d (κk T t) (Rabs T t)) ∧
            (∀ c : ℝ, 0 ≤ c → (∀ t, T ≤ t → t ≤ 2 * T → rp ex lg (Rabs T t) q ≤
              c * (rp ex lg (1 + Wsum ex lg (X T) 0 (S T t) d (κk T t)) q *
                rp ex lg ((1 + Onsum (S T t) (κk T t)) ^ (k - 1)) q)) → D.Rq q T ≤ c * IA T) ∧
            IA T ≤ rp ex lg (IW T) (1 / r) * rp ex lg (IO T) (1 - 1 / r) ∧
            IO T ≤ K₀ * T ∧
            rp ex lg (IW T) (1 / (q * r)) ≤ rp ex lg T (1 / (q * r)) + NOn T + NOff T ∧
            NOn T ≤ K₀ * rp ex lg T (1 / (q * r)) ∧
            (∃ M : ℕ, NOff T ≤ ∑ a ∈ range M, ex θ ^ (a + 1) * NF T a + rp ex lg T (1 / (q * r))) ∧
            ∀ a : ℕ, NF T a ≤ K₁ * rp ex lg T (1 / (q * r)) * ex (-(κ * a * (1 - 1 / w) / (q * r)))
    theorem m3Inputs_of_wired (h : ExpLog ex lg) {κ θ : ℝ} {D : DataB} (hin : M3InputsW ex lg κ θ D) :
        M3Inputs ex lg κ θ D
    theorem theoremUB_wired {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        {κ : ℝ} (hκ : 0 < κ) (hden : DensityK Z β ex lg κ) (hin : DPCInputsWB Z β ex lg κ) :
        TheoremUB Z κ :=
    theorem D_pos_seven_elevenths_wired {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        (hR : RvMUp Z) (hden : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DensityK Z β ex lg (8 / 7 - δ))
        (hin : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DPCInputsWB Z β ex lg (8 / 7 - δ)) :
        ∀ lam : ℝ, 7 / 11 < lam → 0 < Dlam Z lam :=

## Named hypotheses and sources

DPC = Durkan–Pearce-Crump 2610.28527, Lz = Lamzouri 2611.14486 (TeX line numbers as quoted in PROOF.md, the brief and dpclean's LEAN.md; spot-check in 'TeX check' below). L = log T = lg T, X = T^θ, h = 2πφ/log T = φ/sc T. DataB = the analytic quantities of one Step-0′ construction (φ, θ, P, f, η fixed before T): I, M, Hp, S, Bad, mE as in dpclean's Data; Rq q T = ∫_T^{2T}|R_T|^q, Ps s T = ∫_T^{2T}|P(V_φ)|^s. DataB.base = dpclean's Data (R4 := Rq 4, P4 := Ps 4, unused) so that Step3, MainTerm, Step1Pos, Covering are reused verbatim.

| hypothesis | Lean statement (for T ≥ T₀; o(T) = ∀η > 0 ∃T₀ … ≤ ηT) | PROOF.md step | source |
|---|---|---|---|
| Zeros.Ordered, Reflect, ExpLog, hsc, RvMUp | dpclean's (Riemann/Deep/DPCUncond/LEAN.md) | Theorem U′ (multiplicity), Step 2, D(λ) | Lz l.176–177; functional equation; Mathlib exp/log facts; Riemann–von Mangoldt |
| Step3 D.base | \|I − (M + 𝓗)\| ≤ o(T) | Step 3 (unchanged) | DPC l.335–388, 698 |
| MainTerm D.base | S₀ ≤ Σ\|c_n\|², Re M_T ≤ −3δ₀(T/2)Σ\|c_n\|²(1 − η) | Steps 0, 3, 7 (unchanged) | DPC Prop 1 l.262–388; Lemma 6 l.686–702; construction l.743–771 |
| Step1Pos D.base | −Bad − o(T) ≤ I | Step 1; Step 7′ | DPC l.760–765, 774 |
| Covering D.base | mE ≤ (φ/sc T)·𝒩(T) + φ/sc T | Step 7 | DPC l.774, 782 |
| DensityK κ | #{T ≤ γ ≤ 2T, (β − ½) log T ≥ a} ≤ C_d T log T e^{−κa}, all a ∈ ℕ | Theorem U′'s hypothesis; M2–M4, 6e′ | Conrey κ = 8/7 − δ (Baluyot's thesis, Riemann/Deep/Baluyot); Lz l.221–230 |
| Step6ClassesB | ∀u, w > 1 with 1 − 1/u − 1/w > 0, 2uθ < κ ∃C, K and class quantities: (6a) \|𝓗_T\| ≤ Σ_{a≤M(T)} 𝓗_{T,a} + o(T); (a)–(b) 𝓗_{T,a} ≤ K(a+2)^C L⁻¹(e^{2θ(a+1)}Y₁ + Y₂); (d) Y₁ ≤ \|S_a\|^{1−1/w}MY^{1/w}, Y₂ ≤ \|S_a\|^{1/v}MO^{1/u}MY^{1/w}; (c) MY ≤ K(a+2)^C TL, MO ≤ KTL | 6a; 6e′ (a)–(d) | Lemma Z (PROOF.md), Hölder in y and over S_a; M4 (Lz DoubleSumZeros l.251–283, Fujii l.242–249); 6b–6d with all zeros |
| HolderB | ∀q > 2, s > 1, 2/q + 1/s < 1: Bad ≤ K_G (Rq)^{2/q}(Ps)^{1/s} mE^{1−2/q−1/s} | Step 7′ | Hölder (q/2, s, 1/α), G_T ≤ K_G (as DPC l.776–780) |
| MomentM | ∀q ∈ (2, κ/θ): Rq q T ≤ K₁T | M3 | Lemma Z + M1–M3 |
| M3Inputs (fine; replaces MomentM in InputsFineB) | ∀q ∈ (2, κ/θ), r, w > 1, qrθ < κ(1 − 1/w), p = qr: Rq ≤ K_Z(IW)^{1/r}(IO)^{1−1/r}; IO ≤ K₀T; IW^{1/p} ≤ T^{1/p} + NOn + NOff; NOn ≤ K₀T^{1/p}; NOff ≤ Σ_{a<M} e^{θ(a+1)}NF_a + T^{1/p}; NF_a ≤ K₁T^{1/p}e^{−κa(1−1/w)/p} | M1–M3 | Lemma Z at y = 0 + Hölder; Fujii (M1, Step 5b, Lz l.242–249); Minkowski (W ≤ On + Off); Off ≤ Σ_a e^{θ(a+1)}F_a + O(T^{−5}); M2 + Minkowski over m + kernel sum |
| Step5cB | ∀s > 1: Ps s T ≤ K₂T | 5c | M. Riesz DPC l.641; Fujii Lz l.242–249; density Lz l.221–230 |
| DPCInputsB, DPCInputsFineB | ∀(φ, θ), 0 < φ < 1, 0 < θ < min{κ/2, 1}, 0 < φθ < 1 < φ(1+θ) < 2: ∃ DataB with InputsB (InputsFineB) | Step 0′ | DPC l.743–772 (sign condition l.753, order l.772) |
| LemmaZ.term_bound hF, hP; lemmaZ_sum hR, hB | \|F(𝔞)\| ≤ C_I Σ_{ρ⃗∈S^m}(X/𝔞)^{max_i c(ρ_i)}∏κ(ρ_i); Σ_𝔞 C(𝔞)𝔞^{−1/2} ≤ K_pk; \|R_T(t − iy)\| ≤ Σ_τ B_τ + E | Lemma Z (i)–(iv) | explicit formula (IK §5.5; Rodgers Thm 2.1 as cited DPC l.643); Möbius inversion over set partitions; binomial expansion; Fourier inversion; Fubini; integration by parts; Mertens |
| LemmaM.M2_bound hI, hμ, hB, hN | Hölder; meas ≤ \|B_a\|/L; \|B_a\| ≤ CTLe^{−κa}; ∫n_m^{pw} ≤ CT | M2 | Hölder; density bins; Fujii |
| LemmaM4.M4_count hN, hB, hS; M4_moment hI, hS | covering + Hölder; \|B_a\| ≤ CTLe^{−κa}; Σ_{ρ'}Ñ^w ≤ C(\|m\|+2)^w TL; Hölder; Σ_ρ n_m(x_ρ)^{pw} ≤ C'TL | M4 | Lz DoubleSumZeros; 6c–6d's Ñ moments |

## Build

Each module separately, one at a time, lean directly (no lake process), oleans to .scratch/dpcblean/olean (build output, deleted at the end):

    rm .lake && bash .scratch/memcap.sh 290 3600 lake exe cache get        # 14:59–15:01 UTC: peak 243 MB, rc = 0
    cp Riemann/Deep/DPCUncond/Lean/audit/leanpath.txt Riemann/Deep/DPCBeyond/Lean/audit/leanpath.txt     # lake's LEAN_PATH
    # per module m (dir d): pplean's FinSum Model Liminf, dpclean's Params Gaps Tail Phi Step7 TheoremU DLam Hairpin Fine, then
    # ParamsB Rpow Step7B HairpinB SumB DLamB TheoremUB FinalB LemmaZ LemmaM LemmaM4 LemmaM3 FineB LemmaZ2 LemmaZW CheckB WiredB (d = Riemann/Deep/DPCBeyond/Lean):
    LEAN_PATH=.scratch/dpcblean/olean:$(cat Riemann/Deep/DPCBeyond/Lean/audit/leanpath.txt) bash .scratch/memcap.sh 290 1800 \
      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 -o .scratch/dpcblean/olean/$d/$m.olean $d/$m.lean

Scripts and logs: Riemann/Deep/DPCBeyond/Lean/audit/verify.sh (clean rebuild of all 29 modules), mk.sh (one module), verify.log, verify1.log (first full run), peaks.log (every compile of the session), logs/log_<module>.txt, review_prompt.txt, review_out.md, statements.txt. rg over Lean/*.lean: no sorry, admit, axiom, native_decide, decide, skipKernelTC. CheckB.lean: 72 #print axioms lines (WiredB.lean: 4 more), all [propext, Classical.choice, Quot.sound]:

    'DPCBLean.theoremUB' depends on axioms: [propext, Classical.choice, Quot.sound]
    (same line for step0B_nonempty, step0B_signs, step8B_empty, step0B_iff, exists_q_s, exists_u_w, order_of_choices, conrey_threshold,
     conrey_sharp, rp_add, rp_mul, rp_rp, rp_le_rp, rp_lt_rp, rp_le_rp_exp, rp_le_max, step7B_meas, step7B, class_count, class_bound,
     step6B_of_classes, PhiB_small, dpos_of_from, nToInf_of_from, dlam_pos_of_from, cover_of, moments_of, split_c, D_pos_of_kappa,
     theoremUB_conrey, D_pos_seven_elevenths, ex_sum, sum_mul_le, prod_ex_le, prod_ex_le_sup, rp_sup_le_sum, tuple_sum, tuple_bound,
     packet_case1, packet_case2, packet_nofree, pos_part_shift, pos_part_split, term_bound, mono_bound, lemmaZ_sum, M2_bound, geom_le,
     M3_series, M3_exists_w, M4_count, M4_moment, inv_sq_sum, kernel_sum, off_bound, momentM_of, exists_r_w, momentM_of_inputs,
     inputsB_of_fine, theoremUB_fine, D_pos_seven_elevenths_fine, W_split, tuple_case2, term_bound2, term_nofree1,
     term_nofree2, final_shape, lemmaZ_point1, lemmaZ_point2, m3_first_of_lemmaZ, hairpin_pointwise; WiredB: m3Inputs_of_wired,
     inputsFineB_of_wired, theoremUB_wired, D_pos_seven_elevenths_wired)

## Proof route

- (a) ParamsB: interval and sign arithmetic via 1/φ − 1 < θ ⟺ 1 < (1 + θ)φ; q = (2 + κ/θ)/2, 1/s = (1 − 2/q)/2; u = (1 + κ/(2θ))/2, w = u/(u − 1) + 1 (then 1/w < 1 − 1/u < 1 − 2θ/κ); conrey_threshold: δ = min{1, (22φ/7 − 2)/(2φ)}.
- (b) LemmaZ: ex(Σx_i) = ∏ex(x_i) (Finset induction) and Σd_iu_i ≤ DΣu_i ≤ D·lg(X/𝔞), ex monotone; the max is attained (Finset.exists_mem_eq_sup'); tuple_sum: a(ρ_j)∏_ik(ρ_i) = ∏_i[i = j ? ak : k](ρ_i) (Finset.mul_prod_erase), Finset.sum_prod_piFinset, then Σ_j over Fin m; term_bound: (X/𝔞)^{c} ≤ X^{c}, 𝔞^{−1/2−y} ≤ 𝔞^{−1/2}, tuple_bound, packet sum; lemmaZ_sum: m W On^{m−1} ≤ k(1+W)(1+On)^{k−1} per term, #𝒯 terms.
- (c) LemmaM, LemmaM4: rp monotone in the base, (xy)^e = x^e y^e, x^{e}x^{f} = x^{e+f}, x^e ≤ max{1, x} for e ∈ [0, 1]; geometric sums by (Σρ^a)(1 − ρ) = 1 − ρ^M; Σ_{n≤k}1/(n+1)² ≤ 2 − 1/(k+1) by induction (1/(x+1)² ≤ 1/x − 1/(x+1)); LemmaM3: ‖Off‖_p ≤ (K₁e^θ/(1−ρ) + 1)T^{1/p} (off_bound), ‖1+W‖_p ≤ K₃T^{1/p}, IW = (IW^{1/p})^p ≤ K₃^pT, then ∫|R_T|^q ≤ K_Z(K₃^pT)^{1/r}((K₀+1)T)^{1−1/r} = KT; exists_r_w: r = (1 + κ/(qθ))/2 and M3_exists_w at p = qr.
- (d) HairpinB, SumB: class_count as dpclean's 6a counts with slope κ; class_bound: (σTL)^{e₁}(K′BTL)^{1/w} etc. with e₁ + 1/w = 1 and 1/v + 1/u + 1/w = 1, (K′B)^{1/w} ≤ K′B; step6B_of_classes: dpclean's hairpin_sum with (θ, E) = (2θ, 1 − 1/w) and (0, 1/v); PhiB_small: dpclean's Phi_small twice (2θ < κ(1 − 1/w), c·2θ < 1 − 1/w; 0 < κ/v).
- (e) Step7B: δ₁ = δ₀S₀/4 from dpclean's step7_bad; Bad ≤ (|K_G|+1)((|K₁|+1)T)^{2/q}((|K₂|+1)T)^{1/s}mE^α = K_B T^{1−α}mE^α; step7B_meas: T^{1−α} = T·T^{−α}, (c₂T)^α = (δ₁/K_B)T^α ≤ mE^α, strict monotonicity; ε₀ = min{ε₁, ½, c₂/(2H₀)} and c₂TL ≤ H₀(𝒩 + 1) ≤ (c₂/2)TL + H₀ is false for T > 2H₀/c₂ + 1.
- (f) theoremUB: θ, q, s, u, w from (a); step6B_of_classes, PhiB_small with c = (1 − 1/w)/(4θ); dpclean's card_off_le (Step 2), covering (cover_of, copied from dpclean's theoremU), step7B with L = lg, h = φ/sc T, H₀ = 2πφ; DLamB: dpclean's transfer with φ = (φ₀ + min{λ, 1})/2; FinalB: δ from conrey_threshold.

## Fidelity (each Lean object against PROOF.md and the TeX)

| Lean | PROOF.md / source | match |
|---|---|---|
| UFrom, TheoremUB | Theorem U′ | exact: (γ_{n+1} − γ_n)·sc T < φ ⟺ γ_{n+1} − γ_n < 2πφ/log T (dpclean's gapSet); c·T·sc T = (c/2π)T log T; φ ∈ (max{2/(2+κ), ½}, 1) |
| theoremUB, theoremUB_fine | Theorem U′ from Steps 0′–7′ | conditional on DensityK and DPCInputsB (DPCInputsFineB) |
| theoremUB_conrey, D_pos_seven_elevenths(_fine) | "With κ = 8/7 − δ: every φ > 7/11, so D(λ) > 0 for every λ > 7/11" | exact; δ ∈ (0, 8/7) picked per φ or λ (conrey_threshold); inputs at every κ = 8/7 − δ (Conrey's slope holds for each fixed δ > 0; no δ = 0 endpoint is assumed) |
| D_pos_of_kappa | discovery's RH-relevance line: λ₀(κ) = max{2/(2+κ), ½} | 0 < D(λ) for every λ > max{2/(2+κ), ½} |
| step0B_nonempty, step0B_signs, step8B_empty, step0B_iff | Step 0′ ("What changes": θ ∈ (1/φ − 1, min{κ/2, 1}), nonempty iff φ > max{2/(2+κ), ½}; 0 < φθ < 1 < φ(1+θ) < 2, θ < 1); DPC l.753 | exact; emptiness stated for φ > 0 |
| exists_q_s, exists_u_w, order_of_choices | Order of choices: q ∈ (2, κ/θ), s with 2/q + 1/s < 1; u ∈ (1, κ/(2θ)), w > u/(u−1) with κ(1 − 1/w) > 2θ, 1/v = 1 − 1/u − 1/w > 0 | exact; κ(1 − 1/w) > 2θ is proved from w > u/(u − 1) and u < κ/(2θ) |
| conrey_threshold, conrey_sharp | referee C′: "κ = 8/7 − δ gives exactly φ > 7/11" | exact: 2/(2 + 8/7) = 7/11; no φ ≤ 7/11 is above any threshold with δ ∈ [0, 8/7) |
| rp, Rpow lemmas | x^e (Real.rpow) | rp x e = ex(e·lg x) for x > 0, 0 for x ≤ 0: Real.rpow on [0, ∞) for e ≠ 0 (instance not compiled, cap) |
| step7B_meas, step7B | Step 7′ | exact: δ₁T ≤ K T^{2/q+1/s} meas(E_T)^α ⇒ meas(E_T) ≥ c₂T, c₂ = (δ₁/K)^{1/α} independent of ε; the rest of Step 7 (dpclean's step7_bad, covering, contradiction) unchanged |
| HolderB, Step5cB | Step 7′ Hölder; 5c (‖P(V)‖_s ≪ T^{1/s}) | exact (G_T ≤ K_G; K_G may depend on q, s: weaker) |
| MomentM | M3: ∫_T^{2T}\|R_T\|^q ≪ T for q < κ/θ | restricted to q > 2, the range used (reviewer note 6) |
| M3Inputs, off_bound, momentM_of, exists_r_w, momentM_of_inputs | M3: "‖Off‖_p ≪ T^{1/p} Σ_a e^{θ(a+1) − κa(1−1/w)/p} ≪ T^{1/p} whenever pθ < κ (choose w large), and by Lemma Z and Hölder ∫\|R_T\|^q ≪ T" | M3's Hölder split (p = qr, r′ = r/(r − 1)) and Minkowski made explicit; IW = ∫(1+W)^{qr}, IO = ∫(1+On)^{q(k−1)r′}, NOn, NOff, NF_a = p-norms of On, Off, F_a |
| M2_bound | M2: (n^{(a)}_m)^p ≤ 1_{≥1}n_m^p, Hölder, (\|B_a\|/L)^{1−1/w}(CT)^{1/w} | exact: ∫(n^{(a)}_m)^p ≤ CTe^{−κa(1−1/w)} (p-th power moment; its p-th root and Minkowski over m give M3Inputs' NF bound) |
| geom_le, M3_series | M3's bin series | every partial sum ≤ e^θ/(1 − e^{θ−κ(1−1/w)/p}) (convergence of a positive series) |
| M3_exists_w | "choose w large" (M3); "(1 − 1/w)² → κ" (M4, B′) | explicit w = 4κ/(κ − pθ) for 0 ≤ pθ < κ |
| M4_count, M4_moment | M4 (repaired radius (\|m\|+1)/L, referee B′) | exact: ≪ (\|m\|+2)TLe^{−κa(1−1/w)}, then ≪ (\|m\|+2)TLe^{−κa(1−1/w)²} ((\|m\|+2)^{1−1/w} ≤ \|m\|+2) |
| key_ineq, inv_sq_sum, kernel_sum | M4 / B′: Σ_m (1+\|m\|)^{−A}(\|m\|+2) < ∞ for A large | for natural A ≥ 3, over n = \|m\| ∈ ℕ: ≤ 4 (the ℤ-sum ≤ 8) |
| DensityK | Theorem U′'s N(σ, 3T) ≪ T^{1−κ(σ−½)} log T at σ = ½ + a/L, γ ∈ [T, 2T] | implied (N(σ, 2T) ≤ N(σ, 3T)); M2/M4's bins over [T/4, 3T] enter via their own named \|B_a\| bounds |
| Step6ClassesB, class_count | 6a, 6e′ (a)–(d): σ_a = min(2ε + o(1), Ce^{−κa}) | the class bounds hold for the actual \|S_a\| = classSet card; σ_a = min(3ε, C_d e^{−κa}) formed from Step 2's #Z_off ≤ 2εTL (dpclean's card_off_le) and DensityK |
| class_bound | 6e′ (d): 𝓗_{T,a} ≪ T(a+2)^C[e^{2θ(a+1)}σ_a^{1−1/w} + σ_a^{1/v}] | exact exponent arithmetic (constant K(K+1)², (a+2)^{2C}) |
| step6B_of_classes, PhiB_small | 6e′ (e): Σ_a𝓗_{T,a} ≪ Φ(ε)T, Φ(ε) → 0 (A = c log(1/ε), 0 < c < (1−1/w)/(2θ); σ^{1/v} series) | Φ = K₁·Phi(2θ, κ, 1 − 1/w) + K₂·Phi(0, κ, 1/v) with dpclean's Phi (split at A, tail Σ_{a>A}); c = (1 − 1/w)/(4θ); both proved → 0 uniformly in the truncation |
| sum_mul_le, prod_ex_le, prod_ex_le_sup | Lemma Z (iii): \|e^{Σz_iu_i}\| ≤ e^{max(Re z_i)⁺Σu_i} ≤ (X/𝔞)^{max(Re z_i)⁺} | exact on the support (u_i ≥ 0, Σu_i ≤ log(X/𝔞)); ex(Σ) = ∏ex proved |
| rp_sup_le_sum, tuple_sum, tuple_bound | (iv): X^{max c_i} ≤ Σ_i X^{c_i}; Σ_{ρ_1..ρ_m} X^{max} ∏κ_A ≤ m·W·On^{m−1} | exact over every finite set of zeros (full sums = suprema, terms ≥ 0) |
| packet_case1, packet_case2, packet_nofree, pos_part_shift | (iv): 𝔞^{−1/2−y}(X/𝔞)^{…} ≤ 𝔞^{−1/2}X^{…}; 𝔞^{−1/2+y}(X/𝔞)^{y+max d⁺} ≤ X^{y+max d⁺}𝔞^{−1/2}; 𝔞 < X ⇒ 𝔞^{−1/2+y} ≤ X^y𝔞^{−1/2}; (d+y)⁺ ≤ y + d⁺ | exact |
| pos_part_split | 6e′ (a): W_{−y} ≤ On + X^{−y}Off | exact per zero (summation linear) |
| term_bound, mono_bound, lemmaZ_sum | Lemma Z (iii)–(iv) for each partition/binomial term; \|R_T\| ≤ C(1 + W_{−y})(1 + On)^{k−1} | case 1 (R_T(t − iy)); (i)–(ii) and (iii)'s IBP bound named; case 2 (R*) has the same shape with packet_case2, pos_part_shift (W_− ≤ X^yW) |
| Wsum, Onsum, Offsum, W_split | Lemma Z: W_{−y} = Σ_ρX^{(d_ρ−y)⁺}κ_A, On, W = W_{−0}; 6e′ (a): W_{−y} ≤ On + X^{−y}Off | exact over every finite set of zeros |
| tuple_case2, term_bound2, term_nofree1, term_nofree2 | Lemma Z (iv), case 2 (𝔞^{−1/2+y}(X/𝔞)^{y+max d⁺} ≤ X^{y+max d⁺}𝔞^{−1/2}) and terms with no free variable (𝔞 < X ⇒ 𝔞^{−1/2+y} ≤ X^y𝔞^{−1/2}) | exact; the case-2 free exponent y + max d⁺ is (iii) with D = y + max d⁺ ≥ max(d_i + y)⁺ (on a nonempty support X/𝔞 > 1; follow-up re-review note 2) |
| ZPoint1, ZPoint2, final_shape, lemmaZ_point1, lemmaZ_point2 | Lemma Z (i)–(iv), both displays | named inputs (i)–(iii) bundled per point; conclusions exactly |R_T(t − iy)| ≤ C(1 + W_{−y})(1 + On)^{k−1}, |R*(t − iy)| ≤ CX^y(1 + W)(1 + On)^{k−1}, C = NT K_pk C_I k + 1 |
| m3_first_of_lemmaZ, M3InputsW, m3Inputs_of_wired, inputsFineB_of_wired, theoremUB_wired, D_pos_seven_elevenths_wired | M3: "by Lemma Z and Hölder ∫|R_T|^q ≪ T" | Lemma Z at y = 0 + hcmp (monotone, homogeneous integral on the two integrable functions) + hHol (Hölder r, r′); other M3 inputs unchanged |
| hairpin_pointwise | 6e′ (a): |R_T R*|(x − iy) ≤ C[X^y + Off](1 + Off)(1 + On)^{2k} | exact ((1 + On)^{2k} written ((1 + On)^k)², C = K²) |
| dpos_of_from, nToInf_of_from, dlam_pos_of_from | Theorem U′ ⇒ D(λ) (dpclean's dpos_of_U, Repair 7) | dpclean's proof with 2/3 replaced by φ₀ ∈ [0, 1) |

## TeX check (16:03 UTC; e-prints curled to .scratch/dpcblean, deleted after; keyword digest by texcheck.py, no bulk printing)

Lz SmallGapsPaper.tex (575 lines): l.177 ρ_n = β_n + iγ_n; l.223–224 Eq:Conrey N(σ, T) ≪_ε T^{1−(8/7−ε)(σ−1/2)} log T (one log T, uniform in σ ≥ 1/2: DensityK with κ = 8/7 − ε); l.240–246 Lem:Fujii; l.251–257 Lem:DoubleSumZeros. DPC Gaps_for_Dirichlet.tex v2 (1,093 lines): l.262 prop:finite-moments (Prop 1); l.337–341 the shift and its O(e^{−cL²}) sides (Step3); l.619 0 < θ < 1, f vanishing for x ≥ θ; l.641 ∫|F_T|^p ≪ T (5c); l.643 the explicit-formula choice of F̂ (Lemma Z (i)); l.686 lem:proportion-rough-limit (Lemma 6, MainTerm); l.744 "1/φ − 1 < θ < 1"; l.753 "0 < φθ < 1 < φ(1+θ) < 2" (DPCInputsB); l.760–763 the polynomial signs (Step1Pos); l.774 E_T; l.776–779 δT ≤ ∫_{E_T}G_T|R_T|²|P(V_φ)| with meas(E_T)^{1/4} (Hölder (2, 4, 4), replaced by HolderB). All cited lines match their use.

## Review

Reviewer (verdict-only, statement fidelity; PROOF.md, dpclean's reviewed statements, pplean's zero model, the abstraction notes and all new Lean statements pasted inline; prompt and answer in Riemann/Deep/DPCBeyond/Lean/audit/review_prompt.txt, review_out.md): **FAITHFUL**, nine [minor] notes, no [blocking]: (1) conrey_sharp is the parameter obstruction, not necessity of the method's threshold; (2)–(3) Lemma Z's mechanism, packet exponents and summation correct, the explicit formula, IBP and packet summability remain named inputs; (4) M2_bound bounds a p-th power moment (the 1/p enters in M3), M4 correctly loses twice (κ(1 − 1/w)²), kernel_sum over |m| costs a constant; (5) M3Inputs are genuine intermediate hypotheses, not the final moment estimate; (6) MomentM records only 2 < q < κ/θ, the range Theorem U′ uses; (7) 6e′ bookkeeping faithful (factor 2θ, constants before ε, uniform in the truncation); (8) Step 7′ exponents correct, c₂ independent of ε, missing sign assumptions only weaken hypotheses; (9) final theorems' quantifiers correct, δ > 0 with slack, the inputs are non-vacuous; a conditional formalization.

Follow-up re-reviews (verdict-only, new statements only; Riemann/Deep/DPCBeyond/Lean/audit/review2_*, review3_*): first ISSUES, one [blocking] (m3_first_of_lemmaZ's integral functional was assumed monotone on all functions, false for Lean's integral of non-integrable functions), fixed by stating monotonicity, homogeneity (hcmp) and Hölder (hHol) only for the two integrable functions used; then FAITHFUL with four [minor] notes: hcmp and hHol are standard for those functions, IW and IO need local boundedness of W and On (uniform convergence or tail domination of the zero sums), ZPoint1/2 need the truncation/tail realization for the actual polynomial, T₀ > 0.

## Not done

- The named hypotheses are not proved in Lean: Theorem U's DPC construction, Prop 1 and Lemma 6, the slit identity (H), Step 1's positivity, Hölder and the covering bound (measure theory), M. Riesz, Fujii, Lamzouri's DoubleSumZeros, Conrey's density theorem, Lemma Z (i)–(iii) (explicit formula, Möbius inversion over set partitions, binomial expansion, Fourier inversion, Fubini, integration by parts), Mertens, Minkowski, M1, 6b–6d with all zeros, Riemann–von Mangoldt. They need Mathlib's complex analysis, measure theory and Real.log, above the 290 MB cap.
- Lemma Z (follow-up, 16:29–17:05 UTC): both cases are now assembled (LemmaZ2: term_bound2, term_nofree1/2, lemmaZ_point1/2) and Lemma Z is a theorem from its named inputs (i)–(iii) bundled per point in ZPoint1/ZPoint2 (the decomposition into ≤ NT partition/binomial terms, packet sums, free-sum bounds, remainder E ≤ 1). Wired: (1) M3's first inequality — M3InputsW (WiredB) replaces M3Inputs' first clause by ZPoint1 at y = 0 at every t ∈ [T, 2T] plus two named integral facts stated only for the two integrable functions used (hcmp: ∫_T^{2T} monotone and homogeneous on |R_T|^q and (1+W)^q((1+On)^{k−1})^q; hHol: Hölder with r, r′), and m3Inputs_of_wired proves M3Inputs (theoremUB_wired, D_pos_seven_elevenths_wired); (2) 6e′ (a) — hairpin_pointwise proves |R_T R*|(x − iy) ≤ K²(X^y + Off)(1 + Off)(1 + On)^{2k} from ZPoint1, ZPoint2 and W_split. NOT wired: Step6ClassesB's clause (a)–(b) 𝓗_{T,a} ≤ K(a+2)^C L⁻¹(e^{2θ(a+1)}Y₁ + Y₂) is still a named hypothesis; turning it into a theorem from hairpin_pointwise needs inputs not yet named: (6a′) 𝓗_{T,a} ≤ K₀Σ_{ρ∈S_a}J_ρ with J_ρ = ∫_0^{d_ρ}|R_T R*|(x_ρ − iy)Z(x_ρ − iy)dy; (y-int) J_ρ ≤ A∫_0^{d_ρ}X^yZ dy + B∫_0^{d_ρ}Z dy whenever the integrand is ≤ AX^yZ + BZ on [0, d_ρ] (A, B ≥ 0, integrable functions); (b) ∫_0^{d_ρ}X^yZ dy ≤ C(a+1)e^{θ(a+1)}L⁻¹Z̄_ρ and ∫_0^{d_ρ}Z dy ≤ Y_aZ̄_ρ (Hölder in y, ∫_0^{d}X^{ry}dy ≤ X^{rd}/(r log X)); ZPoint1, ZPoint2 at every (x_ρ, y), y ∈ [0, d_ρ]; Y₁, Y₂ defined as Σ_{S_a}(1+On)^{2k}Z̄_ρ, Σ_{S_a}Off²(1+On)^{2k}Z̄_ρ (then Step6ClassesB's Hölder clause (d) refers to these sums). The remaining algebra is elementary ((1 + Off)(C e^{θ(a+1)} + Off) ≤ (2C+2)(e^{2θ(a+1)} + Off²)); it did not fit the 45-minute follow-up. ZPoint1/2 work over a finite set of zeros; applying them to the actual polynomial needs the truncation/tail realization (re-review notes).
- ExpLog's instance (Real.exp, Real.log, tp = 2π) is not compiled (cap), as in dpclean.
- .lake: \`lake exe cache get\` (14:59–15:01 UTC, memcap peak 243 MB, rc = 0) left .lake as a real directory (7.4 GB: Mathlib and dependencies with their cached oleans); kept on purpose (user, follow-up: other Lean agents share the real .lake; do not restore the placeholder). No lake build was run; this run's oleans went to .scratch/dpcblean/olean (build output, deleted at the end).
- Memory cap: during development two compiles of an earlier single-file HairpinB exceeded the cap (time -l footprint 305 and 306 MiB at 15:22 and 15:24 UTC; the 3-s watchdog sampled 289 and 276 MB and did not kill them); the file was split (HairpinB, SumB) and every nlinarith in large contexts removed. One LemmaM compile reached 292 MiB (above the 290 watchdog line, below the 300 cap) before it was split (LemmaM, LemmaM4). All final compiles are below 290 MiB (table above).

## Progress
- 14:59 read brief, PROOF.md, discovery, dpclean LEAN.md and Lean files; started lake exe cache get under memcap 290/3600 (.lake placeholder removed).
- Plan: namespace DPCBLean, files ParamsB, RpowB, LemmaZ, LemmaM, HairpinB, Step7B, TheoremUB, DLamB, CheckB; reuse DPCLean (Tail, Phi, Hairpin, Step7, Gaps, TheoremU, DLam) by import.
- 15:25 compiled (lean -j1, memcap 290/1800, /usr/bin/time -l): deps (pplean FinSum Model Liminf; dpclean Params Gaps Tail Phi Step7 TheoremU DLam Hairpin Fine) 249–279 MiB; ParamsB 273, Rpow 253, Step7B 283, HairpinB 282, SumB 271 MiB.
- 15:22–15:24 CAP BREACH: HairpinB (one file, nlinarith in big contexts) reached footprint 305 then 306 MiB (watchdog sampled 289 and 276 MB, no kill); fixed by removing nlinarith, linarith only, and splitting into HairpinB + SumB (282, 271 MiB).
- 15:28–15:49 compiled TheoremUB 263, DLamB 272, FinalB 253, LemmaZ 272, LemmaM 263 (292 before split), LemmaM4 277, LemmaM3 274, FineB 251, CheckB 251 MiB; CheckB: 62 theorems, all [propext, Classical.choice, Quot.sound].
- 15:50 started clean verification run (Riemann/Deep/DPCBeyond/Lean/audit/verify.sh).
- 15:50–16:02 first clean verification run (verify1.log): all 26 compiles rc = 0, footprints 249–286 MiB (Step7B 286).
- 15:53–15:56 reviewer (statement fidelity, all material inline): FAITHFUL, nine [minor] notes.
- 16:02 Step7B: every linarith in the long proofs restricted (linarith only) for memory margin; 16:03 TeX spot-check of the cited lines (Lz, DPC v2); 16:04 second clean verification run (verify.log) for the final record.
- 16:04–16:17 second clean verification run (verify.log): all 26 compiles rc = 0, footprints 250–282 MiB (Step7B 277), 62 #print axioms lines standard; LEAN.md assembled; oleans, e-prints deleted; .lake left as a real directory (shared with other Lean agents).
- Follow-up 16:29–17:05 UTC: LemmaZ2 (case 2, ZPoint1/2, lemmaZ_point1/2, W_split) 276 MiB; LemmaZW (m3_first_of_lemmaZ, hairpin_pointwise) 274 MiB; WiredB (M3InputsW, m3Inputs_of_wired, theoremUB_wired, D_pos_seven_elevenths_wired). Re-review 1 (new statements, Riemann/Deep/DPCBeyond/Lean/audit/review2_*): ISSUES, one [blocking]: the first version of m3_first_of_lemmaZ asked for an integral functional monotone on all functions (false for Lean's integral of non-integrable functions); fixed by stating monotonicity/homogeneity and Hölder only for the two integrable functions used (hcmp, hHol); re-review 2 (review3_*): FAITHFUL, four [minor] notes (finite moments IW, IO need local boundedness of W, On; truncation/tail realization of ZPoint; T₀ > 0). Third clean verification run (verify.log).

## FORMALIZATIONS.md entry

- Theorem U′ of our discovery [small-gap-seven-elevenths](discoveries/small-gap-seven-elevenths.md) (Durkan–Pearce-Crump's small-gap method past length T^{1/2} without RH) in Lean 4 + Mathlib: DPCBLean.theoremUB(_fine) (𝒩(T) ≥ c(φ)T log T for φ ∈ (max{2/(2+κ), 1/2}, 1)) and DPCBLean.D_pos_seven_elevenths(_fine) (0 < D(λ) for every λ > 7/11, Lamzouri's liminf D, Conrey's κ = 8/7 − δ) from named inputs (Theorem U's DPC construction, Prop 1, Lemma 6, slit identity, Step 1, covering; the density slope; 6e′'s class inequalities; Step 7′'s Hölder; Lemma M3 or its inputs (Lemma Z + Hölder, Fujii, Minkowski, depth bins, M2); 5c; Riemann–von Mangoldt; exp/log facts), with Steps 0′/8′ and the order of choices, Lemma Z's mechanism and summation and Lemma Z itself in both cases from its named inputs (i)–(iii) (wired into M3's first inequality and 6e′ (a)), Lemma M's bookkeeping (M2, M3, M4, kernel sum, M3's assembly), 6e′'s bookkeeping and Φ(ε) → 0, Step 7′ and Theorem U′ ⇒ D(λ) > 0 proved; sorry-free, standard axioms, 2445 lines, every final compile under the 290 MB watchdog (≤ 282 MiB; two development compiles reached 305–306 MiB and were fixed by splitting). Reviewer: FAITHFUL (follow-up re-review: see Review). [LEAN.md](Riemann/Deep/DPCBeyond/LEAN.md)
