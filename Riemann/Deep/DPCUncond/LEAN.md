FORMALIZED: Theorem U (DPCLean.theoremU_fine; for every φ ∈ (2/3, 1) there is c > 0 with 𝒩(T) ≥ c·T·log T/(2π) for all large T) and D(λ) > 0 for every λ > 2/3 with Lamzouri's liminf D (DPCLean.D_pos_two_thirds_fine), from named analytic inputs: DPC's construction with Prop 1 and Lemma 6 (MainTerm), the slit identity (Step3), Step 1's positivity (Step1Pos), Hölder (Holder7), Montgomery–Vaughan (MV5a), M. Riesz + Fujii + density (Step5c), the covering bound (Covering), the Jutila/Conrey density theorem (Density), the 6a–6e class bounds (Step6Classes), Riemann–von Mangoldt's upper bound (RvMUp), pplean's zero model (Zeros.Ordered), the reflection ρ ↦ 1 − ρ̄ (Reflect) and four Mathlib facts about exp/log (ExpLog). Proved, not assumed: Step 2 (#off-line or multiple zeros in [T, 2T] ≤ 2𝒩(T), list and zero-model forms), Steps 0 and 8 (θ ∈ (1/φ − 1, 1/2) exists iff φ > 2/3, and gives 0 < φθ < 1 < φ(1 + θ) < 2), 6e's Hölder bookkeeping (E = (1 − 1/w)(1/2 − 1/(3r)), growth e^{θ(a+1)}) and explicit δ, w, r with θ < (1 − δ)E, 6a's class counts, 6f's class sum and Φ(ε) → 0 (summability proved), 5a's use of θ ≤ 1/2, Step 7's contradiction (𝒩(T) > ε₀TL), and Theorem U ⇒ D(λ) > 0. Sorry-free; axioms propext, Classical.choice, Quot.sound only; no native_decide, no decide, no skipKernelTC.

# DPCUncond (Theorem U) in Lean 4 + Mathlib (dpclean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile). Riemann/Deep/DPCUncond/Lean/, 10 files, 1,459 lines, namespace DPCLean; not imported by Riemann.lean. Imports pplean's zero model and liminf D unchanged (Riemann/Deep/ExplorePointProc/Lean/Model.lean, Liminf.lean, FinSum.lean: PPLean.Zeros, Zeros.Ordered, smallGaps, DPos, Dlam, dpos_iff, ScaleToInf). No Mathlib zeta, Fourier, LSeries, complex analysis, measure theory, Real.log, Real.rpow or π (290 MB cap): integrals are abstract real functions of T, exp/log are abstract (ExpLog), every analytic input is a named hypothesis.

| file | content | lines | footprint (time -l) | memcap reading |
|---|---|---|---|---|
| Gaps | (a) shared, zeroGaps, card_shared_le; Reflect, Exc, excSet, offSet, gapSet (𝒩), card_exc_le, card_off_le, card_gapSet_le | 151 | 254 MiB | 245 MB |
| Params | (b) step0_nonempty, step0_signs, step8_empty, step0_iff, step8_threshold, step8_above; (c) Eexp, exp_mean, holder_bookkeeping, exists_params | 145 | 266 MiB | 248 MB |
| Tail | ExpLog (+ zero, nat_mul, lg_mul, lg_one, lg_inv, le_lg, lt_one), pow_mul_le_one, polygeom, geom_Ioc, geom_tail_le, (d) tail_small | 195 | 279 MiB | 268 MB |
| Phi | (d) Aeps, Phi, first_small, Phi_small | 115 | 265 MiB | 251 MB |
| Hairpin | lg_ex, hairpin_sum (6f's class sum) | 144 | 271 MiB | 0 (between samples) |
| Step7 | (e) step7_bad, step7_holder, step7 | 161 | 272 MiB | 241 MB |
| TheoremU | TheoremU, Data, Step3, MainTerm, Step6, Step1Pos, Holder7, Step5, Covering, Inputs, DPCInputs, theoremU | 166 | 265 MiB | 244 MB |
| DLam | (f) ScaleHalf, RvMUp, scaleToInf_of, scaleHalf_of, gapSet_half_sub, dpos_of_U, nToInf_of_U, dlam_pos_of_U, D_pos_two_thirds | 156 | 279 MiB | 257 MB |
| Fine | deepSet, classSet, Density, Step6Classes, step6_of_classes, MV5a, Step5c, step5_of, InputsFine, DPCInputsFine, inputs_of_fine, theoremU_fine, D_pos_two_thirds_fine, theoremU_seq | 194 | 270 MiB | 208 MB |
| Check | #print axioms (28 theorems) | 32 | 248 MiB | 204 MB |

## Main statements (verbatim; proofs omitted)

    -- Theorem U and D(λ) (Gaps, TheoremU, DLam, Fine); pplean: Z.γ n (n ≥ 1, with multiplicity), Z.N T, Z.sc T = log T/(2π)
    noncomputable def gapSet (Z : Zeros) (φ T : ℝ) : Finset ℕ :=        -- 𝒩(T)
      (Icc 1 (Z.N (2 * T))).filter
        (fun n => T ≤ Z.γ n ∧ Z.γ (n + 1) ≤ 2 * T ∧ (Z.γ (n + 1) - Z.γ n) * Z.sc T < φ)
    def TheoremU (Z : Zeros) : Prop :=
      ∀ φ : ℝ, 2 / 3 < φ → φ < 1 → ∃ c : ℝ, 0 < c ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
        c * T * Z.sc T ≤ ((gapSet Z φ T).card : ℝ)
    theorem theoremU_fine {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β) (hEL : ExpLog ex lg)
        {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T) (hden : Density Z β ex lg)
        (hin : DPCInputsFine Z β ex lg) : TheoremU Z
    theorem D_pos_two_thirds_fine {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
        (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
        (hR : RvMUp Z) (hden : Density Z β ex lg) (hin : DPCInputsFine Z β ex lg) :
        ∀ lam : ℝ, 2 / 3 < lam → 0 < Dlam Z lam                         -- Dlam = Lamzouri's liminf D (pplean)
    theorem theoremU (hZ : Z.Ordered) (hβ : Reflect Z β) (hEL : ExpLog ex lg) (htp : 0 < tp)
        (hsc : ∀ T, 0 < T → Z.sc T * tp = lg T) (hin : DPCInputs Z β ex lg) : TheoremU Z   -- coarse inputs
    theorem D_pos_two_thirds … (hR : RvMUp Z) (hin : DPCInputs Z β ex lg) : ∀ lam, 2 / 3 < lam → 0 < Dlam Z lam
    theorem dpos_of_U (hZ : Z.Ordered) (hU : TheoremU Z) (hs1 : ScaleToInf Z) (hs2 : ScaleHalf Z) (hR : RvMUp Z)
        {lam : ℝ} (hlam : 2 / 3 < lam) : DPos Z lam
    theorem dlam_pos_of_U … : 0 < Dlam Z lam
    theorem theoremU_seq (hU : TheoremU Z) (t : ℕ → ℝ) (ht : ∀ M, ∃ k₀, ∀ k, k₀ ≤ k → M ≤ t k) (h1 : 2 / 3 < φ)
        (h2 : φ < 1) : ∃ c, 0 < c ∧ ∃ k₀, ∀ k, k₀ ≤ k → c * t k * Z.sc (t k) ≤ ((gapSet Z φ (t k)).card : ℝ)

    -- (a) Step 2 (Gaps)
    noncomputable def shared (x : ℕ → ℝ) (n : ℕ) : Finset ℕ := (range n).filter (fun i => ∃ j ∈ range n, j ≠ i ∧ x j = x i)
    noncomputable def zeroGaps (x : ℕ → ℝ) (n : ℕ) : Finset ℕ := (range n).filter (fun i => i + 1 < n ∧ x (i + 1) = x i)
    theorem card_shared_le (x : ℕ → ℝ) (n : ℕ) (hx : ∀ i, i + 1 < n → x i ≤ x (i + 1)) :
        (shared x n).card ≤ 2 * (zeroGaps x n).card
    def Reflect (Z : Zeros) (β : ℕ → ℝ) : Prop :=
      ∀ n : ℕ, 1 ≤ n → β n ≠ 1 / 2 → ∃ m : ℕ, 1 ≤ m ∧ Z.γ m = Z.γ n ∧ β m = 1 - β n
    def Exc (Z : Zeros) (β : ℕ → ℝ) (n : ℕ) : Prop := β n ≠ 1 / 2 ∨ ∃ m : ℕ, 1 ≤ m ∧ m ≠ n ∧ Z.γ m = Z.γ n ∧ β m = β n
    excSet Z β T = {n ∈ Icc 1 N(2T) : T ≤ γ n ≤ 2T, Exc n};  offSet Z β T = {n ∈ Icc 1 N(2T) : T ≤ γ n ≤ 2T, 1/2 < β n}
    theorem card_exc_le (hZ : Z.Ordered) (hβ : Reflect Z β) (hφ : 0 < φ) (T : ℝ) :
        (excSet Z β T).card ≤ 2 * (gapSet Z φ T).card
    theorem card_off_le … : (offSet Z β T).card ≤ 2 * (gapSet Z φ T).card

    -- (b) Steps 0 and 8 (Params)
    theorem step0_nonempty (h1 : 2 / 3 < φ) : ∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < 1 / 2
    theorem step0_signs (h1 : 2 / 3 < φ) (h2 : φ < 1) (h3 : 1 / φ - 1 < θ) (h4 : θ < 1 / 2) :
        0 < φ * θ ∧ φ * θ < 1 ∧ 1 < φ * (1 + θ) ∧ φ * (1 + θ) < 2
    theorem step8_empty (h0 : 0 < φ) (h : φ ≤ 2 / 3) : ¬ ∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < 1 / 2
    theorem step0_iff (h0 : 0 < φ) : (∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < 1 / 2) ↔ 2 / 3 < φ
    theorem step8_threshold (hθ : 0 < 1 + θ) : 1 < φ * (1 + θ) ↔ 1 / (1 + θ) < φ
    theorem step8_above (h0 : 0 ≤ θ) (h1 : θ < 1 / 2) : 2 / 3 < 1 / (1 + θ)

    -- (c) Step 6e (Params)
    noncomputable def Eexp (w r : ℝ) : ℝ := (1 - 1 / w) * (1 / 2 - 1 / (3 * r))
    theorem exp_mean (w r : ℝ) :
        (1 - 1 / w) / 4 + (1 - 1 / w) / 4 + (1 - 1 / w) * (1 - 1 / r) = Eexp w r + Eexp w r + Eexp w r
    theorem holder_bookkeeping {s : ℝ → ℝ} (hs : ∀ x y, s (x + y) = s x * s y) (hs0 : ∀ x, 0 ≤ s x)
        {w r K T B g A₁ A₂ A₃ H : ℝ} (hK : 0 ≤ K) (hT : 0 ≤ T) (hB : 0 ≤ B) (hg : 0 ≤ g) (hA₂ : 0 ≤ A₂) (hA₃ : 0 ≤ A₃)
        (h1 : A₁ ≤ K * T * B * s ((1 - 1 / w) / 4)) (h2 : A₂ ≤ K * T * B * g ^ 3 * s ((1 - 1 / w) / 4))
        (h3 : A₃ ≤ K * T * B * s ((1 - 1 / w) * (1 - 1 / r))) (hH : H ^ 3 ≤ A₁ * A₂ * A₃) :
        H ≤ K * T * B * g * s (Eexp w r)
    theorem exists_params (hθ : θ < 1 / 2) : ∃ δ w r : ℝ, 0 < δ ∧ δ < 1 ∧ 1 < w ∧ 1 < r ∧ θ < (1 - δ) * Eexp w r
        -- δ = x, w = 1/x, r = 1/(3x), x = (1/2 − max θ 0)/4

    -- (d) Step 6f (Tail, Phi, Hairpin)
    structure ExpLog (ex lg : ℝ → ℝ) : Prop where
      add : ∀ x y : ℝ, ex (x + y) = ex x * ex y
      pos : ∀ x : ℝ, 0 < ex x
      mono : StrictMono ex
      ex_lg : ∀ x : ℝ, 0 < x → ex (lg x) = x
    noncomputable def Aeps (lg : ℝ → ℝ) (c ε : ℝ) : ℕ := ⌊c * lg (1 / ε)⌋₊
    noncomputable def Phi (ex lg : ℝ → ℝ) (C : ℕ) (θ κ E c ε : ℝ) (M : ℕ) : ℝ :=
      ((Aeps lg c ε : ℝ) + 2) ^ C * ex θ ^ (Aeps lg c ε + 1) * ex (lg (3 * ε) * E)
        + ∑ a ∈ Ioc (Aeps lg c ε) M, ((a : ℝ) + 2) ^ C * ex θ * ex (θ - κ * E) ^ a
    theorem tail_small (C : ℕ) {q ρ : ℝ} (hq : 0 ≤ q) (hρ0 : 0 ≤ ρ) (hρ1 : ρ < 1) :
        ∀ τ : ℝ, 0 < τ → ∃ N : ℕ, ∀ A M : ℕ, N ≤ A → ∑ a ∈ Ioc A M, ((a : ℝ) + 2) ^ C * q * ρ ^ a ≤ τ
    theorem Phi_small (h : ExpLog ex lg) (C : ℕ) {θ κ E c : ℝ} (hE : 0 ≤ E) (hκ : θ < κ * E) (hc : 0 < c)
        (hcθ : c * θ < E) : ∀ τ : ℝ, 0 < τ → ∃ ε₁ : ℝ, 0 < ε₁ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ ε₁ → ∀ M : ℕ,
          Phi ex lg C θ κ E c ε M ≤ τ
    theorem hairpin_sum (h : ExpLog ex lg) {C : ℕ} {θ κ E c ε K T Cd : ℝ} (M : ℕ) (H s : ℕ → ℝ)
        (hθ : 0 ≤ θ) (hE : 0 ≤ E) (hK : 0 ≤ K) (hT : 0 ≤ T) (hCd : 0 < Cd)
        (hs0 : ∀ a, 0 < s a) (hsε : ∀ a, s a ≤ 3 * ε) (hsd : ∀ a : ℕ, s a ≤ Cd * ex (-(κ * a)))
        (hH : ∀ a : ℕ, H a ≤ K * T * ((a : ℝ) + 2) ^ C * ex θ ^ (a + 1) * ex (lg (s a) * E)) :
        ∑ a ∈ range (M + 1), H a ≤ K * max 1 (ex (lg Cd * E)) * Phi ex lg (C + 1) θ κ E c ε M * T

    -- (e) Step 7 (Step7): all objects real functions of T; h3 … hcov in PROOF.md's order
    theorem step7 {I M Hp S R4 P4 Bad mE L h Nc Off : ℝ → ℝ} {Φ : ℝ → ℕ → ℝ} {H₀ : ℝ}
        (hH₀ : 0 < H₀) (hL : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 1 ≤ L T) (hNc : ∀ T : ℝ, 0 ≤ Nc T)
        (h3 : ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → |I T - (M T + Hp T)| ≤ η * T)
        (hM : ∃ δ₀ S₀ : ℝ, 0 < δ₀ ∧ 0 < S₀ ∧ ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
          S₀ ≤ S T ∧ M T ≤ -3 * δ₀ * (T / 2) * S T * (1 - η))
        (h6 : ∀ ε : ℝ, 0 < ε → ε < 1 → ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T →
          Off T ≤ 2 * ε * T * L T → ∃ m : ℕ, |Hp T| ≤ Φ ε m * T + η * T)
        (hΦ : ∀ τ : ℝ, 0 < τ → ∃ ε₁ : ℝ, 0 < ε₁ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ ε₁ → ∀ m : ℕ, Φ ε m ≤ τ)
        (h2 : ∀ T : ℝ, Off T ≤ 2 * Nc T)
        (h1 : ∀ η : ℝ, 0 < η → ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → -Bad T - η * T ≤ I T)
        (hHol : ∃ KG T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ Bad T ∧ Bad T ^ 4 ≤ KG * R4 T ^ 2 * P4 T * mE T)
        (h5 : ∃ K₁ K₂ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ R4 T ∧ R4 T ≤ K₁ * T ∧ 0 ≤ P4 T ∧ P4 T ≤ K₂ * T)
        (hcov : ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → 0 ≤ mE T ∧ 0 ≤ h T ∧ h T * L T ≤ H₀ ∧ mE T ≤ h T * Nc T + h T) :
        ∃ ε₀ : ℝ, 0 < ε₀ ∧ ∃ T₀ : ℝ, ∀ T : ℝ, T₀ ≤ T → ε₀ * T * L T < Nc T

    -- Steps 5, 6 from finer inputs (Fine)
    theorem step6_of_classes (hEL : ExpLog ex lg) (hθ : 0 ≤ θ) (hden : Density Z β ex lg)
        (h6c : Step6Classes Z β ex lg θ D) : Step6 Z β ex lg θ D
    theorem step5_of (hEL : ExpLog ex lg) (hθ : θ ≤ 1 / 2) (h5a : MV5a ex lg θ D) (h5c : Step5c D) : Step5 D
    theorem inputs_of_fine (hEL : ExpLog ex lg) (hden : Density Z β ex lg) (hin : DPCInputsFine Z β ex lg) :
        DPCInputs Z β ex lg

## Named hypotheses and sources

DPC = Durkan–Pearce-Crump 2610.28527v2, Gaps_for_Dirichlet.tex (curl of https://arxiv.org/e-print/2610.28527v2, 1,093 lines; lines 686–702, 743–756, 766–787 read for this file; download deleted). Lz = Lamzouri 2611.14486, SmallGapsPaper.tex (line numbers as quoted in PROOF.md and pplean). L = log T = lg T, h = 2πφ/log T = φ/sc T. Data = the analytic quantities of one Step-0 construction (φ, θ, q, b₀, P, f, η fixed in DPC's order, l.772): I T = Re ∫_T^{2T} G_T|R_T|²P(V_φ), M T = Re M_T, Hp T = Re 𝓗_T, S T = Σ|c_n|², R4 T = ∫_T^{2T}|R_T|⁴, P4 T = ∫_T^{2T}|P(V_φ)|⁴, Bad T = ∫_{E_T} G_T|R_T|²|P(V_φ)|, mE T = meas E_T.

| hypothesis | Lean statement (all "for T ≥ T₀"; o(T) = ∀ η > 0 ∃ T₀ … ≤ η T) | PROOF.md step | source |
|---|---|---|---|
| Zeros.Ordered (pplean) | γ_n nondecreasing, γ_1 > 0, γ_k ≤ T ⟺ k ≤ N(T) | Theorem U (ordinates with multiplicity) | Lz l.176–177 |
| Reflect | β n ≠ 1/2 ⇒ ∃ m ≥ 1, γ m = γ n, β m = 1 − β n | Theorem U ("an off-line pair gives two equal ordinates"), Step 2 | functional equation; Lz l.177 |
| ExpLog | ex(x + y) = ex x · ex y, ex > 0, ex strictly increasing, ex(lg x) = x for x > 0 | 6f (e^θ, (3ε)^E, log(1/ε)) | Mathlib Real.exp_add, Real.exp_pos, Real.exp_strictMono, Real.exp_log |
| hsc | sc T · tp = lg T for T > 0 (tp = 2π) | Theorem U (2πφ/log T) | pplean's sc T = log T/(2π) |
| RvMUp | N(T) ≤ K T sc T | Theorem U ⇒ D(λ) ("≫ N(T)") | Riemann–von Mangoldt, Lz l.177 |
| Step3 | \|I T − (M T + Hp T)\| ≤ o(T) | Step 3, identity (H) | DPC l.335–388 (Prop 1 after the shift), l.698 (outside [T, 2T]: O(e^{−cL²})) |
| MainTerm | ∃ δ₀ S₀ > 0: S₀ ≤ S T and M T ≤ −3δ₀(T/2) S T (1 − η) for every η > 0, T large | Steps 0, 3, 7 | DPC Prop 1 l.262–388; Lemma 6 l.686–696 (limiting quotient; l.695 denominator/T → positive constant); l.700 (∫G_T\|R_T\|² ~ (T/2)Σ\|c_n\|²); construction l.743–771 |
| Step1Pos | −Bad T − o(T) ≤ I T | Step 1 (Im V counts all zeros), Step 7 | DPC l.774 (outside E_T, Im V ∈ {0, 1}, Re P(V) ≥ 0 by the two-line signs l.760–765) |
| Holder7 | Bad⁴ ≤ KG · R4² · P4 · mE | Step 7 | DPC l.776–780 (Hölder 2, 4, 4) |
| Step5 (coarse) | R4 ≤ K₁ T, P4 ≤ K₂ T, both ≥ 0 | Step 5 (constants independent of ε) | 5a, 5c |
| MV5a (fine) | R4 ≤ (2T + K e^{2θ log T}) B | 5a | Montgomery–Vaughan mean value theorem; Σ\|b_n\|² ≪_η 1 (Mertens) |
| Step5c (fine) | P4 ≤ K₂ T | 5b–5c | M. Riesz DPC l.641; Fujii Lz l.242–249; density Lz l.221–230 |
| Covering | mE ≤ (φ/sc T)·𝒩(T) + φ/sc T | Step 7 ("each consecutive pair with gap < h covers centres of measure ≤ h") | DPC l.774, 782 |
| Step6 (coarse) | ∀ δ ∈ (0,1), w, r > 1 with θ < (1 − δ)E ∃ C, K ≥ 0 ∀ c ∈ (0, E/θ) ∀ ε ∈ (0,1): #offSet ≤ 2εTL ⇒ ∃ m, \|Hp T\| ≤ K Φ(ε) T + o(T) | Step 6 (6a–6f), Repair 6 | Jutila's κ = 1 − δ, Lz l.223–225 |
| Density (fine) | ∀ δ ∈ (0,1) ∃ C_d > 0: #{n ∈ offSet : (β n − 1/2) log T ≥ a} ≤ C_d T log T e^{−(1−δ)a} for all a ∈ ℕ | 6a, 6d, 6f | Jutila 1983 / Conrey (Baluyot), uniform for σ ≥ 1/2, Lz l.221–230 |
| Step6Classes (fine) | ∀ w, r > 1 ∃ C, K, class totals Ha, A₁, A₂, A₃ and a truncation M(T): \|Hp\| ≤ Σ_{a ≤ M(T)} Ha + o(T); Ha³ ≤ A₁A₂A₃; \|S_a\| ≤ sTL ⇒ A₁ ≤ KT(a+2)^C s^{(1−1/w)/4}, A₂ ≤ KT(a+2)^C e^{3θ(a+1)} s^{(1−1/w)/4}, A₃ ≤ KT(a+2)^C s^{(1−1/w)(1−1/r)} | 6a (decomposition; outside [T, 2T] o(1), Step 3), 6b–6e | Poisson representations, sparse kernels (6c), Fujii Lz l.242–249, Lz Lemma DoubleSumZeros l.251–283 |
| DPCInputs / DPCInputsFine | for every (φ, θ) with 0 < φ < 1, 0 < θ < 1/2, 0 < φθ < 1 < φ(1 + θ) < 2: ∃ Data satisfying the inputs above | Step 0 | DPC l.743–772 (sign condition l.753, order of choices l.772) |

## Build

Each module separately, one at a time, lean directly (no lake process), oleans to .scratch/dpclean/olean (build output, deleted at the end):

    cp Riemann/Deep/ExplorePointProc/Lean/audit/leanpath.txt .scratch/dpclean/leanpath.txt        # lake's LEAN_PATH
    # per module m (dir d): pplean's FinSum Model Liminf (d = Riemann/Deep/ExplorePointProc/Lean), then
    # Params Gaps Tail Phi Step7 TheoremU DLam Hairpin Fine Check (d = Riemann/Deep/DPCUncond/Lean):
    LEAN_PATH=.scratch/dpclean/olean:$(cat .scratch/dpclean/leanpath.txt) bash .scratch/memcap.sh 290 1800 \
      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 -o .scratch/dpclean/olean/$d/$m.olean $d/$m.lean

Final run (Lean/audit/verify.sh, 14:01–14:06 UTC): all 13 rc = 0, 0 errors, 0 warnings; peak memory footprint (/usr/bin/time -l) 248–279 MiB (table above; pplean's FinSum 255, Model 264, Liminf 251 MiB), 9–27 s each; no watchdog kill in the whole run. During development the first single-proof Step7 reached 286 MiB and was split (step7_bad, step7_holder: 272 MiB). By lake instead: lake build Riemann.Deep.DPCUncond.Lean.Check (not run: lake's own process adds to the cap). rg over the 10 files: no sorry, admit, axiom, native_decide, decide, skipKernelTC. Check.lean:

    'DPCLean.card_shared_le' depends on axioms: [propext, Classical.choice, Quot.sound]
    (same line for card_exc_le, card_off_le, step0_nonempty, step0_signs, step8_empty, step0_iff, step8_threshold,
     step8_above, exp_mean, holder_bookkeeping, exists_params, polygeom, tail_small, first_small, Phi_small, step7,
     theoremU, dpos_of_U, dlam_pos_of_U, D_pos_two_thirds, hairpin_sum, step6_of_classes, step5_of, inputs_of_fine,
     theoremU_fine, D_pos_two_thirds_fine, theoremU_seq: 28 theorems, all [propext, Classical.choice, Quot.sound])

## Proof route

- (a) card_shared_le, card_exc_le: an entry sharing its value with an entry to its right has a zero gap at i, one to its left a zero gap at i − 1; so the set injects into G ∪ (G + 1) (no run decomposition). In the zero model the partner is the reflected zero (Reflect) or the repeated one; the zero gap lies inside [T, 2T], so it counts in 𝒩(T) (0 · sc T < φ). Result ≤ 2𝒩(T), without PROOF's + O(L).
- (b) Interval and sign arithmetic (div_lt_iff₀, nlinarith).
- (c) holder_bookkeeping: A₁A₂A₃ ≤ (KTB)³ g³ s(e₁)s(e₁)s(e₃) = (KTBg s(E))³ by s(x + y) = s x s y and 2e₁ + e₃ = 3E (ring); cube monotone. exists_params: (1 − x)²(1/2 − x) ≥ 1/2 − 2x > θ.
- (d) polygeom: s = 1 − (1 − t)/(C + 2) has t ≤ s^{C+2} (Bernoulli), s^a(1 + au) ≤ 1 (Bernoulli), a + 2 ≤ K₀(1 + au); so (a + 2)^C t^a ≤ K₀^C s^a; the tail is a geometric tail. first_small: A ≤ c log(1/ε) gives q^{A+1}(3ε)^E ≤ e^θ 3^E (e^{θ − E/c})^A with e^{θ − E/c} < 1 (cθ < E). Phi_small: ε₁ = 1/e^{N/c} makes A(ε) ≥ N. hairpin_sum: split the a-sum at A(ε); a ≤ A uses s_a ≤ 3ε (A + 1 terms, absorbed by the exponent C + 1), a > A uses s_a ≤ C_d e^{−κa}, i.e. s_a^E ≤ C_d^E (e^{−κE})^a.
- (e) step7: δ₁ = δ₀S₀/4; Φ(ε₀) ≤ δ₀S₀/2 and the o(T)'s give I ≤ −(3/8)δ₀S₀T, Bad ≥ δ₁T (step7_bad); Hölder and Step 5 give δ₁⁴TL ≤ Kc H₀(𝒩 + 1) (step7_holder); with 𝒩 ≤ ε₀TL and ε₀ ≤ δ₁⁴/(2KcH₀) this forces TL ≤ 2KcH₀/δ₁⁴, false for large T.
- theoremU: θ from step0_nonempty, signs from step0_signs, (δ, w, r) from exists_params, c = E/(2θ), Φ(ε) := K·Phi, then step7 with L = lg, h = φ/sc T (h·L = 2πφ), Off ≤ 2𝒩 (card_off_le). Fine: step6_of_classes takes s_a = min(3ε, C_d e^{−κa}) (Step 2 and Density bound |S_a|), holder_bookkeeping per class, hairpin_sum; step5_of: e^{2θ log T} ≤ T for θ ≤ 1/2.
- (f) dpos_of_U: φ = (2/3 + min(λ, 1))/2; ScaleHalf (from ExpLog: sc T = sc(T/2) + log 2/(2π)) and ScaleToInf give φ sc T ≤ λ sc(T/2), so gapSet φ (T/2) ⊆ smallGaps λ T; 𝒩_φ(T/2) ≥ c(T/2)sc(T/2) ≥ (c/4)T sc T ≥ (c/4(|K| + 1))N(T) (RvMUp). N(T) → ∞ from Theorem U itself (nToInf_of_U); pplean's dpos_iff gives 0 < Dlam.

## Fidelity (each Lean object against PROOF.md and the TeX)

| Lean | PROOF.md / source | match |
|---|---|---|
| Zeros, Ordered, smallGaps, DPos, Dlam, dpos_iff (pplean, imported) | Lz l.176–183, D(λ) = liminf N(T)⁻¹#{0 < γ_n ≤ T : γ_{n+1} − γ_n ≤ 2πλ/log T} | exact (pplean LEAN.md, reviewer FAITHFUL) |
| gapSet | Theorem U's 𝒩(T) | exact: (γ_{n+1} − γ_n)·sc T < φ ⟺ γ_{n+1} − γ_n < 2πφ/log T for sc T = log T/2π > 0; γ_n ≤ γ_{n+1} from Ordered; the ambient Icc 1 N(2T) removes nothing |
| TheoremU | Theorem U | exact; c·T·sc T = (c/2π)T log T, the 2π absorbed into c |
| D_pos_two_thirds(_fine) | "Hence D(λ) > 0 for every λ > 2/3"; Repair 7 | exact, with Lamzouri's liminf D; transfer on [T/2, T] with 2/3 < φ < min(λ, 1) as in Repair 7 |
| theoremU_seq | brief: T → ∞ by ℕ-indexed sequences | T is real throughout (Theorem U and the liminf D are over real T); the conclusion along any T_k → ∞ is a corollary |
| card_shared_le | Step 2, "runs of length r ≥ 2 contain r − 1 zero gaps and r ≤ 2(r − 1)" | exact (sorted list, entries sharing a value ≤ 2·#zero gaps) |
| Reflect, Exc, card_exc_le, card_off_le | Step 2: #{off-line or multiple, γ ∈ [T, 2T]} ≤ 2εTL + O(L) under 𝒩(T) ≤ εTL | sharper: ≤ 2𝒩(T), no O(L); Reflect asks only for a partner index (implied by the multiset invariance) |
| step0_nonempty, step0_signs, step0_iff, step8_empty | Step 0 ("θ ∈ (1/φ − 1, 1/2), nonempty iff φ > 2/3; 0 < φθ < 1 < φ(1 + θ) < 2"), DPC l.744, 753 | exact; emptiness needs φ > 0 (for φ ≤ 0 the interval is nonempty, irrelevant) |
| step8_threshold, step8_above | Step 8 first bullet | the sign condition φ(1 + θ) > 1 ⟺ φ > 1/(1 + θ) and 1/(1 + θ) > 2/3 for θ < 1/2; the operator identity ⟨1, (I + aS)⁻¹1⟩ = tan(aθ)/a is not formalized |
| Eexp, exp_mean, holder_bookkeeping | 6e: A₁ ≪ Tσ^{(1−1/w)/4}, A₂ ≪ Te^{3θ(a+1)}σ^{(1−1/w)/4}, A₃ ≪ Tσ^{(1−1/w)(1−1/r)}, 𝓗_{T,a} ≤ (A₁A₂A₃)^{1/3} ≪ T(a+2)^C e^{θ(a+1)}σ^E | exact exponent arithmetic; σ^e abstract (any s with s(x + y) = s x s y ≥ 0; ex(lg σ · e) in Fine); one constant K and B = (a+2)^C for all three (PROOF's (a+2), (a+2), (a+2)^C are ≤ (a+2)^C); 𝓗³ ≤ A₁A₂A₃ replaces the cube root |
| exists_params | 6f: "fix w, r so that κE > θ, possible since θ < 1/2 and E ↑ 1/2", κ = 1 − δ | explicit δ, w, r for every θ < 1/2 |
| ExpLog | Real.exp, Real.log | the four Mathlib facts; the instance (Real.exp, Real.log, tp = 2π) is not compiled (Mathlib.Analysis.SpecialFunctions.Log.Basic alone is 326 MB, OddRobin/LEAN.md) |
| Aeps, Phi | Repair 6: A = c log(1/ε), 0 < c < E/θ; 6f: Φ(ε) = (A+2)^C e^{θ(A+1)}(3ε)^E + Σ_{a>A}(a+2)^C e^{θ(a+1)−κEa} | exact terms: q = e^θ, ρ = e^{θ−κE}, e^{θ(a+1)−κEa} = qρ^a, (3ε)^E = e^{log(3ε)·E} (Real.rpow_def_of_pos); A = ⌊c log(1/ε)⌋ (an integer split point); the a-sum truncated at M (the hairpin sum at height T is finite, depths < 1/2) |
| tail_small, first_small, Phi_small | 6f: first term ≪ (log(1/ε) + 2)^C e^θ ε^{E−cθ} → 0, tail ≪ (A+2)^C e^{−(κE−θ)A} → 0 | proved, uniformly in M (equivalent to the series statement); the summability of (a+2)^C ρ^a is proved (Bernoulli), not assumed |
| hairpin_sum | 6f first display, Σ_a 𝓗_{T,a} ≪ T[…] =: Φ(ε)T | proved from per-class bounds and σ_a ≤ min(3ε, C_d e^{−κa}) (6a's σ_a = min(2ε + o(1), Ce^{−κa})) |
| step7 | Step 7 (DPC l.766–787 with Steps 1–6) | exact logic; constants of Step 5 and Hölder are independent of ε by construction (the objects do not depend on ε); ε₀ = min(ε₁, 1/2, δ₁⁴/(2KcH₀)) plays PROOF's min(ε₁, c₂/(4πφ)) |
| Step3 | Step 3: ∫ G_T\|R_T\|²P(V) = M_T + 𝓗_T + o(T) | real parts, implied |
| MainTerm | Steps 0, 3, 7; DPC Lemma 6 l.686–702 | implied: lim_T Re M_T/∫G_T\|R_T\|² ≤ −3δ₀ and ∫G_T\|R_T\|² ~ (T/2)Σ\|c_n\|² give Re M_T ≤ −3δ₀(T/2)Σ\|c_n\|²(1 + o(1)); S₀ > 0 is l.695 (the denominator/T → positive constant) |
| Step6 | Step 6, Repair 6 | the premise is Step 2's output (#Z_off ∩ [T, 2T] ≤ 2εTL), the only use of the contradiction hypothesis in 6a–6f; quantifiers in PROOF's order (κ via δ; w, r) → C, K → c → ε → T |
| Density, Step6Classes, step6_of_classes | 6a–6f | Density: N(σ, 2T) ≪ T^{1−κ(σ−1/2)} log T at σ = 1/2 + a/log T; class bounds stated for every upper bound s of \|S_a\|/(TL), as the 6c–6e proofs use them; 6a's counts, 6e and 6f proved |
| Step1Pos | Step 1, DPC l.774 | implied (DPC's inequality is exact; here − o(T)) |
| Holder7 | DPC l.776–780 | exact, fourth powers instead of fourth roots, KG ≥ sup G_T⁴ |
| Step5, MV5a, Step5c, step5_of | 5a, 5c | MV5a implied by 5a (∫_T^{2T} ≤ ∫_{T/2}^{5T/2}); step5_of proves 5a's ∫\|R_T\|⁴ ≪ T from θ ≤ 1/2 |
| Covering | Step 7, DPC l.774, 782 | sharper than PROOF's 𝒩(T) ≥ c₂T/h − O(L): pairs not counted by 𝒩(T) cover only [T, T + h/2) ∪ (2T − h/2, 2T] |
| DPCInputs, DPCInputsFine | Step 0 | DPC's construction for every admissible (φ, θ); φ, θ, P, f, η fixed before T (l.772) |
| RvMUp, ScaleHalf, ScaleToInf | Riemann–von Mangoldt; log T − log(T/2) = log 2; log T → ∞ | RvMUp named (weaker than RvM: any K); ScaleHalf, ScaleToInf proved from ExpLog and hsc |

Reviewer (verdict-only, statement fidelity; all Lean statements, pplean's zero model, PROOF.md Theorem U, Notation, Steps 0–3, 5–8, Inputs, Repairs and DPC TeX l.686–702, 743–756, 766–787 pasted inline; Lean/audit/review_out.md): FAITHFUL, nine [minor] notes, no [blocking]: a conditional formalization (DPCInputs, DPCInputsFine remain assumptions); Step3, MainTerm faithful (S₀ from the positive denominator limit; the 1 − η form handles the negative main term); Step6, Step6Classes match the class estimates and the order of choices; Step1Pos, Holder7, Step5, MV5a, Step5c faithful (restricting MV to [T, 2T] valid); Covering's additive h valid; Reflect faithful; Density, RvMUp faithful consequences; ExpLog, hsc satisfied by Real.exp, Real.log, tp = 2π but do not themselves fix that normalization; TheoremU, gapSet faithful, 2π absorbed; D_pos_two_thirds(_fine) conclude positivity of the liminf itself.

## Not done

- The named hypotheses are not proved in Lean: DPC's construction, Prop 1 and Lemma 6, the slit identity (H), Step 1's positivity, Hölder and the covering bound (measure theory), Montgomery–Vaughan, M. Riesz, Fujii, the density theorem, the 6a–6e class bounds, Riemann–von Mangoldt. They are analytic and need Mathlib's complex analysis, measure theory and Real.log, which exceed the 290 MB cap.
- ExpLog's instance with Mathlib's Real.exp, Real.log (and tp = 2 * Real.pi) is not written as a compiled file (cap); it is four lemma applications.
- Step 8's operator identity (tan(aθ)/a, the exactness of 1/(1 + θ)) and Step 0's polynomial algebra (q, b₀, the two-line signs of P) are not formalized; they sit inside MainTerm and Step1Pos.
- .lake: left as it is (real directory, another agent's Mathlib cache); no lake command was run and find shows no file under .lake newer than this run's brief. My oleans went to .scratch/dpclean/olean (deleted at the end).

## FORMALIZATIONS.md entry

- Theorem U of our discovery [small-gap-two-thirds](discoveries/small-gap-two-thirds.md) (Durkan–Pearce-Crump's small-gap method without RH) in Lean 4 + Mathlib: DPCLean.theoremU_fine (𝒩(T) ≥ c(φ)T log T for φ ∈ (2/3, 1)) and DPCLean.D_pos_two_thirds_fine (0 < D(λ) for every λ > 2/3, Lamzouri's liminf D) from named inputs (DPC's construction, Prop 1 and Lemma 6; the slit identity; Step 1's positivity; Hölder; the covering bound; Montgomery–Vaughan; M. Riesz and Fujii; the Jutila/Conrey density theorem; the 6a–6e class bounds; Riemann–von Mangoldt; Mathlib's exp/log facts), with Step 2 (#off-line or multiple zeros ≤ 2𝒩(T)), Steps 0 and 8, 6e's Hölder bookkeeping, 6f's class sum and Φ(ε) → 0, 5a's θ ≤ 1/2 step, Step 7's contradiction and Theorem U ⇒ D(λ) > 0 proved; sorry-free, standard axioms, 1,459 lines, every compile under the 290 MB watchdog (≤ 286 MiB; final run ≤ 279 MiB). Reviewer: FAITHFUL. [LEAN.md](Riemann/Deep/DPCUncond/LEAN.md)

Orchestrator (12 Nov 2026): accepted from the build logs (copied from .scratch/dpclean, now deleted, to [Lean/audit/](Lean/audit/): verify.sh, verify.log/out, per-module logs in audit/logs/, the reviewer prompt and verdict; all 13 modules rc = 0, true peaks ≤ 279 MiB, watchdog peaks ≤ 268 MB) and a static scan (no sorry, admit, native_decide, skipKernelTC or axiom declarations; 28 #print axioms lines, standard axioms only). The oleans in .scratch/dpclean/olean were build output (deleted); audit/verify.sh rebuilds them.
