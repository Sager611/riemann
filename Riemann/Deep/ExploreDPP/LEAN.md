FORMALIZED: Lemma 1 (weak duality), monotonicity and Theorem A (the 1/(1+θ)-discrete-sine dual: for 0 < θ ≤ 1 and λ ≤ 1/(1+θ) every f meeting Lamzouri's Thm Main2 hypotheses has J_θ(f) ≤ f^(0); hence no λ ≤ 7/11 in Lamzouri's range 1/2 < θ < 4/7 and no λ ≤ 1/2 at θ = 1) from no hypotheses, including the Fourier series of the periodized triangle; Theorem B (no λ ≤ 0.775 for 0 ≤ θ ≤ 4/7, so Thm Main2's hypotheses cannot be met for λ ≤ 0.775; θ = 1: no λ ≤ 0.595) from the named Arb-certificate hypotheses CertB, CertB1. Lemma 1, Theorem A and Theorem B (i), (ii) also hold for every decaying kernel (f, f^ = O(|x|^(−2)); every real Schwartz f, SmallGap's Lemma T class): weak_duality_decay, theoremA_decay, theoremB_decay, theoremB1_decay. Sorry-free, standard axioms only, no native_decide.

# ExploreDPP in Lean 4 + Mathlib (dpplean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile). Files in Riemann/Deep/ExploreDPP/Lean/: WeakDuality (422 lines: definitions, Poisson on hℤ, Lemma 1, monotonicity), TheoremB (209: certificates, Theorem B), TheoremA (109: Theorem A from TriangleSeries), Triangle (237: 𝓕Δ = sinc², TriangleSeries proved, Theorem A unconditional, corollaries), Literal (25: Lemma 1 with the brief's literal hypotheses), Decay (380: decaying and Schwartz kernels), Check (28: axiom prints); 1,410 lines. No sorry, admit, axiom, native_decide (rg). Certificate side conditions by decide +kernel (kernel reduction, standard axioms).

## Main statements (namespace ExploreDPP)

Conventions: 𝓕 is Mathlib's Fourier transform, 𝓕 f w = ∫ v, 𝐞(−(v w)) • f v = ∫ f(v) e^{−2πivw} dv (Real.fourier_real_eq), Lamzouri's convention (TeX l. 307–309). Real.sinc x = sin x / x, so DERIVATION's sinc(y) = sin(πy)/(πy) is Real.sinc (π * y).

    structure Admissible (lam : ℝ) (f : ℝ → ℝ) : Prop where
      contDiff : ContDiff ℝ 2 f
      hasCompactSupport : HasCompactSupport f
      even : ∀ x, f (-x) = f x
      nonneg : ∀ x, 0 ≤ f x
      fourier_zero : 𝓕 (fun x => (f x : ℂ)) 0 = 1
      fourier_nonpos : ∀ u, lam < |u| → (𝓕 (fun x => (f x : ℂ)) u).re ≤ 0
    def J (θ : ℝ) (f : ℝ → ℝ) : ℝ := f 0 + 2 * ∫ α in (0:ℝ)..θ, f α * α
    def Main2Kernel (θ lam : ℝ) (f : ℝ → ℝ) : Prop := Admissible lam f ∧ 1 < J θ f
    def Main2Hyp (θ lam : ℝ) (f : ℝ → ℝ) : Prop :=
      1 / 2 < θ ∧ θ < 4 / 7 ∧ 0 < lam ∧ lam ≤ 1 ∧ Main2Kernel θ lam f
    def P (h : ℝ) (η : ℕ → ℝ) (α : ℝ) : ℝ :=
      (1 - h) + 2 * h * ∑' n : ℕ, η (n + 1) * Real.cos (2 * π * α * ((n + 1) * h))

    theorem weak_duality {h lam θ : ℝ} {k0 : ℕ} {η : ℕ → ℝ}
        (hh : 0 < h) (hlam : (k0 : ℝ) * h = lam) (hη : Summable η)
        (hcore : ∀ k : ℕ, 1 ≤ k → k < k0 → η k = -1) (hge : ∀ k : ℕ, k0 ≤ k → -1 ≤ η k)
        (hθ0 : 0 ≤ θ) (hθ : θ ≤ 1 / (2 * h))
        (hP : ∀ α, 0 ≤ α → α ≤ 1 / (2 * h) → (if α ≤ θ then α else 0) ≤ P h η α)
        {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ 1
    theorem weak_duality_literal ... (_hk0 : 1 ≤ k0) ... : J θ f ≤ (𝓕 (fun x => (f x : ℂ)) 0).re
        -- same hypotheses plus k0 ≥ 1, conclusion J_θ(f) ≤ f^(0) (Literal.lean)
    lemma Admissible.mono (hf : Admissible lam' f) (h : lam' ≤ lam) : Admissible lam f
    lemma J_mono (hnn : ∀ x, 0 ≤ f x) (hc : Continuous f) {θ' θ : ℝ} (h0 : 0 ≤ θ') (h : θ' ≤ θ) :
        J θ' f ≤ J θ f
    lemma Main2Kernel.mono (hf : Main2Kernel θ' lam' f) (h0 : 0 ≤ θ') (hθ : θ' ≤ θ) (hlam : lam' ≤ lam) :
        Main2Kernel θ lam f

    -- Theorem B (TheoremB.lean)
    def etaOf (L : List ℤ) (k : ℕ) : ℝ := if 1 ≤ k then ((L.getD (k - 1) 0 : ℤ) : ℝ) / 10 ^ 12 else 0
    def Pfin (h : ℝ) (L : List ℤ) (α : ℝ) : ℝ := (1 - h) + 2 * h * ∑ n ∈ Finset.range L.length,
        etaOf L (n + 1) * Real.cos (2 * π * α * ((n + 1) * h))
    def certB : List ℤ := [-1000000000000, …]   -- 235 numerators
    def CertB : Prop :=
      ∀ α : ℝ, 0 ≤ α → α ≤ 140 / 31 → (if α ≤ 4 / 7 then α else 0) ≤ Pfin (31 / 280) certB α
    theorem theoremB (hcert : CertB) {θ lam : ℝ} (hθ0 : 0 ≤ θ) (hθ : θ ≤ 4 / 7)
        (hlam : lam ≤ 31 / 40) {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ 1
    theorem not_main2Hyp (hcert : CertB) {θ lam : ℝ} (hlam : lam ≤ 31 / 40) (f : ℝ → ℝ) :
        ¬ Main2Hyp θ lam f
    theorem no_main2Kernel (hcert : CertB) {θ lam : ℝ} (hθ0 : 0 ≤ θ) (hθ : θ ≤ 4 / 7)
        (hlam : lam ≤ 31 / 40) (f : ℝ → ℝ) : ¬ Main2Kernel θ lam f
    def CertB1 : Prop :=
      ∀ α : ℝ, 0 ≤ α → α ≤ 300 / 119 → (if α ≤ 1 then α else 0) ≤ Pfin (119 / 600) certB1 α
    theorem theoremB1 (hcert : CertB1) {θ lam : ℝ} (hθ0 : 0 ≤ θ) (hθ : θ ≤ 1)
        (hlam : lam ≤ 119 / 200) {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ 1

    -- Theorem A (TheoremA.lean, Triangle.lean)
    def tri (x : ℝ) : ℝ := max 0 (1 - |x|)
    def etaA (a : ℝ) (k : ℕ) : ℝ := -(Real.sinc (π * (k * a))) ^ 2
    def TriangleSeries (a : ℝ) : Prop :=
      ∀ α : ℝ, a + 2 * a * ∑' n : ℕ, Real.sinc (π * ((n + 1) * a)) ^ 2 *
          Real.cos (2 * π * α * ((n + 1) * a)) = ∑' m : ℤ, tri (α - m / a)
    theorem theoremA {θ lam : ℝ} (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hT : TriangleSeries (1 / (1 + θ)))
        (hlam : lam ≤ 1 / (1 + θ)) {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ 1
    theorem triangleSeries {a : ℝ} (ha : 0 < a) : TriangleSeries a
    theorem theoremA' {θ lam : ℝ} (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hlam : lam ≤ 1 / (1 + θ))
        {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ 1
    theorem not_main2Hyp_of_le_7_11 {θ lam : ℝ} (hlam : lam ≤ 7 / 11) (f : ℝ → ℝ) : ¬ Main2Hyp θ lam f
    theorem theoremA_RH {lam : ℝ} (hlam : lam ≤ 1 / 2) {f : ℝ → ℝ} (hf : Admissible lam f) : J 1 f ≤ 1

    -- Decaying kernels (Decay.lean): SmallGap's Lemma T class
    structure AdmissibleDecay (lam : ℝ) (f : ℝ → ℝ) : Prop where
      continuous : Continuous f
      integrable : Integrable f
      f_decay : f =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ)))
      fourier_decay : 𝓕 (fun x => (f x : ℂ)) =O[cocompact ℝ] (fun x : ℝ => |x| ^ (-(2:ℝ)))
      even : ∀ x, f (-x) = f x
      nonneg : ∀ x, 0 ≤ f x
      fourier_zero : 𝓕 (fun x => (f x : ℂ)) 0 = 1
      fourier_nonpos : ∀ u, lam < |u| → (𝓕 (fun x => (f x : ℂ)) u).re ≤ 0
    lemma Admissible.toDecay (hf : Admissible lam f) : AdmissibleDecay lam f
    lemma AdmissibleDecay.of_schwartz (g : 𝓢(ℝ, ℂ)) (hfg : ∀ x, g x = (f x : ℂ)) (heven : ∀ x, f (-x) = f x)
        (hnn : ∀ x, 0 ≤ f x) (h0 : 𝓕 (fun x => (f x : ℂ)) 0 = 1)
        (hneg : ∀ u, lam < |u| → (𝓕 (fun x => (f x : ℂ)) u).re ≤ 0) : AdmissibleDecay lam f
    theorem weak_duality_decay (hypotheses of weak_duality) (hf : AdmissibleDecay lam f) : J θ f ≤ 1
    theorem theoremB_decay (hcert : CertB) {θ lam : ℝ} (hθ0 : 0 ≤ θ) (hθ : θ ≤ 4 / 7)
        (hlam : lam ≤ 31 / 40) {f : ℝ → ℝ} (hf : AdmissibleDecay lam f) : J θ f ≤ 1
    theorem theoremA_decay {θ lam : ℝ} (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hlam : lam ≤ 1 / (1 + θ))
        {f : ℝ → ℝ} (hf : AdmissibleDecay lam f) : J θ f ≤ 1
    theorem theoremB1_decay (hcert : CertB1) {θ lam : ℝ} (hθ0 : 0 ≤ θ) (hθ : θ ≤ 1)
        (hlam : lam ≤ 119 / 200) {f : ℝ → ℝ} (hf : AdmissibleDecay lam f) : J θ f ≤ 1

Proof route (all in Lean): Poisson on hℤ (Mathlib Real.tsum_eq_tsum_fourier_of_rpow_decay applied to f(·/h), with 𝓕 f = O(|u|^{−2}) from Real.fourier_iteratedDeriv for f ∈ C_c²), Σ_ℤ folded by evenness (tsum_int_eq_zero_add_two_mul_tsum_pnat), ∫ f P by integral_tsum_of_summable_integral_norm, P ≥ |α|1_θ on ℝ from the base period by evenness and 1/h-periodicity (round), f^ ≤ 0 at |u| = λ by continuity. Theorem A: 𝓕Δ = sinc² by the fundamental theorem of calculus, 𝓕 of a real even function is real (conjugation + x ↦ −x), Poisson with shift −aα; on [0, (1+θ)/2] only the triangles m = 0, 1 are nonzero (tri_sum_base), then a case split (SA_ge).

## Hypotheses and sources

- CertB (Theorem B): the inequality certified outside Lean by verify_dual.py (python-flint Arb, 360 bits, outward rounding, Taylor cells with |P''| ≤ M2): logs/verify_t47exact_0775.log, PASS, min of P − α·1[α ≤ 4/7] on [0, 140/31] ≥ 0.000204, 90,324 cells, δ = 1/20000. certB = the 235 numerators that verify_dual.py builds from logs/dual_t47_k7X26.json (η_k = −1 for k < 7, else round(e_k·10¹²)/10¹² clamped at −1), regenerated by .scratch/dpplean/cert_nums.py with the same code; 29 entries equal −10¹².
- CertB1 (Theorem B(ii)): logs/verify_t1_0595.log, PASS, min ≥ 0.000239, 100,841 cells; certB1 = 202 numerators from logs/cert_t1k3_X40_lam0.595.json, k0 = 3, h = 119/600.
- Checked in Lean: all certB, certB1 entries ≥ −10¹², and η_k = −1 for 1 ≤ k < k0 (decide +kernel); 7·31/280 = 31/40, 3·119/600 = 119/200, 1/(2h) = 140/31, 300/119 (norm_num). CertB, CertB1 themselves are not checked in Lean (the optional coarse-cell check was not attempted: a degree-235 trigonometric polynomial on ~10⁵ cells without native_decide is out of reach in the time box).
- No other hypotheses: TriangleSeries is proved (triangleSeries), so Theorem A and its corollaries are unconditional.

## Fidelity

Lamzouri 2611.14486, SmallGapsPaper.tex l. 392–400 (local copy .scratch/r442/lamz/src/SmallGapsPaper.tex):

    \begin{thm}\label{Thm:Main2} Let $1/2< \theta<4/7$ be fixed and suppose that there exists a constant
    $0<\lambda\leq 1$ and a nonnegative even compactly supported function $f\in C_c^2(\R)$ such that
    $\widehat{f}(0)=1$, $\widehat{f}(u)\leq 0$ for all $|u|> \lambda$ and
    \begin{equation}\label{Eq:conditionf}
        f(0)+ 2\int_0^{\theta} f(\alpha) \alpha\, d\alpha >1.
    \end{equation}
    Then we have $ D(\lambda)>0.$
    \end{thm}

l. 307–309: \widehat{f}(\xi) = \int_{-\infty}^{\infty} f(u)e^{-2\pi i\xi u }\,du.

| Lean | Source | Match |
|---|---|---|
| Admissible.contDiff, hasCompactSupport | "compactly supported f ∈ C_c²(ℝ)" l. 392 | exact |
| nonneg, even | "nonnegative even" l. 392 | exact (f real-valued, as "nonnegative" implies) |
| fourier_zero | \widehat f(0) = 1 l. 392, convention l. 307–309 = Mathlib 𝓕 | exact |
| fourier_nonpos | \widehat f(u) ≤ 0 for all \|u\| > λ l. 392 | Lean asks only Re 𝓕 f(u) ≤ 0: a weaker hypothesis, so every barrier covers all f Lamzouri allows (for even real f, 𝓕 f is real anyway: fourier_real_of_even) |
| J | (Eq:conditionf) l. 394 | exact |
| Main2Hyp | 1/2 < θ < 4/7, 0 < λ ≤ 1, f-conditions, J > 1 (l. 392–394) | exact |
| Main2Kernel | the f-part of Thm Main2 at (θ, λ) = DERIVATION's "admissible for (θ, λ)" | the brief's "Admissible θ λ f"; Lean splits it into Admissible lam f (θ-free) and Main2Kernel θ lam f |
| P | DERIVATION §1, P(α) = (1 − h) + 2hΣ_{k≥1} η_k cos(2παkh) | exact, k = n + 1 |
| weak_duality | DERIVATION §1 Lemma 1 (repaired: θ ≤ 1/(2h)) | Summable η covers finite support and absolute summability; conclusion J_θ(f) ≤ 1 = f^(0); k0 ≥ 1 dropped (more general; k0 = 0 forces λ = 0, where no f is admissible); 0 ≤ θ added (Lamzouri: θ > 1/2) |
| Admissible.mono, J_mono, Main2Kernel.mono | DERIVATION §1 monotonicity | exact for θ' ≥ 0 (for θ < 0, J_θ = J_{−θ}) |
| CertB, CertB1 | DERIVATION §3, verify_dual.py | same rational η, same P, same interval [0, 1/(2h)], inclusive branch α ≤ θ (the left cell ending at θ certifies P(θ) − θ) |
| theoremB, not_main2Hyp, no_main2Kernel | Theorem B(i) | exact (θ ∈ [0, 4/7] ⊇ Lamzouri's range) |
| theoremB1 | Theorem B(ii), θ = 1, λ ≤ 0.595 | exact, also all 0 ≤ θ ≤ 1 |
| etaA, theoremA, theoremA' | DERIVATION §2 Theorem A, η_k = −sinc²(ka), h = a = 1/(1+θ), k0 = 1 | exact |
| triangleSeries | §2 "Poisson (FT of sinc² is Δ)": aΣ_k sinc²(ka)e^{2πiαka} = Σ_m Δ(α − m/a) | exact, folded to k ≥ 1 |
| tri_sum_base, SA_ge | §2 repaired positivity: one triangle ≤ 1, two overlapping sum to 2 − 1/a ≤ 1 | proved on the base period [0, (1+θ)/2] directly: only m = 0, 1 contribute, S_a = α (α ≤ θ), θ (θ < α ≤ 1), 1 + θ − α (α > 1) |
| not_main2Hyp_of_le_7_11, theoremA_RH | Theorem A's 7/11 (θ < 4/7) and 1/2 (θ = 1) | exact |
| AdmissibleDecay, of_schwartz, theoremB_decay, theoremB1_decay, theoremA_decay | DERIVATION §1 "Schwartz f (Lemma T): same proof"; SmallGap DERIVATION §3 Lemma T (f even, f ≥ 0, Schwartz, f^ ≤ 0 on \|ξ\| ≥ λ, f^(0) > 0, J_θ(f) > f^(0)) | covers every Lemma T kernel after f ↦ f/f^(0) (not written in Lean; J and 𝓕 are linear); Lemma T's closed-set sign condition implies Lean's open one; a real Schwartz f is passed as a complex Schwartz g with real values (the lift 𝓢(ℝ, ℝ) → 𝓢(ℝ, ℂ) is not written); Admissible.toDecay: Admissible ⊆ AdmissibleDecay |

Not formalized: Lemma T's own construction of a C_c² F from a Schwartz f (not needed: the barrier is proved for Schwartz kernels directly), the scaling f ↦ f/f^(0), the probabilistic content of §2 (Palm measures, HKPV) and §4 (collapse).

Reviewer (statement fidelity, Lean statements and TeX l. 392–400 pasted inline; .scratch/dpplean/review_out.md): ISSUES, all minor and none a fidelity error: (1) k0 ≥ 1 omitted, "harmless vacuous extension", restore only for literal fidelity; (2) 0 ≤ θ added, no Lamzouri parameter lost; (3) B, B1 conditional on CertB, CertB1 (as designed); confirmed: Main2 hypotheses faithful, Re ≤ 0 only strengthens the barrier, J ≤ 1 = f^(0), π-scaled sinc, n + 1 indexing, inclusive α ≤ θ branch, θ ≤ 1/(2h) including Theorem A's equality case, 7/11 and 1/2 corollaries. Kept the more general statements; weak_duality_literal (Literal.lean) restores k0 ≥ 1 and the conclusion J_θ(f) ≤ f^(0) verbatim. Second reviewer (Decay statements vs SmallGap's Lemma T; .scratch/dpplean/review2_out.md): ISSUES, qualifications only: B covers Lemma T kernels for λ ≤ 0.775 (the barrier region, as intended), it is a barrier and not Lemma T's construction, B is conditional on CertB; confirmed: normalization f/f^(0), closed vs open sign set in the right direction, complex encoding of real Schwartz f, Admissible ⊂ AdmissibleDecay strictly.

## Build

    rm .lake                      # the placeholder file
    lake exe cache get            # 8,723 cached files decompressed, 47 s
    lake build Riemann.Deep.ExploreDPP.Lean.WeakDuality   # ✔ Built (104 s)
    lake build Riemann.Deep.ExploreDPP.Lean.TheoremB      # ✔ Built (114 s)
    lake build Riemann.Deep.ExploreDPP.Lean.TheoremA      # ✔ Built (102 s)
    lake build Riemann.Deep.ExploreDPP.Lean.Triangle      # ✔ Built (178 s)
    lake build Riemann.Deep.ExploreDPP.Lean.Literal       # ✔ Built (49 s)
    lake build Riemann.Deep.ExploreDPP.Lean.Decay         # ✔ Built (120 s)
    lake env lean Riemann/Deep/ExploreDPP/Lean/Check.lean # exit 0, output below

Each build: "Build completed successfully", no errors, no warnings in our files. Peak lean memory 420–510 MB (top). Check.lean output:

    'ExploreDPP.weak_duality' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.weak_duality_literal' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.Admissible.mono' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.J_mono' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.Main2Kernel.mono' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.theoremB' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.not_main2Hyp' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.no_main2Kernel' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.theoremB1' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.theoremA' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.triangleSeries' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.theoremA'' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.not_main2Hyp_of_le_7_11' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.theoremA_RH' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.weak_duality_decay' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.Admissible.toDecay' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.AdmissibleDecay.of_schwartz' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.theoremB_decay' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.theoremA_decay' depends on axioms: [propext, Classical.choice, Quot.sound]
    'ExploreDPP.theoremB1_decay' depends on axioms: [propext, Classical.choice, Quot.sound]

The modules are not imported by Riemann.lean (default target untouched); build them by name as above.
.lake: left as is (Mathlib cache + our oleans, 7.4 GB); the orchestrator restores the placeholder file.

## FORMALIZATIONS.md entry

- Lamzouri's small-gap criterion has a floor, machine-checked: weak duality (ExploreDPP.weak_duality: a lattice dual with P ≥ |α|·1[|α| ≤ θ] on [0, 1/(2h)] forces J_θ(f) ≤ f^(0) for every f of Thm Main2) and the discrete-sine dual (theoremA': no λ ≤ 1/(1+θ), so none ≤ 7/11 for θ < 4/7 and none ≤ 1/2 at θ = 1), sorry-free, standard axioms, Poisson summation and 𝓕Δ = sinc² in Lean; Theorem B (no λ ≤ 0.775 for θ ≤ 4/7; θ = 1: none ≤ 0.595) from the named Arb certificates CertB, CertB1; weak duality and Theorems A, B also for Schwartz kernels (SmallGap's Lemma T route). [LEAN.md](Riemann/Deep/ExploreDPP/LEAN.md)

Build logs and fidelity reviews: ../logs/lean_build_*.out, ../logs/lean_review*_out.md.
