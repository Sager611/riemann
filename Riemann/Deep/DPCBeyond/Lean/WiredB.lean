import Riemann.Deep.DPCBeyond.Lean.FineB
import Riemann.Deep.DPCBeyond.Lean.LemmaZW

/-!
# dpcbeyond: Lemma M3's first inequality wired to Lemma Z (M3InputsW), and the wired main theorems

M3InputsW: M3Inputs with its first inequality (∫|R_T|^q ≤ K_Z(∫(1+W)^{qr})^{1/r}(∫(1+On)^{q(k−1)r′})^{1−1/r}) replaced
by Lemma Z's NAMED inputs (i)–(iii) (ZPoint1 at y = 0 at every t ∈ [T, 2T]) and two NAMED integral facts (∫_T^{2T}
monotone and positively homogeneous on the two integrable functions |R_T|^q, (1 + W)^q((1 + On)^{k−1})^q; Hölder with
r, r′ = r/(r − 1)); the other five M3 inputs unchanged. m3Inputs_of_wired: M3InputsW ⇒ M3Inputs
(m3_first_of_lemmaZ). theoremUB_wired, D_pos_seven_elevenths_wired: the main results from these inputs.
-/

namespace DPCBLean

open Finset PPLean DPCLean

variable {ex lg : ℝ → ℝ}

/-- M3Inputs with the first inequality from Lemma Z's named inputs (ZPoint1 at y = 0 at every t ∈ [T, 2T]) and two
named integral facts at each T: (hcmp) ∫_T^{2T} is monotone and positively homogeneous on |R_T|^q and
(1 + W)^q((1 + On)^{k−1})^q (Rq q T = ∫|R_T|^q, IA T = ∫(1+W)^q((1+On)^{k−1})^q; both integrable on [T, 2T]);
(hHol) Hölder IA ≤ IW^{1/r} IO^{1−1/r} (IW = ∫(1+W)^{qr}, IO = ∫((1+On)^{k−1})^{qr′}, r′ = r/(r − 1)); the other five
M3 inputs as in M3Inputs. -/
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

/-- **M3's first inequality is a theorem from Lemma Z's named inputs**: M3InputsW ⇒ M3Inputs. -/
theorem m3Inputs_of_wired (h : ExpLog ex lg) {κ θ : ℝ} {D : DataB} (hin : M3InputsW ex lg κ θ D) :
    M3Inputs ex lg κ θ D := by
  intro q r w hq hqκ hr hw hrw
  obtain ⟨CI, Kpk, K₀, K₁, T₀, k, NT, X, S, d, κk, Rabs, IA, IW, IO, NOn, NOff, NF, hCI, hKpk, hk, hK₀,
    hK₁, hT⟩ := hin q r w hq hqκ hr hw hrw
  refine ⟨rp ex lg (NT * (Kpk * CI * k) + 1) q, K₀, K₁, T₀, IW, IO, NOn, NOff, NF, rp_nonneg h _ _, hK₀,
    hK₁, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · intro T hT0
    obtain ⟨hX, hκ, hz, hcmp, hHol, -⟩ := hT T hT0
    exact le_of_le_of_eq (m3_first_of_lemmaZ h hX (by linarith only [hq]) hCI hKpk hk (S T) d (κk T)
      hκ (Rabs T) hz hcmp hHol) (mul_assoc _ _ _).symm
  · intro T hT0
    obtain ⟨-, -, -, -, -, hIO, -⟩ := hT T hT0
    exact hIO
  · intro T hT0
    obtain ⟨-, -, -, -, -, -, hIW, -⟩ := hT T hT0
    exact hIW
  · intro T hT0
    obtain ⟨-, -, -, -, -, -, -, hOn, -⟩ := hT T hT0
    exact hOn
  · intro T hT0
    obtain ⟨-, -, -, -, -, -, -, -, hOff, -⟩ := hT T hT0
    exact hOff
  · intro T hT0
    obtain ⟨-, -, -, -, -, -, -, -, -, hF⟩ := hT T hT0
    exact hF

/-- The named inputs with M3 from Lemma Z (M3InputsW). -/
structure InputsWB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ φ θ : ℝ) (D : DataB) : Prop where
  step3 : Step3 D.base
  main : MainTerm D.base
  step6 : Step6ClassesB Z β ex lg κ θ D.Hp
  step1 : Step1Pos D.base
  holder : HolderB ex lg D
  m3 : M3InputsW ex lg κ θ D
  step5c : Step5cB D
  cover : Covering Z φ D.base

/-- Step 0′ with the wired inputs. -/
def DPCInputsWB (Z : Zeros) (β : ℕ → ℝ) (ex lg : ℝ → ℝ) (κ : ℝ) : Prop :=
  ∀ φ θ : ℝ, 0 < φ → φ < 1 → 0 < θ → θ < min (κ / 2) 1 → 0 < φ * θ → φ * θ < 1 →
    1 < φ * (1 + θ) → φ * (1 + θ) < 2 → ∃ D : DataB, InputsWB Z β ex lg κ φ θ D

theorem inputsFineB_of_wired (hEL : ExpLog ex lg) {Z : Zeros} {β : ℕ → ℝ} {κ : ℝ}
    (hin : DPCInputsWB Z β ex lg κ) : DPCInputsFineB Z β ex lg κ := by
  intro φ θ h1 h2 h3 h4 h5 h6 h7 h8
  obtain ⟨D, inp⟩ := hin φ θ h1 h2 h3 h4 h5 h6 h7 h8
  exact ⟨D, ⟨inp.step3, inp.main, inp.step6, inp.step1, inp.holder, m3Inputs_of_wired hEL inp.m3,
    inp.step5c, inp.cover⟩⟩

/-- **Theorem U′** with M3's first inequality from Lemma Z. -/
theorem theoremUB_wired {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    {κ : ℝ} (hκ : 0 < κ) (hden : DensityK Z β ex lg κ) (hin : DPCInputsWB Z β ex lg κ) :
    TheoremUB Z κ :=
  theoremUB_fine hZ hβ hEL htp hsc hκ hden (inputsFineB_of_wired hEL hin)

/-- **D(λ) > 0 for every λ > 7/11** with M3's first inequality from Lemma Z. -/
theorem D_pos_seven_elevenths_wired {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β)
    (hEL : ExpLog ex lg) {tp : ℝ} (htp : 0 < tp) (hsc : ∀ T : ℝ, 0 < T → Z.sc T * tp = lg T)
    (hR : RvMUp Z) (hden : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DensityK Z β ex lg (8 / 7 - δ))
    (hin : ∀ δ : ℝ, 0 < δ → δ < 8 / 7 → DPCInputsWB Z β ex lg (8 / 7 - δ)) :
    ∀ lam : ℝ, 7 / 11 < lam → 0 < Dlam Z lam :=
  D_pos_seven_elevenths_fine hZ hβ hEL htp hsc hR hden (fun δ h1 h2 => inputsFineB_of_wired hEL (hin δ h1 h2))

end DPCBLean

#print axioms DPCBLean.m3Inputs_of_wired
#print axioms DPCBLean.inputsFineB_of_wired
#print axioms DPCBLean.theoremUB_wired
#print axioms DPCBLean.D_pos_seven_elevenths_wired

