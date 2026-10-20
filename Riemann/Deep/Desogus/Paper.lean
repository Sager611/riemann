import Riemann.Deep.Weil.Defs

/-!
# Desogus, arXiv 2609.20367v3 ("The Three Gates"), laid down and pinned

Paper: localized odd Weil operator `A_{a,-} = C_a - 2|sinh(y/2)⟩⟨sinh(y/2)|` on `L²_odd(-a,a)`;
its quadratic form is `weilQ` (`Riemann/Deep/Weil/Defs.lean`) on real odd `f`; window `L = 2a`.
Integer endpoints `a_k = log k / 2`, i.e. window `log k`.

Chain (all proved below from the cited Props):
* Thm 1.2 (restricted odd Weil criterion; correct, standard)          — `OddCriterion`
* Lemma 1.3/1.4 (zero extension, cofinal endpoints; correct)           — `WeilPosOdd.mono`, `cofinal`
* Certificate 6.64 (`A_7 ≻ 0`, window log 7; plausible, < Mori log 15) — `Base7`
* Thm 8.6 / Cor 8.8 induction step `A_k ≻ 0 ⇒ A_{k+1} ≻ 0` (k ≥ 7)   — `Step`
  ⇒ `desogus_RH`.

Failing step: Thm 8.6's quantitative core (201) with (MASTER-SCALAR+) `g_k ≥ 1.878`:
`⟨H_k f, A_{k+1} H_k f⟩ ≥ g_k |⟨e_k,f⟩|²`. Because `H_k f` minimises the form over all old-core
completions, (201) gives `Q(v) ≥ g_k ⟨e_k,v⟩²` for EVERY odd `v` on window `log (k+1)` (`Thm86Quant`).
Rigorous Arb evaluation (`.scratch/desogus/rigorous.py`) of an explicit odd polynomial `v` on window
`log 8` gives `Q(v)/⟨e_7,v⟩² ≤ 5.7e-7` (N=40), `1.2e-10` (N=80): `Witness7` ⇒ `¬ Thm86Quant 7 g`.
The step's conclusion is not refuted (it holds if RH holds); its proof is.
-/

namespace Riemann.Deep.Desogus

open Riemann.Deep.Weil MeasureTheory

/-- Real-valued odd functions (the paper's test core `D^odd_R`). -/
def OddReal (f : ℝ → ℂ) : Prop := (∀ x, f (-x) = -f x) ∧ ∀ x, (f x).im = 0

/-- `A_{a,-} ⪰ 0` with `L = 2a`: Weil positivity on real odd tests supported in `[-L/2, L/2]`. -/
def WeilPosOdd (L : ℝ) : Prop :=
  ∀ f : ℝ → ℂ, TestFn (L / 2) f → OddReal f → 0 ≤ (weilQ f).re

/-- `A_{a,-} ≻ 0` (strict on non-zero tests). -/
def WeilPosOddStrict (L : ℝ) : Prop :=
  ∀ f : ℝ → ℂ, TestFn (L / 2) f → OddReal f → f ≠ 0 → 0 < (weilQ f).re

theorem weilQ_zero : weilQ (0 : ℝ → ℂ) = 0 := by
  simp [weilQ, poleTerm, archTerm, primeTerm, ft, autocorr]

theorem WeilPosOddStrict.weak {L : ℝ} (h : WeilPosOddStrict L) : WeilPosOdd L := by
  intro f hf ho
  by_cases h0 : f = 0
  · subst h0; rw [weilQ_zero]; simp
  · exact (h f hf ho h0).le

theorem WeilPos.odd {L : ℝ} (h : WeilPos L) : WeilPosOdd L := fun f hf _ => h f hf

/-- Lemma 1.3 (exact zero-extension compression) in form language. -/
theorem WeilPosOdd.mono {L L' : ℝ} (h : L' ≤ L) (hL : WeilPosOdd L) : WeilPosOdd L' :=
  fun f hf ho => hL f (hf.mono (by linarith)) ho

/-- Lemma 1.4 (cofinal-support reduction) at the integer endpoints `log N`, `N ≥ 7`. -/
theorem cofinal (h : ∀ N : ℕ, 7 ≤ N → WeilPosOdd (Real.log N)) (L : ℝ) : WeilPosOdd L := by
  set N : ℕ := max 7 ⌈Real.exp L⌉₊
  have hN : 7 ≤ N := le_max_left _ _
  have hpos : (0 : ℝ) < N := by exact_mod_cast (show 0 < N by omega)
  have hexp : Real.exp L ≤ N :=
    (Nat.le_ceil _).trans (by exact_mod_cast le_max_right 7 ⌈Real.exp L⌉₊)
  exact (h N hN).mono ((Real.le_log_iff_exp_le hpos).mpr hexp)

/-- Theorem 1.2 (restricted odd Weil criterion). Correct (paper's proof checked); cited. -/
def OddCriterion : Prop := (∀ L : ℝ, WeilPosOdd L) → RiemannHypothesis

/-- Certificate 6.64: `A_{a_7,-} ≻ 0`, `a_7 = log 7 / 2` (Arb, LDL pivots ≥ 9.16e-25). -/
def Base7 : Prop := WeilPosOddStrict (Real.log 7)

/-- Corollary 8.8 induction step, via Theorem 8.6: `A_k ≻ 0 ⇒ A_{k+1} ≻ 0`. -/
def Step (k : ℕ) : Prop := WeilPosOddStrict (Real.log k) → WeilPosOddStrict (Real.log (k + 1 : ℕ))

theorem all_endpoints (hb : Base7) (hs : ∀ k : ℕ, 7 ≤ k → Step k) :
    ∀ N : ℕ, 7 ≤ N → WeilPosOddStrict (Real.log N) := by
  intro N hN
  induction N, hN using Nat.le_induction with
  | base => exact_mod_cast (hb : WeilPosOddStrict (Real.log 7))
  | succ k hk ih => exact hs k hk ih

/-- The paper's logic, machine-checked: base + step + odd criterion ⇒ RH. -/
theorem desogus_RH (hb : Base7) (hs : ∀ k : ℕ, 7 ≤ k → Step k) (hc : OddCriterion) :
    RiemannHypothesis :=
  hc (cofinal fun N hN => (all_endpoints hb hs N hN).weak)

/-! ## The failing step: quantitative content of Theorem 8.6 -/

/-- Collar ground profile `sgn(y) e^{-|y|/2}` on `a_k < |y| < a_{k+1}` (paper §6.5: `q_Y(t) = t^{-1/2}`
in Mellin coordinates; Lemmas 6.4, 6.8). Unnormalised. -/
noncomputable def ground (k : ℕ) (y : ℝ) : ℝ :=
  if Real.log k / 2 < |y| ∧ |y| < Real.log (k + 1 : ℕ) / 2 then
    (if 0 < y then 1 else -1) * Real.exp (-|y| / 2) else 0

/-- `‖ground k‖² = 2 (k^{-1/2} - (k+1)^{-1/2})`. -/
noncomputable def groundNormSq (k : ℕ) : ℝ :=
  2 * (Real.exp (-(Real.log k / 2)) - Real.exp (-(Real.log (k + 1 : ℕ) / 2)))

noncomputable def pairing (k : ℕ) (f : ℝ → ℂ) : ℝ := ∫ y, ground k y * (f y).re

/-- (201) + (MASTER-SCALAR+) + harmonic-graph minimality: `Q(v) ≥ g ⟨e_k, v⟩²`, `e_k` unit. -/
def Thm86Quant (k : ℕ) (g : ℝ) : Prop :=
  ∀ f : ℝ → ℂ, TestFn (Real.log (k + 1 : ℕ) / 2) f → OddReal f →
    g * pairing k f ^ 2 ≤ groundNormSq k * (weilQ f).re

/-- Rigorous Arb computation (`.scratch/desogus/rigorous.py`, explicit odd polynomial of degree 79,
`Q(v)/⟨e_7,v⟩² ∈ [1.1016027e-10 ± 3e-18]`), transferred to a smooth test by density of
`C_c^∞(-a,a)` in the form domain (dilate then mollify; the form is continuous for
`∫ log(2+|r|) |v̂|²`). Numerical fact; cited, not an axiom. -/
def Witness7 : Prop :=
  ∃ f : ℝ → ℂ, TestFn (Real.log (8 : ℕ) / 2) f ∧ OddReal f ∧
    groundNormSq 7 * (weilQ f).re < (1 / 1000000) * pairing 7 f ^ 2

theorem thm86_false (hw : Witness7) {g : ℝ} (hg : 1 / 1000000 ≤ g) : ¬ Thm86Quant 7 g := by
  rintro h
  obtain ⟨f, hf, ho, hlt⟩ := hw
  have h1 := h f (by simpa using hf) ho
  nlinarith [sq_nonneg (pairing 7 f)]

/-- Summary: the chain is valid; its quantitative step (paper: `g_7 ≥ 1.878`) is refuted. -/
theorem desogus_gap :
    (Base7 → (∀ k : ℕ, 7 ≤ k → Step k) → OddCriterion → RiemannHypothesis) ∧
      (Witness7 → ¬ Thm86Quant 7 (1878 / 1000)) :=
  ⟨desogus_RH, fun hw => thm86_false hw (by norm_num)⟩

end Riemann.Deep.Desogus
