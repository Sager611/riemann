import Riemann.Deep.QRH.Bridge

/-!
# Dirichlet L-functions right of 7/8: no Siegel zeros (rung 2)

Corollaries of OpenAI's `OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re`
(ported in Riemann/Deep/QRH/, Apache 2.0): every Dirichlet L-function has no zero with
Re s > 7/8, apart from the pole of the principal character at s = 1.

* `LFunction_real_ne_zero`: no real zero β > 7/8 (so no Siegel zeros at all).
* `zeroFree_region_log`: L(s, χ) ≠ 0 for Re s ≥ 1 - 1/(8 log q), q ≥ 3: the classical
  zero-free region, uniform in t and with no exceptional zero.
* `exists_absolute_real_zero_gap`: the exact statement of OpenAI's `SiegelZeros` comparator
  challenge (which OpenAI proved by a separate argument), here with c = log 3 / 8 and a
  two-line proof from the 7/8 theorem.
-/

namespace Fayed.QRHCor

/-- No zero of any Dirichlet L-function with Re s > 7/8 (except the pole of χ = 1 at s = 1). -/
theorem LFunction_ne_zero_of_seven_eighths_lt_re {q : ℕ} [NeZero q]
    (χ : DirichletCharacter ℂ q) {s : ℂ} (hs : 7 / 8 < s.re) (hpole : χ ≠ 1 ∨ s ≠ 1) :
    χ.LFunction s ≠ 0 :=
  OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re χ hs
    (fun h => hpole.elim (fun h1 => h1 h.1) (fun h2 => h2 h.2))

/-- **No Siegel zeros.** No Dirichlet L-function has a real zero β > 7/8. -/
theorem LFunction_real_ne_zero {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) {β : ℝ}
    (hβ : 7 / 8 < β) (hpole : χ ≠ 1 ∨ β ≠ 1) : χ.LFunction (β : ℂ) ≠ 0 :=
  LFunction_ne_zero_of_seven_eighths_lt_re χ (by simpa using hβ)
    (hpole.imp_right fun h h' => h (by exact_mod_cast h'))

/-- A real zero of a nonprincipal L-function satisfies β ≤ 7/8. -/
theorem le_seven_eighths_of_LFunction_eq_zero {q : ℕ} [NeZero q] {χ : DirichletCharacter ℂ q}
    (hχ : χ ≠ 1) {β : ℝ} (hz : χ.LFunction (β : ℂ) = 0) : β ≤ 7 / 8 := by
  by_contra h
  exact LFunction_real_ne_zero χ (not_le.mp h) (Or.inl hχ) hz

lemma one_lt_log_three : 1 < Real.log 3 := by
  rw [Real.lt_log_iff_exp_lt (by norm_num)]
  have := Real.exp_one_lt_d9
  norm_num at this ⊢
  linarith

/-- **Zero-free region of classical shape, uniform in t, no exceptional zero:**
for q ≥ 3, L(s, χ) ≠ 0 whenever Re s ≥ 1 - 1/(8 log q). -/
theorem zeroFree_region_log {q : ℕ} [NeZero q] (hq : 3 ≤ q) (χ : DirichletCharacter ℂ q)
    {s : ℂ} (hs : 1 - 1 / (8 * Real.log q) ≤ s.re) (hpole : χ ≠ 1 ∨ s ≠ 1) :
    χ.LFunction s ≠ 0 := by
  have hlog : 1 < Real.log q :=
    lt_of_lt_of_le one_lt_log_three (Real.log_le_log (by norm_num) (by exact_mod_cast hq))
  have h8 : 1 / (8 * Real.log q) < 1 / 8 := by
    apply one_div_lt_one_div_of_lt (by norm_num)
    linarith
  exact LFunction_ne_zero_of_seven_eighths_lt_re χ (by linarith) hpole

/-- OpenAI's `SiegelZeros` comparator statement (`OAI.SiegelZeros.WeightedTorusJets.
exists_absolute_real_zero_gap`), derived from the 7/8 theorem with c = log 3 / 8. The hypotheses
"primitive" and "real" are not needed. -/
theorem exists_absolute_real_zero_gap :
    ∃ c : ℝ, 0 < c ∧
      ∀ (q : ℕ) [NeZero q], 3 ≤ q →
      ∀ χ : DirichletCharacter ℂ q,
        χ.IsPrimitive → χ ≠ 1 → (∀ a : ZMod q, (χ a).im = 0) →
        ∀ β : ℝ, 0 < β → β < 1 → χ.LFunction (β : ℂ) = 0 →
          c ≤ (1 - β) * Real.log (q : ℝ) := by
  have hl3 : 0 < Real.log 3 := Real.log_pos (by norm_num)
  refine ⟨Real.log 3 / 8, div_pos hl3 (by norm_num), ?_⟩
  intro q _ hq χ _ hχ _ β _ _ hz
  have hβ := le_seven_eighths_of_LFunction_eq_zero hχ hz
  have hlog : Real.log 3 ≤ Real.log q := Real.log_le_log (by norm_num) (by exact_mod_cast hq)
  have h1 : (1 / 8 : ℝ) * Real.log q ≤ (1 - β) * Real.log q :=
    mul_le_mul_of_nonneg_right (by linarith) (by linarith)
  linarith

end Fayed.QRHCor
