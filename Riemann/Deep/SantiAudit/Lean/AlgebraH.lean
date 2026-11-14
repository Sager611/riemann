/-
E1Lean.AlgebraH — the h-algebra of E1 (E1.md §1, §4; SL §5, last paragraph): √ from its two fields,
h > 0 ⟺ k > k*, h increasing in k and in c.
-/
import Riemann.Deep.SantiAudit.Lean.Defs

namespace E1Lean

variable (F : Fns)

/-- The two √ fields determine √ on [0, ∞) (so `F.sqrt` = `Real.sqrt` there). -/
theorem sqrt_eq_of {a s : ℝ} (hs : 0 ≤ s) (h : s * s = a) : F.sqrt a = s := by
  have ha : 0 ≤ a := h ▸ mul_self_nonneg s
  exact (mul_self_inj (F.sqrt_nonneg a) hs).1 ((F.mul_self_sqrt a ha).trans h.symm)

theorem sqrt_le_sqrt' {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) : F.sqrt a ≤ F.sqrt b := by
  by_contra hc
  rw [not_le] at hc
  have := mul_self_lt_mul_self (F.sqrt_nonneg b) hc
  rw [F.mul_self_sqrt a ha, F.mul_self_sqrt b (ha.trans hab)] at this
  linarith

theorem sqrt_lt_of_lt_sq {a B : ℝ} (ha : 0 ≤ a) (hB : 0 ≤ B) (h : a < B * B) : F.sqrt a < B := by
  by_contra hc
  rw [not_lt] at hc
  have := mul_self_le_mul_self hB hc
  rw [F.mul_self_sqrt a ha] at this
  linarith

theorem le_sqrt_of_sq_le {a B : ℝ} (h : B * B ≤ a) : B ≤ F.sqrt a := by
  by_contra hc
  rw [not_le] at hc
  have ha : 0 ≤ a := (mul_self_nonneg B).trans h
  have := mul_self_lt_mul_self (F.sqrt_nonneg a) hc
  rw [F.mul_self_sqrt a ha] at this
  linarith

/-- E1.md §4, SL §1: for k < 1, c < 1: h > 0 ⟺ k > k* = 1 − 2/(2 − c). -/
theorem hc_pos_iff {c k : ℝ} (hk : k < 1) (hc1 : c < 1) : 0 < hc F c k ↔ 1 - 2 / (2 - c) < k := by
  have h2c : 0 < 2 - c := by linarith
  have hA : 0 ≤ (1 - k) * (9 - k - 8 * c) := mul_nonneg (by linarith) (by linarith)
  have key : 1 - 2 / (2 - c) < k ↔ -c < k * (2 - c) := by
    have h2c' : 2 - c ≠ 0 := h2c.ne'
    have e : 1 - 2 / (2 - c) = -c / (2 - c) := by field_simp; ring
    rw [e, div_lt_iff₀ h2c]
  rw [key]
  unfold hc
  have hs := F.sqrt_nonneg ((1 - k) * (9 - k - 8 * c))
  constructor
  · intro h
    have h3 : F.sqrt ((1 - k) * (9 - k - 8 * c)) < 3 + k := by linarith
    have h4 := mul_self_lt_mul_self hs h3
    rw [F.mul_self_sqrt _ hA] at h4
    nlinarith
  · intro h
    have h3 : 0 < 3 + k := by nlinarith
    have h4 : (1 - k) * (9 - k - 8 * c) < (3 + k) * (3 + k) := by nlinarith
    have h5 := sqrt_lt_of_lt_sq F hA h3.le h4
    linarith

/-- h increases in k (SL §5: ∂_k[(1 − k)(9 − k − 8c)] = −(10 − 2k − 8c) < 0 for c < 1). -/
theorem hc_mono_k {c k₁ k₂ : ℝ} (hc1 : c < 1) (hk : k₁ ≤ k₂) (hk2 : k₂ < 1) : hc F c k₁ ≤ hc F c k₂ := by
  unfold hc
  have hA2 : 0 ≤ (1 - k₂) * (9 - k₂ - 8 * c) := mul_nonneg (by linarith) (by linarith)
  have hm : 0 ≤ (k₂ - k₁) * (10 - 8 * c - k₁ - k₂) := mul_nonneg (by linarith) (by linarith)
  have hle : (1 - k₂) * (9 - k₂ - 8 * c) ≤ (1 - k₁) * (9 - k₁ - 8 * c) := by nlinarith
  have := sqrt_le_sqrt' F hA2 hle
  linarith

/-- h increases in c (fixed k < 1, c₁ ≤ c₂ < 1). -/
theorem hc_mono_c {c₁ c₂ k : ℝ} (hk : k < 1) (hcc : c₁ ≤ c₂) (hc2 : c₂ < 1) : hc F c₁ k ≤ hc F c₂ k := by
  unfold hc
  have hA2 : 0 ≤ (1 - k) * (9 - k - 8 * c₂) := mul_nonneg (by linarith) (by linarith)
  have hm : 0 ≤ (1 - k) * (8 * (c₂ - c₁)) := mul_nonneg (by linarith) (by linarith)
  have hle : (1 - k) * (9 - k - 8 * c₂) ≤ (1 - k) * (9 - k - 8 * c₁) := by nlinarith
  have := sqrt_le_sqrt' F hA2 hle
  linarith

end E1Lean
