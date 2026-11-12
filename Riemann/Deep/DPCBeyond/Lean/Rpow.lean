import Riemann.Deep.DPCUncond.Lean.Tail

/-!
# dpcbeyond: abstract real powers from the abstract exp/log

rp ex lg x e := ex (e · lg x) for x > 0, 0 for x ≤ 0: Real.rpow on [0, ∞) for e ≠ 0 (Real.rpow_def_of_pos,
Real.zero_rpow), built on dpclean's ExpLog (Riemann/Deep/DPCUncond/Lean/Tail.lean; Mathlib's Real.log and
Real.rpow exceed the 290 MB cap, Riemann/Deep/OddRobin/LEAN.md). Lemmas: values, x^{e+f}, (xy)^e, (x^e)^f,
monotonicity in the base (e ≥ 0) and in the exponent (x ≥ 1), x^e ≤ max 1 x for e ∈ [0, 1].
-/

namespace DPCBLean

open DPCLean

variable {ex lg : ℝ → ℝ}

/-- x^e on [0, ∞) (e ≠ 0): ex (e · lg x) for x > 0, 0 for x ≤ 0. -/
noncomputable def rp (ex lg : ℝ → ℝ) (x e : ℝ) : ℝ := if 0 < x then ex (e * lg x) else 0

lemma lgex (h : ExpLog ex lg) (x : ℝ) : lg (ex x) = x :=
  h.mono.injective (h.ex_lg _ (h.pos x))

lemma lg_lt (h : ExpLog ex lg) {x y : ℝ} (hx : 0 < x) (hxy : x < y) : lg x < lg y := by
  rcases lt_or_ge (lg x) (lg y) with h1 | h1
  · exact h1
  · have h2 := h.mono.monotone h1
    rw [h.ex_lg _ hx, h.ex_lg _ (lt_trans hx hxy)] at h2
    linarith

lemma lg_le (h : ExpLog ex lg) {x y : ℝ} (hx : 0 < x) (hxy : x ≤ y) : lg x ≤ lg y := by
  rcases hxy.lt_or_eq with h1 | h1
  · exact (lg_lt h hx h1).le
  · rw [h1]

lemma lg_nonneg (h : ExpLog ex lg) {x : ℝ} (hx : 1 ≤ x) : 0 ≤ lg x := by
  have h1 := lg_le h one_pos hx
  rwa [h.lg_one] at h1

lemma lg_nonpos (h : ExpLog ex lg) {x : ℝ} (hx0 : 0 < x) (hx : x ≤ 1) : lg x ≤ 0 := by
  have h1 := lg_le h hx0 hx
  rwa [h.lg_one] at h1

lemma rp_of_pos {x : ℝ} (hx : 0 < x) (e : ℝ) : rp ex lg x e = ex (e * lg x) := ite_eq_left hx

lemma rp_of_nonpos {x : ℝ} (hx : x ≤ 0) (e : ℝ) : rp ex lg x e = 0 := ite_eq_right (not_lt.2 hx)

lemma rp_pos (h : ExpLog ex lg) {x : ℝ} (hx : 0 < x) (e : ℝ) : 0 < rp ex lg x e := by
  rw [rp_of_pos hx]
  exact h.pos _

lemma rp_nonneg (h : ExpLog ex lg) (x e : ℝ) : 0 ≤ rp ex lg x e := by
  rcases lt_or_ge 0 x with hx | hx
  · exact (rp_pos h hx e).le
  · rw [rp_of_nonpos hx]

lemma rp_one (h : ExpLog ex lg) {x : ℝ} (hx : 0 < x) : rp ex lg x 1 = x := by
  rw [rp_of_pos hx, one_mul, h.ex_lg x hx]

lemma rp_zero (h : ExpLog ex lg) {x : ℝ} (hx : 0 < x) : rp ex lg x 0 = 1 := by
  rw [rp_of_pos hx, zero_mul, h.zero]

lemma rp_add (h : ExpLog ex lg) {x : ℝ} (hx : 0 < x) (e f : ℝ) :
    rp ex lg x (e + f) = rp ex lg x e * rp ex lg x f := by
  rw [rp_of_pos hx, rp_of_pos hx, rp_of_pos hx, add_mul, h.add]

lemma rp_mul (h : ExpLog ex lg) {x y : ℝ} (hx : 0 < x) (hy : 0 < y) (e : ℝ) :
    rp ex lg (x * y) e = rp ex lg x e * rp ex lg y e := by
  rw [rp_of_pos (mul_pos hx hy), rp_of_pos hx, rp_of_pos hy, h.lg_mul hx hy, mul_add, h.add]

lemma rp_rp (h : ExpLog ex lg) {x : ℝ} (hx : 0 < x) (e f : ℝ) :
    rp ex lg (rp ex lg x e) f = rp ex lg x (e * f) := by
  rw [rp_of_pos (rp_pos h hx e), rp_of_pos hx, rp_of_pos hx, lgex h]
  congr 1
  ring

lemma rp_neg (h : ExpLog ex lg) {x : ℝ} (hx : 0 < x) (e : ℝ) :
    rp ex lg x (-e) * rp ex lg x e = 1 := by
  rw [← rp_add h hx, neg_add_cancel, rp_zero h hx]

/-- Monotone in the base for e ≥ 0 (every real x ≤ y). -/
lemma rp_le_rp (h : ExpLog ex lg) {x y e : ℝ} (hxy : x ≤ y) (he : 0 ≤ e) :
    rp ex lg x e ≤ rp ex lg y e := by
  rcases lt_or_ge 0 x with hx | hx
  · rw [rp_of_pos hx, rp_of_pos (lt_of_lt_of_le hx hxy)]
    exact h.mono.monotone (mul_le_mul_of_nonneg_left (lg_le h hx hxy) he)
  · rw [rp_of_nonpos hx]
    exact rp_nonneg h _ _

/-- Strictly monotone in the base for e > 0 (0 < x < y). -/
lemma rp_lt_rp (h : ExpLog ex lg) {x y e : ℝ} (hx : 0 < x) (hxy : x < y) (he : 0 < e) :
    rp ex lg x e < rp ex lg y e := by
  rw [rp_of_pos hx, rp_of_pos (lt_trans hx hxy)]
  exact h.mono (mul_lt_mul_of_pos_left (lg_lt h hx hxy) he)

/-- Monotone in the exponent for x ≥ 1. -/
lemma rp_le_rp_exp (h : ExpLog ex lg) {x e f : ℝ} (hx : 1 ≤ x) (hef : e ≤ f) :
    rp ex lg x e ≤ rp ex lg x f := by
  have hx0 : 0 < x := by linarith
  rw [rp_of_pos hx0, rp_of_pos hx0]
  exact h.mono.monotone (mul_le_mul_of_nonneg_right hef (lg_nonneg h hx))

/-- x^e ≤ max 1 x for 0 ≤ e ≤ 1. -/
lemma rp_le_max (h : ExpLog ex lg) {x e : ℝ} (he0 : 0 ≤ e) (he1 : e ≤ 1) :
    rp ex lg x e ≤ max 1 x := by
  rcases lt_or_ge 0 x with hx | hx
  · rcases le_total x 1 with h1 | h1
    · rw [rp_of_pos hx]
      have h2 : e * lg x ≤ 0 := by nlinarith [lg_nonpos h hx h1]
      calc ex (e * lg x) ≤ ex 0 := h.mono.monotone h2
        _ = 1 := h.zero
        _ ≤ max 1 x := le_max_left _ _
    · calc rp ex lg x e ≤ rp ex lg x 1 := rp_le_rp_exp h h1 he1
        _ = x := rp_one h hx
        _ ≤ max 1 x := le_max_right _ _
  · rw [rp_of_nonpos hx]
    exact le_trans zero_le_one (le_max_left _ _)

end DPCBLean

