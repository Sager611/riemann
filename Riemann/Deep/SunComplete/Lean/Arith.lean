import Mathlib.Data.Rat.Floor
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.FieldSimp

/-!
# sunlean · Arith: the rational inequalities behind Proposition 4 and Lemma 5 (NOTE §3)
-/

namespace SunLean

/-- NOTE's `F(a) = Ñ²·a·(n − a − C) − E²·(a + C)(n − a)`. -/
def Fq (Nt E n C a : ℚ) : ℚ := Nt ^ 2 * a * (n - a - C) - E ^ 2 * (a + C) * (n - a)

/-! ### Arithmetic lemmas (plain `ℚ`) -/

/-- Proposition 4 (i): `s(Ñ + E) ≤ n(ε/τ + E)`, multiplied by `τ`. -/
lemma indep_arith {τ E Nt ε n s e : ℚ} (hτ : 0 < τ) (hE : 0 ≤ E) (hε : 0 ≤ ε) (hs0 : 0 ≤ s)
    (hsn : s ≤ n)
    (he : e ≤ ε * s)
    (hL : τ * (Nt * s * s) - n * e ≤ 0 ∨
      (τ * (Nt * s * s) - n * e) ^ 2 ≤ (τ * E) ^ 2 * ((s * (n - s)) * (s * (n - s)))) :
    s * (τ * (Nt + E)) ≤ n * (ε + τ * E) := by
  have hR : 0 ≤ τ * E * (s * (n - s)) :=
    mul_nonneg (mul_nonneg hτ.le hE) (mul_nonneg hs0 (by linarith))
  have hL' : τ * (Nt * s * s) - n * e ≤ τ * E * (s * (n - s)) := by
    rcases hL with h | h
    · linarith
    · by_contra hlt
      push Not at hlt
      have h2 := mul_self_lt_mul_self hR hlt
      linarith [h, h2]
  have hn0 : 0 ≤ n := le_trans hs0 hsn
  have key : s * (τ * Nt * s - n * ε) ≤ s * (τ * E * (n - s)) := by
    have : n * e ≤ n * (ε * s) := mul_le_mul_of_nonneg_left he hn0
    linarith [hL', this]
  rcases hs0.lt_or_eq with hpos | h0
  · have := le_of_mul_le_mul_left key hpos
    linarith [this]
  · rw [← h0, zero_mul]
    exact mul_nonneg hn0 (add_nonneg hε (mul_nonneg hτ.le hE))

/-- Proposition 4 (ii): an empty cut `(A, B)` with `|S| = c` forces `Φ(a, c) ≤ 0`. -/
lemma cut_arith {τ E Nt n a b c e : ℚ} (hτ : 0 < τ) (hNt : 0 < Nt) (ha : 1 ≤ a) (hb : 1 ≤ b)
    (hn : n = a + b + c) (he : e = 0)
    (hL : τ * (Nt * a * b) - n * e ≤ 0 ∨
      (τ * (Nt * a * b) - n * e) ^ 2 ≤ (τ * E) ^ 2 * ((a * (n - a)) * (b * (n - b)))) :
    Nt ^ 2 * a * (n - a - c) ≤ E ^ 2 * (a + c) * (n - a) := by
  have ha0 : 0 < a := by linarith
  have hb0 : 0 < b := by linarith
  have hpos : 0 < τ * (Nt * a * b) := mul_pos hτ (mul_pos (mul_pos hNt ha0) hb0)
  rcases hL with h | h
  · rw [he] at h; linarith
  · have hP : 0 < τ ^ 2 * (a * b) := by positivity
    have h2 : τ ^ 2 * (a * b) * (Nt ^ 2 * (a * b)) ≤
        τ ^ 2 * (a * b) * (E ^ 2 * ((n - a) * (n - b))) := by
      have e1 : (τ * (Nt * a * b) - n * e) ^ 2 = τ ^ 2 * (a * b) * (Nt ^ 2 * (a * b)) := by
        rw [he]; ring
      have e2 : (τ * E) ^ 2 * ((a * (n - a)) * (b * (n - b))) =
          τ ^ 2 * (a * b) * (E ^ 2 * ((n - a) * (n - b))) := by ring
      rw [← e1, ← e2]; exact h
    have h3 := le_of_mul_le_mul_left h2 hP
    have e3 : n - a - c = b := by rw [hn]; ring
    have e4 : n - b = a + c := by rw [hn]; ring
    rw [e3]
    rw [e4] at h3
    linarith [h3]

/-- Proposition 4: the concave quadratic `F` on `[a₁, n/2]`, and the monotonicity `Φ(a,c) ≥ F(a)`. -/
lemma prop4_arith {Nt E n C a c δ₀ : ℚ} (hE : 0 ≤ E) (hNE : E < Nt) (hcC : c ≤ C)
    (ha : δ₀ ≤ a - 1 + c) (h2a : 2 * a + c ≤ n) (hc0 : 0 ≤ c) (ha1 : 1 ≤ a)
    (hF : n / 2 < δ₀ - C + 1 ∨ (0 < Fq Nt E n C (δ₀ - C + 1) ∧ 0 < Fq Nt E n C (n / 2))) :
    E ^ 2 * (a + c) * (n - a) < Nt ^ 2 * a * (n - a - c) := by
  have haa1 : δ₀ - C + 1 ≤ a := by linarith
  have ham : a ≤ n / 2 := by linarith
  have hna : 0 ≤ n - a := by linarith
  have hmono : Fq Nt E n C a ≤ Nt ^ 2 * a * (n - a - c) - E ^ 2 * (a + c) * (n - a) := by
    have e : Nt ^ 2 * a * (n - a - c) - E ^ 2 * (a + c) * (n - a) - Fq Nt E n C a =
        (C - c) * (Nt ^ 2 * a + E ^ 2 * (n - a)) := by unfold Fq; ring
    have h1 : 0 ≤ (C - c) * (Nt ^ 2 * a + E ^ 2 * (n - a)) :=
      mul_nonneg (by linarith)
        (add_nonneg (mul_nonneg (sq_nonneg _) (by linarith)) (mul_nonneg (sq_nonneg _) hna))
    linarith
  have hFa : 0 < Fq Nt E n C a := by
    rcases hF with h | ⟨h1, h2⟩
    · linarith
    · have hκ0 : 0 < (Nt - E) * (Nt + E) := mul_pos (by linarith) (by linarith)
      have hκ : 0 < Nt ^ 2 - E ^ 2 := by linarith [hκ0]
      have hid : ∀ x y : ℚ, Fq Nt E n C x - Fq Nt E n C y =
          (Nt ^ 2 - E ^ 2) * ((x - y) * (n - C - x - y)) := by
        intro x y; unfold Fq; ring
      by_cases hs : a + (δ₀ - C + 1) ≤ n - C
      · have h3 := hid a (δ₀ - C + 1)
        have h4 : 0 ≤ (Nt ^ 2 - E ^ 2) * ((a - (δ₀ - C + 1)) * (n - C - a - (δ₀ - C + 1))) :=
          mul_nonneg hκ.le (mul_nonneg (by linarith) (by linarith))
        linarith
      · push Not at hs
        have h3 := hid a (n / 2)
        have e : (a - n / 2) * (n - C - a - n / 2) = (n / 2 - a) * (a + n / 2 - (n - C)) := by ring
        rw [e] at h3
        have h4 : 0 ≤ (Nt ^ 2 - E ^ 2) * ((n / 2 - a) * (a + n / 2 - (n - C))) :=
          mul_nonneg hκ.le (mul_nonneg (by linarith) (by linarith))
        linarith
  linarith

/-- Lemma 5: (CS) gives `Φ(a, c) > 0` for every admissible `(a, c)` with `c < A'`. -/
lemma lemma5_arith {Nt E n A' Dd a c δ₀ : ℚ} (hE : 0 ≤ E) (hn : 0 < n) (hc0 : 0 ≤ c)
    (hcA : c < A') (hD : Dd ≤ δ₀ + 1) (hDA : A' < Dd) (h2A : 2 * A' < n)
    (hCS : E ^ 2 * Dd < (Dd - A') * Nt ^ 2 * (1 - 2 * A' / n))
    (ha : δ₀ ≤ a - 1 + c) (h2a : 2 * a + c ≤ n) (ha1 : 1 ≤ a) :
    E ^ 2 * (a + c) * (n - a) < Nt ^ 2 * a * (n - a - c) := by
  have hn0 : n ≠ 0 := hn.ne'
  have hCS' : E ^ 2 * Dd * n < (Dd - A') * Nt ^ 2 * (n - 2 * A') := by
    have e : (1 - 2 * A' / n) * n = n - 2 * A' := by field_simp
    have h := mul_lt_mul_of_pos_right hCS hn
    rw [mul_assoc ((Dd - A') * Nt ^ 2), e] at h
    exact h
  have hDA0 : 0 < Dd - A' := by linarith
  have hnA0 : 0 < n - 2 * A' := by linarith
  have hDc0 : 0 < Dd - c := by linarith
  have hprod : (Dd - A') * (n - 2 * A') ≤ (Dd - c) * (n - 2 * c) :=
    mul_le_mul (by linarith) (by linarith) hnA0.le hDc0.le
  have hNt2 : 0 ≤ Nt ^ 2 := sq_nonneg _
  have h1 : (Dd - A') * Nt ^ 2 * (n - 2 * A') ≤ (Dd - c) * Nt ^ 2 * (n - 2 * c) := by
    have := mul_le_mul_of_nonneg_right hprod hNt2
    linarith
  have hK : E ^ 2 * n * c < (Dd - c) * (Nt ^ 2 * (n - 2 * c) - E ^ 2 * n) := by linarith
  have hEnc : 0 ≤ E ^ 2 * n * c := mul_nonneg (mul_nonneg (sq_nonneg E) hn.le) hc0
  have hKpos : 0 < Nt ^ 2 * (n - 2 * c) - E ^ 2 * n := by
    by_contra h
    push Not at h
    linarith [mul_nonneg hDc0.le (neg_nonneg.mpr h)]
  have haD : Dd - c ≤ a := by linarith
  have hB : 0 < a * (Nt ^ 2 * (n - 2 * c) - E ^ 2 * n) - E ^ 2 * n * c := by
    have := mul_le_mul_of_nonneg_right haD hKpos.le
    linarith
  have hna : 0 < n - a := by linarith
  have h2 := mul_pos hna hB
  have h3 : (n - a) * (a * (Nt ^ 2 * (n - 2 * c) - E ^ 2 * n) - E ^ 2 * n * c) ≤
      n * (Nt ^ 2 * a * (n - a - c) - E ^ 2 * (a + c) * (n - a)) := by
    have e : n * (Nt ^ 2 * a * (n - a - c) - E ^ 2 * (a + c) * (n - a)) -
        (n - a) * (a * (Nt ^ 2 * (n - 2 * c) - E ^ 2 * n) - E ^ 2 * n * c) =
        Nt ^ 2 * a * c * (n - 2 * a) := by ring
    have : 0 ≤ Nt ^ 2 * a * c * (n - 2 * a) :=
      mul_nonneg (mul_nonneg (mul_nonneg hNt2 (by linarith)) hc0) (by linarith)
    linarith
  have h4 := lt_of_lt_of_le h2 h3
  by_contra h
  push Not at h
  linarith [mul_nonneg hn.le (sub_nonneg.mpr h)]


end SunLean
