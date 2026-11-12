import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Positivity

/-!
# dpcbeyond (Theorem U′): Steps 0 and 8 and the order of choices (new step (a))

PROOF.md (Riemann/Deep/DPCBeyond), "What changes" and "Order of choices": for κ > 0 and
φ ∈ (max{2/(2+κ), 1/2}, 1) the θ-interval (1/φ − 1, min{κ/2, 1}) is nonempty and every θ in it has
0 < φθ < 1 < φ(1+θ) < 2 (DPC's sign condition, 2610.28527 l.753); q ∈ (2, κ/θ), s with 2/q + 1/s < 1,
u ∈ (1, κ/(2θ)), w > u/(u − 1) with κ(1 − 1/w) > 2θ and 1/v = 1 − 1/u − 1/w > 0 exist; the interval is empty
for 0 < φ ≤ max{2/(2+κ), 1/2}. Conrey's κ = 8/7 − δ: every φ > 7/11 admits some δ ∈ (0, 8/7), no φ ≤ 7/11 does.
-/

namespace DPCBLean

/-- θ > 1/φ − 1 ⟺ 1 < (1 + θ)φ (φ > 0). -/
lemma one_lt_of {φ θ : ℝ} (hφ0 : 0 < φ) (h : 1 / φ - 1 < θ) : 1 < (1 + θ) * φ := by
  have h' : 1 / φ < 1 + θ := by linarith
  rwa [div_lt_iff₀ hφ0] at h'

/-- Step 0′: the θ-interval (1/φ − 1, min{κ/2, 1}) is nonempty above the threshold. -/
theorem step0B_nonempty {κ φ : ℝ} (hκ : 0 < κ) (h1 : max (2 / (2 + κ)) (1 / 2) < φ) :
    ∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < min (κ / 2) 1 := by
  obtain ⟨ha, hb⟩ := max_lt_iff.1 h1
  have hφ0 : 0 < φ := by linarith
  have hA : 2 < φ * (2 + κ) := (div_lt_iff₀ (by linarith)).1 ha
  have hB : 1 / φ < 1 + κ / 2 := by rw [div_lt_iff₀ hφ0]; nlinarith
  have hC : 1 / φ < 2 := by rw [div_lt_iff₀ hφ0]; nlinarith
  have hD : 1 / φ - 1 < min (κ / 2) 1 := lt_min (by linarith) (by linarith)
  exact ⟨(1 / φ - 1 + min (κ / 2) 1) / 2, by linarith, by linarith⟩

/-- Step 0′ signs: every θ in the interval has θ ∈ (0, min{κ/2, 1}) and 0 < φθ < 1 < φ(1+θ) < 2. -/
theorem step0B_signs {κ φ θ : ℝ} (h1 : max (2 / (2 + κ)) (1 / 2) < φ) (h2 : φ < 1)
    (h3 : 1 / φ - 1 < θ) (h4 : θ < min (κ / 2) 1) :
    0 < θ ∧ θ < κ / 2 ∧ θ < 1 ∧ 0 < φ * θ ∧ φ * θ < 1 ∧ 1 < φ * (1 + θ) ∧ φ * (1 + θ) < 2 := by
  have hφ0 : 0 < φ := by have := le_max_right (2 / (2 + κ)) (1 / 2); linarith
  have h5 := one_lt_of hφ0 h3
  have h6 : θ < κ / 2 := lt_of_lt_of_le h4 (min_le_left _ _)
  have h7 : θ < 1 := lt_of_lt_of_le h4 (min_le_right _ _)
  have h8 : 0 < θ := by nlinarith
  have h9 : φ * θ < 1 := by nlinarith
  exact ⟨h8, h6, h7, mul_pos hφ0 h8, h9, by nlinarith, by nlinarith⟩

/-- Step 8′: the θ-interval is empty for 0 < φ ≤ max{2/(2+κ), 1/2}. -/
theorem step8B_empty {κ φ : ℝ} (hκ : 0 < κ) (h0 : 0 < φ) (h : φ ≤ max (2 / (2 + κ)) (1 / 2)) :
    ¬ ∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < min (κ / 2) 1 := by
  rintro ⟨θ, h1, h2⟩
  have h5 := one_lt_of h0 h1
  have h3 : θ < κ / 2 := lt_of_lt_of_le h2 (min_le_left _ _)
  have h4 : θ < 1 := lt_of_lt_of_le h2 (min_le_right _ _)
  rcases le_max_iff.1 h with h6 | h6
  · have h7 : φ * (2 + κ) ≤ 2 := (le_div_iff₀ (by linarith)).1 h6
    nlinarith [mul_lt_mul_of_pos_right h3 h0]
  · nlinarith [mul_lt_mul_of_pos_right h4 h0]

/-- Steps 0′ and 8′: the θ-interval is nonempty iff φ > max{2/(2+κ), 1/2} (φ > 0). -/
theorem step0B_iff {κ φ : ℝ} (hκ : 0 < κ) (h0 : 0 < φ) :
    (∃ θ : ℝ, 1 / φ - 1 < θ ∧ θ < min (κ / 2) 1) ↔ max (2 / (2 + κ)) (1 / 2) < φ := by
  constructor
  · intro hex
    by_contra hn
    exact step8B_empty hκ h0 (not_lt.1 hn) hex
  · exact step0B_nonempty hκ

/-- Order of choices (Step 7′, M3): q ∈ (2, κ/θ) and s > 1 with 2/q + 1/s < 1, for 0 < θ < κ/2. -/
theorem exists_q_s {κ θ : ℝ} (hθ0 : 0 < θ) (hθ : θ < κ / 2) :
    ∃ q s : ℝ, 2 < q ∧ q < κ / θ ∧ 1 < s ∧ 2 / q + 1 / s < 1 := by
  have hk : 2 < κ / θ := by rw [lt_div_iff₀ hθ0]; linarith
  obtain ⟨q, hq⟩ : ∃ q : ℝ, q = (2 + κ / θ) / 2 := ⟨_, rfl⟩
  have hq2 : 2 < q := by rw [hq]; linarith
  have hqk : q < κ / θ := by rw [hq]; linarith
  have hq0 : 0 < q := by linarith
  have ht0 : 0 < 2 / q := by positivity
  have ht1 : 2 / q < 1 := (div_lt_one hq0).2 hq2
  have hx0 : 0 < 1 - 2 / q := by linarith
  refine ⟨q, 2 / (1 - 2 / q), hq2, hqk, (one_lt_div hx0).2 (by linarith), ?_⟩
  rw [one_div_div]
  linarith

/-- Order of choices (6e′ (d), (e)): u ∈ (1, κ/(2θ)) and w > u/(u − 1); then w > 1, 2uθ < κ,
κ(1 − 1/w) > 2θ and 1/v := 1 − 1/u − 1/w > 0, for 0 < θ < κ/2. -/
theorem exists_u_w {κ θ : ℝ} (hθ0 : 0 < θ) (hθ : θ < κ / 2) :
    ∃ u w : ℝ, 1 < u ∧ u < κ / (2 * θ) ∧ u / (u - 1) < w ∧ 1 < w ∧ 2 * u * θ < κ ∧
      2 * θ < κ * (1 - 1 / w) ∧ 0 < 1 - 1 / u - 1 / w := by
  have hκ : 0 < κ := by linarith
  have h2θ : (0 : ℝ) < 2 * θ := by linarith
  have hk : 1 < κ / (2 * θ) := by rw [lt_div_iff₀ h2θ]; linarith
  obtain ⟨u, hu⟩ : ∃ u : ℝ, u = (1 + κ / (2 * θ)) / 2 := ⟨_, rfl⟩
  have hu1 : 1 < u := by rw [hu]; linarith
  have hu2 : u < κ / (2 * θ) := by rw [hu]; linarith
  have hu0 : 0 < u := by linarith
  have hu10 : 0 < u - 1 := by linarith
  have hr0 : 0 < u / (u - 1) := div_pos hu0 hu10
  have hr1 : 1 < u / (u - 1) := (one_lt_div hu10).2 (by linarith)
  have huθ : 2 * u * θ < κ := by
    have h := (lt_div_iff₀ h2θ).1 hu2
    nlinarith
  refine ⟨u, u / (u - 1) + 1, hu1, hu2, by linarith, by linarith, huθ, ?_, ?_⟩
  · have hw : 1 / (u / (u - 1) + 1) < 1 / (u / (u - 1)) :=
      one_div_lt_one_div_of_lt hr0 (by linarith)
    rw [one_div_div, sub_div, div_self hu0.ne'] at hw
    have h3 : κ * (1 / u) < κ * (1 - 1 / (u / (u - 1) + 1)) :=
      mul_lt_mul_of_pos_left (by linarith) hκ
    have h4 : 2 * θ < κ / u := by rw [lt_div_iff₀ hu0]; nlinarith
    rw [mul_one_div] at h3
    linarith
  · have hw : 1 / (u / (u - 1) + 1) < 1 / (u / (u - 1)) :=
      one_div_lt_one_div_of_lt hr0 (by linarith)
    rw [one_div_div, sub_div, div_self hu0.ne'] at hw
    linarith

/-- **Order of choices** (PROOF.md): every φ ∈ (max{2/(2+κ), 1/2}, 1) admits θ, q, s, u, w meeting all
constraints of Steps 0′, 6e′ and 7′. -/
theorem order_of_choices {κ φ : ℝ} (hκ : 0 < κ) (h1 : max (2 / (2 + κ)) (1 / 2) < φ) (h2 : φ < 1) :
    ∃ θ q s u w : ℝ, 1 / φ - 1 < θ ∧ θ < min (κ / 2) 1 ∧ 0 < φ * θ ∧ φ * θ < 1 ∧
      1 < φ * (1 + θ) ∧ φ * (1 + θ) < 2 ∧ 2 < q ∧ q < κ / θ ∧ 1 < s ∧ 2 / q + 1 / s < 1 ∧
      1 < u ∧ u < κ / (2 * θ) ∧ u / (u - 1) < w ∧ 2 * θ < κ * (1 - 1 / w) ∧ 0 < 1 - 1 / u - 1 / w := by
  obtain ⟨θ, h3, h4⟩ := step0B_nonempty hκ h1
  obtain ⟨hθ0, hθκ, -, s1, s2, s3, s4⟩ := step0B_signs h1 h2 h3 h4
  obtain ⟨q, s, q1, q2, q3, q4⟩ := exists_q_s hθ0 hθκ
  obtain ⟨u, w, u1, u2, u3, -, -, u6, u7⟩ := exists_u_w hθ0 hθκ
  exact ⟨θ, q, s, u, w, h3, h4, s1, s2, s3, s4, q1, q2, q3, q4, u1, u2, u3, u6, u7⟩

/-- Conrey's slope κ = 8/7 − δ: every φ > 7/11 is above the threshold for some δ ∈ (0, 8/7). -/
theorem conrey_threshold {φ : ℝ} (h1 : 7 / 11 < φ) :
    ∃ δ : ℝ, 0 < δ ∧ δ < 8 / 7 ∧ max (2 / (2 + (8 / 7 - δ))) (1 / 2) < φ := by
  have hφ0 : 0 < φ := by linarith
  have hg : 0 < 22 / 7 * φ - 2 := by linarith
  obtain ⟨δ, hδ⟩ : ∃ δ : ℝ, δ = min 1 ((22 / 7 * φ - 2) / (2 * φ)) := ⟨_, rfl⟩
  have hd0 : 0 < δ := by rw [hδ]; exact lt_min one_pos (by positivity)
  have hd1 : δ ≤ 1 := by rw [hδ]; exact min_le_left _ _
  have hd2 : δ * φ ≤ (22 / 7 * φ - 2) / 2 := by
    have h := mul_le_mul_of_nonneg_right (show δ ≤ (22 / 7 * φ - 2) / (2 * φ) by
      rw [hδ]; exact min_le_right _ _) hφ0.le
    have e : (22 / 7 * φ - 2) / (2 * φ) * φ = (22 / 7 * φ - 2) / 2 := by
      field_simp
    linarith
  refine ⟨δ, hd0, by linarith, max_lt ?_ (by linarith)⟩
  rw [div_lt_iff₀ (by linarith)]
  nlinarith

/-- No φ ≤ 7/11 is above the threshold for any δ ∈ [0, 8/7): 2/(2 + 8/7 − δ) ≥ 7/11. -/
theorem conrey_sharp {φ δ : ℝ} (h1 : φ ≤ 7 / 11) (hδ0 : 0 ≤ δ) (hδ1 : δ < 8 / 7) :
    ¬ max (2 / (2 + (8 / 7 - δ))) (1 / 2) < φ := by
  intro h
  have h2 : 2 / (2 + (8 / 7 - δ)) < φ := lt_of_le_of_lt (le_max_left _ _) h
  rw [div_lt_iff₀ (by linarith)] at h2
  nlinarith [mul_nonneg (sub_nonneg.2 h1) (show (0 : ℝ) ≤ 22 / 7 - δ by linarith)]

end DPCBLean

