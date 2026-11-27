import Mathlib.Algebra.Order.Chebyshev
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# LeanGale T6 (1/3): the quadratic-form inequality behind `SimplexPartitionBound`

For `a, b : ι → ℝ` (`N = |ι|`, `A = Σ a`, `B = Σ b`):
`psi κ a b = N⟨a,b⟩ − AB + (κ/2)(N|a|² − A² + N|b|² − B²)`.
`psi_nonneg`: `psi κ a b ≥ 0` for `a, b ≥ 0` with a common zero, `0 ≤ κ ≤ 1` and `st ≤ κ²(N − s)(N − t)` for all
naturals `s + t + 1 ≤ N`. Proof: `psi` is concave along `(e_j, −e_j)` (`psi_split`), so the mass of a coordinate
where `a` and `b` are both positive can be moved entirely into `a` or entirely into `b` (induction over the set of
such coordinates); for disjoint supports `S`, `T` Cauchy–Schwarz `A² ≤ |S|·Σa²` and the hypothesis on `|S|, |T|`
finish (`psi_nonneg_disjoint`, `base_real`).
`kappa n` = κₙ = n/(n+2) (n even), √((n−1)/(n+3)) (n odd); `kappa_int`: the hypothesis holds for κₙ and `N = n + 1`.
-/

namespace LeanGale

open Finset

section Psi

variable {ι : Type*} [Fintype ι]

/-- `Ψ_κ(a, b) = N⟨a,b⟩ − AB + (κ/2)(N|a|² − A² + N|b|² − B²)`, `N = |ι|`, `A = Σ a`, `B = Σ b` -/
noncomputable def psi (κ : ℝ) (a b : ι → ℝ) : ℝ :=
  (Fintype.card ι : ℝ) * ∑ k, a k * b k - (∑ k, a k) * ∑ k, b k
    + κ / 2 * ((Fintype.card ι : ℝ) * ∑ k, a k ^ 2 - (∑ k, a k) ^ 2
      + ((Fintype.card ι : ℝ) * ∑ k, b k ^ 2 - (∑ k, b k) ^ 2))

theorem sum_update_id [DecidableEq ι] (a : ι → ℝ) (j : ι) (v : ℝ) :
    ∑ k, Function.update a j v k = ∑ k, a k - a j + v := by
  have h : ∀ k, Function.update a j v k = a k + (if k = j then v - a j else 0) := by
    intro k
    by_cases hk : k = j
    · subst hk; simp
    · simp [hk]
  simp only [h, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  ring

theorem sum_update_sq [DecidableEq ι] (a : ι → ℝ) (j : ι) (v : ℝ) :
    ∑ k, Function.update a j v k ^ 2 = ∑ k, a k ^ 2 - a j ^ 2 + v ^ 2 := by
  have h : ∀ k, Function.update a j v k ^ 2 = a k ^ 2 + (if k = j then v ^ 2 - a j ^ 2 else 0) := by
    intro k
    by_cases hk : k = j
    · subst hk; simp
    · simp [hk]
  simp only [h, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  ring

theorem sum_update_mul [DecidableEq ι] (a b : ι → ℝ) (j : ι) (v w : ℝ) :
    ∑ k, Function.update a j v k * Function.update b j w k = ∑ k, a k * b k - a j * b j + v * w := by
  have h : ∀ k, Function.update a j v k * Function.update b j w k =
      a k * b k + (if k = j then v * w - a j * b j else 0) := by
    intro k
    by_cases hk : k = j
    · subst hk; simp
    · simp [hk]
  simp only [h, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  ring

/-- **Concavity along `(e_j, −e_j)`**: with `α = a j`, `β = b j`,
`(α + β) Ψ(a, b) = α Ψ(a + βe_j, b − βe_j) + β Ψ(a − αe_j, b + αe_j) + (α + β)(N − 1)(1 − κ)αβ`. -/
theorem psi_split [DecidableEq ι] (κ : ℝ) (a b : ι → ℝ) (j : ι) :
    (a j + b j) * psi κ a b =
      a j * psi κ (Function.update a j (a j + b j)) (Function.update b j 0)
        + b j * psi κ (Function.update a j 0) (Function.update b j (a j + b j))
        + (a j + b j) * (((Fintype.card ι : ℝ) - 1) * (1 - κ) * (a j * b j)) := by
  simp only [psi, sum_update_id, sum_update_sq, sum_update_mul]
  ring

/-- the two-block inequality: `−AB + (κ/2)(NP − A² + NR − B²) ≥ 0` from `A² ≤ sP`, `B² ≤ tR`,
`s + t + 1 ≤ N`, `st ≤ κ²(N − s)(N − t)` -/
theorem base_real {N s t κ A B P R : ℝ} (hκ : 0 ≤ κ) (hA0 : 0 ≤ A) (hB0 : 0 ≤ B) (hP : 0 ≤ P)
    (hR : 0 ≤ R) (hs : 0 ≤ s) (ht : 0 ≤ t) (hN : s + t + 1 ≤ N) (hA : A ^ 2 ≤ s * P)
    (hB : B ^ 2 ≤ t * R) (hint : s * t ≤ κ ^ 2 * (N - s) * (N - t)) :
    0 ≤ -(A * B) + κ / 2 * (N * P - A ^ 2 + (N * R - B ^ 2)) := by
  rcases hs.eq_or_lt with hs0 | hs0
  · have hA' : A = 0 := by
      rw [← hs0, zero_mul] at hA
      exact pow_eq_zero_iff (n := 2) (by norm_num) |>.mp (le_antisymm hA (sq_nonneg A))
    have h1 : B ^ 2 ≤ N * R := hB.trans (mul_le_mul_of_nonneg_right (by linarith) hR)
    have h2 : 0 ≤ N * P := mul_nonneg (by linarith) hP
    rw [hA']
    have : 0 ≤ N * P - 0 ^ 2 + (N * R - B ^ 2) := by nlinarith
    nlinarith
  rcases ht.eq_or_lt with ht0 | ht0
  · have hB' : B = 0 := by
      rw [← ht0, zero_mul] at hB
      exact pow_eq_zero_iff (n := 2) (by norm_num) |>.mp (le_antisymm hB (sq_nonneg B))
    have h1 : A ^ 2 ≤ N * P := hA.trans (mul_le_mul_of_nonneg_right (by linarith) hP)
    have h2 : 0 ≤ N * R := mul_nonneg (by linarith) hR
    rw [hB']
    have : 0 ≤ N * P - A ^ 2 + (N * R - 0 ^ 2) := by nlinarith
    nlinarith
  have hNs : 0 < N - s := by linarith
  have hNt : 0 < N - t := by linarith
  have hN0 : 0 ≤ N := by linarith
  have hp : (N - s) * A ^ 2 ≤ s * (N * P - A ^ 2) := by nlinarith
  have hr : (N - t) * B ^ 2 ≤ t * (N * R - B ^ 2) := by nlinarith
  have hp0 : 0 ≤ N * P - A ^ 2 := by
    have : 0 ≤ s * (N * P - A ^ 2) := le_trans (by positivity) hp
    exact (mul_nonneg_iff_of_pos_left hs0).mp this
  have hr0 : 0 ≤ N * R - B ^ 2 := by
    have : 0 ≤ t * (N * R - B ^ 2) := le_trans (by positivity) hr
    exact (mul_nonneg_iff_of_pos_left ht0).mp this
  set p := N * P - A ^ 2 with hpdef
  set r := N * R - B ^ 2 with hrdef
  have h1 : (N - s) * (N - t) * (A ^ 2 * B ^ 2) ≤ s * t * (p * r) := by
    calc (N - s) * (N - t) * (A ^ 2 * B ^ 2) = ((N - s) * A ^ 2) * ((N - t) * B ^ 2) := by ring
      _ ≤ (s * p) * (t * r) := mul_le_mul hp hr (by positivity) (by positivity)
      _ = s * t * (p * r) := by ring
  have h2 : s * t * (A ^ 2 * B ^ 2) ≤ s * t * (κ ^ 2 * (p * r)) := by
    calc s * t * (A ^ 2 * B ^ 2) ≤ κ ^ 2 * (N - s) * (N - t) * (A ^ 2 * B ^ 2) :=
          mul_le_mul_of_nonneg_right hint (by positivity)
      _ = κ ^ 2 * ((N - s) * (N - t) * (A ^ 2 * B ^ 2)) := by ring
      _ ≤ κ ^ 2 * (s * t * (p * r)) := mul_le_mul_of_nonneg_left h1 (by positivity)
      _ = s * t * (κ ^ 2 * (p * r)) := by ring
  have key : A ^ 2 * B ^ 2 ≤ κ ^ 2 * (p * r) := le_of_mul_le_mul_left h2 (by positivity)
  have hsq : (2 * (A * B)) ^ 2 ≤ (κ * (p + r)) ^ 2 := by
    nlinarith [mul_nonneg (sq_nonneg κ) (sq_nonneg (p - r))]
  have hle : 2 * (A * B) ≤ κ * (p + r) :=
    (pow_le_pow_iff_left₀ (by positivity) (by positivity) two_ne_zero).mp hsq
  linarith

/-- **disjoint supports** -/
theorem psi_nonneg_disjoint (κ : ℝ) (hκ : 0 ≤ κ)
    (hint : ∀ s t : ℕ, s + t + 1 ≤ Fintype.card ι →
      (s : ℝ) * t ≤ κ ^ 2 * ((Fintype.card ι : ℝ) - s) * ((Fintype.card ι : ℝ) - t))
    (i₀ : ι) (a b : ι → ℝ) (ha : ∀ k, 0 ≤ a k) (hb : ∀ k, 0 ≤ b k) (hab : ∀ k, a k = 0 ∨ b k = 0)
    (ha0 : a i₀ = 0) (hb0 : b i₀ = 0) : 0 ≤ psi κ a b := by
  classical
  set S := univ.filter (fun k => a k ≠ 0) with hS
  set T := univ.filter (fun k => b k ≠ 0) with hT
  have hX : ∑ k, a k * b k = 0 :=
    Finset.sum_eq_zero fun k _ => by rcases hab k with h | h <;> simp [h]
  have hAS : ∑ k ∈ S, a k = ∑ k, a k := Finset.sum_filter_ne_zero _
  have hPS : ∑ k ∈ S, a k ^ 2 = ∑ k, a k ^ 2 :=
    Finset.sum_filter_of_ne fun k _ h hk => h (by rw [hk]; ring)
  have hBT : ∑ k ∈ T, b k = ∑ k, b k := Finset.sum_filter_ne_zero _
  have hRT : ∑ k ∈ T, b k ^ 2 = ∑ k, b k ^ 2 :=
    Finset.sum_filter_of_ne fun k _ h hk => h (by rw [hk]; ring)
  have hA : (∑ k, a k) ^ 2 ≤ (S.card : ℝ) * ∑ k, a k ^ 2 := by
    rw [← hAS, ← hPS]; exact sq_sum_le_card_mul_sum_sq
  have hB : (∑ k, b k) ^ 2 ≤ (T.card : ℝ) * ∑ k, b k ^ 2 := by
    rw [← hBT, ← hRT]; exact sq_sum_le_card_mul_sum_sq
  have hdisj : Disjoint S T := by
    rw [hS, hT, Finset.disjoint_filter]
    intro k _ h1 h2
    rcases hab k with h | h
    · exact h1 h
    · exact h2 h
  have hsub : S ∪ T ⊆ univ.erase i₀ := by
    intro k hk
    rw [Finset.mem_erase]
    refine ⟨?_, mem_univ k⟩
    rintro rfl
    rcases Finset.mem_union.mp hk with h | h
    · exact (Finset.mem_filter.mp h).2 ha0
    · exact (Finset.mem_filter.mp h).2 hb0
  have hcard : S.card + T.card + 1 ≤ Fintype.card ι := by
    have h1 := Finset.card_le_card hsub
    rw [Finset.card_union_of_disjoint hdisj, Finset.card_erase_of_mem (mem_univ i₀),
      Finset.card_univ] at h1
    have h2 : 0 < Fintype.card ι := Fintype.card_pos_iff.mpr ⟨i₀⟩
    omega
  have hcardR : (S.card : ℝ) + T.card + 1 ≤ Fintype.card ι := by exact_mod_cast hcard
  have := base_real (N := Fintype.card ι) hκ (Finset.sum_nonneg fun k _ => ha k)
    (Finset.sum_nonneg fun k _ => hb k) (Finset.sum_nonneg fun k _ => sq_nonneg (a k))
    (Finset.sum_nonneg fun k _ => sq_nonneg (b k)) (Nat.cast_nonneg _) (Nat.cast_nonneg _) hcardR hA hB
    (hint _ _ hcard)
  unfold psi
  rw [hX]
  linarith

/-- **The quadratic-form inequality**: `Ψ_κ(a, b) ≥ 0` for `a, b ≥ 0` with a common zero `i₀`. -/
theorem psi_nonneg [DecidableEq ι] (κ : ℝ) (hκ0 : 0 ≤ κ) (hκ1 : κ ≤ 1)
    (hint : ∀ s t : ℕ, s + t + 1 ≤ Fintype.card ι →
      (s : ℝ) * t ≤ κ ^ 2 * ((Fintype.card ι : ℝ) - s) * ((Fintype.card ι : ℝ) - t))
    (i₀ : ι) (a b : ι → ℝ) (ha : ∀ k, 0 ≤ a k) (hb : ∀ k, 0 ≤ b k) (ha0 : a i₀ = 0) (hb0 : b i₀ = 0) :
    0 ≤ psi κ a b := by
  have hN : (1 : ℝ) ≤ Fintype.card ι := by
    exact_mod_cast Nat.one_le_iff_ne_zero.mpr (Fintype.card_pos_iff.mpr ⟨i₀⟩).ne'
  suffices H : ∀ F : Finset ι, ∀ a b : ι → ℝ, (∀ k, 0 ≤ a k) → (∀ k, 0 ≤ b k) →
      (∀ k ∉ F, a k = 0 ∨ b k = 0) → a i₀ = 0 → b i₀ = 0 → 0 ≤ psi κ a b from
    H univ a b ha hb (fun k hk => absurd (mem_univ k) hk) ha0 hb0
  intro F
  induction F using Finset.induction_on with
  | empty =>
    intro a b ha hb hab ha0 hb0
    exact psi_nonneg_disjoint κ hκ0 hint i₀ a b ha hb (fun k => hab k (Finset.notMem_empty k)) ha0 hb0
  | insert j F hjF ih =>
    intro a b ha hb hab ha0 hb0
    have hout : ∀ k, k ≠ j → k ∉ F → a k = 0 ∨ b k = 0 := fun k hkj hkF =>
      hab k (by simp [hkj, hkF])
    have h0 : i₀ = j → a j + b j = 0 := by rintro rfl; rw [ha0, hb0, add_zero]
    have h1 : 0 ≤ psi κ (Function.update a j (a j + b j)) (Function.update b j 0) := by
      apply ih
      · intro k
        rcases eq_or_ne k j with rfl | hk
        · rw [Function.update_self]; exact add_nonneg (ha k) (hb k)
        · rw [Function.update_of_ne hk]; exact ha k
      · intro k
        rcases eq_or_ne k j with rfl | hk
        · rw [Function.update_self]
        · rw [Function.update_of_ne hk]; exact hb k
      · intro k hk
        rcases eq_or_ne k j with rfl | hkj
        · right; rw [Function.update_self]
        · rw [Function.update_of_ne hkj, Function.update_of_ne hkj]; exact hout k hkj hk
      · rcases eq_or_ne i₀ j with h | h
        · rw [h, Function.update_self, ← h0 h]
        · rw [Function.update_of_ne h]; exact ha0
      · rcases eq_or_ne i₀ j with h | h
        · rw [h, Function.update_self]
        · rw [Function.update_of_ne h]; exact hb0
    have h2 : 0 ≤ psi κ (Function.update a j 0) (Function.update b j (a j + b j)) := by
      apply ih
      · intro k
        rcases eq_or_ne k j with rfl | hk
        · rw [Function.update_self]
        · rw [Function.update_of_ne hk]; exact ha k
      · intro k
        rcases eq_or_ne k j with rfl | hk
        · rw [Function.update_self]; exact add_nonneg (ha k) (hb k)
        · rw [Function.update_of_ne hk]; exact hb k
      · intro k hk
        rcases eq_or_ne k j with rfl | hkj
        · left; rw [Function.update_self]
        · rw [Function.update_of_ne hkj, Function.update_of_ne hkj]; exact hout k hkj hk
      · rcases eq_or_ne i₀ j with h | h
        · rw [h, Function.update_self]
        · rw [Function.update_of_ne h]; exact ha0
      · rcases eq_or_ne i₀ j with h | h
        · rw [h, Function.update_self, ← h0 h]
        · rw [Function.update_of_ne h]; exact hb0
    rcases (add_nonneg (ha j) (hb j)).eq_or_lt with hz | hpos
    · have haj : a j = 0 := by linarith [ha j, hb j]
      apply ih a b ha hb _ ha0 hb0
      intro k hk
      rcases eq_or_ne k j with rfl | hkj
      · left; exact haj
      · exact hout k hkj hk
    · have hq : 0 ≤ ((Fintype.card ι : ℝ) - 1) * (1 - κ) * (a j * b j) :=
        mul_nonneg (mul_nonneg (by linarith) (by linarith)) (mul_nonneg (ha j) (hb j))
      have : 0 ≤ (a j + b j) * psi κ a b := by
        rw [psi_split κ a b j]
        have := mul_nonneg (ha j) h1
        have := mul_nonneg (hb j) h2
        have := mul_nonneg hpos.le hq
        linarith
      exact (mul_nonneg_iff_of_pos_left hpos).mp this

/-- `Ψ` in shifted coordinates: for `Σ s = Σ t = 0`, `Ψ_κ(s_i − s, t_i − t) = N(⟨s,t⟩ + (κ/2)(|s|² + |t|²))`. -/
theorem psi_shift (κ : ℝ) (s t : ι → ℝ) (i : ι) (hs : ∑ j, s j = 0) (ht : ∑ j, t j = 0) :
    psi κ (fun j => s i - s j) (fun j => t i - t j) =
      (Fintype.card ι : ℝ) * (∑ j, s j * t j + κ / 2 * (∑ j, s j ^ 2 + ∑ j, t j ^ 2)) := by
  have h1 : ∀ u : ι → ℝ, ∑ j, (u i - u j) = (Fintype.card ι : ℝ) * u i - ∑ j, u j := by
    intro u
    rw [Finset.sum_sub_distrib]
    simp
  have h2 : ∀ u v : ι → ℝ, ∑ j, (u i - u j) * (v i - v j) =
      (Fintype.card ι : ℝ) * (u i * v i) - u i * ∑ j, v j - v i * ∑ j, u j + ∑ j, u j * v j := by
    intro u v
    have : ∀ j, (u i - u j) * (v i - v j) = u i * v i - u i * v j - v i * u j + u j * v j :=
      fun j => by ring
    simp only [this, Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum]
    simp
    ring
  have h3 : ∀ u : ι → ℝ, ∑ j, (u i - u j) ^ 2 =
      (Fintype.card ι : ℝ) * u i ^ 2 - 2 * u i * ∑ j, u j + ∑ j, u j ^ 2 := by
    intro u
    have := h2 u u
    simp only [← sq] at this
    rw [this]
    ring
  simp only [psi, h1, h2, h3, hs, ht]
  ring

end Psi

/-! ## κₙ -/

/-- κₙ = 2σₙ² − 1: n/(n+2) for even n, √((n−1)/(n+3)) for odd n -/
noncomputable def kappa (n : ℕ) : ℝ :=
  if n % 2 = 0 then (n : ℝ) / ((n : ℝ) + 2) else Real.sqrt (((n : ℝ) - 1) / ((n : ℝ) + 3))

theorem kappa_nonneg (n : ℕ) : 0 ≤ kappa n := by
  unfold kappa
  split_ifs
  · positivity
  · exact Real.sqrt_nonneg _

theorem kappa_le_one (n : ℕ) : kappa n ≤ 1 := by
  unfold kappa
  split_ifs
  · rw [div_le_one (by positivity)]
    linarith
  · calc Real.sqrt (((n : ℝ) - 1) / ((n : ℝ) + 3)) ≤ Real.sqrt 1 :=
          Real.sqrt_le_sqrt (by rw [div_le_one (by positivity)]; linarith)
      _ = 1 := Real.sqrt_one

/-- **the integer optimisation**: `st ≤ κₙ²(n + 1 − s)(n + 1 − t)` for naturals `s + t ≤ n` -/
theorem kappa_int (n s t : ℕ) (h : s + t ≤ n) :
    (s : ℝ) * t ≤ kappa n ^ 2 * ((n : ℝ) + 1 - s) * ((n : ℝ) + 1 - t) := by
  have hst : (s : ℝ) + t ≤ n := by exact_mod_cast h
  have hs0 : (0 : ℝ) ≤ s := Nat.cast_nonneg s
  have ht0 : (0 : ℝ) ≤ t := Nat.cast_nonneg t
  have hns : 0 ≤ (n : ℝ) + 1 - s := by linarith
  have hnt : 0 ≤ (n : ℝ) - s - t := by linarith
  rcases Nat.even_or_odd n with ⟨m, hm⟩ | ⟨m, hm⟩
  · have hk : kappa n = (n : ℝ) / ((n : ℝ) + 2) := by
      have h2 : n % 2 = 0 := by omega
      simp only [kappa, h2, ↓reduceIte]
    have hpos : (0 : ℝ) < ((n : ℝ) + 2) ^ 2 := by positivity
    rw [hk, div_pow, div_mul_eq_mul_div, div_mul_eq_mul_div, le_div_iff₀ hpos]
    have key : (n : ℝ) ^ 2 * ((n : ℝ) + 1 - s) * ((n : ℝ) + 1 - t) - (s : ℝ) * t * ((n : ℝ) + 2) ^ 2 =
        ((n : ℝ) + 1) * ((n : ℝ) - 2 * s) ^ 2 +
          ((n : ℝ) - s - t) * ((n : ℝ) ^ 2 * ((n : ℝ) + 1 - s) + ((n : ℝ) + 2) ^ 2 * s) := by
      ring
    have e1 : 0 ≤ ((n : ℝ) + 1) * ((n : ℝ) - 2 * s) ^ 2 := by positivity
    have e2 : 0 ≤ ((n : ℝ) - s - t) * ((n : ℝ) ^ 2 * ((n : ℝ) + 1 - s) + ((n : ℝ) + 2) ^ 2 * s) := by
      apply mul_nonneg hnt
      positivity
    linarith
  · have hn1 : (1 : ℝ) ≤ n := by
      have : 1 ≤ n := by omega
      exact_mod_cast this
    have hk : kappa n ^ 2 = ((n : ℝ) - 1) / ((n : ℝ) + 3) := by
      have h2 : ¬n % 2 = 0 := by omega
      simp only [kappa, h2, ↓reduceIte]
      rw [Real.sq_sqrt]
      apply div_nonneg <;> linarith
    have hpos : (0 : ℝ) < (n : ℝ) + 3 := by positivity
    rw [hk, div_mul_eq_mul_div, div_mul_eq_mul_div, le_div_iff₀ hpos]
    -- parity: (2s − n)² ≥ 1
    have hpar : 1 ≤ ((n : ℝ) - 2 * s) ^ 2 := by
      have hnm : (n : ℝ) = 2 * m + 1 := by rw [hm]; push_cast; ring
      have hprod : 0 ≤ ((s : ℝ) - m) * ((s : ℝ) - m - 1) := by
        rcases Nat.lt_or_ge s (m + 1) with h1 | h1
        · have : (s : ℝ) ≤ m := by exact_mod_cast Nat.lt_succ_iff.mp h1
          nlinarith
        · have : (m : ℝ) + 1 ≤ s := by exact_mod_cast h1
          nlinarith
      rw [hnm]
      nlinarith
    have key : ((n : ℝ) - 1) * ((n : ℝ) + 1 - s) * ((n : ℝ) + 1 - t) - (s : ℝ) * t * ((n : ℝ) + 3) =
        (((n : ℝ) - 2 * s) ^ 2 - 1) +
          ((n : ℝ) - s - t) * (((n : ℝ) - 1) * ((n : ℝ) + 1 - s) + ((n : ℝ) + 3) * s) := by
      ring
    have e2 : 0 ≤ ((n : ℝ) - s - t) * (((n : ℝ) - 1) * ((n : ℝ) + 1 - s) + ((n : ℝ) + 3) * s) := by
      apply mul_nonneg hnt
      have : 0 ≤ (n : ℝ) - 1 := by linarith
      positivity
    linarith

end LeanGale
