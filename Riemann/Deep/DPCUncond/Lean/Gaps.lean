import Riemann.Deep.ExplorePointProc.Lean.Model

/-!
# dpcuncond Step 2 (new step (a)): off-line and multiple zeros give zero gaps

List form (card_shared_le): in a sorted finite list x 0 ≤ … ≤ x (n − 1) the entries that share their value with
another entry number at most 2·#{i : x (i + 1) = x i} (a run of r ≥ 2 equal entries holds r − 1 zero gaps,
r ≤ 2(r − 1)). Zero model (card_exc_le, card_off_le): pplean's ordinates γ (with multiplicity, Lamzouri
2611.14486 l.177) plus real parts β with the reflection ρ ↦ 1 − ρ̄ (Reflect); then
#{off-line or multiple zeros with T ≤ γ ≤ 2T} ≤ 2𝒩(T), 𝒩(T) = #gapSet (PROOF Theorem U, Step 2).
-/

namespace DPCLean

open Finset PPLean

/-! ## (a), list form -/

/-- Entries of the sorted list x 0 ≤ … ≤ x (n − 1) that share their value with another entry. -/
noncomputable def shared (x : ℕ → ℝ) (n : ℕ) : Finset ℕ :=
  (range n).filter (fun i => ∃ j ∈ range n, j ≠ i ∧ x j = x i)

/-- Zero gaps: indices i with i + 1 < n and x (i + 1) = x i. -/
noncomputable def zeroGaps (x : ℕ → ℝ) (n : ℕ) : Finset ℕ :=
  (range n).filter (fun i => i + 1 < n ∧ x (i + 1) = x i)

lemma sorted_le {x : ℕ → ℝ} {n : ℕ} (hx : ∀ i, i + 1 < n → x i ≤ x (i + 1)) {i j : ℕ} (hij : i ≤ j)
    (hj : j < n) : x i ≤ x j := by
  induction j, hij using Nat.le_induction with
  | base => exact le_rfl
  | succ j hij ih => exact (ih (by omega)).trans (hx j hj)

/-- **Step 2, list form**: #(entries sharing their value with another entry) ≤ 2·#(zero gaps). -/
theorem card_shared_le (x : ℕ → ℝ) (n : ℕ) (hx : ∀ i, i + 1 < n → x i ≤ x (i + 1)) :
    (shared x n).card ≤ 2 * (zeroGaps x n).card := by
  have hsub : shared x n ⊆ zeroGaps x n ∪ (zeroGaps x n).image (· + 1) := by
    intro i hi
    simp only [shared, mem_filter, mem_range] at hi
    obtain ⟨hin, j, hjn, hji, hxj⟩ := hi
    rcases lt_or_gt_of_ne hji with h | h
    · apply mem_union_right
      rw [mem_image]
      refine ⟨i - 1, ?_, by show i - 1 + 1 = i; omega⟩
      have h1 : x j ≤ x (i - 1) := sorted_le hx (by omega) (by omega)
      have h2 : x (i - 1) ≤ x i := by
        have h3 := hx (i - 1) (by omega)
        rwa [show i - 1 + 1 = i by omega] at h3
      simp only [zeroGaps, mem_filter, mem_range]
      rw [show i - 1 + 1 = i by omega]
      exact ⟨by omega, by omega, by linarith⟩
    · apply mem_union_left
      simp only [zeroGaps, mem_filter, mem_range]
      have h1 : x i ≤ x (i + 1) := hx i (by omega)
      have h2 : x (i + 1) ≤ x j := sorted_le hx (by omega) hjn
      exact ⟨hin, by omega, by linarith⟩
  calc (shared x n).card ≤ (zeroGaps x n ∪ (zeroGaps x n).image (· + 1)).card := card_le_card hsub
    _ ≤ (zeroGaps x n).card + ((zeroGaps x n).image (· + 1)).card := card_union_le _ _
    _ ≤ (zeroGaps x n).card + (zeroGaps x n).card := Nat.add_le_add_left card_image_le _
    _ = 2 * (zeroGaps x n).card := by ring

/-! ## (a), zero model -/

/-- Real parts β n of the zeros, in pplean's ordinate order (Zeros.γ, with multiplicity). The functional
equation ρ ↦ 1 − ρ̄ permutes the zeros with multiplicity: an off-line zero (β ≠ 1/2) has a partner index with the
same ordinate and real part 1 − β. -/
def Reflect (Z : Zeros) (β : ℕ → ℝ) : Prop :=
  ∀ n : ℕ, 1 ≤ n → β n ≠ 1 / 2 → ∃ m : ℕ, 1 ≤ m ∧ Z.γ m = Z.γ n ∧ β m = 1 - β n

/-- The n-th zero is off-line or multiple. -/
def Exc (Z : Zeros) (β : ℕ → ℝ) (n : ℕ) : Prop :=
  β n ≠ 1 / 2 ∨ ∃ m : ℕ, 1 ≤ m ∧ m ≠ n ∧ Z.γ m = Z.γ n ∧ β m = β n

open Classical in
/-- Off-line or multiple zeros with T ≤ γ ≤ 2T. -/
noncomputable def excSet (Z : Zeros) (β : ℕ → ℝ) (T : ℝ) : Finset ℕ :=
  (Icc 1 (Z.N (2 * T))).filter (fun n => T ≤ Z.γ n ∧ Z.γ n ≤ 2 * T ∧ Exc Z β n)

/-- Z_off ∩ {T ≤ γ ≤ 2T} (zeros with β > 1/2): PROOF Step 6a's classes S_a are subsets. -/
noncomputable def offSet (Z : Zeros) (β : ℕ → ℝ) (T : ℝ) : Finset ℕ :=
  (Icc 1 (Z.N (2 * T))).filter (fun n => T ≤ Z.γ n ∧ Z.γ n ≤ 2 * T ∧ 1 / 2 < β n)

/-- 𝒩(T) of Theorem U: n with T ≤ γ_n ≤ γ_{n+1} ≤ 2T and γ_{n+1} − γ_n < 2πφ/log T, unfolded with the zero
model's sc T = log T/(2π): (γ_{n+1} − γ_n)·sc T < φ. (γ_n ≤ γ_{n+1} is Zeros.Ordered.mono; the ambient
Icc 1 N(2T) cuts nothing, since γ_{n+1} ≤ 2T forces n + 1 ≤ N(2T).) -/
noncomputable def gapSet (Z : Zeros) (φ T : ℝ) : Finset ℕ :=
  (Icc 1 (Z.N (2 * T))).filter
    (fun n => T ≤ Z.γ n ∧ Z.γ (n + 1) ≤ 2 * T ∧ (Z.γ (n + 1) - Z.γ n) * Z.sc T < φ)

/-- **Step 2 (zero model)**: #(off-line or multiple zeros with T ≤ γ ≤ 2T) ≤ 2𝒩(T). -/
theorem card_exc_le {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β) {φ : ℝ} (hφ : 0 < φ)
    (T : ℝ) : (excSet Z β T).card ≤ 2 * (gapSet Z φ T).card := by
  have hsub : excSet Z β T ⊆ gapSet Z φ T ∪ (gapSet Z φ T).image (· + 1) := by
    intro n hn
    simp only [excSet, mem_filter, mem_Icc] at hn
    obtain ⟨⟨hn1, hnN⟩, hT1, hT2, hexc⟩ := hn
    obtain ⟨m, hm1, hmn, hγ⟩ : ∃ m : ℕ, 1 ≤ m ∧ m ≠ n ∧ Z.γ m = Z.γ n := by
      rcases hexc with h | ⟨m, hm1, hmn, hγ, -⟩
      · obtain ⟨m, hm1, hγ, hβm⟩ := hβ n hn1 h
        refine ⟨m, hm1, ?_, hγ⟩
        intro hmn
        rw [hmn] at hβm
        apply h
        linarith
      · exact ⟨m, hm1, hmn, hγ⟩
    rcases lt_or_gt_of_ne hmn with h | h
    · apply mem_union_right
      rw [mem_image]
      refine ⟨n - 1, ?_, by show n - 1 + 1 = n; omega⟩
      have h1 : Z.γ m ≤ Z.γ (n - 1) := hZ.mono' hm1 (by omega)
      have h2 : Z.γ (n - 1) ≤ Z.γ n := by
        have h3 := hZ.mono (n - 1) (by omega)
        rwa [show n - 1 + 1 = n by omega] at h3
      simp only [gapSet, mem_filter, mem_Icc]
      rw [show n - 1 + 1 = n by omega]
      refine ⟨⟨by omega, by omega⟩, by linarith, hT2, ?_⟩
      have h4 : Z.γ n - Z.γ (n - 1) = 0 := by linarith
      rw [h4, zero_mul]
      exact hφ
    · apply mem_union_left
      have h1 : Z.γ n ≤ Z.γ (n + 1) := hZ.mono n hn1
      have h2 : Z.γ (n + 1) ≤ Z.γ m := hZ.mono' (by omega) (by omega)
      simp only [gapSet, mem_filter, mem_Icc]
      refine ⟨⟨hn1, hnN⟩, hT1, by linarith, ?_⟩
      have h4 : Z.γ (n + 1) - Z.γ n = 0 := by linarith
      rw [h4, zero_mul]
      exact hφ
  calc (excSet Z β T).card ≤ (gapSet Z φ T ∪ (gapSet Z φ T).image (· + 1)).card := card_le_card hsub
    _ ≤ (gapSet Z φ T).card + ((gapSet Z φ T).image (· + 1)).card := card_union_le _ _
    _ ≤ (gapSet Z φ T).card + (gapSet Z φ T).card := Nat.add_le_add_left card_image_le _
    _ = 2 * (gapSet Z φ T).card := by ring

lemma offSet_sub_excSet (Z : Zeros) (β : ℕ → ℝ) (T : ℝ) : offSet Z β T ⊆ excSet Z β T := by
  intro n hn
  simp only [offSet, excSet, mem_filter] at hn ⊢
  obtain ⟨h1, h2, h3, h4⟩ := hn
  exact ⟨h1, h2, h3, Or.inl h4.ne'⟩

/-- Step 2 as used in 6a: #(Z_off ∩ {T ≤ γ ≤ 2T}) ≤ 2𝒩(T). -/
theorem card_off_le {Z : Zeros} (hZ : Z.Ordered) {β : ℕ → ℝ} (hβ : Reflect Z β) {φ : ℝ} (hφ : 0 < φ)
    (T : ℝ) : (offSet Z β T).card ≤ 2 * (gapSet Z φ T).card :=
  (card_le_card (offSet_sub_excSet Z β T)).trans (card_exc_le hZ hβ hφ T)

lemma card_gapSet_le (Z : Zeros) (φ T : ℝ) : (gapSet Z φ T).card ≤ Z.N (2 * T) := by
  have h := card_le_card (filter_subset
    (fun n => T ≤ Z.γ n ∧ Z.γ (n + 1) ≤ 2 * T ∧ (Z.γ (n + 1) - Z.γ n) * Z.sc T < φ) (Icc 1 (Z.N (2 * T))))
  rw [Nat.card_Icc] at h
  unfold gapSet
  omega

end DPCLean

