import Riemann.Deep.SunComplete.Lean.Sieve

/-!
# sunlean · CaseA1: the case-(A) facts about `H_q` (NOTE §1, §3)

`τ > 0`, `E ≥ 0`; `e(I, I) ≤ ε|I|` on independent sets; `e(A, B) = 0` across an empty cut;
every degree of `H_q` is `≥ φ(q − 1) − ε`.
-/

namespace SunLean

theorem theta_pos (S : Finset ℕ) (hS : ∀ p ∈ S, p.Prime) : 0 < ∏ p ∈ S, (1 - 1 / (p : ℚ)) := by
  induction S using Finset.induction_on with
  | empty => simp
  | insert a S ha ih =>
    rw [Finset.prod_insert ha]
    have h2 : (2 : ℚ) ≤ a := by exact_mod_cast (hS a (Finset.mem_insert_self a S)).two_le
    have h3 : 1 / (a : ℚ) < 1 := (div_lt_one (by linarith)).mpr (by linarith)
    exact mul_pos (by linarith) (ih (fun p hp => hS p (Finset.mem_insert_of_mem hp)))

namespace SieveData
variable {q : ℕ} (d : SieveData q)

theorem core_prime : ∀ p ∈ d.core, p.Prime := by
  intro p hp
  have h : p ∈ primeDivisors (q - 1) := d.hsplit ▸ Finset.mem_union_left _ hp
  unfold primeDivisors at h
  exact (Finset.mem_filter.mp h).2.1

theorem τ_pos (hδ : 0 < d.δ) : 0 < d.τ := mul_pos (theta_pos _ d.core_prime) hδ

theorem WΔ_ge (hδ : 0 < d.δ) : 1 ≤ d.W * d.Δ := by
  have hW : (1 : ℚ) ≤ d.W := one_le_pow₀ (by norm_num)
  have hΔ : (1 : ℚ) ≤ d.Δ := by
    unfold SieveData.Δ
    split_ifs with h
    · exact le_refl _
    · have h1 : (1 : ℚ) ≤ d.siev.card := by exact_mod_cast Nat.one_le_iff_ne_zero.mpr h
      have : 0 ≤ ((d.siev.card : ℚ) - 1) / d.δ := div_nonneg (by linarith) hδ.le
      linarith
  exact one_le_mul_of_one_le_of_one_le hW hΔ

theorem E_nonneg (hδ : 0 < d.δ) {s : ℚ} (hs : 0 ≤ s) : 0 ≤ d.E s := by
  unfold SieveData.E
  have := d.WΔ_ge hδ
  have : 0 ≤ (d.W * d.Δ - 1) * s := mul_nonneg (by linarith) hs
  linarith
end SieveData

theorem epsA_nonneg (q : ℕ) : 0 ≤ epsA q := by
  unfold epsA; split_ifs <;> norm_num

variable {F : Type} [Field F] [Fintype F] [DecidableEq F]

/-- `e(I, I) ≤ ε|I|` for an independent set `I` of `H_q` (NOTE Prop 4 (i)). -/
theorem pairCount_diag (h2 : 2 ≤ Fintype.card F) (hF : FFacts F) {I : Finset F}
    (hI : IndepSet (sumGraph F) I) : (pairCount I I : ℚ) ≤ epsA (Fintype.card F) * I.card := by
  have hsub : (I ×ˢ I).filter (fun p => p.1 + p.2 ∈ primSet F) ⊆
      (I.filter (fun x => x + x ∈ primSet F)).image (fun x => (x, x)) := by
    intro p hp
    rw [Finset.mem_filter, Finset.mem_product] at hp
    obtain ⟨⟨h1, h2⟩, h3⟩ := hp
    have hpe : p.1 = p.2 := by
      by_contra hne
      exact hI p.1 h1 p.2 h2 hne ⟨hne, mem_primSet.mp h3⟩
    rw [Finset.mem_image]
    refine ⟨p.1, Finset.mem_filter.mpr ⟨h1, ?_⟩, Prod.ext rfl hpe⟩
    have h4 := h3
    rw [← hpe] at h4
    exact h4
  have h3 : pairCount I I ≤ (I.filter (fun x => x + x ∈ primSet F)).card :=
    (Finset.card_le_card hsub).trans Finset.card_image_le
  unfold epsA
  split_ifs with hq
  · have hD : I.filter (fun x => x + x ∈ primSet F) = ∅ := by
      rw [Finset.filter_eq_empty_iff]
      intro x _ hx
      have h0 : x + x = 0 := by rw [← two_mul, hF.two_eq_zero hq, zero_mul]
      rw [h0] at hx
      exact not_isPrim_zero h2 (mem_primSet.mp hx)
    rw [hD, Finset.card_empty] at h3
    have h5 : pairCount I I = 0 := by omega
    rw [h5]
    simp
  · have h4 := Finset.card_filter_le I (fun x => x + x ∈ primSet F)
    have h5 : pairCount I I ≤ I.card := le_trans h3 h4
    rw [one_mul]
    exact_mod_cast h5

/-- `e(A, B) = 0` when `A, B` are disjoint with no edge of `H_q` between them. -/
theorem pairCount_cross {A B : Finset F} (hAB : Disjoint A B)
    (hno : ∀ x ∈ A, ∀ y ∈ B, ¬ (sumGraph F).Adj x y) : pairCount A B = 0 := by
  unfold pairCount
  rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
  intro p hp h3
  rw [Finset.mem_product] at hp
  have hne : p.1 ≠ p.2 := by
    intro h
    have hb : p.1 ∈ B := by rw [h]; exact hp.2
    exact Finset.disjoint_left.mp hAB hp.1 hb
  exact hno p.1 hp.1 p.2 hp.2 ⟨hne, mem_primSet.mp h3⟩

/-- Every degree of `H_q` is at least `δ₀ = φ(q − 1) − ε` (NOTE §1). -/
theorem deg_ge (hF : FFacts F) (h2 : 2 ≤ Fintype.card F) (v : F) :
    delta0A (Fintype.card F) ≤ ((Finset.univ.filter ((sumGraph F).Adj v)).card : ℚ) := by
  have hcount : (Finset.univ.filter (fun y => v + y ∈ primSet F)).card = (primSet F).card := by
    have e : Finset.univ.filter (fun y => v + y ∈ primSet F) = (primSet F).image (fun z => z - v) := by
      ext y
      simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_image]
      constructor
      · intro h
        exact ⟨v + y, h, by ring⟩
      · rintro ⟨z, hz, rfl⟩
        have : v + (z - v) = z := by ring
        rw [this]
        exact hz
    rw [e, Finset.card_image_of_injective _ sub_left_injective]
  unfold delta0A epsA
  rw [← hF.card_prim, ← hcount]
  split_ifs with hq
  · have hsub : Finset.univ.filter (fun y => v + y ∈ primSet F) ⊆
        Finset.univ.filter ((sumGraph F).Adj v) := by
      intro y hy
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hy ⊢
      refine ⟨?_, mem_primSet.mp hy⟩
      rintro rfl
      have h0 : v + v = 0 := by rw [← two_mul, hF.two_eq_zero hq, zero_mul]
      rw [h0] at hy
      exact not_isPrim_zero h2 (mem_primSet.mp hy)
    have := Finset.card_le_card hsub
    rw [sub_zero]
    exact_mod_cast this
  · have hsub : Finset.univ.filter (fun y => v + y ∈ primSet F) ⊆
        insert v (Finset.univ.filter ((sumGraph F).Adj v)) := by
      intro y hy
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hy
      by_cases hyv : y = v
      · rw [hyv]; exact Finset.mem_insert_self _ _
      · apply Finset.mem_insert_of_mem
        simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        exact ⟨Ne.symm hyv, mem_primSet.mp hy⟩
    have h1 := Finset.card_le_card hsub
    have h3 := Finset.card_insert_le v (Finset.univ.filter ((sumGraph F).Adj v))
    have h4 : ((Finset.univ.filter (fun y => v + y ∈ primSet F)).card : ℚ) ≤
        ((Finset.univ.filter ((sumGraph F).Adj v)).card : ℚ) + 1 := by
      exact_mod_cast le_trans h1 h3
    linarith

end SunLean
