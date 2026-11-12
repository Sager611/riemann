import Riemann.Deep.EvoAlg.Lean.CorA

/-!
# More examples over `ℚ` (NOTE §0, Corollary A): `K_{1,3}`, `K_{2,3}` are not isomorphic to their
random-walk algebras over `ℚ` (ratios `3`, `3/2` are not rational cubes); `K_{1,8}` is (ratio `8 = 2³`).
-/

namespace EvoAlg

theorem nat_cube_ne' {m k j : ℕ} (h1 : j ^ 3 < k) (h2 : k < (j + 1) ^ 3) : m ^ 3 ≠ k := by
  intro h
  rcases Nat.lt_or_ge j m with hm | hm
  · have := Nat.pow_le_pow_left (by omega : j + 1 ≤ m) 3
    omega
  · have := Nat.pow_le_pow_left hm 3
    omega

/-- A natural number strictly between consecutive cubes is not a rational cube. -/
theorem not_cube_nat (k j : ℕ) (h1 : j ^ 3 < k) (h2 : k < (j + 1) ^ 3) :
    ¬ ∃ t : ℚ, t ^ 3 = k := by
  rintro ⟨t, ht⟩
  have hnum : t.num ^ 3 = k := by
    have := congrArg Rat.num ht
    rw [Rat.num_pow] at this
    simpa using this
  have h3 : t.num.natAbs ^ 3 = k := by
    have := congrArg Int.natAbs hnum
    simpa [Int.natAbs_pow] using this
  exact nat_cube_ne' h1 h2 h3

theorem cbg_connected (m n : ℕ) [NeZero m] [NeZero n] :
    (completeBipartiteGraph (Fin m) (Fin n)).Connected := by
  have hadj : ∀ a b, (completeBipartiteGraph (Fin m) (Fin n)).Adj (Sum.inl a) (Sum.inr b) :=
    fun _ _ => Or.inl ⟨rfl, rfl⟩
  have h : ∀ v, (completeBipartiteGraph (Fin m) (Fin n)).Reachable (Sum.inl 0) v := by
    intro v
    rcases v with a | b
    · exact (hadj 0 0).reachable.trans (hadj a 0).symm.reachable
    · exact (hadj 0 b).reachable
  exact (SimpleGraph.connected_iff_exists_forall_reachable _).mpr ⟨Sum.inl 0, h⟩

/-- For `K_{m,n}` (degrees `n` on the left, `m` on the right, `m ≠ n`): if neither `n/m` nor
`m/n` is a rational cube, then `A(K_{m,n}) ≇ A_RW(K_{m,n})` over `ℚ`. -/
theorem cbg_not_iso (m n : ℕ) [NeZero m] [NeZero n] (hmn : m ≠ n)
    (hl : ∀ a, (completeBipartiteGraph (Fin m) (Fin n)).degree (Sum.inl a) = n)
    (hr : ∀ b, (completeBipartiteGraph (Fin m) (Fin n)).degree (Sum.inr b) = m)
    (hcard : 2 ≤ Fintype.card (Fin m ⊕ Fin n))
    (h1 : ¬ ∃ t : ℚ, t ^ 3 = (n : ℚ) / m) (h2 : ¬ ∃ t : ℚ, t ^ 3 = (m : ℚ) / n) :
    ¬ AIsoRW (completeBipartiteGraph (Fin m) (Fin n)) ℚ := by
  intro h
  rcases (theoremA _ (cbg_connected m n) hcard).mp h with
    ⟨d, hd⟩ | ⟨s, t, d₁, d₂, hb, c, hc, hcube⟩
  · have e1 := hd (Sum.inl 0)
    have e2 := hd (Sum.inr 0)
    rw [hl] at e1
    rw [hr] at e2
    exact hmn (e2.trans e1.symm)
  · have hadj : (completeBipartiteGraph (Fin m) (Fin n)).Adj (Sum.inl 0) (Sum.inr 0) :=
      Or.inl ⟨rfl, rfl⟩
    rcases hb.bip.mem_of_adj hadj with ⟨e1, e2⟩ | ⟨e1, e2⟩
    · have f1 := hb.deg₁ _ e1
      have f2 := hb.deg₂ _ e2
      rw [hl] at f1
      rw [hr] at f2
      subst f1 f2
      exact h1 ⟨c, hcube.symm⟩
    · have f1 := hb.deg₂ _ e1
      have f2 := hb.deg₁ _ e2
      rw [hl] at f1
      rw [hr] at f2
      subst f1 f2
      exact h2 ⟨c, hcube.symm⟩

/-- `A(K_{1,3}) ≇ A_RW(K_{1,3})` over `ℚ` (`3` and `1/3` are not rational cubes). -/
theorem corA_K13 : ¬ AIsoRW (completeBipartiteGraph (Fin 1) (Fin 3)) ℚ :=
  cbg_not_iso 1 3 (by decide) (by decide) (by decide) (by decide)
    (fun ⟨t, ht⟩ => not_cube_nat 3 1 (by norm_num) (by norm_num) ⟨t, by rw [ht]; norm_num⟩)
    (fun ⟨t, ht⟩ => not_cube_nat 9 2 (by norm_num) (by norm_num)
      ⟨3 * t, by rw [mul_pow, ht]; norm_num⟩)

/-- `A(K_{2,3}) ≇ A_RW(K_{2,3})` over `ℚ` (`3/2` and `2/3` are not rational cubes). -/
theorem corA_K23 : ¬ AIsoRW (completeBipartiteGraph (Fin 2) (Fin 3)) ℚ :=
  cbg_not_iso 2 3 (by decide) (by decide) (by decide) (by decide)
    (fun ⟨t, ht⟩ => not_cube_nat 12 2 (by norm_num) (by norm_num)
      ⟨2 * t, by rw [mul_pow, ht]; norm_num⟩)
    (fun ⟨t, ht⟩ => not_cube_nat 18 2 (by norm_num) (by norm_num)
      ⟨3 * t, by rw [mul_pow, ht]; norm_num⟩)

theorem K18_biregular : IsBiregular (completeBipartiteGraph (Fin 1) (Fin 8)) (Set.range Sum.inl)
    (Set.range Sum.inr) 8 1 where
  bip := SimpleGraph.IsBipartiteWith.completeBipartiteGraph (Fin 1) (Fin 8)
  cover v := by
    rcases v with a | b
    · exact Or.inl ⟨a, rfl⟩
    · exact Or.inr ⟨b, rfl⟩
  deg₁ v hv := by
    obtain ⟨a, rfl⟩ := hv
    revert a
    decide
  deg₂ v hv := by
    obtain ⟨b, rfl⟩ := hv
    revert b
    decide

/-- `A(K_{1,8}) ≅ A_RW(K_{1,8})` over `ℚ`: biregular with `d₁/d₂ = 8/1 = 2³`. -/
theorem K18_iso : AIsoRW (completeBipartiteGraph (Fin 1) (Fin 8)) ℚ :=
  aIsoRW_of_biregular _ K18_biregular (c := (2 : ℚ)) two_ne_zero (by norm_num)

end EvoAlg

