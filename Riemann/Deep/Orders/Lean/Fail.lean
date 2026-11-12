import Riemann.Deep.Orders.Lean.Search

/-! The 13 failing orders (NOTE.md §2, Proposition 1): no k-th move can make a new number. -/
namespace Orders

/-- + − × ÷: no legal 2-move prefix, so a perfect N has N ≤ 2. -/
theorem fail_pmxd : ∀ N ms, Perfect (⟨.add, .sub, .mul, .div⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 2) (by decide +kernel) N ms h) (by decide)

/-- + − ÷ ×: no legal 2-move prefix, so a perfect N has N ≤ 2. -/
theorem fail_pmdx : ∀ N ms, Perfect (⟨.add, .sub, .div, .mul⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 2) (by decide +kernel) N ms h) (by decide)

/-- × − + ÷: no legal 2-move prefix, so a perfect N has N ≤ 2. -/
theorem fail_xmpd : ∀ N ms, Perfect (⟨.mul, .sub, .add, .div⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 2) (by decide +kernel) N ms h) (by decide)

/-- × − ÷ +: no legal 2-move prefix, so a perfect N has N ≤ 2. -/
theorem fail_xmdp : ∀ N ms, Perfect (⟨.mul, .sub, .div, .add⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 2) (by decide +kernel) N ms h) (by decide)

/-- − + × ÷: no legal 1-move prefix, so a perfect N has N ≤ 1. -/
theorem fail_mpxd : ∀ N ms, Perfect (⟨.sub, .add, .mul, .div⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 1) (by decide +kernel) N ms h) (by decide)

/-- − + ÷ ×: no legal 1-move prefix, so a perfect N has N ≤ 1. -/
theorem fail_mpdx : ∀ N ms, Perfect (⟨.sub, .add, .div, .mul⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 1) (by decide +kernel) N ms h) (by decide)

/-- − × + ÷: no legal 1-move prefix, so a perfect N has N ≤ 1. -/
theorem fail_mxpd : ∀ N ms, Perfect (⟨.sub, .mul, .add, .div⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 1) (by decide +kernel) N ms h) (by decide)

/-- − × ÷ +: no legal 1-move prefix, so a perfect N has N ≤ 1. -/
theorem fail_mxdp : ∀ N ms, Perfect (⟨.sub, .mul, .div, .add⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 1) (by decide +kernel) N ms h) (by decide)

/-- − ÷ + ×: no legal 1-move prefix, so a perfect N has N ≤ 1. -/
theorem fail_mdpx : ∀ N ms, Perfect (⟨.sub, .div, .add, .mul⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 1) (by decide +kernel) N ms h) (by decide)

/-- − ÷ × +: no legal 1-move prefix, so a perfect N has N ≤ 1. -/
theorem fail_mdxp : ∀ N ms, Perfect (⟨.sub, .div, .mul, .add⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 1) (by decide +kernel) N ms h) (by decide)

/-- ÷ × − +: no legal 5-move prefix, so a perfect N has N ≤ 5. -/
theorem fail_dxmp : ∀ N ms, Perfect (⟨.div, .mul, .sub, .add⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 5) (by decide +kernel) N ms h) (by decide)

/-- ÷ − + ×: no legal 2-move prefix, so a perfect N has N ≤ 2. -/
theorem fail_dmpx : ∀ N ms, Perfect (⟨.div, .sub, .add, .mul⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 2) (by decide +kernel) N ms h) (by decide)

/-- ÷ − × +: no legal 2-move prefix, so a perfect N has N ≤ 2. -/
theorem fail_dmxp : ∀ N ms, Perfect (⟨.div, .sub, .mul, .add⟩ : Order) N ms → N ≤ 5 := fun N ms h =>
  Nat.le_trans (bound_of_reach (k := 2) (by decide +kernel) N ms h) (by decide)

end Orders
