import Riemann.Deep.LeanPrefix.Prop5A

/-!
# Proposition 5, arithmetic core, Steps 2-4 and assembly (prop5_arith)

Step 2: a run R_i ⊆ I (i ≤ d - 2) is an integer progression with step e ≡ αd (mod N),
(L-1)|e| ≤ α - 2, so L e + α ∈ {0, N}. Step 3: L e = β contradicts Star on R_0. Step 4: L e = -α,
λ = α/L: for 1 ≤ t ≤ λ - 1 the run i = t - 1 + ⌊tβ/α⌋ has last point tβ mod α and first point
tβ mod α + α - λ, so tβ mod α ≤ λ - 1; then prop5_finish gives β ≡ 1 (mod α).
-/

namespace Riemann.Deep.LeanPrefix

section prop5

variable {n α β d L : ℕ}

/-- Step 2: a run inside I is an integer progression; its step e satisfies e ≡ αd (mod N) and
|(L-1) e| ≤ α - 2. -/
theorem prop5_run_prog (hα : 3 ≤ α) (hαn : 2 * α < n + 1) (hc : Nat.Coprime (α * d) (n + 1))
    (hdL : d * L = n) (hL : 3 ≤ L) {i : ℕ} (hi : i + 2 ≤ d)
    (hrun : ∀ s, i * L < s → s ≤ (i + 1) * L → Yf n α d s < α) :
    ∃ e : ℤ, ((e : ℤ) : ZMod (n + 1)) = (α : ZMod (n + 1)) * d ∧
      ((L : ℤ) - 1) * e ≤ (α : ℤ) - 2 ∧ -((α : ℤ) - 2) ≤ ((L : ℤ) - 1) * e := by
  have hiL : (i + 1) * L = i * L + L := Nat.succ_mul i L
  have hdn : (i + 1) * L ≤ n := by
    rw [← hdL]; exact Nat.mul_le_mul_right L (by omega)
  have hpos : ∀ s, i * L < s → s ≤ (i + 1) * L → 1 ≤ Yf n α d s := fun s h1 h2 =>
    Yf_pos hc (by omega) (by omega) (by omega)
  obtain ⟨a, ha⟩ : ∃ a : ℤ, a = Yf n α d (i * L + 1) := ⟨_, rfl⟩
  obtain ⟨e, he⟩ : ∃ e : ℤ, e = (Yf n α d (i * L + 2) : ℤ) - Yf n α d (i * L + 1) := ⟨_, rfl⟩
  have hecast : ((e : ℤ) : ZMod (n + 1)) = (α : ZMod (n + 1)) * d := by
    rw [he]; push_cast; rw [Yf_cast, Yf_cast]; push_cast; ring
  have hb1 := hrun (i * L + 1) (by omega) (by omega)
  have hb2 := hrun (i * L + 2) (by omega) (by omega)
  have hp1 := hpos (i * L + 1) (by omega) (by omega)
  have hp2 := hpos (i * L + 2) (by omega) (by omega)
  have hprog : ∀ j, j ≤ L - 1 → (Yf n α d (i * L + 1 + j) : ℤ) = a + j * e := by
    intro j hj
    induction j with
    | zero => rw [Nat.add_zero, ha]; push_cast; ring
    | succ j ih =>
      have ih' := ih (by omega)
      have hu1 := hrun (i * L + 1 + j) (by omega) (by omega)
      have hu2 := hrun (i * L + 1 + (j + 1)) (by omega) (by omega)
      have hq1 := hpos (i * L + 1 + j) (by omega) (by omega)
      have hq2 := hpos (i * L + 1 + (j + 1)) (by omega) (by omega)
      have hX : ((Yf n α d (i * L + 1 + (j + 1)) : ℤ) - Yf n α d (i * L + 1 + j) - e) = 0 := by
        apply int_eq_zero_of_dvd_of_bounds (N := n + 1)
        · rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
          push_cast
          rw [hecast, Yf_cast, Yf_cast]
          push_cast; ring
        · rw [he]; push_cast; omega
        · rw [he]; push_cast; omega
      push_cast
      linear_combination hX + ih'
  have hlast := hprog (L - 1) le_rfl
  rw [show i * L + 1 + (L - 1) = (i + 1) * L by omega] at hlast
  have hl1 := hrun ((i + 1) * L) (by omega) le_rfl
  have hl2 := hpos ((i + 1) * L) (by omega) le_rfl
  have hcastL : ((L - 1 : ℕ) : ℤ) = (L : ℤ) - 1 := by push_cast [Nat.cast_sub (show 1 ≤ L by omega)]; ring
  rw [hcastL] at hlast
  refine ⟨e, hecast, ?_, ?_⟩ <;> rw [ha] at hlast <;> push_cast at hlast ⊢ <;> omega

/-- Step 4: with Δ ≡ -λ (L λ = α), tβ mod α ≤ λ - 1 for 1 ≤ t ≤ λ - 1. -/
theorem prop5_step4 (hN : α + β = n + 1) (hαβ : α < β) (hdL : d * L = n) (hL : 3 ≤ L)
    (hstar : ∀ i s s', i + 2 ≤ d → i * L < s → s ≤ (i + 1) * L → i * L < s' → s' ≤ (i + 1) * L →
      (Yf n α d s < α ↔ Yf n α d s' < α))
    {lam : ℕ} (hlam : L * lam = α) (hΔ : (α : ZMod (n + 1)) * d = -(lam : ZMod (n + 1)))
    {t : ℕ} (ht1 : 1 ≤ t) (ht : t + 1 ≤ lam) : t * β % α + 1 ≤ lam := by
  have hαpos : 0 < α := by
    rcases Nat.eq_zero_or_pos α with h | h
    · rw [h] at hlam; rcases Nat.mul_eq_zero.mp hlam with h1 | h1 <;> omega
    · exact h
  have hlamα : lam ≤ α := by rw [← hlam]; exact Nat.le_mul_of_pos_left lam (by omega)
  obtain ⟨q, hq⟩ : ∃ q, q = t * β / α := ⟨_, rfl⟩
  obtain ⟨ρ, hρ⟩ : ∃ ρ, ρ = t * β % α := ⟨_, rfl⟩
  have hdiv : α * q + ρ = t * β := by rw [hq, hρ]; exact Nat.div_add_mod _ _
  have hρα : ρ < α := by rw [hρ]; exact Nat.mod_lt _ hαpos
  obtain ⟨i, hi⟩ : ∃ i, i = t - 1 + q := ⟨_, rfl⟩
  have hi1 : i + 1 = t + q := by omega
  have hid : (i + 1) * α + ρ = t * (n + 1) := by
    rw [hi1, Nat.add_mul, ← hN, Nat.mul_add, Nat.mul_comm q α, Nat.add_assoc, hdiv]
  -- i + 2 ≤ d
  have hqα : q * α ≤ t * β := by rw [hq]; exact Nat.div_mul_le_self _ _
  have h1 : (t + 1) * α ≤ lam * α := Nat.mul_le_mul_right α ht
  have h2 : t * β ≤ (lam - 1) * β := Nat.mul_le_mul_right β (by omega)
  have h3 : (lam - 1) * β + β = lam * β := by
    rw [← Nat.succ_mul, Nat.succ_eq_add_one, Nat.sub_add_cancel (by omega)]
  have h4 : (n + 1) * lam = lam * α + lam * β := by rw [← hN, Nat.add_mul, Nat.mul_comm α, Nat.mul_comm β]
  have h5 : (i + 2) * α = (t + 1) * α + q * α := by
    rw [show i + 2 = (t + 1) + q by omega, Nat.add_mul]
  have h6 : d * α + lam = (n + 1) * lam := by
    rw [← hlam, ← Nat.mul_assoc, hdL, Nat.succ_mul]
  have hlamβ : lam ≤ β := by omega
  have hi2 : i + 2 ≤ d := by
    have : (i + 2) * α ≤ d * α := by omega
    exact Nat.le_of_mul_le_mul_right this hαpos
  -- the two points of run i
  have hiL : (i + 1) * L = i * L + L := Nat.succ_mul i L
  have hlamZ : (L : ZMod (n + 1)) * lam = α := by
    have := congrArg (Nat.cast : ℕ → ZMod (n + 1)) hlam; push_cast at this; exact this
  have hρZ : ((i : ZMod (n + 1)) + 1) * α + ρ = t * ((n : ZMod (n + 1)) + 1) := by
    have := congrArg (Nat.cast : ℕ → ZMod (n + 1)) hid; push_cast at this; exact this
  rw [n_cast] at hρZ
  have hY2 : Yf n α d ((i + 1) * L) = ρ := by
    apply Yf_eq (by omega)
    push_cast
    linear_combination ((i : ZMod (n + 1)) + 1) * L * hΔ - ((i : ZMod (n + 1)) + 1) * hlamZ - hρZ
  have hY1 : Yf n α d (i * L + 1) = ρ + α - lam := by
    apply Yf_eq (by omega)
    rw [Nat.cast_sub (by omega)]
    push_cast
    linear_combination ((i : ZMod (n + 1)) * L + 1) * hΔ - (i : ZMod (n + 1)) * hlamZ - hρZ
  have hs := hstar i (i * L + 1) ((i + 1) * L) hi2 (by omega) (by omega) (by omega) le_rfl
  rw [hY1, hY2] at hs
  have := hs.mpr hρα
  rw [← hρ]
  omega

/-- Steps 2-4: if some run R_i (i ≤ d - 2) lies in I, then α ∣ n. -/
theorem prop5_case2 (hN : α + β = n + 1) (hα : 3 ≤ α) (hαβ : α < β) (hcop : Nat.Coprime α β)
    (hdL : d * L = n) (hd : 3 ≤ d) (hL : 3 ≤ L)
    (hstar : ∀ i s s', i + 2 ≤ d → i * L < s → s ≤ (i + 1) * L → i * L < s' → s' ≤ (i + 1) * L →
      (Yf n α d s < α ↔ Yf n α d s' < α))
    {i : ℕ} (hi : i + 2 ≤ d) (hrun : ∀ s, i * L < s → s ≤ (i + 1) * L → Yf n α d s < α) :
    α ∣ n := by
  have hc := coprime_ad hN hcop hdL
  have hαn : 2 * α < n + 1 := by omega
  obtain ⟨e, hecast, hb1, hb2⟩ := prop5_run_prog hα hαn hc hdL hL hi hrun
  have hdvd : ((n + 1 : ℕ) : ℤ) ∣ (L : ℤ) * e + α := by
    rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
    push_cast
    rw [hecast]
    linear_combination L_ad (α := α) hdL
  obtain ⟨c, hcdef⟩ := hdvd
  have hLz : (3 : ℤ) ≤ L := by exact_mod_cast hL
  have hNz : 2 * (α : ℤ) < (n : ℤ) + 1 := by exact_mod_cast hαn
  have hlow : -((n : ℤ) + 1) < (L : ℤ) * e + α := by
    rcases le_or_gt 0 e with he | he
    · nlinarith
    · nlinarith
  have hhigh : (L : ℤ) * e + α < 2 * ((n : ℤ) + 1) := by
    rcases le_or_gt 0 e with he | he
    · nlinarith
    · nlinarith
  push_cast at hcdef
  have hc01 : c = 0 ∨ c = 1 := by
    rcases lt_trichotomy c 0 with h | h | h
    · exfalso; nlinarith
    · exact Or.inl h
    · rcases lt_trichotomy c 1 with h' | h' | h'
      · omega
      · exact Or.inr h'
      · exfalso; nlinarith
  rcases hc01 with h0 | h1
  · -- L e = -α: λ = -e
    rw [h0, mul_zero] at hcdef
    have hneg : e < 0 := by nlinarith
    obtain ⟨lam, hlamdef⟩ : ∃ lam : ℕ, (lam : ℤ) = -e := ⟨(-e).toNat, by omega⟩
    have hlam : L * lam = α := by
      have : (L : ℤ) * lam = α := by rw [hlamdef]; linarith
      exact_mod_cast this
    have hΔ : (α : ZMod (n + 1)) * d = -(lam : ZMod (n + 1)) := by
      rw [← hecast, show e = -(lam : ℤ) by omega]; push_cast; ring
    have hlam1 : 1 ≤ lam := by omega
    rcases Nat.lt_or_ge lam 2 with hl | hl
    · have : lam = 1 := by omega
      rw [this, Nat.mul_one] at hlam
      exact ⟨d, by rw [← hdL, ← hlam, Nat.mul_comm]⟩
    · have h3 : 3 * lam ≤ α := by rw [← hlam]; exact Nat.mul_le_mul_right lam hL
      have hβ0 : β % α ≠ 0 := by
        intro h
        have hdv : α ∣ β := Nat.dvd_of_mod_eq_zero h
        have := Nat.Coprime.eq_one_of_dvd hcop hdv
        omega
      have hfin := prop5_finish hl h3 hβ0 (fun t ht1 ht => prop5_step4 hN hαβ hdL hL hstar hlam hΔ ht1 ht)
      have hdm := Nat.div_add_mod β α
      rw [hfin] at hdm
      exact ⟨1 + β / α, by rw [Nat.mul_add, Nat.mul_one]; omega⟩
  · -- L e = β: contradiction on run 0
    rw [h1, mul_one] at hcdef
    have hpos : 0 < e := by nlinarith
    obtain ⟨e', he'⟩ : ∃ e' : ℕ, (e' : ℤ) = e := ⟨e.toNat, by omega⟩
    have hLe : L * e' = β := by
      have : (L : ℤ) * e' = β := by rw [he']; linarith
      exact_mod_cast this
    have he'le : e' + 2 ≤ α := by
      have : ((L : ℤ) - 1) * e' ≤ (α : ℤ) - 2 := by rw [he']; exact hb1
      have h2 : (e' : ℤ) ≤ ((L : ℤ) - 1) * e' := by nlinarith
      omega
    have hY1 : Yf n α d 1 = e' := by
      apply Yf_eq (by omega)
      rw [← he'] at hecast
      push_cast at hecast ⊢
      rw [one_mul, ← hecast]
    have hYL : Yf n α d L = β := by
      apply Yf_eq (by omega)
      rw [← he'] at hecast
      push_cast at hecast
      rw [mul_assoc, ← hecast]
      have := congrArg (Nat.cast : ℕ → ZMod (n + 1)) hLe; push_cast at this; exact this
    have hs := hstar 0 1 L (by omega) (by omega) (by omega) (by omega) (by omega)
    rw [hY1, hYL] at hs
    have := hs.mp (by omega)
    omega

/-- Proposition 5, arithmetic core: Star ⇒ α ∣ n. -/
theorem prop5_arith (hN : α + β = n + 1) (hα : 3 ≤ α) (hαβ : α < β) (hcop : Nat.Coprime α β)
    (hdL : d * L = n) (hd : 3 ≤ d) (hL : 3 ≤ L)
    (hstar : ∀ i s s', i + 2 ≤ d → i * L < s → s ≤ (i + 1) * L → i * L < s' → s' ≤ (i + 1) * L →
      (Yf n α d s < α ↔ Yf n α d s' < α)) :
    α ∣ n := by
  have hc := coprime_ad hN hcop hdL
  by_cases hex : ∃ i s, i + 2 ≤ d ∧ i * L < s ∧ s ≤ (i + 1) * L ∧ Yf n α d s < α
  · obtain ⟨i, s, hi, h1, h2, h3⟩ := hex
    exact prop5_case2 hN hα hαβ hcop hdL hd hL hstar hi
      (fun s' h1' h2' => (hstar i s s' hi h1 h2 h1' h2').mp h3)
  · push Not at hex
    have hnone : ∀ s, 1 ≤ s → s + L ≤ n → ¬ Yf n α d s < α := by
      intro s hs1 hs2 hY
      obtain ⟨i, hi⟩ : ∃ i, i = (s - 1) / L := ⟨_, rfl⟩
      have hdm := Nat.div_add_mod (s - 1) L
      have hml := Nat.mod_lt (s - 1) (show 0 < L by omega)
      rw [← hi] at hdm
      have hiL : (i + 1) * L = i * L + L := Nat.succ_mul i L
      have hLi : L * i = i * L := Nat.mul_comm L i
      have hi2 : i + 2 ≤ d := by
        by_contra hcon
        have : (d - 1) * L ≤ i * L := Nat.mul_le_mul_right L (by omega)
        have e2 : (d - 1) * L + L = d * L := by
          rw [← Nat.succ_mul, Nat.succ_eq_add_one, Nat.sub_add_cancel (by omega)]
        omega
      exact absurd hY (not_lt.mpr (hex i s hi2 (by omega) (by omega)))
    obtain ⟨τ, hτ⟩ := prop5_case1 hα (by omega) hc hdL hd hL hnone
    exact ⟨τ * d, by rw [← hdL, hτ]; ring⟩

end prop5

end Riemann.Deep.LeanPrefix

