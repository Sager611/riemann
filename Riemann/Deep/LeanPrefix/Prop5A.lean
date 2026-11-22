import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Int.GCD
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith

/-!
# Proposition 5, arithmetic core (NOTE §4, Steps 1-4), in x-coordinates

N = n + 1 = α + β, gcd(α, β) = 1, n = d L with d, L ≥ 3. Y(s) = s α d mod N is the s-th point
of the orbit of the rotation by Δ = αd (NOTE Lemma 4: y_s = sΔ); s ↦ Y(s) is a bijection of
[1, n]. The runs are R_i = {Y(s) : iL < s ≤ (i+1)L}; R_{d-1} is the exempt class. I = [1, α-1].
Hypothesis Star: every run R_i with i ≤ d - 2 lies in I or misses I (Y(s) < α ↔ Y(s') < α).
Conclusion: α ∣ n.
-/

namespace Riemann.Deep.LeanPrefix

/-- Y(s) = s α d mod (n + 1). -/
def Yf (n α d s : ℕ) : ℕ := s * (α * d) % (n + 1)

theorem n_cast (n : ℕ) : (n : ZMod (n + 1)) = -1 := by
  have h : ((n + 1 : ℕ) : ZMod (n + 1)) = 0 := ZMod.natCast_self _
  push_cast at h
  linear_combination h

theorem int_eq_zero_of_dvd_of_bounds {N : ℕ} {x : ℤ} (h : (N : ℤ) ∣ x) (h1 : -(N : ℤ) < x)
    (h2 : x < N) : x = 0 := by
  obtain ⟨c, rfl⟩ := h
  rcases lt_trichotomy c 0 with hc | hc | hc
  · nlinarith
  · simp [hc]
  · nlinarith

section arith

variable {n α d : ℕ}

theorem Yf_lt (s : ℕ) : Yf n α d s < n + 1 := Nat.mod_lt _ (Nat.succ_pos n)

theorem Yf_cast (s : ℕ) : ((Yf n α d s : ℕ) : ZMod (n + 1)) = (s : ZMod (n + 1)) * α * d := by
  unfold Yf; rw [ZMod.natCast_mod]; push_cast; ring

theorem Yf_eq {s y : ℕ} (hy : y < n + 1) (h : (s : ZMod (n + 1)) * α * d = y) :
    Yf n α d s = y := by
  have h1 : ((Yf n α d s : ℕ) : ZMod (n + 1)) = (y : ZMod (n + 1)) := by rw [Yf_cast, h]
  rw [ZMod.natCast_eq_natCast_iff', Nat.mod_eq_of_lt (Yf_lt s), Nat.mod_eq_of_lt hy] at h1
  exact h1

theorem exists_inv (h : Nat.Coprime (α * d) (n + 1)) (hn : 1 ≤ n) :
    ∃ ε : ℕ, ε < n + 1 ∧ 0 < ε ∧ (ε : ZMod (n + 1)) * α * d = 1 := by
  obtain ⟨ε, hε, h1⟩ := Nat.exists_mul_mod_eq_one_of_coprime h (by omega)
  refine ⟨ε, hε, ?_, ?_⟩
  · rcases Nat.eq_zero_or_pos ε with h0 | h0
    · rw [h0, Nat.mul_zero, Nat.zero_mod] at h1; omega
    · exact h0
  have h2 : ((α * d * ε % (n + 1) : ℕ) : ZMod (n + 1)) = 1 := by rw [h1]; simp
  rw [ZMod.natCast_mod] at h2
  push_cast at h2
  linear_combination h2

theorem Yf_inj (hc : Nat.Coprime (α * d) (n + 1)) (hn : 1 ≤ n) {s s' : ℕ} (hs : s ≤ n)
    (hs' : s' ≤ n) (h : Yf n α d s = Yf n α d s') : s = s' := by
  obtain ⟨ε, -, -, hε⟩ := exists_inv hc hn
  have h1 : ((Yf n α d s : ℕ) : ZMod (n + 1)) = ((Yf n α d s' : ℕ) : ZMod (n + 1)) := by rw [h]
  rw [Yf_cast, Yf_cast] at h1
  have h2 : (s : ZMod (n + 1)) = s' := by
    linear_combination (ε : ZMod (n + 1)) * h1 - ((s : ZMod (n + 1)) - s') * hε
  rw [ZMod.natCast_eq_natCast_iff', Nat.mod_eq_of_lt (by omega), Nat.mod_eq_of_lt (by omega)] at h2
  exact h2

theorem Yf_zero : Yf n α d 0 = 0 := by simp [Yf]

theorem Yf_pos (hc : Nat.Coprime (α * d) (n + 1)) (hn : 1 ≤ n) {s : ℕ} (h1 : 1 ≤ s) (h2 : s ≤ n) :
    1 ≤ Yf n α d s := by
  by_contra h
  have h0 : Yf n α d s = Yf n α d 0 := by rw [Yf_zero]; omega
  have := Yf_inj hc hn h2 (Nat.zero_le n) h0
  omega

end arith

theorem coprime_ad {n α β d L : ℕ} (hN : α + β = n + 1) (hcop : Nat.Coprime α β)
    (hdL : d * L = n) : Nat.Coprime (α * d) (n + 1) := by
  have h1 : Nat.Coprime α (n + 1) := by rw [← hN]; exact Nat.coprime_self_add_right.mpr hcop
  have h2 : Nat.Coprime d (n + 1) := by
    rw [← hdL, show d * L + 1 = 1 + L * d by ring, Nat.coprime_add_mul_right_right]
    exact Nat.coprime_one_right d
  exact Nat.Coprime.symm (Nat.Coprime.mul_right h1.symm h2.symm)

theorem L_ad {n α d L : ℕ} (hdL : d * L = n) :
    (L : ZMod (n + 1)) * α * d = -(α : ZMod (n + 1)) := by
  have h1 : ((d * L : ℕ) : ZMod (n + 1)) = (n : ZMod (n + 1)) := by rw [hdL]
  push_cast at h1
  rw [n_cast] at h1
  linear_combination (α : ZMod (n + 1)) * h1

section prop5

variable {n α β d L : ℕ}

/-- Step 1: if no run R_i (i ≤ d - 2) meets I, then α ∣ L. -/
theorem prop5_case1 (hα : 3 ≤ α) (hαn : 2 * α < n + 1) (hc : Nat.Coprime (α * d) (n + 1))
    (hdL : d * L = n)
    (hd : 3 ≤ d) (hL : 3 ≤ L) (hnone : ∀ s, 1 ≤ s → s + L ≤ n → ¬ Yf n α d s < α) :
    α ∣ L := by
  have hn3 : 3 * L ≤ n := by rw [← hdL]; exact Nat.mul_le_mul_right L hd
  obtain ⟨ε, hεlt, hε0, hε⟩ := exists_inv hc (by omega)
  obtain ⟨τ, hτ⟩ : ∃ τ, τ = n + 1 - ε := ⟨_, rfl⟩
  have hτcast : (τ : ZMod (n + 1)) = -(ε : ZMod (n + 1)) := by
    rw [hτ, Nat.cast_sub (by omega)]; push_cast; rw [n_cast]; ring
  have key : ∀ j, j * τ ≤ n → Yf n α d (n + 1 - j * τ) = j % (n + 1) := by
    intro j hj
    apply Yf_eq (Nat.mod_lt _ (by omega))
    rw [ZMod.natCast_mod, Nat.cast_sub (by omega)]
    push_cast
    rw [hτcast, n_cast]
    linear_combination (j : ZMod (n + 1)) * hε
  have hYL : Yf n α d (n + 1 - L) = α := by
    apply Yf_eq (by omega)
    rw [Nat.cast_sub (by omega)]
    push_cast
    linear_combination ((α : ZMod (n + 1)) * d) * n_cast n - L_ad (α := α) hdL
  have hτ1 : 1 ≤ τ := by omega
  have hτn : τ ≤ n := by omega
  have step : ∀ j u, u = j * τ → 1 ≤ j → j < α → u ≤ n → u ≤ L - 1 := by
    intro j u hu hj1 hjα hjn
    rw [hu] at hjn ⊢
    have hY := key j hjn
    rw [Nat.mod_eq_of_lt (by omega)] at hY
    have hpos : 1 ≤ n + 1 - j * τ := by omega
    have h1 : ¬ (n + 1 - j * τ + L ≤ n) := fun h => hnone _ hpos h (by rw [hY]; exact hjα)
    have h2 : j * τ ≠ L := by
      intro h; rw [h, hYL] at hY; omega
    omega
  have hall : ∀ j, 1 ≤ j → j < α → j * τ ≤ L - 1 := by
    intro j hj1 hjα
    induction j with
    | zero => omega
    | succ j ih =>
      rcases Nat.eq_zero_or_pos j with h0 | h0
      · subst h0
        exact step 1 (1 * τ) rfl le_rfl hjα (by omega)
      · have h1 := ih h0 (by omega)
        have hb := step 1 (1 * τ) rfl le_rfl (by omega) (by omega)
        refine step (j + 1) ((j + 1) * τ) rfl hj1 hjα ?_
        rw [Nat.succ_mul]; omega
  have hfin : α * τ = L := by
    have h1 := hall (α - 1) (by omega) (by omega)
    have h0 := hall 1 le_rfl (by omega)
    rw [Nat.one_mul] at h0
    have e1 : α * τ = (α - 1) * τ + τ := by
      rw [← Nat.succ_mul, Nat.succ_eq_add_one, Nat.sub_add_cancel (by omega)]
    have hb : α * τ ≤ n := by omega
    have hY : Yf n α d (n + 1 - α * τ) = Yf n α d (n + 1 - L) := by
      rw [key α hb, Nat.mod_eq_of_lt (by omega), hYL]
    have := Yf_inj hc (by omega) (by omega) (by omega) hY
    omega
  exact ⟨τ, hfin.symm⟩

/-- The last step: if tβ mod α ≤ λ - 1 for all 1 ≤ t ≤ λ - 1 (λ ≥ 2, 3λ ≤ α), then β ≡ 1 (mod α). -/
theorem prop5_finish {α β lam : ℕ} (hl2 : 2 ≤ lam) (h3 : 3 * lam ≤ α) (hβ : β % α ≠ 0)
    (h : ∀ t, 1 ≤ t → t + 1 ≤ lam → t * β % α + 1 ≤ lam) : β % α = 1 := by
  have hα : 0 < α := by omega
  obtain ⟨r, hr⟩ : ∃ r, r = β % α := ⟨_, rfl⟩
  have hr1 : 1 ≤ r := by rw [hr]; exact Nat.pos_of_ne_zero hβ
  have hrl : r + 1 ≤ lam := by have := h 1 le_rfl hl2; rw [Nat.one_mul] at this; omega
  by_contra hne
  have hr2 : 2 ≤ r := by omega
  obtain ⟨t0, ht0⟩ : ∃ t0, t0 = (lam + r - 1) / r := ⟨_, rfl⟩
  have hdm := Nat.div_add_mod (lam + r - 1) r
  have hrem := Nat.mod_lt (lam + r - 1) (show 0 < r by omega)
  rw [← ht0] at hdm
  obtain ⟨P, hP⟩ : ∃ P, P = r * t0 := ⟨_, rfl⟩
  rw [← hP] at hdm
  have h1 : lam ≤ P := by omega
  have h2 : P ≤ lam + r - 1 := by omega
  have ht01 : 1 ≤ t0 := by
    rcases Nat.eq_zero_or_pos t0 with h0 | h0
    · rw [h0, Nat.mul_zero] at hP; omega
    · exact h0
  have ht0l : t0 + 1 ≤ lam := by
    by_contra hc
    have : 2 * t0 ≤ P := by rw [hP]; exact Nat.mul_le_mul_right t0 hr2
    omega
  have hmod : t0 * β % α = P := by
    rw [Nat.mul_mod, ← hr, Nat.mod_eq_of_lt (show t0 < α by omega), hP, Nat.mul_comm,
      Nat.mod_eq_of_lt (by rw [← hP]; omega)]
  have := h t0 ht01 ht0l
  omega

end prop5

end Riemann.Deep.LeanPrefix

