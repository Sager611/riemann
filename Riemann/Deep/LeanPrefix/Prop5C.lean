import Riemann.Deep.LeanPrefix.Prop5B
import Riemann.Deep.LeanPrefix.PrimB2
import Riemann.Deep.LeanPrefix.Bridge

/-!
# Proposition 5: primitivity on the odd line when (k+1) ∤ m

prop5_core: (★_d) (StarD, Lemma 3's condition in the cyclic coordinate) for a divisor 1 < d < n
of n forces (k+1) ∣ m. Translation: the point x_s with x_s + 1 = Y(s) = s α d mod (n+1) has
phi x_s = s d mod (n+1) = d + t d - i for s = 1 + t + iL, so its class mod d is d - 1 - i;
StarD on run i ≤ d - 2 is prop5_arith's hypothesis. primitive_of_line: G is primitive
(PrefixRev form, then Mathlib's MulAction.IsPreprimitive).
-/

namespace Riemann.Deep.LeanPrefix

open Equiv

theorem odd_of_dvd_odd {n d : ℕ} (hodd : n % 2 = 1) (h : d ∣ n) : d % 2 = 1 := by
  by_contra hc
  have : 2 ∣ n := Nat.dvd_trans (Nat.dvd_of_mod_eq_zero (by omega)) h
  omega

/-- Proposition 5, Steps 1-4: (★_d) for a divisor 1 < d < n of n forces (k+1) ∣ m. -/
theorem prop5_core {n m k θ d : ℕ} (hk : 2 ≤ k) (hkm : k < m) (hline : n = m + k + 1)
    (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1)) (hθ : (k + 1) * θ % (n + 1) = 1)
    (hd1 : 1 < d) (hdn : d < n) (hdvd : d ∣ n) (hS : StarD n k θ d) : (k + 1) ∣ m := by
  obtain ⟨L, hL⟩ := hdvd
  have hdL : d * L = n := hL.symm
  have hd3 : 3 ≤ d := by have := odd_of_dvd_odd hodd ⟨L, hL⟩; omega
  have hL3 : 3 ≤ L := by
    have hLo := odd_of_dvd_odd hodd ⟨d, by rw [hL, Nat.mul_comm]⟩
    rcases Nat.lt_or_ge L 2 with h | h
    · interval_cases L <;> omega
    · omega
  have hN : (k + 1) + (m + 1) = n + 1 := by omega
  have hθZ : ((k + 1 : ℕ) : ZMod (n + 1)) * θ = 1 := by
    have h1 : (((k + 1) * θ % (n + 1) : ℕ) : ZMod (n + 1)) = ((1 : ℕ) : ZMod (n + 1)) := by rw [hθ]
    rw [ZMod.natCast_mod] at h1; push_cast at h1 ⊢; exact h1
  have hc := coprime_ad hN (Nat.Coprime.symm hg) hdL
  have hstar : ∀ i s s', i + 2 ≤ d → i * L < s → s ≤ (i + 1) * L → i * L < s' →
      s' ≤ (i + 1) * L → (Yf n (k + 1) d s < k + 1 ↔ Yf n (k + 1) d s' < k + 1) := by
    -- the point at time s and its class
    have hiL : ∀ i : ℕ, (i + 1) * L = i * L + L := fun i => Nat.succ_mul i L
    have hrun_le : ∀ i, i + 2 ≤ d → (i + 1) * L ≤ n := by
      intro i hi; rw [← hdL]; exact Nat.mul_le_mul_right L (by omega)
    have hpt : ∀ i s, i + 2 ≤ d → i * L < s → s ≤ (i + 1) * L →
        ∃ x : Fin n, x.val + 1 = Yf n (k + 1) d s ∧ (phi n θ x - 1) % d = d - 1 - i := by
      intro i s hi h1 h2
      have hsn := hrun_le i hi
      have hY1 := Yf_pos hc (by omega) (show 1 ≤ s by omega) (by omega)
      have hY2 := Yf_lt (n := n) (α := k + 1) (d := d) s
      refine ⟨⟨Yf n (k + 1) d s - 1, by omega⟩, by simp only; omega, ?_⟩
      obtain ⟨t, ht⟩ : ∃ t, t = s - 1 - i * L := ⟨_, rfl⟩
      have htL : t ≤ L - 1 := by have := hiL i; omega
      have hs : s = 1 + t + i * L := by omega
      have hphi : phi n θ ⟨Yf n (k + 1) d s - 1, by omega⟩ = s * d % (n + 1) := by
        unfold phi
        simp only
        rw [Nat.sub_add_cancel hY1]
        apply (ZMod.natCast_eq_natCast_iff' _ _ _).mp
        push_cast
        rw [Yf_cast]
        push_cast
        push_cast at hθZ
        linear_combination ((s : ZMod (n + 1)) * d) * hθZ
      have htd : t * d ≤ (L - 1) * d := Nat.mul_le_mul_right d htL
      have hLd : (L - 1) * d + d = n := by
        rw [← Nat.succ_mul, Nat.succ_eq_add_one, Nat.sub_add_cancel (by omega), Nat.mul_comm, hdL]
      have hsd : s * d = (d + t * d - i) + i * (n + 1) := by
        have : s * d = d + t * d + i * n := by rw [hs, ← hdL]; ring
        rw [this, Nat.mul_add, Nat.mul_one]; omega
      rw [hphi, hsd, Nat.add_mul_mod_self_right, Nat.mod_eq_of_lt (show d + t * d - i < n + 1 by omega),
        show d + t * d - i - 1 = (d - 1 - i) + t * d by omega, Nat.add_mul_mod_self_right,
        Nat.mod_eq_of_lt (by omega)]
    intro i s s' hi h1 h2 h1' h2'
    obtain ⟨x, hx, hxc⟩ := hpt i s hi h1 h2
    obtain ⟨y, hy, hyc⟩ := hpt i s' hi h1' h2'
    have hxy := hS x y (by rw [hxc, hyc]) (by rw [hxc]; omega)
    rw [← hx, ← hy]
    omega
  have := prop5_arith hN (by omega) (by omega) (Nat.Coprime.symm hg) hdL hd3 hL3 hstar
  have h2 : (k + 1) ∣ m + (k + 1) := by rw [show m + (k + 1) = n by omega]; exact this
  exact (Nat.dvd_add_left (dvd_refl (k + 1))).mp h2

/-- Proposition 5: on the odd line with gcd(m+1, k+1) = 1 and (k+1) ∤ m, G is primitive. -/
theorem primitive_of_line {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m) (hline : n = m + k + 1)
    (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1)) (hdiv : ¬ (k + 1) ∣ m) :
    IsPrimitive (G n m k) := by
  obtain ⟨θ, hθ⟩ := exists_theta hline (Nat.Coprime.symm hg)
  refine ⟨fun x y => ?_, fun B hB => ?_⟩
  · obtain ⟨t, ht⟩ := P_cycle hline hkm (Nat.Coprime.symm hg) x y
    exact ⟨P n m k ^ t, (G n m k).pow_mem (P_mem_G n m k) t, ht⟩
  · rcases Set.eq_empty_or_nonempty B with h | ⟨x0, hx0⟩
    · left; rw [h]; exact Set.subsingleton_empty
    · rcases block_classes hline hkm hθ hB hx0 with h | h | ⟨d, hd1, hdn', hdn, hmem⟩
      · exact Or.inl h
      · exact Or.inr h
      · exact absurd (prop5_core hk hkm hline hodd hg hθ hd1 hdn' hdn
          (star_of_block hline hkm hθ hodd hB hd1 hdn' hdn hmem)) hdiv

/-- Proposition 5 with Mathlib's primitivity. -/
theorem isPreprimitive_of_line {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m) (hline : n = m + k + 1)
    (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1)) (hdiv : ¬ (k + 1) ∣ m) :
    MulAction.IsPreprimitive (G n m k) (Fin n) :=
  (isPreprimitive_iff_isPrimitive _).mpr (primitive_of_line hk hkm hline hodd hg hdiv)

end Riemann.Deep.LeanPrefix

