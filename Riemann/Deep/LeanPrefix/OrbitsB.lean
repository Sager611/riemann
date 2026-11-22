import Riemann.Deep.LeanPrefix.Orbits
import Riemann.Deep.LeanPrefix.IsCycleP
import Riemann.Deep.LeanPrefix.Bridge

/-!
# Theorem B.1: the orbits of G on the line n = m + k + 1 (any parity)

g = gcd(k+1, m+1). P = r_k r_n r_m is the first return of y ↦ y + (k+1) on Z/(n+1) ∖ 0
(P_firstReturn), and gcd(k+1, n+1) = g, so P^t reaches every point y with y+1 ≡ x+1 (mod g);
r_n maps x+1 to (n+1) - (x+1) ≡ -(x+1). With G_le_cls: the orbit of x is
{y : y+1 ≡ ±(x+1) (mod g)}; there are g/2 + 1 orbits; G is transitive iff g = 1.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv MulAction

/-- cls g x = cls g y iff x+1 ≡ ±(y+1) (mod g). -/
theorem cls_eq_iff {g n : ℕ} (hg : 0 < g) (x y : Fin n) :
    cls g x = cls g y ↔ ((x.val + 1) % g = (y.val + 1) % g ∨ (x.val + 1 + (y.val + 1)) % g = 0) := by
  have ha := Nat.mod_lt (x.val + 1) hg
  have hb := Nat.mod_lt (y.val + 1) hg
  have hs : (x.val + 1 + (y.val + 1)) % g = ((x.val + 1) % g + (y.val + 1) % g) % g := Nat.add_mod _ _ _
  rcases mod_cases ((x.val + 1) % g + (y.val + 1) % g) g (by omega) with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;>
    rw [hs, h2] <;> unfold cls <;> split_ifs <;> omega

theorem cls_le_half {n : ℕ} (g : ℕ) (hg : 0 < g) (x : Fin n) : 2 * cls g x ≤ g := by
  have := Nat.mod_lt (x.val + 1) hg
  unfold cls; split_ifs <;> omega

section line

variable {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
include hline hkm

/-- First return, iterated: some power of P maps x to the point with 1-based value
(x+1 + s(k+1)) mod (n+1), or, when that value is 0, to the next point of the rotation orbit. -/
theorem P_pow_reach (x : Fin n) (s : ℕ) :
    ∃ t : ℕ, ((P n m k ^ t) x).val + 1 =
      if (x.val + 1 + s * (k + 1)) % (n + 1) = 0 then (x.val + 1 + (s + 1) * (k + 1)) % (n + 1)
      else (x.val + 1 + s * (k + 1)) % (n + 1) := by
  have hx := x.isLt
  have hstep : ∀ j : ℕ, (x.val + 1 + (j + 1) * (k + 1)) % (n + 1) =
      ((x.val + 1 + j * (k + 1)) % (n + 1) + (k + 1)) % (n + 1) := by
    intro j
    rw [Nat.add_mod ((x.val + 1 + j * (k + 1)) % (n + 1)), Nat.mod_mod,
      Nat.mod_eq_of_lt (show k + 1 < n + 1 by omega), ← Nat.mod_eq_of_lt (show k + 1 < n + 1 by omega),
      ← Nat.add_mod, Nat.mod_eq_of_lt (show k + 1 < n + 1 by omega)]
    congr 1; ring
  have hstep2 : ∀ j : ℕ, (x.val + 1 + (j + 1 + 1) * (k + 1)) % (n + 1) =
      ((x.val + 1 + j * (k + 1)) % (n + 1) + 2 * (k + 1)) % (n + 1) := by
    intro j
    rw [Nat.add_mod ((x.val + 1 + j * (k + 1)) % (n + 1)), Nat.mod_mod, ← Nat.add_mod]
    congr 1; ring
  induction s with
  | zero =>
    refine ⟨0, ?_⟩
    rw [pow_zero, Perm.one_apply, Nat.zero_mul, Nat.add_zero, Nat.mod_eq_of_lt (by omega),
      ite_eq_right (by omega)]
  | succ s ih =>
    obtain ⟨t, ht⟩ := ih
    by_cases hz : (x.val + 1 + s * (k + 1)) % (n + 1) = 0
    · refine ⟨t, ?_⟩
      rw [ite_eq_left hz] at ht
      have hz1 : (x.val + 1 + (s + 1) * (k + 1)) % (n + 1) = k + 1 := by
        rw [hstep, hz, Nat.zero_add, Nat.mod_eq_of_lt (by omega)]
      rw [ite_eq_right (by omega), ht]
    · refine ⟨t + 1, ?_⟩
      rw [ite_eq_right hz] at ht
      have hP := P_firstReturn hline hkm ((P n m k ^ t) x)
      rw [ht, ← hstep, ← hstep2] at hP
      rw [pow_succ', Perm.mul_apply, hP]

/-- Every y with y+1 ≡ x+1 (mod gcd(k+1, m+1)) is P^t x for some t. -/
theorem P_reach_class (x y : Fin n)
    (h : (x.val + 1) % Nat.gcd (k + 1) (m + 1) = (y.val + 1) % Nat.gcd (k + 1) (m + 1)) :
    ∃ t : ℕ, (P n m k ^ t) x = y := by
  have hx := x.isLt; have hy := y.isLt
  have hgN : Nat.gcd (k + 1) (n + 1) = Nat.gcd (k + 1) (m + 1) := by
    rw [show n + 1 = (m + 1) + (k + 1) by omega, Nat.gcd_add_self_right]
  set g := Nat.gcd (k + 1) (m + 1) with hgdef
  have hgpos : 0 < g := Nat.gcd_pos_of_pos_left _ (by omega)
  have hgk : g ∣ k + 1 := Nat.gcd_dvd_left _ _
  have hgle : g ≤ k + 1 := Nat.le_of_dvd (by omega) hgk
  have hgn : g ∣ n + 1 := by rw [← hgN]; exact Nat.gcd_dvd_right _ _
  obtain ⟨u, -, hu⟩ := Nat.exists_mul_mod_eq_gcd (k := n + 1) (n := k + 1) (by rw [hgN]; omega)
  rw [hgN] at hu
  -- D ≡ y - x (mod n+1), g ∣ D
  set D := y.val + 1 + (n + 1) - (x.val + 1) with hD
  have hgD : g ∣ D := by
    apply Nat.dvd_of_mod_eq_zero
    have e1 : (y.val + 1 + (n + 1)) % g = (x.val + 1) % g := by
      rw [Nat.add_mod, Nat.mod_eq_zero_of_dvd hgn, Nat.add_zero, Nat.mod_mod, h]
    exact Nat.sub_mod_eq_zero_of_mod_eq e1
  obtain ⟨c, hc⟩ := hgD
  obtain ⟨t, ht⟩ := P_pow_reach hline hkm x (u * c)
  have hval : (x.val + 1 + u * c * (k + 1)) % (n + 1) = y.val + 1 := by
    have e2 : u * c * (k + 1) % (n + 1) = D % (n + 1) := by
      rw [show u * c * (k + 1) = (k + 1) * u * c by ring, ← Nat.mod_mul_mod, hu, ← hc]
    rw [Nat.add_mod, e2, ← Nat.add_mod, show x.val + 1 + D = (y.val + 1) + (n + 1) by omega,
      Nat.add_mod_right, Nat.mod_eq_of_lt (by omega)]
  rw [hval, ite_eq_right (by omega)] at ht
  exact ⟨t, Fin.ext (by omega)⟩

theorem rn_add (x : Fin n) : (rev n n x).val + 1 + (x.val + 1) = n + 1 := by
  have hx := x.isLt
  rw [rev_val, ite_eq_left ⟨hx, le_rfl⟩]; omega

/-- Theorem B.1, the orbits: y is in the G-orbit of x iff y+1 ≡ ±(x+1) (mod gcd(k+1, m+1)). -/
theorem mem_orbit_iff_mod (x y : Fin n) :
    y ∈ MulAction.orbit (G n m k) x ↔
      ((y.val + 1) % Nat.gcd (k + 1) (m + 1) = (x.val + 1) % Nat.gcd (k + 1) (m + 1) ∨
        (y.val + 1 + (x.val + 1)) % Nat.gcd (k + 1) (m + 1) = 0) := by
  have hgpos : 0 < Nat.gcd (k + 1) (m + 1) := Nat.gcd_pos_of_pos_left _ (by omega)
  constructor
  · rintro ⟨h, rfl⟩
    have := G_le_cls hline hkm h.2 x
    exact (cls_eq_iff hgpos _ _).mp this
  · rintro (h | h)
    · obtain ⟨t, ht⟩ := P_reach_class hline hkm x y h.symm
      exact ⟨⟨P n m k ^ t, (G n m k).pow_mem (P_mem n m k) t⟩, ht⟩
    · have hgN : Nat.gcd (k + 1) (m + 1) ∣ n + 1 := by
        rw [hline, show m + k + 1 + 1 = (m + 1) + (k + 1) by omega]
        exact Nat.dvd_add (Nat.gcd_dvd_right _ _) (Nat.gcd_dvd_left _ _)
      have hr := rn_add hline hkm x
      have h2 : ((rev n n x).val + 1) % Nat.gcd (k + 1) (m + 1) =
          (y.val + 1) % Nat.gcd (k + 1) (m + 1) := by
        have e1 : ((rev n n x).val + 1 + (x.val + 1)) % Nat.gcd (k + 1) (m + 1) = 0 := by
          rw [hr]; exact Nat.mod_eq_zero_of_dvd hgN
        have h3 := Nat.mod_lt (x.val + 1) hgpos
        have h4 := Nat.mod_lt (y.val + 1) hgpos
        have h5 := Nat.mod_lt ((rev n n x).val + 1) hgpos
        rw [Nat.add_mod] at e1 h
        rcases mod_cases (((rev n n x).val + 1) % Nat.gcd (k + 1) (m + 1) +
            (x.val + 1) % Nat.gcd (k + 1) (m + 1)) (Nat.gcd (k + 1) (m + 1)) (by omega) with
          ⟨_, e3⟩ | ⟨_, e3⟩ <;>
        rcases mod_cases ((y.val + 1) % Nat.gcd (k + 1) (m + 1) +
            (x.val + 1) % Nat.gcd (k + 1) (m + 1)) (Nat.gcd (k + 1) (m + 1)) (by omega) with
          ⟨_, e4⟩ | ⟨_, e4⟩ <;> omega
      obtain ⟨t, ht⟩ := P_reach_class hline hkm (rev n n x) y h2
      refine ⟨⟨P n m k ^ t * rev n n, (G n m k).mul_mem ((G n m k).pow_mem (P_mem n m k) t)
        (rn_mem n m k)⟩, ?_⟩
      show (P n m k ^ t * rev n n) x = y
      rw [Perm.mul_apply]; exact ht

/-- Theorem B.1, the number of orbits: g/2 + 1. -/
theorem card_orbits :
    Nat.card (MulAction.orbitRel.Quotient (G n m k) (Fin n)) = Nat.gcd (k + 1) (m + 1) / 2 + 1 := by
  set g := Nat.gcd (k + 1) (m + 1) with hgdef
  have hgpos : 0 < g := Nat.gcd_pos_of_pos_left _ (by omega)
  let Φ : MulAction.orbitRel.Quotient (G n m k) (Fin n) → Fin (g / 2 + 1) :=
    Quotient.lift (fun x => ⟨cls g x, by have := cls_le_half g hgpos x; omega⟩) (by
      intro a b hab
      have hab' : a ∈ MulAction.orbit (G n m k) b := (MulAction.orbitRel_apply).mp hab
      obtain ⟨h, rfl⟩ := hab'
      apply Fin.ext
      exact G_le_cls hline hkm h.2 b)
  have hinj : Function.Injective Φ := by
    intro a b hab
    induction a using Quotient.inductionOn with
    | h a =>
      induction b using Quotient.inductionOn with
      | h b =>
        apply Quotient.sound
        change MulAction.orbitRel (G n m k) (Fin n) a b
        rw [MulAction.orbitRel_apply, mem_orbit_iff_mod hline hkm]
        have h1 : cls g a = cls g b := congrArg Fin.val hab
        rcases (cls_eq_iff hgpos a b).mp h1 with h2 | h2
        · exact Or.inl h2
        · exact Or.inr h2
  have hsurj : Function.Surjective Φ := by
    intro c
    obtain ⟨x, hx⟩ := cls_surj hline c.val (by have := c.isLt; omega)
    exact ⟨Quotient.mk _ x, Fin.ext hx⟩
  rw [Nat.card_congr (Equiv.ofBijective Φ ⟨hinj, hsurj⟩), Nat.card_eq_fintype_card, Fintype.card_fin]

/-- Theorem B.1: G is transitive iff gcd(k+1, m+1) = 1. -/
theorem isPretransitive_iff_gcd :
    MulAction.IsPretransitive (G n m k) (Fin n) ↔ Nat.gcd (k + 1) (m + 1) = 1 := by
  have hgpos : 0 < Nat.gcd (k + 1) (m + 1) := Nat.gcd_pos_of_pos_left _ (by omega)
  have hgle : Nat.gcd (k + 1) (m + 1) ≤ k + 1 := Nat.le_of_dvd (by omega) (Nat.gcd_dvd_left _ _)
  constructor
  · intro ht
    obtain ⟨h, hh⟩ := ht.exists_smul_eq (⟨0, by omega⟩ : Fin n)
      ⟨Nat.gcd (k + 1) (m + 1) - 1, by omega⟩
    have hmem := (mem_orbit_iff_mod hline hkm (⟨0, by omega⟩ : Fin n)
      ⟨Nat.gcd (k + 1) (m + 1) - 1, by omega⟩).mp ⟨h, hh⟩
    simp only at hmem
    rw [show Nat.gcd (k + 1) (m + 1) - 1 + 1 = Nat.gcd (k + 1) (m + 1) by omega, Nat.mod_self] at hmem
    rcases hmem with h1 | h1
    · by_contra hne
      rw [Nat.mod_eq_of_lt (show 0 + 1 < Nat.gcd (k + 1) (m + 1) by omega)] at h1
      omega
    · rw [Nat.add_mod, Nat.mod_self, Nat.zero_add, Nat.mod_mod] at h1
      by_contra hne
      rw [Nat.mod_eq_of_lt (show 0 + 1 < Nat.gcd (k + 1) (m + 1) by omega)] at h1
      omega
  · intro h1
    refine ⟨fun x y => ?_⟩
    have := (mem_orbit_iff_mod hline hkm x y).mpr (Or.inl (by rw [h1, Nat.mod_one, Nat.mod_one]))
    exact this

end line

end Riemann.Deep.LeanPrefix

