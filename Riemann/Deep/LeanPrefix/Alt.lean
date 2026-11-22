import Riemann.Deep.LeanPrefix.Inputs
import Riemann.Deep.LeanPrefix.Lemma7C
import Riemann.Deep.LeanPrefix.Prop5C
import Riemann.Deep.LeanPrefix.Small23
import Mathlib.Algebra.Field.ZMod

/-!
# Proposition 8: A_n ≤ G on the odd line with gcd(m+1, k+1) = 1 and (k+1) ∤ m

G is primitive (Prop 5) and contains the n-cycle P. If A_n ⊄ G, Jones's theorem (k = 0) gives:
(a) n prime and G affine: r_k ≠ 1 fixes k and k+1, impossible; (b) G ≤ PΓL_d(q): r_k is a
nonidentity involution fixing n - k = m + 1 > n/2 points, contradicting Lemma 7
((d, q) = (2, 4) would give n = 5); (c) n = 11 or 23: G = S_n by kernel computation.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv
open scoped LinearAlgebra.Projectivization

/-- The primitive odd-line triples with n = 11 or 23 generate S_n. -/
theorem small_line_top {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m) (hline : n = m + k + 1)
    (hg : Nat.Coprime (m + 1) (k + 1)) (h : n = 11 ∨ n = 23) : G n m k = ⊤ := by
  rcases h with rfl | rfl
  · have hk4 : k ≤ 4 := by omega
    interval_cases k
    · obtain rfl : m = 8 := by omega
      exact absurd hg (by decide)
    · obtain rfl : m = 7 := by omega
      exact absurd hg (by decide)
    · obtain rfl : m = 6 := by omega
      exact G_11_6_4
  · have hk10 : k ≤ 10 := by omega
    interval_cases k
    · obtain rfl : m = 20 := by omega
      exact absurd hg (by decide)
    · obtain rfl : m = 19 := by omega
      exact absurd hg (by decide)
    · obtain rfl : m = 18 := by omega
      exact G_23_18_4
    · obtain rfl : m = 17 := by omega
      exact absurd hg (by decide)
    · obtain rfl : m = 16 := by omega
      exact G_23_16_6
    · obtain rfl : m = 15 := by omega
      exact absurd hg (by decide)
    · obtain rfl : m = 14 := by omega
      exact absurd hg (by decide)
    · obtain rfl : m = 13 := by omega
      exact absurd hg (by decide)
    · obtain rfl : m = 12 := by omega
      exact G_23_12_10

theorem rk_fix {n k : ℕ} (x : Fin n) (hx : k ≤ x.val) : rev n k x = x :=
  Fin.ext (by rw [rev_val, ite_eq_right (by omega)])

/-- Proposition 8: A_n ≤ G on the odd line with gcd = 1 and (k+1) ∤ m. -/
theorem alt_le_line (hJ : Jones2014_Thm1_2_k0) {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m)
    (hline : n = m + k + 1) (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1))
    (hdiv : ¬ (k + 1) ∣ m) : alternatingGroup (Fin n) ≤ G n m k := by
  by_contra hA
  have hprim := isPreprimitive_of_line hk hkm hline hodd hg hdiv
  have hcyc : ∃ c ∈ G n m k, Perm.IsCycle c ∧ c.support = Finset.univ :=
    ⟨P n m k, P_mem n m k, P_isCycle hline hkm (Nat.Coprime.symm hg), P_support hline hkm⟩
  have hfix : n - k ≤ nfix (rev n k) := fix_rk_card (by omega)
  have hne := rk_ne_one hk (show k ≤ n by omega)
  rcases hJ n (G n m k) hprim hA hcyc with ⟨hp, e, he⟩ | ⟨q, d, F, hF, hFt, hq, hd, hn, e, he⟩ | h
  · -- (a) affine
    obtain ⟨a, b, ha, hab⟩ := he (rev n k) (rk_mem n m k)
    have : Fact (Nat.Prime n) := ⟨hp⟩
    let x1 : Fin n := ⟨k, by omega⟩
    let x2 : Fin n := ⟨k + 1, by omega⟩
    have e1 := hab x1
    have e2 := hab x2
    rw [rk_fix x1 (by simp [x1])] at e1
    rw [rk_fix x2 (by simp [x2])] at e2
    have hne12 : e x1 ≠ e x2 := e.injective.ne (by simp [x1, x2, Fin.ext_iff])
    have ha1 : a = 1 := by
      have h0 : (a - 1) * (e x1 - e x2) = 0 := by linear_combination e2 - e1
      rcases mul_eq_zero.mp h0 with h | h
      · exact sub_eq_zero.mp h
      · exact absurd (sub_eq_zero.mp h) hne12
    have hb0 : b = 0 := by rw [ha1] at e1; linear_combination -e1
    apply hne
    ext z
    have := hab z
    rw [ha1, hb0, one_mul, add_zero] at this
    exact congrArg Fin.val (e.injective this)
  · -- (b) collineation group: Lemma 7
    classical
    have h24 : ¬ (d = 2 ∧ Fintype.card F = 4) := by
      rintro ⟨rfl, h4⟩
      rw [← hq, h4] at hn
      norm_num at hn
      omega
    have := lemma7_invol hd e (rev_rev n k) hne (he _ (rk_mem n m k)) h24
    omega
  · -- (c) n = 11 or 23
    exact hA (by rw [small_line_top hk hkm hline hg h]; exact le_top)

end Riemann.Deep.LeanPrefix

