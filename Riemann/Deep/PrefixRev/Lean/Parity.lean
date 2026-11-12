import Riemann.Deep.PrefixRev.Lean.Defs

/-!
# Parity of prefix reversals (decides A_n versus S_n)

sign r_i = (-1)^(i/2): r_i is a product of i/2 disjoint transpositions.
-/

namespace PrefixRev

open Equiv

/-- Reversal of the block of positions [a, a + i) (identity unless a + i ≤ n). -/
def blkFun (n a i : ℕ) (x : Fin n) : Fin n :=
  if h : a ≤ x.val ∧ x.val < a + i ∧ a + i ≤ n then ⟨2 * a + i - 1 - x.val, by omega⟩ else x

theorem blkFun_val (n a i : ℕ) (x : Fin n) : (blkFun n a i x).val =
    if a ≤ x.val ∧ x.val < a + i ∧ a + i ≤ n then 2 * a + i - 1 - x.val else x.val := by
  unfold blkFun; split_ifs <;> rfl

theorem blkFun_blkFun (n a i : ℕ) (x : Fin n) : blkFun n a i (blkFun n a i x) = x := by
  apply Fin.ext
  rw [blkFun_val, blkFun_val]
  by_cases h : a ≤ x.val ∧ x.val < a + i ∧ a + i ≤ n
  · rw [ite_eq_left h, ite_eq_left (by omega)]; omega
  · rw [ite_eq_right h, ite_eq_right h]

def blk (n a i : ℕ) : Perm (Fin n) where
  toFun := blkFun n a i
  invFun := blkFun n a i
  left_inv := blkFun_blkFun n a i
  right_inv := blkFun_blkFun n a i

theorem blk_val (n a i : ℕ) (x : Fin n) : (blk n a i x).val =
    if a ≤ x.val ∧ x.val < a + i ∧ a + i ≤ n then 2 * a + i - 1 - x.val else x.val :=
  blkFun_val n a i x

theorem blk_small (n a i : ℕ) (hi : i ≤ 1) : blk n a i = 1 := by
  ext x; rw [blk_val, Perm.one_apply]; split_ifs <;> omega

theorem swap_val {n : ℕ} (A B y : Fin n) : (swap A B y).val =
    if y.val = A.val then B.val else if y.val = B.val then A.val else y.val := by
  rw [swap_apply_def]
  by_cases h1 : y = A
  · rw [ite_eq_left h1, ite_eq_left (congrArg Fin.val h1)]
  · rw [ite_eq_right h1, ite_eq_right (fun e => h1 (Fin.ext e))]
    by_cases h2 : y = B
    · rw [ite_eq_left h2, ite_eq_left (congrArg Fin.val h2)]
    · rw [ite_eq_right h2, ite_eq_right (fun e => h2 (Fin.ext e))]

theorem blk_succ_succ (n a i : ℕ) (h : a + i + 2 ≤ n) :
    blk n a (i + 2) = swap (⟨a, by omega⟩ : Fin n) ⟨a + i + 1, by omega⟩ * blk n (a + 1) i := by
  ext x
  have hx := x.isLt
  rw [Perm.mul_apply, swap_val, blk_val, blk_val]
  dsimp only
  split_ifs <;> omega

theorem sign_blk (n i : ℕ) : ∀ a, a + i ≤ n → Perm.sign (blk n a i) = (-1) ^ (i / 2) := by
  induction i using Nat.strong_induction_on with
  | _ i ih =>
    intro a h
    rcases Nat.lt_or_ge i 2 with hi | hi
    · rw [blk_small n a i (by omega), map_one, show i / 2 = 0 by omega, pow_zero]
    · obtain ⟨j, rfl⟩ : ∃ j, i = j + 2 := ⟨i - 2, by omega⟩
      rw [blk_succ_succ n a j (by omega), map_mul, Perm.sign_swap (by rw [Ne, Fin.mk.injEq]; omega),
        ih j (by omega) (a + 1) (by omega), show (j + 2) / 2 = j / 2 + 1 by omega, pow_succ,
        mul_comm]

theorem rev_eq_blk (n i : ℕ) : rev n i = blk n 0 i := by
  ext x; rw [rev_val, blk_val]; split_ifs <;> omega

/-- sign r_i = (-1)^(i/2) for i ≤ n. -/
theorem sign_rev (n i : ℕ) (h : i ≤ n) : Perm.sign (rev n i) = (-1) ^ (i / 2) := by
  rw [rev_eq_blk]; exact sign_blk n i 0 (by omega)

theorem neg_one_pow_eq_one (j : ℕ) : ((-1 : ℤˣ) ^ j = 1) ↔ j % 2 = 0 := by
  rcases Nat.even_or_odd j with h | h
  · rw [h.neg_one_pow]; exact ⟨fun _ => Nat.even_iff.mp h, fun _ => rfl⟩
  · rw [h.neg_one_pow]
    exact ⟨fun e => absurd e (by decide), fun e => absurd e (by have := Nat.odd_iff.mp h; omega)⟩

/-- Parity: G ≤ A_n iff r_n, r_m, r_k are all even iff n/2, m/2, k/2 are all even. -/
theorem G_le_alt_iff {n m k : ℕ} (hm : m ≤ n) (hk : k ≤ n) :
    (∀ g ∈ G n m k, Perm.sign g = 1) ↔ (n / 2 % 2 = 0 ∧ m / 2 % 2 = 0 ∧ k / 2 % 2 = 0) := by
  constructor
  · intro h
    refine ⟨?_, ?_, ?_⟩
    · have := h _ (rn_mem n m k); rwa [sign_rev n n le_rfl, neg_one_pow_eq_one] at this
    · have := h _ (rm_mem n m k); rwa [sign_rev n m hm, neg_one_pow_eq_one] at this
    · have := h _ (rk_mem n m k); rwa [sign_rev n k hk, neg_one_pow_eq_one] at this
  · rintro ⟨h1, h2, h3⟩ g hg
    induction hg using Subgroup.closure_induction with
    | mem x hx =>
      simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hx
      rcases hx with rfl | rfl | rfl
      · rw [sign_rev n n le_rfl, neg_one_pow_eq_one]; exact h1
      · rw [sign_rev n m hm, neg_one_pow_eq_one]; exact h2
      · rw [sign_rev n k hk, neg_one_pow_eq_one]; exact h3
    | one => exact map_one _
    | mul x y _ _ hx hy => rw [map_mul, hx, hy, one_mul]
    | inv x _ hx => rw [map_inv, hx, inv_one]

/-- A subgroup containing A_n and an odd permutation is S_n. -/
theorem eq_top_of_alt {n : ℕ} {H : Subgroup (Perm (Fin n))}
    (hA : ∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ H) (g : Perm (Fin n)) (hg : g ∈ H)
    (hodd : Perm.sign g = -1) : H = ⊤ := by
  rw [Subgroup.eq_top_iff']
  intro σ
  rcases Int.units_eq_one_or (Perm.sign σ) with h | h
  · exact hA σ h
  · have : g⁻¹ * σ ∈ H := hA _ (by rw [map_mul, map_inv, hodd, h]; decide)
    have := H.mul_mem hg this
    rwa [← mul_assoc, mul_inv_cancel, one_mul] at this

/-- On the line n = m + k + 1 with m, k even: all three generators are even iff
n ≡ 1 (mod 4) and 4 | k. -/
theorem parity_line {n m k : ℕ} (hline : n = m + k + 1) (hm2 : m % 2 = 0) (hk2 : k % 2 = 0) :
    (n / 2 % 2 = 0 ∧ m / 2 % 2 = 0 ∧ k / 2 % 2 = 0) ↔ (n % 4 = 1 ∧ k % 4 = 0) := by
  omega

end PrefixRev
