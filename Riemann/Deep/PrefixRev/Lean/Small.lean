import Riemann.Deep.PrefixRev.Lean.IsCycleP
import Riemann.Deep.PrefixRev.Lean.Parity

/-!
# Kernel-checked counterexamples to the literal Conjecture 4 (G = S_n with m + k < n)

A word in the generators is coded by a natural number c (base 3, lowest digit applied first:
0 = r_n, 1 = r_m, 2 = r_k) and a length f. The kernel evaluates the word on every point
(decide +kernel: the kernel itself evaluates; no compiled evaluation) and checks that it is the transposition (0, P 0); then
G = S_n by IsCycleP.G_eq_top_of_swap. Words found by .scratch/prefixlean/words.py.
-/

namespace PrefixRev

open Equiv

/-- Generator number d: 0 = r_n, 1 = r_m, otherwise r_k. -/
def gen (n m k d : ℕ) : Perm (Fin n) := if d = 0 then rev n n else if d = 1 then rev n m else rev n k

def sel (n m k d : ℕ) : ℕ := if d = 0 then n else if d = 1 then m else k

/-- The word coded by (f, c) as a permutation (lowest base-3 digit applied first). -/
def codePerm (n m k : ℕ) : ℕ → ℕ → Perm (Fin n)
  | 0, _ => 1
  | f + 1, c => codePerm n m k f (c / 3) * gen n m k (c % 3)

def revN (i x : ℕ) : ℕ := if x < i then i - 1 - x else x

/-- The same word evaluated on plain naturals (what the kernel computes). -/
def evalCode (n m k : ℕ) : ℕ → ℕ → ℕ → ℕ
  | 0, _, x => x
  | f + 1, c, x => evalCode n m k f (c / 3) (revN (sel n m k (c % 3)) x)

def swapN (a b x : ℕ) : ℕ := if x = a then b else if x = b then a else x

theorem gen_mem (n m k d : ℕ) : gen n m k d ∈ G n m k := by
  unfold gen; split_ifs
  · exact rn_mem n m k
  · exact rm_mem n m k
  · exact rk_mem n m k

theorem codePerm_mem (n m k : ℕ) : ∀ f c, codePerm n m k f c ∈ G n m k
  | 0, _ => (G n m k).one_mem
  | f + 1, c => (G n m k).mul_mem (codePerm_mem n m k f (c / 3)) (gen_mem n m k (c % 3))

theorem gen_val {n m k : ℕ} (hm : m ≤ n) (hk : k ≤ n) (d : ℕ) (x : Fin n) :
    (gen n m k d x).val = revN (sel n m k d) x.val := by
  unfold gen sel revN
  split_ifs <;> rw [rev_val] <;> split_ifs <;> omega

theorem codePerm_val {n m k : ℕ} (hm : m ≤ n) (hk : k ≤ n) :
    ∀ f c (x : Fin n), (codePerm n m k f c x).val = evalCode n m k f c x.val
  | 0, _, _ => rfl
  | f + 1, c, x => by
    rw [codePerm, evalCode, Perm.mul_apply, codePerm_val hm hk f, gen_val hm hk]

theorem codePerm_eq_swap {n m k f c a b : ℕ} (hm : m ≤ n) (hk : k ≤ n) (ha : a < n) (hb : b < n)
    (h : (List.range n).all (fun x => evalCode n m k f c x == swapN a b x) = true) :
    codePerm n m k f c = swap ⟨a, ha⟩ ⟨b, hb⟩ := by
  ext x
  have := List.all_eq_true.mp h x.val (List.mem_range.mpr x.isLt)
  rw [beq_iff_eq] at this
  rw [codePerm_val hm hk, this, swap_val]
  rfl

theorem check_7_4_2 : (List.range 7).all
    (fun x => evalCode 7 4 2 19 231448314 x == swapN 0 3 x) = true := by decide +kernel

theorem check_11_6_4 : (List.range 11).all
    (fun x => evalCode 11 6 4 56 98175891590607452020817016 x == swapN 0 5 x) = true := by decide +kernel

theorem check_13_10_2 : (List.range 13).all
    (fun x => evalCode 13 10 2 67 18456050534072274988258016225154 x == swapN 0 3 x) = true := by decide +kernel

/-- (7, 4, 2): r_7, r_4, r_2 generate S_7 although 4 + 2 < 7. -/
theorem G_7_4_2 : G 7 4 2 = ⊤ := by
  have hP : P 7 4 2 ⟨0, by omega⟩ = ⟨3, by omega⟩ := Fin.ext (by rw [P_val (by omega) (by omega)]; rfl)
  refine G_eq_top_of_swap (by omega) (by omega) (by decide) ⟨0, by omega⟩ ?_
  rw [hP, ← codePerm_eq_swap (by omega) (by omega) (by omega) (by omega) check_7_4_2]
  exact codePerm_mem _ _ _ _ _

/-- (11, 6, 4): r_11, r_6, r_4 generate S_11 although 6 + 4 < 11. -/
theorem G_11_6_4 : G 11 6 4 = ⊤ := by
  have hP : P 11 6 4 ⟨0, by omega⟩ = ⟨5, by omega⟩ :=
    Fin.ext (by rw [P_val (by omega) (by omega)]; rfl)
  refine G_eq_top_of_swap (by omega) (by omega) (by decide) ⟨0, by omega⟩ ?_
  rw [hP, ← codePerm_eq_swap (by omega) (by omega) (by omega) (by omega) check_11_6_4]
  exact codePerm_mem _ _ _ _ _

/-- (13, 10, 2): r_13, r_10, r_2 generate S_13 although 10 + 2 < 13. -/
theorem G_13_10_2 : G 13 10 2 = ⊤ := by
  have hP : P 13 10 2 ⟨0, by omega⟩ = ⟨3, by omega⟩ :=
    Fin.ext (by rw [P_val (by omega) (by omega)]; rfl)
  refine G_eq_top_of_swap (by omega) (by omega) (by decide) ⟨0, by omega⟩ ?_
  rw [hP, ← codePerm_eq_swap (by omega) (by omega) (by omega) (by omega) check_13_10_2]
  exact codePerm_mem _ _ _ _ _

/-- Conjecture 4 of Blanco et al. (arXiv 2511.16959v3), literal form:
2 ≤ k < m < n and m + k < n imply that r_n, r_m, r_k do not generate S_n. -/
def Conjecture4Literal : Prop :=
  ∀ n m k : ℕ, 2 ≤ k → k < m → m < n → m + k < n → G n m k ≠ ⊤

/-- The literal Conjecture 4 is false: (7, 4, 2). -/
theorem not_conjecture4Literal : ¬ Conjecture4Literal := fun h =>
  h 7 4 2 (by omega) (by omega) (by omega) (by omega) G_7_4_2

end PrefixRev
