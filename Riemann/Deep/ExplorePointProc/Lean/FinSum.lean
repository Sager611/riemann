import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Step (iii') of explore40's Theorem 1: the hard-core pair-sum inequality (DERIVATION §4)

Points x 0, ..., x (n-1) with all gaps x (i+1) - x i ≥ λ ≥ 0, an even kernel F with F ≤ 0 on
[U, ∞), U ≤ 3λ, and (H) F x + F (x + y) + A x - A y + β (1 - x) ≤ c for x, y ≥ λ. Then
(0-indexed, gaps g i = x (i+1) - x i)

  Σ_{i,j<n} F (x i - x j) ≤ n F 0 + 2 (n-2) c + 2 β ((x (n-2) - x 0) - (n-2))
                             + 2 (A (g (n-2)) - A (g 0)) + 2 F (g (n-2)),

which is the brief's 1-indexed form (x_{n-1} - x_1, A(g_{n-1}) - A(g_1), F(g_{n-1})).
Pairs at index distance ≥ 3 are ≥ 3λ ≥ U apart (dropped, F ≤ 0); nearest and next-nearest pairs are
bounded by (H) at (g i, g (i+1)): A telescopes, the β-terms sum to β((x (n-2) - x 0) - (n-2)).
pairSum_weighted_le: the same for the weighted kernel F u * w u (w even, ≥ 0, w 0 = 1,
F u * w u ≤ F u + η), with c + 2η in place of c: Lamzouri's weight w → 1 step.
-/

namespace PPLean

open Finset

/-- Hypothesis (H) of Theorem 1 for a kernel F (= f^): on [λ, ∞)². -/
def HH (lam β c : ℝ) (F A : ℝ → ℝ) : Prop :=
  ∀ x y : ℝ, lam ≤ x → lam ≤ y → F x + F (x + y) + A x - A y + β * (1 - x) ≤ c

/-- The pair sum Σ_{i<n} Σ_{j<n} F (x i - x j). -/
def pairSum (F : ℝ → ℝ) (x : ℕ → ℝ) (n : ℕ) : ℝ :=
  ∑ i ∈ range n, ∑ j ∈ range n, F (x i - x j)

lemma pairSum_zero (F : ℝ → ℝ) (x : ℕ → ℝ) : pairSum F x 0 = 0 := by
  simp [pairSum]

/-- Adding the point x k: diagonal F 0 and twice the row Σ_{i<k} F (x k - x i) (F even). -/
lemma pairSum_succ {F : ℝ → ℝ} (heven : ∀ u, F (-u) = F u) (x : ℕ → ℝ) (k : ℕ) :
    pairSum F x (k + 1) = pairSum F x k + 2 * ∑ i ∈ range k, F (x k - x i) + F 0 := by
  have h1 : ∑ i ∈ range k, F (x i - x k) = ∑ i ∈ range k, F (x k - x i) :=
    Finset.sum_congr rfl (fun i _ => by rw [← neg_sub, heven])
  unfold pairSum
  rw [Finset.sum_range_succ]
  simp only [Finset.sum_range_succ, Finset.sum_add_distrib, sub_self]
  rw [h1]
  ring

/-- d consecutive gaps ≥ λ add up: x i + d λ ≤ x (i + d). -/
lemma gap_sum {lam : ℝ} {n : ℕ} {x : ℕ → ℝ} (hgap : ∀ i, i + 1 < n → lam ≤ x (i + 1) - x i) :
    ∀ d i : ℕ, i + d < n → x i + d * lam ≤ x (i + d) := by
  intro d
  induction d with
  | zero => intro i _; simp
  | succ d ih =>
    intro i hi
    have h1 := ih i (by omega)
    have h2 := hgap (i + d) (by omega)
    have h3 : i + (d + 1) = i + d + 1 := by omega
    rw [h3]
    push_cast
    linarith

/-- The new row of the point x (m+2): all but its two nearest left neighbours are ≥ 3λ ≥ U away. -/
lemma row_le {lam U : ℝ} {F : ℝ → ℝ} (hlam : 0 ≤ lam) (hU : U ≤ 3 * lam)
    (hneg : ∀ u, U ≤ u → F u ≤ 0) {x : ℕ → ℝ} {m : ℕ}
    (hgap : ∀ i, i + 1 < m + 3 → lam ≤ x (i + 1) - x i) :
    ∑ i ∈ range (m + 2), F (x (m + 2) - x i) ≤ F (x (m + 2) - x m) + F (x (m + 2) - x (m + 1)) := by
  rw [Finset.sum_range_succ, Finset.sum_range_succ]
  have h0 : ∑ i ∈ range m, F (x (m + 2) - x i) ≤ 0 := by
    apply Finset.sum_nonpos
    intro i hi
    rw [Finset.mem_range] at hi
    apply hneg
    have h1 := gap_sum hgap (m + 2 - i) i (by omega)
    have h2 : i + (m + 2 - i) = m + 2 := by omega
    rw [h2] at h1
    have h3 : (3 : ℝ) ≤ ((m + 2 - i : ℕ) : ℝ) := by
      have : 3 ≤ m + 2 - i := by omega
      exact_mod_cast this
    nlinarith
  linarith

/-- (iii') with n = m + 2 points, by induction on m. -/
theorem pairSum_le_aux {lam U β c : ℝ} {F A : ℝ → ℝ} (hlam : 0 ≤ lam) (hU : U ≤ 3 * lam)
    (heven : ∀ u, F (-u) = F u) (hneg : ∀ u, U ≤ u → F u ≤ 0) (hH : HH lam β c F A)
    (x : ℕ → ℝ) (m : ℕ) (hgap : ∀ i, i + 1 < m + 2 → lam ≤ x (i + 1) - x i) :
    pairSum F x (m + 2) ≤ ((m : ℝ) + 2) * F 0 + 2 * (m : ℝ) * c
      + 2 * β * ((x m - x 0) - m) + 2 * (A (x (m + 1) - x m) - A (x 1 - x 0))
      + 2 * F (x (m + 1) - x m) := by
  induction m with
  | zero =>
    have e : F (x 0 - x 1) = F (x 1 - x 0) := by rw [← neg_sub, heven]
    rw [show (0 : ℕ) + 2 = 0 + 1 + 1 from rfl, pairSum_succ heven, pairSum_succ heven, pairSum_zero]
    simp only [Finset.sum_range_zero, Finset.sum_range_one, Nat.cast_zero, zero_add, sub_self]
    linarith
  | succ m ih =>
    have h1 := ih (fun i hi => hgap i (by omega))
    have h2 := row_le hlam hU hneg (m := m) (x := x) (fun i hi => hgap i (by omega))
    have h3 := hH (x (m + 1) - x m) (x (m + 2) - x (m + 1)) (hgap m (by omega))
      (hgap (m + 1) (by omega))
    have h4 : x (m + 1) - x m + (x (m + 2) - x (m + 1)) = x (m + 2) - x m := by ring
    rw [h4] at h3
    rw [show m + 1 + 2 = (m + 2) + 1 by omega, pairSum_succ heven, show m + 1 + 1 = m + 2 by omega]
    push_cast
    linarith

/-- **Step (iii')** (DERIVATION §4), exact boundary form, n ≥ 2 points, 0-indexed. -/
theorem pairSum_le {lam U β c : ℝ} {F A : ℝ → ℝ} (hlam : 0 ≤ lam) (hU : U ≤ 3 * lam)
    (heven : ∀ u, F (-u) = F u) (hneg : ∀ u, U ≤ u → F u ≤ 0) (hH : HH lam β c F A)
    {n : ℕ} (hn : 2 ≤ n) (x : ℕ → ℝ) (hgap : ∀ i, i + 1 < n → lam ≤ x (i + 1) - x i) :
    pairSum F x n ≤ (n : ℝ) * F 0 + 2 * ((n : ℝ) - 2) * c
      + 2 * β * ((x (n - 2) - x 0) - ((n : ℝ) - 2))
      + 2 * (A (x (n - 1) - x (n - 2)) - A (x 1 - x 0)) + 2 * F (x (n - 1) - x (n - 2)) := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
  have h := pairSum_le_aux hlam hU heven hneg hH x m hgap
  rw [show m + 2 - 2 = m by omega, show m + 2 - 1 = m + 1 by omega]
  push_cast
  linarith

/-- Weighted (iii'): kernel F u * w u with w even, ≥ 0, w 0 = 1 and F u * w u ≤ F u + η;
then (H) holds for it with c + 2η (Lamzouri's weight w → 1 on near pairs; far pairs: F w ≤ 0). -/
theorem pairSum_weighted_le {lam U β c η : ℝ} {F A w : ℝ → ℝ} (hlam : 0 ≤ lam) (hU : U ≤ 3 * lam)
    (heven : ∀ u, F (-u) = F u) (hneg : ∀ u, U ≤ u → F u ≤ 0) (hH : HH lam β c F A)
    (hweven : ∀ u, w (-u) = w u) (hwnn : ∀ u, 0 ≤ w u) (hw0 : w 0 = 1)
    (hη : ∀ u, F u * w u ≤ F u + η)
    {n : ℕ} (hn : 2 ≤ n) (x : ℕ → ℝ) (hgap : ∀ i, i + 1 < n → lam ≤ x (i + 1) - x i) :
    pairSum (fun u => F u * w u) x n ≤ (n : ℝ) * F 0 + 2 * ((n : ℝ) - 2) * (c + 2 * η)
      + 2 * β * ((x (n - 2) - x 0) - ((n : ℝ) - 2))
      + 2 * (A (x (n - 1) - x (n - 2)) - A (x 1 - x 0))
      + 2 * (F (x (n - 1) - x (n - 2)) * w (x (n - 1) - x (n - 2))) := by
  have hGeven : ∀ u, (fun u => F u * w u) (-u) = (fun u => F u * w u) u := by
    intro u; simp only [heven, hweven]
  have hGneg : ∀ u, U ≤ u → (fun u => F u * w u) u ≤ 0 := by
    intro u hu
    have h1 := hneg u hu
    have h2 := hwnn u
    simp only
    nlinarith
  have hGH : HH lam β (c + 2 * η) (fun u => F u * w u) A := by
    intro x y hx hy
    have h1 := hH x y hx hy
    have h2 := hη x
    have h3 := hη (x + y)
    simp only
    linarith
  have h := pairSum_le hlam hU hGeven hGneg hGH hn x hgap
  simp only [hw0, mul_one] at h
  exact h

end PPLean
