import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# Fasi–Higham, m = 3: the real root d of d³ + 3Nd = 2N in (1/2, 2/3)

Proved from the completeness of ℝ (sSup of {x ≤ 2/3 : x³ + 3Nx ≤ 2N}) with explicit polynomial
continuity estimates (below, above), instead of the intermediate value theorem.
-/

namespace FasiLean

theorem F_mono (N x y : ℝ) (hN : 0 ≤ N) (h : x ≤ y) : x ^ 3 + 3 * N * x ≤ y ^ 3 + 3 * N * y := by
  have e : (y ^ 3 + 3 * N * y) - (x ^ 3 + 3 * N * x)
      = (y - x) * ((x + y) ^ 2 / 2 + x ^ 2 / 2 + y ^ 2 / 2 + 3 * N) := by ring
  have : 0 ≤ (y - x) * ((x + y) ^ 2 / 2 + x ^ 2 / 2 + y ^ 2 / 2 + 3 * N) :=
    mul_nonneg (by linarith) (by positivity)
  linarith

theorem below (N s : ℝ) (hN : 1 ≤ N) (hs1 : 1 / 2 ≤ s) (hs2 : s ≤ 2 / 3)
    (h : 2 * N < s ^ 3 + 3 * N * s) : ∃ δ : ℝ, 0 < δ ∧ 2 * N < (s - δ) ^ 3 + 3 * N * (s - δ) := by
  obtain ⟨δ, hδ⟩ : ∃ δ : ℝ, δ = min (1 / 4) ((s ^ 3 + 3 * N * s - 2 * N) / (6 * (N + 1))) :=
    ⟨_, rfl⟩
  have hN1 : 0 < 6 * (N + 1) := by linarith
  have hδ0 : 0 < δ := by
    rw [hδ]; exact lt_min (by norm_num) (div_pos (by linarith) hN1)
  have hδ1 : δ ≤ 1 / 4 := by rw [hδ]; exact min_le_left _ _
  have hδ2 : δ * (6 * (N + 1)) ≤ s ^ 3 + 3 * N * s - 2 * N := by
    have := min_le_right (1 / 4 : ℝ) ((s ^ 3 + 3 * N * s - 2 * N) / (6 * (N + 1)))
    rw [← hδ, le_div_iff₀ hN1] at this
    exact this
  have e : (s ^ 3 + 3 * N * s) - ((s - δ) ^ 3 + 3 * N * (s - δ))
      = δ * (3 * s ^ 2 - 3 * s * δ + δ ^ 2 + 3 * N) := by ring
  have p1 : 0 ≤ (2 / 3 - s) * (s - 1 / 2) := mul_nonneg (by linarith) (by linarith)
  have p2 : 0 ≤ (s - 1 / 2) * δ := mul_nonneg (by linarith) hδ0.le
  have p3 : 0 ≤ (1 / 4 - δ) * δ := mul_nonneg (by linarith) hδ0.le
  have hq : 3 * s ^ 2 - 3 * s * δ + δ ^ 2 ≤ 3 := by linarith
  have : δ * (3 * s ^ 2 - 3 * s * δ + δ ^ 2 + 3 * N) ≤ δ * (3 + 3 * N) :=
    mul_le_mul_of_nonneg_left (by linarith) hδ0.le
  exact ⟨δ, hδ0, by linarith⟩

theorem above (N s : ℝ) (hN : 1 ≤ N) (hs1 : 1 / 2 ≤ s) (hs3 : s < 2 / 3)
    (h : s ^ 3 + 3 * N * s < 2 * N) :
    ∃ δ : ℝ, 0 < δ ∧ s + δ ≤ 2 / 3 ∧ (s + δ) ^ 3 + 3 * N * (s + δ) ≤ 2 * N := by
  obtain ⟨δ, hδ⟩ : ∃ δ : ℝ, δ = min (2 / 3 - s) ((2 * N - (s ^ 3 + 3 * N * s)) / (6 * (N + 1))) :=
    ⟨_, rfl⟩
  have hN1 : 0 < 6 * (N + 1) := by linarith
  have hδ0 : 0 < δ := by
    rw [hδ]; exact lt_min (by linarith) (div_pos (by linarith) hN1)
  have hδ1 : δ ≤ 2 / 3 - s := by rw [hδ]; exact min_le_left _ _
  have hδ2 : δ * (6 * (N + 1)) ≤ 2 * N - (s ^ 3 + 3 * N * s) := by
    have := min_le_right (2 / 3 - s) ((2 * N - (s ^ 3 + 3 * N * s)) / (6 * (N + 1)))
    rw [← hδ, le_div_iff₀ hN1] at this
    exact this
  have e : ((s + δ) ^ 3 + 3 * N * (s + δ)) - (s ^ 3 + 3 * N * s)
      = δ * (3 * s ^ 2 + 3 * s * δ + δ ^ 2 + 3 * N) := by ring
  have q1 : 0 ≤ (2 / 3 - s) * (s - 1 / 2) := mul_nonneg (by linarith) (by linarith)
  have q2 : 0 ≤ (2 / 3 - s) * δ := mul_nonneg (by linarith) hδ0.le
  have q3 : 0 ≤ (1 / 6 - δ) * δ := mul_nonneg (by linarith) hδ0.le
  have hq : 3 * s ^ 2 + 3 * s * δ + δ ^ 2 ≤ 3 := by linarith
  have : δ * (3 * s ^ 2 + 3 * s * δ + δ ^ 2 + 3 * N) ≤ δ * (3 + 3 * N) :=
    mul_le_mul_of_nonneg_left (by linarith) hδ0.le
  exact ⟨δ, hδ0, by linarith, by linarith⟩

theorem exists_root (N : ℝ) (hN : 1 ≤ N) :
    ∃ d : ℝ, 1 / 2 < d ∧ d < 2 / 3 ∧ d ^ 3 + 3 * N * d = 2 * N := by
  have hN0 : 0 < N := by linarith
  let S : Set ℝ := {x | x ≤ 2 / 3 ∧ x ^ 3 + 3 * N * x ≤ 2 * N}
  have hhalf : (1 / 2 : ℝ) ∈ S := ⟨by norm_num, by norm_num; linarith⟩
  have hne : S.Nonempty := ⟨1 / 2, hhalf⟩
  have hbdd : BddAbove S := ⟨2 / 3, fun x hx => hx.1⟩
  have hs1 : 1 / 2 ≤ sSup S := le_csSup hbdd hhalf
  have hs2 : sSup S ≤ 2 / 3 := csSup_le hne (fun x hx => hx.1)
  have hle : sSup S ^ 3 + 3 * N * sSup S ≤ 2 * N := by
    rcases le_or_gt (sSup S ^ 3 + 3 * N * sSup S) (2 * N) with h | h
    · exact h
    · exfalso
      obtain ⟨δ, hδ0, hδ⟩ := below N (sSup S) hN hs1 hs2 h
      obtain ⟨y, hyS, hy⟩ := exists_lt_of_lt_csSup hne (show sSup S - δ < sSup S by linarith)
      have hm := F_mono N (sSup S - δ) y hN0.le hy.le
      linarith [hyS.2]
  have hge : 2 * N ≤ sSup S ^ 3 + 3 * N * sSup S := by
    rcases le_or_gt (2 * N) (sSup S ^ 3 + 3 * N * sSup S) with h | h
    · exact h
    · exfalso
      have hs3 : sSup S < 2 / 3 := by
        rcases lt_or_eq_of_le hs2 with h' | h'
        · exact h'
        · rw [h'] at h; norm_num at h; linarith
      obtain ⟨δ, hδ0, hδ1, hδ2⟩ := above N (sSup S) hN hs1 hs3 h
      have := le_csSup hbdd (show sSup S + δ ∈ S from ⟨hδ1, hδ2⟩)
      linarith
  refine ⟨sSup S, ?_, ?_, le_antisymm hle hge⟩
  · rcases lt_or_eq_of_le hs1 with h | h
    · exact h
    · exfalso; rw [← h] at hge; norm_num at hge; linarith
  · rcases lt_or_eq_of_le hs2 with h | h
    · exact h
    · exfalso; rw [h] at hle; norm_num at hle; linarith

end FasiLean
