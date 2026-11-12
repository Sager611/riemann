import Riemann.Deep.DPCUncond.Lean.Tail
import Mathlib.Algebra.Order.Floor.Semiring

/-!
# dpcuncond Step 6f, part 2 (new step (d)): Φ(ε) → 0 as ε → 0

Phi ex lg C θ κ E c ε M = (A + 2)^C q^{A+1} (3ε)^E + Σ_{A<a≤M} (a + 2)^C q ρ^a with q = e^θ, ρ = e^{θ−κE},
A = A(ε) = ⌊c log(1/ε)⌋ (Repair 6), (3ε)^E = e^{log(3ε)·E}; the a-sum is truncated at M (at height T the
hairpin classes S_a are empty for a ≥ L/2, since every depth is < 1/2). Phi_small: for θ < κE, 0 < c, cθ < E,
E ≥ 0: for every τ > 0 there is ε₁ > 0 with Φ(ε) ≤ τ for all 0 < ε ≤ ε₁ and all M.
-/

namespace DPCLean

open Finset

variable {ex lg : ℝ → ℝ}

/-- A(ε) = ⌊c log(1/ε)⌋ (PROOF.md Repair 6). -/
noncomputable def Aeps (lg : ℝ → ℝ) (c ε : ℝ) : ℕ := ⌊c * lg (1 / ε)⌋₊

/-- Φ(ε) of PROOF 6f, a-sum truncated at M. -/
noncomputable def Phi (ex lg : ℝ → ℝ) (C : ℕ) (θ κ E c ε : ℝ) (M : ℕ) : ℝ :=
  ((Aeps lg c ε : ℝ) + 2) ^ C * ex θ ^ (Aeps lg c ε + 1) * ex (lg (3 * ε) * E)
    + ∑ a ∈ Ioc (Aeps lg c ε) M, ((a : ℝ) + 2) ^ C * ex θ * ex (θ - κ * E) ^ a

/-- The first term of Φ: (A + 2)^C e^{θ(A+1)} (3ε)^E ≤ τ once N ≤ c log(1/ε) (needs cθ < E). -/
lemma first_small (h : ExpLog ex lg) (C : ℕ) {θ E c : ℝ} (hE : 0 ≤ E) (hc : 0 < c)
    (hcθ : c * θ < E) : ∀ τ : ℝ, 0 < τ → ∃ N : ℕ, ∀ ε : ℝ, 0 < ε → (N : ℝ) ≤ c * lg (1 / ε) →
      ((Aeps lg c ε : ℝ) + 2) ^ C * ex θ ^ (Aeps lg c ε + 1) * ex (lg (3 * ε) * E) ≤ τ := by
  intro τ hτ
  have hθE : θ - E / c < 0 := by
    have h1 : θ < E / c := by
      rw [lt_div_iff₀ hc]
      linarith
    linarith
  have ht1 : ex (θ - E / c) < 1 := h.lt_one hθE
  obtain ⟨K, s, u, hK, hs0, hs1, hu, hsu, hb⟩ := polygeom C (h.pos _).le ht1
  have hK₃0 : 0 ≤ ex θ * ex (lg 3 * E) := (mul_pos (h.pos _) (h.pos _)).le
  have hden : 0 < u * τ := mul_pos hu hτ
  obtain ⟨N, hN⟩ := exists_nat_gt (ex θ * ex (lg 3 * E) * K / (u * τ))
  refine ⟨N, fun ε hε hNε => ?_⟩
  have hcℓ0 : 0 ≤ c * lg (1 / ε) := le_trans (Nat.cast_nonneg N) hNε
  have hAle : (Aeps lg c ε : ℝ) ≤ c * lg (1 / ε) := Nat.floor_le hcℓ0
  have hNA : N ≤ Aeps lg c ε := Nat.le_floor hNε
  have hlg3 : lg (3 * ε) = lg 3 - lg (1 / ε) := by
    rw [h.lg_mul (by norm_num) hε, h.lg_inv hε]
    ring
  generalize Aeps lg c ε = A at hAle hNA ⊢
  have hexp : lg (3 * ε) * E ≤ lg 3 * E + A * (-(E / c)) := by
    rw [hlg3]
    have h1 : (A : ℝ) / c ≤ lg (1 / ε) := by
      rw [div_le_iff₀ hc]
      linarith
    have h2 : (A : ℝ) / c * E ≤ lg (1 / ε) * E := mul_le_mul_of_nonneg_right h1 hE
    have h3 : (A : ℝ) * (-(E / c)) = -((A : ℝ) / c * E) := by ring
    rw [h3]
    linarith
  have hex3 : ex (lg (3 * ε) * E) ≤ ex (lg 3 * E) * ex (-(E / c)) ^ A := by
    calc ex (lg (3 * ε) * E) ≤ ex (lg 3 * E + A * (-(E / c))) := h.mono.monotone hexp
      _ = ex (lg 3 * E) * ex (-(E / c)) ^ A := by rw [h.add, h.nat_mul]
  have ht : ex θ * ex (-(E / c)) = ex (θ - E / c) := by rw [sub_eq_add_neg, h.add]
  have hmain : ex θ ^ (A + 1) * ex (lg (3 * ε) * E) ≤ ex θ * ex (lg 3 * E) * ex (θ - E / c) ^ A := by
    calc ex θ ^ (A + 1) * ex (lg (3 * ε) * E)
        ≤ ex θ ^ (A + 1) * (ex (lg 3 * E) * ex (-(E / c)) ^ A) :=
          mul_le_mul_of_nonneg_left hex3 (pow_nonneg (h.pos _).le _)
      _ = ex θ * ex (lg 3 * E) * (ex θ * ex (-(E / c))) ^ A := by rw [pow_succ, mul_pow]; ring
      _ = ex θ * ex (lg 3 * E) * ex (θ - E / c) ^ A := by rw [ht]
  have hpg := hb A
  have h4 := pow_mul_le_one hs0.le hu.le hsu A
  have hNA' : (N : ℝ) ≤ A := by exact_mod_cast hNA
  have hKN : ex θ * ex (lg 3 * E) * K < A * (u * τ) := by
    rw [div_lt_iff₀ hden] at hN
    exact lt_of_lt_of_le hN (mul_le_mul_of_nonneg_right hNA' hden.le)
  have hsA0 : 0 ≤ s ^ A := pow_nonneg hs0.le A
  have hXu : s ^ A * (A * u) ≤ 1 := by nlinarith
  calc ((A : ℝ) + 2) ^ C * ex θ ^ (A + 1) * ex (lg (3 * ε) * E)
      = ((A : ℝ) + 2) ^ C * (ex θ ^ (A + 1) * ex (lg (3 * ε) * E)) := by ring
    _ ≤ ((A : ℝ) + 2) ^ C * (ex θ * ex (lg 3 * E) * ex (θ - E / c) ^ A) :=
        mul_le_mul_of_nonneg_left hmain (by positivity)
    _ = ex θ * ex (lg 3 * E) * (((A : ℝ) + 2) ^ C * ex (θ - E / c) ^ A) := by ring
    _ ≤ ex θ * ex (lg 3 * E) * (K * s ^ A) := mul_le_mul_of_nonneg_left hpg hK₃0
    _ = (ex θ * ex (lg 3 * E) * K) * s ^ A := by ring
    _ ≤ (A * (u * τ)) * s ^ A := mul_le_mul_of_nonneg_right hKN.le hsA0
    _ = τ * (s ^ A * (A * u)) := by ring
    _ ≤ τ * 1 := mul_le_mul_of_nonneg_left hXu hτ.le
    _ = τ := mul_one _

/-- **Φ(ε) → 0 as ε → 0** (PROOF 6f, Repair 6), uniformly in the truncation M. -/
theorem Phi_small (h : ExpLog ex lg) (C : ℕ) {θ κ E c : ℝ} (hE : 0 ≤ E) (hκ : θ < κ * E) (hc : 0 < c)
    (hcθ : c * θ < E) : ∀ τ : ℝ, 0 < τ → ∃ ε₁ : ℝ, 0 < ε₁ ∧ ∀ ε : ℝ, 0 < ε → ε ≤ ε₁ → ∀ M : ℕ,
      Phi ex lg C θ κ E c ε M ≤ τ := by
  intro τ hτ
  obtain ⟨N₁, hN₁⟩ := first_small h C hE hc hcθ (τ / 2) (by linarith)
  have hρ1 : ex (θ - κ * E) < 1 := h.lt_one (by linarith)
  obtain ⟨N₂, hN₂⟩ := tail_small C (h.pos θ).le (h.pos _).le hρ1 (τ / 2) (by linarith)
  refine ⟨1 / ex ((max N₁ N₂ : ℕ) / c), one_div_pos.2 (h.pos _), fun ε hε hε₁ M => ?_⟩
  have hexN : ex ((max N₁ N₂ : ℕ) / c) ≤ 1 / ε := by
    rw [le_div_iff₀ hε]
    rw [le_div_iff₀ (h.pos _)] at hε₁
    linarith
  have hNc : ((max N₁ N₂ : ℕ) : ℝ) ≤ c * lg (1 / ε) := by
    have h1 := h.le_lg hexN
    rw [div_le_iff₀ hc] at h1
    linarith
  have hN1 : (N₁ : ℝ) ≤ c * lg (1 / ε) :=
    le_trans (by exact_mod_cast le_max_left N₁ N₂) hNc
  have hA2 : N₂ ≤ Aeps lg c ε := le_trans (le_max_right N₁ N₂) (Nat.le_floor hNc)
  have hf := hN₁ ε hε hN1
  have ht := hN₂ (Aeps lg c ε) M hA2
  unfold Phi
  linarith

end DPCLean

