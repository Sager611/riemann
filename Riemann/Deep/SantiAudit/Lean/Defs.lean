/-
E1Lean.Defs — definitions of Theorem E1 (Riemann/Deep/SantiAudit/E1.md §1, §4).
Elementary functions are abstract (structure `Fns`): each field is verbatim a Mathlib theorem
(intended instance `Fns.real` in Inst.lean, uncompiled: Real.sqrt, Real.sin, Real.cos, Real.exp, Real.log, Real.rpow).
-/
import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.NormNum

namespace E1Lean

/-- Elementary functions with the Mathlib facts used (field name = Mathlib theorem name, `Real.` namespace). -/
structure Fns where
  sqrt : ℝ → ℝ
  sin : ℝ → ℝ
  cos : ℝ → ℝ
  exp : ℝ → ℝ
  log : ℝ → ℝ
  /-- `rpow T θ = T ^ θ` (Real.rpow); only a symbol for the window H = T^θ. -/
  rpow : ℝ → ℝ → ℝ
  sqrt_nonneg : ∀ x, 0 ≤ sqrt x
  mul_self_sqrt : ∀ x, 0 ≤ x → sqrt x * sqrt x = x
  sin_bound : ∀ x, |x| ≤ 1 → |sin x - (x - x ^ 3 / 6)| ≤ |x| ^ 5 / 100
  cos_bound : ∀ x, |x| ≤ 1 → |cos x - (1 - x ^ 2 / 2)| ≤ |x| ^ 4 * (5 / 96)
  sin_two_mul : ∀ x, sin (2 * x) = 2 * sin x * cos x
  cos_two_mul : ∀ x, cos (2 * x) = 2 * cos x ^ 2 - 1
  sin_sq_add_cos_sq : ∀ x, sin x ^ 2 + cos x ^ 2 = 1
  sin_sub : ∀ x y, sin (x - y) = sin x * cos y - cos x * sin y
  sin_lt : ∀ x, 0 < x → sin x < x
  sin_gt_sub_cube : ∀ x, 0 < x → x - x ^ 3 / 6 < sin x
  exp_pos : ∀ x, 0 < exp x
  exp_zero : exp 0 = 1
  exp_add : ∀ x y, exp (x + y) = exp x * exp y
  exp_one_gt_d9 : 2.7182818283 < exp 1
  quadratic_le_exp_of_nonneg : ∀ x, 0 ≤ x → 1 + x + x ^ 2 / 2 ≤ exp x
  log_le_iff_le_exp : ∀ x y, 0 < x → (log x ≤ y ↔ x ≤ exp y)
  exp_one_lt_d9 : exp 1 < 2.7182818286
  exp_bound_div_one_sub_of_interval : ∀ x, 0 ≤ x → x < 1 → exp x ≤ 1 / (1 - x)

variable (F : Fns)

/-- cot x = cos x / sin x (Mathlib: Real.cot_eq_cos_div_sin). -/
noncomputable def cot (x : ℝ) : ℝ := F.cos x / F.sin x

/-- Wang's density c(θ) = 2 − θ/2 − cot(θ/√2)/√2 (E1.md §1; SL (1)). -/
noncomputable def cW (θ : ℝ) : ℝ := 2 - θ / 2 - cot F (θ / F.sqrt 2) / F.sqrt 2

/-- h as a function of (c, k): (3 + k − √((1 − k)(9 − k − 8c)))/4. -/
noncomputable def hc (c k : ℝ) : ℝ := (3 + k - F.sqrt ((1 - k) * (9 - k - 8 * c))) / 4

/-- h(θ; k) = (3 + k − √((1 − k)(9 − k − 8c(θ))))/4 (E1.md §1; SL (2)). -/
noncomputable def h (θ k : ℝ) : ℝ := hc F (cW F θ) k

/-- Positivity threshold k*(θ) = 1 − 2/(2 − c(θ)) (E1.md §4). -/
noncomputable def kstar (θ : ℝ) : ℝ := 1 - 2 / (2 - cW F θ)

/-- Levinson constant κ = 1 − R⁻¹ log c (SL Thm 1.3, with c = c(P,Q,R,ν)). -/
noncomputable def kappa (c R : ℝ) : ℝ := 1 - R⁻¹ * F.log c

/-- liminf_{T→∞} f T ≥ a, elementary form: ∀ ε > 0, eventually a − ε ≤ f T. -/
def LiminfGe (f : ℝ → ℝ) (a : ℝ) : Prop := ∀ ε > 0, ∃ T₀ : ℝ, ∀ T ≥ T₀, a - ε ≤ f T

/-- The ratio A(T, T^θ)/N(T, T^θ). -/
noncomputable def ratio (A N : ℝ → ℝ → ℝ) (θ : ℝ) : ℝ → ℝ := fun T => A T (F.rpow T θ) / N T (F.rpow T θ)

/-- Real polynomial functions: the closure of constants and x under + and × (= image of Polynomial.eval). -/
inductive IsPoly : (ℝ → ℝ) → Prop
  | const (a : ℝ) : IsPoly (fun _ => a)
  | id : IsPoly (fun x => x)
  | add {f g : ℝ → ℝ} : IsPoly f → IsPoly g → IsPoly (fun x => f x + g x)
  | mul {f g : ℝ → ℝ} : IsPoly f → IsPoly g → IsPoly (fun x => f x * g x)

namespace IsPoly

theorem sub {f g : ℝ → ℝ} (hf : IsPoly f) (hg : IsPoly g) : IsPoly (fun x => f x - g x) := by
  have h := hf.add ((IsPoly.const (-1)).mul hg)
  simp only [neg_one_mul, ← sub_eq_add_neg] at h
  exact h

theorem smul {f : ℝ → ℝ} (a : ℝ) (hf : IsPoly f) : IsPoly (fun x => a * f x) := (IsPoly.const a).mul hf

theorem pow {f : ℝ → ℝ} (hf : IsPoly f) : ∀ n : ℕ, IsPoly (fun x => f x ^ n)
  | 0 => by simpa using IsPoly.const 1
  | n + 1 => by simpa [pow_succ] using (pow hf n).mul hf

theorem comp {f g : ℝ → ℝ} (hf : IsPoly f) (hg : IsPoly g) : IsPoly (fun x => f (g x)) := by
  induction hf with
  | const a => exact IsPoly.const a
  | id => exact hg
  | add _ _ ih1 ih2 => exact ih1.add ih2
  | mul _ _ ih1 ih2 => exact ih1.mul ih2

end IsPoly

end E1Lean
