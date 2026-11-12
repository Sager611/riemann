module

public import Riemann.Deep.Combos.Lean.Lemma1
public import Riemann.Deep.Combos.Lean.Ortho
public import Mathlib.Algebra.Order.BigOperators.Ring.Finset

@[expose] public section

/-!
# Combos NOTE §1, Lemma 2: the interface and the elementary inequalities

Abstract interface (NOTE §0, l.32–45). Mathlib has no ring class groups, so one negative
discriminant D is modelled by an RCData G ι:
* G : a finite abelian group, written additively: G = Cl(O) ≅ I_K(f)/P_{K,Z}(f) (l.32–33);
* ι : the O_K-ideals 𝔞 prime to f (the c_𝔞 of Lemma 2 are supported there, l.79);
* nrm 𝔞 = N𝔞, cls 𝔞 = [𝔞] ∈ G, cont 𝔞 = m(𝔞), prim 𝔞 = 𝔞₀ with 𝔞 = m(𝔞)𝔞₀ (l.37–39);
* characters: Mathlib's AddChar G ℂ, χ(𝔞) := χ([𝔞]) (l.34).
-/

open Finset ComplexConjugate

namespace RingClass

/-- NOTE §0 data for one negative discriminant D (l.32–39); see the module docstring. -/
structure RCData (G ι : Type*) [AddCommGroup G] where
  /-- the discriminant D = f² d_K -/
  D : ℤ
  D_neg : D < 0
  /-- N𝔞 -/
  nrm : ι → ℕ
  /-- [𝔞] ∈ I_K(f)/P_{K,Z}(f) ≅ Cl(O) = G (l.32–33, (F1)) -/
  cls : ι → G
  /-- m(𝔞): the largest m ≥ 1 with 𝔞 ⊂ mO_K (l.37) -/
  cont : ι → ℕ
  /-- 𝔞₀ with 𝔞 = m(𝔞)𝔞₀ (l.37) -/
  prim : ι → ι
  cont_pos : ∀ A, 0 < cont A
  /-- 𝔞₀ is primitive -/
  cont_prim : ∀ A, cont (prim A) = 1
  /-- χ(𝔞) = χ(𝔞₀), because mO_K ∈ P_{K,Z}(f) (l.38–39) -/
  cls_prim : ∀ A, cls (prim A) = cls A
  /-- N(m𝔞₀) = m² N𝔞₀ -/
  nrm_eq : ∀ A, nrm A = cont A ^ 2 * nrm (prim A)
  /-- 𝔞 = m(𝔞)𝔞₀ determines 𝔞 from (m(𝔞), 𝔞₀) -/
  decomp_inj : ∀ A B, prim A = prim B → cont A = cont B → A = B

variable {G ι : Type*} [AddCommGroup G]

/-- (F1)+(F2) as used in Lemma 2 (NOTE l.44–51, l.86–87; Cox Prop 7.20, 7.22, Thm 7.7): for
every class C there is a reduced form Q = (a, b, c) of discriminant D (with [𝔞_Q] = C⁻¹) such that
every primitive 𝔟 prime to f with [𝔟] = C is ideal(x, y) = α𝔞_Q⁻¹ (moved to O_K by (F1)) for some
α = xa + y(−b+√D)/2 with gcd(x, y) = 1 (primitivity: 𝔟 ⊂ mO ⟺ m | gcd(x, y)) and
N𝔟 = Q(x, −y). Only 0 < a, |b| ≤ a ≤ c of "reduced" is recorded (Lemma 1 uses no more). -/
def F12 (d : RCData G ι) : Prop :=
  ∀ C : G, ∃ (a b c : ℤ) (ideal : ℤ × ℤ → ι), 0 < a ∧ |b| ≤ a ∧ a ≤ c ∧
    b ^ 2 - 4 * a * c = d.D ∧
    ∀ B : ι, d.cont B = 1 → d.cls B = C →
      ∃ x y : ℤ, Int.gcd x y = 1 ∧ Qf a (-b) c x y = d.nrm B ∧ ideal (x, y) = B

/-- weighted Cauchy–Schwarz: |Σ z_i|² ≤ (Σ w_i⁻¹)(Σ w_i |z_i|²) -/
lemma normSq_sum_le {κ : Type*} (F : Finset κ) (w : κ → ℝ) (hw : ∀ i ∈ F, 0 < w i) (z : κ → ℂ) :
    Complex.normSq (∑ i ∈ F, z i) ≤ (∑ i ∈ F, (w i)⁻¹) * ∑ i ∈ F, w i * Complex.normSq (z i) := by
  have hre : (∑ i ∈ F, (z i).re) ^ 2 ≤ (∑ i ∈ F, (w i)⁻¹) * ∑ i ∈ F, w i * (z i).re ^ 2 := by
    apply Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul
    · intro i hi; exact (inv_pos.mpr (hw i hi)).le
    · intro i hi; exact mul_nonneg (hw i hi).le (sq_nonneg _)
    · intro i hi
      exact le_of_eq (by rw [← mul_assoc, inv_mul_cancel₀ (hw i hi).ne', one_mul])
  have him : (∑ i ∈ F, (z i).im) ^ 2 ≤ (∑ i ∈ F, (w i)⁻¹) * ∑ i ∈ F, w i * (z i).im ^ 2 := by
    apply Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul
    · intro i hi; exact (inv_pos.mpr (hw i hi)).le
    · intro i hi; exact mul_nonneg (hw i hi).le (sq_nonneg _)
    · intro i hi
      exact le_of_eq (by rw [← mul_assoc, inv_mul_cancel₀ (hw i hi).ne', one_mul])
  have e : ∀ i, w i * Complex.normSq (z i) = w i * (z i).re ^ 2 + w i * (z i).im ^ 2 := by
    intro i; rw [Complex.normSq_apply]; ring
  rw [Complex.normSq_apply, Complex.re_sum, Complex.im_sum,
    Finset.sum_congr rfl (fun i _ => e i), Finset.sum_add_distrib, mul_add]
  linarith

/-- Cauchy–Schwarz: |Σ_{i∈F} z_i|² ≤ #F Σ |z_i|² -/
lemma normSq_sum_le_card {κ : Type*} (F : Finset κ) (z : κ → ℂ) :
    Complex.normSq (∑ i ∈ F, z i) ≤ (F.card : ℝ) * ∑ i ∈ F, Complex.normSq (z i) := by
  have h := normSq_sum_le F (fun _ => (1 : ℝ)) (fun _ _ => one_pos) z
  simp only [inv_one, Finset.sum_const, nsmul_eq_mul, mul_one, one_mul] at h
  exact h

/-- the constant of Lemma 2 (NOTE l.87): h π_C ≤ 30 L² (N + u²), from (F3) h ≤ 3u²L² and
Lemma 1 π_C ≤ 2 + 8N/u² + 4x/u, where u = |D|^{1/4}, x = N^{1/2}. -/
lemma sieve_arith (L h π N x u : ℝ) (hu : 0 < u) (hx2 : x ^ 2 = N)
    (hh : h ≤ 3 * u ^ 2 * L ^ 2) (hπ : π ≤ 2 + 8 * N / u ^ 2 + 4 * x / u) (hπ0 : 0 ≤ π) :
    h * π ≤ 30 * L ^ 2 * (N + u ^ 2) := by
  have hA : 0 ≤ 3 * u ^ 2 * L ^ 2 := by positivity
  have h1 : h * π ≤ (3 * u ^ 2 * L ^ 2) * (2 + 8 * N / u ^ 2 + 4 * x / u) :=
    mul_le_mul hh hπ hπ0 hA
  have h2 : (3 * u ^ 2 * L ^ 2) * (2 + 8 * N / u ^ 2 + 4 * x / u)
      = 3 * L ^ 2 * (2 * u ^ 2 + 8 * N) + 6 * L ^ 2 * (2 * x * u) := by
    field_simp; ring
  have h3 : 2 * x * u ≤ x ^ 2 + u ^ 2 := by
    have := sq_nonneg (x - u); linarith [this]
  have h4 : 0 ≤ L ^ 2 := sq_nonneg L
  have h5 : L ^ 2 * (2 * x * u) ≤ L ^ 2 * (x ^ 2 + u ^ 2) := mul_le_mul_of_nonneg_left h3 h4
  have h6 : L ^ 2 * x ^ 2 = L ^ 2 * N := by rw [hx2]
  have h7 : 0 ≤ L ^ 2 * u ^ 2 := mul_nonneg h4 (sq_nonneg u)
  nlinarith [h1, h2, h5, h6, h7]

end RingClass

