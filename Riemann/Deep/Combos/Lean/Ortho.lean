import Mathlib.Algebra.Group.AddChar
import Mathlib.Basic.Complex.BigOperators
import Mathlib.GroupTheory.OrderOfElement
import Mathlib.Tactic.Linarith

/-!
# Character orthogonality in Bessel form (for Lemma 2, NOTE l.82–83)

For a finite abelian group G (written additively) and any finite set 𝒞 of characters
χ ∈ AddChar G ℂ: Σ_{χ∈𝒞} |Σ_g χ(g) v(g)|² ≤ |G| Σ_g |v(g)|². For 𝒞 = Ĝ the NOTE uses the equality
Σ_χ |Σ_g χ(g)v(g)|² = |G|Σ_g|v(g)|²; Lemma 2 needs only "≤". Proof: first orthogonality
Σ_g χ'(g) conj χ(g) = |G|·[χ = χ'] and 0 ≤ Σ_g |n v(g) − Σ_χ c_χ conj χ(g)|² (Bessel), n = |G|.
|z|² is Complex.normSq z.
-/

open Finset ComplexConjugate

namespace RingClass

variable {G : Type*} [AddCommGroup G] [Fintype G]

lemma char_mul_conj (χ : AddChar G ℂ) (g : G) : χ g * conj (χ g) = 1 := by
  have h1 : χ g ^ Fintype.card G = 1 := by
    rw [← AddChar.map_nsmul_eq_pow, card_nsmul_eq_zero, AddChar.map_zero_eq_one]
  have h2 : Complex.normSq (χ g) ^ Fintype.card G = 1 := by
    rw [← map_pow, h1, map_one]
  have h3 : Complex.normSq (χ g) = 1 :=
    (pow_eq_one_iff_of_nonneg (Complex.normSq_nonneg _) Fintype.card_ne_zero).mp h2
  rw [Complex.mul_conj, h3, Complex.ofReal_one]

/-- g ↦ χ'(g) conj χ(g), again a character -/
def twist (χ χ' : AddChar G ℂ) : AddChar G ℂ where
  toFun g := χ' g * conj (χ g)
  map_zero_eq_one' := by simp
  map_add_eq_mul' a b := by
    simp only [AddChar.map_add_eq_mul, map_mul]; ring

/-- first orthogonality: Σ_g χ'(g) conj χ(g) = |G| [χ = χ'] -/
lemma ortho (χ χ' : AddChar G ℂ) [Decidable (χ = χ')] :
    ∑ g, χ' g * conj (χ g) = if χ = χ' then (Fintype.card G : ℂ) else 0 := by
  classical
  have e : twist χ χ' = 0 ↔ χ = χ' := by
    constructor
    · intro hq
      ext g
      have h1 : χ' g * conj (χ g) = 1 := by
        have := congrArg (fun ψ : AddChar G ℂ => ψ g) hq
        exact this
      have h2 := char_mul_conj χ g
      calc χ g = χ g * (χ' g * conj (χ g)) := by rw [h1, mul_one]
        _ = χ' g * (χ g * conj (χ g)) := by ring
        _ = χ' g := by rw [h2, mul_one]
    · rintro rfl
      ext g
      exact char_mul_conj χ g
  have hs : ∑ g, χ' g * conj (χ g) = ∑ g, twist χ χ' g := rfl
  rw [hs, AddChar.sum_eq_ite]
  split_ifs with h1 h2 h2
  · rfl
  · exact absurd (e.mp h1) h2
  · exact absurd (e.mpr h2) h1
  · rfl

/-- c_χ = Σ_g χ(g) v(g) -/
def coef (v : G → ℂ) (χ : AddChar G ℂ) : ℂ := ∑ g, χ g * v g

/-- Σ_g v(g) conj(u(g)) = Σ_χ |c_χ|², u(g) = Σ_χ c_χ conj χ(g) -/
lemma bessel_A (𝒞 : Finset (AddChar G ℂ)) (v : G → ℂ) :
    ∑ g, v g * conj (∑ χ ∈ 𝒞, coef v χ * conj (χ g))
      = ((∑ χ ∈ 𝒞, Complex.normSq (coef v χ) : ℝ) : ℂ) := by
  have h1 : ∀ g, v g * conj (∑ χ ∈ 𝒞, coef v χ * conj (χ g))
      = ∑ χ ∈ 𝒞, conj (coef v χ) * (χ g * v g) := by
    intro g
    simp only [map_sum, map_mul, Complex.conj_conj, Finset.mul_sum]
    exact Finset.sum_congr rfl (fun χ _ => by ring)
  rw [Finset.sum_congr rfl (fun g _ => h1 g), Finset.sum_comm, Complex.ofReal_sum]
  refine Finset.sum_congr rfl (fun χ _ => ?_)
  rw [← Finset.mul_sum, mul_comm]
  exact Complex.mul_conj (coef v χ)

/-- Σ_g u(g) conj(u(g)) = |G| Σ_χ |c_χ|² -/
lemma bessel_B (𝒞 : Finset (AddChar G ℂ)) (v : G → ℂ) :
    ∑ g, (∑ χ ∈ 𝒞, coef v χ * conj (χ g)) * conj (∑ χ ∈ 𝒞, coef v χ * conj (χ g))
      = (((Fintype.card G : ℝ) * ∑ χ ∈ 𝒞, Complex.normSq (coef v χ) : ℝ) : ℂ) := by
  classical
  have h1 : ∀ g, (∑ χ ∈ 𝒞, coef v χ * conj (χ g)) * conj (∑ χ ∈ 𝒞, coef v χ * conj (χ g))
      = ∑ χ ∈ 𝒞, ∑ χ' ∈ 𝒞, (coef v χ * conj (coef v χ')) * (χ' g * conj (χ g)) := by
    intro g
    simp only [map_sum, map_mul, Complex.conj_conj, Finset.sum_mul_sum]
    exact Finset.sum_congr rfl (fun χ _ => Finset.sum_congr rfl (fun χ' _ => by ring))
  rw [Finset.sum_congr rfl (fun g _ => h1 g), Finset.sum_comm]
  have h2 : ∀ χ ∈ 𝒞, ∑ g, ∑ χ' ∈ 𝒞, (coef v χ * conj (coef v χ')) * (χ' g * conj (χ g))
      = (Fintype.card G : ℂ) * (Complex.normSq (coef v χ) : ℂ) := by
    intro χ hχ
    rw [Finset.sum_comm]
    have h3 : ∀ χ' ∈ 𝒞, ∑ g, (coef v χ * conj (coef v χ')) * (χ' g * conj (χ g))
        = (coef v χ * conj (coef v χ')) * (if χ = χ' then (Fintype.card G : ℂ) else 0) := by
      intro χ' _
      rw [← Finset.mul_sum, ortho χ χ']
    rw [Finset.sum_congr rfl h3]
    simp only [mul_ite, mul_zero]
    simp only [Finset.sum_ite_eq, hχ, ↓reduceIte]
    rw [Complex.mul_conj]
    ring
  rw [Finset.sum_congr rfl h2, ← Finset.mul_sum]
  push_cast
  ring

/-- **Bessel inequality for characters**: Σ_{χ∈𝒞} |Σ_g χ(g)v(g)|² ≤ |G| Σ_g |v(g)|². -/
theorem bessel (𝒞 : Finset (AddChar G ℂ)) (v : G → ℂ) :
    ∑ χ ∈ 𝒞, Complex.normSq (∑ g, χ g * v g) ≤ (Fintype.card G : ℝ) * ∑ g, Complex.normSq (v g) := by
  have hn0 : (0 : ℝ) < Fintype.card G := by exact_mod_cast Fintype.card_pos
  have hA := bessel_A 𝒞 v
  have hB := bessel_B 𝒞 v
  set n : ℝ := (Fintype.card G : ℝ) with hn
  set S : ℝ := ∑ χ ∈ 𝒞, Complex.normSq (coef v χ) with hS
  have hU : ∑ g, Complex.normSq (∑ χ ∈ 𝒞, coef v χ * conj (χ g)) = n * S := by
    have h := congrArg Complex.re hB
    rw [Complex.re_sum, Complex.ofReal_re] at h
    rw [← h]
    refine Finset.sum_congr rfl (fun g _ => ?_)
    rw [Complex.mul_conj, Complex.ofReal_re]
  have hV : ∑ g, (v g * conj (∑ χ ∈ 𝒞, coef v χ * conj (χ g))).re = S := by
    have h := congrArg Complex.re hA
    rwa [Complex.re_sum, Complex.ofReal_re] at h
  have hexp : ∀ g, Complex.normSq ((n : ℂ) * v g - ∑ χ ∈ 𝒞, coef v χ * conj (χ g))
      = n * n * Complex.normSq (v g) + Complex.normSq (∑ χ ∈ 𝒞, coef v χ * conj (χ g))
        - 2 * (n * (v g * conj (∑ χ ∈ 𝒞, coef v χ * conj (χ g))).re) := by
    intro g
    rw [Complex.normSq_sub, Complex.normSq_mul, Complex.normSq_ofReal]
    simp only [mul_assoc, Complex.re_ofReal_mul]
  have hpos : 0 ≤ ∑ g, Complex.normSq ((n : ℂ) * v g - ∑ χ ∈ 𝒞, coef v χ * conj (χ g)) :=
    Finset.sum_nonneg (fun g _ => Complex.normSq_nonneg _)
  rw [Finset.sum_congr rfl (fun g _ => hexp g), Finset.sum_sub_distrib, Finset.sum_add_distrib,
    hU, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, hV] at hpos
  have h3 : n * S ≤ n * (n * ∑ g, Complex.normSq (v g)) := by linarith
  exact le_of_mul_le_mul_left h3 hn0

end RingClass

