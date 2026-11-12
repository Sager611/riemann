import Riemann.Deep.EvoAlg.Lean.CorA2
import Mathlib.Data.Nat.Prime.Defs

/-!
# Corollary A (NOTE §0, first item): cube roots of all primes suffice
-/

namespace EvoAlg

variable {𝕜 : Type*} [Field 𝕜]

/-- If every prime is a cube in `𝕜`, so is every positive integer. -/
theorem cube_of_primes (hp : ∀ p : ℕ, p.Prime → ∃ z : 𝕜, z ^ 3 = p) :
    ∀ n : ℕ, 0 < n → ∃ z : 𝕜, z ^ 3 = n := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro hn
    by_cases h1 : n = 1
    · exact ⟨1, by rw [h1]; simp⟩
    · obtain ⟨p, hpp, m, rfl⟩ := Nat.exists_prime_and_dvd h1
      have hm : 0 < m := Nat.pos_of_ne_zero (fun h => by simp [h] at hn)
      obtain ⟨z, hz⟩ := hp p hpp
      obtain ⟨w, hw⟩ := ih m (lt_mul_of_one_lt_left hm hpp.one_lt) hm
      exact ⟨z * w, by rw [mul_pow, hz, hw, Nat.cast_mul]⟩

/-- **Corollary A** (first item): over every field of characteristic 0 that contains a cube root
of every prime, the survey's statement holds for all finite connected graphs. -/
theorem corA_primes [CharZero 𝕜] (hp : ∀ p : ℕ, p.Prime → ∃ z : 𝕜, z ^ 3 = p) :
    SurveyHolds 𝕜 := by
  refine corA_iff.mpr (fun a b ha hb => ?_)
  obtain ⟨za, hza⟩ := cube_of_primes hp a ha
  obtain ⟨zb, hzb⟩ := cube_of_primes hp b hb
  exact ⟨za / zb, by rw [div_pow, hza, hzb]⟩

end EvoAlg

