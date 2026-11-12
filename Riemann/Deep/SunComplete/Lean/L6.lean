import Riemann.Deep.SunComplete.Lean.L6NT
import Riemann.Deep.SunComplete.Lean.L6Arith

/-!
# sunlean · L6: Lemma 6 (case (A)) without a sieve, for every q ≥ 196·16²·4^15 and every q with
`ω(q − 1) ≥ 16`; the rest is the named hypothesis `Lemma6Amid`

Choice: core = every prime divisor of `q − 1` (`δ = Δ = 1`, `W = 2^r`, `τ = θ ≥ 1/(r + 1)`),
`s = ⌊√q⌋ + 1`, `A' = 1/τ + E`, `D = φ(q − 1)`. (CS) holds once `q ≥ 196(r + 1)²4^r`
(`L6Arith`); `q − 1 ≥ (r + 1)!` gives that for `r ≥ 16`.
-/

namespace SunLean

theorem succ_le_two_pow : ∀ r : ℕ, r + 1 ≤ 2 ^ r
  | 0 => le_refl _
  | r + 1 => by
    rw [pow_succ]
    have := succ_le_two_pow r
    omega

/-- No sieve: the core is every prime divisor of `q − 1`. -/
def noSieve (q : ℕ) : SieveData q :=
  ⟨primeDivisors (q - 1), ∅, Finset.union_empty _, Finset.disjoint_empty_right _⟩

theorem noSieve_δ (q : ℕ) : (noSieve q).δ = 1 := by simp [SieveData.δ, noSieve]

theorem noSieve_τ (q : ℕ) :
    (noSieve q).τ = ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) := by
  unfold SieveData.τ
  rw [noSieve_δ, mul_one]
  rfl

theorem noSieve_E (q : ℕ) (s : ℚ) :
    (noSieve q).E s = (2 ^ (primeDivisors (q - 1)).card - 1) * s + 1 := by
  simp [SieveData.E, SieveData.W, SieveData.Δ, noSieve]

/-- (CS) without a sieve once `q ≥ 196(r + 1)²4^r`, `r = ω(q − 1)`. -/
theorem csa_noSieve (q : ℕ) (hq3 : 3 ≤ q)
    (hbig : 196 * ((primeDivisors (q - 1)).card + 1) ^ 2 * 4 ^ (primeDivisors (q - 1)).card ≤ q) :
    CSA q (noSieve q) ((Nat.sqrt q : ℚ) + 1) := by
  have hq1 : q - 1 ≠ 0 := by omega
  have hnQ : ((q - 1 : ℕ) : ℚ) = (q : ℚ) - 1 := by rw [Nat.cast_sub (by omega)]; simp
  have hθ : 1 ≤ (((primeDivisors (q - 1)).card : ℚ) + 1) *
      ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) :=
    theta_card_bound _ (fun p hp => two_le_of_mem_primeDivisors hp)
  have ha : 14 * ((primeDivisors (q - 1)).card + 1) * 2 ^ (primeDivisors (q - 1)).card ≤
      Nat.sqrt q := by
    rw [Nat.le_sqrt']
    have e : (14 * ((primeDivisors (q - 1)).card + 1) * 2 ^ (primeDivisors (q - 1)).card) ^ 2 =
        196 * ((primeDivisors (q - 1)).card + 1) ^ 2 * 4 ^ (primeDivisors (q - 1)).card := by
      rw [show (4 : ℕ) = 2 ^ 2 by norm_num, ← pow_mul, mul_comm 2 (primeDivisors (q - 1)).card,
        pow_mul]
      ring
    rw [e]; exact hbig
  have hK : 7 * ((primeDivisors (q - 1)).card + 1) * 2 ^ (primeDivisors (q - 1)).card *
      (Nat.sqrt q + 1) + 1 ≤ q := by
    have ht2 : Nat.sqrt q ^ 2 ≤ q := Nat.sqrt_le' q
    have hp1 : 1 ≤ 2 ^ (primeDivisors (q - 1)).card := Nat.one_le_two_pow
    have h1 : 14 * ((primeDivisors (q - 1)).card + 1) * 2 ^ (primeDivisors (q - 1)).card *
        Nat.sqrt q ≤ Nat.sqrt q * Nat.sqrt q := Nat.mul_le_mul_right _ ha
    have h2 : 2 ≤ Nat.sqrt q := by nlinarith
    have h3 : 7 * ((primeDivisors (q - 1)).card + 1) * 2 ^ (primeDivisors (q - 1)).card * 2 ≤
        7 * ((primeDivisors (q - 1)).card + 1) * 2 ^ (primeDivisors (q - 1)).card * Nat.sqrt q :=
      Nat.mul_le_mul_left _ h2
    have h4 : Nat.sqrt q * Nat.sqrt q = Nat.sqrt q ^ 2 := (sq _).symm
    nlinarith
  have hKQ : 7 * ((((primeDivisors (q - 1)).card : ℚ) + 1) *
      (2 ^ (primeDivisors (q - 1)).card * ((Nat.sqrt q : ℚ) + 1))) ≤ (q : ℚ) - 1 := by
    have h : ((7 * ((primeDivisors (q - 1)).card + 1) * 2 ^ (primeDivisors (q - 1)).card *
        (Nat.sqrt q + 1) + 1 : ℕ) : ℚ) ≤ (q : ℚ) := by exact_mod_cast hK
    push_cast at h
    linarith
  have hw1 : (1 : ℚ) ≤ 2 ^ (primeDivisors (q - 1)).card := one_le_pow₀ (by norm_num)
  have hwr : ((primeDivisors (q - 1)).card : ℚ) + 1 ≤ 2 ^ (primeDivisors (q - 1)).card := by
    exact_mod_cast succ_le_two_pow _
  have hε0 : 0 ≤ epsA q := by unfold epsA; split_ifs <;> norm_num
  have hε1 : epsA q ≤ 1 := by unfold epsA; split_ifs <;> norm_num
  have hq3Q : (3 : ℚ) ≤ q := by exact_mod_cast hq3
  have hsq0 : (0 : ℚ) ≤ Nat.sqrt q := Nat.cast_nonneg _
  obtain ⟨c1, c2, c3, c4⟩ := noSieve_arith (n := (q : ℚ) - 1) (q := (q : ℚ))
    (w := 2 ^ (primeDivisors (q - 1)).card) (r1 := ((primeDivisors (q - 1)).card : ℚ) + 1)
    (s := (Nat.sqrt q : ℚ) + 1) (θ := ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)))
    (ε := epsA q) (by ring) (by linarith) hw1 (le_add_of_nonneg_left (Nat.cast_nonneg _)) hwr
    (by linarith) hθ hKQ hε0 hε1
  have hphi : (phi (q - 1) : ℚ) = ((q : ℚ) - 1) * ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) := by
    rw [phi_eq_prod _ hq1, hnQ]
  unfold CSA
  refine ⟨by rw [noSieve_δ]; norm_num, by positivity, ?_,
    1 / (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))) +
      ((2 ^ (primeDivisors (q - 1)).card - 1) * ((Nat.sqrt q : ℚ) + 1) + 1),
    phi (q - 1), ?_, ?_, ?_, c3, ?_⟩
  · have h := Nat.lt_succ_sqrt' q
    have h' : (q : ℚ) < ((Nat.sqrt q + 1 : ℕ) : ℚ) ^ 2 := by exact_mod_cast h
    push_cast at h'
    linarith
  · unfold alphaStarA
    rw [noSieve_τ, noSieve_E]
    exact c1
  · unfold delta0A
    linarith
  · rw [hphi]; exact c2
  · rw [noSieve_E, hphi]; exact c4


/-- (CS) without a sieve once `7·2^r·(⌊√q⌋ + 1) ≤ φ(q − 1)`, `r = ω(q − 1)`. -/
theorem csa_noSieve_phi (q : ℕ) (hq3 : 3 ≤ q)
    (hphi : 7 * 2 ^ (primeDivisors (q - 1)).card * (Nat.sqrt q + 1) ≤ phi (q - 1)) :
    CSA q (noSieve q) ((Nat.sqrt q : ℚ) + 1) := by
  have hq1 : q - 1 ≠ 0 := by omega
  have hnQ : ((q - 1 : ℕ) : ℚ) = (q : ℚ) - 1 := by rw [Nat.cast_sub (by omega)]; simp
  have hθr : 1 ≤ (((primeDivisors (q - 1)).card : ℚ) + 1) *
      ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) :=
    theta_card_bound _ (fun p hp => two_le_of_mem_primeDivisors hp)
  have hfac : ∀ p ∈ primeDivisors (q - 1), 0 ≤ 1 - 1 / (p : ℚ) ∧ 1 - 1 / (p : ℚ) ≤ 1 := by
    intro p hp
    have h2 : (2 : ℚ) ≤ p := by exact_mod_cast two_le_of_mem_primeDivisors hp
    have h3 : 1 / (p : ℚ) ≤ 1 := by rw [div_le_one (by linarith)]; linarith
    have h4 : 0 ≤ 1 / (p : ℚ) := by positivity
    constructor <;> linarith
  have hθ1 : ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) ≤ 1 :=
    Finset.prod_le_one₀ (fun p hp => (hfac p hp).1) (fun p hp => (hfac p hp).2)
  have hθ0 : 0 < ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) := by
    rcases lt_or_ge 0 (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))) with h | h
    · exact h
    · exfalso
      have := mul_le_mul_of_nonneg_left h
        (show (0 : ℚ) ≤ ((primeDivisors (q - 1)).card : ℚ) + 1 by positivity)
      rw [mul_zero] at this
      linarith
  have hphiQ : (phi (q - 1) : ℚ) = ((q : ℚ) - 1) * ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) := by
    rw [phi_eq_prod _ hq1, hnQ]
  have hw1 : (1 : ℚ) ≤ 2 ^ (primeDivisors (q - 1)).card := one_le_pow₀ (by norm_num)
  have hwr0 : ((primeDivisors (q - 1)).card : ℚ) + 1 ≤ 2 ^ (primeDivisors (q - 1)).card := by
    exact_mod_cast succ_le_two_pow _
  have hinv : 1 / (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))) ≤
      ((primeDivisors (q - 1)).card : ℚ) + 1 := by
    rw [div_le_iff₀ hθ0]; linarith
  have hr1 : 1 ≤ 1 / (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))) := by
    rw [le_div_iff₀ hθ0]; linarith
  have hθ : 1 ≤ 1 / (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))) *
      ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) := by
    rw [one_div_mul_cancel hθ0.ne']
  have hKQ : 7 * (1 / (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))) *
      (2 ^ (primeDivisors (q - 1)).card * ((Nat.sqrt q : ℚ) + 1))) ≤ (q : ℚ) - 1 := by
    have h : ((7 * 2 ^ (primeDivisors (q - 1)).card * (Nat.sqrt q + 1) : ℕ) : ℚ) ≤
        (phi (q - 1) : ℚ) := by exact_mod_cast hphi
    rw [hphiQ] at h
    push_cast at h
    rw [show 7 * (1 / (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))) *
        (2 ^ (primeDivisors (q - 1)).card * ((Nat.sqrt q : ℚ) + 1))) =
        7 * 2 ^ (primeDivisors (q - 1)).card * ((Nat.sqrt q : ℚ) + 1) /
          ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)) by ring, div_le_iff₀ hθ0]
    exact h
  have hε0 : 0 ≤ epsA q := by unfold epsA; split_ifs <;> norm_num
  have hε1 : epsA q ≤ 1 := by unfold epsA; split_ifs <;> norm_num
  have hq3Q : (3 : ℚ) ≤ q := by exact_mod_cast hq3
  have hsq0 : (0 : ℚ) ≤ Nat.sqrt q := Nat.cast_nonneg _
  obtain ⟨c1, c2, c3, c4⟩ := noSieve_arith (n := (q : ℚ) - 1) (q := (q : ℚ))
    (w := 2 ^ (primeDivisors (q - 1)).card)
    (r1 := 1 / (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))))
    (s := (Nat.sqrt q : ℚ) + 1) (θ := ∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ)))
    (ε := epsA q) (by ring) (by linarith) hw1 hr1 (le_trans hinv hwr0) (by linarith) hθ hKQ hε0 hε1
  unfold CSA
  refine ⟨by rw [noSieve_δ]; norm_num, by positivity, ?_,
    1 / (∏ p ∈ primeDivisors (q - 1), (1 - 1 / (p : ℚ))) +
      ((2 ^ (primeDivisors (q - 1)).card - 1) * ((Nat.sqrt q : ℚ) + 1) + 1),
    phi (q - 1), ?_, ?_, ?_, c3, ?_⟩
  · have h := Nat.lt_succ_sqrt' q
    have h' : (q : ℚ) < ((Nat.sqrt q + 1 : ℕ) : ℚ) ^ 2 := by exact_mod_cast h
    push_cast at h'
    linarith
  · unfold alphaStarA
    rw [noSieve_τ, noSieve_E]
    exact c1
  · unfold delta0A
    linarith
  · rw [hphiQ]; exact c2
  · rw [noSieve_E, hphiQ]; exact c4

/-- `196(r + 1)²4^r ≤ (r + 1)!` for `r ≥ 16`. -/
theorem bound_fact (r : ℕ) (hr : 16 ≤ r) : 196 * (r + 1) ^ 2 * 4 ^ r ≤ (r + 1).factorial := by
  induction r, hr using Nat.le_induction with
  | base => decide
  | succ r hr ih =>
    have h1 : 4 * (r + 2) ≤ (r + 1) ^ 2 := by nlinarith
    calc 196 * (r + 1 + 1) ^ 2 * 4 ^ (r + 1) = (196 * 4 ^ r * (r + 2)) * (4 * (r + 2)) := by ring
      _ ≤ (196 * 4 ^ r * (r + 2)) * ((r + 1) ^ 2) := Nat.mul_le_mul_left _ h1
      _ = (r + 2) * (196 * (r + 1) ^ 2 * 4 ^ r) := by ring
      _ ≤ (r + 2) * (r + 1).factorial := Nat.mul_le_mul_left _ ih
      _ = (r + 1 + 1).factorial := (Nat.factorial_succ (r + 1)).symm

/-- **Lemma 6 (case (A)), proved part**: every `q ≥ 196·16²·4^15 = 53,876,069,761,024` and every
`q ≥ 3` with `ω(q − 1) ≥ 16` has sieve data with (CS) (no sieve needed). -/
theorem lemma6A_big (q : ℕ) (hq : 3 ≤ q)
    (h : 53876069761024 ≤ q ∨ 16 ≤ (primeDivisors (q - 1)).card) :
    ∃ (d : SieveData q) (s : ℚ), CSA q d s := by
  refine ⟨noSieve q, _, csa_noSieve q hq ?_⟩
  by_cases hr : 16 ≤ (primeDivisors (q - 1)).card
  · have h1 := bound_fact _ hr
    have h2 := fact_le_prod (primeDivisors (q - 1)) (fun p hp => two_le_of_mem_primeDivisors hp)
    have h3 := prod_primeDivisors_le (q - 1) (by omega)
    exact le_trans h1 (le_trans h2 (le_trans h3 (Nat.sub_le q 1)))
  · rcases h with h | h
    · push Not at hr
      have h1 : ((primeDivisors (q - 1)).card + 1) ^ 2 ≤ 16 ^ 2 := Nat.pow_le_pow_left (by omega) 2
      have h2 : 4 ^ (primeDivisors (q - 1)).card ≤ 4 ^ 15 := Nat.pow_le_pow_right (by norm_num) (by omega)
      calc 196 * ((primeDivisors (q - 1)).card + 1) ^ 2 * 4 ^ (primeDivisors (q - 1)).card
          ≤ 196 * 16 ^ 2 * 4 ^ 15 := Nat.mul_le_mul (Nat.mul_le_mul_left _ h1) h2
        _ = 53876069761024 := by norm_num
        _ ≤ q := h
    · exact absurd h hr

/-- **What is left of Lemma 6**, NAMED HYPOTHESIS: (CS) for the prime powers
`30140 ≤ q < 53,876,069,761,024` with `ω(q − 1) ≤ 15` and `φ(q − 1) < 7·2^{ω(q−1)}·(⌊√q⌋ + 1)`
(NOTE Lemma 6 with the sieve). -/
def Lemma6Amid : Prop :=
  ∀ q : ℕ, PrimePow q → 30140 ≤ q → q < 53876069761024 → (primeDivisors (q - 1)).card ≤ 15 →
    phi (q - 1) < 7 * 2 ^ (primeDivisors (q - 1)).card * (Nat.sqrt q + 1) →
    ∃ (d : SieveData q) (s : ℚ), CSA q d s

/-- Lemma 6 from its middle range. -/
theorem lemma6A_of_mid (h : Lemma6Amid) : Lemma6A := by
  intro q hq h1
  by_cases hbig : 53876069761024 ≤ q
  · exact lemma6A_big q (by omega) (Or.inl hbig)
  · by_cases hr : 16 ≤ (primeDivisors (q - 1)).card
    · exact lemma6A_big q (by omega) (Or.inr hr)
    · by_cases hφ : 7 * 2 ^ (primeDivisors (q - 1)).card * (Nat.sqrt q + 1) ≤ phi (q - 1)
      · exact ⟨noSieve q, _, csa_noSieve_phi q (by omega) hφ⟩
      · exact h q hq h1 (by omega) (by omega) (by omega)

end SunLean

#print axioms SunLean.lemma6A_big
#print axioms SunLean.csa_noSieve_phi
#print axioms SunLean.lemma6A_of_mid
