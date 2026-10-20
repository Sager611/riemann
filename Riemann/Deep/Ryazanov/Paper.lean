import Mathlib

/-!
# Ryazanov, arXiv 1808.10774v11 (31 Aug 2026), laid down and pinned

Paper: Báez-Duarte criterion (Theorem D): RH ⟸ `χ = 1_{[1,∞)}` is an `L²((0,∞), t⁻² dt)`-limit of
functions `∑ h_k {t/k}` (class `𝔑`). Sequence (4.1):
`φ*_n = -∑_{k<n!} μ(k) a_k + δ(n) n! a_{n!}`, `a_k(t) = {t/k}`, `δ(n) = ∑_{k<n!} μ(k)/k`.
Proof of Theorem 1: `φ*_n ≡ 1` on `[1,n!)` (correct, by `∑ μ(k)[t/k] = 1`); then
"`φ*_n` are `n!`-periodic because each `a_k`, `k = 1,…,n!`, is `n!`-periodic. Thus `φ*_n ≡ 0` on each
`[n!k, n!k+1)` and coincide with `χ` outside of all these intervals", whence (4.3)
`‖χ - φ*_n‖²_* = ∑ 1/(n!k(n!k+1)) < ζ(2)/n!² → 0`.

Failing step: the periodicity sentence. `a_k` is `n!`-periodic only when `k ∣ n!`.
Proved here for EVERY `n ≥ 3` (Bertrand prime `n!/2 < p < n!`, jump of `φ*_n` at `2p`):
* `akPeriodic_false`   : `a_{n!-1}` is not `n!`-periodic;
* `not_periodic`       : `φ*_n` is not `n!`-periodic;
* `coincide_false`     : `φ*_n ≠ χ` somewhere in `(n!, 2n!)`, outside all the exceptional intervals.
So (4.3) has no basis. Theorem D (true, Báez-Duarte 2003) and the chain are kept as Props.
-/

namespace Riemann.Deep.Ryazanov

open MeasureTheory Filter Topology Finset

/-- `n!`. -/
abbrev N (n : ℕ) : ℕ := Nat.factorial n

/-- Möbius function, in `ℝ`. -/
noncomputable def mu (k : ℕ) : ℝ := (ArithmeticFunction.moebius k : ℝ)

/-- `a_k(t) = {t/k}`. -/
noncomputable def a (k : ℕ) (t : ℝ) : ℝ := Int.fract (t / k)

/-- (4.2) `δ(n) = ∑_{k=1}^{n!-1} μ(k)/k`. -/
noncomputable def delta (n : ℕ) : ℝ := ∑ k ∈ Ico 1 (N n), mu k / k

/-- (4.1) `φ*_n = -∑_{k=1}^{n!-1} μ(k) a_k + δ(n) n! a_{n!}`. -/
noncomputable def phiStar (n : ℕ) (t : ℝ) : ℝ :=
  -(∑ k ∈ Ico 1 (N n), mu k * a k t) + delta n * N n * a (N n) t

/-- `χ = 1_{[1,∞)}`. -/
noncomputable def chi (t : ℝ) : ℝ := if 1 ≤ t then 1 else 0

/-- Class `𝔑`: finite combinations `∑_{k=1}^m h_k a_k`. -/
def InR (φ : ℝ → ℝ) : Prop := ∃ (m : ℕ) (h : ℕ → ℝ), ∀ t, φ t = ∑ k ∈ Ico 1 (m + 1), h k * a k t

/-- Squared distance in `L²_* = L²((0,∞), t⁻² dt)`. -/
noncomputable def distSq (f g : ℝ → ℝ) : ENNReal :=
  ∫⁻ t in Set.Ioi 0, ENNReal.ofReal ((f t - g t) ^ 2 / t ^ 2)

/-- Theorem D (Báez-Duarte 2003, after `t ↦ 1/t`), sufficiency direction. True; cited. -/
def TheoremD : Prop :=
  ∀ φ : ℕ → ℝ → ℝ, (∀ n, InR (φ n)) →
    Tendsto (fun n => distSq chi (φ n)) atTop (𝓝 0) → RiemannHypothesis

/-- The paper's claim: `φ*_n → χ` in `L²_*`. -/
def Converges : Prop := Tendsto (fun n => distSq chi (phiStar n)) atTop (𝓝 0)

lemma one_le_N (n : ℕ) : 1 ≤ N n := Nat.one_le_iff_ne_zero.mpr (Nat.factorial_ne_zero n)

theorem phiStar_mem (n : ℕ) : InR (phiStar n) := by
  refine ⟨N n, fun k => if k < N n then -mu k else delta n * N n, fun t => ?_⟩
  rw [Finset.sum_Ico_succ_top (one_le_N n), phiStar, ← Finset.sum_neg_distrib]
  simp only [lt_irrefl, ite_false]
  congr 1
  refine Finset.sum_congr rfl fun k hk => ?_
  simp only [(Finset.mem_Ico.mp hk).2, ite_true, neg_mul]

/-- The paper's logic, machine-checked: Theorem D + (4.3)-convergence ⇒ RH. -/
theorem ryazanov_RH (hD : TheoremD) (hc : Converges) : RiemannHypothesis :=
  hD phiStar phiStar_mem hc

/-! ## The paper's justification of `Converges`, as Props -/

/-- "each function `a_k`, `k = 1,…,n!`, is `n!`-periodic". -/
def AkPeriodicClaim (n : ℕ) : Prop := ∀ k ∈ Ico 1 (N n + 1), Function.Periodic (a k) (N n)

/-- "functions `φ*_n` are `n!`-periodic". -/
def PeriodicClaim (n : ℕ) : Prop := Function.Periodic (phiStar n) (N n)

/-- "`φ*_n` … coincide with `χ` outside of all `[n!k, n!k+1)`" (the input of (4.3)). -/
def CoincideClaim (n : ℕ) : Prop :=
  ∀ t : ℝ, 1 ≤ t → (∀ k : ℕ, ¬ ((N n * k : ℝ) ≤ t ∧ t < N n * k + 1)) → phiStar n t = 1

/-- The paper's inference "`a_k` periodic ⇒ `φ*_n` periodic" is valid; its premise is not. -/
theorem periodic_of_ak (n : ℕ) (h : AkPeriodicClaim n) : PeriodicClaim n := by
  intro t
  have hN : N n ∈ Ico 1 (N n + 1) := Finset.mem_Ico.mpr ⟨one_le_N n, Nat.lt_succ_self _⟩
  simp only [phiStar]
  rw [h _ hN t]
  congr 2
  refine Finset.sum_congr rfl fun k hk => ?_
  have hk' : k ∈ Ico 1 (N n + 1) := by
    rcases Finset.mem_Ico.mp hk with ⟨h1, h2⟩; exact Finset.mem_Ico.mpr ⟨h1, by omega⟩
  rw [h k hk' t]

/-! ## Arithmetic of the jumps -/

lemma six_le_N {n : ℕ} (hn : 3 ≤ n) : 6 ≤ N n := Nat.factorial_le hn

lemma a_nat {k : ℕ} (hk : 0 < k) (m : ℕ) : a k m = m / k - ((m / k : ℕ) : ℝ) := by
  have hk' : (k : ℝ) ≠ 0 := by exact_mod_cast hk.ne'
  have hm : (m : ℝ) = (m % k : ℕ) + k * (m / k : ℕ) := by exact_mod_cast (Nat.mod_add_div m k).symm
  rw [a, Int.fract_div_natCast_eq_div_natCast_mod]
  field_simp
  linarith

/-- On integers: `φ*_n(m) = ∑_{k<n!} μ(k)⌊m/k⌋ - δ(n) n! ⌊m/n!⌋` (the slopes cancel). -/
lemma phi_nat (n m : ℕ) :
    phiStar n m = ∑ k ∈ Ico 1 (N n), mu k * ((m / k : ℕ) : ℝ) - delta n * N n * ((m / N n : ℕ) : ℝ) := by
  have hN : (0 : ℝ) < N n := by exact_mod_cast one_le_N n
  have hs : ∑ k ∈ Ico 1 (N n), mu k * a k m =
      m * delta n - ∑ k ∈ Ico 1 (N n), mu k * ((m / k : ℕ) : ℝ) := by
    rw [delta, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun k hk => ?_
    rw [a_nat (Finset.mem_Ico.mp hk).1]
    ring
  rw [phiStar, hs, a_nat (one_le_N n)]
  field_simp
  ring

/-- Jump of `φ*_n` at the integer `m+1`. -/
lemma jump (n m : ℕ) :
    phiStar n ((m + 1 : ℕ) : ℝ) - phiStar n m =
      ∑ k ∈ (Ico 1 (N n)).filter (· ∣ m + 1), mu k - (if N n ∣ m + 1 then delta n * N n else 0) := by
  rw [phi_nat, phi_nat, Finset.sum_filter]
  simp only [Nat.succ_div, Nat.cast_add, Nat.cast_ite, Nat.cast_one, Nat.cast_zero, mul_add,
    Finset.sum_add_distrib, mul_ite, mul_one, mul_zero]
  split_ifs <;> ring

lemma sum_mu_divisors {j : ℕ} (hj : 1 < j) : ∑ d ∈ j.divisors, mu d = 0 := by
  have h := congrArg (fun f : ArithmeticFunction ℝ => f j)
    (ArithmeticFunction.coe_moebius_mul_coe_zeta (R := ℝ))
  simp only [ArithmeticFunction.coe_mul_zeta_apply, ArithmeticFunction.intCoe_apply,
    ArithmeticFunction.one_apply_ne (show j ≠ 1 by omega)] at h
  simpa [mu] using h

/-- For `2 ≤ m+1 < n!` the jump is `∑_{d∣m+1} μ(d) = 0`. -/
lemma jump_small (n m : ℕ) (h2 : 1 ≤ m) (hlt : m + 1 < N n) :
    phiStar n ((m + 1 : ℕ) : ℝ) - phiStar n m = 0 := by
  have hf : (Ico 1 (N n)).filter (· ∣ m + 1) = (m + 1).divisors := by
    ext d
    simp only [Finset.mem_filter, Finset.mem_Ico, Nat.mem_divisors]
    constructor
    · rintro ⟨-, hd⟩; exact ⟨hd, by omega⟩
    · rintro ⟨hd, -⟩
      have := Nat.le_of_dvd (by omega) hd
      exact ⟨⟨Nat.pos_of_dvd_of_pos hd (by omega), by omega⟩, hd⟩
  have hnd : ¬ N n ∣ m + 1 := fun h => absurd (Nat.le_of_dvd (by omega) h) (by omega)
  rw [jump, hf, sum_mu_divisors (by omega), ite_eq_right hnd, sub_zero]

/-- For a prime `p ≠ 2` with `n! < 2p` and `p < n!`, the jump at `2p` is `-μ(2p) = -1`. -/
lemma jump_two_p (n m p : ℕ) (hp : p.Prime) (hp2 : p ≠ 2) (hm : m + 1 = 2 * p)
    (h1 : N n < 2 * p) (h2 : p < N n) :
    phiStar n ((m + 1 : ℕ) : ℝ) - phiStar n m = -1 := by
  have hf : (Ico 1 (N n)).filter (· ∣ m + 1) = (m + 1).divisors.erase (m + 1) := by
    ext d
    simp only [Finset.mem_filter, Finset.mem_Ico, Finset.mem_erase, Nat.mem_divisors]
    constructor
    · rintro ⟨⟨-, hd⟩, hdv⟩; exact ⟨by omega, hdv, by omega⟩
    · rintro ⟨hne, hdv, -⟩
      obtain ⟨c, hc⟩ := hdv
      have hc2 : 2 ≤ c := by
        rcases Nat.lt_or_ge c 2 with h | h
        · interval_cases c <;> simp_all
        · exact h
      have : d * 2 ≤ d * c := Nat.mul_le_mul_left d hc2
      refine ⟨⟨?_, by omega⟩, ⟨c, hc⟩⟩
      rcases Nat.eq_zero_or_pos d with h | h
      · subst h; omega
      · exact h
  have hnd : ¬ N n ∣ m + 1 := by
    rintro ⟨c, hc⟩
    rcases Nat.lt_or_ge c 2 with h | h
    · interval_cases c <;> omega
    · have : N n * 2 ≤ N n * c := Nat.mul_le_mul_left _ h; omega
  have hmu : mu (m + 1) = 1 := by
    have hcop : Nat.Coprime 2 p := (Nat.coprime_primes Nat.prime_two hp).mpr (Ne.symm hp2)
    rw [mu, hm, ArithmeticFunction.isMultiplicative_moebius.map_mul_of_coprime hcop,
      ArithmeticFunction.moebius_apply_prime Nat.prime_two, ArithmeticFunction.moebius_apply_prime hp]
    norm_num
  rw [jump, hf, Finset.sum_erase_eq_sub (Nat.mem_divisors_self _ (by omega)),
    sum_mu_divisors (by omega), hmu, ite_eq_right hnd]
  ring

/-- Bertrand: a prime `p` with `n!/2 < p < n!`, `p ≠ 2`. -/
lemma bertrand_p {n : ℕ} (hn : 3 ≤ n) :
    ∃ p, p.Prime ∧ p ≠ 2 ∧ N n + 2 ≤ 2 * p ∧ p < N n := by
  have h6 := six_le_N hn
  have h2 : 2 ∣ N n := Nat.dvd_factorial (by norm_num) (by omega)
  obtain ⟨p, hp, hlt, hle⟩ := Nat.exists_prime_lt_and_le_two_mul (N n / 2) (by omega)
  have hNe : N n = 2 * (N n / 2) := (Nat.mul_div_cancel' h2).symm
  have hpN : p ≠ N n := by
    rintro rfl
    have := (Nat.prime_dvd_prime_iff_eq Nat.prime_two hp).mp h2
    omega
  exact ⟨p, hp, by omega, by omega, by omega⟩

/-! ## The refutations, for every `n ≥ 3` -/

/-- `a_{n!-1}` is not `n!`-periodic: `a(1) = 1/(n!-1) ≠ 2/(n!-1) = a(1+n!)`. -/
theorem akPeriodic_false {n : ℕ} (hn : 3 ≤ n) : ¬ AkPeriodicClaim n := by
  intro h
  have h6 := six_le_N hn
  have hk : N n - 1 ∈ Ico 1 (N n + 1) := Finset.mem_Ico.mpr ⟨by omega, by omega⟩
  have e := h _ hk 1
  have e1 : ((1 : ℝ) + N n) = ((1 + N n : ℕ) : ℝ) := by push_cast; ring
  have e2 : (1 : ℝ) = ((1 : ℕ) : ℝ) := by norm_num
  rw [e1, a, a, e2, Int.fract_div_natCast_eq_div_natCast_mod,
    Int.fract_div_natCast_eq_div_natCast_mod] at e
  have hmod : (1 + N n) % (N n - 1) = 2 := by
    rw [show 1 + N n = 2 + (N n - 1) * 1 by omega, Nat.add_mul_mod_self_left,
      Nat.mod_eq_of_lt (by omega)]
  rw [hmod, Nat.mod_eq_of_lt (show 1 < N n - 1 by omega)] at e
  have hne : ((N n - 1 : ℕ) : ℝ) ≠ 0 := by exact_mod_cast (show N n - 1 ≠ 0 by omega)
  field_simp at e
  norm_num at e

/-- `φ*_n` is not `n!`-periodic: its jump at `2p` is `-1`, at `2p - n!` it is `0`. -/
theorem not_periodic {n : ℕ} (hn : 3 ≤ n) : ¬ PeriodicClaim n := by
  intro hper
  obtain ⟨p, hp, hp2, h1, h2⟩ := bertrand_p hn
  set m := 2 * p - 1 - N n with hmdef
  have hj0 := jump_small n m (by omega) (by omega)
  have hj1 := jump_two_p n (m + N n) p hp hp2 (by omega) (by omega) h2
  have e0 : phiStar n ((m + N n : ℕ) : ℝ) = phiStar n m := by
    rw [← hper (m : ℝ)]; push_cast; ring_nf
  have e1 : phiStar n ((m + N n + 1 : ℕ) : ℝ) = phiStar n ((m + 1 : ℕ) : ℝ) := by
    rw [← hper ((m + 1 : ℕ) : ℝ)]; push_cast; ring_nf
  rw [e0, e1] at hj1
  linarith

lemma not_exceptional {n t : ℕ} (h1 : N n < t) (h2 : t < 2 * N n) (k : ℕ) :
    ¬ ((N n * k : ℝ) ≤ t ∧ (t : ℝ) < N n * k + 1) := by
  rintro ⟨ha, hb⟩
  have ha' : N n * k ≤ t := by exact_mod_cast ha
  have hb' : t < N n * k + 1 := by exact_mod_cast hb
  rcases Nat.lt_or_ge k 2 with h | h
  · interval_cases k <;> omega
  · have : N n * 2 ≤ N n * k := Nat.mul_le_mul_left _ h; omega

/-- `φ*_n` does NOT coincide with `χ` off the intervals `[n!k, n!k+1)`: on `[2p-1, 2p)` and
`[2p, 2p+1)` (both inside `(n!, 2n!)`) the values differ by the jump `-1`. -/
theorem coincide_false {n : ℕ} (hn : 3 ≤ n) : ¬ CoincideClaim n := by
  intro hc
  obtain ⟨p, hp, hp2, h1, h2⟩ := bertrand_p hn
  have h6 := six_le_N hn
  set m := 2 * p - 1 with hmdef
  have hj := jump_two_p n m p hp hp2 (by omega) (by omega) h2
  have ha := hc (m : ℝ) (by exact_mod_cast (show 1 ≤ m by omega))
    (not_exceptional (by omega) (by omega))
  have hb := hc ((m + 1 : ℕ) : ℝ) (by exact_mod_cast (show 1 ≤ m + 1 by omega))
    (not_exceptional (by omega) (by omega))
  rw [ha, hb] at hj
  norm_num at hj

/-- Summary: the chain is valid; every claim the paper uses to get (4.3) is false for all `n ≥ 3`. -/
theorem ryazanov_gap :
    (TheoremD → Converges → RiemannHypothesis) ∧
      ∀ n, 3 ≤ n → ¬ AkPeriodicClaim n ∧ ¬ PeriodicClaim n ∧ ¬ CoincideClaim n :=
  ⟨ryazanov_RH, fun _ hn => ⟨akPeriodic_false hn, not_periodic hn, coincide_false hn⟩⟩

end Riemann.Deep.Ryazanov
