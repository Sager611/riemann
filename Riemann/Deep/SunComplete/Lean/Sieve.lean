import Riemann.Deep.SunComplete.Lean.Defs
import Riemann.Deep.SunComplete.Lean.Arith
import Mathlib.Data.Nat.Prime.Defs
import Mathlib.Data.Finset.Prod
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

/-!
# sunlean · Sieve: NOTE §2 sieve data, case-(A) constants, and the named / certificate hypotheses

Arithmetic definitions copy Mathlib's verbatim but avoid its heavy modules (memory cap):
`phi` = `Nat.totient`, `primeDivisors m` = `Nat.primeFactors m` (`m ≠ 0`), `PrimePow` =
`isPrimePow_nat_iff`'s right side; `FFInst.lean` proves the agreements.
-/

namespace SunLean

/-- Euler's `φ(m)`; verbatim Mathlib's `Nat.totient m = #{a ∈ range m | m.Coprime a}`. -/
def phi (m : ℕ) : ℕ := ((Finset.range m).filter (fun a => m.Coprime a)).card

/-- The prime divisors of `m` (= Mathlib's `Nat.primeFactors m` for `m ≠ 0`). -/
def primeDivisors (m : ℕ) : Finset ℕ := (Finset.range (m + 1)).filter (fun p => p.Prime ∧ p ∣ m)

/-- `q = p^k` with `p` prime, `k ≥ 1` (Mathlib: `isPrimePow_nat_iff`). -/
def PrimePow (q : ℕ) : Prop := ∃ p k : ℕ, p.Prime ∧ 0 < k ∧ p ^ k = q

/-- Sieve data (NOTE §2): the prime divisors of `q − 1` split into a core (product `k`) and
sieving primes `p_1, …, p_s`. -/
structure SieveData (q : ℕ) where
  core : Finset ℕ
  siev : Finset ℕ
  hsplit : core ∪ siev = primeDivisors (q - 1)
  hdisj : Disjoint core siev

namespace SieveData
variable {q : ℕ} (d : SieveData q)
/-- `θ(k) = ∏_{p | k} (1 − 1/p)`. -/
def θ : ℚ := ∏ p ∈ d.core, (1 - 1 / (p : ℚ))
/-- `δ = 1 − Σ_i 1/p_i`. -/
def δ : ℚ := 1 - ∑ p ∈ d.siev, 1 / (p : ℚ)
/-- `Δ = (s − 1)/δ + 2` if `s ≥ 1`, `Δ = 1` if `s = 0`. -/
def Δ : ℚ := if d.siev.card = 0 then 1 else ((d.siev.card : ℚ) - 1) / d.δ + 2
/-- `W(k) = 2^{ω(k)}`. -/
def W : ℚ := 2 ^ d.core.card
/-- `τ = θ(k)·δ`. -/
def τ : ℚ := d.θ * d.δ
/-- `E = (W(k)Δ − 1)·s + 1` for a rational `s ≥ √q` (NOTE: `E = (W(k)Δ − 1)√q + 1`;
Lemma 3 holds with any larger `E`). -/
def E (s : ℚ) : ℚ := (d.W * d.Δ - 1) * s + 1
end SieveData

section
variable {F : Type} [Field F] [Fintype F] [DecidableEq F]

/-- `e(X, Y)` (NOTE Lemma 3, case (A)): ordered pairs `(x, y) ∈ X × Y` with `x + y` primitive
(`x = y` allowed). -/
noncomputable def pairCount (X Y : Finset F) : ℕ :=
  ((X ×ˢ Y).filter (fun p => p.1 + p.2 ∈ primSet F)).card
end

/-- **Lemma 3, case (A)** (NOTE §2; Cohen's prime sieve with Gauss sums), NAMED HYPOTHESIS.
For every finite field `F` (`n = q = card F`, `Ñ = q − 1`), sieve data with `δ > 0`, rational
`s ≥ √q` (`s ≥ 0`, `q ≤ s²`) and `X, Y ⊆ F`:
`e(X,Y) ≥ (τ/n)·[Ñ|X||Y| − E·√(|X|(n−|X|)·|Y|(n−|Y|))]`, in the equivalent squared form
`L ≤ 0 ∨ L² ≤ (τE)²·|X|(n−|X|)·|Y|(n−|Y|)`, `L = τÑ|X||Y| − n·e(X,Y)`. -/
def Lemma3A : Prop :=
  ∀ (F : Type) [Field F] [Fintype F] [DecidableEq F] (d : SieveData (Fintype.card F)) (s : ℚ),
    0 < d.δ → 0 ≤ s → (Fintype.card F : ℚ) ≤ s ^ 2 → ∀ X Y : Finset F,
      d.τ * (((Fintype.card F : ℚ) - 1) * X.card * Y.card) -
          (Fintype.card F : ℚ) * (pairCount X Y : ℚ) ≤ 0 ∨
      (d.τ * (((Fintype.card F : ℚ) - 1) * X.card * Y.card) -
          (Fintype.card F : ℚ) * (pairCount X Y : ℚ)) ^ 2 ≤
        (d.τ * d.E s) ^ 2 * ((X.card * ((Fintype.card F : ℚ) - X.card)) *
          (Y.card * ((Fintype.card F : ℚ) - Y.card)))

/-- `ε` (NOTE §1, case (A)): `1` for odd `q`, `0` for even `q`. -/
def epsA (q : ℕ) : ℚ := if q % 2 = 0 then 0 else 1

/-- `δ₀ = φ(q − 1) − ε` (NOTE §1): a lower bound for every degree of `H_q`. -/
def delta0A (q : ℕ) : ℚ := (phi (q - 1) : ℚ) - epsA q

/-- `α* = n(ε/τ + E)/(Ñ + E)` in case (A) (`n = q`, `Ñ = q − 1`). -/
def alphaStarA (q : ℕ) (d : SieveData q) (s : ℚ) : ℚ :=
  (q : ℚ) * (epsA q / d.τ + d.E s) / ((q : ℚ) - 1 + d.E s)

/-- Proposition 4's criterion in case (A), for one `q` (verify.py `crit_exact('A', q, …)`;
there `s = sqrt_ub(q)`, `C = ⌊α*⌋ − 1`). -/
def CritA (q : ℕ) (d : SieveData q) (s : ℚ) (C : ℕ) : Prop :=
  0 < d.δ ∧ 0 ≤ s ∧ (q : ℚ) ≤ s ^ 2 ∧ d.E s < (q : ℚ) - 1 ∧
    alphaStarA q d s < C + 2 ∧ 1 ≤ delta0A q - C + 1 ∧
    ((q : ℚ) / 2 < delta0A q - C + 1 ∨
      (0 < Fq ((q : ℚ) - 1) (d.E s) q C (delta0A q - C + 1) ∧
        0 < Fq ((q : ℚ) - 1) (d.E s) q C ((q : ℚ) / 2)))

/-- (CS) of Lemma 5 in case (A), with any `A' ≥ α*` and `D ≤ δ₀ + 1`. -/
def CSA (q : ℕ) (d : SieveData q) (s : ℚ) : Prop :=
  0 < d.δ ∧ 0 ≤ s ∧ (q : ℚ) ≤ s ^ 2 ∧ ∃ A' Dd : ℚ, alphaStarA q d s ≤ A' ∧
    Dd ≤ delta0A q + 1 ∧ A' < Dd ∧ 2 * A' < q ∧
    (d.E s) ^ 2 * Dd < (Dd - A') * ((q : ℚ) - 1) ^ 2 * (1 - 2 * A' / q)

/-- The exceptional set `E_A` of Proposition 7 (NOTE §4, verify.py `EA`; 79 values). -/
def EA : List ℕ := [3, 4, 5, 7, 9, 11, 13, 16, 19, 23, 25, 27, 29, 31, 37, 41, 43, 47, 49, 61, 67,
  71, 73, 79, 103, 121, 127, 131, 139, 151, 157, 169, 181, 191, 199, 211, 223, 229, 241, 271, 277,
  283, 307, 313, 331, 337, 361, 379, 421, 463, 541, 547, 571, 631, 661, 691, 841, 859, 967, 991,
  1021, 1051, 1093, 1171, 1231, 1291, 1321, 1381, 1429, 1471, 1531, 1681, 1951, 2311, 2731, 3571,
  4621, 4831, 6091]

/-- **Proposition 7, case (A)**, CERTIFICATE HYPOTHESIS (verify.py `check_prop4`,
`crit_best('A', q)` over all `2^{ω(q−1)}` partitions, exact `Fraction` arithmetic): every prime
power `7 < q < 30140` outside `E_A` passes Proposition 4's criterion. -/
def Prop7A : Prop :=
  ∀ q : ℕ, PrimePow q → 7 < q → q < 30140 → q ∉ EA →
    ∃ (d : SieveData q) (s : ℚ) (C : ℕ), CritA q d s C

/-- **Lemma 6, case (A)**, NAMED HYPOTHESIS (NOTE §4; verify.py `check_lemma6`): (CS) holds for
every prime power `q ≥ 30140` for some sieve data (NOTE's proof takes the core to be the
`t(ω(q−1))` smallest prime divisors; the hypothesis asks only for existence). -/
def Lemma6A : Prop := ∀ q : ℕ, PrimePow q → 30140 ≤ q → ∃ (d : SieveData q) (s : ℚ), CSA q d s

/-- **§5 certificates, case (A)**, CERTIFICATE HYPOTHESIS (certA.json; verify.py
`check_certificates`: She's construction for 62 primes, explicit cycles for 13 prime powers):
every field with `q ∈ E_A`, `q > 7` elements has the Sun–Hou property. -/
def CertA : Prop :=
  ∀ (F : Type) [Field F] [Fintype F], Fintype.card F ∈ EA → 7 < Fintype.card F → SunHou F

/-- The finite-field facts used; Mathlib proves all three (`FFInst.lean`). -/
structure FFacts (F : Type) [Field F] [Fintype F] : Prop where
  primePow : PrimePow (Fintype.card F)
  card_prim : (primSet F).card = phi (Fintype.card F - 1)
  two_eq_zero : Fintype.card F % 2 = 0 → (2 : F) = 0

end SunLean
