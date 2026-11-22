import Mathlib.RingTheory.Polynomial.ShiftedLegendre
import Mathlib.Algebra.Polynomial.FieldDivision
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.Complex.Basic
import Mathlib.Topology.Algebra.Polynomial

/-!
# Christoffel's endpoint bound for polynomials on [0, 1]

deg P < m ⇒ |P(0)|² ≤ m² ∫_0^1 |P|². Proof: shifted Legendre polynomials (Mathlib's
Polynomial.shiftedLegendre, Rodrigues' formula there) are orthogonal on [0,1] with
∫ L̃_n² = 1/(2n+1) and L̃_n(0) = 1 (ours, by iterated integration by parts); the kernel
K_m = Σ_{n<m} (2n+1) L̃_n reproduces P(0) and has ∫ K_m² = m².
-/

open Polynomial MeasureTheory intervalIntegral Nat

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- ∫_0^1 P(x) dx for a complex polynomial evaluated at real points. -/
def pint (P : ℂ[X]) : ℂ := ∫ x in (0 : ℝ)..1, P.eval (x : ℂ)

theorem hasDerivAt_peval (P : ℂ[X]) (x : ℝ) :
    HasDerivAt (fun x : ℝ => P.eval (x : ℂ)) (P.derivative.eval (x : ℂ)) x :=
  (P.hasDerivAt (x : ℂ)).comp_ofReal

theorem continuous_peval (P : ℂ[X]) : Continuous (fun x : ℝ => P.eval (x : ℂ)) :=
  P.continuous.comp Complex.continuous_ofReal

theorem pint_add (P Q : ℂ[X]) : pint (P + Q) = pint P + pint Q := by
  unfold pint; simp only [eval_add]
  exact integral_add ((continuous_peval P).intervalIntegrable _ _)
    ((continuous_peval Q).intervalIntegrable _ _)

theorem pint_sub (P Q : ℂ[X]) : pint (P - Q) = pint P - pint Q := by
  unfold pint; simp only [eval_sub]
  exact integral_sub ((continuous_peval P).intervalIntegrable _ _)
    ((continuous_peval Q).intervalIntegrable _ _)

theorem pint_C_mul (c : ℂ) (P : ℂ[X]) : pint (C c * P) = c * pint P := by
  unfold pint; simp only [eval_mul, eval_C]
  exact integral_const_mul c _

theorem pint_neg (P : ℂ[X]) : pint (-P) = -pint P := by
  unfold pint; simp only [eval_neg]; exact intervalIntegral.integral_neg

theorem pint_C (c : ℂ) : pint (C c) = c := by
  unfold pint; simp

theorem pint_finsetSum {ι : Type*} (s : Finset ι) (f : ι → ℂ[X]) :
    pint (∑ i ∈ s, f i) = ∑ i ∈ s, pint (f i) := by
  unfold pint; simp only [eval_finsetSum]
  exact intervalIntegral.integral_finset_sum fun i _ => (continuous_peval (f i)).intervalIntegrable _ _

theorem pint_ibp (F G : ℂ[X]) :
    pint (F * derivative G) = F.eval 1 * G.eval 1 - F.eval 0 * G.eval 0 - pint (derivative F * G) := by
  unfold pint
  simp only [eval_mul]
  have := intervalIntegral.integral_mul_deriv_eq_deriv_mul (a := (0 : ℝ)) (b := 1)
    (u := fun x : ℝ => F.eval (x : ℂ)) (v := fun x : ℝ => G.eval (x : ℂ))
    (u' := fun x : ℝ => F.derivative.eval (x : ℂ)) (v' := fun x : ℝ => G.derivative.eval (x : ℂ))
    (fun x _ => hasDerivAt_peval F x) (fun x _ => hasDerivAt_peval G x)
    ((continuous_peval _).intervalIntegrable _ _) ((continuous_peval _).intervalIntegrable _ _)
  simpa using this

theorem pint_ibp_iter (G : ℂ[X]) (k : ℕ)
    (hG : ∀ j < k, (derivative^[j] G).eval 0 = 0 ∧ (derivative^[j] G).eval 1 = 0) (F : ℂ[X]) :
    pint (F * derivative^[k] G) = (-1) ^ k * pint (derivative^[k] F * G) := by
  induction k generalizing F with
  | zero => simp
  | succ k ih =>
    rw [Function.iterate_succ_apply' derivative k G, pint_ibp F (derivative^[k] G)]
    have h0 := hG k (Nat.lt_succ_self k)
    rw [h0.1, h0.2, mul_zero, mul_zero, sub_zero, zero_sub,
      ih (fun j hj => hG j (lt_trans hj (Nat.lt_succ_self k))) (derivative F),
      ← Function.iterate_succ_apply derivative k F]
    ring

/-- X^n (1 - X)^n. -/
def rod (n : ℕ) : ℂ[X] := X ^ n * (1 - X) ^ n

theorem rod_vanish (n j : ℕ) (hj : j < n) :
    (derivative^[j] (rod n)).eval 0 = 0 ∧ (derivative^[j] (rod n)).eval 1 = 0 := by
  constructor
  · have hdvd : (X - C (0 : ℂ)) ^ n ∣ rod n := by
      rw [map_zero, sub_zero]; exact dvd_mul_right _ _
    obtain ⟨q, hq⟩ := pow_sub_dvd_iterate_derivative_of_pow_dvd j hdvd
    rw [hq, eval_mul, eval_pow, eval_sub, eval_X, eval_C, sub_self, zero_pow (by omega), zero_mul]
  · have hdvd : (X - C (1 : ℂ)) ^ n ∣ rod n := by
      have e : rod n = (C ((-1 : ℂ) ^ n) * X ^ n) * (X - C 1) ^ n := by
        unfold rod
        have h1 : (1 - X : ℂ[X]) = C (-1) * (X - C 1) := by
          rw [mul_sub, ← C_mul, mul_one, map_neg, map_one, neg_one_mul]; ring
        rw [h1, mul_pow, ← C_pow]; ring
      rw [e]; exact dvd_mul_left _ _
    obtain ⟨q, hq⟩ := pow_sub_dvd_iterate_derivative_of_pow_dvd j hdvd
    rw [hq, eval_mul, eval_pow, eval_sub, eval_X, eval_C, sub_self, zero_pow (by omega), zero_mul]

/-- Shifted Legendre polynomial over ℝ and over ℂ (Mathlib's Polynomial.shiftedLegendre). -/
def sLr (n : ℕ) : ℝ[X] := (shiftedLegendre n).map (Int.castRingHom ℝ)
def sL (n : ℕ) : ℂ[X] := (sLr n).map Complex.ofRealHom

theorem sL_eq_map (n : ℕ) : sL n = (shiftedLegendre n).map (Int.castRingHom ℂ) := by
  unfold sL sLr
  rw [Polynomial.map_map]
  congr 1

theorem sL_eval_ofReal (n : ℕ) (x : ℝ) : (sL n).eval (x : ℂ) = (((sLr n).eval x : ℝ) : ℂ) := by
  unfold sL
  exact Polynomial.eval_map_apply (p := sLr n) (f := Complex.ofRealHom) x

theorem sL_coeff (n k : ℕ) : (sL n).coeff k =
    (((-1 : ℤ) ^ k * n.choose k * (n + k).choose n : ℤ) : ℂ) := by
  rw [sL_eq_map, coeff_map, coeff_shiftedLegendre]; simp

theorem sL_eval_zero (n : ℕ) : (sL n).eval 0 = 1 := by
  rw [← coeff_zero_eq_eval_zero, sL_coeff]; simp

theorem sL_coeff_self (n : ℕ) : (sL n).coeff n = (-1 : ℂ) ^ n * ((2 * n).choose n : ℂ) := by
  rw [sL_coeff]; push_cast; rw [Nat.choose_self, two_mul]; push_cast; ring

theorem sL_coeff_gt (n k : ℕ) (hk : n < k) : (sL n).coeff k = 0 := by
  rw [sL_coeff, Nat.choose_eq_zero_of_lt hk]; simp

theorem sL_degree_lt (n m : ℕ) (h : n < m) : (sL n).degree < m := by
  rw [degree_lt_iff_coeff_zero]
  intro k hk
  exact sL_coeff_gt n k (by omega)

theorem sL_rodrigues (n : ℕ) : C (n ! : ℂ) * sL n = derivative^[n] (rod n) := by
  have h := congrArg (Polynomial.map (Int.castRingHom ℂ)) (factorial_mul_shiftedLegendre_eq n)
  rw [Polynomial.map_mul, ← iterate_derivative_map] at h
  rw [sL_eq_map]
  unfold rod
  convert h using 2
  · simp
  · simp [Polynomial.map_mul, Polynomial.map_pow, Polynomial.map_sub]

/-- Orthogonality: ∫ F L̃_n = 0 for deg F < n. -/
theorem pint_mul_sL_eq_zero (F : ℂ[X]) (n : ℕ) (hF : F.degree < n) : pint (F * sL n) = 0 := by
  have h1 := pint_ibp_iter (rod n) n (rod_vanish n) F
  rw [iterate_derivative_eq_zero_of_degree_lt hF, zero_mul, ← sL_rodrigues] at h1
  have e : F * (C (n ! : ℂ) * sL n) = C (n ! : ℂ) * (F * sL n) := by ring
  rw [e, pint_C_mul] at h1
  have hf : (n ! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos n).ne'
  have h2 : pint (0 : ℂ[X]) = 0 := by unfold pint; simp
  rw [h2, mul_zero] at h1
  exact (mul_eq_zero.mp h1).resolve_left hf

theorem pint_X_pow (a : ℕ) : pint (X ^ a) = 1 / (a + 1) := by
  unfold pint
  simp only [eval_pow, eval_X]
  have : ∫ x in (0 : ℝ)..1, (x : ℂ) ^ a = ((∫ x in (0 : ℝ)..1, x ^ a : ℝ) : ℂ) := by
    rw [← intervalIntegral.integral_ofReal]; push_cast; rfl
  rw [this, integral_pow]; push_cast; ring

/-- Beta integral ∫_0^1 x^a (1 - x)^b = a! b! / (a + b + 1)!. -/
theorem pint_beta (b : ℕ) : ∀ a : ℕ, pint (X ^ a * (1 - X) ^ b) = (a ! * b ! : ℂ) / ((a + b + 1)! : ℂ) := by
  induction b with
  | zero =>
    intro a
    rw [pow_zero, mul_one, pint_X_pow, Nat.factorial_zero, add_zero, Nat.factorial_succ]
    have : (a ! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos a).ne'
    push_cast; field_simp
  | succ b ih =>
    intro a
    set G : ℂ[X] := C ((a + 1 : ℂ)⁻¹) * X ^ (a + 1)
    have hG : derivative G = X ^ a := by
      simp only [G, derivative_C_mul, derivative_X_pow]
      rw [← mul_assoc, ← C_mul]
      have : (a + 1 : ℂ)⁻¹ * ((a + 1 : ℕ) : ℂ) = 1 := by push_cast; field_simp
      rw [this, C_1, one_mul]; simp
    have hF : derivative ((1 - X : ℂ[X]) ^ (b + 1)) = -(C ((b + 1 : ℕ) : ℂ) * (1 - X) ^ b) := by
      rw [derivative_pow, derivative_sub, derivative_one, derivative_X]
      simp
    have h := pint_ibp ((1 - X) ^ (b + 1)) G
    rw [hG, hF] at h
    have e1 : X ^ a * (1 - X : ℂ[X]) ^ (b + 1) = (1 - X) ^ (b + 1) * X ^ a := by ring
    rw [e1, h]
    simp only [G, eval_pow, eval_sub, eval_one, eval_X, sub_self, eval_mul, eval_C, zero_pow
      (Nat.succ_ne_zero _), zero_mul, sub_zero, zero_sub, ne_eq]
    have e2 : -(C ((b + 1 : ℕ) : ℂ) * (1 - X) ^ b) * (C ((a + 1 : ℂ)⁻¹) * X ^ (a + 1)) =
        -(C (((b + 1 : ℕ) : ℂ) * (a + 1 : ℂ)⁻¹) * (X ^ (a + 1) * (1 - X) ^ b)) := by
      rw [C_mul]; ring
    rw [e2, pint_neg, pint_C_mul, ih (a + 1)]
    have h3 : a + 1 + b + 1 = a + (b + 1) + 1 := by ring
    rw [h3]
    have hfa : (a ! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos a).ne'
    have hfb : (b ! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos b).ne'
    have hfab : ((a + (b + 1) + 1) ! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos _).ne'
    rw [Nat.factorial_succ a, Nat.factorial_succ b]
    push_cast
    field_simp
    ring

theorem pint_rod (n : ℕ) : pint (rod n) = ((n !) ^ 2 : ℂ) / ((2 * n + 1)! : ℂ) := by
  unfold rod
  rw [pint_beta]; congr 1
  · ring
  · congr 2; ring

/-- ∫ L̃_n² = 1/(2n+1). -/
theorem pint_sL_sq (n : ℕ) : pint (sL n * sL n) = 1 / (2 * n + 1) := by
  set lead := (sL n).coeff n
  set R := sL n - C lead * X ^ n
  have hR : R.degree < n := by
    rw [degree_lt_iff_coeff_zero]
    intro k hk
    simp only [R, coeff_sub, coeff_C_mul, coeff_X_pow]
    rcases eq_or_lt_of_le hk with h | h
    · rw [← h]; simp [lead]
    · rw [sL_coeff_gt n k h, if_neg (by omega)]; ring
  have hsplit : sL n * sL n = C lead * (X ^ n * sL n) + R * sL n := by
    simp only [R]; ring
  rw [hsplit, pint_add, pint_C_mul, pint_mul_sL_eq_zero R n hR, add_zero]
  -- ∫ X^n L̃_n = (-1)^n ∫ rod
  have h1 := pint_ibp_iter (rod n) n (rod_vanish n) (X ^ n)
  rw [← sL_rodrigues, iterate_derivative_X_pow_eq_smul, Nat.sub_self, pow_zero,
    Nat.descFactorial_self] at h1
  have e : X ^ n * (C (n ! : ℂ) * sL n) = C (n ! : ℂ) * (X ^ n * sL n) := by ring
  rw [e, pint_C_mul] at h1
  have e2 : ((n ! : ℂ) • (1 : ℂ[X])) * rod n = C (n ! : ℂ) * rod n := by
    rw [smul_eq_C_mul, mul_one]
  rw [e2, pint_C_mul, pint_rod] at h1
  have hf : (n ! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos n).ne'
  have hXs : pint (X ^ n * sL n) = (-1) ^ n * (((n !) ^ 2 : ℂ) / ((2 * n + 1)! : ℂ)) := by
    have := mul_left_cancel₀ hf (h1.trans (by ring : (-1 : ℂ) ^ n * ((n ! : ℂ) *
      (((n !) ^ 2 : ℂ) / ((2 * n + 1)! : ℂ))) = (n ! : ℂ) * ((-1) ^ n * (((n !) ^ 2 : ℂ) / ((2 * n + 1)! : ℂ)))))
    exact this
  rw [hXs, show lead = (-1 : ℂ) ^ n * ((2 * n).choose n : ℂ) from sL_coeff_self n]
  have hc : ((2 * n).choose n : ℂ) * ((n !) ^ 2 : ℂ) = ((2 * n)! : ℂ) := by
    have := Nat.choose_mul_factorial_mul_factorial (show n ≤ 2 * n by omega)
    rw [show 2 * n - n = n by omega] at this
    have h2 : (2 * n).choose n * (n !) ^ 2 = (2 * n)! := by rw [← this]; ring
    exact_mod_cast h2
  have hf2 : ((2 * n + 1)! : ℂ) = (2 * n + 1) * ((2 * n)! : ℂ) := by
    rw [Nat.factorial_succ]; push_cast; ring
  have hf3 : ((2 * n)! : ℂ) ≠ 0 := by exact_mod_cast (Nat.factorial_pos _).ne'
  have hodd : (2 * (n : ℂ) + 1) ≠ 0 := by
    have : (0 : ℝ) < 2 * n + 1 := by positivity
    exact_mod_cast this.ne'
  rw [hf2]
  have hsq : ((-1 : ℂ) ^ n) * ((-1 : ℂ) ^ n) = 1 := by rw [← mul_pow]; simp
  calc (-1 : ℂ) ^ n * ((2 * n).choose n : ℂ) * ((-1) ^ n * (((n !) ^ 2 : ℂ) / ((2 * n + 1) * ((2 * n)! : ℂ))))
      = ((-1 : ℂ) ^ n * (-1) ^ n) * (((2 * n).choose n : ℂ) * ((n !) ^ 2 : ℂ)) /
          ((2 * n + 1) * ((2 * n)! : ℂ)) := by ring
    _ = 1 / (2 * n + 1) := by rw [hsq, hc]; field_simp

/-- K_m = Σ_{n<m} (2n+1) L̃_n. -/
def Kpoly (m : ℕ) : ℂ[X] := ∑ n ∈ Finset.range m, C (2 * (n : ℂ) + 1) * sL n

theorem Kpoly_degree_lt (m : ℕ) : (Kpoly m).degree < m := by
  rw [degree_lt_iff_coeff_zero]
  intro k hk
  unfold Kpoly
  rw [finset_sum_coeff]
  refine Finset.sum_eq_zero fun n hn => ?_
  rw [Finset.mem_range] at hn
  rw [coeff_C_mul, sL_coeff_gt n k (by omega), mul_zero]

theorem sL_coeff_self_ne (m : ℕ) : (sL m).coeff m ≠ 0 := by
  rw [sL_coeff_self]
  refine mul_ne_zero (pow_ne_zero _ (by norm_num)) ?_
  exact_mod_cast (Nat.choose_pos (by omega : m ≤ 2 * m)).ne'

theorem two_mul_add_one_ne (m : ℕ) : (2 * (m : ℂ) + 1) ≠ 0 := by
  have : (0 : ℝ) < 2 * m + 1 := by positivity
  exact_mod_cast this.ne'

/-- The reproducing property ∫ P K_m = P(0) for deg P < m. -/
theorem pint_mul_K (m : ℕ) : ∀ P : ℂ[X], P.degree < m → pint (P * Kpoly m) = P.eval 0 := by
  induction m with
  | zero =>
    intro P hP
    have hP0 : P = 0 := by
      ext k; rw [degree_lt_iff_coeff_zero] at hP; simpa using hP k (Nat.zero_le k)
    subst hP0; simp [Kpoly, pint]
  | succ m ih =>
    intro P hP
    set lead := (sL m).coeff m
    have hlead : lead ≠ 0 := sL_coeff_self_ne m
    set c := P.coeff m / lead
    set P' := P - C c * sL m
    have hP' : P'.degree < m := by
      rw [degree_lt_iff_coeff_zero]
      intro k hk
      rw [degree_lt_iff_coeff_zero] at hP
      simp only [P', coeff_sub, coeff_C_mul]
      rcases eq_or_lt_of_le hk with h | h
      · rw [← h]; simp only [c]; field_simp; ring
      · rw [hP k (by omega), sL_coeff_gt m k h]; ring
    have hK : Kpoly (m + 1) = Kpoly m + C (2 * (m : ℂ) + 1) * sL m := by
      unfold Kpoly; rw [Finset.sum_range_succ]
    have hdec : P = C c * sL m + P' := by simp only [P']; ring
    have e : P * Kpoly (m + 1) = C c * (Kpoly m * sL m) + C (c * (2 * (m : ℂ) + 1)) * (sL m * sL m) +
        P' * Kpoly m + C (2 * (m : ℂ) + 1) * (P' * sL m) := by
      rw [hK]; nth_rewrite 1 [hdec]; rw [C_mul]; ring
    rw [e, pint_add, pint_add, pint_add, pint_C_mul, pint_C_mul, pint_C_mul,
      pint_mul_sL_eq_zero (Kpoly m) m (Kpoly_degree_lt m), pint_sL_sq, ih P' hP',
      pint_mul_sL_eq_zero P' m hP']
    have h2 := two_mul_add_one_ne m
    have hev : P.eval 0 = c + P'.eval 0 := by
      simp only [P', eval_sub, eval_mul, eval_C, sL_eval_zero]; ring
    rw [hev]
    field_simp
    ring

theorem pint_K_sq (m : ℕ) : pint (Kpoly m * Kpoly m) = (m : ℂ) ^ 2 := by
  induction m with
  | zero => simp [Kpoly, pint]
  | succ m ih =>
    have hK : Kpoly (m + 1) = Kpoly m + C (2 * (m : ℂ) + 1) * sL m := by
      unfold Kpoly; rw [Finset.sum_range_succ]
    have e : Kpoly (m + 1) * Kpoly (m + 1) = Kpoly m * Kpoly m +
        C (2 * (2 * (m : ℂ) + 1)) * (Kpoly m * sL m) + C ((2 * (m : ℂ) + 1) ^ 2) * (sL m * sL m) := by
      rw [hK, C_mul, C_pow, map_ofNat C 2]; ring
    rw [e, pint_add, pint_add, pint_C_mul, pint_C_mul, ih,
      pint_mul_sL_eq_zero (Kpoly m) m (Kpoly_degree_lt m), pint_sL_sq]
    have h2 := two_mul_add_one_ne m
    push_cast
    field_simp
    ring

/-- The real-valued kernel. -/
def Kr (m : ℕ) (x : ℝ) : ℝ := ∑ n ∈ Finset.range m, (2 * (n : ℝ) + 1) * (sLr n).eval x

theorem Kpoly_eval_real (m : ℕ) (x : ℝ) : (Kpoly m).eval (x : ℂ) = ((Kr m x : ℝ) : ℂ) := by
  unfold Kpoly Kr
  rw [eval_finsetSum]
  push_cast
  refine Finset.sum_congr rfl fun n _ => ?_
  rw [eval_mul, eval_C, sL_eval_ofReal]

theorem continuous_Kr (m : ℕ) : Continuous (Kr m) := by
  unfold Kr
  exact continuous_finsetSum _ fun n _ => continuous_const.mul (sLr n).continuous

theorem integral_Kr_sq (m : ℕ) : ∫ x in (0 : ℝ)..1, (Kr m x) ^ 2 = (m : ℝ) ^ 2 := by
  have h := pint_K_sq m
  unfold pint at h
  simp_rw [eval_mul, Kpoly_eval_real] at h
  have : ∫ x in (0 : ℝ)..1, ((Kr m x : ℝ) : ℂ) * ((Kr m x : ℝ) : ℂ) =
      ((∫ x in (0 : ℝ)..1, (Kr m x) ^ 2 : ℝ) : ℂ) := by
    rw [← intervalIntegral.integral_ofReal]; push_cast; congr 1; ext x; ring
  rw [this] at h
  exact_mod_cast h

/-- **Christoffel's endpoint bound**: deg P < m ⇒ |P(0)|² ≤ m² ∫_0^1 |P|². -/
theorem christoffel (m : ℕ) (P : ℂ[X]) (hP : P.degree < m) :
    ‖P.eval 0‖ ^ 2 ≤ (m : ℝ) ^ 2 * ∫ x in (0 : ℝ)..1, ‖P.eval (x : ℂ)‖ ^ 2 := by
  rcases Nat.eq_zero_or_pos m with hm | hm
  · subst hm
    have hP0 : P = 0 := by
      ext k; rw [degree_lt_iff_coeff_zero] at hP; simpa using hP k (Nat.zero_le k)
    subst hP0; simp
  set f : ℝ → ℝ := fun x => ‖P.eval (x : ℂ)‖
  set g : ℝ → ℝ := fun x => |Kr m x|
  have hf : Continuous f := (continuous_peval P).norm
  have hg : Continuous g := (continuous_Kr m).abs
  set A := ∫ x in (0 : ℝ)..1, f x ^ 2
  set B := ∫ x in (0 : ℝ)..1, f x * g x
  set Cc := ∫ x in (0 : ℝ)..1, g x ^ 2
  have hCc : Cc = (m : ℝ) ^ 2 := by
    simp only [Cc, g, sq_abs]; exact integral_Kr_sq m
  have hCpos : 0 < Cc := by rw [hCc]; positivity
  -- |P(0)| ≤ B
  have h1 : ‖P.eval 0‖ ≤ B := by
    have hrep := pint_mul_K m P hP
    unfold pint at hrep
    simp_rw [eval_mul, Kpoly_eval_real] at hrep
    rw [← hrep]
    refine (intervalIntegral.norm_integral_le_integral_norm (by norm_num)).trans (le_of_eq ?_)
    refine intervalIntegral.integral_congr fun x _ => ?_
    simp only [f, g, norm_mul, Complex.norm_real, Real.norm_eq_abs]
  -- Cauchy–Schwarz by the quadratic trick
  have hq : ∀ t : ℝ, 0 ≤ A - 2 * t * B + t ^ 2 * Cc := by
    intro t
    have hnn : 0 ≤ ∫ x in (0 : ℝ)..1, (f x - t * g x) ^ 2 :=
      intervalIntegral.integral_nonneg (by norm_num) fun x _ => sq_nonneg _
    have hi1 : IntervalIntegrable (fun x => f x ^ 2) volume 0 1 := (hf.pow 2).intervalIntegrable _ _
    have hi2 : IntervalIntegrable (fun x => f x * g x) volume 0 1 := (hf.mul hg).intervalIntegrable _ _
    have hi3 : IntervalIntegrable (fun x => g x ^ 2) volume 0 1 := (hg.pow 2).intervalIntegrable _ _
    have e : ∫ x in (0 : ℝ)..1, (f x - t * g x) ^ 2 = A - 2 * t * B + t ^ 2 * Cc := by
      have e1 : (fun x => (f x - t * g x) ^ 2) =
          fun x => f x ^ 2 - (2 * t) * (f x * g x) + t ^ 2 * g x ^ 2 := by ext x; ring
      rw [e1, intervalIntegral.integral_add (hi1.sub (hi2.const_mul _)) (hi3.const_mul _),
        intervalIntegral.integral_sub hi1 (hi2.const_mul _), intervalIntegral.integral_const_mul,
        intervalIntegral.integral_const_mul]
    rw [← e]; exact hnn
  have hB2 : B ^ 2 ≤ A * Cc := by
    have := hq (B / Cc)
    have e : A - 2 * (B / Cc) * B + (B / Cc) ^ 2 * Cc = A - B ^ 2 / Cc := by field_simp; ring
    rw [e] at this
    have : B ^ 2 / Cc ≤ A := by linarith
    rw [div_le_iff₀ hCpos] at this
    linarith
  have hB0 : 0 ≤ B := (norm_nonneg _).trans h1
  calc ‖P.eval 0‖ ^ 2 ≤ B ^ 2 := pow_le_pow_left₀ (norm_nonneg _) h1 2
    _ ≤ A * Cc := hB2
    _ = (m : ℝ) ^ 2 * A := by rw [hCc]; ring

end Riemann.Deep.LeanNBRate



