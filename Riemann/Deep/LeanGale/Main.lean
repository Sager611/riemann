import Riemann.Deep.LeanGale.Gale
import Riemann.Deep.LeanGale.Family
import Riemann.Deep.LeanGale.P3
import Riemann.Deep.LeanGale.P4
import Riemann.Deep.LeanGale.P7
import Riemann.Deep.LeanGale.P8

/-!
# LeanGale: the main theorems (DISCOVERIES #2, discoveries/gale-r7-r8.md)

* `claim1_n7`, `claim1_n8`: the explicit sets points7.txt (84 points, ℝ⁷) and points8.txt (85 points, ℝ⁸) cannot
  be split into n + 1 parts of diameter ≤ σₙ · D; β₈/D = √(8900856812318928/10000000213244525) = 0.9434435130…
  > σ₇ = 0.9419651451… and β₉/D = √(9023994457279939/10000000157507010) = 0.9499470677… > σ₈ = 0.9486832981…
* `claim2`: for every n = 5…1000 and all unit simplices U ⊂ ℝ^⌊n/2⌋, W ⊂ ℝ^⌈n/2⌉, the family member Vₙ has
  (⌊n/2⌋+1) + (⌈n/2⌉+1) + 3(⌊n/2⌋+1)(⌈n/2⌉+1) points, cannot be split into n + 1 parts of diameter ≤ σₙ · D, and
  β_{n+1}(Vₙ)/D(Vₙ) > σₙ (unit simplices exist in every dimension: `splx_simplex`).
* `gale_fails_le_1000`: Gale's conjecture fails in ℝⁿ for every n = 3…1000 (no hypotheses; n = 3, 4 from our
  sets of discoveries/gale-r3.md and gale-r4-r5.md, kernel-checked certificates).
* `gale_fails`: Gale's conjecture fails in every dimension n ≥ 3, given `Hinrichs2002` (Borsuk's conjecture fails
  in all dimensions > 322) for n > 1000.
* `gale_ball_fails`: the same for the "ball is extremal" form, given also `SimplexPartitionBound`.
-/

namespace LeanGale

open Gale7

/-! ## General lemmas -/

theorem beta_eq {k n : ℕ} (hk : 0 < k) {X : Finset (Fin n → ℝ)} {t : ℝ}
    (hlow : ∀ c : (Fin n → ℝ) → Fin k, ∃ x ∈ X, ∃ y ∈ X, c x = c y ∧ t ≤ sqd x y)
    (hup : ∃ c : (Fin n → ℝ) → Fin k, ∀ x ∈ X, ∀ y ∈ X, c x = c y → sqd x y ≤ t) :
    beta k (XE X) = Real.sqrt t :=
  le_antisymm (beta_le_sqrt hk hup) (sqrt_le_beta hk hlow)

theorem coe_image_toE {n : ℕ} (X : Finset (Fin n → ℝ)) : ((X.image toE : Finset (E n)) : Set (E n)) = XE X := by
  rw [Finset.coe_image]
  rfl

theorem card_image_toE {n : ℕ} (X : Finset (Fin n → ℝ)) : (X.image toE).card = X.card :=
  Finset.card_image_of_injective _ toE_injective

/-- β_{n+1}(V)/diam V > σₙ from σₙ · diam V < β_{n+1}(V) -/
theorem sigma_lt_ratio {n : ℕ} {V : Set (E n)} (hV : Bornology.IsBounded V)
    (h : sigma n * Metric.diam V < beta (n + 1) V) : sigma n < beta (n + 1) V / Metric.diam V := by
  have hle := beta_le_diam (Nat.succ_pos n) hV
  have hD : 0 < Metric.diam V := by
    by_contra hD
    have h0 : Metric.diam V = 0 := le_antisymm (not_lt.mp hD) Metric.diam_nonneg
    rw [h0, mul_zero] at h
    linarith
  rw [lt_div_iff₀ hD]
  exact h

/-- the σ test from integers, as reals -/
theorem sigma_test {n t D : ℕ} (hn : 1 ≤ n) (h : sigmaOK n t D = true) :
    sigma n ^ 2 * ((D : ℝ) / 10 ^ 16) < (t : ℝ) / 10 ^ 16 := sigma_sq_mul_lt' hn h

/-- the exact ratio √(dc²/10¹⁶) / √(D²/10¹⁶) = √(dc²/D²) -/
theorem ratio_sqrt (a b : ℝ) (ha : 0 ≤ a) :
    Real.sqrt (a / 10 ^ 16) / Real.sqrt (b / 10 ^ 16) = Real.sqrt (a / b) := by
  rw [← Real.sqrt_div (div_nonneg ha (by positivity)), div_div_div_cancel_right₀ (by positivity)]

/-- σ < √r from σ² · D < t, r = t/D -/
theorem sigma_lt_sqrt {n : ℕ} {t D : ℝ} (hD : 0 < D) (h : sigma n ^ 2 * D < t) :
    sigma n < Real.sqrt (t / D) := by
  rw [Real.lt_sqrt (sigma_nonneg n), lt_div_iff₀ hD]
  exact h

theorem sqrt_near {r d e : ℝ} (he0 : 0 < e) (he : e ≤ d) (h1 : (d - e) ^ 2 < r) (h2 : r < (d + e) ^ 2) :
    |Real.sqrt r - d| < e := by
  have hde : 0 ≤ d - e := by linarith
  rw [abs_sub_lt_iff]
  constructor
  · have := (Real.sqrt_lt' (by linarith : (0 : ℝ) < d + e)).mpr h2
    linarith
  · have := (Real.lt_sqrt hde).mpr h1
    linarith

/-! ## Claim (1): the sets points7.txt and points8.txt -/

/-- the 84-point set points7.txt ⊂ ℝ⁷ -/
noncomputable def V7 : Finset (E 7) := P7.X7.image toE

theorem hσ7 : sigma 7 ^ 2 * ((10000000213244525 : ℝ) / 10 ^ 16) < (8900856812318928 : ℝ) / 10 ^ 16 := by
  have := sigma_test (n := 7) (t := 8900856812318928) (D := 10000000213244525) (by norm_num) (by decide)
  exact_mod_cast this

theorem V7_key : sigma 7 * Metric.diam (V7 : Set (E 7)) < beta 8 (V7 : Set (E 7)) := by
  rw [V7, coe_image_toE]
  exact sigma_mul_diam_lt_beta (fun x hx y hy => P7.sqd_le x y hx hy) P7.lower hσ7

/-- **Claim (1), ℝ⁷.** The 84 points of points7.txt cannot be split into 8 parts of diameter ≤ σ₇ · D;
β₈/D = √(8900856812318928/10000000213244525) = 0.9434435130… > σ₇ = 0.9419651451… -/
theorem claim1_n7 :
    V7.card = 84 ∧
    (∀ c : E 7 → Fin 8, ∃ k, sigma 7 * Metric.diam (V7 : Set (E 7)) <
      Metric.diam ((V7 : Set (E 7)) ∩ c ⁻¹' {k})) ∧
    beta 8 (V7 : Set (E 7)) / Metric.diam (V7 : Set (E 7)) =
      Real.sqrt (8900856812318928 / 10000000213244525) ∧
    sigma 7 < Real.sqrt (8900856812318928 / 10000000213244525) ∧
    |Real.sqrt (8900856812318928 / 10000000213244525) - 0.9434435130| < 5 / 10 ^ 11 ∧
    |sigma 7 - 0.9419651451| < 5 / 10 ^ 11 := by
  refine ⟨?_, fun c => exists_part_gt (by norm_num) V7_key c, ?_, ?_, ?_, ?_⟩
  · rw [V7, card_image_toE, P7.card_X7]
  · rw [V7, coe_image_toE, beta_eq (by norm_num) P7.lower P7.upper, diam_XE, P7.diam_X7,
      ratio_sqrt _ _ (by norm_num)]
  · have := sigma_lt_sqrt (by norm_num) hσ7
    rwa [div_div_div_cancel_right₀ (by norm_num)] at this
  · exact sqrt_near (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  · have h7 : sigma 7 = Real.sqrt (1 / 2 + Real.sqrt (3 / 5) / 2) := by
      rw [sigma, if_neg (by norm_num)]
      norm_num
    have hs : |Real.sqrt (3 / 5) - 0.774596669241| < 1 / 10 ^ 12 :=
      sqrt_near (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    rw [abs_sub_lt_iff] at hs
    rw [h7]
    apply sqrt_near (by norm_num) (by norm_num)
    · have := hs.2
      norm_num at this ⊢
      linarith
    · have := hs.1
      norm_num at this ⊢
      linarith

/-- the 85-point set points8.txt ⊂ ℝ⁸ -/
noncomputable def V8 : Finset (E 8) := P8.X8.image toE

theorem hσ8 : sigma 8 ^ 2 * ((10000000157507010 : ℝ) / 10 ^ 16) < (9023994457279939 : ℝ) / 10 ^ 16 := by
  have := sigma_test (n := 8) (t := 9023994457279939) (D := 10000000157507010) (by norm_num) (by decide)
  exact_mod_cast this

theorem V8_key : sigma 8 * Metric.diam (V8 : Set (E 8)) < beta 9 (V8 : Set (E 8)) := by
  rw [V8, coe_image_toE]
  exact sigma_mul_diam_lt_beta (fun x hx y hy => P8.sqd_le x y hx hy) P8.lower hσ8

/-- **Claim (1), ℝ⁸.** The 85 points of points8.txt cannot be split into 9 parts of diameter ≤ σ₈ · D;
β₉/D = √(9023994457279939/10000000157507010) = 0.9499470677… > σ₈ = 0.9486832981… -/
theorem claim1_n8 :
    V8.card = 85 ∧
    (∀ c : E 8 → Fin 9, ∃ k, sigma 8 * Metric.diam (V8 : Set (E 8)) <
      Metric.diam ((V8 : Set (E 8)) ∩ c ⁻¹' {k})) ∧
    beta 9 (V8 : Set (E 8)) / Metric.diam (V8 : Set (E 8)) =
      Real.sqrt (9023994457279939 / 10000000157507010) ∧
    sigma 8 < Real.sqrt (9023994457279939 / 10000000157507010) ∧
    |Real.sqrt (9023994457279939 / 10000000157507010) - 0.9499470677| < 5 / 10 ^ 11 ∧
    |sigma 8 - 0.9486832981| < 5 / 10 ^ 11 := by
  refine ⟨?_, fun c => exists_part_gt (by norm_num) V8_key c, ?_, ?_, ?_, ?_⟩
  · rw [V8, card_image_toE, P8.card_X8]
  · rw [V8, coe_image_toE, beta_eq (by norm_num) P8.lower P8.upper, diam_XE, P8.diam_X8,
      ratio_sqrt _ _ (by norm_num)]
  · have := sigma_lt_sqrt (by norm_num) hσ8
    rwa [div_div_div_cancel_right₀ (by norm_num)] at this
  · exact sqrt_near (by norm_num) (by norm_num) (by norm_num) (by norm_num)
  · rw [sigma, if_pos (by norm_num)]
    exact sqrt_near (by norm_num) (by norm_num) (by norm_num) (by norm_num)

/-! ## Claim (2): the family for every n = 5…1000 -/

theorem splitEq (n : ℕ) : n / 2 + (n - n / 2) = n := by omega

/-- reindexing ℝ^k = ℝⁿ along k = n -/
def reidx {k n : ℕ} (h : k = n) (x : Fin k → ℝ) : Fin n → ℝ := fun t => x (Fin.cast h.symm t)

theorem transport {k n : ℕ} (h : k = n) {X : Finset (Fin k → ℝ)} {t D2 : ℝ}
    (hD : ∀ x ∈ X, ∀ y ∈ X, sqd x y ≤ D2)
    (hlow : ∀ c : (Fin k → ℝ) → Fin (k + 1), ∃ x ∈ X, ∃ y ∈ X, c x = c y ∧ t ≤ sqd x y) :
    (∀ x ∈ X.image (reidx h), ∀ y ∈ X.image (reidx h), sqd x y ≤ D2) ∧
    (∀ c : (Fin n → ℝ) → Fin (n + 1), ∃ x ∈ X.image (reidx h), ∃ y ∈ X.image (reidx h),
      c x = c y ∧ t ≤ sqd x y) ∧
    (X.image (reidx h)).card = X.card := by
  subst h
  have e : X.image (reidx (rfl : k = k)) = X := by
    rw [show reidx (rfl : k = k) = id from rfl, Finset.image_id]
  rw [e]
  exact ⟨hD, hlow, rfl⟩

/-- **The family member Vₙ ⊂ ℝⁿ** (NOTE `2): ℝⁿ = ℝ^⌊n/2⌋ ⊕ ℝ^⌈n/2⌉ (the first ⌊n/2⌋ coordinates, then the last
⌈n/2⌉), U, W the vertex vectors of the two regular simplices, parameters family_params.txt row n:
P_i = r_A u_i, Q_j = r_B v_j, A_ij = a u_i − b v_j, B_ij = −c u_i + d v_j, C_ij = e u_i + f v_j. -/
noncomputable def famV (n : ℕ) (U : Fin (n / 2 + 1) → Fin (n / 2) → ℝ)
    (W : Fin (n - n / 2 + 1) → Fin (n - n / 2) → ℝ) : Finset (E n) :=
  ((famSet n U W).image (reidx (splitEq n))).image toE

theorem famV_facts {n : ℕ} (h5 : 5 ≤ n) (h1000 : n ≤ 1000)
    {U : Fin (n / 2 + 1) → Fin (n / 2) → ℝ} {W : Fin (n - n / 2 + 1) → Fin (n - n / 2) → ℝ}
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) :
    (famV n U W).card = (n / 2 + 1) + (n - n / 2 + 1) + 3 * ((n / 2 + 1) * (n - n / 2 + 1)) ∧
    sigma n * Metric.diam (famV n U W : Set (E n)) < beta (n + 1) (famV n U W : Set (E n)) := by
  obtain ⟨hD, hlow, hσ⟩ := rowOK_facts (params_ok h5 h1000) hU hW
  rw [splitEq n] at hσ
  change ∀ x ∈ famSet n U W, ∀ y ∈ famSet n U W, _ at hD
  change ∀ c : (Fin (n / 2 + (n - n / 2)) → ℝ) → Fin (n / 2 + (n - n / 2) + 1),
    ∃ x ∈ famSet n U W, ∃ y ∈ famSet n U W, _ at hlow
  obtain ⟨hD', hlow', hcard⟩ := transport (splitEq n) hD hlow
  refine ⟨?_, ?_⟩
  · rw [famV, card_image_toE, hcard]
    exact card_famSet h5 h1000 hU hW
  · rw [famV, coe_image_toE]
    exact sigma_mul_diam_lt_beta hD' hlow' hσ

/-- **Claim (2).** For every n = 5…1000 and all unit simplices U ⊂ ℝ^⌊n/2⌋, W ⊂ ℝ^⌈n/2⌉ (NOTE `2), Vₙ has
(⌊n/2⌋+1) + (⌈n/2⌉+1) + 3(⌊n/2⌋+1)(⌈n/2⌉+1) points, cannot be split into n + 1 parts of diameter ≤ σₙ · D, and
β_{n+1}(Vₙ)/D(Vₙ) > σₙ; the colouring part is the hand lemma of NOTE `3 (Gale7 `Vx.not_colourable`), the
metric part the exact integer row check `rowOK`. -/
theorem claim2 (n : ℕ) (h5 : 5 ≤ n) (h1000 : n ≤ 1000)
    (U : Fin (n / 2 + 1) → Fin (n / 2) → ℝ) (W : Fin (n - n / 2 + 1) → Fin (n - n / 2) → ℝ)
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) :
    (famV n U W).card = (n / 2 + 1) + (n - n / 2 + 1) + 3 * ((n / 2 + 1) * (n - n / 2 + 1)) ∧
    (∀ c : E n → Fin (n + 1), ∃ k, sigma n * Metric.diam (famV n U W : Set (E n)) <
      Metric.diam ((famV n U W : Set (E n)) ∩ c ⁻¹' {k})) ∧
    sigma n < beta (n + 1) (famV n U W : Set (E n)) / Metric.diam (famV n U W : Set (E n)) := by
  obtain ⟨hc, hk⟩ := famV_facts h5 h1000 hU hW
  exact ⟨hc, fun c => exists_part_gt (Nat.succ_pos n) hk c,
    sigma_lt_ratio (Finset.finite_toSet _).isBounded hk⟩

/-! ## Every dimension -/

theorem hσ3 : sigma 3 ^ 2 * ((10000000211897483 : ℝ) / 10 ^ 16) < (7932058815501876 : ℝ) / 10 ^ 16 := by
  have := sigma_test (n := 3) (t := 7932058815501876) (D := 10000000211897483) (by norm_num) (by decide)
  exact_mod_cast this

theorem hσ4 : sigma 4 ^ 2 * ((10000000179265409 : ℝ) / 10 ^ 16) < (8399124363373578 : ℝ) / 10 ^ 16 := by
  have := sigma_test (n := 4) (t := 8399124363373578) (D := 10000000179265409) (by norm_num) (by decide)
  exact_mod_cast this

/-- **Gale's conjecture fails in ℝⁿ for every n = 3…1000** (no hypotheses). -/
theorem gale_fails_le_1000 : ∀ n : ℕ, 3 ≤ n → n ≤ 1000 → ¬ GaleConjecture n := by
  intro n h3 h1000
  rcases (show n = 3 ∨ n = 4 ∨ 5 ≤ n by omega) with rfl | rfl | h5
  · exact not_gale_of (isBounded_XE _)
      (sigma_mul_diam_lt_beta (fun x hx y hy => P3.sqd_le x y hx hy) P3.lower hσ3)
  · exact not_gale_of (isBounded_XE _)
      (sigma_mul_diam_lt_beta (fun x hx y hy => P4.sqd_le x y hx hy) P4.lower hσ4)
  · exact not_gale_of (Finset.finite_toSet _).isBounded
      (famV_facts h5 h1000 (splx_simplex (by omega)) (splx_simplex (by omega))).2

/-- **Literature input** (A. Hinrichs, "Spherical codes and Borsuk's conjecture", Discrete Math. 243 (2002)
253–256, Theorem 1, p. 253): "There exists a finite set in the unit sphere in R^323 which cannot be partitioned
into 561 sets of smaller diameter. Hence Borsuk's conjecture fails in all dimensions exceeding 322." Borsuk's
conjecture, as stated there: "every bounded set S ⊂ R^d containing at least two points can be partitioned into
at most d + 1 sets of smaller diameter". Formalized: for every d > 322 some bounded S ⊂ ℝ^d with at least two
points has, for every partition into d + 1 parts (the classes of a colouring; empty parts allowed, i.e. at
most d + 1 parts), a part whose diameter is not smaller than diam S. -/
def Hinrichs2002 : Prop :=
  ∀ d : ℕ, 322 < d → ∃ S : Set (E d), Bornology.IsBounded S ∧ S.Nontrivial ∧
    ∀ c : E d → Fin (d + 1), ∃ k, ¬ Metric.diam (S ∩ c ⁻¹' {k}) < Metric.diam S

/-- a Borsuk counterexample is a Gale counterexample (σₙ < 1) -/
theorem not_gale_of_borsuk {d : ℕ} {S : Set (E d)} (hS : Bornology.IsBounded S) (hnt : S.Nontrivial)
    (hB : ∀ c : E d → Fin (d + 1), ∃ k, ¬ Metric.diam (S ∩ c ⁻¹' {k}) < Metric.diam S) :
    ¬ GaleConjecture d := by
  intro hG
  have hpos : 0 < Metric.diam S := by
    obtain ⟨x, hx, y, hy, hxy⟩ := hnt
    exact lt_of_lt_of_le (dist_pos.mpr hxy) (Metric.dist_le_diam_of_mem hS hx hy)
  have h2 : beta (d + 1) S < Metric.diam S := by
    have := hG S hS
    have := sigma_lt_one d
    nlinarith
  haveI : Nonempty (E d → Fin (d + 1)) := ⟨fun _ => 0⟩
  obtain ⟨c, hc⟩ := exists_lt_of_ciInf_lt h2
  obtain ⟨k, hk⟩ := hB c
  exact hk (lt_of_le_of_lt (le_iSup_part S c k) hc)

/-- **Gale's conjecture (1953) fails in every dimension n ≥ 3**: n = 3…1000 by `gale_fails_le_1000`,
n > 1000 by Borsuk counterexamples (`Hinrichs2002`). -/
theorem gale_fails (hHin : Hinrichs2002) : ∀ n : ℕ, 3 ≤ n → ¬ GaleConjecture n := by
  intro n h3
  by_cases h : n ≤ 1000
  · exact gale_fails_le_1000 n h3 h
  · obtain ⟨S, hS, hnt, hB⟩ := hHin n (by omega)
    exact not_gale_of_borsuk hS hnt hB

/-- **Literature input** (the regular-simplex partition of the ball; V. Voronov, "A counterexample to Gale's
conjecture in R^3", arXiv 2610.28828 (unrefereed), main_en.tex l. 96–107): "Central projection of the facets of a
regular n-simplex inscribed in S^{n−1} partitions the sphere into n+1 congruent regular spherical simplices
(taking cones from the center gives the corresponding partition of the ball). … This gives the following upper
bound for n ≥ 2: β_{n+1}(B^n) = β_{n+1}(S^{n−1}) ≤ σ_n" (β there divided by the diameter; here B has diameter 1). -/
def SimplexPartitionBound : Prop :=
  ∀ n : ℕ, 2 ≤ n → beta (n + 1) (Metric.closedBall (0 : E n) (1 / 2)) ≤ sigma n

/-- **The "ball is extremal" form fails in every dimension n ≥ 3** (given the two literature inputs). -/
theorem gale_ball_fails (hHin : Hinrichs2002) (hSB : SimplexPartitionBound) :
    ∀ n : ℕ, 3 ≤ n → ¬ GaleConjectureBall n := by
  intro n h3 hGB
  apply gale_fails hHin n h3
  intro V hV
  exact (hGB V hV).trans (mul_le_mul_of_nonneg_right (hSB n (by omega)) Metric.diam_nonneg)

end LeanGale
