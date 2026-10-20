import Riemann.Deep.Hasse.Height
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Point
import Mathlib.Tactic.FieldSimp

/-!
# Manin's twist and the parallelogram inequality for the height
-/

open Polynomial

namespace Hasse

/-! ## Classical addition formulas over an arbitrary field -/

section General

variable {K : Type*} [Field K]

/-- `x(P+Q) + x(P-Q)` and `x(P+Q)·x(P-Q)` on `Y² = X³ + a₂X² + a₄X + a₆`. -/
lemma sum_prod_X (a₂ a₄ a₆ X Y X₀ Y₀ : K) (hE : Y ^ 2 = X ^ 3 + a₂ * X ^ 2 + a₄ * X + a₆)
    (hE₀ : Y₀ ^ 2 = X₀ ^ 3 + a₂ * X₀ ^ 2 + a₄ * X₀ + a₆) (hX : X ≠ X₀) :
    ((((Y - Y₀) / (X - X₀)) ^ 2 - a₂ - X - X₀) + (((Y + Y₀) / (X - X₀)) ^ 2 - a₂ - X - X₀)) *
        (X - X₀) ^ 2 =
      2 * (X ^ 2 * X₀ + X * X₀ ^ 2 + 2 * a₂ * X * X₀ + a₄ * (X + X₀) + 2 * a₆) ∧
    (((Y - Y₀) / (X - X₀)) ^ 2 - a₂ - X - X₀) * (((Y + Y₀) / (X - X₀)) ^ 2 - a₂ - X - X₀) *
        (X - X₀) ^ 2 =
      (X * X₀ - a₄) ^ 2 - 4 * a₆ * (X + X₀ + a₂) := by
  have hδ : X - X₀ ≠ 0 := sub_ne_zero.mpr hX
  set X₃ := ((Y - Y₀) / (X - X₀)) ^ 2 - a₂ - X - X₀ with hX₃
  set X₄ := ((Y + Y₀) / (X - X₀)) ^ 2 - a₂ - X - X₀ with hX₄
  have e3 : X₃ * (X - X₀) ^ 2 = (Y - Y₀) ^ 2 - (a₂ + X + X₀) * (X - X₀) ^ 2 := by
    rw [hX₃]; field_simp; ring
  have e4 : X₄ * (X - X₀) ^ 2 = (Y + Y₀) ^ 2 - (a₂ + X + X₀) * (X - X₀) ^ 2 := by
    rw [hX₄]; field_simp; ring
  refine ⟨by linear_combination e3 + e4 + 2 * hE + 2 * hE₀, ?_⟩
  set E := X ^ 3 + a₂ * X ^ 2 + a₄ * X + a₆
  set E₀ := X₀ ^ 3 + a₂ * X₀ ^ 2 + a₄ * X₀ + a₆
  set m := (a₂ + X + X₀) * (X - X₀) ^ 2
  have h4 : (X₃ * X₄ * (X - X₀) ^ 2) * (X - X₀) ^ 2 =
      ((X * X₀ - a₄) ^ 2 - 4 * a₆ * (X + X₀ + a₂)) * (X - X₀) ^ 2 := by
    linear_combination (X₄ * (X - X₀) ^ 2) * e3 + ((Y - Y₀) ^ 2 - m) * e4 +
      (Y ^ 2 + E - 2 * Y₀ ^ 2 - 2 * m) * hE + (Y₀ ^ 2 + E₀ - 2 * E - 2 * m) * hE₀
  exact mul_right_cancel₀ (pow_ne_zero 2 hδ) h4

/-- The same formulas after the substitution `X = D u / v`, `X₀ = D T`, `a₂ = c₂ D`,
`a₄ = c₁ D²`, `a₆ = c₀ D³`, divided by `D`. -/
lemma sum_prod_scaled (c₂ c₁ c₀ D u v T Y Y₀ : K) (hD : D ≠ 0) (hv : v ≠ 0)
    (hE : Y ^ 2 = (D * u / v) ^ 3 + c₂ * D * (D * u / v) ^ 2 + c₁ * D ^ 2 * (D * u / v) +
      c₀ * D ^ 3)
    (hE₀ : Y₀ ^ 2 = (D * T) ^ 3 + c₂ * D * (D * T) ^ 2 + c₁ * D ^ 2 * (D * T) + c₀ * D ^ 3)
    (hX : D * u / v ≠ D * T) :
    let X₃ := ((Y - Y₀) / (D * u / v - D * T)) ^ 2 - c₂ * D - D * u / v - D * T
    let X₄ := ((Y + Y₀) / (D * u / v - D * T)) ^ 2 - c₂ * D - D * u / v - D * T
    (u - T * v) ^ 2 * (X₃ / D + X₄ / D) =
      2 * (T * u ^ 2 + T ^ 2 * u * v + 2 * c₂ * T * u * v + c₁ * u * v + c₁ * T * v ^ 2 +
        2 * c₀ * v ^ 2) ∧
    (u - T * v) ^ 2 * (X₃ / D * (X₄ / D)) =
      (T * u - c₁ * v) ^ 2 - 4 * c₀ * v * (u + T * v + c₂ * v) := by
  intro X₃ X₄
  obtain ⟨hs, hp⟩ := sum_prod_X (c₂ * D) (c₁ * D ^ 2) (c₀ * D ^ 3) (D * u / v) Y (D * T) Y₀
    hE hE₀ hX
  constructor
  · calc (u - T * v) ^ 2 * (X₃ / D + X₄ / D)
        = ((X₃ + X₄) * (D * u / v - D * T) ^ 2) * (v ^ 2 / D ^ 3) := by
          field_simp
      _ = 2 * (T * u ^ 2 + T ^ 2 * u * v + 2 * c₂ * T * u * v + c₁ * u * v + c₁ * T * v ^ 2 +
            2 * c₀ * v ^ 2) := by
          rw [hs]; field_simp; ring
  · calc (u - T * v) ^ 2 * (X₃ / D * (X₄ / D))
        = (X₃ * X₄ * (D * u / v - D * T) ^ 2) * (v ^ 2 / D ^ 4) := by
          field_simp
      _ = (T * u - c₁ * v) ^ 2 - 4 * c₀ * v * (u + T * v + c₂ * v) := by
          rw [hp]; field_simp

lemma sum_prod_scaled' (c₂ c₁ c₀ D u v T x Y Y₀ : K) (hD : D ≠ 0) (hv : v ≠ 0)
    (hxuv : x = D * u / v)
    (hE : Y ^ 2 = x ^ 3 + c₂ * D * x ^ 2 + c₁ * D ^ 2 * x + c₀ * D ^ 3)
    (hE₀ : Y₀ ^ 2 = (D * T) ^ 3 + c₂ * D * (D * T) ^ 2 + c₁ * D ^ 2 * (D * T) + c₀ * D ^ 3)
    (hX : x ≠ D * T) :
    (u - T * v) ^ 2 * ((((Y - Y₀) / (x - D * T)) ^ 2 - c₂ * D - x - D * T) / D +
        (((Y + Y₀) / (x - D * T)) ^ 2 - c₂ * D - x - D * T) / D) =
      2 * (T * u ^ 2 + T ^ 2 * u * v + 2 * c₂ * T * u * v + c₁ * u * v + c₁ * T * v ^ 2 +
        2 * c₀ * v ^ 2) ∧
    (u - T * v) ^ 2 * ((((Y - Y₀) / (x - D * T)) ^ 2 - c₂ * D - x - D * T) / D *
        ((((Y + Y₀) / (x - D * T)) ^ 2 - c₂ * D - x - D * T) / D)) =
      (T * u - c₁ * v) ^ 2 - 4 * c₀ * v * (u + T * v + c₂ * v) := by
  subst hxuv
  exact sum_prod_scaled c₂ c₁ c₀ D u v T Y Y₀ hD hv hE hE₀ hX

end General

/-! ## Degree bookkeeping -/

section Degrees

variable {F : Type*} [Field F]

lemma nd_mul {p q : F[X]} {a b : ℕ} (hp : p.natDegree ≤ a) (hq : q.natDegree ≤ b) :
    (p * q).natDegree ≤ a + b := natDegree_mul_le_of_le hp hq

lemma nd_add {p q : F[X]} {n : ℕ} (hp : p.natDegree ≤ n) (hq : q.natDegree ≤ n) :
    (p + q).natDegree ≤ n := natDegree_add_le_of_degree_le hp hq

lemma nd_sub {p q : F[X]} {n : ℕ} (hp : p.natDegree ≤ n) (hq : q.natDegree ≤ n) :
    (p - q).natDegree ≤ n := by
  have := natDegree_sub_le_of_le hp hq; omega

lemma nd_C (c : F) {n : ℕ} : (C c).natDegree ≤ n := by simp

lemma nd_X {n : ℕ} (hn : 1 ≤ n) : (X : F[X]).natDegree ≤ n := natDegree_X_le.trans hn

lemma nd_mono {p : F[X]} {a b : ℕ} (hp : p.natDegree ≤ a) (h : a ≤ b) : p.natDegree ≤ b :=
  hp.trans h

/-- Degree bounds for the three Vieta coefficients. -/
lemma vieta_degrees (c₂ c₁ c₀ : F) (u v : F[X]) (h : ℕ) (hu : u.natDegree ≤ h)
    (hv : v.natDegree ≤ h) :
    ((u - X * v) ^ 2).natDegree ≤ 2 * h + 2 ∧
    (2 * (X * u ^ 2 + X ^ 2 * u * v + 2 * C c₂ * X * u * v + C c₁ * u * v +
        C c₁ * X * v ^ 2 + 2 * C c₀ * v ^ 2)).natDegree ≤ 2 * h + 2 ∧
    ((X * u - C c₁ * v) ^ 2 - 4 * C c₀ * v * (u + X * v + C c₂ * v)).natDegree ≤
      2 * h + 2 := by
  have n2 : (2 : F[X]).natDegree ≤ 0 := by simp
  have n4 : (4 : F[X]).natDegree ≤ 0 := by simp
  have hX1 : (X : F[X]).natDegree ≤ 1 := natDegree_X_le
  have hXv : (X * v).natDegree ≤ h + 1 := by
    have := nd_mul hX1 hv; omega
  have hXu : (X * u).natDegree ≤ h + 1 := by
    have := nd_mul hX1 hu; omega
  refine ⟨?_, ?_, ?_⟩
  · have h1 : (u - X * v).natDegree ≤ h + 1 := nd_sub (nd_mono hu (by omega)) hXv
    have := natDegree_pow_le_of_le 2 h1; omega
  · have hin : (X * u ^ 2 + X ^ 2 * u * v + 2 * C c₂ * X * u * v + C c₁ * u * v +
        C c₁ * X * v ^ 2 + 2 * C c₀ * v ^ 2).natDegree ≤ 2 * h + 2 := ?_
    · exact nd_mono (nd_mul n2 hin) (by omega)
    refine nd_add (nd_add (nd_add (nd_add (nd_add ?_ ?_) ?_) ?_) ?_) ?_
    · exact nd_mono (nd_mul hX1 (natDegree_pow_le_of_le 2 hu)) (by omega)
    · exact nd_mono (nd_mul (nd_mul (natDegree_pow_le_of_le 2 hX1) hu) hv) (by omega)
    · exact nd_mono (nd_mul (nd_mul (nd_mul (nd_mul n2 (nd_C c₂ (n := 0))) hX1)
        hu) hv) (by omega)
    · exact nd_mono (nd_mul (nd_mul (nd_C c₁ (n := 0)) hu) hv) (by omega)
    · exact nd_mono (nd_mul (nd_mul (nd_C c₁ (n := 0)) hX1) (natDegree_pow_le_of_le 2 hv))
        (by omega)
    · exact nd_mono (nd_mul (nd_mul n2 (nd_C c₀ (n := 0)))
        (natDegree_pow_le_of_le 2 hv)) (by omega)
  · refine nd_sub ?_ ?_
    · have h1 : (X * u - C c₁ * v).natDegree ≤ h + 1 :=
        nd_sub hXu (nd_mono (nd_mul (nd_C c₁ (n := 0)) hv) (by omega))
      have := natDegree_pow_le_of_le 2 h1; omega
    · have h1 : (u + X * v + C c₂ * v).natDegree ≤ h + 1 :=
        nd_add (nd_add (nd_mono hu (by omega)) hXv)
          (nd_mono (nd_mul (nd_C c₂ (n := 0)) hv) (by omega))
      exact nd_mono (nd_mul (nd_mul (nd_mul n4 (nd_C c₀ (n := 0))) hv) h1)
        (by omega)

end Degrees

/-! ## The twisted curve over `F(t)` -/

section Twist

open scoped Classical

variable {F : Type*} [Field F]

local notation "ι" => algebraMap F[X] (RatFunc F)

/-- The monic cubic `X³ + c₂X² + c₁X + c₀`. -/
noncomputable def cub (c₂ c₁ c₀ : F) : F[X] := X ^ 3 + C c₂ * X ^ 2 + C c₁ * X + C c₀

lemma cub_monic (c₂ c₁ c₀ : F) : (cub c₂ c₁ c₀).Monic := by
  unfold cub; monicity!

lemma cub_ne_zero (c₂ c₁ c₀ : F) : cub c₂ c₁ c₀ ≠ 0 := (cub_monic c₂ c₁ c₀).ne_zero

lemma cub_natDegree (c₂ c₁ c₀ : F) : (cub c₂ c₁ c₀).natDegree = 3 := by
  unfold cub; compute_degree!

variable (c₂ c₁ c₀ : F)

/-- `D = f(t)` as an element of `F(t)`. -/
noncomputable def Dk : RatFunc F := ι (cub c₂ c₁ c₀)

/-- The variable `t`. -/
noncomputable def Tk : RatFunc F := ι X

/-- Manin's twist `f(t) y² = f(x)`, written as `Y² = X³ + c₂ D X² + c₁ D² X + c₀ D³` with
`X = D x`, `Y = D² y`. -/
noncomputable def twist : WeierstrassCurve (RatFunc F) :=
  ⟨0, ι (C c₂) * Dk c₂ c₁ c₀, 0, ι (C c₁) * Dk c₂ c₁ c₀ ^ 2, ι (C c₀) * Dk c₂ c₁ c₀ ^ 3⟩

lemma Dk_ne_zero : Dk c₂ c₁ c₀ ≠ 0 := RatFunc.algebraMap_ne_zero (cub_ne_zero c₂ c₁ c₀)

lemma Dk_eq : Dk c₂ c₁ c₀ = Tk ^ 3 + ι (C c₂) * Tk ^ 2 + ι (C c₁) * Tk + ι (C c₀) := by
  simp [Dk, Tk, cub]

lemma twist_equation_iff (x y : RatFunc F) : (twist c₂ c₁ c₀).toAffine.Equation x y ↔
    y ^ 2 = x ^ 3 + ι (C c₂) * Dk c₂ c₁ c₀ * x ^ 2 + ι (C c₁) * Dk c₂ c₁ c₀ ^ 2 * x +
      ι (C c₀) * Dk c₂ c₁ c₀ ^ 3 := by
  rw [WeierstrassCurve.Affine.equation_iff]
  simp [twist]

lemma twist_negY (x y : RatFunc F) : (twist c₂ c₁ c₀).toAffine.negY x y = -y := by
  simp [twist]

variable {c₂ c₁ c₀}

lemma four_ne_zero_F (h2 : (2 : F) ≠ 0) : (4 : F) ≠ 0 := by
  have : (4 : F) = 2 * 2 := by norm_num
  rw [this]; exact mul_ne_zero h2 h2

lemma four_ne_zero_X (h2 : (2 : F) ≠ 0) : (4 : F[X]) ≠ 0 := by
  intro h
  apply four_ne_zero_F h2
  have := congrArg (Polynomial.coeff · 0) h
  simpa using this

lemma two_ne_zero_K (h2 : (2 : F) ≠ 0) : (2 : RatFunc F) ≠ 0 := by
  intro h
  apply h2
  apply (algebraMap F (RatFunc F)).injective
  rw [map_ofNat, map_zero, h]

lemma twist_nonsingular {x y : RatFunc F} (h2 : (2 : F) ≠ 0) (hy : y ≠ 0)
    (h : (twist c₂ c₁ c₀).toAffine.Equation x y) : (twist c₂ c₁ c₀).toAffine.Nonsingular x y := by
  rw [WeierstrassCurve.Affine.nonsingular_iff]
  refine ⟨h, Or.inr ?_⟩
  simp only [twist, zero_mul, sub_zero]
  intro hc
  have : (2 : RatFunc F) * y = 0 := by linear_combination hc
  exact (mul_ne_zero (two_ne_zero_K h2) hy) this

lemma P0_equation : (twist c₂ c₁ c₀).toAffine.Equation (Dk c₂ c₁ c₀ * Tk) (Dk c₂ c₁ c₀ ^ 2) := by
  rw [twist_equation_iff]
  linear_combination (Dk c₂ c₁ c₀) ^ 3 * Dk_eq c₂ c₁ c₀

variable (c₂ c₁ c₀)

/-- The point `P₀ = (t, 1)` of the twist (i.e. the identity map of `E`). -/
noncomputable def P0 (h2 : (2 : F) ≠ 0) : (twist c₂ c₁ c₀).toAffine.Point :=
  .some _ _ (twist_nonsingular h2 (pow_ne_zero 2 (Dk_ne_zero c₂ c₁ c₀)) P0_equation)

/-- Height of a point of the twist: projective degree of `x = X / D`. -/
noncomputable def ht : (twist c₂ c₁ c₀).toAffine.Point → ℕ
  | .zero => 0
  | .some x _ _ => htK (x / Dk c₂ c₁ c₀)

variable {c₂ c₁ c₀}

@[simp] lemma ht_zero : ht c₂ c₁ c₀ 0 = 0 := rfl

lemma ht_some {x y : RatFunc F} (h : (twist c₂ c₁ c₀).toAffine.Nonsingular x y) :
    ht c₂ c₁ c₀ (.some x y h) = htK (x / Dk c₂ c₁ c₀) := rfl

@[simp] lemma ht_neg (P : (twist c₂ c₁ c₀).toAffine.Point) :
    ht c₂ c₁ c₀ (-P) = ht c₂ c₁ c₀ P := by
  cases P <;> rfl

lemma ht_P0 (h2 : (2 : F) ≠ 0) : ht c₂ c₁ c₀ (P0 c₂ c₁ c₀ h2) = 1 := by
  rw [P0, ht_some, mul_div_cancel_left₀ _ (Dk_ne_zero c₂ c₁ c₀), Tk, htK_algebraMap,
    natDegree_X]

lemma Tk_ne_Dk_div {x : RatFunc F} (hx : x ≠ Dk c₂ c₁ c₀ * Tk) : x / Dk c₂ c₁ c₀ ≠ Tk := by
  intro h
  apply hx
  rw [← h, mul_div_cancel₀ _ (Dk_ne_zero c₂ c₁ c₀)]

/-- `P₀ + P₀` has height at most 4. -/
lemma ht_double (h2 : (2 : F) ≠ 0) :
    ht c₂ c₁ c₀ (P0 c₂ c₁ c₀ h2 + P0 c₂ c₁ c₀ h2) ≤ 4 := by
  have hD := Dk_ne_zero c₂ c₁ c₀
  have h2K := two_ne_zero_K h2
  have h4K : (4 : RatFunc F) ≠ 0 := by
    have : (4 : RatFunc F) = 2 * 2 := by norm_num
    rw [this]; exact mul_ne_zero h2K h2K
  have hy : Dk c₂ c₁ c₀ ^ 2 ≠ (twist c₂ c₁ c₀).toAffine.negY (Dk c₂ c₁ c₀ * Tk) (Dk c₂ c₁ c₀ ^ 2) := by
    rw [twist_negY]
    intro h
    have : (2 : RatFunc F) * Dk c₂ c₁ c₀ ^ 2 = 0 := by linear_combination h
    exact mul_ne_zero h2K (pow_ne_zero 2 hD) this
  rw [P0, WeierstrassCurve.Affine.Point.add_self_of_Y_ne hy, ht_some,
    WeierstrassCurve.Affine.slope_of_Y_ne rfl hy]
  have hden : Dk c₂ c₁ c₀ ^ 2 - (twist c₂ c₁ c₀).toAffine.negY (Dk c₂ c₁ c₀ * Tk)
      (Dk c₂ c₁ c₀ ^ 2) ≠ 0 := sub_ne_zero.mpr hy
  set ℓ := (3 * (Dk c₂ c₁ c₀ * Tk) ^ 2 + 2 * (twist c₂ c₁ c₀).a₂ * (Dk c₂ c₁ c₀ * Tk) +
        (twist c₂ c₁ c₀).a₄ - (twist c₂ c₁ c₀).a₁ * Dk c₂ c₁ c₀ ^ 2) /
      (Dk c₂ c₁ c₀ ^ 2 - (twist c₂ c₁ c₀).toAffine.negY (Dk c₂ c₁ c₀ * Tk) (Dk c₂ c₁ c₀ ^ 2))
    with hℓdef
  have hℓ : ℓ * 2 = 3 * Tk ^ 2 + 2 * ι (C c₂) * Tk + ι (C c₁) := by
    rw [hℓdef, div_mul_eq_mul_div, div_eq_iff hden, twist_negY]
    simp only [twist]
    ring
  have key : (twist c₂ c₁ c₀).toAffine.addX (Dk c₂ c₁ c₀ * Tk) (Dk c₂ c₁ c₀ * Tk) ℓ /
        Dk c₂ c₁ c₀ =
      ι ((3 * X ^ 2 + 2 * C c₂ * X + C c₁) ^ 2 - 4 * cub c₂ c₁ c₀ * (C c₂ + 2 * X)) /
        ι (4 * cub c₂ c₁ c₀) := by
    rw [div_eq_div_iff hD
      (RatFunc.algebraMap_ne_zero (mul_ne_zero (four_ne_zero_X h2) (cub_ne_zero c₂ c₁ c₀)))]
    simp only [WeierstrassCurve.Affine.addX, twist, map_sub, map_mul, map_add, map_pow,
      map_ofNat]
    rw [← Dk, ← Tk]
    linear_combination (Dk c₂ c₁ c₀ * (2 * ℓ + (3 * Tk ^ 2 + 2 * ι (C c₂) * Tk + ι (C c₁)))) * hℓ
  rw [key]
  refine (htK_le_of_eq_div (mul_ne_zero (four_ne_zero_X h2) (cub_ne_zero c₂ c₁ c₀)) rfl).trans ?_
  refine max_le ?_ ?_
  · unfold cub; compute_degree
  · unfold cub; compute_degree; norm_num

/-- Generic case of the parallelogram inequality. -/
lemma ht_pgram_some (h2 : (2 : F) ≠ 0) {x y : RatFunc F}
    (h : (twist c₂ c₁ c₀).toAffine.Nonsingular x y) (hx : x ≠ Dk c₂ c₁ c₀ * Tk) :
    ht c₂ c₁ c₀ (.some x y h + P0 c₂ c₁ c₀ h2) + ht c₂ c₁ c₀ (.some x y h - P0 c₂ c₁ c₀ h2) ≤
      2 * ht c₂ c₁ c₀ (.some x y h) + 2 := by
  have hD := Dk_ne_zero c₂ c₁ c₀
  rw [sub_eq_add_neg, P0, WeierstrassCurve.Affine.Point.neg_some,
    WeierstrassCurve.Affine.Point.add_of_X_ne hx, WeierstrassCurve.Affine.Point.add_of_X_ne hx,
    ht_some, ht_some, ht_some, WeierstrassCurve.Affine.slope_of_X_ne hx,
    WeierstrassCurve.Affine.slope_of_X_ne hx, twist_negY]
  have hX3 : (twist c₂ c₁ c₀).toAffine.addX x (Dk c₂ c₁ c₀ * Tk)
      ((y - Dk c₂ c₁ c₀ ^ 2) / (x - Dk c₂ c₁ c₀ * Tk)) =
      ((y - Dk c₂ c₁ c₀ ^ 2) / (x - Dk c₂ c₁ c₀ * Tk)) ^ 2 - ι (C c₂) * Dk c₂ c₁ c₀ - x -
        Dk c₂ c₁ c₀ * Tk := by
    simp [twist]
  have hX4 : (twist c₂ c₁ c₀).toAffine.addX x (Dk c₂ c₁ c₀ * Tk)
      ((y - -Dk c₂ c₁ c₀ ^ 2) / (x - Dk c₂ c₁ c₀ * Tk)) =
      ((y + Dk c₂ c₁ c₀ ^ 2) / (x - Dk c₂ c₁ c₀ * Tk)) ^ 2 - ι (C c₂) * Dk c₂ c₁ c₀ - x -
        Dk c₂ c₁ c₀ * Tk := by
    simp [twist]
  rw [hX3, hX4]
  set z := x / Dk c₂ c₁ c₀ with hz
  have hv : ι z.denom ≠ 0 := RatFunc.algebraMap_ne_zero z.denom_ne_zero
  have hxuv : x = Dk c₂ c₁ c₀ * ι z.num / ι z.denom := by
    rw [mul_div_assoc, RatFunc.num_div_denom, hz, mul_div_cancel₀ _ hD]
  have hE : y ^ 2 = x ^ 3 + ι (C c₂) * Dk c₂ c₁ c₀ * x ^ 2 + ι (C c₁) * Dk c₂ c₁ c₀ ^ 2 * x +
      ι (C c₀) * Dk c₂ c₁ c₀ ^ 3 := (twist_equation_iff c₂ c₁ c₀ x y).mp h.1
  have hE₀ := (twist_equation_iff c₂ c₁ c₀ _ _).mp (P0_equation (c₂ := c₂) (c₁ := c₁) (c₀ := c₀))
  obtain ⟨hs, hp⟩ := sum_prod_scaled' (ι (C c₂)) (ι (C c₁)) (ι (C c₀)) (Dk c₂ c₁ c₀) (ι z.num)
    (ι z.denom) Tk x y (Dk c₂ c₁ c₀ ^ 2) hD hv hxuv hE (by rw [← hE₀]) hx
  set x₃ := (((y - Dk c₂ c₁ c₀ ^ 2) / (x - Dk c₂ c₁ c₀ * Tk)) ^ 2 - ι (C c₂) * Dk c₂ c₁ c₀ - x -
    Dk c₂ c₁ c₀ * Tk) / Dk c₂ c₁ c₀ with hx₃
  set x₄ := (((y + Dk c₂ c₁ c₀ ^ 2) / (x - Dk c₂ c₁ c₀ * Tk)) ^ 2 - ι (C c₂) * Dk c₂ c₁ c₀ - x -
    Dk c₂ c₁ c₀ * Tk) / Dk c₂ c₁ c₀ with hx₄
  set A : F[X] := (z.num - X * z.denom) ^ 2 with hA
  set B : F[X] := 2 * (X * z.num ^ 2 + X ^ 2 * z.num * z.denom + 2 * C c₂ * X * z.num * z.denom +
    C c₁ * z.num * z.denom + C c₁ * X * z.denom ^ 2 + 2 * C c₀ * z.denom ^ 2) with hB
  set Cq : F[X] := (X * z.num - C c₁ * z.denom) ^ 2 -
    4 * C c₀ * z.denom * (z.num + X * z.denom + C c₂ * z.denom) with hCq
  have hA0 : A ≠ 0 := by
    refine pow_ne_zero 2 (sub_ne_zero.mpr ?_)
    intro h'
    apply Tk_ne_Dk_div hx
    rw [← hz, ← RatFunc.num_div_denom z, h', map_mul, mul_div_cancel_right₀ _ hv, Tk]
  have hsA : ι A * (x₃ + x₄) = ι B := by
    rw [hA, hB]
    simp only [map_mul, map_sub, map_add, map_pow, map_ofNat]
    rw [← Tk]
    linear_combination hs
  have hpA : ι A * (x₃ * x₄) = ι Cq := by
    rw [hA, hCq]
    simp only [map_mul, map_sub, map_add, map_pow, map_ofNat]
    rw [← Tk]
    linear_combination hp
  refine (htK_add_htK_le hA0 hsA hpA).trans ?_
  obtain ⟨d1, d2, d3⟩ := vieta_degrees c₂ c₁ c₀ z.num z.denom (htK z) (le_max_left _ _)
    (le_max_right _ _)
  exact max_le d1 (max_le d2 d3)

/-- **Parallelogram inequality** for the height on the twist. -/
theorem ht_pgram (h2 : (2 : F) ≠ 0) (P : (twist c₂ c₁ c₀).toAffine.Point) :
    ht c₂ c₁ c₀ (P + P0 c₂ c₁ c₀ h2) + ht c₂ c₁ c₀ (P - P0 c₂ c₁ c₀ h2) ≤
      2 * ht c₂ c₁ c₀ P + 2 := by
  rcases P with _ | ⟨x, y, h⟩
  · show ht c₂ c₁ c₀ (0 + P0 c₂ c₁ c₀ h2) + ht c₂ c₁ c₀ (0 - P0 c₂ c₁ c₀ h2) ≤
      2 * ht c₂ c₁ c₀ 0 + 2
    rw [zero_add, zero_sub, ht_neg, ht_P0, ht_zero]
  · by_cases hx : x = Dk c₂ c₁ c₀ * Tk
    · rcases WeierstrassCurve.Affine.Y_eq_of_X_eq h.1 P0_equation hx with hy | hy
      · have hP : (WeierstrassCurve.Affine.Point.some x y h) = P0 c₂ c₁ c₀ h2 := by
          subst hx hy; rfl
        rw [hP, sub_self, ht_zero, ht_P0]
        have := ht_double (c₂ := c₂) (c₁ := c₁) (c₀ := c₀) h2
        omega
      · have hP : (WeierstrassCurve.Affine.Point.some x y h) = -P0 c₂ c₁ c₀ h2 := by
          subst hx hy; rfl
        rw [hP, neg_add_cancel, ht_zero, ← neg_add', ht_neg, ht_neg, ht_P0]
        have := ht_double (c₂ := c₂) (c₁ := c₁) (c₀ := c₀) h2
        omega
    · exact ht_pgram_some h2 h hx

end Twist

end Hasse
