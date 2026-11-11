import Mathlib.Analysis.Complex.Polynomial.Basic
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.Convex.Topology
import Mathlib.Analysis.Convex.Hull

/-!
# Explicit Kim–Lee wedge, part 1: definitions, geometry, Hurwitz, the Obreschkoff step

DERIVATION.md §2 (parabola, convexity, sector) and §3 (zero-free region ⇒ hyperbolicity), with
Obreschkoff's theorem, the Hadamard partial products and Gauss–Lucas for genus 0 as named hypotheses.
Hurwitz's theorem for real polynomials of bounded degree ((I5), polynomial case) is proved here.
-/

open Polynomial Filter Topology

namespace KimEff

noncomputable section

/-- A real polynomial is hyperbolic if all its complex zeros are real.
(The zero polynomial vanishes at I, so it is not hyperbolic.) -/
def Hyperbolic (p : ℝ[X]) : Prop :=
  ∀ z : ℂ, (p.map (algebraMap ℝ ℂ)).eval z = 0 → z.im = 0

/-- sum over j ≤ d of binom(d, j) b_j X^j. -/
def jensenPoly (b : ℕ → ℝ) (d : ℕ) : ℝ[X] :=
  ∑ j ∈ Finset.range (d + 1), C ((d.choose j : ℝ) * b j) * X ^ j

/-- Jensen polynomial J^{d,n}(X) = sum_{j ≤ d} binom(d,j) γ(n+j) X^j of a real sequence γ. -/
def jensen (γ : ℕ → ℝ) (d n : ℕ) : ℝ[X] := jensenPoly (fun j => γ (n + j)) d

/-- J(P; d)(X) = sum_{j ≤ d} binom(d,j) P^{(j)}(0) X^j for a real polynomial P. -/
def jensenOf (P : ℝ[X]) (d : ℕ) : ℝ[X] := jensenPoly (fun j => (derivative^[j] P).eval 0) d

lemma coeff_jensenPoly (b : ℕ → ℝ) (d j : ℕ) :
    (jensenPoly b d).coeff j = if j ≤ d then (d.choose j : ℝ) * b j else 0 := by
  rw [jensenPoly, finsetSum_coeff]
  simp_rw [coeff_C_mul_X_pow]
  rw [Finset.sum_ite_eq]
  simp only [Finset.mem_range, Nat.lt_succ_iff]

lemma natDegree_jensenPoly_le (b : ℕ → ℝ) (d : ℕ) : (jensenPoly b d).natDegree ≤ d := by
  unfold jensenPoly
  refine natDegree_sum_le_of_forall_le _ _ (fun j hj => ?_)
  exact (natDegree_C_mul_X_pow_le _ _).trans (Nat.lt_succ_iff.mp (Finset.mem_range.mp hj))

/-! ## Hurwitz for real polynomials of bounded degree ((I5), polynomial case) -/

lemma norm_prod_sub_ge (z : ℂ) (s : Multiset ℂ) (hs : ∀ r ∈ s, r.im = 0) :
    |z.im| ^ Multiset.card s ≤ ‖(s.map fun r => z - r).prod‖ := by
  induction s using Multiset.induction_on with
  | empty => simp
  | cons r s ih =>
    have hr : r.im = 0 := hs r (Multiset.mem_cons_self r s)
    have ih' := ih (fun x hx => hs x (Multiset.mem_cons_of_mem hx))
    rw [Multiset.map_cons, Multiset.prod_cons, norm_mul, Multiset.card_cons, pow_succ]
    have h1 : |z.im| ≤ ‖z - r‖ := by
      have := Complex.abs_im_le_norm (z - r)
      simpa [Complex.sub_im, hr] using this
    calc |z.im| ^ Multiset.card s * |z.im| ≤ ‖(s.map fun r => z - r).prod‖ * ‖z - r‖ :=
          mul_le_mul ih' h1 (abs_nonneg _) (norm_nonneg _)
      _ = ‖z - r‖ * ‖(s.map fun r => z - r).prod‖ := mul_comm _ _

/-- A hyperbolic real polynomial satisfies |lc| |Im z|^deg ≤ |p(z)| on ℂ. -/
lemma norm_eval_ge_of_hyperbolic {p : ℝ[X]} (hp : Hyperbolic p) (z : ℂ) :
    |p.leadingCoeff| * |z.im| ^ p.natDegree ≤ ‖(p.map (algebraMap ℝ ℂ)).eval z‖ := by
  set q := p.map (algebraMap ℝ ℂ) with hq
  have hcard : Multiset.card q.roots = q.natDegree := IsAlgClosed.card_roots_eq_natDegree
  have hdeg : q.natDegree = p.natDegree := natDegree_map _
  have hlc : q.leadingCoeff = (p.leadingCoeff : ℂ) := by
    rw [hq, leadingCoeff_map]; rfl
  have hreal : ∀ r ∈ q.roots, r.im = 0 := fun r hr => hp r (isRoot_of_mem_roots hr)
  have heq := C_leadingCoeff_mul_prod_multiset_X_sub_C hcard
  have hev : q.eval z = q.leadingCoeff * (q.roots.map fun r => z - r).prod := by
    conv_lhs => rw [← heq]
    rw [eval_mul, eval_C, eval_multiset_prod, Multiset.map_map]
    congr 2
    exact Multiset.map_congr rfl (fun r _ => by simp)
  rw [hev, norm_mul, hlc, Complex.norm_real, Real.norm_eq_abs, ← hdeg, ← hcard]
  exact mul_le_mul_of_nonneg_left (norm_prod_sub_ge z _ hreal) (abs_nonneg _)

/-- Hurwitz, polynomial case: a coefficientwise limit of hyperbolic real polynomials of degree
at most d whose d-th coefficient is nonzero is hyperbolic. -/
theorem hyperbolic_of_tendsto {p : ℕ → ℝ[X]} {q : ℝ[X]} {d : ℕ}
    (hdeg : ∀ k, (p k).natDegree ≤ d) (hq : q.natDegree ≤ d) (hlead : q.coeff d ≠ 0)
    (hlim : ∀ j, Tendsto (fun k => (p k).coeff j) atTop (𝓝 (q.coeff j)))
    (hhyp : ∀ᶠ k in atTop, Hyperbolic (p k)) : Hyperbolic q := by
  intro z hz
  by_contra him
  have hm : 0 < |z.im| ^ d := pow_pos (abs_pos.mpr him) d
  have hL : 0 < |q.coeff d| := abs_pos.mpr hlead
  have hsum : ∀ r : ℝ[X], r.natDegree ≤ d →
      (r.map (algebraMap ℝ ℂ)).eval z =
        ∑ j ∈ Finset.range (d + 1), ((r.coeff j : ℝ) : ℂ) * z ^ j := by
    intro r hr
    rw [eval_eq_sum_range' (n := d + 1) (natDegree_map_le.trans_lt (Nat.lt_succ_of_le hr))]
    simp [coeff_map]
  have hev : Tendsto (fun k => ((p k).map (algebraMap ℝ ℂ)).eval z) atTop (𝓝 0) := by
    have h : Tendsto (fun k => ∑ j ∈ Finset.range (d + 1), (((p k).coeff j : ℝ) : ℂ) * z ^ j)
        atTop (𝓝 (∑ j ∈ Finset.range (d + 1), ((q.coeff j : ℝ) : ℂ) * z ^ j)) := by
      refine tendsto_finsetSum _ (fun j _ => ?_)
      exact ((Complex.continuous_ofReal.tendsto _).comp (hlim j)).mul_const _
    rw [← hsum q hq, hz] at h
    refine h.congr' (Eventually.of_forall fun k => ?_)
    exact (hsum (p k) (hdeg k)).symm
  have hev' : Tendsto (fun k => ‖((p k).map (algebraMap ℝ ℂ)).eval z‖) atTop (𝓝 0) := by
    simpa using hev.norm
  have h1 : ∀ᶠ k in atTop, |q.coeff d| / 2 < |(p k).coeff d| :=
    ((continuous_abs.tendsto _).comp (hlim d)).eventually (lt_mem_nhds (half_lt_self hL))
  have h2 : ∀ᶠ k in atTop, ‖((p k).map (algebraMap ℝ ℂ)).eval z‖ < |q.coeff d| / 2 * |z.im| ^ d :=
    hev'.eventually (gt_mem_nhds (mul_pos (half_pos hL) hm))
  obtain ⟨k, hk1, hk2, hk3⟩ := (h1.and (h2.and hhyp)).exists
  have hcd : (p k).coeff d ≠ 0 := by
    intro h
    rw [h, abs_zero] at hk1
    linarith [half_pos hL]
  have hnd : (p k).natDegree = d := natDegree_eq_of_le_of_coeff_ne_zero (hdeg k) hcd
  have hlc : (p k).leadingCoeff = (p k).coeff d := by rw [leadingCoeff, hnd]
  have hge := norm_eval_ge_of_hyperbolic hk3 z
  rw [hlc, hnd] at hge
  have hlt : |q.coeff d| / 2 * |z.im| ^ d < |(p k).coeff d| * |z.im| ^ d :=
    mul_lt_mul_of_pos_right hk1 hm
  linarith

/-! ## Geometry (DERIVATION §2.1, §3) -/

/-- The double sector S(δ) = {z : |Im z| ≤ δ |z|} (Kim–Lee, Thm 3). -/
def sector (δ : ℝ) : Set ℂ := {z | |z.im| ≤ δ * ‖z‖}

/-- The parabola P_a = {w : (Im w)^2 ≤ 4a^2(a^2 - Re w)} (DERIVATION (2.1)). -/
def parabola (a : ℝ) : Set ℂ := {w | w.im ^ 2 ≤ 4 * a ^ 2 * (a ^ 2 - w.re)}

/-- W_R = {w ∈ P_a : Re w ≥ -R} (DERIVATION §3). -/
def wRegion (a R : ℝ) : Set ℂ := {w | w ∈ parabola a ∧ -R ≤ w.re}

/-- δ(R) = 2a √(R + a^2) / R (DERIVATION §3). -/
def sectorDelta (a R : ℝ) : ℝ := 2 * a * Real.sqrt (R + a ^ 2) / R

/-- D(R) = R^2 / (4a^2 (R + a^2)) = δ(R)^{-2} (DERIVATION §3). -/
def Dfun (a R : ℝ) : ℝ := R ^ 2 / (4 * a ^ 2 * (R + a ^ 2))

/-- (2.1): |Re z| ≤ a implies z^2 ∈ P_a. -/
lemma sq_mem_parabola {a : ℝ} {z : ℂ} (hz : |z.re| ≤ a) : z ^ 2 ∈ parabola a := by
  have hx : z.re ^ 2 ≤ a ^ 2 := sq_le_sq' (neg_le_of_abs_le hz) (le_of_abs_le hz)
  have h1 : z.re ^ 2 * z.im ^ 2 ≤ a ^ 2 * z.im ^ 2 :=
    mul_le_mul_of_nonneg_right hx (sq_nonneg _)
  have h2 : 0 ≤ a ^ 2 * (a ^ 2 - z.re ^ 2) := mul_nonneg (sq_nonneg a) (sub_nonneg.mpr hx)
  have hre : (z ^ 2).re = z.re ^ 2 - z.im ^ 2 := by rw [sq, Complex.mul_re]; ring
  have him : (z ^ 2).im = 2 * z.re * z.im := by rw [sq, Complex.mul_im]; ring
  show (z ^ 2).im ^ 2 ≤ 4 * a ^ 2 * (a ^ 2 - (z ^ 2).re)
  rw [hre, him]
  nlinarith

/-- (2.1): P_a is convex. -/
lemma convex_parabola (a : ℝ) : Convex ℝ (parabola a) := by
  intro x hx y hy s t hs ht hst
  simp only [parabola, Set.mem_ofPred_eq] at hx hy ⊢
  simp only [Complex.add_re, Complex.add_im, Complex.smul_re, Complex.smul_im, smul_eq_mul]
  have ht' : t = 1 - s := by linarith
  subst ht'
  have key : (s * x.im + (1 - s) * y.im) ^ 2 ≤ s * x.im ^ 2 + (1 - s) * y.im ^ 2 := by
    have h : s * x.im ^ 2 + (1 - s) * y.im ^ 2 - (s * x.im + (1 - s) * y.im) ^ 2
        = s * (1 - s) * (x.im - y.im) ^ 2 := by ring
    nlinarith [mul_nonneg (mul_nonneg hs ht) (sq_nonneg (x.im - y.im))]
  nlinarith [mul_le_mul_of_nonneg_left hx hs, mul_le_mul_of_nonneg_left hy ht]

/-- (2.1): P_a is closed. -/
lemma isClosed_parabola (a : ℝ) : IsClosed (parabola a) := by
  unfold parabola
  exact isClosed_le (Complex.continuous_im.pow 2)
    (continuous_const.mul (continuous_const.sub Complex.continuous_re))

/-- §3: a point of P_a with Re w < -R lies in the sector S(δ(R)). -/
lemma mem_sector_of_parabola {a R : ℝ} (ha : 0 < a) (hR : 0 < R) {w : ℂ}
    (hw : w ∈ parabola a) (hre : w.re < -R) : w ∈ sector (sectorDelta a R) := by
  simp only [parabola, Set.mem_ofPred_eq] at hw
  simp only [sector, Set.mem_ofPred_eq, sectorDelta]
  have hX : 0 < R + a ^ 2 := add_pos hR (pow_pos ha 2)
  have hsq : Real.sqrt (R + a ^ 2) ^ 2 = R + a ^ 2 := Real.sq_sqrt hX.le
  have hδ : 0 ≤ 2 * a * Real.sqrt (R + a ^ 2) / R :=
    div_nonneg (mul_nonneg (mul_nonneg zero_le_two ha.le) (Real.sqrt_nonneg _)) hR.le
  have hspos : 0 < -w.re := by linarith
  have hnorm : -w.re ≤ ‖w‖ := by
    have := Complex.abs_re_le_norm w
    rw [abs_of_neg (by linarith : w.re < 0)] at this
    exact this
  have hδR : 2 * a * Real.sqrt (R + a ^ 2) / R * R = 2 * a * Real.sqrt (R + a ^ 2) :=
    div_mul_cancel₀ _ hR.ne'
  have e : (2 * a * Real.sqrt (R + a ^ 2) / R * -w.re) ^ 2 * R ^ 2
      = 4 * a ^ 2 * (R + a ^ 2) * (-w.re) ^ 2 := by
    calc (2 * a * Real.sqrt (R + a ^ 2) / R * -w.re) ^ 2 * R ^ 2
        = (2 * a * Real.sqrt (R + a ^ 2) / R * R) ^ 2 * (-w.re) ^ 2 := by ring
      _ = 4 * a ^ 2 * Real.sqrt (R + a ^ 2) ^ 2 * (-w.re) ^ 2 := by rw [hδR]; ring
      _ = 4 * a ^ 2 * (R + a ^ 2) * (-w.re) ^ 2 := by rw [hsq]
  have h3 : (a ^ 2 + -w.re) * R ^ 2 ≤ (R + a ^ 2) * (-w.re) ^ 2 := by
    have hpos : 0 ≤ a ^ 2 * (-w.re + R) + R * -w.re :=
      add_nonneg (mul_nonneg (sq_nonneg a) (by linarith)) (mul_nonneg hR.le hspos.le)
    nlinarith [mul_nonneg (by linarith : (0:ℝ) ≤ -w.re - R) hpos]
  have h4 : w.im ^ 2 * R ^ 2 ≤ (2 * a * Real.sqrt (R + a ^ 2) / R * -w.re) ^ 2 * R ^ 2 := by
    rw [e]
    have h5 := mul_le_mul_of_nonneg_left h3 (by positivity : (0:ℝ) ≤ 4 * a ^ 2)
    have h6 := mul_le_mul_of_nonneg_right hw (sq_nonneg R)
    nlinarith
  have h2 : w.im ^ 2 ≤ (2 * a * Real.sqrt (R + a ^ 2) / R * -w.re) ^ 2 :=
    le_of_mul_le_mul_right h4 (pow_pos hR 2)
  calc |w.im| ≤ Real.sqrt ((2 * a * Real.sqrt (R + a ^ 2) / R * -w.re) ^ 2) := Real.abs_le_sqrt h2
    _ = 2 * a * Real.sqrt (R + a ^ 2) / R * -w.re := Real.sqrt_sq (mul_nonneg hδ hspos.le)
    _ ≤ 2 * a * Real.sqrt (R + a ^ 2) / R * ‖w‖ := mul_le_mul_of_nonneg_left hnorm hδ

lemma sectorDelta_pos {a R : ℝ} (ha : 0 < a) (hR : 0 < R) : 0 < sectorDelta a R :=
  div_pos (mul_pos (mul_pos two_pos ha) (Real.sqrt_pos.mpr (add_pos hR (pow_pos ha 2)))) hR

/-- D(R) δ(R)^2 = 1. -/
lemma Dfun_mul_sq {a R : ℝ} (ha : 0 < a) (hR : 0 < R) : Dfun a R * sectorDelta a R ^ 2 = 1 := by
  unfold Dfun sectorDelta
  have hX : 0 < R + a ^ 2 := add_pos hR (pow_pos ha 2)
  have hsq : Real.sqrt (R + a ^ 2) ^ 2 = R + a ^ 2 := Real.sq_sqrt hX.le
  rw [div_pow, mul_pow, mul_pow, hsq]
  field_simp
  norm_num

end

end KimEff

