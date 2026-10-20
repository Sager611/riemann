import Riemann.Deep.FirstZero.EM
import Riemann.Fayed.Defs

/-!
# Fayed.Inputs.em17 (Euler–Maclaurin, k = 4) from Mathlib + the order-8 formula zeta_EM8

At N = ⌊u⌋: zeta_EM8 gives S(N) = ζ(s) + N^{1-s}/(1-s) - EMcorr(N) + O(N^{-σ-3}).
On [N, u] the function hE(x) = x^{1-s}/(1-s) - EMcorr(x) has derivative
s(s+1)(s+2) B₃(x-N)/6 · x^{-s-3}, so moving from N to u costs O(N^{-σ-3}) (mean value).
-/

open Complex Set Filter Topology MeasureTheory

namespace Fayed.InputsProof

open Fayed.FirstZero

noncomputable section

/-- Smooth model of x^{1-s}/(1-s) - EMcorr x on [N, N+1). -/
def hE (s : ℂ) (N : ℝ) (x : ℝ) : ℂ :=
  (x : ℂ) ^ (1 - s) / (1 - s) - ((BP 1 (x - N) : ℂ) * (x : ℂ) ^ (-s) +
    s * (BP 2 (x - N) : ℂ) * (x : ℂ) ^ (-s - 1) +
    s * (s + 1) * (BP 3 (x - N) : ℂ) * (x : ℂ) ^ (-s - 2))

lemma cpow_hd (c : ℂ) {x : ℝ} (hx : 0 < x) :
    HasDerivAt (fun y : ℝ => (y : ℂ) ^ c) (c * (x : ℂ) ^ (c - 1)) x :=
  ((Complex.hasStrictDerivAt_cpow_const (c := c) (ofReal_mem_slitPlane.mpr hx)).hasDerivAt).comp_ofReal

lemma bp_hd {k : ℕ} (hk1 : 1 ≤ k) (hk : k ≤ 7) (N x : ℝ) :
    HasDerivAt (fun y : ℝ => (BP (k + 1) (y - N) : ℂ)) (BP k (x - N) : ℂ) x := by
  have := ((BP_hasDerivAt hk1 hk (x - N)).comp x ((hasDerivAt_id' x).sub_const N)).ofReal_comp
  simpa using this

lemma bp1_hd (N x : ℝ) : HasDerivAt (fun y : ℝ => (BP 1 (y - N) : ℂ)) 1 x := by
  have := (((hasDerivAt_id' x).sub_const N).sub_const (1 / 2 : ℝ)).ofReal_comp
  simpa [BP] using this

lemma hE_hasDerivAt {s : ℂ} (hs : s ≠ 1) (N : ℝ) {x : ℝ} (hx : 0 < x) :
    HasDerivAt (hE s N) (s * (s + 1) * (s + 2) * (BP 3 (x - N) : ℂ) * (x : ℂ) ^ (-s - 3)) x := by
  have hs1 : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm hs)
  have h1 := (cpow_hd (1 - s) hx).div_const (1 - s)
  have a1 := (bp1_hd N x).mul (cpow_hd (-s) hx)
  have a2 := ((bp_hd (k := 1) le_rfl (by norm_num) N x).const_mul s).mul (cpow_hd (-s - 1) hx)
  have a3 := ((bp_hd (k := 2) (by norm_num) (by norm_num) N x).const_mul (s * (s + 1))).mul
    (cpow_hd (-s - 2) hx)
  have := h1.sub ((a1.add a2).add a3)
  convert this using 1
  · funext y; rfl
  · rw [show (1 : ℂ) - s - 1 = -s by ring, show -s - 1 - 1 = -s - 2 by ring,
      show -s - 2 - 1 = -s - 3 by ring]
    field_simp
    ring

lemma abs_BP3_le {t : ℝ} (h0 : 0 ≤ t) (h1 : t ≤ 1) : |BP 3 t| ≤ 1 := by
  simp only [BP]
  have t2 : t ^ 2 ≤ 1 := by nlinarith
  have t3 : t ^ 3 ≤ 1 := by nlinarith
  rw [abs_le]; constructor <;> nlinarith [sq_nonneg t, pow_nonneg h0 3]

lemma Bbar2_eq (x : ℝ) : Bbar 2 x = Int.fract x ^ 2 - Int.fract x + 1 / 6 := by
  simp [Bbar, Polynomial.bernoulli, Finset.sum_range_succ, bernoulli, bernoulli'_two,
    bernoulli'_one, Nat.choose]
  ring

lemma Bbar3_eq (x : ℝ) : Bbar 3 x = Int.fract x ^ 3 - 3 / 2 * Int.fract x ^ 2 + 1 / 2 * Int.fract x := by
  simp [Bbar, Polynomial.bernoulli, Finset.sum_range_succ, bernoulli, bernoulli'_two,
    bernoulli'_one, bernoulli'_three, Nat.choose]
  ring

lemma Bbar1_eq (x : ℝ) : Bbar 1 x = Int.fract x - 1 / 2 := by
  simp [Bbar, Polynomial.bernoulli_one]

lemma cpow_neg_eq_one_div {x : ℝ} (c : ℂ) : (x : ℂ) ^ (-c) = 1 / (x : ℂ) ^ c := by
  rw [cpow_neg, one_div]

/-- EMcorr agrees with the smooth model hE (u ≥ 0, N = ⌊u⌋). -/
lemma EMcorr_eq {u : ℝ} (hu : 0 ≤ u) (s : ℂ) :
    (u : ℂ) ^ (1 - s) / (1 - s) - EMcorr u s = hE s (⌊u⌋₊ : ℝ) u := by
  have hf : Int.fract u = u - (⌊u⌋₊ : ℝ) := by
    rw [Int.fract, natCast_floor_eq_intCast_floor hu]
  have e2 : (-s - 1) = -(s + 1) := by ring
  have e3 : (-s - 2) = -(s + 2) := by ring
  simp only [hE, EMcorr, e2, e3, cpow_neg_eq_one_div, Bbar1_eq, Bbar2_eq, Bbar3_eq, hf, BP]
  push_cast
  ring

lemma S_eq_sum {s : ℂ} (hs0 : s ≠ 0) (u : ℝ) :
    S u s = ∑ n ∈ Finset.range (⌊u⌋₊ + 1), 1 / (n : ℂ) ^ s := by
  have hr : Finset.range (⌊u⌋₊ + 1) = insert 0 (Finset.Icc 1 ⌊u⌋₊) := by
    ext n; simp; omega
  rw [hr, Finset.sum_insert (by simp), S]
  simp only [Nat.cast_zero, zero_cpow hs0, div_zero, zero_add]
  refine Finset.sum_congr rfl fun n _ => ?_
  rw [cpow_neg, one_div]

lemma rpow_mono_neg {N x a : ℝ} (hN : 0 < N) (hx : N ≤ x) (ha : 0 ≤ a) :
    x ^ (-a) ≤ N ^ (-a) := Real.rpow_le_rpow_of_nonpos hN hx (by linarith)

theorem em17_holds : ∀ s : ℂ, 0 < s.re → s.re < 1 → ∃ C : ℝ, ∀ u : ℝ, 1 ≤ u →
    ‖S u s - (riemannZeta s + (u : ℂ) ^ (1 - s) / (1 - s) - EMcorr u s)‖ ≤
      C * u ^ (-(s.re + 3)) := by
  intro s hσ0 hσ1
  have hs1 : s ≠ 1 := by intro h; rw [h] at hσ1; simp at hσ1
  have hs0 : s ≠ 0 := by intro h; rw [h] at hσ0; simp at hσ0
  set σ := s.re with hσ
  set K : ℝ := ‖s‖ * ‖gc s 7‖ + ‖s‖ * (‖gc s 2‖ + ‖gc s 4‖ + ‖gc s 6‖) +
    ‖s * (s + 1) * (s + 2)‖ with hK
  have hK0 : 0 ≤ K := by positivity
  refine ⟨K * 2 ^ (σ + 3), fun u hu => ?_⟩
  set N : ℕ := ⌊u⌋₊ with hNdef
  have hN1 : 1 ≤ N := Nat.le_floor (by exact_mod_cast hu)
  have hNr : (1 : ℝ) ≤ N := by exact_mod_cast hN1
  have hN0 : (0 : ℝ) < N := by linarith
  have hNu : (N : ℝ) ≤ u := Nat.floor_le (by linarith)
  have huN : u < N + 1 := Nat.lt_floor_add_one u
  have hσ3 : 0 ≤ σ + 3 := by linarith
  -- powers
  have pw : ∀ j : ℕ, 2 ≤ j → (N : ℝ) ^ (-σ - 1 - j) ≤ (N : ℝ) ^ (-(σ + 3)) := by
    intro j hj
    apply Real.rpow_le_rpow_of_exponent_le hNr
    have : (2 : ℝ) ≤ j := by exact_mod_cast hj
    linarith
  -- (A) order-8 Euler–Maclaurin at N
  have hA := zeta_EM8 hσ0 hs1 hN1
  rw [← S_eq_sum hs0 u] at hA
  set eA := riemannZeta s - (S u s + (-(N : ℂ) ^ (1 - s)) / (1 - s) + (-(N : ℂ) ^ (-s)) / 2 +
    s * EMPhi s N) with heA
  have hA' : ‖eA‖ ≤ ‖s‖ * ‖gc s 7‖ * (N : ℝ) ^ (-(σ + 3)) := by
    refine hA.trans ?_
    rw [mul_assoc]
    apply mul_le_mul_of_nonneg_left _ (norm_nonneg _)
    have h7 : (N : ℝ) ^ (-(σ + 7)) ≤ (N : ℝ) ^ (-(σ + 3)) :=
      Real.rpow_le_rpow_of_exponent_le hNr (by linarith)
    have hp : 0 ≤ (N : ℝ) ^ (-(σ + 7)) := by positivity
    have : (N : ℝ) ^ (-(σ + 7)) / (σ + 7) ≤ (N : ℝ) ^ (-(σ + 7)) :=
      div_le_self hp (by linarith)
    have hg : 0 ≤ ‖gc s 7‖ := norm_nonneg _
    have hd : ‖gc s 7‖ / 1000000 ≤ ‖gc s 7‖ := div_le_self hg (by norm_num)
    exact mul_le_mul hd (this.trans h7) (by positivity) hg
  -- (B) EMPhi minus its leading term
  set eB := EMPhi s N - 1 / 12 * (N : ℂ) ^ (-s - 1) with heB
  have hB : ‖eB‖ ≤ (‖gc s 2‖ + ‖gc s 4‖ + ‖gc s 6‖) * (N : ℝ) ^ (-(σ + 3)) := by
    have e0 : gf s 0 N = (N : ℂ) ^ (-s - 1) := by simp [gf, gc]
    have : eB = -(1 / 720) * gf s 2 N + 1 / 30240 * gf s 4 N - 1 / 1209600 * gf s 6 N := by
      rw [heB, EMPhi, e0]; ring
    rw [this]
    have n2 := norm_gf (s := s) (j := 2) hN0
    have n4 := norm_gf (s := s) (j := 4) hN0
    have n6 := norm_gf (s := s) (j := 6) hN0
    have p2 := pw 2 (by norm_num)
    have p4 := pw 4 (by norm_num)
    have p6 := pw 6 (by norm_num)
    have g2 : ‖gf s 2 N‖ ≤ ‖gc s 2‖ * (N : ℝ) ^ (-(σ + 3)) := by
      rw [n2]; exact mul_le_mul_of_nonneg_left (by exact_mod_cast p2) (norm_nonneg _)
    have g4 : ‖gf s 4 N‖ ≤ ‖gc s 4‖ * (N : ℝ) ^ (-(σ + 3)) := by
      rw [n4]; exact mul_le_mul_of_nonneg_left (by exact_mod_cast p4) (norm_nonneg _)
    have g6 : ‖gf s 6 N‖ ≤ ‖gc s 6‖ * (N : ℝ) ^ (-(σ + 3)) := by
      rw [n6]; exact mul_le_mul_of_nonneg_left (by exact_mod_cast p6) (norm_nonneg _)
    calc ‖-(1 / 720) * gf s 2 N + 1 / 30240 * gf s 4 N - 1 / 1209600 * gf s 6 N‖
        ≤ ‖gf s 2 N‖ + ‖gf s 4 N‖ + ‖gf s 6 N‖ := by
          refine (norm_sub_le _ _).trans ?_
          refine add_le_add ((norm_add_le _ _).trans (add_le_add ?_ ?_)) ?_ <;>
            (rw [norm_mul]; apply mul_le_of_le_one_left (norm_nonneg _); norm_num)
      _ ≤ _ := by linarith
  -- (C) mean value on [N, u]
  have hC : ‖hE s N u - hE s N N‖ ≤ ‖s * (s + 1) * (s + 2)‖ * (N : ℝ) ^ (-(σ + 3)) := by
    have hmv := norm_image_sub_le_of_norm_deriv_le_segment' (f := hE s N) (a := N) (b := u)
      (C := ‖s * (s + 1) * (s + 2)‖ * (N : ℝ) ^ (-(σ + 3)))
      (f' := fun x => s * (s + 1) * (s + 2) * (BP 3 (x - N) : ℂ) * (x : ℂ) ^ (-s - 3))
      (fun x hx => (hE_hasDerivAt hs1 N (by linarith [hx.1])).hasDerivWithinAt)
      (fun x hx => by
        have hx0 : 0 < x := by linarith [hx.1]
        rw [norm_mul, norm_mul, norm_cpow_eq_rpow_re_of_pos hx0]
        have hb : ‖(BP 3 (x - N) : ℂ)‖ ≤ 1 := by
          rw [Complex.norm_real, Real.norm_eq_abs]
          exact abs_BP3_le (by linarith [hx.1]) (by linarith [hx.2])
        have hr : x ^ (-s - 3).re ≤ (N : ℝ) ^ (-(σ + 3)) := by
          rw [show (-s - 3).re = -(σ + 3) by simp [hσ]; ring]
          exact rpow_mono_neg hN0 hx.1 hσ3
        calc ‖s * (s + 1) * (s + 2)‖ * ‖(BP 3 (x - N) : ℂ)‖ * x ^ (-s - 3).re
            ≤ ‖s * (s + 1) * (s + 2)‖ * 1 * (N : ℝ) ^ (-(σ + 3)) := by gcongr
          _ = _ := by ring)
      u ⟨hNu, le_rfl⟩
    refine hmv.trans ?_
    have : u - N ≤ 1 := by linarith
    have hp : 0 ≤ ‖s * (s + 1) * (s + 2)‖ * (N : ℝ) ^ (-(σ + 3)) := by positivity
    calc ‖s * (s + 1) * (s + 2)‖ * (N : ℝ) ^ (-(σ + 3)) * (u - N)
        ≤ ‖s * (s + 1) * (s + 2)‖ * (N : ℝ) ^ (-(σ + 3)) * 1 :=
          mul_le_mul_of_nonneg_left this hp
      _ = _ := by ring
  -- hE at N
  have hEN : hE s N N = (N : ℂ) ^ (1 - s) / (1 - s) + (N : ℂ) ^ (-s) / 2 -
      s / 12 * (N : ℂ) ^ (-s - 1) := by
    simp only [hE, sub_self, BP]; push_cast; ring
  -- assemble
  have key : S u s - (riemannZeta s + (u : ℂ) ^ (1 - s) / (1 - s) - EMcorr u s) =
      -eA - s * eB - (hE s N u - hE s N N) := by
    have h1 : riemannZeta s + (u : ℂ) ^ (1 - s) / (1 - s) - EMcorr u s =
        riemannZeta s + hE s N u := by
      rw [← EMcorr_eq (by linarith) s]; ring
    rw [h1, heA, heB, hEN]; ring
  rw [key]
  have hpow : (N : ℝ) ^ (-(σ + 3)) ≤ 2 ^ (σ + 3) * u ^ (-(σ + 3)) := by
    have h2N : u ≤ 2 * N := by linarith
    have hm := rpow_mono_neg (by positivity : (0 : ℝ) < u) h2N hσ3
    rw [Real.mul_rpow (by norm_num) hN0.le] at hm
    have h22 : (2 : ℝ) ^ (σ + 3) * 2 ^ (-(σ + 3)) = 1 := by
      rw [Real.rpow_neg (by norm_num), mul_inv_cancel₀ (by positivity)]
    have hpos : 0 ≤ (2 : ℝ) ^ (σ + 3) := by positivity
    calc (N : ℝ) ^ (-(σ + 3)) = 2 ^ (σ + 3) * (2 ^ (-(σ + 3)) * (N : ℝ) ^ (-(σ + 3))) := by
          rw [← mul_assoc, h22, one_mul]
      _ ≤ 2 ^ (σ + 3) * u ^ (-(σ + 3)) := mul_le_mul_of_nonneg_left hm hpos
  have hsB : ‖s * eB‖ ≤ ‖s‖ * ((‖gc s 2‖ + ‖gc s 4‖ + ‖gc s 6‖) * (N : ℝ) ^ (-(σ + 3))) := by
    rw [norm_mul]; exact mul_le_mul_of_nonneg_left hB (norm_nonneg _)
  calc ‖-eA - s * eB - (hE s N u - hE s N N)‖
      ≤ ‖eA‖ + ‖s * eB‖ + ‖hE s N u - hE s N N‖ := by
        refine (norm_sub_le _ _).trans (add_le_add ((norm_sub_le _ _).trans ?_) le_rfl)
        rw [norm_neg]
    _ ≤ K * (N : ℝ) ^ (-(σ + 3)) := by rw [hK]; nlinarith [hA', hsB, hC]
    _ ≤ K * (2 ^ (σ + 3) * u ^ (-(σ + 3))) := mul_le_mul_of_nonneg_left hpow hK0
    _ = K * 2 ^ (σ + 3) * u ^ (-(σ + 3)) := by ring

end

end Fayed.InputsProof
