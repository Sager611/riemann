import Riemann.Deep.DPCBeyond.Lean.Rpow
import Mathlib.Algebra.Order.BigOperators.Ring.Finset

/-!
# dpcbeyond Lemma Z (new step (b)): the mechanism and the summation, with abstract exp/log

PROOF.md Lemma Z (iii)–(iv). Proved here, from no hypotheses beyond ExpLog:
* sum_mul_le, prod_ex_le(_sup) — (iii), the length constraint kept inside the r-fold integral: u_i ≥ 0,
  Σu_i ≤ log(X/𝔞) and Re z_i ≤ D (D = max_i (Re z_i)⁺ ≥ 0) give |e^{Σ z_i u_i}| = ∏ e^{(Re z_i)u_i} ≤ (X/𝔞)^D,
  one factor X^{max d} for the whole r-fold product (not X^{r·max d});
* rp_sup_le_sum — (iv): X^{max_i c_i} ≤ Σ_i X^{c_i};
* tuple_sum, tuple_bound — (iv): Σ_{ρ_1..ρ_m} X^{max_i c(ρ_i)} ∏_i k(ρ_i) ≤ m · W · On^{m−1}, W = Σ_ρ X^{c(ρ)}k(ρ),
  On = Σ_ρ k(ρ), over any finite set of zeros (the full sums are the suprema of these, all terms ≥ 0);
* packet_case1, packet_case2, packet_nofree, pos_part_shift, pos_part_split — (iv)'s packet exponents
  (𝔞^{−1/2∓y}, (X/𝔞)^{…}, (d + y)⁺ ≤ y + d⁺) and 6e′ (a)'s split W_{−y} ≤ On + X^{−y}Off;
* term_bound — (iii)+(iv) for one partition/binomial term with m ≥ 1 free variables, from the NAMED integral
  bound |F(𝔞)| ≤ C_I Σ_{ρ⃗} (X/𝔞)^{max c} ∏ k (Lemma Z (i)–(iii): explicit formula, Fourier inversion, Fubini, IBP)
  and the NAMED packet sum Σ_𝔞 C(𝔞)𝔞^{−1/2} ≤ K_pk (Mertens, (iv)).
Not formalized: the explicit formula (i), the Möbius inversion over set partitions and the binomial expansion (ii),
the integration by parts (iii) (analysis); they enter only through the named bounds of term_bound.
-/

namespace DPCBLean

open Finset DPCLean

variable {ex lg : ℝ → ℝ}

/-- ex turns finite sums into products. -/
lemma ex_sum (h : ExpLog ex lg) {ι : Type*} (s : Finset ι) (x : ι → ℝ) :
    ex (∑ i ∈ s, x i) = ∏ i ∈ s, ex (x i) := by
  classical
  refine Finset.induction_on s ?_ ?_
  · simp [h.zero]
  · intro i s hi ih
    rw [sum_insert hi, prod_insert hi, h.add, ih]

/-- (iii), the length constraint: u_i ≥ 0, Σu_i ≤ ℓ, d_i ≤ D, D ≥ 0 give Σ d_i u_i ≤ D ℓ. -/
lemma sum_mul_le {ι : Type*} (s : Finset ι) (d u : ι → ℝ) {D ℓ : ℝ} (hD : 0 ≤ D)
    (hu : ∀ i ∈ s, 0 ≤ u i) (hd : ∀ i ∈ s, d i ≤ D) (hl : ∑ i ∈ s, u i ≤ ℓ) :
    ∑ i ∈ s, d i * u i ≤ D * ℓ := by
  calc ∑ i ∈ s, d i * u i ≤ ∑ i ∈ s, D * u i :=
        sum_le_sum (fun i hi => mul_le_mul_of_nonneg_right (hd i hi) (hu i hi))
    _ = D * ∑ i ∈ s, u i := by rw [mul_sum]
    _ ≤ D * ℓ := mul_le_mul_of_nonneg_left hl hD

/-- **Lemma Z (iii), the decisive bound**: |e^{Σ z_i u_i}| = ∏ e^{(Re z_i)u_i} ≤ (X/𝔞)^D on the support
(u_i ≥ 0, Σu_i ≤ log(X/𝔞) = lg Y), for any D ≥ 0 with Re z_i ≤ D. -/
theorem prod_ex_le (h : ExpLog ex lg) {ι : Type*} (s : Finset ι) (d u : ι → ℝ) {D Y : ℝ}
    (hY : 0 < Y) (hD : 0 ≤ D) (hu : ∀ i ∈ s, 0 ≤ u i) (hd : ∀ i ∈ s, d i ≤ D)
    (hl : ∑ i ∈ s, u i ≤ lg Y) : ∏ i ∈ s, ex (d i * u i) ≤ rp ex lg Y D := by
  rw [← ex_sum h, rp_of_pos hY]
  exact h.mono.monotone (sum_mul_le s d u hD hu hd hl)

/-- (iii) with D = max_i (Re z_i)⁺. -/
theorem prod_ex_le_sup (h : ExpLog ex lg) {ι : Type*} (s : Finset ι) (hs : s.Nonempty)
    (d u : ι → ℝ) {Y : ℝ} (hY : 0 < Y) (hu : ∀ i ∈ s, 0 ≤ u i) (hl : ∑ i ∈ s, u i ≤ lg Y) :
    ∏ i ∈ s, ex (d i * u i) ≤ rp ex lg Y (s.sup' hs (fun i => max (d i) 0)) := by
  obtain ⟨i₀, hi₀⟩ := hs
  refine prod_ex_le h s d u hY ?_ hu (fun i hi => ?_) hl
  · exact le_trans (le_max_right (d i₀) 0) (le_sup' (fun i => max (d i) 0) hi₀)
  · exact le_trans (le_max_left (d i) 0) (le_sup' (fun i => max (d i) 0) hi)

/-- (iv): X^{max_i c_i} ≤ Σ_i X^{c_i}. -/
lemma rp_sup_le_sum (h : ExpLog ex lg) {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (c : ι → ℝ)
    (X : ℝ) : rp ex lg X (s.sup' hs c) ≤ ∑ i ∈ s, rp ex lg X (c i) := by
  obtain ⟨i, hi, he⟩ := exists_mem_eq_sup' hs c
  rw [he]
  exact single_le_sum (f := fun i => rp ex lg X (c i)) (fun j _ => rp_nonneg h _ _) hi

/-- a(ρ_j) ∏_i k(ρ_i) as one product. -/
lemma mul_prod_eq {ι : Type*} (n : ℕ) (a k : ι → ℝ) (j : Fin (n + 1)) (ρ : Fin (n + 1) → ι) :
    a (ρ j) * ∏ i, k (ρ i) = ∏ i, (if i = j then a (ρ i) * k (ρ i) else k (ρ i)) := by
  rw [← mul_prod_erase univ (fun i => k (ρ i)) (mem_univ j),
    ← mul_prod_erase univ (fun i => if i = j then a (ρ i) * k (ρ i) else k (ρ i)) (mem_univ j)]
  rw [ite_eq_left rfl, prod_congr rfl (fun i hi => ite_eq_right (ne_of_mem_erase hi))]
  ring

/-- Σ_{ρ ∈ S^{n+1}} ∏_i (if i = j then a k else k)(ρ_i) = (Σ a k)(Σ k)^n. -/
lemma tuple_one {ι : Type*} [DecidableEq ι] (S : Finset ι) (a k : ι → ℝ) (n : ℕ) (j : Fin (n + 1)) :
    ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
        ∏ i, (if i = j then a (ρ i) * k (ρ i) else k (ρ i)) =
      (∑ x ∈ S, a x * k x) * (∑ x ∈ S, k x) ^ n := by
  rw [sum_prod_piFinset S (fun i x => if i = j then a x * k x else k x)]
  rw [← mul_prod_erase univ (fun i => ∑ x ∈ S, if i = j then a x * k x else k x) (mem_univ j)]
  simp only
  congr 1
  rw [prod_congr rfl (fun i hi => sum_congr rfl (fun x _ => ite_eq_right (ne_of_mem_erase hi))),
    prod_const, card_erase_of_mem (mem_univ j), card_univ, Fintype.card_fin]
  rfl

/-- (iv), the exchange: Σ_{ρ ∈ S^m} (Σ_j a(ρ_j)) ∏_i k(ρ_i) = m (Σ a k)(Σ k)^{m−1}, m = n + 1. -/
theorem tuple_sum {ι : Type*} [DecidableEq ι] (S : Finset ι) (a k : ι → ℝ) (n : ℕ) :
    ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S), (∑ j, a (ρ j)) * ∏ i, k (ρ i) =
      (n + 1) * (∑ x ∈ S, a x * k x) * (∑ x ∈ S, k x) ^ n := by
  simp_rw [sum_mul]
  rw [sum_comm]
  rw [sum_congr rfl (fun j _ => by simp_rw [mul_prod_eq n a k j]; exact tuple_one S a k n j)]
  rw [sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul]
  push_cast
  ring

/-- **Lemma Z (iv)**: Σ_{ρ_1..ρ_m} X^{max_i c(ρ_i)} ∏_i k(ρ_i) ≤ m · W · On^{m−1} (m = n + 1 ≥ 1, k ≥ 0),
W = Σ_ρ X^{c(ρ)} k(ρ), On = Σ_ρ k(ρ), over a finite set S of zeros. -/
theorem tuple_bound (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] (S : Finset ι) (c k : ι → ℝ)
    (hk : ∀ x, 0 ≤ k x) (X : ℝ) (n : ℕ) :
    ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
        rp ex lg X (univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i) ≤
      (n + 1) * (∑ x ∈ S, rp ex lg X (c x) * k x) * (∑ x ∈ S, k x) ^ n := by
  rw [← tuple_sum S (fun x => rp ex lg X (c x)) k n]
  exact sum_le_sum (fun ρ _ => mul_le_mul_of_nonneg_right
    (rp_sup_le_sum h univ univ_nonempty (fun i => c (ρ i)) X) (prod_nonneg (fun i _ => hk _)))

/-- (iv) case 1 (s = 1/2 + y + it): 𝔞^{−1/2−y}(X/𝔞)^D ≤ 𝔞^{−1/2}X^D for 𝔞 ≥ 1, y, D ≥ 0. -/
lemma packet_case1 (h : ExpLog ex lg) {a X y D : ℝ} (ha : 1 ≤ a) (hX : 0 < X) (hy : 0 ≤ y)
    (hD : 0 ≤ D) : rp ex lg a (-(1 / 2) - y) * rp ex lg (X / a) D ≤
      rp ex lg a (-(1 / 2)) * rp ex lg X D := by
  have ha0 : 0 < a := by linarith
  have h1 : X / a ≤ X := div_le_self hX.le ha
  exact mul_le_mul (rp_le_rp_exp h ha (by linarith)) (rp_le_rp h h1 hD) (rp_nonneg h _ _)
    (rp_nonneg h _ _)

/-- (iv) case 2 (s = 1/2 − y + it): 𝔞^{−1/2+y}(X/𝔞)^{y+D} = X^{y+D}𝔞^{−1/2−D} ≤ X^{y+D}𝔞^{−1/2} (𝔞 ≥ 1, D ≥ 0). -/
lemma packet_case2 (h : ExpLog ex lg) {a X y D : ℝ} (ha : 1 ≤ a) (hX : 0 < X) (hD : 0 ≤ D) :
    rp ex lg a (-(1 / 2) + y) * rp ex lg (X / a) (y + D) ≤
      rp ex lg X (y + D) * rp ex lg a (-(1 / 2)) := by
  have ha0 : 0 < a := by linarith
  have hXa : 0 < X / a := div_pos hX ha0
  have e0 := rp_mul h hXa ha0 (y + D)
  rw [div_mul_cancel₀ X ha0.ne'] at e0
  have hn := rp_neg h ha0 (y + D)
  have e1 : rp ex lg (X / a) (y + D) = rp ex lg X (y + D) * rp ex lg a (-(y + D)) := by
    calc rp ex lg (X / a) (y + D) =
          rp ex lg (X / a) (y + D) * (rp ex lg a (-(y + D)) * rp ex lg a (y + D)) := by
          rw [hn, mul_one]
      _ = (rp ex lg (X / a) (y + D) * rp ex lg a (y + D)) * rp ex lg a (-(y + D)) := by ring
      _ = rp ex lg X (y + D) * rp ex lg a (-(y + D)) := by rw [← e0]
  rw [e1]
  calc rp ex lg a (-(1 / 2) + y) * (rp ex lg X (y + D) * rp ex lg a (-(y + D))) =
        rp ex lg X (y + D) * (rp ex lg a (-(1 / 2) + y) * rp ex lg a (-(y + D))) := by ring
    _ = rp ex lg X (y + D) * rp ex lg a (-(1 / 2) - D) := by
        rw [← rp_add h ha0, show -(1 / 2) + y + -(y + D) = -(1 / 2) - D by ring]
    _ ≤ rp ex lg X (y + D) * rp ex lg a (-(1 / 2)) :=
        mul_le_mul_of_nonneg_left (rp_le_rp_exp h ha (by linarith)) (rp_nonneg h _ _)

/-- (iv), no free variable left: 1 ≤ 𝔞 ≤ X gives 𝔞^{−1/2+y} ≤ X^y 𝔞^{−1/2} (y ≥ 0). -/
lemma packet_nofree (h : ExpLog ex lg) {a X y : ℝ} (ha : 1 ≤ a) (haX : a ≤ X) (hy : 0 ≤ y) :
    rp ex lg a (-(1 / 2) + y) ≤ rp ex lg X y * rp ex lg a (-(1 / 2)) := by
  have ha0 : 0 < a := by linarith
  rw [rp_add h ha0, mul_comm]
  exact mul_le_mul_of_nonneg_right (rp_le_rp h haX hy) (rp_nonneg h _ _)

/-- (iv) case 2's weight: (d + y)⁺ ≤ y + d⁺, so X^{(d+y)⁺} ≤ X^y X^{d⁺} (X ≥ 1, y ≥ 0): W_− ≤ X^y W. -/
lemma pos_part_shift (h : ExpLog ex lg) {X d y : ℝ} (hX : 1 ≤ X) (hy : 0 ≤ y) :
    rp ex lg X (max (d + y) 0) ≤ rp ex lg X y * rp ex lg X (max d 0) := by
  rw [← rp_add h (by linarith)]
  apply rp_le_rp_exp h hX
  rcases le_total d 0 with hd | hd
  · rw [max_eq_right hd]
    exact max_le (by linarith) (by linarith)
  · rw [max_eq_left hd]
    exact max_le (by linarith) (by linarith)

/-- 6e′ (a): a zero with d ≤ y contributes X^0 = 1 to W_{−y}, one with d > y contributes X^{−y}X^d:
X^{(d−y)⁺} ≤ 1 + [d > 0] X^{−y}X^d (X ≥ 1, y ≥ 0); summed: W_{−y} ≤ On + X^{−y}Off. -/
lemma pos_part_split (h : ExpLog ex lg) {X d y : ℝ} (hX : 1 ≤ X) (hy : 0 ≤ y) :
    rp ex lg X (max (d - y) 0) ≤ 1 + (if 0 < d then rp ex lg X (-y) * rp ex lg X d else 0) := by
  have hX0 : 0 < X := by linarith
  rcases le_or_gt d y with hd | hd
  · rw [max_eq_right (by linarith), rp_zero h hX0]
    have h1 : 0 ≤ (if 0 < d then rp ex lg X (-y) * rp ex lg X d else 0) := by
      split_ifs
      · exact mul_nonneg (rp_nonneg h _ _) (rp_nonneg h _ _)
      · exact le_refl 0
    linarith
  · have hd0 : 0 < d := lt_of_le_of_lt hy hd
    rw [max_eq_left (by linarith), ite_eq_left hd0, sub_eq_neg_add, rp_add h hX0]
    have h1 := mul_nonneg (rp_nonneg h X (-y)) (rp_nonneg h X d)
    linarith

/-- **Lemma Z (iii)+(iv) for one term** (case 1, s = 1/2 + y + it) with m = n + 1 free variables: packets 𝔞 ∈ P
with value v 𝔞 ≥ 1 and coefficient C 𝔞 ≥ 0 (stripped of 𝔞^{−s}); NAMED: the free sum after the explicit formula
obeys |F 𝔞| ≤ C_I Σ_{ρ⃗ ∈ S^m} (X/v 𝔞)^{max_i c(ρ_i)} ∏ k(ρ_i) with c ≥ 0 (Lemma Z (i)–(iii), c = (d − y)⁺, k = κ_A),
and Σ_𝔞 C 𝔞 (v 𝔞)^{−1/2} ≤ K_pk (Mertens). Then Σ_𝔞 C 𝔞 (v 𝔞)^{−1/2−y}|F 𝔞| ≤ K_pk C_I m W On^{m−1}. -/
theorem term_bound (h : ExpLog ex lg) {ι α : Type*} [DecidableEq ι] (S : Finset ι) (c k : ι → ℝ)
    (hc : ∀ x, 0 ≤ c x) (hk : ∀ x, 0 ≤ k x) {X y CI Kpk : ℝ} (hX : 0 < X) (hy : 0 ≤ y)
    (hCI : 0 ≤ CI) (n : ℕ) (P : Finset α) (v C F : α → ℝ) (hv : ∀ p ∈ P, 1 ≤ v p)
    (hC : ∀ p ∈ P, 0 ≤ C p)
    (hF : ∀ p ∈ P, |F p| ≤ CI * ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
      rp ex lg (X / v p) (univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i))
    (hP : ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) ≤ Kpk) :
    ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) - y) * |F p| ≤
      Kpk * CI * ((n + 1) * (∑ x ∈ S, rp ex lg X (c x) * k x) * (∑ x ∈ S, k x) ^ n) := by
  set G := (n + 1) * (∑ x ∈ S, rp ex lg X (c x) * k x) * (∑ x ∈ S, k x) ^ n with hG
  have hG0 : 0 ≤ G := by
    rw [hG]
    exact mul_nonneg (mul_nonneg (by positivity) (sum_nonneg (fun x _ =>
      mul_nonneg (rp_nonneg h _ _) (hk x)))) (pow_nonneg (sum_nonneg (fun x _ => hk x)) _)
  have hterm : ∀ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) - y) * |F p| ≤
      C p * rp ex lg (v p) (-(1 / 2)) * (CI * G) := by
    intro p hp
    have hvp := hv p hp
    have hsum : ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
        rp ex lg (X / v p) (univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i) ≤ G := by
      refine le_trans (sum_le_sum (fun ρ _ => mul_le_mul_of_nonneg_right ?_
        (prod_nonneg (fun i _ => hk _)))) (tuple_bound h S c k hk X n)
      have hs0 : 0 ≤ univ.sup' univ_nonempty (fun i => c (ρ i)) := by
        obtain ⟨i, -, he⟩ := exists_mem_eq_sup' univ_nonempty (fun i => c (ρ i))
        rw [he]
        exact hc _
      exact rp_le_rp h (div_le_self hX.le hvp) hs0
    have hF' : |F p| ≤ CI * G := le_trans (hF p hp) (mul_le_mul_of_nonneg_left hsum hCI)
    have ha : rp ex lg (v p) (-(1 / 2) - y) ≤ rp ex lg (v p) (-(1 / 2)) :=
      rp_le_rp_exp h hvp (by linarith)
    exact mul_le_mul (mul_le_mul_of_nonneg_left ha (hC p hp)) hF' (abs_nonneg _)
      (mul_nonneg (hC p hp) (rp_nonneg h _ _))
  calc ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) - y) * |F p| ≤
        ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) * (CI * G) := sum_le_sum hterm
    _ = (∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2))) * (CI * G) := by rw [sum_mul]
    _ ≤ Kpk * (CI * G) := mul_le_mul_of_nonneg_right hP (mul_nonneg hCI hG0)
    _ = Kpk * CI * G := by ring

/-- (iv), the final shape: m W On^{m−1} ≤ k (1 + W)(1 + On)^{k−1} for 1 ≤ m = n + 1 ≤ k (W, On ≥ 0). -/
lemma mono_bound {W On : ℝ} (hW : 0 ≤ W) (hOn : 0 ≤ On) {n k : ℕ} (hnk : n + 1 ≤ k) :
    ((n : ℝ) + 1) * W * On ^ n ≤ k * (1 + W) * (1 + On) ^ (k - 1) := by
  have h1 : On ^ n ≤ (1 + On) ^ n := pow_le_pow_left₀ hOn (by linarith only) n
  have h2 : (1 + On) ^ n ≤ (1 + On) ^ (k - 1) := pow_le_pow_right₀ (by linarith only [hOn]) (by omega)
  have h3 : ((n : ℝ) + 1) ≤ k := by exact_mod_cast hnk
  have h4 : W ≤ 1 + W := by linarith only
  exact mul_le_mul (mul_le_mul h3 h4 hW (Nat.cast_nonneg k)) (le_trans h1 h2) (pow_nonneg hOn n)
    (mul_nonneg (Nat.cast_nonneg k) (by linarith only [hW]))

/-- **Lemma Z, summation over the partition/binomial terms** (case 1, s = 1/2 + y + it): NAMED (Lemma Z (i)–(ii):
explicit formula, Möbius inversion over set partitions, binomial expansion of (Q − Q_pp), Fourier inversion, Fubini)
|R_T(t − iy)| ≤ Σ_{τ∈𝒯} B_τ + E with E = O(T^{−10}) (pole and trivial-zero terms), B_τ the term's packet sum; each
term either has n_τ + 1 ≤ k free variables and B_τ ≤ K_pk C_I (n_τ+1) W On^{n_τ} (term_bound) or none and
B_τ ≤ K_pk C_I. Then |R_T(t − iy)| ≤ #𝒯 · K_pk C_I · k (1 + W)(1 + On)^{k−1} + E, i.e. Lemma Z's
|R_T| ≤ C(1 + W_{−y})(1 + On)^{k−1} (#𝒯 ≤ C(k)). -/
theorem lemmaZ_sum {γ : Type*} (𝒯 : Finset γ) (Bt : γ → ℝ) (nt : γ → ℕ) {Rabs E W On Kpk CI : ℝ}
    {k : ℕ} (hk : 1 ≤ k) (hW : 0 ≤ W) (hOn : 0 ≤ On) (hK : 0 ≤ Kpk * CI)
    (hR : Rabs ≤ ∑ τ ∈ 𝒯, Bt τ + E)
    (hB : ∀ τ ∈ 𝒯, nt τ + 1 ≤ k ∧
      (Bt τ ≤ Kpk * CI * (((nt τ : ℝ) + 1) * W * On ^ nt τ) ∨ Bt τ ≤ Kpk * CI)) :
    Rabs ≤ 𝒯.card * (Kpk * CI * (k * (1 + W) * (1 + On) ^ (k - 1))) + E := by
  have hG : 1 ≤ (k : ℝ) * (1 + W) * (1 + On) ^ (k - 1) := by
    have h1 : (1 : ℝ) ≤ k := by exact_mod_cast hk
    have h2 : (1 : ℝ) ≤ 1 + W := by linarith only [hW]
    have h3 : (1 : ℝ) ≤ (1 + On) ^ (k - 1) := one_le_pow₀ (by linarith only [hOn])
    exact one_le_mul_of_one_le_of_one_le (one_le_mul_of_one_le_of_one_le h1 h2) h3
  have hterm : ∀ τ ∈ 𝒯, Bt τ ≤ Kpk * CI * (k * (1 + W) * (1 + On) ^ (k - 1)) := by
    intro τ hτ
    obtain ⟨hn, h1 | h1⟩ := hB τ hτ
    · exact le_trans h1 (mul_le_mul_of_nonneg_left (mono_bound hW hOn hn) hK)
    · calc Bt τ ≤ Kpk * CI := h1
        _ = Kpk * CI * 1 := (mul_one _).symm
        _ ≤ Kpk * CI * (k * (1 + W) * (1 + On) ^ (k - 1)) := mul_le_mul_of_nonneg_left hG hK
  have hs := sum_le_sum hterm
  rw [sum_const, nsmul_eq_mul] at hs
  linarith only [hR, hs]

end DPCBLean

