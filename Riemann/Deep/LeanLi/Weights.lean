module

public import Riemann.Deep.LeanLi.Smooth
public import Riemann.Deep.LeanLi.LR.Core
public import Mathlib.Analysis.Calculus.BumpFunction.FiniteDimension

/-!
# LeanLi: positive-definite smooth weights (autocorrelations of a flat-top bump)

`v(m) = ψ(m/N)`, `ψ` a smooth bump `= 1` on `[-a, a]`, `= 0` off `(-1/2, 1/2)`;
`w(n) = Σ_m v(m) v(m+n) / Σ_m v(m)²`. Then `Φ(x) = Σ_{|n|≤N} w(|n|) e^{inx} = |Σ_m v(m) e^{imx}|² / Σ v² ≥ 0`.
-/

@[expose] public section

open Real Finset

namespace LeanLi

/-- The square-sum identity behind positive definiteness. -/
lemma autocorr_identity (M : ℕ) (v : ℤ → ℝ) (hv : ∀ k, k ∉ Icc (-(M : ℤ)) M → v k = 0) (x : ℝ) :
    (∑ m ∈ Icc (-(M : ℤ)) M, v m * cos (m * x)) ^ 2 + (∑ m ∈ Icc (-(M : ℤ)) M, v m * sin (m * x)) ^ 2 =
      ∑ n ∈ Icc (-(2 * M : ℤ)) (2 * M), (∑ m ∈ Icc (-(M : ℤ)) M, v m * v (m + n)) * cos (n * x) := by
  set S := Icc (-(M : ℤ)) M
  have h1 : (∑ m ∈ S, v m * cos (m * x)) ^ 2 + (∑ m ∈ S, v m * sin (m * x)) ^ 2 =
      ∑ m ∈ S, ∑ m' ∈ S, v m * v m' * cos ((m' - m : ℤ) * x) := by
    rw [sq, sq, Finset.sum_mul_sum, Finset.sum_mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun m' _ => ?_
    push_cast
    rw [show ((m' : ℝ) - m) * x = m' * x - m * x by ring, cos_sub]
    ring
  have h2 : ∀ m ∈ S, ∑ m' ∈ S, v m * v m' * cos ((m' - m : ℤ) * x) =
      ∑ n ∈ Icc (-(2 * M : ℤ)) (2 * M), v m * v (m + n) * cos (n * x) := by
    intro m hm
    rw [Finset.mem_Icc] at hm
    -- reindex m' = m + n
    have hmap : (Icc (-(2 * M : ℤ)) (2 * M)).map (addLeftEmbedding m) ⊇ S := by
      intro m' hm'
      rw [Finset.mem_Icc] at hm'
      rw [Finset.mem_map]
      exact ⟨m' - m, Finset.mem_Icc.2 ⟨by omega, by omega⟩, by simp⟩
    rw [Finset.sum_subset hmap ?_, Finset.sum_map]
    · refine Finset.sum_congr rfl fun n _ => ?_
      simp [addLeftEmbedding]
    · intro m' _ hm'
      rw [hv m' hm']; ring
  rw [h1, Finset.sum_congr rfl h2, Finset.sum_comm]
  refine Finset.sum_congr rfl fun n _ => ?_
  rw [Finset.sum_mul]

/-- `r(n) = Σ_m v(m) v(m+n)`. -/
noncomputable def autocorr (M : ℕ) (v : ℤ → ℝ) (n : ℤ) : ℝ := ∑ m ∈ Icc (-(M : ℤ)) M, v m * v (m + n)

lemma autocorr_eq_tsum (M : ℕ) (v : ℤ → ℝ) (hv : ∀ k, k ∉ Icc (-(M : ℤ)) M → v k = 0) (n : ℤ) :
    autocorr M v n = ∑' m : ℤ, v m * v (m + n) := by
  unfold autocorr
  rw [tsum_eq_sum (s := Icc (-(M : ℤ)) M) (fun m hm => by rw [hv m hm, zero_mul])]

lemma autocorr_neg (M : ℕ) (v : ℤ → ℝ) (hv : ∀ k, k ∉ Icc (-(M : ℤ)) M → v k = 0) (n : ℤ) :
    autocorr M v (-n) = autocorr M v n := by
  rw [autocorr_eq_tsum M v hv, autocorr_eq_tsum M v hv, ← (Equiv.addRight n).tsum_eq]
  refine tsum_congr fun m => ?_
  simp only [Equiv.coe_addRight]
  rw [show m + n + -n = m by ring, mul_comm]

/-- The cosine form: `Σ_{|n| ≤ 2M} r(n) cos(nx) = r(0) + 2 Σ_{n=1}^{2M} r(n) cos(nx)`. -/
lemma sum_Icc_sym_cos (M : ℕ) (r : ℤ → ℝ) (hr : ∀ n, r (-n) = r n) (x : ℝ) :
    ∑ n ∈ Icc (-(2 * M : ℤ)) (2 * M), r n * cos (n * x) =
      r 0 + 2 * ∑ n ∈ Icc (1 : ℤ) (2 * M), r n * cos (n * x) := by
  set K : ℤ := 2 * M
  have hK : 0 ≤ K := by positivity
  let negE : ℤ ↪ ℤ := ⟨Neg.neg, neg_injective⟩
  have hsplit : Icc (-K) K = (Icc 1 K).map negE ∪ ({0} ∪ Icc 1 K) := by
    ext n
    simp only [Finset.mem_Icc, Finset.mem_union, Finset.mem_map, Finset.mem_singleton]
    constructor
    · intro h
      rcases lt_trichotomy n 0 with hn | hn | hn
      · left; exact ⟨-n, ⟨by omega, by omega⟩, by simp [negE]⟩
      · right; left; exact hn
      · right; right; exact ⟨by omega, by omega⟩
    · rintro (⟨m, hm, rfl⟩ | h | h)
      · simp only [negE, Function.Embedding.coeFn_mk]; constructor <;> omega
      · subst h; constructor <;> omega
      · constructor <;> omega
  have hd1 : Disjoint ((Icc 1 K).map negE) ({0} ∪ Icc 1 K) := by
    rw [Finset.disjoint_left]
    intro n hn hn'
    simp only [Finset.mem_map, Finset.mem_Icc, negE, Function.Embedding.coeFn_mk] at hn
    simp only [Finset.mem_union, Finset.mem_singleton, Finset.mem_Icc] at hn'
    obtain ⟨m, hm, rfl⟩ := hn
    omega
  have hd2 : Disjoint ({0} : Finset ℤ) (Icc 1 K) := by
    rw [Finset.disjoint_singleton_left]; simp
  rw [hsplit, Finset.sum_union hd1, Finset.sum_union hd2, Finset.sum_map, Finset.sum_singleton]
  simp only [negE, Function.Embedding.coeFn_mk, hr, Int.cast_neg, neg_mul, cos_neg, Int.cast_zero,
    zero_mul, cos_zero, mul_one]
  ring

/-- `Σ_{n ∈ [1, N] ⊂ ℤ} g n = Σ_{n ∈ [1, N] ⊂ ℕ} g n`. -/
lemma sum_Icc_int_nat (N : ℕ) (g : ℤ → ℝ) :
    ∑ n ∈ Icc (1 : ℤ) N, g n = ∑ n ∈ Icc 1 N, g (n : ℤ) := by
  have h : (Icc 1 N).map Nat.castEmbedding = Icc (1 : ℤ) N := by
    ext n
    simp only [Finset.mem_map, Finset.mem_Icc, Nat.castEmbedding_apply]
    constructor
    · rintro ⟨m, ⟨h1, h2⟩, rfl⟩; constructor <;> omega
    · rintro ⟨h1, h2⟩; exact ⟨n.toNat, ⟨by omega, by omega⟩, by omega⟩
  rw [← h, Finset.sum_map]; rfl

/-! ## The concrete weights -/

/-- The flat-top bump: `= 1` on `[-a, a]`, `= 0` off `(-1/2, 1/2)`. -/
noncomputable def bumpF (a : ℝ) (ha : 0 < a) (ha' : a < 1 / 2) : ContDiffBump (0 : ℝ) :=
  ⟨a, 1 / 2, ha, ha'⟩

section Bump

variable (a : ℝ) (ha : 0 < a) (ha' : a < 1 / 2)

/-- `v_N(m) = ψ(m/N)`. -/
noncomputable def vN (N : ℕ) (m : ℤ) : ℝ := (bumpF a ha ha') ((m : ℝ) / N)

/-- The weights `w_N(n) = r(n)/r(0)`. -/
noncomputable def wts (N : ℕ) (n : ℕ) : ℝ := autocorr N (vN a ha ha' N) n / autocorr N (vN a ha ha' N) 0

lemma vN_nonneg (N : ℕ) (m : ℤ) : 0 ≤ vN a ha ha' N m := (bumpF a ha ha').nonneg

lemma vN_le_one (N : ℕ) (m : ℤ) : vN a ha ha' N m ≤ 1 := (bumpF a ha ha').le_one

lemma vN_eq_zero {N : ℕ} (hN : 1 ≤ N) {m : ℤ} (hm : (N : ℝ) / 2 ≤ |(m : ℝ)|) : vN a ha ha' N m = 0 := by
  unfold vN
  apply (bumpF a ha ha').zero_of_le_dist
  rw [Real.dist_eq, sub_zero, abs_div, abs_of_pos (by exact_mod_cast hN : (0 : ℝ) < N)]
  show 1 / 2 ≤ |(m : ℝ)| / N
  rw [le_div_iff₀ (by exact_mod_cast hN)]; linarith

lemma vN_eq_one {N : ℕ} {m : ℤ} (hm : |(m : ℝ)| ≤ a * N) (hN : 1 ≤ N) : vN a ha ha' N m = 1 := by
  unfold vN
  apply (bumpF a ha ha').one_of_mem_closedBall
  rw [Metric.mem_closedBall, Real.dist_eq, sub_zero, abs_div, abs_of_pos (by exact_mod_cast hN : (0 : ℝ) < N)]
  show |(m : ℝ)| / N ≤ a
  rw [div_le_iff₀ (by exact_mod_cast hN)]; linarith

lemma vN_supp {N : ℕ} (hN : 1 ≤ N) : ∀ k, k ∉ Icc (-(N : ℤ)) N → vN a ha ha' N k = 0 := by
  intro k hk
  apply vN_eq_zero a ha ha' hN
  rw [Finset.mem_Icc, not_and_or, not_le, not_le] at hk
  have hN' : (0 : ℝ) ≤ N := Nat.cast_nonneg N
  rcases hk with h | h
  · have : (k : ℝ) < -(N : ℝ) := by exact_mod_cast h
    rw [abs_of_neg (by linarith)]; linarith
  · have : (N : ℝ) < k := by exact_mod_cast h
    rw [abs_of_pos (by linarith)]; linarith

lemma autocorr_nonneg (N : ℕ) (n : ℤ) : 0 ≤ autocorr N (vN a ha ha' N) n :=
  Finset.sum_nonneg fun m _ => mul_nonneg (vN_nonneg a ha ha' N m) (vN_nonneg a ha ha' N _)

lemma autocorr_le_zero {N : ℕ} (hN : 1 ≤ N) (n : ℤ) :
    autocorr N (vN a ha ha' N) n ≤ autocorr N (vN a ha ha' N) 0 := by
  have hsupp := vN_supp a ha ha' hN
  rw [autocorr_eq_tsum N _ hsupp, autocorr_eq_tsum N _ hsupp]
  set v := vN a ha ha' N
  have hsum : ∀ k : ℤ, Summable fun m : ℤ => v m * v (m + k) := fun k =>
    summable_of_ne_finset_zero (s := Icc (-(N : ℤ)) N) fun m hm => by rw [hsupp m hm, zero_mul]
  have hsq : ∑' m : ℤ, v (m + n) * v (m + n) = ∑' m : ℤ, v m * v (m + 0) := by
    have := (Equiv.addRight n).tsum_eq (fun m => v m * v m)
    simp only [Equiv.coe_addRight] at this
    rw [this]; simp
  have h0 : Summable fun m : ℤ => v m * v m := (hsum 0).congr fun m => by simp
  have hs2 : Summable fun m : ℤ => v (m + n) * v (m + n) := (Equiv.addRight n).summable_iff.2 h0
  calc ∑' m : ℤ, v m * v (m + n) ≤ ∑' m : ℤ, (v m * v (m + 0) + v (m + n) * v (m + n)) / 2 := by
        refine (hsum n).tsum_le_tsum (fun m => ?_) (((hsum 0).add hs2).div_const 2)
        simp only [add_zero]
        nlinarith [sq_nonneg (v m - v (m + n))]
    _ = ∑' m : ℤ, v m * v (m + 0) := by
        rw [tsum_div_const, (hsum 0).tsum_add hs2, hsq]; ring

lemma sum_vN_ge {N : ℕ} (hN : 1 ≤ N) (f : ℝ → ℝ) (hf : ∀ y, 0 ≤ y → y ≤ 1 → 0 ≤ f y) (hf1 : f 1 = 1) :
    2 * a * N - 1 ≤ ∑ m ∈ Icc (-(N : ℤ)) N, f (vN a ha ha' N m) := by
  set k : ℤ := ⌊a * N⌋
  have hk0 : 0 ≤ k := Int.floor_nonneg.2 (by positivity)
  have hkN : k ≤ N := by
    have : a * N ≤ N := by nlinarith [(Nat.cast_nonneg N : (0:ℝ) ≤ N)]
    exact_mod_cast (Int.floor_le (a * N)).trans this
  have hsub : Icc (-k) k ⊆ Icc (-(N : ℤ)) N := by
    intro m hm; rw [Finset.mem_Icc] at hm ⊢; constructor <;> omega
  have hone : ∀ m ∈ Icc (-k) k, f (vN a ha ha' N m) = 1 := by
    intro m hm
    rw [Finset.mem_Icc] at hm
    have h1 : |(m : ℝ)| ≤ k := by
      rw [abs_le]; constructor
      · have : (-k : ℤ) ≤ m := hm.1; exact_mod_cast (by push_cast at this ⊢; exact this)
      · exact_mod_cast hm.2
    have h2 : (k : ℝ) ≤ a * N := Int.floor_le _
    rw [vN_eq_one a ha ha' (h1.trans h2) hN, hf1]
  calc 2 * a * N - 1 ≤ (2 * k + 1 : ℝ) := by
        have := Int.lt_floor_add_one (a * N)
        linarith
    _ = ∑ m ∈ Icc (-k) k, f (vN a ha ha' N m) := by
        rw [Finset.sum_congr rfl hone, Finset.sum_const, Int.card_Icc, nsmul_eq_mul, mul_one]
        rw [show k + 1 - -k = 2 * k + 1 by ring]
        have : (0 : ℤ) ≤ 2 * k + 1 := by omega
        rw [show (((2 * k + 1).toNat : ℕ) : ℝ) = (((2 * k + 1).toNat : ℤ) : ℝ) by norm_cast,
          Int.toNat_of_nonneg this]
        push_cast; ring
    _ ≤ ∑ m ∈ Icc (-(N : ℤ)) N, f (vN a ha ha' N m) :=
        Finset.sum_le_sum_of_subset_of_nonneg hsub fun m _ _ =>
          hf _ (vN_nonneg a ha ha' N m) (vN_le_one a ha ha' N m)

lemma autocorr_zero_bounds {N : ℕ} (hN : 1 ≤ N) :
    2 * a * N - 1 ≤ autocorr N (vN a ha ha' N) 0 ∧
      autocorr N (vN a ha ha' N) 0 ≤ ∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m := by
  have e : autocorr N (vN a ha ha' N) 0 = ∑ m ∈ Icc (-(N : ℤ)) N, (fun y => y * y) (vN a ha ha' N m) := by
    unfold autocorr; simp
  constructor
  · rw [e]; exact sum_vN_ge a ha ha' hN _ (fun y hy _ => mul_nonneg hy hy) (by norm_num)
  · rw [e]
    refine Finset.sum_le_sum fun m _ => ?_
    have h0 := vN_nonneg a ha ha' N m
    have h1 := vN_le_one a ha ha' N m
    simp only; nlinarith

lemma autocorr_eq_zero_of_ge {N : ℕ} (hN : 1 ≤ N) {n : ℤ} (hn : (N : ℤ) ≤ n) :
    autocorr N (vN a ha ha' N) n = 0 := by
  unfold autocorr
  refine Finset.sum_eq_zero fun m _ => ?_
  by_cases hm : (N : ℝ) / 2 ≤ |(m : ℝ)|
  · rw [vN_eq_zero a ha ha' hN hm, zero_mul]
  · have hm' : |(m : ℝ)| < N / 2 := not_le.1 hm
    have : (N : ℝ) / 2 ≤ |((m + n : ℤ) : ℝ)| := by
      have hn' : (N : ℝ) ≤ n := by exact_mod_cast hn
      push_cast
      rw [abs_lt] at hm'
      rw [le_abs]; left; linarith [hm'.1]
    rw [vN_eq_zero a ha ha' hN this, mul_zero]

lemma autocorr_pos {N : ℕ} (hN : 1 ≤ N) : 0 < autocorr N (vN a ha ha' N) 0 := by
  have h := (autocorr_zero_bounds a ha ha' hN).1
  have h0 : ∀ m, 0 ≤ vN a ha ha' N m := vN_nonneg a ha ha' N
  have hv0 : vN a ha ha' N 0 = 1 := vN_eq_one a ha ha' (by simp; positivity) hN
  have e : autocorr N (vN a ha ha' N) 0 = ∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * vN a ha ha' N m := by
    unfold autocorr; simp
  have : vN a ha ha' N 0 * vN a ha ha' N 0 ≤ ∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * vN a ha ha' N m :=
    Finset.single_le_sum (f := fun m => vN a ha ha' N m * vN a ha ha' N m)
      (fun m _ => mul_nonneg (h0 m) (h0 _)) (Finset.mem_Icc.2 ⟨by omega, by omega⟩)
  rw [hv0, mul_one] at this
  rw [e]; linarith

lemma wts_nonneg (N n : ℕ) : 0 ≤ wts a ha ha' N n := by
  unfold wts
  exact div_nonneg (autocorr_nonneg a ha ha' N n) (autocorr_nonneg a ha ha' N 0)

lemma wts_le_one (N n : ℕ) : wts a ha ha' N n ≤ 1 := by
  rcases Nat.eq_zero_or_pos N with h0 | hN
  · subst h0
    unfold wts
    exact div_le_one_of_le₀ (by
      unfold autocorr; simp only [Nat.cast_zero, neg_zero, Finset.Icc_self, Finset.sum_singleton, zero_add]
      unfold vN; simp) (autocorr_nonneg a ha ha' 0 0)
  unfold wts
  rw [div_le_one (autocorr_pos a ha ha' hN)]
  exact autocorr_le_zero a ha ha' hN n

lemma wts_zero {N : ℕ} (hN : 1 ≤ N) : wts a ha ha' N 0 = 1 := by
  unfold wts; push_cast; exact div_self (autocorr_pos a ha ha' hN).ne'

/-- `r(0) Φ(x) = |Σ v(m) e^{imx}|²`. -/
lemma phi_mul_autocorr {N : ℕ} (hN : 1 ≤ N) (x : ℝ) :
    autocorr N (vN a ha ha' N) 0 * LR.Phi N (wts a ha ha' N) x =
      (∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * cos (m * x)) ^ 2 +
        (∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * sin (m * x)) ^ 2 := by
  have hsupp := vN_supp a ha ha' hN
  rw [autocorr_identity N _ hsupp x]
  set r := fun n : ℤ => autocorr N (vN a ha ha' N) n
  have hr : ∀ n, r (-n) = r n := fun n => autocorr_neg N _ hsupp n
  have e : ∀ n : ℤ, (∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * vN a ha ha' N (m + n)) = r n :=
    fun n => rfl
  simp only [e]
  rw [sum_Icc_sym_cos N r hr x]
  -- truncate the cosine sum at N
  have htrunc : ∑ n ∈ Icc (1 : ℤ) (2 * N), r n * cos (n * x) = ∑ n ∈ Icc (1 : ℤ) N, r n * cos (n * x) := by
    symm
    refine Finset.sum_subset (fun n hn => ?_) fun n hn hn' => ?_
    · rw [Finset.mem_Icc] at hn ⊢; constructor <;> omega
    · rw [Finset.mem_Icc] at hn hn'
      have : (N : ℤ) ≤ n := by omega
      simp only [r]; rw [autocorr_eq_zero_of_ge a ha ha' hN this, zero_mul]
  rw [htrunc, sum_Icc_int_nat]
  have hpos := autocorr_pos a ha ha' hN
  have hPhi : LR.Phi N (wts a ha ha' N) x * r 0 =
      r 0 + 2 * ∑ n ∈ Icc 1 N, r (n : ℤ) * cos ((n : ℝ) * x) := by
    unfold LR.Phi wts
    rw [add_mul, mul_assoc, Finset.sum_mul]
    have h0 : autocorr N (vN a ha ha' N) ((0 : ℕ) : ℤ) / autocorr N (vN a ha ha' N) 0 * r 0 = r 0 := by
      simp only [Nat.cast_zero, r]; field_simp
    rw [h0]
    congr 2
    refine Finset.sum_congr rfl fun n _ => ?_
    simp only [r]; field_simp
  push_cast
  rw [← hPhi]; ring

lemma phi_nonneg (N : ℕ) (x : ℝ) : 0 ≤ LR.Phi N (wts a ha ha' N) x := by
  rcases Nat.eq_zero_or_pos N with h0 | hN
  · subst h0; unfold LR.Phi
    have : Finset.Icc 1 0 = (∅ : Finset ℕ) := by decide
    rw [this, Finset.sum_empty, mul_zero, add_zero]
    exact wts_nonneg a ha ha' 0 0
  have h := phi_mul_autocorr a ha ha' hN x
  have hpos := autocorr_pos a ha ha' hN
  have : 0 ≤ autocorr N (vN a ha ha' N) 0 * LR.Phi N (wts a ha ha' N) x := by
    rw [h]; positivity
  exact (mul_nonneg_iff_of_pos_left hpos).1 this

lemma autocorr_zero_ge_one {N : ℕ} (hN : 1 ≤ N) : 1 ≤ autocorr N (vN a ha ha' N) 0 := by
  have h0 : ∀ m, 0 ≤ vN a ha ha' N m := vN_nonneg a ha ha' N
  have hv0 : vN a ha ha' N 0 = 1 := vN_eq_one a ha ha' (by simp; positivity) hN
  have e : autocorr N (vN a ha ha' N) 0 = ∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * vN a ha ha' N m := by
    unfold autocorr; simp
  have : vN a ha ha' N 0 * vN a ha ha' N 0 ≤ ∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * vN a ha ha' N m :=
    Finset.single_le_sum (f := fun m => vN a ha ha' N m * vN a ha ha' N m)
      (fun m _ => mul_nonneg (h0 m) (h0 _)) (Finset.mem_Icc.2 ⟨by omega, by omega⟩)
  rw [hv0, mul_one] at this
  rw [e]; exact this

lemma autocorr_zero_ge_aN {N : ℕ} (hN : 1 ≤ N) : a * N ≤ autocorr N (vN a ha ha' N) 0 := by
  have h1 := autocorr_zero_ge_one a ha ha' hN
  have h2 := (autocorr_zero_bounds a ha ha' hN).1
  by_cases h : a * N ≤ 1
  · linarith
  · push Not at h; linarith

/-- `Φ(0) ≥ 2aN - 1`. -/
lemma phi_zero_ge {N : ℕ} (hN : 1 ≤ N) : 2 * a * N - 1 ≤ LR.Phi N (wts a ha ha' N) 0 := by
  have h := phi_mul_autocorr a ha ha' hN 0
  simp only [mul_zero, cos_zero, mul_one, sin_zero, Finset.sum_const_zero] at h
  have hpos := autocorr_pos a ha ha' hN
  obtain ⟨hlo, hup⟩ := autocorr_zero_bounds a ha ha' hN
  set Sv := ∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m
  have hSv : 2 * a * N - 1 ≤ Sv := sum_vN_ge a ha ha' hN id (fun y hy _ => hy) rfl
  set r0 := autocorr N (vN a ha ha' N) 0
  have hPhi : LR.Phi N (wts a ha ha' N) 0 = Sv ^ 2 / r0 := by
    rw [eq_div_iff hpos.ne', mul_comm, h]; ring
  rw [hPhi, le_div_iff₀ hpos]
  have hSv0 : 0 ≤ Sv := hpos.le.trans hup
  rcases le_or_gt (2 * a * N - 1) 0 with h' | h'
  · nlinarith [sq_nonneg Sv]
  · nlinarith

/-- The Fourier sum of `v` is the smooth exponential sum. -/
lemma esum_vN_re_im {N : ℕ} (hN : 1 ≤ N) (x : ℝ) :
    (esum (fun n => (bumpF a ha ha') ((n : ℝ) / N)) x).re = ∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * cos (m * x) ∧
    (esum (fun n => (bumpF a ha ha') ((n : ℝ) / N)) x).im = ∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * sin (m * x) := by
  have hsupp := vN_supp a ha ha' hN
  unfold esum
  rw [tsum_eq_sum (s := Icc (-(N : ℤ)) N) (fun m hm => by
    have := hsupp m hm; unfold vN at this; simp only [this]; simp)]
  rw [Complex.re_sum, Complex.im_sum]
  constructor <;> refine Finset.sum_congr rfl fun m _ => ?_ <;> unfold vN
  · rw [show (Complex.I * (((m : ℝ) : ℂ) * (x : ℂ))) = (((m : ℝ) * x : ℝ) : ℂ) * Complex.I by
        push_cast; ring, Complex.re_ofReal_mul, Complex.exp_ofReal_mul_I_re]
  · rw [show (Complex.I * (((m : ℝ) : ℂ) * (x : ℂ))) = (((m : ℝ) * x : ℝ) : ℂ) * Complex.I by
        push_cast; ring, Complex.im_ofReal_mul, Complex.exp_ofReal_mul_I_im]

/-- **Decay of `Φ`**: `Φ(x)(1 + N|x|)^{2K} ≤ C N` for `|x| ≤ π`. -/
theorem phi_decay (K : ℕ) : ∃ C : ℝ, 0 ≤ C ∧ ∀ N : ℕ, 1 ≤ N → ∀ x : ℝ, |x| ≤ π →
    LR.Phi N (wts a ha ha' N) x * (1 + N * |x|) ^ (2 * K) ≤ C * N := by
  obtain ⟨C, hC0, hC⟩ := esum_decay (bumpF a ha ha').contDiff (bumpF a ha ha').hasCompactSupport K
  refine ⟨C ^ 2 / a, by positivity, fun N hN x hx => ?_⟩
  have hN' : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have h1 := hC N hN' x hx
  have hpos := autocorr_pos a ha ha' hN
  have hr0 := autocorr_zero_ge_aN a ha ha' hN
  have hP := phi_mul_autocorr a ha ha' hN x
  obtain ⟨hre, him⟩ := esum_vN_re_im a ha ha' hN x
  set E := esum (fun n => (bumpF a ha ha') ((n : ℝ) / N)) x
  have hnorm : ‖E‖ ^ 2 = (∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * cos (m * x)) ^ 2 +
      (∑ m ∈ Icc (-(N : ℤ)) N, vN a ha ha' N m * sin (m * x)) ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply, hre, him]; ring
  set Φ := LR.Phi N (wts a ha ha' N) x
  set q := (1 + N * |x|) ^ K
  have hq : 0 < q := by positivity
  have hPhi : Φ = ‖E‖ ^ 2 / autocorr N (vN a ha ha' N) 0 := by
    rw [hnorm, ← hP]; field_simp
  have hE : ‖E‖ * q ≤ C * N := h1
  have hE2 : (‖E‖ * q) ^ 2 ≤ (C * N) ^ 2 := pow_le_pow_left₀ (by positivity) hE 2
  rw [hPhi, show (1 + N * |x|) ^ (2 * K) = q ^ 2 by rw [pow_mul']]
  rw [div_mul_eq_mul_div, div_le_iff₀ hpos]
  have hNpos : (0 : ℝ) < N := by linarith
  calc ‖E‖ ^ 2 * q ^ 2 = (‖E‖ * q) ^ 2 := by ring
    _ ≤ (C * N) ^ 2 := hE2
    _ = C ^ 2 / a * N * (a * N) := by field_simp
    _ ≤ C ^ 2 / a * N * autocorr N (vN a ha ha' N) 0 :=
        mul_le_mul_of_nonneg_left hr0 (by positivity)

end Bump

end LeanLi
