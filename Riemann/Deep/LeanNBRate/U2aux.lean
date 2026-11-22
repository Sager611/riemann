import Riemann.Deep.LeanNBRate.U2grp

/-!
# U2 auxiliaries: zero-free gaps, the finite low set, multiplicity bounds, the zero sum
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem finset_avoid (F : Finset ℝ) (x : ℝ) (hx : x ∉ F) : ∃ d > 0, ∀ f ∈ F, d ≤ |f - x| := by
  classical
  induction F using Finset.induction_on with
  | empty => exact ⟨1, one_pos, fun f hf => absurd hf (Finset.notMem_empty f)⟩
  | insert a F ha ih =>
    have hxa : x ≠ a := fun h => hx (h ▸ Finset.mem_insert_self a F)
    have hxF : x ∉ F := fun h => hx (Finset.mem_insert_of_mem h)
    obtain ⟨d, hd, hdF⟩ := ih hxF
    have hpos : 0 < |a - x| := abs_pos.mpr (sub_ne_zero.mpr (Ne.symm hxa))
    refine ⟨min d |a - x|, lt_min hd hpos, fun f hf => ?_⟩
    rcases Finset.mem_insert.mp hf with rfl | hf
    · exact min_le_right _ _
    · exact (min_le_left _ _).trans (hdF f hf)

/-- Absolute ordinates of the zeros up to height b form a finite set. -/
theorem absIm_finite (b : ℝ) : {t : ℝ | ∃ ρ ∈ NTZ, |ρ.im| ≤ b ∧ t = |ρ.im|}.Finite := by
  have h : {t : ℝ | ∃ ρ ∈ NTZ, |ρ.im| ≤ b ∧ t = |ρ.im|} ⊆ (fun ρ : ℂ => |ρ.im|) '' zerosIn (-b - 1) b := by
    rintro t ⟨ρ, hρ, hb, rfl⟩
    refine ⟨ρ, ⟨hρ, ?_, ?_⟩, rfl⟩
    · have := neg_abs_le ρ.im; linarith
    · exact (le_abs_self _).trans hb
  exact ((zerosIn_finite _ _).image _).subset h

/-- **Zero-free gap**: every interval [H, H+1] contains a height Hs with a zero-free annulus. -/
theorem gap_exists (H : ℝ) : ∃ Hs d : ℝ, H ≤ Hs ∧ Hs ≤ H + 1 ∧ 0 < d ∧ d ≤ 1 / 2 ∧
    ∀ ρ ∈ NTZ, |ρ.im| ≤ Hs - d ∨ Hs + d ≤ |ρ.im| := by
  classical
  set F := (absIm_finite (H + 2)).toFinset
  obtain ⟨x, hxI, hxF⟩ := (Set.Ioo_infinite (by linarith : H < H + 1)).exists_notMem_finset F
  obtain ⟨d, hd, hdF⟩ := finset_avoid F x hxF
  refine ⟨x, min d (1 / 2), hxI.1.le, hxI.2.le, lt_min hd (by norm_num), min_le_right _ _,
    fun ρ hρ => ?_⟩
  by_cases hb : |ρ.im| ≤ H + 2
  · have hmem : |ρ.im| ∈ F := (absIm_finite (H + 2)).mem_toFinset.mpr ⟨ρ, hρ, hb, rfl⟩
    have h1 := hdF _ hmem
    have h2 : min d (1 / 2) ≤ |(|ρ.im|) - x| := (min_le_left _ _).trans h1
    rcases le_or_gt (|ρ.im|) x with h | h
    · left; rw [abs_of_nonpos (by linarith)] at h2; linarith
    · right; rw [abs_of_pos (by linarith)] at h2; linarith
  · right
    push Not at hb
    have := min_le_right d (1 / 2 : ℝ)
    linarith [hxI.2]

/-- The finite set of low ordinates (zeros on the critical line, height ≤ b). -/
theorem lowOrd_finite (b : ℝ) : {γ : ℝ | ((1 / 2 : ℂ) + (γ : ℂ) * I) ∈ NTZ ∧ |γ| ≤ b}.Finite := by
  have h : {γ : ℝ | ((1 / 2 : ℂ) + (γ : ℂ) * I) ∈ NTZ ∧ |γ| ≤ b} ⊆ (fun ρ : ℂ => ρ.im) '' zerosIn (-b - 1) b := by
    rintro γ ⟨hγ, hb⟩
    refine ⟨(1 / 2 : ℂ) + (γ : ℂ) * I, ⟨hγ, ?_, ?_⟩, by simp⟩
    · rw [im_half_add]; have := neg_abs_le γ; linarith
    · rw [im_half_add]; exact (le_abs_self _).trans hb
  exact ((zerosIn_finite _ _).image _).subset h

def lowOrd (b : ℝ) : Finset ℝ := (lowOrd_finite b).toFinset

theorem mem_lowOrd {b γ : ℝ} : γ ∈ lowOrd b ↔ ((1 / 2 : ℂ) + (γ : ℂ) * I) ∈ NTZ ∧ |γ| ≤ b := by
  unfold lowOrd; rw [Set.Finite.mem_toFinset]; rfl

theorem nodeSet_low_subset {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) (b : ℝ) :
    nodeSet P {γ | |γ| ≤ b} ⊆ lowOrd b := by
  intro γ hγ
  exact mem_lowOrd.mpr ⟨node_mem_NTZ hP hγ, (mem_nodeSet.mp hγ).2⟩

theorem finset_gap (S : Finset ℝ) : ∃ g₀ > 0, ∀ γ ∈ S, ∀ γ' ∈ S, γ ≠ γ' → g₀ ≤ |γ - γ'| := by
  classical
  induction S using Finset.induction_on with
  | empty => exact ⟨1, one_pos, fun γ h => absurd h (Finset.notMem_empty γ)⟩
  | insert a S ha ih =>
    obtain ⟨g, hg, hgS⟩ := ih
    obtain ⟨d, hd, hdS⟩ := finset_avoid S a ha
    refine ⟨min g d, lt_min hg hd, fun γ hγ γ' hγ' hne => ?_⟩
    by_cases h1 : γ = a
    · have h2 : γ' ∈ S := by
        rcases Finset.mem_insert.mp hγ' with h | h
        · exact absurd (h1.trans h.symm) hne
        · exact h
      rw [h1, abs_sub_comm]; exact (min_le_right _ _).trans (hdS γ' h2)
    · have h3 : γ ∈ S := (Finset.mem_insert.mp hγ).resolve_left h1
      by_cases h4 : γ' = a
      · rw [h4]; exact (min_le_right _ _).trans (hdS γ h3)
      · have h5 : γ' ∈ S := (Finset.mem_insert.mp hγ').resolve_left h4
        exact (min_le_left _ _).trans (hgS γ h3 γ' h5 hne)

/-- Multiplicities grow at most like √|γ| (Titchmarsh 9.2). -/
theorem mult_le_sqrt (hc : Titchmarsh_Thm9_2) : ∃ K : ℝ, 0 ≤ K ∧ ∀ ρ ∈ NTZ,
    (mult ρ : ℝ) ≤ K * Real.sqrt (|ρ.im| + 1) := by
  classical
  obtain ⟨A, T₀, hT₀, hA⟩ := hc
  set b := max T₀ 1 + 1
  set Z := (zerosIn_finite (-b - 1) b).toFinset
  set Ms : ℝ := ∑ ρ ∈ Z, (mult ρ : ℝ)
  refine ⟨max Ms (2 * |A|), le_max_of_le_left (Finset.sum_nonneg fun _ _ => Nat.cast_nonneg _),
    fun ρ hρ => ?_⟩
  have hs1 : 1 ≤ Real.sqrt (|ρ.im| + 1) := Real.one_le_sqrt.mpr (by linarith [abs_nonneg ρ.im])
  by_cases hsmall : |ρ.im| ≤ b
  · have hmem : ρ ∈ Z := (zerosIn_finite _ _).mem_toFinset.mpr
      ⟨hρ, by have := neg_abs_le ρ.im; linarith, (le_abs_self _).trans hsmall⟩
    have h1 : (mult ρ : ℝ) ≤ Ms := Finset.single_le_sum (f := fun ρ => (mult ρ : ℝ))
      (fun _ _ => Nat.cast_nonneg _) hmem
    calc (mult ρ : ℝ) ≤ Ms := h1
      _ ≤ max Ms (2 * |A|) := le_max_left _ _
      _ = max Ms (2 * |A|) * 1 := by ring
      _ ≤ max Ms (2 * |A|) * Real.sqrt (|ρ.im| + 1) := by
          apply mul_le_mul_of_nonneg_left hs1
          exact le_max_of_le_left (Finset.sum_nonneg fun _ _ => Nat.cast_nonneg _)
  · push Not at hsmall
    set t := |ρ.im|
    have hcount : (mult ρ : ℝ) ≤ Ncount (t - 1) t := by
      rcases le_or_gt 0 ρ.im with hpos | hneg
      · have ht : t = ρ.im := abs_of_nonneg hpos
        have := sum_le_Ncount (a := t - 1) (b := t) {ρ} (fun x hx => by
          rw [Finset.mem_singleton] at hx; subst hx
          exact ⟨hρ, by rw [← ht]; linarith, by rw [← ht]⟩)
        rw [Finset.sum_singleton] at this
        exact_mod_cast this
      · have ht : t = -ρ.im := abs_of_neg hneg
        have hc' : conj ρ ∈ NTZ := conj_mem_NTZ hρ
        have := sum_le_Ncount (a := t - 1) (b := t) {conj ρ} (fun x hx => by
          rw [Finset.mem_singleton] at hx; subst hx
          refine ⟨hc', ?_, ?_⟩
          · rw [Complex.conj_im, ← ht]; linarith
          · rw [Complex.conj_im, ← ht])
        rw [Finset.sum_singleton, mult_conj ρ (NTZ_ne_one hρ)] at this
        exact_mod_cast this
    have hmax1 := le_max_left T₀ 1
    have hmax2 := le_max_right T₀ 1
    have hT : T₀ ≤ t - 1 := by simp only [b] at hsmall; linarith
    have h2 := hA (t - 1) hT
    rw [show t - 1 + 1 = t by ring] at h2
    have ht1 : 0 < t - 1 := by linarith
    have hlog : Real.log (t - 1) ≤ 2 * Real.sqrt (t + 1) := by
      have h := Real.log_le_rpow_div ht1.le (by norm_num : (0 : ℝ) < 1 / 2)
      have hs : (t - 1) ^ (1 / 2 : ℝ) ≤ Real.sqrt (t + 1) := by
        rw [← Real.sqrt_eq_rpow]; exact Real.sqrt_le_sqrt (by linarith)
      calc Real.log (t - 1) ≤ (t - 1) ^ (1 / 2 : ℝ) / (1 / 2) := h
        _ = 2 * (t - 1) ^ (1 / 2 : ℝ) := by ring
        _ ≤ 2 * Real.sqrt (t + 1) := by linarith
    have hlog0 : 0 ≤ Real.log (t - 1) := Real.log_nonneg (by simp only [b] at hsmall; linarith)
    calc (mult ρ : ℝ) ≤ Ncount (t - 1) t := hcount
      _ ≤ A * Real.log (t - 1) := h2
      _ ≤ |A| * Real.log (t - 1) := mul_le_mul_of_nonneg_right (le_abs_self A) hlog0
      _ ≤ |A| * (2 * Real.sqrt (t + 1)) := mul_le_mul_of_nonneg_left hlog (abs_nonneg A)
      _ = (2 * |A|) * Real.sqrt (t + 1) := by ring
      _ ≤ max Ms (2 * |A|) * Real.sqrt (t + 1) :=
          mul_le_mul_of_nonneg_right (le_max_right _ _) (Real.sqrt_nonneg _)

theorem norm_sq_NTZ (hRH : RiemannHypothesis) {ρ : ℂ} (hρ : ρ ∈ NTZ) : ‖ρ‖ ^ 2 = 1 / 4 + ρ.im ^ 2 := by
  rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply, NTZ_re_half hRH hρ]; ring

/-- **Summability of the zero sum** (RH, Titchmarsh 9.2). -/
theorem summable_zeroSum (hc : Titchmarsh_Thm9_2) (hRH : RiemannHypothesis) :
    Summable fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2 := by
  obtain ⟨K, hK, hmK⟩ := mult_le_sqrt hc
  have hs := (summable_wF1 hc).mul_left K
  refine Summable.of_nonneg_of_le (fun ρ => by positivity) (fun ρ => ?_) hs
  have hρ := ρ.2
  have hm := hmK _ hρ
  set m : ℝ := (mult (ρ : ℂ) : ℝ)
  set t := |(ρ : ℂ).im|
  have hm0 : 0 ≤ m := Nat.cast_nonneg _
  rw [norm_sq_NTZ hRH hρ]
  have ht0 : 0 ≤ t := abs_nonneg _
  have hsq : (ρ : ℂ).im ^ 2 = t ^ 2 := (sq_abs _).symm
  rw [hsq]
  have hden : 0 < 1 / 4 + t ^ 2 := by positivity
  have hF : Real.sqrt (t + 1) / (1 / 4 + t ^ 2) ≤ F1 t := by
    unfold F1
    have hs0 : 0 ≤ Real.sqrt (t + 1) := Real.sqrt_nonneg _
    have h2 : 1 ≤ Real.sqrt 2 := Real.one_le_sqrt.mpr (by norm_num)
    rcases le_or_gt 1 t with h1 | h1
    · rw [max_eq_left h1]
      have ht2 : 0 < t ^ 2 := by positivity
      calc Real.sqrt (t + 1) / (1 / 4 + t ^ 2) ≤ Real.sqrt (t + 1) / t ^ 2 :=
            div_le_div_of_nonneg_left hs0 ht2 (by linarith)
        _ ≤ 4 * Real.sqrt 2 * Real.sqrt (t + 1) / t ^ 2 := by
            gcongr
            nlinarith
    · rw [max_eq_right h1.le]
      have hs1 : Real.sqrt (t + 1) ≤ 2 := by
        rw [show (2 : ℝ) = Real.sqrt 4 by
          rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]]
        exact Real.sqrt_le_sqrt (by linarith)
      have h22 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)
      have hq : (1 : ℝ) + 1 = 2 := by norm_num
      rw [hq, div_le_iff₀ hden]
      norm_num
      nlinarith [sq_nonneg t]
  calc m ^ 2 / (1 / 4 + t ^ 2) = m * (m / (1 / 4 + t ^ 2)) := by ring
    _ ≤ m * (K * Real.sqrt (t + 1) / (1 / 4 + t ^ 2)) := by gcongr
    _ = m * K * (Real.sqrt (t + 1) / (1 / 4 + t ^ 2)) := by ring
    _ ≤ m * K * F1 t := by gcongr
    _ = K * (m * F1 t) := by ring

/-- The low part of the dual sum is at most the zero sum. -/
theorem low_sum_le_zeroSum (hc : Titchmarsh_Thm9_2) (hRH : RiemannHypothesis) {P : ExpSum} {T : ℝ}
    (hP : ZeroConfig P T) (Tl : Set ℝ) :
    ∑ γ ∈ nodeSet P Tl, ((P.part {γ} γ).M : ℝ) ^ 2 / (γ ^ 2 + 1 / 4) ≤ zeroSum := by
  classical
  set A := nodeSet P Tl
  set e : A → NTZ := fun γ => ⟨(1 / 2 : ℂ) + (γ.1 : ℂ) * I, node_mem_NTZ hP γ.2⟩
  have he : Function.Injective e := by
    intro x y hxy
    have := congrArg (fun ρ : NTZ => (ρ : ℂ).im) hxy
    simp [e] at this
    exact this
  have hM : ∀ γ ∈ A, (P.part {γ} γ).M = nw P γ := fun γ _ => by
    rw [ExpSum.part_M]; unfold nw
    refine Finset.sum_congr rfl fun j _ => ?_
    by_cases h : P.node j = γ <;> simp [h]
  have h1 : ∑ γ ∈ A, ((P.part {γ} γ).M : ℝ) ^ 2 / (γ ^ 2 + 1 / 4) ≤
      ∑ γ ∈ A.attach, (mult ((e γ : NTZ) : ℂ) : ℝ) ^ 2 / ‖((e γ : NTZ) : ℂ)‖ ^ 2 := by
    rw [← Finset.sum_attach A]
    refine Finset.sum_le_sum fun γ _ => ?_
    have h2 : ‖((e γ : NTZ) : ℂ)‖ ^ 2 = γ.1 ^ 2 + 1 / 4 := by
      simp only [e]
      rw [← Complex.normSq_eq_norm_sq, Complex.normSq_apply]; simp; ring
    rw [h2, hM γ.1 γ.2]
    have h3 : (nw P γ.1 : ℝ) ≤ mult ((e γ : NTZ) : ℂ) := by exact_mod_cast nw_le_mult hP γ.2
    have h4 : (0 : ℝ) ≤ nw P γ.1 := Nat.cast_nonneg _
    gcongr
  refine h1.trans ?_
  unfold zeroSum
  rw [← Finset.sum_image (f := fun ρ : NTZ => (mult (ρ : ℂ) : ℝ) ^ 2 / ‖(ρ : ℂ)‖ ^ 2) (g := e)
    (s := A.attach) (fun x _ y _ hxy => he hxy)]
  exact (summable_zeroSum hc hRH).sum_le_tsum _ fun ρ _ => by positivity

theorem zeroSum_nonneg : 0 ≤ zeroSum := tsum_nonneg fun ρ => by positivity

theorem claimFC_nonneg {c : ℝ} (hc : 1 ≤ c) (M₀ : ℕ) {H₀ g₀ : ℝ} (hH : 0 ≤ H₀) (hg : 0 < g₀) (J : ℕ) :
    0 ≤ claimFC c M₀ H₀ g₀ J := by
  unfold claimFC
  have := markovC_nonneg (by linarith : (0 : ℝ) < c) M₀
  positivity

theorem claimFC_mono {c : ℝ} (hc : 1 ≤ c) (M₀ : ℕ) {H₀ g₀ : ℝ} (hH : 0 ≤ H₀) (hg : 0 < g₀) {J J' : ℕ}
    (hJ : J ≤ J') : claimFC c M₀ H₀ g₀ J ≤ claimFC c M₀ H₀ g₀ J' := by
  unfold claimFC
  have := markovC_nonneg (by linarith : (0 : ℝ) < c) M₀
  have hJ' : (J : ℝ) ≤ J' := by exact_mod_cast hJ
  gcongr

end Riemann.Deep.LeanNBRate

