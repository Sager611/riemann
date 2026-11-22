import Riemann.Deep.LeanNBRate.U2parts

/-!
# Zero configurations: node weights, counts, and the RH bounds for high groups
-/

open MeasureTheory Filter Topology Set Complex ComplexConjugate

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- An exponential sum whose nodes are distinct ordinates of nontrivial zeros (height ≤ T), with
multiplicities between 1 and the zero's multiplicity. -/
structure ZeroConfig (P : ExpSum) (T : ℝ) : Prop where
  inj : Function.Injective P.node
  zero : ∀ j, ((1 / 2 : ℂ) + (P.node j : ℂ) * I) ∈ NTZ
  height : ∀ j, |P.node j| ≤ T
  mult_pos : ∀ j, 1 ≤ P.mult j
  mult_le : ∀ j, P.mult j ≤ mult ((1 / 2 : ℂ) + (P.node j : ℂ) * I)

/-- Weight of a node. -/
def nw (P : ExpSum) (γ : ℝ) : ℕ := ∑ j, if P.node j = γ then P.mult j else 0

open Classical in
theorem part_M_eq (P : ExpSum) (g : Finset ℝ) (c : ℝ) : (P.part ↑g c).M = ∑ γ ∈ g, nw P γ := by
  rw [ExpSum.part_M]
  unfold nw
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  by_cases hj : P.node j ∈ (g : Set ℝ)
  · rw [if_pos hj, Finset.sum_ite_eq g (P.node j) (fun _ => P.mult j),
      if_pos (by exact_mod_cast hj)]
  · rw [if_neg hj]
    symm
    refine Finset.sum_eq_zero fun γ hγ => ?_
    rw [if_neg]
    intro h; exact hj (by rw [h]; exact_mod_cast hγ)

theorem nw_node {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) (j : Fin P.J) : nw P (P.node j) = P.mult j := by
  unfold nw
  rw [Finset.sum_eq_single j]
  · rw [if_pos rfl]
  · intro i _ hij
    rw [if_neg]
    intro h; exact hij (hP.inj h)
  · intro h; exact absurd (Finset.mem_univ j) h

theorem nw_of_mem {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {T' : Set ℝ} {γ : ℝ}
    (hγ : γ ∈ nodeSet P T') : ∃ j, P.node j = γ ∧ nw P γ = P.mult j := by
  obtain ⟨⟨j, hj⟩, _⟩ := mem_nodeSet.mp hγ
  exact ⟨j, hj, by rw [← hj]; exact nw_node hP j⟩

theorem nw_pos {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {T' : Set ℝ} {γ : ℝ}
    (hγ : γ ∈ nodeSet P T') : 1 ≤ nw P γ := by
  obtain ⟨j, _, h⟩ := nw_of_mem hP hγ
  rw [h]; exact hP.mult_pos j

theorem nw_le_mult {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {T' : Set ℝ} {γ : ℝ}
    (hγ : γ ∈ nodeSet P T') : nw P γ ≤ mult ((1 / 2 : ℂ) + (γ : ℂ) * I) := by
  obtain ⟨j, hj, h⟩ := nw_of_mem hP hγ
  rw [h, ← hj]; exact hP.mult_le j

theorem node_mem_NTZ {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {T' : Set ℝ} {γ : ℝ}
    (hγ : γ ∈ nodeSet P T') : ((1 / 2 : ℂ) + (γ : ℂ) * I) ∈ NTZ := by
  obtain ⟨j, hj, _⟩ := nw_of_mem hP hγ
  rw [← hj]; exact hP.zero j

theorem node_height {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {T' : Set ℝ} {γ : ℝ}
    (hγ : γ ∈ nodeSet P T') : |γ| ≤ T := by
  obtain ⟨j, hj, _⟩ := nw_of_mem hP hγ
  rw [← hj]; exact hP.height j

theorem im_half_add (γ : ℝ) : ((1 / 2 : ℂ) + (γ : ℂ) * I).im = γ := by simp

/-- Nodes in (a, b] count at most Ncount a b. -/
theorem nodes_count_le {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) (T' : Set ℝ) (a b : ℝ) :
    ∑ γ ∈ (nodeSet P T').filter (· ∈ Ioc a b), nw P γ ≤ Ncount a b := by
  classical
  set S := (nodeSet P T').filter (· ∈ Ioc a b)
  have h1 : ∑ γ ∈ S, nw P γ ≤ ∑ γ ∈ S, mult ((1 / 2 : ℂ) + (γ : ℂ) * I) :=
    Finset.sum_le_sum fun γ hγ => nw_le_mult hP (Finset.mem_filter.mp hγ).1
  have hinj : Set.InjOn (fun γ : ℝ => (1 / 2 : ℂ) + (γ : ℂ) * I) S := by
    intro x _ y _ hxy
    have := congrArg Complex.im hxy
    simpa using this
  rw [← Finset.sum_image hinj] at h1
  refine h1.trans (sum_le_Ncount _ fun ρ hρ => ?_)
  obtain ⟨γ, hγ, rfl⟩ := Finset.mem_image.mp hρ
  obtain ⟨hγS, hγI⟩ := Finset.mem_filter.mp hγ
  exact ⟨node_mem_NTZ hP hγS, by rw [im_half_add]; exact hγI.1, by rw [im_half_add]; exact hγI.2⟩

/-- Nodes in (a, b] with b ≤ 0 count at most Ncount (2a - b) (-a) (mirror image). -/
theorem nodes_count_le_neg {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) (T' : Set ℝ) (a b : ℝ)
    (hab : a < b) :
    ∑ γ ∈ (nodeSet P T').filter (· ∈ Ioc a b), nw P γ ≤ Ncount (-b - (b - a)) (-a) := by
  classical
  set S := (nodeSet P T').filter (· ∈ Ioc a b)
  have h1 : ∑ γ ∈ S, nw P γ ≤ ∑ γ ∈ S, mult (conj ((1 / 2 : ℂ) + (γ : ℂ) * I)) := by
    refine Finset.sum_le_sum fun γ hγ => ?_
    rw [mult_conj _ (NTZ_ne_one (node_mem_NTZ hP (Finset.mem_filter.mp hγ).1))]
    exact nw_le_mult hP (Finset.mem_filter.mp hγ).1
  have hinj : Set.InjOn (fun γ : ℝ => conj ((1 / 2 : ℂ) + (γ : ℂ) * I)) S := by
    intro x _ y _ hxy
    have := congrArg Complex.im hxy
    simpa using this
  rw [← Finset.sum_image hinj] at h1
  refine h1.trans (sum_le_Ncount _ fun ρ hρ => ?_)
  obtain ⟨γ, hγ, rfl⟩ := Finset.mem_image.mp hρ
  obtain ⟨hγS, hγI⟩ := Finset.mem_filter.mp hγ
  refine ⟨conj_mem_NTZ (node_mem_NTZ hP hγS), ?_, ?_⟩
  · simp only [Complex.conj_im, im_half_add]; linarith [hγI.2]
  · simp only [Complex.conj_im, im_half_add]; linarith [hγI.1]

theorem gMin_eq_min' {g : Finset ℝ} (hg : g.Nonempty) : gMin g = g.min' hg :=
  Finset.Nonempty.csInf_eq_min' hg

/-- A group of nodes all of size ≥ H > s₀/2 lies on one side. -/
theorem groups_one_side (P : ExpSum) (T : Set ℝ) {s₀ H : ℝ} (hs₀ : 0 < s₀) (hH : s₀ < 2 * H)
    (hTH : ∀ x ∈ nodeSet P T, H ≤ |x|) {g : Finset ℝ} (hg : g ∈ groups P T s₀) :
    H ≤ gMin g ∨ gMax g ≤ -H := by
  obtain ⟨x, hx, rfl⟩ := mem_groups.mp hg
  set g := grp (nodeSet P T) s₀ x
  have hne : g.Nonempty := ⟨x, self_mem_grp hx⟩
  have hmin := gMin_mem hne
  have hmax := gMax_mem hne
  have hHa := hTH _ (grp_subset x hmin)
  have hHb := hTH _ (grp_subset x hmax)
  by_cases hpos : 0 < gMin g
  · left; rw [abs_of_pos hpos] at hHa; exact hHa
  · push Not at hpos
    by_cases hneg : gMax g < 0
    · right; rw [abs_of_neg hneg] at hHb; linarith
    · push Not at hneg
      exfalso
      have h0 : (0 : ℝ) ∈ comp (nodeSet P T) s₀ x := Icc_subset_comp hs₀ hmin hmax ⟨hpos, hneg⟩
      have := comp_subset_balls hs₀ h0
      simp only [mem_iUnion] at this
      obtain ⟨γ, hγ, hγ0⟩ := this
      have hγH := hTH γ (grp_subset x hγ)
      have : |γ| < s₀ / 2 := abs_lt.mpr ⟨by linarith [hγ0.2], by linarith [hγ0.1]⟩
      linarith

/-- **RH bounds for the high groups** of a zero configuration. -/
theorem hi_group_facts {P : ExpSum} {T : ℝ} (hP : ZeroConfig P T) {H s₀ ε₁ t₁ : ℝ} {Ku : ℕ}
    (hs₀ : 0 < s₀) (hε₁ : 0 < ε₁) (hε₁1 : ε₁ ≤ 1) (hH : t₁ + 2 ≤ H) (ht₁ : 1 ≤ t₁)
    (hcl : ∀ t ≥ t₁, ∀ h, 0 < h → h ≤ ε₁ → (Ncount t (t + h) : ℝ) ≤ ε₁ * Real.log t)
    (hKu : ((Ku : ℝ) + 2) * s₀ ≤ ε₁ / 2) (hKbig : ε₁ * Real.log (T + 2) ≤ Ku) :
    ∀ g ∈ groups P {γ | H < |γ|} s₀,
      1 ≤ (P.part ↑g (gC g)).M ∧ ((P.part ↑g (gC g)).M : ℝ) ≤ ε₁ * Real.log (|gC g| + 2) ∧
      (P.part ↑g (gC g)).M ≤ Ku ∧ gW g ≤ (P.part ↑g (gC g)).M * s₀ ∧ H ≤ |gC g| := by
  classical
  intro g hg
  set Thi : Set ℝ := {γ | H < |γ|}
  set S := nodeSet P Thi
  have hTH : ∀ x ∈ S, H ≤ |x| := fun x hx => (mem_nodeSet.mp hx).2.le
  have hs₀1 : s₀ ≤ 1 / 4 := by
    have : 2 * s₀ ≤ ((Ku : ℝ) + 2) * s₀ := by
      have : (2 : ℝ) ≤ (Ku : ℝ) + 2 := by have := Nat.cast_nonneg (α := ℝ) Ku; linarith
      exact mul_le_mul_of_nonneg_right this hs₀.le
    linarith
  have hsH : s₀ < 2 * H := by linarith
  obtain ⟨x, hx, hgx⟩ := mem_groups.mp hg
  have hne : g.Nonempty := groups_nonempty hg
  set M := (P.part ↑g (gC g)).M
  have hM : M = ∑ γ ∈ g, nw P γ := part_M_eq P g _
  have hnw : ∀ γ ∈ S, 1 ≤ nw P γ := fun γ hγ => nw_pos hP hγ
  have hcard : g.card ≤ M := by
    rw [hM, Finset.card_eq_sum_ones]
    exact Finset.sum_le_sum fun γ hγ => hnw γ (groups_subset hg hγ)
  set h := ((Ku : ℝ) + 2) * s₀
  have hh0 : 0 < h := by positivity
  have hhε : h ≤ ε₁ / 2 := hKu
  have hmin : gMin g ∈ g := gMin_mem hne
  have hmin_height : |gMin g| ≤ T := node_height hP (groups_subset hg hmin)
  set t' := gMin g - s₀ / 2
  have hmin' : (grp S s₀ x).min' ⟨x, self_mem_grp hx⟩ = gMin g := by
    rw [gMin_eq_min' hne]; congr 1 <;> exact hgx.symm
  -- the multiplicity bound with a height-dependent K
  have hside := groups_one_side P Thi hs₀ hsH hTH hg
  have key : ∃ τ : ℝ, t₁ ≤ τ ∧ τ ≤ T + 2 ∧ (τ ≤ |gC g| + 1 ∨ τ ≤ |gMin g| + 1) ∧
      ∑ γ ∈ S.filter (· ∈ Ioc t' (t' + h)), nw P γ ≤ ⌊ε₁ * Real.log τ⌋₊ ∧
      (H ≤ gMin g → τ ≤ |gC g| + 1) := by
    rcases hside with hpos | hneg
    · refine ⟨t', by simp only [t']; linarith, ?_, ?_, ?_, fun _ => ?_⟩
      · have := le_abs_self (gMin g); simp only [t']; linarith
      · left
        have h1 : gMin g ≤ gC g := by unfold gC; linarith [gMin_le (gMax_mem hne)]
        have h2 : gC g ≤ |gC g| := le_abs_self _
        simp only [t']; linarith
      · apply Nat.le_floor
        have hc := nodes_count_le hP Thi t' (t' + h)
        have hb := hcl t' (by simp only [t']; linarith) h hh0 (by linarith)
        calc ((∑ γ ∈ S.filter (· ∈ Ioc t' (t' + h)), nw P γ : ℕ) : ℝ) ≤ (Ncount t' (t' + h) : ℝ) := by
              exact_mod_cast hc
          _ ≤ ε₁ * Real.log t' := hb
      · have h1 : gMin g ≤ gC g := by unfold gC; linarith [gMin_le (gMax_mem hne)]
        have h2 : gC g ≤ |gC g| := le_abs_self _
        simp only [t']; linarith
    · have hmaxv : gMin g ≤ -H := (gMin_le (gMax_mem hne)).trans hneg
      set τ := -t' - 2 * h
      have hτ : t₁ ≤ τ := by simp only [τ, t']; linarith
      refine ⟨τ, hτ, ?_, ?_, ?_, fun hH' => absurd hH' (by linarith)⟩
      · have := neg_abs_le (gMin g); simp only [τ, t']; linarith
      · right
        have : -gMin g ≤ |gMin g| := neg_le_abs _
        simp only [τ, t']; linarith
      · apply Nat.le_floor
        have hc := nodes_count_le_neg hP Thi t' (t' + h) (by linarith)
        have e : -(t' + h) - (t' + h - t') = τ := by simp only [τ]; ring
        have e2 : -t' = τ + 2 * h := by simp only [τ]; ring
        rw [e, e2] at hc
        have hb := hcl τ hτ (2 * h) (by positivity) (by linarith)
        calc ((∑ γ ∈ S.filter (· ∈ Ioc t' (t' + h)), nw P γ : ℕ) : ℝ) ≤ (Ncount τ (τ + 2 * h) : ℝ) := by
              exact_mod_cast hc
          _ ≤ ε₁ * Real.log τ := hb
  obtain ⟨τ, hτ1, hτT, hτc, hcount, hτpos⟩ := key
  have hτpos' : 0 < τ := by linarith
  have hKτ : ⌊ε₁ * Real.log τ⌋₊ ≤ Ku := by
    apply Nat.floor_le_of_le
    calc ε₁ * Real.log τ ≤ ε₁ * Real.log (T + 2) :=
          mul_le_mul_of_nonneg_left (Real.log_le_log hτpos' hτT) hε₁.le
      _ ≤ Ku := hKbig
  have hKh : ((⌊ε₁ * Real.log τ⌋₊ : ℕ) + 2 : ℝ) * s₀ ≤ h := by
    have : ((⌊ε₁ * Real.log τ⌋₊ : ℕ) : ℝ) ≤ Ku := by exact_mod_cast hKτ
    simp only [h]; gcongr
  have hgm := grp_mult_le hs₀ (nw P) hnw hx hKh (by rw [hmin']; exact hcount)
  rw [hgx] at hgm
  have hMK : M ≤ ⌊ε₁ * Real.log τ⌋₊ := by rw [hM]; exact hgm
  have hMKu : M ≤ Ku := hMK.trans hKτ
  -- width
  have hwidth : gW g ≤ M * s₀ := by
    have hw := grp_width hs₀ (x := x) (by rw [hgx]; exact hmin) (by rw [hgx]; exact gMax_mem hne)
      (gMin_le (gMax_mem hne))
    rw [hgx] at hw
    have hc : (g.card : ℝ) ≤ M := by exact_mod_cast hcard
    unfold gW
    have : (g.card : ℝ) * s₀ ≤ M * s₀ := mul_le_mul_of_nonneg_right hc hs₀.le
    have hM0 : (0 : ℝ) ≤ M * s₀ := by positivity
    linarith
  have hcenter : H ≤ |gC g| := groups_center P Thi hs₀ hsH hTH hg
  have hW1 : gW g ≤ 1 := by
    have : (M : ℝ) * s₀ ≤ Ku * s₀ := by
      have : (M : ℝ) ≤ Ku := by exact_mod_cast hMKu
      exact mul_le_mul_of_nonneg_right this hs₀.le
    have : (Ku : ℝ) * s₀ ≤ ((Ku : ℝ) + 2) * s₀ := by nlinarith
    linarith
  refine ⟨?_, ?_, hMKu, hwidth, hcenter⟩
  · rw [hM]
    obtain ⟨y, hy⟩ := hne
    calc 1 ≤ nw P y := hnw y (groups_subset hg hy)
      _ ≤ ∑ γ ∈ g, nw P γ := Finset.single_le_sum (f := nw P) (fun _ _ => Nat.zero_le _) hy
  · -- height-dependent bound
    have hτc' : τ ≤ |gC g| + 2 := by
      rcases hτc with h1 | h1
      · linarith
      · have : |gMin g| ≤ |gC g| + gW g := by
          have e : gMin g = gC g - gW g := by unfold gC gW; ring
          rw [e]
          calc |gC g - gW g| ≤ |gC g| + |gW g| := abs_sub _ _
            _ = |gC g| + gW g := by rw [abs_of_nonneg (gW_nonneg hne)]
        linarith
    calc (M : ℝ) ≤ ⌊ε₁ * Real.log τ⌋₊ := by exact_mod_cast hMK
      _ ≤ ε₁ * Real.log τ := Nat.floor_le (mul_nonneg hε₁.le (Real.log_nonneg (by linarith)))
      _ ≤ ε₁ * Real.log (|gC g| + 2) :=
          mul_le_mul_of_nonneg_left (Real.log_le_log hτpos' hτc') hε₁.le

end Riemann.Deep.LeanNBRate


