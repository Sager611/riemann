import Riemann.Deep.LeanNBRate.ExpSumOps
import Riemann.Deep.LeanNBRate.Groups

/-!
# Configurations: parts of an exponential sum, groups of its nodes
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

namespace ExpSum

open Classical in
/-- The terms of P with nodes in T, nodes shifted by -c (other terms removed: multiplicity 0). -/
def part (P : ExpSum) (T : Set ℝ) (c : ℝ) : ExpSum where
  J := P.J
  node j := if P.node j ∈ T then P.node j - c else 0
  mult j := if P.node j ∈ T then P.mult j else 0
  coef := P.coef

variable (P : ExpSum)

open Classical in
theorem part_eval (T : Set ℝ) (c u : ℝ) : (P.part T c).eval u =
    ∑ j, if P.node j ∈ T then ∑ k ∈ Finset.range (P.mult j),
      P.coef j k * (I * u) ^ k * cexp (I * ((P.node j - c : ℝ) : ℂ) * u) else 0 := by
  unfold eval part
  refine Finset.sum_congr rfl fun j _ => ?_
  by_cases h : P.node j ∈ T
  · simp only [h, if_true]
  · simp only [h, if_false, Finset.range_zero, Finset.sum_empty]

theorem part_eval_shift (T : Set ℝ) (c u : ℝ) :
    cexp (I * c * u) * (P.part T c).eval u = (P.part T 0).eval u := by
  classical
  rw [part_eval, part_eval, Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  by_cases h : P.node j ∈ T
  · simp only [h, if_true, Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    have : cexp (I * c * u) * cexp (I * ((P.node j - c : ℝ) : ℂ) * u) =
        cexp (I * ((P.node j - 0 : ℝ) : ℂ) * u) := by
      rw [← Complex.exp_add]; push_cast; ring_nf
    calc cexp (I * c * u) * (P.coef j k * (I * u) ^ k * cexp (I * ((P.node j - c : ℝ) : ℂ) * u))
        = P.coef j k * (I * u) ^ k * (cexp (I * c * u) * cexp (I * ((P.node j - c : ℝ) : ℂ) * u)) := by
          ring
      _ = _ := by rw [this]
  · simp only [h, if_false, mul_zero]

theorem part_eval_zero_eq (T : Set ℝ) (c : ℝ) : (P.part T c).eval 0 = (P.part T 0).eval 0 := by
  have := P.part_eval_shift T c 0
  simpa using this

open Classical in
theorem part_M (T : Set ℝ) (c : ℝ) :
    (P.part T c).M = ∑ j, if P.node j ∈ T then P.mult j else 0 := rfl

theorem part_node_mem (T : Set ℝ) {c w : ℝ} (hw : 0 ≤ w)
    (h : ∀ j, P.node j ∈ T → |P.node j - c| ≤ w) : ∀ j, (P.part T c).node j ∈ Icc (-w) w := by
  classical
  intro j
  show (if P.node j ∈ T then P.node j - c else 0) ∈ Icc (-w) w
  split_ifs with hj
  · exact abs_le.mp (h j hj)
  · exact ⟨by linarith, hw⟩

theorem part_univ_eval : (P.part univ 0).eval = P.eval := by
  classical
  ext u
  rw [part_eval]
  unfold eval
  refine Finset.sum_congr rfl fun j _ => ?_
  simp only [mem_univ, if_true, sub_zero]

/-- Splitting along a partition of the nodes in T. -/
theorem part_eval_sum (T : Set ℝ) (𝒢 : Finset (Finset ℝ))
    (hcov : ∀ j, P.node j ∈ T → ∃ g ∈ 𝒢, P.node j ∈ g)
    (hsub : ∀ g ∈ 𝒢, ∀ x ∈ g, x ∈ T)
    (hdisj : ∀ g ∈ 𝒢, ∀ g' ∈ 𝒢, ∀ x, x ∈ g → x ∈ g' → g = g') (u : ℝ) :
    ∑ g ∈ 𝒢, (P.part ↑g 0).eval u = (P.part T 0).eval u := by
  classical
  simp_rw [part_eval]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  by_cases hj : P.node j ∈ T
  · obtain ⟨g₀, hg₀, hjg₀⟩ := hcov j hj
    rw [if_pos hj, Finset.sum_eq_single g₀]
    · rw [if_pos (by exact_mod_cast hjg₀)]
    · intro g hg hne
      rw [if_neg]
      intro hjg
      exact hne (hdisj g hg g₀ hg₀ _ (by exact_mod_cast hjg) hjg₀)
    · intro h; exact absurd hg₀ h
  · rw [if_neg hj]
    refine Finset.sum_eq_zero fun g hg => ?_
    rw [if_neg]
    intro hjg
    exact hj (hsub g hg _ (by exact_mod_cast hjg))

end ExpSum

/-! ## Node sets and groups -/

/-- The nodes of P lying in T. -/
def nodeSet (P : ExpSum) (T : Set ℝ) : Finset ℝ := by
  classical exact (Finset.univ.image P.node).filter (· ∈ T)

theorem mem_nodeSet {P : ExpSum} {T : Set ℝ} {x : ℝ} :
    x ∈ nodeSet P T ↔ (∃ j, P.node j = x) ∧ x ∈ T := by
  classical
  unfold nodeSet
  simp [Finset.mem_filter, Finset.mem_image]

/-- The groups of the nodes of P in T at scale s₀. -/
def groups (P : ExpSum) (T : Set ℝ) (s₀ : ℝ) : Finset (Finset ℝ) := by
  classical exact (nodeSet P T).image (grp (nodeSet P T) s₀)

theorem mem_groups {P : ExpSum} {T : Set ℝ} {s₀ : ℝ} {g : Finset ℝ} :
    g ∈ groups P T s₀ ↔ ∃ x ∈ nodeSet P T, grp (nodeSet P T) s₀ x = g := by
  classical
  unfold groups
  rw [Finset.mem_image]

/-- Window of a group: min, max, centre, half-width. -/
def gMin (g : Finset ℝ) : ℝ := sInf (g : Set ℝ)
def gMax (g : Finset ℝ) : ℝ := sSup (g : Set ℝ)
def gC (g : Finset ℝ) : ℝ := (gMin g + gMax g) / 2
def gW (g : Finset ℝ) : ℝ := (gMax g - gMin g) / 2

theorem gMin_le {g : Finset ℝ} {x : ℝ} (hx : x ∈ g) : gMin g ≤ x :=
  csInf_le g.finite_toSet.bddBelow hx

theorem le_gMax {g : Finset ℝ} {x : ℝ} (hx : x ∈ g) : x ≤ gMax g :=
  le_csSup g.finite_toSet.bddAbove hx

theorem gMin_mem {g : Finset ℝ} (hg : g.Nonempty) : gMin g ∈ g := by
  have := Set.Nonempty.csInf_mem (s := (g : Set ℝ)) (by exact_mod_cast hg) g.finite_toSet
  exact_mod_cast this

theorem gMax_mem {g : Finset ℝ} (hg : g.Nonempty) : gMax g ∈ g := by
  have := Set.Nonempty.csSup_mem (s := (g : Set ℝ)) (by exact_mod_cast hg) g.finite_toSet
  exact_mod_cast this

theorem abs_sub_gC_le {g : Finset ℝ} {x : ℝ} (hx : x ∈ g) : |x - gC g| ≤ gW g := by
  have h1 := gMin_le hx
  have h2 := le_gMax hx
  unfold gC gW
  rw [abs_le]; constructor <;> linarith

theorem gW_nonneg {g : Finset ℝ} (hg : g.Nonempty) : 0 ≤ gW g := by
  obtain ⟨x, hx⟩ := hg
  unfold gW; linarith [gMin_le hx, le_gMax hx]

section groups

variable {P : ExpSum} {T : Set ℝ} {s₀ : ℝ}

theorem groups_nonempty {g : Finset ℝ} (hg : g ∈ groups P T s₀) : g.Nonempty := by
  obtain ⟨x, hx, rfl⟩ := mem_groups.mp hg
  exact ⟨x, self_mem_grp hx⟩

theorem groups_subset {g : Finset ℝ} (hg : g ∈ groups P T s₀) : g ⊆ nodeSet P T := by
  obtain ⟨x, _, rfl⟩ := mem_groups.mp hg
  exact grp_subset x

theorem groups_cover {x : ℝ} (hx : x ∈ nodeSet P T) : ∃ g ∈ groups P T s₀, x ∈ g :=
  ⟨grp (nodeSet P T) s₀ x, mem_groups.mpr ⟨x, hx, rfl⟩, self_mem_grp hx⟩

theorem groups_disj {g g' : Finset ℝ} (hg : g ∈ groups P T s₀) (hg' : g' ∈ groups P T s₀)
    {x : ℝ} (hxg : x ∈ g) (hxg' : x ∈ g') : g = g' := by
  obtain ⟨y, _, rfl⟩ := mem_groups.mp hg
  obtain ⟨y', _, rfl⟩ := mem_groups.mp hg'
  rw [← grp_eq_of_mem hxg, ← grp_eq_of_mem hxg']

/-- The windows of distinct groups are s₀-separated. -/
theorem groups_sep (hs₀ : 0 < s₀) {g g' : Finset ℝ} (hg : g ∈ groups P T s₀)
    (hg' : g' ∈ groups P T s₀) (hne : g ≠ g') : gW g + gW g' + s₀ ≤ |gC g - gC g'| := by
  obtain ⟨x, hx, rfl⟩ := mem_groups.mp hg
  obtain ⟨y, hy, rfl⟩ := mem_groups.mp hg'
  have hcne : comp (nodeSet P T) s₀ x ≠ comp (nodeSet P T) s₀ y := by
    intro h; apply hne; ext z; rw [mem_grp, mem_grp, h]
  have hxy : x ≠ y := by intro h; exact hcne (h ▸ rfl)
  have hgne := groups_nonempty hg
  have hgne' := groups_nonempty hg'
  have hab := gMin_le (gMax_mem hgne)
  have hab' := gMin_le (gMax_mem hgne')
  rcases lt_or_gt_of_ne hxy with hlt | hlt
  · have hsep := grp_sep hs₀ hx hy hcne hlt (gMax_mem hgne) (gMin_mem hgne')
    unfold gC gW
    have h0 : 0 ≤ (gMin (grp (nodeSet P T) s₀ y) + gMax (grp (nodeSet P T) s₀ y)) / 2 -
        (gMin (grp (nodeSet P T) s₀ x) + gMax (grp (nodeSet P T) s₀ x)) / 2 := by linarith
    rw [abs_sub_comm, abs_of_nonneg h0]
    linarith
  · have hsep := grp_sep hs₀ hy hx (Ne.symm hcne) hlt (gMax_mem hgne') (gMin_mem hgne)
    unfold gC gW
    have h0 : 0 ≤ (gMin (grp (nodeSet P T) s₀ x) + gMax (grp (nodeSet P T) s₀ x)) / 2 -
        (gMin (grp (nodeSet P T) s₀ y) + gMax (grp (nodeSet P T) s₀ y)) / 2 := by linarith
    rw [abs_of_nonneg h0]
    linarith

theorem part_eval_groups (s₀ : ℝ) (u : ℝ) :
    ∑ g ∈ groups P T s₀, (P.part ↑g 0).eval u = (P.part T 0).eval u := by
  refine P.part_eval_sum T _ (fun j hj => ?_) (fun g hg x hx => ?_) (fun g hg g' hg' x hx hx' => ?_) u
  · exact groups_cover (mem_nodeSet.mpr ⟨⟨j, rfl⟩, hj⟩)
  · exact (mem_nodeSet.mp (groups_subset hg hx)).2
  · exact groups_disj hg hg' hx hx'

end groups

end Riemann.Deep.LeanNBRate

