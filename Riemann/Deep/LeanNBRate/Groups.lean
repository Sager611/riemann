import Mathlib.Topology.Connected.LocallyConnected
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.MeasureTheory.Measure.Lebesgue.Basic
import Mathlib.Order.Interval.Set.Basic

/-!
# [N] §3 Grouping of ordinates at scale s₀

Groups are the nodes lying in one connected component of the s₀/2-neighbourhood of the node set.
Windows of distinct groups are s₀-separated; each window is covered by the s₀/2-balls of its own
nodes (so its length is ≤ #g·s₀); a group's multiplicity is bounded by a short-interval count.
-/

open Set MeasureTheory

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The s₀/2-neighbourhood of a finite set of ordinates. -/
def nbhd (S : Finset ℝ) (s₀ : ℝ) : Set ℝ := ⋃ γ ∈ S, Ioo (γ - s₀ / 2) (γ + s₀ / 2)

/-- The component of x in the neighbourhood. -/
def comp (S : Finset ℝ) (s₀ x : ℝ) : Set ℝ := connectedComponentIn (nbhd S s₀) x

/-- The group of x: the nodes in the component of x. -/
def grp (S : Finset ℝ) (s₀ x : ℝ) : Finset ℝ := by
  classical exact S.filter fun y => comp S s₀ y = comp S s₀ x

variable {S : Finset ℝ} {s₀ : ℝ}

theorem ball_subset_nbhd {γ : ℝ} (hγ : γ ∈ S) : Ioo (γ - s₀ / 2) (γ + s₀ / 2) ⊆ nbhd S s₀ := by
  intro z hz
  simp only [nbhd, mem_iUnion]
  exact ⟨γ, hγ, hz⟩

theorem mem_nbhd (hs₀ : 0 < s₀) {γ : ℝ} (hγ : γ ∈ S) : γ ∈ nbhd S s₀ :=
  ball_subset_nbhd hγ ⟨by linarith, by linarith⟩

theorem ball_subset_comp (hs₀ : 0 < s₀) {γ : ℝ} (hγ : γ ∈ S) :
    Ioo (γ - s₀ / 2) (γ + s₀ / 2) ⊆ comp S s₀ γ :=
  isPreconnected_Ioo.subset_connectedComponentIn ⟨by linarith, by linarith⟩ (ball_subset_nbhd hγ)

theorem self_mem_comp (hs₀ : 0 < s₀) {γ : ℝ} (hγ : γ ∈ S) : γ ∈ comp S s₀ γ :=
  mem_connectedComponentIn (mem_nbhd hs₀ hγ)

theorem comp_eq_of_mem {x y : ℝ} (h : y ∈ comp S s₀ x) : comp S s₀ x = comp S s₀ y :=
  connectedComponentIn_eq h

theorem comp_eq_of_inter {x y z : ℝ} (hx : z ∈ comp S s₀ x) (hy : z ∈ comp S s₀ y) :
    comp S s₀ x = comp S s₀ y := (comp_eq_of_mem hx).trans (comp_eq_of_mem hy).symm

theorem mem_grp {x y : ℝ} : y ∈ grp S s₀ x ↔ y ∈ S ∧ comp S s₀ y = comp S s₀ x := by
  classical
  unfold grp
  rw [Finset.mem_filter]

theorem self_mem_grp {x : ℝ} (hx : x ∈ S) : x ∈ grp S s₀ x := mem_grp.mpr ⟨hx, rfl⟩

theorem grp_subset (x : ℝ) : grp S s₀ x ⊆ S := fun y hy => (mem_grp.mp hy).1

theorem grp_eq_of_mem {x y : ℝ} (hy : y ∈ grp S s₀ x) : grp S s₀ y = grp S s₀ x := by
  ext z
  rw [mem_grp, mem_grp, (mem_grp.mp hy).2]

theorem mem_comp_of_mem_grp (hs₀ : 0 < s₀) {x y : ℝ} (hy : y ∈ grp S s₀ x) : y ∈ comp S s₀ x := by
  obtain ⟨hyS, hc⟩ := mem_grp.mp hy
  rw [← hc]; exact self_mem_comp hs₀ hyS

/-- Nodes in different components are at distance ≥ s₀. -/
theorem sep_of_ne (hs₀ : 0 < s₀) {x y : ℝ} (hx : x ∈ S) (hy : y ∈ S)
    (hne : comp S s₀ x ≠ comp S s₀ y) : s₀ ≤ |x - y| := by
  by_contra h
  push Not at h
  apply hne
  rw [abs_lt] at h
  have h1 : (x + y) / 2 ∈ comp S s₀ x := ball_subset_comp hs₀ hx ⟨by linarith, by linarith⟩
  have h2 : (x + y) / 2 ∈ comp S s₀ y := ball_subset_comp hs₀ hy ⟨by linarith, by linarith⟩
  exact comp_eq_of_inter h1 h2

/-- Distinct components are ordered. -/
theorem comp_lt_of_ne (hs₀ : 0 < s₀) {x y : ℝ} (hx : x ∈ S) (hy : y ∈ S)
    (hne : comp S s₀ x ≠ comp S s₀ y) (hxy : x < y) {a b : ℝ} (ha : a ∈ comp S s₀ x)
    (hb : b ∈ comp S s₀ y) : a < b := by
  by_contra h
  push Not at h
  apply hne
  have hxc := self_mem_comp hs₀ hx
  have hyc := self_mem_comp hs₀ hy
  rcases le_or_gt x b with hxb | hbx
  · have : b ∈ comp S s₀ x :=
      (isPreconnected_connectedComponentIn.Icc_subset hxc ha) ⟨hxb, h⟩
    exact comp_eq_of_inter this hb
  · have : x ∈ comp S s₀ y :=
      (isPreconnected_connectedComponentIn.Icc_subset hb hyc) ⟨hbx.le, hxy.le⟩
    exact comp_eq_of_inter hxc this

/-- Nodes of distinct groups: the one in the lower group is ≥ s₀ below. -/
theorem grp_sep (hs₀ : 0 < s₀) {x y : ℝ} (hx : x ∈ S) (hy : y ∈ S)
    (hne : comp S s₀ x ≠ comp S s₀ y) (hxy : x < y) {a b : ℝ} (ha : a ∈ grp S s₀ x)
    (hb : b ∈ grp S s₀ y) : a + s₀ ≤ b := by
  have hab := comp_lt_of_ne hs₀ hx hy hne hxy (mem_comp_of_mem_grp hs₀ ha) (mem_comp_of_mem_grp hs₀ hb)
  have hne' : comp S s₀ a ≠ comp S s₀ b := by
    rw [(mem_grp.mp ha).2, (mem_grp.mp hb).2]; exact hne
  have := sep_of_ne hs₀ (grp_subset x ha) (grp_subset y hb) hne'
  rw [abs_sub_comm, abs_of_pos (by linarith)] at this
  linarith

/-- A component is covered by the s₀/2-balls of its group. -/
theorem comp_subset_balls (hs₀ : 0 < s₀) {x : ℝ} :
    comp S s₀ x ⊆ ⋃ γ ∈ grp S s₀ x, Ioo (γ - s₀ / 2) (γ + s₀ / 2) := by
  intro z hz
  have hzN : z ∈ nbhd S s₀ := connectedComponentIn_subset _ _ hz
  simp only [nbhd, mem_iUnion] at hzN
  obtain ⟨γ, hγ, hzγ⟩ := hzN
  simp only [mem_iUnion]
  refine ⟨γ, mem_grp.mpr ⟨hγ, ?_⟩, hzγ⟩
  exact (comp_eq_of_inter hz (ball_subset_comp hs₀ hγ hzγ)).symm

/-- Measure comparison: an interval covered by finitely many s₀/2-balls. -/
theorem length_le_card_mul {a b : ℝ} (hab : a ≤ b) (T : Finset ℝ)
    (hcov : Icc a b ⊆ ⋃ γ ∈ T, Ioo (γ - s₀ / 2) (γ + s₀ / 2)) (hs₀ : 0 < s₀) :
    b - a ≤ T.card * s₀ := by
  have h1 : volume (Icc a b) ≤ volume (⋃ γ ∈ T, Ioo (γ - s₀ / 2) (γ + s₀ / 2)) := measure_mono hcov
  have h2 : volume (⋃ γ ∈ T, Ioo (γ - s₀ / 2) (γ + s₀ / 2)) ≤
      ∑ γ ∈ T, volume (Ioo (γ - s₀ / 2) (γ + s₀ / 2)) := measure_biUnion_finset_le T _
  have h3 : ∑ γ ∈ T, volume (Ioo (γ - s₀ / 2) (γ + s₀ / 2)) = T.card * ENNReal.ofReal s₀ := by
    have e : ∀ γ ∈ T, volume (Ioo (γ - s₀ / 2) (γ + s₀ / 2)) = ENNReal.ofReal s₀ := by
      intro γ _; rw [Real.volume_Ioo]; congr 1; ring
    rw [Finset.sum_congr rfl e, Finset.sum_const, nsmul_eq_mul]
  rw [Real.volume_Icc] at h1
  have h4 : ENNReal.ofReal (b - a) ≤ ENNReal.ofReal (T.card * s₀) := by
    rw [ENNReal.ofReal_mul (by positivity), ENNReal.ofReal_natCast]
    exact h1.trans (h2.trans h3.le)
  exact (ENNReal.ofReal_le_ofReal_iff (by positivity)).mp h4

/-- The window of a group lies in its component. -/
theorem Icc_subset_comp (hs₀ : 0 < s₀) {x a b : ℝ} (ha : a ∈ grp S s₀ x) (hb : b ∈ grp S s₀ x) :
    Icc a b ⊆ comp S s₀ x :=
  isPreconnected_connectedComponentIn.Icc_subset (mem_comp_of_mem_grp hs₀ ha)
    (mem_comp_of_mem_grp hs₀ hb)

/-- Width: any two nodes of a group are within #g·s₀. -/
theorem grp_width (hs₀ : 0 < s₀) {x a b : ℝ} (ha : a ∈ grp S s₀ x) (hb : b ∈ grp S s₀ x)
    (hab : a ≤ b) : b - a ≤ (grp S s₀ x).card * s₀ :=
  length_le_card_mul hab _ ((Icc_subset_comp hs₀ ha hb).trans (comp_subset_balls hs₀)) hs₀

/-- Local count: a subinterval [a, b] of a group's window meets ≥ (b - a)/s₀ balls of its nodes. -/
theorem grp_local_count (hs₀ : 0 < s₀) {x a b a' b' : ℝ} (ha : a' ∈ grp S s₀ x)
    (hb : b' ∈ grp S s₀ x) (h1 : a' ≤ a) (h2 : b ≤ b') (hab : a ≤ b) :
    b - a ≤ ((grp S s₀ x).filter fun γ => γ ∈ Ioo (a - s₀ / 2) (b + s₀ / 2)).card * s₀ := by
  classical
  refine length_le_card_mul hab _ ?_ hs₀
  intro z hz
  have hz' : z ∈ comp S s₀ x := Icc_subset_comp hs₀ ha hb ⟨by linarith [hz.1], by linarith [hz.2]⟩
  have := comp_subset_balls hs₀ hz'
  simp only [mem_iUnion] at this ⊢
  obtain ⟨γ, hγ, hzγ⟩ := this
  refine ⟨γ, Finset.mem_filter.mpr ⟨hγ, ⟨by linarith [hz.1, hzγ.2], by linarith [hz.2, hzγ.1]⟩⟩, hzγ⟩

/-- **Multiplicity of a group from a short-interval count.** If the nodes carry multiplicities
m ≥ 1 and every node-count in (t, t + h] (t = min of the group − s₀/2) is ≤ K with (K + 2)s₀ ≤ h,
then the group's total multiplicity is ≤ K. -/
theorem grp_mult_le (hs₀ : 0 < s₀) (m : ℝ → ℕ) (hm : ∀ γ ∈ S, 1 ≤ m γ) {x : ℝ} (hx : x ∈ S)
    {K : ℕ} {h : ℝ} (hKh : ((K : ℝ) + 2) * s₀ ≤ h)
    (hcount : ∑ γ ∈ S.filter (fun γ => γ ∈ Ioc ((grp S s₀ x).min' ⟨x, self_mem_grp hx⟩ - s₀ / 2)
      ((grp S s₀ x).min' ⟨x, self_mem_grp hx⟩ - s₀ / 2 + h)), m γ ≤ K) :
    ∑ γ ∈ grp S s₀ x, m γ ≤ K := by
  classical
  set g := grp S s₀ x
  have hne : g.Nonempty := ⟨x, self_mem_grp hx⟩
  set a := g.min' hne
  set b := g.max' hne
  have ha : a ∈ g := g.min'_mem hne
  have hb : b ∈ g := g.max'_mem hne
  set t := a - s₀ / 2
  have hK0 : (0 : ℝ) ≤ K := by positivity
  rcases le_or_gt b (a + h - s₀) with hcase | hcase
  · -- the whole group lies in (t, t + h]
    have hsub : g ⊆ S.filter (fun γ => γ ∈ Ioc t (t + h)) := by
      intro γ hγ
      refine Finset.mem_filter.mpr ⟨grp_subset x hγ, ⟨?_, ?_⟩⟩
      · have := g.min'_le γ hγ; simp only [t]; linarith
      · have := g.le_max' γ hγ; simp only [t]; linarith
    exact (Finset.sum_le_sum_of_subset hsub).trans hcount
  · exfalso
    have hlen : a + ((K : ℝ) + 1) * s₀ ≤ b := by nlinarith
    have hc := grp_local_count hs₀ ha hb le_rfl hlen (by nlinarith)
    set F := g.filter fun γ => γ ∈ Ioo (a - s₀ / 2) (a + ((K : ℝ) + 1) * s₀ + s₀ / 2)
    have hF : ((K : ℝ) + 1) ≤ F.card := by
      have : (a + ((K : ℝ) + 1) * s₀ - a) = ((K : ℝ) + 1) * s₀ := by ring
      rw [this] at hc
      exact le_of_mul_le_mul_right hc hs₀
    have hsub : F ⊆ S.filter (fun γ => γ ∈ Ioc t (t + h)) := by
      intro γ hγ
      obtain ⟨hγg, hγI⟩ := Finset.mem_filter.mp hγ
      refine Finset.mem_filter.mpr ⟨grp_subset x hγg, ⟨?_, ?_⟩⟩
      · simp only [t]; linarith [hγI.1]
      · simp only [t]; nlinarith [hγI.2]
    have h1 : (F.card : ℝ) ≤ ∑ γ ∈ F, (m γ : ℝ) := by
      rw [Finset.card_eq_sum_ones, Nat.cast_sum]
      exact Finset.sum_le_sum fun γ hγ => by
        exact_mod_cast hm γ (grp_subset x (Finset.mem_filter.mp hγ).1)
    have h2 : ∑ γ ∈ F, (m γ : ℝ) ≤ ∑ γ ∈ S.filter (fun γ => γ ∈ Ioc t (t + h)), (m γ : ℝ) :=
      Finset.sum_le_sum_of_subset_of_nonneg hsub fun _ _ _ => by positivity
    have h3 : ∑ γ ∈ S.filter (fun γ => γ ∈ Ioc t (t + h)), (m γ : ℝ) ≤ K := by exact_mod_cast hcount
    linarith

end Riemann.Deep.LeanNBRate

