import Riemann.Deep.DH.FE

/-!
# Round 81 — filters as Lean theorems

Already proved earlier: σ-blindness (Fayed.R1.fayed_profile_realizable), finite base + tail (Fayed.R45.tail_iff),
symmetry ≠ fixedness (Fayed.R46.zeta_offline_pair), positivity ≡ RH (Fayed.WeilSuff.weil_criterion_weilPos).
New below: the Davenport–Heilbronn filter (from the cited Spira zero), local criteria = zero-free regions,
and a model showing that verification up to any height is consistent with an off-line zero.
Deferred (costly): Beurling counter-models (Diamond–Montgomery–Vorhauer), explicit zero-free region constants.
-/

namespace Fayed.R81

/-- Davenport–Heilbronn filter: the cited Spira zero lies off the critical line, so no argument that
uses only properties shared by DH (Dirichlet series + a ζ-type functional equation) can prove RH. -/
theorem dh_offline (h : Fayed.DH.OffLineZeroClaim) : ∃ s : ℂ, Fayed.DH.DH s = 0 ∧ s.re ≠ 1 / 2 := by
  obtain ⟨s, hs, hre, _⟩ := h
  refine ⟨s, hs, fun h12 => ?_⟩
  rw [h12] at hre
  norm_num [abs_lt] at hre

/-- Local criteria are zero-free regions in disguise: a criterion failing at every zero,
once proved on a region, makes the region zero-free. -/
theorem zero_free_of_criterion {α : Type*} (IsZero C : α → Prop) (R : Set α)
    (hf : ∀ z, IsZero z → ¬ C z) (hr : ∀ z ∈ R, C z) : ∀ z, IsZero z → z ∉ R :=
  fun z hz hzR => hf z hz (hr z hzR)

/-- Finite verification is consistent with failure: for every height T there is a zero set that is
on the line up to height T and still contains an off-line point. -/
theorem finite_check_consistent (T : ℝ) :
    ∃ Z : Set ℂ, (∀ s ∈ Z, |s.im| ≤ T → s.re = 1 / 2) ∧ ∃ s ∈ Z, s.re ≠ 1 / 2 := by
  refine ⟨{(⟨1 / 3, T + 1⟩ : ℂ)}, ?_, ⟨⟨1 / 3, T + 1⟩, rfl, by norm_num⟩⟩
  intro s hs hT
  rw [Set.mem_singleton_iff] at hs
  subst hs
  simp only at hT
  rw [abs_le] at hT
  linarith [hT.2]

end Fayed.R81

#print axioms Fayed.R81.dh_offline
