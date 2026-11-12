module

public import Mathlib.Algebra.Order.BigOperators.Ring.Finset
public import Mathlib.Basic.Complex.Basic
public import Mathlib.Basic.Real.Basic
public import Mathlib.Tactic.Linarith

@[expose] public section

/-!
# Theorems 1 and 1′(a): zero counts, the per-window count, the class split

Zero data of one family (NOTE §0, l.40–41): a finite set X of characters (X = Ĝ) with the trivial
one; win χ t = the distinct zeros ρ of L(s, χ*) with |Im ρ − t| ≤ 1; mult = multiplicity.
N(σ; t) = Σ_χ #{ρ : β ≥ σ, |γ − t| ≤ 1} with multiplicity (NOTE l.40).
Proved here: N(σ;t) ≤ B (|𝒳| + 1) when every window holds ≤ B zeros (NOTE l.162–164,
RCApps2 l.57), and |𝒳| ≤ B_I + B_II from a class I / class II split of one chosen zero per χ
(NOTE l.164–165, RCApps2 l.37–39).
-/

open Finset

namespace RingClass

/-- zero data of one family -/
structure ZeroData (X : Type*) where
  /-- the trivial character -/
  one : X
  /-- the distinct zeros ρ of L(s, χ*) (0 < Re ρ < 1) with |Im ρ − t| ≤ 1 -/
  win : X → ℝ → Finset ℂ
  /-- multiplicity of the zero ρ of L(s, χ*) -/
  mult : X → ℂ → ℕ
  win_im : ∀ χ t ρ, ρ ∈ win χ t → |ρ.im - t| ≤ 1

variable {X : Type*}

/-- N(σ; t) = Σ_χ #{ρ = β + iγ : L(ρ, χ*) = 0, β ≥ σ, |γ − t| ≤ 1}, with multiplicity -/
noncomputable def ZeroData.N [Fintype X] (z : ZeroData X) (σ t : ℝ) : ℕ :=
  ∑ χ, ∑ ρ ∈ (z.win χ t).filter (fun ρ => σ ≤ ρ.re), z.mult χ ρ

/-- 𝒳 = the χ ≠ 1 with a zero in [σ, 1] × [t − 1, t + 1] -/
noncomputable def ZeroData.Xs [Fintype X] (z : ZeroData X) (σ t : ℝ) : Finset X := by
  classical
  exact Finset.univ.filter (fun χ => χ ≠ z.one ∧ ∃ ρ ∈ z.win χ t, σ ≤ ρ.re)

lemma ZeroData.mem_Xs [Fintype X] (z : ZeroData X) (σ t : ℝ) (χ : X) :
    χ ∈ z.Xs σ t ↔ χ ≠ z.one ∧ ∃ ρ ∈ z.win χ t, σ ≤ ρ.re := by
  classical
  unfold ZeroData.Xs
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]

/-- per-window count: N(σ;t) ≤ B (|𝒳| + 1) if every window holds ≤ B zeros -/
lemma ZeroData.N_le [Fintype X] [DecidableEq X] (z : ZeroData X) (σ t B : ℝ)
    (hB : ∀ χ, (∑ ρ ∈ z.win χ t, (z.mult χ ρ : ℝ)) ≤ B) :
    (z.N σ t : ℝ) ≤ B * ((z.Xs σ t).card + 1) := by
  classical
  have hB0 : 0 ≤ B := le_trans (Finset.sum_nonneg (fun _ _ => Nat.cast_nonneg _)) (hB z.one)
  have hsub : ∀ χ, (∑ ρ ∈ (z.win χ t).filter (fun ρ => σ ≤ ρ.re), (z.mult χ ρ : ℝ))
      ≤ if χ ∈ insert z.one (z.Xs σ t) then B else 0 := by
    intro χ
    split_ifs with h
    · exact le_trans (Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
        (fun _ _ _ => Nat.cast_nonneg _)) (hB χ)
    · have h1 : (z.win χ t).filter (fun ρ => σ ≤ ρ.re) = ∅ := by
        rw [Finset.filter_eq_empty_iff]
        intro ρ hρ hσ
        apply h
        rw [Finset.mem_insert]
        right
        rw [z.mem_Xs]
        refine ⟨?_, ρ, hρ, hσ⟩
        intro h2; apply h; rw [h2]; exact Finset.mem_insert_self _ _
      rw [h1, Finset.sum_empty]
  have e1 : (z.N σ t : ℝ) = ∑ χ, ∑ ρ ∈ (z.win χ t).filter (fun ρ => σ ≤ ρ.re), (z.mult χ ρ : ℝ) := by
    unfold ZeroData.N; push_cast; rfl
  rw [e1]
  calc (∑ χ, ∑ ρ ∈ (z.win χ t).filter (fun ρ => σ ≤ ρ.re), (z.mult χ ρ : ℝ))
      ≤ ∑ χ, (if χ ∈ insert z.one (z.Xs σ t) then B else 0) :=
        Finset.sum_le_sum (fun χ _ => hsub χ)
    _ = B * (insert z.one (z.Xs σ t)).card := by
        rw [Finset.sum_ite_mem, Finset.univ_inter, Finset.sum_const, nsmul_eq_mul, mul_comm]
    _ ≤ B * ((z.Xs σ t).card + 1) := by
        apply mul_le_mul_of_nonneg_left _ hB0
        exact_mod_cast Finset.card_insert_le _ _

/-- one chosen zero ρ_χ per χ ∈ 𝒳 (NOTE l.162–163) -/
noncomputable def ZeroData.sel (z : ZeroData X) (σ t : ℝ) (χ : X) : ℂ := by
  classical
  exact if h : ∃ ρ ∈ z.win χ t, σ ≤ ρ.re then Classical.choose h else 0

lemma ZeroData.sel_spec (z : ZeroData X) (σ t : ℝ) (χ : X) (h : ∃ ρ ∈ z.win χ t, σ ≤ ρ.re) :
    z.sel σ t χ ∈ z.win χ t ∧ σ ≤ (z.sel σ t χ).re := by
  classical
  unfold ZeroData.sel
  split_ifs
  exact Classical.choose_spec h

/-- the class split (Lemma 5): if every zero of a χ ≠ 1 with β ≥ σ is of class I or class II, and
every set of characters with selected class-I (class-II) zeros has ≤ B_I (≤ B_II) elements, then
|𝒳| ≤ B_I + B_II -/
lemma ZeroData.Xs_card_le [Fintype X] (z : ZeroData X) (σ t : ℝ) (PI PII : X → ℂ → Prop)
    (hsplit : ∀ χ, χ ≠ z.one → ∀ ρ ∈ z.win χ t, σ ≤ ρ.re → PI χ ρ ∨ PII χ ρ) (BI BII : ℝ)
    (hI : ∀ (S : Finset X) (ρ : X → ℂ), (∀ χ ∈ S, χ ≠ z.one ∧ ρ χ ∈ z.win χ t ∧ σ ≤ (ρ χ).re ∧
      PI χ (ρ χ)) → (S.card : ℝ) ≤ BI)
    (hII : ∀ (S : Finset X) (ρ : X → ℂ), (∀ χ ∈ S, χ ≠ z.one ∧ ρ χ ∈ z.win χ t ∧ σ ≤ (ρ χ).re ∧
      PII χ (ρ χ)) → (S.card : ℝ) ≤ BII) :
    ((z.Xs σ t).card : ℝ) ≤ BI + BII := by
  classical
  have hsel : ∀ χ ∈ z.Xs σ t, χ ≠ z.one ∧ z.sel σ t χ ∈ z.win χ t ∧ σ ≤ (z.sel σ t χ).re := by
    intro χ hχ
    rw [z.mem_Xs] at hχ
    exact ⟨hχ.1, z.sel_spec σ t χ hχ.2⟩
  set SI := (z.Xs σ t).filter (fun χ => PI χ (z.sel σ t χ)) with hSI
  set SII := (z.Xs σ t).filter (fun χ => PII χ (z.sel σ t χ)) with hSII
  have hcover : z.Xs σ t ⊆ SI ∪ SII := by
    intro χ hχ
    obtain ⟨h1, h2, h3⟩ := hsel χ hχ
    rw [Finset.mem_union, hSI, hSII, Finset.mem_filter, Finset.mem_filter]
    rcases hsplit χ h1 _ h2 h3 with h | h
    · exact Or.inl ⟨hχ, h⟩
    · exact Or.inr ⟨hχ, h⟩
  have hI' := hI SI (z.sel σ t) (fun χ hχ => by
    rw [hSI, Finset.mem_filter] at hχ
    obtain ⟨h1, h2, h3⟩ := hsel χ hχ.1
    exact ⟨h1, h2, h3, hχ.2⟩)
  have hII' := hII SII (z.sel σ t) (fun χ hχ => by
    rw [hSII, Finset.mem_filter] at hχ
    obtain ⟨h1, h2, h3⟩ := hsel χ hχ.1
    exact ⟨h1, h2, h3, hχ.2⟩)
  have h4 : (z.Xs σ t).card ≤ SI.card + SII.card :=
    (Finset.card_le_card hcover).trans (Finset.card_union_le _ _)
  have h5 : ((z.Xs σ t).card : ℝ) ≤ SI.card + SII.card := by exact_mod_cast h4
  linarith

end RingClass

