import Mathlib.LinearAlgebra.Projectivization.Cardinality
import Mathlib.Data.Fintype.Pi

/-!
# Lemma 7 (involution form), part A: collineations, fixed points, f² scalar

PG(d-1, F) is Mathlib's ℙ F (Fin d → F). A permutation g of Fin n is a collineation through a
bijection e : Fin n ≃ ℙ F (Fin d → F) (IsCollineationVia) if some σ-semilinear bijection f of F^d
(an element of ΓL_d(F)) maps each line e x into the line e (g x); these are exactly the elements
of PΓL_d(F) acting on the points. If g is an involution and d ≥ 2, then f² = μ·id and σ² = id.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv
open scoped LinearAlgebra.Projectivization

/-- Number of fixed points. -/
def nfix {n : ℕ} (h : Perm (Fin n)) : ℕ := (Finset.univ.filter (fun x => h x = x)).card

/-- g is induced, through e, by a σ-semilinear bijection f of F^d (an element of ΓL_d(F)):
f maps the line e x into the line e (g x). -/
def IsCollineationVia {n d : ℕ} {F : Type} [Field F] (e : Fin n ≃ ℙ F (Fin d → F))
    (g : Perm (Fin n)) : Prop :=
  ∃ (σ : F ≃+* F) (f : (Fin d → F) ≃+ (Fin d → F)),
    (∀ (c : F) (v : Fin d → F), f (c • v) = σ c • f v) ∧
    ∀ (x : Fin n) (v : Fin d → F), v ∈ (e x).submodule → f v ∈ (e (g x)).submodule

section

variable {F : Type} [Field F] {d n : ℕ}

theorem mem_point_iff (p : ℙ F (Fin d → F)) (v : Fin d → F) :
    v ∈ p.submodule ↔ ∃ a : F, a • p.rep = v := by
  rw [Projectivization.submodule_eq, Submodule.mem_span_singleton]

theorem rep_mem (p : ℙ F (Fin d → F)) : p.rep ∈ p.submodule :=
  (mem_point_iff p _).mpr ⟨1, one_smul _ _⟩

theorem eq_mk_of_mem {p : ℙ F (Fin d → F)} {v : Fin d → F} (hv : v ≠ 0) (h : v ∈ p.submodule) :
    p = Projectivization.mk F v hv := by
  obtain ⟨a, ha⟩ := (mem_point_iff p v).mp h
  have ha0 : a ≠ 0 := by rintro rfl; rw [zero_smul] at ha; exact hv ha.symm
  conv_lhs => rw [← Projectivization.mk_rep p]
  refine (Projectivization.mk_eq_mk_iff' F _ _ (Projectivization.rep_nonzero p) hv).mpr ⟨a⁻¹, ?_⟩
  rw [← ha, smul_smul, inv_mul_cancel₀ ha0, one_smul]

theorem point_eq_of_mem {p p' : ℙ F (Fin d → F)} {v : Fin d → F} (hv : v ≠ 0)
    (h : v ∈ p.submodule) (h' : v ∈ p'.submodule) : p = p' := by
  rw [eq_mk_of_mem hv h, eq_mk_of_mem hv h']

theorem mem_mk_self (v : Fin d → F) (hv : v ≠ 0) : v ∈ (Projectivization.mk F v hv).submodule := by
  rw [Projectivization.submodule_mk]; exact Submodule.mem_span_singleton_self v

theorem smul_mem_point {p : ℙ F (Fin d → F)} {v : Fin d → F} (c : F) (h : v ∈ p.submodule) :
    c • v ∈ p.submodule := p.submodule.smul_mem c h

variable {e : Fin n ≃ ℙ F (Fin d → F)} {g : Perm (Fin n)} {σ : F ≃+* F}
  {f : (Fin d → F) ≃+ (Fin d → F)}

theorem f_ne_zero {v : Fin d → F} (hv : v ≠ 0) : f v ≠ 0 := by
  intro h; apply hv; rw [← map_zero f] at h; exact f.injective h

/-- A fixed point's representative is an eigenvector: f v = λ v with λ ≠ 0. -/
theorem fixed_eigen (hf : ∀ (x : Fin n) (v : Fin d → F), v ∈ (e x).submodule → f v ∈ (e (g x)).submodule)
    {x : Fin n} (hx : g x = x) :
    ∃ lam : F, lam ≠ 0 ∧ f (e x).rep = lam • (e x).rep := by
  have h1 := hf x _ (rep_mem (e x))
  rw [hx] at h1
  obtain ⟨a, ha⟩ := (mem_point_iff (e x) _).mp h1
  refine ⟨a, ?_, ha.symm⟩
  rintro rfl
  rw [zero_smul] at ha
  exact f_ne_zero (f := f) (Projectivization.rep_nonzero (e x)) ha.symm

/-- If g is an involution and d ≥ 2, then f² = μ·id (μ ≠ 0) and σ² = id. -/
theorem f_sq_scalar (hd : 2 ≤ d) (hσ : ∀ (c : F) (v : Fin d → F), f (c • v) = σ c • f v)
    (hf : ∀ (x : Fin n) (v : Fin d → F), v ∈ (e x).submodule → f v ∈ (e (g x)).submodule)
    (hg2 : ∀ x, g (g x) = x) :
    ∃ μ : F, μ ≠ 0 ∧ (∀ v, f (f v) = μ • v) ∧ ∀ c, σ (σ c) = c := by
  have hsq : ∀ v : Fin d → F, v ≠ 0 → ∃ a : F, f (f v) = a • v := by
    intro v hv
    obtain ⟨x, hx⟩ : ∃ x, e x = Projectivization.mk F v hv := ⟨e.symm _, e.apply_symm_apply _⟩
    have h1 : v ∈ (e x).submodule := by rw [hx]; exact mem_mk_self v hv
    have h2 := hf (g x) _ (hf x v h1)
    rw [hg2, hx, Projectivization.submodule_mk, Submodule.mem_span_singleton] at h2
    obtain ⟨a, ha⟩ := h2
    exact ⟨a, ha.symm⟩
  choose a ha using hsq
  have hindep : ∀ v w (hv : v ≠ 0) (hw : w ≠ 0), (∀ c : F, v ≠ c • w) → a v hv = a w hw := by
    intro v w hv hw hvw
    have hvw0 : v + w ≠ 0 := by
      intro h; apply hvw (-1); rw [neg_one_smul]; exact eq_neg_of_add_eq_zero_left h
    have e1 := ha (v + w) hvw0
    rw [map_add, map_add, ha v hv, ha w hw, smul_add] at e1
    have e2 : (a (v + w) hvw0 - a v hv) • v = (a w hw - a (v + w) hvw0) • w := by
      rw [sub_smul, sub_smul, sub_eq_sub_iff_add_eq_add, ← e1]; abel
    by_cases hc : a (v + w) hvw0 - a v hv = 0
    · rw [hc, zero_smul] at e2
      rcases smul_eq_zero.mp e2.symm with h | h
      · rw [sub_eq_zero] at hc h; rw [← hc, h]
      · exact absurd h hw
    · exfalso
      apply hvw ((a (v + w) hvw0 - a v hv)⁻¹ * (a w hw - a (v + w) hvw0))
      rw [mul_smul, ← e2, smul_smul, inv_mul_cancel₀ hc, one_smul]
  let i0 : Fin d := ⟨0, by omega⟩
  let i1 : Fin d := ⟨1, by omega⟩
  have hi01 : i0 ≠ i1 := by simp [i0, i1, Fin.ext_iff]
  let u0 : Fin d → F := Pi.single i0 1
  let u1 : Fin d → F := Pi.single i1 1
  have hu0 : u0 ≠ 0 := by
    intro h; have := congrFun h i0; simp [u0] at this
  have hu1 : u1 ≠ 0 := by
    intro h; have := congrFun h i1; simp [u1] at this
  have hu01 : ∀ c : F, u0 ≠ c • u1 := by
    intro c h; have := congrFun h i0
    simp [u0, u1, Pi.single_apply, hi01] at this
  have hall : ∀ v (hv : v ≠ 0), a v hv = a u0 hu0 := by
    intro v hv
    by_cases hm : ∃ c : F, v = c • u0
    · obtain ⟨c, rfl⟩ := hm
      have hv1 : ∀ c' : F, c • u0 ≠ c' • u1 := by
        intro c' h
        have h0 := congrFun h i0
        simp [u0, u1, Pi.single_apply, hi01] at h0
        rw [h0, zero_smul] at hv; exact hv rfl
      rw [hindep _ _ hv hu1 hv1, hindep _ _ hu0 hu1 hu01]
    · push Not at hm
      exact hindep _ _ hv hu0 hm
  refine ⟨a u0 hu0, ?_, ?_, ?_⟩
  · intro h0
    have := ha u0 hu0
    rw [h0, zero_smul] at this
    exact f_ne_zero (f := f) (f_ne_zero (f := f) hu0) this
  · intro v
    by_cases hv : v = 0
    · rw [hv, map_zero, map_zero, smul_zero]
    · rw [ha v hv, hall v hv]
  · intro c
    by_cases hc : c = 0
    · rw [hc, map_zero, map_zero]
    · have hcu : c • u0 ≠ 0 := smul_ne_zero hc hu0
      have e1 := ha (c • u0) hcu
      rw [hσ, hσ, ha u0 hu0, hall _ hcu, smul_smul, smul_smul] at e1
      have hu0' := hu0
      have hμ : a u0 hu0 ≠ 0 := by
        intro h0
        have := ha u0 hu0
        rw [h0, zero_smul] at this
        exact f_ne_zero (f := f) (f_ne_zero (f := f) hu0) this
      have e2 : (σ (σ c) * a u0 hu0 - a u0 hu0 * c) • u0 = 0 := by rw [sub_smul, e1, sub_self]
      rcases smul_eq_zero.mp e2 with h | h
      · have : (σ (σ c) - c) * a u0 hu0 = 0 := by rw [sub_mul, mul_comm c]; exact h
        rcases mul_eq_zero.mp this with h' | h'
        · exact sub_eq_zero.mp h'
        · exact absurd h' hμ
      · exact absurd h hu0

end

end Riemann.Deep.LeanPrefix

