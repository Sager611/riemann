import Riemann.Deep.Orders.Lean.Fail
import Riemann.Deep.Orders.Lean.Small1
import Riemann.Deep.Orders.Lean.Small2
import Riemann.Deep.Orders.Lean.Small3
import Riemann.Deep.Orders.Lean.Small4
import Riemann.Deep.Orders.Lean.Small5
import Riemann.Deep.Orders.Lean.Small6
import Riemann.Deep.Orders.Lean.D_pxmd
import Riemann.Deep.Orders.Lean.D_pxdm
import Riemann.Deep.Orders.Lean.D_pdxm
import Riemann.Deep.Orders.Lean.D_pdmx
import Riemann.Deep.Orders.Lean.D_xpmd
import Riemann.Deep.Orders.Lean.D_xpdm
import Riemann.Deep.Orders.Lean.D_xdpm
import Riemann.Deep.Orders.Lean.D_xdmp
import Riemann.Deep.Orders.Lean.D_dxpm
import Riemann.Deep.Orders.Lean.D_dpxm
import Riemann.Deep.Orders.Lean.D_dpmx

/-! Theorem 1 (exact sets of perfect N) and the Corollary (which orders work). -/
namespace Orders

/-- Theorem 1, + × − ÷: perfect exactly for N ∈ [] or N ≥ 9. -/
theorem thm1_pxmd (N : Nat) : (∃ ms, Perfect (⟨.add, .mul, .sub, .div⟩ : Order) N ms) ↔ N ∈ [] ∨ 9 ≤ N := by
  by_cases h : N < 9
  · rw [small_pxmd N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_pxmd.works N (by omega)⟩

/-- Theorem 1, + × ÷ −: perfect exactly for N ∈ [] or N ≥ 10. -/
theorem thm1_pxdm (N : Nat) : (∃ ms, Perfect (⟨.add, .mul, .div, .sub⟩ : Order) N ms) ↔ N ∈ [] ∨ 10 ≤ N := by
  by_cases h : N < 10
  · rw [small_pxdm N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_pxdm.works N (by omega)⟩

/-- Theorem 1, + − × ÷: perfect exactly for N ∈ []. -/
theorem thm1_pmxd (N : Nat) : (∃ ms, Perfect (⟨.add, .sub, .mul, .div⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_pmxd N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_pmxd N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, + − ÷ ×: perfect exactly for N ∈ []. -/
theorem thm1_pmdx (N : Nat) : (∃ ms, Perfect (⟨.add, .sub, .div, .mul⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_pmdx N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_pmdx N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, + ÷ × −: perfect exactly for N ∈ [] or N ≥ 9. -/
theorem thm1_pdxm (N : Nat) : (∃ ms, Perfect (⟨.add, .div, .mul, .sub⟩ : Order) N ms) ↔ N ∈ [] ∨ 9 ≤ N := by
  by_cases h : N < 9
  · rw [small_pdxm N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_pdxm.works N (by omega)⟩

/-- Theorem 1, + ÷ − ×: perfect exactly for N ∈ [4, 6, 10] or N ≥ 12. -/
theorem thm1_pdmx (N : Nat) : (∃ ms, Perfect (⟨.add, .div, .sub, .mul⟩ : Order) N ms) ↔ N ∈ [4, 6, 10] ∨ 12 ≤ N := by
  by_cases h : N < 12
  · rw [small_pdmx N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_pdmx.works N (by omega)⟩

/-- Theorem 1, × + − ÷: perfect exactly for N ∈ [] or N ≥ 9. -/
theorem thm1_xpmd (N : Nat) : (∃ ms, Perfect (⟨.mul, .add, .sub, .div⟩ : Order) N ms) ↔ N ∈ [] ∨ 9 ≤ N := by
  by_cases h : N < 9
  · rw [small_xpmd N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_xpmd.works N (by omega)⟩

/-- Theorem 1, × + ÷ −: perfect exactly for N ∈ [] or N ≥ 8. -/
theorem thm1_xpdm (N : Nat) : (∃ ms, Perfect (⟨.mul, .add, .div, .sub⟩ : Order) N ms) ↔ N ∈ [] ∨ 8 ≤ N := by
  by_cases h : N < 8
  · rw [small_xpdm N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_xpdm.works N (by omega)⟩

/-- Theorem 1, × − + ÷: perfect exactly for N ∈ []. -/
theorem thm1_xmpd (N : Nat) : (∃ ms, Perfect (⟨.mul, .sub, .add, .div⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_xmpd N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_xmpd N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, × − ÷ +: perfect exactly for N ∈ []. -/
theorem thm1_xmdp (N : Nat) : (∃ ms, Perfect (⟨.mul, .sub, .div, .add⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_xmdp N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_xmdp N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, × ÷ + −: perfect exactly for N ∈ [4, 5, 6, 8, 9, 10] or N ≥ 12. -/
theorem thm1_xdpm (N : Nat) : (∃ ms, Perfect (⟨.mul, .div, .add, .sub⟩ : Order) N ms) ↔ N ∈ [4, 5, 6, 8, 9, 10] ∨ 12 ≤ N := by
  by_cases h : N < 12
  · rw [small_xdpm N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_xdpm.works N (by omega)⟩

/-- Theorem 1, × ÷ − +: perfect exactly for N ∈ [4, 5, 6] or N ≥ 12. -/
theorem thm1_xdmp (N : Nat) : (∃ ms, Perfect (⟨.mul, .div, .sub, .add⟩ : Order) N ms) ↔ N ∈ [4, 5, 6] ∨ 12 ≤ N := by
  by_cases h : N < 12
  · rw [small_xdmp N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_xdmp.works N (by omega)⟩

/-- Theorem 1, − + × ÷: perfect exactly for N ∈ []. -/
theorem thm1_mpxd (N : Nat) : (∃ ms, Perfect (⟨.sub, .add, .mul, .div⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_mpxd N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_mpxd N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, − + ÷ ×: perfect exactly for N ∈ []. -/
theorem thm1_mpdx (N : Nat) : (∃ ms, Perfect (⟨.sub, .add, .div, .mul⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_mpdx N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_mpdx N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, − × + ÷: perfect exactly for N ∈ []. -/
theorem thm1_mxpd (N : Nat) : (∃ ms, Perfect (⟨.sub, .mul, .add, .div⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_mxpd N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_mxpd N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, − × ÷ +: perfect exactly for N ∈ []. -/
theorem thm1_mxdp (N : Nat) : (∃ ms, Perfect (⟨.sub, .mul, .div, .add⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_mxdp N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_mxdp N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, − ÷ + ×: perfect exactly for N ∈ []. -/
theorem thm1_mdpx (N : Nat) : (∃ ms, Perfect (⟨.sub, .div, .add, .mul⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_mdpx N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_mdpx N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, − ÷ × +: perfect exactly for N ∈ []. -/
theorem thm1_mdxp (N : Nat) : (∃ ms, Perfect (⟨.sub, .div, .mul, .add⟩ : Order) N ms) ↔ N ∈ [] := by
  by_cases h : N < 6
  · exact small_mdxp N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_mdxp N ms hp) (by omega), fun h' => ?_⟩
    simp at h'

/-- Theorem 1, ÷ + × −: perfect exactly for N ∈ [2, 3, 4] or N ≥ 9. -/
theorem thm1_dpxm (N : Nat) : (∃ ms, Perfect (⟨.div, .add, .mul, .sub⟩ : Order) N ms) ↔ N ∈ [2, 3, 4] ∨ 9 ≤ N := by
  by_cases h : N < 9
  · rw [small_dpxm N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_dpxm.works N (by omega)⟩

/-- Theorem 1, ÷ + − ×: perfect exactly for N ∈ [2, 3, 4] or N ≥ 12. -/
theorem thm1_dpmx (N : Nat) : (∃ ms, Perfect (⟨.div, .add, .sub, .mul⟩ : Order) N ms) ↔ N ∈ [2, 3, 4] ∨ 12 ≤ N := by
  by_cases h : N < 12
  · rw [small_dpmx N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_dpmx.works N (by omega)⟩

/-- Theorem 1, ÷ × + −: perfect exactly for N ∈ [2, 4, 5, 6, 8, 9, 10] or N ≥ 12. -/
theorem thm1_dxpm (N : Nat) : (∃ ms, Perfect (⟨.div, .mul, .add, .sub⟩ : Order) N ms) ↔ N ∈ [2, 4, 5, 6, 8, 9, 10] ∨ 12 ≤ N := by
  by_cases h : N < 12
  · rw [small_dxpm N h]; exact ⟨Or.inl, fun h' => h'.elim id (fun h' => absurd h' (by omega))⟩
  · exact ⟨fun _ => Or.inr (by omega), fun _ => D_dxpm.works N (by omega)⟩

/-- Theorem 1, ÷ × − +: perfect exactly for N ∈ [2, 4, 5]. -/
theorem thm1_dxmp (N : Nat) : (∃ ms, Perfect (⟨.div, .mul, .sub, .add⟩ : Order) N ms) ↔ N ∈ [2, 4, 5] := by
  by_cases h : N < 6
  · exact small_dxmp N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_dxmp N ms hp) (by omega), fun h' => ?_⟩
    simp only [List.mem_cons, List.not_mem_nil, or_false] at h'; omega

/-- Theorem 1, ÷ − + ×: perfect exactly for N ∈ [2]. -/
theorem thm1_dmpx (N : Nat) : (∃ ms, Perfect (⟨.div, .sub, .add, .mul⟩ : Order) N ms) ↔ N ∈ [2] := by
  by_cases h : N < 6
  · exact small_dmpx N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_dmpx N ms hp) (by omega), fun h' => ?_⟩
    simp only [List.mem_cons, List.not_mem_nil, or_false] at h'; omega

/-- Theorem 1, ÷ − × +: perfect exactly for N ∈ [2]. -/
theorem thm1_dmxp (N : Nat) : (∃ ms, Perfect (⟨.div, .sub, .mul, .add⟩ : Order) N ms) ↔ N ∈ [2] := by
  by_cases h : N < 6
  · exact small_dmxp N h
  · refine ⟨fun ⟨ms, hp⟩ => absurd (fail_dmxp N ms hp) (by omega), fun h' => ?_⟩
    simp only [List.mem_cons, List.not_mem_nil, or_false] at h'; omega

theorem cor_of_works {σ : Order} {R : Prop} (h : ∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms)
    (hr : R) : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms) ↔ R := ⟨fun _ => hr, fun _ => h⟩

theorem cor_of_fails {σ : Order} {R : Prop} (h : ∀ N ms, Perfect σ N ms → N ≤ 5) (hr : ¬R) :
    (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms) ↔ R := by
  constructor
  · rintro ⟨N0, hN⟩
    obtain ⟨ms, hms⟩ := hN (N0 + 6) (by omega)
    have := h _ _ hms; omega
  · intro r; exact absurd r hr

theorem cor_pxmd : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.add, .mul, .sub, .div⟩ : Order) N ms) ↔ (((⟨.add, .mul, .sub, .div⟩ : Order).s3 = .sub ∨ (⟨.add, .mul, .sub, .div⟩ : Order).s4 = .sub) ∧ (⟨.add, .mul, .sub, .div⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨9, D_pxmd.works⟩ (by decide)

theorem cor_pxdm : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.add, .mul, .div, .sub⟩ : Order) N ms) ↔ (((⟨.add, .mul, .div, .sub⟩ : Order).s3 = .sub ∨ (⟨.add, .mul, .div, .sub⟩ : Order).s4 = .sub) ∧ (⟨.add, .mul, .div, .sub⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨10, D_pxdm.works⟩ (by decide)

theorem cor_pmxd : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.add, .sub, .mul, .div⟩ : Order) N ms) ↔ (((⟨.add, .sub, .mul, .div⟩ : Order).s3 = .sub ∨ (⟨.add, .sub, .mul, .div⟩ : Order).s4 = .sub) ∧ (⟨.add, .sub, .mul, .div⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_pmxd (by decide)

theorem cor_pmdx : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.add, .sub, .div, .mul⟩ : Order) N ms) ↔ (((⟨.add, .sub, .div, .mul⟩ : Order).s3 = .sub ∨ (⟨.add, .sub, .div, .mul⟩ : Order).s4 = .sub) ∧ (⟨.add, .sub, .div, .mul⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_pmdx (by decide)

theorem cor_pdxm : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.add, .div, .mul, .sub⟩ : Order) N ms) ↔ (((⟨.add, .div, .mul, .sub⟩ : Order).s3 = .sub ∨ (⟨.add, .div, .mul, .sub⟩ : Order).s4 = .sub) ∧ (⟨.add, .div, .mul, .sub⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨9, D_pdxm.works⟩ (by decide)

theorem cor_pdmx : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.add, .div, .sub, .mul⟩ : Order) N ms) ↔ (((⟨.add, .div, .sub, .mul⟩ : Order).s3 = .sub ∨ (⟨.add, .div, .sub, .mul⟩ : Order).s4 = .sub) ∧ (⟨.add, .div, .sub, .mul⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨12, D_pdmx.works⟩ (by decide)

theorem cor_xpmd : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.mul, .add, .sub, .div⟩ : Order) N ms) ↔ (((⟨.mul, .add, .sub, .div⟩ : Order).s3 = .sub ∨ (⟨.mul, .add, .sub, .div⟩ : Order).s4 = .sub) ∧ (⟨.mul, .add, .sub, .div⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨9, D_xpmd.works⟩ (by decide)

theorem cor_xpdm : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.mul, .add, .div, .sub⟩ : Order) N ms) ↔ (((⟨.mul, .add, .div, .sub⟩ : Order).s3 = .sub ∨ (⟨.mul, .add, .div, .sub⟩ : Order).s4 = .sub) ∧ (⟨.mul, .add, .div, .sub⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨8, D_xpdm.works⟩ (by decide)

theorem cor_xmpd : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.mul, .sub, .add, .div⟩ : Order) N ms) ↔ (((⟨.mul, .sub, .add, .div⟩ : Order).s3 = .sub ∨ (⟨.mul, .sub, .add, .div⟩ : Order).s4 = .sub) ∧ (⟨.mul, .sub, .add, .div⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_xmpd (by decide)

theorem cor_xmdp : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.mul, .sub, .div, .add⟩ : Order) N ms) ↔ (((⟨.mul, .sub, .div, .add⟩ : Order).s3 = .sub ∨ (⟨.mul, .sub, .div, .add⟩ : Order).s4 = .sub) ∧ (⟨.mul, .sub, .div, .add⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_xmdp (by decide)

theorem cor_xdpm : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.mul, .div, .add, .sub⟩ : Order) N ms) ↔ (((⟨.mul, .div, .add, .sub⟩ : Order).s3 = .sub ∨ (⟨.mul, .div, .add, .sub⟩ : Order).s4 = .sub) ∧ (⟨.mul, .div, .add, .sub⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨12, D_xdpm.works⟩ (by decide)

theorem cor_xdmp : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.mul, .div, .sub, .add⟩ : Order) N ms) ↔ (((⟨.mul, .div, .sub, .add⟩ : Order).s3 = .sub ∨ (⟨.mul, .div, .sub, .add⟩ : Order).s4 = .sub) ∧ (⟨.mul, .div, .sub, .add⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨12, D_xdmp.works⟩ (by decide)

theorem cor_mpxd : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.sub, .add, .mul, .div⟩ : Order) N ms) ↔ (((⟨.sub, .add, .mul, .div⟩ : Order).s3 = .sub ∨ (⟨.sub, .add, .mul, .div⟩ : Order).s4 = .sub) ∧ (⟨.sub, .add, .mul, .div⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_mpxd (by decide)

theorem cor_mpdx : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.sub, .add, .div, .mul⟩ : Order) N ms) ↔ (((⟨.sub, .add, .div, .mul⟩ : Order).s3 = .sub ∨ (⟨.sub, .add, .div, .mul⟩ : Order).s4 = .sub) ∧ (⟨.sub, .add, .div, .mul⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_mpdx (by decide)

theorem cor_mxpd : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.sub, .mul, .add, .div⟩ : Order) N ms) ↔ (((⟨.sub, .mul, .add, .div⟩ : Order).s3 = .sub ∨ (⟨.sub, .mul, .add, .div⟩ : Order).s4 = .sub) ∧ (⟨.sub, .mul, .add, .div⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_mxpd (by decide)

theorem cor_mxdp : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.sub, .mul, .div, .add⟩ : Order) N ms) ↔ (((⟨.sub, .mul, .div, .add⟩ : Order).s3 = .sub ∨ (⟨.sub, .mul, .div, .add⟩ : Order).s4 = .sub) ∧ (⟨.sub, .mul, .div, .add⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_mxdp (by decide)

theorem cor_mdpx : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.sub, .div, .add, .mul⟩ : Order) N ms) ↔ (((⟨.sub, .div, .add, .mul⟩ : Order).s3 = .sub ∨ (⟨.sub, .div, .add, .mul⟩ : Order).s4 = .sub) ∧ (⟨.sub, .div, .add, .mul⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_mdpx (by decide)

theorem cor_mdxp : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.sub, .div, .mul, .add⟩ : Order) N ms) ↔ (((⟨.sub, .div, .mul, .add⟩ : Order).s3 = .sub ∨ (⟨.sub, .div, .mul, .add⟩ : Order).s4 = .sub) ∧ (⟨.sub, .div, .mul, .add⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_mdxp (by decide)

theorem cor_dpxm : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.div, .add, .mul, .sub⟩ : Order) N ms) ↔ (((⟨.div, .add, .mul, .sub⟩ : Order).s3 = .sub ∨ (⟨.div, .add, .mul, .sub⟩ : Order).s4 = .sub) ∧ (⟨.div, .add, .mul, .sub⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨9, D_dpxm.works⟩ (by decide)

theorem cor_dpmx : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.div, .add, .sub, .mul⟩ : Order) N ms) ↔ (((⟨.div, .add, .sub, .mul⟩ : Order).s3 = .sub ∨ (⟨.div, .add, .sub, .mul⟩ : Order).s4 = .sub) ∧ (⟨.div, .add, .sub, .mul⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨12, D_dpmx.works⟩ (by decide)

theorem cor_dxpm : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.div, .mul, .add, .sub⟩ : Order) N ms) ↔ (((⟨.div, .mul, .add, .sub⟩ : Order).s3 = .sub ∨ (⟨.div, .mul, .add, .sub⟩ : Order).s4 = .sub) ∧ (⟨.div, .mul, .add, .sub⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_works ⟨12, D_dxpm.works⟩ (by decide)

theorem cor_dxmp : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.div, .mul, .sub, .add⟩ : Order) N ms) ↔ (((⟨.div, .mul, .sub, .add⟩ : Order).s3 = .sub ∨ (⟨.div, .mul, .sub, .add⟩ : Order).s4 = .sub) ∧ (⟨.div, .mul, .sub, .add⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_dxmp (by decide)

theorem cor_dmpx : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.div, .sub, .add, .mul⟩ : Order) N ms) ↔ (((⟨.div, .sub, .add, .mul⟩ : Order).s3 = .sub ∨ (⟨.div, .sub, .add, .mul⟩ : Order).s4 = .sub) ∧ (⟨.div, .sub, .add, .mul⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_dmpx (by decide)

theorem cor_dmxp : (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect (⟨.div, .sub, .mul, .add⟩ : Order) N ms) ↔ (((⟨.div, .sub, .mul, .add⟩ : Order).s3 = .sub ∨ (⟨.div, .sub, .mul, .add⟩ : Order).s4 = .sub) ∧ (⟨.div, .sub, .mul, .add⟩ : Order) ≠ ⟨.div, .mul, .sub, .add⟩) :=
  cor_of_fails fail_dmxp (by decide)

/-- **Corollary.** An order σ of + × − ÷ has perfect sequences for all large N iff − is
third or fourth and σ ≠ ÷ × − + (11 of the 24 orders). -/
theorem corollary (σ : Order) (hσ : σ.IsPerm) :
    (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms) ↔
      ((σ.s3 = .sub ∨ σ.s4 = .sub) ∧ σ ≠ ⟨.div, .mul, .sub, .add⟩) := by
  obtain ⟨a, b, c, d⟩ := σ
  cases a <;> cases b <;> cases c <;> cases d <;>
  first
  | exact absurd hσ (by unfold Order.IsPerm; decide)
  | exact cor_pxmd
  | exact cor_pxdm
  | exact cor_pmxd
  | exact cor_pmdx
  | exact cor_pdxm
  | exact cor_pdmx
  | exact cor_xpmd
  | exact cor_xpdm
  | exact cor_xmpd
  | exact cor_xmdp
  | exact cor_xdpm
  | exact cor_xdmp
  | exact cor_mpxd
  | exact cor_mpdx
  | exact cor_mxpd
  | exact cor_mxdp
  | exact cor_mdpx
  | exact cor_mdxp
  | exact cor_dpxm
  | exact cor_dpmx
  | exact cor_dxpm
  | exact cor_dxmp
  | exact cor_dmpx
  | exact cor_dmxp

end Orders
