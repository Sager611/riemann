import Riemann.Deep.Orders.Lean.Patch

/-! # Theorem 2 scheme: base range by checked witnesses, then N = 3P + c by patches -/
namespace Orders

/-- Base cases: explicit sequences, accepted by the checker. -/
theorem base_of {σ : Order} {N0 K : Nat} (seqs : List (List Move))
    (h : (List.range' N0 K).all (fun N => check σ N (seqs.getD (N - N0) [])) = true) :
    ∀ N, N0 ≤ N → N < N0 + K → ∃ ms, Perfect σ N ms := by
  intro N h1 h2
  rw [List.all_eq_true] at h
  exact ⟨_, (check_iff _ _ _).1 (h N (List.mem_range'.2 ⟨N - N0, by omega, by omega⟩))⟩

/-- Strong induction: every N ≥ N1 is 3(4q + r) + c for an offset c with q ≥ q*_c and
P = 4q + r ≥ N0; the patch for c extends a perfect sequence for P. -/
theorem works_of {σ : Order} (hσ : σ.IsPerm) (r N0 N1 : Nat) (h0 : 1 ≤ N0) (pt : Nat → Patch)
    (offs : List Nat) (hok : ∀ c ∈ offs, patchOK σ r c (pt c) = true)
    (hbase : ∀ N, N0 ≤ N → N < N1 → ∃ ms, Perfect σ N ms)
    (hcov : ∀ N, N1 ≤ N → ∃ c ∈ offs, ∃ q, N = 3 * (4 * q + r) + c ∧ (pt c).qstar ≤ q ∧
      N0 ≤ 4 * q + r) :
    ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms := by
  intro N
  induction N using Nat.strongRecOn with
  | ind N ih =>
    intro hN
    by_cases h : N < N1
    · exact hbase N hN h
    · obtain ⟨c, hc, q, rfl, hq, hP⟩ := hcov N (by omega)
      obtain ⟨ms, hms⟩ := ih (4 * q + r) (by omega) hP
      exact ⟨_, patch_sound hσ (hok c hc) hq hms⟩

end Orders

