import Riemann.Deep.SunComplete.Lean.CaseA1
import Riemann.Deep.SunComplete.Lean.Prop4

/-!
# sunlean · TheoremA: the logical architecture of NOTE §7 for Theorem A

`theoremA_of`: Chvátal–Erdős + Lemma 3 (named hypotheses) + Proposition 7's per-`q` data and the
§5 certificates (certificate hypotheses) + Lemma 6 (named hypothesis) + the finite-field facts
`FFacts` (proved from Mathlib in `FFInst.lean`) ⇒ every finite field with `q > 7` elements has the
Sun–Hou property.
-/

namespace SunLean

variable {F : Type} [Field F] [Fintype F] [DecidableEq F]

/-- The `L3Data` of `H_q` for given sieve data and `s ≥ √q`: `e = pairCount`, `Ñ = q − 1`,
`ε = epsA q`, Lemma 3 from `Lemma3A`. -/
noncomputable def l3dataA (hL3 : Lemma3A) (hF : FFacts F) (h2 : 2 ≤ Fintype.card F)
    (d : SieveData (Fintype.card F)) (s : ℚ) (hδ : 0 < d.δ) (hs : 0 ≤ s)
    (hsq : (Fintype.card F : ℚ) ≤ s ^ 2) : L3Data (sumGraph F) where
  e := pairCount
  τ := d.τ
  E := d.E s
  Nt := (Fintype.card F : ℚ) - 1
  ε := epsA (Fintype.card F)
  hτ := d.τ_pos hδ
  hE := d.E_nonneg hδ hs
  hε := epsA_nonneg _
  hdiag := fun _ hI => pairCount_diag h2 hF hI
  hcross := fun _ _ hAB hno => pairCount_cross hAB hno
  hL3 := fun X Y => hL3 F d s hδ hs hsq X Y

/-- Proposition 4 for `H_q`: `CritA` ⇒ Hamiltonian. -/
theorem ham_of_critA (hCE : ChvatalErdos) (hL3 : Lemma3A) (hF : FFacts F)
    (h3 : 3 ≤ Fintype.card F) {d : SieveData (Fintype.card F)} {s : ℚ} {C : ℕ}
    (hc : CritA (Fintype.card F) d s C) : HasHamCycle (sumGraph F) := by
  obtain ⟨hδ, hs, hsq, hNE, hC, ha1, hF'⟩ := hc
  exact (l3dataA hL3 hF (by omega) d s hδ hs hsq).prop4 hCE (fun v => deg_ge hF (by omega) v)
    h3 hNE C hC ha1 hF'

/-- Lemma 5 for `H_q`: `CSA` ⇒ Hamiltonian. -/
theorem ham_of_CSA (hCE : ChvatalErdos) (hL3 : Lemma3A) (hF : FFacts F)
    (h3 : 3 ≤ Fintype.card F) {d : SieveData (Fintype.card F)} {s : ℚ}
    (hc : CSA (Fintype.card F) d s) : HasHamCycle (sumGraph F) := by
  obtain ⟨hδ, hs, hsq, A', Dd, hA', hD, hDA, h2A, hCS⟩ := hc
  have h3' : (3 : ℚ) ≤ Fintype.card F := by exact_mod_cast h3
  have hNt : (0 : ℚ) ≤ (Fintype.card F : ℚ) - 1 := by linarith
  exact (l3dataA hL3 hF (by omega) d s hδ hs hsq).lemma5 hCE (fun v => deg_ge hF (by omega) v)
    h3 hNt hA' hD hDA h2A hCS

end SunLean

namespace SunLean

/-- **Theorem A** (Sun–Hou, Conjecture 3.5(i)), the architecture of NOTE §7: from Chvátal–Erdős,
Lemma 3, Proposition 7's data, Lemma 6, the §5 certificates and `FFacts F`, every finite field
`F` with `q > 7` elements has a numbering whose cyclic adjacent sums are all primitive. -/
theorem theoremA_of (hCE : ChvatalErdos) (hL3 : Lemma3A) (hP7 : Prop7A) (hL6 : Lemma6A)
    (hCert : CertA) (F : Type) [Field F] [Fintype F] (hF : FFacts F)
    (h7 : 7 < Fintype.card F) : SunHou F := by
  classical
  by_cases hE : Fintype.card F ∈ EA
  · exact hCert F hE h7
  · rw [sunHou_iff (by omega)]
    by_cases hlt : Fintype.card F < 30140
    · obtain ⟨d, s, C, hc⟩ := hP7 _ hF.primePow h7 hlt hE
      exact ham_of_critA hCE hL3 hF (by omega) hc
    · obtain ⟨d, s, hc⟩ := hL6 _ hF.primePow (by omega)
      exact ham_of_CSA hCE hL3 hF (by omega) hc

end SunLean
