module

public import Riemann.Deep.Combos.Lean.SieveC
public import Riemann.Deep.Combos.Lean.ZetaB

@[expose] public section

/-!
# Combos NOTE §1, Lemma 2 (weighted large sieve for ring class characters), NOTE l.76–89

For 0 < ε ≤ 1/2, N ≥ 1 and complex c_𝔞 supported on a finite set S of ideals prime to f with
N𝔞 ≤ N, and every finite set 𝒞 of characters of G = Cl(O) (in particular 𝒞 = Ĝ):
Σ_{χ∈𝒞} |Σ_𝔞 c_𝔞 χ(𝔞)|² ≤ (30/ε) log²(3|D|) (N + |D|^{1/2}) Σ_𝔞 m(𝔞)^{1+2ε} |c_𝔞|²,
from (F1)+(F2) (F12) and (F3) only. |z|² is Complex.normSq z.
-/

open Finset

namespace RingClass

variable {G ι : Type*} [AddCommGroup G]

/-- (F3) (NOTE l.52): h(D) = |Cl(O)| ≤ 3|D|^{1/2} log²(3|D|). -/
def F3 [Fintype G] (d : RCData G ι) : Prop :=
  (Fintype.card G : ℝ) ≤ 3 * |(d.D : ℝ)| ^ (1 / 2 : ℝ) * Real.log (3 * |(d.D : ℝ)|) ^ 2

/-- ζ(1+2ε) ≤ 1 + 1/(2ε) ≤ 1/ε for 0 < ε ≤ 1/2 (NOTE l.89), on finite sets of positive integers -/
lemma zeta_eps (ε : ℝ) (hε : 0 < ε) (hε2 : ε ≤ 1 / 2) (K : Finset ℕ) (hK : ∀ k ∈ K, 0 < k) :
    ∑ k ∈ K, (((k : ℕ) : ℝ) ^ (1 + 2 * ε))⁻¹ ≤ 1 / ε := by
  have h := sum_rpow_neg_le (2 * ε) (by linarith) K hK
  have h3 : 1 ≤ 1 / (2 * ε) := by rw [le_div_iff₀ (by linarith)]; linarith
  have h4 : 1 / ε = 1 / (2 * ε) + 1 / (2 * ε) := by field_simp; ring
  linarith

/-- **Lemma 2** (NOTE l.76–89). -/
theorem lemma2 [Fintype G] [DecidableEq G] [DecidableEq ι] (d : RCData G ι) (hF12 : F12 d)
    (hF3 : F3 d) (ε : ℝ) (hε : 0 < ε) (hε2 : ε ≤ 1 / 2) (N : ℝ) (hN : 1 ≤ N) (S : Finset ι)
    (hS : ∀ A ∈ S, (d.nrm A : ℝ) ≤ N) (c : ι → ℂ) (𝒞 : Finset (AddChar G ℂ)) :
    ∑ χ ∈ 𝒞, Complex.normSq (∑ A ∈ S, c A * χ (d.cls A)) ≤
      30 / ε * Real.log (3 * |(d.D : ℝ)|) ^ 2 * (N + |(d.D : ℝ)| ^ (1 / 2 : ℝ)) *
        ∑ A ∈ S, (d.cont A : ℝ) ^ (1 + 2 * ε) * Complex.normSq (c A) := by
  have hsq : Real.sqrt |(d.D : ℝ)| = |(d.D : ℝ)| ^ (1 / 2 : ℝ) := Real.sqrt_eq_rpow _
  have hF3' : (Fintype.card G : ℝ) ≤
      3 * Real.sqrt |(d.D : ℝ)| * Real.log (3 * |(d.D : ℝ)|) ^ 2 := by
    rw [hsq]; exact hF3
  have h := sieve_core d hF12 (Real.log (3 * |(d.D : ℝ)|)) hF3' (fun m => (m : ℝ) ^ (1 + 2 * ε))
    (fun m hm => Real.rpow_pos_of_pos (by exact_mod_cast hm) _) (1 / ε)
    (fun K hK => zeta_eps ε hε hε2 K hK) 𝒞 N hN S hS c
  beta_reduce at h
  rw [hsq] at h
  calc _ ≤ _ := h
    _ = _ := by ring

end RingClass

