import Riemann.Deep.LeanPrefix.Small

/-!
# Kernel-checked G = S_n for the primitive odd-line triples with n = 23

(23, 18, 4), (23, 16, 6), (23, 12, 10): a word in r_n, r_m, r_k (words from
../PrefixRev/Lean/audit/word_23_*.txt) equals the transposition (0, P 0) on every point (kernel
evaluation, decide +kernel); then G = S_n (IsCycleP.G_eq_top_of_swap). These are the cases n = 23
of Jones's list (c); n = 11 is G_11_6_4 (Small.lean).
-/

namespace Riemann.Deep.LeanPrefix

open Equiv

theorem check_23_12_10 : (List.range 23).all
    (fun x => evalCode 23 12 10 116 2628610569547227077616883639360546291977164134632809607 x == swapN 0 11 x) = true := by decide +kernel

/-- (23, 12, 10): r_23, r_12, r_10 generate S_23. -/
theorem G_23_12_10 : G 23 12 10 = ⊤ := by
  have hP : P 23 12 10 ⟨0, by omega⟩ = ⟨11, by omega⟩ :=
    Fin.ext (by rw [P_val (by omega) (by omega)]; rfl)
  refine G_eq_top_of_swap (by omega) (by omega) (by decide) ⟨0, by omega⟩ ?_
  rw [hP, ← codePerm_eq_swap (by omega) (by omega) (by omega) (by omega) check_23_12_10]
  exact codePerm_mem _ _ _ _ _

theorem check_23_16_6 : (List.range 23).all
    (fun x => evalCode 23 16 6 333 152104814959565062523160982812154116845811772846210820600973877300531702500395065257343064292196751311040108108537702311537304431713405333454853079267869599806 x == swapN 0 7 x) = true := by decide +kernel

/-- (23, 16, 6): r_23, r_16, r_6 generate S_23. -/
theorem G_23_16_6 : G 23 16 6 = ⊤ := by
  have hP : P 23 16 6 ⟨0, by omega⟩ = ⟨7, by omega⟩ :=
    Fin.ext (by rw [P_val (by omega) (by omega)]; rfl)
  refine G_eq_top_of_swap (by omega) (by omega) (by decide) ⟨0, by omega⟩ ?_
  rw [hP, ← codePerm_eq_swap (by omega) (by omega) (by omega) (by omega) check_23_16_6]
  exact codePerm_mem _ _ _ _ _

theorem check_23_18_4 : (List.range 23).all
    (fun x => evalCode 23 18 4 1078 4322447099324449200257585316667407078114590972758878334529595193417877237618188382424330558058464980963477391770033723973103620319830241199658709049768792830696715290268669924207572302039343754142722667291768305536589535579310056068075360650852702266209197472945012555713070882678184092488420113966997872679677459701757107543737344096501831667027334804563144486041382180612531784173364556393889357235212929245241807580762645725315460961464432461437964099090469185168446939195843423005147682008004082969838899370485 x == swapN 0 5 x) = true := by decide +kernel

/-- (23, 18, 4): r_23, r_18, r_4 generate S_23. -/
theorem G_23_18_4 : G 23 18 4 = ⊤ := by
  have hP : P 23 18 4 ⟨0, by omega⟩ = ⟨5, by omega⟩ :=
    Fin.ext (by rw [P_val (by omega) (by omega)]; rfl)
  refine G_eq_top_of_swap (by omega) (by omega) (by decide) ⟨0, by omega⟩ ?_
  rw [hP, ← codePerm_eq_swap (by omega) (by omega) (by omega) (by omega) check_23_18_4]
  exact codePerm_mem _ _ _ _ _

end Riemann.Deep.LeanPrefix
