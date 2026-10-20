import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Nonvanishing
import Riemann.Tree.R007

/-!
# Bridge: OpenAI's quasi-Riemann hypothesis in this repository's vocabulary

`Fayed.R7.QuasiRH δ` (Riemann/Tree/R007.lean) says that ζ has no zero with `1 - δ < Re s`
(`s ≠ 1`); so `QuasiRH (1/2)` is RH, and the zero-free half-plane `Re s > 7/8` is
`QuasiRH (1/8)`. (`QuasiRH (7/8)` would mean no zero right of `1/8`: false, since the critical
line carries zeros.) The proof is OpenAI's `OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re`,
ported under Riemann/Deep/QRH/OAI/ (see NOTICE, REPORT.md).
-/

namespace Fayed.QRH

/-- The zero-free half-plane `Re s > 7/8`, i.e. `QuasiRH (1/8)` in R007's parametrization. -/
theorem quasiRH_seven_eighths : Fayed.R7.QuasiRH (1 / 8) := fun s hs _ =>
  OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re (by norm_num at hs ⊢; linarith)

/-- Corollary: any zero of ζ that would refute Fayed's claimed RH has `Re s ≤ 7/8`;
in particular ζ has no zero with `Re s > 7/8`. -/
theorem zeta_ne_zero_of_seven_eighths_lt_re (s : ℂ) (hs : 7 / 8 < s.re) : riemannZeta s ≠ 0 :=
  OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re (by exact_mod_cast hs)

end Fayed.QRH
