import Riemann.Deep.FasiHigham.Lean.Main
import Riemann.Deep.FasiHigham.Lean.PS

/-!
# Fasi–Higham, m = 3: every Taylor coefficient of 1/q_{k,3} is positive (power-series form)
-/

namespace FasiLean

/-- Fasi–Higham (SIMAX 40 (2019), Conj. 3.2) for m = 3 and every k: all coefficients of the formal
power-series inverse of q_{k,3} are positive. -/
theorem fasi_higham_m3 (k n : ℕ) : 0 < PowerSeries.coeff n (qPS k 3)⁻¹ := by
  rw [coeff_inv]
  exact c_pos k n

end FasiLean
