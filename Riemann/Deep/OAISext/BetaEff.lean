-- oaisext (Riemann/Deep/OAISext/REPORT.md): the zero supremum floored at the contour sentinel 51/100.
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Hecke.ZeroSupremum

namespace OAI

noncomputable section
namespace SevenEighths.HeckeZeroSupremum

/-- `betaEff = max β_* (51/100)`, where `beta = β_*` is the zero supremum of the Hecke family. -/
def betaEff : ℝ := max beta (51 / 100)

theorem fiftyOne_le_betaEff : (51 / 100 : ℝ) ≤ betaEff := le_max_right _ _

theorem beta_le_betaEff : beta ≤ betaEff := le_max_left _ _

theorem betaEff_le_one : betaEff ≤ 1 := max_le beta_le_one (by norm_num)

end SevenEighths.HeckeZeroSupremum

end

end OAI
