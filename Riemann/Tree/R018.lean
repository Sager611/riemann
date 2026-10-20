import Riemann.Tree.R007

/-!
# Round 18 — stepping stones strictly weaker than RH

| id | statement | best known | on path? | status |
|----|-----------|-----------|----------|--------|
| R18.1 | quasi-RH strip Re s > 1-δ | none (Ford 2025 region 1 - c(log t)^{-2/3}(log log t)^{-1/3}) | yes | OPEN stepping stone |
| R18.2 | Lindelöf | Bourgain 13/84 | yes (growth) | OPEN stepping stone |
| R18.3 | density hypothesis | Guth–Maynard 30/13 | partly | OPEN stepping stone |
| R18.4 | > 50% zeros on the line | 41.7% | no | DEAD as stepping stone |
| R18.5 | zero-free region exponent < 2/3 | 2/3 since 1958 | yes | OPEN stepping stone |
| R18.6 | > 19/27 simple zeros unconditionally | conditional only | no | DEAD |
| R18.7 | Λ_dBN < 0.1 | 0.2 | no | DEAD |
| R18.8 | full pair correlation | α ≤ 1 (Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh) | no | DEAD |
-/

namespace Fayed.R18

open Fayed.R7

/-- Quasi-RH strips are nested: a wider strip implies a narrower one. -/
theorem quasiRH_mono {δ δ' : ℝ} (h : δ' ≤ δ) (hq : QuasiRH δ) : QuasiRH δ' :=
  fun s hs hs1 => hq s (by linarith) hs1

/-- RH gives every strip up to width 1/2. -/
theorem quasiRH_of_RH (h : RiemannHypothesis) {δ : ℝ} (hδ : δ ≤ 1 / 2) : QuasiRH δ :=
  quasiRH_mono hδ (RH_iff_QuasiRH_half.1 h)

/-- R18.2: the Lindelöf hypothesis. -/
def Lindelof : Prop :=
  ∀ ε : ℝ, 0 < ε → ∃ C : ℝ, ∀ t : ℝ, 1 ≤ |t| →
    ‖riemannZeta ((1 / 2 : ℂ) + t * Complex.I)‖ ≤ C * |t| ^ ε

end Fayed.R18
