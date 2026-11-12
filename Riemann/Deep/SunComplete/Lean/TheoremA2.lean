import Riemann.Deep.SunComplete.Lean.TheoremA
import Riemann.Deep.SunComplete.Lean.L6b

/-!
# sunlean · TheoremA2: Theorem A with Lemma 6 reduced to its middle range (`Lemma6Amid`)
-/

namespace SunLean

/-- **Theorem A** (Sun–Hou, Conjecture 3.5(i)) from Chvátal–Erdős, Lemma 3, Proposition 7's data,
the middle range of Lemma 6 (`30140 ≤ q < 53,876,069,761,024`, `ω(q − 1) ≤ 15`), the `5
certificates and `FFacts F`. -/
theorem theoremA_of_mid (hCE : ChvatalErdos) (hL3 : Lemma3A) (hP7 : Prop7A)
    (hL6 : Lemma6Amid) (hCert : CertA) (F : Type) [Field F] [Fintype F] (hF : FFacts F)
    (h7 : 7 < Fintype.card F) : SunHou F :=
  theoremA_of hCE hL3 hP7 (lemma6A_of_mid hL6) hCert F hF h7

/-- **Theorem A** with the sharper middle range of Lemma 6 (`Lemma6Amid2`:
`30140 ≤ q < 118,380,036,096`, `ω(q − 1) ≤ 11`, `φ(q − 1) < 7·2^{ω(q−1)}·(⌊√q⌋ + 1)`). -/
theorem theoremA_of_mid2 (hCE : ChvatalErdos) (hL3 : Lemma3A) (hP7 : Prop7A)
    (hL6 : Lemma6Amid2) (hCert : CertA) (F : Type) [Field F] [Fintype F] (hF : FFacts F)
    (h7 : 7 < Fintype.card F) : SunHou F :=
  theoremA_of hCE hL3 hP7 (lemma6A_of_mid2 hL6) hCert F hF h7

end SunLean

#print axioms SunLean.theoremA_of_mid
#print axioms SunLean.theoremA_of_mid2
