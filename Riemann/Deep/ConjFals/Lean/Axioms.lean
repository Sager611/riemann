import Riemann.Deep.ConjFals.Lean.Main
import Riemann.Deep.ConjFals.Lean.TwoAdic

/-! Axiom check for the omegalean development (lim_{a→∞} s_a = 1). -/

#check @OmegaLean.tendsto_s
#check @OmegaLean.explicit_bounds
#check @OmegaLean.le_s
#check @OmegaLean.s_le
#check @OmegaLean.WN_two_adic
#print OmegaLean.W
#print OmegaLean.s
#print OmegaLean.U
#print OmegaLean.L
#print axioms OmegaLean.tendsto_s
#print axioms OmegaLean.explicit_bounds
#print axioms OmegaLean.le_s
#print axioms OmegaLean.s_le
#print axioms OmegaLean.WN_eq
#print axioms OmegaLean.WN_two_adic
