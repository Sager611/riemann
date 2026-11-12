import Riemann.Deep.OddRobin.Lean.Main

namespace OddRobin

theorem c0_eq : c0 = 67 * NN := by norm_num [c0, NN]

theorem candidates_eq : m3 = 3 * NN ∧ m11 = 11 * NN ∧ 5 * m67 = 67 * NN := by
  norm_num [m3, m11, m67, NN]

end OddRobin

#print axioms OddRobin.gap_theorem
#print axioms OddRobin.sharp
#print axioms OddRobin.main
#print axioms OddRobin.c0_eq
#print axioms OddRobin.candidates_eq

