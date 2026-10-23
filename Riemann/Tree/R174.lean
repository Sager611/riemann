/-!
# Round 174 — the OEIS growth conjecture is plausible after all; new claimed counterexamples to audit

* OEIS A399309 and A399861: the sign sequence A191232 is the concatenation of the primes in base 2 (density of 1s
  0.5418 in the first 10⁶ bits; aperiodic, factor complexity 2^k for k ≤ 14). Lyapunov estimates (renormalized
  floating point, three seeds agree): growth 1.169697, 1.155324, 1.152191 at n = 10⁵, 10⁶, 10⁷. The excess over
  Viswanath's 1.1319882 shrinks like 1/log n (each prime contributes a deterministic leading and trailing 1, about
  2/log₂ p of the bits), so the data fit the conjecture. By Copeland–Erdős the sequence is normal in base 2, and for a
  generic point of the Bernoulli measure subadditivity gives limsup |a(n)|^(1/n) ≤ Viswanath's constant; the matching
  lower bound would need more. No counterexample.
* Wormald's conjecture, connected case: 2-edge-cut joins of cubic graphs on 6 and 10 vertices (18 base graphs) all
  split; no counterexample.
* Leads (claimed counterexamples and records, October titles): 2610.01536 (explicit polynomial counterexample to
  Connes' embedding conjecture), 2610.08581 (continuous Aldous–Lyons), 2610.13848 (Nirenberg, closed asymptotic
  curves), 2610.00153 (Yau, bounded mean curvature), 2610.13979 (three conjectures of Luo), 2610.31589 (thrackles
  ≤ 4n/3), 2610.02809 (Euclidean TSP constant), 2610.30027, 2610.08121, 2610.11406, 2610.13888, 2610.24541;
  Zenodo 23584043, 23583993, 23541166, 23497760, 23555563, 23580381.
-/

namespace Fayed.R174
end Fayed.R174
