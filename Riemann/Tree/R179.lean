/-!
# Round 179 — claimed counterexamples checked: two correct, one under-specified (multiple unicast), two not checkable

* Correct (others' results, our own computation): 2610.15145 (Game of Arrows on trees: the spider Sp(1,3,7) has 11
  markable edges, odd, yet Grundy value 0 by exhaustive search over 1,673 reachable states, refuting the parity
  conjecture of Lin, Zandee, Fokkink–Zandee); 2610.27294 (total-current Hohenberg–Kohn: the explicit two-fermion example
  with A and −A has equal density and total current, ∇ × A ≠ 0; mpmath residual < 1.5·10⁻⁴²).
* Under-specified: 2610.10108 (claimed counterexample to the multiple unicast conjecture of network coding): the routing
  side checks (182-vertex PG(2, 9) incidence network, 735 edges, all 157 session distances 5, so flow ≤ 147/157 < 1),
  but the encoding of F₉ and the point and line ordering are not fixed, and under a natural choice the stated property
  d⁻ = d⁺ = 5 fails. Orca audit next.
* Not checkable in 50 s: 2610.31637 (Simon's L² conjecture; infinite-scale construction), 2610.15199 (Rényi area law;
  the entropy algebra checks).
* Agrees on small cases: 2610.01240 (Lehmer's permutation conjecture; 58 signatures, n ≤ 8). Timed out: Pólya annuli
  2610.05146, MDS range 2610.03740.
-/

namespace Fayed.R179
end Fayed.R179
