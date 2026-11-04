/-!
# Round 210 — a Zenodo proof of RH breaks at (7); a scout's book-Ramsey record fails; hunts clean

* Refutation: Cadeddu, "Proof of the Riemann Hypothesis through the limit of an approximated ξ function"
  (Zenodo 10.5281/zenodo.23960915, 4 Nov 2026) breaks at (7): convergence of the approximant to ξ(1−s) gives absolute
  error → 0 only; at a zero the denominator ξ(s(Q)) → 0, so the relative error need not vanish (the paper's own (6) says
  it tends to ∞). No Euler product is used, so the argument would also exclude the off-line zeros of the
  Davenport–Heilbronn function, which exist. Logged in REFUTATIONS.md.
* Not a record: a scout's circulant on Z/254 (55 differences, degree 110) for R(B₆₃, B₆₄) ≥ 255 (2610.36522) has no
  B₆₃, but its complement has books with 93 pages (orchestrator check, .scratch/r210/book.py).
* No counterexample: 2611.03057 (x² + 4y² + 5z² exceptions, window 3.0–3.2·10⁶), 2611.05178 (max-cut, triangle-free
  subcubic, n ≤ 9), 2611.02944 (PF scalar criterion, m ≤ 10⁶), 2611.04530 (Table 2 radius, 100 digits), 2611.05363
  (Krzyż bound, random samples), OEIS A401070 (k ≤ 10⁶). Not testable in two minutes: 2611.02063, 2611.00251.
-/

namespace Fayed.R210
end Fayed.R210
