/-!
# Round 241 — OpenAI's combinatorics proofs audited (no flaw); two scout flaw claims rejected

* Audit (Orca "oaiproofs", Riemann/Deep/OAIProofs/REPORT.md): OpenAI's proofs of Seymour's second-neighbourhood
  conjecture, Saxl's conjecture and the tensor-square conjecture, plus its Hill and Zarankiewicz crossing-number proofs:
  every step read, every lemma with a finite instance retested in stronger forms (e.g. the pruning lemma on all 512
  relations on 3 points, 134M instances): no flaw found. Not our results.
* Rejected (orchestrator checks): 2609.15133 (resistance diameter of line graphs): K₄ − e has Dr(G) = 1 and
  Dr(L(G)) = 2/3 (L(G) is the wheel W₄), strict, as the theorem says; 2609.24239: "[F, 7]_q" has minimum distance 7
  (notation [F, k, δ]_q), so the Etzion–Silberstein bound 6 on the dimension is no contradiction.
* Verified (others' results, correct): 2609.00583 (21-colouring of the plane lattice, separation 5/2), 2610.28541
  (M(20, PG(2,7)) ≤ 21 instance), 2609.14526 (|P₁₀₀| = 2040), 2609.37774 (n = 11 sample).
* Not testable: 2609.37459 (labelled conjecture), Zenodo 23848597 (finite range of a positivity criterion).
-/

namespace Fayed.R241
end Fayed.R241
