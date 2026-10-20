# Minimality of the first certified zero — REPORT

**Result: done, no sorry, no native_decide.** `lake build Riemann.Deep.FirstZero.Min` green.
`#print axioms` = [propext, Classical.choice, Quot.sound] for all three theorems (Riemann/Deep/FirstZero/Min.lean):
- `Fayed.FirstZero.no_zero_below`: ∀ s, ζ s = 0 → 0 < Re s < 1 → 0 < Im s < 14.13 → False.
- `Fayed.FirstZero.first_zero_is_first`: ∃ t₁ ∈ (14.13, 14.14), ζ(1/2+it₁) = 0 and t₁ ≤ Im s for every zero s in the strip with Im s > 0.
- `Fayed.FirstZero.zeros_below_15` (stronger): the only zero in the strip with 0 < Im s ≤ 15 is 1/2+it₁.

## Route (a), with ζ itself, no Γ needed
1. Argument principle for ζ on [1/4, 2] × [−15, 15], with the simple pole at 1 (ported zeta23 `rectangleIntegral'_mul_logDeriv_of_poles`), folded by ζ(s̄) = conj ζ(s): Σ ord − 1 = (1/π)·(arg change along 2 → 2+15i → 1/4+15i → 1/4). See `MinArg.winding`.
2. A zero ρ below 14.13 (or its mirror 1−ρ̄ if Re ρ < 1/4), its conjugate, the certified zero and its conjugate are 4 distinct zeros. So Σ ord ≥ 4, and the arg change must be ≥ 3π (`MinCount.unique_of_cert`).
3. Certified arg change ≤ 3.37 (true value π): 30 pieces, each Log(ζ·c̄) by the fundamental theorem of calculus with an anchor c that keeps Re(ζ·c̄) > 0 along the piece. The right edge is free because Re ζ > 0 there.
4. Re(ζ·c̄) > 0 on each piece comes from a 2nd-order Taylor model of the order-8 Euler–Maclaurin main term, with explicit derivatives, a crude bound on the second derivative, the mean value theorem, and the Euler–Maclaurin remainder (`MinModel.seg_model`). Node values use ball arithmetic as in Cert.lean, generalised to Re s ≠ 1/2 via rational powers n^{-p/q}.

## Files (Riemann/Deep/FirstZero/)
| file | content |
|---|---|
| MinModel.lean | gc derivatives, M/M1/M2 (Euler–Maclaurin main term + derivatives), crude bounds, taylor2, seg_model, arg-from-ball |
| MinArg.lean | logDeriv conj, FTC per piece, TopOK/LeftOK, right edge, `winding` |
| MinCount.lean | four-zero count ⇒ `unique_of_cert`, `no_zero_below_of_cert` |
| MinLib.lean | packaged numeric lemmas (nps_ball, zeta_node, arg_up/lo, leftOK_seg, topOK_seg) |
| MinCert0–5.lean | generated: 31 nodes, 30 pieces (.scratch/minimality/{genlib,plan,gen,drive}.py) |
| Min.lean | chains + final theorems |

## Notes
- Not imported from Riemann.lean (outside edit scope). Orchestrator: add `import Riemann.Deep.FirstZero.Min`.
- Slack: bound 3.37 vs limit 3π ≈ 9.42. Full build of the Min* files ≈ 6 min.
- Cheap follow-ups: the same count shows the first zero is simple (Σ ord < 4), and T can be pushed past 21 to cover the second zero.
