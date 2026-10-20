/-!
# Round 96 — an audit checklist for claimed RH proofs (lessons from ~15 audits)

| test | how to run (≤ 30 min) | catches | example |
|------|------------------------|---------|---------|
| Davenport–Heilbronn | instantiate the paper's hypotheses with DH (certified off-line zero, Fayed.DH.Zero.dh_offline_holds) | FE / symmetry / contour arguments without Euler-product input | R77 contour topology |
| σ-blindness | rerun the decisive step with σ₀ = 0.4 and 0.6 | localisation that cannot see Re s₀ | Fayed (R1, R28) |
| finite base + tail | separate the finite check from the tail; demand a uniform tail bound | numerics + extrapolation | Desogus, Gershon |
| false general lemma | extract every "well-known" lemma; try a counterexample | overgeneralisation | Fayed Thm 1, Wünsche, Ryazanov |
| positivity ≡ RH | check whether the claimed positivity is itself RH-equivalent | circularity | Yakaboylu |
| symmetry ≠ fixedness | does the paper prove ρ = 1 − conj ρ, or only pair zeros? | pairing mistaken for location | Micheas, Covei, Stenger? |
| numerics as proof | list heights/precision; find the step beyond the last height | finite → infinite leaps | Q(A) Zenodo |
| Lean fidelity (ours) | probe x/0, non-summable tsum, binders | vacuous formal statements | Lagarias, HL, Weil archTerm |
-/
