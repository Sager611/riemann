/-!
# Round 450 — node N8 (small gaps), round 19, falsification of a fresh theorem of ours (dpcuncond's Theorem U: Durkan–Pearce-Crump's two-line moments without RH give D(λ) > 0 for every λ > 2/3): the quoted inputs hold in the sources (Jutila's and Conrey's density bounds uniform for σ ≥ 1/2, Lamzouri l.221–230; Fujii's moments for all zeros, l.240–284; DPC uses RH only at l.57, 59–62, 140, 607, 610, 629 and its construction needs only θ < 1); the hairpin weight at 10,000 zero-sampled points equals the random baseline; three referees found only repairable gaps, all closed on re-review (CORRECT, SIGNIFICANT): DISCOVERY small-gap-two-thirds, record 0.768 → 2/3; Banaszczyk balancing, Cartan's lemma and the Combinatorial Nullstellensatz transfer no tool

* Node: N8 (19th round; R448 before). Round type: falsification of a fresh theorem of ours (Riemann/Deep/DPCUncond/PROOF.md,
  agent dpcuncond: unconditionally, for every λ > 2/3 a positive proportion of consecutive ordinates, with multiplicity, have
  γ_{n+1} − γ_n < 2πλ/log T), with the falsification aimed at the off-line hairpin bound and the density input.
* [check] Density input: Lamzouri 2611.14486 TeX l.221–230 states Conrey's N(σ, T) ≪ T^{1−(8/7−ε)(σ−1/2)} log T uniformly for
  σ ≥ 1/2 and cites Jutila's exponent 1 − ε (bibliography l.530–532); Baluyot's thesis confirms uniformity; Step 6a/6f use it
  as |S_a| ≤ C·TL·e^{−κa} for the depth classes d_ρL ∈ [a, a + 1): HOLDS.
* [check] Fujii and the double sum: TeX l.240–284: Fujii's moments count all zeros, h ≤ 100/log T, even moments suffice by
  Cauchy–Schwarz; Lemma DoubleSumZeros' proof gives Step 6c's Σ_γ Ñ_ρ^m ≪ TL(a + 2)^m: HOLDS.
* [check] RH in DPC (2610.28527 TeX): logical RH uses at l.57, 59–62, 140, 607, 610, 629 only; none in Prop 1 (l.262–388) or
  Lemma 6 (l.686–741) after the shift; the construction l.743–771 needs only 0 < φθ < 1 < φ(1 + θ) < 2 and θ < 1: HOLDS.
* [falsify] Hairpin weight at zero-sampled points: Möbius polynomial of length 500, the last 10,000 of 22,491 zeros against
  10,000 random points: mean |R R*(x − iy)| ratios 1.310, 1.307, 1.287 at y = 0.05, 0.1, 0.2 against 1.322, 1.311, 1.286:
  no correlation with the zeros (the rigorous bound needs none: Hölder with sparse kernels): HOLDS.
* [known] A scout's 'FLAW' (Steps 5c and 6b–6d are sketches, R448's obstruction) is a judgment of detail, not an error: the
  orchestrator's referees then checked those steps (below). No prior art (arXiv Oct–Nov 2026 listings, API, Zenodo, HAL).
* Referees (orchestrator, three independent, exact text pasted): Steps 0–3 and 7 GAP repairable (b₀ is a polynomial;
  oppositely signed coincident slits when γ₂ − γ₁ = h need the total signed jump; side endpoints off the shifted ordinates;
  the block endpoint in the reduction); Step 5 CORRECT; Steps 4 and 6 GAP repairable (box sides through on-line zeros,
  truncation of v, depth bins over all off-line zeros, the Ψ moments by cell Hölder–Young–Minkowski, the order of choices).
  dpcuncond repaired all eight; two re-reviews: every gap CLOSED, CORRECT at the level of a refereed proof sketch;
  all referees: SIGNIFICANT. Orchestrator reruns: threshold.py byte-identical, toy_weight.py matches.
* Gate: new theorem → DISCOVERIES (N8, small-gap-two-thirds) + RECORDS (0.768 → 2/3) + TREE N8 (the R448 'heuristic, open'
  clause replaced: DPC needs RH only beyond length T^{1/2}; below 2/3 it needs an unconditional (2+δ)-th moment of rough
  Möbius polynomials longer than T^{1/2}) + a mandatory Lean agent.
* Exploration (tool transfers, tested): Banaszczyk vector balancing (signed phase vectors only after finite sampling),
  Cartan's exceptional-disc lemma (only the triangle-inequality zero-pole bound), Alon's Combinatorial Nullstellensatz (one
  nonzero grid value, not real-rootedness): no transfer.
-/

namespace Fayed.R450
end Fayed.R450
