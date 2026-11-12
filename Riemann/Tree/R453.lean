/-!
# Round 453 — node N8 (critical-line proportion), round 21, can the arena's unconditional zeta23/TEN class reach Chirre–Gonçalves–de Laat's 0.6792? CGL use only Montgomery's |α| ≤ 1 theorem (RH, TeX l.474–478) and F(α) ≥ 0 (unconditional, l.446) with a Cohn–Elkies-type class (f even, f̂ ≥ 0, f eventually ≤ 0; SDP over p(x)e^{−πx²}); Carneiro–Chandee–Littmann–Milinovich (1406.5462) prove Montgomery–Taylor's 1.3274992 = 2^{−1/2}cot 2^{−1/2} + 1/2 optimal for the class R̂ ⊂ [−1, 1], which is the arena's κ₀ = 2 − 1.3274992; zeta23 has no F(α) theorem: off-line zeros enter through complex γ and are bounded by rank inequalities; the arena's first gains came from simple-zero refinements with window-count pressure and pyramid certificates; nobody has imported CGL or CCLM into an unconditional frame: Orca candidate; pluripotential theory, Jacobi forms and Marstrand projections transfer no tool

* Node: N8 (21st round; R452 before). Round type: feasibility check of an Orca candidate (import Chirre–Gonçalves–de Laat's
  semidefinite pair-correlation class into zeta23's unconditional frame), lesson 87: check before spending a slot.
* [check] CGL 1810.08843: inputs are Montgomery's inequality for |α| ≤ 1 (l.241–247; RH through the F-asymptotic l.474–478),
  F(α, T) ≥ 0 (l.446, unconditional); no lower bound for |α| > 1 under RH (GRH adds F ≥ 3/2 − |α| − o(1) on [1, 3/2 − δ],
  l.494–498); class A_LP: even continuous f ∈ L¹, f(0) = f̂(0) = 1, f̂ ≥ 0, f eventually ≤ 0 (l.420–430), f = p(x)e^{−πx²},
  SDP/SOS constraints (l.670–715): the analytic machinery is unconditional; only the F-asymptotic needs RH.
* [known] CCLM 1406.5462 Cor. cor_Mont_Taylor (TeX l.1243–1273): the exact one-delta optimum M(R) ≥ 2^{−1/2}cot(2^{−1/2}) − 1/2
  = 0.3274992…, N*(T) ≤ (1.3274992… + o(1))N(T), extremiser characterised (reproducing-kernel/de Branges optimality, l.284–289,
  496, 1317) for R̂ ⊂ [−1, 1]: Montgomery–Taylor's value is the ceiling of that class and the arena's κ₀ = 2 − 1/c_MT.
* [check] zeta23 port: MV (Hypotheses.lean l.110–114) is the Montgomery–Vaughan Hilbert inequality with gaps δ_r; EF (l.51–56)
  is Weil's pairing W f g = ∫ f̂ ḡ̂ dν_X; off-line zeros enter through gammaOf ρ = (ρ − 1/2)/i and complex u_ρu_ρᵀ, made
  Hermitian by reflection (Tail.lean l.487–515), and are bounded only by rank and positive-index inequalities (Defs.lean
  l.325–333, Assembly/Inputs.lean l.38–55); no F(α) asymptotic appears: a CGL import must be built from EF and MV.
* [check] Arena's first gains over κ₀ (1–2 Oct 2026): attempt-002 (0.6728525), five-point-pyramid (0.67288785), attempt-005,
  five-point-pyramid-2, weighted-five-point-pyramid (0.67298538): simple-zero refinements n = 4 → 6 with window-count
  pressure, then dyadic pyramid certificates and capacity-saturating pair weights; no stated limit.
* [known] No unconditional import of CGL or CCLM exists (arXiv API, Zenodo, our files); unconditional pair-correlation
  proportions: zeta23 2/3 (Lean), arena 0.673639288 (Lean), Schäfer 22837912 0.673425 (superseded); 2610.08965 (mollifier):
  0.6735382 simple on-line, 0.8390046 distinct, 0.678 simple.
* Exploration: Siciak / Bernstein–Walsh (discrete zero sets are polar; no Bernstein–Markov measure), Jacobi heat operator
  (reproduces known low-degree Jensen positivity), Marstrand–Kaufman projections (atomic zero measures have infinite
  energy): no transfer.
-/

namespace Fayed.R453
end Fayed.R453
