/-!
# Round 472 — node N10 (explicit ψ, Mertens and Liouville bounds), round 15, the unseen Zenodo N10 records (no date filter, 140 records, about 18 unseen on topic): Papadopoulos's claimed RH proof treats Möbius values as independent coin flips (REFUTATIONS); Lakhal's RH-conditional remainder bound, Mittermeier's Part 3 reserve reduction (critical state rerun identical), Pavesi's Björner-complex decomposition and Bhattacharjee–Bhattacharya's Lean companion hold; N10 records stand; arithmetic QUE, Drinfeld/Goss zeta and Rogers–Ramanujan transfer no tool

* Node: N10 (15th round; R443 before). Round type: [check]/[falsify] of the unseen Zenodo N10 records, no date filter (lessons 216, 227):
  Zenodo on "Mertens function", "Mertens conjecture", "Chebyshev function", "prime counting function", "Skewes", "prime number
  theorem" explicit, "Ramanujan prime counting inequality", "Rosser" "Schoenfeld", "logarithmic integral", "Mobius function" sum,
  "Liouville function" sum, "Polya conjecture": 140 records, about 18 unseen on topic; 8 texts pre-extracted (lesson 228).
* [check] S. M. Lakhal, Zenodo 23892178 v1.1 (2 Nov 2026), remainder R(N) = ϑ(p_N) − ϑ(x̄) + (π(x̄) − N) ln x̄, x̄ = li⁻¹(N):
  Theorem 2 (R ≥ 0) and Theorem 3 (under RH, R(N) ≤ (1 + η) × printed bound for N ≥ 400, η ≤ 0.09) hold on 1,120 N in [400, 15000]
  (min R = 0.0580436 at N = 13238, max R/bound = 0.33587 at N = 467; direct vs integral form agree to 7·10⁻¹¹): HOLDS.
* [check] R. A. Mittermeier, Zenodo 23377241 v4 (18 Oct 2026), Part 3: Theorem 3.4 (RH ⟺ V_q ≥ 0 at every active event of
  Suzuki's Ψ) uses only the convexity between prime powers and the event structure (scout's symbolic audit); the orchestrator
  reran the critical state at q = 34,186,367 (dense sieve, 256-bit Arb, 18 s, 22 MB): prefix and state identical to the stored
  results, V_q ∈ [0.0214985808236810, 0.0214985808236811]; the 10¹⁰ anchor is Part 2's certificate, imported by hash; the
  infinite tail stays open (no RH claim): HOLDS.
* [check] L. E. Pavesi, Zenodo 22831623 (18 Sep 2026, sixth revision), Björner's complex of squarefree integers: Prop. 4.2
  (index = M(n)) is rank–nullity; Theorem 4.7 (maximal squarefree m ≤ n via primorials) holds for every n ≤ 700 (orchestrator,
  exact) and restates Definition 4.3; the 'critical exponent 1/2' of the kernel decomposition is numerical (§5, Remark 6.2); the
  obstructions concern his own spectral modification; no bound on M(x): HOLDS.
* [falsify] P. C. Papadopoulos, Zenodo 19394456 (3 Apr 2026): the claimed RH proof (l.3744–3817) models the Möbius values of
  squarefree numbers as independent fair 'virtual coin' trials, so its Borel–Cantelli bound concerns a random model, not M(N):
  FLAW → REFUTATIONS; its Table 2 counts are exact (L = 1000: 303 / 304).
* [check] Bhattacharjee–Bhattacharya, Zenodo 24200541 (13 Nov 2026), Lean 4 + Julia companion of 24150355 (R443): 24 finite
  theorems (energy identities, Gram rank obstruction, numeric constants), no sorry, admit, axiom or native_decide; the energy
  exponent 2Θ and any zero location are not formalized; |L(x)| ≤ x/2641 + 8.3√x checked to 3·10⁷, as in COROLLARIES: HOLDS.
* Not read beyond titles: Axler 10474452 (2016 π(x) bounds, superseded by his published ones), Kim–Nguyen data 11550925 (known),
  Caraccioli Abrego 16899990 (framework), Mertens 'falsifications' and frameworks 21575722, 23587181. N10 records stand.
* Exploration: arithmetic QUE (Holowinsky's shifted-convolution sieve bounds Hecke coefficients, no identification with μ or Λ),
  Drinfeld modules/Goss zeta (Newton-polygon distinct slopes fail for y² = x³ + 1 over F₅, P(u) = 1 + 5u²; R341 recorded the
  archimedean gap), Rogers–Ramanujan/Bailey chains (no map from q-series positivity to the Weil form): no transfer.
* Scouts: 8, 7 answered (pavesi timed out; the orchestrator ran its Theorem 4.7 check).
-/

namespace Fayed.R472
end Fayed.R472
