IDEA: Weil window certificates can be split by deflating an explicit near-null ladder (the window cut-offs of Φ, Φ″, Φ⁗, …, where Φ is the Riemann–Pólya kernel whose Fourier transform Ξ kills Weil's distribution unconditionally), which leaves a macroscopically positive, RH-blind remainder and confines all high precision, and every failure we provoked, to a block of size k ≈ c (node N1).

# explore17 — Gaussian free field, Liouville quantum gravity, Kahane's positive-type kernels → RH tree

Plain summary. The field itself (log-correlated fields, chaos, Polyakov determinants, DOZZ) gives nothing new for the tree: its one
positivity tool, Kahane/Junnila–Saksman–Webb "positive on small balls", reproduces at best the known analytic window c = 2. But the
Kahane question "is the window covariance a positive kernel?" led to a structural fact about Weil's form that is useful for N1: the
tiny eigenvalues that make every window certificate need ~17c bits are not mysterious; they come from an explicit family of functions,
and removing that family leaves a form that is comfortably positive. Parts of this picture are Connes's (below); the deflation tests,
the location of the failures and the certificate consequence are, as far as checked, new. No theorem, no record, no RH lever.

Conventions as in ../WeilSimple and ../WeilZeta: window W = [−L/2, L/2], L = ln c, primes n < c, w_n = Λ(n)/√n,
Q(f) = V(L)‖f‖² + ∫₀^L (C_f(0) − C_f(x)) k(x) dx + 2⟨f,cosh(x/2)⟩² − 2⟨f,sinh(x/2)⟩² − 2 Σ_{n<c} w_n C_f(log n),
k(x) = e^{x/2}/sinh x, V(L) = ψ(1/4) − log π + 2 Σ_j c^{−(2j+1/2)}/(2j + 1/2) (derived from φ(t) = Re ψ(1/4 + it/2) − log π by
Re ψ(1/4 + it/2) − ψ(1/4) = ∫₀^∞ (1 − cos tx) k(x) dx). Galerkin space E_N = CCM's cos/sin basis; our float64 matrices agree with
WeilSimple/weil_matrix.py (Arb) to all printed digits at c = 2, 3 (code/t_check.py).

## 1. Survey (one page)
- Kahane 1985 / Junnila–Saksman–Webb, "Decompositions of log-correlated fields with applications" (1808.06838, AAP 2019): a kernel
  is σ-positive if it is a sum of continuous kernels that are positive definite AND pointwise ≥ 0 (needed for Kahane's convexity
  inequality); JSW: log(1/|x − y|) + g(x, y), g ∈ H^{d+ε}, is locally the covariance of an almost ⋆-scale-invariant field plus a
  Hölder field (positivity on small balls). Transfer facts: (i) translation-invariant pieces positive definite on ℝ sum to a Krein
  positive-definite extension, i.e. explore8's dual certificate (ExploreGroupSOS T2), nothing new; (ii) σ-positivity proper fails for
  W on every window: its off-diagonal kernel −k(|u|)/2 − Σ w_n δ_{±log n} is negative near u = 0 and at every prime atom.
- Log-correlated fields / chaos and ζ: Saksman–Webb 1609.00027, Dirichlet L-functions → the same chaos limit (2506.16115), random-model
  maxima (2103.04817, 2404.08575), Ostrovsky 1506.07488 (Bourgade–Kuan–Rodgers zero statistics ↔ chaos total mass, conjectural) and
  1402.6163 (ξ approximated by Mellin transforms of Barnes beta laws): statistical; nothing about where zeros are (R8.2, R36 dead).
- Log-infinite divisibility: Nakamura 1504.03438 (ξ(σ − it)/ξ(σ) is a characteristic function; RH ⟺ "pretended" infinite
  divisibility for 1/2 < σ < 1) and 2306.08317 (RH ⟺ a function is an infinitely divisible characteristic function): restatements.
- Determinants / Polyakov / DOZZ: det′ of flat tori = Epstein zeta (off-line zeros for class number > 1, the Davenport–Heilbronn
  phenomenon), det(Δ − s(1 − s)) of the modular surface contains ζ(2s) (Efrat; resonances ρ/2, Lax–Phillips, Pavlov–Faddeev),
  2112.02771 (det Laplacian of conical spheres ↔ Zamolodchikov's Liouville action), Υ_b zeros lie on lattice cones by construction.
  A self-adjoint Laplacian forces a real spectrum, never RH for its spectral zeta: nothing for N13 (R10.8).
- Closest prior art for §2 T4: Connes 2602.04022 (Letter to Riemann, TeX l. 1156, 1165, 1229): the eigenvectors of QW_λ for its
  "minuscule" eigenvalues are matched numerically by Gram–Schmidt of E(prolate vectors), E(f) lies in the radical of the Weil form,
  and FT(E(h_λ)) → Ξ; CCM 2511.22755 uses the minimal eigenvector. They do not deflate, do not measure the complement and do not
  discuss certificates. Repo: explore8 states "no explicit augmentation ideal exists" (ExploreGroupSOS T1) — corrected here.
- Novelty scan (titles): arXiv math.NT and math.SP listings 2026-10 and 2026-11, show=2000: nothing on deflating Weil window
  certificates; nearest titles 2610.32704 (spectral convergence, assumes RH), 2611.02944 (positivity of Ξ, N6), both seen.

## 2. Transfers, tests, verdicts

### T1 (N1). JSW "positive on small balls": positive piece + norm-dominated remainder. DEAD.
Splits tested (code/run1.py, run1.log; N = 40, both parities, min over parities; exact = λ_min(Q)):
- M0 (fully analytic): exact scale-invariant piece log|t| (λ_min on W_L = 0.34925 − log L even, 1.56376 − log L odd: exact
  scaling, checked c = 1.2 … 100) + inf_t(φ(t) − log|t|) = −3.2044 (t ≈ 0.2) + pole and prime norms.
- M1/M2 (Kahane–Lévy–Khintchine piece ∫₀^L (I − C(x))k dx, PSD on ℝ, + V(L)) with the pole in the remainder.
- M3 (the whole archimedean place, pole included, as the positive piece; primes dominated by the half-overlap bound w_n or by the
  exact Galerkin norm).

| c | half-width | M0 | M1 = M2 | M3 (w_n) | M3 (Galerkin norm) | exact λ_min |
|---|---|---|---|---|---|---|
| 1.1 | 0.048 | −0.50 | +0.84 | +1.02 | +1.02 | 1.02 |
| 1.2 | 0.091 | −1.15 | +0.17 | +0.52 | +0.52 | 0.52 |
| 1.25 | 0.112 | −1.36 | −0.04 | +0.39 | +0.39 | 0.39 |
| 2 | 0.347 | −2.49 | −1.29 | +0.0013 | +0.0013 | 1.3e-3 |
| 2.01 | 0.349 | −3.48 | −1.69 | −0.49 | −0.40 | 1.2e-3 |
| 3 | 0.549 | −3.93 | −2.33 | −0.92 | −0.92 | 5.8e-8 |
| 4 | 0.693 | −5.43 | −3.03 | −1.84 | −1.62 | 8.6e-13 |

Certified half-widths: M0 0.029 (c ≈ 1.06), M1 0.107 (c ≈ 1.24), M3 exactly (ln 2)/2 (c = 2, Connes–Consani/Yoshida) and not one step
further: the first prime atom costs w_2 = 0.49 against λ_min(Q_∞) = 1.3e-3. "Small balls" are smaller than the known window.

### T2 (N1, N13). Kahane σ-positivity of the window covariance; Beurling–Deny structure. DEAD as a proof route (led to T4).
- Exact identity: Q − pole = (Dirichlet form of a symmetric jump process killed outside W: jump density k(|u|)/2 on all of ℝ plus
  atoms of rate w_n at ±log n) − (log π − ψ(1/4) + 2Σ_{n<c} w_n)‖f‖² (the ∫_L^∞ k in V(L) is the killing by jumps longer than L). Markov because Λ(n) ≥ 0 (Davenport–Heilbronn's −f′/f has
  mixed signs). Perron–Frobenius: its ground state is positive and simple (checked: no sign change, c = 1.5 … 4, code/run2.py).
- But the pole-free form is positive only to c ≈ 1.24, and the pole (rank one, positive kernel 2cosh) is anti-Markov in the even
  sector: no Perron–Frobenius argument reaches Q itself.
- Yet the true ground state of Q has no sign change (bulk, |ψ| > 1e-7 max) at c = 1.5, 2, 2.5, 3, 4, 5, 10, 20, so the window Green
  function Q^{−1} ≈ ψ₀ψ₀ᵀ/λ₀ is (numerically) σ-positive. Sturm test (code/run3.py, N = 40, 60 identical; eigenvalues of both
  parities merged): the k-th eigenvector has exactly k sign changes for k ≤ 3 at every c tested (1.5 … 20); beyond that the count
  breaks (k = 4: 6 changes at c = 3, 4; k = 5: 7 at c = 5, 10, 20), so Q^{−1} is not an oscillation kernel (Gantmacher–Krein dead
  beyond order 4). Explanation → T4: ψ₀ ≈ Φ|_W and Φ > 0.

### T3 (N13). Polyakov / zeta-regularised determinants, DOZZ. DEAD (survey §1): every determinant whose zeros contain ρ or ρ/2 is
either ξ itself (Deninger regularised products, circular) or a non-arithmetic spectral zeta with off-line zeros (Epstein).

### T4 (N1; N13 side). The near-null ladder is the cut-off of Weil's unconditional radical. ALIVE (certificate engineering).
Fact (unconditional, classical in Connes's language): Φ(u) = Σ_n (2π²n⁴e^{9u/2} − 3πn²e^{5u/2}) e^{−πn²e^{2u}} has Fourier transform
Ξ, which vanishes at every nontrivial zero, on or off the line; so Q_ℝ(h, Φ^{(j)}) = Σ_ρ H(γ_ρ)(iγ_ρ)^j Ξ(γ_ρ) = 0 for all h
(also for translates and all K∗Φ). On a window only the tails Φ^{(j)}·1_{|u|>L/2} (size ≈ (2πc)^j e^{−πc}) survive:
Q(Φ^{(j)}|_W, h) = −Q_ℝ(Φ^{(j)}1_{W^c}, h). Tests (Arb closed forms → mpmath/float64; code/run4.py, run5.py, run7.py):
- Ground state vs Φ|_W: |cos| = 0.997921 (c = 5), 0.999559 (c = 10); odd ground state vs Φ′|_W: 0.994449, 0.998832.
- Principal angles, bottom-k eigenvectors vs the span of the first k derivatives of the right parity (k ≤ 4, c = 5, 10): k − 1
  cosines ≥ 0.9997, the last 0.90–0.99 (0.97 at c = 10, k = 4).
- Deflation: λ_min of Q on the L²-complement of V_K = span{Φ^{(j)}|_W} (even j for the even block, odd j for the odd block):

| c (N) | block | ladder λ_0 … (rungs < 1e-2) | K = 0 | 1 | 2 | 3 | 4 | … | K = k | complement λ_min |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 (40) | even | 1.2e-17, 6.5e-12, 5.0e-7, 2.5e-3 | ~0 | 6.5e-12 | 4.9e-7 | 2.3e-3 | 0.57 | | 4 | 0.57 |
| 5 (40) | odd | 1.1e-14, 2.0e-9, 6.1e-5 (then 0.10) | ~0 | 1.9e-9 | 5.9e-5 | 8.4e-2 | 0.72 | | 4 | 0.72 |
| 10 (40) | even | 4.7e-43 … 8.0e-3 (9 rungs) | ~0 | … | | | | 8: 3.9e-3 | 9 | 0.29 |
| 10 (40) | odd | 2.1e-39 … 2.0e-4 (8 rungs, then 0.13) | ~0 | … | | | | 8: 7.6e-2 | 9 | 0.34 |
| 20 (60) | even | ≥ 14 rungs below the float64 floor | ~0 | | | | | 16: 8.9e-8, 18: 5.6e-3, 20: 0.085, 22: 0.22 | 24 | 0.31 |
| 20 (60) | odd | ≥ 13 rungs below the float64 floor | ~0 | | | | | 16: 1.1e-6, 18: 3.4e-2, 20: 0.13, 22: 0.24 | 24 | 0.33 |
  Each deflated vector removes exactly one rung (the complement's λ_min equals the next ladder eigenvalue), and once K passes the
  ladder length the remainder is macroscopically positive: K = 4, 9, 24 per parity reach ≥ 0.29 at c = 5, 10, 20 (k(c) ≈ c).
  Stable in N: at c = 10, N = 120 (basis frequencies up to 327) gives 0.40 / 0.57 at K = 10 against 0.41 / 0.58 at N = 40.
  The first rung ratio is ≈ (2πc)² (924 vs 987 at c = 5, 4468 vs 3948 at c = 10), decreasing along the ladder.
- Consequence for certificates (estimate, not yet run): with V = span{Φ^{(j)}|_W}, Q ≥ 0 on the window ⟺ Q_{V⊥} > 0 (margin
  0.3–0.7 here: a double-precision factorisation suffices) and the k × k Schur block S = Q_V − B Q_{V⊥}^{−1} Bᵀ ≥ 0. S's entries are
  ~Q(Φ|_W) ≈ Φ(L/2)² = (2π²c^{9/4}e^{−πc})² up to factors polynomial in c, its bottom is λ_min (ln(1/λ_min)/c = 7.8, 9.75, 11.8 at c = 5, 10, 47; explore8: → 4π), so along the
  bottom direction S needs relative accuracy λ_min/Q(Φ|_W) ≈ e^{−5.5c}: B and Q_V from the tail integrals (no cancellation), and
  Q_{V⊥}^{−1}B by iterative refinement with matrix entries good to ≈ 8c bits (graded error analysis: each ladder direction measured
  against its own size) instead of ≈ 17c bits everywhere. explore8 T1's lemma is not violated: the 17c bits are still resolved, but
  only inside the k × k block S. At c = 47 (ε₄₇ = 1.9e-242): k ≈ 50 per parity, a 53-bit factorisation of the big block plus ~10
  refinement sweeps at ≈ 400 bits, against Chekmachev's order-7170 matrix at ≈ 800 bits (estimate, not run). Untested: the margin of
  Q_{V⊥} in the full space (the lattice method's high-frequency "sporadic bad nodes" lie beyond E_N, N ≤ 120 here) and the rigorous
  error budget.
- N13 side: CCM's missing step one (simple even ground state) is "ψ₀ ≈ Φ|_W > 0 is nodeless", i.e. a property of the ladder block S;
  equally true for Davenport–Heilbronn on its positive windows (Φ_DH is even and ≥ 0, below), so it is not an RH statement.

## 3. Adversarial tests
1. Davenport–Heilbronn type.
   - T1 for DH itself (code/run8.py: odd gamma factor, conductor 5, b_n of −f′/f, no pole): the Kahane piece + norm-dominated
     primes certifies DH up to c ≈ 3.5 (its archimedean part is stronger than ζ's), where DH's exact λ_min is still positive (4.1e-2
     at c = 4, 8.6e-8 at c = 10), far below where its zero 0.8085 + 85.70i can act: a sound bound, no false positivity, and blind
     to ζ versus DH.
   - T2: the Markov structure needs Λ(n) ≥ 0, which DH lacks (b_2 = 0.197, b_3 = −0.312, b_4 = −1.442, …), but it never proved
     anything for ζ either.
   - T4, real DH: Φ_DH(u) = 2e^{3u/2} Σ a_n n e^{−πn²e^{2u}/5} is even (|Φ(u) − Φ(−u)| < 2e-31 max) and positive on [−2.5, 2.5]
     above the 1e-31 noise floor (code/phidh.py): same ladder, so T4 is DH-blind, as a certificate scheme must be.
   - T4, moved zeros (functional-equation symmetric, off-line quartet ±z, ±z̄ added, no Euler product; Q_test = Q + 4 Re F(z)²,
     its own radical P(iD)Φ with P(t) = (t² − z²)(t² − z̄²); code/run6.py):

   | z | c = 5: λ_min(Q_test) | complement of own ladder (K = 2, 4, 6) | c = 10: λ_min | complement (K = 7, 9, 11) |
   |---|---|---|---|---|
   | 10 + 0.3i | −0.025 | −0.020, 0.64, 0.84 | −0.185 | −0.12, 0.46, 0.56 |
   | 30 + 0.3i | −1.5e-5 | −1.4e-5, 0.67, 0.79 | −0.053 | −0.053, 0.37, 0.61 |
   | 60 + 0.3i | +1.1e-17 (not yet seen) | 4.9e-7, 0.56, 0.85 | −1.3e-15 | 2.7e-5, 0.26, 0.54 |

   Failures are detected, appear once the window's ladder reaches the height of the off-line zeros, and always lie inside the
   ladder: after deflation the remainder is positive (0.26–0.85). The remainder is RH-blind; the RH content sits in S.
2. One prime weight perturbed, w_2 → (1 + ε)w_2 (code/run6.py): first-order failure at ε ≈ 3.0e-16 (c = 5) and 8.3e-42 (c = 10)
   (dλ₀/dε = −0.039, −0.057), i.e. along the ladder; the deflated remainder stays at 0.52–0.61 (c = 5) and 0.24–0.35 (c = 10) for
   ε = ±0.1, 0.01. Same location as in test 1.
3. Known windows only? T1 reproduces at most c = 2 (Connes–Consani). T4 is not a window extension: it reorganises how a known
   window is certified; no new window was certified here.

## 4. Ranking (RH relevance × chance × gain)
| transfer | node | RH relevance | chance | gain | verdict |
|---|---|---|---|---|---|
| T4 ladder deflation | N1 (N13) | low: certificate cost, records don't approach RH | medium: works in E_N for c ≤ 10; full-space margin untested | 53-bit factorisation, ≈ 8c-bit assembly, ≈ 17c bits only in a k ≈ c block | ALIVE (engineering) |
| T2 Kahane σ-positivity / Beurling–Deny | N13 | low | nil (pole is anti-Markov; Sturm fails at k = 4–5) | explains CCM's positive ground state only | DEAD |
| T1 JSW small balls | N1 | medium | nil | ≤ c = 2 | DEAD |
| T3 Polyakov / DOZZ | N13 | low | nil | none | DEAD |

Next step for T4 (one run): re-certify c = 16 (Chekmachev: order 133 at 512 bits) with the deflated scheme — Galerkin + tail bound of
Q_{V⊥} in double precision with rigorous rounding, B and Q_V from the tails Φ^{(j)}1_{|u|>L/2} in Arb, S in Arb — and compare cost.

## Suggested NOGO line
- explore17 (Orca, exploration: Gaussian free field, Liouville gravity, Kahane/Junnila–Saksman–Webb positive-type kernels, Polyakov/DOZZ → N1, N13) IDEA (engineering): Kahane pieces + norm-dominated remainder certify only half-width 0.029/0.107/(ln 2)/2 (c = 2 exactly; w_2 = 0.49 ≫ λ_min 1.3e-3), σ-positivity fails (negative jump kernel), W − pole is Markov (Λ ≥ 0) but the pole is anti-Markov; the window ground state is ≈ Φ|_W (cos 0.9996 at c = 10, nodeless; Sturm only for k ≤ 3), and deflating Φ^{(j)}|_W (Weil's unconditional radical, FT Ξ·t^j) removes the near-null ladder rung by rung, leaving λ_min 0.57/0.72 (c = 5), 0.29/0.34 (c = 10), 0.31/0.33 (c = 20, K = 24); off-line-zero and weight-perturbation failures all sit in the k ≈ c ladder block, so certificates need high precision only there (explore8's "no explicit augmentation ideal" corrected); Polyakov/DOZZ/chaos give nothing [ExploreGFF](Riemann/Deep/ExploreGFF/REPORT.md)

## Files (Riemann/Deep/ExploreGFF/)
- REPORT.md (this file); code/gw.py (float64 Galerkin pieces: LK, LOG, V, pole, primes, DH coefficients); code/t_check.py (agreement
  with WeilSimple/weil_matrix.py); code/run1.py (T1 margins); code/run2.py, run3.py (positivity, Sturm counts); code/run4.py (overlaps
  with Φ^{(j)}); code/run5.py, run5_lib.py, run7.py (deflation; exact Φ^{(j)} by the recursion (2z d/dz) on z^a e^{−bz}); code/run6.py
  (adversarial: moved zeros, perturbed w_2); code/run8.py (T1 margins for DH); code/phidh.py (Φ_DH); logs/run1-6.log, logs/run8.log, logs/run7_c20_N60.log, logs/run7_c10_N120.log
  (run2.log's high-precision c = 5, 10 lines are NaN from Arb radius growth in inverse iteration; run3 supersedes them). Uses
  ../WeilSimple/weil_matrix.py (Arb closed forms).
- Run: uv run --no-project --with numpy --with mpmath --with python-flint python3 code/<script> (single-threaded, < 10 min each,
  < 100 MB).

RH relevance: N1 (window-certificate cost; RH content localized in the k ≈ c ladder block), N13 (CCM's ground state ≈ Φ|_W)
