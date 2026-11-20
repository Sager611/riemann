BARRIER: the von Neumann–Wigner rule protects the wrong gap in CCM's step one: the same-sector gap E₂ᵉᵛ/E₁ᵉᵛ (codimension 2) is ≥ 10^4.4 from λ² = 3 on and never binds; the binding gap E₁ᵒᵈᵈ/E₁ᵉᵛ is cross-sector (codimension 1, unprotected) and sets the chain step Δλ² ≈ ln(E₁ᵒᵈᵈ/2E₁ᵉᵛ)/4π ≈ 0.7; odd-sector interlacing E₁ᵒᵈᵈ ≤ E₂ᵉᵛ is certified at 8 windows but has no sector proof, because the Markov part B = Q − pole has a certified parity crossing of its second level at λ² = c× ∈ (5/2, 11/4) (≈ 2.652), and it holds for L(s, χ₅) and Davenport–Heilbronn alike (node N13).

# explore67 (Orca exploration): eigenvalue crossings → N13 (CCM's first missing step)

Orca agent explore67, 20 Nov 2026, 19:08–21:40 UTC. Brief: .scratch/orca/explore67-brief.md. Code: [code/](code/), logs: [logs/](logs/).
Conventions (WeilSimple, CCM 2511.22755 eq. 3.10): window I = [−L/2, L/2], L = ln c, c = λ²; Q = QW_λ = pole + Arch − primes;
pole = 2⟨f, cosh(x/2)⟩² on even f, −2⟨f, sinh(x/2)⟩² on odd f (explore39/57: rank one per sector); B := Q − pole, the Markov part
(killed Lévy process with jump density e^{u/2}/(2 sinh u) plus prime jumps ±ln n of rate Λ(n)n^{−1/2}, Λ(n) ≥ 0; explore17/39).
E_kᵉᵛ, E_kᵒᵈᵈ: k-th eigenvalue of Q on even / odd f; e_k(B) likewise.

## 1. Choice (one tool, one node)
- Tool: interlacing across symmetry sectors in the von Neumann–Wigner classification — two levels of one real-symmetric sector meet
  only in codimension 2, levels of different sectors in codimension 1 — made exact by Aronszajn–Krein for the rank-one pole.
- Node: N13, CCM's step one ("smallest eigenvalue of QW_λ simple with even eigenvector"), i.e. E₁ᵉᵛ < E₁ᵒᵈᵈ and E₁ᵉᵛ < E₂ᵉᵛ.
- Not chosen: (a) N8 via de Bruijn–Newman collisions. A collision is a square-root exceptional point with gap² = 8(t − t_c)
  (explore24); the gap identity d(Δ²)/dt = 8 − 4Δ²S (S > 0) only gives Δ(0)² ≥ Δ(t)² − 8t (Csordas–Smith–Varga), so information flows
  forward in t and no observation at t > 0 prices a collision at t = 0 better than checking t = 0 itself. In the two-parameter real family
  e^{−τD²}cos(hD)Ξ (heat flow × de Bruijn's shift average, both real-zero preserving under RH) a double zero x₀ of Ξ is a diabolical point
  x₀ ± √(2τ + h²) sitting on the corner τ = h = 0 of the admissible region, so no loop encircles it (Berry/Longuet-Higgins sign tests,
  Rellich branches and non-crossing rules are blind to it). The repulsion of the top complex zero by all others (the N2 analogue) is already
  the engine of Raymond's and Hakimoglu-Brown's chains (q′ = −2 − 4hΣ′K_h, DBNHigh, HBDBN).

## 2. Transferred statements (precise)
(S1) Odd-sector interlacing: E₁ᵒᵈᵈ(c) ≤ E₂ᵉᵛ(c) for every c > 1. If true, step one at c ⟺ E₁ᵉᵛ(c) < E₁ᵒᵈᵈ(c), one sector comparison.
(S2) Sector route to (S1) (Aronszajn–Krein + Perron–Frobenius): Qᵉᵛ = Bᵉᵛ + 2|C⟩⟨C| and Qᵒᵈᵈ = Bᵒᵈᵈ − 2|S⟩⟨S| give
  E₂ᵉᵛ ≥ e₂ᵉᵛ(B) and E₁ᵒᵈᵈ ≤ e₁ᵒᵈᵈ(B); so (S1) follows wherever B's second level is odd, e₁ᵒᵈᵈ(B) ≤ e₂ᵉᵛ(B) — the
  "antisymmetric second eigenfunction" property of symmetric Lévy processes killed outside an interval (Bañuelos–Kulczycki 2004 for the
  Cauchy process; stable processes: Kwaśnicki 1012.1133, simplicity 2311.00713; it can fail under perturbation, 2507.19016).
(S3) Hankel route: Qᵉᵛ − Qᵒᵈᵈ = 2H on [0, L/2] with Hankel kernel w(x + y), w = 2cosh(u/2) − ν(u) − Σ Λ(n)n^{−1/2}δ(u − ln n): rank-2
  positive (pole) minus an infinite-rank positive part (ν completely monotone) minus reflections about (ln n)/2 (eigenvalues ±): H has
  infinitely many negative directions, so min–max gives no (S1).

## 3. Tests (logs/: scan_N30.jsonl, scan_N80.jsonl, ritz_*.jsonl, certs.jsonl, bisect.jsonl, adv.jsonl)
Rigour: upper bounds are Arb Rayleigh–Ritz values of the true form on exact dyadic trial vectors in E_N (closed-form Arb entries,
weil_matrix.py; code/ew_ritz.py); lower bounds are WeilSimple's certificate (Chekmachev lattice minorant, Woodbury, Cholesky, rigorous scan
rscan.c) run as code/cert_simple_b.py, which adds one switch WE_NOPOLE=1 (drop the pole column: certify B; in sector defl the deflation
vector is then B's Perron vector). Galerkin scans are floating point at working precision (orientation only).

T1 (S2 falsified, certified). B's second level is odd for small windows and even beyond a parity crossing:
| c = λ² | upper bound | certified lower bound | second level of B |
|---|---|---|---|
| 2 | e₁ᵒᵈᵈ(B) ≤ 0.0904066 | e₂ᵉᵛ(B) ≥ 0.2296 | odd |
| 5/2 | e₁ᵒᵈᵈ(B) ≤ 0.0201624 | e₂ᵉᵛ(B) ≥ 0.03 | odd |
| 11/4 | e₂ᵉᵛ(B) ≤ 0.00926868 | e₁ᵒᵈᵈ(B) ≥ 0.013 | even |
| 3 | e₂ᵉᵛ(B) ≤ 9.66026e−4 | e₁ᵒᵈᵈ(B) ≥ 3.93e−3 | even |
| 4 | e₂ᵉᵛ(B) ≤ 1.17938e−7 | e₁ᵒᵈᵈ(B) ≥ 1.25e−6 | even |
| 5 | e₂ᵉᵛ(B) ≤ 4.41020e−12 | e₁ᵒᵈᵈ(B) ≥ 5.74e−11 | even |
| 7 | e₂ᵉᵛ(B) ≤ 1.93610e−21 | e₁ᵒᵈᵈ(B) ≥ 4.22e−20 | even |
| 10 | e₂ᵉᵛ(B) ≤ 1.77677e−36 (N = 80) | e₁ᵒᵈᵈ(B) ≥ 5.96e−35 (K2 = 401) | even |
Every certificate rigorous (λ_min(I − W) > 0, rigorous scan); every process ≤ 74 MB max RSS (memcap peaks ≤ 20 MB). The 2-dim Ritz bound
is the larger root of det(A − μM) = 0; B's lowest level is the even Perron one (e₁ᵉᵛ(B) ≤ −1.8126 < 0 ≤ E₁ᵒᵈᵈ(Q) ≤ e₁ᵒᵈᵈ(B) at 5/2, Q ≥ 0
certified there), so the second level over all f is odd at 5/2 and even at 11/4: certified parity change in (5/2, 11/4). Sector
eigenvalues are non-increasing in c (min–max); with continuity in c (rescale to a fixed interval; not written out) e₂ᵉᵛ(B) = e₁ᵒᵈᵈ(B) at a
cross-sector crossing (codimension 1, allowed), c× = 2.6519 (Galerkin N = 50; N = 30: 2.6526), level ≈ 0.01978, inside (2, 3) where the only prime term is n = 2. Galerkin
(N = 30 on c ∈ [1.5, 12], N = 80 at 10, 15, 20): odd for c ≤ 2.5, even at every c ≥ 2.75. Mechanism: the jump ±ln 2 rewards even
two-bump profiles with bumps ln 2 apart and penalises their odd partners; without primes (c ≤ 2) the Lévy antisymmetry holds.
So (S2) proves (S1) at most for c < c× and dies there; the 2507.19016 phenomenon (perturbation breaks the second-eigenfunction
property) occurs for ζ's own Markov part, driven by the prime 2.

T2 ((S1) certified at 8 windows): E₁ᵒᵈᵈ ≤ U_odd (Arb Rayleigh quotient) < ε ≤ E₂ᵉᵛ (defl certificate of Q):
| c | 2 | 5/2 | 11/4 | 3 | 4 | 5 | 7 | 10 |
|---|---|---|---|---|---|---|---|---|
| U_odd ≥ E₁ᵒᵈᵈ | 0.0787722 | 1.86351e−3 | 1.57438e−4 | 1.75232e−5 | 5.44116e−10 | 1.29579e−14 | 2.59330e−24 | 9.36624e−40 |
| ε ≤ E₂ᵉᵛ | 0.228 | 0.01197 | 1.49e−3 | 1.62e−4 | 1.0e−8 | 2.99e−13 | 8.86e−23 | 5.1e−38 |
The low spectrum alternates even < odd < even < odd < even at all 35 Galerkin points c ∈ [1.5, 12] (N = 30) and at c = 10, 15, 20
(N = 80): (S1) is true where tested, but by T1 and (S3) it is a ladder property, not a sector one.

T3 (which gap binds; method ceiling of window chains). Galerkin log₁₀ gaps (N = 30 for c ≤ 8; N = 80 at 10, 15, 20):
| c | 1.5 | 2 | 3 | 5 | 8 | 10 | 15 | 20 |
|---|---|---|---|---|---|---|---|---|
| log₁₀(E₂ᵉᵛ/E₁ᵉᵛ), same sector | 1.10 | 2.69 | 4.40 | 5.75 | 6.66 | 7.15 | 7.93 | 7.85 |
| log₁₀(E₁ᵒᵈᵈ/E₁ᵉᵛ), cross-sector | 0.84 | 1.77 | 2.47 | 3.03 | 3.47 | 3.68 | 4.10 | 4.06 |
(Numerical reading: Galerkin values, not certified.) No avoided crossing in the even ladder (no dip of the same-sector gap below 10^4.4 for c ≥ 3; only shallow kinks after prime powers,
e.g. 5.28 → 5.27 across c = 4.25). A chain step c_i → c_{i+1} needs U(c_i) ≈ E₁ᵉᵛ(c_i) below both lower bounds at c_{i+1}; all levels decay
like e^{−4πc} (WeilSimple, CCMFMM), so Δc ≤ ln(min(E₁ᵒᵈᵈ, E₂ᵉᵛ)/2E₁ᵉᵛ)/(4π + o(1)) = ln(E₁ᵒᵈᵈ/2E₁ᵉᵛ)/4π ≈ 0.69 at c = 20 against 1.38
if the same-sector gap bound; BothaChain's 8 steps on [385/12, 227/6] average 0.72. The codimension-2 protected gap never binds; the binding
gap is the unprotected cross-sector one, so non-crossing rules and Hellmann–Feynman gap monotonicity add nothing to monotone chains, and
proving (S1) would only drop WeilSimple's defl certificate (half the work per window), not shorten steps. Lipschitz (Hellmann–Feynman) bounds
on ln E_k in c give the same step, since both sides decay at the same rate 4π.

T4 (adversarial; ew_adv.py with explore19's wmat.py, Galerkin N = 80, N = 90 at c = 31): (S1) holds for L(s, χ₅) at c = 3, 5, 10, 20
(E₁ᵒᵈᵈ/E₂ᵉᵛ = 0.71, 0.51, 0.16, 0.018) and for Davenport–Heilbronn at c = 3, 5, 10, 20 and at its failing window 31
(E₁ᵉᵛ = −6.52e−30 < E₁ᵒᵈᵈ = 7.25e−26 < E₂ᵉᵛ = 3.29e−22), forms with no pole and signed prime weights (no Markov part at all). So (S1) carries
no RH information, like step one itself (explore32). Full alternation is not universal: DH at c = 5 has E₂ᵒᵈᵈ = 1.599 < E₂ᵉᵛ = 1.647.
What tells ζ apart is only the route: ζ alone has a Markov part (Λ(n) ≥ 0), and even that part loses the Lévy parity order at c×.

## 4. Result against the node line (TREE N13)
No new bound or range: the certified windows lie inside WeilSimple's (step one is already certified on [2, 227/6] ∪ {40}). New: (i) a
certified parity crossing of the Markov part's second level at c× ∈ (5/2, 11/4), which kills the Perron–Frobenius/Aronszajn–Krein route
to odd-sector interlacing beyond c×; (ii) (S1) certified at 8 windows (RH-blind); (iii) the identification of the binding gap: CCM's step
one is limited by a codimension-1 (cross-sector) gap, so the von Neumann–Wigner rule, which protects only same-sector gaps, cannot certify
it, and window chains cannot beat Δλ² ≈ ln(E₁ᵒᵈᵈ/2E₁ᵉᵛ)/4π per certificate. Unconditional (no RH used).

## 5. Novelty (checked 20 Nov 2026)
arXiv API: "avoided crossing" ∧ zeta, "level repulsion" ∧ "Riemann zeros", "von Neumann-Wigner" ∧ zeros, "Hellmann-Feynman" ∧ zeta,
"level curvature" ∧ zeta, "Weil quadratic form" ∧ eigenvalue/parity, "second eigenfunction" ∧ Cauchy/Lévy: nothing on Weil's form
(2611.13129, 2610.32704, 2607.24830, 2605.20224 seen; 2507.19016 is the generic perturbation counterexample). lst.py (math.NT, math.SP,
this and last month; Zenodo '"Weil quadratic form" parity'): no match. Repo: explore39 (ExploreFunctional t7_parity) compares even and odd
GROUND states via a reflection identity for c ≤ ρ²; explore57 tabulates Galerkin B spectra (its c = 2 and 5 rows show the two parities,
unremarked); explore32 noted Q's alternation. The crossing c×, its certificate and the binding-gap reading are not in NOGO/TREE/REPORTs.

## Suggested NOGO line
- explore67 (Orca, exploration: eigenvalue crossings — von Neumann–Wigner sector interlacing with the Aronszajn–Krein pole → N13) BARRIER: CCM's step one is bound by the cross-sector gap E₁ᵒᵈᵈ/E₁ᵉᵛ (codimension 1, unprotected; 10^0.8–10^4.1 on λ² ∈ [1.5, 20]), not by the codimension-2 same-sector gap (≥ 10^4.4 from λ² = 3, no avoided crossing), so window chains cannot beat Δλ² ≈ ln(E₁ᵒᵈᵈ/2E₁ᵉᵛ)/4π ≈ 0.7 (BothaChain 0.72); odd-sector interlacing E₁ᵒᵈᵈ ≤ E₂ᵉᵛ holds (certified at λ² = 2, 5/2, 11/4, 3, 4, 5, 7, 10; also χ₅, DH at its failing window: RH-blind) but has no sector proof: the Markov part Q − pole has a certified parity crossing of its second level at λ² ∈ (5/2, 11/4) (≈ 2.652, prime 2 breaks the Lévy antisymmetry), and the parity Hankel difference is indefinite; dBN collision pricing reduces to Csordas–Smith–Varga (a double zero of Ξ is a corner diabolical point x₀ ± √(2τ + h²) of e^{−τD²}cos(hD)Ξ) [ExploreWigner](Riemann/Deep/ExploreWigner/REPORT.md)

## Files and footprints
- code/: weil_matrix.py (WeilSimple's, unchanged), wmat.py (explore19's, unchanged), cert_simple.py → cert_simple_b.py (WeilSimple's
  certificate + the WE_NOPOLE switch, 10 lines marked EW), ew_ritz.py (Arb Ritz upper bounds), ew_scan.py (Galerkin sector spectra of Q
  and B), ew_bisect.py (c×), ew_adv.py (ζ, χ₅, DH), run_certs.sh + certs_list*.txt (the certificate queue).
- logs/: scan_N30.jsonl (35 windows), scan_N80.jsonl, ritz_a/b.jsonl, ritz_c10_N80.jsonl, certs.jsonl (18 runs, 16 rigorous; the two
  c = 10 runs with targets from unconverged N = 30 values failed and were redone with N = 80 targets, K2 = 401), bisect.jsonl, adv.jsonl.
- One compute process at a time, OMP_NUM_THREADS=1, every run under .scratch/memcap.sh 290: memcap peaks ≤ 20 MB, max RSS ≤ 74 MB,
  longest run 85 s (scan_N30). Rerun: bash Riemann/Deep/ExploreWigner/code/run_certs.sh Riemann/Deep/ExploreWigner/code/certs_list.txt
  (needs .scratch/explore67/rscan = cc -O2 Riemann/Deep/WeilSimple/rscan.c -lm).
- Reviewer (verdict only, inline evidence): GAP — (1) and the certificates correct; parity at 5/2 needs e₁ᵉᵛ(B) < e₁ᵒᵈᵈ(B) (added
  above); the exact crossing needs continuity in c (stated as such); which gap binds is a Galerkin reading, not certified (labelled).
- Not done: certificates beyond c = 10; a lower bound for e₂ᵉᵛ(B) on (2, c×) at more than two windows.

RH relevance: N13
