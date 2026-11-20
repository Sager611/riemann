BARRIER: Pólya–Turán real-zero transfer into orthogonal bases cannot certify hyperbolicity of Ξ's approximants — Romik's Meixner–Pollaczek pre-images fail Newton's inequality at every order N ≥ 2 for every positive Φ (ζ: c₂²/c₀c₄ = 0.386 < 1.4), Turán's Hermite pre-image is the de Bruijn–Newman function at time −1 (dead by Λ ≥ 0), and the MP, Hermite and continuous-Hahn approximants have non-real zeros from N = 2 (Arb), as do Laguerre–Pólya toys and Davenport–Heilbronn; Jensen's damping is Maclaurin-extremal (node N6).

# explore63 — special-function positivity → N6 (Jensen hyperbolicity), N2 link

Status: done 20 Nov 2026, 13:37–15:00 UTC. Scripts here, outputs in logs/. All peaks ≤ 40 MB (memcap lines in the logs).

## Transfer

- Field tool: Pólya–Turán real-zero preservation for sums of classical orthogonal polynomials, and positivity of expansion coefficients.
  - Turán (1959, Lemma II, after Pólya 1915): Σ a_k z^k real-rooted ⇒ Σ a_k H_k(z) real-rooted.
  - Romik (1902.06330, §open problems): Σ a_k z^k real-rooted ⇒ Σ (k!/(3/2)_k) a_k f_k(z) real-rooted, f_k = P_k^{(3/4)}(·; π/2) (Meixner–Pollaczek).
  - Romik's expansions: Ξ(t) = Σ (−1)^n c_{2n} f_{2n}(t/2), c_{2n} > 0, c_{2n} ~ 16√2 π^{3/2} √n e^{−4√(πn)}; Ξ(t) = Σ (−1)^n b_{2n} H_{2n}(t).
- Node: N6 (finite hyperbolicity tests for Ξ). The Jensen polynomials are the Taylor-basis case; the orthogonal bases give the Hermite and MP analogues.
- Transferred claims (precise): C1 the MP–Pólya sufficient condition (real-rooted pre-image ⇒ real-rooted approximant) can hold for no N ≥ 2; C2 the Hermite one needs Ξ_{−1/4} ∈ LP (Romik units); C3 the approximants Σ_{n≤N} (−1)^n e_{2n} p_{2n}(t/2) of Ξ in the MP, Hermite and continuous-Hahn bases are hyperbolic for all N (tested, false); Lemma M Jensen's multiplier is extremal.
- Why not (a) N3 or (b) N1: see "Other targets" below.

## Tests (ζ unless stated; scripts here, outputs in logs/)

- Setup checks (logs/mp_lemma_check.log): recurrence (n+1)f_{n+1} = 2x f_n − (n+½) f_{n−1} equals Romik's ₂F₁ definition (err 9e−40); Romik's lemma holds on 200 random real-rooted p (deg 2–10): all images real-rooted; control 1+z² ↦ roots ±0.75i.
- c_{2n}, n ≤ 200 (logs/mp_c2n_zeta.log): two Gauss–Legendre rules agree to 4.9e−45; c₀, c₂, c₄ = 0.480990005928, 0.0194182171985, 0.00203194488774; c_{2n}/(Romik asymptotic) = 0.49, 0.71, 0.78, 0.84 at n = 10, 50, 100, 200 (O(n^{−1/10}) error term); Σ_{n≤200} (−1)^n c_{2n} f_{2n}(t/2) = Ξ(t) to 12 digits at t = 0, 5, γ₁, γ₂, 30, 40.

### C1 — the MP–Pólya route fails at every order N ≥ 2, for every positive Φ (theorem)
- Statement. Let Ξ_Φ(t) = ∫ Φ(u) e^{itu} du, Φ ≥ 0 even, = Σ (−1)^n c_{2n} f_{2n}(t/2). For every N ≥ 2 and every Poisson factor r ∈ (0, 1], G_N(z) = Σ_{n≤N} (−1)^n r^{2n} c_{2n} ((3/2)_{2n}/(2n)!) z^{2n} has non-real zeros, so Romik's lemma never certifies a real-rooted approximant.
- Proof. c_{2n} = ∫ tanh^{2n}(u) cosh^{−3/2}(u) Φ(u) du (u > 0) are Hausdorff moments, so c₂² ≤ c₀c₄ (Cauchy–Schwarz). G_N is a polynomial of degree N in w = z² with coefficients A_n = (−1)^n r^{2n} c_{2n} β_n, β_n = (3/2)_{2n}/(2n)!; real-rooted in z needs N positive roots in w, hence Newton's inequality at k = 1: A₁² ≥ A₀A₂ · 2N/(N−1). With β₁²/β₂ = (15/8)²/(945/384) = 10/7 this reads c₂²/(c₀c₄) ≥ (7/5) N/(N−1) > 1. Contradiction; r cancels.
- ζ: c₂²/(c₀c₄) = 0.3858 against the needed ≥ 1.4.

### C3 — the MP approximants themselves are not hyperbolic (Arb root isolation)
- Ξ_N(t) = Σ_{n≤N} (−1)^n c_{2n} f_{2n}(t/2), degree 2N; all roots isolated by python-flint acb_poly.roots (c-balls of relative radius 1e−44); non-real = imaginary ball excludes 0.
- ζ (logs/mp_approx_roots_zeta.log): non-real zeros for every tested N ≥ 2: N = 2: 4 of 4 (t = ±8.44 ± 2.32i); N = 5: 8/10; 20: 16/40; 50: 24/100; 100: 40/200 (min |Im t| 1.06). Real zeros reproduce γ₁ = 14.13473 from N = 20 and γ₂ = 21.02204 from N = 50 (21.02223 at N = 30); spurious real zeros reach t ≈ 3.8N. (N = 150, 200 at 3000 bits: see the bounded-scale section.)
- LP toys (logs/mp_approx_roots_toys.log; tested N = 2, 3, 5, 10, 20, 40), all with only real zeros: cos(t/2) (Φ = atoms at ±½): non-real from N = 5; sin t / t (Φ = 1 on [−1, 1]): from N = 10; e^{−0.045t²} (Gaussian Φ): from N = 2.
- So MP truncation does not preserve the Laguerre–Pólya class: "Ξ_N hyperbolic for all N" is not implied by RH and is false for ζ at N = 2.

### C1 extended to the whole Meixner–Pollaczek family (logs/mp_lambda.log; mp_lambda.py)
- P_n^{(λ)}(x; π/2), (n+1)P_{n+1} = 2xP_n − (n+2λ−1)P_{n−1}; lemma Σ a_k z^k real-rooted ⇒ Σ (k!/(2λ)_k) a_k P_k^{(λ)} real-rooted: 100 random real-rooted p per λ ∈ {0.25, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 5, 10, 30}, all images real (max |Im| = 0). Romik states only λ = 3/4.
- Expansion: cos(tu) = cosh^{−2λ}(u) Σ (−1)^n tanh^{2n}(u) P_{2n}^{(λ)}(t/2), so c_{2n}^{(λ)} ∝ ∫₀^∞ Φ(u) cosh^{−2λ}(u) tanh^{2n}(u) du: Hausdorff moments for every λ, dilation and Poisson factor.
- Newton at k = 1 forbids a real-rooted pre-image p_N unless R := c₂²/(c₀c₄) ≥ ν(λ)·2N/(N−1), ν(λ) = (2λ+2)(2λ+3)/(12λ(2λ+1)).
  - λ < 1.151 (2ν > 1 ≥ R by Cauchy–Schwarz): fails at every N ≥ 2 for every positive Φ.
  - ζ, λ = 1.25 … 30: R = 0.3834 … 0.3445 < 2ν = 0.943 … 0.356: fails at every N ≥ 2.
  - λ → ∞, any smooth even Φ > 0: R = ⅓(1 + 1/λ + O(λ^{−2})) (Φ''(0) cancels; computed from Gaussian moments) < 2ν = ⅓(1 + 2/λ + O(λ^{−2})).
- Arb: p_N has non-real zeros for every λ on the grid and every N = 2…30 (exactly 2N or 2N − 2 of 2N: at most two real zeros).
- All N large, every λ: c_{2n}^{1/2n} → 1 (Φ > 0 on all of ℝ puts mass at α → 1), so the pre-image has radius of convergence 1 and Erdős–Turán's discrepancy bound spreads the zeros of p_N over the unit circle for every large N.

### C2 — Turán's Hermite route (logs/herm_route.log; herm_route.py)
- Known: Romik's Thm (Hermite–Poisson flow) X_r(t) = Ξ_{(r²−1)/4}(rt), and the pre-image is e^{D²/4}Ξ = Ξ_{−1/4}, the r → 0 end of the rescaled flow X_r(t/r): the Pólya–de Bruijn function at λ = −1/4 in Romik's units (dBN time −1); every dilation s gives time −s²/4 < 0. Rodgers–Tao (Λ ≥ 0): never Laguerre–Pólya.
- Numbers (ζ): b₀, b₂, b₄ = 0.2457, 1.391e−3, 3.669e−6 (Ξ/2 normalization); Newton k = 1 passes for N ≥ 15 (b₂²/(b₀b₄) = 2.145), yet Arb finds non-real zeros of the pre-image p_N at every tested N (2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40), and the Hermite approximants Σ_{n≤N} (−1)^n b_{2n} H_{2n}(t) are far from hyperbolic: real/non-real = 0/4 (N = 2), 0/20 (10), 4/36 (20), 8/72 (40); their real zeros are Ξ's (14.1347, 21.0220, 25.0109).

### Lemma M — Jensen's damping is extremal among finite multiplier sequences (Maclaurin)
- Statement. Let E(z) = Σ_{k≤d} γ_k z^k/k! be real-rooted of degree d with γ_k ≥ 0 and γ₀ = γ₁ = 1 (by Pólya–Schur, the finite multiplier sequences of the first kind, up to z ↦ cz and z ↦ −z). Then γ_k ≤ (d)_k/d^k = ∏_{j<k}(1 − j/d) ≤ e^{−k(k−1)/(2d)} for every k, with equality for all k iff E = (1 + z/d)^d.
- Proof. γ_k = k! e_k(1/a₁, …, 1/a_d) for the roots −a_i; Maclaurin's inequality (e_k/C(d,k))^{1/k} ≤ e₁/d = 1/d gives γ_k ≤ k! C(d,k) d^{−k}.
- Reading for N6. J^{d,0}(x/d) = Σ_k ((d)_k/d^k) γ(k) x^k/k! is exactly the equality case. Every degree-d test that is real-rooted for every Laguerre–Pólya input (the multiplier, i.e. diagonal, class) damps the k-th Taylor coefficient of Ξ at least as much as Jensen's, so it reaches no further than Jensen's √d.
- Not covered: non-diagonal real-rootedness preservers (Borcea–Brändén symbols). The Hermite and MP transforms are two of them (C1, C2).

### Resolution law of the orthogonal approximants (MP: quadratic like Jensen's; continuous Hahn: N^{2/3})
- Romik's c_{2n} ~ e^{−4√(πn)} and |P_{2n}(t/2)| ≈ e^{π|t|/4} against |Ξ(t)| ≈ e^{−π|t|/4}: Ξ_N ≈ Ξ needs 4√(πN) > π|t|/2, i.e. |t| < (8/√π)√N = 4.51√N.
- ζ, measured: the non-real zero closest to the axis (foot of the non-real arc) sits at |Re t| = 8.44, 12.06, 14.67, 21.43, 26.09, 32.32, 39.94, 46.81 for N = 2, 5, 10, 20, 30, 50, 75, 100, against 4.51√N = 6.4, 10.1, 14.3, 20.2, 24.7, 31.9, 39.1, 45.1; inside, the real zeros of Ξ_N are Ξ's (14.13473, 21.02204 from N = 30).
- Independent check (logs/check_approx.log): mpmath polyroots at 300 digits gives the same counts (N = 10: 12 real, 8 non-real; N = 20: 24 real, 16 non-real, non-real at (±14.4513, ±18.6656), (±17.886, ±11.609), (±20.0418, ±6.3253), (±21.4327, ±1.8543)); sign changes on a 0.01 grid give 24 (N = 20) and 76 (N = 50) real zeros; |Ξ₂₀| = 1e−298 at the refined non-real root.
- Davenport–Heilbronn: the same saddle point with Φ_DH ~ exp(−(π/5)e^{2u}) gives c_{2n} ≈ e^{−4√(πn/5)}, so |t| < 2.02√N (heuristic); N = 200 reproduces F(½+it) for t ≤ 20, fails at t = 30 (logs/dh_c2n.log); seeing the first off-line zero ρ = 0.808517182457 + 85.6993484854i (reproduced by findroot) needs N ≈ 1800 by the same rule.
- So MP approximants buy nothing in range (a zero at height t costs N ≈ 0.049t², the T² law of N3/N6 again). Continuous Hahn does better (N ≈ 0.20t^{3/2}, next section) and Taylor partial sums are linear, but only the damped Jensen scheme stays real-rooted.

### C3 at bounded scale and in Romik's third basis
- MP, ζ, c_{2n} to 3.6e−100 (two rules, logs/mp_c2n_zeta_hi.log), 3000-bit isolation (logs/mp_approx_roots_zeta_hi.log): N = 150: 252 real, 48 non-real, foot at |Re t| = 56.50 (edge 55.2); N = 200: 344 real, 56 non-real, foot 65.08 (edge 63.8). N = 300 stopped (single compute slot).
- Continuous Hahn g_n (orthogonal for |Γ(¾+ix)|⁴; (2n+3)g_{n+1} = 8x g_n − (2n+1)g_{n−1}), Romik's d_{2n} (₂F₁ kernel), n ≤ 100, two rules agree to 2.6e−40; Σ_{n≤100} (−1)^n d_{2n} g_{2n}(t/2) = Ξ(t) to 12 digits at t ≤ 50 (logs/hahn_d2n.log).
- Its approximants (logs/hahn_roots_zeta.log): non-real zeros at every N tested: N = 2: 4/4; 10: 16/20; 20: 28/40; 50: 52/100; 75: 72/150. Foot of the non-real arc at |Re t| = 23.94, 30.45, 42.04, 55.35 for N = 20, 30, 50, 75, against the edge 2.96 N^{2/3} = 21.8, 28.6, 40.2, 52.6 predicted from d_{2n} ~ e^{−3(4π)^{1/3} n^{2/3}}. Real zeros = Ξ's (14.13473, 21.02204, 25.01086, 30.42488 at N = 50).
- LP toy cos(t/2) in the same basis (logs/hahn_roots_atom.log): non-real from N = 3 (N = 40: 24 real, 56 non-real).

## Adversarial test (Davenport–Heilbronn; logs/dh_c2n.log, logs/mp_approx_roots_dh.log)
- DH as F(s) = (5/π)^{(s+1)/2} Γ((s+1)/2) f(s), a(n) = 1, κ, −κ, −1, 0 (mod 5), κ = 0.284079043840412; F(½+it) real (|Im/Re| < 3e−40); Φ_DH(u) = 2θ_a(e^{2u})e^{3u/2} even (ratio 1.0 at u = 0.2, 0.5, 0.9) and positive: on [0, 3.6] numerically (Φ_DH(0) = 1.1527), and for all u because the n = 1 term dominates θ_a(x) for x ≥ 1 (the negative terms sum to < 0.6 % of it at x = 1, less beyond), then evenness.
- So DH is also a positive-definite even function: C1 applies verbatim (R = c₂²/(c₀c₄) = 0.4689 < 1.4), and its MP approximants are non-hyperbolic from N = 5 (N = 2, 3 real-rooted; N = 100: 184 real, 16 non-real; real zeros 5.09416, 8.93991 = DH's first critical zeros).
- Verdict: every statement above holds for ζ, DH and LP toys alike. It is a statement about the method (positive Φ, subexponential coefficient decay), not about ζ — as a barrier must be.

## Other targets (why not chosen)
- (a) N3, archimedean part of λ_n: S_∞(n) of Bombieri–Lagarias is explicit, (n/2) log n + O(n); λ_n > 0 is decided by the prime/zero part. A positive Laguerre combination φ gives the prime sum a sign only through Λ ≥ 0, and the pole term ∫φ cancels it to PNT accuracy (pole-anchored, R30/R31). Computing λ_n arithmetically costs Θ(n) digits (Johansson 10⁵, R118).
- (a′) Schoenberg reading of Li (standard, not new): RH ⟹ K_{mn} = λ_m + λ_n − λ_{|m−n|} ⪰ 0, because K is Weil's Gram matrix on Li's test functions 1 − (1 − 1/s)^n (F_m F̃_n = (1 − w^m)(1 − w^{−n})). This is the right PSD object; R395/R482 tested the Toeplitz matrix of λ, which has zero diagonal. It is generic CND, true for any real ordinates (R433).
- (b) N1: Legendre/Jacobi bases do not diagonalise Weil's translation-invariant form on a window. Positive-type kernel decompositions were explore17 (Kahane) and explore27 (design/minorant LP).

## Novelty check (20 Nov 2026)
- arXiv API: "Meixner-Pollaczek" ∧ (xi ∨ zeta ∨ Riemann): 5 hits. 1902.06330 (Romik: the lemma and the open question; the Hermite–Poisson = de Bruijn–Newman theorem is his). 1412.1220 (Inoue: MP partial sums of a meromorphic completed ζ, "almost all zeros real, four complex near the origin" at degree 260 — a different function; ours is entire Ξ with ≈ 0.4N non-real zeros on an arc at the resolution edge). 2512.11405 (expansions of η, ζ in hypergeometric polynomials with zeros on the line; nothing on approximant zeros). 2605.26846 (Stirling polynomials). 0805.0725 (physics).
- arXiv API, no hits bearing on C1–C3 or Lemma M: "orthogonal polynomial" ∧ "Riemann xi" (+2002.12825, physics); "Jensen polynomials" ∧ (Maclaurin ∨ extremal ∨ optimal); "multiplier sequence(s)" ∧ (Jensen ∨ Hermite ∨ Laguerre ∨ Legendre ∨ Meixner ∨ Jacobi) (Legendre/Laguerre/Hermite multiplier-sequence papers 1002.0759, 1108.4662, 1210.7151, 1312.6187, 1404.1631, 1701.02420, 2611.17839: classifications of basis multipliers, nothing on Ξ); "Meixner-Pollaczek" ∧ "real zeros" (2604.25692: interlacing); "Hausdorff moment" ∧ Riemann ∧ xi; Romik ∧ Turán; "partial sums" ∧ "orthogonal polynomial" ∧ "real zeros": 0.
- Monthly listings 2026-10, 2026-11 (math.NT/CA/CV/SP/PR/FA/CO), title regex Meixner|Pollaczek|continuous Hahn|Hermite expan|orthogonal polynomial.*(zeta|xi)|Turán.*(Hermite|xi)|multiplier sequence|Jensen polynomial: 2611.01422, 2611.13932, 2610.15226 (Meixner families; not ζ).
- Zenodo: '"Meixner-Pollaczek" zeta', '"Riemann xi" "orthogonal polynomials"', 'Romik Turan Hermite expansion xi': nothing on Pólya–Turán–Romik routes. Read the abstracts of 24377750/57 (Thorin programme: GGC positivity, not orthogonal expansions), 23941908 (Gegenbauer Hankel determinants for quadrature), 23960915 (claimed RH proof via an approximated ξ, unrelated mechanism).
- Our NOVELTY/NOGO files: no Askey/Gasper/Vietoris/Romik/Meixner entry; Schoenberg only as generic CND (R419, R433).
- Classics: Pólya 1915 and Turán 1959 (Hermite lemma); Maclaurin 1729 (Lemma M is a one-line consequence, likely folklore); Erdős–Turán 1950 (zero discrepancy); Rodgers–Tao 2018 (Λ ≥ 0).

## Verdict
- Barrier, unconditional (C1, C3, Lemma M need no hypothesis; C2 uses Rodgers–Tao). No new bound on a node quantity: the N6 record (d ≤ 1.6144·10²⁵ for all n) is untouched.
- Where the transfer breaks.
  - Pólya's positivity of Φ makes every MP coefficient sequence a Hausdorff moment sequence; Cauchy–Schwarz then breaks Newton's inequality at the first step of every pre-image (C1). Positivity, the field's own tool, is what kills the route.
  - Λ ≥ 0 leaves Ξ no heat-flow room, so Hermite pre-images are never Laguerre–Pólya (C2).
  - Truncations beyond the resolution edge put an arc of spurious zeros off the axis (C3), for Laguerre–Pólya inputs as well. Only damped (multiplier) schemes stay real-rooted, and Jensen's is the least damped (Lemma M).
- Method ceiling, freedoms covered: basis (Hermite, MP for every λ, continuous Hahn), dilation, Poisson factor r, order N. Not covered: other non-diagonal real-rootedness preservers (Borcea–Brändén symbols).
- Certification status. C1 at λ = 3/4: proof (any positive Φ). In λ: proof for λ < 1.151 and asymptotic for λ → ∞; between them a grid (λ = 1.25 … 30, margin ≥ 3 %), not a continuum certificate; every λ for N ≥ N₀(λ) by Erdős–Turán. C3 for ζ: Arb isolation with coefficient balls covering the two-rule quadrature difference, plus an independent 300-digit check; λ-grid and Hermite runs use 45–55-digit point values (their non-real zeros lie ≥ 0.5 from the axis).

| scheme | sees zeros up to | real-rooted for all LP inputs | RH-necessary | ζ outcome |
|---|---|---|---|---|
| Jensen J^{d,n} | d ≍ t² | yes (Pólya–Schur) | yes | hyperbolic, d ≤ 1.6·10²⁵ (record) |
| other finite multipliers | no further than Jensen (Lemma M) | yes | yes | — |
| MP approximants, degree 2N | t ≈ 4.51√N | no (toys from N = 2–10) | no | non-real from N = 2 |
| continuous-Hahn approximants | t ≈ 2.96 N^{2/3} | no (cos from N = 3) | no | non-real from N = 2 |
| Hermite approximants | t ≈ 27 at N = 40 | no | no | 8 real, 72 non-real at N = 40 |
| MP–Pólya pre-image (sufficient) | — | — | sufficient only | Newton fails at every N ≥ 2 |
| Hermite–Pólya pre-image (sufficient) | — | — | needs Λ ≤ −1 | false (Rodgers–Tao) |

## Suggested NOGO line
- explore63 (Orca, exploration: special-function positivity — Pólya–Turán real-zero transfer into Hermite, Meixner–Pollaczek and continuous-Hahn bases (Turán 1959, Romik 1902.06330) → N6, N2) BARRIER: Romik's MP pre-images fail Newton at every order N ≥ 2 for every positive Φ (Hausdorff moments: c₂²/c₀c₄ ≤ 1 < 1.4; ζ 0.386, DH 0.469; every λ < 1.151, grid to λ = 30, and λ → ∞); Turán's Hermite pre-image is de Bruijn–Newman at time −1 (Λ ≥ 0); the MP, Hermite and continuous-Hahn approximants have non-real zeros from N = 2 (Arb; the arc's foot sits at the resolution edge 4.51√N, resp. 2.96N^{2/3}; LP toys and DH alike); Jensen's damping (d)_k/d^k is Maclaurin-extremal among finite multipliers, so N6's T² law stands [ExplorePositivity](Riemann/Deep/ExplorePositivity/REPORT.md)

## Files
- mp_core.py (f_n, ω, c_{2n}, Ξ); mp_lemma_check.py; mp_c2n.py; mp_approx_roots.py; check_approx.py; mp_lambda.py; herm_route.py; dh_c2n.py; hahn_d2n.py; hahn_roots.py; hahn_toys.py; logs/ (coefficient tables c2n_*.txt, d2n_*.txt and every run's output with its memcap peak).

RH relevance: N6 (and N2, through Romik's Hermite–Poisson = de Bruijn–Newman identity).
