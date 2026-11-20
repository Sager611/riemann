NO IDEA: Schrödinger tools re-encode known machinery: the Prüfer angle of Z is Levinson–Montgomery's arg ζ′ count, Birman–Schwinger/Weinstein–Aronszajn are explore39's pinning and the N1 record's own dual check, and phase-space counts 2 ln c·E e^{X_c} are vacuous (side law: the record's BS rank follows that count, ×0.87–1.55 for c = 10…71).

# ExploreBS — explore57: spectral theory of Schrödinger operators → N8, N13, N1 (Orca exploration, 20 Nov 2026)

Field: Birman–Schwinger, Weinstein–Aronszajn/Aronszajn–Donoghue, bound-state counting (Bargmann, Calogero, CLR,
Lieb–Thirring, Berezin–Lieb), Sturm oscillation/comparison and Prüfer angles. Repo check (rg, NOGO + Deep reports):
0 hits for Birman–Schwinger/Weinstein–Aronszajn/Prüfer; Lieb–Thirring/CLR only in R307/R454 (on −Ξ″/Ξ), Sturm in R371,
ExploreGFF (nodal counts k ≤ 3). Scripts: `code/`, outputs: `logs/`.

## 0. Tool → node map (which transfer is new, which re-encodes earlier work)

| tool | node | what it becomes for ζ | earlier work it reduces to |
|---|---|---|---|
| (d) Prüfer angle of Z w.r.t. the θ-oscillator | N8 | φ = 3π/2 − θ − arg ζ′(½+it), exact; crossings = zeros of Z | Levinson–Montgomery 1974, Berndt 1970 (argument principle for ζ′) |
| (d) Prüfer frequency ω | N8 | Z + iZ′/ω = −(L/ω)e^{iθ}(ζ + ζ′/L), L = θ′ − ω: Levinson's V | Levinson 1974 (L = log(T/2π) constant) |
| (d) Sturm comparison per gap | N8 | Picone/Wirtinger per gap | explore53, gapjinf (per-gap class) |
| (a)+(b) BS / W–A, rank-one pole | N13, N1 | Q^ev = B^ev + 2cc*, Q^odd = B^odd − 2ss*; interlacing pins μ₁(B^ev), μ₀(B^odd) within the ladder | explore39 (Lévy pinning), explore17/weildefl (deflation) |
| (a) BS at an off-line point | N1 | Christoffel function of on-line zeros at a = γ − iβ | explore38, explore19 (detection law), morirate |
| (c) Berezin–Lieb / CLR / Weyl | N1 | Σλ₋ ≤ (ℓ/2π)∫σ₋ for every admissible symbol; = 0 iff a nonnegative Krein extension exists | the record: a lattice symbol whose negative nodes a BS check shows not to bind (Chekmachev) |
| (c) Bargmann, Calogero, GGMT | — | need a local operator −d² + V; Weil's form is nonlocal | — |

Chosen transfer for the tests: (d) → N8, because it is the only one not already in the repo; (c) → N1 gets a
quantitative vacuity test; (a)/(b) are recorded as reductions (§3).

## 1. Transfer (d): the Prüfer angle of Hardy's Z → N8

Setting. For t > t* = 6.2898 (θ′(t*) = 0, θ(t*) = −3.531; θ′ > 0 beyond) the θ-oscillator (y′/θ′)′ + θ′y = 0 has solutions cos θ, sin θ. Hardy's Z obeys
(Z′/θ′)′ + θ′Z = R(t) with R the Riemann–Siegel remainder beyond the n = 1 term. Prüfer variables: Z = r sin φ,
Z′/θ′ = r cos φ, φ continuous. Statements (exact unless marked):

- (P1) Z + iZ′/θ′ = −e^{iθ}ζ′(½+it)/θ′ (from Z′ = ie^{iθ}(θ′ζ + ζ′)), hence φ = 3π/2 − θ(t) − arg ζ′(½+it) (mod 2π,
  arg continuous along the line; valid while ζ′(½+it) ≠ 0, i.e. no multiple zero of ζ on the line).
- (P2) φ′ = (θ′(Z′² − ZZ″) + θ″ZZ′)/(Z′² + θ′²Z²) = −θ′ − Re(ζ″/ζ′)(½+it); at a zero of Z, φ′ = θ′ > 0, so φ crosses each
  kπ upward exactly once: #zeros of Z (odd multiplicity) in (a, b] = ⌊φ(b)/π⌋ − ⌊φ(a)/π⌋. With (P1):
  N₀(T) = −(θ(T) + arg ζ′(½+iT))/π + O(1). This is the argument principle for ζ′ on the critical line, i.e. the
  Levinson–Montgomery (1974)/Berndt (1970) count: no new information.
- (P3) Prüfer frequency ω ≠ 0 (any function of t): Z + iZ′/ω = e^{iθ}[(1 − θ′/ω)ζ − ζ′/ω] = −(L/ω)e^{iθ}(ζ + ζ′/L) with
  L = θ′(t) − ω(t). The count (P2) holds for every ω. (P1) is L = 0 (V = ζ′); Levinson's V = ζ + ζ′/L has L = log(T/2π)
  constant, i.e. ω = θ′ − L, so V is analytic and Littlewood's lemma counts its zeros left of the line; Conrey's
  Q(−∂/L)ζ is the higher-order version. Any other ω makes L depend on t and V non-analytic: the extra Prüfer freedom
  is unusable.
- (P4) (Hadamard product of (s − 1)²ζ′(s); constant fixed by N(T) − N₁(T) = (T/2π) log 2 + O(log T), Berndt)
  φ′(t) = ½ log 2 + Σ_{ρ′} (β′ − ½)/((β′ − ½)² + (t − γ′)²) + E(t), E(t) → 0 (heuristic derivation: by Hadamard,
  −Re ζ″/ζ′(½+it) = Σ_{ρ′}(β′ − ½)/|s − ρ′|² + smooth terms from the real zeros of ζ′, the double pole at 1 and the
  constant; −θ′ cancels their ½ log t growth and Berndt's count fixes the constant; E is not bounded explicitly here,
  numerically E ∈ [0.0085, 0.0116] on [30, 95], §3.4).
  RH ⟹ (Speiser) every Poisson term ≥ 0 ⟹ φ′ ≥ ½ log 2 + E(t): Z's Prüfer phase never stalls. Consequence: every gap
  between consecutive zeros is ≤ 2π/log 2 + o(1) ≈ 9.06, far weaker than Goldston–Gonek's (π + o(1))/log log T under RH.
  Routine; not a node advance.
- Gaps: π = ∫_gap φ′ (exact). A long gap needs the ζ′-zero Poisson mass near it to be thin; short gaps need a ζ′ zero
  close to the line (Soundararajan 1998, Zhang 2001, Radziwiłł 2014). Per-gap Sturm/Picone bounds are Wirtinger's class
  (explore53, gapjinf). No non-local inequality comes out: φ′ is a ratio, its mean values over gaps are not available.
- (P4) integrated over [t₀, T] (each Poisson term gives ±π up to edge terms O(log T)) is Levinson–Montgomery:
  N₀(T) = (T/2π) log 2 + N₁⁺(T) − N₁⁻(T) + O(log T), so with Berndt N − N₀ = 2N₁⁻ + O(log T) (N₁^± = zeros of ζ′
  right/left of the line). The Prüfer picture re-derives it in two lines; nothing beyond it.

## 2. Transfers (a)+(b): Birman–Schwinger and Weinstein–Aronszajn → N13, N1 (reductions, no new test needed)

- Exact rank-one structure (wmat.py normalisation, window length ℓ = ln c): Q^ev = B^ev + 2|c⟩⟨c|, Q^odd = B^odd − 2|s⟩⟨s|,
  c = cosh(x/2), s = sinh(x/2) on the window, B = Arch − Prime = E_c − S(c) (explore39: E_c a killed Lévy generator with
  prime jumps, S(c) = 2Σ_{n<c}Λ(n)/√n − ψ(¼) + ln π). The pole is ζ-specific (DH and L(s, χ₅) have none).
- Cauchy/Aronszajn interlacing: λ_k(Q^ev) ∈ [μ_k(B^ev), μ_{k+1}(B^ev)], λ_k(Q^odd) ∈ [μ_{k−1}(B^odd), μ_k(B^odd)].
  Weinstein–Aronszajn: the eigenvalues of Q^ev off σ(B^ev) are the zeros of w(z) = 1 + 2⟨c, (B^ev − z)⁻¹c⟩, so
  Q^ev ≥ 0 ⟺ N(B^ev) = 0, or N(B^ev) = 1 and w(0) ≤ 0 (c generic, 0 ∉ σ(B^ev)); Q^odd ≥ 0 ⟺ B^odd > 0 and
  2⟨s, (B^odd)⁻¹s⟩ ≤ 1 (the rank-one Birman–Schwinger principle).
- Under RH this pins 0 ≤ λ₀(Q^ev) ≤ μ₁(B^ev) ≤ λ₁(Q^ev) and 0 ≤ λ₀(Q^odd) ≤ μ₀(B^odd) ≤ λ₁(Q^odd): the second even and
  first odd decay rates of the Lévy process sit within the Weil ladder (≈ e^{−4πc}) of S(c) — explore39's pinning, now as
  interlacing.
  So a W–A certificate of CCM's step one needs lower bounds on μ₂(B^ev) ≤ λ₂(Q^ev) and on the odd W–A root, both
  ladder-small: no precision is saved by the rank-one splitting.
- Any Birman–Schwinger splitting Q = A − V (A ≥ 0, V ≥ 0): if Qψ = εψ then ‖V^{1/2}A^{−1}V^{1/2}‖ ≥ 1 − ε/⟨ψ, Aψ⟩, so
  the BS margin is ≤ min_k ε_k/⟨ψ_k, Aψ_k⟩ over the ≈ 2c ladder vectors (explore19 counts c per parity): an O(1)
  margin forces A itself to carry the ladder (A ψ_k ≈ ε_k), i.e. a certificate of A ≥ 0 at the same precision
  (explore17's Φ-ladder deflation and weildefl's cost comparison are the instances).
- BS at an off-line point a = γ − iβ (one moved pair, all else on the line) is the Christoffel function of the on-line
  zeros at a (explore38) and gives explore19's detection law c* ≈ 0.072m + 0.73; the rank-two W–A determinant adds no
  freedom. Not repeated here.

## 3. Tests (numbers; scripts in code/, outputs in logs/)

### 3.1 Prüfer (d) on ζ and DH — p1_prufer.py (mpmath siegelz/zeta derivatives, dps 15–20; every run ≤ 16 MB)
- (P1) φ + θ + arg ζ′(½+it) ≡ 3π/2 (mod 2π) to 12 digits at t = 20.5, 37.3, 61.1, 99.9, 143.2, 201.7, 280.4, 1000.3, 5000.7
  (logs/p1_ident.log).
- ζ on [18, 300], step 0.02 (14 101 points, 1146 s): sign changes of Z = 137 = net kπ-crossings of φ = N(300) − N(18)
  (mp.nzeros) = 137; ∫φ′/π = 136.58 (endpoint phases); φ′ > 0 at every point; min φ′ = 0.4445 at t = 18, the left end
  (θ′(18) = 0.5262), against the RH floor ½ log 2 = 0.3466 of (P4); in the bulk φ′ ≈ θ′ ± gap fluctuations.
  logs/p1_zeta_18_300.log.
- DH (Davenport–Heilbronn, f = 5^{−s}Σ a_r ζ(s, r/5), ϑ_DH = Im log Γ(¾ + it/2) + (t/2) log(5/π); logs/p1_dh_18_25.log,
  logs/p1_dh_80_92.log): on [18, 25] φ′ ≥ 0.968;
  on [80, 92] sign changes = kπ-crossings = 6, but φ′ < 0 at 29 grid points, min φ′ = −9.50 at t = 85.70, at DH's
  off-line zero 0.8085 + 85.6993i (explore23's list). p7_dhprime.py locates the cause: DH′ has the zero
  0.4053715 + 85.7051360i, 0.0946 LEFT of the line (Speiser fails; DH's zero pair 0.8085/0.1915 + 85.6993i), whose
  Poisson term at 85.70 is −10.54; with ½ log 2 and ≈ +0.69 from the other DH′ zeros this is the observed −9.50
  (logs/p7_dhprime.log; its other 'roots' at Re s ≈ 51–57 are underflow artefacts of f′ ≈ 5^{−s}, discarded).
  The Prüfer phase runs backwards exactly where Speiser fails: zero verification, nothing more.

### 3.2 Rank-one pole, Cauchy/Aronszajn interlacing, Weinstein–Aronszajn scalars (a)+(b) — p3b_interlace.py
Galerkin, CCM Fourier basis N = 30 (wmat.py of explore19/39), mpmath eigsy at 80 digits (Galerkin values are upper
bounds; c = 10 is not converged at N = 30, only the pattern is used). Lowest three eigenvalues, logs/p3b_interlace.log:

| c | sector | Q | B = Q − pole | W–A scalar at 0 |
|---|---|---|---|---|
| 2 | even | 1.34e−3, 0.661, 1.192 | −1.286, 0.583, 1.176 | w = 1 + 2⟨c,B⁻¹c⟩ = −1.2e−3 |
| 2 | odd | 0.0788, 0.931, 1.383 | 0.0904, 0.933, 1.384 | 1 − 2⟨s,B⁻¹s⟩ = 0.870 |
| 5 | even | 1.22e−17, 6.89e−12, 5.76e−7 | −3.334, 4.41e−12, 4.50e−7 | −8.6e−18 |
| 5 | odd | 1.30e−14, 2.04e−9, 6.26e−5 | 7.46e−10, 3.91e−5, 0.0656 | 4.5e−13 |
| 10 | even | 1.1e−40, 3.2e−34, 2.2e−28 | −5.119, 2.0e−34, 1.7e−28 | −7.3e−41 |
| 10 | odd | 2.6e−37, 2.8e−31, 1.3e−25 | 1.0e−31, 7.2e−26, 8.4e−21 | 7.8e−36 |

Reading: B^ev has exactly one negative eigenvalue (the archimedean–prime well at t ≈ 0, −1.29 … −5.12), the pole lifts
it; μ₁(B^ev) sits between λ₀(Q^ev) and λ₁(Q^ev), and the even W–A scalar is ≈ −0.7 λ₀(Q^ev): every quantity a
Weinstein–Aronszajn or Birman–Schwinger certificate of positivity or of CCM's step one must bound is ladder-small.
The splitting saves no precision (explore39's pinning, as interlacing).

### 3.3 Bound-state counting (c): Weyl / Berezin–Lieb / CLR-type counts are vacuous — p2_berezin.py
Pole-free part B = Toeplitz compression of the symbol σ_c(t) = Ω(t) − X_c(t), Ω = Re ψ(¼ + it/2) − ln π,
X_c(t) = 2Σ_{n<c} Λ(n)n^{−½} cos(t ln n), to Paley–Wiener functions of a window of length ℓ = ln c. Berezin–Lieb
(Tr φ(PAP) ≤ Tr Pφ(A)P, φ convex, φ(0) = 0, φ(x) = x₋) gives the rigorous bound Σ_j λ_j(B)₋ ≤ B(c) := (ℓ/2π)∫(σ_c)₋;
the Weyl count is W(c) := (ℓ/2π)|{σ_c < 0}|. Bohr–Jessen (t ln p equidistributed, Ω ≈ ln(t/2π)) gives
|{σ_c < 0}| ≈ 4π E[e^{X_c}], so W ≈ B ≈ 2ℓE[e^{X_c}], E[e^{X_c}] = Π_p (1/2π)∫exp(2 ln p Σ_k p^{−k/2} cos kθ)dθ ≈ exp(Σ_p ln²p/p)
= exp((½ + o(1)) ln²c). logs/p2_berezin.log (grid step 0.005 over |t| < 2πe^{max X_c + ½}):

| c | ℓ | E[e^{X_c}] | 2ℓE (BJ) | grid W | grid B | truth (Galerkin N = 30): Σλ₋(B), #neg B^ev / B^odd, #neg Q |
|---|---|---|---|---|---|---|
| 2 | 0.693 | 1 | 1.39 | 1.39 | 1.56 | 1.286, 1 / 0, 0 |
| 5 | 1.609 | 2.18 | 7.0 | 7.05 | 7.04 | 3.334, 1 / 0, 0 |
| 10 | 2.303 | 8.13 | 37.4 | 39.3 | 36.8 | 5.119, 1 / 0, 0 |
| 20 | 2.996 | 52.2 | 313 | — | — | 7.21, 1 / 0, 0 |
| 47 | 3.850 | 608.5 | 4686 | — | — | 10.52, 1 / 0, 0 (certified record) |
| 100 | 4.605 | 1.21e4 | 1.11e5 | — | — | — |

(#neg counts at c ≥ 10 from interlacing, since float Galerkin cannot resolve the ladder.) The one real negative
eigenvalue of B is the archimedean well |t| < t* = 6.29 (measure 12.58, the c = 2 grid value; σ_c(0) = −5.37 − X_c(0)); Berezin–Lieb is sharp there
(1.286 vs 1.56 at c = 2). Every other negative excursion of σ_c, where the prime sum beats ln(t/2π) (Bohr–Jessen large
values), binds nothing: positivity is a band-limited cancellation, not smallness. Any phase-space bound with any
admissible symbol σ̃ (same Fourier transform on (−ℓ, ℓ)) vanishes iff σ̃ ≥ 0 can be chosen, i.e. iff a nonnegative
Krein extension exists (Krein 1940), i.e. iff the form is positive: counting tools re-encode the problem. The record
method is exactly this pair: a lattice (periodized) admissible symbol, whose ≈ W(c) negative nodes are shown not to bind
by Birman–Schwinger (§3.5). Bargmann, Calogero and GGMT
need a local operator −d² + V on a half-line; Weil's form is nonlocal in x and a band-limited Toeplitz form in t.

### 3.4 (P4) pointwise — p4_zetaprime.py
Newton from a grid of starts found 30 zeros of ζ′ with 23.3 < γ′ < 133 (first 2.4631619 + 23.29832i; min β′ = 0.6356,
Speiser holds); the argument principle on [0.501, 5] × [1.5, 130] gives 29.0 = the 29 found below 130. With these 30
(and their conjugates) E(t) := φ′(t) − ½ log 2 − Σ_{ρ′} Poisson terms lies in [+0.0085, +0.0116] at t = 30, 35, …, 95,
smooth, rising at both ends (O(1/t) on the left, the ζ′ zeros above 133 on the right), while φ′ itself ranges over
[0.593, 3.059]. So φ′ = ½ log 2 + (Poisson-smoothed ζ′ zeros) to 1 % of its size: the zeros of Z are generated by the
zeros of ζ′, each ζ′ zero at distance δ′ from the line contributing a Poisson bump of mass π and width δ′.
logs/p4_zetaprime.log (284 s, 26 MB).

### 3.5 Side result: the record's Birman–Schwinger rank follows the Bohr–Jessen Weyl count — p5_lattice.py, p6_bsrank.py
Chekmachev's certificate (N1 record; Zenodo 23848597, audited ChekAudit, WeilZeta, ExploreFMM) reduces Q_c ⪰ ε to his
Thm 5.1, P − BBᵀ ⪰ 0 ⟺ λ_max(BᵀP⁻¹B) ≤ 1: the Birman–Schwinger principle for the finite-rank negative part of his
lattice symbol, of order |K−| + |J−| = the number of negative lattice nodes (his Table A.2, p. 34). Against the phase-space
count W(c) = 2 ln c · E[e^{X_c}] of §3.3 (no fit, same formula for every c), logs/p6_bsrank.log:

| c | 5 | 10 | 16 | 20 | 25 | 31 | 45 | 47 | 71 |
|---|---|---|---|---|---|---|---|---|---|
| BS order | 9 | 33 | 133 | 311 | 558 | 1663 | 4424 | 7170 | 38 377 |
| W(c) | 7.0 | 37.4 | 113 | 313 | 496 | 1071 | 4633 | 4686 | 24 895 |
| ratio | 1.29 | 0.88 | 1.18 | 0.99 | 1.13 | 1.55 | 0.95 | 1.53 | 1.54 |

For c = 10…71 (25 windows) the ratio stays in [0.869, 1.553] (geometric mean 1.16, log–log slope 1.036) over three
orders of magnitude; prime c sit high (1.29–1.55), composite c at 0.87–1.19, although c = 30, 31 and c = 45, 47 have the
same atoms and differ only in ℓ (lattice effect, not explained here). Counting negative values of σ_c itself on his two node systems (k < kscan) gives 7, 51, 126,
381, 548 at c = 5, 10, 16, 20, 25 (logs/p5_lattice.log), the same scale. Prediction: c = 100 needs a BS block of order
≈ 1.1·10⁵ (0.97–1.7·10⁵), c = 150 ≈ 9.6·10⁵, c = 200 ≈ 5.6·10⁶, i.e. exp((½ + o(1)) ln²c) — the Bohr–Jessen factor F_S
that aclgh2 found for ACLG's onset. Empirical, about one method's cost; the scan (weildefl's kscan law) dominates anyway.

## 4. Adversarial (iv): does any of it tell ζ from Davenport–Heilbronn or from a positive-definite model?

- Prüfer (P1–P3): identities for every real-on-the-line function with a functional equation (DH: same counts, 6 = 6 on
  [80, 92]); they separate nothing. (P4)'s floor separates ζ from DH (φ′_DH = −9.50 at 85.70) only because Speiser's
  criterion (no ζ′ zero left of the line) is RH itself, locally: it is zero verification, as explore23/explore41 found for
  ordinate statistics.
- Pole interlacing (a)+(b): the pole is ζ's (DH and L(s, χ₅) have none), but the conclusion "every BS/W–A margin is
  ladder-small" holds for DH below its failure c ∈ (30, 31) too (DH even λ₀ = 8.8e−8, 4.5e−18 at c = 10, 20,
  ExploreEntropy's table; 9.0e−8 at c = 10 in logs/p3_polefree_N30.log; explore19's depth slope −4π/5): structure of positive Weil-type forms, not RH information.
- Counting (c): vacuity holds for any positive form whose symbol has negative excursions (any Weil form with
  Euler-product coefficients); W(c) uses the Euler product (independent phases p^{−it}, weights Λ(n)n^{−½}) and so is
  ζ-specific in its numbers, but says nothing about positivity: it is a cost, not a criterion.

## 5. Novelty
- arXiv API: HTTP 429 on every query (web_search approach arxiv ×11, with 15 s retries, and .scratch/axq.py ×4,
  05:08–05:50 UTC); listing scan (.scratch/lst.py, math.NT/CA/CV/SP/FA/PR, Oct–Nov 2026) for Prüfer, Birman–Schwinger,
  Aronszajn, Berezin, Lieb–Thirring, Cwikel, Bargmann, Sturm, Levinson, Speiser, Hardy's Z, zeta′, Weil form, lattice
  minorant, Bohr–Jessen: only off-topic titles plus 2610.32704 (Weil-form spectral convergence, already in TREE) and
  2611.21551 (complex moments of ζ′, not this); Zenodo '"Prüfer" zeta', '"Birman-Schwinger" Weil': nothing relevant.
- Classics covering the transfers: Speiser (1934); Berndt, J. London Math. Soc. 1970 (N₁(T)); Levinson–Montgomery,
  Acta Math. 133 (1974) (zeros of ζ′ and of ζ left of the line, Σ(β′ − ½)); Levinson (1974); Soundararajan (1998),
  Zhang (2001), Radziwiłł (2014) (ζ′ zeros near the line ↔ small gaps); Berezin (1972), Lieb (1973) (Berezin–Lieb);
  Krein's extension theorem; Chekmachev's Thm 5.1 (the BS dual check, 2026). Nothing here is claimed new except the
  empirical size law of §3.5, which is a cost observation, not a discovery (gate: no new theorem, tool or barrier).

## 6. Verdict
- (d) Prüfer/Sturm → N8: no new bound. The Prüfer angle of Z is 3π/2 − θ − arg ζ′ (exact), its crossing count is the
  Levinson–Montgomery/Berndt count, its frequency freedom is Levinson's choice of V, per-gap Sturm comparison is
  Wirtinger's class; the RH floor φ′ ≥ ½ log 2 (P4) is routine and weak (gaps ≤ 2π/log 2 ≈ 9.06). Large gaps would need
  mean values of the ratio φ′ over gaps, which are not available: no non-local inequality comes out.
- (a)+(b) Birman–Schwinger/Weinstein–Aronszajn → N13, N1: re-encodings (explore39's pinning, explore17/weildefl
  deflation, explore38/explore19 detection); interlacing puts every margin into the ladder (§3.2); BS is already the
  record's dual check.
- (c) counting → N1: Weyl/Berezin–Lieb counts 2 ln c·E[e^{X_c}] = exp((½ + o(1)) ln²c) against one true negative
  direction of the pole-free form (the archimedean well, lifted by the pole) and none of Q: vacuous; the record method is
  a lattice symbol plus a Birman–Schwinger proof that its ≈ W(c) negative nodes do not bind (§3.5).
- Side law (empirical, cost only): the record certificate's BS order tracks W(c) within ×0.87–1.55 for c = 10…71;
  c = 100 needs a block of order ≈ 1.1·10⁵ (within reach of ExploreFMM's O(n)-memory structured LDLᵀ; the scan,
  kscan ≈ 1.2·10¹⁶ at c = 100 by weildefl's law, stays the bottleneck).
- No bound, rate or barrier on a node quantity: TREE N8 (67.3639 %, large gaps Λ > 3.34), N1 (c ≤ 47, even c ≤ 71) and
  N13 (step one on λ² ≤ 227/6 ∪ {40}) are unchanged; no method ceiling claimed (nothing to certify). Unconditional except
  the floor (P4), which assumes RH.

## Suggested NOGO line
- explore57 (Orca, exploration: spectral theory of Schrödinger operators — Birman–Schwinger, Weinstein–Aronszajn/Aronszajn–Donoghue, Weyl/Berezin–Lieb/CLR counting, Sturm/Prüfer → N8, N13, N1) NO IDEA: the Prüfer angle of Z against the θ-oscillator is 3π/2 − θ − arg ζ′(½+it) (12 digits), its kπ-crossings are the zeros (137 = 137 on [18, 300]) and its frequency is Levinson's V: the Levinson–Montgomery count; φ′ = ½ log 2 + Poisson terms of the ζ′ zeros (residual 0.009–0.012 on [30, 95]), so RH gives only φ′ ≥ ½ log 2 (gaps ≤ 9.06), and DH's phase runs backwards (φ′ = −9.50 at 85.70, a DH′ zero 0.0946 left of the line): zero verification; per-gap Sturm is Wirtinger (explore53, gapjinf); the rank-one pole interlaces every BS/W–A margin into the ladder (w(0) ≈ −0.7λ₀ at c = 5, 10; explore39's pinning) and BS is already the record's dual check (Chekmachev Thm 5.1); Weyl/Berezin–Lieb counts 2 ln c·E e^{X_c} = exp((½ + o(1)) ln²c) (37, 313, 4686 at c = 10, 20, 47) face one true negative direction (the archimedean well): vacuous; the record is exactly a lattice symbol plus a BS proof that its negative nodes do not bind, and its BS order (his Table A.2) tracks that count within ×0.87–1.55 for c = 10…71 (38,377 vs 24,895 at c = 71), so c = 100 needs ≈ 1.1·10⁵ [ExploreBS](Riemann/Deep/ExploreBS/REPORT.md)

## Files
code/p1_prufer.py (Prüfer scans ζ, DH), p2_berezin.py (Weyl/Berezin–Lieb/Bohr–Jessen), p3_polefree.py and
p3b_interlace.py (Galerkin, pole interlacing, W–A scalars; reuse ExploreFunctional/code/wmat.py), p4_zetaprime.py (ζ′
zeros, (P4) residual), p5_lattice.py (negative lattice nodes of σ_c), p6_bsrank.py (BS order vs W(c)), p7_dhprime.py
(DH′ zero at 85.7); logs/ holds every output (memcap peaks 0–126 MB; one process at a time).

RH relevance: N8 (Prüfer = Levinson's method, no new bound); N1 and N13 (no new bound; a certificate-size law for the N1
record method).
