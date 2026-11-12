IDEA: point-process realizability (a hard-core set of density 1 must have many nearest- and next-nearest-neighbour pairs) augments Lamzouri's small-gap criterion past its certified 0.775 floor: certified D(λ) > 0 for every λ ≥ 0.768, unconditionally (was 0.7851) (node N8)

# explore40: point-process statistics → N8 (Orca exploration, 12 Nov 2026, 09:54–11:23 UTC)

Derivation: [DERIVATION.md](DERIVATION.md). Code: [code/](code/). Certificates: [certs/](certs/). Logs: [logs/](logs/).

## 1. Transfer (bridges measured before choosing)
| bridge | tool | measurement | verdict |
|---|---|---|---|
| N8: ∫_1^b F, multiplicities (RH) | Yamada Var N_L ≥ θ(1−θ) | the [−1, 1] part of F alone gives Var ≥ A(L) = Cin(2πL)/π² > θ(1−θ) for every L > 1 | vacuous (§3a) |
| N8: small gaps, unconditional (Lamzouri Thm Main2, record 0.7851, certified floor 0.775) | realizability of the pair measure of a hard-core set (triangle/hypermetric, exact 1D local realizability, consecutive-gap ground state) | explore34's dual pseudo-lattices violate μ(2λ) ≥ 2μ(λ) − 1 by 0.15 | chosen |
| N1 (Weil as structure-factor positivity), N13 | – | Weil positivity is a test-function statement on f^·conj(f^), no point set to realize | not chosen |
Tool: in physics, a candidate pair correlation is realizable only if it is the pair statistic of an actual point process
(Yamada; Costin–Lebowitz; Kuna–Lebowitz–Speer); for a 1D hard-core gas the pair statistic is built from consecutive gaps (Tonks).
Transfer: Lamzouri's contradiction argument runs on Z_0, the zeros without a short neighbouring gap: an actual set with gaps > λ and
density 1. His step (iii) bounds its off-diagonal pair sum by 0 (f^ ≤ 0 beyond λ). Realizability gives a sharper, negative bound.

## 2. Statements (proofs: DERIVATION.md)
- Theorem 1 (§4). Let 1/2 < θ < 4/7, λ < U ≤ 3λ, f even, f ≥ 0, Schwartz, f^ ≤ 0 on [U, ∞), and A bounded, β ≥ 0, c with
  (H) f^(x) + f^(x + y) + A(x) − A(y) + β(1 − x) ≤ c for all x, y ≥ λ. If J_θ(f) > f^(0) + 2c, then D(λ) > 0.
  A = β = c = 0, U = λ recovers Thm Main2. Proof: Lamzouri's steps (i)–(ii) unchanged; (iii) replaced by summing (H) over consecutive
  gap pairs of Z_0 (pairs three or more apart are ≥ 3λ ≥ U apart, where f^ ≤ 0; A telescopes; the β-term is the mean-gap-1 constraint).
- Corollary (§5). D(λ) > 0 for every λ ≥ 0.768, unconditionally: certificate with f = degree-42 polynomial × Gaussian (22 rational
  Hermite coefficients), U = 3λ, A on 78 nodes; margin J − f^(0) − 2c = +0.0014626 (λ = 0.768), +0.0051 (0.77), +0.0137 (0.775).
- Barrier crossed: explore34 proved Thm Main2 cannot give any λ ≤ 0.775; Theorem 1 gives 0.775 with margin 0.0137 and 0.768.
- Proposition Y (§1). Under RH, every single-window Yamada inequality is implied by Montgomery's theorem (A(L) > θ(1−θ) for L > 1).
- Proposition R (§2). explore34's certified duals (λ = 0.76–0.775, θ = 4/7) are not pair measures of any simple point process.

## 3. Tests
a) Smallest case, Yamada (logs/yamada_margins.log): θ(1−θ) − A(L) = −0.247, −0.067, −0.031, −0.035, −0.098, −0.316 at L = 1, 1.25,
   1.4, 1.5, 1.7, 2; closest approach L ≈ 1.4. GUE Var(1.5) = 0.388 vs 0.25: zeros are far from integer-rigid at the pair level.
b) Smallest realizability case, the triangle μ(2λ) ≥ 2μ(λ) − 1 on explore34's duals (logs/triangle_explore34_duals.log):
   λ = 0.775: μ(λ) = 0.757, μ(2λ) = 0.367, violation +0.147; 0.77: +0.156; 0.76: +0.165; the θ = 1 (RH) duals satisfy it (−0.20, −0.22).
c) Float dual floors with realizability rows (dual_real.py, dual_pattern.py; θ = 4/7, crossings of t* = 0):
| support X | plain | + triangle, packing | exact range-R realizability, R = 20 (range 2.14) | R = 30 (range 3.2) |
|---|---|---|---|---|
| 10 | 0.7632 | 0.7617 | 0.7513 | 0.7437 |
| 18 / 26 | 0.7722 / 0.7758 | – | – | 0.7524 / 0.7562 |
   Lattice step (range 2.2, X = 10): 0.7513 / 0.7517 / 0.7516 at k0 = 7 / 10 / 14. Fit R = 30: 0.7640 − 0.2031/X vs plain 0.7839 − 0.207/X.
   The realizable optimum at λ = 0.75 is a run process: gaps λ (70.7%), 2λ, 1.61, 1.71; two consecutive λ-gaps 42.8%; long gaps isolated.
d) Continuum primal (primal_2gap.py, float): Theorem 1's LP threshold ≈ 0.7678 (D = 20–40, converged; margins +0.0134 at 0.775,
   +0.030 at 0.7851), the same code's plain LP ≈ 0.786 (D = 20); with U < 2λ (nearest neighbours only) no gain. Under RH (θ = 1)
   no gain either: 0.6042 vs plain 0.6046 (D = 20, logs/primal_theta1.log; β = 0), as explore34's θ = 1 duals are triangle-feasible.
e) Certificates (cert_aug.py, exact Sturm + Arb; xcheck_aug.py, mpmath quadrature + direct Fourier transform + 10⁶ random points):
   all PASS, table in DERIVATION §5. Rerun: `cd code && EPS=1/1000 uv run --no-project --with numpy --with python-flint python3
   cert_aug.py ../certs/cert768.json 5` (6 s, 41 MB; c_cert = −0.3478616532); `EPS=1/1000 uv run --no-project --with numpy --with mpmath
   python3 xcheck_aug.py ../certs/cert768.json -0.3478616532` (45 s, 23 MB).
f) At scale, real zeros (22,462 zeros, 100 < t < 2·10⁴, zeros_check.py): D(0.768) = 0.302; with the λ = 0.768 certificate's f the full
   pair sum (1/N)Σ f^ = 0.366 ≥ J = 0.307 (Montgomery-type input holds), above the hard-core ceiling f^(0) + 2c = 0.305: the zeros
   escape the contradiction exactly by having small gaps. Their Z_0 has density 0.439 (energy −0.242 per point vs c = −0.348).

## 4. Adversarial
- Realizability holds for every simple point set; Theorem 1 is a statement about test functions plus Lamzouri's ζ-specific input
  (Conrey's density). Its conclusion is true for the sine process (D(0.768) ≈ 0.32, explore34's HKPV table), Poisson (1 − e^{−0.768} =
  0.536), ζ data (0.302) and Davenport–Heilbronn's on-line zeros (≈ 0.21; explore34: 0.199 at 0.76, 0.219 at 0.7851). It sees
  ordinates only (an off-line pair keeps its ordinates). It does not separate ζ from DH, GUE or Poisson: no RH information.
- What it does separate: the AH-like pseudo-lattices that cap pair-only LPs (explore34) from genuine point processes. Genuine processes
  still cap it: the 7/11-discrete sine DPP has every pair input and no gap below 7/11, so no realizability argument goes below 0.636.

## 5. Novelty
- arXiv math.NT listings 2026-10 and 2026-11 (titles with gaps/consecutive): 2611.14486 (Lamzouri, 0.7892), 2610.22927 (consecutive
  moderate gaps, cluster counting), 2610.26038 (Bondarenko–Heap, RH), 2610.28527 (Durkan–Pearce-Crump), 2610.31744: none uses
  realizability or improves 0.7851. arXiv API (3 queries) and GXL (2 claims: no support; nearest: Aryan 1910.02408, consecutive gaps
  under the AH via Dirichlet polynomials; BGMM 2208.02359) find no hard-core/consecutive-gap bound in a pair-correlation gap argument.
  Lagarias–Rodgers (1905.12123, 1907.03391) build lattice processes mimicking band-limited data (the dual side), not primal bounds.
  Our files: SmallGap/NOVELTY.md (0.7851 not pre-empted on 11 Nov), ExploreDPP (floor 0.775 for Thm Main2 only).
- New as far as checked: Theorem 1, the certificates, Propositions Y and R. Gate: a record on the N8 line (RECORDS) and a tool that
  crosses a certified method barrier (DISCOVERIES candidate); still needs an independent certificate rerun and a full review.

## 6. Reviewer (verdict-only, claims inline, 11:18 UTC)
(1) step (iii′) incl. the mean-gap β-term and dropping pairs ≥ 3 apart: CORRECT, no sign condition on c needed. (2) Lemma T′: GAP as
stated, fixed: h_R even, 0 ≤ h_R ≤ 1, = 1 on [−θ, θ], reserve G compact C² (DERIVATION §4); then F^(0) ≤ f^(0) + 2c_r. (3) statement:
CORRECT with Thm Main2's hypotheses kept; normalise by dividing f, A, β, c together (added). (4) checker: bilinear constant 3M2δ²/8
CORRECT with M2 on [λ, 2U]; GAP: R2, R3 need A's constant tail A(x) = A(U) for x ≥ U, which the certificates use (now stated).
(5) margin not independently re-established by the reviewer (it did not run code): the checker reruns above are the evidence.
Round 2 (11:40 UTC, both closed, mathematics unchanged): (a) Lemma T′ now says why F = f h_R + c_r G, supported far beyond θ, is
admissible: Lamzouri's proof (TeX l.418–440) applies Prop SmallSetZeros only to g = ψ_ε F (ψ_ε = 1 on [−θ, θ], support in
[−θ−ε, θ+ε]) and reaches F ≥ g ≥ 0 through the positive identity Σ F^w = (2/π)∫∫ F(α)|Σ T^{α(ρ−1/2)}/(1 − (ρ − (1/2 + it))²)|²
(DERIVATION §4). (b) c_cert is now printed to 10 digits: −0.3478616532 = R1 = corners −0.3478888283 + 2.718·10⁻⁵ (the old −0.347862 was
its 6-decimal rounding, which read as below R1); the three certificates were rerun after the print-only change, margins identical, logs kept.

## 7. Caveats
- Rests on Lamzouri 2611.14486v1 (unrefereed; audited PASSED by lamzaudit) steps (i)–(ii), as SmallGap does; (iii′) and Lemma T′ are
  ours: one verdict-only reviewer (§6), no referee.
- Checker: exact Sturm (python-flint), Arb values; corner maxima in float64 from Arb upper endpoints (+10⁻¹² slack) and again entirely
  in Arb (FULLARB=1, logs/cert_check_*_fullarb.log): identical to 10⁻⁹; M2 via a Cramér-type third-derivative bound. Margin at 0.768: 1.5·10⁻³.
- B_0 is shifted by 1/1000 in the checker (EPS) to remove two near-double roots of the LP optimum; the certificate is the shifted f.
- Float floors (§3c, d) are not barriers. The 3-gap version (U ≤ 4λ) was not run (memory); the lattice suggests ≈ 0.764.

## 8. Process
One process at a time, single-threaded; LP self-abort at 285 MB (ru_maxrss) plus memcap 290. Peaks: lattice LPs 92–213 MB, continuum
LPs 96–239 MB, checkers 21–47 MB. One LP (dg = 1/200) self-aborted at 285 MB; one (k0 = 20) hit the 900 s alarm. Downloads deleted.

## Suggested NOGO line
- explore40 (Orca, exploration: point-process statistics — Yamada/realizability, hyperuniformity, Palm/hard-core → N8) IDEA: single-window Yamada is vacuous under Montgomery (A(L) = Cin(2πL)/π² > θ(1−θ) for all L > 1: nothing for ∫_1^b F or multiplicities), but explore34's small-gap duals are not point processes (μ(2λ) ≥ 2μ(λ) − 1 fails by 0.15); summing a consecutive-gap inequality over Z_0 (gaps > λ, density 1) augments Lamzouri's Thm Main2 past its certified 0.775 floor: certified D(λ) > 0 for λ ≥ 0.768 (was 0.7851; float floors 0.7678 two-gap, ≈ 0.764 lattice range 3.2; 7/11 the DPP limit; inert at θ = 1); RH-blind [ExplorePointProc](Riemann/Deep/ExplorePointProc/REPORT.md)

RH relevance: N8 (unconditional small gaps: record λ ≥ 0.768 and a tool that crosses the pair-only LP floor); no RH information.

## Files
code/: margins.py (Yamada vs A(L)), yam_dual34.py and tri_check.py (Yamada and triangle on explore34 duals), dual_real.py (triangle/packing
rows), dual_pattern.py (exact range-R realizability), pstruct.py (realizable optimum structure), primal_2gap.py and verify_2gap.py
(continuum LP), lp_cert.py (LP with margins → rational certificate), cert_aug.py (checker), xcheck_aug.py (independent check), proots.py,
zeros_check.py. certs/: cert768.json, cert770.json, cert775.json. logs/: cert_check_*.log, xcheck_768.log, yamada_margins.log,
triangle_explore34_duals.log, scan_X10.log, pat_*.log, primal_scan*.log, struct_R30_X18_lam0.75.json. Lattice LP scripts write logs/ in
the working directory.
