NO IDEA: antenna synthesis maps one-to-one onto Weil window theory, but every transfer is already known (Berry 2012; Connes–Consani 2021; Connes–Consani–Moscovici 2025; Elston, Zenodo 23174403, Oct 2026), depends only on the gamma factor and conductor (Davenport–Heilbronn obeys it until it fails at c ≈ 30.7), and sees off-line zeros only below height ≈ 4.4πc, so it gives N1 a sharper barrier and nothing toward RH.

# ExploreAntenna: antenna arrays and beam synthesis → RH tree (explore19, Orca, 6 Nov 2026)

Window and form as in TREE.md N1: W = [−L/2, L/2], L = ln c, Q(f) = 2⟨f,cosh(x/2)⟩² − 2⟨f,sinh(x/2)⟩² + (1/2π)∫(Re ψ(1/4+it/2) − log π)|F(t)|²dt
− 2Σ_{n<c} Λ(n)n^{−1/2}C_f(log n), F(t) = ∫f e^{itx}dx. By the explicit formula Q(f) = Σ_ρ F(γ_ρ)·conj F(conj γ_ρ), i.e. Σ_γ |F(γ)|² under RH.

## Dictionary (what the field gives)
| antenna / filter design | Weil window |
|---|---|
| aperture of length L, far-field pattern | window W, F(t) = FT of f (entire, type L/2) |
| receivers / prescribed nulls | zeros γ: Q(f) = Σ_γ |F(γ)|² (RH) |
| Woodward–Lawson sampling basis (shifted sincs) | CCM's Fourier basis E_N = span{cos(2πjx/L)} (our Galerkin basis) |
| visible region, Nyquist density L/2π | zeros sparser than L/2π below T* = 2πc, denser above |
| superdirective (supergain) modes, Landau/Chu count | near-null ladder; Landau excess LT*/π − 2N(T*) = 2c − 7/4, i.e. c − 7/8 per parity |
| Chu–Harrington exponential Q, Slepian leakage e^{−2c_s} | depth λ_min ≈ e^{−4πc}·c^5 (prolate leakage at c_s = 2πc) |
| tolerance / sensitivity factor | λ_min moves by ~e^{4πc}× the relative error of one prime weight (WeilDefl: w_2·(1 ± 10⁻⁷⁰) flips c = 16) |
| Schelkunoff / Taylor n̄ null placement | the vector in span{φ_0..φ_{n̄−1}} with F(γ_1..γ_{n̄−1}) = 0 |
| Fejér–Riesz, Carathéodory–Toeplitz, FIR via autocorrelation (Wu–Boyd–Vandenberghe) | Krein: Q ≥ 0 on W ⟺ Weil's distribution on (−L, L) extends to a positive measure; Chekmachev's periodized lattice minorant is a lattice version |
| non-uniform array | prime array: elements at ±log n, weights Λ(n)/√n, n < c |

## Survey (arXiv first; scouts s1–s4, logs in .scratch/explore19)
- Berry, "Riemann zeros in radiation patterns", J. Phys. A 45 (2012) 302001, and II, J. Phys. A 48 (2015) 385203 (no arXiv id): ζ's zeros are the nulls of a
  radiation pattern; truncating the source at |u| ≤ log N keeps the zeros up to t_N (aperture O(log T)). Visible-region picture: known.
- Connes–Consani, arXiv 2106.01715 (TeX l. 181–195, 644–669): the tiny eigenvalues of Weil's form on [S⁻¹, S] come from prolate functions pushed through
  the map E; about μ = S² of them in the even sector, one fewer odd. Superdirective count: known. Connes–Moscovici 2112.05500, CCM 2310.18423: same mechanism.
- Connes–Consani–Moscovici, arXiv 2511.22755 (l. 1324–1326, 1394–1400): prolate defect 1 − χ(λ) ~ C e^{−4πλ² + 9 log λ}; ε_λ tracks it (figure, slope ≈ −12).
  Superdirective depth: known. Chekmachev, Zenodo 23978194: scale Θ(c) ~ c⁵e^{−4πc} and the log-Laplacian edge balance of the minimiser
  (cf. superdirective currents concentrating at the array edge, arXiv 2606.08174). Zenodo 23983773: conditional lower bound e^{−(4π+ε)c}. Zhu 2608.24827, Conj. 12.1.
- Elston, Zenodo 23174403 (v1.1.0, 6 Oct 2026, unrefereed) and github.com/zacharyelston/weil-windows (docs read 6 Nov): the decay law with slope −4π in
  c/q for Dirichlet L (DECAY_LAW.md, P-DL1); the window as a sampling problem with T* = 2πc, Landau 1967, Beurling, Ortega-Cerdà–Seip, Kaiser–Bessel,
  Hermite and E-map prolate trial functions (prolate within 1.05–2.5 of the minimum; RESEARCH_R2_SAMPLING.md, R2_THEOREMS.md); DH's first failure
  c ∈ (30.74404, 30.74497] with certified witnesses (CONTROL_CERTIFICATES.md); the quadruple term −4δ²|F′(γ₀)|², the forward criterion σ₁ = 1 and
  exclusion heights all ≥ 6 orders below 3·10¹² (RESEARCH_R1_EXCLUSION.md). The sampling transfer is essentially this work.
- Time–frequency and superdirectivity: Slepian–Fuchs 1 − λ_n(c) ~ 4√π 8ⁿc^{n+½}e^{−2c}/n!; Bonami–Karoui 1012.3881, 1509.02646; Bonami–Jaming–Karoui
  1804.01257 (explicit bounds, Landau–Widom count); superoscillation cost (Ferreira–Kempf 2003, id per scout quant-ph/0305148, unverified); superdirectivity
  reviews 2602.04121, 2606.08174, Chu-limit modulation 2607.23688. Filter design: Wu–Boyd–Vandenberghe 1996 (not on arXiv); Fejér–Riesz generalisations 2005.11920.
- No arXiv/Zenodo item links antenna synthesis or superdirectivity to Weil positivity by name (arXiv math.NT listings 2026-10, 2026-11 show=2000 scanned
  for weil|prolate|superoscill|antenna|radiation|slepian|paley|sampling|explicit formula|smallest eigenvalue|connes: only 2610.32704, known); but the
  mathematics is the CC/CCM/Elston programme. Our files: no hit for antenna, Dolph, superdirect, Schelkunoff, Fejér–Riesz (brief).

## Transfers, tests, verdicts
Code: Galerkin matrices from closed forms (code/wmat.py, generalising ../WeilSimple/weil_matrix.py to gamma factor (q/π)^{s/2}Γ((s+a)/2):
z0 = 1/4 + a/2, b_k = 2k + 2z0, log(q/π), no pole for DH), eigenvalues by mpmath at 4πc/ln 2 + 96 bits or more. Checks: c = 5 reproduces
WeilSimple's U(5) = 1.2218e-17; DH entries equal an independent x-space quadrature to 1e-41 (logs/check_entries.log).

T1. Superdirective count and depth (Landau/Chu ↔ near-null ladder). Prediction: c − 7/8 near-null directions per parity; depth = prolate leakage at
c_s = 2πc. Test: eigenvalues below 1e-2 (code/count.py, logs/count.log): 4/3, 9/8, 15/14 (even/odd) at c = 5, 10, 16 against c − 7/8 =
4.1, 9.1, 15.1, i.e. ⌊c − 7/8⌋ even and one fewer odd (WeilDefl: c − 1 and c − 2 at c = 16, 20, 25). WeilSimple's U(c), c = 5..20, free fit
ln λ_min = −12.599c + 5.37 ln c + 15.4 (4π = 12.566); with 4π fixed, ln λ_min + 4πc = 5.00 ln c + 15.86 (residual ≤ 0.24; logs/fits.log). It predicts
5.6e-242 at c = 47 (Chekmachev's ε₄₇ = 1.87e-242), 4.6e-372 at c = 71 (his margin 2.4e-373 at ε₇₁) and 6e-530 at c = 100.06 (Wiktorowicz's core
bound 4.58e-530, ../WiktAudit). So a certificate at c needs at least 4πc/ln 2 ≈ 18.1c bits (c = 100: 1810; they used 7300).
Verdict: works, known (CC count; CCM prolate defect; Chekmachev Θ(c); Elston P-DL1).

T2. Visible region and null steering of off-line zeros (Berry + sampling). Test (code/movedpair.py): the xi-like function whose zeros γ_k, γ_{k+1} are
merged into the pair m ± iδ off the line; its form is Q_mod = Q + 4Re F(m+iδ)² − 2F(γ_k)² − 2F(γ_{k+1})² (even or odd f, ± symmetric). First window
c* where the Galerkin λ_min(Q_mod) < 0 (bisection; λ_min is non-increasing in c):
| k | m | δ = 0.3 | δ = 0.03 | δ = 0.003 | m/2π |
|---|---|---|---|---|---|
| 1 | 17.578 | 2.06 ± 0.04 | 2.51 ± 0.04 | 2.95 ± 0.05 | 2.80 |
| 5 | 35.261 | 3.28 ± 0.05 | 3.79 ± 0.06 | 4.36 ± 0.07 | 5.61 |
| 10 | 51.372 | 4.32 ± 0.06 | 5.10 ± 0.07 | | 8.18 |
| 20 | 78.241 | 6.42 ± 0.08 | 7.28 ± 0.09 | | 12.45 |
| 30 | 102.522 | 8.12 ± 0.06 | | | 16.32 |
| 50 | 144.556 | 11.19 ± 0.06 | | | 23.01 |
Fit: c* = 0.0722m + 0.73 (δ = 0.3, residual ≤ 0.12), c* = 0.0790m + 1.07 (δ = 0.03): a window c sees a count-preserving off-line pair only below height
m* ≈ 13.8(c − 0.73) ≈ 4.4πc, about 2.2 T*. Basis check: 1.6× the modes at k = 10 leaves (4.31, 4.375). Mechanism (as in Elston R1): a near-null f with a real
null at m between the pair gets −4δ²|F′(m)|² (the DH witness below has its null at 85.725). Prediction for Davenport–Heilbronn (pair 0.808517 + 85.699348i, conductor 5, slope ×5): c* = 31.7;
for its second pair (0.65083 + 114.16334i): 42.0.
Verdict: the reach law is ours and holds over m = 17.6..144.6, but it is the forward form of Elston's σ₁ = 1 criterion; as a statement about RH it is a
barrier: Chekmachev's c = 47 sees off-line pairs only below ≈ 640, c = 100.06 below ≈ 1370, and matching RH to 3·10¹² needs c ≈ 2·10¹¹.

T3. Taylor n̄ / Schelkunoff null placement versus the Woodward–Lawson Galerkin minimum (code/taylor.py; the brief's basis comparison). The best n̄
gives RQ 4.40e-15 at c = 5 (n̄ = 20) against λ_min 1.22e-17 (×360), and 2.12e-37 at c = 10 (n̄ = 60) against 1.96e-43 (×1.1·10⁶): the loss grows by
about e^{1.6} per unit of c.
Taylor keeps the sinc pedestal's 1/t sidelobes beyond the steered nulls; the good windows are smooth tapers (E-map prolates within 1.05–2.5, Elston;
Kaiser–Bessel trials reach RQ/minimum ≈ 1 at c/q = 2 but 9·10⁴ further out). Verdict: Dolph–Chebyshev/Taylor-type exact null steering is no certificate basis; the antenna-optimal window is the prolate, known.

T4. Fejér–Riesz / Krein dual certificates. Positivity on W ⟺ Weil's distribution on (−L, L) is the Fourier transform of a positive measure (Krein);
an explicit dual must reproduce the atoms at ±log n, which a finite atomic measure cannot, so it needs the zero measure itself (circular), and the
lattice (periodized) version is Chekmachev's minorant. Verdict: known (Krein/de Branges; R034, R040).

## Adversarial tests
(a) Davenport–Heilbronn (code/dh.py, witness.py; no Euler product, Γ((s+1)/2), conductor 5). Galerkin λ_min (even): 4.42e-14 (c = 16), 2.80e-22 (24),
7.80e-27 (28), 6.60e-28 (29, N = 118), 9.57e-29 (30, N = 123), −6.08e-30 (31), −1.92e-29 (32); a second negative eigenvalue from c ∈ (42, 44) (even) and
(44, 46) (odd). Depth slope −2.36..−2.62 per unit c against −4π/5 = −2.51. Rigorous witness at c = 31: Q(f) = [−6.077703566e-30 ± 3.5e-40]‖f‖² for an
explicit 85-term cosine polynomial f (Arb, closed-form entries), its pattern null at 85.725 between the pair. Control L(s, χ₅), χ(2) = i, same gamma
factor and conductor: positive to c = 48 (λ_min 5.3e-5, 2.9e-8, 1.1e-9, 5.5e-11 at c = 16, 32, 40, 48). So every law of T1–T3 holds for DH as for ζ, with
c → c/5; positivity at c < 30.7 is DH-blind; ζ's certificates at c ∈ [31, 47] are past DH's failure. Matches Elston's (30.74404, 30.74497]: not new.
(b) xi with a pair off the line: invisible to every window c < c*(m, δ) (table above). Positivity on any finite window is compatible with off-line
zeros above ≈ 4.4πc; the transfer cannot say anything about RH beyond zero verification to that height.
(c) Restatement? Yes: visible region = explicit formula + prolate concentration (CC, CCM, Elston R2); count/depth = Landau–Widom + Slepian; duals = Krein.

## Suggested NOGO line
- explore19 (Orca, exploration: antenna arrays, superdirectivity, Dolph–Chebyshev/Taylor/Schelkunoff synthesis, Fejér–Riesz/FIR spectral factorization → N1, N8, N13) NO IDEA: the window is an aperture whose pattern F is sampled at the zeros; superdirectivity gives c per parity near-null modes (Connes–Consani) of depth ≈ e^{15.9}c⁵e^{−4πc} (CCM's prolate defect; fitted on c = 5..20, matches the certificates at c = 47, 71, 100, so ≥ 18.1c bits per certificate), Krein/Fejér–Riesz duals are Chekmachev's lattice minorant, Taylor/Schelkunoff null steering loses ×360 (c = 5) to ×1.1·10⁶ (c = 10), and all of it depends only on gamma factor and conductor (DH: depth slope −4π/5, first failure c ∈ (30, 31), as Elston 23174403 found); a merged off-line pair at height m is first seen at c* ≈ 0.072m + 0.73 (δ = 0.3; DH predicted 31.7), so a window sees off-line zeros only below ≈ 4.4πc (c = 47: ≈ 640): N1 records are zero verification to height ~10³ [ExploreAntenna](Riemann/Deep/ExploreAntenna/REPORT.md)

## Caveats
- Thresholds c* and DH's positive side are Galerkin (upper bounds of λ_min; true first failures can only be earlier); the negative side is a real test
  function. The DH witness rests on the standard explicit formula for DH and on the closed forms (checked to 1e-41 against quadrature); it is not reviewed.
- The reach law is empirical (6 heights, 3 values of δ); its constant 13.8 ≈ 4.4π has no derivation here.

## Files (Riemann/Deep/ExploreAntenna/; run from this directory with uv run --no-project --with python-flint --with mpmath --with numpy python3)
- code/wmat.py (forms: zeta, DH, L(s, χ₅) control); code/movedpair.py k δ cmin cmax steps (T2; NMUL env for basis size); code/dh.py [kind=dh|chi5] c...;
  code/witness.py c (Arb witness, pattern nulls); code/check_entries.py (DH entries vs quadrature); code/taylor.py c (T3); code/count.py c... (T1 counts);
  code/fits.py (tables, fits);
  code/batch1.sh, batch2.sh (the runs). logs/: one log per run (mp_k*_d*.log, dh_*.log, chi5_ctrl.log, witness_c31.log, taylor_c*.log, fits.log).
- Scratch: .scratch/explore19/ (STATE.md, working copies).

RH relevance: N1 (barrier: window c sees off-line zeros only below ≈ 4.4πc); no new idea for N8 or N13.
