NO IDEA: gravitational-wave and black-hole mathematics transfers only known structures: Chandrasekhar's isospectral map is the Maass weight shift (it touches only the Gamma factor), quasinormal-mode spectral instability is the Davenport–Heilbronn instability of zeros inside a functional-equation class, a mirrored cavity is Lagarias–Suzuki (it needs only 0 < Re ρ < 1), and chirp matched filtering is Odlyzko–Schönhage/Hiary/Turing; each is blind to Davenport–Heilbronn or restates a known result.

# explore25: gravitational-wave astronomy and black-hole perturbation theory → RH tree
Orca agent explore25, 11 Nov 2026, 12:26–13:20 UTC. Every number below comes from [code/](code/) (Arb through python-flint, mpmath). Outputs are in [outputs.txt](outputs.txt) and derivations in [derivations.md](derivations.md). Each run took under 15 min and under 100 MB, single-threaded.

## Short answer
- **Chandrasekhar–Darboux (N13).** The Regge–Wheeler ↔ Zerilli map is, on the modular surface, the Maass weight shift: φ_{k+2}(s) = φ_k(s)(1 − s + k/2)/(s + k/2). The factor is unimodular on Re s = 1/2, and its "algebraically special" point s = 1 + k/2 is a holomorphic form for k ≥ 2. The resonances ρ/2 are shared by every weight. A finite Darboux chain carries ζ's scattering matrix to that of a ξ with one pair moved off the line, and Davenport–Heilbronn (DH) carries the same structure. **Dead (known).**
- **Spectral instability / pseudospectrum (N13, N8).** The admissible perturbations of an L-function are its functional-equation (FE) class. The test family is F_c = cos c·L(s, χ₁₃) + sin c·G, where G is a DH-type block mod 13. Off-line pairs are born at an exceptional point (two zeros collide) whenever tan c passes a critical value of Z_L/Z_G. Below height 790 every zero stays on the line only for c ∈ (−0.0022, 0.0130) (scan step 0.0005; both ends exact by Newton). The thresholds are set by Lehmer-type pairs of L(s, χ₁₃) (normalized gaps 0.068, 0.086, 0.205). The block alone has about 56 off-line pairs below 800. This is the DH filter, now with numbers. **Dead (barrier; Bombieri–Hejhal, Selberg, Hejhal).**
- **Black hole in a box (N13, N1).** A mirror at cusp height a turns resonances into normal modes, the zeros of a^sΛ(2s) + a^{1−s}Λ(2s − 1). These are real (Lagarias–Suzuki) whenever every zero has 0 < Re ρ < 1. DH's box spectrum is real to height 200, despite its four off-line zeros there. A mod-13 DH-type block has an off-line box mode below height 160 only at its β = 1.18 zero. **Dead (prime-number-theorem level).**
- **Chirps, matched filtering, stationary phase (N8, N10).** These tools are already in zeta computation: Odlyzko–Schönhage, Hiary's quadratic chirp sums, Booker/Büthe explicit-formula Turing, Gabcke. Data on the critical line see a zero pair only through u = δ². Turing's certificate is already optimal in its number of witnesses. **Dead (engineering, already optimal).**
- **Rest.** The overtone law corresponds to the Gamma poles, which cancel the trivial zeros. Regge poles have no continuation in the horocyclic mode. Hellings–Downs positivity is covered by R419. PN/EOB resummation of the Riemann–Siegel remainder has no RH node. **Dead.**

## Survey (arXiv API, abstracts; queries in .scratch/explore25/survey1.sh, survey2.sh)
- **Black-hole perturbation theory.**
  - Chandrasekhar 1975 (RW/Zerilli isospectral, superpartner potentials W² ± W′).
  - Glampedakis–Johnson–Kennefick 1702.06459 (Darboux transformation; key role of algebraically special modes; Kerr analogue).
  - Ould El Hadj–Dolan 2609.17580 (a local defect breaks isospectrality; exceptional points) and 2610.27102.
  - Leaver's continued fractions.
  - Motl gr-qc/0212096 and Motl–Neitzke hep-th/0301173: monodromy gives Re ω → T_H log 3 with equally spaced overtones.
- **Quasinormal-mode (in)stability.**
  - Nollert 1996 (PRD 53, 4397).
  - Jaramillo–Macedo–Al Sheikh 2004.06434: the pseudospectrum shows overtones are unstable, the fundamental mode stable.
  - Review 2308.16227; physical significance 2404.01374.
  - Structured pseudospectrum 2611.00853: admissible perturbations reduce the instability.
  - Transient QNM sums 2406.06685: sums of M modes live ~ log M from non-orthogonality.
  - Green's function = branch cut + QNMs + prompt response: 2603.20490.
- **Data analysis.**
  - Matched filtering, stationary-phase chirps, Owen's template metric (1996).
  - Hellings–Downs (1983): Legendre coefficients ∝ (2l+1)/((l−1)l(l+1)(l+2)), l ≥ 2; NANOGrav 2411.13472.
  - PN Padé/EOB resummation: Damour–Iyer–Sathyaprakash 1998, Buonanno–Damour 1999.
- **Links to ζ found.**
  - 0905.4939 (near-horizon scalar plus area quantization gives the zeros) and 2004.09523 (gauged CPT of the dilation operator, Schwarzschild S-matrix): Berry–Keating xp class (R51).
  - 2608.23528: prime gas → black-hole metric; zeros "cannot play this role".
  - 1103.0834 (zeros and poles of partition functions) and hep-th/9903267 (Landau–Siegel zeros and black-hole entropy): no zero location.
- **ζ side hit by the transfers.**
  - Modular scattering: Faddeev–Pavlov 1972; Lax–Phillips 1976 (R5.5, R58); Gutzwiller 1983 (phase shift arg ζ(1+2ik)); Colin de Verdière 1983 (pseudo-Laplacian); Lagarias–Suzuki math/0412039 (constant term; zeros on the line for y ≥ 1); Phillips–Sarnak (cusp forms dissolve under deformation).
  - Linear combinations: Bombieri–Hejhal 1995 (conditional: almost all zeros on the line); Selberg; Hejhal 1311.4862; 2010.10490 (zeros right of 1/2 + 1/G(T) ~ K₀TG/√log G).
  - Zero computation: Odlyzko–Schönhage 1988; Hiary 0711.5005 and 1002.2443; Booker math/0507502; Büthe 1308.6704; Gabcke 1979; Berry 1995 (Riemann–Siegel high orders); Arias de Reyna 2011.
- **Not found (1990–2026):** any use of Chandrasekhar or Darboux maps, quasinormal-mode pseudospectra, matched filtering, template metrics or EOB resummation on ζ zeros, Z(t) or the Riemann–Siegel remainder.

## Transfers, tests, verdicts

### T1 Chandrasekhar–Darboux isospectrality = Maass weight shift (N13)
- **Claim (D1).**
  - Black hole: V_± = W² ± dW/dr* + β are intertwined by A_± = ∓d/dr* + W. The reflection amplitudes differ by a unimodular rational factor, (μ²(μ²+2) − 12iMω)/(μ²(μ²+2) + 12iMω) with μ² = (l−1)(l+2), up to convention. Its zero or pole is the algebraically special frequency Mω = −iμ²(μ²+2)/12 (−2i for l = 2). QNMs are shared.
  - Modular surface, weight k: φ_k(s) = I_k(s)ζ(2s−1)/ζ(2s) with I_k(s) = ∫((u+i)/|u+i|)^{−k}(1+u²)^{−s}du = (−1)^{k/2}√π Γ(s−1/2)Γ(s)/(Γ(s+k/2)Γ(s−k/2)).
  - The Maass raising operator R_k = iy∂_x + y∂_y + k/2 gives φ_{k+2} = φ_k·C_k with C_k(s) = (1−s+k/2)/(s+k/2). C_k is unimodular on Re s = 1/2, has a zero at s = 1 + k/2 and a pole at −k/2. For k ≥ 2 the zero is the holomorphic Eisenstein series of weight k + 2 (a lowest-weight point = the algebraically special mode). At k = 0 it cancels φ₀'s pole at s = 1 and leaves the quasimodular E₂*. The poles ρ/2 are common to all weights.
- **Test (t1_maass.py).**
  - The closed form matches numerical integration to ≤ 4·10⁻²⁷ for k = 0…6 at three points.
  - C_k matches to 12 digits, and |C_k(1/2+it)| = 1.
  - |φ_k(ρ₁/2 + 10⁻⁸)| = 5.15·10⁷, 5.17·10⁷, 5.25·10⁷ for k = 0, 2, 4: the same pole.
- **(a) DH and moved pair.**
  - Write ξ_mod = ξR with R rational. Then φ_mod = φ·R(2s−1)/R(2s), a unimodular rational factor: a finite complex-seed Crum–Darboux chain. So ζ (RH true) and the moved pair (RH false) lie in one Darboux class.
  - The weight shift acts only on the Gamma factor, and unitarity |φ(1/2+it)| = 1 needs only FE + reality. A formal φ_DH therefore carries the same structure.
- **(b) Restatement.**
  - Maass operators (Maass, Roelcke, Hejhal LNM 1001).
  - Darboux insertion works for finitely many levels only (R308).
  - Lax–Phillips/Pavlov–Faddeev causality ≡ RH (R5.5, R58).
  - Defect-broken isospectrality (2609.17580) is a Phillips–Sarnak deformation: arithmeticity is lost and the ρ/2 move.
- **Dead (known).**

### T2 Quasinormal-mode instability = structured pseudospectrum of an L-function in its FE class (N13, N8): the numeric test
- **Claim (D2).**
  - Class: even primitive characters mod 13, Gamma factor (13/π)^{s/2}Γ(s/2), F(s) = F(1−s), real on the line.
  - Real basis on the line: Z₁₂ (quadratic χ, an L-function), Z_A = Z₃ + Z₉ and Z_B = Z₄ + Z₁₀ (DH-type blocks; Z₉(t) = Z₃(−t) checked to 12 digits).
  - F_c = cos c·L(s, χ₁₃) + sin c·G_A has a zero on the line exactly where r(t) = Z₁₂/Z_A = −tan c.
  - At a critical point of r, the set Im(L/G) = 0 also crosses the line. When tan c passes a negative local minimum value |r*|, two real zeros collide (a square-root exceptional point, EP) and leave the line. They re-enter at negative local maxima.
  - So the structured ε-pseudospectrum {s : (L/G)(s) ∈ [−ε, 0]} is read off the critical values of one real function. This is the scalar analogue of {z : ‖(A − z)⁻¹‖ ≥ 1/ε} restricted to admissible perturbations (2611.00853).
- **Test.** Arb Hardy Z on 32,000 points to height 800, then Newton at 96 bits; residuals < 10⁻²⁵.
  - L(s, χ₁₃): 817 sign changes vs θ(800)/π = 816.24, mean S ≈ 0.1. All zeros are on the line.
  - Block A: 704 sign changes, about 56 off-line pairs below 800. The first is 0.776646 + 13.708061i; 6 of the 44 we located have β > 1, e.g. 1.182959 + 77.19157i.
  - Block B: 782 sign changes, about 17 pairs.
  - Record EP births (lowest height for each threshold):

| sign of c | (t*, \|c*\|) records below 800 |
|---|---|
| c > 0 (negative minima of r, 63 of them) | (13.71, 1.485) (31.79, 1.080) (48.85, 0.870) (77.31, 0.492) (94.17, 0.324) (151.40, 0.178) (203.18, 0.166) (277.41, 0.0920) (604.00, 0.01295) |
| c < 0 (positive maxima of r, 215 of them) | (8.04, 0.188) (15.77, 0.157) (19.17, 0.0956) (47.36, 0.0768) (89.57, 0.0259) (254.25, 0.0189) (573.77, 0.0120) (597.60, 0.00238) (710.14, 0.00223) |

  - **The thresholds are exact.**
    - c = +0.0128 gives only the on-line zero 603.992493; c = +0.0131 gives 0.5096283 + 604.001733i.
    - c = −0.0022 gives on-line zeros only; c = −0.003 gives 0.5161727 + 597.596360i and 0.5220071 + 710.136117i.
    - Square-root EP law: β − 1/2 = 0.0096 at Δc = 1.5·10⁻⁴ and 0.0649 at Δc = 7.05·10⁻³. The ratio is 6.7, against √47 = 6.9.
  - **Lehmer pairs set the stability radius.**
    - The two smallest thresholds sit on the two smallest normalized gaps of L(s, χ₁₃) below 800: 0.068 at 597.6 and 0.086 at 710.1.
    - The c > 0 record at 604.0 combines a 0.205 gap (max |Z₁₂| = 0.105) with Z_A = 8.10.
  - **Whole parameter range.** Scanning c ∈ [−1.570, 1.570] in steps of 0.0005 below height 790, all zeros stay on the line exactly for c ∈ [−0.0020, 0.0125]; every other c shows wrong-sign extrema, the signature of off-line pairs. Negative admixtures are far more unstable: about 25 pairs at c = −0.05 and 59 at c = −0.1, against 1 at +0.05 and 2 at +0.1.
  - **Shape.** Off-line zeros stay a fixed number of mean spacings from the line: (β − 1/2)·log(13t/2π)/2π ≈ 0.11–0.22 for c = 0.5 and 0.32–0.36 for block A, in every 200-band up to 800. This is the 1/log T scale of Selberg and Hejhal. It is not a Nollert–Price branch, where perturbed overtones jump to new open curves, far from the unperturbed overtones, along which the damping grows only logarithmically.
- **(a) DH and moved pair.**
  - The family ends at a DH-type function and is unstable for every c ≠ 0 at large enough height. That instability is the DH filter: a stability argument that uses only the FE class must fail here.
  - Moving a pair off the line is an EP at positive de Bruijn–Newman time (explore24).
- **(b) Restatement.**
  - Bombieri–Hejhal 1995: under GRH and a weak pair-correlation hypothesis, almost all zeros of such combinations lie on the line.
  - Selberg; Hejhal 1311.4862; 2010.10490.
  - EPs: explore24 (T2/T3). R2 (the DH filter).
  - Lehmer pairs play the same role here as in the Csordas–Smith–Varga bounds on Λ.
- **Dead (barrier, made quantitative).** Inside its FE class, RH for L(s, χ₁₃) below height 790 survives only admixtures c ∈ (−0.0022, 0.0130). This window shrinks to {0} as T → ∞, because arbitrarily close Lehmer pairs occur.

### T3 Black hole in a box = Lagarias–Suzuki / pseudo-Laplacian (N13, N1)
- **Claim (D3).**
  - A Dirichlet wall at cusp height a makes the problem self-adjoint (Colin de Verdière's pseudo-Laplacian). The normal modes are the zeros of B_a(s) = a^sΛ(2s) + a^{1−s}Λ(2s−1). All lie on Re s = 1/2 for a ≥ 1, except two real ones iff a > 4πe^{−γ} (Lagarias–Suzuki).
  - Mechanism: Hermite–Biehler (HB) for E(s) = a^sξ(2s). Each zero ρ, grouped with its FE partner 1 − ρ̄, multiplies |E#/E| by (|x−β+id|/|x+β+id|)·(|x−1+β+id|/|x+1−β+id|). This factor is < 1 on Re s > 1/2 when 0 < β < 1 and can exceed 1 when β > 1.
  - Equivalently, Re(ξ′/ξ)(1+iu) = Σ(1−β)/|1+iu−ρ|² > 0: Lagarias' positivity at σ = 1, which holds unconditionally.
  - So the mirror can see only zeros with β ≥ 1. (DH and block A have real Dirichlet coefficients in this normalization, a_n = 2 Re(κχ(n)), so the identity applies to them too.)
- **Test (t3_*.py).**
  - min_u Re(ξ′/ξ)(1+iu) over (0, 200]: ζ 0.0231; moved pair (14.13, 21.02 → 1/2 ± 0.3 + 17.578i) 0.0223; DH 0.0976.
  - DH has off-line zeros 0.80852 + 85.69935i, 0.65083 + 114.16334i, 0.57436 + 166.47931i, 0.72426 + 176.70246i (reproduced), so DH stays positive despite them.
  - Box modes of the entire form a^sξ_X(2s) + a^{1−s}ξ_X(2s−1) on (0, 100] for a = 1, 7.0555, 20 (sign changes / phase count): ζ 79/79.19, 141/141.35, 175/174.50; moved pair the same; DH 130/129.75, 192/191.92, 225/225.07. The phase is monotone, so every box mode is on the line.
  - Mod-13 block A: Re ξ′/ξ < 0 only on u ∈ (76.9, 77.5), around its zero 1.183 + 77.19i. There the box phase dips (−8.05) and the box gets the off-line mode 0.67836243 + 38.547714i.
- **(a) DH-correct but RH-blind.** DH and the moved pair have real box spectra; the mirror can see only β ≥ 1.
- **(b) Restatement.**
  - Lagarias–Suzuki 2006.
  - The HB ladder at the zero-free abscissa (explore22, R323).
  - explore24: absorber/laser points of the modular scattering matrix are the prime number theorem.
- **Dead.**

### T4 Chirps, matched filtering, stationary phase (N8, N10; N2, N6)
- **Claim.**
  - Z(t) = 2Σ n^{−1/2}cos(θ(t) − t log n) + R is a sum of chirps. Chirp n has frequency (1/2π)log(√(t/2π)/n) and enters at its stationary point t = 2πn².
  - The stationary-phase error is the Riemann–Siegel remainder: Gabcke's explicit bounds (1979); Arias de Reyna's coefficients; Berry 1995 (divergent, resurgent).
  - The gravitational-wave toolkit is already used:
    - FFT correlation over a template bank = the Odlyzko–Schönhage multi-evaluation.
    - Quadratic-phase chirp sums = Hiary's O(t^{1/3}) method (0711.5005; amortized 1002.2443).
    - The optimal matched window = explicit-formula Turing (Booker, Büthe, Platt–Trudgian).
- **Off-line-pair matched filter.**
  - Moving γ₁ < γ₂ to γ̄ ± iδ gives Z_mod(t) = Z_red(t)·((t−γ̄)² + u) with u = δ²; u = −g²/4 is the on-line pair. This is linear in u.
  - So any statistic on the line sees only u, with constant Fisher information. Certifying an off-line pair needs precision ~ min(δ², g²/4)·|Z_red|: explore23's δ² law.
  - "Faster than Turing" is impossible in witnesses: sign changes must cover N(T+H) − N(T) zeros, and Turing adds only O(log T). Only the per-evaluation cost can improve, and that is Odlyzko–Schönhage/Hiary.
- **(a) DH.** Z_DH superposes the same chirps (same θ, same log n) with other amplitudes. Its own template fits it, off-line zeros included. The moved pair is flagged only by ζ's own template, which is finite verification.
- **(b) Restatement.** Turing, Booker, Büthe, Odlyzko–Schönhage, Hiary, Gabcke, Berry–Keating; explore23.
- **Dead (no new test needed; T4 adds no numbers).**

### T5 Other items (no test)
- **Overtone law vs N(T) (N8, N7).** The monodromy law (equally spaced overtones) corresponds to the Gamma poles of the archimedean factor. In φ(s) = √π Γ(s−1/2)ζ(2s−1)/(Γ(s)ζ(2s)) they cancel the trivial zeros, so there is no overtone tower. The mean density of nontrivial zeros is the archimedean phase θ(t) (explore23).
- **Regge poles.** The horocyclic Fourier mode n plays the angular momentum. Watson–Sommerfeld over n is Poisson summation, and the Eisenstein coefficients |n|^{s−1/2}σ_{1−2s}(|n|) have no canonical continuation in n. So there are no Regge trajectories of zeros.
- **Completeness and transients.**
  - Green's function = QNMs + branch cut + prompt response is the explicit formula ψ(x) = x − Σx^ρ/ρ − log 2π − ½log(1 − x^{−2}).
  - Long-lived sums of short-lived modes (2406.06685) correspond to Littlewood's Ω±(x^{1/2} log log log x). Both are known.
- **Hellings–Downs / noise-weighted inner products (N1).** Generic positive-definite kernels, DH-blind (R419, explore17).
- **PN/EOB/Padé resummation of the Riemann–Siegel remainder.** Already resurgent (Berry 1995) and not the bottleneck of zero verification.

## Ranking (RH relevance × chance × gain)
| Transfer | Node | RH relevance | Chance | Gain | Verdict |
|---|---|---|---|---|---|
| T2 FE-class pseudospectrum | N13, N8 | barrier only | none for RH | the DH filter with numbers (stable window (−0.0022, 0.0130) below 790) | dead, known mechanism |
| T4 chirps / matched filter | N8, N10 | engineering | low | none: already in Odlyzko–Schönhage, Hiary, Turing | dead |
| T3 mirror / box | N13, N1 | low | 0 | none: PNT-level | dead, Lagarias–Suzuki |
| T1 Chandrasekhar = Maass | N13 | low | 0 | none: Gamma factor only | dead, known |
| T5 rest | N1, N7, N13 | none | 0 | none | dead |

## Suggested NOGO line
- explore25 (Orca, exploration: gravitational-wave astronomy and black-hole perturbation theory — chirp matched filtering, PN/EOB resummation, Chandrasekhar–Darboux isospectrality, quasinormal-mode pseudospectra, mirrors, Hellings–Downs → N13, N8, N10, N1) NO IDEA: Chandrasekhar's RW↔Zerilli map is the Maass weight shift (φ_{k+2} = φ_k(1 − s + k/2)/(s + k/2), Gamma factor only; a finite Darboux chain links ζ to a moved pair); quasinormal-mode instability is DH instability inside an FE class (cos c·L(s, χ₁₃) + sin c·DH-block: off-line pairs born at critical values of Z_L/Z_G; below height 790 all zeros stay on the line only for c ∈ (−0.0022, 0.0130), the ends set by Lehmer pairs; the block alone ≈ 56 pairs below 800); a mirror in the cusp makes modes real using only 0 < Re ρ < 1 (Lagarias–Suzuki; DH's box spectrum real to 200); chirp matched filtering is Odlyzko–Schönhage/Hiary/Turing, and on-line data see a pair only through δ² [ExploreGW](Riemann/Deep/ExploreGW/REPORT.md)

## Files
- [REPORT.md](REPORT.md) (this file), [derivations.md](derivations.md) (D1 Chandrasekhar–Maass, D2 EP births, D3 Hermite–Biehler box), [outputs.txt](outputs.txt) (verbatim results).
- [code/](code/):
  - t1_maass.py (weight-k integral, C_k).
  - t2_grid.py (Arb Hardy Z grid; writes t2_grid.npz, about 13 min).
  - t2_analyze.py, t2_neg.py, t2_look.py, t2_cscan.py (sign counts, critical values, Lehmer pairs, scan over c).
  - t2_newton.py, t2_newton2.py, t2_branch.py (off-line zeros).
  - t3_box.py, t3_dhcal.py, t3_boxA.py, t3_boxA2.py (box spectra, DH calibration).
- Run with `uv run --no-project --with python-flint --with numpy --with mpmath python3 <script>` from code/; t2_* after t2_grid.py.

RH relevance: N13 (barrier only: the FE class is unstable, so any spectral route must use the Euler product exactly; no new idea), also N8, N10, N1 (none).
