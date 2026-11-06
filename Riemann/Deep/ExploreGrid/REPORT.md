NO IDEA: every transfer from grid-cell coding and neural-field theory lands on a known result. The Euler phases t log p mod 2π form a grid code whose capacity, aliasing and random-torus model are Dirichlet/Kronecker, Bagchi recurrence and Bohr–Jessen: correct on Davenport–Heilbronn (DH), but they give zero density only. A neural field with Weil's kernel is the window Gram method, and its Turing mode has a Fourier transform with real zeros automatically (Krein; DH too, even at c = 40, where λ_min = −0.024). Grid optimality is a property of the real axis: the class-number-2 lattice 2x²+xy+2y² beats the square lattice on packing, θ and Z(2), yet its Epstein zeta has the off-line zero 0.75807 + 24.48282i.

# ExploreGrid: grid cells, phase codes, continuous attractors and neural fields → RH tree (explore21, Orca, 6 Nov 2026)

Normalisations. Lattices have covolume 1. Z_Λ(s) = Σ' |v|^{−2s}, θ_Λ(α) = Σ e^{−πα|v|²}. Hexagonal: Z_H(s) = 6(√3/2)^s ζ(s)L(s, χ₋₃).
For the forms of discriminant −15: Z_{2x²+xy+2y²} = ζL(χ₋₁₅) − L(χ₋₃)L(χ₅) and Z_{x²+xy+4y²} = ζL(χ₋₁₅) + L(χ₋₃)L(χ₅). Both identities hold coefficient by coefficient for n ≤ 3000 (logs/t1_verify.log).
The DH function is Σ a(n)n^{−s} with a(n) = 1, κ, −κ, −1, 0 repeating mod 5 (Titchmarsh §10.25). Its Weil form comes from explore19's `wmat.py`: gamma factor Γ((s+1)/2), conductor 5, no pole.
Windows are [−(ln c)/2, (ln c)/2], as in TREE N1.

## Dictionary
| neuroscience | RH side |
|---|---|
| grid module: period λ, phase x/λ mod 1 | prime p: phase t log p mod 2π, period 2π/log p; prime powers p^k are the module's harmonics |
| modular (residue) code with integer periods 1..N; CRT capacity lcm(1..N) | lcm(1..N) = e^{ψ(N)}, so RH ⟺ log lcm(1..N) = N + O(√N log²N) (N10, restated) |
| capacity grows exponentially in the number of modules; first alias | Dirichlet/Kronecker: about log T/1.26 primes localize t ≤ T at 10 % phase tolerance (T1 test); same scale as the Ω-results (Montgomery; Aistleitner–Mahatab–Munsch) |
| ambiguity function Σ wᵢ cos(2πΔ/λᵢ) | Re P(σ+iΔ) ≈ log abs ζ(σ+iΔ): aliases are large values; a zero would be a perfect anti-alias (all p^{−iΔ} ≈ −1) |
| alias = shift under which the code repeats | almost period; Bagchi strong recurrence (≡ RH, NOGO R3.4, R14.6) |
| error correction (decode to the nearest strand; Sreenivasan–Fiete) | strand separation is abs(Σ a_p log p), i.e. gaps between S-units (Tijdeman/Baker); proven, not RH-related |
| toroidal population activity (Gardner et al. 2022) | Bohr compactification T^∞ (dual of ℚ_{>0}^×); the Bohr–Jessen/Bagchi random model (R36) |
| path integration: phases advance linearly in x | t ↦ (t log p)_p is a Kronecker flow; the archimedean phase θ(t) of Hardy's Z is not periodic, so the zeros are not almost periodic in t and recurrence lives only in σ > 1/2 |
| infinitely many weak modules | amplitudes p^{−1/2}, Σ p^{−1} = ∞: Selberg CLT variance (1/2) log log T; no finite set of modules decodes the zeros (Gonek–Hughes–Keating split) |
| place field as a Fourier sum of grid cells (Solstad–Moser–Einevoll 2006) | explicit formula: zero density = smooth − Σ Λ(n)n^{−1/2}·2cos(t log n)/2π (Guinand–Weil); the zero measure as a quasicrystal (R2.6, explore5) |
| sawtooth grid cells {y/k} with integer periods k | Báez-Duarte/Nyman–Beurling with y = 1/x: RH ⟺ the constant is a limit of linear readouts (N5) |
| Amari u_t = −u + w*f(u), linearisation −I + f′W | Weil form Q(f) = ⟨f, Wf⟩ on a window: stability at infinite gain is window positivity (N1; R402 Gram method) |
| Turing instability: ŵ(k) crosses a threshold | Bochner: ŵ of Weil's kernel is Σ δ_γ, positive iff RH; an off-line pair is a growing term 4cos(γ_r x)cosh(γ_i x), seen only on windows c ≳ height/4.4π (explore19) |
| ring attractor: kernel on a circle, eigenvalues = Fourier coefficients | periodized window form: Chekmachev's lattice minorant (the N1 record method) |
| critical (Turing) mode on a bounded domain | ground state of the window form (CCM); its Fourier transform has only real zeros for every kernel (Krein; T2) |
| Mexican-hat kernels | the Fourier transform of a difference of Gaussians has zeros k² = (log C + 2πim)/(b − a), infinitely many non-real; neural kernels give no new Laguerre–Pólya kernels (N6/N2) |
| Evans function of a bump | Fredholm determinant / argument principle (R404, explore18) |
| three-wave resonance, hexagon selection | no resonances among the log p (they are linearly independent), only harmonics k log p: Bohr independence |
| hexagonal optimality: packing, Fisher information, θ (Mathis–Stemmler–Herz; Bétermin; Montgomery 1988) | Z_H factors as ζ·L(χ₋₃) because h(−3) = 1 (Euler product), not because H is optimal (T3) |
| optimal scale ratio √e (Wei–Prentice–Balasubramanian) | radix economy; log p_{i+1}/log p_i → 1, so the prime code is maximally redundant; nothing about zeta |

## Survey (arXiv; the export API returned HTTP 429 all session, so arxiv.org/search and abs pages were used: .scratch/explore21/axsearch.py, absf.py)
- Modular/phase grid code. Fiete–Burak–Brookings, J. Neurosci. 2008 (residue-number-like code). Sreenivasan–Fiete, Nat. Neurosci. 2011 (analog error correction; capacity exponential in the number of modules). Mathis–Herz–Stemmler, PRL 2012 (resolution exponential in neurons). Wei–Prentice–Balasubramanian 1304.0031 (ratio √e; 1.4–1.7 for realistic neurons). Zlokapa–…–Fiete–Tegmark–Chuang 2202.12887 (fault-tolerant networks from grid codes). Klukas–Lewis–Fiete 2020. 2311.04872 (residue numbers in hyperdimensional computing); 2406.18808.
- Lattices and optimality. Mathis–Stemmler–Herz 1411.2136: hexagonal is best in 2D, FCC in 3D, by information theory. Bétermin 2010.08264: the grid-cell Fisher information is built from translated theta functions, and the triangular lattice is *not* maximal for every Gaussian width. Bétermin–Faulhuber 2007.15977 (maximal theta functions). Montgomery, Glasgow Math. J. 1988 (θ_Λ minimal at H for every α); Rankin, Cassels, Ennola, Diananda (Z_Λ(s) minimal at H for s > 0).
- Torus and decoding. Gardner et al., Nature 2022 (grid-module activity lies on a 2-torus); 2510.16216 (decoding by path lifting); 2510.02853; 2501.19262; 2606.27783 (CANN toolkit).
- Neural fields. Amari 1977; Ermentrout–Cowan 1979 and Bressloff–Cowan–Golubitsky–Thomas–Wiener 2001 (hexagons by equivariant bifurcation); Coombes 2005; Bressloff 2012 (reviews); 1803.01380 (Evans functions with oscillatory coupling).
- Epstein zeros. Davenport–Heilbronn 1936 II: h(D) > 1 gives infinitely many zeros in σ > 1. Potter–Titchmarsh 1935; Bombieri–Hejhal; Lamzouri 1907.06387 (N(σ₁, σ₂, T) asymptotics for h ≥ 2); 1811.01613, 1204.6297; Rezvyakova 2411.18492 (positive proportion on the line); McPhedran 1601.01724 and Bétermin–Šamaj–Travěnec 2110.09368 (curves of off-critical zeros for rectangular lattices).
- Links to RH. arxiv.org/search for "grid cells" with zeta, number theory, Diophantine, Kronecker, Riemann zeros or Euler product; "neural field" with Riemann zeta; "Turing pattern" with "Riemann hypothesis"; "phase code" with incommensurate: **0 hits** each. The only overlap found is Bétermin's theta-function grid model, which is about optimality, not zeros.

## T1. The Euler-phase grid code (N11, N10, N7). Best transfer; dead
Transfer. Treat each prime as a grid module with phase t log p mod 2π. Then (i) how many modules localize t, (ii) what aliases are, and (iii) what the population "posterior" Π_p(1 − p^{−s})^{−1} predicts about zeros.

Test (i): capacity (code/t2_gridcode.py a, exact interval sieve, logs/t2_a.log). τ_N(δ) is the first lag > 0 at which all N phases are within δ of 0.

| δ | N = 2 | 4 | 6 | 8 | 10 | 12 | fitted log τ_N per extra prime (theory log(1/2δ)) |
|---|---|---|---|---|---|---|---|
| 0.10 | 17.5 | 171 | 2447 | 5.0·10⁴ | 1.3·10⁵ | 5.0·10⁶ | 1.26 (1.61) |
| 0.05 | 45.6 | 898 | 5.5·10⁴ | 5.5·10⁵ | 1.3·10⁸ | — | 1.86 (2.30) |

So about N ≈ 2 + (log T − 2.9)/1.26 primes (δ = 0.1) separate all heights up to T; at T = 3·10¹² that is 23 primes (p ≤ 83). The slope runs below log(1/2δ), presumably because the aliases have finite width and the log p have near-relations (2¹⁹ ≈ 3¹²). This is Dirichlet's simultaneous-approximation scale. It is the mechanism behind the Ω-theorems, and it localizes t, not the zeros. Zeros at height T need the whole prime sum: the variance Σ p^{−1} diverges.

Test (iii): random torus against t-averages at σ = 0.75 (code/t2_gridcode.py b; t ∈ [1000, 2000), step 0.025, 40,000 points; random model with 6·10⁴ samples and primes < 10⁵; logs/t2_b075.log).

| ε | P(abs ζ < ε): t-avg / random | P(abs DH < ε): t-avg / random | DH: P/ε², t-avg / random |
|---|---|---|---|
| 0.5 | 0.153 / 0.154 | 0.067 / 0.071 | 0.27 / 0.28 |
| 0.2 | 2.8·10⁻³ / 1.7·10⁻³ | 9.1·10⁻³ / 8.6·10⁻³ | 0.23 / 0.21 |
| 0.1 | 0 / 0 | 1.9·10⁻³ / 2.1·10⁻³ | 0.19 / 0.21 |
| 0.03 | 0 / 0 | 1.8·10⁻⁴ / 1.0·10⁻⁴ | 0.19 / 0.11 |

On the grid, min abs ζ = 0.145 and min abs DH = 0.0032. The standard deviation of log abs ζ is 0.670 for the t-average and 0.670 for the random model.
- Adversarial test 1 (DH). Passed. The model gives DH a positive two-dimensional value density at 0, P ≈ 0.2ε², because two Euler products interfere. It gives ζ a left tail lighter than any power. So it does not wrongly "prove" RH for DH. But its theorem-level content is the Bohr–Jessen value distribution (1932) and the density dichotomy: N(σ, T) = O(T) (Bohr–Landau 1914) and o(T) (Carlson 1920) for ζ, against ≫ T zeros in every strip 1/2 < σ₁ < σ < σ₂ < 1 for DH (zeros in σ > 1: Davenport–Heilbronn 1936; in the strip: joint universality, Voronin) and for Epstein zetas with h ≥ 2 (Gonek–Lee; Lamzouri 1907.06387).
- Adversarial test 2 (restatement). Yes. Capacity is Dirichlet/Kronecker. Aliasing is Bagchi recurrence, which DH also has under GRH for χ mod 5 by joint universality; recurrence implies RH only through the Euler-product density theorem (R3.4, R14.6). The random torus is Bohr–Jessen (R36: one deterministic sample). The integer-period residue code is ψ(N). Verdict: dead.

## T2. A neural field whose lateral kernel is Weil's distribution (N1, N13). Dead
Transfer. Linearise Amari's field at a homogeneous state: −I + f′W. At infinite gain it is stable on a window iff the window Weil form is ≥ 0. The Turing critical mode is the window ground state. If Q − λ_min‖·‖² ≥ 0 is singular and λ_min is simple (globally; here the ground state is even), Krein's extension writes it as ∫ abs(F)² dμ with μ ≥ 0. The ground state's F vanishes on supp μ, and any non-real zero w could be flipped to w̄, (z − w̄)/(z − w), giving a second minimizer. So **every zero of the mode's Fourier transform is real, whatever the sign of λ_min.** This is the continuous form of Carathéodory–Toeplitz/Pisarenko.

Test (code/t3_turing.py; even block, Galerkin N = 0.6 c ln c + 20, logs/t3_*.log). The box count is the argument principle on [0, Z] × [−2, 2].

| kernel, c | λ_min | real zeros in (0, Z) / zeros in box | match to the true on-line zeros |
|---|---|---|---|
| ζ, 10 | 5.73·10⁻⁴² | 8 / 8.000 (Z = 45) | all 8 to the bisection resolution (≈ 5·10⁻¹⁵) |
| ζ, 20 | 1.75·10⁻⁷⁹ | 21 / 21.000 (Z = 80) | all 21 zeros below 80 to ≈ 5·10⁻¹⁵ |
| DH, 20 | 4.52·10⁻¹⁸ | 27 / 27.000 (Z = 60) | 5.6·10⁻¹⁷, 5.4·10⁻¹⁵, 1.5·10⁻¹³, …; first 15 to < 2·10⁻³ |
| DH, 40 | **−0.0242** | 69 / 69.000 (Z = 120) | only 10⁻⁴ to 1, with extra zeros (0.90, 2.72, …) and a real zero at 85.729 under the off-line pair 85.6993 ± 0.3085i |

- Adversarial test 1 (DH). The mode property holds for DH both before its failure (c = 20: its first zeros to 5.6·10⁻¹⁷, the phenomenon CCM see for ζ, on DH's scale c/5) and after it (c = 40). Window positivity itself is DH-blind up to c ≈ 30.7 (explore19; Elston). So the transfer proves nothing DH lacks.
- Adversarial test 2 (restatement). Yes: the Gram/Hopfield certificate method (R402), the ring-attractor (periodized) analysis = Chekmachev's minorant, CCM's ground state (N13, WeilSimple), explore19's near-null "superdirective" modes, and explore17's E/I-type splits (W − pole Markov, pole anti-Markov). Mexican hats add no Laguerre–Pólya kernels (Dictionary). Verdict: dead.

## T3. Hexagonal grid optimality and the zeros of lattice zeta functions (N13). Dead
Transfer. Grid cells select H, which optimizes packing, Fisher information and θ. Z_H contains ζ. Does optimality force RH?

Test A (code/t1_lattice.py A, logs/t1_table.log):

| form | D | h | Euler product | packing | θ(1) | Z(2) |
|---|---|---|---|---|---|---|
| x²+xy+y² (H) | −3 | 1 | yes | 0.9069 | 1.159595 | 5.7834 |
| **2x²+xy+2y²**, τ = (−1+i√15)/4 | −15 | 2 | **no** | **0.8112** | **1.171945** | **5.9291** |
| x²+y² (square) | −4 | 1 | yes | 0.7854 | 1.180341 | 6.0268 |
| x²+xy+2y² | −7 | 1 | yes | 0.5937 | 1.221131 | 6.6319 |
| x²+xy+41y² | −163 | 1 | yes | 0.1230 | 2.526573 | 88.80 |

Test B (code/t1_lattice.py B, t1_verify.py, logs/t1_B.log, t1_verify.log). For the rhombic lattice 2x²+xy+2y², the argument principle counts 199 zeros up to T = 200, against 179 sign changes on the line, so 10 off-line pairs. Nine were located. The first is **0.75807089343824 + 24.48282139425657i** (checked with mpmath at 30 digits, abs Z = 2·10⁻³²). Further zeros: 0.5359 + 75.612i, 0.7784 + 98.829i, …
The principal form x²+xy+4y² has 100 zeros against 70 on the line up to T = 100. Its off-line zeros include 0.80001 + 12.03860i and **1.02597261621875 + 61.42279830878940i, with Re s > 1**. The Euler-product ratio there is L(χ₋₃)L(χ₅)/(ζL(χ₋₁₅)) = −1.0006 + 0.0004i over p < 2·10⁶ (exactly −1 is required).
Control: H has 35 zeros, all 35 on the line, up to T = 60.

Test C (code/t1_mellin.py). Write D = θ_Λ − θ_H ≥ 0 (Montgomery). Then M(s) = π^{−s}Γ(s)(Z_Λ − Z_H)(s) = ∫ D α^{s−1} dα is entire with abs M(s) ≤ M(σ). That is the whole content of optimality at complex s. At the off-line zero above: abs M = 1.5·10⁻¹⁶ ≤ M(σ) = 0.0125. At L(χ₋₃)'s first zero 1/2 + 8.0397i: 2.2·10⁻⁵ ≤ 0.0124. The Γ-decay swallows everything.
- Adversarial test 1 (Epstein). Failed by optimality, as it must be. The D = −15 rhombic lattice lies on the arc between the square and hexagonal lattices and beats the square lattice on all three functionals, yet it has off-line zeros from height 24.5. The poor grid codes D = −7 and −163 have Euler products. "Has an Euler product" is not monotone in any grid-code functional.
- Adversarial test 2 (restatement). Yes: Montgomery/Rankin optimality, the Hecke factorization for h = 1, Davenport–Heilbronn II, R073 (Epstein never easier) and R046 (deformation along τ: zeros leave the line only by collision; at H and i every lattice function is stationary, so the zeros move only at second order). Verdict: dead. The low off-line zeros above are an adversarial example, not a claim: off-line Epstein zeros at low height are classical (Potter–Titchmarsh, McPhedran, Bétermin et al.).

## Ranking (RH relevance × chance × gain)
| transfer | RH relevance | chance | gain | why it dies |
|---|---|---|---|---|
| T1 Euler-phase grid code | medium: the only transfer that sees the Euler product | ≈ 0 | none | Bohr–Jessen/Kronecker/Bagchi; density only |
| T2 Weil neural field / Turing mode | high (N1 ≡ RH) | 0 | none | = window Gram method; real-zero mode automatic (Krein) |
| T3 hexagonal optimality | low–medium (Z_H ⊃ ζ) | 0 | none | optimality lives on the real axis; RH-type = class number one |

## Suggested NOGO line
- explore21 (Orca, exploration: grid cells, modular/phase grid codes, continuous attractors, Amari/Wilson–Cowan neural fields, Turing/hexagon selection → N11, N10, N1, N13, N6) NO IDEA: the Euler phases t log p mod 2π are a grid code whose capacity is Dirichlet's (≈ 2 + (log T − 2.9)/1.26 primes localize t at 10 % tolerance, 23 primes at 3·10¹²), whose aliases are Bagchi recurrence (DH recurs too) and whose random torus is Bohr–Jessen (σ = 0.75: ζ none below 0.1 in t-average and model alike, DH P(abs < ε) ≈ 0.2ε²), so it is DH-correct but density-only; a neural field with Weil's kernel is the window Gram method, and its Turing mode has a Fourier transform with only real zeros for every kernel (Krein: DH at c = 40 with λ_min = −0.024 has 69/69 real zeros, one at 85.729 under its pair 85.699 ± 0.309i); hexagonal grid optimality is a real-axis property (abs M(s) ≤ M(σ); slack 10¹⁴ at height 24): the rhombic h = 2 lattice 2x²+xy+2y² beats the square lattice on packing 0.811, θ(1) 1.1719 and Z(2) 5.929, yet its Epstein zeta has the zero 0.75807 + 24.48282i (x²+xy+4y²: 1.02597 + 61.42280i, Re s > 1) [ExploreGrid](Riemann/Deep/ExploreGrid/REPORT.md)

## Files
- code/t1_lattice.py (A: lattice table; B: Epstein zero counts by argument principle and sign changes, off-line zeros by Newton; python-flint), code/t1_verify.py (coefficient identities; mpmath recheck; Euler-ratio check for σ > 1), code/t1_mellin.py (Montgomery–Mellin bound).
- code/t2_gridcode.py (a: alias sieve; b: ζ/DH small values against the Bohr–Jessen model), code/t3_turing.py (window ground states for ζ and DH; uses Riemann/Deep/ExploreAntenna/code/wmat.py).
- logs/: t1_table, t1_B (hex 60, Q2 200, Q1 100), t1_B60 (Q2 rerun with the corrected s = 3 tail check: 12 digits), t1_verify, t1_mellin, t2_a, t2_b075, t3_zeta10, t3_zeta20, t3_dh20, t3_dh40 (.log).
- Run: `uv run --no-project --with python-flint --with mpmath --with numpy python3 code/<script>`; every job was single-threaded and under 160 MB.

RH relevance: N11/N10 (grid code), N1/N13 (neural field), N13 (hexagonal lattice): none advanced.
