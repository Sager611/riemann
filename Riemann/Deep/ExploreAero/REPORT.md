NO IDEA: aerodynamics reduces to the Mellin calculus of the dilation generator (the finite Hilbert transform is i tanh(πH), skew-adjoint, so it adds nothing to Weil's quadratic form; its even partner is the archimedean term itself), and Munk's theorem needs a drag operator that is diagonal in Glauert modes, which Weil's log-order archimedean term is not. The exact Munk restatement of window positivity is the known rank-one pole analysis, and it is DH-blind on every window it can certify. Flutter is the Csordas–Smith–Varga/de Bruijn machinery (explore24).

# ExploreAero: aerodynamics and aeroelasticity → RH tree (explore26, Orca, 11 Nov 2026)

Scripts and outputs: [code/](code/). Derivations: [DERIVATIONS.md](DERIVATIONS.md). Every run took under 2 min and under 100 MB, single-threaded; numerics are double precision (orientation only), validated against WeilSimple's Arb matrices.

## Survey (arXiv API and keyless arXiv search, 11 Nov 2026; API hits in [code/survey_arxiv.txt](code/survey_arxiv.txt))
- **No link to ζ.** Zero arXiv hits for airfoil ∧ zeta, induced drag ∧ zeta, Hilbert transform ∧ "Riemann hypothesis", Koppelman ∧ Pincus. "Hilbert transform" ∧ "explicit formula" returns tomography papers only. The keyless search returned noise.
- **Finite Hilbert transform (FHT).** Koppelman–Pincus 1959 diagonalise it by a Mellin transform: simple absolutely continuous spectrum [−i, i]. Widom 1960 and Okada–Elliott 1991 give the L^p theory: onto with a 1-dimensional kernel for p < 2, injective with codimension-1 range for p > 2, not Fredholm in L². Recent work:
  - Curbera–Okada–Ricker, rearrangement-invariant spaces: 1901.06334, the survey 2310.10228, and 2602.04844 on L^∞.
  - Bertola–Katsevich–Tovbis, multi-interval FHT with commuting differential operators (the Hilbert-transform analogue of prolate functions). 2610.35217 diagonalises linearised EIT through the adjacent-interval FHT, with continuous spectrum.
  - Al-Aifari–Defrise–Katsevich 1312.4727: SVD asymptotics of the truncated Hilbert transform. 2002.12723: inversion of the cosh-weighted FHT.
  - None of these involves arithmetic.
- **Aerodynamics.**
  - Kutta condition: Söhngen 1939 and Tricomi 1951 (inversion with weights √((1∓x)/(1±x))); a variational (Gauss least-constraint) closure, 2601.06495; an adjoint Kutta condition, 2606.06937; unsteady Kutta, 2205.08647.
  - Lifting line and induced drag: Munk 1921 (NACA Rep. 121), Glauert and Multhopp Chebyshev quadrature.
  - Unsteady aerodynamics: Theodorsen 1935 and Wagner/Küssner, with recent fits (2104.15122) and unsteady thin-airfoil theory (2609.05839).
  - Flutter: Zimmerman–Weissenburger 1964 (flutter margin), Nissim 1971 (aerodynamic energy).
- **Nearest number theory.**
  - The xp/dilation picture: Berry–Keating 1999, Connes 1999, the review 2606.24405.
  - Burnol's conductor operator log|x| + log|ξ| for the archimedean explicit formula (math/9902080, math/9809119).
  - Weil windows and CCM ground states: Yoshida 1992, Connes–Consani 2006.13771, CCM 2511.22755, Andrews–Connes–La Gatta 2611.13129.
  - Csordas–Smith–Varga 1994 and Rodgers–Tao.
  - No paper ties the logarithmic Laplacian (Chen–Weth 2019) to Weil's form; the one hit, 2606.04225, is unrelated.

## Transfers, tests, verdicts
| # | transfer (node) | numeric test | adversarial 1: DH / moved pair | adversarial 2: only a known equivalence? | verdict |
|---|---|---|---|---|---|
| T2 | Munk/Glauert: window positivity as minimum drag at unit lift (N1) | Weil form in Legendre and Glauert bases, a ≤ 0.65 | moved pairs pass every window the method can certify | yes: the rank-one pole analysis (Yoshida, Connes–Consani, CCM) | dead (best of three) |
| T1 | Kutta/airfoil operator as a selection rule for H (N13, N5) | the eigen-relation; spectrum of the Legendre–Galerkin truncation | identical for DH (sees only the gamma factor) | yes: Koppelman–Pincus, Berry–Keating xp, the Hilbert pair of R13.6 | dead |
| T3 | flutter margin and energy method for H_t (N2); Theodorsen positivity (N8) | Theodorsen C(k) positive-real check; no zeros on the physical sheet | works for moved pairs (explore24: t_e = 0.1507) | yes: CSV 1994, de Bruijn 1950, Lagarias | dead (pre-empted) |

| rank | transfer | RH relevance | chance | gain |
|---|---|---|---|---|
| 1 | T2 Munk | high (N1 is RH) | ≈ 0: no diagonal structure | engineering only (edge-adapted bases) |
| 2 | T1 Kutta/xp | medium (N13) | 0: skew-adjoint, linear Weyl law | none |
| 3 | T3 flutter | medium (N2) | 0: equals CSV and de Bruijn | none (explore24) |

### T2: Munk and Glauert on Weil's window (N1) — best transfer
- **Reformulation (exact, D2).** On [−a, a] (c = e^{2a}), write Q = 2ℓ² + R with ℓ(f) = ⟨f, cosh(x/2)⟩ and R = archimedean − primes.
  - Even sector: Q ≥ 0 ⟺ R ≥ 0 on ker ℓ and D_min = min{R : ℓ = 1} ≥ −2. In words: the Weil wing's minimum drag at unit lift is at least −2.
  - Odd sector: the pole term is −2⟨f, sinh(x/2)⟩², a "negative lift" (roll). So Q ≥ 0 ⟺ R > 0 and the maximal lift²/drag satisfies 2sᵀR⁻¹s ≤ 1.
  - The Munk optimisers are R⁻¹ℓ and R⁻¹s.
- **Validation ([code/t2a_out.txt](code/t2a_out.txt), [code/odd_conv_out.txt](code/odd_conv_out.txt)).** Our Legendre code converges to WeilSimple's Arb values:

  | c | quantity | ours (Legendre) | WeilSimple (Fourier, N = 40 → 320) |
  |---|---|---|---|
  | 2 | λ₀ even | 1.32932·10⁻³ (K = 32) | 1.33273·10⁻³ → 1.32942·10⁻³ |
  | 2 | λ₀ odd | 0.0731039 (K = 80) | 0.077437 → 0.073690 |
  | 3 | λ₀ even | 5.549·10⁻⁸ | 5.757·10⁻⁸ (N = 40) |

- **Munk's structure is absent ([code/t2_out.txt](code/t2_out.txt)).**
  - Prandtl's drag operator is exactly (π/2)k δ_kl in the Glauert modes sin kθ.
  - Weil's archimedean matrix in the same modes has normalised off-diagonals up to 0.22 at a = 0.3 and 1.21 at a = 0.6. Its diagonal is negative at k = 1 (A₁₁ = −1.53 at a = 0.6) and grows only like log k, against (π/2)k for Prandtl.
  - Diagonal dominance (Gershgorin) holds only at a = 0.1, on an 8-mode span with the tail ignored. It fails from a = 0.2, far inside the analytic windows a ≤ (log 2)/2 of Yoshida and Connes–Consani.
- **Optimum ≈ Munk's loading, but not to the needed accuracy.**

  | a | 0.1 | 0.2 | 0.3 | 0.347 (c = 2) | 0.549 (c = 3) | 0.65 |
  |---|---|---|---|---|---|---|
  | cos(ground state, elliptic) | 0.9982 | 0.9990 | 0.990 | 0.976 | 0.887 | 0.844 |
  | Munk margin 2 + D_min | 2.37 | 0.255 | 1.5·10⁻² | 2.4·10⁻³ | 8.5·10⁻⁸ | 5.9·10⁻¹¹ |

  - The margin tracks λ₀.
  - For a ≥ 0.2 the elliptic loading's drag at unit lift is 0.3–2.5 % above D_min, so Munk's "elliptic is optimal" fails. That excess (0.006–0.05) is 2 times (a = 0.3) to 10⁸ times (a = 0.65) the margin, so no near-elliptic ansatz can certify the sign.
  - Odd sector: the efficiency 2sᵀR⁻¹s rises 2.3·10⁻⁴, 0.141, 0.99068, 0.999329, 0.99999927 at a = 0.1, (log 2)/2, 0.5, (log 3)/2, 0.65. The optimum is ≈ the roll loading x√(1 − x²/a²) (cos 0.998 at a = 0.3).
- **Tip (Kutta-type) behaviour.** The optimum's tip load f₀(a)/f₀(0) falls only slowly with basis size: 0.300, 0.267, 0.242, 0.223, 0.208 for K = 6, 12, 24, 48, 96 at a = 0.1 ([code/t2b_out.txt](code/t2b_out.txt)).
  - This fits logarithmically weak vanishing at the edge: the drag symbol θ′ is of order 0, not Prandtl's order 1, so there is no √ tip law.
  - Engineering consequence: Rayleigh–Ritz upper bounds converge only algebraically. At a = 0.45 the relative change is 5·10⁻⁴ per doubling at K = 96.
  - CCM's odd spaces sin(2πjx/L) vanish at the edges (a Kutta-type constraint), so their odd λ₀ errs like 1/N: 0.0774 → 0.0737 for N = 40 → 320, against 0.07313 already at K = 10 with free edges. Lower-bound certificates are unaffected.
- **Fourier zeros of the optimum → zeta zeros (CCM's phenomenon, known).** At a = 0.65 they are 14.134, 21.022, 25.010, 30.435, 32.975, against γ = 14.135, 21.022, 25.011, 30.425, 32.935.
- **Adversarial 1 (moved pair).** Replace ±γ by ±γ ± iδ. λ₀ first turns negative at these window half-widths a (none = still positive at a = 0.65):

  | moved pair | δ = 0.02 | δ = 0.1 | δ = 0.25 |
  |---|---|---|---|
  | γ₁ | 0.549 | 0.50 | 0.45 |
  | γ₂ | 0.65 | 0.60 | 0.549 |
  | γ₅ | none | none | none |

  Every window a Munk/Gershgorin argument can certify (a ≤ 0.1) is passed by all these moved functions. Windows that see a pair at height γ grow with γ, so this is finite zero verification (the N1 barrier).
- **Adversarial 2.** The reformulation is Weil positivity with the pole term handled as a rank-one perturbation, i.e. a Schur complement: known.

### T1: airfoil operator, Kutta condition and xp (N13; N5)
- **Exact identity.** In u = (1+x)/(1−x), T = (1/π)PV∫·/(t−x) equals i tanh(πH), where H = −i(u d/du + 1/2) is the Berry–Keating dilation generator. Check: T g_τ = i tanh(πτ) g_τ for g_τ = ((1+x)/(1−x))^{iτ}/√(1−x²), to 12 digits at τ = 0, 0.35, 0.7, 1.5 ([code/t1_out.txt](code/t1_out.txt)).
- **The Kutta mode.** 1/√(1−x²) is the generalized eigenvector at τ = 0, i.e. at s = 1/2.
  - The Kutta condition picks the side of the pole of 1/tanh at s = 1/2, and the circulation is its residue (D1).
  - In L^p the index jumps from +1 to −1 at p = 2. This is Beurling's L^p/Nyman mechanism with a single "zero" at s = 1/2; N5 already has the arithmetic version.
  - Discrete Kutta mode: the exact zero eigenvector for odd K has overlap 0.915, 0.897, 0.884 with the projected 1/√(1−x²) at K = 41, 161, 641.
- **Wrong Weyl law.** In the Legendre–Galerkin (Multhopp-type) truncation, τ_k = artanh(μ_k)/π has a linear count. The slope dN/dτ is 4.35, 4.94, 5.73, 6.33, 6.95 for K = 41 … 641 (≈ 0.95 log K + 0.8), and N(0.25, 0.5, 1, 2, 4) = 1, 3, 7, 12, 23 at K = 641. This is box-regularised xp, not Riemann–von Mangoldt's (τ/2π) log(τ/2πe).
- **No contribution to positivity.** T is skew-adjoint, so ⟨f, Tf⟩ = 0 for real f. Its even partner, the same causal kernel 1_{t>0}/(2 sinh(t/2)) with symbol −ψ(1/2 − iτ), is the complex-place archimedean Weil term. ζ's own kernel is e^{x/2}/sinh x = 1/(2 sinh(x/2)) + 1/(2 cosh(x/2)) (D1). So "an airfoil operator that reproduces part of Weil's form" is only the archimedean term itself (Burnol's conductor operator), with the airfoil operator as its Kramers–Kronig conjugate: the R13.6 relation at the archimedean place.
- **Adversarial tests.** Everything here depends only on the gamma factor, so DH and moved-zero ξ behave identically. It only restates known results.

### T3: flutter, energy method, Theodorsen (N2, N8)
- **Flutter margin.** Coupled-mode flutter is the collision of two real modes. For H_t, gap² = 8(t − t_c) is the Zimmerman–Weissenburger margin (a discriminant, linear in the parameter). Its extrapolated onset is the CSV Lehmer-pair bound. explore24 computed it: −7.11203·10⁻⁴ (81-zero ODE) against −7.11206·10⁻⁴ (CSV) for Lehmer's pair at 7005. These are lower bounds on Λ only, moot since Λ ≥ 0.
- **Energy method.** Nissim's aerodynamic energy gives d(Im z)²/dt = −2 for an isolated pair, which is de Bruijn's strip theorem (0.2788 with the Lean half-plane, N2). The p-k and V-g sweeps are the argument principle (R404).
- **Theodorsen.** C(k) = H₁⁽²⁾/(H₁⁽²⁾ + iH₀⁽²⁾) = K₁/(K₀ + K₁)(ik) to 1.3·10⁻¹⁶. On k ∈ (0, 20], 1/2 < Re C ≤ 1, Im C ≤ 0 and |H₀⁽²⁾| < |H₁⁽²⁾|. K₀ + K₁ has 0 zeros on the physical sheet (keyhole winding −3·10⁻¹⁵; [code/t3_out.txt](code/t3_out.txt)). This is a passive, positive-real transfer function; the ξ analogue is Lagarias' Re ξ′/ξ > 0 (explore22).
- **Wagner/Küssner.** 1 − φ is completely monotone (R409, R419).
- **Adversarial tests.** The method also works for moved pairs (explore24: t_e = 0.1507). It restates CSV and de Bruijn.

## Suggested NOGO line
- explore26 (Orca, exploration: aerodynamics: thin-airfoil finite Hilbert transform, Kutta condition, lifting line/Glauert/Munk, Theodorsen/Wagner, flutter → N1, N13, N2, N8) NO IDEA: the airfoil operator is i tanh(πH) of the dilation generator (12 digits; Kutta mode = the s = 1/2 resonance; Multhopp truncations count N ≈ 0.95 log K·τ, not Riemann–von Mangoldt) and skew, so it adds nothing to Weil's form, whose archimedean term is its even Kramers–Kronig partner; Munk recasts window positivity as minimum drag at unit lift ≥ −2 (odd: lift²/drag efficiency ≤ 1, 0.99999927 at a = 0.65), the known rank-one pole analysis, but the log-order drag is far from diagonal in Glauert modes (normalised off-diagonals up to 1.21), Gershgorin works only at a = 0.1, the optimum's tip load vanishes only logarithmically (CCM's edge-vanishing odd basis converges like 1/N), and moved pairs fail only from a ≈ 0.45 (γ₁): DH-blind; flutter = CSV/de Bruijn (explore24), Theodorsen = positive real (explore22) [ExploreAero](Riemann/Deep/ExploreAero/REPORT.md)

## Files
- [DERIVATIONS.md](DERIVATIONS.md): D1 airfoil = i tanh(πH), Kutta, L^p index, the causal-kernel pair; D2 Munk/Glauert and the Weil wing; D3 flutter ↔ H_t.
- [code/t1_airfoil.py](code/t1_airfoil.py) → t1_out.txt: eigen-relation, Legendre–Galerkin spectrum, discrete Kutta mode.
- [code/weilleg.py](code/weilleg.py): Weil form for ζ in even Legendre modes on [−a, a] (x-space archimedean formula).
- [code/t2a_validate.py](code/t2a_validate.py) → t2a_out.txt (and the eigs_quick lines): validation against WeilSimple.
- [code/t2_munk.py](code/t2_munk.py) → t2_out.txt: Munk margin, Glauert structure, ground state, Fourier zeros, moved pairs, Gershgorin.
- [code/t2b_odd.py](code/t2b_odd.py) → t2b_out.txt: odd-sector Munk efficiency, tip value against K. [code/odd_conv.py](code/odd_conv.py) → odd_conv_out.txt: odd convergence against CCM's basis.
- [code/t3_theodorsen.py](code/t3_theodorsen.py) → t3_out.txt: Theodorsen positive-real check and zero count.
- [code/survey_arxiv.txt](code/survey_arxiv.txt): arXiv API query results.
- Run with: uv run --no-project --with mpmath --with numpy --with scipy python3 <script> (python-flint for eigs_quick).

RH relevance: N1 (best transfer, T2); also N13, N2 and N8; no node advanced.
