NO IDEA: Beckner's sharp logarithmic-uncertainty constant ψ(1/4) − ln π = 2θ′(0) is the minimum of Weil's archimedean symbol (Burnol's conductor operator; the killing constant of explore17's Markov form), so its window versions prove archimedean positivity only for half-width < 0.0856 (entropic versions < 0.0122), below Yoshida's (ln 2)/2, no prime-blind bound passes c = 2, and its ground-state prediction (simple, even, nodeless; Perron–Frobenius only for λ² ≤ ρ = 1.3247, the plastic number) holds just as well for L(s, χ₅), Davenport–Heilbronn and moved-pair models.

explore32 (Orca exploration, 11 Nov 2026, 22:31–23:20 UTC). Code [code/](code/), logs [logs/](logs/). Conventions of ../WeilZeta/DERIVATION.md:
supp f ⊂ I = [−a, a], a = (ln c)/2 (CCM: λ² = c), F(t) = ∫ f e^{itx} dx, φ(t) = −ln π + Re ψ(1/4 + it/2) = 2θ′(t),
Q = pole + Arch − primes, pole = 2⟨f, cosh(x/2)⟩² − 2⟨f, sinh(x/2)⟩², Arch = (1/2π)∫ φ|F|², w_n = Λ(n)/√n.

## 1. Transfer

Tools. (a) Beckner, Proc. AMS 123 (1995) 1897: ∫ ln|y| |k|² + ∫ ln|η| |k̂|² ≥ (ψ(1/4) − ln π)‖k‖² on ℝ, k̂(η) = ∫ k e^{2πiyη} dy,
the α-derivative at 0 of the sharp Pitt inequality ∫ |η|^{−α}|k̂|² ≤ C_α ∫ |y|^α |k|², C_α = π^α [Γ((1−α)/4)/Γ((1+α)/4)]²
(the AMS page is behind a bot check; the statement is re-derived in S0 and checked numerically in T1). (b) Hirschman–Beckner
H(|f|²) + H(|f̃|²) ≥ ln(e/2), f̃(ξ) = ∫ f e^{2πixξ}. (c) Maassen–Uffink / Hirschman on a circle: H(p) + H(q) ≥ ln(period).
Bridge read in the source (Burnol, arXiv math/9902080 Thms E2, G1, H2, F2; math/9811040; Haran 1990): the local term of the explicit
formula at a place ν is the spectral analysis of the conductor operator H = ln|x|_ν + ln|y|_ν on L²(K_ν); H is essentially
self-adjoint, bounded below, commutes with dilations, Fourier transform and inversion, and in the Mellin picture is multiplication
by the derivative of Tate's local gamma function on the critical line; at p, on invariant functions, H is bounded with spectrum
[−2 ln p/(√p − 1), 2 ln p/(√p + 1)]. Burnol states only "bounded below" and never cites Beckner (rg of both TeX sources).
Nodes: N1 (analytic window positivity, against Yoshida's (ln 2)/2 and Chekmachev's 1.925), N13 (CCM's first missing step).

## 2. Statements

S0 (identity; re-encoding). For f ∈ L²(ℝ), k(y) = |y|^{−1/2} f(ln|y|) (even, ‖k‖² = 2‖f‖²):
  Arch(f) = ½⟨k, (ln|y| + ln|η|) k⟩ = D‖f‖² + ∫_0^∞ κ(x) ‖f − f(· − x)‖² dx,   κ(x) = e^{x/2}/(e^x − e^{−x}) = Σ_m e^{−(2m+½)x},
  D = ψ(1/4) − ln π = 2θ′(0) = φ(0) = −γ − π/2 − 3 ln 2 − ln π = −5.372183419225666.
Beckner's inequality for even k is "the jump (Dirichlet) part is ≥ 0"; it is sharp because φ − φ(0) = 2∫ κ(x)(1 − cos tx) dx
vanishes only at t = 0 (f ≡ const, k = |y|^{−1/2}, not in L²). Odd k (gamma factor Γ((s+1)/2)) has ψ(3/4) − ln π = −2.2306.
D is exactly the killing constant "log π − ψ(1/4)" of explore17's Markov decomposition of Q − pole (ExploreGFF T2), and
φ(t) = 2θ′(t) < 0 for |t| < 6.29 is the same fact. Locally at p the sharp constant is −2 ln p/(√p − 1) (Burnol F2); these sum to
−∞ over p, so place-by-place Beckner bounds never give a global one.

S1 (N1: window Beckner bound). For supp f ⊂ [−a, a], keeping only the jumps that leave I (κ decreasing, K convex):
  Arch(f) ≥ (D + 2K(a))‖f‖²,   K(a) = ∫_a^∞ κ = Σ_{m≥0} e^{−(2m+½)a}/(2m + ½),   D + 2K(a) = ln(1/(2πa)) − γ + O(a),
positive iff a < 0.085628 (c < 1.18679); odd f with the pole bound −2(sinh a − a): a < 0.085611.
Entropic version (b): Arch/‖f‖² ≥ sup_{β>1} [ln(e/(4a)) − ln ∫ e^{−βφ(2πξ)} dξ]/β (Gibbs), positive iff a < 0.012236 (c < 1.02477).
Lattice version (c) on Chekmachev's identity (f fills half of the period 2L, so H(q) ≥ ln 2):
Σ_k μ(k) q_k ≥ sup_β [ln 2 − ln Σ_k e^{−βμ(k)}]/β = −0.503, −1.401, −2.061 at c = 1.2, 1.5, 2: never positive there.

S2 (N13: Perron–Frobenius, exact range). For c < 2, Q(f) = ∫_I W f² + ½∫∫_{I×I} J(x − y)(f(x) − f(y))² with
J(z) = κ(|z|) − 2cosh(z/2), W(x) = D + ∫_{ℝ∖I} κ(|x − y|) dy + 2∫_I cosh((x − y)/2) dy (pole = 2∫∫ cosh((x − y)/2) f(x) f(y)).
J(z) ≥ 0 ⟺ e^{3|z|} ≤ e^{|z|} + 1 ⟺ e^{|z|} ≤ ρ, ρ = 1.324718 the plastic number (v³ = v + 1). Hence:
  for λ² = c ≤ ρ (λ ≤ 1.150963) the semigroup of CCM's QW_λ is positivity improving (Beurling–Deny: Q(|f|) ≤ Q(f); J > 0 on
  I × I off the diagonal), so by Perron–Frobenius (Reed–Simon XIII.44) and CCM's thmsmallest its smallest eigenvalue is simple
  with a strictly positive, hence even, eigenvector: CCM's first missing step, analytically, on λ ∈ (1, √ρ];
  for c > ρ, f = 1_A − 1_B (small intervals, A − B ⊂ (ln ρ, min(ln c, ln 2)), avoiding every prime distance ln n) gives
  Q(|f|) − Q(f) = −4∫∫_{A×B} J > 0: not positivity preserving; the prime terms Σ_{n<c} w_n(‖f − f(· − ln n)‖² − 2‖f‖²) are
  positive jump atoms and cannot repair a negative jump density. Q − pole is positivity improving for every c (explore17). So
  explore17's "no Perron–Frobenius argument reaches Q itself" holds exactly from λ² = ρ on. Reviewer (verdict only, claims
  inline): S2's five steps CORRECT, one gap (truncate the prime sum to n < c before splitting it), fixed here.

S3 (N13, adversarial content). "The smallest eigenvalue is simple, even, nodeless" holds alike for ζ, L(s, χ₅), DH and moved-pair
models on every positive window tested, and for DH also at a failing window (§4).

## 3. Tests

T1 identities (logs/t1_identities.log, t2_burnol.log). D = 2θ′(0) = φ(0) = closed form to 30 digits; Lévy–Khintchine at
t = 0.5, 1, 5, 14.134725 to 15 digits; killing constant −ln π − γ + ∫(2e^{−x} − 2e^{x/2})/(e^x − e^{−x}) = D. Gaussian e^{−πy²}
(self-dual): −3.8014 = ln(1/2π) − γ − 2 ln 2 > D. Burnol's normalisation, f = (1 − (x/0.5)²)³: Arch = −0.4181825 (x-space)
against ½⟨k, Hk⟩ = −0.4181822 (quadrature of k̂ to |η| = 600), ‖k‖²/‖f‖² = 2.

T2 windows (logs/t3_windows_N40.log, t3b_thresholds.log; CCM basis; Galerkin minima are upper bounds, so exact thresholds are upper
bounds; the odd sine basis vanishes at ±a and converges slowly).

| bound | positive for half-width a < | c < |
|---|---|---|
| global Beckner D | never (−5.372) | — |
| Hirschman–Beckner entropy (b) | 0.012236 | 1.0248 |
| lattice Maassen–Uffink (c) | not at c ≥ 1.2 | — |
| window Beckner S1 | 0.085628 | 1.1868 |
| exact Arch, even | 0.10717 (explore17: 0.107) | 1.2390 |
| exact Arch, odd | 0.3716 (Veltri; ours ≤ 0.3778, N = 80) | ≈ 2.10 |
| Arch + pole (= Q for c ≤ 2), even | ≤ 0.4101 | ≤ 2.2711 |
| Arch + pole, odd | (0.3466, 0.3724] | (2, 2.106] |
| Yoshida / Connes–Consani | 0.3466 | 2 |
| Chekmachev, certified | 1.925 | 47 |

Prime-blind bounds: at c = 2 the Arch + pole margins are 1.330e−3 (even, converged in N) and 0.075 (odd) against the n = 2 atom,
whose norm ln 2/√2 = 0.4901 switches on at c = 2⁺ (Galerkin −0.4901 already at c = 2.1): 368×, so "archimedean minus ‖primes‖"
stops at c = 2 (explore17 found the same). Largest window where a Beckner-type bound beats the prime terms: none; the window bound
is positive only for c < 1.187, where no prime enters.

T3 scale (logs/t3c_scale.log, float64, N = 80; λ_min(Q) from logs/t4_B, c = 47 Chekmachev's ε₄₇).

| c | D + 2K(a) | Arch even / odd | Arch + pole even / odd | primes min even / odd | λ_min(Q) |
|---|---|---|---|---|---|
| 5 | −2.576 | −2.323 / −0.778 | −0.373 / −0.932 | −1.219 / −1.219 | 1.36e−17 |
| 10 | −3.075 | −2.806 / −1.172 | −0.731 / −1.654 | −2.508 / −2.245 | 5.73e−42 |
| 20 | −3.461 | −3.172 / −1.487 | −1.027 / −2.632 | −4.276 / −2.734 | 1.75e−79 |
| 47 | −3.838 | −3.522 / −1.817 | −1.351 / −4.507 | −7.299 / −2.977 | 1.9e−242 |

The archimedean bottom tends to D; separate bounds give ≥ −8.65 at c = 47 against the true +1.9e−242.

T4 ground states (logs/t4_x0.log, t4_A.log, t4_B_zeta_groundstates.log, t4_P_profiles.log).
- Sign change of J at |z| = 0.281200 = ln ρ (S2).
- Q − pole (float64, N = 80): Perron eigenvalue μ₀ = −1.286, −3.334, −5.119, −7.216, −10.525 at c = 2, 5, 10, 20, 47, eigenvector
  nodeless (min/max 0.49–0.84); μ₀ ≈ −2⟨ψ₀, cosh(x/2)⟩² = −1.358, −3.379, −5.141, −7.238, −10.552 (the Perron vector is the pole
  direction); the next even eigenvalues are ≈ 0, as Q ≥ 0 forces (at most one negative eigenvalue).
- Q (mpmath at 4πc log₂e + 128 bits, N = ⌈0.6 c ln c⌉ + 20): even ground state nodeless at c = 2, 3, 5, 10, 20 (f(a/2)/f(0) =
  0.37, 0.13, 0.015, 6.5e−4 at c = 3, 5, 10, 20), odd ground state one-signed on (0, a); overlap with the Perron vector of Q − pole
  0.953, 0.785, 0.657, 0.542, 0.463; spectra alternate even < odd < even (c = 20: 1.75e−79, 9.94e−76, 2.94e−72, 6.44e−69, 1.24e−65,
  1.86e−62). Beyond c = ρ the positivity is a ladder property (explore17: ground state ≈ Φ|_W > 0), not Perron–Frobenius.

## 4. Adversarial

- N1: every bound in S0–S1 depends only on the gamma factor and conductor: the same for every L-function with Γ(s/2) (DH, with
  Γ((s+1)/2) and q = 5, gets ln(5/π) + ψ(3/4) = −0.621 in place of D). No zero information, by construction; window positivity
  below c ≈ 30.7 is DH-blind anyway (explore19).
- N13: L(s, χ₅) and DH have prime weights of both signs (no Markov structure), yet (logs/t4_C_chi5_dh_groundstates.log):

| form | c = 5 even / odd | c = 10 | c = 20 | c = 31 |
|---|---|---|---|---|
| L(s, χ₅) | 6.43e−2 / 0.754 | 1.83e−3 / 0.169 | 8.55e−6 / 1.88e−3 | 4.40e−8 / 2.30e−5 |
| DH | 6.02e−3 / 0.451 | 8.76e−8 / 6.42e−5 | 4.52e−18 / 1.42e−14 | −6.08e−30 / 7.46e−26 |

  all ground states even and nodeless; DH at c = 31 (its form fails, explore19) still has a simple even nodeless ground state.
- Moved pair (explore19's model: γ_k, γ_{k+1} → m ± iδ; code/t5_movedpair.py, logs/t5_movedpair_parity.log): 14 models,
  k ∈ {1, 2, 3, 5, 10, 20}, δ ∈ {0.003, 0.03, 0.3}, first failure at c ∈ [2.03, 7.21] (resolution 1/256, inside explore19's
  brackets): the even sector fails first in all 14 and the ground state is even at all 6 positive windows below each threshold
  (84/84). Parity flips only after failure (k = 5, δ = 0.3, c = 3.53: odd −1.42e−8 < even −3.04e−10, explore19's log). So a
  simple even ground state never flags an off-line pair on a positive window: CCM's step one is RH-neutral; their RH content sits
  in positivity and in the second step.

## 5. Verdict and novelty

- S0 is a re-encoding: Burnol (1998–99) plus Beckner (1995); the Lévy–Khintchine/Markov form is explore17's T2.
- S1: new numbers, weaker than known windows; a relaxation of the exact compressed form cannot beat it, and the exact form already
  stops at 0.107 (even) without the pole, at c = 2 with it (prime-blind).
- S2: corrects explore17's remark; the range λ ≤ 1.151 is inside Kramarenko-Byrd's L ≤ 4/5 and Botha's L ≤ 2 (computer-assisted);
  an analytic proof on a tiny range, not a discovery (no reviewer run; arXiv search found no Beckner–Weil statement).
- S3 extends explore17's DH-blindness of step one to a failing DH window and to moved-pair models.
- No node quantity moves.

## Suggested NOGO line

- explore32 (Orca, exploration: entropic uncertainty — Beckner's logarithmic uncertainty principle and Pitt's inequality, Hirschman–Beckner, Maassen–Uffink → N1, N13) NO IDEA: Beckner's sharp constant ψ(1/4) − ln π = 2θ′(0) = −5.3722 is the bottom of Weil's archimedean symbol (Burnol's conductor operator: Arch(f) = ½⟨k, (ln|y| + ln|η|)k⟩, k = |y|^{−1/2} f(ln|y|)), i.e. the killing constant of explore17's Markov form; on [−a, a] it gives Arch ≥ (D + 2Σ e^{−(2m+½)a}/(2m+½))‖f‖² > 0 only for a < 0.0856 (Hirschman–Beckner entropy a < 0.0122, lattice Maassen–Uffink never) against the exact 0.107 (even) and Yoshida's (ln 2)/2; the c = 2 margin 1.3e−3 against w₂ = 0.49 stops every prime-blind bound; Perron–Frobenius proves CCM's step one exactly for λ² ≤ ρ = 1.3247 (plastic number: κ(z) ≥ 2cosh(z/2) ⟺ e^{3z} ≤ e^z + 1); the simple even nodeless ground state also holds for L(s, χ₅), for DH at its failing window c = 31 (λ_min = −6.1e−30) and for 14 moved-pair models on every positive window (the even sector always fails first), so step one carries no RH information [ExploreEntropy](Riemann/Deep/ExploreEntropy/REPORT.md)

RH relevance: N1 (archimedean window bounds), N13 (CCM's first missing step).

## Files

code/: t1_identities.py (S0 constants, Lévy–Khintchine), t2_burnol.py (Burnol normalisation), t3_windows.py, t3b_thresholds.py,
t3c_scale.py (S1, T2, T3), t4_pf.py x0|A|B|C|P (S2, T4, §4), t5_movedpair.py + t5_run.sh (moved pairs), run_quick.sh (all fast logs);
copied unchanged: weil_matrix.py (../WeilSimple), wmat.py and movedpair.py (../ExploreAntenna/code). Every process single-threaded,
peak RSS ≤ 148 MB (/usr/bin/time -l in each log).
