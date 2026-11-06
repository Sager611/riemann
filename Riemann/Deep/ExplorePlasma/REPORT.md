NO IDEA: plasma stability theory proves zero-free half-planes only for one-sided (Cauchy/Laplace) transforms of single-humped or monotone profiles; for Xi the step to the two-sided transform is Hermite–Biehler for Xi+(−z), which fails (Xi+ has zeros 13.788 − 7.411i, 20.744 − 11.828i, …), and every Gardner/Penrose hypothesis also holds for the Davenport–Heilbronn kernel, a moved-pair xi-like kernel and the de Bruijn–Newman kernels at t < 0, all of which have off-line zeros.

# ExplorePlasma: plasma physics and kinetic theory → RH tree (explore20, Orca, 6 Nov 2026)

Normalisation (code/kern.py): Ξ(z) = ξ(1/2 + iz) = 2∫_ℝ Φ(u)e^{izu}du, Φ(u) = Σ_n (2π²n⁴e^{9u/2} − 3πn²e^{5u/2}) exp(−πn²e^{2u}).
Polymath's kernel is Φ_P(u) = Φ(2u), so Φ_P(u)e^{tu²} ↔ Φ(v)e^{sv²} with t = 4s; every t below is Polymath's (Λ ≤ 0.10756 in TREE N2).
One-sided transform Ξ⁺(z) = 2∫_0^∞ Φ(u)e^{izu}du: Ξ(z) = Ξ⁺(z) + Ξ⁺(−z), and on ℝ, Ξ⁺(x) = Ξ(x)/2 + iS(x) with S(x) = 2∫_0^∞ Φ(u) sin(xu)du.
Plasma: 1D electrostatic, fixed ions, f0 normalised, f̃0(s) = ∫f0(v)e^{−isv}dv; dielectric function
ε(k, kz) = 1 − k⁻²∫f0′(v)/(v − z)dv = 1 + k⁻²∫_0^∞ s f̃0(s) e^{izs}ds (Im z > 0).

## Dictionary
| plasma / kinetic theory | Riemann side |
|---|---|
| equilibrium f0 > 0, single-humped | Φ: positive, even, decreasing on (0, ∞), log-concave ((log Φ)″ ≤ −18.73 on (0, 3], rechecked) |
| characteristic function f̃0; free-streaming density ρ(t) = f̃0(kt) (phase mixing) | Ξ(s)/Ξ(0): the density of the Φ-plasma vanishes exactly at t = γ_n/k |
| dielectric function of the Φ-plasma: 1 + (k²Ξ(0))⁻¹∫_0^∞ sΞ(s)e^{izs}ds | the zeros of Ξ sit inside the integrand, not in the root set |
| Landau continuation into Im z < 0 | stops at Im z = −π/4: Φ is analytic only for abs(Im u) < π/4 (natural boundary of θ); decay rate π/4 = Ξ(t) ~ t^{7/4}e^{−πt/4} |
| Penrose criterion (winding of ε(ℝ)) | argument principle (R404) |
| Gardner: single hump ⇒ no growing mode for any k | monotone kernel ⇒ zero-free half-plane for a one-sided transform (Eneström–Kakeya; T1) |
| Kruskal–Oberman energy ∫abs(h)²/(−f0′/v) + field energy > 0 ⇒ real frequencies | Hilbert–Pólya (R51, R62); Weil's form is its analogue with an indefinite prime term |
| Case–van Kampen modes; Krein signature of the continuum (Hagstrom–Morrison 1002.1039) | Pontryagin index = number of off-line zeros (explore18) |
| Sturrock's wave energy ω ∂D/∂ω of a neutral mode | sgn γ_nΞ′(γ_n) alternates: Ξ is maximally Krein-indefinite (T3) |
| two-stream (double-humped) instability | Pólya's increasing-kernel theorem; de Bruijn's large-t limit (mass at ±u*) (T2) |
| Fokker–Planck (Ornstein–Uhlenbeck) smoothing of f0 | Ξ(e^{−τ}z)·exp(−(1 − e^{−2τ})z²/2): zeros only dilated; the de Bruijn–Newman flow is the Gibbs tilt Φe^{tu²} ("heating"), not collisions |
| Landau roots as zero-collision limits of true eigenvalues (2109.07806) | ζ(2s) zeros as Lax–Phillips/Pavlov–Faddeev resonances (R5.5, R58); viscosity limits (1912.09840) are non-self-adjoint, no line constraint |
| Penrose in a strip ⇒ exponential decay (Volterra / Paley–Wiener; Mouhot–Villani) | Möbius renewal: ζ ≠ 0 on Re s > σ ⇔ M(x) = O(x^{σ+ε}) (N10; renewal equations R378) |
| plasma echo: phase-mixed modes re-cohere | Landau's 1911 formula Σ_{0<γ≤T} x^ρ = −(T/2π)Λ(x) + O_x(log T): the zeros re-cohere at prime powers (explicit formula) |
| analytic / Gevrey "gliding" norms (Mouhot–Villani, Bedrossian–Masmoudi–Mouhot) | no node lacks them: the "equilibrium" Σ_n δ_{log n} has no analyticity strip; continuing 1/ζ is the problem itself |

## Survey (arXiv API, abs/e-print pages, monthly listings)
- Linear theory: Landau 1946; van Kampen 1955 and Case 1959 (singular modes); Penrose, Phys. Fluids 3 (1960) 258 (unstable iff some local minimum v0 of f0 has
  ∫(f0(v) − f0(v0))/(v − v0)² dv > 0, proved by a Nyquist winding); Gardner, Phys. Fluids 6 (1963) 839 (f0 decreasing in energy ⇒ stable, energy proof);
  Hamiltonian/Krein view: Morrison's G-transform, Hagstrom–Morrison 1002.1039 (Krein-like theorems for continuous spectra: instability needs a signature
  change of the continuum, i.e. a local minimum of f0).
- Nonlinear theory: Mouhot–Villani 0904.2760 (analytic/Gevrey data, Penrose condition in a strip, echoes controlled by a gliding analytic radius);
  Bedrossian–Masmoudi–Mouhot 1311.2870; Grenier–Nguyen–Rodnianski 2004.05979; Bedrossian 1605.06841 (echo chains defeat damping in finite regularity);
  Lin–Zeng 1003.3005, 1106.4368 (small BGK waves in H^s, s < 3/2). 2026: Lean formalisation of Mouhot–Villani for Gevrey index s > 1/3 (2609.16801),
  pointwise damping around compactly supported radial equilibria (2610.26600), BGK stability (2606.16194, 2412.07025), 2D unconfined damping (2609.28251).
- Collisions: Villani, Hypocoercivity (math/0609050); Bedrossian 1704.00425 (weak collisions suppress echoes); Landau modes become true eigenmodes in the
  zero-collision limit (2109.07806; Ng–Bhattacharjee–Skiff 1999); resonances as viscosity limits (Galkowski–Zworski 1912.09840).
- Links to ζ/RH: API searches ("Landau damping" with Riemann or "zeta function"; Vlasov with "Riemann hypothesis"; "Penrose criterion" with zeta;
  Penrose with Hermite-Biehler or Nyquist; "van Kampen" with Krein; "plasma dispersion function" with zeros; "Riemann zeros" with plasma; "Landau damping"
  with de Bruijn, Polya, Stieltjes, Hilbert transform) return nothing relevant. Listings 2026-10/11 of math.NT, math.CV, math.CA, math.SP (plasma words)
  and math-ph, physics.plasm-ph 2026-10 (zeta words): only Landau-level items (2610.14108, the Sierra–Townsend line, not Landau damping) and known
  de Bruijn–Newman papers; math.AP and physics.plasm-ph 2026-11 listings did not load (one retry each). Our files: R404's Penrose line, R369's kinetic line,
  explore18's Krein line, R378's renewal line, R91's Eneström–Kakeya row ("not a polynomial").
- Zeta-side classics used: Pólya 1918 (increasing kernels on [0, 1]), 1926 (K_{iz}(a)), 1927 (e^{−u^{2n}}); de Bruijn 1950; Csordas–Norfolk–Varga 1986
  (Φ decreasing); Newman 1976; Rodgers–Tao 1801.05914 (Λ ≥ 0); Lagarias's Herglotz criterion (R399); E. Landau 1911.

## Transfers, tests, verdicts
Code in code/, logs in logs/ (each run single-threaded, peak memory 40–130 MB).

**T1. Gardner/Penrose = monotonicity ⇒ zero-free half-plane for a one-sided transform; the Hermite–Biehler gap (N2, N13).**
Gardner's proof: Im of the boundary value of a Cauchy transform is local (Im C_g(x + i0) = πg(x)); with one sign change (single hump) and Re C < 0 at the
crossing, the boundary curve never meets (0, ∞). For f0 = Φ: (k²Ξ(0))⁻¹∫_0^∞ sΞ(s)e^{izs}ds ≠ −1 on Im z > 0 for every k. Its Fourier twin is the
continuous Eneström–Kakeya theorem: K > 0 decreasing on (0, ∞) ⇒ ∫_0^∞Ke^{izu}du = (i/z)(K(0) + ∫_0^∞K′e^{izu}du) ≠ 0 on Im z > 0, since
abs(∫K′e^{izu}) < ∫abs(K′) = K(0). For K = Φ: Ξ⁺ ≠ 0 on Im z ≥ 0 (grid minimum of abs(Ξ⁺) on [−60, 60] × [0, 20]: 0.0142). This is also the plasma whose f0′ is Ξ:
C_Ξ = 2πiΞ⁺. To reach Ξ itself, write Ξ = 2A with E(z) = Ξ⁺(−z) = A − iB, A = Ξ/2, B = S: if E were Hermite–Biehler, Ξ would have only real zeros.
Test (logs/onesided.log): E is not Hermite–Biehler, in two ways.
(a) Ξ⁺ has zeros in the lower half-plane, so E has zeros in the upper one: 13.78823935 − 7.41080193i (mpmath, 30 digits: abs(Ξ⁺) = 4.7·10⁻³⁷),
20.74424832 − 11.82771449i, 26.83506285 − 15.22958302i, 32.4352 − 18.1320i, 37.7045 − 20.7241i, 42.7294 − 23.0993i, 47.5633 − 25.3121i, 52.2421 − 27.3971i,
56.7912 − 29.3782i (all with 0 < Re < 100, Im > −30), on a curve with Im/Re falling from −0.57 to −0.52.
(b) Even next to the axis: abs(E)² − abs(E#)² ≈ 4y(AB′ − A′B), and at each zero AB′ − A′B = −Ξ′(γ_n)S(γ_n)/2 with S > 0 on (0, ∞) (Φ decreasing; checked
on (0, 60]), so the sign alternates and the inequality fails beside γ2, γ4, γ6, ….
Same for Davenport–Heilbronn (Φ_DH decreasing): Ξ⁺_DH ≠ 0 on the upper grid (minimum 0.0182) and has zeros 7.45665419 − 2.59532437i, 12.23283359 − 4.74352465i,
16.4332 − 6.4509i, … (8 with Re < 80, Im > −14). Verdict: dead. Gardner/Eneström–Kakeya is true for ζ, DH and every decreasing kernel; the missing step is
Hermite–Biehler (R34, R404), false for this E.

**T2. Single hump versus two streams along the de Bruijn–Newman family (N2).** The condition evaluated on Φ (logs/shape.log).
Φ_P e^{tu²} is single-humped (Gardner-stable for all k) iff t < t_b = −(log Φ_P)″(0)/2 = 37.4538; it has a dip at 0 for t > t_b, and the Penrose integral at
the centre, P(t) = ∫(Φ_P(u)e^{tu²} − Φ_P(0))/u² du, satisfies P′(t) = ∫Φ_P e^{tu²}du = 2H_t(0) > 0 (identity checked numerically at t = 0), so P changes sign exactly once, at
t_P = 47.7208 (P(0) = −9.835, P′(0) = 0.12428): for every t > t_P the centre is a local minimum with positive Penrose integral, so the plasma is unstable. Real-rootedness holds exactly for t ≥ Λ ∈ [0, 0.10756].
Davenport–Heilbronn: t_b = 5.1073, t_P = 10.0928, with off-line zeros at t = 0. Moved-pair kernel (below): t_b = 40.76–40.81.
Reading: plasma stability and real zeros are independent: both hold on [Λ, 37.45); stability with non-real zeros for t < 0 (Rodgers–Tao; explicit pair at
t = −2 below) and for DH at t = 0; real zeros with instability for t > 47.72. Pólya's increasing-kernel theorem and de Bruijn's large-t mechanism live on
the two-stream side (mass near ±u*, transform ≈ cosine), Gardner's on the single-hump side; Φ is single-humped by a margin of 37.45 in t.
Barrier (Gaussian damping): every Gardner/Penrose-type hypothesis (positivity, evenness, single hump, log-concavity of K or of K(√u), Penrose's sign)
survives K ↦ Ke^{−au²}, a > 0. H_t has non-real zeros for every t < 0 (Λ ≥ 0), so none of them implies real zeros. Explicit (logs/negt.log): at t = −1.6 the
zeros γ9, γ10 are still real (z_P ≈ 98.30, 99.90); at t = −2 they are the pair z_P = 99.43089514 ± 0.49715483i (abs(H) = 1.5·10⁻²⁹ against 5·10⁻¹⁶ on the line);
at t = −2.4, 99.76784 ± 1.04881i. Verdict: dead; the barrier is folklore (Newman's "RH is barely true") made explicit.

**T3. Case–van Kampen, embedded modes, Krein signature (N13).** The dual plasma f0′ = Ξ (f0 = ∫_{−∞}^v Ξ, a beam edge): its dielectric function
1 − (2πi/k²)Ξ⁺(z) is entire; its Landau roots tend to the zeros of Ξ⁺ (T1) as k → 0. Im ε is local, so its neutral (embedded) modes sit exactly at
phase velocities −γ_n, with k_n² = 2πS(γ_n) = 0.45644, 0.27965, 0.23173, 0.18844, 0.17352, 0.15140, 0.13876, 0.13088, 0.11790, 0.11365 (n = 1–10,
logs/dual.log): the real zeros become embedded eigenvalues of a Vlasov operator, one wavenumber each. Off-line zeros are invisible: for the moved-pair Ξ_off
(δ = 0.25) the modes at −γ1, −γ2 vanish, and the Penrose winding (number of growing modes) at k² = 0.30, 0.25, 0.20, 0.17, 0.15, 0.13, 0.12, 0.11 is
1, 0, 1, 1, 0, 0, 0, 0 for Ξ and 0, 1, 0, 1, 1, 1, 0, 0 for Ξ_off: the parity of #{n: k_n² > k²}, while Ξ has 0 and Ξ_off 2 zeros in the upper half-plane.
Krein/Sturrock: a simple real root ω0 of a real dispersion function D carries wave energy ∝ ω0 D′(ω0); for D = Ξ the signs alternate, so every adjacent pair
may collide and leave the axis (it happens at t = −2 above). A Krein-definite (Gardner-type) system needs poles between roots, i.e. a Herglotz function
such as −Ξ′/Ξ, whose Herglotz property is RH (Lagarias; R399, explore18). Weil's form (N1) is the Kruskal–Oberman analogue only up to its prime term, whose
sign is indefinite; Penrose's lesson (indefinite energy, still stable, decided by winding) maps to "Weil positivity is decided zero by zero", i.e. RH.
Verdict: dead.

**T4. Damping rates, echoes, Gevrey norms (N7, N10).** No numerics. Penrose in a strip ⇒ exponential decay is the Paley–Wiener theorem for Volterra
equations; for the Möbius renewal (1/ζ as the resolvent of the measure −Σ_{n≥2}δ_{log n}) it is the classical equivalence ζ ≠ 0 on Re s > σ ⇔
M(x) = O(x^{σ+ε}). Mouhot–Villani's input (analytic f0 ⇒ the dielectric continues into a strip) has no counterpart: Σ_n δ_{log n} has no analyticity
strip, and continuing 1/ζ past Re s = 1 is the unknown. Echoes are Landau's 1911 formula. Verdict: dead (known equivalences).

## Adversarial tests
| kernel | off-line zeros | positive, even, decreasing on (0, ∞) | log-concave | Gardner/Penrose stable for all k |
|---|---|---|---|---|
| Φ (ζ) | none known (RH to 3·10¹²) | yes | (log Φ)″ ≤ −18.73 | yes (t_b = 37.45) |
| Φ_DH = e^{3u/2}Σ a_n n e^{−πn²e^{2u}/5} | 0.808517 + 85.699348i (abs Λ_DH = 6.7·10⁻³⁵ against 1.9·10⁻²⁹ on the line) | yes | (log Φ_DH)″ ≤ −2.554 | yes (t_b = 5.11) |
| Φ_off, δ = 0.05, 0.25, 0.45 | 17.578382 ± iδ | yes | ≤ −20.38 | yes (t_b = 40.8) |
| Φ_P e^{−2u²} (H_t, t = −2) | 99.43090 ± 0.49715i | yes | yes | yes |
| e^{−u²}(1 + u⁴/4) (R38, Lean) | yes | yes | yes | yes |
- DH kernel (logs/norm.log): a_n = 1, κ, −κ, −1, 0 (n mod 5), κ = 0.284079043840412; evenness (the functional equation, root number 1) holds to 15 digits
  at u = ±0.05, ±0.2, ±0.4; 2∫Φ_DH = Λ_DH(1/2) = 1.43319216877501 and the transform at z = 10 matches. Shapes checked on (0, 3] with exact derivatives
  (grid 0.005), beyond which the n = 1 term dominates; also log-concave in √u (logs/shape.log).
- Moved pair (code/offline_np.py, logs/offline.log): Ξ_off = Ξ·N(z²)/((z² − γ1²)(z² − γ2²)), N(T) = (T − γ_m² + δ²)² + 4γ_m²δ², γ_m = (γ1 + γ2)/2; kernel
  Φ_off = Φ + Aψ1 + Bψ2 with ψ_j(u) = −γ_j⁻¹∫_u^∞ sin(γ_j(w − u))Φ(w)dw (the decaying solution of −ψ″ − γ_j²ψ = Φ, even because Ξ(γ_j) = 0); its transform
  equals Ξ·R to 10⁻¹⁶ at z = 3, 8, 12, 25 and vanishes (10⁻¹⁶, double precision) at γ_m + iδ; shapes checked on (0, 2.4], and on (2.4, 3.2]
  abs(Φ_off/Φ − 1) ≤ 4.2·10⁻⁵.
- Does a transfer "prove" real zeros for DH or the moved pair? T1 and T2 apply to all five kernels alike and prove only one-sided zero-freeness; T3 sees
  only real zeros. None separates ζ from the counterexamples, as required, because none says anything about Ξ's zeros off the line.
- Restatement: T1 = Eneström–Kakeya (R91) + Hermite–Biehler (R34, R404); the Penrose criterion = argument principle (R404); T2's real-zero side = Pólya and
  de Bruijn (N2 line), its barrier = Λ ≥ 0 (R71); T3 = Herglotz (R399), Krein (explore18), Lax–Phillips (R5.5, R58); T4 = N10 equivalences and renewal (R378).

## Ranking (RH relevance × chance × gain)
| transfer | node | RH relevance | chance | gain | verdict |
|---|---|---|---|---|---|
| T1 Gardner/Eneström–Kakeya, Hermite–Biehler gap | N2, N13 | would be RH | ~0: E = Ξ⁺(−z) is not Hermite–Biehler (computed zeros, alternating AB′ − A′B) | none | dead |
| T2 single hump vs two streams, Gaussian-damping barrier | N2 | barrier only | 0 as a route | a clean filter for kernel criteria | dead (barrier) |
| T3 dual plasma, embedded modes, Krein/Sturrock | N13, N1 | low | 0: sees real zeros only; Krein-definite ⇔ Herglotz ⇔ RH | none | dead |
| T4 damping rates, echoes, Gevrey norms | N7, N10 | low | 0 | none | dead (known equivalences) |

## Suggested NOGO line
- explore20 (Orca, exploration: plasma physics and kinetic theory, Landau damping, Penrose/Gardner, Case–van Kampen/Krein, echoes, Fokker–Planck → N2, N13, N7, N10, N1) NO IDEA: Gardner's single-hump stability, read in Fourier variables, gives zero-free half-planes only for one-sided transforms (Eneström–Kakeya: Ξ⁺(z) = 2∫_0^∞Φe^{izu} ≠ 0 on Im z ≥ 0); reaching Ξ = Ξ⁺(z) + Ξ⁺(−z) needs Hermite–Biehler for Ξ⁺(−z), false (Ξ⁺ has zeros 13.788 − 7.411i, 20.744 − 11.828i, …; AB′ − A′B = −Ξ′(γ_n)S(γ_n)/2 alternates); every Gardner/Penrose hypothesis survives Gaussian damping, so Λ ≥ 0 blocks it (H_{−2} has the pair 99.431 ± 0.497i), and the Davenport–Heilbronn and moved-pair kernels are single-humped and log-concave too; Φe^{tu²} turns two-stream only at t_b = 37.45 (Penrose-unstable at 47.72) ≫ Λ; the dual plasma f0′ = Ξ has embedded modes exactly at the real zeros (k_n² = 2πS(γ_n)) and never sees off-line ones [ExplorePlasma](Riemann/Deep/ExplorePlasma/REPORT.md)

## Files
- code/kern.py: Φ and Φ_DH with closed-form derivatives, ξ and Λ_DH (mpmath). code/offline_np.py: the moved-pair kernel Φ_off (numpy).
  code/onesided.py: Ξ⁺ for Φ and Φ_DH, zero finder (chunked; about 75 MB).
- code/run_norm.py, run_shape.py, run_offline.py, run_onesided.py, run_dual.py, run_negt.py → logs/norm.log, shape.log, offline.log, onesided.log,
  dual.log, negt.log (rerun: `uv run --no-project --with mpmath --with numpy python3 code/run_<name>.py`, all six in about 4.5 min).

RH relevance: N2 (Gardner/Penrose kernel criteria on the de Bruijn–Newman family; barrier only), with N13 (Case–van Kampen, Krein) and N7/N10 (damping ↔ zero-free strips) secondary; no node advanced.
