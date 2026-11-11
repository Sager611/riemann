# STATEMENTS: Connes–Moscovici counting law vs Bhattacharjee (node N13)

Sources (local copies in .scratch/prolatefals/, deleted after the round):
- **CM21** = arXiv 2112.05500 v1, the only arXiv version (10 Dec 2021), "Prolate spheroidal operator and Zeta"; e-print TeX `Draft2.tex` (1059 lines). "l.N" = line N of that file.
- **CM22** = PNAS 119(22) e2123174119 (2022), "The UV prolate spectrum matches the zeros of zeta"; JATS/MathML full text from Europe PMC (PMC9295779, fullTextXML). Equation numbers as printed.
- **B26** = D. Bhattacharjee, Zenodo 24127214 v1.0.0 (11 Nov 2026), "The Negative Prolate Spectrum and the Riemann–Siegel Theta Function"; `paper/main.tex` (443 lines). "b.N" = line N.

## Common definitions (all three agree verbatim)

| Object | Definition | CM21 | CM22 | B26 |
|---|---|---|---|---|
| Operator | (W_λξ)(x) = −∂_x((λ²−x²)∂_xξ) + (2πλ)²x²ξ on ℝ | l.260–262 | Eq. 1, 3 | b.56–58 |
| Fourier | F(f)(y) = ∫ f(x) e^{−2πixy} dx | l.307–309 | Eq. before 1 | — |
| W_sa | W_max restricted to L_β: lim_{x→±λ}(λ²−x²)ξ′ = 0; at ±∞ even part: lim(x sin(2πλx)ξ⁺′ − (2πλx cos(2πλx) − sin(2πλx))ξ⁺) = 0, odd part: lim(x cos(2πλx)ξ⁻′ + (2πλx sin(2πλx) + cos(2πλx))ξ⁻) = 0 | l.567–583 | Eqs. 17–19 | b.122–131 |
| Properties | self-adjoint, commutes with F, P_λ, P̂_λ, unique such extension, discrete spectrum unbounded both ways | Thm (label thmwsa) l.590–596 | §1 | b.53 |
| Eigenfunctions | regular at λ from each side; at ∞ ∝ sin(2πλx)/x (even), cos(2πλx)/x (odd) | Cor l.638–646 | — | Lemma 2.1 b.148 |

Negative eigenvalues μ = −E² live on |x| > λ (eigenfunctions vanish on (−λ,λ): b.113–118; Sonin space: CM21 Cor l.720–726). On (λ,∞): ((x²−λ²)u′)′ + (4π²λ²x² + E²)u = 0, u regular at λ, plus the parity condition at ∞ (b.140–146).

Semiclassical count (not a theorem about eigenvalues), CM21 Prop 4.2 l.788–792 = CM22 Prop 3.2 Eq. 24: number of negative eigenvalues ξ with −ξ ≤ E², on even (equally odd) functions, ≈ 2σ(E,λ),
σ(E,λ) ~ (E/2π)(log(E/2π) − 1 + log 4 − 2 log λ) + λ² + o(1), i.e. 2σ = (E/π)(log(2E/(πλ²)) − 1) + 2λ² + o(1).

Dirac operator: D = [[0, ∇+w],[∇−w, 0]], ∇ = p^{1/2}∂ (CM21) resp. p^{1/4}∂p^{1/4} (CM22), p = x²−λ², w any solution of the Riccati equation (CM21 l.838–840, Prop 5.3 l.879–899; CM22 Lemma 4.1, Prop 4.3). CM22 Prop 4.3(ii): Spec D_w = {±√μ : μ ∈ Spec W′⁺_sa}, W′⁺_sa = (I−P_λ)W_sa on even functions ≅ W on L²(λ,∞) with the even conditions; independent of w.

## The two counting laws, side by side (λ = √2, operator 2D)

N(E) := #{eigenvalues of 2D with 0 < Im ≤ E} (CM: "less than E"; the difference is O(1)). Imaginary eigenvalues of 2D are ±2i√(−μ), μ < 0 in Spec W′⁺_sa, so with −(E_n⁺)² the even negative eigenvalues (B26 notation) **N(E) = N₊(E/2)**, t_n = 2E_n⁺.

| | Statement | Proof given |
|---|---|---|
| CM21 Thm 6.1 (l.904–908, eq. ne) and intro (l.271–273) | N(E) ~ (E/2π)(log(E/2π) − 1) + O(1) | N(E) = 2σ(E/2, √2) from the semiclassical Prop 4.2 (l.909–920): a phase-space area, not an eigenvalue count |
| CM22 intro Eq. 2 | N(E) ≈ (E/2π)(log(E/2π) − 1) + O(log E) | — (compatible with both sides) |
| **CM22 Thm 5.1, Eq. 33** | **N(E) ≈ (E/2π)(log(E/2π) − 1) − (log E)/(2π) + O(1)** | Heywood's formula (1.4) of ref. 19 (Heywood 1954) and extensions (refs. 20, 21) after a Liouville transformation W′⁺_sa ≅ S on (0,∞), Sφ = φ″ − qφ, q(y) = −(4π)²cosh²y − ¼(coth²y − 2) (Eq. 34); "a delicate computation which is a refinement of the semiclassical estimate of Proposition 3.2 gives the additional logarithmic contribution −(log E/2π)". No further detail printed. |
| **B26 Thm B** (b.104–113) | ϑ(t_n) = (n + n₀ − 3/8)π + O(t_n^{−2/3}) for some integer n₀; N(E) = ϑ(E)/π + O(1) = (E/2π)(log(E/2π) − 1) + O(1); N(T) − N_ζ(T) = −S(T) + O(1) | Thm A with κ = 4π, E = t/2 (b.377–391) |
| **B26 Thm A** (b.75–86), any λ > 0, κ = 2πλ² | Θ continuous, strictly increasing, Θ(E) = Φ(E) + O(E^{−2/3}), Φ(E) = E log(4E/κ) − E + π/4; −E² ∈ Spec W⁺_sa ⇔ Θ(E) ∈ πℤ, ∈ Spec W⁻_sa ⇔ Θ(E) ∈ π/2 + πℤ; all simple; N_±(E) = (E/π)(log(2E/(πλ²)) − 1) + O(1); E_n log(2E_n/(πλ²)) − E_n = (n + n_± ∓ ¼)π + O(E_n^{−2/3}) | Prüfer angle in x = λ cosh t (b.163–186), Θ(E) = lim(φ − κ cosh t) (b.236–238); criterion Prop 2.5 (b.249); monotone Prop 2.6 (b.261); Bessel J₀ comparison at t = 0 (Lemmas 3.2–3.3, b.307–345); tail by parts (b.347–357); Thm 3.5 (b.359) |

## Same eigenvalues, same normalisation?

Yes. Both count the negative eigenvalues of the even part of W_sa for λ = √2 on (√2, ∞) with the boundary conditions above, through t = 2E (CM22 Prop 4.3(ii) and proof of Thm 5.1; B26 b.377–384). The only change of variable, t = 2E, is linear: it shifts log by a constant and cannot create or remove a (log E)/2π term. The variable change x = λ cosh t (B26) and the Liouville transformation (CM22 Eq. 34) act on x, not on the spectral parameter. So there is no reconciliation: N_CM22(E) − N_B26(E) = −(log E)/2π + O(1) → −∞. Exactly one of CM22 Eq. 33 and CM21 eq. ne (= B26) can hold.

Main terms agree: B26's (E/π)(log(2E/(πλ²)) − 1) equals CM's 2σ(E,λ) minus the constant 2λ²; for λ = √2, t = 2E both give (t/2π)(log(t/2π) − 1).

Also TWW26 (arXiv 2609.32639, l.211–231, 239–246; operator P = −W_1, eigenvalue variable E_TWW = E²) count the same spectrum (both parities): Bohr–Sommerfeld 2I(1 + E_n/4π²) + ¼ + R(E_n) = n ∈ ½ℤ, R(E) = −(1/16π)E^{−1/2}log E + O(E^{−1/2}); Weyl law N₊(E) = 4(√E/2π)(log(√E/2π) − 1 + 2 log 2) + O(1); general μ = (2πλ²)² by Rem. 1.4. Translated: N₊ per parity = (E/π)(log(2E/(πλ²)) − 1) + O(1) = B26 Thm A.

## What is tested (Riemann/Deep/ProlateFals/REPORT.md)

1. log coefficient b in N(t) − (t/2π)(log(t/2π) − 1) = c + b log t: CM22 ⇒ b = −1/(2π) = −0.159155; CM21/B26 ⇒ b = 0.
2. B26 location: ϑ(t_n)/π − n → n₀ − 3/8; exponent of the error.
3. B26 Thm A for other λ and the odd parity: E log(4E/κ) − E + π/4 − Θ(E) → 0.
