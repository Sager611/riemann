# DPC 2610.28527 — re-derivations (audit, 12 Nov 2026)

TeX lines refer to v2 `Gaps_for_Dirichlet.tex`. Scripts in this folder: `opcheck.py` (operator calculus), `largegap.py` (large-gap quadratic forms), `cert.py` (Arb certificate; reads the `Jsym` file that `largegap.py 1000 12 5` writes to `.scratch/dpcaudit/`), `models.py` (model tests). Run with `uv run --no-project --with numpy --with scipy --with python-flint --with mpmath python3 <script>`.

## D1. Reduction of A (L57, L580–583)

If ρ = β+iγ, β ≠ ½, is a zero then so is 1−ρ̄ (ξ_χ(s) = conj ξ_χ(1−s̄), L156): two zeros with ordinate γ, hence γ_{n+1} = γ_n. Infinitely many such ρ at positive heights ⇒ μ = 0. Otherwise all zeros above some Y₀ lie on the line (and, if μ > 0, are simple). Correct, but it rests on the convention "ordinates counted with multiplicity" (L47, L68). Under the convention of the unconditional critical-line literature (Hall, Bredberg, Bui–Milinovich: gaps between consecutive zeros *on the line*) Theorem 1 does not follow. Substantive content of A: μ ≤ ½ whenever all high zeros lie on the line.

## D2. V counts zeros; link to Bondarenko–Heap

Boundary values from Re s > ½: near a zero ½+iγ, arg(ε+i(t−γ)) rises by π; between zeros ξ(½+it) is real, so Im log ξ is constant. Hence V(t) ∈ ℝ + iℤ_{≥0} for a.e. real t, always, with Im V = N(t) = #zeros in (t−h/2, t+h/2). Z := e^{πV} = ξ(½+i(t+h/2))/ξ(½+i(t−h/2)) is real with sign (−1)^N: the quotient of shifted zeta functions used by Bondarenko–Heap 2610.26038.

## D3. Prop 1 (L262–388): every estimate re-derived

- Index sets: integers q ≤ N or reduced m/n, m,n ≤ Y, always (q,Q) = 1 (L241). For q ≠ q′ and integer k: the reduced denominator of kq/q′ is ≤ B = T/L²⁰, so |ω| = |log(kq/q′)| ≥ 1/2B (log(a/b) ≥ log(1+1/b) ≥ 1/2b for a > b).
- Contour: integrand gV^jR(t)R*(t) on [a,b]×[−1,0], a ≈ T/2, b ≈ 5T/2 chosen ≥ T⁻² from shifted zeros; V analytic (½ < Re s < 3/2, Im s > Y₀); boundary limit by the ε-uniform majorant of Lemma 1. Sides: |g| = exp(−(T²−y²)L²/2T²) ≤ e^{−L²/4}, |RR*| ≤ B²‖c‖², |log ξ| ≪ T log²T: negligible. [T/2,5T/2]∖[T,2T]: g ≤ e^{−L²/8}, Lemma 1 bound: negligible.
- On Im t = −1: log L series absolutely convergent, Σ_k|β_r(k)|/k ≤ (C/L)^r. Each term ∫gG^m e^{−itω}: shift only G and the exponential (not V) to Im t = −sgn(ω)H; G is analytic on Re t ∈ [a,b], |Im t| ≤ H (gamma factor has no zeros/poles there), G = iφ + O(1/L) by Stirling (|s| ≍ T, |arg s/2| ≤ 3π/4). |∫| ≪ H(e^{−H|ω|} + e^{−L²/4}e^{−ω}), e^{−H|ω|} ≤ e^{−L¹⁹/4}q′/(kq) since H/2B = L¹⁹/2. With Σ|c_qc_{q′}|q′/q ≤ B²‖c‖²: total o(H‖c‖²).
- Diagonal kq = q′: χ(q)χ̄(kq)χ(k) = 1; intermediate indices stay in the index set; sum = ⟨K^r c,c⟩; ‖|K|‖ = O(1) (Schur test; checked the three cases) bounds the O(1/L) error from G ≈ iφ uniformly in c.
- Verdict: correct. The off-diagonal is super-polynomially small because the Gaussian has width H = T/L while |ω| ≥ L²⁰/2T; length T^{1−o(1)} is genuinely allowed.

## D4. Prop 2 (L404–480)

Squarefree q, q′ force distinct primes coprime to q; μ(nP) = (−1)^rμ(n)² turns −2i into 2i; Mertens gives du/u; uniform first limit of Lemma 2 in v = Σu_i. Limit (2i/π)^r∫∏sin(πφu_i)/u_i∫f(x)f̄(x+Σu)dx = ⟨(A−iφ)^r f,f⟩. ‖(A−iφ)^r‖ ≤ (2φ)^r/r! ⇒ spec A = {iφ}. Correct.

## D5. Lemma 4 (L548–577): the decisive trick

Under the gap hypothesis N ∈ {0,1}. (X+iN)^n: Im E[V^{2j+1}] = Σ_r (−1)^{j−r}C(2j+1,2r)E[X^{2r}1_{N=1}] (top coefficient 2j+1 ≠ 0) and Re E[V^{2j}] = E[X^{2j}] + Σ_{r<j}(−1)^{j−r}C(2j,2r)E[X^{2r}1_{N=1}]: holomorphic moments determine, by triangular induction, every E_T[X^{2j}] and E_T[X^{2r}1_{N=1}], so all absolute moments of V are bounded uniformly in T. Tightness + uniform integrability ⇒ a positive limit measure on the two lines with moments m_j = ⟨A^jf,f⟩/‖f‖². Correct; no absolute-moment input (pair correlation, RH-moment bounds) is needed for A.

## D6. Lemma 3 (L494–546)

ν = push-forward of Ω(dr) = π⁻¹Im⟨(r−C)⁻¹g,g⟩dr, C = e^{πB} (spectrum in the upper half-plane), under r ↦ Log(r+i0)/π. Moments: Riesz–Dunford on the upper half-annulus; (ζ−C*)⁻¹ is analytic there; arcs O((log R)^j/R) and o(1). Uniqueness: the e^{±it/2} combination divided by 2cos(t/2) is analytic in |t| < π ⇒ even moments ≤ C_R(2j)!R^{−2j} ⇒ exponential moments ⇒ ν̂₀(s) + e^{∓s}ν̂₁(s) known ⇒ determinant e^s − e^{−s} ≠ 0. Correct. Numerics (opcheck.py, n = 300 grid): expm(πA) = closed form (L585) to 2e−5; moments of ν (from p₀, p₁) = ⟨A^jg,g⟩, j ≤ 5, φ ∈ {.3,.5,.6,.75}; p₀(0) numeric = sin a sin 2a/(a‖g‖²) (−0.105301 at φ = .6).

## D7. Positivity criterion and sharpness at ½ (ours)

Im(r−C)⁻¹ = (r−C)⁻¹(Im C)(r−C*)⁻¹, Im C = (C−C*)/2i. So the limit measures ν(φ,f) are positive for every f iff Im C ⪰ 0. Kernel: Im C(x,y) = sin a δ(x−y) + a e^{ia(x−y)}e^{ia sgn(x−y)}. With F = e^{−iax}f:

  ⟨Im C f,f⟩ = sin a‖F‖² + a cos a|∫₀¹F|² + a sin a⟨SF,F⟩,  (SF)(x) = i∫₀¹sgn(x−y)F(y)dy,

spec S = {2/((2k+1)π): k ∈ ℤ}, eigenvectors e^{i(2k+1)πx}. Hence φ ≤ ½ ⇒ Im C ⪰ sin(πφ)(1−2φ) ≥ 0: every f gives a positive measure, no contradiction; φ > ½ ⇒ F = e^{−iπx} gives sin a(1−2φ) + 4φ cos a/π < 0. Numerics: λ_min(Im C) = +0.462, 0.000, −0.476, −1.414 at φ = .3, .5, .6, .75. The Möbius-coefficient method proves exactly μ ≤ ½ and cannot prove μ < ½: it needs φ > ½ precisely through sin 2πφ < 0 (A) and sin πφ(1+θ) < 0, θ < 1 (B).

## D8. Theorem 2 (L603–788)

θ ∈ (1/φ−1, 1) makes sin(aθ) > 0 > sin(a(1+θ)). Polynomials are dense in L²(|p₀|dx) (exponential moments), so ∫q²p₀ < 0 for a real q; b₀ with Im b₀(x+i) = Re q(x+i)² exists (triangular in degree), so P = q²+ib₀ has Re P = q² ≥ 0 on the lower line, 0 on the upper line, and Re∫P dν = ∫q²p₀ < 0. Lemma 6 transfers to a fixed rough R_T: Re∫G_T|R_T|²P(V) ≤ −δT. Outside E_T (≥ 2 zeros), Re P(V) ≥ 0; Hölder with ∫|R_T|⁴ ≪ T, ∫|P(V)|⁴ ≪ T (Lemma 5, RH) gives meas E_T ≫ T; each consecutive pair covers ≤ h of centres ⇒ ≫ T log T pairs. Correct. Lemma 5's inputs (Rodgers 2017 tail bound; Re V = Hilbert transform of the zero-count function, L639) are standard under RH. Cosmetic: L665 says "at most ⌊θ/η⌋ prime factors, all exceeding T^η", but supp ρ_η ⊂ (η/2, 2b) gives ⌊2θ/η⌋ and T^{η/2}; harmless.

## D9. Theorem 3 (L790–969)

- Prop 3: numerator move m/n → mp/n and denominator move m/(np) → m/n contribute equally by m ↔ n symmetry; c_{m/n}c_{mp/n}·p^{−½} = (ℓ^{2Ω(mn)}/mn)·ℓ/p·H(m)H(mp)H(n)²; Mertens; log Y/log T → ½ gives sin(πφu/2). Sum/‖c‖² → −4iℓJ/I, so E[N] → φ − (4ℓ/π)J/I. Correct.
- Lemma 7: (log p/log n)_p under ℓ^{2Ω(n)}/n is Poisson–Dirichlet(ℓ²); Ewens moment formula E[Σ_distinct ∏x^{a_U}] = θ^{|P|}Γ(θ)∏Γ(a_U)/Γ(θ+Σa) gives L813–820. Implemented via the exponential formula; sanity e = 1 ⇒ all moments 1 (error 1e−300).
- J_K (L922–931): v = xy, u = x(1−y), Jacobian x; x-integral 1/(ℓ²+r+s+2k+1); y-integral = Beta(ℓ²+r+e·l, 2k+1+e·(k−l)). Re-derived; remainder L934 re-derived.
- Reproduction (largegap.py, prec 1000 bits, I = A_v ⊗ A_t whitened by Cholesky, J̃ = (L_v⁻¹⊗L_tᵀ)𝕎(L_v⁻ᵀ⊗L_t⁻ᵀ), float64 eigvalsh; 225 s, max RSS 267 MB): max J/I and margin 1−φ+(4ℓ/π)J/I at φ = 3.63846723, ℓ = 39327/20000:

| r ≤ | \|j\| ≤ | dim | max J/I | margin |
|---|---|---|---|---|
| 12 | 0 | 13 | 1.027031405460 | −6.716e−2 |
| 12 | 1 | 65 | 1.053623832745 | −5.782e−4 |
| 12 | 2 | 195 | 1.053838841158 | −3.987e−5 |
| 12 | 3 | 455 | 1.053854347236 | −1.045e−6 |
| 12 | 4 | 910 | 1.053854726398 | −9.597e−8 |
| 10 | 5 | 1386 | 1.053854765283 | +1.386e−9 |
| 12 | 5 | 1638 | 1.053854765588851 | **+2.150864e−9** |

  Paper (L953): +2.15088826449e−9. The margin needs the full degree-5 prime-sum space: degree ≤ 4 fails.
- Certificate (cert.py, Arb, prec 1000, 95 s, max RSS 174 MB): b = top eigenvector mapped back to monomials, rounded to exact dyadics; I(b), J₄₀(b) evaluated from the defining formulas in ball arithmetic; Taylor remainder (L934 with K = 40) rem/I = 5.1e−28. Result: J₄₀/I = 1.0538547655888622856 ± 4e−20, rigorous margin **2.15089187038150e−9 ± 6e−25 > 0**. Σ|b|/√I = 7.4617e17 (paper: Σ|b| < 7.46135 at I ≈ 1e−34): same optimizer, extreme cancellation in the monomial basis.

## D10. Model tests

- Alternative Hypothesis (Lagarias–Rodgers): with gaps in ½ℤ and Montgomery's F(α) = |α| on 0 < |α| < 1, F is 2-periodic and the pair density at distance k/2 is w_k = ½ − 2/(π²k²) (k odd), ½ (k even): about 29.7% of gaps equal ½. So AH has μ = ½ and D(½+ε) > 0, agreeing with A and B; by D7 the argument never yields μ < ½. C under AH needs gaps ≥ 4, which AH allows.
- Picket fence Λ(s) = γ(s)+γ(1−s), γ(s) = π^{−s/2}Γ(s/2): functional equation, real on the line, all high zeros on the line at θ(t) ≡ π/2 mod π, so μ = λ = 1. No Euler product ⇒ the prime operator K vanishes ⇒ moments (iφ)^j for every weight ⇒ positive law ν (Stieltjes transform 1/(ζ−e^{iπφ})). models.py: E[V^j] = (iφ)^j (j ≤ 5) and sign Z = (−1)^N exactly. The argument uses the Euler product essentially and claims nothing false here.
- Davenport–Heilbronn (`Riemann/Deep/ExploreHelio/results/dh_zeros.json`: 9 off-line zeros up to 420, e.g. 0.8085+85.699i): off-line pairs share ordinates, so the reduction gives μ = 0, which is true; the analytic core is never reached (no zero-free half-strip, no Λ-supported series for log f).
- Actual zeta zeros (zfinite.py): Prop 1 is a finite-T identity up to e^{−12}-small terms once N is small against H, so it can be tested. Setup: T = 2000, φ = 0.6, h = 0.496, Gaussian centre 3000 with H = 263, and centres t ∈ [1000, 5000]. Z(t) comes from Riemann–Siegel with C₀–C₂ (max error 3.3e−8 against mpmath.siegelz). N(window) is counted by sign changes: 3873 zeros, against 3872.2 from the Riemann–von Mangoldt main term. X = π⁻¹log\|ξ(+)/ξ(−)\|. Prediction: Σ_r C(j,r)γ_{j−r}⟨(K/π)^r c,c⟩/‖c‖², with γ_m the exact Gaussian mean of G^m (γ₁ = −0.1239+0.4866i at this T, not iφ). \|numeric − predicted\| for c = δ₁, Möbius f = 1 and Möbius f = 1−x (N = 20): j=1: 2e−8, 6e−6, 7e−6; j=2: 3e−4, 8e−5, 1e−4; j=3: 7e−4, 8e−4, 8e−4; j=4: 1.5e−3, 5.5e−3, 3.3e−3. Grid steps 0.004, 0.002 and 0.0013 make the j = 4 values scatter on both sides of the prediction (0.118, 0.109, 0.099 against 0.103), so the residual is unbiased quadrature noise from the log singularities. Prop 1 holds numerically at finite T. Neither A nor B has a finite-T statement to test.
