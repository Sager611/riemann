NO IDEA: Girsanov/Esscher change of measure on the Biane–Pitman–Yor law proves Li positivity only as positivity of a probability law (λ₁, λ₂ always, λ₃ whenever λ₁ ≤ 2; an Arb-certified positive symmetric law has λ₃ = −77.30), its Feynman–Kac twin is Polymath15's horizontal-heat mechanism, and the first-passage (Kent/Krein) layer is already Polson's (arXiv 1708.02653, Oct 2026), so nothing moves beyond the node lines (node N3; N2, N13 checked).

# explore56 — stochastic calculus (Girsanov / Esscher tilt, Feynman–Kac, h-transform / first passage) → N3 (N2, N13 checked)
Orca agent explore56, 20 Nov 2026 03:24–06:20 UTC. Scripts here, outputs in logs/. Every process one core, memcap peak ≤ 17 MB.
Prior coverage read first: explore51 (Dyson/zero dynamics, N2), explore39 (h-transform of Weil forms, N13), R10.2 and R342 (BPY positivity allows off-line zeros), R431/R495 (Polson–Sokolov GGC criterion). Synonym scan of NOGO/REPORTs: only those.

## 0. Transfer (one tool: change of measure)
- BPY (math/9912170): Y = √(2/π)·max of a standard Brownian excursion, equivalently Y² = (π/2)∫₀¹|b_s|²ds for a 4-dimensional Brownian bridge (Karhunen–Loève: Y² = (1/π)Σ G_n/n², G_n iid Gamma(2)); then E[Y^s] = 2ξ(s). Checked by hand via E e^{−λY²} = (√(πλ)/sinh√(πλ))² and Mellin (not scripted).
- Tool: the Esscher/Girsanov tilt dQ = Y dP (E Y = 2ξ(1) = 1). FE ξ(s) = ξ(1−s) ⟺ law_Q(1/Y) = law_P(Y) ⟺ under Q_{1/2} ∝ Y^{1/2}P, V = log Y is symmetric (density ∝ Riemann's Φ).
- Statement S1 (exact): log 2ξ(s) = (s−1)·D_s(Q‖P) (Rényi divergence), and λ_n = Σ_{m≤n} C(n,m) κ̃_m/(m−1)!, κ̃_m = (log 2ξ)^{(m)}(1) = cumulants of log Y under Q.
  Hence λ₁ = KL(Q‖P) > 0 and λ₂ = 2KL + Var_Q(log Y) > 0 (BPY 2001; arXiv 2411.08863 for Q(i), Q(√−2)).
- Question transferred: how far does positivity of the tilted law (the whole content of the change of measure) carry Li positivity (N3: RH ⟺ λ_n > 0 ∀n)?

## 1. Statements
- S2 (class, exact): positive symmetric laws ⟷ probability measures ν on [0,∞); the Q-moments are m_k = E_ν[v^k] (k even), E_ν[v^k tanh(v/2)] (k odd). Every such ν gives an entire F(s) = E[Y^s] with F(0) = F(1) = 1, F(s) = F(1−s).
- S3 (λ₃, proved): m₃ ≥ 0 and κ̃₃ = m₃ − 3κ̃₁κ̃₂ − κ̃₁³ give λ₃ = 3κ̃₁ + 3κ̃₂ + κ̃₃/2 ≥ 3κ̃₁ − κ̃₁³/2 + (3/2)κ̃₂(2 − κ̃₁) > 0 whenever λ₁ = κ̃₁ ≤ 2. For ζ, κ̃₁ = 0.0231: λ₃ > 0 by change of measure alone (one coefficient past BPY; λ₃ = 0.2076 is long known).
- S4 (barrier, certified): without the smallness of λ₁ positivity stops at n = 2: law A, M(w) = 0.999 cosh(w/2) + 0.001 cosh(18w), F(s) = M(s−½)/M(½), has λ₁ = 14.375, λ₂ = 80.461, λ₃ = −77.30479 (Arb, 200 bits, li3_certify.py). Its zero s = 0.9232 + 0.1756i (w = 0.42325 + 0.17555i) sits next to s = 1.
- S5 (range given ζ's two cumulants, float evidence only): with κ̃₁, κ̃₂ fixed at ζ's values (0.0230957, 0.0461543) every scanned 3-atom ν (712 feasible) has λ_n > 0 for n ≤ 14, and the minimisers reproduce ζ's own λ₃, λ₄, λ₅ = 0.20764, 0.36880, 0.57560 to 3–4 digits (ζ: 0.207639, 0.368790, 0.575543). Reason: ζ lies 1.06 % inside the class wall m₁ ≤ m₂/2 (tanh(v/2) ≤ v/2), so small λ₁, λ₂ force a concentrated Q, i.e. no low zeros. This is a Li range from 'no zeros below a small height', far below the record n ≤ 106 T².

## 2. Tests (numbers; logs/)
- Smallest case, ζ (li_universal.out): the cumulant formula reproduces λ₁…λ₅ = 0.0230957090, 0.0923457352, 0.2076389206, 0.3687904795, 0.5755427145; Esscher cumulants κ̃₁…κ̃₅ = 0.0230957, 0.0461543, −2.2232·10⁻⁴, −4.4176·10⁻⁴, 1.716·10⁻⁵. So λ₃ = 0.2076 is 99.9 % 'Gaussian' (3κ̃₁ + 3κ̃₂ = 0.20775); κ̃₃ < 0 is the only non-positive input, and S3 absorbs it.
- Class scan (li_universal.out): two-scale laws M = p cosh(aw) + (1−p) cosh(bw), 198 grid points: min λ₁ = 0.1235, min λ₂ = 0.4843 (> 0, as S1 forces), min λ₃ = −77.30, min λ₄ = −888.3, …, min λ₁₀ = −1.0·10⁷: every n ≥ 3 fails somewhere in the class.
- Certificate (li3_certify.out): law A above, Arb balls λ₃ = [−77.3047883197 ± 1.9·10⁻¹¹], λ₁ = [14.3752520159 ± 4.8·10⁻¹¹], λ₂ = [80.4609415086 ± 1.1·10⁻¹¹].
- Bounded scale (li_range2.out): 3-atom laws with ζ's κ̃₁, κ̃₂: 712 feasible, all λ_n > 0, n ≤ 14 (min λ₁₄ = 4.445); two-scale laws matched to ζ's κ̃₁, κ̃₂ at b = 6, 9, 12 give λ₃ = 0.20810, 0.20852, 0.20897 (Arb). A float scan, not a certificate; the ceiling of S5 (largest n forced) is not determined.

## 3. Adversarial (dh_class.out)
- Davenport–Heilbronn (κ = 0.284079043840, conductor 5) is in the class: its theta θ(x) = Σ c_n n e^{−πn²x/5} is self-dual (θ(1/x) = x^{3/2}θ(x), checked to 12 digits) and positive (for x ≥ 1 the n ≥ 2 terms are < 0.33 of the first; x < 1 by self-duality; grid 0.01–100 positive).
- DH: λ₁ = 0.0976361, λ₂ = 0.3888255, λ₃ = 0.8684689, and S3's bound gives λ₃ ≥ 0.84476 from positivity alone, although DH has off-line zeros (0.8085 + 85.70i first). So S1–S3 do not tell ζ from DH; with law A and Epstein-type kernels (R10.2) this is the known 'positive laws allow off-line zeros'. The method sees only that λ₁, λ₂ are small, i.e. that there are no zeros near s = 0, 1.

## 4. Side transfers checked (no new statement)
- Feynman–Kac → N2 (fk_horizontal.out): H_t(z) = (1/8)E[ξ(½ + iz/2 + √(t/2)·G)], G ~ N(0,1), t ≥ 0 (Hubbard–Stratonovich of e^{tu²}); the two sides agree to < 10⁻³⁰ at t = 0, 0.05, 0.1, 0.2 and z = 0, 10, 28.3, 40, 60. So the de Bruijn–Newman flow is the heat semigroup in σ = Re s. Since log|ξ(σ+iT)| ≈ (L/2)|σ − ½| (L = log(T/2π)), the Gaussian weight has two mirror peaks σ = μ ∓ tL/4 (at T = 10⁶, t = 0.1, μ = 0.4 the left peak is at 0.14, predicted 0.10), and a zero needs the two to cancel: this is Polymath15's two-term (σ_eff) approximation, whose ceiling with present inputs is explore51's 4(2θ − 1) / 0.0957. Variation diminishing (the Gaussian kernel is totally positive) fails here: Im ξ(σ + iT) changes sign infinitely often in σ (the Γ phase turns by ≈ π/4 per unit σ), so sign-change counts are infinite.
- Kent's first-passage theorem / Krein strings → N13: pre-empted. Polson, arXiv 1708.02653 (TeX of the Oct 2026 version read): the reciprocal ξ(½)/ξ(½ + √s) is unconditionally the Laplace transform of a self-decomposable law with canonical function Θ(x) = Σ_{γ>0} e^{(ρ−½)²x}; Θ completely monotone ⟺ RH (GGC = first passage of a string); Θ monotone to order H²/80; DH sanity tests. We had re-derived Θ > 0 for all x from 'zeros below height 100 on the line, |β − ½| < ½, N(T) ≤ (T/2π)log(T/2π) for T ≥ 100 (HSW)' before finding it there; dropped.
- Itô–Tanaka / local time: counts real zeros only through Kac's formula, i.e. the first Laguerre inequality (θ′ ≤ 0 for arg(H + iH′)), N6 and Fleischer (ExploreQuant): no transfer.

## 5. Novelty
- arXiv API (logs not kept; titles): 'Riemann'+'Girsanov', 'zeta'+'Cameron–Martin', 'de Bruijn–Newman'+'Brownian', 'Li coefficients'+'probabilistic' (→ 2411.08863, λ₁, λ₂ for Q(i), Q(√−2)), 'Biane'+'excursion' (→ math/9912170), 'Riemann hypothesis'+'Brownian' (→ 1708.02653, 2611.04359 Langevin model: GHS-type susceptibility bound, no Li claim; 0811.1490 Biane on Pólya's ξ*; 2312.12760). Monthly listings Oct–Nov 2026, math.NT/PR/CA/CV (lst.py): no Li/Esscher/Girsanov item. S3 (λ₃ from positivity when λ₁ ≤ 2) not found; it is routine (one moment sign) and below every Li record.

## Suggested NOGO line
- explore56 (Orca, exploration: stochastic calculus — Girsanov/Esscher change of measure on the Biane–Pitman–Yor law → N3; Feynman–Kac → N2, Kent first passage → N13 checked) NO IDEA: change of measure gives Li positivity only as positivity of a law: λ₁ = KL, λ₂ (BPY), λ₃ ≥ 3λ₁ − λ₁³/2 + 1.5κ̃₂(2 − λ₁) > 0 when λ₁ ≤ 2 (ζ, and Davenport–Heilbronn, whose theta is positive: λ₃ = 0.868); an Arb-certified positive symmetric law has λ₃ = −77.30; Feynman–Kac makes H_t the horizontal heat semigroup of ξ (Polymath15's two-term mechanism, explore51's ceiling); the first-passage/Thorin layer is Polson 1708.02653 [ExploreGirsanov](Riemann/Deep/ExploreGirsanov/REPORT.md)

RH relevance: N3 (why probabilistic Li positivity stops after an initial segment, and is Davenport–Heilbronn-blind); N2 and N13 side checks.

## Files
li_universal.py, li3_certify.py (Arb), li_range2.py, dh_class.py, fk_horizontal.py; outputs logs/*.out (memcap peaks in logs/*.err). The brief's NOGO template began 'explore56 (Orca, exploration: anti-concentration …', explore55's field; the line above names this field.
