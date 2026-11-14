# FULLPROOF (robinchaudit): RobinCH's chain, step by step, with source lines

Agent robinchaudit, 14 Nov 2026. RCH = Riemann/Deep/RobinCH/PROOF.md, RHP = Riemann/Deep/RobinHeight/PROOF.md,
CH = Chirre–Helfgott arXiv 2512.15709v1 e-print VonMangoldt.tex (fetched 14 Nov 2026, the only version: tgz dated 18 Dec 2025).
Code: c1 = certify.py, c2 = certify2.py, c3 = certify3.py (RobinCH copies). ✓ = re-derived here; (inh.) = inherited, statement level.

## A. Inputs

| # | used for | check |
|---|---|---|
| H1 RH to 3·10^12 (Platt–Trudgian 3,000,175,332,800) | β = 1/2 for abs(γ) ≤ T0; CH's hypothesis | (inh.) statement |
| H2 half-plane b0 (0.87338… / 7/8) | only zeros abs(γ) > T0, only rows b0 < 1 | (inh.) not audited here |
| H3 HSW N(T) bound, R(t) = 0.1038 log t + 0.2573 log log t + 9.3675 | S_m(T0), S1(T0) | (inh.); integrals ✓ (§E) |
| H4 Σ 1/(ρ(1−ρ)) = 2 + γ − log 4π | low-zero sums | ✓ (classical: 2 Σ 1/ρ) |
| H5 Ingham ψ_1 explicit formula | Lemma 5 / 11_k decomposition | (inh.) |
| H6 Büthe: θ(x) < x ≤ 10^19; x − θ ≤ 1.95√x (1423 ≤ x ≤ 10^19) | η below 10^19, M_low, Lemma 9, ψ−θ ≤ P | (inh.) |
| H7 θ(y) < y log 4 | M_up beyond 10^76 | (inh.) classical |
| H9 Morrill–Platt: Robin for 5041 ≤ n ≤ 10^(10^13.11485) | CA numbers with x_1 < 10^10 (log N < 1.01·10^10) | (inh.) |
| H10 CH Cor. cor:psiexpl | η above 10^19 | ✓ statement vs TeX (§B); proof chain (inh.) CHAudit §2.2 |

## B. H10 against the TeX

- CH l.387–388: "RH holds up to height T" := all zeros with 0 < Im s ≤ T have Re s = 1/2. = RCH l.13–14 ✓.
- CH l.389–397 Cor. cor:psi: RH to T ≥ 10^7, x > max(T, 10^9): abs(ψ(x) − x(π/T)coth(π/T)) ≤ πx/(T−1) + C_T √x. ✓
- CH l.402–406 Cor. cor:psiexpl: "For any x ≥ 1, abs(ψ(x) − x) ≤ π/(3·10^12) x + 113.67 √x". Verbatim = RCH l.10 ✓.
- CH l.4981–4991 (proof): x ≥ 3·10^12 + 3 via Cor. psi with T = 3·10^12 + 1 + π/3 (RH to that height, Platt–Trudgian);
  1 ≤ (π/T)coth(π/T) ≤ 1 + π²/(3T²); π²/(3T²) + π/(T−1) ≤ π/(T−1−π/3) = π/(3·10^12). ✓ (rerun (0): True;
  by hand: π/(T−1−π/3) − π/(T−1) = (π²/3)/((T−1)(T−1−π/3)) > π²/(3T²)).
- CH l.4868–4880 Prop. prop:sagaro (σ = 0): hypothesis RH to T ≥ 10^7 only; x > max(T, 10^9). Proof uses Lemma hardin
  (T′ ∈ [T − 1/2, T]); CHAudit §2.2 rechecked its assembly; nothing above height T is assumed. ✓ (inh.)
- C_T at T = 3·10^12 + 1 + π/3: 113.6688806 ± 5·10^−8 < 113.67 (rerun (0)) ✓.
- RobinCH uses H10 only for x ≥ 10^19 e^(−10^−3) > 3·10^12 + 3 (c3 best_eta3: Büthe only if L.upper ≤ log 10^19 − 10^−3). ✓
- No step uses RH above T0: the only zero inputs are H1 (abs γ ≤ T0), H2 (abs γ > T0, b0 < 1 rows), H3, H4, H5; H10's proof
  assumes RH to 3·10^12 + 2.05 < 3,000,175,332,800. ✓

## C. Lemma 12 (RCH l.22–41)

- ψ − θ = Σ_{k≥2} θ(x^{1/k}); for x ≤ 10^38 each x^{1/k} ≤ 10^19, θ(y) < y (H6) ⇒ 0 ≤ ψ − θ ≤ √x + x^{1/3} + x^{1/4} + (log x/log 2) x^{1/5} =: P_0(x)
  (k ≥ 5 terms: at most log x/log 2 of them, each < x^{1/5}). ✓
- θ − x ≤ ψ − x ≤ πx/3e12 + 113.67√x; x − θ = (x − ψ) + (ψ − θ) ≤ πx/3e12 + 113.67√x + P_0. ⇒ abs(θ − x) ≤ η_CH x. ✓
- η_1 := η_CH(10^19) (factor 1) bounds abs(θ(y)/y − 1) on [10^19, 10^38] (η_CH decreasing for fixed factor: each term decreasing,
  (log x) x^{−4/5} decreasing for x > e^{5/4}). Rerun: η_1 = 3.62631·10^−8. ✓
- 10^38 < x ≤ 10^76: x^{1/2}, x^{1/3} may exceed 10^19 ⇒ θ(x^{1/k}) ≤ (1 + η_1) x^{1/k}; P = (1 + η_1) P_0 (generous for k ≥ 4). ✓
- Jump at 10^38 (factor switch): c3 uses the (1 + η_1) formula whenever L.upper + 10^−3 > log 10^38 (margins4 l. "e1 = ..."). ✓
- Window uniformity: margins4 evaluates η at the ball L − c (c = λ/2) — covers [L_a − c, L_b − c]; η_CH decreasing ⇒ its value
  at L_a − c bounds abs(θ/t − 1) on the whole window [x e^−c, x e^c] of every x in the box. ✓ Büthe branch: window top
  ≤ (L − c).upper + λ < log 10^19 (λ ≤ 9·10^−4). ✓

## D. Reductions (RHP §1–2, used unchanged)

- Lemma 1 (RHP l.27–32): log P(x) = γ + log log x + ε(x)/(x log x) − I(x) − S_2(x), I = ∫_x^∞ ε w. Partial summation + Mertens. ✓
- Lemma 2 (RHP l.36–42): Nicolas at p_k = x ⇐ −I(x) > S_2(x) (D ≥ 0 since u − log(1+u) ≥ 0). ✓ S_2 ≤ 1/(2(x−1)):
  Σ_{m≥2} 1/(m p^m) ≤ 1/(2p(p−1)), Σ_{n>x} 1/(2n(n−1)) = 1/(2⌊x⌋). ✓
- CA numbers (RHP l.47–50): σ(p^a)/σ(p^{a−1}) ≥ p^{1+ε} ⟺ (σ(p^a)/p^a)/(σ(p^{a−1})/p^{a−1}) = 1 + 1/(p + … + p^a) ≥ p^ε ⟺ F(p,a) ≥ ε. ✓
- Lemma 3 (RHP l.52–58): N± maximise σ(n)/n^{1+ε_c} ⇒ log(σ(n)/n) ≤ c + ε_c log n; G(u) = c + ε_c u − log log u has
  G'' = (log u + 1)/(u log u)² > 0 ⇒ max at the ends. ✓
- Lemma 4 (RHP l.59–68): log log log N ≥ log log θ(x) + L_2 from log N ≥ θ(x) + θ(x_2); primes in (x_2, x] have a_p = 1,
  log(1 − p^−2) ≤ −p^−2. ✓ Robin at N_ε ⇐ −I + D < Σ_{x_2<p≤x} p^−2 + L_2 + S_2 (S_2 dropped in the code: safe). ✓
- Lemma 9 (RHP l.118–124): g = 1/(t² log t), ∫ g = E_1(log x_2) − E_1(log x); boundary ε(x)g(x) ≥ −1.95 x^{−3/2}/log x,
  −ε(x_2)g(x_2) ≥ 0; ∫ ε abs(g') ≥ −1.95·3∫ t^{−5/2}/log x_2 = −3.9 x_2^{−3/2}/log x_2. ✓ η_1 piece (RHP l.197–200):
  −η_1(1/(b log b) + 1/(a log a) + 3/(a log a)). ✓
- Theorem B (RHP l.149–154) + H9: CA numbers with x_1 < 10^10 have log N < 1.01·10^10 < 10^13.11485 ln 10 = 2.9996·10^13. ✓

## E. Explicit formula, zero sums (RHP §3, RCH §2)

- Lemma 5 (RHP l.70–88): t − ψ = Σ t^ρ/ρ + log 2π + ½log(1 − t^−2) via Ingham's ψ_1; sign of each piece ✓; b0 = 1 (RHP l.87–88):
  Tonelli with Σ 1/((1−β)|ρ|²) < ∞ (classical zero-free region gives 1 − β ≫ 1/log γ). ✓
- Lemma 6 (RHP l.90–96): two integrations by parts; ∫_x^∞ t h'' = h(x) + x abs(h'(x)); κ = 1 + 2(L+2)/(L(L+1)) from
  x abs(h')/h = (L+2)/(L(L+1)). ✓
- Lemma 7 (RHP l.97–109): Re ρ(1−ρ) = β(1−β) + γ² > 0 ⇒ Σ_{abs γ≤T0} abs(ρ)^−2 ≤ β_0 = 0.0461914; pairing ρ ↔ 1 − ρ̄ fixes γ,
  β ↦ x^{β−1} + x^{−β} convex. ✓
- S_m(T) (RHP l.191–193; c2.Sm): 2[∫_T^∞ t^−m log(t/2π)/2π + R(T)/T^m + m∫_T^∞ R t^{−m−1}], R ≤ (c1+c2+c3) log t (t ≥ e),
  m∫_T^∞ log t·t^{−m−1} = (m log T + 1)/(m T^m). ✓ (rederived). S1(T0) (c1.S1U): 2[N(U)/U + ∫_14^U N/t²], ∫ M/t² =
  (log²(U/2πe) − log²(14/2πe))/4π, ∫ R/t² ≤ (c1+c2+c3)(log 14 + 1)/14. ✓ Value 119.99 vs main term 2·log²(T0/2π)/4π = 115.1.
- Lemma 11_k (RCH l.43–77): every step re-derived:
  (i) 1_{t≥x} = φ + (1_{t≥x} − φ); ω = wφ continuous, piecewise C¹, ω(x') = 0, abs(ω') ≤ abs(w') + w·(k/λ)/t = O(t^−3/log t)
      (density of V ≤ k/λ) ⇒ Lemma 5's termwise argument applies. ✓
  (ii) φ ∈ [1_{t ≥ x'e^λ}, 1_{t ≥ x'}] ⇒ M(x'e^λ) ≤ M_φ ≤ M(x'); ∫ω ≤ 1/(x'L'); abs(V_φ) ≤ x'^−3/L'. ✓
  (iii) t w dt = H(u) du; abs(W) ≤ η_win H(L') m_c; m_c = E abs(V − c) (layer-cake), m_0 = λ/2, m_{λ/2} = (λ/k) MAD_k. ✓
  (iv) (B1'): integrand supported on u ≥ L', abs(e^{su}) H(u) ≤ x'^{β−1} H(L'). ✓
  (v) (B2_k): F_φ = E F(ρ; L'+V) (Fubini, β < 1); F(ρ;y) = −e^{sy}H/s + e^{sy}H'/s² + s^−2 ∫ e^{su}H''; E e^{sV} =
      ((e^{sa}−1)/(sa))^k, abs(e^{sa} − 1) ≤ 2; abs(H(L'+v) − H(L')) ≤ v abs(H'(L')) (H'' > 0); E V = λ/2. ✓
- Zero sum (RCH l.79–85): low zeros abs(F/ρ) ≤ x^{−1/2}h(1+κ/abs ρ)/abs(ρ)² + x'^{−1/2}H(L')m_c/abs(ρ), summed with β_0, β_0/γ_1, S1(T0) ✓;
  high zeros abs(ρ), abs(s) ≥ abs(γ), pairing in x' ✓; code c3.nZlow4/nZhigh4: normalisation √x L·x'^{−1/2} = L e^{c/2} ✓,
  P = √x L (x'^{b0−1} + x'^{−b0})/2 ✓.

## F. Lemma 13 (RCH l.87–95)

- D = u/L − log(1+v) = (u − log(1+u))/L + (v − log(1+v)); z − log(1+z) ≤ z²/2 + abs(z)³/(3(1−abs z)); abs(log(1+u)) ≤ abs(u) + u²
  (abs u ≤ ½) ⇒ abs(v) ≤ η(1+η)/L. Inserting gives exactly RCH l.90 (each bracket term checked by hand). ✓


## G. Code against formulas (c3.margins4, c2.robin_pos, c2.Mlow3/Mup3, c1.box/run)

| term (normalised by n(x) = √x log x) | formula | code | check |
|---|---|---|---|
| M, Nicolas | M(x'e^λ) = M(x e^c) | (L/2).exp()*L*K.Mlow3(L + c) | ✓ |
| M, Robin | M(x') = M(x e^−c) | (L/2).exp()*L*K.Mup3(L − c) | ✓ |
| Mlow3 ≤ 10^38 | [I½(x) − I½(Y1)] − 1.95[I¼(x) − I¼(Y1)] + (1−η_1)[I½(Y1) − I½(Y2)] | same | ✓ (θ(√t) ≥ √t − 1.95 t^¼; ≥ (1−η_1)√t on (Y1, Y2]) |
| Mup3 ≤ 10^38 | base(x) + η_1 base(Y1) + (log 4 − 1)(I½ + I⅓)(Y2) + log 4·J¼(Y2)/log 2 | same | ✓ (H7 beyond 10^76) |
| low zeros | (1 + 1/L)β_0(1 + κ/γ_1) + L e^{c/2} H(L') m_c S1(T0) | nZlow4 | ✓ |
| high zeros | √x L (x'^{b0−1} + x'^{−b0})/2 [H(L')(2k/λ)^k S_{k+2} + abs(H'(L'))(λ/2 S_2 + 2 S_3)] | nZhigh4 | ✓ |
| window | √x L η H(L') (λ/k) MAD_k | nWin4 | ✓ |
| log 2π, V_φ | +log 2π e^{−L/2} e^c L/L' (Robin); −e^{−5L/2} e^{3c} L/L' (Nicolas) | same | ✓ signs |
| S_2 | −√x L/(2(x−1)) (Nicolas); dropped (Robin, ≥ 0) | same | ✓ |
| D (Lemma 13) | √x η²/2 [1 + 2η/(3(1−η)) + (1+η)²/L + 2η(1+η)³/(3L²(1 − η(1+η)/L))] | nDup2 | ✓ |
| Robin positive | n(x_a)[Σ_{x2hi<p≤x_a} p^−2 (Lemma 9 + η_1 piece) + L_2(θ2 ≥ θ(x2lo), θ1 ≥ x_b(1+η))] | robin_pos | ✓ (x_2 increasing in x, enclosed by arb-verified bisection) |
| box | ball ⊇ [La.lower, Lb.upper] | arb(mid, rad) adds rad to the error | ✓ |
| k, λ | k ≥ 1, λ ≤ 9·10^−4 < 10^−3, exact binary floats | choose4, min(·, 9e−4) | ✓ |

No discrepancy found. Ball evaluation encloses each pointwise bound over the box, so m.lower() > 0 proves the margin for every x in the box.

## H. Reruns and independent recomputation (logs in this folder)

- rerun.log: check_all2.py on an unmodified copy (.scratch/robinchaudit/run), memcap 290 MB: ALL PASS, every line identical to
  RobinCH/REPORT.md §4 (box counts, minima, k bounds, exponents 40.5656 / 40.4896 / 35.3789 / 46.7931 / 46.7900 / 46.4480);
  370 s wall, max RSS 87.6 MB.
- dump.log / cover.log: RobinCH's own Level-3 search with every box recorded (3,513 boxes): frontier never skips (0 gaps),
  first box starts at 10^7 / 10^10, last box ends ≥ the claimed endpoint (log10 48.250420, 47.155336, 47.799341, 46.810233).
  Endpoint probes (search continued past the claim): failure at 10^48.252447 (N, 0.87338), 10^47.156207 (R, 0.87338),
  10^47.799393 (N, b0 = 1), 10^46.810759 (R, b0 = 1) = RobinCH's search endpoints (last digit ±1: other box sequence).
- indep.log: indep.py, own code from the formulas (256 bits, monotone terms at box ends), 30 boxes with RobinCH's (k, λ):
  every RobinCH PASS box also passes here, every probe FAIL box also fails; our m ≥ RobinCH's m (ball dependency) by ≤ 2.4%
  on large margins, e.g. min Robin box (b0 = 0.87338, 10^47.1483–10^47.1491): 2.04·10^−6 vs 1.63·10^−6; first Robin box
  (10^10): 0.2152 vs 0.2147; Nicolas at 5·10^8: 1.6024 vs 1.5727. Own (k, λ) on log x ∈ [20.03, 20.04] (OddRobin's Robin box):
  0.468163 = oddraudit F4's 0.468.
- Constants recomputed: S1(T0) = 119.98765, S_2(T0) = 2.95941·10^−12, S_3(T0) = 4.84393·10^−25, β_0 = 0.0461914179322,
  η_1 = 3.6263146·10^−8, η_CH(10^40) = 1.04720·10^−12, C_T = 113.66888; E_1 (flint) = mpmath at 15 digits (two points);
  MAD_1 = 1/4, MAD_2 = 1/3 (by hand), MAD_40 = 1.458561 (normal approximation 1.456731).


## I. Falsification (logs in this folder)

- Lemma 11_k at corners (fals.py → fals.log; method independent of RobinCH/test_b2k.py: F_φ = E F(ρ; L'+V) as a Taylor
  series of G(y) = e^{−sy}F(ρ;y) with exact E_n closed forms against the exact moments E[V^n e^{sV}] of the MGF
  ((e^{sa}−1)/(sa))^k; truncation ≤ 1.7·10^−22): 612 valid cases, k ∈ {1, 2, 40, 57, 120} (certificate uses k ≤ 57),
  γ ∈ {14.134725, 3·10^12 + 1, 10^16}, L ∈ {16.118 (10^7), 111.1 (1.78·10^48)}, β ∈ {0.5, 0.87338…, 0.999}, λγ/k ∈ {π, 2, 4}
  or λ = 10^−3, centred and one-sided, λ ≤ 10^−3 deliberately violated at γ = 14.13 (the proof does not use it).
  (B2_k): never exceeded; tightest 1 − ratio = 6.7·10^−18 (k = 1, γ = 10^16, β = 0.999, λγ/k = π: e^{sa} ≈ −1, where
  abs(e^{sa} − 1) ≤ 2 is sharp — the lemma is optimal there); 60 cases within 10^−6 of the bound, all at λγ/k = π.
  (B1'): worst ratio 0.999995 (k = 120, γ = 14.13, λ = 10^−3; sharp as λγ → 0).
  Method cross-checked against direct quadrature (qcheck.log (a)): agreement 7·10^−29, 4·10^−25, 8·10^−35, and 7·10^−10
  (k = 1, γ = 14.13: quadrature accuracy).
- Lemma 13 (lemma13.log): 144 cases, η ∈ [10^−12, 1/2], u/η ∈ ±{1, 0.99, 0.5}, L ∈ [20, 111.1]: min(1 − D/bound) =
  8.9·10^−15 > 0 (sharp as η → 0, as it must be: both sides ~ η²/(2L)).
- MAD_k (qcheck.log (c)): closed form = quadrature of abs(s − k/2) f_S(s) to 20 digits, k = 3, 12, 40.
- One-sided bounds: every sign in table G checked (M at x'e^λ for Nicolas, at x' for Robin; V_φ ≤ 0 and log 2π∫ω ≥ 0 dropped
  on the safe side; S_2 ≥ 0 dropped for Robin; θ(x_2) from below, θ(x) from above in L_2). No error.
- Box coverage (cover.log): see H. Endpoints (endfine.log, endreach.log): at RobinCH's failure boxes our finer search
  (4 best k, λ × 0.5…2 in 31 steps) also fails for b0 = 0.87338… (N −1.34·10^−6, R −6.55·10^−7): these endpoints are the
  method's limit. For b0 = 1, k = 55 instead of RobinCH's k = 56 gives +0.0027 (N) and +0.0009 (R) at RobinCH's failure
  points, but every point box 0.002 dex further fails: the b0 = 1 claims are conservative by < 0.002 dex.
- RH above T0: none used (B). Robin below 10^10: H9 (D).

## J. Inputs at statement level (TeX e-prints, fetched 14 Nov 2026)

- Platt–Trudgian 2004.09765 main.tex l.62 (theorem): "RH is true up to height 3 000 175 332 800" (title and abstract l.46:
  3·10^12, which alone would not cover CH's T = 3·10^12 + 2.05; the theorem does). = H1 and CH l.402. ✓
- Büthe 1511.02032 PRH2F.tex l.37–39: abs(x − ψ) ≤ 0.94√x (11 < x ≤ 10^19); x − ϑ ≤ 1.95√x (1423 ≤ x ≤ 10^19);
  x − ϑ > 0.05√x (1 ≤ x ≤ 10^19). = H6. ✓
- Hasanalizade–Shen–Wong 2107.06506 main.tex abstract l.165–170, table l.221: abs(N(T) − (T/2π) log(T/2πe)) ≤
  0.1038 log T + 0.2573 log log T + 9.3675, T ≥ e. = H3. ✓
- Not re-checked: H2 (half-plane, Riemann/Deep/QRH; only the b0 < 1 rows), H5 (Ingham), H7, H9 (Morrill–Platt), Robin 1984 /
  Nicolas 1983 (RHP quotes), python-flint/arb (E_1 = −Ei(−z): matches mpmath to 15 digits at z = 25, 55.5).
