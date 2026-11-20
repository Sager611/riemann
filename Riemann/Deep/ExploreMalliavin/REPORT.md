BARRIER: ζ's shift-dependent fourth-order means cannot lift the 1/√3 ceiling: they are carried by a vanishing length (sparse CUE windows plus a background match 𝒮, every normalized microscopic fourth-order mean and the kurtosis with zero density 1/√3, and → 0 with a low-frequency bulk that also copies ζ's one-point law, Selberg's CLT and L¹ ratio), the amplitude-modulation escape X = A·Y is impossible under an exact sinc spectrum, and no finite-chaos model even matches ζ's fourth order (node N8)

# explore62 — Malliavin calculus (fourth-moment / squared-chaos inequalities, Kac–Rice in finite chaos) → N8 (Orca exploration, 20 Nov 2026)

## 0. Repo check
- `rg -n -i 'Malliavin|Nourdin|Peccati|Wiener chaos|fourth.moment theorem|Kac.Rice|Stein.s method' NOGO.md TREE.md Riemann/Deep/*/REPORT.md`: R343 (Kac–Rice: "random-model counts do not transfer"),
  R347 (Stein's method: Berry–Esseen only), R476/explore14 (Beurling–Malliavin, Malliavin's synthesis failure: other Malliavin). No fourth-moment theorem, no squared-chaos inequality.
- Read first: explore60 (ExploreOsc §1 C3: Gaussian witness, 1/√3), ceilaudit7 (CeilAudit7 §B: block mixtures 10.1%; open: G(x); heuristic X = A·Y).

## 1. Transfer (tool, node, statements)
Node N8, explore60's class 𝒮 (all untwisted second-order unfolded means → c(−1)^i sinc^{(i+j)}(u)), unfolded units (maximal zero density 1).
Tool (Malliavin–Peccati school): Malicet–Nourdin–Peccati–Poly, arXiv 1503.02154, Thm 1.1 (TeX checked): F_i in the k_i-th Wiener chaos (F_i ∈ Ker(L + k_i I))
⇒ E[∏F_i²] ≥ ∏E[F_i²], equality iff jointly independent; with hypercontractivity E F⁴ ≤ 9^q (E F²)² in the q-th chaos.
Notation: K_X(x) = E[X(0)²X(x)²]/(E X²)², M4 = K_X(0), G_X = K_X/M4 (ceilaudit7's G).

### 1a. Lemma 1 (independent positive modulation is rigid; refutes ceilaudit7's heuristic)
Statement. Let A ≥ 0 be stationary, not a.s. 0, Y stationary (any law), A and Y independent, X = A·Y (E A² < ∞ follows from E X² < ∞). If the covariance
of X is exactly sinc (spectral measure uniform on [−π, π]), then A is a.s. constant in time. In particular sinc/R̃_A (R̃_A = R_A/R_A(0),
R_A(u) = E A(0)A(u)) is positive definite only for constant A, so ceilaudit7's route "X = A·Y keeps Y's zeros, λ₂(Y) = π²/3 − |R_A″(0)| < π²/3"
never applies, and its "true to first order" is false: for R̃ = (1 + vρ)/(1 + v), FT(sinc/R̃) = (1+v)𝟙 − v(𝟙 * ρ̂) + O(v²) is < 0 just outside [−π, π].
Explicit: A = 1 + εcos(at + U), d = ε²/2: 1/(1 + d cos x) ∝ 1 + 2Σ_{k≥1}(−r)^k cos kx, r = (1 − √(1−d²))/d, so on (π, π + a) FT(sinc/R̃) ∝ Σ_{k≥1, |ξ−ka|≤π}(−r)^k < 0.
Proof. ν = spectral measure of A, μ of Y (μ(ℝ) = 1 after scaling): ν * μ = U := 𝟙_{[−π,π]}dξ/2π. Supports add, so supp ν ⊂ [−a, a], supp μ ⊂ [−b, b]
with a + b = π (a = max supp ν ≥ 0). Atoms of ν: ν({ω}) = E|c_ω|², c_ω = lim (1/2T)∫A(t)e^{−iωt}dt; A ≥ 0 gives |c_ω| ≤ c_0 a.s., so
ν({ω}) ≤ ν({0}) =: m₀ > 0. From ν ≥ m₀δ₀: μ ≤ U/m₀ (density ≤ 1/(2πm₀)). Edge: δ/2π = U([π−δ, π]) ≤ ν([a−δ, a])·μ([b−δ, b]) ≤ ν([a−δ, a])·δ/(2πm₀),
so ν([a−δ, a]) ≥ m₀ for all δ > 0. If a > 0 this gives ν({a}) ≥ m₀, hence |c_a| = c_0 a.s.; since c_a(A(· + s)) = e^{ias}c_a(A), B(t) = A(t)(1 − cos(at + arg c_a))
is a stationary process ≥ 0 with ergodic mean c_0 − |c_a| = 0, so E B(0) = 0 and A(t) = 0 off the lattice {at ≡ −arg c_a mod 2π}: A = 0 a.e., contradiction. So a = 0 and ν = m₀δ₀: A is constant. ∎ (E|A(t) − A(0)|² = 2(R_A(0) − R_A(t)) = 0.)

### 1b. Proposition 2 (transferred tool: squared-chaos inequality; what finite-chaos models can and cannot do)
(i) MNPP Thm 1.1 with d = 2: if every X(t) lies in the q-th Wiener chaos of one Gaussian field, then K_X(x) ≥ 1 for all x (equality iff X(0), X(x)
independent); with hypercontractivity ‖F‖₄ ≤ 3^{q/2}‖F‖₂: G_X(x) = K_X(x)/M4 ≥ 1/M4 ≥ 9^{−q} (exact constants 1/3 for q = 1, 1/15 for q = 2:
second-chaos kurtosis 3 + 12Σλ⁴/(Σλ²)² ≤ 15). Long-block Gaussian mixtures (ceilaudit7): K ≥ Σp_kA_k² ≥ (Σp_kA_k)² = 1 too.
(ii) ζ/CUE violate both, so no finite-chaos (or Gaussian-mixture) model reproduces ζ's fourth order: G_ζ(x) ≈ 3/(π²x²) (x²G = 0.3010, 0.3013, 0.3025
at x = 4, 8, 16, N = 400; logs/cue4.log) is < 1/3 for x > 0.95, < 1/15 for x > 2.14, < 9^{−3} for x > 14.9; and K < 1 at macroscopic shifts:
**K_N(π) = (⌈N/2⌉+1)(⌊N/2⌋+1)/(N+1)² → 1/4** (exact: the symbol 2 − 2cos 2θ splits the Toeplitz matrix into two (2, −1) tridiagonal blocks;
integer check N ≤ 100, e.g. 2601/10201 at N = 100), K_N(φ) < 1 for φ ≳ 1.08 (Fisher–Hartwig 1/(4sin²(φ/2)): φ > π/3). ζ itself at the heights of §2c:
K(4) = 0.354, 0.342, 0.341, 0.388, 0.455 at t = 10⁴…10⁸ (x = 4 is already φ ≈ 2π·4/L there).
(iii) Inside the family the floor does rise: X = Y₁Y₂ (second chaos, Y_i independent Gaussian, λ₂(Y₁) + λ₂(Y₂) = π²/3 forced by U = μ₁ * μ₂) has expected zero
density (√λ₂(Y₁) + √λ₂(Y₂))/π ≥ 1/√3 (Rice; a.s. when both spectra are continuous), equality iff one factor is constant (zeros of a product are the
union). Reviewer's caveat: a line-spectrum factor, μ₁ = (δ_{−2π/3} + δ₀ + δ_{2π/3})/3 with μ₂ = U[−π/3, π/3], has zero-free Y₁ with positive probability
(density 1/(3√3)), but those realizations' own means are not ∝ sinc, so they are not in 𝒮. So the brief's suggested
"H₂-type or products of Gaussians" witnesses cannot test the question: they cannot match ζ's G, and where they live the floor is ≥ 1/√3.


### 1c. Theorem 3 (dilution: microscopic fourth-order means cannot force zeros) — the answer to the brief's question
Unfolding as explore60: F_f(v) = f(t)/√(2θ′(t)), v = θ(t)/π, ⟨·⟩_V = V^{−1}∫_V^{2V}dv. Inputs on f (ζ's values on the right):
(I2) 𝒮: ⟨F^{(i)}(v)F^{(j)}(v+u)⟩_V → c(−1)^i sinc^{(i+j)}(u), i, j ≤ J;
(I4) every normalized microscopic fourth-order mean: ⟨F^{(i₀)}(v)F^{(i₁)}(v+x₁)F^{(i₂)}(v+x₂)F^{(i₃)}(v+x₃)⟩_V/⟨F⁴⟩_V → ζ's limit (fixed x's, i's ≤ J;
contains G(x), R = π²/15 and the full four-point function);  (Ik) ⟨F⁴⟩_V/⟨F²⟩_V² ~ M4_ζ(V) ~ L²/(2π²) (Ingham);
(I1, optional) ⟨|F|⟩_V/⟨F²⟩_V^{1/2} ~ M1_ζ(V) ≍ L^{−1/4} (Heath-Brown, k = 1/2 both bounds).
Statement. There are f satisfying (I2) + (I4) + (Ik) with (1/√3 + o(1))·(θ(2T) − θ(T))/π zeros in [T, 2T], and for each q ∈ (0, 1] with (q/√3 + o(1))·(…);
adding (I1), with o(1)·(…). So the floor over these inputs is 0 (≤ 1/√3 already with a Gaussian background), not above 1/√3: the barrier extends.
Construction. f = background g outside windows, f = λW on windows, W = independent CUE characteristic polynomials of size N = N(V) → ∞ in unfolded units
(one full circle = N units per window; real version W(x) = ∏_j 2sin((2πx/N − θ_j)/2)/√(N+1), ⟨W²⟩ = 1), length fraction p = η², λ² = 1/η, energy fraction η → 0, smooth
C^J junctions of unfolded length 1; background scaled so (1 − p)c_g + η = 1; N chosen with M4_CUE(N) = (N+2)²(N+3)/(12(N+1)) ~ M4_ζ(V) (N ≈ √6 L/π).
- (I2): CUE covariance = normalized Dirichlet kernel → sinc with all derivatives (exact); background ∈ 𝒮; constant 1. Junctions: amplitude-weighted,
  their share of every second- or fourth-order mean (derivatives included: W is a trig polynomial, Bernstein) is O(1/N) of the windows' share → 0.
- (I4), (Ik): pλ⁴ = 1, so the windows contribute M4_CUE(N)·(G_CUE, R_CUE, full CUE four-point function) ~ M4_ζ·(…) by the law of large numbers over
  the ≍ η²V/N → ∞ independent windows, the background O(1) (its derivative fourth moments are O(kurtosis): band-limited Gaussian blocks): ratios →
  the CUE limits, which are ζ's (leading term of the shifted fourth moment, Hughes–Young 2010 / Conrey et al.: arithmetic factor 1 at shifts → 0).
  Variant with ζ itself in the windows: needs window asymptotics; use windows of length ≫ T^{1−δ₀} (δ₀ = power saving of the shifted fourth moment), not T^{0.9}.
- Zeros: (1 − p)·dens(g) + O(p) + O(#junctions) = dens(g) + o(1) (#junctions ≍ η²V/N). Choices of g: the Gaussian with covariance sinc (explore60, Rice 1/√3);
  Gaussian on long blocks of fraction q and ε-constant elsewhere (q/√3, kurtosis 3/q = O(1)); with (I1) these dead blocks already suffice (reviewer):
  q = (π/2)M1_ζ² gives density (π/(2√3))M1_ζ² ≍ L^{−1/2} → 0, kurtosis 3/q = O(L^{1/2}) = o(L²) (0.184 and 9.4 at M1 = 0.45); the thin-band mixture of §2b
  is lower at finite M1 (0.061) but heavier (kurtosis 51 at M1 = 0.45, 886 vs M4_ζ ≈ 5.1·10⁴ at L = 10³). ∎ (Unconditional; ζ enters only through its
  limits: 𝒮 (explore60), the leading shifted fourth moment, Ingham's M4, Heath-Brown's M1. Densities are a.s. along the realization: ergodic background.)
Why it works (chaos view): ζ's fourth cumulant is not small but divergent (κ₄/σ⁴ = M4 − 3 → ∞), the opposite end of the fourth-moment theorem;
a divergent normalized moment is carried by a vanishing length, so any normalized input of order > 2 (any shifts) is diluted away.
Not covered by Theorem 3: unnormalized fourth-order means at macroscopic shifts (x ≍ L, i.e. fixed α in t), where K_ζ = O(1) and can be < 1
(§1b); the Gaussian background adds ≈ 1 there. §1d covers them with ζ itself, under one value-distribution input.

### 1c′. Theorem 3′ (energy in short CUE windows, bulk copied from ζ's own one-point law; adds Selberg's CLT and M1)
Witness: top = independent CUE(N′) windows carrying energy 1 − ε (λ² = (1 − ε)/p), length fraction p → 0, N′ ≈ √(6p)·L/π → ∞ so that
(1 − ε)²M4_CUE(N′)/p = M4_ζ(V); bulk = Φ(Y), Y a stationary Gaussian with spectrum in [−δ, δ] (δ = δ_V → 0), Φ increasing so that Φ(Y(0)) has
ζ's law of F given |F| < b, b = L^{1/4}/log L (i.e. |Z| < L^{3/4}/log L).
- Unconditional tail facts: bulk energy ≤ b⟨|F|⟩ ≍ 1/log L → 0 (Heath-Brown/HRS k = ½); P(|F| > b) ≤ ⟨|F|⟩/b → 0; for 0 < m < 3/2 the part of ⟨|F|^m⟩
  above b is o(⟨|F|^m⟩), ratio ≈ L^{−(m−3/2)²/4} (HRS's unconditional upper bound at k = ¾ with Heath-Brown's lower bounds; m = 1: L^{0.1875}√log L vs L^{1/4}).
  Windows' share of ⟨|F|^m⟩ (m < 2): p^{1 − 3m/4 + m²/8} → 0; their share of the one-point law: total variation O(p).
- So 𝒮, (I4), (Ik) are carried by the windows as in Thm 3; the one-point law of F below b, hence Selberg's CLT and every ⟨|F|^m⟩ with 0 < m < 3/2
  (M1 included), by the bulk; the windows' own one-point law only affects orders ≥ 3/2 other than the matched kurtosis. Negative orders: the copied law
  copies ζ's values, but finiteness for m near −1 needs simple zeros, so they are not claimed (reviewer).
- Zeros: p (windows) + crossings of the level Φ⁻¹(0) by Y, density ≤ δ/π a.s. (Rice; Y ergodic: continuous spectrum) + junctions (2p/N′ per unit) → 0.
  Φ(Y) is not band-limited, but only its zeros (Y's crossings) and its negligible energy enter. Unconditional.
What the inputs cannot see is the bulk's local spectrum: the bulk may oscillate at any frequency in [0, π] (ζ's does, at frequency ≈ π) without moving
any of them; a criterion must weight the bulk at microscopic shifts by length, not energy (mollified moments do exactly that).

### 1d. Proposition 4 (ζ dilutes itself; conditional part marked)
f_a = Z on E_a = {|F| ≥ a⟨F²⟩^{1/2}}; on each run of E_a^c: the constant ±a⟨F²⟩^{1/2} if Z has equal signs at the run's ends (no zero), else a linear
join (one zero); smoothed at scale → 0 for the derivative means. With a = a_T → 0 (e.g. 1/log L):
- unconditional: ⟨Z²𝟙_{E^c}⟩ ≤ a√L⟨|Z|⟩ ≪ aL^{3/4} = o(L) (Heath-Brown) and the filler adds a²L·|E^c| = o(L): all value (underived) second-order
  means of f_a = ζ's + o(1); ⟨Z(t)²Z(t+α)²𝟙_{t∉E}⟩ ≤ a²L⟨Z²⟩ = o(L²) for every shift α, so the two-point fourth-order means K(α), macroscopic and
  unnormalized included, = ζ's + o(L²); normalized microscopic four-point means move by ≤ a√L⟨Z⁴⟩^{3/4} = o(⟨Z⁴⟩) (Hölder).
  Not covered (reviewer): derivative means (the smoothing of the filler amplifies derivatives) and unnormalized macroscopic four-point means other than K(α).
- conditional (Selberg-CLT heuristic for local maxima, not proved): #components of E_a = o(N(T)), hence f_a has o(N(T)) zeros.
Numbers at fixed a in §2c.


## 2. Tests (scripts in this folder, outputs in logs/; every run under memcap 290 MB + /usr/bin/time -l, peaks in the logs)
Cap incidents: a first DH run was killed by memcap at 1336 MB (oscz's 10 MB Hurwitz temporaries grew RSS; fixed by chunk 16, flat 61 MB); a first cuewit run
reached RSS 432 MB unseen by memcap's 2–3 s sampling (Bluestein FFT buffers at a prime-factor length; fixed by 2²¹ points, 139 MB). Logs hold the reruns.
### 2a. Exact CUE data (cue4.py → logs/cue4.log)
- K_N(π) closed form verified as integers for N = 1, 2, 3, 10, 51, 100 (2, 4, 6, 36, 702, 2601) and against slogdet; K_N(π) = 1/2, 4/9, 3/8, …, 2601/10201 → 1/4.
- G_N(x), x = ½, 1, 2, 4, 8, 16: 0.7139/0.2922/0.0731/0.0184/0.0047/0.0012 (N = 100) → 0.7208/0.3009/0.0752/0.0188/0.0047/0.0012 (N = 400) = ceilaudit7's;
  x²G → 3/π² = 0.3040 (Prop 2 thresholds: G < 1/3, 1/15, 9⁻³, 9⁻⁴ beyond x = 0.95, 2.14, 14.9, 44.7).
### 2b. The M1 floor of the background (lpm1.py → logs/lpm1.log)
Continuum of thin-band Gaussian blocks, KKT solution s(ω) = α/(ω + β) (block amplitude A = 1/s², length density s²/π):
| M1 | 0.7978 | 0.70 | 0.60 | 0.45 | 0.30 | 0.20 | 0.10 |
|---|---|---|---|---|---|---|---|
| floor over 𝒮 + M1 | 0.4916 | 0.2296 | 0.1374 | **0.0609** | 0.0217 | 0.0082 | 0.0017 |
| background kurtosis | 3.0 | 7.2 | 15.8 | 51.2 | 202 | 680 | 4510 |
Discrete J = 32 LP (b4.py columns minus the M4 and R rows): 0.0645 at M1 = 0.45 (continuum 0.0609). With M1_ζ ≈ 0.45(26.9/L)^{1/4}: floor 0.061, 0.026,
0.0066, 0.0017 at L = 26.9, 10², 10³, 10⁴, background kurtosis 51, 158, 886, 4230 vs M4_ζ ≈ L²/2π² = 37, 507, 5.1·10⁴, 5.1·10⁶ (ratio → 0, as Thm 3 needs;
at T = 3·10¹² itself the background is not yet dominated). Note M1 = √(2/π) already gives 0.49 < 1/√3: thin-band blocks beat the Gaussian (mean ω/π < rms ω/π).
### 2c. ζ dilutes itself (zdil.py; logs/zdil_zeta.log, 3000 zeros per height, 16 points per unit; Riemann–Siegel vs mpmath siegelz 1–5·10⁻⁶, logs/rschk.log)
f_a of Prop 4 (constant filler), a = τ/rms. "kept" = zeros(f_a)/zeros(Z); dC = max_{u ≤ 8}|⟨f f_u⟩ − ⟨ZZ_u⟩|/⟨Z²⟩; dG = max_{x ∈ {½,1,2,4}}|G_f − G_Z|.
| t₀ (L) | M4 | G(4) | K(4) | a = 0.1: kept / dC / dG | a = 0.2: kept / dC / dG | a = 0.3: kept / dC / dG |
|---|---|---|---|---|---|---|
| first 100 zeros (t ≥ 12, Euler–Maclaurin) | 3.1 | 0.417 | 1.302 | 1.000 / 0.0026 / 0.0000 | 0.880 / 0.024 / 0.0012 | 0.860 / 0.050 / 0.0036 |
| 10⁴ (7.4) | 8.2 | 0.043 | 0.354 | 0.869 / 0.0065 / 0.0001 | 0.696 / 0.033 / 0.0012 | 0.566 / 0.082 / 0.0044 |
| 10⁵ (9.7) | 11.5 | 0.030 | 0.342 | 0.835 / 0.0077 / 0.0001 | 0.638 / 0.035 / 0.0011 | 0.515 / 0.084 / 0.0035 |
| 10⁶ (12.0) | 16.4 | 0.021 | 0.341 | 0.806 / 0.0074 / 0.0001 | 0.594 / 0.037 / 0.0007 | 0.467 / 0.086 / 0.0023 |
| 10⁷ (14.3) | 17.9 | 0.022 | 0.388 | 0.779 / 0.0079 / 0.0001 | 0.571 / 0.037 / 0.0007 | 0.447 / 0.085 / 0.0025 |
| 10⁸ (16.6) | 24.1 | 0.019 | 0.455 | 0.734 / 0.0082 / 0.0001 | 0.550 / 0.036 / 0.0006 | 0.434 / 0.081 / 0.0017 |
| 10⁹ (18.9) | 25.6 | 0.021 | 0.530 | 0.720 / 0.0080 / 0.0001 | 0.526 / 0.037 / 0.0005 | 0.403 / 0.086 / 0.0018 |
| 10¹⁰ (21.2) | 29.2 | 0.029 | 0.844 | 0.700 / 0.0086 / 0.0000 | **0.503** / 0.037 / 0.0003 | 0.383 / 0.086 / 0.0011 |
(10⁹, 10¹⁰: logs/zdil_zeta_hi.log, 257 s, peak 60 MB.) At fixed perturbation the kept fraction falls steadily with T, ≈ −0.18 per unit of log L
(Selberg drift: the moments' carriers {|Z| ≳ rms} shrink); at 10¹⁰ ζ moved by 3.7% in its second-order means, by 1.3% in M4 and by ≤ 0.0003 in G
keeps 50.3% < 1/√3 of its zeros; E_{0.2} (54.9% of the length) carries 99.53% of ⟨Z²⟩ and 100.00% of ⟨Z⁴⟩.


### 2d. Theorem 3's witness at bounded scale (cuewit.py → logs/cuewit.log; 131072 units, 2²¹ points, peak RSS 139 MB)
Gaussian background (exactly flat FFT spectrum on [−π, π]) + K = 295 independent CUE(40) windows, η = 0.300, p = 0.0900, λ = 1.826, cosine junctions:
- zeros per unit 0.6158, predicted (1 − p)/√3 + p = 0.6154; second-order means C(u)/⟨f²⟩ = 0.8970, 0.6261, −0.0159, −0.2116, 0.1280 at u = ¼, ½, 1, 3/2, 5/2
  (sinc: 0.9003, 0.6366, 0, −0.2122, 0.1273).
- Fourth order holds in expectation only: the ensemble formula gives G_f = 0.6995, 0.2764, 0.0721, 0.0211 at x = ½, 1, 2, 4 (exact G_N = 0.7005, 0.2758, 0.0694, 0.0178),
  but the empirical window moments are heavy-tailed (⟨W⁴⟩ per window: median 1.8, max 8408, mean 62.7 against 154.2; moments of moments grow like N⁷),
  so the empirical G_f (0.650, 0.216, 0.081, 0.044 here) needs K ≫ N⁷ windows; Thm 3 has K ≍ η²V/N ≍ T/L ≫ L⁷ (Chebyshev along V = 2^j).
  ζ's own finite-height G (§2c) carries the same sampling caveat.

## 3. Adversarial test (does any of this tell ζ apart from Davenport–Heilbronn?)
- The barrier is about the inputs, and they do not separate ζ from DH: DH (mod 5, κ as explore60) is c₁L(s,χ) + c₂L(s,χ̄) with |c₁| = |c₂|; its fourth
  moment's L⁴ term is |c₁|⁴⟨|L_χ|⁴⟩ + |c₂|⁴⟨|L_χ̄|⁴⟩ (cross terms O(L²)), each with the CUE shape, so DH's normalized G → G_CUE = G_ζ and Theorem 3's
  witness serves DH unchanged; explore60 already has DH in 𝒮 (|a_n|² averaged). At finite height they differ through the cross terms:
  DH at t₀ = 10⁴ (1357 zeros on the line in 1500 units: 0.905 per unit, off-line zeros; logs/zdil_dh.log, oscz copy, |Im|/|Z| ≤ 1.5·10⁻¹¹, peak 44 MB):
  M4 = 4.52 (ζ 8.25), G(½, 1, 2, 4) = 0.527, 0.186, 0.130, 0.203 (ζ 0.582, 0.172, 0.062, 0.043), K(4) = 0.915.
- The dilution acts on DH the same way, only later (smaller kurtosis): a = 0.2 keeps 88.4% of DH's zeros (dC 0.021, dG 0.0010), a = 0.3 keeps 81.2%
  (ζ at the same height: 69.6%, 56.6%); DH's first 96 line zeros (logs/zdil_first100.log): 93.8%, 86.5% (ζ's first 100: 88%, 86%). So nothing here is a statement about ζ's zeros; it is a statement about what the inputs can see.
- Lemma 1 and Prop 2(i, iii) hold for every process with the stated structure; Prop 2(ii) is about ζ/CUE (and DH) data only.

## 4. Where the transfer breaks, and what is left outside the barrier
- The Malliavin–Stein machinery (fourth-moment theorems, Stein bounds on Kac–Rice densities) needs a fixed chaos or bounded κ₄/σ⁴; ζ sits at the other
  end (κ₄/σ⁴ → ∞, G → 0, K < 1 macroscopically): Prop 2 shows no finite-chaos model even reproduces ζ's fourth order, so no chaos-based lower bound on
  zero density can be fed with ζ's moments. The surviving content of the tool is diagnostic: divergent normalized moments ⇒ vanishing carriers ⇒ dilution.
- Inputs not covered (possible escapes): (1) twisted/mollified moments ⟨|ζM|²⟩ (Levinson/Conrey; they flatten the amplitude so the moment's carrier
  has positive length and the bulk's local spectrum becomes visible; the 67.36% record lives here); (2) two-point length-weighted statistics at
  microscopic shifts, e.g. ⟨sgn Z(t)·sgn Z(t+u)⟩ or ⟨(log|Z(t+u)| − log|Z(t)|)²⟩ (they see the bulk's oscillation, i.e. essentially the zeros), and values
  off the line (⟨|ζ(½+α+it)|^{−2k}⟩, 2610.27591 under RH); (3) shifted non-polynomial means such as ⟨|Z(t)||Z(t+x)|⟩/⟨Z²⟩: Thm 3′'s windows match them iff
  ζ's limits equal the CUE's (conjectural for non-integer exponents); (4) macroscopic unnormalized fourth-order means (K(α), α fixed): only Prop 4
  (conditional) covers them; (5) arithmetic inputs (explicit formula, argument principle). Covered now: Selberg's CLT and one-point moments of order
  0 < m < 3/2 (Thm 3′), so a criterion has to weight the bulk at microscopic shifts by length.
- Prop 4's open input, exactly: zeros(f_a) ≤ 1 + #{local maxima of |F| above a}, so it needs "the proportion of the maxima of |Z| between consecutive
  zeros that exceed a·√L tends to 0" for some a = a_T → 0 (a Selberg-type CLT for log|Z| at its relative extrema; Conrey–Ghosh 1985 give only
  the mean of |ζ|² there; Hejhal's CLT for log|ζ′(ρ)| assumes RH and a spacing hypothesis). Not proved here; §2c measures it.


## 5. Novelty (searched 20 Nov 2026)
- arXiv API (axq.py): ti:"squared chaotic" (→ 1503.02154 only); "Wiener chaos" + squares + inequality; zeros + "Wiener chaos" + Kac-Rice (2609.06754,
  1808.02967: random polynomials); "proportion of zeros" + mollifier + necessary/limitation; "critical line" + "large values" + zeros + moments (2605.03665);
  Selberg + CLT + weighted (2507.04150, seen); "multiplicative chaos" + zeros + "critical line" (2506.16115: Dirichlet L-functions as GMC, no zero
  counting); "fourth moment theorem" + zeta: none. Monthly listings 2026-10/11 (lst.py, 7 cats, regex on chaos/Malliavin/fourth moment/mollif/large values/
  critical line/Kac–Rice/zero density/multiplicative chaos/Selberg: 31 titles; closest 2610.14901 (Gibbs measure e^{βH} of a log|ζ| model:
  the tilted measure that carries the moments, no zero count), 2610.27591 (negative moments under RH: an escape of §4)). Zenodo: nothing relevant.
- Status. MNPP (2015) is the tool; Prop 2(i) is its immediate corollary; K_N(π) = (⌈N/2⌉+1)(⌊N/2⌋+1) is Rains' theorem (eigenvalues of U² = two
  independent CUEs of sizes ⌈N/2⌉, ⌊N/2⌋) in disguise. Lemma 1 is an elementary spectral argument (not found stated; it refutes our own heuristic).
  "Moments are carried by the large values" is classical (Keating–Snaith, Soundararajan, Heap–Radziwiłł–Soundararajan, GMC), and that
  unmollified methods give only ≫ T zeros is folklore (Hardy–Littlewood vs Selberg); the method ceiling of Thm 3 for our class 𝒮 + fourth order was
  not found stated. Verdict: a barrier/NOGO entry (routine from classical facts), not a discovery.


## Review (two reviewer runs, 20 Nov)
- Run 1 (§1a–1d): Lemma 1 SOUND (phase typo −arg c_a: fixed). Prop 2 GAP: Rice gives the expected density, a.s. needs ergodic factors (line-spectrum
  counterexample, not in 𝒮): wording fixed. Thm 3 GAP: ζ windows of length T^{0.9} not covered by Hughes–Young's error; junctions need amplitude-weighted
  bounds; derivative fourth moments of the background: fixed by switching to independent CUE windows (ζ windows only as a variant with windows ≫ T^{1−δ₀}),
  amplitude-weighted O(1/N) junction shares, band-limited background blocks. Prop 4 GAP: bounds cover value second-order means and K(α), not derivative
  or general macroscopic four-point means: claim narrowed.
- Run 2 (revised Thm 3): SOUND "for the existential conclusions, choosing dilution sufficiently slowly" (K ≫ N⁷, summable concentration errors;
  κ_g = o(L²)); added its remark that dead blocks with q = (π/2)M1_ζ² already handle (I1).
- Run 3 (Thm 3′): GAP only on negative one-point moments (finiteness needs simple zeros): claim restricted to 0 < m < 3/2; tail exponent
  L^{−(m−3/2)²/4}, windows' TV share O(p), a.s. Rice density by ergodicity of Y and "Φ(Y) is not band-limited" added; no constraint ties the bulk's spectrum to the inputs.

## Suggested NOGO line
- explore62 (Orca, exploration: Malliavin calculus — Malicet–Nourdin–Peccati–Poly squared-chaos inequality, fourth-moment theorem, Kac–Rice in finite chaos → N8 critical-line ceiling) BARRIER: ceilaudit7's open escape (shift-dependent fourth-order means) is closed. Its mechanism X = A·Y (A ≥ 0 independent) is impossible under an exact sinc spectrum (edge-atom lemma; "true to first order" false); finite-chaos models have K ≥ 1 and G ≥ 9^{−q} while ζ/CUE have G ≈ 3/(π²x²) and K_N(π) → 1/4 (exact; Rains), so chaos witnesses cannot match ζ (products of independent Gaussians: density ≥ 1/√3); dilution wins: independent CUE windows of vanishing length carrying the energy plus a background match 𝒮, every normalized microscopic fourth-order mean (G, R, full four-point function, derivatives) and the kurtosis, with zero density 1/√3 (Gaussian background), q/√3 (dead blocks), → 0 with a low-frequency bulk Φ(Y) copying ζ's one-point law below |Z| = L^{3/4}/log L, which adds Selberg's CLT and all one-point moments of order 0 < m < 3/2 (M1 included; HRS tail bounds; unconditional; negative orders only where finite); bounded-scale witness 0.6158 zeros per unit vs 0.6154 predicted; ζ dilutes itself: Z cut to {|Z| ≥ 0.2 rms} (second-order means ±3.7%, M4 −1.3%, G ±0.0003) keeps 69.6% → 50.3% of its zeros from t = 10⁴ to 10¹⁰ (DH 88.4% at 10⁴); reviewer: SOUND after fixes (Thm 3); escapes: mollified moments and other length-weighted microscopic two-point inputs, shifted non-polynomial means, macroscopic K(α), off-line values, arithmetic inputs [ExploreMalliavin](Riemann/Deep/ExploreMalliavin/REPORT.md)

RH relevance: N8 (critical-line proportion; TREE.md record 67.3639%, mollified, unconditional; explore60's ceiling for 𝒮 alone 57.735%; this barrier:
0% for 𝒮 + fourth order + one-point law, not conditional on RH). Unconditional barrier (ζ enters only through explore60's 𝒮-membership, the leading term of its shifted
fourth moment (= CUE), Ingham's M4, Selberg's CLT and the Heath-Brown/HRS fractional-moment bounds): no argument whose inputs on Z are its
second-order means, its normalized microscopic fourth-order means, its kurtosis and its one-point law proves a positive proportion of zeros on the line;
a proof must weight ζ's bulk at microscopic shifts by length, as mollified moments do (Hardy–Littlewood ≫ T unmollified vs Selberg's positive proportion).

## Files
Riemann/Deep/ExploreMalliavin/: cue4.py (exact CUE K_N, G_N, chaos thresholds), lpm1.py (M1 floor, continuum + J-band LP), zdil.py (self-dilution of ζ
and DH; Riemann–Siegel and Euler–Maclaurin), rschk.py (RS vs mpmath), cuewit.py (Thm 3 witness, bounded scale), oscz_copy.py (explore60's oscz.py, Hurwitz chunk 256 → 16 to stop RSS growth);
logs/: cue4.log, lpm1.log, cuewit.log, zdil_zeta.log (10⁴–10⁸), zdil_zeta_hi.log (10⁹, 10¹⁰), zdil_dh.log, zdil_first100.log, rschk.log (memcap peaks ≤ 74 MB inside).

