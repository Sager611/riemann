NO IDEA: every monotone quantity of the de Bruijn–Newman flow runs only in the smoothing direction (t increasing). The complex-side ones (Σ Im z, Σ (Im z)², zero exactly when all zeros are real) are Littlewood's lemma plus the de Bruijn/Planat contraction, and they work the same way for Davenport–Heilbronn, whose pair 0.8085 + 85.70i lands at t ≈ 0.1764. The real-side ones (Perelman-type entropy and Fisher information, Stam, the entropy power inequality, all theorems since arXiv 2602.15822) give only lower bounds on Λ and miss a planted off-line pair by factors of 10² to 10³. Li–Yau/Hamilton log-concavity gives the Turán inequalities at every t (also at t = −3) and for DH as well.

# Exploration explore7: monotonicity formulas of geometric analysis → RH tree (2026-11-05)

Field: Perelman's entropies (F, W, reduced volume, pseudolocality), Huisken's monotonicity and Colding–Minicozzi entropy, Li–Yau and Hamilton Harnack inequalities, and 2015–2026 monotone functionals for zeros of heat flows.
Target: N2 (Λ) and N6 (Laguerre/Turán/Jensen). Not repeated: Almgren/Poon frequency, doubling index, Carleman, Poisson bound for the derivative lever (explore2, [ExploreUC](../ExploreUC/REPORT.md)); Calogero–Moser dynamics (R62, R71).
Conventions: H_t(z) = ∫₀^∞ e^{tu²}Φ(u) cos(zu) du with H₀(z) = Ξ(z/2)/8 (Polymath units). Then ∂_tH = −H'', and simple zeros move by dz_j/dt = 2Σ_{k≠j} 1/(z_j − z_k).

## Survey (arXiv abs and e-print pages, web_search approach arxiv, 24 queries)
- **Geometric analysis.** In one space dimension each monotonicity formula collapses to a known object:
  - Huisken's Gaussian density and Colding–Minicozzi parabolic frequency (2002.11015; ancient caloric functions 1902.01736) become Poon's frequency, which explore2 already tested.
  - Perelman's F and W are Fisher information and the log-Sobolev functional along the heat flow.
  - Li–Yau and Hamilton estimates (generalized: 2312.16181; backward Harnack: 2508.19680; W-entropy: 1303.1242) are lower bounds on the Hessian of log u, that is, log-concavity statements.
- **Finite free information theory.** Garza-Vargas–Srivastava–Stier, arXiv 2602.15822 (Feb 2026), prove these for real-rooted polynomials:
  - monotonicity of the finite free Fisher information Φ_n and entropy χ_n under differentiation;
  - the Stam inequality 1/Φ(p ⊞ q) ≥ 1/Φ(p) + 1/Φ(q);
  - the entropy power inequality N(p ⊞ q) ≥ N(p) + N(q), with N = e^{2χ};
  - Hermite polynomials are the only equality cases.

  Their Lemma 7 shows the score vector is the root velocity of the heat flow, which is exactly the zero ODE above. Hashemi 2604.11922 adds ℓ^p versions (conjectural).
- **Zeros under heat flow.**
  - Hall–Ho 2202.09660 (heat-flow conjecture: transport of the root density until singularities form).
  - Hall–Ho–Jalowy–Kabluchko 2308.11685 and 2304.06665 (heat flow on entire functions of order ≤ 2; GAF zeros move nearly on straight lines).
  - Steinerberger 1811.04844 and Kiselev–Tan 2012.09080 (root density under differentiation).
- **de Bruijn–Newman.**
  - Rodgers–Tao 1801.05914: Λ ≥ 0, a real-side relaxation argument.
  - Dobner 2005.05142: Λ ≥ 0 for the extended Selberg class, via the heat-flowed Dirichlet series Σ e^{(t/4)log²n} n^{−s}.
  - Polymath 1904.12438; Newman–Wu 1901.06596 (survey).
  - Planat 2609.37164: the collective contraction d(y²)/dτ ≤ −2 − 4y²𝒢 for the highest off-line zero, with 𝒢 ≥ −Im(H'/H) off the zero.
  - Planat–Raymond 2611.06245: Λ < 0.156821 by collision exclusion.
  - Raymond, Zenodo 23752103: Λ ≤ 0.107965 (our N2 record). It uses the same contraction, with probe functionals.
- **Laguerre and Pólya-frequency side.** Sra 2611.02944 proves log-concavity (PF₂) and PF_{1.8·10¹²} for the densities of 1/Ξ, and Laguerre inequalities of every order up to height T.
- **Bridge searches.** Queries "Huisken monotonicity zeros analytic functions", "Li–Yau Laguerre–Pólya", "Harnack zeros entire backward heat", "monotone quantity zeros polynomial heat flow imaginary parts" and "de Bruijn–Newman Lyapunov functional": nothing beyond the papers above.

## Why the field stops short
1. **Direction.** Real-rootedness survives the flow toward larger t (de Bruijn) and is generically lost toward smaller t (collisions). So a quantity that is monotone for all configurations and vanishes exactly when every zero is real must decrease in t. An upper bound on Λ then needs that quantity's value at t = 0, which encodes where the zeros of ξ are (zero-free regions), plus a decay rate. Λ ≤ 0 needs the value 0 at t = 0, which is RH itself.
2. **Real-side functionals are one-sided.** Entropy, Fisher information, Stam and the entropy power inequality are defined and proved only where the zeros are real. Applied between t₀ < t₁, they are necessary conditions for real zeros at t₀. They can refute Λ ≤ t₀ but never prove it, which is the Rodgers–Tao direction.
3. **Scale.** ε-regularity and pseudolocality would need the parabolic regularity time at height x, about (4π/log(x/4π))², to exceed the band 0.108 still to be crossed. That holds only up to x ≈ 5·10¹⁷, and the Gaussian density already exceeds the collision value above γ ≈ 6·10⁵ (explore2). Rodgers–Tao's Λ ≥ 0 is this mismatch.
4. **Log-concavity survives the flow at every t.** Multiplying the kernel by e^{tu²} keeps Φ(√w)e^{tw} log-concave in w = u² for every t. So Li–Yau-type inputs (Turán, PF₂) cannot tell t < Λ from t ≥ Λ.

## Transfers and tests
### T1 (N2): Perelman W-entropy → finite free Fisher information, Stam and entropy power inequality
- **Transfer.** Treat H_t as a polynomial p_t = p ⊞_n (Hermite of variance 2t). If p_{t₀} has real zeros, Stam gives 1/Φ_n(p_{t₁}) ≥ 1/Φ_n(p_{t₀}) + 2(n−1)(t₁ − t₀). Here Φ_n(p) = (1/n)Σ_i((2/(n−1))Σ_{j≠i}1/(α_i − α_j))², and the Hermite term is Φ_n(h_τ) = 1/((n−1)τ), confirmed numerically to 7 digits.
  - Use 1 (refute): a violation, seen in data at t₁, would prove non-real zeros at t₀.
  - Use 2 (lower bound): Λ_p ≥ −1/(2(n−1)Φ_n(p₀)).
- **Test** (`t1_stam.py`). The truncated Xi is p = Π_{j≤50}(z² − (2γ_j)²), n = 100.
  - Φ_n decreases along the flow, 7.05·10⁻⁵ → 6.59·10⁻⁵ for t = 0 … 0.2. Stam's relative slack is 0.8–6.3 % and the EPI's 0.015–0.14 %, so both are nearly tight globally.
  - Lower bound. The true first backward collision is at t = −0.3831, at x = −222.78 (the close pair γ = 111.03, 111.87; the isolated-pair estimate −g²/8 gives −0.3571). Stam gives only Λ_p ≥ −71.6 and the EPI ≥ −153, about 190 times too weak.
  - Detection. Off-line pairs were planted at Y = 0.30 and 0.75 (x = 181.3) and at Y = 0.75 (x = 109.4). From data at t₁ = max(0.108, landing + 0.02):
    - Stam would see them only if Φ_n(p_{t₁}) ≥ 1/(2(n−1)t₁). The ratios are 2.1·10⁻³, 9.4·10⁻³ and 9.6·10⁻³, so it is blind by 100–480×.
    - The EPI ratios are 1.4·10³, 5.3·10² and 5.2·10², where ≤ 1 is needed.
- **DH test.** The inequalities follow from the heat flow plus real zeros alone, so they hold equally for H_t^DH on [Λ_DH, ∞). They use nothing that DH lacks.
- **Numeric sanity.** Passed: monotone and consistent on the model, but slack where it matters.
- **Verdict.** Dies, for two reasons:
  - Direction: it gives lower bounds on Λ only.
  - Scale: it is global over n roots. A local, per-unit-length Stam inequality would be Rodgers–Tao's argument and gives Λ ≥ 0, which is known.

### T2 (N6, N2): Li–Yau/Hamilton Harnack → Laguerre/Turán log-concavity
- **Transfer, three facts.**
  1. On the imaginary axis, w(y, t) = H_t(iy) = ∫e^{tu²}Φ cosh(yu) du is a positive eternal solution of w_t = w_yy (Widder form). Hamilton's ancient Harnack inequality, (log w)_yy ≥ 0, is just Cauchy–Schwarz here. The Turán inequalities b_m² > ((2m−1)/(2m+1))b_{m−1}b_{m+1} for b_m = ∫u^{2m}e^{tu²}Φ (Csordas–Norfolk–Varga) say exactly that b_m/Γ(m+½) is log-concave. That follows from log-concavity of e^{tw}Φ(√w), which holds for every t.
  2. On the real axis, the Laguerre ratio ψ = −(log H)'' satisfies ψ_s = ψ'' + 2(log H)'ψ' − 2ψ² in forward-heat time s = −t. The minimum principle keeps ψ ≥ 0 except across a collision, where ψ jumps from +∞ to −∞. So a Li–Yau argument cannot forbid collisions.
  3. Exact hierarchy (derived here; not found in the searched literature). J = |H_t(x+iy)|² = Σ_k L_k(x,t) y^{2k} solves J_t = (J_yy − J_xx)/2. Hence ∂_t L_k = (k+1)(2k+1) L_{k+1} − ½ ∂_x² L_k.
     - The L_k are the generalized Laguerre expressions, for example L₁ = H'² − HH''.
     - In t this system is backward-parabolic, and positivity of L_k is fed only by L_{k+1}.
- **Test** (`t2_harnack.py`).
  - **Turán ratios for ζ** at t = −3, −1, −0.3, 0, 0.3 and m ≤ 25: all > 1. The minimum is 1.0153 (t = −3, m = 25); at m = 1 the ratio rises from 1.0677 to 1.0756 as t increases. Yet H_t has non-real zeros for every t < 0 (Λ ≥ 0).
  - **DH.**
    - Its kernel Φ_DH(v) = e^{3v}Σ a(n) n e^{−πn²e^{4v}/5} is even (functional equation checked), positive, and log-concave in w = v² (max d²/dw² log = −40.7; ζ: −174.5) on the tested range.
    - Its Turán ratios exceed 1 at every t tested (min 1.0206), with larger margins than ζ (m = 1: 1.168–1.229).
  - **Hierarchy identity:** confirmed to 12 digits for k = 0, 1, 2 on a polynomial with a complex pair.
- **DH test.** Fails: DH has the same log-concavity (PF₂) input.
- **Numeric sanity.** Turán holds at t = −3, where zeros are non-real.
- **Verdict.** Dies: blind both to Λ and to DH.

### T3 (N2): Huisken/Colding–Minicozzi on the complex side — the brief's "monotone quantity that vanishes exactly when all zeros are real"
- **Transfer.** Take S(t) = Σ_{Im z_j>0} Im z_j and V(t) = Σ (Im z_j)². Both are zero exactly when all zeros are real, and Λ = inf{t : S(t) = 0}. Exactly,

  dS/dt = −Σ 1/y_j − 2Σ_j Σ_{k real} y_j/|z_j − x_k|² − 2Σ_{j≠k} (y_j + y_k)/|z_j − z̄_k|² < 0,

  with all sums over upper zeros unless marked. A Gaussian (Huisken) weight localises S but adds flux terms. Three readings:
  1. **One-dimensional form.** For a polynomial and Y above all zeros, S = nY/2 − (1/2π)∫[log|p(x+iY)| − log|p(x)|]dx (Littlewood/Jensen). For Ξ the term nY/2 becomes the exact zero count in a window (Turing's method), so evaluating S(0) on a window is RH verification there.
  2. **Rate.** Cauchy–Schwarz gives d(S²)/dt ≤ −2m² (m = number of upper zeros). This never beats de Bruijn's d(y_max²)/dt ≤ −2 in the worst case of one pair. With the real-zero term it becomes Planat's collective contraction, already the engine of Raymond's 0.107965.
  3. **Direction.** Λ ≤ 0 needs S(0) = 0, that is, RH.
- **Test** (`t1_stam.py` (c), `t3_complexside.py`).
  - **Planted models.** S and V decrease strictly. Landing times are 0.0447, 0.2690 and 0.2735, against Y²/2 = 0.0450 and 0.2813; neighbours help by 0.7–4.4 %.
  - **Littlewood identity.** It returns S = 0.600009 against the true 0.6.
  - **DH test.** The DH off-line zero is z = 171.39870 + 0.61703i in Polymath units (s = 0.808517 + 85.699348i); |H| at it is 3·10⁻²² times its value one unit away. Run through its own heat flow (Gauss–Hermite, 28 vs 40 nodes agree to 3·10⁻¹⁶):
    - S = y falls 0.6170 → 0.2315 by t = 0.15 and lands at t* ≈ 0.1764 (7.3 % before de Bruijn's 0.1904).
    - The same pair planted between ζ's zeros 84.7355 and 87.4253 lands at t* ≈ 0.1836.
    - So the functional behaves alike for DH (Λ_DH > 0.17; numerical, not certified) and for a hypothetical ζ pair.
  - **Correction.** The first run's landing fit used a point after landing. The code is fixed and was re-run; both runs are in outputs.txt.
- **Verdict.** The quantity exists, is monotone and is one-dimensional as hoped. But it is Littlewood's lemma plus the de Bruijn/Planat contraction, it is DH-blind, and its value at t = 0 is RH-strength. Dies.

### Argued only
- **Perelman reduced volume and κ-noncollapsing → lower bounds on zero gaps.** Gaps at height x have no uniform lower bound (Lehmer pairs; arbitrarily small normalised gaps are expected), so there is no noncollapsing at a fixed scale. Dies.
- **Hamilton's matrix Harnack for |H| along vertical lines.** J = |H|² is only a subsolution, J_t = J_yy − 2|H'|², and Jensen's K = ∂_y J obeys K_t = K_yy − 2∂_y|H'|². Both give upper bounds only, the wrong sign for Jensen's criterion. Dies.
- **Hall–Ho transport map.** It describes the bulk root density up to the first singularity and is blind to single collisions. Dies.
- **Sturm/Angenent zero numbers of Re or Im H_t(x₀ + iy) along vertical lines.** This is a genuinely forward-parabolic problem in t, but the zero number on the full line is infinite (the phase of H_t rotates as y grows). No usable count. Not pursued.

## Ranking (RH relevance × chance × gain)
| rank | idea | node | relevance | chance | gain | result |
|---|---|---|---|---|---|---|
| 1 | T3: S(t) = Σ Im z_j, monotone, zero exactly when all zeros are real | N2 | high | 0: = Littlewood + Planat contraction; DH lands the same | none beyond Raymond's scheme | dies |
| 2 | T1: Stam / entropy power inequality / Fisher monotonicity (2602.15822) | N2 | medium (lower bounds only) | 0: blind by 10²–10³×; global | Λ ≥ 0 at best (known) | dies |
| 3 | T2: Li–Yau/Hamilton → Turán, Laguerre hierarchy | N6 | low | 0: holds at t = −3 and for DH | none | dies |
| 4 | reduced volume, Hamilton Harnack along vertical lines, Hall–Ho, vertical zero numbers | — | low | 0 | — | argued dead |

## Suggested NOGO line (orchestrator)
- explore7 (Orca exploration: monotonicity formulas — Perelman entropy, Huisken/Colding–Minicozzi, Li–Yau/Hamilton) NO IDEA: monotone quantities of the de Bruijn–Newman flow run only toward larger t.
  - Complex-side ones (Σ Im z, Σ (Im z)², zero exactly when all zeros are real) are Littlewood's lemma plus the de Bruijn/Planat contraction and act alike on Davenport–Heilbronn (its 0.8085 + 85.70i pair lands at t ≈ 0.1764).
  - Real-side ones (Fisher information, entropy, Stam, entropy power inequality; theorems since 2602.15822) give only lower bounds on Λ and miss a planted pair by 10²–10³×.
  - Li–Yau log-concavity gives the Turán inequalities at every t (also t = −3) and for DH.
  - Filter: a quantity that needs real zeros at the earlier time, or holds for every positive log-concave kernel, cannot bound Λ from above.

  [ExploreMonotone](Riemann/Deep/ExploreMonotone/REPORT.md)

## Side notes for node owners
- **N2.** Λ_DH > 0.17 (numerical): the first DH off-line pair is still non-real at t = 0.17 and lands at t ≈ 0.1764 under DH's own flow. This is a cheap lower bound; DH's zeros with Re s > 1 higher up surely push Λ_DH far above 0.5. Not a discovery.
- **N6.** For ζ the Turán ratios increase slowly in t and tend to 1 as t → −∞ (the Gaussian is the equality case). DH's margins are larger than ζ's.
- **N2/N6.** The Laguerre hierarchy ∂_t L_k = (k+1)(2k+1)L_{k+1} − ½∂_x²L_k is a compact way to state how the flow moves Laguerre positivity between orders. It is routine (it follows from J_t = (J_yy − J_xx)/2), so it belongs at most in COROLLARIES.

## Files
- `t1_stam.py` (6.7 s, 40 MB): truncated Xi from 50 zeros; Fisher information, entropy, Stam, EPI, backward collision, planted pairs, V and S monotonicity. It writes `zeros50.json`.
- `t2_harnack.py` (108 s, 40 MB): ζ and DH kernels (evenness, positivity, log-concavity in u²), Turán ratios for t from −3 to 0.3, the Laguerre hierarchy check.
- `t3_complexside.py` (1258 s at nice 19, 45 MB): Littlewood identity; DH off-line zero under the DH heat flow; the same pair planted in ζ.
- `outputs.txt`: all outputs, with the landing-time correction.
- `zeros50.json`: the first 50 ordinates of ζ's zeros (mpmath).
- Runs were single-threaded, under `perl alarm` ≤ 1700 s, one process at a time. Downloads are deleted.

RH relevance: N2 (de Bruijn–Newman); secondarily N6 (Laguerre/Turán).
