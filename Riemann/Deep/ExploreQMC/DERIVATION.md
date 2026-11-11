# explore30 — derivation: Weil's zero sum as a quadrature with nodes at the zero ordinates (node N1)

## 1. Setting
f ∈ C_c^∞(−a, a) (log scale; window c = e^{2a}), F(z) = ∫ f(x) e^{izx} dx, f^*(x) = conj f(−x), g = f * f^*:
supp g ⊂ [−2a, 2a], g(0) = ‖f‖², conj g(−y) = g(y). Weil: W(f) = Σ_ρ ĝ(γ_ρ), γ_ρ = (ρ − 1/2)/i, ĝ(z) = F(z)·conj F(z̄);
RH ⟺ W(f) ≥ 0 for all a and f. Window records: W ≥ ε_c‖f‖² for a = 1.925 (c = 47, all f), 2.131 (c = 71, even f),
claim 2.303 (c ≈ 100).

## 2. Lemma 1 (quadrature form)
For ρ = 1/2 + δ + iγ: ĝ(γ − iδ) = ∫ g(y) e^{δy} e^{iγy} dy. With its partner 1 − ρ̄ (γ + iδ) the pair gives 2K(γ, δ), where
  K(t, δ) := ∫ g(y) cosh(δy) e^{ity} dy = Re[F(t − iδ) conj F(t + iδ)] = |(f cosh δx)^(t)|² − |(f sinh δx)^(t)|²   (real).
Hence  W(f) = Σ_ρ K(Im ρ, Re ρ − 1/2)  (each zero once; K(t, 0) = |F(t)|²).

## 3. Lemma 2 (exactness)
∫_R K(t, δ) dt = 2π g(0) cosh 0 = 2π‖f‖² for every δ, and K(·, δ) has Fourier support in [−2a, 2a] (bandlimited).
So every node kernel has the same integral: moving a zero off the line does not change the "integral", only the node.

## 4. Lemma 3 (Poisson exactness) and the Koksma–Hlawka split
Nodes t_k = t_0 + kΔ with a common δ: Σ_k K(t_k, δ) = (2π/Δ) Σ_m φ(2πm/Δ) e^{2πimt_0/Δ}, φ = g·cosh(δ·). If Δ < π/a only m = 0
survives: Σ_k K = (2π/Δ)‖f‖², the value of on-line zeros of the same density, for every δ. Regular off-line lattices finer than
π/a are invisible to the window. In general
  W(f) = ∫ K(t, 0) dN̄(t) + ∫ K(t, 0) dS(t) + Σ_ρ [K(γ_ρ, δ_ρ) − K(γ_ρ, 0)],   ∫ [K(·, δ) − K(·, 0)] dt = 0,
so off-line zeros act only through the discrepancy (irregularity) of the ordinates; Koksma–Hlawka gives
|∫K dS| ≤ V(K)·sup|S|, too weak to certify (V(K)·sup|S| has the order of the main term D·2π‖f‖²).

## 5. The adversary
Two-bump f = φ((x − a′)/w) + e^{iθ} φ((x + a′)/w), φ(s) = (1 − s²)^8 on [−1, 1], a′ = a − w (supp f = [−a, a]), modulated to
height t_0 (u = t − t_0). For narrow bumps K(t_0 + u, δ) ≈ Ψ(u)[2 + 2cosh(2δa′) cos(2ua′ − θ)], Ψ = |wΦ(wu)|² ≥ 0.
K < 0 exactly where cos(2ua′ − θ) < −sech(2δa′): a fraction p = arccos(sech 2δa′)/π of each period P = π/a′.
Block [t_0 − L, t_0 + L], L = 40/w: every zero of the block is moved into these intervals as an off-line pair with common δ;
counts per period ≈ D·P (D = log(t/2π)/2π), rounded by a second-order error diffusion that keeps S and ∫S near 0
(cluster shift ≤ (1 − q)/2 of the interval). Outside the block: on-line lattice of density D (contributes ≤ 5·10⁻¹⁵ of
D·2π‖f‖²). Then W(f) = Σ_block 2K + Σ_out |F|² < 0. Geometry: gap g = (1 − p)P, count deviation ≈ Dg/2,
∫S excursion ≈ DgP/8 = gP log t/(16π).

## 6. Admissibility (unconditional zero-side information)
(i) RH for 0 < γ ≤ T = 3·10¹² (Platt–Trudgian 2021): the block lies above T.
(ii) no zero with β > B_2 = 0.873381966… (Lean-checked, Richey): δ = 0.37.
(iii) |N(t) − (t/2π)log(t/2πe) − 7/8| ≤ min{Trudgian 2014: 0.112 log t + 0.278 log log t + 2.510; Hasanalizade–Shen–Wong 2022:
     0.1038 log t + 0.2573 log log t + 9.3675; arXiv 2606.31441: 0.0946 log t + 5.746} − 1 (1 reserved for the true S at the
     block ends, which are chosen where S_true ≈ 0 and S_true(left) = S_true(right)).
(iv) |S_1(t)| ≤ 0.059 log t + 2.067 (Trudgian 2011) − 1.
(v) zero density: the block has ≤ 6000 zeros (≤ 2965 pairs) with β ≤ 0.87 near 3·10¹²; explicit density bounds
     (Kadiri–Lumley–Ng, N(σ,T) ≤ C(σ)T^{8(1−σ)/3}log^{5−2σ}T + …) allow ≥ 10⁸ there.
(vi) functional equation: pairs (ρ, 1 − ρ̄) at the same height, conjugates at −γ (|u| ≈ 2t_0, negligible).
(vii) on-line proportion and short-interval results are averages over [T, 2T] or intervals ≥ T^θ; a block of length ≤ 1400 is
     invisible to them.

## 7. Thresholds (logs/qmc1, qmc4, qmc5, qmc7)
- T = 3·10¹², all of (i)–(vii): admissible W < 0 for a ≥ 1.7 (c ≥ 30); not found for a ≤ 1.6 (S_1 binds: 2.94 > 2.76).
  Without the S_1 constraint already from a = 1.2 (c = 11).
- t = 10¹⁰⁰: admissible for a ≥ 1.85.
- t → ∞: the needed ratios are |S|/D → g/2 = 0.559 (a = 1.925, w = 0.03a) against 2π·0.0946 = 0.594, and
  |∫S|/log t → gP/(16π) = 0.037 against 0.059: admissible at every height for a ≥ 1.925; narrow-bump limit (w → 0):
  a ≥ 1.78 by the gap formula g(a) = (1 − arccos(sech 2δa)/π)π/a ≤ 4π·0.0946 = 1.189.
- smaller δ at T = 3·10¹² (balanced): δ = 0.2 at a = 1.925, 2.131, 2.303; δ = 0.1 and 0.05 at a = 2.131, 2.303.

## 8. Consequences
(a) For a ≥ a_1 no proof of W ≥ 0 on [−a, a] can use only (i)–(vii): window certificates beyond a_1 must use the prime side
    (Euler product). a_1(3·10¹²) ≤ 1.7; a_1(T) ≤ 1.925 for every T, relative to the constants in (iii), (iv).
(b) Conversely the certified window c = 47 (Chekmachev; re-certified by explore28) excludes such blocks at every height:
    information about high zeros not implied by (i)–(vii). The excluded blocks are exotic (hundreds of consecutive off-line
    zeros in a fixed phase pattern), so this is not RH progress.
(c) Isolated off-line pairs are invisible to localized f: a count-preserving merged pair at height 15000 lowers W only from
    1.071 to 0.898; DH's nine real pairs (δ up to 0.37) keep min W/(2π‖f‖²) ≥ +0.26 for a ≤ 2.303. They are seen only by the
    low-frequency near-null modes (explore19: reach m* ≈ 4.4πc; DH fails at c ≈ 30.7).
(d) Dependence on the constants: asymptotically this adversary needs leading coefficients ≥ g/(4π) ≈ 0.089 (count) and
    ≥ gP/(16π) ≈ 0.037 (S_1) at a = 1.925; better future bounds would remove it only for large T. At T = 3·10¹² the constant
    terms dominate.

## 9. N5 (examined, not pursued)
d_N² = 1 − bᵀG⁻¹b is the squared worst-case error of the optimal rule with nodes 1/k for the functional ⟨·, χ⟩ in the RKHS of
Vasyunin's kernel, i.e. the posterior variance Var(W(1) | ∫ρ_{1/k} dW, k ≤ N) of Brownian motion; the "target measure" is the
Möbius measure on the nodes themselves (χ = −Σμ(k)ρ_{1/k} formally), so smoothness, fill-distance and weight-positivity theorems
do not apply (signed Möbius weights; explore22's c₁ ≥ 0 floor 0.2598), fooling functions are Burnol's orthogonal vectors
(lower bound C/log N, range-limited, so node selection cannot beat it), and the prediction-theory form is Balazard–Saias–Yor.
Every finite-N statement is DH-blind (explore13: DH's floor appears near N ≈ 10³³⁴). nb2nd's κ̄ is a nonlinear inverse-Gram
functional with no quadrature handle.
