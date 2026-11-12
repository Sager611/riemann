BARRIER: Baernstein's star function, Steiner/circular symmetrisation, Pólya–Szegő, Essén/Hayman–Kennedy and harmonic-measure comparison hold for every subharmonic function, and the t-independent convex u*(σ) = s*·L·max(0, ½ − σ) satisfies all of them together with the line second moments of VΨ while meeting Levinson's bound with equality, so none of them can beat the concavity step (node N8); the only extra structure, zeros as 2π-atoms, buys at most 17% of the bound at T = 10⁴ (an explicit analytic competitor with the same 42 line second moments has right-zero fraction 0.200 against the smeared optimum 0.2407).

# explore45 — symmetrisation and the star function → N8 (Levinson's concavity step)

## Transfer
- Field: potential theory/geometric function theory. Tool: Baernstein's star function u*(σ, θ) = sup_{|E|=θ} ∫_E u(σ+it) dt, subharmonic for subharmonic u; its integral-means inequalities ∫Φ(u) ≤ ∫Φ(v) (Φ convex increasing, v the symmetrised comparison), and the other rearrangement/harmonic-measure inequalities.
- Node N8: Levinson–Conrey: Littlewood on [a, ∞) × [T, 2T], a = ½ − R/L, then ∫ log|VΨ(a+it)| dt ≤ (T/2) log((1/T)∫|VΨ|²) (Jensen), giving κ ≥ 1 − (1/R) log c(P, Q, R). Bridge tested: does subharmonicity of u = log|VΨ| in (σ, t), with the mean squares C(σ) known on a whole range of lines, give a smaller ∫ log|VΨ| than Jensen?

## Statement (rescaled s = ½ + x/L + it; m(x) = mean of u, h(x) = ½ log C(x), ν(x) = fraction of zeros with Re > x)
Theorem B (routine, proof in DERIVATION D1). Constraints: Littlewood (m convex, −m′ = ν, κ_n-normalised), m(+∞) = 0 (leading coefficient 1), m ≤ h on every line. Then
 sup ν(0) = inf_{x<0} h(x)/|x| (Levinson's bound, at its best line),
attained by m*(x) = s*·max(0, −x). The function u*(x, t) := m*(x) is subharmonic (convex in x), has mean m* and mean e^{2u*} = e^{2m*} ≤ C on every line, and its Riesz measure is (s*/2π) dt on x = 0. Every inequality valid for all subharmonic functions in the strip (star function, symmetrisation, Pólya–Szegő, Essén/Hayman–Kennedy, Carleman, harmonic measure, Hadamard–Doetsch three lines) therefore holds for u*, so no combination of them with line second moments beats Levinson. Rearrangement in t also preserves every line mean (the only quantity Littlewood reads), and the star function bounds integral means of convex Φ(u) from above, i.e. moments from zero data, the reverse of the direction the method needs.
What separates log|f| from u*: the Riesz measure of log|f| is a sum of 2π-atoms (zeros) of bounded total density. Test of how much that is worth: T3.

## Tests (all single-threaded, memcap 290 MB; peaks below)
T1 (smallest case). Conrey's formula, P(x) = x, Q(x) = 1 − bx: θ = ½: b = 1.03217, R = 1.24618, κ = 0.347439 (Levinson's 0.3474); θ = 4/7: κ = 0.391130 (conrey_c.py).
T2 (actual VΨ, these parameters, θ = ½, y = 100, T ∈ [10⁴, 1.4·10⁴], dt = 0.01, 4806 zeros, Euler–Maclaurin checked against Arb to 10⁻¹⁰; levin.py, ana.py; peak 207 MB):
| quantity | ζ | DH |
|---|---|---|
| zeros of V right of ½ (winding on [½, 3]) | 0 | 275 (0.057) |
| finite-T Levinson bound ν(0) ≤ h(−R)/(Rκ_n) | 0.2491 (κ-equiv. 0.502) | 0.534 |
| Jensen/concavity loss (h − m)(−R)/(Rκ_n) | 0.1437 | 0.197 |
| Littlewood-weight loss (zeros in (a, ½)) | 0.1056 | 0.261 |
| m_f(0) term (right zeros weighted beyond ½) | −0.0001 | 0.020 |
| Jensen gap J(−R), Var log|f|, J/Var | 0.179, 0.233, 0.77 | 0.245, 0.333, 0.73 |
| Jensen gap without mollifier J_V(−R) | 0.734 | 0.290 |
LP over all 23 lines (only h used): ζ 0.2426 at x = −1 = best single line (Theorem B). Under RH and b < 2 the true ν_V(0) is 0, so all 0.249 is method loss, 58% of it the concavity step. DH is normalised by ζ's count (its own L would be log(5T/2π)).
T3 (quantised inverse problem; invlev.py). Competitor f(z) = p(e^{−z/P_z}), p = Π_{j≤K}(1 − w/w_j), c₀ = 1, K = 1.61 P_z (the zero density of VΨ, 1 + log y/L), period 2πP_z; mean|f|² on line x = Σ|c_k|² e^{−2kx/P_z} (exact on 2K+2 nodes). Maximise M = #{|w_j| ≤ 1} subject to mean|f|² ≤ measured C on every line. Max violation of log S − log C (≤ 0 = feasible):
| profile | P_z | feasible M (ν) | infeasible M: violation | smeared LP |
|---|---|---|---|---|
| 23 lines, x ∈ [−4, 18.9] | 30 | 6 (0.200) | 7: +0.024 | 7.27 (0.2424) |
| 23 lines | 60 | 12 (0.200) | 14: +0.024 | 14.5 |
| 42 lines, x ∈ [−12, 18.9] | 60 | 8, 10, 12 (0.200) | 13: +0.0049, 14: +0.028 | 14.4 (0.2407) |
| 42 lines | 100 | — | 21: +0.0011 (unconverged), 22: +0.0092, 23: +0.023 | 24.1 |
Peaks ≤ 84 MB. Feasible competitors put all M right zeros on the line and a screening layer at x ≈ −2.7 (distance 1.7 from the tight line x = −1).
T4 (independent check, verify.py): the P_z = 60, M = 12 competitor rebuilt by np.poly, 12 roots of np.roots in |w| ≤ 1, max(h_comp − h) = −1.5·10⁻⁴ on all 42 lines; m_comp = 0.200|x| on [−1.4, 0]; Jensen gap 0.030 at x = −1 (true VΨ 0.164), 0.087 at x = 0.
Reading: rigorous, no argument from these 42 line second moments + analyticity + leading coefficient 1 + bandwidth 1.61 proves ν < 0.200 (κ-equivalent > 0.600 at this T; Levinson 0.519); potential-theoretic inequalities stop at 0.2407. Numerical, not proved: the quantised optimum lies in [0.200, ~0.21], strictly below 0.2407.

## Adversarial
Theorem B and T3 use only the profile C(x), so they are blind to ζ versus Davenport–Heilbronn and to GUE/Poisson models with the same density: no RH information, as expected of a barrier. T2 shows the concavity loss is of the same size for ζ (0.144) and DH (0.197); log|VΨ| is not Gaussian (J/Var = 0.77, a complex Gaussian would give J = γ/2 = 0.289 regardless).

## Where it breaks / what could still move N8
- Breaks: all symmetrisation-type inequalities are one-sided upper bounds valid for every subharmonic u; the Levinson LP extremal is subharmonic, so they add nothing (Theorem B).
- Remaining room: the 2π-atomicity of zeros of a bounded-type function. A "quantised Levinson inequality" (sup ν over exponential-type f with c₀ = 1 and line moments ≤ C) lies between 0.200 and 0.2407 here; its proof would need a Jensen-gap lower bound from zeros near a line that survives screening by zeros at comparable distance. Not even the record method: N8's record (0.6736) is TEN-type; the natural target is the certified one-piece GL(2) Levinson–Conrey ceiling 0.07005633 (gl2ceil), which is a ceiling of the Jensen functional only.

## Novelty
arXiv (axq.py): Levinson+concavity, Levinson+arithmetic-geometric, Baernstein+zeta, star function+zeros, Levinson's method+optimal, mollifier+Littlewood+Jensen, long mollifiers, subharmonic+zeta zeros, integral means+Dirichlet series, Mahler measure+Dirichlet polynomial, Littlewood's lemma+critical line: no paper with Theorem B or the inverse problem; hits 2511.06109, 1403.5786 (mollified-function choice, Feng), 1802.10521 (five-twelfths), 2508.11108 (seen) concern mollifier choice or length. Baernstein (Acta 1974; Symmetrization in Analysis, 2019), Essén, Hayman–Kennedy: inequalities for all subharmonic functions, none for Levinson. Theorem B is routine (LP duality); T2's loss decomposition and T3's competitor are measurements, not discoveries.

## Files
code/ (conrey_c.py, levin.py, check.py, ana.py, mkprof.py, invlev.py, verify.py, speed.py; the .npz intermediates were not kept), logs/ip100.log; derivations Riemann/Deep/ExploreSymm/DERIVATION.md (D1–D3).

## Suggested NOGO line
- explore45 (Orca, exploration: symmetrization and the star function — Baernstein, Steiner/Pólya–Szegő, Essén/Hayman–Kennedy, harmonic measure → N8, N7) BARRIER: these inequalities hold for every subharmonic function and the t-independent convex u* = s*·max(0, −x) meets Levinson's bound with the same line second moments, so they cannot beat the concavity step; at T = 10⁴ (θ = ½) the concavity step is 0.144 of Levinson's 0.249 (ζ; DH 0.197), and an analytic competitor with the same 42 line mean squares has right-zero fraction 0.200 (smeared LP 0.2407), so zero-atomicity buys ≤ 17% [ExploreSymm](Riemann/Deep/ExploreSymm/REPORT.md)

RH relevance: N8 (Levinson concavity step; side: N7 zero-density detection uses the same Littlewood + Jensen step).
