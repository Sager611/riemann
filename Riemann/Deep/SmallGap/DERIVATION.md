# SmallGap: derivation

Source: Y. Lamzouri, arXiv 2611.14486v1 (10 Nov 2026), SmallGapsPaper.tex; line numbers below refer to it.
Notation: f^(ξ) = ∫ f(x) e^{−2πixξ} dx; J_θ(f) := f(0) + 2∫_0^θ f(α) α dα; θ0 := 4/7 − 10^{−5} = 399993/700000.

## 1. Reduction (Thm Main2, l.392-400)

If 1/2 < θ < 4/7 is fixed, 0 < λ ≤ 1, and f ∈ C_c²(R) is even, f ≥ 0, f^(0) = 1, f^(u) ≤ 0 for |u| > λ and J_θ(f) > 1, then D(λ) > 0.
Proof skeleton (l.401-449): suppose #E_λ(T) < δN(T) (E_λ: zeros with a neighbouring gap ≤ 2πλ/log T, l.404-406; it contains every off-line and every multiple zero, l.443).
(i) g := ψ_ε f, ψ_ε ∈ C_c^∞(−θ−ε, θ+ε), ψ_ε = 1 on [−θ, θ], ε = (4/7 − θ)/2 (l.413-424); Prop SmallSetZeros (l.331-339) for g on Z_0 = Z∖E_λ gives Σ_{Z_0} g^(…)w ≥ (J_θ(f) + O(δ^{c}))N(T) (l.425-431).
(ii) Σ_{Z_0} f^(…)w = (2/π)∫∫ f(α)|Σ_{Z_0}…|² ≥ the same with g, because f ≥ g ≥ 0 (l.432-436).
(iii) Z_0 is simple, on the line, with gaps > 2πλ/log T, so f^(0) = 1 and f^ ≤ 0 off [−λ, λ] give Σ_{Z_0} f^(…)w ≤ N(T) (l.443-446). Contradiction for small δ (l.448).
Lamzouri's choice (l.452-469): θ = θ0, f = λR^(λ·), R(x) = (sin πx/(πx))²/(1 − x²): J_{θ0} = 1.000100147… at λ = 0.7892.
Reproduced in Arb: 1.0001001476430363596…; this family crosses 1 at λ = 0.7891108 (repro_selberg.py, logs/repro_selberg.log).

## 2. Where θ < 4/7 comes from; can another input be substituted?

Chain: Conrey's density N(σ,T) ≪ T^{1−(8/7−ε)(σ−1/2)} log T, uniform σ ≥ 1/2 (l.221-226) → Lemma ZeroDensity: Σ_{β≥1/2} x^{β−1/2} ≪ N(T) for x ≤ T^{θ'}, θ' < 8/7 (l.227-239) → proof of Prop SmallSetZeros: Hölder with x = T^{2θp²} needs 2θp² < 8/7, i.e. θ < 4(1−ε)²/7 (l.381-389) → Thm Main2 needs θ < 4/7 (l.392). Lemma BGST needs supp ⊂ [−1, 1] (l.322-329).
- Better density exponent: substitutes verbatim. With N(σ,T) ≪ T^{1−c(σ−1/2)} log T uniformly for σ ≥ 1/2, the same lines give Thm Main2 for θ < min(c/2, 1). Only the slope at σ → 1/2+ binds: away from 1/2 classical density bounds already have larger c (e.g. Ingham at σ = 3/4 gives c = 1.6). That slope is c = 2θ_mol for a mollifier of length T^{θ_mol}; Conrey's θ_mol = 4/7 is, to our knowledge, still the longest admissible for ζ. We found no density estimate near σ = 1/2 with slope above 8/7 (arXiv API: mollifier AND zero density, "longer mollifier", zero-density AND critical line; the 2026 density paper 2610.15989 covers σ ∈ [7/10, 4/5] in short strips). Hence no substitution is available today; c = 2 (density hypothesis) would give θ → 1, where our SDP gives λ = 0.60400 at d = 41 (BGMM's 0.6039, Lamzouri l.203).
- 2610.08965 (Knausgård: ≥ 67.353% of zeros simple and on the line): cannot replace Conrey. Prop SmallSetZeros needs the weighted moment Σ_{ρ∈E, β≥1/2} T^{θp(2β−1)} (l.375-389), whose weights grow with β − 1/2. A proportion bounds how many zeros lie off the line (≤ 32.65%), not how far, and E_λ also contains on-line simple zeros with short gaps. Its role (making off-line/multiple zeros harmless) is already played by E_λ (l.443).
- θ itself: any fixed θ < 4/7 is allowed. J grows by 2θf(θ)Δθ ≈ 0.54Δθ, so θ0 → 4/7 gains ≈ 5·10^{−6} in λ (invisible at 4 decimals; not used).

## 3. Truncation: Gaussian-type f suffices (Lemma T)

Lemma T. Let 1/2 < θ < 4/7, 0 < λ ≤ 1, f even, f ≥ 0, Schwartz, f^ ≤ 0 on |ξ| ≥ λ, f^(0) > 0 and J_θ(f) > f^(0). Then some F satisfies Thm Main2 exactly as stated, so D(λ) > 0.
Proof. Scale so f^(0) = 1; put η = J_θ(f) − 1 > 0.
(a) Reserve. R as above, R^(u) = (1 − |u| + sin(2π|u|)/(2π))·1_{|u|≤1} ≥ 0, R^ ∈ C_c² (l.456-462). φ = (1+√5)/2, G(x) := λR^(λx) + (λ/φ)R^(λx/φ). Then G ≥ 0, even, G ∈ C_c², G^(ξ) = R(ξ/λ) + R(φξ/λ), G^(0) = 2. Claim: G^(ξ) ≤ −κ(λ/ξ)^6 for ξ ≥ λ, κ = 16/(81π²φ^8) ≈ 4.26·10^{−4}.
 For |t| ≥ 1, R(t) ≤ 0 and |R(t)| ≥ sin²(πt)/(π²t^4) (R(±1) = 0 by continuity; for |t| > 1 use t²(t² − 1) ≤ t^4). Let y = ξ/λ ≥ 1; then |G^(ξ)| ≥ (sin²πy + sin²πφy)/(π²φ^4y^4). Write M = max(‖y‖, ‖φy‖) (‖·‖: distance to Z), n, m the integers nearest y, φy (n ≥ 1). Then ‖nφ‖ ≤ |nφ − m| ≤ φ‖y‖ + ‖φy‖ ≤ φ²M. For n ≥ 1 and m' nearest nφ: (m' − nφ)(m' + n/φ) = m'² − m'n − n² is a nonzero integer, so ‖nφ‖ ≥ 1/(m' + n/φ) ≥ 1/(√5 n + 1/2) ≥ 1/(3n). With n ≤ 3y/2: M ≥ 2/(9φ²y), and sin²πy + sin²πφy ≥ 4M² ≥ 16/(81φ^4y²). □
(b) Cut-off. h ∈ C_c^∞ even, 0 ≤ h ≤ 1, h = 1 on [−1, 1]; h_R = h(·/R), E_R = f(1 − h_R). By Leibniz, ‖E_R^{(k)}‖_1 ≤ C_k Σ_{j≤k} ∫_{|x|≥R} |f^{(j)}| → 0, so ε_R := max(‖E_R‖_1, (2π)^{−6}‖E_R^{(6)}‖_1) → 0 and |E_R^(ξ)| ≤ ε_R ξ^{−6} for all ξ ≥ λ (λ ≤ 1).
(c) F0 := f h_R + cG with c = η/4 and R so large that ε_R < cκλ^6. F0 is even, ≥ 0, in C_c². For |ξ| ≥ λ: F0^ = f^ − E_R^ + cG^ ≤ 0 + ε_Rξ^{−6} − cκλ^6ξ^{−6} < 0. F0^(0) = 1 − ∫E_R + 2c ∈ (0, 1 + ε_R + 2c]. Since h_R = 1 on [−θ, θ] and G ≥ 0: J_θ(F0) = J_θ(f) + cJ_θ(G) ≥ 1 + η. So J_θ(F0)/F0^(0) ≥ (1 + η)/(1 + 3η/4) > 1, and F := F0/F0^(0) satisfies Thm Main2 (l.392-399). ∎
Remark. Alternatively, Lamzouri's proof needs compact support only for g = ψ_ε f (l.413-431); steps (ii)-(iii) hold for any even f ≥ 0 in L¹ ∩ C² (finite sums over Z_0, all on the line). Lemma T avoids relying on this reading.
Why not perturb f directly: f^ decays like a Gaussian, any compactly supported change of f has a slower-decaying, sign-indefinite transform, so f^ ≤ 0 would fail far out; the reserve G supplies the needed strictly negative tail.

## 4. Certificate scheme (cert.py, cert.json)

Ansatz: f(x) = P(√(2π)x/s) e^{−πx²/s²}, P(u) = Σ_{k=0}^{61} B_k H_{2k}(u) (physicists' Hermite), s = 11/10, B_k ∈ Q (20-digit decimals). Since H_n(√(2π)x)e^{−πx²} has transform (−i)^n H_n(√(2π)ξ)e^{−πξ²}: f^(ξ) = s Q(√(2π)sξ) e^{−πs²ξ²}, Q(v) = Σ(−1)^k B_k H_{2k}(v). P(u) = p(u²), Q(v) = q(v²), deg p = deg q = 61, exact in Q[w].
1. f > 0 on R ⟸ p(0) > 0 and no root of p in [0, ∞): exact Sturm count over Q (primitive integer remainders) = 0; independently, Arb root isolation (flint complex_roots) encloses all 61 roots, the only real one ≈ −0.441.
2. f^ < 0 on |ξ| ≥ λ ⟸ q(w0) < 0 and no root of q in [w0, ∞), w0 ∈ Q, w0 < 2πs²λ² (Arb). Sturm count 0; the only real root of q is ≈ 4.6861382 < w0 ≈ 4.6861393 (the sign change of f^ sits just below λ, as it should at an optimum).
3. f^(0) = s q(0) > 0 (exact). J_{θ0}(f) = p(0) + (s²/π)Σ_n p_n 2^n n! (1 − e^{−X}Σ_{j≤n} X^j/j!), X = πθ0²/s² (incomplete gamma; Arb, 4000 bits).
4. Result: J_{θ0}(f)/f^(0) = 1.0000219019375662111 ± 7·10^{−21} > 1 at λ = 0.7851. Lemma T + Thm Main2: D(0.7851) > 0, hence D(λ) > 0 for all λ ≥ 0.7851 (D is nondecreasing).
Cross-check (xcheck.py, mpmath, 60 digits): f^ by direct quadrature equals the Hermite formula to 1.6·10^{−61} at 7 points; J/f^(0) by quadrature agrees to 20 digits.

## 5. How the certificate was found (sdp.py, make_cert.py)

f = Σ_k a_k ψ_{2k}(x/s), ψ_n the L²-normalised Hermite functions (ψ_n^ = (−i)^nψ_n). Constraints as exact polynomial identities matched in the ψ basis: g = Σ X_ij φ_iφ_j and −ĝ = Σ X'_ij φ_iφ_j + (t² − (sλ)²) Σ Y'_ij φ_iφ_j, φ_i(y) = h_i(√π y), Gram blocks split by parity, PSD; triple products by exact Gauss-Hermite quadrature. Maximise J_θ(f) with f^(0) = 1 (cvxpy + Clarabel, double precision). Even d forces a_d = 0 (leading terms of g and −ĝ have opposite signs), so d is odd.
Rigour margin: the optimum has double zeros; we take 0.98·(optimum at d = 61) + 0.02·(interior point at d = 41 with Gram ⪰ 10^{−6}I), which makes p, −q strictly positive.

Float thresholds λ*(d) (smallest λ with max J > 1, best of s ∈ {1.0, 1.1, 1.2, 1.3}; thresholds.py):

| d | 21 | 31 | 41 | 51 | 61 |
|---|---|---|---|---|---|
| λ*(d) | 0.785901 | 0.785433 | 0.785301 | 0.785172 | 0.785072 |

At λ = 0.7850: max J = 0.99992 (d = 61), 0.99996 (d = 67), 1.000086 (d = 100, solver inaccurate). λ*(d) falls by about 10^{−4} per 10 degrees, so the true extremal value is below 0.78507 and its limit is unknown: we computed no dual (upper) bound for J, hence no floor for this method. Certifying 0.7850 needs d ≳ 80 with better than double precision; with cvxpy that exceeds the 300 MB cap (d = 81: 401 MB).
A compactly supported Selberg-type family (f^ = sin²(πbξ)·Σ_{|n|≤K}[c_n/(bξ − n)² + d_n/(bξ − n)], f ∈ C_c² iff c_0 + 2Σ(c_n + nd_n) = 0; sampled LP, .scratch only) reproduces Lamzouri at K = 1 but reaches only J ≤ 0.9994 at λ = 0.7851 (b ≤ 5, K ≤ 45): worse than the Gaussian family.
