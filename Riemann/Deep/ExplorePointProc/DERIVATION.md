# ExplorePointProc: derivations (explore40, 12 Nov 2026)

Notation as ../ExploreDPP and ../SmallGap: f^(u) = ∫ f(α) e^{−2πiαu} dα; J_θ(f) = f(0) + 2∫_0^θ f(α)α dα; θ0 = 399993/700000 = 4/7 − 10⁻⁵;
unfolded ordinates x = γ log T/2π (mean spacing 1); D(λ) = lower density of consecutive gaps ≤ λ (Lamzouri 2611.14486).
Pair measure of a stationary process of intensity 1 (Palm): μ = Mδ_0 + ν; F = μ^ (Montgomery's F, its δ at 0 = the density part).

## 1. Yamada's condition is vacuous given Montgomery (smallest case)
Integer counts give E[(N−k)(N−k−1)] ≥ 0, i.e. Var N_L ≥ θ(1−θ), θ = frac(L), for the count N_L in a window of length L.
Var N_L = ∫ S(α)(sin πLα/πα)² dα with S = F − δ ≥ 0, S = |α| on [−1, 1] under RH. The [−1, 1] part alone is
A(L) = (2/π²)∫_0^1 sin²(πLα)/α dα = Cin(2πL)/π², Cin(x) = ∫_0^x (1 − cos t)/t dt (increasing).
- L ≤ 1 (k = 0): the condition is E N(N−1) ≥ 0, already implied by ν ≥ 0, M ≥ 1.
- L > 1: A(1) = 0.24699, A(1.445) = 0.278, A(2) = 0.31555; θ(1−θ) > A(1) only for θ ∈ (0.445, 0.555), where A ≥ 0.278 > 1/4.
  So A(L) > θ(1−θ) for every L > 1 (table: logs/yamada_margins.log; closest at L = 1.4: A = 0.27118 vs 0.24).
Conclusion: single-window integrality never constrains F beyond 1; in particular it adds nothing to CCCM's bounds on ∫_1^b F
(2108.09258 Cor. 1.2, Thm 17_25: lower bound 0 for b ≤ 1 + 1/√3 from G_1) or to the Montgomery–Taylor multiplicity bound.
For comparison the GUE variance is A(L) + (1/π²)[1 − cos 2πL + 2πL(π/2 − Si 2πL)] (0.38807 at L = 1.5): far from rigid.

## 2. explore34's small-gap duals are not point processes
Hard-core simple process (gaps ≥ λ), Palm atom μ(λ) := ν({λ}) one-sided. Windows A ∋ λ, −A of length < λ contain ≤ 1 point,
so P_0(points at λ and at −λ) ≥ 2μ(λ) − 1, and each such triple gives a pair at distance 2λ: μ(2λ) ≥ 2μ(λ) − 1 (triangle /
hypermetric inequality; general form E[N_φ(N_φ − 1)] ≥ 0 for integer-valued step functions φ, here φ = 1_A + 1_{−A} − 1_{[−ε,ε]}).
Lattice version (sites kh, μ_k = h(1 + η_k)): μ_p + μ_q − μ_{p+q} ≤ 1, μ_p + μ_q − μ_{|q−p|} ≤ 1, Σ_{k=p}^{p+k0−1} μ_k ≤ 1.
explore34's certified duals (θ = 4/7): μ(λ) = 0.757, μ(2λ) = 0.367 at λ = 0.775: 2μ(λ) − 1 − μ(2λ) = +0.147 (0.156 at 0.77,
0.165 at 0.76): violated. θ = 1 duals satisfy it (−0.20, −0.22). logs/triangle_explore34_duals.log.

## 3. Exact local realizability in the dual LP (float floors)
dual_pattern.py: explore34's lattice LP (maximise t, P(α) ≥ |α|1_θ + t on [0, 1/(2h)], η_k = −1 for k < k0) plus a probability
vector p on hard-core words of length R (ones ≥ k0 sites apart), shift-consistent (P(prefix = u) = P(suffix = u) for every
(R−1)-word u), P(site occupied) = h, and P(X_0 = X_k = 1) = h·μ_k for k0 ≤ k < R. In 1D consistent R-marginals extend to a
stationary (R−1)-Markov process, so this is exactly realizability of the range-(R−1) pair data. t* > 0: witness exists, the
augmented method fails at λ. Crossings λ*(t* = 0), θ = 4/7:
| X | plain | + triangle/packing rows | R = 20 (range 2.14) | R = 30 (range 3.2) |
|---|---|---|---|---|
| 10 | 0.7632 | 0.7617 | 0.7513 | 0.7437 |
| 18 | 0.7722 (explore34) | – | – | 0.7524 |
| 26 | 0.7758 (explore34) | – | – | 0.7562 |
Lattice step at range ≈ 2.2, X = 10: k0 = 7/10/14 → 0.7513/0.7517/0.7516 (k0 = 20 > 900 s). Fit R = 30: λ*(X) ≈ 0.7640 − 0.2031/X
(X = 18 predicted 0.7527, measured 0.7524) against plain 0.7839 − 0.207/X: the floor drops by ≈ 0.020 at every X.
Realizable optimum at λ = 0.75 (X = 18, R = 30, logs/struct_R30_X18_lam0.75.json): Palm gap law λ: 0.707, 2λ: 0.125, 1.607: 0.055,
1.714: 0.113; P(two consecutive λ-gaps) = 0.428; long gaps isolated (an AH-like run structure, now consistent).
These floors are float and mix a restriction (support X) with a relaxation (range R); the theorem below does not use them.

## 4. Theorem 1 (hard-core augmentation of Lamzouri's Thm Main2)
Let 1/2 < θ < 4/7, 0 < λ ≤ 1, λ < U ≤ 3λ. Let f be even, f ≥ 0, Schwartz, f^ ≤ 0 on [U, ∞). Let A: [λ, ∞) → R be bounded,
β ≥ 0, c ∈ R with
 (H) f^(x) + f^(x + y) + A(x) − A(y) + β(1 − x) ≤ c for all x, y ≥ λ.
If J_θ(f) > f^(0) + 2c, then D(λ) > 0.
(A = β = c = 0, U = λ is Thm Main2 itself; its f^ ≤ 0 on [λ, ∞) gives (H).)
Proof. Follow Lamzouri l.401-449 (SmallGap/DERIVATION §1): if #E_λ(T) < δN(T), Z_0 = Z∖E_λ (simple, on the line, consecutive
unfolded ordinates x_1 < … < x_{N_0} with gaps g_i = x_{i+1} − x_i > λ), N_0 ≥ (1 − δ)N; steps (i) and (ii) are unchanged and give
Σ_{γ,γ'∈Z_0} f^((γ−γ')log T/2π)w(γ−γ') ≥ (J_θ(f) + O(δ^c))N (they use only that f is even, ≥ 0 and in C_c², with no condition on
its support: Lemma T′ below).
(iii') Upper bound. Diagonal: N_0 f^(0). Index distance ≥ 3: unfolded distance > 3λ ≥ U, f^ ≤ 0, w > 0: drop. Index distance 1, 2:
w(γ − γ') = 1 − O(u²/log²T) at unfolded distance u, and sup|f^(u)|u² < ∞, so replacing w by 1 costs o(N). Hence
 Σ ≤ N_0 f^(0) + 2Σ_{i ≤ N_0−2}[f^(g_i) + f^(g_i + g_{i+1})] + O(1) + o(N)
   ≤ N_0 f^(0) + 2Σ_i [c − A(g_i) + A(g_{i+1}) − β(1 − g_i)] + o(N)          (by (H) at (x, y) = (g_i, g_{i+1}))
   = N_0(f^(0) + 2c) + 2β(x_{N_0−1} − x_1 − N_0) + O(1) + o(N)                 (A telescopes, A bounded).
x_{N_0−1} − x_1 ≤ T log T/2π = N(1 + O(1/log T)) and N_0 ≥ (1 − δ)N, β ≥ 0: the β-term is ≤ 2β(δ + o(1))N; N_0(f^(0) + 2c) ≤
N(f^(0) + 2c) + δN|f^(0) + 2c|. So J_θ(f) + O(δ^c) ≤ f^(0) + 2c + O(δ) + o(1): contradiction for small δ. □
Normalisation: the criterion is homogeneous; to reach Thm Main2's f^(0) = 1 divide f, A, β and c by f^(0) > 0 together.
Truncation (Lemma T′, as SmallGap Lemma T): h ∈ C_c^∞ even, 0 ≤ h ≤ 1, h = 1 on [−1, 1], h_R = h(·/R) (= 1 on [−θ, θ] for R ≥ 1);
G = λR^(λ·) + (λ/φ)R^(λ·/φ) (Lamzouri's reserve: G ≥ 0, even, C_c², G^(0) = 2); F = f h_R + c_r G, c_r > 0 small, R large. F ≥ 0, even, C_c²;
F^ = f^ − E_R^ + c_r G^ with |E_R^(ξ)| ≤ ε_R ξ^{−6} and G^(ξ) ≤ −κ(λ/ξ)^6 on ξ ≥ λ, so F^ ≤ f^ on [λ, ∞) once ε_R < c_rκλ^6: F^ ≤ 0
on [U, ∞) and (H) holds for F with the same A, β, c; J_θ(F) ≥ J_θ(f) (h_R = 1 on [−θ, θ], G ≥ 0), F^(0) = ∫F ≤ ∫f + 2c_r =
f^(0) + 2c_r (0 ≤ h_R ≤ 1, f ≥ 0). So J_θ(F) > F^(0) + 2c for small c_r.
Why F is admissible although supp F = [−R′, R′] reaches far beyond θ (Lamzouri 2611.14486, TeX l.418–440): Thm Main2 (l.392) asks
for a nonnegative even f ∈ C_c²(R), no support condition, and its proof touches the support only through g = ψ_ε f: l.413–424 take
ψ_ε ∈ C_c^∞((−θ−ε, θ+ε)) even, 0 ≤ ψ_ε ≤ 1, ψ_ε = 1 on [−θ, θ], ε = (4/7 − θ)/2, so g is even, ≥ 0, supported in [−θ−ε, θ+ε];
l.425–429 apply Prop SmallSetZeros to g only: Σ_{Z_0} g^w ≥ (f(0) + 2∫_0^θ f(α)α dα + O(δ^{1/16−7θ/64}))N(T); l.430–436 pass to f by
the identity Σ_{ρ,ρ'∈Z_0} f^(i(ρ−ρ')log T/2π)w(ρ−ρ') = (2/π)∫∫ f(α)|Σ_{ρ∈Z_0} T^{α(ρ−1/2)}/(1 − (ρ − (1/2 + it))²)|² dα dt ≥ the same
with g, as f ≥ g ≥ 0 (the identity needs only f even, integrable, f^ integrable: true for C_c² f); l.437–440 combine. Applied to F
(even, ≥ 0, C_c², h_R = 1 and G ≥ 0 on [−θ, θ]) these lines give Σ_{Z_0} F^w ≥ (J_θ(F) + O(δ^{1/16−7θ/64}))N(T) verbatim; step (iii′)
uses F only through F^ ≤ 0 on [U, ∞) and (H).
Reading: (H) is the 1D ground-state (realizability) dual of "sup over stationary gap processes (g_i ≥ λ, E g = 1) of
E[f^(g_1) + f^(g_1 + g_2)] ≤ c": a hard-core set of density 1 must put pairs at nearest- and next-nearest-neighbour distances,
where the certificates' f^ is strongly negative; Thm Main2 throws this away (bounds the off-diagonal sum by 0).

## 5. Certificates (cert_aug.py; python-flint, exact Q and Arb at 1200 bits)
f(α) = Σ_{n≤21} B_n H_{2n}(√(2π)α/s) e^{−πα²/s²}, s = 11/10, B_n ∈ Q (certs/*.json, B_0 += 1/1000 in the checker, EPS),
f^(u) = s Σ (−1)^n B_n H_{2n}(√(2π)su) e^{−πs²u²}; P(y) = p(y²), Q(v) = q(v²) exact. A piecewise linear on 78 rational nodes
λ + k·dg (dg = (U − λ)/77) and A(x) := A(U) for x ≥ U (constant tail; used in R2, R3); β ∈ Q, β ≥ 0.
1. f > 0: p(0) > 0, Sturm: no root of p in (0, ∞).  2. f^ < 0 on [U, ∞): q(w0) < 0, no root in (w0, ∞), w0 = 2·3.14159265·s²U² < 2πs²U².
3. J_{θ0}(f) = p(0) + (s²/π)Σ_k p_k 2^k k!(1 − e^{−X}Σ_{j≤k}X^j/j!), X = πθ0²/s² (Arb); f^(0) = s q(0).
4. sup of Φ(x, y) = LHS of (H) over [λ, ∞)²: (R1) [λ, U]²: lattice δ = dg/5, Arb upper values of f^ at λ + iδ and 2λ + kδ, max over
   corners + 3M2δ²/8 (bilinear interpolation error; A, β-terms are bilinear-exact on aligned cells); (R2) y ≥ U: A(y) = A(U),
   f^(x + y) < 0: 1D max of f^(x) + A(x) − A(U) + β(1 − x) + M2δ²/8; (R3) x ≥ U: A(x) = A(U), f^(x), f^(x + y) < 0, β(1 − x) ≤
   β(1 − U): A(U) − min A + β(1 − U) (min over the nodes, A piecewise linear). M2 ≥ sup_{[λ,2U]}|f^''| = grid max (step
   1/8000, Arb) + M3/16000, M3 = s(2πs²)^{3/2} Σ|c_n|·0.8163·(2(2n+3))^{3/2} (Cramér |h_n| ≤ 1.0865π^{−1/4} < 0.8163 and
   h_n' = √(n/2)h_{n−1} − √((n+1)/2)h_{n+1}). Corner sums in float64 from Arb upper endpoints, +10⁻¹² slack, and again fully in
   Arb (FULLARB=1): R1 corner maxima −0.347888828296, −0.349620915461, −0.353300378304 (± 5·10⁻¹³), identical to the float pass.
   c_cert = max(R1, R2, R3) + 10⁻¹², R1 binding in all three: R1 = −0.3478888283 + 2.718·10⁻⁵ = −0.3478616532 (λ = 0.768; R2 =
   −0.5640124229, R3 = −0.3481091233), −0.3496209155 + 2.728·10⁻⁵ = −0.3495936313 (0.77; R2 = −0.3497135694, R3 = −0.5688739698),
   −0.3533003783 + 2.736·10⁻⁵ = −0.3532730208 (0.775; R2 = −0.5703300330, R3 = −0.3534744596).
5. margin = J − f^(0) − 2c_cert.
| λ | U | dg | β | M2 ≤ | c_cert | J_{θ0}(f) | f^(0) | margin |
|---|---|---|---|---|---|---|---|---|
| 96/125 = 0.768 | 288/125 = 3λ | 192/9625 | 0.64956 | 4.553 | −0.3478616532 | 0.3068393 | 1.0011 | +0.0014626 |
| 77/100 = 0.77 | 231/100 = 3λ | 1/50 | 0.64792 | 4.547 | −0.3495936313 | 0.3070550 | 1.0011 | +0.0051423 |
| 31/40 = 0.775 | 463/200 | 1/50 | 0.64413 | 4.560 | −0.3532730208 | 0.3083008 | 1.0011 | +0.0137468 |
All PASS, rerun 12 Nov 11:40 UTC after the print-only change to 10 digits: margins identical (logs/cert_check_*.log and
cert_check_*_fullarb.log, ≤ 15 s and ≤ 42 MB each). Independent check (xcheck_aug.py, mpmath 40 digits + float): J by quadrature
0.306839281978626 (λ = 0.768), f^ by direct numerical Fourier transform equal to the Hermite formula at 7 points, max Φ over 10⁶ random
and 15,000 edge points −0.3479745137 ≤ c_cert = −0.3478616532, f^ < 0 on [U, 12], f > 0 on [0, 12] (logs/xcheck_768.log).
Corollary. D(λ) > 0 for every λ ≥ 0.768, unconditionally (inputs as SmallGap: Lamzouri's Prop SmallSetZeros with Conrey's density
θ < 4/7, BGST, Fujii).
Shape (λ = 0.77): f^(λ) = −0.219, min −0.496 at 1.14, f^(1.5) = −0.313, bump +0.011 at 2.135, ≈ 0 at U; J = 0.307, c = −0.350.

## 6. Continuum 2-gap LP (primal_2gap.py, float, sampled; how the certificates were found)
Hermite basis as above (normalised h_{2n}), f ≥ 0 and f^ ≤ 0 sampled, (H) on a dg-grid with A at the nodes, maximise J − 2c, f^(0) = 1.
Plain (A = β = c = 0, U = λ): J(0.7851) = 0.99896 (D = 20), 0.99971 (D = 40). Augmented, U = 2.99λ, dg = 0.02:
| λ | 0.755 | 0.76 | 0.765 | 0.77 | 0.775 | 0.7851 |
|---|---|---|---|---|---|---|
| margin D = 20 | – | −0.0158 | −0.0063 | +0.0029 | +0.0121 | +0.0300 |
| margin D = 30 / 40 | −0.0234 / −0.0233 | −0.0144 / −0.0142 | −0.0051 / – | – | +0.0131 / +0.0134 | – |
Threshold ≈ 0.7678, converged in D. U = 1.99λ (nearest neighbours only): no gain (c = β = 0), so the gain is the next-nearest pair.
θ = 1 (RH, Montgomery's range): plain 0.6046, augmented 0.6042 (D = 20; β = 0, c ≈ +0.003): realizability is nearly inert there.
Next: 3-gap version (U ≤ 4λ, A(x, y) on a 2D grid; memory-bound here) should approach the range-3.2 lattice floor ≈ 0.764.
Hard lower limit for any realizability-augmented pair method: 7/11 (the 7/11-discrete sine DPP is a genuine process, explore34 Thm A).
