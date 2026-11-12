# LostLaw — derivations (lostlaw, Orca, 12 Nov 2026)

Setting = MoriRate (ii′). Even sector. Z ⊂ ℝ symmetric, locally finite (ζ's ordinates under RH, or an RH-consistent model),
Q(u) = Σ_{γ∈Z} m_γ|û(γ)|² on L²(−a, a), a = ½ log μ, û(t) = ∫u e^{itx}. A = Friedrichs operator of Q: compact resolvent,
eigenpairs (E_i, ψ_i), E₁ = λ < E₂ (simple bottom, certified for 2 ≤ μ ≤ 30 and μ = 35, 40 in WeilSimple, CCMFMM). γ₀ ∈ Z simple, e₀ = cos(γ₀x)1_{(−a,a)},
c_i = ⟨e₀, ψ_i⟩, Q′ = Q − 2⟨·, e₀⟩², λ′ = λ_min(Q′), r = λ′/λ, y₀ = γ₀/2πμ, g₂(y) = arccosh(2/|y|) (|y| < 2), 0 (|y| ≥ 2).
Law (MoriRate): log(λ/λ′) = 4g₂(y₀) + o(1).

## 1. Exact reduction to full-form data
Definitions (all of the FULL form Q): θ = 2⟨e₀, A⁻¹e₀⟩ = Σ_i 2c_i²/E_i = sup_u 2⟨u, e₀⟩²/Q(u) ∈ (0, 1];
w₀ = 2c₁²/λ = share of the ground-state energy Q(ψ₁) = λ carried by the samples ±γ₀; ρ = λ Σ_{i≥2} 2c_i²/(E_i(E_i − λ′)).

Proposition 1. Assume Z ∖ {±γ₀} is a uniqueness set for the form domain (true for ζ's ordinates and every model used here; it can fail
for general multisets — reviewer: ordinates ±(n + ½)π with multiplicities (n + 1)², a = 1, lose it after one deletion). If c₁ ≠ 0 then 0 < λ′ < λ and
    1 − θ = r·[w₀/(1 − r) + ρ],   0 ≤ ρ ≤ θλ/(E₂ − λ′) ≤ λ/(E₂ − λ).                                   (1)
Hence
    log(w₀/(1 − θ)) ≤ log(λ/λ′) ≤ log(w₀/(1 − θ)) + log(1 + ρ/w₀) + log(1/(1 − r)).                      (2)
Proof. φ(x) = 2⟨e₀, (A − x)⁻¹e₀⟩ = Σ 2c_i²/(E_i − x) increases on [0, λ) from φ(0) = θ ≤ 1 to +∞ (c₁ ≠ 0), so it has one root of
φ = 1 there, and it is λ′ (Weinstein–Aronszajn; MoriRate Lemma A(c)); λ′ > 0 iff θ < 1, and θ < 1 because Q′ ≥ 0 has trivial
kernel (A′ has compact resolvent, and an eigenfunction with Q′ = 0 would vanish on the uniqueness set Z ∖ {±γ₀}). Subtract: 1 − θ = φ(λ′) − φ(0) = Σ_i 2c_i² λ′/(E_i(E_i − λ′)).
The i = 1 term is λ′·w₀λ/(λ(λ − λ′)) = r w₀/(1 − r); the rest is rρ. ρ ≤ (λ/(E₂ − λ′)) Σ_{i≥2} 2c_i²/E_i ≤ λθ/(E₂ − λ′). (2): from
(1), r = (1 − θ)(1 − r)/(w₀ + (1 − r)ρ); drop (1 − r) ≤ 1 and ρ ≥ 0 for the left inequality, use (1 − r)ρ ≤ ρ for the right. ∎

Corollary 1. Along any sequence with r → 0 and ρ/w₀ → 0 the law is equivalent to w₀/(1 − θ) = e^{4g₂(y₀) + o(1)}: the cost of a
lost sample is the ratio of the ground state's energy share at γ₀ to the extremal share inf_u Q′(u)/Q(u) = 1 − θ.
Numbers (Galerkin, gal.py): ρ/w₀ ≈ 0.011 (gap in (2), μ = 16, 30, first zero); r = λ′/λ ≤ 0.3 for y₀ ≤ 1.9 (log(1/(1 − r)) ≤ 0.36), ≤ 0.05 for y₀ ≤ 1.5;
both factors are exponentially small: ln w₀ = −164.7, ln(1 − θ) = −178.0 at μ = 16, γ₁; ln w₀ = −337.3, ln(1 − θ) = −353.2 at μ = 30.

## 2. MoriRate's sufficient target for the upper half fails at every computed window
MoriRate (ii′) Lemma A(b): λ′ ≥ (1 − θ)λ, so 1 − θ ≥ e^{−4g₂ − o(1)} would give the upper half. By (1), 1 − θ ≤ r(w₀ + ρ)/(1 − r) ≤
r(w₀ + λ/(E₂ − λ))/(1 − r): the target would force w₀ + λ/(E₂ − λ) ≥ e^{−o(1)}, i.e. an order-one ground-state share at a sub-Nyquist
sample or an order-one gap ratio. Certified (theta_rig.py: θ ≥ 2(b·u)²/(uᵀMu) for the Galerkin vector u = M⁻¹b, Arb balls for M and
for b with a certified zero ball): 1 − θ ≤ e^{−X} with X in the table of REPORT.md (X ≈ 2μP(y₀) ≫ 4g₂(y₀)). So Lemma A(b) loses
≈ 2μP(y₀) − 4g₂(y₀) nats (164 at μ = 16, 337 at μ = 30, first zero), growing linearly in μ. Reviewer: finite-μ values do not refute an
estimate with an unspecified o(1); an asymptotic refutation needs 1 − θ(γ₀) ≤ e^{−cμ} for fixed y₀ < 2, which §3's proxy statement suggests
(rate 2P(y₀)) but which is not proved here.

## 3. The rate of 1 − θ: Green potential of the deficit
P(y₀) := ∫_{−2}^{2} g_Ω(t, y₀)(−log|t|) dt, Ω = ℂ ∖ ((−∞, −2] ∪ [2, ∞)), g_Ω its Green function (P(0) = 2π, Mori Thm 4).
Proxy statement (rigorous, Mori's certificate moved from 0 to y₀): for admissible σ (σ ≥ c₀, finite log-moments),
p_{y₀}(σ) := sup_x U_σ(x) − U_σ(y₀) ≥ ∫(U_σ − U_σ(y₀)) dω_{y₀} = ∫ g_Ω(t, y₀) dσ(t) ≥ ∫ g_Ω(t, y₀) dc₀(t) = P(y₀),
ω_{y₀} = harmonic measure of Ω at y₀ (probability on K = ℝ ∖ (−2, 2)), using U_{ω_y} − log|· − y| = g_Ω(·, y) (balayage identity,
MoriRate Lemma B) and Fubini. Spectral check: −log(1 − θ(γ₀)) = 2μP(y₀) + d, d ∈ [−5.3, +1.8] at μ = 16 and 30 for
y₀ ∈ [0.07, 1.51] (−log(1 − θ) from 27.9 to 353.2). So 1 − θ is the spectral counterpart of a virtual node at y₀, and by
Corollary 1, −log w₀ = 2μP(y₀) − 4g₂(y₀) + O(log μ): the ground state's share at a sub-Nyquist sample decays at the same rate.
Not new as a statement about the proxy (it is Mori's Theorem 4(i) argument at another base point); the spectral match is ours.

## 4. What each half needs (consequence of (2))
- Lower half, log(λ/λ′) ≥ 4g₂ − o(1): by (2) it suffices that w₀/(1 − θ) ≥ e^{4g₂ − o(1)}, i.e. a test vector u with
  Q′(u)/Q(u) ≤ w₀ e^{−4g₂ + o(1)} (upper bound for 1 − θ) AND a lower bound for the ground-state share w₀. Or directly a test vector for
  Q′ with Rayleigh quotient ≤ λe^{−4g₂ + o(1)} AND a lower bound for λ. Either way one positivity (certificate) statement about the
  undeleted bottom is needed: test vectors alone give nothing relative to an unknown λ.
- Upper half, log(λ/λ′) ≤ 4g₂ + o(1): a lower bound λ′ ≥ λe^{−4g₂ − o(1)}, a certificate for Q′ sharp to e^{o(1)} relative to the
  undeleted bottom; by (2) equivalently 1 − θ ≥ w₀e^{−4g₂ − o(1)} (a lower bound for the Christoffel-type share, the multiplier side)
  with an upper bound for w₀ (ground-state side).
So each half needs one sampling-side and one ground-state-side estimate, sharp to e^{o(1)}; Mori's machinery gives the sampling
side only to e^{O(log μ)} (§3: d = O(1)–O(log μ) is exactly the unknown part).

## 5. Certificates (certmulti.py) — same inequalities as MoriRate certlow.py, organised for many model forms
Model form Q_S = Q − Σ_{k∈D} 2⟨·, e_{γ_k}⟩² + Σ_{a∈A} 2⟨·, e_a⟩²: deleted zeros = extra B columns √2 e_{γ_k}, added points =
extra U columns e_a with M-weight 2; Q_S − ε ≥ D + UMU* − B′B′* ≥ 0 iff λ_max(B′*P⁻¹B′) ≤ 1 (Chekmachev Thm 5.10, Woodbury),
positivity of Id − W by Lemma 5.5 (approximate Cholesky, Arb residual), as in ChekAudit cert.py / MoriRate certlow.py.
Changes (no new inequality): (a) d_j(ε) = d⁰_j − ε with d⁰_j = ((μ(2j) + 1)/2).lower(), so B weights √(½ − d⁰_j) are ε-free;
(b) Gram(ε) = G₀ + εG₁ + E(ε), G₀ = Φᵀdiag(1/d⁰ − 2)Φ + 2·IP, G₁ = Φᵀdiag(1/d⁰² − 4)Φ + 4·IP (stored at 128 bits as balls),
|E_ab| ≤ ε²(κ n_a n_b + 16 ‖a‖‖b‖), κ = 1/(d_min²(d_min − ε)) + 16, n_a² = Σ_{j∈G} Φ_ja² (exact second-order remainders of
1/(d − ε) and 1/(½ − ε), Cauchy–Schwarz), added as a ball; (c) ε chosen by linear extrapolation of the margin, then the rigorous
check at the returned ε (0.985 of the extrapolated root); (d) inner products of two evaluation columns
⟨cos γx, cos γ′x⟩ = sin((γ−γ′)h)/(γ−γ′) + sin((γ+γ′)h)/(γ+γ′), h = a. Sanity: margins at ε = 0 equal ChekAudit's Δ to all printed
digits (c = 5: 4.002816064e−18; 10: 5.026355931e−44; 16: 7.293623148e−76) and MoriRate certlow's deleted margins (7.4065e−22 at c = 5).
(e) Optional a-posteriori Woodbury (SOLVEMODE=resid, memory): Z̃ ≈ S⁻¹Y by approximate LU, R = Y − SZ̃ in balls; the exact S* ⪰ diag(1/M)
(Gram_UU is a Gram matrix of D⁻¹ > 0), so ‖S*⁻¹‖ ≤ max M and |(Y*ᵀS*⁻¹R*)_ab| ≤ ‖Y_a‖ max M ‖R_b‖, added as a ball to X − YᵀZ̃; norms
from upper endpoints (sqrt of a ball straddling 0 is NaN in Arb). Same certified ε as the standard mode at c = 5, 10 (all digits).
(f) More nodes: K2 = 3 × max(201, ⌈5c²/3⌉) at c = 10 raises ε/Galerkin from 0.61–0.65 to 0.79–0.83 for the full and every model form.
At c = 20 the certificate has order ≈ 1100 (Kp = 788, |G| = 8032), beyond the 300 MB cap for dense Arb matrices (killed at 291 MB in
assembly); ExploreFMM's structured O(n) LDL would be needed.

## 6. Upper bounds beyond the solvable N (tailritz.py) and the finite-μ statements
Rayleigh–Ritz on the orthonormal set {x/|x|} ∪ {φ_j : N₁ < j ≤ N₂}, x = Σ_{k≤N₁} x_kφ_k a Galerkin ground state (exact midpoint
coefficients): H₀₀ = xᵀM x/xᵀx, H₀p = Σ_k M_{j_p k}x_k/|x|, H_pq = M_{j_p j_q} (closed-form entries in Arb, |x| = √(xᵀx) as a ball);
any exact w gives λ ≤ wᵀHw/wᵀw (reviewer: CORRECT). Gains over N₁: 0.018 (μ = 5, 200 → 700), 0.030 (μ = 10, 520 → 1000),
0.026 (μ = 16, 520 → 900) nats; RSS ≤ 222 MB (no N₂ × N₂ solve).
Finite-μ consequences at μ = 10 (F(S) = log(λ/λ_min(Q_S)), certified endpoints, margins from their logarithms):
- F(γ) < 4g₂(y) for γ₄, γ₉, γ₁₀, γ₁₄ (y = 0.484, 0.764, 0.792, 0.968) by ≥ 0.0035, 0.0218, 0.0229, 0.0280: the exact equality fails.
- F({γ₁, γ₁₄, γ₃₅}) ≥ 19.0779 > 19.0203 ≥ F(γ₁) + F(γ₁₄) + F(γ₃₅): strict super-additivity (margin 0.058). For the other multi-sets only
  F(S) > Σ4g₂ is certified (singles of γ₃, γ₂₁, γ₂₇ not certified).
- A resolved pair (γ₁₀ moved to γ₉ + 10⁻¹⁰) has F ∈ [−0.341, −0.008]: it raises λ.
None of these bears on the asymptotic law (an o(1) statement).
