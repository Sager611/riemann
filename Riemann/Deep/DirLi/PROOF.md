# PROOF: explicit Li-coefficient positivity for Dirichlet L-functions from GRH to a finite height

Node N3 (Li / Keiper–Li). Certificate: certify.py (this directory). Sanity checks (not proof): sanity.py.

## 0. Setting

Let χ be a primitive Dirichlet character mod q ≥ 3, a_χ = (1 − χ(−1))/2, and
ξ(s, χ) = (q/π)^{(s+a_χ)/2} Γ((s+a_χ)/2) L(s, χ), entire of order 1. Its zeros are the nontrivial zeros
ρ = β + iγ of L(s, χ); all satisfy 0 < β < 1 (L(1+it, χ) ≠ 0 and the functional equation
ξ(s, χ) = ε(χ) ξ(1 − s, χ̄), |ε(χ)| = 1). Since ξ(s̄, χ̄) = conj ξ(s, χ), the zero multiset is invariant under
ρ ↦ ρ* := 1 − ρ̄ (same ordinate γ, real part 1 − β, same multiplicity). It is not invariant under ρ ↦ ρ̄ unless χ is real.

Li coefficients (Li 1997; Li, Illinois J. Math. 48 (2004); Lagarias, Ann. Inst. Fourier 57 (2007); Omar–Ouni–Mazhouda,
arXiv 1507.03431): with z_ρ = 1 − 1/ρ,
  λ_n(χ) = Σ_ρ [1 − z_ρ^n] := lim_{T→∞} Σ_{|γ|≤T} [1 − z_ρ^n].
Generalized Li criterion: GRH for L(s, χ) ⟺ Re λ_n(χ) ≥ 0 for all n ≥ 1 (then even > 0). For real χ, λ_n(χ) ∈ ℝ;
for complex χ, λ_n(χ̄) = conj λ_n(χ), so Re λ_n(χ) is the natural quantity (the brief's "λ_n(χ) > 0").

Remark (absolute convergence of the real part). Re(1 − z_ρ^n) = n β/|ρ|² + O(n²/|ρ|²) and Σ_ρ 1/|ρ|² < ∞, so
Re λ_n(χ) = Σ_ρ Re(1 − z_ρ^n) converges absolutely; no ordering issue arises below.
Remark (Taylor normalisation, used only in sanity.py). log ξ(1/(1−z), χ) = const + i β_χ/(1−z) + Σ_n c_n z^n/n with
c_n = conj λ_n(χ) (zeros ρ ↦ 1 − ρ̄) plus the term n·iβ_χ from the imaginary linear factor e^{iβ_χ s} of the Hadamard
product; hence Re c_n = Re λ_n(χ). sanity.py checks this against a direct zero sum (χ mod 3, 189 zeros up to 300).

Notation: θ(t) = 2 arctan(1/(2t)) (t ≥ 0, θ(0) = π), g(x) = 2 − 2 cos x, G_n(t) = g(nθ(t)),
t_n = ½ cot(π/(2n)) (so nθ(t_n) = π; G_n decreases from 4 to 0 on [t_n, ∞)), N(t) = N(t, χ) = #{ρ : |γ| ≤ t}
(multiplicity), q_n(t) = 4 sinh²(n a/(2t²)) + 4 n a²/t³.

## 1. Inputs

(I1) Platt, "Numerical computations concerning the GRH", Math. Comp. 85 (2016) 3009–3027, arXiv 1305.3087, final
theorem: for every primitive χ mod q ≤ 400,000, every zero with |γ| ≤ t0(q) has β = 1/2, where
t0(q) = max(10^8/q, 7.5·10^7/q + 200) for even q and max(10^8/q, 3.75·10^7/q + 200) for odd q
(N_χ(t0) in his Turing-method lemma counts zeros with |Im s| ≤ t0, so both signs of γ are covered; his count of
29,565,923,837 L-functions equals Σ_{q≤400000} φ*(q) ≈ (36/π⁴)·(4·10^5)²/2, i.e. every primitive χ, χ and χ̄ alike).
(I2) Bennett–Martin–O'Bryant–Rechnitzer, "Counting zeros of Dirichlet L-functions", Math. Comp. 90 (2021)
1455–1482, arXiv 2005.02989 (not 1812.11592, which is an unrelated paper), Theorem 1.1: N(T, χ) = #{0<β<1, |γ| ≤ T};
for conductor q > 1, T ≥ 5/7, ℓ(T) = log(q(T+2)/(2π)): if ℓ ≤ 1.567 then N(T, χ) = 0; if ℓ > 1.567 then
|N(T, χ) − Ñ(T)| ≤ E(T), Ñ(T) = (T/π) log(qT/(2πe)) − χ(−1)/4, E(T) = 0.22737 ℓ + 2 log(1 + ℓ) − 0.5.
(I3a) Classical: 0 < β < 1, i.e. |β − ½| < a with a = ½.
(I3b) Zero-free half-plane Re s > 7/8 for every Dirichlet L-function (OpenAI, Sep 2026, unrefereed; Lean port
Riemann/Deep/QRH, theorem OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re, used in
Riemann/Deep/QRHCor/Dirichlet.lean). With ρ ↦ 1 − ρ̄: 1/8 ≤ β ≤ 7/8, i.e. |β − ½| ≤ a with a = 3/8.

Facts about E: E is increasing in T and in q; E'(T) = (0.22737 + 2/(1 + ℓ(T)))/(T + 2) ≤ c_E(t1)/T for T ≥ t1, where
c_E(t) = 0.22737 + 2/(1 + ℓ(t)) (decreasing in t and in q). Ñ'(T) = (1/π) log(qT/(2π)), increasing in q.

## 2. Lemma 1 (per zero)

(a) If β = ½: Re(1 − z_ρ^n) = 1 − cos(nθ(|γ|)) = ½ G_n(|γ|).
(b) If ρ = ½ + δ + iγ with 0 < δ ≤ a ≤ ½ (δ < ½), γ ≠ 0, and ρ* = 1 − ρ̄ = ½ − δ + iγ, then
  Re(1 − z_ρ^n) + Re(1 − z_{ρ*}^n) ≥ G_n(|γ|) − q_n(|γ|).
Proof. (a) For ρ = ½ + iγ, γ > 0: z = (−½ + iγ)/(½ + iγ) has |z| = 1 and arg z = arctan(1/(2γ)) + arctan(1/(2γ)) = θ(γ);
γ < 0 gives the conjugate; γ = 0 gives z = −1 = e^{iπ}, θ(0) = π.
(b) Take γ > 0 (γ < 0: z_ρ, z_{ρ*} are conjugated, real parts unchanged). z_ρ = r e^{iφ}, z_{ρ*} = 1/z̄_ρ = r^{−1} e^{iφ}
with r² = ((½−δ)² + γ²)/((½+δ)² + γ²) and φ = f(½ + δ) + f(½ − δ), f(x) = arctan(x/γ). The sum equals
2 − 2 cosh(nℓ) cos(nφ), ℓ = −log r, = g(nφ) − 2(cosh nℓ − 1) cos nφ ≥ g(nφ) − 4 sinh²(nℓ/2).
ℓ = ½ log(1 + 2δ/((½−δ)² + γ²)) ≤ δ/γ² ≤ a/γ². φ − θ(γ) = f(½+δ) + f(½−δ) − 2f(½) = δ² f''(ξ), ξ ∈ (½−δ, ½+δ) ⊂ (0, 1),
|f''(ξ)| = 2ξγ/(ξ² + γ²)² ≤ 2/γ³, so |nφ − nθ| ≤ 2nδ²/γ³ and, as |g'| ≤ 2, g(nφ) ≥ G_n(γ) − 4na²/γ³. ∎
(certify.py part0 tests (b) on 4000 random cases, both signs of γ, a ∈ {1/8,…,1/2}.)

## 3. Corollary (the basic inequality)

Let t0 = t0(q) and a as in (I3a) or (I3b). Every zero with |γ| ≤ t0 is on the line (I1); off-line zeros come in pairs
{ρ, ρ*} (equal multiplicity) with |γ| > t0. Summing Lemma 1 (absolutely convergent sums):
  Re λ_n(χ) ≥ ½ M_n − ½ Q_n,   M_n = Σ_ρ G_n(|γ|) = ∫_{[0,∞)} G_n dN,   Q_n = Σ_{|γ|>t0} q_n(|γ|) = ∫_{(t0,∞)} q_n dN.
M_n uses ordinates only (no RH): every term is ≥ 0, so zeros may be dropped from M_n.

## 4. Lemma 2 (Stieltjes with an error band)

Let N = Ñ + R on [α, ∞) with |R| ≤ E, E nondecreasing, Ñ ∈ C¹, N nondecreasing right-continuous.
(a) H ≥ 0 continuous monotone on [α, β]: ∫_{(α,β]} H dN ≥ ∫_α^β H Ñ' dt − 2 max(H(α), H(β)) E(β).
(b) H decreasing to 0, H·N → 0 and H·E → 0 at ∞, E ∈ C¹: ∫_{(α,∞)} H dN ≥ ∫_α^∞ H Ñ' dt − 2H(α)E(α) − ∫_α^∞ E'H dt,
    and ∫_{(α,∞)} H dN ≤ ∫_α^∞ H (Ñ' + E') dt + 2H(α)E(α).
Proof. (a) ∫ H dN = [HN]_α^β − ∫ N dH = ∫ H dÑ + H(β)R(β) − H(α)R(α) − ∫ R dH, |∫ R dH| ≤ E(β)|H(β) − H(α)|; both monotone
cases give ≥ −2 max(H(α),H(β)) E(β). (b) ∫_{(α,∞)} H dN = −H(α)N(α) + ∫_α^∞ N |H'| dt; insert N = Ñ + R, |R| ≤ E, and
∫_α^∞ E|H'| = E(α)H(α) + ∫_α^∞ E'H. ∎ (For a general monotone H read |H'|dt as the measure −dH; every H used here is
smooth.) The constant −χ(−1)/4 of Ñ cancels in every identity (only Ñ' enters), so the
bounds hold for both parities; only the integer roundings of §5C use ±¼ explicitly.

## 5. Gain: lower bounds for ½M_m, for all integers m in a block [n1, n2]

Common facts: G_m(t) ≤ (mθ(t))² ≤ m²/t² (g(x) ≤ x², θ(t) ≤ 1/t). Substitution u = mθ(t), t(u) = ½cot(u/(2m)),
dt = −(1/(4m)) csc²(u/(2m)) du ≥ (m/u²)|du| in absolute value, t(u) = (m/u)·x cot x with x = u/(2m), and
x cot x ≥ 1 − x²/2 ≥ e^{−x²} for 0 < x ≤ π/4 (series x cot x = 1 − x²/3 − Σ_{k≥2} 2ζ(2k)(x/π)^{2k}, tail ≤ 0.025x⁴ ≤ x²/6).
Hence where log(q t/(2π)) ≥ 0: Ñ'(t)|dt| ≥ (m/π u²)(log(qm/(2π)) − log u − x²) du.
K0 = ∫_0^π (1−cos u)/u² du = Si(π) − 2/π = 1.215317…, K1 = ∫_0^π (1−cos u) log u/u² du = −0.0516076…,
∫_0^∞ (1−cos u)/u² du = π/2, ∫_U^∞ (1−cos u)/u² du ≤ 2/U, K1∞ = ∫_0^∞ (1−cos u) log u /u² du = (π/2)(1 − γ_E)
(Mellin transform; the integrand is ≥ 0 for u ≥ 1, so ∫_0^U ≤ K1∞ for U ≥ 1).

A (zeros with |γ| > t_m; m ≥ 3, t_{n1} > 5/7, ℓ(t_{n1}) > 1.567, q t_{n1} > 2π): by Lemma 2(b) with H = G_m (decreasing
from G_m(t_m) = 4), and ∫ E'G_m ≤ c_E(t_m) m²/(2t_m²),
  ½M_m ≥ (m/π)[K0(log(qm/(2π)) − π²/(4m²)) − K1] − 4E(t_m) − c_E(t_m) m²/(4t_m²).
Block: the main term increases in m once its bracket is > 0 at n1 (certified); E(t_m) ≤ E(t_{n2}); c_E(t_m) ≤ c_E(t_{n1});
m/t_m = 2m tan(π/(2m)) decreases in m.
B (zeros with |γ| > t_low = √n1 ≤ t_{n1}; used for n1 ≥ 10^4): on [t_low, ∞), G_m is monotone on the pieces cut at the
points τ_k with mθ(τ_k) = kπ (k = 1, …, ⌊U_m/π⌋, U_m = mθ(t_low)); each bounded piece ends at ≤ t_m and costs ≤ 8E(t_m)
(Lemma 2(a), G ≤ 4), the last piece [t_m, ∞) costs 8E(t_m) + c_E(t_m)m²/(2t_m²) (Lemma 2(b)). With x ≤ θ(t_low)/2 ≤ 1/(2t_low):
  ½M_m ≥ (m/π)[(π/2 − 2/U_m)(log(qm/(2π)) − 1/(4t_low²)) − (π/2)(1 − γ_E)] − 4E(t_m)(U_m/π + 2) − c_E(t_m) m²/(4t_m²).
Block: certified U_{n1} > 4/π and bracket(n1) > 0 make the main term increasing in m; costs taken at n2 (E, U) and n1 (c_E).
For n1 comparable to m this is asymptotically (m/2)(log m + log(q/2π) − 1 + γ_E), the GRH main term of λ_m(χ) (Lagarias).
C (single n, small n): s0 = max(t_n, 5/7) (t_n < 5/7 only for n ≤ 2; then G_n is decreasing on [s0, ∞) as well).
By Fubini, Σ_{|γ|>s0} G_n(|γ|) = ∫_{s0}^∞ (N(s) − N(s0)) |G_n'(s)| ds. On a geometric grid s0 < s1 < … < sK:
  ½M_n ≥ ½ Σ_k max(0, N⁻(s_k) − N⁺(s0)) (G_n(s_k) − G_n(s_{k+1})) + ½ max(0, N⁻(s_K) − N⁺(s0)) G_n(s_K),
with integer bounds N⁻(s) = max(0, ⌈Ñ₀(s) − ¼ − E(s)⌉) if ℓ(s) > 1.567 (else 0) and N⁺(s0) = 0 if ℓ(s0) ≤ 1.567, else
⌊Ñ₀(s0) + ¼ + E(s0)⌋ (Ñ₀ = Ñ without the parity constant; ±¼ covers both parities).

## 6. Loss: upper bound for ½Q_n

q_n is decreasing in t and increasing in n. Lemma 2(b) (upper form) with α = t0:
  ½Q_n ≤ ½[∫_{t0}^∞ q_n(t) w(t) dt + 2q_n(t0)E(t0)],   w(t) = (1/π) log(qt/(2π)) + c_E(t0)/t ≥ Ñ' + E' (> 0).
The integral is bounded by an upper Riemann sum Σ q_n(p_k)(W(p_{k+1}) − W(p_k)) (W' = w) on a partition of [t0, p_K] that is
fine (steps 0.004 in x = na/(2t²)) where x ∈ [X − 6, X], X = na/(2t0²), plus the tail ∫_{p_K}^∞ ≤ C n²a² ∫ w/t⁴ + 4na² ∫ w/t³,
C = (sinh x_K/x_K)² (sinh x/x increasing), with ∫_{T}^∞ log(ct)/t⁴ = (3L+1)/(9T³), ∫_T^∞ log(ct)/t³ = (2L+1)/(4T²), L = log(cT).
For X < 10^{−3} only the tail formula from t0 is used.

## 7. Theorem and certification

Theorem. For a primitive χ mod q (3 ≤ q ≤ 400,000) and integers m in [n1, n2]: if gain(n1, n2) − loss(n2) > 0 then
Re λ_m(χ) > 0 (gain = max(A, B) of §5, loss of §6 at n2 ≥ m since Q_m ≤ Q_{n2}). Single n via C.
certify.py covers [1, N] by regime C for n ≤ n_C (until A certifies a single n) followed by a chain of blocks
[n_i + 1, n_{i+1}], each a certified ball-arithmetic sign; floats only choose the block ends. The chain stops when the float
version of the same inequality fails at n_i + 1 or the relative step falls below 10^{−9}.

q-blocks (moduli q1 ≤ q ≤ q2 of one parity class): every q-dependence above is monotone. The gain uses q1 in log(qm/2π) and
in c_E, q2 in E; the loss uses q2 in w and E(t0), q1 in c_E, and t0 = t0(q2) ≤ t0(q) (t0 decreasing in q within a parity
class; zeros of L(s, χ) with |γ| > t0(q) are among those with |γ| > t0(q2), and q_n ≥ 0). N⁻ uses q1 in Ñ₀ and q2 in E;
N⁺ uses q2. So one chain certifies all moduli of the block.

## 8. What is not used / caveats

No zeros are computed in the proof. Neither Trudgian's 2015 Dirichlet zero-count bound (reported error, see
Hasanalizade–Shen–Wong 2022 fn. 2, following BMOR) nor any zeta input enters. Variant a = 3/8 rests on OpenAI's unrefereed
7/8 theorem (Lean-checked port in this repo); variant a = 1/2 uses refereed inputs only (I1, I2, I3a).
