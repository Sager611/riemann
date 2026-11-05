# explore6 derivations (Möbius disjointness → RH tree)

Notation: x = 2^λ; s_2(n) = binary digit sum; F_g(α) = Σ_{n<x} g(n) e(α s_2(n)) = Σ_{k=0}^{λ} c_k z^k with z = e(α);
M(x) = Σ_{n≤x} μ(n); Θ = exponent of M (|M(x)| = x^{Θ+o(1)}); B_2 = (777 − √1081)/852 = 0.8733819666.

## D1. Bernstein–Walsh / two-constants inequality on the arc (used in T1)
Arc A_δ = {e^{iθ}: 2πδ ≤ θ ≤ 2π − 2πδ} (all α with ‖α‖ ≥ δ), domain Ω = Ĉ ∖ A_δ, which contains z = 1 and ∞.
For a polynomial P of degree ≤ λ, u = log|P| − λ g_Ω(·, ∞) is subharmonic on Ω (bounded at ∞) and equals log|P| on A_δ, so

  (D1.1)  log|P(1)| ≤ λ g_Ω(1, ∞) + ∫_{A_δ} log|P| dω_1        (two-constants form)
  (D1.2)  log|P(1)| ≤ λ g_Ω(1, ∞) + log sup_{A_δ}|P|          (Bernstein–Walsh sup form),

with Riesz mass at the zeros: log|P(1)| = λ g_Ω(1,∞) + ∫ log|P| dω_1 − Σ_{P(ζ)=0, ζ∈Ω} g_Ω(1, ζ). Equality in (D1.1) iff all
zeros of P lie on A_δ.

**Conformal map.** w = i(1−z)/(1+z) sends e^{iθ} ↦ tan(θ/2), so A_δ ↦ {|w| ≥ a}, a = tan(πδ), 1 ↦ 0, ∞ ↦ −i. Then ζ = w/a,
v = (1+ζ)/(1−ζ) (slit ↦ (−∞, 0]), u = √v (Ω ↦ right half-plane), 1 ↦ u = 1, ∞ ↦ p = e^{−iβ}, β = (π − 2πδ)/2.
- Green function: g_Ω(1, ∞) = log|(1 + p̄)/(1 − p)| = log cot(β/2) = **log((1 + sin πδ)/cos πδ)**.
- Harmonic measure at 1: dω_1 = dt/π on u = i tan t, t ∈ (−π/2, π/2); the boundary point is θ(t) = 2 arctan(−a sec 2t)
  (t and −t are the two sides of the arc). Used for quadrature in t1_analyze.py.
- Check: P = (1+z)^λ has its only zero at −1 ∈ A_δ, so (D1.1) is an equality; t1_analyze.py confirms 1.0000 for
  δ = 0.002 … 0.3, and the sup form exceeds it by exactly λ log_2(1 + sin πδ) (base-2 units).

## D2. The g ≡ 1 polynomial and the closed loop (T1)
- For g ≡ 1: F_1(α) = (1 + e(α))^λ, |F_1(α)| = x |cos πα|^λ = x^{1 − c(α)}, c(α) = −log_2|cos πα| ≈ 7.12 α² (α small).
- For μ write c_k = M·b_k + r_k with b_k = binom(λ,k)/2^λ, so F_μ = M·P_0 + R, P_0 = ((1+z)/2)^λ, R(1) = 0 (exact).
- Exact identity: μ(2m) = −μ(m), s_2(2m) = s_2(m) for odd m, so c_k(λ) = Σ_{m odd, 2^{λ−1} ≤ m < 2^λ, s_2(m)=k} μ(m).
  With the matching random model (independent signs on odd squarefree m, density 8/π², class size binom(λ−2, k−2)) the
  residuals r_k have rms z ≈ 0.64–1.06 for λ = 15 … 27 (t1_out.txt): the classes look random; no hidden rigidity.
- In exponent units, (D1.1) with a proven bound |F_μ(α)| ≤ x^{1−σ(α)} on A_δ gives
  **Θ_out ≤ 1 − ∫_{A_δ}(σ(α) − c(α)) dω_1(α)**   (because g_Ω(1,∞)/ln 2 = ∫ c dω_1 exactly, by D1 equality for (1+z)^λ).
- The frequency-0 component M·P_0 alone has |M·P_0(α)| = x^{Θ − c(α)}, so any valid bound has σ(α) ≤ 1 − Θ + c(α) wherever
  it dominates R, and then Θ_out ≥ Θ: the transfer returns the exponent fed in, never less. With the sup form (D1.2) the
  output is Θ + log_2(1 + sin πδ) (strictly worse).
- Type I/II (Gelfond, Mauduit–Rivat, Bourgain) proofs bound the m = 1 Type I term Σ_{k≤x} f(k) — the g ≡ 1 sum — in
  absolute value, so their saving satisfies σ(α) ≤ c(α) + o(1): plugged into D2 they give Θ_out ≥ 1.
- Breaking the loop needs |F_μ(α)| below |M|·x^{−c(α)}, i.e. knowledge of M (or cancellation of μ in the m-sum of Type I
  terms at frequency 0, which is a Mertens-type estimate at s = 1 − c(α) + 2πij/ln 2). With the B_2 half-plane as input
  that cancellation is available and the loop returns B_2.
- RH-sufficient form (true but stronger than RH, open even under GRH): if for each fixed δ > 0,
  sup_{‖α‖≥δ}|F_μ(α)| ≪ x^{1/2+ε}, then Θ ≤ 1/2 + log_2((1+sin πδ)/cos πδ) for every δ, hence RH.
- Numbers to beat B_2 via the sup form: σ(δ) > (1 − B_2) + c(δ) + log_2(1 + sin πδ) = 0.172, 0.354, 0.588 at δ = 0.01,
  0.05, 0.10, against the method cap c(δ) = 0.0007, 0.018, 0.072.

## D3. Correlations and Mertens (T2)
- (D3.1) M(X)² = C(0) + 2 Σ_{1≤h<X} C_X(h), C_X(h) = Σ_{n≤X−h} μ(n)μ(n+h), C(0) = Σ_{n≤X} μ(n)² ≈ (6/π²)X.
- (D3.2) Σ_{x≤X} M(x)² = D(0) + 2 Σ_{h≥1} D(h), D(h) = Σ_n μ(n)μ(n+h)(X+1−n−h).
- (D3.3) Short windows: Σ_x (M(x+H) − M(x))² = Σ_{|h|<H}(H − |h|) C(h) (Good–Churchhouse predicts ≈ (6/π²)H per x).
- Triangle route: M(X)² ≤ C(0) + 2Σ_h |C_X(h)|. If |C_X(h)| ≪ X^{1−η+ε} uniformly for h < X then Θ ≤ 1 − η/2;
  square-root Chowla (η = 1/2) gives Θ ≤ 3/4; beating B_2 needs η > 2(1 − B_2) = 0.2532 uniformly up to h ≍ X.
  Measured at X = 2^22: √(C(0) + 2Σ|C|) = X^{0.7355} for μ and for random signs alike, while |M(X)| = X^{0.356}.
- Fixed-shift theorems (MRT averaged, Tao log, Pilatte, OpenAI ordinary two-point, all log-power) control C_X(h) for h
  fixed or h ≤ H = o(X); the identity needs signed cancellation over h up to X.

## D4. Rotations under a zero-free half-plane for all Dirichlet L (T3)
Assume L(s,χ) ≠ 0 for Re s > Θ for every Dirichlet χ (true here with Θ = B_2 + 2.2·10^{−254}). Then 1/L(s,χ) ≪ (q(|t|+2))^ε on
Re s ≥ Θ + ε (Borel–Carathéodory + convexity), and Σ_{m≤y} μ(m)χ(m) ≪ y^{Θ+ε}q^ε (Perron).
- Dirichlet with Q = x/R: α = a/q + β, q ≤ x/R, |β| ≤ R/(qx).
- q ≤ R (major): n = dm, d = (n,q); expand e(am/q') over χ mod q' (|τ(χ)| ≤ √q'); partial summation in e(βdm):
  S(α) ≪ q^{1/2+ε}(1 + |β|x) x^{Θ+ε} ≪ R x^{Θ+ε}.
- q > R (minor; Vaughan-type bound for μ): S(α) ≪ (x q^{−1/2} + x^{4/5} + (qx)^{1/2}) log^4 x ≪ x R^{−1/2} + x^{4/5}.
- R = x^{2(1−Θ)/3} gives **sup_α|Σ_{n≤x} μ(n)e(nα)| ≪ x^{(2+Θ)/3+ε}** = x^{0.95779} at Θ = B_2 (GRH: 5/6; Baker–Harman's
  sharper GRH treatment gives 3/4). The family contains α = 0, so its exponent is ≥ Θ: no route back to N7.
- Universality cap (Kátai–BSZ, Daboussi, Montgomery–Vaughan use only multiplicativity and |g| ≤ 1): g(p) = e(−pφ),
  g(p^k) = 0, gives Σ_{n≤x} g(n)e(nφ) = π(x)(1 + o(1)) (measured 1.0013 at 10^7), so such arguments cap at x/log x even
  for the golden rotation, where μ itself measures ≈ 0.4–0.9 √x.
