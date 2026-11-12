# ExploreDPP: derivations (explore34, 12 Nov 2026)

Notation: f^(u) = ∫ f(α) e^{−2πiαu} dα; Δ(α) = (1 − |α|)_+; sinc x = sin(πx)/(πx); J_θ(f) = f(0) + 2∫_0^θ f(α) α dα.
Lamzouri 2611.14486 Thm Main2 (l. 392-400; SmallGap/DERIVATION.md §1): if 1/2 < θ < 4/7, f ∈ C_c²(R) even, f ≥ 0, f^(0) = 1,
f^(u) ≤ 0 for |u| > λ and J_θ(f) > 1, then D(λ) > 0. SmallGap's Lemma T: same with Schwartz f. Call such f admissible for (θ, λ).
Both reduce to the LP: is sup J_θ(f) − f^(0) over admissible f positive? θ = 1 is the RH (Montgomery) case of the same functional
(SmallGap's SDP at θ = 1 gives 0.60400 = CGdL's 0.6039).

## 1. Weak duality (lattice dual)
Data: h > 0, k0 ≥ 1 with k0 h = λ, K ≥ k0, reals η_1..η_K with η_k = −1 (k < k0), η_k ≥ −1 (k ≥ k0); η_k := 0 for k > K.
σ := Σ_{k≠0} h(1 + η_|k|) δ_{kh} ≥ 0, supported on |u| ≥ λ. P(α) := (1 − h) + 2h Σ_{k=1}^K η_k cos(2παkh) (even, period 1/h).
Lemma 1. If P(α) ≥ |α|·1[|α| ≤ θ] on [0, 1/(2h)] and θ ≤ 1/(2h), then J_θ(f) ≤ f^(0) for every admissible f.
Proof. f ∈ C_c² gives f^(u) = O(u^{−2}), so Poisson summation h Σ_k f^(kh) = Σ_m f(m/h) holds, and
 ∫ f^ dσ = Σ_m f(m/h) − h∫f + ∫ f(α)·2hΣ_k η_k cos(2παkh) dα = f(0) + Σ_{m≠0} f(m/h) + ∫ f (P − 1).
Subtract J_θ(f) − f^(0) = f(0) + ∫ f(α)(|α|1_θ − 1) dα:
 ∫ f^ dσ − (J_θ(f) − f^(0)) = Σ_{m≠0} f(m/h) + ∫ f(α)(P(α) − |α|1_θ(α)) dα ≥ 0,
since f ≥ 0 and P − |α|1_θ ≥ 0 pointwise on R: on [−1/(2h), 1/(2h)] by hypothesis and evenness; for |α| > 1/(2h) ≥ θ the
target is 0 and P(α) = P(α − m/h) with α − m/h in the base period, where P ≥ |·|1_θ ≥ 0. θ = 1/(2h) is allowed: [−θ, θ] is then the
closed base period, and its endpoints, shared with the neighbouring periods, satisfy P ≥ 0 there too (only θ ≤ 1/(2h) is used, to
keep the support of 1_θ inside one period up to its endpoints; the integral ignores null sets anyway).
σ ≥ 0 sits on |u| ≥ λ, where f^ ≤ 0 (at |u| = λ by continuity), so ∫ f^ dσ ≤ 0. Hence J_θ(f) − f^(0) ≤ 0. □
Monotonicity: admissible for (θ', λ') with θ' ≤ θ, λ' ≤ λ implies admissible for (θ, λ) and J_θ' ≤ J_θ; so one certificate kills all
smaller θ and λ. Schwartz f (Lemma T): same proof.
Reading: σ is a pseudo pair-correlation measure with hard core λ; G = Σ_{m≠0} δ_{m/h} + P − |α|1_θ ≥ 0 says its structure factor P
is ≥ Montgomery's |α| on [−θ, θ] and ≥ 0 elsewhere. Any stationary simple point process of intensity 1 with S(α) = |α| on [−θ, θ]
and no two points closer than λ gives such a certificate (with continuous instead of lattice support).

## 2. Theorem A: the determinantal (discrete sine) dual, floor ≥ 1/(1+θ)
a-discrete sine process (Lagarias–Rodgers 1905.12123, Prop. discrete_sine_exists; Macchi–Soshnikov for 0 < a ≤ 1): the DPP on aZ
with K(x, y) = a·sinc(x − y). HKPV applies: counts are sums of independent Bernoulli(eigenvalues of K restricted).
Pair measure (Palm, intensity 1): σ_a = Σ_{k≠0} a(1 − sinc²(ka)) δ_{ka}, i.e. Lemma 1 with h = a, k0 = 1, η_k = −sinc²(ka), K = ∞.
Poisson (FT of sinc² is Δ): a Σ_{k∈Z} sinc²(ka) e^{2πiαka} = Σ_m Δ(α − m/a), so P(α) = 1 − Σ_m Δ(α − m/a) =: S_a(α).
For |α| ≤ 1/a − 1 only m = 0 contributes: S_a(α) = |α|. Elsewhere, since 1/a ≥ 1, at most two triangles are nonzero at any α:
one triangle contributes at most 1, and two overlapping triangles sum to exactly 2 − 1/a ≤ 1 on their overlap; so S_a ≥ 0
(S_a(m/a) = 0 for m ≠ 0).
So with a = 1/(1+θ): P = |α| on [−θ, θ], P ≥ 0 elsewhere; Lemma 1 applies since θ ≤ (1+θ)/2 = 1/(2a) for θ ≤ 1, with equality
at θ = 1 (a = h = 1/2), the case Lemma 1 now allows (K = ∞ is fine: the series converges absolutely).
Theorem A. For 0 < θ ≤ 1 and λ ≤ 1/(1+θ) no f is admissible with J_θ(f) > f^(0). θ = 1: λ ≤ 1/2 (the Alternative Hypothesis
obstruction, Lagarias–Rodgers); θ < 4/7: λ ≤ 7/11 = 0.63636.
Probabilistic content: the randomly shifted 1/(1+θ)-discrete sine is a genuine stationary simple point process (a mixture of DPPs)
with Montgomery's structure factor on |α| < θ and every gap ≥ 1/(1+θ); its Palm gap law (logs/hkpv_gaps.log): P(gap = a) = 0.5046,
P(gap = 2a) = 0.4225 for a = 7/11; 0.2974, 0.4526 for a = 1/2.

## 3. Theorem B: certified lattice duals beyond the DPP
LP (dual_lp.py, dual_lp2.py): maximise t subject to P(α_i) − |α_i|1_θ ≥ t on an active set of α_i ∈ [0, 1/(2h)] (cutting
planes from a 0.0002-grid, P evaluated in chunks), −1 ≤ η_k ≤ 20, HiGHS dual simplex in float64; the float solution is only a
candidate. Certificate check, verify_dual.py (python-flint arb/acb, 360 bits, outward rounding throughout; every quantity below is a
ball and only its lower endpoint is used):
- η: rationals with denominator 10^12; η_k = −1 exactly for k < k0 and clamped to ≥ −1 otherwise (asserted), so the weights
  w_k = h(1 + η_k) ≥ 0 and no atom lies in (−λ, λ); h = λ/k0 and θ are exact rationals.
- Grid: all multiples iδ in [0, 1/(2h)] plus θ and 1/(2h), sorted; consecutive points bound closed cells that cover [0, 1/(2h)]
  completely. Since θ is a grid point, no cell straddles θ.
- Cells left of θ (right endpoint ≤ θ) check g = P − α; cells right of θ check g = P. At θ the target |α|1_θ drops from θ to 0, so
  G = P − |α|1_θ jumps up there (a jump, not a kink): the binding lower branch P(θ) − θ is the closed right endpoint of the left cell
  [θ − δ', θ] and is certified by that cell.
- Cell minorant: on [a0, a0 + δ'] (δ' the cell length), g(a0 + x) ≥ g(a0) + g'(a0)x − M2x²/2 ≥ min(g(a0), g(a0) + g'(a0)δ') − M2δ'²/2,
  with M2 = 2hΣ|η_k|(2πkh)² ≥ sup|P''| (the concave quadratic is minimised at an endpoint); g, g' from P, P' computed by Horner in
  acb at z = e^{2πiαh} (360 bits absorb the rectangular-ball growth (|Re z| + |Im z|)^K).
- PASS iff the minimum lower endpoint over all cells is > 0.
- θ = 4/7 (covers every θ < 4/7, hence Thm Main2's whole range), λ = 0.775, h = 31/280, k0 = 7, K = 235 (η on (0, 26]):
  min lower bound of G on [0, 1/(2h)] = 0.000204 > 0 (90,324 cells, δ = 1/20000, M2 = 217,389; 34 s, 78 MB).
- θ = 4/7, λ = 0.77, h = 77/1000, k0 = 10, K = 234: min ≥ 0.00193 (64,937 cells, M2 = 74,094).
- θ = 4/7 and θ = 399993/700000, λ = 0.76, h = 19/375, k0 = 15, K = 197: min ≥ 0.0031679 (49,344 cells, M2 = 17,370.8).
- θ = 1, λ = 0.595, h = 119/600, k0 = 3, K = 202 (η on (0, 40]): min ≥ 0.000239 (100,841 cells, δ = 1/40000, M2 = 480,797).
- θ = 1, λ = 0.59, h = 59/300, k0 = 3, K = 153: min ≥ 0.00361 (50,848 cells); λ = 0.55, h = 1/20, k0 = 11, K = 200: min ≥ 0.0322.
- Not certified: θ = 4/7, λ = 0.78, h = 13/50, K = 269 (float t* = 2.2·10⁻⁴ on a 0.0002 grid): FAIL, min −0.0021 at α ≈ 1.4769.
Theorem B. (i) For every θ ≤ 4/7 and λ ≤ 0.775 no f is admissible with J_θ(f) > f^(0): Lamzouri's Thm Main2 (and SmallGap's Lemma T)
cannot prove D(λ) > 0 for any λ ≤ 0.775. (ii) θ = 1: no λ ≤ 0.595 (the RH pair-correlation floor is above the AH value 1/2).
Upper ends (the method works above them): 0.7851 for θ < 4/7 (SmallGap, certified; its float SDP reaches 0.7849, uncertified) and
0.6039 for θ = 1 (CGdL, published; our float SDP of the J_1 LP: 0.6040). So the floors lie in [0.775, 0.7851] and [0.595, 0.6039].
Floor vs support (float): λ*(X) ≈ 0.7839 − 0.207/X for θ = 4/7 (X = 10…40, residuals ≤ 2·10⁻⁴) and ≈ 0.600 − 0.17/X for θ = 1
(X = 30, 40); appending the kink-exact tail of ψ_h (closed form) does not change λ*(X), so the 1/X loss comes from truncating the
AH-like pseudo-lattice, not from the kink of |α| at 0.
Shape of the λ = 0.775 dual: atoms 0.757 (6.8 × the uniform weight h) at u = λ, 0.37 at 2λ, 0.43 at 3.14λ, 0.44 at 4.14λ, no
mass on [0.89, 1.44] and [2.77, 3.10]: an AH-like pseudo-lattice with spacing ≈ λ (the λ = 0.76 dual: atom 0.737 at λ, empty
[0.81, 1.42]). θ = 1, λ = 0.55: atom 0.421 at 0.55, empty to
1.05, mass at 1.10-1.20.

## 4. What the DPP tool cannot do (the collapse)
A lower bound on D(λ), or on simple-zero or gap proportions, that uses the zeros only through band-limited correlations must also
hold for every point process with the same band-limited correlations. The a-discrete sine DPPs (a ≤ 1/2 for the full Rudnick–Sarnak
range, Lagarias–Rodgers; a ≤ 1/(1+θ) for the pair range θ) are such processes: determinantal, simple, with every HKPV property
(Bernoulli-sum counts, negative association, real-rooted count PGFs) and no gaps below a. So HKPV/strong Rayleigh/Soshnikov
structure adds no constraint that separates the sine process from them. The band-limited bounds of this type (Beurling–Selberg,
CCLM, CGdL, Thm Main2) are LP bounds, and the DPP enters only on their dual side, as a certificate that such a method fails.
Arguments with other inputs (mollifiers, moments of ζ, resonance) are not covered by this.
