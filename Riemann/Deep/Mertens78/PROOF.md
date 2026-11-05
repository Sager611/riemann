# PROOF — explicit power-saving bounds for the Mertens function from the zero-free half-plane Re s > 7/8

Agent mertens78, 2026-11-05. Node N10 (explicit ψ, π − li and Mertens bounds). Companion: REPORT.md (results, novelty, certificates).
Notation: s = σ + it; M(x) = Σ_{n≤x} μ(n); B = 7/8; H = {s : Re s > B}; P_x(u) = x/(π(x² + u²)) (Poisson kernel of a half-plane at height x > 0; ∫P_x = 1).

## 0. Inputs

- (I1) ζ(s) ≠ 0 for Re s > 7/8 (OpenAI 2026; sorry-free Lean port in Riemann/Deep/QRH, theorem Fayed.QRH.zeta_ne_zero_of_seven_eighths_lt_re). Nothing else about zeros is used.
- (I2) Rademacher, Math. Z. 72 (1959) 192–204, Theorem 4 (as quoted in Kadiri–Lumley–Ng, arXiv 2101.12263, Lemma "bndLschi"; checked in their TeX): for 0 < η ≤ 1/2 and −η ≤ σ ≤ 1 + η,
  |ζ(s)| ≤ 3 |1+s|/|1−s| · (|1+s|/2π)^{(1+η−σ)/2} · ζ(1+η).
- (I3) For σ > 1: |ζ(s)| ≤ ζ(σ) and |1/ζ(s)| = |Σ μ(n) n^{−s}| ≤ ζ(σ)/ζ(2σ).
- (I4) Herglotz–Nevanlinna representation of positive harmonic functions on a half-plane (e.g. Axler–Bourdon–Ramey, Harmonic Function Theory, 2nd ed., Ch. 7, positive harmonic functions on the upper half-space; rotated): every harmonic v ≥ 0 on H has the form
  v(B + x + iy) = c·x + ∫_ℝ P_x(y − w) dν(w)   (x > 0),  with c ≥ 0 and ν a positive Borel measure, ∫ dν(w)/(1+w²) < ∞.
- (I5) Rigorous ball arithmetic (Arb via python-flint 0.9.0) for the numerical certificates in §3.

## 1. A lower bound for |ζ| on Re s = σ₁ > 7/8 (upper bound for 1/ζ)

Let Z(s) = (s − 1)ζ(s)/s. By (I1), (s − 1)ζ(s) is entire, equals 1 at s = 1 and has no zero in H; s ≠ 0 on H. So u(s) := log|Z(s)| is harmonic on H.

**Lemma A1 (harmonic majorant).** Fix 0 < η ≤ 1/2, κ := (1 + η − B)/2 = (1/8 + η)/2 and
  h(s) := A + log|1+s| − log|s| + κ log|1+s|,   A := log(3ζ(1+η)) + κ log(8/15).
Then h is harmonic on Re s > 0 and u ≤ h on H.

Proof. Harmonicity: log|1+s|, log|s| are harmonic on Re s > 0. For B < σ ≤ 1 + η, (I2) gives
log|Z(s)| ≤ log(3ζ(1+η)) + log|(1+s)/s| + α(ℓ − log 2π), with α = (1+η−σ)/2 ∈ [0, κ) and ℓ = log|1+s| ≥ log(1+σ) > log(15/8).
The map α ↦ α(ℓ − log 2π) − κ(ℓ − log(15/8)) is affine; at α = 0 it is −κ(ℓ − log(15/8)) ≤ 0, at α = κ it is κ(log(15/8) − log 2π) < 0. So α(ℓ − log 2π) ≤ κ(ℓ − log(15/8)) = κ log(8|1+s|/15), i.e. u ≤ h.
For σ > 1 + η, (I3) and |s − 1| ≤ |s| give u ≤ log ζ(1+η) < h(s), because log 3 > 0, |1+s|/|s| > 1 and 8|1+s|/15 > 1. ∎

So v := h − u is a nonnegative harmonic function on H.

**Lemma A2 (certificates).** Let x₁ > 0 and let λ = Σ_j λ_j δ_{(x_j, y_j)} be a finite signed combination of points with heights x_j > 0 such that
  (C1) Σ_j λ_j P_{x_j}(w − y_j) ≥ P_{x₁}(w) for every real w;   (C2) Σ_j λ_j x_j ≥ x₁.
Then every harmonic v ≥ 0 on H satisfies, for every t ∈ ℝ,  v(B + x₁ + it) ≤ Σ_j λ_j v(B + x_j + i(t + y_j)).

Proof. By (I4), Σ_j λ_j v(B + x_j + i(t+y_j)) − v(B + x₁ + it) = c(Σ_j λ_j x_j − x₁) + ∫ [Σ_j λ_j P_{x_j}(t + y_j − w) − P_{x₁}(t − w)] dν(w). With w = t + w′ and evenness of P, the bracket is the left side minus the right side of (C1) at w′, hence ≥ 0; the first term is ≥ 0 by (C2) and c ≥ 0. (Finite sum of convergent integrals.) ∎

*Harnack case.* λ = ρ δ_{(x₀, 0)} with ρ = x₀/x₁ > 1 satisfies (C1): ρ x₀/(x₀² + w²) ≥ x₁/(x₁² + w²) ⟺ x₀²(x₁² + w²) ≥ x₁²(x₀² + w²) ⟺ (x₀² − x₁²)w² ≥ 0; and (C2): ρx₀ = x₀²/x₁ ≥ x₁. (This is the half-plane Harnack inequality v(B+x₁+it) ≤ (x₀/x₁) v(B+x₀+it); it also follows from the Poisson formula on the discs of radius R centred at B + R + it, R → ∞, without (I4).)

**Theorem A3.** Let σ₁ > B, x₁ = σ₁ − B, and let λ satisfy (C1), (C2) with every point on a line σ_j = B + x_j > 1. Put s_j = σ_j + i(t + y_j), s₁ = σ₁ + it and, for Re s = σ > 1,
  L₀(s) = log(ζ(σ)/ζ(2σ)) + log|s/(s−1)|  (≥ −u(s) by (I3)),   L₁(s) = log ζ(σ) + log|(s−1)/s|  (≥ u(s) by (I3)).
Then for every t ∈ ℝ
  log|1/ζ(σ₁+it)| ≤ Σ_{λ_j>0} λ_j [h(s_j) + L₀(s_j)] − Σ_{λ_j<0} |λ_j| [h(s_j) − L₁(s_j)] − h(s₁) + log|(s₁−1)/s₁|  =: log Φ(t).

Proof. Lemma A2 with v = h − u: −u(s₁) ≤ Σ_j λ_j (h(s_j) − u(s_j)) − h(s₁). For λ_j > 0 use −u(s_j) ≤ L₀(s_j); for λ_j < 0, λ_j(h − u)(s_j) = −|λ_j| h(s_j) + |λ_j| u(s_j) ≤ −|λ_j|(h(s_j) − L₁(s_j)). Finally log|1/ζ(s₁)| = −u(s₁) + log|(s₁−1)/s₁|. ∎

**Corollary A4 (polynomial form).** Let Λ = Σ_j λ_j, c_j = |y_j|, a = (Λ − 1)κ, T₁ ≥ max(1, 2 max_j c_j), and
  log C_∞ = (Λ − 1)A + Σ_{λ_j>0} λ_j log(ζ(σ_j)/ζ(2σ_j)) + Σ_{λ_j<0} |λ_j| log ζ(σ_j) + R/T₁,
  R = σ₁ + Σ_{λ_j>0} λ_j[(1+κ)(c_j + 1 + σ_j) + 2c_j + 2] + Σ_{λ_j<0} |λ_j|[(3 + 2κ)c_j + σ_j].
Then |1/ζ(σ₁ + it)| ≤ Φ(|t|) ≤ C_∞ |t|^a for |t| ≥ T₁. (Integrating (C1) gives Λ ≥ 1, so a ≥ 0.)

Proof (t ≥ T₁; t ↦ −t is symmetric). Elementary: |t + y_j| ≥ t − c_j ≥ t/2; log(t + d) ≤ log t + d/t; log(t − c) ≥ log t − 2c/t for t ≥ 2c.
- λ_j > 0: log|1+s_j| ≤ log(t + c_j + 1 + σ_j) ≤ log t + (c_j + 1 + σ_j)/t; log|s_j| ≥ log|t + y_j| ≥ log t − 2c_j/t; so h(s_j) ≤ A + κ log t + [(1+κ)(c_j+1+σ_j) + 2c_j]/t. Also log|s_j/(s_j−1)| ≤ log(1 + 1/|s_j − 1|) ≤ 1/(t − c_j) ≤ 2/t.
- λ_j < 0: log|1+s_j| ≥ log t − 2c_j/t, log|s_j| ≤ log t + (c_j + σ_j)/t, so −h(s_j) ≤ −A − κ log t + [2(1+κ)c_j + c_j + σ_j]/t; and log|(s_j−1)/s_j| ≤ 0 (σ_j > 1/2).
- s₁: log|1+s₁| ≥ log t, log|s₁| ≤ log t + σ₁/t, log|(s₁−1)/s₁| ≤ 0, so −h(s₁) + log|(s₁−1)/s₁| ≤ −A − κ log t + σ₁/t.
Summing, the coefficient of κ log t is Σλ_j − 1 = Λ − 1, the constants give log C_∞ − R/T₁, the 1/t terms give R/t ≤ R/T₁. ∎

Remark (where the exponent comes from). For the Harnack certificate, a = (x₀/x₁ − 1)κ > ((1 − σ₁)/(σ₁ − B))/16, so θ* = (σ₁ + a)/(1 + a) > 0.92 for every Harnack choice (the infimum 0.92 is approached at σ₁ = 0.9, σ₀ → 1, η → 0: with d = σ₁ − B, 1 − θ* < 16d(1/8 − d)/(15d + 1/8) ≤ 0.08, maximal at d = 1/40). Signed certificates (C1) found by linear programming reduce Λ − 1 at the price of larger Σ|λ_j| (hence larger C_∞); we certify down to θ = 0.8992. In principle every θ > 7/8 is reachable explicitly (the Borel–Carathéodory + three-circles argument of Titchmarsh §14.2 gives log|1/ζ| ≤ c(log t)^α, α < 1, with computable c), but with far worse constants.

## 2. From 1/ζ to M(x)

**Lemma B1.** For y > 0 and c > 1: M₁(y) := ∫₀^y M(u) du = Σ_{n≤y} μ(n)(y − n) = (1/2πi) ∫_{(c)} y^{s+1} ds/(s(s+1)ζ(s)).
Proof. (1/2πi)∫_{(c)} z^s ds/(s(s+1)) = (1 − 1/z)·[z ≥ 1] for z > 0 (standard; residues at 0 and −1 for z > 1, shift right for z < 1). Insert 1/ζ(s) = Σ μ(n)n^{−s} (absolutely convergent on Re s = c, and ∫|ds|/|s(s+1)| < ∞), take z = y/n. ∎

**Lemma B2 (contour shift).** For 0.89 ≤ σ₁ < 1 (all cases used here; for 7/8 < σ₁ < 0.89 the same holds by the qualitative bound 1/ζ(s) ≪_ε |t|^ε on Re s ≥ 7/8 + ε, Lean: Riemann/Deep/PNTError/InvZeta.lean): M₁(y) = (1/2π) ∫_ℝ y^{s+1} dt/(s(s+1)ζ(s)), s = σ₁ + it, whenever the certificate of Corollary A4 has a < 1.
Proof. 1/ζ is analytic on H (its zero at s = 1 is harmless: no residue), and 0, −1 ∉ H. Fix σ₀ = 1.02, η = 0.05 (κ = 0.0875). For σ ∈ [σ₁, σ₀), Theorem A3 with the Harnack certificate at height x = σ − B (ρ = x₀/(σ − B) ≤ x₀/(σ₁ − B)) and Corollary A4 give |1/ζ(σ + iT)| ≤ C |T|^{a′} for |T| ≥ 1, uniformly in σ, with a′ = (0.145/(σ₁ − 7/8) − 1)·0.0875 ≤ 0.76 for σ₁ ≥ 0.89; for σ ≥ σ₀, |1/ζ| ≤ ζ(σ₀)/ζ(2σ₀). So the horizontal integrals over [σ₁, c] ± iT are O(T^{a′−2}) → 0 (a′ < 2 would suffice). The line integral on Re s = σ₁ converges absolutely by Corollary A4 for the certificate used (a < 1, asserted in the scripts). ∎

**Lemma B3 (unsmoothing).** For 0 < h < 2x: |M(x) − (1/h)∫_{x−h/2}^{x+h/2} M(u)du| ≤ h/4 + 1.
Proof. |M(u) − M(x)| ≤ #{n ∈ ℤ between x and u} ≤ |u − x| + 1; average over u. ∎

**Lemma B4.** For x > 0, 0 < h < 2x, X = x + h/2, τ = 2X/h, s = σ₁ + it:
  |M(x)| ≤ h/4 + 1 + (X^{σ₁}/2π) J(τ),   J(τ) := ∫_ℝ |1/ζ(s)| · min(1/|s|, τ/(|s||s+1|)) dt.
Proof. (1/h)∫_{x−h/2}^{x+h/2} M = (M₁(x+h/2) − M₁(x−h/2))/h; apply B2 to both, and bound |(x+h/2)^{s+1} − (x−h/2)^{s+1}| ≤ min(|s+1| h X^{σ₁}, 2X^{σ₁+1}) (the first from ∫ of (s+1)w^s over [x−h/2, x+h/2]). With B3. ∎

**Lemma B5.** With C_∞, a, T₁ from Corollary A4, assuming 0 < a < 1 (true for every certificate used: 0.09 ≤ a ≤ 0.19, asserted in batch.py), and J₀ := 2∫₀^{T₁} |1/ζ(σ₁+it)| dt/|σ₁+it|, for every τ > 0:
  J(τ) ≤ J₀ + 2C_∞ (1/a + 1/(1−a)) τ^a.
Proof. The integrand is even in t. On [0, T₁] use min(·) ≤ 1/|s|. On [T₁, ∞) use |s|, |s+1| ≥ t and |1/ζ| ≤ C_∞t^a: ∫_{T₁}^∞ t^a min(1/t, τ/t²)dt equals (τ^a − T₁^a)/a + τ^a/(1−a) if τ ≥ T₁ and τT₁^{a−1}/(1−a) ≤ τ^a/(1−a) if τ < T₁; both are ≤ τ^a(1/a + 1/(1−a)). ∎

**Theorem B6 (assembly).** Assume 0 < a < 1. Let θ ∈ (0, 1) with σ₁ + a(1 − θ) − θ ≤ 0 (i.e. θ ≥ θ* := (σ₁ + a)/(1 + a)), c_h > 0, and x_a with x_a^{1−θ} ≥ c_h. Put
  C := c_h/4 + x_a^{−θ} + K(x_a)·J₀/(2π) + Q(x_a)·2C_∞(1/a + 1/(1−a))/(2π),
  K(x) = x^{σ₁−θ}(1 + (c_h/2)x^{θ−1})^{σ₁},  Q(x) = x^{σ₁+a(1−θ)−θ}(1 + (c_h/2)x^{θ−1})^{σ₁}(2/c_h + x^{θ−1})^a.
Then |M(x)| ≤ C x^θ for all x ≥ x_a, and |M(x)| ≤ max(C, x_a^{1−θ}) x^θ for all x ≥ 1.
Proof. Take h = c_h x^θ in B4 (h ≤ x < 2x since x^{1−θ} ≥ c_h). Then X^{σ₁}/x^θ = K(x) and X^{σ₁}τ^a/x^θ = Q(x) (use τ = 2x^{1−θ}/c_h + 1). By B5, |M(x)|/x^θ ≤ c_h/4 + x^{−θ} + K(x)J₀/2π + Q(x)·2C_∞(1/a+1/(1−a))/2π. Each of x^{−θ}, K, Q is nonincreasing in x (the exponents σ₁ − θ and σ₁ + a(1−θ) − θ are ≤ 0, and θ − 1 < 0). The second claim: |M(x)| ≤ x ≤ x_a^{1−θ} x^θ for 1 ≤ x ≤ x_a. ∎

Optimisation used (not part of the proof): the limit of the bound as x → ∞ at θ = θ* is c_h/4 + D c_h^{−a} with D = C_∞(1/a + 1/(1−a))2^a/π, minimal at c_h = (4aD)^{1/(1+a)}.

## 3. Certificates (what the scripts check)

Scripts in Riemann/Deep/Mertens78/ (copied from .scratch/mertens78/): cert.py (Arb), run.py (assembly), lp.py/mkcert.py (LP search for λ; floating point, untrusted), chk.py.
For each parameter set the following are computed or checked in ball arithmetic (precision 120 bits):
1. (C1), (C2). Harnack: analytic (above). LP certificates (finitely many decimal masses): (C2) directly; (C1) on |w| ≥ Y by the bracket (Y = the first of 2 max c_j + 50, 10³, …, 10⁷ where the bracket is positive)
   w²·π·[LHS − RHS] ≥ Σ_{λ_j>0} λ_j x_j w²/((|w| + c_j)² + x_j²) + Σ_{λ_j<0} λ_j x_j w²/(|w| − c_j)² − x₁,
   which is nondecreasing in |w| > max c_j (each fraction w²/((w+c)²+d²) increases, each w²/(w−c)² decreases), checked > 0 at |w| = Y; and on [−Y, Y] by adaptive bisection with ball evaluation of LHS − RHS (every ball must have positive lower end; start: 4000 equal pieces of [−Y₀, Y₀], Y₀ = min(Y, 2 max c_j + 50), and pieces [w, 1.002w] beyond).
2. κ, A, Λ, a, log C_∞ (Corollary A4 with T₁ = 1000).
3. J₀ = 2∫₀^{1000} |1/ζ(σ₁+it)|/|σ₁+it| dt by Arb ζ on t-balls of radius 0.01 (t < 100) and 0.005 (t ≥ 100): each ball gives a rigorous lower bound for |ζ|, and the integral is bounded by Σ (sup on ball) × length.
4. θ = θ* rounded up to 4 decimals; c_h from the optimisation; x_a = 10^{⌈log₁₀(k c_h)/(1−θ)⌉} for several k; C of Theorem B6 at x_a (upper end of the ball).

## 4. Self-review (two passes)

Pass 1 (signs and kernels):
- The factor 1/s: Lemma B4 keeps 1/|s| and 1/|s(s+1)| exactly; the 1/(s(s+1)) kernel gives absolute convergence, so no truncated Perron error and no horizontal segments at finite T remain.
- Horizontal segments: only in B2, where they vanish because a′ ≤ 0.76 < 2 uniformly on [σ₁, c].
- Near t = 0: 1/ζ is analytic at s = 1 (zero of 1/ζ), so moving from Re s = c > 1 to σ₁ < 1 crosses no pole; J₀ uses true values of |1/ζ(σ₁+it)| on [0, 1000] (max ≤ 5.7 for σ₁ ≥ 0.89 there, certified).
- Majorant: Rademacher's |1−s| denominator is cancelled by the factor (s − 1) of Z; at s = 1 the bound is continuous. The region σ > 1 + η is covered separately.
- L₀ uses |1/ζ| ≤ ζ(σ)/ζ(2σ) (not 1/ζ(σ)); L₁ uses |ζ| ≤ ζ(σ). Signs of λ_j: positive masses need −u ≤ L₀, negative masses need u ≤ L₁ (Theorem A3 proof).
- Herglotz (C2) term: c(Σλ_jx_j − x₁) ≥ 0 needs (C2); checked for each certificate.
Pass 2 (asymptotics and monotonicity):
- Hypotheses 0 < a < 1 (B5, B6) and σ₁ ≥ 0.89 (B2) were missing in the first draft; added after a reviewer pass (reviewer verdict on the draft: "FLAWED" only for these missing hypotheses; "the signed-weight inequalities are correct; B6's monotonicity claims are correct in the intended parameter range").
- Corollary A4's 1/t error bounds each use only |t + y_j| ≥ t − c_j ≥ t/2 (T₁ ≥ 2 max c_j asserted in code).
- Lemma B5 holds for all τ (both cases), so Theorem B6 needs no lower bound on τ.
- Theorem B6's monotonicity: Q contains x^{σ₁+a(1−θ)−θ}; the code asserts σ₁ + a(1−θ) − θ ≤ 0 in ball arithmetic.
- The bound |M(x)| ≤ x for x ≥ 1 is trivial (|μ| ≤ 1), used only for the "all x ≥ 1" form.


## 5. Certified results (results.json; every line produced by batch.py in ball arithmetic)

Columns: certificate (H = Harnack, ML = LP on lines σ ∈ {1.005, 1.01, 1.02, 1.05, 1.1, 1.3, 2}, S = LP on σ = 1.01), σ₁, η, Λ, a, θ (= θ* rounded up), C valid for x ≥ x_a, C′ valid for all x ≥ 1, crossover with Ramaré's 0.013x/log x.

| name | σ₁ | η | Λ | a | θ | C (x ≥ x_a) | C′ (x ≥ 1) | beats 0.013x/log x from |
|---|---|---|---|---|---|---|---|---|
| H093 (Harnack, σ₀ = 1.02) | 0.93 | 0.1 | 2.6364 | 0.18409 | 0.9409 | 6.10e5 (x ≥ 1e163) | 7.91e5 | 1e175.8 |
| ML093g5 | 0.93 | 0.1 | 2.1707 | 0.13171 | 0.9382 | 1.50e6 (x ≥ 1e160) | 1.84e6 | 1e174.0 |
| ML092g8 | 0.92 | 0.1 | 2.2051 | 0.13557 | 0.9296 | 1.53e8 (x ≥ 1e169) | 1.85e8 | 1e181.4 |
| ML092g8 | 0.92 | 0.05 | 2.2051 | 0.10544 | 0.9277 | 6.22e8 (x ≥ 1e172) | 7.54e8 | 1e185.2 |
| ML091g12 | 0.91 | 0.03 | 2.4622 | 0.11332 | 0.9192 | 2.94e12 (x ≥ 1e200) | 3.54e12 | 1e211.9 |
| S090g20 | 0.90 | 0.01 | 3.0255 | 0.13672 | 0.9121 | 2.75e20 (x ≥ 1e275) | 3.35e20 | 1e287.0 |
| S090g40 | 0.90 | 0.005 | 2.5473 | 0.10058 | 0.9092 | 1.03e26 (x ≥ 1e326) | 1.23e26 | 1e340.0 |
| S089g400 | 0.89 | 0.005 | 2.9526 | 0.12692 | 0.9024 | 5.60e110 (x ≥ 1e1173) | 6.70e110 | 1e1190.1 |
| S089g3000 | 0.89 | 0.005 | 2.3944 | 0.09063 | 0.8992 | 8.20e381 (x ≥ 1e3826) | 1.01e382 | 1e3847.6 |

Consequences |M(x)| ≤ x^{θ′} for x ≥ X (X = max(x_a, C^{1/(θ′−θ)}), best row):

| θ′ | 0.95 | 0.94 | 0.93 | 0.92 | 0.91 | 0.90 |
|---|---|---|---|---|---|---|
| X | 1e395 (ML092g8, η 0.05) | 1e600 (ML091g12) | 1e1142 (S090g20) | 1e2409 (S090g40) | 1e14573 (S089g400) | 1e477393 (S089g3000) |

J₀ ≤ 18.3–18.9 and max_{0≤t≤1000} |1/ζ(σ₁+it)| ≤ 4.72 (σ₁ = 0.93), 4.93 (0.92), 5.15 (0.91), 5.40 (0.90), 5.68 (0.89), all certified.
