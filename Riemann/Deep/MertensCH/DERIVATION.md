# DERIVATION — |M(x)| ≤ c x^θ for all x ≥ 1 from RH to height 3·10^12, the half-plane Re s > 7/8 and Chirre–Helfgott

Agent mertensch, 2026-11-11. Companion of REPORT.md. Builds on Riemann/Deep/Mertens78/PROOF.md (cited as M78; Lemmas A1–A4, B1–B2 are used unchanged).
Notation: s = σ + it; M(x) = Σ_{n≤x} μ(n); Q(x) = #{n ≤ x squarefree}; R(x) = Q(x) − 6x/π²; Z(s) = (s−1)ζ(s)/s; u = log|Z|.

## 0. Inputs

- (I1) ζ(s) ≠ 0 for Re s > 7/8 (Lean-checked, Riemann/Deep/QRH; M78 (I1)).
- (I2) Rademacher's convexity bound (M78 (I2)): |ζ(s)| ≤ 3|1+s|/|1−s|·(|1+s|/2π)^{(1+η−σ)/2} ζ(1+η) for 0 < η ≤ 1/2, −η ≤ σ ≤ 1+η.
- (I3) |ζ(s)| ≤ ζ(σ), |1/ζ(s)| ≤ ζ(σ)/ζ(2σ) for σ > 1.
- (I4) Platt–Trudgian (Bull. LMS 53 (2021) 792–797): every zero ρ of ζ with 0 < Im ρ ≤ T_PT := 3 000 175 332 800 has Re ρ = 1/2.
- (I5) Chirre–Helfgott, arXiv 2511.14736v1, Cor. "mertensimple" (TeX l.389–395): |M(x)| ≤ 3x/(π·10^10) + 11.39√x for all x ≥ 1. Used as a black box (unrefereed).
- (I6) The nine Mertens78 certificates (M78 §2, §5; Riemann/Deep/Mertens78/all.json, c_*.json): each gives a line σ_H ∈ [0.89, 0.93], an exponent 0 < a < 1 and, by M78 Cor. A4 with T₁ := T′, |1/ζ(σ_H + it)| ≤ C_∞(T′) |t|^a for |t| ≥ T′.
- (I7) Arb ball arithmetic (python-flint 0.9).

Fixed parameters: T′ = 3·10^12, T₁ = 1000, σ_L = 0.85, σ₀′ = 1.02, η′ = 0.1, β = 3/(2π²). (T′ + σ₀′ − 1/2 < T_PT, asserted.)

## 1. Where Chirre–Helfgott's linear term comes from (why the literal route fails)

Their Theorem 1.1 (TeX l.246–290): (1/x^{1−σ})Σ_{n≤x} a_n n^{−σ} = a_∞ ι_{δ,σ} + δ Σ_{|Im ρ|≤T} Res(w_{δ,σ}(s)A(s)x^{s−1}) + ε, δ = π/2T, |ι| ≤ δ, w_{δ,σ}(s) = coth(δ(s−σ)) − tanh(δ(1−σ)).
- The weight is band-limited: 1_{[−1,1]}(t)·w(1 + iTt) is a Fourier transform; the smoothed sum uses A(s) only on 1 + i[−T, T].
- ι is the L¹ distance between the sharp cutoff and the band-limited approximant, multiplied by sup|a_n|. Their own remark (TeX l.5148): "1/T comes from the crude estimate |Σ_{N<n≤N+H} μ(n)| ≤ H applied to H ≫ N/T". So the zeros above T are paid for by a trivial bound for μ on a window of relative width ≈ π/2T around x.
- The continuation w(s) is periodic in t (period 2T) and does not decay, so ∫_{|t|>T} |w(s)x^{s−1}/ζ(s)| dt diverges on every line Re s = B′; and the window sum Σμ(n)(1_{n≤x} − φ(n/x)) has a kernel with a jump at n = x and only (log)^{−2} decay, so its Mellin transform cannot be moved off Re s = 0 either. The zeros above T can only be used after smoothing at a scale h ≪ x, i.e. with the M78 kernel. That is what §2–§3 do.
- CH used T = 10^10 + 1 (Platt's residue sums, TeX l.3240–3290, 11.35051 = δΣ|coth(δρ)/ζ′(ρ)|), not T_PT.

## 2. Lemmas

**Lemma H0 (squarefree unsmoothing).** |R(y)| ≤ 3√y for y > 0. For x ≥ 1, 0 < h ≤ x, X = x + h/2:
  |M(x) − (1/h)∫_{x−h/2}^{x+h/2} M(u) du| ≤ βh + 6√X,  β = 3/(2π²).
Proof. Q(y) = Σ_{d≤√y} μ(d)⌊y/d²⌋, so R(y) = −yΣ_{d>√y} μ(d)/d² − Σ_{d≤√y} μ(d){y/d²}; with Σ_{d≥m} d^{−2} < 1/m² + 1/m ≤ 2/m (m = ⌊√y⌋ + 1 > √y): |R(y)| ≤ 2√y + √y (for 0 < y < 1, |R(y)| = 6y/π² ≤ 3√y). For u > x, |M(u) − M(x)| ≤ Q(u) − Q(x) = (6/π²)(u − x) + R(u) − R(x) ≤ (6/π²)|u − x| + 6√X; the same for u < x. Average over u ∈ [x − h/2, x + h/2]: (6/π²)(h/4) = βh. ∎

**Lemma H1 (RH strip).** Let κ′ = (1/2 + η′)/2 = 0.3, A′ = log(3ζ(1+η′)) + κ′ log(2/3) (> 0), h′(s) = A′ + (1+κ′)log|1+s| − log|s|. Then u ≤ h′ on Re s > 1/2. For 1/2 < σ < σ₀′, t ≥ T₁ with t + σ₀′ − 1/2 < T_PT:
  |1/ζ(σ + it)| ≤ C_L(σ) t^{a_L(σ)},  Λ_L(σ) = (2σ₀′ − σ − 1/2)/(σ − 1/2),  a_L = (Λ_L − 1)κ′,
  log C_L(σ) = (Λ_L − 1)A′ + Λ_L log(ζ(σ₀′)/ζ(2σ₀′)) + [σ + Λ_L((1+κ′)(1+σ₀′) + 2)]/T₁.
Proof. u ≤ h′: M78 Lemma A1 with B = 1/2 (same proof; ℓ = log|1+s| > log(3/2), the affine map α ↦ α(ℓ − log 2π) − κ′(ℓ − log(3/2)) is ≤ 0 at α = 0 and α = κ′; for σ > 1 + η′, u ≤ log ζ(1+η′) < h′). v′ = h′ − u ≥ 0 is harmonic on the open disc D of centre s₀ = σ₀′ + it and radius ρ = σ₀′ − 1/2: D ⊂ {Re s > 1/2, |Im s| < T_PT}, where Z has no zeros by (I4) (and Z(1) = 1). Harnack in D at distance r = σ₀′ − σ: v′(σ + it) ≤ ((ρ + r)/(ρ − r)) v′(s₀) = Λ_L v′(s₀). So −u(σ + it) ≤ Λ_L(h′(s₀) − u(s₀)) − h′(σ + it) ≤ Λ_L(h′(s₀) + L₀(s₀)) − h′(σ + it) by (I3), and log|1/ζ| = −u + log|(s−1)/s|. This is the inequality of M78 Theorem A3 for the single mass λ = Λ_L at (σ₀′, 0); the elementary estimates of M78 Corollary A4 with c_j = 0 (valid for t ≥ T₁ ≥ 1) turn it into C_L(σ)t^{a_L(σ)}. ∎
Certified values (Arb, logs/cert_ch.log): Λ_L(0.85) = 1.9714, a_L = 0.29143, log C_L = 10.048.

**Lemma H2 (contour).** Let σ_H be the line of a certificate (I6), σ_L = 0.85 < σ_H, and Γ = {σ_L + it: |t| ≤ T′} ∪ {σ ± iT′: σ_L ≤ σ ≤ σ_H} ∪ {σ_H + it: |t| ≥ T′}. For y > 0, M₁(y) := ∫₀^y M = (1/2πi)∫_Γ y^{s+1} ds/(s(s+1)ζ(s)).
Proof. M78 Lemma B1 gives the formula on Re s = c > 1. 1/ζ(s) y^{s+1}/(s(s+1)) is analytic on the closed region between Re s = c and Γ: for |t| ≤ T′ it lies in Re s ≥ σ_L > 1/2, zero-free by (I4) (1/ζ is analytic at s = 1); for |t| ≥ T′ in Re s ≥ σ_H > 7/8 (I1). Cauchy on [σ_L, c] × [−T′, T′]; for |t| ≥ T′ the horizontal integrals at height ±U → ∞ vanish exactly as in M78 Lemma B2 (σ_H ≥ 0.89). ∎

**Lemma H3 (one x).** For x ≥ 2·10^8, 0 < h ≤ x, X = x + h/2, τ = 2X/h:
  |M(x)| ≤ βh + 6√X + (X^{σ_L}/2π)(J_L0 + 2C_L T′^{a_L}/a_L) + (σ_H − σ_L)C_L* T′^{a_L} X^{σ_H} min(1/T′, τ/T′²)/π + (C_∞(T′)/π) X^{σ_H} G(τ, T′),
where a_L, C_L are at σ_L; J_L0 := 2∫₀^{T₁} |1/ζ(σ_L+it)|/|σ_L+it| dt ≤ 19.5963 (Arb, M78 J0_direct at σ = 0.85, 95 000 balls); log C_L* is log C_L(σ_L) with the first σ in the bracket replaced by σ_H; and
  G(τ, T′) = τT′^{a−1}/(1−a) (τ ≤ T′),  (τ^a − T′^a)/a + τ^a/(1−a) (τ ≥ T′);  G(τ, T′) ≤ τ^a(1/a + 1/(1−a)) for all τ > 0.
Proof. Lemma H0 and (1/h)∫M = (M₁(x+h/2) − M₁(x−h/2))/h with Lemma H2 and |(x+h/2)^{s+1} − (x−h/2)^{s+1}| ≤ min(|s+1|hX^σ, 2X^{σ+1}) (M78 B4): |M(x)| ≤ βh + 6√X + (1/2π)∫_Γ X^{Re s}|1/ζ(s)| k(s)|ds|, k(s) = min(1/|s|, τ/(|s||s+1|)); the integrand is even in t.
- Re s = σ_L, |t| ≤ T′: k ≤ 1/|s|; on [0, T₁] this is J_L0/2, on [T₁, T′] Lemma H1 gives ≤ C_L∫t^{a_L−1} ≤ C_L T′^{a_L}/a_L; both signs of t: factor 2.
- Horizontal segments: |s|, |s+1| ≥ T′, so k ≤ min(1/T′, τ/T′²); X^σ ≤ X^{σ_H}; by Lemma H1, |1/ζ(σ + iT′)| ≤ C_L(σ)T′^{a_L(σ)} ≤ C_L* T′^{a_L(σ_L)} because Λ_L decreases in σ, A′ > 0, log(ζ(σ₀′)/ζ(2σ₀′)) > 0, σ ≤ σ_H, T′ ≥ 1; two segments of length σ_H − σ_L.
- Re s = σ_H, |t| ≥ T′: (I6) and |s|, |s+1| ≥ t: (1/2π)·2∫_{T′}^∞ X^{σ_H}C_∞ t^a min(1/t, τ/t²)dt = (C_∞/π)X^{σ_H}G(τ, T′) (split at t = τ if τ > T′). The bound G ≤ τ^a(1/a + 1/(1−a)): for τ ≤ T′, (τ/T′)^{1−a} ≤ 1; for τ ≥ T′ drop −T′^a/a. ∎

**Lemma H4 (a range).** Fix θ with σ_L ≤ θ < 1, a certificate, c_h > 0 and x_k < x_{k+1} with c_h ≤ x_k^{1−θ}, x_k ≥ 2·10^8. Take h = c_h x^θ (so h ≤ x). Then |M(x)| ≤ c_k x^θ on [x_k, x_{k+1}] with
  c_k = βc_h + 6 x_k^{1/2−θ}r^{1/2} + x_k^{σ_L−θ}r^{σ_L}(J_L0 + 2C_LT′^{a_L}/a_L)/2π + (σ_H−σ_L)C_L*T′^{a_L} D min(1/T′, τ₁/T′²)/π + (C_∞/π) D G(τ₁, T′),
  r = 1 + (c_h/2)x_k^{θ−1},  D = x_*^{σ_H−θ} r^{σ_H} (x_* = x_k if σ_H ≤ θ, else x_{k+1}),  τ₁ = 2x_{k+1}^{1−θ}/c_h + 1.
Proof. Divide Lemma H3 by x^θ: X/x = 1 + (c_h/2)x^{θ−1} decreases in x; x^{1/2−θ}, x^{σ_L−θ} decrease; τ = 2x^{1−θ}/c_h + 1 increases and G, min(1/T′, τ/T′²) increase in τ. Each term is a product of monotone factors, each bounded at its worst endpoint. ∎

**Lemma H5 (final range).** If moreover σ_H + a(1−θ) − θ ≤ 0 and σ_H ≤ θ, then |M(x)| ≤ c_fin x^θ for all x ≥ X₁ (c_h ≤ X₁^{1−θ}), where c_fin = βc_h + 6X₁^{1/2−θ}r^{1/2} + X₁^{σ_L−θ}r^{σ_L}(J_L0 + 2C_LT′^{a_L}/a_L)/2π + (σ_H−σ_L)C_L*T′^{a_L}X₁^{σ_H−θ}r^{σ_H}/(πT′) + (C_∞/π)(1/a + 1/(1−a))X₁^{σ_H+a(1−θ)−θ} r^{σ_H}(2/c_h + X₁^{θ−1})^a, r at X₁.
Proof. Lemma H3 with G ≤ τ^a(1/a + 1/(1−a)), τ^a = x^{a(1−θ)}(2/c_h + x^{θ−1})^a and min(·) ≤ 1/T′: every term divided by x^θ is nonincreasing in x (M78 Theorem B6 argument). ∎

**Theorem.** For each θ in the table of REPORT.md, |M(x)| ≤ c x^θ for every x ≥ 1, with c = max(c_CH, max_k c_k, c_fin), where c_CH = max(3/(π10^10) + 11.39, 3X₀^{1−θ}/(π10^10) + 11.39X₀^{1/2−θ}) covers [1, X₀] by (I5) (the ratio CH(x)/x^θ is a sum of two exponentials in log x, hence convex, so its maximum on [1, X₀] is at an endpoint), the ranges [x_k, x_{k+1}] tile [X₀, X₁] (Lemma H4) and Lemma H5 covers [X₁, ∞).
The ranges are contiguous, x₀ = X₀ and x_K = X₁, and X₁ = X₀ when there are none (asserted). The parameters (X₀, X₁, every x_k = 10^{n_k/10^6} with steps 0.05 in log₁₀ x up to 10^3000 and 1 beyond, certificate and c_h per range, final certificate and c_h) are chosen by floating-point search; every term is evaluated in Arb (cert_ch.py); c is rounded up to 3 digits and re-checked against every piece in Arb.

## 3. What limits the result

- At the hand-over x ≈ X₀ (≈ 10^245) the optimal τ ≈ 4·10^10, far below T′; the hybrid stays effective until τ ≈ T′, i.e. x^{1−θ} grows by ≈ 70, then the line σ_H carries all of t ≤ τ again and the constant is M78's. For x → ∞ nothing changes: θ ≥ θ* of a half-plane certificate is still required.
- The RH-strip pieces are negligible: with the logged constants (e^{20.35} for J_L0 + 2C_LT′^{a_L}/a_L, C_L* ≈ e^{10.05}, T′^{a_L} ≈ 10^{3.64}) one gets T_d ≤ 10^{−13} for x ≥ 10^{240} and T_e ≤ 10^{−11} T_c, so the bound is βh + (C_∞/π)X^{σ_H}G(τ, T′) up to a relative 10^{−10} (estimate; the certificate includes both terms anyway). A larger verified height T would help as T^{−(1−a)/2} in the hybrid window.

## 4. Routine rerun with the Lean-checked half-plane B′ = 0.87338196661320752 (bprime.py)

- (I1′) ζ(s) ≠ 0 for Re s > B′ (TREE N7, Lean). Lemma A1 holds with B′: κ = (1 + η − B′)/2, A = log(3ζ(1+η)) + κ log(1/(1+B′)) (same proof, ℓ > log(1 + B′)).
- Translation: give every certificate the lines σ_j′ = σ_j − D and σ₁′ = σ₁ − D, D = 7/8 − B′ = 0.00161803338679248. The heights σ_j′ − B′ = σ_j − 7/8 and the shifts y_j are unchanged, so (C1) and (C2) are literally the same inequalities (checked anew for H093, ML093g5, ML092g8). Lemma A2 on Re s > B′, Theorem A3 and Corollary A4 then hold with the new κ, A and ζ-values at σ_j′ (all σ_j′ > 1).
- Contour shift (M78 B2) with B′: Harnack at σ₀ = 1.02, η = 0.05 gives the uniform exponent a′ = ((1.02 − B′)/(σ₁′ − B′) − 1)(1.05 − B′)/2 < 2 for every σ₁′ used (asserted).
- Effect: Λ is unchanged, κ grows by D/2, σ₁ drops by D; θ* drops by about 0.0014. The constants change through ζ(σ_j′) and J₀ at σ₁′ (recomputed, Arb). This is the "plain routine gain" asked for in the brief.

## 5. Self-review

- Kernel bound: |(x+h/2)^{s+1} − (x−h/2)^{s+1}| ≤ min(|s+1|hX^σ, 2X^{σ+1}) needs σ ≥ 0; Γ has σ ≥ σ_L = 0.85.
- Γ avoids zeros: |t| ≤ T′ < T_PT with σ ≥ σ_L > 1/2; |t| ≥ T′ with σ ≥ σ_H > 7/8. No residue: 1/ζ is analytic at s = 1; s = 0, −1 lie left of Γ.
- Harnack disc of Lemma H1: centre σ₀′ + it, radius σ₀′ − 1/2 = 0.52, used only for T₁ ≤ t ≤ T′, so the disc lies in 1/2 < Re s, |Im s| < T′ + 0.52 < T_PT.
- Monotonicity in Lemma H4: each factor is evaluated at its worst endpoint; τ₁ uses x_{k+1}; D uses x_{k+1} when σ_H > θ (happens for θ < 0.93).
- Rounding: c is the Arb upper bound of the maximum of all pieces, rounded up to 3 digits and re-compared in Arb (cert_ch.py, certify()).
- Reviewer pass (reviewer subagent, 2026-11-11, on §0–§5 and certify()): FLAWED for two missing hypotheses only, θ < 1 in Lemmas H4–H5 and "no gap when there are no ranges"; both are now stated and asserted. H0, H1, H3 and the implementation of H4–H5 (endpoints, factor 1/π, G and its bound, rounding) were confirmed.
- Not re-derived here: M78 Lemmas A1–A4, B1–B2 (reviewed in M78) and the (C1) bisection checks of the LP certificates (Mertens78/verify.log); the black boxes (I4), (I5).
