NO DISCOVERY: Stenger's §4 step (Lemmas 4.2/4.3) is a false complex Rolle theorem; a standard objection, not new.

# Audit: Stenger, arXiv 1708.01209v19 (July 2024; also Adv. Pure Math. 2023)

**Verdict: the proof fails at §4, the one sentence used in both Lemma 4.2 and Lemma 4.3.** Everything in §2 is correct (or has harmless slips).

## Chain
1. (1.2)/(2.1) G(z) = ∫₀^∞ y^{z−1}/(1+e^y) dy = (1−2^{1−z})Γ(z)ζ(z). **Correct** (checked with mpmath to 8 digits).
2. Lemma 2.6(i): G and ζ have the same zeros in D. **Correct.** Slips: the (−1)^n in (2.9) should not be there; the (ii) argument uses "G(z) = G(1−z)" but the functional equation has the factor K. Harmless, since the zeros are isolated anyway.
3. Lemmas 2.9 and 2.10: κ(σ,x) = e^{σx}/(1+exp eˣ) > 0, so ∂_σ^{2m}G(σ) = ∫x^{2m}κ > 0 on (0,1), and ∂_σ^{2m+1}G(σ) = −∫₀^∞x^{2m+1}Δ < 0 on (0,1/2]. **Correct.** The garbled Taylor line (2.15) is fixable: δ(1/2,y) = Σ_{n≥1}(y^{1−2n} − y^{2n−1})/n! > 0. Numbers for m ≤ 4: `.scratch/stenger/num.py`.
4. **§4 (fails).** Quote, Lemma 4.3 I.1: "each of G(s ± iτ₁) vanishes with multiplicity ν₁+1 … G vanishes on ℓ[s+iτ₁, s−iτ₁] with total multiplicity 2ν₁+2. Equivalently, G^{(2ν₁+2)}(s) vanishes with multiplicity one, but which is not allowed by Lemma 2.10." Lemma 4.2 uses the same sentence at s = 1/2.

## Why step 4 is false
- It needs Rolle's theorem for a complex-valued function on a segment, and that theorem is false: e^z − 1 vanishes at 0 and 2πi, but its derivatives never vanish.
- Rolle applied to Re G (even in t) gives nothing new. It only gives an odd t-derivative vanishing at t = 0, which Lemma 2.10(iii) already says is identically 0. The count is also off by one: N zeros give the (N−1)-th derivative, not the N-th.
- Even a zero of G^{(N)} somewhere on the segment would not touch Lemma 2.10, which only speaks about the real point s.
- **It proves too much.** With ν = 0 the step forbids simple zeros too, including those on L (Lemma 4.2 at s = 1/2). So it would rule out the first zero 1/2 + 14.1347i, and with it the Hardy–Littlewood result the paper itself cites as [16].

## Davenport–Heilbronn test: the argument fails it
DH = Σaₙn^{−s} with coefficients (1, κ, −κ, −1, 0) repeating mod 5 (κ ≈ 0.284), real on ℝ. Γ·DH = ∫y^{s−1}θ(y)dy, where θ = (u+κu²−κu³−u⁴)/(1−u⁵), u = e^{−y}, factors as u(1−u)(1+(1+κ)u+u²)/(1−u⁵) > 0. That is the same positivity as 1/(1+e^y). Numerically, the sign table of Lemma 2.10 holds for Γ·DH (σ = 0.2, 0.5, 0.8; derivatives 0..5). Nothing in §4 uses the Euler product, so §4 would "prove" RH for DH, which has a certified zero at 0.808517 + 85.699348i.

## Lean: `Riemann/Deep/Stenger/Paper.lean`
Compiles with `lake env lean`; no `sorry`; axioms are only propext, Classical.choice and Quot.sound.
- `EvenPos H` is Lemma 2.10(i). `StepOn a H` is the §4 step (a zero at σ+iτ with σ ≤ a forces some H^{(2m)}(σ) = 0, m ≥ 1). `StepHalf` is the paper's range (0, 1/2].
- `chain`: EvenPos + Step ⇒ no zero off the real axis at all (it proves too much).
- `not_stepHalf_cosh`: **unconditional counterexample.** cosh(z − 3/10) is real on ℝ, has all even derivatives cosh(σ − 3/10) > 0, and vanishes at 3/10 + iπ/2.
- `stenger_refuted : EvenPos G → ¬ StepHalf G`, via the certified first zeta zero `Fayed.FirstZero.first_zero`.
- `theta_pos`: the DH kernel is positive. `dh_test : EvenPos GD → ¬ StepOn (9/10) GD`, via `Fayed.DH.Zero.offLineZero_holds`.
- Not formalized: EvenPos for G and for Γ·DH. Both are true (the kernel is positive), and the numbers above check them; the Mellin-integral proof is not in Lean.

## Prior critiques
Earlier versions (v9: "Schwarz reflection, Rolle's theorem and the Cauchy–Riemann equations") were dismissed informally, e.g. a 2017 Zhihu thread that points to Davenport–Heilbronn. The complex-Rolle failure is folklore. No formal refutation was found; it is not a significant result. The author was not contacted.
