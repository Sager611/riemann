# Quasi-RH for Grössencharacters of Q(√−3) and CM elliptic curves (2026-10-21)

Status: a new statement with a complete proof **relative to** OpenAI, "The Quasi-Riemann Hypothesis" (5 Oct 2026,
github.com/openai/math, preprints/The-Quasi-Riemann-Hypothesis-October-5-2026, file build/paper2.tex; below "paper2",
numbering of paper2.pdf). Like every current record, it is conditional on that unrefereed preprint.

## 1. Statement
K = Q(√−3), O = Z[ω], primary = ≡ 1 (mod 3), α(n) = n/|n| on primary n (a Hecke character mod 3, infinity type z/|z|).

**Theorem G.** For every unitary Hecke character ψ of K, of finite or infinite order and any infinity type,
L_K(s, ψ) ≠ 0 for Re s > 11/12 (for ψ = N^{it} the pole at s = 1 − it is allowed).

- **G1 (CM elliptic curves).** If E/Q has CM by an order of K (j(E) ∈ {0, 54000, −12288000}; e.g. y² = x³ + D for every
  D ∈ Q^×, y² = x³ + 1, the Fermat cubic x³ + y³ = 1), then L(E, s) ≠ 0 for Re s > 17/12. Proof: by Deuring,
  L(E, s) = L_K(s, ψ_E) with ψ_E algebraic of weight 1, so L(E, s) = L_K(s − 1/2, ψ_E N^{−1/2}) with a unitary ψ.
- **G2 (CM newforms).** A newform f of weight k ≥ 2 with CM by K has L(f, s) ≠ 0 for Re s > (k − 1)/2 + 11/12 (same proof).
- **G3 (routine).** For fixed m ≠ 0, Σ_{N p ≤ x} α(p)^m log N p ≪_{m,ε} x^{11/12+ε} (explicit formula): power-saving
  angular equidistribution of Eisenstein primes in fixed sectors.

Before this, fixed zero-free half-planes were known only for finite-order Hecke L-functions (OpenAI 7/8; Liu
2610.12234, 2610.26835; Xu–Sang 2610.15291; Cureton–Thorner 2610.00200: all finite order). Novelty search (scout,
2026-10-21: arXiv API "quasi-Riemann", web): no stated zero-free half-plane for Grössencharacters, CM elliptic curves
or any GL(2) L-function.

## 2. Reduction (Step 0)
Every unitary ψ is ν0·α^m·N^{it} with ν0 of finite order and m ∈ Z (a unitary character of C^× is
z ↦ (z/|z|)^m|z|^{2it}; dividing by α^m N^{it} leaves trivial infinity type). N^{it} only shifts s, and
L(s, ψ̄) = conj L(s̄, ψ) moves zeros to conjugate points. So it suffices to treat **ν = α^m ν0 with m ≠ 1**
(for m = 1 apply the result to ν̄, which has m = −1). Put **k := 1 − m ≠ 0**.

## 3. Proof: paper2 with α replaced by its k-th power
paper2's dual coefficient is a_ξ(n) = ᾱ(n)γ2(n)ξ(n) (TeX line 807). Define a^[k]_ξ := ᾱ^k γ2 ξ = α^m·a_ξ.
1. **§3 (zero-free from the mean square)** uses ν only as a unimodular multiplicative coefficient: (3.2) and the Mellin
   argument hold verbatim; L_K(s, ν) is entire for ν ≠ 1 (Hecke). So Prop. 3.1 for ν gives Theorem G.
2. **§4.** Apply (4.15) to ν0 and multiply by α^m(z1)·conj α^m(z2): the right side becomes Σ c_ξ a^[k]_ξ(z1)·conj a^[k]_ξ(z2).
   (4.6) holds for a^[k] (α is completely multiplicative on primary elements). Hence Prop. 4.5 holds with the family
   E^[k] = (4.9) with a^[k]. The finite ray class group now carries only ν0, G, R.
3. **§5 and §7: the family is closed with the same k.** (5.3) becomes T^[k](X;Ψ) with ᾱ(n)^k and ᾱ(b)^{3k}; (5.8) holds
   with μ(h)ᾱ(h)^{3k}Ψ(h)³ (complete multiplicativity). By (7.3), a^[k](u1)·conj a^[k](u2)·γ(χ_{u1}·conj χ_{u2}) =
   α^m(u1)·conj α^m(u2)·μ(u1)μ(u2)(ξG)(u1/u2), so the Möbius polynomial (7.1) gets the factor α^m(n); then the identity
   of Lemma 7.2, multiplied by α^m(z1)·conj α^m(z2), returns a^[k]. Everything else in §§5, 7 and App. B uses only
   |a| ≤ 1, (4.6), Lemma 4.4 and Lemma B.2.
4. **§6 and App. A.2: the one new input** is Prop. 6.2 / Lemmas 6.3–6.4 for T^[k]. Lemma 6.1 holds with ᾱ(nb³)^k.
   In A.2 replace ∂_z̄ by D_k := ∂_z̄^k (k ≥ 1; for k ≤ −1 use ∂_z^{|k|} and complex-conjugate everything below).
   - (a) Modes: D_k ĕ(ℓz) = (2πi ℓ̄)^k ĕ(ℓz), and ∫_0^∞ v^{μ−1}K_{1/3}(av)dv = 2^{μ−2}a^{−μ}Γ((μ+1/3)/2)Γ((μ−1/3)/2)
     (DLMF 10.43.19). With J_k(s) := ∫_0^∞ D_kΘ_Ψ(z,v)|_{z=0} v^{2s+k−2} dv one gets, for Re s > 1,
     J_k(s) = c_k (27/(2π)²)^s Γ(s + k/2 + 1/6)Γ(s + k/2 − 1/6) T^[k](s, Ψ), c_k ≠ 0 (k = 1 is (A.16)).
   - (b) **Lemma A (chain rule).** With (A.9), z' = −δ'/c − z̄/(c²(v²+|z|²)), v' = v/(N(c)(v²+|z|²)). As power series in
     (z, z̄): z' − z'_0 = −z̄/(c²v²)·(1 + |z|²/v²)^{−1}, z̄' − z̄'_0 = O(|z|), v' − v'_0 = O(|z|²). The z⁰z̄^k coefficient of
     Φ(z', v') comes only from the leading term of (z' − z'_0)^k, so for every real-analytic Φ(w, u) (θ∘H or its
     conjugate, whichever (A.8) uses) ∂_z̄^k[Φ(g^{−1}(z + a/c, v))]|_{z=0} = (−(cv)^{−2})^k (∂_w^k Φ)(−δ'/c, 1/(N(c)v)).
     (k = 1 is the paper's "∂z'/∂z̄ = −1/(c²v²), no height-derivative term".) For k ≥ 1 the horizontal derivative kills
     the constant terms (3^{5/2}/2)v^{2/3} of all cusp expansions (6.8), so J_k is entire, exactly as in A.2.
   - (c) Substituting v ↦ (N(c_h)v)^{−1} gives the functional equation
     J_k(s) = (−1)^k Σ_h c_F(h)κ(g_{1,h})(c̄_h/c_h)^k N(c_h)^{1−2s} J^∨_{h,k}(1 − s), with
     J^∨_{h,k}(s) := (A.17) with ∂_z → ∂_z^k and v^{2s−1} → v^{2s+k−2}; (c̄_h/c_h)^k is unimodular (k = 1: the paper's α(c_h)²),
     = c'_k (27/(2π)²)^s Γ(s + k/2 + 1/6)Γ(s + k/2 − 1/6) Σ_ℓ d_σ(ℓ) α(ℓ)^k ĕ(−δ'_h ℓ/c_h) N(ℓ)^{−s}.
     Stirling, Phragmén–Lindelöf and the contour shift run as in A.2 (T^[k] entire, no poles crossed).
   - (d) Result: (6.9) with d(ℓ)α(ℓ) replaced by d(ℓ)α(ℓ)^k and V♯ by
     V♯_k(x) = (1/2πi)∫_(0) V̂*(−t) Γ(k/2+2/3+t)Γ(k/2+1/3+t)/(Γ(k/2+2/3−t)Γ(k/2+1/3−t)) ((2π)⁴x/27)^{−t} dt (k = 1: (6.6));
     |C| ≪_k 1. The finite-place analysis (A.3)–(A.12) involves only the finite-order Ψ, so Lemma 6.3 is unchanged up
     to unimodular scalars. Lemma 6.4: coefficient bounds unchanged; (6.12) holds with k-dependent constants, since the
     first numerator pole is at t = −k/2 − 1/3 ≤ −5/6 and on −1/4 ≤ Re t ≤ A the Gamma quotient is ≪_k (1+|Im t|)^{4 Re t}.
   - (e) Lemma 6.6 and Prop. 5.2 are verbatim (A_{u,m} gets the unimodular factor α^k; |A_{u,m}| ≤ 27).
5. Prop. 5.1 (induction), Prop. 3.1 and §3 then give Theorem G. Constants depend on ν (through m, ν0, S), which the
   continuation argument of §3 allows.

**Why m = 1 needs the conjugation.** For k = 0 there is no horizontal derivative: J_0 keeps the constant terms
v^{2/3}, the radial cubic Gauss-sum series has the Kubota–Patterson pole, and the resulting main terms are not covered by
the proof of Prop. 5.2. Conjugating ν turns k = 0 into k = 2.

**Checks.** paper2 §§3–7 and App. A, B read in full (pdftotext; conjugation signs from the TeX source).
`angular_check.py` (mpmath, 3 s): the Bessel–Mellin identity holds to 1e−16; Lemma A holds to 1e−31 for k = 2, 3
at a generic cusp map. Not machine-checked.

## 4. Record hunt (exponent level; no improvement of B♭ = 0.874886647)
- Combining the two OpenAI methods is already done in the 7/8 paper: its (17.87) is paper2's Prop. 3.1, and its long-
  witness branch e(r) = (1+5r)/6 (Lemma 17.6) is paper2's sixth-power trick. paper2 alone is capped at 11/12; even a
  perfect mean square for H ≥ D^h gives only 6/7 (h = 6/7). Prime amplifiers in paper2 lose: slot cost 2z, gain z.
- **New lever:** the P_a^{1/6} term of the additive Gram bound (7/8 paper, Prop. 15.2). If it were P_a^θ, Liu's
  system gives B = 0.8724 (θ = 1/8), 0.8695 (1/12), 0.8637 (0), below the 13/15 floor, which assumes θ = 1/6
  (`gram_theta.py`, 30 s, 30 MB; the coarse grid gives 0.874917 at θ = 1/6). **No-go:** the exceptional frequencies
  k = a⁶e resonate exactly when the column character can be made principal. That needs ξ⁶ = 1 and an unweighted Gauss
  modulus, and the principal signal needs both (ξ(a)⁶ = 1 gives the ζ(6z) pole, §7.2; the unweighted s-sum gives
  ζ(w) at w = 1). A ξ of order 5, or a Möbius or nonprincipal weight on s, removes the loss and also kills the signal.
  Open lead: apply Cauchy–Schwarz only off the resonant subspace and bound the projection of the theta rows onto it.
- Other families: the sextic family is forced. Möbius absorption μγ₋₁ ∝ ᾱγ₂ needs n = 6 with the cubic theta.
  Quartic: g(χ^{−1}) = unit/g(χ), not a theta coefficient. Quadratic: self-dual, circular. A different base field gives
  the same exponents (Liu Part III). There is no positivity step to swap. Averaging in Z gives cancellation over bad
  rows only for separated zero ordinates; at heights ≤ Z^τ aligned zeros are the worst case.
- ζ′: ζ′(s) ≠ 0 for Re s < 1 − B♭, |t| ≥ t0 (Levinson–Montgomery with log ζ ≪ (log t)^{θ<1} right of B♭): routine.
- Xu–Sang refine the 7/8 *boundary* (to 0.8748867, same families). They keep the decay of rescaled prime tuples
  (ℓ up to 1/5 with no reflected excess), average the numerator over the detected class with the amplified fourth
  moment, and license κ in two stages. Their Cor. 9.1 is the Möbius-sum consequence.
