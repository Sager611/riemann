# Quasi-RH for Hecke characters of Q(√−3) with nonzero angular type (2026-10-21, repaired and novelty-checked)

**Status: KNOWN, and superseded.** The statement below, its elliptic-curve case and its proof route were already on
Hexagon on 13 Oct 2026:
- Rubinstein–Restrepo Borrero, Hexagon 2610.00220v1, Part I: OpenAI's 11/12 argument with ∂_z̄²θ gives
  L(E_d, s) ≠ 0 for Re s > 17/12, E_d: y² = x³ + d (our G1 for j = 0, same route). Part II: L(s, ᾱη) ≠ 0 for
  Re s > 7/8 for every finite-order η of K, so L(E,s) ≠ 0 for Re s > 11/8 for every E with CM by K.
- Hu, Hexagon 2610.00216v1, Thm 1.2: L_F(s, ψ) ≠ 0 for Re s > 7/8 for every unitary Hecke character ψ of every
  number field, infinite order included (pole only for ψ = N^{it}, at 1 + it); §9: CM motives at w/2 + 7/8.
Hu's 7/8 implies Theorem G, G1, G2 and G3 below with 11/12 replaced by 7/8. What remains is a corrected,
independent write-up of the general-k step of the 11/12 route (R–RB §17 lists the general-k step for their 7/8 route
as "what would still have to be written"; Hu reaches all types by another route). Everything is relative to unrefereed
preprints: OpenAI, "The Quasi-Riemann Hypothesis" (5 Oct 2026, github.com/openai/math, build/paper2.tex; "paper2",
numbering of paper2.pdf).

## 1. Statement
K = Q(√−3), O = Z[ω], λ = √−3, primary = ≡ 1 (mod 3), α(n) = n/|n| on primary n (a Hecke character mod 3).
Every unitary Hecke character ψ of K has a unique form ψ((a)) = ν0((a))·α(a)^m·N(a)^{it} (a primary, prime to the
conductor) with ν0 of finite order, **angular type m ∈ Z**, t ∈ R.

**Theorem G (m ≠ 0).** If m ≠ 0, L_K(s, ψ) is entire and L_K(s, ψ) ≠ 0 for Re s > 11/12.
- m = 0 is paper2 Thm 1.1 for ψN^{−it}, shifted by it. If moreover ν0 is principal, L_K(s, ψ) = ζ_K(s − it) times a
  finite Euler product: a simple pole at s = 1 + it (the draft said 1 − it), no zeros in Re s > 11/12; paper2 §3
  handles it meromorphically (L^S·M_W = Ŵ holds on Re s > 11/12 off the pole).
- Scope: finite-order characters already give GL(2) cases (L_K(s,χ) = L(s, θ_χ), θ_χ a dihedral weight-1 form; paper2
  and Liu) and N^{it} twists give trivial infinite-order ones. The only content beyond paper2 is m ≠ 0, i.e. CM forms of
  weight ≥ 2 and CM elliptic curves. The draft's "first for any GL(2)" was false.
- **G1.** E/Q with CM by an order of K (j ∈ {0, 54000, −12288000}): L(E,s) = L_K(s − 1/2, ψ), m = ±1 (Deuring), so
  L(E,s) ≠ 0 for Re s > 17/12. Known: R–RB Part I (j = 0, 17/12); R–RB Part II and Hu (all three j, 11/8).
- **G2.** A newform of weight w ≥ 2 with CM by K: L(f,s) = L_K(s − (w−1)/2, ψ), m = ±(w−1), so L(f,s) ≠ 0 for
  Re s > (w−1)/2 + 11/12. Known stronger: Hu, (w−1)/2 + 7/8.
- **G3 (routine).** Σ_{N𝔭≤x} α(𝔭)^m log N𝔭 ≪_ε x^{11/12+ε} for 0 < |m| ≤ x (explicit formula; the zero-free region is
  uniform in m), hence Eisenstein primes in sectors with error x^{11/12+ε}. Known stronger via Hu (7/8).

## 2. Reduction (Step 0, corrected)
L(s, ψ) = L(s − it, ψN^{−it}), and L(s, ψ̄) = conj L(s̄, ψ) with ψ̄ of type −m. So it suffices to treat
**ν = α^m ν0 with m ≥ 1, and k := 1 + m ≥ 2**. The draft had k = 1 − m: a sign error (see §3.2).

## 3. Proof: paper2 with ᾱ replaced by ᾱ^k
paper2's dual coefficient is a_ξ = ᾱγ₂ξ (TeX l. 807). Put a^[k]_ξ := ᾱ^k γ₂ ξ = ᾱ^m a_ξ.
1. **§3.** The Mellin argument is verbatim: L^S(s, ν) is entire for m ≠ 0 (Hecke), the Euler factors at S are nonzero
   on Re s > 0, so Prop. 3.1 (mean square (3.2)) for ν gives Theorem G. Constants may depend on m, ν0, S.
2. **§4, (4.15).** Expanding |Σ μνχ_n(u)W|² puts ν̄(z₁)ν(z₂) = ᾱ^m(z₁)α^m(z₂)·ν̄0(z₁)ν0(z₂) on the left of (4.15).
   Apply (4.15) to ν0 (expand ν̄0(t)G(t⁻¹)) and multiply by ᾱ^m(z₁)·conj ᾱ^m(z₂): the right side becomes
   Σ c_ξ a^[k]_ξ(z₁)·conj a^[k]_ξ(z₂), k = 1 + m. (The factor ᾱ(z₁)α(z₂) inside a_ξ(z₁)·conj a_ξ(z₂) comes from (4.5) and
   (4.4).) (4.6) and |a^[k](v)χ_v(he⁵)|² = 1_{(v,h)=1} hold (ᾱ^k completely multiplicative), so Prop. 4.5 reduces
   Prop. 3.1 to (4.12) for E^[k] := (4.9) with a^[k].
3. **§5, §7: closed with the same k.** (5.3) becomes T^[k](X; Ψ) with ᾱ(n)^k ᾱ(b)^{3k}; (5.8) holds with
   μ(h)ᾱ(h)^{3k}Ψ(h)³, and β₀ in Lemma 5.3 gets ᾱ(b)^{3k}. (7.3) becomes a^[k](u₁)·conj a^[k](u₂)·γ(χ_{u₁}·conj χ_{u₂})
   = ᾱ^m(u₁)α^m(u₂)μ(u₁)μ(u₂)(ξG)(u₁/u₂). So p_y(n) in (7.1) gets ᾱ^m(n), and the t-factors cancel. In Lemma 7.2 the
   first column then carries ξ₁ᾱ^m, and its identity times ᾱ^m(z₁)α^m(z₂) returns a^[k]_{ξ'}. Lemma 7.3, Lemma 4.4 and
   App. B use only |a| ≤ 1 and (4.6); α does not occur in §7 or App. B.
4. **§6 and App. A.2: the one new input** is Prop. 5.2 for T^[k]. Lemma 6.1 holds with ᾱ(nb³)^k. Replace ∂_z̄ by ∂_z̄^k.
   - (a) With ∂_z̄^k ĕ(ℓz) = (2πiℓ̄)^k ĕ(ℓz), DLMF 10.43.19 and ᾱ(λ⁻³)^k·i^k = 1:
     J_k(s) := ∫_0^∞ ∂_z̄^kΘ_Ψ|_{z=0} v^{2s+k−2} dv = (3^{5/2}/4)(27/(2π)²)^s Γ(s+k/2+1/6)Γ(s+k/2−1/6) T^[k](s,Ψ)
     for Re s > 1 (k = 1 is (A.16)).
   - (b) **Lemma A (chain rule).** With (A.9), z' − z'₀ = −z̄/(c²v²)·(1 + |z|²/v²)⁻¹, z̄' − z̄'₀ = O(|z|) and
     v' − v'₀ = O(|z|²), each a power series in z and z̄. Only the leading term of (z' − z'₀)^k contributes to the z⁰z̄^k
     coefficient, so for Φ = θ̄∘H_h, ∂_z̄^k[Φ(g⁻¹(z + a/c, v))]_{z=0} = (−(cv)⁻²)^k (∂_w^kΦ)(−δ'/c, 1/(N(c)v)).
     For k ≥ 1 this kills the constant terms (3^{5/2}/2)v^{2/3} of (6.8), so J_k is entire, as in A.2.
   - (c) Substituting v ↦ (N(c_h)v)⁻¹ and using (A.8):
     J_k(s) = (−1)^k Σ_h c_F(h)·**conj κ(g_{1,h})**·**ᾱ(c_h)^{2k}**·N(c_h)^{1−2s} J^∨_{h,k}(1 − s), with
     J^∨_{h,k}(s) := ∫_0^∞ ∂_z^k{θ̄(H_h(z,v))}|_{z=−δ'_h/c_h} v^{2s+k−2} dv
     = (i^k/4)**(2π)^{−2s}**Γ(s+k/2+1/6)Γ(s+k/2−1/6) Σ_{0≠ℓ∈λ⁻⁴O} d_σ(ℓ)α(ℓ)^k ĕ(−δ'_hℓ/c_h)N(ℓ)^{−s}.
     k = 1 gives (A.18), (A.19). Three draft errors are fixed here. The draft had κ; it must be conj κ, the conjugate
     multiplier of (A.8), which (A.12) combines with C_{p,j} to give B_{p,1} = χ_p³, the quadratic factor needed for the
     large sieve. The draft had α(c_h)²; it must be ᾱ(c_h)² = c̄_h/c_h. The draft had (27/(2π)²)^s; it must be (2π)^{−2s}.
     Stirling, Phragmén–Lindelöf and the contour shift run as in A.2 (T^[k] is entire, so no poles are crossed).
   - (d) With t = 1/2 − s: (6.9) holds with d(ℓ)α(ℓ)^k, C_k = ((−i)^k/81)ᾱ(c)^{2k}φ̂(h₀)·conj κ₀·∏(as in A.2), |C_k| ≤ 1/81,
     V♯_k(x) = (1/2πi)∫_(0) V̂*(−t) Γ(k/2+2/3+t)Γ(k/2+1/3+t)/(Γ(k/2+2/3−t)Γ(k/2+1/3−t)) ((2π)⁴x/27)^{−t} dt
     (k = 1: (6.6)). The (2π), 27 and N(c) powers do not depend on k. The first numerator pole is at
     t = −k/2 − 1/3 ≤ −4/3, so (6.12) holds with k-dependent constants (the Gamma quotient is ≪_k (1+|Im t|)^{4 Re t}).
     The finite places ((A.3)–(A.12), Lemma 6.3) are unchanged: the derivative acts only at ∞.
   - (e) Lemma 6.6: A_{u,m'}(n,b) gets the factor α(uλ^{m'}nb³)^k and still |A| ≤ 27; (6.13) allows arbitrary
     coefficients. The proof of Prop. 5.2 is verbatim.
5. Prop. 5.1 (induction), Prop. 3.1 and §3 then give Theorem G.

**Agreement with R–RB Part I (k = 2).** Their (8.4) is Lemma A at k = 2 (factor (cv)⁻⁴). Their reflected kernel
Γ(4/3+t)Γ(5/3+t)/(Γ(4/3−t)Γ(5/3−t)) is V♯₂. Their one-mode Mellin carries (2π)^{−2s}, their check (B.16) uses conj κ,
and their coefficient is "cubic Gauss sums times (n̄/|n|)²" = a^[2]. This matches every repaired point.

**Why m < 0 is conjugated.** Then k = 1 + m ≤ 0. For k = 0 (m = −1) there is no horizontal derivative: J₀ keeps the
constant terms v^{2/3}, the radial cubic Gauss-sum series has the Kubota–Patterson pole, and the proof of Prop. 5.2 does
not cover its main terms. Negative k would need ∂_z^{|k|} formulas that were never derived. Conjugation avoids both.

**Checks.** paper2 §§3–7 and App. A.2 were re-read line by line against the TeX (numbering checked in the pdf).
`angular_check.py` (mpmath; 16 s, 22 MB) checks four things, each to ≤ 1e−16:
(1) the Bessel–Mellin identity;
(2) Lemma A for k = 1–4, and that the constant term is killed;
(3) the one-mode forward and dual Mellin constants (c_k = 3^{5/2}/4, dual factor (i^k/4)(2π)^{−2s});
(4) the one-mode functional equation in (c) for k = 1, 2, 3.
None of this is machine-checked.

**Novelty search (2026-10-21).** Sources: arXiv HTML search (the API was rate-limited); the Hexagon API; OpenAI's list
of 749 preprints (7/8, 11/12, Landau–Siegel, cubic first moment, Gaussian primes, modularity); Liu 2610.12234 and
2610.26835; Xu–Sang 2610.15291; Cureton–Thorner (Hexagon 2610.00200v2); FSW 2610.18491; the web; the evidence API.
- Liu, Xu–Sang, Cureton–Thorner and FSW: finite order only (Liu's angular Λ is a tool, not a target).
- OpenAI's cubic-first-moment preprint has angular prime sums but no half-plane.
- Hits: Hu 2610.00216v1 and R–RB 2610.00220v1 (both Hexagon, 13 Oct 2026).

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
