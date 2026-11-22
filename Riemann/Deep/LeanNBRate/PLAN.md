# LeanNBRate — PLAN (lemma DAG, Lean interfaces)

Target: discoveries/nyman-beurling-continuous-rate.md, Claim (a)–(c). Proof: DbarNB/DERIVATION.md [N], ExploreDbar/REPAIR.md [R],
ExploreDbar/DERIVATION.md [X], audited chain DbarAudit/FULLPROOF.md [F]. Restructured below for Lean; every deviation from [N]/[R]
is a re-proof, not an input. Notation: L = log(1/λ), ℓ = L − √L, s = ½ + it on the line, w(t) = 1/(2π(¼+t²)),
‖F‖² = ∫|F|²w, N_ℓ(p) = ∫_0^ℓ(|p′|² + ¼|p|²) + ½|p(0)|² + ½|p(ℓ)|².

## 0. Statement (Defs.lean, Main.lean)
- volPos = volume.restrict (Ioi 0); L2pos = Lp ℝ 2 volPos; nbFun θ x = Int.fract (θ/x); chiVec = class of 1_(0,1).
- nbSet λ = {f | ∃ θ ∈ Icc λ 1, f =ᵐ nbFun θ}; nbDist λ = infDist chiVec (span ℝ (nbSet λ)).
- nymanSet λ = {f | ∃ n c θ, (∀k, θ k ∈ Icc λ 1) ∧ Σ c k θ k = 0 ∧ f =ᵐ Σ c k nbFun (θ k)}; nymanDist λ = infDist chiVec (nymanSet λ).
- NTZ = {ρ // ζ ρ = 0 ∧ 0 < re ρ < 1}; mult ρ = analyticOrderNatAt riemannZeta ρ; S = Σ' m²/‖ρ‖².
- MainClaim: RH → Summable ∧ Tendsto (log(1/λ)·nbDist²) (𝓝[>]0) (𝓝 S) ∧ ((∀ρ, mult ρ = 1) → S = 2+γ−log 4π)
  ∧ Tendsto (log(1/λ)·nymanDist²) (𝓝[>]0) (𝓝 S).

## 1. Assembly (Main.lean)
- lower: Burnol (input, complex L²) + nbDist_R = nbDist_C (ours) ⇒ liminf ≥ S; nymanDist ≥ nbDist (nymanSet ⊆ span).
- upper: U1 + U2 ⇒ limsup L·nymanDist² ≤ S; nbDist ≤ nymanDist.  value: Davenport input + RH.  summable: zero-count input.

## 2. U1 reduction (Mellin, closure, construction)  — [R](3i)(3ii)(3iii)(3v), [N]§6, [X]§7
- Closure.lean: G entire, |G s| ≤ C(1+|s|)^K e^{a max(0,−σ)}, a ≤ L − ½ ⇒ nymanDist(e^{−L}) ≤ ‖1 − ζG‖ + |G 1|·‖ζ‖.
  (mollify G·φ̂_δ^K, L² Paley–Wiener [input Rudin 19.3] ⇒ Laplace of h ∈ L²(0,b), b ≤ L; Bochner integral of
  ρ_{e^{-u}} − e^{-u}ρ_1 lies in closure(nymanSet); Mellin–Plancherel; dominated convergence δ → 0.)
- Mellin.lean: M(ρ_θ)(s) = −θ^s ζ(s)/s (0<σ<1); Mellin–Plancherel for the L¹∩L² functions used.
- Window.lean: Φ (spectrum [0,2/T₁]), bounds; Kernel B; Trivial.lean: R_tr; Contour.lean: J (divided-difference contour
  sums over good-height rectangles), tiling independence; Entire.lean: G := (H − J)/ζ − R_tr extends to an entire
  function; Growth.lean: (G2); LineBounds.lean: ‖J‖² ≤ e^{−cL²}, ‖ζR_tr‖², ‖1−Φ‖² ≤ C/L², |G 1| small.
- U1: nymanDist(e^{−L})² ≤ (1 + L^{−1/2})·Ustar + C·L^{−3/2}.

## 3. Dual form (Dual.lean) — [N]§1
- E_Z = span{u^k e^{iγu}}; Ustar Z ℓ := sup |p 0|²/N_ℓ p ≤ 1; p* representer; μ = ν_{p*};
  E_low(s) = ∫ e^{−u(s−½)} dμ interpolates Z with multiplicities, ‖E_low‖² = Ustar; |E_low z|² ≤ Ustar·N_ℓ(e_z).
- OU: ∫ e^{−ixt} w(t) dt = e^{−|x|/2}; ‖∫e^{−u(s−½)}dμ‖² = ∫∫e^{−|u−v|/2}dμdμ̄; representer identity.

## 4. U2 interpolation asymptotics (RH)
- ClaimF.lean (finite Z_lo, ℓ → ∞): |p(0)|² ≤ ((S_lo+ε)/ℓ)N(p). Christoffel: deg q < m ⇒ |q 0|² ≤ (m²/ℓ)∫_0^ℓ|q|²
  (explicit dual polynomial, Lagrange); polynomial Nikolskii/Markov constants by finite-dim norm equivalence.
- Groups.lean: groups at scale s₀ = L^{−3/4}; M_g ≤ M*(|γ_g|) from the RH short-interval count (input Littlewood).
- Multiplier.lean (replaces [N] Lemma A): r ∈ E_g (M nodes in [−w,w], W = max(w, KM²/ℓ)), Θ = e^{iθx}P(θ) near nodes,
  |Λ_r(Θ)| ≤ C^{M+deg P}·sup_{[−2W,2W]}|P|·sup_{[0,ℓ]}|r| (Fourier multiplier + Turán–Nazarov tails). Gives r′ bound,
  Newton coefficients c_k bound.
- TuranJets.lean: Turán–Nazarov (input, pure exponentials) ⇒ jets by forward-difference limits.
- LemmaB.lean, Newton.lean (Newton–Hermite expansion; Hermite–Genocchi input or proof), LemmaLambda.lean,
  LemmaD.lean (constant-C almost orthogonality suffices), LemmaE.lean, HighR.lean (R″ high part), U2.lean.

## 5. Literature inputs (LitInputs.lean; each a def … : Prop with citation + verbatim statement)
Burnol 2002 (math/0103058, Thm, TeX l.197–206); Friedland–Yomdin 1107.0039 Thm 1.1 (Nazarov); Goldston–Gonek math/0511092
(Littlewood, RH: N(t+h)−N(t) − (h/2π)log(t/2π) ≪ log t/log log t); Titchmarsh 9.7 (good heights); Titchmarsh 9.2/9.4
(zero counts) or zeta23's proofs; Rudin RCA 19.3 (Paley–Wiener); Davenport ch.12 (Σ Re 1/ρ); convexity |ζ| bounds in −1≤σ≤2.

## 6. Order of work
Defs → LitInputs → Main skeleton (sorry interfaces) → high-zero core (Multiplier, TuranJets, LemmaB, Newton, Λ, D, E) →
ClaimF → Dual → Contour/Entire/Growth/LineBounds → Mellin/Closure → U1/U2 → remove sorries → reviewers.
