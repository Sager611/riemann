# hmean working notes (appended in ≤60-line blocks)

## 0. Setup (03:25 UTC)
- S_B = Σ_{0<γ≤T} ζ′(ρ)²B(ρ)², B = Σ_{k≤y} μ(k)P(log(y/k)/log y)k^{-s}, y = T^θ. H: S_B ~ (TL³/2π)(u² − 1/12), u = ½ + θ∫P.
- Literature (scouts, .scratch/hmean/scout_k2.md, scout_lit.md): θ = 0 case Σζ′(ρ)² is CONJECTURED (RH) by Hughes–Pearce-Crump,
  arXiv 2509.07792 conj:secmom, leading term (1/2π)∫(L³/6 + …)dt; = our recipe at θ = 0 (1/6). No paper proves any k ≥ 2 holomorphic moment.
- Tools: Durkan–Hughes–Pearce-Crump 2601.18025: Σ_{T<γ≤2T} χ(ρ)X^ρ for all X ≥ 1 (three ranges, error O(T^{1/2}log²T) + edge terms); χ only, not χ².
  BHB 1302.5018: θ < ½ forced by Lemma 2, M_{ν,3} ≪ y^{1/3}T^{5/6+ε} (large moduli q ≤ y), via Heath-Brown's generalised Vaughan identity.

## 1. Exact identities at zeros (no RH)
- ζ(s) = χ(s)ζ(1−s) ⇒ ζ′(s) = χ′(s)ζ(1−s) − χ(s)ζ′(1−s); at any nontrivial zero ζ(1−ρ) = 0, so ζ′(ρ) = −χ(ρ)ζ′(1−ρ).
- Hence S_B = −Σ_γ χ(ρ)ζ′(ρ)ζ′(1−ρ)B(ρ)²; under RH = −Σ_γ χ(ρ)|ζ′(ρ)|²B(ρ)² (twisted second moment, twist χ(ρ)B(ρ)²).
- Also at zeros, ζ(ρ) = 0 in the AFE gives χ(ρ) ≈ −Σ_{n≤τ}n^{-ρ} (τ = γ/2π), hence ζ′(ρ) ≈ Ẽ(ρ) := Σ_{n≤τ}(ℓ − log n)n^{-ρ}, ℓ = log τ:
  S_B ≈ Σ_γ ζ′(ρ)·[Ẽ(ρ)B(ρ)²] = a CGG/BHB first moment Σζ′(ρ)A(ρ) with A of length T^{1+2θ} ≥ T, for EVERY θ ≥ 0.

## 2. Contour structure (left edge)
- Σ_γ G(ρ)w(γ) = (1/2πi)[∫_(c) − ∫_(1−c)](ζ′/ζ)GW; right edge and the χ′/χ part are O(T^{-A}) for G = ζ′²B² (no n = 1 coefficient).
- Left edge: ζ′(s)² = χ(s)²[Σ_a(log a − ℓ)a^{s−1}]², (ζ′/ζ)(1−s) = −ΣΛ(m)m^{s−1}; χ(s)²X^s has its stationary point at t₀ = 2π√X,
  value ∝ e(2√X). So S_B ≈ −(e^{iπ/4}/√2) Σ_h c_h h^{-3/4} Σ_{a,b,m} Λ(m)(log a−ℓ)(log b−ℓ)(abm)^{-1/4} w(t₀) e(2√(abm/h)), t₀ = 2π√(abm/h).
- Exact squares abm/h = k² alone give only O(T^{1/2}L³); the main term is a resonance (Voronoi-type) of the 3-fold convolution with e(2√·).
- For fixed prime m ≤ T^{1−δ}, Poisson in (a,b) resonates only at ν | m (ν ∈ {1, m}) and those terms are α-only / β-only in the shifted
  version (they die under ∂α∂β). The swap (pseudo-variance) term needs ν | a with ν ∤ m: it comes from the smooth part of Λ (composite m),
  i.e. from Poisson in m, available only where Λ(m) can be replaced by its mean.
- For m ≳ T the explicit formula for Σ_m Λ(m)e(2√(abm/h))g(m) returns a sum over zeros at height ≈ T of χ(ρ)²(ab/h)^{iγ}…: the method maps
  S_B to itself (checked: it reproduces Σ_γ [χ²-dual AFE part of ζ′(ρ)²]). No RH-only estimate breaks the loop.

## 3. Modulus localization (recipe, leading order) — kappa.py
- Ẽ_K(s) = Σ_{n≤K}(ℓ − log n)n^{-s} = (1/2πi)∫_(c) (K^z/z)(ℓ + ∂_β)ζ(s+β)|_{β=z} dz.  F_B(K) := Σ_γ ζ′(ρ)Ẽ_K(ρ)B(ρ)²; F_B(τ) ≈ S_B.
- Recipe: S_B(α,β) := Σ_γ ζ(ρ+α)ζ(ρ+β)B(ρ)² ≈ (T/2π)[ℓ + 2λ + (ζ′/ζ)(1+α) + (ζ′/ζ)(1+β) + D(α,β)],
  D = [β²g(α) − α²g(β)]/(αβ(β−α)), g(x) = e^{−xℓ}(1 − xλ)² = Σ g_j x^j, λ = μℓ, μ = θ∫P.
- ∂_α at α = 0 (the ±1/α cancel): Φ(β) = 1/β² + g₁/β + g₂ − g(β)/β² = −Σ_{j≥3} g_j β^{j−2} (entire).
- (ℓ + ∂)Φ = [−2/β³ + (ℓ − g₁)/β² + ℓg₁/β + ℓg₂] + e^{−βℓ}[2h/β³ − h′/β²], h = (1 − βλ)².
  Perron with K = τ^κ: polynomial part closes left (residue at 0); exponential part carries (K/τ)^z, closes right for K < τ: 0.
- Result: F_B(τ^κ) ~ (T/2π)ℓ³ f(κ), f(κ) = 1/6 + μ + μ² + μ(1 − κ)² + (1 − κ)³/3   (0 ≤ κ ≤ 1).
  Checks: f(1) = u² − 1/12 (H); f(0) = ½ + 2μ + μ² = recipe for ℓΣζ′(ρ)B(ρ)² (θ=0: Fujii's ½); the same recipe gives BHB's S₁ ~ (TL²/2π)(½ + θ∫P).
- Meaning: moduli n ∈ (τ^κ, τ] of the length-T twist carry μ(1−κ)² + (1−κ)³/3 (units ℓ² per zero); at κ = ½: μ/4 + 1/24, i.e. 25 % (θ = 0)
  to 22 % (θ = ½) of H's main term. The twist Ẽ_K B² has length τ^κ y²: BHB's first-moment range (length < T^{1/2}) is κ + 2θ < ½.
- Numerics, zeros 50001–10⁵ (L_rms = 9.11; logs/kappa_low.log; 59 s, memcap peak 70 MB): F(1) = exact S_B to 4 digits (ζ′(ρ) ≈ Ẽ_τ(ρ)).
  Tail F(κ) − F(1) vs prediction: θ = 0: κ = .5: .0436/.0417, .6: .0223/.0213, .7: .0096/.0090, .8: .0028/.0027;
  θ = ½ (μ = .257): κ = .5: .1103/.1060, .6: .0644/.0625, .7: .0318/.0322. Small κ deviates by O(1/L) (κ=0, θ=0: .562 vs .5).

## 4. Reviews and fixes (03:47–03:57)
- Reviewer 1 (main term): GAP;ROUTINE. Fixed: ζ′ ≈ Ẽ_τ only as a numerical reading (sharp AFE endpoint); exact use is ζ′ = Ẽ_K + χG_{τ/K}(1−ρ) + err
  (accurate at K ≍ √τ, checked to 4 digits); f(κ) for 0 < κ < 1.
- Reviewer 2 (obstruction): GAP;ROUTINE. Fixed: everything except A2 labelled heuristic; conclusion "these routes do not prove H"; ≥ 1/24 needs ∫P ≥ 0.

## 5. Three-piece split (kappa_low2/hi3): S_B = P₀ + 2X₁ + X₂ at K = M = √τ. Pure χ² piece X₂ = .160 (ℓ 9.1) → .313 (ℓ 24.5) at θ = ½,
  S_B = .530 → .485: the χ² part is large and growing; pieces cancel. A bound on the χ² part cannot suffice. Third height (zeros6 tail,
  ℓ ≈ 12.1, 2·10⁴ zeros) running: logs/kappa_z6.log.

## 6. Dead end checked (heuristic): lower bound Σ(Im Y)² via projection on short D(ρ): low-frequency coefficients of ẼB are
  ≈ (P′(1)/θ − P(1))Λ(N) (from Σ_{k|N}μ(k)P(1 − log k/log y) ≈ P′(1)Λ(N)/log y and Σ_{nk=N}log n μ(k) = Λ(N)); zero for BHB's P at θ = ½.
