# LiDens working notes (appended as we go; REPORT.md assembled at the end)

## Step 1: re-derivation of Li200k's per-zero inequality (falsification target)
- z_ρ = 1 − 1/ρ. ρ = 1/2 + δ + iγ, ρ* = 1 − ρ̄ = 1/2 − δ + iγ (same height γ > 0, same multiplicity by ξ(s) = ξ(1 − s), ξ(s̄) = conj ξ(s)).
- |z_ρ|² = (γ² + (1/2 − δ)²)/(γ² + (1/2 + δ)²) =: r², |z_ρ*| = 1/r; ℓ := −log r = ½ log(1 + 2δ/(γ² + (1/2 − δ)²)) ≤ δ/(γ² + (1/2 − δ)²) ≤ δ/γ².
- arg z_ρ = arg z_ρ* = φ = arctan((1/2 − δ)/γ) + arctan((1/2 + δ)/γ); θ(γ) = 2 arctan(1/(2γ)) = φ at δ = 0;
  h(x) = arctan(x/γ): h(c+δ) + h(c−δ) − 2h(c) = δ²h″(ξ), |h″(ξ)| ≤ 2ξ/γ³ ≤ 2/γ³ (ξ ≤ 1), so |φ − θ| ≤ 2δ²/γ³.
- Quartet ρ, ρ*, ρ̄, ρ̄*: Σ[1 − z^n] = 4 − 2(r^n + r^−n)cos nφ = Σ_{ρ, ρ*} [2 − 2cosh(nℓ)cos nφ].
  Per upper zero: 2 − 2cosh(nℓ)cos nφ ≥ g(nφ) − 2(cosh nℓ − 1) ≥ g(nθ) − 2n|φ − θ| − 4sinh²(nℓ/2) ≥ g(nθ(γ)) − ψ_n(γ, δ),
  ψ_n(t, a) := 4 sinh²(n a/(2t²)) + 4 n a²/t³ (Li200k's q_n with a(t) ≡ a), g(x) = 2 − 2cos x, |g′| ≤ 2.
- So λ_n ≥ M_n − Σ_{γ > T0, β ≠ 1/2} ψ_n(γ, |β − 1/2|) (on-line zeros: exact g-term, zero loss; γ ≤ T0 on the line by Platt–Trudgian).
- Monotonicity (analytic): ψ_n(t, a) increases in a ≥ 0 and in n, decreases in t > 0. For a' > a ≥ 0, f = ψ_n(t,a') − ψ_n(t,a):
  with s = n/(2t²), d/ds[sinh²(a's) − sinh²(as)] = a' sinh(2a's) − a sinh(2as) > 0 (x ↦ x sinh(2xs) increasing), so f decreases
  in t and increases in n; 4n(a'² − a²)/t³ likewise. Hence every Stieltjes step below has a nonnegative, decreasing integrand.

## Step 1b: what KLN counts (TeX 2101.12263)
- l.121: N(σ,T) = #{ρ: ζ(ρ) = 0, 0 < γ < T, σ < β < 1} (strict both); l.259–261: Littlewood's lemma (Titchmarsh 9.9.1) counts with multiplicity.
- l.193–217 (Thm 1.1 = Lemma 4.14 with C1, C2): for σ > 1/2 + d/log H0, T ≥ H0, H ∈ [1002, H0):
  N(σ,T) ≤ (T − H) log T/(2πd) · log(1 + C1 (log kT)^{2σ}(log T)^{4(1−σ)} T^{8(1−σ)/3}/(T − H)) + C2/(2πd) (log T)²  (KLN-NsigT1)
  ≤ A(σ)(log kT)^{2σ}(log T)^{5−4σ}T^{8(1−σ)/3} + B(σ)(log T)², A = C1/(2πd), B = C2/(2πd)  (KLN-NsigT2).
- l.1690: Table 1 gives A, B "for σ = σ0" (not for all σ ≥ σ0). Lemma 4.14 with a row's parameters holds for every σ > smin = 1/2 + d/log H0
  = 0.5117 with A(σ) = exp(c0 + c1σ), B(σ) = q0 − q1σ (Ingham4 cert.py class KP; c1 ≈ 5.0–5.1 > 0): we evaluate each row at σ itself.
- Rows at H0 = 3e12 (build_kln_H0('3e12'), logs/kln_rows.log): A(σ0) = 8.5122, 8.9580, 9.4220, 9.9044, B = 3.37737, 3.30437, 3.22157,
  3.13878 for σ0 = 0.85..0.88 — reproduces the brief. Row 0.60 not admissible at this H0.
- Count needed: zeros with |β − 1/2| > a_j and T0 < γ ≤ t = 2·N(σ_j, t⁺) (σ_j = 1/2 + a_j; β > σ_j plus reflections β < 1 − σ_j),
  ≤ 2·KLN(σ_j, t) (bound continuous in t, valid for t ≥ H0 = 3e12 < T0); N(σ_j, T0⁺) = 0 by RH to T0.

## Step 2: the modified bound (implemented in code/cert.py)
Slabs s_0 = s* < s_1 < … < s_J = B⁺ (exact dyadics; B⁺ ≥ θ₁₄₁ = 0.873345856106718576801184796… resp. 436666397/500000000), a_j = s_j − 1/2.
Every zero has 1 − θ ≤ β ≤ θ (half-plane + reflection), so |β − 1/2| ≤ a_J. Slab cap a(ρ) = a* = a_0 if |β − 1/2| ≤ a_0, a_{j+1} if a_j < |β − 1/2| ≤ a_{j+1}.
Per zero (T0 < γ ≤ t_K): ψ(γ,|β−½|) ≤ ψ(γ,a(ρ)) = ψ(γ,a*) + Σ_j f_j(γ)·1[|β−½| > a_j], f_j = ψ(·,a_{j+1}) − ψ(·,a_j) ≥ 0, ↓ in t.
Per zero (γ > t_K): ψ(γ,|β−½|) ≤ ψ(γ,a*) + ψ(γ,a_B). Hence
  Σ_{γ>T0} ψ_n(γ,|β−½|) ≤ Qall_n(a*) + Σ_j S_j(n) + R_n,
  Qall_n(a*) = Li200k Q_hi with a ≡ a* (verify.py l.126–145: 2ψ(T0)E(T0) + upper sum of ∫_{T0}^{100T0} ψ w + tail_bound),
  C_j(t) := #{T0 < γ ≤ t : |β−½| > a_j} = 2N(s_j, t⁺) (KLN count; β > s_j plus reflections β < 1 − s_j; none below T0),
  C_j(t) ≤ U_j(t) := min(2·KLN(s_j, t), N⁺(t)), N⁺(t) = Ñ(t) − Ñ(T0) + E(t) + E(T0) ≥ N(t) − N(T0⁻) (HSW),
  Σ_{T0<γ≤t_K, |β−½|>a_j} f_j(γ) = f_j(t_K)C_j(t_K) + ∫_{T0}^{t_K} C_j(t)(−f_j′(t))dt  (C_j(T0) = 0)
      ≤ S_j(n) := f_j(t_K)U_j(t_K) + Σ_k U_j(t_{k+1})(f_j(t_k) − f_j(t_{k+1}))   (C_j ↑, f_j ↓),
  R_n = 2ψ(t_K,a_B)E(t_K) + tail_bound(n, a_B, t_K) ≥ Σ_{γ>t_K} ψ(γ,a_B) (Li200k's Stieltjes step from t_K).
KLN(s,t) = (t−H) log t/(2πd)·log(1 + 2πd A(s)(log kt)^{2s}(log t)^{4(1−s)}t^{8(1−s)/3}/(t−H)) + B(s)(log t)² (Lemma 4.14 (bnd-N@H0), TeX l.1548–1566;
  H = H0 − 1, H0 = 3e12 ≤ T0, k = 1), for each slab the row (0.85–0.88) minimising it at t = T0; '--plain' uses (bnd-Nasymp) (the brief's form).
Grid: t_k − T0 = 1e8·1.05^k up to 0.02T0, then ×1.01 to 1.3T0, ×1.03 to t_K = 2T0 (exact integers; K = 172). Blocks: Q at the block's top
  (ψ_n ↑ in n per zero), gain M_lo(n1,n2) of Li200k; chain 1 = Li200k's own chain with a = a_B (89 … 106.255746 T0²), chain 2 continues.
Departures from the brief's literal step (2), both rigorous: (i) the count is min(2·KLN, total count N⁺) instead of 2·KLN; (ii) KLN's
  log-form (bnd-N@H0) instead of (bnd-Nasymp); '--plain --nomin' is the brief's literal bound. s* is free (slabs below 0.85 use row 0.85
  at σ < 0.85, valid since Lemma 4.14 holds for σ > 1/2 + d/log H0 = 0.5117).

