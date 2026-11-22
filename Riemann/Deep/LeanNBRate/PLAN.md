# LeanNBRate — PLAN (lemma DAG, Lean interfaces) — revised 22 Nov 2026 (session 2)

Target: discoveries/nyman-beurling-continuous-rate.md, Claim (a)–(c). Proof sources: DbarNB/DERIVATION.md [N],
ExploreDbar/REPAIR.md [R], ExploreDbar/DERIVATION.md [X], DbarAudit/FULLPROOF.md [F]. Restructured for Lean; every
deviation from [N]/[R] is a re-proof, not an input. L = log(1/λ), ℓ = L − √L, s₀ = L^{−3/4}, T₂ ∈ [L⁶, L⁶+4] good height,
w(t) = 1/(2π(¼+t²)), Q_ℓ(p) = ∫_0^ℓ(|p′|² + ¼|p|²) (≤ N_ℓ(p), endpoint terms dropped: they are ≥ 0).

## 0. Statement (Defs.lean; Main.lean to do)
MainClaim: RH → Summable ∧ Tendsto (log(1/λ)·nbDist²) (𝓝[>]0) (𝓝 S) ∧ ((∀ρ, mult ρ = 1) → S = 2+γ−log 4π) ∧ same for nymanDist.

## 1. Assembly (Main.lean, to do)
- lower: Burnol input (complex L², Ioc) ⇒ nbDist ≥ nbDistC (real span embeds isometrically) ⇒ liminf ≥ S; nymanDist ≥ nbDist.
- upper: limsup L·nymanDist² ≤ S from U1 + U2 (real part of a complex Nyman approximant is a real one).
- value: Davenport + RH (Re 1/ρ = ½/|ρ|²). summable: Titchmarsh 9.2 count.

## 2. U2 (DONE pieces marked ✓)
✓ ExpSum (jets), turan_expSum, Nikolskii, turan_middle, sup_sq_le, Kernel/Multiplier/MLBound, Bernstein, Markov (w=0).
✓ LemmaB (one group, γ² gain), LemmaLambda (frequency leakage, polynomial decay order N), tail_bound,
✓ LemmaD crude-count form: (1 − ((1+η⁻¹)G+1)2πδ)Σ_g∫|Q_g|² ≤ (1+η)∫|ΣQ_g|² (G = #groups ≤ #zeros ≤ poly(L)).
✓ Groups.lean (components of the s₀/2-nbhd: separation, width ≤ #g·s₀, grp_mult_le from a short-interval count),
✓ Config.lean (ExpSum.part T c, nodeSet, groups, gC/gW, groups_sep, part_eval_groups),
✓ HighR.highR: (1 − Cδ)|p_hi(0)|² ≤ (Σ_g 16(2κ)^{4(M_g−1)}/c_g²)/ℓ·(1+η')∫|p_hi′ψ|².
✓ Legendre.christoffel: deg P < m ⇒ |P(0)|² ≤ m²∫_0^1|P|²;  LowPart.christoffel_expSum (all nodes 0).
✓ LowPart: ibp_osc_bound, oubil (B(f,g) = ∫_0^ℓ f′ḡ′ + ¼fḡ), cross_bound (generic, sup bounds in, 1/|γ₁−γ₂| out).
- TODO LowPart: diagonal Q(e^{iγu}r) ≥ (γ²+¼ − C/ℓ)∫|r|²; Claim F core: Q(p_lo) ≥ (1 − C/ℓ)Σ_j(γ_j²+¼)∫|r_j|²,
  |p_lo(0)|² ≤ (S_lo/ℓ)·Q(p_lo)/(1 − C/ℓ) (norm_sum_sq_le_of_le + christoffel_expSum).
- TODO LemmaE′ (no endpoint terms): |B(p_lo,p_hi)| ≤ θ_L (Q(p_lo)·Σ_g∫|Q_g|²)^{1/2}, via cross_bound per (j,g)
  with group sups (Lemma B by-product, sup_sq_le, Bernstein η = s₀), C-S over g with Σ_g e^{βM_g}/Ω_g² ≤ C(H).
- TODO U2.lean: Q(p) ≥ Q(p_lo) + ∫|p_hi′ψ|² − 2|B| ⇒ |p(0)|² ≤ ((S_lo + ε + S_hi)/ℓ)Q(p)/(1−θ). Inputs: GG ⇒ M_g ≤
  A log t/log log t (grp_mult_le); Titchmarsh 9.2 ⇒ #zeros ≤ C T log T, Σ_ρ|γ|^{−3/2} < ∞ (tails → 0).
- TODO Dual.lean: Hilbert space W = L²(0,ℓ)² × ℂ², T f = (f′, f/2, f(0)/√2, f(ℓ)/√2); f(0) = ⟨Tf, e⟩; p* := preimage of
  orthogonal projection of e on T(E_Z) ⇒ ⟨f,p*⟩_N = f(0) on E_Z. E_low(s) := ⟨e^{−(s−½)·}, p*⟩_N (entire, explicit),
  interpolates at Z with multiplicity (Z symmetric: needs mult(conj ρ) = mult ρ), ‖E_low‖²_w = N(p*) = p*(0) (OU kernel:
  ∫e^{−ixt}w = e^{−|x|/2}, reproducing identity), |E_low(s)|² ≤ N(p*)·N(e_s) (C-S).

## 3. U1 (to do) — design decided
- Φ: Laplace transform of a Gevrey bump on [0, 2/T₁], ∫φ = 1; |Φ| ≤ 1 on the line, ‖1−Φ‖²_w ≪ L⁻², |Φ| ≤ (CN²/(δ₁|t|))^N.
- Tiles Q_n = [−½,3/2]×[t_n,t_{n+1}] (good heights from Titchmarsh 9.7, gaps ≥ 2 by choosing t_n ∈ [3n,3n+1]) + mirrors.
  K_n(s) := (2πi)⁻¹∮_{∂Q_n} H(w)B(s−w)ζ[s,w]/ζ(w)dw, H = ΨΦ, B(x) = ((1−e^{−ax})/(ax))², a = (L−1)/2.
  Cauchy formula on rectangles (port PNT+/zeta23 ResidueTheoremAtOrigin' for ∮1/(w−s) = 2πi) ⇒ near any high zero
  H − ΣK = −ζ·ΣJ ⇒ G₀ := (H − ΣK)/ζ − R_tr has nonnegative meromorphic order everywhere; G := toMeromorphicNFOn G₀.
- R_tr(s) = Σ_k r_k e^{−(L−1)(s+2k)}/(s+2k), r_k = H(−2k)/ζ′(−2k), ζ′(−2k) from Mathlib's functional equation.
- Growth for PW: maximum principle on rectangles [−¾,7/4]×[t_n,t_{n+1}] (K-form on good heights, J-form on σ = −¾, 7/4);
  compactness for |t| ≤ T₂+4; σ > 7/4 Dirichlet bound; σ < −¾ functional equation + Gamma bound (input).
- PW: prove by contour shift (Cauchy on rectangles + Fourier inversion) for mollified G (decay (1+|s|)^{−2}).
- Closure: Riemann sums Σh(u_k)Δ(ρ_{e^{−u_k}} − e^{−u_k}ρ₁) ∈ nymanSet; dominated convergence on the Mellin side
  (needs Mellin–Plancherel for these functions: L¹∩L² Plancherel input if Mathlib lacks it) ⇒ nymanDist ≤ ‖1 − ζA‖_w.

## 4. Literature inputs (LitInputs.lean)
Burnol 2002; Nazarov/Friedland–Yomdin Turán; Goldston–Gonek 2007; Titchmarsh 9.2, 9.7, (2.1.5); Rudin 19.3 (may be
dropped if PW is proved); Davenport ch.12. Likely additions for U1: convexity |ζ|,|ζ′| ≪ (2+|t|)^B in −1 ≤ σ ≤ 2;
|1/ζ| bounded on σ = −¾ (|t| ≥ 1) via functional equation/Stirling; L¹∩L² Plancherel (Rudin 9.13) if needed.

