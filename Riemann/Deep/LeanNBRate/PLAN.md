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
✓ LowPart.claimF_core (Claim F core), LemmaE.lemmaE (Lemma E′, no endpoint terms), Counting (cluster_bound, shells),
✓ ZeroConfig/GroupSums/U2grp/U2aux/U2asym/U2/U2props/U2main/U2final: U2final.u2 (RH, inputs TN, GG, T9.2):
  ∀ε>0 ∃ℓ₀ ∀ℓ≥ℓ₀ ∀T≤ℓ⁷ ∀P (ZeroConfig P T: nodes = distinct ordinates of zeros on the line, 1 ≤ mult ≤ m_ρ):
  ‖P.eval 0‖² ≤ (zeroSum+ε)/ℓ·Re oubil ℓ P P′ P P′ (= Q_ℓ(P)). Parameters: ℓ = q⁴, s₀ = q⁻³, η = q⁻², η′ = q⁻⁴,
  N = 40, multiplicity cutoff Ku = ⌊30ε₁ log q⌋+1, split height Hs in a zero-free gap.
- TODO Dual.lean: Hilbert space W = L²(0,ℓ)² × ℂ², T f = (f′, f/2, f(0)/√2, f(ℓ)/√2); f(0) = ⟨Tf, e⟩; p* := preimage of
  orthogonal projection of e on T(E_Z) ⇒ ⟨f,p*⟩_N = f(0) on E_Z. E_low(s) := ⟨e^{−(s−½)·}, p*⟩_N (entire, explicit),
  interpolates at Z with multiplicity (Z symmetric: needs mult(conj ρ) = mult ρ), ‖E_low‖²_w = N(p*) = p*(0) (OU kernel:
  ∫e^{−ixt}w = e^{−|x|/2}, reproducing identity), |E_low(s)|² ≤ N(p*)·N(e_s) (C-S).

## 3. U1 (to do) — design revised 22 Nov (session 3): single infinite contours, no tiles, no Stirling
Constants: L large, ℓ = L − √L, T₁ = L², T₂ ∈ [L⁶, L⁶+1] good height (T9.7), T₂′ ∈ [T₂−2, T₂−1] good height.
- Dual (Dual.lean): V = ℂ-span of u^k e^{iγu} (γ = ordinates of zeros with |γ| ≤ T₂, both signs, k < m_ρ), functions ℝ→ℂ.
  N(f,g) = ∫_0^ℓ(f′ḡ′ + ¼fḡ) + ½(f(0)ḡ(0) + f(ℓ)ḡ(ℓ)); definite on V (N(f)=0 ⇒ f=0 on [0,ℓ] ⇒ f=0, identity thm
  for real-analytic f). Type synonym EZ with InnerProductSpace.ofCore; Riesz (toDual) ⇒ p* with N(f,p*) = f(0) ∀f∈V.
  p* is an ExpSum with ZeroConfig ⇒ U2: p*(0) = N(p*) ≤ (S+ε)/ℓ. E_low(s) := N(e_s, p*), e_s(u) = e^{−(s−½)u}.
  Facts: entire; E_low^{(k)}(ρ) = δ_{k0} (k < m_ρ) ⇒ order(1−E_low, ρ) ≥ m_ρ; |E_low(w)|² ≤ N(e_w)·N(p*) (C–S);
  ‖E_low‖²_w = N(p*) via E_low(½+it) = ∫e^{−itu}dν̄*, ν* = (−p*″+¼p*)du + (½p*(0)−p*′(0))δ₀ + (½p*(ℓ)+p*′(ℓ))δ_ℓ,
  ∫e^{−ixt}w dt = e^{−|x|/2} (Fourier inversion of e^{−|x|/2}), Green identity ∫e^{−|u−v|/2}dν*(v) = p*(u).
- Window (Window.lean): φ = twin(2/T₁)/∫twin ≥ 0 on [0,2/T₁], Φ(s) = ∫φ(u)e^{−(s−½)u}du. |Φ(½+it)| ≤ 1,
  |1−Φ(½+it)| ≤ 2|t|/T₁, Gevrey decay |Φ(σ+it)| ≤ C e^{−c√(|t|/T₁)} for −2 ≤ σ ≤ 3 (IBP n times, n ≈ √(|t|/T₁)).
- H := (1 − E_low)Φ. Kernel κ(x) = e^{−(L−1)x}. D(s,w) := dslope ζ w s (= ζ[s,w]).
  Ω = [−½, 3/2]×[T₂,∞), Ω̄ its conjugate. K(s) = (2πi)⁻¹∫_{∂Ω}H(w)κ(s−w)D(s,w)/ζ(w)dw (+ Ω̄), analytic on ℂ∖{1}.
  J(s) = (2πi)⁻¹∫_{∂Ω}H(w)κ(s−w)/((s−w)ζ(w))dw, analytic off ∂Ω. Cauchy (rectangles [−½,3/2]×[T₂,R], R→∞, H decays;
  Goursat (Mathlib) + ∮1/(w−s) = 2πi (port PNT+ ResidueCalcOnRectangles)) ⇒ K = ζJ + H·1_Ω.
  Trivial zeros: R_tr(s) = Σ_k r_k e^{−(L−1)(s+2k)}/(s+2k), r_k = (H−K)(−2k)/ζ′(−2k); ζ′(−2k) = Λ(2k+1)π^{−k}(−1)^k k!/2
  (ζ = Λ/Γℝ, (1/Γ)′(−k) = (−1)^k k!), |r_k| super-exponentially small.
  G₀ := (H − K)/ζ − R_tr; removable singularities at zeros (inside Ω: = −J − R_tr; low zeros: order(H) ≥ m_ρ;
  −2k: simple pole cancelled; s = 1: 1/ζ → 0) ⇒ G entire (limUnder update / toMeromorphicNFOn).
- Line bounds: 1 − ζG = (1−Φ) + ΦE_low + K + ζR_tr on σ = ½; |K(s)| ≤ (|ζ(s)|+1)·tiny (|H| ≤ e^{−cL²} on ∂Ω,
  |1/ζ| ≤ ζ(3/2) on σ = 3/2, ≤ 4πζ(3/2) on σ = −½ (exact: ζ = Λ/Γℝ, Λ(s)=Λ(1−s), |Γ(¾+iy)| = |(−¼+iy)Γ(−¼+iy)|),
  ≤ CT₂^A on the bottom side; D(s,w) near the crossing via ζ′ bound from (2.1.5): |ζ(z)| ≤ |z|(1/Re z + 1/(1−Re z))).
  ∫|ζ(½+it)|²w < ∞ from the Mellin identity for ρ(1/x). ⇒ ‖1−ζG‖²_w ≤ (1+L^{−1/2})N(p*) + O(L^{−3/2}).
- Growth (for PW; constants may depend on L): σ ≥ 7/4: H/ζ − J, |1/ζ| ≤ ζ(7/4); inside Ω off ∂Ω: −J; compact
  |t| ≤ T₂+2, bounded σ: continuity; near ∂Ω sides: PL in vertical strips [−3/2, ½] (boundary: σ = −3/2 good line,
  σ = ½ inside Ω / compact) and right side analog; left half-plane: PL in strips between good lines σ = ½−n
  (|1/ζ(½−n+it)| = π^n/Π_{j<n}|s/2+j|·|1/ζ(½+n−it)|, elementary), crude a priori bounds (|Γ(z)| ≤ Γ(Re z), reflection,
  duplication), max principle on discs around −2k. Result |G(s)| ≤ C(1+|s|)^N e^{(L−1)max(0,−σ)}.
- PW (Rudin 19.3) on f(z) = e^{−i(a/2)z}·(G·M^{N+2})(c − iz), M(s) = (1−e^{−δs})/(δs); support [−ε, a+ε]; shift by e^{−εs}
  ⇒ G̃ = ĝ, g ∈ L² supported in [0, L]; DCT as δ, ε → 0.
- Closure (Mellin.lean): ‖χ − Σc_kρ_{θ_k}‖² = (2π)⁻¹∫|1 + ζ(s)Σc_kθ_k^s|²/|s|² (Plancherel L¹∩L² input + mellin_eq_fourier
  + (2.1.5) ⇒ M[ρ_θ](s) = −θ^sζ(s)/s); continuous g: Riemann sums + DCT; L² g: L¹-density; constraint via G̃ − G̃(1)·1
  (θ = 1 term), |G(1)| tiny (= −J(1) − R_tr(1)).

## 4. Literature inputs (LitInputs.lean)
Burnol 2002; Nazarov/Friedland–Yomdin Turán; Goldston–Gonek 2007; Titchmarsh 9.2, 9.7, (2.1.5); Rudin 19.3 (PW);
Davenport ch.12. Addition planned: Rudin 9.13 (Plancherel for L¹∩L², Mathlib has only the abstract L² isometry).

