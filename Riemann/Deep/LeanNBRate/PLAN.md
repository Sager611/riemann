# LeanNBRate — PLAN (lemma DAG, Lean interfaces) — revised 22 Nov 2026 (session 5: DONE, Main.nb_rate proved)

Target: discoveries/nyman-beurling-continuous-rate.md, Claim (a)–(c). Proof sources: DbarNB/DERIVATION.md [N],
ExploreDbar/REPAIR.md [R], ExploreDbar/DERIVATION.md [X], DbarAudit/FULLPROOF.md [F]. Restructured for Lean; every
deviation from [N]/[R] is a re-proof, not an input. L = log(1/λ), ℓ = L − √L, s₀ = L^{−3/4}, T₂ ∈ [L⁶, L⁶+4] good height,
w(t) = 1/(2π(¼+t²)), Q_ℓ(p) = ∫_0^ℓ(|p′|² + ¼|p|²) (≤ N_ℓ(p), endpoint terms dropped: they are ≥ 0).

## 0. Statement (Defs.lean; Main.lean ✓ nb_rate : NBInputs → RH → claim)
MainClaim: RH → Summable ∧ Tendsto (log(1/λ)·nbDist²) (𝓝[>]0) (𝓝 S) ∧ ((∀ρ, mult ρ = 1) → S = 2+γ−log 4π) ∧ same for nymanDist.

## 1. Assembly ✓ (Assembly.mainClaim_of_upper, Final5.upper_bound, Main.nb_rate)
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

## 3. U1 ✓ — design v4 (session 4): explicit g by inverse Laplace, no PW, no growth bounds, no entire G
(✓ session 5: KJ1–5 identity theorem, ContourBd/KBd/RtrBd/GF1 bounds, LineE Plancherel identity, HConc1–3 concrete H, Final1–5)
Constants: L = log(1/λ) large, ℓ = L − √L (dual/U2), T₁ = L², n ≈ L, a = 2/(nT₁), T₂ ∈ [L⁶, L⁶+1] good height (T9.7,
σ ∈ [−1,2]); Ω = [−1,2]×[T₂,∞), Ω̄ mirror; line c = 3; κ(x) = e^{−Lx}; w = 1/(2π(¼+t²)).
✓ Dual1–5: W = L²(0,ℓ)²×ℂ² (Temb), Riesz representer p* (ExpSum, ZeroConfig), N(p*) ≤ (S+ε)/ℓ (U2), E_low(s) = N(e_s,p*)
  entire, Green form E_low(s) = ∫e^{−(s−½)u}conj(gDens) + conj(gA) + e^{−(s−½)ℓ}conj(gB), |E_low| ≤ √N(e_s)√N(p*).
- Window: M_a(x) = ∫_0^1 e^{−axv}dv, Φ(s) = M_a(s−½)^n: |Φ| ≤ 1, |1−Φ| ≤ |t|/T₁ on the line; |Φ| ≤ (2e^{a|σ−½|}/(a|t|))^n.
- H = (1−E_low)Φ (vanishes to order m_ρ at low zeros: E_low^{(k)}(ρ) = N((−u)^k e_ρ, p*) = δ_{k0}); F = H/ζ.
- g(u) = (2π)⁻¹∫F(3+iy)e^{u(3+iy)}dy (continuous); Ĝ(s) = ∫_0^L g e^{−us} = (2π)⁻¹∫F(z)(1−κ(s−z))/(s−z)dy (Fubini).
- Contour (σ_s > 3), Ψ₀ = F/(s−·), Ψ = Fκ(s−·)/(s−·): VI(Ψ₀,3) = 2πiF(s) (strip [3,c'], pole s, c'→∞);
  VI(Ψ,3) = VI(Ψ,2) = VI(Ψ,−1) + Rect[−1,2]×[−T₂,T₂] (= 0: removable at low zeros, 1) + UpperU − LowerU (PNT+
  DiffVertRect) and VI(Ψ,−1) = VI(Ψ,½−n) + 2πiΣ_{−2k>½−n}A_k, A_k = H(−2k)κ(s+2k)/((s+2k)ζ′(−2k)); n→∞ ⇒
  Ĝ = F − R_tr − J, R_tr = ΣA_k, J = (2πi)⁻¹(UpperU − LowerU)(Ψ).
- K(s) = (2πi)⁻¹(UpperU−LowerU)(w ↦ H(w)κ(s−w)dslope ζ w s/ζ(w)), analytic on ℂ∖{1} (parametric integral +
  Cauchy estimates); K = ζJ off Ω̄∪Ω̄c (Cauchy–Goursat on the half-strip); identity theorem on ℂ∖({1}∪{−2k})
  (countable complement, path-connected) ⇒ on the line 1 − ζĜ = (1−Φ) + ΦE_low + K + ζR_tr.
- Line bounds: ‖1−Φ‖²_w ≤ 4/T₁; ‖ΦE_low‖²_w ≤ ‖E_low‖²_w = N(p*) (ŵ = e^{−|x|/2} by Fourier inversion; Green:
  ∫e^{−|u−v|/2}dν = p*(v) by two FTCs); |K| ≤ tiny·(1+|ζ|) (split dslope = ζ(s)/(s−w) − ζ(w)/(s−w) for |s−w| ≥ 1/8,
  |ζ′| ≤ C|z| near the line from (2.1.5) + Cauchy); |ζR_tr| ≤ |ζ|·C√L e^{−(5/4)√L}; Nyman: c₀ = Ĝ(1), |Ĝ(1)| tiny.
- ζ/Γ facts: |Γ(z)| ≤ Γ(Re z); |Γ(½+it)|² = π/cosh πt; PL on [½,3/2] ⇒ |Γ(z)| ≤ √(2π)(|z|+1)e^{−π|Im z|/2};
  1/ζ(s) = 2(2π)^{−s}Γ(s)cos(πs/2)/ζ(1−s) (Re s < 0, s ∉ −ℕ) ⇒ polynomial |1/ζ| on lines/strips left; zeros in Re < 0
  are −2k, ζ′(−2k) = (−1)^k(2k)!ζ(2k+1)/(2(2π)^{2k}); |1/ζ| ≤ ζ(σ) (σ > 1, Möbius L-series).
- Rect.lean: port PNT+ (ResidueTheoremInRectangle, DiffVertRect_eq_UpperLowerUs); Strip lemma (finite simple poles,
  decay, T→∞), UpperU vanishing (holomorphic half-strip, decay), removable-point extension for Cauchy–Goursat.
- Closure: nymanDist_sq_le_laplace with g̃ = −g, c₀ = Ĝ(1) (constraint exact), L·D_N(e^{−L})² ≤ L(√N(p*) + small)².

## 4. Literature inputs (LitInputs.lean)
Burnol 2002; Nazarov/Friedland–Yomdin Turán; Goldston–Gonek 2007; Titchmarsh 9.2, 9.7, (2.1.5); Rudin 19.3 (PW);
Davenport ch.12. Addition planned: Rudin 9.13 (Plancherel for L¹∩L², Mathlib has only the abstract L² isometry).

