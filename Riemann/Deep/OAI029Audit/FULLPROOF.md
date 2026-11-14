# FULLPROOF — zero-exclusion chain of 029 Theorem 1.2 (map, with TeX line numbers)
Source: openai/math preprints/Primitive-roots-for-every-admissible-integer-base-October-4-2026/build/sections/*.tex (fetched 14 Nov 2026).
PDF numbering (main.tex order, theorem counter per section): §3=02-arithmetic (L3.1–3.4), §4=02a-coefficient-family (P4.1),
§5=03-theta (P5.1, C5.2, L5.3–5.5), §6=04-reflection (L6.1–6.4), §7=05-poisson (L7.1, P7.2, P7.3), §8=06-zero-free (L8.1, P8.2), App. A=10-whittaker-calculation.

## Statement (01-introduction.tex l.38–45)
Thm 1.2: F **cyclotomic** number field containing μ12; every finite-order Hecke character η of F: L_F(s,η) has no zeros in Re s > 1 − 10^−6.
(Wording note: our TREE line says "fields containing μ12"; the TeX requires F cyclotomic. 06-zero-free l.752–764 says "in the scope of Thm 1.2".)

## §8 (06-zero-free.tex) — the endgame, checked line by line
- L8.1 (l.37–145) conductor large sieve: Σ_ω |Σ b_n ω(n)|² ≪ (N+U²)(NU)^ε Σ|b|²; Gauss-sum expansion, orthogonality ≤ Σ_d |T_{d/q}|², then dual of L7.1 with V = C0·U. Checked: correct (standard).
- (8.3) conservative convexity |L(σ+it,ω)| ≪ U^{1/2+ε}(2+|t|)^{A_F}, .003 ≤ σ ≤ 4 (Tate). Standard.
- P8.2 (l.165–367) zero-density/detector: D_U mollifier length U², detector e^{−q/U^10}; shift to Re v = 1/2 − β:
  exponent 10(1/2−β)+3/2 ≤ −3.46 (β ≥ .996) ✓; dyad N ∈ (U²/2, U^11]: (N+U²)N^{1−2σ0} ≪ N^{.008} ≤ U^{.088} ✓;
  #detected ≪ U^{.088+ρ+ε} ≤ U^{1/2} (ρ < .1) ✓. Reciprocal bound: disks at 2+it, R2=1.0035 (left edge .9965 > .996 ✓),
  r2=1.003, r0=.999 (Re ≥ 1.001 ✓), r1=1.002 (reaches .998 ✓); θ = log(1.002/.999)/log(1.003/.999) ≈ 0.7505 < 1 ✓. Correct.
- Scales (from (8.x) contour-exponent, l.377–388): a = 9/10, b = 133/1000, h0 = 1/10 (solved from C*(ξ)=ξ−8/15 and the printed margins).
  E*(σξ,σw,σz) = a/2 + σξ − 1 + h0σz + b(σw − 1); C*(1) = 7/15.
- Principal row u=1 (l.420–470): I_1 = c0 f_η + R_w + R_z; R_w at (σ1,.95,ζ1): (1+h0)δ1 − .05b = −.00664989 ✓;
  R_z at (σ1,1,.16): δ1 + h0(.16−1/6) = −.000666566… ✓; δ1 = 10^−7. ξ stays at Re > 1: no zero input. ✓
- Nonprincipal rows (l.478–633): contour (σ1,.003,.4): E* − C*(1) = 10^−7 + 7/300 − .132601 = −.1092676… ✓;
  U ≤ Z^.09 (α=1): +.099 → −.0102676 ✓; exceptional E(U), Z^.09<U≤Z^.10001 (α=1/2): +.060006 → −.0492616 ✓;
  good rows moved to Re ξ = .998 via P8.2: −.0020001 + 7/300 − .132601 + .110011 = −377/300000 = −.00125667 ✓;
  U > Z^.10001: Re ξ = Re w = 4, Re z = R0 = 350000: 3.849 + .1R0 + .10001(1−R0) = .44901 ✓ (7/15 − .44901 = .017657 ✓).
  Total ≪ Z^{C*(1) − .0004} after losses < 10^−4 ✓ (slack .00126 − .0001 − .0004 = .00076).
- Completion (l.647–764): reflection exponent (2M0 − 1 − b)/2 = 933/2000 = 7/15 − 1/6000 (M0 = 1033/1000 ✓ exact);
  with P7.2 loss 10^−5: f_η(Z) ≪ Z^{7/15 − 1/20000} (1/6000 − 10^−5 > 1/20000 ✓); f_η = inverse Mellin of
  A(ξ) = e^{(ξ−5/6)²} H_η(ξ)/L^S(ξ,η) in the variable C*(ξ) = ξ − 8/15; T(s) converges for Re s > 1 − 1/20000 ✓;
  identity theorem ⇒ 1/L^S holomorphic there provided H_η holomorphic and zero-free on Re s > 1 − 1/20000 (P7.3 second neighbourhood).
  The width 10^−6 is a rounding-down of 1/20000 = 5·10^−5; the theorem's real input is the 1/6000 margin of the reflection bound.

## §7 (05-poisson.tex) — probe, Poisson, Euler product
- (7.4) probe I(X,Y,Z) l.92–102: ∫_{F_S^1/U_S^6} Y^{-1}Σ_s W1(q_s/Y)(q_sX)^{-1/2}Σ_m q_s^{-1/2}g_s(−m)W0(m/(t0ℓ_{q_sX}))B_m(Z)dt0.
- L7.1 l.118–192 (additive large sieve on R ⊂ F_S, (N+V²)): separation |d/s − d'/s' + r|_S ≥ V^{-2}, majorant Φ with δ = c·min(N^{-1/d_F},V^{-2/d_F}). Correct.
- P7.2 l.194–246: additive factor ℓ²-mass ≪ (XY)^{-1}(XY+Y²)/Y ≪ Y^{-1} (b < a); Cauchy–Schwarz with P4.1 ⇒ I ≪ Z^{(2M0−1−b)/2+ε}. Correct.
- Poisson l.248–389: finite transform (7.6) re-derived (m = v + Aj, no coprimality needed) ✓; multiplier (q_sX)^{-1/2}(X/q_A)q_s^{1/2} = X^{1/2}/q_A ✓ (L3.1: R self-dual, covolume 1);
  orbit unfolding over U_S^6 (free on h ≠ 0; quotient measure) gives F(s,A,h)V0(Xq_h/q_A) ✓; h = ur^6 bijection (L3.3) ✓;
  Mellin with t = ξ+z−1 gives (7.8) and exponent E* ✓. Normalization constants multiply every term of I alike: only c0 ≠ 0 matters.
- P7.3 l.407–766: local factor re-derived by hand (t = 6d, 6d+1, k = 0, 1; t ≡ 3, 4 vanish) ✓; closed form (7.14) ✓;
  region exponents −1.001, −1.908, −2.948 ✓; δ = 1/4000 ✓. Numerics (euler_check.py, euler_factor.py): C_{t,k}(J) brute force vs (7.12),
  136 cases at p = 13, 37, max rel diff 1.6e−12; closed Euler factor vs series 4.4e−16 (12 random points); H(E1=0)−1 formula 4e−10;
  planted errors (τ1→τ5; ϑ→−ϑ) fail (1.95; 7e−2).

## §4–§6 — reflected second moment (02a-coefficient-family.tex, 04-reflection.tex)
- B_m(Z) (4.2) l.41–49; P4.1 l.101–121: Σ_{h~P,d~D}|B_m|² ≪ Z^ε P^{1/2}(D + D²P/Z) ⇒ Σ|B_m|² ≪ Z^{2M0−1+ε}.
- L6.1 l.222–328: θ_υ chosen with S-Whittaker 1_U·(angular derivative Bessel); constant term killed by angular L^{-1} (6.6) ✓;
  frequency ν = γA forced by [νy] = 1, [y] ∈ B, E ∩ B = 1, L3.2 injectivity ✓; unfolding Σ_γ∫_{Ω_l} → ∫_{X_{q_A l}} (ℓ_T homomorphism,
  no Jacobian) ✓; L5.5 turns V(q_A T/Z) into P_G(q_A/Z) ✓. Exact identity with (4.2) via (6.1) and (6.2) ✓.
- L6.2 l.339–431: re-derived: torus factor χ(t)^4χ(v)^{4(e−2)}, spherical χ(t)^{2e}χ(v)^{−4e}W_e ⇒ (6.9) ✓; j=1,e=0: χ(t)^{8−1}, ω_p = Σχ³ψ/τ1 ✓.
- L6.3 l.625–910: ν0 = εa'c'n'^3/(df)²; local units at p|d, p|f, q|c' ✓; pairing χ_p(q)χ_q(p)^{−4} = χ_p(q)^3 (L3.2 reciprocity) ✓;
  reflected unfolding v = γν0/y (inversion preserves sheet measures) ✓; profile argument ZQ/(q_d²q_f²) ✓; Fubini majorant (6.24) ✓.
  Normalization self-check: reflected coefficients λ(c')/√q_c' with |λ| ≤ 1, no q_d power: matches a GL2 twist functional equation
  (dual length cond²/Z, root number of modulus 1).
- L6.4 l.941–1025 vs Goldmakher–Louvel arXiv 1112.1642v2 Def. 1 + Thm 1.1 (quoted in REPORT): family Ξ_c of F(√c)/F, c = q_S,
  order 2, reciprocity with C ≡ 1 (allowed generators), conductor of Ξ_{c1}Ξ̄_{c2} exactly c1c2 when cl2 agree ✓; rows/columns squarefree ✓.
- P4.1 proof l.1027–1228: tail (Φ ≪ x^{−N}) ✓; retained: Fourier separation (5.18), GL per τ, Minkowski over C-dyads, n', a' ✓;
  Σ_f ΠE_p ≪ q_h^{δ1} ✓; final P^{1/2}D ≤ Z^{M0} ≤ Z^{2M0−1}, P^{1/2}D²P/Z ≪ Z^{2M0−1} ✓.

## §3 (02-arithmetic.tex) — arithmetic inputs, checked line by line
- L3.1 l.56–110: R self-dual, covolume 1, Poisson on R and aR (Weil, Basic NT) ✓.
- L3.2 l.134–235: dim V_ℓ = 2s_f + 2r, dim U_S/U_S^ℓ = r + s_f (μ_ℓ ⊂ F, all places above 6 in S) ✓; injectivity via trivial class group of R ✓;
  H_S alternating since −1 = ζ12^6 ✓; E Lagrangian, isotropic complement B, CRT V = V_2 ⊕ V_3 ✓. Exact reciprocity (3.4) for allowed
  generators from the product formula and isotropy of B ✓; (3.6) χ_a(ε) = H_S(ε,a)^{-1} ✓.
- L3.3 l.285–385: Kummer characters, conductor ≤ C·N rad(u), u = 1 only principal, distinctness, bijection (u, r) ↦ [ur^6] ✓.
- L3.4 l.411–625: G2G4 = Q, G2² = J(ρ,ρ)G4, (Q^{-1/2}G2)³ = Q^{-1/2}J ✓; J(ρ,ρ) = −π_+ (J ≡ −1 mod 3, |J|² = Q, q | J) ✓;
  Davenport–Hasse lift ⇒ γ2(p)³ = −Λ(p), Λ = Λ_{Q(ω)}∘N angular ±1 ✓; γ3 via Weil index product formula, a character of B ✓;
  J(χ,χ) = χ̄(4)J(φ,χ) ⇒ G1G4 = χ̄(4)G2G3 ⇒ γ1γ2 = −ΛG ✓; CRT (3.9) ✓. Numerics: degree-one primes 13..97 (checks.log).

## §5 + App. A — theta input
- P5.1 = KP84 specialised (cover c = 2 on PGL2; exceptional datum |x³/y³|^{1/6}); C5.2 l.186–231 ✓ (unique c n³ decomposition, γ2 by CRT).
- (5.?) C_S = H_S² (03-theta l.233–251) from C_p = (x,y)_{6,p}² and Hilbert reciprocity ✓.
- Complex Bessel: Laplace/Gauss integral ⇒ ρK_{1/3}(cρ) ✓; bounds min(|y|^{1/3}, |y|^{-A}) ✓.
- Finite S-places (incl. above 2, 3): Jacquet exponents Q^{-1}, Q^{-2} per triple step (KP Thm I.2.9(e), no tameness hypothesis)
  ⇒ |W(y)| ≪ |y|^δ, δ < 1/3, at every finite place ✓ — this settles kusasyu36's wild-place concern (need only δ > 0).
- L5.3 finite cutoff ✓; L5.4 sheet profiles (weighted sup at v0 × weighted integrals elsewhere) ✓; L5.5 Mellin
  ∫t^{s+1}K_{1/3}(c√t)dt/t = 2^{2s+1}c^{-2s-2}Γ(s+5/6)Γ(s+7/6) ✓, K̂(−s) = e^{s²}/V̂(s) entire ✓ ⇒ exact P_G ✓.
- App. A l.183–495 re-derived: matrix identity wn(H) = n(H^{-1})diag(−H^{-1},H)n̄(H^{-1}) ✓; shell factor Q^{-2e/3} ✓;
  shell (A.?) Q^{l/3}I_{e,l}e_{2l−e} ✓; I_{e,l} table ✓; J_{3k} = Q^{-k}(e0 + Q^{-1/6}g4e2), J_{3k+1} = Q^{-k-1/2}g2(…), J_{3k+2} = 0 ✓;
  unit rule χ(u)^{2e} ✓. Orientation ι is a free choice (the two choices give complex-conjugate thetas); used consistently.
- Global Weyl lift (kusasyu36 gap 3): in KP's cocycle model the rational splitting is γ ↦ (γ_v, 1)_v (Π_v σ_v = 1 on rational pairs by
  Hilbert reciprocity; κ_v(γ) = 1 for almost all v), i.e. the product of the identity sections used locally; any residual discrepancy
  would be a fixed or per-row element of μ3, unimodular, invisible in Σ|B_m|².

## External inputs at statement level (sources: Numdam PMIHES_1984__59__35_0, PMIHES_1985__62__203_0; arXiv 1112.1642v2)
- KP84 Prop 0.1.3 (p.44): |n|_p = 1 ⇒ canonical lift = s on H ∩ K, W, N+ ∩ K ✓ (029 uses it only at exterior, tame places).
- KP84 Thm I.2.9 (p.72): V0(ω') unique irreducible sub of V(w0ω'), quotient of V(ω'); (e) Jacquet module = ind(w0ω'); (f) K'-fixed vector ✓.
- KP84 Cor I.3.6 (p.79): dim Wh = 1 iff r = n−1 & 2(c+1) ≡ 0 (mod n), or r = n; "this assumption [|n|_p = 1] is not needed" ✓ (r=2,n=3,c=2).
- KP84 Thm II.2.1 (p.118): global Θ ≅ ⊗ local V0(ω_v) ✓. Thm II.2.5 (p.122): (a) tame dims > 0 ⇒ global coefficient ≠ 0; (b) tame dims = 1 ⇒
  dim = 1 at every place ✓.
- KP85 erratum (p.203): pp.119–120 assertion false (Suzuki); c(η) = lim C_S(η)T(S), T(S) = Π(1+q_v^{-1})…; Thms II.2.2–II.2.3, p.130 modified.
  029 App. A states exactly this and never uses the withdrawn assertion (normalises each exterior functional at φ°_v) ✓.
- Goldmakher–Louvel Def. 1 + Thm 1.1 (quoted in REPORT): 029 L6.4's verification of (1) order 2, (2) reciprocity with C ≡ 1, (3) conductor of
  Ξ_{c1}Ξ̄_{c2} exactly c1c2 when cl2 agree, all match ✓; bound (M+N)(MN)^ε over squarefree ideals of I(c) ✓.
