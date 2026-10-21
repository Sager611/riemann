# CMQuartic: class group characters of CM fields — relative large sieve, zero density, primes in every class for almost all fields

Date 2026-10-21. Brief: .scratch/orca/cmquartic-brief.md. Paper proofs, not refereed, no Lean. References: [C] = [Combos NOTE](../Combos/NOTE.md)
(Lemma 3, Lemma 5, (F4)–(F6)), [D2] = [Density2 NOTE](../Density2/NOTE.md) (Proposition 6, Proposition 7, Lemma 4), [R] = [RCApps NOTE](../RCApps/NOTE.md) (§2),
[R2] = [RCApps2 NOTE](../RCApps2/NOTE.md) (Proposition 1, Theorems 2′ and 4). Scratch (.scratch/cmquartic/): cmcount.gp + run2.gp (class counts, Lemmas 2–3),
reldet.gp (Lemma 1), leastp.gp (least primes per class), params_cm.py → params_cm.out (Theorem B parameters), axq*.py → axq*.out (arXiv), gxl.py → gxl.out,
src/ (sources of Thorner–Zaman 1604.01750 and Habegger–Pazuki 1506.02485). Unconditional unless marked.

## Verdicts

| Question | Answer | Status |
|---|---|---|
| usable large sieve over ideals of a quartic (any) CM field K, F = K⁺ fixed? | **yes**: weights on the O_F-content of 𝔞, characters trivial on ideals from F (≤ h_F cosets cover Ĉl_K); constant C_F ε⁻¹log^{2n−1}(3d_K)(N + √d_K), essentially sharp (Lemma 3). Key step: ββ′ᶜ − βᶜβ′ ∈ N_{K/F}(𝔠)𝔇_{K/F} (Lemma 1), so a class holds ≤ 2^{n−1}(1 + (4π)ⁿe^{nr_F}d_F N/√d_K) O_F-primitive ideals of norm ≤ N (Lemma 2) | new; resolves the Density2 lead. The Z-primitive version of [D2] fails for CM fields, by d_K^{1/4−o(1)} (Remark 2.2) |
| zero density, each K | Theorem A: Σ_{χ∈Ĉl_K} N(σ;t,χ) ≤ C(F,ε)(d_K Tⁿ)^{2(1−σ)/(3−2σ)+ε} for 0.6 ≤ σ ≤ 1 | new statement; Pascadi–Thorner (Thm 1.3) gives d_K^{3(1−σ)/(3−2σ)} |
| zero density, family | Proposition C: over the CM fields K/F of the δ₀-rough family with N𝔡_{K/F} ∈ (X/2, X], the non-genus characters have Σ N(σ,T) ≤ C(XT^{(3n+2)/5})^{10(1−σ)/(3−2σ)+ε} | new (Humphries–Thorner over F, fed with Theorem A's second moment) |
| application | Theorem B: for every A > 3/2, all but O(X^{1−c}) fields of that family have ≥ H/(16h_K log 2x) primes 𝔓 of degree one in every ideal class with N𝔓 ∈ (x, x+H], for all x ≥ D^A, H ≥ xD^{−κ} (D = N𝔡_{K/F} ≍ d_K). Corollary B1: some β ∈ O_K has prime norm ≤ 2d_K^A; every β ∈ O_K∖F has N(β) ≥ d_K/(4ⁿd_F²) | new. Per field before: d_K^{694} (Thorner–Zaman 2017), GRH d_K^{1+ε} |
| where it stops | all CM K/F (small genus conductors need a Q-aspect log-free density with Deuring–Heilbronn over F), F varying (constants), A < 3/2 (same barrier as [R2] §4), fields with no quadratic subfield ([D2] Prop 7) | §7 |

## 0. Notation and inputs

- F totally real, n := [F:Q], fixed; ε_F > 1 a fundamental unit if n = 2; h_F, h_F⁺ = |Cl_F⁺| (narrow class group). All constants C_F, C(F,…)
  are effective and may depend on F.
- K/F **CM**: [K:F] = 2, K totally imaginary; c = the nontrivial element of Gal(K/F). For every embedding σ: K → C, σ∘c = σ̄. For each real
  place v of F fix σ_v: K → C above v; then N_{K/Q}β = Π_v|σ_vβ|², and θ_v(β) := arg σ_vβ ∈ R/2πZ.
- 𝔡 = 𝔡_{K/F}, 𝔇 = 𝔇_{K/F} (N_{K/F}𝔇 = 𝔡), D := N_{F/Q}𝔡, d_K = d_F²D (tower formula). η = η_{K/F}: the quadratic Hecke character of F of K.
  For a fractional ideal 𝔠 of K, 𝔠𝔠ᶜ = N_{K/F}(𝔠)O_K; for an ideal 𝔤 of F, N_{K/Q}(𝔤O_K) = N_{F/Q}(𝔤)². N := N_{K/Q} on K unless marked.
- Cl_K, h_K; j: Cl_F → Cl_K extension of ideals; **𝒵_K** := {χ ∈ Ĉl_K : χ∘j = 1} = the dual of Cl_K/j(Cl_F) ("F-trivial" characters). Ĉl_K is
  the disjoint union of |j(Cl_F)| ≤ h_F cosets ψ₀𝒵_K.
- **Content.** For an integral ideal 𝔞 of K, 𝔤(𝔞) := the largest integral ideal 𝔤 of F with 𝔤O_K | 𝔞 (it exists: 𝔤O_K | 𝔞 iff
  e(𝔓/𝔭)v_𝔭(𝔤) ≤ v_𝔓(𝔞) for all 𝔓 | 𝔭). 𝔞 is **O_F-primitive** if 𝔤(𝔞) = O_F. Every 𝔞 is uniquely 𝔞 = 𝔤(𝔞)O_K·𝔞₀ with 𝔞₀ O_F-primitive, and
  N_{F/Q}(𝔤(𝔞))² ≤ N𝔞. Weighted norm: ‖c‖²_ε := Σ_𝔞 N_{F/Q}(𝔤(𝔞))^{1+2ε}|c_𝔞|².
- N(σ;t,χ) := #{ρ = β+iγ : L(ρ,χ) = 0, β ≥ σ, |γ−t| ≤ 1}, N(σ,T;χ) := #{ρ : β ≥ σ, |γ| ≤ T}, with multiplicity; T := |t|+3.
  e(σ) := 2(1−σ)/(3−2σ).

Facts.
- (a) For χ ∈ 𝒵_K and integral 𝔞: χ(𝔞) = χ(𝔞₀), and χ̄(𝔞) = χ(𝔞ᶜ) because 𝔞𝔞ᶜ = N_{K/F}(𝔞)O_K.
- (b) Every χ ∈ Ĉl_K is primitive with conductor O_K; Λ(s,χ) = d_K^{s/2}Γ_C(s)ⁿL(s,χ) satisfies Λ(s,χ) = W(χ)Λ(1−s,χ̄), |W(χ)| = 1 (Hecke).
  All χ have the same gamma factor and conductor d_K; the analytic conductor at height t is ≍_n d_K(|t|+3)^{2n}.
- (c) h_K ≤ C(n)√d_K log^{2n−1}(3d_K): h_K = w_K√d_K κ_K/((2π)ⁿR_K) (class number formula, r₂ = n), κ_K = res_{s=1}ζ_K ≤
  (e log d_K/(2(2n−1)))^{2n−1} (Louboutin 2000), R_K ≥ 0.2 (Friedman 1989), w_K ≤ C(n).

Inputs (cited).
- (F4) [IK Prop 5.7]: a primitive Hecke L-function of K (of F) has ≤ c₀log(d_K(|t|+3)^{2n}) (≤ c₀log(d_FN𝔣(|t|+3)ⁿ)) zeros with |γ−t| ≤ 1.
- (F5) [IK Thm 5.3, Prop 5.4]: the AFE for χ ∈ Ĉl_K∖{1} at s = ½+iu, with A = d_K and V_s(y) = (2πi)⁻¹∫_{(2)}y^{−w}e^{w²}(Γ_C(s+w)/Γ_C(s))ⁿdw/w;
  Stirling gives |(y∂_y)^rV_s(y)| ≤ C(n,r,B)min(1, (y/(|u|+3)ⁿ)^{−B}).
- (I6) [IK Thm 5.11; [R] (I6)]: for w ∈ C_c^∞([x,x+H]), x ≥ 2, W(s) := ∫w(u)u^{s−1}du and a primitive Hecke character χ of K or of F:
  Σ_𝔞 Λ(𝔞)χ(𝔞)w(N𝔞) = δ_{χ=1}W(1) − Σ_ρ W(ρ) + R_χ, |R_χ| ≤ C(n)(H/x + x^{−1/2})log(A(χ)x).
- (HT) Humphries–Thorner, Compos. Math. 160 (2024), arXiv 2103.05634, over the number field F (their paper works over any number field; read
  in .scratch/rcapps/src/2103.05634): let m ∈ {1,2}, 𝒮 a finite set of cuspidal automorphic representations of GL_m(A_F) (unitary central
  character trivial on R_{>0}) satisfying GRC, Q ≥ max_{π∈𝒮}C(π) (C(π) = d_F^m N𝔮_π Π_{v|∞}Π_j(3+|μ_{j,π}(v)|)), T ≥ 3, ε > 0 and
  V := Q^{m+ε}T^{m²n/2+1+ε}|𝒮|. For 0 ≤ σ ≤ 1 and every Y ≥ 3 with log Y ≍ log V,
  Σ_{π∈𝒮}N_π(σ,T) ≪_{ε,F} V^ε[Y^{2(1−σ)} + Y^{1/2−σ}(V·I₂)^{1/2} + Y^{1−σ}],  I₂ := Σ_{π∈𝒮}∫_{−T}^{T}|L^{ur}(½+iv,π)|²dv,
  and (convexity) I₂ ≪_ε V^εQ^{m/2}T^{mn/2+1}|𝒮|. This is HT's inequality (bound_int) in §5 with both μ-integrals bounded by their Cor 4.4
  (cor:MVT) and θ_m = 0 (their Thm 4.1 uses the best exponent towards GRC for π ∈ 𝒮): [R] (I3) over F.
- (TZ) Thorner–Zaman, Algebra Number Theory 11 (2017), arXiv 1604.01750, Theorem 3.2 (arXiv numbering; source read): for a number field M and a
  congruence class group H (narrow ray classes), Σ_{χ mod H}N(σ,T,χ) ≪ (e^{O(n_M)}D_M²Q_HT^{n_M+2})^{81(1−σ)} for ½ ≤ σ < 1 and
  T ≥ max(n_M^{5/6}D_M^{−4/(3n_M)}Q_H^{−4/(9n_M)}, 1), absolute constants. Used once: M := H_F⁺ (narrow Hilbert class field of F), H = all ideals,
  i.e. for the single fixed function ζ_M = Π_{ψ∈Ĉl_F⁺}L(s,ψ): N_M(σ,T) ≤ C_M T^{B_M(1−σ)} (T ≥ 3), B_M := 81([M:Q]+2).
- (ZFR) For the fixed field M there are effective c_M > 0, β_M < 1 such that every zero of ζ_M has β ≤ β_M if |γ| ≤ 3 and
  β ≤ 1 − c_M/log(|γ|+3) if |γ| ≥ 3 (classical: Landau; e.g. Lagarias–Montgomery–Odlyzko 1979; Stark 1974 for an effective β_M).

## 1. Relative lattice lemmas

**Lemma 1 (relative determinant).** Let 𝔠 be a fractional ideal of K and β, β′ ∈ 𝔠 with β′/β ∉ F. Then Δ := ββ′ᶜ − βᶜβ′ ≠ 0, Δᶜ = −Δ,
Δ ∈ N_{K/F}(𝔠)𝔇, and

  Π_v 2|σ_vβ||σ_vβ′||sin(θ_v(β) − θ_v(β′))| = |N_{K/Q}Δ|^{1/2} ≥ N𝔠·√d_K/d_F.

*Proof.* Δ = 0 ⟺ (β′/β)ᶜ = β′/β ⟺ β′/β ∈ F; Δᶜ = −Δ is clear. Put x := ββ′ᶜ ∈ 𝔠𝔠ᶜ = N_{K/F}(𝔠)O_K, so Δ = x − xᶜ. *Claim:* if 𝔤 is a
fractional ideal of F and x ∈ 𝔤O_K, then x − xᶜ ∈ 𝔤𝔇. It suffices to check this at every prime 𝔭 of F. O_{K,𝔭} := O_K ⊗ O_{F,𝔭} is free of rank 2
over the DVR O_{F,𝔭} and 1 is primitive in it (π_𝔭⁻¹ ∉ O_{K,𝔭}), so O_{K,𝔭} = O_{F,𝔭}[y] for some y; then 𝔇O_{K,𝔭} = (f′(y)) = (y − yᶜ) for the
minimal polynomial f of y (Neukirch, ANT III (2.4)). With 𝔤O_{F,𝔭} = gO_{F,𝔭} write x = g(a + by), a, b ∈ O_{F,𝔭}; then x − xᶜ = gb(y − yᶜ) ∈ g𝔇O_{K,𝔭}.
So (Δ) ⊆ N_{K/F}(𝔠)𝔇 and |N_{K/Q}Δ| ≥ N_{K/Q}(N_{K/F}(𝔠)O_K)·N_{K/Q}(𝔇) = (N𝔠)²·D = (N𝔠)²d_K/d_F². Finally σ_v∘c = σ̄_v gives
σ_vΔ = z − z̄ with z = σ_vβ·conj(σ_vβ′), so |σ_vΔ| = 2|Im z| = 2|σ_vβ||σ_vβ′||sin(θ_v(β) − θ_v(β′))|, and |N_{K/Q}Δ| = Π_v|σ_vΔ|². □

For K imaginary quadratic (F = Q) this is "det(β, β′) ∈ covol·Z", the input of [D2] Lemma 1. The point is that the CM involution gives a
determinant over O_F, and the different makes its norm as large as √d_K, not just N𝔠.

**Lemma 2 (O_F-primitive ideals in a class).** There is r_F ≥ 0 depending only on F (for n = 2: e^{r_F} = ε_F) such that for every CM field K
over F, every C ∈ Cl_K and every N ≥ 1,

  π_C(N) := #{𝔞 ∈ C integral, O_F-primitive, N𝔞 ≤ N} ≤ 2^{n−1}(1 + (4π)ⁿe^{nr_F}d_F·N/√d_K).

For n = 2: π_C(N) ≤ 2 + 2(4π)²ε_F²d_F N/√d_K (= 2 + 4134N/√d_K for F = Q(√5)).

*Proof.* Fix an integral 𝔠 ∈ C⁻¹. The integral 𝔞 ∈ C are the (β)𝔠⁻¹ with 0 ≠ β ∈ 𝔠, β unique up to O_K^×, N𝔞 = Nβ/N𝔠, and 𝔞 ⊆ 𝔭O_K ⟺ β ∈ 𝔭𝔠.
(i) *Proportional primitive generators are unit multiples.* Let 𝔞 = (β)𝔠⁻¹, 𝔞′ = (β′)𝔠⁻¹ be O_F-primitive with λ := β′/β ∈ F. Write λO_F = 𝔪𝔫⁻¹ with
𝔪, 𝔫 coprime integral. Then 𝔫O_K·𝔞′ = 𝔪O_K·𝔞; if 𝔭 | 𝔫, no prime above 𝔭 divides 𝔪O_K, so 𝔭O_K | 𝔞, contradicting primitivity. So 𝔫 = O_F,
symmetrically 𝔪 = O_F, and λ ∈ O_F^×.
(ii) *Shape normalization.* ℓ(β) := (log|σ_vβ|² − n⁻¹log Nβ)_v lies in H₀ := {t ∈ Rⁿ : Σt_v = 0}. For u ∈ O_F^×, N_{K/Q}u = N_{F/Q}(u)² = 1 and
ℓ(uβ) = ℓ(β) + λ(u), λ(u) := (2log|σ_vu|)_v. By Dirichlet, Λ_F := λ(O_F^×) is a lattice of full rank n−1 in H₀ and ker λ = {±1}. Fix a bounded
fundamental domain 𝒫 for Λ_F in H₀ and put r_F := sup_{t∈𝒫}max_v t_v (n = 2, 𝒫 = [−½,½)·λ(ε_F): r_F = log ε_F). For each integral 𝔞 ∈ C pick one
generator β_𝔞 of 𝔞𝔠 with ℓ(β_𝔞) ∈ 𝒫. If 𝔞 ≠ 𝔞′ are O_F-primitive then β_𝔞′/β_𝔞 ∉ F: otherwise it is a unit u by (i), λ(u) ∈ (𝒫 − 𝒫) ∩ Λ_F = {0},
u = ±1 and 𝔞 = 𝔞′. If N𝔞 ≤ N then |σ_vβ_𝔞|² = (Nβ_𝔞)^{1/n}e^{ℓ_v(β_𝔞)} ≤ (N·N𝔠)^{1/n}e^{r_F} for every v.
(iii) *Angle boxes.* Let k ≥ 1 be the least integer with kⁿ ≥ (4π)ⁿe^{nr_F}d_F N/√d_K, and cut (R/2πZ)ⁿ into kⁿ boxes Π_v[2πj_v/k, 2π(j_v+1)/k).
If the angle vectors (θ_v(β_𝔞))_v, (θ_v(β_𝔞′))_v of two ideals counted by π_C(N) lay in one box, then |sin(θ_v(β_𝔞) − θ_v(β_𝔞′))| < 2π/k for all v, and
Lemma 1 and (ii) would give N𝔠√d_K/d_F ≤ Π_v 2|σ_vβ_𝔞||σ_vβ_𝔞′||sin(…)| < Π_v 2(N·N𝔠)^{1/n}e^{r_F}(2π/k) = (4π/k)ⁿe^{nr_F}N·N𝔠 ≤ N𝔠√d_K/d_F.
So each box holds at most one,
π_C(N) ≤ kⁿ ≤ (1 + 4πe^{r_F}(d_FN/√d_K)^{1/n})ⁿ ≤ 2^{n−1}(1 + (4π)ⁿe^{nr_F}d_FN/√d_K). □

Remarks.
- 2.1 (size). The O_F-primitive ideals have Dirichlet series ζ_K(s)/ζ_F(2s), so their average number per class with N𝔞 ≤ N is
  ~ κ_KN/(ζ_F(2)h_K) = (2π)ⁿR_KN/(w_Kζ_F(2)√d_K) ≍_F N/√d_K. Lemma 2 is sharp up to a factor depending on F. Numerically (§6) the largest class
  count is 175 to 1.1·10⁴ times smaller than the bound (the constant (4π)ⁿe^{nr_F}d_F is crude).
- 2.2 (Z-primitivity fails; why [D2] stopped). Call 𝔞 Z-primitive if 𝔞 ⊄ mO_K for all m ≥ 2 (the notion of [D2] §5). For a prime 𝔭 of F of
  degree one, 𝔭O_K is Z-primitive, has norm p², and lies in j(Cl_F). So some class of j(Cl_F) (≤ h_F classes) contains ≫_F √N/log N Z-primitive
  ideals of norm ≤ N, against 1 + N/√d_K. Any large sieve over Ĉl_K (or 𝒵_K, on which all these ideals give χ = 1) with weight 1 on Z-primitive
  ideals has constant ≫_F h_K√N/log N, which at N = √d_K exceeds N + √d_K by d_K^{1/4−o(1)}. This is the quartic analogue of [D2] Proposition 7
  (loss d^{1/3}/log²d for pure cubic fields): the ideals of a subfield crowd into few classes. For CM fields the cure is to weight by the
  O_F-content and use 𝒵_K; this is what Lemma 2 allows.
- 2.3 (principal class). If β ∈ O_K∖F then β − βᶜ = Δ(β,1) ∈ 𝔇∖{0} and |σ_v(β − βᶜ)| ≤ 2|σ_vβ|, so N(β) ≥ 4^{−n}|N(β − βᶜ)| ≥ d_K/(4ⁿd_F²). So every
  prime ideal of degree one in the principal class has norm ≥ d_K/(4ⁿd_F²): exponent 1 is necessary in Theorem B.
- 2.4 ([D2] Proposition 8). Lemma 2 with 𝔞 = 𝔤O_K𝔞₀ and Cauchy–Schwarz over 𝔤 gives the L² pair bound (P_K, Δ) of [D2] §5 with
  Δ ≪_{F,ε} (M d_K)^ε for every CM field over F. So [D2]'s conditional Proposition 8 holds unconditionally for these fields.

## 2. The weighted large sieve

**Lemma 3.** For 0 < ε ≤ ½, N ≥ 1 and complex c_𝔞 supported on N𝔞 ≤ N,

  Σ_{χ∈𝒵_K}|Σ_𝔞 c_𝔞χ(𝔞)|² ≤ ζ_F(1+2ε)·h_K·2^{n−1}(1 + (4π)ⁿe^{nr_F}d_FN/√d_K)·‖c‖²_ε ≤ (B_F/ε)log^{2n−1}(3d_K)(N + √d_K)‖c‖²_ε,

with B_F depending only on F. The same holds for every coset ψ₀𝒵_K (apply it to c_𝔞ψ₀(𝔞)). Both terms are attained up to F-constants
(c = δ_{𝔞=O_K} gives |𝒵_K| ≥ h_K/h_F; c = 1 on all O_F-primitive ideals of norm ≤ N gives ≥ (κ_K/ζ_F(2) + o(1))N, from χ = 1).

*Proof.* Write 𝔞 = 𝔤O_K·𝔟 with 𝔟 O_F-primitive. By Fact (a), Σ_𝔞 c_𝔞χ(𝔞) = Σ_𝔟 B_𝔟χ(𝔟) with B_𝔟 := Σ_𝔤 c_{𝔤𝔟}, and by Cauchy–Schwarz
|B_𝔟|² ≤ ζ_F(1+2ε)Σ_𝔤 N𝔤^{1+2ε}|c_{𝔤𝔟}|², so Σ_𝔟|B_𝔟|² ≤ ζ_F(1+2ε)‖c‖²_ε. Orthogonality on Cl_K/j(Cl_F) gives
Σ_{χ∈𝒵_K}|Σ_𝔟B_𝔟χ(𝔟)|² = |𝒵_K|Σ_{C̄}|Σ_{[𝔟]∈C̄}B_𝔟|² ≤ |𝒵_K|·max_{C̄}π_{C̄}(N)·Σ_𝔟|B_𝔟|², where π_{C̄}(N) ≤ |j(Cl_F)|·max_Cπ_C(N) and
|𝒵_K||j(Cl_F)| = h_K. Insert Lemma 2. The second form uses ζ_F(1+2ε) ≤ C_F/ε and Fact (c). □

**Lemma 3′ (coefficient sums; [C] Lemma 3 with the content weight).** For k ≥ 1, ε > 0, M ≥ 2:
Σ_{N𝔞≤M}N_{F/Q}(𝔤(𝔞))^{1+2ε}d_k(𝔞)²/N𝔞 ≤ C(F,k,ε)M^{2ε}, and the same sum without 1/N𝔞 is ≤ C(F,k,ε)M^{1+2ε}.
*Proof.* 𝔞 = 𝔤O_K𝔟 with N𝔞 = N𝔤²N𝔟, so N𝔤^{1+2ε}/N𝔞 ≤ M^εN𝔤^{−1}N𝔟^{−1}. d_k(𝔞) ≤ d_k(𝔤O_K)d_k(𝔟) and d_k(𝔤O_K) ≤ d_k(𝔤)² (prime by
prime: split d_k(𝔭^a)², inert d_k(𝔭^a), ramified d_k(𝔓^{2a}) = binom(2a+k−1,k−1) ≤ binom(a+k−1,k−1)², since every k-tuple of sum 2a is a sum of
two k-tuples of sum a). The number of ideals of norm m is ≤ d_{2n}(m) (in K) and ≤ d_n(m) (in F), so the sum is
≤ M^εΣ_{N𝔤≤√M}d_k(𝔤)⁴/N𝔤·Σ_{N𝔟≤M}d_{k²}(𝔟)/N𝔟 ≤ M^ε(1+log M)^{C(n,k)}. The second sum is ≤ M times the first. □

## 3. Zero density for each field

**Theorem A.** Let F be a fixed totally real field of degree n. For every CM field K over F, every ε > 0, t ∈ R (T = |t|+3):
- (a) Σ_{χ∈Ĉl_K}N(σ;t,χ) ≤ C(F,ε)(d_KTⁿ)^{e(σ)+ε} for 0.6 ≤ σ ≤ 1, e(σ) = 2(1−σ)/(3−2σ);
- (b) Σ_{χ∈Ĉl_K}N(σ;t,χ) ≤ c₀h_K log(d_KT^{2n}) ≤ C_F√d_K log^{2n}(3d_KT) for every σ.

So Σ_χ #{ρ : β ≥ σ, |γ| ≤ T} ≤ C(F,ε)(d_KTⁿ)^{e(σ)+ε}T for T ≥ 3. The family Ĉl_K has h_K = d_K^{1/2±o(1)} members (Brauer–Siegel, F fixed)
and conductors d_K: a thin family, and (a) is |Ĉl_K|^{4(1−σ)/(3−2σ)+o(1)}, the exponent of [R2] Theorem 1′ for imaginary quadratic fields.

*Proof.* (b) is (F4) with Fact (c). For (a), [R2] Proposition 1 (which runs [D2] Proposition 6) with this dictionary; every step not listed
is verbatim.
1. 𝔔 := d_K. The Z-trivial family is 𝒵_K; "Z-trivial" is replaced by "F-trivial", and the two facts used in [D2] Prop 6 (χ(m𝔞) = χ(𝔞) and
   χ̄(𝔞) = χ(𝔞̄)) are Fact (a): χ(𝔤O_K𝔞) = χ(𝔞) and χ̄(𝔞) = χ(𝔞ᶜ). Ĉl_K is the union of ≤ h_F cosets ψ₀𝒵_K; treat each and add.
2. f = 1: every χ is primitive with A(χ) = d_K (Fact (b)), so E_χ = 1 and the 𝔢 | f^∞ splitting of [C] Lemma 4, step 3, is void.
3. (LS) is Lemma 3 with B = B_F log^{2n−1}(3d_K) (Prop 6 and [R2] Prop 1 use only B ≤ 𝔔^{o(1)}) and Λ := √d_K = 𝔔^{1/2}. The coefficient sums are
   Lemma 3′ in place of [C] Lemma 3.
4. AFE (F5): the gamma factor is Γ_C(s)ⁿ for every χ, so in [D2] Prop 6, step 2, U := |u|+3, the dyadic weights are
   μ_j = C(B₁)min(1, (2^{j−1}/(√d_K Uⁿ))^{−B₁}), and the dual sum is handled by χ̄(𝔞) = ψ̄₀(𝔞)ξ(𝔞ᶜ) (ξ ∈ 𝒵_K) and 𝔞 ↦ 𝔞ᶜ. The result is

   Σ_{χ∈ψ₀𝒵_K, χ≠1}|L(½+iu,χ)M_X(½+iu,χ)|² ≤ C(F,ε)d_K^ε(X d_K^{1/2}(|u|+3)ⁿ)^{1+ε}  (1 ≤ X ≤ d_K).   (2.1′)

5. Zero detector [C] Lemma 5: unchanged (L(s,χ) is entire for χ ≠ 1). Zero counts per window: (F4).
6. [R2] Prop 1, class II: by (2.1′) with X = 1, Σ_{χ≠1}|L(½+iv,χ)|² ≤ C(F,ε)d_K^{1/2+ε₀}(|v|+3)^{n+ε₀}; by Lemma 3 and Lemma 3′ the mollifier of length
   V := √d_K has Σ_χ|M_V(½+iv,χ)|² ≤ C(F,ε)d_K^{1/2+ε₀}. Hence |𝒳_II| ≤ C(F,ε)Y^{1/2−σ}d_K^{1/2+ε₀}T^{n/2+ε₀}. Class I is verbatim (blocks of length
   M ≥ V, where Lemma 3 has constant ≤ 3(B/ε₀)M) and gives |𝒳_I| ≤ C(F,ε)(d_KT)^{ε₀}Z^{2−2σ+2ε₀}, Z = 2Y log²Y.
7. Choose Y := (d_KTⁿ)^{1/(3−2σ)} (≥ V since 3−2σ ≤ 2). Then Y^{2−2σ} = Y^{1/2−σ}(d_KTⁿ)^{1/2} = (d_KTⁿ)^{e(σ)}. Rename ε. □

The Carlson version ([D2] Prop 6 itself) gives the exponent 2σ(1−σ) in d_K; (a) is smaller for every σ ∈ (½,1).

*Comparison.* Pascadi–Thorner (arXiv 2508.14888, Thm 1.3; any finite set of GL₁ forms over any number field), with base field K, 𝒮 = Ĉl_K,
Q ≍ d_K, |𝒮| = h_K = d_K^{1/2+o(1)}: (d_K^{−2}Q³|𝒮|⁴)^{(1−σ)/(3−2σ)+o(1)} = d_K^{3(1−σ)/(3−2σ)+o(1)}. Theorem A has exponent 2(1−σ)/(3−2σ): smaller by
the factor 2/3 for every σ. Kowalski–Michel type and log-free bounds (Weiss; Thorner–Zaman) are weaker for fixed σ < 1.

## 4. Zero density over a family of CM fields

Family. Fix δ₀ ∈ (0,1). 𝒦(X) = 𝒦_F^{δ₀}(X) := the CM fields K/F with X/2 < D ≤ X, 2 ∤ 𝔡_{K/F}, and N𝔭 > D^{δ₀} for every prime 𝔭 | 𝔡_{K/F}
("δ₀-rough"; for δ₀ ≥ ½ this means 𝔡_{K/F} is prime). Then 𝔡 is squarefree with ω(𝔡) < 1/δ₀; a given 𝔡 belongs to ≤ |Cl_F(𝔡∞)[2]| ≤ C(F,δ₀)
fields, so |𝒦(X)| ≤ C(F,δ₀)X. Prime 𝔡 alone give ≫_F X/log X fields: K = F(√−π) with π ≫ 0, (π) prime, π ≡ −1 mod 4O_F (then −π is a square
mod 4, so 𝔡_{K/F} = (π)); such (π) form a ray class mod 4∞ (Hecke). For K, let Ĝ_g(K) := {χ ∈ Ĉl_K : χᶜ = χ}
(genus characters) and Ĝ_c(K) := Ĉl_K∖Ĝ_g(K).

**Proposition C.** For X, T ≥ 3, 0 ≤ σ ≤ 1 and ε > 0:
Σ_{K∈𝒦(X)}Σ_{χ∈Ĝ_c(K)}N(σ,T;χ) ≤ C(F,δ₀,ε)(X T^{(3n+2)/5})^{10(1−σ)/(3−2σ)+ε}.

*Proof.* [R] Proposition 1 over F. For χ ∈ Ĝ_c(K), π_χ := AI_K^F(χ) is a cuspidal automorphic representation of GL₂(A_F) (χ ≠ χᶜ) with
L(s,π_χ) = L(s,χ); at each real place of F it is the weight-one representation D₁ (parameters 0, 1); its Satake parameters are
{χ(𝔓), χ(𝔓ᶜ)} at split 𝔭 and {±χ(𝔭O_K)^{1/2}} at inert 𝔭, so GRC holds. C(π_χ) = 12ⁿd_F²D ≤ Q := 12ⁿd_F²X. π_χ ≅ π_{χ′} (χ′ ∈ Ĉl_{K′}) forces
π_χ ≅ π_χ ⊗ η_{K′/F}; a cuspidal π has at most three such quadratic twists, and inside one K, π_χ ≅ π_{χ′} iff χ′ ∈ {χ, χᶜ}: each π occurs ≤ 6
times. 𝒮 := {π_χ}; |𝒮| ≤ Σ_K h_K ≤ C X^{3/2}log^{2n−1}(3X) (Fact (c)). Apply (HT) with m = 2: V = Q^{2+ε}T^{2n+1+ε}|𝒮| ≤ C X^{7/2+2ε}T^{2n+1+ε}.
L^{ur}(s,π_χ) = L(s,χ)Π_{𝔭|𝔡}(1 − χ(𝔓_𝔭)N𝔭^{−s}) (𝔓_𝔭² = 𝔭O_K), so |L^{ur}(½+iv,π_χ)| ≤ 2^{1/δ₀}|L(½+iv,χ)|, and (2.1′) with X = 1 gives
I₂ ≤ 4^{1/δ₀}Σ_K C(F,ε)d_K^{1/2+ε}T^{n+1+ε} ≤ C X^{3/2+2ε}T^{n+1+ε} (HT's convexity bound is |𝒮|Q·T^{n+1} ≍ X^{5/2}). Choose
Y := (V·X^{3/2}T^{n+1})^{1/(3−2σ)} (log Y ≍ log V). Then Y^{1/2−σ}(V·I₂)^{1/2} ≤ X^{O(ε)}Y^{2(1−σ)}, Y^{1−σ} ≤ Y^{2(1−σ)}, and
Y^{2(1−σ)} = (X^{5+O(ε)}T^{3n+2+O(ε)})^{2(1−σ)/(3−2σ)}. Rename ε. □

## 5. Primes in every ideal class, for almost all fields

**Lemma 4 (genus characters).** Let K/F be CM with 2 ∤ 𝔡. Let 𝒢(K) be the group of finite-order Hecke characters ψ of F that are unramified at
every finite 𝔭 ∤ 𝔡 and whose restriction to O_𝔭^× is trivial or the quadratic character of (O_F/𝔭)^× for 𝔭 | 𝔡. Then η ∈ 𝒢(K); ψ ↦ ψ∘N_{K/F} maps
𝒢(K) onto Ĝ_g(K) with kernel {1, η}; and L(s,ψ∘N_{K/F}) = L(s,ψ)L(s,ψη) (primitive characters). The conductor of ψ ∈ 𝒢(K) is Π_{𝔭∈S(ψ)}𝔭 (times
real places), S(ψ) := {𝔭 | 𝔡 : ψ ramified at 𝔭}; ψ is unramified at all finite primes iff ψ ∈ 𝒢₀ := Ĉl_F⁺. |𝒢(K)| ≤ h_F⁺2^{ω(𝔡)}.

*Proof.* (onto) Let χ ∈ Ĝ_g(K) and M := H_K^{ker χ} (H_K the Hilbert class field). As χ∘c = χ, c(ker χ) = ker χ, so M/F is Galois, and c acts
trivially on Gal(M/K) ≅ Cl_K/ker χ (𝔞ᶜ𝔞⁻¹ ∈ ker χ). So Gal(M/F) is a central extension of Z/2 by a cyclic group, hence abelian. M/K is unramified,
so an inertia group of M/F meets Gal(M/K) trivially and has order ≤ 2; it is nontrivial only above places ramified in K/F. Extend χ to a character
of Gal(M/F); by class field theory it is a Hecke character ψ of F with χ = ψ∘N_{K/F}, unramified outside 𝔡∞, and at 𝔭 | 𝔡 its restriction to
O_𝔭^× has order ≤ 2, hence is 1 or the Legendre character (𝔭 odd). (into) For ψ ∈ 𝒢(K): at 𝔓 | 𝔭 | 𝔡, N(O_𝔓^×) = (O_𝔭^×)² (ramified quadratic,
𝔭 odd), on which ψ_𝔭 is trivial; elsewhere ψ is unramified. So ψ∘N ∈ Ĉl_K, and (ψ∘N)ᶜ = ψ∘N. (kernel) ψ∘N = 1 iff ψ is trivial on the norms from
the idele class group of K, i.e. ψ ∈ {1, η}. (L-functions) Ind_K^F(ψ∘N) = ψ ⊕ ψη. At 𝔭 | 𝔡 exactly one of ψ, ψη is unramified, and its Euler factor
is that of ψ∘N at 𝔓_𝔭 (a uniformizer of K_𝔓 has norm −π_𝔭u, a local norm, on which ψ_𝔭 = (ψη)_𝔭). (rest) A character of the narrow ray class group
mod 𝔡∞ that is trivial on all O_𝔭^× is a character of Cl_F⁺; so 𝒢(K)/𝒢₀ embeds in Π_{𝔭|𝔡}{1, Legendre}. □

So, for K ∈ 𝒦(X), the genus L-functions are the fixed L(s,ψ), ψ ∈ Ĉl_F⁺ (the factors of ζ_M, M = H_F⁺), and the L(s,ψ) with ψ ∈ 𝒢(K)∖𝒢₀,
of conductor norm > D^{δ₀}. No genus factor of small nontrivial conductor occurs; this is why the family is δ₀-rough (§7.1).

**Theorem B.** Let F be a totally real field of degree n, δ₀ ∈ (0,1) and A > 3/2. There are effective c, κ > 0 and X₀ (depending on F, δ₀, A)
such that for X ≥ X₀ all K ∈ ∪_{Y≤X}𝒦_F^{δ₀}(Y) outside a set of size ≤ X^{1−c} satisfy: for every C ∈ Cl_K, every x ≥ D^A and every
H ∈ [xD^{−κ}, x],

  #{𝔓 ∈ C prime : N𝔓 = p prime, x < p ≤ x+H} ≥ H/(16h_K log 2x).   (5.1)

For n = 2 and A ≥ 1.55 one may take κ = min(10⁻³, δ₀/4, κ₀) and c = min(4·10⁻³, δ₀/3), where κ₀ = κ₀(F,A) > 0 is the effective constant of
step (e).
Since x ≥ d_K^A implies x ≥ D^A, (5.1) holds for x ≥ d_K^A and H ≥ xd_K^{−κ/2} (D ≥ d_F²).

*Proof.* [R] §2 with the changes of [R2] Theorem 2′ and the following simplifications (no Siegel-zero factor ν_D occurs). The statement for A
implies it for every A′ > A; take 3/2 < A ≤ 1.6.

*Parameters.* σ₁(A) := 3/2 − 1/A; φ(σ) := (10/(3−2σ) − A)(1−σ), decreasing on [0.8,1], φ(σ₁) = 2(2−A) < 1. Fix σ₀ ∈ [0.8, σ₁) with φ(σ₀) < 1;
ε′ := (1−φ(σ₀))/4; s₁ := (1−σ₀)(A − 2/(3−2σ₀)) > 0; κ := min(κ₀, ε′/(n+2), s₁/(n+2), δ₀/(n+2)); ε₀ := ε₁ := κ/10; σ_g := 0.995; σ_j := σ₀ + jε′/A
(0 ≤ j ≤ J, J minimal with σ_J ≥ 1−ε′/A); T_X := 2X^{κ+ε₁}. params_cm.py checks every inequality below (n = 2, 3, 4, 6; A = 1.55, σ₀ = 0.84:
φ(σ₀) = 0.96412, ε′ = 0.00897, s₁ = 0.00558; κ = 10⁻³ passes).

*Exceptional sets* (X a power of 2, K ∈ 𝒦(X)).
- E₁(X): some j has N^c_K(σ_j, T_X) := Σ_{χ∈Ĝ_c(K)}N(σ_j,T_X;χ) > X^{A(1−σ_j)−2ε′}. Markov and Proposition C (ε = ε₀), with φ(σ_j) ≤ 1−4ε′ and
  T_X^{(3n+2)·2(1−σ)/(3−2σ)} ≤ T_X^{2(3n+2)/7} (σ ≥ 0.8): #E₁(X) ≤ CJX^{1−2ε′+2ε₀+(2(3n+2)/7)(κ+ε₁)} ≤ CJX^{1−ε′}.
- E₂(X): some ψ ∈ 𝒢(K)∖𝒢₀ has a zero with β > σ_g, |γ| ≤ T_X. For dyadic Z ≥ (X/2)^{δ₀}, let Ψ(Z) be the set of ψ ∈ ∪_K𝒢(K)∖𝒢₀ with conductor
  norm in (Z, 2Z]. Each is determined by its squarefree conductor 𝔣 up to Ĉl_F⁺ (two such characters with the same local components on all O_𝔭^×
  differ by a character of Cl_F⁺), so |Ψ(Z)| ≤ C_F Z, and ψ ∈ 𝒢(K) needs 𝔣 | 𝔡_{K/F}, which holds for ≤ C(F,δ₀)X/Z fields K ∈ 𝒦(X). (HT) with
  m = 1 and its convexity bound for I₂ (V ≍ Z^{2+ε}T^{n/2+1+ε}, Y = (Q^{1/2}|𝒮|T^{n/2+1}V)^{1/(3−2σ)}):
  Σ_{ψ∈Ψ(Z)}N_ψ(σ_g,T_X) ≪ (Z^{7/2}T_X^{n+2})^{2(1−σ_g)/(3−2σ_g)+ε} ≤ Z^{0.035+ε}T_X^{0.01(n+2)}. So, with ε := δ₀/100 and κ ≤ δ₀/(n+2),
  #E₂(X) ≤ Σ_Z C(X/Z)Z^{0.035+ε}T_X^{0.01(n+2)} ≤ CX^{1−0.965δ₀+ε+0.011(n+2)κ} ≤ X^{1−0.9δ₀}.
- Fields with D < D₀ are put into the exceptional set. Summing over dyadic X: c := min(ε′, 0.9δ₀)/2 works (for A ≥ 1.55: c = min(4·10⁻³, δ₀/3)).

*Decomposition.* Fix K ∈ 𝒦(X)∖(E₁ ∪ E₂) with D ≥ D₀, C ∈ Cl_K, x ≥ D^A, H ∈ [xD^{−κ}, x], T₁ := 2D^{κ+ε₁} ≤ T_X. Take w as in [R] §2.3
(w₀ ∈ C_c^∞((0,1)), 0 ≤ w₀ ≤ 1, w₀ = 1 on [¼,¾], w(u) = w₀((u−x)/H)): W(1) ≥ H/2, |W(ρ)| ≤ Hx^{β−1}, |W(ρ)| ≤ C_mHx^{β−1}(2x/(H|ρ|))^m. Let
S_C := Σ_{𝔞∈C}Λ_K(𝔞)w(N𝔞) and ψ_w(ψ) := Σ_𝔤Λ_F(𝔤)ψ(𝔤)w(N𝔤) (ideals 𝔤 of F). Orthogonality and Lemma 4 (each genus χ = ψ∘N contributes
S_χ(w) = ψ_w(ψ) + ψ_w(ψη), and χ̄(C) = ψ̄(N_{K/F}𝔞_C) = (ψη)‾(N_{K/F}𝔞_C) for 𝔞_C ∈ C prime to 𝔡) give

  h_K S_C = Σ_{ψ∈𝒢(K)} ε_ψ(C)ψ_w(ψ) + Σ_{χ∈Ĝ_c(K)} χ̄(C)S_χ(w),  |ε_ψ(C)| = 1, ε_1(C) = 1.   (5.2)

Apply (I6) to each term: the pole of ζ_F (ψ = 1) gives W(1); the rest is −Σ_ρW(ρ) + R. Bound the parts:
- (a) zeros of χ ∈ Ĝ_c(K) with β ≥ σ₀, |γ| ≤ T₁: as K ∉ E₁ and x^{σ−1} ≤ 2^AX^{−A(1−σ)}, E_a ≤ HJ2^AX^{−ε′} ([R] §2.4(a)).
- (b) zeros of all χ ∈ Ĉl_K with β < σ₀, |γ| ≤ T₁ (this includes the genus factors, whose zeros are zeros of L(s,ψ∘N)): per window |γ−k| ≤ 1,
  Theorem A(a) on [0.6,σ₀] gives x^{σ−1}log x·N(σ;k) ≤ C D^{−g(σ)+ε₀}(|k|+3)^{ne(σ)+nε₀}(A log D + 10), g(σ) := (1−σ)(A − 2/(3−2σ)) ≥ s₁ (both factors are
  positive and decreasing on [0.6,σ₀]); Theorem A(b) handles [½,0.6] and β < ½ with D^{1/2−0.4A+o(1)} ≤ D^{−0.1}. With ne(σ) ≤ 4n/9:
  E_b ≤ HD^{−s₁+ε₀+(1+4n/9+ε₀)(κ+ε₁)+o(1)} + HD^{−0.1+1.1(κ+ε₁)}.
- (c) zeros of any character with |γ| > T₁: each window holds ≤ c₀h_K log(d_K(|k|+3)^{2n}) zeros (all of Ĉl_K); with m ≥ 2/ε₁, E_c ≤ HD^{−1}.
- (d) zeros of ψ ∈ 𝒢(K)∖𝒢₀ with β ≥ σ₀, |γ| ≤ T₁: as K ∉ E₂, β ≤ σ_g, and there are ≤ |𝒢(K)|·CT₁log(DT₁) of them: E_d ≤ HD^{−0.005A+κ+ε₁+o(1)}.
- (e) zeros of ψ ∈ 𝒢₀ = Ĉl_F⁺ with β ≥ σ₀, |γ| ≤ T₁: these are zeros of the fixed ζ_M (M = H_F⁺, multiplicities added). Put L := log(T₁+3),
  y := log x ≥ A log D, so y/L ≥ A/(2(κ+ε₁)) for D ≥ D₀. By (ZFR) and (TZ) (N_M(σ,T₁) ≤ C_Me^{B_ML(1−σ)}):
  Σ x^{β−1} ≤ N_M(0,3)x^{β_M−1} + x^{σ₀−1}N_M(σ₀,T₁) + y∫_{c_M/L}^{1/2}C_Me^{−u(y−B_ML)}du ≤ o(1) + 2C_M exp(−c_My/(2L)) ≤ o(1) + 2C_M exp(−c_MA/(4(κ+ε₁))),
  using y ≥ 2B_ML. With κ₀ := min(c_MA/(4.4 log(64C_M)), A/(4.4B_M)) this is ≤ 1/32 + o(1): E_e ≤ H/16 for D ≥ D₀.
- (f) remainders: Σ|R| ≤ C(h_K + |𝒢(K)|)(H/x + x^{−1/2})log(d_Kx) ≤ HD^{−1} (x ≥ D^{3/2}, H ≥ xD^{−κ}).
*Conclusion.* (a), (b), (d) are ≤ HD^{−η} for some η = η(F,δ₀,A) > 0 by the choice of κ (params_cm.py: margins 0.0076 and 0.0045 at A = 1.55,
κ = 5·10⁻⁴). So h_KS_C ≥ H/2 − H/16 − 4HD^{−η} ≥ H/4 for D ≥ D₀. Ideals 𝔞 with N𝔞 = p^f, f ≥ 2 (primes of degree ≥ 2 and prime powers)
contribute ≤ C(n)(Hx^{−1/2} + log x)log^{2n+1}x ≤ H/(16h_K) to S_C (√x ≥ D^{3/4}, h_K ≤ D^{1/2+o(1)}, H/h_K ≥ D^{A−κ−1/2−o(1)}). So
Σ log p ≥ 3H/(16h_K) over the 𝔓 ∈ C with N𝔓 = p ∈ (x, x+H], which gives (5.1). □

**Corollary B1.** Let K be as in Theorem B (outside the exceptional set).
- Every ideal class of K contains a prime 𝔓 with N𝔓 = p prime, p ≤ 2D^A ≤ 2d_K^A; and ≫ x/(h_K log x) of them in (x, 2x] for every x ≥ D^A.
- There is β ∈ O_K with N_{K/Q}(β) prime and ≤ 2d_K^A (principal class). For every CM field K over F, every β ∈ O_K∖F has N(β) ≥ d_K/(4ⁿd_F²)
  (Remark 2.3), and an element of F never has prime norm. So the least prime norm from O_K lies in [d_K/(4ⁿd_F²), 2d_K^A].
- Form version: for every integral 𝔞 ∈ C⁻¹, the relative norm form β ↦ N_{K/F}(β) on the rank-2 O_F-module 𝔞 (totally positive definite)
  takes a value with N_{K/F}(β)O_F = 𝔭·N_{K/F}(𝔞), 𝔭 a prime of F of degree one with N𝔭 ≤ 2d_K^A (take 𝔓 ∈ C from (5.1) and 𝔓𝔞 = (β)).

**Corollary B2 (conditional; routine).** Assume (Hu) (Hu, Hexagon 2610.00216: L(s,ψ) ≠ 0 for Re s > 7/8, every Hecke character of every number
field), resp. (L) (Re s > B♭ = 0.874887). Then for every CM field K over F, every ideal class contains a prime of degree one with norm
≤ C(F,ε)d_K^{8/5+ε}, resp. d_K^{2/(3−2B♭)+ε} = d_K^{1.5997+ε}. *Proof.* [C] Theorem 2 with Theorem A in place of [C] Theorem 1 ([R2] §5:
x^σN(σ;k) ≤ C x d_K^{−(1−σ)(A−2/(3−2σ))+ε}(|k|+3)ⁿ on [0.6, B] and 2/(3−2σ) ≤ 2/(3−2B)), Fact (c) for the error terms, and Lemma 4 is not needed (all
zeros, genus ones included, satisfy β ≤ B). For COROLLARIES.md, not DISCOVERIES.md.

## 6. Numerical checks (all pass)

- reldet.gp: for K = F(√−π) with (F,π) = (Q(√5), y+24), (Q(√5), 7), (Q(√2), 2y+23), (Q(√3), 2y+25), (Q(√10), 4y+27), (Q(√5), 3y+4) (y² = 5, 2, 3, 10;
  for Q(√5), y² = y+1): N(𝔇_{K/F}) = d_K/d_F² in every case, and for 300 random pairs β, β′ ∈ O_K each, Δ = ββ′ᶜ − βᶜβ′ ∈ 𝔇_{K/F} (1800/1800);
  |NΔ|^{1/2}/(√d_K/d_F) is an integer ≥ 1 (Lemma 1).
- run1.gp + run2.gp (cmcount.gp): 10 fields with d_K from 1225 to 1.49·10⁸ (h_K up to 408, Cl_K up to Z/204 × Z/2, h_F up to 2): for
  N ∈ {√d_K/4, √d_K, 4√d_K, N_max} (N_max ≤ 6000) the largest number of O_F-primitive ideals of norm ≤ N in one class is ≤ 1/175 of the Lemma 2 bound, and
  h_K·max_Cπ_C(N)/(N + √d_K) ≤ 0.61 (the large sieve constant of Lemma 3 is ≍ N + √d_K). Z-primitive counts are larger (up to 1.6×), as Remark 2.2 predicts.
- leastp.gp: least prime of degree one in every class, 7 fields (d_K ≤ 8·10⁶): the largest is ≤ 0.003d_K, attained by the principal class in
  6 of 7 cases, and there it is 1.04–1.97 times the lower bound d_K/(16d_F²) of Remark 2.3 (e.g. d_K = 8008400: 21169 vs 20021).
- params_cm.py → params_cm.out: φ(σ₁(A)) = 2(2−A); all inequalities of Theorem B for A = 1.55, σ₀ = 0.84 and n ∈ {2,3,4,6}.

## 7. Where the method stops (pinned)

1. **All CM fields K/F, not only δ₀-rough ones.** If 𝔡_{K/F} has prime factors of small norm, Lemma 4 produces genus factors L(s,ψ) of small
   nontrivial conductor, up to 2^{ω(𝔡)} of them, one possibly with a Siegel zero shared by a positive proportion of fields. [R] handles the
   analogue (χ_d, d | D) with the Q-aspect log-free density and Deuring–Heilbronn repulsion of Thorner–Zaman (arXiv 2208.11123) for all Dirichlet
   characters of conductor ≤ Q. Over F we have only (TZ) for one congruence class group at a time; summing it over 2^{ω(𝔡)} moduli is not log-free,
   and the smooth part of 𝔡 exceeds D^{θ} for a positive proportion of 𝔡. Needed: a Q-aspect log-free zero density with Deuring–Heilbronn for Hecke
   characters of F (not found; An, arXiv 2004.14410, proves log-free densities for automorphic families over number fields; whether it
   covers this case, with the repulsion, was not checked). With it, Theorem B extends to all K/F with N𝔡 ≤ X (with a factor ν_K as in [R]).
2. **F varying ("all quartic CM fields with d_K ≤ X").** Lemma 2's constant (4π)ⁿe^{nr_F}d_F is ε_F²d_F for n = 2, up to e^{c√d_F}. Cutting 𝒫 into unit
   cubes replaces e^{nr_F} by C^n(1 + vol 𝒫) ≍ R_F (sketch), so Theorem A holds with a factor d_F^{O(1)}. But (HT) and (TZ) are used with constants
   depending on F, so Theorem B is for each fixed F. Since the quartic CM fields with d_K ≤ X and d_F > X^η number ≪ X^{1−η+o(1)}, an F-uniform form
   of (HT) would give the full family.
3. **Non-CM fields.** Lemma 1 holds verbatim for every quadratic extension K/F (Δ ∈ N_{K/F}(𝔠)𝔇_{K/F}); the count then needs hyperbolic
   coordinates at real places of F that split in K (as in [D2] Lemma 3) and the relative units (rank r_K − r_F), and the threshold stays √d_K
   because h_KR_K ≍ √d_Kκ_K. Sketch only, not proved here. For fields with **no quadratic subfield** (cubic, S₄/A₄ quartic) there is no involution:
   [D2] Proposition 7 (pure cubic) shows the Z-primitive count fails, and no substitute is known.
4. **A < 3/2.** At A = 3/2, σ = 5/6 the per-field sum 1 + e(σ), the family exponent 10(1−σ)/(3−2σ) and the Markov bar 1 + A(1−σ) meet at 5/4, exactly
   as in [R2] §4 (the CM family has the same numerology: X fields, X^{1/2} characters each, conductor X). Lowering 3/2 needs a family large sieve
   below the Rankin–Selberg threshold X^{7/2}.
5. **Every K, unconditionally.** Theorem A loses d_K^ε near σ = 1 (not log-free), so it gives no least-prime bound for an individual field ([D2] §6.1).

## 8. Novelty audit (2026-10-21)

- arXiv API, 46 distinct queries in three sweeps (axq.out, axq2.out, axq3.out; failed queries rerun until answered): CM field(s) × large sieve /
  zero density / least prime / prime ideal × class; quartic CM × prime / class group; class group characters × CM / large sieve; ideal class ×
  least prime (ideal); every ideal class × prime; Hilbert class field × least prime / Chebotarev × almost all; weight one / dihedral × Hilbert
  modular × zero density; CM extension × zero density; totally imaginary quadratic extension × prime; relative norm / binary Hermitian forms ×
  primes; primitive ideals × CM field; anticyclotomic × zero density; Linnik × CM field; least prime ideal × almost all / number field;
  zero density × class group L-functions / Hilbert modular; class group / ideal classes / Grössencharacters × large sieve; au:Thorner × class group.
- Relevant hits, none implying Lemmas 2–3, Theorem A, Proposition C or Theorem B:
  - Habegger–Pazuki (1506.02485, source read): for CM K over a fixed totally real F, the proportion of classes in a coset of Cl_K with an ideal of
    norm ≤ ε√d_K is ≪ ε^{1/2}, via subconvexity for class group L-functions. An average over classes; no per-class bound, no large sieve.
  - Sardari (1802.06193): every class of a number field has a prime ≪ h_K²log²d_K, under GRH. Thorner–Zaman (1604.01750, Thm 1.1 with L = H_K):
    ≪ d_K^{694} unconditionally, for every K.
  - Pascadi–Thorner (2508.14888, Thm 1.3): GL₁ families over any field: d_K^{3(1−σ)/(3−2σ)} for Ĉl_K (§3).
  - An (2004.14410), Pierce–Turnage-Butterbaugh–Wood (1709.09637), Lemke Oliver–Thorner–Zaman (2012.14422), Lemke Oliver–Smith (2405.08383):
    Chebotarev for almost all fields of fixed Galois group; Hilbert class fields have degree h_K → ∞ (not covered).
  - Deshouillers–Gun–Ramaré–Sivaraman (2210.11051): every ray class contains a product of three primes of small norm. Van Order (1903.06722):
    moments of class group twists over CM K/F, no zeros. Tomé (2307.16775): Hecke L-values of CM extensions.
- GXL evidence API, 3 claims (gxl.out): no supporting passage for a large sieve or zero density for class group characters of CM fields, nor for
  primes in every class of almost all CM fields (closest: Habegger–Pazuki, Sardari).
- Verdict: Lemmas 1–3, Theorem A, Proposition C, Theorem B and Corollary B1 are new as stated. The new ingredient is Lemma 1 (with Lemma 2); the
  rest is the project's machinery ([C], [D2], [R], [R2]) with (HT) and (TZ) over F.

## 9. Dead ends (NOGO)

- Z-primitive weights for CM fields (the [D2] §5 setting): fail by d_K^{1/4−o(1)} (Remark 2.2).
- Ordering primitive vectors by one angle ([D2] Lemma 1) does not extend to n ≥ 2 angles; boxes plus Lemma 1 replace it (sharp up to F-constants).
- All CM K/F: blocked by small genus conductors without a Q-aspect log-free density + Deuring–Heilbronn over F (§7.1).
- Uniformity in F, A₀ < 3/2, every K unconditionally: §7.2, §7.4, §7.5.
