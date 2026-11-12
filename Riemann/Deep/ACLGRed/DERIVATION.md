# ACLGRed derivation notes (aclgred, 2026-11-12)

Source: arXiv e-print 2611.13129 v1 (only version on the abs page, 12 Nov 2026), prolate-weil-fine15.tex, 1527 lines,
md5 bbd15bcd…; same file as ../ACLGFals and ../ACLGH2. "l." = line of that file. W = [1/λ, λ] (the window), L² = L²(W, d*u).
N_λ(<c) := number of eigenvalues of QW_λ below c, both ι-sectors.

## 1. The chain of Sec. 8 (statement, lines, hypotheses, verdict)

| # | statement | lines | uses | verdict |
|---|---|---|---|---|
| C1 | Thm TF: S ⊇ {∞} ∪ {p ≤ λ²}, f real C^∞, supp f ⊂ W, h = f*f*: Σ_{v∈S} W_v(h) = −Tr(θ(h)(1−P−P̂)); if Mf(±½) = 0, QW_λ(f,f) = Tr(θ(f)(1−P−P̂)θ(f)*) | 1119–1130 | [Co99], [CM] L. 2.35, Thm 2.36, §8 | cited; sign and cyclicity rechecked: correct |
| C2 | Lemma shift: QW_λ(f,f) = 2 log σ‖f‖² + Tr(θ(f)(1−P_σ−P̂_σ)θ(f)*) | 1133–1148 | C1 | correct: Tr θ(h)1_(1,t] = h(1) log t (kernel h(u/v)), same for P̂ via F_S^{-1}θ(h)F_S = θ(h̃); h(1) = ‖f‖² |
| C3 | Lemma H1 (Range P_t ∩ Range P̂_t = 0), H2 (κ_n ∈ (0,1), O(n^{-1/2})), H0 (θ(h)P_Son trace class), Halmos (eq. halmos), eq. Kn, eq. decompK | 1178–1220 | Coburn, [CCqi], [CCMsl] | Kn rechecked: on V_n, 1−P−P̂ = 1_{V_n} − \|e⟩⟨e\| − \|f⟩⟨f\| has eigenvalue −κ on e+f, +κ on e−f; Gram inverse correct |
| C4 | (H2′): Σ_{n>n0} κ_n(t)‖θ(f)a_n‖² ≤ ε_t(n0)Λ(f), ε_t(n0) → 0, uniformly in f | 1236–1239 | hypothesis (open) | — |
| C5 | Lemma transport; α_nm = ⟨f, e^λ_m * e_n*⟩, β_nm = ⟨f, f^λ_m * f_n*⟩ | 1265–1275 | supp f ⊂ W | correct (⟨g,h⟩ = (g*h*)(1), commutative group) |
| C6 | Prop pairexp: QW_λ = 2 log σ‖f‖² + ‖P_Son θ(f)*‖²_HS + Σ_{n≤n0} B_n + R_{n0}; B_n ≥ Σ_m w_n(κ_n − κ^λ_m)(\|α\|²+\|β\|²); R_{n0} ≥ −ε_σ(n0)Λ(f) | 1278–1296 | C2, C3, C5, (H2′) | correct; w_n(κ^λ_m − κ_n) ≤ w_n(1−κ_n) = κ_n/(1+κ_n) ≤ ½ |
| C7 | Lemma archdom: Λ(f) ≤ QW_λ(f,f) + C_0λ‖f‖², C_0 absolute | 1306–1316 | eqs. Psi, Wp, WR | correct: \|W_02\| ≤ 2(λ−λ^{-1})‖f‖² (Cauchy–Schwarz, ∫_W u d*u = λ−λ^{-1}); \|ΣW_p\| ≤ 2‖f‖²Σ_{p^m≤λ²}log p·p^{-m/2} ≈ 4λ‖f‖² |
| C8 | Thm condpos: ε_σ(n0) ≤ min(1, log σ/(C_0λ)), M = {m: κ_m(λσ) > κ_{n0}(σ)} ⇒ QW_λ ≥ (log σ)/(1+ε)‖f‖² ≥ ½ log σ‖f‖² on W_λ(σ,n0), codim ≤ 2 + 2(n0+1)\|M\| | 1318–1330 | C6, C7, (H2′) | correct: (1+ε)QW ≥ (2 log σ − εC_0λ)‖f‖² ≥ log σ‖f‖²; codim count correct |
| C9 | Lemma schur: V ⊕ W orthogonal, D ≥ d‖w‖², d > 0 ⇒ ε_λ ≥ min(0, λ_min(𝔖)), 𝔖(v) = A(v) − ‖D^{-1/2}B*v‖²; Q ≥ 0 ⇔ 𝔖 ≥ 0 | 1348–1355 | — | correct (completing the square; d only makes D^{-1} bounded) |
| C10 | Prop monotone: ε_λ non-increasing; ε_λ ≥ −η(λ), η → 0 ⇒ QW_λ ≥ 0 ∀λ ⇔ RH | 1357–1362 | Weil criterion ([CCsel] (4)) | correct; only a sequence λ_k → ∞ with η(λ_k) → 0 is needed (ε_λ monotone) |
| C11 | Split V = V_low ⊕ V_tr ⊕ Cρ, (R1), (R2), bounds (schurred), Weyl ⇒ ε_λ ≥ −max(α,η) − β | 1364–1374 | C8, C9, (R1), (R2) | GAP G1, G2; GAP F2 (sec. 3, SIGNIFICANT); size claims G3, F1 |
| C12 | Prop reduction: (H2′)+(R1)+(R2), τ = λ^{-A}, A with α,β → 0, λ_min(𝔖\|_{V′}) ≥ −η(λ) → 0 ⇒ RH | 1376–1378 | C11, C10 | valid implication after G1 fix (η, β re-established on V ⊖ V_low); its (R2) contradicts Fact count (F2) |
| C13 | "block of dimension O(Aσ² log λ) + 2" | 1364, 1380 | Slepian (archimedean) | false: G3 (archimedean log slip), F1, F2, sec. 4 (semilocal) |

(R1) l. 1366: QW(g,h) = Ψ(g* * h) (rechecked: g**h(x) = ∫g(u)h(ux)d*u = F(x)), so θ(φ)Rad ⊂ Rad for test functions φ: correct.

## 2. Gaps of the reduction step (C11)

G1 (orthogonality). V_low, V_tr, Cρ are spans of different conditioning vectors: the sum is direct, not orthogonal. "λ_min(𝔖) ≥
−max(α,η) − β by Weyl's inequality" (l. 1374) needs V = V_low ⊕ V′ orthogonal. For a direct sum with cosine c of the minimal angle,
only λ_min(𝔖) ≥ −(max(α,η) + β)/(1−c). Fix: V′ := V ⊖ V_low (same dimension); the hypothesis of Prop. reduction must be placed on
this V′. ROUTINE.

G2 ((R2) as derived, l. 1367). "k = K − r, K ∈ Rad(QW), ‖r‖ ≤ Cτ^{1/2}(log λ)^C‖k‖, so that |QW_λ(k,v)| = |QW(r,v)| ≤ …‖k‖(QW_λ(v,v)
+ C_0λ‖v‖²)^{1/2}". The L² bound does not give the bilinear bound: E(S_0) ⊂ Rad(QW) and ME(f)(it) = ζ(½+it)Mf(it); multiplication
by a function nonzero a.e. has dense range, so Rad(QW) is L²-dense and dist_{L²}(k, Rad) = 0 for every k. The operative content of
(R2) is the bilinear bound itself (it needs the form size of r, e.g. by estimating pole, archimedean and prime parts of QW(r,v)
separately in absolute value, as [CK] do for QW(r_λ,r_λ)). A Cauchy–Schwarz route for r outside W would need the global QW to be
bounded below, which fails without RH (an off-line zero β+iγ gives QW(g,g)/‖g‖² ≥ −C T^{2β}log T only, unbounded as supp g grows).
So (R2) must be stated as the bilinear bounds; ROUTINE as a reformulation, not circular if proved term by term.

G3 (archimedean count, l. 1364). With the paper's own factor e^{−π²/(2 log(8πλ²σ²))} per index, the number of indices with
δ_m(λσ) ∈ (λ^{−A}, 1) is A log λ · 2 log(8πλ²σ²)/π² = (4A/π²) log²λ (1 + O(1/log λ)), not O(A log λ); Landau–Widom gives the same
order ((1/π²) log(2πλ²σ²) log(1/τ) for all n). So dim V_tr = O(Aσ² log²λ) in the archimedean case. Slip; no logical effect. Exact rungs (sec. 8) fit Landau–Widom:
even rungs in (1e-10, 1e-2) = 5, 6, 6 at t² = 20, 40, 80 against (1/2π²) log(2πt²) log(1e8) = 4.5, 5.2, 5.8; the per-even-index step in
log δ is 2π²/log(2πt²) (3.2 at t² = 80), not the paper's π²/(2 log(8πt²)) (0.65): constant off by ≈ 5, order log²λ unchanged.

## 3. F2: (R2) on V_low as defined contradicts the eigenvalue count (GAP, SIGNIFICANT; reviewer, sec. 9)

Definition (l. 1364): V_low = span of the conditioning vectors 1_W(e^λ_m * e_n*), 1_W(f^λ_m * f_n*), n ≤ n0, attached to modes m with
δ_m(λσ) ≤ τ (scale λσ, not λ).

(a) (R2) ⇒ dim V_low ≤ N_λ(≤ θ_λ), θ_λ = Cτ(log λ)^C. Proof: (R2) gives |QW_λ(k,k′)| ≤ θ_λ‖k‖‖k′‖ on V_low, so QW_λ ≤ θ_λ on V_low;
min–max: the (dim V_low)-th eigenvalue is ≤ θ_λ.

(b) dim V_low ≥ M_low := #{m: δ_m(λσ) ≤ τ}. Archimedean model (S = {∞}, e_n = √2 u^{1/2}h_{n,σ}(u), e^λ_m = √2 u^{1/2}h_{m,λσ}(u)),
n = 0 alone: for x ≤ λ, c_m(x) := (e^λ_m * e_0*)(x) = 2√x ∫_0^σ h_{m,λσ}(xw) h_{0,σ}(w) dw (substitute u = xw; the cut u ≤ λσ is
inactive for x ≤ λ), √x times an entire function of x. If Σ c_m·c_m(x) = 0 on W, then H := Σ c_m h_{m,λσ} (even, entire) satisfies
∫_0^σ H(xw)h_{0,σ}(w)dw = 0 for all x; with H(z) = Σ a_k z^{2k}: a_k μ_k = 0, μ_k = ∫_0^σ w^{2k}h_{0,σ}(w)dw > 0 (h_{0,σ} > 0, ground
state), so H = 0 and all c_m = 0. Hence the window restrictions are linearly independent. Semilocal pair: it contains the
archimedean near-self-dual directions (l. 1198, [CCMsl] Prop. 4.1: η_S(h_{m,t}) ∈ Range P_t, F_Sη_S(h) = η_S(ĥ) ∈ Range P̂_t), so its
M_low is at least the archimedean one (min–max, up to the factor ‖η_S‖ ≤ ∏_{p∈S}(1−p^{-1/2})^{-1} against δ^{1/2}); independence of
their window convolutions is the analogue of the moment argument (not proved here for S ≠ {∞}).

(c) Counts. Archimedean: M_low = #{even m ≤ 4λ²σ² with δ_m(λσ) ≤ τ} = 2λ²σ² − O(A log²λ) (Slepian, G3). Fact count (l. 742–749,
λ² ≤ 150) and ../ACLGFals (λ² = 2…20): N_λ(< 10^{-2}) = 2λ² − 3. With θ_λ < 10^{-2}: (R2) needs 2λ²σ² − O(A log²λ) ≤ 2λ² − 3,
false for σ² > 1 + O(A log²λ/λ²). Ratio: dim V_low/N_λ(θ) ≥ σ²(1 − o(1)); for σ = 2, a factor 4.

(d) Mechanism. A mode deep at scale λσ occupies the phase-space box [−λσ,λσ]²; window functions resolve Mellin frequency t ≈ 2πxy;
the near-null (minuscule) part of QW_λ is the undersampled band |t| ≤ γ = 2πλ² (deficit density (1/π)log(γ/|t|), total 2λ², l. 1390).
Deep scale-λσ modes with phase-space support at |xy| > λ² (number ≈ 2λ²(σ² − 1 − 2 log σ), the "corners") feed Mellin frequencies above γ
(oversampled) into V_low, and the radical element θ(e_n*)K^λ_m of (R1) is not concentrated in W: it is supported on (0, ∞) (supp e^λ_m =
(0, λσ], supp e_n* = [1/σ, ∞)), with u^{∓1/2} tails (coefficients ∝ ∫h_{m,λσ}, h_{n,σ}(0)). So "in the radical up to the same tails as in
(bridge)" (l. 1366) is false for σ > 1: the tails of (bridge) are alias tails of size δ_4^{1/2} for E-images of scale-λ prolates
supported in [−λ,λ]; here the part outside W is not small. Theorem upperK (l. 573–577), the [CK] input that (R2) is to extend (l. 1369),
is for H_K(λ), prolates at scale λ: it covers at most the N_λ(θ) ≈ 2λ² directions of the scale-λ ladder.

## 4. F1: a rigorous lower bound for the transition-and-pole block

Given (H2′): Thm condpos gives QW_λ ≥ ½ log σ on V^⊥ ∩ form domain (Remark (e)), so min–max gives dim V ≥ N_λ(< ½ log σ).
Given (R2) on some V_low ⊂ V: dim V_low ≤ N_λ(≤ θ_λ) (3a). With V′ := V ⊖ V_low (G1):
   dim V′ ≥ N_λ(< ½ log σ) − N_λ(≤ θ_λ) = #{eigenvalues of QW_λ in (θ_λ, ½ log σ)}.
This holds for every admissible split. The paper's claim dim V′ = O(Aσ² log λ) + 2 is therefore equivalent to: QW_λ has only
O(Aσ² log λ) eigenvalues in (θ_λ, ½ log σ). Fact cond (l. 879–891): the non-minuscule edge μ_λ ≈ 2/λ → 0, so for every fixed σ the
non-minuscule eigenvalues enter (θ, ½ log σ) as λ grows; how many is the size of V′ (sec. 5, counts). Galerkin eigenvalues are min–max
upper bounds, so Galerkin counts below a level are lower bounds for N_λ.

## 5. Sizes with superpolynomial n0 (semilocal count of ../ACLGH2)

ACLGH2 (semiclassical, not proved): n0 ≥ N_{S,σ}(κ♮) ≥ Σ_{κ>κ♮}κ² ≈ 2σ²F_S, \|M\| ≥ 2λ²σ²F_S, codim ≥ 8λ²σ⁴F_S², log F_S = Σ_{p≤λ²}(log p)²/p
+ O(1) = 2 log²λ + O(log λ), so codim = exp((1+o(1)) log²λ²) = λ^{(4+o(1)) log λ}.
- Which estimate needs polynomial dimensions? None of C8–C12 as implications: Lemma schur, Weyl and α, β are dimension-free once (R2)
  is assumed uniformly on V_low. Explicitly (my bookkeeping): ‖D^{-1/2}B*k‖² ≤ θ_λ(1 + 2C_0λ/log σ)‖k‖²/… so α = O(τλ^{1+o(1)}/log σ);
  β ≤ 3Cτ^{1/2}(log λ)^C(1 + 2C_0λ/log σ) sup_{v∈V′}‖v‖_form/‖v‖ with ‖v‖²_form ≤ (QW_λ(v,v)/‖v‖² + C_0λ)‖v‖² = O(λ)‖v‖² if the
  conditioning vectors have Mellin content below exp(λ^{O(1)}) (the semiclassical horizon is 2πt²exp(Σ_{p≤λ²} 2 log p/(√p−1)) = exp((4 + o(1))λ),
  so log of it is O(λ)); hence β = O(τ^{1/2}λ^{3/2+o(1)}): A > 3 suffices, independent of dim V.
- Polynomial dimensions are needed only on the route the paper names for (R2) (l. 1369–1370): per-vector [CK] estimates transferred to
  the span. Theorem upperK carries a factor K (Cauchy–Schwarz over K prolates, eq. formcomp); transferring ‖r_j‖ ≤ ε_0‖k_j‖ from N spanning
  vectors to their span costs (N/g)^{1/2}, g the least eigenvalue of the normalized Gram matrix. With N = 2(n0+1)M_low ≥ 4σ²F_S·M_low:
  ε_0(N/g)^{1/2} → 0 needs τ^{1/2}F_S^{1/2} → 0 (times g^{-1/2}·poly), F_S^{1/2} = exp(log²λ + O(log λ)) = λ^{log λ + O(1)}; the
  inequality that breaks for every fixed A is λ^{−A/2}·(4σ²F_S)^{1/2} → ∞; threshold A > 2 log λ + O(1), i.e. τ ≤ exp(−2 log²λ − O(log λ))
  ≈ F_S^{-1}. Theorem upperK's factor K (not √K) with K ≍ N gives the same order, τ·F_S → 0.
- Lowering τ to exp(−c log²λ) is consistent with Thm condpos (τ does not enter it) and with (R2) as a hypothesis (only the threshold moves);
  it does not shrink V′: V_tr grows as τ falls (archimedean O(cσ² log³λ)), and the semilocal thin-pocket modes (κ ≈ √(A/2π) ∈ (κ_{n0},
  1−τ), the bulk of the 2t²F_S count) are in V_tr for every τ ≤ 1 − O(1). So dim V′ is quasi-polynomial whatever τ.
- Remark condpos(b), l. 1335: the archimedean sizes are right (n0 = 4σ² + O(log λ), \|M\| = 2λ²σ²(1+o(1)), codim O(σ⁴λ²)); the semilocal
  sentence "both polynomial in λ under a polynomial decay of the angles" is wrong in its natural reading: the angles do decay polynomially
  for each S (Lemma H2, κ_n = O(n^{-1/2})), but the constant C_S in κ_n ≤ C_S n^{-1/2} must satisfy C_S² ≥ κ♮²N_S(κ♮) ≳ 2t²F_S κ♮², with
  κ♮ ≍ 1/λ (needed: ε ≤ log σ/(C_0λ), ε ≈ A_R κ_{n0}, ACLGH2 sec. 3.3), i.e. quasi-polynomial (semiclassical). True only if "polynomial
  decay" means κ_n(t) ≤ λ^{O(1)}n^{-a} uniformly in S ∋ p ≤ λ², which the semiclassical count contradicts.

## 6. Circularity (iv)

None besides Sec. 6 (Remark l. 1035). C1–C10 use no sign of QW (Remark condpos(a)). The hypothesis on 𝔖|_{V′} is necessary for RH
(RH ⇒ QW_λ ≥ 0 ⇒ 𝔖 ≥ 0 on V, Lemma schur), so it is a reformulation of RH on V′, not circular. Two places where a hidden sign could
enter: (R2)'s second bound |QW_λ(k,k′)| ≤ Cτ(log λ)^C does not follow from its first (v = k′ gives only Cτ^{1/2}λ^{1/2}(log λ)^C); it
needs QW_λ(k,k′) = QW(r,r′) and the form size of r, r′ (G2). F2 uses only QW_λ ≤ θ′ on V_low with θ′ → 0, so either bound suffices.
(R2)'s factor (QW_λ(v,v) + C_0λ‖v‖²)^{1/2} is unconditional for v ∈ W (Lemma archdom ⇒ QW_λ + C_0λ ≥ Λ ≥ 0); for r outside W only a
term-by-term (absolute-value) bound is admissible (G2). The sentence "in the computed range their conditioned Schur complement is positive"
(l. 1374) is numerical evidence, not used in a proof.

## 7. Weakest extra hypothesis (v)

Given (H2′), RH ⇔ λ_min(𝔖_λ on all of V) ≥ −o(1) along a sequence (Lemma schur, Prop monotone): the reduction's content is only the split.
The weakest split hypothesis that closes it and is consistent with F1/F2:
 (R2*) V_low ⊂ V^+ := V + E_λ(H_{K_τ}(λ)) (adding conditioning vectors keeps Thm condpos), V_low = E_λ(H_{K_τ}(λ)) ⊕ odd analogue, the
 scale-λ ladder up to δ_{4K}(λ) ≤ τ (≤ N_λ(θ) dimensions), with the bilinear bound sup|QW_λ(k,v)|/(‖k‖‖v‖_form) ≤ λ^{−3/2−ε};
 (S*) λ_min(𝔖_λ|_{V^+ ⊖ V_low}) ≥ −η(λ), η → 0.
(R2*) is the uniform-in-index [CK] estimate the paper asks for, at scale λ, plus the Gram bound g_K(λ) ≥ λ^{−C} (≥ 0.151 for λ² ≤ 50,
../ACLGFals 2.3). (S*) lives on dim(V^+ ⊖ V_low) ≥ #{eigenvalues of QW_λ in (θ_λ, ½ log σ)} directions (sec. 4), at most quasi-polynomially
many under the ACLGH2 count.

## 8. Computation: QW_λ on V_low (code/vlow.py, logs/vlow_*.json)

Conditioning vectors 1_W(e^λ_m * e_n*) for the archimedean pair (exact pair for λ² ≤ 2, where no prime enters and, l. 1388, "all the
hypotheses of Theorem condpos are known"; a model of the semilocal pair for λ² = 5, 10), via the formula of 3(b); prolates by the
Bouwkamp matrix in double precision; QW_λ = Galerkin matrix of ../WeilSimple/weil_matrix.py on E_N (N = 120; N = 240 at λ² = 2 agrees
to 4 digits); "dim" = number of Gram eigenvalues above cut × largest (normalized vectors), QW_λ compressed to that subspace of V_low.
| case | M_low = #{m: δ_m(λσ) ≤ 1e-10} | vectors | dim at cut 1e-4 / 1e-8 / 1e-11 | QW eigenvalues < 1e-2 | ≥ 0.05 | largest | N_λ(<1e-2) |
|---|---|---|---|---|---|---|---|
| λ² = 2, σ = 2, n = 0 | 11 | 11 | 11 / 11 / 11 (least Gram eigenvalue 7.8e-4) | 1 / 1 / 1 | 10 / 10 / 10 | 3.03 | 1 |
| λ² = 2, σ = 2, n ≤ 14 | 11 | 88 | 11 / 13 / 15 | 1 | 10 / 12 / 14 | 3.54 | 1 |
| λ² = 2, σ = 3, n = 0 | 30 | 30 | 11 / 21 / 30 | 1 | 10 / 20 / 29 | 4.61 | 1 |
| λ² = 5, σ = 2, n = 0 (model) | 34 | 34 | 13 / 22 / 30 | 3 / 6 / 7 | 9 / 15 / 23 | 3.47 | 7 |
| λ² = 5, σ = 2, n ≤ 14 (model) | 34 | 272 | 32 / 38 / 39 | 7 | 25 / 31 / 32 | 3.85 | 7 |
| λ² = 10, σ = 2, n = 0 (model) | 73 | 73 | 15 / 31 / 47 | 4 / 9 / 12 | 10 / 22 / 35 | 3.75 | 17 |
At λ² = 2, σ = 2 the eleven V_low eigenvalues are 0.0013, 0.074, 0.661, 0.927, 1.19, 1.38, 1.56, 1.78, 1.94, 2.80, 3.03, i.e. V_low is
close to the span of the lowest eigenvectors of QW_{√2} (Galerkin: 0.0013, 0.0739, 0.661, 0.927, 1.19, 1.38, …): one near-null
direction, ten of energy 0.07–3. Projection loss onto E_N ≤ 6% (λ² = 2), ≤ 27% (λ² = 10): the λ² = 5, 10 rows are model values.
Exact rung counts (code/prolcount.py): #{even m: δ_m(t) ≤ 1e-10} = 5, 14, 34, 53, 73, 83, 152 at t² = 5, 10, 20, 30, 40, 45, 80; Galerkin
N_λ(<1e-2) = 7, 17, 29, 37, 47, 56 at λ² = 5, 10, 16, 20, 25, 30. So for σ = 2: 34 vs 7, 73 vs 17, 152 vs 37 (ratio ≈ σ²).

## 9. Reviewer (recorded) and status

Reviewer verdicts on the decisive arguments (brief: .scratch/aclgred/reviewer_brief.md): G1 GAP, ROUTINE (also directness of the spans is
not established; η, β must be re-established on V ⊖ V_low); F1 CORRECT, ROUTINE (Galerkin interval counts diagnostic, not certified);
F2 GAP, SIGNIFICANT (missing: a rigorous bound N_λ(≤θ) ≤ 2λ² + o(λ²) at the moving threshold; finite numerics do not supply it; the
global-support mechanism alone does not prove loss of near-radicality); G3 GAP, ROUTINE (spacing must hold down to λ^{-A}; keep the cutoff
δ < 1 − κ_{n0}(σ)); Prop. reduction: the block argument survives with orthogonal blocks; "vacuous" not established. Hence F2 is reported
as an inconsistency with the paper's own count (Fact count), not as a proof that (R2) is false.
Route to a rigorous refutation: the hypothesis set implies RH, so it suffices to prove RH ⇒ N_λ(≤θ_λ) < #{m: δ_m(λσ) ≤ τ} for large λ.
Under RH, QW_λ(f,f) = Σ_γ |Mf(iγ)|² for every f ∈ L²(W) (Ψ(h) = Σ_ρ Mh(ρ − ½), pole term included), Mf ∈ PW_{log λ}; the needed statement is a Landau-type deficit bound: the
near-null space of the zero-sampling form on PW_{log λ} has dimension (1/π)∫_0^{2πλ²}(log λ² − log(t/2π))dt·(1 + o(1)) = 2λ²(1 + o(1)).

## 10. The τ that avoids F2

#{m: δ_m(λσ) ≤ τ} ≤ N_λ(θ) ≈ 2λ² needs τ ≤ δ_{4λ²}(λσ) = exp(−4πλ²σ² I(σ_{4λ²}) + O(log λ)) (uniform asymptotics, l. 486),
exponentially small in λ², not λ^{−A} or exp(−c log²λ). With such τ, V_low is (approximately) the scale-λ ladder, α, β = O(τ^{1/2}λ^{3/2+o(1)})
are tiny, and V_tr contains the scale-λσ modes 4λ² ≲ m ≲ 4λ²σ²: dim V′ ≥ M − N_λ(θ) ≈ 2λ²(σ² − 1) for n = 0 alone (archimedean model,
moment argument, with the count), consistent with F1.

## 11. Follow-up: can the count in F2 be proved from RH? (aclgred, 12 Nov, second session)

Answer: not as things stand (an open obstacle, not an impossibility proof). Under RH the near-null count is bounded BELOW by the multiplicity excess of the zeros
(Prop. L, proved). So the upper bound F2 needs, N_λ(≤ θ_λ) < (2 − ε)λ²σ², would prove that at most O(λ²) of the ≈ 2λ² log λ zeros
up to 2πλ² are multiple. That is open even under RH. The deep-mode count of (b) is proved for the semilocal pair (11.5); injectivity
is proved only for the archimedean pair.

### 11.1 Inputs
(I1) Explicit formula in the normalization (Psi): for f ∈ C_c^∞((1/λ, λ)) (complex), h = f * f*, Ψ(h) = Σ_ρ Mh(ρ − ½) (zeros with
multiplicity); Mh(s) = Mf(s)·conj(Mf(−s̄)). Under RH, ρ = ½ + iγ: QW_λ(f,f) = Σ_γ m_γ |F(γ)|², F(t) := Mf(it), γ distinct ordinates,
m_γ multiplicities (the paper uses this at l. 1374: "the net value QW_λ(ρ,ρ) is the sum over the zeros").
(I2) Variables: y = log u ∈ [−ℓ, ℓ], ℓ = log λ, g(y) = f(e^y), F(t) = ∫ g(y)e^{ity}dy, ‖f‖² = ‖g‖² = (1/2π)∫|F|²; |F(t)|² ≤ 2ℓ‖g‖².
(I3) Concentration operator A (time [−ℓ_1, ℓ_1], angular band [−T, T]) on L²(−ℓ_1, ℓ_1), kernel sin(T(x−y))/(π(x−y)), ⟨Ag,g⟩ =
(1/2π)∫_{−T}^{T}|G|². tr A = 2ℓ_1T/π. tr A² = ∫_{−2ℓ_1}^{2ℓ_1}(2ℓ_1 − |u|)sin²(Tu)/(π²u²)du, so tr(A − A²) ≤ 2/π² + (2/π²)∫_0^{2ℓ_1T}sin²v/v dv
≤ (2/π²) log(2ℓ_1T) + 3/π² (2ℓ_1T ≥ 1). Since Σ_{μ≤½}μ ≤ 2 tr(A − A²): #{μ_k > ½} ≥ tr A − 2tr(A − A²) ≥ 2ℓ_1T/π − (4/π²)log(2ℓ_1T) − 6/π².
(I4) Zeros (unconditional, Riemann–von Mangoldt with S(T) = O(log T)): N(T) = (T/2π)log(T/2πe) + O(log T); N(T+H) − N(T) ≤
(H/2π)log(T+H) + C_N log(T+H). N_d(T) = number of distinct ordinates in (0, T]; E(T) := N(T) − N_d(T) ≥ 0 (multiplicity excess).
(I5) Band-limited plateau: ω ∈ C_c^∞(−1,1) even, ω(0) = 1, ω̌(s) = ∫ω(x)e^{−isx}dx; ψ_δ(s) = (δ/2π)ω̌(δs) (∫ψ_δ = 1, ‖ψ_δ‖_1 = B_ω,
∫_{|s|>R}|ψ_δ| ≤ c_k(δR)^{1−k}); Φ := 1_{[−T−H, T+H]} * ψ_δ. Then Φ is the Fourier transform of φ_Φ = 2π·(sin((T+H)y)/(πy))·(inverse
transform of ψ_δ), a C_c^∞ function supported in (−δ, δ); |Φ − 1| ≤ ε_1 := c_k(δH)^{1−k} on [−T, T]; |Φ(t)| ≤ c_k(δ(|t| − T − H))^{1−k}
for |t| > T + H.

### 11.2 Proposition L (RH)
Let ℓ_1 = ℓ − δ, T ≥ 2, H ≥ 1, k ≥ 3, ε_1 ≤ ¼, Σ_tail := Σ_{|γ|>T+2H} m_γ|Φ(γ)|², and 4ℓΣ_tail ≤ θ(1 − ε_1)². Then
   N_λ(≤ θ) ≥ 2ℓ_1T/π − (4/π²)log(2ℓ_1T) − 6/π² − 2N_d(T + 2H).
Proof. E := span of the eigenfunctions of A with μ > ½ (dim by I3); E′ := {g ∈ E: G(γ) = 0 for the 2N_d(T+2H) distinct ordinates
|γ| ≤ T + 2H}, dim E′ ≥ dim E − 2N_d(T+2H). For g ∈ E′ let F := GΦ and f the window function with Mf(i·) = F: its log-time function
g * φ_Φ is C^∞ with support in [−ℓ_1 − δ′, ℓ_1 + δ′] ⊂ (−ℓ, ℓ), so f ∈ C_c^∞((1/λ, λ)), the core where (I1) holds; g ↦ f is linear and
injective (Φ ≢ 0 is entire). By (I1), (I2): QW_λ(f,f) = Σ_{|γ|>T+2H} m_γ|G(γ)|²|Φ(γ)|² ≤ 2ℓ‖g‖²Σ_tail, and ‖f‖² = (1/2π)∫|G|²|Φ|² ≥
(1 − ε_1)²⟨Ag,g⟩ ≥ ½(1 − ε_1)²‖g‖². So QW_λ ≤ θ on a subspace of dimension dim E′; min–max. ∎
Tail: Σ_tail ≤ Σ_{j≥1} 2[N(T+(j+2)H) − N(T+(j+1)H)]·c_k²(δjH)^{2−2k} ≤ C_k (δH)^{2−2k} H log(T + H) (k ≥ 3).

### 11.3 Corollary L (RH), at the deficit maximum T_0 = 2πλ²
Take T = T_0, H = λ, δ = λ^{−1/2}, k = A + 3, θ_λ ≥ λ^{−A}. Then ε_1 ≤ c_kλ^{−(k−1)/2}, 4ℓΣ_tail ≤ C_kλ^{2−k}log²λ ≤ θ_λ for λ ≥ λ_0(A), and
2ℓ_1T_0/π = 4λ²ℓ − 4λ^{3/2}, 2N_d(T_0 + 2λ) ≤ 2N_d(T_0) + O(λ log λ), 4λ²ℓ = 2N(T_0) + 2λ² + O(log λ) (I4). Hence
   N_λ(≤ θ_λ) ≥ 2λ² + 2E(2πλ²) − O(λ^{3/2}).
(i) Fact count's lower half is a theorem under RH: N_λ(≤ θ_λ) ≥ 2λ² − O(λ^{3/2}) for every θ_λ ≥ λ^{−A}.
(ii) Where the upper bound fails: any bound N_λ(≤ θ_λ) ≤ 2λ² + B_λ gives E(2πλ²) ≤ ½B_λ + O(λ^{3/2}). The requested 2λ² + o(λ²) would give
E(T) = o(T) on T = 2πλ²; F2 at fixed σ needs N_λ(≤ θ_λ) < (2 − ε)λ²σ² and would give E(T) ≤ (σ² − 1)T/(2π)(1 + o(1)), i.e. all but a
fraction ≍ (σ² − 1)/log T of the zeros up to T simple. Known under RH: Montgomery's Σ_{0<γ≤T} m_γ (with multiplicity) ≤ (4/3 + o(1))N(T)
gives E(T) ≤ Σ(m − 1) ≤ Σ(m² − m)/2 ≤ (1/6 + o(1))N(T), and later improvements (Cheer–Goldston, Bui–Heath-Brown, ≈ 70% simple) keep
E(T) ≍ N(T). The gap is a factor log T. So F2 cannot be closed from RH by counting; any closure by counting proves a new
simple-zeros theorem.
(iii) Close pairs count too: let G vanish at one ordinate per cluster of diameter ≤ ρ. Then |G(γ′)| ≤ ρ·sup|G′| ≤ ρℓ(2ℓ)^{1/2}‖g‖ at the
other members, and N_λ(≤ θ_λ) ≥ 2λ² + 2E_ρ(2πλ²) − O(λ^{3/2}) as soon as 4B_ω²ρ²ℓ³E_ρ ≤ θ_λ(1 − ε_1)²/2. Here E_ρ is the number of zeros within ρ of a
lower one. So an upper bound also needs E_ρ(T_0) = o(λ²) at ρ ≈ λ^{−A/2−1}(log λ)^{−2}, a statement on extremely small gaps (Montgomery's
pair correlation covers Fourier support [−1,1] only and does not give it).

### 11.4 What the upper bound would need beyond (ii), (iii)
A sufficiency (sampling) estimate on the oversampled band t > T_0: functions of PW_{log λ} that are small at all zeros above T_0 have
little mass there, modulo o(λ²) dimensions. Beurling–Landau-type sufficient conditions need maximal gaps below the Nyquist spacing π/log λ.
The RH gap bound (π + o(1))/log log t (Goldston–Gonek) is far larger, and the zeros' local count deviates by S(t) = O(log t/log log t).
So no available theorem supplies it. Unlike (ii), this is a missing method, not a proven obstruction.

### 11.5 Status of (b) for the semilocal pair
(b1) proved: #{m: δ_m(t) ≤ τ} ≥ #{even n: δ_n^{arch}(t) ≤ τ_1}, τ_1 := τ² b_S²/(6B_S²), b_S = ∏_{p∈S}(1 + p^{−1/2})^{−1},
B_S = ∏_{p∈S}(1 − p^{−1/2})^{−1}. Inputs: η_S(h) ∈ Range P_t for h supported in [−t,t] ([CCMsl] Prop. 4.1, l. 1094–1095); F_Sη_S = η_S F_∞ (l. 1095);
‖η_S(h)‖² = (1/2π)∫|ζ_S(½+it)|²|M[u^{1/2}h](it)|²dt (M E_S(h)(s) = ζ_S(s + ½)M[u^{1/2}h](s), S-integer sum), so b_S‖h‖ ≤ ‖η_S(h)‖ ≤ B_S‖h‖
(norms up to the even-function factor). Proof: H_0 ⊕ H_2 = span of h_{n,t}, n ≡ 0, 2 mod 4, δ_n ≤ τ_1 (orthogonal in L²(R)); for h = h_0 + h_2 put
x = η_S(h) ∈ Range P_t and y = η_S(ĥ_0 − ĥ_2) = F_S η_S(h_0) − F_S η_S(h_2) ∈ F_S(Range P_t) = Range P̂_t. Then ‖x − y‖ ≤ B_S(‖h_0 − ĥ_0‖ + ‖h_2 + ĥ_2‖)
≤ B_S(4τ_1)^{1/2}‖h‖ ≤ (B_S/b_S)(4τ_1)^{1/2}‖x‖ < τ‖x‖. Here ‖h_0 − ĥ_0‖² = Σ\|c_n\|²((1 − χ_n)² + 1 − χ_n²) = Σ\|c_n\|²·2δ_n ≤ 2τ_1‖h_0‖²
by (leak), (chi) (same for h_2 + ĥ_2).
So ‖P̂_t x‖ ≥ (1 − τ)‖x‖ on η_S(H_0 ⊕ H_2), and min–max for P_tP̂_tP_t gives the count (η_S injective). For S = {∞} ∪ {p ≤ λ²σ},
log(B_S/b_S) = 2Σ_{p∈S}p^{−1/2} + O(1) = O(λ/log λ), so log(1/τ_1) = O(λ). Uniform rung asymptotics (l. 486, [Du]: log(1/δ_n) = 2γI(σ_n)
+ O(log γ) for σ_n ≤ σ_0 < 1, 2n + 1 ≈ (2γ/π)J(σ_n)) give δ_n ≤ e^{−cλ²} ≤ τ_1 for all even n ≤ (γ/π)J(σ_0). That is λ²σ²J(σ_0)(1 − o(1))
of them, and J(σ_0) → 2 as σ_0 → 1. Hence M_low ≥ (2 − ε)λ²σ² for every ε > 0, λ ≥ λ_0(ε), for the semilocal pair.
(b2) not proved for S ≠ {∞}: dim V_low ≥ M_low needs injectivity of g ↦ (1_W(g * e_n*))_{n ≤ n0} on the deep spectral subspace. The moment
argument of 3(b) uses analyticity of the archimedean pair vectors. The semilocal ones are Range P_t eigenvectors built with U_fin
(dilations by p^k with weights ≍ p^{−|k|/2}), and their analyticity is unknown; the dilated archimedean transforms do not sum
analytically off the real axis.

### 11.6 What is rigorous now
Theorem R (archimedean pair, or semilocal pair under (b2)). If (R2) holds on V_low as defined (l. 1364), then for every ε > 0 and λ ≥ λ_0(ε):
   N_λ(≤ θ_λ) ≥ dim V_low ≥ (2 − ε)λ²σ².
ACLG's hypotheses imply RH (Prop. reduction + Prop. monotone), and under RH Corollary L gives N_λ(≤ θ_λ) ≥ 2λ² + 2E(2πλ²) − O(λ^{3/2}). So the hypothesis
set forces a near-null space σ²(1 − ε) times the deficit count 2λ² of l. 1390 and of Fact count. It is inconsistent if and only if
the near-null count stays below (2 − ε)λ²σ² for large λ. That is a theorem about the zeros (simplicity at the level E(T) = O(T),
small gaps at scale λ^{−A/2−1}, and a sampling estimate above T_0), open under RH. So F2 stays a GAP, now with its exact missing input,
and is not upgraded to INCONSISTENT HYPOTHESES.

### 11.7 Reviewer (second run, brief logs/reviewer_brief2.md)
Prop. L CORRECT, SIGNIFICANT (regime H ≥ 1, 0 < δ < ℓ, ε_1 < 1, k > 3/2). Corollary L CORRECT, ROUTINE (integer k ≥ A + 3).
Consequence GAP, ROUTINE: the count ⇒ E-bound implication is correct (Montgomery's sum counted with multiplicity); "best known" is a literature
claim; a necessary condition is an open obstacle, not impossibility. (b1) CORRECT, ROUTINE (given its inputs, including uniform asymptotics
down to τ_1). Status GAP, not INCONSISTENT HYPOTHESES: E(T) = O(T) is necessary, not sufficient, for the upper count.
