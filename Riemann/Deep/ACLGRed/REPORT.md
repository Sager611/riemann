GAP: Prop. reduction (l. 1364–1378). After the orthogonality fix (V′ := V ⊖ V_low, routine) the implication holds, but (R2) on V_low = conditioning vectors of modes deep at scale λσ forces N_λ(≤θ_λ) ≥ dim V_low ≥ (2 − ε)λ²σ² (archimedean pair; for the semilocal pair the deep-mode count is proved and injectivity is open), against the paper's own count N_λ ≈ 2λ². This does not become a contradiction under RH as things stand: under RH, N_λ(≤θ_λ) ≥ 2λ² + 2E(2πλ²) − O(λ^{3/2}) (E = multiplicity excess of the zeros, proved here), so the needed cap first requires E(T) = O(T) (known: E ≤ (1/6 + o(1))N(T)) and a sampling estimate. Repaired with an exponentially small τ, the transition-and-pole block has ≥ #{eigenvalues of QW_λ in (θ, ½ log σ)} ≈ 2λ²(σ² − 1) directions, not O(Aσ² log λ) + 2.

# ACLGRed: refereeing the reduction "(H2′) + (R1) + (R2) + Schur sign ⇒ ε_λ ≥ −η(λ) ⇒ RH" of Andrews–Connes–La Gatta (arXiv 2611.13129) (aclgred, 2026-11-12)

Source: e-print TeX prolate-weil-fine15.tex, 1527 lines, md5 bbd15bcd…, v1 (the abs page lists only v1, 9 Nov 2026); "l." = its lines.
Long derivations: [DERIVATION.md](DERIVATION.md). Inputs read first: ../ACLGFals/STATEMENTS.md, ../ACLGFals/REPORT.md, ../ACLGH2/REPORT.md.
Notation: W = [1/λ, λ], N_λ(<c) = number of eigenvalues of QW_λ below c (both ι-sectors), θ_λ = Cτ(log λ)^C (the (R2) scale).

## 1. The chain of Sec. 8
| step | statement | TeX | hypotheses | verdict |
|---|---|---|---|---|
| Thm TF | S ⊇ {∞} ∪ {p ≤ λ²}, Mf(±½) = 0 ⇒ QW_λ(f,f) = Tr(θ(f)(1−P−P̂)θ(f)*) | 1119–1130 | [Co99], [CM] Thm 2.36 | cited, signs rechecked: correct |
| Lemma shift | QW_λ = 2 log σ‖f‖² + Tr(θ(f)(1−P_σ−P̂_σ)θ(f)*) | 1133–1148 | Thm TF | correct (trace h(1)log t + h(1)log w recomputed) |
| H1, H2, H0, Halmos, Kn, decompK | pair structure; 1−P_t−P̂_t = P_Son + Σκ_n(\|b_n⟩⟨b_n\| − \|a_n⟩⟨a_n\|) | 1178–1220 | Coburn, [CCqi], [CCMsl] | correct (eigenvalues ∓κ on e ± f rechecked) |
| (H2′) | Σ_{n>n0}κ_n‖θ(f)a_n‖² ≤ ε_t(n0)Λ(f), ε_t(n0) → 0 | 1236–1239 | open | — |
| Prop pairexp | QW_λ = 2 log σ‖f‖² + ‖P_Son θ(f)*‖²_HS + Σ_{n≤n0}B_n + R_{n0}, B_n ≥ −(κ_n/(1+κ_n))Σ_{κ^λ_m>κ_n}(\|α\|²+\|β\|²) | 1278–1296 | (H2′) for R_{n0} | correct |
| Lemma archdom | Λ(f) ≤ QW_λ(f,f) + C_0λ‖f‖² | 1306–1316 | — | correct (pole ≤ 2(λ−λ^{-1})‖f‖², primes ≈ 4λ‖f‖²) |
| Thm condpos | ε_σ(n0) ≤ min(1, log σ/(C_0λ)) ⇒ QW_λ ≥ ½ log σ‖f‖² on W_λ(σ,n0), codim ≤ 2 + 2(n0+1)\|M\| | 1318–1330 | (H2′) | correct |
| (R1) | Rad(QW) is a convolution ideal | 1366 | — | correct for test functions (QW(g,h) = Ψ(g**h) rechecked) |
| (R2) | V_low near-radical: \|QW_λ(k,v)\| ≤ Cτ^{1/2}(log λ)^C‖k‖‖v‖_form, \|QW_λ(k,k′)\| ≤ θ_λ‖k‖‖k′‖ | 1367, 1369 | open | GAP G2 (L² form vacuous); GAP F2 (contradicts Fact count) |
| Lemma schur | ε_λ ≥ min(0, λ_min(𝔖)); Q ≥ 0 ⇔ 𝔖 ≥ 0 on V | 1348–1355 | V ⊥ W, D ≥ d > 0 | correct |
| block bounds | 𝔖\|_{V_low} ≥ −α, 𝔖\|_{V′} ≥ −η, ‖𝔖_{V_low,V′}‖ ≤ β ⇒ λ_min(𝔖) ≥ −max(α,η) − β (Weyl) | 1370–1374 | (R2), V′ = V_tr ⊕ Cρ | GAP G1 (needs V′ ⊥ V_low) |
| Prop reduction | (H2′)+(R1)+(R2), τ = λ^{-A}, λ_min(𝔖\|_{V′}) ≥ −η(λ) → 0 ⇒ ε_λ ≥ 0 ∀λ ⇒ RH | 1376–1378 | above | valid implication after G1 (η, β re-proved on V ⊖ V_low) |
| Prop monotone | ε_λ non-increasing; ε_λ ≥ −η(λ), η → 0 ⇒ QW_λ ≥ 0 ∀λ ⇔ RH | 1357–1362 | Weil criterion | correct; η(λ_k) → 0 on one sequence λ_k → ∞ suffices |
| size claim | V′ "of dimension O(Aσ² log λ) + 2" | 1364, 1380 | Slepian (archimedean) | false: G3, F1, F2, sec. 4 |

## 2. Findings on the reduction step
- G1 (orthogonality; GAP, ROUTINE). V_low, V_tr, Cρ are spans of overlapping families of conditioning vectors: the sum is direct at best,
  not orthogonal, and Weyl's inequality on the 2×2 block form needs ‖x+y‖² = ‖x‖² + ‖y‖². For a direct sum with minimal-angle cosine c
  one only gets λ_min(𝔖) ≥ −(max(α,η) + β)/(1 − c). Fix: V′ := V ⊖ V_low; η and β must then be proved on that complement.
- G2 ((R2) as derived, l. 1367). "‖r‖ ≤ Cτ^{1/2}(log λ)^C‖k‖, so that |QW(r,v)| ≤ …" is a non sequitur: E(S_0) ⊂ Rad(QW) and
  ME(f)(it) = ζ(½+it)Mf(it), so Rad(QW) is L²-dense and dist_{L²}(k, Rad) = 0 for every k. The content of (R2) is the bilinear bounds.
  The second bound does not follow from the first (v = k′ gives only Cτ^{1/2}λ^{1/2}(log λ)^C). For r outside W a Cauchy–Schwarz bound
  would need the global QW to be bounded below, false without RH. So the bounds must come term by term (pole, archimedean, primes in
  absolute value), as [CK] do for QW(r_λ, r_λ).
- F2 ((R2) vs V_low; GAP, SIGNIFICANT). V_low (l. 1364) = span of 1_W(e^λ_m * e_n*), 1_W(f^λ_m * f_n*), n ≤ n0, δ_m(λσ) ≤ τ.
  (a) (R2) ⇒ QW_λ ≤ θ_λ on V_low ⇒ (min–max) dim V_low ≤ N_λ(≤ θ_λ).
  (b) dim V_low ≥ #{m: δ_m(λσ) ≤ τ}: for the archimedean pair, (e^λ_m * e_0*)(x) = 2√x ∫_0^σ h_{m,λσ}(xw)h_{0,σ}(w)dw on x ≤ λ; a vanishing
  combination on W gives ∫_0^σ H(xw)h_{0,σ}(w)dw ≡ 0, H = Σc_m h_{m,λσ} = Σa_k z^{2k}, so a_k ∫_0^σ w^{2k}h_{0,σ} = 0, and h_{0,σ} > 0
  forces H = 0. The semilocal pair has at least these directions (l. 1198, [CCMsl] Prop. 4.1); independence of their convolutions is the
  analogue, not proved here.
  (c) The inequality that fails: #{even m: δ_m(λσ) ≤ τ} = 2λ²σ² − O(A log²λ) ≤ N_λ(≤ θ_λ) = 2λ² − 3 (Fact count, l. 742–749:
  below 10^{-10} each sector of QW_λ has one eigenvalue fewer than U_λ). Instances (σ = 2, exact rungs): 34 vs 7 (λ² = 5),
  73 vs 17 (λ² = 10), 152 vs 37 (λ² = 20). Both sides are exact rungs or the paper's own count; only the persistence of Fact count for large λ is assumed.
  (d) Mechanism: deep scale-λσ modes past the scale-λ plunge (4λ² < m ≲ 4λ²σ²) are not near-radical on W. The radical element θ(e_n*)K^λ_m
  of (R1) lives on (0, ∞), so "in the radical up to the same tails as in (bridge)" (l. 1366) fails for σ > 1. And Theorem upperK, the [CK]
  input (R2) is to extend (l. 1369), covers H_K(λ), i.e. scale λ, ≈ 2λ² directions.
  (e) Avoiding F2 needs τ ≤ δ_{4λ²}(λσ) = exp(−4πλ²σ²I(σ_{4λ²}) + O(log λ)) (l. 486), exponentially small, not λ^{−A}.
- F1 (lower bound for V′; CORRECT, ROUTINE). Thm condpos + min–max: dim V ≥ N_λ(< ½ log σ). (R2) + min–max: dim V_low ≤ N_λ(≤ θ_λ). Hence
  for every admissible split dim V′ ≥ #{eigenvalues of QW_λ in (θ_λ, ½ log σ)}. Fact cond (l. 879–891, edge μ_λ ≈ 2/λ → 0) puts
  non-minuscule eigenvalues below every fixed ½ log σ for large λ; the Galerkin counts of sec. 3 grow like λ² for σ = 3, 4 on λ² ≤ 30. With F2's τ, dim V′ ≥ #{m: δ_m(λσ) ≤ τ} − N_λ(θ) ≈ 2λ²(σ² − 1) (archimedean model, n = 0 alone).
- G3 (archimedean count; GAP, ROUTINE). With the paper's own factor e^{−π²/(2 log(8πλ²σ²))} per index, the rungs with δ_m(λσ) ∈
  (λ^{−A}, 1 − κ_{n0}(σ)) number (2/π²) A log λ · log(8πλ²σ²) = (4A/π²) log²λ(1 + o(1)), not O(A log λ); dim V_tr = O(Aσ² log²λ). Landau–Widom,
  (1/π²) log(2πt²) log(1/τ) over all n, fits the exact rungs (even rungs in (1e-10, 1e-2): 5, 6, 6 at t² = 20, 40, 80 against 4.5, 5.2, 5.8)
  and gives a per-even-index step of 2π²/log(2πt²) in log δ (3.2 at t² = 80) against the paper's π²/(2 log(8πt²)) (0.65): their constant
  is off by a factor ≈ 5, the log²λ order stands.

## 3. Computations
QW_λ on V_low (code/vlow.py; archimedean pair: exact for λ² ≤ 2, where l. 1388 says all hypotheses of Thm condpos are known; a model
for λ² = 5, 10; Galerkin QW_λ of ../WeilSimple/weil_matrix.py, N = 120, N = 240 at λ² = 2 equal to 4 digits; dim = Gram eigenvalues above
cut × largest; eigenvalues of QW_λ compressed to that subspace of V_low):
| case | #{δ_m(λσ) ≤ 1e-10} | vectors | dim (cut 1e-4 / 1e-8 / 1e-11) | eig < 1e-2 | eig ≥ 0.05 | largest | N_λ(<1e-2) |
|---|---|---|---|---|---|---|---|
| λ² = 2, σ = 2, n = 0 | 11 | 11 | 11 / 11 / 11 | 1 | 10 | 3.03 | 1 |
| λ² = 2, σ = 2, n ≤ 14 | 11 | 88 | 11 / 13 / 15 | 1 | 10 / 12 / 14 | 3.54 | 1 |
| λ² = 2, σ = 3, n = 0 | 30 | 30 | 11 / 21 / 30 | 1 | 10 / 20 / 29 | 4.61 | 1 |
| λ² = 5, σ = 2, n ≤ 14 (model) | 34 | 272 | 32 / 38 / 39 | 7 | 25 / 31 / 32 | 3.85 | 7 |
| λ² = 10, σ = 2, n = 0 (model) | 73 | 73 | 15 / 31 / 47 | 4 / 9 / 12 | 10 / 22 / 35 | 3.75 | 17 |
At λ² = 2, σ = 2, V_low is 11-dimensional (least normalized Gram eigenvalue 7.8e-4) with QW eigenvalues 0.0013, 0.074, 0.661, 0.927,
1.19, 1.38, 1.56, 1.78, 1.94, 2.80, 3.03: one near-null direction, ten of energy 0.07–3, against (R2)'s Cτ(log λ)^C. Projection loss onto
E_N ≤ 6% at λ² = 2, ≤ 27% at λ² = 10.
Eigenvalue counts for F1 (code/counts2.py on full Galerkin spectra, N = 240, 300 at λ² = 25, 30; diagnostic, not certified):
| λ² | N_λ(<1e-2) | (1e-2, ½ log 2) | (1e-2, ½ log 3) | (1e-2, ½ log 4) | (1e-2, 1) |
|---|---|---|---|---|---|
| 5 | 7 | 1 | 1 | 2 | 5 |
| 10 | 17 | 1 | 2 | 5 | 12 |
| 16 | 29 | 2 | 4 | 7 | 20 |
| 20 | 37 | 2 | 5 | 9 | 24 |
| 25 | 47 | 1 | 6 | 16 | 32 |
| 30 | 56 | 3 | 9 | 24 | 39 |
Exact rung counts #{even m: δ_m(t) ≤ 1e-10} (code/prolcount.py, mpmath): 5, 11, 14, 30, 34, 53, 73, 83, 152 at t² = 5, 8, 10, 18, 20, 30,
40, 45, 80 (t² = 8, 18 confirm the 11 and 30 of the λ² = 2 rows).

## 4. Size question with superpolynomial n0
- ACLGH2's count (semiclassical): n0 ≥ 2σ²F_S, \|M\| ≥ 2λ²σ²F_S, codim ≥ 8λ²σ⁴F_S², log F_S = Σ_{p≤λ²}(log p)²/p + O(1) = 2 log²λ + O(log λ):
  codim = exp((1 + o(1)) log²λ²).
- Prop. reduction as an implication is dimension-free: Lemma schur, Weyl and α, β do not see dim V once (R2) is uniform on V_low.
  Bookkeeping (mine): α = O(τλ^{1+o(1)}/log σ), β = O(τ^{1/2}λ^{3/2+o(1)}) if V′ vectors have Mellin content below exp(λ^{O(1)})
  (horizon 2πt²exp(Σ_{p≤λ²} 2 log p/(√p − 1)) = exp((4 + o(1))λ), so Λ ≤ O(λ)‖v‖²); A > 3 suffices.
- Polynomial dimensions are needed on the route to (R2) the paper names (l. 1369–1370): Theorem upperK's factor K (eq. formcomp) and the
  transfer from N spanning vectors to their span, factor (N/g)^{1/2}, N = 2(n0+1)M_low ≥ 4σ²F_S·M_low. The inequality that breaks for every
  fixed A: λ^{−A/2}(4σ²F_S)^{1/2} → ∞, since F_S^{1/2} = λ^{log λ + O(1)}; threshold A > 2 log λ + O(1), i.e. τ ≤ exp(−2 log²λ) ≈ F_S^{-1}.
- τ = exp(−c log²λ) is consistent with Thm condpos (τ does not enter it) and with (R2) as a hypothesis, but it does not shrink V′:
  V_tr grows as τ falls, the thin-pocket modes (κ ≈ √(A_pocket/2π) ∈ (κ_{n0}, 1 − τ), the bulk of 2t²F_S) are in V_tr for every τ, and
  F2 needs τ exponentially small. So dim V′ ≳ 2λ²(σ² − 1) (archimedean pair, given the count) and is quasi-polynomial under the semilocal count.
- Remark condpos(b) (l. 1335): archimedean sizes correct (n0 = 4σ² + O(log λ), \|M\| = 2λ²σ²(1 + o(1)), codim O(σ⁴λ²)). The semilocal sentence
  ("both polynomial in λ under a polynomial decay of the angles") is not correct as a size claim. The angles do decay polynomially for each S
  (Lemma H2: κ_n = O(n^{−1/2})), but κ_n ≤ C_S n^{−1/2} with C_S² ≥ κ♮²N_S(κ♮) ≳ 2t²F_S/(A_Rλ)², quasi-polynomial. It holds only if
  "polynomial decay" means constants λ^{O(1)} uniformly in S ∋ p ≤ λ², which the semiclassical count contradicts.

## 5. Circularity
None besides Sec. 6 (their Remark l. 1035). Thm TF through Thm condpos use no sign of QW (Remark condpos(a), checked). Lemma schur is an
identity; the sign hypothesis on V′ is implied by RH (RH ⇒ 𝔖 ≥ 0 on V), so Prop. reduction reformulates RH on V′ rather than assuming it.
The one hidden-positivity risk is G2: a Cauchy–Schwarz bound for QW(r, v) with r outside W would need the global form bounded below.

## 6. Weakest extra hypothesis
Given (H2′), RH ⇔ λ_min(𝔖_λ on V) ≥ −o(1) along a sequence (Lemma schur + Prop monotone); all content is in the split. The weakest split
consistent with F1, F2:
(R2*) V_low := E_λ(H_{K_τ}(λ)) ⊕ odd analogue, the scale-λ ladder with δ_{4K}(λ) ≤ τ, added to the conditioning space (V⁺ = V + V_low keeps
Thm condpos), with sup |QW_λ(k,v)|/(‖k‖‖v‖_form) ≤ λ^{−3/2−ε} on it: the uniform-in-index [CK] estimate at scale λ plus g_K(λ) ≥ λ^{−C}
(the Gram constant is ≥ 0.151 for λ² ≤ 50, ../ACLGFals 2.3);
(S*) λ_min(𝔖_λ|_{V⁺ ⊖ V_low}) ≥ −η(λ), η → 0, on a block of dimension ≥ #{eigenvalues of QW_λ in (θ_λ, ½ log σ)}.
Another repair: σ − 1 ≍ log²λ/λ² shrinks the excess 2λ²(σ² − 1) to O(log²λ), at margin ½ log σ ≍ log²λ/λ² and with (H2′) needed at
ε_σ(n0) ≤ log σ/(C_0λ) ≍ log²λ/λ³ (larger n0); the semilocal thin-pocket modes still make V′ quasi-polynomial (sec. 4).
A rigorous refutation of (R2) as defined needs only RH ⇒ N_λ(≤ θ_λ) < #{m: δ_m(λσ) ≤ τ} (the hypothesis set implies RH). Under RH,
QW_λ(f,f) = Σ_γ|Mf(iγ)|² for every f ∈ L²(W) (the pole term is part of Ψ), with Mf ∈ PW_{log λ}. So this is a Landau-type deficit bound: the near-null space of zero
sampling on PW_{log λ} has dimension (1/π)∫_0^{2πλ²} log(2πλ²/t) dt (1 + o(1)) = 2λ²(1 + o(1)) (the deficit of l. 1390). Sec. 11: this upper bound
would imply E(T) = O(T) for the multiplicity excess of the zeros, which is open under RH; only the matching lower bound is proved.

## 7. What the earlier audits imply for each hypothesis
| hypothesis | role in Sec. 8 | status |
|---|---|---|
| (H2′) | tail R_{n0} in Thm condpos | open, untestable as stated (ACLGFals); angle form of (H2ii) consistent at fixed S, T(n0,R) ≈ A_Rκ_{n0}, but n0 ≥ 2σ²F_S (ACLGH2, semiclassical) |
| (R1) | radical ideal | proved for test φ; for φ = e_n* (not compactly supported) continuity is needed (routine) |
| (R2) | α, β | Gram part holds (≥ 0.151, λ² ≤ 50, ACLGFals); L² form vacuous (G2); on V_low as defined contradicts Fact count (F2) |
| Schur sign on V′ | η | necessary for RH; block ≥ #eig(θ, ½ log σ) (F1); untested |
| Fact count | used in F2 | ACLG λ² ≤ 150; ACLGFals λ² = 2…20; recounted λ² = 5…30 here |
| LA, Lemma gap, QW_λ ≥ 0 | Sec. 6 only | not used by Sec. 8; ACLGFals: no failure on λ² ≤ 47 |

## 8. Reviewer (one run, decisive arguments with TeX quoted inline; brief logs/reviewer_brief.md, which cites l. 1200 for the line at l. 1198)
G1: GAP, ROUTINE (directness not established either; re-prove η, β on V ⊖ V_low). F1: CORRECT, ROUTINE (Galerkin interval counts
diagnostic, not certified). F2: GAP, SIGNIFICANT (missing a rigorous bound N_λ(≤ θ) ≤ 2λ² + o(λ²) at the moving threshold; the
global-support mechanism alone does not prove loss of near-radicality). G3: GAP, ROUTINE (keep the cutoff δ < 1 − κ_{n0}(σ); spacing
needed down to λ^{−A}). Prop. reduction: survives with orthogonal blocks; "vacuous" not established. Hence the first line reads GAP, not
FALSE STEP. The λ² = 2 computation (sec. 3) was done after the review.

## 9. Novelty and prior work
arXiv listings math.NT and math.FA 2026-11 (fetched 12 Nov): only 2611.13129 itself; no comment or reply. Our files: ACLGFals, ACLGH2
(hypotheses and sizes; ACLGH2 wrote "Prop. reduction's bookkeeping assumes polynomial dimensions", refined in sec. 4: only the per-vector
route to (R2) does), R419/R432/Harvest6 list the hypotheses only. Candidates for the orchestrator: REFUTATIONS (gap in a claimed
conditional proof: F2 + G1) or NOGO (N13: "(R2) at scale λσ contradicts the minuscule count; the transition block is ≥ 2λ²(σ² − 1)").
The Landau-type deficit bound of sec. 6 would turn F2 into a theorem.

## 10. Verdict in plain words
The logic of Andrews–Connes–La Gatta's reduction is almost right: each lemma checks, and the final implication holds once one block is
made orthogonal to the other. What fails is the bookkeeping that makes it look like a reduction. Their step (R2) asks that a whole family
of conditioning vectors be almost invisible to the Weil form, but that family is about σ² times larger than the number of almost-invisible
directions the Weil form has (their own count, about 2λ²). In the smallest case, where everything is known, only one of its eleven
directions is almost invisible; the other ten have energies from 0.07 to 3. Repaired, the small leftover block that should decide RH
has at least about 2λ²(σ² − 1) dimensions, not a logarithmic number, and quasi-polynomially many if the semilocal angle count holds. So
the route does not shrink RH to a small finite check. Turning the contradiction with their count into a proof needs more than RH:
under RH the number of almost-invisible directions is at least 2λ² plus twice the number of repeated zeros up to height 2πλ² (proved
here). So capping it would prove that almost all zeros are simple, which is not known.

## 11. Follow-up (second session): can RH give the count that F2 needs?
Question: does RH imply N_λ(≤ θ_λ) ≤ 2λ² + o(λ²) (θ_λ = Cτ(log λ)^C, τ = λ^{−A})? If so, F2 becomes a contradiction: ACLG's hypotheses imply RH,
and (R2) forces N_λ ≥ dim V_low ≈ 2λ²σ². Answer: not as things stand. The count first requires E(T) = O(T), beyond known results, and a sampling estimate (an open
obstacle, not an impossibility proof). Proofs: DERIVATION.md sec. 11.
- Prop. L / Corollary L (RH; proved). Inputs: the explicit formula QW_λ(f,f) = Σ_γ m_γ|Mf(iγ)|² on C_c^∞((1/λ,λ)); the trace bound
  #{μ > ½} ≥ 2ℓ_1T/π − (4/π²)log(2ℓ_1T) − 6/π² for time-band concentration; Riemann–von Mangoldt with S(T) = O(log T); a band-limited
  plateau Φ ∈ PW_δ. Functions GΦ, with G concentrated on |t| ≤ T and vanishing at the distinct ordinates |γ| ≤ T + 2H, give
     N_λ(≤ θ_λ) ≥ 2λ² + 2E(2πλ²) − O(λ^{3/2}),   E(T) = N(T) − N_d(T) (multiplicity excess),
  at T = 2πλ², H = λ, δ = λ^{−1/2}. By-product: under RH the lower half of Fact count is a theorem, N_λ(≤ θ_λ) ≥ 2λ² − O(λ^{3/2}).
- Where the requested upper bound fails. Any bound N_λ(≤ θ_λ) ≤ 2λ² + B_λ implies E(2πλ²) ≤ ½B_λ + O(λ^{3/2}). The bound F2 needs at fixed σ,
  (2 − ε)λ²σ², would imply E(T) = O(T): all but O(T) of the ≈ (T/2π)log T zeros up to T simple. Under RH the best known is
  E(T) ≤ (1/6 + o(1))N(T) (Montgomery's Σ_{0<γ≤T} m_γ ≤ (4/3 + o(1))N(T); improvements to ≈ 70% simple keep E ≍ N), a factor log T short.
  Close pairs count the same way: zeros within ρ ≈ λ^{−A/2−1}(log λ)^{−2} of another add to the lower bound, and pair correlation does not
  control them. An upper bound would also need a sampling estimate above T_0 = 2πλ², where RH gap bounds ((π + o(1))/log log t) exceed
  the Nyquist spacing π/log λ; no such theorem exists. So the obstruction is exact for multiplicities and close pairs; for the sampling
  estimate the method is missing.
- (b) for the semilocal pair. (b1) proved: #{m: δ_m(λσ) ≤ τ} ≥ #{even n: δ_n^{arch}(λσ) ≤ τ²b_S²/(6B_S²)} ≥ (2 − ε)λ²σ². Proof by min–max on
  η_S(deep prolates), with [CCMsl] Prop. 4.1, F_Sη_S = η_SF_∞ and b_S ≤ |ζ_S(½+it)| ≤ B_S, log(B_S/b_S) = O(λ/log λ). (b2), injectivity of
  V_low's restriction map, is proved only for the archimedean pair (moment argument); for S ≠ {∞} it needs analyticity of the semilocal pair
  vectors, unknown.
- Theorem R (archimedean pair; semilocal pair under (b2)): (R2) on V_low ⇒ N_λ(≤ θ_λ) ≥ (2 − ε)λ²σ² for λ ≥ λ_0(ε). ACLG's hypotheses
  therefore force a near-null space σ²(1 − ε) times the deficit count 2λ² (l. 1390, Fact count). They are inconsistent exactly when the
  near-null count stays below that, which is a theorem about the zeros (E(T) = O(T), no clusters at scale λ^{−A/2−1}, sampling above T_0),
  open under RH. F2 stays GAP; it is not upgraded to INCONSISTENT HYPOTHESES.
- Reviewer (second run, decisive argument inline, brief logs/reviewer_brief2.md): Prop. L CORRECT, SIGNIFICANT (regime H ≥ 1, 0 < δ < ℓ,
  ε_1 < 1, k > 3/2); Corollary L CORRECT, ROUTINE (integer k ≥ A + 3); Consequence GAP, ROUTINE (the implication count ⇒ E-bound is correct, with
  Montgomery's sum counted with multiplicity, Σ_distinct m_γ² ≤ (4/3 + o(1))N; but "best known" is a literature claim, not proved from the inputs,
  and a necessary condition shows an open obstacle, not that RH cannot imply the count); (b1) CORRECT, ROUTINE (given CCMsl Prop. 4.1, the
  intertwining and uniform rung asymptotics down to τ_1); status: GAP (open), not INCONSISTENT HYPOTHESES. E(T) = O(T) is necessary for the
  fixed-σ upper count but not sufficient; an actual spectral upper-count proof remains to be found.

## Process note
One run broke the 300 MB cap: vlow.py's first version peaked at 1.36 GB for 12 s (/usr/bin/time; a 181 × 360 000 Legendre table built
at once). It was rewritten row by row; every later run stayed ≤ 68 MB (spectra 49–68 MB, prolcount 43 MB, vlow ≤ 60 MB under a 0.5-s
RSS poller with kill at 270 MB). One compute process at a time, single-threaded, nice 19 for the spectra loop.

RH relevance: N13 (a conditional route to RH)

## Files
- DERIVATION.md (chain, G1–G3, F1, F2 with proofs, sizes, circularity, (R2*)/(S*), computations, reviewer record), this REPORT.md.
- code/: vlow.py (QW_λ on V_low), counts2.py (eigenvalue counts), prolcount.py (rung counts), run_spec.sh (full Galerkin spectra via
  ../ACLGFals/code/spec_small.py), mon.sh (RSS-capped runner).
- logs/: vlow_*.json/.time, spec_all.jsonl, spec_c*.time, counts_table.md, prolcount*.jsonl/.time, reviewer_brief.md, reviewer_brief2.md.

Orchestrator note (12 Nov 2026): an independent reviewer of §11 (DERIVATION.md): Prop. L and Corollary L CORRECT but ROUTINE (standard concentration plus finite rank; state 0 < δ < ℓ, 2ℓ_1T ≥ 1, integer k ≥ A + 3); (ii) correct; two wordings overclaim: 'E(T) ≍ N(T)' should read 'known bounds give E(T) ≤ (1/6 + o(1))N(T)' (no lower bound for E is known), and 'F2 cannot be closed from RH by counting' should read 'closing F2 by counting would prove an almost-all-simple theorem' (necessity, not impossibility). Recorded as COROLLARIES (Prop. L) and NOGO/TREE N13 (the (R2) count), not as a discovery. Orchestrator rerun of prolcount.py at t² = 8, 18, 20: 11, 30, 34, identical.
