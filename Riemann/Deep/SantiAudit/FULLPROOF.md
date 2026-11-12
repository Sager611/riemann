# SantiAudit — full referee notes (Santibáñez-Leal, Zenodo 23134787 v0.03, 17 pp., 4 Oct 2026)

## 1. Theorem chain (exact statements from the PDF)
Notation: N(T,H) zeros T<γ≤T+H (mult.); S(T,H) simple zeros with β=1/2; O(T,H) distinct sign changes of Z(t) in (T,T+H]; H=T^θ, L=log T.
- (1) Wang [14]: c(θ) = 2 − θ/2 − (1/√2)cot(θ/√2), > 0 iff θ > θ0 = 0.550193964744.
- (2) Parity [4, Prop 2.5] (author's own earlier paper): any lower density k for O/N gives S/N ≥ h(θ;k) = (3+k−√((1−k)(9−k−8c(θ))))/4.
  Derivation (ours): from (q−s)(1−o) ≥ 2(1−s)², q ≤ 2−c, o ≥ k, solve the quadratic in u=1−s. h>0 ⟺ k > 1−2/(2−c(θ)) (checked by hand).
- (3)-(5) νA = min{1/2, 2θ−1, (5θ−2)/6} (elementary, Sec. 3); νB = min{1/2, (17/33)(2θ−1)} (Bettin–Chandee [1]); 0<ν<ν*=max(νA,νB).
- Thm 1.1 (D+MV): θ∈[0.527,1): liminf S/N ≥ max{0, c, (c+2κθ)/3, h(θ;κθ)} > 0; detector ν=499/10000, κ ≥ 71589990975836346694543352023356087639/(2·10^39) ≈ 0.0357950; h(0.527;κ) > 0.0005947542001.
- Thm 1.2 (D): ∫ w|Vψ(σ0+it)|² = c(P,Q,R,ν)ŵ(0) + O(H/L), V = Q(−(1/L)d/ds)ζ, ψ = Σ_{h≤y} μ(h)P(log(y/h)/log y) h^{−s−1/2+σ0}, σ0 = 1/2−R/L, y=T^ν, w smooth, =1 on [T,T+H], supp [T−Δ,T+H+Δ], Δ=H/L, w^(j) ≪ Δ^{−j};
  c = 1 + (1/ν)∫∫(w_R(v)P'(u) + ν w_R'(v)P(u))² du dv, w_R(v)=e^{Rv}Q(v).
- Thm 1.3 (D): if Q(x)+Q(1−x) ≡ nonzero const, liminf O/N ≥ κ = 1 − (1/R)log c.
- Sec. 2 (Young's short proof [15] localized): valid for ν < θ−1/2 only (off-diagonal hm≠kn killed by parts: |log(hm/kn)| ≥ 1/(2√(hkmn)), √(hkmn) ≤ y T^{(1+ε)/2}).
- Sec. 3: Gaussian weight e^{−(t−T)²/H²}; exact formula (10) = 2π√(k/h)Σσ_{α,−β}(n)K_β(nk/h) − residue; Estermann FE (12) per h=gq,k=gp; dual Dirichlet polynomial in (pq/n)^{it}, |t| ≲ T/H; Lemma 3.1 mean square; Cauchy–Schwarz → relative errors TM/H² + TM³/H^{5/2}.

## 2. Heuristic cross-check of the claimed range (ours)
- Trivial bound of the Voronoi dual (additive twist e(−pn/q), primal bump n ∈ qT/(2πp) ± qH/(2πp), dual length ≈ pqT/H², Bessel amplitude):
  per (g,p,q) ≈ √(T/H)/g; total ≈ M²√(T/H); relative M²T^{1/2}/H^{3/2} ⟹ ν < (3θ−1)/4 (= Steuding's T^{1/3}M^{4/3} range; at θ=1 gives BCHB's 1/2).
- At θ=0.527: (3θ−1)/4 = 0.14525 > 2θ−1 = 0.054 > ν = 0.0499 > θ−1/2 = 0.027. So the claimed range is inside the classical trivial-dual range: the statement is not over-reaching; only the proof route is new.

## 3. Referee of the chain (line by line; ✓ = checked by us)
### 3a. Detector constants (MV) — recomputed exactly (arb_check.py, python-flint, 4000 bits; log arb_check.log)
- Data: CAOS source zip (Zenodo 23135354) EXP-010 frozen-parameters.json, set ν=499/10000, R=19321/1000, r=7709/10000, 100 rationals x_j.
- Exact: deg Q = 201, Q(0)=1, Q(1)=0, Q(y)+Q(1−y) ≡ 1; P(0)=0, P(1)=1, ∫PP' = 1/2. c = 1/2 + (C_P/ν)∫e^{2Rv}Q² + νB_P∫e^{2Rv}(RQ+Q')² (✓ from ∫ww' = −1/2).
- c = 123210507.3277167487001823621…; κ = 0.0357949954879181733472716760116780438195 (paper's rational lower bound 0.035794995487918173347271676011678043819500: ✓, margin 3e-41).
- c(0.527) = −0.0723758394637755709…; onset k* = 1−2/(2−c) = 0.03492408958…; h(0.527;κ) = 0.000594754200102962… > 0.0005947542001 ✓ (margin 3e-15).
- ν = 0.0499 < νA(0.527) = 2θ−1 = 27/500 ✓; but ν > θ−1/2 = 0.027: at θ=0.527 the claim needs Sec. 3 (Sec. 2 alone is insufficient).
### 3b. Sec. 4 (Levinson counting, Thm 1.3) ✓
- Lemma 4.1 (χ(s)ζ^{(k)}(1−s) = (−1)^k(D+λ)^kζ, E_β, βZ = 2Re(ωṼ)) re-derived ✓. (21)-(22) ✓ (λ−L = O(1) on t ≍ T).
- Lemma 4.2: argument principle on [1/2,σ1]×[T,T+H] for F = Ṽ/∏(s−1/2−it_j)^{k_j}; Δφ = Δϑ + Δ_brt arg Ṽ − πK − 2πN_> (sign checked); level-crossing count O ≥ |Δφ|/π − 1 − J ✓; horizontal sides O(L) (Backlund–Jensen) ✓.
- Thm 1.3: Littlewood on [σ0,σ1]; horizontal O(L); right edge O(1)+O(H/L); concavity + |a+b|² ≤ (1+L^{−1/2})|a|²+(1+L^{1/2})|b|² + Thm 1.2 + (22) ✓; N_Ṽ ≤ (HL/4πR)log c + o(HL); O/N ≥ 1 − (1/R)log c ✓. Short window costs only O(L²) in counts vs N ≍ HL: harmless.
### 3c. Sec. 5 (transfer to simple zeros) ✓
- [4, Prop 2.5] (Zenodo 22860012 v0.07): (Q−S)(N−O) ≥ 2(N−S)² for every conjugation-invariant finite multiset, Q = ‖A‖²_HS = Σ m_z m_w K(z−w)², A = Σ m_z|v_z⟩⟨v_z̄|, v_z = η e^{−2πizu}, K = (η²)^.
  Re-derived: A self-adjoint, tr A = N; Bessel Q ≥ Σα_j²; Gram–Schmidt ranges U (dim ≤ d = r+k), M (dim ≤ S, Σα ≤ S), L (α ≤ 0: e ⊥ v_x, pair term 2m(|g⟩⟨g|−|h⟩⟨h|)); x² ≥ 2tx−t², 2x−1, 2tx ⟹ Q ≥ 2tN−(2t−1)S−t²d (t ≥ 1); t=(N−S)/d ≥ 2; N−O ≥ 2d ✓. Off-line zeros allowed: no RH.
- Same Q as Wang's S_K(I) (Wang (both-counts)); Wang Thm 2.2 + (eq:unweighted): S_K/N → C(f) for f=η², supp η ⊂ (−λ/2,λ/2), every 0<λ<θ<1 ✓ (range covers θ=0.527).
- R(A,b) smaller root ✓; R(2−c,1−k) = h ✓; monotonicity (c' = ½cot²(θ/√2), ∂_k) ✓; branches c and (c+2k)/3 re-derived (per-ordinate m² ≥ 2m − 1_{m=1}; m² ≥ 2m + 2·1_odd − 3·1_{m=1}) ✓.
### 3d. Sec. 3 (extended moment, ν < 2θ−1 and (5θ−2)/6) — re-derived ✓
- (10) Gaussian twisted moment = 2π√(k/h)Σσ_{α,−β}(n)K_β(nk/h) − 2πζ(α+β)G(1−α)(h/k)^{1/2−α}: re-derived (ζ(1/2+β−it) = χ(1−w)ζ(w), χ(1−w) = Mellin of 2cos(2πx), x = (kn/h)e^{2y/H}); residue term ≪ e^{−T²/H²} (negligible) ✓.
- K_β on x ≍ T: ε=+1 branch ≪ e^{−(T+2πx)²/H²}; ε=−1 branch = e(−x)·Gaussian bump at x = T/2π, width H; primal n ∈ qT/(2πp) ± O(qH/p), additive phase e(−pn/q) ✓; Taylor remainder (T/H²)^J ✓.
- (12) = classical Estermann FE with σ_{a,b}(n) = n^{−b}σ_{b−a}(n) (bars on Ā lost in pdftotext; "inverse additive phase" in text) ✓.
- (14) prefactor 4π/(gpq)·(q/p)^{−z/2}n^{−z/2}, Γ(u−d/2)Γ(u+d/2)cos(πu)(pq/4π²n)^u Ã(1−u−z/2): re-derived from (12) with F̃(s) = (q/p)^sÃ(s), u = 1−s−z/2 ✓.
- (15) L¹ multiplier HT^{−σ}(T/H)^{2σ} ✓; Lemma 3.1 (spacing δ = 1/(16PQN), collision R = 2N·maxτ) ✓; (17) by Cauchy–Schwarz on |t| ≤ K0 = T/H ✓; ‖c‖₂ ≈ (PQ)^{σ−1/2}N^{1/2−σ}, R^{1/2} ≈ N^{1/2} ✓.
- σ = 1/2, N ≤ N0 = TPQ/H²: (H/g)a0^{1/2}N^{1/2}(1+PQN/K0)^{1/2}; N > N0: σ = 1/2+D gains exactly (N0/N)^D (σ-dependence is H(N0/N)^σ) ✓; Σ over P,Q ≤ M/g, g ⟹ (18) H T^ε(TM/H² + TM³/H^{5/2}) ✓.
- Exponents at θ−η = 52699/100000, ν = 499/10000: 1−2θ+ν = −51/12500, 1−5θ/2+3ν = −6711/40000 (recomputed exactly) ✓. p=1 or q=1: N0 ≤ T^{1+ν−2θ} < 1 ✓. K-Bessel branch O(H(M²/T)^σ) ✓.
- Main term = Estermann residues (H√π/g)[ζ(1+α+β)q^{−1−α}p^{−1−β} + (T/2π)^{−α−β}ζ(1−α−β)q^{−1+β}p^{−1+α}] = Young's diagonal term per (h,k); Young's Lemma 7 evaluation needs only ν < 1/2 ✓.
- Compact windows (3.7): w_D = φ_{s0} * Σ_{j<D}(−s0²/4)^j w^{(2j)}/j!, ŵ_D = ŵ e^{−u}Σ_{j<D}u^j/j! (u = s0²ξ²/4), ‖w−w_D‖∞ ≪ (s0/Δ)^{2D} = T^{−2Dη}L^{2D}; ∫|mixture| = O(H) ✓; costs η in θ (needs ν < ν*(θ−η)) ✓.
### 3e. Onsets (onsets.py/onsets.log; exact κ of the 8 frozen detectors)
| moment range used | least θ proved (best frozen detector) | optimal-shape heuristic κ = 0.71733ν |
|---|---|---|
| Sec. 2 only, ν < θ−1/2 (Young localized) | 0.5339856 (ν = 0.0339) | 0.53399 |
| paper's νA (Sec. 3) | 0.5264286 (ν = 0.0499) — paper states 0.527 ✓ | 0.52573 |
| trivial dual / Steuding range ν < (3θ−1)/4 (not proved here) | 0.5031740 (ν = 0.0999) | 0.4938 → every θ > 1/2 |

## 4. Falsification (smallest failing parameter per lemma; adversarial tests)
- Lemma 2.1 (Young localized, off-diagonal killed by parts): needs ν < θ−1/2 (√(hkmn) ≤ yT^{(1+ε)/2} < Δ); at θ = 0.527 it allows only ν < 0.027 < 0.0499 — Sec. 2 alone stops at θ = 0.5339856 (ν = 0.0339 detector). Not misused: the paper invokes Sec. 3 there.
- (18) with ν = 0.0499: fails for θ ≤ (1+ν)/2 = 0.52495 (TM/H² = T^{1+ν−2θ}); (5θ−2)/6 > ν needs θ > 0.45988. Binding bottom constraint is the parity onset, not (18).
- Detector: h(θ; 0.0357950) > 0 ⟺ θ > 0.5264286 (onsets.log). Paper's 0.527 sits 0.0006 above: consistent.
- Wang's input: c(θ) < 0 below θ0 = 0.5501939647, but the pair-sum asymptotic S_K(I) = (C(f)+o(1))HL/2π is proved for every 0<λ<θ<1 (Wang Lemma (eq:unweighted), Prop (F-pointwise): errors O(H + xL³ + …), x ≤ T^λ). Used inside its range: Thm 5.1 needs the q-bound, not c > 0.
- Parity product (37) and Lamzouri S ≥ 2N−Q and 3S ≥ 2N−Q+2O: 120,000 random conjugation-invariant multisets (real points mult 1–3, conjugate pairs heights 1e−7…3, mult 1–2), kernels K = (η²)^ for indicator and optimal cosine densities on (−λ/2,λ/2), λ = 0.527: min normalized slack 0 (−4e−16 rounding; equality at orthogonal configurations); targeted {x0, z, z̄} family: min slack +1.4e−5 → 0⁺ as Im z → 0 (double-point limit). No violation (parity_fals.py/.log, 19 MB).
- Exact formula (10) numerically (check10.py/.log, mpmath, T=40, H=6, α=0.03+0.02i, β=0.05−0.01i): (h,k)=(1,1): rel. diff 1.4e−11; (2,1): 5.0e−10; (1,2): 2.0e−12. ✓
- Hidden hypotheses: none. RH not used (Lamzouri/parity allow non-real points; Wang via BGST unconditional explicit formula); no pair-correlation conjecture (Wang's localized Montgomery formula is proved for α < θ); no simplicity assumption (O counts distinct sign changes; S is derived).
- Minor: Thm 5.1 writes η ∈ C_c, Wang's lemma needs C_c^∞ (density argument, harmless); pdftotext drops the bars on Ā in (12) (text says "inverse additive phase").

## 5. Reviewer (most delicate step: the dual bound (17)→(18), and our sharper bound)
- Item A (paper's (18) from the block bound, exponents −51/12500, −6711/40000): CORRECT.
- Item B (absolute x-space bound of the same dual via the Hankel form: |D(g,p,q,n)| ≪ (H/g)(pq)^{−3/4}(nT)^{−1/4}n^ε(1+√(n/N0))^{−A} ⟹ total ≪ T^ε M²(T/H)^{1/2}, relative T^{1/2+ε}M²/H^{3/2}, admissible iff ν < (3θ−1)/4): CORRECT; (3θ−1)/4 > min(2θ−1,(5θ−2)/6) on (1/2,1). Pointwise the B majorant is not always smaller (needs M < (T/H)^{1/2}); the admissible range is larger everywhere.

## 6. Lead E1 (ours, beyond the audited claim): every fixed θ > 1/2
Claim. For every fixed θ ∈ (1/2,1): liminf S(T,T^θ)/N(T,T^θ) ≥ h(θ;κ1) ≥ h(1/2;κ1) = 0.0055353 > 0, κ1 = 0.0844779079697484032032 (P,Q of the frozen ν = 999/10000 detector re-used at ν = 3/25, R = 177/20; nu12.log, exact rationals + Arb 4000 bits; k*(1/2) = 0.0766235550685757017).
Proof sketch (only §3.5 changes; everything else as audited):
1. Exact transformation (10)–(12), (14) as in the paper ((10) checked numerically to 1e−10).
2. Dual term in x-space (paper §3.6 Hankel symbols): D(g,p,q,n) = (2π b_{gq}b_{gp}/(gq))·(c/p)∫A(ξ)𝒥_d(X(ξ))dξ·σ(n)e(∓p̄n/q), X = 4π√(nξ/(pq)) ≍ λ = √(nT/(pq)) ≥ T^{(1−2ν)/2}; 𝒥_d = √(2/(πX))Σ_± e^{±iX}a_±(X,d), X^k∂^k a_± = O_k(1).
3. A(ξ) = ξ^{−β+j}H^{−2j}Φ_j((T−2πξ)/H), Φ_j Schwartz (Fourier transform of the Gaussian × (e^{2y/H}−1−2y/H)^j): non-oscillating bump of width H at T/2π.
4. X' ≍ λ/T, X^{(k)} ≪ X/ξ^k; K integrations by parts: |∫A𝒥_d| ≪ H λ^{−1/2}(1+Hλ/T)^{−K} = H(pq/(nT))^{1/4}(1+√(n/N0))^{−K}, N0 = pqT/H².
5. Σ_n d(n)n^{−1/4}(1+√(n/N0))^{−K} ≪ T^ε N0^{3/4} (also if N0<1) ⟹ per (g,p,q): ≪ T^ε g^{−1}(T/H)^{1/2}; Σ over p,q ≤ M/g, g: ≪ T^ε M²(T/H)^{1/2}.
6. Relative error T^{1/2+2ν−3θ/2+ε}: power-small iff ν < (3θ−1)/4 (Steuding's classical line, now with shifts, any fixed Q). Other constraints unchanged: θ > 1/2 (Taylor (T/H²)^J), ν < 1/2 (K-branch, Young's main term), compact windows at θ−η.
7. ν = 3/25 < (3θ−1)/4 for all θ > 37/75 = 0.4933; κ1 > k*(θ) for all θ ≥ 1/2 (k* decreasing) ⟹ h(θ;κ1) > 0 on (1/2,1) via Thm 1.3 + Thm 5.1 (unchanged).
Values h(θ;κ1): θ=0.5: 0.0055353; 0.51: 0.016352; 0.527: 0.034045; 0.55: 0.056689; 0.6: 0.10138 (lead_h.log). With the frozen ν = 0.0999 detector alone: θ > 0.5031740.
Status: step 5 bound reviewed CORRECT (reviewer, Item B). Not yet: full write-up of steps 2–4 uniformly in shifts/profiles, independent second check, AGENTS.md novelty gate.
Novelty pointers: arXiv math.NT listings 2026-10 (854 titles) and 2026-11 (280) show no short-interval simple-critical result; Zenodo concept 22984154 latest = v0.03; Wang 2609.07918v2 (7 Oct) states only θ > θ0; Steuding (2002) had the range (3θ−1)/4 with linear Q (θ ≥ 0.552).
