# DERIVATION — the z = 0 clause of Liu's M2 (Thm 12.16) as used in Remark 12.18

Sources: Liu 2610.26835v2 e-print number_fields_v4.tex ("l."); OpenAI 7/8 paper.tex at openai/math 1581d4cb ("o.l."); our Lean port Riemann/Deep/QRH/OAI/NumberTheory/DirichletL ("Lean"). Script: exc_orders.py (log: logs/exc_orders.log).

## 0. Claim and its M2 coordinates
- Remark 12.18 (l. 7836–7848): F ⊃ μ12 cyclotomic, ψ_k(n) = ν(n)χ_n(k), q_k ≤ U, ψ_k nonprincipal, X = U^{1/2}: Σ_k |X^{-1/2} Σ_n ψ_k(n) W(q_n/X)|^4 ≪ U^{1+ε}.
- In M2 coordinates: Z = U, m = 1, q = 0, M = m + q = 1, n1 = n2 = 1/2, A = n1 + n2 = M. Since A > 5M/6 this is a *centered* input.
- The proof is a transcription of OpenAI §18 (o.l. 12477–14984) to F: same threshold A ≤ 5M/6 (o.l. 12594–12599), L = M/4 (o.l. 12941–12944), v_p(h') ≡ 4 mod 6 (o.l. 14356–14357), maximum A − M at v = L (o.l. 14765–14767).

## 1. Path of the z = 0 induction
- Rank R = 2{2(N+1)t + 2b(M) + c} + p (l. 7622–7640): t = 0 (no prime sums) precedes t = 1, so a plain input only calls plain inputs; no amplifier (w = w_o = ℓ = 0, l. 6774–6776), no prime sum, no zero-free hypothesis. OpenAI's paper has the same order ("all zero-slot bands completed before the positive-slot bands", o.l. 12592). The Lean port does not: Energy/NaturalZeroStage.lean:26 actual_zero_stage takes PositiveAt at the child width and hbeta.
- Reflection (Lemma 12.6, l. 6149): a plain of length n > M/2 is replaced by length M − n + ξ.
- Centering (l. 7650–7660, o.l. 12941–12971): P = S(X1)S(X2), P_comp = S(Y1)S(Y2) with Y1 = Z^{L}, Y2 = Z^{A−L}, Y1Y2 = X1X2, L = M/4; |P|² ≤ 2|P − P_comp|² + 2|P_comp|². After reflecting its longer factor P_comp has total length L + M − (A − L) < 2M/3 + ξ < 5M/6: uncentered, same width, smaller rank. P − P_comp goes through the transforms.

## 2. The two transforms (plain case)
- T1, Poisson in the rows k (l. 6576–6750). Expanded square over whole products U = l1l2, V. Origin: O(Z^{m+ε}) (UV powerful). Nonzero h: first bridge (l. 6628)
  H/(X q_e √q_r) Σ_{h≠0} G_ξ(r,h) Σ_{(a,b)=1} A_C(a)Ā_D(b)(q_aq_b)^{-1/2} τ_C(a)τ̄_D(b) G(a,h)Ḡ(b,−h) K1(Hq_h/(q_eq_rq_aq_b)),
  dual length K0 = 2A − c − d + R + E − m. Cauchy–Schwarz in h; target N_C ≪ Z^{A−c+q̃+B_c−s0+ε}; outer exponent M + ½{B_c + B_d − (c + d − 2p − R)} ≤ M (brace ≤ 0, l. 6741–6750; LiuAudit exponents_m2).
- Positive enlargement: h-range K → K + g, g = J_+ + σ, J = d − c + (K0 − K) (plain case). This is the only lossy step; it is what lowers the width.
- T2, Poisson in h (l. 6830–6981). Second bridge (l. 6850): Z^{K+g−a0} Σ_{u,v} B(u)B̄(v)τ1(u)τ̄1(v)(q_uq_v)^{-1/2} Σ_j F(u,v;j) K2(Z^{K+g}q_j/(q_uq_v)).
  The Gauss sums cancel: F(D2a, E2b; j) = F(D2,E2;j) χ_a(j) χ̄_b(−j) on the coprime locus (Lemma 12.10; LiuAudit corr_q12, 165/165 over Q(ζ12)). So the child is Σ_j |Σ_a B(D2a)τ1(a)χ_a(j) q_a^{-1/2}|², again a plain product moment of the same convolution coefficients, rows j = G_cV_id h'.
  Width (l. 6922–6927): m' = 2a0 − K − g − g2 − V, q' = q̃ + t2 + V, M' = M + J − g − g2 + t2 ≤ M − σ (t2 ≤ g2).
  j = 0: excess A − M − D0 − B_c + g ≤ A − M + σ + δ1 (plain case; 5σ/3 with amplifier, l. 6870–6880). Child allowance Δ' = b2 − p2 + B_c ≥ 0 (l. 6966).
- Model check (one line): with c = d = R = E = 0, a0 = A, K = 2A − m: normalised moment Ê(m) ≤ Ê(m − g) + terminal, i.e. T2∘enlarge∘T1 maps width m to width m − g with the same lengths; reflection rebalances; depth ≤ 4M/σ; at width ≤ ρ absolute counting.

## 3. Exceptional rows (l. 7146–7327) — the order-dependent step
- Identification (l. 7157–7174): at an i = 1 common prime of the nonunit partition G_cV_id ∋ p²; a row inducing a character of the fixed S-supported group needs 2 + v_p(h') ≡ 0 mod 6, so v_p(h') ≥ 4 (Kummer: the inducing character is unramified outside S iff all exterior valuations lie in prescribed classes mod 6).
- Count (l. 7176–7192): h' = (fixed sixth-power-free part, norm ≥ Z^{2f}) × (sixth power) × (finitely many S-unit classes mod U_S^6): #h' ≪ Z^{(m'−2f)/6+ε0} ≤ Z^{(m'−f)/6+δ2/6+ε0}, f = 2v1.
- Raw paired volume per row: Z^{a0−b2+2θ_N+ε0} (each side Z^{α_j/2}, α_j = a0 − c2, a0 − d2).
- Excess = count + volume − (M' + Δ'). exc_orders.py (1) derives it symbolically for any order n (count exponent 1/n):
  E_x = A − (1 − 1/n)M − F1^{(n)} − F2^{(n)},
  F1^{(n)} = c/n + (1−1/n)D0 + 2w/n + q̃/n + w_o + B_c − (1−1/n)g + ℓ, F2^{(n)} = 2b2 − (1−1/n)g2 − p2 + t2 + V/n + f/n.
  At n = 6 this is exactly the printed excess and F1, F2 at l. 7203–7214 (check: True).
- Lower bounds (exc_orders.py (2)): F2 ≥ 2b2/3, equality at the i = 1 nonunit line (7−5)/6 + 1/3 (the +1/3 is the v_p(h') ≥ 4 sparsity); F1 ≥ (2/3)(c + w) − 3σ, the 2/3 attained in the J < 0 branch.
- Uncentered (A ≤ 5M/6): E_x ≤ 3σ + O(ξ + ε0).
- Centered (l. 7257–7327): the common Mellin mean (Lemma 12.3) cancels the leading products because X1'X2' = Y1'Y2'; saving (r − r̃1)_+, r = (L − v)_+, v = c + w + min(c2,d2). With F1 + F2 ≥ 2v/3 − 3σ:
  E_x ≤ A − 5M/6 + sup_v{−2v/3 − (L − v)_+} + 3σ = A − 5M/6 − 2L/3 + 3σ = (A − M) + 3σ at L = M/4.
  Tight: the sup is attained at v = L, and −2L/3 = −M/6 exactly cancels the sextic excess A − 5M/6 at A = M. Feasible window: L ≥ M/4 (saving), L < 5M/12 (comparison must land in the uncentered range after reflection); exc_orders.py (3): best slack M/9, paper's choice slack 0.
- Hence the whole Lindelöf strength of Remark 12.18 sits in one balance: sixth-power count exponent 1/6 against centering saving 2L/3. All four exceptional ingredients (Kummer identification, sixth-power count, raw volume, Mellin mean) are standard and were rechecked here at the level of exponents.

## 4. Order-n bookkeeping (consistency check, not a proof)
- Without centering the recursion gives E_x ≤ A − (1 − 1/n)M ≤ M/n at A = M, i.e. Σ ≪ U^{1+1/n}: 7/6 (n = 6), 4/3 (n = 3). These are exactly the balanced exponents of de Faveri's large sieve for n-th order characters (arXiv 2610.04045 Thm 1.1: (AB)^ε(A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n})), which Xu–Sang (2610.17675) show sharp for sextic squarefree families. So uncentered M2 = optimal general-coefficient large sieve; centering is the one step that uses the smooth convolution structure, and Remark 12.18 cannot follow from any large sieve with arbitrary coefficients.
- Cubic transfer: needed saving M/3; comparison constraint L < M/3; F2^{(3)} ≥ b2 (v_p(h') = 1 allowed, f = v1); F1^{(3)} ≥ (11/15)c in the J < 0 branch if B_c is unchanged. Best saving min(κ3,1)L < (11/45)M: deficit ≥ 4M/45; even κ3 = 1 closes only at the boundary L = M/3. Agrees with Remark 12.18 ("not decided"). Not rederived for n = 3: B_c and the T1 brace, the local Gauss tables (i ∈ {1, 3, 4}), the dual-origin count.

## 5. F-specific inputs on the z = 0 path
| input | l. | Q(√−3) model | status |
|---|---|---|---|
| positive row lift on G/Γ (Lemma 12.1) | 5776 | rows = elements z ∈ O, finite units; Lean Radial (Moments/RadialEligibleEnergy.lean:17) | prose, standard |
| orbit Poisson with coprimality (Prop 8.1) | 3494 | lattice Poisson | LiuAudit ✓ |
| compact dual support (Lemma 9.1) | 3614 | same | LiuAudit ✓ |
| allowed generators, reciprocity (Lemma 3.4) | 1908 | primary elements | LiuAudit ✓ (standard; not tested numerically) |
| Gauss sums, correlation (Lemmas 12.8–12.10) | 6319–6552 | o.l. §18 Gauss identities | LiuAudit ✓n over Q(ζ12) |
| unit twists Υ_S ⊂ Θ, exception invariance | 6726, 7593 | six units | prose |
| Kummer identification, sixth-power count mod U_S^6 | 7157–7192 | o.l. 14337–14369 | rechecked (§3) |
| uniform mean with coprimality (Lemma 12.3) | 5917 | o.l. §18 | LiuAudit ✓ |
| reflection (Lemma 12.6) | 6149 | Dirichlet/Hecke FE | LiuAudit ✓ |
| returned-family measures, enlargements, δ_edge (Lemma 12.13, App. C.4) | 7407, 12741 | Lean width induction (Energy/) | G1: prose only |

## 6. What could still fail
Only the G1 prose: (a) the positive enlargement of the h-rows to a smooth ball when the unit rank is ≥ 1 (the ball is a union of unit orbits on G/Γ, not a lattice box); (b) common Fourier separation of K1, K2 on fixed boxes with polynomial height losses; (c) δ_edge ≤ C*(η + ξ + ε0) for the weighted measures. For F = Q(√−3) all three are the Lean-checked constructions of Energy/ (modulo the hbeta hypothesis, REPORT §1). No step of the z = 0 path was found false.
