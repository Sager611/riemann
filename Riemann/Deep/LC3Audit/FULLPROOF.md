# FULLPROOF — referee proofs for lincomb3 (DERIVATION §3–§5, §8, §9)

Referee: lc3audit, 12 Nov 2026. Sources read in TeX: Lamzouri arXiv 2610.26841v2 (IsolatedZeros2.tex, "Lam l.n"),
arXiv 2610.08965 (knausgard839.tex, "K l.n"). DERIVATION = Riemann/Deep/LinComb3/DERIVATION.md. N_T := (T/2π) log T.

## 0. Inputs (exact form used) and what is load-bearing

I1 (Lam l.873–881, Lem:SelbergCLT). χ_1 ≠ χ_2 primitive mod q, same parity: (log|L_1(½+it)|, log|L_2(½+it)|)/√(½ log log T),
 t uniform on [T, 2T], → (X_1, X_2) i.i.d. N(0,1). Cited as a special case of BH Thm B plus Selberg's density bound (l.878–880).
I1N (DERIVATION 3.1). The same for N functions jointly (BH Thm B for N functions). Not read (Duke 1995, no open copy).
I2 (Lam l.898–908, Lem:BHLem9; BH Lemma 7, from which Lamzouri derives it, is cited at l.896 but not stated, and is not used here). δ small fixed: for > (1−δ)T(log T)/M integers ℓ there is A_ℓ ⊂ [y_ℓ, y_{ℓ+1}],
 meas A_ℓ > (1−δ)M/log T, such that for j = 1, 2 and all t ∈ A_ℓ: N(t + (1−√δ)M/log T, χ_j) − N(t, χ_j) = (1−√δ)M/2π + O(√log M/δ²)
 (Eq:StabilityZeros) and (1/M)∫_0^M |log|L_j(½+it+iu/log T)| − log|L_j(½+it)|| du ≪ √log M/δ² (Eq:BH6.10).
I3 (Lam l.856–861). N_F(T) ~ N_T; N(2T, χ) − N(T, χ) = N_T + O(T).
I4 (Lam l.968–975, dyadic Lem:GaPa). r ∈ L¹ real, even, supp ⊂ [−1, 1], Lipschitz at 0, w(u) = 4/(4−u²):
 Σ_{T<γ,γ'≤2T} r̂(i(ρ−ρ') log T/2π) w(ρ−ρ') = (T/2π) log(qT)(r(0) + 2∫_0^1 α r(α) dα + o_{r,q}(1)), all zeros, with multiplicity.
I5 (K l.517–535, l.1081–1104, l.1149–1235). Lemma lem:threshold (linear algebra), (eq:key2), the certificate (eq:certificate2) with
 (x, b) = (8722, 4389)·10⁻⁶ "for all small ε and all finite sets Y of n real points", C(f) < C'_0 = 6639737/5000000.
I6 ((P), DERIVATION 9.1(iv)). Mean square ∫_T^{2T}(S(t+h,χ) − S(t,χ))² dt ≪_χ T log(2 + h log T), 0 < h ≤ 1 (Selberg 1946, Fujii 1974).
 Not read: AMS serves Trans. AMS 196 (1974) only behind a browser check; Selberg's 1946 memoir has no open copy.

R0 (load-bearing reduction, new in this audit). If δ_1 = … = δ_N = δ then 1 − Σ_j min(1 − δ, 1/N) = max(Nδ − (N−1), 0).
 So for equal inputs the capacity term never helps, and for N ∈ {2, 3} with δ = s₀ (1 − s₀ = 0.32646 < 1/3) the headline
 3s₀ − 2 and 2s₀ − 1 follow from the capacity-free inequality count ≥ (1 − ε)N_T − Σ_j |H_j(T)| (§1, P8).
 Consequence: DERIVATION 3.4 (dominance shares 1/N), 4.6 (capacity) and input I1N are NOT load-bearing for the record;
 the record needs only I1 for each pair, I2 for each function, I3, and Lemma M (I4–I6). 3.4/4.6 are still checked below.

## 1. Proposition N (DERIVATION §1–§4), complete proof

P1 (real function). e^{2iθ_j} = ε̄_j and Λ(s, χ) = ε_χ Λ(1−s, χ̄) (Lam l.832–837) give e^{iθ_j}Λ(s, χ_j) = conj(e^{iθ_j}Λ(1 − s̄, χ_j)),
 so Λ_j is real on Re s = ½ and Z_j(t) = e^{i(ϑ(t)+θ_j)} L_j(½+it) is real (|γ(½+it)| > 0). G = Σ b_j Z_j is real-analytic, G ≢ 0 (F ≢ 0:
 a_F(n) = Σ b_j e^{iθ_j} χ_j(n) is not identically 0 by independence of characters). ord_{t₀} G = ord_{½+it₀} F. If G(a)G(b) < 0 then (a, b)
 contains a zero of G of odd order (finitely many zeros on [a, b]; G changes sign exactly at odd-order zeros). ✓ (DERIVATION 1.1–1.3, 2.)
P2 (parameters; repairs DERIVATION 3.3). Lamzouri's I2 is for a pair. Put δ' = δ/(N−1), σ = √δ'. Apply I2 with δ' to (χ_1, χ_j), 2 ≤ j ≤ N;
 G := ∩_j G_j, A_ℓ := ∩_j A_ℓ^{(j)}: #([1, T log T/M] ∖ G) < δT(log T)/M, meas A_ℓ > (1 − δ)M/log T, and for t ∈ A_ℓ, all j ≤ N:
 N(t + (1−σ)M/log T, χ_j) − N(t, χ_j) = (1−σ)M/2π + O_N(√log M/δ²), and (Eq:BH6.10) with O_N(√log M/δ²).
 [DERIVATION 3.3 says "δ/N in place of δ" but keeps windows (1 − √δ)M and starts √δM/log T: the window of I2 is tied to its δ.
 Using σ in both places (below) is the consistent form. Slip, no effect on any constant.]
P3 (separation; DERIVATION 3.2). E := ∪_{i<k} E_ik, E_ik := {t ∈ [T, 2T] : |log|L_i| − log|L_k|| < δ√(log log T)} at ½+it.
 I1 for (χ_i, χ_k) and the continuous map (x, y) ↦ (x − y)/√2: (log|L_i| − log|L_k|)/√(log log T) → N(0, 1); {|x| < δ} has null boundary,
 so meas E_ik = (P(|X| < δ) + o(1))T ≤ (0.7979δ + o(1))T and meas E ≤ C(N,2)(0.7979δ + o(1))T. Uses only the bivariate I1. ✓
P4 (good blocks). ℓ good :⇔ ℓ ∈ G and J_ℓ := (y_ℓ, y_ℓ + σM/log T) ∩ A_ℓ ∩ E^c ≠ ∅. For ℓ ∈ G, meas((y_ℓ, y_ℓ + σM/log T) ∩ A_ℓ) ≥ (σ − δ)M/log T
 ≥ σM/(2 log T) when δ ≤ σ/2 (⇔ δ ≤ 1/(4(N−1))). A non-good ℓ ∈ G has (y_ℓ, y_ℓ + σM/log T) ∩ A_ℓ ⊂ E; these sets are disjoint in ℓ, so
 #(G ∖ good) ≤ 2 meas(E) log T/(σM) ≪_N (δ/σ)T(log T)/M. Bad ℓ in total: ≪_N √δ T(log T)/M. ✓ (Lam l.917–922.)
P5 (stable logs). Good ℓ: pick t_ℓ ∈ J_ℓ. (Eq:BH6.10) at t_ℓ ∈ A_ℓ makes every log|L_j(½+it_ℓ)| finite (finite integral), and Markov at level
 λ_T = (log log T)^{1/3} gives B_ℓ ⊂ [0, M] with e_T := meas([0, M] ∖ B_ℓ) ≤ N·M·O_N(√log M/δ²)/λ_T → 0 (T → ∞; M, δ fixed) and
 |log|L_j(½+it_ℓ+iu/log T)| − log|L_j(½+it_ℓ)|| ≤ λ_T for all u ∈ B_ℓ, j ≤ N. ✓ (Lam l.924–935.)
P6 (dominance; DERIVATION 4.3). t_ℓ ∉ E: j_ℓ := argmax_j log|L_j(½+it_ℓ)| is unique and beats every rival by ≥ δ√(log log T). For u ∈ B_ℓ,
 i ≠ j_ℓ: the gap stays ≥ δ√(log log T) − 2λ_T ≥ (δ/2)√(log log T) (T ≥ T_0(M)). On the line Z_i/Z_{j_ℓ} is real with modulus |L_i/L_{j_ℓ}|, so
 G = b_{j_ℓ}Z_{j_ℓ}(1 + η), η real, |η| ≤ (N−1)(max|b|/min|b|)e^{−(δ/2)√(log log T)} < 1: sign G = sign(b_{j_ℓ}Z_{j_ℓ}) ≠ 0 on t_ℓ + B_ℓ/log T. ✓
P7 (one block; DERIVATION 4.4). I_ℓ := (t_ℓ, t_ℓ + (1−σ)M/log T] ⊂ (y_ℓ, y_{ℓ+1}) ⊂ (T, 2T] (end blocks discarded); disjoint. k_ℓ := number of zeros of
 L_{j_ℓ} with ordinate in I_ℓ (multiplicity) = (1−σ)M/2π + O_N(√log M/δ²) (P2 at t_ℓ). Ordinates z_1 ≤ … ≤ z_{k_ℓ}; τ := (log log T)^{−1/4}.
 An interior z_m (2 ≤ m ≤ k_ℓ − 1) with both gaps > 2πτ/log T is a simple zero on the line (else z_m repeats: multiplicity, or 1 − ρ̄ ≠ ρ).
 Each open gap (z_{m−1}, z_m) ⊂ I_ℓ has u-length > 2πτ > e_T (T large: τ = (log log T)^{−1/4} ≫ e_T ≍ (log log T)^{−1/3}), so meets t_ℓ + B_ℓ/log T; pick one point per such gap.
 Z_{j_ℓ} has no zero inside a gap (no ordinate there), so a constant sign on it, and opposite signs on the two gaps at a simple zero z_m.
 By P6 G has opposite signs at the two chosen points around z_m; by P1 an odd-order zero of F lies between them. Consecutive chosen points
 bound disjoint intervals, so distinct m give distinct zeros. A non-isolated interior zero has a gap ≤ 2πτ/log T to a neighbour in the global
 order of ordinates (I_ℓ is an interval), so lies in H_{j_ℓ}(T) (Lam l.864, with multiplicity: a repeated ordinate has gap 0).
 Hence #{odd zeros of F in I_ℓ} ≥ k_ℓ − 2 − #(H_{j_ℓ}(T) ∩ I_ℓ). ✓ (Lam l.941–957.)
P8 (sum; capacity-free form). With δ = M^{−1/3} (σ ≤ M^{−1/6}, √log M/(δ²M) = √log M·M^{−1/3}) and M ≥ M_0(ε, N), T ≥ T_0(ε, M):
 Σ_{good ℓ}(k_ℓ − 2) ≥ (1 − O_N(√δ))(T log T/M)((1−σ)M/2π − O_N(√log M/δ²) − 2) ≥ (1 − ε/2)N_T, so, the I_ℓ being disjoint,
  #{odd zeros of F, t ∈ (T, 2T]} ≥ (1 − ε/2)N_T − Σ_j Σ_{good ℓ: j_ℓ=j} #(H_j ∩ I_ℓ) ≥ (1 − ε/2)N_T − Σ_j |H_j(T)|.   (P8)

P9 (capacity; DERIVATION 3.4, 4.6; needs I1N, not needed for the record by R0). S_j := {t ∈ [T, 2T] : log|L_j| > max_{i≠j} log|L_i| at ½+it}.
 The normalised vector lies in the open cone {x_j > max_{i≠j} x_i}, boundary ⊂ ∪_{i≠j}{x_i = x_j} (Gaussian-null), so by I1N and portmanteau
 meas S_j = (P(X_j = max_i X_i) + o(1))T = (1/N + o(1))T (exchangeable, ties null). For good ℓ with j_ℓ = j, P6 gives
 t_ℓ + (B_ℓ ∩ [0, (1−σ)M])/log T ⊂ I_ℓ ∩ S_j, of measure ≥ ((1−σ)M − e_T)/log T, disjoint in ℓ. So #{good ℓ: j_ℓ = j} ≤ (1/N + o(1))T log T/((1−σ)M − e_T)
 and Σ_{good ℓ: j_ℓ=j} #(H_j ∩ I_ℓ) ≤ Σ_{good ℓ: j_ℓ=j} k_ℓ ≤ (1/N + o(1))N_T(1 + O_N(√log M/(δ²M))) ≤ (1/N + ε/2)N_T (M ≥ M_0, T ≥ T_0).
 With P8: count ≥ (1 − ε/2)N_T − Σ_j min(|H_j(T)|, (1/N + ε/2)N_T), which implies Proposition N as stated. ✓ Order: ε → M(ε, N), δ = M^{−1/3} → T.

## 2. Theorem N (DERIVATION §5) and Theorem N^R (DERIVATION §8)

N1. For fixed c > 0 and T large, τ < c, so a c-isolated zero of L_j in (T, 2T] is τ-isolated and not in H_j(T):
 |H_j(T)| ≤ N(2T, χ_j) − N(T, χ_j) − B_j(T, c) ≤ (1 − δ_j + o(1))N_T (I3). Insert in Prop N; sum over the half-open dyadic blocks
 (T/2^{k+1}, T/2^k] above T_0 (Σ_k (T/2^{k+1}) log(T/2^{k+1}) = T log T − O(T)); divide by N_F(T) ~ N_T (I3); ε → 0. ✓ Multiplicity bookkeeping
 as DERIVATION 5.4 (H_j, k_ℓ, N count with multiplicity; each odd zero of F counted once): ✓.
R1 (definition). Y*_j(T, κ) := simple zeros ½+iy of L_j on the line, T < y ≤ 2T, with no other zero of L_j ON the line (any multiplicity) at
 |y' − y| < 2πκ/log T. δ^R_j(κ) := liminf_T |Y*_j(T, κ)|/N_T. H^R_j(T) := (zeros of L_j in (T, 2T], with multiplicity) ∖ Y*_j(T, κ), so
 |H^R_j(T)| = N(2T, χ_j) − N(T, χ_j) − |Y*_j(T, κ)| ≤ (1 − δ^R_j(κ) + o(1))N_T (I3). A c-isolated zero lies in Y*_j(T, c), so δ^R_j(c) ≥ δ_j(c).
R2 (one block). Good ℓ, j = j_ℓ, y ∈ Y*_j ∩ I_ℓ with [y − πκ/log T, y + πκ/log T] ⊂ I_ℓ. U_± := the open halves of (y − πκ/log T, y + πκ/log T).
 The real zeros of Z_j are exactly the ordinates of the zeros of L_j on the line (an off-line zero β+iγ, β ≠ ½, is not a zero of t ↦ L_j(½+it));
 U_± contain none (their points are < πκ/log T < 2πκ/log T from y), so Z_j has constant sign on U_− and on U_+, opposite (y simple).
 Each U_± has u-length πκ > e_T for T ≥ T_0(κ, M), so contains a point of t_ℓ + B_ℓ/log T; by P6 and P1, F has an odd-order zero in
 (y − πκ/log T, y + πκ/log T). Two points of Y*_j are ≥ 2πκ/log T apart (each is an on-line zero), so these open intervals are disjoint.
 Points of Y*_j ∩ I_ℓ failing the inclusion lie within πκ/log T of an end of I_ℓ: at most one per end. Hence
 #{odd zeros of F in I_ℓ} ≥ |Y*_j ∩ I_ℓ| − 2 = k_ℓ − #(H^R_j ∩ I_ℓ) − 2. ✓ (This is DERIVATION 8.2; the "κ/2-neighbourhood" there is this one.)
R3. P8, P9 and N1 verbatim with H^R_j: liminf_T N_0^odd(F; T)/N_F(T) ≥ 1 − Σ_j min(1 − δ^R_j(κ), 1/N) (needs I1N) and, capacity-free,
 ≥ Σ_j δ^R_j(κ) − (N − 1) (needs only I1, I2, I3). ✓ Theorem N^R holds. Order: κ fixed first (T_0 depends on κ), then ε, M, δ, T.
R4 (sanity of the mechanism, falsification attempt). The only property of Y* used is: the sign of Z_j flips at y and nowhere else within
 πκ/log T. Off-line zeros of L_j near y are irrelevant to Z_j's sign and also to P6 (P6 uses |L_j| on the line only, through I2). Off-line
 zeros do enter k_ℓ and I2's stability, but only as counted members of H^R_j. No hidden use of on-line-ness of the neighbours: ✓.


## 3. Lemma M (DERIVATION §9), complete proof

Setting (K l.753–778, read): f even probability density on I = [−½, ½]; f_ε = η_ε², η_ε real, even, smooth, supported inside I, ∫f_ε = 1,
 K_ε = f̂_ε (real, even, Schwartz, K_ε(0) = 1), |K_ε − K| ≤ d_ε → 0 (eq:smooth). F_z(u) = η_ε(u)e^{−2πiuz}, ⟨F_z, F_w⟩ = K_ε(w − z̄).
 Z = finite multiset, invariant under z ↦ z̄ with equal multiplicities ν_z; A := Σ_z ν_z F_z⟨F_{z̄}, ·⟩ on L²(I), finite rank, A* = A,
 tr A = Σν_z = |Z| =: N, tr A² = ‖A‖²_HS = Σ ν_zν_w K_ε(z − w)². For real y, P_y := F_y⟨F_y, ·⟩ ⪰ 0, tr P_y = 1, tr(P_uP_v) = K_ε(u−v)² =: g(u−v).
M1 (lem:threshold, K l.517–535; re-proved). P ⪰ 0, Q Hermitian with ≤ p positive eigenvalues, A = P + Q, c ∈ ℝ ⇒
 tr A² ≥ 2c tr A − 2c tr P + tr φ_c(P) − c²p. Proof: Z := Q_−; tr((P−Z)Q_+) = tr(PQ_+) ≥ 0, so tr A² ≥ tr(P−Z)² + tr Q_+²;
 tr Q_+² ≥ 2c tr Q_+ − c²p (x² ≥ 2cx − c² on ≤ p eigenvalues); remains tr(P−Z)² + 2c tr Z ≥ tr φ_c(P), i.e. with X := P − cI,
 tr X_+² − 2tr(XZ) + tr Z² ≥ 0, true as tr(XZ) ≤ tr(X_+Z) (tr(X_−Z) ≥ 0) and tr(X_+ − Z)² ≥ 0. ✓ Pure linear algebra.
M2 ((eq:key2) for ANY such Z; K l.1072–1082). Let Y = the real points of multiplicity 1 (n of them), P := Σ_{y∈Y} P_y, Q := A − P:
 multiple real points w give m_wP_w ⪰ 0 (rank 1, trace m_w ≥ 2); a pair {z, z̄}, z ∉ ℝ, gives ν(vw* + wv*) = (ν/2)((v+w)(v+w)* − (v−w)(v−w)*)
 (v = F_z, w = F_{z̄}): ≤ 1 positive eigenvalue, trace ν(⟨w,v⟩ + ⟨v,w⟩) = 2ν. n_+ is subadditive (Weyl), so Q has ≤ p positive eigenvalues,
 p := #multiple real points + #pairs, and N − n = tr Q ≥ 2p. M1 with c = 2: tr A² ≥ 4(N − n) − 4p + tr φ_2(P) ≥ 2(N − n) + tr φ_2(U_ε),
 U_ε := (K_ε(y_i − y_j)) = V*V, P = VV* (Ve_i = F_{y_i}) have the same nonzero spectrum and φ_2(0) = 0. Nothing about zeros is used. ✓
M3 (certificate; K l.1084–1090, 1149–1235; Lean GramEnergyStatement, anc/lean/Simple673/Spec.lean, proved in Gram/Main.lean from
 BetaBound and TransferStatement, the latter from the native_decide search). For all small ε and EVERY finite set Y of n distinct reals:
 tr φ_2(U_ε) − n ≥ xn − b·spread(Y) − ω(ε)n − C_2, (x, b) = (8722, 4389)·10⁻⁶, ω(ε) = C_1d_ε ≥ 0, C_2 = 7x + Ω.
 Hypotheses, each checked in the TeX/Lean for arbitrary real sets: (a) |U_ij − K(y_i−y_j)| ≤ d_ε (eq:smooth); (b) U_SS ⪯ 2I on
 θ-separated S (lem:maj with w = m_εf_ε ≤ f ≤ g, Λ = 1.78942 < 2·m_ε for small ε; Lean SeparatedBound); K² ≥ β = 0.05106 ≥ x on [−θ, θ];
 U_ε is a PSD Gram matrix of distinct points (Lean: Function.Injective y, G.PosSemidef). x < 1, b ≥ 0. None refers to zeros. ✓
 The Lean headline (K l.2009–2011) is about ζ only: simple_critical_zeros_673 : ∀ ε > 0, ∃ T₀, ∀ T ≥ T₀,
 (1669159/2478195 − ε)·Ncount 0 T ≤ N0simple 0 T, axioms propext, Classical.choice, Quot.sound + the second search's native_decide.
M4 (merging; DERIVATION 9.1(i)–(ii)). Fix κ ∈ (0, 1/16], κ' := √κ. Real points of Z: distinct y_a with multiplicities m_a. Graph: a ~ b iff
 |y_a − y_b| < κ. On ℝ its components are runs of consecutive points (if a < y < b with a, b in one component, the path from a to b
 has a step p < y < p' with p' − p < κ, so y is joined to p or p'). Merge each component C with |C| ≥ 2 and diam C ≤ κ' into ȳ_C := Σ_C m_ay_a/m_C,
 multiplicity m_C := Σ_C m_a ≥ 2. Z' := result. (i) Z' conjugation-invariant, |Z'| = N, tr A' = N. (ii) ȳ_C ∈ [min C, max C]; every other real
 point lies outside [min C, max C] (runs) at distance ≥ κ from C, so ȳ_C is a new distinct point; the multiplicity-1 real points of Z' are
 exactly Y* ⊔ L, Y* := multiplicity-1 points forming singleton components, L := multiplicity-1 points of components with diam > κ'. ✓

M5 (energy perturbation; DERIVATION 9.1(iii)). A' − A = Σ_C D_C, D_C := Σ_{a∈C} m_a(P_{ȳ_C} − P_{y_a}) (off-line terms identical in A, A').
 tr(D_CD_{C'}) = Σ_{a∈C,b∈C'} m_am_b M(s_a, t_b; Δ), s_a := y_a − ȳ_C, t_b := y_b − ȳ_{C'}, Δ := ȳ_C − ȳ_{C'},
 M(s,t;Δ) := g(Δ) − g(Δ−t) − g(Δ+s) + g(Δ+s−t) = −∫_0^s∫_0^t g''(Δ+σ−τ)dτdσ (∂_σ∂_τ g(Δ+σ−τ) = −g'').
 Taylor: g''(Δ+w) = g''(Δ) + g'''(Δ)w + R(w), |R(w)| ≤ ½w² sup_{|v−Δ|≤|w|}|g⁗(v)|. ∫∫1 = st, ∫∫(σ−τ) = (s²t − st²)/2, and
 Σ_a m_as_a = Σ_b m_bt_b = 0 (centroids) kill Σ m_am_b st and Σ m_am_b(s²t − st²)/2. |s_a|, |t_b| ≤ κ' (centroid inside the run), |σ−τ| ≤ 2κ',
 |∫∫R| ≤ ½·ψ(Δ)·|s||t|(2κ')² ≤ 2ψ(Δ)κ'⁴, ψ(Δ) := sup_{|v−Δ|≤2κ'}|g⁗(v)|. As |Δ − (y_a − y_b)| ≤ 2κ', ψ(Δ) ≤ ψ̃(y_a − y_b),
 ψ̃(u) := sup_{|v−u|≤4κ'}|g⁗(v)| ≤ Ψ₁(u) := sup_{|v−u|≤1}|g⁗(v)| (κ ≤ 1/16). Hence ‖A' − A‖²_HS = Σ_{C,C'} tr(D_CD_{C'})
 ≤ 2κ'⁴ Σ_{a,b real} m_am_b Ψ₁(y_a − y_b) = 2κ² Σ_{a,b} m_am_bΨ₁(y_a − y_b). ✓ (κ'⁴ = κ².) Exact identity, no asymptotics.
M6 (pair counts (P) from I6; DERIVATION 9.1(iv)). Zeros of L(s, χ) in the block, ordered pairs with multiplicity (ρ = ρ' included),
 normalised distance d = |γ−γ'| log T/2π. Claim: #{d ≤ R} ≤ C_χ R N_T for all R ≥ 1, T ≥ T_0(χ). Proof. h := 2πR/log T.
 (1 ≤ R ≤ log T/6π, so 3h ≤ 1.) n(γ) := #{ρ' : |γ'−γ| ≤ h} ≤ h⁻¹∫_{γ−2h}^{γ−h}(N(t+3h) − N(t))dt ([γ−h, γ+h] ⊂ (t, t+3h] there);
 Σ_ρ n(γ) ≤ h⁻¹∫_{T−2h}^{2T}(N(t+3h) − N(t))²dt; N(t+3h) − N(t) = (3h/2π)log(qt/2π) + S(t+3h) − S(t) + O(1/t) gives
 ≤ h⁻¹(2T(3R + o(1))² + 2C_χT log(2 + 6πR) + O(T)) ≪_χ R N_T. (R > log T/6π: N(t + D) − N(t) ≪ (D + 1)log t with D = h ≤ 2πR/log T
 gives ≪ N_T(2πR + log T) ≪ R N_T.) ✓ Then, real points only (a sub-count), Ψ₁ Schwartz (g = K_ε² Schwartz):
 Σ_{a,b} m_am_bΨ₁(y_a − y_b) ≤ Σ_{k≥0} sup_{[k,k+1)}Ψ₁ · C_χ(k+1)N_T =: C_χΨ_εN_T, Ψ_ε < ∞, uniform in κ ≤ 1/16.
 So ‖A' − A‖_HS ≤ c_εκ√N_T, c_ε := (2C_χΨ_ε)^{1/2}. (DERIVATION's c_εκ√N.) ✓ [Only I6 here is quoted.]
M7 (long chains; 9.1(vi)). y ∈ C, diam C > κ', consecutive gaps < κ: one side of y carries > κ'/2 of the run, so (y, y ± κ'/2] holds
 n points with (n+1)κ > κ'/2, n ≥ 1/(2√κ) − 1 ≥ 1/(4√κ) (κ ≤ 1/16). Each y ∈ L has ≥ 1/(4√κ) partners at distance ≤ ½:
 |L| ≤ 4√κ·#{ordered pairs, d ≤ 1} ≤ 4√κ·C_χN_T. ✓ (DERIVATION: 8C_χ√κN, with its (1+R) convention.)
M8 (conclusion; 9.2). ‖A'‖_HS ≤ ‖A‖_HS + ‖A' − A‖_HS (A, A' finite-rank self-adjoint on L²(I)) and the energy (M9):
 tr A'² ≤ (√(C(f_ε) + o(1)) + c_εκ)²N. M2 for Z' and M3 for Y' := Y* ⊔ L (n' points, spread(Y') ≤ spread(all real points) ≤ N + o(N),
 b ≥ 0, ω ≥ 0): (√(C(f_ε)+o(1)) + c_εκ)²N ≥ 2(N − n') + (1 + x)n' − bN − ω(ε)N − C_2 − o(N), i.e.
 (1 − x)n' ≥ (2 − (√C(f_ε) + c_εκ)² − b − ω(ε))N − o(N); |Y*| ≥ n' − |L|. Given η: choose ε (C(f_ε) → C(f) < C'_0, ω(ε) → 0), then κ
 (c_εκ and 4C_χ√κ small): liminf |Y*|/N ≥ s₀ − η, s₀ = (2 − C'_0 − b)/(1 − x) = 1669159/2478195 (exact rational, recomputed). ✓
M9 (energy for fixed χ; 9.3). I4 with r = R_ε := f_ε*f_ε and r = R_ε'' (real, even, smooth, supp ⊂ [−1, 1]: admissible) and
 (R_ε − R_ε''/(4L²))^(x) = (1 + π²x²/L²)K_ε(x)² = K_ε(x)²/w(ρ−ρ') at x = i(ρ−ρ')L/2π (L = log T): Σ K_ε(z−z')² = (C(f_ε) + o(1))N_T,
 C(f_ε) = R_ε(0) + 2∫_0^1αR_ε(α)dα = ∫f_ε² + ∬|s−t|f_εf_ε, exactly as K l.791–805 for ζ. Conjugation invariance: Λ(1−ρ̄, χ) = ε_χ conj Λ(ρ, χ) = 0. ✓
 Y* (block-internal isolation) vs R1 (isolation from all on-line zeros): differ by O(log T) zeros near T, 2T. ✓ Lemma M holds.


## 4. Numerical falsification attempts (logs/lemmaM_hs.log; lemmaM_hs.py; peak 189 MB, 91 s, one thread)

Setup: 800 unfolded ζ zeros (t ≈ 18000–18600, Riemann/Deep/ExploreHelio/results/zeta_zeros.npy), two windows: a C^∞ bump η² and
 2610.08965's second window smoothed at the edges (ε = 0.05). Planted per κ: tight pairs, spread clusters (diam ≤ κ', gaps < κ, up to
 1/√κ points), double points, long chains (extent 2κ'). ‖A'−A‖²_HS computed exactly on the Fourier side, Σw_iw_jg(p_i−p_j) =
 ∫R(α)|Σw_ie(−αp_i)|²dα (R = f_ε*f_ε, 4096 composite Gauss nodes), cross-checked in x-space on 25 clusters (agreement ≤ 1.6·10⁻⁸ relative, 12 runs).
F1 (termwise centroid bound M5): max over all cluster pairs of |tr(D_CD_{C'})|/(2ψ(Δ)m_Cm_{C'}κ'⁴) = 6.7·10⁻⁴ … 7.3·10⁻⁴ in all 8 planted
 runs (κ = 1/16, 0.02, 0.005, 0.001; both windows). Never > 1. The value matches the leading term computed by hand:
 tr(D_CD_{C'}) = −¼g⁗(Δ)V_CV_{C'} + O(κ'⁶), V_C := Σ_a m_as_a²; for an evenly filled cluster V_C ≈ m_Cκ'²/12, ratio 1/1152 = 8.7·10⁻⁴.
F2 (sums): E2 ≤ B2 ≤ B2' in every run; ‖A'−A‖_HS/(κ√N) ∈ [0.09, 0.90] against the proved constant (2ΣΨ₁/N)^{1/2} ∈ [17, 99]
 (ratio ≤ 0.0092; both grow as κ falls because the planted pair density per point grows, pairs(d ≤ 1)/N = 6.6 → 58).
 Unplanted zeros at κ = 0.2 and 0.3 (natural clusters: 7 and 15 merged): termwise ratio ≤ 6.6·10⁻⁴, E2 ≤ B2.
F3 (centroid is a refinement, not load-bearing): merging into min C instead inflates E2 by factors 40 to 2.4·10⁴ and obeys the first-order bound
 κ'²Σ_{C,C'}m_Cm_{C'}sup_{|v−Δ|≤2κ'}|g''(v)| in all runs (|M(s,t;Δ)| ≤ |s||t|sup|g''|). Under (P) that bound is ≤ κ C' N, so
 ‖A'−A‖_HS = O(κ^{1/2}√N) → 0 already suffices for M8; centroids only improve the rate to O(κ√N).
F4 (combinatorics M4(ii), M7, triangle): simple points of Z' = Y* ⊔ L exactly and all positions distinct (all runs); |L| ≤ 4√κ·#pairs(d ≤ 1);
 |√tr A'² − √tr A²| ≤ ‖A'−A‖_HS. Unplanted zeros: no component with ≥ 2 points at κ ≤ 1/16 (GUE repulsion), Z' = Z.
F5 (q = 7 numerics of lincomb3, logs only, nothing rerun): the per-block inequality (count ≥ isolated − 2) holds in 63 + 114 + 127 = 304
 strongly dominated blocks; the "isolated" there is isolation from on-line zeros (analyse.py isolated_frac on sign changes of Z_j), i.e. the
 Y* notion of Theorem N^R. Dominance shares for b = (1,1,1): 0.329/0.330/0.341. For b = (1,−2,½), (3,1,−1) the shares are 0.247/0.640/0.113 and
 0.654/0.172/0.174: at t ≤ 2400, log|b_j| is comparable to √(½ log log T) ≈ 1, so 1/N emerges only as log|b|/√(log log T) → 0.
 Not a counterexample (3.4 is asymptotic and, by R0, not load-bearing), but the capacity o(1) is numerically large.


## 5. Findings, load-bearing structure, what survives

G1 (slip, repaired in P2) DERIVATION 3.3 applies BH Lemma 9 "with δ/N in place of δ" but keeps the windows (1 − √δ)M and starts √δM/log T;
 the window of Lemma 9 is tied to its δ. Use σ = √(δ/(N−1)) for both: every bound of §4 holds with O_N constants. No constant changes.
G2 (structure) R0: with equal inputs the capacity term is inert; the record (N = 2, 3) needs neither 3.4/4.6 nor BH Thm B for N functions.
 3.4/4.6 are nevertheless correct given I1N (P9).
G3 (precision) (P) holds uniformly: #{d ≤ R} ≤ C_χRN_T for all R ≥ 1 (M6; for R > log T/6π from N(t+1) − N(t) ≪ log t). DERIVATION's
 "T ≥ T_0(R)" plus "uniform for 1 ≤ R ≤ log T" says the same. Its |L| ≤ 8C√κN is M7's 4C√κN in a (1+R) convention.
G4 (wording) "Lean-checked certificate": the abstract Gram bound is a Lean theorem (GramEnergyStatement, any finite injective real family,
 PSD Gram matrix; proved from BetaBound and TransferStatement, the latter resting on the native_decide axiom of the second search).
 The Lean headline is about ζ; Lemma M itself is a written argument outside Lean (M4–M8 above).
G5 (inputs read only as quoted) I6: AMS serves Fujii (Trans. AMS 196, 1974) only behind a browser check; Selberg 1946 has no open copy.
 BH Thm B / Lemma 9 (Duke 1995) and Garunkštis–Paliulionytė Lemma 3.1 (arXiv API returned nothing) were checked in Lamzouri's statements
 (Lam l.873–908, 968–975), as the brief allows. I1N (BH Thm B for N functions) is not load-bearing (G2).
G6 (numerics) No counterexample. The capacity o(1) is numerically large at t ≤ 2400 (F5); harmless (G2).
No false step found in Proposition N, Theorem N^R or Lemma M.

What each bound needs:
 N = 3: 3s₀ − 2 = 17029/826065 = 0.0206146 (= 51087/2478195): Thm N^R capacity-free (I1 pairwise, I2, I3) + Lemma M (I4, I5, I6).
 N = 2: 2s₀ − 1 = 860123/2478195 = 0.3470764: same inputs.
 Fallback N = 3: 3p − 2 = 0.0174739 (p = 0.672491294): Thm N capacity-free + IsoZero's c-isolated bound for L(s, χ); no Lemma M, no
 2610.08965, no I6 (but Wu's mollified moment and the hinge bound of IsoZero instead; I4 is used by both routes).
Comparison with discoveries/isolated-zeros.md (N8): Lemma M gives ≥ s₀ − η simple zeros on the line isolated from the OTHER ON-LINE zeros at a
 non-explicit κ(η); isolated-zeros gives 0.6724 isolated from ALL zeros at explicit c ≤ 0.00476. Neither implies the other.

