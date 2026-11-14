# FULLPROOF (bothaaudit): proofs behind REPORT.md

Notation: c = λ², L = log λ = (ln c)/2 (Botha's window, his "L"), Q = QW_λ (CCM 2511.22755 sec. 3; = Botha's Q_L, R415), A_λ its operator;
λ_k^even, λ_k^odd its min-max values on even/odd functions. E_N = span of the orthonormal bases of ../CCMFMM/weil_matrix.py
(even φ_j = N_j cos(w_j x), j = 0..N; odd ψ_j = (2/L')^{1/2} sin(w_j x), j = 1..N; L' = ln c, w_j = 2πj/L'); T_N = (Q(b_i, b_j)) in closed form.

## 1. The entry identity (gen.py; checked entrywise against weil_matrix.py, gen.py test: every entry overlaps, mid difference ≤ 2^-597)
weil_matrix.py's off-diagonal entry (i ≠ j) is v = ps P_i P_j + ar − 2 pr with ps = 2 (even), −2 (odd), P_j = pole_vec(j), and
- even: ar = −2 N_i N_j s (w_j S_j − w_i S_i)/(w_i² − w_j²), pr = Σ_k ω_k N_i N_j s (w_j sin(w_j ℓ_k) − w_i sin(w_i ℓ_k))/(w_i² − w_j²);
- odd:  ar = −2 (2/L') s (w_i S_j − w_j S_i)/(w_i² − w_j²),  pr = Σ_k ω_k (2/L') s (w_i sin(w_j ℓ_k) − w_j sin(w_i ℓ_k))/(w_i² − w_j²);
s = (−1)^{i+j}, ℓ_k = ln n_k, ω_k = Λ(n_k) n_k^{-1/2} over prime powers n_k < c. With R_j = S_j + σ_j, σ_j = Σ_k ω_k sin(w_j ℓ_k):
  ar − 2pr = −2 N_i N_j s (w_j R_j − w_i R_i)/(w_i² − w_j²)  (even),   −2 (2/L') s (w_i R_j − w_j R_i)/(w_i² − w_j²)  (odd).
Put x_j = w_j², and a_j = N_j (−1)^j, u_j = −2 w_j R_j (even); a_j = (2/L')^{1/2} (−1)^j w_j, u_j = −2 R_j/w_j (odd). Then in both sectors
  T_ij = ps P_i P_j + a_i a_j (u_j − u_i)/(x_i − x_j)   (i ≠ j),   T_ii = d_i := weil_matrix.py's diagonal.
(odd: a_i a_j (u_j − u_i) = (2/L') s w_i w_j (−2R_j/w_j + 2R_i/w_i) = −2(2/L') s (w_i R_j − w_j R_i).) The even j = 0 row is the case w_0 = 0, u_0 = 0. ∎
So T_N is "diagonal + rank one + Cauchy-like" (the particle form of ../ExploreFMM); deig.c builds it from 5n generator mid-points.

## 2. Lemma G (rigorous Galerkin upper bounds; point.py)
For exact-dyadic y_1, y_2 ∈ R^{N+1} (even) and z_1 ∈ R^N (odd), let H = (y_a^T T_N y_b), G = (y_a^T y_b) (2×2), all enclosed in Arb from the generator
balls of sec. 1 (no rounding of T_N: every entry is evaluated from the balls). Then
  λ_1^even ≤ U := H_11/G_11,   λ_2^even ≤ U2 := largest root of det(H − tG) = 0,   λ_1^odd ≤ Uo := z_1^T T_N z_1/z_1^T z_1.
Proof: f = Σ y_j φ_j lies in the form domain (E_N ⊂ domain, CCM sec. 3.3) with Q(f) = y^T T_N y and ‖f‖² = y^T y (orthonormal basis).
Courant–Fischer for A_λ (bounded below, compact resolvent: CCM sec. 3.3, Botha sec. 2.3): λ_1 ≤ Q(f)/‖f‖² for any f ≠ 0 of the parity, and
λ_2^even ≤ max over the 2-dimensional even space span(f_1, f_2) of the Rayleigh quotient = largest generalised eigenvalue of (H, G)
Hypotheses: y_1, z_1 ≠ 0 and y_1, y_2 linearly independent, i.e. G_11 > 0 and det G > 0 (reviewer's point, sec. 5). Verified by Arb
semantics at every point: a quotient by a ball containing 0, or a square root of a ball meeting (−∞, 0), is infinite or NaN, and every printed
U, U2, Uo is finite (y_1, y_2 are Gram–Schmidt orthonormal up to rounding, G ≈ Id). The vectors are deig.c's float eigenvector
approximations; rigour does not depend on their quality. The hybrid step (sec. 3) uses only U (one nonzero vector). ∎

## 3. Lemma M and the hybrid step test
(Monotonicity, Botha lem:mono = ../WeilSimple/DERIVATION.md sec. 5) every λ_k^σ is non-increasing in c. Hence a bound computed at
c_q ≤ e^{2L} is an upper bound at L: λ(e^{2L}) ≤ λ(c_q) ≤ U(c_q). Our grid takes c_q = ⌊e^{2L}·10^6⌋/10^6 (Python Decimal, 50 digits).
Hybrid step: if U(c_q(L_i)) < min(μ_odd, μ_l2)(L_{i+1}), with Botha's certified lower bounds μ (his Theorem V, imported), then for every
L ∈ [L_i, L_{i+1}]: λ_1^even(L) ≤ λ_1^even(L_i) ≤ U < μ(L_{i+1}) ≤ min(λ_1^odd, λ_2^even)(L_{i+1}) ≤ min(λ_1^odd, λ_2^even)(L)
(Botha cor:step with our U in place of his fr_sweep U(L_i)). This removes fr_sweep from the trust base; his μ remain imported.
Consistency test (falsification of his μ): μ_odd(L) ≤ Uo(c_q(L)) and μ_l2(L) ≤ U2(c_q(L)) must hold; a violation refutes his record.

## 4. The falsification scan (orientation only; no rigour claimed)
deig.c: T_N from the generator mid-points (sec. 1), dense LDLᵀ in place (packed rows; FLINT arb_approx_dot, mid-points only; no pivoting:
T_N is positive definite here, every pivot printed positive, "neg 0"), block inverse iteration V ← T_N^{-1}V with Gram–Schmidt (3 even, 2 odd
vectors, 4 steps; the per-step change of log10 diag(H) is ≤ 1e-7 after step 1 because e2/e1 ≥ 1.4e9 and o2/o1 ≥ 5e8 in every run), Rayleigh–Ritz matrix
H = (Lᵀ Y)ᵀ D (Lᵀ Y) (no cancellation: D > 0), eigenvalues by mpmath (post.py). Precision p ≈ log2(1/e1) + 160 bits: the LDL backward error
n·2^-p·‖T_N‖ is ≤ 2^-140·e1. Validation: c = 30, N = 240: e1, e2, e3, o1, o2 equal those of a dense python-flint arb_mat computation
(scan.py) to 12 digits, and e1 = CCMFMM's certified U(30). A crossing in the true spectrum would need o1/e1 or e2/e1 → 1; the truncated
ratios are converged to ≤ 5 % (REPORT sec. 2) and stay ≥ 5.5e4.

## 5. The chain beyond c = 30 (our own certificate; ../CCMFMM/DERIVATION.md, unchanged)
ccm/ holds ../CCMFMM's chain_ccm.py, run_ccm.py, gstate.py, prep_ccm.py, o1fix.py, weil_matrix.py, slt.c, rscan.c, copied unchanged except
prep_ccm.py's search path for ../WeilZeta/cert_zeta.py (one extra entry) and, from the second step on (19:33 UTC), chain_ccm.py's ladder beyond
c = 30: one rung ν = ε/10 instead of ν = ε/10, 10ε (time; any rung with the right inertia and ν > κη certifies the level ε, the second rung only
brackets λ_min(Id − W)). Binaries rebuilt with CCMFMM's commands (REPORT.md there, l. 170–171).
Rules (CCMFMM sec. 0, Lemmas P', I', R, B'): at c_{i+1}, both certificates (odd sector; deflated even sector with g from gstate.py at NG)
at ε = L(c_{i+1}) = 2.00001·U(c_i), U(c_{i+1}) = Q(g)/‖g‖² in Arb, scan premise rscan.c with K2 = max(401, ⌈5c²/3⌉); monotonicity covers
[c_i, c_{i+1}]. Start: CCMFMM's chain end c = 30, U(30) = 3.27353749418e-150. Step heuristics (logs/chain_int.jsonl): o1 at integers from our
float scan, scaled by 0.6 (the truncation overestimate seen at c = 30.57); they choose step sizes only and carry no rigour.

## 6. The s(L) test (scheck.py)
C_b = Σ_{log n < 2b} ω_n (τ_{log n} + τ_{−log n}) on L²(−b, b) is bounded self-adjoint and commutes with parity, so sup ⟨C_b h, h⟩/‖h‖² is the
larger of its even and odd suprema, each ≥ the top eigenvalue of its compression to E_N (Galerkin, from below). For f = Σ y_j b_j,
⟨C_b f, f⟩ = 2 Σ_n ω_n C_f(log n) = yᵀ(2 PR_N)y with PR_N = Σ_n ω_n (C_ij(log n)), C_ij(x) = ∫ b_i(y) b_j(y − x) dy symmetric in i, j for
these bases (weil_matrix.py's closed forms, 0 ≤ log n < 2b). So Ritz_N ≤ sup ≤ s is necessary; Ritz_N > s would refute (V1). Double precision
suffices: the entries are O(1) and the observed margins (≥ 0.0044) exceed the rounding error by > 10 orders. b = ⌈2^24 L⌉/2^24 (= L here).
