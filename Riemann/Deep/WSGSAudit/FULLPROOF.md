# FULLPROOF: the WeilSimple chain, step by step, with source lines (wsgsaudit, 2026-11-20)

Sources: CCM = arXiv 2511.22755 TeX (fetched 2026-11-20); WS = Riemann/Deep/WeilSimple/ (read-only); WZ = Riemann/Deep/WeilZeta/.
Notation: c = λ², L = ln c = 2 log λ, I = [-L/2, L/2], x = log u; f real of one parity (complex f splits, Q(u+iv) = Q(u)+Q(v)).

## Step 1. The form (CCM eq. (3.10), TeX l. 535-542)
QW_λ(f) = ∫|f^(t)|² 2θ'(t) dt/2π + 2Re(f^(i/2) conj f^(-i/2)) - Σ_{1<n≤λ²} Λ(n)<f|T(n)f>, f^(s) = ∫ f(u)u^{-is}d*u (l. 433-436),
θ(t) = -(t/2)log π + Im log Γ(1/4 + it/2) (l. 459-461), so 2θ'(t) = φ(t) := -log π + Re ψ(1/4 + it/2).
1a. Pole: f^(±i/2) = <f, e^{±x/2}> ⇒ 2Re(...) = 2<f,cosh(x/2)>² - 2<f,sinh(x/2)>²: +2<f,cosh>² on even f, -2<f,sinh>² on odd f
    (WS weil_matrix.py pole_vec + block `ps`, l. 87-95). Also = W_{0,2}^# of CCM l. 491 (checked: ∫_{-L}^{L} C(y)2cosh(y/2) = <f,cosh>²-<f,sinh>²).
1b. Archimedean: Gauss ψ(z) = -γ + ∫_0^∞ (e^{-u} - e^{-zu})/(1-e^{-u}) du, u = 2x, Plancherel (1/2π)∫F1 F2 cos(tx) dt = C12(x):
    A(f1,f2) = -(log π+γ)<f1,f2> + ∫_0^∞ [2e^{-x}<f1,f2> - 2e^{x/2}C12(x)]/(e^x - e^{-x}) dx  (WS DERIVATION §2).
    Audit identity: A = -W_R with CCM's own matrix formula (weinfty, l. 711): the difference is
    ||f||²[log 4 + log((c-1)/(c+1)) - 2∫_0^L dx/(e^x+1) - log(1-c^{-2})] = 0 (each log cancels). Fubini: |F|² ~ t^{-2}, kernel ~ log t: absolutely convergent.
1c. Primes: (f^**f)(n) = C_f(log n) is supported in [-L, L], continuous, C_f(±L) = 0, so n ≤ c ⇔ n < c (atom n = c contributes 0);
    weil_matrix.py l. 24-28 keeps n < c (fmpq comparison, exact). Weight Λ(n)n^{-1/2}, both signs ⇒ -2Σ w_n C12(log n). ✓
1d. Basis: CCM U_n = L^{-1/2}e^{2πinx/L} on [0,L] (l. 355-357), V_n = U_n(log(λu)) (l. 562-564) = (-1)^n L^{-1/2} e^{i w_n x};
    QW_λ^N = restriction to span{V_n : |n| ≤ N} (l. 812-813); γV_j = V_{-j} (l. 839) is f(x) → f(-x). Even part = span{cos w_j x, j ≤ N},
    odd part = span{sin w_j x, 1 ≤ j ≤ N}: exactly WS's φ_j, ψ_j (orthonormal). Eigenvalues are basis-independent ⇒ even-simplicity
    (CCM Def. even-simple, l. 852) of WS's two blocks ⇔ of CCM's matrix τ. ✓
1e. Numerics (logs/indep_c5.log, own code indep.py): 25 entries (even j ≤ 9, odd j ≤ 9) by x-space Gauss-Legendre quadrature with
    numerically computed autocorrelation vs weil_matrix Arb closed forms: max diff 2.5e-14 (double precision); lattice form of Step 2
    with my own μ(k) on cos²(πx/L) and sin(ax)+sin(2ax)/2 vs the matrix: 6.674893046841396e-5 vs ...583e-5 (truncated tail 2e-18), odd equal.
    A first 40-digit attempt (mpmath tanh-sinh, killed by memcap at 360 MB after 3 entries) matched (0,0), (0,1), (0,2) to 6e-40.
    CCM's own t-integral ∫φ|F|²/2π for the (0,0) entry at c = 5 (t ≤ 4000 + asymptotic tail): 2.9e-11 from the closed form (logs/tcheck_c5.log).
    Rough regime (f = φ_0, jumps at ∂I, F(r_k) ~ 1/k): lattice sum to K plus tail vs T[0,0]: 4.2e-9 (K = 2·10^4), 1.9e-10 (K = 10^5) (logs/latrough_c5.log).

## Step 2. Lattice identity on all of L²(I) (WZ DERIVATION §3, Chekmachev Thm 3.1; code WS cert_simple.py l. 47-79)
For f ∈ L²(I), r_k = πk/L, as an identity in (-∞, +∞]:  Q(f) = pole(f) + (1/2L) Σ_{k∈Z} μ(k)|F(r_k)|²,
μ(k) = -log π + Re ψ(1/4 + ir_k/2) + 2(-1)^k Σ_n c^{-b_n} b_n/(b_n²+r_k²) - 2Σ_{m<c} w_m cos(r_k ln m), b_n = 2n+1/2.
Audit: the (-1)^k term is the restriction of the kernel -Σ_n e^{-b_n|x|} to [-L, L] at r_k (∫_L^∞ e^{-bx}cos(r_k x)dx = (-1)^k c^{-b} b/(b²+r_k²)), re-derived;
the singular part is Σ_k |F(r_k)|²(1 - cos r_k x)/x ≥ 0 termwise (Tonelli). SymZeta.mu = own mpmath μ to 2e-37 (8 nodes, k ≤ 12345);
the identity reproduces the closed-form matrix (Step 1e). Q is lsc and bounded below on L², = CCM's closed form (E is a core, CCM l. 569-576);
since Q_closure(f) = liminf_{E∋g_n→f} Q(g_n) ≥ Q(f) by lsc, every lower bound proved for (3.10) on L² holds for A_λ.
Sectors (f even: φ_j = N_j cos w_j x, e_k = cos(πkx/L), k odd; f odd: sin, j ≥ 1): Σ_{k∈Z}|F(r_k)|² = 2L||f||² (Parseval on [-L, L]) and
Σ_j <f,φ_j>² = ||f||² give (1/L)Σ_{k odd≥1}<f,e_k>² = ||f||²/2, hence for every t1 and ε the identity
  Q(f) - ε||f||² = Σ_j [(μ(2j)+t1)/2 - ε]<f,φ_j>² + (1/L)Σ_{k odd≥1} (μ(k) - t1)<f,e_k>² ± 2<f,v>²   (v = cosh(x/2) even, sinh odd).
Closed forms <φ_j,e_k>, <φ_j,v>, <e_k,v>, ||v||² (cert_simple.py l. 301-311) = quadrature to 2e-31 (logs/indep_c5.log).

## Step 3. Minorant, Woodbury, positivity, deflation (cert_simple.py l. 259-455)
3a. Drop k ∉ K (needs μ(k) ≥ t1 = 1 for every odd k ∉ K) and use dl_j ≥ deff = t1/2 - ε for j ∉ E (needs μ(2j) ≥ 0): l. 259-290.
    Weights are Arb lower bounds (l. 283-287: w_k = ((μ-1)/L).lower(), dl_j = (...).lower()): lowering a coefficient only weakens the minorant.
3b. Infinitely many nodes: for k > ks, μ(k) ≥ monotone(ks+1) = h(r) - σ2 - 2(tailsum + trem) ≥ t1 (l. 76-89): h(r) = -log π + Re ψ(1/4+ir/2)
    is increasing (∂_y Re ψ(x+iy) = Σ 2y(n+x)/((n+x)²+y²)² > 0), b/(b²+r²) decreasing, |2Σw cos| ≤ σ2 = 2Σ|w|. Same for even j with threshold 0.
    Nodes in (K2, ks] by rscan.c: exact 64-bit phases, 2^16 table + Taylor, block bounds H ≤ h, T ≥ 2tailsum + 2trem (monotone), a-priori
    error E (l. 152-156) with E < margin/10 asserted; unflagged ⇒ μ(k) ≥ μ_c - E ≥ thr + 1e-6 - E; flagged nodes are simply kept in K/E (l. 259).
    Empirical: rscan_dbg.c (rscan.c + print) vs Arb at c = 19, 4204 nodes up to k = 2·10^9: max error 1.53e-9 = 0.385 E (logs/rscan_test.log).
3c. P0 = deff Id + Σ_{j∈G}(dl_j - deff)φ_jφ_j^T (G = {dl_j > 0}), P = P0 + Σ_{K+} w_k e_k e_k^T (+2vv^T even) (+M gg^T defl), N = BB^T with columns
    √(deff-dl_j)φ_j (dl_j ≤ 0), √(-w_k)e_k (w_k ≤ 0) (+√2 v odd). Then Q - ε||·||² (+M<·,g>²) ≥ P - N on L² (bounded operators).
    P0^{-1} = (1/deff)Id + Σ_G (1/dl_j - 1/deff)φφ^T; Cauchy sums by partial fractions (l. 352-354; distinct y: k² > 0 > y_v).
    g-entries (l. 345-350): <g, P0^{-1}b> = Σ_l (g_l/p_l)<φ_l, b>, p_l = dl_l (l ∈ G) else deff, l ≤ NG ≤ K2 asserted in E (l. 297). ✓
3d. P - BB^T ≥ 0 ⟺ B^T P^{-1} B ≤ Id; Woodbury W = X - Y^T S^{-1} Y, S = D^{-1} + U^T P0^{-1} U > 0 (D > 0); for any exact Z,
    Y^T S^{-1}Y ≥ Z^TY + Y^TZ - Z^TSZ ⇒ W ≤ W_up; code l. 375-397 forms Id - W_up = Id - X + (Z^TY)^T - Z^T(SZ - Y) blockwise (checked index by index;
    each column block of X is overwritten only after its own use).
3e. Positivity (l. 403-455): L = approximate Cholesky (exact dyadic, lower triangle only, arb_mat zero-initialized), E = IW - LL^T - s0 I enclosed
    in Arb; max row sum ≤ s0 - ν ⇒ λ_min(M) ≥ s0 - ‖E_M‖ ≥ ν > 0 for every symmetric M in the ball matrix (Gershgorin/∞-norm). ✓
3f. Deflation: Q(f) - ε||f||² + M<f,g>² ≥ 0 ∀ even f ⇒ Q(f) ≥ ε||f||² for even f ⊥ g ⇒ λ_2^e ≥ ε (every 2-dim subspace meets g^⊥).
    Falsification (logs/defl_falsify.log, c = 5, ε = 2.443577e-17 > U(5)): normal g passes; g ⊥ ξ_N (<g,ξ> = 1e-57) FAILS (margin -7.2e-18);
    M = 1e-40 FAILS (-7.6e-18), as they must (f = ξ_N gives Q(f) - ε||f||² < 0).

## Step 4. Logic (WS DERIVATION §3, §5; run_windows.py l. 58-74; chain.py l. 27-50)
4a. Parity: Q(γf) = Q(f) term by term (<γf,cosh> = <f,cosh>, <γf,sinh> = -<f,sinh>; φ even; ±log n both summed), so cross terms vanish,
    Q(f_e + f_o) = Q(f_e) + Q(f_o), A_λ and every E_N commute with γ.
4b. Dimension count (no discreteness needed): m = min(L_odd, L_2) > U. Ran 1_{(-∞,m)}(A_λ) is γ-invariant; an odd vector in it has
    Q(f) < m||f||² ≤ L_odd||f||², contradicting (b); its even part has dim ≤ 1 (a 2-dim subspace meets g^⊥, contradicting 3f) and ≥ 1
    (Q(g) ≤ U||g||² < m||g||²). So the bottom of σ(A_λ) is a simple eigenvalue μ_λ ≤ U with even eigenvector. ✓ (reviewer: CORRECT)
4c. Truncations: E_N = span{V_n : |n| ≤ N} (CCM l. 812) ⊂ form domain; compressions have per-sector min-max values ≥ the full ones, so
    λ_1^o(N) ≥ L_odd, λ_2^e(N) ≥ L_2, and λ_1^e(N) ≤ U once g ∈ E_N (N ≥ NG). Hence QW_λ^N is even-simple (CCM Def. even-simple, l. 852) for every
    N ≥ NG; Thm 5.10 (finmain, l. 1087; called "Theorem 5.10" at l. 1392) needs exactly this: its normalization δ_N(ξ) = 1 is derived from
    even-simplicity in CCM (l. 856-864). Our φ_j/ψ_j are an orthonormal basis of the same E_N, so nothing else is needed. ✓
4d. Monotonicity: L²(I_c) ⊂ L²(I_{c'}) (c < c', extension by 0) preserves parity and Q (atoms c ≤ n < c' meet C_f only where C_f = 0), so every
    per-sector min-max value is non-increasing in c. Chain step: for c ∈ [c_i, c_{i+1}], λ_1^e(c) ≤ Q(g_i)/||g_i||² = U(c_i) < ε_{i+1}
    ≤ λ_1^o(c_{i+1}) ≤ λ_1^o(c) and ≤ λ_2^e(c_{i+1}) ≤ λ_2^e(c); then 4b at c (with min-max values) gives the simple even ground state. ✓
4e. Logs (chaincheck.py, exact Decimal/Fraction): 37 contiguous steps 2 → 20, Uprev of step i+1 = U_up of step i (string-equal), Uprev < ε, both
    certificates rigorous; 19 integer windows U < min(L_odd, L_2). 52 windows (19 integers + 38 chain points, 5 shared) have U(c) < a certified
    lower bound at some c' ≥ c with g ∈ E_NG: N₀ = NG ∈ [30, 160]. table.py (copy) agrees: "covered interval [2, 20], 37 steps".
    Ratios min(L)/U at the integer windows: 7.68 (c = 2) to 151 (c = 20). Note: chain points use ε = L_odd = L_2 = 2.00001 U(c_i).

## Step 5. OEIS A401172 (bracket.py, abracket.py, a401172.py, uimprove.py; logs/a401172_table.md)
5a. Definition (OEIS entry %N, fetched 2026-11-20): a(n) = ⌊-log10 κ(n)⌋, κ(n) = inf Q(f)/||f||² over nonzero smooth f supported in
    (-log(n)/2, log(n)/2), Q the zeta Weil form (%C formula = Step 1). Data n = 2..9: 2,7,12,17,22,27,32,37 (keyword "more").
5b. κ(n) = bottom of A_λ at c = n: lower bounds hold on all of L²; for the upper bound U (g ∈ E, nonzero at ∂I) one needs C_c^∞(I°) to be a form
    core: dilation f(x/s) → f (∫φ(u/s)|F|² → ∫φ|F|² by dominated convergence, φ(u/s) - φ(u) bounded) and mollification. Not stated by WS; holds.
5c. Lower bound: plain even certificate at ε (sector even, no deflation) ⇒ λ_1^e ≥ ε; odd: λ_1^o ≥ L_odd(n) > U(n) ≥ ε (int.jsonl), so κ(n) ≥ ε.
    Upper: rigorous Rayleigh quotients (int.jsonl, uimprove.jsonl; Arb, exact dyadic g). Decision: ⌊-log10 U⌋ = ⌊-log10 L⌋ = a ⇒ 10^{-a-1} < L ≤ κ ≤ U ≤ 10^{-a}.
5d. Rechecked from the logs: a(10..13, 15, 16, 18, 19) = 42, 48, 53, 58, 69, 74, 85, 90; n = 2..9 equal the OEIS data; n = 14, 17, 20 open.
    n = 10, 16, 19 rest on L = 1.00001·10^{-a-1} (abracket.jsonl, rigorous); n = 16 rerun here: margin 2.95443483618e-76, identical.
    Cross-check: OEIS %C cites Q(g)/||g||² ≤ 8.62132e-49 at n = 11 (Wiktorowicz): inside our [2.7411e-49, 9.45177e-49], no conflict.

## Step 6. Reruns (copy /tmp/wsgs_copy, one process at a time, RSS by /usr/bin/time -l, all < 290 MB) and falsification
| run | result vs log | RSS MB |
|---|---|---|
| c = 2, rerun.sh | U_up 0.00133518326736, L_2 2.970683e-2, L_odd 1.025550e-2 identical; margins 0.3404/0.0837 (log 0.3402/0.0821: log used K2 = 201, rerun 401) | 67 |
| step 29/4 → 31/4 | both rigorous at ε = 7.350532e-29; U(31/4) = 1.08221871520e-31; margins identical | 74 |
| step 16 → 133/8 | both rigorous at ε = 7.486021e-75; U(133/8) = 1.99822754250e-78; margins identical | 157 |
| c = 9, 12 (CCM's λ = 3, √12), run_windows from scratch incl. NG search | identical U, L_2, L_odd, margins; NG = 80 ≤ 120 | 101 |
| c = 19 (L = 1.472), defl and odd at ε = 5.809267e-89 | both rigorous; U(19) = 4.22270760017e-91 identical; odd margin identical | 268 |
| OEIS n = 16 decision (even, ε = 1.00001e-75) | rigorous, margin 2.95443483618e-76 identical | 154 |
| table.py on logs/chain.jsonl | 37 records ok, covered [2, 20] | - |
Not rerun (above the cap): integer windows 18, 20 and chain steps 37/2 → 20 (312-676 MB); the orchestrator reran 16 → 20 on 6 Nov 2026.
Falsification (all behave as a sound certificate must; logs/negctl.log, defl_falsify.log):
| test | expected | got |
|---|---|---|
| odd, ε = 1.02 λ_1^o(E_60) at c = 5; ε = 1.02 λ_1^o(E_80) at c = 31/4 (truncated eigenvalues are upper bounds) | fail | fail (margins -3.5e-15, -2.5e-28) |
| defl, ε = 1.02 λ_2^e(E_N), same windows | fail | fail (-8.0e-13, -8.8e-26) |
| odd, ε = 0.5 λ_1^o(E_N) (tightness) | - | pass: the minorant is within a factor 2 of the truth |
| stretched step 2 → 3 at ε = 2.00001 U(2) | fail | fail (-3.3e-3) |
| odd at ε = 2.00001 U(2), c = 9/4, 7/3, 19/8, 29/12 / 5/2, 21/8, 11/4, 3 (λ_1^o(E_60) = 1.3e-2 ... 3.6e-3 / 1.8e-3 ... 1.6e-5) | pass ⇒ ε ≤ λ^N | pass / fail, consistent at all 8 |
| c = 7 odd, K2 = 201 / 401 | fail / pass | fail (-2.66e-24) / pass (1.0e-24) |
| defl at c = 5, ε = 2.443577e-17 > U(5): g normal / g ⊥ ξ_N / M = 1e-40 | pass / fail / fail | pass / fail / fail |
| rscan.c rounding, c = 19, k ≤ 2·10^9 | ≤ E | 0.385 E |
