AUDITED SOUND UP TO: one unrepeated large run. The authors ran the 7300-bit core (K = 1088), the 192-bit assembly through M = 3,678,208 modes and the 918,464-pivot LDL once, and deposited only SHA-256 hashes and receipts of their core files and input matrices; a re-run needs about 6 GB of memory and 1-3 core-days per sector. All three gaps the authors declared are closed: the binary128 run plus our double-double division lemma, and our exact enumeration of the prime-norm breakpoints, which proves rho_W <= 11.101037276776648 < 11.101037276778 (the value used).

# WiktAudit: is Wiktorowicz-Wiktorowicz's Weil positivity certificate at half-width 2.3029 rigorous?

Audited: R. Wiktorowicz, B. SI Wiktorowicz, "An Unconditional Exponential Upper Bound for the Windowed Weil Quadratic Form, with Computer-Assisted Positivity Certificates", Zenodo 23869032 v3, 1 Nov 2026 (concept 23447302). Its supplement ZENODO_v3_n100.zip has code (kod/), receipts (wyniki/) and rho_w_cert_n100.json. The method is that of rStation paper #60 (bramble_trace, "Full-source Schur certificates for Weil positivity at support radius two", 26 Oct 2026, a = 2). That paper uses the Green identity of rStation #45, Lemma 3.1 (H. Li, a = (log 11)/2).

RH relevance: N1 (Weil positivity on windows; RH iff positivity on every window).

## 1. Verdict in short
- **Claim.** Let a = 23029/10000, so e^{2a} = 100.063. The certificate claims inf sigma(A_a) > 0 in both parity sectors, hence Q(f) > 0 for every nonzero complex f supported in [-a, a]. Q is exactly our functional (ChekAudit sec. 1) with c = e^{2a}. The certificate gives no explicit eps (sec. 2, step 10). It covers every smaller window, so kappa(n) > 0 for 2 <= n <= 100 (OEIS A401172).
- **Mathematics.** No false or unjustified step. I re-derived every reduction (sec. 2), including the a-dependent pieces that the authors carried over from the a = 2 text.
- **Their gap (i), no error bound for double-double division: closed.**
  - The binary128 run certifies both sectors by itself: IEEE correctly rounded arithmetic, same input SHA-256, margins 0.0173269704163 (even) and 0.0170198143082 (odd).
  - The double-double run is also covered: + and * are Joldes-Muller-Popescu 2017, Algorithms 6 and 10 (Theorems 3.1 and 5.1). Their division (QD's accurate_div) has relative error < 10.0001 u^2 < 2^-102.6, which we prove in DDDIV.md. The budget assumes 2^-95.
- **Their gaps (ii) and (iii), breakpoint completeness by count only and the rho_W generator not deposited: closed by recomputation.**
  - rhow_exact.c builds the operator from scratch. It enumerates all breakpoints exhaustively with exact integer labels; no floating-point comparison decides anything.
  - Level 6 has 14,951,282 cells, exactly their count. Schur test with h = B^5 1: the proved bound is ‖B_p‖ <= ‖B‖ <= 11.101037276776648, below the 11.101037276778 used. 124 s, 75 MB.
- **Further independent checks.** All match to the stated precision.
  - Our x-space Weil form (no Green identity) against their matrix entries.
  - Their pipeline code at small size against dense Schur complements.
  - All a-dependent constants, the omitted tail, the assembly error and the LDL budget, recomputed from the receipts.
- **Remaining.** The big Arb and LDL computations were executed once and cannot be repeated from the deposit (sec. 5). They are rigorous by construction (Arb balls; a proved backward-error bound) if the code ran as deposited. We read that code line by line and ran it at small size.
- **Consequence for N1.** If accepted, the window record becomes half-width 2.3029 (c = 100.06), positivity only. Among certificates with an explicit eps, Chekmachev's c <= 47 stays the record.
- **Novelty pass.** rStation #60 already claims a = 2 (c = e^4 = 54.6 > 47) on 26 Oct 2026, before Chekmachev. TREE.md has no rStation entry (sec. 6).

## 2. The inequality chain
H = (a/2) E A E with E_qq = (-1)^q, in the parity bases on (0, a): even sqrt(2/a) cos(q pi x/a) (q = 0: 1/sqrt(a)); odd sqrt(2/a) sin((q+1/2) pi x/a). Put k = q (even) or q + 1/2 (odd), lambda_q = k^2. A positive scalar and a diagonal unitary keep the sign of the form.

Labels used below:
- (A): rigorous as deposited (Arb or exact).
- (B): floating point, closed by a rigorous error analysis.
- (C): needed recomputation; done here.
- (D): unverifiable.

| # | step (statement) | where | label | our check |
|---|---|---|---|---|
| 1 | Q on C_c^inf(-a,a) is semibounded and closable; it splits orthogonally by parity, each sector unitarily equals L^2(0,a); smooth tests are a form core; zero extension gives every smaller window | #60 sec. 2, #45 sec. 2 | A (math) | re-derived. Form = ChekAudit's Q_c (pole 2 Re(fhat(i/2) conj fhat(-i/2)), primes n < e^{2a}) |
| 2 | Green identity: A = D + B_p + sum_j r_j (x) r_j + eps p (x) p. D = diag(Re psi(1/4 + i pi k/(2a)) - log pi); r_j >= 0 rank one, c_j = 2j + 1/2; B_p prime operator (shifts by log m and the reflection log m - x); p_eps = e^{x/2} + eps e^{-x/2} (negative rank one when odd) | #45 Lemma 3.1, #60 (2.5)-(2.6) | A (math) | re-derived (digamma expansion (3.3), 2c R_c - E_c = r (x) r). Bounded parts: responses (constant plus Carleman kernel), primes, pole |
| 3 | Displacement rank 2: H_qr = (u_q v_r - v_q u_r)/(lambda_q - lambda_r), plus the exact diagonal (3.5)-(3.6) and Arb tail balls (3.7) with e^{-2a c_j} | source.py l. 31-87 | A | **x-space check: our own Galerkin matrix (Weil regularisation, closed-form autocorrelations, no Green identity) equals their H to <= 8.4e-72 relative, both parities, indices up to 100, at a = 2.3029** (xcheck_both.log) |
| 4 | Core q < K = 1088 at 7300 bits: exact dyadic T, ‖T M0 T^T - I‖_inf <= eta < 1, so M0 > 0 and M0 >= (1-eta)/‖T‖_F^2 = 4.58e-530 (even), 2.72e-524 (odd). Certified solves a = M0^{-1} u, b = M0^{-1} v (residual/gap) | core.py l. 60-102; #60 sec. 4 | A as deposited, not re-run | logic re-derived; the same code at K = 24 matches dense Arb solves (test_pipeline_*.log) |
| 5 | Post-core S = H_{>=K} - C^T M0^{-1} C has generators (U,V) = T(lambda_q)(u,v) and diagonal S_qq = H_qq - g^T T^T J T' g. The sums use a 1400-bit convolution (q < 69,632), then 96 moments with remainder (6.3) as balls | pipeline.py l. 9-87, tail.py l. 62-74; #60 (5.1)-(5.6), (6.2)-(6.3) | A as deposited | identities re-derived (displacement Schur complement; (5.6) and moment remainders checked). **Small run: U, V, S_qq balls contain the dense values (rel. 1e-53 to 1e-55 for q >= 100) in both the convolution and the moment branch, both parities** |
| 6 | Prime norm rho_W >= ‖B_p‖: Schur test sup B h/h with h = B^5 1 for the positive dominating operator B | tail.py l. 13-57, rho_w_cert_n100.json | was D (summary JSON only), now C done, A | **rhow_exact.c: 14,951,282 cells, every comparison exact, max ratio 11.101037276776603, proved <= 11.101037276776648 < 11.101037276778** (sec. 4) |
| 7 | Distant block q >= N = 919,552: S_{>=N} >= diag(t_q), t_q = (a/2) d_k - (a/2) rho_W - reserve. The reserve 0.227 covers the trace of the core covariance (0.13896 even, 0.13803 odd), the odd pole tail and the dropped PSD responses and even pole. t_N = 1.0434 > 0, t_q increasing | tail.py l. 97-154; #60 (6.4)-(6.10) | A | t_N, t_M, J2, J4, abs(r_q) < 0.591, pole tail recomputed in mpmath, equal to 50+ digits (check_receipts.log). The covariance needs the core: A as deposited |
| 8 | Finite Schur matrix on [K, N): calS = PSP - C T^{-1} C^* = calS_M - R_M. calS_M keeps all columns N <= q < M exactly (192-bit Arb, displacement rank 4, (7.3)-(7.4)) and is stored as binary64 centres with operator error (7.5) <= 7.19e-10 / 7.46e-10, using H_* = (1 + log(N-K))/K | pipeline.py l. 89-147; #60 (7.1)-(7.5) | A as deposited | (7.2)-(7.5) re-derived. **Small run: binary entries equal the dense Schur complement to binary64.** (7.5) recomputed from the receipt |
| 9 | Omitted tail ‖R_M‖ <= (sqrt(SU) vM sqrt(J4) + sqrt(SV)[aM sqrt(J2) + (aM+bM) sqrt(J4)])^2 / (t_M [1 - (k_{N-1}/k_M)^2]^2) = 0.28267 (even) / 0.28298 (odd) | tail.py l. 156-161; #60 (7.6)-(7.8) | A | Hilbert-Schmidt and Minkowski steps re-derived; value recomputed from the receipt fields (SU, SV, aM, bM, vM come from the run) |
| 10 | Point LDL of S_input - 0.30 I. Structured backward error <= n[64uG(1+C/p0)H_* + 4u(D + C^2/p0)] + 1e-12 (local residuals zero-extended and telescoped; row-sum norm; maxima monitored in the run) = 1.0e-12 (binary128) and 6.4e-12 (double-double). Hence calS >= (0.30 - residual - assembly - omitted) I = 0.017327 / 0.017020 > 0. Two square completions ((9.1) with T >= t_N > 0 and C bounded; (9.2) with M0 > 0) give H >= eps I for some eps > 0 in both sectors | structured_ldl.cpp l. 26-39 (double-double), 76-123 (LDL), finish.py l. 8-61; #60 sec. 8-9 | B, closed | (8.3)-(8.5) re-derived, including factor rounding (first-column error <= 9uG/abs(Delta lambda), cross term 18uGC/p0, generator term 16uG(1+C/p0)). Budget recomputed (check_receipts.log). Binary128 is IEEE; double-double per DDDIV.md. Explicit eps would need bounds on ‖M0^{-1} C_core‖ and ‖C‖, which are not deposited |

First false or unjustified step: none.

## 3. The three self-declared gaps
- **(i) Double-double division.** Closed twice.
  - The binary128 run alone suffices. GCC's __float128 + - * / (libgcc soft-fp) are correctly rounded. The code uses only those, plus exact conversions from binary64. The shift 0.30 is parsed by strtoflt128; its error, below 1e-33, lies inside the 1e-12 slack. Unit roundoff 2^-113, and finish.py uses 2^-112.
  - For double-double, operator+ and operator* are exactly JMP 2017 Algorithms 6 and 10 (checked line by line against hal-01351529v3), with bounds 3u^2 + 13u^3 and 7u^2.
  - Our lemma bounds the three-digit QD division by 10.0001 u^2. Every operation is therefore within 2^-102.6, far inside the 2^-95 budget (DDDIV.md).
  - Their 27-28-digit agreement between the two runs is consistent with this, but our proof does not use it.
- **(ii) Completeness of the rho_W breakpoints, checked by count only.** Closed.
  - Our breakpoint set is closed under all 105 preimage maps by construction; by induction, h_k = B^k 1 is constant on every cell.
  - Our count, 14,951,282, equals the count of positions they checked. Their certificate lists 14,951,361 cells, which they explain as 79 window edges.
  - Minimum cell width 1.07e-11, as theirs.
- **(iii) rho_W generator not deposited.** Closed: rhow_exact.c, gen_consts.py and rhow_bound.py are a complete independent generator and checker. Our level-4 value 11.159163993665569 equals their quoted "order-4 value 11.159163994".

## 4. What we recomputed (all single-threaded, under 130 MB)
| check | result | file |
|---|---|---|
| Schur test for ‖B‖, exact labels (e,P,Q) = e a + log(P/Q). Comparisons: equal e exact in 128-bit; otherwise against integer enclosures of 2^56 e^{a} and 2^56 e^{2a} from 256-bit Arb balls; 0 undecided | levels 1-6: 43, 793, 11964, 148367, 1607080, 14951282 cells; max ratio 11.101037276776603; min h_5 = 112230.997 | rhow_exact.c, consts.h (gen_consts.py), rhow_exact6.log |
| rounding allowance: weights = smallest double above the Arb upper end; <= 35 positive products summed in round-to-nearest; Higham's bound for sums of positive terms and one division | factor 1 + 4.0e-15, bound 11.101037276776648, slack 1.35e-12 | rhow_bound.py, rhow_bound.log |
| first float version (double positions, TOL 1e-13, 0 near-ties) | identical numbers | rhow.c, weights.h, rhow6.log |
| own x-space Galerkin matrix against their source rows at a = 2.3029 | max relative difference 8.4e-72 (both parities; indices 0-5 and {0,1,7,30,63,100}) | xcheck_both.py/.log, xcheck_even.py/.log (uses ChekAudit galerkin_x.py) |
| their core.py and pipeline.py at K = 24, N = 512, M = 4096 against dense Arb Schur complements | U, V, S_qq enclosed in both branches; assembled binary64 entries equal the dense calS_M | test_pipeline.py, test_pipeline_even.log, test_pipeline_odd.log |
| a-dependent constants and the final budget, mpmath only: V_h = 2.8423, E_h = 645.81, t_N, t_M, J2M, J4M, r bound, odd pole tail, omitted (7.7), assembly error (7.5), residual (8.5), margins (dd and f128, both sectors) | all equal the receipts to 50+ digits | check_receipts.py/.log |

Notes:
- The negative controls in the deposit (rho_W x 0.99; list without 89) stop at the binding gate in tail.py before any computation. They test the gate, not the margin.
- The margin vanishes if rho_W is understated by about 1.2% (through t_M in the omitted tail). Our proved bound is 1.35e-12 below the value used.

## 5. What remains: exactly which computation
- **Core** (step 4): 1088 x 1088 Arb matrices at 7300 bits, about 1.15 GB each, several alive at once (M, L^{-1}, T, T M, T M T^T). Our estimate is >= 6 GB; their run took 3433 s with 16 Arb threads. Above our 1.5 GB limit. Its output a23029_*_core.json is not deposited (SHA-256 in sha_duze_pliki.txt).
- **Source, transfer and assembly** (steps 5, 7-9) through M = 3,678,208 need that core JSON. Their run took about 1280 s per sector. The receipt numbers that come from it: covariance 0.13896 / 0.13803, SU, SV, aM, bM, vM and the generator maxima.
- **LDL** (step 10) on the two 36.7 MB binaries a23029_*_pipeline.json.bin. These are not deposited (SHA-256 only). Cost: 918,464 pivots, about 4.2e11 pair updates. Their run took 20 core-hours (dd) and 54 core-hours (binary128) per sector.
- None of these needs the 369 GB machine specifically. A workstation with about 8 GB and one to three core-days per sector could repeat everything with the deposited code. Until someone does, trust rests on one execution of published, readable code. That code is internally consistent: dd and f128 runs agree on the same SHA-256, and every receipt identity we could test holds.

## 6. Novelty pass (6 Nov 2026)
- **arXiv.** python3 .scratch/lst.py 'Weil|positivity|quadratic form' math.NT for Nov and Oct 2026: 19 title hits, no windowed positivity certificate. 2610.32704 (RH-conditional spectral convergence) and 2611.02944 (PF orders of the Xi densities) are unrelated to window records.
- **Zenodo.** "Weil quadratic form", 15 newest: nothing beyond a = 2.3029. Related: 23983773 (conditional lower bound e^{-(4 pi + eps) n}), 23939977 (Markov representation), 23847003 (Chekmachev data).
- **Concept 23447302.** Three versions only: v1 23447303 (18 Oct, upper bound), v2 23673520 (26 Oct, kappa(23) >= 2^-68468), v3 23869032 (1 Nov, this certificate). No v4.
- **HAL.** "Weil positivity" and "Weil quadratic form": nothing relevant.
- **rStation** (all 113 riemann_hypothesis papers, via its API):
  - #45 (Li, a = (log 11)/2 = 1.199, 23 Oct);
  - #60 (bramble_trace, a = 2, 26 Oct; long-double LDL with u = 2^-63, 4-step prime graph rho = 8.920292);
  - nothing claims more than a = 2 except this deposit.
- **Our files.** No mention of rStation in TREE.md, DISCOVERIES.md or the NOVELTY files.
- **For the orchestrator:**
  - the N1 entry should list rStation #60 (half-width 2, c = 54.6, all f; same method; not audited by us) as a claim beyond Chekmachev's 1.925;
  - this deposit, at half-width 2.3029, is now audited sound up to the single unrepeated run.

## 7. Files (Riemann/Deep/WiktAudit/)
- REPORT.md (this file) and DDDIV.md (double-double error lemma).
- rhow_exact.c, gen_consts.py, consts.h, rhow_bound.py, rhow_exact6.log, rhow_bound.log: the exact rho_W certificate.
- rhow.c, weights.py, weights.h, rhow6.log: the first float version.
- xcheck_both.py/.log and xcheck_even.py/.log: the formula check.
- test_pipeline.py and test_pipeline_{even,odd}.log: the algebra check.
- check_receipts.py/.log: the constants and budget.
- Run:
  - cc -O2 -ffp-contract=off -o rhow_exact rhow_exact.c -lm && ./rhow_exact 6 > rhow_exact6.log && python3 rhow_bound.py rhow_exact6.log
  - uv run --no-project --with python-flint python3 gen_consts.py
  - uv run --no-project --with mpmath python3 check_receipts.py <wyniki>/tor1/results_f128
  - the two Python checks need the deposit's kod/ directory as their first argument.
- Downloads (deposit, rStation texts, JMP PDF) were deleted after use. State: .scratch/wiktaudit/STATE.md.
