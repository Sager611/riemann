AUDIT INCOMPLETE — no counterexample and no false step in Botha's Pieces III–IV; independently certified (our upper and lower bounds) on λ² ∈ [2, 385/12 (= 32.0833)] plus the windows 35, 40, and with our rigorous upper bounds replacing his fr_sweep on all of λ² ∈ [30, e⁴], his Theorem V lower bounds μ remain imported trust on (385/12, e⁴] except at 35, 40.

# BothaAudit — Botha's simple even ground state on L ∈ (1.7006, 2] (λ² ∈ (30, e⁴]): falsification, audit, independent extension (bothaaudit, 14 Nov 2026)

## 1. Short answer
Botha's Theorem 1 survives every test we could run on L ∈ (1.7006, 2], and part of it is now independent of his trust base:
| part of (30, e⁴] | lower bounds (λ_2^even, λ_1^odd) | upper bound (λ_1^even) | status |
|---|---|---|---|
| [30, 385/12 (= 32.0833)] | ours (CCMFMM chain continued, sec. 5) | ours | independently certified |
| 35, 40 (single windows) | ours (CCMFMM) | ours | independently certified (CCMFMM, before this round) |
| (385/12 (= 32.0833), e⁴] | Botha's Theorem V (certify.py), imported | ours (point.py, sec. 4) | certified conditional on his μ only |
| whole range | — | — | float scan: no crossing, gap ratio ≥ 5.5e4 (sec. 2); all 62 of his μ below our upper bounds (sec. 4a) |
Before this round the stretch (30, e⁴] \ {35, 40} rested on two imported pieces (his μ and his fr_sweep U); fr_sweep is now redundant.

## 1b. Object
Botha, Zenodo 23852617 v2.0.0 (1 Nov 2026; earlier version v1.1.1 = 23832221; no later version): Theorem 1, the Weil form QW_λ of
Connes–Consani–Moscovici (arXiv 2511.22755, = his Q_L with L = log λ, R415) has a simple, even ground state for every L ∈ (0, 2].
Pieces: I–II (0, 1.19]; III [1.19, 1.4035] (Elston's certifier + deflation); IV [1.4035, 2] (rh_60-route certifier code/rh60route/certify.py,
53 windows, odd and deflated-even lower bounds μ, upper bounds U(L) from code/upper/fr_sweep.py, step test cor:step).
Our coverage before this round: λ² ∈ [2, 30] (WeilSimple [2, 20], CCMFMM [20, 30]) and the windows λ² = 35, 40.

## 2. Falsification scan (float level): no crossing, no near-degeneracy on λ² ∈ [30, 54.6]
Grid: c = 30 and Botha's 31 windows c_q = ⌊e^{2L}·10^6⌋/10^6 = 30.570758, …, 54.598150 (step ≈ 0.8 in c; logs/grid.txt). Per point: generators in Arb
(gen.py), the truncated even block (n = N + 1) and odd block (n = N) built from them and factored in C (deig.c: dense LDLᵀ in place, FLINT
arb_approx_dot at 670–1107 bits, ≈ 3.3·(150 + 5.46(c − 30)) + 160), block inverse iteration with 3 (even) and 2 (odd) vectors, Rayleigh–Ritz; then
point.py's rigorous Rayleigh quotients (sec. 4). N = max(320, c²/5) (320 … 600). Validation: at c = 30, N = 240 the C path gives U = 3.27353749418e-150,
CCMFMM's certified U(30) to all 12 digits, and the Python Arb path (scan.py) agrees to 12 digits for e1, e2, e3, o1, o2. No negative pivot anywhere.
Results (logs/points.jsonl, falsify.out):
| quantity | min over the grid | max | trend |
|---|---|---|---|
| o1/e1 | 5.52e4 (c = 30.57) | 1.78e5 (c = 54.60) | rising overall (wiggles from the N rule) |
| e2/e1 | 1.58e9 (c = 30.57) | 1.65e10 (c = 54.60) | rising |
| relative gap (o1 − e1)/o1 | 0.99998 | — | ≥ 10⁻³ by 4.7 orders |
Truncation adequacy at both grid ends: the ratios move ≤ 5 % while N doubles (c = 30.570758: N = 240/320/420/520 → o1/e1 = 5.24/5.52/5.43/5.50e4,
e2/e1 = 1.50/1.58/1.60/1.61e9; c = 54.598: N = 500/650/800 → 1.76/1.79/1.74e5 and 1.65/1.67/1.66e10). The eigenvalues themselves converge slowly
and together (e1 at c = 30.57: 3.00/2.25/1.99/1.83e-153; Botha's rigorous U = 1.554e-153 from a larger basis), so the scan decides the ordering,
not the values. Botha's own float levels agree: his mu_ref (= his upper bound for λ_1^odd) at L = 2 is 3.24e-278; our N = 500/650/800 values
6.11/4.73/4.09e-278 extrapolate (o∞ + C N^{-2.23}, three-point fit) to ≈ 3.0e-278. Verdict of the scan: no candidate counterexample.

## 3. Audit of Pieces III–IV (code and records read; rh60 v1 and rh45 v1 fetched read-only from rstation.ai)
- **Piece III** ((1.19, 1.4035], λ² ∈ (10.8, 16.6]): the claim is covered independently by ../WeilSimple (λ² ∈ [2, 20]); his Elston + deflation
  certificate itself is not re-audited here (R415: unaudited beyond L = 1.19).
- **Piece IV on [1.4035, 1.7006]** (λ² ∈ [16.6, 30]): covered independently by WeilSimple + CCMFMM. The audit below concerns (1.7006, 2]:
  31 windows L_i = 1.71002197265625, …, 2 (data/piece4/glue: g_<L>_odd.json, g_<L>_l2.json, upper_<L>_{0,1}.jsonl).
- **Imported method.** rh60 v1 (bramble_trace, rStation, 26 Oct 2026) proves only positivity, inf σ(A_a) > 0 for a ≤ 2 (its Thm 1.1), by full-source
  Schur certificates: exact core elimination (displacement rank 2), interval assembly of the intermediate Schur matrix, a moment/Abel bound for the
  distant block, a structured backward-error LDL. rh45 v1 (lihongmin) supplies the Green-response identity and the bases. Botha's Theorem V is his
  own shifted (μ) and deflated (ϑ⟨Φ, ·⟩²) general-L version; supplement/PROOFS.md P1–P23 re-derives every constant (written 31 Oct with Claude,
  reviewed only by AI agents, no human check). PROOFS §11 lists what it does not prove: stage S's cell evaluation (Elston's schur_upper, the
  prime-norm bound s(L) of (V1)), the implementation of Elston's integer Cholesky (V2), the operator setup (closability, min-max), Arb itself.
- **Structure of one certificate** (record fields, L = 2 odd): core K = 512 with M = A_core − μ ≻ 0 by exact integer Cholesky (m_core = 1.79e-278,
  half the float λ_min(M) = 3.58e-278), certified core solves at prec_core = 3392 bits (solve_err_a = 1.5e-357: needed, cond(M) ~ 1e278);
  intermediate block [K, N), N = 69632, Schur complement S_input − δ_0 by Lemma SR (exact integers, ε_E = 6.4e-16, δ_0 = 35/256); far block
  q ≥ N by t′_q = H(κ_q) − μ − s − τ_N ≥ t′_N = 1.022 (s = 7.98, τ_N = 0.763); omitted tail r_M = 0.135 (M = 4N); final margin V5 = 0.00217.
  Every certified number is an Arb end or an exact rational (MAP.md). The size N is forced by the far block: H(κ_N) must exceed s + τ_N ≈ 8.7.
- **The levels μ are not fitted to the step test.** μ_odd = (2/3)·mu_ref with mu_ref = his own fr_sweep upper bound for λ_1^odd (par 1, K = 640:
  8.472e-149 at L_1 = 1.71002, 3.238e-278 at L = 2); μ_l2 = (1/4 … 1/1024)·lam2hat (first passing rung; 1/1024 at L = 1.89587, a weak but valid rung).
- **Step rule.** cor:step with lem:mono (= our chain rule, ../WeilSimple sec. 5): U(L_i) < min(μ_odd, μ_l2)(L_{i+1}) gives simple-even on [L_i, L_{i+1}].
  Recomputed in Arb from his records (compare2.py, last column; independent of his check_steps.py): every step in (1.7006, 2] passes, ratios 2.20–2.49 (min
  2.201 into L = 1.87665; his Piece IV minimum 2.0986 lies below 1.7006), plus 470 and 449 at the refinement windows 1.91913, 1.92365.
- **Rule E** (N raised by 1024 at 1.92364501953125 after the registered λ_2 rungs failed at V4) is post hoc, but not needed: λ_2^ev(1.92365) ≥
  λ_2^ev(1.93262) ≥ μ_l2(1.93262) = 1.302e-236 (lem:mono), so min(μ_odd(1.92365) = 1.232e-236, 1.302e-236)/U(1.91913) = 448 > 1; the next step
  uses only bounds at 1.93262. Closure valid (also with our U, sec. 4).
- **His certificate runs cannot be rerun under the 290 MB cap**: record field peak_gb = 1.49–3.43 GB per window (N = 12288–69632, 300–2630 s on a
  64-vCPU VM). We reran what fits: the record checker (sec. 3a) and an independent recomputation of every bound we can bound from the other side.

## 3a. Rerun of Botha's own scripts (under the cap)
- code/rh60route/check_steps.py data/piece4/glue data/piece4/windows.txt --stamp 7a3c6b28…b9e8: "ALL STEPS PASS (53 comparisons): simple-even
  on [1.4035, 2] (minimum step ratio 2.099)"; with --l2-from 1.92364501953125=1.9326171875 (no rule E): the same line; 35 MB, < 15 s
  (logs/botha_check_steps*.out). This validates his records against his step rule; it recomputes no window (his own caveat).
- certify.py / sr.c (Theorem V) and fr_sweep.py (K = 640, prec 4000: t = 1125 s and > 225 MB for the matrix alone) do not fit the cap:
  not rerun. In their place, sec. 4 recomputes the upper side independently and tests every lower bound against it.

## 4. Independent upper side: our rigorous U replaces fr_sweep on all of (1.7006, 2] (point.py, compare2.py; proofs FULLPROOF.md secs. 1–3)
At every grid point point.py evaluates, in Arb from the generator balls (no rounded matrix), the Rayleigh quotient of deig's exact-dyadic even
vector y_1 (U ≥ λ_1^even), the 2×2 Rayleigh–Ritz pencil on span(y_1, y_2) (U2 ≥ λ_2^even) and the odd quotient (Uo ≥ λ_1^odd), at c_q < e^{2L}
(checked in Arb), so by monotonicity they bound the eigenvalues at L itself. Peaks: ≤ 44 MB per point (deig ≤ 41 MB), ≤ 50 s per point.
(a) Consistency of every Botha lower bound with our upper bounds (a ratio < 1 would refute his record): Uo/μ_odd ∈ [2.20, 3.11], U2/μ_l2 ∈ [3.37, 6.39]
    (1365 at L = 1.89587, his weak 1/1024 rung), all 31 windows: no contradiction. Our U exceeds his U(L) by 1.43–2.15 (smaller basis: N ≤ 600).
(b) Hybrid step test, his μ with OUR U: U_ours(c_q(L_i)) < min(μ_odd, μ_l2)(L_{i+1}) holds for all 30 consecutive pairs in (1.7006, 2], ratios
    1.071 (L_i = 1.83643, N = 320) … 1.641, and 261, 266 into the refinement windows; first step from c = 30 (CCMFMM's certified U(30) =
    3.27353749418e-150) into L_1 = 1.71002: ratio 17.25. Hence, by cor:step (FULLPROOF sec. 3), simple-even holds on all of
    L ∈ [1.7006, 2] (λ² ∈ [30, e⁴]) with fr_sweep removed from the trust base; Botha's Theorem V lower bounds μ remain imported.
    Rule E: the step into 1.92365 needs μ_l2(1.92365) > U_ours(1.91913) = 4.642e-239; the rule-E record (4.572e-232) or, without rule E,
    μ_l2(1.93262) = 1.302e-236 by lem:mono (ratio 280) both pass, so rule E is not needed here either.
(c) Reviewer (sec. 5): GAP only in the wording "any vectors work" (Lemma G needs y_1, z_1 ≠ 0 and y_1, y_2 independent; verified by Arb at every
    point, FULLPROOF sec. 2); the hybrid step itself (one nonzero vector) CORRECT, including the rounding direction of c_q and the j = 0 row.

## 5. Our own certificate beyond λ² = 30: CCMFMM's chain continued (logs/chain.jsonl; method FULLPROOF sec. 5)
Same rules as ../CCMFMM (both certificates at ε = 2.00001·U(c_i), deflated even and odd sectors, Lemma B′ backward error, proved node scan),
started from its certified endpoint U(30) = 3.27353749418e-150:
| c_i → c_{i+1} | U(c_i) | ε = L(c_{i+1}) | U(c_{i+1}) | certified ν (defl, odd) | n (defl / odd) | NG | bits | ok | peak RSS MB (g; prep, LDL defl; odd) | s |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 → 92/3 | 3.27353749418e-150 | 6.547108e-150 | 6.75090730136e-154 | 6.547108e-149, 6.547108e-151 | 3332 / 3331 | 320 | 896 | ok | 111; [83, 16]; [72, 15] | 361.6 |
| 92/3 → 251/8 | 6.75090730136e-154 | 1.350188e-153 | 1.00300264109e-157 | 1.350188e-154, 1.350188e-154 | 4450 / 4449 | 320 | 960 | ok | 115; [None, 34]; [100, 32] | 611.1 |
| 251/8 → 385/12 | 1.00300264109e-157 | 2.006015e-157 | 1.72477599362e-161 | 2.006015e-158, 2.006015e-158 | 4763 / 4762 | 320 | 960 | ok | 106; [99, 31]; [95, 28] | 975.8 |

covered: [30, 385/12] = [30, 32.083333]; ok steps 3, failed attempts 0; max RSS 115 MB
Hence (Courant–Fischer + monotonicity, ../WeilSimple sec. 5): for every λ² ∈ [30, 385/12] = [30, 32.0833] (L ∈ [1.7006, 1.7341]) QW_λ has a simple,
even ground state; with WeilSimple and CCMFMM, λ² ∈ [2, 32.0833]. Steps of 0.67–0.71 in c; 6–16 min each on the shared machine. The chain was stopped
at 20:01 UTC for time during the step from 385/12. Single windows λ² = 47, 54: not attempted (time; CCMFMM's λ² = 40 window took 1783 s, and the
dense gstate block at NG ≈ 550–750 would approach the cap).

## 6. Botha's stage-S input s(L) (V1, imported), falsification test (scheck.py, logs/scheck.out)
(V1) needs s ≥ sup ⟨C_b h, h⟩/‖h‖² for C_b = Σ_{log n < 2b} Λ(n) n^{-1/2} (τ_{log n} + τ_{−log n}) on L²(−b, b) (= minus the prime part of Q). Our
Galerkin values (top eigenvalue of 2·PR_N on our even and odd bases, N = 200/400/800, double precision) rise with N and stay below his s at all
31 windows, margin s − Ritz₈₀₀ ≥ 0.0044 (L = 1.71002: 5.63092/5.63155/5.63187 vs s = 5.63654; L = 2: 7.97326/7.97400/7.97429 vs s = 7.97984).
No violation; the near-agreement (0.07 %) also confirms the normalisation. This is a consistency test, not a proof of (V1).

## 7. Reviewer (one reviewer subagent, on the most delicate new step: Lemma G + the hybrid step, values pasted inline)
Verdict GAP, on wording only: Lemma G's "any vectors work" needs y_1, z_1 ≠ 0 and y_1, y_2 independent (else det(H − tG) ≡ 0); closed in
FULLPROOF sec. 2 (Arb would return an infinite or NaN ball otherwise; every bound printed is finite). Everything else CORRECT: exact-dyadic vectors
from floats are legitimate trial vectors; the pencil root bounds λ_2^even; the j = 0 even row matches the closed form; rounding c_q down is the
right direction; the hybrid step uses only the one-vector U and proves simple-even on [log(30)/2, 2] given Botha's μ.

## 8. Does Botha's Theorem 1 hold on L ∈ (1.7006, 2]?
- As a claim: yes, nothing contradicts it. No crossing or near-degeneracy (sec. 2); every one of his 62 lower bounds in the range lies below our
  rigorous upper bound (sec. 4a); his step test reruns (sec. 3a) and survives replacing his upper bounds by ours (sec. 4b); rule E is not needed.
- Independently certified (our lower and upper bounds only): λ² ∈ [2, 385/12 (= 32.0833)] (WeilSimple [2, 20], CCMFMM [20, 30], this round [30, 385/12 (= 32.0833)]),
  plus CCMFMM's windows λ² = 35, 40.
- Hybrid (our U, his μ): all of λ² ∈ [30, e⁴]; fr_sweep is no longer needed anywhere in (1.7006, 2].
- Still imported trust on λ² ∈ (385/12 (= 32.0833), e⁴] (except 35, 40): his lower bounds μ, i.e. certify.py and sr.c implementing Theorem V (PROOFS
  P1–P23, AI-written and AI-reviewed), stage S's s(L) (Elston's schur_upper; our Galerkin test, sec. 6, finds no violation), Elston's integer
  Cholesky for (V2), and the operator setup.

## 9. Cost and memory (one process at a time, single-threaded, background priority; every process under memcap 290 MB)
| job | peak (maxrss / memcap) | wall time |
|---|---|---|
| point.py per window (gen + deig even/odd + Arb Rayleigh quotients), N = 320–600 | ≤ 44 MB (deig ≤ 41 MB) | 12–50 s |
| truncation study c = 54.598, N = 800 (deig) | 69 MB | 48 s per sector |
| CCMFMM chain step (gstate NG = 320; prep; slt LDL, n ≈ 3300–4500) | gstate 111, prep 72–83, LDL 15–16 MB | 6–15 min |
| Botha's check_steps.py (both variants) | 35 MB | < 15 s |
| Botha's certify.py per window (his record field peak_gb) | 1.47–3.43 GB | not run (cap) |

## 10. Files (Riemann/Deep/BothaAudit/)
- Scan and upper bounds: scan.py (FastWeilMat, Python Arb cross-check), gen.py (generators + entrywise test), deig.c (float LDL + inverse
  iteration; build: cc -O2 -o deig deig.c -lflint -lgmp -lmpfr -lm), post.py, one.sh, scan_grid.sh, mkgrid.py, point.py (rigorous U, U2, Uo),
  scan2.sh, compare2.py (his records vs ours; hybrid steps; his step ratios), falsify.py, scheck.py (sec. 6), weil_matrix.py (copy of CCMFMM's).
- Chain: ccm/ (CCMFMM scripts, FULLPROOF sec. 5), chain.sh, table_chain.py; binaries and runs lived in .scratch/bothaaudit (deleted).
- compare2.py and scheck.py read Botha's zip unpacked to .scratch/bothaaudit/b (fetch: Zenodo API record 23852617; deleted after the run).
- Logs: logs/points.jsonl (31 windows), logs/points_mem.txt, logs/compare2.out, logs/falsify.out, logs/conv.jsonl, conv30.jsonl, scan_N240.jsonl,
  scan_float.jsonl (float-only first pass), val.jsonl, logs/chain.jsonl / chain.out, logs/botha_check_steps*.out, logs/scheck.out, logs/grid.txt.

RH relevance: N13
