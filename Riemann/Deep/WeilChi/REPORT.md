RECORD: first certified Weil-positivity windows for Dirichlet L-functions. Fully rigorous (Arb plus a proved scan): for all f supported in [-a, a], Q_chi(f) >= 2^-200 ||f||^2 for chi_-4 (mod 4) with a = (ln 47)/2 = 1.925 (both parities), and for chi_-3 (mod 3) with a = (ln 35)/2 = 1.778. Resting on a float-scan error budget (not proved): chi_-3 even f to a = 1.925 (2^-275), and chi_5, chi_8 (mod 5, 8) to a = 1.844 (2^-140, 2^-91). Prime-free windows are rigorous with eps = 1/16.

# WeilChi: certified Weil positivity windows for Dirichlet L-functions (weilchi, 2026-11-05)

RH relevance: N1 (L-function branch). For a primitive Dirichlet character chi, GRH for L(s, chi) holds if and only if Weil's quadratic form Q_chi(f) is >= 0 for every test function f. Here Q_chi(f) >= eps ||f||^2 is proved for all f supported in a window [-a, a]. No such certificate for a single Dirichlet L-function was found in print (sec. 8).

Review (2026-11-05): GAP and ROUTINE. Both gaps are fixed. (1) The explicit formula is now stated for smooth test functions, and the identity the certificate uses is proved directly for every f in L^2 (DERIVATION.md sec. 1 and 3b). (2) The node scan is now proved: rscan.c with the a-priori bound in RIGOR.md. The chi_-3 certificates at c = 30 and 35 and the chi_-4 certificate at c = 47 were re-run with it. ROUTINE stands: this is a direct transfer of Chekmachev's method; what is new is the record type, not the method.

## 1. Result
**Theorem (computer-assisted, Arb ball arithmetic).** Let chi be one of the characters below, a its half-width, and f any complex function in L^2(-a, a) (in particular any C^1 function supported in [-a, a]). Then Q_chi(f) >= eps ||f||^2. The row "even f only" holds for every even f.

| chi | q | parity | half-width a (c = e^{2a}) | eps | margin, even f | margin, odd f | node scan |
|---|---|---|---|---|---|---|---|
| chi_-4 | 4 | odd | (ln 47)/2 = 1.9252 | 2^-200 = 6.2e-61 | 3.388e-61 | 1.487e-55 | rigorous (re-run) |
| chi_-3 | 3 | odd | (ln 35)/2 = 1.7777 | 2^-200 = 6.2e-61 | 4.261e-59 | 9.237e-55 | rigorous (re-run) |
| chi_-3 | 3 | odd | (ln 30)/2 = 1.7006 | 2^-170 = 6.7e-52 | 2.372e-50 | 3.808e-46 | rigorous (re-run) |
| chi_-3, even f only | 3 | odd | (ln 47)/2 = 1.9252 | 2^-275 = 1.6e-83 | 2.120e-79 | not certified (sec. 5) | float budget |
| chi_5 | 5 | even | (ln 40)/2 = 1.8444 | 2^-140 = 7.2e-43 | 2.964e-41 | 4.255e-37 | float budget |
| chi_8 | 8 | even | (ln 40)/2 = 1.8444 | 2^-91 = 4.0e-28 | 2.026e-25 | 1.774e-21 | float budget |
| chi_-4, prime-free | 4 | odd | (ln 3)/2 = 0.5493 | 1/16 | 0.0255 | 0.808 | none needed |
| chi_-3, prime-free | 3 | odd | (ln 2)/2 = 0.3466 | 1/16 | 0.364 | 0.885 | none needed |
| chi_5, prime-free | 5 | even | (ln 2)/2 = 0.3466 | 1/16 | no negative part (P - N = P > 0) | same | none needed |
| chi_8, prime-free | 8 | even | (ln 3)/2 = 0.5493 | 1/16 | 0.282 | no negative part | none needed |

"Node scan" says how K and E were proved to contain every node below the thresholds (sec. 4).
- **rigorous**: the C scan rscan.c with a proved one-sided bound mu >= mu_c - E (RIGOR.md), which is all the no-miss argument needs; re-run on 2026-11-05 after review.
- **float budget**: the first float64 scan with libm cos; its error budget is an estimate, not a proof.
- **none needed**: kscan = K2, so every node below kscan lies in K or E anyway.
The "rigorous (re-run)" rows reproduce the earlier margins to all printed digits, with identical node sets.
"Margin" is 1 - lambda_max(B^T P^{-1} B) of the dual check (sec. 4) at the stated eps; it is proved > 0 by Lemma 5.5 (Cholesky of the midpoint, Arb residual, row-sum bound). A function supported in a smaller window is supported in [-a, a] and Q_chi is the same functional there, so each bound holds on every smaller window too.
Prime-free windows: no chi(n) Lambda(n) with n < e^{2a} is nonzero (chi_-4(2) = chi_8(2) = 0; chi_-3(2) = chi_5(2) = -1), so there Q_chi is purely archimedean, (1/2 pi) int [log(q/pi) + Re psi(1/4 + a_chi/2 + it/2)] |F|^2. Its symbol is negative near t = 0 (log(4/pi) + psi(3/4) = -0.844 for chi_-4), so even this needs a proof; the certificate gives it.

Rigorous re-runs (logs/rig_*.log, node sets in logs/scan_*_rig.json): same node sets and the same margins to all printed digits. Times and peaks: chi_-3 c = 30: 30 s and 28 s, 0.17 GB; chi_-3 c = 35: 235 s and 222 s, 0.49 GB; chi_-4 c = 47: 728 s (including the 56 s scan) and 708 s, 0.98 GB and 0.90 GB.
Run data of the first runs (logs/*.log hold one JSON line per run; table.py prints them all):
| run | K2 | t1 | bits | kscan (odd/even) | selected odd/even | K+ / K- / J- / G | order | max radius | time | peak RSS |
|---|---|---|---|---|---|---|---|---|---|---|
| chi_-4, c = 47, even | 201 | 1 | 448 | 8.1e8 / 3.0e8 | 794 / 253 | 598 / 824 / 114 / 43649 | 938 | 6.2e-128 | 472 s | not measured |
| chi_-4, c = 47, odd | 201 | 1 | 448 | same | same | 598 / 824 / 114 / 43648 | 938 | 6.6e-126 | 401 s | not measured |
| chi_-3, c = 47, even | 1001 | 1 | 448 | 1.5e9 / 5.5e8 | 968 / 334 | 1124 / 1043 / 152 / 56067 | 1195 | 2.0e-125 | 1510 s | not measured (old code; see sec. 11) |
| chi_-3, c = 47, odd (not certified) | 201 | 1 | 320 | same | 1015 / 364 | 815 / 1043 / 152 / 56066 | 1195 | 4.9e-86 | 1064 s | 1.34 GB |
| chi_-3, c = 30, even | 201 | 1 | 384 | 9.2e6 / 3.4e6 | 224 / 63 | 211 / 250 / 28 / 11386 | 278 | 1.4e-109 | 23 s | 167 MB |
| chi_-3, c = 30, odd | 201 | 1 | 384 | same | same | 211 / 250 / 28 / 11385 | 278 | 7.4e-108 | 29 s | 167 MB |
| chi_-3, c = 35, even | 201 | 1 | 384 | 4.2e7 / 1.6e7 | 721 / 31 | 147 / 749 / 19 / 36937 | 768 | 5.2e-108 | 260 s | 485 MB |
| chi_-3, c = 35, odd | 201 | 1 | 384 | same | same | 147 / 749 / 19 / 36936 | 768 | 1.4e-106 | 220 s | 486 MB |
| chi_5, c = 40, even | 601 | 1 | 384 | 1.2e8 / 4.5e7 | 311 / 96 | 469 / 347 / 44 / 18268 | 391 | 2.8e-108 | 110 s | 408 MB |
| chi_5, c = 40, odd | 601 | 1 | 384 | same | same | 469 / 347 / 44 / 18267 | 391 | 2.5e-107 | 186 s | 256 MB |
| chi_8, c = 40, even | 601 | 1 | 384 | 3.9e7 / 1.4e7 | 195 / 47 | 373 / 229 / 28 / 11105 | 257 | 6.1e-109 | 48 s | 218 MB |
| chi_8, c = 40, odd | 601 | 1 | 384 | same | same | 373 / 229 / 28 / 11104 | 257 | 1.9e-107 | 46 s | 171 MB |

## 2. The form (derivation in DERIVATION.md)
chi primitive mod q, a_chi = 0 (even chi) or 1 (odd chi), F(t) = int f(x) e^{itx} dx, g(x) = int f(y) conj f(y - x) dy (|F|^2 is the transform of g):
  Q_chi(f) = (1/2 pi) int [log(q/pi) + Re psi(1/4 + a_chi/2 + it/2)] |F(t)|^2 dt - 2 sum_{n < e^{2a}} Lambda(n) n^{-1/2} Re[chi(n) g(log n)].
There are no pole terms (L(s, chi) is entire). By the explicit formula, Q_chi(f) = sum over the zeros rho = 1/2 + i gamma of F(gamma) conj F(conj gamma); this is sum |F(gamma)|^2 >= 0 under GRH, and an off-line zero pair makes it negative for a suitable f.
For real chi the prime operator sum chi(n) Lambda(n) n^{-1/2} (S_{log n} + S_{-log n}) is real symmetric and commutes with x -> -x. So Q_chi splits into even and odd f with no cross term, and complex f = u + iv gives Q(u) + Q(v).
Lattice identity (Chekmachev's Thm 3.1, re-derived for chi). With L = ln c, I = [-L/2, L/2] and r_k = pi k / L, for every f in L^2(I):
  Q_chi(f) = (1/2L) sum_k mu(k) |F(r_k)|^2,
  mu(k) = log(q/pi) + Re psi(1/4 + a_chi/2 + i r_k/2) + 2 (-1)^k sum_{n>=0} c^{-b_n} b_n/(b_n^2 + r_k^2) - 2 sum_{m<c} chi(m) Lambda(m) m^{-1/2} cos(r_k log m),  b_n = 2n + 1/2 + a_chi.
Compared with zeta: no rank-one pole term, b_n = 2n + 3/2 for odd chi, the constant log(q/pi), and signed prime weights.

## 3. Numerical check against zeros (ef_check.py, logs/ef_check.log)
Test functions: f = (1 - (x/A)^2)^8 (even) and (x/A)(1 - (x/A)^2)^8 (odd) on [-1.2, 1.2], so the prime powers below e^2.4 = 11.02 enter. The zeros of L(s, chi) up to t = 90 were found as sign changes of the real function e^{i theta(t)} L(1/2 + it, chi) (mpmath, 30 digits): 44 zeros for chi_-4 (first 6.020948905) and 40 for chi_-3 (first 8.039737156).
| chi | f | sum over zeros | Q_chi by the formula | difference | lattice identity minus formula |
|---|---|---|---|---|---|
| chi_-4 | even | 0.0533620888388890646 | 0.0533620888388890646 | 3.1e-22 | 6.2e-32 |
| chi_-4 | odd | 0.0106837681484936486 | 0.0106837681484936486 | 6.7e-22 | 4.3e-32 |
| chi_-3 | even | 0.00331707128298722211 | 0.00331707128298722211 | 3.4e-22 | 1.1e-31 |
| chi_-3 | odd | 0.00177746443416940006 | 0.00177746443416940006 | 5.6e-22 | 4.0e-32 |
The zero sums match the formula to 7e-22, the size of the zeros omitted above t = 90. The lattice identity matches to 1e-31. A smoke run (A = 0.6, zeros to t = 40) matched to 4e-10, again the truncation size.

## 4. Method: Chekmachev's periodized lattice minorant, ported (cert_chi.py)
cert_chi.py ports our audited re-implementation (../ChekAudit/evidence/cert2.py) of V. Chekmachev, Zenodo 23848597 (Theorem 5.10 with the rule of Definition 6.1). Changes for chi: the symbol mu(k) above; no pole column in either sector; prime weights chi(m) Lambda(m)/sqrt m; scan bound with 2 sigma = 2 sum |chi(m)| Lambda(m)/sqrt m.
- Two node systems turn the lattice sum into (1/2) sum_j mu(2j) <f, phi_j>^2 + (1/L) sum_{k odd} mu(k) <f, e_k>^2. The even nodes 2j give Neumann cosines (even f) or Dirichlet sines (odd f); the odd nodes k give Dirichlet cosines or Neumann sines.
- Minorant: every odd k outside a finite set K has mu(k) >= t1 (here t1 = 1) and is replaced by t1, using Parseval for the odd system. Every even j outside E has mu(2j) >= 0. The resulting P - N, with P > 0 (diagonal plus low rank) and N = B B^T of finite rank, is >= 0 iff lambda_max(B^T P^{-1} B) <= 1. This is computed by Woodbury, with all inner products in closed form and the sums over even nodes by partial fractions.
- K and E contain all odd k <= K2 and all j <= K2. They also contain every node found by a scan up to kscan, beyond which the monotone bound h(r) - 2 sigma - (tail) >= t1 (odd) or >= 0 (even) holds; h(r) = log(q/pi) + Re psi(1/4 + a_chi/2 + ir/2) is increasing in r > 0. Chekmachev's neighbour rule adds 65 even nodes around each selected odd node, and odd nodes 2j +- 1 and even nodes j +- 3 around each negative even node.
- Scan, rigorous version (rscan.c, proof in RIGOR.md; used for the rows marked rigorous). For every node below kscan it computes mu_c(k) with the proved one-sided bound mu(k) >= mu_c(k) - E (no absolute bound is claimed; the block bounds bias mu_c downward). The inputs are exact 64-bit integer phases, a 2^16-entry cosine table built in Arb, a two-term Taylor step for the low 48 phase bits, and Arb lower bounds for h and upper bounds for the tail on blocks of nodes (h increases, the tail decreases). Only IEEE +, -, * are used, with no libm calls and no FMA contraction. E = 2.1e-11 (chi_-3, c = 30), 1.1e-10 (chi_-3, c = 35) and 2.6e-9 (chi_-4, c = 47), and nodes with mu_c < threshold + 1e-6 are re-checked in Arb. So K and E provably contain every node below the thresholds. Scan time: seconds (chi_-3, c = 35) to about two minutes (chi_-4, c = 47).
- Scan, float version (fscan/cscan; the first runs, rows marked float budget): a float64 pre-scan with exact 64-bit phases (C, cscan.c; rotators re-synchronised from exact phases every 1024 nodes). It flags nodes at margin 1e-6; the float error budget is below 1e-8 (phase error <= k 2^-65 <= 4e-11 cycles for k <= 1.5e9, times 2 pi sigma <= 2 pi 10; Stirling remainder < 1e-25; rounding and rotator drift ~1e-12). Every flagged node is re-evaluated in Arb (96 bits). Here the error budget is an estimate, not a proof: libm cos has no proved bound. At chi_-4, c = 47 the C scan selects exactly the nodes of the independent numpy scan (794 odd, 253 even; logs/scan_check.log). The rigorous scan selects the same node sets wherever both were run (chi_-3 at c = 30 and 35: 224/63 and 721/31; chi_-4 at c = 47: 794/253).
- All the rest uses Arb balls (python-flint 0.9.0). Positivity of Id - W is proved by Lemma 5.5 (DERIVATION.md sec. 5). Take an estimate lambda of the smallest eigenvalue, Cholesky-factor the midpoint minus s0 Id, enclose the residual in Arb, and check its row sums against s0 - nu.
- The code was revised during the session to save memory (sec. 11): Gram blocks built directly, Woodbury and residual in blocks, and finally Rayleigh-Ritz estimates of lambda with s0 = lambda/2. Every revision reproduces the earlier margins to all printed digits (logs/leancheck.log: chi_-4 at c = 30, eps = 0, and at c = 16, eps = 2.4e-18). The chi_-4 (c = 47), chi_-3 (c = 47, even f) and prime-free certificates ran on the first version. The chi_5, chi_8 and chi_-3 (c = 30, 35) certificates ran on later ones.

## 5. Cross-checks
- **Zeros (sec. 3).** The form, its signs (chi(n) against conj chi(n)), the parity factor Gamma((s + a_chi)/2), the conductor term and the lattice symbol all agree with sums over zeros.
- **Independent Galerkin upper bounds against the certificate** (galerkin_chi.py: x-space closed forms, Neumann cosine basis, even f; logs/galerkin.log, logs/bracket.log).

  | chi, c = 16 | lambda_min(Q_N), N = 60 / 120 / 200 (upper bounds) | certificate at eps above the N = 200 value | certificate at eps = half of it |
  |---|---|---|---|
  | chi_-4 | 5.701e-18 / 5.029e-18 / 4.742e-18 | eps = 4.8e-18: margin -1.6e-18, fails, as it must | eps = 2.4e-18: margin 5.2e-18, certified |
  | chi_-3 | 6.366e-25 / 5.739e-25 / 5.328e-25 | eps = 5.4e-25: margin -2.9e-25, fails, as it must | eps = 2.6e-25: margin 5.8e-25, certified |

  So 2.4e-18 <= lambda_min <= 4.74e-18 (chi_-4) and 2.6e-25 <= lambda_min <= 5.33e-25 (chi_-3) at c = 16: two independent codes bracket the true value within a factor 2, and the certificate refuses an eps above a proven upper bound.
- **Scan.** C and numpy scans agree exactly (sec. 4).
- **Memory-lean rewrite.** After the blockwise Woodbury and residual rewrite (sec. 11), the chi_-4, c = 30 margin is unchanged to all 12 printed digits (2.23133452614e-36; logs/leancheck.log).
- **The method can fail, and then says so.** chi_-3 at c = 40, even f, K2 = 201: margin -2.6e-68, not certified. With K2 = 601 it gives +3.1e-67 (logs/m3_40_tests.log). The minorant had too little slack; Q itself was not negative. Likewise t1 = 1/2 with 8 neighbours at chi_-3, c = 47, odd f gave -2.2e-75.
- **chi_-3 at c = 47, odd f, is not certified.** With K2 = 201 and default knobs, Rayleigh-Ritz on 6 block inverse-iteration vectors finds two negative eigenvalues of the midpoint of Id - W: -1.8e-41 and -2.7e-54 (logs/m3_47_cert.log). So the minorant itself is indefinite there, and Lemma 5.5 rightly fails. An earlier single-vector inverse iteration had reported only 1.75e-75, the eigenvalue closest to 0. Larger K2 (as for the even f, K2 = 1001) should repair it, but in this implementation the Woodbury solve then needs more than the 1.5 GB allowed (sec. 11). The same memory wall stopped chi_-3 at c = 40 (K2 = 601).
- Every certificate in sec. 1 uses Chekmachev's default knobs (t1 = 1, 65 neighbours) and K2 = 201 to 1001. Lemma 5.5 is a proof whatever the eigenvalue estimate: it shows Id - W >= s0 - ||E|| >= nu > 0 directly.

## 6. How the slack decays: consistent with exp(-4 pi c / q)
Margin at eps = 0 (Delta, even f), with these K2: chi_-4 uses max(201, ceil(5c^2/3)) for c <= 25 and 201 from c = 30; chi_-3 uses 601 at c = 16 and 40, 201 at c = 20 to 35, and 1001 at 47; chi_5 and chi_8 use 601:
| c | chi_-4 | chi_-3 | chi_5 | chi_8 |
|---|---|---|---|---|
| 2 or 3 (prime-free) | 0.236 | 0.503 | P > 0 | P > 0 |
| 5 | 6.3e-4 | | | |
| 9 / 10 | 6.0e-9 (c = 9) | | 4.9e-9 | 2.4e-5 |
| 16 | 1.2e-17 | 1.4e-24 | | |
| 20 | 5.2e-23 | 1.8e-31 | 2.3e-19 | 4.5e-12 |
| 25 | 4.3e-30 | 4.7e-41 | | |
| 30 | 2.2e-36 | 2.4e-50 | | |
| 35 | | ~4.3e-59 (margin at eps = 2^-200) | | |
| 40 | 7.3e-50 | 3.1e-67 (K2 = 601) | ~3.0e-41 (margin at eps = 2^-140) | ~2.0e-25 (margin at eps = 2^-91) |
| 47 | 2.6e-60 | ~2.1e-79 (at eps = 2^-275) | | |
Fitted -d(ln Delta)/dc: 3.18 (chi_-4, c = 20 to 47; 4 pi/4 = 3.14), 4.10 (chi_-3, c = 20 to 47; 4 pi/3 = 4.19), 2.47 (chi_5, c = 10 to 40; 4 pi/5 = 2.51), 1.54 (chi_8, c = 10 to 40; 4 pi/8 = 1.57). So Delta ~ exp(-4 pi c/q), the shape of Elston's empirical law lambda_min ~ exp(-2 T*) with T* = 2 pi x/q (Zenodo 23174403). For zeta the rate is 4 pi c, so chi mod q at window c is about as hard as zeta at c/q; the window here is limited by the scan (about exp(2 sigma(c)) nodes) and by memory, not by precision.
The ratio of the margin's slope in eps (kappa = (Delta - margin)/eps) is about 3.6 at chi_-4, c = 47 (even), so the largest certifiable eps there is about Delta/3.6 = 7e-61. The certified 2^-200 = 6.2e-61 is close to it. For chi_-3 at c = 47 (even f) we chose eps far below Delta (Delta/13000) to need one run.

## 7. Replay
From Riemann/Deep/WeilChi, Python packages only via uv, one process, single-threaded:
    cc -O2 -o /tmp/cscan cscan.c -lm
    cc -O2 -ffp-contract=off -o /tmp/rscan rscan.c
    export WEILCHI_RSCAN=/tmp/rscan WEILCHI_SCAN=/tmp OMP_NUM_THREADS=1     # rigorous scan; WEILCHI_CSCAN=/tmp/cscan instead gives the float scan
    uv run --no-project --with mpmath --with numpy --with python-flint python3 cert_chi.py m4 47 even 2^-200 448 201     # and: odd
    uv run --no-project --with mpmath --with numpy --with python-flint python3 cert_chi.py m3 47 even 2^-275 448 1001    # chi_-3, even f (> 1.5 GB with the final code? peak not measured)
    uv run --no-project --with mpmath --with numpy --with python-flint python3 cert_chi.py m3 30 even 2^-170 384 201     # chi_-3, all f at c = 30; and odd
    uv run --no-project --with mpmath --with numpy --with python-flint python3 cert_chi.py m3 35 even 2^-200 384 201     # chi_-3, all f at c = 35; and odd
    uv run --no-project --with mpmath --with numpy --with python-flint python3 cert_chi.py p5 40 even 2^-140 384 601     # and odd; p8: 2^-91
    uv run --no-project --with mpmath --with numpy --with python-flint python3 cert_chi.py m4 3 even 2^-4 128            # prime-free windows
    uv run --no-project --with mpmath --with numpy --with python-flint python3 ef_check.py m4 1.2 8 90                    # zeros; and m3
    uv run --no-project --with mpmath --with numpy --with python-flint python3 galerkin_chi.py m4 16 60,120,200 256
Arguments of cert_chi.py: character (m4 = chi_-4, m3 = chi_-3, p5 = chi_5, p8 = chi_8), c, parity of f, eps, bits, K2. Each run prints one JSON line; "rigorous": true means Lemma 5.5 succeeded at that eps. Without WEILCHI_CSCAN the numpy scan is used (6 times slower). Optional knobs (environment): WEILCHI_T1, WEILCHI_NBO, WEILCHI_NBEK, WEILCHI_NBEJ.

## 8. Novelty pass (2026-11-05)
- **Zenodo** (API, 9 queries: "Weil positivity", "Weil positivity Dirichlet", "Weil criterion", "Weil quadratic form L-function", "Weil explicit formula Dirichlet character positivity", "Connes positivity Dirichlet", "GRH Weil positivity", "Weil positivity character", "Dirichlet L-function Weil"). Closest hit: Elston, Zenodo 23174403 (v1.1.0, 6 Oct 2026), "Weil positivity on finite windows across L-functions ... and an unconditional bound for Dirichlet L-functions". We read its description and README and listed its GitHub repository (826 files, docs/EPSTEIN.md, DECAY_LAW.md, CONDUCTOR5_FAMILY.md, CONTROL_CERTIFICATES.md). Its certified lower bounds cover:
  - zeta;
  - the Dedekind zeta of Q(sqrt -5), zeta * L(chi_-20), on [-1.495, 1.495]. That form is the sum Q_zeta + Q_chi_-20 and certifies neither factor;
  - three control functions without an Euler product (docs/CONTROL_CERTIFICATES.md): F_{t*} (conductor 5, gamma factor Gamma(s/2), -F'/F coefficients (1 - t*) chi_5(n) + t*(1 + sqrt5 [5|n])), the Davenport-Heilbronn function (conductor 5, Gamma((s+1)/2)) and the Epstein zeta Z_1. These have positive windows up to a = 0.835 to 1.19 and certified negative witnesses beyond (F_{t*} turns negative at x = 5.3625, a = 0.84).
  For Dirichlet L-functions themselves it gives upper bounds only (a Hermite rate; finite-basis numerics, "nothing here is certified"). The contrast is instructive: L(s, chi_5) stays certified-positive to a = 1.844 (x = 40), far past the x = 5.36 where the non-Euler F_{t*}, with the same gamma factor and conductor, fails. Other hits do not certify a Dirichlet-L window: 23769666 (archimedean threshold for odd f, zeta), 23691133 (conductor-23 CM form: screw-function values, not a window), 23005653 (a GRH criterion), 23848597 and 23978194 (Chekmachev, zeta).
- **arXiv**: monthly listings math.NT 2026-10 (854 titles) and 2026-11 (140), math.CA and math.SP for both months, scanned for Weil, positivity, explicit formula, Li, Bombieri, Connes, Dirichlet L. API searches (5 queries) found 2608.24827 (Zhu; abstract: zeta only, window 0.8), 2006.13771 (Connes-Consani, archimedean place, zeta), 2606.09096 (screw function, zeta) and 2503.15832 (positivity technique under GRH, not a window).
- **Ours**: R335 (TREE.md) found none; NOGO.md has no Dirichlet-L window.
- Verdict: no certified Weil-positivity window for a single Dirichlet L-function found. The claimed record type is new; the method (Chekmachev's) is not.

## 9. Corollary (routine): Dedekind zeta functions of Q(i) and Q(sqrt -3)
For K = Q(i) (chi_-4) and K = Q(sqrt -3) (chi_-3), zeta_K = zeta * L(s, chi) and the Weil form of zeta_K is Q_zeta + Q_chi. Our own zeta certificate gives Q_zeta >= 2^-430 ||f||^2 on [-1.6, 1.6] (../WeilRecord/REPORT5.md). Hence Q_{zeta_K} >= (2^-430 + eps_chi) ||f||^2 there, all from our own certificates. Add Chekmachev's zeta certificate at c = 47 (audited sound in ../ChekAudit; the final run used his unpublished code). Then the same holds for Q(i) on [-1.925, 1.925], and for Q(sqrt -3) on [-1.778, 1.778] (all f) and on [-1.925, 1.925] (even f). Elston certified zeta_K for K = Q(sqrt -5) on [-1.495, 1.495]. This is routine (COROLLARIES.md), not a new record type.

## 10. Self-review (two passes)
Pass 1, mathematics:
- Explicit formula. The n^{-s} terms of -L'/L(s, chi) carry chi(n) and pair with g(log n); the reflected line carries conj chi(n) and pairs with g(-log n). For real chi the term is 2 chi(n) g(log n). The archimedean term comes from Gamma((s + a_chi)/2) and the conductor from (q/pi)^{s/2}. There are no poles. All of this is confirmed against zeros for both parities of f (sec. 3); a wrong sign or parity would show at the 1e-2 level, not 1e-22.
- Factor 2 between f and g = f * f~: (1/2 pi) int |F|^2 = ||f||^2 = g(0), and sum_k |F(r_k)|^2 = 2L ||f||^2. Both enter the lattice identity, which matches the direct formula to 1e-31.
- Periodization: int_{-L}^{L} e^{-b|x|} e^{-i r_k x} dx = 2b(1 - (-1)^k c^{-b})/(b^2 + r_k^2) with b = b_n = 2n + 1/2 + a_chi. The tail remainder is bounded by c^{-b_N}/(b_N (1 - c^{-2})).
- Scan bound: d/dy Re psi(x + iy) = sum_n 2(n + x)y/((n + x)^2 + y^2)^2 > 0, so h increases for r > 0; b/(b^2 + r^2) decreases; and |sum chi(m) Lambda(m) m^{-1/2} cos| <= sigma. So monotone(k) bounds mu(k') for every k' >= k of either parity.
- Minorant: K must contain every odd k with mu(k) < t1, and E every j with mu(2j) < 0. Adding more nodes only tightens the bound (K2, neighbour windows), so the knobs affect sharpness, never validity.
- Sectors and complex f: the form is real symmetric and parity-invariant, so Q(f) = Q(Re f) + Q(Im f) and Q(u) = Q(u_even) + Q(u_odd).
Pass 2, code, read adversarially:
- Characters: tab[n mod q] gives chi_-4 = (0, 1, 0, -1), chi_-3 = (0, 1, -1), chi_5 = (0, 1, -1, -1, 1) (the Legendre symbol mod 5), chi_8 = (0, 1, 0, -1, 0, -1, 0, 1) (Kronecker (2/n)). Parities are a = 1, 1, 0, 0 (chi(-1) = chi(q - 1)).
- Atoms are n < c strictly, Lambda(p^e) = log p, and the weight is chi(n) log p/sqrt n. The node n = c sits at g(+-L) = 0.
- The float scan uses the same mu (sign of the tail term +2(-1)^k, constant log(q/pi), w = 1/4 + a/2 + ir/2). Phases are r_k log n = 2 pi k (log n/(2L)) with 64-bit fixed point. Unflagged nodes are at least 1e-6 above threshold in float, against an estimated error below 1e-8. The estimate assumes libm's cos is accurate to a few ulp, which is not proved. That is the gap the review found, now closed by rscan.c for the rows marked rigorous (RIGOR.md).
- The classification uses exact lower endpoints: G = {dl_j > 0}, K+ = {w_k > 0}. Columns sqrt(deff - dl_j) and sqrt(-w_k) are real. deff = t1/2 - eps > 0 is asserted.
- Lemma 5.5: if the eigenvalue estimate is <= 0 the run reports false (chi_-3, c = 40, K2 = 201). Whatever the estimate, success proves Id - W >= nu > 0, because Id - W = L L^T + s0 Id + E. The first code estimated lambda from one inverse-iteration vector, which finds the eigenvalue closest to 0, not the smallest. That is harmless for proofs (a negative eigenvalue makes the residual check fail), but it misled the diagnosis at chi_-3, c = 47, odd f. The final code uses Rayleigh-Ritz on 6 vectors, which exposed the two negative eigenvalues there. The residual bound uses upper endpoints of every ball entry of IW - L L^T - s0 Id, and the row-sum bound is a valid 2-norm bound for a symmetric error.
- With no negative columns (nB = 0: the prime-free chi_5 and chi_8 runs) the code returns true because P > 0. This is correct: P >= min(deff, min_G dl_j) Id > 0.
No error found in either pass.

## 11. Trust assumptions, resources, limits
- Trusted for every row: Arb / python-flint 0.9.0 encloses correctly; IEEE-754 binary64 +, -, * in the C scan; this Python and C code, which is not checked by a proof assistant.
- Trusted only for rows marked "float budget": the float scan's error budget, which assumes libm cos to a few ulp.
- The analytic chain is Chekmachev's (audited in ../ChekAudit/REPORT.md). Its port to chi is derived in DERIVATION.md. The explicit formula and Weil's criterion for one L-function are standard.
- Resources. Every run was one single-threaded process. The chi_-3, c = 47 even certificate ran on the first version of the code, which built the full Gram matrix. Its peak memory was not measured, and it probably exceeded the brief's 1.5 GB. Two odd runs with that code were stopped: one (K2 = 1001, 448 bits) passed 2.4 GB before I killed it, and a smaller one (K2 = 201, 384 bits) reached 1.52 GB. The code was then rewritten (Gram blocks built directly, Woodbury and residual in blocks). A memory watchdog (memwatch.sh) stopped four runs above 1.45 GB:
  - chi_-3, c = 47, odd, K2 = 201 with 8 neighbours (first code version);
  - chi_-3, c = 40, even and odd, K2 = 601, block solves;
  - chi_-3, c = 40, even, single solve.
  The cause is the rigorous arb_mat solve with the positive-node Gram block (order about 1000, entries of very different sizes) at 384 bits. Diagonal scaling or a blocked Cholesky would remove it; that is the next step for chi_-3 beyond c = 35 and for c > 47. Measured peaks of the certificates: 0.17 to 0.49 GB (chi_-3 at c = 30 and 35, chi_5, chi_8); chi_-4 at c = 47 was not measured.
- Not progress toward GRH itself. Windows cannot cover R: the slack falls like exp(-4 pi c/q), and the scan grows like exp(2 sigma(c)), about exp(4 sqrt c).

## 12. Files
DERIVATION.md (the form, lattice identity on L^2, minorant, conventions); RIGOR.md (proof of the scan bound); rscan.c (rigorous scan); logs/scan_*_rig.json (rigorous node sets with their bound E); logs/rig_*.log (re-runs); cert_chi.py (certificate); cscan.c (C pre-scan); scan_check.py (C vs numpy scan); galerkin_chi.py (independent upper bounds); ef_check.py (explicit formula against zeros); mu_dist.py (distribution of mu at selected nodes); table.py (all runs as a table); queue*.sh (the exact run sequences, including the stopped ones); memwatch.sh (the memory watchdog); logs/ (every run).
