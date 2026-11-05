NO RECORD: c = 47 was not re-certified and nothing beyond it was tried: its matrix step (order 7170 plus a positive block of order 6063, at about 1000 bits) needs about 35 GB of memory and about a day on one core, against this run's 1.5 GB and 3 hours. New instead: a full rigorous scan of all 3.49e10 nodes at c = 47 (and 3.45e10 at c = 45) returns exactly Chekmachev's node sets, the one premise of his record that his own checker only samples; our own zeta code certifies windows up to c = 25 (half-width (ln 25)/2 = 1.609, all f).

# WeilZeta: independent zeta certificates and the full c = 45, 47 node scans (weilzeta, 2026-11-05)

RH relevance: N1 (Weil positivity: RH holds iff Q(f) >= 0 on every window).

## 1. Results in short
- **No new record.** The N1 record stays Chekmachev's (Zenodo 23848597, 1 Nov 2026): Q(f) >= eps_c ||f||^2 for all f supported in [-(ln 47)/2, (ln 47)/2] (half-width 1.925) and for even f up to c = 71 (2.131).
- **The scan premise of the c = 45 and c = 47 certificates is now checked in full, independently.** His Theorem 5.10 needs every node outside the certificate's finite sets to have mu(k) >= 1 (odd k) and mu(2j) >= 0 (even). Our C scan has a proved rounding bound (RIGOR.md: E = 1.8e-7 against a flag margin of 1e-5); every flagged node was re-checked in Arb. Result:

  | c | nodes scanned (odd k <= kscan, even 2j <= kscan_e) | selected beyond K2 (odd / even) | Chekmachev's retained nodes with mu < 1 / mu < 0 | missing / extra | largest selected odd k, even j | time, peak memory |
  |---|---|---|---|---|---|---|
  | 47 (K2 = 3683) | 34,865,511,020 (kscan 5.10e10, 1.88e10) | 6047 / 2172 | 6047 / 2172 | 0 / 0 | 10,263,762,037; 897,154,762 (equal) | 658 s, 73 MB |
  | 45 (K2 = 3375) | 34,471,726,347 (kscan 5.04e10, 1.85e10) | 2381 / 4562 | 2381 / 4562 | 0 / 0 | 284,419,129; 2,088,261,922 (equal) | 708 s, 77 MB |

  His retained sets also contain every odd k <= K2 and every j <= K2, and our symbol equals his published 40-digit values at 2000 (c = 47) and 600 (c = 45) random retained nodes up to k = 3e9, to 1.0e-38, i.e. to his truncation (logs/chek_spot.log). So condition (i) of his Theorem 5.10 holds for his published c = 45 and c = 47 data, checked by code that shares nothing with his. His own checker (Zenodo 23986313, 5 Nov 2026) does the full scan only for c <= 31: "For c = 45, 47 and 71 a full scan would take days, so the symbol and the matrix are verified completely, while the nodes are checked on a sample or interval by interval." With this scan, every premise of the c = 47 certificate has a complete check. The final matrix inequality at c = 47 still rests on his own two implementations (his Rust code and his python-flint checker); we did not re-run it (sec. 5).
- **Our own certificates (Arb plus the proved scan, sec. 2).** cert_zeta.py ports WeilChi's pipeline to zeta. It reproduces the independent ChekAudit certificates at c = 5, 10, 16 to all 12 printed digits, matches Chekmachev's set sizes (K+, K-, bad even nodes, G, matrix order) exactly at c = 20, 25, and at his own eps gives his margins (his manifests store a rigorous lower bound about 1e-6 below the eigenvalue). Largest window certified by our code, both parities: c = 25 (half-width (ln 25)/2 = 1.609, all f).

## 2. Certificates
Theorem (computer-assisted, Arb ball arithmetic). For each row, Q(f) >= eps ||f||^2 for every f in L^2(-a, a) of the stated parity, a = (ln c)/2. Both parities at the same c give all complex f with the smaller of the two eps (the form is real symmetric and commutes with x -> -x): Q(f) >= eps ||f||^2 for all f with eps = 7.447e-18 at c = 5, 1.5128e-75 at c = 16, 2^-327 = 3.7e-99 at c = 20 and 2^-416 = 5.9e-126 at c = 25. A bound at c holds on every smaller window.

| c | half-width a | parity of f | eps | margin (ours) | Chekmachev, same eps | bits | order / U | time (s) | peak RSS (GB) | status |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.8047 | even | 7.447489666e-18 | 4.00280318265e-19 | 4.002799e-19 (eps 7.4474896668e-18) | 256 | 9 / 94 | 2 | 0.05 | rigorous (Arb + proved scan) |
| 5 | 0.8047 | odd | 5.693186289e-15 | 8.28212749750e-16 | 8.282119e-16 (eps 5.6931862892e-15) | 256 | 9 / 93 | 1 | 0.05 | rigorous (Arb + proved scan) |
| 10 | 1.1513 | even | 9.133835606e-44 | 5.02635510181e-45 | 5.026350e-45 (eps 9.1338356067e-44) | 384 | 33 / 132 | 2 | 0.06 | rigorous (Arb + proved scan) |
| 16 | 1.3863 | even | 1.512800554e-75 | 7.29362304210e-77 | 7.293616e-77 (eps 1.5128005543e-75) | 512 | 133 / 280 | 22 | 0.12 | rigorous (Arb + proved scan) |
| 16 | 1.3863 | odd | 2.105612047e-71 | 2.99318593786e-72 | 2.993183e-72 (eps 2.1056120472e-71) | 512 | 133 / 279 | 17 | 0.14 | rigorous (Arb + proved scan) |
| 20 | 1.4979 | even | 2^-327 | 3.10571228010e-97 |  | 448 | 311 / 789 | 83 | 0.70 | rigorous (Arb + proved scan) |
| 20 | 1.4979 | odd | 2^-327 | 6.70422678283e-92 |  | 448 | 311 / 788 | 81 | 0.68 | rigorous (Arb + proved scan) |
| 25 | 1.6094 | even | 2^-416 | 5.21043141861e-124 |  | 576 | 558 / 1058 | 295 | 1.08 | rigorous (Arb + proved scan) |
| 25 | 1.6094 | odd | 4.0546561810437569560e-119 | 8.07878305192e-120 | 8.078775e-120 (eps 4.0546561810e-119) | 576 | 558 / 1057 | 414 | 1.08 | rigorous (Arb + proved scan) |

"margin" = smallest eigenvalue estimate of Id - W_up at that eps (the dual check of Chekmachev's Thm 5.10); Lemma 5.5 then proves Id - W_up >= nu > 0 (Cholesky of the midpoint minus s0 = margin/2, Arb residual below s0/2: "rigorous": true in every row; the logs give nu and the residual bound). "order / U" are the orders of Id - W_up and of the positive block. Every row used the rigorous scan (rscan.c, RIGOR.md). Time and peak RSS (maximum resident set size) from /usr/bin/time -l, one single-threaded process; under memory pressure macOS compresses pages, so RSS can understate the footprint (seen in the c = 29 runs, logs/rss.log). Full run data: python3 table.py (logs/*.log).
c = 29 (order 1084, U 1205) was attempted at Chekmachev's eps (3.3887641037680338350e-145 even, 1.5329466444386175069e-140 odd). The even run exceeded 1.5 GB before finishing, at 640 and at 576 bits, and was stopped; the odd run was then cancelled (logs/z29.time, logs/memwatch.log). Not certified. Rows c = 5, 10, 16 used ChekAudit's eps (Chekmachev's, truncated to 10 digits); c = 25 odd used his exact eps.

## 3. The form and the explicit-formula check (DERIVATION.md)
Q(f) = 2 |<f, cosh(x/2)>|^2 - 2 |<f, sinh(x/2)>|^2 + (1/2 pi) int [Re psi(1/4 + it/2) - log pi] |F(t)|^2 dt - 2 sum_{n < e^{2a}} Lambda(n) n^{-1/2} Re g(log n),
with F(t) = int f e^{itx} and g = f * conj f(-.). The first two terms are the poles at s = 0, 1. By the Guinand-Weil formula Q(f) = sum_rho F(gamma) conj F(conj gamma); RH iff Q >= 0 on C_c^infty.
Check against zeros (ef_check_zeta.py, logs/ef_check_zeta.log): f = (1 - (x/1.2)^2)^8 and (x/1.2)(1 - (x/1.2)^2)^8 on [-1.2, 1.2] (prime powers below e^2.4 = 11.02), 25 zeros from mpmath.zetazero up to t = 90.

| f | sum over zeros | Q by the formula: pole + arch - primes | difference | omitted zeros above t = 90 (estimate) | lattice identity minus formula |
|---|---|---|---|---|---|
| even | 7.925427240828247e-9 | 7.925427240828625e-9 = 1.0533736 - 0.93160091 - 0.12177269 | 3.8e-22 | 3.6e-22 | 3.0e-31 |
| odd | 4.936784172559485e-8 | 4.936784172559506e-8 = -0.0010485641 - 0.0080237652 + 0.0090723786 | 2.1e-22 | 3.5e-22 | 1.2e-32 |

The terms are of size 1 and cancel to 1e-8, so a wrong sign or factor in any of them would show at the 1e-2 level. The certificate's lattice identity matches the formula to 1e-31.

## 4. Method (what is new relative to WeilChi and ChekAudit)
- cert_zeta.py: Chekmachev's periodized lattice minorant (Thm 5.10 with the rule of Def. 6.1) in python-flint 0.9.0 balls. From WeilChi's cert_chi.py: the rigorous scan and the blockwise layout. From ChekAudit's cert2.py: the zeta symbol and the rank-one pole column (+2 cosh(x/2) in the positive part for even f, -2 sinh(x/2) in the negative part for odd f). Closed forms: DERIVATION.md sec. 4.
- Woodbury step by Chekmachev's Lemma 5.3: W <= W_up = X - Y^T Z - Z^T Y + Z^T S Z for ANY Z (S > 0), with Z an approximate floating solve; the loss is quadratic in the error of Z. This replaces WeilChi's rigorous arb_mat solve, which was its memory wall. Column blocks keep the memory near X + Y + S + Z.
- Block Cholesky in floating point (panel solves and trailing updates as arb_mat products); the residual is still enclosed in Arb, so this changes speed, not rigor.
- rscan2.c: the scan with a proved early exit (RIGOR.md, last section). Atoms by decreasing weight; a node is rejected once the partial sum plus the remaining weights proves mu >= threshold + margin. On average 1.6 of the 21 atoms per node at c = 45 and 47, about 19 ns per node. On c = 20 and 25 it selects exactly the nodes of rscan.c (logs/scancmp_sets.log). At c = 45, 47 the flag margin is 1e-5, because the phase rounding term 2 pi k e 2^-64 reaches 1.8e-7 at k = 5e10.
- Validation chain: c = 5, 10, 16 equal ChekAudit's independent certificates to all printed digits (logs/val.log); every code revision reproduces them (logs/val2.log, logs/val3.log); set sizes at c = 20, 25 equal Chekmachev's manifests.

## 5. Why c = 47 was not certified here (cost model)
Chekmachev's c = 47 data: Id - W_up of order 7170 (6301 odd negative nodes, 869 bad even nodes), positive block 6063 plus the pole column, |G| = 371,000, 1200 bits; margin 6.5e-244 at eps_47 = 1.87e-242. Any certificate must resolve about 1e-244 against entries of size 1, so it needs at least about 850 bits.
- Memory: X, Y, S, Z hold 7170^2 + 2 * 6064 * 7170 + 6064^2 = 1.75e8 balls, at 900-1200 bits about 180-220 bytes each: 32-38 GB.
- Time on one core, extrapolated from our runs: Cauchy sums over G, 13234 columns x 371,000 terms = 4.9e9 terms (3.5 h in Python at 2.6 us per term); Gram assembly 8.8e7 entries (about 45 min); approximate solve and products about 9e11 multiply-adds at 1000 bits plus a Cholesky of 1.2e11 (about 10 h). About a day in all.
- Scaling seen here: c = 20: order 311 / U 789, 2 min, 0.5-0.7 GB; c = 25: 558 / 1058, 5-7 min, 1.08 GB; c = 29: 1084 / 1205, did not fit: the even run passed the 1.45 GB watchdog after 17 min at 640 bits and after 13.5 min at 576 bits (max RSS 1.54 and 1.50 GB) and was stopped. c = 30 (U = 2739) and c = 31 (order 1663, U = 1628) do not fit in 1.5 GB in this layout.
- Beyond the record: all f at c = 53 or more needs a matrix larger than at c = 47 (the scan itself, about 1e11 nodes, would take about 30 min with rscan2); even f at c > 71 needs a scan beyond 1.7e13 nodes (3-4 days on one core with rscan2) and a matrix of order beyond 38,377. Neither fits; nothing beyond c = 47 was attempted.

## 6. Replay (one row in under 10 minutes)
From Riemann/Deep/WeilZeta; Python packages only via uv; one single-threaded process:
    cc -O2 -ffp-contract=off -o /tmp/rscan rscan.c
    export WZ_RSCAN=/tmp/rscan WZ_SCAN=/tmp WZ_CHOL_BS=64 OMP_NUM_THREADS=1
    uv run --no-project --with mpmath --with numpy --with python-flint python3 cert_zeta.py 20 even 2^-327 448        # 2 min, 0.5 GB: "rigorous": true, margin 3.10571228010e-97
    uv run --no-project --with mpmath --with numpy --with python-flint python3 cert_zeta.py 16 odd 2.105612047e-71 512   # 10 s: margin 2.99318593786e-72
    uv run --no-project --with mpmath python3 ef_check_zeta.py 1.2 8 90                                                  # 1 min: zeros against the formula
The c = 47 scan (11 min, 73 MB) and the comparison (needs truncated-weil-form-certificates-en-v1.1.zip, 212 MB, from Zenodo 23847003):
    cc -O2 -ffp-contract=off -o /tmp/rscan2 rscan2.c
    WZ_MARGIN=1e-5 WZ_RSCAN2=/tmp/rscan2 WZ_SCAN=/tmp uv run --no-project --with mpmath --with numpy --with python-flint python3 scan_zeta.py 47 3683 24
    python3 cmp_chek.py /tmp/scan_zeta_47_3683_rig.json <dir>/certificates/c47_even_eps/nodes_odd.csv <dir>/certificates/c47_even_eps/nodes_even.csv 3683
Arguments of cert_zeta.py: c, parity of f, eps (0, 2^-k or a decimal), bits, optional K2 (default max(201, ceil(5c^2/3)), Chekmachev's). One JSON line per run; "rigorous": true means Lemma 5.5 succeeded at that eps. The selected node lists of both full scans are in logs/scan_zeta_47_3683_rig.json and logs/scan_zeta_45_3375_rig.json.

## 7. Novelty pass (2026-11-05)
- python3 .scratch/lst.py 'Weil|positivity|explicit formula' math.NT --zenodo 'Weil positivity': math.NT listings for 2026-10 and 2026-11 (12 title hits) and Zenodo. Relevant:
  - Zenodo 23986313 (V. Chekmachev, 5 Nov 2026), "Independent verification of certificates of lower bounds for the truncated Weil form": his python-flint checker of his own published data (symbol balls, scan, matrix). Complete for c <= 31 (c = 31 even: 33 min on 18 processes); at c = 45, 47, 71 the node scan is sampled or done interval by interval. Our full scans fill that gap for c = 45 and 47.
  - Chekmachev's paper 23848597 is still v1.0, the only version of concept 23848596. His data record has v1.1 (23847003). His companion paper 23978194 moved to v2.0 on 5 Nov 2026; it is theory (boundary balance r_c), claims no new window, and says boundedness of r_c is open.
  - arXiv 2610.32704, "Spectral convergence associated with Weil's quadratic form": assumes RH; eigenvalues of the localized form converge to zero ordinates. No certified window.
  - arXiv 2611.02944, "On some positivity properties of the Riemann Xi function": title-level hit, not a Weil-window certificate.
- Our files: ../ChekAudit (audit; c <= 16 re-certified; float scan of k <= 2e8 at c = 47 only, 95% / 98% of the selected nodes), ../WeilChi (Dirichlet L-functions), ../WeilRecord (our earlier 1.55 and 1.6).
- Verdict: no record. The full c = 45, 47 scans are a new independent verification of a published certificate, not new mathematics. Our certificates up to c = 25 (half-width (ln 25)/2 = 1.609, all f) replicate known results with independent code.

## 8. Self-review
- Mathematics. The explicit formula is checked against zeros for both parities; the pole signs (+cosh for even f, -sinh for odd f), the constant -log pi, b_n = 2n + 1/2 and the weights Lambda(n)/sqrt n all enter that check. The lattice identity on L^2 and the minorant are Chekmachev's, audited SOUND in ../ChekAudit and re-derived for chi in ../WeilChi.
- Lemma 5.3, new in our code path: (Z - S^{-1}Y)^T S (Z - S^{-1}Y) >= 0 with S = D^{-1} + U^T P0^{-1} U > 0 (D > 0: the positive K+ weights and the pole weight 2). So W <= W_up for every Z, and Id - W_up > 0 implies Id - W > 0. Switching to it left the c = 16 and c = 20 margins unchanged to 12 digits.
- Scan. The early exit rejects only nodes whose proved lower bound already clears threshold + margin; the others get the full rscan.c test. Flagged at c = 47: 6049 odd, 2173 even; the Arb re-check (96 bits) kept 6047 and 2172. At c = 45: 2381 and 4566 flagged, 2381 and 4562 kept.
- Our kscan (c = 47: 50,977,469,937 odd, 18,753,563,152 even) is slightly below his (50,977,471,206; 18,753,563,601): our monotone bound uses the tail sum itself (WeilChi DERIVATION sec. 5), his uses Lemma 4.2. Both are valid; beyond our kscan no node can be selected, so the comparison covers all k.
- The comparison uses his even-sector node files; per his README the node sets of the two sectors coincide at c = 45, 47 (his odd-sector manifests have the same k_scan values).
- A bug found and fixed: memwatch.sh (copied from ../WeilChi) used the caller's variable p inside its ancestor test, so it killed the parent queue script instead of the oversized process. Here that stopped one queue; the oversized c = 29 run (640 bits, 1.54 GB) was then stopped by hand, and for about 45 s it overlapped the c = 45 scan (two processes at once). Fixed in this folder (memwatch.sh header). ../WeilChi/memwatch.sh has the same bug.
No error found in the mathematics or in the certificate code.

## 9. Trust assumptions
Arb / python-flint 0.9.0 encloses correctly; IEEE-754 binary64 +, -, * in rscan.c and rscan2.c (compiled with -ffp-contract=off); this Python and C code, which no proof assistant has checked. The scan comparison reads Chekmachev's published node files as data; our scans do not depend on them.

## 10. Files
DERIVATION.md (zeta form, explicit formula, lattice identity, minorant, Lemma 5.3); RIGOR.md (scan error bound, early exit); cert_zeta.py (certificate); rscan.c, rscan2.c (scans); scan_zeta.py (chunked full scan); cmp_chek.py and chek_spot.py (comparison with Chekmachev's data); ef_check_zeta.py (explicit formula against zeros); table.py, report_table.py (tables); queue*.sh, chain*.sh, memwatch.sh, rsslog.sh (the exact run sequences); logs/ (every run, including the stopped one).
