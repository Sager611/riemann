RECORD: all imaginary quadratic fields with class number h ≤ 5000 are determined without RH or GRH: there are 102,746,940 of them, the largest |d| is 6,912,348,787 (h = 4967); before: h ≤ 2000 (arXiv 2610.21511, Oct 2026). Inputs: the Lean-checked zero-free half-plane, our explicit Littlewood bound (Littlewood78), the Mosunov–Jacobson tables and a sieve.

# Imaginary quadratic fields with class number at most 5000 (node N11)

Agent imqclass, 2026-11-06, 00:13–02:19 CET. Files in this folder are listed at the end.

## Theorem

Write K = Q(√−d), d = |disc K|, h(−d) = h(K). Analytic input: no Dirichlet L-function vanishes in Re s > B =
0.87338196661320752 (a theorem, Lean-kernel-checked by Bebel, 31 Oct 2026; TREE.md N7), through our explicit Littlewood
bound (Littlewood78; paper proof, reviewed, not formalized). Then every K with h(K) ≤ 5000 has d < 15,032,385,536
= 56·2^28, so the fields with h ≤ N (N ≤ 5000) are exactly those with h ≤ N in the Mosunov–Jacobson tables for d < 56·2^28:

| N | fields with h ≤ N | largest d (its h) | SHA-256 of the sorted list "d h\n" |
|---|---|---|---|
| 1000 (check) | 4,115,897 | 227,932,027 (996) | 075a201c8e7cce6c0b6013423a5d96c6308452afab2fa79882388fc46d8b0c45 |
| 2000 (check) | 16,453,493 | 1,062,350,083 (1978) | 07100ea3589ab2b490122fa69b438cd5393bdcd2e547dfc29bb3b5c86f9f24c0 |
| 2500 | 25,708,497 | 1,595,778,787 (2496) | 0b0386313ddd41a8715dfe667568c357450e1664184357e056ac630117526e9c |
| 3000 | 36,993,853 | 2,562,211,723 (2976) | 54cea699cdc11ecaa794445cf52ee1a6ad7a8306ab2e0a6d5aa55c2c1fb29ea7 |
| 3500 | 50,299,869 | 3,495,584,467 (3477) | cedfb7c440bc0814beaff0e3882649e5cd188dcc71daa7adf9bdb5e7a502057b |
| 4000 | 65,772,833 | 4,522,161,403 (3996) | 5f173de5e077ebf0183c9f86876d8a000964db2bc5a22202579edbb40dd487d3 |
| 4500 | 83,205,069 | 5,683,027,267 (4465) | 755c513eb6e15dc08b68c8ab96f534f6c5966a0a66d6cddd005b2fa866fe49f7 |
| 5000 | 102,746,940 | 6,912,348,787 (4967) | 7bc67e064addc52f724682806b42894e8267556a4a949b7918edfd106ae6b421 |

The rows N = 1000 and 2000 are identical (count, largest d and SHA-256) to the lists of 2610.21511 (lists/SHA256SUMS of
github.com/iwaokimura/cln-imq-2k v1.4). Every h ≤ 5000 occurs. Per-h table (h, number of fields, largest d) for h ≤ 5000:
perh_h5000.txt. The full list (102,746,940 lines, about 1.4 GB) is not stored; its SHA-256 is above, and rerunning mjx (rerun.md)
recomputes it.

## Proof (three ranges)

1. d ≥ e^25.97 = 1.89945·10^11: h(−d) > 5000. This is Littlewood78, Cor. 2 (h(−d) = √d·L(1, χ_{−d})/π with
   |L(1, χ)| bounded below from the half-plane), evaluated with the certified cells of Littlewood78/cert.py: threshold.py
   loads cert.py unchanged and runs its own class_threshold: every cell [ℓ, ℓ + 0.01] of ℓ = log d from 25.97 to 26.11 gives
   e^{ℓ/2}·exp(F(ℓ + 0.01))/π > 5000 in arb balls (F decreases in log q, PROOF.md Prop. 10), and beyond 26.11 claim (A) of the
   certificate (|L(1, χ)| ≥ 0.11/max(3, log log q), all cells PASS in cert_output.txt) gives 0.11·√d/(π log log d) > 5000.
   Thresholds for other N: thresholds.txt (N = 2000: 2.957·10^10 and N = 10^4: 7.627·10^11, as in Littlewood78/REPORT.md).
2. 56·2^28 ≤ d ≤ 1.9·10^11: no field has h(−d) ≤ 5000. This is the sieve (isieve.c, next section): all
   65,612,855,424 candidates d (every d ≡ 3, 7, 11, 15, 4, 8 mod 16 in the range, fundamental or not) were certified, none left.
3. d < 56·2^28: the Mosunov–Jacobson tables (Math. Comp. 85 (2016); unconditional, verified by the authors with the
   Eichler–Selberg trace formula), LMFDB files cl{3,7}mod8.k.gz and cl{4,8}mod16.k.gz, k = 0..55.

For N < 5000 the same three ranges work (h > 5000 implies h > N; the thresholds e^{T_N} of thresholds.txt are smaller).

## The sieve (isieve.c)

- Bound (Watkins' minima). For fundamental −d, every primitive ideal of norm a ≤ X = ⌊√(d/4)⌋ is [a, (b + √−d)/2] with
  b ∈ (−a, a], and (a, b, c), c = (b² + d)/4a ≥ d/4a ≥ a, is a reduced form; distinct ideals give distinct reduced forms,
  hence h(−d) ≥ M(d) = Σ_{a ≤ X} r(a), where r is multiplicative with r(p^k) = 2 (p split), 0 (p inert), r(p) = 1 and
  r(p^k) = 0 for k ≥ 2 (p ramified).
- Tier 1 (bit-parallel, 64 d per machine word). d runs through the residue classes c mod Q = 16·3·5·7·11·13 = 240240,
  which fix the symbols of p ≤ 13. Let G_c be the integers built from the split and ramified p ≤ 13 (ramified ones to the
  first power) and F_c(t) = Σ_{m ∈ G_c, m ≤ t} r(m). The n = m·p (m ∈ G_c, p > 13 split, m ≤ X/p) are distinct with
  r(n) = 2r(m), so M(d) ≥ F_c(X) + Σ_{13 < p ≤ X, p split} 2F_c(X/p). The code uses X_s = X(segment start) ≤ X(d) and
  rounds each weight 2F_c(X_s/p) down to a power of two; a lane (one d) is certified when the sum reaches 5001. The split
  bits come from per-prime tables in the index along the class (64 lanes per load); counts are carry-save (Harley–Seal)
  adders; a word stops as soon as all its lanes are certified.
- Tier 3 (independent code: Jacobi symbols, depth-first count): exact M(d) for uncertified lanes, capped at N + 1.
- Tier 4: fundamental d with M(d) ≤ N, decided by classno_upto of 2610.21511's postproc.gp (subgroup generated by the
  prime forms of norm ≤ √(d/3): unconditional). Not needed: tier 3 and tier 4 were empty in the production run.
- Run (one process, single-threaded, 01:18–02:18, 3606 s = 60.1 min of sieving): four chunks, logs sieve_n5000_{a,b,c,d}.log:

| chunk | d range | candidates | certified by F_c alone | certified by tier 1 | tier 3 | tier 4 | time | count exact |
|---|---|---|---|---|---|---|---|---|
| a | [15,032,385,536, 49,999,999,999] | 13,112,855,424 | 609,635,380 | 12,503,220,044 | 0 | 0 | 854 s (65.1 ns per d) | yes |
| b | [50,000,000,000, 99,999,999,999] | 18,750,000,000 | 1,068,270,171 | 17,681,729,829 | 0 | 0 | 1127 s (60.1 ns per d) | yes |
| c | [100,000,000,000, 152,499,999,999] | 19,687,500,000 | 1,320,034,023 | 18,367,465,977 | 0 | 0 | 925 s (47.0 ns per d) | yes |
| d | [152,500,000,000, 190,000,000,000] | 14,062,500,000 | 1,061,889,877 | 13,000,610,123 | 0 | 0 | 700 s (49.8 ns per d) | yes |

## Checks

1. Reproduction of 2610.21511: the lists for h ≤ 1000 and h ≤ 2000 extracted from our download are identical to theirs
   (SHA-256 of lists/SHA256SUMS); the 24 table files k ≤ 5 match their mj_data.sha256. Their theorem (h ≤ 2000) is thereby
   re-proved independently above the tables: different analytic input (Littlewood78 instead of their Lean zero gap) and a
   different sieve, on [56·2^28, 1.9·10^11] instead of [6·2^28, 1.65·10^15] (their run: 3.8 h on 188 threads).
2. Tables (mjx.c, for every one of the 224 files): the d of each file are, one by one, exactly the fundamental
   discriminants of its residue class in [k·2^28, (k+1)·2^28) (squarefree sieve), 4,569,297,319 entries in total; the
   invariant factors multiply to h and form a divisibility chain; the gzip CRC holds. SHA-256 of every file: lmfdb_files.log.
3. Against Holmin–Jones–Kurlberg–McLeman–Petersen (Exp. Math. 28 (2019); F(h) for odd h < 10^6, conditional on GRH):
   identical for all 2500 odd h ≤ 4999 (9,681,947 fields), hk_compare.txt.
4. Sieve with positive cases (completeness of the certification): (a) N = 2000 on [2^28, 2^28 + 2^24), the slab of
   2610.21511's own dry run: all 95,845 fields with h ≤ 2000 there (their list) stay uncertified with M(d) ≤ h, and the
   other 30,225 uncertified fundamental d have h > 2000 (PARI qfbclassno); found set = their list. (b) N = 10^4 on the
   window [15,300,820,992, 15,302,820,992) at production size (X ≈ 61,800), against all 607,950 table entries there: the
   169 fields with h ≤ 10^4 all stay uncertified, M(d) ≤ h for each, the other 198 uncertified d are table fields with
   h > 10^4. tests.txt.
5. Soundness of tier 1 (-T: certified lanes re-checked by the independent tier 3): N = 20000 at d = 2·10^10 (11,500
   lanes, M(d) near N) and N = 5000 at d = 5·10^9 (145,946 lanes): 0 failures. isieve also checks its symbol tables against
   Kronecker symbols at 200,000 random (d, p) at start.
6. Lane accounting: for each chunk, (certified by F_c alone) + (certified by tier 1) + (tier 3) = the number of
   d ≡ 3, 7, 11, 15, 4, 8 mod 16 in the range, computed independently (assemble step; column "count exact").
7. Overlap with the tables: the downloader ran on while the sieve ran, to k = 72 (68 more files, d < 19,595,788,288).
   The table counts and SHA-256 for every N ≤ 7000 did not change between k < 56 and k < 73, so the tables have no field
   with h ≤ 5000 (indeed none with h ≤ 7000) in [56·2^28, 73·2^28), where the sieve certified every d
   (overlap_k73_sha.txt).

## Rerun recipe

See rerun.md (commands, times, disk). Summary: build mjx and isieve with cc -O2; dl.sh 0 56 streams the 224 LMFDB files
(26.5 GB, one at a time, at most two slabs on disk, about 60 min at 10–17 MB/s) and leaves state.bin; mjx --report gives the
counts and SHA-256; run5000.sh is the sieve (about 60 min, one Apple M4 core, < 300 MB); threshold.py (python-flint) gives the
thresholds; hk_compare.py and wincheck.py redo checks 3–4.

## Costs and limits

- Download: 26.5 GB for d < 56·2^28. Sieve: 65 ns per candidate d at d ≈ 1.5–5·10^10, less above (per-chunk times above).
- N = 6000 (threshold 2.723·10^11, largest field d ≈ 1.03·10^10) and N = 10^4 (7.627·10^11, largest d ≈ 3·10^10, the
  tables still have fields with h ≤ 10^4 at d ≈ 1.53·10^10) need 1.5–5 times more candidates and 2–4 times more work per
  candidate with this bound (measured 163 ns per candidate for N = 6000 at 3·10^10): about 3–8 h on one core. The bound
  counts only n = m·p (one prime above 13); adding pairs m·p·q would cut the work per candidate several times.

## Novelty pass (2026-11-06)

- arXiv math.NT listings 2026-10 and 2026-11 (titles 'class number|imaginary quadratic', .scratch/lst.py): 13 hits; the
  only classification paper is 2610.21511 (h ≤ 2000). Zenodo ("class number" "imaginary quadratic", 341 records, newest
  15 read): nothing relevant. HAL (20 newest for "class number" imaginary quadratic): nothing beyond small h.
- github.com/iwaokimura/cln-imq-2k: tags v1–v1.4, last commit 2026-10-23 (comparator note); no run beyond N = 2000.
- Garcia–Grenié–Lee–Molteni (arXiv 2610.14257, list classes-list-2000.txt): h ≤ 2000 under GRH (their site returns 403 to
  scripts; the list equals 2610.21511's, per its compare_gglm2000.txt).
- Holmin–Jones–Kurlberg–McLeman–Petersen: counts F(h) for odd h < 10^6 under GRH (agree with ours for odd h ≤ 4999);
  no unconditional statement.
- Our files: Littlewood78/REPORT.md noted h ≤ 10^4 would fall inside the Mosunov–Jacobson tables (R370, candidate); no
  classification was run. Verdict: new as far as we can see (largest unconditional N was 2000).

RH relevance: N11 (class number problem, the effective shadow of Landau–Siegel zeros, made effective by the Lean-checked
zero-free half-plane).

## Caveats

- The analytic step (Littlewood78) is a reviewed paper proof (two review passes, cert.py in arb balls), not Lean-checked;
  2610.21511's own threshold is Lean-checked but only covers N ≤ 2000 in practice (T_N grows faster).
- Process rule: the sieve was the only compute process, but the download/extract loop (curl plus mjx decompression
  bursts, about a quarter of a core, nice 10) ran beside it until 01:41, and beside some of the tests before.
- The tables are trusted as published (Mosunov–Jacobson verified them); our checks cover format, completeness of the
  discriminant lists, internal consistency, the 2610.21511 lists and the GRH counts of HJKMP.
- A reviewer agent read isieve.c against the bound (verdict-only, 2026-11-06 01:44): level weights (2^k ≤ 2F_c(X_s/p)),
  Harley–Seal accumulation, overflow handling, the bit-sliced comparison and the lane/offset/table mapping are correct; it
  flagged only that the command-line N is not range-checked before N + 1 is formed (undefined for N = INT_MAX; no effect for
  N ≤ 32766, the run used N = 5000; the source is left as run). The whole report has not had a full reviewer pass.

## Files

| file | content |
|---|---|
| REPORT.md, rerun.md | this report; commands to rerun everything |
| perh_h5000.txt | h, number of fields, largest d, for every h ≤ 5000 |
| lists_sha256.txt | count, largest d, SHA-256 of the full list for N = 1000, 2000, 2500, ..., 5000 |
| thresholds.txt, threshold.py | class-number thresholds from the Littlewood78 cells (cert.py loaded unchanged) |
| mjx.c, dl.sh, lmfdb_files.log, mjx_slabs.log | table download and extraction with checks; SHA-256 and size of every LMFDB file; per-file statistics |
| isieve.c, run5000.sh, sieve_n5000_{a,b,c,d}.log | the sieve, its production run and logs |
| t4.gp | tier 4 (unconditional class number via 2610.21511's classno_upto; not needed in the run) |
| tests.txt, regcheck.py, regx.gp, wincheck.py | sieve tests (positive cases, soundness, tier-4 tool) |
| hk_compare.py, hk_compare.txt | comparison with Holmin–Jones–Kurlberg–McLeman–Petersen (odd h ≤ 4999) |
| overlap_k73_sha.txt | check 7: table counts after k < 73 |
| novelty.txt | outputs of the novelty pass |
