AUDIT PASSED: every item of class-number-10000 holds as written, given its two stated inputs (the published Mosunov–Jacobson tables and the refereed Littlewood78 bound); no false step, no counterexample, one test-design weakness and two wording slips.

# ImqAudit: line-by-line referee of class-number-10000 (Riemann/Deep/ImqClass2), 19 Nov 2026

Agent imqaudit (Orca), node N11. Read-only web (cln-imq-2k v1.4 postproc.gp, 71 lines). Nobody contacted. No Lean, no
commits; Riemann/Deep/ImqClass2/, ImqClass/ and discoveries/ untouched. Chain with source line numbers: FULLPROOF.md.
Every compute process ran alone, single-threaded, under .scratch/memcap.sh 290 (peaks below; max 218 MB, isieve2).

## Verdict per item (claim of discoveries/class-number-10000.md, verbatim pieces)

| Item | Claim | Holds? | Evidence |
| --- | --- | --- | --- |
| (1) | for fundamental −d, h(−d) ≥ M(d) = Σ_{a ≤ ⌊√(d/4)⌋} r(a) | yes | FULLPROOF Step 1 (reduced forms (a, b, (b²+d)/4a), b ∈ (−a, a], c > a; primitive); hcount.c: M = M2 ≤ h on 300 + 528 + 33 + 127 d |
| (2) | d ≡ 3, 7, 11, 15, 4, 8 mod 16 cover every fundamental d; non-fundamental handled | yes | Step 1d; isieve2.c l.322-329 (classes), l.148-151 (fundamental() before tier 4) |
| (3) | isieve2.c bound, weights, bit-sliced exactness | yes | Step 2 (all of l.89-505 read); window test: tier 4 = own exact {M ≤ N} at N = 10⁴ and 3·10⁴ |
| (4) | 275,512,500,000 candidates, 38 chunks, no gap, 33 tier 4 | yes | logcheck.py: lanes = own count per chunk, consecutive [2.8·10¹⁰, 762,700,000,000], bad 0, t4 = 33 |
| (5) | t4.gp (Kimura's classno_upto) unconditional; 33 tier-4 → only d = 28,154,696,707 (h = 9819) | yes | Step 3b (code read); own reduced-form count: 32 have h ≥ 10,050, M agrees on all 33 |
| (6) | d ≥ e^27.36: h > 10⁴ (claim (A) beyond 27.53, 17 cells) | yes | Step 4: F decreasing in log d (−kp W, kp > 0); threshold.py rerun identical; e^27.36 = 7.62600·10¹¹ < sieve top |
| (7) | d < 111·2²⁸: the tables, extracted with checks | yes (input trusted) | Step 5: mjx.c checks read; mjx_slabs.log: 444 files, k = 0..110 contiguous, 9,056,999,962 entries |
| (8) | 410,810,155 fields; largest d 28,154,696,707 (h = 9819); every h ≤ 10⁴ occurs | yes | perh_h10000.txt: all 10⁴ h present, total 410,810,155, 18 rows recomputed; 127 sampled (h, largest d) exact by hcount |

## Falsification (all passed; files in this folder, raw outputs summarized)

- hcount.c (own, independent of isieve2.c): exact h by counting reduced forms (b ≤ √(d/3), divisors of (b² + d)/4 from a
  root sieve), M = number of reduced forms with a ≤ X, and M2 = Σ r(a) by a multiplicative enumeration. Validated against
  PARI qfbclassno(−d, 1) on 300 random fundamental d ∈ [10⁵, 10⁹] (seed 20261119): 0 mismatches; M = M2 ≤ h for all;
  and (consistency only, GRH-conditional) against PARI quadclassunit on 23 sieve-range samples d > 3·10¹¹: 0 mismatches.
- 33 tier-4 d of chunk x5: hcount M equals x5_t4.txt's M for all 33; h ≤ 10⁴ only for 28,154,696,707 (9819); others ≥ 10,050.
- 528 sieve-range samples (gen_adv.py, seed 1119): 400 hardest-class d (d ≡ 3 mod 8 and (−d|p) = −1 for every odd
  p ≤ 47, log-uniform on [2.8·10¹⁰, 7.627·10¹¹]), 88 d next to the range ends, chunk edges, 111·2²⁸ and e^27.36, 40 d with
  2, 3, 5, 7, 11, 13 split. All fundamental; M = M2 ≤ h; min h 11,684 (d = 33,781,251,043, M = 11,091); M ≥ 10,001 for all
  (so no tier-1-certified lane with M ≤ 10⁴); M/h ∈ [0.949, 0.958] (advcheck.py). Top regime (phase B cut at
  PT = 250000 < X): 1500 d ∈ [5·10¹¹, 7.627·10¹¹] with d ≡ 3 mod 8 and every odd p ≤ 59 inert (gen_adv2.py, seed 7;
  1498 fundamental): exact min M = 41,354 (d = 502,660,741,507), none ≤ 10⁴.
- Window test, production binary (isieve2.c sha256 be745424…, cc -O3 -mcpu=apple-m4) on [28,154,600,000, 28,154,800,000]
  (75,000 lanes, 60,780 fundamental d, own exact M for each by hcount MONLY): N = 10⁴: tier 4 = {28,154,696,707} = own
  {M ≤ N}; N = 3·10⁴ (12,152 tier-3 lanes): tier 4 = own {M ≤ N} exactly (8045 d, equal M values), 0 certified lanes
  with M ≤ N, nearest certified lanes M = N + 2 (wincmp.py). A tier-1 over-count anywhere in the window would show here.
  Second window, production regime with pairs and phase B (nA 1292, nP 15,999): [131,000,000,000, 131,000,040,000], N = 3·10⁴,
  12,158 fundamental d: tier 4 = own {M ≤ N} (8 d, equal M), 0 certified lanes with M ≤ N, nearest M = N + 181; peak 140 MB.
- Sieve rerun: one production segment of chunk c08, [130,904,138,240, 132,872,184,319] (738,017,280 lanes, -T 20000):
  lanes = certF 18,743,296 + cert1 719,273,984, t3 0, bad 0; 46 s, peak 218 MB.
- Per-h table: all h = 1..10⁴ present (count ≥ 1), total 410,810,155; the 18 rows (count, largest d) of lists_sha256.txt
  recomputed from perh_h10000.txt: identical; 127 h (120 random, seed 1119, plus 1..5, 9819, 10⁴): the listed largest d is
  fundamental with exactly that h (hcount; perhsample.py).


## Code audit notes (isieve2.c, threshold.py, mjx.c, t4.gp)

- Bound (Step 2): a = m·n (m 13-smooth, n 13-rough) is unique; the n used (1; p^k with 17 ≤ p ≤ Xs/17; pq with
  17 ≤ p < q, pq ≤ Xs; p with Xs/17 < p ≤ min(Xs, 250000)) are distinct, so the sum is ≤ M(d); every rounding (round2,
  cap 32767) is downward; Xs = X(segment start) ≤ X(d), F_c nondecreasing. Phase B runs over pmed = primes ≥ 17
  (l.284-287), so p ≥ 17 there too.
- Exactness: T_p bit x = [−QS·x nonzero square mod p], built directly to bit p + 64VW + 64 (periodic, no wrap logic);
  lane j of class c at x = c·QS⁻¹ + j mod p (l.223, l.222; QS⁻¹ verified l.297); fetch's two-step shift is defined for
  shift 0. Harley–Seal: o + 2t + 4f + 8e + 16·c16 = item count per chunk (≤ 256; c16 has 5 bits, ≤ 16); addsh ripple-adds
  count·w' over 16 planes, any bit at plane ≥ 16 sets a sticky overflow (sum ≥ 2¹⁶ > 10,001); geq is an exact MSB-first
  comparison. Items are masked by the live lanes; tail lanes by GVM.
- Pairs: posT[t] = number of pairs with ⌊Xs/pq⌋ ≥ t (prefix, products sorted); posT[tmax + 1] = 0 since 17·19 = 323 is the
  least product; each emitted run [posT[tb+1], posT[t+1]) has one t-range of equal round2(4F(t)); every pair index < nA
  checked (l.373), so the stored phase-A bits exist (a group leaving phase A early never enters the pairs).
- Lanes: jlo/jhi (l.408-409) give exactly the d ≡ c mod QS in the segment; segments [LO + i·LSEG, …] are contiguous to
  HI; chunk logs contiguous (logcheck.py). Lesson 277: the scans start at 28,000,000,000 < 29,796,335,616 (tables) and
  end at 762,700,000,000 > e^27.36 = 762,600,393,219.6 (decimal check), so the three ranges overlap at both joins.
- threshold.py loads cert.py's definitions only (imports, constants up to the first say(), functions) and reruns
  class_threshold; cert.py sha256 009ce50f… unchanged; rerun output identical (thr10000.txt). Cells: Step 4.
- mjx.c: Step 5. The completeness check is two-sided (each d must be the next fundamental d; none may remain at the
  file end); the four classes are merged with a strict order check; gzip errors and truncation abort.
- t4.gp/classno_upto: Step 3b; unconditional (generation by prime classes of norm ≤ √(d/3); exact coset enumeration with
  normalized reduced keys); it returns min(h, N + 1). Used on fundamental d only.


## Weaknesses and slips (none changes a statement)

- W1 (test design, no proof step uses it): isieve2's -T K tests lanes with d·2654435761 ≡ 0 mod K, i.e. d ≡ 0 mod K. For
  16 | K no lane of the six classes is ever tested (ImqAudit's -T 20000 run: "tested 0"); even K tests only d ≡ 8 mod 16,
  always with K's primes ramified. In tests.txt the K = 3000 and 9000 runs (l.3-4, l.32) tested only d ≡ 8 mod 16, and
  the K = 2000 run (T7, l.61) was aborted; the odd-K runs (1499, 19, 49, 4499, 5) carry the evidence (its caveat says
  so). ImqAudit's window tests check every lane instead.
- S1 (formula wording): ImqClass2/REPORT.md l.68 and isieve2.c l.13 write phase B as Σ_{Xs/17 < p ≤ min(Xs, PT)}
  without "p ≥ 17"; as written the display is false for tiny Xs (reviewer: d = 23, Xs = 2 gives 5 > M = 3). The code
  sums over pmed = primes ≥ 17 only (l.284-287), so the implemented bound is right; production Xs ≥ 83,666.
- S2 (condensed wording): the page says the sieve "certified M(d) ≥ 10,001 for all … candidates except 33 tier-4
  entries"; precisely 177 lanes needed the exact tier 3 (144 with M ≥ 10,001, 33 tier 4), as ImqClass2/REPORT.md l.50-51 and l.154 say.

## Reruns and reproductions

- Production slice: chunk c08's first four segments [125,000,000,000, 132,872,184,319] rerun with the production build:
  the log line "seg 4 … X=180903 nA 1292 nP 15993 lanes 2952069120 certF 74973184 cert1 2877095936 (A 2439912960,
  pairs 413885440, B 23297536) t3 0" is identical to sieve_logs/p_c08.log; peak 208 MB.
- threshold.py 10000 (python-flint): "10000 27.36 7.627e11 27.53", identical to thresholds.txt.
- Holmin et al. (GRH, consistency only): class_group_orders.txt (6.9 MB, sha256 8b9551b8…, as in ImqClass2, deleted
  after use) against perh_h10000.txt: 5000 odd h ≤ 9999, 0 differences, 36,132,505 fields both (hkcheck.py, own script).
- Not rerun (forbidden or out of scope): the 53 GB table download and mjx extraction (the 18 SHA-256 rows), the full
  sieve (3.5 h). Their correctness rests on the code audit (Steps 2, 5), the chunk and slab logs, the reruns above and the
  cross-checks (rows ≤ 5000 = ImqClass, ≤ 2000 = 2610.21511, odd h = HJKMP, 127 largest d exact).


## Reviewer (verdict only, the delicate step pasted inline: Claim L = h ≥ M(d) with r at p = 2, and the sieve inequality)

GAP in ImqAudit's paste, not in the code: the pasted phase-B sum lacked "p ≥ 17" (counterexample d = 23, Xs = 2); with
p ≥ 17 (as implemented, l.284-287) the inequality is a valid lower bound; Claim L correct including p = 2; even the literal
display certifies nothing below Xs = 221 (it is ≤ 8224 < 10,001 there), so the conclusion h > 10⁴ holds. Also noted a typo
in the paste (r(1) = 1). Resulting status: CORRECT (S1 above records the wording).

## What survives if an input fails

- Littlewood78 wrong: range (i) fails, so completeness is lost; what survives is "every field with h ≤ 10⁴ and
  d ≤ 7.627·10¹¹" (tables + sieve); the literature's h ≤ 2000 (2610.21511) does not depend on our files.
- Tables wrong: nothing below 2.98·10¹⁰ is ours; the cross-checks (ImqClass rows, 2610.21511 to 2000, HJKMP for odd h,
  127 exact largest d) agree with the published data but cannot exclude an error inside it.
- The sieve and tier 4 hold independently of both (Steps 1-3 plus the falsification above).

## Dependents (rg -n -i 'class-number-10000|imqclass2|410,810,155' …)

RECORDS.md l.20 (the record); TREE.md l.15 (N11 node); discoveries/class-number-5000.md l.3 (superseded pointer);
Riemann/Tree/R399.lean (inputs check), R435.lean (Bhattacharjee row comparison with perh_h10000.txt); Riemann/Deep/
HPAudit/REPORT.md l.136 and LittlewoodAudit/REPORT.md (listed as dependent of the half-plane). All unaffected (PASSED);
suggested page fixes: S1 (add p ≥ 17 to the phase-B display), S2 optional.


## Files (Riemann/Deep/ImqAudit/)

| file | content |
| --- | --- |
| REPORT.md, FULLPROOF.md | this report; the chain with source line numbers |
| hcount.c | own exact h (reduced forms), M and M2 = Σ r(a); MONLY=1 for M2 only (cc -O2 -o hcount hcount.c -lm) |
| gen_small.py, pari_h.gp, pari_qcu.gp | 300 random fundamental d ≤ 10⁹ and PARI qfbclassno(−d, 1) (validation of hcount); quadclassunit at d > 3·10¹¹ |
| gen_adv.py, gen_adv2.py, advcheck.py | 528 sieve-range samples (hardest classes, range ends, split-heavy), 1500 top-range hardest classes; checks |
| winlist.py, wincmp.py | all fundamental d of a window; tier 4 against own exact M |
| logcheck.py, slabcheck.py | the 38 chunk logs (coverage, lane counts) and mjx_slabs.log (444 files) |
| perhsample.py, hkcheck.py | per-h table consistency, 127-h sample; HJKMP comparison |
| logs/ | isieve2 logs (c08 segments, windows), t4/adv/perh/small hcount outputs, PARI output, threshold rerun, STATE.md |

Process: every run under memcap 290 MB; peaks: isieve2 208–218 MB (c08), 140 MB (window 2), 52 MB (window 1, N = 3·10⁴);
hcount ≤ 41 MB; PARI 17 MB; threshold.py too short-lived for memcap's 2-3 s sampling (cert.py states < 100 MB). Downloads (postproc.gp, class_group_orders.txt) deleted.

RH relevance: N11

