UPDATED (7 records improve: N10 Mertens |M(x)| ≤ 9.69·10^4 x^0.9382, 1.52·10^4 x^0.9409; N10/N11 ExplProg exponent 7.8978 → 7.8956; N10 Ingham3 n^k thresholds, k = 3–6; N11 PrimRootNum exponents 26.85 / 23.69 / 7.896; N2 Λ ≤ 0.107550869585; N3 Li 106.2557 T² (routine); N10 conditional Ingham4 x^0.841. Unchanged at printed precision: Robin and Nicolas ranges, v_2(n) ≤ 72, odd Robin, Littlewood/class numbers, Kim–Lee wedge, Ingham4 x^0.856, zero density)

# HP141 — rerun of certificate-backed half-plane dependents with PR 141's θ₁₄₁

Orca agent hp141, 2026-11-14, from 19:13 UTC. θ₁₄₁ = exact rational of Riemann/Tree/R466.lean l.10–12 (127/128-digit numerator/denominator)
= 0.873345856106718576801184796065…, 3.611·10⁻⁵ below the old 0.87338196661320752. Helper: theta141.py (TH exact Fraction; TH_UP17/25/40
decimal roundings UP; TH_FLOAT_UP = smallest double ≥ TH). A zero-free half-plane Re s > θ stays true for larger θ, so every decimal here is
rounded up. Originals untouched; copies under HP141/<folder>/; every run under .scratch/memcap.sh 290 ≤1800 plus /usr/bin/time -l.

## Summary table (old = audited state, C3 = 9.4926, B = 0.87338196661320752; new = θ₁₄₁; peak = memcap top-MEM / time -l max RSS)

| folder | script | record | old | new | peak MB | s |
|---|---|---|---|---|---|---|
| MertensCH | cert_ch.py --bp (+ bprime.py) | N10 \|M(x)\| ≤ c x^θ, all x ≥ 1 | 9.77·10^4 (0.9382), 1.54·10^4 (0.9409) | **9.69·10^4, 1.52·10^4** | 51/69 | 93 |
| ExplProg | cert.py all | N10/N11 primes in progressions, Linnik, non-residues, L′/L | exponent 7.8978 | **7.8956** (all constants kept) | 14/33 | 23 |
| Ingham3 | cert.py | N10 prime in (n^k, (n+1)^k] | k = 3..6: e^627.01, e^117.21, e^57.75, e^36.25 | **e^626.82, e^117.17, e^57.73, e^36.23** (k = 7, all n: same) | 186/201 | 345 |
| PrimRootNum | cert.py all | N11 least prime primitive root | (1.956 log p)^26.86, ^23.70, ^7.898 | **^26.85, ^23.69, ^7.896** | 148/— | 320 |
| DBNB2 | make_flow_b2.py, verify_arcs_141.py | N2 de Bruijn–Newman Λ | 0.107559925318 | **0.107550869585** (exact arc verifier PASS) | 23/—; verifier 44/60 | 143 + 347 |
| Li200k | verify.py | N3 Li range (routine) | 106.245188 T² | **106.255746 T²** | 22/— | 79 |
| Ingham4 | cert.py part cs | N10 conditional (CS) short intervals | x^0.842; n₀(5, 4, 3) = e^45.584, e^91.1775, e^468.4233 | **x^0.841**; e^45.568, e^91.14, e^468.15 | 148/— | 502 |
| Ingham4 | cert.py part main | N10 prime in (x, x + x^0.856] | 0.856 (0.855 fails, R 1.0005) | 0.856 (0.855 fails, R 1.0018) | 89/104 | 113 |
| RobinCH | certify3.py 4; record_theta141.py | N4 Robin ≤ 10^(10^46.79), Nicolas p_k ≤ 1.78·10^48 | endpoints 10^48.252447 / 10^47.156206 | 10^48.252564 / 10^47.156287 (printed: same) | 19/37 | 85 + 212 |
| OddRobin | transfer.py | N4 odd Robin ≤ 10^(10^46.79) | 10^(10^46.7931) (1.43e47) | 10^(10^46.7940) (1.433e47; printed: same) | 21/40 | 132 |
| Robin2Adic | cert.py, check.py | N4 v_2(n) ≤ 72 | E* 8.7461707·10⁻²³, K = 72 | 8.7364609·10⁻²³, K = 72 | 25/40 | 37 + 191 |
| Littlewood78 | cert.py | N11 \|L(1, χ)\| ≥ 0.11/max(3, log log q); h(−d) thresholds | limit 0.11694; tails 0.115659…0.116933 | limit 0.11697; tails 0.115692…0.116966 (stated: same) | 16/33 | 216 |
| KimEff | certify.py | N6 Kim–Lee wedge, n ≥ 998,616,377,771, 6.45 n²/log² n | 3.6/(4a²) ≥ 6.4555883 | 6.4568372 (stated: same) | 0/33 | 7 |
| DensityB, B2, Q | verify.py | N7 density (formula already in COROLLARIES) | ALL PASS | ALL PASS | 13/28 | 8, 8, 32 |

## Draft RECORDS.md rows (orchestrator to place; nothing edited outside HP141)
- (N10, half-plane) Explicit power-saving Mertens bound: |M(x)| ≤ 9.69·10^4 x^0.9382 and ≤ 1.52·10^4 x^0.9409 for all x ≥ 1 (14 Nov: PR 141's
  θ₁₄₁ in place of 0.87338196661320752; was 9.77·10^4, 1.54·10^4). Script HP141/MertensCH/cert_ch.py --bp (+ bprime.py), log
  HP141/MertensCH/logs/cert_ch_bp.log; 200-bit twin logs/cert_ch_p200_bp.log.
- (N10/N11, half-plane) Explicit primes in progressions from Re s > θ₁₄₁: |ψ(x;q,a) − x/φ(q)| ≤ 1.5 x^θ₁₄₁ log²(qx) (x ≥ 2); least prime
  p ≡ a (q) ≤ (q log q)^7.8956, ≤ (0.72 q log q)^7.8956 for q ≥ 10⁴; a prime n ≤ (log q)^7.8956 with χ(n) ∉ {0, 1}; |L′/L(1, χ)| ≤
  7.8956 log log q + 10.53 (was 7.8978 everywhere, B = 0.8733819666132076). Script HP141/ExplProg/cert.py, log logs/cert_output.txt
  (320-bit twin logs/cert_prec320_output.txt). Also updates COROLLARIES l.35 (|L′/L|) and l.36 (Miller bases (log n)^7.8956).
- (N10, half-plane) Primes between consecutive powers: for every n ≥ e^626.82 (k = 3; was e^627.01), e^117.17 (k = 4; was e^117.21),
  e^57.73 (k = 5; was e^57.75), e^36.23 (k = 6; was e^36.25); k = 7 for every n ≥ 1 unchanged. Script HP141/Ingham3/cert.py, log
  logs/cert_run.log (C3 = 9.4926).
- (N2, half-plane) de Bruijn–Newman: Λ ≤ 0.107550869585 (Raymond's chain, Zenodo 23752103, restarted at PR 141's θ₁₄₁; was ours
  0.107559925318 at 0.87338196661320752). Scripts HP141/DBNB2/make_flow_b2.py + verify_arcs_141.py (author's exact arc verifier, flow
  streamed), logs logs/make_flow_theta141.log, logs/verify_arcs_141.log (PASS, 220,162 arcs); float model logs/float_check.log. N2 lives in
  TREE.md (no RECORDS row yet): the TREE N2 record line is the place, or a new RECORDS row in this wording.
- (N11, half-plane) Least prime primitive root: g*(p) ≤ (1.956 log p)^26.85 for every prime p ≥ 3 (was 26.86); (0.814 log p)^23.69 for
  p ≥ 5 with ω(p − 1) ≤ 10⁴ (was 23.70); (κ_ω log p)^7.896 (was 7.898). Script HP141/PrimRootNum/cert.py all (imports HP141/ExplProg/cert.py),
  log logs/cert_run.log. (PrimRootNum is not yet a RECORDS row; TREE/REPORT wording for the orchestrator.)
- Routine, not RECORDS (COROLLARIES/TREE wording): Li λ_n > 0 for n ≤ 106.2557 T² (Li200k, N3; the record is Fabbian's 77.03 T²);
  conditional Ingham4 x^0.841, n₀(5) = e^45.568, n₀(4) = e^91.14, n₀(3) = e^468.15 (on Chourasiya–Simonič).

## Spot-checks (task 4)
- Spot-check 1 (ExplProg, second precision): cert_prec320.py = cert.py with ctx.prec 320 instead of 200 (logs/cert_prec320_output.txt,
  17 s, 17 MB): no CERT FAIL; all 56 result lines structurally identical and every printed ball overlaps its 200-bit twin; decisive margins
  identical: (iii) 0.001213, (iv) 10.525, (ii) 0.007355 (0.72 q log q), (i) 1.48538 ≤ 1.5.
- Spot-check 2 (MertensCH headline, second precision): cert_ch_p200.py + bprime_p200.py = the same pipeline at 200 bits instead of 120
  (logs/cert_ch_p200_bp.log, logs/cert_ch_p200_results_bp.json; 44 s, 51 MB): θ = 0.9382 gives c = 969e2, log₁₀ c = 4.986323777050766,
  X0 = 10^242.80, X1 = 10^333.65, 1817 ranges, final ML092g8_e05_Bp — identical to the 120-bit run.
- Spot-check 3 (DBNB2): independent float model, see the DBNB2 section (agrees with the exact producer to 3·10⁻¹⁶).
- arith_check.py (arith_check.log): exact Fraction checks of every rounding in the drafted statements (exponents 7.8956, 26.85, 23.69,
  Ingham3 thresholds, 6.456, upward roundings of θ₁₄₁): ALL PASS.

## Per-folder results

Common: copies get STCAudit's C3 = 9.4926 correction (RECORDS audit note, 14 Nov) wherever HSW's misprinted 9.3675 appears, so "old"
is the audited current state. Direction: every script needs an UPPER bound for the zero-free abscissa (Re s > B zero-free; any B ≥ θ₁₄₁ is
valid); exact rationals are used as Arb balls arb(p)/q containing θ₁₄₁, decimals are rounded up.

### MertensCH (N10): |M(x)| ≤ c·x^θ for all x ≥ 1 — IMPROVES (9.77·10^4 → 9.69·10^4 at 0.9382; 1.54·10^4 → 1.52·10^4 at 0.9409)
- bprime.py with B′ = θ₁₄₁ rounded up at 10⁻¹⁷ (0.87334585610671858; every Mertens78 certificate translated by −(7/8 − B′), Decimal-exact
  shift 0.00165414389328142); cert_ch.py unchanged (imports this bprime; Mertens78 read-only; logs/JL0.json copied: σ_L = 0.85, B-free).
- cert_ch.py --bp 0.9409 0.9382 (logs/cert_ch_bp.log, 93 s, 51 MB; results logs/cert_ch_results_bp.json): θ = 0.9409: c = 152e2
  (log₁₀ 4.1818; was 154e2, 4.1875); θ = 0.9382: c = 969e2 (log₁₀ 4.9863; was 977e2, 4.9899); final certificate ML092g8_e05_Bp both times;
  c is printed rounded up (10^4.1818 = 15,199 ≤ 15,200; 10^4.9863 = 96,895 ≤ 96,900). Other θ of the original table not rerun (time).

### ExplProg (N10/N11 explicit progressions): IMPROVES (exponent 7.8978 → 7.8956)
- cert.py all with B = A(Fraction(θ₁₄₁)) (exact ball); the T* check is replaced by B < 0.87338196661320752 (logs/cert_output.txt: 23 s,
  14 MB, no CERT FAIL; compare.py → logs/compare.txt). 1/(1 − B) = 7.89551742454 (was 7.89781…, printed 7.8978); printed upward: 7.8956.
- All stated constants survive with B = θ₁₄₁: (i) max ratio 1.48538 ≤ 1.5 (x ≥ 2), 0.1731 (x ≥ 64), 0.0323 (x ≥ e^10); (ii) margins 0.2121
  (q log q, q ≥ 10²), 0.007355 (0.72 q log q, q ≥ 10⁴); (iii) worst log-margin 0.001213 > 0 (was 0.008514); (iv) 10.525 (was 10.526) ≤ 10.53.
- New statement: (i) |ψ(x;q,a) − x/φ(q)| ≤ 1.5 x^θ₁₄₁ log²(qx), x ≥ 2; (ii) least prime p ≡ a (q) ≤ (q log q)^7.8956 (q ≥ 2), ≤ (0.72 q log q)^7.8956
  (q ≥ 10⁴); (iii) a prime n ≤ (log q)^7.8956 with χ(n) ∉ {0, 1}; (iv) |L′/L(1, χ)| ≤ 7.8956 log log q + 10.53.

### Ingham3 (N10): a prime between n^k and (n+1)^k — IMPROVES for k = 3, 4, 5, 6 (k = 7 for every n ≥ 1 unchanged)
- cert.py with B = θ₁₄₁ rounded up at 10⁻¹⁷, assert B ≥ θ₁₄₁ ball (replaces B > T*), C3 = 9.4926 (logs/cert_run.log, logs/cert.log,
  logs/cert_result.json). Old = original log = STCAudit's 9.4926 rerun (identical thresholds).
- Main pipeline (KLN, RH to 3·10^12): k = 3: n₀ = e^627.0100 → e^626.8167 (x₀ = e^1881.03 → e^1880.45); k = 4: e^117.2100 → e^117.1700;
  k = 5: e^57.7460 → e^57.7300; k = 6: e^36.2417 → e^36.2283; k = 7: every n ≥ 100 certified (max R 0.3466, was 0.3482), n ≤ 100 direct.
- New statement (rounded up): a prime between n^k and (n+1)^k for every n ≥ e^626.82 (k = 3), e^117.17 (k = 4), e^57.73 (k = 5), e^36.23 (k = 6).

### PrimRootNum (N11): least prime primitive root g*(p) — IMPROVES (exponents 26.86 → 26.85, 23.70 → 23.69, 7.898 → 7.896)
- cert.py all, importing HP141/ExplProg/cert.py (B = exact θ₁₄₁ ball); sieve.py, tail.py, small.py, grid.py, logs/sieve_best.json copied
  (try 1 failed on the missing small.py: logs/cert_run_try1_missing_small.log). logs/cert_run.log: 320 s, 148 MB; kappa_table.txt, tail_output.txt.
- final: g*(p) ≤ (489/250 log p)^26.844759 (was 26.852415) ⇒ (1.956 log p)^26.85 for every p ≥ 3 (binding: tail, margin 3.7·10⁻⁴);
  Theorems B/B′: (1.4312 log p)^(3/(1−B)) and (0.814 log p)^(3/(1−B)), 3/(1−B) = 23.686552 (was 23.693307) ⇒ 23.69; small primes ≤ 10⁷ rechecked
  with exponent 23.686552 (worst ratio 0.1083). Theorem A: (κ_ω log p)^7.8956; 407 of the 10⁴ tabulated κ_ω rise by one unit in the 3rd
  digit (e.g. κ₄₄ 4650 → 4660), each with a smaller total bound for log p ≥ log p_ω. Corollaries g*(p) < p (p ≥ 10^41.032) and < √p − 2
  (p ≥ 10^97.547): unchanged.

### Li200k (N3, outside the brief's rg list: it writes a = B − 1/2 = 0.37338196661320752): Li range — IMPROVES vs the audited value
- verify.py with variant F replaced by G: a = θ₁₄₁ − 1/2 exact (arb ball), C3 = 9.4926 (logs/verify.log, 79 s, 22 MB).
- G: λ_n > 0 for all n ≤ 956,413,498,941,653,124,095,410,176 = 106.255746 T² (T = 3,000,175,332,800); F with the old a and C3 = 9.4926
  (STCAudit/dep_Li200k_94926.log): 956,318,469,361,119,813,265,522,688 = 106.245188 T² (TREE N3 prints 106.25, made with 9.3675; COROLLARIES
  R466: "≈ 106.24 T², not recertified"). Variants A, B, C (B-free) reproduce 78.713642, 79.738933, 105.774404 T².
- Routine (COROLLARIES/TREE N3 wording), not RECORDS: Fabbian's 77.03 T² (RH to height T only) stays the record.

### DBNB2 (N2): de Bruijn–Newman Λ ≤ 0.107559925318 — IMPROVES to Λ ≤ 0.107550869585 (producer + exact arc verifier PASS)
- make_flow_b2.py (Raymond's 67-row flow, exact FLINT rationals, author's 2^−96 outward rounding) with B0 = exact θ₁₄₁, inputs hash-checked
  from Zenodo 23752103 (Proof_and_Replay.zip, sha256 bf5c7b3d… OK; only the plan, the source and ACTUAL_67_ROW_ARC_INPUTS.json kept, in
  .scratch/hp141/dbnb2): logs/make_flow_theta141.log (143 s, 23 MB), logs/FLOW_THETA141_SUMMARY.json.
- Endpoint T + q_T/2 = 0.1075508695842417 (220,162 arcs, T = 110081/1024000, q_T = 9.9786·10⁻⁵; headline ceiling 21510173917/2·10^11 =
  0.107550869585) vs 0.10755992531713451 (220,181 arcs): gain 9.056·10⁻⁶, slope 0.2508 (REPORT's dΛ/db0 = 0.2507); startup gap
  9.3096·10⁻¹³ > 5·10⁻¹³ and kernel cap margin positive (both in the summary).
- The author's arc verifier as patched by DBNB2 (verify_arcs_b2.py) logged 676 MB max RSS: rerun streamed instead (UPDATE below).
- Independent float model (float_check.py, from the TeX formulas, not the producer; logs/float_check.log, 24 s, 9 MB), b0 = smallest
  double ≥ θ₁₄₁: scheme 1 (certificate model) endpoint 0.10755086958424454 with the same 220,162 arcs (producer: 0.1075508695842417,
  Δ = 3·10⁻¹⁶); scheme 2 (RK4) 0.1075508349; old b0 rows reproduce 0.10755992531713356 / 0.1075598907.
- UPDATE 20:11 UTC — exact arc verifier PASS, so DBNB2 is certified: patch_verifier141.py = DBNB2's patch_verifier.py (mode A: same
  constant swaps, read from the θ₁₄₁ flow; initial square (2θ₁₄₁ − 1)² + 10⁻¹² exact) plus one I/O-only change in the author's
  verify_arcs_fresh.py: the flow is streamed record by record (StreamFlow) instead of json.loads of all 220k records. verify_arcs_141.py
  (copied here; logs/verify_arcs_141.log, 347 s, memcap 44 MB, max RSS 60 MB vs 676 MB): "PASS": true, 220,162 arcs, 5,776,946 exact
  checks, endpoint 2130264443492687883657300572721/19807040628566084398385987584000 = 0.1075508695842417 < headline 21510173917/2·10^11
  (= 0.107550869585) < endpoint + 10⁻¹²; flow sha256 ea2ab380…, witness trace sha256 88ca31c0… (423 MB, deleted; regenerable);
  logs/EXACT_ARC_VERIFICATION_FRESH_V1.json kept; verifier141.diff = author's verify_arcs_fresh.py → verify_arcs_141.py (constants
  from the flow, mode-A tube switch, StreamFlow; no mathematical check touched). Not rerun: mode B (tube kept), compact_arc_witnesses.py and the author's 15-stage
  replay; per DBNB2/REPORT the other stages do not involve q or b0 (source certificates q-free; kernel/own-pair certificates hold on [0, cap]).

### Ingham4 (N10): a prime in (x, x + x^0.856] for every x ≥ 1 — UNCHANGED; conditional (CS) part IMPROVES to x^0.841
- cert.py part main with B = θ₁₄₁ rounded up at 10⁻¹⁷, assert B ≥ θ₁₄₁ ball (replaces B > T*), C3 = 9.4926, and the "0.855 must fail"
  assert turned into a report (logs/part_main_run.log, 113 s, 89 MB; logs/part_main.json). θ = 0.856: CERTIFIED, max R 0.9673 at
  log x = 191.82 (was 0.9711; tail 0.0141). θ = 0.855: still FAILS (R = 1.0018 at log x = 187.74; the C3 correction and θ₁₄₁ act in opposite
  directions). Parts a0, delta1, cs (conditional CS thresholds n₀(k), x^0.842) not rerun.
- Part cs (conditional on the unrefereed Chourasiya–Simonič density, its constant corrected; "0.841 must fail" assert turned into a report;
  logs/part_cs_run.log, 502 s, 148 MB): θ = 0.841 now CERTIFIED (max R 0.9975 at log x = 174.20; was FAILS, R = 1.0005) ⇒ conditional
  x^0.841 (was x^0.842); k = 6 every n ≥ 1 unchanged; n₀(5) = e^45.584 → e^45.568, n₀(4) = e^91.1775 → e^91.14, n₀(3) = e^468.4233 → e^468.15.
  Conditional, so not a RECORDS row (TREE/Ingham4 REPORT wording for the orchestrator). Part a0/delta1 (comparison variants) not rerun.

### RobinCH (N4): Robin 5041 ≤ n ≤ 10^(10^46.79), Nicolas p_k ≤ 1.78·10^48 — UNCHANGED at printed precision
- certify3.py 4 <θ₁₄₁ p/q> (endpoint search, "Level 3" = code level 4; logs/endpoint_level4_theta141.log, 85 s, 19 MB):
  Nicolas 10^48.252447 → 10^48.252564, Robin (largest prime of CA) 10^47.156206 → 10^47.156287; old b0 rerun in the same 9.4926 copy
  (logs/endpoint_level4_old.log) reproduces the original endpoints exactly.
- record_theta141.py (claims as STCAudit/dep_robinch_record.py; logs/record_theta141.log, 212 s, 19 MB): θ₁₄₁ passes 1.788e48 / 1.433e47
  ⇒ Robin to 10^(10^46.7940) (old 4-digit claims 1.788e48 / 1.432e47 ⇒ 10^(10^46.7937)); the printed 3-digit claims pass unchanged.
- Gain: +1.2·10⁻⁴ (Nicolas) and +8·10⁻⁵ (Robin) in log₁₀ x; at the record's 3-digit/2-decimal precision nothing moves. Optional finer
  RECORDS wording: 10^(10^46.7937) → 10^(10^46.7940), Nicolas p_k ≤ 1.788·10^48 (unchanged). Recommendation: no RECORDS change.
- code level 3 (wrong level, aborted; logs/endpoint_codelevel3_theta141_aborted.log) kept only as a trace.

### OddRobin (N4): odd n, N < n ≤ 10^(10^46.79) with the half-plane — UNCHANGED at printed precision
- transfer.py (default b0 = exact θ₁₄₁ p/q; imports HP141/RobinCH, C3 = 9.4926) 5e8 1.433e47 (logs/transfer_theta141_1.433e47.log,
  132 s, 21 MB): PASS, 1288 boxes, min margin 2.6·10⁻⁷ ⇒ all odd N < n ≤ 10^(10^46.7940) (original: 1.43e47 ⇒ 10^(10^46.7931); the old
  b0's Robin endpoint 1.4329e47 allows at most 1.432e47 ⇒ 10^(10^46.7937)). Printed 10^(10^46.79): unchanged.

### Robin2Adic (N4): Robin for n > 5040 with v_2(n) ≤ 72 — UNCHANGED
- cert.py with B_BEBEL = θ₁₄₁ rounded up (0.87334585610671858), B_UP dyadic above it, assert B ≥ θ₁₄₁ ball replaces B ≥ T*; HSW 9.4926;
  imports HP141/Ingham4/cert.py (KLN rows only, B-independent). logs/cert_run.log: 37 s, 23 MB.
- KLN Lemma 4.14, RH(3·10^12): E* ≤ 8.7364609·10⁻²³ = 2^−73.2773 → K = 72 (was 8.7461707·10⁻²³ = 2^−73.276, K = 72); CS variants K = 75
  (unchanged). K = 73 needs E* < 2^−74 (k_of_L.py: L ≥ 111.84); the half-plane is not the binding input (RobinCH's L = 108.57 is).

### Littlewood78 (N11): |L(1, χ)| ≥ 0.11/max(3, log log q), class numbers — UNCHANGED (limit constant only, 4th digit)
- cert.py, B = θ₁₄₁ rounded up at 10⁻²⁵ (0.8733458561067185768011848; arb ball above θ₁₄₁), cert_output.txt, logs/cert_run.log (216 s, 16 MB);
  compare.py lists every changed key line against the original.
- Unchanged: claim (A) 0.11/max(3, log log q); LamB = 3.13 for 0.115 (0.01 grid) and LamB2 = 4.0 for 0.116 (0.05 grid); every class-number
  threshold (h(−d) > 10⁴ for d ≥ e^27.36, …; 0.01 cells), so the N11 class-number record (h ≤ 10,000) is untouched (its search bound binds).
- Moved: limit (π²/6)e^−γ(1 − B) 0.1169398 → 0.1169732 (printed 0.11694 → 0.11697); tails 0.115659 → 0.115692 (log log q ≥ 40),
  0.116424 → 0.116457 (≥ 100), 0.116933 → 0.116966 (≥ 10⁴); upper side e^γ/(1 − B) 14.06650 → 14.06249. Tightest (B)-cell margin
  4.7·10⁻⁴ → 7.3·10⁻⁴: 3.12 still fails the 0.01 grid (the search starts at the first passing multiple).

### KimEff (N6): Kim–Lee wedge, J^{d,n} hyperbolic for n ≥ 998,616,377,771, d ≤ 6.45 n²/log² n — UNCHANGED at printed precision
- certify.py with a = θ₁₄₁ − 1/2 rounded up at 10⁻²⁵ (0.3733458561067185768011848; a = half-width of Ξ's zero strip, larger a weaker),
  logs/certify.log (7 s): ALL CHECKS PASS; n_* line, c = 0.55 constant C = 0.87160 and the two-thirds containment (d > 1.0471·10^8) identical.
- Moved: 3.6/(4a²) ≥ 6.4555883 → 6.4568372 (stated constant 6.450 kept; 6.456 would now be available, not adopted); crossover
  D_n > D_all from n ≥ 4.6974173·10^13 → 4.6969477·10^13; Phi_up, E_up at the 4th digit.
- Related formula-only Jensen update (Chasse sector lemma, d ≤ 1.61440·10^25 for all n) is already in COROLLARIES (R466).

### DensityB (N7 zero density), DensityB2, DensityQ: exact Fraction checks with B = θ₁₄₁
- DensityB/verify.py: 24 PASS, ALL PASS (the "B above Liu's contact" check is inverted on purpose: θ₁₄₁ < B_2, PR 141). The record's
  formula is B-generic; its θ₁₄₁ values (σ* = 0.856069, A(0.87) ≤ 1.692204) are already in COROLLARIES (R466): UNCHANGED here.
- DensityB2/verify.py (NO RECORD folder): 73 PASS, ALL PASS.
- DensityQ/verify.py: ALL PASS (32 s, 13 MB). Its record (4(1−σ)/(6σ−3) on (5/6, 7/8]) needs no half-plane (REPORT l.18): UNCHANGED.

### Skipped (one line each; lesson 157: logged peaks read first)
| folder/script | reason |
|---|---|
| DBNB2/evidence/patch_verifier.py (as is) | its verifier logged 708,952,064 B max RSS (evidence/A_verify_run.log) > 290 MB; replaced by patch_verifier141.py (same swaps + streamed flow, 60 MB; DBNB2 section) |
| STCAudit/dep_robinch_record.py | audit rerun of RobinCH's claims; superseded here by HP141/RobinCH/record_theta141.py (same code path, θ₁₄₁) |
| STCAudit/dep_robinheight_record.py, RobinHeight/check_all.py | RobinHeight's record (10^(10^38.60)) is superseded by RobinCH's 10^(10^46.79); no record can move |
| RobinCH/check_all2.py | parts (2)/(3) = record_theta141.py's claims; Level 2-CH is superseded by Level 3; other parts half-plane-free |
| OddRobin/check_all.py | gap.py, dfs_gap.py, oddca.py do not use the half-plane; its transfer.py is rerun (below) |
| R2AAudit/falsify.py | audit falsification of Robin2Adic's lemmas (BB only sets test ranges); not a record; K unchanged (Robin2Adic rerun) |
| N8Fam/cfam.py | BQ only in the 'Q' form-factor model, which Johns's F₁ supersedes (REPORT l.96); the conditional record uses F₁ |
| Harvest6/k_est.py | harvest float estimate, superseded by Robin2Adic/cert.py (rerun) |
| Exp48/cert.py | NO RECORD (known, Jaswal); B only in side check (2d) 1/(4(1−B)) < 11.2, true a fortiori (θ₁₄₁: 1.9739) |
| BVIrr/params.py | comparison row a = 1/(1−T) only (θ₁₄₁: 7.8955 vs 7.8978); the result uses Cureton–Thorner's 31/32 |
| DensityB/repro/antedb_check.py, explore.py | float exploration / needs an external ANTEDB clone; DensityB's record is certified by verify.py (rerun) |
| B2Below/infimum.py, levers.py | mechanism analysis of PR 137's own parameter point (2183367/2500000); no record depends on it |
| MertensCH/bprime.py (main mode) | its module is what cert_ch.py --bp imports (rerun with θ₁₄₁); the stand-alone Mertens78-only table is routine |
| ExplProg/{sanity,lsanity,ksanity,xcheck}.py | float sanity scripts; the certificate is cert.py (rerun) |

## Footprints (every run; .scratch/memcap.sh 290 MB cap, one process at a time, single thread)
| log | memcap peak MB (top MEM) | max RSS MB (time -l) | wall s | rc |
|---|---|---|---|---|
| DBNB2/logs/float_check.log | 9 rc=0 | 17 | 24.38 | 0 |
| DBNB2/logs/make_flow_theta141.log | 23 rc=0 | 39 | 143.30 | 0 |
| DBNB2/logs/patch_verifier141.log | 0 rc=0 | 46 | 8.28 | 0 |
| DBNB2/logs/verify_arcs_141.log | 44 rc=0 | 60 | 347.14 | 0 |
| DensityB/logs/verify.log | 0 rc=0 | 28 | 7.54 | 0 |
| DensityB2/logs/verify.log | 0 rc=0 | 28 | 7.53 | 0 |
| DensityQ/logs/verify.log | 13 rc=0 | 27 | 31.61 | 0 |
| ExplProg/logs/cert_output.txt | 14 rc=0 | 33 | 22.89 | 0 |
| ExplProg/logs/cert_prec320_output.txt | 17 rc=0 | 36 | 17.91 | 0 |
| Ingham3/logs/cert_run.log | 186 rc=0 | 201 | 345.06 | 0 |
| Ingham4/logs/part_cs_run.log | 148 rc=0 | 152 | 502.38 | 0 |
| Ingham4/logs/part_main_run.log | 89 rc=0 | 104 | 113.43 | 0 |
| KimEff/logs/certify.log | 0 rc=0 | 33 | 6.72 | 0 |
| Li200k/logs/verify.log | 22 rc=0 | 37 | 83.80 | 0 |
| Littlewood78/logs/cert_run.log | 16 rc=0 | 33 | 215.78 | 0 |
| MertensCH/logs/cert_ch_bp.log | 51 rc=0 | 69 | 93.15 | 0 |
| MertensCH/logs/cert_ch_p200_bp.log | 51 rc=0 | 69 | 44.10 | 0 |
| OddRobin/logs/transfer_theta141_1.433e47.log | 21 rc=0 | 40 | 131.61 | 0 |
| PrimRootNum/logs/cert_run.log | 148 rc=0 | 159 | 320.01 | 0 |
| PrimRootNum/logs/cert_run_try1_missing_small.log | 84 rc=1 | 118 | 87.87 | 1 |
| Robin2Adic/logs/cert_run.log | 23 rc=0 | 37 | 36.77 | 0 |
| Robin2Adic/logs/check_run.log | 25 rc=0 | 40 | 191.00 | 0 |
| RobinCH/logs/endpoint_level4_old.log | 19 rc=0 | 34 | 70.53 | 0 |
| RobinCH/logs/endpoint_level4_theta141.log | 19 rc=0 | 37 | 84.61 | 0 |
| RobinCH/logs/record_theta141.log | 19 rc=0 | 34 | 212.13 | 0 |

## Not done (precisely)
- MertensCH: only the two record rows θ = 0.9409, 0.9382 rerun; the other 18 θ of the original table and the T′ = 1000 baseline not rerun.
- Ingham4: parts a0 and delta1 (comparison variants) not rerun. DBNB2: mode B (tube kept), compact_arc_witnesses.py and the author's
  15-stage replay not rerun (b0-free stages). No reviewer verdict was requested (records, no discovery); the streamed-verifier change
  (DBNB2/patch_verifier141.py, I/O only) is the one place a reviewer could look before the N2 line is placed.
- Files: copies with patch.diff per folder, theta141.py, footprints.py/.md, arith_check.py/.log; downloads (Zenodo zip, flow, 423 MB
  witness trace) deleted; .scratch/hp141 holds only STATE.md and the run/queue/patch scripts.

RH relevance: N7 (records of dependents)
