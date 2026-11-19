AUDIT PASSED (one erratum): every N5 sign and value claim holds for all N ≤ 1280 — a certified rerun of the final code (memory-lean copy, under the cap) and independent mpmath, direct-integration and double checks reproduce it — but "ball width < 10⁻²³" is false for 62 N (max 3.14·10⁻²³ at N = 1068; true bound < 3.2·10⁻²³), which changes no conclusion.

# nbaudit: line-by-line audit of the N5 record (RECORDS.md N5 line, Riemann/Deep/NBRate/)

Audited claim (RECORDS.md l.31): "d_N² for every N ≤ 1280 (Vasyunin's Gram summed exactly, ball width < 10⁻²³); d_N² log N − C (C = 2 + γ − log 4π) changes sign 15 times on [2, 1280], first below C at N = 41, below on [584, 1280] (was N ≤ 128 rigorous, Zenodo 24088197; d_300² interval, Zenodo 21898513; 2611.03924 in double precision)." Chain with source lines: [FULLPROOF.md](FULLPROOF.md).

## Items
| # | Item | Verdict | Evidence |
|---|---|---|---|
| 1 | d_N² = Báez-Duarte's distance in L²(0,∞), c free | holds | BD math/0202141 l.44, 56-66 (H = L₂(0,∞), ρ_a = ρ(1/(ax)), χ = 1_(0,1]); Burnol l.168-176 same space; code = V1 |
| 2 | Vasyunin Gram closed form | holds | direct exact integration (no cotangents) = Vasyunin to 2·10⁻³⁹ on all 210 pairs ≤ 20 (incl. non-coprime: homogeneity) and 24 random coprime pairs ≤ 640 (b2); Zenodo 21898513's digamma Gram (V2 form, n ≤ 300) = ours to 8.6·10⁻¹⁵ (b3) |
| 3 | "summed exactly" | holds | cot(πm/k) as 220-bit Arb balls → 2⁻¹⁴⁰ integers, cotangent sums accumulated exactly (integer-valued float products, partial sums < 2⁴⁷), error ≤ h·(max rad + 2⁻¹⁴¹), h = ⌊(k−1)/2⌋ (s2 l.28-55) |
| 4 | Linear algebra | holds | d_N² = 1 − bᵀG_N⁻¹b; any c: d_N² = Q(c) − rᵀG_N⁻¹r ∈ [Q − ‖r‖²/λ_min, Q]; exact integer Q and r, ε_G, ε_b propagated (s2 l.140-180) |
| 5 | λ_min(G_N) ≥ 4.768·10⁻⁷, N ≤ 1280 | holds | Weyl + interlacing; float error of E: coefficients 3u, 9u ≤ 12u used; reviewer CORRECT; rerun reproduces ‖E‖_F ≤ 5.487·10⁻¹⁴ |
| 6 | ball width < 10⁻²³ | **fails as written** | exact widths of all 1280 intervals (c1, c2): max 3.138·10⁻²³ (N = 1068), 62 N ≥ 10⁻²³ (first N = 225), median 1.1·10⁻²⁴. NBRate's "≤ 6.4·10⁻²⁴ for every N" extrapolated its 18 logged N (those widths reproduce exactly). Correct: < 3.2·10⁻²³ |
| 7 | C = 2 + γ − log 4π | holds | = Σ_{Re ρ=½} m(ρ)/|ρ|² under RH (BCF l.79); role: Burnol liminf d_N² log N ≥ Σ m_ρ²/|ρ|² (Bu l.204) ≥ C under RH, and the BBLS/BCF asymptotic d_N² ~ C/log N (BCF l.75-83); the finite sign pattern needs no hypothesis; Arb value 0.04619141793 |
| 8 | 15 sign changes on [2,1280] | holds | certified at N = 41, 49, 206, 213, 220, 283, 287, 372, 381, 383, 400, 412, 418, 525, 584 (c2); every N decided, tightest margin 2.045·10⁻⁶ (N = 411) |
| 9 | first below C at N = 41 | holds | certified + at every N ∈ [2,40], − at 41: d_41² ∈ 0.012326491948693913786970[00, 11], so d_41² log 41 − C = −4.1610·10⁻⁴ (width 3.5·10⁻²⁵) |
| 10 | below C on [584,1280] | holds | certified − at all 697 N; on (640,1280] margin ≥ 6.80·10⁻⁴ (minimum at N = 1280) |
| 11 | every N ≤ 1280, no gaps | holds | loop N = 1..NMAX, no undecided N; our exact interval file has all 1280 N (c2) |
| 12 | "was" clause | imprecise | Zenodo 24088197 is withdrawn by its author ("SUPERSEDED & WITHDRAWN", files restricted): unverifiable. Zenodo 21898513 and 2611.03924 compute the CONSTRAINED distance (Σc_k/k = 0, L²(0,1), = NBRate's V2 ≥ V1), not V1: their values agree with our V2 (below) but are not the same quantity |

## Falsification runs (all under 300 MB at the end; logs/)
| Run | What | Result |
|---|---|---|
| a1 + a2 | verbatim s1/s2 copies, OUT → /tmp, NMAX 640 (lesson 148); a2 adds exact-width tracking and an exact lo/hi dump | all 16 logged rows identical to NBRate's 1280 log to every printed digit; same 16 sign runs on [2,640] (15 changes); peak 114 MB top, 142 MB RSS |
| c1 | memory-lean copy of s2 (reused buffers, int32 limbs, μ = 2⁻²¹ given, blockwise λ certificate, ramp skipped); byte-identical lo/hi to a2 at 640 | full certified N ≤ 1280: λ_min ≥ 4.768371·10⁻⁷ (‖E‖_F ≤ 5.487·10⁻¹⁴, as logged), the 18 logged widths reproduce exactly, no undecided N, same 16 runs; 276 MB RSS, 242 MB top |
| a3 | verbatim s3_v2 copy (certified V2), NMAX 640 | V2 at n = 40, 80, 160, 320, 640 identical to NBRate's certified values (e.g. 0.007023264380134277 at 640); 106 MB RSS |
| b1 | own mpmath code from Vasyunin's formula (full cot sums, no symmetry trick, 60 digits), N ≤ 100 | 100/100 values inside the certified intervals; 2 sign changes on [2,100] (41, 49); V2 at n = 40, 80 equals NBRate's certified V2 to all 16 printed digits |
| b2 | Gram by direct exact integration ∫{u/h}{u/k}u⁻²du (pieces + Hurwitz-ζ tail), no Vasyunin | entries = Vasyunin to 2.3·10⁻³⁹ (pairs ≤ 20), 3.5·10⁻³⁹ (24 pairs ≤ 640), = s1's double Gram to 1.4·10⁻¹⁷; d_N², N ≤ 20, inside the certified intervals |
| b3 | own double sweep, all N ≤ 1280 (not rigorous) | same 16 runs; max deviation from certified midpoints 5.9·10⁻¹⁴; d_N² non-increasing; |Δ²(d_N² log N)| ≤ 6.1·10⁻⁵ on (640,1280] (smooth) |
| b4 | BCF ramp I_N, own double code, N ≤ 1280 | within 10⁻¹³ of a2's certified ramp bounds (N ≤ 640); window means 3.5480, 3.5496, 3.5478 as in COROLLARIES l.51 |
| c2 | checks on the 1280 exact intervals | no gaps; 981 of 1279 N below C (as reported); lo_{N+1} ≤ hi_N for all N; 1280- and 640-run intervals intersect at every N ≤ 640; all 16 REPORT table rows agree to their printed digits |

Comparisons (V2 = constrained distance, the quantity the earlier works computed):
| n | ours V2 (b3 double; b1 mpmath) | NBRate s3 certified V2 | 2611.03924 Table tab:dn | Zenodo 21898513 |
|---|---|---|---|---|
| 40 | 0.01273583539133116; 0.012735835391331124 | 0.01273583539133112 | 1.274e−2 | – |
| 160 | 0.009356604331762 | 0.009356604331758 | 9.357e−3 | – |
| 300 | 0.008002939956659 (V1: 0.007994899430933) | – | – | [8.00293986, 8.00294005]·10⁻³ ∋ ours |
| 640 | 0.007023264380137 | 0.007023264380134 | 7.023e−3 | – |
| 1280 | 0.006366197363451 | 0.006366197363433 | 6.366e−3 | – |
Zenodo 21898513's Gram cache (n = 2..300) equals our V2 Gram TᵀGT to 8.6·10⁻¹⁵. Zenodo 24088197: withdrawn, files restricted.

## What stays unrechecked
- Vasyunin's formula itself (not its implementation, which b1/b3 confirm) for all 1280² pairs: checked to 40 digits on 234 pairs (b2) and to 10⁻¹⁴ on all pairs ≤ 300 (Zenodo 21898513's digamma Gram, V2 form; R461: Zenodo 22831678's digamma form to 10⁻¹⁶); larger pairs rest on Vasyunin's theorem (1995, not read, no TeX).
- BCF-ramp bounds (COROLLARIES l.51, "3.5478 on [800, 1280]"): certified only to N ≤ 640 here (a2; c1 skipped them); b4 (own double code) lies within 10⁻¹³ of a2's certified ramp at all N ≤ 640 (max 5.6·10⁻¹⁵) and reproduces the window means 3.5480 [200,400], 3.5496 [400,800], 3.5478 [800,1280] and I_1280 = 0.0756 ≈ 12 d²_1280 (double, not certified above 640). Certified V2 (s3_v2) rerun verbatim at NMAX 640 (a3): n = 40, 80, 160, 320, 640 reproduce to every digit; its n = 1280 value only checked in double (b3, 1.8·10⁻¹⁴).
- IEEE semantics of the BLAS/LAPACK used (exact integer products < 2⁵³ assumed of numpy matmul): checked indirectly (b1 containment at N ≤ 100, b3 agreement).

## Process notes
- MEMORY: my first b3 attempt was killed by memcap at 543 MB (top); its debug rerun reached 711 MB ru_maxrss while memcap reported 20 MB (memcap under-measures some uv-launched runs; /usr/bin/time and in-script maxrss are reliable). c1's first 1280 attempt hit my 280 MB in-process watchdog (323 MB RSS, 241 MB top). Cause in all cases (and of NBRate's 567 MB): per-iteration O(N²) temporaries (outer products, solve_triangular's slice copy) retained by macOS malloc; reused buffers fix it with byte-identical output.
- Reviewer (λ_min float-error step, pasted inline with exact formulas): first run timed out at 170 s with recap CORRECT (coefficients ≈ 3u on |G_f| + μ, 9u/(1−8u) on ‖Y_i‖‖Y_j‖; Y exact); shorter rerun: **CORRECT** (conservative 10u/(1−9u) < 12u; underflow absorbed by 12uμ). Robustness: ‖E‖_F = 5.5·10⁻¹⁴ ≪ μ = 4.8·10⁻⁷, so even a 10⁶-fold error in it leaves λ_min ≥ 4.2·10⁻⁷ and every sign unchanged.

## What survives, dependents (rg -n -i 'nbrate|NBRate|d_N²|1280' TREE.md COROLLARIES.md RECORDS.md DISCOVERIES.md NOGO.md discoveries/)
- Survives: every d_N² digit, the certified sign of d_N² log N − C at every N ∈ [2,1280], 15 changes, first below at 41, below on [584,1280], 981/1279 below; width bound corrected to < 3.2·10⁻²³ (immaterial: tightest sign margin 2·10⁻⁶).
- RECORDS.md l.31 (N5): erratum "ball width < 10⁻²³" → "< 3.2·10⁻²³"; the "was" clause should say the earlier results are for the constrained distance (V2) and that Zenodo 24088197 is withdrawn (orchestrator: we do not edit RECORDS).
- NBRate REPORT "Resulting ball widths: ≤ 6.4·10⁻²⁴ for every N ≤ 1280 (V1)": false, max 3.14·10⁻²³ (V2's ≤ 8.7·10⁻²⁴ covers only its 6 n).
- TREE.md l.9 (N5: 15 changes, 41, [584,1280], d_640² log 640 = 0.0453405057), TREE.md l.497 + Riemann/Tree/R461.lean l.10 (d_100² = 0.0101887625, d_320², d_640²), ExploreChristoffel REPORT l.31 (12-digit agreement): hold.
- COROLLARIES.md l.51 ("certified Gram values give 3.5478 on [800, 1280]"): Gram confirmed; the three window means reproduce (b4, double; certified ramp agreement to N ≤ 640). NB2nd, NOGO l.603/628: own computations, no use of the width.
- New artefact: logs/dN2_certified_1280.txt, the certified enclosure of d_N² for every N ≤ 1280 (NBRate's rig_V1.npz is gone with .scratch/nbrate).

## Files
Riemann/Deep/NBAudit/: REPORT.md; FULLPROOF.md (chain with source line numbers); code/ a1_s1copy.py, a2_s2copy.py, a3_s3copy.py (verbatim copies, paths → /tmp/nbaudit), c1_s2lean.py (lean certified copy), b1_vas_mp.py, b2_direct.py, b3_dbl.py, b4_ramp.py, c2_check.py, c3_table.py; logs/ (run logs with /usr/bin/time and memcap peaks; dN2_certified_1280.txt). Run order: a1 640 → a2 640 → a3 640 → b1 100 → b2 → b3 1280 → c1 1280 −21 → c2 → c3 → b4 1280 (uv run --no-project --with python-flint --with numpy --with scipy --with mpmath).

RH relevance: N5
