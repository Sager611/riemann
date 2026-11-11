IDEA: Chekmachev's Woodbury–Cholesky matrix is an N-body interaction matrix (kernel 1/(y − y′), particles at y = k², 4j², two charges each: displacement rank 2), so a structured LDLᵀ with a backward-error bound decides his window certificates in O(n) memory: it re-certifies the N1 record c = 47 for all f below 125 MB per process (dense: ≈ 35 GB), while H-matrix compression gains at most ≈ 3× at certificate accuracy (node N1).

# ExploreFMM — N-body / hierarchical-matrix tools → Weil positivity certificates (explore28, 11 Nov 2026)

RH relevance: N1 (Weil positivity on windows; the c = 47 record's matrix stage, blocked by memory in WeilZeta, WeilSimple and WeilDefl).

## 1. The transfer
| tool from N-body / boundary-integral computing | what it becomes for Chekmachev's certificate | result here |
|---|---|---|
| particle representation: a kernel matrix is never stored, only points and charges | the whole Woodbury + Cholesky stage is the 1/(y − y′) interaction of n particles at y = k², 4j², y_v with two charges each (Lemma P): displacement rank 2 | **exact, O(n) memory** |
| structured elimination on the charges (Schur/GKO recursion, the fast-direct-solver view of FMM matrices) | inertia of J″ = G + diag(δ) decides Id − W > ν (Lemma I) | works with pivoting + generator balancing |
| (c) a-posteriori certification of a compressed factorisation | backward-error bound accumulated per step, O(n) (Lemma B, after Wiktorowicz's point LDL) | rigorous; loses 17–33 digits instead of ~300 |
| (b) H/HODLR matrices (numerical low rank) | off-diagonal ranks measured at certificate accuracy | **barrier**: 1.1–1.2× less storage at c = 25, 29; ≈ 3× (≈ 7 GB) at c = 47 |
| (a) FMM for the Cauchy sums F(y) = Σ_G γ_j/(y − 4j²) | 4.6·10⁹ terms at c = 47; a factor-16 scale-separated expansion leaves 23.7 % of the pairs | not built (≤ 4×); an exact-integer kernel (Lemma F) gives 5.7× |
Proofs: DERIVATION.md. Prior art: the displacement-rank-2 form H_qr = (u_q v_r − v_q u_r)/(λ_q − λ_r) and a structured point LDL with a proved backward-error bound are the core of Wiktorowicz–Wiktorowicz's c ≈ 100 pipeline (their x-space Galerkin matrix; ../WiktAudit/REPORT.md rows 3, 5, 10). New here: the same structure in Chekmachev's lattice-minorant matrix (P0⁻¹-Gram of e_k, φ_j, v) and its inertia form, which removes his 35 GB wall.

## 2. Statements (DERIVATION.md)
- **Lemma P (particle form).** With F(y) = Σ_{j∈G} γ_j/(y − 4j²), γ_j = (1/dl_j − 1/deff)A_j²: e_k ↦ (y, u, w) = (k², β_k, β_kF(k²)); v ↦ (y_v, β_v, β_vF(y_v) + κ/deff); φ_j (dl_j ≤ 0) ↦ (4j², 0, A_j/deff). Then ⟨c_a, P0⁻¹c_b⟩ = (u_aw_b − w_au_b)/(y_a − y_b) for a ≠ b.
- **Lemma I (inertia).** J″(ν) = G + diag(δ), δ_U = 1/D_a, δ_B = −(1 − ν)/s_a: Id − W > ν·Id ⟺ In J″(ν) = (n_U, n_B, 0) (Haynsworth; congruences without square roots).
- **Lemma B (backward-error certificate).** The recursion run from exact points at each step, with local errors bounded by radii and H_a = Σ_b 1/|y_a − y_b|, gives η ≥ ‖A_0 − J_s‖ for a matrix A_0 whose exact pivots are the computed ones. #pos = n_U ⟹ λ_min(Id − W) ≥ ν − κη; #pos ≥ n_U + 1 ⟹ λ_min(Id − W) < ν + κη, κ = (1 + ‖Y_s‖²/(σ(σ − η)))·max_B s_a2^{−2e_a}, σ = min_U δ_a2^{2e_a} ≤ λ_min(S_s).
- **Lemma F (exact-integer Cauchy sums).** k² − 4j² is an exact integer, so each term is two fixed-point truncated divisions with error < 1 unit; radii added.

## 3. Measurements (i): numerical ranks (rank.c; HODLR, points sorted by y, J″ Jacobi-scaled, ACA with exact generator updates)
Max rank of the off-diagonal blocks per level, at tolerance τ relative to max|J″|; storage = Σ r(|I| + |J|) + dense leaves (32), against the dense lower triangle n(n + 1)/2.
| c (even) | n | level 0 block / ranks at τ = 1e-16, 1e-32, 1e-64, 1e-100, margin | level 1 block / ranks | storage / dense at 1e-16 | at margin accuracy |
|---|---|---|---|---|---|
| 10 | 165 | 82 / 31, 51, 74, 81, 81 (1e-150) | 41 / 24, 37, 41, 41, 41 | 1.11 | 1.96 (no gain) |
| 16 | 413 | 206 / 37, 62, 95, 114, 114 (1e-150) | 103 / 30, 52, 82, 102, 103 | 0.65 | 1.52 |
| 25 | 1616 | 808 / 34, 62, 110, 152, 178 (1e-124) | 404 / 42, 75, 131, 185, 217 | 0.25 | 0.91 |
| 29 | 2289 | 1144 / 35, 65, 113, 159, 211 (1e-146) | 572 / 44, 82, 144, 205, 271 | 0.19 | 0.81 |
(c = 16 odd: as even within 2 ranks.) Ranks grow ≈ 1.1–1.4 per decade of τ (Zolotarev-type geometric decay of Cauchy blocks), so a binary64 H-matrix would be 5× smaller, but the certificate needs τ ≈ margin. Extrapolation to c = 47 (τ = 1e-244, n = 13,234): top levels ≈ 330, 410, 350, then full rank below block size ≈ 400: ≈ 3.0·10⁷ stored numbers against 8.8·10⁷, ≈ 7 GB at 1408 bits. **Numerical low rank (H-matrices) cannot bring c = 47 under 300 MB; the exact displacement structure needs 5n = 66,170 numbers.**

FMM for the Cauchy sums at c = 47 (371,000 sources, 12,365 targets): share of target-source pairs within a factor 2^R of each other: 23.7 % (R = 4), 34.9 % (6), 45.2 % (8), 54.6 % (10). A scale-separated multipole (p ≈ prec/R terms) removes the rest, so ≤ 4× at best; at 1000+ bits linear-distance FMM needs p ≈ 630 terms at separation 1/3 and does not pay. Lemma F instead: 0.42 μs/term at 1024 bits, 0.46 at 1408 (WeilZeta's Python: 2.6 μs).

## 4. Certificates (ii): two-sided enclosures of λ_min(Id − W)
Point mode (Lemma B) at ν = m(1 − 10⁻⁹) passes and at ν = m(1 + 10⁻⁹) fails (rigorously), m = the margin of WeilZeta/WeilDefl; so λ_min(Id − W) = m(1 ± 10⁻⁹) ± κη, and Q(f) ≥ ε‖f‖² on the sector. Same Chekmachev node sets (rigorous scan rscan.c, WeilZeta's code), same ε.
| c | parity | ε | bits | n = n_U + n_B | m (WeilZeta/WeilDefl) | κη | ball-mode bracket (bisection) | RSS: prep (Python) / LDL (C) | sums / LDL time |
|---|---|---|---|---|---|---|---|---|---|
| 10 | even | 9.133835606e-44 | 384 | 165 = 132 + 33 | 5.02635510181e-45 | 6.5e-99 | [5.02635510180509, 5.02635510180755]e-45 | 45 / 3.4 MB | 0.1 s / 0.0 s |
| 16 | even | 1.512800554e-75 | 512 | 413 = 280 + 133 | 7.29362304210e-77 | 1.5e-137 | [7.29362304210000, 7.29362304210356]e-77 | 65 / 3.6 MB | 1.8 / 0.2 s |
| 16 | odd | 2.105612047e-71 | 512 | 412 = 279 + 133 | 2.99318593786e-72 | 1.3e-132 | [2.99318593785708, 2.99318593785854]e-72 | 47 / 3.7 MB | 1.6 / 0.2 s |
| 20 | even | 2^-327 | 448 | 1100 = 789 + 311 | 3.10571228010e-97 | 1.5e-114 | [3.10571228010455, 3.10571228010607]e-97 | 66 / 4.3 MB | 11 / 1.3 s |
| 20 | odd | 2^-327 | 448 | 1099 = 788 + 311 | 6.70422678283e-92 | 2.9e-110 | fails at 448 bits; (6.704226, 6.704227)e-92 at 768 | 50 / 4.3 MB | 12 / 1.3 s |
| 25 | even | 2^-416 | 896 | 1616 = 1058 + 558 | 5.21043141861e-124 | 1.2e-246 | [5.21043141796, 5.21043141861]e-124 | 70 / 5.2 MB | 10 / 5.3 s |
| 25 | odd | 4.0546561810437569560e-119 | 896 | 1615 = 1057 + 558 | 8.07878305192e-120 | 4.0e-238 | [8.07878305091, 8.07878305192]e-120 | 57 / 5.2 MB | 10 / 4.6 s |
| 29 | even | 3.3887641037680338350e-145 | 1024 | 2289 = 1205 + 1084 | 1.42611291013e-146 | 2.7e-275 | [1.42611290942, 1.42611291013]e-146 | 82 / 5.8 MB | 44 / 10 s |
| 29 | odd | 1.5329466444386175069e-140 | 1024 | 2288 = 1204 + 1084 | 2.25792756276e-141 | 1.6e-275 | fails at 1024 bits; (2.2579275, 2.2579276)e-141 at 2048 | 75 / 6.2 MB | 46 / 9 s |
The set sizes (n_U, order, |G|) equal WeilZeta's and Chekmachev's at every c. WeilZeta's dense run: 0.06–1.08 GB up to c = 25 and did not fit 1.5 GB at c = 29; WeilDefl's deflated run: c = 29 in 0.27 GB and 850 s.

## 5. Ball mode and its barrier (interval wrapping)
Running the same recursion directly on Arb balls is also rigorous (each pivot ball excludes 0), but radii are re-amplified at every step. Digits lost beyond the inherent cancellation (precision − |log₁₀ pivot_min| − |log₁₀ worst relative radius|): 8 (c = 10), 11 / 37 (16 even / odd), 18 / 65 (20), 34 / 95 (25), 61 / ≈ 314 (29). Raising the working precision of slt alone does not help (input radii are amplified); only more precise inputs do. Without pivoting and generator balancing it fails already at c = 16 (step 141 of 412, odd). Point mode loses 17–33 digits at every c tested, so it is the method for c = 47.

## 6. Test for a new bound (iii): c = 47
Same pipeline (run47.sh): WeilZeta's full rigorous c = 47 scan (all 3.49·10¹⁰ nodes, ../WeilZeta/logs/scan_zeta_47_3683_rig.json), Chekmachev's node rules, ε = 2⁻⁸⁰³ = 1.8719·10⁻²⁴² (his ε₄₇ = 1.874·10⁻²⁴²; ours is 0.1 % smaller, so the statement is 0.1 % weaker).
| sector | bits | n = n_U + n_B, |G| | prep (Python) | Cauchy sums (3 chunks) | point LDL at ν | pivots | η, κ, κη | result |
|---|---|---|---|---|---|---|---|---|
| even | 1408 | 13,234 = 6064 + 7170, 371,000 (= Chekmachev's manifest: K+ 6063, K− 6301, bad even 869) | 372 s, 117 MB | 761 + 733 + 770 s, 101 MB, 0.46 μs/term | 575 s, 123 MB, ν = 3·10⁻²⁴⁴ | 6064 positive, 7170 negative, none zero | 4.8e-291, 3.5e13, 1.7e-277 | **λ_min(Id − W) ≥ 3·10⁻²⁴⁴ − 1.7·10⁻²⁷⁷ > 0** |
| odd | 1152 | 13,233 = 6063 + 7170, 371,000 (bad even 868) | 274 s, 95 MB | 556 + 482 + 621 s, 89 MB, 0.36 μs/term | 404 s, 106 MB, ν = 10⁻²⁵⁰ | 6063 positive, 7170 negative, none zero | 1.1e-291, 2.5e21, 2.7e-270 | **λ_min(Id − W) ≥ 10⁻²⁵⁰ − 2.7·10⁻²⁷⁰ > 0** (last pivot −9.0e-237: the odd margin is far larger) |
η ≈ 10⁻²⁹¹ at both precisions is a floor from double underflow in the bound (radii below 10⁻³⁰⁸ are rounded up), a valid overestimate; scaling the radii would lower it by ≈ 100 digits. A rank.c run on the c = 47 particles crashed (its rook search has no NaN guard); the c = 47 H-matrix numbers above stay an extrapolation.
**Theorem (computer-assisted; Chekmachev's minorant, our code).** For Weil's form Q of ζ (../WeilZeta/DERIVATION.md sec. 2), Q(f) ≥ 2⁻⁸⁰³‖f‖² for every f ∈ L²(−a, a) (both parities, hence all complex f), a = (ln 47)/2 = 1.925. This re-certifies the N1 record (Chekmachev, Zenodo 23848597) independently of his Rust and python-flint code: the scan is WeilZeta's, the matrix stage is ours, and every process stayed below 125 MB. It is not a new window. Beyond c = 47: all f at c = 53 needs WeilZeta's scan (≈ 10¹¹ nodes, ≈ 30 min with rscan2) and ≈ 10¹⁰ Cauchy terms (≈ 1.5 h at 0.46 μs) plus an O(n²) LDL (n ≈ 2·10⁴, ≈ 25 min), all within the memory cap; not run (time).

## 7. Adversarial test (iv)
explore27's moved pair (γ₁, γ₂ → 17.58 ± 0.3i, zero count kept; its exact Δμ(r) = Σ ±2 sin((γ − r)L)/(γ − r) added to μ in Arb, ζ's node sets) at c = 10, ε = 9.133835606e-44 and ε = 0: the inertia is rigorously wrong in both sectors (even: 2 positive B pivots at ν = 0, 1 at ν = −0.1, correct at ν = −1; odd: 1 at ν = 0 and −0.1, correct at −1), so λ_min(Id − W_moved) ∈ (−1, −0.1) and nothing is certified, consistent with explore27's float λ_min = −0.16 / −0.18. The linear algebra adds no false positives; all RH content stays in Chekmachev's minorant.

## 8. Memory and time
| c | dense (WeilZeta): RSS, time | deflated (WeilDefl): RSS, time | particle form (here): RSS (Python prep / C), time |
|---|---|---|---|
| 25 | 1.08 GB, 295–414 s | 0.15–0.17 GB, 284–311 s | 57–70 MB / 8.5 MB, 15 + 10 + 5 s |
| 29 | > 1.5 GB (stopped) | 0.26–0.27 GB, 843–850 s | 75–82 MB / 17 MB, 40 + 45 + 10 s |
| 47 | ≈ 35 GB, ≈ 1 day (estimate) | ≈ 1.3 GB, 10–12 h (estimate) | even 117 MB / 123 MB, 6 + 38 + 10 min; odd 95 MB / 106 MB, 5 + 28 + 7 min (measured, one core) |
The prep (Chekmachev's scan reuse, μ at all 371,000 j, node rules) is WeilZeta's code in a streaming mode; the matrix stage is O(n) in memory and O(n²) in time; the remaining cost at c = 47 is the 4.6·10⁹-term Cauchy sum.

## 9. Honest assessment
- New for the tree, not new mathematics: Lemma P (Chekmachev's matrix has displacement rank 2), the inertia form (Lemma I) and a backward-error certificate on it (Lemma B, Wiktorowicz's technique). The certified values are Chekmachev's; the gain is cost (O(n) memory) and independence of his code.
- Barriers measured: H-matrix compression at certificate accuracy (≤ 1.2× at c ≤ 29, ≈ 3× at c = 47); FMM for the Cauchy sums at 1000+ bits (≤ 4×); interval wrapping in ball mode (≈ 300 digits at c = 29 odd).
- Beyond c = 47 the wall is no longer memory: all f at c = 53 needs Chekmachev's scan (≈ 10¹¹ nodes, ≈ 30 min with rscan2) plus a Cauchy sum of ≈ 10¹⁰ terms (≈ 1–2 h) and an O(n²) LDL; even f at c > 71 needs ≥ 1.7·10¹³ scan nodes (weildefl's wall). The lattice minorant still cannot reach Wiktorowicz's c ≈ 100.
- Novelty: arXiv search (2 queries, approach arxiv) returned nothing relevant; prior art is Wiktorowicz–Wiktorowicz (displacement rank 2, structured point LDL, proved backward error; ../WiktAudit). The c = 47 result re-checks a published certificate; it is an audit result, not a discovery.
- Unchecked: the scan premise is WeilZeta's proved rscan (c ≤ 29) and WeilZeta's full rigorous c = 47 scan (../WeilZeta/logs/scan_zeta_47_3683_rig.json), reused here; Arb/FLINT 3.6 and python-flint 0.9 are trusted.

## Suggested NOGO line
- explore28 (Orca, exploration: fast multipole method and hierarchical matrices — particle representation, Schur/GKO structured elimination, H/HODLR compression, a-posteriori backward-error certification → N1) IDEA (engineering; the structure is Wiktorowicz's): Chekmachev's Woodbury–Cholesky matrix is the 1/(y − y′) interaction of particles at y = k², 4j² with two charges each (displacement rank 2), so the inertia of an O(n)-memory structured LDLᵀ with a backward-error bound decides Id − W > ν: two-sided λ_min enclosures equal WeilZeta/WeilDefl margins at c = 10–29 (both parities) in ≤ 82 MB, and c = 47 for all f is re-certified (ε = 2⁻⁸⁰³) below 125 MB per process (dense ≈ 35 GB); numerical low rank fails at certificate accuracy (HODLR storage 0.8–0.9 of dense at c = 25, 29, ≈ 7 GB at c = 47), FMM for the Cauchy sums gains ≤ 4× at 1000+ bits, ball-arithmetic elimination wraps (≈ 300 digits lost at c = 29 odd) while the backward-error bound loses 17–33; beyond 47 the walls are the scan and 10¹⁰-term Cauchy sums, not memory [ExploreFMM](Riemann/Deep/ExploreFMM/REPORT.md)

## Files
- prep_slt.py: Chekmachev's node sets, weights and Cauchy sources (WeilZeta's cert_zeta.py symbol and scan; LEAN=1 streams at c = 47); writes base.idx (particles) and base.src (sources).
- slt.c: Cauchy sums (Lemma F, chunkable: SLT_RANGE/SLT_FOUT/SLT_FIN), particle dump (-g/-l), structured LDL in ball mode (default) or point mode (SLT_POINT=1, Lemma B), bisection (-b lo hi steps).
- rank.c: HODLR ranks (sec. 3). val.sh, val_pt.sh, run47.sh: the exact run sequences. logs/: every JSON line and /usr/bin/time record.
- Replay (from this folder; flint via Homebrew; Python only via uv): cc -O2 -o slt slt.c -lflint -lgmp -lmpfr -lm -I/opt/homebrew/include -L/opt/homebrew/lib; cc -O2 -ffp-contract=off -o rscan ../WeilZeta/rscan.c; WZ_RSCAN=$PWD/rscan WZ_SCAN=/tmp uv run --no-project --with mpmath --with numpy --with python-flint python3 prep_slt.py 29 odd 1.5329466444386175069e-140 1024 /tmp/c29o; SLT_POINT=1 ./slt /tmp/c29o 1024 2.257927560502e-141 2.257927565018e-141 (≈ 2 min, < 90 MB: first line passes, second fails). c = 47: copy ../WeilZeta/logs/scan_zeta_47_3683_rig.json to scan/, then ./run47.sh even 2^-803 1408 3e-244 (≈ 55 min, < 125 MB) and ./run47.sh odd 2^-803 1152 1e-250 (≈ 40 min, < 110 MB); every process stays under its 1790 s alarm. The c = 47 runs used a byte-identical build for the sums and the current slt.c for the LDL.

RH relevance: N1.
