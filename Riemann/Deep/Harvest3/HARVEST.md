CANDIDATES: 10 (best: arXiv 1803.02056 + 2605.26020, node N11, class groups of exponent 4 and 8 classified unconditionally: our explicit Littlewood bound empties the "at most one missing field")

Scope: arXiv monthly listings math.NT/CA/CV/SP for 2026-01 to 2026-06 (3966 entries; 2161 math.NT with cross-lists); 2569 abstracts read through the API (all math.NT plus title-keyword hits in CA/CV/SP); node-keyword scoring and a hand pass gave 169 node-touching papers, of which 27 are already ours (dropped unread) and 142 were read in TeX (`.scratch/harvest3/dl_scan.py`, Harvest 2's phrase list = the brief's list plus 14 similar phrases; 5 had no source). 106 had a match, 60 a strong one (187 strong lines, 716 "sufficiently large"/"absolute constant" lines; hits kept in `.scratch/harvest3/hits/`). A second abstract pass for GRH-, Siegel-zero- or Siegel–Walfisz-conditional statements found 2605.26020, 2605.02426 and 2601.21901, read in TeX with the cited source 1803.02056. Hexagon and Zenodo not scanned (Harvests 1 and 2). Each hit checked against TREE, COROLLARIES, DISCOVERIES, NOGO, REFUTATIONS, FORMALIZATIONS, discoveries/, Riemann/Tree and Riemann/Deep, and against Harvest 1 and Harvest 2. Harvest agent harvest3, 2026-11-06, 05:28–06:32 UTC.

Main finding: the January–June core-node papers (N1–N9) are mostly ours already (all 27 drops) or asymptotic. What is left comes from before October 2026: results whose only obstruction was a Siegel zero, GRH-sized least primes or the classical zero-free region. Our Littlewood78 bound, explprog's numeric Linnik bound and the Lean-checked half-plane remove these obstructions. The two best finds (rows 1 and 3) fell outside the node-keyword selection and came from the conditional-statement pass.

## Ranked candidates

| # | id (source line) | node | statement in the source | what is missing | inputs (explicit where?) | h | gain | risk |
|---|---|---|---|---|---|---|---|---|
| 1 | arXiv 1803.02056 Thm l.176–182, Lemma Siegel l.312–321, Lemma Chen l.338–352, §6 l.600–616 (Elsenhans–Klüners–Nicolae, Acta Arith.); arXiv 2605.26020 Thm l.258–270, Thm l.229–247, l.666–668 (Aiken–Ross–Shvydkoy–Xue, May 2026) | N11 | EKN: "Assuming the non-existence of Siegel-zeros, our computation found all fields with exponents 2,4,8. Without assuming any unproven statement … at most one missing field" (56, 203, 778 fields; largest \|D\| = 430,950,520). 2605.26020: the E = 2, 4, 8 tables over all discriminants and the classification of D whose CM points all lie on the boundary of the fundamental domain hold "conditionally on the non-existence of Siegel zeros". Its conductor step is unconditional (l.666–668) | the one possible exception: a field whose h(D) falls below Lemma Chen's √\|D\|/(7.732π log\|D\|) (EKN ran the search with that bound, l.611–613) | Littlewood78: h(−d) ≥ 0.0350√d/max(3, log log d) for every fundamental −d < −4 (Riemann/Deep/Littlewood78/REPORT.md l.1, l.35; Lean-checked half-plane). It beats Lemma Chen's 0.0412√d/log d for d ≥ 34 and Lemma Siegel's 0.0767√d/log d for d > e^{6.6} (checked here) | 2–3 | unconditional classification of imaginary quadratic fields with class group of exponent 4 or 8 (and 2); 2605.26020's E = 2, 4, 8 tables and its boundary-CM theorem become unconditional | low. OpenAI's QRH note (5 Oct) completed exponent 2 (idoneal numbers) only (Riemann/Deep/Combos/NOTE.md l.413). Check FSW 2610.18491, 2610.21511, Zenodo 23547020 and 23607492 |
| 2 | arXiv 2601.15564 l.138–142, l.150–153, l.156 (Jan 2026); 2602.22502 (RH-conditional analogue); R381 | N10 | Ingham: a prime lies in (n³, (n+1)³) for all sufficiently large n. "still an open problem as for whether [it] holds for all n ≥ 1" (l.142). The explicit route "appears futile": n ≥ exp(exp(33.3)) (Dudek 2016), exp(exp(32.537)) (Cully-Hugill thesis). Computation: n³ ≤ 1.649·10^40 | a single-exponential threshold. Every n is known only for k ≥ 8 (2610.20058, our ShortPrimes); k = 4…7 have no usable thresholds | Lean-checked half-plane B = 0.87338196661320752 (TREE N7); RH to 3·10^12 (Platt–Trudgian); Kadiri–Lumley–Ng 2018 Table 1, explicit N(σ,T) (checked in Riemann/Deep/RamRepair/NOTE.md §3.7); Hasanalizade–Shen–Wong N(T); Chirre–Helfgott 2512.15709 for n ≤ kT₀/(2π) | 6 | n₀(3) drops from exp(exp(32.537)) to about exp(600–800) (my model with a KLN-type C₁ ∈ [1, 10³]: the zeros near Re s = B contribute x^{−(1−B)/9}); first explicit n₀(k) for k = 4…7; stretch goal: every n for k = 7 if the density bound closes at x ≈ e^{202} | medium. R381 recorded "k = 7 would need a half-plane below 6/7: KNOWN", which holds only when the half-plane is used without density. Buthe78 REPORT l.89: KLN beats the trivial count only for σ ≳ 0.85 at T = 3·10^12, so the k = 7 stretch probably fails. Visser, Johnston or Cully-Hugill may get there first (no Oct–Nov arXiv title on 6 Nov) |
| 3 | arXiv 2605.02426 Thm main3 l.127–128 (GRH), abstract l.103, l.118–124, Rem. l.252–262, l.329 (in a comment block) (May 2026) | N11 (N10) | every n > 24 is a prime plus a non-squarefree integer: proved for odd n, and for all n under GRH. Unconditionally it holds only for large n, via the ineffective Siegel–Walfisz theorem. "cannot be achieved … unless an explicit version of the Siegel–Walfisz theorem becomes available" | even n > 8·10^9 | explprog: least prime ≡ a (mod m) ≤ (0.72 m log m)^{7.8978} for m ≥ 10⁴ (TREE N11; Riemann/Deep/ExplProg); the paper's check for 24 < n ≤ 8·10^9 (l.143). My check: for every prime q ≤ 200, each reduced class mod q² contains a prime ≤ 6,734,249 (worst q = 197), and (0.72 q² log q²)^{7.8978} < ∏_{p<q} p for 200 < q < 3000 (`.scratch/harvest3/check_nsf.py`) | 1 | the paper's open problem solved unconditionally. Take q, the least prime not dividing n, and a prime p ≡ n (mod q²) with p < n | low. Routine (COROLLARIES). Check v2 and Zenodo |
| 4 | arXiv 2605.03349 Thm l.41–43, l.49–51; Lemma siegel l.82–85, used at l.590–591 | N11 | for primes q ≥ q₀(ε), every reduced class mod q has n, m ≤ q^{5/2+ε} with λ(n) = −1, λ(m) = +1; "the threshold … is uncomputable" (Siegel) | a numeric q₀(ε) | Littlewood78 L(1,χ_q) ≥ 0.11/max(3, log log q) in place of Siegel; explicit large sieve (Montgomery–Vaughan); explicit Pólya–Vinogradov for the O(√q) of Lemma hyperbola l.97–103 | 5 | first numeric threshold (8-page proof; every ≪_ε to be tracked) | low–medium; modest significance |
| 5 | arXiv 2603.16023 eqs l.547–549, Rem. l.611, Discussion l.1753–1757 | N9 | explicit \|Δ_k(x)\| bounds for 3 ≤ k ≤ 9 from 0.611/0.566 t^{1/6} log t and explicit 2nd/4th moments; "an improvement in any of these bounds … would directly lead to slight improvement in the exponent of x" | our sharper ζ inputs | Explicit611 (0.568 t^{1/6} log t for t ≥ 3, 0.4821 for t ≥ 8.97·10^17); SubWeyl 54.6 t^{27/164} (t ≥ 3) | 3 | explicit Δ_k with the sub-Weyl exponent and smaller constants; a small gain | low; routine |
| 6 | 2605.21876 l.318, 429, 701, 766, 784; 2604.11941 l.2719; 2603.11196 l.209, 1114, 6149 | N11 | constants ineffective only through Siegel or Siegel–Walfisz (reverse Goldbach, unconditional simultaneous non-vanishing, shifted-prime laws) | computable constants | COROLLARIES "Effective Siegel–Walfisz" (2610.21511 gap); explprog | 1–2 each | effective, not numeric | low gain |
| 7 | arXiv 2603.17334 Thm 2 l.149–155, l.163, Rem. l.1066–1068 | N8 | unconditional Θ = ϑ = 2.36·10⁻⁶ for r-gaps between zeros; "It is possible to sharpen these bounds with more work" (Korolev's inequality refined by their Prop. on moments of S(t+h) − S(t)) | the refinement | the paper's §§3–5 (explicit Tsang); STT22 | 4 | larger explicit Θ; numeric only | low gain |
| 8 | arXiv 2603.01711 Thm LBLD l.101–112, l.181 | N9 | K_α ≫ (Cα² log α)^{−α²} in the lower bound for large deviations of log\|ζ(1/2+it)\|; "We make no effort to specify the absolute constant C" | numeric C | the Arguin–Bailey off-line machinery (AB24), mean values of Dirichlet polynomials | 6+ | first explicit C; weak RH link | high |
| 9 | arXiv 2606.22851 l.6225, l.7139 | off-tree (N10, computation) | δ in D(x, log^A x) > 1 + δ is "in principle explicitly computable"; C₀ "effectively computable" | numbers | the paper's own lemmas | 4 | first explicit δ, C₀ | low relevance |
| 10 | arXiv 2601.21642 l.274–275 | off-tree (N12) | in O_k(q^{−η}), the dependence on k is not explicit (Fouvry–Kowalski–Michel–Sawin complexity bounds) | explicit k-dependence | the FKMS conductor/complexity bounds | 6+ | explicit k-dependence | high; off the RH path |

## Mini-briefs (top 5)

**1. Class groups of exponent 4 and 8, unconditionally (N11).** Goal: show that EKN's lists (56, 203 and 778 fundamental discriminants with exponent 2, 4, 8) are complete, that 2605.26020's E = 2, 4, 8 tables over all discriminants are complete, and that its boundary-CM theorem (l.229–247) holds.
Inputs: Littlewood78 (REPORT.md l.1, l.35); EKN 1803.02056 (Lemma Siegel l.312–321: 0.655/(πe) = 0.0767; Lemma Chen l.338–352: 1/(7.732π) = 0.0412; Thm upp_bound l.401–413: N₈ = 58, d₅₈ ≤ 3.17·10^110; §6 l.600–616: the search used Lemma Chen); 2605.26020 l.258–270, l.666–668.
Method: Arb proof that 0.0350/max(3, log log d) > 0.0412/log d for all d ≥ 10⁶ (case split at d = e^{e³}), so no field violates Lemma Chen and EKN's single exception cannot exist; then 2605.26020's unconditional conductor algorithm. Optional: rerun EKN's Rédei-matrix search for E = 4 (k ≤ 24) or 2605.26020's code.
Certificate: Arb inequality script; if a search is rerun, the counts 203/778 and the largest \|D\| = 430,950,520.
Verify: reviewer checks that EKN use the class-number bound only to cap \|D\| for each number of prime factors (no other Siegel input), and that Littlewood78 covers every fundamental D (D ≡ 0 mod 4 and 8 included). Novelty: OpenAI QRH note (5 Oct), FSW 2610.18491, 2610.21511, Zenodo 23547020/23607492, Oct–Nov listings 'exponent|idoneal|class group' (none on 6 Nov), Zenodo '"class group" AND exponent'.

**2. Explicit Ingham for cubes and k-th powers (N10).** Goal: x₀(k) such that every (n^k, (n+1)^k) with n^k ≥ x₀(k) contains a prime, for k = 3…7 (target n₀(3) ≤ e^{1000}). Test whether n₀(7) ≤ 7T₀/(2π) = 3.34·10^12; if so, every n works for k = 7.
Inputs: Bebel's B (TREE N7); T₀ = 3·10^12; KLN 2018 Table 1 (RamRepair NOTE §3.7); HSW N(T); Chirre–Helfgott 2512.15709 (low range); Dudek 2016 (Funct. Approx.) and the Cully-Hugill thesis as the template for the smoothed explicit formula; 2601.15564 Thm compthm.
Method: a smoothed explicit formula for θ on (x, x+h] with h = k x^{1−1/k}. Zeros below T₀ lie on the line. For zeros above T₀, bound ∫ x^{σ−1} dN(σ,T): use KLN for σ ≥ σ₁ and the trivial count below (σ₁ < 1 − 1/k). Optimise the kernel order and T; prove monotonicity in x once.
Certificate: Arb (python-flint) on a log-x grid. Sanity check: Dudek's double-exponential threshold must come back when the classical region replaces B.
Verify: reviewer re-derives the truncation, ψ-versus-θ and trivial-zero terms. Adversarial focus: the σ ≈ 0.85 band at T ≈ 10^13. Novelty: Visser 2610.20058, Lee 2602.14340 (R381), the Cully-Hugill thesis, listings 'cube|consecutive powers|Ingham|Legendre', Zenodo '"consecutive cubes"'.

**3. A prime plus a non-squarefree integer, for every n > 24 (N11).** Goal: unconditional proof of 2605.02426 Thm main3.
Inputs: explprog's numeric Linnik bound (TREE N11; Riemann/Deep/ExplProg); the paper's computation for 24 < n ≤ 8·10^9 (l.143); Rosser–Schoenfeld lower bounds for θ(x).
Method: let q be the least prime not dividing n; a prime p ≡ n (mod q²) with p < n gives q² \| n − p. For q ≤ 200, least primes in all reduced classes mod q² are ≤ 6,734,249 (`check_nsf.py`). For q > 200, n ≥ ∏_{p<q} p > (0.72 q² log q²)^{7.8978} (checked for q < 3000; beyond via θ(x) > x(1 − 1/log x)).
Certificate: exact sieve to 10⁷, rerun; Arb for the q > 200 inequality.
Verify: reviewer confirms explprog covers every reduced class and every modulus, and reruns the n ≤ 8·10^9 check (or an independent check to 10⁹). Novelty: the paper's v2, Zenodo 'square-free AND prime', Oct–Nov listings 'square-?free' (none on 6 Nov).

**4. Liouville sign changes with a numeric threshold (N11).** Goal: a numeric q₀(ε) in 2605.03349 Thm main, for ε ∈ {1, 1/2, 1/4, 1/10}.
Inputs: Littlewood78 (in place of Lemma siegel, l.82–85, used at l.590–591); explicit large sieve (Montgomery–Vaughan); explicit Pólya–Vinogradov; the paper's Lemmas comb, periodicity, dominant, exceptional (l.116–470).
Method: track each ≪_ε constant. The end step needs q^{1−ε/9}·C < Ŵ(0)·0.11 q/max(3, log log q) − C′√q.
Certificate: Arb table of q₀(ε); brute force over primes q ≤ 2000 (λ in every class up to q^{5/2+ε}).
Verify: reviewer checks every lemma's constant. Novelty: the paper's v2 and Matomäki–Teräväinen (products in progressions).

**5. Explicit Δ_k from our ζ bounds (N9).** Goal: rerun 2603.16023's explicit Δ_k(x), 3 ≤ k ≤ 9, with 0.568 (and 0.4821) t^{1/6} log t and with 54.6 t^{27/164}, re-choosing T for the sub-Weyl exponent.
Inputs: 2603.16023 l.547–549, l.611, l.1753–1757; Riemann/Deep/Explicit611 and the subweyl REPORT.
Method: replace the pointwise bounds, re-optimise T, ε and ε₁, and find where the sub-Weyl exponent wins.
Certificate: Arb tables like the paper's Tables (fourth- and second-moment cases).
Verify: reviewer (routine). Novelty: Tudzi 2025 (k = 5, 6), Bellotti–Yang 2024 (k ≥ 30), the paper's v2.

## Dropped hits
- Already ours (27, not read): 2602.05614, 2605.03665, 2602.03626, 2605.27833, 2606.16097, 2602.06199, 2601.22413, 2602.04022, 2602.14340, 2603.21490, 2606.12376, 2605.20224, 2605.30767, 2606.22562, 2604.24505, 2606.27516, 2602.20313, 2604.25579, 2606.09096, 2602.15905, 2606.27323, 2601.12495, 2603.03029, 2604.00693, 2606.31441, 2602.05731, 2606.17494. Read but cited in our Deep files: 2606.30567 (RCApps2, CM3 NOVELTY; ineffective via Siegel), 2601.04189 (RSFam), 2605.09282 (N8Fam; no source), 2605.31356 (ExploreStatMech), 2604.17532 (RSFam, ChebST).
- Conditional on GRH, GLH, RH or other hypotheses: 2603.25612 (GLH), 2601.21901 (GLH for m > 2 via Burgess-type sums, l.381–412; ERH Kummer lemma l.756–772), 2602.22502 (RH), 2603.20093 (GRH + LI), 2604.05733 (RH), 2601.15610 (RH), 2606.10242, 2606.25094, 2603.21256, 2606.01514 (GRH), 2605.24363 (θ = ∞ hypothesis), 2603.28104 (narrow-box hypothesis), 2604.11941 main theorem (GRH), 2605.17825 (GRH part; unconditional K = 8 needs a research project), 2605.19566 (U = N^{4/49}: the worst σ ≈ 3/4 lies below B, so the half-plane removes only the exp(log^{2/3+ε} N) factor).
- Already explicit, or limited by inputs other than zeros: 2603.17334 (row 7 only), 2604.19681, 2606.07785, 2606.04677 (M < 1.78·10^32 is bound by the digit-discrepancy rate c₄ = log 10/6, Rem. l.155–163, not by zeros), 2603.10356 and 2602.22720 (sieve; squares need √x intervals), 2601.15564 Thm almostthm (P₂ between all cubes, done).
- Off-tree or research-level: 2604.06848 (optimal exponent, l.257), 2603.12159 (monotonicity of C_{2d+1}), 2603.23336 (Bohr's problem), 2602.07690 (Beurling), 2602.07558, 2606.12907, 2604.23661, 2606.02407 (Barnes analogue of Lindelöf, l.1270), 2601.18687, 2606.29555, 2606.07312 (density of log ζ on the line, open), 2601.18408 (qualitative Brauer–Siegel descent), 2603.24548 (a limit relation), 2605.00695.

## Flags for our files
- Riemann/Tree/R381.lean l.12–15: "k = 7 would need a half-plane below 6/7 < B_2: KNOWN" holds only for the half-plane alone. With an explicit zero-density estimate the k = 3…7 thresholds become single-exponential, and every n for k = 7 is still open for our tools (row 2).
- Littlewood78's page (discoveries/littlewood-half-plane.md) can cite row 1 as an application once it is done. Harvest 2 already flagged its "before" line (Zenodo 23400116, 23459046).
- 2603.16023 uses 2602.05614's 0.611/0.566 (its "R2026"). Our Explicit611 values improve them, and nothing in that paper is wrong.

## Papers scanned
2601.00292 nosrc · 2601.00660 hit · 2601.00910 hit · 2601.02623 hit · 2601.03984 hit · 2601.04189 hit
2601.06292 hit · 2601.09276 none · 2601.11052 none · 2601.12601 hit · 2601.13855 none · 2601.15448 none
2601.15564 hit · 2601.15610 hit · 2601.16193 hit · 2601.17079 hit · 2601.18025 none · 2601.18408 hit
2601.18687 hit · 2601.19177 hit · 2601.19429 hit · 2601.21642 hit · 2601.21901 hit · 2602.01409 none
2602.05788 hit · 2602.06304 hit · 2602.07558 hit · 2602.07690 hit · 2602.08985 hit · 2602.09797 hit
2602.10497 none · 2602.12371 none · 2602.13058 none · 2602.19132 none · 2602.20917 hit · 2602.22127 hit
2602.22167 none · 2602.22502 hit · 2602.22720 hit · 2602.23100 hit · 2603.00768 nosrc · 2603.01711 hit
2603.02118 hit · 2603.05609 hit · 2603.07641 none · 2603.10241 hit · 2603.10356 hit · 2603.10746 hit
2603.11196 hit · 2603.12159 hit · 2603.15197 none · 2603.16023 hit · 2603.17334 hit · 2603.18802 hit
2603.20093 hit · 2603.21005 hit · 2603.21256 hit · 2603.21506 hit · 2603.21535 none · 2603.21739 none
2603.22124 hit · 2603.23250 hit · 2603.23336 hit · 2603.24313 none · 2603.24548 hit · 2603.25076 none
2603.25612 hit · 2603.25814 nosrc · 2603.25915 hit · 2603.25961 hit · 2603.26507 hit · 2603.28104 none
2604.02306 hit · 2604.02614 hit · 2604.02960 hit · 2604.03051 none · 2604.05733 hit · 2604.06848 hit
2604.07474 hit · 2604.08776 none · 2604.11205 hit · 2604.11941 hit · 2604.17532 hit · 2604.17832 none
2604.19681 hit · 2604.21195 hit · 2604.22058 none · 2604.22272 none · 2604.23517 hit · 2604.23661 hit
2604.23900 hit · 2604.24481 none · 2604.24754 hit · 2604.25969 hit · 2604.26564 hit · 2605.00695 hit
2605.02426 hit · 2605.03349 hit · 2605.05180 none · 2605.07869 hit · 2605.09282 nosrc · 2605.09807 hit
2605.13715 hit · 2605.14871 hit · 2605.17012 hit · 2605.17825 hit · 2605.18676 hit · 2605.19105 hit
2605.19566 hit · 2605.21476 none · 2605.21876 hit · 2605.22421 none · 2605.23014 hit · 2605.24363 hit
2605.26020 hit · 2605.29470 nosrc · 2605.29991 hit · 2605.31356 hit · 2606.01115 none · 2606.01384 hit
2606.01514 hit · 2606.02407 hit · 2606.04512 none · 2606.04677 hit · 2606.04744 hit · 2606.05224 hit
2606.07312 hit · 2606.07329 hit · 2606.07785 hit · 2606.09818 hit · 2606.10242 hit · 2606.11813 none
2606.12484 none · 2606.12907 hit · 2606.15635 hit · 2606.22851 hit · 2606.25094 hit · 2606.25849 hit
2606.25903 hit · 2606.26440 hit · 2606.28657 hit · 2606.29555 hit · 2606.30567 hit · 2606.30603 hit
2606.30618 hit
Also read as an input (outside the window): 1803.02056 (EKN, 2018).
("hit" = at least one phrase match, strong or weak; "nosrc" = no TeX source (HTTP 404 from arxiv.org and export.arxiv.org).)

RH relevance: N11, N10, N9, N8
