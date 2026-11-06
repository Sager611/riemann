RECORD: for three distinct primitive Dirichlet characters of one modulus and parity and any nonzero real b_j, at least 51087/2478195 = 2.061% of the zeros of b_1e^{iθ_1}L(s,χ_1) + b_2e^{iθ_2}L(s,χ_2) + b_3e^{iθ_3}L(s,χ_3) lie on the critical line with odd multiplicity, unconditionally (was 4.91·10⁻⁹, Dousselin 2025); for two functions 860123/2478195 = 34.707% (was 34.498%, ours; Bombieri–Hejhal's 34.5% needs GRH); for four or more functions the method gives nothing without a new lemma (stated below)

## Claim
Let q ≥ 3, χ_1, …, χ_N distinct primitive characters mod q with χ_j(−1) all equal, b_j nonzero reals, e^{2iθ_j} = i^a τ̄(χ_j)/√q, F = Σ_j b_j e^{iθ_j} L(s, χ_j) (functional equation, no Euler product). Then
 liminf_{T→∞} N_0^{odd}(F; T)/N_F(T) ≥ 1 − Σ_j min(1 − s₀, 1/N),  s₀ = 1669159/2478195 = 0.6735382,
where N_0^{odd} counts distinct zeros on the line of odd multiplicity. So N = 3: ≥ 3s₀ − 2 = 51087/2478195 = 0.0206146; N = 2: ≥ 2s₀ − 1 = 860123/2478195 = 0.3470764; N ≥ 4: 0.
Independent, weaker route on our earlier reviewed inputs only (stage 1): N = 3: ≥ 3p − 2 ≥ 0.0174738 with p = 0.672491294, the isolated-zero proportion of each L(s, χ_j) (IsoZero).
Three new pieces (DERIVATION §4, §8, §9): (1) Lamzouri's unconditional Bombieri–Hejhal proposition for N functions, with a capacity refinement (each L_j dominates on 1/N of the blocks);
(2) only isolation from the other zeros ON the line matters, because off-line zeros do not change the sign of the Hardy function; (3) Lemma M: merging κ-clusters of real points into multiple points
changes the Hilbert–Schmidt norm of the zero operator by O(κ√N), so the Lean-checked simple-zero certificate of arXiv 2610.08965 (67.3538%) also counts simple zeros κ-isolated from the other on-line zeros, for each fixed L(s, χ).

## Before
- N = 3: Dousselin (arXiv 2311.10285 = Forum Math. 37 (2025) 1871–1912, Selberg's method made explicit): κ_F ≥ 4.91·10⁻⁹ (N = 3), 3.68·10⁻⁹ (N = 4), 2.94·10⁻⁹ (N = 5), 2.161·10⁻⁶/(N log N) for N large.
  Correction to the brief and to R401: N ≥ 3 was not "only Selberg's non-explicit proportion"; Selberg (1998) proved c/N², Dousselin made it explicit. New/old for N = 3: 4.2·10⁶.
- N = 2: ours 0.3449825 (discoveries/linear-combinations-line.md); Lamzouri 2610.26841v2 Thm 1.4: 31%; Bombieri–Hejhal (Duke 1995, Thm C): 34.5% under GRH; 100% under GRH plus well-spacing.
- Lamzouri's Remark 1.5 needs δ_1 + δ_2 > 1; for N = 3 one needs isolated proportions above 2/3, and his 0.6552 gives 3·0.6552 − 2 < 0: N = 3 was out of reach of 2610.26841's constants.

## Method
1. Proposition N (DERIVATION §3–4): Lamzouri's proof (TeX l.909–964) line by line for N functions: joint Selberg CLT for N functions (BH Thm B; χ_iχ̄_k non-principal), separation set over all pairs
   (measure O(N²δT)), BH Lemma 9 per function, dominance of the argmax over N − 1 rivals, sign changes forced by isolated zeros of the dominant L_j. Count ≥ (1 − ε)(T/2π) log T − Σ_j (bad zeros of L_j in its own
   dominance blocks); capacity: those blocks carry (1/N + o(1)) of the zeros of L_j.
2. Theorem N^R (§8): the forced sign changes need only simple on-line zeros with no other on-line zero within 2πκ/log T.
3. Lemma M (§9): 2610.08965's proof of its Thm thm:line-general uses the zeros only through an abstract inequality (valid for every conjugation-invariant multiset), a certificate for every finite real set,
   the energy tr A² ≤ (C(f_ε) + o(1))N and spn ≤ N + o(N). Merge each κ-chain of real points of diameter ≤ √κ into its centroid (total multiplicity): A' − A involves real points only; centroids kill the
   second- and third-order terms of the mixed differences of K_ε², so ‖A' − A‖_HS ≤ c_ε κ √N (classical pair counts, Selberg's mean square of S(t+h, χ) − S(t, χ)); tr A'² ≤ (√tr A² + c_εκ√N)²;
   long chains are O(√κ N). Hence for every η > 0 some κ > 0 gives δ^R ≥ s₀ − η. For fixed χ the energy is Garunkštis–Paliulionyte's unconditional pair correlation over all zeros (as stated by Lamzouri, l.968–975),
   applied to R_ε = f_ε * f_ε and R_ε'' exactly as in 2610.08965's prop:energy.
4. Barrier (§7): from per-function counts, block counts and dominance measures alone the counting step is sharp (bad zeros of L_j may sit in L_j's own dominance blocks), so N ≥ 4 needs some δ_j > 3/4
   or a joint lemma.

## Missing lemma for N ≥ 4 (and to go beyond Σδ_j − (N − 1))
Non-concentration of the non-isolated zeros of L_j on the blocks where L_j dominates: limsup |H_j ∩ Ŝ_j|/|H_j| ≤ θ_j (Ŝ_j = blocks where |L_j| is largest, equivalently Selberg's polynomials P_j(γ) = Re Σ_{p≤X} χ_j(p)p^{−½−iγ} largest).
Then the proportion is ≥ 1 − Σ_j min(θ_j(1 − δ_j), 1/N); N = 4 needs θ < 0.766, N = 5: 0.613, N = 10: 0.306; independence (θ_j = 1/N) gives ≥ s₀ = 67.35% for every N.
Route: a pair-correlation theorem weighted by smooth functions of (P_1(γ), …, P_N(γ)) (Montgomery/BGST/Garunkštis–Paliulionytė twisted by short Dirichlet polynomials at the zeros); not found in the literature.
Numerically (q = 7, t ≤ 2400) the non-isolated zeros of L_j are under-represented in L_j's dominance blocks (0.18 vs 0.25 share), i.e. no adversarial concentration.

## Certificate
cd Riemann/Deep/LinComb3 && uv run --no-project --with python-flint --with mpmath python3 cert.py
(about 8 s, one thread; logs/cert.log). Stage 1 in Arb (imports IsoZero/cert.py hinge_p, read-only): p(10⁻⁵) ≥ 0.672491294, 3p − 2 ≥ 0.017473881 (asserted > 0.0174737), table N = 2…6, ceilings,
paper-input check 3p_A − 2 < −0.034. Stage 2 in exact rationals: s₀ = (2 − 6639737/5000000 − 4389·10⁻⁶)/(1 − 8722·10⁻⁶) = 1669159/2478195 (asserted equal to 2610.08965's value),
1 − 3(1 − s₀) = 51087/2478195, 1 − 2(1 − s₀) = 860123/2478195, N = 4: 0 (asserted).

## Checks
- Hypotheses of Lamzouri's proposition for N functions (DERIVATION §1–3): common gamma factor (same modulus, primitive, same parity) — the only use of "same modulus"; phases θ_j (e^{iθ_j}ε_j = e^{−iθ_j}); G real;
  N_F(2T) − N_F(T) = (T/2π) log T + O(T) (BH (3.12), Dousselin (eq:N)); CLT orthogonality including conjugate pairs χ, χ̄; dyadic blocks throughout; multiplicities (H_j, k_ℓ, N_F with multiplicity,
  zeros of F counted once). No GRH, no well-spacing: BH Thm B and Lemmas 7, 9 only as Lamzouri states them unconditionally (Selberg's density theorem).
- The N = 3 headline does not use the capacity step (1 − s₀ < 1/3); capacity only shows N ≥ 4 gives exactly 0.
- Reviewer (verdict-only, pasted summaries, logs/review.log): N-function extension SOUND; capacity step GAP (fixed-δ loss) → repaired (the factors 1 − √δ cancel; order ε → M → T);
  on-line isolation SOUND; centroid cancellation and HS triangle inequality SOUND; off-line terms cancel in A' − A, SOUND; pair counts from Selberg's mean square SOUND; application of 2610.08965's chain to the
  merged multiset: FLAW as first stated (b, ω ≥ 0; merged points must not coincide) → both hold (TeX l.1084; centroid in [min C, max C], ≥ κ from other points), added; fixed-χ energy: GAP (evenness of f) →
  η_ε is even (TeX l.759), added. A compressed whole-chain summary drew GAP ("off-line interactions in the HS norm"), answered by the narrow check above (off-line terms are identical in A and A').
- Numerics (q = 7, three odd characters, t ≤ 2400, three choices of b; numerics/, logs/analyse.log, logs/argcount.log): exact zero counts by the argument principle (integers to 10⁻⁴); zeros on the line with
  odd multiplicity 0.90, 0.70, 0.70 of all zeros (off-line 10–30%); dominance shares 0.33 each; the per-block inequality of Proposition N holds in all 304 strongly dominated blocks; nothing contradicts the bounds.
- DERIVATION constants rg-checked against cert.py output (p, 3p − 2, 2p − 1, 2/3 margin, 4 − 3C_MT, p_A, s₀, C'_0, x, b, 51087/2478195, 860123/2478195, θ thresholds, 2/√(2π)).

## Novelty
arXiv listings math.NT/CA/CV Oct–Nov 2026 (lst.py 'linear combination|Davenport.Heilbronn|critical line|Epstein|Selberg class|Bombieri.Hejhal'): only 2610.26841 (N = 2). arXiv API (all dates): Dousselin 2311.10285,
Lamzouri–Lee 2010.10490 (zeros near the line), Rezvyakova 2501.00551/2411.18492 (degree 2, Epstein: positive proportion, not explicit), Hejhal 1311.4862; 2609.15329 (Selberg's method optimised, ζ only).
Lamzouri 2610.26841: v2 is current. Zenodo '"linear combinations" AND "L-functions"': nothing. HAL: Dousselin's paper and thesis (tel-05317244, 2025), hal-05441218 (no proportion). GXL evidence API (2 claims):
nothing beyond Dousselin for N ≥ 3. Our files: only the N = 2 side result; R401's "non-explicit only" is out of date (Dousselin). Full log: logs/novelty.log.
Risk: nobody has run Selberg's method with the 2609.15329 optimisation on linear combinations; it could beat 4.91·10⁻⁹ by far, but it is not in the literature.

## Caveat
- Stage 2 rests on arXiv 2610.08965 (unrefereed; its local inequality is a 535-million-box computer proof, Lean native_decide; our audit Riemann/Deep/ZetaZeros found no flaw), on Lamzouri's dyadic form of
  Garunkštis–Paliulionytė (as in our N = 2 side result), on Selberg's mean square for S(t+h, χ) − S(t, χ) (quoted, not re-derived), and on BH Thm B and Lemmas 7, 9 as cited by Lamzouri (not re-read).
  Lemma M is new and checked only by verdict-only reviewer runs on pasted pieces; a specialist should read DERIVATION §9.
- Stage 1 (0.0174738) rests only on the inputs of the reviewed N = 2 side result (IsoZero) plus Proposition N.
- liminf only, fixed q, χ_j, b_j; no uniformity in q. The κ in Lemma M is not explicit; it need not be: Theorem N^R works with any fixed κ > 0, because the exceptional set in each Bombieri–Hejhal block has measure tending to 0.
- The barrier (§7.2) is about the method's inputs, not the zeros.

RH relevance: N12, N8 (critical-line proportions for functions with a functional equation and no Euler product: the Davenport–Heilbronn family; and, as a by-product of Lemma M, ≥ 67.35% of the zeros of
ζ and of each fixed L(s, χ) are simple, on the line and isolated from the other zeros on the line at some scale κ > 0)

## Files
- Riemann/Deep/LinComb3/REPORT.md (this), DERIVATION.md (every step and constant), cert.py (certificate), logs/cert.log, logs/novelty.log, logs/review.log
- numerics/zgrid.py (Hardy Z grids, writes Z_q7_n*.npy in the working directory), numerics/analyse.py, numerics/argcount.py; logs/zgrid.log, logs/analyse.log, logs/argcount.log
