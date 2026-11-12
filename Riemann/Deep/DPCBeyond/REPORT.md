THEOREM: unconditionally D(λ) > 0 for every λ > 7/11 = 0.6364 (was 2/3): Durkan–Pearce-Crump's two-line method runs at every mollifier length θ < κ/2 = 4/7, because the fourth moment can be replaced by a pointwise bound of R_T through the zeros (Lemma Z) plus Conrey's density slope κ = 8/7.

# REPORT — dpcbeyond (Orca, N8): beyond 2/3 in the unconditional two-line method

## Short answer
- Yes for θ ∈ (½, 4/7): an unconditional substitute for the fourth moment exists, and it is not a mean value theorem. **Lemma Z** (PROOF.md): for t ∈ [T/2, 3T], 0 ≤ y ≤ ½,
  |R_T(t − iy)| ≤ C(1 + W_{−y})(1 + On)^{k−1}, |R*(t − iy)| ≤ C X^{y}(1 + W)(1 + On)^{k−1}, X = T^θ,
  On = Σ_ρ κ_A(L|γ − t|) (smoothed local zero count), W = Σ_ρ X^{(β−½)⁺} κ_A(L|γ − t|). Proof: multivariate explicit formula for the prime factors of n, with the cutoff f(Σ log q_i/L) kept inside the r-fold integral, so Σu_i ≤ log X and an r-fold product costs X^{max d}, not X^{r·max d}. DPC's RH route (l.665–670, |e_r| ≪ (1 + |S_1|)^r) loses exactly this.
- Consequences (Lemma M): ∫_T^{2T}|R_T|^q ≪ T for every q < κ/θ (candidate (a), with δ = κ/θ − 2 > 0 for θ < 4/7); uniform integrability (candidate (b)) follows; Step 7 needs only q > 2 (Step 7′). Hairpins (candidate (c), Step 6e′): bounded pointwise by Lemma Z instead of Hölder with ‖R‖₄, cost e^{2θa} per depth bin against e^{−κa}: threshold 2θ < κ.
- Hence Theorem U′: every φ > max{2/(2 + κ), ½} works (θ < min{κ/2, 1}; the ½ and 1 caps bind only if κ > 2); κ = 8/7 − δ gives 7/11. With Jutila's κ = 1 − δ it gives back 2/3: the whole gain is Conrey's mollifier length 4/7 > ½, imported through the density theorem.

## The three candidates
| Candidate | Status | Where |
|---|---|---|
| (a) ∫\|R_T\|^{2+δ} ≪ T, θ > ½ | PROVED for θ < κ/2 with 2 + δ < κ/θ (Lemma Z + density + Fujii) | PROOF.md Lemma M (M3) |
| (b) uniform integrability of \|R_T\|² | follows from (a); Step 7 needs only (a) with any q > 2 | Step 7′ |
| (c) new allocation in 6e | pointwise Lemma Z on the slits; Hölder only over zeros, never on R_T | Step 6e′ |
| fourth moment ∫\|R_T\|⁴ ≪ T, θ > ½ | not proved and not needed; Lemma Z gives it only for θ < κ/4 = 2/7 (MV: ½) | — |

## Obstruction side (why 7/11 is this method's end with Conrey's κ)
- Heuristic, not proved (referee C): the estimates show that θ < κ/2 suffices, not that it is necessary. Heuristically a zero of depth d = a/L makes R_T ≈ X^{d}·(polynomial loss) near γ (the mollifier residue; N1 shows I_j ≈ e^{θa}, while the lossy route gives e^{jθa}). If so, a zero configuration consistent with Conrey's theorem (isolated off-line zeros, N(½ + a/L) ≍ TL e^{−κa}) would give ∫|R_T|^q ≍ T Σ_a e^{(qθ − κ)a}/poly(a), so density information alone would not yield (a) with q ≥ κ/θ, and (a) with q > 2 would need θ < κ/2.
- The hairpin of such a zero is ≈ ∫_0^{d} X^{d−y}X^{d+y} dy = d X^{2d} (R_T damped, R* amplified at the mirror point): Σ_a TL e^{−κa}(a/L)e^{2θa} needs 2θ < κ. Heuristically these are terms of the exact identity, not artefacts of Hölder, so a sign/cancellation between hairpins or a better uniform slope κ near σ = ½ would be needed to move the limit λ₀(κ) = max{2/(2 + κ), ½} (κ = 2, i.e. the density hypothesis near ½, would give DPC's RH value ½).
- Mean value routes for θ > ½ (what the brief proposed first) are blocked independently: R_T² has length T^{2θ} > T, and ∫|R_T²|² ≪ T is equivalent to a Selberg-type variance of almost-primes in intervals (x, x + x/T], x ≤ T^{2θ}, with square-root cancellation — out of reach unconditionally; Huxley's large values cover only V ≥ T^{θ/4}, MV's fourth moment only V ≥ T^{θ−½}; the range 1 ≪ V ≤ T^{θ−½} is exactly where Lemma Z (zeros) is needed.


## Numerics (logs/, all under memcap 290 MB, single thread; peaks 0–146 MB)
- N1 n1_mechanism.py: planted zero, depth a/L, j copies at t = γ (θ = 0.55, η = 0.1): I_j/e^{θa} decays polynomially (a = 40: I_2 = 2.5·10⁶, I_1² = 3.5·10¹², e^{θa} = 3.6·10⁹). The cutoff keeps X^{max d}.
- N2 n2_tails.py: DPC's R_T at T = 10⁴, 10⁶, θ = 0.45, 0.55, 0.65, 2·10⁵ random t: E|R|²/Σ|c|² = 1.000, E|R|⁴/(E|R|²)² = 1.4–1.8, max |R|²/Σ|c|² ≤ 10.8, L²-tail mass above 8Σ|c|² ≤ 0.0008. No concentration at θ > ½ (at these heights all zeros are on the line, the on-line regime of Lemma Z).
- N3 n3_identity.py: the j = 1, 2 multivariate explicit formula behind Lemma Z with 482 zeta zeros (mpmath) near t = 1000: |LHS − RHS| ≤ 4.5·10⁻¹² (j = 1), 9.8·10⁻¹³ (j = 2).

## Review
- Review 1 (logs/review1.txt): Lemma Z GAP (block prime q vs packet q^h at s = ½ − y + it), repaired; Lemma M, Step 6e′, Step 7′ CORRECT.
- Review 2 (logs/review2.txt): Lemma Z GAP (prime-power packets carry ρ_η(h log q/L)/h), repaired as proposed, 'no further obstruction in (i)–(iv)'; adversarial reviewer: NO FLAW FOUND (mechanism sound; no contradiction with known theorems — the direction is density ⇒ moments; asked to confirm that the fourth moment was Theorem U's only θ < ½ input: PROOF.md for Theorem U says so at Step 5a, and a text search finds θ < ½ only in Steps 0, 5a, 6e–6f, 7 and in Step 8's discussion).
- Review 3 (logs/review3.txt): Step 6e′ including the all-zero discrete moments CORRECT (wording repairs applied); audit of Theorem U's Steps 0, 3, 7 at θ ∈ [½, 4/7): NO OTHER θ<1/2 USE.
- Review 4 (logs/review4.txt): third pass on the repaired Lemma Z: CORRECT (notation fix applied).
- Independent referees (orchestrator, 12 Nov 2026): A (Lemma Z) CORRECT, SIGNIFICANT; B (Lemma M, 6e′) GAP repairable, SIGNIFICANT — M4 radius (m+1)/L → (|m|+1)/L, repaired; C (Step 7′, order of choices, assembly) GAP repairable, SIGNIFICANT — range now φ ∈ (max{2/(2+κ), ½}, 1), θ < min{κ/2, 1} (unchanged for κ = 8/7), and the sharpness claim weakened to a heuristic; both repaired in PROOF.md (Status lists them); CLOSED re-reviews B′ and C′: both CLOSED (orchestrator, 12 Nov 2026).
- Not done: line-by-line referee of the whole chain at θ ∈ (½, 4/7); Lean (excluded by the brief, mandatory for a discovery).

## Novelty (logs/novelty_axq.txt)
- arXiv API (8 queries: mollifier moments unconditional; explicit formula + mollifier + zero density; small gaps; uniform integrability of Dirichlet polynomials; Dirichlet polynomial moments + zero density; positive proportion of gaps; Durkan/Pearce-Crump; large values + mollifier) and the math.NT listings 2026-10 (854 entries) and 2026-11 (280 entries), titles filtered for gaps/zeros/mollifier/density/large values: no unconditional small-gap proportion below 2/3 other than ours (Theorem U), no zero-based pointwise or moment bound for mollifiers longer than T^{1/2}. Closest: DPC 2610.28527 l.643–670 (the explicit formula for one prime sum under RH, then |e_r| ≪ (1 + |S_1|)^r), 2610.18923 (large values for θ ∈ (½, 4/7), different problem), 2611.05000 (long Dirichlet polynomial approximations to (ζ'/ζ)^k, under RH). Zenodo/HAL not searched this round.
- Lever: λ₀(κ) = max{2/(2 + κ), ½} for any uniform density slope κ near σ = ½ with only a log T loss; the recent 67% critical-line papers (2610.08965, Lamzouri 2609.02882, Alpöge–Furman 2608.13637) were not checked for a better κ.

## Files
Riemann/Deep/DPCBeyond/: REPORT.md, PROOF.md (Theorem U′: Lemma Z, Lemma M, Steps 6e′, 7′), DERIVATION.md (D1–D6), n1_mechanism.py, n2_tails.py, n3_identity.py, logs/ (n1, n2, n3, review1, review2, novelty_axq).

## Suggested NOGO line
- dpcbeyond (Orca, N8 small gaps, 12 Nov 2026): the fourth-moment cap θ < ½ of the unconditional two-line method is removable (Lemma Z: R_T bounded through the zeros, cost X^{max d} not X^{kd}); the new wall is θ < min{κ/2, 1}, λ > max{2/(2+κ), ½} = 7/11 with Conrey: heuristically one isolated off-line zero of depth d costs X^{2d} in its hairpin and X^{qd} in the q-th moment (sufficiency proved, necessity not); past 7/11 needs a uniform density slope κ > 8/7 near σ = ½ or cancellation between hairpins.

RH relevance: N8

