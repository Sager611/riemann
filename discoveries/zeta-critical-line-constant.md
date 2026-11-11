# Explicit bound for ζ on the critical line: |ζ(1/2+it)| ≤ 0.568 t^{1/6} log t (node N9)

**Claim.** For all t ≥ 3, |ζ(1/2+it)| ≤ 0.568 t^{1/6} log t. For all t ≥ 8.97·10^17, |ζ(1/2+it)| ≤ 0.40913 t^{1/6} log t + 3.01756 t^{1/6} − 2.920 ≤ 0.4821 t^{1/6} log t. For t ≥ 10^10 the constant is 0.5598. Taking the lemmas of 2602.05614 exactly as stated (one lemma's t ≥ 5·10^6 kept), the first bound holds with 0.5782.

**Before.** 0.611 for t ≥ 3 and 0.566 for t ≥ 8.97·10^17 (Revers, arXiv 2602.05614 Thm 1, Feb 2026, unrefereed); 0.618 for t ≥ 3 (Hiary–Patel–Yang, J. Number Theory 256 (2024), arXiv 2207.02366). Nothing better on arXiv, Zenodo or HAL as of 6 Nov 2026. The method cannot go below 0.5075 for t ≥ 3 (2602.05614 l. 96–101).

**Method.** Same lemmas and method as 2602.05614: a Riemann–Siegel main sum, a trivial part 2Σ_{n<4K} n^{−1/2}, and blocks of length K bounded by van der Corput / Weyl differencing (Lemmas 2.3–2.9). The paper's uniform relaxations are removed:
- Each block gets its own constant μ_r and its own differencing parameter c_r.
- The last, partial block is bounded by its true length (the paper's own remark, l. 2045–2049).
- Each of 39,394 intervals t ∈ [2πn₁², 2π(n₁+1)²) up to 10^10 gets its own integer K ∈ [t^{1/3}, t^{0.35}], instead of 10 hand-tuned regions.
- The Weyl sums are exact for small n₁ and for the first tail blocks.
- The tail split √(a+b) ≤ √a + √b is replaced by a monotone sum–integral comparison.
- ζ on [3, 300] is evaluated in Arb, not quoted from Hiary.

0.611 binds where the trivial and van der Corput bounds meet (t ≈ 3.3·10^7). 0.568 binds at t ≈ 2.30·10^6 (n₁ = 605, K = 133).

**Certificate.** certify.py (python-flint 0.9.0 Arb, 96 bits, about 6.5 min on one thread, 49 MB): parts A (balls on [3, 300]), B (exact trivial bound, n₁ ≤ 500), C (39,394 intervals to 10^10), E (tail t ≥ 10^10), F (t ≥ 8.97·10^17). Every lemma hypothesis is asserted on every interval (K³ > t, K − 1 < t^{0.35}, R ≥ 4, t ≥ 2012, 0 < c < 3, M − 1 < 3W^{1/3}). The orchestrator reran it on a copy (6 Nov 2026): "ALL PARTS OK", output identical to run_final.log. A script reproduces all 40 region coefficients of the paper to every printed digit.

**Review.** Two reviewers read the paper's TeX at the cited lines.
- Items 6–9 (exact block sums, tail integral, the t ≥ 5·10^6 hypothesis, the low range): CORRECT. The proof of Lemma 2.7 uses t ≥ 5·10^6 only to get R ≥ 4. The asinh bounds must point the right way; the code uses the exact asinh at R and the upper bound at r₁.
- Items 1–5 (per-block constants, last block, per-interval K, exact Weyl sums): each valid. One "gap" asked whether K/t^{1/3} ≤ 1 + t₀^{−1/3} holds for the freely chosen K. The orchestrator checked the code: the tail parts E/F use K = ⌈t^{1/3}⌉, and the free K appears only in part C, which is checked interval by interval in balls. The second "gap" (heights below 5·10^6) is item 8, judged correct by the other reviewer.

**Novelty.** The agent searched arXiv listings Oct–Nov 2026 (math.NT, math.CA), the arXiv API (Revers, sub-Weyl, explicit van der Corput), HAL and Zenodo. The orchestrator repeated the listing scan and ran Zenodo phrase queries and an OpenAlex search: the only related item is a verification note for Hiary–Patel–Yang (Zenodo 22952308, constant 0.618). Revers is still v1.

**Caveat.** Uses Lemmas 2.5–2.9 of 2602.05614 (unrefereed, v1), audited correct by revaudit (11 Nov 2026); the same certificate also passes with Hiary–Patel–Yang's refereed Lemmas 2.5/2.6 in their place (max ratio 0.56797, t ≥ 8.97·10^17: 0.482061; [RevAudit](../Riemann/Deep/RevAudit/REPORT.md)), so no unrefereed constant is needed. Python certificate, not formalized.

**RH relevance:** N9 (moments, Lindelöf, subconvexity): the explicit Weyl-type bound on the critical line feeds explicit zero-density and zero-counting estimates; this lowers its constant by 7% for all t ≥ 3 and by 15% for t ≥ 8.97·10^17.

**Files.** [Riemann/Deep/Explicit611/](../Riemann/Deep/Explicit611/) (REPORT.md, certify.py, sanity_authors.py, crosscheck_tail.py, consistency.py, logs).
