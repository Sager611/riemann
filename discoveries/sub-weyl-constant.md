# Explicit sub-Weyl bound: |ζ(1/2+it)| ≤ 54.6 t^{27/164} (node N9)

**Claim.** For all real t ≥ 3, |ζ(1/2+it)| ≤ 54.6 t^{27/164}. For t ≥ e^2000 the constant is 52.36; on ranges it is smaller (35.58 on [e^75, e^100], 45.08 on [e^175, e^200], 50.42 on [e^275, e^300]). The bound is sharper than every published Weyl-type bound (c t^{1/6} log t forms) for t ≥ e^{88.7}, and than our 0.568 t^{1/6} log t page's large-t form for t ≥ e^{101.3}.

**Before.** 66.7 t^{27/164} for t ≥ 3 (Patel–Yang, arXiv 2302.13444 = J. Number Theory 262 (2024), Thm 1.1); 307.098 (Patel 2021). Nothing better on arXiv, Zenodo or HAL as of 6 Nov 2026.

**Error found in the previous proof.** Patel–Yang's Lemma 3.5 (TeX l. 876, 921–923) prints Σ_{k≤K} a_k^{−α} < μ₂(α) h^{αK} t^{−α/2}, but the last equality drops a factor (2πh)^α; with α = 17/328 their C₂ piece is too small by (2πh)^{17/328} ≈ 1.10, their own parameters give 68.54, and 82 of their 163 table rows (log t from 380 to 785) exceed 66.7. Our 54.6 implies their Theorem 1.1, so the statement stands (REFUTATIONS).

**Method.** Patel–Yang's lemmas (1.2 Weyl differencing, 1.4 Yang's k-th derivative test, 2.1–2.3 Poisson/B-process, 3.1 Riemann–Siegel, the display in the proof of Lemma 3.4, TeX l. 826–831) with two changes: no square-root split in Lemma 3.4 (TeX l. 860), the C₁ and C₂ pieces summed inside one square root per block (worth 65.8 → 54.2 at log t = 540); and per-block parameters θ₁, η₁ with, per block, the better of the full Lemma 3.4 bound and the 4th-derivative test. The block domain is derived without the threshold θ₂t^{7/17} (REPORT, "Block domain (derivation)"). Several typo-level slips of the paper are repaired (REPORT).

**Certificate.** certify.py (python-flint arb, 128 bits; 241 s, 61 MB): [3, e^75] from Hiary–Patel–Yang's 0.618 t^{1/6} log t (53.983 at e^75); [e^75, e^2000] by 325 interval evaluations, every term a monomial bounded by its exact supremum over the interval, every lemma hypothesis asserted per interval; maximum 54.59967 at log t ∈ [648, 650]; tail t ≥ e^2000 by a closed formula (52.358). The orchestrator reran it on a copy (6 Nov 2026): ALL PARTS OK, identical maximum. crosscheck_terms.py matches each term to Lemma 3.4 (largest log ratio 4.5·10⁻¹³); sanity_paper.py reproduces the paper's 163 rows and the μ₂ correction.

**Review.** The agent's reviewer (three rounds, evidence inline): COMPLETE. Orchestrator's reviewer: first GAP (the extension to blocks with a ≤ θ₂t^{7/17} was asserted, not derived); after the derivation was added, CORRECT (threshold dependence removed, every certificate block covered including the last, the 4th-derivative alternative holds without the threshold, C₁ and C₂ summed in one square root); SIGNIFICANT if validated (about an 18% lower constant; same exponent; incremental method).

**Novelty.** lst.py (math.NT, math.CA, Oct–Nov 2026 titles: no explicit sub-Weyl paper), Zenodo "sub-Weyl" (3 unrelated hits), arXiv API ("27/164", "sub-Weyl": only 2302.13444 is explicit; 2609.00537 improves the B-process and says it could refine Patel–Yang but gives no constant), OpenAlex (6 citing papers all use 66.7, none notes the slip), HAL (nothing), our files (R376 candidate only).

**Caveat.** Conditional on the cited lemmas of the refereed paper (checked against arXiv v1; the journal host blocks scripts) and on python-flint/arb.

**RH relevance:** N9 (explicit sub-Weyl/subconvexity constants on the critical line; inputs to explicit zero-density and S(T) bounds that cite 66.7).

**Files.** [Riemann/Deep/SubWeyl/](../Riemann/Deep/SubWeyl/) (REPORT.md, certify.py, run_final.log, review.log, crosscheck_terms.py, sanity_paper.py).
