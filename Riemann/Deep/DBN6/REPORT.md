NO DISCOVERY: the bound to beat under the 7/8 theorem is not 0.13963 but D. J. Raymond's Λ ≤ 0.116397726408291 (Zenodo 10.5281/zenodo.23272725, 10 Oct 2026, restated in Planat arXiv 2610.18214). Every route in the brief models at 0.1218–0.1392. The one lever that could go lower, to about 0.115, needs a uniform bound on the vertical derivative of S that nobody has certified.

# De Bruijn–Newman below 0.13963 under the 7/8 theorem (2026-10-21). Details: NOTE.md; scripts: evidence/. Model numbers are not certificates.

**1. The bar is 0.116398.**
- Raymond, "Retuned translated probes and the conditional bound Λ < 0.116398", v1.0, 10 Oct 2026. It assumes no zeta zero with Re s > 7/8, the same hypothesis as Hakimoglu-Brown. It is computer-assisted and unrefereed.
- What it rests on: his own wall at X₀ = 6000000185827 (|H_t/B_t| > 3/4 on x = X₀, 0 ≤ y ≤ 1, t ≤ 41/250); Platt–Trudgian's RH check up to 3000175332800; his source floors for all frequencies on 608 time cells; and 57 exact kernel certificates.
- Planat's arXiv 2610.18214 (13 Oct) restates the bound and ships the source files as ancillary data.
- Our floating-point replay of the chain from those files gives 0.11639772640829064, matching his exact ceiling to 15 digits (evidence/raymond_model.out).
- "Our own certified computations" (allowed by the brief) have the same status as his. So any c ≥ 0.1164 is not a record.
- The NOGO entry that treats 0.13963 as the bar should cite 0.116398 instead.

**2. Routes in the brief.** These come from a model of Hakimoglu-Brown's chain that uses his formulas with X as a parameter. It reproduces his table: C₀ 5.2904 vs 5.2902, λ₀ 1.723 vs 1.718, k₀ 8.670 vs 8.638, final 0.13943 vs 0.13963.

| Route | Wall used | Model bound |
|---|---|---|
| (a) finer steps, best (λ, μ) per step | Polymath at 6·10¹⁰, height 1/5 | 0.13921 |
| (b) own box wall at 6·10¹² (not computed) | height 1/5 | 0.12408 |
| (c) Polymath Table 1 wall, X = 5·10¹²+194858 (peer-reviewed) | y ≥ 0.16733, t ≤ 0.186 | 0.12350 |
| (c) Raymond's line wall (lands at height 0) | unrefereed | 0.12179 |

Why all stay above 0.1164: Hakimoglu-Brown compares at height p = 61/20, because his approximation is only proved for 3 ≤ y ≤ 6. Even at 6·10¹² this caps the rate at k ≈ 11.3–13.6. Raymond's probes sit at heights 0.69–1.88 and reach rates of 11–24.

**3. Levers on Raymond's chain** (model run on his exact source data):
- *Denser kernel catalogue* (830 kernels, 10 shifts, ρ step 0.005): 0.1163974618. The gain is 3·10⁻⁷, so this lever is dead.
- *Vertical derivative.* This is Hakimoglu-Brown's Lemma 10 applied to translated probes; Planat piloted it only at the first cell, with μ = 1/500.
  - My scan raises λ from 1.99 to 2.26 at y/h ≤ 0.3 and to 2.38 at y/h ≤ 0.2.
  - Model bounds for an assumed uniform bound G on the mean ∂_hS at the probes: G = ψ₂(σ_eff)/4 gives 0.11233; the phase-free heat-deformed ideal G gives 0.11474; ideal G + 0.1 gives 0.11505.
  - A record near 0.115 is therefore plausible, but only with a certified, all-frequency G at heights 0.69–1.88 for t ≤ 0.117. That is a second-derivative counterpart of Raymond's whole source producer. Planat says it is not certified, and existing floors cannot be differentiated to get it.

**What blocks:** (1) the bar, 0.116398, is below every route the brief asks for; (2) the only lever below it needs that new physical certificate, which is days of work and was not attempted.

**Not done:** no barrier computation (route b cannot give a record) and no new certificates (no bound is claimed). Raymond's 0.116398 deposit has not been audited; we audited only the premises of his 0.149 result.

**Novelty scan:** arXiv math.NT, CV, GM and CA for 2026-09 and 2026-10 (evidence/novelty.txt). The Zenodo API returned 403 today; this morning's cached search shows nothing after 23272725.

**Resources:** one nice'd process at a time, about 20 CPU-minutes in total, peak RSS 139 MB, scratch up to 50 MB (deleted).
