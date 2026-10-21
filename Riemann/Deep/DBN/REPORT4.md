NO DISCOVERY: Raymond's premises A1–A3 behind Planat's Λ ≤ 0.148612 are written, computer-assisted proofs, not heuristics. Every identity and inequality I re-derived holds, and every new numerical test agrees. I found no false or unjustified premise.

# Audit: premises A1–A3 of arXiv:2610.18214, from D. J. Raymond's Zenodo deposits
**Deposits.** Found with the Zenodo search `creators.name:Raymond AND Bruijn`; his other records are unrelated physics.
- **22626161** (8 Sep 2026): Λ < 0.163297 via the Polymath criterion. Planat does not use it.
- **22867547** (21 Sep): Λ < 0.149396 from the classical strip. Its *Analytic Foundations* (Part A, §2–12) contains the proofs of A1–A3.
- **23272725** (10 Oct): a retuned source and 44 new kernels.
  - It proves Λ < 0.116398 only under the hypothesis that no zeta zero has Re s > 7/8.
  - Its Appendix B fixes the Gaussian bound.
  - It ships the `SOURCE_RETUNED` floors that Planat uses.

## The premises as Planat states them, and verdicts
**A1 (source floors).** For each provider, cell and shift, the average of −Im H_t'/H_t at x ± rh, at height h, exceeds ℓ for all t in the cell and all x ≥ X (for example ℓ = 6.321162 at cell 0).
- **Verdict: proved**, computer-assisted.
- **The proof chain:**
  - Raymond extends the Riemann–Siegel error to Re s ∈ [−1, 2] and gets derivatives by Cauchy's estimate on a disk with frozen cutoff (Foundations §4–5).
  - Then come an exact moment-cone identity (§6), Taylor time bounds (§7) and integral envelopes for every future cutoff (§8–9).
  - Then a surplus decomposition and a paired-phase lemma (Correlated Probes §2–6), and the positivity c_t ≥ 0 (Supplement §3).
- I re-derived each of these; all hold.

**A2 (wall and localisation).** |H_t/B_t| > 3/4 on x = X for 0 ≤ y ≤ 1 and t ≤ 41/250. With Platt–Trudgian's verified RH, a homotopy argument then rules out nonreal zeros with |Re z| ≤ X.
- **Verdict: proved** (Foundations §2–3).
- The ingredients are prefix sums > 79/100, an Abel-summation transfer to all y, and Appendix B's error bound 2/369 < 1/180.

**A3 (first contact).** Built from Polymath Thm 1.5(i) (compactness for t > 0), Prop 3.1 (zero motion, with Hermite splitting ruling out a multiple contact) and de Bruijn's strip bound.
- **Verdict: proved** from published theorems.
- The constant in Thm 1.5 is ineffective, but only compactness is used, so this is harmless.

## New independent checks (`.scratch/raymond/work/RESULTS.txt`)
- **Masses.** I recomputed them in numpy: Euler-mollified head, product tail, reflected branch, and the actual cutoff with new terms. This covered cells 320, 450 and 595 in all 6 classes, plus cells 400 and 500 in classes 0–4.
  - Every value lies inside Raymond's enclosures, with slack ≥ 5.5e-4 below and ≥ 1.9e-3 above, and slack ≥ 0.020 on the weighted mass.
- **Floor formula (9.2)**, recomputed from the stored budgets with my own Euler cap, for all 3,312 budgets (2 providers × cells 320–595 × 6 classes):
  - Each stored floor is 1.0e-7 below my value (deliberate rounding).
  - The exact cone margin is ≥ 1.0e-7 > 0, and A+ ≤ 0.294 < 4/5.
- **True function at t = 0** (mpmath at height X/2):
  - |H_0/B_0| on the wall is 80.99, 23.72, 11.77, 7.26 and 5.07 at y = 0, ¼, ½, ¾ and 1, against the claimed > 0.75.
  - At cell 0, the probe average is 6.4356 against the floor 6.321162.
- **Positivity lemma:** c_t(n) ≥ 0 for every n ≤ 690988 at t = 1e-4, 0.05, 0.1 and 0.152.
- **Toy check of the zero-motion identity.** A toy polynomial under ∂_t P = −P'' gives dq/dt = −13.14954 numerically, against −13.14964 from the identity. −Im P'/P equals Σ K_h exactly.
- **From REPORT3.** The floors lie below the phase infimum at X (margins 0.025–0.30), the wall prefix is 0.800597, all 57 kernels pass, and the replay gives 0.14861178932171.

## Errors Raymond already corrected himself (none affects Planat)
- **Polymath constant:** 5/6 → 11/12.
- **September wall:** used a narrower historical error field; it is replaced by the broad bound.
- **September "fixed displacement" explanation** for the sharp continuum kernel is false. With r² = 4/15 held fixed at w = 1/10, the self term is 2100/1321, above the cap 30/19 (Supplement §8).
  - Planat's history uses fixed shifts only (`uses_continuum_shift: false`; 19,023 steps over 6 fixed r).

## Remaining trust, not proof
- **Not regenerated independently:** the native production of all 608 cells (two lanes of one Arb/FLINT library) and the per-cell paired credits. I checked only the one-point budgets of cells 320–595.
- **Unrefereed:** the analytic text is unrefereed and AI-assisted.
- **Tightest margins:** the cone conditions (1e-7, tight by design) and cell 0 against the phase infimum (0.025).

**Files** (`.scratch/raymond/`):
- Small files and unzipped sources from the three deposits.
- `work/{masses,masses2,floorcheck,wall_true,toy_checks}.py`.
- CPU about 5 min, niced; peak 0.68 GB.
