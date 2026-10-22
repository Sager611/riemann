DISCOVERY: Λ ≤ 0.147373586543538, from Planat's unchanged producer and auditor (PASS; 18,864 steps; all 9 of his negative controls rejected). The inputs are his sources with our 786 certified kernels and with source floors raised to new certified paired-probe floors. Premises that remain: A2 (wall), A3 (first contact) and the analytic half of A1 (Raymond's model-to-function transfer). Our certificate replaces only the numerical moment-cone and paired-credit step of A1.

# Raising the A1 floors (2026-10-21)
**Previous record.** Our own Λ ≤ 0.148409 (REPORT5), on the same premises; Planat's posted bound is 0.148612. The new bound is lower by 0.00104.

**Where the slack was.** A1 bounds S̄, the average of −Im H'/H at the two probe points x ± rh. Raymond bounds each point by a worst-case moment cone. He credits only constant-label pairs with n ≤ 1024 for the fact that both points cannot be worst at once. Yet every label's phase at the two points differs by exactly e^{±i(rh/2)λ} (reflected labels to within 1e-11), so every pair of labels is constrained.

**New bound** (for each cell, height, shift and frequency subrange):
- For each pair of labels, the two-probe moment form has a closed-form minimum over the relative phase: −|α e^{−iθ} + β e^{iθ}|.
- A split grid then bounds (g₊ + g₋)/2, the mean of the two probes' Re(V/Q): if g₊ ≥ c_k then g₋ ≤ 2c − c_{k+1}, each step proved by Ψ = H_{c'}(z₋) + μH_{c_k}(z₊) ≥ 0.
- The paired Euler term is an exact one-dimensional supremum; it increases with ρ (Poisson averaging).
- Labels that switch on between N₀ and N are explicit linked labels; a cutoff index present at one probe only is charged as unlinked. For s ≥ 8/5 we use Raymond's envelope (8.1)–(8.3).
- Floor = L/4 − 10⁻⁶ − Raymond's physical error (9.1, with our own A₊, B₊) − paired Euler term − c/2. The cell floor is the minimum over 7 frequency subranges.

**Rigor.**
- Every transcendental value comes from Arb. Bin masses are upper bounds from exact integer counting and log-convexity. They use Raymond's lemma that head coefficients are positive and nondecreasing in time.
- All other arithmetic is interval arithmetic with IEEE outward rounding. Every certificate was verified when it was produced.
- A separate re-run re-verified 48 of 48 randomly sampled floors that the producer consumed.
- Tests: interval arithmetic against exact rationals; masses never below direct sums (about 1% above); the Euler bound never below a fine grid; on 14 toy systems the certificate is never below a brute-force maximum.

**Gains.** 10,150 of 12,160 floors raised. Mean gain by cell range:

| cells | 0–99 | 100–199 | 200–299 | 300–399 | 400–607 |
|---|---|---|---|---|---|
| gain | +0.06 (retuned provider only) | +0.08 | +0.12 | +0.12 | +0.10 |

Example: cell 500, r = 3/5 goes from 5.9044 to 5.9997. Early cells keep the higher old floors of the enhanced provider.

**Limits.**
- The relaxation lets every pair of labels pick its own phase. Real phases are multiplicative (χ(pq) = χ(p)χ(q)), and using that would lower Φ, the probe average of Re(V/Q), from about 0.85 to about 0.48 at cell 500 (about 0.18 more in floors). That needs a new argument.
- Losses inside the method: grid steps (about 0.001) and mass bins (about 1%).

**Novelty.** A record only inside this conditional framework; no external result below 0.1486 is known (REPORT5 search). The certificates are computer-generated; the analytic premises are unrefereed.

**Cost.** About 1.6 CPU-hours (3 niced processes, at most 0.25 GB each). Scratch peaked under 10 MB.

**Files** (EVIDENCE2/):
- floors_dbn5.json and certs_dbn5.json.gz: the floors and their certificate parameters.
- classical_dbn5.json.gz and its AUDIT; 07_planat_run_dbn5.txt; logs 08–10.
- scripts/ and README.md: how to rerun.
