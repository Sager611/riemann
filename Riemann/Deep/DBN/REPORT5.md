DISCOVERY (conditional on the premises A1–A3 of Planat's posted Λ ≤ 0.148612, arXiv:2610.18214): Λ ≤ 0.148408913066829 is the bound backed by Planat's unchanged producer and auditor, with the new kernels proved by separate exact certificates. The sharper Λ ≤ 0.148368783506957 uses two-shift mixtures and is backed only by our superset of his auditor (all his checks plus mixture and kernel-geometry checks). Evidence: EVIDENCE/README.md.

# Pushing Planat's classical-strip comparison below 0.148612 (2026-10-21)
**What sets 0.148612.** The bound is the time for q' = −2 − c·q to take q from 1 to 10⁻⁸. Here c = (2λ/h)(ℓ − M/h), where ℓ is a source floor (A1), λ a kernel constant, M the own-pair cap and h the height. Float replays (`sens.py`, `q0sens.py`) give these sensitivities:
- **Floors.** dU/dℓ = −0.016 per unit floor, spread evenly over the first three quarters of the run. Raising floors to the measured torus infimum (REPORT3) would give 0.14556 (a quarter of the way: 0.14767).
- **RH height.** Doubling the RH height X would give 0.14576. Starting at q₀ = 0.9752 instead of 1 would give 0.14697.
- **Kernel catalogue (the free lever).** Floors exist for 10 shifts r, but Raymond's 57 kernels give a family over ρ (where y ≤ ρh) for only 2 of them. Early and middle cells sit at y/h ≈ 0.44–0.46, where the best λ drops steeply in ρ (r = 0.6: 1.980 at ρ = 0.42, 1.890 at ρ = 0.45). So the comparison was forced onto a coarse ρ or a worse shift.

## New certified inputs (exact rationals; no new transcendental input)
1. **786 single kernels.** These cover all 10 declared shifts, with ρ spaced 10⁻⁴ on [0.40, 0.48] and more coarsely below; λ = numerical optimum − 3·10⁻⁶.
   - Each is proved by exact integer Bernstein subdivision (`certify_kernel.py`; 9,510 leaves) and again by Planat's own algorithm (`check_geometry.py`; 13 kernels need its depth cap raised from 30 to 34).
   - The polynomial identity 2K_y/y − λJ = 4P/(Q*D₄) was re-derived in sympy with residual 0 (`derive.py`).
   - The own-pair caps are checked exactly.
2. **Run with Planat's unchanged producer and auditor.** The inputs are his refined sources with these kernels appended; the floors are byte-identical (checked). U = **0.148408913066829** in 18,997 steps; the audit passes and all 9 negative controls are rejected (`run_phase1.py`).
3. **345 two-shift mixtures** (same provider, height and cell), one per cell for cells 0–344. The kernel is 2K_y ≥ (y/h)(λ₁J_{r₁} + λ₂J_{r₂}), with rate (2/h)Σλ_k(ℓ_k − M_k/h); A1 is used only at its own provider, height, cell and shifts.
   - The degree-9 polynomial identity checks with residual 0 (`derive_mix.py`).
   - The exact Bernstein certifier (`certify_mix.py`) reduces correctly to the single kernels, and its sympy and dictionary polynomials agree.
   - A direct float check of the kernel inequality on every mixture gives a minimum ratio of 1.0000015.
   - My own auditor (`audit_mix.py`, imports no producer code) re-certifies every mixture and passes: U = **0.148368783506957** in 18,992 steps (9,012 of them on mixtures). All 6 negative controls are rejected.
4. **Gains.** Planat's own Weyl-tail gain was 0.0000218. The finer catalogue gives 0.0002029 and the mixtures a further 0.0000401, for **0.0002430** in total.

## Status and novelty
- **Premises.** Wall, contact argument, heights, floors, cells, seed q₀ = 1.000001² and terminal 10⁻⁸ are all unchanged. A1–A3 are the same premises Planat assumes; REPORT4 judged them proved (computer-assisted, unrefereed).
- **Literature.** As of 2026-10-21, the arXiv listing and a Zenodo search show nothing below 0.148612 after 2610.18214. Raymond's last deposit is dated 10 Oct. The 0.13963 and 0.116398 bounds assume no zeta zeros with Re s > 7/8, so they are not comparable.
- **Limit of this lever.** Comparison-layer gains are now nearly used up. A per-cell mixture LP shows no gain after cell 345. Per-substep mixtures and 64 substeps per cell each add about 10⁻⁶.

## Leads (not done)
1. **Split-abscissa start, about −0.0014 (to about 0.1470).** The zero-free region (R₀ = 5.558691) gives q ≤ 0.9752 for zeros with x ≤ Xe², but not as x → ∞. The fix is to treat x ≥ Xe² separately using its higher class floors (those classes sit about 0.4–1.2 above class 0 at late cells). It needs three things:
   - a wall at the split for y ≳ 0.9 and t ≲ 0.02;
   - class floors for the early cells, from Raymond's 123 MB `Proof_and_Replay.zip`;
   - a short written multi-region contact lemma.
2. **Sharper source floors, up to −0.0029.** These would come from branch-and-bound over small-prime phases in place of the cone relaxation, which leaves a slack of about 0.3 at late cells.

**Cost.** About 35 CPU-min, one niced process (17 min of it a single sympy simplify); peak RSS about 0.1 GB; scratch 172 MB, of which the venv is 170 MB.

**Files** (`dbn4/`): scripts, `catalogue.json` (786 kernels), `mixed_kernels.json` (345), both histories (`.json.gz`), and the audits. To rerun, work in `.scratch/dbn4` with its `venv` and run `gen_catalogue.py`, `run_phase1.py`, `produce_mix.py` and `audit_mix.py`.
