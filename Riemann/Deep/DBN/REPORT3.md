NO DISCOVERY: no false step in Planat's Λ ≤ 0.148612 (arXiv:2610.18214). It is an explicitly conditional implication, with premises A1–A3 imported from D. J. Raymond's unrefereed Zenodo deposits. Every independent test passed.

# Audit: Λ ≤ 0.148611789321711 (Planat, arXiv:2610.18214 v1, 13 Oct 2026)
**How it gets below REPORT2's ≈0.1705.** It does not use Polymath's criterion (Theorem 1.2) at all. It uses Raymond's collective-contraction comparison (Zenodo 22867547, 21 Sep 2026, Λ < 0.149396; retuned release 23272725, 10 Oct 2026).
- **Exact identity** at a simple highest zero x+iy, with q = y²: q' = −2 − 4y·Σ_other K_y(x−a, b), where K ≥ 0. I re-derived it from Polymath Prop 3.1.
- **Probe.** Lower-bound the average S̄ of −Im H_t'/H_t at x ± rh. The probe height h ∈ [0.69, 2.18] lies far above the zeros, so the series converges absolutely (amplitude exponent ≥ 1.334).
- **Kernels.** 57 rational kernels give 2K_y ≥ λ(y/h)·J. Hence q' < −2 − cq with c = (2λ/h)(ℓ − M/h), where ℓ is a floor on S̄ valid for all x ≥ X.
- **Comparison.** An exact ODE runs from q0 = 1.000002 (the classical strip) down to 1e−8: 19,023 steps, rates 10.2 to 22.2. De Bruijn then adds 5e−9.
- **Consequence.** No zero-free region at small y is needed, so REPORT2's phase-torus limit (≈0.144) does not apply.

## Inputs and rigor
- **RH height.** Only X/2 for X = 6000000185827 (Gomila's X). The Platt–Trudgian margin is 175,239,886; nothing above Platt–Trudgian is used.
- **Wall (A2).** |H_t/B_t| > 3/4 on x = X for 0 ≤ y ≤ 1 and t ≤ 0.164; y = 0 is included, so no real zero crosses X.
  - Source: 4·690988 prefix checks (3.1) > 79/100, plus errors ≤ 1/100 + 1/180, in Raymond's Arb at 256 and 384 bits.
- **Floors (A1).** 608 time cells × 2 providers × 6 frequency classes, computed in Raymond's two Arb/FLINT lanes (one shared library).
  - They come from the moment cone C|Q|² − Re(V·conj Q) ≥ (1+A)[C(1+A)−B] − D_C; I checked the identity.
  - Planat reassembled the stored evidence; he did not regenerate it.
- **Error terms.** Polymath Prop 6.1 (stated for real σ) and Prop 6.3, re-derived for Re s ∈ [−1, 2]. The remainder is ≤ 4·exp(−(1+y)L/4 − tL²/16 + 1/10); totals are ≤ 2e−6 in value and 2e−4 in derivative.
- **Contact compactness (A3).** Polymath Thm 1.5(i) (zeros with x ≥ e^{C/t} are real) and Hermite splitting (Prop 3.1(ii)), both cited correctly. Kernels and caps are Bernstein-certified; the scalar history is in exact rationals.
- **New in this paper.**
  - Both providers and all 57 kernels in the classical strip give 0.14863355.
  - A finite Weyl tail on cells 320–595 subtracts a further 0.0000218. It uses Yang's Lemma 2.1 plus the A-process; its constants match Yang (arXiv:2301.03165v2).

## First unjustified step (none false)
Premise A1, "S̄_{t,h_{p,i},r}(x) > ℓ_{p,i,r} (t ∈ I_i, x ≥ X)", is assumed, not proved; for example ℓ = 6.321162 at cell 0. A2 and A3 are assumed in the same way. Planat says so himself: "their complete mathematical validation has not been independently completed in this audit". So this is a conditional result, not an error.

## Independent checks (all pass; .scratch/dbn3/RESULTS.txt)
- **Wall (3.1)**, with my own 160-bit phases: the minimum prefix is 0.800597 (at J = 1, t = 0.0205) against the claimed > 0.79. Full sums are 31.5–39.7.
- **Floors against the phase infimum.** I minimised S̄ with L-BFGS over all 55,874 prime phases plus the reflected phase, using the two-branch model at X.
  - Margin (infimum minus floor) by cell: 0: 0.025; 25: 0.030; 100: 0.050; 200: 0.124; 300: 0.293; 500: 0.300; 594: 0.253.
  - At cell 0, an exact Euler decoupling gives 6.34557, against the optimiser's 6.34592.
  - No floor exceeds what the phases allow.
- **Height extension.** I computed |H/B − f| from heat-kernel ξ integrals at x = 2·10³ and 2·10⁴, for y from 0.69 to 2.41 and t = 0, 0.05, 0.15. It is 0.09–0.16 times the claimed bound.
- **Kernels.** A grid check of all 57: the slack tends to 0 only as d → ∞ for λ = 2, where ρ² ≤ 1 − 3r² holds. All self caps are valid.
- **Scalar replay** in 50-digit mpmath: every restart is at least the exact endpoint (excess ≤ 2⁻⁹⁰), and U = 0.14861178932171. Floors, rates and eligibility all match.
- **Polymath slip (already known).** In the proof of Polymath Prop 6.3, "1/(12(T'/2−0.33)) ≤ (1/3)/(4(T−3))" is off by a factor of 2. Raymond's release already fixes it (5/6 → 11/12); it is harmless.

**Not checked:** Raymond's native production of the 608 cells; the activation and paired-credit algebra beyond reading; Planat's 192- and 256-bit tail run; the 7/8-conditional appendix.
**Lead:** with each consumed floor raised to its measured phase infimum (same kernels and heights), the replay gives about 0.1458. The slack sits in the cone relaxation at late cells, so a sharper rigorous source bound (for example branch-and-bound over small-prime phases) is the route to a record in this framework.
**Files** (.scratch/dbn3/): `wall_check.py`, `torus_min.py`, `euler_inf.py`, `height_check.py`, `kernel_check.py`, `hist.py`, `ideal_replay.py`, `RESULTS.txt`; sources in `src/`, `ray1/`, `ray2/` and `yang/`. CPU about 20 min; peak 0.40 GB.
