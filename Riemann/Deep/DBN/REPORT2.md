NO DISCOVERY: the Polymath criterion at the Platt–Trudgian height cannot plausibly beat the lowest posted bound, 0.148612 (Planat, arXiv:2610.18214). Its best certificate type reaches about 0.1705; its ideal phase-uniform limit is about 0.144–0.145. No complete certificate below 0.1779532546 was produced.

# Lower de Bruijn–Newman bound inside the Polymath framework (Theorems 1.2 and 1.3)
**The bar (posted, Oct 2026):** 0.148612 (Planat, arXiv:2610.18214, 13 Oct, under stated premises); 0.149345 (Planat, 2609.37164 v2); 0.149396 (Raymond); 0.15737 and 0.15766 (Gordon lanes); 0.172422 (Kornaus); 0.17795 (Planat, relative to Gomila); 0.1787854 (Gomila). Conditional: 0.13963 (arXiv:2610.12803, assumes no zeta zeros with Re s > 7/8). X is pinned by Platt–Trudgian (X ≤ 6.00035e12), so the first window is N0 = 690988.

## Findings
1. **Gomila's lane has slack.** His finite lane uses the plain triangle weight |b|+g|a|, which gives a floor of 7.9137e-7 at N0 (reproduced in Arb). Polymath's own improved triangle inequality (proof of claim (b)), with weight max(|b−ga|, c|b+ga|), gives **0.0769953** at the same row (Arb, 128 bits). This weight is equally rigorous for pointwise non-vanishing: dist(Ef, (−∞,0]) ≤ |Ef|.
2. **First-window frontier.** This is a necessary condition, explored in float64 with floor ≥ 1e-6, as min t0+y0²/2 over y0. Five-prime Euler mollifier (2,3,5,7,11; "P11"): 0.17167. Seven primes (P17): 0.17091. Nine primes (P23): 0.17050. The optimum is at y0² ≈ 0.005–0.006 and t0 ≈ 0.1675–0.169. Gains shrink with each added prime. Rescaling the mollifier coefficients does not help: λp = s·bp is best at s = 1. The floor rises with N and with y (spot checks up to N = 2e6 and y² = 0.66), so the first window at y0 is the binding one.
3. **Arb-certified first-window floors** (improved weight):
   - (t0, y0) = (0.158, 0.1809), P11: floor **0.0052723**. Candidate Λ ≤ 0.174362405. Gomila's barrier covers this row (it needs y0 ≥ 0.1809 and t0 ≤ 0.16125).
   - (t0, y0²) = (0.1677, 0.006), P23: floor **0.0025779**. Candidate Λ ≤ 0.1707. This needs a new barrier: for t in [0.1543, 0.1677] down to y ≈ 0.077, and at all heights for t > 0.16125.
   - The approximation error at these rows is about 3–4e-7.
4. **Ideal limit (Bohr-torus obstruction).** Any certificate built only from coefficient sizes (all triangle and mollifier types) must hold for every choice of prime phases. Over all 55k prime phases up to N0, I minimised Φ = |F| − g|G|. The function f has a torus zero wherever Φ < 0.
   - Φ < 0 was found at Λ = 0.1445 (y0² = 0.001 and 0.005), 0.145 (0.002, 0.01), 0.146 (0.02), 0.1475 (0.035), 0.151 (0.05) and 0.157 (0.08).
   - So no such certificate exists below about **0.144**. The record sits about 0.004 above that edge.
   - At Λ = 0.1485 (y0² = 0.002) the torus minimum is Φ ≈ 0.21, but the best floor (P23) is −0.30. Beating the record needs a new certificate type that recovers most of that 0.5 gap over windows N0 to about 5·N0. No such method is known.
   - At Gomila's row the torus minimum is 0.485, against a certified floor of 0.097.
5. **Planat's collective contraction** adds at most about 0.001–0.002 at y0 ≈ 0.18, and nothing at small y0. It cannot close the gap.

## Cost (laptop, one niced process)
- **This run:** about 45 min CPU; peak 103 MB (Arb), 31 MB (torus optimiser); scratch 2 MB.
- **Full certificate for 0.174362 (relative to Gomila's barrier):**
  - grouped Arb floors on the bottom edge y = y0 for N0…N1, with the improved weight;
  - left, top and right edges via the argument principle, as in Polymath claim (b);
  - an all-N tail bound (Polymath's large-N lemma or Gomila's tail lemma);
  - estimate: 2–4 CPU-hours plus about a day of code and audit.
  - Not run: it is not a record, since 0.172422 and lower bounds are posted.
- **0.1707 row:** additionally needs a new barrier (Gomila's used 883 time prisms). Estimate: tens of CPU-hours.
- **Below 0.148612:** unrealistic now.

## Leads
- A branch-and-bound over small-prime phases, plus grouping by the large prime factor, could approach the torus limit. This is the only route toward a record in this framework.
- The 0.148612–0.172422 lanes were not audited here. If they fell, rows 0.1707 and 0.1744 would be the best bounds in the classical framework.

**Files** (`Riemann/Deep/DBN/dbn2/`): `floor.c` float64 explorer, `floor_arb.c` Arb certifier, `bohr.c` torus optimiser, `scan*.py` and `bohr_cont*.py` drivers (run next to the compiled binaries), `RESULTS.txt` all numbers, `arb_rows.log` Arb output. Scratch: `.scratch/dbn2/`.
