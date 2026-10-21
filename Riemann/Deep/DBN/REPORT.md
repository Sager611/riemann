DISCOVERY: in arXiv:2609.37164 (Planat, v1 and v2), the lower-barrier certificate (Theorem "Directed lower barrier after tau_0") evaluates the wrong function: its node values use F1 + gamma*conj(F2) instead of Polymath's F1 + gamma*F2. All four printed edge minima match that slip to within 2e-12; the true values differ by up to 0.35. Not fatal: corrected, the margin is still 0.312, so the Lambda <= 0.177953 chain stands. Gomila's Lambda <= 0.1787854: no gap found.

# Audit: Lambda <= 0.1787854 (Gomila) and arXiv:2609.37164 (Planat)
**Sources.**
- Gomila: github.com/judegomila/dbn-lambda-01787854-candidate-audit; primary paper dated 20 Aug 2026 (not 2025, not on arXiv, not peer reviewed).
- arXiv:2609.37164 is by **Michel Planat**, not Gomila (v1 29 Sep, v2 7 Oct 2026). It imports Gomila's certificate and claims Lambda <= 0.1779532546 ("audit-relative").
- Inputs: Polymath arXiv:1904.12438v2 (Theorems 1.2 and 1.3) and Platt–Trudgian 2021. Texts are in `.scratch/dbn/{gom,pm15,planat_v1,planat_v2}`.

## Gomila's argument, reconstructed
- **Criterion.** Polymath Theorem 1.2 at X = 6000000185827, t0 = 129/800, y0^2 = 87677/2500000, giving Lambda <= t0 + y0^2/2 = 0.1787854.
- **(i) RH verified to X/2** = 3000000092913.5, below the Platt–Trudgian height 3000175332800 (margin 1.75e8).
- **(ii) Zero-free region:** H_t0 != 0 for x >= X + sqrt(1 - y0^2), y0 <= y <= sqrt(1 - 2t0). Windows N = 690988..3840000 use a triangle lower bound on |f_t| with Euler mollifiers ({2,3,5,7,11} down to {2,3}); smallest floor 7.91366e-7 at N = 690988. All N >= 3.84e6 use a separate tail lemma (D < 0.999721).
- **(iii) Barrier:** the closed box [X, X+1] x [0.1809, 1], t in [0, t0], certified in 883 time prisms; allowance 0.00125 (actual error <= 3.57e-4); smallest margin 0.5198.
- **Error terms:** Polymath Theorem 1.3, H/B = f_t + O(e_A + e_B + e_C0), with the two-leg f_t (not the older A + B - C, E1–E3 form). e_C0 dominates: about 2.335e-7 at the weakest row.
- No heuristic step found. Fragile point: the finite lane rests on a mollified mass of 0.9999986 < 1 (absolute margin 5.6e-7).

## Pinned error (Planat section 4: R_w = [X, X+1] x [0.08, 0.1809], tau in [tau0, tau_w])
- **Printed** node minima (bottom, top, left, right): 1.3453859649899646, 1.6143370467511939, 2.1748966688549862, 1.4259781421906720.
- **Conjugated leg** reproduces them, differences +1e-14, +1.4e-12, -4e-13, +1.9e-12 (`variants2.py`, variant e). **True** Polymath (14) gives 1.2302290964, 1.6997759691, 2.1357351192, 1.7726685868. The bottom-edge value at x = X + 8/12, tau_w equals a direct mpmath sum of all 690988 terms to 1e-15.
- **Which is right:** at x = 1e6, y = 0.08 my f stays within Polymath's error bound of H/B (H from a heat-kernel integral of xi): 0.0095 against a bound of 0.0133. The conjugated variant misses by 3.35.
- **Consequence:** the printed margin 0.42726 is unproven. With true nodes and Planat's own bounds: 1.23023 - 0.64128 - 0.27684 - 1.5e-6 = 0.3121 > 0; a fine scan gives true min Re f about 1.2302 on the boundary. The barrier and Lambda <= 0.177953 relative to Gomila survive. This is a broken proof step (REFUTATIONS.md material), not a counterexample.
- **Novelty:** not found reported (web search; Gomila's issues). The other v2 lanes (Gordon 0.15737, Raymond 0.149345, Kornaus-hybrid 0.1718/0.16836) were not audited and may share this evaluator.

## Strongest checks (independent code; all agree)
- **Criterion text:** Polymath Theorem 1.2 matches Gomila's (i)-(iii) exactly. Its proof (Prop. 3.3) also gives Planat's height bound M(tau0) < y0.
- **Weakest row:** my numpy code gives L_690988 = 7.913667e-7 (stored 7.91366e-7) and E = 2.33493e-7 (stored bound 2.33495e-7).
- **Window and tail:** x* lies in the first window, margins 5.38e6 and 1.20e7. At the tail cutoff: D = 0.9997198 < 0.999721, M_max = 1.608289, (1-D)/M_max = 1.742e-4.
- **Barrier mesh minima reproduced:** prism 1 4.2783620904; prism 883 1.5431990755993192 (stored ...189515). So Gomila's evaluator uses the correct second leg; the conjugated variant is off by 6e-6.
- **Planat's algebra:** pair kernel, bridge numerator (c^2 - r^2)(c^2 - r^2 - 4) and transport kernel checked exactly in sympy; (1/2)Re(M'/M) = 6.72293909. True L_0.99 near X is 6.45–6.49 (G >= 4.4) against the certified > 3.679. The two-envelope comparison holds.
- **Not re-run:** the full 3.15M-row sweep, every prism's derivative bounds, the 7688 Taylor coefficients.

**Files** (`.scratch/dbn/`): `floor.py` (floor); `feval.py` and `theta_build.py` (f_t at X, Arb phases); `validate_f.py` (H/B check); `check_gomila_mesh.py`; `variants2.py` and `scan_lower.py` (Planat nodes); `misc_checks.py` and `sym_tail.py`. CPU about 8 min; peak 0.9 GB (transient).
