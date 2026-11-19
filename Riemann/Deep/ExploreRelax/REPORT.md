BARRIER: rows of off-line zero pairs that obey the heat-flow zero ODE (extreme case: the heat flow of cos(πρz) − cosh(πρy₀)) pass all 67 of Raymond's certified sources while reaching the line only at t = 0.1018, and pass every phase-uniform σ_eff > 1 input until (16/L²)·log cosh(L(2θ₁₄₁ − 1)/4) = 0.0957, so the record 0.10755 can drop by at most 5.3% on its sources and 11% overall (node N2).

# explore51 — Dyson Brownian motion / local relaxation flow → N2 (de Bruijn–Newman)
Orca agent explore51. Started 14 Nov 2026 21:45 UTC; the session stalled for 5 days after 20 minutes; work resumed 19 Nov 18:53 UTC. Scripts in this folder, outputs in logs/. All processes: one core, peak ≤ 64 MB RSS (/usr/bin/time -l; memcap peaks ≤ 47 MB).

## 0. Transfer
- Tool: the heat-flow zero dynamics (Csordas–Smith–Varga; Polymath 15 Prop 3.1), read as β = ∞ Dyson Brownian motion. Two DBM notions are tested: Erdős–Schlein–Yau's local equilibrium (relaxation of the environment of a particle) and Bourgade–Erdős–Yau–Yin's coupling of two configurations that have the same density.
- Normalisation: H_t(z) = ∫ e^{tu²}Φ(u)cos(zu)du, ∂_t H = −∂_z²H, H_0(z) = ξ(1/2 + iz/2)/8. Zeros obey ż_k = 2Σ'_j 1/(z_k − z_j). An on-line zero is z = 2γ. Zero density at height Z: ρ(Z) = log(Z/4π)/(4π). At Raymond's wall X = 6000000185827: ρ = 2.139978, spacing s = 0.467295, and πρ = L/4 = 6.722939 with L = log(T₀/2π) = 26.891756.
- Top zero z = x + iy, q = y²: q' = −2 − 4yΣ_other (y − b_j)/|z − z_j|², with b_j = Im z_j. The −2 is the conjugate (de Bruijn). Raymond's chain bounds the sum below by a probe functional J (signed Poisson combinations at heights h > y, with σ_eff = (1+h)/2 + tL/8 > 5/4) and certifies J ≥ S(t) for ζ from the Dirichlet series, uniformly in x ≥ X.
- Why ESY convexity does not transfer: on the real line the flow is gradient ascent of 2Σ log|x_i − x_j|, a convex log-gas. There ESY/BEYY relaxation holds, with conductances 2/(x_i − x_j)²; Rodgers–Tao use it for Λ ≥ 0, and so does Kastoris (Zenodo 23842668). Off the line the flow is a saddle flow: ẋ = +∂_x Re W, ẏ = −∂_y Re W, with W = 2Σ log(z_i − z_j). Nothing drives complex configurations to the real equilibrium; they have a slow invariant manifold, the row below.

## 1. Statement
- **Exact slow solution.** Take an infinite row of pairs x_m ± iy with pair spacing a = 2/ρ (every zero off the line, density ρ). It is the zero set of f(z) = cos(πρz) − cosh(πρy). Under ∂_t f = −f'', the row becomes f_t = e^{(πρ)²t}cos(πρz) − cosh(πρy₀), so cosh(πρy(t)) = cosh(πρy₀)e^{−(πρ)²t}. Equivalently ẏ = −πρ·coth(πρy). It reaches the line at
  t_B = log cosh(πρy₀)/(πρ)² = (16/L²)·log cosh(Ly₀/4) = 4y₀/L − 16 log 2/L² + O(e^{−Ly₀/2}/L²).
- **Worst case for the sum.** For a given density, a zero at the top height contributes the least: its kernel is 2/(d² + 4y²) per zero, against 2/(d² + y²) for a real zero. So the row descends at half the mean-field speed of a pair surrounded by real zeros: πρ against 2πρ for y ≫ s. This kernel-level fact is Raymond's "main mass ≤ 1/y" (DBNAudit §1). The row realises it as an actual solution.
- **Input-consistency.** Replace ζ's zeros in a window above T₀ by such a row; 64 pairs (128 zeros) already give t_B to 2·10⁻⁶ (§2). The row's signature in a probe at height h is one Fourier mode at frequency πρ, of relative size 1 − [tanh(πρ(h−y)/2) + tanh(πρ(h+y)/2)]/2 ≈ e^{−πρ(h−y)}. ζ's own phase-uniform band for the same probe is (1/2)Σ_p log p·p^{−σ_eff}. It always contains the signature (§4): margin ≥ 3.15 along the whole trajectory for every single-height probe, hence for every positive untranslated combination. The row also stays below the σ_eff = 1 line y = 1 − tL/4 (above it, Dirichlet-series data could exclude zeros): that line moves at L/4 = πρ ≤ πρ·coth(πρy) and starts higher.
- **Barrier.** Take any argument whose information above T₀ is (i) the strip |Im z| ≤ y₀ at t = 0 and (ii) phase-uniform σ_eff > 1 data that the row satisfies. Such an argument cannot prove Λ < t_B: the modified configuration satisfies the same inputs and is still off the line at every t < t_B. This holds for every comparison quantity, Lyapunov function or probe design the row passes. The barrier is model-level: the modified function is a Hadamard product with ζ's zeros outside the window, not a function with an Euler product; it only has to meet the inputs. The class contains Raymond's chain: his 67 signed, translated designs pass the periodic row with margin +0.81 (§4, raycheck).
- Values at T₀ = 3·10¹² (logs/analysis.out):

| strip y₀ | t_B (ceiling) | best proved |
|---|---|---|
| classical β < 1, y₀ = 1 | 0.133409 | 0.156821 (Planat–Raymond), 0.1486 (Planat, given A1–A3) |
| 7/8, y₀ = 3/4 | 0.096223 | 0.107965 (Raymond), 0.13963 (Hakimoglu-Brown) |
| θ₁₄₁, y₀ = 0.746691712 | **0.095731** | **0.107551** (record, HP141) |

- KKL made explicit: no density-level proof of Polymath's Theorem Zero(i) ("x ≥ e^{C/t} ⇒ real") has C < 4(2θ − 1); with the classical strip, C < 4. For the θ₁₄₁ ceiling to reach 0.09, 0.08, 0.05 or 0.02, RH must be verified to T₀ ≈ 10^13.4, 10^15.2, 10^25.0 or 10^64.0.

## 2. Tests: the smallest case, then scale (row_exact.py, logs/row_exact.out)
- Method: put N zeros per period P. Then F = Π 2 sin(π(z − z_j)/P) is a trigonometric polynomial, and the flow multiplies each Fourier mode by e^{(2πm/P)²t}. This is exact (no ODE step). Its zeros obey ż_k = 2Σ (π/P)cot(π(z_k − z_j)/P), which is the zero ODE of the P-periodic continuation. Hitting time = the first t with N sign changes per period (bisection to 10⁻⁷; grid 2¹⁷).
- Smallest case: N = 2, one pair per period a = 2s. Trig flow 0.0957314033; formula 0.0957314031.
- Check against an ODE: N = 4 (one pair and two real zeros), RK4 with dt = 10⁻⁶. max Im z at t = 0.02: ODE 0.5481656508, trig 0.5481656508 (agree to 3·10⁻¹⁵). Hitting time 0.073178.
- Scale: N = 400 zeros per period at ρ(X). k pairs at height y₀ = 2θ₁₄₁ − 1 replace 2k consecutive real zeros, each pair at the midpoint of the two zeros it replaces. Environments: lattice, CUE(400) (two seeds), and ζ zeros n = 1001…1400 (mpmath.zetazero, unfolded by θ(t)/π + 1, rescaled to spacing s). Hitting times:

| k pairs | 1 | 2 | 4 | 8 | 16 | 64 | 200 (all) |
|---|---|---|---|---|---|---|---|
| lattice | 0.074832 | 0.081686 | 0.091564 | 0.094956 | 0.095622 | 0.095730 | 0.095732 |
| CUE seed 1 | 0.076551 | 0.085075 | 0.095111 | 0.101237 | 0.104068 | 0.104243 | 0.104988 |
| CUE seed 2 | 0.070818 | 0.077436 | 0.086584 | 0.091131 | 0.097067 | 0.101181 | 0.104933 |
| ζ statistics | 0.077444 | 0.083811 | 0.092421 | 0.095823 | 0.096580 | 0.097206 | 0.097508 |

- Reading:
  - A single pair among real zeros (ESY's "local equilibrium") reaches the line by 0.071–0.077, against de Bruijn's y₀²/2 = 0.278774. Even certified knowledge that every neighbour is real would give only about 0.075 at this height. The reason: y₀ = 1.6s is not ≫ s, so the mean-field speed 2πρ = 13.4 is never reached. With a gap g ≈ s around the pair, ẏ ≈ −1/y − 4ρ·arctan(y/g), about −10 at y₀.
  - A periodic row converges to t_B with 2k ≳ 30 zeros (k = 16: 0.095622).
  - Irregular rows are slower: their slowest pair sits in a local gap (up to 0.104988). §4 shows the probes exclude them.

## 3. Raymond's chain against the solutions (analysis.py, logs/analysis.out)
- Comparison path: the DBNB2 run (b0 = 0.87338196661…, endpoint 0.107559925), damping.txt. The rows are rerun at the same y₀ = 2b0 − 1. q = (max Im z)²:

| t | chain | periodic row | CUE row s1 | CUE row s2 | ζ-stat row | one pair (CUE) |
|---|---|---|---|---|---|---|
| 0.0125 | 0.452508 | 0.439189 | 0.451012 | 0.450827 | 0.441307 | 0.394228 |
| 0.025 | 0.358836 | 0.334818 | 0.355649 | 0.355368 | 0.338628 | 0.262578 |
| 0.05 | 0.204695 | 0.168120 | 0.198443 | 0.198163 | 0.173860 | 0.087292 |
| 0.075 | 0.091965 | 0.055927 | 0.084657 | 0.084471 | 0.061345 | 0.003205 |
| 0.1 | 0.016452 | 0 | 0.010623 | 0.010500 | 0 | 0 |
| hit | 0.107560 | 0.095742 | 0.104999 | 0.104944 | 0.097519 | 0.076558 |

- Falsification test passed: every solution stays below the certified path at all times. The CUE rows come within 0.0015 in q at t = 0.0125. So early on the chain is nearly sharp for GUE-irregular configurations.
- Rates. Write the row's rate as r_row(q) = (2πρ√q·coth(πρ√q) − 2)/q. Raymond's band rates (7/8 run, DBNAudit §1) are 0.874, 0.880, 0.870, 0.861, 0.860, 0.817, 0.826 and 0.775 of r_row at the top of each q band, from 0.5625 down to 10⁻⁴. No chain with inputs the row satisfies can have a rate above r_row(q) at any q. The room left is 12–22% in rate and 0.0118 in Λ (11%).


## 4. Input-consistency: which configurations the probes can see (analysis.py, analysis2.py, analysis3.py)
- ζ's band for a height-h probe at t = 0: J_ζ(x) = πρ − (1/2)Σ_n Λ(n)n^{−(1+h)/2}cos((x/2)log n) (explicit formula; mean πρ = L/4). By Kronecker the phase-uniform infimum is πρ − Pen(h), with Pen(h) = (1/2)(−ζ′/ζ((1+h)/2)), so every certified S satisfies S ≤ πρ − Pen. A configuration passes a probe if its relative deficit 1 − min_x J/πρ is at most Pen/πρ.
- At t = 0, for the rows at y₀ (logs/analysis.out §d):

| h | Pen/πρ | periodic row | CUE row s1 | ζ-stat row | CUE, all real |
|---|---|---|---|---|---|
| 1.05 | 2.93 | 0.115 | 0.288 | 0.235 | 0.152 |
| 1.36 | 0.373 | 0.016 | 0.160 | 0.095 | 0.124 |
| 2.0 | 0.112 | 2.2e−4 | 0.093 | 0.020 | 0.093 |
| 2.5 | 0.0649 | 7.6e−6 | **0.078** | 0.0076 | **0.077** |
| 4.0 | 0.0215 | 3.2e−10 | **0.049** | 0.0011 | **0.048** |
| 6.0 | 0.0075 | 4e−16 | **0.030** | 2.5e−4 | **0.029** |

- Whole trajectory. Use only primes p ≤ N₀ = √(x/4π) = 690,989 (a lower bound for Pen) and the H_t exponent (1+h)/2 + (t/4)log(x/4πp) (Polymath's b_n^t). Probe heights run down to σ_eff = 1.02, at t ≤ 0.09. The periodic row passes every single-height probe with margin min Pen/deficit = 3.15 (closed form, analysis2) and 3.35 (exact J_t = −Im F_t′/F_t from the evolved coefficients, analysis3).
- Irregular rows with GUE fluctuations at every scale fail at h ≥ 2.5. So does an all-real CUE configuration. ζ's zeros are prime-rigid at large scales (the counting function smoothed over ≳ 2 is governed by p = 2, 3, 5, …; Berry's saturation), and the chain's high probes (Raymond's row 0 uses heights 2.5–6) enforce exactly that.
- Rigidified CUE rows (CUE pair midpoints with every displacement mode of wavelength > l_c pair spacings removed), checked along the whole trajectory:

| l_c | hit (seeds 1, 2, 3) | min Pen/deficit |
|---|---|---|
| 8 | 0.099287, 0.100225, 0.099886 | 1.56, 1.75, 1.57 |
| 12 | 0.101807, 0.099493, 0.100594 | 1.56, 1.66, 1.36 |
| 16 | 0.101419, 0.102754 | 1.28, 1.02 |
| 24 | 0.101554 | 0.80 (fails at h = 6) |

- So sources built from single-height probes cannot prove Λ < 0.1018 (margin 1.56). That is 5.3% below the record.
- **Raymond's actual sources (raycheck.py, logs/raycheck.out).** Data: his 67 certified rows, verification/arcs/ACTUAL_67_ROW_ARC_INPUTS.json from Zenodo 23752103 (Proof_and_Replay.zip, sha256 bf5c7b3d…; 442 KB extracted, then deleted). Each row gives signed translated nodes (h, v, w), jets on rows 0, 16 and 44, and the affine source S_i(t). For each configuration, J_i(x, t) was computed from the evolved coefficients: G = −Im F_t′/F_t matches the closed form to 10⁻¹⁰. The minimum over 8192 x and three times per row of J_i − S_i, with v of either sign:

| configuration | hit | min J − S | at the top zero |
|---|---|---|---|
| periodic row | 0.095732 | +0.812 | +0.812 |
| CUE row l_c = 8 | 0.099287 | +0.371 | +0.453 |
| CUE row l_c = 12 | 0.101807 | **+0.249** | +0.303 |
| CUE row l_c = 16 (seed 2) | 0.102754 | −0.051 (t = 0.0975) | +0.185 |

- So a configuration that satisfies every one of Raymond's source inequalities (and the strip, and lies above X) takes 0.1018. His chain, or any re-optimised comparison path on these sources, cannot prove Λ < 0.1018. To go below, new source certificates must exclude such rows; the floor for any phase-uniform design stays 0.0957.
- The run stopped at its 1700 s cap (peak 45 MB) after 4 of 7 configurations, slowed by the jet image sums.
- pathcheck.py: all three rigidified rows stay below his certified path at every tabulated time (e.g. l_c = 12: q(0.05) = 0.1881 against 0.2047). This is a consistency check of his kernel lemma; it passed.
- Residual risk (designs other than Raymond's). The row's signature in a design Σ_l w_l P_{h_l} is at most Σ_l |w_l|·2πρ·e^{−πρ(h_l − y)}·cosh-factor. That is exponentially small in the probe–row distance; for Raymond's row-0 design (9 nodes, heights 2.5–6, |w| ≤ 16.9) at t = 0 it is ≲ 10⁻³ of the mean. The primes p ≈ N₀ (v_p = log p/2 ≈ πρ) give every design a band contribution ≥ (1/π)∫e^v|R̂(v)|dv near πρ, of ratio ≈ W·e^{πρ(1−y₀)}/(πρ) ≈ 1.6–3 to the signature when R̂ is smooth there. A design whose transform has a peak narrower than ~0.6 at πρ could in principle see the row. An LP over all signed designs failed numerically (HiGHS status 4) and is not used.

## 5. Adversarial (brief item iv)
- Is this only a lemma about the ODE? It is a statement about ζ's chain: the modified ζ (same zeros below T₀, same strip, a row in one window above T₀) has the same density and passes the chain's inputs. A density-matched model is exactly what the chain cannot tell apart from ζ. That is the content of the barrier.
- Does it tell ζ apart from Davenport–Heilbronn? Yes: the ceiling uses RH to T₀ and the Euler-product strip θ₁₄₁, and DH has neither (its off-line pair 0.8085 + 85.70i sits at low height; Λ_DH > 0.2306).
- The positive direction (a relaxation lemma) does not survive: the row is itself a solution of the ODE and consistent with the inputs, so no estimate valid for every configuration can beat t_B.
- Chain validity: the certified comparison path dominates every solution tested (§3). The kernel inequality forces this: a solution that is slower and has J ≥ S along its path would refute the kernel lemma. None was found. The irregular rows that come within 0.0015 of the path fail the probes at h ≥ 2.5 (§4).


## 6. Novelty (brief item v; 19 Nov 2026)
- arXiv API queries (logs/novelty.txt): "de Bruijn-Newman" with upper / limit / barrier / heuristic / dynamics / Dyson / Calogero / local / trigonometric / verification, and "heat flow" with zeros / trigonometric polynomial. Monthly listings Oct–Nov 2026 (math.NT/CA/CV/SP/PR/FA/CO; titles matching Bruijn|Newman|heat flow): only seen items (2611.06245, 2610.12803, 2610.18214, 2610.27343, 2609.37164).
- Zenodo "de Bruijn-Newman" (167 records, newest read):
  - Kastoris 23842668 (31 Oct) re-proves Rodgers–Tao's local-equilibrium step through a force heat equation with conductances 2/(x_i − x_j)²: the real-line DBM tool, used for Λ ≥ 0.
  - Moki & Julio 24318337 (17 Nov): degree-2 bounds only.
  - Neither gives a ceiling.
- Classics:
  - Polymath 15 (1904.12438 TeX): Theorem Zero has an ineffective C; its "gas, liquid, solid" remark describes the cooling picture with no worst case.
  - Rodgers–Tao (1801.05914): lower bound only.
  - Planat 2609.37164: collective contraction, no ceiling.
  - Raymond 23752103: the kernel-level main-mass ceiling 1/y, not a configuration.
  - Ki–Kim–Lee and Csordas–Smith–Varga: qualitative.
  - The trigonometric example cos z − c with Λ = log c is elementary and surely classical.
- Our files: DBN/REPORT2 gives the Polymath-criterion limit ≈ 0.144–0.145 for one method. It is consistent with this report: the classical-strip row ceiling 0.1334 lies below it and binds every density-level method. No row, configuration ceiling or log cosh formula appears in NOGO, TREE, Riemann/Tree or Riemann/Deep.
- The brief's "Dobner 2005.13434" is a quantum-computing paper (wrong ID). Dobner's Newman-conjecture work concerns Λ ≥ 0.
- Verdict: the input-consistent row as a ceiling for ζ's chains, the formula (16/L²)log cosh(L(2θ−1)/4) with C ≥ 4(2θ − 1) in Theorem Zero, and the measured 11% gap to the record look new. The single-probe ceiling 0.1018 is a model computation. No reviewer verdict yet.

## 7. What would go lower
- The ceiling moves only with the inputs. Raise RH verification: t_B = 0.09 needs T₀ ≈ 10^13.4 (today 10^12.5). Sharpen the strip: dt_B/dθ ≈ 8/L = 0.30. Or add information that sees a row above T₀: a short-window zero-density bound at β ≈ θ, or height-specific evaluations of H_t near the line, which phase-uniform data cannot supply.
- Inside the chain's own class, the largest possible gain is 0.0118, and probably at most ≈ 0.006 for single-height-style sources (§4). Done here (raycheck.py): against his sources the ceiling is ≥ 0.1018, so re-optimising his arcs, stopping rule or comparison ODE can gain at most 0.0057. Only new sources can go further, and none can go below 0.0957.
- Relaxation/DBM: even certified local equilibrium (every neighbour of the top pair real) gives only ≈ 0.071–0.077 at T₀. Λ ≤ 0 remains RH itself (t_B → 0 only as L → ∞ or θ → 1/2).

## Suggested NOGO line
- explore51 (Orca, exploration: local relaxation flow / Dyson Brownian motion → N2) BARRIER: rows of off-line pairs (heat flow of cos(πρz) − cosh(πρy₀), ẏ = −πρ·coth(πρy)) pass the chain's inputs. A rigidified row passes all 67 of Raymond's certified sources and takes 0.1018, so his chain cannot go below 0.1018. For any phase-uniform σ_eff > 1 data from θ₁₄₁ and RH to 3·10¹², the floor is (16/L²)log cosh(L(2θ−1)/4) = 0.0957. Theorem Zero needs C ≥ 4(2θ−1). Certified local equilibrium would give only ≈ 0.075 [ExploreRelax](Riemann/Deep/ExploreRelax/REPORT.md)

## Files
- row_exact.py: exact periodic heat flow (trigonometric polynomial), the N = 2 and N = 4 checks, the N = 400 environments, hitting times. Output logs/row_exact.out and .json.
- zeta_env.py: 400 ζ zeros (n = 1001–1400), unfolded. Output logs/zeta_unfolded.json, zeta_env.out.
- analysis.py: the ceiling table, the comparison with Raymond's path and rates, single-probe visibility at t = 0. Output logs/analysis.out and .json.
- analysis2.py: the periodic row against the band at t > 0; rigidified rows at t = 0. Output logs/analysis2.out and .json.
- analysis3.py: consistency along the whole trajectory from the evolved coefficients. Output logs/analysis3.out and .json.
- raycheck.py <rows67.json>: configurations against Raymond's 67 certified sources. Output logs/raycheck.out. To get the input: download Zenodo 23752103 Proof_and_Replay.zip, then unzip -p it 'Proof_and_Replay/Lambda_Peer_Review/verification/arcs/ACTUAL_67_ROW_ARC_INPUTS.json'.
- pathcheck.py: rows against his certified path. Output logs/pathcheck.out.
- logs/novelty.txt: the novelty queries and their hits.
- Replay (one core, < 64 MB each): uv run --no-project --with numpy [--with mpmath] python3 Riemann/Deep/ExploreRelax/<script>.py, from the repo root; zeta_env.py first (~5 min). The others take 10 s–3 min.

RH relevance: N2

## Orchestrator review (19 Nov 2026, orchestrator 45)
- Rerun of analysis2.py on a copy: identical (overall min ratio 3.15; rigidified rows identical), memcap peak 48 MB.
- Reviewer: GAP, ROUTINE. Sampled source checks (8192 x, three t per row) do not certify uniform bounds in x ≥ X and t; single-height margins do not cover signed designs (the narrow-peak case stays open); the window splice (strip, Hadamard structure, tail and source bounds) and the spliced collision time are not certified; "data the configuration satisfies" is a conditional restriction, not a ceiling for all phase-uniform methods; (a) is the classical cos(az) − c example. A method floor would need a certified feasible witness for the complete relaxation.
- Gate: not a discovery; NOGO line and a TREE N2 barrier note (model ceiling).
