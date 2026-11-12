THEOREM: unconditionally, for every λ > 2/3 a positive proportion of consecutive zeta-zero ordinates (with multiplicity) satisfy γ_{n+1} − γ_n < 2πλ/log T, i.e. D(λ) > 0 (Durkan–Pearce-Crump's Thm B argument with polynomial length T^θ, θ < 1/2).

# dpcuncond — Durkan–Pearce-Crump's small-gap theorem without RH

Short answer: the off-line zeros are harmless. Each one adds an explicit "hairpin" term to DPC's moments with weight T^{θ(β−1/2)}, not T^{2θ(β−1/2)}. With Jutila's density slope 1 − δ, those terms vanish as the share of off-line zeros → 0, for every θ < 1/2. The cap θ < 1/2 comes from the unconditional fourth moment of R_T (Montgomery–Vaughan), the only RH-free replacement for DPC's Rodgers-based Lemma 5. With θ → 1/2, DPC's sign condition φ(1+θ) > 1 gives every λ > 2/3. This is the exact limit of the method at θ < 1/2.

Full argument: PROOF.md (Steps 0–8, inputs table).

## Answers to the brief

| Item | Result |
|---|---|
| (1) Contradiction hypothesis | If 𝒩(T) ≤ εTL, then off-line and multiple zeros in [T, 2T] number ≤ 2εTL + O(L) (runs of equal ordinates). With log ξ continued horizontally from Re s = +∞, Im V_φ = N(t + h/2) − N(t − h/2) counts **all** zeros, so DPC's integer count (l.774–785) needs no RH. The price: the lower-half-plane continuation of V has vertical slits |
| (2) One off-line zero | Slits at γ ∓ h/2, depth d = β − 1/2, V jumps ±2i. Exact identity (H), PROOF Step 3: real-line moment = DPC's Im t = −1 main term (unchanged, Euler series) + sides + Σ over merged cut points x of −i∫_0^{depth} G_T R_T R*(x − iy)[P(V_L + J_x(y)) − P(V_L)]dy, with J_x = ±2i for a single zero. Weight R_T R*(x − iy) = M(½ + y + ix)·conj M(½ − y + ix): continuous mean ≍ T^{θ(β−½)}·S(0)·sinh(a)/a. This is a heuristic; the rigorous factor e^{θ(a+1)} comes from PROOF 6c–6e. Growth in j: polynomial in V_L (degree j − 1), no further T-power. Clusters: all coincident shifted cuts are merged, including oppositely signed ones when two ordinates differ by exactly h. The total jump is 2i(m₊ − m₋), and P is applied to the total, never to separate increments; other zeros enter only through V_L |
| (3) Absolute moments | ∫|R_T|⁴ ≪ T by Montgomery–Vaughan (θ < 1/2, roughness bounds Σ|b_n|²). ∫|Im V|^p ≪ T by Fujii (Lz l.242–249, all zeros). ∫|Re V|^p ≪ T: Re V = H[F̃_T] + O(1), where F̃_T is the zero count with each off-line zero Poisson-smeared at scale d. Bound it by M. Riesz + Fujii + any density slope |
| (4) Assembly | THEOREM above. Inputs: DPC (audited correct), Montgomery–Vaughan, Fujii, Lamzouri's DoubleSum argument, M. Riesz, Jutila's slope 1 − δ (or Conrey's 8/7). Hairpin sum ≪ T Σ_a (a+2)^C e^{θ(a+1)} min(ε, e^{−κa})^{E}, E ↑ 1/2, so θ < κ/2 is needed |

## Comparison

| λ such that D(λ) > 0 | Hypothesis | Source |
|---|---|---|
| **> 2/3 = 0.6667** | **none** | **this (θ < 1/2)** |
| ≥ 0.768 | none | explore40 (Lamzouri + hard-core inequality), record before |
| ≥ 0.7851 / 0.7892 | none | smallgap / Lamzouri 2611.14486 |
| 0.775 | none | floor of Lamzouri's criterion (small-gap-floor) |
| 7/11 = 0.6364 | none, heuristic | R448 heuristic for θ < 4/7. Here it would need θ ∈ (1/2, 4/7), which the fourth moment blocks |
| > 0.6878 / > 0.6553 / > 0.6039 | RH | Soundararajan / Wu / BGMM (Lz l.186) |
| > 1/2 | RH | DPC Thm B |

So the unconditional bound now beats Soundararajan's RH bound and lies 0.0114 above Wu's.

## Exact limit and the obstruction beyond 2/3

- On L²[0, θ], Im e^{πA_φ} ⪰ 0 iff φ(1+θ) ≤ 1. The key identity is ⟨1, (I + aS)⁻¹1⟩ = tan(aθ)/a, checked numerically to 1e−6. So DPC's g is optimal, and length T^θ yields exactly λ > 1/(1+θ). θ = 1 gives 1/2 (DPCAudit D7).
- Off-line zeros allow every θ < κ/2: 4/7 with Conrey, i.e. 7/11, and θ < 1 under the density hypothesis. The binding constraint is θ < 1/2, from fourth moments of R_T and of R*(· − iY) on a horizontal line.
- OBSTRUCTION for λ < 2/3: the method needs ∫_T^{2T}|R_T|^{2+δ} ≪ T for rough Möbius polynomials longer than T^{1/2}. No unconditional bound of this kind is known; DPC get it under RH from Rodgers' zero moments. Montgomery–Vaughan gives only (T + T^{2θ})Σ|b_n|².

## Toy-model tests (scripts in this folder, logs in logs/)

| Test | Setup | Outcome |
|---|---|---|
| Identity (H), toy_identity.py → toy_identity.json | F(s) = planted off-line pairs × ∏(s − ½ − iγ_k) × ζ(s+1); Q ∈ {1, 12}; j = 1, 2, 3; one zero and a cluster of three (d = 0.25, 0.1, 0.4) | residual 1e−7 to 1.6e−4 with hairpins (quadrature noise); 0.03 to 38 without them; peak 22 MB, 783 s |
| Weight, toy_weight.py → toy_weight.json | Möbius R of length Q = 100, 1000, 5000; 4000 random x; a = y log Q ∈ [0, 3] | E\|RR*\|/√(S₊S₋) = 0.86–1.03 and E\|RR*\|/S(0) ≈ sinh(a)/a, so the weight is T^{θ(β−½)}; E\|RR*\|/S₊ falls like e^{−a}, which rules out T^{2θ(β−½)}; peak 18 MB, 8 s |
| Threshold, threshold.py → threshold.json | discretised Im e^{πA_φ} on [0, θ], θ = 0.45, 0.5, 4/7, 1 | λ_min changes sign at 1/(1+θ) (sign only; the printed scale is arbitrary); tan(aθ)/a identity to 1e−6; peak 18 MB |

## Novelty

lst.py on arXiv math.NT, math.CA, math.CV, math.SP, math.PR, math.FA, math.CO for 2026-10 and 2026-11 (pattern: small gaps / gaps between / positive proportion / consecutive zeros / spacing / pair correlation) found 13 titles. None is an unconditional version of DPC. 2610.22927 concerns moderate gaps (lower bounds), 2610.31744 biases, 2610.26038 is RH. arXiv API 'small gaps' + zeros: nothing new beyond 2608.07399 (Siegel zeros). Three arXiv web searches: none. Zenodo ("small gaps" zeta; "Durkan" "Pearce-Crump"): none. Repo files (NOGO, DISCOVERIES, RECORDS, Harvest6, R448): R448 lists this as open, with heuristic 7/11.

## Review

Self-review: one full-proof reviewer timed out with no output. Three parallel reviewers (no tools, argument pasted inline) then checked Steps 1–4, Step 5 and Steps 6–7.

| Block | Verdict | Points raised | Repair (PROOF.md) |
|---|---|---|---|
| 1–4 | GAP | (H) sign correct; Im t = −1 main term unaffected; Step 4's weight is a continuous average, not a bound at the zero ordinates (marked FATAL for Step 4 read alone); overlapping cuts need the total jump | Step 4 reworded as the heuristic weight, made rigorous by Step 6; cuts counted once per ordinate with jump 2i·m(y) |
| 5 | GAP, repairable | Poisson tails in 5c; M. Riesz needs m > 1; state the density bound as N(½+δ, 3T) ≪ TL e^{−cδL} | Young's inequality with ‖K_j‖₁ = 1; m ≥ 2; density stated, bin j = 0 by the total count |
| 6–7 | GAP, repairable; "no fatal defect in the sparse-kernel or exponent argument" | ξ̂ zero-free only at heights [T/2, 5T/2]; R*_Y tail keeps its growth factor; q, w > 1; order of choices | localised boxes plus harmonic measure of the sides; tail bound absorbed in A₂; E = (1−1/w)(1/2 − 1/(3r)) confirmed; order (δ, w, r) → ε → T |

All points are folded in.

**Orchestrator review (orchestrator 33).** Three reviewers and five falsification scouts found no fatal error.

| Reviewer | Scope | Verdict |
|---|---|---|
| A | Steps 0–3, 7 | GAP, repairable; SIGNIFICANT if closed |
| B | Step 5 | CORRECT; SIGNIFICANT |
| C | Steps 4, 6 | GAP, repairable; SIGNIFICANT if closed |

The scouts confirmed four things:
- Lz l.221–230: Conrey's and Jutila's density estimates are uniform for σ ≥ ½.
- Lz l.240–284: Fujii's bound counts all zeros, for h ≤ 100/log T.
- DPC uses RH only at l.57, 59–62, 140, 607, 610 and 629, never inside Prop 1 (l.262–388) or Lemma 6 (l.686–741) after the shift.
- DPC's construction (l.743–771) needs only 0 < φθ < 1 < φ(1+θ) < 2 and θ < 1.

A numerical test of the hairpin weight at points sampled at actual zeros (Q = 500, 10,000 zeros) matches the random-point weight.

All eight repairs are applied and listed at the end of PROOF.md ("Repairs (orchestrator 33 review)"):

| # | Step | Repair |
|---|---|---|
| 1 | 0 | b₀ is a real polynomial of degree 2 deg q + 1, from DPC l.757–759; θ < 1 only (l.744, 753) |
| 2 | 3 | coincident shifted cuts merged, oppositely signed included; total jump 2i(m₊ − m₋); tip circles 0; sides a, b ≥ c/L from shifted ordinates |
| 3 | 5a | Montgomery–Vaughan main term on [T/2, 5T/2] is 2T |
| 4 | 6b | box sides ≥ c/L from shifted ordinates (V̂ uniformly T^{O(1)} there); v truncated to I |
| 5 | 6a, 6d | Ñ_ρ kept inside the polynomial factor; bins over all of Z_off; Ψ moments by h-cells, Hölder, Young, Minkowski |
| 6 | 6f | order: κ, E (κE > θ), A = c log(1/ε) with c < E/θ, ε, T |
| 7 | Thm U ⇒ D(λ) | 2/3 < φ < min(λ, 1), apply on [T/2, T] |
| 8 | 4 | continuous-mean heuristic; rigorous factor from 6c–6e |

## Resources

One process at a time, single-threaded, all under memcap.sh 290. Peaks: 22, 18, 18, 0 MB. No worktrees, no commits, no Lean. Downloads deleted at the end.

RH relevance: N8 (zero spacing).
