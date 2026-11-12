AUDIT PASSED — Theorem U (D(λ) > 0 for every λ > 2/3, unconditionally) survives a line-by-line referee: complete proofs of Steps 5c and 6b–6f (FULLPROOF.md) close every step; seven cosmetic corrections, no error.

# uaudit — line-by-line referee of Theorem U (dpcuncond, PROOF.md)

## Short answer

The proof is correct. FULLPROOF.md gives every inequality of Steps 5c and 6b–6f, and states what each constant depends on: φ, θ, f, η, P and the auxiliary exponents m, p, q, r, w, δ. No constant depends on ε, on T or on the zero configuration.

- **Slit identity (H).** It holds with merged signed jumps in all 9 adversarial configurations, to 10⁻¹¹. The rule PROOF.md rejects (separate increments per zero) fails on every coincident cut, so the merging is necessary and is the right rule.
- **DPC's main term.** After the contour shift it uses no RH, for every θ < 1 and DPC's rough Möbius coefficients.
- **Uniformity.** Step 5 never uses the contradiction hypothesis.
- **Exponents.** E = (1 − 1/w)(½ − 1/(3r)) and the order of choices both check.
- **Reviewer.** CORRECT on the most delicate chain: the box representation of V̂, then the Δ and Ψ moments, then the bound for A₃.

## Brief items

| Item | Verdict | Where (FULLPROOF.md) |
|---|---|---|
| (1) Complete proofs of 5c, 6b, 6c, 6d, 6e (+ 6f) | done; no false step | §1–§5 |
| (i) Identity (H), orientation, clustered and coincident cuts | correct: Ω_δ is simply connected; walking b → a goes down the right side of each slit and up the left, giving −i∫(F_R − F_L); P is applied to the total jump 2i(m₊ − m₋); tips, intermediate tips and corners contribute O(r\|log r\|^D) → 0; numerics below | §6 |
| (ii) DPC's main term after the shift: no RH; θ < ½; rough coefficients | correct: the only zero uses are Lemma 1 and the move at l.341. After it, only the Euler series on Re s = 3/2 and the zero-free completion factor G (shifted to Im t = ∓H) appear. Support T^θ < T/L^{20} (l.702). Lemma 6 uses Mertens and Laplace transforms, never Lemma 5 | §6 |
| (iii) Step 5 uniform in ε | yes: Step 5 never uses Step 2; c₂ in Step 7 is ε-free | §6 |
| (iv) Poisson/box representations with log boundary singularities | correct: shrink to (x₁, x₂) × (−1, −δ) and pass to the limit by dominated convergence with the uniform log majorant; top density ≤ P_y (maximum principle); sides ≤ (4/π)e^{−π(T/8−1)} | §2 (b4) |
| (v) Sparse-kernel Hölder 6c: near and far parts, overlaps, tails | correct (Lemma C): near part by the power mean with N_ρ; far part by dyadic shells and Minkowski; Hölder on S_a with F2d; Poisson tails ≪ T^{3θ/2}σ_a(a+1)²/L | §3 |
| (vi) Δ and Ψ moments 6d | correct, after the Ñ window correction (2 below); the near bound needs the two logs of each ρ' paired | §2 (b5), §4 |
| (vii) Exponents 6e–6f, order of choices | correct: (A₁A₂A₃)^{1/3} ≪ T(a+2)^{D+2}e^{θ(a+1)}σ_a^E; Φ(ε) ≪ (log 1/ε)^{D+3}(ε^{E−cθ} + ε^{c(κE−θ)}) → 0 | §5 |
| Steps 1, 3, 7: hidden assumptions | none: Step 1's arg count is exact (on the line (s − ρ)(s − 1 + ρ̄) = −d² − u² < 0); Step 7 uses only 4th moments | §6 |

## Corrections (all cosmetic; FULLPROOF §8)

1. 6b: I = [7T/8 − 1, 17T/8 + 2], not +1. 5c is needed on I, not only on [T, 2T]; it holds on any subinterval of [3T/4, 9T/4].
2. 6c/6d: the window of Ñ_ρ (4Y_a + h) misses near ρ' with |γ' − γ_ρ| < 5h/2 when a = 0 and φ > 0.42. Use 4Y_a + 3h; the moments are unchanged.
3. 6d: the constant 32Ψ_ρ becomes (2/φ)Ψ_ρ, and (a+2)^C becomes (a+2)^{m+1}.
4. 6b: the near bound |ℓ| ≤ C + Σ log⁺(h/|t − w|) holds only with the logs of each ρ' paired. Unpaired, each log would cost log(1/h) ≍ log L.
5. Step 1: "1 − ρ̄ adds 0" is a jump count; the continuous arguments of the pair cancel.
6. Step 3: a zero with slit in [a, b] but γ ∉ [T/2, 5T/2] is not in Z_off. It contributes o(1); alternatively take a ≥ T/2 + 1.
7. Step 3, slits with γ ∉ [T, 2T]: read the log majorant in y, using |s − ρ''| ≥ |y − d''|. In that reading it is uniform in the configuration; read in x it would not be.

## Falsification numerics (scripts here, logs in logs/)

| Test | Setup | Outcome |
|---|---|---|
| (H): slit_identity.py | X = e^{αs}∏(s − ρ)^m, complex cubic P, Möbius-type R (length 30), Gaussian weight. 9 configurations: single zero; γ₂ − γ₁ = h with equal depths (J ≡ 0) and unequal depths; two off-line zeros at one ordinate; a double off-line zero; an on-line zero at the off-line ordinate; an on-line shifted ordinate at a slit top; a mixed cluster of 6; deep zeros near the sides | residual 2.5·10⁻¹³–1.5·10⁻¹¹; corner cases 5.6·10⁻⁹ and 1.0·10⁻⁹ (tanh-sinh); 3·10⁻⁴–50 without hairpins; separate increments 0.16, 10.6, 6.6 on the coincident cuts. Peak 50 MB (138 MB mpmath) |
| Lemma C on actual zeros: sparse_kernel.py | zeros in [2000, 4000] (1960); subsets of 5% and 30% (random, tightest local clusters, contiguous block); a = 0, 2, 5; mixed or full depths; q = 2, 4 (72 runs) | ‖K‖_q^q/((a+2)^qσT) ≤ 0.113 even with σ instead of σ^{1−1/w}; ‖K‖_q^q ≤ 1.01× the near-part bound; max N_ρ ≤ 5. Peak 50 MB |
| 5c on actual zeros: vmoments.py | T = 2000, φ = 0.7: the far part Re V_exact − H[F̃] (mpmath log\|ξ\|, 120 points), and moments over [T, 2T] | far part: mean −0.106, sd 0.114, max 0.315 (O(1) as claimed). ∫\|Re V\|^m/T = 0.295, 0.384 (m = 2, 4); ∫\|Im V\|^m/T = 0.601, 0.811, 1.652 (m = 2, 4, 6); M. Riesz ratio ‖Re V‖_m/‖F̃‖_m = 0.50, 0.70, below Pichorides' cot(π/2m) = 1, 2.41 |
| 5c with planted off-line pairs: vmoments_planted.py | J = [T, T + 400]; 5% or 30% of zeros replaced by pairs at d = 1/L, 3/L, 0.45 | per unit length, m = 2/4/6: \|Re V\| 0.31/0.41/1.28 unplanted, 0.30–0.38/0.38–0.48/1.15–1.28 planted. Smeared Im V̂ 0.57/0.72/1.35 → at most 0.82/1.46/3.73. Sharp Im V (pairs twice) at most 1.10/3.98/23.2. All stay O(1) per unit length. Peak 50 MB |

Numerics test identities and constants at T ≤ 4000, not asymptotics. No test came close to breaking a claimed inequality.

## Inputs (TeX re-read)

- **Lz l.221–230.** Conrey's N(σ, T) ≪ T^{1−(8/7−ε)(σ−½)}log T holds uniformly for σ ≥ ½ (full proof in Baluyot's thesis). Selberg's exponent is ¼ and Jutila's 1 − ε. Jutila's suffices: with κ = 1 − δ and E ↑ ½, κE > θ is possible for every θ < ½.
- **Lz l.240–249.** Fujii: ∫_0^T(N(t + h) − N(t))^{2k} ≪ T for 0 ≤ h ≤ 100/log T. N counts all zeros with multiplicity. Fujii's 1974 range was corrected in 1981 (Lz footnote).
- **Lz l.251–284.** DoubleSumZeros: 2π/L-cells, ab^q ≤ a^{q+1} + b^{q+1}, then Fujii with Cauchy–Schwarz. FULLPROOF F2c–F2d use the same argument at windows O((a+2)/L).
- **DPC's RH uses.** l.57–62, 140, 158 (Lemma 1's hypothesis), 233 (the Section 3 standing hypothesis, used at l.335–341), 607, 610, 629. R450's list omitted l.158 and l.233; these are exactly the uses Steps 1 and 3 replace.

## Review

Reviewer (gpt-6-astra, no tools, chain pasted inline): **CORRECT** on 6b–6d and A₃. It checked the analyticity window of ξ̂ against the boxes, the Poisson representation with top-edge log singularities, the harmonic-measure bounds, the far and near Δ bounds (the slit owner is among the near ρ', and its endpoint log is integrable), the Ψ chain (TLe^{−κa'/4}) and the Jensen–Hölder assembly. Its one remark — use min(σ_a, C) — is already in the definition of σ_a (FULLPROOF §5, 6f).

## Residual risks

- DPC 2610.28527 is unrefereed. I re-checked that its post-shift main term and Lemma 6 are RH-free, not their algebra, which DPCAudit audited under RH.
- Fujii and the density theorem are used as Lamzouri quotes them.
- Lean: not done (mandatory per DISCOVERIES rule; out of this brief).

## Files

Riemann/Deep/UAudit/: REPORT.md, FULLPROOF.md, slit_identity.py, sparse_kernel.py, vmoments.py, vmoments_planted.py, logs/ (slit_identity*.json/log, sparse_kernel.json/log, vmoments.log, vmoments_planted.json/log).

## Resources

One compute process at a time, single-threaded, all under memcap.sh 290. Peaks: 50, 138, 50 and 50 MB; vmoments.py was killed by me after m = 6 (see its log). No worktrees, no commits, no Lean. Downloads deleted.

RH relevance: N8
