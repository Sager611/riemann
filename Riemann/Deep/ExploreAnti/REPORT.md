IDEA: under RH, replacing the complex Cauchy–Schwarz step of Conrey–Ghosh–Gonek by the sharp real (2D) second-moment small-ball bound doubles the effective mollifier length, κ* ≥ 1 − (1+2θ)⁻³ → 7/8 instead of 1 − (1+θ)⁻³ → 19/27, provided the holomorphic mean value Σ_γ ζ′(ρ)²M(ρ)² has its ratios-recipe main term (Conjecture H) (node N8)

# explore55 — anti-concentration through moments → N8 (simple zeros, RH branch)

## Transfer
- Tool (small-ball theory, "anti-concentration through moments"): for a random vector W ∈ ℝ^d with m = EW,
  P(W ≠ 0) ≥ mᵀ(E WWᵀ)⁻¹m, sharp given first and second moments (Paley–Zygmund at level 0 / Christoffel–Markov–Stieltjes
  with linear tests). For complex Y, the complex form P(Y ≠ 0) ≥ |EY|²/E|Y|² is the circular special case: it loses whenever
  Y is improper, E(Y − EY)² ≠ 0 (the signal-processing "widely linear" gain, Picinbono–Chevalier 1995).
- Node N8, RH branch: CGG (Proc. LMS 1998, GLH) and Bui–Heath-Brown (arXiv 1302.5018, RH only): κ* ≥ 19/27 from
  N*(T) ≥ |Σ_γ ζ′(ρ)B(ρ)|²/Σ_γ|ζ′(ρ)B(ρ)|², B(s) = Σ_{k≤y} μ(k)P(log(y/k)/log y)k⁻ˢ, y = T^θ, θ < ½, P(0) = 0, P(1) = 1.
- Why ζ′(ρ)B(ρ) is improper: under RH ζ′(ρ) = −i e^{−iϑ(γ)}Z′(γ) (checked to 5·10⁻⁸ at 8 zeros), so ζ′(ρ)² = −χ(ρ)|ζ′(ρ)|²,
  and the phase χ(ρ) = e^{−2iϑ(γ)} correlates with |ζ′(ρ)|² (zeros sit between Gram points): Σ_γ ζ′(ρ)² is of full size TL³.
  Equivalently, Re(ζ′(ρ)B(ρ)) = ½ζ′(ρ)ψ(ρ) with the two-piece mollifier ψ(s) = B(s) − χ(1−s)B(1−s) (Bui–Conrey–Young-type χ-twisted piece).

## Statement
- Lemma (proved, no RH): Y_γ = ζ′(ρ)B(ρ) vanishes at multiple zeros, so by Cauchy–Schwarz over simple zeros
  N*(T) ≥ (Σ_{0<γ≤T} Re Y_γ)² / Σ_{0<γ≤T} (Re Y_γ)², and Σ(Re Y)² = ½S₂ + ½Re S_B, S₂ = Σ|Y|² (RH), S_B = Σ_γ Y_γ² (holomorphic).
- Known (BHB, l. 240 of the TeX, there P(1) = 1; S₁ is linear, S₂ quadratic in P): S₁ = Σ Y ~ (TL²/2π)u, S₂ ~ (TL³/2π)A, L = log(T/2π),
  p = P(1): u = p/2 + θ∫P, A = u² + p²/12 + ∫P′²/(12θ) (= p²/3 + θp∫P + θ²(∫P)² + ∫P′²/(12θ)).
- Conjecture H (new, holomorphic): for 0 < θ < ½, S_B = Σ_{0<γ≤T} ζ′(ρ)²B(ρ)² = (TL³/2π)(u² − p²/12) + o(TL³)
  (bilinear: Σ_γ Y_PY_Q ~ (TL³/2π)(u_Pu_Q − p_Pp_Q/12)). Equivalently the pseudo-variance E(Y − EY)² = −p²L²/12.
- Theorem (conditional on RH and H): N*(T) ≥ (κ₂(θ) − ε)N(T), κ₂(θ) = u²/(u² + ∫P′²/(24θ)) maximised over P:
  κ₂(θ) = 1 − (1+2θ)⁻³, attained by P = [(½+θ)x − θx²/2]/(½ + θ/2) (θ = ½: P = (4x − x²)/3); θ → ½: κ* ≥ 7/8 = 0.875.
  Complex CGG with the same inputs: 1 − (1+θ)⁻³ (19/27 at θ = ½, BHB's P = −θx² + (1+θ)x). So κ₂(θ) = κ_CGG(2θ).
- Mechanism: Var Y = p²/12 + ∫P′²/(12θ) but Var Re Y = ∫P′²/(24θ), Var Im Y = p²/12 + ∫P′²/(24θ) (units L²): the p²/12 and half of
  the mollifier term sit in the direction orthogonal to the mean, which the real test ignores.
- Ceiling over real-part tests of the family: over all P with P(0) = 0 (any P(1)), equivalently all real-linear combinations of the
  Re Y_{P_i} (Σc_i Re Y_{P_i} = Re Y_{Σc_iP_i}), the best is exactly 1 − (1+2θ)⁻³ (Riesz representer of P ↦ P(1)/2 + θ∫P in
  ⟨P,Q⟩ = ∫P′Q′ is R = (½+θ)x − θx²/2, sup u²/∫P′² = ∫R′² = (3 + 6θ + 4θ²)/12; θ = ½: P = (4x − x²)/3, u = 7/9, ∫P′² = 28/27, κ₂ = 7/8):
  exact sympy algebra in [exact.py](exact.py), [logs/exact.log](logs/exact.log).
- Not proved (referee, Block 3): the ceiling over all real-linear combinations of Re/Im of several Y_P additionally needs
  Im Σ_γ Y_P·conj(Y_Q) = o(TL³); the diagonal S₂ and H do not determine these antisymmetric mixed covariances. Evidence only
  ([mixcov.py](mixcov.py), [logs/mixcov.log](logs/mixcov.log)): heuristically the recipe density for Σ_γ ζ′(ρ)ζ′(1−ρ)B_P(ρ)B_Q(1−ρ) is real
  for real shifts, so its main term is real and (Hermitian symmetry) equals the polarisation of S₂; data at L = 9.39 (zeros 50001–10⁵,
  P ∈ {x, x², x³, P_BHB, P₂}): |Im Σ Y_P conj(Y_Q)|/(NL²) ≤ 1.2·10⁻⁴ (Re ≈ 0.8), |Cov(Re Y_P, Im Y_Q)| ≤ 1.6·10⁻⁴ L², and adding the Im parts
  to the J = 3 real test moves it by ≤ 2·10⁻⁶ (θ = 0.35: 0.891332 → 0.891333; θ = 0.45: 0.923217 → 0.923219); at 10¹² (L = 24.48, 10⁴ zeros):
  ≤ 1.3·10⁻³, ≤ 2.1·10⁻³ L² and ≤ 4.5·10⁻⁵ (θ = 0.35: 0.835436 → 0.835481; θ = 0.45: 0.883755 → 0.883764); general-p algebra: [exact_p.py](exact_p.py), [logs/exact_p.log](logs/exact_p.log).
  Not covered either: other test families (Feng's Λ^{*k} pieces, ζ″-pieces, B(1−ρ)-twists; the last adds ≤ 0.004 numerically).
- Threshold: H for some θ_H > ¼ already beats 19/27 (κ₂(¼) = 19/27 exactly); H at θ_H → ½ gives 7/8. If θ could reach 1: 26/27.

## Why H should hold (ratios recipe, leading order)
- Σ_γ G(ρ) = (2πi)⁻¹∮(ζ′/ζ)G; right edge ≈ 0; left edge: ζ′/ζ(s) = χ′/χ(s) − ζ′/ζ(1−s). Write B(s) = (2πi)⁻¹∫y^z K(z)/ζ(s+z)dz,
  K(z) = Σ_j p_j j!/((log y)^j z^{j+1}), G = ζ(s+α)ζ(s+β)B(s)², and apply the CFZ ratios recipe to
  ζ(s+α)ζ(s+β)ζ(1−s+δ)/(ζ(s+z₁)ζ(s+z₂)ζ(1−s+γ)) (arithmetic factor 1 at leading order). The χ′/χ term and the identity term are
  sums of functions of α alone or β alone and die under ∂_α∂_β; the swaps α↔δ, β↔δ give
  D(α,β) = [β²g(α) − α²g(β)]/(αβ(β−α)), g(x) = e^{−xℓ}(p − xλ)², ℓ = log(t/2π), λ = θL∫P, and ∂_α∂_β D|₀ = p²ℓ³/6 + pℓ²λ + ℓλ².
  Hence S_B ~ (TL³/2π)(p²/6 + θp∫P + θ²(∫P)²) = (TL³/2π)(u² − p²/12), u = p/2 + θ∫P.
- Same recipe, checks: Σ_γ ζ′(ρ) ~ (T/4π)L² (classical) and BHB's S₁ ~ (TL²/2π)(p/2 + θ∫P) exactly (sympy, [exact.py](exact.py), general p: [exact_p.py](exact_p.py)).
- Proof route (not done): S_B is a first moment of a degree-2 object at the zeros; the swap main term comes from the χ²-stationary
  points t₀ = 2π√(n/h) with n/h a square; the errors need cancellation in Σ Λ-weighted sums with phase e(2√(n/h)), h ≤ y²,
  i.e. Bui–Heath-Brown-type generalised Vaughan work. The Landau–Gonek extension with χ (arXiv 2601.18025) is the nearest tool.

## Tests (single-threaded, under memcap 290 MB; incidents in [logs/MEMORY_NOTE.txt](logs/MEMORY_NOTE.txt): one uncapped 436 MB peak (fixed), four runs killed at the cap, one at 297 MiB)
Z′(γ) by Riemann–Siegel ([rs.py](rs.py)), validated against mpmath at 8 zeros (rel. err. 2.2·10⁻³ at γ₁, ≤ 1.2·10⁻⁵ beyond zero 1000,
[logs/rs_validate.log](logs/rs_validate.log)); at 10¹² phases reduced mod 2π in mpmath, |Z(γ)| ≤ 7.5·10⁻⁷ at all 10⁴ zeros ([moll_hi.py](moll_hi.py)).
Odlyzko's tables (zeros1, zeros6 tail, zeros3 near 10¹²). P = x unless stated; u, A, B in units L² (B rotated so EY > 0; |Im B| ≤ 0.006 where recorded).

| zeros | L | θ | CGG (complex) | 2D (real) | u | A | B | B − u² (pred −1/12 = −0.0833) |
|---|---|---|---|---|---|---|---|---|
| 1–100 | 3.63 | ¼ | 0.7535 | 0.9231 | 0.656 | 0.572 | 0.362 | −0.069 |
| 1–100 | 3.63 | ½ | 0.8917 | 0.9757 | 0.829 | 0.770 | 0.637 | −0.049 |
| 50001–10⁵ | 9.39 | ¼ | 0.6105 | 0.8177 | 0.6663 | 0.7272 | 0.3586 | −0.0854 |
| 50001–10⁵ | 9.39 | ½ | 0.7856 | 0.9275 | 0.8054 | 0.8258 | 0.5730 | −0.0757 |
| 1901053–2001052 | 12.10 | ¼ | 0.5795* | 0.7901* | 0.6668 | 0.7673 | 0.3582 | −0.0864 |
| 1901053–2001052 | 12.10 | ½ | 0.7599* | 0.9126* | 0.8013 | 0.8449 | 0.5622 | −0.0799 |
| 10¹²+1–10¹²+10⁴ | 24.48 | ¼ | 0.5299 | 0.7451 | 0.6398 | 0.7725 | 0.3263 | −0.0831 |
| 10¹²+1–10¹²+10⁴ | 24.48 | ½ | 0.7237 | 0.8924 | 0.7726 | 0.8248 | 0.5131 | −0.0839 |
| asymptotic (recipe) | ∞ | ¼ | 0.4839 | 0.7009 | 0.625 | 0.8073 | 0.3073 | −0.0833 |
| asymptotic (recipe) | ∞ | ½ | 0.6923 | 0.8710 | 0.75 | 0.8125 | 0.4792 | −0.0833 |
(* = recomputed from u, A, B for P = x; the optimal polynomials give, at L = 9.39 / 12.10: CGG(P_BHB) 0.8044 / 0.7782 and 2D(P₂) 0.9326 / 0.9183
at θ = ½; three-polynomial optima at L = 24.48, θ = ½: 0.7432 vs 0.8995.) Logs: [moll_v1](logs/moll_v1.log), [fit_z1](logs/fit_z1.log),
[fit_z6](logs/fit_z6.log), [hi_t025](logs/hi_t025.log), [hi_t050](logs/hi_t050.log).
- Reading: B − u² sits at −1/12 within 0.009 from L = 9.4 on and within 0.0006 at L = 24.5 (θ = ¼ and ½), while u and A still drift
  by O(1/L) toward the BHB constants; the finite-T bounds overshoot the asymptotics because Var Y is still below its limit.
- Unmollified ([logs/unmoll.log](logs/unmoll.log)): E ζ′(ρ)/L = 0.573, 0.554, 0.545, 0.520 (L = 6.0, 8.0, 9.4, 24.5; → ½); E ζ′(ρ)²/L² = 0.202, 0.174, 0.172,
  0.189 − 0.057i (the last noisy: 10⁴ zeros, heavy tails): positive, near the recipe's 1/6; pseudo-variance −0.125 (L = 9.4), −0.081 (L = 24.5).

## Independent model check (CUE / function-field limit, [cue.py](cue.py))
L(u) = Π(1 − a_j u), a_j eigenvalues of Haar U(N), zeros u_j = ā_j, "ζ′(ρ)" ↔ u L′(u), mollifier Σ_{k≤θN}P((K−k)/K)h_k(a)u^k (1/L = Σh_k u^k);
N plays L. 60,000–80,000 zeros per row ([logs/cue_v1.log](logs/cue_v1.log), [logs/cue_opt.log](logs/cue_opt.log)):
- P = x, θ = ½, N = 10…200: u = 0.739–0.753 (pred 0.75), A = 0.804–0.816 (0.8125), B − u² = −0.082 ± 0.002 (−1/12), CGG 0.679–0.696 (0.6923),
  2D 0.862–0.872 (0.8710); θ = ¼: u = 0.598–0.628 (0.625), B − u² = −0.079…−0.084, 2D 0.633–0.706 (0.7009).
- Optimal polynomials, N = 50 / 200: θ = ½: CGG(P_BHB) 0.7059 / 0.7050 (19/27 = 0.7037), 2D(P₂) 0.8752 / 0.8751 (7/8);
  θ = ¼: 0.4819 / 0.4954 (1 − (4/5)³ = 0.4880), 0.6959 / 0.7103 (19/27). The model has no arithmetic, so it checks the recipe algebra, not H for ζ.

## Control: Milinovich–Ng's J₋₁ bound (arXiv 1106.1160) gains nothing ([mncontrol.py](mncontrol.py), [logs/mncontrol.log](logs/mncontrol.log))
Their J₋₁ ≥ |Σ B(1−ρ)/ζ′(ρ)|²/Σ|B(ρ)|² is "half the conjectured value" because θ/(1+θ) < ½ for θ < 1, not because of impropriety:
adding the conjugate features moves the captured fraction of J₋₁ only 0.5719 → 0.5785 (zeros 10⁴–2·10⁴, θ = ½) and 0.5311 → 0.5363
(zeros 5·10⁴–10⁵); empirical J₋₁/((3/π³)ΔT) = 0.990, 0.997 (Gonek's conjecture). So the gain is specific to statistics with a full-size
pseudo-moment (ζ′(ρ)·mollifier), as the recipe predicts.

## Adversarial test (Davenport–Heilbronn, [dh_test.py](dh_test.py), [logs/dh_test.log](logs/dh_test.log))
DH (κ = 0.2840…, explore23's dhlib), 501 on-line zeros in [5, 640] (|Z_DH(γ)| ≤ 1.1·10⁻¹²), f′(ρ) = −i e^{−iϑ}Z′_DH, mollifier = truncated
Dirichlet inverse of DH's coefficients; ζ at the same height (370 zeros, L = 4.62):
- θ = ¼: DH CGG 0.7216 → 2D 0.8688, B − u² = −0.0720; ζ 0.7117 → 0.8960, −0.0719. θ = ½: DH 0.8424 → 0.9392, −0.0612; ζ 0.8687 → 0.9641, −0.0505.
- With ζ's Möbius mollifier DH degrades (θ = ½: 0.4161 → 0.5935, Var 1.13 vs 0.10).
- Verdict: the impropriety and the 2D gain are generic to functions with a real Z-function (functional equation, zeros on the line);
  they do not tell ζ from DH. The ζ-specific inputs are only the mollified main terms (Möbius = Euler product) and RH for S₂.
  This is an N8 (simple-zero) tool, not RH progress; DH violates RH, so the conditional statement is vacuous for it.

## Smallest case
First 100 zeros (Riemann–Siegel Z′, rel. err. ≤ 2.2·10⁻³; L = 3.63): CGG 0.7535 / 0.8917 → 2D 0.9231 / 0.9757 (θ = ¼ / ½); both trivially ≤ 1 (all simple),
B − u² = −0.069 / −0.049 (already negative, far from −1/12 at this L). No finite-T statement is claimed: all computed zeros are simple
(Platt–Trudgian); the content is the asymptotic main term, so certification here is the exact algebra (exact.py), not interval arithmetic.

## Further tests of H
- Bilinear H at 10¹² ([hi_cross.py](hi_cross.py), [logs/hi_cross.log](logs/hi_cross.log)): for P ∈ {x, x², x³, P_BHB, P₂}, all 15 pairs
  B_ij − u_iu_j lie in [−0.0865, −0.0828] (θ = 0.35) and [−0.0850, −0.0840] (θ = 0.45), vs −1/12 = −0.0833 (needed for the Riesz ceiling);
  2D(P₂) / CGG(P_BHB) = 0.8319 / 0.6388 (θ = 0.35) and 0.8819 / 0.7124 (θ = 0.45); asymptotic 0.7965 / 0.5936 and 0.8542 / 0.6720.
- P(1)-dependence (sharp test of −P_i(1)P_j(1)/12): for P = x − x² (P(1) = 0) the same data give B − u² = (−0.0838) − 2(−0.0835) + (−0.0844)
  = −0.0012 (θ = 0.35) and +0.0002 (θ = 0.45), predicted 0: the pseudo-variance is carried by P(1), as the recipe says.
- Falsification attempt of H beyond BHB's range ([falsH.py](falsH.py), [logs/falsH.log](logs/falsH.log), zeros 50001–10⁵, P = x):
  B − u² = −0.0757, −0.0726, −0.0712 at θ = ½, ¾, 1 (mollifier² up to length T²); drift 0.012 at L = 9.4, the size of the
  lower-order terms seen at θ = ½: no break found. 2D vs CGG there: 0.9275/0.7856, 0.9642/0.8658, 0.9798/0.9085.

## Comparison with the node line (TREE.md N8)
- Unconditional record (67.3639 %, TEN certificates) untouched: the method needs RH to evaluate S₂ (off the line Σ|Y|² is not a contour sum, and
  ψ = B − χ(1−s)B(1−s) is large at off-line zeros since |χ(1−ρ)| = (γ/2π)^{β−½}).
- RH branch: 19/27 = 0.7037 (CGG with GLH; Bui–Heath-Brown with RH alone) is the benchmark quoted in the N8 line; this gives
  7/8 = 0.875 conditional on RH + H, and > 19/27 from H at any θ_H > ¼. It is a reduction (new tool + one open mean value), not a theorem.

## Novelty (adversarial)
- arXiv API ([.scratch/axq.py]): "simple zeros"+mollifier / +"19/27" / +"real part" / +"Riemann hypothesis"; "discrete mean value";
  "Conrey, Ghosh and Gonek"; "negative moment"+derivative; "proportion of simple zeros"; pseudo-variance+zeta. Hits read: BHB 1302.5018
  (TeX: complex Cauchy only, no real part), Milinovich–Ng 1106.1160 (TeX), 0706.1763, 2311.13554 (Σζ(ρ+α)X(ρ)Y(1−ρ), one ζ),
  2106.03005 and 2305.14253 (Σζ^{(n)}(ρ), Shanks), 2509.07788/07792 (unmollified holomorphic products at zeros, conjectural),
  2601.18025 (Σχ(ρ)X^ρ), 1806.01959, 2109.15311, 2310.07360, 2501.14545: none uses the real (widely-linear) Cauchy–Schwarz or
  Σζ′(ρ)²B(ρ)². Monthly listings math.NT 2026-10/11 (lst.py, simple zero|derivative|mollif|discrete moment|Shanks): 9 titles, none.
  Zenodo '"simple zeros" zeta': unrelated top hits. Repo NOVELTY.md files and NOGO/TREE: nothing. CGG 1998 itself not read
  (not on arXiv; BHB §2 summarises it with the complex form). Verdict: not found; an expert check of the CGG/Bui literature is still needed.

## Risks and what is missing
- H is the whole content. The recipe and the CUE model agree exactly; ζ data support it to 6·10⁻⁴ at L = 24.5. A proof needs the χ²
  stationary-phase main term (n/h a square) and cancellation in Λ-weighted sums with phase e(2√(n/h)) for h ≤ T^{2θ}; Levinson-type
  experience (Conrey's θ > ½ needed Kloosterman sums) says ranges are where such proofs stop: θ_H > ¼ is the target.
- The gain is generic to real Z-functions (DH, CUE), so it is an N8 tool with no RH content of its own.
- Scope: every CGG-type simple-zero argument for an L-function with a functional equation inherits the gain under GRH, since
  L′(ρ)² = −ε χ_L(ρ)|L′(ρ)|² at critical zeros (Dirichlet, Dedekind (2310.07360), GL(2) families): same Lemma, same holomorphic input.
- Next brief (Orca, worker): prove H for θ < θ_H inside Bui–Heath-Brown's framework (contour, Gonek-type lemma for χ², generalised
  Vaughan), first rigorously for small θ, then push θ_H past ¼; falsify on the way at the 10¹² window (θ ∈ (¼, ½), P(1) ≠ 1, two-polynomial cross terms B₁₂ = u₁u₂ − 1/12).

## Suggested NOGO line
- explore55 (Orca, exploration: anti-concentration through moments — the sharp real second-moment small-ball bound P(W ≠ 0) ≥ mᵀ(EWWᵀ)⁻¹m, i.e. widely-linear Paley–Zygmund → N8) IDEA: Conrey–Ghosh–Gonek / Bui–Heath-Brown's complex Cauchy–Schwarz ignores that Y = ζ′(ρ)B(ρ) is improper (ζ′(ρ)² = −χ(ρ)|ζ′(ρ)|² under RH); the real test doubles the effective mollifier length, κ* ≥ 1 − (1+2θ)⁻³ → 7/8 instead of 1 − (1+θ)⁻³ → 19/27, under RH plus Conjecture H: Σ_γ ζ′(ρ)²B(ρ)² ~ (TL³/2π)(u² − 1/12) (ratios recipe; CUE model 0.8751 vs 7/8; ζ zeros: pseudo-variance −0.0839/−0.0831 at 10¹² vs −1/12); exact ceiling over real-part tests of the family (with Im parts it also needs Im Σ_γ Y_P·conj(Y_Q) = o(TL³): unproved, data ≤ 1.2·10⁻⁴); H at any θ > ¼ beats 19/27; generic to real Z-functions (DH −0.072), so no RH content; Milinovich–Ng's J₋₁ gains ≤ 0.007; H unproved [ExploreAnti](Riemann/Deep/ExploreAnti/REPORT.md)

## Files
[rs.py](rs.py) (Riemann–Siegel Z′), [moll.py](moll.py), [moll_hi.py](moll_hi.py) (10¹² window), [fit.py](fit.py), [unmoll.py](unmoll.py),
[falsH.py](falsH.py), [cue.py](cue.py), [mncontrol.py](mncontrol.py), [dh_test.py](dh_test.py) (uses ExploreConcentration/code/dhlib.py),
[exact.py](exact.py) (sympy), [NOTES.md](NOTES.md) (working notes), [logs/](logs/) (all outputs, memcap peaks, [MEMORY_NOTE.txt](logs/MEMORY_NOTE.txt)).
Data: Odlyzko's zeros1, zeros6 (tail), zeros3 (deleted from .scratch after the runs; re-download from www-users.cse.umn.edu/~odlyzko/zeta_tables).

RH relevance: N8 (RH-conditional proportion of simple zeros); the mechanism holds for DH and CUE, so it does not bear on RH itself.
