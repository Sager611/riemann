AUDIT PASSED: (a) [R2] Propositions 1 and 2 are correct as stated (Prop 1 is not load-bearing for RCApps5; Prop 2 is F₂ for D > 0), (b) the Humphries–Thorner inequality [R] uses is exactly HT §5 (bound_int) with Cor 4.4 (cor:MVT) at θ₂ = 0, which GRC for the dihedral, CM and λ = ¼ Maass families allows, and [R] uses it within its hypotheses, and (c) every [R4] step that RCApps5 §6 and the indefinite-forms corollary use is correct; the only defects are cosmetic constants (C(π_χ) ≤ 16D, not 12D, for odd narrow characters).

# RC5Audit follow-up: [R2] Props 1–2, the Humphries–Thorner input, [R4] for D > 0 and sectors (2026-11-12, orchestrator 31)

Long derivations: [DERIVATION2.md](DERIVATION2.md). Scripts and outputs: [checks/](checks/) (`fampos.gp`, `fampos*.out`). References: [R] = [RCApps NOTE](../RCApps/NOTE.md), [R2] = [RCApps2 NOTE](../RCApps2/NOTE.md), [R4] = [RCApps4 NOTE](../RCApps4/NOTE.md), [R5] = [RCApps5 NOTE](../RCApps5/NOTE.md), [D2] = [Density2 NOTE](../Density2/NOTE.md), [C] = [Combos NOTE](../Combos/NOTE.md); HT = Humphries–Thorner, arXiv 2103.05634 (e-print of 13 Mar 2024 = Compositio Math. 160 (2024) 1041–1072), TeX line numbers.

## Verdict table

| Step | Lines | Inputs | Verdict |
|---|---|---|---|
| (a) [R2] Prop 1: \|ℱ\| ≤ 2B(1+Λ), V = 𝔔^{1/2}, Y = (𝔔T)^{1/(3−2σ)} ≥ V, Z | R2 28–39 | (LS), [D2] Prop 6 step 4 | correct |
| (a) Prop 1 detector with mollifier length V | R2 38–39 | [C] Lemma 5 (refereed by RCAudit), [D2] Prop 6 step 3 | correct |
| (a) Prop 1 class II, first moment, Cauchy–Schwarz over all of ℱ | R2 40–47 | [D2] (2.1) at X = 1; (LS) and [C] Lemma 3 for M_V | correct: \|𝒳_II\| ≪ Y^{1/2−σ}𝔔^{1/2+ε}T^{1/2+ε} |
| (a) Prop 1 class I without k-th powers | R2 49–54 | (LS) at M = 2^{j+1} > V ≥ Λ/2, Sobolev | correct: \|𝒳_I\| ≪ (𝔔T)^εZ^{2−2σ+ε} |
| (a) Prop 1 conclusion, e(σ) = 2(1−σ)/(3−2σ) | R2 55–56 | — | correct: Y^{2−2σ} = Y^{1/2−σ}(𝔔T)^{1/2} = (𝔔T)^{e(σ)}; e(σ) < 2σ(1−σ) ⟺ (2σ−1)(1−σ) > 0 |
| (a) Theorem 1′(a), ring class characters | R2 58–61 | [C] Lemmas 2, 3, 4, 5 (all refereed) | correct; supersedes [C] Thm 1's exponent 2σ(1−σ) on [0.6, 1] |
| (a) Theorem 1′(b), (c) | R2 62–64 | [D2] Lemma 4 (M > 0), Lemma 5, (2.1) | arguments checked (Fejér + sector count; hyperbolic count), correct; (d) not checked; the per-D bounds themselves are not used by [R5] (their inputs [D2] Lemmas 4–5 are, in Theorems 2^∠, 2⁺) |
| (a) [R2] Prop 2 (D > 0 family density, F₂) | R2 112–122 | (I3), [D2] (2.1) for Ĝ⁺ | correct: AI(χ) cuspidal ⟺ χ² ≠ 1, Maass ¼ tempered, GRC; cosmetic: C(π_χ) ≤ 16D (odd type μ = (1,1)), not 12D |
| (a) Use in [R5] Thm 3 step 2′ | R5 159–163 | [R] Prop 1 at β > 9/10 and σ > 9/10 | correct |
| (a) Use in [R5] Thm 3 step 3 | R5 164 | [C] Thm 1(b) only | correct; [R2] Prop 1 not used |
| (a) Use in [R5] Thm 3⁺ | R5 209–215 | [R2] Prop 2 (F₂), h⁺ ≤ 5D^{1/2}log²(3D) per window | correct |
| (a) Use in [R5] Thm 3^∠ | R5 216–223 | [R] Prop 1′, V·I₂′ ≍ X⁵M⁶T⁵ | correct |
| (b) HT large sieve and the role of θ_n | HT 716–727, 779–800 | (LRS_2) for pairs in 𝒮 | θ_n enters only via the ramified Euler factors and (preconvex); θ₂ = 0 is legitimate when every π ∈ 𝒮 satisfies GRC |
| (b) HT cor:MVT (published Cor 4.4) | HT 842–869 | Gallagher, thm:main_theorem_1 | matches [R]'s V = Q^{2+ε}T^{3+ε}\|𝒮\| (n = 2, F = Q, π′ = 1) |
| (b) HT eqn:bound_int (§5) | HT 875–979 | detector, boxes of height 2(log Y)² | matches [R] (I3): Y enters only through log Y ≍ log X, κ = 0 for π′ = 1; HT then bound I₂ by preconvexity, [R] keeps I₂ |
| (b) [R] use: Prop 1, Prop 1′, [R2] Prop 2 | R 57–67, 167–169; R2 112–122 | GRC, Q ≥ C(π), log Y ≍ log V, T ≥ 2log²Y | within hypotheses; uniform (constants depend on n, n′, [F:Q], ε only); cosmetic: boxes reach \|v\| ≤ T + log²Y (use 2T) |
| (c) [R4] Lemma 1⁺ (Grössencharacters λ_k) | R4 31–36 | — | correct |
| (c) [R4] Lemma 2⁺(i), unfolding the units | R4 38–55 | Poisson in u, ŵ = ν̂² ≥ 0 | correct (re-derived Σ_kŵ(k/R)e(ku/R) = RΣ_jw(u−jR)) |
| (c) [R4] Lemma 2⁺(ii), h⁺R⁺ = √Δ L(1,χ)Π(1−(d_K/p)/p) | R4 43–44, 56–58 | Dirichlet's formula | correct |
| (c) [R5] Lemma P⁺ | R5 175–182 | Lemma 2⁺(i) for Δ_e | correct |
| (c) [R4] Theorem 1⁺ at w = 1: Lemma 5⁺, squares, Lemma 6⁺, §3.5 | R4 67–94; R5 184–191 | [R3] §3.6–§3.8 | correct (Lemma 6⁺ re-derived: ∫f = √M·v·K⁺(v), W smooth via w = w₀(x²), hypothesis X ≥ 128c√MX^δ) |
| (c) [R4] Theorem 2⁺ | R4 96–103; R5 206–208 | [D2] Lemma 5, (2.1) for Ĝ⁺ | correct |
| (c) [R4] Theorem 3⁺, H ≤ x | R4 105–118; R5 209–215 | [R2] Thm 4 dictionary | correct; [R5] keeps H ∈ [xD^{−κ}, x] (the [R4] statement says only H ≥ xD^{−κ}) |
| (c) [R2] Corollary 5 (indefinite forms, a² − nb²) | R2 147–156 | narrow class ↔ proper equivalence | correct: N(ua + v(−b+√D)/2) = aQ(u, −v) |
| (c) [R4] Lemma 2^∠; [R5] Lemma P^∠, Thm 1^∠ | R4 123–144; R5 193–204 | — | correct |
| (c) [R4] Theorem 2^∠ | R4 146–152 | [D2] Lemma 4 (M > 0) for class II | correct (Lemma 4's argument checked) |
| (c) [R4] / [R5] Theorem 3^∠ | R4 154–170; R5 216–223 | [R] Prop 1′ | correct: 4Aδ + 2κ + ε₁ + ε₂ ≤ 6δ |

## Findings (none breaks a proof)

1. **[R2] Proposition 1 holds**, so the project's best per-D density for ring class characters is e(σ) = 2(1−σ)/(3−2σ) on 0.6 ≤ σ ≤ 1 (Theorem 1′(a)), better than Combos Theorem 1's 2σ(1−σ) by the factor σ(3−2σ) (at most 9/8, at σ = ¾). Its inputs are exactly the refereed Combos Lemmas 2–5. RCApps5 does not use it.
2. **Humphries–Thorner.** [R] quotes HT's intermediate inequality (§5 bound_int with Cor 4.4), not their Theorem 1.2. The quote is correct; θ₂ = 0 is justified because θ_n enters only through (LRS_2) for pairs in 𝒮, and every family used satisfies GRC. ([R] adds that Pascadi–Thorner 2508.14888 remove the condition; not re-checked here, and not needed.)
3. **Cosmetic.** C(π_χ) ≤ 16D for odd narrow characters ([R2] says 12D); HT's box argument integrates over |v| ≤ T + (log Y)², not T; [D2] (2.1) is stated on Re s = ½ and used at ½+δ (routine, as for D < 0). All only change constants.
4. **Not re-refereed:** [R4] §3–§4 (a_E signs, sharpening; unused), [R2] Theorem 1′(d), [D2] Theorem C.

## Numerical test, D > 0 (`fampos.gp`, PARI 2.19)

D > 0 analogue of `fam.c`. The family is the fundamental D ∈ (X/2, X]. Narrow class groups come from bnrinit(bnf, [1, [1,1]]), and each class from bnrisprincipal of one prime above each split or ramified p ≥ 19, with [𝔭̄] = −[𝔭]. The detector b is that of `fam.c` (σ = 5/6, V = √X/2, admissible n). LHS = Σ_D h⁺Σ_C|S_C|² = Σ_DΣ_{χ∈Ĝ⁺}|Σbλ_χ|². R_fam := LHS/(Σ_D h⁺·Σ|b|²) normalises by the actual family size; R_X := LHS/(X^{3/2}Σ|b|²) is the Theorem 1⁺ normalisation.

| X | Σ_D h⁺ (Σ_D h⁺R⁺) | N | R_fam killed | R_fam pole | R_X killed | R_X pole | χ = 1 share, killed |
|---|---|---|---|---|---|---|---|
| 1000 | 454 (5420) | 128 / 256 / 512 | 1.13 / 1.44 / 1.41 | 15.3 / 24.1 / 39.4 | 0.016 / 0.021 / 0.020 | 0.22 / 0.35 / 0.57 | 29–50 % |
| 2000 | 1033 (16109) | 128 / 512 / 1024 | 1.34 / 1.52 / 1.26 | 15.1 / 41.4 / 69.3 | 0.015 / 0.018 / 0.015 | 0.17 / 0.48 / 0.80 | 29–51 % |
| 4000 | 2268 (46207) | 256 / 1024 / 2048 | 1.15 / 1.77 / 1.33 | 22.6 / 60.5 / 110.9 | 0.010 / 0.016 / 0.012 | 0.20 / 0.54 / 0.99 | 23–56 % |
| 8000 | 5290 (131782) | 256 / 1024 / 4096 | 1.37 / 2.10 / 1.23 | 19.7 / 48.2 / 179.3 | 0.010 / 0.016 / 0.009 | 0.15 / 0.36 / 1.33 | 32–59 % |

- Theorem 1⁺'s normalisation: the pole-killed R_X stays at 0.009–0.021 for N from 2.9√X to X/2 (at N = X/8: 0.016, 0.017, 0.016, 0.016 for X = 1000 … 8000; at N = X/2 it falls, 0.020 → 0.009). The pole-present R_X grows almost linearly in N (0.15 → 1.33, ×9 for N ×16, at X = 8000) and with X at fixed N/X = ½ (0.57, 0.80, 0.99, 1.33, slope ≈ 0.41).
- Family-size normalisation: the pole-killed sum is 1.1–2.1 times Σh⁺·Σ|b|² for N ≥ 2.9√X, the family size of the narrow class characters, far below Theorem 1⁺'s X^{3/2} (Σ_D h⁺ = 0.45–0.66X here). Its slow rise at N = X/8 (1.13 → 2.10) is the trivial character, about half of the sum. With a mollifier of length V = √X/2 that term can be as large as X·V·Σ|b|² ≍ X^{3/2}Σ|b|² (heuristic scale), which Theorem 1⁺ allows. At N ≤ 2√X it is 3–8 (χ = 1 carries 75–90 %, as for D < 0).
- Memory 24 MB, 241 s for X = 8000; gp under a perl fork-and-kill timer.

RH relevance: N7, N11. These inputs feed the "every A > 1" least-prime line (records only, no RH node moves). Together with RC5Audit's first report, the analytic chain of RCApps5 Theorems 3, 3⁺ and 3^∠ is now refereed: the family sieve, the density, F₂, the HT input, and the D > 0 and sector reductions, down to the refereed Combos Lemmas and the published HT inequality. Only margin-checked: [R] §2's exceptional-set and Siegel bookkeeping ((I4)–(I6), E₂, §2.4(e), §2.5).
