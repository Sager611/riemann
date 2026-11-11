AUDIT PASSED: Lamzouri arXiv 2609.02882v2 (latest): Prop 2.1, Lemma 3.1 (= BGST Lemma 5, same hypotheses) and Lemma 3.2 re-derived line by line, no gap; Thm 1.1 follows; constants recomputed in Arb; no counterexample in 27,000 random multisets and 288 adversarial local searches; one harmless misquote of Lemma 3.2 in 2610.26841.

Auditor: lamz2audit (Orca), 11 Nov 2026. Sources: TeX e-prints 2609.02882 v1 and v2, BGST 2306.04799 v1, 2610.26841 v2, 2610.08965 v1, 2608.13637 v2, 1002.4127, six citing papers; Lean repo AxiomMath/ZetaZerosV2 @ 1d7edeb. Line numbers are TeX lines of 2609.02882 v2 unless stated.

## 1. Statements used

| # | Statement (TeX) | Hypotheses | Used by |
|---|---|---|---|
| S1 | Prop 2.1 (Pro:FourierHilbert, l.303–323). For a finite, nonempty, conjugation-invariant multiset Z and K = FT(η²): #(simple real) ≥ 2\|Z\| − ΣK(z−s)² (l.306); #distinct ≥ (3/2)\|Z\| − (1/2)ΣK² (l.311); #simple + #real (with mult.) ≥ 3\|Z\| − ΣK² (l.316); if ΣK² ≤ A\|Z\|, 1 ≤ A < 2: #(simple or real, with mult.) ≥ (5+2√2−2A)/(3+2√2)·\|Z\| (l.320) | λ > 0; η ∈ L²(ℝ) real, even, supp ⊂ (−λ, λ), FT(η²)(0) = 1; no smoothness | Thm 1.1; 2610.26841 re-proves its variants (l.284–313); IsoZero through them and DERIVATION §3.5 (κ = 2 − C(f)); 2609.07918, 2610.22927, Zenodo 23002514 |
| S2 | Lemma 3.1 (Lem:BGST, l.529–535) = BGST Lemma 5 (2306.04799 l.349–358): Σ_{0<γ,γ′≤T} FT(f)(i(ρ−ρ′)log T/2π)·w(ρ−ρ′) = (T/2π) log T·(f(0) + 2∫₀¹ αf(α)dα + O_f(1/√log T)), w(u) = 4/(4−u²), all zeros with multiplicity | f ∈ L¹ real, even, supp ⊂ [−1, 1], Lipschitz at 0. Quoted with exactly BGST's hypotheses; only O → O_f differs | Lemma 3.2; 2610.08965 prop:energy (l.780–815) directly; SmallGap (2611.14486) directly |
| S3 | Lemma 3.2 (Lem:ConstructionKernel, l.583–626): for every ε > 0 there exists a real even η_ε ∈ C_c^∞((−½, ½)) with K_ε(0) = 1 and Σ_{ρ,ρ′} K_ε(i(ρ−ρ′)log T/2π)² = (C_{η_ε} + O_ε(1/√log T))(T/2π) log T, \|C_{η_ε} − C_MT\| < ε, C_MT = ½ + cot(1/√2)/√2 | none beyond Lemma 3.1 | Thm 1.1; quoted by 2610.26841 l.670–701 (see F1); IsoZero DERIVATION §1.2 (l.28–29); ARTZ 2610.22927 (l.200–253); Dirichlet analogue 2610.26841 l.976–997 ("GaPa + the proof of Lemma 3.2") → our linear-combinations pages |
| S3′ | Form G, what the proof (l.655–853) actually shows: for every real even η ∈ C_c^∞((−½, ½)) with ∫η² = 1, ΣK² = (C(η²) + O_η(1/√log T))(T/2π) log T, C(f) = Q(0) + 2∫₀¹αQ = ∫f² + ∬\|u−v\|f(u)f(v), Q = f∗f | as S3 | IsoZero (f_0 cutoffs for the headline, cos(wx) cutoffs in column H+w); 2610.08965 proves the same form itself |
| S4 | Thm 1.1 (Thm:Main, l.211–236): liminf N₀ˢ/N ≥ C₀ = 0.6725007, N_d/N ≥ C₁ = 0.8362504, N_{s∪0}/N ≥ C₂ = 0.8876200, (N_s + N₀)/2N ≥ C₁ | unconditional | the 2 − C_MT = 0.672501 of IsoZero and TREE N8 |
| — | Not used downstream: Rem 2.2 (l.325–332, ΣK² ≥ \|Z\|), Rem 3.3 (l.856–860, optimality C_η ≥ C_MT via CCLM Cor. 14), the comparison with AF Thm 5.7 (l.627–654) | | |

"K1" in the TREE N8 line is not a result of 2609.02882: in 2610.08965 it is the bibliography key of Knausgård's first paper, arXiv 2609.33043 (bib l.2438–2440), cited with BGST for the energy (prop:energy, "[BGST, Lemma 5], [K1, Proposition 5]"). 2610.08965 uses nothing of 2609.02882 logically: it re-proves the energy formula and the operator splitting in-line (l.780–830; "in [Lam, Section 3] the weight is removed by a modification of the test function as well").

## 2. Proof check

P1. Prop 2.1 (proof l.333–514), every line re-derived:
- l.341 K(z − s̄) = ⟨f_z, f_s⟩ (η real); l.345 ‖f_x‖ = 1 for real x; l.349 ‖g_z‖² − ‖h_z‖² = Re⟨f_z, f_z̄⟩ = K(0) = 1.
- l.358 ΣK(z−s)² = ΣK(z−s̄)² (conjugation invariance with multiplicities) = ‖F‖²; this is a square, not a modulus, and the sum is real.
- l.369 F = Σm_x f_x⊗f_x + 2Σm_z(g_z⊗g_z − h_z⊗h_z).
- l.381 linear independence: an exponential polynomial vanishing on a set of positive measure vanishes; (g_z, h_z) is an invertible change of (f_z, f_z̄); so D_U = r+k, D_V = n+r+k, D_W = n+r+2k.
- l.384–399 Φ̄(u) = Φ(−u) holds for f_x (x real), g_z and h_z; all inner products are real; Gram–Schmidt has real coefficients; Ψ_j is orthonormal; Bessel's inequality applies.
- l.403–409 α_j is real, and each ⟨·, ψ_j⟩² ≥ 0.
- l.422 Σ_U α ≥ Σ_{ℓ>n} m_x + 2Σm_z ≥ 2(r+k); l.430 Σ_{V∖U} α ≤ n; l.435 α_j ≤ 0 on W∖V; l.439 Σ_W α = \|Z\|.
- (2.1): from a²+4 ≥ 4a, a²+1 ≥ 2a and α ≤ 0. (2.2): the proof gives even n+r+k on the left. (2.3): the coefficient check Σ_U − Σ_{V∖U} − 2r − 2p − 4q + n ≥ 0 holds.
- (2.4): t = 2+√2 is the root of t² − 4t + 2 = 0 with b₂ = t²/4 > 0. Then b₁ = 5/2+√2−A, b₂ = 3/2+√2, b₃ = 1/2+√2, and the q-coefficient 4b₂ − t² = 0. Result (5+2√2−2A)/(3+2√2).

P2. Lemma 3.2 (proof l.655–853), re-derived:
- f₀ = cos(√2x)/(√2 sin(1/√2)) on [−½, ½] is positive (√2/2 < π/2) and has ∫f₀ = 1.
- η_δ = ψ_δ√f₀/√A_δ ∈ C_c^∞((−½, ½)), since √f₀ is smooth on supp ψ_δ ⋐ (−½, ½).
- Q_δ = f_δ∗f_δ ∈ C_c^∞((−1, 1)), so Q_δ and Q_δ″ satisfy Lemma 3.1.
- r_{δ,T} = Q_δ − Q_δ″/(4 log²T) (l.742) has FT = (1 + π²z²/log²T)K², which equals K²/w at z = i(ρ−ρ′)log T/2π (l.777–790).
- The Q″ term is O_δ(N/log²T).
- C_δ → C(f₀) by L¹ ∩ L² convergence (\|u+v\| ≤ 1 on the support). C(f₀) = C_MT (Euler–Lagrange: 2f₀ + 2∫\|u−v\|f₀ is constant on [−½, ½], so C(f₀) = ½ + f₀(½); recomputed three ways, §3).

P3. Thm 1.1 (proof l.862–969):
- Z_T = {i(ρ−½)log T/2π}. Conjugation corresponds to ρ ↦ 1−ρ̄, which has the same multiplicity by the functional equation. An element is real iff β = ½, z − s = i(ρ−ρ′)log T/2π, and N(T) ~ (T/2π) log T.
- One routine step is not written: for (2.4) take A = C_{η_ε} + ε for T ≥ T₀(ε). A ≥ 1 holds by Rem 2.2, and A < 2 since C_MT < 2 − 2ε.

P4. BGST, the load-bearing analytic input:
- Lemma 5 follows from Thm 1 as printed (BGST l.360–371): split at loglogT/logT; the tail is ≤ ∫\|r\|/log T.
- Proof of Thm 1, scanned (BGST l.186–330):
  - Lemma 1: the σ ↦ 1−σ subtraction algebra is correct.
  - Lemma 3 (F(x,T) = (2/π)∫\|Σ…\|²dt, hence F ≥ 0 with off-line zeros): ∫dt/((1+t²)(1+(t+a)²)) = 2π/(4+a²) for \|Im a\| < 1 (residues). The shift by iδ crosses no pole, since \|δ\|, \|δ′\| < ½, and the re-indexing ρ′ ↦ 1−ρ̄′ is valid.
  - Lemma 4 (truncation at Z = T log²T, x^{δ+δ′} ≤ x^{1−2η}): consistent.
  - Assembly with the Vinogradov–Korobov region and Goldston–Montgomery's R(x,T): F(α) = T^{−2α}(log T + O(1)) + α + O(1/√log T) follows.
  - Not re-derived: GM87's unconditional R(x,T) and the VK zero-free region (classical, refereed).

P5. Conrey:
- 2609.02882 does not use Conrey's mean value; Con89 appears only in the introduction.
- For IsoZero, Conrey's Theorem 2 as quoted by Bui–Conrey–Young (1002.4127 TeX l.335–344) carries these hypotheses: θ₁ < 4/7, Q(0) = 1, Q(x) + Q(1−x) = const, P₁(0) = 0, P₁(1) = 1, σ₀ = ½ − R/L, R bounded positive. IsoZero DERIVATION §3.1 (l.80–83) quotes it with all of them (Q ≡ 1, c₁ formula identical). Conrey's Crelle paper itself was not read (as in R395).

Flags (none affects a theorem):
- F1 (2610.26841, not 2609.02882). l.622–701 fixes δ = ε, takes any admissible cutoff ψ_ε and quotes Lemma 3.2 as \|C_{η_ε} − C_MT\| < ε. That is false for steep cutoffs.
  - Indicator limit: D(δ) = C − C_MT = −δ + (cot(√2(½−δ)) − cot(1/√2))/√2, and D(δ) − δ has derivative csc²(√2(½−δ)) − 2 ≥ csc²(1/√2) − 2 = 0.3695 > 0. So D(δ) > δ for every 0 < δ < ½, with D/δ → 1.3695.
  - Smooth ψ with transition width ≤ δ/2 also violates it: D/δ = 1.05–1.78 at δ = 0.003–0.09 (cutoff.log).
  - Harmless: 2610.26841 uses only C_{η_ε} = C_MT + O(ε) (l.816). That holds uniformly in ψ.
    - With g = f_ε − f₀ (∫g = 0), the first variation of C at f₀ vanishes, so C − C_MT = ∫g² + ∬\|u−v\|gg.
    - Upper bound: ∬\|u−v\|gg ≤ 0. On the boundary layer \|g\| ≤ f₀ ≤ f₀(½−ε), and on the core g = O(ε). Hence C − C_MT ≤ 2ε·f₀(½−ε)² + O(ε²) = cot²(1/√2)·ε + O(ε²) = 1.3695ε + O(ε²). The indicator is the extremal case.
    - Lower bound: C − C_MT ≥ −(∫\|g\|)² = −O(ε²).
  - 2609.02882's own statement is existential and correct.
- F2 (ours). The TREE N8 line reads "Lamzouri 2609.02882, its K1"; K1 is Knausgård 2609.33043.
- F3. v1 printed C_MT = 1.3274992963206… (rounded); v2 has …205… (truncated, correct). C_MT = 1.32749929632058835.
- F4. The (2.4) hypothesis 1 ≤ A < 2 is stronger than needed: b₁ > 0 needs only A < 5/2 + √2.
- F5. The step noted in P3 is unwritten but routine.

## 3. Falsification tests (scripts here, logs in logs/; all single-threaded, < 60 MB, < 1 min each)

| Test | Result |
|---|---|
| consts.py (Arb 256 bits) | C_MT = 1.32749929632058835426, C₀ = 0.672500703679411646, C₁ = 0.836250351839705823, C₂ = 0.887620008173354340, (2+C₀)/3 = 0.890833567893; the claims 67.25%, 83.62%, 88.76% and (from C₂) 1 − C₂ < 0.1124 hold. C(f₀) = C_MT by three routes (Euler–Lagrange closed form; Arb integral; mpmath Q₀(0) + 2∫αQ₀ with Q₀ = f₀∗f₀ by double quadrature): agree to 1e−22. Indicator-cutoff D(δ) > δ certified at δ = 0.1, …, 1e−6 |
| cutoff.py | the paper's gentle cutoffs: C_δ − C_MT = 0.78–0.90δ > 0, converging to C_MT from above (consistent with CCLM); steep cutoffs violate F1's quote |
| fourier.py | FT(f″) = (2πiz)²FT(f) at random complex z, and FT(r_T)·w = K² at z = i(ρ−ρ′)log T/2π with \|Re(ρ−ρ′)\| < 1 and T up to 10³⁰: relative errors ≤ 2e−31 |
| propfals.py | Prop 2.1 (four inequalities) and Rem 2.2 on 27,000 random conjugation-invariant multisets: 6 kernels (Montgomery–Taylor f₀, cos 2.5x, box and Fejér with λ = ½, 2, 3); spaced, clustered and wide positions; multiplicities 1–3; Im from 1e−9 to 3. Min slack/\|Z\| ≥ −4.4e−16 (float noise at exact equality cases); ΣK² real to 2e−16. Adversarial local search on (2.1), 24 starts × 1500 steps per kernel: min ΣK² − (n+4k) = 1.2e−8 … 6.6e−4 > 0, reached as an off-line pair collapses onto the line (the double-zero equality case). Arb recheck of every tightest case: slacks ≥ 0, or balls containing 0 at exact equality |
| propfals2.py | 16 starts × 9 settings, pairs held at Im z ≥ y₀: min slack of (2.1) = 0.0129, 0.224, 2.50 at y₀ = 0.02, 0.1, 0.3; of (2.2) ≥ 1.005; of (2.3) = 0.0096, 0.167, 2.54; all Arb-confirmed |
| earlier, LamzAudit N3 | BGST/Montgomery main term on Odlyzko's 2·10⁶ zeros: consistent within the predicted 1/log T drift |

No violation anywhere.

## 4. Literature

Versions:
- 2609.02882: v1 (2 Sep 2026) and v2 (8 Sep 2026), the latest. v1→v2 adds (2.3), (2.4), Rem 2.2 and Appendix A; Lemmas 3.1 and 3.2 are unchanged except O → O_f and one digit of C_MT; the proof of Lemma 3.2 lost only commented-out text. No erratum or withdrawal.
- BGST 2306.04799: v1 only (7 Jun 2023); published Acta Arith. 214 (2024) 357–376; the published text was not read.

Alpöge–Furman 2608.13637:
- v2 (19 Aug 2026) is the latest. Proof found by Claude, verified by the authors, formalized in the zeta23 Lean library.
- Its energy input is Montgomery's prime-side second moment (Thm thm:HS, the "Theorem 5.7" of l.641), an independent derivation of what Lemma 3.2 gives (window ψ, dyadic range).
- No reported error. 2610.08965 (l.1968–1990) notes that AF §7.2's period-256 "optimal" law (p₀ = 0.68182868746, not public) is not the least share. That concerns AF's barrier section, not its theorem.

2610.08965 (v1 only):
- Re-proves the energy formula from BGST Lemma 5 with the same R − R″/4L² trick (prop:energy).
- This is a re-derivation, not a correction. No paper corrects Lemma 3.2.

Citing works (OpenAlex, 11 Nov 2026): 14 for 2609.02882 and 25 for AF; TeX contexts scanned, none reports a flaw. Re-uses:
- 2609.07918 re-proves Lemma 3.2 in short intervals with a local pair correlation.
- 2610.22927 uses it for fixed δ.
- 2610.23467 recalls the weight removal.
- 2610.24612 is a human exposition of AF.

Zenodo:
- Fabbian 23027645 sharpens (2.1) to ΣK² ≥ x + (\|Z\| − x)²/N_m. This is consistent and follows in three lines from Lamzouri's ranges (re-derived here): with s = Σ_U α ≥ M = \|Z\| − x and D_U = N_m, Cauchy–Schwarz gives Σ_U α² ≥ s²/N_m, and the other two ranges give ≥ 2(\|Z\| − s) − x. The sum s²/N_m + x + 2M − 2s increases for s ≥ N_m, and M ≥ 2N_m, so its minimum is at s = M: x + M²/N_m ≥ 2\|Z\| − x.
- Fabbian 23198632 (barrier: no certificate in its class beats C₂) is consistent.
- Zeraoulia 22286993 and Eskelinen 23002514 build on Prop 2.1.
- No flaw reported.

Formal proof, AxiomMath/ZetaZerosV2, commit 1d7edeb (3 Oct 2026, "make main theorem unconditional"):
- Challenge/Basic.lean states the four parts of Prop 2.1 and all four parts of Thm 1.1 with no hypotheses: Mathlib riemannZeta, open strip, multiplicity = analyticOrderNatAt, counts as in the paper. The statements are faithful and not vacuous: a zero count of 0 would make them false.
- kernelConstruction (ZetaZeros/Zeta/Kernel.lean:996) is Lemma 3.2 in its existential form, with o(1) in place of O_ε(1/√log T).
- pairCorrelationFormula (PairCorrelation/FromLandau.lean:52) proves BGST Lemma 5 for test functions vanishing on 1 − ε ≤ \|x\| and Lipschitz at 0 globally (enough for Q_δ and Q_δ″). riemannVonMangoldt is also proved.
- Checked here:
  - The GitHub Actions build of 1d7edeb succeeded.
  - No sorry, axiom, native_decide, implemented_by, extern or unsafe in the 100 Lean files outside the challenge stubs.
  - The 12 directly imported PrimeNumberTheoremAnd modules (@c19697c) have no live sorry.
- Not checked: the comparator run (permitted axioms propext, Quot.sound, Classical.choice) is claimed in the README and was not re-run (no Lean builds here).

## 5. Dependency map

| Our entry | Input from 2609.02882 | Fails if Lemma 3.2 failed | Fails if Prop 2.1 failed | Fails if BGST Lemma 5 failed |
|---|---|---|---|---|
| discoveries/isolated-zeros.md, RECORDS l.17, TREE N8 "isolated zeros ≥ 0.6724" | Lemma 3.2, form G (via 2610.26841 l.670–701) | yes | no (it uses 2610.26841's own Pro:IsolatedSMALL, audited in IsoZero/AUDIT.md) | yes |
| linear-combinations-line.md (0.3449), RECORDS l.21 (34.71%, 2.06%), TREE N8/N12 items | Dirichlet analogue (Garunkštis–Paliulionytė + the proof of Lemma 3.2) | yes (the transferred argument) | no | its Dirichlet analogue (GaPa) |
| 67.353% record (2610.08965), ZetaZeros audit, N8Fam, distinct 0.83900, simple 0.67800 | none (re-proved in-line) | no | no | yes |
| SmallGap λ ≥ 0.7851 (2611.14486) | none (Lemma BGST only, l.322–329, per LamzAudit/REPORT.md) | no | no | yes |
| GapLambda, ConjFals, GL(2) 7.005%, Selberg N ≥ 4 | none | no | no | no |

What survives each failure:
- Without Lemma 3.2, which cannot fail independently (it is two applications of Lemma 3.1 plus integration by parts, and Axiom's formal proof builds, comparator run claimed): 2610.08965 and AF.
- Without BGST: AF's prime-side route, which carries the same analytic content derived independently.
- Without both: PRZZ 0.417293 (on the line) and 0.407511 (simple, on the line), and Ki–Lee 70% distinct. BGST's own 61.7% (Thm 2) needs a thin-box hypothesis.

Suggested updates (for the orchestrator):
- TREE N8: "its K1" → "Knausgård 2609.33043 (K1)".
- isolated-zeros.md caveat: "Lemma 3.2 refereed (LamzAudit2): sound; machine-checked in AxiomMath/ZetaZerosV2 (comparator claimed)".
- IsoZero/AUDIT.md: add F1.

## Review

Two verdict-only reviewer runs (inline evidence, logs/review.log):
- R1, on F1: (a)–(c) SOUND. (d) FLAW as written: the uniform constant was quoted with max f₀ instead of the boundary value. It is corrected above to cot²(1/√2) = 1.3695, with the lower bound −O(ε²) added; harmlessness confirmed.
- R2, on the chain P1–P3: SOUND conditional on BGST Lemma 5. Its two literal flags were slips in my inline summary, not in the paper: the reality property is stated only for f_x, g_z and h_z (l.384), and the sums are over 0 < γ, γ′ ≤ T (l.529–535, 583–626).

## 6. Referee verdict (plain words)

Lamzouri's argument is correct:
- The key proposition is pure linear algebra in a Hilbert space and checks line by line.
- Every inequality survived randomized and adversarial float tests, with ball-arithmetic rechecks of the tightest cases.
- The only analytic input is the refereed BGST pair-correlation theorem, quoted with its exact hypotheses. Its proof is a faithful unconditional transcription of Montgomery's, resting on a positivity identity I re-derived.
- The step that removes the weight, Lemma 3.2, is elementary and right, and the constant 2 − C_MT = 0.6725007 is right.
- Axiom's Lean library claims a complete hypothesis-free formal proof of the whole theorem, including the pair-correlation input. Its build passes and its sources contain no sorry or axiom, though I did not run their axiom checker.

The only defect found is elsewhere: Lamzouri's later isolated-zeros paper restates Lemma 3.2 with a quantitative "< ε" that fails for steep cutoffs. Only an O(ε) bound is used there, so nothing breaks. Our isolated-zero record therefore rests on a sound input.

Residual trust:
- BGST's classical inputs: GM87's mean value and the VK zero-free region.
- Conrey's Crelle theorem as quoted by BCY (IsoZero only).
- Axiom's comparator run.

## Files
Riemann/Deep/LamzAudit2/: REPORT.md (this), consts.py, cutoff.py, fourier.py, propfals.py, propfals2.py, logs/{consts,cutoff,fourier,propfals,propfals2,review}.log. Run each with uv run --no-project --with python-flint --with mpmath --with numpy python3 <script>.
