AUDIT PASSED: RCApps5 NOTE §2–5 are correct as stated (Lemma P; Corollary P; Theorem 1, the family large sieve at size X^{3/2} for every N ≤ X^{1−4δ}, for the pole-killed detector sequences only; Theorem 2, Σ_D N_D(σ) ≪ X^{2−σ+ε} for ½ < σ ≤ 9/10; Theorem 3, every A > 1), and no counterexample appears for X ≤ 64000, N ≤ X.

# RC5Audit: referee of RCApps5 (Lemma P, family large sieve, density on average over D, least primes for A > 1), 2026-11-12

Brief: `.scratch/orca/rc5audit-brief.md`. Long derivations: [DERIVATION.md](DERIVATION.md). Scripts and outputs: [checks/](checks/). Line numbers refer to [RCApps5 NOTE](../RCApps5/NOTE.md); [R3] = [RCApps3 NOTE](../RCApps3/NOTE.md), [R] = [RCApps NOTE](../RCApps/NOTE.md).

## Verdict table

| Step | NOTE lines | Inputs | Verdict |
|---|---|---|---|
| §1 diagnosis: F₁ = (13−8σ)/(7−2σ) and F₂ = 10(1−σ)/(3−2σ) cross at σ_c = 5 − √69/2 | 22–46 | [R3] Thm 2, [R] Prop 1 | correct (σ_c = 0.846688, 6/(7−2σ_c) = 1.130662) |
| Lemma P (a): only equal f-parts survive; (c, f) = 1 | 63–66 | p ∥ M, p² \| Δ ⇒ no γ | correct |
| Lemma P (b): γ ↦ γ/e, r_Δ(e²M′) = r_{Δ_e}(M′) | 67–70 | Δ_e ≡ Δ (8) | correct |
| Lemma P (c): r_{Δ_e}(M) = (w_e/2)[1](T_M), B_e = (w_e/2)Σ_C\|S_C\|² | 71–77 | [R3] Lemma 1, P_{K,Z}(f_e) | correct; re-checked with own code on 231 non-maximal orders |
| Lemma P statement, "the term e = 1 is Σ_C\|S_C\|²" | 61 | — | cosmetic: needs \|Δ\| > 4; harmless (Φ = 0 at Δ = −3, −4 for X ≥ 32) |
| Corollary P: dual identity, LHS ≤ CX^{1/2}log X·𝒬′(b) | 79–89 | h ≪ X^{1/2}log X, Φ ≥ 0 | correct (cosmetic: [C] (F3) has log², absorbed in X^δ) |
| Remarks (1)–(3) | 91–95 | — | correct |
| Thm 1: 𝒬′ = ½Σ_c Q_c with Q_c = [R3]'s Q_{c,1} | 106–108 | supp W ⊂ [1/(4ℓ²), 4] = [1/16, 4] (ℓ := 2, [R3] l.79) | correct |
| Thm 1 (i): Lemma 6 (Poisson in T) at w = 1 | 110–112 | [R3] Lemma 6 | correct: c√M ≤ 4N, 256X^{1−3δ} ≤ X; Lemma 6 re-derived (main term, local sums, error) |
| Thm 1 (ii): M = n₂n₂′ odd, cube-free | 113 | [R3] Lemma 7 | correct |
| Thm 1 (iii): 2^ate ≤ E₀ = 23N/(c√X), uniform in M | 114 | S(u) = 0 for u ≥ √8 | correct |
| Thm 1 (iv): E₀² ≤ L ⟺ N ≤ 4cX/529 | 115–116 | [HB95] Thm 1 | correct |
| Thm 1: y = 0 part, secondary part, main part | 117–120 | [R3] Lemma 5(a), §3.6, §3.7, Lemma 9, Jutila | correct; no step uses w > 1 or coprimality to a conductor; (Z) kills the pole also under the restrictions (m, κ) = 1 and the twists ψ(·/q), and for ∂_{t′}b (W^♯ satisfies (Z)) |
| Thm 1 conclusion, X^{3/2+14δ}N^{1−2σ} | 120–122 | Σd(n)⁴\|b\|² ≪ X^δN^{1−2σ} ([R3] l.149; uses N ≤ X^{1−4δ}) | correct; only for the detector sequences b_{j,t′}, ∂_{t′}b_{j,t′} (l.100); false for arbitrary b (numerics below) |
| Thm 2, class II | 132–134 | [R3] Lemma 4(II), Lemma 3(c), [D2] (2.1) | correct: exponent 2 − σ + 6δ; (2.1) is stated on Re s = ½, used at ½+δ (inherited from [R3], same proof, routine) |
| Thm 2, class I | 134–136 | Thm 1, Sobolev, Φ̂_β decay | correct: \|t′\| ≤ X^δ, 2Z ≤ X^{1−4δ}, X^{3/2}V^{1−2σ} = 2^{2σ−1}X^{2−σ} |
| Thm 2 remark: min{2−σ, 10(1−σ)/(3−2σ)}, crossing (√41−3)/4 | 138–141 | [R] Prop 1 | correct |
| Thm 3, step 1′ (δ = (A−1)/8000, κ = δ/4, c = δ/2) | 157–158 | [R3] step 1 constraints | correct |
| Thm 3, step 2′ (Markov, φ ≤ 1 − 8ε′) | 159–163 | Thm 2 (T_X + 3 ≤ X^δ/4), [R] Prop 1 | correct: X^{1−6ε′+62δ} = X^{1−5.38ε′}; F₂ decreasing, F₂(9/10) = 5/6 |
| Thm 3, step 3 (β < σ₀) | 164 | [C] Thm 1(b) | correct: 4Aδ + κ + ε₁ ≤ 5.06δ ≤ 6δ |
| Thm 3, step 4 and [R] §2.5 margins | 165–168 | [R] §2.2–§2.5 | correct: ε′ = 100δ > δ/2, x ≥ \|D\|^{1+2δ}, c ≤ min{ε′, 0.505δ} |
| "Where A > 1 is used" | 170–171 | — | correct |
| §6 reductions (Lemmas P⁺, P^∠; Thms 1⁺, 1^∠, 2⁺, 2^∠, 3⁺, 3^∠) | 175–223 | [R4] | reductions checked (γ/e keeps γ ≫ 0, u_γ, arg; Poisson conditions); [R4] itself not re-refereed |

## Findings (none breaks a proof)

1. **Scope of Theorem 1.** It holds only for the detector sequences, as the NOTE states (l.100) and the discovery file says. For arbitrary b the family constant is ≍ X·N: the trivial character gives about N/√X times the family size, genus characters about XN/q (b = (·/q) on the D with q* | D). Measured below.
2. **Cosmetic.** l.61 lacks |Δ| > 4; Corollary P cites h ≤ CX^{1/2}log X ([C] (F3) has log²); [D2] (2.1) is stated on Re s = ½ but used at ½+δ (inherited from [R3]); [R3] uses F_{κ,q,s} with q = T₁ in (S) and q = T₁T₅ in the triple; §6 reuses ℓ (pole killer ℓ = 2 vs sector exponent ℓ). All absorbed in X^{O(δ)} or notational.
3. **Effective but astronomical X₀.** Lemma 6 needs X^{3δ} ≥ 256, so X₀ ≥ 2^{64000/(3(A−1))} (A = 1.1: 2^{213333}). The theorems are asymptotic; no numerical least-prime bound follows.
4. **Side observation (not claimed in NOTE).** With all genus characters removed, the family constant stays at the family size for arbitrary b up to N = X/2 (λ₁ = 0.88–0.99 × the largest diagonal entry). This suggests a general-coefficient non-genus family sieve. Lead only.

## Numerical falsification tests (all passed)

Own code ([checks/fam.c](checks/fam.c)): fundamental D < 0, X/2 < |D| ≤ X, classes = reduced forms, v_C(n) = r_F(n)/2, b on admissible n (prime factors ≥ 19), Parseval Σ_χ|Σbλ_χ|² = hΣ_C|S_C|², genus characters from assigned characters. Checks: Σ_D h(D) = 5444 at X = 2000 matches PARI qfbclassno; Σ_C v_C(n) = Π_{p|n}(1+(D/p)) on 2000 random (D, n) per run, 0 failures; genus orthogonality on every D. Max memory 259 MB; all runs single-threaded under alarms.

- **Lemma P** ([checks/lemmaP_mine.py](checks/lemmaP_mine.py)): direct (T, y) counts against reduced primitive forms of every Δ/e². 231 orders Δ ∈ [−2600, −7], all non-maximal, all with a nonzero e > 1 block, 14 with Δ_e ∈ {−3, −4}. Max relative deviation 2.7·10⁻¹⁵; no negative block.
- **Theorem 1 on its own sequences** (X = 2000, 4000, …, 64000; N = 2^j from 32 up to ≈ X; σ ∈ {0.51, 0.55, 5/6}; V ∈ {2, √X/2}; t′ ∈ [0, 40], including 14.1347 and 21.022; b and ∂_{t′}b). Normalised by the family size, R1 := LHS/(Σ_D h(D)·Σ|b|²):

| Sequence | R1 (all X, N) | Slope in N on [√X, X^0.95], X = 2000 → 32000 | Slope in X (2000–32000) at N/X = 1/64, 1/16, 1/4, 1 |
|---|---|---|---|
| pole-killed, σ = 5/6, V = √X/2 | 0.81–1.91 (> 1.1 only at N ≲ 2√X) | −0.21 → −0.09 | −0.23, +0.04, −0.02, +0.05 |
| pole-killed, σ = 0.55 | 0.78–1.76 | −0.14 → −0.05 | −0.24, +0.05, −0.01, +0.05 |
| pole-killed, σ = 0.55, V = 2 | 0.76–1.61 | −0.11 → −0.01 | −0.19, 0.00, 0.00, +0.01 |
| pole-killed, t′ = 5 | 1.00–1.06 | −0.01 → 0.00 | 0.00, −0.02, 0.00, 0.00 |
| pole-killed, ∂_{t′}b | 0.79–1.78 | −0.15 → −0.06 | −0.24, +0.05, −0.01, +0.05 |
| pole present (no ℓ-term), σ = 5/6 | 1.35–53.8 | +0.67 → +0.72 | +0.16, +0.08, +0.28, +0.35 |
| pole present, σ = 0.9, V = 2 | 1.31–88.2 | +0.68 → +0.83 | +0.25, +0.42, +0.49, +0.51 |

  At X = 64000 (N = 512, 2048, 8192, 32768) every pole-killed sequence gives R1 = 0.95–1.10, while the pole sequence gives 2.05, 2.99, 11.6, 45.7. t′ scan at X = 16000: max R1 = 1.045. The pole-killed sums grow with exponent 0 in X and ≤ 0 in N, as Theorem 1 says (LHS ≍ X^{3/2}‖b‖²).
- **The dual form of Corollary P** (`fam.c` mode 2; all orders with |Δ| ∈ [X/4, 2X], my own smooth Φ, Lemma P blocks Σ_e (w_e/2)Σ_C|S_C^{(Δ/e²)}(b_e)|²; X = 2000–32000, N up to ≈ X). 𝒬′(b)/(X·Σ|b|²) for the pole-killed sequences (σ = 5/6, 0.55, ∂_{t′}b) is 0.45–0.63 for every N ≥ 2.1√X (0.67–1.21 at N = 1.4√X and 2.0√X), flat in N and X. So the proof's intermediate bound Σ_c|Q_c| ≪ X^{1+ε}‖b‖² is already visible at X ≤ 32000, far below the range X^{3δ} ≥ 256 where its Poisson step applies. This reproduces RCApps5's famq.c (≈ 0.5) with different Φ and support (prime factors ≥ 19 here, prime to 3 there). With the pole present the ratio grows to 41.6 (X = 32000, N = 32768).
- **Pole handling.** With (Z) the trivial character carries 4–15 % of the sum for N ≥ 4√X (18–56 % only for N ≲ 2√X, where a_V removes most of b; still bounded). Without (Z) it carries up to 99.4 %, and the sum grows like N/√X relative to the family size: the term √X|Σb|² of Remark (3).
- **Large-sieve constants** ([checks/eigs.py](checks/eigs.py)). λ₁ is the top eigenvalue of Σ_DΣ_χλ_χλ_χ^* on admissible n ∈ [N, 2N), for X = 2000–16000 and N = X/32 … X/2, in units of Σ_D h(D):

| X | N | all χ | χ ≠ 1 | non-genus χ | largest diagonal |
|---|---|---|---|---|---|
| 2000 | 1024 | 10.62 | 3.30 | 1.38 | 1.49 |
| 4000 | 2048 | 15.05 | 4.50 | 1.37 | 1.56 |
| 8000 | 4096 | 21.14 | 5.92 | 1.32 | 1.37 |
| 16000 | 8192 | 29.65 | 8.11 | 1.31 | 1.33 |

  Measured exponents:
  - all χ: 0.65 → 0.85 in N, 0.35–0.49 in X at fixed N/X = 1/16 … 1/2. So λ₁ ≍ X·N, a factor N/√X above C(X^{3/2} + N)X^ε (the pole).
  - χ ≠ 1: 0.37 → 0.63 in N, 0.19–0.43 in X (genus characters).
  - non-genus: 0.14 → 0.10 in N, tracking the largest diagonal entry (a divisor effect); −0.03 to +0.04 in X.

  Projecting out ψ(n)(n/q) for q ≤ N/40 and ψ mod 8 (the directions of [R3]'s main term 𝔐) brings the full constant down to 1.05–2.04. This matches the proof: an arbitrary b can be large only along these twists, and the detector sequences are small along all of them (Lemma 9).

## Literature check

- `lst.py 'large sieve|class group|binary quadratic form|quadratic forms|least prime|ring class|ideal class|zero density|density hypothesis|Linnik'` returned 24 arXiv title hits (math.NT/CA/CV/SP/PR/FA/CO, 2026-10 and 2026-11) and, on Zenodo, '"large sieve" "class group"' (404 loose matches; newest 15 read by title), '"zero density" "class group"' (415), '"class group characters"' (0) and '"least prime" AND "quadratic form"' (0). None is a large sieve or zero-density result for class group characters on average over D. Abstracts read: 2610.20554 (a large sieve for sums of two squares, Ostmann's problem) and Zenodo 24093928 (split primes and the exponent of class groups, conditional on unrefereed half-planes); RCAudit had already read 2610.04045 (de Faveri), 2610.15989 and 2610.32470.
- Eight arXiv searches (large sieve × class group characters × average over D; zero density × class group L-functions; least prime × ideal class × almost all D; Duke–Kowalski families; Pascadi–Thorner; Ditchen; Heath-Brown–Pierce; density hypothesis × discriminants). Hits:
  - Pascadi–Thorner 2508.14888: large sieve for arbitrary GL_n families, independent of Ramanujan; constant of type N + Q²|S| per the [R3] NOVELTY note, here ≥ X^{7/2}.
  - Ditchen 1312.1502: Bombieri–Vinogradov over D; least prime |D|^{20/3+ε} outside a density-zero set.
  - Heath-Brown–Pierce 1409.3177: the pair-count device, used for ℓ-torsion averages.
  - Duke–Kowalski (Invent. Math. 139, 2000; not on arXiv): automorphic families with conductor losses.
  
  None states a family-size sieve for class group characters or Σ_D N_D(σ) ≪ X^{2−σ+ε}.
- Contradiction? Only for arbitrary coefficients: the trivial and genus characters give an elementary counterexample (measured above), and RCApps5 does not claim that case. GRH implies Theorem 2 outright, so no unconditional lower bound can contradict it.
- Verdict: no prior statement and no contradiction found, consistent with RCApps5 NOVELTY2.

## Reviewer

Three `reviewer` runs on the decisive argument, each with the NOTE lines pasted verbatim ([checks/reviewer_brief.md](checks/reviewer_brief.md), [reviewer_brief2.md](checks/reviewer_brief2.md), [reviewer_brief3.md](checks/reviewer_brief3.md)).
- **Run 1** (98d18730; RCApps5 l.50–122, [R3] Lemma 1, Lemma 6, §3.6):
  - Q1 Lemma P: "CORRECT, SIGNIFICANT".
  - Q2 Corollary P: "CORRECT, SIGNIFICANT".
  - Q3 w = 1 reduction: "CORRECT, ROUTINE … No conductor-coprimality or w > 1 obstruction appears".
  - Q4 "GAP, SIGNIFICANT in the supplied verification": pole killing by (Z) "cannot be checked without the detector definitions, (Z), and Lemma 9".
  - Side notes: l.61 needs |Δ| > 4; the norm comparison Σd⁴|b|² ≪ X^δ‖b‖² was not in the excerpt. "No counterexample to Lemma P or the w = 1 reduction found."
- **Run 2** (277484c4; adds [R3] l.121–150 and l.276–305):
  - G1, pole killing under the restrictions and for ∂_{t′}b: "CLOSED (SIGNIFICANT)".
  - G3, |Δ| > 4: "CLOSED (ROUTINE)".
  - G2, the norm comparison: "OPEN (ROUTINE)", because "absorption into X^δ needs a range such as N ≤ X^C". That range is a hypothesis of Theorem 1 (l.99, "V/4 ≤ N ≤ X^{1−4δ}").
- **Run 3** (79479199; adds the statements [R3] l.138–149 and RCApps5 l.99–102):
  - "the N-range fixes absorption", leaving one point: the support [N/(4ℓ²), 4N] = [N/16, 4N] "without a bound on ℓ". [R3] l.79 fixes it: "Fix the prime ℓ := 2."
  - The run states: "With that established, the displayed argument closes G2."
- Net: every point raised is closed by a quoted NOTE line; ROUTINE. No GAP or FLAW is reported, and none was found by my own derivation (DERIVATION D1–D5).

RH relevance: N7, N11. Theorem 2 is a zero-density estimate on average over a GRH family (an N7 tool on the N11 family), and this audit referees the step that takes the least-prime line from 1.1307 to every A > 1. Records only; no RH node moves. Not re-refereed here: [R] Proposition 1 (the Humphries–Thorner detector, used for β > 9/10) and [R4] (D > 0, sectors).
