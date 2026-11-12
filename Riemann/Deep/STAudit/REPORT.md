AUDIT PASSED: the Sato–Tate-in-every-ideal-class chain (STLemP Theorem 1′, Lemma L, Lemma W, Theorem MS, Theorem 2₁, Proposition F₂, thresholds 1, (23 − 3√5)/8, 3; ChebST's 11/3, (M+3)/3 and minorants; RSFam's sign machinery) is correct as stated (one routine omission in Theorem MS's S₂ paragraph, repaired with no change), and nothing breaks at small sizes.

# STAudit: referee of STLemP, ChebST and RSFam as used by [discoveries/sato-tate-classes.md](../../../discoveries/sato-tate-classes.md) (2026-11-12)

Brief: `.scratch/orca/staudit-brief.md`. Derivations: [DERIVATION.md](DERIVATION.md) (D1–D8). Scripts, outputs, reviewer briefs: [checks/](checks/). Line numbers: [S] = [STLemP NOTE](../STLemP/NOTE.md), [CS] = [ChebST NOTE](../ChebST/NOTE.md), [RS] = [RSFam NOTE](../RSFam/NOTE.md), [R3] = [RCApps3 NOTE](../RCApps3/NOTE.md). Refereed inputs taken as given: [RC5Audit REPORT](../RC5Audit/REPORT.md), [REPORT2](../RC5Audit/REPORT2.md) (Lemma P, Corollary P, [R5] Theorem 1 at w = 1, Lemma P⁺, Theorem 1⁺, [R4]).

## Verdict table

| Step | Lines | Inputs | Verdict |
|---|---|---|---|
| Class 𝔅, Theorem 1′ (family sieve for sequences with controlled quadratic twists) | S 56–91 | [R3] §3.4–§3.8 at w = 1 (refereed), [R3] Lemma 9 → (B2) | correct: b enters [R3] only through its support, Σd⁴\|b\|² (from (B1)) and §3.7's B_{κ,q}(1), B′_{κ′,q}(1); (B2) bounds both (characters real, so b̄ is covered); E^θ by dyadic partial summation (D1) |
| Lemma L (ring-class lift) | S 95–110 | Lemma P (refereed), h(p₀²D) = h(D)(p₀ − (D/p₀)) | correct (D2) |
| Corollary L (beyond the conductor) | S 112–118 | Lemma L, Theorem 1′ | correct: X′ ∈ [N^{1+5δ}, 4N^{1+5δ}], N ≤ X′^{1−4δ} |
| Lemma L⁺ (D > 0) | S 120–124 | [R4] Lemma 2⁺, [R5] Lemma P⁺, Theorem 1⁺ (refereed) | correct as a reduction |
| Lemma H (quadratic-twist mean values) | S 134–151 | IK Thm 5.3, Heath-Brown 1995 | correct, standard: θ = 0 (f), θ = ½ (Sym²f) |
| Lemma W(i) (twisted detectors ∈ 𝔅) | S 153–177 | L(s, Sym^m f ⊗ ψ̃) entire for m ≤ 2, every ψ̃ (f non-CM); Lemma H | correct (D3) |
| Lemma W(ii) (AFE pieces ∈ 𝔅), incl. the 2026-10-22 N = 1 repair | S 157–158, 178–184 | — | correct |
| Theorem MS: Σ_{D,χ}\|L(½+δ+iv, Sym^m f ⊗ θ_χ)\|² ≪ X^{3(m+1)/4+60δ} | S 188–220 | Lemma L, Theorem 1′, Lemma W(ii), AFE | correct; finding 1 (S₂ outer weight s^{−½+δ} omitted; routine; bound unchanged) (D4) |
| Detector, class II jointly over (D, χ), class I | S 224–236 | [CS] Lemmas F1, F2 | correct |
| Theorem 2₁: Σ_{D,χ}N₁(σ) ≪ X^{2−σ+40δ} | S 238–246 | Theorem MS, Lemma W(i), Theorem 1′ (p₀ = 1) | correct (D5) |
| Proposition F₂: X^{15(5−4σ)/(16(2−σ))}, ¾ ≤ σ ≤ 13/14 | S 252–268 | Theorem MS (m = 2), Corollary L per block | correct; exact (`thr_check.py`) and optimiser (`thr_opt2.py`: E* = F₂) |
| Proposition T¹ (one-sided, every A > 1) | S 276–294 | Theorem 2₁ (σ ≤ 9/10), [RS] Prop F (σ > 9/10), [RS] §5 | correct: margins re-derived (D6) |
| Proposition T² (signs, every A > (23 − 3√5)/8) | S 296–306 | Prop F₂, [RS] Prop F (σ > 13/14) | correct: max f = f((1+3√5)/44) = (23 − 3√5)/8 exactly; σ₀ = 0.75448 ≥ ¾ |
| m = 3: free Y in [CS] Prop F′, Proposition T₃ (A > 3) | S 312–325 | [CS] Prop F′, F″; PT Cor 1.2 | correct, conditional on the PT preprint (as stated) |
| Theorem D′, Corollary D′1 | S 327–341 | [RS] Theorem D, minorants x(x+2)/8, (α+x)/(2+α) | correct |
| Sato–Tate table, discrepancy 10/(3A−3), A(J) ≤ max{11/3, 25.2/L + 1.7} | S 343–348; CS §1, §4, §5 | [CS] Theorems 1, 2, minorants (1.1) | correct; the A(J) bound holds with room (finding 3) |
| [CS] Prop F′/F″, (LS1)–(LS3), routes H/H′ (11/3, (M+3)/3) | CS 165–200, 252–351 | HT Thm 4.1, PT Thm 1.1 and Cor 1.2 (TeX) | correct: per-D ≤ X^{3/2}, PT on [X^{3/2}, X^{11/4}], HT above; 11/3 is the method's limit (S 370–373) |
| [RS] Prop F (route F, X^{(m+10)e(σ)}) and §5 bookkeeping (E₁, E₂, (a)–(f)) | RS 229–240, 316–375 | HT bound_int + cor:MVT, BTZ, (J5) (TeX) | correct as used (only near σ = 1 in T¹, T²) |
| Imported statements | — | PT 2508.14888v2, HT 2103.05634v2, BTZ 1804.06402v3, Thorner 2404.06482v2 (J5), read in TeX | as cited (D7); PT is a preprint |
| D > 0, narrow classes | S 350–356 | [R4], [R5] §6 (refereed), [RS] §6, [CS] §6 | reductions correct; not re-derived in full |
| Limits, GRH paragraph | S 358–381 | — | correct |

## Findings (none changes a statement)

1. **Theorem MS, S₂ (S 217–220), routine.** "S₂ is treated identically" omits that the outer Cauchy–Schwarz weight becomes d_{2m+2}(s)s^{−½+δ}. With 𝔖₂ ≤ (4P)^δ𝔖 and N ≤ 4P/s, the S₂ family sum is ≤ C(4P)^{3δ}(log X)^{C}X′^{3/2+21δ}, and (8P₀)^{3δ}\|ε_m\|² ≤ C (δ ≤ (m+1)/24): the S₂ bound equals the S₁ bound. The NOTE's own count ("N^{4δ}") reaches the same X^{(2m+2)δ}. Reviewers: B "GAP — SIGNIFICANT … the claimed exponent survives", C (on the repair) "GAP — ROUTINE … MS's exponent is unchanged".
2. **Discovery page, cosmetic.** "With published inputs only, 11/3 becomes 5": without PT, 𝒜(3) = 10/3, 𝒜(4) = 11/3 still (route F), 𝒜(M) = (M+7)/3 for 5 ≤ M ≤ 8, 5 for 8 ≤ M ≤ 12, (M+3)/3 above ([CS] l.37). "No bound with any power of |D| was known" means unconditionally (under GRH the method itself gives every A > 1 for every D, S 375–378).
3. **A(J) bound is loose, not wrong.** For 4812 intervals (`minor.py`) the least n with μ_ST(J) > λ_n(α) + λ_n(β) has nL ≤ 10.0 (≈ 3π, at an end of [0, π]), against the 37.7 used; M*(J) ≲ 20/L + 2.
4. **Uniformity.** Every constant depends only on (f, m, δ, A, J); D enters only through 𝔇(X). X₀ is astronomical (Lemma 6 needs X′^{3δ} ≥ 256, δ ≤ (A−1)/8000 or smaller): no numerical bound follows (as for RCApps5).

## Numerical tests (E = 11a1, f = f_E; every process < 290 MB under `memcap.sh`, single-threaded)

**(a) L-values** ([checks/msq.gp](checks/msq.gp), [msq2.gp](checks/msq2.gp); PARI 2.19 lfun from the Euler product of E/K ⊗ χ, conductor 121D², root number −(D/11); lfuncheckfeq ≤ −51 bits on every L-function; trivial χ equals L(E,s)L(E^D,s); 11 ∤ D; 170 discriminants, 1259 L-functions, |D| ≤ 600, s = ½ + δ + it analytic). S := Σ_DΣ_χ\|L\|² per dyadic window, divided by the family size Σh:

| window | #D | Σh | (δ,t) = (0,0) | (0,1.5) | (.05,0) | (.05,1.5) | (.15,0) | (.15,1.5) |
|---|---|---|---|---|---|---|---|---|
| (37.5, 75] | 11 | 35 | 11.85 | 4.83 | 10.03 | 3.54 | 7.36 | 2.32 |
| (75, 150] | 20 | 93 | 7.70 | 12.36 | 6.21 | 8.34 | 4.44 | 4.61 |
| (150, 300] | 44 | 300 | 8.21 | 12.61 | 6.28 | 7.99 | 4.36 | 4.19 |
| (300, 600] | 84 | 821 | 12.79 | 14.59 | 9.17 | 8.76 | 5.73 | 4.39 |

- Measured exponents of S in X over the four windows: 1.58, 2.02, 1.50, 1.92, 1.42, 1.80 (Σh: 1.53); over the last two windows 2.09, 1.66, 2.00, 1.58, 1.85, 1.52 (Σh: 1.45); t = 4 (|D| ≤ 400 only, windows (100, 200] → (200, 400]): 1.63–1.85.
- S/Σh fits a + b log X better than cX^{1/2} in all six columns (RSS 1.7–19 against 4.9–52). The X^{1/2} growth of the per-D bound (×1.41 per window) is absent at t = 1.5 (12.36 → 12.61 → 14.59; 4.61 → 4.19 → 4.39).
- Limits: at |D| ≤ 600 the AFE length ≈ 2.5|D|(\|t\|+3)²/9 still exceeds the window's family size, and half the D have root number −1 (L(½) = 0 for all χ), so the central columns fluctuate. Consistent with X^{3/2+o(1)}; no contradiction. Sym² L-values (degree 6, conductor 121²|D|³) were not cheap; Sym² is tested in (b).

**(b) The mechanism of Theorem MS: family sums of AFE pieces** ([checks/fam2.c](checks/fam2.c), [fam3.c](checks/fam3.c); own code, all fundamental D < 0 in (X/2, X], classes = reduced forms, r_C(n) = R_Q(n)/2 from lattice points; Σh checked against PARI: 5444 at X = 2000). b(n) = coef(n)n^{−σ}φ(n/N), n squarefree with primes ≥ 19, φ a bump on [1, 2]; R1 := F/G with F = Σ_D hΣ_C\|S_C(b)\|² and G its pure diagonal Σ_D hΣ_n\|b(n)\|²r_D(n):

| X | Σh | N/X | λ_E n^{−½} | λ_{Sym²E} n^{−½} | a_Vλ_E n^{−0.55} | a_Vλ_{Sym²} n^{−0.8} | pole (coef 1) |
|---|---|---|---|---|---|---|---|
| 4000 | 15217 | 1/16 … 2 | 1.010–1.055 | 1.044–1.117 | | | 2.24 → 42.7 |
| 8000 | 43175 | 1/16 … 2 | 1.014–1.038 | 1.007–1.097 | 1.018–1.043 | 1.016–1.138 | 2.82 → 59.9 |
| 16000 | 122446 | 1/16 … 2 | 1.003–1.022 | 1.032–1.076 | | | 3.60 → 84.4 |
| 32000 | 346139 | 1/16 … 2 | 1.003–1.018 | 1.025–1.044 | 1.014–1.022 | 1.029–1.066 | 4.76 → 119.4 |

- Twisted sequences: the family sum equals its diagonal within 0–14 % up to N = 2X (twice the conductor), exponent 0.00 in X and in N; F itself grows like Σh ∝ X^{1.50}. Per D, max F_D/G_D ≤ 5.7 (pole sequence: up to 690). This is Lindelöf on average for the double family at the level of the AFE pieces, which is what Theorem MS proves.
- Beyond the conductor ([checks/fam3_beyond.out](checks/fam3_beyond.out)): N/X = 4, 8, 16 at X = 2000, 4000 and N/X = 16, 32 at X = 1000 (N = 32000 ≈ X^{3/2}, the Sym² AFE length): twisted R1 = 1.00–1.20, pole R1 = 61–340. The proved Sym² bound X^{9/4} (via the lift) is not sharp here; the data sit at the family size, as STLemP §8 also found.
- Pole sequence: R1 grows like N (slopes 0.86–1.0) and like X^{0.49} at fixed N/X (21.5, 30.2, 42.4, 59.8 at N = X): the √X\|Σb\|² term that class 𝔅 excludes. The test detects the failure mode it is meant to.

**(c) Signs in every class** ([checks/sign.c](checks/sign.c); 𝔭 above p located by qfbprimeform + reduction, 𝔭̄ in the inverse class; a_p from PARI, p < 2·10⁶). For each D, p^± := max over classes of the least p ∤ 11D with a prime above p in the class and ±a_E(p) > 0:

| \|D\| | #D | max log p⁺/log\|D\| | max log p⁻/log\|D\| | a_p > −√p/2 | p^± > 2\|D\|^{2.0365} | min 2\|D\|^{2.0365}/max(p⁺, p⁻) | a_p > −√p/2 needs p > 2\|D\| |
|---|---|---|---|---|---|---|---|
| 3, 4 (by hand) | 2 | p⁺ = 13, 5 | p⁻ = 7, 17 | p = 13, 5 | 0 | 1.44 (D = −3), 1.98 (D = −4) | 1 (D = −3) |
| 5–200 | 60 | 1.961 (−8) | 1.632 (−23) | 1.611 (−7) | 0 | 2.34 (−8) | 12 |
| 201–1000 | 243 | 1.393 | 1.369 | 1.364 | 0 | 97 (−335) | 14 |
| 1001–5000 | 1219 | 1.293 | 1.320 | 1.182 | 0 | 297 (−1079) | 30 |
| 5001–20000 | 4555 | 1.218 | 1.201 | 1.158 | 0 | 2848 (−5879) | 19 |
| 20001–60000 | 12159 | 1.129 | 1.132 | 1.105 | 0 | 17430 (−20135) | 3 |
| 60001–150000 | 27360 | 1.106 | 1.123 | 1.104 | 0 | 57250 (−76079) | 3 |

- The sign claim holds for every one of the 45,598 discriminants with |D| ≤ 150000; tightest D = −3 (needs p = 13, bound 18.7), then D = −4 and D = −8. The worst exponent falls from 1.96 to 1.11–1.12, toward the GRH value 1.
- One-sided (α = ½) at A = 1 exactly: a few D per range need p slightly above 2|D| (3 of 27,360 at the top); allowed (A > 1, exceptional set, astronomical X₀), and the worst exponent 1.10 is falling.

## Literature check

- `lst.py 'Sato|ideal class|class group|binary quadratic|Rankin|least prime' --zenodo '"Sato-Tate" "class group"'`: 24 arXiv title hits (2026-10, 2026-11), Zenodo 345 loose matches (newest 15 by title). Abstracts read: 2610.36464 (Amri) and 2610.42022 (El Kaouni), effective Chebotarev–Sato–Tate, both conditional on GRH-type/motivic conjectures; 2610.32470 (Shintani lifts), 2610.18491 (quasi-RH, small generators). Zenodo 24093928 was already read by RC5Audit. None concerns Sato–Tate conditions on primes in ideal classes of varying fields, or family mean values over D and χ.
- 12 arXiv searches (approach "arxiv"; 2 refused as irrelevant) and 2 GXL Evidence API claims (the double-family mean square; signs in every class): no supporting or refuting source. GXL's two "refutes" are mis-matches (0907.4747, the quadratic-twist family L(½, f⊗χ_{8d}); math/0510197, a per-field count with an error term).
- Related, not overlapping: Michel–Venkatesh (per field, first moment, non-vanishing); Nordentoft 2105.09130 (wide moments per field via Waldspurger); Blomer–Brumley 2009.07093 (GRH, per field); Jana 2012.07817 (GL(n)×GL(n), different family); Munsch–Toma 2503.05269 (quadratic character sums); Ditchen 1312.1502 (primes in forms on average). Duke–Friedlander–Iwaniec (class group L-functions, per field) was not found on arXiv.
- Verdict: nothing states Σ_{|D|~X}Σ_χ\|L(½+δ+it, f⊗θ_χ)\|² ≪ X^{3/2+O(δ)}, the lift, or the thresholds; nothing contradicts them (they are GRH-consistent, and no unconditional lower bound beyond the family size is known). Consistent with STLemP NOVELTY.md.

## Reviewer

Three `reviewer` runs, briefs with the NOTE and [R3] lines pasted verbatim ([checks/reviewer_A.md](checks/reviewer_A.md), [reviewer_B.md](checks/reviewer_B.md), verdicts in [reviewer_verdicts.md](checks/reviewer_verdicts.md)):
- A (Theorem 1′ against [R3] §3.5–§3.7): "VERDICT: CORRECT, ROUTINE … (B2) is a hypothesis to be verified for each sequence: this is a routine conditional extension, not a new cancellation theorem".
- B (Theorem MS given Theorem 1′ and Lemma L): "GAP — SIGNIFICANT" for the S₂ outer weight, with "Thus the claimed exponent survives; the theorem is not refuted … The significant Lindelöf on average conclusion follows after the explicit S₂ repair".
- C (the repair of finding 1): "GAP — ROUTINE … the repair works with mechanical clarifications; MS's exponent is unchanged" (its clarifications, N ≤ 4P/s and δ ≤ (m+1)/24, are in D4).
- Net: one routine bookkeeping gap, closed; no FLAW.

## Not re-derived here

[RS] (R1)–(R4), Lemmas A, B (read, standard), [CS] Lemmas M1–M5 (classical Chebyshev–Markov–Stieltjes; their consequence checked numerically), the Rankin–Selberg local theory (J1), (J2), the D > 0 Theorem MS beyond its reduction, and the sector variant (not claimed by the page).

RH relevance: N11. The chain is a GRH-family zero-density result on average (density hypothesis on average for f ⊗ θ_χ) with least-prime and Sato–Tate consequences; it now stands refereed down to RC5Audit's inputs and the published HT/BTZ/Thorner theorems (plus the PT preprint for M ≥ 3). Records only; no RH node moves.
