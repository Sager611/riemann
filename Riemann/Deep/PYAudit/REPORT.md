AUDIT PASSED: every lemma our 54.6 t^{27/164} certificate takes from Patel–Yang 2302.13444v1 (1.2, 1.4 = Yang 2301.03165 Lem 2.3–2.5, 2.1–2.3, 3.1 = Hiary Lem 2.1/Gabcke, and the Lemma 3.4 display l. 826–831) is true as used. Each was re-derived line by line, and every decimal constant is rounded the safe way. About 6.2·10³ new numeric cases, plus 4.5·10⁴ checks of the r-sum inequality, found 0 violations. New findings: 5 proof gaps (G1–G5), each closed without changing a constant, and 7 misprints. One dependence: 54.6 needs the recursion reading of Yang's proof (the other reading of the printed δ₃ gives 55.12). The only error that changes a number is still our Lemma 3.5 μ₂ slip.

# PYAudit: referee report on the lemmas of Patel–Yang 2302.13444 (node N9, load-bearing input of the 54.6 record)

Sources: 2302.13444v1 e-print `main.tex` (line numbers "l."); Yang 2301.03165 v1 (10 Jan 2023, the version PY cite; "Y1 l.") and v2 (16 Jan 2024 = J. Math. Anal. Appl. 534 (2024) 128124; "Y2 l."); Hiary 1507.01261v4 `explicitCorput.tex` ("H l."); Gabcke via RevAudit (`Riemann/Deep/RevAudit/REPORT.md` §3); Dhiman–Kadiri–Quesada-Herrera 2609.00537 ("DKQ l."). Scripts and logs are in this directory. R438 had already falsification-tested Lemma 1.4 (k = 3, 4, 5, log and polynomial phases), Lemma 2.3 (log blocks) and Lemma 2.2; those tests are not repeated.

## 1. Statements, hypotheses, use in `SubWeyl/certify.py`

**Reading of S_f.** Yang defines S_f(a, N) := |Σ_{a<n≤a+N} e(f(n))| (Y1 and Y2 l. 101–102). PY's own definition (l. 78) drops |·|, but Lemmas 1.2 and 1.4 are Yang's, word for word. So S_f is the modulus. Lemma 1.2 also holds with Re S_{g_r}, because its proof is Cauchy–Schwarz over N + q − 1 shifted windows.

| PY lemma (l.) | Statement (S = modulus) | Hypotheses | Use in certify.py |
|---|---|---|---|
| 1.2 `weyl_differencing` (100–108) = Y1 Lem 2.3 (317–327) | S_f(a,N)² ≤ (N−1+q)(N/q + (2/q)Σ_{r=1}^{q−1}(1−r/q)S_{g_r}(a,N−r)), g_r(x) = f(x+r) − f(x) | a, N integers; any integer q > 0 | display l. 825–831 → `interval_bound` 'aba' (F1·F2); `tail_bound` C1…E7b |
| 1.3 `poisson_summation_formula` (110–125) = Patel 2021 Thm 2.31 | Poisson formula with 40/√π·λ₂^{−1/2} + A₁log(β−α+4) + A₂(b−a)λ₂^{1/5}λ₃^{1/5} + A₃ | — | **not used** (PY replace it by Lemma 2.3, l. 179) |
| 1.4 `kth_deriv_test` (129–152) = Y1 Lem 2.4 (330–401) + 2.5 (403–536) | S_f(a,N) ≤ A_k h^{2/K}Nλ_k^{1/(2K−2)} + B_kN^{1−2/K}λ_k^{−1/(2K−2)}, K = 2^{k−1}, every η > 0 | a, N integers, N > 0; f ∈ C^k, f^{(k)} monotone; λ_k ≤ \|f^{(k)}\| ≤ hλ_k on (a, a+N], h > 1 | `A3B3`, `AkBk`: k = 5 → c2, c3, c7a, c7b (h = HC5·H⁹·κ); k = 4 → `fth`, tail D3/D4 (h = H⁴) |
| 2.1 `titchmarsh_lem_47` (212–338) | \|Σ_{a<n≤b}e(f(n)) − Σ_{α−½<m<β+½}∫_a^b e(f(x)−mx)dx\| ≤ (3/π)log(β−α+2) + 4 | f′ continuous, decreasing; α = f′(b), β = f′(a) | inside 2.3 (→ E6, the log term) |
| 2.2 `titchmarsh_lem_46` (340–465) | \|∫_a^b e(f) − e(f(c)−⅛)/\|f″(c)\|^{½}\| ≤ (2·3^{2/3}/π^{2/3})h₃^{1/3}λ₂^{−1}λ₃^{1/3} + (1/π)(1/\|f′(a)\| + 1/\|f′(b)\|) | f″ < 0; λ₂ ≤ \|f″\| ≤ h₂λ₂, \|f‴\| ≤ h₃λ₃ on (a,b); f′(c) = 0, c ∈ [a,b] | inside 2.3 (→ c5b, log term) |
| 2.3 `imp_poisson_summation_formula` (467–507) | \|Σe(f(n)) − Σ_{α<ν≤β}e(f(x_ν)−νx_ν−⅛)/\|f″(x_ν)\|^{½}\| ≤ 4.686λ₂^{−½} + (2·3^{2/3}/π^{2/3})h₂h₃^{1/3}(b−a)λ₃^{1/3} + (5/π)log(β−α+2) + 6 | f ∈ C³, f′ decreasing on [a,b]; λ₂ ≤ \|f″\| ≤ h₂λ₂; \|f‴\| ≤ h₃λ₃ | c4 = (8/3)·4.686√π H^{3/2}; c5a, c5b; E6 = 6 + (5/π)log 2 |
| 3.1 `riemann_siegel_formula` (513–521) = H Lem 2.1 (274–310) | \|ζ(½+it)\| ≤ 2\|Σ_{n≤⌊√(t/2π)⌋}n^{−½−it}\| + 1.48t^{−¼} + 0.127t^{−¾} | t ≥ 200 | `rs` in `interval_bound` (L₀ ≥ 75) |
| 3.4 `exponential_sum_lemma` (636–658), proof 660–861 | display l. 826–831: \|S_f(a,b−a)\|² ≤ (h−1+q/a)((h−1)a²/q + (1800/2911)c₃a^{4/5}(qt)^{11/30} + (28800/54481)c₄a^{9/20}(qt)^{61/120} + (8/3)E₄a^{5/2}(tq)^{−½} + (9/14)E₅a^{2/3}(tq)^{1/3} + E₆a) | a < b ≤ ha, 1 < h ≤ 2, integer q ≥ 1 (SubWeyl REPORT A) | Part C per block (E₅(a), repairs 3a, 3b); Part T in stated form |
| 3.2, 3.3, 3.5 (557, 576, 863) | S₁, S₂, S₃ bounds | | Part T (`tail_bound`) only, with the SubWeyl corrections |

## 2. Proof checks

### Yang 2301.03165 (source of 1.2, 1.4)
- **Lem 2.1 (second-derivative test, Y1 l. 130–287): correct.** k ≤ Nhλ₂ + ε₁ − ε₂, y_j − x_j ≤ 2Δ/λ₂, Kuzmin–Landau 2/(πΔ) and boundary cases checked. Y1 misprints: "f > 0" should read f″ > 0; "S_f(a, y₁ − a)" should read y₀. In J(ε), Case 1 counts 2/(πΔ₀) where H has 1/(πΔ₀), and Case 2 uses −1/(πΔ₀) in J(Δ₀); both are conservative. The Case-3 display (Y1 l. 277–281) drops H = −1/(πΔ₀), so its "< 2/√(πλ₂) − ½" is false as printed; with H, J < 1/(πΔ₀) + ½ ≤ 2/(πΔ₀). Y2 (l. 165–303) fixes the f″ typo and the J(ε) bookkeeping; "y₁" is still there (Y2 l. 210). The conclusion holds.
- **Lem 2.3 (Weyl): correct** (Cauchy–Schwarz, N + q − 1 windows).
- **Lem 2.4 (k = 3): correct.** Checked:
  - the trivial cases λ₃ ≥ λ₀ and λ₃ ≤ λ₁ = (η/(κN))³, with 32κ(1+κ) = 3√πη^{3/2} and δ₃ = √(1+κ);
  - the main case: Lem 2.1 on g_r with λ₂ = rλ₃, where g_r″ is monotone because f‴ is;
  - the r-sums 32/(15√π), 1/3 and 32/(3√π), and q = ⌈ηλ₃^{−1/3}⌉;
  - the final factor 1 + κ.
- **Lem 2.5 (k ≥ 4): correct**, with one gap that the trivial bound closes:
  - Edge N ≤ 2336: it needs 2^{25/3}(32/(3√(5π)))² = 2336.258 > 2336 and 64(K−1)⁴ − (2K−1)(2K−3)(4K−3)(4K−5) = 20K² − 40K + 19 > 0. Both hold.
  - Edge λ ≤ λ₀(k): it needs (B₃/δ₃)⁴·9πη/1024 = 1, which holds. Y1 l. 454 prints the reciprocal of B₃/δ₃; Y2 fixes it.
  - Main case: α_j, β_j, D_j ≤ δ_j (from q/N ≤ 2(9πη/1024)^{1/J}N^{2/J−1}), the exponent 1/(4J−2), 2^{7/12}α_j^{½} = 2^{19/12}(J−1)/√((2J−1)(4J−3)) and (2β_j)^{½} = 2^{3/2}(J−1)/√((2J−3)(4J−5)) all check.
  - **G1 (Y1 l. 500, Y2 l. 526):** "λ^{−1/(2J−1)} < q ≤ 2λ^{−1/(2J−1)}" needs λ_{j+1} ≤ 1. When λ_{j+1} > 1 the trivial bound gives the claim, because A_{j+1} ≥ δ_jh^{−1/J} > h^{−1/J}.
- **Reading of the recursion.** PY l. 136–150 (like Y1 l. 403–415) define δ₃ twice: once in A₃, B₃ and once by the δ_k formula at k = 3. Yang's proof settles which is meant (D_j ≤ δ_j, Y1 l. 505, Y2 l. 531):
  - A₃ and B₃ carry the third-derivative δ₃ = √(½ + ½√(1 + (3/8)√πη^{3/2})).
  - The step A₃ → A₄ uses δ₃ = √(1 + 2·2337^{−½}(9πη/1024)^{¼}).
  - Y2 l. 434–441 states this explicitly (J := 2^{j−1}).
  - `certify.py` uses this reading. An independent mpmath implementation written from the proof agrees with it to 3·10⁻³⁸ (`akbk_check.py`).
  - The other reading inflates A₅ by 1–4 % and B₅ by 2–7 %. It raises the certificate to 55.11965 on log t ∈ [648, 650], against 54.59967 (`reading_check.py`), and the tail to 52.86.
  - **So 54.6 needs the proof's reading, which is the correct one.** The other reading is merely more conservative.
- **Y1 → Y2 (reviewer revision):** no constant of §2 changed. The diff is notation, the J index, the B₃ reciprocal, the J(ε) bookkeeping and a new remark. (R438 said "v1 only"; v2 exists.)

### PY §2
- **Lem 2.1: correct.** Checked: Euler–Maclaurin; the sawtooth Fourier series; the digamma bounds; S₁₂ via m₀ and n = ⌊m−β⌋. The total is (10 log 2 + 1 + γ + 9/32)/π + 1 = 3.798 < 4, and the step γ + 2 = 2.5772 ≤ ½ + 3 log 2 = 2.5794 holds.
  - **G2 (l. 219, 335):** the normalisation −½ < α ≤ ½ fails at α = ½. There the m = 0 integral is in the identity but outside (α−½, β+½). Normalising to −½ ≤ α < ½ instead makes every step hold (≤ at α = −½).
  - **G3 (l. 244–254):** the second mean-value theorem is applied to a complex integrator. The bound 2\|g(a)\| + 2\|g(b)\| still holds: split at the sign change of the monotone g, and on each piece \|∫g dφ\| ≤ 2 max\|g\|.
- **Lem 2.2: correct**, with one gap. Checked: Rogers' bound; the Fresnel tail \|∫_X^∞e^{−iu}u^{−½}du\| ≤ X^{−½} (by contour rotation); I₄–I₆; the choice of δ.
  - **G4 (l. 366–370, 411, 455–464):** when c ± δ ∉ [a, b] the proof integrates over [c−δ, c+δ] and uses \|f‴\| ≤ h₃λ₃ and f″ < 0 outside (a, b). Repair: extend f past a and b with f″ held constant (C^{2,1}, Taylor with integral remainder). The extension keeps f″ ≤ −λ₂ and \|f‴\| ≤ h₃λ₃ a.e., so the statement is true. In PY's use, g_r lives on (0, ∞) anyway.
  - Misprints: l. 452 has λ₃^{−1/3} for λ₃^{1/3}. The chosen δ balances the two terms but does not minimise them (valid). l. 186–190 quotes Rogers' bound without "f′ monotone", which is false without it; every use has f′ monotone.
- **Lem 2.3: correct**, with one gap in the proof.
  - Constants: 4.686 = 2 + 2·1.343, where 1.343 ≥ 1.342169 is the quadratic Kershner value; √(2π)·1.342169 = 3.36431 matches PY's recomputed γ₀ = 3.3643 (l. 194 footnote). The T₃ term is (number of ν, which is < β−α ≤ (b−a)h₂λ₂) times the Lemma 2.2 term. Also 5/π = 2/π + 3/π and 6 = 2 + 4.
  - **G5 (l. 490–499):** the S₃ chain "Σ1/(ν−α) < Σ_{1≤n<β−α−½}1/(n−½) < ψ(β−α) − ψ(½)" is false as printed (`s3_examples.py`):
    - The n-range must be n < β−α. At α = 0.4, β = 3.6 the left side is 2.6763 > 2.6667.
    - The ψ bound fails too. At α = 0.495, β = 3.505 the sum is 3.0439 > ψ(3.01) − ψ(½) = 2.8902. At α = ½ − 10⁻⁷, β = 3/2 + 10⁻⁷ the sum is 2.0000 > 1.3863.
    - Correct version: each sum is ≤ ψ(N+½) − ψ(½) with N = ⌈β−α⌉ − 1 (the sharp sup, approached as ν₁−α and β−ν_N → ½⁺). So \|S₃\| ≤ (2/π)(ψ(N+½) − ψ(½)). The supremum of this minus (2/π)log(β−α) − 2 is −0.72676, approached as β−α → 1⁺ (−0.75 as β−α → ∞; `s3_examples.py`).
    - For β−α ≤ 1, S₃ is empty. The printed (2/π)log(β−α) + 2 is then negative when β−α < e^{−π}, so the final step must use log(β−α+2) there, which it can.
    - The lemma's constants stand.
    - Misprint: l. 490 has "a ≤ ν ≤ b" for α ≤ ν ≤ β.
- **Lem 3.1 = Hiary Lem 2.1 (H l. 274–310): the same statement** (Hiary writes n^{−½+it}, the conjugate). Hiary uses Gabcke p. 9 with \|R₀(t)\| < 0.127t^{−¾} for t ≥ 200, and β = cos(π/8), so β(2π)^{¼} = 1.462717 < 1.48. RevAudit §3 checked Gabcke's 2015 edition, Satz 3.2.2(b): t-scale, t ≥ 200, MATCH, and max\|R₀\|t^{¾} = 0.12229 at 856 heights. Not repeated here.

### PY §3 (Lemma 3.4 up to the display; Lemmas 3.2, 3.3, 3.5)
- **Re-derived and confirmed:**
  - x_ν, and the closed form of \|φ_r⁽⁵⁾\| at l. 694–697 (it matches the numerical 5th derivative to 15 digits, `falsify.py phi5`);
  - x ≤ ¼;
  - (8x³+36x²+60x+35)/(x+2)^{7/2} ∈ [6704/2187, 35/(8√2)] and decreasing; x·P′/P ∈ [−0.0094, 0], so φ⁽⁵⁾ is monotone;
  - λ₅ = (53632√2π⁴/729)a⁹/(tr)⁴ and h₅ = (76545√2/107264)h⁹ (exact);
  - (N_bound); c₁, c₂, b₁, b₂; the weight ≤ π^{½}h^{3/2}a^{3/2}(tr)^{−½}; c₃, c₄;
  - λ₂ = tr/(πh³a³), h₂ = h³, λ₃ = 3tr/(πh⁴a⁴), h₃ = h⁴; T₁, and T₃ ≤ (6/π)h³(h−1)(tr)^{⅓}a^{−⅓};
  - log(2+x) ≤ log 2 + 0.893x^{⅓} (max 0.892413 at x = 31.6);
  - the r-sum factors 1800/2911, 28800/54481, 8/3, 9/14, 1, with the inequality Σ_{r<q}(1−r/q)r^s ≤ q^{1+s}/((1+s)(2+s)) proved (concavity for 0 < s ≤ 1, monotonicity for −1 < s ≤ 0) and checked to q = 10⁵;
  - the display at l. 826–831, every exponent of C₁…E₃, and the exponents of Lemma 3.5.
  - 0.4750 ≥ 0.474985 and 0.2531 ≥ 0.253040, so both are rounded the safe way.
- **New misprints (harmless):**
  - l. 725 gives "N = b − a" for the length of the ν-sum, which is β − α. Repairs 3a and 3b of the SubWeyl REPORT supply the "+1" and the ν-range below α.
  - l. 746 and 776 cite Lemma 2.2 where Lemma 2.3 is applied.
  - l. 858 drops the factor E₅, which the statement's E₂ contains.
  - l. 547 "a₀′ ≤ a_K" is reversed (a_K ≤ a₀′). This is harmless, because the S₂ block bounds are maxima over prefixes.
- **Known items (SubWeyl REPORT), each confirmed:**
  - μ₂ (l. 876, 921–923): a factor (2πh)^α is dropped;
  - h₃ has θ₂ for θ₃ (l. 581); h₀ (l. 864); C₀ (l. 560);
  - the E₁ exponent (l. 849) and the q₀ exponent (l. 838);
  - the last-block hypothesis (l. 535 vs 637);
  - A₄(η₂, h₃), which is conservative;
  - μ₄'s subtracted term (l. 891), invalid only for log t ≳ 1.7·10⁴.

## 3. Falsification (new; complements R438)
| Test (`falsify.py` mode / script) | What | Cases | Max LHS/RHS |
|---|---|---|---|
| `consts.py` | every decimal constant above, high precision | — | all rounded the safe way |
| `weyl` | Lemma 1.2: log, cubic, quadratic and √x phases; N ≤ 1500; q ∈ {1, 2, 3, N, N+5} plus random q ≤ 3N | 2788 | 1.000000000000 (equality at q = 1, \|S\| = N) |
| `weyl` | r-sum inequality, 15 exponents used, q ≤ 10⁵ | — | 0.9999999999 |
| `lem21` | Lemma 2.1 (not in R438): log, quadratic and cubic phases; α pinned at ±½ (the G2 corner); Gauss–Legendre integrals | 238 | 0.191 (at α = −½) |
| `poissoncorner` | Lemma 2.3 at the G5 configuration (ν₁−α, β−ν_N → ½⁺; α, β next to integers), λ₂ ∈ [0.003, 2] | 2992 | 0.195 |
| `phi5` | Lemma 1.4, k = 5, on the actual phase φ_r of Lemma 3.4 (40 digits, 2400 ≤ N ≤ 12000), 30 η ∈ [0.01, 6.5]; the second set is in Yang's main case for the largest η | 36 + 36 | 0.012, 0.016 |
| `block` | the certificate's per-block inequality (l. 826–831 + repairs 3a/3b, constants from certify.py) against exact block sums, t ≤ 3·10⁹, min over q and η | 151 | 0.046 (\|S\|²/RHS) |

There were no violations, so nothing needed a ball-arithmetic recheck. The bounds are loose by large factors, so these tests catch only gross errors; the line-by-line checks of §2 carry the audit.

## 4. Literature
- **2302.13444:**
  - arXiv: still v1 (27 Feb 2023); published as J. Number Theory 262 (2024), doi 10.1016/j.jnt.2024.03.011. The journal text was not fetched (the publisher returns 403 to scripts).
  - OpenAlex: 8 citing works (2212.06867 RNT, 2412.15470 MCOM, JMAA 127426, JMAA 2026 Phragmén–Lindelöf, JNT 2024 log-free density, PAMS Deuring–Heilbronn, 2609.22624, a 2026 Springer chapter). None reports an error; no corrigendum was found.
- **Yang 2301.03165:** v2 = JMAA 534 (2024) 128124. OpenAlex lists 17 citing works; none is an erratum.
- **Pre-emption:**
  - `lst.py 'sub-?Weyl|27/164|exponent pair|van der Corput|explicit.*zeta|…' math.NT math.CA --zenodo '"sub-Weyl"'` (2026-10, 2026-11) gives 2 title hits, both on the van der Corput property of sets (2610.07765, 2611.05267). The Zenodo hits are unrelated.
  - An arXiv search for "explicit sub-Weyl 27/164" finds only 2302.13444. No newer explicit 27/164 constant exists, so 54.6 is not pre-empted.
- **Related work:**
  - Arias de Reyna 2407.02094 gives an independent explicit d-th derivative test of another shape; it is not an erratum.
  - 2609.00537 (DKQ) refines Lemma 2.3 (NOGO l. 529: no gain for us).
  - **Side remark:** DKQ's (S3-bound) (DKQ l. 1317–1320) reuses PY's G5 step and sharpens the constant 2 to 1.251. That is not valid:
    - In PY's ν-range, β−α → 1⁺ gives \|S₃\| → 4/π = 1.2732 > 1.251.
    - In DKQ's own range (0 < f′(b) < 1, ν = 1 … ⌊f′(a)⌋−1), the ν = 1 term is 1/(π(1−f′(b))). At f′(b) = 0.99, f′(a) = 5.5 this gives S₃ = 32.91 > 2.21.
    - Whether DKQ's Corollary itself survives was not checked; it is outside our chain.
- **Patel's thesis (Ohio State 2021):** no open copy was found. Lemma 1.3 is not used, and the r-sum inequality it is cited for is proved directly above.

## 5. Dependency map
| Input | Used by | If it failed | What survives |
|---|---|---|---|
| Lemma 1.2 (Weyl; elementary Cauchy–Schwarz) | 54.6 (Parts C, T); PY 66.7; also HPY's and our 0.568's own A-process steps | every van der Corput bound falls | convexity-type bounds only |
| Lemma 1.4 (Yang) | 54.6: ABA³B blocks (A₅, B₅ ≈ 94 % of the bound), 4th-derivative blocks (A₄, B₄); PY 66.7; Yang's critical-strip bounds | 54.6 and 66.7 fall | HPY 0.618 t^{1/6}log t and our 0.568 (their own second- and third-derivative tests; RevAudit: 0.568 holds with HPY's lemmas, max 0.56797) |
| Lemma 1.4 *reading* | 54.6 at log t ∈ [638, 650] | the other reading gives 55.12 at width-2 intervals | a constant ≈ 55.1–55.2 (not certified) |
| Lemmas 2.1–2.3 | 54.6 (c4, c5a, c5b, E6: E-terms < 10⁻¹⁷ of the bound); 66.7; DKQ | the E-terms change; a gross failure of the B-process would void 54.6 and 66.7 | HPY 0.618 and our 0.568 (no B-process) |
| Lemma 3.1 (Hiary/Gabcke) | every Riemann–Siegel-based explicit bound: 54.6, 66.7, 0.568, HPY | all fall | only bounds from other approximate functional equations |
| Lemma 3.5 μ₂ (error) | PY's 66.7 only | already refuted (REFUTATIONS l. 28) | 54.6 implies 66.7 |
| HPY Thm 1.1 (not audited here) | Part H, t ∈ [3, e^{75}] (53.983) | Part H must be redone | Parts C and T |
Entries that rest on these inputs: RECORDS l. 12 (54.6, `discoveries/sub-weyl-constant.md`), REFUTATIONS l. 28 (μ₂; arithmetic, independent of the lemmas), TREE N9, NOGO l. 529. RECORDS l. 11 (0.568) shares only the Riemann–Siegel input (Lemma 3.1) with this chain.

**Residual, not re-derived:**
- Kershner's sharp constant for general f (the quadratic value was recomputed);
- Rogers' first-derivative bound and the generalised Kuzmin–Landau lemma (classical and refereed);
- the numerics behind Gabcke's 0.127 (RevAudit tested it at 0.12229);
- HPY Thm 1.1;
- the JNT text.

## Verdict (plain words)
The paper's tools are sound. Every estimate our sub-Weyl record borrows from Patel and Yang, and from Yang's earlier paper, was checked step by step and holds. The paper has more typos than usual, and in five places a proof step is wrong or incomplete as written, but each one has a short fix that changes none of the numbers we use. Its decimal constants are all rounded the safe way, and about six thousand new numerical tests found no case where a lemma fails. One point needs care: the paper prints an ambiguous formula for one constant of the derivative test. Our certificate uses the reading that Yang's own proof gives; the other reading would push our 54.6 to about 55.1. The one real error in the paper is still the factor of about 1.10 in its last lemma, which we found before.

## Files
- `consts.py`: high-precision checks of every decimal constant. Log: `logs/consts.log`.
- `falsify.py`: the tests weyl, lem21, phi5, block and poissoncorner. Logs: `logs/{weyl,lem21,phi5,phi5_main,block,poissoncorner}.log`.
- `akbk_check.py`: A_k, B_k from Yang's proof, compared with certify.py, plus both misreadings. Log: `logs/akbk.log`.
- `reading_check.py`: the certificate under the other reading. Log: `logs/reading.log`.
- `s3_examples.py`: the G5 counterexamples and the DKQ remark. Log: `logs/s3_examples.log`.
