AUDIT INCOMPLETE: nothing failed; 2512.15709 (ψ) checks out; open for 2511.14736 (Mertens): Platt's private residue sums above height 1.9·10^6 (they set only the √x constant 11.39) and a line-by-line re-derivation of its technical lemmas (§§2–5, §7; tested end to end).

# Referee audit of Chirre–Helfgott 2511.14736v1 and 2512.15709v1 (node N10, load-bearing inputs)

Agent chaudit, 11 Nov 2026, 14:17–15:31 UTC. Read-only audit: TeX e-prints (line numbers below are those of the arXiv e-print `Moebius.tex`, `VonMangoldt.tex`; quotes checked against the TeX, not pdftotext). Scripts and logs in this directory (`logs/`). No Lean built, nothing committed.

## 0. Answer first

| question | answer |
|---|---|
| Is either bound false? | No violation anywhere: exact M, m, ψ to 10^9, published M(10^k) for k ≤ 25, and explicit-formula tests at T = 10^3, 10^4, 10^5 that a wrong formula fails by 4–75×. |
| Are the numbers right? | Every checkable constant reproduced: Platt's sums at T = 10^6 + 1 and T = 1893193.5 to 10 digits each (3.5 million zeros), the 1/ζ minima on Im s = T to 9 digits, C_T = 113.66888, the zero-sum constants of Prop. vihuela, the four extrema of Lemma pernic, all seven (c1, c2) pairs of Lemma edesmo (exactly). |
| Is the method right? | The two extremal approximants do what is claimed (L1 errors, one-sidedness, decay, numerically). The analytic core of 2512.15709 (contour shift, Graham–Vaaler specialisation, main proposition, Cor. 1.2(a)) has an independent sorry-free Lean formalisation (PNT+ `IEANTN/CH2` + Daugherty–Ryan–Ward, Zenodo 23501708), for 10^7 ≤ T ≤ 1.1·10^7. Nothing comparable exists for 2511.14736. |
| What is open? | 2511.14736 only: Platt's sums at T = 10^7 + 1 … 10^10 + 1 are private computations (our extrapolation agrees to 0.04–0.11 %; the constant 11.39 needs 5·10^−5); its technical lemmas (TV norms, contour shift, square-free R(x) propositions) were read and tested end to end, not re-derived line by line. |
| Do our records change? | No. Fallbacks are in §5 if either paper ever fails. |
| Misprints | 6, none load-bearing (§3). |

## 1. arXiv 2511.14736v1 "Optimal bounds for sums of bounded arithmetic functions" (18 Nov 2025, 44 pp.)

### 1.1 Statements we use
- Thm. mainthmA (l.246–289): A(s) = Σ a_n n^−s meromorphic on ℂ, a_∞ = sup|a_n| < ∞, T ≥ 4π, A(s)T^s bounded on a ladder S. For σ ∈ ℝ, x > e²T: x^{σ−1}Σ_{n≤x} a_n n^−σ = a_∞ι + δ Σ_{ρ∈Z_A(T)∪{σ}} Res_ρ(w_{δ,σ}(s)A(s)x^{s−1}) + ε, δ = π/2T, |ι| ≤ δ tanhc((σ−1)δ), w_{δ,σ}(s) = coth(δ(s−σ)) − tanh(δ(1−σ)), |ε| ≤ (π/4)(a_∞/L + a_∞/L² + I)/T² + 2a_∞/x.
- Cor. mertens (l.342–380): zeros simple to T, σ ≥ −1, x ≥ e²T, max_{r≤1} 1/|ζ(r ± iT)| ≤ log²x ⇒ |Σ_{n≤x} μ(n)n^−σ − (δ Σ_{|γ|≤T} w_{δ,σ}(ρ)x^{ρ−σ}/ζ′(ρ) + 1/ζ(σ))| ≤ (π/2)x^{1−σ}/(T−1) + 2x^−σ.
- Cor. mertensimple (l.390–396): |M(x)| ≤ 3x/(π·10^10) + 11.39√x, |m(x)| ≤ 3/(π·10^10) + 11.39/√x, x ≥ 1. Our MertensCH uses the M part (x ≤ 10^242.85 … 10^3840).

### 1.2 Proof check (chain of Cor. mertensimple)
1. x < 3.61·10^17: eq. garadiol of Cor. mertensimplon (l.3484–3487, T = 10^9, C₂ = 9.758736 = 2·4.87936778767100 rounded up) plus the arithmetic at l.4814–4817: rechecked (1.4753·10^−9 < 2.7149·10^−9). ✓
2. x ≥ 3.61·10^17: Cor. mertosimp (l.4791–4807) = Cor. metamert (l.4744–4789) with T = 10^10 + 1, c = 0.0134: C = 2·5.675256… ≤ 11.350514 plus 2.9c ≤ 0.0389 gives 11.389414; the slack 0.000586√x absorbs |1/ζ(0)| = 2 for x ≥ 1.2·10^7. e³(cT)² = 3.607·10^17 ≤ 3.61·10^17. |w_{δ,0}| ≤ |coth δs| on Re s = ½ via eq. nipaul (l.3445–3447, identity checked), so no tanh correction is needed. ✓
3. Prop. mainsqfprop (l.4652–4735) ← Prop. summsec3 (l.1363–1385) ← Prop. sumwizsqf (l.1094–1177) ← Lemma thirno (l.1052–1086): read line by line; the factor 6/π² in the L1 term (hence 3/(πT) = (6/π²)·π/2T, the 1/ζ(2) gain) follows from Stieltjes integration against Q = 6x/π² + R; the identity d|e^{−πy/T}f| = … (l.1132–1137) checked. ✓
4. Thm. mainthmA ← Prop. gendau (l.2789–2889): residue normalisation re-derived (2π/T)·Res(Φ·A x^s) = δ sgn(1−σ) Res(w A x^s) with Φ from eq. girasol3; L1 term (2π/T)·tanh(λ/4)/λ = tanh(δ(σ−1))/(σ−1). ✓ Lemmas basicfour/edge (l.713–801), sumwiz (l.869–951), adartail (l.851–866): read, correct.
5. Extremal input Prop. carlitfou (l.1996–2093): only the attained error of the given φ_λ matters for correctness (optimality does not). Numerically (extremal_check.py): ‖φ̂_λ − I_λ‖₁ on [−150,150] = tanh(λ/4)/λ − (1.7·10^−4 … 4·10^−5), within the tail ≤ 2.7·10^−4, for λ = 0.5, 2, 5, −1, −0.01; sup_{|y|≥20} y²|φ̂_λ − I_λ| ≤ 0.01989 = 1/16π (Lemma arborio's κ, sharp as λ → 0). ✓
6. Square-free constant c = 0.0134 (Cor. bataclar l.4634–4639 ← Lemmas nopgik l.4589–4605, gopnik l.4607–4633, MOST's |R(x)| ≤ 1.12543x^{1/4} to 10^18): no circularity (both lemmas use Cor. mertensimplon, not mertosimp). All coefficients of nopgik/gopnik recomputed from Props. qmediano/gould (0.772103, 0.860845, 8.98413, 0.999; 6.51·10^−11, 2.5222, 1.042, 1.086, 1.99·10^−5, 6.01·10^7, final 0.84491, 0.00161, 3.9·10^−5). The (c1, c2) pairs of Lemma edesmo are exactly the optimal ones (edesmo.c, full periods up to (30030)²). Props. gould/qmediano themselves (l.4356–4583) read, not re-derived. Sensitivity: any c ≤ 0.0136 keeps 11.39 (the garadiol branch covers x ≤ 1.2·10^18, so the hand-over e³(cT)² may move).
7. Not re-derived line by line: Lemmas artanor (TV constants c1 < 17/5, c2 < 2/9), kolobrz/Prop. shiftbounded (contour shift, the π/8T² term), paulin/pommedupe (trivial zeros), zetinvbound. They feed only the x/T², O(a_∞) and x^−2 terms; the paper fits these into the gap between π/2T and π/2(T−1) ((7/4)(π/4) + 0.45e^−4 = 1.383 < π/2), and for M(x) any O(1) excess is absorbed by 11.39√x above 10^16. The end-to-end tests (1.4) exercise them.

### 1.3 Computational inputs and our recomputation
| input | source | ours |
|---|---|---|
| δΣ_{0<γ<T}\|coth δρ\|/\|ζ′(ρ)\| (Lemma platt109, l.3285–3319), T = 10^6 + 1: 2.66161277991001 | D. Platt, private, interval arithmetic | 2.661612780157 (Odlyzko's 2,001,052 zeros; ζ′ by Arb below 50,000th zero, vectorised Riemann–Siegel above, max rel. dev. from Arb 8·10^−8) |
| same, T = 1893193.5: 2.85417779533422 | Platt | 2.854177795581 (Odlyzko's zeros plus 1,498,948 found by Riemann–Siegel sign changes on a 0.05 grid with a close-pair search, 307 pairs inside one cell; total exactly 3,500,000 = N(1893193.5) as quoted in Remark (a) to lem:chamlu; positions and \|ζ′\| agree with Arb to 10^−9 and 10^−8 at five sampled indices) |
| same, T = 10^7+1, 10^8+1, 10^9, 10^10+1: 3.36904620179490, 4.10963816503581, 4.87936778767100, 5.675256 | Platt | not reachable (10^10 needs 3·10^10 zeros). Extrapolation from our zeros with J(t) = Σ_{γ≤t}1/\|ζ′\| ≈ (t/2π)(K₁L^{1/4} + K₂L^{−3/4}) (HKO shape, J/((t/2π)L^{1/4}) = 0.959 … 0.998 for t = 10^4 … 1.1·10^6): fit on [10^4, 2·10^5]: −0.007, +0.008, +0.024, +0.040 % at 10^7+1, 10^8+1, 10^9, 10^10+1 (and −0.023 % against our exact 10^6+1 value); fit on [2·10^5, 1.13·10^6]: +0.034, +0.062, +0.088, +0.112 %; one-term fit −0.09 … −1.4 %. Plausible; not a proof. |
| max_{−1/64≤σ≤1} 1/\|ζ(σ+iT)\| (Lemma chamlu, l.3220–3268) | authors, Arb/MPFI bisection | grid of 1200 σ + refinement (chamlu.py): 0.591342108, 0.551087905, 0.579764045, 0.669256578, 0.536165862 (at σ = 0.52029) for T = 1893193.5, 10^7+1, 10^8+1, 10^9, 10^10+1 = paper to 9 digits; row T = 10^6+1 (0.894198297) is the value at T = 10^6 (true max at 10^6+1: 0.835322), a mislabel in the safe direction; σ < −1/64 samples ≤ 0.0016 |
| \|M(x)\| ≤ 0.570591√x, 33 ≤ x ≤ 10^16 (Lemma hurst) | Hurst, Math. Comp. 2018 | max \|M(x)\|/√x = 0.567105 (x = 199) on [33, 10^9] |
| \|m(x)\| ≤ 0.569449/√x, 3 ≤ x ≤ 10^14 | Helfgott's book, code on request | max \|m(x)\|√x = 0.461713 on [3, 10^9] |
| RH and simple zeros to 3·10^12 | Platt–Trudgian, Bull. LMS 2021 | refereed; not rechecked |

Where 11.35051 comes from: 2 × Platt's one-sided sum at T = 10^10 + 1 (about 3.2·10^10 zeros, each with ζ′(ρ)), computed by D. Platt in interval arithmetic and reported to 7 digits (lower T to 15 digits); unpublished, no data or code released. There is no tail: the weight vanishes for |γ| > T, so the whole sum is Platt's computation and its rigour is his. Margin: 11.39 fails only if Platt's one-sided sum at 10^10 + 1 exceeds 5.675256 by 0.00029 (5·10^−5 relative). The leading term 3/(π·10^10) does not depend on Platt's sums at all.

### 1.4 Sanity and falsification tests
- Exact M(n), m(n) for n ≤ 10^9 (sieve.c): max |M(x)|/(3x/(π10^10) + 11.39√x) = 0.0878, max |m(x)|/(…) = 0.124 (both at x = 1); M(10^9) = −222. OEIS A084237 (M(10^k), k ≤ 23) and Hurst 2607.07566 (M(10^24) = 7189337839, M(10^25) = −258560632948): ratios ≤ 0.030 for k ≥ 1 (k = 10). Classical bounds give no conflict: the leading term is of the order 1/T expected from zeros up to T (Prop. counterex, l.1569, shows π/2T is attained for some A(s)); it beats Ramaré's 0.013x/log x on about [3·10^8, e^(1.4·10^8)] and Cohen–Dress–El Marraki's x/4345 from 2.5·10^9 on.
- Explicit-formula test of Cor. mertens (eft.py) with Odlyzko zeros, complex ζ′(ρ), exact M, m on 3·10^4–4.6·10^4 points of [e²T, 10^9] (x = n and x → (n+1)^−):

| T | zeros | max \|M − P₀\| / bound | max \|M − P₀\| / δx | max \|m − P₁\| / bound | rms(M − P₀)/rms(M) where δx < √x | same with zero sum of wrong sign |
|---|---|---|---|---|---|---|
| 10^3 | 649 | 0.391 | 0.452 | 0.391 | 9.6 % | fails, 4.6× bound |
| 10^4 | 10,142 | 0.319 | 0.369 | 0.319 | 2.8 % | fails, 14× |
| 10^5 | 138,069 | 0.328 | 0.380 | 0.328 | 1.8 % | fails, 45× |

(condition max 1/|ζ(r ± iT)| ≤ log²x holds: 1.16, 3.02, 0.52 against log²x ≥ 79.) The formula reproduces M(x) to 2–10 % and never leaves |ι| ≤ δ.

## 2. arXiv 2512.15709v1 "Optimal bounds for sums of non-negative arithmetic functions" (17 Dec 2025, 52 pp.)

### 2.1 Statements we use
- Thm. mainthmB (l.313–353): a_n ≥ 0, A(s) simple pole at 1, bounded on a ladder and an admissible contour C; for x > T: explicit formula with weights (π/T)coth(π(s−σ)/T) on real poles, (2π/T) Im Σ ω⁺_{T,σ}(ρ)x^{ρ−1}Res, error O*((π/T)Σ_real Res + 2I_{+,C}/T² + (2π/T) Re Σ θ_{T,1}(ρ)x^{ρ−1}Res).
- Cor. psi (l.389–397): RH to T ≥ 10^7, x > max(T, 10^9): |ψ(x) − x(π/T)coth(π/T)| ≤ πx/(T−1) + C_T√x, C_T = (1/2π)log²(T/2π) − (1/6π)log(T/2π).
- Cor. psiexpl (l.402–406): |ψ(x) − x| ≤ πx/(3·10^12) + 113.67√x for x ≥ 1. Used by RobinCH, OddRobin, RamRepair (§4).

### 2.2 Proof check (chain of Cor. psiexpl)
1. Proof l.4981–4999: x < 3·10^12 + 3 by Lemma pernic; else Cor. psi with T = 3·10^12 + 1 + π/3 (RH verified to 3,000,175,332,800); π²/3T² + π/(T−1) ≤ π/(T−1−π/3) checked; C_T = 113.668881 ≤ 113.67. ✓
2. Prop. sagaro (l.4868–4949): T′ ∈ [T − ½, T] from Lemma hardin; small terms 0.0001 + 0.0098 < 0.01 of the 1.01/6π reserve; coth-shift bounds; π/2T² + π/(T − ½) < π/(T−1). Arithmetic rechecked. ✓
3. Prop. vihuela (l.3948–4007): the combination |Im Σa| + |Re Σb| ≤ max_ξ Σ|ω⁺ + ξiθ| x^{−1/2} (l.3400–3402) is valid. Lemma adar (l.3801–3918), salmon (l.3921–3940), thonny (l.3451–3544), tritura (l.3550–3640), Lehmanmodern (l.5800–5856), oldie (l.5687–5706), brut, cathinv: read line by line, correct. Constants recomputed (vihuela_check.py, zsums.py): C₁ = 0.168938533, C₂ = 0.164184385 (2C₁ = 0.337877 ≥ 0.337876, 2(C₁ − C₂) = 0.0095083 ≥ 0.0095); err₁/t₀ = 8.161·10^−4 + 1.4433·10^−4|σ−½|; err₂/T = 1.404·10^−6 + 3.97·10^−6|σ−½|; 2Σ_{γ≤2·10^4}1/γ = 10.3193172 (= (1/2π)log²(t₀/2π) − 0.0343588, paper "0.03435…" truncated). Assembled bound vs claim: 31.6688 ≤ 31.6904 (T = 10^7), 55.7135 ≤ 55.7519 (10^9), 113.6144 ≤ 113.6546 (3·10^12). Tritura's integral inequality holds with slack 2.1–4.0. ✓
4. Actual zero sums (2π/T)Σ_{0<γ≤T}|ω⁺_{T,σ} + ξiθ_{T,1}| for σ ∈ {0,1}, ξ = ±1 at T = 10^3, 10^4, 10^5, 10^6, 1.13·10^6 lie below the vihuela bound by 0.150, 0.210, 0.267, 0.323, 0.326 (outside its range T ≥ 10^7; the margin grows). ✓
5. Lemma hardin (l.4802–4838) and Prop. coronidis (l.4734–4793): assembly rechecked (39/5 + 68/log T < 51/4 and 37 + 311/log T < 60 for T ≥ 10^6; coronidis' final inequality at x = 15: −0.061 < 0). Integral lemmas (adiaro, moruno, ranadi, arles, adioso) read, not re-derived; formalised for T ≤ 1.1·10^7 by DRW.
6. Graham–Vaaler input Prop. grahamvaaler (l.1645–1693), numerically for ν = 0.5, 2, 5: F⁺ ≥ I_ν and F⁻ ≤ I_ν everywhere on the grid (to 10^−13), F⁺(0) = 1, F⁻(0) = 0, L1 errors equal 1/(1−e^−ν) − 1/ν and 1/ν − 1/(e^ν − 1) within the tail. ✓
7. Lemma pernic (l.4955–4980) reproduced to 10^9 (sieve2.c, sieve3.c): max (ψ(x)−x)/√x = 0.7905927544 at 110102617, min on [10^4, 10^9] = −0.7509024437 at 36917099^−, and for √x·(Σ_{n≤x}Λ(n)/n − log x + γ): 0.7878013324 at 110102617, −0.7585825520 at 1423^−.
8. Explicit-formula test of Thm. mainthmB for A = −ζ′/ζ, σ = 0 (eft.py): ψ(x)/x stays inside the sandwich main ± W, W = π/T − (2π/T)Re Σθ_{T,1}x^{ρ−1}: max |ψ/x − main|/W = 0.586, 0.686, 0.736 for T = 10^3, 10^4, 10^5; the error reaches 0.58, 0.82, 0.97 of π/T (the claimed optimal constant shows); wrong-sign control fails 4.1×, 20×, 75×.
9. Sanity: max |ψ(x) − x|/(πx/(3·10^12) + 113.67√x) = 0.0124 on [1, 10^9]; Büthe's refereed |ψ(x) − x| ≤ 0.94√x (11 < x ≤ 10^19, Math. Comp. 2018) is 120× inside the bound.
10. Where the constants come from: C_T is analytic (N(T) via Rosser's |Q(t)| ≤ 0.137 log t + 0.443 log log t + 1.588); the only data are RH to 3·10^12 (Platt–Trudgian), the zeros below 2·10^4 (we reproduced 2Σ1/γ) and Lemma pernic's sieve (reproduced to 10^9).

### 2.3 Independent formalisation
PrimeNumberTheoremAnd (`IEANTN/CH2/CH2.lean`, `CH2_part1.lean`, 10,810 lines) formalises Lemma 5.1 (contour shift), Prop. 5.2 (Graham–Vaaler specialisation), Lemma 4.3, Props. 2.3–2.4 of this paper; Daugherty–Ryan–Ward (Zenodo 23501708, 20 Oct 2026, `CHRoute`, 17,395 lines) drop the sorried upstream corollaries, replace the sorried decay lemmas by their own, add Prop. 6.1 (σ = 0) and the ζ-specific lemmas, and prove `cor_1_2_a_explicit`: for 10^7 ≤ T ≤ 1.1·10^7, x ≥ 10^9, RH to T, no zeros with |Im s| ≤ 2, a zero gap near T and a numerical zero-weight inequality (margin 1/10, discharged by FLINT for T = 10^5, 10^7), |ψ(x) − x(π/T)coth(π/T)| ≤ πx/(T−1) + C_T√x. Their axiom audit (237 CH2 lemmas) lists only propext, Classical.choice, Quot.sound. We did not rebuild it (no Lean in this audit). It replaces Prop. vihuela by computation, so the T = 3·10^12 instance still rests on vihuela, checked in 2.2.3.

## 3. Misprints and small gaps (none load-bearing)
1. 2512.15709 Lemma pernic eq. malina (l.4964–4967): the bounds −0.7585825520 ≤ … ≤ 0.787 hold for √x·(Σ_{n≤x}Λ(n)/n − log x + γ); the factor √x is missing. Its range "1 ≤ x ≤ 10^12 + 3" also falls short of the 3·10^12 + 3 needed in the proof of Cor. psiexpl (second claim); the proof text (l.4979) says the computation covers [10^10, 3·10^12 + 3]. The ψ part (to 10^13) is complete.
2. 2512.15709 Prop. coronidis (l.4792): "the left side increases for x ≥ e" should read "decreases"; as written the argument would not extend from x = 15.
3. 2512.15709 eq. espresso (l.3863) and l.3867: "=" for "≤".
4. 2511.14736 Lemma chamlu: the row T = 10^6 + 1 gives the maximum for T = 10^6 (0.894198296); at 10^6 + 1 the maximum is 0.835322, so the stated bound stays valid.
5. 2511.14736 Lemma edesmo proof (l.4342–4350): ν is defined without n (should be max |Q_q(n) − c₁n|), and the listed ν equal c₂ for q ≤ 5 but c₂/2 for q ≥ 7. The pairs in the statement are exactly right.
6. 2511.14736 Cor. mertensimplon proof (l.3513–3566): the step from Cor. mertens to "+4" sits in a comment block in the e-print; it is the obvious |1/ζ(0)| + 2.

## 4. Literature
- arXiv abs pages (11 Nov 2026): both still v1 only, no journal reference, no comments beyond page counts. arXiv API search (Chirre AND Helfgott): only these two plus 2005.03162. No erratum found.
- Citations (OpenAlex): 2511.14736 (W4416561291): Lee–Leong 2208.06141 (later version), Zuniga-Alterman 2610.17720. 2512.15709 (W7116070227): Fiori 2609.22624; Daugherty–Ryan–Ward Zenodo 23475183 and 23501708 (the formalisation above; the only independent check found).
- Code: none public. Platt's residue sums and the 1/ζ minima are private; Helfgott's m(x) code "available upon request"; sagetex is loaded but unused.

## 5. Dependency map
| our entry | uses | if 2511.14736 fails | if 2512.15709 fails |
|---|---|---|---|
| MertensCH record \|M(x)\| ≤ 9.77·10^4x^0.9382, 1.54·10^4x^0.9409 (RECORDS N10, discoveries/mertens-power-saving.md, TREE N10) | Cor. mertensimple for x ≤ X₀ (leading term 3/(π10^10) binding there; 11.39 never binding) | full failure: back to Mertens78's 1.84·10^6x^0.9382, 7.91·10^5x^0.9409; if only the square-free step fails, Cor. mertensimplon (πx/(2·10^10) + 11.350514√x, needs only Thm. mainthmA + Platt + Hurst) multiplies every constant by at most π²/6 ≈ 1.645, θ unchanged (estimate, not re-certified) | – |
| RobinCH: Robin to 10^(10^46.79), Nicolas to p_k ≤ 1.78·10^48; 10^(10^46.44) without half-plane (RECORDS N4, discoveries/robin-nicolas-ranges.md, TREE N4, R353, R402) | Cor. psiexpl above 10^19 (both πx/3·10^12 and 113.67√x enter) | – | back to RobinHeight's 10^(10^38.71) (own pointwise Lemma 10_k) |
| OddRobin range above largest prime 5·10^8 (discoveries/odd-robin.md) | RobinCH | – | range falls to what RobinHeight's inputs give (not recomputed); the criterion itself and the range to e^(5.0·10^8) need neither paper |
| RamRepair Thm A to e^753.98 (NOVELTY.md) | 2512.15709 | – | e^672.84 with refereed inputs only |
| Buthe78 "NO RECORD" verdict, NOGO line on buthe78, TREE N10 "Büthe-type ranges to 1.6·10^27" | 2512.15709 as the comparator | – | verdict reopens: Büthe–Johnston + 7/8 (ψ to 5.19·10^26) would need re-evaluation |
| NOGO R307 (CH + Fiori–Jaskari), R379, R389, R428, Harvest3/5 | quoted only | – | – |

Dependents of one input count once: every N4 range above 10^19 rests on Cor. psiexpl alone; the N10 Mertens record rests on Cor. mertensimple alone.

## 6. Referee verdict in plain words
Both papers hold up. Nothing we could compute contradicts them, and every number we could recompute came out the same, twice to ten digits. The prime-counting paper (2512.15709) is in the best shape: its core argument has been checked by a computer proof assistant by another group, and the one step that proof leaves to hand calculation we redid line by line. The Möbius paper (2511.14736) is consistent everywhere we looked, but two things stay on trust: the sums over zeta zeros up to height ten billion, which only D. Platt has computed, and the long chain of technical lemmas, which we tested as a whole rather than step by step. Our records that use these papers stand; if either paper ever breaks, the fallbacks in §5 apply. Both are still unrefereed preprints and should be cited as such.

## 7. Suggested index line
chaudit (N10, audit): Chirre–Helfgott 2511.14736 and 2512.15709 refereed by us: no error; all checkable constants reproduced (Platt's sums at 10^6 + 1 and 1893193.5 to 10 digits); ψ paper's core Lean-checked by PNT+/DRW; open: Platt's private sums above 1.9·10^6; 6 harmless misprints.

## 8. Files
`zp_rs.py` (ζ′ at 2,001,052 zeros), `zsums.py` (Platt and vihuela sums), `platt_extrap.py`, `zeros2.py` (zeros to 1893193.5), `chamlu.py`, `eft.py` (explicit-formula tests), `extremal_check.py`, `vihuela_check.py`, `sieve.c`, `sieve2.c`, `sieve3.c`, `edesmo.c`, `misc_checks.py`, `zsum2.py`; outputs in `logs/`. Inputs deleted after use: the two e-prints, Odlyzko's `zeros6.gz` (https://www-users.cse.umn.edu/~odlyzko/zeta_tables/), the DRW release; to rerun, download `zeros6.gz` to `.scratch/chaudit/dl/`.
