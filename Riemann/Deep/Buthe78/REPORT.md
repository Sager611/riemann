NO RECORD: Chirre–Helfgott (arXiv 2512.15709, Dec 2025, Cor. "psiexpl") prove |ψ(x) − x| ≤ πx/(3·10^12) + 113.67√x for all x ≥ 1, which already gives Büthe's bounds up to about 1.6·10^27 (π(x) up to about 1.1·10^27) and never uses zeros above height 3·10^12, so the 7/8 half-plane cannot improve it. The 7/8 half-plane does enlarge the Büthe–Johnston method 1.30-fold (certified here: ψ to 5.19·10^26, π to 3.68·10^26; Johnston 2022 stated 1.101·10^26, not 2.169·10^25 as the brief assumed).

# Büthe's partial-RH range and the 7/8 half-plane (node N10)

Agent buthe78, 2026-11-05. Files: this report, `verify.py` (ball-arithmetic certificate, 5–10 s).

## 1. Short answer

1. **The brief's baseline was out of date.** Büthe's condition 4.92·√(x/log x) ≤ T with T = 3·10^12 gives 2.169·10^25, but Johnston (arXiv 2109.02249, Ramanujan J. 2022, Thm 3.2 and Cor 3.3) already proved all four bounds (ψ, θ, Π, π) for x ≤ 1.101·10^26.
2. **The 7/8 half-plane helps Büthe's method by a factor of 1.30 in x.** Büthe's smoothing kernel is the Logan function, which is not band-limited. Zeros above c/ε enter through the tail (e:rem1) with weight x^{Re ρ}/sinh c. With Re ρ ≤ 7/8 the kernel parameter can drop from c ≈ ½ log x + D to c ≈ ⅜ log x + D. That shrinks the |ψ − ψ_{c,ε}| term, which grows like √c, at fixed height c/ε = T.
3. **Re-optimizing Büthe's own tools helps more than 7/8 does.** The two tools are a one-sided smoothing shift α > 0 in his Prop. 4 and a stepwise zero sum from his Lemma 3. Büthe uses both in his Theorem 1, but neither he nor Johnston used them in the √x-range theorem. Without 7/8 they give a factor of about 2.5 (π: 2.82·10^26); with 7/8 the range reaches 3.68·10^26 for π.
4. **None of this is a record.** Chirre–Helfgott use band-limited extremal weights (Graham–Vaaler / Beurling–Selberg) and get error πx/T. Their explicit corollary implies ranges about 3× larger than our best (table below). Their method sees no zeros above T, so 7/8 adds nothing there.

## 2. Ranges (upper end X of "bound holds for lower limit ≤ x ≤ X"), T = 3·10^12

| Method | uses 7/8 | ψ (√x log²x/8π) | θ (same) | Π (√x log x/8π) | π (same) | status |
|---|---|---|---|---|---|---|
| Büthe 2016, Thm 2 (4.92 condition) | no | 2.169e25 | 2.169e25 | 2.169e25 | 2.169e25 | published |
| Johnston 2022, Cor 3.3 | no | 1.101e26 | 1.101e26 | 1.101e26 | 1.101e26 | published (proof has a slip, claim true: see §6) |
| B0: Büthe's split, parameters re-optimized | no | 1.99e26 | 1.86e26 | 1.40e26 | 1.28e26 | certified here |
| M0: B0 with the 7/8 tail (the R295 idea) | yes | 2.60e26 | 2.43e26 | 1.82e26 | 1.67e26 | certified here |
| M1: α > 0 and stepwise zero sum | no | 3.98e26 | 3.77e26 | 3.01e26 | 2.82e26 | certified here |
| M2: M1 with the 7/8 tail | yes | **5.19e26** | **4.92e26** | **3.92e26** | **3.68e26** | certified here |
| Chirre–Helfgott 2025, Cor. psiexpl (+ ψ−θ bound) | no | ≈1.66e27 | ≈1.56e27 | ≈1.18e27 | ≈1.09e27 | implied (one line, float; not stated in their paper) |
| Chirre–Helfgott 2025, their stated condition 2π²√x/log x ≤ T | no | 8.2e25 | – | – | – | stated |

Gain from 7/8 alone: M0/B0 = 1.30–1.31 and M2/M1 = 1.30–1.31, for all four statements. For π, the brief's estimate (c = ⅜ log x + 5 with Büthe's ε) gives 3.9·10^25. That is correct, but it is below Johnston's value.

Tradeoffs:

| | Büthe-type (Logan kernel, M2) | Chirre–Helfgott (band-limited) |
|---|---|---|
| zeros above T | enter, weight ≈ x^{β}e^{−c}; 7/8 helps | never enter; 7/8 irrelevant |
| smoothing error at x = 3.68e26 | ≈ 0.66·√x log x (≈ 2.0·πx/T) | πx/T (optimal for support [−T, T]) |
| zero sum at x = 3.68e26 | ≈ 1.67·√x log x (smaller: effective height ≈ T/√c) | (1/2π)log²(T/2π)·√x ≈ 1.86·√x log x |
| refereed inputs | yes, apart from 7/8 (Lean-checked here) | preprint v1, not yet refereed |

## 3. Error terms of Büthe's proof (TeX of 1410.7015, `ppnt.tex`, `zsum.tex`, `psum.tex`, `explicit-formula.tex`)

Parameters: c (Logan width), ε (time support), a (split point a·c/ε); Büthe takes c = ½ log x + 5, ε = log^{3/2}x/(8√x), a = √(2/c), α = 0.

| Term | Source | Bound | Dependence |
|---|---|---|---|
| constant | Prop 1 (e:exp-form) and the 0.024 bound | Θ(2): \|C₁\| = 1.861, \|Σ a(ρ)/ρ\| ≤ 0.024, 8ε\|log ε\| | needs ε small; no c dependence; holds for every zero |
| 𝓔₁ tail, \|γ\| > c/ε | Prop 3 (e:rem1), via Lemma 2 | 0.16 (x+1)/sinh c · e^{0.71√(cε)} log(3c) log(c/ε) | ≈ 0.32 e^{−(c−½log x)}√x·log(3c)·½log x; this forces c ≥ ½ log x + O(1); unconditional (Re ρ < 1) |
| 𝓔₂, a c/ε < \|γ\| ≤ c/ε | Prop 3 (e:rem2) | (1+11cε)/(πca²) log(c/ε) cosh(c√(1−a²))/sinh c · √x | ≈ e^{−k}/(2πk)·½ log x·√x with a² = 2k/c; RH needed up to c/ε |
| 𝓔₃, 0 < \|γ\| ≤ a c/ε | Lemma 3 (ii) | (√x/2π) log²(ac/(2πε)) | the main term; RH needed up to a c/ε |
| 𝓔₄, \|ψ − ψ_{c,ε}\| main part | Prop 4 (α = 0), Prop 5 | e^{2ε} log(e^ε x)·2εx\|ν_c(0)\|/log B, \|ν_c(0)\| = I₁(c)/(2 sinh c) ≈ (2πc)^{−1/2} | ∝ εx/√c: grows with ε, shrinks with c |
| 𝓔₅, prime powers | Prop 4 | e^{2ε} log(e^ε x)(2.01 ε√x + ½ log log 2x²) | negligible (≈ 10^{−11}·√x log x) |
| RH height | text before (e:zsum1-bound) | c/ε ≤ T | with Büthe's choice c/ε = (4 + 40/log x)√(x/log x) ≤ 4.92√(x/log x) |

## 4. Modified proof (each change next to the original label)

The hypotheses become: RH for 0 < Im ρ ≤ T = 3·10^12 (Platt–Trudgian, arXiv 2004.09765), and, in M0/M2, ζ(s) ≠ 0 for Re s > 7/8 (OpenAI, Oct 2026; Lean-checked in `Riemann/Deep/QRH`). We prove the strong form |ψ(y) − y| < √y log y (log y − C)/(8π) for 2.169·10^25 ≤ y ≤ X'. Here C = 0 (ψ), C_θ (θ), 2 (Π), 2 + C_θ (π), with C_θ(y) = 8π(y^{1/2} + y^{1/3} + (log y/log 2) y^{1/4})/(√y log y) ≈ 0.42. Below 2.169·10^25 we use Büthe's theorem with T = 3·10^12, as Johnston's Prop 3.1 does.

1. **Parameters** (replaces "c = ½log x + 5, ε = log^{3/2}x/(8√x)"): c = θ₀ log y + D with θ₀ = ⅜ (7/8 used) or ½ (not used), ε = c/T, so c/ε = T exactly. This replaces the requirement "c/ε ≤ 4.92√(x/log x) ≤ T". D, α and k are fixed per statement (verify.py, CASES).
2. **(e:rem1)/(e:zsum1-bound) → 7/8 tail.** Lemma 2 bounds Σ_{\|γ\|>c/ε} \|ℓ(z_ρ)\|/\|γ\| for all zeros. The zeros ρ and 1 − ρ̄ have conjugate arguments z, z̄ and equal \|ℓ\|. Pairing them gives Σ \|a(ρ)x^ρ/ρ\| ≤ max_β ½(x^β + x^{1−β})·0.32(…) = **0.16 (x^{7/8} + x^{1/8})/sinh c · e^{0.71√(cε)} log(3c) log(c/ε)**, since β ↦ x^β + x^{1−β} is convex and β ∈ [1/8, 7/8]. Büthe has x + 1 here. With c = ⅜ log y + D this is ≈ 0.32 e^{−D}√y·log(3c)·log T, the same size as before.
3. **(e:zsum2-bound)+(e:zsum3-bound) → stepwise sum** (M1, M2; M0 keeps Büthe's pair with a² = 2k/c, k ≈ 0.88 instead of 1). All zeros with \|γ\| ≤ T lie on the line, \|a(ρ)\| ≤ ℓ(γ) (since λ = ℓ(i/2) ≥ 1), and ℓ(t) = sinh(c√(1−(t/T)²))/(√(1−(t/T)²) sinh c) decreases on [0, T]. With a grid t₁ = 10^6 < t₂ < … < t_n < T < t_{n+1} = T + 1 (400 geometric blocks), Lemma 3 (ii) and (i) give
   Σ_{0<\|γ\|≤T} \|a(ρ)x^ρ/ρ\| ≤ 2√x [ (1/4π)log²(t₁/2π) + Σ_i ℓ(t_i)((1/4π)(log²(t_{i+1}/2π) − log²(t_i/2π)) + 5 log t_i/t_i) ].
   This is Büthe's piecewise lemma in `chebb.tex`, with the error term taken at the left end t_i, as Lemma 3 states. It needs RH only up to T and no Lemma 5.5 of Büthe 2018.
4. **(e:psum-bound) → Prop 4 with α ∈ (0,1)** (M1, M2). Take x = e^{±αε}y in Prop 4. Then
   |ψ(y) − y| ≤ αε e^{αε} y + |Σ_ρ a(ρ)x^ρ/ρ| + 2 + A(x,c,ε,α),
   with A = e^{2ε} log(e^ε x)[2εx|ν_c(α)|/log B + 2.01ε√x + ½ log log 2x²] and B = εxe^{−ε}|ν_c(α)|/(2μ_c(α)).
   Here μ_c(α) = ∫_α^1 η_c and |ν_c(α)| = ∫_α^1 (u − α)η_c(u)du, with η_c(u) = c I₀(c√(1−u²))/(2 sinh c). Both are computed from the exact series Σ_n (c²/4)^n/(n!)² ∫_α^1 (1−u²)^n(…)du, with a recursion and a tail bound. At α = 0 this gives |ν_c(0)| = I₁(c)/(2 sinh c), as it should.
   The optimum is α ≈ 0.10–0.12. Example at c = 27: |ν_c| falls from 0.0757 to 0.0312, and at the optimum the shift term equals the smoothing term.
5. **ψ → θ** (Büthe: "ψ(x) − ψ(√x) ≤ θ(x) ≤ ψ(x)"; that first inequality is reversed, see §6). We use ψ(y) − θ(y) = Σ_{k≥2} θ(y^{1/k}) < y^{1/2} + y^{1/3} + (log y/log 2)y^{1/4}, from θ(t) < t for t ≤ 1.39·10^17 (Platt–Trudgian 2016). This gives the constant C_θ above.
6. **Partial summation** (unchanged): it needs the strong ψ bound (C = 2) or the strong θ bound (C = 2) on [5000, x]. The constants at 5000, |Π(5000) − li(5000) − (ψ(5000) − 5000)/log 5000| = 0.4786 and |π(5000) − li(5000) − (θ(5000) − 5000)/log 5000| = 4.9130, are both < √5000/(4π) = 5.6270 (recomputed in verify.py).
7. **Side conditions**, checked in every box: ε < 10^{−4} (covers ε < 10^{−3}, 10^{−2} and 1/10, and log x > 2/|log ε|), c > 3, B > 1, a·T > 5000 (M0/B0), x > 100.

Example at the π endpoint of M2 (y = 3.68·10^26, c = 27.14, α = 0.1163), in units of √y log y: target (L − C)/8π = 2.3380; zeros ≤ T 1.6701; tail > T 0.0099; |ψ − ψ_{c,ε}| 0.3280; shift 0.3299; margin 3.2·10^{−5}.

## 5. Certificate

`verify.py` uses python-flint 0.9 (arb balls, 160 bits). It covers log y ∈ [log 2.169·10^25, log X'] by boxes of width 0.02, bisected when needed (minimum width 10^{−7}). It encloses the normalized margin, target minus every term above, with log y, c, ε and the shift δ ∈ [−αε, αε] as balls. ℓ is decreasing in c, so the stepwise sum is evaluated at c_lo. It certifies all 16 cases (B0, M0, M1, M2 × ψ, θ, Π, π) in 4.1–7.6 s single-threaded; every box has a positive lower bound.

Checks run in scratch:
- **Controls.** Just past the float optimum the certificate fails, as it must: M2 π at 3.70e26, M2 ψ at 5.25e26, B0 π at 1.30e26.
- **Cross-check.** Two independent float models (mpmath, unnormalized, besseli) agree with verify.py's margins to 10^{−11}.
- **Parameter search.** Coordinate search over (D, k) or (D, α), with ε = c/T optimal at the end point; scratch files opt2.py and opt3.py, deleted.

## 6. Findings on the literature

- **Johnston 2022 (arXiv 2109.02249 v2), proof of Thm 3.2.** The bound 𝓔₁ ≤ 0.000032√x log x log log x divides by log x where it should divide by log log x: 0.0000316 = 0.16(x+1)/sinh c·log(3c)/(√x log x)·½ at x = 2.169e25, D = 6. The true constant is 0.000436, about 14× larger.
  - His later iterations lower D to 2.34. There the true tail is 0.068·√x log x against the 0.005 he used, and his final parameters give margin −0.021 at 1.101e26. His last iteration therefore fails as written.
  - **The claim itself is true.** Re-optimized parameters (B0: c = ½ log x + 4.79, a² = 1.764/c, ε = c/T) prove all four bounds to 1.28e26 (π) and further, so Cor 3.3 stands. This is a repairable slip, not a refutation.
- **Büthe 2016, end of proof of Thm 2.** "ψ(x) − ψ(√x) ≤ θ(x)" is false; the correct direction is θ(x) ≤ ψ(x) − ψ(√x) (x = 8: ψ(8) − ψ(√8) = 6.04 > θ(8) = 5.35). This is harmless. Büthe 2018 (arXiv 1511.02032) gives |ψ(t) − t| ≤ 0.94√t for t ≤ 10^19, and with ψ − θ < 1.4262√t (Rosser–Schoenfeld) and θ(t) < t that yields the strong θ bound on [5000, 10^19]. The same 0.94√t also replaces Büthe's Theorem 1 tables on [e^18, e^44]. Theorem 1 lacks the log(3c) of Lemma 2 (corrected by Bhattacharjee 2024, see Riemann/Deep/RamRepair/NOVELTY.md); we never use it.
- **Büthe 2016, piecewise lemma in `chebb.tex`.** The error is taken at T₁ (20π log T₁/T₁), but Lemma 3 gives it at T₀. This is minor; verify.py uses T₀.
- **Improved half-planes.** Liu 2610.12234 (0.874957), Xu–Sang 2610.15291 (0.8748867, conditional on OpenAI's inputs), Liu 2610.26835 (0.874887) and Chang–Wang 2611.00924 (0.87446, conditional) would lower c by at most 0.0005·log y ≈ 0.03, a negligible effect. We kept 7/8 (Lean-checked here).
- **Zero-density tail ("optionally").** At height 3·10^12, explicit N(σ,T) bounds (Kadiri–Lumley–Ng type, log^{5−2σ}T factors) beat the trivial count only for σ ≳ 0.85. Our estimate is ≲ 5% more range, and none for the band-limited route. Not done.

## 7. Novelty check

- **arXiv math.NT listings.** 2026-09 (812 entries), 2026-10 (854) and 2026-11 (140 so far), scanned by title keywords. Abstracts read for every hit on partial RH, explicit PNT, half-planes or prime counting: 2609.23222 Fiori–Jaskari (7/8 version of the classical bound, a different statement), 2610.12234, 2610.15291, 2610.26835, 2611.00924, 2610.18491, 2610.39120, 2609.02713. None extends the Büthe/Johnston range.
- **arXiv search** (web_search, arXiv mode) and **arXiv author search** (Chirre, Helfgott) found Johnston 2109.02249 and **Chirre–Helfgott 2512.15709**. We read the TeX of the latter. Their corollary under T = 2π²√x/log x states Schoenfeld's ψ bound only to 8.2e25. Their Cor. psiexpl (all x ≥ 1, T = 3·10^12) implies it to about 1.66e27, which supersedes everything here.
- **Zenodo** (phrase searches): nothing on this range. Related but different are Daugherty–Ryan–Ward 23475184 (Lean ψ bounds from 7/8, thresholds 2^256 and up), Mittermeier 23544568 (Suzuki screw function) and liqsweep/Padhi (Ramanujan's inequality).
- **Own files:** RamRepair NOVELTY (Büthe method + 7/8 for Ramanujan, no range statement); R295 (missed both Johnston and Chirre–Helfgott).

**Verdict: no record.** What is new but minor:
1. The 7/8 half-plane improves non-band-limited (Logan-kernel) partial-RH bounds by 1.3× in range.
2. Büthe's α > 0 and stepwise zero sum carry over to his √x theorem; certified ranges for ψ, θ, Π, π up to 3.98e26 / 3.77e26 / 3.01e26 / 2.82e26 use refereed inputs only.
3. Johnston's 𝓔₁ slip, and a check that his corollary still holds.

Suggested NOGO line: "7/8 half-plane cannot extend partial-RH prime bounds: the optimal band-limited weights (Chirre–Helfgott 2025) never see zeros above T; it helps only Büthe's Logan kernel (×1.3, still 3× short)."

## 8. RH relevance

Node N10 (explicit ψ and π − li bounds). The ranges are tied to the RH height: X scales like T² up to log factors, and 7/8 does not change that exponent. The useful lever is the RH height T (R286: no extension beyond 3·10^12 found), not the half-plane.
