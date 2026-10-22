# B♭ = 0.874886647: the stage map, what binds, and attempts to lower it (2026-10-22)

Brief: .scratch/orca/bflat-brief.md. Verdict: [REPORT.md](REPORT.md), NO DISCOVERY.

## 0. Records scan (22 Oct 2026, 03:10 UTC)
- **arXiv.** Listings: math.NT 2026-10 (625 entries with cross-lists), math.NT new (54), math.CV and math.GM 2026-10. The 2026-11 listing is not open yet. API searches of abstracts: "quasi-Riemann" (12), "zero-free half-plane" (7). Quasi-RH constants found: Liu 2610.12234 (0.874957019), Xu–Sang 2610.15291 (0.8748867), Liu 2610.26835 v1 (B♭). Nothing lower.
- **Zenodo.** Searches: "quasi-Riemann" (32 records), "zero-free half-plane" (31), "cubic theta" (4). Constants found, all at least B♭:
  - Chakhkiev 0.87496 (23221313); Mori 0.874957020 (23367242);
  - Newton 0.87488664962, with a barrier of 0.874886649617 inside his fixed κ-range (23414581). Liu cites it;
  - Daugherty–Ryan–Ward, hypothetical 0.87494137 (23462320).
- **Hexagon** (hexagonmath.org API, all 342 works): only 7/8 (Hu 2610.00216), 31/32 (Cureton–Thorner 2610.00200), the CM curves paper (2610.00220) and an implications paper (2610.00178).
- **openai/math:** 749 preprint folders, last commit 12 Oct. Quasi-RH preprints of 30 Sep (7/8) and 5 Oct only.

## 1. The stage map (Liu 2610.26835, TeX l. 8832–10003)
**Setup.** Balanced lengths l_x = (1−b−ℓ)/2, l_y = (1+b−ℓ)/2, h = (1+b+3ℓ)/2; C0 = −1/4 + b/6 + 5ℓ/4. The direct bound matched at B(ℓ) = 11/12 − ℓ/4. A row class has norm U = Z^d, zero class δ = 2a−1 ∈ [1/8, 5/6], and amplitude x = g/δ ∈ [0, 1/2].

**Step.** Input: a proved bound β_* ≤ B. Put κ = 2B − 1, which is what M2 needs: β_* ≤ (1+κ)/2. The output is B′ = B(ℓ) for every (b, ℓ) with (b, ℓ, κ) admissible:
- (A1) b > 0, 1/6 ≤ ℓ ≤ 1/5, 0.7497 ≤ κ ≤ 3/4, B(ℓ) > 437/500, l_x > ℓ, l_y − ℓ > 11b/6.
- (A2) 6κℓ > 3/4, 21/100 < ℓ/h < 1/(6κ), h/2 < l_y < 1/2.
- (A3) E(δ, x) < 0 on [1/8, 5/6] × [0, 1/2], where
  E = C0 + δ/2 + δx{ℓ + (h/4 − l_y/2)/(3κ)} − (3h/4)(1 − R_κ(δ, x)).
- (A4) max{C0, C0 + A_sm/8} < −1/250, C0 + (1/2 + 3ℓ/2)/50 < −1/100, 13h/75 − l_y/2 + 1/50 < −1/20, with A_sm = 1/2 + (ℓ/2)(1/2 + 2ℓ/h) − h/2.

**Where E comes from.** Each term is one input. Write g = δx.
- **Count.** R_κ = min over t ∈ [1, 3/2] of max(R_short(t), R_long(t)), where R_short(t) = max over r of min{A_I(r), A_S(t−r)}.
  - A_I(r) = 1 − δr − g(1−r): M1, the Möbius second moment, which allows primes up to r + 2z ≤ 1.
  - A_S(m) = 1 − 2δm − 2g(1−2m)/(6κ): M2, the plain fourth moment, which allows primes up to 2m + 6κz ≤ 1.
  - R_long(t) = 1/6 + (5/6 − δ)t: sixth-power amplification, Σ|M_u(U^r)|² ≪ U^{(1+5r)/6}.
  - Result: R_κ = 1 − δ + (5/6−δ)δP/(2J), with D = 3 − (1+2c)x, P = (2−2cx)(1−x), J = (5/6−δ)D + δP and c = 1/(3κ). This is Liu's (eq:v-refined-count); bflat_exact.py derives it symbolically from the three inputs.
- **Numerator.** Hölder plus M2 give U^{(3R+1)/4 − θ/2}, with shift θ = g/(6κ).
- **Contour.** Factors Z^{a−B}, Y^{−1/2−θ}, prime sums Z^{gℓ}, and the reflection cost U^θ.

**Closed form** (Liu, Lemmas v:lem:refined-critical and v:lem:refined-family).
- At x = 1/2, ∂E/∂b is proportional to cr(δ,κ) = (3−324κ²)δ² + (756κ²−54κ−5)δ + 40κ − 300κ².
- Its root δ_c(κ) ∈ (23/50, 12/25) is the critical class. There E = λ(ℓ − ℓ*(κ)) does not depend on b, and b*(κ) solves ∂_δE = 0.
- Hence the step yields every B′ > Φ(B) := F♭(2B−1), with F♭(κ) = 11/12 − ℓ*(δ_c(κ), κ)/4, and nothing at or below it.

**Fixed point.** F♭(κ) = (1+κ)/2 is the same as ℓ* = 5/3 − 2κ.
- The resultant in δ is −10800κ²(6κ−1)²(15κ−2)(42κ+5)·(16686κ⁴ − 27405κ³ + 9936κ² + 1290κ − 275). In B this is 2·(133488B⁴ − 376596B³ + 384534B² − 167541B + 26231), Liu's quartic, reproduced exactly.
- The quartic changes sign on κ ∈ [7497732940799, 7497732940800]·10⁻¹³ (exact rational check).
- B♭ = 0.87488664703996413094…, κ♭ = 0.74977329407992826189….

**Iterates from 7/8.** 0.8748869894, 0.8748866481, 0.874886647043, 0.8748866470399736. Errors: 3.4e-7, 1.0e-9, 3.1e-12, 9.4e-15. Φ′(B♭) = 0.00302.

**κ is not what stops the iteration.** F♭ is nearly flat:
- F♭(1) = 0.875169: with no zero-free input the system is worse than 7/8.
- F♭(3/4) = 0.874887.
- F♭(1/2) = 0.874312 (formal).
The iteration stops where this almost constant curve crosses (1+κ)/2.

## 2. What binds at B♭ (bflat_exact.py, bflat_checks.py)
**Tight constraints.** (A3) is tight at exactly one class, (δ_c, x) = (0.4693448587, 1/2), as a double root in δ. Also tight: κ = 2β_* − 1 and B′ = B(ℓ). The point is ℓ* = 0.1671200785, b* = 0.0865349593, h = 0.7939476, l_y = 0.4597074.

**At the critical class:**
- **Count.** R = 0.5903323535; the density hypothesis would give 0.5306551. Three estimates are tight at once:
  - M1 at r* = 0.7457, with primes z_M = 0.1272;
  - M2 at m* = 0.4183, with primes z_S = 0.0363;
  - amplification at t* = 1.1640.
- **Numerator.** K = (3R+1)/4 = 0.6927493, θ = 0.0521652.
- **Amplitude cap.** ∂E/∂x = +0.0324, so the cap x ≤ 1/2 is active. The cap says prime sums are at most P^{δ/2}, which follows from the zero-free rectangle.
- ∂E/∂ℓ = 1.089.

**Sensitivity** (bflat_sens.py; same pipeline, one constant changed):

| constant (source) | published | dB♭/d(constant) | example |
|---|---|---|---|
| amplification exponent (OpenAI Lemma 17.6) | 5/6 | +0.00737 | 4/5 → 0.8746276; 2/3 → 0.8732353 |
| M1 prime cost (Lemma 17.1: r + 2z ≤ 1) | 2 | +0.00208 | 7/4 → 0.8742420; 3/2 → 0.8731969 |
| M2 prime cost (Lemma 18.1: n1 + n2 + 6κz ≤ M) | 6 | +0.000189 | 5 → 0.8746582; 3 → 0.8736975 |
| numerator moment order (Hölder) | 4 | −0.00349 | sixth moment: 0.869351 (QRHImprove) |
| amplitude cap x ≤ | 1/2 | +0.00747 | a theorem, not a choice |
| licensing κ = 2B − 1 | — | Φ′ = 0.003 | κ = 1/2: 0.874312 (formal) |

**Slack at B♭.** Exact rational slacks were computed at a rational point within 10⁻²⁰ of the tangent point. All are positive except OpenAI's 2/27 bound, which applies only for κ ≥ 3/4 and is not used by the refinements.

| constraint | value | bound | slack |
|---|---|---|---|
| row count c/12: selected plain length z_S ≤ 1/(18κ) | 0.03634 | 0.07410 | 0.0378 |
| Xu–Sang: 1/(18κ) ≤ 1/13 | 0.07410 | 0.07692 | 0.0028 |
| OpenAI plain capacity 2/27 (κ ≥ 3/4 only; unused here) | 0.07410 | 0.07407 | −2.2e-5 |
| prime supply ℓ/h against the worst class need (δ = 1/8) | 0.17734 | 0.21049 | 0.033 |
| 13/15 floor: ℓ ≤ 1/5, κ ≥ 11/15 | 0.16712, 0.74977 | 0.2, 0.73333 | 0.033, 0.016 |
| Lemmas 18.2/18.3 | needed (M2 in its centred range: 2m* + z_S = 0.8728 > 5/6) | no constant | — |
| M1 side condition 2r + 8z ≤ 3 | 2.5086 | 3 | 0.49 |
| detection ranges for t, r, m | | | ≥ 0.08 |
| numerator gate R ≤ (6κ − 2g)ℓ/h | 0.5903 | 0.8481 | 0.26 |
| (A1) ℓ ≥ 1/6; (A2) ℓ/h > 21/100; 6κℓ > 3/4 | | | 4.5e-4; 4.9e-4; 1.8e-3 |
| (A4a) small-δ range (coarse bound R ≤ 1 − 2δ/3) | −0.004179 | −0.004 | 1.8e-4; true margin with the exact count: 0.0157 |

## 3. Attempts to beat B♭ (all fail; numbers from bflat_checks.py)
1. **More stages.** Φ is a 0.003-contraction toward B♭, so no finite chain from 7/8 goes below it (Liu, Thm v:thm:refined-optimality; our iterates agree).
2. **Averaged or pointwise numerator, chosen class by class** (outside Liu's theorem). At ℓ = ℓ* + ε the best b still leaves max over classes = +1.089ε. The pointwise endpoint at δ_c is +3.96e-4. No gain.
3. **Recombining M1, M2 and amplification.**
   - Hölder between the two short bounds: their convex mixtures never beat the minimum.
   - Several cutoffs t at once: the worst witness pair can be taken separately for each t, so the bound is the minimum over t, which is already used.
   - M2 with a trivial second factor (|S_m·1·Q|²): worse by mδ(1 − x/(3κ)) > 0.
   - Prime sums alone (M1 at r = 0): U^{1−g} = U^{0.765}.
   - Splitting the long witness M_t by Heath-Brown's identity into M_a·M_a·S_b (a ≤ t/2, b ≤ 1/2): at best 0.6657 > 0.5903.
   - Amplifying with selected primes in the enlarged family (H = U^{t+2z}): the exponent grows by (5/3 − 2g)z.
   - The sextic large sieve with general coefficients (U + N + (UN)^{2/3}): mean-square exponent 1.443 against 1.137 for amplification at t*, and worse than M1 for r ≤ 1.
   - Halász–Montgomery with prime polynomials: needs length ≥ U^{R/(1−β_*)} = U^{4.7}.
4. **Constants inside the lemmas.**
   - **5/6.** The multipliers a with q_a ≤ (H/U)^{1/6} are all the sixth powers there are, and H = max(U, D^{1+c}) is the best choice. Other multipliers (u ↦ uc², uc³) change the character.
   - **6κ.** It is the pointwise slot bound |Q|² ≤ Z^{κz} (OpenAI (3.5)) set against the budget M/6. That budget is left by the uncentred threshold 5M/6, which comes from the sixth-power form of the Θ-rows ((3.6); greedy children κd_z ≤ F_act/6). Daugherty–Ryan–Ward (Zenodo 23462320) certify the 6 as optimal for this accounting. Even κ → 1/2 gives only 0.874312.
   - **The M1 cost 2.** Slots enter OpenAI's §17 recursion as arbitrary bounded mark coefficients, through the invariants F−M−Q−z0 ≥ c and 4F−3M−6z0 ≥ c (l. 9296–9338, 9808, 12208). No step bounds a prime sum pointwise, so β_* cannot enter. A lower cost needs a new M1.
   - **Lemmas 18.2/18.3.** They are κ-free Poisson identities with no exponent to tune.
5. **Unequal lengths, smaller θ, other detectors, and the 13/15 floor.** No gain; see QRHImprove REPORT.md §2 and REPORT3.md. The floor binds only for ℓ > 1/5.

## 4. What a record needs
Any one of these would lower B♭; each is a new theorem for the thin sextic family:
- a count below U^{0.590332} for the class (δ, x) = (0.4693, 1/2), that is, a zero-density exponent below 2.2249 at σ = 0.7347 (QRHImprove REPORT2: the best known for ζ there is 2.2477);
- an amplification exponent below 5/6, an M1 prime cost below 2, or an M2 prime cost below 6κ;
- a numerator average below U^{(3R+1)/4 − θ/2} on the class, for example from a sixth moment.

None is in the literature (REPORT2 checked de Faveri–Pandey, Guth–Maynard, Chen–Gupta–Li and Heath-Brown). After any such gain, the next constraint to tighten is Liu's coarse small-δ range (A4a). It is only a convenience bound: with the exact count its margin is 0.0157.

## Sources, scripts, cost
- **Sources.** Liu arXiv 2610.26835 (number_fields_v2.tex: §v:sec:refined, l. 1360–1560 and 8832–10003); Liu 2610.12234 (combine_v2-2.tex: §§19–20, l. 1416–2437, and §29); Xu–Sang 2610.15291; OpenAI 7/8 paper.tex (§17 l. 9209–12476, §18 l. 12477–14984); Zenodo 23462320 and 23414581.
- **Scripts:**
  - bflat_exact.py: symbolic count, critical class, resultant, fixed point. 2 s, 56 MB.
  - bflat_sens.py: stage map and sensitivities. 4 s CPU, 96 MB.
  - bflat_checks.py: slacks (exact rationals) and the attempts. 3 s, 41 MB.
- **Cost.** No Lean builds, no commits. All downloads were in .scratch/bflat/, which is now deleted.
