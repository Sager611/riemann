RECORD: Robin for 5041 <= n <= 10^(10^46.79), Nicolas for p_k <= 1.78*10^48 (b0 = 0.87338196661320752)

# Robin and Nicolas ranges with the Chirre–Helfgott pointwise bound (node N4)

Agent robinch, 2026-11-05. Base: Riemann/Deep/RobinHeight (agent robinh). Files in this folder:
- PROOF.md: the new lemmas (12: pointwise θ from Chirre–Helfgott; 11_k: k-box ramp; 13: sharper second-order term).
- certify3.py: the certificates (levels 2, 3, 4) and the endpoint search; certify.py, certify2.py: unchanged copies used as libraries.
- check_all2.py: reruns every inequality behind this report (ALL PASS, 460 s on one core, 80 MB).
- test_b2k.py: numerical test of Lemma 11_k (448 cases, PASS).

## 1. Short answer

1. **The brief's change works as predicted.** Replacing our own pointwise bound η ≈ 4·10^−11 (Lemma 10_k) by the Chirre–Helfgott bound η ≈ π/(3·10^12) = 1.05·10^−12, and changing nothing else, moves the largest prime by a factor of 71 (Robin) and 55 (Nicolas):
   - Robin for 5041 ≤ n ≤ 10^(10^40.56) (largest prime 8.47·10^40; before 1.19·10^39), Nicolas for p_k ≤ 7.17·10^41 (before 1.31·10^40), b0 = 0.87338196661320752.
2. **With a cheap pointwise bound, a smoother cutoff pays.** The tail integral's cutoff at x was smoothed by a linear ramp. A ramp built from k boxes (the distribution function of a sum of k uniform variables), centred at x, damps the zeros above 3·10^12 by a factor (2k/(λγ))^k, while its window still costs only η times the ramp's mean absolute deviation. With a sharper bound for Robin's second-order term this gives:
   - **Robin for 5041 ≤ n ≤ 10^(10^46.79)** (largest prime 1.43·10^47), **Nicolas for p_k ≤ 1.78·10^48** (b0 = 0.87338196661320752). That is about 10^8 times RobinHeight's largest prime (1.19·10^39).
3. **The half-plane now matters little.** With no zero-free half-plane at all, the same method gives Robin to 10^(10^46.44) and Nicolas to 6.30·10^47. Chirre–Helfgott needs only RH up to 3·10^12 + 1 + π/3 (Platt–Trudgian verified 3,000,175,332,800), and the k-box ramp damps the zeros above 3·10^12 whatever their real part.
4. **Our own pointwise lemmas are gone.** Lemmas 10 and 10_k are no longer used; the pointwise bound is Chirre–Helfgott (a preprint) plus Büthe below 10^19. New own lemmas: 11_k and 13 (short, reviewed, numerically tested) and 12 (two lines).
5. **What limits X now.** Pointwise knowledge of θ(x) at relative precision π/T0. The window costs about √x · (π/T0) · (λ/k)·MAD_k with λ T0 ≈ 90 and k ≈ 30–60, and Robin's second-order term costs √x η²/2; both reach the margin near x ≈ 10^47–10^48, i.e. roughly x ≈ (T0²/c)² with c ≈ 13 (a heuristic, not a theorem). Possible routes to a further gain: RH to a greater height, or a Beurling–Selberg-type treatment of the corner of the tail kernel itself.

## 2. Ranges

All rows of this work are certified by check_all2.py in arb ball arithmetic (192 bits), from 10^7 (Nicolas; below: direct check of all 664,579 primorials) and from 10^10 (Robin; below: Morrill–Platt).

| method | zero input above T0 = 3·10^12 | Nicolas p_k ≤ | Robin: largest prime ≤ | Robin: all n ≤ |
|---|---|---|---|---|
| **Level 3 (centred k-box ramp + CH)** | Re ρ ≤ 0.87338… | **1.78·10^48** | **1.43·10^47** | **10^(10^46.79)** |
| Level 3 | Re ρ ≤ 7/8 | 1.76·10^48 | 1.42·10^47 | 10^(10^46.79) |
| Level 3 | none (β < 1) | 6.30·10^47 | 6.46·10^46 | 10^(10^46.44) |
| Level 2-CH (brief: only η changed) | Re ρ ≤ 0.87338… | 7.17·10^41 | 8.47·10^40 | 10^(10^40.56) |
| Level 2-CH | Re ρ ≤ 7/8 | 5.99·10^41 | 7.11·10^40 | 10^(10^40.48) |
| Level 2-CH | none | 3.57·10^36 | 5.51·10^35 | 10^(10^35.37) |
| before: RobinHeight Level 2 (Lemma 10_k) | Re ρ ≤ 0.87338… | 1.31·10^40 | 1.19·10^39 | 10^(10^38.71) |
| before: RobinHeight Level 2 | none | 1.42·10^34 | 2.18·10^33 | 10^(10^32.98) |
| Polak 2026 (Zenodo, sharp cutoff) | none | – | 1.65·10^23 | 10^(10^22.85) |
| Morrill–Platt 2018 (refereed) | – | – | 3.0·10^13 | 10^(10^13.11) |

The "all n ≤" column is 10^(10^Y) with Y = log10(θ(X_R)/log 10), θ(X_R) ≥ X_R(1 − η) (since log N_ε ≥ θ(X_R)), Y rounded down to two decimals (check_all2.py prints Y to four).
Intermediate step (search log only, not in check_all2.py): the one-sided k-box ramp (λ/2 window cost) reaches Nicolas 2.44·10^46 and Robin 3.08·10^45 (b0 = 0.87338…).

Tradeoffs:

| | Level 2-CH | Level 3 |
|---|---|---|
| window cost (normalised) | √x η λ/2 | √x η (λ/k) MAD_k ≈ √x η λ·0.23/√k |
| cost of a zero above T0 | x^(β−1)(2/λ)/γ³ | x^(β−1)(2k/(λγ))^k/γ² (+ tiny terms) |
| λ T0 at the end of the range | ≈ 7000 | ≈ 90 (k ≈ 30–60) |
| Robin second-order term | η²/log x (RobinHeight) | η²/(2 log x)(1 + …) (Lemma 13) |
| what limits X | window vs zeros: X ~ (c/(S3 η))^(1/b0) | window + D: X ~ (T0²/c)², almost independent of b0 |
| new own lemmas | none (Lemma 12 is two lines) | 11_k, 13 |
| status | certified; changes only η | certified; reviewer: formulas and code OK, two wording errors fixed |

## 3. The change

### 3.1 Chirre–Helfgott (H10)

- Paper: A. Chirre, H. A. Helfgott, "Optimal bounds for sums of non-negative arithmetic functions", arXiv 2512.15709v1 (17 Dec 2025, 52 pp.). Status: preprint; no journal reference on arXiv (checked 5 Nov 2026).
- Statement used, Corollary cor:psiexpl, TeX line 404 (with l.403 "For any $x\geq 1$,"):
  "For any $x\geq 1$, \[|\psi(x)-x|\leq \frac{\pi}{3\cdot 10^{12}}\cdot x + 113.67 \sqrt{x},\]"
- It is ψ, not θ. Its proof (TeX l.4981–4991) applies Corollary cor:psi = Proposition prop:sagaro (σ = 0, TeX l.4868–4880) with T = 3·10^12 + 1 + π/3 for x ≥ 3·10^12 + 3, and a brute-force lemma below.
- Hypothesis of Prop. prop:sagaro: "the Riemann hypothesis holds up to height T ≥ 10^7", which the paper defines (l.387–388) as all non-trivial zeros with 0 < Im s ≤ T having real part 1/2; range x > max(T, 10^9). No zero-free region and no zero-density input.
- Our H1 (Platt–Trudgian: RH to 3,000,175,332,800) implies it. We use H10 only for x ≥ 10^19·e^(−10^−3), where it rests on Prop. prop:sagaro.
- check_all2.py (0) rechecks the two numerical steps of the corollary's proof in arb: C_T = (1/2π) log²(T/2π) − (1/6π) log(T/2π) = 113.66888… < 113.67, and π²/(3T²) + π/(T − 1) ≤ π/(3·10^12).

### 3.2 Where η enters (Lemma 12)

- η_CH(x) = π/(3·10^12) + 113.67/√x + P(x)/x, with P(x) = √x + x^(1/3) + x^(1/4) + (log x/log 2) x^(1/5) (times 1 + η_1 above 10^38) our existing bound for ψ − θ. Then |θ(x) − x| ≤ η_CH(x)·x in both directions: θ − x ≤ ψ − x, and x − θ = (x − ψ) + (ψ − θ).
- η := min(Büthe 1.95/√x when the box and its window lie below 10^19, η_CH). It is used exactly where Lemma 10_k was: the window term, Robin's second-order term D(x), and θ(x) ≤ x(1 + η) in Robin's L_2.
- η_1 (θ on [10^19, 10^38]: prime-power zone 10^38–10^76, Robin's prime-square sum, θ(x_2) ≥ x_2(1 − η_1)) is now η_CH(10^19) = 3.6263·10^−8 instead of Lemma 10's 2.70·10^−8. This removes Lemma 10 completely; rerunning Level 2-CH (b0 = 0.87338…) with Lemma 10's η_1 gives the same endpoints to six decimals (10^41.855569 and 10^40.928055).
- At the end of the range η = 1.04720·10^−12 (η_CH(10^40)), against 3.6·10^−11 from Lemma 10_k at 10^40.

### 3.3 Level 3 (Lemma 11_k, Lemma 13)

- Ramp: φ(t) = Φ(log t − log x + λ/2), Φ the distribution function of V = U_1 + … + U_k, U_i uniform on [0, λ/k]. Mean λ/2, so the ramp is centred at x.
- Window: |W| ≤ η H(L′) E|V − λ/2|, E|V − λ/2| = (λ/k) MAD_k, MAD_k = (2/(k+1)!) Σ_{j≤k/2} (−1)^j C(k, j)(k/2 − j)^(k+1) exactly (≈ 0.8 of the standard deviation λ/√(12k)).
- Zeros above T0: |F_φ(ρ)| ≤ x′^(β−1)[H(L′)(2k/(λ|s|))^k/|s| + |H′(L′)|(λ/(2|s|) + 2/|s|²)], s = ρ − 1, because the Fourier–Laplace transform of V is ((e^(sa) − 1)/(sa))^k with |e^(sa) − 1| ≤ 2.
- Zeros below T0, the prime-power integral M, the log 2π and trivial-zero terms: as in Lemma 11, with x replaced by x′ = x e^(−λ/2) or x e^(λ/2) on the safe side.
- Lemma 13: D(x) ≤ (η²/(2L))[1 + 2η/(3(1 − η)) + (1 + η)²/L + …], asymptotically half of RobinHeight's η²/L(1 + (1 + η)²/L).
- Per box: k ∈ {1, …, 120} and λ minimise window + zeros (float), then the margin is certified in arb for the best two k and three λ.

## 4. Certificates

- Command: uv run --no-project --with python-flint --with mpmath --with numpy python3 check_all2.py (in this folder).
- Output (run 5 Nov 2026, 19:51–19:59 CET):

    (0) PASS: C_T = [113.6688806 +/- 4.72e-8] < 113.67: True; pi^2/(3T^2) + pi/(T-1) <= pi/3e12: True; T = 3e12 + 1 + pi/3 <= 3000175332800 (Platt-Trudgian): True
    (1) PASS: Nicolas holds for all 664579 primorials with p_k <= 9999991 (min of margin*sqrt(p)log p = 1.9517 at p = 5); 7s
        eta_1 = eta_CH(1e19) = [3.62631e-8 +/- 4.61e-14]; eta_CH(1e40) = [1.04720e-12 +/- 1.31e-18]
    (2) Level 2-CH, b0 = 0.87338196661320752: Nicolas on [1e7, 7.17e41]: PASS (373 boxes, min 5.04e-06, k <= 0); Robin on CA with largest prime in [1e10, 8.47e40]: PASS (1212 boxes, min 8.53e-06, k <= 0) => Robin for 5041 <= n <= 10^(10^40.5656); 77s
    (2) Level 2-CH, b0 = 7/8: Nicolas on [1e7, 5.99e41]: PASS (354 boxes, min 0.000528, k <= 0); Robin on CA with largest prime in [1e10, 7.11e40]: PASS (1177 boxes, min 4.57e-06, k <= 0) => Robin for 5041 <= n <= 10^(10^40.4896); 176s
    (2) Level 2-CH, b0 = 1: Nicolas on [1e7, 3.57e36]: PASS (301 boxes, min 0.000671, k <= 0); Robin on CA with largest prime in [1e10, 5.51e35]: PASS (1125 boxes, min 3.3e-06, k <= 0) => Robin for 5041 <= n <= 10^(10^35.3789); 221s
    (3) Level 3, b0 = 0.87338196661320752: Nicolas on [1e7, 1.78e48]: PASS (404 boxes, min 0.00147, k <= 40); Robin on CA with largest prime in [1e10, 1.43e47]: PASS (1161 boxes, min 1.63e-06, k <= 40) => Robin for 5041 <= n <= 10^(10^46.7931); 298s
    (3) Level 3, b0 = 7/8: Nicolas on [1e7, 1.76e48]: PASS (404 boxes, min 0.00109, k <= 40); Robin on CA with largest prime in [1e10, 1.42e47]: PASS (1189 boxes, min 8.78e-06, k <= 40) => Robin for 5041 <= n <= 10^(10^46.7900); 379s
    (3) Level 3, b0 = 1: Nicolas on [1e7, 6.30e47]: PASS (419 boxes, min 9.35e-06, k <= 57); Robin on CA with largest prime in [1e10, 6.46e46]: PASS (1171 boxes, min 2.03e-06, k <= 57) => Robin for 5041 <= n <= 10^(10^46.4480); 460s
    ALL PASS total 460s

- Runtime 460 s wall (296 s CPU) on one core, peak memory 80 MB (/usr/bin/time -l).
- Endpoint searches (certify3.py 2 / 3 / 4; adaptive boxes from width 0.25 down to 10^−6 in log x): Level 2-CH 10^41.855569 / 10^40.928055 (0.87338…), 10^41.778132 / 10^40.852343 (7/8), 10^36.553270 / 10^35.741733 (none); Level 3 10^48.252447 / 10^47.156206, 10^48.247204 / 10^47.152610, 10^47.799392 / 10^46.810759 (Nicolas / Robin). The claims above are these endpoints rounded down to three digits.
- test_b2k.py (mpmath, 448 cases: β ∈ {0.5, 0.875, 0.97}, γ ∈ {30, 300, 3000}, L ∈ {23, 70}, k ∈ {1, 2, 4, 8, 16, 32}, one-sided and centred, 2k/(λγ) ∈ {0.2, 0.3, 0.6, 2}): worst ratio true/bound 0.9954 for (B2_k) and 0.9725 for (B1′); the closed forms for F and MAD_k agree with direct quadrature. PASS.

## 5. Self-review

Pass 1 (code against formulas, by me): every term of margins3/margins4 matches PROOF.md; Level 2-CH differs from certify2.py only in η, η_1 and the λ grid. Found and fixed: λ was capped at the float 10^−3, a hair above the lemma's 10^−3 (RobinHeight had the same edge); now 9·10^−4.

Pass 2 (adversarial, the brief's list):
1. Chirre–Helfgott's hypothesis is RH up to height 3·10^12 + 1 + π/3, implied by H1 (3,000,175,332,800). It uses no zero-free region, so the b0 = 1 rows need no half-plane (H2) at all: only H1, H3–H10 (RobinHeight's inputs plus Chirre–Helfgott).
2. ψ → θ: the upper bound θ − x ≤ ψ − x and the lower bound x − θ ≤ (x − ψ) + (ψ − θ) are both used; upper bounds of θ enter L_2 and D, lower bounds enter L_2 (θ(x_2)), log N_ε ≥ θ(X_R), and M (θ(√t)).
3. Ranges: Chirre–Helfgott holds for every x ≥ 1, so there is no gap; Büthe is used only when the window stays below 10^19; Lemma 12's ψ − θ bound is valid to 10^76 and we stop at 10^49.
4. Window uniformity: η is evaluated at the left end of each window (ball L − λ/2), the factor 1 + η_1 is switched on for any box or window reaching beyond 10^38 (η_CH jumps up by 3.6·10^−27 there).
5. Robin's margins compare n(x_a)·(lower bound of the positive terms) with the sup over the box of n(x)·(negative terms); valid since n increases.

Independent reviewer (read-only subagent, PROOF.md sections 1–3 and certify3.py): Lemma 11_k (both bounds, Fubini, signs, the transform bound), the centred decomposition and the M endpoints in the code, the pairing and normalisation, Lemma 13's derivation and nDup2, and Lemma 12's inequalities are OK. It flagged two wording errors, both fixed: Lemma 13's "twice as large" holds only as η → 0, and η_CH is decreasing only for a fixed factor (jump at 10^38; the code already used the larger factor).

Second reviewer (REPORT.md and PROOF.md against the check_all2.py output): all ranges round down correctly except one exponent (Level 2-CH, 7/8: 40.49 → 40.48); the suggested DISCOVERIES line rounded up and left out the half-plane input; Chirre–Helfgott's height is 3·10^12 + 1 + π/3, not 3·10^12; C_T was mistyped in PROOF.md; novelty statements must read "found in these searches". All fixed; it confirmed the two PROOF.md repairs.

## 6. Novelty pass (5 Nov 2026)

- arXiv math.NT listings 2026-10 (866 entries) and 2026-11 (141), keyword scan of titles (scout): only Molnar 2610.12816 (a parameterised Robin-type inequality, no range in n).
- arXiv search (arXiv mode): nothing beyond 2510.23889 (structure of a least CA counterexample) and older claimed proofs.
- Zenodo phrase searches ("Robin's inequality", "Robin inequality", "colossally abundant", "Nicolas inequality", Chirre AND Helfgott, primorial AND Riemann): Polak 21808589 (10^(10^22.85)) is still the only deposit with a verified Robin range found in these searches; Fabbian 23782485 (uniform Mertens bound); Daugherty–Ryan–Ward 23501709 (Lean-formalised Chirre–Helfgott argument for Dusart's bound on Σ 1/p, x ≥ 10372) and 23475184 (Lean Chebyshev bounds from the 7/8 half-plane): neither gives a Robin or Nicolas range. Others are claimed proofs or unrelated.
- Our files: RobinHeight (10^(10^38.71)); R334/R344 records. No work found in these searches combines Chirre–Helfgott with Robin or Nicolas.
- Verdict: new record (Level 2-CH: routine substitution of a preprint bound into our method; Level 3: new smoothing step for this problem).

Suggested lines for the orchestrator:
- TREE.md N4 record: "Robin verified for 5041 ≤ n ≤ 10^(10^46.79), Nicolas for p_k ≤ 1.78·10^48 (RH to 3·10^12 + 1 + π/3 via Chirre–Helfgott's pointwise ψ bound, the 0.87338… half-plane, centred k-box smoothing; robinch); without any half-plane 10^(10^46.44) and 6.30·10^47; before: RobinHeight 10^(10^38.71), Polak 10^(10^22.85)".
- DISCOVERIES.md: "Robin's inequality verified to 10^(10^46.79) and Nicolas' to p_k ≤ 1.78·10^48 (Chirre–Helfgott ψ bound, k-box smoothed tail, 0.87338 half-plane); 10^(10^46.44) without the half-plane (N4)."

## 7. Caveats

- Chirre–Helfgott is a preprint (v1, Dec 2025). Every row now depends on it above 10^19 (Level 2-CH and Level 3 alike).
- Lemmas 11_k, 12, 13 are ours and unpublished: one mathematical reviewer pass (formulas and code OK), one consistency pass, and the numerical tests of Lemma 11_k.
- The certificates trust python-flint (arb); E_1 is computed as −Ei(−z), as in RobinHeight.
- The half-plane inputs (7/8, 0.87338…) are no longer what drives the range; the b0 = 1 rows need no zero-free half-plane: RH to 3·10^12 + 1 + π/3 (with Chirre–Helfgott), the N(T) bound, Büthe, Morrill–Platt and classical facts.

RH relevance: N4 (Robin's criterion: verified range of an RH-equivalent inequality)
