RECORD: a prime between n^7 and (n+1)^7 for every n ≥ 1 (the smallest exponent known to work for all n drops from 8 to 7); a prime between n^k and (n+1)^k for every n ≥ e^{627.01} (k = 3; before exp(e^{32.537})), n ≥ e^{117.21} (k = 4), n ≥ e^{57.75} (k = 5), n ≥ e^{36.25} (k = 6) (k = 4, 5, 6: first explicit thresholds).

# Ingham3 — primes between consecutive k-th powers from the Lean-checked half-plane and explicit zero density (node N10)

## Claim
**Theorem.** (a) For every integer n ≥ 1 the interval (n^7, (n+1)^7) contains a prime.
(b) For k = 3, 4, 5, 6 and every integer n with log(n^k) ≥ 1881.03, 468.84, 288.73, 217.45 respectively, (n^k, (n+1)^k) contains a prime. So n₀(3) = e^{627.01} ≈ 2.03·10^{272}, n₀(4) = e^{117.21} ≈ 8.0·10^{50}, n₀(5) = e^{57.746} ≈ 1.2·10^{25}, n₀(6) = e^{36.242} ≈ 5.5·10^{15}.

Inputs only: ζ(s) ≠ 0 for Re s > B₂ + 2.2·10^{−254} (Bebel, Lean; B = 0.8733819666132076 used); RH up to 3·10^{12} (Platt–Trudgian 2021); Kadiri–Lumley–Ng's explicit zero density (JMAA 465 (2018), Lemma 4.14, evaluated by us at every σ; its subconvexity input holds by Hiary–Patel–Yang 2024); N(t) ≤ t log t (from Hasanalizade–Shen–Wong 2022). No RH above 3·10^{12}, no unrefereed zero count, and the Chirre–Helfgott low range is not needed.

## Before
- Cubes: n ≥ exp(e^{33.3}) (Dudek 2016), exp(e^{32.892}) (Cully-Hugill, J. Number Theory 247 (2023)), exp(e^{32.537}) (Cully-Hugill's thesis); computation to n³ ≤ 1.649·10^{40} (arXiv 2601.15564, which calls all n "open" and the explicit route "futile").
- All n: k ≥ 90 (Cully-Hugill–Johnston 2402.04272), k ≥ 86 (Lee 2602.14340), k ≥ 8 from the 7/8 half-plane (Visser 2610.20058; our ShortPrimes). Our R381: "k = 7 would need a half-plane below 6/7", true only without density.
- k = 4, 5, 6: only "sufficiently large".

## Method (DERIVATION.md)
1. x = n^k, δ = k log(1+1/n); weight w(u) = Φ(log(u/x)/δ) with Φ the triangle 1_{[0,1/2]} * 1_{[0,1/2]}, supported exactly on [n^k, (n+1)^k]. Explicit formula Σ Λ(n)w(n) = W(1) − Σ_ρ W(ρ) − ∫w/(u³−u), W(s) = δx^sΦ̃(δs): absolutely convergent, no truncation. A prime exists if Z + PP/M₀ + 1/(x³−x) < 1 (Z: weighted zero sum, PP: prime powers in the interval).
2. Kernel: |Φ̃(δρ)|/∫Φ ≤ G(δ|γ|), G = (S(λ/4) + 2δ_max min(1/2, 1/λ))², S the non-increasing envelope of |sinc|.
3. Pairing ρ ↔ 1 − ρ̄: Z ≤ 2[x^{−0.485}·S + log x ∫_{0.515}^{B} x^{σ−1}F(σ)dσ]; S (all zeros, weight x^{−1/2}) via N(t) ≤ t log t; F(σ) (zeros with β > σ, all above T₀) via KLN, F(σ) ≤ A(σ)δ^{−a}D^{5−2σ}Mom(…) + …, a = 8(1−σ)/3.
4. KLN at every σ: for fixed parameters log C₁ is affine in σ; we reimplemented all of KLN's constants (U, V, C₁…C₈, b₁…b₁₂) and reproduce all 20 rows of their Table 1 (A and B within 0.002), then use the 20 parameter sets at all σ in 144 σ-cells.
5. Monotonicity in x: cells [e^{L_a}, e^{L_b}] use x-powers at L_a (log x·x^{σ−1} decreases for log x > 1/(1−σ)) and the kernel at δ(e^{L_b}); beyond L* an explicit majorant R_an(L), proved non-increasing, is checked once.
6. k = 7 below the cells (n ≤ 100): direct prime search.

Why it works: the zeros that bind lie at β ∈ [0.83, B], heights 3·10^{12}–10^{13}; there KLN allows ≤ 2·10^{10}–5·10^{11} zeros, and x^{β−1} ≤ x^{−0.1266}. For k = 7 the worst case is log x ≈ 205.5 (n ≈ 5.6·10^{12}): Z ≤ 0.155. For k = 3 the exponent margin is only (1 − 8/9)(1 − B) = 0.0141 per unit of log x, hence e^{627}.

## Certificate (one command)
uv run --no-project --with python-flint --with mpmath --with sympy python3 Riemann/Deep/Ingham3/cert.py
(Arb, 128 bits; 193 s, 199 MB peak, one thread, one-shot (nothing left running); writes logs/cert.log and logs/cert_result.json.) It checks B > T*, the sinc envelope constant, KLN admissibility and Table 1 reproduction, then every cell and every tail; output "CERTIFIED" per k. k = 7: 3592 cells on log x ∈ [32.236, 240] (max bound 0.1546), tail from 240 (0.244), n ≤ 100 direct. k = 3…6: 1651–2226 cells from log x₀(k), tails 0.0117–0.0876.

## Checks
- Independent float quadrature (check.py; exact KLN shape, best row per σ, true envelope, no moment tables): certificate / direct = 1.011–1.227 at 10 points; logs/check.log.
- δ = 1 in every KLN row (Theorem 1.1 is printed with δ ≥ 1, Lemma 4.14 and Table 1 use δ ≈ 0.303): (a) still certified (max 0.348); log x₀ = 1947.12, 491.86, 305.60, 231.99 (n₀(3) = e^{649.04}).
- Sanity: Ford's zero-free region instead of B in the same pipeline gives n₀(3) = exp(e^{32.90}) (Cully-Hugill: exp(e^{32.892})); the double exponential comes back.
- σ ≈ 0.85 band at T ≈ 10^{13}: KLN gives N(0.85, 10^{13}) ≤ 1.06·10^{11} against the trivial 1.54·10^{13}; crossover only at σ ≈ 0.78; the certificate never uses the trivial count; σ < 0.83 contributes < 7·10^{−4} at the k = 7 worst point.
- Kernel bound tested on 6·10⁴ random points (max ratio 0.9999999995 < 1). Ramaré's b₁, b₂ (inside KLN) rechecked to 10⁷.
- Truncation: none. ψ versus θ: prime powers inside the interval only. Trivial zeros: ≤ W(1)/(x(x²−1)).

## Novelty
NOVELTY.md: arXiv listings Oct–Nov 2026, arXiv API, web search, Zenodo, HAL, GXL evidence API, our files. Nothing on k = 7 for all n or on a single-exponential cube threshold. One unaccepted HAL claim (Diouf 2020, elementary "primes between all consecutive cubes") is not audited; 2601.15564 (2026) treats the problem as open.

## Caveat
- Rests on Bebel's half-plane (Lean-kernel-checked, unrefereed) and on KLN (refereed). KLN's Table 1 used Hiary's 0.63, whose proof had an error; Hiary–Patel–Yang's refereed 0.618 restores it. KLN cite Ramaré's lemmas from a paper whose final density bound had a mistake (KLN footnote); b₃, b₄ were not rechecked by us.
- KLN Theorem 1.1 says δ ≥ 1; we follow their Lemma 4.14 and Table 1 (δ > 0). With δ = 1 everywhere the result (a) stands and log x₀(k) grows by 3.5% (k = 3) to 6.7% (k = 6).
- The explicit formula is the classical one (Davenport ch. 17); no Lean formalization here.
- For k = 3…6 the gap between computation/Chirre–Helfgott (n ≲ kT₀/2π, or n³ ≤ 1.649·10^{40}) and n₀(k) stays open; k = 6 for all n would need about 6× fewer zeros near Re s = B at heights 3·10^{12}–10^{13} than KLN allows (our bound is 5.9 at n = 6T₀/2π).

RH relevance: N10

## Files
- REPORT.md (this), DERIVATION.md (all steps and constants), NOVELTY.md
- cert.py (the certificate), check.py (independent float cross-check, breakdown, kernel test)
- logs/cert.log, logs/cert_result.json, logs/check.log
