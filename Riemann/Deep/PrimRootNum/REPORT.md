RECORD: from the Lean-checked zero-free half-plane (B = 0.8733819666132076), every prime p ≥ 3 has a prime primitive root g*(p) ≤ (1.956 log p)^{3.4/(1−B)} ≤ (1.956 log p)^{26.86}; if p ≥ 5 and ω(p−1) ≤ 10⁴ (every p ≥ 5 below e^{104403}), g*(p) ≤ (0.814 log p)^{3/(1−B)} ≤ (0.814 log p)^{23.70}; and g*(p) ≤ (κ_ω log p)^{7.898} with κ_ω tabulated for ω = ω(p−1) ≤ 10⁴.

## Claim
Assume H: ζ(s) and every Dirichlet L(s, χ) are zero-free on Re s > B. H is a theorem (Bebel, QuasiRiemannTracker, Lean-kernel-checked:
Re s > B₂ + 2.2·10⁻²⁵⁴, B₂ = (777 − √1081)/852 < B; TREE N7). Let g*(p) be the least prime that is a primitive root mod p and ω = ω(p − 1).

1. **Theorem A (per ω).** For every ω ≤ 10⁴ and every prime p with ω(p − 1) = ω: g*(p) ≤ (κ_ω log p)^{1/(1−B)}, 1/(1−B) = 7.8978, with
   κ₁ = 3.53 (p − 1 a power of 2), κ₂ = 11.1 (e.g. Sophie Germain p = 2l + 1), κ₅ = 54.1, κ₁₀ = 216, κ₂₀ = 868, κ₅₀ = 6130,
   κ₁₀₀ = 2.72·10⁴, κ₁₀₀₀ = 4.32·10⁶, κ₁₀⁴ = 6.85·10⁸ (all 10⁴ values in logs/kappa_table.txt; κ_ω ≤ 0.095 p_ω² for
   ω ≥ 100 and ≤ 0.070 p_ω² for ω ≥ 1000, p_ω the ω-th prime: Martin's shape (ω log ω)² with a number).
2. **Theorem B (exponent 3/(1−B) = 23.693).** For every prime p with ω(p − 1) ≤ 10⁴ — in particular every p < e^{104403} —
   g*(p) ≤ (1.4312 log p)^{3/(1−B)}, and g*(p) ≤ (0.814 log p)^{3/(1−B)} if moreover p ≥ 5.
3. **Theorem C (the rest).** If ω(p − 1) > 10⁴ (then p > e^{104403}): g*(p) ≤ (1.955 log p)^{3.4/(1−B)} ≤ (1.955 log p)^{26.86}.
4. **Theorem D (every prime).** For every prime p ≥ 3: g*(p) ≤ (1.956 log p)^{26.86}.
5. **Corollary E.** g*(p) < √p − 2 for every prime p ≥ 10^{97.547}; g*(p) < p for every prime p ≥ 10^{41.032} (cross-check, below).

Table (log₁₀ of the bound for the worst admissible ω and for typical ω):

| p ≈ | max ω(p−1) | worst ω (Thm A) | Thm B (0.814) | Thm D (1.956) | ω = 2 | ω = 5 |
|---|---|---|---|---|---|---|
| 10¹⁰ | 10 | 10^29.2 | 10^30.2 | 10^44.4 | 10^19.0 | 10^24.4 |
| 10¹⁰⁰ | 53 | 10^49.0 | 10^53.9 | 10^71.3 | 10^26.9 | 10^32.3 |
| 10¹⁰⁰⁰ | 350 | 10^71.0 | 10^77.5 | 10^98.1 | 10^34.8 | 10^40.2 |
| 10¹⁰⁰⁰⁰ | 2584 | 10^94.1 | 10^101.2 | 10^125.0 | 10^42.7 | 10^48.1 |

## Before
- Unconditionally, no bound below p was known for g*(p): "It is not known whether ĝ(p) < p for all p, or even for all sufficiently
  large p" (Cohen–Trudgian, arXiv 1602.02440). Our own Linnik theorem (Riemann/Deep/ExplProg) gives only g*(p) ≤ (p log p)^{7.898}.
  For the least primitive root g(p) (not prime): g(p) < p^{5/8} for p > 3 (arXiv 2206.14193), g(p) < √p − 2 for p > 1.77·10⁵⁴ (2608.20431).
- Under GRH: ĝ(p) ≤ ((8/5)(2^{ω(p−1)} − 1) log p)² for p ≥ 10⁹ and ĝ(p) < √p − 2 for p > 2791 (McGown–Treviño–Trudgian, 1508.05182);
  Shoup/Martin g*(p) ≪ (ω log ω)⁴ log² p (no constants).
- From a zero-free rectangle: Martin, Acta Arith. 80 (1997), Prop. 3, g*(q) ≪_σ (ω² log ω c₀⁻¹ log q)^{1/(1−σ)} (no constants).
  October 2026 Zenodo notes (dr vengy: 23341593, 23352868, 23394711, 23423890, 23621868) apply it to OpenAI's 7/8: g*(p) ≪ (log p)²⁴,
  "No numerical value of C is asserted". Zenodo 23623384: g*(p) < p for p ≥ 10⁴¹ (later 1.4·10³³), QRH-conditional, unrefereed,
  plain sieve at x = p; no polylogarithmic numeric bound.

## Method
DERIVATION.md has every step. In short:
- Detection with Λ(m)w(m/x) over prime powers m (r^k generates G ⇒ r does), ExplProg's kernel and lemmas: every χ ≠ χ₀ mod p is
  primitive, so R = max|Σ Λχw| ≤ (x^B Z(L) + x^{1/2}Z_L(L))/2 + W(0)/(1 − e^{−2Y}) and Σ_{p∤m} Λw ≥ xW(1) − … (Lemma 1.2).
- **Shifted sieve** (Iwaniec 1978 Lemma 1, as used by Martin; re-derived as Lemma S2): any lower-bound sieve weights on the first ω
  primes transport to the actual prime factors of p − 1, with main term ≥ the first-ω-primes main term T_ω and cost C_ω. This makes the
  worst case p − 1 = P_ω·k and reduces everything to one finite sieve problem per ω.
- **ω ≤ 10⁴:** Rosser's linear-sieve weights at level y (any prefix-closed set with even boundary works, Lemma S1), enumerated exactly
  (sieve.py; main terms as exact integers in units 2⁻⁶⁴); the best of 79 levels per ω; (C − 1)/T ≈ 0.14 p_ω², optimal y ≈ p_ω^{1.9–2.0}.
  Then one Arb inequality per ω at L = log(P_ω + 1), extended to all larger L with frozen kernel constants and supremum
  majorants R̂ for the error side (Theorems A, B).
- **ω > 10⁴:** a sieve with an analytic main term: exact Möbius below 50, Buchstab at every prime q ≤ z, Selberg's Λ² upper bounds
  at level (z^{2.3}/q)^{1/2}; G_q(ξ) bounded below three ways (log min(ξ, q); first-order inclusion–exclusion with Rosser–Schoenfeld;
  Rankin), numerically for q ≤ 10⁵ and uniformly beyond; Abel summation with explicit Mertens errors; 205 cells in log z plus an
  exponential-decay tail give T ≥ 0.1014·Π_{p≤z}(1 − 1/p) for every z ≥ 104743 (Theorem C).
- p ≤ 10⁷: exact g*(p) (664578 primes).

## Certificate
One command (≈ 4 min, one thread, < 0.3 GB; python-flint Arb; every inequality a certified ball sign, else abort):

    cd Riemann/Deep/PrimRootNum && uv run --no-project --with python-flint --with mpmath --with numpy python3 cert.py all

Output: logs/cert_output.txt (also logs/kappa_table.txt, logs/tail_output.txt). Inputs re-used by path from ExplProg/cert.py.

## Checks
- sanity.py: Lemma S1 by brute force over all n | P_ω (ω ≤ 13, 11 levels): 0 violations; Lemma S2's main-term inequality in exact
  rationals for 200 random prime sets per (ω, y), 40 pairs: all ≥ 0; sieve.py's bulk enumeration against a plain DFS (ω ≤ 200):
  identical counts and main terms (logs/sanity_output.txt).
- Exact g*(p), 3 ≤ p ≤ 10⁷: max 211 (p = 4022911), records listed in DERIVATION sec. 6; all ≤ 211 < (0.814 log 5)^{3/(1−B)} = 601.44 (Arb), so the bound holds for 5 ≤ p ≤ 10⁷ (float scan worst
  ratio of logs 0.108 at p = 5).
- adv.py: least primes ≡ 1 mod P_m (p − 1 primorial-like, m ≤ 45, p up to 10⁸¹), Sophie Germain, Fermat and 2^a3^b + 1 primes:
  g*(p) ≤ (κ_ω log p)^{7.9} everywhere (worst log ratio 0.083); Rosser–Schoenfeld R1–R5 checked at every prime ≤ 10⁷ (logs/adv_output.txt).
- Three focused referee passes (read-only reviewer subagents on DERIVATION secs. 1–2–4, 5 lemmas, 5–7). Findings, all fixed in the text
  (the code already did the right thing in the first three): endpoint tests (4.4) and (5.6) stated with explicit supremum majorants R̂_m, R̂₁ rather
  than a monotonicity claim for Ro (the code always used them); (5.3) only for ξ > x₀; Rankin candidates only with r_c < 1; u₁ written as 11.55926… not 11.5593; small-p bound
  now proved by max g* = 211 < (0.814 log 5)^{3/(1−B)} (Arb) plus monotonicity instead of a float scan; exponents written with ≤ and
  rounded up; decimal thresholds made sufficient (10^{41.032}, 10^{97.547}); normalised u ≥ 60 bounds displayed. A fourth pass on
  tail.py/cert.py against DERIVATION sec. 5 found only interval-comparison hazards (max/min selections between overlapping balls,
  ψ's branch tests on balls); fixed by comparing exact endpoints and deciding branches on exact rationals, then re-run (same numbers).
- Self-review (DERIVATION sec. 9): no GRH, no zero density, no unrefereed zero count; all roundings upward (κ, c) or downward (τ);
  monotonicity claims argued; smallest p (3, 5, 7: conductor < 9 handled by Lemma Z at log 9); p = 3 excluded only from the 0.814 form.

## Novelty
NOVELTY.md. Searched arXiv math.NT/CO listings 2026-10 and 2026-11 (lst.py), Zenodo ('"least prime primitive root"',
'"prime primitive root"', 'primitive root explicit bound'), the seven dr vengy notes (descriptions and full text: no constant
computed; 23623384 only g*(p) < p), Martin 1997 (TeX), GXL evidence API (2 claims), Hexagon, OpenAI's primitive-roots preprint (our
audit), Johnston–Naik 2610.08569 (different question), ArtinHalfPlane, our repo. No numeric polylogarithmic bound for g*(p) exists;
under GRH the explicit bounds carry 2^{ω(p−1)}. Corollary E's g*(p) < p for p > 10⁴¹ reproduces Zenodo 23623384 (QRH-conditional)
and is a cross-check, not a claim; g*(p) < √p − 2 for p ≥ 10^{97.547} has no GRH-free predecessor.

## Caveat
- Conditional only on H, which is Lean-kernel-checked (tracker) but not refereed. Not refereed; no external reviewer yet.
- The exponent 3/(1−B) = 23.693 is proved for ω(p − 1) ≤ 10⁴ only; for ω(p − 1) > 10⁴ (p > e^{104403}) the exponent is 3.4/(1−B) = 26.852.
  **Blocking lemma** for 23.693 everywhere: an explicit lower bound T_ω ≥ c/(log p_ω)² for Rosser's lower-bound weights at level
  y ≍ p_ω² valid for every ω (an explicit form of Iwaniec's Lemma 2 / the beta-sieve f(s) near s = 2); here it is computed exactly only
  for ω ≤ 10⁴, and the analytic Buchstab–Selberg substitute needs level p_ω^{2.3} (with our lower bounds for Selberg's G its main
  term stays positive only for s ≳ 2.25).
- Constants are far from optimal: sieve levels come from a coarse grid, the weights are Rosser's (not optimised), and the
  character-sum constant is ExplProg's (72% of it from zeros below height 5/7); real g*(p) are tiny (≤ 211 below 10⁷).
- g(p) ≤ g*(p), so every bound also holds for the least primitive root g(p); no better exponent for g(p) falls out (sums over
  integers rather than primes would need growth bounds for L(s, χ) inside the strip, which give only p^ε-type results).
- Rosser–Schoenfeld 1962 inequalities are used as published (ranges x ≥ 286, θ: x ≥ 41) and checked numerically only up to 10⁷.

RH relevance: N11

## Files
- REPORT.md (this file), DERIVATION.md (all steps, constants table sec. 10), NOVELTY.md.
- cert.py (certificate: core, sieve tables, Theorems A/B, small primes, tail, uniform c, corollaries), sieve.py (exact sieve weights),
  grid.py (choice of levels → logs/sieve_best.json), tail.py (sec. 5), small.py (exact g*(p)), sanity.py (Lemmas S1, S2), adv.py.
- logs/: cert_output.txt, cert_time.txt, kappa_table.txt, tail_output.txt, sieve_best.json, grid_run.txt, sanity_output.txt, adv_output.txt.

Suggested DISCOVERIES line: "Least prime primitive root g*(p) ≤ (1.956 log p)^{26.86} for every prime; (0.814 log p)^{23.70} when ω(p−1) ≤ 10⁴ ([report](Riemann/Deep/PrimRootNum/REPORT.md))".
