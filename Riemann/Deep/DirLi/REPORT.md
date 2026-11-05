RECORD: every primitive Dirichlet character χ mod q ≤ 400,000 has Re λ_n(χ) > 0 for all 1 ≤ n ≤ 26·t0(q)², t0(q) = Platt's GRH height (q ≤ 1000: n ≤ 41·10^16/q²; q = 3: n ≤ 6.05·10^16), from refereed inputs; with the 7/8 half-plane n ≤ 35·t0(q)² (q ≤ 1000: 55·10^16/q²; q = 3: 8.16·10^16) — first explicit certified Li ranges for Dirichlet L-functions

# dirli report (node N3, Li / Keiper–Li coefficients for Dirichlet L-functions), 2026-11-05

## Short answer

- **What is new.** Explicit, certified positivity ranges for the generalized Li coefficients
  λ_n(χ) = Σ_ρ [1 − (1 − 1/ρ)^n] of Dirichlet L-functions. Nothing explicit existed: Omar–Ouni–Mazhouda (arXiv 1507.03431)
  tabulate n ≤ 40 from 10^4 zeros assuming RH (not certified) and their Prop. 1 ("n ≤ T² for T ≥ T₁ sufficiently large") has
  no explicit T₁. Our ranges come from Platt's GRH verification to height t0(q) (10^8/q for q ≤ 125,000) plus an explicit
  zero count, and reach N(q) ≈ 27–54 t0(q)² (refereed inputs) and 37–73 t0(q)² (adding the 7/8 half-plane), slowly
  decreasing in q (41–54 and 56–73 for q ≤ 1000).
- **Headline numbers (certified, every primitive χ of the modulus, both parities, complex χ via Re λ_n):**

| q | t0(q) (Platt) | N(q), refereed inputs (a = 1/2) | N/t0² | N(q), with 7/8 half-plane (a = 3/8) | N/t0² |
|---|---|---|---|---|---|
| 3 | 3.33·10^7 | 60,469,357,461,824,346 ≈ 6.05·10^16 | 54.42 | 81,573,243,404,167,256 ≈ 8.16·10^16 | 73.42 |
| 4 | 2.5·10^7 | 33,613,820,289,736,842 ≈ 3.36·10^16 | 53.78 | 45,352,014,948,651,417 ≈ 4.54·10^16 | 72.56 |
| 5 | 2·10^7 | 21,313,957,270,270,402 ≈ 2.13·10^16 | 53.28 | 28,760,439,636,336,494 ≈ 2.88·10^16 | 71.90 |
| 7 | 1.43·10^7 | 10,721,266,235,936,412 | 52.53 | 14,469,688,200,660,378 | 70.90 |
| 8 | 1.25·10^7 | 8,161,869,755,212,833 | 52.24 | 11,016,303,953,316,466 | 70.50 |
| 100 | 10^6 | 46,556,556,033,802 | 46.56 | 62,943,150,794,887 | 62.94 |
| 1000 | 10^5 | 413,010,283,643 | 41.30 | 559,492,859,836 | 55.95 |

  All 749 moduli 3 ≤ q ≤ 1000 with primitive characters were certified one by one (results_q1000.txt); N(q) and N(q)/t0(q)²
  decrease monotonically in q, so **for all q ≤ 1000: Re λ_n(χ) > 0 for n ≤ 41.30·10^16/q² (a = 1/2) and
  n ≤ 55.94·10^16/q² (a = 3/8)**. Moduli 1000 < q ≤ 400,000: 604 blocks [q1, q2] (q2 ≤ 1.02 q1, odd moduli and
  moduli ≡ 0 mod 4 separately; all 299,250 moduli covered), one certified N per block (results_blocks.txt), e.g.

| moduli (class) | t0(q2) | N, a = 1/2 | N, a = 3/8 |
|---|---|---|---|
| 1001–1021 (odd) | 97,943 | 395,722,709,413 | 536,085,524,385 |
| 9993–10191 (odd) | 9,812.6 | 3,455,030,662 | 4,693,301,959 |
| 99597–101587 (odd) | 984.4 | 29,283,017 | 39,936,837 |
| 123844–126316 (≡ 0 mod 4) | 793.75 | 18,683,349 | 25,492,899 |
| 307988–314144 (≡ 0 mod 4) | 438.7 | 5,395,272 | 7,372,374 |
| 398419–399999 (odd) | 293.75 | 2,319,928 | 3,173,606 |
| 398420–400000 (≡ 0 mod 4) | 387.5 | 4,155,533 | 5,680,214 |

  Uniformly, **for every q ≤ 400,000 and every primitive χ mod q: Re λ_n(χ) > 0 for n ≤ 26.12·t0(q)² (a = 1/2) and
  n ≤ 35.73·t0(q)² (a = 3/8)** (minimum of N_block/t0(q1)² over the blocks, at odd q ∈ [301949, 307985]; t0 decreases in q
  within a parity class, so N(q) ≥ N_block ≥ c·t0(q1)² ≥ c·t0(q)²; the per-q minimum for q ≤ 1000 is 41.30 resp. 55.95).
- **Sharpness.** The gain bound is asymptotically the GRH main term (n/2)(log n + log(q/2π) − 1 + γ), checked against
  directly computed λ_n(χ) (sanity.py); the range is set by the worst case "every zero above t0 at distance a from the line".
  R315's uncertified estimates (6.1·10^16, 4.1·10^11, 4.1·10^6 for q = 3, 1000, 400,000) are confirmed
  (6.05·10^16, 4.13·10^11, 4.16·10^6).
- **RH relevance: N3 — none beyond the inputs.** As for zeta (Li200k, Fabbian), finite Li ranges are a reformulation of
  "GRH to height t0 + strip width": N ≈ (t0²/a)·log(…), no information about zeros above t0.

## Theorem

Let χ be a primitive character mod q, 3 ≤ q ≤ 400,000, t0(q) Platt's height, a = 1/2 (zeros above t0 only known to lie in
0 < β < 1) or a = 3/8 (7/8 half-plane). Then Re λ_n(χ) > 0 for 1 ≤ n ≤ N_a(q), with N_a(q) as certified in
results_q1000.txt (q ≤ 1000) and results_blocks.txt (q > 1000; one value per block of moduli of one parity class).
(For real χ, λ_n(χ) is real; for complex χ, λ_n(χ̄) = conj λ_n(χ). GRH for L(s, χ) ⟺ Re λ_n(χ) ≥ 0 for all n.)

## Proof (summary; full text PROOF.md)

1. Re λ_n(χ) = Σ_ρ Re(1 − z_ρ^n), z_ρ = 1 − 1/ρ, converges absolutely. A zero on the line contributes ½G_n(|γ|),
   G_n(t) = 2 − 2cos(nθ(t)), θ(t) = 2arctan(1/2t). Off-line zeros pair as ρ, ρ* = 1 − ρ̄ (same ordinate; the zeros of a
   complex χ are not conjugation-symmetric, so pairs, not quartets), and a pair contributes ≥ G_n(|γ|) − q_n(|γ|),
   q_n(t) = 4sinh²(na/2t²) + 4na²/t³ (Li200k's per-zero bound, re-derived). Hence Re λ_n ≥ ½M_n − ½Q_n with
   M_n = Σ_ρ G_n(|γ|) (ordinates only) and Q_n = Σ_{|γ|>t0} q_n(|γ|).
2. N(t, χ) = #{|γ| ≤ t}: Bennett–Martin–O'Bryant–Rechnitzer (BMOR) Thm 1.1, |N − (t/π)log(qt/2πe) + χ(−1)/4| ≤
   0.22737ℓ + 2log(1+ℓ) − 0.5, ℓ = log(q(t+2)/2π) > 1.567, t ≥ 5/7. Stieltjes integration with this error band (the parity
   constant cancels) bounds M_n from below in three regimes: per n for small n (integer counts, both parities), zeros
   above t_n (n < 10^4), and all zeros above √n cut into monotone pieces (n ≥ 10^4); Q_n from above by an upper Riemann sum.
3. Block chains [n_i + 1, n_{i+1}] with a certified sign gain(n_i+1, n_{i+1}) − loss(n_{i+1}) > 0 (python-flint balls;
   floats only choose block ends) cover [1, N]. q-blocks use worst-case q in each monotone ingredient and t0 of the
   largest modulus in the block (t0 decreases in q within a parity class of q).

## Inputs

| input | used for | status |
|---|---|---|
| Platt, Math. Comp. 85 (2016), arXiv 1305.3087, final theorem | GRH for all zeros of height at most t0(q) (both signs), every primitive χ mod q ≤ 400,000; t0 = max(10^8/q, 7.5·10^7/q + 200) (q even), max(10^8/q, 3.75·10^7/q + 200) (q odd) | refereed |
| BMOR, "Counting zeros of Dirichlet L-functions", Math. Comp. 90 (2021), **arXiv 2005.02989** | N(t, χ) error | refereed |
| classical 0 < β < 1 | a = 1/2 | classical |
| OpenAI 7/8 half-plane for every Dirichlet L-function (Lean port Riemann/Deep/QRH; QRHCor/Dirichlet.lean) + functional equation | a = 3/8 | unrefereed, Lean-checked |

Note: the orchestrator's id arXiv 1812.11592 is "A Geometric Theory of Higher-Order Automatic Differentiation"; BMOR is
2005.02989. Trudgian's 2015 Dirichlet zero count (reported error, Hasanalizade–Shen–Wong 2022 fn. 2) is not used.
Li200k's zeta certificate used Trudgian 2014 for zeta (same reported error); per the orchestrator it should be re-run
with HSW 2022 Cor. 1.2 (outside this brief).

## Certificate (Riemann/Deep/DirLi/certify.py; python-flint 0.9.0, 128-bit balls, one thread)

- part0: the pair inequality on 4000 random off-line zeros (both signs of γ, |γ| ∈ [10, 10^4], a ∈ {1/8,…,1/2}); min slack 9.9·10^−6.
- small: q = 3, 4, 5, 7, 8, both variants, 3.6 s, 41 MB peak.
- q1000: all 749 moduli 3 ≤ q ≤ 1000, both variants, ≈ 8 min CPU (one process; wall time longer because the machine slept).
- blocks: moduli 1000 < q ≤ 400,000, 604 blocks of ratio ≤ 1.02 per parity class, both variants, 291 s; coverage of all
  199,500 odd and 99,750 even moduli checked. Consistency: blocks on (900, 1000] give N_block ≤ min per-q N in every block
  (within 3·10^−5), and one-modulus blocks (q = 999, 1000) reproduce the per-q values exactly.
- Bug found and fixed during the run: a relative-step stop rule ended the chain right after the small-n regime for 8 moduli
  (15, 16, 17, 75, 76, 77, 615, 617), giving tiny but true N; the rule now applies only for n > 10^5; all q rerun.

## Sanity checks (Riemann/Deep/DirLi/sanity.py; not part of the proof; 36 s)

- Normalisation: zeta λ_1 = 0.0230957089661, λ_2 = 0.0923457352280 from the Taylor series of log ξ at s = 1.
- Re λ_n(χ) for n ≤ 300 from the Taylor series of log ξ(s, χ) at s = 1 (complex χ: the real part is unaffected by the
  imaginary linear Hadamard factor), every primitive χ mod 3, 4, 5, 7, 8, 11, 12 (22 characters): our certified lower bound
  is below Re λ_n(χ) for every n (max ratio 0.835); Re λ_300(χ) is within 0.8% of (n/2)(log n + log(q/2π) − 1 + γ)
  (Lagarias' constant; Omar–Ouni–Mazhouda's c_χ with log(q/π) is off by (log 2)/2).
- χ mod 3: 189 zeros up to height 300 from sign changes of Z(t, χ) (count matches BMOR's main term 378.8 for both signs);
  Σ_{|γ|≤300}(1 − cos nθ) plus the expected tail reproduces the Taylor λ_n(χ) to 10^−4 relative for n ≤ 5 (λ_1 = 0.0566150).
- A first sanity run flagged "bound above λ_54"; the cause was the sanity script's default series truncation (ctx.cap = 10),
  not the bound; with the full series length the check passes.

## Novelty check

arXiv API (Li coefficients/criterion with Dirichlet, L-functions, Keiper, height, zero-free), arXiv monthly listings
math.NT and math.CV 2026-10 and 2026-11 (no Li-coefficient titles; 2610.38870 is Jensen–Pólya, 2610.13565 Appell
polynomials), Zenodo API ("Li coefficients"/"Li criterion"/"Keiper-Li" with Dirichlet; Fabbian; GRH Li): Fabbian 22839423
(zeta only), Fabbian 23027727 (general ξ, analyses the shape of Brown's 2005 converse-direction proofs, no ranges),
Pavesi 22944511 (τ-moments for Dirichlet L, no ranges), Chao 23370435 (zeta, Platt–Trudgian height, no computation).
Omar–Ouni–Mazhouda 1507.03431 (tables n ≤ 40 under RH, non-explicit Prop. 1). Brown 2005 (J. Number Theory 111) relates
Li-type positivity and zero-free regions for a general class containing Dirichlet L-functions; its printed proofs rest on a
lemma found defective (Fabbian 23027727); no explicit Dirichlet numbers are known to us.
GXL evidence API: daily quota exhausted (HTTP 429), not run. Not found elsewhere.

## Review

Self-review twice (PROOF.md §2–§7: pair inequality, absolute convergence, multiplicities and negative ordinates, Lemma 2
both directions, the factor 2 in g = 2(1 − cos) against the halving, piece count, monotonicity in n and in q).
Reviewer subagent (verdict only, proof text inline): "VERDICT: CORRECT" for (1) pair inequality, (2) Re λ_n ≥ ½M_n − ½Q_n,
(3) Stieltjes lemma and parity cancellation, (4) gains A, B, C, (5) loss and q-block worst cases; remarks: state Lemma 2
with the measure −dH for general monotone H (all uses are smooth), and the asymptotic sharpness of B holds for n1 ~ m.

## Classification and suggested repo edits (orchestrator)

- Method routine (Fabbian/Li200k transferred to Dirichlet L with BMOR); numbers new. Suggested: DISCOVERIES side result or
  COROLLARIES, "explicit Li positivity for all primitive Dirichlet L mod q ≤ 400,000; q ≤ 1000: n ≤ 41·10^16/q²
  (55·10^16/q² with 7/8); q = 3: 6.05·10^16 / 8.16·10^16".
- TREE N3 record line: add "Dirichlet L: Re λ_n(χ) > 0 for n ≤ 26 t0(q)² for all primitive χ mod q ≤ 400,000
  (q = 3: 6.05·10^16; 35 t0(q)² and 8.16·10^16 with 7/8), from Platt's GRH heights (dirli, R###)".

## Files

- Riemann/Deep/DirLi/REPORT.md — this report.
- Riemann/Deep/DirLi/PROOF.md — derivation.
- Riemann/Deep/DirLi/certify.py — certificate (modes part0, small, q1000, blocks).
- Riemann/Deep/DirLi/sanity.py — numerical checks.
- Riemann/Deep/DirLi/results_q1000.txt, results_blocks.txt — certified N per modulus / per block.
