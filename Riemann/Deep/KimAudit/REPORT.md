AUDIT PASSED — items 1–3 of jensen-kim-lee-wedge and its a = 1/2 line hold as written: every step re-derived, the cited inputs match their sources, the Arb certificate reruns byte-identically and certifies exactly the stated bounds, and no independent recomputation or falsification test found a discrepancy.

# kimaudit: line-by-line referee of jensen-kim-lee-wedge (N6)

Orca agent kimaudit, 19 Nov 2026 (22:21–23:25 UTC). Object: [discoveries/jensen-kim-lee-wedge.md](../../../discoveries/jensen-kim-lee-wedge.md),
Riemann/Deep/KimEff/ (DERIVATION.md, certify.py, logs). The chain with source line numbers: [FULLPROOF.md](FULLPROOF.md). No Lean
(kimefflean not audited); KimEff/ untouched.

## Verdict per item

| item | claim | verdict | evidence |
|---|---|---|---|
| 1 | Θ^{(n)} ≠ 0 on W_{R_n} = P ∩ {Re w ≥ −R_n}, n ≥ n_* | holds | D §2, §4–§9 re-derived (FULLPROOF §2, §4–§9); certify rerun identical; own majorants equal at u = 12 and monotone on 5766 points; reviewer on §8.5: CORRECT; at n = 60 and 1000, W_{R_n} has no zero |
| 2 | J^{d,n} hyperbolic for d ≤ D_n, hence for d ≤ 6.45 n²/log² n | holds | D §3 (parabola → sector, Kim–Lee Thm 3, Hurwitz) re-derived; KL Thm 3 = (I3) verbatim; D_n/(6.45 n²/log² n) ≥ 1.0155 on n ∈ [n_*, 10^297] |
| 3 | n ≥ max(n_*, C d^c), four (c, C) | holds | closed form, 40,000-point grid sup and 2404-point direct sampling: D_n ≥ d in every case |
| a = 1/2 | d ≤ 3.59 n²/log² n | holds | certify_half rerun identical; ratio ≥ 1.0175 |
| numbers | n_* = 998,616,377,771; D_n > D_all for n ≥ 4.6974173·10^13 | correct | Z(12) = 998,616,377,770.9088; exact crossover 46,974,172,465,412.6 |

## Cited inputs
- (I1) Riemann's integral (Titchmarsh 2.16.1–2, Φ = 2Σ…): ξ(1/2 + z) = 4∫φ cosh(zt) ✓; γ(m) = 32 m! Mt(m)/(2m)! equals
  8 m! ξ^{(2m)}(1/2)/(2m)! numerically (m = 0, 1, 2; mpmath ξ) and Θ(w) = 8ξ(1/2 + √w) at w = −50 + 10i ✓.
- (I2) a = 0.37338196661320752 = 0.87338196661320752 − 1/2: TREE.md l.11 (N7) records the Lean-checked half-plane
  Re s > 0.87338196661320752 (Bebel, 31 Oct) and the later, smaller Lean-checked bounds (B_2 + 2.2·10⁻²⁵⁴, B_2, 0.8733468 (PR 137),
  θ₁₄₁ = 0.87334585… (PR 141)); each implies (I2) by the symmetry ρ ↦ 1 − ρ̄. Refereed as an input by HPAudit; not re-audited here.
- (I3) Kim–Lee Thm 3 (KL TeX l.171–178) and its Corollary (l.183–189) are D l.42–45 verbatim: real P, Q; Z(P) ⊂ S(δ) = {|Im z| ≤ δ|z|};
  Q hyperbolic; deg Q ≤ δ^{−2} ⇒ P(D)Q hyperbolic; J(P;d)(z) = z^d (P(D)z^d)(1/z). KL cite Obreschkoff 1961 for the proof (not seen
  here). Random test: 0/800 failures inside, 173/800 in the control one degree outside; D's sharpness example confirmed.
- (I4) DLMF 10.21.15, 10.16.9, 10.21.3 ✓ (ν = n − 1/2 ≥ 1/2); S1 = 1/b, S2 = 1/(b²(b+1)) are the Rayleigh sums ✓. (I5) Hadamard
  (order 1/2 ⇒ genus 0) and Hurwitz ✓.

## Code audit (certify.py)
- Each majorant (C l.62–85, l.104) equals the D formula; re-implemented from D's text (mono.py): same values at u = 12 to 10 digits.
- Directed rounding: up()/lo() step a decimal until the Arb comparison holds; checks are whole-ball strict; 256 bits ✓.
- Certified: peak bracket (M1, Q2), local interval (M3a/b), phase including the largest |Im w| (vmax = 2a√(R_n + a²) in Φ_up),
  side tails (E_up), crude tail (Hc concavity, Hcp_low ≥ 1, c6 > 1/(2n), logratio), FINAL ≥ 1.8747279 ✓.
- Thresholds (lesson 277): the proof of (i) has no scan: one evaluation at u_* = 12 ≤ u(n_*) = 12 + 4.4·10⁻¹⁴, then analytic
  monotonicity (D §9, checked by hand and on the grid); H2 and H3 at n_* exactly; D(n) uses a certified upper bound of u(n), so R_n and
  D_n are lower bounds; the crossover bisection keeps a certified point ✓.
- Rerun on a copy in /tmp: both logs byte-identical, "ALL CHECKS PASS" (logs/certify_rerun.log, logs/certify_half_rerun.log).

## Reviewer (most delicate step: D §8.5 crude tail + §8.6 assembly, formulas pasted inline)
- Run 1: (A)–(H) each correct (the 2π² cancellation, 1/|Hc'(U)| ≤ 1, the e^{h(u2)} normalisation, U²Wmax/(4bn) ≤ θU/4 + the vmax
  term, −LG/n ≤ θ/4 + θ²/(32(1−θ)), −log η ≤ log(5n/K)/2, F' = ([2(ρ−1)u − 1]e^{2y} + 1 + 2ρ)/u² > 0, d(logratio)/dn < 0);
  verdict GAP, caused by my prompt only: I wrote c6 ≥ 0.786 (D l.186: 0.78582827) and omitted φ ≥ φ_1 and the bracket constants.
- Run 2 (closure, D's exact statements): CORRECT. Its precision note (c6 from the rounded F and G_up alone gives 0.78582819) is moot:
  certify.py computes c6 as one Arb ball; mono.py gives c6 = 0.7858282718.

## Remarks (non-blocking; no claim changes)
1. D l.207 "(≥ 4.448)" is the margin printed by C l.122 (value − 1); the value is 5.448 ≥ 1. Wording only.
2. D l.67 assumes R > 2a²; not needed (δ ≥ 1 gives D(R) ≤ 1, trivially fine).
3. C l.157 prints the non-certified side of the two-thirds bisection ("needs n ≥ 3.41725e37"); a comparison number, not in the
   theorem; recomputed 3.4172529·10^37.
4. Analytic steps outside the code, checked by hand: the k ≥ 12 tail of ε0, 4n/u + 9 ≤ 5n (−log η bound), every §9 direction.
5. What the result still rests on: (I2) the Lean-checked half-plane (HPAudit's referee, not redone), (I3) Obreschkoff as stated by
   Kim–Lee (original 1961 paper not seen; random test consistent, statement sharp), and Arb (python-flint) correctness.


## Dependents (rg -i 'kimeff|kim-lee-wedge|jensen-kim|6\.45 n'): none needs a change
TREE.md l.10 (N6 line, explicit Kim–Lee wedge entry); RECORDS.md l.22 (record, n ≥ 9.99·10^11, d ≤ 6.45 n²/log² n); FORMALIZATIONS.md
l.11 (kimefflean, partial Lean with (i) as a named hypothesis); Riemann/Tree/R423.lean (the round), R413.lean (no effective Kim 1996);
Riemann/Deep/HP141/KimEff (rerun with a = θ₁₄₁ − 1/2: unchanged at printed precision); the page itself. DISCOVERIES.md and
COROLLARIES.md: no entry.


## Falsification (own mpmath code; FULLPROOF §12)
- n_*, D_n, crossover, (c, C), comparisons recomputed from D's formulas (indep.py): all agree (numbers in the table above).
- §9 monotone majorants re-implemented from D's text (mono.py), 5766 grid points, both a: no violation, minima at u = 12.
- Zero-free parabola at moderate n from the Taylor coefficients γ(n+m) (theta_n.py, argument principle): n = 60: W_{R_n} zero-free,
  first zero at 7.09 R_n; n = 1000 (300 digits): W_{R_n} zero-free (reliable), W_{2R_n} zero-free (indicative). No violation.
- Mixture identity at n = 60 (mixture.py): Taylor = Bessel mixture to 28 digits on W_{R_n}.
- J^{d,n} for n = 3, 5, 8, 10 at d = ⌊D_n⌋, ⌊D_n⌋ + 1, 2⌊D_n⌋ (jensen_small.py): 12/12 hyperbolic (outside too, as expected: these d are
  inside the all-n range); γ normalisation equals mpmath's ξ.
- Kim–Lee Thm 3 (obresch.py): 0/800 random failures inside, 173/800 in the control, sharp example confirmed.

## Files (Riemann/Deep/KimAudit/)
REPORT.md, FULLPROOF.md; indep.py, mono.py, theta_n.py, jensen_small.py, obresch.py, mixture.py; logs/ (indep, mono, mono_half,
theta_n60, theta_n1000 (+ run1), jensen_small, obresch, mixture, certify_rerun, certify_half_rerun; .err files hold memcap peaks,
all ≤ 48 MB). Every script: uv run --no-project --with mpmath python3 <script> under .scratch/memcap.sh 290.

RH relevance: N6
