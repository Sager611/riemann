AUDIT PASSED: every item of the isolated-zeros claim holds as written — the chain (Lamzouri's Hilbert-space inequality with explicit errors, plus the new hinge lemma from Conrey's published mean value) re-derived step by step against the TeX, all constants recomputed independently, the certified values at every stated threshold reproduced, and no counterexample found.

Auditor: isoaudit (Orca), 19 Nov 2026, 21:49–22:16 UTC. Target: discoveries/isolated-zeros.md, Riemann/Deep/IsoZero/ (DERIVATION.md, cert.py). Proof chain with source lines: [FULLPROOF.md](FULLPROOF.md) §A–H.

## Items (verbatim from the page) and status
| item | status | evidence |
|---|---|---|
| For every c ≤ 0.00476, liminf B_0(T,c)/N(T) ≥ 0.6724, unconditionally | HOLDS | certified 0.6724005185 at c = 0.00476 (margin 5.2e-7); B_0 non-increasing in c (FULLPROOF F7, G4); inputs unconditional (E1, E2) |
| B_0 = zeros 0 < γ ≤ T whose two neighbouring ordinates are > 2πc/log T away (Lamzouri Thm 1.1) | HOLDS | [La] l.211–212 (the page omits "n ≥ 2"; harmless) |
| such zeros are simple and on the line; counted with multiplicity | HOLDS | ρ and 1 − ρ̄ share the ordinate (A1); [La] l.587–600 |
| c = 0.001: certified 0.672487 | HOLDS | certified 0.6724872882; float 0.672487289 (F2) |
| ≥ 0.672 for c ≤ 0.01106 | HOLDS | 0.6720007083; float c_max 0.01106798 |
| ≥ 0.67 for c ≤ 0.02485 | HOLDS | 0.6700018526; c_max 0.02485920 |
| ≥ 0.66 for c ≤ 0.05495 | HOLDS | 0.6600054000; c_max 0.05496151 |
| Lamzouri's 0.6552 for every c ≤ 0.06427 | HOLDS | 0.6552053893; c_max 0.06427959 |
| ≥ 0.5247 at c = 1/6 (paper 0.50875) | HOLDS | 0.5247019978 with the cos(wx) kernel, w² = 3.75 (f_0 gives 0.5167056); rests on the general form of La's Lemma 3.2 (LamzAudit2 S3′ "Form G"), not on the f_0-only statement quoted in [La] l.670–700 |
| Lamzouri: 0.6552 for an unstated c | HOLDS | [La] l.213–215; his 0.655202908 and 0.508751056 reproduced (E5) |

## What was checked
1. Derivation line by line ([D] §0–4, 6, 7; §5 via c = 1/6; §8 not in scope): every hypothesis, range, constant, normalisation and inequality direction against [La] l.258–616, 640–700, 740–822 and BCY 1002.4127 l.205–240, 337–347 (FULLPROOF A–E). Lamzouri's (E) checked including the signed term −h⊗h (handled by Bessel, l.393–402). The finite inequality (F) behind every value is re-derived (C5); the hinge lemma (D1–D6) and its two independent normalisation checks: R → 0 reproduces Conrey's announced 0.0806 and BCHB's 0.0845, large R reproduces Conrey's density exponent 8/7 (E3).
2. c-dependence enters only through λ_C (K(2c), K(4c)), θ̄ = π²c²/2 and ρ_max(c), α_ρ(c); C(f) and B(ρ) are c-free. Thresholds need only one certified value each, by monotonicity of B_0 in c (F7).
3. Code (G): cert.py certifies exactly (H); outward rounding sound (one immaterial 1e-16-level float blemish); side conditions asserted; thresholds decided by exact certified comparisons.
4. Independent numbers (F): own mpmath code, K by re-derived closed forms checked against quadrature, C(f) by double quadrature; all table values and thresholds reproduced to ≤ 1e-9.
5. Falsification (H): 128,000 + 10,500 random/adversarial multisets at corner parameters (c up to 1/4 = the 4c ≤ 1 limit, three kernels): no violation of (F); first 10⁵ real zeros: B_0/N ≥ 0.9967 for every stated c.
6. Reviewer on the most delicate step (hinge lemma, pasted inline): GAP on presentation only — see Issues 1–2; no step or number changes.

## Issues (none load-bearing)
1. [D] l.133 defines α_ρ for 0 < ρ ≤ ρ_max and argues finiteness by "h(ρ) = 0 is approached from below"; that is no proof at ρ = ρ_max for general continuous h (reviewer's counterexample). Harmless: H is real-analytic, and cert.py uses only ρ ≤ ρ_max − 0.02, where α_ρ is certified finite. Fix: state 0 < ρ < ρ_max.
2. Pairing and Littlewood details (reviewer): pairs are ρ ↔ 1 − ρ̄ (as [D] §0 says); Jensen with T − 1; zeros on the rectangle's boundary by limits. o(T)-level.
3. logs/cert.log prints "p(c0) >= …" rounded to nearest (%.6f), overstating the certified value by up to 5e-7 (0.672401 for 0.6724005185). The claims are unaffected; the exact values are in logs/certvals.log.
4. The c = 1/6 value uses a non-f_0 kernel, so it needs Lamzouri's pair-correlation lemma in its general form (shown by its proof, LamzAudit2 S3′), not as quoted in 2610.26841.
5. The method's ceiling 2 − C_MT = 0.672501 (page Caveat) checked: p ≤ 1 − (C − 1)/(1 − λ) ≤ 2 − C.

## Dependents (all survive)
RECORDS.md l.23 (N8 record); TREE.md l.12 (N8 node: isolated zeros and the linear-combination result); discoveries/isolated-zeros.md; discoveries/linear-combinations-line.md (2p − 1 ≥ 0.344982508 uses this p(c) at c = 10⁻⁴; its transfer to L(s,χ) via Wu's Prop. 3.1, [D] §8, is outside this audit); discoveries/linear-combinations-three.md (imports the hinge bound); Riemann/Tree/R392, R395, R406.

## Files
FULLPROOF.md (chain, source lines, all checks); indep_lib.py, indep.py (mpmath recomputation, logs/indep.log, 331 s, peak 49 MB); certvals.py (exact certified values via the unmodified cert.py, logs/certvals.log, 11 s, 15.8 MB footprint); fintest_lib.py, fintest.py, fintest2.py (finite-inequality falsification, logs/fintest.log, logs/fintest2.log); zcheck.py (real zeros, logs/zcheck.log; Odlyzko zeros1, 1.8 MB, downloaded to .scratch and deleted; ran 2 s, below memcap's first sample); logs/review.log (reviewer verdict and resolution).

RH relevance: N8
