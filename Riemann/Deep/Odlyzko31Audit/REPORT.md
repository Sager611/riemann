AUDIT PASSED: every mathematical and numerical item of discriminant-bounds-strip holds as written, given Cureton–Thorner's Lean theorem (and Bebel's for the abelian constants) as verified in HPAudit; only three wording slips (record caveat, table description) were found, none touching a bound.

# odlaudit: line-by-line audit of "unconditional root-discriminant bounds beyond the Odlyzko–Poitou limits"

Scope: discoveries/discriminant-bounds-strip.md, Riemann/Deep/Odlyzko31/ (DERIVATION.md, certify.py, kernels.json, certified.json, logs). Read-only on both. Chain with line numbers: [FULLPROOF.md](FULLPROOF.md). Code and logs: this directory. Started 19:44 UTC, 19 Nov 2026.

## Verdict per item

| Item (page "Claim") | Verdict | Evidence (this audit) |
|---|---|---|
| Explicit formula (P), Γ_R and Γ_C terms, pole term E = Φ(0)+Φ(1), primes | holds | re-derived by hand (FULLPROOF Step 1); matches Odlyzko 1990 survey eq. (2.3) term by term; identity to 5e-13…2.5e-11 for ℚ, ℚ(√−3), ℚ(√5), cubic d = −23 in independent PARI/GP code |
| H(15/32) from Cureton–Thorner | holds (input as audited) | Lean `main` read in Challenge/DedekindZetaChallenge.lean @09afbb4: existential continuation f on ℂ∖{1}, so it covers 31/32 < Re s ≤ 1 (not the junk-valued LSeries); FE reflection gives β ≥ 1/32 |
| Lemma A (FT of cosh(ax)/cosh(wx) > 0, \|a\| < w) | holds | closed form re-derived; 256-bit quadrature at 5 (a, w, t): diff ≤ 1.4e-76 |
| Lemma B (G, H autocorrelations, Ŝ ≥ 0, 0 ≤ S ≤ 1) | holds | both autocorrelations re-derived; FT = \|ĥ\|²/‖h‖² at 4 t each, diff ≤ 7e-80 |
| Prop C: Re Φ(ρ) ≥ 0 under H(w), F ≥ 0 drops primes, bound (T) | holds | convolution of nonnegative FTs; boundary \|a\| = w gives FT = bŜ(bγ) ≥ 0; reviewer: CORRECT |
| Signature and degree monotonicity (r1 = 0 value bounds all; value at n holds for all larger n) | holds | I_r(F) ≤ π/2 needs only F ≥ 0; E > 0 |
| Prop D: admissible F ≤ sech(wx); F_b → limit | holds | Bochner (φ = F·cosh(w·) continuous, integrable, φ̂ ≥ 0 ⇒ \|φ\| ≤ φ(0) = 1); dominated convergence; reviewer: CORRECT |
| A_c(15/32) = 23.3878, A_r(15/32) = 65.6439 | holds | mpmath 256-bit: 23.3878942309, 65.6439656548 (classical check: log 2, π/2 − 1, 22.3816160955, 60.8395403239 exact) |
| Abelian 26.9662, 84.2999 (w_ab = 0.3733819666132076 ≥ 0.87338196661320752 − 1/2) | holds (input as audited) | 26.9662525665, 84.2999165937 |
| rd > 22.38169 (n ≥ 2732), > 60.84021 (tr, n ≥ 1618), > 22 (n ≥ 1581), > 60 (tr, n ≥ 1223) | holds | recomputed at 256 bits: 22.3816960154, 60.8402198645, 22.0002274248, 60.0003991246 (each above its decimal; Arb balls < 1e-9 wide); scans start at n = 2 (doubling + bisection), certified exactly at the stated n |
| Certified table n = 2…10⁶ (103 degrees, every signature n ≤ 12) | holds | 40 random entries + 7 special lines recomputed at 256 bits: 0 failures, max rel. diff to certified.json floats 3.9e-15; rerun of certify.py (rerun/) |
| Finitely many fields with rd < 23.3878 (tr: < 65.6439) | holds | bound → A > threshold for large n, then Hermite |
| Falsification by fields | no counterexample | 23 minimal fields n ≤ 8, 42 LMFDB minima, 121,588 quadratic fields, their Hilbert class fields (degree 2h, rd = √\|D\|, \|D\| ≤ 2·10⁵), 149,999 cyclotomic fields and 149,997 real cyclotomic subfields (m ≤ 2·10⁵): tightest margins ℚ(√5) 0.543%, ℚ(√−3) 0.642%, Hilbert class field of ℚ(√−71) (n = 14) 3.53% |
| Side claims: "Odlyzko's 1976 table reaches 22 only at n = 8862"; abelian "elementary conductor count is stronger for large n" | hold | data/odlyzko1976_tables12.csv: unc_tc 21.845 at n = 4840, 22.021 at n = 8862; own count of primitive (even) characters by conductor: product of the n smallest conductors beats the ab table from n = 220 (20.943 vs 20.377) and, totally real, from n = 1332 (73.117 vs 72.347), as DERIVATION §5 states |
| Table sanity | consistent | every ct entry < Prop D limit A_r^{r1/n}A_c^{2r2/n} (min gap 0.2%); ≥ our classical recomputation entrywise (equal only at the formal r1 = 0 entries n = 2, 3); ≤ Odlyzko's 1976 GRH Table 1 at all 166 common rows (equal only at n = 2, r1 = 0: 1.721) ; the data file's GRH column itself spot-checked with own w = 0 code (grh_spot.py, logs/grh_spot.log): 19.2786, 59.1745, 29.1030, 5.7378 vs file 19.268, 59.141, 29.094, 5.734 at (100, 0), (100, 100), (1000, 0), (8, 0), i.e. the file holds Odlyzko's slightly rounded-down values |

## (2) Code audit of certify.py (line refs to Riemann/Deep/Odlyzko31/certify.py)
- Exact inputs: w as fmpq ('15/32', l.20-22); b parsed to an exact rational from its decimal string (l.82). All 235 ct entries in kernels.json carry wq = '15/32' (none smaller); ab entries carry 3733819666132076/10^16 ≥ 0.87338196661320752 − 1/2.
- Shapes and derivatives (l.24-36): S_G, S_H and dS_G, dS_H re-derived by hand; correct.
- I_c (l.50): [δ, 2b] by acb.integral; [0, δ] enclosed in ±Mδ with M = sup_{[0,δ]}|F'| from a ball evaluation of F' = S'(x/b)/b·sech(wx) − S(x/b)·w·tanh(wx)·sech(wx) (correct derivative), using |1 − F| ≤ Mx and 2 sinh(x/2) ≥ x; tail −log tanh(b/2) = ∫_{2b}^∞ dx/(2 sinh(x/2)) (F = 0 there). Correct.
- I_r (l.51): tail π/2 − 2 arctan tanh(b/2) = ∫_{2b}^∞ dx/(2 cosh(x/2)) (d/du 2 arctan tanh(u/2) = sech u). Correct. E (l.52): 4∫_0^{2b} F cosh(x/2). Correct.
- Rigor of acb.integral: integrands are meromorphic (poles at 2πik, iπ(k+½), iπ(k+½)/w, none on [δ, 2b]); a ball meeting a pole evaluates to a non-finite ball, so the analytic flag needs no extra check; non-convergence would show as a wide ball, rejected by rd.rad() > 1e-9 (l.88).
- Rounding: all arithmetic in Arb balls (outward by construction); arb(mid, rad) with arb arguments adds |rad| to mid's own radius (conservative). floor3 and THRESHOLD/CROSSOVER use certain Arb comparisons against exact rationals resp. Arb 4πe^γ, 4πe^{1+γ}.
- Formula (l.86): v = (r1 log A + (n − r1) log B − E)/n = (T)/n with n − r1 = 2r2. Correct.
- Limits (l.57-66): [0, δ] piece ≤ w²δ²/4 from 1 − sech y ≤ y²/2; tails in [0, −log tanh 50] etc. Correct.
- Rerun on a copy (.scratch/odlaudit/rerun, python-flint, 15 MB peak): ALL CERTIFIED, certified.json byte-identical (logs/certify_rerun.log).
- Independent recomputation (consts_mp.py, own mpmath code, tails integrated numerically, not by the closed forms): logs/consts_mp.log, logs/consts_mp_special.log.

## (3) Falsification
- (i) Identity: identity.gp, identity2.gp (logs/identity.log, logs/identity2.log). PARI's lfunzeros on the product ζ_K of ℚ(√−3) misses close pairs of zeros from its two factors (residual 7.5e-9 resp. 1.8e-6 at T = 250); with zeros of ζ and L(s, χ₋₃) taken separately the residual is 5.2e-13. A tool artifact, not a formula error. Control: omitting the r1 term shifts ℚ(√5) by 0.888.
- (ii) Fields: fields_test.py (logs/fields_test.log); class numbers from PARI quadclassunit (classno.gp; GRH-conditional in PARI, test data only). No field violates any certified bound. A known field cannot falsify a correct table, since all table values lie below Odlyzko's GRH bounds and no field is known to violate those; the table sanity rows above are the stronger test.
- (iii) Lean scope: the brief's trap (lesson 136) is avoided; statement fetched from GitHub (raw, 09afbb4) and Mathlib@0493595e DedekindZeta.lean l.49-50 (LSeries definition) read. 31/32 elaborates in ℝ.
- Sharpness control (posit.py, logs/posit.log): with the certified crossover kernel, Re Φ(½ + a + iγ) ≥ 0 for every |a| ≤ 15/32 on a fine γ-grid (minimum −2.8e-15 at a = 15/32, rounding), but a hypothetical zero with β = 0.975 (a = 0.475) at γ ≈ 0.28 would contribute −0.25. The table rests on the Cureton–Thorner strip at exactly its stated width; it would not survive a zero with β > 31/32.

## Wording slips (do not affect any bound)
1. Page Caveat and Odlyzko31/REPORT.md l.84: "for n ≤ 6 the gain over our classical recomputation is at most 0.1%". False for most signatures: ct/classical − 1 = 0.110% (n = 6, r1 = 0), 0.158% (4, 4), 0.245% (5, 5), 0.353% (6, 6). ≤ 0.1% holds only for (n, r1) ∈ {(2,0), (2,2), (3,1), (3,3), (4,0), (4,2), (5,1)}. The record is claimed only from n = 8 on, so the record scope is unaffected (it becomes conservative).
2. "Best classical values known to us": Odlyzko's 1990 survey (§2, after eq. (2.3)) says the best unconditional bounds use L. Tartar's kernel g (refs B10, B11); NOVELTY.md and the Caveat name Diaz y Diaz but not Tartar. At n = 8 the certified ct value 5.661 beats Odlyzko's "5.656 with better kernels" by only 0.08%. The small-n record (8 ≤ n < 2732) is therefore thin and unverified against Tartar/Diaz y Diaz values; n ≥ 2732 (tc) and n ≥ 1618 (tr) are records outright, since every classical admissible kernel stays below 4πe^γ, 4πe^{1+γ} (Prop D at w = ½).
3. Minor: "certified table covers n = 2 … 10⁶" means 103 tabulated degrees, extended to all n by degree monotonicity.

## Residual risks and what this audit did not check
- The proofs of the inputs: Cureton–Thorner (5509 Lean modules) and Bebel (abelian constants only) are taken as audited in Riemann/Deep/HPAudit/REPORT.md, parts (B) and (D). That audit checked the Comparator statement match, the three standard axioms, and acceptance by three kernels, all from third-party CI logs, with no local kernel replay. If either input fell, the bounds would fall back to the classical ones (22.3816, 60.8395), which this audit also recomputed.
- The functional equation, Hadamard theory and the explicit formula are classical and not formalised. They were checked by hand, against Odlyzko 1990 eq. (2.3), and numerically in four fields.
- Poitou 1977 and Diaz y Diaz 1980 were not consulted (no open copy was fetched). Tartar's kernel values were not computed. This affects only the small-n record wording (slip 2), not validity.
- Novelty was not re-audited: Odlyzko's survey (§6, "Open Problem 6.3") discusses how partial knowledge of zeros helps but gives no strip-based constants, which is consistent with NOVELTY.md.

## Review
One reviewer subagent, given Steps 1-5 inline (the Lean statement, the strip derivation via the identity theorem and the functional equation, Lemmas A and B, the positivity convolution including |a| = w, monotonicity, and the Bochner extremality with the limit): CORRECT. Its added notes: at |a| = w the transform 2πδ is a positive measure; take n → ∞ before b → ∞ (or a slow diagonal b(n)).

## Files
consts_mp.py (mpmath constants, Lemma A/B, thresholds, random entries), identity.gp and identity2.gp (explicit formula identity), classno.gp + fields_test.py (falsification fields, table sanity), posit.py (sharpness control), FULLPROOF.md, logs/. grh_spot.py (GRH data-file spot check). Footprints (memcap peaks): consts_mp 174 MB, grh_spot 160 MB, fields_test 88 MB, posit 22 MB, gp scripts ≤ 115 MB stack, certify rerun 15 MB. classno.gp writes .scratch/odlaudit/quadh.csv (deleted after use; rerun classno.gp before fields_test.py).

RH relevance: N12


