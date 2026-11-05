VERDICT: CORRECT
NOVELTY: NEW

# Independent review of PROOF.md (diagonal Fasi–Higham positivity, all odd m)

Reviewer: fhrev, 2026-11-05. I read PROOF.md, REPORT.md and all the code (fh2lib, good, goodsmall, mktab, assembly, large_tau,
smalln, identities, verify, the sanity scripts, ../FasiHigham/arb_cert.py). I re-derived every lemma by hand, re-ran
verify.py, and tested the theorem and the intermediate bounds with independent code. I found no gap. All findings
below are cosmetic or concern the sanity tooling; none affects the proof. Scripts and logs are in `.scratch/fhrev/`.

## Verdict in one paragraph

The exact layer is right. That covers the Padé remainder with E_j > 0, p = θ_m(x)/θ_m(0), X = πI_ν/K_ν (which needs
e^{iπν} = −i for odd m), the normalisations with κ, the left-half-plane formula and the reflection identity, and the
contour decomposition of Lemma 2. The analytic layer is right too: the Liouville–Green bounds of Lemma 1, the
saddle-point lower bound of Lemmas 3/3′, and the crude, Neg and I_L upper bounds. The finite checks are done on
intervals in ball arithmetic. In τ they use t-ranges, concavity of Φ, and linearity in t of the crude rates. In r they
use r*-boxes and covering cells. In angle they use chunked ball suprema. Every "valid for all ν ≥ 203.5" claim rests on
an explicit monotonicity that the code checks. The cases n ≤ 12, {n ≥ 8, n < 60ν} and {n ≥ 60ν} cover every n. The
cases m = 1, odd 3 ≤ m ≤ 201 and odd m ≥ 203 cover every odd m. verify.py passes. Independent exact and ball
computations agree with the theorem and with every intermediate bound I tested.

## Findings

1. **§1 (R), (I1), (I2), (LG): re-derived, correct.** (Severity: none.)
   - (R): z^{2m+1}∫e^{tz}Φ = (−1)^m(e^z q − p), with E_j = m!(m+j)!/(j!(2m)!(2m+j+1)!).
   - (I1): K_{m+½}(x)x^ν = √(π/2)e^{−x}θ_m(x) (DLMF 10.49.12), and θ_m(z/2)/θ_m(0) = p(z).
   - (I2): K_ν(xe^{iπ}) = e^{−iπν}K_ν − iπI_ν, and −iπe^{iπν} = (−1)^mπ = −π for odd m.
   - (LG): e^z/p = Λe^{ν(ζ+Ψ_e)}(1+ζ²)^{1/4}/(1+h_K) and X = κe^{2νη}(1+h_I)/(1+h_K). Here
     (2πν)^{−1/2}(2ν/π)^{1/2} = 1/π, and κ makes h_I(0) = 0.
   - Useful closed form: 1 + h_K(ζ) = √ν ν^{−ν} θ_m(νζ) e^{ν(Ψ_e−ζ)} (1+ζ²)^{1/4}. So 1 + h_K continues to all of
     ℂ \ ±i[1,∞). The reflection identity therefore holds for every Re w > 0, not only |w| < 1. h_I is even
     because I_ν(x)x^{−ν}, Ψ_e and (1+ζ²) are even.

2. **§2 Lemma 1 (LG error bounds): correct. This is the complex-variable LG theorem.** (Severity: cosmetic.)
   - The Volterra equation h = ∫K ψ(1+h) has kernel K = (1−e^{−2ν(v−ξ)})/(2ν). On ξ-progressive paths |K| ≤ 1/ν.
   - Picard iteration gives |h| ≤ exp(∫|ψ||dξ|/ν) − 1. By hand: ψ dξ/dζ = ζp³(15p²−3)/12 = 2U1′(ζ), so the bound is
     exp(2V(U1)/ν) − 1.
   - This is Olver, *Asymptotics and Special Functions*, Ch. 6 §11 (Thm 11.1), with error-control function
     F = −∫ψ dξ and V(F) = 2V(U1). Its hypotheses hold:
     - f = (1+ζ²)/ζ² and g = −1/(4ζ²) are holomorphic in the simply connected domain {Re ζ > 0}.
     - f ≠ 0 there; the turning points ±i lie on the boundary.
     - ξ = η is single-valued.
     - The reference points are ζ → ∞ for K and ζ = 0 (ξ → −∞) for I.
     - Every point used is reached by a progressive path with finite variation:
       rays (d Re η/dρ = Re s/ρ > 0); rightward horizontals (arg s ∈ (0, arg ζ) in the first quadrant, also for
       |ζ| > 1); and [0, iy] for y < 1.
   - Identification with K_ν and I_ν follows from recessiveness at the reference point. Points on Re ζ = 0 follow by
     continuity.
   - Uniformity in ν: V does not depend on ν, so ε = e^{2V/ν} − 1 decreases in ν.
   - Nit: §2 says the factor 2 is "lost against Olver's sharper constant". The bound used is exactly the
     complex-variable theorem. The observed ratio |h|/bound ≈ ½ reflects the leading term |h| ≈ V/ν (the
     real-variable constant).

3. **§2 λ-grid: validated against true values.** (Severity: none.) Test: `lgtest.py` and `lgextra.py`.
   - h_K is computed from the exact θ_m form above; h_I from Arb's `bessel_i`.
   - All 2448 grid cells, 9 points each (22 032 points), at ν0 = 203.5: max |h_I|/ε_I = 0.494 and
     max |h_K|/ε_K = 0.430.
   - The disc |ζ| ≤ 0.3, including Re ζ < 0: ≤ 0.49.
   - 12 points beyond the grid, up to |ζ| = 1.85 (these are used directly by good.py's Cauchy regions): ≤ 0.49.
   - lg_box returns no bound near ζ = i, so a missing cell makes the run fail rather than pass silently.

4. **§4 Lemmas 3/3′: re-derived step by step, correct.** (Severity: cosmetic, plus one condition worth stating.)
   - Checked: the critical point r_ex (F′ ≥ (B − λc₂/ν)/ρ, with B increasing); the Gaussian constants M4/(8β²),
     (5/24)M3²β^{1/2}β₋^{−7/2} and the Mills tail; the purely imaginary cubic term; the (F1) deficit with the
     decreasing-in-φ chunk evaluation; the horizontal-segment monotonicity in t and in ρ; the ν- and a-monotonicity
     conditions; the n-monotonicity of every term in Lemma 3′ (quotients of affine functions); and ν_eff ≤ ν.
   - Cosmetic: multiplying E by √(b_h/b_l) is unnecessary, since (2πβ)^{−1/2} ≥ (2πβ_hi)^{−1/2} already. It is harmless.
   - Unstated condition: assembly reads λ_W0 from the grid on [0.8r_lo, min(1.25r_hi, 1.74)]. This needs
     |r_ex − r*| ≤ 0.2r* and r_ex ≤ 1.74. Mechanically checked (`mech.py`): max Δ/r_lo = 0.103, and
     max r_hi + Δ = 1.6015. Recomputed e values equal the table entries. Fix: state the condition in §4.

5. **§5–§6 assembly: the l-partition is exhaustive, and the bad-term bounds are correct.** (Severity: none.)
   - The loop covers l ≤ ⌊τ_b/2⌋ (N_l ≥ 0). Neg covers l ≥ ⌊τ_a/2⌋+1 (N_l ≤ −1, since boxes are cut at even
     integers). The extra crude term at τ_b = 2L only adds to the bad sum.
   - Arc, segment and horizontal crude pieces: re-derived. y^{−N} ≤ r_c^{−N}, ∫dy/y gives the log log term, and
     dζ/ζ gives the factor 1/|ζ| on the horizontals.
   - Neg geometric series: re-derived. Both the "big" and "small" kap cases are valid for all ν ≥ ν0.
   - I_L: re-derived.
   - ν-monotonicity: the code requires ν0·gap > 1, which is stronger than the ½ stated in §6.

6. **§7 Proposition 2: re-derived, correct.** (Severity: none.)
   - Re-derived: dρ/dl = −2s²/(s+aρ); the identity Φ = r(1−log r) + (1−r)g + τg; I0 = 1.0464 and g_min = 0.40154
     (spot values 1.0475 and 0.4019 at r = ½).
   - The classes far / good (l ≥ 5) / start / Neg / I_L are exhaustive, and l̂ ≥ 25.
   - The counts and the τ- and ν-derivative bounds are correct.
   - Consistency re-runs of large_tau.py with τ0 = 100, 200 and 1000 all pass.

7. **§8 smalln.py: the positivity criterion is slightly too weak.** (Severity: minor, code only; no gap.)
   - It accepts coefficients that are ≥ 0 in x = m − 203, so a zero constant term would allow c_n(203) = 0.
   - All 12 numerator and denominator constant terms are in fact > 0 (`smalln_const.py`), and c_n(203) > 0 exactly.
   - Fix: require the constant terms to be > 0.

8. **Coverage and off-by-one: correct. REPORT.md describes the case split wrongly.** (Severity: presentation.)
   - The ranges are n ≤ 12 (Prop 3) ∪ {n ≥ 8, τ < 60} (Prop 1) ∪ {τ ≥ 60} (Prop 2). They overlap at 8 ≤ n ≤ 12;
     the boxes are half-open, and Prop 2 starts exactly at τ = 60.
   - In m: m = 1; odd 3 ≤ m ≤ 201 (verify.py `arb_ok` demands exactly range(3, 202, 2), and my re-run gives 100/100
     POS); odd m ≥ 203 ⟺ ν ≥ 203.5. N_l = n − l(2m+1) is an integer.
   - REPORT.md says the proof covers "n ≤ 2m through 2m ≤ n ≤ Cm² up to n ≥ Cm²". That is not the proof's split. The
     analytic proof has no Cm² threshold: Prop 2 covers all n ≥ 60ν. The ≈0.73m² threshold belongs only to the
     partial-fraction certificate for m ≤ 201.
   - Fix: reword REPORT.md as "n ≤ 12; 8 ≤ n < 60ν; n ≥ 60ν".

9. **sanity.py (S1): a comment has the direction backwards.** (Severity: cosmetic.)
   - The comment says "gf at nu0 is smaller". It is larger: at ν > ν0 the true G is G* − ½log(ν/ν0) + l/(12ν0) −
     l/(12ν).
   - S1 runs at m = 203, where ν = ν0, so it is unaffected. My checks below use the correct ν-dependence.

10. **check_neg_IL.py (S3) is only valid for m up to about 250.** (Severity: sanity tooling only; not used by the proof.)
    - It truncates the E-series at 1100 terms. At m = 401 (|z| up to 683) it reports false "violations", for example
      n = 20: Neg numeric 2136 against bound 32.6.
    - With X = πI_ν/K_ν taken from Arb Bessel functions (`negtest.py`), the m = 203 values reproduce S3
      (−5.77, 37.44, 606.32), and at m = 401 every Neg bound holds. The six n tested go up to 30 000 (τ = 74.7),
      with log-margins of 25 to 127.
    - Fix: use the Bessel form, or EJ ≳ 3|z|max.

11. **verify.py: the step C success test is weak.** (Severity: cosmetic.) It only counts " boxes " lines. Validity is
    still enforced downstream: assembly.good_factor and large_tau re-check coverage and e < 0.9. Fix: assert those
    conditions in step C as well.

12. **Robustness.** The proof has comfortable margins. The worst bad/good ratio is 5.5e−18, and the minimum rate gap is
    0.183 (ν0·gap ≥ 37). Every Good-bound check below sits above the bound by a log-margin ≥ 0.044.
    Residual risk is the usual one for a computer-assisted proof: correctness of python-flint/Arb and of about
    1300 lines of Python. I reviewed that code; it is not formally verified.

## Independent checks (all passed)

| check | what | result |
|---|---|---|
| verify.py --sanity (run on a copy in .scratch/fhrev/run) | steps A–G and S1–S4 | ALL STEPS PASSED (1298 s); regenerated tab203.json is byte-identical to the committed file |
| exact rationals (python-flint fmpq_series, `exact.py exact`) | m = 203 (n ≤ 6000); m = 205, 251, 301, 351, 401 (n ≤ 5000) | every c_n > 0. The proof's lower bound (best Good l, correct ν-dependence) is ≤ the true value for every n ≥ 8; min log(true/LB) = 0.0437–0.0461, smallest at l = 0, r* ≈ 0.2, then near τ ≈ 2.0–2.5 (0.050). The margin profile does not depend on ν. |
| Arb balls (`exact.py arb`) | m = 203, n ≤ 40 000 (τ ≤ 98); m = 401, n ≤ 30 000 (τ ≤ 75), across the τ = 60 seam | every c_n > 0; LB ≤ true for sampled n, min margin 0.045 / 0.052 |
| Lemma 3 by quadrature, independent integrand (θ_m and Arb I_ν, no E-series; `goodtest.py`) | 61 cases at m = 203, 67 at m = 401, 10 at m = 1001, 9 at m = 203 with l = 100, 300, 1000. Includes the weakest table rows: l ≥ 5 at r* ≈ 0.1165 (e = 0.49), l = 1 at r* ≈ 0.031, rows straddling y_J, and r* ≈ 1.59 | every J_l > 0 and above the bound. Min margin 0.045; at the weak rows 0.40–2.4. Reproduces S2 (log J = 615.3595). |
| crude bounds (`crudetest.py`) | m = 401: start terms l = 1–5 with N = 0–40, and far terms | all hold; min log-margin 0.97 |
| Neg bounds (`negtest.py`) | m = 203 and m = 401 | all hold (see finding 10) |
| LG grid (`lgtest.py`, `lgextra.py`) | see finding 3 | all hold, ratio ≤ 0.494 |
| table sanity (`mech.py`) | r_ex range against the grid edge; e recomputed | see finding 4 |

## Novelty (2026-11-05)

No earlier proof of the diagonal case, for all odd m or for all m ≥ any bound, was found. The sources below are listed
newest first. Downloads are in `.scratch/fhrev/nov/`.

- **Goodchild, Zenodo 23919117.**
  - The Zenodo API shows a single version, v1 of 3 Nov 2026.
  - The GitHub repo JacobGoodchild/findformula3 has 81 commits; the last push was 3 Nov 13:02 UTC.
  - Its results/pade_positivity/README says: "(1) is a theorem for each listed m, not for all m (the general diagonal
    statement is still open)".
  - The paper's sec_pade.tex lists "Is 1/q_{m,m} positive for every odd m?" as an open problem.
  - It certifies odd m ≤ 361 individually. A numerical scan reaches m ≤ 601 with n ≤ 20000.
  - Goodchild's newer repos (findformula4, formwork) contain no Padé or Fasi files.
- **Scardino, arXiv 2610.39051.** Only v1 exists (28 Oct 2026). It proves k = 0 and eventual positivity, and says the
  diagonal is open. The evidence API retrieves exactly these passages: "The diagonal case is precisely the case that our
  counterexamples leave untouched … remains an open problem".
- **arXiv monthly listings.** math.CA, math.NA and math.CV for 2026-10 and 2026-11, 1530 entries in total. The only
  relevant hit is 2610.39051. Unrelated: 2610.25457 (Sheffer, dual Padé), 2610.09251 (repeated-pole approximation of
  e^z), 2610.39453. arXiv API searches (Fasi AND Higham; Padé AND denominator AND coefficients; "reverse Bessel
  polynomial"; Padé AND "matrix exponential") find nothing new.
- **Citations of Fasi–Higham 2019** (doi 10.1137/18M1228876). OpenAlex lists 30 citing works, Semantic Scholar 32 with
  contexts, Crossref 26. The only one on positivity is Scardino.
- **Fasi, PhD thesis (Manchester 2019), Conj. 5.2.** No proof: "none of them led us to a proof of Conjecture 5.2".
- **MathWorld** "Fasi-Higham Positivity Conjecture": "The diagonal case remains open".
- **Other searches.** MathOverflow and Math.SE (StackExchange API): nothing relevant; MO 41226 is about the numerator.
  Zenodo searches (Fasi AND Padé; Padé AND exponential; "reverse Bessel"): only Goodchild. Evidence API, 2 claims: only
  Scardino. MathDB was behind a Cloudflare wall (on 4 Nov fhnov saw it marked "Open").

Hence NEW. Positivity for m ≤ 361 was already certified, one m at a time, by Goodchild (and ours for m ≤ 201). The new
content is the uniform proof for all odd m ≥ 203, and with it the first proof for every odd m ≥ 363.
