RECORD: Lambda > 3.2648 unconditionally for consecutive zeros on the critical line (was 3.18, Bui-Milinovich)

# gaplambda: large gaps between consecutive critical zeros (node N8)

## Claim
Let Λ = limsup (t_{n+1} − t_n)/(2π/log t_n), t_n the ordinates of zeros of ζ on the critical line. Then
**Λ > 3.2648**, unconditionally. Bui and Milinovich's criterion (arXiv 1410.3635v3: Wirtinger's inequality plus
the Bettin-Bui-Li-Radziwiłł twisted fourth moment, ϑ < 1/4) gives Λ > κ whenever h1 = c0/(κ²c1) > 1. We take
r = 1, ϑ = 2499/10000, v = 131063/100000 and an explicit degree-14 amplifier polynomial P, and certify
h1 = 1.00005692466 ± 2.5·10⁻¹² at κ = 4081/1250 = 3.2648 in Arb (`cert_bm.py`). Under RH the bound also holds
for all zeros, but there the record is 3.638 (Durkan-Pearce-Crump 2026).

## Before: which bounds need RH
| Bound | Who | Hypothesis | Source (TeX line) |
|---|---|---|---|
| Λ ≥ √(11/2) = 2.345 | Hall 2002, L⁴ Wirtinger | none (critical zeros) | 2511.14415 l.116, l.477 |
| Λ > 2.63 / 2.76 | Hall 2005 / Bredberg 2011 | none | 2511.14415 l.116 |
| **Λ > 3.18** | Bui-Milinovich (1410.3635, QJM 2018) | **none** for critical zeros; RH only for all zeros (λ) | 1410.3635v3 l.53, l.70; 2610.28527 l.81 |
| Λ > 2.64 | Bui-Hall-Subirà Jorge (2511.14415) | none | l.118-122 |
| λ > 3.18 | Bui-Milinovich | RH | 2510.14309 l.423 |
| **λ > 3.63846723** | Durkan-Pearce-Crump (2610.28527, Oct 2026), all primitive Dirichlet L incl. ζ | **RH** | 2610.28527 l.45, l.76 |

## Method
1. **BM's reduction, re-derived** (DERIVATION_BM.md §2): F = e^{ivt𝓛}ζ(½+it)ζ(½+it+iκπ/𝓛)M(½+it). Under Λ ≤ κ
   its zeros are at most (1+o(1))κπ/𝓛 apart. Wirtinger and BM's mean values (Thm smoothed, l.300-313) then give
   Λ > κ if h1 = c0/(κ²c1) > 1 (l.327-330).
2. **Amplifier.** BM used a degree-2 P (l.333). c0 and c1 are quadratic forms in P, so the best P of each degree
   solves a generalised eigenproblem (float search, `search/bm/`). Largest κ by degree (float scan, then
   certified at the marked degrees):

| deg P | 2 | 4 | 6 | 8 | 10 | 11 | 12 | 14 |
|---|---|---|---|---|---|---|---|---|
| κ_max ≈ | 3.183 | 3.228 | 3.246 | 3.2548 | 3.2597 | 3.2614 | 3.2628 | 3.2648 |
| certified at κ = | 3.18 | | | 3.254 | 3.258 | 3.261 | 3.2628 | **3.2648** |

3. **Certificate** (DERIVATION_BM.md §3):
   - t3 and t4 in closed form; the polynomial factor cancels the E-function denominators, leaving four terms.
   - The x-integrals are exact (Taylor-expanded phases with rigorous remainder, Dirichlet integrals).
   - (t1, t2) use a 20-point Gauss-Legendre rule with a rigorous Bernstein-ellipse remainder.
   - Everything runs in Arb balls (192 bits) from exact rationals.
   - Float quadrature is unreliable at high degree: a degree-12 P from a 14-point float scan certifies only to
     h1 = 0.9975 (float said 1.0024). Every reported value is from the certificate.

## Certificate
`cd Riemann/Deep/GapLambda && uv run --no-project --with python-flint python3 cert_bm.py` (≈ 2.5 min, ≈ 32 MB;
output `logs/cert_bm.txt`):
- (a) BM's published numerics (ϑ = 1/4, κ = 3.18, v = 1.26, P = 1 − 5.8x + 6.4x²): h1 = 1.0021778 > 1.0002 (l.337).
- (b) The same P at ϑ = 2499/10000: h1 = 1.00201577199 > 1, a rigorous Λ > 3.18 with ϑ strictly inside BBLR's
  range (the first certificate of their theorem).
- (c) Record: deg-14 P, κ = 3.2648, h1 ∈ [1.00005692466 ± 2.5·10⁻¹²], H = c0 − κ²c1 ∈ [2.4580659·10⁻⁹ ± 4.5·10⁻¹⁷] > 0.
- ALL CHECKS PASS.

## Checks
- Ball values equal the float quadrature where it converges (degree 2: c0 = 1.29667198837·10⁻⁵ vs 1.296672·10⁻⁵;
  degree 8: h1 = 1.004200 at κ = 3.25 in both).
- A Monte Carlo of BM's raw printed 9-fold integrand (l.307-311) matches the float quadrature (z = 0.61, 0.27).
- Im c0 and Im c1 contain 0 (the main terms are real), and c0, c1 > 0.
- The remainders are tiny: Taylor ≤ 6·10⁻²², t-Gauss ≤ 10⁻²⁴ before the factor (κπ)⁻²; ‖P‖_∞ is bounded rigorously
  (Taylor bound on 400 subintervals).
- The search used 14-, 16- and 18-point float scans: 14 points fail at degree 12 (certified 0.9975), while 16 and
  18 points match the certificate (degree 12: float 0.998547 at κ = 3.2645 vs the certified slope).

## Novelty
- No unconditional improvement of 3.18 exists in the sources searched: BM (doi 10.1093/qmath/hax047) has 22 citing
  works in OpenAlex, none improving it.
- BHSJ (Nov 2025, l.116) and Durkan-Pearce-Crump (Oct 2026, l.81) still call 3.18 the best Wirtinger value.
- The Oct-Nov 2026 arXiv listings, Zenodo and HAL show nothing; details in NOVELTY.md.
- The method is BM's own; the gain comes from a higher-degree polynomial plus a certificate (BM gave numerics only).

## Caveat
- Taken as printed: BM's Lemma main and Theorem smoothed (main terms c_j, l.300-313, l.370-384) and BBLR's twisted
  fourth moment (l.500-517). We re-derived the reduction, the c_j at the shifts and every evaluation step
  (DERIVATION_BM.md), not those mean-value theorems.
- The certificate proves the main-term inequality h1 > 1 rigorously; the passage to Λ > κ is BM's argument
  (o(1) terms as T → ∞).
- Not yet reviewed by a referee-style reviewer; the orchestrator should run one (DERIVATION_BM.md §3 is the core).
- The ceiling of the method is near 3.266-3.27 (degree steps 0.0017, 0.0013, then about 0.001 per degree).
  Larger r showed no gain (r = 2, degree 8: h1(3.25) = 0.9986).
- Repo flag (kept): Riemann/Tree/R421.lean l.13 mislabels 3.18 as RH-conditional.

RH relevance: N8 (zero spacing; Montgomery's Λ = ∞ direction); low to moderate (record for a classical
unconditional constant).

## Files
- **Certificate:** cert_bm.py (BM criterion) and logs/cert_bm.txt.
- **Derivations:** DERIVATION_BM.md (reduction, re-derived vs printed, why any P is admissible).
- **Float search:** search/bm/ (bmfloat.py, bmscan.py, bmopt.py, bmcheck.py, bmmc.py, bmtest.py) and logs/bm/
  (scans with JSON, cert_tries/ with every certificate attempt, including the failed 14-point degree-12 one).
- **L⁴ route:** cert.py, core.py, rcheck.py, DERIVATION.md, search/, logs/ (section below).
- NOVELTY.md

## L⁴ route (Bui-Hall-Subirà Jorge 2511.14415): certified Λ > 2.907, not a record
- Optimising P and r in BHSJ's L⁴-Wirtinger argument (Thm Wi l.507-516, amplified fourth moment Thm 1.1,
  ϑ < 1/8) certifies Λ > 2.907 with a degree-11 P (r = 1, ϑ = 1249/10000); `cert.py`, DERIVATION.md.
- Per degree, r = 1: 2.6410 (P = 1), 2.7435, 2.7993, 2.8328, 2.8548, 2.8694, 2.8800, 2.8883, 2.8947, 2.8998,
  2.9038, 2.9071 (degree 11). It levels off near 2.94; r = 2, 3 are worse (2.6168, 2.5749 at P = 1).
- Their printed u = 0.6, v = 0.4 (l.684) gives 2.6395, not > 2.64. The theorem holds with u = 0.5942, v = 0.4038
  (certified 2.641).
- Checks: Hall's √(11/2) is reproduced exactly; a Monte Carlo of the printed 16-fold integrand matches.
- These are quartic forms in P, so there is no eigenproblem; the search is multi-start BFGS on well-conditioned
  tensors.

## History
- Brief 1: the L⁴ route (above).
- Follow-up 1: BM's input is BBLR's twisted fourth moment, ϑ < 1/4 (l.246, l.259, l.301, l.371, l.500-517).
  No longer input applies: Conrey's 4/7 is a second moment, Alpöge-Furman and Lamzouri are pair-correlation
  arguments, and BHSJ Thm 1.1 is BBLR applied to A². The float estimate Λ > 3.25 came from optimising P.
- Follow-up 2: the certificate above.
