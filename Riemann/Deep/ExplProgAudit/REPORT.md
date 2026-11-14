AUDIT PASSED: every step of DERIVATION.md and cert.py checks out, and independent 256-bit code reproduces every claimed constant (1.5, 0.174, 0.0324, 0.0076; 0.72, 0.45331; 7.8978; 10.53); items (1)–(4) hold as written, conditional on H and the cited refereed inputs.

# ExplProgAudit — referee report on discovery linnik-numeric (Riemann/Deep/ExplProg), 14 Nov 2026

Agent explaudit (Orca). Read-only web (arXiv e-prints only); nobody contacted; no Lean; Riemann/Deep/ExplProg/ and discoveries/ untouched. Line-by-line chain with D:/C: line numbers: FULLPROOF.md. Every process ≤ 95 MB (memcap).

## Verdict per item (assume H: no zeros of ζ or any L(s,χ) with Re s > B = 0.8733819666132076)

| Item | Holds as written? | Decisive evidence |
| --- | --- | --- |
| (1) ψ error ≤ 1.5 x^B log²(qx), x ≥ 2; 0.174 (x ≥ 64), 0.0324 (x ≥ e¹⁰), limit 0.0076 | yes | max over pieces 1.48535 (q = 1, x ≈ 2; exact value 1.48478), cells 0.17309 (q = 3, x ≈ 64–70), 0.03229 (X ≥ 10), tail 0.014204 (X ≥ 60), 0.0075898 (X ≥ 10⁶); indep. recomputed |
| (2) p(a,q) ≤ (q log q)^7.8978 (q ≥ 2), (0.72 q log q)^7.8978 (q ≥ 10⁴), (φ(q) log q)^7.8978 (q ≥ 3); c_∞ = 0.453311 | yes | tail margins 21.22% (c₂ = 1, q ≥ 100), 0.7455% (0.72, q ≥ 10⁴), φ-version 21.10%; sieve q < 100 (worst ratio 0.426); indep. recomputed |
| (3) χ(n) ∉ {0,1} for a prime n ≤ (log q)^7.8978, q ≥ 3 | yes | trivial n ≤ q − 1 up to e^25.6 (tightest q = 3: 2 ≤ 2.10), tail c₃ = 1 from log q ≥ 25 (margin 38.87%) |
| (4) \|L′/L(1,χ)\| ≤ 7.8978 log log q + 10.53, χ primitive | yes | worst cell C = 10.5263 at q = 3 (cert 10.526 ± 3·10⁻⁴) |
| HP141 θ₁₄₁ rows (7.8956) | change nothing but B | diff: only l.42 (B) and the core check/print l.47, 49; log has no CERT FAIL |

Exponent: 1/(1 − B) = 7.8977691664599; rounding up to 7.8978 is safe because every base (q log q, c₂ q log q, φ(q) log q, log q, log log q) exceeds 1 (resp. 0) in its range.

## What was checked
- **A. Derivation (all of DERIVATION.md).** Explicit formula for compactly supported weights (residues xW(1), −x^ρW(ρ), −x^{−m}W(−m) with T(χ) = {2,4,…}, {0,2,…}, {1,3,…}; the even-character residue −W(0); (1.1)); imprimitive reduction (1.2); orthogonality (1.3) with Ω over conductors 3 ≤ q* | q; pairing ρ ↔ 1 − ρ̄ incl. β = 1/2 and multiplicities; Lemma Z (BMOR → N⁺ via conductor monotonicity, concavity, the ℓ ≤ 1.567 branch via R(t;9) ≥ 0; Stieltjes with real zeros); all antiderivatives by differentiation; Lemma Zζ, Zζ'; kernels K1, K2 (all constants); Theorems 2, 2', 3 (tail lemma monotonicity; f(q) increasing for q ≥ 16); Theorem 1 (6.1–6.4, coverage); Theorem 4 (EF1, Bernoulli, Dusart). No false step.
- **B. Sources (TeX).** BMOR 2005.02989 (only version = Math. Comp. 2021) Thm 1.1 l.170–175 and N(T,χ) l.157–165: verbatim as used. Dusart 1002.0442 (Σ_{p≤x} log p/p ≤ log x + E + 0.2/log x + 0.2/log²x, x ≥ 2974): as used. Bennett et al. 1802.00085 v3: the page's exp(8√q log³q) is their all-q corollary (0.03 for q ≥ 58). H: HPAudit PASSED (every modulus and χ, ζ as q = 1).
- **Low zeros at height 0 (brief question).** Legitimate: BMOR is used only at T ≥ 5/7; below, only N(t) ≤ N(5/7), and those zeros get weight G(0) = n₀ = max G.
- **C. Code (all 679 lines of cert.py).** Every formula = its derivation line; all inputs exact rationals (rg list in FULLPROOF.md); balls at 200 bits, every proof inequality a certified pos(); floats only choose δ (then exact rationals) and cell edges; undecided comparisons abort (one immaterial exception, NOTE 2). Cells: X-slabs 4.16 → ≥ 60 contiguous, L-cells log 3 → ≥ X₁ contiguous, q = 1, 2 separate; x ∈ [64, e^4.16) inside the exact range (x < 66); tail X ≥ 60 for 0 ≤ L ≤ X; x < q separately. No gap. Theorem 3 uses Ω₁/L (no 1/q), correctly.
- **D. Independent recomputation (indep.py, mpmath 256 bits; zero sums by quadrature, not cert's closed forms).** All §0–§3 constants, all Theorem 2/2'/3 table margins (to every printed digit), 21 Theorem-1 cells at 256 bits and 133 at 64 bits (all extremes, X ≈ 10, 20, 60 boundaries, 120 random) (indep ≤ cert, |rel. diff| ≤ 9.2·10⁻⁴), the 6.4 tail (0.0142041, 0.011127, 0.0078753, 0.0075898) and direct non-majorised T₁…T₆ at X ≥ 60 (max 0.008495), Theorem 4 cells, the exact range. kcheck.py: kernel majorants K1, K2, Thm-4 on 4000 random zeros (max ratio 0.992 ≤ 1).
- **E. Falsification.** q = 1, x = 2 (1.48478, forced); q = 2, 3, 4 exact ψ to 10⁶ (max R 0.568, 0.119, 0.078); least primes q ≤ 300 and the sieve/analytic boundary q = 95–105 (ratios ≤ 0.16); n_χ for q ≤ 400 (tightest q = 3); |L′/L(1,χ)| for q = 3, 4, 5, 7, 8 (≤ 0.64); sanity.py at reduced size (N = 10⁷, q ≤ 1500: failures 0). Real characters: nothing but H excludes an exceptional zero; one real zero with β > B would break Lemma P.
- **F. Reviewer** (Lemma P + Lemma Z, BMOR branches, low zeros, real zeros, multiplicities): CORRECT.

## Notes (none affects a claim)
1. D:158 justifies freezing k₁, k₂ at Y_m by monotonicity; the simpler valid reason is G ≤ G_m pointwise, so Lemma Z applies to G_m directly.
2. zsum_coeffs/zsum_zeta: an undecided ball comparison lo > 5/7 (or > 14) integrates from 5/7; error O(radius²) ≈ 10⁻¹¹⁰.
3. 0.174, 0.0324, 0.0147, 0.0112, 0.0079, 0.0076 and 10.53 are read off printed ball upper bounds, not pos()-checked against the decimal; margins 9·10⁻⁴, 1.1·10⁻⁴, 1.0·10⁻⁴, 6.7·10⁻⁵, 2.5·10⁻⁵, 9.8·10⁻⁶, 3.7·10⁻³ ≫ float rounding.
4. Page "limit constant 0.45331" is c_∞ = 0.4533109… rounded at the 6th digit; a description of the limit, not a bound.
5. D:184 "the tail lemma of §4 applies verbatim": not verbatim (Ω₁ has no factor 1/q in Theorem 3); cert.py uses the right form.
6. D:214–224 rounds constants up (0.1137 ≥ 0.113685, 0.0196 ≥ 0.274/14, 0.00979 ≥ 0.137/14); cert.py uses the exact values.
7. Page 'Before': exp(8√q log³q) is Bennett et al.'s all-q corollary; they also give exp(0.03√q log³q) for q ≥ 58. The comparison stands.

## Not checked (scope)
- H itself (HPAudit PASSED; Lean-checked, unrefereed: the page's caveat stands). Proofs of BMOR, Rosser 1941 Thm 19, RS62 Thms 9/12/15, Dusart, Davenport's Σ Re(1/ρ) (refereed; statements checked: BMOR and Dusart in TeX, the others as standard citations). Arb/python-flint correctness (trusted). Independent recomputation covers 133 of the 7307 cells (all extremes and boundaries), not every cell.

## Files
FULLPROOF.md (chain, D:/C: line numbers, per-step OK/NOTE); indep.py (modes const lemz linnik nonres pnt ptail lprime small); kcheck.py; falsify.py; certcells.py (cert's own per-cell values); logs/ (indep_*.txt, kcheck_output.txt, falsify_output.txt, sanity_small_output.txt). Commands: cd Riemann/Deep/ExplProgAudit && uv run --no-project --with mpmath python3 indep.py <mode> (each ≤ 4 min, ≤ 100 MB); uv run --no-project --with numpy --with mpmath python3 falsify.py.

RH relevance: N10, N11
