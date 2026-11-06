RECORD: unconditional root-discriminant bounds beyond the 1977 Odlyzko–Poitou limits. Every number field K of degree n with r1 real places has rd_K = |d_K|^{1/n} ≥ (65.6439 − o(1))^{r1/n} (23.3878 − o(1))^{2r2/n} (was 4πe^{1+γ} = 60.8395 and 4πe^γ = 22.3816), from Cureton–Thorner's Lean-checked strip Re ρ ≤ 31/32. Explicitly, rd_K > 22.38169 for every K of degree ≥ 2732 and rd_K > 60.84021 for every totally real K of degree ≥ 1618. The classical method stays below 4πe^γ (resp. 4πe^{1+γ}) at every degree. The values are certified in Arb, with a table for n = 2 … 10^6.

# odlyzko31 report (N12, N7; 6 Nov 2026)

## Claim
1. **Asymptotic constants.** Under H(w), meaning every zero of ζ_K satisfies |Re ρ − 1/2| ≤ w, we get rd_K ≥ (A_r(w) − o(1))^{r1/n} (A_c(w) − o(1))^{2r2/n}, where A_c(w) = exp(γ + log 8π − I_c(w)), A_r(w) = A_c(w)·exp(π/2 − I_r(w)), and I_c(w), I_r(w) are the integrals of (1 − sech wx)/(2 sinh(x/2)) and (1 − sech wx)/(2 cosh(x/2)) over (0, ∞). Within the explicit-formula method these constants are the best possible: every admissible kernel satisfies F ≤ sech(wx) (DERIVATION §3). The certified values (Arb enclosures) are:

| zero information | w | totally real A_r | all signatures A_c |
|---|---|---|---|
| none (classical: Odlyzko, Poitou) | 1/2 | 60.839540 | 22.381616 |
| **Cureton–Thorner 31/32, every K (Lean-checked)** | 15/32 | **65.643966** | **23.387894** |
| Bebel's Dirichlet half-plane, **abelian K** (Lean-checked); every K given Liu's unrefereed B₂ | 0.3733819666132076 | 84.299917 | 26.966253 |
| Liu arXiv 2610.26835 B♭ (unrefereed), every K | 0.374886647040 | 83.949939 | 26.903236 |
| GRH | 0 | 215.332516 | 44.763232 |

2. **Crossover degrees, all certified.** rd_K > 22.38169 > 4πe^γ for every K of degree n ≥ 2732, and rd_K > 60.84021 > 4πe^{1+γ} for every totally real K of degree n ≥ 1618. Equivalently, a field with rd_K ≤ 4πe^γ has degree < 2732; classically no degree bound existed for such fields. Smaller thresholds come earlier: rd_K > 22 for every K of degree n ≥ 1581, against n = 8862, the first degree where Odlyzko's 1976 Table 2 exceeds 22. Totally real rd_K > 60 for n ≥ 1223. Both are certified (THRESHOLD lines).
3. **Table** ([TABLES.md](TABLES.md), 693 table entries). Selected entries; each holds for degree n and every larger degree. The "all" rows use formal r1 = 0, which is valid for every signature:

| n | Odlyzko 1976 uncond. (all / tot. real) | classical, ours | **31/32, ours** | abelian (B), ours | Odlyzko 1976 GRH |
|---|---|---|---|---|---|
| 8 | 5.646 / 10.568 | 5.651 / 10.573 | **5.661 / 10.630** | 5.689 / 10.787 | 5.734 / 11.036 |
| 20 | 9.779 / 21.642 | 9.783 / 21.653 | **9.844 / 21.986** | 10.004 / 22.909 | 10.270 / 24.594 |
| 100 | 16.454 / 41.728 | 16.459 / 41.741 | **16.760 / 43.304** | 17.623 / 48.128 | 19.268 / 59.141 |
| 1000 | 20.895 / 55.966 | 20.898 / 55.974 | **21.571 / 59.307** | 23.704 / 70.760 | 29.094 / 107.548 |
| 10^4 | 22.049 / 59.746 | 22.050 / 59.750 | **22.913 / 63.934** | 25.784 / 79.081 | 34.768 / 141.218 |
| 10^6 | 22.348 / 60.691 | 22.366 / 60.788 | **23.340 / 65.457** | 26.708 / 83.041 | 39.923 / 176.415 |
| ∞ | 22.3816 / 60.8395 | same | **23.3878 / 65.6439** | 26.9662 / 84.2999 | 44.7632 / 215.3325 |

   The table also gives mixed signatures for every r1 when n ≤ 12, e.g. n = 8: 5.661, 6.560, 7.656, 8.994, 10.630 for r1 = 0, 2, 4, 6, 8. For any (n, r1), [tables.csv](tables.csv) lists (b, log A, log B, E) per kernel, and rd_K ≥ exp((r1 log A + 2r2 log B − E)/n).
4. **Finiteness.** Only finitely many number fields have rd_K < 23.3878 (classically < 22.3816), and only finitely many totally real fields have rd_K < 65.6439 (classically < 60.8395).
5. **Abelian fields (ii).** Bebel's half-plane gives the "abelian" column: limits 26.9662 and 84.2999, and 22.3825 for n ≥ 482 (all fields) and 60.8673 for totally real n ≥ 285. Caveat: for abelian K the conductor–discriminant formula alone gives rd_K ≥ (product of the n smallest primitive conductors)^{1/n}, which grows like √n. Among the tabulated degrees it is stronger than the B-bound from n = 220 (all fields) and n = 1332 (totally real); see column elem in TABLES.md. So the abelian column is the best bound we know only for moderate n.
6. **Ghitza–Yamauchi (iii): no.** In arXiv 2509.00635 (TeX l. 397, 489, 521, 567), GRH enters only through Odlyzko's Table 1. They need:
   - rd ≥ 27.328 for n ≥ 660, any signature;
   - rd ≥ 32 for some degree bound (n < 4800), any signature;
   - rd ≥ 81 (n < 280) and rd ≥ 128 (n < 4800) for totally real fields.

   Their fields K = Q̄^{ker ρ̄} are Galois with nonabelian group, so only 31/32 applies unconditionally, and its suprema 23.3879 and 65.6440 lie below every threshold, at every degree. None of Theorem A(a)–(c), Cor. B, Cor. E or Cor. F(2) becomes unconditional.
   - Under Liu's unrefereed B₂ the totally real limit 84.2999 exceeds 81: rd > 81 for totally real n ≥ 35,448 (certified THRESHOLD line). That would start their length-2, p = 3 iteration, but only after rerunning their Magma search over n < 35,448. We have not done this. The all-signature thresholds stay out of reach (26.9663 < 27.328 < 32).
   - Side remark: their GRH use in Prop. divisibility(b) is unnecessary where it is applied. In Prop. length1, n < 80 < 660 already excludes the nonsolvable case. In Prop. length2(b), K is totally real and Odlyzko's unconditional Table 2 gives rd ≥ 54.122 > 27 for n ≥ 600. This changes no theorem's status.

## Before
Odlyzko (1976 tables; Acta Arith. 1976) and Poitou (Sém. DPP 1976/77) give unconditional limits 4πe^γ = 22.3816 and 4πe^{1+γ} = 60.8395, the suprema of the method when zeros may lie anywhere in 0 < Re s < 1. Finite-n tables: Odlyzko 1976 (Table 2) and Diaz y Diaz 1980 (n ≤ 100; not consulted, see Caveat). Under GRH the limits are 44.763 and 215.333. Before October 2026 no zero-free strip of fixed width was known for every ζ_K.

## Method
- Poitou's explicit formula (P) (DERIVATION §1) gives log|d_K| = n(γ + log 8π − I_c(F)) + r1(π/2 − I_r(F)) − E(F) + Σ_ρ Φ(ρ) + 2Σ log N𝔭 · N𝔭^{−m/2} F(m log N𝔭), with E(F) = 4∫_0^∞ F(x) cosh(x/2) dx.
- Kernels are F(x) = S(x/b)·sech(wx), where S ∈ {G, H} are Odlyzko's 1976 shapes, normalised autocorrelations of cos(πt/2) and 1 + cos(πt) on [−1, 1].
  - F ≥ 0 drops the primes.
  - Under H(w), each zero contributes Re Φ(ρ) = FT[S(·/b)·cosh(a·)/cosh(w·)](γ) ≥ 0 with |a| ≤ w. This uses Lemma A: cosh(ax)/cosh(wx) has the positive Fourier transform (2π/w)cos(πa/2w)cosh(πt/2w)/(cosh(πt/w) + cos(πa/w)).
  - Optimality: positive definiteness forces F(x)cosh(wx) ≤ 1.
- The brief proposed kernels ∫e^{−(s−1/2)|x|}dμ(s), μ ≥ 0 on s > 31/32. This fails (DERIVATION §4):
  - Mass on s ≤ 1 makes E(F) infinite. Equivalently, −ζ_K'/ζ_K(s) = 1/(s − 1) − γ_K + … cancels the "positive" pole term 2/(1 − s).
  - Mixtures over s > 1 are no better than a single s (Stark), whose limits are 11.19 and 22.38.
- Parameters: b per (n, r1) comes from a float optimiser (make_kernels.py), using w' = w; a free w' ≥ w was never better (64 runs, logs/opt_wp.log). Every entry is then certified in Arb (certify.py). Within one kernel the bound increases with n, so each entry holds for all larger degrees.

## Certificate (one command)
`cd Riemann/Deep/Odlyzko31 && uv run --no-project --with python-flint python3 certify.py kernels.json`
It runs in 35 s and prints 4 `CROSSOVER … PASS` and 3 `THRESHOLD … PASS` lines, the limit enclosures and `ALL CERTIFIED` for 699 kernels and 700 entries ([logs/certify.log](logs/certify.log)). `python3 tables.py` then regenerates TABLES.md and tables.csv from certified.json.

## Checks
- **Formula (P), three independent routes** ([logs/checks.log](logs/checks.log)):
  1. With Odlyzko's kernels it reproduces all 322 (b, A, B, E) rows of his 1976 Tables 3 and 4. Ours exceed his "lower estimates" by 3·10⁻⁵ to 10⁻³ relative. Two GRH-table E entries (b = 15, 70) lie below his own closed form, a rounding slip or typo in his table that does not matter here.
  2. For K = ℚ, (P) holds as an identity with 300 Riemann zeros: residual −4.4·10⁻¹⁴ for H(x/2)sech(15x/32), and −3.8·10⁻⁹ against a tail bound of 8.6·10⁻⁹ for G(x/1.6).
  3. The Gamma terms, computed in t-space with mpmath's digamma, match Poitou's x-space form to 10⁻¹³ and 10⁻¹⁰.
- **Per-n reproduction.** Our classical (w = 1/2) and GRH (w = 0) recomputations match Odlyzko's Tables 2 and 1 at every listed n, slightly higher because we optimise b and integrate exactly.
- **Kernel class** ([logs/kernopt.log](logs/kernopt.log)). Optimising over nonnegative kernels h = q², with q a 3-term cosine series, gains at most 0.1% over H (n = 8) and ≤ 0.002% for n ≥ 1000. It reproduces Odlyzko's remark "5.656 unconditionally for n ≥ 8 with better kernels" (classical: 5.6563). So for n ≥ 8 the 31/32 gain (+0.18% at n = 8, +1.8% at 100, +3.9% at 10⁴) exceeds anything a better classical kernel of this type gives.
- **Known fields** ([logs/known_fields.log](logs/known_fields.log)). We took the smallest-|d| field per (degree ≤ 10, signature) and the smallest abelian ones from the LMFDB API: 42 fields, with 6 queries lost to rate limiting. All satisfy the certified 31/32 bounds, and the abelian ones also satisfy the B bounds. The tightest margins are 0.5% for ℚ(√5) and 0.6% for ℚ(√−3), and 1.1–1.3% for the cubic fields of discriminant −23 and 49.
- **Ghitza–Yamauchi values** were read in their TeX: l. 397 (27.328), l. 437–460 (Table 2 uses), l. 489 (n < 4800 from 32), l. 521 (n < 280 from 81), l. 567 (n < 4800 from 128).

## Novelty
[NOVELTY.md](NOVELTY.md):
- arXiv math.NT listings for October and November 2026;
- arXiv API (11 queries, including every paper mentioning Cureton or quasi-Riemann);
- the TeX of Cureton–Thorner, Frei–Sofos–Widmer and 2610.20600;
- Zenodo (7 phrases) and HAL (4).

None derives discriminant bounds from a zero-free strip.

## Caveat
- The improvement is a direct substitution of a 2026 input into a 1977 method; reviewers may call it routine. The new content is the numbers and the optimality of sech(wx) for general w.
- The inputs carry different status:
  - Cureton–Thorner's theorem is Lean-checked; their paper is unrefereed.
  - The explicit formula (P) is classical, not Lean-checked. Our checks confirm its constants numerically.
  - Every table entry is Arb-certified given (P).
- Diaz y Diaz's 1980 optimised unconditional tables were not available. For n ≤ 6 our gain over our own classical recomputation is ≤ 0.1%, so a better classical kernel might match it there; those degrees are fully known anyway. For 8 ≤ n < 2732 the record is relative to the best classical values known to us: Odlyzko 1976 and our optimised classical recomputation, whose kernel-class check reproduces Odlyzko's 5.656 at n ≥ 8. For n ≥ 2732 (all signatures) and n ≥ 1618 (totally real) it holds against any classical table, since the classical method stays below 4πe^γ and 4πe^{1+γ}.
- The abelian column is superseded for large n by the elementary conductor count (Claim 5).

RH relevance: N12, N7. A zero-free strip for every Dedekind zeta function (L-function analogue, N12; zero-free regions, N7) moves the classical discriminant bounds strictly between the unconditional limit and the GRH limit. The extremal kernel sech(wx) shows the bounds now depend continuously on the strip width: as w → 0 they tend to the GRH constants.

## Files
- [DERIVATION.md](DERIVATION.md): the formula, positivity, optimality, why not the s-family, abelian count, certification details.
- [NOVELTY.md](NOVELTY.md), [TABLES.md](TABLES.md), [tables.csv](tables.csv), certified.json, kernels.json.
- Scripts:
  - kern.py: float evaluation;
  - opt.py: w' study;
  - make_kernels.py: parameter choice;
  - **certify.py**: the certificate;
  - tables.py;
  - checks.py: the formula checks;
  - kernopt.py: kernel-class study;
  - known_fields.py.
- data/odlyzko1976_tables12.csv and data/odlyzko1976_tables34.csv (Odlyzko's public 1976 tables); data/lmfdb_min.json.
- logs/.
