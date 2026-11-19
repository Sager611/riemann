AUDIT PASSED: Fodge's Appendix C holds as written — Lemmas C.2–C.4 and Prop. C.1 re-derived line by line (the explicit B-process of C.3 included; reviewer CORRECT on its stationary-phase step), every constant recomputed in Arb, 0 of ≈16,000 intermediate-inequality tests, 32,374 C.4 block sums, 44 heights of C.1 bookkeeping and 397 C.5 sums exceeded; C.5–C.8 and C.7's proof and code also hold, with only C.7's full 421-interval table not rerun (first, middle and last intervals rerun, identical).

# FodgeAudit — referee report on Appendix C of Fodge, Zenodo 24325329 v40 (17 Nov 2026), nodes N9 → N7

Source: `zfr_exact_stechkin_v40_2026-11-17.pdf` (36 pp., pdftotext lines l., page images read at 130 dpi), code bundle `code_exact_stechkin_v31.zip`. Inputs re-read: HPY arXiv 2207.02366 TeX (Lemma 2.1, l. 150–163). Full step-by-step chain with line numbers: [FULLPROOF.md](FULLPROOF.md). Earlier checks not redone: R482 (Lemma C.3 on 26,500 exact sums, max 0.547), R484 (N7 certificate), RevAudit (Revers 2602.05614).

## 1. Line-by-line referee (statements as printed)

| Item | Statement (pdf) | Verdict | Remarks |
|---|---|---|---|
| C.2 l. 1979–1995 | weighted Weyl differencing, any w_j ≥ 0 | CORRECT | re-derived (index ranges, Cauchy–Schwarz, conjugate pairing needs r(m) ≥ 0) |
| C.3 l. 2012–2136 | \|Σ_{n<L′−m}e(g(n))\| ≤ µ′L(m/W)^{1/2} + 4.524(W/m)^{1/2} + 1.0004 + 2.24m/W, N ≥ 10⁴, L ≤ N/75 | CORRECT | all 7 steps re-derived (FULLPROOF §C.3): Euler/Poisson with ρ vanishing at half-integers; E-frequencies 4(πλ)^{−1/2}; two partial integrations off the stationary range; explicit stationary phase (u = yq^{1/2}, τ ≤ 1/(2π²U³), b1′ = u^{−2}∫s a″, the five a_yy terms ⇒ Ξ); Euler–Maclaurin for the main terms; π/sin closed form for the boundary sums |Σ_{j∉E}B_j| ≤ ½ |
| C.4 l. 2137–2177 | \|Σ_{n<L′}e(f(n))\|² ≤ (L/W^{1/3}+c)(αL+βW^{2/3}) | CORRECT | flat weights only; the four m-sums; term-by-term passage to α, β |
| C.1 proof l. 2178–2283 | 0.3080532 t^{1/6}log t + 2.405897 t^{1/6} + 506.6, t ≥ 2.9·10¹² | CORRECT | RS input = HPY Lemma 2.1 verbatim (t ≥ 200, 1.48t^{−1/4} + 0.127t^{−3/4}); blocks, µ′ = ((N+L)/(N′+L))^{3/2}φ(L/N′), ᾱ, β̄, e1, e0, κ, κ″, C′, I∞ all re-derived |
| C.5 l. 2284–2297 | C.3 for general g | CORRECT | C.3's proof uses only the listed properties; C⁴ needed for a_yy |
| C.6 l. 2298–2334, C.8 l. 2344–2390 | second/third differencing | CORRECT (proof); constants see §3 | phase derivatives, λ, h, ε3, ε4, µ″, µ3 re-derived |
| C.7 l. 2335–2343, proof l. 2391–2447 | U(t) on 421 intervals of log t ∈ [28.69, 87.25] | CORRECT (proof and code logic) | monotonicity in s = t^{1/3} re-derived for every term incl. µ′, µ″, µ3; code mirrors the proof; coverage contiguous |

Inequality directions: all checked. Two places where the printed form is weaker than what the argument gives (harmless): F(β)−F(α) ≤ (1−h^{−1/2})λ^{−1/2} ≤ (h^{1/2}−1)λ^{−1/2}; Σ(1−m/M)m = (M²−1)/6 ≤ M²/6. One place needing the right argument: u_b ≥ (β−j)/(hλ)^{1/2} holds via u_b²/2 = ∫_0^{β−j} v dv/ψ″ (the naive y_b-bound gives only (β−j)/(hλ^{1/2})). One sign subtlety: C′ = log((1+u⁺)/(√(2π)ν)) = −0.855 < 0, so b2 decreases in α∞; the code evaluates b2 as one interval, so the direction is right.

Independent reviewer (subagent `reviewer`, pure reasoning, the stationary-phase step of C.3 pasted inline with the exact formulas: (a) the decomposition I_j = a0e(ψ(x_j)+1/8) + B_j(b) − B_j(a) + R_j and the R_j bound, (b) u_b, |u_a| lower bounds, (c) b1′ = u^{−2}∫s a″ and a″ = a(a_y² + a a_yy), (d) the Ξ bound, (e) the count of stationary j and the 0.00027 total): **VERDICT: CORRECT** on all five; for (b) it used |ψ′|² ≤ 2hλ(ψ − ψ(x_j)) = hλu², equivalent to ours.

## 2. Constants recomputed (Arb, python-flint, 256 bits; `consts.py`, `consts2.py`, logs/)

| Constant | Stated ingredients | Recomputed | Printed | Margin |
|---|---|---|---|---|
| h (C.3) | ((1+1/75)·10⁴/(10⁴−½))³ | 1.0406918 | ≤ 1.0407 | 8.2·10⁻⁶ |
| ε3 | 3/(10⁴−½) | 3.000150·10⁻⁴ | ≤ 3.0002·10⁻⁴ | 5·10⁻⁹ |
| 8/√π + (√h−1)/2 | h = 1.0407 | 4.5235902 | 4.524 | 4.1·10⁻⁴ |
| h^{5/2}ΞL²/(4π) | L/N′ ≤ (1/75)·10⁴/(10⁴−½) | 2.63619·10⁻⁴ | 0.00027 | 6.4·10⁻⁶ |
| 1 + ε3/8 + amplitude | | 1.0003011 | 1.0004 | 9.9·10⁻⁵ |
| 7ζ(3)h^{3/2}/π² + 21ζ(3)h/(2π²) | | 2.2360152 | 2.24 | 4.0·10⁻³ |
| µ′ = ((N+L)/(N′+L))^{3/2}φ(L/N′) | identity | |diff| ≤ 7·10⁻⁷² on 2000 random (N, L) | — | exact |
| u⁺, µ′⁺, κ, α∞, β̄, e1, e0, C′ (C.1) | proof formulas | 0.0125738148, 1.0094669225, 1.0000195170, 1.2508161940, 9.2745371544, 1122.67312, 395.244733, −0.8551498 | as printed | all ≤ printed |
| I∞ | ∫_{π^{1/3}ν}^∞(ᾱ^{1/2}−α∞^{1/2})dy/y | 0.0529032 (acb.integral to e⁴⁰ + tail) | ≤ 0.0531080 (code: upper sums) | 2.0·10⁻⁴ |
| b1 | κα∞^{1/2}/(3π^{1/6}) | 0.308053152008 | 0.3080532 | 4.8·10⁻⁸ |
| b2 | 4√ν + (2κ/π^{1/6})(α∞^{1/2}C′ + I∞) | 2.4055585 (exact I∞); 2.40589694 (code's I∞ bound) | 2.405897 | 3.4·10⁻⁴; 6·10⁻⁸ |
| S | π^{1/6}e1(1+u⁺)²/(2πα∞)^{1/2} + e0κ″((1/6)log T1 + C′)/(π^{1/6}α∞^{1/2}T1^{1/6}) | 506.575070 | 506.6 | 0.025 |
| variant ν = 1.15, t ≥ e^{27.4} | same | 0.3080562827, 2.4285868, 485.309030 | 0.3080571, 2.428866, 485.31 | 8.2·10⁻⁷, 2.8·10⁻⁴, 9.7·10⁻⁴ |

## 3. Code audit

- `appendix_c/appC_sp_consts.py` (59 lines; = `appC_sp_consts_v2.py` with T1 fixed): every formula matches the proof (FULLPROOF); mpmath `iv` at 40 digits, outward rounding; ζ(3) replaced by 1.2020570 > ζ(3) (all constants increase in it); h replaced by 1.0407 after asserting the bound; N0 = νT1^{1/3} is a valid lower bound for every block's N; u⁺, N/N′, µ′⁺ evaluated at the worst N; I∞ by left-endpoint upper sums of the decreasing integrand on a 4000-step geometric grid to 10⁷ plus the tail A/(4√c α∞^{1/2}·10⁷), evaluated at lo(y) (larger) and with interval widths (y_{k+1} − y_k) — valid; b2 computed as one interval, so C′ < 0 is handled; asserts (log T1)/6 + C′ > 1 and RS(T1) < 2. Reran by R482 (identical log); my flint recomputation agrees to all printed digits (b1 identical to 12 digits; b2 differs only through the exact I∞).
- Coverage in t: C.1 is analytic and uniform for t ≥ T1 — every t-dependence enters through monotone quantities evaluated at T1 (N0, u⁺, T1^{−1/3} in e1, ((1/6)log t + C′)t^{−1/6} decreasing beyond T1, RS remainder); no scan, no gap (lesson 277 n/a).
- `appendix_c2/propC1p_v6.py` (Prop. C.7, 150 lines, python-flint balls at 80 bits): block k taken iff it can occur for t ≤ t_b (s ≥ 2πν_a²ρ^{2k}), evaluated at s_k = max(s_a, lo(2πν_a²ρ^{2k})) with y_k, l_k as balls from the ceilings; asserts l_k/y_k ≤ 1/75; methods gated by lower bounds N ≥ 10⁴ (C.2+C.3), 5·10⁴ (C.6), 2·10⁴ (C.8); µ′, µ3 in cancellation-free closed forms (checked algebraically); shift groups take each increasing term at the group top and each decreasing term at its bottom, with exact group weight sums (parabola closed forms exact: `consts2.py`); free integers chosen in floating point (admissible); U rounded up and asserted > ball. Matches the proof.
- Coverage of Prop. C.7: `run_table_v6.py` chains intervals through identical 4-decimal endpoint strings; `prof_table_v6.txt` has 421 contiguous intervals 28.6900 → 87.2500 (no gap, no overlap); the N7 certificate needs C.7 only on [ℓ1, 86.45] and uses C.1 beyond 86.45 (l. 1318–1322).
- C.5 constants as used by C.6 (N ≥ 5·10⁴: h ≤ 1.0544517, 4.5269492, 1.0005007, 2.2716010 ≤ printed 4.52695, 1.00051, 2.27161 ≤ used 4.527, 1.0006, 2.272) and C.8 (N ≥ 2·10⁴: h ≤ 1.0686019, 4.5303827, 1.0008414, 2.3083408 ≤ printed 4.53039, 1.00085, 2.30835 ≤ used 4.531, 1.001, 2.31): recomputed in Arb, all hold, worst at the smallest N (rechecked at 2N). Parabola-weight closed forms (`parab_forms.py`: W1, r(m), group sums Σ_{m=a}^{b} r(m)): exact, 0 mismatches (M < 60 all m; 40 random M ≤ 400).

## 4. Falsification of the intermediate inequalities of the C.3 proof (`falsify.py`, logs/gridA*.jsonl, `summarize.py`)

Exact data: phase g(x) = (t/2π)(log(N+x+m) − log(N+x)) with Taylor coefficients at the midpoint in mpmath (40 digits; g(x0) mod 1, g′(x0) − j exact), the rest in double; exact sums in mpmath; I_j = ∫_a^b e(ψ_j) by composite 20-node Gauss–Legendre (≤ 2 cycles per panel); x_j in closed form. 300 cases: N = 10⁴, L = 133 = ⌊N/75⌋ (the corner), λ = m/W ∈ {10⁻⁴, 10⁻³, 0.01, 0.05, 0.2, 0.5, 1, 1.1, 2} × (L′, m) ∈ {(133,1), (133,2), (133,66), (133,131), (133,132), (40,7), (3,1), (2,1)} × {untuned; β = k+½+10⁻⁶ (stationary point 10⁻⁶ from the excluded zone at b); α = k−½−10⁻⁶ (same at a); ε = β − j_β = ½ − 10⁻⁶ (largest π/sin term)}, plus the actual first block of Prop. C.1 at t = T1 (N = 13548, L = L′ = 170, m = 1, 2, 3, λ = 0.358, 0.715, 1.073) × the four modes. Per case: every exceptional j, every nonstationary j within 20 of [α, β], every stationary j (when > 40 of them: the 15 nearest each end + 10 interior; the main-term sum M1 always over all). ≈ 600 exceptional, ≈ 12,300 nonstationary, ≈ 3,300 stationary instances.

| Check (displayed inequality) | Max LHS/RHS | Where |
|---|---|---|
| E1 \|I_j\| ≤ 4(πλ)^{−1/2}, j ∈ E | 0.5577 | L′−m = 2, λ = 0.5, ε → ½ |
| N1 4π²\|R_j\| ≤ g″(a)/\|α−j\|³ + g″(b)/\|β−j\|³ + ∫\|(g″/ψ′³)′\| (nonstationary; TV from a grid = lower estimate of the RHS) | 0.999992 | b−a = 1, λ = 10⁻³: sharp (the two boundary terms align when ψ′ ≈ half-odd and the phase is nearly linear); never exceeded |
| N2 the same with hλ, ε3/(2d²), hλ/d³ | 0.6593 | L′ = 2, λ = 0.01 |
| S1 \|R_j\| ≤ (a0/2π²)(u_b^{−3}+\|u_a\|^{−3}) + (1/2π)∫\|b1′\| (stationary) | 0.4998 | λ = 1.1, L′ = 40, m = 7 |
| S2 per-j form h^{3/2}λ/(2π²)((β−j)^{−3}+(j−α)^{−3}) + h^{3/2}Ξ(b−a)/(4πλ) | 0.4893 | C.1 block, λ = 1.073 |
| S3 u_b ≥ (β−j)/(hλ)^{1/2}, \|u_a\| ≥ (j−α)/(hλ)^{1/2}, a0 ≤ λ^{−1/2} | 0.99991 | sharp (g″(x_j) → λ), never exceeded |
| B1 \|Σ_{j∉E}B_j(b)\|, \|Σ_{j∉E}B_j(a)\| ≤ ½ (closed form, itself checked against 4·10⁵-term symmetric sums) | 0.5756 | λ = 1.1, ε → ½ |
| M1 Σ_stat g″(x_j)^{−1/2} ≤ µ′L(m/W)^{1/2} + (h^{1/2}−1)/(2λ^{1/2}) | 0.9915 | λ = 1.1, L′ = 133, m = 1 (nearly sharp, as expected) |
| D1 identity: S − Σ_E I_j − main − boundary − Σ_window R_j vs the N2 tail bound outside the window | 0.0304 (there residual 1.1·10⁻⁶ against tail bound 3.7·10⁻⁵) | 246 non-subset cases: the whole decomposition (Poisson, e(1/8), signs of B_j) is exact |
| C3 the lemma itself | 0.2691 | consistent with R482's 0.547 (other corners) |

No ratio exceeds 1. Two intermediate bounds are attained asymptotically (N1, S3) and one nearly (M1): the proof has no slack there but no error either.

## 5. Falsification of C.4, the C.1 bookkeeping, and C.5–C.8

- **C.1 block by block** (`bookkeeping.py blocks`, mpmath 30 digits, exact integer blocks N_k, L_k): 44 heights t ∈ {T1, T1(1+10⁻⁷), 3·10¹², 10¹³ … 10⁴⁰} ∪ 30 log-random in [T1, 10⁸T1] (315–1166 blocks each). Every per-block claim of the proof holds at every block (max ratios, each must be ≤ 1): 1/(u_kr) 1 (equality when 80 | N_k, allowed), u_k/u⁺ 0.99951, 10⁴/N_k 0.738, 75L_k/N_k 0.943, y0/x_k ≤ 1, µ′_k/µ′⁺ 0.9999915, α_k/ᾱ(x_k) 0.9999941, β_k/β̄ ≤ 1 (equality only in the limit x_k → y0), E_k/((e1x_k + e0)t^{−1/3}) 0.99606, u/√(1+u)/(κlog(1+u)) 0.9999870, u/(κ″log(1+u)) 0.9999970, N_end/((1+u⁺)n1) 0.99984. The sum head + Σ_k 2N_k^{−1/2}((L_k/x_k + c)(α_kL_k + β_kx_k²))^{1/2} + RS(t) against b1t^{1/6}log t + b2t^{1/6} + S: ratio 0.8915 at T1, rising to 0.9990 at 10²⁵, 0.99958 at 10³⁵, 0.99928 at 10⁴⁰ — never above 1. The closing estimates are asymptotically tight (relative slack ≈ 4·10⁻⁴ at 10³⁵), so the printed constants cannot be lowered much without re-doing the integral comparison; no error.
- **C.4 on exact block sums** (`bookkeeping.py c4 100`): 32,374 blocks (every block up to n1, the last one cut) at 100 random t ∈ [T1, 3T1], every prefix length L′ ≤ L: max |Σ_{n<L′}e(f(n))|²/RHS = 0.0488.
- **C.5 as used by C.6 / C.8** (`c5test.py`, mpmath 45 digits): exact sums of second differences (N = 5·10⁴, L = 666, constants 4.527/1.0006/2.272, main term µ″L(mk/V)^{1/2}), 199 cases, λ = mk/V log-uniform in [3·10⁻⁴, 2.5]: max 0.3746; third differences (N = 2·10⁴, L = 266, 4.531/1.001/2.31, µ3L(mkj/V3)^{1/2}), 198 cases: max 0.3486.
- **Prop. C.7**:
  - Rerun of `propC1p_v6.py 50.0 50.05` (python-flint, 175 s, peak 38 MB): U = 14.24639, 1042 blocks (single 380, double 414, triple 248), parts head 0.25298, triple 0.71613, double 4.55374, single 8.72304, remainder 0.00048 — identical to the bundle's `parts50.log` and to the table and the text (p.35).
  - Rerun of the first interval `propC1p_v6.py 28.6900 28.8900` (11 s): U = 13.38796, 287 blocks all single-differenced, head 4.74396, single 8.62722, remainder 0.01676 — identical to the table's first line (coverage starts at 28.69 ≤ log T1 = 28.6957).
  - Rerun of the last interval `propC1p_v6.py 87.0000 87.2500` (270 s, peak 30 MB): U = 18.79375, 1542 blocks (single 547, double 779, triple 216) — identical to the table's last line.
  - The key new step of C.7's proof (B_k/s² non-increasing in s for fixed y, l and integers), tested numerically (`c7mono.py`): 360 random (y, l, M, K, J) for each of C.2+C.3, C.6, C.8 bounds, 8 values of s ∈ [s_k, 50s_k]: B/s² strictly decreasing everywhere (largest relative change between consecutive s: −1.3·10⁻⁶).

## 6. What holds, what survives, dependents

- Holds as written: Lemma C.2, Lemma C.3 (every step and constant), Lemma C.4, Proposition C.1 (0.3080532 t^{1/6}log t + 2.405897 t^{1/6} + 506.6 for t ≥ 2.9·10¹²) and its variant (0.3080571, 2.428866, 485.31 for t ≥ e^{27.4}), Lemma C.5, Lemmas C.6 and C.8 (constants included), Proposition C.7 (proof and code logic; numerically one interval rerun, §7).
- No step fails, so nothing needs a fallback. For the record, the paper's own fallbacks: the N7 certificate holds at R₀ = 4.6892 with C.7 as stated and at 4.7079 with C.1 alone (l. 1573–1578); both beat Bellotti–Trudgian–Yang's 4.896. The N7 record R₀ = 4.6892 survives this audit as far as Appendix C is concerned.
- Sharpness found (not errors): N1 (the two-partial-integrations remainder) and S3 (u_b, |u_a|, a0 bounds) are attained asymptotically; M1 (Euler–Maclaurin for the main terms) to 0.9915; the C.1 closing estimates to 0.9996 of the block-by-block sum at t = 10³⁵.
- Our dependents (`rg -n -i 'fodge|24325329'` over TREE.md COROLLARIES.md RECORDS.md REFUTATIONS.md NOGO.md DISCOVERIES.md discoveries Riemann/Tree): TREE.md l. 11 (N7 record line: R₀ = 4.6892, via Prop. C.1), l. 13 (N9 explicit record: Prop. C.1), l. 518 (R482 index), l. 520 (R484 index); Riemann/Tree/R482.lean, R484.lean. None in RECORDS.md, COROLLARIES.md, REFUTATIONS.md, NOGO.md, DISCOVERIES.md, discoveries/. The N9 line's "the B-process proof, l. 2021–2136, not refereed line by line" can now read "refereed line by line (FodgeAudit)".

## 7. Unchecked residue

- The full Prop. C.7 table (421 intervals, ≈ 2.3 h at the author's speed, longer here) was not rerun; the first, [50, 50.05] and the last interval were (§5, identical); the rest rests on the audited code and the bundle's log (`prof_table_v6.txt`, used by R484).
- Gabcke's 0.127 t^{−3/4} under HPY Lemma 2.1: imported, as in RevAudit (its residual).
- Correctness of mpmath `iv` and Arb is assumed.

RH relevance: N9, N7
