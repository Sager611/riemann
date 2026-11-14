CORRECTED: C (Table 2) ×1.78–1.90 with both flaws fixed, ×1.81–2.05 with every slip found: 3·10^12 < T ≤ e^29: 37341.72 → 66414.87 (76414.26); e^500 < T ≤ e^1000: 6.23·10^6 → 1.17·10^7 (1.18·10^7); e^3000 < T ≤ e^481958: 5.76·10^8 → 1.06·10^9; B unchanged

# BellottiFix: Bellotti's explicit log-free zero density (arXiv 2405.12545v2, JNT 269 (2025)) with the flaws corrected

## Short answer
- Flaw (i) is real but sits in c3, not c4. Each distance band holds representatives on both sides of γ_r, so (4.11) (δ < |γ_r−γ_s| ≤ 1)
  must double: c3 → 2c3, and c3 is 82–93% of c2+c3+c4. (4.13)/c4 carries an unexplained factor 2.4. The two-sided count needs
  only 2 + 4δe^{1/6−π/2}/I ≤ 2.120, so the printed c4 is valid. Our earlier "c4 lacks a factor 2, +15.7%" (R160, R162, REFUTATIONS) is wrong.
- Flaw (ii) is real: |G(½+it)| ≤ d5 W/log W drops signs. The certified absolute-value bound 9.75707 W/log²W is smaller than the
  printed one at every row's W (log W ∈ [8.83, 14.90]): c5 drops 5–40% and C drops 0.6–5.9%.
- Falsification found three more slips, all making C bigger. (3.10) prints (u+w)^4 where (v+w)^4 belongs (its code does the same):
  +13% at T0 = 3·10^12, ≥ 0.5% while log T0 ≤ 33. (4.9) assumes X^y/y decreasing, which can fail in rows 26–38: c2 → (5F_A−4)c2, F_A ≤ 1.1362.
  The code's d3 is too small by a factor log²T0: d4 +0.15–1.1%. Smaller slips: (c)–(e) in NOTE §3.
- Bellotti's claim "sharpest explicit zero-density estimate uniformly for σ ∈ [α0, 1]" fails after the correction. By his own
  Table 1 margins, the corrected bound is weaker than Kadiri–Lumley–Ng / Johnston–Yang at (σ, T) = (0.993, 10^13), (0.990, e^46.2)
  and (0.985, e^90): −54…−75%. It still wins at the other 5 points (+24…+95%).
- Rafiq's "7 unreproduced rows" are not a real discrepancy. Rows 32 and 34–37 are printed to 3 significant figures, rounded up.
  Only rows 20–21 differ beyond rounding, and they are printed 0.27–0.29% above the code (the safe side). His dps and quadrature worry
  is harmless: dps 10 vs 30 changes C by < 5e-10.

## Corrected Table 2
"both flaws" = c3 → 2c3, the (4.9) repair needed by the corrected count (c2 → (5F_A−4)c2, F_A = 1 in rows 1–25), and c5 with the
certified |G| bound. "all fixes" = both flaws plus (3.10) (v+w)^4, d3·log²T0, √2 in b1, d1 at VW, e^{1/(6·0.447 log T)}, and an
honest partial summation for d2 (NOTE §3). Rounded up as Bellotti does (2 decimals; 3 significant figures above 10^6). Both
corrected columns come from code/indep.py at 50 digits; the corrected code matches them to 5e-10.

| # | T0 | T1 | alpha0 | B | C printed | C his code | C both flaws fixed | ratio | C all fixes | ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 3·10^12 | e^29 | 0.9927 | 14.160 | 37341.72 | 37341.72 | 66414.87 | 1.779 | 76414.26 | 2.046 |
| 2 | e^29 | e^30 | 0.9927 | 14.084 | 37564.18 | 37564.18 | 66849.95 | 1.780 | 75333.28 | 2.005 |
| 3 | e^30 | e^31 | 0.9925 | 13.588 | 40014.59 | 40014.59 | 71240.48 | 1.780 | 75998.39 | 1.899 |
| 4 | e^31 | e^32 | 0.9923 | 13.108 | 42840.01 | 42840.01 | 76305.96 | 1.781 | 79273.09 | 1.850 |
| 5 | e^32 | e^33 | 0.9921 | 12.647 | 46015.02 | 46015.02 | 82030.96 | 1.783 | 84129.86 | 1.828 |
| 6 | e^33 | e^34 | 0.9919 | 12.226 | 49270.71 | 49270.71 | 87885.92 | 1.784 | 89588.55 | 1.818 |
| 7 | e^34 | e^35 | 0.9918 | 11.856 | 52249.66 | 52249.66 | 93268.04 | 1.785 | 94790.68 | 1.814 |
| 8 | e^35 | e^36 | 0.9916 | 11.505 | 55524.80 | 55524.80 | 99202.68 | 1.787 | 100687.03 | 1.813 |
| 9 | e^36 | e^37 | 0.9914 | 11.134 | 59715.72 | 59715.72 | 106734.91 | 1.787 | 108282.25 | 1.813 |
| 10 | e^37 | e^38 | 0.9912 | 10.782 | 64228.93 | 64228.93 | 114908.23 | 1.789 | 116564.75 | 1.815 |
| 11 | e^38 | e^39 | 0.9910 | 10.475 | 68446.65 | 68446.65 | 122531.22 | 1.790 | 124303.64 | 1.816 |
| 12 | e^39 | e^40 | 0.9909 | 10.217 | 71907.32 | 71907.32 | 128782.01 | 1.791 | 130645.33 | 1.817 |
| 13 | e^40 | e^41 | 0.9907 | 9.951 | 76267.80 | 76267.80 | 136683.00 | 1.792 | 138682.14 | 1.818 |
| 14 | e^41 | e^42 | 0.9906 | 9.720 | 80004.48 | 80004.48 | 143445.17 | 1.793 | 145552.46 | 1.819 |
| 15 | e^42 | e^43 | 0.9905 | 9.501 | 83811.09 | 83811.09 | 150399.98 | 1.795 | 152619.73 | 1.821 |
| 16 | e^43 | e^44 | 0.9903 | 9.259 | 89209.33 | 89209.33 | 160180.30 | 1.796 | 162584.28 | 1.823 |
| 17 | e^44 | e^45 | 0.9901 | 9.020 | 95233.68 | 95233.68 | 171139.66 | 1.797 | 173749.88 | 1.824 |
| 18 | e^45 | e^46.2 | 0.9900 | 8.843 | 99227.30 | 99227.30 | 178422.24 | 1.798 | 181157.94 | 1.826 |
| 19 | e^46.2 | e^50 | 0.9899 | 8.634 | 104761.04 | 104761.04 | 188633.01 | 1.801 | 191523.19 | 1.828 |
| 20 | e^50 | e^60 | 0.9892 | 7.936 | 128806.85 | 128461.39 | 232021.62 | 1.801 | 235786.02 | 1.831 |
| 21 | e^60 | e^70 | 0.9876 | 6.661 | 194323.75 | 193755.46 | 352330.44 | 1.813 | 358526.59 | 1.845 |
| 22 | e^70 | e^80 | 0.9860 | 5.809 | 278045.07 | 278045.07 | 508089.00 | 1.827 | 517548.14 | 1.861 |
| 23 | e^80 | e^90 | 0.9850 | 5.216 | 370655.73 | 370655.73 | 680293.43 | 1.835 | 693330.58 | 1.871 |
| 24 | e^90 | e^100 | 0.9850 | 4.831 | 425721.47 | 425721.47 | 781856.52 | 1.837 | 796256.77 | 1.870 |
| 25 | e^100 | e^110 | 0.9850 | 4.513 | 488901.14 | 488901.14 | 896954.68 | 1.835 | 912981.94 | 1.867 |
| 26 | e^110 | e^120 | 0.9850 | 4.264 | 545744.21 | 545744.21 | 1.01·10^6 | 1.834 | 1.02·10^6 | 1.866 |
| 27 | e^120 | e^130 | 0.9850 | 4.032 | 629490.27 | 629490.27 | 1.16·10^6 | 1.838 | 1.18·10^6 | 1.869 |
| 28 | e^130 | e^140 | 0.9850 | 3.855 | 694045.43 | 694045.43 | 1.28·10^6 | 1.838 | 1.30·10^6 | 1.868 |
| 29 | e^140 | e^150 | 0.9850 | 3.696 | 771373.78 | 771373.78 | 1.42·10^6 | 1.837 | 1.44·10^6 | 1.866 |
| 30 | e^150 | e^160 | 0.9850 | 3.572 | 825913.07 | 825913.07 | 1.52·10^6 | 1.839 | 1.55·10^6 | 1.867 |
| 31 | e^160 | e^170.2 | 0.9850 | 3.448 | 909966.07 | 909966.07 | 1.68·10^6 | 1.840 | 1.70·10^6 | 1.867 |
| 32 | e^170.2 | e^500 | 0.9850 | 3.337 | 1.12·10^6 | 1.12·10^6 | 2.07·10^6 | 1.839 | 2.10·10^6 | 1.867 |
| 33 | e^500 | e^1000 | 0.9850 | 2.152 | 6.23·10^6 | 6.23·10^6 | 1.17·10^7 | 1.867 | 1.18·10^7 | 1.884 |
| 34 | e^1000 | e^1500 | 0.9850 | 1.820 | 2.12·10^7 | 2.12·10^7 | 3.94·10^7 | 1.856 | 3.96·10^7 | 1.867 |
| 35 | e^1500 | e^2000 | 0.9850 | 1.684 | 6.47·10^7 | 6.47·10^7 | 1.23·10^8 | 1.896 | 1.24·10^8 | 1.904 |
| 36 | e^2000 | e^2500 | 0.9850 | 1.610 | 1.73·10^8 | 1.73·10^8 | 3.28·10^8 | 1.890 | 3.29·10^8 | 1.899 |
| 37 | e^2500 | e^3000 | 0.9850 | 1.577 | 2.56·10^8 | 2.56·10^8 | 4.83·10^8 | 1.883 | 4.85·10^8 | 1.891 |
| 38 | e^3000 | e^481958 | 0.9850 | 1.551 | 5.76·10^8 | 5.76·10^8 | 1.06·10^9 | 1.829 | 1.06·10^9 | 1.835 |

## The corrections (proofs in NOTE.md §1–2)
- (4.11), l.843–849, per r: Σ_{d≤1} X^{Re z}|Γ(z)| ≤ 1.0067 X^{2δ} δ^{-1}(½ + 2H_{⌊1/δ⌋} + 2.5(F_A−1)), with H_{⌊1/δ⌋} ≤ 1.12 log(1/ν_i(T1))
  checked in every row. Printed: ½ + H_{1/δ−1}. Reason: the k-th nearest representative on each side has d ≥ (k−1)δ.
- (4.13), l.858–871, per r: Σ_{d>1} ≤ X^{2δ}√(2π)(2I/δ + 4e^{1/6−π/2}) ≤ X^{2δ}√(2π)·2.4·I/δ (printed), since 2 + 4(1−α0)e^{1/6−π/2}/I ≤ 2.1199.
- l.933: |G(½+it)| ≤ Σ|θ_jθ_k|[j,k]^{-1/2} ≤ c²P·W/log²W ≤ 9.75707 W/log²W (arb), where [j,k]^{-1/2} = (jk)^{-1/2}Σ_{r|(j,k)}h(r),
  h(p) = √p−1, c = sup S(y)/√y ≤ 2.574121012 (orchestrator's reviewer, 14 Nov: rounding repaired; Arb ball upper end 2.5741210108618) and P ≤ 1.4725205, so c²P < 9.757067. A cruder 16P ≤ 23.5603 would raise C by 3.5–17% instead.
- Reviewer, first pass: A (c3), B (c4 valid), D ((v+w)^4) CORRECT; C GAP only because the numerical constants were not yet certified.
  cert.py (arb) then certified them. Second pass: E (F_A) and F (certification) both GAP, on stated bounds only. E: F_A(0.576) =
  1.13615 > 1.1361, now stated as ≤ 1.1362, and the step fails only for some δ. F: propagate the float error, and quote c²P ≤ 9.757067.
  Both are fixed in NOTE §1–2. The computations used the exact F_A and 9.7572, so no number changes.

## Spot checks at 50 digits (code/indep.py, built from the TeX formulas rather than his code; logs/spot50.json)
- Row 1 (3·10^12, e^29]: code-like 37341.716328 = his code. Corrected: 66414.865421 (c3 4445.43, c5 0.03913). All fixes: 76414.258109
  (c1 0.98122 → 0.92256 from (v+w)^4).
- Row 33 (e^500, e^1000]: his code 6229979.24. Corrected: 11632252.85 (F_A 1.0966, c2 30905 → 45831, c3 ×2). All fixes: 11737811.63.
- Over all 38 rows, the code-like mode matches his code to 2.9e-10. Every row passes its checks: F4need ≤ 2.12 < 2.4,
  H ≤ 1.12 log(1/ν), and c1² − c5 ≥ 0.81 > 0.

## RECORDS-ready lines (drafts)
- (N7) Explicit log-free zero density for ζ: N(σ,T) ≤ C T^{B(1−σ)} for σ ≥ α0 on Bellotti's 38 T-ranges (3·10^12 < T ≤ e^481958).
  It holds with his B and α0 and with C ×1.81–2.05 the printed values. Examples: C = 76414.27, B = 14.160, α0 = 0.9927 for
  T ≤ e^29; C = 1.06·10^9, B = 1.551, α0 = 0.985 for T ≤ e^481958. The printed C (JNT 269 (2025)) is not proved. [BellottiFix]
- (N7) Against Kadiri–Lumley–Ng / Johnston–Yang, the corrected estimate loses at (0.993, 10^13), (0.990, e^46.2) and (0.985, e^90).
  It keeps margins of 24–95% at (0.993, e^46.2), (0.99/0.993, e^90) and (0.99/0.993, e^170.2). No other explicit log-free estimate
  for ζ is known: Thorner–Zaman 2208.11123 covers Dirichlet L-functions with far larger constants. arXiv searches turned up nothing
  newer (2412.02068 Carlson, 2507.15184 Ingham and 2609.23222 Fiori–Jaskari are not log-free; the last uses the KLN form).

## Effect on our lines that cite Bellotti
- REFUTATIONS.md l.9 (Bellotti): replace "factor about 2 in c4 … 15.7%". New text: (4.11) undercounts, so c3 ×2; c4 stands (2.4);
  the G-bound proof is invalid but the corrected bound is smaller; (3.10) (v+w)^4 adds up to 13%; (4.9) monotonicity fails in rows 26–38;
  corrected C ×1.81–2.05 (this table). "Theorem not disproved" stays: corrected constants exist.
- REFUTATIONS.md l.10 and Riemann/Deep/Ramanujan (Padhi): "inherits the factor-2 undercount of Bellotti (4.13)" is unsupported,
  because (4.13) holds as printed. Padhi's far term 2√(2π)·2.4·I′ uses a one-sided I′, and its 2.4 may cover both sides as Bellotti's does.
  Our q = 1.03–1.07 came from doubling on top of 2.4. This needs a re-audit, not done here. Open points: what Padhi's I′ is exactly;
  whether his N3 (Bellotti's tables, understated ×1.8–2.0) binds near log x = 875; his G-step is harmless with 9.757 W/log²W (wL ∈ [15, 24]).
- TREE.md R160, R162 (c4 / 15.7%): superseded by this report. R164: unchanged in substance.
- RamRepair (Ramanujan N_R, ours): unaffected. It doubled the far term and used 16P, both conservative, and it never used N3.
- Other Bellotti papers we cite are not touched: Bellotti–Kerr (COROLLARIES l.49), the Piltz database (Combos), Bellotti–Wong
  (STCAudit), BTY 2603.21490 (R279, RamRepair).

## Files and runs
- code/constants_orig.py: his constants.py, verbatim. code/constants_fix.py: a corrected copy with every change marked FIX
  (modes fix, fix16, r162, both, c3only, c5only). code/harness.py runs one copy per table row. code/indep.py: the independent
  50-digit implementation. code/prodP.py, sqfree.py, cert.py: the G constants. code/table.py, kln.py, mono.py, split.py: summaries.
  code/runall.sh runs everything.
- logs/: orig.json (31/38 exact, the other 7 explained above), fix*.json, indep_all.json, spot50.json, table.md, kln.log, mono.log,
  cert.log, plus memcap and /usr/bin/time peaks: every run < 141 MB RSS (cap 290 MB), one process at a time.


RH relevance: N7
