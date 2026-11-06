# NOVELTY — selbergn4 (explicit Selberg-method proportion for linear combinations of N ≥ 4 Dirichlet L-functions)

Checked 2026-11-06, 13:10–13:30 UTC. Verdict: no source states an explicit proportion for N ≥ 4 functions above Dousselin's 3.70·10⁻⁹ (N = 4), and no source has an explicit c/N.

## Searches
1. arXiv monthly listings (this and last month), `python3 .scratch/lst.py 'linear combination|Selberg class|critical line|Dirichlet L-functions|Selberg.s method|mollif' math.NT math.CV math.CA
   --zenodo '"linear combinations" "critical line"'`: 14 title hits, relevant: 2610.26841 (Lamzouri, N = 2, 31%), 2610.08965 (ζ, 67.35% simple on the line), 2610.23467 (ζ, 8/9 simple or on line),
   2610.36885 (large values of Hardy's function, Selberg class; abstract read: no zero proportions). Nothing on N ≥ 3 functions. Zenodo: 2850 records, top 15 unrelated (flood models,
   RH "proofs", Hardy's function means), none on linear combinations.
2. arXiv API (all dates, `.scratch/axq.py`): au:Dousselin → 2311.10285 (the source), 2408.01401, 2303.05890 (class numbers; no 2026 follow-up);
   abs:"linear combination" AND "critical line" AND "zeros" → 2610.26841, 2508.11108 (short mollifiers for ζ), 2311.10285, 2010.10490 (Lamzouri–Lee, zeros near the line), 2007.04317,
   1602.06328, 1311.4862 (Hejhal), 1301.5144, 0704.3448, 0102007; abs:"Selberg" AND "linear combination" AND "proportion" → 2610.26841, 2311.10285 only;
   abs:"Selberg's method" AND "critical line" → 2609.15329 only (ζ, N = 1: 6.78% via Zhuravlev's Littlewood–Jensen form; abstract and intro read); abs:"Davenport-Heilbronn" AND "critical line" →
   five papers, none with a proportion.
3. HAL (api.archives-ouvertes.fr, q = Dousselin): the paper (hal-04302034, 2023), the thesis tel-05317244 (2025-06-19, "Corps cubiques, méthode de ramollissement et nombre de classes"),
   two class-number papers. Thesis Théorème I.1.2 (text extracted locally, lines 553–580): κ_F ≥ 7.41·10⁻⁹, 4.93·10⁻⁹, 3.70·10⁻⁹, 2.95·10⁻⁹, 1.47·10⁻⁹, 1.46·10⁻¹⁰, 1.44·10⁻¹¹ for
   N = 2, 3, 4, 5, 10, 100, 1000 (slightly above the paper's 7.38, 4.91, 3.68, 2.94, 1.46·10⁻⁹, 1.45·10⁻¹⁰, 1.43·10⁻¹¹), same 2.161·10⁻⁶/(N log N) large-N form (valid only for log N ≫ 862).
   So "Before" = thesis values.
4. GXL evidence API, one claim ("Selberg mollifier method gives an explicit positive proportion … of order c/N or c/(N log N)"): 246 passages; 1 "supports" (Rezvyakova 2501.00551, degree 2,
   non-explicit), 2 "refutes" (1505.00623, 2401.02813: Selberg's result does not cover other combinations), neutral: 2311.10285, 2610.26841, 2504.08322 (CLT for linear combinations of
   log L), 1206.1679, 1609.03738 — no explicit bound for N ≥ 3 beyond Dousselin.
5. Our files: DISCOVERIES.md / discoveries/: linear-combinations-line.md (N = 2, 34.50%), linear-combinations-three.md (N = 3, 2.061%; N = 2, 34.707%; "N ≥ 4 gives nothing");
   R406 proposed exactly this round; no earlier SelbergN4 claim.

## Relation to known results
- Selberg (1998, "Old and new conjectures…", and Collected Papers): κ_F ≥ c/N² (non-explicit) and the remark that c/(N log N) is attainable. Dousselin made c/(N log N) explicit.
  The c/N obtained here (κ_F ≥ 4.32·10⁻⁴/N for all N ≥ 11) is, as far as these searches show, not stated anywhere.
- For N = 2, 3 our LinComb3 bounds (34.707%, 2.061%; different method) are far better than the present Selberg-method values (1.62·10⁻⁴, 1.16·10⁻⁴); the record claim is for N ≥ 4.
- 2609.15329 (N = 1, ζ): a different, stronger form of Selberg's method (Zhuravlev: Littlewood's lemma + Jensen on the analytic half) with exact diagonal evaluation; not applied to
  linear combinations. Transferring it to the dominance sets S_j is the natural next step (REPORT, Caveat).
