# Fourth moment below κ = 3/4: working notes (2026-10-21)

Brief: `.scratch/orca/kappa-brief.md`. Verdict: [REPORT.md](REPORT.md), NO REFUTATION.

## Sources
- **OpenAI 7/8:** `preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex` (16,677 lines), fetched from github.com/openai/math via the GitHub API and raw.githubusercontent.
- **Xu–Sang 2610.15291v7:** `seven_eighths_and_beyond_v7.tex` (2,894 lines).
- **Liu 2610.12234:** `combine_v2-2.tex` (6,338 lines).
- **Liu 2610.26835:** `number_fields_v2.tex` (11,410 lines).
- The three arXiv sources come from the e-prints. All downloads lived in .scratch/kappa/, which is now deleted.
- **Lean:** our port Riemann/Deep/QRH/OAI/NumberTheory/DirichletL (openai/math commit 1581d4cb).

## 1. Lemma 18.1 and its proof (OpenAI §18, l. 12477–14984)
- **Statement** (l. 12531–12579): 3/4 ≤ κ ≤ 1. The z > 0 case requires n1 + n2 + 6κz ≤ M and, when κ < 1, β_* ≤ (1+κ)/2. The mesh is "uniform for κ ∈ [3/4,1]".
- **Supporting lemmas:**
  - Lemma 18.2, "Common coefficient under complete extraction" (l. 13953–13996): no κ.
  - Lemma 18.3, "Masked rectangle cancellation" (l. 14545–14680): no κ.
- **Search.** Every occurrence of κ, 3/4, 7/2, M/21, M/20 and 2M/9 in §18 was found with rg. Each is in the table below.
- **z elsewhere.** Outside the table, the slot length z occurs only in mask erasure (l. 12651, where z decreases) and z' ≤ z (l. 13831).

| Line(s) | Use of κ | At 13/18 ≤ κ < 3/4 |
|---|---|---|
| 12833–12898 (3.4)–(3.5) | Contour shift to s_κ = (1+κ)/2 + e. Bound Z^{κz+ε}. Losses "uniform in [3/4,1]". | OK. Lemma logarithmic-control (l. 1531) holds for every b ∈ [1/2,1], and s_κ ≥ 31/36. Prime powers need s_κ ≥ 1/2. |
| 12929–12936 (3.6) | z ≤ M/(6κ) ≤ 2M/9, and κz ≤ M/6 | The 2M/9 part is false (it becomes 3M/13), but label (3.6) is never cited again. Only κz ≤ M/6 is used, for the terminal loss ρ/6. |
| 12947–12950 | 6κ−1 ≥ 7/2 ⇒ z < M/21 < M/20, and A − z > 2L | 7/2 and M/21 are false. z < M/20 ⟺ κ ≥ 13/18 (since (6κ−1)z < M/6). A − z > 47M/60 > M/2. |
| 12958–12984 | A_comp ≤ 3M/2 − A + 2z + ξ ≤ 23M/30 + ξ; plus (6κ−1)z it is ≤ 14M/15 + ξ; margin M/15, with ξ ≤ ρ/30 | Unchanged: these lines already use M/20 only. Tight at 13/18. |
| 14214–14280 (3.14), (3.16) | F_act = (A − M' + (6κ−1)z)_+ ≤ 6(w+ℓ) + …; greedy cost κd_z ≤ F_act/6 + κη | OK. Needs only 6κ−1 ≥ 0. |
| 14768 (2.19) | (3.9) ⇒ A ≤ M | OK for κ ≥ 1/6 |
| 14824–14827 (2.1j) | 0 ≤ 6κ−1 ≤ 5 | OK for κ ∈ [1/6,1] |
| 14885–14889 (2.1i) | Choices of ξ, η, ε0, "uniform for κ ∈ [3/4,1]" | OK. The choices do not involve κ; ξ < ρ/30 uses the margin M/15. |
| 14905–14907 | Principal exponent κz, contour cost 2ez | OK |

- **General threshold.** In the comparison, write q = 6κ−1. The bounds become A_comp ≤ (2/3 + 1/(3q))M + ξ and A_comp + qz ≤ (5/6 + 1/(3q))M + ξ, so the margin is (q−2)M/(6q). It is positive for every κ > 1/2 and equals M/15 at q = 10/3.
- **Why 13/18.** It is exactly the point where OpenAI's printed constants (M/20, 23/30, 14/15, M/15) survive unchanged.
- **No use below 11/15.** Going lower is useless: Liu 2610.26835 (l. 10205–10213) shows B ≥ 13/15 for this method, hence κ ≥ 11/15 > 13/18.
- **§19 uses κ, but not as part of the lemma.** OpenAI's row count (κ = 3/4 + 2Δ, l. 15100; plain capacity ≤ 2/27, l. 15382–15384) uses κ ≥ 3/4 only to check that the selection fits the prime supply. The refinements redo this count themselves.

## 2. Lean port
- **The induction.** The plain moment is the Energy/ width induction:
  - CertifiedBand: Energy/CappedWidthInduction.lean:51;
  - certified_bands and terminal_certificate: Energy/CertifiedExistence.lean:21, 42.
- **Hypotheses:** hκ0 : 3/4 ≤ κ, the zero-free hypothesis 2β−1 ≤ κ, and 51/100 ≤ β.
- **The mesh.** It is reserve/(κ+1) (Energy/WidthSchedule.lean:19), so it is uniform on bounded κ.
- **Instantiation.** κ = kappaPlain D = 3/4 + 2Δ with Δ = β_* − 7/8 > 0. See Moments/DetectorPlainMomentParameters.lean:24, Detector/FinalAssemblyCertifiedBands.lean:32 and Detector/FinalAssemblyUnconditional.lean:20.
- **Consumption of 3/4 ≤ κ.** Found by a declaration scan (every lemma with the hypothesis, plus the lemmas it is passed to). Only these use it substantively:
  - Energy/ReferenceLowBranchGeometry.lean:24, slot_budgets: 7/2·ell ≤ (6κ−1)ell ⇒ ell ≤ M/21. reflected_low (l. 30) then gives 16M/21 + ξ and 13M/14 + ξ.
  - Energy/ReferenceState.lean:77, balanced_reference_geometry: ell ≤ M/21, which gives M/2 ≤ log(X1X2) and the comparison scales. balanced_reflection_margins (l. 98) then gives 23M/30 and 14M/15.
  - Moments/ReflectionRetainedLength.lean:63, positive_slot_width_drop: z ≤ M/21 ⇒ 23M/30 + ξ and 14M/15 + ξ. This is the paper's (comparison-uncentered-margin).
  - Energy/LiveClippingDefect.lean:188, robust_slot_budgets: ell ≤ M/21 + 2e/7. It feeds robust_reflected_low (l. 210: 16M/21 + ξ + 4e/7 and 13M/14 + ξ + 11e/7) and the reserve ξ + 11e/7 ≤ M/14, with e = ξ = M0/56 (l. 263–270, 350–388).
  - Energy/OriginalHighReflectionGeometry.lean:16, ReferenceLowWindow.lean:30 and the …DeletedGeometry gates use those results and the hypothesis xi ≤ width/14.
- **Every other consumption** is κ ≥ 0, κ ≥ 1/6 or 6κ−1 ≥ 0:
  - FirstSourceParameters, Natural*Stage, Capped*, ZeroGrowthPhysical, PositiveBalancedAdmission and PositiveHighSource.balanced_volume_bound;
  - the Moments/*ChildCapacity lemmas take only 0 ≤ κ.
- **Outside the lemma:** Hecke/DetectorRowCountCrossing.lean:106, plain_capacity_bound ((1−2m)/(6κ) ≤ 2/27). This is OpenAI's own row count, and it fails at κ = 13/18 (1/13 > 2/27). The refinements do not use it.
- **kappa_margins.py** checks the robust gate by an exact rational LP (vertex enumeration):
  - κ = 3/4: 16/21 and 13/14 with no slack, and 131/168 and 53/56 at e = ξ = 1/56;
  - κ = 13/18: 23/30 and 14/15, and 659/840 ≤ 5/6 and 799/840 ≤ 1.
  So a Lean version of the extension needs only constant changes: M/21 → M/20, 16/21 → 23/30, 13/14 → 14/15 and reserve M/14 → M/15 (this needs ξ + 8e/5 ≤ M/15, which holds since 13/280 < 1/15).

## 3. What the refinements need and how they justify it
- **Xu–Sang.**
  - Need: two stages. Stage 0 uses κ0 = 3/4 and is licensed by 7/8; it gives β_* ≤ B0 = 69991/80000. Stage 1 uses κ1 = 29991/40000 with (1+κ1)/2 = B0 (Prop. "parameters", l. 1624–1660; admissibility at l. 1366). So the extension is needed at the single value κ1 = 0.749775.
  - Justification: Lemma kappa-extension (l. 485–500) and App. A (l. 1867–1971). They check the prime estimate, the floor (3M/13), the comparison (10/3, M/20, 47M/60, 23M/30, 14M/15, reserve M/15), the child repair (F/6 + η) and the support error (≤ 5). The claim that the transforms and exceptional rows do not depend on κ is verified in §1 above.
  - Their row count (l. 915–1077) uses κ ≥ 13/18 only in "selected plain length ≤ 1/(18κ) ≤ 1/13". The claim "c ≤ 6/13 ⇒ 3/5 ≤ r_* ≤ 1, 1/3 ≤ t − r_* ≤ 1/2" is true; it in fact holds for c ≤ 1. The numerator (l. 1079–1260) applies the moment with 2m + 6κz_p ≤ 1. Both are consistent.
- **Liu 2610.12234.**
  - Justification: §28 (l. 4333–4621), Prop. "iv-plain", proved for [13/18,1]. It contains:
    - the prime bound with s_κ ≥ 31/36;
    - (iv-width-floor) z ≤ 3M/13, κz ≤ M/6;
    - 6κ−1 ≥ 10/3, z < M/20, A − z > 47M/60;
    - (iv-comparison-margins) 23M/30 and 14M/15 with ξ_ref ≤ ρ/30;
    - (iv-greedy-child);
    - 0 < 6κ−1 ≤ 5.
  - Its "No assertion is made below κ = 13/18" (l. 4621) is consistent.
  - Its count, Prop. "iv-count" (l. 4650–4745), uses z_S ≤ c/12 < 1/5 with c = 2/(3κ).
  - Need: κ ≈ 2·0.874957 − 1 ≈ 0.74991, and a rational triple with 13/18 < κ < 3/4 (l. 5905).
- **Liu 2610.26835.**
  - Justification: Part II reproves the moment over a general F (Thm v:thm:M2-assembled, l. 6971–6981: "κ ≥ 13/18 ensures 6κ−1 ≥ 10/3"). The comparison at l. 5814–5832 uses z ≤ (M−A)/q_κ < M/20 and A_comp + q_κz ≤ (5/6 + 1/(3q_κ))M + ξ ≤ 14M/15 + ξ. The terminal bound is κz ≤ M/6 (l. 7064).
  - Need: κ♭ = 2B♭ − 1 = 0.749773, reached by iteration from 3/4.
- **Verdict.** All three are proofs by modification of the source. Their list of κ-dependent steps matches every κ-use found above, in both the TeX and the Lean. None is machine-checked.

## 4. Not checked here
- OpenAI's proof for κ ≥ 3/4 itself: it is Lean-checked through our port (Riemann/Deep/QRH/REPORT.md), and was taken as given.
- The refinements' other steps: the averaged numerator beyond its use of the moment, their parameter certificates, and Liu's general-field Part II beyond the moment's κ-steps.

## Cost
- kappa_margins.py: 0.5 s and 16 MB peak (pure stdlib). Everything else was rg and reading.
- No Lean builds; no commits.
