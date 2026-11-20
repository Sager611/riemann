INCOMPLETE: every row named in the brief, and every body named in the follow-up, is now stated one level down with inputs marked KNOWN/NEW/PLUMBING, and the sources are located: the multiplicity m ≍ √x is the explicit rank choice RankVal, the inert gain is proved for every D7 datum, and the m-step kernel chain is asserted in laws_true_fam. 8 finite probes pass with no discrepancy, but the NEW rows that live on D7/TC data have no finite instance I could build, so they are unprobed (reviewer, logs/reviewer.md part 2).

# AlpRows: the HypsAlt theta-lattice rows on the closed route of Alpöge's Lean RH archive (alprows, 20 Nov 2026, 20:32–21:30 UTC)

Archive .scratch/alp/lean/src/ (READ-ONLY; no Lean, no lake). Paths relative to it; TPS = Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries. Details: ROWS.md (§0–§4 first pass, §5 one level down), PROBES.md (P1–P8). dist = alpdep name-BFS cone distance; −1 can be a BFS miss (§5g).

## Short answer
- **Plumbing around NEW rows.** HypsAlt (TPS/Hyps.lean:22-90) is assembled by nonempty_hyps_of_allRows_alt from allRows_of_open0_alt (TPS/LinkROfIsLinked.lean:50-86). That theorem receives the theta-lattice ties hPA, hLamP, hDivG, hDcPtG and selJ verbatim from OpenParameterStatements (GeometricWeakTies.lean:39-51). Those rows are proved in open_orbit_of_unramified_exponents_alt: hDivG ← bk_law_entries_alt → at_hDivG_main_alt, from published representation-theory inputs plus the BkLaw of the built points (ROWS §5a).
- **The multiplicity m ≍ √x is a parameter choice (§5h).** FamP.sizes is filled by szRow_numAt from szRow_sizesAt / sizes_sizeN_sizeM_…, which compute it from Val.ranks = RankVal. The window m ≥ x^{1/2−ε} in ChainTied's RestSize (built by linkCore, LossesSmallShares.lean:281-289) comes from WinF ← size_windows_of_thick_density ← Fm.sizes, a field of the family structure FamP (ParConsts.lean:44; SizeM, ElevenLe.lean:27-30: m within x^{±szTh} of √(x·δ1)). No theorem proves it for lattices. Families with these sizes come from classSel / member_bundles_large_unramified_exponents_alt (§5c).
- **The inert-prime gain is a theorem, not a hypothesis.** LinkG is first asserted for a family by linkG_fam_noKer, through link_of_local_closed. Its local input loc_proved (AdmissibleFamilyNoKernel.lean:284-286) proves SecondTheoremReduced D for every D7 datum, so the gain follows from the D7 lattice definitions. Whether D7 data are actual theta-lift lattices is carried by the ties (§5f).
- **The m-step theta-kernel chain is on the cone.** For built data it is first asserted in laws_true_fam (FivePointLawsTrueChains.lean:46-63, dist 17), from I.Ok ∧ I.JetsIn. Its ingredients:
  - ChTh.Ok of the theta chain of length m−2, via ChainTheta.chTh_ok_of_colIn' (dist 18);
  - kernel bounds at conjugate (inert) and level prime places, from JetMoments (NEW) and InducedSimpleLowestAlcove (KNOWN, Jantzen II.5.6).
  It feeds the row aP10 (Conclusion5Sharp). This corrects AlpDep, which put the chain-of-kernels theorem off the cone: only the twin Formalism.L2.Ker.chTh_ok is unused (§5g).
- **The cone misses dot-notation calls.** FamGLmin.link8, linkGLOBAL7_min, linkGLOBAL7_items and linkCore have dist −1, but they are on the route (LinkGl.lean:30 calls hGLm.link8). So "dist −1" is not proof of non-use.

## 1. Rows and marks
| Row | Mark | Probe |
|---|---|---|
| HypsAlt assembly; allRows_of_open0_alt; LinkEta, LinkG producers; linkDIVt_of_rowsT / closed_linkDIVt…; FamGLmin.link8 → link7 → linkGLOBAL7_min | PLUMBING | — |
| hP SiegelMassRatioCM; hG; hU; LevelIndexAtPlace; InducedSimpleLowestAlcove; published inputs of allRows/open_orbit | KNOWN | P1–P3, P5 |
| hn4 ClassNumberFormula₂ | KNOWN in kind, NEW as stated | P5, P6 |
| LinkDIV_at bookkeeping; LinkG_at gain row; Win_at windows | NEW | P4, P7, P8 PASS |
| OpenParameterStatements rows (DivRows3G ties, …), DefiniteSideCondition/RowsU, FamGLmin fields (slope inequalities, gapped multiplicity, integrality), linkCore, SecondTheoremReduced/Loc, JetMoments, laws_true_fam, LinkS, LinkIRR, SmallLevelRankEq, DesignExists | NEW | none: no finite instance without a D7/TC datum |

## 2. Probes (PROBES.md; exact arithmetic; memcap peaks ≤ 47 MB)
- P1–P3: hP's right side is exactly |L_F(1−M,χ^M)|/2^d, the Gan–Hanke–Yu ratio (30/30 cases, 5 fields, M = 3..8). Over Z[ζ12] the class number is one up to rank 3 and R(3) = 1/36 = #Aut(I_2)/#Aut(I_3). IsCMUnram excludes every imaginary quadratic field.
- P4: in DIV, the minimizer is m = x^{1/2}.
- P5: the index factors of ClassNumberFormula₂, by brute force.
- P6: the ClassNumberFormula₂ right side is an integer in 1548/1548 cases. Open point: the GL_k index at primes split in F but inert in E/F.
- P7: in the gain row, demand/cap ≤ 1/3. With the theta slope as gain, the partner weights must supply a/(p+1).
- P8: the Win_at windows are satisfiable from x ≈ 10^{2061–6213} and empty below.

## 3. Review (logs/reviewer.md)
- First pass: INCOMPLETE because of the unopened bodies; this continuation opens them.
- Second pass, on the continuation claims R1–R5: AGREE on R2 (gain derived), R3 (chain on the cone) and R4 (BFS misses dot calls). CANNOT JUDGE on R1 (sizes provenance), fixed by §5h. DISAGREE on R5 (ROWS MAPPED) because the NEW D7/TC rows are unprobed, so the verdict is INCOMPLETE.
- Fixes applied: smallest IsCMUnram fields corrected; P1 rerun with 5 fields; P2 claim narrowed; P4 now says it checks exponents only.

## 4. Limits
- Not opened (two levels below the named bodies):
  - hDivG_linked_of_mIstar2;
  - member_bundles_large_unramified_exponents_alt and classSel bodies;
  - linkG_fam_noKer body and stmtE_of_pl/restChain;
  - linkDIV_at_words;
  - the supplier of I.Ok ∧ I.JetsIn for laws_true_fam;
  - the level structure of cls 0 (P6 open point).
- NEW rows without finite instances are stated, not tested.

## Files
- REPORT.md, ROWS.md (§0–§5), PROBES.md (P1–P8).
- code/: lookup.py, p1_mass.py (P1–P3), p2_rows.py (P4–P6), p7_gain.py (P7), p8_win.py (P8).
- logs/: lookup1-18.log, p1_mass.log, p2_rows.log, p7_gain.log, p8_win.log, path_laws_true_fam.log, reviewer.md.

