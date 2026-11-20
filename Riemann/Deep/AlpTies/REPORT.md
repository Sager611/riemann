INCOMPLETE: the theta-lattice ties are proved at good built points that carry actual theta-kernel data, and no quantifier at any screened point is vacuous, but the m-step chain laws_true_fam bounds its kernels only at the split level prime T_u, not at inert primes, and three things stay unopened: whether the kernel is nonzero, the existence theorem for good points, and the supply of 'admitted' (where the inert gain is decided).

# AlpTies: the theta-lattice ties of Alpöge's Lean archive (alpties, 20 Nov 2026, 21:26–22:05 UTC)

Archive .scratch/alp/lean/src/ (READ-ONLY; read only, no Lean, no lake). Paths are relative to it. GI = Development/AbelianVarietyTheory/Proofs/GlobalInterpolation; TPS = Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries. Details: TIES.md §1–§11, PROBES.md A–E. dist = alpdep name-BFS cone distance (code/lk.py, logs/lk1-17.log).

## Short answer
- **The ties are proved, at points built from theta data.**
  - hPA, hLamP, hDivG, hDcPtG and selJ are rows on the slot of good built points. They are proved in open_orbit_of_unramified_exponents_alt (TPS/OpenOrbitOfUnramifiedExponents.lean:39-69) from named published inputs, the BkLaw field and MIstar2.
  - A good built point is by definition the position of a ColW datum (ColIn + Ok + Laws) together with a Built structure (TIES §3, §4).
  - hDivG = DivRows3G (Formalism/LocalToGlobal/Divisibility/TupUsed.lean:112-126) ties the D7 lattices and periods to the indefinite realisation: HLatTie of RR, θ-periods of lift i+1 = jW(θ_{i,t}·lPrj·x), local indices and η-exponents, Bk = column form.
- **No vacuity found at the screened points (§3).** TC data exist (InEx.lg_colIn, dist 16; the author's model step_colW). Every chain lattice has full rank (TC.Laws.lat). m ≥ 2 always, and m ≍ √x on the route (a parameter, alprows §5h). The Used0 guards are removed by TupAllUsed. The ⊥/0 arguments of colsBk only fix a type index. BkLaw is a non-vacuous identification (the archive's StateChecks satisfy it and break it). The families over good points are built by classSel (dist 8).
- **Correction: the m-step chain is a level-prime statement.**
  - KerBar and KerP hold at places w over T = B.T u (OverT w := resChar w = T).
  - Built.Split.nonconj forces P ≠ P̄ for every P ∈ ST, so T_u splits in E/E⁺ (probe C).
  - KerBar ('conjugatePrimePlaces') is the T^{f⌊(t+1)/2⌋}-divisibility of each kernel step at P̄; KerP is integrality at P up to T^{f N_i + Σ sDiv}; MomAt is likewise at T_u.
  - So laws_true_fam states nothing at inert primes. Alprows' gloss "conjugate (inert)" misreads it (TIES §6).
- **The real content of the kernel conditions sits in KerEq.kmap (probe D).** The zero kernel satisfies the other five conditions. kmap forces k to realise the chain's Kmap, which is built from the jets at 0 of the theta function thAt at pairs of lattice points (TIES §8). Whether thAt ≠ 0 at built data was not checked.
- **The inert-prime gain is decided outside the ties.** It is LinkG, derived by loc_proved from 'admitted D', a 40-clause conjunction (TableLineDegreeFive.lean:234-240; it includes FrobeniusExpansionCongruence, Herr and Fontaine clauses and ZeroGainAbove). Its supplier is traced in the follow-up line below. The divisibility itself is LinkG_at's LocLe: NuLat ⊂ P^{⌊eτ_w⌋}·MramL at inert w (DataDIVDefs.lean:12-55).
- **Follow-up (TIES §11).** admitted D is proved for every D7 by second_theorem_closed_everywhere (PlaneCountTheoremWithConstants.lean:1121-1135, dist 13), from Loc6 (Weyl dimension, branching and Pieri rules: KNOWN; Model.*/Joints.* lemmas: NEW) and the extremal projector. The inert-place content is one D-independent local theorem, aux₃ (TableLineDegreeFive.lean:30, which ignores D: an inequality over unramified quadratic DVR settings with p ≥ 13), and is the next falsification target. No lemma forces thAt or Kmap to be nonzero; the only nonvanishing condition is ChIdx.Ok.nonvan, local factors for steps i < m − 2. A zero kernel would make only Conclusion5Sharp (aP10, a denominator cost bound) trivially true, which is harmless. So the verdict is unchanged.
- **Marks.**
  - KNOWN: InducedSimpleLowestAlcove = Jantzen, RAGS II.5.6 (probe A: 0 counterexamples in 111 cases; Sym^p of GL_2 is never simple).
  - NEW (in-archive, no literature match found here, no adversarial novelty check run): JetMoments/MomAt (classical analogue: Katz's p-adic CM periods at ordinary primes), the KerBar/KerP/KerUK kernel bounds, the five ties and their proofs.
  - Received identification: BkLaw. PLUMBING: Sees, the slot/row assembly.

## 1. The objects in standard notation (TIES §1, §2, §5, §8)
- CM field E/E⁺, word u, ranks N_i = k + i(t+1), i = 0..m.
- Chain of steps S_i: decorated lattice classes (deco.L ⊂ E^{N_i}, full lattices), finite stabilisers, unitary groups via gam.
- V_i = ∇(λ_i), the polynomial model over LL of the weight λ_i = lamTrue u i.
- Kernel step i → i+1: k = kAt = KerT.ker (Definitions/GlobalTwisting/BelowFiniteNonempty.lean:20-23), a theta-series kernel on Hom(E^{N_i}, E^{N_{i+1}}) with test function Ξ. Its complex shadow is Kmap = (unit ρ-factors)·(divDe ∘ jets at 0 of thAt ∘ sharp).
- The dual pair is the unitary pair (U(V_i), U(V_{i+1})); the Weil representation is split by psiWeil P ξ0 at P ∈ ST (LatAtLattice.lean:75-80). The model is polynomial/Fock (jets at 0; the archive's seesaw is a Fock-sum identity, AlpCore CORE §6).
- Places:
  - T_u: a level prime, split in E/E⁺, with ST = vT u (one prime of each conjugate pair).
  - w ↔ P (seenP = true) carries KerP; w ↔ P̄ carries KerBar.
  - The inert primes of ASK enter only through NuDef/LinkG.
- D7 lattices: D.div (RR, Nu u, Bk u) are data of the later stage Z.V (DivStage, Definitions/AutomorphicMembership/NuInRange.lean:10-30), tied to the chain by the isometric embedding nuIn u with range ramSub of lat 0. Which lattice Nu u is, and whether OneMapIntoCC holds, are fixed by the rows.

## 2. Instantiation trace on the route (TIES §3, §4, §7c, §9)
- allRows_of_open0_alt (TPS/LinkROfIsLinked.lean:50) → six_rows_at_member_alt (TPS/SixRowsAtMember.lean:11-24): rows on the slot sAnd (IsGoodE 𝒮) True, with hC := goodTwistCores_over_trueColumns.
- → row_level_valuations (GI/LevelPrimeColumnValuationsBundlePointsOverColumn.lean:11-19) → closing_true_cols w.I w.ok hv.jets … → … → laws_true_fam (dist 17).
- So the TC datum is w.I of the ColW attached to a good point. I.Ok = w.ok. I.JetsIn comes from the Built structure (colView_of_built, GI/ColViewOfBuilt.lean:13; Built fields at Definitions/AutomorphicMembership/NonzeroPlacesOverPrimes.lean:147-170).
- The built datum is dat I l Z u a := Inh.ofOk … (GI/MoveStep.lean:26-28); its kernel is kAt … (:44-45); Split comes from SplitPts.exists_split (:61-67).
- The route proves that the slot is inhabited (classSel, CharacterBoxes/ClassSel.lean:14-29 → memberTop_final_contextParameters), non-constructively; no explicit instance is exhibited. ExBL only gives bundles at points that already exist (TIES §9).

## 3. Non-vacuity screen
| Screen | Finding | Evidence |
|---|---|---|
| TC hypothesis Ok | satisfiable; Laws, Law0 derived from Ok | InEx.lg_colIn (ExistsColumnInputAssembled.lean:328); Model.step_colW (NotConditionOfFalse.lean:25-33) |
| JetsIn | rationality of normalised kernel jets (ColIn.JetsIn, JetsInIffStepField.lean:109-125), read off Built | colView_of_built |
| lattices | full rank (Submodule.IsLattice) | TC.Laws.lat (LatAtLattice.lean:124) |
| m | two_le_m; i < m ranges over ≥ 2 values; m ≍ √x on the route | Ranks; alprows §5h |
| Used0 guards | removed under TupAllUsed | TupUsed.lean:73-110 |
| colsBk placeholders | type index only | colData2 … B := C.B (DivisionCompatibleToSub2Indef2.lean:92-100) |
| BkLaw | non-vacuous identification | StateChecks.exists_bkLaw / exists_not_bkLaw (NonzeroPlacesOverPrimesOfLifts.lean:104-113) |
| kernel conditions | k = 0 passes five of six; KerEq.kmap carries the content | PROBES D; DataKernels.lean:63-65 |
| kernel ≠ 0, Nu u ≠ 0 | NOT checked | — |

## 4. One level down (TIES §7)
- open_orbit_of_unramified_exponents_alt (dist 7):
  - hPA, hDcPtG ← linkedAt_*_prmR (dist 8; no published input): NEW.
  - tie, selJ, hMi, hLamP ← theorems with published inputs (Howe unramified, Siegel–Weil, doubling, Harish-Chandra, Casselman–Wallach, Zelevinsky, Flath, Dixmier–Malliavin; in-archive statements not compared): NEW with KNOWN-named inputs.
  - hDivG ← bk_law_entries_alt → at_hDivG_main_alt → hDivG_linked_of_mIstar2 → divisionRows_pairing_vanishes_of_rational_pairing_law (GI/BkLawEntries.lean:40-48, dist 12) + pairing_definite_base_vanish_of_miStar2 (dist 12): NEW.
- laws_true_fam (dist 17):
  - hyp_of_colIn (18): extraction.
  - equivariance laws, normalised theta family, chTh_ok_of_colIn', trueChainKernel_eq_pointKernel_at_representatives, kAt_uk (18): NEW.
  - KerBar/KerP bounds (18) from jetMoments (19; := jetMoments_of_unitPeriod, 20): NEW.
  - weylModule_simple_closedLowestAlcove (18) = InducedSimpleLowestAlcove: KNOWN (Jantzen II.5.6).
- Nothing STRONGER than the literature was identified; the NEW statements have no literature counterpart that I could match.

## 5. Probes (PROBES.md)
- A PASS: the alcove condition is the right one. Sufficient: 0 counterexamples in 111 cases. Sharp for GL_2 at n = p. Not necessary: 16 cases are simple outside the alcove.
- B toy: CM jets at split p are integral up to p^{O(log a)} with no growth; at inert p they gain Σ a/p^{odd}. This is consistent with MomAt being an integrality statement at split T_u.
- C PASS: the Split axioms force P ≠ P̄. Two methods agree in 134/134 cases. ST is nonempty only for T split in E/E⁺.
- D: the zero-kernel screen (logic).
- E: no finite theta-lift instance of DivRows3G fits; BkLaw is checked logically and via the StateChecks.

## 6. Review (logs/reviewer.md)
- The independent reviewer AGREES with INCOMPLETE.
- It DISAGREED with two claims as first worded:
  - 'no vacuity' (C2): ColW nonemptiness alone does not exhibit a joint instance;
  - 'BkLaw = plumbing' (C3).
- AGREE on C4 (Jantzen). CANNOT JUDGE on C1 in full and on C5.
- Fixes applied: C1 reworded to 'stated only at level-T places'; C2 reworded to 'no degeneracy found', with the classSel existence and ExB semantics added; C3 replaced by the StateChecks evidence (non-vacuous identification); NEW qualified as 'in-archive, no novelty check'.

## 7. Limits (open, in order of decisiveness)
1. Nonvanishing of thAt jets (Kmap ≠ 0) and dim Nu u = dim ramSub > 0 at built data. Reduced in TIES §10: Kmap ≠ 0 follows if the test function xiAt (an XiEntry/XiCpl1/xiLast product of the decoration indicator phiT and local factors) is nonzero and the coset-sum theta xiV does not annihilate it.
2. The proof of the uniform local theorem aux₃ (second_theorem_closed_everywhere, TIES §11a), and NuDef/NprimeDef as stage constraints.
3. The bodies of classSel/memberTop_final_contextParameters and compatibleSystem_weilBlock_ne_alt (existence).
4. JetMoments' proof (unitPeriodTranslates_moments_finite, MomEq, MomLat) and a faithful finite MomAt instance.
5. The in-archive 'Published' statements versus the literature.

## Files
- REPORT.md, TIES.md (§1–§11), PROBES.md (A–E).
- code/: lk.py (multi-regex name/cone lookup), p_alcove.py (A), p_jets.py (B), p_split.py (C).
- logs/: lk1-17.log, p_alcove.log, p_jets.log, p_split.log, reviewer.md.

