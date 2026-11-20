# TIES — the theta-lattice ties of Alpöge's archive, unfolded (alpties, 20 Nov 2026)

Archive .scratch/alp/lean/src/ (READ-ONLY). Paths relative to it. GI = Development/AbelianVarietyTheory/Proofs/GlobalInterpolation; TPS = Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries. dist = alpdep name-BFS cone distance (code/lk.py, logs/lk*.log).

## 1. laws_true_fam (GI/FivePointLawsTrueChains.lean:46-63, dist 17)
- Statement (:46-50): ∀ B : Mem.Base, X : PreX B.toPre, I : TC B X, I.Ok → I.JetsIn → ∀ l : I.Laws, l0 : I.Law0 l, Z : TC.Later B (I.cols l), u : B.Word, a : I.ι, ∃ Bd : Built.Dat (I.rec7 l Z) u with Bd.ch = tw (I.tw l l0 Z) u a, ∀ i ≤ R.m Nonempty (Split (Bd.J.A i) Bd.J.T), Bd.Sees, ∀ i < R.m ∃ k : Bd.Ker i with KerEq ∧ KerDil ∧ KerUK ∧ KerBar ∧ KerP, and Bd = dat I l Z u a.
- Witness (:58): Bd := dat I l Z u a (explicit construction; the last conjunct pins it), split u a h i, sees l Z u a; kernel = the placeholder of kerEq (:59), KerUK ← kernelCondition_at_row_family (:36-44: kAt … (I.idx a) I.A u i).
- Chain: hth (:54-56) = ChTh.Ok of chTh … I.A u (B.R.m − 2) for every a, u, via ChainTheta.chTh_ok_of_colIn'. KerBar ← trueChainKernels_bound_conjugatePrimePlaces kernelJets_levelPlace_integralMoments weylModule_simple_closedLowestAlcove …; KerP ← trueChainKernels_bound_levelPrimePlaces (same two inputs) (:60-62).
- Inputs: hyp_of_colIn (hH.ana gives hcc hδ hχ), normalisedFamily_satisfies_equivarianceLaws (hE), normalisedThetaFamily_vectors_lattices_thetaZero (.1 = hV), trueChainKernel_eq_pointKernel_at_representatives (KerEq).

## 2. The built datum Bd : Built.Dat D u (Definitions/GlobalTwisting/KerIsTerm.lean, dist 17)
- Dat (:42-58) over a D7 datum D and word u: a number field E (StarRing), embeddings emb : Tau → (E →+* LL), I : ChainIn E LL Tau D.R (the chain of steps I.S i, i = 0..m), J : AtIn I (lattice data J.A i : TDat, integer J.T, finite prime set J.ST), Q : PerIn I (periods), Z : FmIn J places, proofs hI hJ hQ, laws l, lf, hT.
- Bd.ch : Ch D u (:62) = the chain built by forms B.J B.Q … (the forms of the chain); KerEq.kmap ties a kernel to Bd.ch.C.Kmap.
- Split (A : TDat S ST) (T : ℕ) (:23-37): a setoid 'off' on each data set S.D b, invariant under S.G b, with SameOff ST (lattices agree off ST ∪ conj ST), injectivity of the decoration on classes, surjectivity onto Ok decorations, up/dn/dia moves preserve classes, gam_surj (every unitary γ ∈ GL_n(E) is hit), T_mem (T ∈ P), T_nsq (T ∉ P²), nonconj (Q ≠ P̄ for P, Q ∈ ST).
- nonconj with Q = P forces P ≠ P̄ for every P ∈ ST: the level primes in ST are SPLIT in E/E⁺ (P ≠ P̄), never inert.
- Sees (:82-84): for every place w over T (OverT) and P ∈ ST, emb(tauOf w P)(x or x̄ by seenP w) ∈ O_w ⇔ v_P(x) ≥ 0: w induces P (seenP w = true) or P̄ (false).
- Ker i (:86-87) := ∀ b b', D_i b → D_{i+1} b' → (V_i →ₗ[LL] V_{i+1}), V_i = VV LL (lam i) (the algebraic model of the weight-λ_i representation over LL).
- KerEq k (:91-96): equiv (G_i × G_{i+1}-equivariance through rho) and kmap: Bd.ch.C.Kmap i c c' (ev_i F) = ev_{i+1}(k … F) (k realises the chain's kernel map on evaluations).
- KerDil (:98-99): k invariant under the 'dia' moves (diamond operators at P ∈ ST).
- KerUK (:101-103): Σ_{d₁ above d' at P} k(d, d₁) F = T^{⌊(t+1)/2⌋} Σ_{d₀ above d at P} k(d₀, d') F (a Hecke-type U_P relation).
- KerBar (:105-108): at w over T with seenP w = false (w ↔ P̄): k maps valAt w i (Lsh_i) into T^{f⌊(t+1)/2⌋}·valAt w (i+1) (Lsh_{i+1}).
- KerP (:110-113): at w over T with seenP w = true (w ↔ P): T^{f N_i + Σ_τ sDiv i τ}·k maps valAt w i (Lsh_i) into valAt w (i+1) (deco_{i+1}).L.
- Vacuity screen: k = 0 satisfies equiv, KerDil, KerUK, KerBar (F' = 0), KerP (0 ∈ L). Only KerEq.kmap pins k: it forces Kmap i c c' (ev_i F) = ev_{i+1}(k F). So the kernel conditions have content iff Bd.ch.C.Kmap is the actual theta kernel and nonzero (§3).
- Places: KerBar and KerP live at the places over T (the split level primes of ST and their conjugates), not at the inert primes of ASK.


## 3. Instantiation trace: where the TC datum I of laws_true_fam comes from (route, top-down)
- allRows_of_open0_alt (TPS/LinkROfIsLinked.lean:50) → period_entries_tie_and_level_valuations_at_built (TPS/SixRowsAtMemberStatement.lean:12) := six_rows_at_member_alt (TPS/SixRowsAtMember.lean:11-24): each of the six rows holds as Lg1 KPx (IsBuiltE 𝒮) (IsGoodE 𝒮) (row at slot S), witnessed with S = sAnd (IsGoodE 𝒮) True (hT := largeBaseProperty_true, :20).
- The only bridge to theta data is hC : OverCols (sAnd (IsGoodE 𝒮) True) := (goodTwistCores_over_trueColumns 𝒮 True).2 (:18-19): every good built point has a column view.
- OverCols Q := ∀ K X f T c, Q … → HasCols … (Definitions/GlobalTwisting/AndLeft.lean:24); HasCols := positionOf (ColView w Z c.O.tw c.O.nch) (:22); ColView (:15-20) = jets : ColIn.JetsIn w.I, rowInt : RowInt w.I w.l Z, tw_eq : tw = w.tw Z, tup : ∃ TupTied … with tp.xi = w.I.xi.
- row_level_valuations (GI/LevelPrimeColumnValuationsBundlePointsOverColumn.lean:11-19): for Q with OverCols Q, via property_at_constructedPoints_of_trueColumnsView (GI/PropertyAtConstructedPointsOfTrueColumnsViewStatement.lean:8; := Fc.Cv.view_at_built), at every point β with β.Is HasCols the property Conclusion5Sharp D u ∧ TwistedColumnBoundAtT is obtained from closing_true_cols w.I w.ok hv.jets hs w.l Z w.l0 O.X u (:18).
- So the TC datum of laws_true_fam is w.I, the column datum of the ColW structure w attached to a built point β; I.Ok = w.ok (a field of ColW), I.JetsIn = hv.jets (a field of the view).


## 4. Existence of TC data (non-vacuity of I.Ok)
- step_colW (Definitions/AutomorphicMembership/NotConditionOfFalse.lean:25-33): from SiegelMassRatioCM (hP, KNOWN), for large bases with KPx, Nonempty (ColW B f T); the ColIn I and I.Ok come from InEx.lg_colIn K (:27, :32), and TC.Laws and Law0 are DERIVED from I.Ok alone (TrueColsTest.laws a1, TrueColsTest.law0 a1, :33). Laws.lat (LatAtLattice.lean:124) makes every chain lattice deco.L a full lattice (Submodule.IsLattice) of E^{N_i}: nonzero when N_i ≥ 1.
- exists_colW, not_step_false, not_stepW_false (:74-89) prove that the Step conditions are not False given HostExQ, LevelPrime, SiegelMassRatioCM, RayField, PrincipalIdealSt, CMCurveGoodOutsidePrimeSt: the author's own non-vacuity certificate for ColW (dist −1: off the name cone, but its premise step_colW is the producer of the ColW data).
- Good built points are by definition positions of ColW data: IsGoodE 𝒮 := positionOf (GoodE 𝒮), GoodE = BuiltE ∧ GoodO, BuiltE = ∃ b : Built 𝒮 K w Z c, EntC b.ker (ThetaChr.lean:80-91); positionOf P := ∃ w : ColW, Z, c' with P and the point = ⟨w.C, Z, c'⟩ (PositionOfIntro.lean:61-65).
- ColIn (TrueColumnChoicesHeckeCharReferenceEqOne.lean:165-201) = choices over a Base B: frame C, Hecke characters Ch (B.R.m of them), prime selection sel, finite prime set S, finite index type ι of character tuples tup (injective, tup mu0 = 1), Hecke characters xi, local constants cst, chain indices idx, ambient J. ColIn.Ok (:219-251): DiagSpaces, Ch.Ok, Ch.Tied (to μ0, η), Ch.WeilTied (to the Weil splitting character ξ0), level primes T u ∈ S, ramified primes ∈ S, xi reads tup, xi(muStar) = C.del, sharp cst, idx Ok, J.Ok, duStar.


## 5. The five tie rows of OpenParameterStatements (TPS/GeometricWeakTies.lean:39-51; defs TPS/CountsIrreducibilityConstants.lean:70-93)
- Each row is LinkedAt p r := Lg1 KPx (IsBuiltE p.s) (IsGoodE p.s) (fun S => r S (BundleLinkage p)) (TPS/RowGood.lean:38-39), Lg1 KP Q Q' r := ∃ R, LargeBaseCondition KP Q R ∧ r (sAnd Q' R) (Definitions/GlobalTwisting/AllRows.lean:15-16): the row holds on the slot of good built points (R = True in every use found).
- hPA: Row (AtAll PrimitivityBundleCondition) = PrimRows3 b.B.Q.Rl b.B.Q.F ∧ b.P.A3.fam3.TypesSep3 (:74-75).
- hLamP: RowsHoldLam3 b.P (:77-78). hDcPtG: DcPtG b.P (:83-84). selJ: Nonempty (IndefiniteRealisationJoin b.R b.S b.B.Q) (:86-87).
- hDivG: DivRows3G b.P (:80-81); DivRows3G (Formalism/LocalToGlobal/Divisibility/TupUsed.lean:112-126) =
  - sharpRR, baseRR: HLatTie of the D7 row-side lattice RR (sharp/base) with the indefinite realisation's RR via XE (an isometry tie of hermitian lattices);
  - divisionRows_rationalLift: RationalLiftSpanAtUsedTuples;
  - thTheta (ThThetaG, :61-64): θ-period of lift i+1, thetaP (tup t) (i+1) (jMt … (eMt i t x)) = jW i (θ_{i,t} · lPrj_{i,t} · x), for every i : Fin m and used tuple t;
  - num (NumTieG, :94-105): local index Ij, exponent jw at w, and etaJ (i+1) w of D equal the indefinite side's Iw, jw(HT, Tl), and etaJ (or 0 off the twist set Tw);
  - colsBk (:124-125): D.div.Bk u (eNu x) (eRR y) = (colData2 P.C ⊥ ⊥ 0).B (wU u) x y. Checked: colData2 … .B := C.B (DivisionCompatibleToSub2Indef2.lean:92-100), so the ⊥/0 arguments only fix a type index; not a degeneracy.
- Guards: ThThetaG, DcPtG, NumTieG.Ij/jw hold only for tuples with Used0; thTheta_of_G, dcPt_of_G, numTie_of_G (TupUsed.lean:73-110) remove the guard under TupAllUsed P (on the route from tup_all_used_of_realisation, alprows ROWS §3).
- BkLaw D := RationalMapValues D → D.BuRational (Definitions/AutomorphicMembership/NonzeroPlacesOverPrimesLifts.lean:107), a field bk of the parameter structure St (:109-121). RationalMapValues (Formalism/SharpRestriction/Divisibility/WeightDegreeSizes.lean:1264-1265): every complex pairing value Bu u (jN φ) (jEE x) lies in embLL(LL). BuRational (Formalism/LocalRingCardinality/RhoDetPrimeEqMuTildeMulChiPow.lean:364-365): Bu u (jN φ)(jEE x) = embLL (Bk u φ x). So BkLaw says: if the complex pairing is LL-rational, the formal map Bk is its LL-descent. With Bk defined by descent (BkDef, :1267-1269) it is automatic: PLUMBING (an identification, not a theorem about lattices).
- HasBare (NonzeroPlacesOverPrimesLifts.lean:122-123, dist 14) := ∃ St with IsEmpty Loc and all η j v = 0: a 'bare' parameter with no local data and zero exponents exists as a notion on the cone (use not traced).


## 6. Places: the m-step kernel chain is a split-level-prime statement; the inert gain sits elsewhere
- OverT w := resChar w = X.A.T (Definitions/GlobalTwisting/LevelChainSharpColumnValuationOwnChainIff.lean:79). Built.Split.nonconj (KerIsTerm.lean:36) forces P ≠ P̄ for P ∈ ST. So KerBar (seenP w = false: w ↔ P̄, the 'conjugate prime places') and KerP (w ↔ P, the 'level prime places') live over the split level prime T_u = B.T u (probe C). MomAt likewise assumes [Fact (T u).Prime] and works at seenBarC/seenC places over T u.
- KerBar: per step, k maps the local lattice at P̄ into T^{f⌊(t+1)/2⌋}·(next lattice), with N_{i+1} − N_i = t + 1. KerP: T^{f N_i + Σ_τ sDiv i τ}·k is integral at P. KerUK is a U_P-type relation with factor T^{⌊(t+1)/2⌋}, the δ^{1/2} normalisation of a split-place Hecke action.
- None of these concern the inert primes p of ASK. The inert-prime gain is not in laws_true_fam.
- The inert gain is LinkG (alprows §5f): linkG_fam_noKer (Formalism/LocalToGlobal/GroupLattices/AdmissibleFamilyNoKernel.lean:241-271) via loc_proved (:284-286) = criterion_restricted_of_closed D (stmtE D) for every D7.
  - Correction to alprows: SecondTheoremReduced D (PlaneCountTheoremWithConstants.lean:1020-1025) is 'admitted D … → aux₃ bounds → ConclusionE D'. admitted (TableLineDegreeFive.lean:234-240) is a 40-clause conjunction (adm1…adm26, FrobeniusExpansionCongruence, HerrExtensionClause, HerrCocycleExtensionClass, FontaineFormalProperties, ZeroGainAbove, LevelPrimesProviso, …). So the gain is derived from admitted, not from the D7 definitions alone.
  - The family-level supply of admitted (via SecondTheoremAtFamily ← second_theorem_closed_form ← Loc6 local_model_results, :43-60; and LocalDataSystem's FormalGroupSystem with per-lift lower bounds low u w ≤ Σ_{j=1}^{m} Σ_v Fn v j − Σ_j dG j, :79-91) was not traced.
- Inert-place lattice shape (Formalism/LocalRingCardinality/DataDIVDefs.lean:12-55):
  - MramL D u = (lat 0 of the chain).L pulled back along nuIn;
  - NuDef: NuLat u ⊂ MramL at every w ∤ T, with equality unless w is inert and not dropped, and colength lenMN·e_w;
  - NprimeDef: Nu u = Iu·NuLat with Iu = ∏_w P_w^{−⌊e_w·sg u w⌋};
  - LinkG_at asks LocLe at level ⌊ur(τ_w)⌋, i.e. NuLat ⊂ P^{⌊e τ_w⌋}·MramL locally: the divisibility claim at inert w.
- DivStage (Definitions/AutomorphicMembership/NuInRange.lean:10-30): the D7 lattices D.div (RR, Nu u, Bk) are data of a later stage Z.V. They are tied to the chain only by nuIn u (injective, range = (C.fm u).ramSub, isometric at every place) and by jN, jEE (basis to basis). So which lattice Nu u is, and whether Bk maps it into cc (OneMapIntoCC, DataArith.lean:180), is fixed by the rows, not by the TC datum.


## 7. One level down (task 3). dist from logs/lk*.log. Marks: KNOWN (reference) / NEW / PLUMBING
### 7a. open_orbit_of_unramified_exponents_alt (TPS/OpenOrbitOfUnramifiedExponents.lean:39-69, dist 7)
| Input | file:line, dist | Statement (short) | Mark |
|---|---|---|---|
| hNE | :42-43 | extensionWithList → Nonempty extensionWithLaw, by projection x.1 | PLUMBING |
| hB (BkLaw) | :44, field q.x.nx.2.1.2.down of the point's parameter extension | RationalMapValues → BuRational (§5) | PLUMBING (received identification) |
| hPA | linkedAt_hPA_prmR, TPS/HermWCellLastMulMulVecEqOfInWbStatement.lean:15, 8 | PrimRows3 ∧ TypesSep3 at good points; no published input | NEW |
| hDcPtG | linkedAt_hDcPtG_prmR, TPS/RowEDcPtLinkedStatement.lean:15, 8 | DcPtG: embLL(Dc.Pt) = perCol (period columns) at used tuples | NEW |
| tie | linkedAt_tie_of_nonempty_extensionWithLaw, TPS/AtTieStatement.lean:15, 8 | Row_tie_alt; inputs DiscreteSeriesCover, HarishChandraModules, CasselmanWallach, WeilRepresentationUnitaryPair, DiscreteSeriesExists, HarishChandraFinite, CuspRapidDecay | NEW, KNOWN-named inputs |
| selJ | entries_combined_real_alt, TPS/JoinLinkedRealizationStatement.lean:15, 8 | Nonempty IndefiniteRealisationJoin; same published inputs + hPA hDc hTie | NEW |
| hMi | main_theorem_published_inputs_third_alt_alt, TPS/MainTheoremPublishedInputsThird.lean:34, 8 | MIstar2 rows; inputs incl. howeUnram_holds, WeilSiegelFormulaUnitary, SiegelWeilUnitaryConvergent(2), DoublingMethod, UnramifiedDoublingGL, StrongApproxSU | NEW, KNOWN-named inputs (Howe, Siegel–Weil, Piatetski-Shapiro–Rallis doubling) |
| hDivG | bk_law_entries_alt (TPS/AtHDivGOfMiStatement.lean:20, 8) := at_hDivG_main_alt (TPS/AtHDivGOfMi.lean:66-67, 9) := hDivG_linked_of_mIstar2 (:60-65, 10) | DivRows3G at good points, from hDivG_of_mIstar2 (:49-58): divisionRows_pairing_vanishes_of_rational_pairing_law (GI/BkLawEntries.lean:40-48, 12; inputs hC, Holds, BkLaw, BtVanish) and pairing_definite_base_vanish_of_miStar2 (Development/AbelianVarietyTheory/Proofs/DivisorProofs/PairingDefiniteBaseVanishOfMiStar2.lean:10-17, 12; MIstar2 ⇒ BtVanish) | NEW |
| hLamP | real_group_inputs_combined_alt_alt, TPS/RealGroupInputsCombined.lean:39, 8 | RowsHoldLam3; inputs HarishChandraAnalyticity, DiscreteSeriesClassification, CuspTraceClass, DixmierMalliavin, Flath, ZelevinskyPolynomialRing, GeometricLemma, GreenCuspidal | NEW, KNOWN-named inputs |
- Ties are proved at bundles built from witnesses (BundleWitness q b; the bundle's P equals realPackageOfRealInput w.y.y, AtHDivGOfMi.lean:52-58). The published inputs are Props in GlobalTheory.Published. Their in-archive statements were not compared with the literature (AlpAudit / AlpGal scope).
### 7b. laws_true_fam (GI/FivePointLawsTrueChains.lean:46-63, dist 17)
| Input | file:line, dist | Statement (short) | Mark |
|---|---|---|---|
| hyp_of_colIn | Definitions/GlobalTwisting/JetsAllOfColIn.lean:96-105, 18 | from I.Ok: analytic data (curve ideal cc, discriminant 12th root δ, χ) | PLUMBING-type extraction |
| normalisedFamily_satisfies_equivarianceLaws | Proofs/IwasawaDecomposition/NormalisedFamilySatisfiesEquivarianceLawsStatement.lean:7, 18 | equivariance of the normalised family | NEW |
| normalisedThetaFamily_vectors_lattices_thetaZero | GI/NormalisedThetaFamilyVectorsLatticesThetaZeroStatement.lean:7, 18 | vectors, lattices and θ(0) data of the normalised theta family | NEW |
| ChainTheta.chTh_ok_of_colIn' | GI/ChThOkOfColIn.lean:148-153, 18 | ChTh.Ok of the theta chain of length m − 2 | NEW |
| trueChainKernel_eq_pointKernel_at_representatives | GI/TrueChainKernelEqPointKernelAtRepresentativesStatement.lean:7, 18 | KerEq.kmap: chain Kmap = evaluation of kAt | NEW (carries all kernel content, PROBES D) |
| kerDil, split, sees | GI/MoveStep.lean:53-71 | diamond invariance; Split from SplitPts.exists_split; Sees from Inh.sees_ofOk | NEW / NEW / PLUMBING |
| kernelCondition_at_row_family (KerUK) | GI/FivePointLawsTrueChains.lean:36-44 | U_P relation via UpLaw.kAt_uk, moveStep_fam | NEW |
| trueChainKernels_bound_conjugatePrimePlaces (KerBar) | GI/TrueChainKernelsBoundConjugatePrimePlacesStatement.lean:7, 18 | T^{f⌊(t+1)/2⌋}-divisibility at P̄ over the split level prime | NEW |
| trueChainKernels_bound_levelPrimePlaces (KerP) | GI/TrueChainKernelsBoundLevelPrimePlacesStatement.lean:7, 18 | bounded T-denominators at P | NEW |
| jetMoments (JetMoments/MomAt) | GI/KernelJetsLevelPlaceIntegralMoments.lean:13-17, 19 := jetMoments_of_unitPeriod (GI/JetMomentsOfUnitPeriod.lean:19, 20) | p-adic integrality of kernel jet sums at T_u relative to a unit-period table | NEW as stated; classical analogue Katz (Ann. Math. 104 (1976); Invent. Math. 49 (1978)): p-adic periods at ordinary CM points |
| weylModule_simple_closedLowestAlcove / InducedSimpleLowestAlcove | GI/WeylModuleSimpleClosedLowestAlcoveStatement.lean:7, 18; PublishedTheorems/GlobalTheta/WeylDualWeylLatticesAgreeOfLowestAlcove.lean:19-25, 6 | ∇(λ) simple when ⟨λ+ρ, α∨⟩ ≤ p | KNOWN: Jantzen, RAGS, II.5.6 (PROBES A) |
### 7c. Where the TC data come from (route)
- colView_of_built (GI/ColViewOfBuilt.lean:13, 13) → overCols_built (:20, 13) → overCols_goodE (GI/GoodTwistCoresOverTrueColumns.lean:15, 11): the ColView (JetsIn, RowInt, tw_eq, tup) is read off the Built structure (Definitions/AutomorphicMembership/NonzeroPlacesOverPrimes.lean:147-170: pin, coh, g, col : ColumnDataSplitPrimes, what₂, rx, n, nx, sx, ker0, ker : KerT, div : DivT, th : ThT, obj, ox). PLUMBING.
- ColIn with Ok: InEx.lg_colIn (Proofs/IwasawaDecomposition/ExistsColumnInputAssembled.lean:328, 16). NEW (existence).
- Good built points of the families: member_bundles_large_unramified_exponents_alt (TPS/MemberBundlesLargeUnramifiedExponents.lean:34-63) → compatibleSystem_weilBlock_ne_alt (TPS/ExistsWeilBlockNeAtMember.lean:41, 9) for ExBL; member_exists_baseLattice_alt for MemberE. NEW, not opened.


## 8. What the chain kernel is (KerEq.kmap target)
- Bd.ch.C = ChainIn.chain (Definitions/AutomorphicMembership/PredMemGamRanks.lean:442-463), unf i c c' := fnOf ∘ I.jet i c c' ∘ polOf (:439-440).
- For the built datum, I = TC.chIn u a (LatAtLattice.lean:68-73) with jet i c c' := jetMap (jetC cD (th u a i …)) 0, th := thAt I.Ch I.sel (tup a) (cst a) (idx a) I.A u i b b' p p' : VCx → ℂ (:64-66).
- So Kmap i = (unit-modulus ρ-factors)·(divDe ∘ jets at 0 of the theta function thAt at the pair of points (p, p') ∘ sharp) (DataKernels.lean:63-65).
- k = kAt = KerT.ker … (xiAt …) is its algebraic incarnation (KerEq.kmap). AtIn: T := B.T u, ST := (B.primes X).vT u, psi P := psiWeil P.1 ξ0 (LatAtLattice.lean:75-80): the level prime of word u and the Weil-representation splitting at P.
- Whether thAt ≠ 0 at built data, which is what makes the kernel bounds non-trivial (PROBES D), was not checked. ColumnChain.unf is a free field (DataKernels.lean:13-35), so nonvanishing is not built into the structure.


## 9. After review (logs/reviewer.md): corrections and extra evidence
- BkLaw is not automatic. The archive's own checks: StateChecks.exists_bkLaw (Definitions/AutomorphicMembership/NonzeroPlacesOverPrimesOfLifts.lean:104-108) makes BkLaw hold for every D7 by choosing Bk := the LL-descent of Bu (withBk D b); StateChecks.exists_not_bkLaw (:110-113) breaks it for another choice b when Nu u, RR ≠ 0 and RationalMapValues holds; bkLaw_of_bkDef (:31). So BkLaw is a non-vacuous identification constraint (Bk = LL-descent of the complex period pairing), received as a field of the point's parameters. It is not a lattice theorem; mark: received identification (not 'mere plumbing').
- In hDivG the antecedent RationalMapValues is discharged by divisionPairing_rational_on_rationalPoints' (Proofs/IwasawaDecomposition/DivisionPairingRationalOnRationalPointsStatement.lean:15, dist 13), used at GI/BkLawEntries.lean:47.
- ExBL KP Q Q' Lk := ∃ R, LargeBaseCondition KP Q R ∧ ExB (sAnd Q' R) Lk (Definitions/GlobalTwisting/EventuallyAllScaleContextsAndOfBoth.lean:71). ExB (ValidForConstructedPointClassesEventuallyWorking.lean:259) says that every point in the slot has, for each c, a linked bundle. It does NOT assert that points exist.
- Existence of families over good points: UpperTheory.Top.classSel (Development/AutomorphicArithmetic/CharacterBoxes/ClassSel.lean:14-29, dist 8) builds the class cls (sAnd Q' R) Lk E and its members via memberTop_final_contextParameters from SiegelMassRatioCM, VolTie, and BaseCond ∧ FirstG ∧ ExportF ∧ SecondG. So the route proves (non-constructively; body below memberTop not opened) that the slot used by the rows is inhabited by the families the endgame uses. No explicit instance is exhibited.
- Wording fixes:
  - 'no inert-prime gain in laws_true_fam' means: its kernel conditions are stated only at places over T u (split in E/E⁺). Indirect consequences at inert primes are not excluded by this reading.
  - 'no vacuity' means: no degenerate datum found at the screened points (TC existence, lattice rank, m, Used0 guards, colsBk placeholders, BkLaw). Nondegeneracy (thAt jets ≠ 0, Nu u ≠ 0) is not proved here.
  - NEW = in-archive statement with no literature match found in this audit; no adversarial novelty check was run.


## 10. Nonvanishing of the chain kernel, reduced (late check, 21:59–22:02 UTC; logs/lk15-16.log)
- thAt u i b b' p p' := thOf … N N' (xiAt … i b b' p p') (Definitions/AutomorphicMembership/XiAtPenultimateEqXiCpl1.lean:274-276). thOf is linear in the test function, with thOf … 0 = 0 (:268-271); thAt = 0 for i ≥ m (thAt_of_le, :282-284).
- The test function is xiAt (:157-162):
  - XiEntry for i + 2 < m;
  - XiCpl1 for i + 2 = m;
  - xiLast for i + 1 = m;
  - 0 beyond m.
  XiEntry (Formalism/LocalToGlobal/Columns/FiniteSumNonzeroLargeDepth.lean:1116-1119) = phiT(d, d', x)·∏_{P ∈ rowSG} XiAt(local factor at P), where phiT ∈ {0, 1} is the indicator that x maps the decoration of p into that of p' compatibly at every P ∈ ST (:646-648).
- thOf = thC = xiV … (vol = 1) (Formalism/LocalToGlobal/KernelLattices/XiIndCClass.lean:239): a finite sum over coset representatives (repsOf LamC/AC/BC) of the translated θ, i.e. a lattice theta kernel.
- So Kmap ≠ 0 at a built step ⇐ some x has phiT = 1 and all local factors ≠ 0, and xiV does not annihilate it. Neither was decided here; this is open limit 1.


## 11. Follow-up (22:02–22:08 UTC; logs/lk17.log): who supplies 'admitted', and is the kernel forced nonzero?
### 11a. 'admitted D' is proved for EVERY D7; the inert-place content is one D-independent local theorem
- Producer: second_theorem_closed_everywhere h70 h D (Formalism/SharpRestriction/GrowthBounds/PlaneCountTheoremWithConstants.lean:1121-1135, dist 13). For an arbitrary D it builds the whole 40-clause admitted D and applies h.p14R D (Loc6.p14R = loc_proved) to get ConclusionE D. It discharges the aux hypotheses with derived_statements_with_constants and split_theorem_with_constants.
- Route: linkG_fam_noKer ← h14 : SecondTheoremAtFamily (:975-977 = constants ∧ OrbitValuationBound Fm C) ← second_theorem_closed_form (:1740-1748, dist 12) ← second_theorem_absolute_constant (:1137-1142, dist 12) ← second_theorem_closed_everywhere applied to ONE datum D0. The family assembly family_second_bound_of_theorem (:941, dist 13) gives Adm6F.
- Key fact: aux₃ (_D : D7) (C : ℝ) (TableLineDegreeFive.lean:30-48, dist 13) ignores D. It is a uniform local statement: ∀ p ≥ 13 and every UnramifiedQuadraticSetting p (DVRs W ⊂ Qq, O ⊂ Kf, involution σ) with orbit and weight hypotheses, a valuation inequality holds. So "the gain holds for every family" is this single local theorem plus family bookkeeping.
- Inputs of second_theorem_closed_everywhere, marked:
| Input | Source | Mark |
|---|---|---|
| tool_statements_hold h70 D; plane_count_theorem_with_constants h70 h D | PlaneCountTheoremWithConstants.lean | NEW |
| TopChoiceBounds.orbit_positions_admitted_holds, orbit_counts_admitted_holds, symbol_space_orbit_holds | TopChoiceBounds | NEW |
| ~25 Joints.*_gives_orbit_count, valuation_transport_of_lemma, space_statement_of_lemmas, transport_statement_holds | Joints | NEW (assembly into orbit counts) |
| Loc6 Joints.Lit.weyl_dimension_formula_index_sets, branching_rule_count_form, pieri_rule_count_form | local_model_results (Formalism/LocalToGlobal/GroupLattices/AdmissibleFamilyNoKernel.lean:43-60, dist 12) | KNOWN (Weyl dimension formula, GL branching, Pieri rule) |
| Loc6 WittPI.wittIdentity_holds | same | KNOWN-type (Witt identity), statement not compared |
| Loc6 Model.dualInclusionLattices, bePerfect, lPOneEvery, flatTop, dualHeightBasis, hlSpans, testSide, fewBadOrbits, products_sections_integral; SelectionLemmas.Pf.secondPartClosed; Joints.joint_types_corollary_sentences, coarse_units_lemma | same | NEW |
| h70 : ExtremalProjectorTools = ExtremalProjector.extremal_projector_statement | AdmissibleFamilyNoKernel.lean:40-41 | KNOWN-named (extremal projector of Asherova–Smirnov–Tolstoy), statement not compared |
| h.p14R = loc_proved (admitted → ConclusionE for every D) | AdmissibleFamilyNoKernel.lean:284-286 | NEW (closed form, stmtE) |
- The candidate named in §6 (LocalDataSystem's FormalGroupSystem with per-lift bounds Σ_{j≤m} Fn, :79-91) is NOT the supplier. system_of_local_data has dist −1 (logs/lk12.log).
- Consequence for falsification: the inert-prime local content is a finite-hypothesis local inequality (aux₃) over unramified quadratic DVR settings with p ≥ 13, independent of any lattice or theta datum. It is the natural next probe target: build small instances of UnramifiedQuadraticSetting and test the inequality.
### 11b. Nonvanishing of thAt / Kmap: not forced on the cone
- rg over the 9 files opened in §1–§10 finds no lemma asserting thAt ≠ 0, Kmap ≠ 0, kAt ≠ 0 or xiAt ≠ 0.
  - The only match is the converse phiT_of_xiAt_ne_zero (xiAt x ≠ 0 → phiT = 1; XiAtPenultimateEqXiCpl1.lean:244).
  - A cone-wide name regex found no such theorem either (logs/lk17.log; its 22 hits are unrelated names).
- The only nonvanishing on the cone is ChIdx.Ok.nonvan (Formalism/LocalToGlobal/Columns/FiniteSumNonzeroLargeDepth.lean:1027). It is a field of I.Ok through idx_ok, with hm : R.m = m2 + 2:
  - for every step i < m − 2, index b and row prime P ∈ rowSG, some al has FiniteSumNonzeroLargeDepth (:988-989), i.e. a LOCAL finite sum E1at ≠ 0 for all large depths A', B';
  - this concerns the local factors of XiEntry at row primes only;
  - it excludes neither phiT ≡ 0 nor global cancellation, and says nothing for the last two steps (XiCpl1, xiLast).
- If the zero kernel realised Kmap, the route conclusion made vacuous is Conclusion5Sharp (Formalism/SharpRestriction/GlobalSchurBounds/ProductFirstZerothParts.lean:315-318), the first conjunct of aP10:
  - it states: ∀ w over T and c with v_w(c) = n ≥ e_w·m·ET, Into (c • colT) Nlam (S_m w), a denominator bound for the m-step column map colT;
  - with colT = 0 it holds trivially.
- That vacuity is harmless for the route. aP10 is an upper bound on denominators (a cost). ETDef (:320-322) sets tau = gain = F = 0 at the places over T, so the level prime carries no gain anyway.
- Nonvanishing would have to come from the rows that demand it (LinkDIV_at's BesselId: Bu u a nonzero multiple of besS; RRFull/NuFull spans; the ranks du_u in LinkS), which lie outside laws_true_fam.
- Net: laws_true_fam would remain true and useful as a cost bound even for a zero kernel. Its content is integrality at the split level prime, not a gain.

