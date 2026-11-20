# OUTLINE — proof cone of Alpöge's riemannHypothesis (alpmath, 20 Nov 2026)

Paths relative to `.scratch/alp/lean/src/` (read-only copy of lean.tar.gz, sha256 9631b4f0…). "D/AA" = `Development/AutomorphicArithmetic`. Plain-language reading of each Lean statement; timeline item numbers refer to remark.pdf p. 3–4.
File-level cone: the import closure of `D/AA/CharacterBoxes/MainConclusion.lean` is 175,178 of the 175,191 non-Mathlib modules (logs/import_cone.log), so "on the cone" below means reached by name from the root, read by hand.

## L0 root
- `Main.lean:4` riemannHypothesis : Mathlib RiemannHypothesis := rhOfAll; `Main.lean:12` the four Dirichlet forms (Statements.lean, Mathlib-only).
- `Supplements/Kernels/MainTheorems.lean:7-16`: allChi := UpperTheory.FinalTheorems.everyCharacter; rhOfAll := Formalism.VW.rh_of_every allChi; critAll, claimAll from allChi.
- Main also imports `MainTheoremsAll/Critical` → SmallClassNumbers.{ClassNumberThree, LenEqZeroIff, OnePoint, CorGKOf}: built, but NOT in the import closure of MainTheorems (off the RH cone).

## L1 GRH for all Dirichlet L ⇒ RH (routine)
- `D/AA/CharacterBoxes/EveryCharacter.lean:13` EveryCharacter: L(s,χ) ≠ 0 for every N, χ mod N, Re s > 1/2, except the pole (χ = 1, s = 1). KNOWN-shape (GRH for Dirichlet L).
- `Formalism/GammaFactorNonvanishing/CritLine.lean:29` rh_of_every, `:11` critLine; `D/AA/CharacterBoxes/RiemannHypothesisOfEveryCharacter.lean:211,215` (ζ ≠ 0 on Re s > 1/2 from χ = 1, then the functional equation). Routine.
- `D/AA/CharacterBoxes/ClaimStrip.lean:59` strip: half-plane ⇒ strip by reflection (prim_reflect :42, gamma_ne_zero :35). Routine.

## L2 the top assembly
- `D/AA/CharacterBoxes/MainConclusion.lean:10` everyCharacter := everyCharacter_of_row galoisLocalGlobalStrong_recParameter.
- `MainConclusion.lean:7-8` everyCharacter_of_row h := every_character_of_literature_of_vWin ⟨h, recTExists, h⟩ (vwin_of_row h).
- Inputs: (a) Row0b = Galois representations for regular algebraic polarized cuspidal π on GL_a over CM fields with local–global compatibility (`D/AA/Foundations/GaloisRepresentations/TypeConditions/WeightsConjBis.lean:69`; Chenevier–Harris, Shin, Caraiani, BLGGT: KNOWN); (b) RecTExists = local Langlands for GL_n with its properties (`D/AA/Foundations/LocalLanglands/LocalCentralUnits/RecTExists.lean:20`, proved there; Harris–Taylor, Henniart, Scholze: KNOWN); (c) VWin (below).

## L3 the three branches
- B1 (literature bundle) `D/AA/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean:215` structure OpenLiteratureStatements := {caraianiLocalGlobalDual : Row, localLanglandsExists : RecTExists, galoisRepOfPolarizedPairLit : Row}; 196 declarations there turn it into ~190 "Published" statements: 152 ignore O and cite in-archive proofs (*_holds, e.g. `:219` logDerOne_holds, `:221` ellipticMainTheoremCM_holds, `:223` heckeTateFE), the rest route O's three fields (`:283` O.localLanglandsExists; Row at `:286,:460,:618`) (logs/lit_bundle.log).
- `…EveryCharacterOfLiteratureOfVWin.lean` (last theorem) every_character_of_literature_of_vWin O hV := UpperTheory.Lit.every_character_of_theta_integrals_local_matching (published_results_eigen O) (published_results_first O) (theta integrals agree, local matching) (topCovolSlim2_use_alt_alt …: covolume bound) (recT) hV.
- B2 (Galois representations) `D/AA/Foundations/GaloisRepresentations/IntegralStructures/RowSrc.lean:6` galoisLocalGlobalStrong_recParameter := galoisRepresentations_of_thetaLift tLiftK (`…/GaloisRepresentationsOfThetaLift.lean:6`) := galoisLocalGlobalStrong_of_geometry_lift geoOne_holds tLiftK (`…/IntegralStructures/Development.lean:30`).
- B3 (Vandiver weakening) `D/AA/AutomorphicVarieties/ModularityArgument/ThresholdArguments/DenseThresholdCases/VwinOfRow.lean:4` vwin_of_row h := automorphic_window_of_inputs (OpenImageCondition.mk h) baseChangeN h (`…/TopMainZj.lean:71`).

## L4 inside B1 (analytic engine; NEW)
- `…/ThirdPlaceSeries/SiegelEisensteinSecondTermHypotheses.lean:69-72` every_character_of_theta_integrals_local_matching P F hA hM hX hL hV := every_character_of_green_v_win (green_of_regular P.greenRegularCharacters) hV (nonempty_hyps_definite_side_of_closed_bound … `:44`).
- `…/ThirdPlaceSeries/TopPlug.lean:16` every_character_of_green_v_win hG hV h := imprimitive (zero_free_all_of_green_v_win hG hV h); same shape as `D/AA/CharacterBoxes/EveryCharacter.lean:16` plugAll.
- Green (`Development/AbelianVarietyTheory/Definitions/NonGaloisDescent/LabelDatum.lean:108`): Green's parametrisation of cuspidal characters of GL_d(F_g) by irreducible polynomials (degree ∏(g^i − 1), values on elliptic elements); from GreenRegularCharacters (`…/PublishedTheorems/GaloisCohomology₂/GreenRegularCharacters.lean:13`). KNOWN (Green 1955). Probes G1–G8.
- `Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/ExplicitFormulaErrorOfHostConditions.lean:14,18` ZeroFreeOn Q / ZeroFreeAll: every primitive ψ (PChar) has L(s, ψ) ≠ 0 for Re s > 1/2. `D/AA/CharacterBoxes/Imprimitive.lean:5` imprimitive: primitive ⇒ all χ (Euler factors at p ∣ N; routine).
- `D/AA/CharacterBoxes/ZeroFreeAllOfTopG.lean:26` zeroFreeAll_of_topG (la = 3/2 in `EveryCharacter.lean:10`): ZeroFreeAll from TopG (`:12`: for la > 1 and every requirement G with HostExW, if each family in the class satisfies HypV then the class is FamilyClassConnected) and HostExW for every ψ.
- `…/TotientGaloisClasses/TopWin.lean:26,30,90` zeroFreeOn_of_tiedWin_cyc: if for every family with a Chebyshev inequality on windows, Λ.Gm(x) ≤ x^{1/2 + (eD + eL + 7)ε} for large x in the windows (`:26`), and such families exist at window scale la for all small ε (`:30`), then no zeros right of 1/2. Plain reading: a square-root Chebyshev bound on geometric windows, uniform in ε → 0, gives GRH (Landau-type; ROUTINE given the bound). Timeline 5–6 (Guinand–Weil sign, W(1_S) ≥ 0, "> 1 − η").
- `D/AA/CharacterBoxes/DivinNOfHypV.lean:10` HypV: the host fields' Sha lengths satisfy VVle with B = √b; `:13` FamilyClassConnected: for all small ε there is a family whose windows are ChainTied (the tied chain of lifts; timeline 7–10, 14–16). This is the claimed new mathematics (Stepanov determinant with multiplicity from m ≍ √X theta lifts): NEW; no finite consequence isolated (all statements are asymptotic in x with existential constants).

## L4 inside B2 (Galois representations; KNOWN statements)
- `…/GaloisRepresentations/IntegralStructures/Development.lean` galoisLocalGlobal_of_geometry_lift geo up := rest (coreAll_of (seedFinOdd_of_geoOne geo) (congFinOdd_of_ex wtLacEx_holds) up): Galois representations by congruences/limits from seeds of odd rank.
- GeoOne (`…/IntegralStructures/GeometricFlagBases/HodgeTrAt.lean:79`): for automorphic U on a definite unitary group U(a), a odd, a Galois representation ρ with local–global compatibility at good places (Kottwitz, Harris–Taylor, Shin: KNOWN); proved by geoOne_holds (`…/GeometricFlagBases/WeightOneGeometricRealization.lean:11`).
- TLiftK (`…/ThetaParseval/HalfIndexLattices/TLiftK.lean:12`): every automorphic U on U(n) with weights in range lifts by theta to U(d) for all d ≥ d₀ (stable-range non-vanishing: Rallis, J.-S. Li: KNOWN); proved by tLift_holds_k (`…/ThetaKernels/TLiftK.lean:4`). Timeline 8 ("iterated theta lifts U(X₀) → … → U(X_m)").

## L4 inside B3 (Vandiver weakening; NEW, consistent)
- VWin (`D/AA/CharacterBoxes/ExistsPageN.lean:44-49`): ∀ primitive ψ, ∀ ε ∈ (0, 1/52), ∀ 2-power page e(x) with 2^e ≍ x^{2+8ε}, ∃ J, X₀: for X ≥ X₀ and every WinPair (q ≡ 3 mod 8, d ≡ 7 mod 16 primes, q ∣ d+1, every odd p ≤ Y^{ε/5} divides d+1, log q ≍ Y^{ε/5}, log d ≍ Y^{ε/5} log²Y) there is Y ∈ [X^{1/J}, √X] with: for every subfield E ⊂ Q(ζ_N) cut out by a ≡ ±1 mod 2^{e+2}, ψ(a) = 1, (a/q) = (a/d) = 1 (IsChiField `:19`), Σ_l log l Σ_{0<e′<x⁴⁰¹} len Sha¹(G_{E,S_l}, Q_l/Z_l(−e′)) ≤ exp(x^ε/2) (VVle `:12`; Sha, shaLen `Formalism/GaloisCohomology/Statement.lean:75,78`). By Poitou–Tate and Rost–Voevodsky these Sha are duals of H²(O_E[1/l], Z_l(1+e′)) ⊇ K-theory l-parts (K_{2e′}(O_E) ⊗ Z_l). Timeline 17 and footnote 6.
- Consistency: the bound is far above what the Lichtenbaum/Brauer–Siegel size of these groups allows to be violated (log-size ≲ t²[E:Q] log N ≈ x⁸⁰²·exp(x^{ε/5} log²x) ≪ exp(x^ε/2)); no finite case is asserted. NEW-and-consistent, unfalsifiable by finite computation.
- Footnote 6's "vanishing of K-groups of Z": the Lean statements are `D/AA/AutomorphicVarieties/SmallClassNumbers/ClassNumberThree.lean:16` Bound (∃ C, A: Σ_{l ∈ Ls} log l·len_l(e) ≤ C(e+1)^A over primes l ≥ 5 with (l−1) ∤ 2e, e even, over Q) and `…/LenEqZeroIff.lean:65` Large (∃ C, A: Sha¹(Q, Q_l/Z_l(−e)) = 0 when log l > C(e+1)^A); proved at `…/CorGKOf.lean:27` bound, `LenEqZeroIff.lean:70` large_of. Main build only (off the RH cone). For e even these groups are dual to K_{2e}(Z) ⊗ Z_l = K_{4k}(Z)_l: Large says K_{4k}(Z) has no l-torsion for log l > C(2k+1)^A, i.e. |K_{4k}(Z)| ≤ exp(poly(k)) up to its small primes. Consistent with Vandiver for p < 2^31 (Hart–Harvey–Ong 2017; Kurihara/Soulé: Vandiver at p ⇔ no p-torsion in K_{4k}(Z) for all k), with K₄(Z) = 0 (Rognes 2000) and K₈(Z) = 0 (Dutour Sikirić, Elbaz-Vincent, Kupers, Martinet 2019, "Voronoi complexes in higher dimensions, cohomology of GL_N(Z) for N ≥ 7 and the triviality of K₈(Z)"); constants existential, so no table can contradict it. "Cyclicity up to 100" has no Lean counterpart that we found (rg Vandiver/irregular/KGroup over the declaration index: none).

## Timeline of remark.pdf ↔ Lean nodes
| remark item | Lean node (file:line) | status |
|---|---|---|
| 1 Stepanov determinant; 4 W(1_S) ≥ 0 via a determinant; 7 GGP periods with interlacing weights | FamilyClassConnected / ChainTied (`D/AA/CharacterBoxes/DivinNOfHypV.lean:13`), TopG (`ZeroFreeAllOfTopG.lean:12`); WeylBranching (interlacing, `…/PublishedFoundations/FiniteIntertwining.lean:129`, probe W1) | NEW (engine); interlacing input KNOWN |
| 2 Dp ≡ 0 on a CM curve, p inert | CMCurveGoodOutsidePrimeSt (literature list, `…EveryCharacterOfLiteratureOfVWin.lean:277`) | KNOWN-type input (not probed) |
| 3 qN = M k + R, sign-only Weil/Deligne | (Weil-bound tree `Development/WeilBound/…`; not traced) | not traced |
| 5–6 Guinand–Weil sign; W(1_(0,1] + X^η 1_(0,1/X]) ≥ 0 ⇒ zero-free > 1 − η | ZeroFreeOn via Chebyshev inequality on windows (`Formalism/…/TotientGaloisClasses/TopWin.lean:26,90`) | ROUTINE given the window bound |
| 8 iterated theta lifts U(X₀) → … → U(X_m) | TLiftK (`…/ThetaParseval/HalfIndexLattices/TLiftK.lean:12`), tLiftK (`…/ThetaKernels/TLiftK.lean:4`) | KNOWN (stable range) |
| 9 Chowla–Selberg, Weyl, Siegel mass | `GroundTheory/HeightsAndTheta/Heights/ChowlaSelberg/Assembly.lean:32-39` (probes A1–A8); HilbertClassFieldJ/CM (probes H1–H4); SiegelMassRatioCM in `D/AA/CharacterBoxes/ClassSel.lean:14` | KNOWN, probed |
| 11, 19 Iwasawa tower, Behrend | IsChiField/PageN/WinPair (`ExistsPageN.lean:19,27,39`): 2-power cyclotomic page 2^e ≍ x^{2+8ε}, auxiliary primes q, d | NEW construction |
| 12 seesaw, Rallis | RallisInnerProductConvergent, SiegelWeilUnitary* (literature list) | KNOWN inputs |
| 13 Bost slope inequality | covolume branch topCovolSlim2_use_alt_alt (`…EveryCharacterOfLiteratureOfVWin.lean:812`), PublishedHodgeBundleTheorems (arithmetic Hilbert–Samuel, Moriwaki) | KNOWN inputs |
| 16 relative trace formula | QuadraticClassFieldTheory, HeckeLEntireAndNonzeroAtOne (TwistedTraceFormula folder; probes Q1, L1) | KNOWN inputs |
| 17 Vandiver weakening (footnote 6) | VWin (`ExistsPageN.lean:49`), vwin_of_row (`VwinOfRow.lean:4`); over Q: Bound, Large (SmallClassNumbers, Main build only) | NEW, consistent |
| 20 GRH | EveryCharacter (`EveryCharacter.lean:13`) | target |
| (Galois input) | Row0b (`WeightsConjBis.lean:69`), GeoOne (`HodgeTrAt.lean:79`), Green (probes G1–G8) | KNOWN |
Depth reached by hand: L0 Main → L1 EveryCharacter → L2 everyCharacter_of_row → L3 three branches → L4 engine/TLiftK/GeoOne/VWin → L5 definitions (TopWin, DivinNOfHypV, Statement.lean Sha). Depth 5.
