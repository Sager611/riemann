CORE MAPPED — the core is a Chebyshev/Stepanov-type UPPER bound for a weighted sum over inert prime powers (ASK ≤ MAIN + O(x^(1/2+ε))), whose 'multiplicity' is m ≍ √x theta-lift steps each contributing its own inert-prime divisibility of a period while a Bost slope inequality plus Chowla–Selberg periods bound the archimedean size with relative loss 1/m + m/x, and a windowed Landau theorem turns the one-sided bound into zero-freeness; tests A–C passed; mapped at statement level only (the proof that actual lifts satisfy the per-step hypotheses StepHyp is not located). Follow-up (CORE.md '## 9. StepHyp and VWin discharge'): the geometric per-step hypothesis (theta-kernel lattices with rational jets) is discharged by stepHyp_of_la and chained by chTh_ok (Formalism/LocalToGlobal/KernelLattices/KernelLocalFactors/ThDistC.lean:1871, 488-512), the definite-side chain hypotheses by nonempty_hyps_definite_side_of_closed_bound (SiegelEisensteinSecondTermHypotheses.lean:44), and VWin by vwin_of_row → v_window_of_literal_automorphy_n → v_window_of_automorphy_csd_ss_mono_odd (automorphy lifting, WindowBoundOfLiteratureLevel.lean:33), with the passage from the kernel chain to DIV's per-lift numbers still untraced.

# AlpCore report (alpcore, 20 Nov 2026, 17:54–18:40 UTC)

Archive paths are relative to .scratch/alp/lean/src/ (read-only). Full map: CORE.md (§8 = 1-page sketch).

## 1. The core, bottom-up
| Piece | Statement | file:line |
|---|---|---|
| (c) endgame | one-sided bound Gm ≤ C·y^θ on windows [Y, Y^(1+δ)] for arbitrarily large Y ⇒ L(s,ψ) ≠ 0 for Re s > θ; la = 3/2; θ ↓ 1/2 | Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/ZeroFreeOfRpowWindow.lean:11; Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/TopWin.lean:58-95; .../ExplicitFormulaErrorOfHostConditions.lean:11-41 (Gm, wt) |
| (c) target | Gz x = Σ Λ(n)(x/n−1) − x log x + (1+γ)x (remark item 5) | Formalism/StatementOC.lean:11-21 |
| arithmetic | ASK = Σ_{p inert} 2 log p Σ_{i odd-class} (a/p^i − w)^+; ASK − MAIN = character average of Gz/Gchi (proved, colour_balance) | Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean:38-110; ChainSlopeWindows/StronglyIrreducible.lean:695 |
| (a) Stepanov in Arakelov form | DivData (lattices RR, Nu u, integral maps Nu u → RR^∨); DIV, SlopesN, SlopeIneq, Count ⇒ mean ASK ≤ mean Main + ErrX | Formalism/LocalRingCardinality/DataArith.lean:139; Formalism/RankGapBounds/MeanBoundSecondForm.lean:37-200 |
| (a) multiplicity | ThmDiv over J = Fin m lifts: per-lift Γ-cost minus inert divisibility ≤ ph(a_ρ) + O(b·La); sum gives m·ph(a_t) in DIV | Formalism/SizeWindowBrackets/SBrDiv.lean:103-160, :165, :223-245 |
| (b) chain | Ranks (m even, N(j) = k + j(t+1), sigW j = (k+j−1, m−j)); chain datum Com (characters μ, η); per-step Rallis/positivity/seesaw statements | Formalism/LocalRingCardinality/Params.lean:6-70; Development/AbelianVarietyTheory/Definitions/LatticeFrames/ChainParSucc.lean:19-36; .../Proofs/DecoratedLatticeProofs/ThetaLiftChain*Statement.lean, DoubleLiftFockSumEqSumIteratedLiftsStatement.lean:14 |
| (b) √x | RestSize: x^(1/2−ε) ≤ m; explicit m ≍ x^(1/2−ε/16), n ≍ x^(2+8ε), b = exp(x^ε), a = xb | Formalism/SizeWindowBrackets/SChain.lean:60-75, SizesTextAt |
| cost (item 14) | error bracket nT·aE(1+m²/k) + m(2nmb x^ε + …) ⇒ relative 1/m + m/x | Formalism/SizeWindowBrackets/ErrorFormsWeakWindowsFalse.lean:133 |
| Chowla–Selberg | logOm = γ + log(π/v) − L'/L(1,χ_K), K = Q(√−d) | Development/AbelianVarietyTheory/Proofs/ModuliProofs/PrincipalRootReZeroImSq.lean:25-60 |
| GGP interlacing | HCParam u/v interlacing, gaps (t+3)/2; strictInterlace_ganGrossPrasad_big_small | Params.lean:48-70; .../Proofs/TwistedOperators/StrictInterlaceGanGrossPrasadBigSmallStatement.lean:7 |
| weak Vandiver | VWin: Σ_l log l Σ_{e<x^401} length Ш¹(G_{E,S}, Q_l/Z_l(−e)) ≤ exp(x^ε/2) for E ⊂ Q(ζ_N) | Development/AutomorphicArithmetic/CharacterBoxes/ExistsPageN.lean:12-49; Formalism/GaloisCohomology/Statement.lean:54-78 |

Classical vs new: Landau-type endgame, orthogonality, Stirling, Chowla–Selberg, Rallis/seesaw, supersingular divisibility (height-2 formal group) are classical in kind; new are the integral, nonvanishing m-step iteration with per-step divisibility (StepHyp), the windowed Landau formulation, the host-existence conditions (§3) and the proof of VWin.

## 2. Falsification tests (own code; every run under .scratch/memcap.sh 290)
| Test | What | Result |
|---|---|---|
| A (code/testA_cm_divisibility.py, logs/testA_K1000.log, peak 14 MB) | remark item 2 on y²=4x³−4x (Z[i]) and y²=4x³−4 (Z[ω]): v_p of the a-th derivative of ℘ − z^(−2) at 0, a ≤ 2000, exact rationals | inert p: v_p ≥ a/p + a/p³ + … − 0.96·log_p a for every a; the full v_p(a!) fails linearly; split p: v_p ≈ 0. PASS: inert-only, odd-power divisibility, exactly the classes ASK counts |
| B (code/testB_cost_bookkeeping.py, logs/testB.log) | explicit SizesTextAt parameters in log-space | RestSize, CostWindows, SizeWindows hold from log x ≈ 2.4e3/4.2e3/1.0e4 (ε = 0.04/0.02/0.01); relative error x^(−1/2+O(ε)); balance 1/m vs m/x ⇒ m ≍ √x; a per-lift cost ∝ m would give m ≍ x^(1/3) and exponent 2/3 (= remark item 14). PASS (documents Lean-checked arithmetic) |
| C (code/testC_main_normalisation.py, logs/testC.log, peak 55 MB) | one-character toy: ASK − MAIN = Gz − Gchi_K on real primes, K = Q(√−d), d = 7…199, y ≤ 2·10⁶ | observed |ASK − MAIN|/√y ≤ 1.15; flipping the sign of L'/L(1,χ) in MAIN drifts to 2400·√y. PASS: MAIN is the true main term; the core bound is GRH-scale, not vacuous, not contradicted |
No counterexample. Junk values: none involved (all tested quantities finite, positive arguments).

## 3. Least plausible steps (ranked; reviewer verdict in logs/reviewer_core.md)
1. StepHyp for the actual chain (SBrDiv.lean:223-245): every one of m ≍ √x lifts must carry its own inert-prime divisibility with Γ-cost exceeding it by ≤ ph(a) + O(b·La), uniformly. New, unexplained, and its proof site was not located. This is where Stepanov-over-Z would fail if it fails.
2. VWin (proved by vwin_of_row, …/DenseThresholdCases/VwinOfRow.lean:4): weak Vandiver for ~x^(2+8ε)-degree abelian fields and ~x^401 twists; beyond the literature (reviewer's first choice).
3. Host conditions (CORE.md §9): EAtH.bdG asks ‖G(θ,x)‖ ≤ x^(1/2+ε) at the point x for every 'Bounded' character θ ∈ (Xp∖Kept) ∪ Xm (ExplicitFormulaErrorOfHostConditions.lean:144-150; Analytic.lean:88-136). Plausible only if all of them carry the auxiliary q-/d-quadratic characters (conductor ≫ x^(O(1/ε)), chosen by averaging); Xm = nontrivial on the d-sign iod (ZsmulGeneratorEqZeroOfProjection.lean:104) fits, the 'first factor' A0 was not identified. Supplied by supplyChi (CharacterBoxes/SupplyChi.lean:8) from VWin.
4. Windowed Landau (ZeroFreeOfRpowWindow.lean:11): plausible, Lean-proved via lemmaL; specialised, not a textbook statement.
Reviewer's sign objection (ASK − MAIN is a SIGNED combination, real zeros of L(s,χ_K) enter with the wrong sign) is answered in the Lean by §9: the target Gm has positive weights on KEPT = ⟨ψ⟩ and the negatively signed characters are 'Bounded', controlled by host conditions (item 3).

## 4. Consistency: Euler product
Used essentially: primVal is multiplicative (FirstPrimeNeSecondPrime.lean:65); ASK and DIV are per-prime divisibility at inert primes (Formalism/LocalRingCardinality/DataDIVDefs.lean:45-55 InertK/NuDef); Frobenius classes come from a real abelian extension (HostFacts :47; TieArith, Formalism/SharpRestriction/GlobalSchurBounds/TieGlobal.lean:84-96). Davenport–Heilbronn functions have no Frobenius classes or CM periods, so the method does not apply to them; no inconsistency.

## 5. Against the Lean Zulip chain (CORE.md §10)
Agrees: RH ⇐ EveryCharacter ⇐ every_character_of_theta_integrals_local_matching (…/ThirdPlaceSeries/SiegelEisensteinSecondTermHypotheses.lean:69) with hV = VWin (i), hX = CovolumeBoundOfX (ii), P/F/hA/hM theta integrals and hL = RecTExists (iii); internal Galois representations via OpenLiteratureStatements (EveryCharacterOfLiteratureOfVWin.lean:215-218). Differs: the summary lists inputs only; this report locates the mechanism that turns them into zeros (inert-prime upper bound + windowed Landau), and reads (ii) as the auxiliary-construction (covolume) side of the slope method (not traced).

## 6. Not located / next
- Proof that actual lifts satisfy StepHyp (per-step divisibility, Γ-cost kerPlace DataDIVDefs.lean:74); the remark's explicit 'Stepanov determinant' / GGP-period determinant (names in Development/WeilBound/CoreStatements/JacobiBlockDeterminants/ are Gram/Jacobi minors, not obviously it).
- Proof of host existence (admissibleClassesTwistedChar) and of VWin.
- Next test with teeth: instantiate StepHyp's per-step balance phν − fa ≤ ph(a) + eE for a 1-step lift in the smallest Ranks (k = 6, m = 2, t ≥ 11) with real CM data, once its proof site is found.

## Files
- Riemann/Deep/AlpCore/REPORT.md (this), Riemann/Deep/AlpCore/CORE.md (map §0-§11, sketch §8).
- Riemann/Deep/AlpCore/code/: importgraph.py, graphq.py (import graph; every file reachable from Main), kwdecls.py (keyword census of 1.06 M declarations), testA_cm_divisibility.py, testB_cost_bookkeeping.py, testC_main_normalisation.py.
- Riemann/Deep/AlpCore/logs/: testA_K1000.log, testB.log, testC.log, reviewer_core.md.
- Scratch (not for commit): .scratch/alpcore/ (files.txt, decls.txt = 202 MB index of all 1.06 M declarations, reusable by alpaudit/alpmath, delete when done; mods.txt, edges.bin, nlines.bin, reach.txt, STATE.md).
