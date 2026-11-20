FAITHFUL-KNOWN: the formal AdZero is the Bloch–Kato adjoint Selmer group H¹_f(L+, ad r) (B_cris at l, unramified elsewhere), not an ordinary one, and the AdjointLiftingTheorem is a special case of Thorne, arXiv:2207.04925, Theorem 5.1 (RACP π over a CM field, r|G_{F(ζ_{p^∞})} irreducible ⇒ H¹_f(F+, ad r) = 0, no residual hypothesis); its ordinary core OrdZero follows routinely, and every known non-étale test point is excluded by a named hypothesis.

# AlpAdj: the in-archive AdjointLiftingTheorem (DTC/AutomorphicPointGalois.lean:35-39): faithful and new, known, or a discrepancy?

Archive .scratch/alp/lean/src/ (READ-ONLY; no Lean, no lake). AV = Development/AutomorphicArithmetic/AutomorphicVarieties; DTC = AV/ModularityArgument/ThresholdArguments/DenseThresholdCases; OZ = AV/TheoremsMuLinked/OrdZero.lean; UODR = AV/WeilGroupCharacters/UniversalOrdinaryDeterminantRing.lean. Details and quotes: ADJ.md; numerics: PROBES.md.

## Short answer
- AdZero V ρ (OZ:68-72) = every cocycle of G_{L+,S} in the Clozel–Harris–Taylor adjoint ad~ρ (R(c)X = −P⁻¹XᵀP, AV/WeakShimuraTheory/EigenbasisLifts/IsTanZero.lean:31-34) that satisfies IsBK (unramified at w ∤ l; trivial in V ⊗ B_cris at w | l; AV/AutomorphicForms/BlochKatoSelmerVanishes.lean:135-140, AV/PadicHodgeTheory/FiniteCohomology.lean:112-114) is a coboundary: H¹_f(L+, ad r) = 0. B is Fontaine's B_cris via LocTied (Formalism/LocalRingCardinality/AdicCompletionLocal/LocTied.lean:11-30).
- The ORDINARY statement is the companion core DenseThresholdMainCase (DTC/AutomorphicPointGalois.lean:30-34): OrdZero (OZ:62-64) = Greenberg's fixed-weight ordinary Selmer group (upper triangular at w | l with diagonal unramified, OZ:22-27; unramified at St; relaxed at v0). The archive derives AdZero from it (DTC/ThFinal.lean:18). AlpAut's AUT.md §6 called AdZero 'ordinary': that reading was wrong; both objects are standard.
- Hypotheses ⊇ Thorne's: L CM; 2 ≤ n ≤ 4; AutCsd0 = cuspidal, regular algebraic, conjugate self-dual π with Frobenius matching (Foundations/GaloisRepresentations/AdmissibleResidues/IsAutomorphicOverImagQuad.lean:35-36), so ρ ≅ r_{π,ι}; TypeI = ρ|G_{L(ζ_{l^∞})} irreducible (OZ:14, G' = ker χ_cyc, UODR:73). Thorne (TeX l.400): 'no condition on the residual representation (which might even be trivial)'. The residual diagonal datum, the Hida setting and B enter only the hypotheses; the conclusion is about ρ alone.
- The brief's premise is false: Newton–Thorne (arXiv:1912.11265) need the p-ADIC image r(G_{F(ζ_{p^∞})}) 'enormous', Thorne 2022 only irreducibility; only Allen (arXiv:1411.7661) needs residual adequacy.
- OrdZero is a routine corollary of Thorne: ordinary deformations are de Rham (weights strictly decrease up the flag), Free forbids ψ_i = ψ_j·χ_cyc^{−1} so H¹_g = H¹_f at l, and purity makes v0 harmless (ADJ.md §8).
- No counterexample: weight one fails Arith/IsRegAlg/Free/Par; the split-multiplicative exceptional zero fails Free; CM-by-L, Yoshida, Saito–Kurokawa and U(2)×U(1) points fail TypeI and cuspidality; the trace part dies under fixed weight; the Euler characteristic is balanced (§2).

## 1. Definitions (standard notation; ADJ.md §1-§4)
| Name | file:line | Meaning |
|---|---|---|
| HidaSetting, Lev, frame | AV/WeilGroupCharacters/ContinuousResP.lean:73-162 | l odd; L/L+ unramified at all finite places; v | l, St, v0 split; Γ = G_{L+,S}, G₀ = G_{L,S}, c, ε = χ_cyc, Gv = G_{L_{w(v)}} |
| Resid, V0ok | UODR:109-117; AV/WeilGroupCharacters/ArithmeticPointClassicality.lean:14-16 | n pairwise distinct residual characters, polarized; unramified at v0, equal on D_{v0} |
| UniversalOrdinaryDeterminantRing B | UODR:138-161 | universal polarized ordinary Chenevier determinant, local characters ψ_{v,i} free (Hida family over LamObj, :168-176), minimal at St, det free up to polarization |
| x.Arith, Par | ArithmeticPointClassicality.lean:55,72; OZ:56-58 | regular locally algebraic weights; parallel regular weights at each v | l |
| IsAutCsdA | AV/ModularityArgument/AutomorphyLifting/IsAutCsdA.lean:11-12 | AutCsd0 (RACSDC, attached) ∧ WitA (local conditions, satisfiable: AV/ModularityArgument/AutomorphyLifting/WeakOfAtt.lean:27-31; LeviSmoothLocalization/WinUpDown.lean:38-41) |
| TypeI, N3, Free | OZ:14-20 | ρ|G' irreducible; regular semisimple elements dense in G'; ∃ g ∈ G_w ∩ G' with ψ_i(g)/ψ_j(g) of infinite order |
| AdZero | OZ:68-72 | H¹_f(L+, ad~ρ) = 0 (Bloch–Kato) |
| OrdZero | OZ:31-34, 62-64 | fixed-weight Greenberg ordinary H¹ of ad~ρ = 0 |

## 2. Test points (brief task 2; ADJ.md §7)
| Point | Excluded by | file:line |
|---|---|---|
| (a) weight-one dihedral RM (Bellaïche–Dimitrov; P1: 14 forms at N ≤ 150, e.g. N = 39, RM by Q(√13), Satake ratios roots of unity at l = 17, 23, 29, 43) | h.Incr (HT weights (0,0)), IsRegAlg, Free, Par | ArithmeticPointClassicality.lean:55,72; IsAutomorphicOverImagQuad.lean:36; OZ:19-20, 56-58 |
| (b) split multiplicative at l (P2: y² − 11xy + 3y = x³, N = 2118, a_3 = +1, 3-torsion; Tate normal form t = −25, N = 4495, a_5 = +1) | Free (ψ_1/ψ_2 = χ_cyc, equal on G_w ∩ G'); also ψ_i ≠ ψ_j·χ in OrdVan2Alt | OZ:19-20; DTC/OrdVan2.lean:19 |
| (c) CM by L or by K ⊂ L | TypeI (ρ|G_L = ψ ⊕ ψ^c), cuspidality | OZ:14; IsAutomorphicOverImagQuad.lean:36 |
| (d) Yoshida, Saito–Kurokawa, U(2)×U(1), Eisenstein | TypeI, cuspidality | same |
| (e) trace part | scalars = η_{L/L+}; fixed weight leaves Hom(Cl_L, Q̄_l)^{c=−1} = 0; free weight would give rank d (P3: Z_l-ranks 2, 2, 3 = d + 1 for Q(√−6), Q(√−14), Q(ζ_8)) | OZ:27; IsTanZero.lean:34 |
| (f) Greenberg–Wiles count | balanced: d·n(n−1)/2 − d·n(n−1)/2 = 0 (P4) | — |
- Even the points excluded only by Free are not counterexamples: at the Tate curve the fixed-weight Greenberg line equals H¹_f(G_l, ad⁰) (unit Kummer line; c_q dies), and H¹_f vanishes (Thorne). Free is a hypothesis of the archive's proof, not of the truth of the statement.
- Non-vacuity: a non-CM E/Q with a rational l-torsion point, good ordinary at l, base-changed to a CM field L with L/L+ unramified, l and the bad places split, ζ_l ∉ L, and v0 with N v0 ≡ 1 mod l, satisfies every hypothesis (residual 1 ⊕ ω).

## 3. The proof one level down (brief task 3; ADJ.md §6)
- adjoint_lifting_of_cohomology (DTC/ThFinal.lean:18, cone dist 11) = adjointLiftingTheorem_of_locPrints_of_core (DTC/AdLOf.lean:104, dist 12) (locPrints_of_litN hH.litN, DTC/HVFinalU.lean:12-13 := ⟨coh.lfK, ltK, leK, lf, lt, le, bk⟩: local finiteness, Tate duality, Euler characteristic, Bloch–Kato local; used for the passage OrdZero ⇒ AdZero at l, body of AdLOf.lean:104 not traced) (core_alt_bis hH hC, ThFinal.lean:17, dist 11).
- Core: ordinary_vanishing_two_of_printed_hida (DTC/TOrdVan.lean:159-192, dist 12) proves OrdVan2Alt (DTC/OrdVan2.lean:13-20, dist 13) by Hida theory (hida_x_of_printed_hida, DTC/THidaX.lean:17-23, dist 13: unitary Hecke data, a Hecke ideal J of the ordinary determinant ring through which every automorphic point factors) and patching at auxiliary places chosen from regular semisimple elements of ρ(G') (RegDet q sg, AV/TheoremsMuLinked/OrdinaryHeckeAlgebraIntegral.lean:156-157; Pa.CohP, AV/PatchingModules/TaylorWilesObstructionBound.lean:23; core_o_of_literal_automorphy_n, DTC/TCoreO.lean:129, dist 13). This is the Newton–Thorne 2019 mechanism (Taylor–Wiles places with regular semisimple Frobenius, TeX of 2207.04925 l.392, 398), not Skinner–Wiles and not residual adequacy. Mark: KNOWN (same statement; Thorne 2022 is stronger: irreducibility only, no N3, no ordinarity).
- Inputs (projections found in the 23,970 theorem files of its proof cone, logs/ordvan_projections.log, an over-approximation): hH.litN with aut.{art, berger, gk, js, gal, chA, bc, chW, mw, fin}, coh.{lt, lf, le, pt, gw, gf, ge, uo, bk}, loc.{bz, zel, geo, llc, conv, ...}; quat, smoR, gkR, admR, galR, conjC, ssR, descU. All are classical foundations (local Langlands and Bernstein–Zelevinsky, Jacquet–Shalika, Galois representations of polarized pairs, Berger, Poitou–Tate, Greenberg–Wiles, Bloch–Kato local); aut.car, aut.ctemp, aut.pur do not occur. None is an adjoint-vanishing, R = T or automorphy-lifting statement, so no 'literature' Prop can have dropped an adequacy or big-image hypothesis: the vanishing is re-proved in the archive. Marks: KNOWN (foundations). Leaf bodies not opened.

## 4. Probes (PROBES.md)
- P0 TeX of 2207.04925: Theorem 5.1 and l.400 quoted; P1 weight one; P2 exceptional zero; P3 trace part (ray class groups); P4 Euler characteristic. Peaks ≤ 64 MB RSS.

## 5. Review (logs/reviewer.md)
- Pass 1 (inline E1-E9): INCOMPLETE, blocking 'B-trivial = crystalline not shown'; also Schur for the quantifier over P and the Frobenius normalization. Fixed: LocTied gives an equivariant, Frobenius-compatible ring isomorphism B ≅ Per.Bcris (Fontaine's construction from A_inf); Schur and the Harris–Taylor normalization added (ADJ.md §8). Pass 2: AGREE.

## 6. Consequences
- AlpAut's input (A) (AUT.md §6-§7: 'adjoint-Selmer vanishing at every automorphic point of a residually diagonal family, n ≤ 4, over CM fields' marked as beyond the literature) is KNOWN: Thorne 2022 Theorem 5.1 for H¹_f; OrdZero as a corollary. The archive's novelty on that route must sit in (B) promodularity through residually reducible deformation spaces (PMExAlt, DTC/PotentialModularityTower.lean:38-66) and (C) the seeds. Not a discovery; candidate COROLLARIES/NOGO wording is for the orchestrator.
- No interpretation discrepancy: AdZero and OrdZero are standard objects; only the commissioning description (H¹_ord for AdZero; residual-image requirement of NT) was wrong.

## 7. Limits
- Statement-level reading plus top-level proof terms; leaf bodies (Hida data, patching, Per.Bcris construction bodies) not opened; no build (forbidden).
- Literature: Thorne and Newton–Thorne read in the arXiv TeX; Allen from its abstract. Theorem number 5.1 computed from the TeX counters (label thm_vanishing_over_CM_field).

## Files
- Riemann/Deep/AlpAdj/REPORT.md (this file), ADJ.md (definitions, literature match, proof, test points, review fixes), PROBES.md (P0-P4).
- code/probes.py (PARI probes P1-P4), code/reach_lit.py (reachable literature Props), code/reach_files.py (theorem files of a proof cone).
- logs/P1.log, logs/P1_rerun.log (with /usr/bin/time), logs/P2.log, logs/P3.log, logs/P4.log, logs/reach_ordvan2.log, logs/ordvan_projections.log, logs/reach_files.err, logs/reviewer.md.
- Scratch (not for the repo): .scratch/alpadj/ (STATE.md, tex.py, texl.py, thmno.py, 2207.04925.src, ordvan_files.txt).

