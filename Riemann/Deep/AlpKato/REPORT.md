KNOWN: both targets are implied by published theorems — KatoTwoW asserts only the RATIONAL Bloch–Kato Selmer vanishing H¹_f(Q, V_f(e)) = 0 = H¹_f(Q, V_f̄(e)) for an odd-weight newform (Kato 2004), and LemmaK's bound log #H²(O_L[1/l], Z_l(e+1)) ≤ (2D)^{(2en)^C} follows from Bayer-Fluckiger–Emery–Houriet 2014 plus Quillen–Lichtenbaum; no counterexample.

# AlpKato: the in-archive theorems KatoTwoW and LemmaK of Alpöge's Lean archive — restricted, new, or false?

Archive: .scratch/alp/lean/src/ (READ-ONLY; source read, no Lean, no lake). Paths relative to it; DTC = Development/AutomorphicArithmetic/AutomorphicVarieties/ModularityArgument/ThresholdArguments/DenseThresholdCases, AV = Development/AutomorphicArithmetic/AutomorphicVarieties. Unfolded definitions with file:line: KATO.md. Probes: PROBES.md.

## Short answer
- KatoTwoW is KNOWN. Its BlochKatoSelmerVanishes (AV/AutomorphicForms/BlochKatoSelmerVanishes.lean:144) says every continuous 1-cocycle with values in V = (Q̄_l)² (Mathlib's PadicAlgCl, norm topology) that is unramified away from l and crystalline at l is a coboundary: the rational H¹_f(Q, V) = 0, i.e. FINITENESS of the integral Selmer group, not its vanishing. Under Used (KatoTwoOfUseN.lean:11-41) det W = θ and θ = χ_K (ker θ = G_K, θ² = 1), so f has weight 2e+1 and nebentypus χ_K, W ≅ V_f(e) and tw W θ = W ⊗ det W⁻¹ ≅ W^∨ ≅ V_f̄(e). For odd weight every critical value is non-zero (non-vanishing on Re s = (k+1)/2 plus the functional equation; Kato 2004 Thm 13.5(1)), and Kato 2004 (Astérisque 295, Thm 14.2: L(f, r) ≠ 0, 1 ≤ r ≤ k−1 ⇒ Sel(Q, W_f(r)) finite) gives both vanishings. Used is far stronger than needed. AlpDep REPORT §4b's "vanishing where Kato gives only finiteness" was a misreading of the definition.
- LemmaK is KNOWN, not NEW. Unfolded (LemmaK.lean:11-17, Statement.lean:54-73, SeparatedDvrSecondCohomologyMem.lean:68-102): ∃ C ∀ prime l, totally complex L ∋ μ_l (μ₄ if l = 2), e ≥ 1: H²_cts(Gal(L_{S_l}/L), Z_l(e+1)) is finite and log # ≤ (2D)^{(2en)^C}, S_l the primes over l (fixed by LData.hS). By Quillen–Lichtenbaum (Rost–Voevodsky; Rognes–Weibel at l = 2) this is #(K_{2e}(O_L) ⊗ Z_l); Bayer-Fluckiger–Emery–Houriet (arXiv:1406.5953, Thm 1.1) bound log card_ℓ K_{2e}(O_L)_tors by (2e·n)^α D^β, and the small-torsion loss is ≤ ℓ^{rank} with rank ≤ log₂ h_L + n. Together they give LemmaK with any C ≥ 11 (P2). Soulé 2003 alone does not (his bound is exponential in √D and his Thm 3 excludes exactly the primes l ≤ n+1 that LemmaK is about); Emery 2013 alone does not (K₂ only, inverts l | w_F).
- No counterexample. The probes found where each statement would break — an integral Selmer reading (CM form of weight 15 and level 7, l = 23 split, 23 | L^{alg}), or an LData with S larger than S_l (log #H² unbounded) — and the archive's definitions exclude both.
- Where the strength comes from: KatoTwoW from an in-archive Euler system theorem (zetaUr) plus the literature field LitCoh; LemmaK from an in-archive polynomial-in-D count of perfect-form cells (counting relative to a basis of common minimal vectors, plus HNF classes), the same improvement over Soulé that BEH 2014 published.

## 1. KatoTwoW (…/SphericalDecomposition/KatoTwoW.lean:12-16)
| Item | Finding | Evidence |
|---|---|---|
| Coefficients | rational: V = (Q̄_l)², Qbar l = PadicAlgCl l; no archive instance declaration changes its topology or that of the Galois groups (logs/topology_instances.log) | KATO.md B2 |
| Local conditions | away from l: c restricted to inertia is a coboundary (H¹_ur); at l: c ⊗ 1 is a coboundary in V ⊗_{Q_l} Bcr (H¹_f); Bcr = choice witness tied by a ring isomorphism eB to the archive's Per.Bcris (LocTied.lean:11-) | BlochKatoSelmerVanishes.lean:135-144, FiniteCohomology.lean:112-125 |
| Representation | Att forces tw W cyc^(−e) ≅ V_f, so W ≅ V_f(e), pure of weight 0; det W = θ = χ_K, θ² = 1 ⇒ tw W θ ≅ W^∨ (P5) | Newform.lean:43-48, KatoTwoOfUseN.lean:30 |
| Literature | Kato 2004 Thm 14.2 (L(f, r) ≠ 0, 1 ≤ r ≤ k−1 ⇒ Bloch–Kato Selmer group of the critical twist finite); odd weight ⇒ every critical value non-central and non-zero (functional equation + non-vanishing on Re s = (k+1)/2) | KATO.md C |
| Archive route | twoS ← VanishOrdL ← LitCoh + zetaUr (Euler system ES W A, ordinary, auxiliary g with l ∤ level) + SrcS | DTC/KatoTwoWOfPublishedH.lean:38-49, SphAtTw.lean:16-30, VanishOrdL.lean:13-26, ZetaUr.lean:18-21 |
| Probe | integral reading refuted in principle by f_7 = θ(ψ^{14}) (weight 15, level 7, CM by Q(√−7)), l = 23: L(f_7, 8)·π⁶/Ω^{14} = 7²·23/(2^{10}·3·5); Used excludes it (level 7 = N·7 forces N = 1, odd; CM by K violates irr) | PROBES.md P1 |

## 2. LemmaK (AV/GeneratorCountTuples/LemmaK.lean:11-17; ∃ CK at DTC/KFinal.lean:11)
| Item | Finding | Evidence |
|---|---|---|
| H2 | continuous 2-cocycles modulo d1 of continuous 1-cochains: Tate's continuous cohomology; Z_l with l-adic topology | SeparatedDvrSecondCohomologyMem.lean:68-102 |
| S | fixed: LData.hS : S = Above L l; GS = G_L/closure⟨inertia outside S⟩; chi = cyclotomic | Statement.lean:54-64, Global.lean:90-145 |
| Primes | HasMuQ ⇒ (l−1) | n: only l ≤ n+1 occur (P4) | LemmaK.lean:11 |
| Literature | finiteness Soulé 1979; size: BEH 2014 + Quillen–Lichtenbaum + rank ≤ log₂ h_L + n ⇒ LemmaK for every C ≥ 11 | PROBES.md P2 |
| Not sufficient | Soulé 2003 (exp(√D) shape; p ≥ n+2 only), Emery 2013 (K₂, l ∤ 6w_F) | PROBES.md P2 |
| Archive route | SFinal (Poitou–Tate, local finiteness, local Tate duality) ← ShaK (C = 5·(max 24 9 + 2)+11, NumThree.lean:216) ← SNum ← STors ← torsAll (C = 24) ← KT.tors_le (Gabber-type bound on the perfect-form complex) + KT.orb/nearcoord (orbits ≤ (2D·nN)^{12(p+1)(nN)³N}) | KATO.md A3-A4 |
| Probes | S enlarged ⇒ log #H² ≥ Σ v_l(Nv^e − 1) log l, unbounded (Q(i), l = 2: Σ v₂ = 235,237 over 78,437 primes of norm ≤ 10⁶): hS is load-bearing; small l: no hypothesis fails | PROBES.md P3, P4 |

## 3. Consequences for the route (AlpDep REPORT §4b)
- Two of count I's three key inputs (count_i_of_printed_bound_k_kato_automorphy, DTC/LinkedUnlE.lean:49-53), KatoTwoW and LemmaK, are known theorems in the form the archive states them; the third, AutomorphyUnramifiedCoefficientsAlt (automorphyUnramified_of_litN_of_gdEnd, DTC/CohInp.lean:18), was not examined. Whatever strength count I has lies in that input or in the combination (the count theorem's body was not opened).
- LemmaK also feeds ThmBab (with Poitou–Tate, local duality, Ferrero–Washington: all KNOWN), so ThmBab's inputs are all KNOWN; its window-level packaging was not re-examined here.
- Unaffected: cases II/III (IsAut4: automorphy of every constituent for every odd l, AlpDep §4b NEW) and case V (ThmS6Core_s, not judged).

## 4. Method and limits
- Read the statements and definitions; for LemmaK also the top-level proof steps (LemmaKOfSteps.lean, ExponentLe.lean, HermiteNormalForms/{Assembly,ColsMinSmul,Nearcoord,SmulStrictMono,SmulNeZeroLat}.lean); proof bodies below that level not opened. No Lean: instance resolution and the kernel check were not reproduced; my reading that V carries the norm topology rests on the absence of overriding instances (rg) and Mathlib's defaults.
- Residual citation risk: Kato 2004 Thm 14.2 was checked only through Kim's restatement. Used.level (AV/AutomorphicForms/BlochKatoSelmerVanishes.lean:89-115) allows l | N: W is then potentially ordinary at l with finite-order inertia characters of conductor dividing N|d_K| (LevelAt :102-109). Whether Thm 14.2 covers such l was not re-read in the original; the archive's own route first twists W to an auxiliary form g with l ∤ level(g) and a_l(g) a unit (SrcL), the good ordinary case.
- Not audited: whether the archive's concrete Per.Bcris is Fontaine's B_cris. (Continuity of W is part of Used: Level.cont, BlochKatoSelmerVanishes.lean:112.)
- Literature read in full text: Soulé 2003 (emis.de/ft/9497), Emery 2013 (arXiv:1303.6132), Bayer-Fluckiger–Emery–Houriet (arXiv:1406.5953). Kato 2004 Thm 14.2 (Sel(Q, W_f(r)) finite when L(f, r) ≠ 0) and Thm 13.5(1) (L(f, r) ≠ 0 for critical r unless k is even and r = k/2) checked as restated in C.-H. Kim, arXiv:2203.12157, Cor. 9.3 with its proof and Rem. 9.4; Rost–Voevodsky, Rognes–Weibel and Tate's K₂(Z[i]) = K₂(Z[ζ₃]) = 0 cited from memory.
- Compute: one process at a time under .scratch/memcap.sh 290 (peaks ≤ 49 MB, logs/).

## 5. Review (two reviewer subagents, no tools, evidence inline; logs/reviewer.md)
- KatoTwoW reviewer: C1 (rational reading) AGREE; C2 UNSURE (normalization det W = χ_K, non-vanishing argument, theorem numbers, scope in l), verdict INCOMPLETE. Fixed: det W = θ = χ_K derived explicitly; non-vanishing stated as boundary line + functional equation; theorem numbers checked against Kim's restatement (arXiv:2203.12157, Rem. 9.4, proof of Cor. 9.3). Open: the scope of Kato's Thm 14.2 at l | N was not re-read in the original (§4).
- LemmaK reviewer: every component AGREE (Quillen–Lichtenbaum comparison, finiteness, l ≤ n+1, loss ≤ ℓ^rank, rank via Kummer–Brauer, class-number bound); C5 UNSURE until BEH's truncation convention and explicit exponents are verified, verdict INCOMPLETE. Fixed: both verified in the full text of arXiv:1406.5953 (Thm 1.1 with explicit exponents, card_ℓ = order of A modulo the subgroup generated by elements of order ≤ ℓ) and quoted in PROBES.md P2 with the C = 11 derivation.
- After the fixes the verdict stays KNOWN for both. The one open citation point could at most make KatoTwoW "KNOWN except possibly at l | N", where it would be new but predicted by Bloch–Kato; it cannot produce a counterexample.

## Files
- Riemann/Deep/AlpKato/REPORT.md (this file), KATO.md (unfolded definitions, proof skeletons, literature table), PROBES.md (P1-P6).
- code/p1_cm7.py, code/p1b_cm7_weights.py, code/p1c_cm7_weights_pslq.py (CM forms of level 7, critical L-values, PSLQ), code/p2_bounds.py (BEH/Soulé vs kBound), code/p3_sfixed.py (localization sums), code/p5_twist_small_l.py (dual twist, excluded and admissible primes).
- logs/p1_cm7.log, logs/p1b_cm7_weights.log, logs/p1c_cm7_weights_pslq.log, logs/p2_bounds.log, logs/p3_sfixed.log, logs/p5_twist_small_l.log, logs/topology_instances.log (rg for topology instances), logs/reviewer.md (prompts, replies, fixes).
- Scratch (not for the repo): .scratch/alpkato/ (STATE.md, reviewer_prompts.md, downloaded PDFs and their text: soule2003, emery, houriet, kim).
