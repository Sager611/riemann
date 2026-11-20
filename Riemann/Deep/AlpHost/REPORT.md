KNOWN: PairEx and ThinH are formalized classical analytic number theory with matching ranges — the one good pair (q, d) comes from the standard 'almost all quadratic twists' argument at sum lengths polylogarithmic in the conductor (quadratic reciprocity + Gallagher–Linnik PNT in APs with Deuring–Heilbronn, Landau–Page, Siegel; Montgomery's zero density near σ = 1; Gallagher–Sobolev for the continuum of y), and ThinH is the de la Vallée Poussin PNT for characters with exceptional zeros excluded by RealBox; no circular input and no false statement found (reviewer: KNOWN-CONFIRMED after 4 rounds).

# AlpHost: the hosts of L. Alpöge's Lean archive (PairEx, pair_share, ThinH) — known averaging, new, or false? (alphost, 20 Nov 2026, 21:00–22:10 UTC)

Archive src = .scratch/alp/lean/src (READ-ONLY; no Lean, no lake); AA = src/Development/AutomorphicArithmetic; JCS = src/Formalism/GL3MeanSquares/GenusMeanAsymptotics/JacobiCharacterSums. Statements in standard notation and the full proof tree with file:line: HOSTS.md; numerics: PROBES.md. BFS distance to riemannHypothesis from .scratch/o66_q.py: every declaration below is on the cone (T = theorem).

## Short answer
- The average-to-individual step is classical and harmless: PairGood asks GRH-size bounds |Gchi(χ, y)| ≤ y^{1/2+ε/2} only for y ∈ [Y, Y^A] while the conductor is ≈ e^{Y^{ε/5}log²Y}, i.e. for character sums of length (log q)^{O(1)}; for such short sums 'square-root cancellation for all but a vanishing share of the family' follows from orthogonality over the family (the sums are shorter than √q), not from any GRH-type input. GRH itself would give these bounds for every pair (margin y^{ε/10}).
- The zero-free boxes (Re s > 1 − 1/400 up to height ≈ conductor, also for j = 1 where the height exceeds d) cost nothing: the archive's ZeroDensityBox (≤ C Q^{1/4} bad primitive characters of conductor ≤ Q, box height Q) is far weaker than Montgomery's (Q³)^{2/399} log^{14} Q, against family sizes ≥ e^{0.9 z} and e^{0.9 z'}.
- ThinH needs y ≥ e^{x^{ε/3}} with log-conductor ≤ 2x^{ε/5}log²x: the classical saving exp(−c log y/(λ + √log y)) ≳ exp(−c x^{2ε/15}/log²x) beats the required e^{−x^{ε/8}} because 2/15 > 1/8; real characters with q | cond or d | cond need RealBox (Siegel's bound would be useless at this conductor size), and the others have conductor ≤ 8B and are handled by L(1, χ) ≠ 0 + finiteness, so c = 1 is uniform in B.
- The target ψ enters only as M = ψ.N; the GRH-size bounds are demanded only for ψ^k·χ_q^i·χ_d^j with (i, j) ≠ (0, 0) (boundedTw_key), the kept ψ^k get only ThinH. No input of pair_share or ThinH uses anything about zeros of ψ.
- Consequence for the RH tree: node (I2) of Alpöge's chain (alpchain) is classical; the load-bearing new mathematics, if any, is in (I1) VWin and (I3) the theta-lift lattice data. Checking hosts cannot confirm or refute RH.

## (1) The statements (HOSTS.md §1–5, 7)
- PairGood M ε Y X q d (AA/CharacterBoxes/PairGood.lean:19-27), PairEx (:29-31), RealBox (:33-35), ThinH (:37-40); WinPair (AA/CharacterBoxes/ExistsPageN.lean:27-37): q ≡ 3 (8), log q ∈ (z/2, z]; d ≡ 7 (16), log d ∈ (z log²Y/2, z log²Y], d ≡ −1 mod q and mod every odd p ≤ z (so χ_d(p) = +1 for all p ≤ z and p = q; PROBES (a1) confirms L'/L(1, χ_d) ≈ −Σ_{p ≤ z} log p/(p − 1)); z = Y^{ε/5}.
- An.gChi (src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/Analytic.lean:14-24) is the SMOOTHED signed sum −Σ_n Λ(n)χ(n) min(1, y/n) (for χ = 1: G(y) = −y∫_y^∞(ψ(u) − u)u^{−2}du), not ψ(y, χ); An.Induces (:26-27) = 'the primitive character inducing f'; logDerivativeAtOne = L'/L(1, χ); ChebyshevCharacterBound (:145-147), BoundedCharsPointwise (:154-156), ExceptionalCharacterKept (:158-159); IsCharE = primitive characters of the host group G0 × MA, conductor | QS·q·d.

## (2) PairEx down to its analytic inputs (HOSTS.md §6, §8)
| Input | file:line | statement (ranges) | status, reference |
|---|---|---|---|
| QH.meanSqK | AA/MeanSquareBounds/meanSqK.lean:51-54 | Σ_{q ∈ P}\|Σ_{n ≤ Y^{400K}} Λ(n)min(1, y/n)κ(n)(n/q)\|² ≤ C y log y \|P\|, all \|κ\| ≤ 1, y ∈ [Y, Y^K], P = primes ≡ a0 (m0) in (√U, U], log U ∈ [Y^{ε/6}, Y^{ε/5}log²Y], Y ≥ Y0(ε, K) | KNOWN: derived in HOSTS §8 from reciprocity, Chebyshev, Landau–Page and PNTAPRelW (sums shorter than √q, conductors s ≤ Y^{800K} = (log U)^{O(K/ε)}); same scheme as moments of L(1, χ) over quadratic families (Granville–Soundararajan 2003) [not re-read]; no large sieve for real characters needed |
| QH.twistK | AA/MeanSquareBounds/twistK.lean:148-151 | ∃ n* ∀ squarefree odd s ≤ Y^A, (s, m0) = 1, s ≠ n*: \|Σ_{q ∈ P}(q/s)\| ≤ Y^{−B}\|P\| | KNOWN: PNTAPRelW (relative form) + CRT + Landau–Page with a common cutoff m0 Z² |
| PNTAPRelW | JCS/MeanSquareHypWIff.lean:38-46; proved JCS/PNTFinal.lean:22 [15] | relative error A0(e^{−c3 log X/log M} + e^{−√(c3 log X)}) to the exceptional-sensitive main term (li X − χ1(b) li X^{β1})/φ(M), log X ≥ C0 log M log log M | KNOWN: Linnik's three principles (log-free density, Deuring–Heilbronn, zero-free region), Gallagher 1970 (Invent. Math. 11, 329–339), Iwaniec–Kowalski 2004 ch. 18 [theorem number of the relative form not re-read] |
| ExcBound | src/Formalism/LandauZeroFreeBounds/LogDerivRegionChiOf.lean:50-53; proved src/Formalism/LandauZeroFreeBounds/SiegelWalfiszAPOfSiegelLower.lean:32 [15] from siegelLower (FinalSiegel.lean:17) [16] | β ≤ 1 − c(ε) q^{−ε} for real χ ≠ 1 | KNOWN: Siegel 1935 (Davenport ch. 21); used with ε = 1 only, weaker than the effective 1 − c/(√q log²q) |
| QH.bad_share | AA/MeanSquareBounds/BadShare.lean:251-254 | #{q: ∃ y ∈ [Y, X], \|S\| > y^{1/2+τ}/2} ≤ η\|P\| from MeanSq | KNOWN: Markov + Gallagher–Sobolev (\|S(y)\|² ≤ 2\|S(Y)\|² + 2Y∫\|∂S\|², ∂S is again an S, :118-148) over dyadic blocks |
| zero_density_box | JCS/CountZAll.lean:13 (def JCS/CountSt.lean:7-11) [17] | ≤ C Q^{1/4} primitive χ, cond ≤ Q, with a zero in Re s > 1 − 1/400, \|Im s\| ≤ Q | KNOWN: Montgomery 1971 (LNM 227, Thm 12.2) gives O(Q^{6/399} log^{14} Q); proof in JCS/CountZ*.lean is Montgomery–Halász (mollifier ⌈Q⌉⁴, polynomials to ≈ Q⁸, duality) |
| count_a | JCS/MSqCount.lean:346-350 [15] | \|PPset(m0, a0, X)\| ≥ X^{9/10} if log X ≥ (1/4)log²x·log m0 | KNOWN: PNTAPRelW + ExcBound(ε = 1) |
| box_of_S, tailK, kappaC | AA/MeanSquareBounds/BoxOfS.lean:28-31, TailKExp.lean:92-95, kappaC.lean:9-11 | Gchi = −S(Z, y) − y·tail_Z; \|Y^K tail\| ≤ Y^ε and \|L'/L(1)\| ≤ C log Y from the box | KNOWN: contour shift into the box with Landau's lemma (LandauLocal3, LogDerivRegionChiOf.lean:57-71) |
| kappa4 | AA/MeanSquareBounds/Kappa4Loglog.lean:126-128 [11] | box + log Q ≤ 3Y^{ε/5}log²Y ⇒ \|L'/L(1, χ)\| ≤ 4 log Y | KNOWN in kind ('zero-free strip ⇒ L'/L(1) ≪ δ^{−1} log log q'); constant ≈ 80ε/(1 − θ) < 4 since ε < 1/52 (reviewer: even δ' = 1/800 gives 160ε < 4) |

## (2b) Counting: all conditions at once (QH.stage, AA/MeanSquareBounds/StagePow.lean:41-46; pair_share PairShare.lean:26-31)
- Stage 1, q ∈ PPset(8, 3, e^z): bad q = (some ξ mod 8M has a zero in the box) ∪ (some ξ has max_{y ∈ [Y, Y^K]}|S| > y^{1/2+τ}/2). First set ≤ C(8M e^z)^{1/4} (distinct q give distinct primitive characters: conductor = cond ξ · q with q > F ≥ cond ξ); second ≤ Σ_ξ (η/2R)|P| (bad_share per ξ, R = number of ξ); |P| ≥ e^{0.9z} (count_a); stage_pow (:29-31): C U^{3/8} ≤ (η/2)U^{9/10}. The continuum y ∈ [Y, Y^A] is inside bad_share (derivative bound per dyadic block). The L'/L(1) bound and the tail need no averaging (kappa4/kappaC/tailK from the box). Box heights: j = 0: N0 e^z = 8M e^z = height(M, ε, Y, 0) exactly.
- Stage 2, for each good q: d ∈ PPset(16·oddM(z, q), a0, e^{z'}), z' = z log²Y (ExistsCoprimePrimePairSetSecondEqPPset.lean:56-58), log m0 ≈ 2z ≤ 4z'/log²Y (admitted), primes of m0 ≤ z or = q > Y^{400K}; Ξ = {ξ·χ_q^a}; box height N0 e^{z+z'} = height(M, ε, Y, 1): the j = 1 height exceeds d, which costs nothing because ZeroDensityBox allows T = Q: bad d ≤ C(8M e^{z+z'})^{1/4} ≪ e^{0.9z'} ≤ |P_d|.
- pair_exists (PairExists.lean:9-12) takes η = 1/2: one q in G1, then one d in G2. Rows ⇒ WinPair (pairEx.lean:139-143). All quantifiers close; Y1 = Y1(M, ε, A), non-uniform and ineffective (Siegel), which is harmless downstream: winPairs uses PairEx at M = ψ.N, A = J for one fixed target ψ and ε (AA/CharacterBoxes/WinPairs.lean:12).

## (3) Circularity (HOSTS.md §7)
- Bounded(Level.cyc ψ).KEPT = twists ψ^k χ_q^i χ_d^j, (i, j) ≠ (0, 0) (boundedTw_key, AA/CharacterBoxes/BoundedTw.lean:5-9); KEPT = ⟨ψ⟩ (ExplicitFormulaErrorOfHostConditions.lean:105-121). The pair (q, d) is chosen by averaging over q and d with ξ = ψ^k as an arbitrary bounded coefficient (meanSqK) and a family-wide zero count (zero_density_box); kappa4/tailK use only the twisted character's own box; ThinH uses de la Vallée Poussin + RealBox + small_zeros. None uses a zero-free region or a bound for ψ^k itself, nor (i, j) = (0, 0). Verdict: no circularity.

## (4) Probes (PROBES.md; code/, logs/)
| Probe | Result |
|---|---|
| (a) 7833 primes d ≤ 10^6 (+194 WinPair-like d ≡ 7 (16), ≡ −1 (105)), ξ = 1 and χ_{−4} | max_y\|ψ(y, ξχ_d)\|/√(y log y) on [10², 10^5]: median ≈ 2.0, max 3.8; margin y^{1/104} ≈ 1.05–1.12: 0 good d (expected; NOT a counterexample, PairEx is ∃ Y1); \|L'/L(1)\| ≤ 0.56 log 10² (Riesz approximation, uncertified) |
| (a3) archive regime Z² < d: 70435 primes d ∈ (10^6, 2·10^6], Z = 10³ | family mean of \|S\|² = exact square-pairing diagonal within 1–7 % (M = 1) and 1 % (M = 4): the meanSq_abs mechanism |
| (a4)/(d) scale | pointwise margin needs log y ≳ 430; Markov step needs Y1 ≳ 10^{540} (M = 1, ε = 1/52); Y1 grows like φ(8M)^{2/ε} and is ineffective |
| (b) certified argument principle (arb balls), 17 real characters, conductors 3–24 (T = 8·cond ≤ 192) and 163 (T = 40) | 0 zeros in Re s > 1 − 1/400 for all; positive control at σ0 = 0.4 finds the known zeros (2, 2, 4); 1009 not certified (time cap) |
| (b2) j = 1 height vs density | Montgomery count e^{0.015 log Q}, archive C Q^{1/4}, family ≥ e^{0.9 log Q (1 − o(1))}: the height beyond d is free |
| (c) ThinH | classical saving x^{2ε/15}/log²x vs x^{ε/8}: holds from log x ≈ 1.7·10^5 (ε = 1/52); a Siegel zero would kill it, so RealBox is necessary for q- or d-divisible real characters and is supplied by PairGood.boxes; others have conductor ≤ 8B; G(y)/√y = 0.048 → 0.009 for y = 10³ → 10^7 (principal case well defined) |

## Reviewer (logs/reviewer.md; 4 rounds, briefs .scratch/alphost/reviewer_task*.md)
- Round 1 GAP(1, 7): my summary had dropped meanSqK's Y0(ε, K) and the principal-character definition of gChi; both RESOLVED in round 2 by the exact Lean statements. Round 2 asked for a classical derivation of MeanSq (now HOSTS.md §8). Round 3 GAP(5): normalization of the PNT error by |P| when an exceptional character of conductor | m0 thins P. Round 4: RESOLVED by the relative form of PNTAPRelW (Σ_b main ≤ 2|P|), PNTAPRelW CLASSICAL (quantitative Deuring–Heilbronn + log-free density give the factor min{1, (1 − β1) log X}); the single n* needs a common cutoff m0Z² and quantitative absorption, not per-modulus uniqueness (wording adopted). Final: KNOWN-CONFIRMED.

## What is and is not shown
- Shown: every analytic input of PairEx and ThinH, read at statement level (bodies not checked, no Lean build), is a classical theorem or a routine consequence, with ranges that fit; the counting closes for all conditions at once; the auxiliary family avoids the target.
- Not shown: that the Lean bodies are correct (kernel check not run, per Adrian's order); exact theorem numbers of the relative PNT form and of Granville–Soundararajan's moment computation [not re-read]; anything about (I1) VWin or (I3) the lattice data, where alpchain located the possible new mathematics.

## Discoveries
None. Nothing for DISCOVERIES.md, RECORDS.md or REFUTATIONS.md: (I2) is classical. Tree consequence (for the orchestrator): close the 'hosts' branch of the Alpöge audit as KNOWN; the audit's remaining load is (I3) LinkDIV_at/DIVinN and (I1) VWin.

## Next checks (ranked)
1. (I3): open closed_linkDIVt_of_closed_bridge_kernel and DivisionKernelConstantBounds (alpchain next check 1).
2. (I1): VWin through automorphy lifting (VVle's Sha-length bound √(exp(x^ε))).
3. Optional: read the bodies of zero_density_box (JCS/CountZ*.lean) and pnt_progressions_relative_weak_holds for axioms/sorry-free status (alpaudit's axiom scan covers this if it ran on the cone).

## Files
- Riemann/Deep/AlpHost/REPORT.md — this report.
- Riemann/Deep/AlpHost/HOSTS.md — the host statements in standard notation (§1–5), PairEx proof tree with file:line and status (§6), target entry (§7), classical derivation of MeanSq (§8).
- Riemann/Deep/AlpHost/PROBES.md — probes (a) twist distributions and the mean-square mechanism, (b) certified zero-free boxes and the j = 1 density count, (c) ThinH exponents and the principal case, (d) scale and uniformity of Y1.
- Riemann/Deep/AlpHost/code/stmts.py — statement extractor (declaration headers up to ':=') for a Lean directory; outputs .scratch/alphost/msb_stmts.txt, jcs_stmts.txt.
- Riemann/Deep/AlpHost/code/probe_a.py (numpy; modes win, sample), probe_a2.py (numpy; Z² < d regime), probe_b.py (python-flint; argument principle; args D-list, T cap, σ0), probe_cd.py (numpy; exponent tables, Y1 scale, G(y)).
- Riemann/Deep/AlpHost/logs/probe_a_win.log, probe_a_sample.log, probe_a2.log, probe_b_control.log, probe_b.log, probe_b2.log, probe_cd.log (memcap peaks ≤ 41 MB), reviewer.md (4 rounds), reviewer_round1_raw.md.
- .scratch/alphost/STATE.md, reviewer_task.md … reviewer_task4.md (reviewer briefs), msb_stmts.txt, jcs_stmts.txt (extracted statements), flint_api.py, fix_tables.py (helpers).

