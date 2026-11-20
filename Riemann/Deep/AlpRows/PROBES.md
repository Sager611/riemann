# PROBES — finite checks of the HypsAlt rows (alprows, 20 Nov 2026)

All code in code/, logs in logs/; exact rationals (fractions/sympy) unless stated; memcap peaks in the logs. Archive paths relative to .scratch/alp/lean/src/.

## P1. hP = StandardLattice.SiegelMassRatioCM: normalization of the right side (code/p1_mass.py, logs/p1_mass.log, peak 14 MB)
- Row (StandardLattice.lean:281-285; RHS StatementMass.lean:21-43): mass(I_M)/mass(I_{M−1}) = massRatioRightSide K M for IsCMUnram K, M ≥ 3.
- Instances: IsCMUnram biquadratic fields (P3): K = Q(ζ12) (F = Q(√3)), Q(√−3,√−7) (Q(√21)), Q(√−2,√−3) (Q(√6)), Q(i,√−7) (Q(√7)), Q(√−3,√−5) (Q(√15)); M = 3..8. (First run had 3 fields; the reviewer noted Q(i,√−7) is smaller than Q(√−3,√−5); rerun with 5.)
- Lean side computed literally: Euler product over the prime ideals P of O_K (splitting type from the three quadratic subcharacters; P̄ = P iff complex conjugation ∈ decomposition group), eulerFac as defined, primes < 2·10⁵ (M = 3) or < 2·10⁴, vInv = |disc K|^{1/4}/2, d = 2.
- Independent side: |L_F(1−M, χ^M)|/2^d, L_F(s,χ) = L(s,χ_a)L(s,χ_b) (M odd), ζ(s)L(s,χ_F) (M even), from generalized Bernoulli numbers B_{M,χ}.
- Result: 30/30 agree to ≤ 7.6e−15 (truncation of the Euler product): e.g. Q(ζ12): 1/36, 23/240, 5/12, 1681/504, 427/12, 257543/480 for M = 3..8. So the Lean RHS is exactly the consecutive-ratio of the Gan–Hanke–Yu / Prasad mass formula for unimodular hermitian lattices with hyperspecial stabilizers everywhere. PASS.

## P2. hP: consistency with automorphism groups (same script)
- Aut(I_M) over Z[ζ12] is the monomial group μ12^M ⋊ S_M (order 12^M M!): norm-one vectors are the 12M coordinate vectors (a totally positive α ∈ Z[√3] has Tr α ≥ 2, so Σ x_i x̄_i = 1 forces one x_i x̄_i = 1, and x x̄ = 1 ⇒ x ∈ μ12).
- mass(I_1) = 1/12 directly (h_K = 1; 2+√3 = N(1+ζ12), so the genus of ⟨1⟩ has one class), = the GHY value 2|L_F(0,χ)|/4 = 2·(1/6)/4.
- GHY masses 1/12, 1/288, 1/10368 = 1/#Aut(I_n) for n = 1, 2, 3 ⇒ class number one up to n = 3, and the Lean ratio R(3) = 1/36 = #Aut(I_2)/#Aut(I_3). A normalization too small by a factor 2 (e.g. /2^{d+1}) would give mass(I_3) < 1/#Aut(I_3), impossible (a factor too large is not excluded by this lower bound; reviewer). PASS (consistency; independence limited to n = 1 and the Aut count).

## P3. IsCMUnram (StatementMass.lean:44-49) at the brief's fields
- Q(i): θ − θ̄ = 2b·i ∈ (1+i) for all θ = a+bi; Q(√−3), Q(√−7): θ − θ̄ = b√d_K ∈ (√d_K). So IsCMUnram fails for every imaginary quadratic field (the ramified prime of K/Q is ramified in K/K⁺ = K/Q): the mass row hP, and through the host (D.host.cmStar, cmUnram) the GenusMassAt/ClassNumberFormula₂ rows, cannot be instantiated over Q(i), Q(√−3), Q(√−7). Check of the criterion |d_a d_b| = d_F (⇔ |disc K| = d_F²): true for the three fields of P1, false for the control Q(i,√5) (80 ≠ 5). PASS (definitional).


## P4. LinkDIV_at error bookkeeping: the 1/m + m/x balance (code/p2_rows.py, logs/p2_rows.log)
- Row: DIV conjunct of LinkDIV_at (Formalism/SharpRestriction/GlobalTheoremInputs.lean:83-85) with MAIN u (Formalism/LocalRingCardinality/DataParams₂.lean:232-234), B0 = sizeBudget = f·m·b·x^ε, Nc = budgetNc (DataArith.lean:200-206), SizesTextAt sizes (n ≍ x^{2+8ε}, k = x^{9/2+10ε}, |TT| ≥ x^{−ε/8}n/2, #T2 ≥ x^{−ε}n, a = xb, b = exp(x^ε); f = R.f = 2n (Formalism/LocalRingCardinality/Params.lean:30), taken as f = n (constant factors do not move the minimizer), Cabs = 2 as used by link_gl; constants C set to 1, so P4 checks exponents only, not the actual constants or the producer linkage).
- Relative error (single-lift C|TT|aLa + per-lift C m(B0+Nc)x^{Cabs ε}) / (m|TT| a log a), m = x^μ: minimized at μ = 0.50 for x = 10²⁴, 10⁴⁸ (ε = 0.01, 0.02; μ = 0.55 at x = 10¹² from log factors), minimum ≈ x^{−1/2+O(ε)} (log10 = −11.9 at x = 10²⁴). So the rows reproduce the remark's m ≍ √x: the single-lift term is the 1/m, the per-lift budget B0 ∝ m makes the m/x. PASS.

## P5. idxColLoc and LevelIndexAtPlace (finite ingredients of ClassNumberFormula₂; same script)
- idxColLoc g k = #GL_k(Z/g)/#{unitriangular} (GapAtLeast.lean:56-69) = ∏_{i≤k}(g^i − 1); LevelIndexAtPlace (StandardLattice.lean:446-451) #GL_k(F_q) = #Stab(φ)·(q^k − 1). Brute force over all k×k matrices, q ∈ {2,3,5}, k ≤ 3 (q = 5: k ≤ 2): 8/8 agree. PASS (KNOWN, orbit–stabilizer).

## P6. ClassNumberFormula₂ (hn4; StandardLattice.lean:453-457): integrality of #cls = massI(k)·idxPi(g0,g1)·(T^k−1)^f
- Host E = Q(ζ12), F = Q(√3), f = 2; masses from P2 (1/12, 1/288, 1/10368, 23/2488320, 23/5971968, 38663/3009871872 for k = 1..6); g0 < g1, T primes split in F (≡ ±1 mod 12), GBig (2fk < g_i − 1), k = 3..6: 1548 cases, 0 non-integral right sides. PASS (weak: integrality only).
- Open reading point (not a discrepancy): SplitsF asks only that g split in the totally real F; for g ≡ 11 mod 12 (split in F, inert in E/F) the local unitary group is U_k(F_{g²}/F_g), whose order differs from GL_k(F_g), while idxPi uses the GL_k index. Whether the level structure of (D.mem u).C.cls 0 makes this right was not opened.

## Not probed (no finite content reachable without the realisation data)
- DivRows3G (theta-lattice ties), SmallLevelRankEq, LinkS, LinkIRR, LinkEta: identities/inequalities between fields of the family data D = Fm.D x i (D7) and its indefinite realisation; no small instance can be built without constructing D7 objects. LinkG_at's inert gain demand Σ_t vonMangoldtClassSum ≈ a·p/(p²−1) (wPart = 0) versus the theta slope a/(p²−1): see AlpToy (not repeated).


## P7. LinkG_at inert gain row: internal feasibility and the partner weight it needs (code/p7_gain.py, logs/p7_gain.log)
- Row (GlobalTheoremInputs.lean:129-134) per lift j and t at an inert place over p, with vonMangoldtClassSum exactly as defined (DataParams₂.lean:240-243: i < Nat.log p a + 2, odd i, (a/p^i − wPart)^+).
- Demand with wPart = 0 (≈ a·p/(p²−1)) vs the per-lift cap 4a/p (|TT_w| = |TT|): ratio ≤ 1/3 for a ∈ {10³, 10⁶, 10⁹}, p ≤ 29 (max at p = 2). The row's own bounds are compatible. PASS.
- If the gain per lift is the theta-data slope a/(p²−1) (AlpToy, Lean's rho), the row needs partner weights Σ_{i odd} wPart ≥ a/(p+1) (e.g. p = 3, a = 10⁶: 250000 = a/(p+1)), i.e. the weights w must absorb exactly the gam − rho gap AlpToy found. Not a discrepancy (wPart is a free field of the data, tied to w(partner) in ASK/MAIN), but it locates where the note's a/p + a/p³ + ⋯ and the integral theta slope a/p² + ⋯ are reconciled: in the partner weights of the gain row.


## P8. Win_at (the window row that carries m ≥ x^{1/2−ε} into RestSize; QuadraticLogDerivCongr.lean:231-246) at the SizesTextAt sizes (code/p8_win.py, logs/p8_win.log)
- Conjuncts with known sizes (m-window, x·m ≤ k, k ≤ n²x^{1/2−ε}, k·b' ≤ a ⇔ k ≤ X', 2f(x+1)e^{−c6 x^{ε/8}} ≤ 1/x, f = 2n): all hold from log10 x = 3626 (ε = 0.01, c6 = 1), 6213 (ε = 0.01, c6 = 10⁻³), 2061 (ε = 1/60, c6 = 1), 3618 (ε = 1/60, c6 = 10⁻³); at x = 10¹⁰ the k ≤ X' and f-exp conjuncts fail. Consistent (WinF only asks eventual truth); the route's windows are empty below x ≈ 10²⁰⁰⁰. PASS.
- Newly opened NEW rows without a finite instance (need a D7/TC datum): OpenParameterStatements rows, RowsU, SecondTheoremReduced/Loc, FamGLmin fields, JetMoments, laws_true_fam. Not probed. InducedSimpleLowestAlcove is KNOWN; its condition λ_i − λ_j + (j−i) ≤ p is exactly the closed bottom alcove (N = 2: Sym^d simple iff d ≤ p−1).

