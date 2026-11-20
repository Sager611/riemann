# AlpHost — reviewer rounds (independent 'reviewer' subagent; briefs in .scratch/alphost/reviewer_task*.md)

## Round 1 (statements, counting, references, circularity pasted inline)
- 1 MeanSq: DISAGREE as stated — missing large-Y threshold (a resonating κ = χ_{q0} gives |S(q0)|² ≈ y² log²Z). FIX: the Lean statement has it (QH.meanSqK: ∃ C ∀ ε K ∃ Y0 ∀ Y ≥ Y0); with |P| ≥ U^{9/10} the resonance costs y² log²Z/|P| ≤ Y^{2K}(400K log Y)² e^{−0.9Y^{ε/6}}.
- 2 PNTAPRelW: AGREE (Gallagher + Deuring–Heilbronn; zeros outside a narrow exceptional region absorbed by smaller c3).
- 3 bad_share: AGREE (given MeanSq, Y large in τ). 4 ZeroDensityBox: AGREE (Q^{6/399} log^{14} Q ≪ Q^{1/4}).
- 5 kappa4: AGREE asymptotically (even δ' = 1/800 gives (160ε + o(1)) log Y < 4 log Y); the tail needs a separate shift near 1 − 1/400 (tailK does this).
- 6 counting: AGREE. 7 ThinH: DISAGREE — principal character (pole of L'/L(1), ψ ~ u). FIX: gChi(1, y) := G(y) = Σ Λ(n)(y/n − 1) − y(log y − 1 − γ) (Analytic.lean:14-24) and psi0(1, u) = ψ(u) − u (src/Formalism/PeriodIntegrals/ExplicitFormula.lean:9-11): main term subtracted; G(y) = −y∫_y^∞ (ψ(u) − u)u^{−2}du.
- 8 circularity: AGREE. 9 numerics: AGREE (interpretation). Round-1 verdict: GAP(1, 7).
## Round 2 (exact Lean statements for items 1 and 7)
- Item 1 RESOLVED; item 7 RESOLVED (sign and constant of the tail identity correct; 2ε/15 > ε/8 rate works). New request: a derivation of MeanSq from classical theorems. Verdict: GAP(classical provenance of MeanSq).
## Round 3 (derivation HOSTS.md §8)
- Step (5) fails as written: Gallagher's error is O(Li(U)/φ(m0)·E), and an exceptional χ* of conductor | m0 with χ*(a0) = 1 can make |P| ≪ Li(U)/φ(m0). Verdict: GAP(5).
## Round 4 (the RELATIVE form of PNTAPRelW, JCS/MeanSquareHypWIff.lean:38-46)
- (i) RESOLVED: π(U; M, b) ≥ (1 − ρ)main(e, b, U) gives Σ_b main ≤ 2|P|, so Σ_b|π − main| ≤ 2ρ|P|; an exceptional character induced from m0 cancels in the twist. (ii) PNTAPRelW CLASSICAL: quantitative Deuring–Heilbronn + log-free density supply the factor min{1, (1 − β1) log X}. Note: the single n* needs a common cutoff Q = m0 Z² and quantitative absorption of other exceptional terms, not uniqueness per modulus (wording adopted in REPORT).
- VERDICT: KNOWN-CONFIRMED.
