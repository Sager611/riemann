/-!
# Round 531 — node N5 (Nyman–Beurling/Báez-Duarte), round 24, falsification of nbdir's COROLLARY (the continuous Nyman–Beurling rate for primitive Dirichlet L-functions under GRH): the constant, the trivial zero s = 0 of even χ, the central zero, novelty and new numerics all hold; three explorations, no transfer

* Node: N5 (24th round; R520 before). Round type: falsification of our newest result's deltas from its audited predecessor (R501 pattern): nbdir (COROLLARIES.md, 20 Nov 2026, routine, three referees CORRECT) transfers the audited ζ theorem (DbarNB) to L(s,χ).
* (1) cconst: from the Hadamard product of ξ(s,χ) and Re B(χ) = −Σ Re(1/ρ): Σ m_ρ Re(1/ρ) = ½log(q/π) + ½ψ((1+a)/2) + Re L′/L(1,χ), and Re(1/ρ) = 1/(2|ρ|²) on GRH, so C_χ is as printed; 30 digits: C_{χ₋₄} = 0.1555679799235859 (closed form of L′(1,χ₋₄)), C_{χ₅} = 0.1565569527, C_{χ₈} = 0.2354315618. HOLDS.
* (2) trivial: the residue integrand of (c2) is F₀ = H/L (DbarNB §5 (iv)–(v)), not H/(sL): the factor 1/s belongs to the Mellin transform L(s,χ)/s of the target and to the orthogonality step, so at s₀ = 0 for even χ the pole is simple, r₀ = H(0)/L′(0,χ), and ‖L R_tr‖² ≪ e^{−(2a+1)(√L−2)} stands; s = 0 is excluded from the sums of §5 and §7. HOLDS.
* (3) central: ρ = ½ enters once (self-reflected, counted once in the Z-groups), with weight m²/|ρ|² = 4m² in both bounds; the equivalence 'limit = C_χ ⟺ GRH + simple zeros' survives a central zero (4m² against 4m). HOLDS.
* (4) numer (scout timed out; orchestrator ran nbdir_D.py from a copy with its output paths moved to /tmp, peak 40 MB): real even χ mod 8, NT = 16384, h = 0.0015625: L·D² = 0.2763, 0.2098, 0.1735, 0.1681, 0.1771, 0.1767, 0.1859, 0.1952 at L = 2, 3, 4, 5, 6, 8, 10, 12, below C₈ = 0.2354 and increasing from L = 8 (as for χ₅: 0.1510 at L = 16 against C₅ = 0.1566; the dip at small L is shared). HOLDS.
* (5) sota: no published continuous rate or GRH upper bound of this order for Dirichlet L or the Selberg class (closest: de Roton 2006 Thm 2, lower bound; arXiv 1608.07887, 1101.1199, 1705.09918). HOLDS (not found).
* Explorations: Jensen–Carleman as a stability estimate (inverse problems, N8): weighted zero counting only, blind to Re ρ (r = 15 reproduces the ±14.1347i pair); Nelson–Aalen and Aalen's intensity model (survival analysis, N9): risk-set assumptions fail for deterministic zero counts; p-adic Weierstrass preparation (Dwork, N12): no correspondence with complex zeros. No transfer.
* Verdict: every probe of nbdir's deltas holds; the COROLLARY stands as recorded.
-/

namespace Fayed.R531
end Fayed.R531
