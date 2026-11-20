NG = open('NOGO.md').read(); T = open('TREE.md').read()
line = ("- explore59 (Orca, exploration: complex analysis: Hörmander's L² method for ∂̄ — pole removal with the division isometry v ↦ v/ζ → N5, continuous Nyman–Beurling distance D(λ)) IDEA, not a proof: "
 "explicit Mittag-Leffler corrections with band-limited kernels (∂̄ as the heuristic) reduce, under RH and a standard contour lemma, D(λ)² ≤ (1 + L^{−1/2})·U_{Z(L⁶)}(L − √L) + O(L^{−3/2}) to the cost U of interpolating 1 at the zeros below L⁶ in the one-sided Ornstein–Uhlenbeck Paley–Wiener space (kernel K(γ,γ) = (¼+γ²)L + 1), with no 1/ζ′(ρ); "
 "an m-fold zero costs m²/(L|ρ|²) at the spectral edge (Claim F: Legendre P_k(−1)² = 1, Burnol's m² from above); RH + tail domination (Lemma R′, open) ⇒ lim D(λ)² log(1/λ) = Σ m_ρ²/|ρ|², RH + sup_L L·U < ∞ (Lemma R″, open) ⇒ D(λ)² ≪ 1/log(1/λ) (was (log log)^{5/2+ε}/√log type); "
 "referees: GAP (∂̄ transfer, spectral membership; repaired 6 of 7), then GAP in the high-zero contour isolation and line estimate (repairable by the divided difference (ζ(s) − ζ(w))/(s − w)), ROUTINE as a reduction; RH-blind (Davenport–Heilbronn's on-line part gives the same 0.1836); "
 "no effect on d_N, whose remaining step is explore38's integer-dilation penalty = o(1/log N); proof agent dbarnb (R506) [ExploreDbar](Riemann/Deep/ExploreDbar/REPORT.md)")
t_old = "matches our Vasyunin entries to 10⁻¹⁶ (R461)."
t_new = t_old + (" Under RH, explore59 (20 Nov 2026, Orca exploration, IDEA; [ExploreDbar](Riemann/Deep/ExploreDbar/REPORT.md)) bounds D(λ)² ≤ (1 + L^{−1/2})·U_{Z(L⁶)}(L − √L) + O(L^{−3/2}), L = log(1/λ), by the cost U of interpolating 1 at the zeros below L⁶ in the one-sided Ornstein–Uhlenbeck Paley–Wiener space (Mittag-Leffler pole removal with band-limited kernels; no hypothesis on ζ′(ρ); an m-fold zero costs m²/(L|ρ|²) at the band edge), "
 "so RH + an open uniform interpolation lemma (R″: sup_L L·U < ∞; R′: tail domination) would give D(λ)² ≪ 1/log(1/λ) and lim D(λ)² log(1/λ) = Σ m_ρ²/|ρ|² (BBLS–Burnol's conjecture for the continuous problem; called conjectural in Burnol math/0112254, R506); "
 "two referees: GAP, then GAP only in the high-zero contour isolation and line estimate (repairable), ROUTINE as a reduction; toy Dirichlet polynomials approach Σ m²/|ρ|² from below (R506); proof agent dbarnb running.")
assert 'explore59 (Orca' not in NG and T.count(t_old) == 1
if not NG.endswith('\n'): NG += '\n'
open('NOGO.md', 'w').write(NG + line + '\n'); open('TREE.md', 'w').write(T.replace(t_old, t_new, 1)); print('ok')
