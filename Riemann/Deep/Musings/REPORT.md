NO DISCOVERY: Nadim (arXiv 2605.30767) claims no proof of RH; its one real lemma (§3, contours from the critical line never meet off it) is correct (proved in Lean); the scout lead is false; only small errors in §2 and §5.

# Audit: Nadim, "Musings on the Riemann Hypothesis" (arXiv 2605.30767, May 2026)

**Verdict: there is no failing step, because the paper proves nothing about RH and says so.** Its conclusion: "We have presented a few basic ideas in the hope that they might inspire new approaches". The program stops at an unknown conformal map (§4) and a sufficient condition the author says "do[es] not appear to hold" (§5, §6).

## Chain (z = s − 1/2, ξ = f + ig)
1. §2: f is even and g is odd in each of x and y; a root in one quadrant gives roots in all four. **Correct.** Typo: "Im z⁴ = x³y − xy³" should be 4x³y − 4xy³ (`paper_z4_false`).
2. §3 lemma: "the blue and red curves that emanate horizontally from the vertical axis can never intersect with each other or with themselves". The reason given: they would enclose a region with f = 0 or ∂f/∂n = 0 on its boundary, so f ≡ 0. **Correct.** Lean proves a sharper form: Im ξ² = 2fg is harmonic and vanishes on the boundary, so it vanishes inside, so ξ² is constant (`no_enclosed_region`, `sec3_chain`). The only unformalized input is the picture fact that two such curves plus a piece of the critical line bound a region (the `Encloses` Prop).
3. §4 half-strip model: f = Σ bₙ sin(ny)cosh(nx). Correct for finite sums. The examples with off-axis zeros (b₂ = −1/4; sin y + 0.2 sin 4y) all involve new contours that do not start on the axis, which matches §3.
4. §5 Taylor criterion: if (−1)ᵏF⁽²ᵏ⁾(y) ≥ 0 for all k, then f(x,y) ≥ F(y) for all x. **Correct** (`taylor_criterion`). **But it never applies to ξ:** as x → ∞, arg Γ(s/2) grows like (y/2)·log x, so Re ξ(1/2+x+iy) changes sign infinitely often on every horizontal line y ≠ 0. Numerically there are 7, 8, 13 and 16 sign changes for x in [0,100] at y = 17.5, 20, 48.9, 81 (`signs.py`, not rigorous). So the criterion fails at every y where F(y) ≠ 0 (`criterion_fails`), and "if F(y) had this property everywhere, the Riemann Hypothesis would hold" says nothing. The paper only says the criterion is "not clear" / "do[es] not appear to hold".
5. §5 error: the "different harmonic function" Σ sin((2k+1)y)cosh((2k+1)x)/(2k+1)³ diverges for every x > 0, because its terms are unbounded at y = π/2 (`sec5_series_not_summable`). So that boundary problem has no solution at all; the paper only doubts "finite truncation".

## Scout lead: refuted
The lead said §2 claims Re = 0 and Im = 0 curves cannot cross off the line, and the §3 example contradicts it. In fact the claim is in §3, it covers only curves that start on the axis, and it is a theorem. The off-axis crossings in the paper's §4 examples come from curves that do not start on the axis, which the paper itself points out.

## Lean (`Riemann/Deep/Musings/Paper.lean`: builds, sorry-free, axioms propext/choice/Quot.sound)
`im_z4`, `paper_z4_false`, `im_sq_zero_on_closure` (maximum principle applied to exp(±iξ²)), `setOf_im_eq_zero_not_mem_nhds`, `no_enclosed_region` (uses the open mapping and identity theorems), `Encloses`, `sec3_chain`, `taylor_criterion`, `criterion_fails`, `sec5_series_not_summable`. Axiom check: `.scratch/musings/Axioms.lean`. Not yet imported in `Riemann.lean` (orchestrator's job).

## Novelty (web_search; no contact with the author)
- No critique of 2605.30767 was found; only listings (arXiv, alphaXiv, arxiv.gg, papers.cool).
- Mapping the zero curves of Re ζ and Im ζ is old work (Arias de Reyna, "X-ray of Riemann's zeta-function", 2003).
- The Im ξ² argument and the Γ-argument growth are standard. Nothing here is significant.

## Files
`.scratch/musings/main.tex` (source), `signs.py` (sign changes of Re ξ), `Axioms.lean`.
