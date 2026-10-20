# Fayed.Inputs: which fields are proved

Build: `lake build Riemann.Deep.Inputs.All` (green). Every theorem below uses only the axioms [propext, Classical.choice, Quot.sound]. No sorry.

| field | status | Lean name | file |
|-------|--------|-----------|------|
| leibniz44 | proved | `Fayed.InputsProof.leibniz44_holds` | Leibniz.lean |
| em17 | proved | `Fayed.InputsProof.em17_holds` | EM.lean |
| convexity | proved (before this work) | `Fayed.Convexity.convexity_holds` | ../Convexity/Convexity.lean |
| strip, reflect | proved (before this work) | `Fayed.strip_holds`, `Fayed.reflect_holds` | ../../Fayed/StripReflect.lean |
| perron33 | open, kept as a named hypothesis | `Fayed.InputsProof.Perron33` | All.lean |

Main results (All.lean):
- `inputs_of_perron33 : Perron33 → Inputs`, and `inputs_iff_perron33 : Inputs ↔ Perron33`
- `RH_of_theorem1C_perron33 : Perron33 → Theorem1C → RiemannHypothesis`

The unconditional `inputs_holds` / `RH_of_theorem1C` are NOT reached, because perron33 is still open.

How the proofs go:
- em17: at N = ⌊u⌋ the existing order-8 Euler–Maclaurin formula `zeta_EM8` (FirstZero/EM.lean) gives S(N). On [N,u], x^{1-s}/(1-s) − EMcorr(x) has derivative s(s+1)(s+2)B₃(x−N)/6·x^{−s−3}. The mean value inequality then gives O(N^{−σ−3}) ≤ 2^{σ+3}·O(u^{−σ−3}).
- leibniz44: substitute u = eʳ, so the integrand A(y)e^{rW(y)}ζ(W(y)) is jointly continuous on ℝ². P(r,w) = ∫₀^w has continuous partial derivatives (dominated differentiation, the fundamental theorem of calculus, continuity of parametric primitives), so Mathlib's `hasStrictFDerivAt_uncurry_coprod` applies. Then compose with r ↦ (r, ±e^{8r}) and with log.

Plan for perron33 (not done; estimated 1500+ lines):
1. Rectangle [1−σ−η, 2]×[−T,T] with T = u⁸. The only pole is z = 1−s, with residue u^{1−s}/(1−s); z = 0 lies outside the rectangle. Residue machinery on rectangles exists in .scratch/fm/.../zeta23/FromPNTPlus/ResidueCalcOnRectangles.lean and still has to be ported.
2. Horizontal sides: this needs |ζ| ≤ C·T^{η/2+ε} uniformly for Re w ∈ [1−η, 3]. The MaxMod/Bootstrap work gives it for Re < 1; Re ∈ [1, 1+δ] still needs a uniform version.
3. On Re z = 2, expand ζ as its Dirichlet series and use the truncated Perron bound |P(y) − 1_{y>1}| ≤ y²/(πT|log y|).
4. Crude bound for the sum of errors: if u ∈ Dγ then |u−n| ≥ γ, so 1/|log(u/n)| ≤ max(u,n)/γ. The total error is then O(u³/T) = O(u^{−5}), with no need to split the sum.
