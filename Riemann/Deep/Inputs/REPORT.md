# Fayed.Inputs: all six fields are proved

Build: `lake build Riemann.Deep.Inputs.Perron` (green). Every theorem below uses only the axioms [propext, Classical.choice, Quot.sound]. No sorry, no new axioms, no named hypotheses left.

| field | status | Lean name | file |
|-------|--------|-----------|------|
| leibniz44 | proved | `Fayed.InputsProof.leibniz44_holds` | Leibniz.lean |
| em17 | proved | `Fayed.InputsProof.em17_holds` | EM.lean |
| convexity | proved (before this work) | `Fayed.Convexity.convexity_holds` | ../Convexity/Convexity.lean |
| strip, reflect | proved (before this work) | `Fayed.strip_holds`, `Fayed.reflect_holds` | ../../Fayed/StripReflect.lean |
| perron33 | proved | `Fayed.InputsProof.perron33_holds` | Perron.lean (+ PerronZeta, PerronKernel, PerronShift) |

Main results:
- All.lean: `inputs_iff_perron33 : Inputs ↔ Perron33`, `RH_of_theorem1C_perron33`.
- Perron.lean: `perron33_holds : Perron33`, `inputs_holds : Inputs`, `RH_of_theorem1C : Theorem1C → RiemannHypothesis`.

So Fayed's reduction is now unconditional on its literature inputs. RH still does not follow: Theorem 1C (the L'Hôpital rule on a disconnected domain) is false (Riemann/Fayed/Theorem1.lean); that is the failing step of the paper.

Note: `Riemann.lean` does not import `Riemann.Deep.Inputs.Perron` yet (outside this task's edit scope), so a plain `lake build` does not build it; add `import Riemann.Deep.Inputs.Perron` there.

How the proofs go:
- em17: at N = ⌊u⌋ the existing order-8 Euler–Maclaurin formula `zeta_EM8` (FirstZero/EM.lean) gives S(N). On [N,u], x^{1-s}/(1-s) − EMcorr(x) has derivative s(s+1)(s+2)B₃(x−N)/6·x^{−s−3}. The mean value inequality then gives O(N^{−σ−3}) ≤ 2^{σ+3}·O(u^{−σ−3}).
- leibniz44: substitute u = eʳ, so the integrand A(y)e^{rW(y)}ζ(W(y)) is jointly continuous on ℝ². P(r,w) = ∫₀^w has continuous partial derivatives (dominated differentiation, the fundamental theorem of calculus, continuity of parametric primitives), so Mathlib's `hasStrictFDerivAt_uncurry_coprod` applies. Then compose with r ↦ (r, ±e^{8r}) and with log.
- perron33 (≈920 lines, four files):
  1. PerronZeta.lean, `norm_zeta_le`: ‖ζ(σ+it)‖ ≤ 13(1+log|t|)|t|^{max(1−σ,0)} for 1/2 ≤ σ ≤ 3, |t| ≥ 1, from the first-order formula `Zeta0EqZeta` (ported PNT+, WeilEF/Z23/FromPNTPlus/ZetaBounds.lean) at N = ⌊|t|⌋. This replaces the uniform convexity bound: it is weaker, but enough.
  2. PerronKernel.lean, `perron_kernel` (Apostol 11.12 with explicit constant): ‖(1/2π)∫_{−T}^{T} e^{L(c+iy)}/(c+iy) dy − [L>0]‖ ≤ 2e^{cL}/(π|L|T). L > 0: rectangle [−T,c]×[−T,T], residue 1 at 0 (`ResidueTheoremOnRectangleWithSimplePole'`); L < 0: rectangle [c,c+T]×[−T,T], Cauchy (`HolomorphicOn.vanishesOnRectangle`). Far side bounded by e^{−|L|T} ≤ 1/(|L|T), so no limit R → ∞ is needed.
  3. PerronShift.lean, `shift_identity`: rectangle [1−σ−η, 2]×[−T,T] for ζ(s+z)u^z/z; the only pole inside is z = 1−s (residue u^{1−s}/(1−s), via `riemannZeta_residue_one` + `ResidueOfTendsTo`); z = 0 is outside since 1−σ−η > 0. `horiz_bound`: with T = u⁸ each horizontal side is ≤ 52(1+log 2u⁸)u^E/u⁸, E = max(1−σ+7η, 2) < 4 (because 7η < 3+σ when η < min(σ,1−σ)).
  4. Perron.lean: on Re z = 2, ζ(s+z)u^z/z = Σ(n+1)^{−s}e^{L_n z}/z, L_n = log u − log(n+1); sum and integral swap (`integral_tsum_of_summable_integral_norm`). On D_γ, |u−n| ≥ γ for all integers n, so 1/|L_n| ≤ (u+n+1)/γ and the Perron error sum is ≤ 4u³ζ(1+σ)/(πγT) = O(u^{−5}). Total error ≤ C·u^{−4+ε} with C, U independent of M₀.
