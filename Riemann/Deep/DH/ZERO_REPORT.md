# DH off-line zero — REPORT

**Result: done at the full 1e-5 precision, no sorry, no native_decide.** `lake build Riemann.Deep.DH.ZeroFinal` green.
`Fayed.DH.Zero.offLineZero_holds : Fayed.DH.OffLineZeroClaim` and `Fayed.DH.Zero.dh_offline_holds : ∃ s, DH s = 0 ∧ s.re ≠ 1/2`
(via `Fayed.R81.dh_offline`); `#print axioms` = [propext, Classical.choice, Quot.sound]. The zero lies within 5e-6 of z_A = 0.8085171825 + 85.6993484854 i (mpmath: 0.80851718246 + 85.69934848538 i).

## Route (not a, b or c: Newton–Kantorovich + Banach fixed point)
1. **Continuation** (ZeroRepr): DH = Σ_{m<M} blocks + 5^{-s} Σ_j c_j [A(M+j/5) + s·J_j(s)] on Re s > 0, s ≠ 1. Equal for Re s > 1 by per-class telescoping; identity theorem on {Re s > 0, s ≠ 1}. J_j (shifted sawtooth integral) is analytic as a Mellin transform of a bounded function (`mellin_differentiableAt_of_isBigO_rpow`), so no parametric-integral work.
2. **Order-8 Euler–Maclaurin** (ZeroEM): shifted copy of FirstZero/EM.lean; ‖DH − Gem M‖ ≤ Ebd, Gem in closed form (n^{-s} and an explicit rational Qf).
3. **Newton lemma** (ZeroNewton): F entire, ‖F z0‖ ≤ e0, difference quotient over a step k within ed of d, ‖F''‖ ≤ M2 ⇒ z ↦ z − F z/d contracts the closed r-ball ⇒ a zero. Cauchy estimates give F'' from sup |F| and the difference quotient of the Euler–Maclaurin error.
4. **Certificate** (M = 60, 244 values n^{-s} per point, 2 points z_A and z_A + 1e-5 i): logs by log(1−x) series, n^{-σ} = exp((1−σ)log n/2)²/n, e^{iθ} degree 29, products for composites. Ball radius 1.8e-16; Euler–Maclaurin bound 3.6e-6 at z_A (true error 4.6e-9).
5. Numbers: |DH(z_A)| ≤ 3.62e-6, |d| ≥ 1.2323, ed ≤ 9.7e-5, |DH| ≤ 157.6 on radius 0.74 ⇒ |DH''| ≤ 1151, q = 0.0141, r = 5e-6 (box needs < 9.5e-6).

## Files (Riemann/Deep/DH/)
| file | content |
|---|---|
| ZeroRepr.lean | sawtooth integral, Mellin analyticity, telescoping, DH_eq_Grep |
| ZeroEM.lean | J_EM8, DH_sub_Gem_le, Gem_eval |
| ZeroNewton.lean | deriv2_le, diff_quot_le, newton_zero (pure Mathlib) |
| ZeroLib.lean / ZeroBounds.lean | log/exp/cis/W/Gmain balls; κ enclosure, Gem_ball, Ebd_le, crude bounds, Gem analytic |
| ZeroCertLog/A/B.lean | generated (.scratch/dhzero/gen.py full); ~10 s each |
| ZeroFinal.lean | assembly (.scratch/dhzero/final.py + genfinal.py) |

## Notes
- Not imported from Riemann.lean (outside edit scope): add `import Riemann.Deep.DH.ZeroFinal`. Not committed. Whole chain builds in about 3 min.
