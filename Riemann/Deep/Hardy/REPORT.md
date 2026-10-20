- 2026 start: Kernel.lean (theta kernel Hk, strip, evenness, bounds)
- Repr.lean green: Λ(1/2+it) = (1/2)·𝓕Hk(t/4π), Λ even & real on critical line
- Shift.lean green: contour shift, Λ(1/2+it)e^{αt} = (1/2)𝓕(Hk(·-2αi)), exponential decay
- Moments.lean green: J m α = 2π(-2i)^m Hk^(m)(-2αi) for |α|<π/4 (Fourier inversion + diff under ∫)
- Limit.lean green: Hk^(m) = Rk^(m) - Ck^(m), Ck^(2n) = 16^-n Ck, cusp bound for Rk, Rk^(m)(-2αi) → 0 as α→π/4⁻
- Hardy.lean green: Fayed.Hardy.hardy proved, sorry-free

## Final summary (DONE)
Theorem: Fayed.Hardy.hardy : Set.Infinite {t : ℝ | riemannZeta ((1/2 : ℂ) + t * Complex.I) = 0}
#print axioms: [propext, Classical.choice, Quot.sound]. No sorry/admit/axiom in Riemann/Deep/Hardy/.
Build: lake build Riemann.Deep.Hardy.Hardy green (≈2 min for the last module). 1254 lines, 6 files.

Route (Hardy 1914 / Titchmarsh §10.2, kernel normalisation Hk(u) = e^{-u/4}θ(e^{-u}) - 2cosh(u/4)):
1. Kernel.lean — θ(y) = jacobiTheta₂ 0 (iy); Hk holomorphic on |Im u| < π/2, even (theta FE),
   real on ℝ, ‖Hk u‖ ≤ C e^{-|Re u|/4} on closed sub-strips.
2. Repr.lean — Λ(s) = ½·mellin(θ - 1 - x^{-1/2})(s/2) on 0<Re s<1 (from Mathlib's f_modif minus
   two indicator Mellin transforms), then mellin_eq_fourier: Λ(1/2+it) = ½·𝓕Hk(t/4π).
   Λ(1/2+it) is even (FE) and real (Hk real).
3. Shift.lean — rectangle-Cauchy contour shift on a horizontal strip; gives
   Λ(1/2+it)e^{αt} = ½·𝓕(Hk(·-2αi))(t/4π) and ‖Λ(1/2+it)‖ ≤ C e^{-γ|t|}, γ < π/4.
4. Moments.lean — J m α = ∫Λ(1/2+it)t^m e^{αt}dt: Fourier inversion at 0 gives J 0 α = 2πHk(-2αi);
   differentiation under ∫ + induction gives J m α = 2π(-2i)^m Hk^{(m)}(-2αi), |α|<π/4.
5. Limit.lean — Hk^{(2n)} = Rk^{(2n)} - 16^{-n}Ck; at the cusp τ = -1 the shift θ(τ) = θ_{1/2}(τ+1)
   plus jacobiTheta₂_functional_equation gives ‖Rk u‖ ≤ (K/δ)e^{-1/(72δ)} on circles of radius δ/2
   about (-π/2+δ)i; Cauchy estimates ⇒ Rk^{(m)}(-2αi) → 0 as α → π/4⁻.
6. Hardy.lean — finitely many zeros ⇒ L(t)=Λ(1/2+it) has constant sign s for |t|>T (IVT + evenness);
   for large n, s·Re J(2n)α ≥ m₀(T+1)^{2n} - A·T^{2n} > 0 uniformly in α∈[0,π/4); but
   s·Re J(2n)α = 2π·s(-4)^n·(Re Rk^{(2n)}(-2αi) - 16^{-n}2cos(α/2)) < 0 for s(-1)^n>0, α near π/4.

Module names differ from the brief's plan (Xi/Fourier/Limit/Hardy): Xi ↔ Repr (Lam_conj, Lam_neg),
Fourier ↔ Repr+Shift+Moments. No Mathlib infrastructure was missing.
Not touched: Riemann.lean (orchestrator may add import Riemann.Deep.Hardy.Hardy), TREE.md, git.
Lint noise only: push_neg deprecation, one unused variable warning.
Significance: known theorem (Hardy 1914); formal infrastructure, not a discovery.
