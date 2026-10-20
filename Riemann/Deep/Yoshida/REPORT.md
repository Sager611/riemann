# Yoshida deep dive: Weil positivity on small windows (Lean)

**Verdict (round 3): window 0.55, up from 1/3; not yet log 2 ≈ 0.693.** Sorry-free, no `native_decide`,
axioms = propext, Classical.choice, Quot.sound:
- `Fayed.Yoshida.weilPos_le_55 : ∀ L ≤ 11/20, WeilPos L` (f supported in [−0.275, 0.275]).
- `weilPos_le_half` (L ≤ 1/2) and `weilPos_le_third` (now a corollary) in `Windows.lean`.
- Build: `lake build Riemann.Deep.Yoshida.Windows` ≈ 8 min wall (Half ≈ 2.5 min, W55 ≈ 3.5 min of kernel
  arithmetic), resident ≈ 6.9 GB (mostly mapped Mathlib). Not imported by Riemann.lean. `Third.lean` still
  builds on its own but is no longer needed.

## What changed: the cap is gone
Round 2 bounded `∫ |F|² t` (t = deficit of the weight `W = 2πμ + J` below a level s) by `|F|² ≤ 2aN`.
That is a trace bound and dies at L ≈ 0.355 whatever the cells. Round 3 bounds it by the Taylor
moments of F (a Galerkin step with a rigorous tail):
1. `F(r) = Σ_{m<K} r^m B_m + E(r)`, `B_m = ∫ f(u)(iu)^m/m!`, `|E| ≤ 2(|r|a)^K/K! ∫|f|` (Mathlib `exp_bound'`).
2. `∫ t |P|² = xᵀHx + yᵀHy` exactly (x, y = Re, Im of B; H = Hankel matrix of the moments of the step t).
3. Bessel: `N ≥ 2(d·x + e·y) − dᵀĜd − eᵀĜe` for all real d, e (from `0 ≤ ∫|f − ḡ|²`, g a polynomial).
4. With d = Λx: if `M = α(Λ + Λᵀ − ΛᵀĜΛ) − H ⪰ 0` then `∫ t|P|² ≤ αN`; the remainder adds β = 8a^{2K+1}·Mom_{2K}/K!².
5. `2π Re Q ≥ N (2πs − (1+ε)α − (1+1/ε)β) ≥ 0`.
`M ⪰ 0` is certified by `UᵀMU` diagonally dominant with `UW = I` (rational U, W).

Numerically this is close to sharp: the true lowest eigenvalue (W units) is 0.033 at L = 1/2, 0.0075 at
0.6, 0.0013 at log 2; the moment bound with exact W gives 0.025 / 0.005 / 0.0002 (`galerkin.py`).

## Reflective checking (new)
All numerics run as rational arithmetic inside the kernel (`decide +kernel`), with soundness lemmas:
digamma partial sums rounded down to 2⁻⁴⁰ plus a telescoping tail bound (`MuTail.two_pi_mu_ge_tel`),
γ < 0.57725 from a rounded H_{2^14} (`euler_lt'`), sin/cos by reduction mod 2π (π to 20 digits) and
alternating Taylor bounds, cells on a 1/256 grid with integer moment sums. About 0.1 s per cell.

## Files
| file | content |
|---|---|
| `Moment.lean` | Taylor split, polynomial cell integrals, Bessel inequality, `weilPos_of_moment` |
| `MuTail.lean` | telescoping digamma tail, sharper γ |
| `Reflect.lean` | rational checkers + soundness, `weilPos_of_qcert` |
| `Half.lean`, `W55.lean` | GENERATED certificates: a = 1/4 (950 cells, K = 14, margin 0.090 in 2π·W units); a = 11/40 (1696 cells, K = 16, s = 0.9, margin 0.031) |
| `Windows.lean` | `weilPos_half`, `weilPos_55`, `weilPos_le_55` |
| `.scratch/yoshida3/gen3.py` | generator (`python3.11 gen3.py 11/40 9/10 1/128 0 40 16 W55 150`); `qcheck.py` mirrors the Lean checker exactly |

## Why it stops at 0.55, and the next rung (numerics in `.scratch/yoshida3/`)
At L = 0.6 the certified margin is −0.003 (2π units) with 1/256 cells; the ideal moment bound is +0.020.
`diag.py` splits the loss: per-cell W lower bound ≈ 0.014 (J via a global Lipschitz constant, μ only at
the left end of a cell), μ accuracy ≈ 0.005 (M = 20 terms; the h₂ tail term is dropped when r < √2).
Planned fixes, in order of value:
1. **Affine cells** (t linear on each cell): μ(r) − μ(lo) ≥ (r − lo)·σ from the term-wise differences
   (no derivatives); J(r) ≥ J(p) + J′(p)(r − p) − C₂(r − p)²/2 from the cos inequality inside
   J = ∫ 2cosh(x/2) cos(rx). The loss becomes second order, so cells can be about 1/64 wide.
2. Sharper μ: tail factor c = 1 − 1/(2(M+5/4)²) for every r, M = 40, log π < 1.14473 (error 1e−4).
3. Larger K plus finer cells near the top eigenfunction. Expected reach: 0.6 is clear, 0.65 is tight
   (ideal 0.008); log 2 − ε needs μ to about 1e−5 and is not in reach with this method as it stands.
