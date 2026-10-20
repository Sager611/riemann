# Yoshida deep dive: Weil positivity on small windows (Lean)

**Verdict (round 4): the whole Yoshida range is formalized.** `Fayed.Yoshida.weilPos_lt_log2 : ∀ L < log 2,
WeilPos L` (`Log2.lean`), sorry-free, no `native_decide`, axioms = propext, Classical.choice, Quot.sound.
`weilPos_le_60` and `weilPos_le_65` are corollaries. This is Yoshida's 1992 theorem, proved by our own
numerical route, not by his argument.

## Round 4: what changed
1. **Affine cells** (`Affine.lean`, `Chord.lean`). The deficit is `t = A + B|r|` per cell. Two chord bounds:
   `2πμ` is concave in `u = (r/2)²` (every digamma term `−c/(c²+u)` is), and `y ↦ cos(yx) − x²y²/2` is
   concave, so `J(r) ≥ chord − (16/3)a³cosh a·θ(1−θ)(hi−lo)²`. The loss per cell is second order:
   `δ = (mH − mL)(hi − lo)/(4(hi + lo)) + E_j(hi − lo)²`. Cells 1/64 wide now beat round 3's 1/256 cells.
2. **Wider support** (`weilQ_re_eq_W'`). For `f` supported in `[−a, a]` with `2a < log 2`, the prime sum is
   still zero, but the pole kernel `2cosh(x/2)` may be cut at any `2a₁ ≥ 2a`, giving the weight `W(a₁, ·)`;
   Taylor and Bessel steps also hold with `a₁`. So **one** certificate at a rational `a₁ ≥ log 2 / 2` covers
   every window below `log 2`. The certificate is at `a₁ = 1733/5000 = 0.3466` (log 2 / 2 = 0.346574).
3. Sharper constants: `log π < 1.14474` (`log_pi_lt'`); `2πμ` and `J` are bounded once per grid point
   (shared by neighbouring cells, `chkA`).

| certificate | cells | K | s | margin (2π units, α ≈ 6.9) |
|---|---|---|---|---|
| `WLog2` (a₁ = 0.3466) | 1179 (1172 affine) up to R = 23.6 | 24 | 1.1 | 0.0038 |
| test only: a = 0.3 (L = 0.6) | 742 | 16 | 0.8 | 0.0188 (deleted, superseded) |

The float model (`.scratch/yoshida4/aff.py`) puts the ideal margin at log 2 near 0.0054; rounding and cell
loss eat 30%. K = 32, s = 1.3 gives 0.0050 if a later change needs headroom.

## Cost
`lake build Riemann.Deep.Yoshida.Log2`: `WLog2` alone ≈ 19 min wall standalone (443 s user), 52 min wall when
the machine was busy; peak footprint ≈ 300 MB, resident ≈ 5–7 GB (mapped Mathlib). Not imported by
`Riemann.lean`. Generator: `cd .scratch/yoshida4 && python3.11 gen4.py 1733/5000 11/10 24 6 WLog2 150` (2 s).

## Yoshida 1992 / Connes–Consani 2021 (goal 2)
Not formalized and no longer needed: the full range `L < log 2` now follows from one rational certificate.
Yoshida's proof works with the same prime-free quadratic form but argues by operator bounds (the form is a
compact perturbation of a positive operator, checked numerically on a finite piece); Connes–Consani (arXiv
2006.13771) prove positivity on `[2^{-1/2}, 2^{1/2}]` via prolate spheroidal functions. Both rest on the same
fact our certificate checks: the smallest eigenvalue of the prime-free form stays positive (≈ 0.0013 in W units)
up to `log 2`. Beyond `log 2` the prime 2 enters and this method, as is, stops.

## Beyond log 2

DISCOVERY: first machine-checked proof of Weil positivity beyond `log 2`: `Fayed.Yoshida.weilPos_le_75 : ∀ L ≤ 3/4,
WeilPos L` (`Beyond.lean`; the prime 2 enters); the inequality itself was known only numerically (Zhu 1.6, Mori
`log 15`, our Arb 7/4). Sorry-free, no `native_decide`, axioms = propext, Classical.choice, Quot.sound.

**How.** `BeyondCore.lean`: for every window `2a < log 3` only `n = 2` enters, and
`Re Q(f) = (1/2π)∫|F|² W3` with `W3 = 2πμ + J(b,·) − √2 log 2·cos(r log 2) + Σ d_k cos(r x_k)`, any `b ≥ a` and any
free terms with `|x_k| > 2a` (they integrate to zero against `|F|²`: the autocorrelation vanishes there).
`weilPos_of_affine3` is the affine moment certificate for `W3`; it covers every `L ≤ 2a₁` at once.
1. **Free cosines (the key step).** With the bare prime weight the level-set moment bound needs `s ≥ 1.2` and
   deficit bumps up to `r ≈ 64` (Taylor order ≈ 80, ~10⁴ cells). Twelve free terms `x = jℓ/8`, `j = 9..20`,
   chosen by L-BFGS (`.scratch/yoshida5/opt2.py`), cut the bumps to `r ≈ 29` and keep 72% of the true margin
   (by Krein-type duality such terms make the bound exact in the limit).
2. **`ℤ` balls** (`BeyondBall.lean`, `BeyondReflect.lean`): the kernel does ~1 ms per `ℚ` operation but ~0.05 ms
   per `ℤ` one, so `2πμ`, `cos`, `sin` and the harmonics `cos(jθ)` (angle addition) are checked in fixed point
   `2⁻⁶²` with proved rounding: 16 ms and 115 ms per grid point instead of ~1 s.
3. **Certificate** `BeyondW75.lean` (`gen5.py`): `a = 3/8`, `s = 6/5`, 3440 cells (3204 affine) to `r = 196.7`,
   `K = 40`, margin 0.0016 in 2π units (true value 0.0031; float model 0.0017).

| window | true λ_min (W units) | ideal margin with free terms (2π units) | status |
|---|---|---|---|
| 0.70 | ≈ 0.0012 | ≈ 0.005 | covered by 0.75 |
| 0.75 | 4.9e-4 | 0.0022 | **proved** |
| 0.80 | 1.8e-4 | ≈ 0.0009 (without free terms) | next rung, tight |
| 1.0 | 9.3e-7 | — | needs μ to ~1e-7 and a new precision regime |

**Cost.** `lake build Riemann.Deep.Yoshida.Beyond`: 3.7 h wall on a loaded machine (7600 s user), peak RSS 6.9 GB.
Not imported by `Riemann.lean`. Most kernel time is the moment check `momCheckA`: `(Hi : ℤ) ^ p` unfolds by
recursion for all 3440 cells up to `p = 82`. Next: `Nat.pow` (GMP in the kernel) and skipping zero cells there;
then window 0.8 with more free terms and finer cells near `r = 0`.
Prior art: Mathlib, openai/math, facebookresearch/atlas-lean have no Weil positivity windows (checked 2026-10-20).

---

# Round 3 (kept for reference)

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
