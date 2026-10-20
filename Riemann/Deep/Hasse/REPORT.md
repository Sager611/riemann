NO DISCOVERY: Hasse's theorem is already formalized in Lean 4 (Meta's ATLAS library, github.com/facebookresearch/atlas-lean, file MathlibExt/AlgebraicGeometry/EllipticCurve/HasseBound.lean, added 2026-10-05, Manin's twist, all characteristics, no sorry). Ours is an independent second proof.

# Hasse's theorem (RH for elliptic curves over finite fields), sorry-free

Build: `lake build Riemann.Deep.Hasse.Hasse` (green, no warnings). Axioms of every main theorem:
propext, Classical.choice, Quot.sound only. No native_decide. 1187 lines in 5 files.
Each file compiles in about 6 s with peak memory about 2.6 GB. Nothing runs at idle.

## What is proved (namespace Hasse)
- Rung 3, general case, odd characteristic (F finite, ringChar F ≠ 2):
  - `sum_quadraticChar_sq_le`: (∑ₓ χ(x³+c₂x²+c₁x+c₀))² ≤ 4q, for every monic cubic, singular ones included.
  - `hasse`: |#{y²=x³+c₂x²+c₁x+c₀} + 1 − (q+1)| ≤ 2√q.
  - `hasse_weierstrass`: ((q+1) − (#affine points + 1))² ≤ 4q for any Mathlib WeierstrassCurve.
  - `hasse_elliptic`: ((q+1) − Nat.card W.toAffine.Point)² ≤ 4q for W elliptic (the same statement as ATLAS).
- Rung 2, the families: `hasse_x3_add` (y²=x³+D), `hasse_x3_sub` (y²=x³−Dx), plus exact counts
  `card_x3_add_of_mod_three` (q ≡ 2 mod 3 ⇒ q affine points) and
  `card_x3_sub_of_mod_four` (q ≡ 3 mod 4 ⇒ q affine points).
- Rung 1, y²=x³−x over ZMod p: `rung1_bound` (p odd prime) and `rung1_exact` (p ≡ 3 mod 4 ⇒ a_p = 0).

## Method (Manin 1956, written with upper bounds only)
Twist f(t)y² = f(x) over K = F(t), written as a Mathlib Weierstrass curve with X = f(t)x and Y = f(t)²y.
It has the points P₀ = (t,1) and P₁ = (t^q, f(t)^((q−1)/2)). ht(P) is the degree of x(P) in lowest terms.
1. Height.lean: `htK_add_htK_le`. If A(x₃+x₄) = B and A·x₃x₄ = C, then ht x₃ + ht x₄ ≤ max deg(A,B,C).
   Proof: explicit Bezout identity plus a leading-coefficient argument.
2. Twist.lean: `ht_pgram`. ht(P+P₀) + ht(P−P₀) ≤ 2 ht P + 2 for every P, using Mathlib's group law.
   Inputs: the classical formulas for x(P+Q)+x(P−Q) and x(P+Q)·x(P−Q), and a doubling bound.
3. Count.lean: ht P₁ = q, and ht(P₁ ± P₀) ≤ q + 1 ∓ S. The common factor ∏(X−a)^(1±χ(f a))
   cancels from numerator and denominator (Frobenius f(t)^q = f(t^q), Euler's criterion).
4. Main.lean: a sequence with second differences ≤ 2 stays below n² − Sn + q (`le_parabola`).
   Heights are ≥ 0, and q is odd, so S² ≤ 4q (`sq_le_of_parabola_nonneg`).
Only upper bounds are used, so the proof never needs the exact parallelogram law, separability, or a coprimality
argument. That is why singular cubics come for free. A small simplification, not a new result.

## Not done
- Characteristic 2 (Artin–Schreier twist; poles at h(t)=0 complicate step 3). ATLAS has it.
- Exact CM count a_p = ±2a with p = a²+b², p ≡ 1 mod 4 (needs quartic Jacobi sums). Hasse–Weil zeta / F_{q^n}.

## Files
Riemann/Deep/Hasse/{Height,Twist,Count,Main,Hasse}.lean; scratch copies in .scratch/hasse/.
Not imported from Riemann.lean; TREE.md, NOGO.md, DISCOVERIES.md untouched; no commits (orchestrator's call).
