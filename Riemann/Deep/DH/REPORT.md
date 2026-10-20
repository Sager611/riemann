# DH/FE — Davenport–Heilbronn functional equation (done, no sorry)

File: `Riemann/Deep/DH/FE.lean` (namespace `Fayed.DH`). `lake build Riemann.Deep.DH.FE` green; axioms of all main theorems = [propext, Classical.choice, Quot.sound].

1. `χ`: explicit mod-5 character `0,1,i,-i,-1` (`χ 2 = I`, odd). `χ_isPrimitive`, `χinv_isPrimitive` (conductor ∣ 5 prime, χ ≠ 1).
2. `κ = (√(10-2√5)-2)/(√5-1)`; proved `κ_eq : κ = 2 sin(2π/5) - (√5+1)/2`. `DH s = a·L(χ,s) + b·L(χ̄,s)`, `a=(1-iκ)/2`, `b=(1+iκ)/2`.
3. Root numbers from explicit Gauss sums (`τ(χ) = ζ+iζ²-iζ³-ζ⁴`, `ζ = (√5-1)/4 + i sin(2π/5)`): `a_mul_rootNumber : a·ε(χ) = b`, `b_mul_rootNumber_inv : b·ε(χ̄) = a` (ring identity via `linear_combination`, coefficients from a sympy Gröbner reduction).
   - `completedDH_one_sub`: `completedDH (1-s) = completedDH s` for ALL s, `completedDH s = 5^(s/2)(a Λ(χ,s) + b Λ(χ̄,s))` (Mathlib completed L, entire).
   - `DH_functional_equation`: `5^((1-s)/2) Γℝ(2-s) DH(1-s) = 5^(s/2) Γℝ(s+1) DH(s)` off Gamma poles.
   - `classicalDH_one_sub`: `Λ(1-s)=Λ(s)`, `Λ(s) = (5/π)^(s/2) Γ((s+1)/2) DH(s)`, off poles (Lean's Γ is junk 0 at poles, so hypotheses needed there only).
4. `coeff n = a χ(n) + b χ̄(n)` takes values `1, κ, -κ, -1, 0` (n ≡ 1,2,3,4,0 mod 5): `coeff_values`, `coeff_real`, `coeff_im`; `DH_eq_LSeries` for Re s > 1.
- Off-line zero near 0.808517 + 85.699348 i: cited def Prop `OffLineZeroClaim` only (Spira 1994), not proved.
- Not done: conj-symmetry `DH(s̄) = conj DH(s)` of the continuation (only coefficient realness).
