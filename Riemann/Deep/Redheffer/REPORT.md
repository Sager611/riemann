# Redheffer: det R_n = M(n) — DONE

- File: `Riemann/Deep/Redheffer/Det.lean` (namespace `Fayed.Redheffer`), no sorry; axioms [propext, Classical.choice, Quot.sound].
- `det_redheffer (n) (hn : 0 < n) : (redheffer n).det = ∑ k ∈ Icc 1 n, μ k`.
- Statement as briefed is FALSE at n = 0 (empty det = 1, M(0) = 0) → added `hn : 0 < n`; `det_redheffer_zero : (redheffer 0).det = 1`.
- Route (Cramer, not the det lemma): R = updateCol D 0 𝟙; det D = 1 (upper unitriangular); E·D = 1 with E = Möbius matrix; cramer D 𝟙 = E *ᵥ 𝟙; entry 0 = Σ_{j<n} μ(j+1).
- Key Mathlib: `Matrix.cramer_apply`, `Matrix.mulVec_cramer`, `Matrix.mulVec_mulVec`, `Matrix.det_of_isUpperTriangular`,
  `ArithmeticFunction.coe_mul_zeta_apply` + `moebius_mul_coe_zeta` (Σ_{d|c} μ d = [c=1]), `Fin.sum_univ_eq_sum_range`, `Finset.sum_Icc_succ_top`.
- Helper lemmas: `det_divMat`, `sum_divisors_moebius`, `sum_divisors_dvd_moebius` (Σ_{m|b, a|m} μ(m/a) = [a=b]), `moebMat_mul_divMat`, `redheffer_eq_updateCol`.
- Verify: `lake build Riemann.Deep.Redheffer.Det` → green, no warnings expected besides #print axioms info.
- Gap: none. Not imported by `Riemann.lean` (brief forbade editing other files).
