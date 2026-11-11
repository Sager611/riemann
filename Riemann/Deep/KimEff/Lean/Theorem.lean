import Riemann.Deep.KimEff.Lean.Step
import Riemann.Deep.KimEff.Lean.Arith

/-!
# Explicit Kim–Lee wedge, part 4: the main theorem for an abstract sequence

For any real sequence γ, any family F n (playing Θ^{(n)}) and a > 0: the named inputs imply
(ii) J^{d,n} hyperbolic for n ≥ n_* and d ≤ D_n, and, for a = 0.37338196661320752,
(iii) J^{d,n} hyperbolic for d ≤ 6.45 n^2/log^2 n and for n ≥ max(n_*, 1.1853e-10 d).
-/

namespace KimEff

noncomputable section

open Real

/-- (i), DERIVATION §4–§9 with the Arb certificate certify.py: for n ≥ n_* and u = u(n)
(u > 0, Z(u) = n), Θ^{(n)} has no zero w with (Im w)^2 ≤ 4a^2(a^2 - Re w), Re w ≥ -R_n. -/
def ZeroFreeRegion (a : ℝ) (F : ℕ → ℂ → ℂ) : Prop :=
  ∀ n : ℕ, nStar ≤ n → ∀ u : ℝ, 0 < u → Zfun u = n → ∀ w ∈ wRegion a (Rfun n u), F n w ≠ 0

/-- The named inputs (cited theorems and the unformalized analytic steps; sources in LEAN.md). -/
structure Inputs (γ : ℕ → ℝ) (F : ℕ → ℂ → ℂ) (a : ℝ) : Prop where
  obreschkoff : ObreschkoffKimLee
  hadamard : ∀ n, HadamardApprox γ F n
  gaussLucas : ∀ n, GaussLucasGenus0 F n
  zerosSq : ZerosSq a F
  gammaPos : ∀ m, 0 < γ m
  zeroFree : ZeroFreeRegion a F

/-- (ii): J^{d,n} is hyperbolic for n ≥ n_*, u = u(n) and d ≤ D_n = D(R_n). -/
theorem jensen_hyperbolic_ii {γ : ℕ → ℝ} {F : ℕ → ℂ → ℂ} {a : ℝ} (ha : 0 < a)
    (H : Inputs γ F a) {n : ℕ} (hn : nStar ≤ n) {u : ℝ} (hu : 0 < u) (hZ : Zfun u = n)
    {d : ℕ} (hd : (d : ℝ) ≤ Dfun a (Rfun n u)) : Hyperbolic (jensen γ d n) := by
  have hn1 : (1 : ℝ) ≤ n := by
    have h : (1 : ℕ) ≤ n := le_trans (by norm_num [nStar]) hn
    exact_mod_cast h
  have hR : 0 < Rfun n u := by
    unfold Rfun
    exact div_pos (mul_pos (by norm_num) (pow_pos (by linarith) 2)) (pow_pos (by linarith) 2)
  exact hyperbolic_of_zeroFree_W ha hR H.obreschkoff (H.hadamard n) (H.gaussLucas n) H.zerosSq
    H.gammaPos (H.zeroFree n hn u hu hZ) hd

/-- (iii): J^{d,n} is hyperbolic for n ≥ n_* and d ≤ 6.45 n^2/log^2 n (a = 0.37338196661320752). -/
theorem jensen_hyperbolic_iii {γ : ℕ → ℝ} {F : ℕ → ℂ → ℂ} (H : Inputs γ F a0) {n : ℕ}
    (hn : nStar ≤ n) {d : ℕ} (hd : (d : ℝ) ≤ 6.45 * (n : ℝ) ^ 2 / log n ^ 2) :
    Hyperbolic (jensen γ d n) := by
  have hn' : (nStar : ℝ) ≤ n := by exact_mod_cast hn
  obtain ⟨u, hu, hZ, hub⟩ := exists_u hn'
  exact jensen_hyperbolic_ii a0_pos H hn hu hZ (hd.trans (le_Dfun_of_le_log hn' hu hub))

/-- (iii), power form (c, C) = (1, 1.1853e-10): n ≥ max(n_*, C d) suffices. -/
theorem jensen_hyperbolic_power_one {γ : ℕ → ℝ} {F : ℕ → ℂ → ℂ} (H : Inputs γ F a0) {n d : ℕ}
    (hn : nStar ≤ n) (hd : 1.1853e-10 * (d : ℝ) ≤ n) : Hyperbolic (jensen γ d n) :=
  jensen_hyperbolic_iii H hn (le_of_power_one (by exact_mod_cast hn) hd)

end

end KimEff

