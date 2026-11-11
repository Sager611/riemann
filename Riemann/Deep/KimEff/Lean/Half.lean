import Riemann.Deep.KimEff.Lean.Xi
import Mathlib.NumberTheory.LSeries.Dirichlet

/-!
# Explicit Kim–Lee wedge, part 7: the classical-strip variant a = 1/2

With a = 1/2 the zero-location input (I2) is Mathlib's ζ(s) ≠ 0 for Re s > 1, so it is no
longer a hypothesis; (i) must then be the a = 1/2 zero-free region (certify.py 0.5,
logs/certify_half.log), and (iii) reads d ≤ 3.59 n^2/log^2 n.
-/

namespace KimEff

noncomputable section

open Real

/-- (I2) for a = 1/2, from Mathlib: ζ(s) ≠ 0 for Re s > 1. -/
theorem halfPlane_half : HalfPlane (1 / 2) := fun s hs =>
  riemannZeta_ne_zero_of_one_lt_re (by linarith)

/-- H3 for a = 1/2: 3.59 n^2/log^2 n ≤ D(R(n, u)) when u > 0, u + 1 ≤ (log n)/2, n ≥ n_*. -/
theorem le_Dfun_half {n u : ℝ} (hn : (nStar : ℝ) ≤ n) (hu : 0 < u)
    (hub : u + 1 ≤ log n / 2) : 3.59 * n ^ 2 / log n ^ 2 ≤ Dfun (1 / 2) (Rfun n u) := by
  have hn' := nStar_le hn
  have hn0 : 0 < n := lt_of_lt_of_le (by norm_num) hn'
  have hlog := log_ge_27 hn
  have hl : 0 < log n := by linarith
  have hl2 : 0 < log n ^ 2 := pow_pos hl 2
  have hnh : 0 < n - 1 / 2 := by linarith
  have hR : 0 < Rfun n u := by
    unfold Rfun
    exact div_pos (mul_pos (by norm_num) (pow_pos hnh 2)) (pow_pos (by linarith) 2)
  have hRlow : 3.6 * (n - 1 / 2) ^ 2 / log n ^ 2 ≤ Rfun n u := by
    unfold Rfun
    have h1 : (u + 1) ^ 2 ≤ (log n / 2) ^ 2 := pow_le_pow_left₀ (by linarith) hub 2
    calc 3.6 * (n - 1 / 2) ^ 2 / log n ^ 2 = 9 / 10 * (n - 1 / 2) ^ 2 / (log n / 2) ^ 2 := by
          field_simp
          ring
      _ ≤ 9 / 10 * (n - 1 / 2) ^ 2 / (u + 1) ^ 2 :=
          div_le_div_of_nonneg_left (by positivity) (pow_pos (by linarith) 2) h1
  have h4 : (0 : ℝ) < 4 * (1 / 2) ^ 2 := by norm_num
  have hmain : 3.59 * n ^ 2 / log n ^ 2
      ≤ (3.6 * (n - 1 / 2) ^ 2 / log n ^ 2 - (1 / 2) ^ 2) / (4 * (1 / 2) ^ 2) := by
    rw [div_le_div_iff₀ hl2 h4, sub_mul, div_mul_cancel₀ _ hl2.ne']
    have hls := log_sq_le (by linarith : (1:ℝ) ≤ n)
    have hnn : 998616377771 * n ≤ n * n := mul_le_mul_of_nonneg_right hn' hn0.le
    nlinarith
  calc 3.59 * n ^ 2 / log n ^ 2
      ≤ (3.6 * (n - 1 / 2) ^ 2 / log n ^ 2 - (1 / 2) ^ 2) / (4 * (1 / 2) ^ 2) := hmain
    _ ≤ (Rfun n u - (1 / 2) ^ 2) / (4 * (1 / 2) ^ 2) :=
        div_le_div_of_nonneg_right (by linarith) h4.le
    _ ≤ Dfun (1 / 2) (Rfun n u) := sub_div_le_Dfun (by norm_num) hR

/-- Main theorem for Riemann's ξ with the classical strip (a = 1/2): no zero-location
hypothesis; (ii) J^{d,n} hyperbolic for d ≤ D_n and (iii) for d ≤ 3.59 n^2/log^2 n, n ≥ n_*. -/
theorem xi_jensen_hyperbolic_half (hO : ObreschkoffKimLee)
    (hHad : ∀ n, HadamardApprox xiGamma ThetaDeriv n)
    (hGL : ∀ n, GaussLucasGenus0 ThetaDeriv n) (hpos : ∀ m, 0 < xiGamma m)
    (hi : ZeroFreeRegion (1 / 2) ThetaDeriv) :
    (∀ n : ℕ, nStar ≤ n → ∀ u : ℝ, 0 < u → Zfun u = n → ∀ d : ℕ,
        (d : ℝ) ≤ Dfun (1 / 2) (Rfun n u) → Hyperbolic (jensen xiGamma d n)) ∧
    (∀ n : ℕ, nStar ≤ n → ∀ d : ℕ, (d : ℝ) ≤ 3.59 * (n : ℝ) ^ 2 / Real.log n ^ 2 →
        Hyperbolic (jensen xiGamma d n)) := by
  have H : Inputs xiGamma ThetaDeriv (1 / 2) :=
    ⟨hO, hHad, hGL, zerosSq_xi (by norm_num) halfPlane_half, hpos, hi⟩
  refine ⟨fun n hn u hu hZ d hd => jensen_hyperbolic_ii (by norm_num) H hn hu hZ hd, ?_⟩
  intro n hn d hd
  have hn' : (nStar : ℝ) ≤ n := by exact_mod_cast hn
  obtain ⟨u, hu, hZ, hub⟩ := exists_u hn'
  exact jensen_hyperbolic_ii (by norm_num) H hn hu hZ (hd.trans (le_Dfun_half hn' hu hub))

end

end KimEff

