import Riemann.Deep.DPCBeyond.Lean.LemmaM
import Riemann.Deep.DPCBeyond.Lean.TheoremUB

/-!
# dpcbeyond Lemma M3's assembly: ∫_T^{2T}|R_T|^q ≪ T for q < κ/θ from named norm inequalities

PROOF.md M3: at p = qr (r > 1, pθ < κ(1 − 1/w)), Lemma Z at y = 0 and Hölder give
∫|R_T|^q ≤ K_Z (∫(1+W)^{qr})^{1/r}(∫(1+On)^{q(k−1)r'})^{1/r'}; M1 gives ∫(1+On)^{s} ≪ T and ‖On‖_p ≪ T^{1/p};
Minkowski with W ≤ On + Off gives ‖1+W‖_p ≤ T^{1/p} + ‖On‖_p + ‖Off‖_p; Off ≤ Σ_a e^{θ(a+1)}F_a + O(T^{−5}) gives
‖Off‖_p ≤ Σ_{a<M} e^{θ(a+1)}‖F_a‖_p + T^{1/p}; M2 (+ Minkowski over m and the kernel sum) gives
‖F_a‖_p ≤ K₁T^{1/p}e^{−κa(1−1/w)/p}. All five are NAMED (M3Inputs); proved: the bin sum (off_bound, via M3_series),
the exponent arithmetic (momentM_of) and the choice of r, w (exists_r_w), hence MomentM (momentM_of_inputs).
-/

namespace DPCBLean

open Finset DPCLean

variable {ex lg : ℝ → ℝ}

/-- M3's bin sum: ‖F_a‖_p ≤ K₁τe^{−κa(1−1/w)/p} for all a and pθ < κ(1 − 1/w) give
Σ_{a<M} e^{θ(a+1)}‖F_a‖_p ≤ K₁τ e^θ/(1 − e^{θ−κ(1−1/w)/p}). -/
lemma off_bound (h : ExpLog ex lg) {θ κ w p K₁ τ : ℝ} (hp : 0 < p) (hpθ : p * θ < κ * (1 - 1 / w))
    (hK₁ : 0 ≤ K₁) (hτ : 0 ≤ τ) (NF : ℕ → ℝ)
    (hF : ∀ a : ℕ, NF a ≤ K₁ * τ * ex (-(κ * a * (1 - 1 / w) / p))) (M : ℕ) :
    ∑ a ∈ range M, ex θ ^ (a + 1) * NF a ≤
      K₁ * τ * (ex θ / (1 - ex (θ - κ * (1 - 1 / w) / p))) := by
  have e : ∀ a : ℕ, ex θ ^ (a + 1) * ex (-(κ * a * (1 - 1 / w) / p)) =
      ex (θ * (a + 1) - κ * a * (1 - 1 / w) / p) := by
    intro a
    rw [← h.nat_mul, ← h.add]
    congr 1
    push_cast
    ring
  calc ∑ a ∈ range M, ex θ ^ (a + 1) * NF a ≤
        ∑ a ∈ range M, K₁ * τ * ex (θ * (a + 1) - κ * a * (1 - 1 / w) / p) := by
        apply sum_le_sum
        intro a _
        rw [← e a]
        have h1 := mul_le_mul_of_nonneg_left (hF a) (pow_nonneg (h.pos θ).le (a + 1))
        have e2 : ex θ ^ (a + 1) * (K₁ * τ * ex (-(κ * a * (1 - 1 / w) / p))) =
            K₁ * τ * (ex θ ^ (a + 1) * ex (-(κ * a * (1 - 1 / w) / p))) := by ring
        linarith only [h1, e2]
    _ = K₁ * τ * ∑ a ∈ range M, ex (θ * (a + 1) - κ * a * (1 - 1 / w) / p) := by rw [mul_sum]
    _ ≤ K₁ * τ * (ex θ / (1 - ex (θ - κ * (1 - 1 / w) / p))) :=
        mul_le_mul_of_nonneg_left (M3_series h hp hpθ M) (mul_nonneg hK₁ hτ)

/-- **M3's assembly** at p = qr: the five named inequalities give ∫_T^{2T}|R_T|^q ≤ K T. -/
theorem momentM_of (h : ExpLog ex lg) {θ κ w q r KZ K₀ K₁ T₀ : ℝ} (hq : 0 < q) (hr : 1 < r)
    (hpθ : q * r * θ < κ * (1 - 1 / w)) (hKZ : 0 ≤ KZ) (hK₀ : 0 ≤ K₀) (hK₁ : 0 ≤ K₁)
    {Rq IW IO NOn NOff : ℝ → ℝ} {NF : ℝ → ℕ → ℝ}
    (hZ : ∀ T : ℝ, T₀ ≤ T → Rq T ≤ KZ * rp ex lg (IW T) (1 / r) * rp ex lg (IO T) (1 - 1 / r))
    (hIO : ∀ T : ℝ, T₀ ≤ T → IO T ≤ K₀ * T)
    (hIW : ∀ T : ℝ, T₀ ≤ T →
      rp ex lg (IW T) (1 / (q * r)) ≤ rp ex lg T (1 / (q * r)) + NOn T + NOff T)
    (hOn : ∀ T : ℝ, T₀ ≤ T → NOn T ≤ K₀ * rp ex lg T (1 / (q * r)))
    (hOff : ∀ T : ℝ, T₀ ≤ T → ∃ M : ℕ,
      NOff T ≤ ∑ a ∈ range M, ex θ ^ (a + 1) * NF T a + rp ex lg T (1 / (q * r)))
    (hF : ∀ T : ℝ, T₀ ≤ T → ∀ a : ℕ,
      NF T a ≤ K₁ * rp ex lg T (1 / (q * r)) * ex (-(κ * a * (1 - 1 / w) / (q * r)))) :
    ∃ K T₁ : ℝ, ∀ T : ℝ, T₁ ≤ T → Rq T ≤ K * T := by
  have hp : 0 < q * r := mul_pos hq (by linarith only [hr])
  have hpθ' : q * r * θ < κ * (1 - 1 / w) := hpθ
  have hr0 : 0 < r := by linarith only [hr]
  have hir : 0 ≤ 1 / r := by positivity
  have hir1 : 1 / r ≤ 1 := by rw [div_le_one hr0]; exact hr.le
  have hir' : 0 ≤ 1 - 1 / r := by linarith only [hir1]
  obtain ⟨S, hS⟩ : ∃ S : ℝ, S = ex θ / (1 - ex (θ - κ * (1 - 1 / w) / (q * r))) := ⟨_, rfl⟩
  have hlt : θ < κ * (1 - 1 / w) / (q * r) := by
    rw [lt_div_iff₀ hp, mul_comm]
    exact hpθ'
  have hρ1 : ex (θ - κ * (1 - 1 / w) / (q * r)) < 1 := h.lt_one (by linarith only [hlt])
  have hS0 : 0 ≤ S := by
    rw [hS]
    exact div_nonneg (h.pos θ).le (by linarith only [hρ1])
  obtain ⟨K₃, hK₃⟩ : ∃ K₃ : ℝ, K₃ = 2 + K₀ + K₁ * S := ⟨_, rfl⟩
  have hK₃0 : 0 < K₃ := by
    rw [hK₃]
    have h1 := mul_nonneg hK₁ hS0
    linarith only [h1, hK₀]
  have hK₀1 : 0 < K₀ + 1 := by linarith only [hK₀]
  refine ⟨KZ * rp ex lg (rp ex lg K₃ (q * r)) (1 / r) * rp ex lg (K₀ + 1) (1 - 1 / r),
    max T₀ 1, fun T hT => ?_⟩
  have hT1 : T₀ ≤ T := le_trans (le_max_left _ _) hT
  have hT0 : 0 < T := lt_of_lt_of_le one_pos (le_trans (le_max_right _ _) hT)
  obtain ⟨τ, hτ⟩ : ∃ τ : ℝ, τ = rp ex lg T (1 / (q * r)) := ⟨_, rfl⟩
  have hτ0 : 0 < τ := by rw [hτ]; exact rp_pos h hT0 _
  obtain ⟨M, hM⟩ := hOff T hT1
  have hsum := off_bound h hp hpθ' hK₁ hτ0.le (NF T) (by have h1 := hF T hT1; rw [← hτ] at h1; exact h1) M
  rw [← hS] at hsum
  have hW : rp ex lg (IW T) (1 / (q * r)) ≤ K₃ * τ := by
    have h1 := hIW T hT1
    have h2 := hOn T hT1
    rw [← hτ] at h1 h2 hM
    have e : K₃ * τ = τ + K₀ * τ + (K₁ * τ * S + τ) := by rw [hK₃]; ring
    linarith only [h1, h2, hM, hsum, e]
  have hIW' : IW T ≤ rp ex lg K₃ (q * r) * T := by
    rcases le_or_gt (IW T) 0 with h0 | h0
    · exact le_trans h0 (mul_nonneg (rp_nonneg h _ _) hT0.le)
    · calc IW T = rp ex lg (rp ex lg (IW T) (1 / (q * r))) (q * r) := by
            rw [rp_rp h h0, one_div_mul_cancel hp.ne', rp_one h h0]
        _ ≤ rp ex lg (K₃ * τ) (q * r) := rp_le_rp h hW hp.le
        _ = rp ex lg K₃ (q * r) * T := by
            rw [rp_mul h hK₃0 hτ0, hτ, rp_rp h hT0, one_div_mul_cancel hp.ne', rp_one h hT0]
  have hIO' : IO T ≤ (K₀ + 1) * T := by
    have h1 := hIO T hT1
    have h2 : K₀ * T ≤ (K₀ + 1) * T := mul_le_mul_of_nonneg_right (by linarith only) hT0.le
    linarith only [h1, h2]
  have r1 : rp ex lg (IW T) (1 / r) ≤ rp ex lg (rp ex lg K₃ (q * r)) (1 / r) * rp ex lg T (1 / r) := by
    rw [← rp_mul h (rp_pos h hK₃0 _) hT0]
    exact rp_le_rp h hIW' hir
  have r2 : rp ex lg (IO T) (1 - 1 / r) ≤ rp ex lg (K₀ + 1) (1 - 1 / r) * rp ex lg T (1 - 1 / r) := by
    rw [← rp_mul h hK₀1 hT0]
    exact rp_le_rp h hIO' hir'
  have eT : rp ex lg T (1 / r) * rp ex lg T (1 - 1 / r) = T := by
    rw [← rp_add h hT0, show 1 / r + (1 - 1 / r) = 1 by ring, rp_one h hT0]
  calc Rq T ≤ KZ * rp ex lg (IW T) (1 / r) * rp ex lg (IO T) (1 - 1 / r) := hZ T hT1
    _ = KZ * (rp ex lg (IW T) (1 / r) * rp ex lg (IO T) (1 - 1 / r)) := by ring
    _ ≤ KZ * ((rp ex lg (rp ex lg K₃ (q * r)) (1 / r) * rp ex lg T (1 / r)) *
        (rp ex lg (K₀ + 1) (1 - 1 / r) * rp ex lg T (1 - 1 / r))) :=
        mul_le_mul_of_nonneg_left (mul_le_mul r1 r2 (rp_nonneg h _ _)
          (mul_nonneg (rp_nonneg h _ _) (rp_nonneg h _ _))) hKZ
    _ = KZ * rp ex lg (rp ex lg K₃ (q * r)) (1 / r) * rp ex lg (K₀ + 1) (1 - 1 / r) *
        (rp ex lg T (1 / r) * rp ex lg T (1 - 1 / r)) := by ring
    _ = KZ * rp ex lg (rp ex lg K₃ (q * r)) (1 / r) * rp ex lg (K₀ + 1) (1 - 1 / r) * T := by
        rw [eT]

/-- Order of choices for M3: for 0 < q < κ/θ there are r > 1 and w > 1 with qrθ < κ(1 − 1/w). -/
theorem exists_r_w {θ κ q : ℝ} (hθ : 0 < θ) (hq : 0 < q) (hqκ : q < κ / θ) :
    ∃ r w : ℝ, 1 < r ∧ 1 < w ∧ q * r * θ < κ * (1 - 1 / w) := by
  have hqθ : q * θ < κ := by rwa [lt_div_iff₀ hθ] at hqκ
  have hqθ0 : 0 < q * θ := mul_pos hq hθ
  have hk : 1 < κ / (q * θ) := by rw [lt_div_iff₀ hqθ0]; linarith only [hqθ]
  obtain ⟨r, hr⟩ : ∃ r : ℝ, r = (1 + κ / (q * θ)) / 2 := ⟨_, rfl⟩
  have hr1 : 1 < r := by rw [hr]; linarith only [hk]
  have hr2 : r * (q * θ) < κ := by
    have h1 : r < κ / (q * θ) := by rw [hr]; linarith only [hk]
    rwa [lt_div_iff₀ hqθ0] at h1
  have hκ : 0 < κ := by linarith only [hqθ, hqθ0]
  obtain ⟨w, hw1, hw2, -⟩ := M3_exists_w (θ := θ) (p := q * r) hκ (by positivity)
    (by linarith only [hr2])
  exact ⟨r, w, hr1, hw1, by linarith only [hw2]⟩

/-- The named inputs of M3 (Lemma Z + Hölder, M1, Minkowski with W ≤ On + Off, the bins of Off, M2) for every
q ∈ (2, κ/θ) and every r, w > 1 with qrθ < κ(1 − 1/w). -/
def M3Inputs (ex lg : ℝ → ℝ) (κ θ : ℝ) (D : DataB) : Prop :=
  ∀ q r w : ℝ, 2 < q → q < κ / θ → 1 < r → 1 < w → q * r * θ < κ * (1 - 1 / w) →
    ∃ (KZ K₀ K₁ T₀ : ℝ) (IW IO NOn NOff : ℝ → ℝ) (NF : ℝ → ℕ → ℝ), 0 ≤ KZ ∧ 0 ≤ K₀ ∧ 0 ≤ K₁ ∧
      (∀ T : ℝ, T₀ ≤ T → D.Rq q T ≤ KZ * rp ex lg (IW T) (1 / r) * rp ex lg (IO T) (1 - 1 / r)) ∧
      (∀ T : ℝ, T₀ ≤ T → IO T ≤ K₀ * T) ∧
      (∀ T : ℝ, T₀ ≤ T →
        rp ex lg (IW T) (1 / (q * r)) ≤ rp ex lg T (1 / (q * r)) + NOn T + NOff T) ∧
      (∀ T : ℝ, T₀ ≤ T → NOn T ≤ K₀ * rp ex lg T (1 / (q * r))) ∧
      (∀ T : ℝ, T₀ ≤ T → ∃ M : ℕ,
        NOff T ≤ ∑ a ∈ range M, ex θ ^ (a + 1) * NF T a + rp ex lg T (1 / (q * r))) ∧
      (∀ T : ℝ, T₀ ≤ T → ∀ a : ℕ,
        NF T a ≤ K₁ * rp ex lg T (1 / (q * r)) * ex (-(κ * a * (1 - 1 / w) / (q * r))))

/-- **Lemma M3 from its named inputs**: MomentM (∫_T^{2T}|R_T|^q ≪ T for every q ∈ (2, κ/θ)). -/
theorem momentM_of_inputs (h : ExpLog ex lg) {κ θ : ℝ} (hθ : 0 < θ) {D : DataB}
    (hin : M3Inputs ex lg κ θ D) : MomentM κ θ D := by
  intro q hq hqκ
  obtain ⟨r, w, hr, hw, hrw⟩ := exists_r_w hθ (by linarith only [hq]) hqκ
  obtain ⟨KZ, K₀, K₁, T₀, IW, IO, NOn, NOff, NF, hKZ, hK₀, hK₁, hZ, hIO, hIW, hOn, hOff, hF⟩ :=
    hin q r w hq hqκ hr hw hrw
  exact momentM_of h (by linarith only [hq]) hr hrw hKZ hK₀ hK₁ hZ hIO hIW hOn hOff hF

end DPCBLean

