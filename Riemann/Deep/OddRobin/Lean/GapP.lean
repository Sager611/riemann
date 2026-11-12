import Riemann.Deep.OddRobin.Lean.NumB
import Riemann.Deep.OddRobin.Lean.NumC

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

/-! ## (1a), (1e), (1f): per-prime analysis -/

lemma p3 (hL : LogLike L) (a : ℕ) :
    0 ≤ G L 3 4 - G L 3 a ∧ (G L 3 4 - G L 3 a < 13 / 10000 → a = 4 ∨ a = 5) := by
  have hp : Nat.Prime 3 := by norm_num
  have i3 := inc_3_3_ge hL
  have i4 := inc_3_4_le hL
  have i5 := inc_3_5_le hL
  have g4 : G L 3 4 = G L 3 3 + inc L 3 3 := G_succ hL hp 3
  have g5 : G L 3 5 = G L 3 4 + inc L 3 4 := G_succ hL hp 4
  have g6 : G L 3 6 = G L 3 5 + inc L 3 5 := G_succ hL hp 5
  rcases le_or_gt 4 a with h | h
  · have up : G L 3 a ≤ G L 3 4 := G_le_up hL hp (by linarith) h
    refine ⟨by linarith, fun hc => ?_⟩
    by_contra hne
    have h6 : 6 ≤ a := by omega
    have i6 : inc L 3 6 ≤ inc L 3 4 := inc_anti hL hp (by norm_num)
    have up6 : G L 3 a ≤ G L 3 6 := G_le_up hL hp (by linarith) h6
    linarith
  · have i2 : inc L 3 3 ≤ inc L 3 2 := inc_anti hL hp (by norm_num)
    have dn : G L 3 a ≤ G L 3 3 := G_le_down hL hp (b := 2) (by linarith) (by omega)
    exact ⟨by linarith, fun hc => by exfalso; linarith⟩

lemma low3 (hL : LogLike L) : 11 / 10000 < G L 3 4 - G L 3 5 := by
  have g5 : G L 3 5 = G L 3 4 + inc L 3 4 := G_succ hL (by norm_num) 4
  have := inc_3_4_le hL
  linarith

lemma p5 (hL : LogLike L) (a : ℕ) :
    0 ≤ G L 5 3 - G L 5 a ∧ (G L 5 3 - G L 5 a < 13 / 10000 → a = 2 ∨ a = 3) := by
  have hp : Nat.Prime 5 := by norm_num
  have i1 := inc_5_1_ge hL
  have i2 := inc_5_2_ge hL
  have i3 := inc_5_3_le hL
  have g2 : G L 5 2 = G L 5 1 + inc L 5 1 := G_succ hL hp 1
  have g3 : G L 5 3 = G L 5 2 + inc L 5 2 := G_succ hL hp 2
  have g4 : G L 5 4 = G L 5 3 + inc L 5 3 := G_succ hL hp 3
  rcases le_or_gt 3 a with h | h
  · have up : G L 5 a ≤ G L 5 3 := G_le_up hL hp (by linarith) h
    refine ⟨by linarith, fun hc => ?_⟩
    by_contra hne
    have h4 : 4 ≤ a := by omega
    have i4 : inc L 5 4 ≤ inc L 5 3 := inc_anti hL hp (by norm_num)
    have up4 : G L 5 a ≤ G L 5 4 := G_le_up hL hp (by linarith) h4
    linarith
  · have dn : G L 5 a ≤ G L 5 2 := G_le_down hL hp (b := 1) (by linarith) (by omega)
    refine ⟨by linarith, fun hc => ?_⟩
    by_contra hne
    have h1 : a ≤ 1 := by omega
    have i0 : inc L 5 1 ≤ inc L 5 0 := inc_anti hL hp (by norm_num)
    have dn1 : G L 5 a ≤ G L 5 1 := G_le_down hL hp (b := 0) (by linarith) (by omega)
    linarith

lemma low5 (hL : LogLike L) : 75 / 100000 < G L 5 3 - G L 5 2 := by
  have g3 : G L 5 3 = G L 5 2 + inc L 5 2 := G_succ hL (by norm_num) 2
  have := inc_5_2_ge hL
  linarith

lemma p7 (hL : LogLike L) (a : ℕ) :
    0 ≤ G L 7 2 - G L 7 a ∧ (G L 7 2 - G L 7 a < 13 / 10000 → a = 2) := by
  have hp : Nat.Prime 7 := by norm_num
  have i1 := inc_7_1_ge hL
  have i2 := inc_7_2_le hL
  have g2 : G L 7 2 = G L 7 1 + inc L 7 1 := G_succ hL hp 1
  have g3 : G L 7 3 = G L 7 2 + inc L 7 2 := G_succ hL hp 2
  rcases le_or_gt 2 a with h | h
  · have up : G L 7 a ≤ G L 7 2 := G_le_up hL hp (by linarith) h
    refine ⟨by linarith, fun hc => ?_⟩
    by_contra hne
    have h3 : 3 ≤ a := by omega
    have i3 : inc L 7 3 ≤ inc L 7 2 := inc_anti hL hp (by norm_num)
    have up3 : G L 7 a ≤ G L 7 3 := G_le_up hL hp (by linarith) h3
    linarith
  · have i0 : inc L 7 1 ≤ inc L 7 0 := inc_anti hL hp (by norm_num)
    have dn : G L 7 a ≤ G L 7 1 := G_le_down hL hp (b := 0) (by linarith) (by omega)
    exact ⟨by linarith, fun hc => by exfalso; linarith⟩

lemma p11 (hL : LogLike L) (a : ℕ) :
    0 ≤ G L 11 1 - G L 11 a ∧ (G L 11 1 - G L 11 a < 13 / 10000 → a = 1 ∨ a = 2) := by
  have hp : Nat.Prime 11 := by norm_num
  have i0 := inc_11_0_ge hL
  have i1 := inc_11_1_le hL
  have i2 := inc_11_2_le hL
  have g1 : G L 11 1 = G L 11 0 + inc L 11 0 := G_succ hL hp 0
  have g2 : G L 11 2 = G L 11 1 + inc L 11 1 := G_succ hL hp 1
  have g3 : G L 11 3 = G L 11 2 + inc L 11 2 := G_succ hL hp 2
  rcases le_or_gt 1 a with h | h
  · have up : G L 11 a ≤ G L 11 1 := G_le_up hL hp (by linarith) h
    refine ⟨by linarith, fun hc => ?_⟩
    by_contra hne
    have h3 : 3 ≤ a := by omega
    have i3 : inc L 11 3 ≤ inc L 11 1 := inc_anti hL hp (by norm_num)
    have up3 : G L 11 a ≤ G L 11 3 := G_le_up hL hp (by linarith) h3
    linarith
  · have ha : a = 0 := by omega
    subst ha
    exact ⟨by linarith, fun hc => by exfalso; linarith⟩

lemma low11 (hL : LogLike L) : 89 / 100000 < G L 11 1 - G L 11 2 := by
  have g2 : G L 11 2 = G L 11 1 + inc L 11 1 := G_succ hL (by norm_num) 1
  have := inc_11_1_le hL
  linarith

lemma inc_mid0 (hL : LogLike L) {p : ℕ} (h13 : 13 ≤ p) (h61 : p ≤ 61) :
    13 / 10000 < inc L p 0 := by
  have hs : σ 1 (p ^ 0) = 1 := by simp [isMultiplicative_sigma.map_one]
  rw [inc_val p 0 1 hs]
  have hp0 : (13 : ℝ) ≤ p := by exact_mod_cast h13
  have hp1 : (p : ℝ) ≤ 61 := by exact_mod_cast h61
  have e1 : (1 : ℝ) + 1 / ((p : ℝ) * ((1 : ℕ) : ℝ)) = 1 + 1 / (p : ℝ) := by simp
  rw [e1]
  have h1 : L (62 / 61 : ℝ) ≤ L (1 + 1 / (p : ℝ)) := by
    refine hL.mono _ _ (by norm_num) ?_
    have : 1 / (61 : ℝ) ≤ 1 / (p : ℝ) := one_div_le_one_div_of_le (by linarith) hp1
    linarith
  have h2 : L (p : ℝ) ≤ L 61 := hL.mono _ _ (by linarith) hp1
  have h2' : L 13 ≤ L (p : ℝ) := hL.mono _ _ (by norm_num) hp0
  have h3 : eps L * L (p : ℝ) ≤ eps L * L 61 := mul_le_mul_of_nonneg_left h2 (eps_pos hL).le
  have h4 := mul_le_mul (eps_le hL) (L61_le hL) (by linarith [L61_ge hL]) (by norm_num)
  linarith [L_u61_ge hL]

lemma inc_mid1 (hL : LogLike L) {p : ℕ} (hp : p.Prime) (h13 : 13 ≤ p) :
    inc L p 1 < -(13 / 10000) := by
  have hs : σ 1 (p ^ 1) = p + 1 := by
    rw [sigma_one_apply_prime_pow hp, Finset.sum_range_succ, Finset.sum_range_one]
    ring
  rw [inc_val p 1 (p + 1) hs]
  have hp0 : (13 : ℝ) ≤ p := by exact_mod_cast h13
  have h1 : L (1 + 1 / ((p : ℝ) * ((p + 1 : ℕ) : ℝ))) ≤ L (183 / 182 : ℝ) := by
    refine hL.mono _ _ (by positivity) ?_
    have hq : (182 : ℝ) ≤ (p : ℝ) * ((p + 1 : ℕ) : ℝ) := by push_cast; nlinarith
    have : 1 / ((p : ℝ) * ((p + 1 : ℕ) : ℝ)) ≤ 1 / 182 := one_div_le_one_div_of_le (by norm_num) hq
    linarith
  have h2' : L 13 ≤ L (p : ℝ) := hL.mono _ _ (by norm_num) hp0
  have h4 := mul_le_mul (eps_ge hL) ((L13_ge hL).trans h2') (by norm_num) (eps_pos hL).le
  linarith [L_u182_le hL]

lemma pmid (hL : LogLike L) {p : ℕ} (hp : p.Prime) (h13 : 13 ≤ p) (h61 : p ≤ 61) (a : ℕ) :
    0 ≤ G L p 1 - G L p a ∧ (G L p 1 - G L p a < 13 / 10000 → a = 1) := by
  have i0 := inc_mid0 hL h13 h61
  have i1 := inc_mid1 hL hp h13
  have g1 : G L p 1 = G L p 0 + inc L p 0 := G_succ hL hp 0
  have g2 : G L p 2 = G L p 1 + inc L p 1 := G_succ hL hp 1
  rcases le_or_gt 1 a with h | h
  · have up : G L p a ≤ G L p 1 := G_le_up hL hp (by linarith) h
    refine ⟨by linarith, fun hc => ?_⟩
    by_contra hne
    have h2 : 2 ≤ a := by omega
    have i2 : inc L p 2 ≤ inc L p 1 := inc_anti hL hp (by norm_num)
    have up2 : G L p a ≤ G L p 2 := G_le_up hL hp (by linarith) h2
    linarith
  · have ha : a = 0 := by omega
    subst ha
    exact ⟨by linarith, fun hc => by exfalso; linarith⟩

lemma p67 (hL : LogLike L) (a : ℕ) :
    0 ≤ G L 67 0 - G L 67 a ∧ (G L 67 0 - G L 67 a < 13 / 10000 → a = 0 ∨ a = 1) := by
  have hp : Nat.Prime 67 := by norm_num
  have i0 := inc_67_0 hL
  have i1 := inc_67_1_le hL
  have g1 : G L 67 1 = G L 67 0 + inc L 67 0 := G_succ hL hp 0
  have g2 : G L 67 2 = G L 67 1 + inc L 67 1 := G_succ hL hp 1
  have up : G L 67 a ≤ G L 67 0 := G_le_up hL hp (by linarith) (Nat.zero_le a)
  refine ⟨by linarith, fun hc => ?_⟩
  by_contra hne
  have h2 : 2 ≤ a := by omega
  have i2 : inc L 67 2 ≤ inc L 67 0 := inc_anti hL hp (by norm_num)
  have up2 : G L 67 a ≤ G L 67 2 := G_le_up hL hp (by linarith) h2
  linarith

lemma p71 (hL : LogLike L) (a : ℕ) :
    0 ≤ G L 71 0 - G L 71 a ∧ (G L 71 0 - G L 71 a < 13 / 10000 → a = 0 ∨ a = 1) := by
  have hp : Nat.Prime 71 := by norm_num
  have i0 := inc_71_0_le hL
  have i1 := inc_71_1_le hL
  have g1 : G L 71 1 = G L 71 0 + inc L 71 0 := G_succ hL hp 0
  have g2 : G L 71 2 = G L 71 1 + inc L 71 1 := G_succ hL hp 1
  have up : G L 71 a ≤ G L 71 0 := G_le_up hL hp (by linarith) (Nat.zero_le a)
  refine ⟨by linarith, fun hc => ?_⟩
  by_contra hne
  have h2 : 2 ≤ a := by omega
  have i2 : inc L 71 2 ≤ inc L 71 0 := inc_anti hL hp (by norm_num)
  have up2 : G L 71 a ≤ G L 71 2 := G_le_up hL hp (by linarith) h2
  linarith

lemma low71 (hL : LogLike L) : 1 / 1000 < G L 71 0 - G L 71 1 := by
  have g1 : G L 71 1 = G L 71 0 + inc L 71 0 := G_succ hL (by norm_num) 0
  have := inc_71_0_le hL
  linarith

lemma inc_big0 (hL : LogLike L) {p : ℕ} (h73 : 73 ≤ p) : inc L p 0 < -(13 / 10000) := by
  have hs : σ 1 (p ^ 0) = 1 := by simp [isMultiplicative_sigma.map_one]
  rw [inc_val p 0 1 hs]
  have hp0 : (73 : ℝ) ≤ p := by exact_mod_cast h73
  have e1 : (1 : ℝ) + 1 / ((p : ℝ) * ((1 : ℕ) : ℝ)) = 1 + 1 / (p : ℝ) := by simp
  rw [e1]
  have h1 : L (1 + 1 / (p : ℝ)) ≤ L (74 / 73 : ℝ) := by
    refine hL.mono _ _ (by positivity) ?_
    have : 1 / (p : ℝ) ≤ 1 / 73 := one_div_le_one_div_of_le (by norm_num) hp0
    linarith
  have h2' : L 73 ≤ L (p : ℝ) := hL.mono _ _ (by norm_num) hp0
  have h4 := mul_le_mul (eps_ge hL) ((L73_ge hL).trans h2') (by norm_num) (eps_pos hL).le
  linarith [L_u73_le hL]

lemma pbig (hL : LogLike L) {p : ℕ} (hp : p.Prime) (h73 : 73 ≤ p) (a : ℕ) :
    0 ≤ G L p 0 - G L p a ∧ (G L p 0 - G L p a < 13 / 10000 → a = 0) := by
  have i0 := inc_big0 hL h73
  have g1 : G L p 1 = G L p 0 + inc L p 0 := G_succ hL hp 0
  have up : G L p a ≤ G L p 0 := G_le_up hL hp (by linarith) (Nat.zero_le a)
  refine ⟨by linarith, fun hc => ?_⟩
  by_contra hne
  have h1 : 1 ≤ a := by omega
  have i1 : inc L p 1 ≤ inc L p 0 := inc_anti hL hp (by norm_num)
  have up1 : G L p a ≤ G L p 1 := G_le_up hL hp (by linarith) h1
  linarith


end OddRobin
