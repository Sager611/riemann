import Riemann.Deep.OddRobin.Lean.Final

namespace OddRobin

open ArithmeticFunction
open scoped ArithmeticFunction.sigma

variable {L : ℝ → ℝ}

/-! ## (1d) log log log is steeper than eps on [log N, log c0] -/

lemma Dstep (hL : LogLike L) {x X T T2 e : ℝ} (hX : 1 < L X) (hXx : L X ≤ L x) (hT : L x ≤ T)
    (hT2 : L T ≤ T2) (he0 : 0 ≤ e) (he : e * (T * T2) ≤ 1) :
    L (L (L X)) + (e * L x - e * L X) ≤ L (L (L x)) := by
  have hm1 : 0 < L X := by linarith
  have hl1 : 0 < L x := by linarith
  have hm2 : 0 < L (L X) := hL.pos_of_one_lt hX
  have hl2m : L (L X) ≤ L (L x) := hL.mono _ _ hm1 hXx
  have hl2 : 0 < L (L x) := by linarith
  have hT0 : 0 < T := by linarith
  have hl2T : L (L x) ≤ T2 := (hL.mono _ _ hl1 hT).trans hT2
  have s1 : (L x - L X) / L x ≤ L (L x) - L (L X) := by
    have h := hL.one_sub_inv_le (L x / L X) (div_pos hl1 hm1)
    rw [hL.map_div hl1 hm1, inv_div] at h
    rw [sub_div, div_self hl1.ne']
    linarith
  have s2 : (L (L x) - L (L X)) / L (L x) ≤ L (L (L x)) - L (L (L X)) := by
    have h := hL.one_sub_inv_le (L (L x) / L (L X)) (div_pos hl2 hm2)
    rw [hL.map_div hl2 hm2, inv_div] at h
    rw [sub_div, div_self hl2.ne']
    linarith
  have hd : 0 ≤ L x - L X := by linarith
  have h0 : 0 ≤ L (L x) - L (L X) := by linarith
  have s1' : L x - L X ≤ (L (L x) - L (L X)) * L x := by
    rw [div_le_iff₀ hl1] at s1; linarith
  have hprod : L x * L (L x) ≤ T * T2 := mul_le_mul hT hl2T hl2.le hT0.le
  have hT20 : 0 < T2 := by linarith
  have hTT : 0 < T * T2 := mul_pos hT0 hT20
  have s3 : (L x - L X) / (T * T2) ≤ (L (L x) - L (L X)) / L (L x) := by
    rw [div_le_div_iff₀ hTT hl2]
    nlinarith [mul_le_mul_of_nonneg_right s1' hl2.le, mul_le_mul_of_nonneg_left hprod h0]
  have s4 : e * (L x - L X) ≤ (L x - L X) / (T * T2) := by
    rw [le_div_iff₀ hTT]
    nlinarith [mul_le_mul_of_nonneg_left he hd]
  linarith

/-! ## (1b)–(1d): the total cost of a violator is below 13/10000 -/

lemma total_cost (hL : LogLike L) {E g : ℝ} (hE : 0 < E) (hLE : L E = g)
    (hg : (5772145 / 10000000 : ℝ) ≤ g) {n : ℕ} (hn0 : n ≠ 0) (h1 : NN < n) (h2 : n < c0)
    (hv : ¬ ORf L E n) :
    ∑ p ∈ n.primeFactors ∪ P61, (G L p (NN.factorization p) - G L p (n.factorization p))
      < 13 / 10000 := by
  rw [Finset.sum_sub_distrib, decomp hL NN_ne_zero (NN_pf_sub.trans Finset.subset_union_right),
    decomp hL hn0 Finset.subset_union_left]
  have hN0 : (0 : ℝ) < (NN : ℝ) := by norm_num [NN]
  have hn' : (NN : ℝ) < n := by exact_mod_cast h1
  have hc' : (n : ℝ) < c0 := by exact_mod_cast h2
  have hl1 : L (NN : ℝ) ≤ L n := hL.mono _ _ hN0 hn'.le
  have hl2 : L (n : ℝ) ≤ L c0 := hL.mono _ _ (by linarith) hc'.le
  have hNg := LNN_gt_one hL
  have hll : 0 < L (L (n : ℝ)) := hL.pos_of_one_lt (by linarith)
  have hv' : L E - L 2 + L (L (L (n : ℝ))) ≤ L ((σ 1 n : ℝ) / n) :=
    not_lt.1 (fun h => hv ((ORf_iff_log hL hE (Nat.pos_of_ne_zero hn0) hll).2 h))
  have hσn : (0 : ℝ) < (σ 1 n : ℝ) := by exact_mod_cast sigma_pos_iff.2 (Nat.pos_of_ne_zero hn0)
  have hσN : (0 : ℝ) < (σ 1 NN : ℝ) := by
    exact_mod_cast sigma_pos_iff.2 (Nat.pos_of_ne_zero NN_ne_zero)
  rw [hL.map_div hσn (by linarith), hLE] at hv'
  have hD := Dstep hL hNg hl1 (hl2.trans (Lc0_le hL)) (LT_le hL) (eps_pos hL).le (eps_TT2 hL)
  have hδ := delta0_lt hL hg
  rw [hL.map_div hσN hN0] at hδ
  have e : (L (σ 1 NN : ℝ) - (1 + eps L) * L NN) - (L (σ 1 n : ℝ) - (1 + eps L) * L n)
      = (L (σ 1 NN : ℝ) - L NN) - (L (σ 1 n : ℝ) - L n) + (eps L * L n - eps L * L NN) := by ring
  rw [e]
  linarith


end OddRobin
