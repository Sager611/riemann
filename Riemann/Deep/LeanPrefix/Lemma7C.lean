import Riemann.Deep.LeanPrefix.Lemma7B

/-!
# Lemma 7 (involution form): a nonidentity involution of PΓL_d(q) on the points of
PG(d-1, q), d ≥ 2, (d, q) ≠ (2, 4), fixes at most half of the points

NOTE Lemma 7 (here for involutions, which is what Proposition 8 uses for r_k).
f² = μ·id and σ² = id (Lemma7A). Linear case: every fixed point lies in V_λ or V_{-λ};
|V_±| ≤ q^{d-1}, |V_λ| |V_{-λ}| ≤ q^d. Semilinear case: q = q₀², every fixed point has a
representative in W = {v : f v = λ₀ v} (quadratic Hilbert 90), |W|² ≤ q^d.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv
open scoped LinearAlgebra.Projectivization

theorem arith_lin {q Q P : ℕ} (hq : 2 ≤ q) (hQ : q ≤ Q)
    (h : P + 1 ≤ Q ∨ (3 ≤ q ∧ P + 2 ≤ Q + q)) : 2 * P + 1 ≤ q * Q := by
  rcases h with h | ⟨h3, h⟩
  · have : 2 * Q ≤ q * Q := Nat.mul_le_mul_right Q hq
    omega
  · have h1 : 3 * Q ≤ q * Q := Nat.mul_le_mul_right Q h3
    have h2 : (q - 2) * q ≤ (q - 2) * Q := Nat.mul_le_mul_left _ hQ
    have e1 : (q - 2) * Q + 2 * Q = q * Q := by
      rw [← Nat.add_mul, Nat.sub_add_cancel (by omega)]
    have e2 : (q - 2) * q + 2 * q = q * q := by
      rw [← Nat.add_mul, Nat.sub_add_cancel (by omega)]
    have h4 : 3 * q ≤ q * q := Nat.mul_le_mul_right q h3
    nlinarith

theorem arith_semi {q0 d P W : ℕ} (hq0 : 2 ≤ q0) (hd : 2 ≤ d) (h22 : ¬ (d = 2 ∧ q0 = 2))
    (hW : W * W ≤ q0 ^ d * q0 ^ d) (hP : P + 1 ≤ W) :
    2 * (q0 + 1) * P + 1 ≤ q0 ^ d * q0 ^ d := by
  have hR : W ≤ q0 ^ d := Nat.mul_self_le_mul_self_iff.mp hW
  have hR2 : 2 * q0 + 1 ≤ q0 ^ d := by
    rcases Nat.lt_or_ge d 3 with h | h
    · have hd2 : d = 2 := by omega
      have hq3 : 3 ≤ q0 := by
        by_contra hc; exact h22 ⟨hd2, by omega⟩
      rw [hd2, pow_two]
      have := Nat.mul_le_mul_right q0 hq3
      omega
    · have h1 : q0 ^ 3 ≤ q0 ^ d := Nat.pow_le_pow_right (by omega) h
      have h2 : 4 * q0 ≤ q0 ^ 3 := by
        rw [show q0 ^ 3 = q0 * q0 * q0 by ring]
        have : 4 ≤ q0 * q0 := by nlinarith
        exact Nat.mul_le_mul_right q0 this
      omega
  have h3 : 2 * (q0 + 1) * P ≤ (q0 ^ d + 1) * (q0 ^ d - 1) :=
    Nat.mul_le_mul (by omega) (by omega)
  have h4 : (q0 ^ d + 1) * (q0 ^ d - 1) + 1 = q0 ^ d * q0 ^ d := by
    have : 1 ≤ q0 ^ d := Nat.one_le_pow _ _ (by omega)
    obtain ⟨R, hR'⟩ : ∃ R, q0 ^ d = R + 1 := ⟨q0 ^ d - 1, by omega⟩
    rw [hR', Nat.add_sub_cancel]; ring
  omega

section

variable {F : Type} [Field F] [Fintype F] [DecidableEq F] {d n : ℕ}

theorem inj_smul_points (e : Fin n ≃ ℙ F (Fin d → F)) (Fx : Finset (Fin n)) (C : Finset F)
    (hC : ∀ c ∈ C, c ≠ 0) (vx : Fin n → Fin d → F)
    (hvx : ∀ x ∈ Fx, vx x ≠ 0 ∧ vx x ∈ (e x).submodule) :
    Set.InjOn (fun p : Fin n × F => p.2 • vx p.1) (Fx ×ˢ C : Finset (Fin n × F)) := by
  intro p hp p' hp' hpp
  simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe] at hp hp'
  simp only at hpp
  have h0 : p.2 • vx p.1 ≠ 0 := smul_ne_zero (hC _ hp.2) (hvx _ hp.1).1
  have hx : p.1 = p'.1 := by
    apply e.injective
    apply point_eq_of_mem h0 (smul_mem_point _ (hvx _ hp.1).2)
    rw [hpp]; exact smul_mem_point _ (hvx _ hp'.1).2
  have hc : p.2 = p'.2 := by
    rw [hx] at hpp
    have : (p.2 - p'.2) • vx p'.1 = 0 := by rw [sub_smul, hpp, sub_self]
    rcases smul_eq_zero.mp this with h | h
    · exact sub_eq_zero.mp h
    · exact absurd h (hvx _ hp'.1).1
  exact Prod.ext hx hc

theorem card_nonzero : (Finset.univ.filter (fun c : F => c ≠ 0)).card = Fintype.card F - 1 := by
  rw [Finset.filter_ne' Finset.univ (0 : F), Finset.card_erase_of_mem (Finset.mem_univ _),
    Finset.card_univ]

/-- Lemma 7, involution form. -/
theorem lemma7_invol (hd : 2 ≤ d) (e : Fin n ≃ ℙ F (Fin d → F)) {g : Perm (Fin n)}
    (hg2 : ∀ x, g (g x) = x) (hg1 : g ≠ 1) (hcoll : IsCollineationVia e g)
    (h24 : ¬ (d = 2 ∧ Fintype.card F = 4)) : 2 * nfix g ≤ n := by
  obtain ⟨σ, f, hσ, hf⟩ := hcoll
  obtain ⟨μ, hμ0, hμ, hσ2⟩ := f_sq_scalar hd hσ hf hg2
  have hq2 : 2 ≤ Fintype.card F := Fintype.one_lt_card
  have hcard : Fintype.card F ^ d - 1 = n * (Fintype.card F - 1) := by
    have := Projectivization.card F (Fin d → F)
    rw [Nat.card_eq_fintype_card, card_vec, Nat.card_congr e.symm, Nat.card_eq_fintype_card,
      Fintype.card_fin, Nat.card_eq_fintype_card] at this
    exact this
  obtain ⟨Q, hQ⟩ : ∃ Q, Q = Fintype.card F ^ (d - 1) := ⟨_, rfl⟩
  have hqQ : Fintype.card F ^ d = Fintype.card F * Q := by
    rw [hQ, ← pow_succ', Nat.sub_add_cancel (by omega)]
  have hQq : Fintype.card F ≤ Q := by
    rw [hQ]; calc Fintype.card F = Fintype.card F ^ 1 := (pow_one _).symm
      _ ≤ _ := Nat.pow_le_pow_right (by omega) (by omega)
  obtain ⟨Fx, hFx⟩ : ∃ Fx : Finset (Fin n), Fx = Finset.univ.filter (fun x => g x = x) := ⟨_, rfl⟩
  have hnfix : nfix g = Fx.card := by rw [hFx]; rfl
  have hfixed : ∀ x ∈ Fx, g x = x := fun x hx => by rw [hFx] at hx; exact (Finset.mem_filter.mp hx).2
  have hev : ∀ x : Fin n, ∃ lam : F, g x = x → lam ≠ 0 ∧ f (e x).rep = lam • (e x).rep := by
    intro x
    by_cases hx : g x = x
    · obtain ⟨lam, h1, h2⟩ := fixed_eigen hf hx
      exact ⟨lam, fun _ => ⟨h1, h2⟩⟩
    · exact ⟨0, fun h => absurd h hx⟩
  choose lam hlam using hev
  have hnotscalar : ∀ c : F, ∃ v, f v ≠ c • v := by
    intro c
    by_contra hc
    push Not at hc
    apply hg1
    ext x
    simp only [Perm.one_apply]
    congr 1
    apply e.injective
    exact point_eq_of_mem (f_ne_zero (f := f) (Projectivization.rep_nonzero (e x)))
      (hf x _ (rep_mem (e x))) (by rw [hc]; exact smul_mem_point c (rep_mem (e x)))
  have hrep : ∀ x ∈ Fx, (e x).rep ≠ 0 ∧ (e x).rep ∈ (e x).submodule :=
    fun x _ => ⟨Projectivization.rep_nonzero _, rep_mem _⟩
  -- the eigen-sets S c = {v : f v = c • v}
  have hSsub : ∀ c : F, ∀ a ∈ Finset.univ.filter (fun v => f v = c • v),
      ∀ b ∈ Finset.univ.filter (fun v => f v = c • v), a - b ∈ Finset.univ.filter (fun v => f v = c • v) := by
    intro c a ha b hb
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha hb ⊢
    rw [map_sub, ha, hb, smul_sub]
  rcases Finset.eq_empty_or_nonempty Fx with hE | ⟨x1, hx1⟩
  · rw [hnfix, hE, Finset.card_empty]; omega
  have hx1' := hfixed x1 hx1
  by_cases hlin : ∀ c, σ c = c
  · ------------------------------------------------------------------ linear case
    have hSsmul : ∀ c : F, ∀ (t : F), ∀ a ∈ Finset.univ.filter (fun v => f v = c • v),
        t • a ∈ Finset.univ.filter (fun v => f v = c • v) := by
      intro c t a ha
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha ⊢
      rw [hσ, hlin, ha, smul_comm]
    obtain ⟨l, hl⟩ : ∃ l, l = lam x1 := ⟨_, rfl⟩
    have hsq : ∀ x, g x = x → lam x * lam x = μ := by
      intro x hx
      have h1 := hμ (e x).rep
      rw [(hlam x hx).2, hσ, hlin, (hlam x hx).2, smul_smul] at h1
      have h2 : (lam x * lam x - μ) • (e x).rep = 0 := by rw [sub_smul, h1, sub_self]
      rcases smul_eq_zero.mp h2 with h | h
      · exact sub_eq_zero.mp h
      · exact absurd h (Projectivization.rep_nonzero _)
    have hpm : ∀ x, g x = x → lam x = l ∨ lam x = -l := by
      intro x hx
      have : (lam x - l) * (lam x + l) = 0 := by
        rw [hl]; linear_combination hsq x hx - hsq x1 hx1'
      rcases mul_eq_zero.mp this with h | h
      · exact Or.inl (sub_eq_zero.mp h)
      · exact Or.inr (eq_neg_of_add_eq_zero_left h)
    obtain ⟨Sp, hSp⟩ : ∃ S : Finset (Fin d → F), S = Finset.univ.filter (fun v => f v = l • v) := ⟨_, rfl⟩
    obtain ⟨Sm, hSm⟩ : ∃ S : Finset (Fin d → F), S = Finset.univ.filter (fun v => f v = (-l) • v) := ⟨_, rfl⟩
    have h0p : (0 : Fin d → F) ∈ Sp := by rw [hSp]; simp
    have h0m : (0 : Fin d → F) ∈ Sm := by rw [hSm]; simp
    have hpq : Sp.card * Fintype.card F ≤ Fintype.card F * Q := by
      obtain ⟨w, hw⟩ := hnotscalar l
      have := card_mul_q_le Sp (by rw [hSp]; exact hSsub l) (by rw [hSp]; exact hSsmul l)
        (w := w) (by rw [hSp]; simpa using hw)
      rw [card_vec, hqQ] at this; exact this
    have hmq : Sm.card * Fintype.card F ≤ Fintype.card F * Q := by
      obtain ⟨w, hw⟩ := hnotscalar (-l)
      have := card_mul_q_le Sm (by rw [hSm]; exact hSsub (-l)) (by rw [hSm]; exact hSsmul (-l))
        (w := w) (by rw [hSm]; simpa using hw)
      rw [card_vec, hqQ] at this; exact this
    have hpQ : Sp.card ≤ Q := by
      rw [Nat.mul_comm] at hpq; exact Nat.le_of_mul_le_mul_left hpq (by omega)
    have hmQ : Sm.card ≤ Q := by
      rw [Nat.mul_comm] at hmq; exact Nat.le_of_mul_le_mul_left hmq (by omega)
    have hrep1 : (e x1).rep ∈ Sp := by rw [hSp]; simp [(hlam x1 hx1').2, hl]
    have hpq' : Fintype.card F ≤ Sp.card :=
      q_le_card Sp (by rw [hSp]; exact hSsmul l) hrep1 (Projectivization.rep_nonzero _)
    -- the injection
    have hmaps : Set.MapsTo (fun p : Fin n × F => p.2 • (e p.1).rep)
        (Fx ×ˢ Finset.univ.filter (fun c : F => c ≠ 0) : Finset (Fin n × F)) (Sp.erase 0 ∪ Sm.erase 0 : Finset _) := by
      intro p hp
      simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, Finset.mem_filter,
        Finset.mem_univ, true_and] at hp
      have hx := hfixed _ hp.1
      have h0 : p.2 • (e p.1).rep ≠ 0 := smul_ne_zero hp.2 (Projectivization.rep_nonzero _)
      have hfv : f (p.2 • (e p.1).rep) = lam p.1 • (p.2 • (e p.1).rep) := by
        rw [hσ, hlin, (hlam _ hx).2, smul_comm]
      simp only [Finset.coe_union, Set.mem_union, Finset.mem_coe, Finset.mem_erase]
      rcases hpm _ hx with h | h
      · left; refine ⟨h0, ?_⟩; rw [hSp]; simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        rw [hfv, h]
      · right; refine ⟨h0, ?_⟩; rw [hSm]; simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        rw [hfv, h]
    have hinj := Finset.card_le_card_of_injOn _ hmaps
      (inj_smul_points e Fx _ (fun c hc => (Finset.mem_filter.mp hc).2) (fun x => (e x).rep) hrep)
    rw [Finset.card_product, card_nonzero] at hinj
    have hP : Fx.card * (Fintype.card F - 1) + 1 ≤ Q ∨
        (3 ≤ Fintype.card F ∧ Fx.card * (Fintype.card F - 1) + 2 ≤ Q + Fintype.card F) := by
      by_cases hll : l = -l
      · left
        have hSS : Sm = Sp := by rw [hSm, hSp, ← hll]
        rw [hSS, Finset.union_self, Finset.card_erase_of_mem h0p] at hinj
        omega
      · -- two eigen-sets, Sp ∩ Sm = {0}
        have hl0 : l ≠ 0 := by rw [hl]; exact (hlam x1 hx1').1
        have h3 : 3 ≤ Fintype.card F := by
          have hs : ({0, l, -l} : Finset F).card = 3 := by
            rw [Finset.card_insert_of_notMem, Finset.card_pair hll]
            simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
            exact ⟨Ne.symm hl0, by intro h; apply hl0; exact neg_eq_zero.mp h.symm⟩
          rw [← hs]; exact Finset.card_le_univ _
        have hcap : ∀ v ∈ Sp, v ∈ Sm → v = 0 := by
          intro v hvp hvm
          rw [hSp] at hvp; rw [hSm] at hvm
          simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hvp hvm
          have : (l - -l) • v = 0 := by rw [sub_smul, ← hvp, ← hvm, sub_self]
          rcases smul_eq_zero.mp this with h | h
          · exact absurd (sub_eq_zero.mp h) hll
          · exact h
        have hprod : Sp.card * Sm.card ≤ Fintype.card F * Q := by
          have := card_mul_le_of_add_inj Sp Sm 1 (by
            intro a ha b hb a' ha' b' hb' hab
            rw [one_smul, one_smul] at hab
            have h1 : a - a' = b' - b := by rw [sub_eq_sub_iff_add_eq_add, hab, add_comm]
            have hma : a - a' ∈ Sp := by rw [hSp] at ha ha' ⊢; exact hSsub l _ ha _ ha'
            have hmb : b' - b ∈ Sm := by rw [hSm] at hb hb' ⊢; exact hSsub (-l) _ hb' _ hb
            rw [h1] at hma
            have h2 := hcap _ hma hmb
            rw [h2] at h1
            exact ⟨sub_eq_zero.mp h1, (sub_eq_zero.mp h2).symm⟩)
          rw [card_vec, hqQ] at this; exact this
        rcases Nat.lt_or_ge Sm.card 2 with hm1 | hm2
        · left
          have hm1' : Sm.card = 1 := by
            have := Finset.card_pos.mpr ⟨0, h0m⟩; omega
          have := Finset.card_union_le (Sp.erase 0) (Sm.erase 0)
          rw [Finset.card_erase_of_mem h0p, Finset.card_erase_of_mem h0m, hm1'] at this
          omega
        · right
          refine ⟨h3, ?_⟩
          obtain ⟨v, hv, hv0⟩ : ∃ v ∈ Sm, v ≠ 0 := by
            by_contra hc
            push Not at hc
            have : Sm ⊆ {0} := fun v hv => Finset.mem_singleton.mpr (hc v hv)
            have := Finset.card_le_card this
            rw [Finset.card_singleton] at this; omega
          have hmq' : Fintype.card F ≤ Sm.card :=
            q_le_card Sm (by rw [hSm]; exact hSsmul (-l)) hv hv0
          have := Finset.card_union_le (Sp.erase 0) (Sm.erase 0)
          rw [Finset.card_erase_of_mem h0p, Finset.card_erase_of_mem h0m] at this
          have hsum : Fintype.card F * (Sp.card + Sm.card) ≤ Fintype.card F * (Q + Fintype.card F) := by
            have h5 : (Sp.card - Fintype.card F) * (Sm.card - Fintype.card F) ≥ 0 := Nat.zero_le _
            nlinarith
          have hsum' : Sp.card + Sm.card ≤ Q + Fintype.card F :=
            Nat.le_of_mul_le_mul_left hsum (by omega)
          omega
    have := arith_lin hq2 hQq hP
    have h6 : 2 * (Fx.card * (Fintype.card F - 1)) + 1 ≤ n * (Fintype.card F - 1) + 1 := by
      rw [← hcard, hqQ]; omega
    rw [hnfix]
    have h7 : (2 * Fx.card) * (Fintype.card F - 1) ≤ n * (Fintype.card F - 1) := by
      rw [Nat.mul_assoc]; omega
    exact Nat.le_of_mul_le_mul_right h7 (by omega)
  · ------------------------------------------------------------------ semilinear case
    push Not at hlin
    obtain ⟨ω, hω⟩ := hlin
    obtain ⟨F0, hF0⟩ : ∃ F0 : Finset F, F0 = Finset.univ.filter (fun c : F => σ c = c) := ⟨_, rfl⟩
    have hq : Fintype.card F = F0.card ^ 2 := by rw [hF0]; exact card_eq_sq_fixed hσ2 hω
    have hq02 : 2 ≤ F0.card := by
      have hs : ({0, 1} : Finset F) ⊆ F0 := by
        intro c hc; rw [hF0]; simp only [Finset.mem_insert, Finset.mem_singleton] at hc
        rcases hc with rfl | rfl <;> simp
      have := Finset.card_le_card hs
      rwa [Finset.card_pair zero_ne_one] at this
    obtain ⟨l0, hl0⟩ : ∃ l, l = lam x1 := ⟨_, rfl⟩
    have hnorm : ∀ x, g x = x → lam x * σ (lam x) = μ := by
      intro x hx
      have h1 := hμ (e x).rep
      rw [(hlam x hx).2, hσ, (hlam x hx).2, smul_smul] at h1
      have h2 : (σ (lam x) * lam x - μ) • (e x).rep = 0 := by rw [sub_smul, h1, sub_self]
      rcases smul_eq_zero.mp h2 with h | h
      · rw [mul_comm]; exact sub_eq_zero.mp h
      · exact absurd h (Projectivization.rep_nonzero _)
    have hW : ∀ x : Fin n, ∃ w : Fin d → F,
        g x = x → w ≠ 0 ∧ w ∈ (e x).submodule ∧ f w = l0 • w := by
      intro x
      by_cases hx : g x = x
      · have hlx := (hlam x hx).1
        obtain ⟨c, hc0, hc⟩ := hilbert90_quad hσ2 hω (ν := l0 * (lam x)⁻¹) (by
          rw [map_mul, map_inv₀]
          have h1 := hnorm x hx
          have h2 := hnorm x1 hx1'
          rw [← hl0] at h2
          have hσl : σ (lam x) ≠ 0 := by rw [Ne, map_eq_zero]; exact hlx
          field_simp
          rw [← h1] at h2
          linear_combination h2)
        refine ⟨c • (e x).rep, fun _ => ⟨smul_ne_zero hc0 (Projectivization.rep_nonzero _),
          smul_mem_point c (rep_mem _), ?_⟩⟩
        rw [hσ, hc, (hlam x hx).2, smul_smul, smul_smul]
        congr 1
        field_simp
      · exact ⟨0, fun h => absurd h hx⟩
    choose w hw using hW
    obtain ⟨W, hWdef⟩ : ∃ W : Finset (Fin d → F), W = Finset.univ.filter (fun v => f v = l0 • v) := ⟨_, rfl⟩
    have h0W : (0 : Fin d → F) ∈ W := by rw [hWdef]; simp
    have hmaps : Set.MapsTo (fun p : Fin n × F => p.2 • w p.1)
        (Fx ×ˢ F0.filter (fun c : F => c ≠ 0) : Finset (Fin n × F)) (W.erase 0 : Finset _) := by
      intro p hp
      simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, Finset.mem_filter] at hp
      have hx := hfixed _ hp.1
      have hσp : σ p.2 = p.2 := by rw [hF0] at hp; exact (Finset.mem_filter.mp hp.2.1).2
      simp only [Finset.coe_erase, Set.mem_diff, Finset.mem_coe, Set.mem_singleton_iff]
      refine ⟨?_, smul_ne_zero hp.2.2 (hw _ hx).1⟩
      rw [hWdef]; simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      rw [hσ, hσp, (hw _ hx).2.2, smul_comm]
    have hinj := Finset.card_le_card_of_injOn _ hmaps
      (inj_smul_points e Fx (F0.filter (fun c : F => c ≠ 0)) (fun c hc => (Finset.mem_filter.mp hc).2) w
        (fun x hx => ⟨(hw x (hfixed x hx)).1, (hw x (hfixed x hx)).2.1⟩))
    rw [Finset.card_product, Finset.filter_ne', Finset.card_erase_of_mem (by rw [hF0]; simp),
      Finset.card_erase_of_mem h0W] at hinj
    have hWW : W.card * W.card ≤ F0.card ^ d * F0.card ^ d := by
      have := card_mul_le_of_add_inj W W ω (by
        intro a ha b hb a' ha' b' hb' hab
        have h1 : a - a' = ω • (b' - b) := by
          rw [smul_sub, sub_eq_sub_iff_add_eq_add, hab, add_comm]
        have hma : a - a' ∈ W := by rw [hWdef] at ha ha' ⊢; exact hSsub l0 _ ha _ ha'
        have hmb : b' - b ∈ W := by rw [hWdef] at hb hb' ⊢; exact hSsub l0 _ hb' _ hb
        rw [hWdef] at hma hmb
        simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hma hmb
        rw [h1, hσ, hmb, smul_smul, smul_smul] at hma
        have h2 : (σ ω * l0 - l0 * ω) • (b' - b) = 0 := by rw [sub_smul, hma, sub_self]
        have hl00 : l0 ≠ 0 := by rw [hl0]; exact (hlam x1 hx1').1
        have h3 : σ ω * l0 - l0 * ω ≠ 0 := by
          rw [mul_comm, ← mul_sub]; exact mul_ne_zero hl00 (sub_ne_zero.mpr hω)
        rcases smul_eq_zero.mp h2 with h | h
        · exact absurd h h3
        · rw [h, smul_zero, sub_eq_zero] at h1
          exact ⟨h1, (sub_eq_zero.mp h).symm⟩)
      rw [card_vec, hq, ← pow_mul, show 2 * d = d + d by ring, pow_add] at this
      exact this
    have h22 : ¬ (d = 2 ∧ F0.card = 2) := by
      rintro ⟨h1, h2⟩; apply h24; refine ⟨h1, ?_⟩; rw [hq, h2]; norm_num
    have hPW : Fx.card * (F0.card - 1) + 1 ≤ W.card := by
      have := Finset.card_pos.mpr ⟨0, h0W⟩
      omega
    have := arith_semi hq02 hd h22 hWW hPW
    have hcard' : F0.card ^ d * F0.card ^ d - 1 = n * ((F0.card + 1) * (F0.card - 1)) := by
      have e1 : F0.card ^ d * F0.card ^ d = Fintype.card F ^ d := by
        rw [hq, ← pow_mul, ← pow_add, two_mul]
      have e2 : Fintype.card F - 1 = (F0.card + 1) * (F0.card - 1) := by
        rw [hq]
        obtain ⟨r, hr⟩ : ∃ r, F0.card = r + 1 := ⟨F0.card - 1, by omega⟩
        rw [hr]
        have : (r + 1) ^ 2 = (r + 1 + 1) * r + 1 := by ring
        rw [this, Nat.add_sub_cancel, Nat.add_sub_cancel]
      rw [e1, hcard, e2]
    rw [hnfix]
    have h7 : (2 * Fx.card) * ((F0.card + 1) * (F0.card - 1)) ≤ n * ((F0.card + 1) * (F0.card - 1)) := by
      have e1 : (2 * Fx.card) * ((F0.card + 1) * (F0.card - 1)) = 2 * (F0.card + 1) * (Fx.card * (F0.card - 1)) := by ring
      rw [e1]; omega
    exact Nat.le_of_mul_le_mul_right h7 (Nat.mul_pos (by omega) (by omega))

end

end Riemann.Deep.LeanPrefix

