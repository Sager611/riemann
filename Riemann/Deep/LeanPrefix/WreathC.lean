import Riemann.Deep.LeanPrefix.WreathB
import Riemann.Deep.LeanPrefix.Bridge

/-!
# Theorem B.3: n odd, α = k + 1 ∣ n: interval blocks and G ≅ D_{2α} ≀ D_{2q}, q = n/α

wreath_range: the image of D_{2α} ≀ D_{2q} under its imprimitive action, transported along the
block coordinates, is G. Hence G ≃* D_{2α} ≀ D_{2q} and |G| = (2α)^q · 2q. The intervals
[iα, (i+1)α) are blocks; G is not primitive.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv DihedralGroup

theorem wreath_range {n m k q α : ℕ} [NeZero q] [NeZero α] (hn : n = q * α) (hm : m = (q - 1) * α)
    (hk : k + 1 = α) (hq : 3 ≤ q) (hα : 3 ≤ α) (hodd : α % 2 = 1) :
    (wreathToPerm hn).range = G n m k := by
  have hN : wreathToPerm hn (wN α q) = rev n n := wreathToPerm_eq hn _ _ (block_rn hn hq hα)
  have hM : wreathToPerm hn (wM α q) = rev n m := wreathToPerm_eq hn _ _ (block_rm hn hm hq hα)
  have hK : wreathToPerm hn (wK α q) = rev n k := wreathToPerm_eq hn _ _ (block_rk hn hk hq hα)
  apply le_antisymm
  · rintro _ ⟨w, rfl⟩
    have hH := wreath_eq_top hq hα hodd ((G n m k).comap (wreathToPerm hn))
      (by show wreathToPerm hn (wN α q) ∈ G n m k; rw [hN]; exact rn_mem n m k)
      (by show wreathToPerm hn (wM α q) ∈ G n m k; rw [hM]; exact rm_mem n m k)
      (by show wreathToPerm hn (wK α q) ∈ G n m k; rw [hK]; exact rk_mem n m k)
    have : w ∈ (G n m k).comap (wreathToPerm hn) := by rw [hH]; trivial
    exact this
  · apply G_le
    · rw [← hN]; exact ⟨_, rfl⟩
    · rw [← hM]; exact ⟨_, rfl⟩
    · rw [← hK]; exact ⟨_, rfl⟩

theorem G_mulEquiv_wreath {n m k q α : ℕ} [NeZero q] [NeZero α] (hn : n = q * α)
    (hm : m = (q - 1) * α) (hk : k + 1 = α) (hq : 3 ≤ q) (hα : 3 ≤ α) (hodd : α % 2 = 1) :
    Nonempty (G n m k ≃* Wreath α q) :=
  ⟨(MulEquiv.subgroupCongr (wreath_range hn hm hk hq hα hodd).symm).trans
    (MonoidHom.ofInjective (wreathToPerm_injective hn hq hα)).symm⟩

/-- Fibres of a function b form blocks of any subgroup of fibreStab b. -/
theorem isBlockOf_fibre {n : ℕ} {β : Type} (b : Fin n → β) (H : Subgroup (Perm (Fin n)))
    (hH : H ≤ fibreStab b) (c : β) : IsBlockOf H {x | b x = c} := by
  intro g hg
  by_cases hex : ∃ x, b x = c ∧ b (g x) = c
  · obtain ⟨x, hx, hgx⟩ := hex
    left
    have hst := hH hg
    ext z
    constructor
    · rintro ⟨y, hy, rfl⟩
      have : b y = b x := by rw [hy, hx]
      show b (g y) = c
      rw [← hgx]; exact (hst y x).mp this
    · intro hz
      refine ⟨g⁻¹ z, ?_, by simp⟩
      show b (g⁻¹ z) = c
      have hgg : g (g⁻¹ z) = z := by simp
      have h1 : b (g (g⁻¹ z)) = b (g x) := by rw [hgg, hgx]; exact hz
      rw [← hx]; exact (hst _ x).mpr h1
  · right
    rw [Set.disjoint_left]
    rintro _ ⟨y, hy, rfl⟩ hgy
    exact hex ⟨y, hy, hgy⟩

/-- Theorem B.3 on the line: n odd, (k+1) ∣ n. -/
theorem line_B3 {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m) (hline : n = m + k + 1)
    (hodd : n % 2 = 1) (hdiv : (k + 1) ∣ n) :
    ¬ MulAction.IsPreprimitive (G n m k) (Fin n) ∧
    (∀ i : ℕ, MulAction.IsBlock (G n m k) {x : Fin n | x.val / (k + 1) = i}) ∧
    (∃ e : Fin n ≃ ZMod (n / (k + 1)) × ZMod (k + 1),
      (wreathPerm (k + 1) (n / (k + 1))).range.map e.symm.permCongrHom.toMonoidHom = G n m k) ∧
    Nonempty (G n m k ≃* Wreath (k + 1) (n / (k + 1))) ∧
    Nat.card (G n m k) = (2 * (k + 1)) ^ (n / (k + 1)) * (2 * (n / (k + 1))) := by
  obtain ⟨q, hq⟩ := hdiv
  have hqdef : n / (k + 1) = q := by rw [hq, Nat.mul_div_cancel_left _ (by omega)]
  rw [hqdef]
  have hq2 : 2 ≤ q := by
    rcases Nat.lt_or_ge q 2 with h | h
    · interval_cases q <;> omega
    · exact h
  have hqodd : q % 2 = 1 := by
    by_contra h
    have : 2 ∣ n := by rw [hq]; exact Dvd.dvd.mul_left (Nat.dvd_of_mod_eq_zero (by omega)) _
    omega
  have hq3 : 3 ≤ q := by omega
  have hαodd : (k + 1) % 2 = 1 := by
    by_contra h
    have : 2 ∣ n := by rw [hq]; exact Dvd.dvd.mul_right (Nat.dvd_of_mod_eq_zero (by omega)) _
    omega
  have : NeZero q := ⟨by omega⟩
  have : NeZero (k + 1) := ⟨by omega⟩
  have hn : n = q * (k + 1) := by rw [hq, Nat.mul_comm]
  have hm : m = (q - 1) * (k + 1) := by
    rw [Nat.sub_mul, Nat.one_mul, ← hn]; omega
  have hblocks : ∀ i : ℕ, MulAction.IsBlock (G n m k) {x : Fin n | x.val / (k + 1) = i} := by
    intro i
    rw [isBlock_iff_isBlockOf]
    exact isBlockOf_fibre (fun x : Fin n => x.val / (k + 1)) _ (G_le_intervals hline hkm ⟨q, hq⟩) i
  have hiso := G_mulEquiv_wreath hn hm rfl hq3 (by omega) hαodd
  have hperm : ∃ e : Fin n ≃ ZMod q × ZMod (k + 1),
      (wreathPerm (k + 1) q).range.map e.symm.permCongrHom.toMonoidHom = G n m k :=
    ⟨blockEquiv hn, by rw [← MonoidHom.range_comp]; exact wreath_range hn hm rfl hq3 (by omega) hαodd⟩
  refine ⟨?_, hblocks, hperm, hiso, ?_⟩
  · intro hp
    have htriv := hp.isTrivialBlock_of_isBlock (hblocks 0)
    rcases htriv with h | h
    · have h0 : (⟨0, by omega⟩ : Fin n) ∈ {x : Fin n | x.val / (k + 1) = 0} := by
        show (0 : ℕ) / (k + 1) = 0; exact Nat.div_eq_of_lt (by omega)
      have h1 : (⟨1, by omega⟩ : Fin n) ∈ {x : Fin n | x.val / (k + 1) = 0} := by
        show (1 : ℕ) / (k + 1) = 0; exact Nat.div_eq_of_lt (by omega)
      have := h h0 h1
      simp [Fin.ext_iff] at this
    · have hmem : (⟨n - 1, by omega⟩ : Fin n) ∈ {x : Fin n | x.val / (k + 1) = 0} := by
        rw [h]; trivial
      simp only [Set.mem_setOf_eq] at hmem
      have : k + 1 ≤ n - 1 := by omega
      have := (Nat.le_div_iff_mul_le (show 0 < k + 1 by omega)).mpr (show 1 * (k + 1) ≤ n - 1 by omega)
      omega
  · obtain ⟨e⟩ := hiso
    rw [Nat.card_congr e.toEquiv, card_wreath]

end Riemann.Deep.LeanPrefix

