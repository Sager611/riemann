import Riemann.Deep.LeanPrefix.Alt
import Riemann.Deep.LeanPrefix.B2

/-!
# Which three prefix reversals generate S_n below the diagonal: the full answer
(discoveries/prefix-reversals.md, Claim; NOTE.md Theorems A and B)

r_i = rev n i reverses the first i positions of Fin n (0-based: x ↦ i - 1 - x for x < i);
G n m k = ⟨r_n, r_m, r_k⟩ (Subgroup.closure); positions are x.val + 1 (1-based) in B.1.
Literature inputs (hypotheses): LZZW2026_Thm3_1 (arXiv 2610.30892, Thm 3.1) and
Jones2014_Thm1_2_k0 (Bull. Aust. Math. Soc. 89 (2014), Thm 1.2, k = 0); see Inputs.lean.
-/

namespace Riemann.Deep.LeanPrefix

open Equiv

theorem zmod_cast_iff (g a b : ℕ) : ((a : ℤ) ≡ (b : ℤ) [ZMOD g] ↔ a % g = b % g) := by
  unfold Int.ModEq
  rw [← Int.natCast_mod, ← Int.natCast_mod, Nat.cast_inj]

theorem zmod_neg_cast_iff (g a b : ℕ) : ((a : ℤ) ≡ -(b : ℤ) [ZMOD g] ↔ (a + b) % g = 0) := by
  rw [Int.modEq_iff_dvd, show -(b : ℤ) - a = -((a + b : ℕ) : ℤ) by push_cast; ring, dvd_neg,
    Int.natCast_dvd_natCast, Nat.dvd_iff_mod_eq_zero]

theorem mod4_iff (n k : ℕ) : (n ≡ 1 [MOD 4] ∧ 4 ∣ k) ↔ (n % 4 = 1 ∧ k % 4 = 0) := by
  unfold Nat.ModEq; rw [Nat.dvd_iff_mod_eq_zero]

/-- On the odd line with gcd = 1, (k+1) ∤ m and not (n ≡ 1 mod 4, 4 ∣ k): G = S_n. -/
theorem line_top (hJ : Jones2014_Thm1_2_k0) {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m)
    (hline : n = m + k + 1) (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1))
    (hdiv : ¬ (k + 1) ∣ m) (h4 : ¬ (n % 4 = 1 ∧ k % 4 = 0)) : G n m k = ⊤ := by
  have hA := alt_le_line hJ hk hkm hline hodd hg hdiv
  have hA' : ∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ G n m k :=
    fun σ h => hA (Equiv.Perm.mem_alternatingGroup.mpr h)
  have hev := even_of_coprime hline hodd hg
  have hnot : ¬ (∀ g ∈ G n m k, Perm.sign g = 1) := by
    rw [G_le_alt_iff (by omega) (by omega)]; omega
  obtain ⟨g, hgG, hsg⟩ : ∃ g ∈ G n m k, Perm.sign g ≠ 1 := by
    by_contra hc
    apply hnot
    intro g hg'
    by_contra hs
    exact hc ⟨g, hg', hs⟩
  exact eq_top_of_alt hA' g hgG ((Int.units_eq_one_or _).resolve_left hsg)

theorem alternatingGroup_ne_top {n : ℕ} (hn : 2 ≤ n) : alternatingGroup (Fin n) ≠ ⊤ := by
  intro h
  have : swap (⟨0, by omega⟩ : Fin n) ⟨1, by omega⟩ ∈ alternatingGroup (Fin n) := by rw [h]; trivial
  rw [Equiv.Perm.mem_alternatingGroup, Perm.sign_swap (by simp [Fin.ext_iff])] at this
  exact absurd this (by decide)

/-- **Main theorem** (discoveries/prefix-reversals.md, Claim), from the two literature inputs.

1. Theorem A: for 2 ≤ k < m < n with m + k < n, ⟨r_n, r_m, r_k⟩ = S_n iff n is odd,
   m + k = n - 1, gcd(m+1, k+1) = 1, (k+1) ∤ m and not (n ≡ 1 mod 4 and 4 ∣ k).
2. On the line m + k = n - 1 (2 ≤ k < m), with g = gcd(m+1, k+1), α = k + 1:
   B.1 the orbits are {y : y ≡ ±x (mod g)} (positions 1-based), there are g/2 + 1 of them,
   and G is transitive iff g = 1;
   B.2 n even, g = 1: transitive, imprimitive, with a block of n/2 points (two blocks);
   B.3 n odd, α ∣ n: imprimitive, the intervals of length α are blocks, G is the imprimitive
   wreath product D_{2α} ≀ D_{2n/α} in its action on ZMod (n/α) × ZMod α (up to a bijection with
   Fin n), G ≅ D_{2α} ≀ D_{2n/α} (Wreath α (n/α)) and |G| = (2α)^{n/α}·2(n/α);
   B.4 n odd, g = 1, α ∤ n: primitive, A_n ≤ G, G = A_n iff n ≡ 1 (mod 4) and 4 ∣ k,
   and G = S_n otherwise. -/
theorem prefix_reversals (hOff : LZZW2026_Thm3_1) (hJ : Jones2014_Thm1_2_k0) :
    (∀ n m k : ℕ, 2 ≤ k → k < m → m < n → m + k < n →
      (G n m k = ⊤ ↔ (Odd n ∧ m + k = n - 1 ∧ Nat.gcd (m + 1) (k + 1) = 1 ∧ ¬ (k + 1) ∣ m ∧
        ¬ (n ≡ 1 [MOD 4] ∧ 4 ∣ k)))) ∧
    (∀ n m k : ℕ, 2 ≤ k → k < m → m + k = n - 1 →
      (∀ x y : Fin n, y ∈ MulAction.orbit (G n m k) x ↔
          ((y : ℤ) + 1 ≡ (x : ℤ) + 1 [ZMOD (Nat.gcd (m + 1) (k + 1))] ∨
            (y : ℤ) + 1 ≡ -((x : ℤ) + 1) [ZMOD (Nat.gcd (m + 1) (k + 1))])) ∧
      Nat.card (MulAction.orbitRel.Quotient (G n m k) (Fin n)) = Nat.gcd (m + 1) (k + 1) / 2 + 1 ∧
      (MulAction.IsPretransitive (G n m k) (Fin n) ↔ Nat.gcd (m + 1) (k + 1) = 1) ∧
      (Even n → Nat.gcd (m + 1) (k + 1) = 1 →
        MulAction.IsPretransitive (G n m k) (Fin n) ∧ ¬ MulAction.IsPreprimitive (G n m k) (Fin n) ∧
        ∃ B : Set (Fin n), MulAction.IsBlock (G n m k) B ∧ 2 * B.ncard = n) ∧
      (Odd n → (k + 1) ∣ n →
        ¬ MulAction.IsPreprimitive (G n m k) (Fin n) ∧
        (∀ i : ℕ, MulAction.IsBlock (G n m k) {x : Fin n | x.val / (k + 1) = i}) ∧
        (∃ e : Fin n ≃ ZMod (n / (k + 1)) × ZMod (k + 1),
          (wreathPerm (k + 1) (n / (k + 1))).range.map e.symm.permCongrHom.toMonoidHom = G n m k) ∧
        Nonempty (G n m k ≃* Wreath (k + 1) (n / (k + 1))) ∧
        Nat.card (G n m k) = (2 * (k + 1)) ^ (n / (k + 1)) * (2 * (n / (k + 1)))) ∧
      (Odd n → Nat.gcd (m + 1) (k + 1) = 1 → ¬ (k + 1) ∣ n →
        MulAction.IsPreprimitive (G n m k) (Fin n) ∧ alternatingGroup (Fin n) ≤ G n m k ∧
        (G n m k = alternatingGroup (Fin n) ↔ (n ≡ 1 [MOD 4] ∧ 4 ∣ k)) ∧
        (¬ (n ≡ 1 [MOD 4] ∧ 4 ∣ k) → G n m k = ⊤))) := by
  refine ⟨fun n m k hk hkm hmn hsum => ?_, fun n m k hk hkm hl => ?_⟩
  · -- Theorem A
    rw [mod4_iff]
    constructor
    · intro htop
      by_cases hl : m + k < n - 1
      · exact absurd htop (hOff n m k hk hkm hmn hl)
      · have hline : n = m + k + 1 := by omega
        obtain ⟨h1, h2, h3, h4⟩ := line_necessary hk hkm hline htop
        exact ⟨Nat.odd_iff.mpr h1, by omega, h2, h3, h4⟩
    · rintro ⟨h1, h2, h3, h4, h5⟩
      exact line_top hJ hk hkm (by omega) (Nat.odd_iff.mp h1) h3 h4 h5
  · have hline : n = m + k + 1 := by omega
    have hgc : Nat.gcd (k + 1) (m + 1) = Nat.gcd (m + 1) (k + 1) := Nat.gcd_comm _ _
    have hdivm : ∀ d, d ∣ n ↔ d ∣ m + (k + 1) := fun d => by rw [show m + (k + 1) = n by omega]
    refine ⟨fun x y => ?_, ?_, ?_, ?_, ?_, ?_⟩
    · -- B.1, the orbits
      rw [mem_orbit_iff_mod hline hkm, hgc,
        show (y : ℤ) + 1 = ((y.val + 1 : ℕ) : ℤ) by push_cast; ring,
        show (x : ℤ) + 1 = ((x.val + 1 : ℕ) : ℤ) by push_cast; ring,
        zmod_cast_iff, zmod_neg_cast_iff]
    · -- B.1, the number of orbits
      rw [card_orbits hline hkm, hgc]
    · -- B.1, transitivity
      rw [isPretransitive_iff_gcd hline hkm, hgc]
    · -- B.2
      intro hev hg
      exact line_B2 hk hkm hline (Nat.even_iff.mp hev) (by rw [Nat.Coprime, hgc]; exact hg)
    · -- B.3
      intro hodd hdiv
      exact line_B3 hk hkm hline (Nat.odd_iff.mp hodd) hdiv
    · -- B.4
      intro hodd hg hdiv
      have hodd' := Nat.odd_iff.mp hodd
      have hdivm' : ¬ (k + 1) ∣ m := by
        intro h; apply hdiv; rw [hdivm]; exact Nat.dvd_add h dvd_rfl
      have hprim := isPreprimitive_of_line hk hkm hline hodd' hg hdivm'
      have hA := alt_le_line hJ hk hkm hline hodd' hg hdivm'
      have hAn : (n % 4 = 1 ∧ k % 4 = 0) → G n m k = alternatingGroup (Fin n) := by
        intro h4
        apply le_antisymm _ hA
        intro g hg'
        exact Equiv.Perm.mem_alternatingGroup.mpr (G_le_alt_of_mod4 hline hkm h4 g hg')
      have hSn : ¬ (n % 4 = 1 ∧ k % 4 = 0) → G n m k = ⊤ :=
        line_top hJ hk hkm hline hodd' hg hdivm'
      rw [mod4_iff]
      refine ⟨hprim, hA, ⟨fun h => ?_, hAn⟩, hSn⟩
      by_contra h4
      exact alternatingGroup_ne_top (n := n) (by omega) (h ▸ hSn h4)

end Riemann.Deep.LeanPrefix

