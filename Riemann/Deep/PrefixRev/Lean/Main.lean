import Riemann.Deep.PrefixRev.Lean.BlocksB

/-!
# Theorem A and Theorem B.4 from named hypotheses

Named hypotheses (Props, never axioms):
* OffLineLZZW: 2610.30892 (Lei-Zhao-Zhu-Wu), Theorem 3.1, first sentence (TeX l. 221-223):
  m + k < n - 1 ⇒ ⟨r_n, r_m, r_k⟩ < S_n (for 2 ≤ k < m < n, TeX l. 109-111).
* Prop5Primitive: the NOTE's Proposition 5 (n odd, line, gcd = 1, (k+1) ∤ m ⇒ G primitive).
  OPEN OBLIGATION: our own result, not formalized here.
* JonesLemma7: a primitive group of degree n containing an n-cycle and a nonidentity element
  with more than n/2 fixed points contains A_n. This is Jones, Bull. Aust. Math. Soc. 89 (2014),
  Thm 1.2, case k = 0 (arXiv 1209.5169, TeX l. 54-72), combined with the NOTE's Lemma 7
  (collineation fixed-point bound) and the fixed-point counts in AGL_1(p), L_2(11), M_11, M_23.
-/

namespace PrefixRev

open Equiv

/-- B is a block of H: every h ∈ H maps B onto B or onto a set disjoint from B. -/
def IsBlockOf {n : ℕ} (H : Subgroup (Perm (Fin n))) (B : Set (Fin n)) : Prop :=
  ∀ h ∈ H, (⇑h '' B = B) ∨ Disjoint (⇑h '' B) B

/-- H is primitive: transitive, and every block is trivial. -/
def IsPrimitive {n : ℕ} (H : Subgroup (Perm (Fin n))) : Prop :=
  (∀ x y : Fin n, ∃ h ∈ H, h x = y) ∧
    ∀ B : Set (Fin n), IsBlockOf H B → B.Subsingleton ∨ B = Set.univ

/-- Named hypothesis: 2610.30892, Theorem 3.1, first sentence. -/
def OffLineLZZW : Prop :=
  ∀ n m k : ℕ, 2 ≤ k → k < m → m < n → m + k < n - 1 → G n m k ≠ ⊤

/-- Named hypothesis (OPEN OBLIGATION, our Proposition 5): primitivity on the odd line. -/
def Prop5Primitive : Prop :=
  ∀ n m k : ℕ, 2 ≤ k → k < m → n = m + k + 1 → n % 2 = 1 → Nat.Coprime (m + 1) (k + 1) →
    ¬ (k + 1) ∣ m → IsPrimitive (G n m k)

/-- Named hypothesis: Jones 2014 Thm 1.2 (k = 0) + Lemma 7 (+ fixed points in the
sporadic and affine cases). -/
def JonesLemma7 : Prop :=
  ∀ (n : ℕ) (H : Subgroup (Perm (Fin n))), IsPrimitive H →
    (∃ c ∈ H, ∀ x y : Fin n, ∃ t : ℕ, (c ^ t) x = y) →
    (∃ h ∈ H, h ≠ 1 ∧ n < 2 * (Finset.univ.filter (fun x => h x = x)).card) →
    ∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ H

/-- r_k fixes the n - k points x ≥ k. -/
theorem fix_rk_card {n k : ℕ} (hk : k ≤ n) :
    n - k ≤ (Finset.univ.filter (fun x : Fin n => rev n k x = x)).card := by
  have h := Finset.card_le_card_of_injOn (s := (Finset.univ : Finset (Fin (n - k))))
    (t := Finset.univ.filter (fun x : Fin n => rev n k x = x))
    (fun i => (⟨k + i.val, by have := i.isLt; omega⟩ : Fin n))
    (by
      intro i _
      have hi := i.isLt
      rw [Finset.mem_coe, Finset.mem_filter]
      refine ⟨Finset.mem_univ _, ?_⟩
      apply Fin.ext
      rw [rev_val]
      dsimp only
      rw [ite_eq_right (by omega)])
    (by
      intro i _ j _ hij
      have := congrArg Fin.val hij
      exact Fin.ext (by simp only at this; omega))
  rwa [Finset.card_univ, Fintype.card_fin] at h

theorem rk_ne_one {n k : ℕ} (hk : 2 ≤ k) (hkn : k ≤ n) : rev n k ≠ 1 := by
  intro h
  have h0 : (rev n k ⟨0, by omega⟩).val = k - 1 := by
    rw [rev_val, ite_eq_left (show (0 : ℕ) < k ∧ k ≤ n from ⟨by omega, hkn⟩)]
    dsimp only
    omega
  rw [h, Perm.one_apply] at h0
  dsimp only at h0
  omega

theorem P_mem_G (n m k : ℕ) : P n m k ∈ G n m k :=
  (G n m k).mul_mem ((G n m k).mul_mem (rk_mem n m k) (rn_mem n m k)) (rm_mem n m k)

theorem even_of_coprime {n m k : ℕ} (hline : n = m + k + 1) (hodd : n % 2 = 1)
    (hg : Nat.Coprime (m + 1) (k + 1)) : m % 2 = 0 ∧ k % 2 = 0 := by
  by_contra h
  have h2 : 2 ∣ Nat.gcd (m + 1) (k + 1) := Nat.dvd_gcd (by omega) (by omega)
  rw [hg] at h2
  exact absurd (Nat.le_of_dvd (by omega) h2) (by omega)

/-- B.4, first part: on the odd line with gcd = 1 and (k+1) ∤ m, A_n ⊆ G. -/
theorem line_alt_le (hP5 : Prop5Primitive) (hJ : JonesLemma7) {n m k : ℕ} (hk : 2 ≤ k)
    (hkm : k < m) (hline : n = m + k + 1) (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1))
    (hdiv : ¬ (k + 1) ∣ m) : ∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ G n m k := by
  apply hJ n (G n m k) (hP5 n m k hk hkm hline hodd hg hdiv)
  · exact ⟨P n m k, P_mem_G n m k, P_cycle hline hkm (Nat.Coprime.symm hg)⟩
  · refine ⟨rev n k, rk_mem n m k, rk_ne_one hk (by omega), ?_⟩
    have := fix_rk_card (n := n) (k := k) (by omega)
    omega

/-- B.4: on the odd line with gcd = 1 and (k+1) ∤ m: G = A_n if n ≡ 1 (mod 4) and 4 | k,
and G = S_n otherwise. -/
theorem line_B4 (hP5 : Prop5Primitive) (hJ : JonesLemma7) {n m k : ℕ} (hk : 2 ≤ k)
    (hkm : k < m) (hline : n = m + k + 1) (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1))
    (hdiv : ¬ (k + 1) ∣ m) :
    ((n % 4 = 1 ∧ k % 4 = 0) → ∀ σ : Perm (Fin n), σ ∈ G n m k ↔ Perm.sign σ = 1) ∧
    (¬ (n % 4 = 1 ∧ k % 4 = 0) → G n m k = ⊤) := by
  have hA := line_alt_le hP5 hJ hk hkm hline hodd hg hdiv
  have hev := even_of_coprime hline hodd hg
  constructor
  · intro h4 σ
    exact ⟨fun hσ => G_le_alt_of_mod4 hline hkm h4 σ hσ, hA σ⟩
  · intro h4
    have hnot : ¬ (∀ g ∈ G n m k, Perm.sign g = 1) := by
      rw [G_le_alt_iff (by omega) (by omega)]; omega
    obtain ⟨g, hgG, hsg⟩ : ∃ g ∈ G n m k, Perm.sign g ≠ 1 := by
      by_contra hc
      apply hnot
      intro g hg
      by_contra hs
      exact hc ⟨g, hg, hs⟩
    have hsg' : Perm.sign g = -1 := (Int.units_eq_one_or _).resolve_left hsg
    exact eq_top_of_alt hA g hgG hsg'

/-- Theorem A (from the named hypotheses): for 2 ≤ k < m < n with m + k < n,
⟨r_n, r_m, r_k⟩ = S_n iff n is odd, m + k = n - 1, gcd(m+1, k+1) = 1, (k+1) ∤ m and
not (n ≡ 1 mod 4 and 4 | k). -/
theorem theoremA (hOff : OffLineLZZW) (hP5 : Prop5Primitive) (hJ : JonesLemma7) {n m k : ℕ}
    (hk : 2 ≤ k) (hkm : k < m) (hmn : m < n) (hsum : m + k < n) :
    G n m k = ⊤ ↔ (n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧
      ¬ (n % 4 = 1 ∧ k % 4 = 0)) := by
  constructor
  · intro htop
    by_cases hl : m + k < n - 1
    · exact absurd htop (hOff n m k hk hkm hmn hl)
    · have hline : n = m + k + 1 := by omega
      obtain ⟨h1, h2, h3, h4⟩ := line_necessary hk hkm hline htop
      exact ⟨h1, by omega, h2, h3, h4⟩
  · rintro ⟨h1, h2, h3, h4, h5⟩
    exact (line_B4 hP5 hJ hk hkm (by omega) h1 h3 h4).2 h5

/-- Theorem A, the direction that needs no Prop 5 and no Jones: G = S_n ⇒ the conditions
(only the published off-line result is assumed). -/
theorem theoremA_necessary (hOff : OffLineLZZW) {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m)
    (hmn : m < n) (hsum : m + k < n) (htop : G n m k = ⊤) :
    n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧
      ¬ (n % 4 = 1 ∧ k % 4 = 0) := by
  by_cases hl : m + k < n - 1
  · exact absurd htop (hOff n m k hk hkm hmn hl)
  · have hline : n = m + k + 1 := by omega
    obtain ⟨h1, h2, h3, h4⟩ := line_necessary hk hkm hline htop
    exact ⟨h1, by omega, h2, h3, h4⟩

end PrefixRev
