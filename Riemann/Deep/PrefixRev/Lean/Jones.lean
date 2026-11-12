import Riemann.Deep.PrefixRev.Lean.PrimB2
import Mathlib.Data.Nat.Prime.Defs

/-!
# Jones's theorem (k = 0) and Lemma 7 as two separate named hypotheses

Jones, Primitive permutation groups containing a cycle, Bull. Aust. Math. Soc. 89 (2014),
Theorem 1.2 (arXiv 1209.5169, TeX l. 54-72): G primitive of degree n, A_n ⊄ G, G contains a
cycle fixing k = 0 points (an n-cycle) ⇒ (a) C_p ≤ G ≤ AGL_1(p), n = p prime; or
(b) PGL_d(q) ≤ G ≤ PΓL_d(q), n = (q^d - 1)/(q - 1), d ≥ 2, q a prime power; or
(c) G = L_2(11), M_11 or M_23 with n = 11, 11 or 23.

JonesK0 Coll states this with: (a) weakened to "nonidentity elements fix ≤ 1 point" (true in
AGL_1(p)); (b) through a predicate Coll q d H whose intended meaning is "PGL_d(q) ≤ H ≤ PΓL_d(q)
acting on the points of PG(d-1, q)"; (c) weakened to the fixed-point bounds 3 (L_2(11), M_11 on
11 points: M_11 is sharply 4-transitive) and 7 (M_23 on 23 points). Lemma7 Coll is the NOTE's
Lemma 7 for the same predicate. Both are Props (hypotheses), never axioms.
-/

namespace PrefixRev

open Equiv

/-- Number of fixed points. -/
def nfix {n : ℕ} (h : Perm (Fin n)) : ℕ := (Finset.univ.filter (fun x => h x = x)).card

/-- Named hypothesis: Jones 2014, Theorem 1.2, case k = 0 (see the module doc). -/
def JonesK0 (Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop) : Prop :=
  ∀ (n : ℕ) (H : Subgroup (Perm (Fin n))), IsPrimitive H →
    ¬ (∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ H) →
    (∃ c ∈ H, ∀ x y : Fin n, ∃ t : ℕ, (c ^ t) x = y) →
    (n.Prime ∧ ∀ h ∈ H, h ≠ 1 → nfix h ≤ 1) ∨
    (∃ q d : ℕ, (∃ p e : ℕ, p.Prime ∧ 0 < e ∧ q = p ^ e) ∧ 2 ≤ d ∧ n * (q - 1) = q ^ d - 1 ∧
      Coll q d H) ∨
    (n = 11 ∧ ∀ h ∈ H, h ≠ 1 → nfix h ≤ 3) ∨ (n = 23 ∧ ∀ h ∈ H, h ≠ 1 → nfix h ≤ 7)

/-- Named hypothesis: the NOTE's Lemma 7 (a nonidentity collineation of PG(d-1, q),
(d, q) ≠ (2, 4), fixes at most n/2 points). -/
def Lemma7 (Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop) : Prop :=
  ∀ (n q d : ℕ) (H : Subgroup (Perm (Fin n))), Coll q d H → 2 ≤ d → n * (q - 1) = q ^ d - 1 →
    ¬ (d = 2 ∧ q = 4) → ∀ h ∈ H, h ≠ 1 → 2 * nfix h ≤ n

/-- B.4, first part, from JonesK0 and Lemma7: A_n ⊆ G on the odd line, gcd = 1, (k+1) ∤ m. -/
theorem line_alt_le_jones {Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop}
    (hP5 : Prop5Primitive) (hJ : JonesK0 Coll) (hL : Lemma7 Coll) {n m k : ℕ} (hk : 2 ≤ k)
    (hkm : k < m) (hline : n = m + k + 1) (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1))
    (hdiv : ¬ (k + 1) ∣ m) : ∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ G n m k := by
  by_contra hA
  have hfix : n - k ≤ nfix (rev n k) := fix_rk_card (by omega)
  have hne := rk_ne_one hk (show k ≤ n by omega)
  have hrk := rk_mem n m k
  have hev := even_of_coprime hline hodd hg
  rcases hJ n (G n m k) (hP5 n m k hk hkm hline hodd hg hdiv) hA
      ⟨P n m k, P_mem_G n m k, P_cycle hline hkm (Nat.Coprime.symm hg)⟩ with
    ⟨_, h⟩ | ⟨q, d, _, hd, hn, hC⟩ | ⟨h11, h⟩ | ⟨h23, h⟩
  · have := h _ hrk hne; omega
  · have h24 : ¬ (d = 2 ∧ q = 4) := by
      rintro ⟨rfl, rfl⟩
      have : n * 3 = 15 := by simpa using hn
      omega
    have := hL n q d (G n m k) hC hd hn h24 _ hrk hne
    omega
  · have := h _ hrk hne; omega
  · have := h _ hrk hne; omega

/-- Theorem A from: the published off-line result, the arithmetic core of Prop 5 (open),
Jones's theorem (k = 0) and Lemma 7. -/
theorem theoremA_jones {Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop}
    (hOff : OffLineLZZW) (hC : Prop5Core) (hJ : JonesK0 Coll) (hL : Lemma7 Coll) {n m k : ℕ}
    (hk : 2 ≤ k) (hkm : k < m) (hmn : m < n) (hsum : m + k < n) :
    G n m k = ⊤ ↔ (n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧
      ¬ (n % 4 = 1 ∧ k % 4 = 0)) := by
  constructor
  · exact theoremA_necessary hOff hk hkm hmn hsum
  · rintro ⟨h1, h2, h3, h4, h5⟩
    have hline : n = m + k + 1 := by omega
    have hA := line_alt_le_jones (prop5_of_core hC) hJ hL hk hkm hline h1 h3 h4
    have hev := even_of_coprime hline h1 h3
    have hnot : ¬ (∀ g ∈ G n m k, Perm.sign g = 1) := by
      rw [G_le_alt_iff (by omega) (by omega)]; omega
    obtain ⟨g, hgG, hsg⟩ : ∃ g ∈ G n m k, Perm.sign g ≠ 1 := by
      by_contra hc
      apply hnot
      intro g hg
      by_contra hs
      exact hc ⟨g, hg, hs⟩
    exact eq_top_of_alt hA g hgG ((Int.units_eq_one_or _).resolve_left hsg)

/-- B.4, the alternating branch, from the same hypotheses: if n ≡ 1 (mod 4) and 4 | k then
G = A_n (the S_n branch is theoremA_jones). -/
theorem line_B4_jones {Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop}
    (hC : Prop5Core) (hJ : JonesK0 Coll) (hL : Lemma7 Coll) {n m k : ℕ} (hk : 2 ≤ k)
    (hkm : k < m) (hline : n = m + k + 1) (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1))
    (hdiv : ¬ (k + 1) ∣ m) (h4 : n % 4 = 1 ∧ k % 4 = 0) :
    ∀ σ : Perm (Fin n), σ ∈ G n m k ↔ Perm.sign σ = 1 := fun σ =>
  ⟨fun hσ => G_le_alt_of_mod4 hline hkm h4 σ hσ,
    line_alt_le_jones (prop5_of_core hC) hJ hL hk hkm hline hodd hg hdiv σ⟩

end PrefixRev
