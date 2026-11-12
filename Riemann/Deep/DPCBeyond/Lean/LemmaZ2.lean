import Riemann.Deep.DPCBeyond.Lean.LemmaZ

/-!
# dpcbeyond Lemma Z, case 2 and the pointwise theorem from the named inputs (i)–(iii)

term_bound2, term_nofree1, term_nofree2: (iii)+(iv) per term (case 2: packets 𝔞^{−1/2+y}, free part
(X/𝔞)^{y + max d⁺}, i.e. (iii) with D = y + max d⁺ ≥ max (d_i + y)⁺, prod_ex_le; no free variable: 𝔞 ≤ X).
ZPoint1, ZPoint2: Lemma Z's NAMED inputs (i)–(iii) at one point (the decomposition into ≤ NT partition/binomial
terms, the packet sums, the free-sum bounds, the O(T^{−10}) remainder E ≤ 1). lemmaZ_point1, lemmaZ_point2:
**Lemma Z** from them: |R_T(t − iy)| ≤ K(1 + W_{−y})(1 + On)^{k−1}, |R*(t − iy)| ≤ K X^y (1 + W)(1 + On)^{k−1},
K = NT·K_pk·C_I·k + 1, W = W_{−0}. W_split: W_{−y} ≤ On + X^{−y}Off (6e′ (a)).
-/

namespace DPCBLean

open Finset DPCLean

variable {ex lg : ℝ → ℝ}

/-- W_{−y} = Σ_ρ X^{(d_ρ − y)⁺} κ(ρ) (W = W_{−0}). -/
noncomputable def Wsum (ex lg : ℝ → ℝ) {ι : Type*} (X y : ℝ) (S : Finset ι) (d κ : ι → ℝ) : ℝ :=
  ∑ x ∈ S, rp ex lg X (max (d x - y) 0) * κ x

/-- On = Σ_ρ κ(ρ). -/
noncomputable def Onsum {ι : Type*} (S : Finset ι) (κ : ι → ℝ) : ℝ := ∑ x ∈ S, κ x

/-- Off = Σ_{d_ρ > 0} X^{d_ρ} κ(ρ). -/
noncomputable def Offsum (ex lg : ℝ → ℝ) {ι : Type*} (X : ℝ) (S : Finset ι) (d κ : ι → ℝ) : ℝ :=
  ∑ x ∈ S.filter (fun x => 0 < d x), rp ex lg X (d x) * κ x

/-- **6e′ (a)**, summed: W_{−y} ≤ On + X^{−y}Off (X ≥ 1, y ≥ 0, κ ≥ 0). -/
theorem W_split (h : ExpLog ex lg) {ι : Type*} (S : Finset ι) (d κ : ι → ℝ) (hκ : ∀ x, 0 ≤ κ x)
    {X y : ℝ} (hX : 1 ≤ X) (hy : 0 ≤ y) :
    Wsum ex lg X y S d κ ≤ Onsum S κ + rp ex lg X (-y) * Offsum ex lg X S d κ := by
  unfold Wsum Onsum Offsum
  rw [mul_sum, sum_filter, ← sum_add_distrib]
  apply sum_le_sum
  intro x _
  have h1 := mul_le_mul_of_nonneg_right (pos_part_split h (d := d x) hX hy) (hκ x)
  split_ifs at h1 ⊢ with hd
  · have e : (1 + rp ex lg X (-y) * rp ex lg X (d x)) * κ x =
        κ x + rp ex lg X (-y) * (rp ex lg X (d x) * κ x) := by ring
    linarith only [h1, e]
  · have e : (1 + 0) * κ x = κ x := by ring
    rw [add_zero]
    linarith only [h1, e]

/-- (iii)+(iv), one tuple in case 2: 𝔞^{−1/2+y}((X/𝔞)^{y+M}Q) ≤ 𝔞^{−1/2}(X^y(X^M Q)). -/
lemma tuple_case2 (h : ExpLog ex lg) {a X y M Q : ℝ} (ha : 1 ≤ a) (hX : 0 < X) (hM : 0 ≤ M)
    (hQ : 0 ≤ Q) : rp ex lg a (-(1 / 2) + y) * (rp ex lg (X / a) (y + M) * Q) ≤
      rp ex lg a (-(1 / 2)) * (rp ex lg X y * (rp ex lg X M * Q)) := by
  have h1 := packet_case2 h (y := y) ha hX hM
  rw [rp_add h hX] at h1
  have e1 : rp ex lg a (-(1 / 2) + y) * (rp ex lg (X / a) (y + M) * Q) =
      (rp ex lg a (-(1 / 2) + y) * rp ex lg (X / a) (y + M)) * Q := by ring
  have e2 : rp ex lg a (-(1 / 2)) * (rp ex lg X y * (rp ex lg X M * Q)) =
      (rp ex lg X y * rp ex lg X M * rp ex lg a (-(1 / 2))) * Q := by ring
  rw [e1, e2]
  exact mul_le_mul_of_nonneg_right h1 hQ

/-- **Lemma Z (iii)+(iv) for one term, case 2** (s = 1/2 − y + it, R*(t − iy)): NAMED |F 𝔞| ≤
C_I Σ_{ρ⃗∈S^m}(X/𝔞)^{y + max_i c(ρ_i)}∏k(ρ_i) (c = d⁺ ≥ 0) and Σ_𝔞 C 𝔞 v^{−1/2} ≤ K_pk give
Σ_𝔞 C 𝔞 v^{−1/2+y}|F 𝔞| ≤ K_pk C_I X^y m W On^{m−1}. -/
theorem term_bound2 (h : ExpLog ex lg) {ι α : Type*} [DecidableEq ι] (S : Finset ι) (c k : ι → ℝ)
    (hc : ∀ x, 0 ≤ c x) (hk : ∀ x, 0 ≤ k x) {X y CI Kpk : ℝ} (hX : 0 < X) (hCI : 0 ≤ CI) (n : ℕ)
    (P : Finset α) (v C F : α → ℝ) (hv : ∀ p ∈ P, 1 ≤ v p) (hC : ∀ p ∈ P, 0 ≤ C p)
    (hF : ∀ p ∈ P, |F p| ≤ CI * ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
      rp ex lg (X / v p) (y + univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i))
    (hP : ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) ≤ Kpk) :
    ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) + y) * |F p| ≤
      Kpk * CI * (rp ex lg X y * ((n + 1) * (∑ x ∈ S, rp ex lg X (c x) * k x) * (∑ x ∈ S, k x) ^ n)) := by
  obtain ⟨G, hG⟩ : ∃ G : ℝ, G = (n + 1) * (∑ x ∈ S, rp ex lg X (c x) * k x) * (∑ x ∈ S, k x) ^ n :=
    ⟨_, rfl⟩
  rw [← hG]
  have hG0 : 0 ≤ G := by
    rw [hG]
    exact mul_nonneg (mul_nonneg (by positivity) (sum_nonneg (fun x _ =>
      mul_nonneg (rp_nonneg h _ _) (hk x)))) (pow_nonneg (sum_nonneg (fun x _ => hk x)) _)
  have hXy : 0 ≤ rp ex lg X y := rp_nonneg h _ _
  have hterm : ∀ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) + y) * |F p| ≤
      C p * rp ex lg (v p) (-(1 / 2)) * (CI * (rp ex lg X y * G)) := by
    intro p hp
    have hvp := hv p hp
    have hsum : rp ex lg (v p) (-(1 / 2) + y) * ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
        rp ex lg (X / v p) (y + univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i) ≤
        rp ex lg (v p) (-(1 / 2)) * (rp ex lg X y * G) := by
      rw [mul_sum]
      have h2 : ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S), rp ex lg (v p) (-(1 / 2)) *
          (rp ex lg X y * (rp ex lg X (univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i))) =
          rp ex lg (v p) (-(1 / 2)) * (rp ex lg X y * ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (n + 1) => S),
            rp ex lg X (univ.sup' univ_nonempty (fun i => c (ρ i))) * ∏ i, k (ρ i)) := by
        rw [mul_sum, mul_sum]
      have h3 := tuple_bound h S c k hk X n
      rw [← hG] at h3
      refine le_trans (sum_le_sum (fun ρ _ => ?_)) (le_of_eq_of_le h2
        (mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left h3 hXy) (rp_nonneg h _ _)))
      have hs0 : 0 ≤ univ.sup' univ_nonempty (fun i => c (ρ i)) := by
        obtain ⟨i, -, he⟩ := exists_mem_eq_sup' univ_nonempty (fun i => c (ρ i))
        rw [he]
        exact hc _
      exact tuple_case2 h hvp hX hs0 (prod_nonneg (fun i _ => hk _))
    have h4 := mul_le_mul_of_nonneg_left (hF p hp) (rp_nonneg h (v p) (-(1 / 2) + y))
    have h5 : rp ex lg (v p) (-(1 / 2) + y) * |F p| ≤ CI * (rp ex lg (v p) (-(1 / 2)) * (rp ex lg X y * G)) := by
      refine le_trans h4 (le_of_eq_of_le (by ring) (mul_le_mul_of_nonneg_left hsum hCI))
    calc C p * rp ex lg (v p) (-(1 / 2) + y) * |F p| = C p * (rp ex lg (v p) (-(1 / 2) + y) * |F p|) := by
          ring
      _ ≤ C p * (CI * (rp ex lg (v p) (-(1 / 2)) * (rp ex lg X y * G))) :=
          mul_le_mul_of_nonneg_left h5 (hC p hp)
      _ = C p * rp ex lg (v p) (-(1 / 2)) * (CI * (rp ex lg X y * G)) := by ring
  calc ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) + y) * |F p| ≤
        ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) * (CI * (rp ex lg X y * G)) := sum_le_sum hterm
    _ = (∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2))) * (CI * (rp ex lg X y * G)) := by rw [sum_mul]
    _ ≤ Kpk * (CI * (rp ex lg X y * G)) :=
        mul_le_mul_of_nonneg_right hP (mul_nonneg hCI (mul_nonneg hXy hG0))
    _ = Kpk * CI * (rp ex lg X y * G) := by ring

/-- A term with no free variable, case 1: Σ_𝔞 C v^{−1/2−y}|F| ≤ K_pk C_I when |F| ≤ C_I (y ≥ 0, v ≥ 1). -/
theorem term_nofree1 (h : ExpLog ex lg) {α : Type*} {y CI Kpk : ℝ} (hy : 0 ≤ y) (hCI : 0 ≤ CI)
    (P : Finset α) (v C F : α → ℝ) (hv : ∀ p ∈ P, 1 ≤ v p) (hC : ∀ p ∈ P, 0 ≤ C p)
    (hF : ∀ p ∈ P, |F p| ≤ CI) (hP : ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) ≤ Kpk) :
    ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) - y) * |F p| ≤ Kpk * CI := by
  calc ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) - y) * |F p| ≤
        ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) * CI := by
        apply sum_le_sum
        intro p hp
        exact mul_le_mul (mul_le_mul_of_nonneg_left (rp_le_rp_exp h (hv p hp) (by linarith only [hy]))
          (hC p hp)) (hF p hp) (abs_nonneg _) (mul_nonneg (hC p hp) (rp_nonneg h _ _))
    _ = (∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2))) * CI := by rw [sum_mul]
    _ ≤ Kpk * CI := mul_le_mul_of_nonneg_right hP hCI

/-- A term with no free variable, case 2: 𝔞 ≤ X gives Σ_𝔞 C v^{−1/2+y}|F| ≤ K_pk C_I X^y. -/
theorem term_nofree2 (h : ExpLog ex lg) {α : Type*} {X y CI Kpk : ℝ} (hy : 0 ≤ y) (hCI : 0 ≤ CI)
    (P : Finset α) (v C F : α → ℝ) (hv : ∀ p ∈ P, 1 ≤ v p) (hvX : ∀ p ∈ P, v p ≤ X)
    (hC : ∀ p ∈ P, 0 ≤ C p) (hF : ∀ p ∈ P, |F p| ≤ CI)
    (hP : ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) ≤ Kpk) :
    ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) + y) * |F p| ≤ Kpk * (CI * rp ex lg X y) := by
  calc ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2) + y) * |F p| ≤
        ∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2)) * (CI * rp ex lg X y) := by
        apply sum_le_sum
        intro p hp
        have h1 := packet_nofree h (hv p hp) (hvX p hp) hy
        have h2 := mul_le_mul h1 (hF p hp) (abs_nonneg _)
          (mul_nonneg (rp_nonneg h _ _) (rp_nonneg h _ _))
        have h3 := mul_le_mul_of_nonneg_left h2 (hC p hp)
        have e1 : C p * rp ex lg (v p) (-(1 / 2) + y) * |F p| =
            C p * (rp ex lg (v p) (-(1 / 2) + y) * |F p|) := by ring
        have e2 : C p * (rp ex lg X y * rp ex lg (v p) (-(1 / 2)) * CI) =
            C p * rp ex lg (v p) (-(1 / 2)) * (CI * rp ex lg X y) := by ring
        linarith only [h3, e1, e2]
    _ = (∑ p ∈ P, C p * rp ex lg (v p) (-(1 / 2))) * (CI * rp ex lg X y) := by rw [sum_mul]
    _ ≤ Kpk * (CI * rp ex lg X y) :=
        mul_le_mul_of_nonneg_right hP (mul_nonneg hCI (rp_nonneg h _ _))

/-- Lemma Z's NAMED inputs (i)–(iii) at one point, case 1 (Rabs = |R_T(t − iy)|; zero set S with depths d and kernel
weights κ = κ_A(L|γ − t|)): (i)–(ii) |R| ≤ Σ_{τ∈𝒯}Σ_{𝔞∈P_τ} C 𝔞^{−1/2−y}|F| + E with #𝒯 ≤ NT and E ≤ 1 (pole and
trivial-zero terms); packets 1 ≤ v, 0 ≤ C, Σ C v^{−1/2} ≤ K_pk (Mertens); each term has n_τ + 1 ≤ k free variables with
(iii) |F| ≤ C_I Σ_{ρ⃗∈S^{n_τ+1}}(X/𝔞)^{max_i (d_i − y)⁺}∏κ, or no free variable and |F| ≤ C_I. -/
def ZPoint1 (ex lg : ℝ → ℝ) {ι : Type*} [DecidableEq ι] (X y CI Kpk : ℝ) (k NT : ℕ) (S : Finset ι)
    (d κ : ι → ℝ) (Rabs : ℝ) : Prop :=
  ∃ (𝒯 : Finset ℕ) (nt : ℕ → ℕ) (P : ℕ → Finset ℕ) (v C F : ℕ → ℕ → ℝ) (E : ℝ), 𝒯.card ≤ NT ∧ E ≤ 1 ∧
    Rabs ≤ ∑ τ ∈ 𝒯, ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) - y) * |F τ p| + E ∧
    ∀ τ ∈ 𝒯, nt τ + 1 ≤ k ∧ (∀ p ∈ P τ, 1 ≤ v τ p ∧ 0 ≤ C τ p) ∧
      ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2)) ≤ Kpk ∧
      ((∀ p ∈ P τ, |F τ p| ≤ CI * ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (nt τ + 1) => S),
          rp ex lg (X / v τ p) (univ.sup' univ_nonempty (fun i => max (d (ρ i) - y) 0)) *
            ∏ i, κ (ρ i)) ∨
        ∀ p ∈ P τ, |F τ p| ≤ CI)

/-- Case 2 (Rabs = |R*(t − iy)|, packets 𝔞^{−1/2+y}): as ZPoint1 with the free bound (X/𝔞)^{y + max_i d_i⁺} and, for
terms without free variable, packets 𝔞 ≤ X. -/
def ZPoint2 (ex lg : ℝ → ℝ) {ι : Type*} [DecidableEq ι] (X y CI Kpk : ℝ) (k NT : ℕ) (S : Finset ι)
    (d κ : ι → ℝ) (Rabs : ℝ) : Prop :=
  ∃ (𝒯 : Finset ℕ) (nt : ℕ → ℕ) (P : ℕ → Finset ℕ) (v C F : ℕ → ℕ → ℝ) (E : ℝ), 𝒯.card ≤ NT ∧ E ≤ 1 ∧
    Rabs ≤ ∑ τ ∈ 𝒯, ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) + y) * |F τ p| + E ∧
    ∀ τ ∈ 𝒯, nt τ + 1 ≤ k ∧ (∀ p ∈ P τ, 1 ≤ v τ p ∧ 0 ≤ C τ p) ∧
      ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2)) ≤ Kpk ∧
      ((∀ p ∈ P τ, |F τ p| ≤ CI * ∑ ρ ∈ Fintype.piFinset (fun _ : Fin (nt τ + 1) => S),
          rp ex lg (X / v τ p) (y + univ.sup' univ_nonempty (fun i => max (d (ρ i) - 0) 0)) *
            ∏ i, κ (ρ i)) ∨
        ∀ p ∈ P τ, v τ p ≤ X ∧ |F τ p| ≤ CI)

lemma final_shape {c NT KC Q E Rabs : ℝ} (hc : c ≤ NT) (hKC : 0 ≤ KC) (hQ : 1 ≤ Q) (hE : E ≤ 1)
    (hR : Rabs ≤ c * (KC * Q) + E) : Rabs ≤ (NT * KC + 1) * Q := by
  have h2 : c * KC * Q ≤ NT * KC * Q :=
    mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right hc hKC) (by linarith only [hQ])
  have e1 : c * (KC * Q) = c * KC * Q := by ring
  have e2 : (NT * KC + 1) * Q = NT * KC * Q + Q := by ring
  linarith only [hR, h2, e1, e2, hE, hQ]

/-- **Lemma Z, case 1, from its named inputs**: |R_T(t − iy)| ≤ (NT K_pk C_I k + 1)(1 + W_{−y})(1 + On)^{k−1}. -/
theorem lemmaZ_point1 (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] {X y CI Kpk : ℝ} {k NT : ℕ}
    {S : Finset ι} {d κ : ι → ℝ} {Rabs : ℝ} (hX : 0 < X) (hy : 0 ≤ y) (hCI : 0 ≤ CI) (hKpk : 0 ≤ Kpk)
    (hk : 1 ≤ k) (hκ : ∀ x, 0 ≤ κ x) (hz : ZPoint1 ex lg X y CI Kpk k NT S d κ Rabs) :
    Rabs ≤ (NT * (Kpk * CI * k) + 1) * ((1 + Wsum ex lg X y S d κ) * (1 + Onsum S κ) ^ (k - 1)) := by
  obtain ⟨𝒯, nt, P, v, C, F, E, hcard, hE, hR, hT⟩ := hz
  have hW : 0 ≤ Wsum ex lg X y S d κ := sum_nonneg (fun x _ => mul_nonneg (rp_nonneg h _ _) (hκ x))
  have hOn : 0 ≤ Onsum S κ := sum_nonneg (fun x _ => hκ x)
  have hB : ∀ τ ∈ 𝒯, nt τ + 1 ≤ k ∧
      (∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) - y) * |F τ p| ≤
        Kpk * CI * (((nt τ : ℝ) + 1) * Wsum ex lg X y S d κ * Onsum S κ ^ nt τ) ∨
       ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) - y) * |F τ p| ≤ Kpk * CI) := by
    intro τ hτ
    obtain ⟨hn, hvC, hP, hfree | hnone⟩ := hT τ hτ
    · exact ⟨hn, Or.inl (term_bound h S (fun x => max (d x - y) 0) κ (fun x => le_max_right _ _) hκ
        hX hy hCI (nt τ) (P τ) (v τ) (C τ) (F τ) (fun p hp => (hvC p hp).1)
        (fun p hp => (hvC p hp).2) hfree hP)⟩
    · exact ⟨hn, Or.inr (term_nofree1 h hy hCI (P τ) (v τ) (C τ) (F τ) (fun p hp => (hvC p hp).1)
        (fun p hp => (hvC p hp).2) hnone hP)⟩
  have hZ := lemmaZ_sum 𝒯 (fun τ => ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) - y) * |F τ p|) nt
    hk hW hOn (mul_nonneg hKpk hCI) hR hB
  have hQ : 1 ≤ (1 + Wsum ex lg X y S d κ) * (1 + Onsum S κ) ^ (k - 1) :=
    one_le_mul_of_one_le_of_one_le (by linarith only [hW]) (one_le_pow₀ (by linarith only [hOn]))
  have hc : (𝒯.card : ℝ) ≤ NT := by exact_mod_cast hcard
  have e : (𝒯.card : ℝ) * (Kpk * CI * (k * (1 + Wsum ex lg X y S d κ) * (1 + Onsum S κ) ^ (k - 1))) + E =
      (𝒯.card : ℝ) * (Kpk * CI * k * ((1 + Wsum ex lg X y S d κ) * (1 + Onsum S κ) ^ (k - 1))) + E := by
    ring
  rw [e] at hZ
  exact final_shape hc (mul_nonneg (mul_nonneg hKpk hCI) (Nat.cast_nonneg k)) hQ hE hZ

/-- **Lemma Z, case 2, from its named inputs**: |R*(t − iy)| ≤ (NT K_pk C_I k + 1) X^y (1 + W)(1 + On)^{k−1}. -/
theorem lemmaZ_point2 (h : ExpLog ex lg) {ι : Type*} [DecidableEq ι] {X y CI Kpk : ℝ} {k NT : ℕ}
    {S : Finset ι} {d κ : ι → ℝ} {Rabs : ℝ} (hX : 1 ≤ X) (hy : 0 ≤ y) (hCI : 0 ≤ CI) (hKpk : 0 ≤ Kpk)
    (hk : 1 ≤ k) (hκ : ∀ x, 0 ≤ κ x) (hz : ZPoint2 ex lg X y CI Kpk k NT S d κ Rabs) :
    Rabs ≤ (NT * (Kpk * CI * k) + 1) *
      (rp ex lg X y * ((1 + Wsum ex lg X 0 S d κ) * (1 + Onsum S κ) ^ (k - 1))) := by
  obtain ⟨𝒯, nt, P, v, C, F, E, hcard, hE, hR, hT⟩ := hz
  have hX0 : 0 < X := by linarith only [hX]
  have hXy : 1 ≤ rp ex lg X y := by
    have h1 := rp_le_rp_exp h hX hy
    rwa [rp_zero h hX0] at h1
  have hW : 0 ≤ Wsum ex lg X 0 S d κ := sum_nonneg (fun x _ => mul_nonneg (rp_nonneg h _ _) (hκ x))
  have hOn : 0 ≤ Onsum S κ := sum_nonneg (fun x _ => hκ x)
  have hCI' : 0 ≤ CI * rp ex lg X y := mul_nonneg hCI (rp_nonneg h _ _)
  have hB : ∀ τ ∈ 𝒯, nt τ + 1 ≤ k ∧
      (∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) + y) * |F τ p| ≤
        Kpk * (CI * rp ex lg X y) * (((nt τ : ℝ) + 1) * Wsum ex lg X 0 S d κ * Onsum S κ ^ nt τ) ∨
       ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) + y) * |F τ p| ≤ Kpk * (CI * rp ex lg X y)) := by
    intro τ hτ
    obtain ⟨hn, hvC, hP, hfree | hnone⟩ := hT τ hτ
    · refine ⟨hn, Or.inl ?_⟩
      have h1 := term_bound2 h S (fun x => max (d x - 0) 0) κ (fun x => le_max_right _ _) hκ hX0 hCI
        (nt τ) (P τ) (v τ) (C τ) (F τ) (fun p hp => (hvC p hp).1) (fun p hp => (hvC p hp).2) hfree hP
      refine le_of_le_of_eq h1 ?_
      unfold Wsum Onsum
      ring
    · exact ⟨hn, Or.inr (term_nofree2 h hy hCI (P τ) (v τ) (C τ) (F τ) (fun p hp => (hvC p hp).1)
        (fun p hp => (hnone p hp).1) (fun p hp => (hvC p hp).2) (fun p hp => (hnone p hp).2) hP)⟩
  have hZ := lemmaZ_sum 𝒯 (fun τ => ∑ p ∈ P τ, C τ p * rp ex lg (v τ p) (-(1 / 2) + y) * |F τ p|) nt
    hk hW hOn (mul_nonneg hKpk hCI') hR hB
  have hQ1 : 1 ≤ (1 + Wsum ex lg X 0 S d κ) * (1 + Onsum S κ) ^ (k - 1) :=
    one_le_mul_of_one_le_of_one_le (by linarith only [hW]) (one_le_pow₀ (by linarith only [hOn]))
  have hQ : 1 ≤ rp ex lg X y * ((1 + Wsum ex lg X 0 S d κ) * (1 + Onsum S κ) ^ (k - 1)) :=
    one_le_mul_of_one_le_of_one_le hXy hQ1
  have hc : (𝒯.card : ℝ) ≤ NT := by exact_mod_cast hcard
  have e : (𝒯.card : ℝ) * (Kpk * (CI * rp ex lg X y) *
      (k * (1 + Wsum ex lg X 0 S d κ) * (1 + Onsum S κ) ^ (k - 1))) + E =
      (𝒯.card : ℝ) * (Kpk * CI * k * (rp ex lg X y *
        ((1 + Wsum ex lg X 0 S d κ) * (1 + Onsum S κ) ^ (k - 1)))) + E := by ring
  rw [e] at hZ
  exact final_shape hc (mul_nonneg (mul_nonneg hKpk hCI) (Nat.cast_nonneg k)) hQ hE hZ

end DPCBLean

