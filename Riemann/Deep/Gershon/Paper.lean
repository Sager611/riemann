import Mathlib

/-!
# Gershon, "The De Bruijn–Newman Constant Is Zero" (Preprints.org 202604.1513 v1, 22 Apr 2026)

Paper's chain to RH (Λ = 0), §11 (Parts I–III are not used for the global claim; the paper itself
says, after Thm 1.1, that log-concavity does not imply Laguerre–Pólya, cf. `Riemann/Tree/R038.lean`):
* `γ_k = ∫₀^∞ Φ(u) u^{2k}/(2k)! du`, `D_r(n) = det (γ_{n+i-j})_{i,j<r}` (Def. 11.1).
* Thm 11.1 (Edrei–Schoenberg): `D_r(n) > 0 ∀ r n ⇒ Ξ ∈ LP` (correct, cited) — `EdreiSchoenberg`.
* Thm 11.3: regions A (r ≤ 2), B (3 ≤ r ≤ 50), C2 (r ≥ 51, n ≤ 99), C1 (r ≥ 51, n ≥ 100)
  cover the plane — `gershon_RH` (proved).
* Region C1 = Prop. 11.12 case n ≥ 100: the unitarity ratio `R(n) = S/(a n) + ε_n`, with
  `S = Σ_s C_s` computed at n = 30 and transferred to every n ≥ 100 by **Lemma 11.6**:
  `|C_s(n) - C_s^∞| ≤ K_s δⁿ`, `Σ K_s < 2`, `δ = 0.452`, where `C_s(n) = n² |log Θ_s(n)|`.

Failing step: Lemma 11.6 is false. Arb (`.scratch/gershon/arb_c1.py`) encloses
`C_1(50) - C_1(30) = 0.0331832194541477343… ± 7e-23`, while Lemma 11.6 allows `≤ 2(δ³⁰+δ⁵⁰) ≈ 9e-11`
(`lemma116_false`). Its proof is also void: it expands the Taylor coefficients of the ENTIRE
`g(z) = Σ γ_k z^k` as `γ_m = R₁ρ₁^m + R₂ρ₂^m + O(δ₃^m ρ₁^m)` (that is the expansion of `1/g`);
for an entire series this forces `R₁ = 0` (`expansion_forces_R1_zero`), and `K_s = O(|R₂/R₁|^s)`
is meaningless. The "analytical alternative" (Cor. 11.2) assumes the zeros `-t_m²` of `g` are real
and ordered (i.e. RH) and proves only `D_r(n) > 0 ∀ r ⇒ unitarity`, the converse of what C1 needs.
-/

namespace Riemann.Deep.Gershon

open Filter Topology Real MeasureTheory

/-- Riemann–Jacobi kernel `Φ(u) = Σ_{m≥1} (2π²m⁴e^{9u} - 3πm²e^{5u}) exp(-πm²e^{4u})`. -/
noncomputable def Phi (u : ℝ) : ℝ :=
  ∑' m : ℕ, (2 * π ^ 2 * ((m : ℝ) + 1) ^ 4 * exp (9 * u) - 3 * π * ((m : ℝ) + 1) ^ 2 * exp (5 * u)) *
    exp (-π * ((m : ℝ) + 1) ^ 2 * exp (4 * u))

/-- Gershon's coefficient sequence `γ_k = m_{2k}/(2k)!` (any u-scaling gives the same `L`, `Θ`). -/
noncomputable def gam (k : ℕ) : ℝ := (∫ u in Set.Ioi (0 : ℝ), Phi u * u ^ (2 * k)) / (2 * k).factorial

/-- Contiguous Toeplitz minor `D_r(n) = det (γ_{n+i-j})_{0 ≤ i,j < r}`, `γ_k = 0` for `k < 0`. -/
noncomputable def D (γ : ℕ → ℝ) (r n : ℕ) : ℝ :=
  (Matrix.of fun i j : Fin r => if (j : ℕ) ≤ n + i then γ (n + i - j) else 0).det

/-- `D_r(n) > 0` for all `r ≥ 1`, `n ≥ 0` (the paper's `PF_∞`). -/
def PF (γ : ℕ → ℝ) : Prop := ∀ r n : ℕ, 1 ≤ r → 0 < D γ r n

/-- Thm 11.1 (Edrei–Schoenberg; correct, cited). -/
def EdreiSchoenberg : Prop := PF gam → RiemannHypothesis

/-- Positivity of `D` on a region of the `(r, n)` plane. -/
def PosOn (P : ℕ → ℕ → Prop) : Prop := ∀ r n : ℕ, 1 ≤ r → P r n → 0 < D gam r n

def RegionA : Prop := PosOn fun r _ => r ≤ 2
def RegionB : Prop := PosOn fun r _ => 3 ≤ r ∧ r ≤ 50
def RegionC2 : Prop := PosOn fun r n => 51 ≤ r ∧ n ≤ 99
def RegionC1 : Prop := PosOn fun r n => 51 ≤ r ∧ 100 ≤ n

/-- Thm 11.3 + Cor 11.3, machine-checked: the four regions + Edrei–Schoenberg ⇒ RH. -/
theorem gershon_RH (ha : RegionA) (hb : RegionB) (hc2 : RegionC2) (hc1 : RegionC1)
    (hes : EdreiSchoenberg) : RiemannHypothesis := by
  refine hes fun r n hr => ?_
  rcases le_or_gt r 2 with h | h
  · exact ha r n hr h
  rcases le_or_gt r 50 with h' | h'
  · exact hb r n hr ⟨h, h'⟩
  rcases le_or_gt n 99 with h'' | h''
  · exact hc2 r n hr ⟨h', h''⟩
  · exact hc1 r n hr ⟨h', h''⟩

/-! ## Level-1 quantities (Def. 11.2) and Lemma 11.6 at s = 1 -/

/-- `L_1(n) = γ_n² / (γ_{n-1} γ_{n+1})`. -/
noncomputable def L1 (γ : ℕ → ℝ) (n : ℕ) : ℝ := γ n ^ 2 / (γ (n - 1) * γ (n + 1))

noncomputable def eps (γ : ℕ → ℝ) (n : ℕ) : ℝ := 1 - 1 / L1 γ n

/-- `Θ_1(n) = ε(n)² / (ε(n-1) ε(n+1))`. -/
noncomputable def Theta1 (γ : ℕ → ℝ) (n : ℕ) : ℝ := eps γ n ^ 2 / (eps γ (n - 1) * eps γ (n + 1))

/-- `C_1(n) = n² |log Θ_1(n)|`. -/
noncomputable def C1 (γ : ℕ → ℝ) (n : ℕ) : ℝ := (n : ℝ) ^ 2 * |Real.log (Theta1 γ n)|

/-- Lemma 11.6 at `s = 1`: `|C_1(n) - C_1^∞| ≤ K_1 δⁿ` for all `n` (paper: `K_1 ≤ Σ K_s < 2`). -/
def Lemma116 (γ : ℕ → ℝ) (Cinf K δ : ℝ) : Prop := ∀ n : ℕ, |C1 γ n - Cinf| ≤ K * δ ^ n

/-- Arb certificate (`.scratch/gershon/arb_c1.py`, 400 bits, rigorous tails):
`C_1(50) - C_1(30) ∈ 0.0331832194541477343 ± 7e-23`. Numerical fact; cited, not an axiom. -/
def Witness : Prop := 3 / 100 < C1 gam 50 - C1 gam 30

/-- The failing lemma is false. -/
theorem lemma116_false (hw : Witness) {Cinf K : ℝ} (hK : K ≤ 2) :
    ¬ Lemma116 gam Cinf K (452 / 1000) := by
  intro h
  have h30 := h 30
  have h50 := h 50
  have hsmall : 2 * ((452 / 1000 : ℝ) ^ 30 + (452 / 1000) ^ 50) < 3 / 100 := by norm_num
  have hKx : K * ((452 / 1000 : ℝ) ^ 30 + (452 / 1000) ^ 50) ≤
      2 * ((452 / 1000 : ℝ) ^ 30 + (452 / 1000) ^ 50) :=
    mul_le_mul_of_nonneg_right hK (by positivity)
  have h1 := le_abs_self (C1 gam 50 - Cinf)
  have h2 := neg_abs_le (C1 gam 30 - Cinf)
  unfold Witness at hw
  nlinarith

/-! ## Lemma 11.6's proof premise: a geometric expansion of entire coefficients -/

/-- Coefficients of an entire power series decay faster than every geometric sequence. -/
theorem entire_coeff_div_geom {γ : ℕ → ℝ} (hent : ∀ x : ℝ, Summable fun m => γ m * x ^ m)
    (ρ : ℝ) : Tendsto (fun m => γ m / ρ ^ m) atTop (𝓝 0) := by
  refine (hent ρ⁻¹).tendsto_atTop_zero.congr fun m => ?_
  rw [inv_pow, div_eq_mul_inv]

/-- Paper (proof of Lemma 11.6): `γ_m = R₁ρ₁^m + R₂ρ₂^m + E_m`, `|E_m| ≤ B (δ₃ρ₁)^m`, `ρ₂ < ρ₁`,
`δ₃ < 1`. For the entire `g` (paper, Prop. 11.5) this forces `R₁ = 0`, so the "common factor ρ₁ⁿ"
and `K_s = O(|R₂/R₁|^s)` are void. -/
theorem expansion_forces_R1_zero {γ : ℕ → ℝ} (hent : ∀ x : ℝ, Summable fun m => γ m * x ^ m)
    {R₁ R₂ ρ₁ ρ₂ δ₃ B : ℝ} (hρ₁ : 0 < ρ₁) (hρ₂ : 0 ≤ ρ₂) (h12 : ρ₂ < ρ₁) (hδ0 : 0 ≤ δ₃)
    (hδ : δ₃ < 1)
    (hexp : ∀ m : ℕ, |γ m - (R₁ * ρ₁ ^ m + R₂ * ρ₂ ^ m)| ≤ B * (δ₃ * ρ₁) ^ m) : R₁ = 0 := by
  have h0 := entire_coeff_div_geom hent ρ₁
  have hq : Tendsto (fun m : ℕ => R₂ * (ρ₂ / ρ₁) ^ m) atTop (𝓝 0) := by
    simpa using (tendsto_pow_atTop_nhds_zero_of_lt_one (div_nonneg hρ₂ hρ₁.le)
      ((div_lt_one hρ₁).2 h12)).const_mul R₂
  have hd : Tendsto (fun m : ℕ => B * δ₃ ^ m) atTop (𝓝 0) := by
    simpa using (tendsto_pow_atTop_nhds_zero_of_lt_one hδ0 hδ).const_mul B
  have hsum : Tendsto (fun m : ℕ => |R₂ * (ρ₂ / ρ₁) ^ m| + B * δ₃ ^ m) atTop (𝓝 0) := by
    simpa using hq.abs.add hd
  have hlim : Tendsto (fun m : ℕ => γ m / ρ₁ ^ m) atTop (𝓝 R₁) := by
    refine tendsto_iff_dist_tendsto_zero.2 (squeeze_zero (fun _ => dist_nonneg) (fun m => ?_) hsum)
    have hpos : 0 < ρ₁ ^ m := pow_pos hρ₁ m
    have he := hexp m
    rw [mul_pow] at he
    have key : γ m / ρ₁ ^ m - R₁ - R₂ * (ρ₂ / ρ₁) ^ m =
        (γ m - (R₁ * ρ₁ ^ m + R₂ * ρ₂ ^ m)) / ρ₁ ^ m := by
      rw [div_pow]; field_simp; ring
    have hb : |γ m / ρ₁ ^ m - R₁ - R₂ * (ρ₂ / ρ₁) ^ m| ≤ B * δ₃ ^ m := by
      rw [key, abs_div, abs_of_pos hpos, div_le_iff₀ hpos]
      linarith
    rw [Real.dist_eq]
    have := abs_add_le (γ m / ρ₁ ^ m - R₁ - R₂ * (ρ₂ / ρ₁) ^ m) (R₂ * (ρ₂ / ρ₁) ^ m)
    simp only [sub_add_cancel] at this
    linarith
  exact tendsto_nhds_unique hlim h0

/-! ## Second, independent check: the unitarity ratio (61) at n = 100

Telescoping (53): `μ_r = μ_1 + Σ_{1≤s<r} log Θ_s`, and (paper, Cor. 11.2) `μ_r(n) → 2 log(t_{n+1}/t_n)`.
So `Σ_s |log Θ_s(n)| ≥ μ_1 - μ_∞`: the ratio `R(n)` is at least `1 - μ_∞/μ_1`, a zero-gap quantity, not
`S/(a n)`. Arb (`arb_r100.py`): `μ_1(100) = 0.01512954850896…`, `t_100 = 236.5242296658…`,
`t_101 = 237.7698204809…`, `μ_∞/μ_1 = 0.69432490…`, so `R(100) ≥ 0.30567…` against the claimed `< 0.149`. -/

/-- Eq. (53), summed. -/
def Telescoping (mu lth : ℕ → ℝ) : Prop :=
  ∀ r : ℕ, 1 ≤ r → mu r = mu 1 + ∑ s ∈ Finset.Ico 1 r, lth s

/-- Prop. 11.12 / (61): every partial drag `Σ_{s<r} |log Θ_s|` is at most `R · μ_1`. -/
def Claim61 (mu lth : ℕ → ℝ) (R : ℝ) : Prop :=
  ∀ r : ℕ, ∑ s ∈ Finset.Ico 1 r, |lth s| ≤ R * mu 1

theorem claim61_false {mu lth : ℕ → ℝ} {muinf : ℝ} (ht : Telescoping mu lth)
    (hlim : Tendsto mu atTop (𝓝 muinf)) (hw : muinf < 4 / 5 * mu 1) (hmu : 0 < mu 1) :
    ¬ Claim61 mu lth (149 / 1000) := by
  intro hc
  have hev : ∀ᶠ r in atTop, (851 / 1000) * mu 1 ≤ mu r := by
    filter_upwards [eventually_ge_atTop 1] with r hr
    have h1 := ht r hr
    have h2 := hc r
    have h3 : -(∑ s ∈ Finset.Ico 1 r, lth s) ≤ ∑ s ∈ Finset.Ico 1 r, |lth s| := by
      rw [← Finset.sum_neg_distrib]
      exact Finset.sum_le_sum fun s _ => neg_le_abs (lth s)
    linarith
  have := ge_of_tendsto hlim hev
  linarith

/-- Summary: the chain is valid; its load-bearing transfer lemma (11.6) is refuted. -/
theorem gershon_gap :
    (RegionA → RegionB → RegionC2 → RegionC1 → EdreiSchoenberg → RiemannHypothesis) ∧
      (Witness → ∀ Cinf K : ℝ, K ≤ 2 → ¬ Lemma116 gam Cinf K (452 / 1000)) :=
  ⟨gershon_RH, fun hw _ _ hK => lemma116_false hw hK⟩

end Riemann.Deep.Gershon

#print axioms Riemann.Deep.Gershon.gershon_gap
#print axioms Riemann.Deep.Gershon.expansion_forces_R1_zero
#print axioms Riemann.Deep.Gershon.claim61_false
