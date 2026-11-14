FORMALIZED: Theorem E1 — for every θ ∈ (1/2, 1), liminf S/N ≥ h(θ; κ₁) ≥ h(1/2; κ₁) > 0.0055 — as a sorry-free Lean theorem (`E1Lean.theoremE1`) from the named inputs H1 (E1-M moment), H2 (SL Thm 1.3), H3 (SL Thm 5.1) and H4 (Arb certificate for c(P₁, Q₁, 177/20, 3/25)); H3 is itself derived (`sl51`, `theoremE1_parity`) from the parity product per window, Wang's pair-sum limit and the optimal-density limit; all algebra and numerics proved from Mathlib's elementary-function theorems (fields of `Fns`).

# E1 Lean formalization (e1lean, Orca, 12–14 Nov 2026)

Files: `Riemann/Deep/SantiAudit/Lean/` — Defs, Detector, AlgebraH, AlgebraC, Numeric, Kappa, Arch, Parity, Bounds2, Check (compiled); Inst (uncompiled intended instance, 54 lines). Namespace `E1Lean`. Standard axioms only; no sorry, admit, axiom, decide, native_decide.

## Short answer
| Part of the brief | Status | Lean |
|---|---|---|
| (0) definitions h, c, k*, κ (E1.md §1, §4) | done | `hc`, `h`, `cW`, `kstar`, `kappa`, `cot` (Defs) |
| (1) architecture: E1 from H1–H4 | proved | `theoremE1`; non-vacuity `inputs_consistent` |
| (2) h > 0 ⟺ k > k*; h ↑ in k, c; c ↑ and < 1; ν, exponent bookkeeping | proved | `hc_pos_iff`, `h_pos_iff`, `hc_mono_k`, `hc_mono_c`, `cW_lt`, `cW_lt_one`, `h_mono_theta`, `nu_range_iff`, `exponent_iff` |
| (3) κ₁ − k*(1/2) > 0, h(1/2; κ₁) > 0.0055 from H4 | proved, with two-sided enclosures | `kappa_bounds`, `kappa1_gt_kstar`, `h_half_gt`; `e1_enclosures`: c(1/2) ∈ [−0.16596709, −0.16596366], κ₁ ∈ [0.08447231, 0.08448814], h(1/2; κ₁) ∈ (0.0055, 0.005543], κ₁ − k*(1/2) ∈ (0.0078, 0.00787] (rational Taylor bounds, no numeric axiom) |
| (4) finite parity part of SL Thm 5.1 | partly | proved: SL Thm 5.1's deduction from Prop 2.5 + Wang + optimum (`sl51`), R(2 − c, 1 − k) = h (`Rroot_eq`), Lamzouri bound ⟹ parity product (`parity_of_lamzouri`); not: the Lamzouri bound itself (Hilbert space, Gram–Schmidt) |
| detector admissibility | proved | `detector_admissible`: P₁, Q₁ polynomials, P(0)=0, P(1)=1, Q(0)=1, Q(1)=0, Q(y)+Q(1−y) ≡ 1 (for any parameters x_j) |

## Named hypotheses and their sources
| Name | Statement (Lean: Arch.lean, Parity.lean) | Source | Kind |
|---|---|---|---|
| H1 | ∀ θ ∈ (1/2,1), 0 < ν < (3θ−1)/4: SL Thm 1.2's conclusion holds for every admissible (R, P, Q) | E1.md §1 (Thm E1-M), §2 Steps 1–5 (new; reviewer CORRECT, E1.md §7) | analytic, new |
| H2 | θ ∈ (1/2,1), 0 < ν < 1/2, (R,P,Q) admissible, Q(x)+Q(1−x) ≡ C ≠ 0, Thm 1.2 at (θ,ν) ⟹ liminf O/N ≥ κ(P,Q,R,ν) | SL Thm 1.3 (Zenodo 23134787 v0.03, pdftotext -layout l. 125–129; proof §4.3, l. 814–836) | analytic, attributed |
| H3 | θ ∈ (1/2,1), 0 ≤ k < 1, liminf O/N ≥ k ⟹ liminf S/N ≥ h(θ; k) | SL Thm 5.1 (l. 844–845; proof l. 854–869) | analytic, attributed |
| SL51Inputs (replaces H3) | counts N > 0, S ≥ 0, O ≤ N eventually; (Q_f − S)(N − O) ≥ 2(N − S)² eventually, ∀ f; limsup Q_f/N ≤ C(f); ∀ δ > 0 ∃ f, C(f) ≤ 2 − c(θ) + δ | Riemann–von Mangoldt; SL ref [4] Prop 2.5 (Zenodo 22860012 v0.07; SL l. 861); Wang arXiv 2609.07918 Thm 2.2 + (eq:unweighted) (SL l. 861–862); cosine optimum, λ ↑ θ (SL l. 868–869) | finite + analytic |
| H4 | 3302.286923204427565638 ≤ c(P₁, Q₁, 177/20, 3/25) ≤ 3302.286923204427565640 | e1_constants.py/.log (python-flint Arb, 4000 bits): c = [3302.286923204427565639012 ± 4.36e−22] | certificate |
| `Fns` fields | √, sin, cos, exp, log facts (below) | Mathlib theorems of this toolchain (names, binders checked in .lake/packages/mathlib) | Mathlib |

`Fns` fields = Mathlib theorems (namespace Real): sqrt_nonneg, mul_self_sqrt, sin_bound (|sin x − (x − x³/6)| ≤ |x|⁵/100), cos_bound (≤ |x|⁴·5/96), sin_two_mul, cos_two_mul, sin_sq_add_cos_sq, sin_sub, sin_lt, sin_gt_sub_cube (x − x³/6 < sin x), exp_pos, exp_zero, exp_add, exp_one_gt_d9, exp_one_lt_d9, quadratic_le_exp_of_nonneg, exp_bound_div_one_sub_of_interval (e^x ≤ 1/(1 − x), 0 ≤ x < 1), log_le_iff_le_exp; plus `rpow` (T^θ, a bare symbol). cot := cos/sin (Real.cot_eq_cos_div_sin). Instance: Inst.lean (`Fns.real`, uncompiled: Trigonometric.Bounds, ExponentialBounds, Log.Basic, Pow.Real exceed the 300 MB cap).

## Main statements (verbatim; proofs in the files)
```lean
theorem theoremE1 (h1 : H1 I) (h2 : H2 F I) (h3 : H3 F I) (h4 : H4 I) (θ : ℝ) (hθ1 : 1 / 2 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₁ F I)) ∧ h F (1 / 2) (κ₁ F I) ≤ h F θ (κ₁ F I) ∧
      (0.0055 : ℝ) < h F (1 / 2) (κ₁ F I)
theorem theoremE1_parity (I : Inputs) (h1 : H1 I) (h2 : H2 F I) (hin : ∀ θ : ℝ, 1 / 2 < θ → θ < 1 → SL51Inputs F I θ)
    (h4 : H4 I) (θ : ℝ) (hθ1 : 1 / 2 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₁ F I)) ∧ h F (1 / 2) (κ₁ F I) ≤ h F θ (κ₁ F I) ∧
      (0.0055 : ℝ) < h F (1 / 2) (κ₁ F I)
theorem e1_enclosures (I : Inputs) (h4 : H4 I) :
    (-16596708521464472031 / 100000000000000000000 : ℝ) ≤ cW F (1 / 2) ∧ cW F (1 / 2) ≤ (-16596366764648141727 / 100000000000000000000 : ℝ) ∧
    (37379 / 442500 : ℝ) ≤ κ₁ F I ∧ κ₁ F I ≤ 6231 / 73750 ∧
    (0.0055 : ℝ) < h F (1 / 2) (κ₁ F I) ∧ h F (1 / 2) (κ₁ F I) ≤ 0.005543 ∧
    (0.0078 : ℝ) < κ₁ F I - kstar F (1 / 2) ∧ κ₁ F I - kstar F (1 / 2) ≤ 0.00787
theorem kappa1_gt_kstar (h4 : H4 I) {θ : ℝ} (hθ1 : 1 / 2 ≤ θ) (hθ2 : θ ≤ 1) :
    0 < κ₁ F I - kstar F (1 / 2) ∧ kstar F θ < κ₁ F I
theorem inputs_consistent : H1 trivialInputs ∧ H2 F trivialInputs ∧ H3 F trivialInputs ∧ H4 trivialInputs
theorem h_half_gt {κ : ℝ} (hk : (37379 / 442500 : ℝ) ≤ κ) (hk1 : κ < 1) : (0.0055 : ℝ) < h F (1 / 2) κ
theorem h_half_le {κ : ℝ} (hk : κ ≤ 6231 / 73750) : h F (1 / 2) κ ≤ 0.005543
theorem kstar_half_lt {κ : ℝ} (hk : (37379 / 442500 : ℝ) ≤ κ) (hk1 : κ < 1) : kstar F (1 / 2) < κ
theorem kappa_bounds {c : ℝ} (hlo : (3302286923204427565638 / 1000000000000000000 : ℝ) ≤ c)
    (hhi : c ≤ 3302286923204427565640 / 1000000000000000000) :
    (37379 / 442500 : ℝ) ≤ kappa F c (177 / 20) ∧ kappa F c (177 / 20) < 1
theorem kappa_upper {c : ℝ} (hlo : (3302286923204427565638 / 1000000000000000000 : ℝ) ≤ c) :
    kappa F c (177 / 20) ≤ 6231 / 73750
theorem exp_810242 : (3302286923204427565640 / 1000000000000000000 : ℝ) ≤ F.exp (810242 / 100000)
theorem exp_Llo_le : F.exp (8 + 2557 / 25000) ≤ (3302286923204427565638 / 1000000000000000000 : ℝ)
theorem le_log_of_exp_le {c L : ℝ} (hc : 0 < c) (hle : F.exp L ≤ c) : L ≤ F.log c
theorem cW_half_ge : (-16596708521464472031 / 100000000000000000000 : ℝ) ≤ cW F (1 / 2)
theorem cW_half_le : cW F (1 / 2) ≤ (-16596366764648141727 / 100000000000000000000 : ℝ)
theorem trig_xlo : (34623307866193455119 / 100000000000000000000 : ℝ) ≤ F.sin (353553 / 1000000) ∧ F.cos (353553 / 1000000) ≤ (46907426160882655471 / 50000000000000000000 : ℝ)
theorem trig_xhi : (93814832600643526557 / 100000000000000000000 : ℝ) ≤ F.cos (3535534 / 10000000) ∧ F.sin (3535534 / 10000000) ≤ (3462336234640537867 / 10000000000000000000 : ℝ) ∧ 0 < F.sin (3535534 / 10000000)
theorem sqrt2_bounds : (141421356237 / 100000000000 : ℝ) ≤ F.sqrt 2 ∧ F.sqrt 2 ≤ 141421356238 / 100000000000
theorem hc_pos_iff {c k : ℝ} (hk : k < 1) (hc1 : c < 1) : 0 < hc F c k ↔ 1 - 2 / (2 - c) < k
theorem h_pos_iff {θ k : ℝ} (h0 : 0 < θ) (h1 : θ ≤ 1) (hk : k < 1) : 0 < h F θ k ↔ kstar F θ < k
theorem hc_mono_k {c k₁ k₂ : ℝ} (hc1 : c < 1) (hk : k₁ ≤ k₂) (hk2 : k₂ < 1) : hc F c k₁ ≤ hc F c k₂
theorem hc_mono_c {c₁ c₂ k : ℝ} (hk : k < 1) (hcc : c₁ ≤ c₂) (hc2 : c₂ < 1) : hc F c₁ k ≤ hc F c₂ k
theorem cot_sub_gt {x y : ℝ} (hx : 0 < x) (hxy : x < y) (hy : y ≤ 1) : y - x < cot F x - cot F y
theorem cW_lt {θ₁ θ₂ : ℝ} (h1 : 0 < θ₁) (h12 : θ₁ < θ₂) (h2 : θ₂ ≤ 1) : cW F θ₁ < cW F θ₂
theorem cW_lt_one {θ : ℝ} (h0 : 0 < θ) (h1 : θ ≤ 1) : cW F θ < 1
theorem h_mono_theta {θ₁ θ₂ k : ℝ} (h1 : 0 < θ₁) (h12 : θ₁ ≤ θ₂) (h2 : θ₂ ≤ 1) (hk : k < 1) :
    h F θ₁ k ≤ h F θ₂ k
theorem nu_range_iff (θ : ℝ) : (3 / 25 : ℝ) < (3 * θ - 1) / 4 ↔ 37 / 75 < θ
theorem nu2_range_iff (θ : ℝ) : (1249 / 10000 : ℝ) < (3 * θ - 1) / 4 ↔ 3749 / 7500 < θ
theorem exponent_iff (θ ν : ℝ) : 1 / 2 + 2 * ν - 3 * θ / 2 < 0 ↔ ν < (3 * θ - 1) / 4
theorem range_lt_half {θ : ℝ} (h : θ < 1) : (3 * θ - 1) / 4 < 1 / 2
theorem sqrt_eq_of {a s : ℝ} (hs : 0 ≤ s) (h : s * s = a) : F.sqrt a = s
theorem detector_admissible : IsPoly P₁ ∧ IsPoly Q₁ ∧ P₁ 0 = 0 ∧ P₁ 1 = 1 ∧ Q₁ 0 = 1 ∧ Q₁ 1 = 0 ∧
    ∀ y, Q₁ y + Q₁ (1 - y) = 1
theorem Rroot_eq (c k : ℝ) : Rroot F (2 - c) (1 - k) = hc F c k
theorem sl51 (I : Inputs) {ι : Type} (Qf : ι → ℝ → ℝ → ℝ) (C : ι → ℝ) {θ k : ℝ}
    (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hk1 : k < 1)
    (hcount : ∃ T₀ : ℝ, ∀ T ≥ T₀, 0 < I.N T (F.rpow T θ) ∧ 0 ≤ I.S T (F.rpow T θ) ∧
      I.O T (F.rpow T θ) ≤ I.N T (F.rpow T θ))
    (hpar : ∀ f, ∃ T₀ : ℝ, ∀ T ≥ T₀, 2 * (I.N T (F.rpow T θ) - I.S T (F.rpow T θ)) ^ 2 ≤
      (Qf f T (F.rpow T θ) - I.S T (F.rpow T θ)) * (I.N T (F.rpow T θ) - I.O T (F.rpow T θ)))
    (hwang : ∀ f, ∀ ε > 0, ∃ T₀ : ℝ, ∀ T ≥ T₀, ratio F (Qf f) I.N θ T ≤ C f + ε)
    (hopt : ∀ δ > 0, ∃ f, C f ≤ 2 - cW F θ + δ)
    (hO : LiminfGe (ratio F I.O I.N θ) k) : LiminfGe (ratio F I.S I.N θ) (h F θ k)
theorem H3_of_inputs (I : Inputs) (hin : ∀ θ : ℝ, 1 / 2 < θ → θ < 1 → SL51Inputs F I θ) : H3 F I
theorem parity_of_lamzouri {Q N S O d : ℝ} (hd : 0 ≤ d) (hNS : d ≤ N - S) (hNO : 2 * d ≤ N - O)
    (hL : ∀ t : ℝ, 1 ≤ t → 2 * t * N - (2 * t - 1) * S - t ^ 2 * d ≤ Q) :
    2 * (N - S) ^ 2 ≤ (Q - S) * (N - O)
```

## Definitions (verbatim; `variable (F : Fns) (I : Inputs)`; docstrings in the files)
```lean
noncomputable def cot (x : ℝ) : ℝ := F.cos x / F.sin x
noncomputable def cW (θ : ℝ) : ℝ := 2 - θ / 2 - cot F (θ / F.sqrt 2) / F.sqrt 2
noncomputable def hc (c k : ℝ) : ℝ := (3 + k - F.sqrt ((1 - k) * (9 - k - 8 * c))) / 4
noncomputable def h (θ k : ℝ) : ℝ := hc F (cW F θ) k
noncomputable def kstar (θ : ℝ) : ℝ := 1 - 2 / (2 - cW F θ)
noncomputable def kappa (c R : ℝ) : ℝ := 1 - R⁻¹ * F.log c
def LiminfGe (f : ℝ → ℝ) (a : ℝ) : Prop := ∀ ε > 0, ∃ T₀ : ℝ, ∀ T ≥ T₀, a - ε ≤ f T
noncomputable def ratio (A N : ℝ → ℝ → ℝ) (θ : ℝ) : ℝ → ℝ := fun T => A T (F.rpow T θ) / N T (F.rpow T θ)
noncomputable def qsum : List ℝ → ℕ → ℝ → ℝ
  | [], _, _ => 0
  | a :: as, j, y => a * (cheb (2 * j + 1) (1 - 2 * y) - (1 - 2 * y)) + qsum as (j + 1) y
noncomputable def Qfam (xs : List ℝ) (y : ℝ) : ℝ := 1 - y + qsum xs 1 y
noncomputable def Sfun (u : ℝ) : ℝ :=
  u + u ^ 3 / 6 + u ^ 5 / 120 + u ^ 7 / 5040 + u ^ 9 / 362880 + u ^ 11 / 39916800 + u ^ 13 / 6227020800
noncomputable def Pfam (r x : ℝ) : ℝ := Sfun (r * x) / Sfun r
noncomputable def Q₁ : ℝ → ℝ := Qfam xs₁
noncomputable def P₁ : ℝ → ℝ := Pfam (7707 / 10000)
def Adm (R : ℝ) (P Q : ℝ → ℝ) : Prop := 0 < R ∧ IsPoly P ∧ IsPoly Q ∧ P 0 = 0 ∧ P 1 = 1 ∧ Q 0 = 1
def MomentRange (θ ν : ℝ) : Prop := ∀ R P Q, Adm R P Q → I.Moment θ ν R P Q
def H1 : Prop := ∀ θ ν : ℝ, 1 / 2 < θ → θ < 1 → 0 < ν → ν < (3 * θ - 1) / 4 → MomentRange I θ ν
def H2 : Prop := ∀ (θ ν R : ℝ) (P Q : ℝ → ℝ), 1 / 2 < θ → θ < 1 → 0 < ν → ν < 1 / 2 → Adm R P Q →
  (∃ C : ℝ, C ≠ 0 ∧ ∀ x, Q x + Q (1 - x) = C) → MomentRange I θ ν →
  LiminfGe (ratio F I.O I.N θ) (kappa F (I.cdet P Q R ν) R)
def H3 : Prop := ∀ θ k : ℝ, 1 / 2 < θ → θ < 1 → 0 ≤ k → k < 1 →
  LiminfGe (ratio F I.O I.N θ) k → LiminfGe (ratio F I.S I.N θ) (h F θ k)
def H4 : Prop := (3302286923204427565638 / 1000000000000000000 : ℝ) ≤ I.cdet P₁ Q₁ (177 / 20) (3 / 25) ∧
  I.cdet P₁ Q₁ (177 / 20) (3 / 25) ≤ 3302286923204427565640 / 1000000000000000000
noncomputable def κ₁ : ℝ := kappa F (I.cdet P₁ Q₁ (177 / 20) (3 / 25)) (177 / 20)
noncomputable def Rroot (A b : ℝ) : ℝ := (4 - b - F.sqrt (b * (b + 8 * (A - 1)))) / 4
def SL51Inputs (I : Inputs) (θ : ℝ) : Prop :=
  ∃ (ι : Type) (Qf : ι → ℝ → ℝ → ℝ) (C : ι → ℝ),
    (∃ T₀ : ℝ, ∀ T ≥ T₀, 0 < I.N T (F.rpow T θ) ∧ 0 ≤ I.S T (F.rpow T θ) ∧
      I.O T (F.rpow T θ) ≤ I.N T (F.rpow T θ)) ∧
    (∀ f, ∃ T₀ : ℝ, ∀ T ≥ T₀, 2 * (I.N T (F.rpow T θ) - I.S T (F.rpow T θ)) ^ 2 ≤
      (Qf f T (F.rpow T θ) - I.S T (F.rpow T θ)) * (I.N T (F.rpow T θ) - I.O T (F.rpow T θ))) ∧
    (∀ f, ∀ ε > 0, ∃ T₀ : ℝ, ∀ T ≥ T₀, ratio F (Qf f) I.N θ T ≤ C f + ε) ∧
    (∀ δ > 0, ∃ f, C f ≤ 2 - cW F θ + δ)
```
`Inputs` = ⟨N, S, O : ℝ → ℝ → ℝ; Moment : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop; cdet : (ℝ → ℝ) → (ℝ → ℝ) → ℝ → ℝ → ℝ⟩ (abstract; meanings in Arch.lean docstrings). `IsPoly` = inductive closure of constants and x under +, × (= polynomial functions ℝ → ℝ). `cheb` = Chebyshev T_n by T₀ = 1, T₁ = x, T_{n+2} = 2xT_{n+1} − T_n. `xs₁` = the 100 x_j of params_exp010.json (set ν = 999/10000), verbatim.

## Fidelity
SL = Santibáñez-Leal, Zenodo 23134787 v0.03 (`pdftotext -layout`, line numbers l.).
| Lean | E1.md / SL | Check |
|---|---|---|
| `cW θ = 2 − θ/2 − cot(θ/√2)/√2`, `cot = cos/sin` | E1.md §1; SL (1) l. 54–56 | ✓ verbatim; cot = Real.cot (Real.cot_eq_cos_div_sin); `sqrt_eq_of`: the two √ fields force F.sqrt = Real.sqrt on [0, ∞) |
| `h θ k = hc (cW θ) k = (3 + k − √((1 − k)(9 − k − 8c)))/4` | E1.md §1; SL (2) l. 58–63 | ✓ |
| `kstar θ = 1 − 2/(2 − cW θ)` | E1.md §4; SL l. 66 | ✓ |
| `kappa c R = 1 − R⁻¹ log c`; `κ₁ = kappa (cdet P₁ Q₁ (177/20) (3/25)) (177/20)` | SL Thm 1.3 l. 128; E1.md §4 | ✓ |
| `LiminfGe f a` = ∀ ε > 0, eventually a − ε ≤ f T | liminf_{T→∞} f ≥ a | ✓ (extended-real liminf; for the bounded ratios = Mathlib's a ≤ liminf f atTop) |
| `ratio A N θ T = A(T, T^θ)/N(T, T^θ)` | S(T,T^θ)/N(T,T^θ), O/N | ✓ (x/0 = 0 immaterial: N > 0 eventually) |
| `Qfam xs₁`, `Pfam (7707/10000)`, `Sfun` | SL §6 l. 878–887; params_exp010.json (set ν = 999/10000); e1_constants.py | ✓ same family, data, index shift (x_j ↦ T_{2j+1}, j ≥ 1), degree-13 sinh Taylor S |
| `Adm`, `detector_admissible` | SL Thm 1.2 l. 98–100 (P(0)=0, P(1)=1, Q(0)=1); Thm 1.3 l. 125–126 (Q(x)+Q(1−x) nonzero constant) | ✓ proved for all x_j (SL l. 883–886) |
| `Inputs.Moment` (opaque) | SL Thm 1.2's display l. 115–123 for every admissible weight w | named; docstring only |
| `Inputs.cdet` (opaque) | c(P,Q,R,ν), SL l. 120–123 | named; only H2 and H4 mention it, with the same arguments |
| `IsPoly` | "real polynomial" (SL l. 98–100) | ✓ closure of constants, x under +, × |

Hypotheses against their sources:
- H1 = Theorem E1-M (E1.md §1): SL Thm 1.2 for 0 < ν < (3θ − 1)/4, θ ∈ (1/2, 1), all admissible (R, P, Q). Faithful.
- H2 is SL Thm 1.3 at every (θ, ν) where Thm 1.2's conclusion holds (`MomentRange θ ν`), not only for ν in SL's range (5): E1 uses it at ν = 3/25, outside (5) near θ = 1/2. This is the audited content of SL §4 (FULLPROOF.md §3b; E1.md §3 row "Lemma 4.1, 4.2, Thm 1.3 | none": the proof uses only Thm 1.2's conclusion and (22) = Thm 1.2 for Q = 1, 1 + x^j, all admissible). The extra premises ν < 1/2, θ ∈ (1/2, 1) only weaken H2. Stronger than SL's printed statement exactly in this way.
- H3 = SL Thm 5.1 (l. 844–845) with the extra premises θ ∈ (1/2, 1), k < 1 (weaker; E1.md §3 row "κ < 1, c < 1").
- SL51Inputs: count form of SL l. 861 ((q_T − s_T)(1 − o_T) ≥ 2(1 − s_T)² times N²), eventually in T, for every test function (Prop 2.5 holds for every finite conjugation-invariant multiset, so this is weaker than the source); only the upper half of Wang's q_T → C(f); "C(f) → 2 − c(λ), then λ → θ" (l. 868–869) as ∀ δ > 0 ∃ f, C(f) ≤ 2 − c(θ) + δ (continuity of c in λ is inside this input). Proved from them: `sl51` (= H3), with the ε-bookkeeping and the smaller root R(A, b) of SL l. 864–867.
- H4: the Arb ball of e1_constants.log rounded outward to 18 decimals.
- Brief's "c increasing given cot decreasing" is not enough: c′ > 0 needs x + cot x decreasing (cot′ = −1 − cot²); Lean proves y − x < cot x − cot y on 0 < x < y ≤ 1 (`cot_sub_gt`) from sin_sub, sin_lt, sin_gt_sub_cube.
- Non-vacuity: `inputs_consistent` gives, for every `F`, an `Inputs` satisfying H1–H4 (N = S = O = 1), so the hypotheses are not contradictory; their truth for ζ is the analytic content.

Numerics (Lean-proved, rational): √2 ∈ [1.41421356237, 1.41421356238]; sin, cos at 353553/10⁶ from Taylor at 353553/(8·10⁶) and three doublings; cot(1/(2√2)) ∈ [2.70958180, 2.70958664]; c(1/2) ∈ [−0.165967085215, −0.165963667646] (Arb: −0.165963850365); e^{8.10242} ≥ c_hi (e > 2.7182818283, e^x ≥ 1 + x + x²/2) and e^{8.10228} ≤ c_lo (e < 2.7182818286, e^x ≤ 1/(1 − x) on 64 pieces); κ₁ ∈ [37379/442500, 6231/73750] = [0.08447232, 0.08448814] (Arb: 0.08447791); h(1/2; κ₁) ∈ (0.0055, 0.005543] (Arb: 0.0055353); κ₁ − k*(1/2) ∈ (0.0078, 0.00787] (Arb: 0.0078544). The decimals labelled Arb are e1_constants.log values, not Lean claims.

Not formalized (named inputs only): E1-M's proof (Mellin–Parseval, Bessel/Hankel, integration by parts, Estermann FE, Young's main term); SL §4 (Littlewood, argument principle); Wang's theorem; the Lamzouri bound Q ≥ 2tN − (2t − 1)S − t²d (Hilbert space; its consequence `parity_of_lamzouri` is proved); the double integral c(P, Q, R, ν) (opaque; H4 is the Arb certificate); the count asymptotics.

## Build
Mathlib: `rm .lake && bash .scratch/memcap.sh 290 3000 lake exe cache get` (12 Nov 22:20 UTC; the machine was suspended 12 Nov 22:38 → 14 Nov 17:23 UTC and memcap's wall clock killed the run on resume, TIME); rerun `bash .scratch/memcap.sh 290 1500 lake exe cache get` (14 Nov 17:24): 7555 already decompressed, 1360 decompressed, no downloads, peak 125 MB, rc = 0. LEAN_PATH = `lake env printenv LEAN_PATH` (Riemann/Deep/SantiAudit/Lean/audit/leanpath.txt). Imports: Mathlib.Basic.Real.Basic (Data.Real.Basic is deprecated in this Mathlib), Tactic.Linarith, Positivity, Ring, FieldSimp, NormNum, LinearCombination. Each module separately, one at a time, lean directly (no lake build), oleans to .scratch/e1lean/olean (build output, deleted at the end):

    LEAN_PATH=.scratch/e1lean/olean:$(cat Riemann/Deep/SantiAudit/Lean/audit/leanpath.txt) OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 1800 \
      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 \
      -o .scratch/e1lean/olean/Riemann/Deep/SantiAudit/Lean/$m.olean Riemann/Deep/SantiAudit/Lean/$m.lean

Scripts: Lean/audit/mk.sh (one module), verify.sh (clean rebuild), gendet.py (xs₁ from params_exp010.json), gennum.py (Numeric's rational enclosures; Bounds2's trig_xhi block was generated the same way at x = 3535534/10⁷, exporting the opposite ends). Final clean rebuild (14 Nov, system load ≈ 15, background QoS), every module rc = 0:

| module | wall | footprint (time -l) | memcap peak (top MEM) |
|---|---|---|---|
| Defs | 29 s | 243 MiB | 173 MB |
| Detector | 45 s | 253 MiB | 244 MB |
| AlgebraH | 63 s | 252 MiB | 235 MB |
| AlgebraC | 106 s | 281 MiB | 272 MB |
| Numeric | 53 s | 262 MiB | 240 MB |
| Kappa | 53 s | 248 MiB | 229 MB |
| Arch | 62 s | 250 MiB | 231 MB |
| Parity | 105 s | 284 MiB | 257 MB |
| Bounds2 | 72 s | 274 MiB | 249 MB |
| Check | 37 s | 239 MiB | 207 MB |

Peaks over the whole session: footprint ≤ 294 MiB (one early unsplit Algebra.lean, top MEM 274 MB; split into AlgebraH/AlgebraC), top MEM ≤ 274 MB; all under the 300 MB cap.

## Static scan
`rg -n -w 'sorry|admit|native_decide|skipKernelTC|decide' Riemann/Deep/SantiAudit/Lean/*.lean` and `rg -n '^\s*axiom\s' Riemann/Deep/SantiAudit/Lean/*.lean`: no match (rc = 1).

## #print axioms (Check.lean, 39 lines, all identical)

    'E1Lean.theoremE1' depends on axioms: [propext, Classical.choice, Quot.sound]
    (same line for kappa1_gt_kstar, inputs_consistent, detector_admissible, hc_pos_iff, h_pos_iff, hc_mono_k, hc_mono_c, cot_sub_gt, cW_lt, cW_lt_one, h_mono_theta, nu_range_iff, nu2_range_iff, exponent_iff, range_lt_half, sqrt_eq_of, sqrt2_bounds, trig_xlo, cW_half_ge, exp_810242, kappa_bounds, h_half_gt, kstar_half_lt, Rroot_eq, root_step, per_T, ratio_parity, sl51, H3_of_inputs, theoremE1_parity, parity_of_lamzouri, trig_xhi, cW_half_le, le_log_of_exp_le, exp_Llo_le, kappa_upper, h_half_le, e1_enclosures)

## Review (statement fidelity, verdict-only reviewer subagents; prompts and answers in Lean/audit/)
- Pass 1 (E1.md §1 statement, SL Thm 1.2/1.3/5.1 and the Prop 2.5 step pasted with all Lean definitions, `Fns`, `Inputs`, H1–H4, main statements): **ISSUES** — (1) conditional reduction, (2) abstract `Inputs`/`rpow`, (7) H4 is a hypothesis: by the brief's design (named inputs); (3) H2 extends SL Thm 1.3 beyond SL's range (5): flagged in Fidelity (audited content of SL §4); (4) parity required for every T: fixed, now eventual; (5) `xs₁` not auditable from the paste: checked mechanically, its 100 rationals equal params_exp010.json's x_j (set ν = 999/10000) in order, P₁ uses that set's r_P = 7707/10000, and e1_constants.py (H4) reads the same set; (6) numerics weaker than E1.md's decimals: two-sided enclosures added (`e1_enclosures`).
- Pass 2 (same paste + these replies): **FAITHFUL**.
- Pass 3 (`e1_enclosures` against E1.md §4's Arb values): **FAITHFUL** — all Arb values inside the bounds, directions correct; conditional on H4.

## Disk and .lake
.lake was the placeholder file; deleted and filled by `lake exe cache get` (allowed by the brief). Left as it is: a real directory with Mathlib and its dependencies (7.1 GB); the orchestrator restores the placeholder. No `lake build` was run. Build output .scratch/e1lean/olean and the downloads (SL PDF, its text, the Zenodo JSON) deleted at the end.

## FORMALIZATIONS.md entry
- Theorem E1 (simple critical zeros in every window (T, T + T^θ], θ > 1/2; N8) formalized as an architecture: `E1Lean.theoremE1` — liminf S/N ≥ h(θ; κ₁) ≥ h(1/2; κ₁) > 0.0055 for every θ ∈ (1/2, 1) — from named inputs H1 (E1-M moment), H2 (SL Thm 1.3), H3 (SL Thm 5.1, itself derived from the parity product, Wang's pair sum and the optimal density: `theoremE1_parity`), H4 (Arb certificate); algebra, detector admissibility, SL Thm 5.1's deduction and the E1.md §4 enclosures proved; sorry-free, standard axioms, Mathlib real functions as named fields (instance uncompiled). [LEAN.md](Riemann/Deep/SantiAudit/LEAN.md)
