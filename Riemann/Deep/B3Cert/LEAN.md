FORMALIZED: Theorem E1′ — for every θ ∈ (371/750, 1), liminf S/N ≥ h(θ; κ₂) ≥ h(371/750; κ₂) > 3.35·10⁻⁵ — as a sorry-free Lean theorem (`E1PLean.theoremE1P`) from the named inputs H1′ (E1-M′ moment, θ ∈ (1/3, 1)), H2′ (SL Thm 1.3 on θ ∈ (0, 1)), H3′ (SL Thm 5.1 on θ ∈ (0, 1)) and H4′ (Arb certificate for c(P₁, Q₁, 44/5, 121/1000)); H1′ is itself derived from the Gaussian-form moment plus SL's compact-window step (`H1P_of_gauss`, `theoremE1P_gauss`), H3′ from SL Thm 5.1's inputs (`H3P_of_inputs`, `theoremE1P_parity`, `theoremE1P_fine`); the edge algebra, all E1′ constants (two-sided enclosures) and the detector threshold θ_det ∈ (0.4946365, 0.4946366) are proved from Mathlib's elementary-function theorems (E1Lean's `Fns` fields).

# E1′ Lean formalization (b3lean, Orca, 20 Nov 2026)

Files: `Riemann/Deep/B3Cert/Lean/` — Trig, Exp, Num, Arch, Window, Parity, DetTrig, Det, Check (all compiled; 9 files, 757 lines). Namespace `E1PLean`; imports E1Lean (`Riemann/Deep/SantiAudit/Lean/`, read-only, unchanged) and reuses its definitions verbatim: `Fns`, `cot`, `cW`, `hc`, `h`, `kstar`, `kappa`, `LiminfGe`, `ratio`, `IsPoly`, `Inputs`, `Adm`, `MomentRange`, `P₁`, `Q₁`, `SL51Inputs`, and its lemmas (`detector_admissible`, `h_mono_theta`, `h_pos_iff`, `hc_mono_k/c`, `cW_lt`, `cW_lt_one`, `cot_sub_gt`, `exponent_iff`, `sqrt2_bounds`, `dbl`, `exp_nat`, `le_log_of_exp_le`, `sl51`). Intended instance of `Fns`: E1Lean's Inst.lean (`Fns.real`, uncompiled; no new field). Standard axioms only; no sorry, admit, axiom, decide, native_decide, skipKernelTC.

## Short answer
| Part of the brief | Status | Lean |
|---|---|---|
| (0) reuse E1Lean's definitions; same detector | done | imports; `check_det.py`: E1Lean's 100 x_j ≡ params_exp010.json set ν = 999/10000 (in order, exact rationals), r_P = 7707/10000 ≡ json; b3c_arb.py selects that set, reads its r_P, same index shift T_{2j+1}, j ≥ 1 (and its κ(3/25, 177/20) reproduces E1's κ₁ to 30 digits) |
| (1) architecture: E1′ from H1′–H4′ | proved | `theoremE1P`; variants `theoremE1P_gauss` (H1′ → Gaussian moment + compact windows), `theoremE1P_parity` (H3′ → SL51Inputs), `theoremE1P_fine` (both); non-vacuity `inputs_consistentP`, `gauss_consistent` |
| compact-window arithmetic | proved | `eta_edge` (η = (θ − 371/750)/2 ⟹ θ′ = θ − η > 371/750 ⟹ 121/1000 < (3θ′ − 1)/4), `eta_step` (every ν > 0), `H1P_of_gauss` |
| (2) algebra | proved | `nu_edge_iff` (121/1000 < (3θ − 1)/4 ⟺ θ > 371/750), `nu_edge_eq` (equality at 371/750), `exponent_edge` (via `exponent_iff`); h(θ; κ₂) ≥ h(371/750; κ₂) by E1Lean's `h_mono_theta` (already θ-free: 0 < θ₁ ≤ θ₂ ≤ 1; nothing to widen) |
| (3) κ₂ − k*(371/750) > 0, h(371/750; κ₂) > 3.35·10⁻⁵ from H4′ | proved, two-sided | `kappa2_gt_kstar`, `e1p_enclosures`: c(371/750) ∈ [−0.1857718405, −0.18577184], κ₂ ∈ [0.0850387602, 0.0850387611], k*(371/750) ∈ [0.0849914142, 0.0849914146], h(371/750; κ₂) ∈ (0.000033545, 0.000033547], κ₂ − k*(371/750) ∈ (0.000047345, 0.000047347] |
| (4) detector threshold | proved | `detector_threshold`: k*(0.4946366) < κ₂ < k*(0.4946365); h(θ; κ₂) > 0 on [0.4946366, 1], ≤ 0 on (0, 0.4946365]; 0.4946366 < 371/750 (moment edge binds) |
| (4) H4′ from exact rationals | not done | c's closed form A·e^{2R} + B needs c(P, Q, R, ν) defined (double integral; measure theory is over the memory cap) and e to ~10⁻²² (Fns has e to 10⁻¹⁰): not light |

## Named hypotheses and their sources
| Name | Statement (Lean: Arch.lean, Window.lean, E1Lean Parity.lean) | Source | Kind |
|---|---|---|---|
| H1′ (`H1P`) | ∀ θ ∈ (1/3, 1), 0 < ν < (3θ − 1)/4: SL Thm 1.2's conclusion for every admissible (R, P, Q) | REPORT 'Moment form' (E1-M′); FULLPROOF §§1–5 (Lemma B, Lemma R, §4 branches, §5 windows); reviewer passes 1–4 + referee CORRECT | analytic, new |
| H1G | ∀ θ′ ∈ (1/3, 1), 0 < ν < (3θ′ − 1)/4: Gaussian-weighted moment G θ′ ν R P Q (width T^θ′, uniform in centres T′ ≍ T), every admissible (R, P, Q) | FULLPROOF §§2–4 and §5 l. 108–109 ("the only input is the Gaussian moment at width s0 = T^θ′ … whenever ν < (3θ′ − 1)/4 (θ′ > 1/3)") | analytic, new (replaces H1′) |
| HCW | ∀ 0 < η < θ < 1, ν > 0, admissible (R, P, Q): G (θ − η) ν R P Q ⟹ SL Thm 1.2's conclusion at window T^θ | SL l. 730–758 (profile-free compact windows); FULLPROOF §5 l. 106–111 ("θ − η > 1/2" of SL l. 731 not needed) | analytic, attributed (replaces H1′) |
| H2′ (`H2P`) | θ ∈ (0, 1), 0 < ν < 1/2, (R, P, Q) admissible, Q(x) + Q(1 − x) ≡ C ≠ 0, Thm 1.2 at (θ, ν) ⟹ liminf O/N ≥ κ(P, Q, R, ν) | SL Thm 1.3 (Zenodo 23134787 v0.03, l. 125–129); θ-free per FULLPROOF §5 l. 112–113 (Lemma 4.2/Thm 1.3, l. 823–928) | analytic, attributed |
| H3′ (`H3P`) | θ ∈ (0, 1), 0 ≤ k < 1, liminf O/N ≥ k ⟹ liminf S/N ≥ h(θ; k) | SL Thm 5.1 (l. 844–845; proof l. 854–869); θ-free per FULLPROOF §5 l. 113–116 (Wang 2609.07918 Thm local-pair, TeX l. 231–232, "Let 0 < λ < θ < 1 be fixed") | analytic, attributed |
| SL51Inputs on (0, 1) (replaces H3′) | E1Lean's `SL51Inputs F I θ` (counts, parity product, Wang's pair-sum limit, optimum) for every θ ∈ (0, 1) | as E1Lean's LEAN.md (Riemann–von Mangoldt; SL ref [4] Prop 2.5; Wang Thm 2.2; cosine optimum) | finite + analytic |
| H4′ (`H4P`) | 3138.997965019943478277 ≤ c(P₁, Q₁, 44/5, 121/1000) ≤ 3138.997965019943478278 | code/b3c_arb.py, logs/b3c_arb.log (python-flint Arb, 4000 bits): c = [3138.99796501994347827760614809 ± 1.36e−27] ⊂ the interval | certificate |
| `Fns` fields | unchanged from E1Lean (√, sin, cos, exp, log facts; `rpow` symbol) | Mathlib theorems (E1Lean LEAN.md) | Mathlib |

## Main statements (verbatim; `variable (F : Fns) (I : Inputs)`; proofs in the files)
```lean
theorem nu_edge_iff (θ : ℝ) : (121 / 1000 : ℝ) < (3 * θ - 1) / 4 ↔ 371 / 750 < θ
theorem nu_edge_eq : (3 * (371 / 750 : ℝ) - 1) / 4 = 121 / 1000
theorem exponent_edge (θ : ℝ) : 1 / 2 + 2 * (121 / 1000 : ℝ) - 3 * θ / 2 < 0 ↔ 371 / 750 < θ
theorem theoremE1P (h1 : H1P I) (h2 : H2P F I) (h3 : H3P F I) (h4 : H4P I) (θ : ℝ) (hθ1 : 371 / 750 < θ)
    (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₂ F I)) ∧ h F (371 / 750) (κ₂ F I) ≤ h F θ (κ₂ F I) ∧
      (0.0000335 : ℝ) < h F (371 / 750) (κ₂ F I)
theorem kappa2_gt_kstar (h4 : H4P I) {θ : ℝ} (hθ1 : 371 / 750 ≤ θ) (hθ2 : θ ≤ 1) :
    0 < κ₂ F I - kstar F (371 / 750) ∧ kstar F θ < κ₂ F I
theorem e1p_enclosures (h4 : H4P I) :
    (-1857718405 / 10000000000 : ℝ) ≤ cW F (371 / 750) ∧ cW F (371 / 750) ≤ (-18577184 / 100000000 : ℝ) ∧
    (850387602 / 10000000000 : ℝ) ≤ κ₂ F I ∧ κ₂ F I ≤ 850387611 / 10000000000 ∧
    (849914142 / 10000000000 : ℝ) ≤ kstar F (371 / 750) ∧ kstar F (371 / 750) ≤ 849914146 / 10000000000 ∧
    (0.000033545 : ℝ) < h F (371 / 750) (κ₂ F I) ∧ h F (371 / 750) (κ₂ F I) ≤ 0.000033547 ∧
    (0.000047345 : ℝ) < κ₂ F I - kstar F (371 / 750) ∧ κ₂ F I - kstar F (371 / 750) ≤ 0.000047347
theorem inputs_consistentP : H1P trivialInputsP ∧ H2P F trivialInputsP ∧ H3P F trivialInputsP ∧ H4P trivialInputsP
theorem eta_edge {θ : ℝ} (hθ : 371 / 750 < θ) :
    0 < (θ - 371 / 750) / 2 ∧ (θ - 371 / 750) / 2 < θ ∧ 371 / 750 < θ - (θ - 371 / 750) / 2 ∧
      (121 / 1000 : ℝ) < (3 * (θ - (θ - 371 / 750) / 2) - 1) / 4
theorem eta_step {θ ν : ℝ} (hν0 : 0 < ν) (hν : ν < (3 * θ - 1) / 4) :
    0 < (θ - (4 * ν + 1) / 3) / 2 ∧ (θ - (4 * ν + 1) / 3) / 2 < θ ∧
      1 / 3 < θ - (θ - (4 * ν + 1) / 3) / 2 ∧ ν < (3 * (θ - (θ - (4 * ν + 1) / 3) / 2) - 1) / 4
theorem H1P_of_gauss (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) (hG : H1G G) (hW : HCW I G) : H1P I
theorem theoremE1P_gauss (F : Fns) (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) (hG : H1G G) (hW : HCW I G)
    (h2 : H2P F I) (h3 : H3P F I) (h4 : H4P I) (θ : ℝ) (hθ1 : 371 / 750 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₂ F I)) ∧ h F (371 / 750) (κ₂ F I) ≤ h F θ (κ₂ F I) ∧
      (0.0000335 : ℝ) < h F (371 / 750) (κ₂ F I)
theorem gauss_consistent : H1G (fun _ _ _ _ _ => True) ∧ HCW trivialInputsP (fun _ _ _ _ _ => True)
theorem H3P_of_inputs (I : Inputs) (hin : ∀ θ : ℝ, 0 < θ → θ < 1 → SL51Inputs F I θ) : H3P F I
theorem theoremE1P_parity (I : Inputs) (h1 : H1P I) (h2 : H2P F I) (hin : ∀ θ : ℝ, 0 < θ → θ < 1 → SL51Inputs F I θ)
    (h4 : H4P I) (θ : ℝ) (hθ1 : 371 / 750 < θ) (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₂ F I)) ∧ h F (371 / 750) (κ₂ F I) ≤ h F θ (κ₂ F I) ∧
      (0.0000335 : ℝ) < h F (371 / 750) (κ₂ F I)
theorem theoremE1P_fine (I : Inputs) (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) (hG : H1G G) (hW : HCW I G)
    (h2 : H2P F I) (hin : ∀ θ : ℝ, 0 < θ → θ < 1 → SL51Inputs F I θ) (h4 : H4P I) (θ : ℝ) (hθ1 : 371 / 750 < θ)
    (hθ2 : θ < 1) :
    LiminfGe (ratio F I.S I.N θ) (h F θ (κ₂ F I)) ∧ h F (371 / 750) (κ₂ F I) ≤ h F θ (κ₂ F I) ∧
      (0.0000335 : ℝ) < h F (371 / 750) (κ₂ F I)
theorem cW_detp_ge : (-185884824698 / 1000000000000 : ℝ) ≤ cW F (4946366 / 10000000)
theorem cW_detm_le : cW F (4946365 / 10000000) ≤ (-185885200141 / 1000000000000 : ℝ)
theorem kstar_anti {θ₁ θ₂ : ℝ} (h1 : 0 < θ₁) (h12 : θ₁ ≤ θ₂) (h2 : θ₂ ≤ 1) : kstar F θ₂ ≤ kstar F θ₁
theorem kstar_det_bounds : kstar F (4946366 / 10000000) ≤ 850387096 / 10000000000 ∧
    (850388666 / 10000000000 : ℝ) ≤ kstar F (4946365 / 10000000)
theorem detector_threshold (h4 : H4P I) :
    kstar F (4946366 / 10000000) < κ₂ F I ∧ κ₂ F I < kstar F (4946365 / 10000000) ∧
    (∀ θ : ℝ, 4946366 / 10000000 ≤ θ → θ ≤ 1 → 0 < h F θ (κ₂ F I)) ∧
    (∀ θ : ℝ, 0 < θ → θ ≤ 4946365 / 10000000 → h F θ (κ₂ F I) ≤ 0) ∧ (4946366 / 10000000 : ℝ) < 371 / 750
theorem cW_edge_ge : (-1857718405 / 10000000000 : ℝ) ≤ cW F (371 / 750)
theorem cW_edge_le : cW F (371 / 750) ≤ (-18577184 / 100000000 : ℝ)
theorem kappa2_bounds {c : ℝ} (hlo : (3138997965019943478277 / 1000000000000000000 : ℝ) ≤ c) (hhi : c ≤ (3138997965019943478278 / 1000000000000000000 : ℝ)) :
    (850387602 / 10000000000 : ℝ) ≤ kappa F c (44 / 5) ∧ kappa F c (44 / 5) < 1
theorem kappa2_upper {c : ℝ} (hlo : (3138997965019943478277 / 1000000000000000000 : ℝ) ≤ c) : kappa F c (44 / 5) ≤ 850387611 / 10000000000
theorem h_edge_gt {κ : ℝ} (hk : (850387602 / 10000000000 : ℝ) ≤ κ) (hk1 : κ < 1) :
    (0.000033545 : ℝ) < h F (371 / 750) κ
theorem h_edge_le {κ : ℝ} (hk : κ ≤ 850387611 / 10000000000) : h F (371 / 750) κ ≤ 0.000033547
theorem kstar_edge_bounds : (849914142 / 10000000000 : ℝ) ≤ kstar F (371 / 750) ∧
    kstar F (371 / 750) ≤ 849914146 / 10000000000
theorem sq_lo {y l l' : ℝ} (h : l ≤ F.exp y) (h0 : 0 ≤ l) (e : l' ≤ l * l) : l' ≤ F.exp (2 * y)
theorem sq_up {y u u' : ℝ} (h : F.exp y ≤ u) (e : u * u ≤ u') : F.exp (2 * y) ≤ u'
theorem exp_Lhi2 : (3138997965019943478278 / 1000000000000000000 : ℝ) ≤ F.exp (8 + 258294547 / 5000000000)
theorem exp_Llo2 : F.exp (8 + 51658903 / 1000000000) ≤ (3138997965019943478277 / 1000000000000000000 : ℝ)
theorem base_sin {a : ℝ} (ha : 0 < a) (ha1 : a ≤ 1) :
    a - a ^ 3 / 6 ≤ F.sin a ∧ F.sin a ≤ a - a ^ 3 / 6 + a ^ 5 / 100
theorem cos_two {b sl sh cl' ch' : ℝ} (hs : sl ≤ F.sin b ∧ F.sin b ≤ sh) (h0 : 0 ≤ sl)
    (e3 : cl' ≤ 1 - 2 * sh ^ 2) (e4 : 1 - 2 * sl ^ 2 ≤ ch') : cl' ≤ F.cos (2 * b) ∧ F.cos (2 * b) ≤ ch'
theorem trig_xlo2 : (2677290321021843573842093297 / 7812500000000000000000000000 : ℝ) ≤ F.sin (13991286177 / 40000000000) ∧ F.cos (13991286177 / 40000000000) ≤ (939447389341493878758616448191 / 1000000000000000000000000000000 : ℝ)
theorem trig_xhi2 : (18788947786460058914970169321 / 20000000000000000000000000000 : ℝ) ≤ F.cos (87445538607 / 250000000000) ∧ F.sin (87445538607 / 250000000000) ≤ (13707726445781940956107589929 / 40000000000000000000000000000 : ℝ) ∧ 0 < F.sin (87445538607 / 250000000000)
theorem trig_det_lo : (171336594019346940584443320563 / 500000000000000000000000000000 : ℝ) ≤ F.sin (349760894081 / 1000000000000) ∧ F.cos (349760894081 / 1000000000000) ≤ (469727337451832017944062175461 / 500000000000000000000000000000 : ℝ)
theorem trig_det_hi : (234863674778899578598089897037 / 250000000000000000000000000000 : ℝ) ≤ F.cos (174880411687 / 500000000000) ∧ F.sin (174880411687 / 500000000000) ≤ (171336560831794961546611225969 / 500000000000000000000000000000 : ℝ) ∧ 0 < F.sin (174880411687 / 500000000000)
```

## Definitions (verbatim; new in E1PLean; `variable (F : Fns) (I : Inputs)`; docstrings in the files)
```lean
def H1P : Prop := ∀ θ ν : ℝ, 1 / 3 < θ → θ < 1 → 0 < ν → ν < (3 * θ - 1) / 4 → MomentRange I θ ν
def H2P : Prop := ∀ (θ ν R : ℝ) (P Q : ℝ → ℝ), 0 < θ → θ < 1 → 0 < ν → ν < 1 / 2 → Adm R P Q →
  (∃ C : ℝ, C ≠ 0 ∧ ∀ x, Q x + Q (1 - x) = C) → MomentRange I θ ν →
  LiminfGe (ratio F I.O I.N θ) (kappa F (I.cdet P Q R ν) R)
def H3P : Prop := ∀ θ k : ℝ, 0 < θ → θ < 1 → 0 ≤ k → k < 1 →
  LiminfGe (ratio F I.O I.N θ) k → LiminfGe (ratio F I.S I.N θ) (h F θ k)
def H4P : Prop := (3138997965019943478277 / 1000000000000000000 : ℝ) ≤ I.cdet P₁ Q₁ (44 / 5) (121 / 1000) ∧
  I.cdet P₁ Q₁ (44 / 5) (121 / 1000) ≤ 3138997965019943478278 / 1000000000000000000
noncomputable def κ₂ : ℝ := kappa F (I.cdet P₁ Q₁ (44 / 5) (121 / 1000)) (44 / 5)
noncomputable def trivialInputsP : Inputs :=
  ⟨fun _ _ => 1, fun _ _ => 1, fun _ _ => 1, fun _ _ _ _ _ => True,
    fun _ _ _ _ => 6277995930039886956555 / 2000000000000000000⟩
def H1G (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) : Prop :=
  ∀ θ ν : ℝ, 1 / 3 < θ → θ < 1 → 0 < ν → ν < (3 * θ - 1) / 4 → ∀ R P Q, Adm R P Q → G θ ν R P Q
def HCW (G : ℝ → ℝ → ℝ → (ℝ → ℝ) → (ℝ → ℝ) → Prop) : Prop :=
  ∀ θ η ν : ℝ, 0 < η → η < θ → θ < 1 → 0 < ν → ∀ R P Q, Adm R P Q → G (θ - η) ν R P Q → I.Moment θ ν R P Q
```
Reused from E1Lean unchanged (verbatim in Riemann/Deep/SantiAudit/LEAN.md): `cot x = cos x / sin x`, `cW θ = 2 − θ/2 − cot(θ/√2)/√2`, `hc c k = (3 + k − √((1 − k)(9 − k − 8c)))/4`, `h θ k = hc (cW θ) k`, `kstar θ = 1 − 2/(2 − cW θ)`, `kappa c R = 1 − R⁻¹ log c`, `LiminfGe`, `ratio`, `Inputs` (N, S, O, Moment, cdet), `Adm`, `MomentRange`, `IsPoly`, `cheb`, `Qfam`, `Pfam`, `Sfun`, `xs₁`, `Q₁ = Qfam xs₁`, `P₁ = Pfam (7707/10000)`, `SL51Inputs`.

## Fidelity
REPORT = Riemann/Deep/B3Cert/REPORT.md, FP = FULLPROOF.md (l. = file lines). SL = Santibáñez-Leal: "pdftotext l." = Zenodo 23134787 v0.03 (E1Lean's numbering), "TeX l." = main.tex of 23135354 (FP's numbering).
| Lean | Source | Check |
|---|---|---|
| `theoremE1P`: for 371/750 < θ < 1, LiminfGe (S/N at T^θ) (h θ κ₂) ∧ h(371/750; κ₂) ≤ h(θ; κ₂) ∧ 0.0000335 < h(371/750; κ₂) | REPORT l. 8: "For every fixed θ ∈ (371/750, 1): liminf S(T, T^θ)/N(T, T^θ) ≥ h(θ; κ) ≥ h(371/750; κ) = 3.3546·10⁻⁵ > 0" | ✓; 3.35·10⁻⁵ is the brief's bound; `e1p_enclosures` adds h(371/750; κ₂) ∈ (0.000033545, 0.000033547] ∋ 3.354573606·10⁻⁵ (Arb, FP §8 l. 175). REPORT's '= 3.3546·10⁻⁵' and 'κ = 0.08503876039176576024' are roundings (≈) of the Arb balls h = 3.354573606042305490255·10⁻⁵ ± 1.5·10⁻³⁰, κ = 0.0850387603917657602445070238 ± 5·10⁻³² (reviewer pass 1); Lean keeps κ₂ and h symbolic and proves only inequalities |
| `κ₂ = kappa (cdet P₁ Q₁ (44/5) (121/1000)) (44/5)` | REPORT l. 10: κ = κ(P, Q, R = 44/5, ν = 121/1000), (P, Q) = SL's frozen EXP-010 detector; FP l. 115–116 | ✓; same detector as E1Lean (check_det.py, Short answer (0)) |
| S, N, h, c, k* | REPORT l. 9–10; FP l. 114 (k*(θ) = 1 − 2/(2 − c(θ))) | ✓ E1Lean's definitions, unchanged |
| H1′ | REPORT l. 11–12 (Moment form E1-M′: "holds for every θ ∈ (1/3, 1) and 0 < ν < (3θ − 1)/4") | ✓ same range; "any fixed real Q with Q(0) = 1" with R > 0, P(0) = 0, P(1) = 1 = E1Lean's `Adm`; as E1Lean's H1, the compact-window step is inside |
| H1G + HCW | FP l. 106–110: windows w_D from Gaussians of width s0 = T^{θ−η}; "The only input is the Gaussian moment at width s0 = T^θ′, θ′ := θ − η, which §§2–4 give whenever ν < (3θ′ − 1)/4 (θ′ > 1/3)", uniform in centres T′ ≍ T; SL's "θ − η > 1/2" (TeX l. 731) dropped | ✓ H1G = §§2–4's output on θ′ ∈ (1/3, 1); HCW = the profile-free step, premises 0 < η < θ < 1 only; D (fixed, from θ, η, ν) is inside HCW; `eta_edge` = FP l. 109–110's choice η = (θ − 371/750)/2; `H1P_of_gauss` uses η = (θ − (4ν + 1)/3)/2 for general ν |
| H2′ | SL Thm 1.3 (pdftotext l. 125–129); FP l. 111–112 ("(Emean) → Lemma 4.2 and Thm 1.3 (TeX l. 823–928: N(T, H) = (HL/2π)(1 + o(1)) and Backlund–Jensen O(L) hold for H = T^θ, every θ > 0)"), l. 116 ("Every input of this chain other than Thm 1.2 is θ-free on (0, 1)") | E1Lean's H2 with θ ∈ (1/2, 1) widened to (0, 1); otherwise identical: used wherever Thm 1.2's conclusion holds for all admissible weights (Q = 1, 1 + x^j included: (Emean), FP l. 111); premise ν < 1/2 = SL Lemma 2.2 (REPORT l. 16); as in E1Lean, stronger than SL's printed range (5) in ν, weaker in the extra premises |
| H3′ | SL Thm 5.1 (pdftotext l. 844–845; TeX l. 932–963); FP l. 112–114 (Wang's Thm local-pair "Let 0 < λ < θ < 1 be fixed", c(θ) defined on (0, 1), c′ > 0), l. 116 | E1Lean's H3 with θ ∈ (1/2, 1) widened to (0, 1); extra premise k < 1 (weaker) as before |
| SL51Inputs on (0, 1) | E1Lean's definition (counts, parity product per window, Wang's limit, optimum) | unchanged; `H3P_of_inputs` = E1Lean's `H3_of_inputs` with 1/2 < θ replaced by 0 < θ (`sl51` needs only 0 < θ ≤ 1, k < 1) |
| H4′ | REPORT l. 36–37 (step 5); FP §8 l. 172 ("c = 3138.99796501994347827760614809 ± 1.4·10⁻²⁷"); logs/b3c_arb.log l. 3 ("[3138.99796501994347827760614809 +/- 1.36e-27]") | ✓ ball ⊂ [3138.997965019943478277, …278] (outward rounding to 18 decimals; ball ends …60614673, …60614945) |
| `nu_edge_iff`, `nu_edge_eq` | REPORT l. 32, 38 ("(3·371/750 − 1)/4 = 121/1000 exactly"); FP l. 109 | ✓ |
| `detector_threshold` | REPORT l. 39, FP §8 l. 177–178: θ_det = 0.494636567599178 ± 7·10⁻¹⁷ < 371/750 | ✓ θ_det ∈ (0.4946365, 0.4946366) as sign change of κ₂ − k* (no IVT needed: k* antitone, `kstar_anti`) |

What is proved and what is named: proved (sorry-free): the architecture (E1′ from H1′–H4′ and its three finer variants), the compact-window arithmetic, the edge algebra, every numeric constant of FP §8 used by the Theorem (from H4′ and the `Fns` fields), the detector threshold enclosure, non-vacuity. Named analytic hypotheses: H1′ (or H1G + HCW), H2′, H3′ (or SL51Inputs on (0, 1)). Certificate hypothesis: H4′. Not formalized: Lemma B, Lemma R, the §4 branches (the content of H1G), SL §4 (H2′), Wang's theorem and the Lamzouri bound (H3′; E1Lean proves the deduction `sl51`), the double integral c(P, Q, R, ν) (opaque; H4′ certifies it), the count asymptotics.

What changed from E1Lean: θ-ranges H1 (1/2, 1) → H1′ (1/3, 1); H2, H3, SL51Inputs-range (1/2, 1) → (0, 1); H4 (R = 177/20, ν = 3/25) → H4′ (R = 44/5, ν = 121/1000); E1Lean's H1 contained the compact-window step, here it is also split off (H1G + HCW, `H1P_of_gauss`); numerics 150× tighter (below); new `detector_threshold`. No E1Lean lemma carried θ > 1/2 except H1–H3 and `H3_of_inputs` (`h_mono_theta`, `h_pos_iff`, `cW_lt`, `cW_lt_one`, `cot_sub_gt`, `sl51` need only 0 < θ ≤ 1); those four are replaced by new definitions/lemmas, E1Lean's files are untouched.

## Numerics (Lean-proved, rational; generated by Riemann/Deep/B3Cert/Lean/audit/gen_trig.py, gen_exp.py, gen_det.py, mk_lean.py, mk_num.py, mk_det.py; all bounds checked first in exact Python Fractions, check_final.py)
- Precision plan (brief's warning: h − 3.35·10⁻⁵ ≈ 4.6·10⁻⁸; ∂h/∂k = 0.709, ∂h/∂c = 0.297 at the edge): one-step Taylor at x ≈ 0.35 is useless (~10⁻⁴); E1Lean's 3 doublings with cos_bound would give ~3·10⁻⁶ in c. Used: sin at x/32 from sin_gt_sub_cube and sin_bound (relative error ≤ (x/32)⁴/120 ≈ 1.2·10⁻¹⁰), cos(x/32) = 1 − 2sin²(x/64) (`cos_two`; error ~(x/64)⁶), five doublings with E1Lean's `dbl` (cos 2a = 1 − 2sin²a), every bound rounded outward to 30 decimals.
- x = (371/750)/√2 via E1Lean's √2 ∈ [1.41421356237, 1.41421356238]: x_lo = 13991286177/(4·10¹⁰) ≤ x ≤ x_hi = 87445538607/(2.5·10¹¹); cot monotone (`cot_sub_gt`). Fractions: c(371/750) ∈ [−0.185771840418, −0.185771840061] (width 3.6·10⁻¹⁰; Arb −0.185771840139); stated rounded: [−0.1857718405, −0.18577184].
- κ₂: e⁸ from e ∈ (2.7182818283, 2.7182818286) (`exp_one_gt_d9`, `exp_one_lt_d9`, `exp_nat`); e^f ≥ 8 squarings (`sq_lo`) of 1 + y + y²/2 (y = f/2⁸); e^g ≤ 18 squarings (`sq_up`) of 1/(1 − y) (y = g/2¹⁸). log c(P₁, Q₁, 44/5, 121/1000) ∈ [8.051658903, 8.0516589094] (L_lo = log c − 5.6·10⁻⁹, L_hi = log c + 8.5·10⁻¹⁰), κ₂ ∈ [0.0850387602, 0.0850387611] (Fractions: [κ − 9.6·10⁻¹¹, κ + 6.3·10⁻¹⁰]).
- h: `sqrt_lt_of_lt_sq`, `le_sqrt_of_sq_le` at the rounded ends with E1Lean's monotonicity `hc_mono_k`, `hc_mono_c`; Fractions: h(c_lo, κ_lo) = 3.354549·10⁻⁵ > 3.3545·10⁻⁵, h(c_hi, κ_hi) = 3.354628·10⁻⁵ ≤ 3.3547·10⁻⁵; κ₂ − k*: (4.73457·10⁻⁵, 4.73468·10⁻⁵) inside the stated (0.000047345, 0.000047347].
- Threshold: c(0.4946366) ≥ −0.185884824698, c(0.4946365) ≤ −0.185885200141 (two more chains, DetTrig.lean), so k*(0.4946366) ≤ 0.0850387096 < 0.0850387602 ≤ κ₂ ≤ 0.0850387611 < 0.0850388666 ≤ k*(0.4946365) (margins 5.1·10⁻⁸, 1.06·10⁻⁷).
- Decimals labelled Arb are b3c_arb.log values, not Lean claims.

## Build
Mathlib: `rm .lake && bash .scratch/memcap.sh 290 1500 lake exe cache get` (20 Nov 05:42 UTC; 52 s, peak 246 MB, rc = 0). LEAN_PATH = `lake env printenv LEAN_PATH` (Riemann/Deep/B3Cert/Lean/audit/leanpath.txt). Imports: E1Lean's modules only (Numeric, Bounds2, Parity) plus our own; no new Mathlib import. Each module separately, one at a time, lean directly (no lake build), oleans to .scratch/b3lean/olean (deleted at the end):

    LEAN_PATH=.scratch/b3lean/olean:$(cat Riemann/Deep/B3Cert/Lean/audit/leanpath.txt) OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 1800 \
      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 \
      -o .scratch/b3lean/olean/$d/$m.olean $d/$m.lean

Scripts: Riemann/Deep/B3Cert/Lean/audit/mk.sh <dir> <Module> (one module), verify.sh (clean rebuild into a fresh olean dir, E1Lean's 9 modules then ours, order Trig, Exp, Num, Arch, Window, Parity, DetTrig, Det, Check). Final clean rebuild (20 Nov 06:05–06:15 UTC, system load ≈ 10 from other agents, background QoS), every module rc = 0:

| module | wall | footprint (time -l) | memcap peak (top MEM) |
|---|---|---|---|
| E1Lean Defs | 15 s | 243 MiB | 170 MB |
| E1Lean Detector | 15 s | 252 MiB | 204 MB |
| E1Lean AlgebraH | 19 s | 252 MiB | 233 MB |
| E1Lean AlgebraC | 41 s | 278 MiB | 264 MB |
| E1Lean Numeric | 16 s | 265 MiB | 248 MB |
| E1Lean Kappa | 20 s | 248 MiB | 188 MB |
| E1Lean Arch | 23 s | 249 MiB | 214 MB |
| E1Lean Parity | 49 s | 286 MiB | 274 MB |
| E1Lean Bounds2 | 38 s | 273 MiB | 251 MB |
| Trig | 17 s | 258 MiB | 189 MB |
| Exp | 74 s | 254 MiB | 240 MB |
| Num | 126 s | 273 MiB | 248 MB |
| Arch | 45 s | 255 MiB | 233 MB |
| Window | 27 s | 247 MiB | 207 MB |
| Parity | 9 s | 239 MiB | 0 MB (exited between samples) |
| DetTrig | 9 s | 255 MiB | 235 MB |
| Det | 30 s | 267 MiB | 249 MB |
| Check | 18 s | 238 MiB | 144 MB |

All compiles of the session (38 runs, peaks.log + verify.log): footprint ≤ 286 MiB, top MEM ≤ 274 MB (both E1Lean's Parity); E1PLean files alone: footprint ≤ 273 MiB, top MEM ≤ 262 MB; one failed compile (Num, a reduced-literal mismatch with Exp's c_hi, 125 s, 262 MB; fixed by writing c_hi unreduced). Scripts (Python, exact Fractions, single-threaded) each ended before memcap's first sample (peak 0 MB). Every run under memcap 290 MB; all below the 300 MB cap.

## Static scan
`rg -n -w 'sorry|admit|native_decide|skipKernelTC|decide' Riemann/Deep/B3Cert/Lean/*.lean` and `rg -n '^\s*axiom\s' Riemann/Deep/B3Cert/Lean/*.lean`: no match (rc = 1). No errors or warnings in any compile log.

## #print axioms (Check.lean, 37 lines, all identical)

    'E1PLean.theoremE1P' depends on axioms: [propext, Classical.choice, Quot.sound]
    (same line for theoremE1P_gauss, theoremE1P_parity, theoremE1P_fine, kappa2_gt_kstar, e1p_enclosures, detector_threshold, inputs_consistentP, gauss_consistent, H1P_of_gauss, H3P_of_inputs, eta_edge, eta_step, nu_edge_iff, nu_edge_eq, exponent_edge, cW_edge_ge, cW_edge_le, kappa2_bounds, kappa2_upper, h_edge_gt, h_edge_le, kstar_edge_bounds, exp_Lhi2, exp_Llo2, sq_lo, sq_up, base_sin, cos_two, trig_xlo2, trig_xhi2, trig_det_lo, trig_det_hi, cW_detp_ge, cW_detm_le, kstar_anti, kstar_det_bounds)

## Review (statement fidelity, verdict-only reviewer subagents; prompts and answers in Riemann/Deep/B3Cert/Lean/audit/review_prompt.txt, review1_out.md, review2_out.md, copied below in substance)
- Pass 1 (REPORT l. 1–17, 36–39, FP l. 105–117, 179–181, b3c_arb.log l. 3, E1Lean's definitions, all new definitions with docstrings, main statements, disclosed choices (a)–(e)): **ISSUES** — one item: REPORT's displayed decimals (κ = 0.08503876039176576024, h(371/750; κ) = 3.3546·10⁻⁵) are roundings of the Arb balls, stated as equalities; "Lean correctly retains the exact symbolic κ₂ and h … Do not strengthen Lean"; otherwise faithful (ranges, H2P's disclosed extension, H3P on (0, 1), H4P ⊇ ball, HCW, η arithmetic, conclusion). Repair: fidelity row above (REPORT.md is outside this task's edit scope: orchestrator to write ≈ there).
- Pass 2 (same paste + pass 1 + reply): **FAITHFUL** — "No Lean statement correction required under the stated disclosures."

## Disk and .lake
.lake was the placeholder file; deleted and filled by `bash .scratch/memcap.sh 290 1500 lake exe cache get` (20 Nov 05:42 UTC: cloned the packages, decompressed Mathlib's cache; 52 s, memcap peak 246 MB, rc = 0). Left as it is: a real directory with Mathlib and its dependencies (~6.5 GB); the orchestrator restores the placeholder. No `lake build` was run. Build output .scratch/b3lean/olean, olean_v deleted at the end; generators, logs and review files kept in Riemann/Deep/B3Cert/Lean/audit/.

## FORMALIZATIONS.md entry
- Theorem E1′ (simple critical zeros in every window (T, T + T^θ], θ > 371/750, below √T; N8) formalized as an architecture: `E1PLean.theoremE1P` — liminf S/N ≥ h(θ; κ₂) ≥ h(371/750; κ₂) > 3.35·10⁻⁵ for every θ ∈ (371/750, 1) — from named inputs H1′ (E1-M′, θ ∈ (1/3, 1); also split into the Gaussian moment + SL's compact windows), H2′ (SL Thm 1.3), H3′ (SL Thm 5.1, also from its inputs), H4′ (Arb certificate); edge algebra, all FULLPROOF §8 constants (two-sided, ~10⁻⁹) and the detector threshold θ_det ∈ (0.4946365, 0.4946366) proved; sorry-free, standard axioms, E1Lean's `Fns`; reviewer FAITHFUL (pass 2). [LEAN.md](Riemann/Deep/B3Cert/LEAN.md)
