PARTIAL: Lemma 1 and Lemma 2 (Combos NOTE §1) are proved in Lean from (F1)+(F2) and (F3) alone; for Theorem 1′(a) (RCApps2 Prop 1, exponent e(σ) = 2(1−σ)/(3−2σ), which supersedes Theorem 1's 2σ(1−σ)) the per-window count, the class split, every exponent balance and both 𝔔-cases (`large_case`, `small_case`) are proved from (F3), (F4), Lemma 5 and the two class bounds as named hypotheses. Left: the 25-line final combination `thm1'a` (written, not compiled: the 290 MB watchdog killed five tries of two layouts at 302–322 MB), Lemma 3, proofs of Lemmas 4–5 and of the class I/II analytic chains, and Theorem 1(a) itself.

# Combos ring class zero density in Lean 4 + Mathlib (rclean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile). Files: `Riemann/Deep/Combos/Lean/`, 26 modules, 1,500 lines. No `sorry`, `admit`, `axiom` or `native_decide` (checked with rg). Namespace `RingClass`.

## What is proved

| NOTE item | Lean | Inputs |
|---|---|---|
| Lemma 1 (l.67–75) | `lemma1` (LatticeA–C, Lemma1) | none |
| Lemma 2 (l.76–89), constant 30/ε | `lemma2` (Ortho, SieveA–C, ZetaA–B, Lemma2) | `F12` = (F1)+(F2), `F3` |
| orthogonality step (l.82–83) | `bessel`: Σ_{χ∈𝒞}\|Σ_g χ(g)v(g)\|² ≤ \|G\|Σ_g\|v(g)\|² | Mathlib `AddChar.sum_eq_ite` |
| ζ(1+2ε) ≤ 1/ε (l.89) | `sum_rpow_neg_le`, `zeta_eps` | none (telescoping) |
| per-window count, class split (NOTE l.162–165; RCApps2 l.33–39, 57) | `ZeroData.N_le`, `ZeroData.Xs_card_le` (Count) | none |
| exponent bookkeeping, Thm 1′(a) (RCApps2 l.35–61) | `balance`, `V_le_Y`, `classII_exp`, `classI_exp`, `window_log`, `final_exp`, `large_arith`, `small_arith` (Expo1p*, Thm1pL) | none |
| Thm 1′(a), 𝔔 ≥ 𝔔₀ and 𝔔 < 𝔔₀ | `large_case`, `small_case` (Thm1pLg, Thm1pSm) | `RCZeros.Hyp` |
| Thm 1′(a) | `thm1'a`: draft `Lean/audit/Thm1pD.lean.draft`, `Thm1pMain.lean.draft`, not compiled | `RCZeros.Hyp` |

## Main statements

```lean
theorem lemma1 (a b c : ℤ) (ha : 0 < a) (hb : |b| ≤ a) (hac : a ≤ c) (N : ℝ) (hN : 1 ≤ N) :
    (PQ a b c N : ℝ) ≤ 2 + 8 * N / Real.sqrt |((b : ℝ) ^ 2 - 4 * a * c)|
      + 4 * Real.sqrt N / Real.sqrt (Real.sqrt |((b : ℝ) ^ 2 - 4 * a * c)|)
-- PQ a b c N = ncard {(x,y) : ℤ² | gcd x y = 1 ∧ a x² + b x y + c y² ≤ N}  (finite: PQset_finite)

structure RCData (G ι) [AddCommGroup G] where   -- G = Cl(O), ι = O_K-ideals prime to f
  D : ℤ; D_neg : D < 0; nrm : ι → ℕ; cls : ι → G; cont : ι → ℕ; prim : ι → ι
  cont_pos : ∀ A, 0 < cont A; cont_prim : ∀ A, cont (prim A) = 1; cls_prim : ∀ A, cls (prim A) = cls A
  nrm_eq : ∀ A, nrm A = cont A ^ 2 * nrm (prim A); decomp_inj : ∀ A B, prim A = prim B → cont A = cont B → A = B
def F12 d := ∀ C : G, ∃ a b c (ideal : ℤ × ℤ → ι), 0 < a ∧ |b| ≤ a ∧ a ≤ c ∧ b^2 − 4ac = d.D ∧
  ∀ B, d.cont B = 1 → d.cls B = C → ∃ x y, gcd x y = 1 ∧ Qf a (−b) c x y = d.nrm B ∧ ideal (x, y) = B
def F3 d := (card G : ℝ) ≤ 3 * |D| ^ (1/2) * log (3|D|) ^ 2
theorem lemma2 (d : RCData G ι) (hF12 : F12 d) (hF3 : F3 d) (ε) (hε : 0 < ε) (hε2 : ε ≤ 1/2) (N) (hN : 1 ≤ N)
    (S : Finset ι) (hS : ∀ A ∈ S, nrm A ≤ N) (c : ι → ℂ) (𝒞 : Finset (AddChar G ℂ)) :
    ∑ χ ∈ 𝒞, normSq (∑ A ∈ S, c A * χ (cls A)) ≤
      30 / ε * log (3|D|) ^ 2 * (N + |D| ^ (1/2)) * ∑ A ∈ S, (cont A : ℝ) ^ (1 + 2ε) * normSq (c A)

structure ZeroData (X) where one : X; win : X → ℝ → Finset ℂ; mult : X → ℂ → ℕ; win_im : ρ ∈ win χ t → |Im ρ − t| ≤ 1
def N z σ t := ∑ χ, ∑ ρ ∈ (win χ t).filter (σ ≤ re ·), mult χ ρ ;  def Xs z σ t := {χ | χ ≠ one ∧ ∃ ρ ∈ win χ t, σ ≤ re ρ}
lemma N_le (hB : ∀ χ, ∑ ρ ∈ win χ t, mult χ ρ ≤ B) : N σ t ≤ B * (#Xs σ t + 1)
lemma Xs_card_le (hsplit : ∀ χ ≠ one, ∀ ρ ∈ win χ t, σ ≤ re ρ → PI χ ρ ∨ PII χ ρ)
    (hI : ∀ S ρ, (∀ χ ∈ S, χ ≠ one ∧ ρ χ ∈ win χ t ∧ σ ≤ re (ρ χ) ∧ PI χ (ρ χ)) → #S ≤ BI) (hII : same, PII, BII) :
    #Xs σ t ≤ BI + BII

structure RCZeros (X) extends ZeroData X where
  absD : ℝ; absD_ge : 3 ≤ absD; A : X → ℝ; A_ge : 1 ≤ A χ; A_le : A χ ≤ absD
  detI : ℝ → ℝ → X → ℂ → ℝ    -- |Σ_{N𝔞>V} a_𝔞χ(𝔞)N𝔞^{−ρ}e^{−N𝔞/Y}|
  detII : ℝ → X → ℝ → ℝ       -- ∫|L(½+i(γ+u),χ)M_V(½+i(γ+u),χ)|e^{−|u|}du
structure RCZeros.Hyp d c₀ c₁ CI CII Q0 : Prop where     -- 𝔔 = 9|D|, V = 𝔔^{1/2}, T = |t| + 3
  F3 : card X ≤ 3 absD^{1/2} log²(3 absD)
  F4 : ∀ χ t, ∑ ρ ∈ win χ t, mult χ ρ ≤ c₀ log(A χ (|t|+3)²);  c₀_nonneg : 0 ≤ c₀
  L5 : ∀ σ ∈ [0.6,1], ∀ V Y, 2 ≤ V → V ≤ Y → ∀ χ ≠ one, ∀ t, ∀ ρ ∈ win χ t, σ ≤ re ρ →
         1/4 ≤ detI V Y χ ρ ∨ c₁ Y^{σ−1/2} ≤ detII V χ (Im ρ)
  classII : ∀ ε₀ ∈ (0,1/10], σ ∈ [0.6,1], t, Y ≥ V, S, ρ with class-II zeros in the window:
         #S ≤ CII ε₀ · Y^{1/2−σ} 𝔔^{1/2+ε₀} T^{1/2+ε₀}
  classI : ∀ ε₀ ∈ (0,1/10], Q0 ε₀ ≤ 𝔔, σ ∈ [0.6,1], t, Y ≥ V, S, ρ with class-I zeros:
         #S ≤ CI ε₀ · (𝔔T)^{ε₀} (2Y log²Y)^{2−2σ+2ε₀}
  CI_nonneg, CII_nonneg
lemma large_case (h : d.Hyp c₀ c₁ CI CII Q0) (h0 : 0.6 ≤ σ) (h1 : σ ≤ 1) (he : 0 < ε₀) (he1 : ε₀ ≤ 1/10)
    (hQ : Q0 ε₀ ≤ 9 absD) : N σ t ≤ c₀ (2/ε₀) (CI ε₀ · 8/ε₀² + CII ε₀ + 1) · (9 absD T)^{ee σ + 5ε₀}
lemma small_case (h) (h0) (h1) (he) (hQ : 9 absD < Q0 ε₀) :
    N σ t ≤ c₀ (2/ε₀) |3 (Q0 ε₀/9)^{1/2} log²(3 Q0 ε₀/9) + 1| · (9 absD T)^{ee σ + 5ε₀}
lemma final_exp : (9 absD T)^{e + ε/2} ≤ 81 (absD T)^{e + ε}     -- 1 ≤ absD, T; e ≤ 1; 0 < ε ≤ 1
-- draft, not compiled:
theorem thm1'a (c₀ c₁) (CI CII Q0) (ε) (hε : 0 < ε) (hε1 : ε ≤ 1) : ∃ C, ∀ (X : Type) [Fintype X] [DecidableEq X]
    (d : RCZeros X), d.Hyp c₀ c₁ CI CII Q0 → ∀ σ t, 0.6 ≤ σ → σ ≤ 1 → N σ t ≤ C * (absD * (|t|+3)) ^ (ee σ + ε)
```

The draft proof: C := 81·c₀(2/ε₀)(CI·8/ε₀² + CII + 1 + |K_S|) at ε₀ = ε/10, case split on 𝔔₀(ε₀) ≤ 9|D|, `large_case`/`small_case`, `fin_step`, `final_exp` (5ε₀ = ε/2). All its ingredients are compiled; only this glue is not.

## Named hypotheses and sources

| Name | Lean | Source |
|---|---|---|
| (F1)+(F2) | `F12` | NOTE l.44–51; Cox Prop 7.20, 7.22, Thm 7.7 (𝔟 ∩ O = α𝔞_Q⁻¹, N𝔟 = Q(x,−y), 𝔟 ⊂ mO ⟺ m \| gcd(x,y)) |
| (F3) | `F3`, `Hyp.F3` | NOTE l.52–53 (Cox Thm 7.24) |
| (F4) | `Hyp.F4` with `A_le`, `A_ge` | NOTE l.54–55, l.35 (IK Prop 5.7; A(χ*) ≤ \|D\|) |
| Lemma 5 | `Hyp.L5` | NOTE l.134–149 (mollifier length V in place of X) |
| class II bound | `Hyp.classII` | RCApps2 l.40–48 (Lemma 4 at mollifier length 1, Lemma 2, Lemma 3, Cauchy–Schwarz over χ) |
| class I bound | `Hyp.classI` | RCApps2 l.50–55 (tail beyond Z, dyadic blocks, removal of β, Lemma 2, Sobolev, Lemma 3; 𝔔 ≥ 𝔔₀(ε₀)) |

Lemma 4 enters only through `classII`. The integrals of Lemma 5 (II) and of the class II step are abstracted into `detII`: MeasureTheory.Integral.Bochner (573 MB of oleans) does not fit the memory cap.

## Build

Each module compiled alone, `lean -j1 -o` to `.scratch/rclean/olean` under `bash .scratch/memcap.sh 290 1500` (`Lean/audit/mk.sh`; LEAN_PATH from `lake env`, computed once). Peak top MEM of the successful builds (MB): LatticeA 122, LatticeB 277, LatticeC 277, Lemma1 124, Ortho 253, SieveA 238, SieveB 119, SieveC 120, ZetaA 120, ZetaB 120, Lemma2 120, Count 120, Expo1pA 121, Expo1pB 120, Expo1pC1 266, Expo1pC2 120, Expo1pC3 286, Expo1pC 120, Expo1pD 120, Thm1pI 273, Thm1pL 222, Thm1pW 190, Thm1pLg 274, Thm1pSm 244. Every successful build: exit 0, no errors; one linter warning (unused binder he0 in final_exp, Expo1pD). Killed attempts (watchdog, peaks 294–327 MB) are in `Lean/audit/peaks.log`; each was followed by a split or lighter tactics (`linarith only`, no `nlinarith` with big contexts), never by an unchanged rerun except `thm1'a` (five tries, two layouts) and the final restricted-Y recompile (killed, reverted, see Reviewer).

`#print axioms` (CheckL12.lean, CheckBK.lean, Lean/audit/CkLg.lean, CkSm.lean): `lemma1`, `bessel`, `sieve_core`, `lemma2`, `ZeroData.N_le`, `ZeroData.Xs_card_le`, `balance`, `V_le_Y`, `classII_exp`, `classI_exp`, `window_log`, `final_exp`, `large_arith`, `small_arith`, `large_case`, `RCZeros.window`, `fin_step`, `small_case`: each "depends on axioms: [propext, Classical.choice, Quot.sound]".

## Fidelity

- `lemma1`: l.69–75. P_Q(N) is the ncard of the honest set (proved finite). Hypothesis |b| ≤ a ≤ c, 0 < a instead of "reduced" (weaker; only 3a² ≤ \|D\| is used). Roots by `Real.sqrt`: \|D\|^{−1/2} = 1/√\|D\|, N^{1/2} = √N, \|D\|^{−1/4} = 1/√√\|D\|.
- `RCData`: l.32–39. ι is "O_K-ideals prime to f"; cls = [𝔞] via (F1); cont = m(𝔞), prim = 𝔞₀; `cls_prim` is l.38–39; `nrm_eq`, `decomp_inj` encode 𝔞 = m(𝔞)𝔞₀ (unique). Characters: Mathlib `AddChar G ℂ` (G written additively).
- `F12`: what l.86–87 uses from (F1)+(F2): each primitive 𝔟 of class C is α𝔞_Q⁻¹ with gcd(x,y) = 1, N𝔟 = Q(x,−y); Q = (a,b,c), Q(x,−y) = Qf a (−b) c x y.
- `lemma2`: l.77–79, for every finite set 𝒞 of characters (𝒞 = Ĝ is the NOTE; Mathlib's Fintype instance on AddChar lives in a module too heavy for the cap). Orthogonality is used as Bessel's inequality (≤), which is all the proof needs; \|z\|² = Complex.normSq z. Support: Finset S with N𝔞 ≤ N.
- `ZeroData`: win χ t = distinct zeros of L(s,χ*) with \|γ − t\| ≤ 1, mult = multiplicity; N is l.40.
- `RCZeros.Hyp`: constants c₀, c₁, CI, CII, Q0 are arguments, fixed before D (∃ C before ∀ X d in `thm1'a`): C(ε) is uniform in D. c₁ > 0 is not required (weaker). `Hyp.F3` uses X = Ĝ (\|Ĝ\| = h).
- Parameters: RCApps2 l.35–36 exactly: ε₀ = ε/10, 𝔔 = \|D\|(M₀+3)² with M₀ = 0, V = 𝔔^{1/2}, Y = (𝔔T)^{1/(3−2σ)}, Z = 2Y log²Y; window bound log(\|D\|T²) instead of log(16𝔔T²) (smaller). The σ-range 0.6 ≤ σ ≤ 1 enters through L5 (Lemma 5), V ≤ Y (3 − 2σ ≤ 2) and p = 2 − 2σ + 2ε₀ ∈ [0, 1].
- `thm1'a` is stated for ε ≤ 1 (larger ε follow, since \|D\|T ≥ 1).

What the interface abstracts away: the ring class group, ideals and their arithmetic (only the RCData fields), Hecke L-functions, their zeros and multiplicities (only `win`, `mult`), the analytic conductor (only `A`, `A_le`), and the two Lemma 5 quantities (`detI`, `detII`: arbitrary reals). Nothing ties `win` to actual zeros; the theorem is a reduction: any data satisfying the named inputs obeys the bound.

## Reviewer

One reviewer run (statements and NOTE/RCApps2 text inline; full verdict `Lean/audit/review.md`): (A) Lemma 1 FAITHFUL; (B) Lemma 2 + F12/F3/RCData FAITHFUL ("F12 is a consequence, not a stronger arithmetic assertion"); (C) counting FAITHFUL; (E) exponent lemmas FAITHFUL; (F) FAITHFUL as a conditional reduction (uncompiled draft); (D) ISSUE: `classI`/`classII` quantify over every Y ≥ V, the paper displays them at Y = (𝔔T)^{1/(3−2σ)}. OVERALL: ISSUES.
Response: for class II the paper's derivation has no Y on its right side and holds for every Y ≥ V; for class I it absorbs log⁴(3𝔔Y) into (𝔔T)^{ε₀}, so for Y ≫ (𝔔T)^{O(1)} the all-Y field is stronger than the paper. The proofs use only the chosen Y. A restricted version (fields at Y = (𝔔T)^{1/(3−2σ)} only) is written (`Lean/audit/Thm1pI.restrictedY.lean.draft`, `Thm1pLg.restrictedY.lean.draft`; Thm1pLg changes by dropping two arguments) but its rebuild was killed by the watchdog (317–319 MB), so the compiled files keep the all-Y fields. Open: recompile the restricted version, then `thm1'a`.

## .lake

The repo-root `.lake` was a real directory (shared Mathlib cache): reused, not deleted, no `lake exe cache get`, no `lake build`; our oleans went to `.scratch/rclean/olean`. Left as is; the orchestrator restores the placeholder.

## FORMALIZATIONS.md entry

- Combos ring class zero density ([discovery](discoveries/ring-class-zero-density.md)) in Lean 4 + Mathlib, PARTIAL: Lemma 1 (primitive lattice points) and Lemma 2 (weighted large sieve, constant (30/ε)log²(3\|D\|)(N+\|D\|^{1/2})) proved from (F1)+(F2), (F3) only; for Theorem 1′(a) (RCApps2, exponent 2(1−σ)/(3−2σ)) the per-window count, class split, exponent balances and both 𝔔-cases proved from (F3), (F4), Lemma 5 and the class I/II bounds as named hypotheses; final glue `thm1'a` and the restricted-Y fix not compiled (memory cap). Reviewer: FAITHFUL except class bounds quantified over all Y. [LEAN.md](Riemann/Deep/Combos/LEAN.md)
