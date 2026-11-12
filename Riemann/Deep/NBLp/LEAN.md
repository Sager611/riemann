PARTIAL: Proposition N (Theorem A (ii) ⇒ (i), PROOF.md §3) is proved in Lean from Mathlib alone (NBLp.propN, NBLp.zeroFree_of_chiInW; axioms propext, Classical.choice, Quot.sound), and Lemma M (§4) is proved; the §4–§5 toolkit (Mellin–Plancherel, Lemmas T, PW, Γ(b), H, A, P) is stated as named hypotheses. Left: Theorem A (i) ⇒ (ii) (§7 A1–A7) and the corollaries. Their Lean drafts (draft/Sec7.lean, draft/A4.lean) were never compiled: the orchestrator's 290 MB per-process cap (08:36 UTC) is lower than what the imports alone need (Mathlib.Analysis.SpecialFunctions.Pow.Real 334 MB, Mathlib.MeasureTheory.Function.LpSeminorm.Defs 343 MB, Mathlib.NumberTheory.LSeries.RiemannZeta 342 MB).

# NBLp in Lean 4 + Mathlib (nblplean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile.toml). Source: PROOF.md (Orca agent nblp, 5 Nov 2026), discovery [nyman-beurling-lp](../../../discoveries/nyman-beurling-lp.md).

## Files

| File | Lines | Status |
|---|---|---|
| Basic.lean | 354 | checked (07:57–08:11 UTC). Sections A–F of NymanBeurling/Suff.lean copied verbatim, plus isPreconnected_U from Mertens/ZeroFree.lean, so it needs Mathlib only: ζ(s) = s/(s−1) − sF(s) on Re s > 0, s ≠ 1 (zeta_eq); M(s) = ∫₀^∞ {1/x}x^{s−1}dx = 1/(1−s) + F(s) (M_eq); ζ(s) = 0 ⇒ ∫₀¹ ρ_θ x^{s−1} = 0 (integral_rho_eq_zero) |
| PropN.lean | 167 | checked. PROOF.md §0 definitions; Proposition N |
| Toolkit.lean | 179 | checked. §4–§5 named hypotheses; Lemma M proved |
| draft/Sec7.lean | 521 | NOT checked. §7 objects, [C8] facts, A1, A5, A6 (L^p part), A7 |
| draft/A4.lean | 332 | NOT checked. A4 from A3, A4(b), Lemma H, Lemma Γ(b), Plancherel |

No sorry, admit, axiom or native_decide in any file.

## Main Lean statements (checked)

~~~lean
namespace NBLp
-- PropN.lean
def rho (k : ℕ) (x : ℝ) : ℝ := Int.fract (1 / (k * x))                                  -- ρ_k
def fk (k : ℕ) (x : ℝ) : ℝ := Int.fract (1 / (k * x)) - (k : ℝ)⁻¹ * Int.fract (1 / x)   -- f_k
abbrev vol01 : Measure ℝ := volume.restrict (Ioc 0 1)
def approx (N : ℕ) (c : ℕ → ℂ) (x : ℝ) : ℂ := 1 - ∑ k ∈ Finset.Icc 2 N, c k * (fk k x : ℂ)
def ChiInW (p : ℝ) : Prop := ∀ ε > 0, ∃ N c, eLpNorm (approx N c) (ENNReal.ofReal p) vol01 < ENNReal.ofReal ε
def ZeroFree (θ : ℝ) : Prop := ∀ s : ℂ, θ < s.re → riemannZeta s ≠ 0
lemma fk_eq_rhoNB (k x) : fk k x = rhoNB (1 / k) x      -- the page's ρ_a(x) = {a/x} − a{1/x}, a = 1/k
theorem propN {p : ℝ} (hp : 1 < p) {w : ℂ} (hz : riemannZeta w = 0) (hw : 1 / p < w.re) : ¬ ChiInW p
theorem zeroFree_of_chiInW {p : ℝ} (hp : 1 < p) (h : ChiInW p) : ZeroFree (1 / p)
-- Toolkit.lean
def L2w (σ) (F : ℝ → ℂ) : ℝ := ∫ x in Ioi 0, ‖F x‖ ^ 2 * x ^ (2σ − 1)                   -- ‖F‖²_{L²_σ}
def MemL2w (σ) (F) : Prop := AEStronglyMeasurable F (volume.restrict (Ioi 0)) ∧ IntegrableOn (‖F‖² x^{2σ−1}) (Ioi 0)
def chiP (x) : ℂ := if x ≤ 1 then 1 else 0                                                -- χ₊
structure MellinToolkit (HM : ℝ → (ℝ → ℂ) → (ℝ → ℂ) → Prop) : Prop   -- HM σ F G: F ∈ L²_σ, 𝓜_σF = G a.e.
  mem : HM σ F G → MemL2w σ F
  plancherel : HM σ F G → Integrable (‖G‖²) ∧ L2w σ F = (2π)⁻¹ ∫ ‖G t‖²
  add, smul (c : ℂ), congr (F =ᵐ F', G =ᵐ G')
  l1 : MemL2w σ F → IntegrableOn (x^{σ−1} F) (Ioi 0) → HM σ F (t ↦ ∫ x in Ioi 0, x^{σ+it−1} F x)
def LemmaT HM, LemmaPW HM, LemmaGb, LemmaH, LemmaA, LemmaP : Prop      -- table below
def gam s := ½ s (s−1) Gammaℝ s;  def Gr ε s := gam (s+ε) / gam s;  def xi s := ½ s (s−1) Λ₀(s) + ½
def PN N s := ∑ k ∈ Icc 1 N, μ(k) (1 − k/N) k^{−s};  def invZeta s := if s = 1 then 0 else (riemannZeta s)⁻¹
theorem lemmaM {p r a : ℝ} (hp1 : 1 ≤ p) (hr0 : 0 < r) (hr : 1/p = 1/2 + 1/r) (ha : a < 1/p)
    {F : ℝ → ℂ} (hF : AEStronglyMeasurable F vol01) :
    eLpNorm F (ofReal p) vol01 ≤ ofReal ((1 − r(a − ½))^{−1/r}) * eLpNorm (x ↦ x^{a−½} F x) 2 vol01
~~~

## Named hypotheses and sources (Toolkit.lean)

| Name | Lean form | PROOF.md | Source |
|---|---|---|---|
| MellinToolkit HM | membership, Plancherel ‖F‖²_{L²_σ} = (2π)⁻¹∫\|G\|², linearity, a.e. congruence, the L¹ formula 𝓜_σF(t) = ∫F x^{s−1} | §4, l. 118–121 | [C7] Plancherel; true for HM σ F G := F ∈ L²_σ ∧ G = 𝓜_σF a.e. |
| LemmaT HM | (a) 0 < σ < 1, k ≥ 1: x^{σ−1}ρ_k ∈ L¹, ρ_k ∈ L²_σ, ‖ρ_k‖² = k^{−2σ}‖ρ_1‖², ‖ρ_1‖² ≤ 1/(2σ) + 1/(2−2σ), 𝓜_σρ_k = −k^{−s}ζ(s)/s; (b) σ > 0: 𝓜_σχ₊ = 1/s | l. 123–130 | [C1] (Titchmarsh 2.1.5); its integral formula is Basic.lean's M_eq + mellin_comp_mul_left |
| LemmaPW HM | Ψ holomorphic on an open U ⊇ {a ≤ Re s ≤ b}, sup_σ ∫\|Ψ(σ+it)\|² < ∞, 𝓜_bF = Ψ(b+i·) ⇒ 𝓜_aF = Ψ(a+i·) | l. 139–152 | Paley–Wiener-type shift, proved in PROOF.md |
| LemmaGb | 0 < ε ≤ 1, Re s ≥ ½, segment [s, s+ε] in {Re w < 1} ∪ {\|Im w\| ≥ 1}: \|Γr_ε(s)\| ≤ e^{5ε}(2+\|s\|)^{ε/2}; if Re(s+ε) < 1: \|Γr_ε(s)\|² ≤ e^{8ε}(2+\|t\|)^ε | l. 155–171, Lemma Γ(b) | [C3] Binet, coth series (DLMF 5.9.1, 4.36.3); Γ(a) not stated separately |
| LemmaH | all nontrivial zeros β ≤ θ*, ε > 0, Re s ≥ θ* − ε/2, s ≠ 1, Re s > 0 ⇒ \|ξ(s)\| ≤ \|ξ(s+ε)\| and (s+ε ≠ 1 ⇒ \|ζ(s)\| ≤ \|ζ(s+ε)\|\|Γr_ε(s)\|) | l. 172–182 | [C2] Hadamard product (Davenport ch. 12) |
| LemmaA | ζ ≠ 0 on Re s > θ (½ ≤ θ < 1), δ ∈ (0,(1−θ)/2], η' > 0 ⇒ ∃C, σ ≥ θ+2δ: \|ζ\|⁻¹ ≤ C(2+\|t\|)^{η'}, and \|ζ\| ≤ C(2+\|t\|)^{η'} if \|t\| ≥ 1 | l. 183–207 | [C4] Borel–Carathéodory, [C5] three circles, [C6] (Titchmarsh 4.11) |
| LemmaP | the same zero-free hypothesis, η ∈ (0,1), δ, η' ∈ (0,½], σ − η ≥ θ+2δ ⇒ \|P_N(s) − 1/ζ(s)\| ≤ C N^{−η}(2+\|t\|)^{η'}; σ − η > 1 ⇒ ≤ (1/(1−η) + 1/η)N^{−η} | l. 208–225 | Riesz–Perron formula, proved in PROOF.md |

Lemma M (l. 131–138) is proved, not assumed. Prop N uses none of the hypotheses.

## Compile commands and output

Each module compiled singly, oleans to .scratch/nblplean/olean (not .lake/build), with Lean/audit/mk.sh:

~~~
O=$PWD/.scratch/nblplean/olean
perl -e 'alarm 1500; exec @ARGV' lake env sh -c "LEAN_PATH=\"$O:\$LEAN_PATH\" exec lean -o $O/Riemann/Deep/NBLp/Lean/$n.olean Riemann/Deep/NBLp/Lean/$n.lean"
~~~

| Module | Output | EXIT | Peak (8 s sampler) |
|---|---|---|---|
| Basic | (none) | 0 | 343 MB |
| PropN | 'NBLp.propN' depends on axioms: [propext, Classical.choice, Quot.sound] | 0 | 456 MB |
| Toolkit | 'NBLp.lemmaM' depends on axioms: [propext, Classical.choice, Quot.sound] | 0 | 487 MB |

No errors, no sorry warnings. These ran under the brief's original memory rule ("330–450 MB accepted for lean"), before the orchestrator's correction to 290 MB; they are not re-runnable under the new cap.

## What blocked §7

Peak memory footprint of lean on a file with only the imports and example : True := trivial (direct lean, /usr/bin/time -l, or the memcap watchdog):

| Import | Peak |
|---|---|
| Mathlib.Data.Real.Basic | 223 MB |
| Mathlib.Analysis.SpecialFunctions.Pow.Real | 334 MB |
| Mathlib.MeasureTheory.Function.LpSeminorm.Defs | 343 MB |
| LpSeminorm.TriangleInequality + Function.Floor + Measure.Lebesgue.Basic | 318 MB (killed) |
| Mathlib.Analysis.SpecialFunctions.Gamma.Basic | 324 MB (killed) |
| Mathlib.NumberTheory.LSeries.RiemannZeta | 342 MB (killed) |
| Mathlib.Analysis.SpecialFunctions.Gamma.Deligne | 333 MB (killed) |

On top of the imports, simp, field_simp, norm_num, linarith and nlinarith build the full simp set (another 130–260 MB); ring, positivity, gcongr, simp only stay flat. ChiInW needs eLpNorm, x^s needs Pow.Real, ζ needs LSeries.RiemannZeta: every §7 file exceeds 290 MB on import. A cap-compliant rewrite (zeta as a parameter Z, tactic-light proofs: .scratch/nblplean/light/) was abandoned for the same reason. Two 10 s probes (Pow.Real, LpSeminorm.Defs) peaked at 334 and 343 MB; the watchdog sampled 285 and 278 MB and did not kill them. The watchdog through lake env samples the lake process (about 120 MB) and can miss lean's peak, so its "peak 120 MB" for the light definitions file is not evidence of compliance. Every probe or attempt the watchdog killed (imports above; three attempts at the light toolkit file) went above 290 MB for a few seconds before the kill (313–408 MB).

## §7 drafts (NOT checked by Lean)

- draft/Sec7.lean: §7 objects φ_ε, E_ε, Φ_ε, H_ε, H_{ε,N}, Φ_{ε,N}, R_ε = ζ(s)·invZeta(s+ε), Ψ_ε and r, ε₀, C_Γ, K_p, A_p (l. 247–257, §1); [C8] facts from Mathlib's ζ (1/ζ(1+ε) ≤ ε via Σn^{−σ} ≥ ∫₁^∞ x^{−σ}; ζ decreasing on (1,∞)); A1's identities on (0,1]; A5 (Lemma M + B(ε) ≤ K_pε ⇒ ‖χ + H_ε‖_p ≤ A_pε^{1−1/p} + ε); A6 (‖Φ_ε − Φ_{ε,N}‖²_{L²_σ''} ≤ CN^{−ε} and Lemma P at 1+ε ⇒ ‖H_{ε,N} − H_ε‖_p ≤ C'N^{−ε/2}); A7 (bound_of_A5_A6: N₀(p,ε); chiInW_of_eps: N → ∞ for fixed ε, then ε → 0); chiInW_of_A4_A6. First compile: 6 errors, all fixed; the second compile was killed by the supervisor at 361 MB.
- draft/A4.lean: A4 (B(ε) ≤ K_pε) from A3 on σ' = 1/p − ε/2, Lemma H, Lemma Γ(b), the toolkit's Plancherel identity and A4(b) (JEval: J converges and equals 2π[R(2σ')/(2σ') − 1/((2σ'−1)ζ(1+ε))], l. 280–291), via \|1−R\|² = 1 − 2Re R + \|R\|², ∫(σ²+t²)⁻¹ = π/σ and e^{8ε}y^ε ≤ 1 + εe⁴(8 + log y)y^{1/2}. Never compiled.
- To check the drafts (about 350–500 MB per compile, one at a time, after Basic, PropN, Toolkit): lake env sh -c "LEAN_PATH=\"$O:\$LEAN_PATH\" lean -j1 -o $O/Riemann/Deep/NBLp/Lean/draft/Sec7.olean Riemann/Deep/NBLp/Lean/draft/Sec7.lean", then the same for draft/A4.
- Not drafted: A2, A3 (Lemma PW applied to Ψ_ε; needs the removable zero of 1/ζ at 1), A4(b) itself (residues, or Parseval plus the Dirichlet series of φ_ε at 2σ'), A6's L²_σ'' bound (Lemmas T, H, A, Γ(b), P), the corollaries.

## Fidelity

- ρ_k, f_k, χ, V, W_p, d_{N,p}: PROOF.md l. 10–17. ChiInW p is χ ∈ W_p = closure of span_ℂ{f_k : k ≥ 2}: every element of the span is Σ_{k=2}^N c_k f_k (zero padding), so closure membership is "∀ε ∃N, c: ‖χ − Σ c_k f_k‖_p < ε", i.e. d_{N,p} → 0 (l. 17, 36). Coefficients complex as in l. 17; (0,1] instead of (0,1) changes a null set.
- The page's (0,1) form: fk_eq_rhoNB proves f_k = ρ_{1/k} with ρ_a(x) = {a/x} − a{1/x} (l. 11), definitionally Suff.lean's rhoNB.
- propN: Proposition N (l. 105–114) for every real 1 < p (l. 105: 1 < p < ∞), any zero w₀ with Re w₀ > 1/p; Re w₀ < 1 is derived (Mathlib's riemannZeta_ne_zero_of_one_le_re), not assumed. zeroFree_of_chiInW is Theorem A (ii) ⇒ (i) (l. 33–35); (i) is ZeroFree (1/p) (l. 34), harmless at s = 1 (riemannZeta 1 ≠ 0).
- Báez-Duarte [C10] (l. 82, arXiv math/0202141): RH ⟺ χ ∈ W_2. ChiInW 2 is that statement's W_2 membership; propN at p = 2 is its necessity half (a zero with Re > ½ gives ¬ChiInW 2), the same direction as NymanBeurling/Suff.lean's RH_of_NBApprox for the continuous family.
- Toolkit: as in the table. Reviewer caveats: the docstring's "^" is PROOF.md's angular Fourier transform ĥ(t) = ∫h(u)e^{−iut}du ([C7], l. 75), not Mathlib's 2π-normalized 𝓕 (then 𝓜_σF(t) = 𝓕g(t/2π)); the toolkit's 1/(2π) and its L¹ formula are consistent with each other. LemmaA at s = 1 bounds a finite junk value of ‖ζ(1)‖⁻¹, absorbed by the existential constant. Lemma Γ(a) is not stated (only Γ(b) is used in §7).

Reviewer (subagent, statement fidelity): Parts 1–2 (the checked files) FAITHFUL, with the caveats above; part 3 (drafts) not certifiable, definitions and §7 statements not supplied (and the drafts are unchecked anyway). Verdict: Lean/audit/reviewer_verdict.md.

## .lake

.lake was a real directory (other agents' shared Mathlib cache, 7.4 GB); left exactly as found. No lake exe cache get, no placeholder restore; oleans went to .scratch/nblplean/olean.

## FORMALIZATIONS.md entry

- NBLp (our discovery [nyman-beurling-lp](discoveries/nyman-beurling-lp.md)), PARTIAL: Proposition N (a zero with Re ρ > 1/p ⇒ 1_(0,1) ∉ L^p-closure of span{ρ_{1/k} : k ≥ 2}, i.e. Theorem A (ii) ⇒ (i)) proved from Mathlib alone (NBLp.propN, standard axioms) and Lemma M proved, §4–§5 lemmas as named hypotheses; §7 (i) ⇒ (ii) drafted, unchecked (imports exceed the 290 MB cap). Reviewer: checked statements FAITHFUL. [LEAN.md](Riemann/Deep/NBLp/LEAN.md)



