import Riemann.Deep.ExplorePointProc.Lean.Theorem1
import Riemann.Deep.ExplorePointProc.Lean.PL
import Riemann.Deep.ExplorePointProc.Lean.Herm

/-!
# The certificate cert768 and the corollary: D(λ) > 0 for every λ ≥ 0.768 (DERIVATION §5)

cert768 (certs/cert768.json, checked by code/cert_aug.py with EPS = 1/1000, m = 5): λ = 96/125,
θ = 399993/700000 = 4/7 − 10⁻⁵, U = 3λ = 288/125, s = 11/10, β = 416188506333/640725853462,
f(α) = Σ_{n≤21} B_n H_{2n}(√(2π)α/s) e^{−πα²/s²}, f^(u) = s Σ (−1)^n B_n H_{2n}(√(2π)su) e^{−πs²u²},
A piecewise linear on the 78 nodes λ + k·dg, dg = 192/9625, A(x) = A(U) = 0 for x ≥ U.
Each checked fact is one named hypothesis below (exact statement; checker lines in the docstrings).
Rational constants (from .scratch/pplean/cert_exact.out): c_cert = −0.3478616532233634 (binary64) ≤
c0 = −347861653223/10¹²; J_θ(f) ≥ J_lo = 0.306839281979 − 3.75·10⁻¹³ ≥ 153419640989/(5·10¹¹);
f^(0) = 1.0011 − 3·10⁻¹⁶ ≤ 10011/10000; min_k A_k = A_0 = −378347742845/758342134484.
Proved here: (H) on [λ, ∞)² from R1, R2 and R3 (R3 derived), the margin, Theorem 1 at λ = 96/125, monotonicity.
-/

namespace PPLean

/-- cert_aug.py l.43–46: p(0) > 0 and, by an exact Sturm sequence, no root of p in (0, ∞), where
f(α) = p(2πα²/s²) e^{−πα²/s²}: f > 0 on [0, ∞). -/
def CertPos (f : ℝ → ℝ) : Prop := ∀ α : ℝ, 0 ≤ α → 0 < f α

/-- l.47–52: q(w0) < 0 and no root of q in (w0, ∞), w0 = 2·3.14159265·s²U² < 2πs²U², where
f^(u) = s q(2πs²u²) e^{−πs²u²}: f^ < 0 on [U, ∞), U = 288/125. -/
def CertNeg (C : TestCalc) (f : ℝ → ℝ) : Prop := ∀ u : ℝ, 288 / 125 ≤ u → C.hat f u < 0

/-- l.53–66 (Arb, 1200 bits): J_θ(f) = 0.306839281979 ± 3.75·10⁻¹³, θ = 399993/700000. -/
def CertJ (C : TestCalc) (f : ℝ → ℝ) : Prop :=
  (153419640989 / 500000000000 : ℝ) ≤ C.J (399993 / 700000) f

/-- l.67: f^(0) = s·q(0), an exact rational 1.0010999999999997… ≤ 1.0011. -/
def CertF0 (C : TestCalc) (f : ℝ → ℝ) : Prop := C.hat f 0 ≤ 10011 / 10000

/-- The certificate's A (l.98–103, 129): piecewise linear through the 78 nodes (values in certs/cert768.json,
all in [A_0, 0]) and A(x) = A(U) = 0 for x ≥ U; so A(U) = 0 and A_0 ≤ A ≤ 0 on [λ, ∞). -/
def CertA (A : ℝ → ℝ) : Prop :=
  (∀ y : ℝ, 288 / 125 ≤ y → A y = 0) ∧
    ∀ y : ℝ, 96 / 125 ≤ y → (-378347742845 / 758342134484 : ℝ) ≤ A y ∧ A y ≤ 0

/-- Φ(x, y) = f^(x) + f^(x + y) + A(x) − A(y) + β(1 − x), the left side of (H). -/
noncomputable def Phi (C : TestCalc) (f A : ℝ → ℝ) (x y : ℝ) : ℝ :=
  C.hat f x + C.hat f (x + y) + A x - A y + 416188506333 / 640725853462 * (1 - x)

/-- R1 = [λ, U]² (l.98–105, 126–127): Arb upper values at the lattice corners (δ = dg/5), max −0.3478888283,
plus the bilinear remainder 3M₂δ²/8 = 2.718·10⁻⁵, M₂ ≥ sup_{[λ,2U]} |f^''| (l.77–96); c_cert = R1 + 10⁻¹². -/
def CertR1 (C : TestCalc) (f A : ℝ → ℝ) : Prop :=
  ∀ x y : ℝ, 96 / 125 ≤ x → x ≤ 288 / 125 → 96 / 125 ≤ y → y ≤ 288 / 125 →
    Phi C f A x y ≤ -347861653223 / 1000000000000

/-- R2 (l.128): y ≥ U drops f^(x + y) < 0 and A(y) = A(U): the one-dimensional bound
max_{x∈[λ,U]} f^(x) + A(x) − A(U) + β(1 − x) + M₂δ²/8 = −0.5640124229. -/
def CertR2 (C : TestCalc) (f A : ℝ → ℝ) : Prop :=
  ∀ x : ℝ, 96 / 125 ≤ x → x ≤ 288 / 125 →
    C.hat f x + A x - A (288 / 125) + 416188506333 / 640725853462 * (1 - x)
      ≤ -347861653223 / 1000000000000

/-- (H) on [λ, ∞)² with c = c0: R1 (assumed), R2 (assumed, 1D), R3 = x ≥ U (proved: A(x) = 0, f^(x), f^(x + y) < 0,
β(1 − x) ≤ β(1 − U), −A(y) ≤ −A_0; R3 = −0.3481091233 ≤ c0, l.129). -/
lemma certHH {C : TestCalc} {f A : ℝ → ℝ} (hneg : CertNeg C f) (hA : CertA A)
    (hR1 : CertR1 C f A) (hR2 : CertR2 C f A) :
    HH (96 / 125) (416188506333 / 640725853462) (-347861653223 / 1000000000000) (C.hat f) A := by
  intro x y hx hy
  obtain ⟨htail, hrange⟩ := hA
  rcases le_or_gt (288 / 125) x with hxU | hxU
  · have h1 := hneg x hxU
    have h2 := hneg (x + y) (by linarith)
    have h3 := htail x hxU
    have h4 := (hrange y hy).1
    have h5 : (416188506333 / 640725853462 : ℝ) * (1 - x)
        ≤ 416188506333 / 640725853462 * (1 - 288 / 125) :=
      mul_le_mul_of_nonneg_left (by linarith) (by norm_num)
    have h6 : -(-378347742845 / 758342134484 : ℝ) + 416188506333 / 640725853462 * (1 - 288 / 125)
        ≤ -347861653223 / 1000000000000 := by norm_num
    rw [h3]
    linarith
  · rcases le_or_gt (288 / 125) y with hyU | hyU
    · have h1 := hneg (x + y) (by linarith)
      have h2 := htail y hyU
      have h3 := htail (288 / 125) le_rfl
      have h4 := hR2 x hx hxU.le
      rw [h3] at h4
      rw [h2]
      linarith
    · exact hR1 x y hx hxU.le hy hyU.le

/-- The margin J − f^(0) − 2c0 ≥ 153419640989/(5·10¹¹) − 1.0011 + 2·0.347861653223 = 0.001462588424 > 0. -/
lemma cert_margin {C : TestCalc} {f : ℝ → ℝ} (hJ : CertJ C f) (hF0 : CertF0 C f) :
    C.hat f 0 + 2 * (-347861653223 / 1000000000000) < C.J (399993 / 700000) f := by
  unfold CertJ at hJ
  unfold CertF0 at hF0
  have h : (10011 / 10000 : ℝ) + 2 * (-347861653223 / 1000000000000) < 153419640989 / 500000000000 := by
    norm_num
  linarith

/-- D(96/125) > 0 from Theorem 1 and the certificate. -/
theorem dpos_768 {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z)
    (h12 : Steps12 Z C) (hFC : FourierCc2 C) (hT : LemmaT' C) {f A : ℝ → ℝ}
    (hS : C.Sch f) (he : ∀ x, f (-x) = f x) (hpos : CertPos f) (hneg : CertNeg C f)
    (hJ : CertJ C f) (hF0 : CertF0 C f) (hA : CertA A) (hR1 : CertR1 C f A) (hR2 : CertR2 C f A) :
    DPos Z (96 / 125) := by
  have hfn : ∀ x, 0 ≤ f x := by
    intro x
    rcases le_or_gt 0 x with hx | hx
    · exact (hpos x hx).le
    · rw [← he x]
      exact (hpos (-x) (by linarith)).le
  have hAb : ∃ K : ℝ, ∀ y, 96 / 125 ≤ y → |A y| ≤ K := by
    refine ⟨378347742845 / 758342134484, fun y hy => abs_le.2 ⟨?_, ?_⟩⟩
    · have := (hA.2 y hy).1
      linarith
    · have := (hA.2 y hy).2
      linarith
  exact theorem1 hZ hsc hRvM h12 hFC hT (θ := 399993 / 700000) (U := 288 / 125)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    hS he hfn (fun u hu => (hneg u hu).le) hAb (by norm_num) (certHH hneg hA hR1 hR2)
    (cert_margin hJ hF0)

/-- **Corollary** (DERIVATION §5): D(λ) > 0 for every λ ≥ 0.768 (monotonicity of D in λ). -/
theorem dpos_ge_768 {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered) (hsc : ScaleToInf Z) (hRvM : RvM Z)
    (h12 : Steps12 Z C) (hFC : FourierCc2 C) (hT : LemmaT' C) {f A : ℝ → ℝ}
    (hS : C.Sch f) (he : ∀ x, f (-x) = f x) (hpos : CertPos f) (hneg : CertNeg C f)
    (hJ : CertJ C f) (hF0 : CertF0 C f) (hA : CertA A) (hR1 : CertR1 C f A) (hR2 : CertR2 C f A) :
    ∀ lam : ℝ, 96 / 125 ≤ lam → DPos Z lam := fun _ hlam =>
  (dpos_768 hZ hsc hRvM h12 hFC hT hS he hpos hneg hJ hF0 hA hR1 hR2).mono hlam

/-- CertA holds for the certificate's own A (PL.lean: last node 0, all nodes in [A_0, 0]). -/
lemma certA_Acert : CertA Acert := ⟨Acert_tail, Acert_range⟩

/-- **Corollary, exact form.** Named: the Sturm facts CertP, CertQ, the Arb facts CertJ, CertR1, CertR2 (for the
certificate's own piecewise-linear A), the Hermite–Gaussian identities HermRep, f Schwartz and even, and the
analytic inputs of Theorem 1. Proved: f > 0, f^ < 0 on [U, ∞), f^(0) ≤ 1.0011, A's tail and range, R3, the
margin, Theorem 1 at λ = 96/125 and the monotonicity of D. -/
theorem dpos_ge_768_exact {Z : Zeros} {C : TestCalc} (hZ : Z.Ordered) (hsc : ScaleToInf Z)
    (hRvM : RvM Z) (h12 : Steps12 Z C) (hFC : FourierCc2 C) (hT : LemmaT' C) {f : ℝ → ℝ}
    (hS : C.Sch f) (he : ∀ x, f (-x) = f x) (hrep : HermRep f (C.hat f)) (hP : CertP) (hQ : CertQ)
    (hJ : CertJ C f) (hR1 : CertR1 C f Acert) (hR2 : CertR2 C f Acert) :
    ∀ lam : ℝ, 96 / 125 ≤ lam → DPos Z lam :=
  dpos_ge_768 hZ hsc hRvM h12 hFC hT hS he (fun α _ => pos_of_rep hrep hP α)
    (fun u hu => neg_of_rep hrep hQ u hu) hJ (fh0_of_rep hrep) certA_Acert hR1 hR2

end PPLean
