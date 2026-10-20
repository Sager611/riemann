-- Ported from https://github.com/openai/math (commit 1581d4cbdd69a0faf34b2089614fe2a814316a6a), file lean/OAI/NumberTheory/DirichletL/Detector/GramNumerator.lean.
-- Apache License 2.0 (Riemann/Deep/QRH/LICENSE, NOTICE). Modified: import paths rewritten to Riemann.Deep.QRH.*.
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Detector.LowGramFixedCoefficient
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Detector.GramMean

namespace OAI

noncomputable section
open scoped Classical
namespace SevenEighths.ProbeGramCommon
open CanonicalQuadraticSieve CanonicalRowCompletion CompletedGauss ConcreteTraceCRT
local notation "O" => ActualEisensteinCubic.O
local notation "λ₀" => ConcretePrimeRowBridge.goodLambda

def numeratorExtension (u : Oˣ) (a b : ℕ) (r : O) (hr : Supported (Ideal.span {r})) : O→*ℂ :=
  numeratorBadTwist u a b r hr*movingNumeratorRow r hr

lemma numeratorExtension_norm (u : Oˣ) (a b : ℕ) (r : O) (hr : Supported (Ideal.span {r})) (n : O) :
    ‖numeratorExtension u a b r hr n‖≤1 := by
  rw [numeratorExtension,MonoidHom.mul_apply,norm_mul]
  exact (mul_le_of_le_one_left (norm_nonneg _) (numeratorBadTwist_norm _ _ _ _ _ _)).trans (movingNumeratorRow_norm _ _ _)

lemma numeratorExtension_periodic (u : Oˣ) (a b : ℕ) (r : O) (hr : Supported (Ideal.span {r})) :
    CanonicalCoefficientClass.FactorsModulo (Ideal.span {(72:O)}*Ideal.span {r})
      (numeratorExtension u a b r hr) := by
  intro x y hxy
  change numeratorBadTwist u a b r hr x*movingNumeratorRow r hr x=
    numeratorBadTwist u a b r hr y*movingNumeratorRow r hr y
  rw [numeratorBadTwist_periodic _ _ _ _ _ x y (Ideal.mul_le_left hxy),
    movingNumeratorRow_periodic _ _ x y (Ideal.mul_le_right hxy)]

lemma numeratorExtension_primary (u : Oˣ) (a b : ℕ) (r k n : O)
    (hr : Supported (Ideal.span {r})) (hpr : λ₀^2∣r-1)
    (hk : k=u.val*λ₀^a*(2:O)^b*r)
    (hn : Supported (Ideal.span {n})) (hpn : λ₀^2∣n-1) :
    numeratorExtension u a b r hr n=idealRowHom k (Ideal.span {n}) := by
  rw [hk]
  exact (idealRowHom_factor_numerator u a b r n hr hpr hn hpn).symm

lemma numeratorExtension_neg_primary (u : Oˣ) (a b : ℕ) (r k n : O)
    (hr : Supported (Ideal.span {r})) (hpr : λ₀^2∣r-1)
    (hk : k=u.val*λ₀^a*(2:O)^b*r)
    (hn : Supported (Ideal.span {n})) (hpn : λ₀^2∣n-1) :
    numeratorExtension (-u) a b r hr n=idealRowHom (-k) (Ideal.span {n}) := by
  apply numeratorExtension_primary _ _ _ _ _ _ hr hpr _ hn hpn
  simp only [Units.val_neg,hk]
  ring

lemma numerator_good_norm_le (u : Oˣ) (a b : ℕ) (r k : O)
    (hk : k=u.val*λ₀^a*(2:O)^b*r) (hk0 : k≠0) :
    Ideal.absNorm (Ideal.span {r})≤Ideal.absNorm (Ideal.span {k}) := by
  apply Nat.le_of_dvd
  · exact Nat.pos_of_ne_zero (Ideal.absNorm_eq_zero_iff.not.mpr (Ideal.span_singleton_eq_bot.not.mpr hk0))
  · apply map_dvd Ideal.absNorm
    rw [Ideal.span_singleton_dvd_span_singleton_iff_dvd]
    exact ⟨u.val*λ₀^a*(2:O)^b,by rw [hk];ring⟩

end SevenEighths.ProbeGramCommon
end

end OAI
