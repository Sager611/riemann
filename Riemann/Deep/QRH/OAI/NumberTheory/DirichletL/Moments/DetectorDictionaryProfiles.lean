-- Ported from https://github.com/openai/math (commit 1581d4cbdd69a0faf34b2089614fe2a814316a6a), file lean/OAI/NumberTheory/DirichletL/Moments/DetectorDictionaryProfiles.lean.
-- Apache License 2.0 (Riemann/Deep/QRH/LICENSE, NOTICE). Modified: import paths rewritten to Riemann.Deep.QRH.*.
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Moments.DetectorDictionary
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Hecke.InverseAmplificationRowwise
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Hecke.DetectorDyadicCutoff

namespace OAI

noncomputable section
open scoped Classical ContDiff ComplexConjugate
namespace SevenEighths.CenteredMomentDetectorDictionary
open HeckeFamily HeckeInverseAmplification HeckeDetectorCoefficientTransfer
open HeckeDetectorRowwisePolynomial HeckeDetectorDyadicProfiles

lemma logTest_eq_iterate (W : ℝ→ℂ) (n : ℕ) : logTest W n=(logProfile^[n]) W := by
  induction n with
  | zero => rfl
  | succ n ih => simp only [logTest,Function.iterate_succ_apply',ih]

lemma oriented_logTest (reverse : Bool) (W : ℝ→ℂ) (n : ℕ) :
    orientedProfile reverse ((logProfile^[n]) W)=logTest (orientedProfile reverse W) n := by
  rw [←logTest_eq_iterate]
  cases reverse
  · rfl
  · induction n with
    | zero => rfl
    | succ n ih =>
      simp only [logTest,orientedProfile,ite_true] at ih ⊢
      funext x
      rw [←logProfile_conjugate]
      exact congrArg (fun f : ℝ→ℂ=>logProfile f x) ih

theorem detector_profile_support (reverse : Bool) (n : ℕ) (σ t : ℝ) :
    Function.support (twistProfile (orientedProfile reverse ((logProfile^[n]) positiveAnnular)) σ t)
      ⊆Set.Icc (1/4) (9/4) := by
  apply (oriented_twist_support reverse _ σ t).trans
  rw [←logTest_eq_iterate]
  exact (logTest_support positiveAnnular n).trans positiveAnnular_support

theorem detector_profile_smooth (reverse : Bool) (n : ℕ) (σ t : ℝ) :
    ContDiff ℝ ∞ (twistProfile (orientedProfile reverse ((logProfile^[n]) positiveAnnular)) σ t) := by
  apply twistProfile_smooth _ σ t (1/4) (9/4) (by norm_num)
  · cases reverse <;> simp only [orientedProfile,ite_true]
    · rw [←logTest_eq_iterate]
      exact (logTest_support positiveAnnular n).trans positiveAnnular_support
    · rw [conjugate_profile_support,←logTest_eq_iterate]
      exact (logTest_support positiveAnnular n).trans positiveAnnular_support
  · have hw := logTest_smooth positiveAnnular n (1/4) (9/4) (by norm_num)
      positiveAnnular_support positiveAnnular_smooth
    rw [logTest_eq_iterate] at hw
    cases reverse
    · exact hw
    · exact conjugate_profile_smooth _ hw

end SevenEighths.CenteredMomentDetectorDictionary

end

end OAI
