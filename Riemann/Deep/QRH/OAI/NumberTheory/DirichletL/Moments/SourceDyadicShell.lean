-- Ported from https://github.com/openai/math (commit 1581d4cbdd69a0faf34b2089614fe2a814316a6a), file lean/OAI/NumberTheory/DirichletL/Moments/SourceDyadicShell.lean.
-- Apache License 2.0 (Riemann/Deep/QRH/LICENSE, NOTICE). Modified: import paths rewritten to Riemann.Deep.QRH.*.
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Moments.SourceLowerSupport
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Moments.SecondActiveDyadic

namespace OAI

noncomputable section
open scoped Classical BigOperators

namespace SevenEighths.CenteredMomentSourceDyadicShell
open CenteredMomentSectorLocalization CenteredMomentSourceLowerSupport
local notation "O" => ActualEisensteinCubic.O

lemma residual_window_bounds (C I:Ideal O) (hC:C≠0) (L U:ℝ) (n:ℤ)
    (hL:L≤(Ideal.absNorm (C*I):ℝ)) (hU:(Ideal.absNorm (C*I):ℝ)≤U)
    (hw:dyadicWeight n (Ideal.absNorm I:ℝ)≠0) :
    L/(Ideal.absNorm C:ℝ)≤dyadicScale n ∧
      dyadicScale n≤4*U/(Ideal.absNorm C:ℝ) := by
  have hc:(0:ℝ)<Ideal.absNorm C:=by
    exact_mod_cast Nat.pos_of_ne_zero (Ideal.absNorm_eq_zero_iff.not.mpr hC)
  have hs:=dyadicWeight_support n hw
  simp only [map_mul,Nat.cast_mul] at hL hU
  constructor
  · apply (div_le_iff₀ hc).mpr
    exact hL.trans (by nlinarith [hs.2])
  · apply (le_div_iff₀ hc).mpr
    nlinarith [hs.1]

lemma raw_window_lower (C I:Ideal O) (hC:C≠0) (a H:ℝ) (ha:0<a) (n:ℤ)
    (hL:a*H≤(Ideal.absNorm (C*I):ℝ))
    (hw:dyadicWeight n (Ideal.absNorm I:ℝ)≠0) :
    H/(Ideal.absNorm C:ℝ)≤dyadicScale n/a := by
  have hc:(0:ℝ)<Ideal.absNorm C:=by
    exact_mod_cast Nat.pos_of_ne_zero (Ideal.absNorm_eq_zero_iff.not.mpr hC)
  have hs:=(dyadicWeight_support n hw).2.le
  simp only [map_mul,Nat.cast_mul] at hL
  apply (div_le_div_iff₀ hc ha).mpr
  nlinarith

end SevenEighths.CenteredMomentSourceDyadicShell

end

end OAI
