-- Ported from https://github.com/openai/math (commit 1581d4cbdd69a0faf34b2089614fe2a814316a6a), file lean/OAI/NumberTheory/DirichletL/Dictionary/InverseUniformScale.lean.
-- Apache License 2.0 (Riemann/Deep/QRH/LICENSE, NOTICE). Modified: import paths rewritten to Riemann.Deep.QRH.*.
import Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Moments.DetectorDictionaryUniformTests

namespace OAI

noncomputable section
open scoped Classical SchwartzMap
namespace SevenEighths.DetectorDictionaryInverseUniform
open HeckeInverseAmplification

def scaleCLM : 𝓢(ℝ,ℂ)→L[ℝ]𝓢(ℝ,ℂ) :=
  (-1/2:ℝ) • ContinuousLinearMap.id ℝ (𝓢(ℝ,ℂ))-
    (SchwartzMap.smulLeftCLM ℂ (fun x : ℝ=>x)).comp (SchwartzMap.derivCLM ℝ ℂ)

lemma scaleCLM_apply (W : 𝓢(ℝ,ℂ)) (x : ℝ) : scaleCLM W x=scaleProfile W x := by
  simp only [scaleCLM,sub_apply,smul_apply,
    ContinuousLinearMap.id_apply,ContinuousLinearMap.comp_apply,
    SchwartzMap.smulLeftCLM_apply_apply Function.HasTemperateGrowth.id',
    SchwartzMap.derivCLM_apply,Complex.real_smul,scaleProfile]
  push_cast
  ring

end SevenEighths.DetectorDictionaryInverseUniform

end

end OAI
