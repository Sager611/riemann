module

import all Riemann.Deep.Combos.Lean.Count
import all Riemann.Deep.Combos.Lean.Expo1pA
import all Riemann.Deep.Combos.Lean.Expo1pB
import all Riemann.Deep.Combos.Lean.Expo1pC1
import all Riemann.Deep.Combos.Lean.Expo1pC2
import all Riemann.Deep.Combos.Lean.Expo1pC3
import all Riemann.Deep.Combos.Lean.Expo1pC
import all Riemann.Deep.Combos.Lean.Expo1pD
import all Riemann.Deep.Combos.Lean.Thm1pL
import all Riemann.Deep.Combos.Lean.Thm1pK
import all Riemann.Deep.Combos.Lean.Thm1pI
import all Riemann.Deep.Combos.Lean.Thm1pW
import all Riemann.Deep.Combos.Lean.Thm1pLg
import all Riemann.Deep.Combos.Lean.Thm1pSm
import all Riemann.Deep.Combos.Lean.Thm1pD
import all Riemann.Deep.Combos.Lean.Thm1pMain
public meta import Lean.Elab.Command

open Lean Elab Command

/-! axiom audit of the Theorem 1′(a) chain without descending into Mathlib proofs: with the private
parts of our modules loaded (import all), walk every RingClass constant reachable from thm1'a (types
and values); report RingClass constants without a value, and every external constant that is
sorryAx / Lean.ofReduceBool / Lean.trustCompiler or comes from outside Init/Std/Lean/Mathlib/Batteries/Qq/Aesop. -/

run_cmd do
  let env ← getEnv
  let mut todo : Array Name := #[`RingClass.thm1'a]
  let mut seen : NameSet := {}
  let mut ext : NameSet := {}
  let mut noval : Array Name := #[]
  let mut nrc : Nat := 0
  for _ in [0:200000] do
    match todo.back? with
    | none => break
    | some c =>
      todo := todo.pop
      if seen.contains c then continue
      seen := seen.insert c
      let own := c.getRoot == `RingClass || (match env.getModuleIdxFor? c with
        | some i => ((env.header.moduleNames[i.toNat]?).getD .anonymous).getRoot == `Riemann
        | none => false)
      if own then
        nrc := nrc + 1
        match env.find? c with
        | some ci =>
          for d in ci.type.getUsedConstants do todo := todo.push d
          match ci with
          | .thmInfo v => for d in v.value.getUsedConstants do todo := todo.push d
          | .defnInfo v => for d in v.value.getUsedConstants do todo := todo.push d
          | .opaqueInfo v => for d in v.value.getUsedConstants do todo := todo.push d
          | .inductInfo _ | .ctorInfo _ | .recInfo _ | .quotInfo _ => pure ()
          | .axiomInfo _ => noval := noval.push c
        | none => noval := noval.push c
      else
        ext := ext.insert c
  let bad := ext.toList.filter fun n =>
    n == `sorryAx || n == `Lean.ofReduceBool || n == `Lean.trustCompiler ||
    (match env.getModuleIdxFor? n with
     | some i =>
       let m := (env.header.moduleNames[i.toNat]?).getD .anonymous
       !([`Init, `Std, `Lean, `Mathlib, `Batteries, `Qq, `Aesop].contains m.getRoot)
     | none => true)
  logInfo m!"own constants walked (RingClass and private auxiliaries): {nrc}; external constants: {ext.toList.length}; RingClass constants without value: {noval}; flagged external constants: {bad}"
