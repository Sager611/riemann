PORTED: OpenAI's proof that zeta has no zero with Re s > 7/8 builds sorry-free on our Lean/Mathlib; axioms exactly propext, Classical.choice, Quot.sound.

What builds (Lean v4.35.0-rc2, mathlib 0653561):
- lake build Riemann.Deep.QRH.OAI.NumberTheory.DirichletL.Nonvanishing: green (7059 jobs).
- lake build Riemann.Deep.QRH.Bridge: green. #print axioms gives [propext, Classical.choice, Quot.sound] for
  OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re, OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re,
  Fayed.QRH.quasiRH_seven_eighths and Fayed.QRH.zeta_ne_zero_of_seven_eighths_lt_re (check file: .scratch/qrh2/ax.lean).
- No sorry, admit, native_decide, axiom, unsafe, implemented_by, extern or custom syntax anywhere in the port.

What was ported:
- Import closure of OAI.NumberTheory.DirichletL.Nonvanishing at openai/math 1581d4cb: 2924 files, 486,490 lines
  (longest import chain 404 files), under OAI/. Every file's git blob hash matches upstream before editing.
- Third-party: PrimeNumberTheoremAnd.Wiener closure (4 files) and RellichKondrachov (15 files) under Ext/, with
  OpenAI's lean4341 compatibility patches applied. Licences: all Apache 2.0 (LICENSE, NOTICE).
- Upstream pins Lean v4.34.1 + mathlib d13f23b7; ours is 69 mathlib commits away.

Bridge (Bridge.lean):
- R007's QuasiRH δ means "no zero with Re s > 1 - δ" (QuasiRH (1/2) is RH). The 7/8 half-plane is therefore
  QuasiRH (1/8); the brief's "QuasiRH (7/8)" would claim no zero right of 1/8, which is false. The theorem
  Fayed.QRH.quasiRH_seven_eighths proves QuasiRH (1/8). Corollary zeta_ne_zero_of_seven_eighths_lt_re: no zero of
  zeta with Re s > 7/8, so any counterexample to Fayed's RH claim lies in 1/2 < Re s ≤ 7/8.

Local changes beyond import paths (all keep the proofs; NOTICE points here; headers not edited, to avoid a rebuild):
1. 33 multi-line set literals in 19 files: {X} -> {(X)}. Lean 4.35 parses a set literal whose continuation line
   has application arguments left of the brace as set-builder notation. List: .scratch/qrh2/fix_braces_log.json.
2. congr 1; ext v; simp only [..] -> congr 1; all_goals (ext v; simp only [..]) in Moments/FirstSectors,
   Moments/CommonAllocationSum and one more (.scratch/qrh2/fix_congr_log.json): congr 1 now closes the goal.
3. Reflection/SectorData: ULift.ext x y h -> ULift.ext h (x, y are now implicit).
4. PNT's Fourier.lean/Sobolev.lean wrapped in namespace QRHPNT, Wiener.lean opens it: the same names exist in
   Riemann/Deep/WeilEF/Z23/FromPNTPlus, so without this the two ports could not be imported together.

Costs and caveats:
- Build: about 3.5 h of wall time at 4 Lean threads; peak was ~4 parallel Lean processes.
- Disk: oleans 1.4 GB, plus Lake's C files in .lake/build/ir/Riemann/Deep/QRH (~8 GB, not needed for checking).
- Riemann.lean does not import the port (not in the default build); add import Riemann.Deep.QRH.Bridge to include it.
- Upstream source fetched by curl from raw.githubusercontent.com; scripts and closure lists in .scratch/qrh2/.
