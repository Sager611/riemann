STATIC AUDIT PASSED: source-level screens of all 175,187 files and hand review of every metaprogram found no kernel bypass, no unfaithful statement and no source-level trigger pattern of a known Lean v4.33.1 kernel bug (the 4.34.1 reference-count-overflow class cannot be screened statically); this is not a validation of the proof, since compilation, the elaborated-statement comparison, the axiom closure and an independent kernel replay (with a checker that rejects the #15908 exploit) remain open.

# Alpöge RH Lean archive: static trust-boundary audit (alpaudit, 20 Nov 2026, 17:24–18:35 UTC)

Input: `lean.tar.gz` sha256 9631b4f0…c7fc and `remark.pdf` sha256 3fc97784…eb8d (re-hashed here; equal to Hertaspace's SHA256SUMS). Toolchain `leanprover/lean4:v4.33.1`; Mathlib `0df444a360ea`. No Lean, no lake (300 MB cap). Scope: tasks (1) bypass hunt, (2) statement fidelity, (3) kernel-risk profile, (6) external replays. Tasks (4) falsification probes and (5) remark→module map moved to alpmath (Riemann/Deep/AlpMath/) by orchestrator note.

## Short answer

If the build succeeds, the five theorems of `src/Main.lean` mean what they say, and nothing in the source can make the kernel accept less than a proof, as far as a complete static scan of all 175,187 files can tell: no unchecked-declaration API, no `debug.skipKernelTC`, no private Lean-core API opened, no global notation or term macro, no root alias, no redefinition of core names a known v4.33.1 kernel bug needs. Whether the proofs typecheck is NOT established here: that needs the build (Hertaspace: 24,167/183,903 jobs at 18:10 UTC) and the axiom check and kernel replays.

## Method

- Import cone: rg dump of all 909,822 `import` lines, BFS from Main (code/cone2.py): the cone is **every one of the 175,187 project modules** (incl. BuildAid's 5) plus 348 external modules (344 Mathlib, 3 Lean: `Lean`, `Lean.Elab.Command`, `Lean.LibrarySuggestions.Basic`; 1 Batteries: `Batteries.Tactic.OpenPrivate`). So every finding below is on the cone.
- Seven rg -j1 passes over src (import dump, scan2-scan7) (ASCII-mode regexes, memcap peaks ≤ 5 MB), classified by code/classify.py and small scripts; counts in logs/classify.log. All meta code was then read by hand: the 4 `elab`s, 5 `run_cmd`s, 32 `run_tac`s, 7 files importing Lean directly.
- Lean v4.33.1 sources read: src/Lean/Environment.lean, AddDecl.lean, CoreM.lean, Replay.lean, Elab/Declaration.lean, Elab/MutualDef.lean, src/kernel/inductive.cpp.
- Pinned Mathlib read: RiemannZeta.lean (byte-identical to .scratch/alp/RiemannZeta_pin.lean), DirichletContinuation.lean, DirichletCharacter/Basic.lean, Nonvanishing.lean, lake-manifest.json.
- Archive layout: only 3 non-.lean files (lakefile.toml, lake-manifest.json, lean-toolchain), no hidden files, no src/Mathlib, src/Lean, src/Init or src/Batteries shadow trees; every inherited package rev in lake-manifest.json equals Mathlib's own manifest at the pin (batteries 4488d40d, aesop 3448c0bc, Qq 92c15be1, plausible b7eb3304, LeanSearchClient 5f4d51b8, importGraph 16f02aa7, proofwidgets 4be2e3d5, Cli 6130a478).

## (2) Statement fidelity (pinned Mathlib 0df444a3)

Mathlib (quotes ≤ 5 lines each; RiemannZeta.lean = .scratch/alp/RiemannZeta_pin.lean):
- RiemannZeta.lean:121 `def riemannZeta := hurwitzZetaEven 0` (analytic continuation; Dirichlet series only for Re s > 1, :207 `zeta_eq_tsum_one_div_nat_cpow`).
- RiemannZeta.lean:185-187 `def RiemannHypothesis : Prop := ∀ (s : ℂ) (_ : riemannZeta s = 0) (_ : ¬∃ n : ℕ, s = -2 * (n + 1)) (_ : s ≠ 1), s.re = 1 / 2`; trivial zeros :173 `riemannZeta_neg_two_mul_nat_add_one (n : ℕ) : riemannZeta (-2 * (n + 1)) = 0` (same expression, so exactly the trivial zeros are excluded).
- DirichletContinuation.lean:61 `noncomputable def LFunction (χ : DirichletCharacter ℂ N) (s : ℂ) : ℂ := ZMod.LFunction χ s` (`[NeZero N]` in scope; docstring :57-58: agrees with `LSeries χ` for Re s > 1); :67-69 `LFunction_modOne_eq {χ : DirichletCharacter ℂ 1} : LFunction χ = riemannZeta`; :92 `differentiableAt_LFunction … (hs : s ≠ 1 ∨ χ ≠ 1)` (pole only for χ = 1 at s = 1).
- DirichletCharacter/Basic.lean:246 `conductor χ := sInf (conductorSet χ)`, :292 `def IsPrimitive : Prop := conductor χ = n`.

Junk values: ζ(1) and L(1_N, 1) are junk in Mathlib; every statement excludes them (s ≠ 1, `χ ≠ 1 ∨ s ≠ 1`, or Re s < 1). N = 0 (DirichletCharacter ℂ 0 = characters of ℤˣ, LFunction not meaningful) is excluded by `[NeZero N]`.

| Main theorem (Main.lean line) | type | intended theorem | faithful? |
|---|---|---|---|
| riemannHypothesis (4) | Mathlib `RiemannHypothesis` | RH | yes: Mathlib's own definition, root constant, no capture (TABLE c1-c7) |
| riemannZeta_eq_zero_re (6) | `RiemannZetaZerosOnLine` (Statements.lean:3-4) | RH | yes: identical to RiemannHypothesis up to binder names; `-2 * ((n : ℂ) + 1)` is what Mathlib's `-2 * (n + 1)` elaborates to (ℕ cast at the leaf, binop%); the proof term is `rhOfAll : RiemannHypothesis`, so the kernel itself must find the two defeq |
| LFunction_ne_zero_of_half_lt_re (8) | `DirichletLNonvanishingRightOfHalf` (6-8) | GRH for all Dirichlet L-functions, zero-free form Re s > 1/2 | yes: for imprimitive χ, L(χ,s) = L(χ*,s)∏_{p∣N}(1 − χ*(p)p^{−s}) and the factors vanish only on Re s = 0, so this is GRH for χ*; N = 1 gives ζ(s) ≠ 0 for Re s > 1/2, s ≠ 1 (LFunction_modOne_eq), which with the functional equation is RH; only the pole (χ = 1, s = 1) is excluded |
| LFunction_eq_zero_re (10) | `DirichletLZerosOnLine` (10-12) | GRH, strip form, all χ | yes: strip 0 < Re s < 1 removes the Euler-factor zeros (Re s = 0) and the trivial zeros (Re s ≤ 0) of imprimitive χ; nontrivial zeros lie in the strip (no zeros on Re s = 1, Mathlib Nonvanishing.lean `LFunction_ne_zero_of_one_le_re`), so strip form = GRH |
| LFunction_eq_zero_re_of_isPrimitive (12) | `PrimitiveDirichletLZerosOnLine` (14-16) | GRH, strip form, primitive χ | yes (weaker than line 10) |

Half-plane vs strip: line 8 is stronger than RH-on-strip by itself (it also asserts non-vanishing on 1/2 < Re s, including Re s ≥ 1 where it is Mathlib's theorem); both are the standard GRH. q = 1: the only character mod 1 is 1 and LFunction_modOne_eq turns lines 8/10 into statements about ζ.


## (1) Bypass hunt (details: TABLE.md a-d)

- Unchecked-declaration paths of Lean v4.33.1 (TABLE b1-b16): every public one has 0 archive hits; every private one would need `open private … from Lean.*`, and all 220 `open private` clauses target project modules or 3 Mathlib private lemmas (Mathlib.LinearAlgebra.SymplecticGroup). No `debug.*` option, no `#eval`, no `initialize`, no `unsafe`, no `implemented_by`/`extern`, no file IO.
- All meta code in the archive was enumerated and read: Imports.lean:5-8 (`kexactB`), OwnStruct.lean:7-64 (`kcopy`, `#kcopy_probe`), CutCptCopy.lean:14-91 (`cpt_copy_thms`, used once, CutMainCpt.lean:945), BuildAid/D1-D5.lean:4 (`run_cmd` deny-list entries), 32 `run_tac` blocks (goal-level). The elaborator's type check is skipped by `kexactB`/`kcopy` (unchecked `MVarId.assign`), but never the kernel's: an ill-typed assignment makes the declaration fail at `addDecl`. `cpt_copy_thms` adds renamed copies with `Lean.addDecl` (checked).
- 93 syntax-level declarations (91 macro/elab + 2 unif_hint): 26 global (22 tactic macros, 4 elabs above), 11 scoped, 56 local. No global notation, infix, term macro, `syntax`, `macro_rules` or `elab_rules`; no `@[macro]`/`@[term_elab]`/`@[command_elab]`/`@[tactic]` attribute. None can re-parse Main.lean (ASCII, no `open`, no tactic block).
- Name capture: impossible at the root (duplicate root constants fail at import; root aliases: none); all 36 project decls named re/IsPrimitive/LFunction/Complex.* are namespaced; MainTheorems.lean's `open Complex DirichletCharacter` cannot weaken Main because Main restates each type and the kernel checks it.

## (3) Kernel-risk profile of the cone (= all 175,187 modules)

| feature | count | where / examples | risk |
|---|---|---|---|
| `decide` (kernel evaluates a Decidable instance) | 6,926 lines, 2,465 files (5,159 `by decide`) | densest: Formalism/LocalRingCardinality/SpanPowerMembership/Primes.lean (416 lines, ~790 uses: Lucas primality certificates via windowed `a ^ e % n = r` facts, :12 `lucas_list`, :26 `pow_mod_of_two_pow_window`); OverFacts.lean 87; 0 decide lines over String/Char | relies on the kernel's GMP-accelerated Nat literals (`Nat.pow`, `Nat.mod`, `Nat.mul`, `Nat.ble`, …); see #15908 below |
| numerals ≥ 10^19 (≥ 20 digits; superset of ≥ 2^64) | 639 in 4 files | Primes.lean 627; ExampleLayerExistsGenerators.lean 6 (longest 60 digits, :228 `levelPrimes`, :264 `T0`, :271 `e0`); CyclotomicIsCyclotomicExtension.lean 4 | GMP literal path; independent check of these certificates is alpmath's (4) |
| String / Char literals | 1,905 / ≤ 57 lines | notation/macro atoms, error messages in meta code | not kernel-evaluated (no decide over String/Char) |
| unchecked goal assignment (`kexactB`, term macros `kc%` `kq%` `kd%` `copyDec%` `copyIsV%`, `kcopy`) | kexactB 278, kc% 400, kcopy 9, copyDec% 4, in 72 files | e.g. PlaceModelImp.lean:27 `kc% e => (by kexactB e)`; GoodPartEq.lean:17 `tr_via` | kernel is the only checker of these subterms: more kernel defeq work, no soundness loss |
| `cpt_copy_thms` (addDecl of renamed copies) | 1 use (CutMainCpt.lean:945) | CutCptCopy.lean:68-82 | kernel-checked |
| nested inductives | 0 of 120 `inductive` decls (heuristic: no constructor argument applies List/Array/Option/Multiset/Finset/Vector/Prod/Sum/Subtype/Sigma/Set/Thunk to the type itself) | — | #16490/#16491 not exercised |
| mutual inductives | 1 block | PrincipalizationCertificate.lean:10 | standard |
| `_nested` references | 83 substrings, 0 tokens beginning `_nested` | e.g. `hasGen_nested` (ExistsPrimitiveNormalisation.lean:123) | none |
| direct `Quot` use | 2,606 lines: Quot.out 1,822, Quot.mk 751, Quot.out_eq 416, Quot.sound 160, Quot.lift 55, Quot.ind 38 | Mathlib API | standard (Quot.sound is one of the 3 standard axioms) |
| well-founded recursion | termination_by 20, decreasing_by 19 (27 files); `partial` 2 (meta code only) | — | WF defs are irreducible by default; whether a decide/rfl forces the kernel through WellFounded.fix was not checked |
| hand-built `Declaration`s | `defnDecl`/`thmDecl` 1 each (CutCptCopy), `inductDecl`/`axiomDecl`/`opaqueDecl` 0 | — | no crafted inductives (the #15704/#16490 entry point) |
| `prelude` modules / redefinition of core names (Nat.add…Nat.log2, outParam, semiOutParam, optParam, autoParam, String.mk, Char.ofNat, Eq, Quot, Nat) | 0 (scan5) | every module imports Mathlib, so Init's definitions are always present | preconditions of #15704 and #15908 absent |
| `--tstack=400000` (lakefile.toml:4) | — | worker-thread stack size (KiB) | crash avoidance for deep recursion on large (kexactB/kcopy-built) terms; no soundness effect |

### Known Lean kernel bugs vs v4.33.1 (leanprover/lean4 issues, read 18:10–18:20 UTC via gh api)

| issue | state (UTC) | what | in v4.33.1? | precondition | archive |
|---|---|---|---|---|---|
| #15704 kernel soundness: `consume_type_annotations` + `elim_only_at_universe_zero` | closed completed 2026-10-17 | large elimination wrongly granted to a Prop inductive when binder domains use outParam/semiOutParam/optParam/autoParam whose definitions differ from their syntactic meaning | yes (v4.33.1 inductive.cpp:222, 227 call `consume_type_annotations`) | nonstandard definitions of those four names (fresh env / `prelude`) | 0 `prelude`, 0 redefinitions, Init present everywhere: not reachable |
| #15908 axiom-free proof of False: native `Nat.add` reduction trusts the name | closed completed 2026-10-24 | crafted export file with a user `Nat.add`; Lean default + paranoid kernels and nanoda accept False; lean4lean, con-leche, con-ron reject | yes | a constant named `Nat.add` (etc.) that is not Init's: export-file / `prelude` only | not reachable from these sources; matters for export-based replays (see (6)) |
| #16490 `_nested_1` auxiliary name referenceable | open, filed 2026-11-16 | anonymous/macro-scoped nested type ⇒ aux name escapes the `_nested` prefix check; reporter: no inconsistency derived | yes (inductive.cpp:1105 `mk_unique_name(*g_nested + J_name)`, :1222-1231 prefix test) | nested inductive with anonymous or macro-scoped name + user reference to the aux name | 0 nested inductives, 0 `_nested*` tokens, 0 `inductDecl`: not exercised |
| #16491 nested constructor loses binder info/names | open, filed 2026-11-16 | alpha-equivalent nested occurrences restored to the first one; reporter: not unsound | yes (inductive.cpp:1087 `if (p.first == Iparams)`) | two alpha-equivalent nested occurrences | 0 nested inductives; alpha-only, cannot change Main's types |
| #16489 `add_inductive` accepts a constructor type with a loose bound variable, stores a different type | closed 2026-11-17 | kernel inductive check | yes (fixed one day after 4.36.0, per Deicyde NOTES.md:103) | hand-built `Declaration` (the elaborator never emits loose bvars) | 0 `inductDecl`: not reachable |
| PRs #15241, #15289, #15288 (4.34.1 runtime: frozen/overflowed reference counts, `sharecommon_quick_fn`) | merged 2026-09-23/24 | > INT_MAX references to one shared subterm can become a use-after-free; 4.34.1 notes: potentially a proof of False with deliberately constructed input (≥ 18 GB RAM) | yes (not backported to 4.33.1, Deicyde NOTES.md:100-103) | a term with > 2^31 references to one subterm | not excluded statically, but implausible for modules that build in ~2 min each; kernels not built on the Lean runtime (nanoda, con-*, lean4lean) are unaffected |
| #15100, #15101 `leanchecker` accepts mismatched constructor keys / arbitrary quotInfo entries | closed 2026-09-09 | export-replay checker bugs | n/a to sources | crafted export | relevant only to replays |
| #8554 soundness bug in {Expr, Level}.data | closed 2025-05-30 | — | fixed long before v4.33.1 | — | — |

Title search `soundness` (25 hits) plus label:bug `kernel False` since 2025-06 (15 hits); not exhaustive.

## (6) External builds and replays (read-only, gh api)

| project | state | UTC |
|---|---|---|
| Hertaspace (update) | `latest_job 25654/183903`, built_ok 9,511, with warnings 8,468, failed 1 (the .ir cache error); `verification/comparator.json`: five theorem names, permitted_axioms [propext, Quot.sound, Classical.choice], `enable_nanoda: false` (comparator uses Lean's own kernel only) | status 18:21 |
| Deicyde/riemann-hypothesis (Jack McCarthy) | byte-exact mirror + Comparator under Lean v4.36.0 with external kernels leanchecker-paranoid, nanoda, con-leche, con-ron, plus lean4lean master 67f63b0c (built with v4.33.1) as sixth checker; smoke test on a Mathlib-sized export: five accept, lean4lean `--paranoid` rejects (`unknown constant 'Decidable.intro'`, a checker limitation) (verification/NOTES.md:115-139); notes that a pass by the 4.33.1 kernel alone (= a successful build) is limited evidence because of the unported 4.34.1 runtime fixes | pushed 18:01:42 |
| shosonoda (update) | build on RIKEN's RIKYU (aarch64, 144 cores/node), then a Comparator rebuild; README: 753,569 internal imports, acyclic, depth 202, every module imported by Main. **Our cross-check (code/graph_stats.py on our own import dump): 753,569 internal imports (753,261 distinct), acyclic, Main the only unimported module, longest chain 203 modules = depth 202: identical** | README read 18:31 |
| Lean Zulip #general "Discussion: Alpoge Proves the Riemann Hypothesis" | login needed; not read | — |
| Hertaspace/alpoge-rh-lean-verification | build `latest_job 24167/183903`, built_ok 8,678, with warnings 7,814, `failed 1` (Statements.lean:1:0 could not read Mathlib `…/Bicategory/Extension.ir`, "invalid header": an I/O error in the Mathlib cache, Statements compiles standalone per their VERIFICATION.md:72-78); axiom check empty (pending); comparator (tag v4.33.0) + clean-room statement check + environment scans + replay by Lean paranoid/default kernels, nanoda 0.4.17, con-leche, con-ron (bundled with Lean 4.35.0-rc2, export format 3.1.0) planned (VERIFICATION.md:86-102) | status 18:10; last commit 18:10:49 |
| brunoloff/rh-lean-progress | progress.json: phase build, files_done 18,744/175,187, lines 2.94 M/17.37 M, failed 0, axioms null, replay null, ETA 2026-11-21 04:47 | updated 18:20:41 |
| rohanarun/alpoge-rh-replay | crowd-sourced per-module replay (175,187 module tasks, 2,858 packets); "Full replay in progress. No full RH verification receipt" | last commit 18:24:27 |
| shosonoda/lean-rh-alpoge | independent build-report repo; README describes the method (Mathlib from cache; comparator sandbox, axiom check, kernel replay); no result in the README | last commit 17:18:04 |
| TomasOrtega/RH, PyJudge/…, 1337isnot1337/…, ShashiQubit/…, Reen-Ai-H/…, adamrsophie/RiemannHypothesis, harttraveller/riemann, realRoc/… | mirrors, guides, explainers; no checks | 17:19–18:16 |

Replay caveat from #15908: Lean's own default and paranoid kernels and nanoda accepted a crafted export proving False; for an export-based replay to be independent evidence, at least one of lean4lean / con-leche / con-ron (which reject it) should accept the export, or the Lean kernels should include the 24 Oct fix (4.35.0-rc2 should; not verified here). Correction (R538): no such fix exists: #15908 was closed as not a kernel bug, and Comparator's exact comparison of Nat.add and the other primitives rejects such exports.

## Sorry-producing tactics and elaboration entry points

- scan7 (all files): 0 `sorry`, `admit`, `sorryAx`, `exact?`, `apply?`, `rw?`, `aesop?`, `polyrith`, `says`, `extract_goal`, `proof_wanted`; `stop` 4 hits, none the tactic (field `stop : κ n = tau` LemmaDSkeletonUnit.lean:280, projections `hp.stop`/`H.stop`); `hint` 22,623 hits, hypothesis names except ≤ 6 bare-line uses, and Mathlib's `hint` closes a goal with a real proof or fails; `unseal` 39 (elaborator transparency only). Error-recovery sorries come with an error, which fails the build; an explicit `sorry` would only warn (7,814 modules built with warnings at 18:10), hence the axiom closure check stays required.
- How executable elaboration code was found: (i) every command that runs code (`#eval`, `run_cmd`, `run_elab`, `run_meta`, `run_tac`, `initialize`, `elab`, `elab_rules`, `macro_rules`, `syntax`, `macro`, attributes `@[macro|term_elab|command_elab|tactic|*_parser|init|env_linter|csimp|implemented_by|extern]`) was scanned and every hit read; (ii) every function typed in a meta monad (MetaM, CoreM, TacticM, CommandElabM, TermElabM, AttrM, ImportM, EIO, BaseIO): only OwnStruct.lean and CutCptCopy.lean; (iii) direct imports of Lean core meta modules: 7 files (Imports.lean, OwnStruct.lean, BuildAid/D1-D5). Untyped `do` code would evade (ii) but must be invoked from an entry point in (i), all of which were read. Mathlib's and core's own elaborators, tactics, attributes and deriving handlers are trusted scope.

## Independent reviewer (subagent, 18:27 UTC) and fixes

Verdict: "qualified static-screen pass; proof validation incomplete". Points and fixes: (1) narrow the headline → verdict line rewritten; (2) implicit sorry not excluded by lexical scans → scan7 added above, axiom closure listed as open; (3) compare elaborated statements → open (comparator / clean-room check); (4) unchecked `assign` and checked `addDecl` are not bypasses → agreed, as stated; (5) artifact provenance → trusted scope below; (6) known-bug exclusions overstated → reworded as "checked known trigger patterns"; aliases are not constants (the #15908 GMP path keys on constant names), no command macro or elab generates an inductive, the nested-inductive test is textual; (7) lexical negatives are not an allowlist → entry-point enumeration documented above; (8) build counters prove progress only, replay needs a checker shown to reject the #15908 exploit → agreed.

## Trusted, not checked here

Lean v4.33.1 binaries and kernel; Mathlib 0df444a3 and its 8 pinned dependencies (and any prebuilt cache, e.g. Hertaspace's `lake exe cache get`); the build machines. A static audit cannot authenticate artifacts, flags or search paths of a particular build.

## Open items

1. The build of all 183,903 jobs (Hertaspace 24,167 at 18:10 UTC; brunoloff 18,744/175,187 files at 18:20 UTC).
2. Axiom closure of the five theorems ⊆ {propext, Classical.choice, Quot.sound}.
3. Comparator and clean-room check of the elaborated statements (incl. binders and instances in Statements.lean).
4. Independent kernel replay including lean4lean, con-leche or con-ron (#15908; and the 4.34.1 runtime fixes are absent from 4.33.1): Hertaspace's comparator.json has `enable_nanoda: false`, so its comparator run alone is Lean-kernel evidence; Deicyde's setup runs nanoda, con-leche, con-ron and lean4lean. The 4.35.0-rc2 / 4.36.0 kernels' inclusion of the 17 and 24 Oct fixes not verified here.
5. Residual static gaps: nested-inductive test is a textual heuristic over `inductive` declarations (recursive `structure`s not scanned); whether a `decide`/`rfl` makes the kernel unfold well-founded recursion not checked.
6. Tasks (4) falsification probes (incl. the Lucas certificates in Primes.lean and the 60-digit constants) and (5) remark→module map: with alpmath.

## Files

- `Riemann/Deep/AlpAudit/REPORT.md` (this), `TABLE.md` (rows a1-a93, b1-b16, c1-c7, d1-d14).
- `code/cone.py` (first cone attempt, killed: O(n²) bug), `code/cone2.py` (import cone from rg dump), `code/graph_stats.py` (internal imports, acyclicity, depth; logs/graph_stats.log, peak 251 MB), `code/classify.py` (scan classifier), `code/table_syntax.py` (TABLE section a).
- `logs/cone.log` (cone counts), `logs/cone_external_imports.txt` (348 external modules), `logs/classify.log` (counts), `logs/syntax93.tsv` (file:line, scope, atom).
- `code/aux/` small checks: an4/an5 (export namespaces, nested inductives), an6/an10 (name capture + namespace walk), an7 (Quot, numerals, decide), an11 (manifest vs Mathlib pin, non-ASCII in statement files); they read the scratch dumps.
- The 16 MB list of the 175,187 cone modules was in `.scratch/alpaudit/cone_files.txt` (not committed; deleted after the audit, regenerable by code/cone2.py).
- Scratch (not in repo, deleted after the audit): `.scratch/alpaudit/` scan dumps scan2–scan7, rg pattern files, Lean/Mathlib source copies (`web/`), issue bodies, external READMEs (`ext/`), STATE.md.
