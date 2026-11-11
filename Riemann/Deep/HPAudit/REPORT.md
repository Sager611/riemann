AUDIT PASSED: (A) ζ(s) ≠ 0 and L(s, χ) ≠ 0 for Re s > 7/8, every Dirichlet χ (OpenAI); (B) L(s, χ) ≠ 0 for Re s > 0.87338196661320752 and for Re s > B_2 + 2.16·10⁻²⁵⁴, every modulus and χ (Bebel, tracker PR 101, PR 103); (C) the same for Re s > B_2 = (777 − √1081)/852 exactly (Richey, PR 129; includes ζ as q = 1); (D) ζ_K(s) ≠ 0 and finite-order Hecke L ≠ 0 for Re s > 31/32, every number field K (Cureton–Thorner, AxiomMath@09afbb4). Each Lean statement says what we cite, in Mathlib's own definitions. Each proof passed a Comparator statement match against a challenge built from Mathlib alone, used only propext, Classical.choice and Quot.sound, and was accepted by three independently written kernels: Lean, nanoda and con-ron. Those checks ran in third-party CI; we read their logs and matched hashes but did not re-run a kernel ourselves.

# HPAudit: referee report on the four Lean-checked zero-free half-planes (node N7), 11 Nov 2026

Agent hpaudit (Orca). Read-only GitHub/web access (gh api GET, curl -o); nobody was contacted. Command outputs are in logs/. Clones and downloads were in .scratch/hpaudit/ and were deleted afterwards; STATE.md is kept, and the scripts are in tools/.
Method, per input: (1) locate the code and print the exact statement; trace every name to its definition. (2) Look for axioms and red flags: scan every source file with comments and strings stripped (logs/scan-summary.txt). (3) Read the CI evidence and check its hashes against the sources. (4) Compare the statement with the literature and with known barriers. Part 5 is the dependency map. No local kernel replay was run. The brief requires one only if the CI evidence is missing or unclear, and here it is complete for all four inputs. A replay would also not fit the limits: the solution exports are 2.2 GB, and nanoda alone peaked at 6.3 GB on the 7/8 proof (tracker receipt pr69 README).

## Definitions behind all four statements (logs/mathlib-defs.txt)
- **DirichletCharacter.LFunction** (Mathlib, LSeries/DirichletContinuation.lean @d13f23b7). It is defined as ZMod.LFunction, the Hurwitz-zeta analytic continuation. It equals the Dirichlet series for Re s > 1 (LFunction_eq_LSeries), is differentiable everywhere except s = 1 for χ = 1 (differentiableAt_LFunction), and LFunction_modOne_eq gives LFunction = riemannZeta for q = 1. By the identity theorem it is the true L(s, χ) on ℂ∖{1}.
  - An imprimitive χ gives the imprimitive L-function. Its extra Euler factors vanish only on Re s = 0, so it has the same zeros in Re s > 0 as the primitive one.
  - riemannZeta is Mathlib's. No 'def riemannZeta' exists anywhere in Riemann/ or in the tracker sources.
- **NumberField.dedekindZeta** (Mathlib @0493595e, l.49) is defined as an LSeries. It returns the junk value 0 wherever the series diverges, i.e. for Re s ≤ 1. A naive statement 'dedekindZeta K s ≠ 0 for Re s > 31/32' would therefore be false in Lean. Statement (D) avoids this trap (see D).

## (A) OpenAI, Re s > 7/8 (load-bearing via Riemann/Deep/QRH)
- **Code.** openai/math @1581d4cb, lean/OAI/NumberTheory/DirichletL/Nonvanishing.lean. 1581d4cb is the tracker's pin fd4aeeb2 plus one workflow file (.github/workflows/lock-prs.yml), so the Lean is identical. Our port is Riemann/Deep/QRH (2944 files) and the bridge is Riemann/Deep/QRH/Bridge.lean.
- **Exact statements.**
  - theorem OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re {s : ℂ} (hs : (7 / 8 : ℝ) < s.re) : riemannZeta s ≠ 0
  - theorem OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re {q : ℕ} [NeZero q] (χ : DirichletCharacter ℂ q) {s : ℂ} (hs : (7 / 8 : ℝ) < s.re) (hpole : ¬ (χ = 1 ∧ s = 1)) : _root_.DirichletCharacter.LFunction χ s ≠ 0
  - Ours: Fayed.QRH.quasiRH_seven_eighths : Fayed.R7.QuasiRH (1 / 8), where QuasiRH δ := ∀ s : ℂ, 1 − δ < s.re → s ≠ 1 → riemannZeta s ≠ 0 (R007 l.25, namespace Fayed.R7).
- **Definition trace.** All names come from Mathlib. The hypotheses are satisfiable and not vacuous. The statements cover all heights, every q ≥ 1 and every χ, primitive or not. The ζ theorem even holds at s = 1, because Mathlib sets riemannZeta 1 = (γ − log 4π)/2 ≠ 0.
- **Axioms.**
  - Our QRH/REPORT.md records [propext, Classical.choice, Quot.sound] for the four declarations. Its check file (.scratch/qrh2/ax.lean) has since been deleted.
  - Three external runs confirm the same three axioms: the tracker receipt pr69 Comparator run, Reilly (Zenodo 23505836, 34 #print axioms) and Daugherty–Ryan–Ward (Zenodo 23461828).
- **Red flags.** The scan of all 2944 port files finds no sorry, admit, axiom, native_decide, ofReduceBool, implemented_by, extern, unsafe, opaque, debug option, run_cmd or elab. It finds 4 'decide +kernel' proofs (kernel-evaluated, no extra axiom) and one local 'def LFunction' for Hecke L-functions of ℚ(√−3) that no statement we use touches. Verdict: nothing can affect the kernel-checked statement.
- **Independent checks** (logs/A-independent-replays.txt):
  - Tracker receipt pr69-ac76d5ef (23 Oct): the OpenAI baseline at 7/8, all moduli, passed Comparator, con-ron (122,039 declarations), nanoda and Lean.
  - Reilly, Zenodo 23505836: Comparator plus nanoda (3a24072) over all 189 Comparator configurations of openai/math; QuasiRiemannHypothesis ran at adc7f12, whose import closure is unchanged. The challenges QuasiRiemannHypothesis and DirichletSevenEighths match Thm 1.1 of the 30 Sep manuscript.
  - Daugherty–Ryan–Ward, Zenodo 23461828: rebuilt all 2926 modules in a separate tree and replayed them in Lean's kernel in a fresh environment. Their statement audit against Mathlib definitions agrees with ours.
  - Other replays: nealmcb (Comparator on all 416 configurations) and davegoldblatt (both kernels).
  - Our port matches the triple-checked tracker baseline: 2901 of 2924 files are identical after normalising headers and imports. The other 23 differ only by the documented proof-script edits (19 set-literal parentheses, 3 congr → all_goals, 1 ULift.ext). See logs/A-port-vs-tracker-baseline.txt.
- **Literature.**
  - OpenAI's manuscript of 30 Sep 2026, "The Quasi-Riemann Hypothesis: A Zero-Free Half-Plane Re s > 7/8", Thm 1.1 covers all Dirichlet L and the Hecke L of ℚ(√−3). The Lean statement covers the Dirichlet and ζ parts we use.
  - The text is AI-generated. The openai/math README says "the writeup for the Re(s) > 11/12 zero-free region ... was human edited". Liu (2610.12234, 2610.26835) and our own audits re-derived the argument with varied parameters.
  - No doubt or erratum anywhere. openai/math history.md withdraws only the Weil-class papers. The arXiv and Zenodo searches and the tracker, openai/math and AxiomMath issues show nothing against this result. The closest is GENESIS R300.1 (Zenodo 23482824, an AI-team review): "strongly supported but not conclusively verified", because that team did not reproduce the chain itself.
- **Barriers.** 7/8 lies above the method floors 13/15 (R109) and 167/192 (Liu's free-M floor, tracker issue 87), so it is consistent with both.
- **Verdict: holds.**

## (B) Bebel: Re s > 0.87338196661320752 (PR 101) and Re s > T = B_2 + 2.16·10⁻²⁵⁴ (PR 103)
- **Code.** github.com/GTimG/QuasiRiemannTracker:
  - PR 101: submissions/2026-10-31-positive-gap-contact-return, head 20a4c9a5, merged 1 Nov 09:53Z.
  - PR 103: 2026-11-01-closest-fraction-contact, head cae2e523, merged 1 Nov 12:23Z.
  - Each Candidate.bound imports 4 direct and 3 transitive earlier submissions: openai-baseline, hailey-collet-1252270103, arvi-repeated-marks-on-sil-face, arvi-selected-leaf, sil-density-split, density-split-fine-chain and arvi-density.
  - The src trees on main @375bf5f equal the admitted trees (PR 103's src tree 884a35c3 = record source_tree). No commit has touched these src dirs since admission.
- **Exact statement.** theorem Candidate.bound {q : ℕ} [NeZero q] (χ : _root_.DirichletCharacter ℂ q) {s : ℂ} (hs : (θ : ℝ) < s.re) (hpole : ¬ (χ = 1 ∧ s = 1)) : _root_.DirichletCharacter.LFunction χ s ≠ 0.
  - PR 101: θ = 5458637291332547/6250000000000000, which is B_2 + 7.10·10⁻²⁰.
  - PR 103: θ = 2804966284…654341/3211614609…872271 (127/127 digits), which is B_2 + 2.16·10⁻²⁵⁴.
  - The tracker checks it as QRH.canonical_bound. Its challenge comes from verifier/upstream-challenge.lean (sha pinned in pins.json) with the bound rendered from integers.
  - Our regeneration of PR 103's Challenge.lean (sha e188be23…) equals the CI's challenge-src hash.
- **Definition trace.** As for (A): Mathlib only, every modulus and every character, all heights. Any shadowing (for example a candidate's QRH.DirichletCharacter) would change the elaborated statement, and Comparator would then reject it.
- **Axioms and red flags.** The CI's permitted list is the standard three. The CI compile logs print '#print axioms Candidate.bound' as [propext, Classical.choice, Quot.sound], and none of the 1650 (PR 101) or 1651 (PR 103) per-module logs mentions sorry; all 492 #print axioms lines in each run are the standard three (logs/B-pr10{1,3}-print-axioms.txt). The scan of the 473 own files plus 8,531 dependency files finds no sorry, admit, axiom, native_decide, ofReduceBool, extern or debug option. The other hits are explained in scan-summary.txt and are benign (decide +kernel, MetaM certificate tactics, a dependency audit via run_cmd, one custom theorem command through Lean.addDecl, which runs the kernel).
- **Independent checks** (logs/B-*, BC-*): the PR-admission workflow ran on self-hosted runners.
  - PR 101: run 41175022316, accepted 31 Oct 22:45Z; judge.log shows con-ron accepted 135,905 declarations, then nanoda, then Lean, then "Your solution is okay!".
  - PR 103: run 41294260719, accepted 1 Nov 11:03Z; con-ron 135,931 declarations, then nanoda and Lean.
  - Both runs passed the positive, mismatched-statement and ill-typed controls.
  - The 2923 Mathlib/OpenAI modules restored from a signed build cache do not weaken this: the kernels replay the whole exported environment, not the cache.
  - PR 127 (darkhipo, AI): an outside rebuild of PR 103 on macOS with the default kernel found the standard three axioms (offered and declined).
- **Literature.**
  - The claims in PR 101 and PR 103 match the Lean statements.
  - B_2 is the contact point of the Dirichlet comparison scheme, written up by Liu, Hexagon 2611.00082. We audited that note in Hexagon2: sound given 2610.26835v2.
  - No public doubt.
- **Barriers.** Consistent: B_2 > 167/192 > 13/15.
- **Verdict: holds.**

## (C) Richey: Re s > B_2 = (777 − √1081)/852 exactly (PR 129)
- **Code.** submissions/2026-11-10-exact-contact-point, 215 files, head 515dbe0c, merge commit 535397cd (11 Nov 11:13Z, by nielstron). It has the same dependencies as (B), and the src tree 23015664 equals the record's.
- **Exact statements.**
  - Candidate.bound as in (B), with hs : ((777 − Real.sqrt 1081) / 852 : ℝ) < s.re.
  - Candidate.bound_le : (777 − √1081)/852 ≤ T, where T is PR 103's fraction. It only ranks the entry on the leaderboard; we computed T − B_2 = 2.16·10⁻²⁵⁴ > 0.
  - Both are checked as QRH.canonical_bound and QRH.canonical_bound_le. We regenerated Challenge.lean (sha c89bb4d3…) and Solution.lean (7fc4ffd7…); both equal the run's contract (logs/C-pr129-Challenge-Solution.lean.txt). With q = 1 the statement covers ζ.
- **Axioms and red flags.** Same as (B). The CI log of Candidate.Main gives Candidate.bound and Candidate.bound_le as [propext, Classical.choice, Quot.sound]. None of the 1629 compile logs mentions sorry, and all 485 #print axioms lines are the standard three (logs/C-pr129-print-axioms.txt).
- **Independent check.** Run 42995393841, accepted 11 Nov 10:45Z. judge.log shows con-ron accepted 135,951 declarations, then nanoda, then Lean, then "Your solution is okay!"; the controls behaved (logs/C-pr129-*). The record is on branch admission-records (report sha 34c32c4c…).
- **Literature.**
  - Liu's 2611.00082 claims B_2 for every finite-order Hecke L over every number field.
  - The Lean result covers only Hecke L over ℚ, i.e. ζ and Dirichlet L. Richey says (PR 129 comment) that he converted his proof from a stronger formalisation for unitary Hecke L. That formalisation is unpublished and unchecked.
  - Our pages say "Lean-checked" only for ζ and Dirichlet L; for Hecke L over other fields they cite Liu as a human proof. This is right: no mislabel was found (search over TREE, COROLLARIES, RECORDS, DISCOVERIES and discoveries/).
- **Barriers.** It reaches the contact B_2 and claims nothing below it. This is consistent with BeyondB2 (going below B_2 needs a phase-sensitive estimate).
- **Verdict: holds.**

## (D) Cureton–Thorner: Re s > 31/32 for every number field
- **Code.**
  - github.com/AxiomMath/DedekindZetaZeroFree @09afbb4 (13 Oct 2026, Lean 4.35.0-rc3, Mathlib 0493595e): 5528 .lean files. The proof closure is 5509 modules, including 1658 vendored class-field-theory modules; outside the repository it imports only Mathlib, not Tau Ceti.
  - Archived as PalomarArchive/AxiomMath--DedekindZetaZeroFree--42e46c822b86, tag palomar/PALOMAR-2026-10-13-000010-v2/09afbb4; also Zenodo 23373237 (R399).
- **Exact statements** (Challenge/DedekindZetaChallenge.lean):
  - theorem DedekindZetaZeroFree.main (K : Type) [Field K] [NumberField K] : ∃ f : ℂ → ℂ, DifferentiableOn ℂ f ({1}ᶜ : Set ℂ) ∧ (∀ s, 1 < s.re → f s = NumberField.dedekindZeta K s) ∧ ∀ s, 31 / 32 < s.re → s ≠ 1 → f s ≠ 0
  - theorem DedekindZetaZeroFree.heckeCorollary (K) [Field K] [NumberField K] (m : Ideal (𝓞 K)) (hm : m ≠ ⊥) (χ : RayClassCharacter K m) : ∃ f, DifferentiableOn ℂ f {1}ᶜ ∧ (¬ χ.IsPrincipal → Differentiable ℂ f) ∧ (∀ s, 1 < s.re → HasSum (fun I => χ I * (absNorm I : ℂ) ^ (−s)) (f s)) ∧ ∀ s, 31/32 < s.re → (s ≠ 1 ∨ ¬ χ.IsPrincipal) → f s ≠ 0
  - Library version used by our pages: Main.lean, dedekindZeta_ne_zero, with the same statement.
- **Definition trace.**
  - Because Mathlib's dedekindZeta is an LSeries, the statement asserts an f that is holomorphic on the connected set ℂ∖{1} and agrees with ζ_K on Re s > 1. Such an f is exactly the continuation of ζ_K, so the statement is exactly Thm 1.1 and is not vacuous.
  - The coefficient 31/32 elaborates in ℝ, because s.re : ℝ forces ℝ under binop%.
  - RayClassCharacter is a structure defined in the challenge itself: multiplicative on ideals, zero exactly on ⊥ and on ideals not coprime to 𝔪, and equal to 1 on (α) for every totally positive integral α ≡ 1 mod 𝔪. Its instances are exactly the characters of the narrow ray class group mod 𝔪∞, so the corollary covers every finite-order Hecke character, primitive or not. It is inhabited (the principal character), so the corollary is not vacuous.
  - The library defines RayClassCharacter and IsPrincipal in the same words (Arithmetic/RayClasses/Character.lean l.12–31), and Comparator confirmed the two definitions are identical.
  - The statement covers every K : Type (every number field up to isomorphism) and all heights.
  - Shadowing scan over 57,131 declarations: the only collision is DedekindZetaZeroFree.re, a RingEquiv. It cannot capture s.re, which is field notation on Complex. No global instance among the 1415 is on ℂ, on ℝ or on a bare type variable.
- **Axioms and red flags.**
  - sorry appears only twice, both in the challenge template (the CI build log has exactly these 2 'declaration uses sorry' warnings).
  - There is no axiom, native_decide, ofReduceBool, extern, implemented_by, unsafe, elab or debug option.
  - The 7 opaque declarations are Prop-valued with ':= by' proofs; the kernel type-checks them.
  - The 6 run_cmd calls only generate equation and congruence lemmas.
  - The repository's own CI ran lake build only: it set leanchecker, nanoda and axiom-audit to false.
- **Independent check** (logs/D-palomar-*): the Palomar registry's run 38381813019 ('Verify submission ikz6pd89w5zq', existing id PALOMAR-2026-10-13-000010, 13 Oct 10:13–13:35Z) on commit 09afbb4 reports status pass, stage complete.
  - It ran Comparator on main and heckeCorollary, permitting only the standard three axioms, against a canonical challenge of trust level high whose only dependency is Mathlib.
  - con-ron accepted 206,506 declarations, then nanoda and the Lean kernel, then "Your solution is okay!".
  - The report's sha256 values for Challenge, Solution, comparator.json, lake-manifest and lakefile equal our files at 09afbb4.
  - Tracker PR 114, which would admit this result, has not run any check. Its note "no Comparator or independent kernel replay has been run" refers to the tracker's own port, not to this run.
- **Literature.**
  - Hexagon 2610.00200 (v1 and v2, both 12 Oct) gives the abstract: ζ_K ≠ 0 on Re s > 31/32 for every K, and finite-order Hecke L wherever holomorphic, "proved and formalized in Lean using AxiomProver". This matches the Lean statements exactly.
  - Hexagon labels the paper "a mix of human-written and AI-generated text; Human understanding: some parts".
  - No public doubt.
  - Stronger claims for the same family exist but are human-only, not Lean-checked: Liu, Hexagon 2611.00082, B_2 for every finite-order Hecke L over every number field (audited sound in Hexagon2); before it 2610.00655v2, 0.874073; Hu, 2610.00216v1, 7/8 for unitary Weil L.
- **Barriers.** Different method (metaplectic, Kazhdan–Patterson), and 31/32 lies far above every known floor.
- **Verdict: holds.**

## Common-mode risks (why "holds" is not certainty)
- **Kernels.** (B), (C), (D) and the tracker's check of (A) used the same three binaries (lean bf8d54e4…, nanoda_bin 8241c5e6…, con-ron 4e5616d9…). Only (A) was also replayed with other builds: nanoda 3a24072 and Lean 4.34.1. A shared kernel bug has happened before: Lean PR 14807 (merged 18 Aug 2026) fixed a proof of False that nanoda also accepted. Every toolchain used here postdates that fix.
- **Shared code.** (B) and (C) import OpenAI's baseline, and (D) adapts openai/math analytic code (adc7f124). A logical error is excluded by the kernels, not by independent proofs.
- **Infrastructure.** We trust the tracker's self-hosted runners and the Palomar runner to have run what they logged. Every log is bound by hashes to the exact sources we read.
- **Missing human proof.** No input has a refereed human proof. The texts for (A) and (D) are AI-generated, and (B) and (C) are Lean-only (B_2 itself is explained in Liu's note). For a kernel-checked theorem this matters only for understanding, not for truth.

## Dependency map (part 5)
Inputs:
- **H7/8** = (A).
- **HB** = (B) or (C). Any one of PR 101, 103 or 129 implies every constant our pages use (0.87338196661320752 and 0.8733819666132076); HB implies H7/8.
- **HCT** = (D).
- **Not Lean inputs:** Liu's B_2 for Hecke L over number fields (2611.00082) and Liu's B♭ (2610.26835). Both are human proofs, unrefereed, audited by us.

Scenarios:
- (A) fails alone: nothing falls, since HB implies 7/8.
- (B) and (C) fail: HB-entries fall back to their H7/8 versions, listed in the second column below.
- (D) fails: HCT-entries keep only their abelian-field parts, which follow from HB or A.
- All fail: everything falls back to the classical zero-free regions.

| entry (file) | input | if only H7/8 holds | if HCT fails |
|---|---|---|---|
| Littlewood bound 0.11/max(3, log log q) (littlewood-half-plane) | HB | 0.105 (page) | unaffected |
| Kopp–Lagarias 1.7 (unit-generated-orders) | HB via Littlewood (needs ≥ 0.037) | holds | unaffected |
| Class numbers h ≤ 10,000 (class-number-10000, RECORDS) | HB via Littlewood | needs a rerun with 0.105 | unaffected |
| Numeric Linnik, exponent 7.8978 (linnik-numeric) | HB | exponent 8 (stated, not run) | unaffected |
| Least prime primitive root (prime-primitive-root); unit-group cover 2.1067 log Q + 58 (unit-group-cover) | HB | weaker constants, not computed | unaffected |
| Primes between seventh powers (primes-between-powers) | HB + RH to 3·10¹² + KLN | k = 8 (Visser; R323) | unaffected |
| Prime in (x, x + x^0.856] (primes-short-intervals) | HB + KLN | x^{7/8} (R323) | unaffected |
| Robin to 10^(10^46.79), Nicolas to p_k ≤ 1.78·10⁴⁸ (RECORDS) | HB | 1.76·10⁴⁸; with no half-plane 10^(10^46.44), 6.30·10⁴⁷ | unaffected |
| Mertens power saving (RECORDS) | H7/8 | holds | unaffected |
| Kim–Lee wedge, a = 0.37338… (RECORDS); Jensen d ≤ 1.6141·10²⁵ (N6) | HB | wedge with a = 3/8 (not recomputed); Jensen d ≤ (16T² − 4)/9 | unaffected |
| Λ ≤ 0.107559925318 (N2); Li 106.25 T² (N3); Λ_χ ≤ 0.278828186 (COR) | HB | Raymond 0.107964910079; 105.77 T²; 9/32 | unaffected |
| Zero density below B, Lindelöf at Re s > B, Möbius series in Hardy spaces, differenced ξ, the bound on L′/L(1, χ) by 7.8978 log log q, Miller bases (log n)^{7.8978}, class groups of exponent 8 (COR) | HB | their 7/8 versions (exponent 7.8978 becomes 8; class groups: holds) | unaffected |
| Linnik L ≤ 7/3, NB in L^{8/7}, family density 5/6 < σ ≤ 7/8, Gauss mean values, Dirichlet Li (COR) | H7/8 | holds | unaffected |
| NB L^p for Dirichlet L, p ≤ 1.14497 (DISC N5) | HB | p ≤ 8/7 | unaffected |
| Chebotarev least prime (0.49 log 3Δ_L)^32 (chebotarev-numeric A) | HCT | n/a | falls; part (B), exponent 7.8978, rests on Liu's human Hecke B_2 |
| Random polynomials irreducible without RH (random-polynomials-without-rh) | HCT (BV Prop. 9) | n/a | falls; the Legendre corollary survives with HB (2^{7.9π(d)}) |
| Root-discriminant limits 23.3878 / 65.6439 (discriminant-bounds-strip) | HCT (w = 15/32) | n/a | general K falls; abelian K keeps 26.9662 / 84.2999 via HB |
| Λ_K ≤ 225/512; NB L^p for ζ_K, p ≤ 32/31; koiran; serrenum; detfactor (b); BS-explicit Hecke part; harvest4 items (COR) | HCT | n/a | general K falls; abelian K keeps Λ_K ≤ 0.278828186 and p ≤ 1/B via HB |
| Artin L half-plane; least prime represented by forms, exponent 1.7498; chebotarev (B) (COR, page) | Liu human, not Lean | — | unaffected (not Lean inputs) |

TREE node lines: N2, N3, N4, N5, N6, N7, N9, N10, N11 and N12 cite these inputs as listed above. FORMALIZATIONS l.8–9 (PNTError, QRHCor) import the H7/8 port; l.11 (kimeff) takes HB as a named hypothesis. The DISCOVERIES line "Consequences of the Lean-checked zero-free half-planes (… audit pending)" can now cite this audit; that edit is the orchestrator's.

## Referee verdict
All four inputs are sound theorems at the standard that formal verification provides.
- Each statement says, in Mathlib's correct analytic-continuation definitions, exactly what our pages cite. The one trap, Mathlib's Dedekind zeta being 0 outside convergence, is handled correctly by AxiomMath's existential formulation.
- Each proof passed a statement match against a challenge built from Mathlib alone, used only the three standard axioms, and was accepted by three separately written proof checkers.
- What remains is trust in those checkers, which share binaries across (B), (C), (D), and in the two CI services.
- Nothing here has a human-refereed proof. (A) also has two independent replays on other kernel builds.
- The roughly 15 dependent entries stand as stated. The DISCOVERIES consequences line may drop "audit pending".

Files: logs/{A-*, B-*, BC-*, C-*, D-*, mathlib-defs.txt, scan-summary.txt}; scripts in tools/ (README.md there).
