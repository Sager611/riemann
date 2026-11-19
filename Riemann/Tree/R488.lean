/-!
# Round 488 — node N7 (zero verification and zero-free regions), round 20, Lean [check] of fresh formalizations: Mercuri–Mehta's 10⁷ comparator targets are faithful batch statements (the reduction is not comparator-checked); Munawer's Littlewood equivalence is faithful and hypothesis-free, his de la Vallée Poussin region is known (PrimeNumberTheoremAnd, Aug 2026), his forward Nyman–Beurling is our direction, his RH 'criteria' are conditional; Hanke v0.6.0 formalizes Davenport; DKW, Killip–Simon sum rules and Lloyd's theorem transfer no tool

* Node: N7 (20th round; R484 before). Round type: Lean [check] (statement fidelity, axioms, sorry, where certificates live;
  lessons 136, 214), texts pre-extracted in .scratch/r488/txt by orchestrator 44/45 (Zenodo 15–19 Nov).
* [check] Mercuri–Mehta, Zenodo 24247904 (15 Nov; concept 24090489, their 10⁷ run, T = 10,001,072): COMPARATOR.md plus targets for
  all 165,741 batches (orchestrator: tar listing, 1 Euler–Maclaurin + 165,740 Riemann–Siegel; last batch BatchSpec 56 10001047967 1000
  10001072606 1000 ends at 10,001,072.606). BatchSpec (Zeta/Zeros/Certificates/Basic.lean l.105–108, tag v1-1e7) = ∃ a
  CertifiedZeroIntervals with n intervals and the stated end points; CertifiedZeroIntervals (Zeta/Zeros/Interval.lean l.186–190): each
  interval holds t with riemannZeta (1/2 + I t) = 0, intervals chained I.snd ≤ J.fst: a faithful 'n distinct zeros on the line'
  statement. Config: permitted axioms propext, Quot.sound, Classical.choice; nanoda disabled. As COMPARATOR.md says, the reduction
  (seams plus Turing counting ⇒ riemannHypothesisUpTo_10001072) is not comparator-checked and the archived oleans carry no proof
  terms: HOLDS (scout mm1e7 timed out; orchestrator's check). Our 10⁸ check (R456) is unchanged; N7 line updated.
* [check] Munawer, Zenodo 24353571 (Littlewood's equivalence, v3.0.2): riemannHypothesis_iff_mertens_isBigO
  (MathlibPR_LittlewoodConverse.lean l.1536–1542) states _root_.RiemannHypothesis ↔ ∀ ε > 0, M(x) = O(x^{1/2+ε}) with
  M(x) = Σ_{n ≤ ⌊x⌋} μ(n) (Mathlib's moebius); only proved lemmas as Prop arguments; no sorry, axiom or native_decide outside
  comments; VerifyNoTraps.lean prints standard axioms (not rebuilt by us; Mathlib unpinned): HOLDS (scout littlewood), a formalization
  of a classical result.
* [check] Munawer, Zenodo 24353574 (de la Vallée Poussin region): deLaValleePoussin_zero_free_region_unconditional is ∃ c > 0
  (witness min(1/8, c₀)); the parametric version has c = 27/(256 C_pole³ C_lip⁴ C_growth); clean scan (scout dlvp). KNOWN:
  PrimeNumberTheoremAnd's StrongPNT.lean ZeroInequality (∃ E ∈ (0, 1/14), every zero with |γ| ≥ 2 has β ≤ 1 − E/log|γ|; PR #1635,
  4 Aug 2026; Strong PNT proved there 23 Sep 2026; orchestrator). Its CriticalLineZeros.lean sign-change lemma takes an arbitrary Z.
* [check] Munawer, Zenodo 24353573 (forward Nyman–Beurling): MathlibPR_NymanBeurling.lean l.1308–1329 is density ⇒ RH, the direction
  of ours (RH_of_NBApprox, Riemann/Deep/NymanBeurling): KNOWN. Paper L (Zenodo 24375767, v3.0.25): conditional reductions only
  (riemannHypothesis_of_logarithmic_plaque_bound from a plaque bound; WaveletRayleighRHCriterion.lean l.122–124 proves a finite-N
  Rayleigh inequality, not an RH criterion); no RH claim, so nothing to refute (scout nbtitch; orchestrator read §4 of the paper).
* [check] Munawer, RiemannBlueprint v3.0.14 (Zenodo 24375747): 'axiom-free' holds (both sorry hits are in doc comments; no axiom,
  native_decide, implemented_by, opaque or extern); every RH statement is conditional (PurePlaqueRHCriterion.lean l.216, 237:
  TitchmarshMellinProperty plus a plaque Prop); nothing beyond classical results (scout blueprint). Newer versions 24378371 and
  24378373 (19 Nov 19:44) → next round.
* [check] Hanke, analytic-number-theory-lean v0.6.0 (Zenodo 24353815): Davenport Chs. 8, 13, 15, 16, 19 (N(T) with an existential
  O(log T), Page's theorem with explicit 1 − 1/(800 log z), PNT in progressions with the exceptional zero, c₀ = 1/16200); sorry only in
  the Comparator Challenge files, by design; the author says the new statements were not re-read against the books (orchestrator).
* Zenodo sweep from 19 Nov (orchestrator, o44_zsweep.py): 9 records; Prodromov done (R487); 24357107 (Báez-Duarte 'spectral
  criteria'), 24364735 (claimed Goldbach proof), 24369317 (Jacobi-symbol sums), 24356722, 24359794, 24378143 and the two Munawer
  versions pre-extracted in .scratch/r489/txt → next round.
* Exploration: the Dvoretzky–Kiefer–Wolfowitz inequality (zeros are deterministic, no i.i.d. law to sample), Killip–Simon sum rules
  (need a Jacobi matrix whose spectrum is the zeros: the missing Hilbert–Pólya input), Lloyd's theorem and Krawtchouk interlacing
  (no analytic carrier for Li's Jensen polynomials): no transfer.
* Scouts: 8, 7 answered (mm1e7 timed out; its check done by the orchestrator).
-/

namespace Fayed.R488
end Fayed.R488
