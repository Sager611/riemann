/-!
# Round 439 — node N12 (function-field and Selberg-class analogues), round 14, audit of the Lean claims behind the Booker–Lee–Zubrilina record: theorem13 (their Thm 1.3, levels N ≤ 36 except 29, 31) and every release endpoint check out with only the standard axioms in the authors' verification records (the scout's 'draft, not certified' came from stale doc comments); one more spanning check (N = 5, Q = 6, rank 6) holds; Cossaboom's GL₂ density bound authenticated (+o(1), unitary cuspidal π); degree-2 conductors 2–4 still open; Bézier/Bernstein bases, ensemble data assimilation and Guttman/Rasch scaling transfer no tool

* Node: N12 function-field and Selberg-class analogues (14th round; R435 before). Round type: audit the Lean and computational claims behind
  a node record (AGENTS.md: Lean axioms and sorries, statement vs literature).
* [audit] Booker–Lee–Zubrilina 2611.07921 v1 (Farmer–Koutsoliotas–Lemurell for holomorphic gamma factors, trivial character, N ≤ 36 except
  29, 31): code at github.com/arbooker/twistless_converse (3,372 .lean files). The scout reported FLAW: lean/QuantitativeArithmeticSmall-
  LevelAllCasesV1.lean l. 10–11, 40–42 and SmallLevelTwentySixOriginalModularV1.lean l. 10–11 call theorem13 a draft, 'not certified until its
  dependency closure passes'. Orchestrator check: ReleaseCheck.lean prints the axioms of theorem13 (l. 9); lean/verification/logs/Quantitative-
  ArithmeticSmallLevelAllCasesV1.log: theorem13, modular, circ, exact_count depend on [propext, Classical.choice, Quot.sound]; level_twenty_six
  the same; V4_ALLCASES_CHAPTERS_CLUSTER_CHECKED_20261106.json lists theorem13 among required_declarations_verified (2,044 modules, 47,018
  declarations, exit 0). So the comments are stale and the claim HOLDS on the authors' records (cluster runs, > 6 GiB; not rerun here;
  fidelity of the Lean statement to the TeX not checked).
* [falsify] BLZ's spanning criterion (TeX l. 213–241): N = 5, Q = 6 has 6 legal rows and exact rank 6 (SymPy and mod 1,000,003): HOLDS;
  levels 29 and 31 are excluded only because the longer proofs of Thm 1.3 were not written (TeX l. 124–127), not by the criterion.
* [check] Cossaboom 2609.27715 Thm mainautomorphic (TeX l. 143–145): N_π(σ, T) ≤ T^{(5/2)(1−σ)+o(1)} uniformly in 1/2 ≤ σ ≤ 1 for unitary
  cuspidal π on GL₂(𝔸_ℚ): node line said '+ε' and 'GL₂/Q', corrected to the exact form; 2408.00617 l. 146–148: A_F(η) = max(d_F, 2) for
  η ≤ 1/10, i.e. density hypothesis for degree ≤ 2 in Re s ≥ 0.9, as quoted: HOLDS.
* [check] degree-2 Selberg class with conductor 2, 3, 4 (2503.02492 l. 165): 2110.00311 needs functional equations for twists, 2009.12329
  covers conductor 1; nothing newer: still open. [check] state of the art since R435: no new arXiv item after 4 Nov; Zenodo 24111210 and
  24020618 only restate function-field spectra: HOLDS.
* Exploration (tool transfers, tested): Bernstein-basis positivity certifies J^{d,n} > 0 only on an interval (d = 2, n = 0: roots ±6.5788i),
  an ensemble Kalman update needs an observation model of the zeros and its PSD covariance carries no RH information, Guttman/Rasch
  scaling needs response matrices: no transfer.
-/

namespace Fayed.R439
end Fayed.R439
