/-!
# Round 227 — Gale's conjecture fails for every n ≥ 5 (verified, reviewed; novelty check running); hunts clean

* Gale (Orca "gale7", Riemann/Deep/Gale7/): certified counterexamples in R⁷ (84 points, β₈/D = 0.9434435130 > σ₇ =
  0.9419651451) and R⁸ (85 points, 0.9499470677 > σ₈), and an explicit 8-parameter family Vₙ (two regular simplices in
  orthogonal subspaces) whose long-distance graph has no proper (n+1)-colouring by a hand proof, with β_{n+1}/D > σₙ
  checked exactly for every n = 5…1000; with Borsuk counterexamples (n ≥ 64) Gale fails for every n ≥ 5, and with our
  earlier n = 3, 4 sets for every n ≥ 3. Orchestrator re-ran verify.py (R⁷, R⁸), family.py (996 dimensions) and the
  lemma checks: all pass. Reviewer: colouring lemma sound, combination valid, significant. Orca novelty check running.
* No counterexample: 2610.30583 (its explicit counterexample checked), 2611.04862 (Lemma 2.1(2), 162,886 cases),
  2610.40291 (lens minimiser, random starts), 2611.03191 (Gaussian family). Not testable: 2611.04730, 2611.05228,
  2611.03169 (unspecified constant), Zenodo 23964169 (computer-assisted Krzyż n = 7 claim; 54 GB logs omitted).
-/

namespace Fayed.R227
end Fayed.R227
