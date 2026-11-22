BLOCKED-MEMCAP
UTC: 2026-11-22 11:35
sorry count: 0 (no proof files: blocked before the statement compiles)
files: Riemann/Deep/LeanNBRate/: STATUS.md; audit/probe.sh (rerun from the repo root), audit/fp.py (RSS and physical-footprint sampler), audit/probe.log
main theorem: none (planned Riemann.Deep.LeanNBRate.nb_rate: (a) the limit under `(hRH : RiemannHypothesis)`, (b) the value 2 + γ − log 4π when all nontrivial zeros are simple, (c) Nyman's constrained variant)
literature inputs: none yet (candidates per the brief: Burnol math/0103058 TeX l.169–206; Turán–Nazarov; Titchmarsh §9.7, §14.13; Riemann–von Mangoldt; Mellin–Plancherel; Hadamard product)
NOT-FORMALIZED: all of (a)–(c) (blocked)
last reviewer verdict: none
next step: Adrian decides (Options); then nudge this tab.

## Measured (audit/probe.log, 2026-11-22 11:33 UTC; every run under .scratch/memcap.sh 290)
The statement needs `riemannZeta` and `RiemannHypothesis`, defined only in Mathlib.NumberTheory.LSeries.RiemannZeta (l.121, l.185); also MeasureTheory.Lp (LpSpace.Basic), Real.eulerMascheroniConstant (Harmonic.EulerMascheroni), Metric.infDist (HausdorffDistance).

| run | memcap.sh (ps RSS) | fp.py: max RSS / physical footprint |
|---|---|---|
| `lake env true` | killed at 497–660 MB | — |
| `lake env lean`, import RiemannZeta (module or not; PROBE.txt) | killed at 390–661 MB (lake's own RSS) | — |
| `lean --threads=1` direct, Init only | passed between samples | 379 / 75 MB |
| same, `module`, Init only | passed | 230 / 74 MB |
| direct, `module`, import Lean | killed at 364–464 MB | 491 / 136 MB |
| direct, `module`, import Mathlib.Order.Basic | killed at 312–433 MB | 312 / 98 MB |
| direct, import RiemannZeta | killed at 315–385 MB after 1.5 s | 402 / 159 MB |
| direct, `module`, import RiemannZeta | killed at 304–362 MB after 1.6 s | 376 / 158 MB |
| direct, `module`, LpSpace.Basic / EulerMascheroni / HausdorffDistance | killed at 343 / 464 / 572 MB | footprint 146 / 178 / 174 MB |

Every kill came before the import finished: the true peak cannot be measured under the cap. Module files gave no gain before the kill. Physical footprint, 12 Nov (NBLp/LEAN.md, non-module, top MEM or time -l): RiemannZeta import 342 MB (killed), Gamma.Basic 324 MB (killed); NBLp's Basic.lean (RiemannZeta + LSeries.Dirichlet) compiled at 343 MB under the earlier 450 MB rule; simp, field_simp, norm_num, linarith add 130–260 MB.

Cause: memcap.sh caps `ps -o rss`, which counts the clean, file-backed mmapped .olean pages; the brief's own check (`top -l 1 -pid <pid> -stats mem`, physical footprint) is 2.3–5× smaller here. Under memcap.sh no import of Lean or Mathlib and no lake run passes, so T2–T4 are blocked the same way (LeanLi/STATUS.md: same finding).

## Options (Adrian)
| option | T1 | cost |
|---|---|---|
| memcap.sh on physical footprint, cap 290 MB | still blocked: RiemannZeta's import is 342 MB footprint (non-module); module mode unmeasured | none |
| T1's per-process cap raised to ~600 MB footprint or ~1 GB RSS (estimate: import ≈ 343 MB footprint + tactics 130–260 MB) | unblocked; lean direct, no lake | memory |
| cap kept | stays BLOCKED-MEMCAP | — |
