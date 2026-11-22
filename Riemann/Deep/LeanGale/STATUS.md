BLOCKED-MEMCAP
UTC: 2026-11-22T11:44:41Z
sorry count: 0 (no Lean files: stopped at the import probe, COMMON.md rule)
files: STATUS.md, PLAN.md, tools/fpcap.py (footprint watchdog), tools/probe.sh (import probes), tools/closure.py (olean closure sizes), tools/hom.py (H-homomorphism test)
main theorem: none yet (planned `LeanGale.gale_fails` : ∀ n ≥ 3, ¬ GaleConjecture n)
literature inputs: none yet (planned: Voronov arXiv 2610.28828 for n = 3, 4; a Borsuk-counterexample source for n > 1000, see PLAN.md)
NOT-FORMALIZED: everything (blocked before the first Lean file)
last reviewer verdict: none
next step: Adrian decides (options below); then PLAN.md

## Blocker
The faithful statement needs `EuclideanSpace ℝ (Fin n)`, defined only in `Mathlib.Analysis.InnerProductSpace.PiL2` (l. 114, abbrev of `PiLp 2`). Importing it alone exceeds the 290 MiB cap: killed at 293 MiB (as `module`) / 296 MiB (plain); extrapolated true peak ≈ 310–317 MiB (module), more plain. Even `Mathlib.Analysis.Normed.Lp.PiLp` (the metric of PiLp 2) hits the cap (module: killed at 290.x MiB).

Probe (tools/probe.sh; `lean -j1`, file = one import; physical footprint, lifetime max, MiB):

| import | plain | module |
|---|---|---|
| Mathlib.Order.Basic | 131 | 117 |
| Mathlib.Data.Real.Basic | 211 | 190 |
| Mathlib.Topology.MetricSpace.Bounded (Metric.diam) | 253 | 232 |
| Mathlib.Analysis.Real.Sqrt | 270 | 246 |
| Mathlib.Analysis.InnerProductSpace.Basic | killed (315) | 281 |
| Mathlib.Analysis.Normed.Lp.PiLp | killed (312) | killed (290) |
| Mathlib.Analysis.InnerProductSpace.PiL2 (EuclideanSpace) | killed (296) | killed (293) |
| Mathlib.NumberTheory.LSeries.RiemannZeta | killed (304) | killed (293) |

Run-to-run noise ±10 MiB; `MallocNanoZone=0`, `MallocMaxMagazines=1`, `-DElab.async=false`, `-s 1024` save ≤ 4 MiB.

## Note for the orchestrator: .scratch/memcap.sh measures RSS, not the order's metric
- `top -stats mem` (the order's check) = "Physical memory footprint" (man top). memcap.sh reads `ps -o rss`, which counts Lean's clean, shared, file-backed olean mmaps: `import Mathlib.Order.Basic` has max RSS 842 MB but footprint 125 MiB; `import Mathlib.Analysis.Real.Sqrt` RSS 1.78 GB, footprint 271 MiB. Under memcap.sh every Mathlib import is killed within a second, and `lake env …` itself is killed (lake RSS 623 MB, it elaborates the lakefiles): PROBE.txt's "661 MB" is RSS.
- tools/fpcap.py: same interface as memcap.sh (`uv run --no-project python3 tools/fpcap.py 290 1800 cmd …`), kills the tree when any process's lifetime-max physical footprint (proc_pid_rusage) exceeds the cap. Call `lean` directly with LEAN_PATH (tools/probe.sh builds it) instead of `lake env lean`.
- By footprint RiemannZeta is still over the cap (killed at 293 MiB as module, 304 plain).

## Options for Adrian

| option | allows | statement | cost/risk |
|---|---|---|---|
| A | one statement module importing PiL2 up to ≈ 330 MiB footprint; every other module < 290 | exact: EuclideanSpace ℝ (Fin n), Metric.diam | one process 7–14 % over the order, short (import + one small proof) |
| B | nothing; all < 290 | on `Fin n → ℝ` with √Σ(xᵢ−yᵢ)² and a finite-set diameter (Gale7 style); bridge to EuclideanSpace stated in a docstring, not compiled | not Mathlib's EuclideanSpace/Metric.diam: fails COMMON.md DoD 1 |
| C | nothing; `PiLp 2 (fun _ : Fin n => ℝ)` (= EuclideanSpace by abbrev) | near-exact | blocked too: PiLp's import alone reaches the cap |
