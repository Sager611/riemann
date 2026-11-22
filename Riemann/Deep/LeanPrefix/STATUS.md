BLOCKED-MEMCAP
UTC: 2026-11-22 11:42
Sorry count: 0 (no Lean files: blocked before the first compile)
Files: Riemann/Deep/LeanPrefix/: STATUS.md, PLAN.md (statement sketch, lemma DAG, routes), memprobe.txt, tools/ (probe.sh, fpsamp.py, cone.py)
Main theorem: none yet (planned Riemann.Deep.LeanPrefix.prefix_reversals, PLAN.md)
Literature inputs: none yet (planned: OffLineLZZW = arXiv 2610.30892 Thm 3.1, first sentence; Jones 2014 Thm 1.2 (k = 0) only on route J)
NOT-FORMALIZED: the whole Claim (blocked)
Last reviewer verdict: none
Next step: Adrian decides (Options); then nudge this tab with the decision.

## Blocker (memprobe.txt; rerun: bash Riemann/Deep/LeanPrefix/tools/probe.sh)
- The faithful statement needs alternatingGroup (SpecificGroups/Alternating.lean:73), MulAction.IsPreprimitive (GroupAction/Primitive.lean:91), DihedralGroup and SemidirectProduct (wreath product); A_n ≤ G needs Jordan's 3-cycle theorem (GroupAction/Jordan.lean) or Jones's input with Projectivization. Cone of the statement imports: 2,880 modules, 555 MB .olean.
- Under .scratch/memcap.sh 290 no Mathlib import compiles: every import above, and Perm.Sign alone (the 12 Nov floor, 232-246 MB then), was killed at 308-1174 MB RSS within 1.2-3.7 s, before the import finished; physical footprint at the kill 125-162 MB. Only `import Init` passes (155 MB).
- Cause: memcap.sh caps `ps -o rss`, which counts clean, file-backed, mmapped .olean pages; the brief's check `top -l 1 -pid <pid> -stats mem` is the physical footprint, 2.5-4.6x smaller here. With 68% free memory nothing is reclaimed, so RSS ≈ all touched .olean pages; on 12 Nov the same imports passed (RSS ≈ footprint). Same finding: LeanLi, LeanNBRate STATUS.md.
- Footprint estimate (12 Nov calibration): statement imports ~320 MB, Jones route ~395 MB, plus 30-260 MB for tactics: above 290 MB even on the footprint metric.

## Options (Adrian)
| Option | T4 | Cost |
|---|---|---|
| Keep memcap.sh 290 (RSS) | impossible: no Mathlib import, no lake run passes | none |
| Cap the footprint (top MEM) at 290 MB (e.g. LeanLi/tools/memcap_fp.sh) | likely still blocked: statement imports ~320 MB est.; possible only with a lighter restatement (sign.ker for alternatingGroup, IsBlock form of primitivity, wreath product as an explicit subgroup of Perm, Jordan's theorem reproved): deviates from the brief | none |
| Footprint cap ~450 MB for T4's Lean runs | feasible as briefed (route C: explicit 3-cycle + Mathlib's Jordan; no CFSG input) if a uniform 3-cycle construction is found; else route J | ~450 MB per process, RSS ~1 GB (clean pages) |
| Footprint cap ~600 MB | also route J (Jones 2014 input + our Lemma 7 over Projectivization, GaloisField) | ~600 MB per process |
