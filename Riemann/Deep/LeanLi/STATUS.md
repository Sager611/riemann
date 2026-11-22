BLOCKED-MEMCAP
UTC: 2026-11-22 11:40
Sorry count: 0 (no proof files; blocked before the first compile)
Files: Riemann/Deep/LeanLi/: STATUS.md, probe/Import.lean, logs/memprobe.log, tools/ (env.sh, probe.sh, fp.c, fpmem.c, cone.py, memcap_fp.sh)
Main theorem: none yet (planned Riemann.Deep.LeanLi.li_remainder_laws: the five parts of discoveries/li-remainder.md)
Literature inputs: none yet (planned: Li 1997 and Bombieri–Lagarias 1999 for Li's formula; Riemann–von Mangoldt, Titchmarsh Thm 9.4; Montgomery 1973, pair correlation under RH; Selberg 1946, ∫S²)
NOT-FORMALIZED: all five parts (blocked); part 5's sentence "a = 2 first fails at n = 5083" stays NOT-FORMALIZED in any case (no certified λ_n computation for n ≤ 5083 fits the cap)
Last reviewer verdict: none
Next step: Adrian decides (Options); then nudge this tab with the decision.

## Blocker (measured 2026-11-22 11:32 UTC; logs/memprobe.log; rerun: bash Riemann/Deep/LeanLi/tools/probe.sh)
- The faithful statement needs Mathlib's `RiemannHypothesis` and `riemannZeta`, defined in `Mathlib.NumberTheory.LSeries.RiemannZeta` (line 185); its import cone: 4,596 modules, 946 MB of exported .olean.
- `lean --threads=1` called directly (no lake), `module` file with `public import Mathlib.NumberTheory.LSeries.RiemannZeta`: memcap.sh killed it at 319 MB RSS before the import finished (physical footprint 138 MB at that moment); non-module file: killed at 340 MB (footprint 133 MB). The true peak cannot be measured under the cap.
- No Mathlib import fits memcap.sh 290: `Mathlib.Order.Basic` killed at 430 MB RSS (footprint 118 MB), `Mathlib.Data.Real.Basic` at 596 MB (footprint 159 MB); `import Init` alone peaks at 220 MB RSS (footprint 67 MB). This holds for every Mathlib-based Lean task.
- lake alone: `lake env true` reaches 620–661 MB RSS (footprint 103–120 MB; it loads core oleans for the packages' lakefile.lean). PROBE.txt's 661 MB is lake's own peak, not Lean's. `lake env lean` and `lake build` cannot run under memcap.sh 290 whatever the imports; tools/env.sh sets LEAN and LEAN_PATH without lake.
- Cause: memcap.sh caps RSS (`ps -o rss`), which counts the clean, file-backed, mmapped .olean pages (reclaimable page cache). Adrian's check `top -l 1 -pid <pid> -stats mem` reports the physical memory footprint (man top), 2.4–3.8× smaller in these runs.
- memcap.sh polls every 0.5 s: runs overshot to 596 MB RSS before the kill.

## Options (Adrian)
| Option | T3 | Cost |
|---|---|---|
| Keep memcap.sh 290 (RSS) | impossible: no Mathlib import compiles | none |
| Cap the physical footprint at 290 MB (= `top -stats mem`): tools/memcap_fp.sh, drop-in for memcap.sh | possible if the RiemannZeta import plus elaboration stays under 290 MB footprint (138 MB at the RSS kill, import unfinished). First step: `. Riemann/Deep/LeanLi/tools/env.sh; bash Riemann/Deep/LeanLi/tools/memcap_fp.sh 290 1800 "$LEAN" --threads=1 Riemann/Deep/LeanLi/probe/Import.lean` | RSS ≳ 1 GB per Lean process, mostly clean reclaimable pages |
| Raise the RSS cap for Lean (≳ 1 GB, unmeasured) | possible | RSS growth under memory pressure |

## Scope if unblocked (preliminary)
- Reusable: LiRate Lean (Theorem 1's core and Corollaries 1–4 from named hypotheses).
- Ours, to prove: Lemma A (explicit zero sum, constant 7/4); Theorem 2 (Lemma A, Schur's test; Selberg's ∫S² as input); Theorem 4 (smoothed mean 7/4, φ ∈ C_c^∞(0, ∞)); Lemma 0 (¬RH ⇒ |R(n)| ≥ e^{cn} infinitely often, via log ξ(1/(1 − z))); the Fourier estimates of the ψ∗ψ weights; the trend transform; the τ_T asymptotics; the instantiation of LiRate's abstract zero family by ζ's zeros with multiplicity.
- To define faithfully: λ_n as the symmetric limit Σ_{|γ|≤T} over nontrivial zeros with analytic multiplicity; Merlini–Sala–Sala's λ_tiny exactly as in 1903.06233.
- Mathlib has no Hadamard product, zero counting, Stirling in vertical strips or explicit formula: large new analytic Lean even with literature inputs (rough estimate ≥ 5,000 lines).
