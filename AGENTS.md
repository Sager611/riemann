# riemann — rules

Goal: Prove or disprove the Riemann Hypothesis by endless Lean proof-tree search (≥8 scouts/round, NOGO log, commit per round, never stop).

## Work loop (never stop)
1. Pick frontier node in `TREE.md`.
2. Branch: ≥8 parallel scouts; directions = some straightforward + many wild/out-of-box/original.
3. Write every branch in Lean (`Riemann/Tree/R###/`): claim, counterexample or reduction; `sorry` only for open obligations, marked.
4. Dead ends → `NOGO.md`, few words each.
5. Hard branch → Orca agent: `pi --provider anthropic --model claude-opus-5-5 --thinking max` (user, 2026-10-20), tight full brief; close its tab when done.
5b. Before any formalization: check Mathlib, openai/math, facebookresearch/atlas-lean, zeta23 trees yourself (GitHub API); scouts' "none known" is unreliable (Hasse was in ATLAS).
5c. Disk low: truncate `.lake/build/ir/**/*.setup.json` to 0 bytes (no rebuild; QRH port had 8.9 GB). Never delete `.c` files (triggers rebuild).
6. `lake build` green → git commit and push per round.

## Discoveries
One-liner + link in `DISCOVERIES.md`; tell user; `reviewer` judges significance; scouts check online novelty; insignificant or known → remove.

## Constraints
- Never ask permission; decide, act, log.
- Orca agents: no new worktrees (.lake ≈ 6 GB, disk tight); own subdir; orchestrator commits.
- Few parallel Lean compiles; scratch in `.scratch/`. Lean LSP off (`.pi-lens.json`): each worker loads Mathlib (~3 GB).
- `rg`, never `grep`. No Docker. No vitest.
- Compile: `lake env lean <file>`; all: `lake build`.

## Layout
`Riemann/Fayed/` paper; `Riemann/Tree/` branches; `TREE.md` index; `NOGO.md`; `DISCOVERIES.md`.
