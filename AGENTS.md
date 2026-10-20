# riemann — rules

Goal: Fayed arXiv 2209.01890 (v55) "proof" of RH → Lean 4 + Mathlib; pin failing step; then fix it by endless proof-tree search.

## Work loop (never stop)
1. Pick frontier node in `TREE.md`.
2. Branch: ≥8 parallel scouts; directions = some straightforward + many wild/out-of-box/original.
3. Write every branch in Lean (`Riemann/Tree/R###/`): claim, counterexample or reduction; `sorry` only for open obligations, marked.
4. Dead ends → `NOGO.md`, few words each.
5. Hard branch → Orca agent: `pi --provider anthropic --model claude-opus-5-5 --thinking high`, tight full brief; close its tab when done.
6. `lake build` green → git commit per round.

## Discoveries
One-liner + link in `DISCOVERIES.md`; tell user; `reviewer` judges significance; scouts check online novelty; insignificant or known → remove.

## Constraints
- Never ask permission; decide, act, log.
- Orca agents: no new worktrees (.lake ≈ 6 GB, disk tight); own subdir; orchestrator commits.
- Few parallel Lean compiles; scratch in `.scratch/`.
- `rg`, never `grep`. No Docker. No vitest.
- Compile: `lake env lean <file>`; all: `lake build`.

## Layout
`Riemann/Fayed/` paper; `Riemann/Tree/` branches; `TREE.md` index; `NOGO.md`; `DISCOVERIES.md`.
