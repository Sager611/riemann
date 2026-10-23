# riemann — rules

Goal: Prove or disprove the Riemann Hypothesis by endless Lean proof-tree search (≥8 scouts/round, NOGO log, commit per round, never stop).
Priority (user, 2026-10-20): hunt NEW discoveries (new theorems, records, counterexamples, false steps); formalizing known results is secondary. DISCOVERIES.md strict: new theorems, records, counterexamples to established statements only. Broken claimed proofs → REFUTATIONS.md; routine/conditional deductions → COROLLARIES.md; formalizations → FORMALIZATIONS.md.
Disk (user, 2026-10-21): keep on disk only what is actively worked on; delete scratch, downloads, caches and `.lake` of finished work immediately (re-fetch Mathlib with `lake exe cache get` when a build is needed).

## Work loop (never stop)
1. Pick frontier node in `TREE.md`.
2. Branch: ≥8 parallel scouts; directions = some straightforward + many wild/out-of-box/original.
3. Write every branch in Lean (`Riemann/Tree/R###/`): claim, counterexample or reduction; `sorry` only for open obligations, marked.
4. Dead ends → `NOGO.md`, few words each.
5. Never ultraworker (user, 2026-10-21): all edit/build/research work beyond scouts → Orca agent: `taskpolicy -b nice -n 19 pi --no-lens --provider anthropic --model claude-opus-5-5 --thinking max` with OMP/OPENBLAS/MKL/VECLIB_MAXIMUM/NUMEXPR thread counts = 1 (user, 2026-10-20; harness is pi again: user uninstalled pig on 2026-10-22, never reinstall it; CPU limits after the user's CPU pause, 2026-10-23: low-power cores, one single-threaded job per agent; pi-lens off: its pyright/opengrep servers cost up to ~600 MB per agent), tight full brief; close its tab when done.
5b. Before any formalization: check Mathlib, openai/math, facebookresearch/atlas-lean, zeta23 trees yourself (GitHub API); scouts' "none known" is unreliable (Hasse was in ATLAS).
5c. Disk low: truncate `.lake/build/ir/**/*.setup.json` to 0 bytes (no rebuild; QRH port had 8.9 GB). Never delete `.c` files (triggers rebuild). Disk emergency (<1 GB): stop our heavy jobs; check `/var/folders/*/*/X/*code_sign_clone` (stale Chrome clones), other projects' `/tmp/wt-*`, swap; our QRH oleans (1.6 GB) are deletable (rebuild costs hours).
5e2. Search down (Exa out of credits, keyless rate-limited, 2026-10-21): scouts may use only `web_search({approach:"arxiv"})` inside try/catch. pdftotext drops overlines/accents: confirm notation-level flaws in the TeX (`https://arxiv.org/e-print/<id>`). Under heavy load run scouts at concurrency 4. Doc-only Tree files: check with plain `lean <file>` (no lake, no Mathlib). Write `.lean` files via bash/python, never the write/edit tools: pi-lens then spawns `lake serve`, which clones Mathlib into `.lake` despite `.pi-lens.json`. `.lake` is a placeholder FILE that blocks stray lake/LSP Mathlib clones; delete it only for a real build.
5d. Scout timeouts can orphan their python/bash (ppid 1, cwd = repo); check `ps -Ao pid,ppid,etime,args | awk '$2==1'` and kill ours (a 15 h mpmath loop was found).
5e. Orca agents stall silently on "Response was truncated": if no Working line, read the terminal tail and nudge ("continue; write long derivations to files"); an idle agent may need --enter twice. Hung tool calls show "Elapsed Nh" in the terminal: kill the process; briefs must demand timeouts on scripts.
6. `lake build` green → git commit, then `git push origin main` after EVERY commit (user, 2026-10-21; repo is public: never commit secrets).

## Discoveries
`DISCOVERIES.md`: one short line per entry (≤ 25 words) linking to `discoveries/<name>.md` (statement, before, method, review, novelty, files; user 2026-10-23); tell user; `reviewer` judges significance; scouts check online novelty; insignificant or known → remove. Records: also scan arXiv monthly listings (`https://arxiv.org/list/<cat>/<YYYY-MM>?show=2000`, math.NT/CV/GM, this and last month) and Zenodo; a missed 8 Oct paper (Λ ≤ 0.13963) voided our de Bruijn–Newman entry.

## Constraints
- Never ask permission; decide, act, log.
- Orca agents: no new worktrees (.lake ≈ 6 GB, disk tight); own subdir; orchestrator commits.
- Few parallel Lean compiles; scratch in `.scratch/`. Lean LSP off (`.pi-lens.json`): each worker loads Mathlib (~3 GB).
- `rg`, never `grep`. No Docker. No vitest.
- Compile: `lake env lean <file>`; all: `lake build`.

## Layout
`Riemann/Fayed/` paper; `Riemann/Tree/` branches; `TREE.md` index; `NOGO.md`; `DISCOVERIES.md`.
