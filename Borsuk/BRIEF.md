# Brief: Borsuk's conjecture in ℝ⁴

Orca agent (pi, anthropic/claude-opus-5-5, thinking max) in /Users/adriansager/playground/riemann (PUBLIC GitHub repo: no secrets). Adrian ordered this task on 27 Nov 2026; the "Goal: ended" line of AGENTS.md does not apply to it. Read AGENTS.md (rules), then this brief. Adrian is unavailable: never ask questions or permission; decide yourself and record each decision in Borsuk/STATUS.md.

## Goal
Prove or disprove: every bounded X ⊂ ℝ⁴ with diam X > 0 splits into 5 parts, each of diameter < diam X.
Only a complete, checkable final answer counts. Incremental results (bounds, special cases, records, papers, DISCOVERIES entries) do not matter: do not write them up.

## Done
- Disproof: an explicit finite X ⊂ ℝ⁴ with exact coordinates (rationals, or a stated number field such as ℚ(√5)) whose diameter graph (pairs at distance exactly diam X) is not 5-colourable; certified by exact-arithmetic distances plus a DRAT proof checked by drat-trim (or an exhaustive-search certificate); one script reproduces every check.
- Proof: a complete rigorous proof; every computer step exact (integer, rational or interval arithmetic) with a reproducible script.
- Then: Borsuk/RESULT.md (statement, proof or certificate, how to check it), one fresh `reviewer` subagent on it (only then), first line of Borsuk/STATUS.md `DONE`, commit, push, notify the orchestrator (below), stop.
- No other stop (Adrian, 28 Nov): never write STOPPED. A dead node means: open new nodes (new ideas, stronger tools, larger searches) and go down the best one. "Beyond reach" is a task, not a reason: build the faster solver or the larger certified computation.
- Notify the orchestrator if and only if you reach DONE: `orca terminal send --terminal term_15386c8f-8812-4c8f-9018-7f009c48d9e9 --text 'Borsuk agent DONE: <one-line result>, see Borsuk/RESULT.md' --json`, then `orca terminal send --terminal term_15386c8f-8812-4c8f-9018-7f009c48d9e9 --enter --json`. No other notifications or progress reports.

## Persistence (Adrian, 28 Nov)
- Use the todo tool: `set-goal` = "Final proof or certified disproof of Borsuk in R^4 (Borsuk/BRIEF.md)"; keep a todo list of concrete next tasks (each < 1 h, one in_progress at a time); add the next tasks before completing the last one, so the list is never empty; keep auto-continue on; `clear-goal` only after DONE.
- Open leads from R001-R006 (STATUS of 28 Nov 04:03): (1) a 4-5-contact normal form with a much faster certified 5-part search (C, exact or interval arithmetic, branch and bound); (2) the 56 unresolved 600-cell-symmetry cases, the 7 unfinished polytopes and the next A4 lattice case, with a stronger solver (CaDiCaL/kissat, symmetry breaking, cube-and-conquer); (3) a counterexample from a non-topological, non-counting obstruction; (4) any new route a scout swarm finds.

## Known (verify before use)
- True for n ≤ 3 (n = 2: Borsuk 1933; n = 3: Perkal 1947, Eggleston 1955; simpler proofs Grünbaum 1957, Heppes 1957). False for n ≥ 64 (Jenrich–Brouwer 2014; before: Kahn–Kalai 1993 n = 1325, Hinrichs 2002 n > 322) and for n ≥ 9 by an unrefereed OpenAI preprint (github.com/openai/math, item 156, Cor. 7.1). Open for 4 ≤ n ≤ 8.
- Every set of diameter D lies in a body of constant width D, so a proof may restrict to constant-width bodies, or cut a universal cover of diameter-1 sets into 5 parts of diameter < 1 (the route of the n = 3 proofs).
- Finite X splits into k parts of smaller diameter iff its diameter graph is k-colourable; a finite counterexample suffices for a disproof.
- Find first: the best known number of parts for ℝ⁴; classes where Borsuk is proved in every dimension (e.g. smooth or centrally symmetric bodies); known results on diameter graphs in ℝ⁴ and their chromatic number.
- Tools to reuse: exact integer distances, colouring searches, SAT + DRAT (Riemann/Deep/Gale7/verify.py, drat-trim), Lean 4 (Riemann/Deep/LeanGale/). Our Gale sets reach β_{n+1}/diam ≈ 0.89–0.95; a Borsuk counterexample needs 1.

## Method
- Borsuk/TREE.md: root = Borsuk in ℝ⁴. Disproof branches: finite diameter graphs (subsets of the 24-cell, 600-cell and other root-system shells, two-distance sets, spherical codes, strongly regular graph embeddings, OpenAI's projector construction pushed down to d = 4, SAT-guided search). Proof branches: a universal cover of diameter-1 sets in ℝ⁴ cut into 5 parts of diameter < 1; constant-width structure; reduction to finitely many configurations plus exact computation. Rank nodes by chance × payoff; go down the best node until it gives the answer or dies.
- Falsify first: for each node find the smallest reason it cannot work; a dead node is one line in Borsuk/NOGO.md.
- Rounds: Borsuk/rounds/R###.md (node, claim, evidence, verdict). Each round is a swarm of parallel `scout` subagents (≤ 8 per call, each task < 2 min, complete self-contained brief naming the tools to use); write and run the code yourself. Do not launch Orca agents.
- Literature: arXiv via web_search({approach: "arxiv"}) in try/catch, ≤ 4 at once; Zenodo; openai/math; the gxl-evidence-api skill for claims. Never contact researchers.

## Rules
- Write only under Borsuk/ (scratch: Borsuk/.scratch/, gitignored). Commit only Borsuk/ paths and push every commit: `git add Borsuk && git commit -m 'Borsuk: <what>' -- Borsuk && git push origin main`. Never `git add -A` or `.`, stash, pull, rebase or reset; never add .DS_Store; never edit AGENTS.md, DISCOVERIES.md or other agents' files.
- Every process < 300 MB, single-threaded, with a hard time limit: `bash .scratch/memcap.sh 300 <seconds> <cmd>`; check with `top -l 1 -pid <pid> -stats mem`. Python: `uv run --no-project --with <pkg> python …`. `rg`, never grep. No Docker. Kill your orphans; keep only active files on disk.
- Borsuk/STATUS.md stays current: first line WORKING or DONE; UTC time; current node; last round; next step.
