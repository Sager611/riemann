# riemann — rules

> [!IMPORTANT]
> Editing rules: ≤ 100 lines excluding the Changelog. Maximum deletion, minimum addition: each change leaves the fewest words with the densest information. Rules only; reasons, dates and anecdotes go to the Changelog. Every change appends one Changelog row.

Goal: prove or disprove RH by endless Lean proof-tree search; never stop.

## Priority

- New discoveries (theorems, records, counterexamples, false steps) advancing a named node of the RH approach tree at the top of `TREE.md`: RH criteria/equivalents, quantitative zero or prime-error records, L-function analogues, barriers to approaches.
- Exploration share: 3 of 8 scouts + 1 Orca agent on remote fields, no RH link needed upfront; judged by whether it yields a new idea for an RH node.
- Formalizing known results is secondary.
- Formalizing a new discovery is mandatory.

## Work loop

1. Pick one tree node (seed the tree if missing); rank by RH relevance × chance × gain.
2. ≥ 8 scouts on it: state of the art, computation (extend ranges/constants), adversarial (counterexamples to lemmas the approach needs), wild.
3. Record the round in `Riemann/Tree/R###.lean` (claim, counterexample or reduction; `sorry` only for marked open obligations); dead ends → one line in `NOGO.md`.
4. Work beyond scouts (edit, build, research) → Orca agent with a tight full brief; close its tab when done.
5. Check green → commit → `git push origin main`, every commit. Public repo: no secrets.

## Discoveries

- `DISCOVERIES.md`: one line (≤ 25 words) per entry → `discoveries/<name>.md` (statement, before, method, review, novelty, RH node advanced, files); off-tree results under "Side results"; tell the user.
- Entry needs own verification, a `reviewer` verdict and an adversarial novelty check (arXiv monthly listings `https://arxiv.org/list/<cat>/<YYYY-MM>?show=2000`, this and last month; Zenodo; our NOVELTY.md files). Known or insignificant → remove.
- Broken claimed proofs → `REFUTATIONS.md`; routine/conditional deductions → `COROLLARIES.md`; formalizations → `FORMALIZATIONS.md`.

## Layout

`Riemann/Fayed/` paper; `Riemann/Tree/` rounds; `TREE.md` tree + round index; `NOGO.md`; `DISCOVERIES.md`.

## Changelog


| Timestamp (UTC)  | Who          | What                                                                                                                                     | Why                                                                        |
| ---------------- | ------------ | ---------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| 2026/10/20 01:28 | orchestrator | Created: rules for checking Fayed's RH proof in Lean                                                                                     | Project start                                                              |
| 2026/10/20 16:10 | orchestrator | Orca agents: `--thinking max`                                                                                                            | Adrian                                                                     |
| 2026/10/20 17:15 | orchestrator | Check external Lean libraries before formalizing                                                                                         | Scouts' "none known" missed Hasse in ATLAS                                 |
| 2026/10/20 18:27 | orchestrator | Goal: prove or disprove RH; push each round                                                                                              | Adrian                                                                     |
| 2026/10/20 21:34 | orchestrator | Truncate `*.setup.json` when disk is low                                                                                                 | QRH port used 8.9 GB                                                       |
| 2026/10/20 23:47 | orchestrator | Priority: new discoveries over formalizations                                                                                            | Adrian                                                                     |
| 2026/10/21 01:28 | orchestrator | DISCOVERIES: new mathematics only; formalizations → FORMALIZATIONS.md                                                                    | Adrian                                                                     |
| 2026/10/21 01:36 | orchestrator | Never ultraworker; use Orca agents                                                                                                       | Adrian                                                                     |
| 2026/10/21 01:41 | orchestrator | Refutations and corollaries out of DISCOVERIES                                                                                           | Adrian                                                                     |
| 2026/10/21 05:48 | orchestrator | Kill orphaned scout processes                                                                                                            | A 15 h mpmath loop was found                                               |
| 2026/10/21 09:32 | orchestrator | Nudge agents stalled on truncated responses                                                                                              | Silent stalls                                                              |
| 2026/10/21 09:40 | orchestrator | Idle agents may need Enter twice                                                                                                         | Nudges ignored                                                             |
| 2026/10/21 17:50 | orchestrator | Push after every commit                                                                                                                  | Adrian                                                                     |
| 2026/10/21 18:27 | orchestrator | Kill hung tool calls ("Elapsed Nh")                                                                                                      | Hung calls                                                                 |
| 2026/10/21 22:05 | orchestrator | Disk emergency playbook                                                                                                                  | Disk nearly full                                                           |
| 2026/10/21 23:16 | orchestrator | Keep only active work on disk                                                                                                            | Adrian                                                                     |
| 2026/10/21 23:45 | orchestrator | Write `.lean` via bash                                                                                                                   | pi-lens spawned `lake serve`, which cloned Mathlib                         |
| 2026/10/22 00:47 | orchestrator | Scouts: arXiv-only web search; TeX for notation flaws                                                                                    | Exa out of credits, keyless search throttled                               |
| 2026/10/22 01:03 | orchestrator | `.lake` placeholder file                                                                                                                 | Stray Mathlib clones                                                       |
| 2026/10/22 01:30 | orchestrator | Novelty: scan arXiv monthly listings                                                                                                     | A missed 8 Oct paper (Λ ≤ 0.13963) voided our de Bruijn–Newman entry       |
| 2026/10/23 01:30 | orchestrator | DISCOVERIES: one line per entry, details in `discoveries/`                                                                               | Adrian                                                                     |
| 2026/10/23 03:35 | orchestrator | Orca harness: pig → pi                                                                                                                   | Adrian uninstalled pig                                                     |
| 2026/10/23 05:47 | orchestrator | Orca: pi-lens off, low-power cores, single-threaded                                                                                      | Adrian's CPU pause; pi-lens cost \~600 MB per agent                        |
| 2026/11/05 05:01 | orchestrator | Scope discoveries to the RH approach tree                                                                                                | Adrian: discoveries were not moving toward RH                              |
| 2026/11/05 05:06 | orchestrator | Exploration share: 3 of 8 scouts + 1 Orca agent on remote fields                                                                         | Adrian: RH may need ideas from remote fields                               |
| 2026/11/05 05:16 | orchestrator | Terse rewrite; callout; Changelog; check own NOVELTY files                                                                               | Adrian: file too fluffy; the sofa "record" had missed Tomar's better bound |
| 2026/11/05 05:19 | orchestrator | Moved the Orca agents, Scouts, Lean and disk sections to the workspace daydreaming topic `2026-11-05 - riemann-orca-scouts-lean-disk.md` | Adrian: operating notes belong in daydreaming                              |
| 2026/11/11 12:50 | Adrian       | Priority: formalizing a new discovery is mandatory; markdown reformatted                                                                 | Adrian                                                                     |
