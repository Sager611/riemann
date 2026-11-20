READY TO APPLY: all 14 anchors are unique in today's TREE.md, 10 node sentences now match the nodes' current records and pending items, the five pending checks are settled (all become ready COROLLARIES lines except class-group exponents 6 and 7, which stay conditional), and 'python3 Riemann/Deep/PostRH/apply.py --write --date "<d Mon YYYY>"' passes every assertion in its dry run; the five new COROLLARIES lines still want a reviewer.

# postrh — consequences of the claimed Lean RH + Dirichlet GRH (ROOT, TREE.md l.18)
Orca agent postrh, 20 Nov 2026 18:37–(see end) UTC. Drafts only; no repo file outside Riemann/Deep/PostRH/ edited; nothing committed.

## Answer first
- What the Lean statements give (CONSEQUENCES.md §1): RH for ζ and GRH for every Dirichlet L (all moduli, imprimitive and principal included; .scratch/alp/lean/src/Statements.lean l.3–16), hence GRH for ζ_K of abelian K, GLH for Dirichlet L, no Siegel zeros. Nothing else.
- Becomes unconditional (CONDITIONAL.md §A, §B; DRAFT.md §2): Λ = 0; λ_n > 0; Robin/Lagarias/Nicolas; d_N → 0 and liminf d_N² log N ≥ C; every J^{d,n}; ≥ 19/27 simple zeros; DPC's gaps (< 1/2 + ε, λ > 3.63846723); Lindelöf, Harper's moments, Simonič's explicit bound; Schoenfeld's ψ, θ, π − li bounds; CMS primes in [x, x + (22/25)√x log x]; Soundararajan's Mertens bound; LLS's explicit Linnik (φ(q) log q)², non-residue (log q)², Littlewood and class-number bounds; Grosswald; Bach's Miller witnesses; 26 GRH- or RH-conditional results in our own files (CONDITIONAL.md §A; elite primes, Alladi–Johnson, LPS PageRank, Legendre-sum moments, quadratic character sums, Elkies primes, Gauss mean values, super-positivity for every real χ, our Li-remainder, odd-Robin, Nyman–Beurling and small-gap theorems, …).
- Moot (CONDITIONAL.md §C, 28 rows): our Λ, Li, Robin, Jensen, Linnik, Littlewood (q ≥ 10^10), x^0.856, k-th-power, small- and large-gap, Dirichlet half-plane and density records; the 67.36% simple-zero record (by 19/27); |ζ(1/2+it)| records only for t ≥ 10^50.81 (Fodge's; ours from 10^46.92 and 10^48.09).
- Stays conditional (§D, §E): non-abelian Dedekind ζ, effective Chebotarev, Artin's conjecture (Hooley's HRH), ring/Hilbert-class Hecke L (incl. Hecke L of Q(√−3) that are not Galois-invariant), weight-one/Artin L, elliptic-curve and GL(2) L, everything needing simple zeros (BCF's d_N² ~ C/log N, the constant C in our continuous rate), well-spacing (Bombieri–Hejhal 100%), Conjecture H (our 7/8).
- Re-aim (OPEN.md): 1 GRH beyond Dirichlet (smallest case ζ_K for the cubic field of discriminant −23, i.e. L(s, ψ) for a cubic class-group character of ℚ(√−23)); 2 simple zeros beyond 19/27; 3 pair correlation; 4 GL(2)/Selberg class; 5 Hilbert–Pólya; 6 Mertens; 7 Artin; 8 explicit records; 9 Legendre.

## Method and checks
- Inputs read-only: TREE.md l.5–18 (sed), COROLLARIES/RECORDS/DISCOVERIES, discoveries/, Exp48 and DPCAudit REPORTs (rg); the brief's rg pattern found 56 + 7 + 125 + 11 matching lines; 13 scouts classified them (.scratch/postrh/s_chunk00–12.md); I re-checked every UNC and corrected 7 scout verdicts (CONDITIONAL.md §F).
- Literature: 26 arXiv e-prints fetched to .scratch/postrh/tex/ (1801.05914, 1302.5018, 1810.08843, 2610.28527, 0908.2008, 2010.13307, 1305.4618, 1410.7015, 1708.04122, 0705.0723, math/0008177, math/0202141, 1902.07321, math/0404394, 2005.05142, 1309.3595, 1508.05182, 1803.02056, math/0412262, math/0305114, 1211.5191, 2006.13771, 2412.17199, 2304.13605, 2509.24119, 2504.15658), 2 more by API metadata only (0810.3587, math/0612106); every quoted hypothesis line verified mechanically: code/quotecheck.py, 72/72 quotes found at the stated line ± 3 (logs/quotecheck.log).
- Computation: code/n9cross.py (mpmath): Simonič's RH bound crosses our 0.568 and 0.4821 t^{1/6} log t at log₁₀ t = 46.92, 48.09 and Fodge's bound at 50.81, staying below afterwards (sampled to 10^100; logs/n9cross.log). memcap peaks: TeX fetches 21 MB and 13 MB; n9cross.py and quotecheck.py ended before the first 2 s sample (reported 0 MB).
- Falsification first: 10 traps checked explicitly (CONDITIONAL.md §E), e.g. Hooley's HRH (Moree math/0412262 l.555), Heath-Brown's elliptic GRH (math/0305114 l.64), BCF's hidden simplicity (1211.5191 l.249), Hecke L of the abelian field ℚ(√−3), Legendre.


## Caveats (falsification of this report)
- Coverage: our files fully (the brief's pattern); the literature only on the node lines and only what I could quote from arXiv TeX in the time limit. Non-arXiv sources are named, not quoted: Robin 1984, Nicolas 1983, Schoenfeld 1976 (via Büthe), Bach 1990, Li 1997, Bombieri–Lagarias 1999, Montgomery 1973, GGOS 2000.
- UNC? rows wait on one named check each (DRAFT.md §3): Dobner's normalization vs COROLLARIES.md:18; EKN/ARSX class-group computations; GGOS's GRH = Dirichlet-only; L-functions of Mangerel, Gaillard–Jindal, Martin, Brun 2504.15658.
- "Moot" means superseded as a frontier result; explicit constants that RH does not make explicit (ε_c in N1, the all-x Mertens inequality, the uniform polylog ĝ(p) bound, |L(1, χ)| for q < 10^10, |ζ(1/2+it)| below 10^50.81) stay records.
- Everything is conditional on ROOT's verification, which this agent did not touch (alpaudit, alpmath, alpcore).


## Reviewer
- Independent reviewer (fresh context, 46 rows A1–A24, B1–B22 pasted inline with their quotes): INCOMPLETE, 42 CORRECT, 4 GAP. All four repaired in CONDITIONAL.md: A10 (super-positivity only for non-principal real χ; χ = 1 excluded), A12 (primitive χ mod q ≥ 3; conductor one is A14), A15 (mean square between (c₁ − o(1))N^{3/2}√log N, c₁ = 0.0557905…, and O(N^{3/2}√log N·log log N), not an exact order), B15 (Li coefficients of the completed primitive ξ(s, χ), not of raw imprimitive L). Also fixed: cross-references A21, A22 → B18, B20; scope kept as in the sources (A3 ≍ not ~, A4 even m, A7 the paper's curve assumptions, A9 only the stated exponents, B9 Büthe's thresholds use finite verification, B20 −q fundamental).
- The reviewer also judged A13 (EKN), B3 (GGOS's GRH is Dirichlet) and B16 (Dobner) CORRECT; I keep them UNC? because those residual checks were asserted, not shown (EKN l.488–499 takes small discriminants from tables and an undescribed Magma run).
- Delta review (second fresh reviewer: repaired A10, A12, A15, B15 and the new A25–A27, added after the first review from their TeX): 6 CORRECT, 1 GAP: A15's rounded 0.0558 exceeded the proved c₁ = (3/5)^{1/2}(4/5)²/(2π√2) = 0.05579056… (LiAudit/FULLPROOF.md F4); repaired to the exact c₁ in CONDITIONAL.md and CONSEQUENCES.md.

## Masullo addition (orchestrator 64, 19:20 UTC)
- A. M. Masullo's claim (18:18 UTC, unverified) of GRH for the Hecke L-functions of every number field gives GRH for every ζ_K and zero-free Re s > 1/2 for every Artin L (Brauer induction), not Artin holomorphy, non-CM GL(n ≥ 2) or non-CM elliptic curves. CONDITIONAL.md §D and §E now carry a Masullo column (Alpöge verdicts stay the main column): D 16 discharged (5 with a [check]), D25 partly, D26 only for Hecke-type families, 8 not; E 6 discharged, 4 not. DRAFT.md §4 and OPEN.md (line after the table) say which items move: OPEN item 1 shrinks to non-CM GL(2) and beyond (smallest case L(E, s), E = 11a1, or L(Δ, s)), item 7 (Artin) would become a theorem.
- Not reviewed by a subagent (added in the 15-minute window).

## Files
- Riemann/Deep/PostRH/REPORT.md: this report.
- Riemann/Deep/PostRH/CONSEQUENCES.md: the five Lean statements, derived D1–D4, the not-given list; one consequence line per node N1–N13 and ROOT, with sources and the records made moot.
- Riemann/Deep/PostRH/CONDITIONAL.md: §A our UNCONDITIONAL rows (27), §B literature UNCONDITIONAL rows (22, TeX quotes), §C MOOT (28), §D our STAYS CONDITIONAL (26), §E traps (10), §F scout verdicts corrected; §D and §E carry the Masullo column.
- Riemann/Deep/PostRH/OPEN.md: 9 open problems ranked by relevance × chance × gain, each with a falsification target and a first round; one line on Masullo's claim.
- Riemann/Deep/PostRH/DRAFT.md: §1 14 TREE.md insertions (unique anchors, sentences ≤ 25 words), §2 COROLLARIES block (38 lines, ≤ 30 words each), §3 pending checks, §4 the items Masullo's claim would move.
- Riemann/Deep/PostRH/code/n9cross.py, logs/n9cross.log: heights where Simonič's explicit RH bound beats the explicit t^{1/6} bounds.
- Riemann/Deep/PostRH/code/quotecheck.py, logs/quotecheck.log: mechanical check of 72 quoted hypothesis lines.
- Scratch, not deliverables: .scratch/postrh/STATE.md, SCOUT_BRIEF.md, chunk00–12.txt and m_*.txt (scout inputs), s_chunk00–12.md (scout tables), review_rows.md, tex/ (26 e-prints, lit_hits*.md), lit*.py.

## Recheck 21:15–21:40 UTC (postrh2, 20 Nov 2026)
- Previous verdict line (postrh, 19:25 UTC), superseded by line 1: READY (14 node consequences, 45 UNCONDITIONAL, 31 STAYS CONDITIONAL, 9 open problems): if Alpöge's five Lean statements are verified, every conditional result that needs only RH for ζ or GRH for Dirichlet L becomes a theorem (5 more after one named check each), while everything needing non-abelian, Hecke, Artin or GL(2) GRH, or simple zeros, stays conditional.
- Anchors (§1): all 14 occur exactly once in TREE.md (python str.count = 1 and rg -F -c -e = 1; TREE.md at HEAD 1259790d, 582 lines; node lines 5–18). None of the 14 sentences is present yet.
- What changed on the node lines since the 19:03 UTC drafts (word diff against 2b6f9b7d, the last TREE.md commit before 19:03 UTC): N2 +explore66 (Ważewski blocks, barrier), N8 +explore66 (Conley index barrier), N13 +explore67 (eigenvalue-crossing barrier), ROOT +Sutherland/alpdep/postrh/R538 text and its tag '[check]' → '[check: every decisive outside check passed except the Lean kernels' verdict …]'. N1, N3–N7, N9–N12 byte-identical.
- The ROOT line's tag is the bracket right after the anchor ('— [check: …] Claimed proof …'): an insertion right after the anchor would push the tag behind our sentence, so the ROOT row now inserts after the closing ']' of the leading tag (apply.py finds the tag at run time; the tag text is never edited).

- Word count rule used: tokens containing a letter or digit (≥, >, <, +, ~, → not counted); every §1 sentence now has ≤ 25 (N1 has 24 such words, 27 whitespace tokens).
- §1 sentences changed (old → new, why), all against the node's current text (sed -n '<k>p' TREE.md | cut, keyword windows):
  - N2: dropped '[check]' and added Λ_K = 0 for abelian K and 'Λ_χ records moot' (Dobner check settled, §3; COROLLARIES R324's Λ_χ ≤ 0.2788 becomes moot).
  - N3: + 'PR 144's included' (node l.7 lists 'PR 144's B, pending: 110.4775 T²', a pending Li-range record).
  - N4: 'range records moot' → 'range, 2-adic, t-free records moot' (node l.8: robin2adic v_2(n) ≤ 72/75 record, 21-free, 25-free 'candidate pending review'); 'n ≥ 5041' → 'n > 5040'.
  - N7: + 'QRT PR 144' and 'verification heights' (node l.11 opens 'Pending: QuasiRiemannTracker PR 144 … Re s > 0.873332794', ends with the formal/non-formal zero-verification records).
  - N8: names the records ROOT beats (node l.12: 67.3639% simple and on the line, distinct on the line 0.673639; ours: small gaps 7/11 (dpcbeyond), large gaps 3.367 (explore53)) and adds ≥ 84.86% distinct (CGdL under Dirichlet GRH, §3); attribution 'Bui–Heath-Brown' left to COROLLARIES for length.
  - N9: + '13/84 moot' (node l.13 record μ(1/2) ≤ 13/84); 10^50.81 rechecked against logs/n9cross.log (Fodge crossover log10 t = 50.8061; Fodge is still the node's t^{1/6} record above 2.9·10^12).
  - N10: + '7/8' (node l.14 record ψ(x) = x + O(x^{7/8+ε}), M(x) = O(x^{7/8+ε})); x^0.856 (ingham4) and Mertens exponent 0.9382 (mertensch) still the node's records; the explicit all-x Mertens inequality stays (CONDITIONAL.md:70).
  - N11: + 'Kerr's Siegel-zero exclusion superseded' (node l.15's record: Landau–Siegel zeros excluded, Kerr 2610.21099; GRH implies it).
  - N12: 'our GL(2) records stand' → 'our linear-combination and root-discriminant records stand' (node l.16 has no GL(2) record of ours; its records of ours are lincomb3, selbergn4, odlyzko31).
  - ROOT: bare paths → markdown links; placed after the tag (above).
  - Unchanged, still matching: N1 (window records), N5 (C − 6·10⁻¹² and p ≤ 1.14497 records covered by 'liminf ≥ C' and 'L^p for p < 2'), N6 (d ≤ 1.6144·10^25 and the Dirichlet degree records), N13 (explore67's barrier concerns CCM step one; simple ground state still open).

- §3, the five pending checks, settled from the e-print TeX (.scratch/postrh2/tex/<id>/; file and line numbers in DRAFT.md §3):
  - B16 Λ_χ, Λ_K (Dobner 2005.05142 paper.tex): UNCONDITIONAL-UNDER-ROOT → §2. Same heat flow as COROLLARIES.md:18 (ξ^F_t((1+iz)/2) = ∫e^{tu²}Φ_F(u)e^{izu}du, l.236–240; γ unique up to scaling, l.244; Λ_F ≥ 0, l.119); Λ = 0 is scale-invariant anyway. The line also makes our abelian Λ_K and Λ_χ bounds (COROLLARIES.md:18, 28, 57, 58) moot; the all-field Λ_K ≤ 225/512 (line 28) stays, as it covers non-abelian K.
  - A13 exponents 3, 5 (EKN 1803.02056; ARSX 2605.26020): UNCONDITIONAL-UNDER-ROOT → §2. EKN's Thm 1 (l.148–175) is stated without hypothesis for |D| ≤ 3.1·10^20, ERH only bounds |D| (Bach–Sorenson split prime for the quadratic field, l.250–258: 9.7·10^10, 2.3·10^20, l.268–275); ARSX l.262–267, 666. Exponents 6, 7: STAYS CONDITIONAL, on a search nobody has run (ERH bounds 2.5·10^25, 3.9·10^30 > 3.1·10^20; ARSX l.268–269).
  - B3 CGdL 1810.08843: UNCONDITIONAL-UNDER-ROOT (l.257 defines GRH as 'for Dirichlet L-functions'; GGOS's input l.494–497 under it). 0.6845 itself is below 19/27, so no line; but the same paper's GRH results N_d(T) ≥ (0.8486 + o(1))N(T) (l.302–310) and the 0.9350 Dirichlet average (l.359–362), absent from postrh's rows, are now one §2 line and N8's '≥ 84.86% distinct'. Two arXiv searches (approach arxiv) found no later improvement (not exhaustive; 2501.14545 is a box-condition 2/3).
  - Martin thm2 (2509.24119): UNCONDITIONAL-UNDER-ROOT → §2: ERH 'for quadratic Dirichlet characters' (l.174) enters thm2 (l.188) only through the exponent-2 and -3 lists (l.620, 665). Found on the way: A27's §2 line (finiteness, Thm thm:fin) was already unconditional (COROLLARIES harvest4 (c), chebexp split prime); replaced by the thm2 line.
  - D26 Brun B < 2.1594 (Dunn 2504.15658): UNCONDITIONAL-UNDER-ROOT → §2. GRH enters only at l.187–197 via Johnston–Starichkova 2208.01229 Lemma 5.1 (l.621–638), proved from Grenié–Molteni's Chebotarev for ℚ(ζ_q)/ℚ (1709.07609 assumes RH for Dedekind zeta functions; ℚ(ζ_q) and its subfields are abelian, so these are products of Dirichlet L). Moduli [d₁, d₂] ≤ 2 lie outside Lemma 5.1's 3 ≤ q; Schoenfeld's RH bound (≤ √x log x/8π) covers them.
- DRAFT.md §2 is now 42 lines: 38 + B16, B3 (0.8486, 0.9350), D26 (Brun), A13 (exponents 3, 5), with A27 replaced by Martin's thm2. These five lines are mine and have no reviewer verdict yet. Orchestrator 67 (R546, 22:00 UTC): all five checked against their TeX (Dobner l.110–119, 194–217; CGdL l.302–310, 359–363; Dunn l.43, 187–247 with Johnston–Starichkova l.621–638 and Grenié–Molteni l.180; EKN Thms 1–2 and 2605.26020 Thm 1; Martin l.174, 188, 620, 665): all HOLD.
- apply.py (Riemann/Deep/PostRH/apply.py, python3, stdlib only): parses DRAFT.md §1 (table rows) and §2 ('- If ROOT holds' lines), so later DRAFT edits flow through; exact-substring edits only, each edit point asserted unique; ROOT's sentence goes after its leading [ … ] tag (bracket-depth parse), and it asserts the tag stays first and unchanged; asserts: every anchor count = 1, no sentence present, no 'If ROOT holds (Alpöge, verified' heading and no block line in COROLLARIES.md, TREE.md line count unchanged, files unchanged since read (sha) — all before the first write; --write writes TREE.md then COROLLARIES.md via temp file + os.replace; --date sets the heading's date (default: run date, UTC); --preview DIR writes the would-be files to a scratch dir only.
- QA: dry run ALL ASSERTIONS PASSED (log below, logs/apply_dryrun.log). Preview (.scratch/postrh2/preview) diffed by .scratch/postrh2/previewcheck.py: exactly lines 5–18 change, each by one contiguous insertion right after its anchor (ROOT: after '…page' below)] '), the ROOT line still starts with its tag, COROLLARIES.md gains the heading and 42 lines as a pure append. One guard test passed --write together with --preview .; apply.py refused it before any write branch ('FAIL … Nothing written'), and git status shows TREE.md and COROLLARIES.md unchanged. apply.py was never run with --write otherwise. Negative test (.scratch/postrh2/negtest.py, dry run on copies in a scratch mirror): a node line already carrying 'If ROOT holds' is refused (exit 1, nothing written). Negative test (.scratch/postrh2/negtest.py, dry run on copies in a scratch mirror): a node line already carrying 'If ROOT holds' is refused (exit 1, nothing written). Watchdog peaks: 0 MB reported (runs < 1 s, between samples; stdlib only).
- Caveats: (1) TREE.md moves under us (the ROOT line grew 20,286 → 22,232 chars during this recheck); apply.py re-checks at run time, the log is a snapshot (TREE.md sha 1f0f30a14736d2e5, 21:39 UTC; every dry run since 21:32 passed on each new TREE.md). (2) If the orchestrator rewrites a node's opening words (anchor), that row fails loudly; fix the anchor in DRAFT.md §1. (3) Pass --date with the verification date. (4) EKN do not describe how superfluous candidates are removed (l.755); completeness needs only the superset, so this is residual, not a gap.


### Dry run, verbatim (logs/apply_dryrun.log, part 1 of 2)
```text
apply.py DRY RUN  2026-11-20 21:39 UTC
DRAFT.md sha 1c78bbad737396f3; TREE.md sha 1f0f30a14736d2e5 (582 lines); COROLLARIES.md sha 5d2b8b6eb133b1c0 (81 lines)
TREE.md: 14 insertions, +2042 chars, line count unchanged
  l.5   N1    count 1, after anchor: + 'If ROOT holds: Q(f) ≥ 0 for every f and every Dirichlet χ (Weil), > 0 for compac'
  l.6   N2    count 1, after anchor: + 'If ROOT holds: Λ = 0 (Rodgers–Tao 1801.05914), Λ_χ = 0 for real primitive χ, Λ_K'
  l.7   N3    count 1, after anchor: + 'If ROOT holds: λ_n > 0 for all n and Re λ_n(χ) > 0 for every primitive χ (Li; La'
  l.8   N4    count 1, after anchor: + "If ROOT holds: Robin for n > 5040, Lagarias' H_n form, Nicolas for all primorial"
  l.9   N5    count 1, after anchor: + 'If ROOT holds: d_N → 0; liminf d_N² log N ≥ C; L^p density for p < 2; d_N² ~ C/l'
  l.10  N6    count 1, after anchor: + 'If ROOT holds: every J^{d,n} hyperbolic (Pólya), also for Ξ_t, t ≥ 0, and Ξ_χ of'
  l.11  N7    count 1, after anchor: + 'If ROOT holds: no Dirichlet-L zero in Re s > 1/2; QRT PR 144, Dirichlet half-pla'
  l.12  N8    count 1, after anchor: + 'If ROOT holds: all zeros on the line, ≥ 19/27 simple, ≥ 84.86% distinct; DPC gap'
  l.13  N9    count 1, after anchor: + "If ROOT holds: Lindelöf (13/84 moot), GLH for Dirichlet L, Harper's moments; Sim"
  l.14  N10   count 1, after anchor: + "If ROOT holds: Schoenfeld's bounds, a prime in [x, x + (22/25)√x log x] for x ≥ "
  l.15  N11   count 1, after anchor: + "If ROOT holds: GRH for every Dirichlet L and abelian ζ_K; Kerr's Siegel-zero exc"
  l.16  N12   count 1, after anchor: + 'If ROOT holds: nothing beyond degree 1; GL(2), elliptic-curve and Selberg-class '
  l.17  N13   count 1, after anchor: + "If ROOT holds: no operator follows; QW_λ ≥ 0, but CCM's simple ground state and "
  l.18  ROOT  count 1, after ROOT tag '[check: every decisive outside check pas'…: + 'If verified: RH and Dirichlet GRH proved; consequences in [CONSEQUENCES](Riemann'
COROLLARIES.md: append heading '## If ROOT holds (Alpöge, verified 20 Nov 2026)' and 42 lines (+8252 chars)
  + '- If ROOT holds: Λ = 0 (Rodgers–Tao [1801.05914](https://arxiv.org/abs/1801.0591'
  + '- If ROOT holds: Λ_χ = 0 for every primitive real χ and Λ_K = 0 for every abelia'
  + "- If ROOT holds: Weil's Q(f) ≥ 0 for every test function and every Dirichlet χ ("
  + '- If ROOT holds: λ_n > 0 for all n and Re λ_n(χ) > 0 for every primitive χ, all '
  + '- If ROOT holds: our Li-remainder theorems: mean square between N^{3/2}√log N an'
  + "- If ROOT holds: Robin for every n ≥ 5041 and Lagarias' σ(n) ≤ H_n + e^{H_n} log"
  + '- If ROOT holds: σ(n) < (e^γ/2) n log log n for every odd n > 277093800961598968'
  + '- If ROOT holds: d_N → 0 ([math/0202141](https://arxiv.org/abs/math/0202141) l.4'
  + '- If ROOT holds: natural-number L^p Nyman–Beurling density for every 1 < p < 2, '
  + '- If ROOT holds: D(λ)² log(1/λ) → Σ_ρ m_ρ²/|ρ|² (continuous Nyman–Beurling), and'
  + '- If ROOT holds: every Jensen polynomial J^{d,n} of ξ is hyperbolic (Pólya, [190'
  + '- If ROOT holds: at least 19/27 of the zeros of ζ are simple (Bui–Heath-Brown [1'
  + '- If ROOT holds: at least 84.86% of the zeros of ζ are distinct, and on average '
  + '- If ROOT holds: a positive proportion of normalized gaps < 1/2 + ε, and λ > 3.6'
  + '- If ROOT holds: log|ζ(1/2+it)| ≤ ((log 2)/2 + o(1)) log t/log log t (Chandee–So'
  + '- If ROOT holds: |ζ(1/2+it)| ≤ exp((0.5 + 6.992/((log t)^{0.252} log log t)) log'
  + '- If ROOT holds: ∫_T^{2T}|ζ(1/2+it)|^{2k} dt ≪_k T log^{k²} T for every fixed k '
```
### Dry run, part 2 of 2
```text
  + "- If ROOT holds: Schoenfeld's |ψ(x) − x| ≤ (√x/8π)log²x (x > 59) and |π(x) − li("
  + '- If ROOT holds: a prime in [x, x + (22/25)√x log x] for x ≥ 4 (CMS [1708.04122]'
  + "- If ROOT holds: Brun's constant B < 2.1594 (Dunn [2504.15658](https://arxiv.org"
  + '- If ROOT holds: M(x) ≪ √x exp((log x)^{1/2}(log log x)^{5/2+ε}) (Soundararajan '
  + '- If ROOT holds: the least prime p ≡ a mod q is ≤ (φ(q) log q)² for every q > 3 '
  + '- If ROOT holds: the least quadratic non-residue mod a prime q ≥ 5 is below (log'
  + '- If ROOT holds: for q ≥ 3000 the primes ≤ (log q + B(q))² generate (ℤ/q)* (LLS '
  + '- If ROOT holds: 1/|L(1,χ)| ≤ (12e^γ/π²)(log log q − log 2 + 1/2 + …) and h(−q) '
  + "- If ROOT holds: EKN's lists of the 17 and 27 imaginary quadratic fields with cl"
  + '- If ROOT holds: every odd composite n has a strong Miller–Rabin witness ≤ 2 log'
  + "- If ROOT holds: ĝ(p) < √p − 2 for p > 2791 and Grosswald's g(p) < √p − 2 for p "
  + '- If ROOT holds: ĝ(p) ≤ ((8/5)(2^{ω(p−1)} − 1) log p)² for every prime p ≥ 10⁹ ('
  + '- If ROOT holds: elite primes E(x) ≪ x^{3/4}(log x)^{1/2} ([2611.19610](https://'
  + '- If ROOT holds: #{l : C(k,l) ≥ 0} ≪ φ(k)/log k for the Alladi–Johnson constants'
  + "- If ROOT holds: Hua's PageRank bound Δ_q(1 − c/log q) ≍ 1/N_q on Lubotzky–Phill"
  + '- If ROOT holds: 2s-th moments of short Legendre-symbol sums over primes with er'
  + '- If ROOT holds: long quadratic character sums with the log₃ gain, √X/x ([2509.0'
  + '- If ROOT holds: least Elkies prime ≤ (2 log 4q + 4)² and its counting asymptoti'
  + '- If ROOT holds: Kim [2412.19975](https://arxiv.org/abs/2412.19975) Thm 1.2, who'
  + '- If ROOT holds: the conditional exponents for Gauss mean values and the sphere '
  + '- If ROOT holds: Λ^{(m)}(1/2, χ) ≥ 0 for every m ≥ 0 and every non-principal pri'
  + '- If ROOT holds: the near-null spectrum statement for the windowed Weil form (CO'
  + '- If ROOT holds: every large even N (effective N₀) is a + b with λ(a) = λ(b) = −'
  + "- If ROOT holds: Gaillard–Jindal's GRH reduction of SSR-SLP to one dimension ([2"
  + "- If ROOT holds: Martin's complete lists of degree-2 and -3 rationality fields o"
ALL ASSERTIONS PASSED
DRY RUN: nothing written. Apply with: python3 Riemann/Deep/PostRH/apply.py --write [--date '<d Mon YYYY>']
exit 0; memcap: peak 0 MB rc=0
```

### Files (postrh2)
- Riemann/Deep/PostRH/apply.py: the one-command apply script (dry run by default; --write; --preview DIR).
- Riemann/Deep/PostRH/logs/apply_dryrun.log: its dry-run output (pasted above).
- DRAFT.md: §1 10 sentences updated, ROOT row placed after the tag; §2 42 lines; §3 settled with file and line numbers.
- CONDITIONAL.md: a closing note on the five settled rows.
- Scratch, not deliverables: .scratch/postrh2/ (STATE.md, tex/ with the seven e-prints' .tex/.bbl, anchors.py, nodediff.py, kw.py, draft1.py, draft2.py, previewcheck.py, negtest.py); the preview was deleted after the check (recreate: apply.py --preview DIR).
