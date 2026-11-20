READY (14 node consequences, 45 UNCONDITIONAL, 31 STAYS CONDITIONAL, 9 open problems): if Alpöge's five Lean statements are verified, every conditional result that needs only RH for ζ or GRH for Dirichlet L becomes a theorem (5 more after one named check each), while everything needing non-abelian, Hecke, Artin or GL(2) GRH, or simple zeros, stays conditional.

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
