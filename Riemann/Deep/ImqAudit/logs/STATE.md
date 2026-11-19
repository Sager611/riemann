# imqaudit STATE
Start 2026-11-19 21:12:46 UTC. Stop-rule writeup at 22:32, hard end 22:42.

- [ ] read page, REPORT, rerun, checks
- [ ] referee M(d) lower bound
- [ ] code audit isieve2.c, threshold.py, mjx.c, t4.gp
- [ ] falsify: exact h for sample d
- [ ] rerun one slice
- [ ] per-h table sample
- [ ] reviewer on delicate step
- [ ] REPORT.md, FULLPROOF.md
21:50 window test done; t4 33 hcount ok; adv 528 ok; c08 seg4 rerun ok
21:54 reviewer: GAP in my paraphrase only (phase-B p>=17 omitted; code has it: pmed starts at 17); CORRECT for the code; small-Xs argument given
21:55 c08 seg1-4 rerun identical to p_c08.log seg 4 line; HJKMP 0 diffs; REPORT sections 1-3 written; now window2 at 1.31e11 N=30000
21:59 REPORT, FULLPROOF drafted; remaining: final read-through, cleanup
22:06 adv2 top-range 1500 hard d: min M 41354; all checks done; final cleanup next
22:08 DONE: AUDIT PASSED; REPORT.md, FULLPROOF.md, code, logs in Riemann/Deep/ImqAudit/; scratch cleaned
