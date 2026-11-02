NO RECORD (2026-10-23, stopped 2026-11-02): pushing past half-width 1.55 needs an out-of-core rewrite; the existing code peaked at 4-5 GB at L = 1.55 (grams ~6750 s, certificate ~3300 s CPU), far above the 300 MB per-process limit. Agent notes follow.

# weilrec5 STATE (Orca agent; brief .scratch/orca/weilrec5-brief.md)
Goal: certify Weil positivity half-width > 1.55 (aim 1.6) under: 1 single-threaded proc, <=1 h/run, <300 MB.
## Log
- 2026-10-23 03:10 start. Read PROOF.md, REPORT3.md, cert4.py, grams_c.c, run4 logs. M4 (4P+6E), we run on E-cores (taskpolicy -b).
- 1.55 cost: grams ~6750 s user, cert 3300 s user, 4-5 GB peak -> existing code cannot run under 300 MB at any L. Need lean out-of-core rewrite.
