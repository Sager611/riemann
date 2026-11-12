# STAudit checks (2026-11-12; E = 11a1; all runs under `.scratch/memcap.sh 290`, single-threaded; peak 208 MB (msq.gp))
- `msq.gp` → `msq.out` (|D| ≤ 400), `msq2.gp` → `msq2.out` (401 ≤ |D| ≤ 600): columns D, h, #coefficients, lfuncheckfeq bits, [Σ_χ|L(1+δ+it)|² for (δ,t) = (0,0),(0,1.5),(0,4),(.05,0),(.05,1.5),(.05,4),(.15,0),(.15,1.5),(.15,4)] (msq2: t ∈ {0, 1.5} only), ms. Arithmetic normalisation (centre 1). 11 | D skipped. Analysis: `msq_an.py`, `msq_an2.py`.
- `fam2.c`, `fam3.c` → `fam2.out`, `fam3.out`, `fam3_beyond.out`: usage `./fam3 X N/X ...` with `ap.txt` from `apdump.gp`; R1 = F/G (family sum over its pure diagonal); modes 0 λ_E, 1 λ_{Sym²E}, 2 pole, 3 a_Vλ_E n^{−0.55}, 4 a_Vλ_{Sym²} n^{−0.8}.
- `sign.c` → `sign.out`: usage `./sign Dlo Dhi 2.0365` (needs `ap.txt`).
- `thr_check.py` (23 exact sympy checks, 0 failures), `thr_opt2.py` (free-parameter optimiser for m = 2: 2.0365), `minor.py` (Christoffel minorants, 4812 intervals).
- `reviewer_A.md`, `reviewer_B.md`: reviewer briefs (quotes verbatim); `reviewer_verdicts.md`: the three verdicts.
