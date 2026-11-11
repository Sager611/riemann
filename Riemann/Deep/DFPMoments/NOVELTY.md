# NOVELTY: optimised moment exponents from de Faveri–Pandey (2610.01035 v1)

Checked 6 Nov 2026 and again 11 Nov 2026 (logs: `logs/novelty_*.log`).

| Source | Query | Result |
|---|---|---|
| arXiv listings math.NT, math.CA (+ math.CV on 11 Nov), 2026-10 and 2026-11, `?show=2000` | titles ~ sixth moment, twelfth moment, large value, moments of the Riemann zeta, power moments | 5 hits: 2610.01035 itself; 2610.18923 (Dirichlet polynomials, trilinear sums; abstract: short-interval primes, not zeta moments); 2610.22107, 2610.24958, 2610.36885 off topic |
| arXiv API | id:2610.01035; au:Faveri; "de Faveri" AND Pandey; au:Pandey AND zeta; abs "sixth moment" AND zeta; abs "twelfth moment" | v1 only; no follow-up by the authors; nothing newer than 2610.01035 on sixth or twelfth moments |
| OpenAlex | works citing W7215832989 (= 2610.01035) | 2 records, both 2506.16766 (log-power savings in E_2(T) and kth moments, 8 ≤ k ≤ 12; M_12 ≪ T²(log T)^15; no power saving, no optimisation of DFP) |
| HAL | "sixth moment" AND zeta | one 2023 record (multiple integral formulas), unrelated |
| Zenodo | "sixth moment" "Riemann zeta"; "Faveri" zeta | nothing relevant |
| ANTEDB (teorth/expdb) | commits since 25 Sep 2026; latest 20 PRs and issues; literature.py, zeta_moments.tex, zeta_large_values.tex | last commit 26 Sep 2026; no mention of de Faveri–Pandey; no exponent pair newer than TTY25 |
| Our files | rg 445121, 37/30, de Faveri, sixth/twelfth moment in TREE, NOGO, DISCOVERIES, COROLLARIES, Riemann/Deep | R280, R327, R343, R386 recomputed the printed exponents; R396, R412 list "optimise δ_k and ν" as a harvest candidate, not done; DensityB, QRHImprove, BFlat cite DFP for other purposes |

Verdict: the optimised table and curve for 6 < A < 12 and the explicit range ν = 131/288 are not in the literature we can
see. They are an optimisation the authors explicitly left undone (TeX l.74: "We do not optimize either the savings δ_k or
the range ν"), using only published bounds; expect it to be regarded as routine.
