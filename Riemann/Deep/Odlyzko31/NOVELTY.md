# odlyzko31: novelty check (6 Nov 2026)

Question: has anyone derived unconditional root-discriminant lower bounds (Odlyzko-Poitou type) from the 2026 zero-free strips (Cureton-Thorner's 31/32 for every zeta_K, the Dirichlet half-planes, OpenAI's 7/8 for zeta, Liu's Hecke half-planes)?

## Searches
| Source | Query | Result |
|---|---|---|
| arXiv math.NT listings 2026-10, 2026-11 (`.scratch/lst.py`, titles) | discriminant, Odlyzko, Poitou, root discriminant | 4 title hits: 2610.13295, 2610.16505 (Hecke modules of ternary lattices), 2610.19597 (Erdős-Odlyzko-Sárközy), 2610.26510 (density of A4-quartic discriminants); none relevant |
| arXiv API (all years, newest first) | all:"root discriminant"; abs:Odlyzko AND abs:discriminant; abs:Cureton; abs:"discriminant bounds"; abs:Poitou AND abs:discriminant; abs:"quasi-Riemann"; abs:"quasi-Riemann" AND abs:discriminant; abs:"Odlyzko bound(s)"; abs:"zero-free" AND abs:"root discriminant"; abs:"Dedekind zeta" AND abs:"zero-free strip"; abs:discriminant AND abs:Dedekind AND abs:unconditional | no paper on discriminant lower bounds from a zero-free strip; 2026 hits are applications of the strips to other problems (2611.06705, 2610.37179, 2610.26835, 2610.24900, 2610.20600, 2610.18491, 2610.15989, 2610.20058, 2611.02280) |
| TeX of the closest 2026 papers | Cureton-Thorner (Hexagon 2610.00200v2, the source of the 31/32 strip); Frei-Sofos-Widmer 2610.18491 (applications of quasi-RH to class groups); 2610.20600 (small defining polynomials; cites Cureton-Thorner) | no discriminant lower bound; "Odlyzko" occurs only as Lagarias-Odlyzko (Chebotarev); 2610.20600 uses Minkowski's second theorem only |
| Zenodo API | "discriminant bounds"; "Odlyzko bounds"; "root discriminant" "quasi-Riemann"; "discriminant" "Cureton"; "Odlyzko" "zero-free"; "root discriminant" "unconditional"; "discriminant bounds" "number field" | nothing relevant (unit-distance, Fermat, Malle-type records; Zenodo 23605077 "Zero-Free Boundary Strips ..." has no discriminant content) |
| HAL API | "root discriminant"; "Odlyzko"; "discriminant bounds"; "quasi-Riemann" | nothing after 2021 on discriminant bounds |
| Our files | rg Odlyzko/Poitou/root discriminant in Riemann/Deep, DISCOVERIES, COROLLARIES, NOVELTY files | only Harvest4 candidate 8 (this task) |

## Before (the comparison baseline)
- Unconditional: Odlyzko (1976 tables; Invent. Math. 1975, Acta Arith. 1976), Poitou (Sém. Delange-Pisot-Poitou 1976/77): rd_K >= (60.8395 - o(1))^{r1/n} (22.3816 - o(1))^{2 r2/n}; constants 4 pi e^{1+gamma}, 4 pi e^gamma, the suprema of the method when zeros may lie anywhere in the critical strip (DERIVATION §3, w = 1/2). Diaz y Diaz (Publ. Math. Orsay 1980) tabulated optimised unconditional values for n <= 100 (not available to us; Odlyzko's own remark, "D^{1/n} >= 5.656 for n >= 8 with better kernels", is reproduced by our classical kernel optimisation, 5.6563).
- GRH: 44.7632, 215.3325 (Odlyzko, Serre).
- Zero-free regions near s = 1 (de la Vallée Poussin type, Stark) do not change the asymptotic constants; no earlier unconditional strip of fixed width for every zeta_K existed before Cureton-Thorner (Oct 2026).

## Verdict
No prior or concurrent derivation found. The improvement is a direct consequence of a 2026 input (a "routine" substitution in the explicit-formula method), but the numbers (limits 23.3878 / 65.6439, crossovers 2732 / 1618, the table) are new. The optimality statement for general strip width w (every admissible kernel satisfies F <= sech(wx), so A_c(w), A_r(w) are the method's exact limits) was not found stated anywhere; for w = 1/2 and w = 0 it is the classical description of the Odlyzko-Poitou and Serre constants.
