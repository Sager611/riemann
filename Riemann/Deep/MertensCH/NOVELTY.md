# NOVELTY — mertensch (2026-11-11)

Question: is there an explicit bound |M(x)| ≤ c x^θ (θ < 1, all x ≥ 1) better than our Mertens78 family, or the hybrid method for M(x)?

| source | query / scope (fetched 2026-11-06 and 2026-11-11) | result |
|---|---|---|
| arXiv math.NT listings 2026-10 and 2026-11 (lst.py) | titles ~ Mertens, Möbius, reciprocal of zeta, 1/zeta, zero-free, half-plane, explicit (33 hits) | no explicit M(x) power saving; neighbours: 2610.12234 (half-planes for QRH), 2610.26835 (Hecke L half-planes), 2610.32012 (|ζ_K(1+it)|), 2610.09567 and 2610.11487 (Möbius in short intervals, not explicit), 2611.10499 (unitary Möbius, Goldbach type) |
| arXiv API (axq.py) | abs:"Mertens function" AND (explicit OR bound); abs:Mertens AND abs:explicit; abs:"Mertens function" AND abs:"zero-free" | only 2208.06141 (sub-power savings), 2511.14736 (CH), 2608.01498 (Mertens sums, seen in R379), 2607.09797, 2607.04366 (not bounds for M); last query empty |
| arXiv versions | 2511.14736, 2512.15709 | both still v1 (2025-11-18, 2025-12-17) |
| HAL API | "Mertens function"; Möbius AND explicit AND bound | hal-05371887 (= CH 2511.14736), hal-05838582 (Bi–Hou counterexample bound), nothing else relevant |
| Zenodo API | "Mertens function" | 23700212 (counterexample bound), others unrelated; no explicit power saving |
| our files | DISCOVERIES.md, discoveries/mertens-power-saving.md, Mertens78/REPORT.md §5, TREE N10 | the record being improved is ours (Mertens78, 5 Nov 2026) |

Method novelty: moving the line of integration to Re s = σ_L for |t| ≤ T where RH is verified, and keeping a zero-free region above T, is standard for ψ(x) (e.g. Büthe 2016/2018, Fiori–Kadiri–Swidinsky 2023, Fiori–Jaskari 2609.23222). For M(x) it had no power-saving partner before Mertens78, so the combination is new but not a new technique. Chirre–Helfgott's weight is band-limited; their linear term cannot be removed with a zero-free half-plane (DERIVATION.md §1); our bound replaces their method by the hybrid contour beyond x ≈ 10^245, it does not refine their theorem.

Verdict: new as explicit constants (improves our own record family at every θ in [0.8992, 0.9409]); the method is a standard hybrid. The part due to RH to 3·10^12 is 8× at θ = 0.9409, 6.1× at 0.9382, 2.1× at 0.93, at most 1.42× for θ ≤ 0.9277; the rest (Chirre–Helfgott below 10^245, re-optimised smoothing, squarefree unsmoothing, certificates moved to the Lean half-plane 0.87338…) is routine.
