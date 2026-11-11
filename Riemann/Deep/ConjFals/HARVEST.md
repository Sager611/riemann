# ConjFals harvest — numerically testable conjectures on the RH tree (arXiv 2023-2026)

Method: 98 arXiv API queries (abstract "conjecture(s)/conjectured" + node keywords; title queries on node objects) → 441 papers since 2023 (scripts and query lists in `harvest/`; the downloaded abstracts and e-prints were deleted after the run); 105 TeX e-prints scanned for conjecture/question/problem environments and conjecture phrases (`harvest/scan.py`); Zenodo/HAL keyword queries (AI-generated RH notes; nothing testable beyond the kernel-level claims below); GXL evidence API for the Redheffer item. Ranked by (chance of a counterexample in reach) × (RH relevance). TeX lines refer to the e-print's main file.

| # | Source (TeX line) | Statement | Node | Status here |
|---|---|---|---|---|
| 1 | Jerby 2406.06548 conj:rep (l. 313, 777); 2310.14415 Conj "Repulsion G-B-G" (l. 266) and Conj G-B-G (l. 753); 2610.21952 (l. 323: "to n ~ 2.5·10^7") | g_n bad ((−1)^n Z(g_n) < 0), g_{n±1} good ⇒ μ(g_n) = \|Z'(g_n)/Z(g_n)\| > 4 | N8 (Gram/Rosser zero counting), N13 (proposed repulsion mechanism) | tested, holds (REPORT §1) |
| 2 | Johnston–Leong–Tudzi 2408.04143 Conj newcon2 (l. 196), Conj (l. 1142), Sun's Conj (l. 128) | \|W_a(x)\| < x^{log₂ a} for real a ≥ 2, x ≥ x₀(a), W_a = Σ_{n≤x}(−a)^{Ω(n)}; lim_{a→∞} s_a = 1; \|W(x)\| < x for x ≥ 3078 | N10 (Mertens/Pólya-type sums) | tested, holds; lim s_a = 1 proved (REPORT §2) |
| 3 | Barrett–Jarvis (LAA 1992), restated 2511.13627 (l. 524) | "small" eigenvalues of the Redheffer matrix have negative real part | N10 (det R_n = M(n)) | restatement false at n = 300; real conjecture disproved by Cardon, arXiv 0809.0076 (REPORT §3) |
| 4 | Farr–Pauli–Saidak (Funct. Approx. 64, 2021), restated 2301.11209 (l. 523–525) | \|C_α(a)\| ≤ 2\|e^{α(log w − 1/w)}\|, w = W₀(iα/2π), C_α(a) = γ_α(a) − log^α a/a | N3 (Stieltjes constants) | tested, holds (REPORT §4) |
| 5 | Michałowski 2602.20313 Open Problem 1 (l. 820) | K(u) = Φ(\|u\|) is PF₄ for all configurations (PF₅ fails) | N2 (de Bruijn–Newman kernel) | tested, no violation (REPORT §5) |
| 6 | Henry 2601.00631 Conj (l. 286) = 2601.09276 Conj (l. 113) | 0 < x < 1 ⇒ π cot πx + x < ζ(x) − ψ(x) < π cot πx + (γ+½)x + γ − ½ | none (ζ on (0,1)) | tested, holds (REPORT §6) |
| 7 | Spira, restated 2405.12557 (l. 164) | all non-trivial zeros of Z_{[t/2]}(t) are real | N8 | statement vague; Lehmer pair test only (REPORT §7) |
| 8 | 2511.02106 Remark (l. 1343) | Lemma x₂ < √(2x) (proved for κ ≥ 3/2) for all κ > 1 | N4 | tested, holds (REPORT §8) |
| 9 | Ihara via 2402.13127 (l. 239); 2407.09113 (l. 532) | (a₀−ε) log m < γ_{Q(ζ_m)} < (a₁+ε) log m; κ(q) = o(1) | N11 | asymptotic, not finitely falsifiable |
| 10 | 2608.24827 Conj (l. 1926) | −ln λ_min(L) ~ 2π² N(T*)/ln N(T*), T* = 2πe^{2L} | N1 | asymptotic; c = 47 check: prediction 542.6 vs certified 556.6 |
| 11 | 2311.04860 ELI (l. 159) | \|Σ ℓ_γ γ\| ≥ C_ε e^{…}, C_ε unspecified | zeros (LI) | not falsifiable |
| 12 | 2602.03408 (l. 841–1207); 2606.22562 (l. 1155–1281); 2409.17708 (l. 218); 2312.00211; 2312.12760; 2609.33794 (l. 420); 2310.14415 (l. 384, 1354, 1444) | limits "except countably many a", RH-conditional statements, equivalences, open questions, path-existence claims | — | not finite statements |

Found already refuted by others (not ours): Berndt–Straub's S_{2n} conjecture (2610.18065, χ₁₆₃); Barrett–Jarvis (Cardon 0809.0076); proved: Ramanujan-type polynomial conjecture (2306.10283).

Skipped as done here: 2511.13691 (R268), 2602.15905 (R213), 2609.14295 (R332), 2610.12816 (R177), 2611.13606 (R422).
