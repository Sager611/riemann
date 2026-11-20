# DRAFT — for the orchestrator, if ROOT (TREE.md l.18) is verified

Exact-substring insertions; nothing here edits repo files. Every anchor below occurs exactly once in TREE.md (rg -F -c, 20 Nov 2026 19:03 UTC; rechecked by postrh2 21:2x UTC at HEAD 1259790d, REPORT.md "Recheck"). Apply with Riemann/Deep/PostRH/apply.py (dry run by default). Sources and verdicts: CONDITIONAL.md (rows A*, B*), CONSEQUENCES.md.

## 1. TREE.md node lines: insert "<sentence> " immediately after the anchor (ROOT: after the anchor's leading [ … ] tag, which stays untouched)
| anchor (exact) | sentence (≤ 25 words) |
|---|---|
| `- N1 Weil/Connes positivity — ` | If ROOT holds: Q(f) ≥ 0 for every f and every Dirichlet χ (Weil), > 0 for compactly supported f ≠ 0; window records become quantitative only. |
| `- N2 de Bruijn–Newman Λ — ` | If ROOT holds: Λ = 0 (Rodgers–Tao 1801.05914), Λ_χ = 0 for real primitive χ, Λ_K = 0 for abelian K (Dobner 2005.05142); our Λ, Λ_χ records moot. |
| `- N3 Li/Keiper — ` | If ROOT holds: λ_n > 0 for all n and Re λ_n(χ) > 0 for every primitive χ (Li; Lagarias math/0404394); range records moot, PR 144's included. |
| `- N4 Robin/Lagarias — ` | If ROOT holds: Robin for n > 5040, Lagarias' H_n form, Nicolas for all primorials, our odd Robin beyond N₀; range, 2-adic, t-free records moot. |
| `- N5 Nyman–Beurling/Báez-Duarte — ` | If ROOT holds: d_N → 0; liminf d_N² log N ≥ C; L^p density for p < 2; d_N² ~ C/log N needs BCF's hypotheses. |
| `- N6 Jensen-polynomial hyperbolicity — ` | If ROOT holds: every J^{d,n} hyperbolic (Pólya), also for Ξ_t, t ≥ 0, and Ξ_χ of real primitive χ; degree records moot. |
| `- N7 Zero-free regions and zero-density constants — ` | If ROOT holds: no Dirichlet-L zero in Re s > 1/2; QRT PR 144, Dirichlet half-planes, densities, verification heights moot; Hecke/non-abelian ζ_K regions stay live. |
| `- N8 Critical-line proportion — ` | If ROOT holds: all zeros on the line, ≥ 19/27 simple, ≥ 84.86% distinct; DPC gaps < 1/2 + ε, λ > 3.638 beat records 67.36%, 7/11, 3.367; simplicity open. |
| `- N9 Moments, Lindelöf, subconvexity — ` | If ROOT holds: Lindelöf (13/84 moot), GLH for Dirichlet L, Harper's moments; Simonič's explicit RH bound beats every explicit t^{1/6} bound for t ≥ 10^50.81. |
| `- N10 Explicit ψ and π − li bounds — ` | If ROOT holds: Schoenfeld's bounds, a prime in [x, x + (22/25)√x log x] for x ≥ 4; 7/8, x^0.856 and Mertens-exponent records moot. |
| `- N11 GRH families and Siegel zeros — ` | If ROOT holds: GRH for every Dirichlet L and abelian ζ_K; Kerr's Siegel-zero exclusion superseded; LLS's explicit Linnik, non-residue, Littlewood; non-abelian, Artin, Hecke stay open. |
| `- N12 Function-field and Selberg-class analogues — ` | If ROOT holds: nothing beyond degree 1; GL(2), elliptic-curve and Selberg-class GRH stay open; our linear-combination and root-discriminant records stand. |
| `- N13 Hilbert–Pólya and spectral approaches — ` | If ROOT holds: no operator follows; QW_λ ≥ 0, but CCM's simple ground state and self-adjointness stay open. |
| `- ROOT (RH itself) — ` | If verified: RH and Dirichlet GRH proved; consequences in [CONSEQUENCES](Riemann/Deep/PostRH/CONSEQUENCES.md); loop re-aims at GRH beyond Dirichlet ([OPEN](Riemann/Deep/PostRH/OPEN.md)). |


## 2. COROLLARIES.md block (one line per UNCONDITIONAL result; each conditional on ROOT's verification; A4, A20 omitted as already unconditional or superseded, B3's 0.6845 superseded by 19/27, A12 + A14 and A18 + A19 share a line; postrh2 added B16, B3 (0.8486, 0.9350), D26, A13 and replaced A27's finiteness line, already unconditional, by Martin's complete lists)
- If ROOT holds: Λ = 0 (Rodgers–Tao [1801.05914](https://arxiv.org/abs/1801.05914) Thm, l.129, with RH ⟺ Λ ≤ 0, l.68); our Λ ≤ 0.1074 moot (postrh B1).
- If ROOT holds: Λ_χ = 0 for every primitive real χ and Λ_K = 0 for every abelian field K (Dobner [2005.05142](https://arxiv.org/abs/2005.05142) l.113, 119: Λ_F ≥ 0 on S♯; heat flow l.236–240, as in R324; our Λ_χ and abelian Λ_K bounds, COROLLARIES R324 and later, moot; the all-field Λ_K ≤ 225/512 stays; postrh B16).
- If ROOT holds: Weil's Q(f) ≥ 0 for every test function and every Dirichlet χ (Weil's criterion), Q(g∗g̃) > 0 for compactly supported g ≠ 0 (ĝ has exponential type, so O(T) zeros, against N(T) ~ (T/2π) log T; postrh CONSEQUENCES N1; R542).
- If ROOT holds: λ_n > 0 for all n and Re λ_n(χ) > 0 for every primitive χ, all n (Lagarias [math/0404394](https://arxiv.org/abs/math/0404394) l.1096; postrh B15).
- If ROOT holds: our Li-remainder theorems: mean square between N^{3/2}√log N and N^{3/2}√log N·log log N, limsup \|R(n)\|/(n log n)^{1/4} ≥ 0.289 ([page](discoveries/li-remainder.md); postrh A15).
- If ROOT holds: Robin for every n ≥ 5041 and Lagarias' σ(n) ≤ H_n + e^{H_n} log H_n for every n ≥ 1 ([math/0008177](https://arxiv.org/abs/math/0008177) l.92–97, 171–176; postrh B12).
- If ROOT holds: σ(n) < (e^γ/2) n log log n for every odd n > 277093800961598968947975375 ([page](discoveries/odd-robin.md); postrh A16).
- If ROOT holds: d_N → 0 ([math/0202141](https://arxiv.org/abs/math/0202141) l.47) and liminf d_N² log N ≥ 2 + γ − log 4π (Burnol math/0103058 Thm 1.3: ≥ Σ m_ρ²/|ρ|² ≥ Σ m_ρ/|ρ|² = 2 + γ − log 4π under RH; postrh B13; R542).
- If ROOT holds: natural-number L^p Nyman–Beurling density for every 1 < p < 2, for ζ, every Dirichlet L and every abelian ζ_K ([page](discoveries/nyman-beurling-lp.md); abelian ζ_K by NBLpK's Theorem A_K, COROLLARIES R417 line; postrh A17).
- If ROOT holds: D(λ)² log(1/λ) → Σ_ρ m_ρ²/\|ρ\|² (continuous Nyman–Beurling), and the analogue for primitive χ mod q ≥ 3 (DISCOVERIES N5; COROLLARIES nbdir; postrh A12, A14).
- If ROOT holds: every Jensen polynomial J^{d,n} of ξ is hyperbolic (Pólya, [1902.07321](https://arxiv.org/abs/1902.07321) l.277), and of Ξ_χ for primitive real χ (postrh B14, CONSEQUENCES N6).
- If ROOT holds: at least 19/27 of the zeros of ζ are simple (Bui–Heath-Brown [1302.5018](https://arxiv.org/abs/1302.5018) Thm, l.73–75, "Assuming RH"; postrh B2).
- If ROOT holds: at least 84.86% of the zeros of ζ are distinct, and on average 93.50% of the zeros of primitive Dirichlet L-functions are simple (Chirre–Gonçalves–de Laat [1810.08843](https://arxiv.org/abs/1810.08843) l.257 'GRH for Dirichlet L-functions', l.302–310, 359–362; its 0.6845 simple is below 19/27; postrh B3).
- If ROOT holds: a positive proportion of normalized gaps < 1/2 + ε, and λ > 3.63846723 (DPC [2610.28527](https://arxiv.org/abs/2610.28527) l.62, l.76; postrh A18–A19).
- If ROOT holds: log\|ζ(1/2+it)\| ≤ ((log 2)/2 + o(1)) log t/log log t (Chandee–Soundararajan [0908.2008](https://arxiv.org/abs/0908.2008) l.100–103; postrh B6).
- If ROOT holds: \|ζ(1/2+it)\| ≤ exp((0.5 + 6.992/((log t)^{0.252} log log t)) log t/log log t) for t ≥ 10^40 (Simonič [2010.13307](https://arxiv.org/abs/2010.13307) Thm main; postrh B7).
- If ROOT holds: ∫_T^{2T}\|ζ(1/2+it)\|^{2k} dt ≪_k T log^{k²} T for every fixed k ≥ 0 (Harper [1305.4618](https://arxiv.org/abs/1305.4618) l.87; postrh B8).
- If ROOT holds: Schoenfeld's \|ψ(x) − x\| ≤ (√x/8π)log²x (x > 59) and \|π(x) − li(x)\| ≤ (√x/8π)log x (x > 2657) (Büthe [1410.7015](https://arxiv.org/abs/1410.7015) l.58, 453; B9).
- If ROOT holds: a prime in [x, x + (22/25)√x log x] for x ≥ 4 (CMS [1708.04122](https://arxiv.org/abs/1708.04122) l.368), hence between consecutive k-th powers, k ≥ 3 (routine; B10).
- If ROOT holds: Brun's constant B < 2.1594 (Dunn [2504.15658](https://arxiv.org/abs/2504.15658) Thm l.43; its only GRH input is Johnston–Starichkova [2208.01229](https://arxiv.org/abs/2208.01229) Lemma 5.1, l.621–638: Grenié–Molteni for ℚ(ζ_q)/ℚ, a product of Dirichlet L; postrh D26).
- If ROOT holds: M(x) ≪ √x exp((log x)^{1/2}(log log x)^{5/2+ε}) (Soundararajan [0705.0723](https://arxiv.org/abs/0705.0723) Thm 1; Balazard–de Roton [0810.3587](https://arxiv.org/abs/0810.3587); postrh B11).
- If ROOT holds: the least prime p ≡ a mod q is ≤ (φ(q) log q)² for every q > 3 (LLS [1309.3595](https://arxiv.org/abs/1309.3595) l.217–220; postrh B17).
- If ROOT holds: the least quadratic non-residue mod a prime q ≥ 5 is below (log q)² (LLS [1309.3595](https://arxiv.org/abs/1309.3595) l.158–160; postrh B18).
- If ROOT holds: for q ≥ 3000 the primes ≤ (log q + B(q))² generate (ℤ/q)* (LLS [1309.3595](https://arxiv.org/abs/1309.3595) l.135–145: one lies outside each proper subgroup; B19, A24).
- If ROOT holds: 1/\|L(1,χ)\| ≤ (12e^γ/π²)(log log q − log 2 + 1/2 + …) and h(−q) ≥ (π/(12e^γ))√q(…)^{−1} (−q fundamental) for q ≥ 10^10 (LLS l.243–275; postrh B20).
- If ROOT holds: EKN's lists of the 17 and 27 imaginary quadratic fields with class group of exponent 3 and 5 are complete, and so are ARSX's exponent-3 and -5 tables over all negative discriminants ([1803.02056](https://arxiv.org/abs/1803.02056) Thm 1 l.148–177: computed to 3.1·10^20 without hypothesis, above the ERH bounds 9.7·10^10, 2.3·10^20 of l.268–275; [2605.26020](https://arxiv.org/abs/2605.26020) l.262–267, 666; postrh A13).
- If ROOT holds: every odd composite n has a strong Miller–Rabin witness ≤ 2 log² n (Bach 1990, Math. Comp. 55; postrh B22).
- If ROOT holds: ĝ(p) < √p − 2 for p > 2791 and Grosswald's g(p) < √p − 2 for p > 409 (MTT [1508.05182](https://arxiv.org/abs/1508.05182) l.45; B21).
- If ROOT holds: ĝ(p) ≤ ((8/5)(2^{ω(p−1)} − 1) log p)² for every prime p ≥ 10⁹ (MTT [1508.05182](https://arxiv.org/abs/1508.05182) l.76; postrh A23).
- If ROOT holds: elite primes E(x) ≪ x^{3/4}(log x)^{1/2} ([2611.19610](https://arxiv.org/abs/2611.19610) Thm B; GRH only via Dirichlet character sums; postrh A1).
- If ROOT holds: #{l : C(k,l) ≥ 0} ≪ φ(k)/log k for the Alladi–Johnson constants ([2610.23660](https://arxiv.org/abs/2610.23660) Thm grh; postrh A2).
- If ROOT holds: Hua's PageRank bound Δ_q(1 − c/log q) ≍ 1/N_q on Lubotzky–Phillips–Sarnak graphs ([2610.00109](https://arxiv.org/abs/2610.00109) Thm main-results (ii), GRH for L(s, χ_D) only; postrh A3).
- If ROOT holds: 2s-th moments of short Legendre-symbol sums over primes with error Q^{1/2+o(1)} ([2604.23661](https://arxiv.org/abs/2604.23661) Thm l.700–705; postrh A5).
- If ROOT holds: long quadratic character sums with the log₃ gain, √X/x ([2509.07651](https://arxiv.org/abs/2509.07651) Thm 1.3, f = 1; postrh A6).
- If ROOT holds: least Elkies prime ≤ (2 log 4q + 4)² and its counting asymptotic, curve assumptions as there (Lippmann–McGown [2311.17231](https://arxiv.org/abs/2311.17231) Thms 1–2; postrh A7).
- If ROOT holds: Kim [2412.19975](https://arxiv.org/abs/2412.19975) Thm 1.2, whose only assumption (l.79) is a Dirichlet zero-free half-plane (postrh A8).
- If ROOT holds: the conditional exponents for Gauss mean values and the sphere problem ([2611.10606](https://arxiv.org/abs/2611.10606) Thm thm:cond; [2610.29636](https://arxiv.org/abs/2610.29636) Thm thm:qrh; postrh A9).
- If ROOT holds: Λ^{(m)}(1/2, χ) ≥ 0 for every m ≥ 0 and every non-principal primitive real χ (Goldfeld–Huang 2017 Thm 1.2; postrh A10).
- If ROOT holds: the near-null spectrum statement for the windowed Weil form (COROLLARIES aclgred entry, "under RH"; postrh A11).
- If ROOT holds: every large even N (effective N₀) is a + b with λ(a) = λ(b) = −1 (Mangerel [2412.17199](https://arxiv.org/abs/2412.17199) Thm l.127, Dirichlet GRH; postrh A25).
- If ROOT holds: Gaillard–Jindal's GRH reduction of SSR-SLP to one dimension ([2304.13605](https://arxiv.org/abs/2304.13605) l.360; Chebotarev only in ℚ(√a); postrh A26).
- If ROOT holds: Martin's complete lists of degree-2 and -3 rationality fields of weight-2 CM newforms with trivial nebentypus, 23 quadratic and 18 cubic ([2509.24119](https://arxiv.org/abs/2509.24119) Thm thm2 l.188; ERH for quadratic characters, l.174, enters only through the exponent-2 and -3 lists, l.620, 665; finiteness for every d is already unconditional, COROLLARIES harvest4 (c); postrh A27).

## 3. Pending checks: settled (postrh2, 20 Nov 2026 21:3x UTC; TeX in .scratch/postrh2/tex/, lines of the arXiv e-print source)
- B16 Λ_χ, Λ_K (Dobner 2005.05142 paper.tex): UNCONDITIONAL-UNDER-ROOT, moved to §2. ξ^F_t((1+iz)/2) = ∫ e^{tu²}Φ_F(u)e^{izu} du (l.236–240), the Rodgers–Tao flow with s = (1+iz)/2, as COROLLARIES.md:18 uses (zeros of H_0^χ in \|Im z\| ≤ 2(B − 1/2), de Bruijn Λ ≤ D²/2); γ unique up to scaling (l.244); Λ_F ≥ 0 (l.119); L(s, χ) real primitive and ζ_K lie in S♯ (l.194–216). Λ = 0 is scale-invariant anyway.
- A13 exponents 3, 5 (EKN 1803.02056 main.tex; ARSX 2605.26020 boundary_CM_pts.tex): UNCONDITIONAL-UNDER-ROOT, moved to §2. EKN Thm 1 (l.148–175) is stated without hypothesis for \|D\| ≤ 3.1·10^20 (split-prime enumeration 4p^c = x² + \|D\|y², l.745–763, and a sieve for no split prime ≤ 193, l.785–794); ERH enters only through Bach–Sorenson's split prime for the quadratic field (Thm Bach_Bound l.250–258), bounds 9.7·10^10 (c = 3) and 2.3·10^20 (c = 5) (l.268–275); l.177. ARSX l.262–267: E = 2, 3, 4, 5, 8 complete under ERH, only input EKN's lists (l.666). Residual: EKN's removal of superfluous candidates (l.755) is not described.
- A13 exponents 6, 7: STAYS CONDITIONAL, on a finite search: EKN's ERH bounds 2.5·10^25, 3.9·10^30 (l.273–274) exceed the 3.1·10^20 computed; ARSX l.268–269 make E = 6, 7 complete only if no fundamental D has exponent 6 in [3.1·10^20, 2.5·10^25], resp. 7 in [3.1·10^20, 3.9·10^30].
- B3 CGdL 1810.08843 (LPBandZETAV_17.tex): UNCONDITIONAL-UNDER-ROOT. l.257 defines 'GRH' as 'the Generalized Riemann Hypothesis for Dirichlet L-functions'; GGOS's input F(x, T) ≥ 3/2 − \|x\| − o(1) (l.494–497) is under it. Its 0.6845 (l.280–283) is below 19/27 (B2), so no line; the same paper's GRH corollaries N_d(T) ≥ (0.8486 + o(1))N(T) (l.302–310, with BHB's 19/27, l.296–300) and N_{Φ,s}(Q) ≥ (0.9350 + o(1))N_Φ(Q) (Thm thmsimplezerosofdir, l.359–362) were not in postrh's rows and now are (§2).
- A27 Martin thm2 (2509.24119 CM-rat.tex): UNCONDITIONAL-UNDER-ROOT, moved to §2 (replacing A27's finiteness line, already unconditional by COROLLARIES harvest4 (c)). ERH is 'for quadratic Dirichlet characters' (l.174); thm2 (l.188) needs the exponent-2 list (Weinberger, l.620; also unconditional, COROLLARIES Exp48) and the exponent-3 list (EKN, l.665), exponents dividing d = 2, 3 (cor:exp l.371); the rest is Größencharacter algebra and Sage (l.264).
- D26 Brun B < 2.1594 (2504.15658 main.tex): UNCONDITIONAL-UNDER-ROOT, moved to §2. 'Assume GRH' (l.43) enters only at l.187–197 via Johnston–Starichkova 2208.01229 Lemma 5.1 (main.tex l.621–638: \|π(x; q, a) − li(x)/φ(q)\| ≤ c_π √x log x for 3 ≤ q ≤ √x, from Grenié–Molteni's Chebotarev for ℚ(ζ_q)/ℚ (JNT 200, 2019 = arXiv 1709.07609, abstract: 'the Riemann hypothesis for the Dedekind zeta functions'; ℚ(ζ_q) and its subfields are abelian, so their ζ are products of Dirichlet L)). Residual: moduli [d₁, d₂] ≤ 2 fall outside Lemma 5.1's q ≥ 3; Schoenfeld's RH bound covers them.

## 4. If Masullo's GRH for all Hecke L (claim of 18:18 UTC) is verified
- Moves to discharged (CONDITIONAL.md Masullo column): D1–D3, D5–D7, D9–D12, D15–D19, D21, E1–E3, E5–E7 (effective Chebotarev, Hooley's Artin conjecture, Odlyzko–Poitou–Serre, Weinberger, Breuillard–Varjú, ring-class A = 1, Λ_K = 0, weight-one Λ_f); D25, D26 partly; stays: D4, D8, D13, D20, D22–D24, D27, E4, E8–E10 (simple zeros, non-CM GL(n ≥ 2), Conjecture H, well-spacing).
