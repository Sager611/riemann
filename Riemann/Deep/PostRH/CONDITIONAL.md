# CONDITIONAL — RH- and GRH-conditional results on our files and node lines, after the Lean statements L1–L5

Given (CONSEQUENCES.md): L1–L5 = RH for ζ and GRH for every Dirichlet L (D1), hence GRH for ζ_K of abelian K (D2), GLH for Dirichlet L (D3), no Siegel zeros (D4). Verdicts: **UNC** = the conditional result follows in full; **UNC?** = follows if the bracketed check holds; **MOOT** = an unconditional record or entry of ours is superseded; **STAYS** = a hypothesis remains (named). Quotes are verbatim from the cited line; "our file" quotes are from the entry, TeX quotes from the arXiv e-print (copies in .scratch/postrh/tex/).
Method: rg of the brief's pattern over COROLLARIES.md (56 lines), RECORDS.md + DISCOVERIES.md (7), discoveries/ (125 lines, 37 files), TREE.md l.5–17 (11 lines); 13 scouts classified the windows (.scratch/postrh/s_chunk00–12.md); every UNC below re-checked by me against the quoted line; scout errors are listed in §F.

## A. Our files: UNCONDITIONAL (the hypothesis concerns only ζ or Dirichlet L)
| # | id | result | hypothesis, quoted | why it follows |
|---|---|---|---|---|
| A1 | COROLLARIES.md:63 | elite primes E(x) ≪ x^{3/4}(log x)^{1/2} (arXiv 2611.19610 Thm B) | "using GRH only through Lemma lem:grh (l.337–342: \|Σ_{p≤x} χ(p)\| ≤ c₀√x log(qx) for every non-principal χ mod q)" | Dirichlet χ mod q (D1) |
| A2 | COROLLARIES.md:70 | #{l : C(k,l) ≥ 0} ≪ φ(k)/log k for Alladi–Johnson constants (2610.23660 Thm grh) | "using GRH only through θ(x,χ) ≪ x^{1/2}(log kx)² (l.379–398)" | χ mod k Dirichlet |
| A3 | COROLLARIES.md:71 | PageRank on LPS graphs: Δ_q(1 − c/log q) ≍_{p,I} 1/N_q for large admissible q (Hua 2610.00109 Thm main-results (ii)) | "proves it under GRH for L(s, χ_D), D < 0" | quadratic Dirichlet |
| A4 | COROLLARIES.md:72 | even m ≥ 2: #{x < p/2 : {x^m/p} > 1/2} = p/4 + O_m(√p log log p) (2607.21259 Thm estimate_GRH) | "GRH enters only through Montgomery–Vaughan's Σ_{n≤x} χ(n) ≪ √p log log p" | Dirichlet; already unconditional at full strength in our entry (no gain) |
| A5 | COROLLARIES.md:76 | 2s-th moments of short Legendre sums over primes with Q^{1/2+o(1)} (Munsch–Shparlinski–Sun 2604.23661 Thm l.700–705) | "GRH enters only through Lemma GRH (l.824–831: Σ_{t≤T}(t/m) ≤ T^{1/2}m^{o(1)}" | (·/m) Dirichlet |
| A6 | COROLLARIES.md:79 | long quadratic character sums with the log₃ gain, √X/x (2509.07651 Thm 1.3, f = 1) | "using it only through Darbar–Maiti's lemma (Σ_{\|d\|≤X} χ_d(n) = main + O(X^{1/2+ε}…))" | d ↦ χ_d(n) is a Dirichlet character in d; 2512.24147 Thm 1.1 only as far as its proof is right ("as written its proof covers only x ≤ X^{δ/2}", same entry): PARTIAL |
| A7 | COROLLARIES.md:39(a) | least Elkies prime ≤ (2 log 4q + 4)², N_E(X) = π(X)/2 + O(√X …), under the paper's own curve assumptions (Lippmann–McGown 2311.17231 Thms 1–2) | "GRH enters only through Prop. charsum (l. 148–173); Bebel's Dirichlet half-plane gives the character sums" | Legendre symbols of t² − 4q: Dirichlet |
| A8 | COROLLARIES.md:40(a) | Kim 2412.19975 Thm 1.2 | "its Assumption (l. 79) is exactly a zero-free half-plane Re s > 1 − c₂ for every Dirichlet L-function" | L3 gives every c₂ ≤ 1/2 |
| A9 | COROLLARIES.md:49 | Gauss mean values and the sphere problem with the exponents stated in 2611.10606 Thm thm:cond and 2610.29636 Thm thm:qrh (no stronger 'GRH exponents' claimed) | "assume only that no Dirichlet L-function" (vanishes in the stated half-plane) | L3 |
| A10 | COROLLARIES.md:50 | super-positivity Λ^{(m)}(1/2, χ) ≥ 0, every m ≥ 0, now for every non-principal primitive real χ, conductor q ≥ 3 (was q ≤ 400,000); χ = 1 excluded (ξ is not of this form at 1/2; reviewer) | "Goldfeld–Huang 2017 Thm 1.2 needs only no zero of Λ(s, χ) with 1/2 < σ < 1" | L3; the conductor range becomes MOOT |
| A11 | COROLLARIES.md:54 | near-null spectrum of the windowed Weil form QW_λ | "Near-null spectrum of the windowed Weil form under RH" | L1 |
| A12 | COROLLARIES.md:81 | lim D_χ(λ)² log(1/λ) = Σ_ρ m_ρ²/\|ρ\|²-type limit for every primitive χ mod q ≥ 3 (conductor one is A14) | "GRH for L(s,χ) implies lim_{λ→0} D_χ(λ)² log(1/λ)" | D1; equality with the simple-zero constant STAYS (simple zeros) |
| A13 | COROLLARIES.md:37 + Exp48/REPORT.md l.28–30 | EKN's and ARSX's lists of imaginary quadratic fields with class-group exponent 3, 5 (6, 7) complete | "odd exponents … need a small split prime (ERH via Bach–Sorenson)"; EKN TeX 1803.02056 l.254: "Let K be a quadratic number field … Assume the extended Riemann hypothesis. Then there exists a split prime p ≤ (1.881 log(\|D\|) + 2·0.34 + 5.5)²" | ζ_K, K quadratic: abelian (D2). UNC? [check: EKN/ARSX class-group computations below the bound unconditional] |
| A14 | DISCOVERIES.md:10, discoveries/nyman-beurling-continuous-rate.md:1 | continuous Nyman–Beurling rate D(λ)² log(1/λ) → Σ_ρ m_ρ²/\|ρ\|² | "Under RH alone" | L1 |
| A15 | DISCOVERIES.md:14, discoveries/li-remainder.md l.7, 46, 48 | Li remainder: Σ_{n≤N} R(n)² ≥ (c₁ − o(1)) N^{3/2}(log N)^{1/2}, c₁ = (3/5)^{1/2}(4/5)²/(2π√2) = 0.0557905…, and Σ_{N<n≤2N} R(n)² ≪ N^{3/2}(log N)^{1/2} log log N (the latter with Montgomery's pair-correlation theorem, itself RH-only; li-remainder.md l.6), limsup \|R(n)\|/(n log n)^{1/4} ≥ 0.289, pointwise (0.202 + o(1))√n log n/√(log log n), exponent in [1/4, 1/2] | "Under RH limsup \|R(n)\|/(n log n)^{1/4} ≥ 0.289"; "Pointwise only (0.202 + o(1))·√n log n/√(log log n) is proved, under RH" | L1 |
| A16 | DISCOVERIES.md:8, discoveries/odd-robin.md:3 | σ(n) < (e^γ/2) n log log n for every odd n > 277093800961598968947975375 | "RH ⟺ σ(n) < (e^γ/2) n log log n for every odd n > N" | L1, direction ⟹ |
| A17 | discoveries/nyman-beurling-lp.md:3, 15, 21; COROLLARIES.md:46 | L^p (natural-number) Nyman–Beurling density for every 1 < p < 2, for ζ, every Dirichlet L and every abelian ζ_K | "ζ(s) ≠ 0 on the open half-plane Re s > 1/p"; "L(s, χ) ≠ 0 on Re s > 1/p"; "zero-freeness on Re s > 1/p is the only hypothesis" | L1, L3, D2; p-ranges MOOT |
| A18 | discoveries/small-gap-two-thirds.md:7, small-gap-seven-elevenths.md:5, small-gap-floor.md:3 | positive proportion of normalized gaps < λ for every λ > 1/2 (DPC); 0.6039, 0.6553, 0.6878 superseded | "Under RH: λ > 0.6878 (Soundararajan), 0.6553 (Wu), 0.6039 (Bui–Goldston–Milinovich–Montgomery), 1/2 (Durkan–Pearce-Crump, arXiv 2610.28527, Thm B)"; TeX 2610.28527 l.62: "Assume the Riemann Hypothesis for the Riemann zeta-function" | L1; our 7/11, 2/3 MOOT |
| A19 | DPCAudit/REPORT.md l.13, TREE.md:12 | large gaps λ > 3.63846723 | TeX 2610.28527 l.76: "Assume the Riemann Hypothesis for $L(s,\chi)$. Then $\lambda>3.63846723$." | D1; ζ included: 2610.28527 l.45 "Let $\chi$ have conductor $Q$, with $Q=1$ corresponding to the Riemann zeta-function"; 3.367 MOOT |
| A20 | discoveries/linear-combinations-line.md:7 | ≥ 34.5% of zeros of two-term Dirichlet combinations on the line (Bombieri–Hejhal 1995 Thm C) | "Under GRH alone: Bombieri–Hejhal (Duke 1995, Theorem C), 34.5%" | constituents Dirichlet; superseded by ours 34.707% |
| A21 | discoveries/linnik-numeric.md:9 | polylogarithmic least non-residue | "Polylogarithmic non-residue bounds: GRH (Bach)" | Dirichlet; explicit (log q)² by LLS (B18) |
| A22 | discoveries/littlewood-half-plane.md:5, 7, 17 | Littlewood's 1/\|L(1,χ)\| ≤ (12e^γ/π² + o(1)) log log q | "Explicit 1/log log q bounds existed only under GRH (Lamzouri–Li–Soundararajan 2015, arXiv 1309.3595; Languasco–Trudgian 2011.08348; Lumley 1804.09850" (page l.5) | Dirichlet; explicit for q ≥ 10^10 by LLS (B20) |
| A23 | discoveries/prime-primitive-root.md:12 | ĝ(p) ≤ ((8/5)(2^{ω(p−1)} − 1) log p)² for p ≥ 10⁹ | "Under GRH: ĝ(p) ≤ ((8/5)(2^{ω(p−1)} − 1) log p)² for p ≥ 10⁹ (McGown–Treviño–Trudgian, arXiv 1508.05182)"; TeX l.76, characters mod p l.86 | Dirichlet mod p |
| A24 | discoveries/unit-group-cover.md:5 | every (ℤ/q)* generated by the primes ≤ (log q + B(q))², q ≥ 3000 | "under GRH a sketch covers every q with polylogarithmically many primes"; rigorous via LLS Thm l.135–145 (primes outside any proper subgroup) | Dirichlet mod q |
| A25 | COROLLARIES.md:39(c) | every sufficiently large even N (effective N₀) is a + b with λ(a) = λ(b) = −1 (Mangerel 2412.17199 Thm, l.127) | TeX l.86: "Assuming the Generalised Riemann Hypothesis for Dirichlet $L$-functions (GRH)"; Lemma lem:GRHEp l.368: "Assume that the product $L$-function $\prod_{\chi \pmod{N}, \chi \neq \chi_0} L(s,\chi)$ is zero-free" | Dirichlet mod N |
| A26 | COROLLARIES.md:44(b) | randomized polynomial-time reduction of SSR-SLP to its one-dimensional variant at GRH strength (Gaillard–Jindal 2304.13605 Thm l.360) | TeX l.423: effective Chebotarev (Serre Thm 4) "for the splitting field of $X^2 - a$" | ℚ(√a): abelian (D2) |
| A27 | COROLLARIES.md:44(c) | finitely many degree-d rationality fields of weight-k CM newforms (Martin 2509.24119 Thm thm:fin, l.176) | TeX l.174: "subject to the Extended Riemann Hypothesis (ERH) for quadratic Dirichlet characters" | quadratic Dirichlet; its Thm thm2 (complete d = 2, 3 lists) also uses EKN's lists: UNC? with A13 |


## B. Literature on our node lines: UNCONDITIONAL (TeX hypothesis line quoted)
| # | node | result | source, hypothesis quoted | L-functions |
|---|---|---|---|---|
| B1 | N2 | Λ = 0 | Rodgers–Tao 1801.05914 l.129 "One has $\Lambda \geq 0$." (unconditional); l.68 "The Riemann hypothesis is the equivalent to the assertion $\Lambda \leq 0$" | ζ |
| B2 | N8 | ≥ 19/27 of zeros simple | Bui–Heath-Brown 1302.5018 l.73–75 "Assuming RH we have κ* ≥ 19/27"; abstract l.26–30: CGG had it "under the additional assumption of the Generalised Lindelöf Hypothesis" | ζ (CGG: + GLH for Dirichlet L, D3) |
| B3 | N8 | N_s(T) ≥ (0.6792 + o(1))N(T); with GGOS's GRH input F(x,T) ≥ 3/2 − \|x\| − o(1) on 1 ≤ \|x\| ≤ 3/2 − δ (l.494–497), N_s(T) ≥ (0.6845 + o(1))N(T) | Chirre–Gonçalves–de Laat 1810.08843 Corollary l.276–283: "Assuming RH, we have … 0.6792 … Assuming GRH, we have … 0.6845" | ζ; 0.6845: UNC? [check GGOS's GRH is for Dirichlet L only]; both below B2 |
| B4 | N8 | ≥ c_ε T log T gaps < (1/2 + ε)·mean | Durkan–Pearce-Crump 2610.28527 l.62 "Assume the Riemann Hypothesis for the Riemann zeta-function" | ζ |
| B5 | N8 | λ > 3.63846723 | 2610.28527 l.76 "Assume the Riemann Hypothesis for $L(s,\chi)$" | one Dirichlet L |
| B6 | N9 | log\|ζ(1/2+it)\| ≤ ((log 2)/2 + o(1)) log t/log log t | Chandee–Soundararajan 0908.2008 l.100 "Assume RH." (Thm zetabound, l.103) | ζ |
| B7 | N9 | \|ζ(1/2+it)\| ≤ exp((0.5 + 6.992/((log t)^{0.252} log log t)) log t/log log t), t ≥ 10^40 | Simonič 2010.13307 Thm main l.123 "Assume the Riemann Hypothesis." | ζ |
| B8 | N9 | ∫_T^{2T}\|ζ(1/2+it)\|^{2k} dt ≪_k T log^{k²} T | Harper 1305.4618 l.87 "Assume the Riemann Hypothesis is true, and let $k \geq 0$ be fixed." | ζ |
| B9 | N10 | \|ψ(x) − x\|, \|θ(x) − x\| ≤ (√x/8π) log² x (x > 59, x > 599), \|π(x) − li(x)\| ≤ (√x/8π) log x (x > 2657) | Büthe 1410.7015 l.58–60: Schoenfeld's bound "which is implied by the Riemann hypothesis"; l.453–454 table; original Schoenfeld 1976 Math. Comp. 30 (not arXiv; ψ from x ≥ 73.2 there); the smaller thresholds also use the finite RH verification recorded by Büthe (reviewer) | ζ |
| B10 | N10 | a prime in [x, x + (22/25)√x log x], x ≥ 4 | Carneiro–Milinovich–Soundararajan 1708.04122 l.368 "Assume the Riemann hypothesis. Then, for $x\ge 4$, there is always a prime number in the interval" | ζ |
| B11 | N10 | M(x) ≪ √x exp((log x)^{1/2}(log log x)^{14}) | Soundararajan 0705.0723 l.62 "Theorem 1 Assume RH." (l.64); Balazard–de Roton 0810.3587: (log log x)^{5/2+ε} | ζ |
| B12 | N4 | Robin for every n ≥ 5041; Lagarias' H_n form for every n ≥ 1 | Lagarias math/0008177 l.173–176 "Robin's criterion states that the Riemann hypothesis is true if and only if σ(n) < e^γ n log log n for all n ≥ 5041"; l.92–97 Problem E; Thm th1 "Problem E is equivalent to the Riemann hypothesis" (l.118) | ζ |
| B13 | N5 | 1 ∈ closure of B_nat in L²(0,1), d_N → 0 | Báez-Duarte math/0202141 l.47 "Riemann hypothesis is equivalent to the stronger statement that $\chi\in\overline{\Bnat}$" | ζ |
| B14 | N6 | every J^{d,n}_γ hyperbolic | GORZ 1902.07321 l.277 "he showed that the RH is equivalent to the hyperbolicity of all Jensen polynomials" | ζ |
| B15 | N3 | λ_n > 0 (ζ) and Re λ_n(χ) > 0 for primitive χ, Li coefficients of the completed ξ(s, χ) (not raw imprimitive L, whose Euler factors add zeros on Re s = 0; reviewer); R(n) = O(√n log n) in Lagarias' normalization | Lagarias math/0404394 l.1096 "We can now state a general version of Li's criterion." (automorphic; GL(1) = Dirichlet) | ζ, Dirichlet |
| B16 | N2 | Λ_χ = 0 for primitive real χ; Λ_K = 0 for abelian K | Dobner 2005.05142 l.118–120 "$\Lambda_F \geq 0$ for every $F \in\mathcal{S}^\sharp$." + D1/D2 | S♯ (l.194ff: Dirichlet series, finite-order continuation, functional equation with Γ-factors, Re μ_i ≥ 0) contains L(s, χ), χ real primitive (root number 1), and every ζ_K; Λ_F via ξ^F_t, l.113. UNC? [check only that COROLLARIES.md:18's Λ_χ uses the same heat-flow normalization] |
| B17 | N11 | P(a,q) ≤ (φ(q) log q)², q > 3 | LLS 1309.3595 l.217–220 "Assume GRH. Let $q>3$ … $P(a,q)\leq (\phi(q)\log q)^2$" ; GRH there is for Dirichlet L: l.576 "Suppose that GRH holds for $L(s,\chi)$" | Dirichlet |
| B18 | N11 | least quadratic non-residue mod prime q ≥ 5 below (log q)² | LLS l.158–160 "Assume GRH. If $q\ge 5$ is prime, the least quadratic non-residue $\pmod q$ lies below $(\log q)^2$" | Dirichlet |
| B19 | N11 | least prime ∤ q outside a proper subgroup H ≤ (log q + B(q))², q ≥ 3000 | LLS l.135 "Assume GRH. Let $q\geq 3000$ be an integer, and let $H$ be a proper subgroup" (bound l.140) | Dirichlet |
| B20 | N11 | 1/\|L(1,χ)\| ≤ (12e^γ/π²)(log log q − log 2 + 1/2 + 1/log log q + 14 log log q/log q), primitive χ, q ≥ 10^10; h(−q) ≥ (π/(12e^γ))√q(…)^{−1} for −q a fundamental discriminant | LLS l.243–248, l.272–275 | Dirichlet |
| B21 | N11 | ĝ(p) < √p − 2 for p > 2791; g(p) < √p − 2 for p > 409 (Grosswald) | McGown–Treviño–Trudgian 1508.05182 l.45 "Assume GRH."; l.86 "Let $\chi(n)$ denote a non-principal Dirichlet character modulo $p$" | Dirichlet mod p |
| B22 | N11 | every odd composite n has a strong witness ≤ 2 log² n | Bach 1990, Math. Comp. 55 (not on arXiv; our COROLLARIES.md:36 "under GRH 2 log² n (Bach 1990)"); arXiv-checkable form: B19 with H ⊇ strong liars | Dirichlet mod n |


## C. MOOT: our unconditional records and entries superseded (by A/B; CONSEQUENCES.md has the node view)
| id | record or entry | superseded by |
|---|---|---|
| RECORDS.md:7, TREE N1 | window positivity of Weil's form (as RH evidence) | Q ≥ 0 for all f (N1); explicit ε_c not implied, stays as a quantitative record |
| RECORDS.md:9–10 | Robin/Nicolas ranges, Robin for v₂(n) ≤ 72 | B12 |
| RECORDS.md:11 | Mertens exponent 0.9382 (exponent only) | B11; the explicit all-x inequality STAYS (no explicit RH bound) |
| RECORDS.md:12 | Λ ≤ 0.107410439495 | B1 |
| RECORDS.md:13 | explicit Linnik (q log q)^{7.8956} | B17 |
| RECORDS.md:15 | prime between consecutive k-th powers (k ≥ 3) beyond e^{626.82} etc. | B10 for n ≥ 2 (routine) |
| RECORDS.md:16 | g*(p) ≤ (1.956 log p)^{26.85} | PARTIAL: B21/A23 below p, Shoup–Martin ≪ (ω log ω)⁴ log² p non-explicit; uniform explicit polylog STAYS ours |
| RECORDS.md:17 | (x, x + x^0.856] | B10 (shorter for x ≥ 4; x < 4 by hand) |
| RECORDS.md:18 | \|L(1,χ)\| ≥ 0.11/max(3, log log q) | B20 for q ≥ 10^10 (≈ 0.46/log log q); q < 10^10 STAYS ours |
| RECORDS.md:20 | explicit log-free zero density for ζ | L1 (N(σ,T) = 0, σ > 1/2) |
| RECORDS.md:21–22 | \|ζ(1/2+it)\| ≤ 0.568 t^{1/6} log t; 54.6 t^{27/164} | B7 beats 0.568, 0.4821 and Fodge's bound for t ≥ 10^46.92, 10^48.09, 10^50.81 (logs/n9cross.log); below, the records stand |
| RECORDS.md:25–26, COROLLARIES.md:7, :22 | Jensen ranges (ζ; real χ ≤ 400,000) | B14; D1 |
| RECORDS.md:28 | large gaps Λ > 3.367 | A19/B5 |
| RECORDS.md:30 | small gaps λ > 7/11 (and 2/3) | A18/B4 |
| RECORDS.md:37, COROLLARIES.md:8, :20 | Li ranges (ζ; Dirichlet n ≤ 26.12 t0(q)²) | B15 |
| COROLLARIES.md:9 | liminf d_N² log N ≥ C − 6·10⁻¹² | C exactly (N5) |
| COROLLARIES.md:10, nyman-beurling-lp.md | L^{8/7}; p ≤ 1.14497 | A17 |
| COROLLARIES.md:3 | Linnik L ≤ 7/3 | B17 (P(a,q) ≤ (φ(q) log q)²) |
| COROLLARIES.md:6 | Lindelöf on σ > 7/8 | D3 |
| COROLLARIES.md:12, :15 | effective Siegel–Walfisz, effective Duke (via the exceptional-zero gap) | input only: D4 removes the exceptional zero (ψ(x;q,a) = x/φ(q) + O(√x log² x) under GRH); Duke's rate γ unchanged (needs subconvexity, not GRH-Dirichlet) |
| COROLLARIES.md:18 | Λ_χ ≤ 0.2788 | B16 (UNC?) |
| COROLLARIES.md:28 (abelian part) | Λ_K ≤ 225/512 | B16 for abelian K; non-abelian STAYS |
| COROLLARIES.md:35, :36 | \|L′/L(1,χ)\| ≤ 7.8956 log log q + 10.53; Miller bases ≤ (log n)^{7.8956} | "GRH gives 2 log log q + O(1)" (entry); B22 |
| COROLLARIES.md:50 (range) | super-positivity for conductor ≤ 400,000 | A10 |
| TREE N7, COROLLARIES.md:14 (Dirichlet part), :30 | PR 144's 0.873332794, θ₁₄₁, B₂ for Dirichlet L; family densities | L3 |
| TREE N8 | 67.3639% simple on the line; ceilings 31/46, 153/227 | B2 (19/27) |
| TREE N9 | μ(1/2) ≤ 13/84 | D3 |
| DISCOVERIES.md:15–16 | Lamzouri-criterion floor, hard-core crossing (as small-gap records) | A18 (true statements about methods, no longer frontier) |


## D. Our files: STAYS CONDITIONAL (remaining hypothesis named)
| # | id | result | remaining hypothesis | Masullo (if his Hecke GRH holds) |
|---|---|---|---|---|
| D1 | COROLLARIES.md:42, :65, :66 | Serre open image, effective isogeny, isogeny primes with GRH-strength bounds (Mayle–Wang etc.) | GRH for Artin/Dedekind L of division fields K(E[ℓ]) (non-abelian) | discharged: effective Chebotarev in the Galois division fields needs GRH for their ζ_L (Lagarias–Odlyzko) [check: no Artin-holomorphy step] |
| D2 | COROLLARIES.md:44(a) | Gelander–Slutsky lattice bound with Landau's 13√x(log D_k + d log x) | GRH for ζ_k, general number field k | discharged: GRH for ζ_k, every k |
| D3 | COROLLARIES.md:40(e) | Fuselier et al. (2306.03051) | GRH through Wesolowski's Thm 6.4 (class-group/quaternion L, not Dirichlet) | discharged [check]: Wesolowski's GRH is for Hecke L of imaginary quadratic orders (class-group expansion), not modular L |
| D4 | COROLLARIES.md:51 | second-order term of the BCF ramp | simple zeros + Σ_{\|γ\|≤T}\|ζ′(ρ)\|⁻² ≪ T^{3/2−δ} (BCF 1211.5191 l.48) | not discharged: simple zeros and a moment bound, not GRH |
| D5 | COROLLARIES.md:61; TREE N2 | Λ_f for weight-one newforms; Λ_K for non-abelian K | GRH for 2-dim Artin L (Deligne–Serre) and non-abelian ζ_K | discharged: Λ_K ≤ 0 from GRH for ζ_K (Dobner: = 0); weight-one L(s, ρ_f): Brauer induction + Deligne–Serre holomorphy + functional equation |
| D6 | COROLLARIES.md:67 | generalized Euclidean pure orders with GRH bounds (Weinberger route) | GRH for Kummer fields Q(ζ_n, ε^{1/n}) | discharged: Weinberger needs GRH for ζ of the Kummer fields |
| D7 | COROLLARIES.md:68 | Koblitz–Zywina least prime divisors | GRH for division fields of E | discharged: GRH for ζ of the division fields [check: no Artin-holomorphy step] |
| D8 | COROLLARIES.md:69 | divisor moments of F(a_f(p)) | GRH for ζ of Galois-representation fields + symmetric-power zero-free regions (Hypothesis A) | not discharged: symmetric-power L of non-CM f (GL(n), n ≥ 3); its Dedekind part is discharged |
| D9 | COROLLARIES.md:73 | KMR Selmer distributions with ERH errors | ERH in effective Chebotarev ("ERH enters only through effective Chebotarev (l.485–490)"); fields not abelian over Q as far as stated [check] | discharged: effective Chebotarev needs only GRH for ζ_L |
| D10 | COROLLARIES.md:74 | generalised Euclid–Mullin with GRH thresholds | Chebotarev in Kummer fields | discharged: ζ of the Kummer fields |
| D11 | COROLLARIES.md:75 | largest prime factor of a_f(p) + a_g(p) at GRH strength | Chebotarev in mod-h Galois-representation fields | discharged: ζ of the mod-h Galois-representation fields [check: no Artin-holomorphy step] |
| D12 | COROLLARIES.md:78 | effective Brauer–Siegel with log log for Artin L at GRH strength | GRH for ζ_K over general K/k | discharged: GRH for ζ_K over every K/k |
| D13 | COROLLARIES.md:80, :81 (equality part); DISCOVERIES.md:10 (constant C) | D(λ)² log(1/λ) → C = 2 + γ − log 4π | all zeros simple ("if and only if RH holds and every nontrivial zero of ζ is si[mple]") | not discharged: simplicity of the zeros of ζ |
| D15 | discoveries/chebotarev-numeric.md:9; RECORDS.md:14 | N𝔭 ≤ (4 log Δ_L + 2.5 n_L + 5)² (Bach–Sorenson) | GRH for ζ_L of non-abelian L (abelian L/Q: UNC) | discharged: Bach–Sorenson assume GRH for ζ_L |
| D16 | discoveries/cm-least-primes.md:14 | least prime d_K^{1+ε} for every field | GRH for ζ of the Galois closure (abelian: UNC) | discharged: ζ of the Galois closure |
| D17 | discoveries/discriminant-bounds-strip.md:9, 11, 25; RECORDS.md:19 | Odlyzko–Poitou–Serre limits 44.7632, 215.3325 | GRH for ζ_K, K non-abelian (abelian K: UNC) | discharged: GRH for every ζ_K |
| D18 | discoveries/least-primes-forms.md:11; DISCOVERIES.md:11 | ring-class least primes with A = 1; ring-class densities | GRH for ring-class Hecke L (dihedral over Q) | discharged: ring-class characters are Hecke characters of K; ring-class densities become moot (N = 0 for σ > 1/2) |
| D19 | discoveries/random-polynomials-without-rh.md:3, 9 | Breuillard–Varjú irreducibility probability | GRH for ζ of Q[x]/(P), Galois group S_n | discharged: ζ of Q[x]/(P) |
| D20 | discoveries/sato-tate-classes.md:20, 23, 41 | Sato–Tate in ideal classes at the GRH threshold; Blomer–Brumley's per-D form | GRH for Hecke/Rankin–Selberg L over imaginary quadratic fields | not discharged: twists L(f ⊗ ψ) of a non-CM newform are GL(2) × GL(1)_K, automorphic [check the page's exact L-functions] |
| D21 | discoveries/class-number-10000.md:5 | HJKMP counts for odd h < 10⁶ "under GRH" | GRH for class-group verification (Hecke L of K, dihedral) [check]; GGLM's h ≤ 2000 MOOT (unconditional, 2610.21511) | discharged: class-group certification needs GRH for Hecke L of K (Bach) |
| D22 | TREE.md:5 | Mori's sharpness of the windowed Weil eigenvalue | zero separation | not discharged: zero separation |
| D23 | TREE.md:12 | κ* → 7/8 | RH + Conjecture H (Σ_γ ζ′(ρ)²B(ρ)² asymptotic) | not discharged: Conjecture H |
| D24 | TREE.md:12, linear-combinations-line.md:7 | 100% of zeros of Dirichlet combinations on the line | Bombieri–Hejhal's well-spacing hypothesis | not discharged: well-spacing |
| D25 | TREE.md:15 | Omega results for infinite-order Hecke characters; 2505.11226 (m = 1), 2502.11987 | ERH for Grössencharacters; non-abelian or modular L | partly: Grössencharacter ERH discharged [check: Masullo covers infinite-order characters]; Artin parts zero-free; modular parts not |
| D26 | TREE.md:15 | GRH savings over ≥ X^{1/4} family members (NOGO); Brun's constant B < 2.1594 (2504.15658 l.43 "Assume GRH.") | families not identified in our files; Brun: UNC? [check: its GRH input, Johnston 2025 Lemma 5.1 (l.197), is Dirichlet-only] | discharged only if the family is Dedekind/Hecke/Artin [check]; Brun is Dirichlet already (UNC?) |
| D27 | TREE.md:16, RECORDS.md:29, :35, DISCOVERIES.md:17 | degree-2 Selberg class, GL(2) gaps and critical-line proportions | GRH for GL(2) L-functions | not discharged: non-CM GL(2); CM and weight-one subcases discharged |


Masullo flags (orchestrator 64, 19:20 UTC): A. M. Masullo's unverified claim (18:18 UTC, 940 pp.) of GRH for the Hecke L-functions of every number field gives GRH for every ζ_K and, by Brauer induction (integer powers of Hecke L, all holomorphic and nonzero in Re s > 1/2), no zero of any Artin L in Re s > 1/2; not Artin holomorphy, not GL(n) automorphic (n ≥ 2, non-CM), not non-CM elliptic curves. Tally: D 16 discharged (5 with a [check]), D26 only if its family is Hecke-type, D25 partly, 8 not (D4, D8, D13, D20, D22, D23, D24, D27); E 6 discharged (E1–E3, E5–E7), 4 not (E4, E8–E10). Results that also need Artin holomorphy (e.g. Murty–Murty–Saradha's effective Chebotarev) stay; the [check] flags name where that matters.

## E. Traps checked explicitly (each STAYS, with the hypothesis line)
| # | trap | hypothesis, quoted | why L1–L5 do not give it | Masullo (if his Hecke GRH holds) |
|---|---|---|---|---|
| E1 | Artin's primitive-root conjecture (Hooley 1967) | Moree math/0412262 l.555–556: "Hooley proved in 1967 that if the Riemann Hypothesis holds for the number fields $\Q(\zeta_k,g^{1/k})$ with $k$ squarefree, then" (HRH, l.575–576) | Gal(Q(ζ_k, g^{1/k})/Q) ⊂ ℤ/k ⋊ (ℤ/k)*, non-abelian for k ≥ 3 (k = 3: S₃); ζ of these fields contains 2-dim Artin factors | discharged: HRH is GRH for ζ of Q(ζ_k, g^{1/k}), Dedekind |
| E2 | effective Chebotarev for non-abelian extensions | discoveries/chebotarev-numeric.md:9: "Under GRH, N𝔭 ≤ (4 log Δ_L + 2.5 n_L + 5)²" (Bach–Sorenson) | needs GRH for ζ_L of the Galois closure; only abelian L/Q reduce to Dirichlet L (residue classes mod the conductor, B17) | discharged: GRH for ζ_L of the Galois closure |
| E3 | Hecke L of non-abelian (ring-class, Hilbert-class) fields | discoveries/least-primes-forms.md:11: "GRH improves the least prime represented in every ring class to constant A = 1"; Bach's class-group generators ≤ 12 log²\|d\| (Bach 1990) likewise need GRH for the class-group L-functions | ring class fields are generalized dihedral over Q; their characters give irreducible 2-dim Artin L. Not a trap: class NUMBER computation, which needs "the assumption of GRH for $L(s,\chi_{-q})$" (LLS 1309.3595 l.270): Dirichlet, UNC | discharged: Hecke L of K; Bach's generators follow |
| E4 | elliptic-curve L | Heath-Brown math/0305114 l.63–65: "Assume that the $L$-functions of all the curves $E_{r,s}$ satisfy the Riemann Hypothesis. Then … $\leq 2+o(1)$" | GL(2) L(E, s) not given | not discharged: L(E, s) of non-CM curves (GL(2)); CM curves only |
| E5 | Hecke L of an abelian field that are not Dirichlet | TREE N11: "zero-free Re s > 11/12 for finite-order Hecke L-functions of Q(√−3)" | Galois-invariant ψ = χ∘N: L(ψ,s) = L(χ,s)L(χχ₋₃,s), UNC; non-invariant ψ: Ind ψ irreducible (weight-one CM form), STAYS, though Q(√−3) is abelian | discharged: Hecke L of ℚ(√−3) |
| E6 | weight-one forms / 2-dim Artin | COROLLARIES.md:61: "L(s, f) = L(s, ρ_f) (Deligne–Serre) is, by Brauer induction," | Brauer induction gives quotients of Hecke L of non-abelian fields | discharged: zero-free on Re s > 1/2 by Brauer induction; entire by Deligne–Serre; functional equation |
| E7 | Dedekind ζ of S_n fields (random polynomials) | discoveries/random-polynomials-without-rh.md:3 (Breuillard–Varjú's "Dedekind-zeta RH theorem") | Q[x]/(P) has Galois group S_n | discharged: Dedekind ζ of S_n fields |
| E8 | hidden simplicity | Bettin–Conrey–Farmer 1211.5191 l.249: "Assume the Riemann hypothesis and that all the zeros of $\zeta(s)$ are simple." | simplicity is open (OPEN.md #2) | not discharged: simplicity |
| E9 | Legendre's conjecture | (n+1)² − n² = 2n + 1 ≍ 2√x, but B10 needs (22/25)√x log x | RH-strength intervals are too long by log x | not discharged: not a GRH question |
| E10 | Montgomery's pair-correlation conjecture | Montgomery 1973 proves F(α) only for \|α\| ≤ 1 under RH | the conjecture (F = 1 for \|α\| ≥ 1) is beyond RH | not discharged: beyond any GRH |

## F. Falsification log (scout verdicts corrected)
- COROLLARIES.md:9, :20, :22 were "UNC" (finite-range records): MOOT, the full statements follow (C exactly; Re λ_n(χ) > 0 for all n; J^{d,n} for all d).
- TREE.md:15 Littlewood was "MOOT, superseded by the unconditional 0.11/max(3, log log q)": reversed; ours is superseded for q ≥ 10^10 by LLS l.243–248 (≈ 0.4618/(log log q + …)).
- discoveries/small-gap-*.md (chunk10) and li-remainder.md l.12, 46, 48 (chunk07) were "STAYS?(check), unclear": UNC, the pages say "Under RH" (ζ only), DPC TeX l.62.
- COROLLARIES.md:37, unit-group-cover.md:5 were "STAYS?(check)": UNC? (EKN TeX l.254, quadratic field) and UNC (LLS l.135, Dirichlet mod q).
- COROLLARIES.md:72 "UNC": true but no gain (our entry is already full strength).
- Our own trap: DPC Thm C assumes RH "for L(s,χ)", one Dirichlet L, not ζ alone; covered by D1 either way.

