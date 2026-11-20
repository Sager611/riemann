# DRAFT — for the orchestrator, if ROOT (TREE.md l.18) is verified

Exact-substring insertions; nothing here edits repo files. Every anchor below occurs exactly once in TREE.md (rg -F -c, 20 Nov 2026 19:03 UTC). Sources and verdicts: CONDITIONAL.md (rows A*, B*), CONSEQUENCES.md.

## 1. TREE.md node lines: insert "<sentence> " immediately after the anchor
| anchor (exact) | sentence (≤ 25 words) |
|---|---|
| `- N1 Weil/Connes positivity — ` | If ROOT holds: Q(f) ≥ 0 for every f and every Dirichlet χ, > 0 for f ≠ 0 (Weil); window records become quantitative only. |
| `- N2 de Bruijn–Newman Λ — ` | If ROOT holds: Λ = 0 (Rodgers–Tao 1801.05914), Λ_χ = 0 for real primitive χ (Dobner 2005.05142, [check]); our Λ records moot. |
| `- N3 Li/Keiper — ` | If ROOT holds: λ_n > 0 for all n and Re λ_n(χ) > 0 for every primitive χ (Li; Lagarias math/0404394); range records moot. |
| `- N4 Robin/Lagarias — ` | If ROOT holds: Robin for all n ≥ 5041, Lagarias' H_n form, Nicolas for all primorials, our odd Robin beyond N₀; range records moot. |
| `- N5 Nyman–Beurling/Báez-Duarte — ` | If ROOT holds: d_N → 0; liminf d_N² log N ≥ C; L^p density for p < 2; d_N² ~ C/log N needs BCF's hypotheses. |
| `- N6 Jensen-polynomial hyperbolicity — ` | If ROOT holds: every J^{d,n} hyperbolic (Pólya), also for Ξ_t, t ≥ 0, and Ξ_χ of real primitive χ; degree records moot. |
| `- N7 Zero-free regions and zero-density constants — ` | If ROOT holds: no zero of any Dirichlet L with Re s > 1/2; Dirichlet half-planes and densities moot; Hecke/non-abelian ζ_K regions stay live. |
| `- N8 Critical-line proportion — ` | If ROOT holds: all zeros on the line, ≥ 19/27 simple (Bui–Heath-Brown), DPC gaps λ > 3.638 and < 1/2 + ε; simplicity open. |
| `- N9 Moments, Lindelöf, subconvexity — ` | If ROOT holds: Lindelöf and GLH for Dirichlet L, Harper's moments; Simonič's explicit RH bound beats every explicit t^{1/6} bound for t ≥ 10^50.81. |
| `- N10 Explicit ψ and π − li bounds — ` | If ROOT holds: Schoenfeld's bounds, a prime in [x, x + (22/25)√x log x] for x ≥ 4; x^0.856 and Mertens-exponent records moot. |
| `- N11 GRH families and Siegel zeros — ` | If ROOT holds: GRH for every Dirichlet L and abelian ζ_K; LLS's explicit Linnik, non-residue, Littlewood; non-abelian, Artin, Hecke stay open. |
| `- N12 Function-field and Selberg-class analogues — ` | If ROOT holds: nothing beyond degree 1; GL(2), elliptic-curve and Selberg-class GRH stay open; our GL(2) records stand. |
| `- N13 Hilbert–Pólya and spectral approaches — ` | If ROOT holds: no operator follows; QW_λ ≥ 0, but CCM's simple ground state and self-adjointness stay open. |
| `- ROOT (RH itself) — ` | If verified: RH and Dirichlet GRH proved; consequences in Riemann/Deep/PostRH/CONSEQUENCES.md; loop re-aims at GRH beyond Dirichlet (OPEN.md). |


## 2. COROLLARIES.md block (one line per UNCONDITIONAL result; each conditional on ROOT's verification; A4, A20, B3 omitted as already unconditional or superseded, A12 + A14 and A18 + A19 share a line)
- If ROOT holds: Λ = 0 (Rodgers–Tao [1801.05914](https://arxiv.org/abs/1801.05914) Thm, l.129, with RH ⟺ Λ ≤ 0, l.68); our Λ ≤ 0.1074 moot (postrh B1).
- If ROOT holds: Weil's Q(f) ≥ 0 for every test function and every Dirichlet χ, Q(f) > 0 for f ≠ 0 (Weil's criterion; postrh CONSEQUENCES N1).
- If ROOT holds: λ_n > 0 for all n and Re λ_n(χ) > 0 for every primitive χ, all n (Lagarias [math/0404394](https://arxiv.org/abs/math/0404394) l.1096; postrh B15).
- If ROOT holds: our Li-remainder theorems: mean square between N^{3/2}√log N and N^{3/2}√log N·log log N, limsup \|R(n)\|/(n log n)^{1/4} ≥ 0.289 ([page](discoveries/li-remainder.md); postrh A15).
- If ROOT holds: Robin for every n ≥ 5041 and Lagarias' σ(n) ≤ H_n + e^{H_n} log H_n for every n ≥ 1 ([math/0008177](https://arxiv.org/abs/math/0008177) l.92–97, 171–176; postrh B12).
- If ROOT holds: σ(n) < (e^γ/2) n log log n for every odd n > 277093800961598968947975375 ([page](discoveries/odd-robin.md); postrh A16).
- If ROOT holds: d_N → 0 ([math/0202141](https://arxiv.org/abs/math/0202141) l.47) and liminf d_N² log N ≥ 2 + γ − log 4π exactly (Burnol math/0103058 Thm 1.3; postrh B13).
- If ROOT holds: natural-number L^p Nyman–Beurling density for every 1 < p < 2, for ζ, every Dirichlet L and every abelian ζ_K ([page](discoveries/nyman-beurling-lp.md); postrh A17).
- If ROOT holds: D(λ)² log(1/λ) → Σ_ρ m_ρ²/\|ρ\|² (continuous Nyman–Beurling), and the analogue for primitive χ mod q ≥ 3 (DISCOVERIES N5; COROLLARIES nbdir; postrh A12, A14).
- If ROOT holds: every Jensen polynomial J^{d,n} of ξ is hyperbolic (Pólya, [1902.07321](https://arxiv.org/abs/1902.07321) l.277), and of Ξ_χ for primitive real χ (postrh B14, CONSEQUENCES N6).
- If ROOT holds: at least 19/27 of the zeros of ζ are simple (Bui–Heath-Brown [1302.5018](https://arxiv.org/abs/1302.5018) Thm, l.73–75, "Assuming RH"; postrh B2).
- If ROOT holds: a positive proportion of normalized gaps < 1/2 + ε, and λ > 3.63846723 (DPC [2610.28527](https://arxiv.org/abs/2610.28527) l.62, l.76; postrh A18–A19).
- If ROOT holds: log\|ζ(1/2+it)\| ≤ ((log 2)/2 + o(1)) log t/log log t (Chandee–Soundararajan [0908.2008](https://arxiv.org/abs/0908.2008) l.100–103; postrh B6).
- If ROOT holds: \|ζ(1/2+it)\| ≤ exp((0.5 + 6.992/((log t)^{0.252} log log t)) log t/log log t) for t ≥ 10^40 (Simonič [2010.13307](https://arxiv.org/abs/2010.13307) Thm main; postrh B7).
- If ROOT holds: ∫_T^{2T}\|ζ(1/2+it)\|^{2k} dt ≪_k T log^{k²} T for every fixed k ≥ 0 (Harper [1305.4618](https://arxiv.org/abs/1305.4618) l.87; postrh B8).
- If ROOT holds: Schoenfeld's \|ψ(x) − x\| ≤ (√x/8π)log²x (x > 59) and \|π(x) − li(x)\| ≤ (√x/8π)log x (x > 2657) (Büthe [1410.7015](https://arxiv.org/abs/1410.7015) l.58, 453; B9).
- If ROOT holds: a prime in [x, x + (22/25)√x log x] for x ≥ 4 (CMS [1708.04122](https://arxiv.org/abs/1708.04122) l.368), hence between consecutive k-th powers, k ≥ 3 (routine; B10).
- If ROOT holds: M(x) ≪ √x exp((log x)^{1/2}(log log x)^{5/2+ε}) (Soundararajan [0705.0723](https://arxiv.org/abs/0705.0723) Thm 1; Balazard–de Roton [0810.3587](https://arxiv.org/abs/0810.3587); postrh B11).
- If ROOT holds: the least prime p ≡ a mod q is ≤ (φ(q) log q)² for every q > 3 (LLS [1309.3595](https://arxiv.org/abs/1309.3595) l.217–220; postrh B17).
- If ROOT holds: the least quadratic non-residue mod a prime q ≥ 5 is below (log q)² (LLS [1309.3595](https://arxiv.org/abs/1309.3595) l.158–160; postrh B18).
- If ROOT holds: for q ≥ 3000 the primes ≤ (log q + B(q))² generate (ℤ/q)* (LLS [1309.3595](https://arxiv.org/abs/1309.3595) l.135–145: one lies outside each proper subgroup; B19, A24).
- If ROOT holds: 1/\|L(1,χ)\| ≤ (12e^γ/π²)(log log q − log 2 + 1/2 + …) and h(−q) ≥ (π/(12e^γ))√q(…)^{−1} (−q fundamental) for q ≥ 10^10 (LLS l.243–275; postrh B20).
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
- If ROOT holds: finitely many degree-d rationality fields of weight-k CM newforms (Martin [2509.24119](https://arxiv.org/abs/2509.24119) Thm thm:fin, ERH for quadratic characters; postrh A27).

## 3. Pending checks (UNC? in CONDITIONAL.md; insert only after the check)
- Λ_χ = 0 (primitive real χ), Λ_K = 0 (abelian K): Dobner's S♯ contains both (2005.05142 l.194ff); check only that COROLLARIES.md:18's Λ_χ uses Dobner's heat-flow normalization (B16).
- EKN/ARSX exponent-3, 5 (6, 7) lists complete: their class-group computations below the bound must be unconditional (A13).
- 0.6845 simple zeros (CGdL with GGOS): GGOS's GRH must be Dirichlet-only (B3).
- Martin's complete d = 2, 3 lists (2509.24119 Thm thm2): with the EKN check (A27, A13).
- Brun's constant B < 2.1594 (2504.15658): Johnston 2025 Lemma 5.1's GRH must be Dirichlet-only (D26).


## 4. If Masullo's GRH for all Hecke L (claim of 18:18 UTC) is verified
- Moves to discharged (CONDITIONAL.md Masullo column): D1–D3, D5–D7, D9–D12, D15–D19, D21, E1–E3, E5–E7 (effective Chebotarev, Hooley's Artin conjecture, Odlyzko–Poitou–Serre, Weinberger, Breuillard–Varjú, ring-class A = 1, Λ_K = 0, weight-one Λ_f); D25, D26 partly; stays: D4, D8, D13, D20, D22–D24, D27, E4, E8–E10 (simple zeros, non-CM GL(n ≥ 2), Conjecture H, well-spacing).
