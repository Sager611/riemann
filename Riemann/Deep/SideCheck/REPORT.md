# SideCheck: review and novelty check of two side records (agent sidecheck, 2026-11-05)

(1) VERDICT: CORRECT (ξ(K₆) ≥ 33 holds, and the definition matches every literature form); NOVELTY: NEW (nothing above 25 is documented in any source we could open; one unopened source: Xu–Shao–Su–Li 2009, Table 2)

- **Own verification.** I wrote my own checker, `code/xi_check.py`, using a different method from `Riemann/Deep/Xi6/verify.py` (which brute-forces every 6-set). Mine computes, for each pair of colours, the exact clique number of the union graph (Bron–Kerbosch search on bitmasks) and counts its K₆ and K₅ subgraphs. On `K33_Z11x3.txt`: the matrix is symmetric with entries 0/1/2; the colour classes have 198, 154 and 176 edges. The unions {0,1}, {0,2} and {1,2} all have clique number 5, with no K₆ (they have 1155, 2255 and 748 K₅'s, all multiples of 11, matching the order-11 symmetry). So every K₆ sees all three colours.
- **Controls** (`code/xi_controls.py`):
  - The W/X tables in Xi6/REPORT §2 rebuild the file exactly, and the colouring is invariant under the shift x ↦ x+1.
  - The cubic-residue K₃₁ colouring of arXiv 2611.03878 fails for k = 6 (clique numbers 6, 6, 6) and passes for k = 7.
  - 20 of 20 random colourings of K₃₃ are rejected.
  - Only 44 of the 1056 single-edge recolourings stay valid.
  - Also valid: K55_circulant_K7 for k = 7 (clique numbers 6, 6, 6), and K26 and K28 for k = 6.
- **Definitions** (checked in TeX/PDF). In every form, colours are labels: any relabelling is equivalent, and the diagonal case has no ordering issue. The condition is "every K₆ uses all three colours", not "no monochromatic K₆" (that would be the 3-colour Ramsey number R₃(6), a different and much larger number). The forms:
  - Almási–Axenovich, arXiv 2512.08053 (TeX l.166): ξ(H) is the largest n such that some 3-colouring of K_n gives every copy of H all three colours. Our statement is exactly this.
  - 2611.03878 (l.77) uses the same definition.
  - Set-colouring Ramsey number: R_{3,2}(H) = ξ(H)+1 (2512.08053, l.244–245); R(n; r, s) in Conlon–Fox–He–Mubayi–Suk–Verstraëte (arXiv 2206.11371). Translation: give each edge the list of the two colours it does not have; then "a K₆ whose lists share a colour" means "a K₆ missing that colour".
  - Chung–Liu 1978: r³₂(K₆,K₆,K₆), as stated in Budden–DeJonge (arXiv 2610.41971, l.107).
  - Harborth–Möller 1999: R_{2,3}(K₆).
  - Budden–Stender–Zhang 2017: R³₂(K₆), with R³₂(K₃) = 5 and R³₂(K₄) = 10, i.e. ξ + 1.
  - Munemasa–Shinohara 2019 (J. Indones. Math. Soc. 25:146): complementary Ramsey number R̄(6,6,6) = least n such that every 3-colouring has a colour class of independence number ≥ 6, i.e. ξ + 1.
  - Xu–Shao–Su–Li 2009 (Graphs Combin. 25:863): multigraph Ramsey number f⁽²⁾(6,6,6), per Munemasa–Shinohara.
  - Erdős–Gyárfás generalized Ramsey number: f(33,6,3) = 3.
  - So the result can be stated as ξ(K₆) ≥ 33, or equivalently R(6;3,2) = R³₂(K₆) = R̄(6,6,6) ≥ 34.
- **Prior bounds.**
  - The grid colouring gives (k−1)² = 25. This is Chung–Liu 1978 Thm 2.6, i.e. R̄(6,6,6) ≥ 26 (Munemasa–Shinohara §7 credits Chung–Liu Thm 2.6 with R̄(5,5,5) ≥ 17, which is (k−1)²+1 for k = 5).
  - 2611.03878 claims ξ(K₆) ≥ 31 with a wrong colouring (our REFUTATIONS); it names 25 as the previous bound.
  - Almási–Axenovich (Dec 2025) give exact or bracketed values only for K₃, K₄, K₅ (l.197).
  - Budden–DeJonge (30 Oct 2026) list the known r³₂ values; none is for K₆.
  - Budden–Stender–Zhang 2017 (Zenodo 10636187, PDF read) has nothing for K₆ ("little is known for n ≥ 5").
  - DS1.18 (24 Apr 2026, §9.1) excludes weakened and set-colouring Ramsey numbers.
- **Gallai papers cannot preempt it.** Gallai colourings are colourings with no three-coloured triangle. Wang–Zhang 2026, "Weakened Gallai–Ramsey numbers for K6", is paywalled. My exact recursion over Gallai partitions (`code/gallai_xi.c`) shows the largest Gallai colouring in which every K₆ sees three colours has 11 vertices. The recursion reproduces 8 published values of Budden–DeJonge (for example their gr³₂(K₄,K₄,K₄) = 5 and gr³₂(K₄,K₅,K₆) = 8). So no Gallai-type construction reaches 33.
- **Not opened** (paywalled or not public):
  - Xu–Shao–Su–Li 2009, Table 2: small lower bounds for f⁽²⁾ = R̄. It contains R̄(5,5,5) ≥ 20 (quoted by Munemasa–Shinohara and by arXiv 1605.05074). Bustamante–Stein (1510.05190) say it "studies r_{3,2}(K_t) for small t". Its K₆ entry, if any, is unknown. Springer's page shows a bot challenge.
  - Budden, *Weakened Ramsey Theory* (CRC, June 2026), §3.1 (only the chapter abstract was read).
  - Xu–Liang–Luo, *Ramsey Theory* (De Gruyter 2018), ch. 16 "Set-Coloring Ramsey Theory".
  - Song–Mao 2025 (SSRN 5380593).
  - The KIT theses of Ortlieb (2017), Lee (2022) and Almási (2023).

  The papers that cite or survey these sources since then (2512.08053, 2611.03878, 2610.41971) do not report anything above 25.
- **Small fixes for Xi6/REPORT.md** (do not affect the record):
  - "previous proved bound 25" should read: "best documented lower bound 25 (grid); the Xu–Shao–Su–Li 2009 table was not seen".
  - §7 compares ξ(K₇) ≥ 55 with "31 published" (Budden–Stender–Zhang). The grid already gives ξ(K₇) ≥ 36 (and ξ(K₈) ≥ 49 > 43), so the baseline for K₇ is 36.
  - The upper bound 103 was not re-checked here.

(2) VERDICT: CORRECT (C_41b ≤ 19/10, unconditional, computer-assisted; the proof and the certificate checked independently); NOVELTY: NEW (no proved bound below 2.2195 anywhere; only non-rigorous or conditional claims, listed below)

- **Definition.** The 41b page (fetched 5 Nov 2026): a connected rigid planar shape (as in 41a) that passes a unit-width Z-shaped corridor with the corners sufficiently far apart. PROOF.md handles this formulation, with corners ≥ 10·diam(S) apart. It also handles Romik's "right and left turn from the same placement" (arXiv 1606.08111 TeX §1.3; Kallus–Romik 1706.06630 §5). Connectedness is used for the bound on the x-extent; 41a and Kallus–Romik (l.81) both require it.
- **Proof review.** I re-derived every step by hand:
  - Sofa frame: the second corner equals Rot_π(L₀)+c, so its hallway angle is π−θ. The left-turn hallway M(L₀) equals Rot_{−π/2}(L₀)+(0,1), angle −π/2−θ.
  - Strip bounds: in the vertical segment |cos θ| < 1/T; at the end |sin θ| < 1/T. For times on the midline, D ≥ 10·diam and diam > 1 put the whole sofa inside the vertical strip.
  - Every case of Step 3: k constant with k ≤ −1, 0, 1, ≥ 2 (with m ≥ 0 or m ≤ −1), and k changing along the motion. Each yields one of the pairs {0,2}, {0,3}, {1,2}, {1,3}. The same-placement formulation yields {0,2}, {0,3}, {2,3}, {1,3}, or Q₃ alone.
  - Step 4: the x-extent bound, (1+|cos|)/|sin| + (1+|sin|)/|cos| = 29/6 ≤ 155/32, and the clamping (both directions, null-set changes only).
  - Step 5: the expanded arms contain C(ψ,u) for every u in the box; inclusion–exclusion over {Pt, Qt, Pt∩Qt}.

  No gap found. One harmless overstatement: uniqueness of k needs only ε < π/2, not ε < π/4.
- **Certificate.**
  - `verify.py` (sha256 0f6890c0…e6e0, alarm 1800 s, single thread) printed CERTIFIED in 240 s, with the same per-case leaf counts and maxima as OptConst/REPORT.md.
  - Independent re-check (`code/sofa_prep.py` + `code/sofa_recheck.c`): my own tree traversal, using the split rule stated in PROOF.md, and my own area routine (vertical-slab integration, exact for piecewise-linear cross-sections; floating point). It consumes every tree bit and gives the same leaf counts: 2800, 513503, 3006, 2626, 2100, 2787, 831, 833. The per-case maxima agree with verify.py to 12 digits (largest 1.899999822499 in case {0,2}). The smallest margin to 1.9 is 1.8·10⁻⁷, far above rounding error.
- **Conventions** (`code/sofa_convention.py`). An explicit ambidextrous lens (area 1.2284) was moved through the Z-corridor. The sofa-frame hallway equals C(ψ,u) at all 366,000 sampled point/time pairs, with ψ = −θ at corner 1 and π−θ at corner 2. The hallway angles fill [0°, 90°] and [180°, 270°], i.e. windows 0 and 2. The lens fits the certificate's corridors at φ and π+φ.
- **Trying to break it** (`code/sofa_attack.c`). A multistart local maximisation of the exact relaxation found suprema of 1.82976 for {0,2} (the agent reported 1.8298), 1.45 for the other pairs and 1.3812 for a single window. All are below 1.9, consistent with the certificate.
- **Prior and competing claims.**
  - The 41b page lists 2√2 and 2.2195 (Baek 2024, which the 41a page marks as announced and unverified). Unconditional published: 2.37 (Kallus–Romik 2018, via 41a).
  - Kallus–Romik's TeX (end of §5) only proposes the ambidextrous extension.
  - Baek's arXiv 2411.19826 and 2406.10725 contain nothing on the ambidextrous problem (the TeX mentions it only in a bibliography abstract and in a Mathematica file).
  - O. Mohammed, Zenodo 23667275 (26 Oct 2026, with a GitHub repo) calls U₆ = 1.6906 … U₃₀ = 1.6505 "rigorous upper bounds". His own `src/upper_bound.py` says these are heuristic SciPy maxima, i.e. lower estimates ("to claim a proven ceiling we would need branch-and-bound"). His `src/interval_bound.py` only certifies small balls around his optimum. His model also assumes each turn sweeps [0, π/2], which is unproved for this problem. So these are not upper bounds.
  - V. Bonfioli, Zenodo 22734294 (13 Sep 2026): a conditional local stability bound near Romik's sofa, under a curvature hypothesis (RC) and further conditions. Not a global bound.
- **Caveats.** Computer-assisted, unrefereed, written by AI agents; no human has reviewed it yet.

## Sources searched (all 2026-11-05 UTC)
- **arXiv TeX read:**
  - 2512.08053, 2611.03878, 2206.11371 (set-colouring definitions, prior bounds).
  - 2610.41971 (Budden–DeJonge, 30 Oct 2026).
  - 1406.2050 (Munemasa–Shinohara).
  - 1605.05074, 1510.05190 (how they cite Xu–Shao–Su–Li).
  - 1706.06630 (Kallus–Romik), 1606.08111 (Romik).
  - 2411.19826 and 2406.10725 (Baek), 2608.11206.
- **arXiv API** (all categories, all dates): "weakened Ramsey", "set-coloring/colouring Ramsey", "chromatic Ramsey number", "complementary Ramsey", "multigraph Ramsey", "(t-1)-chromatic", Erdős–Gyárfás, generalized Ramsey f(n,p,q), polychromatic, balanced colouring, Gallai-Ramsey K_6; sofa (129 hits), ambidextrous (40), au:Baek_Jineon (18); abstracts of 2610.41971, 2610.24484, 2610.26061, 2610.07546, 2610.29584, 2609.35079, 2512.16599, 2505.20652, 2602.03865.
- **arXiv monthly listings**, show=2000, keyword scan of all titles: math.CO 2026-10 (1460 entries) and 2026-11 (274); math.MG 2026-10 (256) and 2026-11 (38); math.OC 2026-10 (1124) and 2026-11 (193).
- **Surveys and pages:** Radziszowski DS1.18 PDF (24 Apr 2026); the teorth 41a and 41b pages.
- **GitHub teorth/optimizationproblems** (read only): searches "sofa", "ambidextrous", "41b" (#25, #114, #132, #164, #223); the 100 latest PRs and issues (up to #308, 4 Nov); the commit history of constants/41b.md.
- **GitHub OmarMohammedSTEM/Ambidextrous-Sofa-Problem** (read only): README, UNIFIED_PAPER.md, src/upper_bound.py, src/interval_bound.py, data/upper_bound_results.json.
- **Zenodo API:**
  - Queries: "moving sofa", ambidextrous AND sofa, title:ambidextrous, sofa AND "upper bound", "weakened Ramsey", "set-coloring/set coloring Ramsey", "2-chromatic Ramsey", polychromatic AND Ramsey, "Ramsey turnaround", Almási Axenovich.
  - Records opened: 23667275, 21195095, 22734294, 20259705, 10636187 (PDF).
- **Bibliographic APIs:**
  - OpenAlex: works citing Xu–Shao–Su–Li (8), Harborth–Möller (13), Chung–Liu (19) and Jakhar–Moun (2), plus the abstract of the Budden book chapter.
  - Crossref: Harborth–Möller, Chung–Liu, Xu–Shao–Su–Li, Jakhar–Moun, the Budden book.
  - Semantic Scholar: citation contexts of Xu–Shao–Su–Li and Harborth–Möller.
- **Journal PDF:** Munemasa–Shinohara (jims-a.org).
- **GXL evidence API:** 3 claims (a K₆ bound above 25; an ambidextrous bound below 2.2195; constructive multigraph bounds). Nothing beyond the sources above.
- **Our files:** rg over Riemann/Deep/**/*.md NOVELTY files, REFUTATIONS, TREE and NOGO. Only R237 (the 2611.03878 refutation) is related.

## Files
- `code/`: independent checkers, all standard library or plain C, single-threaded, a few seconds each:
  - `xi_check.py` and `xi_controls.py` for (1);
  - `sofa_prep.py` + `sofa_recheck.c`, `sofa_convention.py` and `sofa_attack.c` for (2);
  - `gallai_xi.c` for the Gallai bound.

  Build with `cc -O2 -o X X.c -lm`, then run `python3 sofa_prep.py <cert.json> <dir>` and `./sofa_recheck <dir>`.
