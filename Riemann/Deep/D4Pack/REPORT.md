REFUTATION: Lemmas radial-new (TeX 1032–1040) and shell-new (TeX 985–1002) are false, so Thms local-general/main-general (Δ₄ = π²/16, TeX 500–503, 615–618) are unproved; certificate: exact packings E1–E5 (examples.py).

**What fails** (arXiv 2609.25120 v2, 24 Sep 2026)
- Radial reduction (TeX 1032–1040; used at 693–697, 1352–1356, 6684): pulling every neighbour within 2√2 in to distance 2 gives a lower bound that is often not a packing. The paper then treats only "contact configurations" (TeX 1723: "Every local configuration considered in this paper … is a contact configuration") with |W| ≤ 24. Neighbours at distance in (2, 2√2) whose directions are under 60° apart are never handled.
  - E2/E3: a 49-ball packing (48 neighbours in F4 directions). True cell 19.216 (rational variant E3: exactly 737327/40000). Its corner has 48 directions, inner products up to 0.7071, and volume 6.594. That is below 7.799, the paper's floor for every contact configuration (TeX 2144).
  - E4: two adjacent contacts are replaced by one centre at 2.3094 between them. The packing is exact (min dist² = 2). True cell 407/48; the corner is not a contact configuration (max ⟨u,u'⟩ = 0.577).
  - E5: the single-deviation hypothesis holds (24 root-direction neighbours at 2.3, one non-root at 2.8). The corner has volume 7.941 < 8 although the cell is 13.987. So the proof of Thm local-uncond (TeX 6684–6690) breaks too; its statement may still be true.
- Shell localisation (TeX 985–1002): a centre at distance ≥ 2√2 can cut the cell. E1: delete one root (cell 25/3) and add a centre at 2.9 on its ray; the exact volume is 1328453/160000 = 8.3028.
- So the general bound needs exactly the case nobody has proved: cells with non-contact neighbours within 2√2. That case is the core of the open 24-cell conjecture (Musin 2018).

**Not new.** The authors' GitHub (creelie/d4-voronoi-cells) v1.4.0 release notes, dated 24 Sep 2026 12:28 UTC, say "Lemma 2.7 (shell localisation) and Lemma 2.9 (radial reduction) … are false as stated", with the same 2.9 and 49-ball examples. arXiv v2 went up 3 h later (15:41 UTC) still carrying them. The repo's current paper (main, 23 Oct 2026) states the general case as an open conjecture ("Crowded near-contacts") and claims only partial cases. arXiv still shows only v2.

**What holds** (contact configurations, all neighbours within 2√2 at distance 2)
- Constants recomputed with mpmath (check_constants.py, 15 OK): the covering table (20 rows, decreasing), A_*, 8−A_* = 0.0928555703, A_4, ω(60°) = 0.0014454086, ω₄(60°), r_*, r_4, 2r_23, and the D4 pair and triangle counts.
- m = 23 certificate (repo file multi_cap/continuation_out/certificate_d8.txt): F₀…F₈ exactly positive definite, f_k ≥ 0, and B = 0.0929000002 > 8−A_*. Condition (C) at 337 360 random admissible triples plus 36 root triples: minimum slack 0.0100 (×1000 scale), no violation. This is a sample, not a proof; we did not rerun their branch-and-bound.
- m = 24 does not silently assume uniqueness. TeX 2493–2511 uses the rational kernel of de Laat–Leijenhorst–de Muinck Keizer (arXiv 2404.18794, a computer-assisted uniqueness proof), so uniqueness is claimed proved, not open. Their data (4TU, mirrored in the repo) has MD5 02acd5270f7b3fa799abdeb5291706fd, matching TeX 2893: 331 text files and 127 dense blocks of "p//q" rationals. 6 sampled blocks are exactly positive definite.
- Cap inequality (TeX 6546): over 1804 directions the largest cap is 0.33333330 ≤ 1/3 (cap_check.py).
- Direct search for a cell below 8 (search.py): 52 projected descents from perturbed D4, D4 with holes, 600-cell subsets and random shells (21–30 neighbours). The smallest volume found was 8.0223.

**Minor:** ω₃(60°,60°,60°) is 6.04e-5, not 5.65e-5 (TeX 4856); our quadrature agrees with Monte Carlo to 0.1%. So the deletion three-point bound is 8.14052, not 8.140848 (TeX 4853). This does not affect any proof.

