# D4Pack audit notes: arXiv 2609.25120 v2 ("Sphere Packing in Dimension 4 and the Twenty-Four-Cell Conjecture")

Line numbers refer to D4.tex of the v2 e-print (10669 lines). Frames: "unit" means unit balls, with the D4 cell of volume 8; "native" means the D4 roots ±e_i±e_j are the contacts, packing distance √2 and volume × 4 = unit volume.

## Sources
- e-print: https://arxiv.org/e-print/2609.25120. The abs page (fetched 23 Oct 2026) lists only v1 (20 Sep) and v2 (24 Sep 15:41 UTC).
- Supplement: https://github.com/creelie/d4-voronoi-cells (main last pushed 23 Oct 2026; releases v1.0.0 … v1.4.0, v1.12.0, v3). v2 cites Zenodo 10.5281/zenodo.22880335, which is v1.3.0 (190 MB); the concept DOI is 10.5281/zenodo.22766562.
- Release v1.4.0 notes (24 Sep 2026 12:28 UTC) say: "v1.3.0 asserted the general bound without that hypothesis. It did so through Lemma 2.7 (shell localisation) and Lemma 2.9 (radial reduction), and both are false as stated", citing the 2.9-centre and 49-ball examples (repo independent_verification/proof_gaps/). The repo paper on main has a new abstract ("… unless 24 to 28 other centres lie within √6 of its centre, not all at distance 2") and a new Conjecture "Crowded near-contacts".
- LLM24 data: doi:10.4121/74ce1c25-6fca-4680-8a36-e9c18e7e9594, mirrored in the repo as third_party/llm24-certificate/LasserreSphericalCodes.zip.part-00 and part-01.

## Task 1: the quantitative steps of the main chain
| step | TeX | verdict |
|---|---|---|
| lem:refcell-vol (24-cell, volume 8) | 879–895 | ✓ exact (voro4.py self-test) |
| lem:hmin-new (support function ≥ 1, = 1 at roots) | 946–953 | ✓ |
| lem:shell-new (only distances in [2,2√2) matter) | 985–1002 | ✗ E1 |
| lem:root-aligned-new | 1007–1028 | ✓ as stated |
| lem:radial-new ("minimum at the all-contact corner") | 1032–1040 | ✗ E2–E5 |
| "every local configuration … is a contact configuration" | 1723 | ✗ (same gap) |
| prop:polar-form (iii), \|W\| ≤ 24 by Musin | 1727–1745 | ✓ for contact configurations only |
| lem:radial-form, thm:covering-bound, table, cor:m22 | 2016–2148 | ✓ all 20 rows recomputed |
| prop:extendable, cor:remaining | 3386–3420 | ✓ (logic) |
| lem:cap-count (sin²r ≥ (k−1)/2k) | 4746–4770 | ✓ proof re-derived |
| prop:truncated, A_4 = 7.647558 | 4781–4810 | ✓ |
| deletion three-point bound 8.140848, ω₃ = 5.65e-5 | 4853–4856 | ✗ minor: 8.14052, 6.04e-5 |
| thm:strict-reduction, A_* = 7.907144430 | 4870–4888 | ✓ |
| lem:certificate (SDP gives the pair bound) | 5033–5070 | ✓ proof re-derived |
| thm:certificate, B = 0.0929000002 | 5081–5097 | ✓ positivity and B exact; (C) sampled only |
| thm:m24 ← thm:twentyfour-points, lem:las2, lem:gram, prop:verified, lem:root-lattice | 2487–2830 | rests on the LLM24 kernel; data sample checked |
| prop:cap-reformulation, lem:crosspolytope, thm:cap-inequality (vol Q = 25/3, cap ≤ 1/3) | 6316–6560 | ✓ 25/3 exact; cap numerics |
| proof of thm:local-uncond via the corner | 6681–6690 | gap (E5) |

## Task 2: numbers
check_constants.log: 15 OK and 2 FAIL. Both failures are the ω₃ / 8.140848 item above. Our Λ₃ quadrature agrees with a 2×10⁷-point Monte Carlo to within 0.1% at three radii.

## Task 3: explicit cells (examples.py; volumes in the unit frame)
- E1 (native): roots α ≠ α₀, plus 1.45·α₀. Volume exactly 1328453/160000 = 8.30283 (deletion alone: 25/3). The new centre is at 2.9 > 2√2.
- E2: 48 centres at 1/sin 22.5° = 2.6131, in the 24 root directions and the 24 directions ±e_i, (±½)⁴. True cell 19.216162; corner 6.593944 (48 directions, max inner product 0.7071).
- E3 (native, rational): (5/4)α and (39/20)σ for σ ∈ {±e_i, (±½)⁴}, minimum squared distance 2.0525. Volume exactly 737327/40000 = 18.4332; corner 6.593944.
- E4 (native): roots except (1,1,0,0) and (1,0,1,0), plus (4/3,2/3,2/3,0). Volume exactly 407/48 = 8.4792; corner 8.231729, max inner product 0.5774.
- E5 (native): (23/20)α and (99/50)e₁. Volume exactly 174835089/12500000 = 13.9868; corner 7.941125 < 8.
- Search (search.py, seeds 1–3, 52 starts): projected gradient descent with an exact-geometry volume gradient. Best volumes per start family: perturbed D4 8.0223, D4 with holes 8.5234, 600-cell subsets 9.3592, random shells 8.9313. Nothing below 8.

## Task 4: the m = 24 step
It is explicit, not silent. TeX 2493–2511: "The theorem was first proved by de Laat, Leijenhorst and de Muinck Keizer [LLM24, Thm. 5.3] … No statement of [LLM24] is used; what is taken from that work is the certificate itself". Uniqueness of the 24-point kissing configuration is therefore claimed proved: by LLM24 (arXiv 2404.18794), computer-assisted. The brief's "open" is outdated. The intro claim "no computation that has not been checked" (TeX 298–300) still leans on that external kernel and on the authors' own scripts in the repo.

## Task 5: external data
- m = 23: certificate_d8.txt (IEEE doubles, f₀…f₈ and F₀…F₈).
  - Our check (cert23_check.py): exact LDL shows every F_k positive definite; B = 0.09290000022287 against 8 − A_* = 0.09285557029365.
  - Condition (C): minimum slack 1.175e-2 at 337 360 random admissible triples and 9.99e-3 at the 36 root triples; the value at (−1/16)³ is 2.12e-2. The margin is thin (4.4e-5 in B); their interval branch-and-bound was not rerun.
- m = 24: the reassembled zip (152 065 368 bytes) has MD5 02acd5270f7b3fa799abdeb5291706fd, as stated in TeX 2893.
  - Contents: 698 entries, 331 .txt files under 4_24, 127 dense blocks, rationals written p//q (Julia).
  - Exact LDL on 6 blocks (orders 1, 1, 1, 6, 31, 2): all positive definite. The zonal matrices and the polynomial identities were not rerun.

## Scripts (python3 + numpy + mpmath; no scipy)
- voro4.py: exact and float 4-polytope volume (vertex enumeration plus pulling triangulation).
- check_constants.py: run as check_constants.py D4.tex.
- examples.py, search.py, cap_check.py.
- cert23_check.py: takes the path to certificate_d8.txt.

## Process notes and leads
- One Monte Carlo run (removed) peaked at 947 MB, above the 300 MB rule. Every other run stayed under 170 MB and 2 min.
- Lead for a future audit: the repo-only version claims that the cells alone bound the density of every 4D packing by 0.63477, "below the three-point bound". That number is not on arXiv and is unverified here.

