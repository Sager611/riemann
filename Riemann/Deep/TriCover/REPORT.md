RECORD: 1.0633022984 with 6003 copies, exact check (was 250000000/235141779 ≈ 1.0631883499 with 303 copies, arXiv 2609.06558); limit of such chains proves λ_simplex ≥ 1.0633083237

# Covering factor of non-separable families of homothetic triangles (2026-10-22)
- **Problem.** A family of triangles is non-separable if no line that misses all of them has triangles on both sides. Its covering factor is the size of the smallest triangle covering the union, divided by the sum of the sizes. λ_simplex is the supremum over all such families; known: ≥ 1.0631883 (303 copies), ≤ 3/2.
- **Result A.** config.txt: 6003 homothetic copies, covering factor 500000000000000/470233160167623 = 1.063302298421. That is +1.139·10⁻⁴ over the old record. verify.py checks it exactly in 0.2 s with the standard library.
- **Result B (theorem).** λ_simplex ≥ 1.0633083237. Chains with equal steps at the paper's continuum parameters exist for every large number of copies, and their covering factors tend to the paper's "candidate" 1.06330832375. The paper had left open whether finite chains reach it. Proof in NOTE.md §4; limit.py certifies the numbers with interval arithmetic.
- **Same size as the paper.** With 303 copies the same construction gives 1.063188861057, above the paper's 1.063188349868; their rounding to 10⁻⁹ cost 5·10⁻⁷.
- **How.** The paper's cyclic chain recurrence, optimized for 2000 bridges per chain, recomputed with 60 digits, then every coordinate is rounded down (which enlarges each triangle) so that all 12003 required segment contacts hold strictly. The check is exact: every contact is tested in rationals, and non-separability is proved by forced 2-colouring from each colouring of the three supports.
- **Cross-check.** An exact projection test (the union is connected in every direction) agrees on the 18- and 63-copy versions.
- **Tried, no gain.** Free contact points, broken symmetry and swapped segment relations: no descent from the chain optimum. Random global search for 4 and 6 copies: weaker than the known values. The continuum optimum sits on a boundary that this pattern cannot cross cheaply.
- **Open.** Beating 1.0633083 needs a new contact pattern. The bridges use at most 53% of their projected length on real gaps, so a better pattern may exist.
- **Significance.** Result A only pushes the paper's own method further. Result B turns the paper's numerical candidate into a theorem; that is the stronger statement.
- **Files.** Riemann/Deep/TriCover/: NOTE.md, config.txt, verify.py, gen.py (numpy, scipy, mpmath), limit.py (mpmath).

