# Record: ξ(K₆) ≥ 33 (a 3-colouring of K₃₃ in which every K₆ sees all three colours)

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** The edges of K₃₃ can be coloured red, blue and green so that every copy of K₆ uses all three colours; equivalently the union of any two colour classes is K₆-free. Hence ξ(K₆) ≥ 33 (Almási–Axenovich parameter), the complementary / set-colouring Ramsey number R̄(6,6,6) = R(6; 3, 2) ≥ 34, and ξ(K₇) ≥ 33.

**Before.** 25 (grid colouring; Chung–Liu 1978, Thm 2.6). arXiv 2611.03878 (2 Nov 2026) claimed ξ(K₆) ≥ 31 with a cubic-residue colouring of K₃₁ that is wrong (186 two-coloured K₆; our REFUTATIONS).

**Method.** SAT and exhaustive search over symmetric colourings: circulant colourings exist only up to n = 26; a block-circulant colouring on Z₁₁ × {0, 1, 2} (cyclic shift symmetry) gives n = 33; n = 34 and 36 in the searched symmetric families are UNSAT.

**Certificate.** [K33_Z11x3.txt](../Riemann/Deep/Xi6/K33_Z11x3.txt) and [verify.py](../Riemann/Deep/Xi6/verify.py): all C(33,6) = 1,107,568 six-sets checked, none two-coloured. Orchestrator re-check with networkx: each two-colour union has clique number 5.

**Review.** Independent review (agent sidecheck): CORRECT; the definition matches every literature form (complementary, weakened, set-colouring Ramsey numbers).

**Novelty.** Nothing above 25 documented in the sources opened (Chung–Liu, Harborth–Möller, Munemasa–Shinohara 1406.2050, Conlon et al., Almási–Axenovich 2512.08053, 2611.03878, the EJC dynamic survey); one source not seen: Xu–Shao–Su–Li, Graphs Combin. 25 (2009), Table 2 (it has R̄(5,5,5) ≥ 20). Provisional ([check](../Riemann/Deep/SideCheck/REPORT.md)).

**RH relevance.** None (side result).

**Files.** [report](../Riemann/Deep/Xi6/REPORT.md), [colouring](../Riemann/Deep/Xi6/K33_Z11x3.txt), [verify.py](../Riemann/Deep/Xi6/verify.py), [review](../Riemann/Deep/SideCheck/REPORT.md)

