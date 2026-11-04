# Theorem: the Fasi–Higham positivity conjecture holds for m = 3 and every k

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** For every k ≥ 0, all Taylor coefficients of 1/q_{k,3}(z) at 0 are positive, where q_{k,m}(z) = Σ_j (−1)^j C(m,j) (k+m−j)!/(k+m)! z^j is the denominator of the (k,m) Padé approximant to e^z. This is the m = 3 case of Fasi and Higham's conjecture (SIAM J. Matrix Anal. Appl. 40 (2019), Conj. 3.2: positivity for all odd m and all k).

**Before.** Proved only for k = 0 (all odd m) and m = 1; for fixed odd m ≥ 3 only for large k ("eventual positivity"); the conjecture fails for m = 5, 7, 9 and k ≫ m (Scardino, arXiv 2610.39051, 28 Oct 2026). For m = 3 there was only numerical evidence (k < 1500, Goodchild, Zenodo 23919117, 3 Nov 2026).

**Method.** With N = k + 1 and y = N − z, (k+1)(k+2)(k+3) q_{k,3} = y³ + 3Ny + 2N. Its real root gives the zero z₀ = N + d with 1/2 < d < 2/3; the two complex zeros are farther from 0 (|z_±|² − z₀² = 3N(1 − d) > 0), the zeros are simple, and positivity of every coefficient follows from the partial-fraction expansion once 2|R_±| < R₀ for the residues, which reduces to one exact polynomial identity and g(d) = −3d³ + 12d² − 30d + 12 < 0 on [1/2, 2/3). Tight: 1 − 2|R_±|/R₀ ≈ 1/(3(k+1)). A short pen-and-paper proof.

**Certificate.** [verify.py](../Riemann/Deep/FasiHigham/verify.py) (standard library + mpmath) checks the identity and the sign argument (part A), and re-checks the computational side results (parts B, C); re-run by the orchestrator on 4 Nov 2026: ALL OK.

**Review.** Correct (all steps, including simplicity of the zeros and the n = 0 case); significant: the whole m = 3 case of an open conjecture.

**Novelty.** In no source found: Fasi–Higham (2018/2019) and Fasi's thesis (no case proved), Scardino 2610.39051 (no m = 3 statement), Goodchild's Zenodo record (m = 3 numerical only), MathWorld and MathDB status pages (open); checked 4 Nov 2026 ([novelty report](../Riemann/Deep/FasiHigham/NOVELTY.md)). Our computer-assisted diagonal and first-failure data (odd m ≤ 201) are not new: Goodchild's record goes to m = 361.

**Files.** [report](../Riemann/Deep/FasiHigham/REPORT.md), [note with the proof](../Riemann/Deep/FasiHigham/NOTE.md), [verify.py](../Riemann/Deep/FasiHigham/verify.py)
