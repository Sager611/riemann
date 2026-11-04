# Record: kissing number in dimension 38 is at least 624020

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** τ₃₈ ≥ 624020: there are 624020 unit vectors in R³⁸ with pairwise angles of at least 60°.

**Before.** τ₃₈ ≥ 611522 (A. Kravatskiy, github.com/alexlegeartis/KissingNumbers, 28 Oct 2026; Henry Cohn's table, page of 2 Nov 2026). Upper bound 26496684 (Leijenhorst–de Laat).

**Method.** Kravatskiy's configuration in R³⁸ = R²⁴ ⊕ R¹⁴ keeps 589560 caps, 1596 axis points and 6656 latitude points; its "equator" was Q. Lin's 6855 lines (arXiv 2611.03682) of the norm-48 shell of the Leech lattice L, found by local search. We replace it by all 13104 lines of norm-48 vectors of M = {x ∈ L : Fx ≡ 0 mod 3}, a sublattice of index 3⁶ with no norm-32 vector (found through the Eisenstein structure L/θL = F₃¹²). In M any two of these vectors differ and sum to vectors of norm ≥ 48, so |v·w| ≤ 24, i.e. angles ≥ 60°; against caps the same holds inside L; latitude points are handled by Cauchy–Schwarz, axis points are orthogonal. Count: 611522 − 13710 + 26208 = 624020 (+12498). M itself is known (Nebe's extremal 3-modular lattice of dimension 24, 1998); its use as an equator is new.

**Certificate.** [verify.py](../Riemann/Deep/Kissing/verify.py) (Python standard library, exact integers, about 4 minutes): rebuilds the Golay code and the Leech lattice, checks that none of the 196560 minimal vectors lies in M, that the 13104 lines are distinct, of norm 48 and in M, and brute-forces all 85,850,856 pairs (max |v·w| = 24). Kravatskiy's own two verifiers, with only the equator file swapped ([replay_kravatskiy.py](../Riemann/Deep/Kissing/replay_kravatskiy.py)), print "VERIFIED: tau_38 >= 624020" and "VERIFIED independently". Both re-run by the orchestrator on 4 Nov 2026.

**Review.** Correct; significant (a new record in Cohn's kissing-number table).

**Novelty.** No source has τ₃₈ > 611522 or a 38D equator with more than 6855 lines: Cohn's table, the KissingNumbers repository up to commit 1f9e23d (4 Nov 2026, 21:04 UTC), Lin's repository and Zenodo record, PackingStar and other kissing repositories, arXiv math.MG/math.CO lists for Oct–Nov 2026, Zenodo; checked 4 Nov 2026, 20:57–21:20 UTC ([novelty report](../Riemann/Deep/Kissing/NOVELTY.md)).

**Files.** [report](../Riemann/Deep/Kissing/REPORT.md), [note](../Riemann/Deep/Kissing/NOTE.md), [equator lines](../Riemann/Deep/Kissing/equator_lines_13104.txt), [F (mod 3)](../Riemann/Deep/Kissing/F3.txt), [verify.py](../Riemann/Deep/Kissing/verify.py)
