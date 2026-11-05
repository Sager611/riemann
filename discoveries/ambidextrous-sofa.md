# Record: the ambidextrous moving sofa constant is at most 19/10

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** Every connected rigid planar shape that can be moved through a unit-width Z-shaped hallway (a right turn and, far away, a left turn) has area at most 19/10; the same holds for Romik's formulation (right and left turn from one placement). This is constant 41b of Tao et al.'s optimization-constants list.

**Before.** Upper 2.2195 on the list page (from Baek's announced, unverified proof for the ordinary sofa, 41a), 2.37 published (Kallus–Romik 2018); lower 1.64495521 (Romik 2016).

**Method.** Kallus–Romik's finite-angle containment bound extended to two hallways (they proposed this and left it open): orientation lift plus intermediate value theorem plus Gerver's strip bounds reduce arbitrary rigid motions to two full quadrant windows of hallway angles; one rational hallway angle per window (arctan(20/21)); eight cases, each closed by branch and bound with every leaf bound recomputed exactly.

**Certificate.** [verify.py](../Riemann/Deep/OptConst/verify.py) with [sofa_certificate.json](../Riemann/Deep/OptConst/sofa_certificate.json) prints CERTIFIED (240 s, one thread); an independent tree traversal and area routine by the reviewer give the same leaf counts.

**Review.** Independent review (agent sidecheck): CORRECT, no gap (one harmless overstatement); attack by multistart local maximisation of the relaxation finds suprema ≤ 1.8298 < 1.9.

**Novelty.** No proved bound below 2.2195 anywhere (list page, its GitHub repository, Kallus–Romik follow-ups, Baek's papers, arXiv math.MG/math.OC listings 2018–2026, Zenodo; 5 Nov 2026); only non-rigorous or conditional claims ([check](../Riemann/Deep/SideCheck/REPORT.md)).

**RH relevance.** None (side result).

**Files.** [report](../Riemann/Deep/OptConst/REPORT.md), [proof](../Riemann/Deep/OptConst/PROOF.md), [verify.py](../Riemann/Deep/OptConst/verify.py), [review](../Riemann/Deep/SideCheck/REPORT.md)

