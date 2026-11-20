/-!
# Round 510 — node N10 (Mertens/Möbius), round 17, [check] of Liflandsky arXiv 2607.09797 (9 Jul 2026; found by a node-keyword API query, not in our files): its criterion Φ(e^{−x}) = Σ μ(n)e^{−nx} = O(x^{−1/2}) ⇒ RH and simple zeros is correct and KNOWN (Báez-Duarte 2005, Thm 4.1, which the paper itself cites, l. 1347); the multiplicity and square-summability propositions and the Fejér endpoint characterization (under RH + simple zeros) hold; the residue and explicit-formula numerics did not fit a scout; Mermin–Wagner, Papoulis sampling and Beale–Kato–Majda transfer no tool

* Node: N10 (17th round; R501 before). Round type: [check] of an unseen on-node paper from node-keyword arXiv API queries over the last six months (lesson 350; the N4, N6, N2, N3 and N7 queries found nothing new on-node), plus three exploration scouts.
* (1) crit: thm:criterion (TeX l. 190–200, proof l. 1169–1253): O(x^{−1/2}) gives G(s) = Γ(s)/ζ(s) = O(1/η) at σ = 1/2 + η, so no zeros right of the line and poles of order m ≤ 1 (G(ρ + η) ≍ η^{−m}). Correct; KNOWN: Báez-Duarte, 'Möbius convolutions and the Riemann hypothesis', IJMMS 2005, Thm 4.1, which the paper's own remark (l. 1347) says specializes exactly to this transform.
* (2) mult: prop:multlog (l. 1254–1288, logarithmic losses bound multiplicities) and prop:residueL2 (l. 1290–1340, square-summability of the residues from the endpoint bound): Mellin exponents and Abel diagonal limits check. HOLDS.
* (3) resid (double poles of Γ(s)x^{−s}/ζ(s) at s = −2m): no answer within the scout's limit (a 30-digit contour quadrature). NOT TESTABLE.
* (4) numer (thm:main at x = 1, 3 against Σ_{n ≤ 4000} μ(n)e^{−nx} with 60 zeros): no answer within the scout's limit. NOT TESTABLE.
* (5) fejer: thm:fejer (l. 1006–1017): under RH and simple zeros (S), the endpoint bound O(x^{−1/2}) ⟺ the Fejér condition (H_F); globally RH + (S) + (H_F), not RH alone; the Fejér kernel K_T(v) = (1 − cos Tv)/(πTv²) has mass 1 and a triangular transform (checked). HOLDS.
* Exploration: (a) the Mermin–Wagner infrared bound → N13: forbids order only for a hypothetical translation-invariant Goldstone dispersion in d ≤ 2 (R048 tried statistical mechanics): no transfer; (b) Papoulis' generalized (two-channel) sampling → N13: needs Ξ band-limited, which it is not: no transfer; (c) the Beale–Kato–Majda criterion → N13: no canonical nonlinear evolution or enstrophy for the zeros (Kato's criterion ruled out in R452): no transfer.
-/

namespace Fayed.R510
end Fayed.R510
