# Which three prefix reversals generate the symmetric group (below the diagonal): a complete answer

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** Let r_i reverse the first i positions of a permutation of n letters. For 2 ≤ k < m < n with m + k < n, the reversals r_n, r_m, r_k generate all of S_n if and only if n is odd, m + k = n − 1, gcd(m+1, k+1) = 1, (k+1) does not divide m, and not (n ≡ 1 mod 4 and 4 | k). On the line m + k = n − 1, with g = gcd(m+1, k+1) and α = k + 1, the full structure is determined: the orbits are the classes {x : x ≡ ±c (mod g)}, ⌊g/2⌋ + 1 of them, so the group is transitive iff g = 1; for even n and g = 1 it is imprimitive with two blocks; for odd n and g = 1 it is an imprimitive wreath product D_{2α} ≀ D_{2n/α} when α divides n, and otherwise a primitive group containing A_n (A_n exactly when n ≡ 1 mod 4 and 4 | k).

**Before.** Conjecture 4 of Blanco et al. (arXiv 2511.16959) claimed that m + k < n always gives a proper subgroup; arXiv 2610.30892 (Oct 2026) proved this off the line m + k = n − 1 with n odd, and the necessity of gcd(m+1, k+1) = 1 was known. The conjecture's literal form is false on that line (smallest counterexamples (7, 4, 2), (11, 6, 4), (13, 10, 2)).

**Method.** The product r_k r_n r_m is the first-return map of x ↦ x + α on Z/(n+1) minus 0, hence an n-cycle when gcd = 1; block systems must be residue classes, and a lap-counting argument decides imprimitivity; Jones's classification of primitive groups with an n-cycle and a fixed-point bound for collineations give A_n ⊆ G; parity decides A_n versus S_n.

**Certificate.** verify.py (30 s) checks every claim for all n ≤ 101, including A_n ⊆ G without the classification (explicit Jordan p-cycles) and group orders by sympy for n ≤ 25.

**Review.** Correct; significant as the exact classification of the open line (the bare counterexamples alone are minor).

**Novelty.** No source states the classification of the boundary line ([novelty](../Riemann/Deep/PrefixRev/NOVELTY.md)).

**Lean.** PROVED (LeanPrefix, 22 Nov 2026; Lean 4 + Mathlib, sorry-free, standard axioms): `Riemann.Deep.LeanPrefix.prefix_reversals` ([Main.lean](../Riemann/Deep/LeanPrefix/Main.lean)) states the whole Claim, from two literature inputs as hypotheses: arXiv 2610.30892 Thm 3.1 (unrefereed) and Jones 2014 Thm 1.2 (k = 0); reviewer FULLY-PROVED-AND-FAITHFUL ([STATUS](../Riemann/Deep/LeanPrefix/STATUS.md)).

**Files.** [report](../Riemann/Deep/PrefixRev/REPORT.md), [note](../Riemann/Deep/PrefixRev/NOTE.md), [verify.py](../Riemann/Deep/PrefixRev/verify.py)
