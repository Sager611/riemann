# Record: a quaternary covering code with K₄(10,3) ≤ 768

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** There is a code of 768 words of length 10 over a 4-letter alphabet such that every one of the 4¹⁰ words is within Hamming distance 3 of a codeword. So K₄(10,3) ≤ 768.

**Before.** K₄(10,3) ≤ 832 (Bertolo, di Nasso, Di Pasquale, Jurcovich, Santisi 2004, in G. Kéri's covering-code tables); lower bound 313.

**Method.** The octacode (256 words, K₄(8,2) ≤ 256; Maharaj, arXiv 2610.03760, Sept 2026) amalgamated on its last coordinate with a 12-word code of length 3 found by a set-cover search. The textbook amalgamated-direct-sum theorem does not apply (the octacode is not normal), so the search was needed.

**Certificate.** verify.py checks the covering radius exhaustively in 0.1 s (standard library only); an independent C program agrees word by word.

**Review.** Certificate correct; a modest but genuine record (one table entry, −7.7%).

**Novelty.** Kéri's tables, newer tables and repositories, and this year's papers list 832 or nothing lower ([novelty](../Riemann/Deep/CovCodes/NOVELTY.md)).

**Files.** [code](../Riemann/Deep/CovCodes/codes/K4_10_3_768.txt), [verify.py](../Riemann/Deep/CovCodes/verify.py)
