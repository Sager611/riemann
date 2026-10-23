# Which cyclic orders of the four operations make every number: a complete answer

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** Start from the number 2. Each move combines two numbers already made, using the next operation of a fixed cyclic order of +, ×, −, ÷ (results positive, divisions exact), and must make a new number at most N; a sequence is perfect if after N − 1 moves every number 1, …, N has been made. An order has perfect sequences for all large N if and only if − is third or fourth in it and it is not ÷ × − +: 11 of the 24 orders. The exact set of perfect N is determined for every order.

**Before.** arXiv 2610.10095 (Oct 2026) proved that the order + × − ÷ works for every N ≥ 9 and asked which other orders work (TeX l. 921–925).

**Method.** Thirteen orders fail by short arguments about the first few moves. For the eleven working orders, a uniform bulk construction plus 132 head and tail patches, each checked symbolically for every q ≥ 2, turns a perfect sequence for P into one for 3P + c; induction starts at N ≤ 53, with exhaustive search below.

**Certificate.** verify.py (6 s): exact search for all orders up to N = 60, symbolic checks of all patches, the constructions checked for every N ≤ 400 and at five large N.

**Review.** Correct as a computer-assisted proof, like the paper's own; significant as a complete answer to the question.

**Novelty.** Nothing found beyond the paper itself ([novelty](../Riemann/Deep/Orders/NOVELTY.md)).

**Files.** [report](../Riemann/Deep/Orders/REPORT.md), [note](../Riemann/Deep/Orders/NOTE.md), [patches](../Riemann/Deep/Orders/PATCHES.md), [verify.py](../Riemann/Deep/Orders/verify.py)
