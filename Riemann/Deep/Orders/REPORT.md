DISCOVERY: An order of + × − ÷ (start 2) has perfect sequences for all large N iff − is 3rd or 4th and the order is not ÷×−+: 11 of 24 orders; exact N-sets for all 24.

Question: arXiv 2610.10095 (TeX l. 921–925). Proved (NOTE.md, Theorem 1), exact perfect-N sets:
- N ≥ 9: + × − ÷ (the paper), + ÷ × −, × + − ÷.  N ≥ 10: + × ÷ −.  N ≥ 8: × + ÷ −.
- ÷ + × −: 2, 3, 4, N ≥ 9.  ÷ + − ×: 2, 3, 4, N ≥ 12.  + ÷ − ×: 4, 6, 10, N ≥ 12.
- × ÷ − +: 4, 5, 6, N ≥ 12.  × ÷ + −: 4–6, 8–10, N ≥ 12.  ÷ × + −: 2, 4–6, 8–10, N ≥ 12.
- Fail: ÷ × − + only 2, 4, 5; ÷ − + ×, ÷ − × + only 2; ten orders never (− first, or +− / ×− first).

How:
- Failures: short hand proofs, decided within the first five moves (NOTE.md §2).
- Working orders: from a perfect P with move P a multiplication, one uniform bulk
  (× 3(P−j), ± 1, ÷ by 2; Lemma 2) plus a 12–30-move head/tail patch per offset c gives
  3P + c. 132 patches, values linear in q (P = 4q + r), each checked for all q ≥ q* ≤ 2 by
  linear witnesses (Lemma 3); induction from N1 ≤ 53; exact search below N1.
- Small N: exhaustive search, cross-checked by plain brute force for N ≤ 14.

Checks: verify.py (about 6 s, peak about 100 MB): exact search for all orders, N ≤ 60;
symbolic and numeric checks of all patches; the construction for every N ≤ 400 and five
large N. One-off runs: construction for every N ≤ 2500 and random N < 200000; direct
search (no construction) for every 12 ≤ N ≤ 160.

Scout table corrected: N = 11 is not perfect for + ÷ − ×, × ÷ + −, ÷ × + − (N0 = 12, not
10, 8, 8); the sporadic small N above were missing.

Side result: the paper's order has a 15-move patch for c = 3 (3P + 3 from P, for every
perfect P ≡ 2 mod 4), where the paper found none and used c = 15 (23 moves).

Status: complete, computer-assisted (patch data plus checker), like the paper's own proof.
Novelty: an arXiv search found only the paper itself; a wider check is still needed.
