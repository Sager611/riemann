# ShortPrimes — a prime in (x, x + x^{7/8}] for every real x ≥ 1 (node N10; routine corollary, R323)

**Claim.** For every real x ≥ 1 the interval (x, x + x^{7/8}] contains a prime.

**Inputs.** (I1) ζ(s) ≠ 0 for Re s > 7/8 (OpenAI 2026; Lean port [QRH](../QRH/)). (I2) Fiori–Jaskari, arXiv 2609.23222 Thm 2 (unrefereed): if every zero has Re ρ ≤ 7/8, then |ψ(x) − x| < E(x) = 0.04934 x^{3/4} log² x + 0.02636425 x^{7/8} for x > e³. (I3) ψ(y) − θ(y) < 1.4262 √y for y > 0 (Rosser–Schoenfeld 1962). (I4) Every gap between consecutive primes below 4·10^18 is ≤ 1476 (Oliveira e Silva–Herzog–Pardi, Math. Comp. 2014).

**Proof.** Put h = x^{7/8}. As ψ − θ is nondecreasing, θ(x+h) − θ(x) > h − E(x) − E(x+h) − 1.4262 √(x+h) = h(1 − r(x)). Every term of r(x) is decreasing in L = log x for L ≥ 16, and r(e^{41.56}) = 0.99992 < 1 (`short.py`, mpmath, 40 digits; r = 1 at L = 41.55892), so the claim holds for x ≥ e^{41.56} = 1.1202·10^18. For 4186.14 ≤ x < 4·10^18, (I4) gives p_{n+1} − x ≤ 1476 ≤ x^{7/8} (1476^{8/7} = 4186.13). For 2 ≤ x < 4186.14 the worst case is x = p_n; a sieve to 10^7 finds no gap p_{n+1} − p_n > p_n^{7/8}. For 1 ≤ x < 2, x + x^{7/8} ≥ 2. ∎

**Review.** Reviewer: the argument is correct given I1–I4 as stated; routine, not a discovery. Asymptotically c·x^{7/8} works for any c > 2·0.02636425 = 0.0527 (explicit cutoff grows as c ↓).

**Comparison.** Log-power intervals are shorter for moderate x: Trudgian's x/(111 log² x) until about 10^{49}, Dusart-type x/(5000 log² x) until about 10^{64}. Fiori–Jaskari state no short-interval result (TeX checked, R323).

**RH relevance:** N10 (explicit prime-error bounds): an explicit power-saving prime interval valid for all x ≥ 1, from the 7/8 half-plane; RH would give about √x log x.

Files: `short.py` (cutoff, sieve), `out.txt`.
