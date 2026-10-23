# The Sun–Hou conjecture holds for every finite field

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim (unconditional, complete).** For every prime power q > 7 there is a numbering a₁, …, a_q of F_q in which all q cyclic adjacent sums a_i + a_{i+1} are primitive elements. This is printed Conjecture 3.5(i) of Z.-W. Sun, with Q.-H. Hou (Nanjing Univ. J. Math. Biquarterly 36 (2019)). It is false for q = 3, 5, 7.

**Before.** Proved for q > 18,888,871 and for primes p > 61 (She, arXiv 2609.19114); Sun checked primes below 545. The prime powers 8 ≤ q ≤ 18,888,871 were open.

**Method.**
- Chvátal–Erdős (κ ≥ α) with a Hoffman-type independence bound and an exact cut inequality.
- Cohen's prime sieve with Gauss sums lowers the threshold from 18,888,871 to 6,091.
- The remaining fields: She's three-consecutive-primitive-roots construction or explicit certified cycles.
- An independent C check of every q < 30,140.

**Not claimed.** The same method proves the corrected Conjecture 3.8 (a₀ + a_i a_{i+1} primitive for all (q, a₀) ≠ (11, 2)), but that classification was posted first by an anonymous Zenodo preprint (10.5281/zenodo.23307915).

**Status.** Paper proof plus certificates; not refereed; not in Lean.

**Review.** Theorem A correct. One missing hypothesis (u ≠ 0) in a lemma used only for Conjecture 3.8, fixed. Certificates spot-checked independently (q = 9, 16, 25, 27).

**Novelty.** She's paper and the Zenodo preprint leave it open (the latter announces a "separate release").

**Files.** [report](../Riemann/Deep/SunComplete/REPORT.md), [note](../Riemann/Deep/SunComplete/NOTE.md), [verify.py](../Riemann/Deep/SunComplete/verify.py)
