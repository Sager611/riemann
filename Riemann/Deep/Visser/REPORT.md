REFUTATION: §3.2 step "π(n²+n)>π(n²) is equivalent to π(x+x^{1/2})>π(x)" (TeX l.222), used to claim a prime in [x,x+x^{1/2+1/m}] "for arbitrary m>2 and 1≤x≤4.9e27" (l.230), is false. Worse, Theorem 2 itself is false for every m in (4.7078,15.831] at x just above 7 (m=15.831, x=7.5: [7.5,10.61] has no prime; RH plays no role). The intermediate claim stays unproved for 1.014e20<x≤4.97e27 when m>ln x/ln(2+x^{-1/2}) (≈66–92).

Paper: M. Visser, "Primes in short intervals assuming the Riemann hypothesis", arXiv 2610.27830v1 (20 Oct 2026). Checked against the TeX source; numbers come from check.py in this directory (18 MB, 3 s).
## 1. The step and where it is used
- Upper Oppermann gives a prime in (x,x+√x) only when x=n². π(x+√x)=π(x) on [3,3.21), [7,8.15), [13,13.35), [23,24.09), [31,31.40) and [113,116.22) (integers 3, 7, 8, 13, 23, 24, 31, 113–116).
- Used in: Theorem 2 (abstract l.111, statement l.156, proof l.220–241); Corollary and Theorem 3 (l.248–270, which is Theorem 2 with x=y^M, m=2M/(M−2)); Conclusions l.284 ("full Oppermann lets m→∞" is also false: m→∞ needs a prime in [x,x+√x], which fails at x=8). Theorem 1 does not use it.
- The same equivalence appears in the companion paper arXiv 2610.20058 ("Primes between 8th powers?", TeX l.478–481). It does no harm there, because x≥1e8 and x^{7/8} is far larger than 2√x+1.
## 2. Exact truth of Theorem 2 (real x, as the paper uses it: CMS in §2, x=y^M in Theorem 3)
- Just above a prime p with gap g to the next prime, the statement fails iff m > 1/(ln g/ln p − 1/2): p=7 (g=4) for m>4.70778; p=3 for m>7.6377; p=23 for m>13.997.
- For m=15.831 the failure set is exactly (3,3.1066)∪(7,7.8164)∪(23,23.1345). The Corollary "[y^M,y^M+y^{M−1}] contains a prime" fails for 2.289≤M<3.4772 (same x-sets, with y=x^{1/M}). Even for integer x, the "all m>2" step fails (x=114, m=25).
## 3. Rescue check (best published data)
- The exhaustive maximal-gap search now reaches at least 1.014e20 (85th maximal gap: 1854, after the prime 101412319996363309069; R. Smith 2026; OEIS A002386/A005250 b-files). Up to there, every gap after a prime p≥127 is <√p (records from 523 on satisfy g≤0.934√p).
- Sorenson–Webster (Res. Number Theory 11 (2025) 4) verified BOTH halves of Oppermann for n≤7.05e13. That is one prime in each of (n²,n²+n) and (n²+n,(n+1)²), not every gap. It gives a prime in (x,x+2√x+1] for x≤4.97e27 (the upper half alone gives only 3√x+2).
- CMS (assuming RH) gives a prime in [x,x+0.88√x ln x]. For x≤4.97e27 that interval is at least 40√x, so it is weaker than Sorenson–Webster. No explicit RH bound gets below 2√x.
- Theorem 2's own range (m≤15.831): every x≥23.14 is proved (gap data up to 1.014e20; Sorenson–Webster up to 4.97e27, since x^{1/m}≥18>2; CMS beyond). The only bad x are the three small intervals above, and there the theorem is false, not merely unproved.
- Intermediate claim (all m>2): false on the six intervals in §1, true on [116.22,1.014e20], open on 1.014e20<x≤4.97e27 for m>ln x/ln(2+x^{-1/2}) (66.5 rising to 92.0). The brief's ln x/ln 3 (≈40–58) assumed upper-half data only and gap data only up to 2^64.
## 4. What correct versions can say (assuming RH)
- Theorem 2': [x,x+x^{1/2+1/m}] contains a prime for all real x≥1 iff 0<m≤m₀=1/(ln4/ln7−1/2)=4.70778 (sharp). For real x≥23.14, or for all integers x≥1, the paper's range 0<m≤15.831 holds.
- Theorem 3 is true as stated but needs a new proof: a generalized Andrica bound q^{1/M}−p^{1/M}≤0.5104 (M=2.29, maximum at p=7) up to 1.014e20, then Sorenson–Webster, then CMS. Keeping the factor M in (y+1)^M−y^M≥M·y^{M−1}, the same inputs give every M≥2.2253.
## 5. Minor
- §2 states CMS for p_n≥3 (l.172–178), but it is false at p=3 (gap 2 > 0.88√3 ln3 = 1.67), so [x,x+0.88√x ln x] has no prime for x in (3,3.183). This is harmless, since CMS is used only for x≥4.9e27. Typos: "[y^M,y+y^{M−1}]" should be "[y^M,y^M+y^{M−1}]" (l.263–268); "m≥8/3" should be "m≤8/3" (l.289).
## 6. Novelty
arXiv shows v1 only (comment "Brief note; 6 pages", no journal ref). An arXiv search and a web search found no public note about this flaw.
