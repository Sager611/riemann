COUNTEREXAMPLES: 1811.10503 Conj 4.9 (Q and F_5: n=2, A={0,1,2}; F_11: n=3,4, A={0,1,2,5}; F_13: n=4, A={0,1,3,4}); 1208.2683 Conj 2.1 (n=9, 0<alpha<0.0023); 1208.2683 Conj 3.1 ((A,B)=(3,-26), n=5); 1211.1588 Conj 4.30(iii) (m=1, n=2); 1809.07766 Conj 6.8 (6.15) (n=7); 1402.6641 Conj 2.25 (n=5, true count of prime ideals); misprint-level: 2405.03626 Conj 3.9 (p=7), 1405.0290 Conj 4.4(ii) (n=8), 1211.1588 Conj 2.6 (4006). Scout leads 1809.07766 Conj 6.1, 6.6(ii), 1402.6641 Conj 3.5(i): MISREAD.

# Sun conjecture audit (sunaudit)

All statements read from the latest arXiv TeX (`e-print`, final versions: 1811.10503v6, 1208.2683v9, 1211.1588v29, 1809.07766v10, 1402.6641v12, 2405.03626v8, 1405.0290v5, 1601.04988v3, 2606.18234v5). Every row below is re-checked by `verify.py` (stdlib only, 0.2 s, 17 MB: ALL PASS).

| paper | conjecture | instance | proof / certificate | novelty |
|---|---|---|---|---|
| 1811.10503 (J. Algebraic Combin.) | 4.9: \|S(A)\| >= min{p(F), (\|A\|-n)n(n+1)/2 + n(n^2-1)/6 + 1}, S(A) = {sum k a_k, a_k distinct in A}, p(F) > n+1, \|A\| >= n+delta_{n,3} | Q or F_5, n=2, A={0,1,2}: S(A)={1,2,4,5}, 4 < 5. F_7 n=2 A={0,1,2,4}: 6<7. F_11 n=3 and n=4, A={0,1,2,5}: 9<11. F_13 n=4, A={0,1,3,4}: 10<11 | enumeration of all injective n-tuples. In Q only 3-term APs fail (searched A in {0..11}, \|A\|<=6, n<=6); in F_p (p<=13, all A up to affine maps; p=13, n>6 only \|A\|<=9) failures only for n=2,3,4. Also refutes Remark 4.9's group version (G=Z/p) | not in the 2 citing works (Semantic Scholar); Evidence API: nothing. New |
| 1208.2683 (Shangri-La proc.) | 2.1: for any alpha>0, (1/n) sum_{k<=n} p_k^alpha < p_n^alpha/(alpha+1) for all n >= 2ceil(alpha)^2+ceil(alpha)+6 | alpha=0.001, n=9: LHS-RHS = +2.32e-7. Fails at n=9 for all 0<alpha<0.0023364; only failure on a grid alpha in [0.001,6], n <= threshold+400 | 50-digit Decimal; first order in alpha: theta(23)/9 > log 23 - 1, i.e. 223092870 e^9 = 1.8077e12 > 23^9 = 1.8012e12 | Sun's Remark 2.1: checked only alpha=2..700. Not in 49 citing works (titles) nor Evidence API. New |
| 1208.2683 | 3.1 (2nd part): A>1, B!=0, A^2>4B, (A>2 or B>=-9), u_0=0,u_1=1,u_{n+1}=Au_n-Bu_{n-1}: (u_{n+1}^{1/(n+1)}/u_n^{1/n})_{n>=4} strictly decreasing | (A,B)=(3,-26), n=5: u_5,u_6,u_7 = 1459, 9135, 65339 give r_5 < r_6. Largest failing B: A=3: -26, A=4: -49, 5: -81, 6: -122, ..., 18: -1375 (all at n=5) | exact: u_6^70 < u_5^42 u_7^30 | Hou-Sun-Wen (1208.3903) proved only the Fibonacci case; Mansour-Shattuck (1410.6943) only "for n >= N(seq)"; Wang-Zhu (1303.5595) do not treat u_n. New |
| 1211.1588 (CANT 2017) | 4.30(iii): p^{(m+1)}_{n+1}^{1/(n+1)}/p^{(m+1)}_n^{1/n} < p^{(m)}_{n+1}^{1/(n+1)}/p^{(m)}_n^{1/n} < 1 for all m,n | m=1, n=2: 11^{1/3}/5^{1/2} > 5^{1/3}/3^{1/2} | exact: 11^2*3^3 = 3267 > 5^5 = 3125. Only failure for m=1 n<34833, m=2 n<3717, m=3 n<518, m=4 n<97 | not in 11 citing works. New |
| 1809.07766 (FFA 59, 2019) | 6.8 (6.15): det[floor(i^2 j^2/n)]_{1<=i,j<=(n-1)/2} != 0 iff n=9 or n prime > 7, n = 3 mod 4 | n=7: det[[0,0,1],[0,2,5],[1,5,11]] = -2 != 0 | exact det. Only failure for odd n<160; "prime > 3" version holds n<160 (likely typo) | Evidence API: only Sun's paper; no later work found |
| 1402.6641 (7th China-Japan sem.) | 2.25: for n>4 some k<=n has #(prime ideals of Z[i] with norm <= kn) a prime = 1 mod 4 | n=5: counts at 5,10,15,20,25 are 3,4,6,8,8. Also fails n=6,7,13,17,21 (n<=300) | brute-force Gaussian-prime enumeration. Remark 2.25 formula pi(sqrt x)+#{sqrt x<p<=x, p != 3 mod 4} counts each split prime once (true: twice); with that formula the claim holds n<=300 | statement false as written; error in Remark 2.25. No OEIS/paper. New, minor |
| 2405.03626v8 | 3.9 (p = 3 mod 4): det[x+((j+k)/p)+((j-k)/p)+(j/p)y+(k/p)z+(jk/p)w]_{0..(p-3)/2} = (-1)^{(h(-p)+1)/2} p^{(p-7)/4}(p-2m_p w)x | p=7: det = -7x - 8z + 4xw - 4yz (at (0,0,1,0): -8, formula 0) | exact det. Without the (k/p)z term the formula holds for all p = 3 mod 4, p<110, with Sun's m_7..m_59 = 2,1,-3,-1,3,1,0,8 | misprint; insignificant |
| 1405.0290 (Nanjing Univ. J.) | 4.4(ii): n>7 = k+m with p = phi(k)+phi(m/2)-1 prime, 2 a primitive root mod p | n=8: m even gives p=2,2,2 | exhaustive | misprint: OEIS A235987 (Sun) states phi(n-k)/2; that version fails only n<=5 (n<=3000). Known |
| 1211.1588 | 2.6: 2C ⊇ {4006,4007,...}, C={x>0: 2x-3, 2x+3 prime} | 4006 ∉ C+C (4006=2*2003, 4009=19*211) | exhaustive; only miss in [4000,20000] | off by one: Sun's OEIS A187759 comment lists 4006 as exception. Known |
| 1809.07766 | 6.1 (scout lead "fails p=5,11,13") | holds: p=5: s=0,t=1,rhs=1; p=13: 3+9 vs 2; p=11 is not 1 mod 4 | holds for all p = 1 mod 4, p<3000 | MISREAD |
| 1809.07766 | 6.6(ii) (6.12) (scout lead "fails p=17,41,73") | holds: p=17: 5 inversions vs 3 | holds for all odd p<3000 | MISREAD |
| 1402.6641 | 3.5(i) (scout lead "fails n=6") | 6 = 2+4: p_2,p_4,p_5 = 3,7,11 | every s in 5..20000 represented, 3,4 not | MISREAD |

## Checked, no failure (literal statements, smallest allowed cases)
- 1811.10503: 4.1 (all n-subsets of Z/5..Z/13, Z/25, Z/5^2, Z/35, Z/5^3, Z/11:Z/5, Z/31:Z/5, Heis(F_5); Z/49, Z/7^2 n<=4 all, n=5 sampled), 4.2 (n<=21), 4.3-4.5 (n<=300; 4.4(ii) fails exactly at n=7 as stated), 4.6 (n<=10), 4.7 (n<=10), 4.8 (n<=14), 4.10 (all subsets of Z/4..Z/12; Z/3^2, Z/3xZ/9, Z/5^2, Z/3^3 sizes 4..8 and full; the excluded Z/2^2xZ/2, Z/2xZ/4, Z/2xZ/6, Z/4^2 do fail).
- 1601.04988 Conj 1.2: every non-cyclic group Z/2^2 ... Z/4xZ/12, Z/7^2, Z/2^3xZ/4 (exhaustive DFS).
- 2606.18234: Conj 1.1(i) n<=6, w_2(2)=5, w_2(3)=9 (exhaustive); Conj 1.2 s_1(3..5)=2n+1 exhaustive, s_1(6)>=13, s_2(3)>=13.
- 2405.03626: 3.1-3.8, 3.9 (p = 1 mod 4), 4.1-4.9 (p<60, 7 points each), existential forms 3.7(ii), 3.8(ii), 4.6 (p<80), 4.10 (p<80), 5.1-5.8 (odd n<90).
- 1809.07766: 6.2 (all a, p<200), 6.3, 6.4, 6.5(i), 6.6, 6.8 (6.14) (p<400, n<160).
- 1405.0290: 2.1-2.4, 3.1-3.9, 4.1(i) (p,q<=600; the 15 exceptional pairs are exactly right), 4.2, right-triangle, three consecutive roots (exceptions 13,29,61 exact), gap, S_p bound and g-quad (p<=200), g+-1/g mod p^2 (p<=400, exceptions exact), 4.4(i), 5.1-5.3 and primitive-prime-divisor conjectures (n<=150; q(n) n<=600, last exception 203).
- 1208.2683: 2.2 (n<2*10^4), 2.3 (q<=45, n<=3000), 2.4 (primes <2*10^6: only a=-2, n=13), 2.5, 2.6, 2.8-2.11 (E, N_0, N exactly as claimed), 2.12-2.17, 3.1 (Fibonacci part), 3.2-3.13 (n<=400; N(r), M(r) tight).
- 1211.1588: 2.2-2.25, 3.1-3.27, 3.29-3.35, 3.37-3.39, 3.40(i), 3.43, 3.44, 4.2, 4.3(i) coincidence lists (exact), 4.17(i), 4.18(i) (n<=1500, some n<=300-800).
- 1402.6641: 2.1-2.24, 3.1-3.6, 3.8-3.15, 3.17-3.19, 3.21-3.25, 4.1-4.5, 4.6-4.10 (n<=1000 or less; 4.4(ii) n<=462).
- Not testable finitely (asymptotic, "infinitely many", "every rational"): skipped.

## Files
- `verify.py`: re-checks all rows (run `python3 verify.py`).
- `audit/`: the search scripts (stdlib Python; run inside `audit/`), one per paper (c1811*, c1601, c2606, c2405*, c1809*, c1405*, c1208*, c1211*, c1402*).
