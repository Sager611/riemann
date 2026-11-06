# The double-double LDL run: per-operation error bounds (closes gap (i) in writing)

Setting: binary64, round to nearest, u = 2^-53, no overflow; a double-word (DW) number x = (xh, xl) is normalized, |xl| <= ulp(xh)/2 <= u|xh|. Code: `structured_ldl.cpp` lines 26-39 (`two_sum`, `qtwo_sum`, `operator+ - * /`), compiled with `-ffp-contract=off` (no hidden FMA), so every line is the stated rounded operation.

## Addition and multiplication: published theorems
- `operator+` (l. 29-31) is line by line Algorithm 6 (AccurateDWPlusDW) of Joldes, Muller, Popescu, "Tight and rigorous error bounds for basic building blocks of double-word arithmetic", ACM TOMS 44(2) 2017, art. 15res (hal-01351529v3): 2Sum, 2Sum, RN, Fast2Sum, RN, Fast2Sum. Theorem 3.1: relative error <= 3u^2/(1-4u) < 3u^2 + 13u^3.
- `operator*` (l. 34-35) is line by line Algorithm 10 (DWTimesDW1, QD's mul): 2Prod by `fma`, t1 = RN(xh yl), t2 = RN(xl yh), c2 = RN(t1+t2), c3 = RN(c1+c2), Fast2Sum. Theorem 5.1: relative error <= 7u^2/(1+u)^2 < 7u^2.
- Subtraction is `a + (-b)`: the same algorithm on an exactly negated DW.

## Division (QD `accurate_div`, not in JMP 2017): Lemma
`operator/` (l. 36-39): q1 = RN(ah/bh); r = a - b*(q1,0); q2 = RN(rh/bh); r' = r - b*(q2,0); q3 = RN(r'h/bh); (sh,sl) = Fast2Sum(q1,q2); return (sh,sl) + (q3,0).

**Lemma.** For normalized DW a, b with b != 0 and no underflow, the result R satisfies |R - A/B| <= (10u^2 + 103u^3)|A/B| < 10.0001 u^2 |A/B| < 2^-102.67 |A/B|, where A = ah+al, B = bh+bl.

**Proof.**
1. A = ah(1+alpha), B = bh(1+beta), |alpha|,|beta| <= u. So q1 = (A/B)(1+e1) with |e1| <= (1+u)^2/(1-u) - 1 < 3.01u.
2. P1 = B q1 exactly; the computed product is P1(1+m1), |m1| < 7u^2 (Thm 5.1, with yl = 0 allowed). The computed r is (A - P1(1+m1))(1+n1), |n1| < 3u^2+13u^3 (Thm 3.1; its error is relative to the exact sum of its two DW inputs). Since A - P1 = -A e1, |r| <= 3.02u|A|.
3. As in 1, q2 = (r/B)(1+e2), |e2| < 3.01u, so |q2| <= 3.03u|A/B|. As in 2, the computed r' = (r - P2(1+m2))(1+n2) with P2 = B q2, |r'| <= 3.02u|r| <= 9.13u^2|A|, and q3 = (r'/B)(1+e3), |q3| <= 9.2u^2|A/B|.
4. Exact identities: A - P1 = r/(1+n1) + P1 m1 and r - P2 = r'/(1+n2) + P2 m2. Dividing by B:
   A/B = q1(1+m1) + [q2(1+m2) + q3/((1+e3)(1+n2))]/(1+n1).
   Hence |A/B - (q1+q2+q3)| <= 7u^2(1+3.01u)|A/B| + 3.03u(10.1u^2)|A/B| + 9.2u^2(3.02u)|A/B| <= (7u^2 + 80u^3)|A/B|.
5. |q2| <= 3.05u|q1|, so Fast2Sum(q1,q2) is exact (Dekker) and returns a normalized DW equal to q1+q2. The final AccurateDWPlusDW adds relative error < 3u^2+13u^3 (Thm 3.1). Total: (10u^2 + 103u^3)|A/B|. QED.

## Consequence for the budget (8.5)
Every scalar operation of the dd LDL (four products, three additions and one division per off-diagonal entry; one division, products and subtractions per update) satisfies fl(x op y) = (x op y)(1+delta) with |delta| < 10.0001u^2 < 2^-102.6. `finish.py` l. 50-53 uses the declared unit roundoff 2^-(96-1) = 2^-95, which is 2^7.6 times larger. Underflow contributes absolute errors of order 2^-1074 per operation, far inside the fixed 1e-12 slack. So the dd margins 0.017327 / 0.017020 are proved, given correct execution.

The binary128 run needs none of this: GCC's `__float128` +,-,*,/ (libgcc soft-fp) are IEEE 754 correctly rounded, unit roundoff 2^-113 <= 2^-112 used, and it certifies both sectors on the same input files (same SHA-256) with margins 0.0173269704163 / 0.0170198143082.
