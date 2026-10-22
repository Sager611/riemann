# Check of Visser, arXiv 2610.27830v1, Thm 2/3 (sec 3.2). Pure python, ~40 MB, <1 min.
import math, signal
signal.alarm(300)
N = 10**6
s = bytearray([1])*(N+1); s[0] = s[1] = 0
for i in range(2, int(N**0.5)+1):
    if s[i]: s[i*i::i] = bytearray(len(range(i*i, N+1, i)))
P = [i for i in range(N+1) if s[i]]; del s
G = list(zip(P, P[1:]))
# maximal prime gaps (p, g): OEIS A002386/A005250 b-files, 85 terms, exhaustive to >= 1.014e20
REC = [(2, 1), (3, 2), (7, 4), (23, 6), (89, 8), (113, 14), (523, 18), (887, 20), (1129, 22), (1327, 34), (9551, 36), (15683, 44), (19609, 52), (31397, 72), (155921, 86), (360653, 96), (370261, 112), (492113, 114), (1349533, 118), (1357201, 132), (2010733, 148), (4652353, 154), (17051707, 180), (20831323, 210), (47326693, 220), (122164747, 222), (189695659, 234), (191912783, 248), (387096133, 250), (436273009, 282), (1294268491, 288), (1453168141, 292), (2300942549, 320), (3842610773, 336), (4302407359, 354), (10726904659, 382), (20678048297, 384), (22367084959, 394), (25056082087, 456), (42652618343, 464), (127976334671, 468), (182226896239, 474), (241160624143, 486), (297501075799, 490), (303371455241, 500), (304599508537, 514), (416608695821, 516), (461690510011, 532), (614487453523, 534), (738832927927, 540), (1346294310749, 582), (1408695493609, 588), (1968188556461, 602), (2614941710599, 652), (7177162611713, 674), (13829048559701, 716), (19581334192423, 766), (42842283925351, 778), (90874329411493, 804), (171231342420521, 806), (218209405436543, 906), (1189459969825483, 916), (1686994940955803, 924), (1693182318746371, 1132), (43841547845541059, 1184), (55350776431903243, 1198), (80873624627234849, 1220), (203986478517455989, 1224), (218034721194214273, 1248), (305405826521087869, 1272), (352521223451364323, 1328), (401429925999153707, 1356), (418032645936712127, 1370), (804212830686677669, 1442), (1425172824437699411, 1476), (5733241593241196731, 1488), (6787988999657777797, 1510), (15570628755536096243, 1526), (17678654157568189057, 1530), (18361375334787046697, 1550), (18470057946260698231, 1552), (18571673432051830099, 1572), (20733746510561442863, 1676), (68068810283234182907, 1724), (101412319996363309069, 1854)]
L = REC[-1][0]
def root(f, lo, hi):
    for _ in range(200):
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return hi
def fails(th):  # real x in (p, x*) with no prime in [x, x+x^th]
    return [(p, round(root(lambda x: x + x**th - q, p, q), 4)) for p, q in G if p + p**th < q]
# 1. Thm 2: sharp m for all real x >= 1; failure sets
t, p0, q0 = max((math.log(q-p)/math.log(p), p, q) for p, q in G if p > 2)
print('m0 =', 1/(t-0.5), 'binding pair', (p0, q0))
for m in (15.831, 2*2.29/0.29):
    print('m=%.4f real-x failures' % m, fails(0.5+1/m),
          'integer failures', [x for p, q in G for x in range(p+1, q) if x + x**(0.5+1/m) < q])
# 2. the 'equivalent assertion' pi(x+sqrt x) > pi(x)
print('pi(x+sqrt x)=pi(x) on [p,x*):', [(p, round(root(lambda x: x+math.sqrt(x)-q, p, q), 4)) for p, q in G if p+math.sqrt(p) < q])
assert all(g <= math.sqrt(p) for p, g in REC if p >= 523)  # => gaps < sqrt p for 127 <= p <= L
# 3. CMS restated for x >= 3
h = lambda x: 0.88*math.sqrt(x)*math.log(x)
print('CMS interval empty on', [(p, round(root(lambda x: x+h(x)-q, p, q), 4)) for p, q in G if p >= 3 and p+h(p) < q])
# 4. Thm 3: generalized Andrica to L (sieve + records in force)
for M in (2.2253, 2.29):
    a = max(q**(1/M)-p**(1/M) for p, q in G)
    b = max((g/M)*p**(1/M-1) for p, g in REC if p >= 1000)
    print('M=%.4f max q^(1/M)-p^(1/M): %.4f (p<=1e6), <= %.1e beyond' % (M, a, b))
# 5. constants: x** = 4.9e27 (SW Oppermann both halves to n = 7.05e13 gives prime in (x, x+2sqrt x+1])
X = 4.9e27; lx = math.log(X)
print('Visser m(x**) =', lx/(math.log(lx)+math.log(0.88)))
print('M_min keeping factor M =', root(lambda M: math.log(M)+(0.5-1/M)*lx-math.log(0.88*lx), 2.01, 3.0))
for x in (float(L), X):
    print('x=%.3e open for m > ln x/ln(2+x^-1/2) = %.2f ; x^(1/15.831) = %.1f' % (x, math.log(x)/math.log(2+x**-0.5), x**(1/15.831)))
