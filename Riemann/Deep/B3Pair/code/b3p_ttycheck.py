# b3pair: exact check of TTY 2501.16779 Theorem new-exp-pair via Lemma beta-duality (alpha in [0,1/2], l-k >= 1/2),
# with and without Bourgain-derived rows of TTY Table beta-table (TeX l.437-482). Exact rationals.
from fractions import Fraction as Fr
import itertools, sys
F = Fr
def A(p): k,l = p; return (k/(2*k+2), (k+l+1)/(2*k+2))
def B(p): k,l = p; return (l-F(1,2), k+F(1,2))
# classical pairs: A/B words (no BB) of length <= D applied to (0,1)
D = int(sys.argv[1]) if len(sys.argv) > 1 else 12
pairs = {(F(0),F(1))}
front = {(F(0),F(1))}
for d in range(D):
    new = set()
    for p in front:
        new.add(A(p)); new.add(B(p))
    new -= pairs
    pairs |= new; front = new
# keep only lines minimal at some grid point of [0,1/2] (a subset of valid bounds: conservative for VALID)
grid = [i/4000 for i in range(2001)]
keep = set()
for a in grid:
    keep.add(min(pairs, key=lambda p: float(p[0]) + float(p[1]-p[0])*a))
cl = [(k, l-k, F(0), F(1,2), 'cl') for (k,l) in keep]   # beta <= k + (l-k) a on [0,1]
print('classical envelope lines kept:', len(cl))
# TTY Table beta-table rows: (intercept, slope, a0, a1, source); range [a0,a1)
T = [
 (F(1,20), F(3,4), F(0), F(1,4), 'HB k=5'),
 (F(0), F(19,20), F(1,4), F(890,3277), 'HB k=5'),
 (F(89,2706), F(2243,2706), F(890,3277), F(199,716), 'Hux17.1'),
 (F(1,66), F(235,264), F(120,419), F(754,2579), 'Hux17.1'),
 (F(9,217), F(1389,1736), F(754,2579), F(251324,841245), 'Bou:AD'),
 (F(2371,43205), F(52209,69128), F(251324,841245), F(861996,2811205), 'TY:A'),
 (F(13,146), F(47,73), F(861996,2811205), F(87,275), 'Hux17.1'),
 (F(11,244), F(191,244), F(87,275), F(423,1295), 'Hux17.1'),
 (F(89,1282), F(454,641), F(423,1295), F(227,601), 'Hux17.1'),
 (F(29,280), F(173,280), F(227,601), F(12,31), 'Hux17.1'),
 (F(1,32), F(103,128), F(12,31), F(1508,3825), 'Hux17.1'),
 (F(18,199), F(521,796), F(1508,3825), F(62831,155153), 'Bou:D'),
 (F(569,2800), F(1053,2800), F(62831,155153), F(143,349), 'Hux19.2'),
 (F(491,5530), F(1812,2765), F(143,349), F(263,638), 'Hux19.2'),
 (F(113,1345), F(897,1345), F(263,638), F(1673,4038), 'Hux19.2'),
 (F(2,9), F(1,3), F(1673,4038), F(5,12), 'Bou3.18'),
 (F(1,12), F(2,3), F(5,12), F(3,7), 'Bou3.18'),
 (F(13,84), F(1,2), F(3,7), F(1,2)+F(1,10**9), 'Bou pair'),
]
def check(pair, rows, tag):
    k,l = pair; m = l-k
    Lf = lambda a: k + m*a
    lines = rows
    cand = {F(0), F(1,2)}
    for r in lines: cand |= {r[2], r[3]}
    allL = [(r[0], r[1]) for r in lines] + [(k, m)]
    for (c1,s1),(c2,s2) in itertools.combinations(set(allL), 2):
        if s1 != s2:
            a = (c2-c1)/(s1-s2)
            if 0 <= a <= F(1,2): cand.add(a)
    cand = sorted(c for c in cand if 0 <= c <= F(1,2))
    worst = None
    # evaluate on each open subinterval's closure using lines valid on the open subinterval
    for a0, a1 in zip(cand, cand[1:]):
        mid = (a0+a1)/2
        val = [r for r in lines if r[2] <= mid < r[3]]
        if not val:
            worst = (F(10**9), a0, a1, None); break
        for a in (a0, a1):
            e = min(r[0]+r[1]*a for r in val)
            d = e - Lf(a)
            if worst is None or d > worst[0]:
                src = min(val, key=lambda r: r[0]+r[1]*a)[4]
                worst = (d, a0, a1, src)
    ok = worst[0] <= 0
    print(f"{tag}: pair ({k}, {l}) s=2k+l={2*k+l}={float(2*k+l):.6f}: max(env-L)={worst[0]} ({float(worst[0]):.3e}) on [{float(worst[1]):.5f},{float(worst[2]):.5f}] src={worst[3]} -> {'VALID' if ok else 'NOT SHOWN'}")
    return ok
new = [(F(89,1282),F(997,1282)), (F(652397,9713986),F(7599781,9713986)), (F(10769,351096),F(609317,702192)), (F(89,3478),F(15327,17390))]
print(f"classical A/B words of length <= {D}: {len(pairs)} pairs")
for p in new:
    check(p, T + cl, 'full table + classical')
nobou = [r for r in T if not r[4].startswith('Bou') and not r[4].startswith('TY')]
for p in new:
    check(p, nobou + cl, 'no Bourgain/TY rows, classical pairs')
# without the Trudgian-Yang A-row: Bourgain-derived PAIRS as full-range lines (an exponent pair bounds beta on all of [0,1])
bou_pairs = [(F(9,217), F(1461,1736)-F(9,217), F(0), F(1,2)+F(1,10**9), 'Bou AD pair'), (F(18,199), F(593,796)-F(18,199), F(0), F(1,2)+F(1,10**9), 'Bou D pair'),
             (F(13,84), F(55,84)-F(13,84), F(0), F(1,2)+F(1,10**9), 'Bou pair')]
noTY = [r for r in T if not r[4].startswith('TY')] + bou_pairs
check(new[0], noTY + cl, 'no TY row, Bourgain pairs full range')
check(new[0], [r for r in T if not r[4].startswith('TY') and not r[4].startswith('Bou')] + bou_pairs + cl, 'no TY row, Bourgain rows only as pair lines')
# also: which rows are tight for the s-optimal pair at alpha = 1/3
k,l = new[0]; print('L(1/3) =', k+(l-k)*F(1,3), '= s/3 with s =', 2*k+l)
for r in T:
    if r[2] <= F(1,3) < r[3]: print('row valid at 1/3:', r[4], r[0]+r[1]*F(1,3))
