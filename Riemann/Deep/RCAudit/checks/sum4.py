import re, math, json, sys
rows = []
for line in open('l4g_out.txt'):
    nums = [float(t) for t in re.findall(r'-?\d+\.?\d*(?:E-?\d+)?', line.replace(' E', 'E'))]
    # [D, dK, f, h, X1,X2,X3, u1,u2,u3, res(9), time]
    D, dK, f, h = nums[0:4]; Xs = nums[4:7]; us = nums[7:10]; res = nums[10:19]; t = nums[19]
    rows.append(dict(D=int(D), dK=int(dK), f=int(f), h=int(h), Xs=Xs, us=us, res=res, t=t))
print(len(rows), 'D; total time', round(sum(r['t'] for r in rows), 1), 's')
def R(r, iu, ix):  # ratio to |D|^{1/2} X (|u|+3)
    n = -r['D']; return r['res'][iu*3 + ix] / (n**0.5 * r['Xs'][ix] * (abs(r['us'][iu]) + 3))
fund = [r for r in rows if r['f'] == 1 and r['h'] > 1]; nonf = [r for r in rows if r['f'] > 1 and r['h'] > 1]
for name, S in (('fundamental', fund), ('non-fundamental', nonf)):
    print('==', name, len(S))
    for iu in range(3):
        for ix in range(3):
            v = [(R(r, iu, ix), r['D']) for r in S]
            w = [(R(r, iu, ix) / math.log(3*-r['D'])**3, r['D']) for r in S]
            print(f"  u={S[0]['us'][iu]:g} X={S[0]['Xs'][ix]:g}: max S/(|D|^.5 X (u+3)) = {max(v)[0]:.3g} (D={max(v)[1]}); max of that / log^3(3|D|) = {max(w)[0]:.3g} (D={max(w)[1]})")
    # growth by decade, X=1,u=0
    import collections
    B = collections.defaultdict(list)
    for r in S: B[int(math.log10(-r['D']))].append(r)
    for k in sorted(B):
        xs = B[k]
        a = [R(r, 0, 0) for r in xs]; b = [R(r, 0, 0) / math.log(3*-r['D'])**3 for r in xs]; c = [R(r, 0, 2) for r in xs]
        print(f"  |D| decade 10^{k}: n={len(xs)}  X=1,u=0: max R={max(a):.3g} median R={sorted(a)[len(a)//2]:.3g}  max R/log^3={max(b):.3g}  | X=20: max R={max(c):.3g}")
# log-log slope fit of R (X=1,u=0) over fundamental D with |D|>=100
import statistics
xs = [math.log(-r['D']) for r in fund if -r['D'] >= 100]; ys = [math.log(R(r, 0, 0)) for r in fund if -r['D'] >= 100]
mx, my = statistics.mean(xs), statistics.mean(ys)
sl = sum((a-mx)*(b-my) for a, b in zip(xs, ys)) / sum((a-mx)**2 for a in xs)
ys3 = [math.log(R(r, 0, 0) / math.log(3*-r['D'])**3) for r in fund if -r['D'] >= 100]; my3 = statistics.mean(ys3)
sl3 = sum((a-mx)*(b-my3) for a, b in zip(xs, ys3)) / sum((a-mx)**2 for a in xs)
print('fundamental |D|>=100: slope d log R / d log|D| =', round(sl, 3), '; after dividing by log^3(3|D|):', round(sl3, 3))
