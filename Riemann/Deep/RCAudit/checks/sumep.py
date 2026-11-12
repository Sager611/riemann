import json, math, collections
r = json.load(open('ep_big.json')) + json.load(open('ep_big2.json'))
r = [x for x in r if x['h'] > 1]
B = collections.defaultdict(list)
for x in r: B[round(math.log10(-x['D']) * 4) / 4].append(x)
print('log10|D|  n  mean R=S/|D|^.5   mean R/(L1^2 log^3|D|)   mean R/(L1*log^3|D|)  mean S/(h*sum_C 1/a * log^2(sqrt|D|/(8pi)))')
prev = None
for k in sorted(B):
    xs = B[k]; n = len(xs)
    R = sum(x['R'] for x in xs)/n
    a = sum(x['R']/(x['L1']**2*math.log(-x['D'])**3) for x in xs)/n
    b = sum(x['R']/(x['L1']*math.log(-x['D'])**3) for x in xs)/n
    c = sum(x['S']/(x['h']*x['sum_inv_a']*math.log(math.sqrt(-x['D'])/(8*math.pi))**2) for x in xs if -x['D'] > 1e4)/max(1, sum(1 for x in xs if -x['D'] > 1e4))
    sl = '' if prev is None else f"  local slope dlogR/dlog|D| = {(math.log(R)-math.log(prev[1]))/((k-prev[0])*math.log(10)):.3f}"
    print(f'{k:5.2f} {n:3d}  {R:12.4g}  {a:10.4g}  {b:10.4g}  {c:8.4g}{sl}')
    prev = (k, R)
