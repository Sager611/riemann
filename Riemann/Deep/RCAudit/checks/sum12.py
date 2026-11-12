import json, math, sys
for fn in sys.argv[1:]:
    r = json.load(open(fn)); print('==', fn, len(r), 'rows; D range', min(x['D'] for x in r), max(x['D'] for x in r))
    def top(key, k=3):
        s = sorted(r, key=lambda x: -x[key])[:k]
        return '; '.join(f"{x[key]:.4g} (D={x['D']},f={x['f']},h={x['h']},N={x['N']})" for x in s)
    print('Lemma1 max P_Q(N)/bound:', top('L1'))
    for e in (0.5, 0.1, 0.01):
        print(f'Lemma2 eps={e}: max K/claimed:', top(f'R{e}', 2))
        print(f'   max K/(N+|D|^.5):', top(f'K{e}', 2))
    print('proof quantity h*max pi_C/(N+|D|^.5) max:', top('hpi', 2), '| vs 30 log^2(3|D|) >=', round(30*math.log(9)**2, 1))
    # growth: bin by |D| and N-type
    import collections
    B = collections.defaultdict(list)
    for x in r:
        n = -x['D']; t = 'N=sqrt' if x['N'] <= int(n**0.5) else ('N=|D|' if x['N'] == n else ('N=4|D|' if x['N'] == 4*n else ('N=|D|^1.5' if x['N'] == int(n**1.5) else 'N=|D|^2')))
        bn = 10**int(math.log10(n)) if n >= 10 else 1
        B[(t, bn)].append(x)
    print('bin(Ntype,|D| decade): max K0.1/(N+sqrt), max K0.1/log(3|D|), max unweighted/(N+sqrt), median unweighted')
    for k in sorted(B):
        xs = B[k]; u = sorted(x['unw'] for x in xs)
        print(' ', k, len(xs), round(max(x['K0.1'] for x in xs), 3), round(max(x['K0.1']/math.log(3*-x['D']) for x in xs), 3), round(max(x['unw'] for x in xs), 2), round(u[len(u)//2], 2))
