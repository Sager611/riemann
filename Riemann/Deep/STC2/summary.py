# summary.py -- collect stc2 certificates; check A1 = ceil_5(C1) >= C1 and A3 = ceil_3(Cz + sup) against the printed Table 2
import json, re, math
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 80
KJ = json.load(open('Riemann/Deep/STC2/kernels.json'))
T2 = {'L5': ('0.20176', 0.70360002, 0.66536440, '1.369', 0.13949090, '0.844'), 'L20': ('0.10823', 1.46807534, 0.70207071, '2.171', -0.18278649, '1.286'),
      'L60': ('0.09473', 2.10426949, 0.72894156, '2.834', -0.15847575, '1.946'), 'L150': ('0.09349', 2.25578346, 0.73031685, '2.987', -0.15757634, '2.099'),
      'C400': ('0.08399', 4.46374163, 0.73756551, '5.202', -0.14582673, '4.318'), 'D6000': ('0.08382', 4.58951035, 0.73573753, '5.326', -0.14641036, '4.444')}
cz = json.load(open('Riemann/Deep/STC2/cz_arb_1e08.json'))
edge = {}
for ln in open('Riemann/Deep/STC2/certK2_edge_all.log'):
    m = re.match(r'(\w+): \(K\) CERTIFIED .* (\{.*\})$', ln.strip())
    if m: edge[m.group(1)] = json.loads(m.group(2))
box = set(re.findall(r'^(\w+): \(K\) CERTIFIED on 0<=b<=1/2, t != 0 \(2-D boxes', open('Riemann/Deep/STC2/certK2_box_all.log').read() + open('Riemann/Deep/STC2/certK2_L20_box.log').read(), re.M))
print('| line | (K) edge: min scaled H on b=1/2 / U(b,0)-pi / tail Q | (K) 2-D | Cz own <= | printed | supF N own <= | printed | supF S own <= | printed | A1 >= C1 | A3N own -> printed | A3S own -> printed |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|')
allok = True
for nm, (a1s, czp, sNp, a3n, sSp, a3s) in T2.items():
    c = KJ[nm]['c'][0]; C1 = arb(c)/(2*arb.pi()); A1 = arb(Fr(a1s).numerator)/Fr(a1s).denominator
    okA1 = bool(A1 >= C1) and bool(A1 - arb(1)/100000 < C1)
    eN = json.load(open('Riemann/Deep/STC2/elem_%s_N.json' % nm)); eS = json.load(open('Riemann/Deep/STC2/elem_%s_S.json' % nm))
    czu = cz[nm]['Cz_upper']; vN = czu + eN['sup_upper']; vS = czu + eS['sup_upper']
    okN = vN <= float(a3n) and czu <= czp and eN['sup_upper'] <= sNp; okS = vS <= float(a3s) and eS['sup_upper'] <= sSp
    rp = lambda x: math.ceil(x*1000 - 1e-9)/1000
    okP = abs(rp(czp + sNp) - float(a3n)) < 1e-9 and abs(rp(czp + sSp) - float(a3s)) < 1e-9
    e = edge[nm]; allok = allok and okA1 and okN and okS and okP and nm in box
    print('| %s | %.4e / %.4e / %.3e | %s | %.8f | %.8f | %.8f | %.8f | %.8f | %.8f | %s (C1 = %.7f) | %.5f -> %s | %.5f -> %s |' % (
        nm, min(e['b+1'], e['b-1']), e['a'], min(e['c+1'], e['c-1']), 'CERTIFIED' if nm in box else '-', czu, czp, eN['sup_upper'], sNp,
        eS['sup_upper'], sSp, 'ok' if okA1 else 'FAIL', float(C1.mid()), vN, a3n + ('' if okN else ' FAIL'), vS, a3s + ('' if okS else ' FAIL')))
print('printed A3 = ceil_3(printed Cz + printed sup) for all lines:', 'yes' if allok else 'CHECK')
print('ALL CHECKS PASS' if allok else 'SOME CHECK FAILED')

